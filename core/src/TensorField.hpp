//---------------------------------------------------------------------------//
// \file TensorField.hpp
// \author Stuart Slattery
// \brief Tensor field declaration.
//---------------------------------------------------------------------------//

#ifndef FOOD_TENSORFIELD_HPP
#define FOOD_TENSORFIELD_HPP
 
#include <string>
#include <vector>

#include "Types.hpp"
#include "TypeTraits.hpp"
#include "Unit.hpp"
#include "Domain.hpp"
#include "TensorTemplate.hpp"

#include <iMesh.h>
#include <iBase.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include <Tpetra_Map.hpp>

#include <Intrepid_FieldContainer.hpp>

namespace FOOD
{

template<class Scalar>
class TensorField
{

  public:

    //@{
    //! Typedefs.
    typedef unsigned long int                        OrdinalType;
    typedef Teuchos::Comm<int>                       Communicator_t;
    typedef Teuchos::RCP<const Communicator_t>       RCP_Communicator;
    typedef Teuchos::RCP<Domain>                     RCP_Domain;
    typedef Teuchos::RCP<TensorTemplate>             RCP_TensorTemplate;
    typedef Teuchos::RCP<Unit>                       RCP_Unit;
    typedef Tpetra::Map<OrdinalType>                 Map_t;
    typedef Teuchos::RCP<const Map_t>                RCP_Map;
    typedef Intrepid::FieldContainer<Scalar>         MDArray;
    typedef int                                      ErrorCode;
    //@}

  private:

    // The communicator this field is defined on.
    RCP_Communicator d_comm;

    // The degrees of freedom represented by this field. Stored in a
    // multidimensional array. Shallow copy of tag data.
    MDArray d_dofs;

    // Tpetra map for the degrees of freedom represented by this field.
    RCP_Map d_dof_map;

    // The domain this field is defined on.
    RCP_Domain d_domain;

    // The entity type this field is defined on.
    std::size_t d_entity_type;

    // The entity topology this field is defined on.
    std::size_t d_entity_topology;

    // The coordinate system for physical field coordinates.
    std::size_t d_coord_type;

    // The tensor template for this field.
    RCP_TensorTemplate d_tensor_template;

    // The units for this field.
    RCP_Unit d_unit;

    // The name of this field.
    std::string d_name;

    // Degrees of freedom tag on the mesh.
    iBase_TagHandle d_dof_tag;

  public:

    // Constructor.
    TensorField( RCP_Communicator comm,
		 RCP_Domain domain,
		 const int entity_type,
		 const int entity_topology,
		 const int coord_type,
		 RCP_TensorTemplate tensor_template,
		 RCP_Unit unit,
		 const std::string &name );

    // Destructor.
    ~TensorField();

    // Attach this field to tag data.
    void attachToTagData( iBase_TagHandle dof_tag,
			  ErrorCode &error );

    // Attach this field to array data and tag the mesh.
    void attachToArrayData( Teuchos::ArrayRCP<Scalar> dof_array,
			    int storage_order,
			    ErrorCode &error );

    // Evaluate the degrees of freedom of this field at a set of coordinates
    // in a particular entity.
    Teuchos::ArrayView<Scalar> 
    evaluateDF( iBase_EntityHandle entity,
		Teuchos::Tuple<double,3> coords,
		int is_param );

    // Evaluate gradient of the degrees of freedom of this field at a set of
    // coordinates in a particular entity. 
    Teuchos::ArrayView<Scalar> 
    evaluateGradDF( iBase_EntityHandle entity,
		    Teuchos::Tuple<double,3> coords,
		    int is_param );

    // Evaluate the Hessian of the degrees of freedom of this field at a set
    // of coordinates in a particular entity. 
    Teuchos::ArrayView<Scalar> 
    evaluateHessianDF( iBase_EntityHandle entity,
		       Teuchos::Tuple<double,3> coords,
		       int is_param );

    //! Get the communicator this tensor field is defined on.
    RCP_Communicator getTensorFieldComm() const
    { return d_comm; }

    //! Get a view of all the degrees of freedom for this field.
    Teuchos::ArrayView<Scalar> getTensorFieldDFView()
    { return Teuchos::ArrayView<Scalar>( d_dofs.getData() ); }

    //! Get a const view of all the degrees of freedom for this field.
    Teuchos::ArrayView<const Scalar> getTensorFieldDFConstView() const
    { return Teuchos::ArrayView<const Scalar>( d_dofs.getData() ); }

    // Get a component view of the degrees of freedom for this field.
    Teuchos::ArrayView<Scalar> 
    getTensorFieldComponentView( int component );

    // Get a const component view of the degrees of freedom for this field.
    Teuchos::ArrayView<const Scalar>
    getTensorFieldConstComponentView( int component ) const;

    // Get const degrees of freedom for a particular entity in the domain.
    Teuchos::ArrayRCP<const Scalar> 
    getTensorFieldConstEntDF( iBase_EntityHandle entity,
			      ErrorCode &error ) const;

    // Get const degrees of freedom for an array of entities in the
    // domain. Returned implicitly interleaved. 
    Teuchos::ArrayRCP<const Scalar> 
    getTensorFieldConstEntArrDF( iBase_EntityHandle *entities,
				 int num_entities,
				 ErrorCode &error ) const;

    //! Get the Tpetra map for the degrees of freedom.
    RCP_Map getTensorFieldDFMap() const
    { return d_dof_map; }

    //! Get the domain this field is defined on.
    RCP_Domain getTensorFieldDomain() const
    { return d_domain; }

    //! Get the entity type this field is defined on.
    int getTensorFieldEntityType() const
    { return d_entity_type; }

    //! Get the entity topology this field is defined on.
    int getTensorFieldEntityTopology() const
    { return d_entity_topology; }

    //! Get the coordinate system for physics field coordinates.
    int getTensorFieldCoordType() const
    { return d_coord_type; }

    //! Get the tensor template for this field.
    RCP_TensorTemplate getTensorFieldTemplate() const
    { return d_tensor_template; }

    //! Get the unit for this field.
    RCP_Unit getTensorFieldUnit() const
    { return d_unit; }

    //! Get the name of this field.
    const std::string& getTensorFieldName() const
    { return d_name; }

    //! Get the degrees of freedom tag on the mesh this field is associated
    //! with.
    iBase_TagHandle getTensorFieldDFTag() const
    { return d_dof_tag; }

  private:

    // Map the degrees of freedom.
    void mapDF();
}; 

} // end namespace FOOD

#include "TensorField_Def.hpp"

#endif // end FOOD_TENSORFIELD_HPP

//---------------------------------------------------------------------------//
// end TensorField.hpp
//---------------------------------------------------------------------------//

