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

namespace FOOD
{

template<class ScalarType_T>
class TensorField
{

  public:

    //@{
    //! Typedefs.
    typedef ScalarType_T                             ScalarType;
    typedef int                                      OrdinalType;
    typedef Teuchos::Comm<int>                       Communicator_t;
    typedef Teuchos::RCP<const Communicator_t>       RCP_Communicator;
    typedef Teuchos::RCP<Domain>                     RCP_Domain;
    typedef Teuchos::RCP<TensorTemplate>             RCP_TensorTemplate;
    typedef Teuchos::RCP<Unit>                       RCP_Unit;
    typedef Tpetra::Map<OrdinalType>                 Tpetra_Map_t;
    typedef Teuchos::RCP<const Tpetra_Map_t>         RCP_Tpetra_Map;
    typedef Teuchos::ArrayView<ScalarType>           View;
    typedef Teuchos::ArrayView<const ScalarType>     ConstView;
    typedef int                                      ErrorCode;
    //@}

  private:

    // The communicator this field is defined on.
    RCP_Communicator d_comm;

    // The degrees of freedom represented by this field. Interleaved component
    // storage.
    Teuchos::ArrayRCP<ScalarType> d_dofs;

    // Tpetra map for the degrees of freedom represented by this field.
    RCP_Tpetra_Map d_dof_map;

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
		 int entity_type,
		 int entity_topology,
		 int coord_type,
		 RCP_TensorTemplate tensor_template,
		 RCP_Unit unit,
		 const std::string &name );

    // Destructor.
    ~TensorField();

    // Attach this field to tag data.
    void attachToTagData( iBase_TagHandle dof_tag,
			  ErrorCode &error);

    // Attach this field to array data and tag the mesh.
    void attachToArrayData( Teuchos::ArrayRCP<ScalarType> dof_array,
			    int storage_order,
			    ErrorCode &error);

    // Evaluate the degrees of freedom of this field at a set of coordinates
    // in a particular entity.
    Teuchos::ArrayRCP<ScalarType> evaluateDF( iBase_EntityHandle entity,
					      Teuchos::Tuple<double,3> coords);

    // Evaluate gradient of the degrees of freedom of this field at a set of
    // coordinates in a particular entity. 
    Teuchos::ArrayRCP<ScalarType> 
    evaluateGradDF( iBase_EntityHandle entity,
		    Teuchos::Tuple<double,3> coords);

    // Evaluate the Hessian of the degrees of freedom of this field at a set
    // of coordinates in a particular entity. 
    Teuchos::ArrayRCP<ScalarType> 
    evaluateHessianDF( iBase_EntityHandle entity,
		       Teuchos::Tuple<double,3> coords);

    //! Get a view of all the degrees of freedom for this field.
    View getTensorFieldDFView()
    { return View(d_dofs); }

    //! Get a const view of all the degrees of freedom for this field.
    ConstView getTensorFieldDFConstView() const
    { return View(d_dofs); }

    // Get a component view of the degrees of freedom for this field.
    View getTensorFieldComponentView(int component);

    // Get a const component view of the degrees of freedom for this field.
    ConstView getTensorFieldConstComponentView(int component) const;

    // Get const degrees of freedom for a particular entity in the domain.
    ConstView getTensorFieldConstEntDF( iBase_EntityHandle entity,
					ErrorCode &error );

    // Get const degrees of freedom for an array of entites in the
    // domain. Returned implicitly interleaved. 
    ConstView
    getTensorFieldConstEntArrDF( Teuchos::ArrayRCP<iBase_EntityHandle> entities );

    //! Get the Tpetra map for the degrees of freedom.
    RCP_Tpetra_Map getTensorFieldDFMap() const
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

