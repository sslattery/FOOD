//---------------------------------------------------------------------------//
// \file TensorField.hpp
// \author Stuart Slattery
// \brief Tensor field declaration.
//---------------------------------------------------------------------------//

#ifndef FOOD_TENSORFIELD_HPP
#define FOOD_TENSORFIELD_HPP
 
#include <string>

#include "Types.hpp"
#include "TypeTraits.hpp"
#include "Unit.hpp"
#include "Domain.hpp"
#include "DFuncKernel.hpp"
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
    typedef DFuncKernel<Scalar>                      DFuncKernel_t;
    typedef Teuchos::RCP<DFuncKernel_t>              RCP_DFuncKernel;
    typedef Teuchos::RCP<TensorTemplate>             RCP_TensorTemplate;
    typedef Teuchos::RCP<Unit>                       RCP_Unit;
    typedef Tpetra::Map<OrdinalType>                 Map_t;
    typedef Teuchos::RCP<const Map_t>                RCP_Map;
    typedef Intrepid::FieldContainer<Scalar>         MDArray;
    typedef int                                      ErrorCode;
    typedef iBase_EntityHandle                       EntityHandle;
    typedef iBase_TagHandle                          TagHandle;
    //@}

  private:

    // The communicator this field is defined on.
    RCP_Communicator d_comm;

    // The degrees of freedom represented by this field. Stored in a
    // multidimensional array access wrapper. Shallow copy of tag data.
    MDArray d_dofs;

    // Tpetra map for the degrees of freedom represented by this field.
    RCP_Map d_dof_map;

    // The domain this field is defined on.
    RCP_Domain d_domain;

    // The distribution function kernel to be used to evaluate this field.
    RCP_DFuncKernel d_dfunckernel;

    // The coordinate system for physical field coordinates.
    std::size_t d_coord_type;

    // The tensor template for this field.
    RCP_TensorTemplate d_tensor_template;

    // The units for this field.
    RCP_Unit d_unit;

    // The name of this field.
    std::string d_name;

    // Degrees of freedom tag on the mesh.
    TagHandle d_dof_tag;

  public:

    // Constructor.
    TensorField( RCP_Communicator comm,
		 RCP_Domain domain,
		 RCP_DFuncKernel dfunckernel,
		 const int coord_type,
		 RCP_TensorTemplate tensor_template,
		 RCP_Unit unit,
		 const std::string &name );

    // Destructor.
    ~TensorField();

    // Attach this field to tag data.
    void attachToTagData( TagHandle dof_tag,
			  ErrorCode &error );

    // Attach this field to array data and tag the mesh.
    void attachToArrayData( Teuchos::ArrayRCP<Scalar> dof_array,
			    int storage_order,
			    ErrorCode &error );

    // Evaluate the degrees of freedom of this field at a set of coordinates
    // in a particular entity.
    void evaluateDF( const EntityHandle entity,
		     const MDArray &coords,
		     const int is_param,
	             MDArray &dfunc_values );

    // Evaluate gradient of the degrees of freedom of this field at a set of
    // coordinates in a particular entity. 
    void evaluateGradDF( const EntityHandle entity,
			 const MDArray &coords,
			 const int is_param,
			 MDArray &dfunc_values );

    // Evaluate the Hessian of the degrees of freedom of this field at a set
    // of coordinates in a particular entity. 
    void evaluateHessianDF( const EntityHandle entity,
			    const MDArray &coords,
			    const int is_param,
			    MDArray &dfunc_values );

    //! Get the communicator this tensor field is defined on.
    RCP_Communicator getComm() const
    { return d_comm; }

    //! Get the degrees of freedom for this field.
    MDArray getDF()
    { return d_dofs; }

    //! Get const degrees of freedom for this field.
    MDArray getConstDF()
    { return d_dofs; }

    //! Get a view of all the degrees of freedom for this field.
    Teuchos::ArrayView<Scalar> getDFView()
    { return d_dofs.getData()(); }

    //! Get a const view of all the degrees of freedom for this field.
    Teuchos::ArrayView<const Scalar> getDFConstView() const
    { return d_dofs.getData()(); }

    // Get degrees of freedom for a particular entity in the domain.
    MDArray getEntDF( EntityHandle entity, 
		      ErrorCode &error ) const;

    // Get degrees of freedom for an array of entities in the
    // domain. Returned implicitly interleaved. 
    MDArray getEntArrDF( EntityHandle *entities,
			 int num_entities,
			 ErrorCode &error ) const;

    //! Get the Tpetra map for the degrees of freedom.
    RCP_Map getDFMap() const
    { return d_dof_map; }

    //! Get the domain this field is defined on.
    RCP_Domain getDomain() const
    { return d_domain; }

    //! Get the distribution function kernel this field is defined on.
    RCP_DFuncKernel getDFuncKernel() const
    { return d_dfunckernel; }

    //! Get the coordinate system for physics field coordinates.
    int getCoordType() const
    { return d_coord_type; }

    //! Get the tensor template for this field.
    RCP_TensorTemplate getTensorTemplate() const
    { return d_tensor_template; }

    //! Get the unit for this field.
    RCP_Unit getUnit() const
    { return d_unit; }

    //! Get the name of this field.
    const std::string& getName() const
    { return d_name; }

    //! Get the degrees of freedom tag on the mesh this field is associated
    //! with.
    TagHandle getDFTag() const
    { return d_dof_tag; }

  private:

    // Map the degrees of freedom with globally unique ID's.
    void mapDF();
}; 

} // end namespace FOOD

#include "TensorField_Def.hpp"

#endif // end FOOD_TENSORFIELD_HPP

//---------------------------------------------------------------------------//
// end TensorField.hpp
//---------------------------------------------------------------------------//

