//---------------------------------------------------------------------------//
// \file TensorField.hpp
// \author Stuart Slattery
// \brief Tensor field definition.
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
    typedef Teuchos::Comm<OrdinalType>               Communicator_t;
    typedef Teuchos::RCP<const Communicator_t>       RCP_Communicator;
    typedef Teuchos::RCP<Domain>                     RCP_Domain;
    typedef Teuchos::RCP<TensorTemplate>             RCP_TensorTemplate;
    typedef Teuchos::RCP<Unit>                       RCP_Unit;
    typedef Tpetra::Map<OrdinalType>                 Tpetra_Map_t;
    typedef Teuchos::RCP<const Tpetra_Map_t>         RCP_Tpetra_Map;
    typedef int                                      ErrorCode;
    //@}

  private:

    // The communicator this field is defined on.
    RCP_Communicator d_comm;

    // The degrees of freedom represented by this field. Interleaved storage.
    Teuchos::ArrayRCP<ScalarType> d_dofs;

    // Tpetra map for the degrees of freedom represented by this field.
    RCP_Tpetra_Map d_dof_map;

    // The domain this field is defined on.
    RCP_Domain d_domain;

    // The entity type this field is defined on.
    int d_entity_type;

    // The entity topology this field is defined on.
    int d_entity_topology;

    // The coordinate system for physical field coordinates.
    int d_coord_type;

    // The tensor template for this field.
    RCP_TensorTemplate d_tensor_template;

    // The units for this field.
    RCP_Unit d_unit;

    // The name of this field.
    std::string d_name;

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
    ErrorCode attachToTagData( iBase_TagHandle dof_tag );

    // Attach this field to array data.
    ErrorCode attachToArrayData( Teuchos::ArrayView<ScalarType> dof_array );

    //! Get a view of all the degrees of freedom for this field.
    Teuchos::ArrayView<ScalarType> getTensorFieldDFView()
    { return Teuchos::ArrayView<ScalarType>(d_dofs); }

    //! Get a const view of all the degrees of freedom for this field.
    const Teuchos::ArrayView<const ScalarType> getTensorFieldDFConstView() const
    { return Teuchos::ArrayView<ScalarType>(d_dofs); }

    //! Get a component view of the degrees of freedom for this field.
    Teuchos::ArrayView<ScalarType> getTensorFieldComponentView(int component);

    //! Get a const component view of the degrees of freedom for this field.
    const Teuchos::ArrayView<const ScalarType> 
    getTensorFieldConstComponentView(int component) const;

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
}; 

} // end namespace FOOD

#include "TensorField_Def.hpp"

#endif // end FOOD_TENSORFIELD_HPP

//---------------------------------------------------------------------------//
// end TensorField.hpp
//---------------------------------------------------------------------------//

