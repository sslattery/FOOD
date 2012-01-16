//---------------------------------------------------------------------------//
// \file TensorField.hpp
// \author Stuart Slattery
// \brief Tensor field definition.
//---------------------------------------------------------------------------//

#ifndef FOOD_TENSORFIELD_HPP
#define FOOD_TENSORFIELD_HPP
 
#include <string>

#include "Types.hpp"
#include "Unit.hpp"
#include "Domain.hpp"
#include "TensorTemplate.hpp"

#include <iMesh.h>
#include <iBase.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>

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
    typedef Teuchos::RCP<Domain>                     RCP_Domain;
    typedef Teuchos::RCP<TensorTemplate>             RCP_TensorTemplate;
    typedef Teuchos::RCP<Unit>                       RCP_Unit;
    typedef Teuchos::ArrayRCP<ScalarType>            DOFArray;
    typedef Tpetra::Map<OrdinalType>                 Tpetra_Map_t;
    typedef Teuchos::RCP<Tpetra_Map_t>               RCP_Tpetra_Map;
    typedef int                                      ErrorCode;
    //@}

  private:

    // The degrees of freedom represented by this field.
    DOFArray d_dofs;

    // Tpetra map for the degrees of freedom represented by this field.
    RCP_Tpetra_Map d_dof_map;

    // The domain this field is defined on.
    RCP_Domain d_domain;

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
    TensorField( RCP_Domain domain,
		 int coord_type,
		 RCP_TensorTemplate tensor_template,
		 RCP_Unit unit,
		 const std::string &name );

    // Destructor.
    ~TensorField();

    // Attach this field to tag data.
    ErrorCode attachToTagData( iBase_TagHandle dof_tag,
			       ScalarType untagged_values );

    // Attach this field to array data.
    ErrorCode attachToArrayData( DOFArray dof_array,
				 int storage_order );

    //! Get the degrees of freedom for this field.
    DOFArray getTensorFieldDF() const
    { return d_dofs; }

    //! Get the Tpetra map for the degrees of freedom.
    RCP_Tpetra_Map getTensorFieldMap() const
    { return d_dof_map; }

    //! Get the domain this field is defined on.
    RCP_Domain getTensorFieldDomain() const
    { return d_domain; }

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

