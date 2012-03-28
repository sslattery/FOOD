//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file TensorField.hpp
 * \author Stuart Slattery
 * \brief Tensor field declaration.
 */
//---------------------------------------------------------------------------//

#ifndef FOOD_TENSORFIELD_HPP
#define FOOD_TENSORFIELD_HPP
 
#include <string>

#include "FieldTypes.hpp"
#include "TypeTraits.hpp"
#include "Unit.hpp"
#include "Domain.hpp"
#include "DFuncKernel.hpp"
#include "TensorTemplate.hpp"

#include <iMesh.h>
#include <iBase.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ArrayRCP.hpp>

namespace FOOD
{

template<class Scalar>
class TensorField
{

  public:

    //@{
    //! Typedefs.
    typedef unsigned long int                        OrdinalType;
    typedef Teuchos::RCP<Domain>                     RCP_Domain;
    typedef DFuncKernel<Scalar>                      DFuncKernel_t;
    typedef Teuchos::RCP<DFuncKernel_t>              RCP_DFuncKernel;
    typedef Teuchos::RCP<TensorTemplate>             RCP_TensorTemplate;
    typedef Teuchos::RCP<Unit>                       RCP_Unit;
    //@}

  private:

    // The degrees of freedom represented by this field. Stored in a
    // multidimensional array access wrapper. Shallow copy of tag data.
    Teuchos::ArrayRCP<Scalar> d_dofs;

    // The domain this field is defined on.
    RCP_Domain d_domain;

    // The distribution function kernel to be used to evaluate this field.
    RCP_DFuncKernel d_dfunckernel;

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
    TensorField( RCP_Domain domain,
		 RCP_DFuncKernel dfunckernel,
		 RCP_TensorTemplate tensor_template,
		 const std::string &name );

    // Destructor.
    ~TensorField();

    // Attach this field to tag data.
    void attachToTagData( iBase_TagHandle dof_tag,
			  int &error );

    // Attach this field to array data and tag the mesh.
    void attachToArrayData( const Teuchos::ArrayRCP<Scalar> &dof_array,
			    const int storage_order,
			    int &error );

    // Evaluate the degrees of freedom of this field at a set of coordinates
    // in a particular entity.
    void evaluateDF( const iBase_EntityHandle entity,
		     const double coords[3],
		     const int is_param,
	             Teuchos::ArrayRCP<Scalar> &dfunc_values );

    // Evaluate gradient of the degrees of freedom of this field at a set of
    // coordinates in a particular entity. 
    void evaluateGradDF( const iBase_EntityHandle entity,
			 const double coords[3],
			 const int is_param,
			 Teuchos::ArrayRCP<Scalar> &dfunc_values );

    //! Get the degrees of freedom for this field.
    Teuchos::ArrayRCP<Scalar> getDF()
    { return d_dofs; }

    //! Get const degrees of freedom for this field.
    Teuchos::ArrayRCP<Scalar> getConstDF() const
    { return d_dofs; }

    //! Get a view of all the degrees of freedom for this field.
    Teuchos::ArrayView<Scalar> getDFView()
    { return d_dofs(); }

    //! Get a const view of all the degrees of freedom for this field.
    Teuchos::ArrayView<const Scalar> getDFConstView() const
    { return d_dofs(); }

    // Get degrees of freedom for a particular entity in the domain.
    Teuchos::ArrayRCP<Scalar> getEntDF( iBase_EntityHandle entity ) const;

    // Get degrees of freedom for an array of entities in the
    // domain. Returned implicitly interleaved. 
    Teuchos::ArrayRCP<Scalar> getEntArrDF( iBase_EntityHandle *entities,
					   int num_entities ) const;

    //! Get the domain this field is defined on.
    RCP_Domain getDomain() const
    { return d_domain; }

    //! Get the distribution function kernel this field is defined on.
    RCP_DFuncKernel getDFuncKernel() const
    { return d_dfunckernel; }

    //! Get the tensor template for this field.
    RCP_TensorTemplate getTensorTemplate() const
    { return d_tensor_template; }

    //! Set the unit for this field.
    void setUnit( RCP_Unit unit )
    { d_unit = unit; }

    //! Get the unit for this field.
    RCP_Unit getUnit() const
    { return d_unit; }

    //! Get the name of this field.
    const std::string& getName() const
    { return d_name; }

    //! Get the degrees of freedom tag on the mesh this field is associated
    //! with.
    iBase_TagHandle getDFTag() const
    { return d_dof_tag; }
}; 

} // end namespace FOOD

#include "TensorField_Def.hpp"

#endif // end FOOD_TENSORFIELD_HPP

//---------------------------------------------------------------------------//
// end TensorField.hpp
//---------------------------------------------------------------------------//

