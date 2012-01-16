//---------------------------------------------------------------------------//
// \file TensorField_Def.hpp
// \author Stuart Slattery
// \brief Tensor field definition.
//---------------------------------------------------------------------------//

#ifndef FOOD_TENSORFIELD_DEF_HPP
#define FOOD_TENSORFIELD_DEF_HPP

namespace FOOD
{

/*!
 * \brief Constructor.
 */
template<class ScalarType>
TensorField<ScalarType>::TensorField(RCP_Domain domain,
				     int coord_type,
				     RCP_TensorTemplate tensor_template,
				     RCP_Unit unit,
				     const std::string &name )
    : d_domain(domain)
    , d_coord_type(coord_type)
    , d_tensor_template(tensor_template)
    , d_unit(unit)
    , d_name(name)
{ /* ... */ }

/*!
 * \brief Destructor.
 */
template<class ScalarType>
TensorField<ScalarType>::~TensorField()
{ /* ... */ }

}

#endif // end FOOD_TENSORFIELD_DEF_HPP

//---------------------------------------------------------------------------//
// end TensorField_Def.hpp
//---------------------------------------------------------------------------//

