//---------------------------------------------------------------------------//
// \field DFuncKernel_Def.hpp
// \author Stuart Slattery
// \brief Distribution function kernel definition.
//---------------------------------------------------------------------------//

#ifndef FOOD_DFUNCKERNEL_DEF_HPP
#define FOOD_DFUNCKERNEL_DEF_HPP

#include "BasisFactory.hpp"

#include <Teuchos_ParameterList.hpp>

namespace FOOD
{

/*!
 * \brief Constructor.
 */
template<class ScalarType>
DFuncKernel<ScalarType>::DFuncKernel( const int entity_topology,
				      const int discretization_type,
				      const int basis_degree,
				      const int basis_operator_type )
    : d_basis(0)
{ /* ... */ }

/*!
 * \brief Destructor.
 */
template<class ScalarType>
DFuncKernel<ScalarType>::~DFuncKernel()
{ /* ... */ }

} // end namespace FOOD

#endif // end FOOD_DFUNCKERNEL_DEF_HPP

//---------------------------------------------------------------------------//
// end DFuncKernel_Def.hpp
//---------------------------------------------------------------------------//

