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
				      const int basis_operator_type,
    				      const int basis_degree )
    : d_basis(0)
{
    BasisFactory<ScalarType,MDArray> basis_factory;
    d_basis = basis_factory.create( entity_topology,
				    discretization_type,
				    basis_operator_type,
				    basis_degree );
}

/*!
 * \brief Destructor.
 */
template<class ScalarType>
DFuncKernel<ScalarType>::~DFuncKernel()
{ /* ... */ }

/*!
 * \brief Evaluate the degrees of freedom for this kernel at a specific
 * location. Coordinates are parametric. 
 */
void DFuncKernel<ScalarType>::evaluateDF( MDArray &dfunc_values, 
					  const MDArray &coords )
{

}


} // end namespace FOOD

#endif // end FOOD_DFUNCKERNEL_DEF_HPP

//---------------------------------------------------------------------------//
// end DFuncKernel_Def.hpp
//---------------------------------------------------------------------------//

