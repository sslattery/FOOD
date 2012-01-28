//---------------------------------------------------------------------------//
// \field DFuncKernel_Def.hpp
// \author Stuart Slattery
// \brief Distribution function kernel definition.
//---------------------------------------------------------------------------//

#ifndef FOOD_DFUNCKERNEL_DEF_HPP
#define FOOD_DFUNCKERNEL_DEF_HPP

#include "BasisFactory.hpp"
#include "CellTopologyFactory.hpp"

#include <Intrepid_FunctionSpaceTools.hpp>

namespace FOOD
{

/*!
 * \brief Constructor.
 */
template<class Scalar>
DFuncKernel<Scalar>::DFuncKernel( const int entity_topology,
				  const int discretization_type,
				  const int basis_operator_type,
				  const int basis_degree )
    : d_entity_topology(entity_topology)
    , d_discretization_type(discretization_type)
    , d_basis_operator_type(basis_operator_type)
    , d_basis(0)
    , d_cell_topology(0)
{
    BasisFactory<Scalar,MDArray> basis_factory;
    d_basis = basis_factory.create( entity_topology,
				    discretization_type,
				    basis_operator_type,
				    basis_degree );

    CellTopologyFactory topo_factory;
    d_cell_topology = topo_factory.create( entity_topology );
}

/*!
 * \brief Destructor.
 */
template<class Scalar>
DFuncKernel<Scalar>::~DFuncKernel()
{ /* ... */ }

/*!
 * \brief Evaluate the degrees of freedom for this kernel at a specific
 * location.
 *
 * \param dfunc_values The evaluated distribution function
 * values. ArrayDim[basis_carindality][num_points].
 * \param coords The coords to evaluate the distribution functions
 * at. ArrayDim[num_points][dimension].
 */
template<class Scalar>
void DFuncKernel<Scalar>::evaluateDF( MDArray &values_at_coords, 
				      const MDArray &coords )
{
    d_basis->getValues( values_at_coords, 
			coords, 
			Intrepid::OPERATOR_VALUE );
}

/*!
 * \brief Evaluate the gradient of the degrees of freedom for this kernel at a
 * specific location.
 *
 * \param dfunc_values The evaluated distribution function gradient.
 * values. ArrayDim[basis_carindality][num_points].
 * \param coords The coords to evaluate the distribution functions
 * at. ArrayDim[num_points][dimension].
 */
template<class Scalar>
void DFuncKernel<Scalar>::evaluateGradDF( MDArray &values_at_coords, 
					  const MDArray &coords )
{
    d_basis->getValues( values_at_coords, 
			coords, 
			Intrepid::OPERATOR_GRAD );
}

} // end namespace FOOD

#endif // end FOOD_DFUNCKERNEL_DEF_HPP

//---------------------------------------------------------------------------//
// end DFuncKernel_Def.hpp
//---------------------------------------------------------------------------//

