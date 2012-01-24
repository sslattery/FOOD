//---------------------------------------------------------------------------//
// \field DFuncKernel_Def.hpp
// \author Stuart Slattery
// \brief Distribution function kernel definition.
//---------------------------------------------------------------------------//

#ifndef FOOD_DFUNCKERNEL_DEF_HPP
#define FOOD_DFUNCKERNEL_DEF_HPP

#include "BasisFactory.hpp"
#include "CellTopologyFactory.hpp"

#include <Shards_CellTopology.hpp>

#include <Intrepid_DefaultCubatureFactory.hpp>

namespace FOOD
{

/*!
 * \brief Constructor.
 */
template<class Scalar>
DFuncKernel<Scalar>::DFuncKernel( const int entity_topology,
				  const int discretization_type,
				  const int basis_operator_type,
				  const int basis_degree,
				  const int cubature_degree )
    : d_basis(0)
    , d_cubature(0)
{
    BasisFactory<Scalar,MDArray> basis_factory;
    d_basis = basis_factory.create( entity_topology,
				    discretization_type,
				    basis_operator_type,
				    basis_degree );

    CellTopologyFactory cell_topo_factory;
    Teuchos::RCP<shards::CellTopology> cell_topo = 
	cell_topo_factory.create( entity_topology );
    
    Intrepid::DefaultCubatureFactory<Scalar> cubature_factory;
    d_cubature = cubature_factory.create( *cell_topo, cubature_degree );
}

/*!
 * \brief Destructor.
 */
template<class Scalar>
DFuncKernel<Scalar>::~DFuncKernel()
{ /* ... */ }

/*!
 * \brief Evaluate the degrees of freedom for this kernel at a specific
 * location. Coordinates are parametric (defined in the reference cell).
 *
 * \param dfunc_values The evaluated distribution function
 * values. ArrayDim[basis_cardinality][num_points].
 * \param coords The coords to evaluate the distribution functions
 * at. ArrayDim[num_points][xyz].
 */
template<class Scalar>
void DFuncKernel<Scalar>::evaluateDF( MDArray &dfunc_values, 
				      const MDArray &coords )
{
    d_basis.getValues( dfunc_values, coords, Intrepid::OPERATOR_VALUE );
}

} // end namespace FOOD

#endif // end FOOD_DFUNCKERNEL_DEF_HPP

//---------------------------------------------------------------------------//
// end DFuncKernel_Def.hpp
//---------------------------------------------------------------------------//

