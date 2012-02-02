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
DFuncKernel<Scalar>::DFuncKernel( const int eval_type,
				  const int eval_topology,
				  const int dof_entity_type,
				  const int dof_entity_topology,
				  const int coordinate_type,
				  const int discretization_type,
				  const int basis_function_space,
				  const int basis_degree )
    : d_eval_type(eval_type)
    , d_eval_topology(eval_topology)
    , d_dof_entity_type(dof_entity_type)
    , d_dof_entity_topology(dof_entity_topology)
    , d_coordinate_type(coordinate_type)
    , d_discretization_type(discretization_type)
    , d_basis_function_space(basis_function_space)
    , d_basis(0)
    , d_cell_topology(0)
{
    BasisFactory<Scalar,MDArray> basis_factory;
    d_basis = basis_factory.create( eval_topology,
				    discretization_type,
				    basis_function_space,
				    basis_degree );

    CellTopologyFactory topo_factory;
    d_cell_topology = topo_factory.create( eval_topology );
}

/*!
 * \brief Destructor.
 */
template<class Scalar>
DFuncKernel<Scalar>::~DFuncKernel()
{ /* ... */ }

/*!
 * \brief Evaluate the basis for this kernel at a specific
 * location.
 *
 * \param dfunc_values The evaluated distribution function
 * values. ArrayDim[basis_carindality][num_points].
 * \param coords The coords to evaluate the distribution functions
 * at. ArrayDim[num_points][dimension].
 */
template<class Scalar>
void DFuncKernel<Scalar>::evaluateValueBasis( MDArray &values_at_coords, 
					      const MDArray &coords )
{
    d_basis->getValues( values_at_coords, 
			coords, 
			Intrepid::OPERATOR_VALUE );
}

/*!
 * \brief Evaluate the gradient of the basis for this kernel at a
 * specific location.
 *
 * \param dfunc_values The evaluated distribution function gradient.
 * values. ArrayDim[basis_carindality][num_points].
 * \param coords The coords to evaluate the distribution functions
 * at. ArrayDim[num_points][dimension].
 */
template<class Scalar>
void DFuncKernel<Scalar>::evaluateGradBasis( MDArray &values_at_coords, 
					     const MDArray &coords )
{
    assert( FOOD_HGRAD == d_basis_function_space );
    d_basis->getValues( values_at_coords, 
			coords, 
			Intrepid::OPERATOR_GRAD );
}

/*!
 * \brief Evaluate the divergence of the basis for this kernel at a
 * specific location.
 *
 * \param dfunc_values The evaluated distribution function divergence.
 * values. ArrayDim[basis_carindality][num_points].
 * \param coords The coords to evaluate the distribution functions
 * at. ArrayDim[num_points][dimension].
 */
template<class Scalar>
void DFuncKernel<Scalar>::evaluateDivBasis( MDArray &values_at_coords, 
					    const MDArray &coords )
{
    assert( FOOD_HDIV == d_basis_function_space );
    d_basis->getValues( values_at_coords, 
			coords, 
			Intrepid::OPERATOR_DIV );
}

/*!
 * \brief Evaluate the curl of the basis for this kernel at a
 * specific location.
 *
 * \param dfunc_values The evaluated distribution function curl.
 * values. ArrayDim[basis_carindality][num_points].
 * \param coords The coords to evaluate the distribution functions
 * at. ArrayDim[num_points][dimension].
 */
template<class Scalar>
void DFuncKernel<Scalar>::evaluateCurlBasis( MDArray &values_at_coords, 
					     const MDArray &coords )
{
    assert( FOOD_HCURL == d_basis_function_space );
    d_basis->getValues( values_at_coords, 
			coords, 
			Intrepid::OPERATOR_CURL );
}

/*!
 * \brief Transform evaluated basis values to physical frame.
 */
template<class Scalar>
void DFuncKernel<Scalar>::transformValue( MDArray &transformed_eval,
					  const MDArray &basis_eval )
{
    if ( d_basis_function_space == FOOD_HGRAD )
    {
	Intrepid::FunctionSpaceTools::HGRADtransformVALUE<Scalar,MDArray,MDArray>( 
	    transformed_eval, basis_eval );
    }
    else if ( d_basis_function_space == FOOD_HDIV )
    {
	// Intrepid::FunctionSpaceTools::HDIVtransformVALUE<Scalar,MDArray,MDArray>( 
	//     transformed_eval, basis_eval );
    }
    else if ( d_basis_function_space == FOOD_HCURL )
    {
	// Intrepid::FunctionSpaceTools::HCURLtransformVALUE<Scalar,MDArray,MDArray>( 
	// transformed_eval, basis_eval );
    }
}

/*!
 * \brief Transform evaluated basis operator values to physical frame.
 */
template<class Scalar>
void DFuncKernel<Scalar>::transformOperator( MDArray &transformed_eval,
					     const MDArray &basis_eval )
{
    // if ( d_basis_function_space == FOOD_HGRAD )
    // {
    // 	Intrepid::FunctionSpaceTools::HGRADtransformGRAD<Scalar,MDArray,MDArray>( 
    // 	    transformed_eval, basis_eval );
    // }
    // else if ( d_basis_function_space == FOOD_HDIV )
    // {
    // 	Intrepid::FunctionSpaceTools::HDIVtransformDIV<Scalar,MDArray,MDArray>( 
    // 	    transformed_eval, basis_eval );
    // }
    // else if ( d_basis_function_space == FOOD_HCURL )
    // {
    // 	Intrepid::FunctionSpaceTools::HCURLtransformCURL<Scalar,MDArray,MDArray>( 
    // 	    transformed_eval, basis_eval );
    // }
}

} // end namespace FOOD

#endif // end FOOD_DFUNCKERNEL_DEF_HPP

//---------------------------------------------------------------------------//
// end DFuncKernel_Def.hpp
//---------------------------------------------------------------------------//

