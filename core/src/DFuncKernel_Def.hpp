//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \field DFuncKernel_Def.hpp
 * \author Stuart Slattery
 * \brief Distribution function kernel definition.
 */
//---------------------------------------------------------------------------//

#ifndef FOOD_DFUNCKERNEL_DEF_HPP
#define FOOD_DFUNCKERNEL_DEF_HPP

#include <cassert>

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
    d_cell_topology = topo_factory.create( eval_topology, d_basis->getCardinality() );
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
void DFuncKernel<Scalar>::evaluateValueBasis( MDArray &dfunc_values, 
					      const MDArray &coords )
{
    d_basis->getValues( dfunc_values, 
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
void DFuncKernel<Scalar>::evaluateGradBasis( MDArray &dfunc_grad_values, 
					     const MDArray &coords )
{
    assert( FOOD_HGRAD == d_basis_function_space );
    d_basis->getValues( dfunc_grad_values, 
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
void DFuncKernel<Scalar>::evaluateDivBasis( MDArray &dfunc_div_values, 
					    const MDArray &coords )
{
    assert( FOOD_HDIV == d_basis_function_space );
    d_basis->getValues( dfunc_div_values, 
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
void DFuncKernel<Scalar>::evaluateCurlBasis( MDArray &dfunc_curl_values, 
					     const MDArray &coords )
{
    assert( FOOD_HCURL == d_basis_function_space );
    d_basis->getValues( dfunc_curl_values, 
			coords, 
			Intrepid::OPERATOR_CURL );
}

/*!
 * \brief Transform evaluated basis values to physical frame.
 */
template<class Scalar>
void DFuncKernel<Scalar>::transformValue( MDArray &transformed_eval,
					  const MDArray &reference_coords,
					  const MDArray &cell_nodes,
					  const MDArray &basis_eval )
{
    MDArray jacobian( 1, 
		      reference_coords.dimension(0),
		      reference_coords.dimension(1),
		      reference_coords.dimension(1) );

    MDArray jacobian_det( 1, reference_coords.dimension(0) );

    MDArray jacobian_inv( 1, 
			  reference_coords.dimension(0),
			  reference_coords.dimension(1),
			  reference_coords.dimension(1) );

    switch ( d_basis_function_space )
    {
	case FOOD_HGRAD:

	    Intrepid::FunctionSpaceTools::HGRADtransformVALUE<Scalar,MDArray,MDArray>( 
		transformed_eval, basis_eval );

	    break;

	case FOOD_HDIV:

	    Intrepid::CellTools<double>::setJacobian( 
		jacobian, 
		reference_coords,
		cell_nodes,
		*d_cell_topology );

	    Intrepid::CellTools<double>::setJacobianDet( jacobian_det, jacobian );

	    Intrepid::FunctionSpaceTools::HDIVtransformVALUE<Scalar,MDArray,MDArray>( 
		transformed_eval, jacobian, jacobian_det, basis_eval );

	    break;

	case FOOD_HCURL:

	    Intrepid::CellTools<double>::setJacobian( 
		jacobian, 
		reference_coords,
		cell_nodes,
		*d_cell_topology );

	    Intrepid::CellTools<double>::setJacobianInv( jacobian_inv, jacobian );

	    Intrepid::FunctionSpaceTools::HCURLtransformVALUE<Scalar,MDArray,MDArray>( 
		transformed_eval, jacobian_inv, basis_eval );

	    break;

	default:
	    
	    assert( d_basis_function_space == FOOD_HGRAD ||
		    d_basis_function_space == FOOD_HDIV  ||
		    d_basis_function_space == FOOD_HCURL );
    }
}

} // end namespace FOOD

#endif // end FOOD_DFUNCKERNEL_DEF_HPP

//---------------------------------------------------------------------------//
// end DFuncKernel_Def.hpp
//---------------------------------------------------------------------------//

