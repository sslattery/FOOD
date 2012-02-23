//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \field IntrepidKernel_Def.hpp
 * \author Stuart Slattery
 * \brief Distribution function kernel definition for Intrepid basis
 * implemenentations.
 */
//---------------------------------------------------------------------------//

#ifndef FOOD_INTREPIDKERNEL_DEF_HPP
#define FOOD_INTREPIDKERNEL_DEF_HPP

#include <cassert>

#include "BasisFactory.hpp"

#include <Teuchos_Tuple.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include <Intrepid_FunctionSpaceTools.hpp>

namespace FOOD
{

/*!
 * \brief Constructor.
 */
template<class Scalar>
IntrepidKernel<Scalar>::IntrepidKernel( RCP_IntrepidBasis intrepid_basis,
					const int entity_topology, 
					const int discretization_type,
					const int function_space_type )
    : d_intrepid_basis( intrepid_basis )
{
    this->b_cardinality = d_intrepid_basis->getCardinality();
    this->b_degree = d_intrepid_basis->getDegree();
    this->b_topology = entity_topology;
    this->b_cn_type = FOOD_SHARDSCN;
    this->b_coord_type = FOOD_CARTESIAN;
    this->b_discretization_type = discretization_type;
    this->b_function_space_type = function_space_type;
}

/*!
 * \brief Destructor.
 */
template<class Scalar>
IntrepidKernel<Scalar>::~IntrepidKernel()
{ /* ... */ }

/*!
 * \brief Evaluate the value of the distribution function kernel at a given
 * set of parametric coordinates.  
 */
template<class Scalar>
void IntrepidKernel<Scalar>::dfuncValue( Scalar **values, 
					 int *num_values,
					 const double param_coords[3] )
{
    num_values = this->b_cardinality;

    Teuchos::Tuple<int,2> value_dimensions = { num_values, 1 };
    MDArray dfunc_values( value_dimensions );

    MDArray coords(1,3);
    coords(0,0) = param_coords[0];
    coords(1,0) = param_coords[1];
    coords(2,0) = param_coords[2];

    d_basis->getValues( dfunc_values, 
			coords, 
			Intrepid::OPERATOR_VALUE );

    values = dfunc_values->getData()->get();
}

/*!
 * \brief Evaluate the operator value of a distribution function kernel at a
 * given set of parametric coordinates. The operator is defined by the
 * function space. (i.e. FOOD_HDIV space returns the divergence of the
 * distribution function kernel.) 
 */
template<class Scalar>
void IntrepidKernel<Scalar>::dfuncOperator( Scalar **values, 
					    int *num_values,
					    const double param_coords[3] )
{
    switch( this->b_function_space_type )
    {
	case FOOD_HGRAD:

	    num_values = 3*this->b_cardinality;

	    Teuchos::Tuple<int,3> grad_dimensions = { this->b_cardinality, 1, 3 };
	    MDArray dfunc_grad( grad_dimensions );

	    MDArray coords(1,3);
	    coords(0,0) = param_coords[0];
	    coords(1,0) = param_coords[1];
	    coords(2,0) = param_coords[2];
    
	    d_intrepid_basis->getValues( dfunc_grad, 
					 coords, 
					 Intrepid::OPERATOR_GRAD );

	    values = dfunc_grad->getData()->get();
	    break;

	case FOOD_HDIV:

	    num_values = 3*this->b_cardinality;

	    Teuchos::Tuple<int,3> div_dimensions = { this->b_cardinality, 1, 3 };
	    MDArray dfunc_div( div_dimensions );

	    MDArray coords(1,3);
	    coords(0,0) = param_coords[0];
	    coords(1,0) = param_coords[1];
	    coords(2,0) = param_coords[2];

	    d_intrepid_basis->getValues( dfunc_div, 
					 coords, 
					 Intrepid::OPERATOR_DIV );
	    break;

	case FOOD_HCURL:

	    num_values = 3*this->b_cardinality;

	    Teuchos::Tuple<int,3> curl_dimensions = { this->b_cardinality, 1, 3 };
	    MDArray dfunc_curl( curl_dimensions );

	    MDArray coords(1,3);
	    coords(0,0) = param_coords[0];
	    coords(1,0) = param_coords[1];
	    coords(2,0) = param_coords[2];

	    d_intrepid_basis->getValues( dfunc_curl_values, 
					 coords, 
					 Intrepid::OPERATOR_CURL );
	    break;
    }
}

/*!
 * \brief Transform a point from the physical frame in a physical cell to
 * the reference frame of the reference cell for the distribution function
 * kernel topology. 
 */
template<class Scalar>
void IntrepidKernel<Scalar>::transformPoint( 
    double param_coords[3],
    const double physical_coords[3],
    const iMesh_Instance mesh,
    const iBase_EntityHandle physical_cell )
{
    int error = 0;

    iBase_EntityHandle *element_nodes = 0;
    int element_nodes_allocated = 0;
    int element_nodes_size = 0;
    iMesh_getEntAdj( mesh,
		     physical_cell,
		     iBase_VERTEX,
		     &element_nodes,
		     &element_nodes_allocated,
		     &element_nodes_size,
		     &error );
    assert( iBase_SUCCESS == error );

    TopologyTools::MBCN2Shards( element_nodes, 
				element_nodes_size,
				this->b_topology );

    int coords_allocated = 0;
    int coords_size = 0;
    double *coord_array = 0;
    iMesh_getVtxArrCoords( mesh,
			   element_nodes,
			   element_nodes_size,
			   iBase_INTERLEAVED,
			   &coord_array,
			   &coords_allocated,
			   &coords_size,
			   &error );
    assert( iBase_SUCCESS == error );

    Teuchos::Tuple<int,3> cell_node_dimensions = { 1, element_nodes_size, 3 };
    MDArray cell_nodes( Teuchos::Array<int>(cell_node_dimensions), 
			coord_array );

    free( element_nodes );
    free( coord_array );

    MDArray reference_point( 1, 3 );
    MDArray coords( 1, 3 );
    coords(0,0) = physical_coords[0];
    coords(0,1) = physical_coords[1];
    coords(0,2) = physical_coords[2];
 
    Intrepid::CellTools<double>::mapToReferenceFrame( 
	reference_point,
	coords,
	cell_nodes,
	d_intrepid_basis->getBaseTopology(),
	0 );

    param_coords[0] = reference_point(0,0);
    param_coords[1] = reference_point(0,1);
    param_coords[2] = reference_point(0,2);
}

/*!
 * \brief Transform the value of the distribution function kernel at a
 * given set of parametric coordinates back to the physical frame for the
 * given physical cell. 
 */
template<class Scalar>
void IntrepidKernel<Scalar>::transformValue( 
    Scalar **transformed_values, 
    const Scalar *values,
    const int num_values,
    const double param_coords[3],
    const iBase_EntityHandle physical_cell )
{
    MDArray jacobian( 1, 1, 3, 3 );
    MDArray jacobian_det( 1, 1 );
    MDArray jacobian_inv( 1, 1, 3, 3 );

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
		d_intrepid_basis->getBaseCellTopology() );

	    Intrepid::CellTools<double>::setJacobianDet( jacobian_det, jacobian );

	    Intrepid::FunctionSpaceTools::HDIVtransformVALUE<Scalar,MDArray,MDArray>( 
		transformed_eval, jacobian, jacobian_det, basis_eval );

	    break;

	case FOOD_HCURL:

	    Intrepid::CellTools<double>::setJacobian( 
		jacobian, 
		reference_coords,
		cell_nodes,
		d_intrepid_basis->getBaseCellTopology() );

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

/*!
 * \brief Transform the operator value of a distribution function kernel
 * at a given set of parametric coordinates back to the physical frame for
 * the given physical cell. The operator is defined by the function
 * space. (i.e. FOOD_HDIV space returns the divergence of the distribution
 * function kernel.)  
 */
template<class Scalar>
void IntrepidKernel<Scalar>::transformOperator( 
    Scalar **transformed_values,
    const Scalar *values,
    const int num_values,
    const double param_coords[3],
    const iBase_EntityHandle physical_cell )
{
    int error = 0;

    iBase_EntityHandle *element_nodes = 0;
    int element_nodes_allocated = 0;
    int element_nodes_size = 0;
    iMesh_getEntAdj( mesh,
		     physical_cell,
		     iBase_VERTEX,
		     &element_nodes,
		     &element_nodes_allocated,
		     &element_nodes_size,
		     &error );
    assert( iBase_SUCCESS == error );

    TopologyTools::MBCN2Shards( element_nodes, 
				element_nodes_size,
				this->b_topology );

    int coords_allocated = 0;
    int coords_size = 0;
    double *coord_array = 0;
    iMesh_getVtxArrCoords( mesh,
			   element_nodes,
			   element_nodes_size,
			   iBase_INTERLEAVED,
			   &coord_array,
			   &coords_allocated,
			   &coords_size,
			   &error );
    assert( iBase_SUCCESS == error );

    Teuchos::Tuple<int,3> cell_node_dimensions = { 1, element_nodes_size, 3 };
    MDArray cell_nodes( Teuchos::Array<int>(cell_node_dimensions), 
			coord_array );

    free( coord_array );
    free( element_nodes );

    MDArray jacobian( 1, 1, 3, 3 );
    MDArray jacobian_det( 1, 1 );
    MDArray jacobian_inv( 1, 1, 3, 3 );
    MDArray reference_point( 1, 3 );
    reference_point(0,0) = param_coords[0];
    reference_point(0,1) = param_coords[1];
    reference_point(0,2) = param_coords[2];

    switch ( d_basis_function_space )
    {
	case FOOD_HGRAD:

	    Intrepid::CellTools<double>::setJacobian( 
		jacobian, 
		reference_point,
		cell_nodes,
		d_intrepid_basis->getCellTopology() );

	    Intrepid::CellTools<double>::setJacobianInv( jacobian_inv, 
							 jacobian );

	    MDArray transformed_eval( 1, this->b_cardinality, 1, 3 );
	    num_values = 3*this->b_cardinality;
	    Intrepid::FunctionSpaceTools::HGRADtransformGRAD<double>( 
		transformed_eval, jacobian_inv, basis_eval );
	    break;

	case FOOD_HDIV:
	    
	    break;

	case FOOD_HCURL:
	    
	    break;
    }

    values = transformed_eval->getData()->get();
}

/*!
 * \brief Evaluate the distribution function using function coefficients
 * and physical frame distribution function kernel values.
 */
template<class Scalar>
void IntrepidKernel<Scalar>::evaluate( Scalar **function_values,
				       const Scalar* coeffs,
				       const int num_coeffs,
				       const Scalar* dfunc_values,
				       const int num_dfunc_values )
{
    int dim1 = this->b_cardinality;
    int dim2 = num_dfunc_values / dim1;

    assert( num_coeffs == dim1 );
    assert( num_dfunc_values == dim1 );

    MDArray component_values( 1, 1, 1 );
    MDArray component_coeffs( 1, dim1 );

    for ( int m = 0; m < dim1; ++m )
    {
	component_coeffs(0,m) = coeffs[m];
    }

    MDArray physical_eval( 1, dim1, 1, dim2 );
    Intrepid::FunctionSpaceTools::evaluate<Scalar>( component_values,
						    component_coeffs, 
						    physical_eval );

    for ( int p = 0; p < coords.dimension(0); ++p )
    {
	for ( int d = 0; d < coords.dimension(1); ++d )
	{
	    dfunc_values(0,p,n,d) = component_values(0,p,d);
	}
    }
}

} // end namespace FOOD

#endif // end FOOD_INTREPIDKERNEL_DEF_HPP

//---------------------------------------------------------------------------//
// end IntrepidKernel_Def.hpp
//---------------------------------------------------------------------------//

