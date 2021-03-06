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

#include "Exception.hpp"
#include "TopologyTools.hpp"

#include <Teuchos_Tuple.hpp>

#include <Intrepid_FunctionSpaceTools.hpp>
#include <Intrepid_CellTools.hpp>

namespace FOOD
{

/*!
 * \brief Constructor.
 */
template<class Scalar>
IntrepidKernel<Scalar>::IntrepidKernel( RCP_IntrepidBasis intrepid_basis,
					const int entity_type, 
					const int entity_topology, 
					const int discretization_type,
					const int function_space_type )
    : d_intrepid_basis( intrepid_basis )
{
    this->b_cardinality = d_intrepid_basis->getCardinality();
    this->b_degree = d_intrepid_basis->getDegree();
    this->b_type = entity_type;
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
void IntrepidKernel<Scalar>::dfuncValue( Teuchos::ArrayRCP<Scalar> &values, 
					 const double param_coords[3] )
{
    MDArray coords(1,3);
    coords(0,0) = param_coords[0];
    coords(0,1) = param_coords[1];
    coords(0,2) = param_coords[2];

    if ( this->b_function_space_type == FOOD_HGRAD )
    {
	MDArray grad_values( this->b_cardinality, 1 );

	d_intrepid_basis->getValues( grad_values, 
				     coords, 
				     Intrepid::OPERATOR_VALUE );

	values = grad_values.getData();
    }
    else if ( this->b_function_space_type == FOOD_HDIV )
    {
	MDArray div_values( this->b_cardinality, 1, 3 );

	d_intrepid_basis->getValues( div_values, 
				     coords, 
				     Intrepid::OPERATOR_VALUE );

	values = div_values.getData();
    }
    else if ( this->b_function_space_type == FOOD_HCURL )
    {
	MDArray curl_values( this->b_cardinality, 1, 3 );

	d_intrepid_basis->getValues( curl_values, 
				     coords, 
				     Intrepid::OPERATOR_VALUE );

	values = curl_values.getData();
    }
    else
    {	    
	testPrecondition( this->b_function_space_type == FOOD_HGRAD ||
			  this->b_function_space_type == FOOD_HDIV  ||
			  this->b_function_space_type == FOOD_HCURL,
			  "Invalid function space type" );
    }
}

/*!
 * \brief Evaluate the operator value of a distribution function kernel at a
 * given set of parametric coordinates. The operator is defined by the
 * function space. (i.e. FOOD_HDIV space returns the divergence of the
 * distribution function kernel.) 
 */
template<class Scalar>
void IntrepidKernel<Scalar>::dfuncOperator( Teuchos::ArrayRCP<Scalar> &values, 
					    const double param_coords[3] )
{
    MDArray coords(1,3);
    coords(0,0) = param_coords[0];
    coords(0,1) = param_coords[1];
    coords(0,2) = param_coords[2];

    if ( this->b_function_space_type == FOOD_HGRAD )
    {
	MDArray dfunc_grad( this->b_cardinality, 1, 3 );
    
	d_intrepid_basis->getValues( dfunc_grad, 
				     coords, 
				     Intrepid::OPERATOR_GRAD );

	values = dfunc_grad.getData();
    }
    else if ( this->b_function_space_type == FOOD_HDIV )
    {
	MDArray dfunc_div( this->b_cardinality );

	d_intrepid_basis->getValues( dfunc_div, 
				     coords, 
				     Intrepid::OPERATOR_DIV );
	    
	values = dfunc_div.getData();
    }
    else if ( this->b_function_space_type == FOOD_HCURL )
    {
	MDArray dfunc_curl( this->b_cardinality, 1, 3 );

	d_intrepid_basis->getValues( dfunc_curl,
				     coords, 
				     Intrepid::OPERATOR_CURL );

	values = dfunc_curl.getData();
    }
    else
    {
	testPrecondition( this->b_function_space_type == FOOD_HGRAD ||
			  this->b_function_space_type == FOOD_HDIV  ||
			  this->b_function_space_type == FOOD_HCURL,
			  "Invalid function space type" );
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

    Teuchos::Tuple<int,3> cell_node_dimensions;
    cell_node_dimensions[0] = 1;
    cell_node_dimensions[1] = element_nodes_size;
    cell_node_dimensions[2] = 3;
    MDArray cell_nodes( Teuchos::Array<int>(cell_node_dimensions), 
			coord_array );

    MDArray reference_point( 1, 3 );
    MDArray coords( 1, 3 );
    coords(0,0) = physical_coords[0];
    coords(0,1) = physical_coords[1];
    coords(0,2) = physical_coords[2];
 
    Intrepid::CellTools<double>::mapToReferenceFrame( 
	reference_point,
	coords,
	cell_nodes,
	d_intrepid_basis->getBaseCellTopology(),
	0 );

    param_coords[0] = reference_point(0,0);
    param_coords[1] = reference_point(0,1);
    param_coords[2] = reference_point(0,2);

    free( element_nodes );
    free( coord_array );
}

/*!
 * \brief Transform the value of the distribution function kernel at a
 * given set of parametric coordinates back to the physical frame for the
 * given physical cell. 
 */
template<class Scalar>
void IntrepidKernel<Scalar>::transformValue( 
    Teuchos::ArrayRCP<Scalar> &transformed_values, 
    const Teuchos::ArrayRCP<Scalar> &values,
    const double param_coords[3],
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

    Teuchos::Tuple<int,3> cell_node_dimensions;
    cell_node_dimensions[0] = 1;
    cell_node_dimensions[1] = element_nodes_size;
    cell_node_dimensions[2] = 3;
    MDArray cell_nodes( Teuchos::Array<int>(cell_node_dimensions), 
			coord_array );

    MDArray jacobian( 1, 1, 3, 3 );
    MDArray jacobian_det( 1, 1 );
    MDArray jacobian_inv( 1, 1, 3, 3 );
    MDArray reference_coords( 1, 3 );
    reference_coords(0,0) = param_coords[0];
    reference_coords(0,1) = param_coords[1];
    reference_coords(0,2) = param_coords[2];

    if ( this->b_function_space_type == FOOD_HGRAD )
    {
	Teuchos::Tuple<int,3> transformed_grad_dimensions;
	transformed_grad_dimensions[0] = 1;
	transformed_grad_dimensions[1] = this->b_cardinality;
	transformed_grad_dimensions[2] = 1;
	MDArray transformed_grad( transformed_grad_dimensions,
				  transformed_values );

	Teuchos::Tuple<int,2> grad_value_dimensions;
	grad_value_dimensions[0] = this->b_cardinality;
	grad_value_dimensions[1] = 1;

	MDArray basis_grad( grad_value_dimensions, values );

	Intrepid::FunctionSpaceTools::HGRADtransformVALUE<Scalar,MDArray,MDArray>( 
	    transformed_grad, basis_grad );
    }
    else if ( this->b_function_space_type == FOOD_HDIV )
    {
	Intrepid::CellTools<double>::setJacobian( 
	    jacobian, 
	    reference_coords,
	    cell_nodes,
	    d_intrepid_basis->getBaseCellTopology() );

	Intrepid::CellTools<double>::setJacobianDet( jacobian_det, jacobian );

	Teuchos::Tuple<int,4> transformed_div_dimensions;
	transformed_div_dimensions[0] = 1;
	transformed_div_dimensions[1] = this->b_cardinality;
	transformed_div_dimensions[2] = 1;
	transformed_div_dimensions[3] = 3;
	MDArray transformed_div( transformed_div_dimensions, 
				 transformed_values );

	Teuchos::Tuple<int,3> div_value_dimensions;
	div_value_dimensions[0] = this->b_cardinality;
	div_value_dimensions[1] = 1;
	div_value_dimensions[2] = 3;
	MDArray basis_div( div_value_dimensions, values );
	Intrepid::FunctionSpaceTools::HDIVtransformVALUE<Scalar,MDArray,MDArray>( 
	    transformed_div, jacobian, jacobian_det, basis_div );
    }
    else if ( this->b_function_space_type == FOOD_HCURL )
    {
	Intrepid::CellTools<double>::setJacobian( 
	    jacobian, 
	    reference_coords,
	    cell_nodes,
	    d_intrepid_basis->getBaseCellTopology() );

	Intrepid::CellTools<double>::setJacobianInv( jacobian_inv, jacobian );

	Teuchos::Tuple<int,4> transformed_curl_dimensions;
	transformed_curl_dimensions[0] = 1;
	transformed_curl_dimensions[1] = this->b_cardinality;
	transformed_curl_dimensions[2] = 1;
	transformed_curl_dimensions[3] = 3;
	MDArray transformed_curl( transformed_curl_dimensions,
				  transformed_values );

	Teuchos::Tuple<int,3> curl_value_dimensions;
	curl_value_dimensions[0] = this->b_cardinality;
	curl_value_dimensions[1] = 1;
	curl_value_dimensions[2] = 3;
	MDArray basis_curl( curl_value_dimensions, values );
	Intrepid::FunctionSpaceTools::HCURLtransformVALUE<Scalar,MDArray,MDArray>( 
	    transformed_curl, jacobian_inv, basis_curl );
    }
    else
    {	    
	testPrecondition( this->b_function_space_type == FOOD_HGRAD ||
			  this->b_function_space_type == FOOD_HDIV  ||
			  this->b_function_space_type == FOOD_HCURL,
			  "Invalid function space type" );
    }

    free( coord_array );
    free( element_nodes );
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
    Teuchos::ArrayRCP<Scalar> &transformed_values,
    const Teuchos::ArrayRCP<Scalar> &values,
    const double param_coords[3],
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

    Teuchos::Tuple<int,3> cell_node_dimensions;
    cell_node_dimensions[0] = 1;
    cell_node_dimensions[1] = element_nodes_size;
    cell_node_dimensions[2] = 3;
    MDArray cell_nodes( Teuchos::Array<int>(cell_node_dimensions), 
			coord_array );

    MDArray jacobian( 1, 1, 3, 3 );
    MDArray jacobian_det( 1, 1 );
    MDArray jacobian_inv( 1, 1, 3, 3 );
    MDArray reference_point( 1, 3 );
    reference_point(0,0) = param_coords[0];
    reference_point(0,1) = param_coords[1];
    reference_point(0,2) = param_coords[2];

    if ( this->b_function_space_type == FOOD_HGRAD )
    {
	Intrepid::CellTools<double>::setJacobian( 
	    jacobian, 
	    reference_point,
	    cell_nodes,
	    d_intrepid_basis->getBaseCellTopology() );

	Intrepid::CellTools<double>::setJacobianInv( jacobian_inv, 
						     jacobian );

	MDArray transformed_grad( 1, this->b_cardinality, 1, 3 );
	Teuchos::Tuple<int,3> grad_dimensions;
	grad_dimensions[0] = this->b_cardinality;
	grad_dimensions[1] = 1;
	grad_dimensions[2] = 3;
	MDArray basis_grad( grad_dimensions, values );
	Intrepid::FunctionSpaceTools::HGRADtransformGRAD<double>( 
	    transformed_grad, jacobian_inv, basis_grad );

	transformed_values = transformed_grad.getData();
    }
    else if ( this->b_function_space_type == FOOD_HDIV )
    {

    }
    else if ( this->b_function_space_type == FOOD_HCURL )
    {

    }
    else
    {
	testPrecondition( this->b_function_space_type == FOOD_HGRAD ||
			  this->b_function_space_type == FOOD_HDIV  ||
			  this->b_function_space_type == FOOD_HCURL,
			  "Invalid function space type" );
    }

    free( coord_array );
    free( element_nodes );
}

/*!
 * \brief Evaluate the distribution function using function coefficients
 * and physical frame distribution function kernel values.
 */
template<class Scalar>
void IntrepidKernel<Scalar>::evaluate( 
    Teuchos::ArrayRCP<Scalar> &function_values,
    const Teuchos::ArrayRCP<Scalar> &coeffs,
    const Teuchos::ArrayRCP<Scalar> &dfunc_values )
{
    int dim1 = this->b_cardinality;
    testPrecondition( dim1 == (int) coeffs.size(),
	"Function coefficients size does not match basis cardinality" );
    
    MDArray function_coeffs( 1, dim1 );
    for ( int m = 0; m < dim1; ++m )
    {
	function_coeffs(0,m) = coeffs[m];
    }
    MDArray basis_eval( 1, dim1, 1 );
    for ( int i = 0; i < dim1; ++i )
    {
	basis_eval( 0, i, 0 ) = dfunc_values[i];
    }

    Teuchos::Tuple<int,2> function_dimensions;
    function_dimensions[0] = 1;
    function_dimensions[1] = 1;
    MDArray function_eval( function_dimensions, function_values );
    Intrepid::FunctionSpaceTools::evaluate<Scalar>( function_eval,
						    function_coeffs, 
						    basis_eval );
}

/*!
 * \brief For a given point in the physical frame of a cell, compute the
 * determinant of the Jacobian.
 */
template<class Scalar>
void IntrepidKernel<Scalar>::jacobianDet( 
    Scalar &determinant,
    const double param_coords[3],
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

    Teuchos::Tuple<int,3> cell_node_dimensions;
    cell_node_dimensions[0] = 1;
    cell_node_dimensions[1] = element_nodes_size;
    cell_node_dimensions[2] = 3;
    MDArray cell_nodes( Teuchos::Array<int>(cell_node_dimensions), 
			coord_array );

    MDArray reference_coords( 1, 3 );
    reference_coords(0,0) = param_coords[0];
    reference_coords(0,1) = param_coords[1];
    reference_coords(0,2) = param_coords[2];

    MDArray jacobian( 1, 1, 3, 3 );
    Intrepid::CellTools<double>::setJacobian( 
	jacobian, 
	reference_coords,
	cell_nodes,
	d_intrepid_basis->getBaseCellTopology() );

    MDArray jacobian_det( 1, 1 );
    Intrepid::CellTools<double>::setJacobianDet( jacobian_det, jacobian );
    determinant = jacobian_det(0,0);

    free( element_nodes );
    free( coord_array );
}

} // end namespace FOOD

#endif // end FOOD_INTREPIDKERNEL_DEF_HPP

//---------------------------------------------------------------------------//
// end IntrepidKernel_Def.hpp
//---------------------------------------------------------------------------//

