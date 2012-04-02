//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file IntrepidQuadrature_Def.hpp
 * \author Stuart Slattery
 * \brief Quadrature definition for Intrepid cubature rule implementations.
 */
//---------------------------------------------------------------------------//

#ifndef FOOD_INTREPIDQUADRATURE_DEF_HPP
#define FOOD_INTREPIDQUADRATURE_DEF_HPP

#include "CellTopologyFactory.hpp"

#include <Teuchos_Tuple.hpp>

#include <Shards_CellTopology.hpp>

#include <Intrepid_CellTools.hpp>
#include <Intrepid_FunctionSpaceTools.hpp>

namespace FOOD
{
/*!
 * \brief Constructor.
 */
template<class Scalar>
IntrepidQuadrature<Scalar>::IntrepidQuadrature( 
    RCP_IntrepidCubature intrepid_cubature,
    const int entity_type, 
    const int entity_topology )
    : d_intrepid_cubature( intrepid_cubature )
{
    this->b_num_points = d_intrepid_cubature->getNumPoints();
    this->b_dimension = d_intrepid_cubature->getDimension();
    this->b_type = entity_type;
    this->b_topology = entity_topology;
}

/*!
v * \brief Destructor.
 */
template<class Scalar>
IntrepidQuadrature<Scalar>::~IntrepidQuadrature()
{ /* ... */ }

/*!
 * \brief Get the quadrature rule.
 */
template<class Scalar>
void IntrepidQuadrature<Scalar>::getQuadratureRule( 
    Teuchos::ArrayRCP<Scalar> &coordinates,
    Teuchos::ArrayRCP<Scalar> &weights ) const
{
    MDArray points( this->b_num_points, this->b_dimension );
    MDArray point_weights( this->b_num_points );

    d_intrepid_cubature->getCubature( points, point_weights );

    coordinates = points.getData();
    weights = point_weights.getData();
}

/*!
 * \brief Integrate over a cell.
 */
template<class Scalar>
void IntrepidQuadrature<Scalar>::integrate(
    Teuchos::ArrayRCP<Scalar> &integrated_values,
    const Teuchos::ArrayRCP<Scalar> &values,
    const int cardinality,
    const iMesh_Instance mesh,
    const iBase_EntityHandle physical_cell )
{
    // get cubature
    MDArray cubature_points( this->b_num_points, this->b_dimension );
    MDArray cubature_weights( this->b_num_points );
    d_intrepid_cubature->getCubature( cubature_points, cubature_weights );

    // compute jacobian
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

    CellTopologyFactory topo_factory;
    Teuchos::RCP<shards::CellTopology> cell_topo = 
	topo_factory.create( this->b_topology, element_nodes_size );

    MDArray jacobian( 
	1, this->b_num_points, this->b_dimension, this->b_dimension );
    Intrepid::CellTools<Scalar>::setJacobian( 
	jacobian, cubature_points, cell_nodes, *cell_topo );
	
    MDArray jacobian_det( 1, this->b_num_points );
    Intrepid::CellTools<Scalar>::setJacobianDet( jacobian_det, jacobian );

    MDArray weighted_measure( 1, this->b_num_points );
    Intrepid::FunctionSpaceTools::computeCellMeasure<Scalar>( 
	weighted_measure, jacobian_det, cubature_weights );

    Teuchos::Tuple<int,3> function_dimensions;
    function_dimensions[0] = 1;
    function_dimensions[1] = cardinality;
    function_dimensions[2] = this->b_num_points;
    MDArray function( function_dimensions, values );

    MDArray weighted_function( 
	1, cardinality, this->b_num_points );
    Intrepid::FunctionSpaceTools::multiplyMeasure<Scalar>( 
	weighted_function, weighted_measure, function );
   
    MDArray integrated_function( 1 );
    Intrepid::FunctionSpaceTools::integrate<Scalar>( 
	integrated_function, function, weighted_function, Intrepid::COMP_CPP );

    integrated_values = integrated_function.getData();
}

} // end namespace FOOD

#endif // end FOOD_INTREPIDQUADRATURE_DEF_HPP

//---------------------------------------------------------------------------//
// end IntrepidQuadrature_Def.hpp
//---------------------------------------------------------------------------//

