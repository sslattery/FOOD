//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file TopologyTools.cpp
 * \author Stuart Slattery
 * \brief TopologyTools definition
 */
//---------------------------------------------------------------------------//

#include <vector>
#include <cassert>

#include "Exception.hpp"
#include "CellTopologyFactory.hpp"
#include "TopologyTools.hpp"

#include <Teuchos_ENull.hpp>
#include <Teuchos_RCP.hpp>

#include <Shards_CellTopology.hpp>

#include <Intrepid_FieldContainer.hpp>
#include <Intrepid_CellTools.hpp>

namespace FOOD
{

/*! 
 * \brief Get the number of linear nodes for a particular iMesh topology.
 */
int TopologyTools::numLinearNodes( const int entity_topology )
{
    int num_nodes = 0;
    
    switch( entity_topology )
    {
	case iMesh_POINT:
	    num_nodes = 1;
	    break;

	case iMesh_LINE_SEGMENT:
	    num_nodes = 2;
	    break;

	case iMesh_TRIANGLE:
	    num_nodes = 3;
	    break;

	case iMesh_QUADRILATERAL:
	    num_nodes = 4;
	    break;

	case iMesh_TETRAHEDRON:
	    num_nodes = 4;
	    break;

	case iMesh_HEXAHEDRON:
	    num_nodes = 8;
	    break;

	default:
	    num_nodes = 0;
    }

    return num_nodes;
}

/*!
 * \brief Reorder a list of element nodes from MBCN ordering to Shards
 * ordering. 
 */
void TopologyTools::MBCN2Shards( iBase_EntityHandle *element_nodes, 
				 const int num_nodes,
				 const int entity_topology )
{
    std::vector<iBase_EntityHandle> temp_nodes( num_nodes );

    switch( entity_topology )
    {
	case iMesh_LINE_SEGMENT:

	    switch( num_nodes )
	    {
		default:
		    for ( int n = 0; n < num_nodes; ++n )
		    {
			temp_nodes[n] = element_nodes[n];
		    }
		    break;
	    }

	    break;

	case iMesh_TRIANGLE:

	    switch( num_nodes )
	    {
		default:
		    for ( int n = 0; n < num_nodes; ++n )
		    {
			temp_nodes[n] = element_nodes[n];
		    }
		    break;
	    }

	    break;

	case iMesh_QUADRILATERAL:

	    switch( num_nodes )
	    {
		default:
		    for ( int n = 0; n < num_nodes; ++n )
		    {
			temp_nodes[n] = element_nodes[n];
		    }
		    break;

	    }

	    break;

	case iMesh_TETRAHEDRON:

	    switch( num_nodes )
	    {
		default:
		    
		    for ( int n = 0; n < num_nodes; ++n )
		    {
			temp_nodes[n] = element_nodes[n];
		    }
		    break;
	    }

	    break;

	case iMesh_HEXAHEDRON:

	    switch( num_nodes )
	    {
		case 27:

		    for ( int n = 0; n < 20; ++n )
		    {
			temp_nodes[n] = element_nodes[n];
		    }
		    temp_nodes[20] = element_nodes[26];
		    temp_nodes[21] = element_nodes[24];
		    temp_nodes[22] = element_nodes[25];
		    temp_nodes[23] = element_nodes[23];
		    temp_nodes[24] = element_nodes[21];
		    temp_nodes[25] = element_nodes[20];
		    temp_nodes[26] = element_nodes[22];
		    break;

		default:

		    for ( int n = 0; n < num_nodes; ++n )
		    {
			temp_nodes[n] = element_nodes[n];
		    }
		    break;
	    }
	    break;

	default:

	    for ( int n = 0; n < num_nodes; ++n )
	    {
		temp_nodes[n] = element_nodes[n];
	    }
    }

    for ( int n = 0; n < num_nodes; ++n )
    {
	element_nodes[n] = temp_nodes[n];
    }
}

/*!
 * \brief Point in volume query on an entity.
 */
bool TopologyTools::pointInVolume( const iMesh_Instance mesh,
				   const iBase_EntityHandle entity,
				   const double coords[3] )
{
    int error = 0;
    int topology = 0;
    iMesh_getEntTopo( mesh, entity, &topology, &error );
    assert( iBase_SUCCESS == error );

    iBase_EntityHandle *element_nodes = 0;
    int element_nodes_allocated = 0;
    int element_nodes_size = 0;
    iMesh_getEntAdj( mesh,
		     entity,
		     iBase_VERTEX,
		     &element_nodes,
		     &element_nodes_allocated,
		     &element_nodes_size,
		     &error );
    assert( iBase_SUCCESS == error );

    MBCN2Shards( element_nodes, element_nodes_size, topology );

    int num_linear_nodes = numLinearNodes( topology );

    CellTopologyFactory topo_factory;
    Teuchos::RCP<shards::CellTopology> cell_topo = 
	topo_factory.create( topology, num_linear_nodes );

    std::vector<iBase_EntityHandle> linear_nodes;
    for ( int n = 0; n < num_linear_nodes; ++n )
    {
	linear_nodes.push_back( element_nodes[n] );
    }

    int coords_allocated = 0;
    int coords_size = 0;
    double *coord_array = 0;
    iMesh_getVtxArrCoords( mesh,
			   &linear_nodes[0],
			   num_linear_nodes,
			   iBase_INTERLEAVED,
			   &coord_array,
			   &coords_allocated,
			   &coords_size,
			   &error );
    assert( iBase_SUCCESS == error );

    Teuchos::Tuple<int,3> cell_node_dimensions;
    cell_node_dimensions[0] = 1;
    cell_node_dimensions[1] = num_linear_nodes;
    cell_node_dimensions[2] = 3;
    Intrepid::FieldContainer<double> cell_nodes( 
	Teuchos::Array<int>(cell_node_dimensions), coord_array );

    Intrepid::FieldContainer<double> find_coords(1,3);
    find_coords(0,0) = coords[0];
    find_coords(0,1) = coords[1];
    find_coords(0,2) = coords[2];
    Intrepid::FieldContainer<double> reference_points(1,3);
    Intrepid::CellTools<double>::mapToReferenceFrame( reference_points,
						      find_coords,
						      cell_nodes,
						      *cell_topo,
						      0 );

    bool return_val = Intrepid::CellTools<double>::checkPointInclusion( 
	&reference_points[0],
	3,
	*cell_topo);

    free( element_nodes );
    free( coord_array );
    cell_topo = Teuchos::null;

    return return_val;
}

} // end namespace FOOD

//---------------------------------------------------------------------------//
// end TopologyTools.cpp
//---------------------------------------------------------------------------//
