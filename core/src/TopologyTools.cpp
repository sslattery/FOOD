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

#include "TopologyTools.hpp"

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
		case 2:
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
		case 3:
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
		case 4:
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
		case 4:
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
		case 8:

		    for ( int n = 0; n < num_nodes; ++n )
		    {
			temp_nodes[n] = element_nodes[n];
		    }
		    break;

		case 20:

		    break;

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

} // end namespace FOOD

//---------------------------------------------------------------------------//
// end TopologyTools.cpp
//---------------------------------------------------------------------------//
