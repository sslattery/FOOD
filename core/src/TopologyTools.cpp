//---------------------------------------------------------------------------//
// \file TopologyTools.cpp
// \author Stuart Slattery
// \brief TopologyTools definition
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
void TopologyTools::MBCN2Shards( EntityHandle *element_nodes, 
				 const int num_nodes,
				 const int entity_topology )
{
    std::vector<EntityHandle> temp_nodes( num_nodes );

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

/*!
 * \brief Get the coordinates of the reference cell for the given topology.
 */
void TopologyTools::getReferenceCoords( 
    Intrepid::FieldContainer<double> &ref_coords,
    const int num_nodes,
    const int entity_topology )
{
    switch( entity_topology )
    {
	case iMesh_LINE_SEGMENT:

	    switch( num_nodes )
	    {
		case 2:

		    ref_coords(0,0) = 0.0;
		    ref_coords(0,1) = 0.0;
		    ref_coords(0,2) = 0.0;

		    ref_coords(1,0) = 1.0;
		    ref_coords(1,1) = 0.0;
		    ref_coords(1,2) = 0.0;

		    break;
	    }

	    break;

	case iMesh_TRIANGLE:

	    switch( num_nodes )
	    {
		case 3:

		    ref_coords(0,0) = 0.0;
		    ref_coords(0,1) = 0.0;
		    ref_coords(0,2) = 0.0;

		    ref_coords(1,0) = 1.0;
		    ref_coords(1,1) = 0.0;
		    ref_coords(1,2) = 0.0;

		    ref_coords(2,0) = 0.0;
		    ref_coords(2,1) = 1.0;
		    ref_coords(2,2) = 0.0;

		    break;

		case 6:

		    ref_coords(0,0) = 0.0;
		    ref_coords(0,1) = 0.0;
		    ref_coords(0,2) = 0.0;

		    ref_coords(1,0) = 1.0;
		    ref_coords(1,1) = 0.0;
		    ref_coords(1,2) = 0.0;

		    ref_coords(2,0) = 0.0;
		    ref_coords(2,1) = 1.0;
		    ref_coords(2,2) = 0.0;

		    ref_coords(3,0) = 0.5;
		    ref_coords(3,1) = 0.0;
		    ref_coords(3,2) = 0.0;

		    ref_coords(4,0) = 0.5;
		    ref_coords(4,1) = 0.5;
		    ref_coords(4,2) = 0.0;

		    ref_coords(5,0) = 0.0;
		    ref_coords(5,1) = 0.5;
		    ref_coords(5,2) = 0.0;

		    break;
	    }

	    break;

	case iMesh_QUADRILATERAL:

	    switch( num_nodes )
	    {
		case 4:

		    ref_coords(0,0) = -1.0;
		    ref_coords(0,1) = -1.0;
		    ref_coords(0,2) =  0.0;

		    ref_coords(1,0) =  1.0;
		    ref_coords(1,1) = -1.0;
		    ref_coords(1,2) =  0.0;

		    ref_coords(2,0) =  1.0;
		    ref_coords(2,1) =  1.0;
		    ref_coords(2,2) =  0.0;

		    ref_coords(3,0) = -1.0;
		    ref_coords(3,1) =  1.0;
		    ref_coords(3,2) =  0.0;

		    break;

		case 9:

		    ref_coords(0,0) = -1.0;
		    ref_coords(0,1) = -1.0;
		    ref_coords(0,2) =  0.0;

		    ref_coords(1,0) =  1.0;
		    ref_coords(1,1) = -1.0;
		    ref_coords(1,2) =  0.0;

		    ref_coords(2,0) =  1.0;
		    ref_coords(2,1) =  1.0;
		    ref_coords(2,2) =  0.0;

		    ref_coords(3,0) = -1.0;
		    ref_coords(3,1) =  1.0;
		    ref_coords(3,2) =  0.0;

		    ref_coords(4,0) =  0.0;
		    ref_coords(4,1) = -1.0;
		    ref_coords(4,2) =  0.0;

		    ref_coords(5,0) =  1.0;
		    ref_coords(5,1) =  0.0;
		    ref_coords(5,2) =  0.0;

		    ref_coords(6,0) =  0.0;
		    ref_coords(6,1) =  1.0;
		    ref_coords(6,2) =  0.0;

		    ref_coords(7,0) = -1.0;
		    ref_coords(7,1) =  0.0;
		    ref_coords(7,2) =  0.0;

		    ref_coords(8,0) =  0.0;
		    ref_coords(8,1) =  0.0;
		    ref_coords(8,2) =  0.0;

		    break;
	    }

	    break;

	case iMesh_TETRAHEDRON:

	    switch( num_nodes )
	    {
		case 4:

		    ref_coords(0,0) = 0.0;
		    ref_coords(0,1) = 0.0;
		    ref_coords(0,2) = 0.0;

		    ref_coords(1,0) = 1.0;
		    ref_coords(1,1) = 0.0;
		    ref_coords(1,2) = 0.0;

		    ref_coords(2,0) = 0.0;
		    ref_coords(2,1) = 1.0;
		    ref_coords(2,2) = 0.0;

		    ref_coords(3,0) = 0.0;
		    ref_coords(3,1) = 0.0;
		    ref_coords(3,2) = 1.0;

		    break;

		case 10:

		    ref_coords(0,0) = 0.0;
		    ref_coords(0,1) = 0.0;
		    ref_coords(0,2) = 0.0;

		    ref_coords(1,0) = 1.0;
		    ref_coords(1,1) = 0.0;
		    ref_coords(1,2) = 0.0;

		    ref_coords(2,0) = 0.0;
		    ref_coords(2,1) = 1.0;
		    ref_coords(2,2) = 0.0;

		    ref_coords(3,0) = 0.0;
		    ref_coords(3,1) = 0.0;
		    ref_coords(3,2) = 1.0;

		    ref_coords(4,0) = 0.5;
		    ref_coords(4,1) = 0.0;
		    ref_coords(4,2) = 0.0;

		    ref_coords(5,0) = 0.5;
		    ref_coords(5,1) = 0.5;
		    ref_coords(5,2) = 0.0;

		    ref_coords(6,0) = 0.0;
		    ref_coords(6,1) = 0.5;
		    ref_coords(6,2) = 0.0;

		    ref_coords(7,0) = 0.0;
		    ref_coords(7,1) = 0.0;
		    ref_coords(7,2) = 0.5;

		    ref_coords(8,0) = 0.5;
		    ref_coords(8,1) = 0.0;
		    ref_coords(8,2) = 0.5;

		    ref_coords(9,0) = 0.0;
		    ref_coords(9,1) = 0.5;
		    ref_coords(9,2) = 0.5;
	    }

	    break;

	case iMesh_HEXAHEDRON:

	    switch( num_nodes )
	    {
		case 8:

		    ref_coords(0,0) = -1.0;
		    ref_coords(0,1) = -1.0;
		    ref_coords(0,2) = -1.0;

		    ref_coords(1,0) =  1.0;
		    ref_coords(1,1) = -1.0;
		    ref_coords(1,2) = -1.0;

		    ref_coords(2,0) =  1.0;
		    ref_coords(2,1) =  1.0;
		    ref_coords(2,2) = -1.0;

		    ref_coords(3,0) = -1.0;
		    ref_coords(3,1) =  1.0;
		    ref_coords(3,2) = -1.0;

		    ref_coords(4,0) = -1.0;
		    ref_coords(4,1) = -1.0;
		    ref_coords(4,2) =  1.0;

		    ref_coords(5,0) =  1.0;
		    ref_coords(5,1) = -1.0;
		    ref_coords(5,2) =  1.0;

		    ref_coords(6,0) =  1.0;
		    ref_coords(6,1) =  1.0;
		    ref_coords(6,2) =  1.0;

		    ref_coords(7,0) = -1.0;
		    ref_coords(7,1) =  1.0;
		    ref_coords(7,2) =  1.0;

		    break;

		case 20:

		    ref_coords(0,0) = -1.0;
		    ref_coords(0,1) = -1.0;
		    ref_coords(0,2) = -1.0;

		    ref_coords(1,0) =  1.0;
		    ref_coords(1,1) = -1.0;
		    ref_coords(1,2) = -1.0;

		    ref_coords(2,0) =  1.0;
		    ref_coords(2,1) =  1.0;
		    ref_coords(2,2) = -1.0;

		    ref_coords(3,0) = -1.0;
		    ref_coords(3,1) =  1.0;
		    ref_coords(3,2) = -1.0;

		    ref_coords(4,0) = -1.0;
		    ref_coords(4,1) = -1.0;
		    ref_coords(4,2) =  1.0;

		    ref_coords(5,0) =  1.0;
		    ref_coords(5,1) = -1.0;
		    ref_coords(5,2) =  1.0;

		    ref_coords(6,0) =  1.0;
		    ref_coords(6,1) =  1.0;
		    ref_coords(6,2) =  1.0;

		    ref_coords(7,0) = -1.0;
		    ref_coords(7,1) =  1.0;
		    ref_coords(7,2) =  1.0;

		    ref_coords(8,0) =  0.0;
		    ref_coords(8,1) = -1.0;
		    ref_coords(8,2) = -1.0;

		    ref_coords(9,0) =  1.0;
		    ref_coords(9,1) =  0.0;
		    ref_coords(9,2) = -1.0;

		    ref_coords(10,0) =  0.0;
		    ref_coords(10,1) =  1.0;
		    ref_coords(10,2) = -1.0;

		    ref_coords(11,0) = -1.0;
		    ref_coords(11,1) =  0.0;
		    ref_coords(11,2) = -1.0;

		    ref_coords(12,0) = -1.0;
		    ref_coords(12,1) = -1.0;
		    ref_coords(12,2) =  0.0;

		    ref_coords(13,0) =  1.0;
		    ref_coords(13,1) = -1.0;
		    ref_coords(13,2) =  0.0;

		    ref_coords(14,0) =  1.0;
		    ref_coords(14,1) =  1.0;
		    ref_coords(14,2) =  0.0;

		    ref_coords(15,0) = -1.0;
		    ref_coords(15,1) =  1.0;
		    ref_coords(15,2) =  0.0;

		    ref_coords(16,0) =  0.0;
		    ref_coords(16,1) = -1.0;
		    ref_coords(16,2) =  1.0;

		    ref_coords(17,0) =  1.0;
		    ref_coords(17,1) =  0.0;
		    ref_coords(17,2) =  1.0;

		    ref_coords(18,0) =  1.0;
		    ref_coords(18,1) =  1.0;
		    ref_coords(18,2) = -1.0;

		    ref_coords(19,0) = -1.0;
		    ref_coords(19,1) =  0.0;
		    ref_coords(19,2) = -1.0;

		    break;

		case 27:

		    ref_coords(0,0) = -1.0;
		    ref_coords(0,1) = -1.0;
		    ref_coords(0,2) = -1.0;

		    ref_coords(1,0) =  1.0;
		    ref_coords(1,1) = -1.0;
		    ref_coords(1,2) = -1.0;

		    ref_coords(2,0) =  1.0;
		    ref_coords(2,1) =  1.0;
		    ref_coords(2,2) = -1.0;

		    ref_coords(3,0) = -1.0;
		    ref_coords(3,1) =  1.0;
		    ref_coords(3,2) = -1.0;

		    ref_coords(4,0) = -1.0;
		    ref_coords(4,1) = -1.0;
		    ref_coords(4,2) =  1.0;

		    ref_coords(5,0) =  1.0;
		    ref_coords(5,1) = -1.0;
		    ref_coords(5,2) =  1.0;

		    ref_coords(6,0) =  1.0;
		    ref_coords(6,1) =  1.0;
		    ref_coords(6,2) =  1.0;

		    ref_coords(7,0) = -1.0;
		    ref_coords(7,1) =  1.0;
		    ref_coords(7,2) =  1.0;

		    ref_coords(8,0) =  0.0;
		    ref_coords(8,1) = -1.0;
		    ref_coords(8,2) = -1.0;

		    ref_coords(9,0) =  1.0;
		    ref_coords(9,1) =  0.0;
		    ref_coords(9,2) = -1.0;

		    ref_coords(10,0) =  0.0;
		    ref_coords(10,1) =  1.0;
		    ref_coords(10,2) = -1.0;

		    ref_coords(11,0) = -1.0;
		    ref_coords(11,1) =  0.0;
		    ref_coords(11,2) = -1.0;

		    ref_coords(12,0) = -1.0;
		    ref_coords(12,1) = -1.0;
		    ref_coords(12,2) =  0.0;

		    ref_coords(13,0) =  1.0;
		    ref_coords(13,1) = -1.0;
		    ref_coords(13,2) =  0.0;

		    ref_coords(14,0) =  1.0;
		    ref_coords(14,1) =  1.0;
		    ref_coords(14,2) =  0.0;

		    ref_coords(15,0) = -1.0;
		    ref_coords(15,1) =  1.0;
		    ref_coords(15,2) =  0.0;

		    ref_coords(16,0) =  0.0;
		    ref_coords(16,1) = -1.0;
		    ref_coords(16,2) =  1.0;

		    ref_coords(17,0) =  1.0;
		    ref_coords(17,1) =  0.0;
		    ref_coords(17,2) =  1.0;

		    ref_coords(18,0) =  1.0;
		    ref_coords(18,1) =  1.0;
		    ref_coords(18,2) = -1.0;

		    ref_coords(19,0) = -1.0;
		    ref_coords(19,1) =  0.0;
		    ref_coords(19,2) = -1.0;

		    ref_coords(20,0) =  0.0;
		    ref_coords(20,1) =  0.0;
		    ref_coords(20,2) =  0.0;

		    ref_coords(21,0) =  0.0;
		    ref_coords(21,1) =  0.0;
		    ref_coords(21,2) = -1.0;

		    ref_coords(22,0) =  0.0;
		    ref_coords(22,1) =  0.0;
		    ref_coords(22,2) =  1.0;

		    ref_coords(23,0) = -1.0;
		    ref_coords(23,1) =  0.0;
		    ref_coords(23,2) =  0.0;

		    ref_coords(24,0) =  1.0;
		    ref_coords(24,1) =  0.0;
		    ref_coords(24,2) =  0.0;

		    ref_coords(25,0) =  0.0;
		    ref_coords(25,1) = -1.0;
		    ref_coords(25,2) =  0.0;

		    ref_coords(26,0) =  0.0;
		    ref_coords(26,1) =  1.0;
		    ref_coords(26,2) =  0.0;

		    break;
	    }
	    break;

	default:
	    
	    for ( int i = 0; i < ref_coords.dimension(0); ++i )
	    {
		for ( int j = 0; j < ref_coords.dimension(1); ++j )
		{
		    ref_coords(i,j) = 0.0;
		}
	    }
    }
}

} // end namespace FOOD

//---------------------------------------------------------------------------//
// end TopologyTools.cpp
//---------------------------------------------------------------------------//
