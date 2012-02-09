//---------------------------------------------------------------------------//
// \file TopologyTools.cpp
// \author Stuart Slattery
// \brief TopologyTools definition
//---------------------------------------------------------------------------//

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

} // end namespace FOOD

//---------------------------------------------------------------------------//
// end TopologyTools.cpp
//---------------------------------------------------------------------------//
