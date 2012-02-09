//---------------------------------------------------------------------------//
// \file TopologyTools.hpp
// \author Stuart Slattery
// \brief TopologyTools definition.
//---------------------------------------------------------------------------//

#ifndef FOOD_TOPOLOGYTOOLS_HPP
#define FOOD_TOPOLOGYTOOLS_HPP

#include "Types.hpp"

#include <iBase.h>
#include <iMesh.h>

namespace FOOD
{

namespace TopologyTools
{

/*! 
 * \brief Get the number of linear nodes for a particular iMesh topology.
 */
int numLinearNodes( const int entity_topology )
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
    }

    return num_nodes;
}

} // end namespace TopologyTools

} // end namepsace FOOD

#endif // end FOOD_TOPOLOGYTOOLS_HPP

//---------------------------------------------------------------------------//
// end TopologyTools.hpp
//---------------------------------------------------------------------------//
