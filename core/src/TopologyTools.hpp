//---------------------------------------------------------------------------//
// \file TopologyTools.hpp
// \author Stuart Slattery
// \brief TopologyTools declaration.
//---------------------------------------------------------------------------//

#ifndef FOOD_TOPOLOGYTOOLS_HPP
#define FOOD_TOPOLOGYTOOLS_HPP

#include "Types.hpp"

#include <iBase.h>
#include <iMesh.h>

#include <Intrepid_FieldContainer.hpp>

namespace FOOD
{

namespace TopologyTools
{

// Get the number of linear nodes for a particular iMesh topology.
int numLinearNodes( const int entity_topology );

// Reorder a list of element nodes from MBCN ordering to Shards ordering.
void MBCN2Shards( iBase_EntityHandle *element_nodes, 
		  const int num_nodes,
		  const int entity_topology );

// Get the coordinates of the reference cell for the given topology.
void getReferenceCoords( Intrepid::FieldContainer<double> &ref_coords,
			 const int num_nodes,
			 const int entity_topology );

} // end namespace TopologyTools

} // end namepsace FOOD

#endif // end FOOD_TOPOLOGYTOOLS_HPP

//---------------------------------------------------------------------------//
// end TopologyTools.hpp
//---------------------------------------------------------------------------//
