//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file TopologyTools.hpp
 * \author Stuart Slattery
 * \brief TopologyTools declaration.
 */
//---------------------------------------------------------------------------//

#ifndef FOOD_TOPOLOGYTOOLS_HPP
#define FOOD_TOPOLOGYTOOLS_HPP

#include "DiscretizationTypes.hpp"

#include <iBase.h>
#include <iMesh.h>

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

// Point in volume query on an entity.
bool pointInVolume( const iMesh_Instance mesh,
		    const iBase_EntityHandle entity,
		    const double coords[3] );

} // end namespace TopologyTools

} // end namepsace FOOD

#endif // end FOOD_TOPOLOGYTOOLS_HPP

//---------------------------------------------------------------------------//
// end TopologyTools.hpp
//---------------------------------------------------------------------------//
