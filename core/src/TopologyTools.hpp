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

namespace FOOD
{

namespace TopologyTools
{

// Get the number of linear nodes for a particular iMesh topology.
int numLinearNodes( const int entity_topology );

} // end namespace TopologyTools

} // end namepsace FOOD

#endif // end FOOD_TOPOLOGYTOOLS_HPP

//---------------------------------------------------------------------------//
// end TopologyTools.hpp
//---------------------------------------------------------------------------//
