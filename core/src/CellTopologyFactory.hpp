//---------------------------------------------------------------------------//
// \file CellTopologyFactory.hpp
// \author Stuart Slattery
// \brief Factory method declaration for cell topology data.
//---------------------------------------------------------------------------//

#ifndef FOOD_CELLTOPOLOGYFACTORY_HPP
#define FOOD_CELLTOPOLOGYFACTORY_HPP

#include <Teuchos_RCP.hpp>

#include <Shards_CellTopology.hpp>

namespace FOOD
{

class CellTopologyFactory
{
  public:

    // Consructor.
    CellTopologyFactory();

    // Destructor.
    ~CellTopologyFactory();

    // Factor method.
    Teuchos::RCP<shards::CellTopology> create( const int entity_topology );
};

} // end namespace FOOD

#endif // end FOOD_CELLTOPOLOGYFACTORY_HPP

//---------------------------------------------------------------------------//
// end CellTopologyFactory.hpp
//---------------------------------------------------------------------------//
