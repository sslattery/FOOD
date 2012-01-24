//---------------------------------------------------------------------------//
// \file CellTopologyFactory.cpp
// \author Stuart Slattery
// \brief Factory method defintion for cell topology data.
//---------------------------------------------------------------------------//

#include "CellTopologyFactory.hpp"

#include <iMesh.h>

namespace FOOD
{

/*!
 * \brief Constructor.
 */
CellTopologyFactory::CellTopologyFactory()
{ /* ... */ }

/*!
 * \brief Destructor.
 */
CellTopologyFactory::~CellTopologyFactory()
{ /* ... */ }

/*!
 * \brief Factory method.
 */
Teuchos::RCP<shards::CellTopology> 
CellTopologyFactory::create( const int entity_topology )
{
    Teuchos::RCP<shards::CellTopology> new_topology;

    switch( entity_topology )
    {
	// case iMesh_POINT:
	    
	//     new_topology = Teuchos::rcp( 
	// 	new shards::CellTopology(
	// 	    shards::getCellTopologyData< shards::Node<> >() ) );
	//     break;

	case iMesh_LINE_SEGMENT:
	    
	    new_topology = Teuchos::rcp( 
		new shards::CellTopology(
		    shards::getCellTopologyData< shards::Line<> >() ) );
	    break;

	case iMesh_TRIANGLE:
	    
	    new_topology = Teuchos::rcp( 
		new shards::CellTopology(
		    shards::getCellTopologyData< shards::Triangle<> >() ) );
	    break;

	case iMesh_QUADRILATERAL:
	    
	    new_topology = Teuchos::rcp( 
		new shards::CellTopology(
		    shards::getCellTopologyData< shards::Quadrilateral<> >() ) );
	    break;

	case iMesh_TETRAHEDRON:
	    
	    new_topology = Teuchos::rcp( 
		new shards::CellTopology(
		    shards::getCellTopologyData< shards::Tetrahedron<> >() ) );
	    break;

	case iMesh_HEXAHEDRON:
	    
	    new_topology = Teuchos::rcp( 
		new shards::CellTopology(
		    shards::getCellTopologyData< shards::Hexahedron<> >() ) );
	    break;
    }

    return new_topology;
}

}

//---------------------------------------------------------------------------//
// end CellTopologyFactory.cpp
//---------------------------------------------------------------------------//

