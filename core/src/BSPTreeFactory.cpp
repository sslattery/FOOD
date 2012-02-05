//---------------------------------------------------------------------------//
// \file BSPTreeFactory.cpp
// \author Stuart Slattery
// \brief BSPTreeFactory definition
//---------------------------------------------------------------------------//

#include "BSPTreeFactory.hpp"

namespace FOOD
{

/*! 
 * \brief Constructor.
 */
BSPTreeFactory::BSPTreeFactory()
{ /* ... */ }

/*! 
 * \brief Destructor.
 */
BSPTreeFactory::~BSPTreeFactory()
{ /* ... */ }

/*! 
 * \brief Factory method.
 */
BSPTreeFactory::RCP_BSPTree BSPTreeFactory::create( const int tree_type,
						    RCP_Domain domain, 
						    const int entity_type,
						    const int entity_topology )
{
    RCP_BSPTree new_tree;

    switch( tree_type )
    {
	case OCTREE:

	    new_tree = Teuchos::rcp( 
		new Octree( domain, entity_type, entity_topology ) );
	    break;

	case KDTREE:
	    new_tree = Teuchos::rcp( 
		new KDTree( domain, entity_type, entity_topology ) );
	    break;
    }

    return new_tree;
}

} // end namespace FOOD

//---------------------------------------------------------------------------//
// end BSPTreeFactory.cpp
//---------------------------------------------------------------------------//
