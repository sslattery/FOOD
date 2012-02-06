//---------------------------------------------------------------------------//
// \file BSPTreeFactory.hpp
// \author Stuart Slattery
// \brief BSPTreeFactory declaration
//---------------------------------------------------------------------------//

#ifndef FOOD_BSPTREEFACTORY_HPP
#define FOOD_BSPTREEFACTORY_HPP

#include "Domain.hpp"
#include "BSPTree.hpp"
#include "Octree.hpp"
#include "KDTree.hpp"

#include <Teuchos_RCP.hpp>

namespace FOOD
{

class BSPTreeFactory
{
  public:
   
    //@{
    //! Typedefs.
    typedef Teuchos::RCP<Domain>               RCP_Domain;
    typedef Teuchos::RCP<BSPTree>              RCP_BSPTree;
    //@}

    //! Constructor.
    BSPTreeFactory();

    //! Destructor.
    ~BSPTreeFactory();

    //! Factory method.
    RCP_BSPTree create( const int tree_type,
			RCP_Domain domain, 
			const int entity_type,
			const int entity_topology );
}; 

} // end namespace FOOD

#endif // end FOOD_BSPTREEFACTORY_HPP

//---------------------------------------------------------------------------//
// end BSPTreeFactory.hpp
//---------------------------------------------------------------------------//
