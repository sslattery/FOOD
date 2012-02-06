//---------------------------------------------------------------------------//
// \file BSPTree.hpp
// \author Stuart Slattery
// \brief BSPTree protocol declaration.
//---------------------------------------------------------------------------//

#ifndef FOOD_BSPTREE_HPP
#define FOOD_BSPTREE_HPP

#include <iBase.h>

#include <Teuchos_Describable.hpp>

#include <Intrepid_FieldContainer.hpp>

namespace FOOD
{

class BSPTree : public Teuchos::Describable
{

  public:

    //@{
    //! Typedefs.
    typedef Intrepid::FieldContainer<double>          MDArray;
    //@}

    //! Tree type enumerations.
    enum TreeType {
	TreeType_MIN = 0,
	OCTREE = TreeType_MIN,
	KDTREE,
	TreeType_MAX = KDTREE
    };

    //! Constructor.
    BSPTree()
    { /* ... */ }

    //! Destructor.
    virtual ~BSPTree()
    { /* ... */ }

    //! Build the BSPTree.
    virtual void buildTree() = 0;

    //! Search the BSPTree for a point.
    virtual bool findPoint( iBase_EntityHandle &found_in_entity,
			    const MDArray &coords ) = 0;
}; 

} // end namespace FOOD

#endif // end FOOD_BSPTREE_HPP

//---------------------------------------------------------------------------//
// end BSPTree.hpp
//---------------------------------------------------------------------------//
