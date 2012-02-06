//---------------------------------------------------------------------------//
// \file KDTree.hpp
// \author Stuart Slattery
// \brief KDTree declaration.
//---------------------------------------------------------------------------//

#ifndef FOOD_KDTREE_HPP
#define FOOD_KDTREE_HPP

#include "PointQuery.hpp"
#include "Domain.hpp"
#include "BSPTree.hpp"

#include <iBase.h>
#include <iMesh.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayView.hpp>

#include <Intrepid_FieldContainer.hpp>

namespace FOOD
{

class KDTreeNode
{
  public:
    iBase_EntitySetHandle node_set;
    KDTreeNode* parent;
    KDTreeNode* child1;
    KDTreeNode* child2;
    double bounding_box[6];
    bool is_leaf;
    int axis;
    double median;
 
    KDTreeNode()
       : node_set(0)
       , parent(0)
       , child1(0)
       , child2(0)
       , is_leaf(false)
       , axis(0)
       , median(0.0)
    { /* ... */ }

    ~KDTreeNode()
    { /* ... */ }
};

class KDTree : public BSPTree
{
  public:

    //@{
    //! Typedefs.
    typedef Teuchos::RCP<Domain>                      RCP_Domain;
    typedef Intrepid::FieldContainer<double>          MDArray;
    //@}

  private:

    // The domain this octree is generated for.
    RCP_Domain d_domain;

    // The entity type this tree is built on.
    std::size_t d_entity_type;

    // The entity topology this tree is built on.
    std::size_t d_entity_topology;

    // The root node of the tree.
    KDTreeNode* d_root_node;

  public:

    // Constructor.
    KDTree( RCP_Domain domain, 
	    const int entity_type,
	    const int entity_topology );

    // Destructor.
    ~KDTree();

    // Build the tree.
    void buildTree();

    // Locate a point.
    bool findPoint( iBase_EntityHandle &found_in_entity,
		    const MDArray &coords );

  private:

    // Build a tree node.
    void buildTreeNode( KDTreeNode* node );

    // Given a point, find its leaf node in the tree.
    KDTreeNode* findLeafNode( KDTreeNode* node, 
			      iBase_EntityHandle &nearest_neighbor,
			      const MDArray &coords );

    // Search a node for a point.
    bool findPointInNode( KDTreeNode* node,
			  iBase_EntityHandle &nearest_neighbor,
			  const MDArray &coords,
			  iBase_EntityHandle &found_in_entity );

    // Get the bounding box of a set of entities.
    void getEntSetBox( iBase_EntitySetHandle entity_set, 
		       double bounding_box[6] );

    // Determine if a point is inside a bounding box.
    bool isPointInBox( const double box[6], const MDArray &coords );

    // Determine if a point is inside the entities adjacent to another point.
    bool isPointInAdj( iBase_EntityHandle point, 
		       const MDArray &coords,
		       iBase_EntityHandle &found_in_entity );

    // Determine if an entity is inside a bounding box.
    bool isEntInBox( const double box[6], iBase_EntityHandle entity );

    // Slice a box along the specified axis and return the value for slicing.
    void sliceBox( KDTreeNode* node );

    // Compute the median of a range of values.
    double median( Teuchos::ArrayView<double>::iterator begin,
		   Teuchos::ArrayView<double>::iterator end );
};

} // end namespace FOOD

#endif // end FOOD_KDTREE_HPP

//---------------------------------------------------------------------------//
// end KDTree.hpp
//---------------------------------------------------------------------------//
