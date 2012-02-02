//---------------------------------------------------------------------------//
// \file KDTree.hpp
// \author Stuart Slattery
// \breif KDTree definition.
//---------------------------------------------------------------------------//

#ifndef FOOD_KDTREE_HPP
#define FOOD_KDTREE_HPP

#include "PointQuery.hpp"
#include "Domain.hpp"

#include <iBase.h>
#include <iMesh.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_ArrayView.hpp>

#include <Intrepid_FieldContainer.hpp>

namespace FOOD
{

class KDTreeNode
{
    iBase_EntitySetHandle node_set;
    Teuchos::RCP<KDTreeNode> child1;
    Teuchos::RCP<KDTreeNode> child2;
    Teuchos::Tuple<double,6> bounding_box;
    bool is_leaf;
 
    KDTreeNode()
       : node_set(0)
       , child1(0)
       , child2(0)
       , is_leaf(false)
    { /* ... */ }

    ~KDTreeNode()
    { /* ... */ }
};

class KDTree
{

  private:

    // The domain this octree is generated for.
    RCP_Domain d_domain;

    // The entity type this tree is built on.
    std::size_t d_entity_type;

    // The entity topology this tree is built on.
    std::size_t d_entity_topology;

    // The root node of the tree.
    RCP_Node d_root_node;

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
    void buildTreeNode( RCP_Node node );

    // Search a node for a point.
    bool findPointInNode( RCP_Node node,
			  iBase_EntityHandle &found_in_entity,
			  const MDArray &coords );

    // Get the bounding box of a set of entities.
    void getEntSetBox( iBase_EntitySetHandle entity_set, Box &bounding_box );

    // Determine if a point is inside a bounding box.
    bool isPointInBox( const Box &box, const MDArray &coords );

    // Determine if an entity is inside a bounding box.
    bool isEntInBox( const Box &box, iBase_EntityHandle entity );

    // Slice a box along the specified axis and return the resulting boxes.
    void sliceBox( RCP_Node node, const int axis, const double axis_coord );
};

} // end namespace FOOD

#endif // end FOOD_KDTREE_HPP

//---------------------------------------------------------------------------//
// end KDTree.hpp
//---------------------------------------------------------------------------//
