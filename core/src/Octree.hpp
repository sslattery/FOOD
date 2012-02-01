//---------------------------------------------------------------------------//
// \file Octree.hpp
// \author Stuart Slattery
// \breif Octree definition.
//---------------------------------------------------------------------------//

#ifndef FOOD_OCTREE_HPP
#define FOOD_OCTREE_HPP

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

class OctreeNode
{
  public:
    iBase_EntitySetHandle node_set;
    Teuchos::RCP<OctreeNode> parent;
    Teuchos::RCP<OctreeNode> child1;
    Teuchos::RCP<OctreeNode> child2;
    Teuchos::Tuple<double,6> bounding_box;
    int cutting_axis;
    bool is_leaf;
 
    OctreeNode()
       : node_set(0)
       , parent(0)
       , child1(0)
       , child2(0)
       , cutting_axis(0)
       , is_leaf(false)
    { /* ... */ }

    ~OctreeNode()
    { /* ... */ }
};

class Octree
{

  public:

    //@{
    //! Typedefs.
    typedef Teuchos::RCP<Domain>                      RCP_Domain;
    typedef Teuchos::RCP<OctreeNode>                  RCP_Node;
    typedef Intrepid::FieldContainer<double>          MDArray;
    typedef Teuchos::Tuple<double,6>                  Box;
    //@}

  private:

    // Cutting axis enum.
    enum CuttingAxis {
	CuttingAxis_MIN = 0,
	X_Axis = CuttingAxis_MIN,
	Y_Axis,
	Z_Axis,
	CuttingAxis_MAX = Z_Axis 
    };

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
    Octree( RCP_Domain domain, 
	    const int entity_type,
	    const int entity_topology );

    // Destructor.
    ~Octree();

    // Build the tree.
    void buildTree();

    // Locate a point.
    bool findPoint( iBase_EntityHandle *found_in_entity,
		    const MDArray &coords );

  private:

    // Build a tree node.
    void buildTreeNode( RCP_Node node );

    // Search a node for a point.
    bool findPointInNode( RCP_Node node,
			  iBase_EntityHandle *found_in_entity,
			  const MDArray &coords );

    // Get the bounding box of a set of entities.
    void getEntSetBox( iBase_EntitySetHandle entity_set,
		       Box &bounding_box );

    // Determine if a point is inside a bounding box.
    bool isPointInBox( const Box &box,
		       const MDArray &coords );

    // Determine if an entity is inside a bounding box.
    bool isEntInBox( const Box &box,
		     iBase_EntityHandle entity );

    // Slice a box along the specified axis and return the resulting boxes.
    void sliceBox( const Box &parent_box,
		   Box &child_box1,
		   Box &child_box2,
		   const int cutting_axis );
};

} // end namespace FOOD

#endif // end FOOD_OCTREE_HPP

//---------------------------------------------------------------------------//
// end Octree.hpp
//---------------------------------------------------------------------------//
