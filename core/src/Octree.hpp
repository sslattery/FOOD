//---------------------------------------------------------------------------//
// \file Octree.hpp
// \author Stuart Slattery
// \brief Octree declaration.
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
    Teuchos::Tuple< Teuchos::RCP<OctreeNode>, 8 > children;
    Teuchos::Tuple<double,6> bounding_box;
    bool is_leaf;
 
    OctreeNode()
       : node_set(0)
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
    typedef iBase_EntityHandle                        EntityHandle;
    //@}

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
    Octree( RCP_Domain domain, 
	    const int entity_type,
	    const int entity_topology );

    // Destructor.
    ~Octree();

    // Build the tree.
    void buildTree();

    // Locate a point.
    bool findPoint( EntityHandle &found_in_entity,
		    const MDArray &coords );

  private:

    // Build a tree node.
    void buildTreeNode( RCP_Node node );

    // Search a node for a point.
    bool findPointInNode( RCP_Node node,
			  EntityHandle &found_in_entity,
			  const MDArray &coords );

    // Get the bounding box of a set of entities.
    void getEntSetBox( iBase_EntitySetHandle entity_set, Box &bounding_box );

    // Determine if a point is inside a bounding box.
    bool isPointInBox( const Box &box, const MDArray &coords );

    // Determine if an entity is inside a bounding box.
    bool isEntInBox( const Box &box, EntityHandle entity );

    // Slice a node box into 8 children boxes.
    void sliceBox( RCP_Node node );
};

} // end namespace FOOD

#endif // end FOOD_OCTREE_HPP

//---------------------------------------------------------------------------//
// end Octree.hpp
//---------------------------------------------------------------------------//
