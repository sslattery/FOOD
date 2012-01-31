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

#include <Intrepid_FieldContainer.hpp>

namespace FOOD
{

struct OctreeNode
{
    iBase_EntitySetHandle node_set = 0;
    Teuchos::RCP<OctreeNode> parent = 0;
    Teuchos::RCP<OctreeNode> child1 = 0;
    Teuchos::RCP<OctreeNode> child2 = 0;
    Teuchos::Tuple<double,6> bounding_box;
    int cutting_plane = 0;
    bool is_leaf = false;
};

class Octree
{

  public:

    //@{
    //! Typedefs.
    typedef Teuchos::RCP<Domain>                      RCP_Domain;
    typedef Teuchos::RCP<OctreeNode>                  RCP_Node;
    typedef Intrepid::FieldContainer<double>          MDArray;
    typedef Teuchos::Tuple<int,6>                     Box;
    //@}

  private:

    // The domain this octree is generated for.
    RCP_Domain d_domain;

    // The root node of the tree.
    RCP_Node d_root_node;

  public:

    // Constructor.
    Octree( RCP_Domain domain );

    // Destructor.
    ~Octree();

    // Build the tree for a specific entity topology.
    void buildTree( const int entity_topology );

    // Locate a point. Return false if we didnt find it.
    bool findPoint( iBase_EntityHandle *found_in_entity,
		    const MDArray &coords );

  private:

    // Get the bounding box of a set of entities.
    void getEntSetBox( const iBase_EntitySetHandle entity_set,
		       Box &bounding_box );

    // Determine if a point is inside a bounding box.
    bool isPointInBox( const Box &box,
		       const MDArray &coords );

    // Determine if an entity is inside a bounding box.
    bool isEntInBox( const Box &box,
		     const iBase_EntityHandle );

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
