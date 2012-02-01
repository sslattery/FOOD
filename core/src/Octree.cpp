//---------------------------------------------------------------------------//
// \file Octree.cpp
// \author Stuart Slattery
// \breif Octree declaration.
//---------------------------------------------------------------------------//

#include <cassert>
#include <algorithm>

#include "Octree.hpp"

namespace FOOD
{

/*! 
 * \brief Constructor.
 */
Octree::Octree( RCP_Domain domain, 
		const int entity_type,
		const int entity_topology )
    : d_domain(domain)
    , d_entity_type(entity_type)
    , d_entity_topology(entity_topology)
    , d_root_node( Teuchos::rcp(new OctreeNode) )
{
    d_root_node->node_set = d_domain->getMeshSet();
    d_root_node->cutting_axis = X_Axis;
    getEntSetBox( d_root_node->node_set, d_root_node->bounding_box );   
}

/*!
 * \brief Destructor.
 */
Octree::~Octree()
{ /* ... */ }

/*! 
 * \brief Build the tree.
 */
void Octree::buildTree()
{
    buildTreeNode( d_root_node );
}

/*!
 * \brief Locate a point. Return false if we didn't find it.
 */
bool Octree::findPoint( iBase_EntityHandle *found_in_entity,
			const MDArray &coords )
{
    return findPointInNode( d_root_node, found_in_entity, coords );
}

/*! 
 * \brief Build a tree node.
 */
void Octree::buildTreeNode( RCP_Node node )
{
    int error = 0;

    // See if we have more than 1 node. If not, we're at a leaf.
    int num_node_elements = 0;
    iMesh_getNumOfTopo( d_domain->getMesh(),
			node->node_set,
			d_entity_topology,
			&num_node_elements,
			&error );
    assert( iBase_SUCCESS == error );

    if ( num_node_elements > 1 )
    {
	// Create the children.
	RCP_Node child1 = Teuchos::rcp( new OctreeNode );
	RCP_Node child2 = Teuchos::rcp( new OctreeNode );

	// Define the new cutting axis.
	int cutting_axis = 0;
	if ( node->cutting_axis == 0 )
	{
	    cutting_axis = 1;
	}
	else if ( node->cutting_axis == 1 )
	{
	    cutting_axis = 2;
	}

	// Create the new bounding boxes and assign them to the children.
	sliceBox( node->bounding_box,
		  child1->bounding_box,
		  child2->bounding_box,
		  cutting_axis );

	// Create the entity sets in the child nodes.
	iMesh_createEntSet( d_domain->getMesh(),
			    1,
			    &(child1->node_set),
			    &error );
	assert( iBase_SUCCESS == error );

	iMesh_createEntSet( d_domain->getMesh(),
			    1,
			    &(child2->node_set),
			    &error );
	assert( iBase_SUCCESS == error );

	// Add the elements in the parent set to the child sets.
	int node_elements_allocated = 0;
	int node_elements_size = 0;
	iBase_EntityHandle *node_elements;
	iMesh_getEntities( d_domain->getMesh(),
			   node->node_set,
			   d_entity_type,
			   d_entity_topology,
			   &node_elements,
			   &node_elements_allocated,
			   &node_elements_size,
			   &error );
	assert( iBase_SUCCESS == error );

	for (int n = 0; n < node_elements_size; ++n )
	{
	    if ( isEntInBox( child1->bounding_box, node_elements[n] ) )
	    {
		iMesh_addEntToSet( d_domain->getMesh(),
				   node_elements[n],
				   child1->node_set,
				   &error );
		assert( iBase_SUCCESS == error );
	    }
	    else if ( isEntInBox( child2->bounding_box, node_elements[n] ) )
	    {
		iMesh_addEntToSet( d_domain->getMesh(),
				   node_elements[n],
				   child2->node_set,
				   &error );
		assert( iBase_SUCCESS == error );
	    }
	}

	free( node_elements );

	// Set relationships.
	node->child1 = child1;
	node->child2 = child2;
	child1->parent = node;
	child2->parent = node;

	// Recurse.
	buildTreeNode( node->child1 );
	buildTreeNode( node->child2 );
    }
    else
    {
	std::cout << "FOUND A LEAF " << std::endl;
	node->is_leaf = true;
    }
}

/*!
 * \brief Search a node for a point. If its a leaf node, return the element we
 * found it in.
 */
bool Octree::findPointInNode( RCP_Node node,
			      iBase_EntityHandle *found_in_entity,
			      const MDArray &coords )
{
    int error = 0;
    bool return_val = false;

    if ( !node->is_leaf )
    {
	// See which child the point is in.
	if ( isPointInBox( node->child1->bounding_box, coords ) )
	{
	    return_val = findPointInNode( node->child1,
					  found_in_entity,
					  coords );
	}
	else if ( isPointInBox( node->child2->bounding_box, coords ) )
	{
	    return_val = findPointInNode( node->child2,
					  found_in_entity,
					  coords );
	}
    }

    else
    {
	int node_elements_allocated = 0; 
	int node_elements_size = 0; 
	iMesh_getEntities( d_domain->getMesh(),
			   node->node_set,
			   d_entity_type,
			   d_entity_topology,
			   &found_in_entity,
			   &node_elements_allocated,
			   &node_elements_size,
			   &error );
	assert( iBase_SUCCESS == error );

	if ( node_elements_size > 0 )
	{
	    if ( PointQuery::point_in_ref_element( d_domain->getMesh(),
						   found_in_entity[0],
						   coords ) )
	    {
		return_val = true;
	    }
	}
    }

    return return_val;
}

/*!
 * \brief Get the bounding box of a set of entities.
 */
void Octree::getEntSetBox( iBase_EntitySetHandle entity_set,
			   Box &bounding_box )
{
    int error = 0;
    
    int num_set_vertices = 0;
    iMesh_getNumOfTopo( d_domain->getMesh(),
			entity_set,
			iMesh_POINT,
			&num_set_vertices,
			&error );
    assert( iBase_SUCCESS == error );

    iBase_EntityHandle *set_vertices = 0;
    int set_vertices_allocated = 0;
    int set_vertices_size = 0;
    iMesh_getEntities( d_domain->getMesh(),
		       entity_set,
		       iBase_VERTEX,
		       iMesh_POINT,
		       &set_vertices,
		       &set_vertices_allocated,
		       &set_vertices_size,
		       &error );
    assert( iBase_SUCCESS == error );

    int set_coords_allocated = 3*set_vertices_size;
    int set_coords_size = 0;
    double *coords = 0;
    iMesh_getVtxArrCoords( d_domain->getMesh(),
			   set_vertices,
			   set_vertices_size,
			   iBase_BLOCKED,
			   &coords,
			   &set_coords_allocated,
			   &set_coords_size,
			   &error );
    assert( iBase_SUCCESS == error );

    Teuchos::ArrayView<double> set_coords( coords, set_coords_size );

    Teuchos::ArrayView<double>::const_iterator set_x_it_begin 
	= set_coords.begin();

    Teuchos::ArrayView<double>::const_iterator set_x_it_end 
	= set_x_it_begin + set_vertices_size;

    Teuchos::ArrayView<double>::const_iterator set_y_it_begin 
	= set_x_it_end;
    Teuchos::ArrayView<double>::const_iterator set_y_it_end 
	= set_y_it_begin + set_vertices_size;
    
    Teuchos::ArrayView<double>::const_iterator set_z_it_begin 
	= set_y_it_end;
    Teuchos::ArrayView<double>::const_iterator set_z_it_end 
	= set_z_it_begin + set_vertices_size;

    bounding_box[0] = *(std::min_element( set_x_it_begin, set_x_it_end ));
    bounding_box[1] = *(std::max_element( set_x_it_begin, set_x_it_end ));
    bounding_box[2] = *(std::min_element( set_y_it_begin, set_y_it_end ));
    bounding_box[3] = *(std::max_element( set_y_it_begin, set_y_it_end ));
    bounding_box[4] = *(std::min_element( set_z_it_begin, set_z_it_end ));
    bounding_box[5] = *(std::max_element( set_z_it_begin, set_z_it_end ));

    free( set_vertices );
    free( coords );
}

/*!
 * \brief Determine if a point is inside a bounding box.
 */
bool Octree::isPointInBox( const Box &box,
			   const MDArray &coords )
{
    bool return_val = false;
    
    if ( coords(0,0) >= box[0] &&
	 coords(0,0) <= box[1] &&
	 coords(0,1) >= box[2] &&
	 coords(0,1) <= box[3] &&
	 coords(0,2) >= box[4] &&
	 coords(0,2) <= box[5] )
    {
	return_val = true;
    }

    return return_val;
}

/*!
 * \brief Determine if an entity is inside a bounding box. The entire element
 * must be in the box.
 */
bool Octree::isEntInBox( const Box &box,
			 iBase_EntityHandle entity )
{
    int error = 0;
    bool return_val = true;

    iBase_EntityHandle *element_nodes = 0;
    int element_nodes_allocated = 0;
    int element_nodes_size = 0;
    iMesh_getEntAdj( d_domain->getMesh(),
		     entity,
		     iBase_VERTEX,
		     &element_nodes,
		     &element_nodes_allocated,
		     &element_nodes_size,
		     &error );
    assert( iBase_SUCCESS == error );

    int coords_allocated = element_nodes_size*3;
    int coords_size = 0;
    double *coord_array = 0;
    iMesh_getVtxArrCoords( d_domain->getMesh(),
			   element_nodes,
			   element_nodes_size,
			   iBase_INTERLEAVED,
			   &coord_array,
			   &coords_allocated,
			   &coords_size,
			   &error );
    assert( iBase_SUCCESS == error );
    
    int n = 0;
    while ( n < element_nodes_size && !return_val )
    {
	if ( !( coord_array[3*n]   >= box[0] &&
		coord_array[3*n]   <= box[1] &&
		coord_array[3*n+1] >= box[2] &&
		coord_array[3*n+1] <= box[3] &&
		coord_array[3*n+2] >= box[4] &&
		coord_array[3*n+2] <= box[5] ) )
	{
	    return_val = false;
	}	
	++n;
    }

    free( element_nodes );
    free( coord_array );

    return return_val;
}

/*!
 * \brief Slice a box along the specified axis and return the resulting boxes.
 */
void Octree::sliceBox( const Box &parent_box,
		       Box &child_box1,
		       Box &child_box2,
		       const int cutting_axis )
{
    if ( cutting_axis == X_Axis )
    {
	child_box1[0] = parent_box[0];
	child_box1[1] = (parent_box[1] - parent_box[0]) / 2 + parent_box[0];
	child_box1[2] = parent_box[2];
	child_box1[3] = parent_box[3];
	child_box1[4] = parent_box[4];
	child_box1[5] = parent_box[5];

	child_box2[0] = (parent_box[1] - parent_box[0]) / 2 + parent_box[0];
	child_box2[1] = parent_box[1];
	child_box2[2] = parent_box[2];
	child_box2[3] = parent_box[3];
	child_box2[4] = parent_box[4];
	child_box2[5] = parent_box[5];
    }
    else if ( cutting_axis == Y_Axis )
    {
	child_box1[0] = parent_box[0];
	child_box1[1] = parent_box[1];
	child_box1[2] = parent_box[2];
	child_box1[3] = (parent_box[3] - parent_box[2]) / 2 + parent_box[2];
	child_box1[4] = parent_box[4];
	child_box1[5] = parent_box[5];

	child_box2[0] = parent_box[0];
	child_box2[1] = parent_box[1];
	child_box2[2] = (parent_box[3] - parent_box[2]) / 2 + parent_box[2];
	child_box2[3] = parent_box[3];
	child_box2[4] = parent_box[4];
	child_box2[5] = parent_box[5];
    }
    else
    {
	child_box1[0] = parent_box[0];
	child_box1[1] = parent_box[1];
	child_box1[2] = parent_box[2];
	child_box1[3] = parent_box[3];
	child_box1[4] = parent_box[4];
	child_box1[5] = (parent_box[5] - parent_box[4]) / 2 + parent_box[4];

	child_box2[0] = parent_box[0];
	child_box2[1] = parent_box[1];
	child_box2[2] = parent_box[2];
	child_box2[3] = parent_box[3];
	child_box2[4] = (parent_box[5] - parent_box[4]) / 2 + parent_box[4];
	child_box2[5] = parent_box[5];
    }
}


} // end namespace food

//---------------------------------------------------------------------------//
// end Octree.cpp
//---------------------------------------------------------------------------//
