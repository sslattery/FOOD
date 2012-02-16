//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file Octree.cpp
 * \author Stuart Slattery
 * \brief Octree defintion.
 */
//---------------------------------------------------------------------------//

#include <cassert>
#include <algorithm>
#include <vector>

#include "Octree.hpp"

namespace FOOD
{

/*! 
 * \brief Constructor.
 */
Octree::Octree( iMesh_Instance mesh,
		iBase_EntitySetHandle mesh_set, 
		const int entity_type,
		const int entity_topology )
    : d_mesh(mesh)
    , d_mesh_set(mesh_set)
    , d_entity_type(entity_type)
    , d_entity_topology(entity_topology)
    , d_root_node( Teuchos::rcp(new OctreeNode) )
{ /* ... */ }

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
    d_root_node->node_set = d_mesh_set;
    getEntSetBox( d_root_node->node_set, d_root_node->bounding_box );

    buildTreeNode( d_root_node );
}

/*!
 * \brief Locate a point. Return false if we didn't find it.
 */
bool Octree::findPoint( iBase_EntityHandle &found_in_entity,
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

    // Create the children.
    for ( int i = 0; i < 8; ++i )
    {
	node->children[i] = Teuchos::rcp( new OctreeNode );

	iMesh_createEntSet( d_mesh,
			    1,
			    &(node->children[i]->node_set),
			    &error );
	assert( iBase_SUCCESS == error );
    }

    // Create the new bounding boxes and assign them to the children.
    sliceBox( node );

    // Add the elements in the parent set to the child sets.
    std::vector<iBase_EntityHandle> root_list;
    int node_elements_allocated = 0;
    int node_elements_size = 0;
    iBase_EntityHandle *node_elements;
    iMesh_getEntities( d_mesh,
		       node->node_set,
		       d_entity_type,
		       d_entity_topology,
		       &node_elements,
		       &node_elements_allocated,
		       &node_elements_size,
		       &error );
    assert( iBase_SUCCESS == error );

    bool element_found = false;
    for (int n = 0; n < node_elements_size; ++n )
    {
	element_found = false;

	for ( int m = 0; m < 8; ++m )
	{
	    if ( !element_found )
	    {
		if ( isEntInBox( node->children[m]->bounding_box, 
				 node_elements[n] ) )
		{
		    iMesh_addEntToSet( d_mesh,
				       node_elements[n],
				       node->children[m]->node_set,
				       &error );
		    assert( iBase_SUCCESS == error );

		    if ( node != d_root_node )
		    {
			iMesh_rmvEntFromSet( d_mesh,
					     node_elements[n],
					     node->node_set,
					     &error );
			assert( iBase_SUCCESS == error );
		    }

		    element_found = true;
		}
	    }

	}

	if ( !element_found )
	{
	    root_list.push_back( node_elements[n] );
	}
    }

    free( node_elements );

    // Special case for the root node.
    if ( node == d_root_node )
    {
	iBase_EntitySetHandle new_root_set;
	iMesh_createEntSet( d_mesh,
			    1,
			    &new_root_set,
			    &error );
	assert( iBase_SUCCESS == error );

	iMesh_addEntArrToSet( d_mesh,
			      &root_list[0],
			      (int) root_list.size(),
			      new_root_set,
			      &error );
	assert( iBase_SUCCESS == error );

	node->node_set = new_root_set;
    }

    // See if we have any child entities.
    Teuchos::Tuple<int,8> num_child_ents;
    int total_child_ents = 0;
    for (int i = 0; i < 8; ++i )
    {
	iMesh_getNumOfTopo( d_mesh,
			    node->children[i]->node_set,
			    d_entity_topology,
			    &num_child_ents[i],
			    &error );
	assert( iBase_SUCCESS == error );

	total_child_ents += num_child_ents[i];
    }
    
    // Recurse.
    if ( total_child_ents == 0 )
    {
	node->is_leaf = true;
    }
    else
    {
	for (int i = 0; i < 8; ++i)
	{
	    if ( num_child_ents[i] == 0 )
	    {
		node->children[i]->is_leaf = true;
	    }
	    else 
	    {
		buildTreeNode( node->children[i] );
	    }
	}
    }
}

/*!
 * \brief Search a node for a point. If its a leaf node, return the element we
 * found it in.
 */
bool Octree::findPointInNode( RCP_Node node,
			      iBase_EntityHandle &found_in_entity,
			      const MDArray &coords )
{
    int error = 0;
    bool return_val = false;

    // First check at the node level.
    iBase_EntityHandle *node_elements = 0;
    int node_elements_allocated = 0; 
    int node_elements_size = 0; 
    iMesh_getEntities( d_mesh,
		       node->node_set,
		       d_entity_type,
		       d_entity_topology,
		       &node_elements,
		       &node_elements_allocated,
		       &node_elements_size,
		       &error );
    assert( iBase_SUCCESS == error );

    int i = 0;
    if ( node_elements_size > 0 )
    {
	while ( i < node_elements_size && !return_val )
	{
	    if ( PointQuery::pointInRefElement( d_mesh,
						node_elements[i],
						coords ) )
	    {
		return_val = true;
		found_in_entity = node_elements[i];
	    }
	    ++i;
	}
    }

    free( node_elements );

    // If we found the element or we're at a leaf then we're done. Otherwise,
    // recurse through the children if this point is in their bounding box.
    if ( !return_val && !node->is_leaf )
    {
	int j = 0;
	while ( j < 8 && !return_val )
	{
	    if ( isPointInBox( node->children[j]->bounding_box, coords ) )
	    {
		return_val = findPointInNode( node->children[j],
					      found_in_entity,
					      coords );
	    }
	    ++j;
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
    iMesh_getNumOfTopo( d_mesh,
			entity_set,
			iMesh_POINT,
			&num_set_vertices,
			&error );
    assert( iBase_SUCCESS == error );

    iBase_EntityHandle *set_vertices = 0;
    int set_vertices_allocated = 0;
    int set_vertices_size = 0;
    iMesh_getEntities( d_mesh,
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
    iMesh_getVtxArrCoords( d_mesh,
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
    iMesh_getEntAdj( d_mesh,
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
    iMesh_getVtxArrCoords( d_mesh,
			   element_nodes,
			   element_nodes_size,
			   iBase_INTERLEAVED,
			   &coord_array,
			   &coords_allocated,
			   &coords_size,
			   &error );
    assert( iBase_SUCCESS == error );
    
    int n = 0;
    while ( n < element_nodes_size && return_val )
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
void Octree::sliceBox( RCP_Node node )
{
    // child 0
    node->children[0]->bounding_box[0] = node->bounding_box[0];
    node->children[0]->bounding_box[1] = (node->bounding_box[1] - node->bounding_box[0]) / 2 + node->bounding_box[0];
    node->children[0]->bounding_box[2] = node->bounding_box[2];
    node->children[0]->bounding_box[3] = (node->bounding_box[3] - node->bounding_box[2]) / 2 + node->bounding_box[2];
    node->children[0]->bounding_box[4] = node->bounding_box[4];
    node->children[0]->bounding_box[5] = (node->bounding_box[5] - node->bounding_box[4]) / 2 + node->bounding_box[4];

    // child 1
    node->children[1]->bounding_box[0] = (node->bounding_box[1] - node->bounding_box[0]) / 2 + node->bounding_box[0];
    node->children[1]->bounding_box[1] = node->bounding_box[1];
    node->children[1]->bounding_box[2] = node->bounding_box[2];
    node->children[1]->bounding_box[3] = (node->bounding_box[3] - node->bounding_box[2]) / 2 + node->bounding_box[2];
    node->children[1]->bounding_box[4] = node->bounding_box[4];
    node->children[1]->bounding_box[5] = (node->bounding_box[5] - node->bounding_box[4]) / 2 + node->bounding_box[4];

    // child 2
    node->children[2]->bounding_box[0] = (node->bounding_box[1] - node->bounding_box[0]) / 2 + node->bounding_box[0];
    node->children[2]->bounding_box[1] = node->bounding_box[1];
    node->children[2]->bounding_box[2] = (node->bounding_box[3] - node->bounding_box[2]) / 2 + node->bounding_box[2];
    node->children[2]->bounding_box[3] = node->bounding_box[3];
    node->children[2]->bounding_box[4] = node->bounding_box[4];
    node->children[2]->bounding_box[5] = (node->bounding_box[5] - node->bounding_box[4]) / 2 + node->bounding_box[4];

    // child 3
    node->children[3]->bounding_box[0] = node->bounding_box[0];
    node->children[3]->bounding_box[1] = (node->bounding_box[1] - node->bounding_box[0]) / 2 + node->bounding_box[0];
    node->children[3]->bounding_box[2] = (node->bounding_box[3] - node->bounding_box[2]) / 2 + node->bounding_box[2];
    node->children[3]->bounding_box[3] = node->bounding_box[3];
    node->children[3]->bounding_box[4] = node->bounding_box[4];
    node->children[3]->bounding_box[5] = (node->bounding_box[5] - node->bounding_box[4]) / 2 + node->bounding_box[4];

    // child 4
    node->children[4]->bounding_box[0] = node->bounding_box[0];
    node->children[4]->bounding_box[1] = (node->bounding_box[1] - node->bounding_box[0]) / 2 + node->bounding_box[0];
    node->children[4]->bounding_box[2] = node->bounding_box[2];
    node->children[4]->bounding_box[3] = (node->bounding_box[3] - node->bounding_box[2]) / 2 + node->bounding_box[2];
    node->children[4]->bounding_box[4] = (node->bounding_box[5] - node->bounding_box[4]) / 2 + node->bounding_box[4];
    node->children[4]->bounding_box[5] = node->bounding_box[5];

    // child 5
    node->children[5]->bounding_box[0] = (node->bounding_box[1] - node->bounding_box[0]) / 2 + node->bounding_box[0];
    node->children[5]->bounding_box[1] = node->bounding_box[1];
    node->children[5]->bounding_box[2] = node->bounding_box[2];
    node->children[5]->bounding_box[3] = (node->bounding_box[3] - node->bounding_box[2]) / 2 + node->bounding_box[2];
    node->children[5]->bounding_box[4] = (node->bounding_box[5] - node->bounding_box[4]) / 2 + node->bounding_box[4];
    node->children[5]->bounding_box[5] = node->bounding_box[5];

    // child 6
    node->children[6]->bounding_box[0] = (node->bounding_box[1] - node->bounding_box[0]) / 2 + node->bounding_box[0];
    node->children[6]->bounding_box[1] = node->bounding_box[1];
    node->children[6]->bounding_box[2] = (node->bounding_box[3] - node->bounding_box[2]) / 2 + node->bounding_box[2];
    node->children[6]->bounding_box[3] = node->bounding_box[3];
    node->children[6]->bounding_box[4] = (node->bounding_box[5] - node->bounding_box[4]) / 2 + node->bounding_box[4];
    node->children[6]->bounding_box[5] = node->bounding_box[5];

    // child 7
    node->children[7]->bounding_box[0] = node->bounding_box[0];
    node->children[7]->bounding_box[1] = (node->bounding_box[1] - node->bounding_box[0]) / 2 + node->bounding_box[0];
    node->children[7]->bounding_box[2] = (node->bounding_box[3] - node->bounding_box[2]) / 2 + node->bounding_box[2];
    node->children[7]->bounding_box[3] = node->bounding_box[3];
    node->children[7]->bounding_box[4] = (node->bounding_box[5] - node->bounding_box[4]) / 2 + node->bounding_box[4];
    node->children[7]->bounding_box[5] = node->bounding_box[5];
}


} // end namespace food

//---------------------------------------------------------------------------//
// end Octree.cpp
//---------------------------------------------------------------------------//
