//---------------------------------------------------------------------------//
// \file KDTree.cpp
// \author Stuart Slattery
// \brief KDTree declaration.
//---------------------------------------------------------------------------//

#include <cassert>
#include <algorithm>
#include <vector>
#include <iterator>

#include "KDTree.hpp"

namespace FOOD
{

/*! 
 * \brief Constructor.
 */
KDTree::KDTree( RCP_Domain domain, 
		const int entity_type,
		const int entity_topology )
    : d_domain(domain)
    , d_entity_type(entity_type)
    , d_entity_topology(entity_topology)
    , d_root_node( Teuchos::rcp(new KDTreeNode) )
{ /* ... */ }

/*!
 * \brief Destructor.
 */
KDTree::~KDTree()
{ /* ... */ }

/*! 
 * \brief Build the tree.
 */
void KDTree::buildTree()
{
    d_root_node->node_set = d_domain->getMeshSet();
    getEntSetBox( d_root_node->node_set, d_root_node->bounding_box );

    buildTreeNode( d_root_node );
}

/*!
 * \brief Locate a point. Return false if we didn't find it.
 */
bool KDTree::findPoint( iBase_EntityHandle &found_in_entity,
			const MDArray &coords )
{
    RCP_Node starting_node = findLeafNode( d_root_node, coords );
    return findPointInNode( starting_node, found_in_entity, coords );
}

/*! 
 * \brief Build a tree node.
 */
void KDTree::buildTreeNode( RCP_Node node )
{
    int error = 0;

    // Create the children.
    node->child1 = Teuchos::rcp( new KDTreeNode );
    node->child2 = Teuchos::rcp( new KDTreeNode );
    node->child1->parent = node;
    node->child2->parent = node;
    if ( node->axis == 0 )
    {
	node->child1->axis = 1;
	node->child2->axis = 1;
    }
    else if ( node->axis == 1 )
    {
	node->child1->axis = 2;
	node->child2->axis = 2;
    }

    iMesh_createEntSet( d_domain->getMesh(),
			1,
			&(node->child1->node_set),
			&error );
    assert( iBase_SUCCESS == error );

    iMesh_createEntSet( d_domain->getMesh(),
			1,
			&(node->child2->node_set),
			&error );
    assert( iBase_SUCCESS == error );

    // Create the new bounding boxes and assign them to the children.
    double balance = sliceBox( node );

    // Add the vertices in the parent set to the child sets.
    std::vector<iBase_EntityHandle> root_list;
    int node_vertices_allocated = 0;
    int node_vertices_size = 0;
    iBase_EntityHandle *node_vertices;
    iMesh_getEntities( d_domain->getMesh(),
		       node->node_set,
		       d_entity_type,
		       d_entity_topology,
		       &node_vertices,
		       &node_vertices_allocated,
		       &node_vertices_size,
		       &error );
    assert( iBase_SUCCESS == error );

    int node_coords_allocated = 3*node_vertices_size;
    int node_coords_size = 0;
    double *coords = 0;
    iMesh_getVtxArrCoords( d_domain->getMesh(),
			   node_vertices,
			   node_vertices_size,
			   iBase_BLOCKED,
			   &coords,
			   &node_coords_allocated,
			   &node_coords_size,
			   &error );
    assert( iBase_SUCCESS == error );

    bool vertex_found = false;
    for (int n = 0; n < node_vertices_size; ++n )
    {
	vertex_found = false;

	if ( coords[3*n + node->axis] < balance )
	{
	    iMesh_addEntToSet( d_domain->getMesh(),
			       node_vertices[n],
			       node->child1->node_set,
			       &error );
	    assert( iBase_SUCCESS == error );

	    vertex_found = true;
	}
	else if ( coords[3*n + node->axis] > balance )
	{
	    iMesh_addEntToSet( d_domain->getMesh(),
			       node_vertices[n],
			       node->child2->node_set,
			       &error );
	    assert( iBase_SUCCESS == error );

	    vertex_found = true;
	}

	if ( node != d_root_node && vertex_found )
	{
	    iMesh_rmvEntFromSet( d_domain->getMesh(),
				 node_vertices[n],
				 node->node_set,
				 &error );
	    assert( iBase_SUCCESS == error );
	}

	if ( !vertex_found && node == d_root_node )
	{
	    root_list.push_back( node_vertices[n] );
	}
    }

    free( node_vertices );
    free( coords );

    // Special case for the root node.
    if ( node == d_root_node )
    {
	iBase_EntitySetHandle new_root_set;
	iMesh_createEntSet( d_domain->getMesh(),
			    1,
			    &new_root_set,
			    &error );
	assert( iBase_SUCCESS == error );

	iMesh_addEntArrToSet( d_domain->getMesh(),
			      &root_list[0],
			      (int) root_list.size(),
			      new_root_set,
			      &error );
	assert( iBase_SUCCESS == error );

	node->node_set = new_root_set;
    }

    // See if we have any child entities. If more than 2, recurse.
    int total_child_ents = 0;
    int num_child_ents = 0;

    iMesh_getNumOfTopo( d_domain->getMesh(),
			node->child1->node_set,
			d_entity_topology,
			&num_child_ents,
			&error );
    assert( iBase_SUCCESS == error );
    total_child_ents += num_child_ents;
    if ( num_child_ents < 3 )
    {
	node->child1->is_leaf = true;
    }
    else
    {
	buildTreeNode( node->child1 );
    }

    iMesh_getNumOfTopo( d_domain->getMesh(),
			node->child2->node_set,
			d_entity_topology,
			&num_child_ents,
			&error );
    assert( iBase_SUCCESS == error );
    total_child_ents += num_child_ents;
    if ( num_child_ents < 3 )
    {
	node->child2->is_leaf = true;
    }
    else
    {
	buildTreeNode( node->child2 );
    }
    
    if ( total_child_ents == 0 )
    {
	node->is_leaf = true;
    }
}

/*!
 * \brief Given a point, find its leaf node in the tree.
 */
KDTree::RCP_Node KDTree::findLeafNode( RCP_Node node, const MDArray &coords )
{
    RCP_Node return_node;
    if ( !node->is_leaf )
    {
	if ( isPointInBox( node->child1->bounding_box, coords ) )
	{
	    return_node = findLeafNode( node->child1, coords );
	}
	else if ( isPointInBox( node->child2->bounding_box, coords ) )
	{
	    return_node = findLeafNode( node->child2, coords );
	}
    }
    return return_node;
}

/*!
 * \brief Search a node for a point. Return the element we found it in.
 */
bool KDTree::findPointInNode( RCP_Node node,
			      iBase_EntityHandle &found_in_entity,
			      const MDArray &coords )
{
    int error = 0;
    bool return_val = false;

    // First check at the local level.
    int node_vertices_allocated = 0;
    int node_vertices_size = 0;
    iBase_EntityHandle *node_vertices;
    iMesh_getEntities( d_domain->getMesh(),
		       node->node_set,
		       d_entity_type,
		       d_entity_topology,
		       &node_vertices,
		       &node_vertices_allocated,
		       &node_vertices_size,
		       &error );
    assert( iBase_SUCCESS == error );

    iBase_EntityHandle *node_adj_elements = 0;
    int node_adj_elements_allocated = 0; 
    int node_adj_elements_size = 0; 
    int *offset = 0;
    int offset_allocated = 0; 
    int offset_size = 0; 
    iMesh_getEntArrAdj( d_domain->getMesh(),
			node_vertices,
			d_entity_type,
			d_entity_topology,
			&node_adj_elements,
			&node_adj_elements_allocated,
			&node_adj_elements_size,
			&offset,
			&offset_allocated,
			&offset_size,
			&error );
    assert( iBase_SUCCESS == error );

    int i = 0;
    if ( node_adj_elements_size > 0 )
    {
	while ( i < node_adj_elements_size && !return_val )
	{
	    if ( PointQuery::point_in_ref_element( d_domain->getMesh(),
						   node_adj_elements[i],
						   coords ) )
	    {
		return_val = true;
		found_in_entity = node_adj_elements[i];
	    }
	    ++i;
	}
    }

    free( node_adj_elements );
    free( node_vertices );

    // If we found an element then we're done, otherwise recurse back up the
    // tree until we hit the root node. At that point we're done.
    if ( !return_val && node != d_root_node )
    {
	return_val = findPointInNode( node->parent,
				      found_in_entity,
				      coords );
    }

    return return_val;
}

/*!
 * \brief Get the bounding box of a set of entities.
 */
void KDTree::getEntSetBox( iBase_EntitySetHandle entity_set,
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
bool KDTree::isPointInBox( const Box &box,
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
bool KDTree::isEntInBox( const Box &box,
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
double KDTree::sliceBox( RCP_Node node )
{
    int error = 0;

    // Get the node vertices.
    int node_vertices_allocated = 0;
    int node_vertices_size = 0;
    iBase_EntityHandle *node_vertices;
    iMesh_getEntities( d_domain->getMesh(),
		       node->node_set,
		       d_entity_type,
		       d_entity_topology,
		       &node_vertices,
		       &node_vertices_allocated,
		       &node_vertices_size,
		       &error );
    assert( iBase_SUCCESS == error );

    int node_coords_allocated = 3*node_vertices_size;
    int node_coords_size = 0;
    double *coords = 0;
    iMesh_getVtxArrCoords( d_domain->getMesh(),
			   node_vertices,
			   node_vertices_size,
			   iBase_BLOCKED,
			   &coords,
			   &node_coords_allocated,
			   &node_coords_size,
			   &error );
    assert( iBase_SUCCESS == error );

    // Find the median point in the cutting node->axis and make the new boxes.
    Teuchos::ArrayView<double> node_coords( coords, node_coords_size );
    Teuchos::ArrayView<double>::const_iterator node_it_begin;
    Teuchos::ArrayView<double>::const_iterator node_it_end;
    double cutting_val = 0.0;
    if ( node->axis == 0 )
    {
	node_it_begin = node_coords.begin();
	node_it_end = node_coords.begin() + node_vertices_size;
	
	cutting_val = median( node_it_begin, node_it_end );

	node->child1->bounding_box[0] = node->bounding_box[0];
	node->child1->bounding_box[1] = node->bounding_box[0] + cutting_val;
	node->child1->bounding_box[2] = node->bounding_box[2];
	node->child1->bounding_box[3] = node->bounding_box[3];
	node->child1->bounding_box[4] = node->bounding_box[4];
	node->child1->bounding_box[5] = node->bounding_box[5];

	node->child2->bounding_box[0] = node->bounding_box[0] + cutting_val;
	node->child2->bounding_box[1] = node->bounding_box[1];
	node->child2->bounding_box[2] = node->bounding_box[2];
	node->child2->bounding_box[3] = node->bounding_box[3];
	node->child2->bounding_box[4] = node->bounding_box[4];
	node->child2->bounding_box[5] = node->bounding_box[5];
    }
    else if ( node->axis == 1 )
    {
	node_it_begin = node_coords.begin() + node_vertices_size;
	node_it_end = node_coords.begin() + 2*node_vertices_size;

	cutting_val = median( node_it_begin, node_it_end );

	node->child1->bounding_box[0] = node->bounding_box[0];
	node->child1->bounding_box[1] = node->bounding_box[1];
	node->child1->bounding_box[2] = node->bounding_box[2];
	node->child1->bounding_box[3] = node->bounding_box[2] + cutting_val;
	node->child1->bounding_box[4] = node->bounding_box[4];
	node->child1->bounding_box[5] = node->bounding_box[5];

	node->child2->bounding_box[0] = node->bounding_box[0];
	node->child2->bounding_box[1] = node->bounding_box[1];
	node->child2->bounding_box[2] = node->bounding_box[2] + cutting_val;
	node->child2->bounding_box[3] = node->bounding_box[3];
	node->child2->bounding_box[4] = node->bounding_box[4];
	node->child2->bounding_box[5] = node->bounding_box[5];
    }
    else 
    {
	node_it_begin = node_coords.begin() + 2*node_vertices_size;
	node_it_end = node_coords.end();

	cutting_val = median( node_it_begin, node_it_end );

	node->child1->bounding_box[0] = node->bounding_box[0];
	node->child1->bounding_box[1] = node->bounding_box[1];
	node->child1->bounding_box[2] = node->bounding_box[2];
	node->child1->bounding_box[3] = node->bounding_box[3];
	node->child1->bounding_box[4] = node->bounding_box[4];
	node->child1->bounding_box[5] = node->bounding_box[4] + cutting_val;

	node->child2->bounding_box[0] = node->bounding_box[0];
	node->child2->bounding_box[1] = node->bounding_box[1];
	node->child2->bounding_box[2] = node->bounding_box[2];
	node->child2->bounding_box[3] = node->bounding_box[3];
	node->child2->bounding_box[4] = node->bounding_box[4] + cutting_val;
	node->child2->bounding_box[5] = node->bounding_box[5];
    }

    free( node_vertices );
    free( coords );
    return cutting_val;
}

/*
 * \brief Compute the median of a range of values.
 */
double KDTree::median( Teuchos::ArrayView<double>::const_iterator begin,
		       Teuchos::ArrayView<double>::const_iterator end )
{
    std::size_t middle = std::distance( begin , end ) / 2;
    std::nth_element( begin, begin+middle, end );
    return *(begin + middle);
}

} // end namespace food

//---------------------------------------------------------------------------//
// end KDTree.cpp
//---------------------------------------------------------------------------//
