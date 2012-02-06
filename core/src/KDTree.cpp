//---------------------------------------------------------------------------//
// \file KDTree.cpp
// \author Stuart Slattery
// \brief KDTree definition.
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
    , d_root_node( new KDTreeNode )
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
    bool return_val = false;

    if ( isPointInBox( d_root_node->bounding_box, coords ) )
    {
	iBase_EntityHandle nearest_neighbor;
	KDTreeNode* starting_node = findLeafNode( d_root_node, 
						  nearest_neighbor, 
						  coords );

	return_val = findPointInNode( starting_node, 
				      nearest_neighbor, 
				      coords,
				      found_in_entity );
    }

    return return_val;
}

/*! 
 * \brief Build a tree node.
 */
void KDTree::buildTreeNode( KDTreeNode* node )
{
    int error = 0;

    // Create the children.
    node->child1 = new KDTreeNode;
    node->child2 = new KDTreeNode;
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
    sliceBox( node );

    // Add the vertices in the parent set to the child sets.
    std::vector<iBase_EntityHandle> root_list;
    int node_vertices_allocated = 0;
    int node_vertices_size = 0;
    iBase_EntityHandle *node_vertices;
    iMesh_getEntities( d_domain->getMesh(),
		       node->node_set,
		       iBase_VERTEX,
		       iMesh_POINT,
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
			   iBase_INTERLEAVED,
			   &coords,
			   &node_coords_allocated,
			   &node_coords_size,
			   &error );
    assert( iBase_SUCCESS == error );

    bool vertex_found = false;
    for (int n = 0; n < node_vertices_size; ++n )
    {
	vertex_found = false;

	if ( coords[3*n + node->axis] < node->median )
	{
	    iMesh_addEntToSet( d_domain->getMesh(),
			       node_vertices[n],
			       node->child1->node_set,
			       &error );
	    assert( iBase_SUCCESS == error );

	    vertex_found = true;
	}
	else if ( coords[3*n + node->axis] > node->median )
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

    // See if we have any child vertices. If more than 2, recurse.
    int total_child_verts = 0;
    int num_child_verts = 0;

    iMesh_getNumOfTopo( d_domain->getMesh(),
			node->child1->node_set,
			iMesh_POINT,
			&num_child_verts,
			&error );
    assert( iBase_SUCCESS == error );
    total_child_verts += num_child_verts;
    if ( num_child_verts < 3 )
    {
	node->child1->is_leaf = true;
    }
    else
    {
	buildTreeNode( node->child1 );
    }

    iMesh_getNumOfTopo( d_domain->getMesh(),
			node->child2->node_set,
			iMesh_POINT,
			&num_child_verts,
			&error );
    assert( iBase_SUCCESS == error );
    total_child_verts += num_child_verts;
    if ( num_child_verts < 3 )
    {
	node->child2->is_leaf = true;
    }
    else
    {
	buildTreeNode( node->child2 );
    }
    
    if ( total_child_verts == 0 )
    {
	node->is_leaf = true;
    }
}

/*!
 * \brief Given a point, find its leaf node in the tree.
 */
KDTreeNode* KDTree::findLeafNode( KDTreeNode* node, 
				  iBase_EntityHandle &nearest_neighbor,
				  const MDArray &coords )
{
    int error = 0;

    KDTreeNode* return_node;
    if ( !node->is_leaf )
    {
	if ( coords(0,node->axis) < node->median )
	{
	    return_node = findLeafNode( node->child1, nearest_neighbor, coords );
	}
	else if ( coords(0,node->axis) > node->median )
	{
	    return_node = findLeafNode( node->child2, nearest_neighbor, coords );
	}
	else
	{
	    return_node = node;
	}
    }
    else
    {
	return_node = node;
    }

    if ( return_node == node )
    {
	double distance = 1.0e99;
	int node_vertices_allocated = 0;
	int node_vertices_size = 0;
	iBase_EntityHandle *node_vertices;
	iMesh_getEntities( d_domain->getMesh(),
			   node->node_set,
			   iBase_VERTEX,
			   iMesh_POINT,
			   &node_vertices,
			   &node_vertices_allocated,
			   &node_vertices_size,
			   &error );
	assert( iBase_SUCCESS == error );

	if ( node_vertices_size > 0 )
	{
	    double local_distance = 0;
	    double x = 0;
	    double y = 0;
	    double z = 0;
	    for ( int i = 0; i < node_vertices_size; ++i )
	    {
		iMesh_getVtxCoord( d_domain->getMesh(),
				   node_vertices[i],
				   &x,
				   &y,
				   &z,
				   &error );
		assert( iBase_SUCCESS == error );

		local_distance = (coords(0,0) - x)*(coords(0,0) - x) +
				 (coords(0,1) - y)*(coords(0,1) - y) +
				 (coords(0,2) - z)*(coords(0,2) - z);

		if ( local_distance < distance )
		{
		    distance = local_distance;
		    nearest_neighbor = node_vertices[i];
		}
	    }
	}

	free( node_vertices );
    }

    return return_node;
}

/*!
 * \brief Search a node for a point. Return the element we found it in.
 */
bool KDTree::findPointInNode( KDTreeNode* node,
			      iBase_EntityHandle &nearest_neighbor,
			      const MDArray &coords,
			      iBase_EntityHandle &found_in_entity )
{
    int error = 0;
    bool entity_found = false;

    // First check at the local level.
    int node_vertices_allocated = 0;
    int node_vertices_size = 0;
    iBase_EntityHandle *node_vertices;
    iMesh_getEntities( d_domain->getMesh(),
		       node->node_set,
		       iBase_VERTEX,
		       iMesh_POINT,
		       &node_vertices,
		       &node_vertices_allocated,
		       &node_vertices_size,
		       &error );
    assert( iBase_SUCCESS == error );

    double x = 0;
    double y = 0;
    double z = 0;
    iMesh_getVtxCoord( d_domain->getMesh(),
		       nearest_neighbor,
		       &x,
		       &y,
		       &z,
		       &error );
    assert( iBase_SUCCESS == error );

    double distance = (coords(0,0) - x)*(coords(0,0) - x) +
		      (coords(0,1) - y)*(coords(0,1) - y) +
		      (coords(0,2) - z)*(coords(0,2) - z);
	
    if ( node_vertices_size > 0 )
    {
	double local_distance = 0;
	for ( int i = 0; i < node_vertices_size; ++i )
	{
	    entity_found = isPointInAdj( node_vertices[i],
					 coords,
					 found_in_entity );
	    
	    if ( !entity_found )
	    {
		iMesh_getVtxCoord( d_domain->getMesh(),
				   node_vertices[i],
				   &x,
				   &y,
				   &z,
				   &error );
		assert( iBase_SUCCESS == error );

		local_distance = (coords(0,0) - x)*(coords(0,0) - x) +
				 (coords(0,1) - y)*(coords(0,1) - y) +
				 (coords(0,2) - z)*(coords(0,2) - z);

		if ( local_distance < distance )
		{
		    distance = local_distance;
		    nearest_neighbor = node_vertices[i];
		}
	    }
	}
    }

    free( node_vertices );

    // Check to see if there could be points on the other side of the
    // splitting plane. If we're at the root node then we're done
    if ( node != d_root_node && !entity_found)
    {
	double distance_to_plane = (node->parent->median - coords(0,node->parent->axis))*
				   (node->parent->median - coords(0,node->parent->axis));
	if ( distance > distance_to_plane )
	{
	    if ( node == node->parent->child1 )
	    {
		findPointInNode( node->parent->child2, 
				 nearest_neighbor, 
				 coords,
				 found_in_entity );
	    }
	    else 
	    {
		findPointInNode( node->parent->child1, 
				 nearest_neighbor, 
				 coords,
				 found_in_entity);
	    }
	}
    
	findPointInNode( node->parent, nearest_neighbor, coords, found_in_entity );
    }

    return entity_found;
}

/*!
 * \brief Get the bounding box of a set of entities.
 */
void KDTree::getEntSetBox( iBase_EntitySetHandle entity_set,
			   double bounding_box[6] )
{
    double grow = 1.0;
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

    bounding_box[0] = *(std::min_element( set_x_it_begin, set_x_it_end ))-grow;
    bounding_box[1] = *(std::max_element( set_x_it_begin, set_x_it_end ))+grow;
    bounding_box[2] = *(std::min_element( set_y_it_begin, set_y_it_end ))-grow;
    bounding_box[3] = *(std::max_element( set_y_it_begin, set_y_it_end ))+grow;
    bounding_box[4] = *(std::min_element( set_z_it_begin, set_z_it_end ))-grow;
    bounding_box[5] = *(std::max_element( set_z_it_begin, set_z_it_end ))+grow;

    free( set_vertices );
    free( coords );
}

/*!
 * \brief Determine if a point is inside the entities adjacent to another
 * point. 
 */
bool KDTree::isPointInAdj( iBase_EntityHandle point, 
			   const MDArray &coords,
			   iBase_EntityHandle &found_in_entity )
{
    bool return_val = false;
    int error = 0;

    iBase_EntityHandle *adj_elements = 0;
    int adj_elements_allocated = 0; 
    int adj_elements_size = 0; 
    iMesh_getEntAdj( d_domain->getMesh(),
    		     point,
    		     d_entity_type,
    		     &adj_elements,
    		     &adj_elements_allocated,
    		     &adj_elements_size,
    		     &error );
    assert( iBase_SUCCESS == error );

    int i = 0;
    if ( adj_elements_size > 0 )
    {
    	while ( i < adj_elements_size && !return_val )
    	{
    	    if ( PointQuery::point_in_ref_element( d_domain->getMesh(),
    						   adj_elements[i],
    						   coords ) )
    	    {
    		return_val = true;
    		found_in_entity = adj_elements[i];
    	    }
    	    ++i;
    	}
    }
 
    free( adj_elements );

    return return_val;
}

/*!
 * \brief Determine if a point is inside a bounding box.
 */
bool KDTree::isPointInBox( const double box[6],
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
bool KDTree::isEntInBox( const double box[6],
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
void KDTree::sliceBox( KDTreeNode* node )
{
    int error = 0;

    // Get the node vertices.
    int node_vertices_allocated = 0;
    int node_vertices_size = 0;
    iBase_EntityHandle *node_vertices;
    iMesh_getEntities( d_domain->getMesh(),
		       node->node_set,
		       iBase_VERTEX,
		       iMesh_POINT,
		       &node_vertices,
		       &node_vertices_allocated,
		       &node_vertices_size,
		       &error );
    assert( iBase_SUCCESS == error );

    int node_coords_allocated = 0;
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
    Teuchos::ArrayView<double>::iterator node_it_begin;
    Teuchos::ArrayView<double>::iterator node_it_end;
    if ( node->axis == 0 )
    {
	node_it_begin = node_coords.begin();
	node_it_end = node_coords.begin() + node_vertices_size;
	
	node->median = median( node_it_begin, node_it_end );

	node->child1->bounding_box[0] = node->bounding_box[0];
	node->child1->bounding_box[1] = node->bounding_box[0] + node->median;
	node->child1->bounding_box[2] = node->bounding_box[2];
	node->child1->bounding_box[3] = node->bounding_box[3];
	node->child1->bounding_box[4] = node->bounding_box[4];
	node->child1->bounding_box[5] = node->bounding_box[5];

	node->child2->bounding_box[0] = node->bounding_box[0] + node->median;
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

	node->median = median( node_it_begin, node_it_end );

	node->child1->bounding_box[0] = node->bounding_box[0];
	node->child1->bounding_box[1] = node->bounding_box[1];
	node->child1->bounding_box[2] = node->bounding_box[2];
	node->child1->bounding_box[3] = node->bounding_box[2] + node->median;
	node->child1->bounding_box[4] = node->bounding_box[4];
	node->child1->bounding_box[5] = node->bounding_box[5];

	node->child2->bounding_box[0] = node->bounding_box[0];
	node->child2->bounding_box[1] = node->bounding_box[1];
	node->child2->bounding_box[2] = node->bounding_box[2] + node->median;
	node->child2->bounding_box[3] = node->bounding_box[3];
	node->child2->bounding_box[4] = node->bounding_box[4];
	node->child2->bounding_box[5] = node->bounding_box[5];
    }
    else 
    {
	node_it_begin = node_coords.begin() + 2*node_vertices_size;
	node_it_end = node_coords.end();

	node->median = median( node_it_begin, node_it_end );

	node->child1->bounding_box[0] = node->bounding_box[0];
	node->child1->bounding_box[1] = node->bounding_box[1];
	node->child1->bounding_box[2] = node->bounding_box[2];
	node->child1->bounding_box[3] = node->bounding_box[3];
	node->child1->bounding_box[4] = node->bounding_box[4];
	node->child1->bounding_box[5] = node->bounding_box[4] + node->median;

	node->child2->bounding_box[0] = node->bounding_box[0];
	node->child2->bounding_box[1] = node->bounding_box[1];
	node->child2->bounding_box[2] = node->bounding_box[2];
	node->child2->bounding_box[3] = node->bounding_box[3];
	node->child2->bounding_box[4] = node->bounding_box[4] + node->median;
	node->child2->bounding_box[5] = node->bounding_box[5];
    }

    // Cleanup.
    free( node_vertices );
    free( coords );
}

/*
 * \brief Compute the median of a range of values.
 */
double KDTree::median( Teuchos::ArrayView<double>::iterator begin,
		       Teuchos::ArrayView<double>::iterator end )
{
    std::size_t middle = std::distance( begin , end ) / 2;
    std::nth_element( begin, begin+middle, end );
    return *(begin + middle);
}

} // end namespace food

//---------------------------------------------------------------------------//
// end KDTree.cpp
//---------------------------------------------------------------------------//
