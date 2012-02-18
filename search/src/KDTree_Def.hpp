//---------------------------------------------------------------------------//
/*! 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file KDTree_Def.hpp
 * \author Stuart Slattery
 * \brief KDTree definition.
 */
//---------------------------------------------------------------------------//

#ifndef FOOD_KDTREE_DEF_HPP
#define FOOD_KDTREE_DEF_HPP

#include <cassert>
#include <cmath>
#include <algorithm>
#include <vector>
#include <iterator>

#include "KDTree.hpp"
#include "TopologyTools.hpp"
#include "PointQuery.hpp"

namespace FOOD
{

/*! 
 * \brief Constructor.
 */
template<int DIM>
KDTree<DIM>::KDTree( iMesh_Instance mesh, 
		     iBase_EntitySetHandle mesh_set, 
		     const int entity_type,
		     const int entity_topology )
    : d_mesh(mesh)
    , d_mesh_set(mesh_set)
    , d_entity_type(entity_type)
    , d_entity_topology(entity_topology)
    , d_large(HUGE_VAL)
    , d_points(0)
    , d_ptindx(0)
    , d_rptindx(0)
    , d_nodes(0)
{ /* ... */ }

/*!
 * \brief Destructor.
 */
template<int DIM>
KDTree<DIM>::~KDTree()
{ /* ... */ }

/*! 
 * \brief Build the tree.
 */
template<int DIM>
void KDTree<DIM>::buildTree()
{
    int error = 0;

    // Get the domain linear element vertices and their coordinates.
    iBase_EntityHandle *elements = 0;
    int elements_allocated = 0;
    int elements_size = 0;
    iMesh_getEntities( d_mesh,
		       d_mesh_set,
		       d_entity_type,
		       d_entity_topology,
		       &elements,
		       &elements_allocated,
		       &elements_size,
		       &error );
    assert( iBase_SUCCESS == error );

    for ( int n = 0; n < elements_size; ++n )
    {
	iBase_EntityHandle *element_nodes = 0;
	int element_nodes_allocated = 0;
	int element_nodes_size = 0;
	iMesh_getEntAdj( d_mesh,
			 elements[n],
			 iBase_VERTEX,
			 &element_nodes,
			 &element_nodes_allocated,
			 &element_nodes_size,
			 &error );
	assert( iBase_SUCCESS == error );

	int num_linear_nodes = 
	    TopologyTools::numLinearNodes( d_entity_topology );
	for ( int i = 0; i < num_linear_nodes; ++i )
	{
	    d_points.push_back( element_nodes[i] );
	}

	free( element_nodes );
    }

    free( elements );

    int coords_allocated = 0;
    int coords_size = 0;
    double *coords;
    iMesh_getVtxArrCoords( d_mesh,
			   &d_points[0],
			   (int) d_points.size(),
			   iBase_BLOCKED,
			   &coords,
			   &coords_allocated,
			   &coords_size,
			   &error );
    assert( iBase_SUCCESS == error );

    // Setup the point indices.
    d_ptindx.resize( d_points.size() );
    d_rptindx.resize( d_points.size() );
    for ( int k = 0; k < (int) d_points.size() ; ++k )
    {
	d_ptindx[k] = k;
    }

    // Compute the number of tree nodes.
    int m = 1;
    for ( int n = (int) d_points.size() ; n; n >>=1 )
    {
	m <<= 1;
    }
    int num_nodes = 2*d_points.size() - (m >> 1);
    if ( m < num_nodes )
    {
	num_nodes = m;
    }
    num_nodes--;
    d_nodes.resize(num_nodes);

    // Create the root node.
    Point<DIM> lo(-d_large, -d_large, -d_large);
    Point<DIM> hi(d_large, d_large, d_large);

    d_nodes[0] = Teuchos::rcp( 
	new KDTreeNode<DIM>( lo, hi, 0, 0, 0, 0, d_points.size()-1 ) );

    // Setup task lists for up to 2^50 points.
    std::vector<int> taskparent(50);
    std::vector<int> taskdim(50);

    // Partition the tree.
    int jbox = 0;
    taskparent[1] = 0;
    taskdim[1] = 0;
    int nowtask = 1;
    int tparent, tdim, ptlo, pthi, np, kk;
    int *hp;
    double *cp;
    while ( nowtask )
    {
	tparent = taskparent[nowtask];
	tdim = taskdim[nowtask--];
	ptlo = d_nodes[tparent]->ptlo;
	pthi = d_nodes[tparent]->pthi;
	hp = &d_ptindx[ptlo];
	cp = &coords[tdim*d_points.size()];
	np = pthi - ptlo + 1;
	kk = (np-1) / 2;

	(void) partition( kk, hp, np, cp );

	// Create the children.
	hi = d_nodes[tparent]->hi;
	lo = d_nodes[tparent]->lo;
	hi.x[tdim] = coords[tdim*d_points.size() + hp[kk]];
	lo.x[tdim] = coords[tdim*d_points.size() + hp[kk]];

	d_nodes[++jbox] = Teuchos::rcp( 
	    new KDTreeNode<DIM>( 
		d_nodes[tparent]->lo, hi, tparent, 0, 0, ptlo, ptlo+kk ) );

	d_nodes[++jbox] = Teuchos::rcp( 
	    new KDTreeNode<DIM>( 
		lo, d_nodes[tparent]->hi, tparent, 0, 0, ptlo+kk+1, pthi ) );

	d_nodes[tparent]->child1 = jbox-1;
	d_nodes[tparent]->child2 = jbox;

	if ( kk > 1 )
	{
	    taskparent[++nowtask] = jbox-1;
	    taskdim[nowtask] = (tdim+1) % DIM;
	}
	if ( np - kk > 3 )
	{
	    taskparent[++nowtask] = jbox;
	    taskdim[nowtask] = (tdim+1) % DIM;
	}
    }

    // Create reverse indices.
    for ( int j = 0; j < (int) d_points.size(); ++j )
    {
	d_rptindx[ d_ptindx[j] ] = j;
    }

    // Cleanup.
    free( coords );
}

/*!
 * \brief Locate the nearest neighbor point in the mesh.
 */
template<int DIM>
void KDTree<DIM>::nearestNeighbor( const std::array<double,3> &coords,
				   iBase_EntityHandle &nearest_neighbor )
{
    Point<DIM> search_point( coords[0], coords[1], coords[2] );
    int nearest_idx = nearest( search_point );
    nearest_neighbor = d_points[nearest_idx];
}

/*!
 * \brief Get the element a point is located in. Return false if we didn't
 * find anything.
 */
template<int DIM>
bool KDTree<DIM>::getElement( const std::array<double,3> &coords,
			      iBase_EntityHandle &element )
{
    element = 0;
    Point<DIM> search_point( coords[0], coords[1], coords[2] );
    return findElement( search_point, element );
}

/*! 
 * \brief Calculate the distance between a point and a box.
 */
template<int DIM>
double KDTree<DIM>::dist( const Point<DIM> &p1, const Point<DIM> &p2 )
{
    double d = 0.0;
    for ( int i = 0; i < DIM; ++i )
    {
	d += (p2.x[i] - p1.x[i])*(p2.x[i] - p1.x[i]);
    }

    double sqr_d = pow( d, 0.5 );

    return sqr_d;
}

/*!
 * \brief Calculate the distance between a point and an iMesh Point.
 */
template<int DIM>
double KDTree<DIM>::dist( iBase_EntityHandle p1, const Point<DIM> &p2 )
{
    int error = 0;

    double p1x[3] = {0.0, 0.0, 0.0};
    iMesh_getVtxCoord( d_mesh,
		       p1,
		       &p1x[0],
		       &p1x[1],
		       &p1x[2],
		       &error );
    assert( iBase_SUCCESS == error );

    double d = 0.0;

    for (int i = 0; i < DIM; ++i )
    {
	d += (p2.x[i] - p1x[i])*(p2.x[i] - p1x[i]);
    }

    double sqr_d = pow( d, 0.5 );

    return sqr_d;
}

/*! 
 * \brief Calculate the distance between a point and a box.
 */
template<int DIM>
double KDTree<DIM>::dist( const Box<DIM> &b, const Point<DIM> &p)
{
    double d = 0.0;

    for (int i = 0; i < DIM; ++i )
    {
	if ( p.x[i] < b.lo.x[i] ) d += (p.x[i]-b.lo.x[i])*(p.x[i]-b.lo.x[i]);
	if ( p.x[i] > b.hi.x[i] ) d += (p.x[i]-b.hi.x[i])*(p.x[i]-b.hi.x[i]);
    }

    double sqr_d = pow( d, 0.5 );

    return sqr_d;
}

/*!
 * \brief Partition an array by index.
 */
template<int DIM>
int KDTree<DIM>::partition( const int k, int *indx, int n, double *arr )
{
    int l = 0;
    int ir = n - 1;
    int mid = 0;
    int i = 0;
    int j = 0;
    int ia = 0;
    double a;

    for (;;)
    {
	if ( ir <= l+1 )
	{
	    if ( ir == l+1 && arr[indx[ir]] < arr[indx[l]] )
	    {
		std::swap(indx[l],indx[ir]);
	    }
	    return indx[k];
	}
	else
	{
	    mid = (l+ir) >> 1;
	    std::swap(indx[mid],indx[l+1]);

	    if ( arr[indx[l]] > arr[indx[ir]] )
	    {
		std::swap(indx[l],indx[ir]);
	    }
	    if ( arr[indx[l+1]] > arr[indx[ir]] )
	    {
		std::swap(indx[l+1],indx[ir]);
	    }
	    if ( arr[indx[l]] > arr[indx[l+1]] )
	    {
		std::swap(indx[l],indx[l+1]);
	    }

	    i = l+1;
	    j = ir;
	    ia = indx[l+1];
	    a = arr[ia];

	    for (;;)
	    {
		do i++; while (arr[indx[i]] < a);
		do j--; while (arr[indx[j]] > a);
		if (j < i) break;
		std::swap(indx[i],indx[j]);
	    }

	    indx[l+1] = indx[j];
	    indx[j] = ia;

	    if ( j >= k)
	    {
		ir = j - 1;
	    }
	    if ( j <= k )
	    {
		l = i;
	    }
	}
    }
}

/*!
 * \brief Locate the node a point is in.
 */
template<int DIM>
int KDTree<DIM>::locate( Point<DIM> p )
{
    int nb = 0;
    int jdim = 0;
    int d1 = 0;
    
    while ( d_nodes[nb]->child1 )
    {
	d1 = d_nodes[nb]->child1;
	if ( p.x[jdim] <= d_nodes[d1]->hi.x[jdim] ) nb = d1;
	else nb = d_nodes[nb]->child2;
	jdim = (jdim+1) % DIM;
    }
    return nb;
}

/*!
 * \brief Locate the node a point is in based on index.
 */
template<int DIM>
int KDTree<DIM>::locate( int ip )
{
    int jh = d_rptindx[ip];
    int nb = 0;
    int d1 = 0;
    
    while ( d_nodes[nb]->child1 )
    {
	d1 = d_nodes[nb]->child1;
	if ( jh <= d_nodes[d1]->pthi ) nb = d1;
	else nb = d_nodes[nb]->child2;
    }
    return nb;
}

/*!
 * \brief Return the index of the point that is the nearest neighbor to a
 * given point.
 */
template<int DIM>
int KDTree<DIM>::nearest( Point<DIM> p )
{
    int task[50];
    double d = 0.0;
    double dnrst = d_large;
    int nrst = 0;

    // Find the starting node.
    int k = locate(p);
    
    // Find the nearest neighbor in that starting node.
    for ( int i = d_nodes[k]->ptlo; i <= d_nodes[k]->pthi; ++i )
    {
	d = dist( d_points[d_ptindx[i]] , p );
	if ( d < dnrst )
	{
	    nrst = d_ptindx[i];
	    dnrst = d;
	}
    }

    // Traverse back up the tree searching for a closer neighbor.
    task[1] = 0;
    int ntask = 1;
    while ( ntask )
    {
	k = task[ntask--];
	if ( dist( d_nodes[k], p ) < dnrst )
	{
	    if ( d_nodes[k]->child1 )
	    {
		task[++ntask] = d_nodes[k]->child1;
		task[++ntask] = d_nodes[k]->child2;
	    }
	    else 
	    {
		for ( int i = d_nodes[k]->ptlo; i <= d_nodes[k]->pthi; ++i )
		{
		    d = dist( d_points[d_ptindx[i]], p );

		    if ( d < dnrst )
		    {
			nrst = d_ptindx[i];
			dnrst = d;
		    }
		}
	    }
	}
    }

    return nrst;
}

/*!
 * \brief Determine if a point is inside the elements adjacent to a point. 
 */
template<int DIM>
bool KDTree<DIM>::pointInAdjElements( Point<DIM> p,
				      iBase_EntityHandle point,
				      iBase_EntityHandle &element )
{
    int error = 0;
    bool return_val = false;

    std::array<double,3> coords = { p.x[0], p.x[1], p.x[2] };

    iBase_EntityHandle *adj_elements = 0;
    int adj_elements_allocated = 0;
    int adj_elements_size = 0;
    iMesh_getEntAdj( d_mesh,
		     point,
		     d_entity_type,
		     &adj_elements,
		     &adj_elements_allocated,
		     &adj_elements_size,
		     &error );
    assert( iBase_SUCCESS == error );

    for ( int i = 0; i < adj_elements_size; ++i )
    {
	if ( !return_val )
	{
	    return_val = PointQuery::pointInRefElement( d_mesh,
							adj_elements[i],
							coords );
	    if ( return_val )
	    {
		element = adj_elements[i];
	    }
	}
    }

    free( adj_elements );

    return return_val;
}

/*!
 * \brief Find the element a point resides in.
 */
template<int DIM>
bool KDTree<DIM>::findElement( Point<DIM> p, iBase_EntityHandle &element )
{
    bool point_found = false;

    int task[50];
    double d = 0.0;
    double dnrst = d_large;

    // Find the starting node.
    int k = locate(p);
    
    // Find the nearest neighbor in that starting node.
    for ( int i = d_nodes[k]->ptlo; i <= d_nodes[k]->pthi; ++i )
    {
	if ( !point_found )
	{
	    point_found = pointInAdjElements( p,
					      d_points[d_ptindx[i]], 
					      element );

	    d = dist( d_points[d_ptindx[i]] , p );

	    if ( d < dnrst )
	    {
		dnrst = d;
	    }
	}
    }

    // Traverse back up the tree searching for a closer neighbor.
    task[1] = 0;
    int ntask = 1;
    while ( ntask && !point_found )
    {
	k = task[ntask--];
	if ( dist( *d_nodes[k], p ) < dnrst )
	{
	    if ( d_nodes[k]->child1 )
	    {
		task[++ntask] = d_nodes[k]->child1;
		task[++ntask] = d_nodes[k]->child2;
	    }
	    else 
	    {
		for ( int i = d_nodes[k]->ptlo; i <= d_nodes[k]->pthi; ++i )
		{
		    if ( !point_found )
		    {
			point_found = pointInAdjElements( p,
							  d_points[d_ptindx[i]],
							  element );

			d = dist( d_points[d_ptindx[i]], p );

			if ( d < dnrst )
			{
			    dnrst = d;
			}
		    }
		}
	    }
	}
    }

    return point_found;
}

} // end namespace food

#endif // end FOOD_KDTREE_DEF_HPP

//---------------------------------------------------------------------------//
// end KDTree_Def.hpp
//---------------------------------------------------------------------------//
