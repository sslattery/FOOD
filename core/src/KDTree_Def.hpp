//---------------------------------------------------------------------------//
// \file KDTree_Def.hpp
// \author Stuart Slattery
// \brief KDTree definition.
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
KDTree<DIM>::KDTree( RCP_Domain domain, 
		     const int entity_type,
		     const int entity_topology )
    : d_domain(domain)
    , d_entity_type(entity_type)
    , d_entity_topology(entity_topology)
    , d_large(1.0e99)
    , d_num_points(0)
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
    EntityHandle *elements = 0;
    int elements_allocated = 0;
    int elements_size = 0;
    iMesh_getEntities( d_domain->getMesh(),
		       d_domain->getMeshSet(),
		       d_entity_type,
		       d_entity_topology,
		       &elements,
		       &elements_allocated,
		       &elements_size,
		       &error );
    assert( iBase_SUCCESS == error );

    for ( int n = 0; n < elements_size; ++n )
    {
	EntityHandle *element_nodes = 0;
	int element_nodes_allocated = 0;
	int element_nodes_size = 0;
	iMesh_getEntAdj( d_domain->getMesh(),
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
	    ++d_num_points;
	}

	free( element_nodes );
    }

    free( elements );

    int coords_allocated = 0;
    int coords_size = 0;
    double *coords;
    iMesh_getVtxArrCoords( d_domain->getMesh(),
			   &d_points[0],
			   d_num_points,
			   iBase_BLOCKED,
			   &coords,
			   &coords_allocated,
			   &coords_size,
			   &error );
    assert( iBase_SUCCESS == error );
    assert( d_num_points * 3 == coords_size );

    // Setup the point indices.
    d_ptindx.resize( d_num_points );
    d_rptindx.resize( d_num_points );
    for ( int k = 0; k < d_num_points; ++k )
    {
	d_ptindx[k] = k;
    }

    // Compute the number of tree nodes.
    int m = 1;
    for ( int n = d_num_points; n; n >>=1 )
    {
	m <<= 1;
    }
    int num_nodes = 2*d_num_points - (m >> 1);
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
	new KDTreeNode<DIM>( lo, hi, 0, 0, 0, 0, d_num_points-1 ) );

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
	cp = &coords[tdim*d_num_points];
	np = pthi - ptlo + 1;
	kk = (np-1) / 2;

	(void) partition( kk, hp, np, cp );

	// Create the children.
	hi = d_nodes[tparent]->hi;
	lo = d_nodes[tparent]->lo;
	hi.x[tdim] = coords[tdim*d_num_points + hp[kk]];
	lo.x[tdim] = coords[tdim*d_num_points + hp[kk]];

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
    for ( int j = 0; j < d_num_points; ++j )
    {
	d_rptindx[ d_ptindx[j] ] = j;
    }

    // Cleanup.
    free( coords );
}

/*!
 * \brief Locate the nearest neighbor in the mesh.
 */
template<int DIM>
void KDTree<DIM>::nearestNeighbor( const MDArray &coords,
				   EntityHandle &nearest_neighbor )
{
    Point<DIM> search_point( coords(0,0), coords(0,1), coords(0,2) );
    int nearest_idx = nearest( search_point );
    nearest_neighbor = d_points[nearest_idx];
}

/*!
 * \brief Get the element a point is located in. Return false if we didn't
 * find anything.
 */
template<int DIM>
bool KDTree<DIM>::getElement( const MDArray &coords,
			      EntityHandle &element )
{
    element = 0;
    Point<DIM> search_point( coords(0,0), coords(0,1), coords(0,2) );
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
double KDTree<DIM>::dist( EntityHandle p1, const Point<DIM> &p2 )
{
    int error = 0;

    double p1x[3] = {0.0, 0.0, 0.0};
    iMesh_getVtxCoord( d_domain->getMesh(),
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
		task[++ntask] = d_nodes[k]->chidl2;
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
				      EntityHandle point,
				      EntityHandle &element )
{
    int error = 0;
    bool return_val = false;

    MDArray coords(1,3);
    coords(0,0) = p.x[0];
    coords(0,1) = p.x[1];
    coords(0,2) = p.x[2];

    EntityHandle *adj_elements = 0;
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

    for ( int i = 0; i < adj_elements_size; ++i )
    {
	if ( !return_val )
	{
	    return_val = PointQuery::pointInRefElement( d_domain->getMesh(),
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
bool KDTree<DIM>::findElement( Point<DIM> p, EntityHandle &element )
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
