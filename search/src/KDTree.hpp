//---------------------------------------------------------------------------//
/*! 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file KDTree.hpp
 * \author Stuart Slattery
 * \brief KDTree declaration.
 */
//---------------------------------------------------------------------------//

#ifndef FOOD_KDTREE_HPP
#define FOOD_KDTREE_HPP

#include <vector>

#include "TopologyTools.hpp"

#include <iBase.h>
#include <iMesh.h>

#include <Teuchos_RCP.hpp>

namespace FOOD
{

//---------------------------------------------------------------------------//
/*! 
 * \brief Point structure for the KDTree.
 */
//---------------------------------------------------------------------------//
template<int DIM>
class Point
{
  public:

    // Coordinates.
    double x[DIM];

    // Copy constructor.
    Point( const Point &p )
    {
	for ( int i = 0; i < DIM; ++i) x[i] = p.x[i];
    }
    
    // Assignment operator.
    Point& operator= ( const Point &p )
    {
	for ( int i = 0; i < DIM; ++i ) x[i] = p.x[i];
	return *this;
    }

    // Comparison operator.
    bool operator== ( const Point &p ) const
    {
	for ( int i = 0; i < DIM; ++i )
	{
	    if ( x[i] != p.x[i] ) return false;
	}
	return true;
    }

    // Coordinate Constructor.
    Point( double x0 = 0.0, double x1 = 0.0, double x2 = 0.0 )
    {
	x[0] = x0;
	if (DIM > 1) x[1] = x1;
	if (DIM > 2) x[2] = x2;
    }
};

//---------------------------------------------------------------------------//
/*!
 * \brief Box structure for the KDTree.
 */
//---------------------------------------------------------------------------//
template<int DIM>
class Box
{
  public:

    Point<DIM> lo;
    Point<DIM> hi;
    
    Box()
    { /* ... */ }

    Box( const Point<DIM> &_lo, const Point<DIM> &_hi)
	: lo(_lo)
	, hi(_hi)
    { /* ... */ }

    ~Box()
    { /* ... */ }
};

//---------------------------------------------------------------------------//
/*!
 * \brief Node structure for the KDTree.
 */
//---------------------------------------------------------------------------//
template<int DIM>    
class KDTreeNode : public Box<DIM>
{
  public:
    int parent;
    int child1;
    int child2;
    int ptlo;
    int pthi;
 
    KDTreeNode()
	: parent(0)
	, child1(0)
	, child2(0)
	, ptlo(0)
	, pthi(0)
    { /* ... */ }

    KDTreeNode( Point<DIM> _lo,
		Point<DIM> _hi,
		int _parent,
		int _child1,
		int _child2,
		int _ptlo,
		int _pthi )
	: Box<DIM>( _lo, _hi )
	, parent(_parent)
	, child1(_child1)
	, child2(_child2)
	, ptlo(_ptlo)
	, pthi(_pthi)
    { /* ... */ }

    ~KDTreeNode()
    { /* ... */ }
};

//---------------------------------------------------------------------------//
/*!
 * \brief KDTree data structure.
 *
 * This class constructs a KDTree of the vertices in the provided mesh
 * set. Nearest neighbor searches are provided for all types and topologies
 * equal to and below the those provided to the constructor. Offered for 1, 2,
 * and 3 dimensional meshes.
 */
//---------------------------------------------------------------------------//
template<int DIM>
class KDTree
{
  public:

    //@{
    //! Typedefs.
    typedef Teuchos::RCP< KDTreeNode<DIM> >           RCP_Node;
    //@}

  private:

    // The mesh instance this tree is generated for.
    iMesh_Instance d_mesh;

    // The mesh set this tree is generated for.
    iBase_EntitySetHandle d_mesh_set;

    // The entity type this tree can search for.
    std::size_t d_entity_type;

    // The entity topology this tree can search for.
    std::size_t d_entity_topology;

    // Large value to get the tree started.
    double d_large;

    // Points in the domain corresponding to linear entities.
    std::vector<iBase_EntityHandle> d_points;

    // Point indices.
    std::vector<int> d_ptindx;

    // Reverse point indices.
    std::vector<int> d_rptindx;

    // Tree nodes.
    std::vector<RCP_Node> d_nodes;

  public:

    // Constructor.
    KDTree( iMesh_Instance mesh, 
	    iBase_EntitySetHandle mesh_set,
	    const int entity_type,
	    const int entity_topology );

    // Destructor.
    ~KDTree();

    // Build the tree.
    void buildTree();

    // Locate the nearest neighbor point in the mesh set (linear entity points
    // only).
    void nearestNeighbor( const double coords[3],
			  iBase_EntityHandle &nearest_neighbor );

    // Get the element a point is located in. Return false if the point is not
    // in the domain of the mesh set.
    bool getElement( const double coords[3],
		     iBase_EntityHandle &element );

  private:

    // Calculate the distance between two points.
    double dist( const Point<DIM> &p1, const Point<DIM> &p2 );

    // Calculate the distance between a point and an iMesh Point.
    double dist( iBase_EntityHandle p1, const Point<DIM> &p2 );

    // Calculate the distance between a point and a box.
    double dist( const Box<DIM> &b, const Point<DIM> &p );

    // Partition an array by index.
    int partition( const int k, int *indx, int n, double *arr );

    // Locate the node a point is in.
    int locate( Point<DIM> p );

    // Locate the node a point is in based on index.
    int locate( int ip );

    // Return the index of the point that is nearest neighbor to a given
    // point. 
    int nearest( Point<DIM> p );

    // Determine if a point is inside the elements adjacent to a point.
    bool pointInAdjElements( Point<DIM> p,
			     iBase_EntityHandle point,
			     iBase_EntityHandle &element );

    // Find the element a point resides in.
    bool findElement( Point<DIM> p, iBase_EntityHandle &element );
};

} // end namespace FOOD

#include "KDTree_Def.hpp"

#endif // end FOOD_KDTREE_HPP

//---------------------------------------------------------------------------//
// end KDTree.hpp
//---------------------------------------------------------------------------//
