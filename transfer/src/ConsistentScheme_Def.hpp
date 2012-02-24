//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file ConsistentScheme.hpp
 * \author Stuart Slattery
 * \brief Consistent finite element interpolation scheme definition.
 */
//---------------------------------------------------------------------------//

#ifndef FOOD_CONSISTENTSCHEME_DEF_HPP
#define FOOD_CONSISTENTSCHEME_DEF_HPP

#include <Teuchos_ArrayRCP.hpp>

namespace FOOD
{

/*!
 * \brief Constructor.
 */
template<class Scalar>
ConsistentScheme<Scalar>::ConsistentScheme(RCP_TensorField dof_domain, 
					   RCP_TensorField dof_range )
    : d_dof_domain( dof_domain )
    , d_dof_range( dof_range )
    , d_kdtree(0)
{ /* ... */ }

/*!
 * \brief Destructor.
 */
template<class Scalar>
ConsistentScheme<Scalar>::~ConsistentScheme()
{ /* ... */ }

/*!
 * \brief Setup for interpolation.
 */
template<class Scalar>
void ConsistentScheme<Scalar>::setup()
{
    int error = 0;

    // Get the range mesh vertices for interpolation.
    iBase_EntityHandle *range_vertices = 0;
    int range_vertices_allocated = 0;
    int range_vertices_size = 0;
    iMesh_getEntities( d_dof_range->getDomain()->getMesh(),
		       d_dof_range->getDomain()->getMeshSet(),
		       iBase_VERTEX,
		       iMesh_POINT,
		       &range_vertices,
		       &range_vertices_allocated,
		       &range_vertices_size,
		       &error );
    assert( iBase_SUCCESS == error );

    int coords_allocated = range_vertices_size*3;
    int coords_size = 0;
    double *coord_array = 0;
    iMesh_getVtxArrCoords( d_dof_range->getDomain()->getMesh(),
			   range_vertices,
			   range_vertices_size,
			   iBase_INTERLEAVED,
			   &coord_array,
			   &coords_allocated,
			   &coords_size,
			   &error );
    assert( iBase_SUCCESS == error );

    // Setup the kdtree to search the domain mesh with range vertices.
    d_kdtree = Teuchos::rcp( 
	new KDTree<3>( d_dof_domain->getDomain()->getMesh(),
		       d_dof_domain->getDomain()->getMeshSet(),
		       d_dof_domain->getDFuncKernel()->getEntityType(),
		       d_dof_domain->getDFuncKernel()->getEntityTopology() ) );
    d_kdtree->buildTree();

    // Generate a mapping for interpolation.
    double local_coords[3] = { 0.0, 0.0, 0.0 };
    iBase_EntityHandle found_entity = 0;
    for ( int n = 0; n < range_vertices_size; ++n )
    {
	found_entity = 0;

	local_coords[0] = coord_array[3*n];
	local_coords[1] = coord_array[3*n+1];
	local_coords[2] = coord_array[3*n+2];

	if ( d_kdtree->getElement( local_coords, found_entity ) )
	{
	    d_map.insert( std::pair<iBase_EntityHandle,iBase_EntityHandle>( 
			      range_vertices[n], found_entity ) );
	}
    }
    
    // Cleanup.
    free( range_vertices );
    free( coord_array );
}

/*!
 * \brief Perform value interpolation of the degrees of freedom from the
 * domain to the range.
 */
template<class Scalar>
void ConsistentScheme<Scalar>::transferValueDF()
{
    int error = 0;

    // Get the range mesh vertices for interpolation.
    iBase_EntityHandle *range_vertices = 0;
    int range_vertices_allocated = 0;
    int range_vertices_size = 0;
    iMesh_getEntities( d_dof_range->getDomain()->getMesh(),
		       d_dof_range->getDomain()->getMeshSet(),
		       iBase_VERTEX,
		       iMesh_POINT,
		       &range_vertices,
		       &range_vertices_allocated,
		       &range_vertices_size,
		       &error );
    assert( iBase_SUCCESS == error );

    int coords_allocated = range_vertices_size*3;
    int coords_size = 0;
    double *coord_array = 0;
    iMesh_getVtxArrCoords( d_dof_range->getDomain()->getMesh(),
			   range_vertices,
			   range_vertices_size,
			   iBase_INTERLEAVED,
			   &coord_array,
			   &coords_allocated,
			   &coords_size,
			   &error );
    assert( iBase_SUCCESS == error );

    // Do the interpolation.
    int num_components = d_dof_domain->getTensorTemplate()->getNumComponents();
    int dof_size = range_vertices_size*num_components;
    Teuchos::ArrayRCP<double> interpolated_vals(dof_size, 0.0);
    Teuchos::ArrayRCP<double> ent_vals;
    for ( int p = 0; p < range_vertices_size; ++p )
    {
	double local_coords[3] = { coord_array[3*p],
				   coord_array[3*p+1],
				   coord_array[3*p+2] };

	if ( d_map[ range_vertices[p] ] )
	{
	    ent_vals = interpolated_vals.persistingView(p, num_components);
	    d_dof_domain->evaluateDF( d_map[ range_vertices[p] ],
				      local_coords,
				      false,
				      ent_vals );
	}
    }
    d_dof_range->attachToArrayData( interpolated_vals,
				    iBase_INTERLEAVED,
				    error );
    assert( iBase_SUCCESS == error );

    // Cleanup.
    free( range_vertices );
    free( coord_array );
}

} // end namespace FOOD

#endif // end FOOD_CONSISTENTSCHEME_DEF_HPP

