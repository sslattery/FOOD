//---------------------------------------------------------------------------//
// \file FEMInterpolate.hpp
// \author Stuart Slattery
// \brief Finite element interpolation declaration.
//---------------------------------------------------------------------------//

#ifndef FOOD_FEMINTERPOLATE_DEF_HPP
#define FOOD_FEMINTERPOLATE_DEF_HPP

namespace FOOD
{

/*!
 * \brief Constructor.
 */
template<class Scalar>
FEMInterpolate<Scalar>::FEMInterpolate(RCP_TensorField dof_domain, 
				       RCP_TensorField dof_range )
    : d_dof_domain( dof_domain )
    , d_dof_range( dof_range )
    , d_kdtree(0)
{ /* ... */ }

/*!
 * \brief Destructor.
 */
template<class Scalar>
FEMInterpolate<Scalar>::~FEMInterpolate()
{ /* ... */ }

/*!
 * \brief Setup for interpolation.
 */
template<class Scalar>
void FEMInterpolate<Scalar>::setup()
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
	new KDTree<3>( d_dof_domain->getDomain(),
		       d_dof_domain->getDFuncKernel()->getEvalType(),
		       d_dof_domain->getDFuncKernel()->getEvalTopology() ) );
    d_kdtree->buildTree();

    // Generate a mapping for interpolation.
    MDArray local_coords(1,3);
    iBase_EntityHandle found_entity = 0;
    int num_found = 0;
    for ( int n = 0; n < range_vertices_size; ++n )
    {
	found_entity = 0;

	local_coords(0,0) = coord_array[3*n];
	local_coords(0,1) = coord_array[3*n+1];
	local_coords(0,2) = coord_array[3*n+2];

	if ( d_kdtree->getElement( local_coords, found_entity ) )
	{
	    ++num_found;
	    d_map.insert( std::pair<iBase_EntityHandle,iBase_EntityHandle>( 
			      range_vertices[n], found_entity ) );
	}
    }
    
    std::cout << " PERCENT FOUND " << (double) num_found / range_vertices_size << std::endl;
    // Cleanup.
    free( range_vertices );
    free( coord_array );
}

/*!
 * \brief Perform value interpolation of the degrees of freedom from the
 * domain to the range.
 */
template<class Scalar>
void FEMInterpolate<Scalar>::interpolateValueDF()
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
    MDArray local_vals( 1, 1, num_components );
    MDArray local_coords( 1, 3 );
    for ( int p = 0; p < range_vertices_size; ++p )
    {
	local_coords(0,0) = coord_array[3*p];
	local_coords(0,1) = coord_array[3*p+1];
	local_coords(0,2) = coord_array[3*p+2];

	if ( d_map[ range_vertices[p] ] )
	{
	    d_dof_domain->evaluateDF( d_map[ range_vertices[p] ],
				      local_coords,
				      false,
				      local_vals );

	    for ( int n = 0; n < num_components; ++n )
	    {
		interpolated_vals[p*num_components + n] = local_vals(0,0,n);
	    }
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

#endif // end FOOD_FEMINTERPOLATE_DEF_HPP

