//---------------------------------------------------------------------------//
// \file cxx_main.cpp
// \author Stuart Slattery
// \brief Driver for tet part interpolation example.
//---------------------------------------------------------------------------//

#include <cassert>
#include <string>
#include <iostream>
#include <map>

#include <Domain.hpp>
#include <TensorTemplate.hpp>
#include <DFuncKernel.hpp>
#include <TensorField.hpp>
#include <PointQuery.hpp>

#include <iBase.h>
#include <iMesh.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ENull.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include <Intrepid_FieldContainer.hpp>

int main(int argc, char* argv[])
{
    typedef Intrepid::FieldContainer<double> MDArray;

    int error;

    // The tensor template and dfunckernel can be shared by both fields.
    Teuchos::RCP<FOOD::TensorTemplate> tensor_template = Teuchos::rcp(
	new FOOD::TensorTemplate(0, 1, FOOD::FOOD_REAL, Teuchos::null) );

    Teuchos::RCP< FOOD::DFuncKernel<double> > dfunckernel =
	Teuchos::rcp( new FOOD::DFuncKernel<double>( iMesh_TETRAHEDRON,
						     iBase_VERTEX,
						     iMesh_POINT,
						     FOOD::FOOD_CARTESIAN,
						     FOOD::FOOD_FEM,
						     FOOD::FOOD_HGRAD,
						     1 ) );

    // Set up the fine mesh.
    iMesh_Instance fine_mesh;
    iMesh_newMesh("", &fine_mesh, &error, 0);
    assert( iBase_SUCCESS == error );

    iBase_EntitySetHandle fine_root_set;
    iMesh_getRootSet( fine_mesh, &fine_root_set, &error );
    assert( iBase_SUCCESS == error );

    std::string fine_mesh_filename = "tagged_fine_tet_part.vtk";
    iMesh_load( fine_mesh, 
		fine_root_set, 
		&fine_mesh_filename[0], 
		"", 
		&error,
		(int) fine_mesh_filename.size(),
		0 );

    Teuchos::RCP<FOOD::Domain> fine_domain = Teuchos::rcp(
	new FOOD::Domain(fine_mesh, fine_root_set) );

    FOOD::TensorField<double> fine_field( getDefaultComm<int>(),
					  fine_domain,
					  dfunckernel,
					  FOOD::FOOD_CARTESIAN, 
					  tensor_template,
					  Teuchos::null,
					  "FINE_FIELD" );

    std::string fine_tag_name = "domain";
    iBase_TagHandle fine_tag;
    iMesh_getTagHandle( fine_domain->getMesh(),
			&fine_tag_name[0],
			&fine_tag,
			&error,
			(int) fine_tag_name.size() );
    TEST_ASSERT( iBase_SUCCESS == error );

    fine_field.attachToTagData( fine_tag, error );
    TEST_ASSERT( iBase_SUCCESS == error );

    // Set up the coarse mesh.
    iMesh_Instance coarse_mesh;
    iMesh_newMesh("", &coarse_mesh, &error, 0);
    assert( iBase_SUCCESS == error );

    iBase_EntitySetHandle coarse_root_set;
    iMesh_getRootSet( coarse_mesh, &coarse_root_set, &error );
    assert( iBase_SUCCESS == error );

    std::string coarse_mesh_filename = "tagged_coarse_tet_part.vtk";
    iMesh_load( coarse_mesh, 
		coarse_root_set, 
		&coarse_mesh_filename[0], 
		"", 
		&error,
		(int) coarse_mesh_filename.size(),
		0 );

    Teuchos::RCP<FOOD::Domain> coarse_domain = Teuchos::rcp(
	new FOOD::Domain(coarse_mesh, coarse_root_set) );

    FOOD::TensorField<double> coarse_field( getDefaultComm<int>(),
					    coarse_domain,
					    dfunckernel,
					    FOOD::FOOD_CARTESIAN, 
					    tensor_template,
					    Teuchos::null,
					    "COARSE_FIELD" );

    std::string coarse_tag_name = "range";
    iBase_TagHandle coarse_tag;
    iMesh_getTagHandle( coarse_domain->getMesh(),
			&coarse_tag_name[0],
			&coarse_tag,
			&error,
			(int) coarse_tag_name.size() );
    TEST_ASSERT( iBase_SUCCESS == error );

    coarse_field.attachToTagData( coarse_tag, error );
    TEST_ASSERT( iBase_SUCCESS == error );

    // Generate the mapping for the interpolation.
    iBase_EntityHandle *domain_elements = 0;
    int domain_elements_allocated = 0;
    int domain_elements_size = 0;
    iMesh_getEntities( fine_domain->getMesh(),
		       fine_domain->getMeshSet(),
		       iBase_REGION,
		       iMesh_TETRAHEDRON,
		       &domain_elements,
		       &domain_elements_allocated,
		       &domain_elements_size,
		       &error );
    TEST_ASSERT( iBase_SUCCESS == error );

    iBase_EntityHandle *range_vertices = 0;
    int range_vertices_allocated = 0;
    int range_vertices_size = 0;
    iMesh_getEntities( coarse_domain->getMesh(),
		       coarse_domain->getMeshSet(),
		       iBase_VERTEX,
		       iMesh_POINT,
		       &range_vertices,
		       &range_vertices_allocated,
		       &range_vertices_size,
		       &error );
    TEST_ASSERT( iBase_SUCCESS == error );

    int coords_allocated = range_vertices_size*3;
    int coords_size = 0;
    double *coord_array = 0;
    iMesh_getVtxArrCoords( coarse_domain->getMesh(),
			   range_vertices,
			   range_vertices_size,
			   iBase_INTERLEAVED,
			   &coord_array,
			   &coords_allocated,
			   &coords_size,
			   &error );
    assert( iBase_SUCCESS == error );

    std::map<iBase_EntityHandle,iBase_EntityHandle> range_to_domain;
    bool point_mapped = false;
    MDArray local_coords(1,3);
    int m = 0;
    for ( int n = 0; n < range_vertices_size; ++n )
    {
	point_mapped = false;
	m = 0;
	while ( m < domain_elements_size && !point_mapped )
	{
	    local_coords[0][0] = coord_array[3*n];
	    local_coords[0][1] = coord_array[3*n+1];
	    local_coords[0][2] = coord_array[3*n+2];

	    if ( point_in_ref_element( fine_domain->getMesh(),
				       domain_elements[m],
				       local_coords ) )
	    {
		range_to_domain.insert(
		    std::pair<iBase_EntityHandle,iBase_EntityHandle>( 
			range_vertices[n], domain_elements[m] ) );
		point_mapped = true;
	    }
	    ++m;
	}
    }
    
    // perform the interpolation
    Teuchos::ArrayRCP<double> interpolated_vals(range_vertices_size);
    MDArray local_vals(1,1);
    for ( int p = 0; p < range_vertices_size; ++p )
    {
	local_coords[0][0] = coord_array[3*p];
	local_coords[0][1] = coord_array[3*p+1];
	local_coords[0][2] = coord_array[3*p+2];

	if ( range_to_domain( range_vertices[p] ) )
	{
	    fine_field.evaluateDF( range_to_domain( range_vertices[p] ),
				   local_coords,
				   false,
				   local_vals );
	    interpolated_vals[p] = local_vals(0,0);
	}
	else
	{
	    interpolated_vals[p] = 0.0;
	}
    }

    coarse_field.attachToArrayData( interpolated_vals,
				    iBase_INTERLEAVED,
				    error );
    assert( iBase_SUCCESS == error );

    std::string interp_file = "interpolated_coarse_part.vtk";
    iMesh_save( coarse_domain->getMesh(),
		coarse_domain->getMeshSet(),
		&interp_file[0],
		"",
		&error,
		(int) interp_file.size(),
		0 );
    assert( iBase_SUCCESS == error );

    // cleanup
    free( domain_elements );
    free( range_vertices );
    free( coord_array );

    return 0;
}
