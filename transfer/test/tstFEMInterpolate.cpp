//----------------------------------*-C++-*----------------------------------//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
//
/*!
 * \file   mesh/test/tstFEMInterpolate.cpp
 * \author Stuart Slattery
 * \brief  FEMInterpolate class unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <cassert>

#include <Types.hpp>
#include <Quantity.hpp>
#include <Unit.hpp>
#include <Domain.hpp>
#include <DFuncKernel.hpp>
#include <TensorTemplate.hpp>
#include <TensorField.hpp>
#include <FEMInterpolate.hpp>

#include <iMesh.h>
#include <iBase.h>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include <Intrepid_FieldContainer.hpp>

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//

template<class Ordinal>
Teuchos::RCP<const Teuchos::Comm<Ordinal> > getDefaultComm()
{
#ifdef HAVE_MPI
    return Teuchos::DefaultComm<Ordinal>::getComm();
#else
    return Teuchos::rcp(new Teuchos::SerialComm<Ordinal>() );
#endif
}

void create_hex_mesh(iMesh_Instance &mesh)
{
    // Create the mesh instance.
    int error;
    iMesh_newMesh("", &mesh, &error, 0);
    assert( iBase_SUCCESS == error );

    // Create a tag for the vertex elements.
    iBase_TagHandle vertex_tag;
    std::string vertex_tag_name = "vertex_tag";
    iMesh_createTag( mesh,
		     &vertex_tag_name[0],
		     1,
		     iBase_DOUBLE,
		     &vertex_tag,
		     &error,
		     (int) vertex_tag_name.size());
    assert( iBase_SUCCESS == error );

    // Generate vertices.
    int num_i = 11;
    int num_j = 11;
    int num_k = 11;
    int dx = 1.0;
    int dy = 1.0;
    int dz = 1.0;
    int idx = 0;
    std::vector<double> coords(3*num_i*num_j*num_k, 0.0);
    for ( int k = 0; k < num_k; ++k )
    {
	for ( int j = 0; j < num_j; ++j )
	{
	    for ( int i = 0; i < num_i; ++i )
	    {
		idx = i + num_i*j + num_i*num_j*k;
		coords[3*idx]     = i*dx;
		coords[3*idx + 1] = j*dy;
		coords[3*idx + 2] = k*dz;
	    }
	}
    }

    iBase_EntityHandle *vertices;
    int vertices_allocated = 0;
    int vertices_size = 0;
    iMesh_createVtxArr( mesh,
			(int) coords.size() / 3,
			iBase_INTERLEAVED,
			&coords[0],
			(int) coords.size(),
			&vertices,
			&vertices_allocated,
			&vertices_size,
			&error);
    assert( iBase_SUCCESS == error );
    assert( vertices_allocated == num_i*num_j*num_k );
    assert( vertices_size == num_i*num_j*num_k );

    // Tag the vertices.
    std::vector<double> vertex_tag_data( (int) coords.size() / 3, 0.0 );
    double data = 0.0;
    std::vector<double>::iterator vertex_data_iterator;
    for (vertex_data_iterator = vertex_tag_data.begin();
	 vertex_data_iterator != vertex_tag_data.end();
	 ++vertex_data_iterator )
    {
	*vertex_data_iterator = data;
	data += 1.0;
    }

    iMesh_setDblArrData( mesh,
			 vertices,
			 vertices_size,
			 vertex_tag,
			 &vertex_tag_data[0],
			 (int) vertex_tag_data.size(),
			 &error);
    assert( iBase_SUCCESS == error );

    // Create a vertex set and add the vertex array.
    iBase_EntitySetHandle vertex_set;
    iMesh_createEntSet(mesh, 1, &vertex_set, &error);
    assert( iBase_SUCCESS == error );

    iMesh_addEntArrToSet(mesh, vertices, vertices_size, vertex_set, &error);
    assert( iBase_SUCCESS == error );

    // Create a hex set.
    iBase_EntitySetHandle hex_set;
    iMesh_createEntSet(mesh, 1, &hex_set, &error);
    assert( iBase_SUCCESS == error );

    // Create a tag for the hex elements.
    iBase_TagHandle hex_tag;
    std::string hex_tag_name = "hex_tag";
    iMesh_createTag( mesh,
		     &hex_tag_name[0],
		     1,
		     iBase_DOUBLE,
		     &hex_tag,
		     &error,
		     (int) hex_tag_name.size());
    assert( iBase_SUCCESS == error );

    // Generate hexahedrons from vertices, tag them, and add them to the hex
    // set. 
    iBase_EntityHandle connectivity[8];
    int hex_creation_status = 0;
    double hex_data = 0.0;
    for (int k = 0; k < num_k - 1; ++k)
    {
	for (int j = 0; j < num_j - 1; ++j)
	{
	    for (int i = 0; i < num_i - 1; ++i)
	    {
		idx = i + num_i*j + num_i*num_j*k;
		connectivity[0] = vertices[ idx ];
		connectivity[1] = vertices[ idx + 1 ];
		connectivity[2] = vertices[ idx + 1 + num_i ];
		connectivity[3] = vertices[ idx +     num_i ];
		connectivity[4] = vertices[ idx +             num_i*num_j ];
		connectivity[5] = vertices[ idx + 1 +         num_i*num_j ];
		connectivity[6] = vertices[ idx + 1 + num_i + num_i*num_j ];
		connectivity[7] = vertices[ idx +     num_i + num_i*num_j ];

		iBase_EntityHandle hexahedron;
		iMesh_createEnt( mesh,
				 iMesh_HEXAHEDRON,
				 connectivity,
				 8,
				 &hexahedron,
				 &hex_creation_status,
				 &error);
		assert( iBase_SUCCESS == error );
		assert( iBase_NEW == hex_creation_status );

		iMesh_setDblData( mesh,
				  hexahedron,
				  hex_tag,
				  hex_data,
				  &error);
		assert( iBase_SUCCESS == error );
		
		hex_data += 1.0;

		iMesh_addEntToSet(mesh,
				  hexahedron,
				  hex_set,
				  &error);
		assert( iBase_SUCCESS == error );
	    }
	}
    }
}

void create_tet_mesh( iMesh_Instance &mesh )
{
    // Setup iMesh instance
    int error;
    iMesh_newMesh("", &mesh, &error, 0);
    assert( iBase_SUCCESS == error );

    // Create a tag for the vertex elements.
    iBase_TagHandle vertex_tag;
    std::string vertex_tag_name = "vertex_tag";
    iMesh_createTag( mesh,
		     &vertex_tag_name[0],
		     1,
		     iBase_DOUBLE,
		     &vertex_tag,
		     &error,
		     (int) vertex_tag_name.size());
    assert( iBase_SUCCESS == error );

    // Generate vertices.
    int num_i = 11;
    int num_j = 11;
    int num_k = 11;
    int dx = 1.0;
    int dy = 1.0;
    int dz = 1.0;
    int idx = 0;
    std::vector<double> coords(3*num_i*num_j*num_k, 0.0);
    for ( int k = 0; k < num_k; ++k )
    {
	for ( int j = 0; j < num_j; ++j )
	{
	    for ( int i = 0; i < num_i; ++i )
	    {
		idx = i + num_i*j + num_i*num_j*k;
		coords[3*idx]     = i*dx;
		coords[3*idx + 1] = j*dy;
		coords[3*idx + 2] = k*dz;
	    }
	}
    }

    iBase_EntityHandle *vertices;
    int vertices_allocated = 0;
    int vertices_size = 0;
    iMesh_createVtxArr( mesh,
			(int) coords.size() / 3,
			iBase_INTERLEAVED,
			&coords[0],
			(int) coords.size(),
			&vertices,
			&vertices_allocated,
			&vertices_size,
			&error);
    assert( iBase_SUCCESS == error );
    assert( vertices_allocated == num_i*num_j*num_k );
    assert( vertices_size == num_i*num_j*num_k );

    // Tag the vertices.
    std::vector<double> vertex_tag_data( (int) coords.size() / 3, 0.0 );
    double data = 0.0;
    std::vector<double>::iterator vertex_data_iterator;
    for (vertex_data_iterator = vertex_tag_data.begin();
	 vertex_data_iterator != vertex_tag_data.end();
	 ++vertex_data_iterator )
    {
	*vertex_data_iterator = data;
	data += 1.0;
    }

    iMesh_setDblArrData( mesh,
			 vertices,
			 vertices_size,
			 vertex_tag,
			 &vertex_tag_data[0],
			 (int) vertex_tag_data.size(),
			 &error);
    assert( iBase_SUCCESS == error );

    // Create a vertex set and add the vertex array.
    iBase_EntitySetHandle vertex_set;
    iMesh_createEntSet(mesh, 1, &vertex_set, &error);
    assert( iBase_SUCCESS == error );

    iMesh_addEntArrToSet(mesh, vertices, vertices_size, vertex_set, &error);
    assert( iBase_SUCCESS == error );

    // Create a tet set.
    iBase_EntitySetHandle tet_set;
    iMesh_createEntSet(mesh, 1, &tet_set, &error);
    assert( iBase_SUCCESS == error );

    // Create a tag for the tet elements.
    iBase_TagHandle tet_tag;
    std::string tet_tag_name = "tet_tag";
    iMesh_createTag(mesh,
		    &tet_tag_name[0],
		    3,
		    iBase_DOUBLE,
		    &tet_tag,
		    &error,
		    (int) tet_tag_name.size());
    assert( iBase_SUCCESS == error );

    // Generate tetrahedrons from vertices, tag them, and add them to the tet
    // set. Decompose each cube into 5 tetrahedrons.
    iBase_EntityHandle connectivity[4];
    int tet_creation_status = 0;
    double tet_data = 0.0;
    double data_arr[3] = {tet_data, tet_data, tet_data};
    int v0 = 0;
    int v1 = 0;
    int v2 = 0;
    int v3 = 0;
    int v4 = 0;
    int v5 = 0;
    int v6 = 0;
    int v7 = 0;
    for (int k = 0; k < num_k - 1; ++k)
    {
	for (int j = 0; j < num_j - 1; ++j)
	{
	    for (int i = 0; i < num_i - 1; ++i)
	    {
		// Starting corner vertex index.
		idx = i + num_i*j + num_i*num_j*k;
		v0 = idx;
		v1 = idx + 1;
		v2 = idx + 1 + num_i;
		v3 = idx +     num_i;
		v4 = idx +             num_i*num_j;
		v5 = idx + 1 +         num_i*num_j;
		v6 = idx + 1 + num_i + num_i*num_j;
		v7 = idx +     num_i + num_i*num_j;

		// Tetrahedron 1.
		connectivity[0] = vertices[ v0 ];
		connectivity[1] = vertices[ v1 ];
		connectivity[2] = vertices[ v3 ];
		connectivity[3] = vertices[ v4 ];

		iBase_EntityHandle tetrahedron_1;
		iMesh_createEnt( mesh,
				 iMesh_TETRAHEDRON,
				 connectivity,
				 4,
				 &tetrahedron_1,
				 &tet_creation_status,
				 &error);
		assert( iBase_SUCCESS == error );
		assert( iBase_NEW == tet_creation_status );

		data_arr[0] = tet_data;
		data_arr[1] = tet_data;
		data_arr[2] = tet_data;
		iMesh_setDblArrData( mesh,
				     &tetrahedron_1,
				     1,
				     tet_tag,
				     data_arr,
				     3,
				     &error);
		assert( iBase_SUCCESS == error );
		
		tet_data += 1.0;

		iMesh_addEntToSet( mesh,
				   tetrahedron_1,
				   tet_set,
				   &error);
		assert( iBase_SUCCESS == error );

		// Tetrahedron 2.
		connectivity[0] = vertices[ v1 ];
		connectivity[1] = vertices[ v2 ];
		connectivity[2] = vertices[ v3 ];
		connectivity[3] = vertices[ v6 ];

		iBase_EntityHandle tetrahedron_2;
		iMesh_createEnt(mesh,
				iMesh_TETRAHEDRON,
				connectivity,
				4,
				&tetrahedron_2,
				&tet_creation_status,
				&error);
		assert( iBase_SUCCESS == error );
		assert( iBase_NEW == tet_creation_status );

		data_arr[0] = tet_data;
		data_arr[1] = tet_data;
		data_arr[2] = tet_data;
		iMesh_setDblArrData( mesh,
				     &tetrahedron_2,
				     1,
				     tet_tag,
				     data_arr,
				     3,
				     &error);
		assert( iBase_SUCCESS == error );
		
		tet_data += 1.0;

		iMesh_addEntToSet( mesh,
				   tetrahedron_2,
				   tet_set,
				   &error);
		assert( iBase_SUCCESS == error );

		// Tetrahedron 3.
		connectivity[0] = vertices[ v6 ];
		connectivity[1] = vertices[ v5 ];
		connectivity[2] = vertices[ v4 ];
		connectivity[3] = vertices[ v1 ];

		iBase_EntityHandle tetrahedron_3;
		iMesh_createEnt( mesh,
				 iMesh_TETRAHEDRON,
				 connectivity,
				 4,
				 &tetrahedron_3,
				 &tet_creation_status,
				 &error);
		assert( iBase_SUCCESS == error );
		assert( iBase_NEW == tet_creation_status );

		data_arr[0] = tet_data;
		data_arr[1] = tet_data;
		data_arr[2] = tet_data;
		iMesh_setDblArrData( mesh,
				     &tetrahedron_3,
				     1,
				     tet_tag,
				     data_arr,
				     3,
				     &error);
		assert( iBase_SUCCESS == error );
		
		tet_data += 1.0;

		iMesh_addEntToSet(mesh,
				  tetrahedron_3,
				  tet_set,
				  &error);
		assert( iBase_SUCCESS == error );

		// Tetrahedron 4.
		connectivity[0] = vertices[ v4 ];
		connectivity[1] = vertices[ v7 ];
		connectivity[2] = vertices[ v6 ];
		connectivity[3] = vertices[ v3 ];

		iBase_EntityHandle tetrahedron_4;
		iMesh_createEnt( mesh,
				 iMesh_TETRAHEDRON,
				 connectivity,
				 4,
				 &tetrahedron_4,
				 &tet_creation_status,
				 &error);
		assert( iBase_SUCCESS == error );
		assert( iBase_NEW == tet_creation_status );

		data_arr[0] = tet_data;
		data_arr[1] = tet_data;
		data_arr[2] = tet_data;
		iMesh_setDblArrData( mesh,
				     &tetrahedron_4,
				     1,
				     tet_tag,
				     data_arr,
				     3,
				     &error);
		assert( iBase_SUCCESS == error );
 		
		tet_data += 1.0;

		iMesh_addEntToSet( mesh,
				   tetrahedron_4,
				   tet_set,
				   &error);
		assert( iBase_SUCCESS == error );

		// Tetrahedron 5.
		connectivity[0] = vertices[ v3 ];
		connectivity[1] = vertices[ v1 ];
		connectivity[2] = vertices[ v6 ];
		connectivity[3] = vertices[ v4 ];

		iBase_EntityHandle tetrahedron_5;
		iMesh_createEnt( mesh,
				 iMesh_TETRAHEDRON,
				 connectivity,
				 4,
				 &tetrahedron_5,
				 &tet_creation_status,
				 &error);
		assert( iBase_SUCCESS == error );
		assert( iBase_NEW == tet_creation_status );

		data_arr[0] = tet_data;
		data_arr[1] = tet_data;
		data_arr[2] = tet_data;
		iMesh_setDblArrData( mesh,
				     &tetrahedron_5,
				     1,
				     tet_tag,
				     data_arr,
				     3,
				     &error);
		assert( iBase_SUCCESS == error );
		
		tet_data += 1.0;

		iMesh_addEntToSet( mesh,
				   tetrahedron_5,
				   tet_set,
				   &error);
		assert( iBase_SUCCESS == error );
	    }
	}
    }
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEUCHOS_UNIT_TEST( FEMInterpolate, constructor_test )
{
}

//---------------------------------------------------------------------------//
//                        end of tstFEMInterpolate.cpp
//---------------------------------------------------------------------------//
