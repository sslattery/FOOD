//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/test/tstiMesh.cpp
 * \author Stuart Slattery
 * \brief  iMesh unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include <iMesh.h>
#include <iBase.h>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>

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

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEUCHOS_UNIT_TEST( iMesh, hex_mesh_test )
{
    // Setup iMesh instance
    int error;
    iMesh_Instance mesh;
    iMesh_newMesh("", &mesh, &error, 0);
    TEST_ASSERT( iBase_SUCCESS == error );

    // Create a tag for the vertex elements.
    iBase_TagHandle vertex_tag;
    std::string vertex_tag_name = "vertex_tag";
    iMesh_createTag(mesh,
		    &vertex_tag_name[0],
		    1,
		    iBase_DOUBLE,
		    &vertex_tag,
		    &error,
		    (int) vertex_tag_name.size());
    TEST_ASSERT( iBase_SUCCESS == error );

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
    iMesh_createVtxArr(mesh,
		       (int) coords.size() / 3,
		       iBase_INTERLEAVED,
		       &coords[0],
		       (int) coords.size(),
		       &vertices,
		       &vertices_allocated,
		       &vertices_size,
		       &error);
    TEST_ASSERT( iBase_SUCCESS == error );
    TEST_ASSERT( vertices_allocated == num_i*num_j*num_k );
    TEST_ASSERT( vertices_size == num_i*num_j*num_k );

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

    iMesh_setDblArrData(mesh,
			vertices,
			vertices_size,
			vertex_tag,
			&vertex_tag_data[0],
			(int) vertex_tag_data.size(),
			&error);
    TEST_ASSERT( iBase_SUCCESS == error );

    // Create a vertex set and add the vertex array.
    iBase_EntitySetHandle vertex_set;
    iMesh_createEntSet(mesh, 1, &vertex_set, &error);
    TEST_ASSERT( iBase_SUCCESS == error );

    iMesh_addEntArrToSet(mesh, vertices, vertices_size, vertex_set, &error);
    TEST_ASSERT( iBase_SUCCESS == error );

    // Create a hex set.
    iBase_EntitySetHandle hex_set;
    iMesh_createEntSet(mesh, 1, &hex_set, &error);
    TEST_ASSERT( iBase_SUCCESS == error );

    // Create a tag for the hex elements.
    iBase_TagHandle hex_tag;
    std::string hex_tag_name = "hex_tag";
    iMesh_createTag(mesh,
		    &hex_tag_name[0],
		    1,
		    iBase_DOUBLE,
		    &hex_tag,
		    &error,
		    (int) hex_tag_name.size());
    TEST_ASSERT( iBase_SUCCESS == error );

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
		iMesh_createEnt(mesh,
				iMesh_HEXAHEDRON,
				connectivity,
				8,
				&hexahedron,
				&hex_creation_status,
				&error);
		TEST_ASSERT( iBase_SUCCESS == error );
		TEST_ASSERT( iBase_NEW == hex_creation_status );

		iMesh_setDblData(mesh,
				 hexahedron,
				 hex_tag,
				 hex_data,
				 &error);
		TEST_ASSERT( iBase_SUCCESS == error );
		
		hex_data += 1.0;

		iMesh_addEntToSet(mesh,
				  hexahedron,
				  hex_set,
				  &error);
		TEST_ASSERT( iBase_SUCCESS == error );
	    }
	}
    }

    // Get the root set to write to file.
    iBase_EntitySetHandle root_set;
    iMesh_getRootSet(mesh, &root_set, &error);
    TEST_ASSERT( iBase_SUCCESS == error );

    // Write mesh to file.
    std::string filename = "hex_mesh_test.vtk";
    iMesh_save(mesh,
	       root_set,
	       &filename[0],
	       "",
	       &error,
	       (int) filename.size(),
	       0);
    TEST_ASSERT( iBase_SUCCESS == error );
}

//---------------------------------------------------------------------------//
//                        end of tstiMesh.cpp
//---------------------------------------------------------------------------//
