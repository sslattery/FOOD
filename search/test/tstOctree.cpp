//----------------------------------*-C++-*----------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file   mesh/test/tstOctree.cpp
 * \author Stuart Slattery
 * \brief  Octree class unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include <Octree.hpp>

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

void create_hex_mesh(iMesh_Instance &mesh)
{
    // Create the mesh instance.
    int error;
    iMesh_newMesh("", &mesh, &error, 0);
    assert( iBase_SUCCESS == error );

    // Generate vertices.
    int num_i = 4;
    int num_j = 4;
    int num_k = 4;
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

    // Generate hexahedrons from vertices and add them to the hex set. 
    iBase_EntityHandle connectivity[8];
    int hex_creation_status = 0;
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

		iMesh_addEntToSet(mesh,
				  hexahedron,
				  hex_set,
				  &error);
		assert( iBase_SUCCESS == error );
	    }
	}
    }
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEUCHOS_UNIT_TEST( Octree, tree_build_and_search_test )
{
    int error;
    iMesh_Instance mesh;
    create_hex_mesh(mesh);

    iBase_EntitySetHandle root_set;
    iMesh_getRootSet(mesh, &root_set, &error);
    TEST_ASSERT( iBase_SUCCESS == error );
    
    FOOD::Octree octree( mesh, root_set, iBase_REGION, iMesh_HEXAHEDRON );
    octree.buildTree();

    double coords1[3];
    coords1[0] = 1.5;
    coords1[1] = 1.5;
    coords1[2] = 1.5;

    double coords2[3];
    coords2[0] = -1.4;
    coords2[1] = 2.6;
    coords2[2] = 7.34;

    iBase_EntityHandle found_hex = 0;

    TEST_ASSERT( octree.findPoint( found_hex, coords1 ) );

    TEST_ASSERT( !octree.findPoint( found_hex, coords2 ) );
}

//---------------------------------------------------------------------------//
//                        end of tstOctree.cpp
//---------------------------------------------------------------------------//
