//----------------------------------*-C++-*----------------------------------//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
//
/*!
 * \file   mesh/test/tstDomain.cpp
 * \author Stuart Slattery
 * \brief  Domain class unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include <DiscretizationTypes.hpp>
#include <Domain.hpp>

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

TEUCHOS_UNIT_TEST( Domain, constructor_test )
{
    int error;
    iMesh_Instance mesh;
    iMesh_newMesh("", &mesh, &error, 0);
    TEST_ASSERT( iBase_SUCCESS == error );

    iBase_EntitySetHandle root_set;
    iMesh_getRootSet(mesh, &root_set, &error);
    TEST_ASSERT( iBase_SUCCESS == error );
    
    FOOD::Domain domain(mesh, root_set, FOOD::FOOD_MBCN);
    TEST_ASSERT( domain.getMesh() == mesh );
    TEST_ASSERT( domain.getMeshSet() == root_set );
}

TEUCHOS_UNIT_TEST( Domain, mesh_iterator_test )
{
    int error;
    iMesh_Instance mesh;
    iMesh_newMesh("", &mesh, &error, 0);
    TEST_ASSERT( iBase_SUCCESS == error );

    iBase_EntityHandle *vertices;
    int vertices_allocated = 0;
    int vertices_size = 0;
    double coords[] = { 0,0,0, 1,1,1, 2,2,2, 3,3,3, 4,4,4, 5,5,5 };
    iMesh_createVtxArr(mesh,
		       6,
		       iBase_INTERLEAVED,
		       coords,
		       18,
		       &vertices,
		       &vertices_allocated,
		       &vertices_size,
		       &error);
    TEST_ASSERT( iBase_SUCCESS == error );
    TEST_ASSERT( vertices_allocated == 6 );
    TEST_ASSERT( vertices_size == 6 );

    iBase_EntitySetHandle vertex_set;
    iMesh_createEntSet(mesh, 1, &vertex_set, &error);
    TEST_ASSERT( iBase_SUCCESS == error );

    iMesh_addEntArrToSet(mesh, vertices, vertices_size, vertex_set, &error);
    TEST_ASSERT( iBase_SUCCESS == error );

    FOOD::Domain domain(mesh, vertex_set, FOOD::FOOD_MBCN);

    iBase_EntityIterator domain_vertex_iterator = 0; 
    error = domain.initEntIter(iBase_VERTEX, 
			       iMesh_POINT, 
			       &domain_vertex_iterator);
    TEST_ASSERT( iBase_SUCCESS == error );

    iBase_EntityIterator mesh_vertex_iterator = 0;
    iMesh_initEntIter(mesh, 
		      vertex_set, 
		      iBase_VERTEX,
		      iMesh_POINT,
		      0,
		      &mesh_vertex_iterator,
		      &error);
    TEST_ASSERT( iBase_SUCCESS == error );

    iBase_EntityArrIterator domain_vertex_arr_iterator = 0; 
    error = domain.initEntArrIter(iBase_VERTEX, 
				  iMesh_POINT, 
				  1, 
				  &domain_vertex_arr_iterator);
    TEST_ASSERT( iBase_SUCCESS == error );

    iBase_EntityArrIterator mesh_vertex_arr_iterator = 0;
    iMesh_initEntArrIter(mesh, 
			 vertex_set, 
			 iBase_VERTEX,
			 iMesh_POINT,
			 6,
			 0,
			 &mesh_vertex_arr_iterator,
			 &error);
    TEST_ASSERT( iBase_SUCCESS == error );
}

TEUCHOS_UNIT_TEST( Domain, domain_retrieval_test )
{
    int error = 0;
    iMesh_Instance mesh;
    iMesh_newMesh("", &mesh, &error, 0);
    TEST_ASSERT( iBase_SUCCESS == error );

    iBase_EntityHandle *vertices = 0;
    int vertices_allocated = 0;
    int vertices_size = 0;
    double coords[] = { 0,1,2,3,4,5, 0,1,2,3,4,5, 0,1,2,3,4,5};
    iMesh_createVtxArr(mesh,
		       6,
		       iBase_BLOCKED,
		       coords,
		       18,
		       &vertices,
		       &vertices_allocated,
		       &vertices_size,
		       &error);
    TEST_ASSERT( iBase_SUCCESS == error );
    TEST_ASSERT( vertices_allocated == 6 );
    TEST_ASSERT( vertices_size == 6 );

    iBase_EntitySetHandle vertex_set;
    iMesh_createEntSet(mesh, 1, &vertex_set, &error);
    TEST_ASSERT( iBase_SUCCESS == error );

    iMesh_addEntArrToSet(mesh, vertices, vertices_size, vertex_set, &error);
    TEST_ASSERT( iBase_SUCCESS == error );

    FOOD::Domain domain(mesh, vertex_set, FOOD::FOOD_MBCN);

    iBase_EntityArrIterator domain_vertex_arr_iterator = 0; 
    error = domain.initEntArrIter(iBase_VERTEX, 
				  iMesh_POINT, 
				  6, 
				  &domain_vertex_arr_iterator);
    TEST_ASSERT( iBase_SUCCESS == error );

    iBase_EntityHandle *domain_vertices = 0;
    int domain_vertices_allocated = 6;
    int domain_vertices_size = 0;
    int domain_vertices_have_data = 0;
    iMesh_getNextEntArrIter(mesh,
			    domain_vertex_arr_iterator,
			    &domain_vertices,
			    &domain_vertices_allocated,
			    &domain_vertices_size,
			    &domain_vertices_have_data,
			    &error);
    TEST_ASSERT( iBase_SUCCESS == error );
    TEST_ASSERT( domain_vertices_allocated == 6);
    TEST_ASSERT( domain_vertices_size == 6);
    TEST_ASSERT( domain_vertices_have_data != 0);

    double *domain_coords = 0;
    int domain_coords_allocated = 3*domain_vertices_size;
    int domain_coords_size = 0;
    iMesh_getVtxArrCoords(mesh,
			  domain_vertices,
			  domain_vertices_size,
			  iBase_INTERLEAVED,
			  &domain_coords,
			  &domain_coords_allocated,
			  &domain_coords_size,
			  &error);
    TEST_ASSERT( iBase_SUCCESS == error );
    TEST_ASSERT( domain_coords_allocated == 18 );
    TEST_ASSERT( domain_coords_size == 18 );

    for (int i = 0; i < 6; ++i)
    {
	TEST_ASSERT( domain_coords[3*i]   == (double) i );
	TEST_ASSERT( domain_coords[3*i+1] == (double) i );
	TEST_ASSERT( domain_coords[3*i+2] == (double) i );
    }
}

//---------------------------------------------------------------------------//
//                        end of tstDomain.cpp
//---------------------------------------------------------------------------//
