//----------------------------------*-C++-*----------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file   mesh/test/tstPointQuery.cpp
 * \author Stuart Slattery
 * \brief  PointQuery class unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <cassert>

#include <PointQuery.hpp>

#include <iMesh.h>
#include <iBase.h>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_ArrayRCP.hpp>

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

TEUCHOS_UNIT_TEST( PointQuery, hex_query )
{
    // Create a hex-8 element.
    int error;
    iMesh_Instance mesh;
    iMesh_newMesh( "", &mesh, &error, 0);
    TEST_ASSERT( iBase_SUCCESS == error );

    double vtx_coords[24] = { 0,0,0, 1,0,0, 1,1,0, 0,1,0,
			      0,0,1, 1,0,1, 1,1,1, 0,1,1 };
    int num_verts = 8;
    int new_coords_size = 24;
    int new_vertex_handles_allocated = 8;
    int new_vertex_handles_size = 0;
    iBase_EntityHandle *vertex_handles = 0;
    iMesh_createVtxArr( mesh,
			num_verts,
			iBase_INTERLEAVED,
			vtx_coords,
			new_coords_size,
			&vertex_handles,
			&new_vertex_handles_allocated,
			&new_vertex_handles_size,
			&error );
    TEST_ASSERT( iBase_SUCCESS == error );

    int status = 0;
    iBase_EntityHandle hex_element;
    iMesh_createEnt( mesh,
		     iMesh_HEXAHEDRON,
		     vertex_handles,
		     new_vertex_handles_size,
		     &hex_element,
		     &status,
		     &error );  
    TEST_ASSERT( iBase_SUCCESS == error );

    double point1[3] = { 0.5, 0.5, 0.5 };
    TEST_ASSERT( FOOD::PointQuery::pointInRefElement( 
		     mesh, hex_element, point1 ) == true );

    double point2[3] = { 0.25, 1.5, -0.5 };
    TEST_ASSERT( FOOD::PointQuery::pointInRefElement( 
		     mesh, hex_element, point2 ) == false );
}

TEUCHOS_UNIT_TEST( PointQuery, quadratic_hex_query )
{
    // Create a hex-27 element.
    int error;
    iMesh_Instance mesh;
    iMesh_newMesh( "", &mesh, &error, 0);
    TEST_ASSERT( iBase_SUCCESS == error );

    double vtx_coords[81] = { 0,0,0, 1,0,0, 1,1,0, 0,1,0,
			      0,0,1, 1,0,1, 1,1,1, 0,1,1, // linear nodes
			      0.5,0,0, 1,0.5,0, 0.5,1,0, 0,0.5,0,
			      0,0,0.5, 1,0,0.5, 1,1,0.5, 0,1,0.5,
			      0.5,0,1, 1,0.5,1, 0.5,1,1, 0,0.5,1, // quad 1
			      0.5,0.5,0.5, // centroid
			      0.5,0.5,0, 0.5,0.5,1, 0,0.5,0.5, 1,0.5,0.5,
			      0.5,0,0.5, 0.5,1,0.5 }; // quad 2
    int num_verts = 27;
    int new_coords_size = 81;
    int new_vertex_handles_allocated = 27;
    int new_vertex_handles_size = 0;
    iBase_EntityHandle *vertex_handles = 0;
    iMesh_createVtxArr( mesh,
			num_verts,
			iBase_INTERLEAVED,
			vtx_coords,
			new_coords_size,
			&vertex_handles,
			&new_vertex_handles_allocated,
			&new_vertex_handles_size,
			&error );
    TEST_ASSERT( iBase_SUCCESS == error );

    int status = 0;
    iBase_EntityHandle hex_element;
    iMesh_createEnt( mesh,
		     iMesh_HEXAHEDRON,
		     vertex_handles,
		     new_vertex_handles_size,
		     &hex_element,
		     &status,
		     &error );  
    TEST_ASSERT( iBase_SUCCESS == error );
    TEST_ASSERT( iBase_NEW == status );

    double point1[3] = { 0.5, 0.5, 0.5 };
    TEST_ASSERT( FOOD::PointQuery::pointInRefElement( 
		     mesh, hex_element, point1 ) == true );

    double point2[3] = { 0.25, 1.5, -0.5 };
    TEST_ASSERT( FOOD::PointQuery::pointInRefElement( 
		     mesh, hex_element, point2 ) == false );
}

//---------------------------------------------------------------------------//
//                        end of tstPointQuery.cpp
//---------------------------------------------------------------------------//
