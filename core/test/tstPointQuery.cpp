//----------------------------------*-C++-*----------------------------------//
/*!
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

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEUCHOS_UNIT_TEST( PointQuery, hex_query )
{
    typedef Intrepid::FieldContainer<double> MDArray;

    // Create a hex mesh.
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

    MDArray point1(1,3);
    point1(0,0) = 0.5;
    point1(0,1) = 0.5;
    point1(0,2) = 0.5;
    TEST_ASSERT( FOOD::PointQuery::point_in_ref_element( 
		     mesh, hex_element, point1 ) == true );

    MDArray point2(1,3);
    point2(0,0) = 0.25;
    point2(0,1) = 1.5;
    point2(0,2) = -0.5;
    TEST_ASSERT( FOOD::PointQuery::point_in_ref_element( 
		     mesh, hex_element, point2 ) == false );
}

//---------------------------------------------------------------------------//
//                        end of tstPointQuery.cpp
//---------------------------------------------------------------------------//