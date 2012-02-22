//----------------------------------*-C++-*----------------------------------//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
//
/*!
 * \file   mesh/test/tstTensorField.cpp
 * \author Stuart Slattery
 * \brief  TensorField class unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <cassert>

#include <DiscretizationTypes.hpp>
#include <FieldTypes.hpp>
#include <Quantity.hpp>
#include <Unit.hpp>
#include <Domain.hpp>
#include <DFuncKernel.hpp>
#include <TensorTemplate.hpp>
#include <TensorField.hpp>

#include <iMesh.h>
#include <iBase.h>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ENull.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include <Intrepid_FieldContainer.hpp>

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//

bool softEqualityCheck( double d1, double d2 )
{
    double epsilon = 1.0e-12;
    double diff = d1 - d2;
    diff *= diff;
    diff = pow( diff, 0.5 );
    if ( diff < epsilon ) return true;
    else return false;
}
	 

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

TEUCHOS_UNIT_TEST( TensorField, constructor_test )
{
    // Create the domain for this field.
    int error;
    iMesh_Instance mesh;
    iMesh_newMesh("", &mesh, &error, 0);
    TEST_ASSERT( iBase_SUCCESS == error );

    iBase_EntitySetHandle root_set;
    iMesh_getRootSet(mesh, &root_set, &error);
    TEST_ASSERT( iBase_SUCCESS == error );
    
    Teuchos::RCP<FOOD::Domain> domain = Teuchos::rcp(
	new FOOD::Domain(mesh, root_set, FOOD::FOOD_MBCN) );

    // Create the quantity for this field.
    Teuchos::Tuple<int,7> numerator;
    Teuchos::Tuple<int,7> denominator;
    for (int i = 0; i < 7; ++i)
    {
	numerator[i] = i;
	denominator[i] = 6 - i;
    }

    Teuchos::RCP<FOOD::Quantity> quantity = Teuchos::rcp(
	new FOOD::Quantity(numerator, denominator, "FOO_QUANTITY") );

    // Create the units for this field.
    Teuchos::RCP<FOOD::Unit> unit = Teuchos::rcp(
	new FOOD::Unit(quantity, 1.4, 4.3, "FOO_UNIT") );

    // Create the tensor template for this field.
    Teuchos::RCP<FOOD::TensorTemplate> tensor_template = Teuchos::rcp(
	new FOOD::TensorTemplate(0, 1, FOOD::FOOD_REAL, quantity) );

    // Create a distribution function for this field.
    Teuchos::RCP< FOOD::DFuncKernel<double> > dfunckernel =
	Teuchos::rcp( new FOOD::DFuncKernel<double>( iBase_REGION,
						     iMesh_HEXAHEDRON,
			                             iBase_FACE,
						     iMesh_QUADRILATERAL,
						     FOOD::FOOD_CARTESIAN,
						     FOOD::FOOD_FEM,
						     FOOD::FOOD_HGRAD,
						     FOOD::FOOD_SHARDSCN,
						     1 ) );

    // Create the field and check basic accessors.
    FOOD::TensorField<double> field( getDefaultComm<int>(),
				     domain,
				     dfunckernel,
				     FOOD::FOOD_CARTESIAN, 
				     tensor_template,
				     unit,
				     "FOO_FIELD" );

    TEST_ASSERT( field.getComm() == getDefaultComm<int>() );
    TEST_ASSERT( field.getDomain() == domain );
    TEST_ASSERT( field.getCoordType() == FOOD::FOOD_CARTESIAN );
    TEST_ASSERT( field.getTensorTemplate() == tensor_template );
    TEST_ASSERT( field.getUnit() == unit );
    TEST_ASSERT( field.getName() == "FOO_FIELD" );
}

TEUCHOS_UNIT_TEST( TensorField, dof_hex_mesh_vertex_tag_test )
{
    // Create a hex mesh.
    int error;
    iMesh_Instance mesh;
    create_hex_mesh( mesh );

    // Generate the domain for the field on the root set.
    iBase_EntitySetHandle root_set;
    iMesh_getRootSet(mesh, &root_set, &error);
    TEST_ASSERT( iBase_SUCCESS == error );
    
    Teuchos::RCP<FOOD::Domain> domain = Teuchos::rcp(
	new FOOD::Domain(mesh, root_set, FOOD::FOOD_MBCN) );

    // Create the quantity for this field.
    Teuchos::Tuple<int,7> numerator;
    Teuchos::Tuple<int,7> denominator;
    for (int i = 0; i < 7; ++i)
    {
	numerator[i] = i;
	denominator[i] = 6 - i;
    }

    Teuchos::RCP<FOOD::Quantity> quantity = Teuchos::rcp(
	new FOOD::Quantity(numerator, denominator, "VERTEX_QUANTITY") );

    // Create the units for this field.
    Teuchos::RCP<FOOD::Unit> unit = Teuchos::rcp(
	new FOOD::Unit(quantity, 1.4, 4.3, "VERTEX_UNIT") );

    // Create the tensor template for this field. The hex vertices are tagged
    // with a scalar field.
    Teuchos::RCP<FOOD::TensorTemplate> tensor_template = Teuchos::rcp(
	new FOOD::TensorTemplate(0, 1, FOOD::FOOD_REAL, quantity) );

    // Create a distribution function kernel for the field.
    Teuchos::RCP< FOOD::DFuncKernel<double> > dfunckernel =
	Teuchos::rcp( new FOOD::DFuncKernel<double>( iBase_REGION,
						     iMesh_HEXAHEDRON,
						     iBase_VERTEX,
						     iMesh_POINT,
						     FOOD::FOOD_CARTESIAN, 
						     FOOD::FOOD_FEM,
						     FOOD::FOOD_HGRAD,
						     FOOD::FOOD_SHARDSCN,
						     1 ) );

    // Create the field.
    FOOD::TensorField<double> field( getDefaultComm<int>(),
				     domain,
				     dfunckernel,
				     FOOD::FOOD_CARTESIAN, 
				     tensor_template,
				     unit,
				     "VERTEX_FIELD" );

    // Get the vertex tag to attach to the field.
    std::string tag_name = "vertex_tag";
    iBase_TagHandle vertex_tag;
    iMesh_getTagHandle( mesh,
			&tag_name[0],
			&vertex_tag,
			&error,
			(int) tag_name.size() );
    TEST_ASSERT( iBase_SUCCESS == error );

    // Attach the field to the tag.
    field.attachToTagData( vertex_tag, error );
    TEST_ASSERT( iBase_SUCCESS == error );

    // Test the tag attachment.
    int myRank = getDefaultComm<int>()->getRank();
    int mySize = getDefaultComm<int>()->getSize();
    int num_vertices = 1331;

    TEST_ASSERT( (int) field.getDFMap()->getGlobalNumElements()
		 == num_vertices*mySize );
    TEST_ASSERT( (int) field.getDFView().size() == num_vertices );
    TEST_ASSERT( (int) field.getDFConstView().size() == num_vertices );

    for (int i = 0; i < num_vertices; ++i)
    {
	TEST_ASSERT( (int) field.getDFMap()->getNodeElementList()[i]
		     == myRank*mySize + i );
	TEST_ASSERT( field.getDFView()[i] == (double) i );
	TEST_ASSERT( field.getDFConstView()[i] == (double) i );
    }
}

TEUCHOS_UNIT_TEST( TensorField, dof_tet_mesh_region_tag_test )
{
    // Create a tet mesh.
    int error;
    iMesh_Instance mesh;
    create_tet_mesh( mesh );

    // Generate the domain for the field on the root set.
    iBase_EntitySetHandle root_set;
    iMesh_getRootSet(mesh, &root_set, &error);
    TEST_ASSERT( iBase_SUCCESS == error );
    
    Teuchos::RCP<FOOD::Domain> domain = Teuchos::rcp(
	new FOOD::Domain(mesh, root_set, FOOD::FOOD_MBCN) );

    // Create the quantity for this field.
    Teuchos::Tuple<int,7> numerator;
    Teuchos::Tuple<int,7> denominator;
    for (int i = 0; i < 7; ++i)
    {
	numerator[i] = i;
	denominator[i] = 6 - i;
    }

    Teuchos::RCP<FOOD::Quantity> quantity = Teuchos::rcp(
	new FOOD::Quantity(numerator, denominator, "TET_QUANTITY") );

    // Create the units for this field.
    Teuchos::RCP<FOOD::Unit> unit = Teuchos::rcp(
	new FOOD::Unit(quantity, 1.4, 4.3, "TET_UNIT") );

    // Create the tensor template for this field. The tet volumes are tagged
    // with a 3-vector field.
    Teuchos::RCP<FOOD::TensorTemplate> tensor_template = Teuchos::rcp(
	new FOOD::TensorTemplate(1, 3, FOOD::FOOD_REAL, quantity) );

    // Create a distribution function kernel for the field.
    Teuchos::RCP< FOOD::DFuncKernel<double> > dfunckernel =
	Teuchos::rcp( new FOOD::DFuncKernel<double>( iBase_REGION,
						     iMesh_TETRAHEDRON,
						     iBase_REGION,
						     iMesh_TETRAHEDRON,
						     FOOD::FOOD_CARTESIAN,
						     FOOD::FOOD_FEM,
						     FOOD::FOOD_HGRAD,
						     FOOD::FOOD_SHARDSCN,
						     1 ) );

    // Create the field.
    FOOD::TensorField<double> field( getDefaultComm<int>(),
				     domain,
				     dfunckernel,
				     FOOD::FOOD_CARTESIAN, 
				     tensor_template,
				     unit,
				     "TET_FIELD" );

    // Get the vertex tag to attach to the field.
    std::string tag_name = "tet_tag";
    iBase_TagHandle vertex_tag;
    iMesh_getTagHandle( mesh,
			&tag_name[0],
			&vertex_tag,
			&error,
			(int) tag_name.size() );
    TEST_ASSERT( iBase_SUCCESS == error );

    // Attach the field to the tag.
    field.attachToTagData( vertex_tag, error );
    TEST_ASSERT( iBase_SUCCESS == error );

    // Test the tag attachment.
    int myRank = getDefaultComm<int>()->getRank();
    int mySize = getDefaultComm<int>()->getSize();
    int num_tets = 5000;
    int num_dof = num_tets*tensor_template->getNumComponents();

    TEST_ASSERT( (int) field.getDFMap()->getGlobalNumElements() ==
		 num_dof*mySize );
    TEST_ASSERT( (int) field.getDFView().size() == num_dof );
    TEST_ASSERT( (int) field.getDFConstView().size() == num_dof );

    int i1 = 0;
    int i2 = 0;
    int i3 = 0;
    for (int i = 0; i < num_tets; ++i)
    {
	i1 = 3*i;
	i2 = 3*i + 1;
	i3 = 3*i + 2;
	TEST_ASSERT( (int) field.getDFMap()->getNodeElementList()[i1]
		     == myRank*mySize + i1 );
	TEST_ASSERT( (int) field.getDFMap()->getNodeElementList()[i2]
		     == myRank*mySize + i2 );
	TEST_ASSERT( (int) field.getDFMap()->getNodeElementList()[i3]
		     == myRank*mySize + i3 );

	TEST_ASSERT( field.getDFView()[i1] == (double) i );
	TEST_ASSERT( field.getDFView()[i2] == (double) i );
	TEST_ASSERT( field.getDFView()[i3] == (double) i );

	TEST_ASSERT( field.getDFConstView()[i1] == (double) i );
	TEST_ASSERT( field.getDFConstView()[i2] == (double) i );
	TEST_ASSERT( field.getDFConstView()[i3] == (double) i );
    }
}

TEUCHOS_UNIT_TEST( TensorField, dof_hex_mesh_region_array_test )
{
    typedef Intrepid::FieldContainer<double> MDArray;

    // Create a hex mesh.
    int error;
    iMesh_Instance mesh;
    create_hex_mesh( mesh );

    // Generate the domain for the field on the root set.
    iBase_EntitySetHandle root_set;
    iMesh_getRootSet(mesh, &root_set, &error);
    TEST_ASSERT( iBase_SUCCESS == error );
    
    Teuchos::RCP<FOOD::Domain> domain = Teuchos::rcp(
	new FOOD::Domain(mesh, root_set, FOOD::FOOD_MBCN) );

    // Create the quantity for this field.
    Teuchos::Tuple<int,7> numerator;
    Teuchos::Tuple<int,7> denominator;
    for (int i = 0; i < 7; ++i)
    {
	numerator[i] = i;
	denominator[i] = 6 - i;
    }

    Teuchos::RCP<FOOD::Quantity> quantity = Teuchos::rcp(
	new FOOD::Quantity(numerator, denominator, "HEX_QUANTITY") );

    // Create the units for this field.
    Teuchos::RCP<FOOD::Unit> unit = Teuchos::rcp(
	new FOOD::Unit(quantity, 1.4, 4.3, "HEX_UNIT") );

    // Create the tensor template for this field. The hex vertices are tagged
    // with a scalar field.
    Teuchos::RCP<FOOD::TensorTemplate> tensor_template = Teuchos::rcp(
	new FOOD::TensorTemplate(0, 1, FOOD::FOOD_REAL, quantity) );

    // Create a distribution function kernel for the field.
    Teuchos::RCP< FOOD::DFuncKernel<double> > dfunckernel =
	Teuchos::rcp( new FOOD::DFuncKernel<double>( iBase_REGION,
						     iMesh_HEXAHEDRON,
						     iBase_REGION,
						     iMesh_HEXAHEDRON,
						     FOOD::FOOD_CARTESIAN,
						     FOOD::FOOD_FEM,
						     FOOD::FOOD_HGRAD,
						     FOOD::FOOD_SHARDSCN,
						     1 ) );

    // Create the field.
    FOOD::TensorField<double> field( getDefaultComm<int>(),
				     domain,
				     dfunckernel,
				     FOOD::FOOD_CARTESIAN, 
				     tensor_template,
				     unit,
				     "HEX_FIELD" );

    // Generate a degrees of freedom array to attach to the field.
    int num_hex = 1000;
    Teuchos::ArrayRCP<double> dof_array(num_hex, 0.0);
    Teuchos::ArrayRCP<double>::iterator dof_iterator;
    double data = 0.0;
    for (dof_iterator = dof_array.begin(); 
	 dof_iterator != dof_array.end(); 
	 ++dof_iterator)
    {
	*dof_iterator = data;
	data += 1.0;
    }

    // Attach the field to the array.
    field.attachToArrayData( dof_array, iBase_INTERLEAVED, error );
    TEST_ASSERT( iBase_SUCCESS == error );

    // Test the array attachment.
    int myRank = getDefaultComm<int>()->getRank();
    int mySize = getDefaultComm<int>()->getSize();

    TEST_ASSERT( (int) field.getDFMap()->getGlobalNumElements()
		 == num_hex*mySize );
    TEST_ASSERT( (int) field.getDFView().size() == num_hex );
    TEST_ASSERT( (int) field.getDFConstView().size() == num_hex );

    for (int i = 0; i < num_hex; ++i)
    {
	TEST_ASSERT( (int) field.getDFMap()->getNodeElementList()[i]
		     == myRank*mySize + i );
	TEST_ASSERT( field.getDFView()[i] == (double) i );
	TEST_ASSERT( field.getDFConstView()[i] == (double) i );
    }

    // Check that the mesh got tagged.
    iBase_TagHandle field_tag = field.getDFTag();

    iBase_TagHandle test_tag;
    std::string tag_name = "HEX_FIELD";
    iMesh_getTagHandle( mesh,
			&tag_name[0],
			&test_tag,
			&error,
			(int) tag_name.size());
    TEST_ASSERT( iBase_SUCCESS == error );
    TEST_ASSERT( test_tag == field_tag );

    iBase_EntityHandle *dof_entities = 0;
    int entities_allocated = num_hex;
    int entities_size = 0;
    iMesh_getEntities( mesh,
		       root_set,
		       iBase_REGION,
		       iMesh_HEXAHEDRON,
		       &dof_entities,
		       &entities_allocated,
		       &entities_size,
		       &error );
    TEST_ASSERT( iBase_SUCCESS == error );

    int tag_values_allocated = num_hex*sizeof(double);
    int tag_values_size = 0;
    Teuchos::ArrayRCP<double> field_tag_data(num_hex, 0.0);
    iMesh_getArrData( mesh,
		      dof_entities,
		      num_hex,
		      field_tag,
		      &field_tag_data,
		      &tag_values_allocated,
		      &tag_values_size,
		      &error );
    TEST_ASSERT( iBase_SUCCESS == error );
    TEST_ASSERT( tag_values_allocated == tag_values_size );

    for (int i = 0; i < num_hex; ++i)
    {
	TEST_ASSERT( field_tag_data[i] == dof_array[i] );
	TEST_ASSERT( iBase_SUCCESS == error );
    }
}

TEUCHOS_UNIT_TEST( TensorField, hex_evaluation_test )
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

    // Generate the domain for the field on the root set.
    iBase_EntitySetHandle root_set;
    iMesh_getRootSet(mesh, &root_set, &error);
    TEST_ASSERT( iBase_SUCCESS == error );
    
    Teuchos::RCP<FOOD::Domain> domain = Teuchos::rcp(
	new FOOD::Domain(mesh, root_set, FOOD::FOOD_MBCN) );

    // Create the quantity for this field.
    Teuchos::Tuple<int,7> numerator;
    Teuchos::Tuple<int,7> denominator;
    for (int i = 0; i < 7; ++i)
    {
	numerator[i] = i;
	denominator[i] = 6 - i;
    }

    Teuchos::RCP<FOOD::Quantity> quantity = Teuchos::rcp(
	new FOOD::Quantity(numerator, denominator, "HEX_QUANTITY") );

    // Create the units for this field.
    Teuchos::RCP<FOOD::Unit> unit = Teuchos::rcp(
	new FOOD::Unit(quantity, 1.4, 4.3, "HEX_UNIT") );

    // Create the tensor template for this field. The hex vertices are tagged
    // with a scalar field.
    Teuchos::RCP<FOOD::TensorTemplate> tensor_template = Teuchos::rcp(
	new FOOD::TensorTemplate(0, 1, FOOD::FOOD_REAL, quantity) );

    // Create a distribution function kernel for the field.
    Teuchos::RCP< FOOD::DFuncKernel<double> > dfunckernel =
	Teuchos::rcp( new FOOD::DFuncKernel<double>( iBase_REGION,
						     iMesh_HEXAHEDRON,
						     iBase_VERTEX,
						     iMesh_POINT,
						     FOOD::FOOD_CARTESIAN, 
						     FOOD::FOOD_FEM,
						     FOOD::FOOD_HGRAD,
						     FOOD::FOOD_SHARDSCN,
						     1 ) );

    // Create the field.
    FOOD::TensorField<double> field( getDefaultComm<int>(),
				     domain,
				     dfunckernel,
				     FOOD::FOOD_CARTESIAN, 
				     tensor_template,
				     unit,
				     "HEX_FIELD" );

    // Attach the field to array data. These are nodal values but they are
    // bound to the hex, so we tag the hex with them.
    Teuchos::ArrayRCP<double> hex_dof1(8, 6.54);
    field.attachToArrayData( hex_dof1, iBase_INTERLEAVED, error );
    TEST_ASSERT( iBase_SUCCESS == error );

    // Check that we can get the DOF on the hex.
    MDArray dof_coeffs(1,8,1);
    dof_coeffs = field.getEntDF( hex_element, error );
    TEST_ASSERT( iBase_SUCCESS == error );
    for ( int n = 0; n < dfunckernel->getBasisCardinality(); ++n )
    {
	TEST_ASSERT( dof_coeffs(0,n,0) == 6.54 );
    }

    // Evaluate the basis at a set of coordinates in the hex element.
    MDArray eval_coords1(1,3);
    eval_coords1(0,0) = 0.5;
    eval_coords1(0,1) = 0.5;
    eval_coords1(0,2) = 0.5;

    MDArray dfunc_values1(1,1,1);
    field.evaluateDF( hex_element, eval_coords1, false, dfunc_values1 );
    TEST_ASSERT( dfunc_values1(0,0,0) == 6.54 );


    // Create new DOFs and attach again;
    Teuchos::ArrayRCP<double> hex_dof2(8, 0.0);
    hex_dof2[4] = 1.0;
    hex_dof2[5] = 1.0;
    hex_dof2[6] = 1.0;
    hex_dof2[7] = 1.0;
    field.attachToArrayData( hex_dof2, iBase_INTERLEAVED, error );
    TEST_ASSERT( iBase_SUCCESS == error );

    // Evaluate the second basis at a set of coordinates in the hex element.
    MDArray eval_coords2(2,3);
    eval_coords2(0,0) = 0.5;
    eval_coords2(0,1) = 0.5;
    eval_coords2(0,2) = 0.5;

    eval_coords2(1,0) = 0.75;
    eval_coords2(1,1) = 0.75;
    eval_coords2(1,2) = 0.75;

    MDArray dfunc_values2(1,2,1);
    field.evaluateDF( hex_element, eval_coords2, false, dfunc_values2 );
    TEST_ASSERT( dfunc_values2(0,0,0) == 0.5 );
    TEST_ASSERT( dfunc_values2(0,1,0) == 0.75 );
}

TEUCHOS_UNIT_TEST( TensorField, hex_vector_evaluation_test )
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

    // Generate the domain for the field on the root set.
    iBase_EntitySetHandle root_set;
    iMesh_getRootSet(mesh, &root_set, &error);
    TEST_ASSERT( iBase_SUCCESS == error );
    
    Teuchos::RCP<FOOD::Domain> domain = Teuchos::rcp(
	new FOOD::Domain(mesh, root_set, FOOD::FOOD_MBCN) );

    // Create the tensor template for this field. The hex vertices are tagged
    // with a 3-vector field.
    Teuchos::RCP<FOOD::TensorTemplate> tensor_template = Teuchos::rcp(
	new FOOD::TensorTemplate(1, 3, FOOD::FOOD_REAL, Teuchos::null) );

    // Create a distribution function kernel for the field.
    Teuchos::RCP< FOOD::DFuncKernel<double> > dfunckernel =
	Teuchos::rcp( new FOOD::DFuncKernel<double>( iBase_REGION,
						     iMesh_HEXAHEDRON,
						     iBase_VERTEX,
						     iMesh_POINT,
						     FOOD::FOOD_CARTESIAN, 
						     FOOD::FOOD_FEM,
						     FOOD::FOOD_HGRAD,
						     FOOD::FOOD_SHARDSCN,
						     1 ) );

    // Create the field.
    FOOD::TensorField<double> field( getDefaultComm<int>(),
				     domain,
				     dfunckernel,
				     FOOD::FOOD_CARTESIAN, 
				     tensor_template,
				     Teuchos::null,
				     "HEX_FIELD" );

    // Attach the field to array data. These are nodal values but they are
    // bound to the hex, so we tag the hex with them.
    Teuchos::ArrayRCP<double> hex_dof1(24, 6.54);
    field.attachToArrayData( hex_dof1, iBase_INTERLEAVED, error );
    TEST_ASSERT( iBase_SUCCESS == error );

    // Check that we can get the DOF on the hex.
    MDArray dof_coeffs(1,8,3);
    dof_coeffs = field.getEntDF( hex_element, error );
    TEST_ASSERT( iBase_SUCCESS == error );
    for ( int n = 0; n < dfunckernel->getBasisCardinality(); ++n )
    {
	for ( int c = 0; c < (int) tensor_template->getNumComponents(); ++c )
	{
	    TEST_ASSERT( dof_coeffs(0,n,c) == 6.54 );
	}
    }

    // Evaluate the basis at a set of coordinates in the hex element.
    MDArray eval_coords1(1,3);
    eval_coords1(0,0) = 0.5;
    eval_coords1(0,1) = 0.5;
    eval_coords1(0,2) = 0.5;

    MDArray dfunc_values1(1,1,3);
    field.evaluateDF( hex_element, eval_coords1, false, dfunc_values1 );
    TEST_ASSERT( dfunc_values1(0,0,0) == 6.54 );
    TEST_ASSERT( dfunc_values1(0,0,1) == 6.54 );
    TEST_ASSERT( dfunc_values1(0,0,2) == 6.54 );

    // Create new DOFs and attach again;
    Teuchos::ArrayRCP<double> hex_dof2(24, 0.0);
    for ( int i = 12; i < 24; ++i )
    {
	hex_dof2[i] = 1.0;
    }
    field.attachToArrayData( hex_dof2, iBase_INTERLEAVED, error );
    TEST_ASSERT( iBase_SUCCESS == error );

    // Evaluate the second basis at a set of coordinates in the hex element.
    MDArray eval_coords2(2,3);
    eval_coords2(0,0) = 0.5;
    eval_coords2(0,1) = 0.5;
    eval_coords2(0,2) = 0.5;

    eval_coords2(1,0) = 0.75;
    eval_coords2(1,1) = 0.75;
    eval_coords2(1,2) = 0.75;

    MDArray dfunc_values2(1,2,3);
    field.evaluateDF( hex_element, eval_coords2, false, dfunc_values2 );
    TEST_ASSERT( dfunc_values2(0,0,0) == 0.5 );
    TEST_ASSERT( dfunc_values2(0,0,1) == 0.5 );
    TEST_ASSERT( dfunc_values2(0,0,2) == 0.5 );
    TEST_ASSERT( dfunc_values2(0,1,0) == 0.75 );
    TEST_ASSERT( dfunc_values2(0,1,1) == 0.75 );
    TEST_ASSERT( dfunc_values2(0,1,2) == 0.75 );
}

TEUCHOS_UNIT_TEST( TensorField, quadratic_hex_evaluation_test )
{
    typedef Intrepid::FieldContainer<double> MDArray;

    // Create a hex-27 element. MBCN ordering.
    int error;
    iMesh_Instance mesh;
    iMesh_newMesh( "", &mesh, &error, 0);
    TEST_ASSERT( iBase_SUCCESS == error );

    double vtx_coords[81] = { 0,0,0, 1,0,0, 1,1,0, 0,1,0,
			      0,0,1, 1,0,1, 1,1,1, 0,1,1, // linear nodes
			      0.5,0,0, 1,0.5,0, 0.5,1,0, 0,0.5,0,
			      0,0,0.5, 1,0,0.5, 1,1,0.5, 0,1,0.5,
			      0.5,0,1, 1,0.5,1, 0.5,1,1, 0,0.5,1, // 20
			      0.5,0,0.5, 1,0.5,0.5, 0.5,1,0.5, 0,0.5,0.5,
			      0.5,0.5,0, 0.5,0.5,1, 0.5,0.5,0.5 }; // 27

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

    // Generate the domain for the field on the root set.
    iBase_EntitySetHandle root_set;
    iMesh_getRootSet(mesh, &root_set, &error);
    TEST_ASSERT( iBase_SUCCESS == error );
    
    Teuchos::RCP<FOOD::Domain> domain = Teuchos::rcp(
	new FOOD::Domain(mesh, root_set, FOOD::FOOD_MBCN) );

    // Create the tensor template for this field. The hex vertices are tagged
    // with a scalar field.
    Teuchos::RCP<FOOD::TensorTemplate> tensor_template = Teuchos::rcp(
	new FOOD::TensorTemplate(0, 1, FOOD::FOOD_REAL, Teuchos::null) );

    // Create a distribution function kernel for the field.
    Teuchos::RCP< FOOD::DFuncKernel<double> > dfunckernel =
	Teuchos::rcp( new FOOD::DFuncKernel<double>( iBase_REGION,
						     iMesh_HEXAHEDRON,
						     iBase_VERTEX,
						     iMesh_POINT,
						     FOOD::FOOD_CARTESIAN, 
						     FOOD::FOOD_FEM,
						     FOOD::FOOD_HGRAD,
						     FOOD::FOOD_SHARDSCN,
						     2 ) );

    // Create the field.
    FOOD::TensorField<double> field( getDefaultComm<int>(),
				     domain,
				     dfunckernel,
				     FOOD::FOOD_CARTESIAN, 
				     tensor_template,
				     Teuchos::null,
				     "HEX_FIELD" );

    // Attach the field to array data. These are nodal values but they are
    // bound to the hex, so we tag the hex with them.
    Teuchos::ArrayRCP<double> hex_dof1(27, 6.54);
    field.attachToArrayData( hex_dof1, iBase_INTERLEAVED, error );
    TEST_ASSERT( iBase_SUCCESS == error );

    // Check that we can get the DOF on the hex-27.
    MDArray dof_coeffs(1,27,1);
    dof_coeffs = field.getEntDF( hex_element, error );
    TEST_ASSERT( iBase_SUCCESS == error );
    for ( int n = 0; n < dfunckernel->getBasisCardinality(); ++n )
    {
	TEST_ASSERT( dof_coeffs(0,n,0) == 6.54 );
    }

    // Evaluate the basis at a set of coordinates in the hex-27 element.
    MDArray eval_coords1(1,3);
    eval_coords1(0,0) = 0.5;
    eval_coords1(0,1) = 0.5;
    eval_coords1(0,2) = 0.5;

    MDArray dfunc_values1(1,1,1);
    field.evaluateDF( hex_element, eval_coords1, false, dfunc_values1 );
    TEST_ASSERT( dfunc_values1(0,0,0) == 6.54 );

    // Create new DOFs and attach again;
    Teuchos::ArrayRCP<double> hex_dof2(27, 0.0);
    hex_dof2[4] = 1.0;
    hex_dof2[5] = 1.0;
    hex_dof2[6] = 1.0;
    hex_dof2[7] = 1.0;
    hex_dof2[12] = 0.5;
    hex_dof2[13] = 0.5;
    hex_dof2[14] = 0.5;
    hex_dof2[15] = 0.5;
    hex_dof2[16] = 1.0;
    hex_dof2[17] = 1.0;
    hex_dof2[18] = 1.0;
    hex_dof2[19] = 1.0;
    hex_dof2[20] = 0.5;
    hex_dof2[22] = 1.0;
    hex_dof2[23] = 0.5;
    hex_dof2[24] = 0.5;
    hex_dof2[25] = 0.5;
    hex_dof2[26] = 0.5;
    field.attachToArrayData( hex_dof2, iBase_INTERLEAVED, error );
    TEST_ASSERT( iBase_SUCCESS == error );

    // Evaluate the second basis at a set of coordinates in the hex element.
    MDArray eval_coords2(2,3);
    eval_coords2(0,0) = 0.5;
    eval_coords2(0,1) = 0.5;
    eval_coords2(0,2) = 0.5;

    eval_coords2(1,0) = 0.75;
    eval_coords2(1,1) = 0.75;
    eval_coords2(1,2) = 0.75;

    MDArray dfunc_values2(1,2,1);
    field.evaluateDF( hex_element, eval_coords2, false, dfunc_values2 );
    TEST_ASSERT( dfunc_values2(0,0,0) == 0.5 );
    TEST_ASSERT( dfunc_values2(0,1,0) == 0.75 );
}

TEUCHOS_UNIT_TEST( TensorField, hex_grad_eval_test )
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

    // Generate the domain for the field on the root set.
    iBase_EntitySetHandle root_set;
    iMesh_getRootSet(mesh, &root_set, &error);
    TEST_ASSERT( iBase_SUCCESS == error );
    
    Teuchos::RCP<FOOD::Domain> domain = Teuchos::rcp(
	new FOOD::Domain(mesh, root_set, FOOD::FOOD_MBCN) );

    // Create the tensor template for this field. The hex vertices are tagged
    // with a scalar field.
    Teuchos::RCP<FOOD::TensorTemplate> tensor_template = Teuchos::rcp(
	new FOOD::TensorTemplate(0, 1, FOOD::FOOD_REAL, Teuchos::null) );

    // Create a distribution function kernel for the field.
    Teuchos::RCP< FOOD::DFuncKernel<double> > dfunckernel =
	Teuchos::rcp( new FOOD::DFuncKernel<double>( iBase_REGION,
						     iMesh_HEXAHEDRON,
						     iBase_VERTEX,
						     iMesh_POINT,
						     FOOD::FOOD_CARTESIAN, 
						     FOOD::FOOD_FEM,
						     FOOD::FOOD_HGRAD,
						     FOOD::FOOD_SHARDSCN,
						     1 ) );

    // Create the field.
    FOOD::TensorField<double> field( getDefaultComm<int>(),
				     domain,
				     dfunckernel,
				     FOOD::FOOD_CARTESIAN, 
				     tensor_template,
				     Teuchos::null,
				     "HEX_FIELD" );

    // Attach the field to array data. These are nodal values but they are
    // bound to the hex, so we tag the hex with them.
    Teuchos::ArrayRCP<double> hex_dof1(8, 6.54);
    field.attachToArrayData( hex_dof1, iBase_INTERLEAVED, error );
    TEST_ASSERT( iBase_SUCCESS == error );

    // Check that we can get the DOF on the hex.
    MDArray dof_coeffs(1,8,1);
    dof_coeffs = field.getEntDF( hex_element, error );
    TEST_ASSERT( iBase_SUCCESS == error );
    for ( int n = 0; n < dfunckernel->getBasisCardinality(); ++n )
    {
	TEST_ASSERT( dof_coeffs(0,n,0) == 6.54 );
    }

    // Evaluate the basis at a set of coordinates in the hex element.
    MDArray eval_coords1(1,3);
    eval_coords1(0,0) = 0.5;
    eval_coords1(0,1) = 0.5;
    eval_coords1(0,2) = 0.5;

    MDArray dfunc_values1(1,1,1,3);
    field.evaluateGradDF( hex_element, eval_coords1, false, dfunc_values1 );
    TEST_ASSERT( softEqualityCheck( dfunc_values1(0,0,0,0) ,0.0 ) );
    TEST_ASSERT( softEqualityCheck( dfunc_values1(0,0,0,1) ,0.0 ) );
    TEST_ASSERT( softEqualityCheck( dfunc_values1(0,0,0,2) ,0.0 ) );

    // Create new DOFs and attach again;
    Teuchos::ArrayRCP<double> hex_dof2(8, 0.0);
    hex_dof2[4] = 1.0;
    hex_dof2[5] = 1.0;
    hex_dof2[6] = 1.0;
    hex_dof2[7] = 1.0;
    field.attachToArrayData( hex_dof2, iBase_INTERLEAVED, error );
    TEST_ASSERT( iBase_SUCCESS == error );

    // Evaluate the second basis at a set of coordinates in the hex element.
    MDArray eval_coords2(2,3);
    eval_coords2(0,0) = 0.5;
    eval_coords2(0,1) = 0.5;
    eval_coords2(0,2) = 0.5;

    eval_coords2(1,0) = 0.75;
    eval_coords2(1,1) = 0.75;
    eval_coords2(1,2) = 0.75;

    MDArray dfunc_values2(1,2,1,3);
    field.evaluateGradDF( hex_element, eval_coords2, false, dfunc_values2 );
    TEST_ASSERT( softEqualityCheck( dfunc_values2(0,0,0,0), 0.0 ) );
    TEST_ASSERT( softEqualityCheck( dfunc_values2(0,0,0,1), 0.0 ) );
    TEST_ASSERT( softEqualityCheck( dfunc_values2(0,0,0,2), 1.0 ) );
    TEST_ASSERT( softEqualityCheck( dfunc_values2(0,1,0,0), 0.0 ) );
    TEST_ASSERT( softEqualityCheck( dfunc_values2(0,1,0,1), 0.0 ) );
    TEST_ASSERT( softEqualityCheck( dfunc_values2(0,1,0,2), 1.0 ) );
}

//---------------------------------------------------------------------------//
//                        end of tstTensorField.cpp
//---------------------------------------------------------------------------//
