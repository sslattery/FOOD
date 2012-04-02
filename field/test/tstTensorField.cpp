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
#include <DFuncKernelFactory.hpp>
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
	new FOOD::TensorTemplate(0, 1, FOOD::FOOD_REAL) );

    // Create a distribution function for this field.
    FOOD::DFuncKernelFactory<double> kernel_factory;
    Teuchos::RCP< FOOD::DFuncKernel<double> > dfunckernel = 
	kernel_factory.create( iBase_REGION,
			       iMesh_HEXAHEDRON,
			       FOOD::FOOD_FEM,
			       FOOD::FOOD_HGRAD,
			       1 );

    // Create the field and check basic accessors.
    FOOD::TensorField<double> field( domain,
				     dfunckernel,
				     tensor_template,
				     "FOO_FIELD" );

    field.setUnit( unit );

    TEST_ASSERT( field.getDomain() == domain );
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

    // Create the tensor template for this field. The hex vertices are tagged
    // with a scalar field.
    Teuchos::RCP<FOOD::TensorTemplate> tensor_template = Teuchos::rcp(
	new FOOD::TensorTemplate(0, 1, FOOD::FOOD_REAL) );

    // Create a distribution function kernel for the field.
    FOOD::DFuncKernelFactory<double> kernel_factory;
    Teuchos::RCP< FOOD::DFuncKernel<double> > dfunckernel = 
	kernel_factory.create( iBase_REGION,
			       iMesh_HEXAHEDRON,
			       FOOD::FOOD_FEM,
			       FOOD::FOOD_HGRAD,
			       1 );

    // Create the field.
    FOOD::TensorField<double> field( domain,
				     dfunckernel,
				     tensor_template,
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
    int num_vertices = 1331;

    TEST_ASSERT( (int) field.getDFView().size() == num_vertices );
    TEST_ASSERT( (int) field.getDFConstView().size() == num_vertices );

    for (int i = 0; i < num_vertices; ++i)
    {
	TEST_ASSERT( field.getDFView()[i] == (double) i );
	TEST_ASSERT( field.getDFConstView()[i] == (double) i );
    }
}

TEUCHOS_UNIT_TEST( TensorField, hex_evaluation_test )
{
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
	new FOOD::TensorTemplate(0, 1, FOOD::FOOD_REAL) );

    // Create a distribution function kernel for the field.
    FOOD::DFuncKernelFactory<double> kernel_factory;
    Teuchos::RCP< FOOD::DFuncKernel<double> > dfunckernel = 
	kernel_factory.create( iBase_REGION,
			       iMesh_HEXAHEDRON,
			       FOOD::FOOD_FEM,
			       FOOD::FOOD_HGRAD,
			       1 );

    // Create the field.
    FOOD::TensorField<double> field( domain,
				     dfunckernel,
				     tensor_template,
				     "HEX_FIELD" );

    // Attach the field to array data. These are nodal values but they are
    // bound to the hex, so we tag the hex with them.
    Teuchos::ArrayRCP<double> hex_dof1(8, 6.54);
    field.attachToArrayData( hex_dof1, iBase_INTERLEAVED, error );
    TEST_ASSERT( iBase_SUCCESS == error );

    // Check that we can get the DOF on the hex.
    Teuchos::ArrayRCP<double> dof_coeffs(8);
    dof_coeffs = field.getEntDF( hex_element );
    for ( int n = 0; n < dfunckernel->getCardinality(); ++n )
    {
	TEST_ASSERT( dof_coeffs[n] == 6.54 );
    }

    // Evaluate the basis at a set of coordinates in the hex element.
    double eval_coords1[3] = { 0.5, 0.5, 0.5 };
    Teuchos::ArrayRCP<double> dfunc_values1(1);
    field.evaluateDF( hex_element, eval_coords1, false, dfunc_values1 );
    TEST_ASSERT( dfunc_values1[0] == 6.54 );

    // Create new DOFs and attach again;
    Teuchos::ArrayRCP<double> hex_dof2(8, 0.0);
    hex_dof2[4] = 1.0;
    hex_dof2[5] = 1.0;
    hex_dof2[6] = 1.0;
    hex_dof2[7] = 1.0;
    field.attachToArrayData( hex_dof2, iBase_INTERLEAVED, error );
    TEST_ASSERT( iBase_SUCCESS == error );

    // Evaluate the second basis at a set of coordinates in the hex element.
    double eval_coords2[3] = { 0.5, 0.5, 0.5 };
    double eval_coords3[3] = { 0.75, 0.75, 0.75 };

    Teuchos::ArrayRCP<double> dfunc_values2(1);
    Teuchos::ArrayRCP<double> dfunc_values3(1);
    field.evaluateDF( hex_element, eval_coords2, false, dfunc_values2 );
    field.evaluateDF( hex_element, eval_coords3, false, dfunc_values3 );
    TEST_ASSERT( dfunc_values2[0] == 0.5 );
    TEST_ASSERT( dfunc_values3[0] == 0.75 );
}

// TEUCHOS_UNIT_TEST( TensorField, hex_vector_evaluation_test )
// {
//     // Create a hex mesh.
//     int error;
//     iMesh_Instance mesh;
//     iMesh_newMesh( "", &mesh, &error, 0);
//     TEST_ASSERT( iBase_SUCCESS == error );

//     double vtx_coords[24] = { 0,0,0, 1,0,0, 1,1,0, 0,1,0,
// 			      0,0,1, 1,0,1, 1,1,1, 0,1,1 };
//     int num_verts = 8;
//     int new_coords_size = 24;
//     int new_vertex_handles_allocated = 8;
//     int new_vertex_handles_size = 0;
//     iBase_EntityHandle *vertex_handles = 0;
//     iMesh_createVtxArr( mesh,
// 			num_verts,
// 			iBase_INTERLEAVED,
// 			vtx_coords,
// 			new_coords_size,
// 			&vertex_handles,
// 			&new_vertex_handles_allocated,
// 			&new_vertex_handles_size,
// 			&error );
//     TEST_ASSERT( iBase_SUCCESS == error );

//     int status = 0;
//     iBase_EntityHandle hex_element;
//     iMesh_createEnt( mesh,
// 		     iMesh_HEXAHEDRON,
// 		     vertex_handles,
// 		     new_vertex_handles_size,
// 		     &hex_element,
// 		     &status,
// 		     &error );  
//     TEST_ASSERT( iBase_SUCCESS == error );

//     // Generate the domain for the field on the root set.
//     iBase_EntitySetHandle root_set;
//     iMesh_getRootSet(mesh, &root_set, &error);
//     TEST_ASSERT( iBase_SUCCESS == error );
    
//     Teuchos::RCP<FOOD::Domain> domain = Teuchos::rcp(
// 	new FOOD::Domain(mesh, root_set, FOOD::FOOD_MBCN) );

//     // Create the tensor template for this field. The hex vertices are tagged
//     // with a 3-vector field.
//     Teuchos::RCP<FOOD::TensorTemplate> tensor_template = Teuchos::rcp(
// 	new FOOD::TensorTemplate(1, 3, FOOD::FOOD_REAL) );

//     // Create a distribution function kernel for the field.
//     FOOD::DFuncKernelFactory<double> kernel_factory;
//     Teuchos::RCP< FOOD::DFuncKernel<double> > dfunckernel = 
// 	kernel_factory.create( iBase_REGION,
// 			       iMesh_HEXAHEDRON,
// 			       FOOD::FOOD_FEM,
// 			       FOOD::FOOD_HGRAD,
// 			       1 );

//     // Create the field.
//     FOOD::TensorField<double> field( domain,
// 				     dfunckernel,
// 				     tensor_template,
// 				     "HEX_FIELD" );

//     // Attach the field to array data. These are nodal values but they are
//     // bound to the hex, so we tag the hex with them.
//     Teuchos::ArrayRCP<double> hex_dof1(24, 6.54);
//     field.attachToArrayData( hex_dof1, iBase_INTERLEAVED, error );
//     TEST_ASSERT( iBase_SUCCESS == error );

//     // Check that we can get the DOF on the hex.
//     Teuchos::ArrayRCP<double> dof_coeffs(24);
//     dof_coeffs = field.getEntDF( hex_element );
//     for ( int n = 0; n < 24; ++n )
//     {
// 	TEST_ASSERT( dof_coeffs[n] == 6.54 );
//     }

//     // Evaluate the basis at a set of coordinates in the hex element.
//     double eval_coords1[3] = { 0.5, 0.5, 0.5 };

//     Teuchos::ArrayRCP<double> dfunc_values1(3);
//     field.evaluateDF( hex_element, eval_coords1, false, dfunc_values1 );
//     TEST_ASSERT( dfunc_values1[0] == 6.54 );
//     TEST_ASSERT( dfunc_values1[1] == 6.54 );
//     TEST_ASSERT( dfunc_values1[2] == 6.54 );

//     // Create new DOFs and attach again;
//     Teuchos::ArrayRCP<double> hex_dof2(24, 0.0);
//     for ( int i = 12; i < 24; ++i )
//     {
// 	hex_dof2[i] = 1.0;
//     }
//     field.attachToArrayData( hex_dof2, iBase_INTERLEAVED, error );

//     // Evaluate the second basis at a set of coordinates in the hex element.
//     double eval_coords2[3] = { 0.5, 0.5, 0.5 };
//     double eval_coords3[3] = { 0.75, 0.75, 0.75 };

//     Teuchos::ArrayRCP<double> dfunc_values2;
//     Teuchos::ArrayRCP<double> dfunc_values3;
//     field.evaluateDF( hex_element, eval_coords2, false, dfunc_values2 );
//     field.evaluateDF( hex_element, eval_coords3, false, dfunc_values3 );

//     std::cout << dfunc_values2[0] << std::endl;
//     std::cout << dfunc_values2[1] << std::endl;
//     std::cout << dfunc_values2[2] << std::endl;
//     std::cout << dfunc_values3[0] << std::endl;
//     std::cout << dfunc_values3[1] << std::endl;
//     std::cout << dfunc_values3[2] << std::endl;
//     TEST_ASSERT( dfunc_values2[0] == 0.5 );
//     TEST_ASSERT( dfunc_values2[1] == 0.5 );
//     TEST_ASSERT( dfunc_values2[2] == 0.5 );
//     TEST_ASSERT( dfunc_values3[0] == 0.75 );
//     TEST_ASSERT( dfunc_values3[1] == 0.75 );
//     TEST_ASSERT( dfunc_values3[2] == 0.75 );
// }

// TEUCHOS_UNIT_TEST( TensorField, quadratic_hex_evaluation_test )
// {
//     // Create a hex-27 element. MBCN ordering.
//     int error;
//     iMesh_Instance mesh;
//     iMesh_newMesh( "", &mesh, &error, 0);
//     TEST_ASSERT( iBase_SUCCESS == error );

//     double vtx_coords[81] = { 0,0,0, 1,0,0, 1,1,0, 0,1,0,
// 			      0,0,1, 1,0,1, 1,1,1, 0,1,1, // linear nodes
// 			      0.5,0,0, 1,0.5,0, 0.5,1,0, 0,0.5,0,
// 			      0,0,0.5, 1,0,0.5, 1,1,0.5, 0,1,0.5,
// 			      0.5,0,1, 1,0.5,1, 0.5,1,1, 0,0.5,1, // 20
// 			      0.5,0,0.5, 1,0.5,0.5, 0.5,1,0.5, 0,0.5,0.5,
// 			      0.5,0.5,0, 0.5,0.5,1, 0.5,0.5,0.5 }; // 27

//     int num_verts = 27;
//     int new_coords_size = 81;
//     int new_vertex_handles_allocated = 27;
//     int new_vertex_handles_size = 0;
//     iBase_EntityHandle *vertex_handles = 0;
//     iMesh_createVtxArr( mesh,
// 			num_verts,
// 			iBase_INTERLEAVED,
// 			vtx_coords,
// 			new_coords_size,
// 			&vertex_handles,
// 			&new_vertex_handles_allocated,
// 			&new_vertex_handles_size,
// 			&error );
//     TEST_ASSERT( iBase_SUCCESS == error );

//     int status = 0;
//     iBase_EntityHandle hex_element;
//     iMesh_createEnt( mesh,
// 		     iMesh_HEXAHEDRON,
// 		     vertex_handles,
// 		     new_vertex_handles_size,
// 		     &hex_element,
// 		     &status,
// 		     &error );  
//     TEST_ASSERT( iBase_SUCCESS == error );
//     TEST_ASSERT( iBase_NEW == status );

//     // Generate the domain for the field on the root set.
//     iBase_EntitySetHandle root_set;
//     iMesh_getRootSet(mesh, &root_set, &error);
//     TEST_ASSERT( iBase_SUCCESS == error );
    
//     Teuchos::RCP<FOOD::Domain> domain = Teuchos::rcp(
// 	new FOOD::Domain(mesh, root_set, FOOD::FOOD_MBCN) );

//     // Create the tensor template for this field. The hex vertices are tagged
//     // with a scalar field.
//     Teuchos::RCP<FOOD::TensorTemplate> tensor_template = Teuchos::rcp(
// 	new FOOD::TensorTemplate(0, 1, FOOD::FOOD_REAL) );

//     // Create a distribution function kernel for the field.
//     FOOD::DFuncKernelFactory<double> kernel_factory;
//     Teuchos::RCP< FOOD::DFuncKernel<double> > dfunckernel = 
// 	kernel_factory.create( iBase_REGION,
// 			       iMesh_HEXAHEDRON,
// 			       FOOD::FOOD_FEM,
// 			       FOOD::FOOD_HGRAD,
// 			       1 );

//     // Create the field.
//     FOOD::TensorField<double> field( domain,
// 				     dfunckernel,
// 				     tensor_template,
// 				     "HEX_FIELD" );

//     // Attach the field to array data. These are nodal values but they are
//     // bound to the hex, so we tag the hex with them.
//     Teuchos::ArrayRCP<double> hex_dof1(27, 6.54);
//     field.attachToArrayData( hex_dof1, iBase_INTERLEAVED, error );
//     TEST_ASSERT( iBase_SUCCESS == error );

//     // Check that we can get the DOF on the hex-27.
//     Teuchos::ArrayRCP<double> dof_coeffs(27);
//     dof_coeffs = field.getEntDF( hex_element );
//     for ( int n = 0; n < dfunckernel->getCardinality(); ++n )
//     {
// 	TEST_ASSERT( dof_coeffs[n] == 6.54 );
//     }

//     // Evaluate the basis at a set of coordinates in the hex-27 element.
//     double eval_coords1[3] = { 0.5, 0.5, 0.5 };
//     Teuchos::ArrayRCP<double> dfunc_values1(1);
//     field.evaluateDF( hex_element, eval_coords1, false, dfunc_values1 );
//     TEST_ASSERT( dfunc_values1[0] == 6.54 );

//     // Create new DOFs and attach again;
//     Teuchos::ArrayRCP<double> hex_dof2(27, 0.0);
//     hex_dof2[4] = 1.0;
//     hex_dof2[5] = 1.0;
//     hex_dof2[6] = 1.0;
//     hex_dof2[7] = 1.0;
//     hex_dof2[12] = 0.5;
//     hex_dof2[13] = 0.5;
//     hex_dof2[14] = 0.5;
//     hex_dof2[15] = 0.5;
//     hex_dof2[16] = 1.0;
//     hex_dof2[17] = 1.0;
//     hex_dof2[18] = 1.0;
//     hex_dof2[19] = 1.0;
//     hex_dof2[20] = 0.5;
//     hex_dof2[22] = 1.0;
//     hex_dof2[23] = 0.5;
//     hex_dof2[24] = 0.5;
//     hex_dof2[25] = 0.5;
//     hex_dof2[26] = 0.5;
//     field.attachToArrayData( hex_dof2, iBase_INTERLEAVED, error );
//     TEST_ASSERT( iBase_SUCCESS == error );

//     // Evaluate the second basis at a set of coordinates in the hex element.
//     double eval_coords2[3] = { 0.5, 0.5, 0.5 };
//     double eval_coords3[3] = { 0.75, 0.75, 0.75 };

//     Teuchos::ArrayRCP<double> dfunc_values2(1);
//     Teuchos::ArrayRCP<double> dfunc_values3(1);
//     field.evaluateDF( hex_element, eval_coords2, false, dfunc_values2 );
//     field.evaluateDF( hex_element, eval_coords3, false, dfunc_values3 );

//     TEST_ASSERT( dfunc_values2[0] == 0.5 );
//     TEST_ASSERT( dfunc_values2[1] == 0.5 );
//     TEST_ASSERT( dfunc_values2[2] == 0.5 );
//     TEST_ASSERT( dfunc_values3[0] == 0.75 );
//     TEST_ASSERT( dfunc_values3[1] == 0.75 );
//     TEST_ASSERT( dfunc_values3[2] == 0.75 );
// }

// TEUCHOS_UNIT_TEST( TensorField, hex_grad_eval_test )
// {
//     // Create a hex mesh.
//     int error;
//     iMesh_Instance mesh;
//     iMesh_newMesh( "", &mesh, &error, 0);
//     TEST_ASSERT( iBase_SUCCESS == error );

//     double vtx_coords[24] = { 0,0,0, 1,0,0, 1,1,0, 0,1,0,
// 			      0,0,1, 1,0,1, 1,1,1, 0,1,1 };
//     int num_verts = 8;
//     int new_coords_size = 24;
//     int new_vertex_handles_allocated = 8;
//     int new_vertex_handles_size = 0;
//     iBase_EntityHandle *vertex_handles = 0;
//     iMesh_createVtxArr( mesh,
// 			num_verts,
// 			iBase_INTERLEAVED,
// 			vtx_coords,
// 			new_coords_size,
// 			&vertex_handles,
// 			&new_vertex_handles_allocated,
// 			&new_vertex_handles_size,
// 			&error );
//     TEST_ASSERT( iBase_SUCCESS == error );

//     int status = 0;
//     iBase_EntityHandle hex_element;
//     iMesh_createEnt( mesh,
// 		     iMesh_HEXAHEDRON,
// 		     vertex_handles,
// 		     new_vertex_handles_size,
// 		     &hex_element,
// 		     &status,
// 		     &error );  
//     TEST_ASSERT( iBase_SUCCESS == error );

//     // Generate the domain for the field on the root set.
//     iBase_EntitySetHandle root_set;
//     iMesh_getRootSet(mesh, &root_set, &error);
//     TEST_ASSERT( iBase_SUCCESS == error );
    
//     Teuchos::RCP<FOOD::Domain> domain = Teuchos::rcp(
// 	new FOOD::Domain(mesh, root_set, FOOD::FOOD_MBCN) );

//     // Create the tensor template for this field. The hex vertices are tagged
//     // with a scalar field.
//     Teuchos::RCP<FOOD::TensorTemplate> tensor_template = Teuchos::rcp(
// 	new FOOD::TensorTemplate(0, 1, FOOD::FOOD_REAL) );

//     // Create a distribution function kernel for the field.
//     FOOD::DFuncKernelFactory<double> kernel_factory;
//     Teuchos::RCP< FOOD::DFuncKernel<double> > dfunckernel = 
// 	kernel_factory.create( iBase_REGION,
// 			       iMesh_HEXAHEDRON,
// 			       FOOD::FOOD_FEM,
// 			       FOOD::FOOD_HGRAD,
// 			       1 );

//     // Create the field.
//     FOOD::TensorField<double> field( domain,
// 				     dfunckernel,
// 				     tensor_template,
// 				     "HEX_FIELD" );

//     // Attach the field to array data. These are nodal values but they are
//     // bound to the hex, so we tag the hex with them.
//     Teuchos::ArrayRCP<double> hex_dof1(8, 6.54);
//     field.attachToArrayData( hex_dof1, iBase_INTERLEAVED, error );
//     TEST_ASSERT( iBase_SUCCESS == error );

//     // Check that we can get the DOF on the hex.
//     Teuchos::ArrayRCP<double> dof_coeffs(8);
//     dof_coeffs = field.getEntDF( hex_element );
//     for ( int n = 0; n < dfunckernel->getCardinality(); ++n )
//     {
// 	TEST_ASSERT( dof_coeffs[n] == 6.54 );
//     }

//     // Evaluate the basis at a set of coordinates in the hex element.
//     double eval_coords1[3] = { 0.5, 0.5, 0.5 };
//     Teuchos::ArrayRCP<double> dfunc_values1(3);
//     field.evaluateGradDF( hex_element, eval_coords1, false, dfunc_values1 );
//     TEST_ASSERT( softEqualityCheck( dfunc_values1[0] ,0.0 ) );
//     TEST_ASSERT( softEqualityCheck( dfunc_values1[1] ,0.0 ) );
//     TEST_ASSERT( softEqualityCheck( dfunc_values1[2] ,0.0 ) );

//     // Create new DOFs and attach again;
//     Teuchos::ArrayRCP<double> hex_dof2(8, 0.0);
//     hex_dof2[4] = 1.0;
//     hex_dof2[5] = 1.0;
//     hex_dof2[6] = 1.0;
//     hex_dof2[7] = 1.0;
//     field.attachToArrayData( hex_dof2, iBase_INTERLEAVED, error );
//     TEST_ASSERT( iBase_SUCCESS == error );

//     // Evaluate the second basis at a set of coordinates in the hex element.
//     double eval_coords2[3] = { 0.5, 0.5, 0.5 };
//     double eval_coords3[3] = { 0.75, 0.75, 0.75 };

//     Teuchos::ArrayRCP<double> dfunc_values2(1);
//     Teuchos::ArrayRCP<double> dfunc_values3(1);
//     field.evaluateDF( hex_element, eval_coords2, false, dfunc_values2 );
//     field.evaluateDF( hex_element, eval_coords3, false, dfunc_values3 );

//     TEST_ASSERT( dfunc_values2[0] == 0.5 );
//     TEST_ASSERT( dfunc_values2[1] == 0.5 );
//     TEST_ASSERT( dfunc_values2[2] == 0.5 );
//     TEST_ASSERT( dfunc_values3[0] == 0.75 );
//     TEST_ASSERT( dfunc_values3[1] == 0.75 );
//     TEST_ASSERT( dfunc_values3[2] == 0.75 );
// }

TEUCHOS_UNIT_TEST( TensorField, integral_test )
{
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
	new FOOD::TensorTemplate(0, 1, FOOD::FOOD_REAL) );

    // Create a distribution function kernel for the field.
    FOOD::DFuncKernelFactory<double> kernel_factory;
    Teuchos::RCP< FOOD::DFuncKernel<double> > dfunckernel = 
	kernel_factory.create( iBase_REGION,
			       iMesh_HEXAHEDRON,
			       FOOD::FOOD_FEM,
			       FOOD::FOOD_HGRAD,
			       1 );

    // Create the field.
    FOOD::TensorField<double> field( domain,
				     dfunckernel,
				     tensor_template,
				     "HEX_FIELD" );

    // Attach the field to array data. These are nodal values but they are
    // bound to the hex, so we tag the hex with them.
    Teuchos::ArrayRCP<double> hex_dof1(8, 6.54);
    field.attachToArrayData( hex_dof1, iBase_INTERLEAVED, error );
    TEST_ASSERT( iBase_SUCCESS == error );

    // Integrate the field over the hex.
    field.integrateDF();

    // Check that the integral was tagged to the mesh.
    std::string tag_name = "HEX_FIELD_INTEGRAL";
    iBase_TagHandle integral_tag;
    iMesh_getTagHandle( mesh,
			&tag_name[0],
			&integral_tag,
			&error,
			(int) tag_name.size() );
    TEST_ASSERT( iBase_SUCCESS == error );
    
    // Check the tag data.
    double data = 0.0;
    iMesh_getDblData( mesh,
		      hex_element,
		      integral_tag,
		      &data,
		      &error );
    TEST_ASSERT( iBase_SUCCESS == error );
    TEST_ASSERT( data == 6.54 );
    std::cout << "DATA " << data << std::endl;
}

//---------------------------------------------------------------------------//
//                        end of tstTensorField.cpp
//---------------------------------------------------------------------------//
