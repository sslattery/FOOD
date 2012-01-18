//----------------------------------*-C++-*----------------------------------//
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

#include <Types.hpp>
#include <Quantity.hpp>
#include <Unit.hpp>
#include <TensorTemplate.hpp>
#include <TensorField.hpp>

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
	new FOOD::Domain(mesh, root_set) );

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
	new FOOD::TensorTemplate(0, 1, FOOD::REAL, quantity) );

    // Create the field and check basic accessors.
    FOOD::TensorField<double> field( getDefaultComm<int>(),
				     domain,
				     iBase_VERTEX,
				     iMesh_POINT,
				     FOOD::CARTESIAN, 
				     tensor_template,
				     unit,
				     "FOO_FIELD" );

    TEST_ASSERT( field.getTensorFieldDomain() == domain );
    TEST_ASSERT( field.getTensorFieldEntityType() == iBase_VERTEX );
    TEST_ASSERT( field.getTensorFieldEntityTopology() == iMesh_POINT );
    TEST_ASSERT( field.getTensorFieldCoordType() == FOOD::CARTESIAN );
    TEST_ASSERT( field.getTensorFieldTemplate() == tensor_template );
    TEST_ASSERT( field.getTensorFieldUnit() == unit );
    TEST_ASSERT( field.getTensorFieldName() == "FOO_FIELD" );
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
	new FOOD::Domain(mesh, root_set) );

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
	new FOOD::TensorTemplate(0, 1, FOOD::REAL, quantity) );

    // Create the field.
    FOOD::TensorField<double> field( getDefaultComm<int>(),
				     domain,
				     iBase_VERTEX,
				     iMesh_POINT,
				     FOOD::CARTESIAN, 
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

    TEST_ASSERT( (int) field.getTensorFieldDFMap()->getGlobalNumElements()
		 == num_vertices*mySize );
    TEST_ASSERT( (int) field.getTensorFieldDFView().size() == num_vertices );
    TEST_ASSERT( (int) field.getTensorFieldDFConstView().size() == num_vertices );

    for (int i = 0; i < num_vertices; ++i)
    {
	TEST_ASSERT( (int) field.getTensorFieldDFMap()->getNodeElementList()[i]
		     == myRank*mySize + i );
	TEST_ASSERT( field.getTensorFieldDFView()[i] == (double) i );
	TEST_ASSERT( field.getTensorFieldDFConstView()[i] == (double) i );
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
	new FOOD::Domain(mesh, root_set) );

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
	new FOOD::TensorTemplate(1, 3, FOOD::REAL, quantity) );

    // Create the field.
    FOOD::TensorField<double> field( getDefaultComm<int>(),
				     domain,
				     iBase_REGION,
				     iMesh_TETRAHEDRON,
				     FOOD::CARTESIAN, 
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
    int num_dof = num_tets*tensor_template->getTensorTemplateNumComponents();

    TEST_ASSERT( (int) field.getTensorFieldDFMap()->getGlobalNumElements() ==
		 num_dof*mySize );
    TEST_ASSERT( (int) field.getTensorFieldDFView().size() == num_dof );
    TEST_ASSERT( (int) field.getTensorFieldDFConstView().size() == num_dof );

    int i1 = 0;
    int i2 = 0;
    int i3 = 0;
    for (int i = 0; i < num_tets; ++i)
    {
	i1 = 3*i;
	i2 = 3*i + 1;
	i3 = 3*i + 2;
	TEST_ASSERT( (int) field.getTensorFieldDFMap()->getNodeElementList()[i1]
		     == myRank*mySize + i1 );
	TEST_ASSERT( (int) field.getTensorFieldDFMap()->getNodeElementList()[i2]
		     == myRank*mySize + i2 );
	TEST_ASSERT( (int) field.getTensorFieldDFMap()->getNodeElementList()[i3]
		     == myRank*mySize + i3 );

	TEST_ASSERT( field.getTensorFieldDFView()[i1] == (double) i );
	TEST_ASSERT( field.getTensorFieldDFView()[i2] == (double) i );
	TEST_ASSERT( field.getTensorFieldDFView()[i3] == (double) i );

	TEST_ASSERT( field.getTensorFieldDFConstView()[i1] == (double) i );
	TEST_ASSERT( field.getTensorFieldDFConstView()[i2] == (double) i );
	TEST_ASSERT( field.getTensorFieldDFConstView()[i3] == (double) i );
    }
}

TEUCHOS_UNIT_TEST( TensorField, dof_hex_mesh_region_array_test )
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
	new FOOD::Domain(mesh, root_set) );

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
	new FOOD::TensorTemplate(0, 1, FOOD::REAL, quantity) );

    // Create the field.
    FOOD::TensorField<double> field( getDefaultComm<int>(),
				     domain,
				     iBase_REGION,
				     iMesh_HEXAHEDRON,
				     FOOD::CARTESIAN, 
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

    TEST_ASSERT( (int) field.getTensorFieldDFMap()->getGlobalNumElements()
		 == num_hex*mySize );
    TEST_ASSERT( (int) field.getTensorFieldDFView().size() == num_hex );
    TEST_ASSERT( (int) field.getTensorFieldDFConstView().size() == num_hex );

    for (int i = 0; i < num_hex; ++i)
    {
	TEST_ASSERT( (int) field.getTensorFieldDFMap()->getNodeElementList()[i]
		     == myRank*mySize + i );
	TEST_ASSERT( field.getTensorFieldDFView()[i] == (double) i );
	TEST_ASSERT( field.getTensorFieldDFConstView()[i] == (double) i );
    }

    // Check that the mesh got tagged.
    iBase_TagHandle field_tag = field.getTensorFieldDFTag();

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
    }
}

//---------------------------------------------------------------------------//
//                        end of tstTensorTemplate.cpp
//---------------------------------------------------------------------------//
