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
    iMesh_createTag(mesh,
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
    iMesh_createVtxArr(mesh,
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

    iMesh_setDblArrData(mesh,
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
    iMesh_createTag(mesh,
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
		iMesh_createEnt(mesh,
				iMesh_HEXAHEDRON,
				connectivity,
				8,
				&hexahedron,
				&hex_creation_status,
				&error);
		assert( iBase_SUCCESS == error );
		assert( iBase_NEW == hex_creation_status );

		iMesh_setDblData(mesh,
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

    // Create the tensor template for this field.
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
    field.attachToTagData( vertex_tag );

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

//---------------------------------------------------------------------------//
//                        end of tstTensorTemplate.cpp
//---------------------------------------------------------------------------//
