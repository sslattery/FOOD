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

    // Create the tensor template for this field.
    Teuchos::RCP<FOOD::TensorTemplate> tensor_template = Teuchos::rcp(
	new FOOD::TensorTemplate(0, 1, FOOD::REAL, quantity) );

    // Create the units for this field.
    Teuchos::RCP<FOOD::Unit> unit = Teuchos::rcp(
	new FOOD::Unit(quantity, 1.4, 4.3, "FOO_UNIT") );

    // Create the field and check basic accessors.
    FOOD::TensorField<double> field( domain, 
				     FOOD::CARTESIAN, 
				     tensor_template,
				     unit,
				     "FOO_FIELD" );

    TEST_ASSERT( field.getTensorFieldDomain() == domain );
    TEST_ASSERT( field.getTensorFieldCoordType() == FOOD::CARTESIAN );
    TEST_ASSERT( field.getTensorFieldTemplate() == tensor_template );
    TEST_ASSERT( field.getTensorFieldUnit() == unit );
    TEST_ASSERT( field.getTensorFieldName() == "FOO_FIELD" );
}

TEUCHOS_UNIT_TEST( TensorField, dof_tag_test )
{
    // Create a single quad for the domain.
    int error;
    iMesh_Instance mesh;
    iMesh_newMesh("", &mesh, &error, 0);
    TEST_ASSERT( iBase_SUCCESS == error );

    iBase_EntityHandle *vertices;
    int vertices_allocated = 0;
    int vertices_size = 0;
    double coords[] = { 0,0,0, 1,0,0, 1,1,0, 1,0,0 };
    iMesh_createVtxArr(mesh,
		       4,
		       iBase_INTERLEAVED,
		       coords,
		       12,
		       &vertices,
		       &vertices_allocated,
		       &vertices_size,
		       &error);
    TEST_ASSERT( iBase_SUCCESS == error );
    
    int quad_creation_status = 0;
    iBase_EntityHandle quad;
    iMesh_CreateEnt(mesh,
		    iMesh_QUADRILATERAL,
		    vertices,
		    4,
		    &quad,
		    &quad_creation_status,
		    &error);
    TEST_ASSERT( iBase_SUCCESS == error );
    TEST_ASSERT( iBase_NEW == quad_creation_status );

    iBaseEntitySetHandle quad_set;
    iMesh_CreateEntSet(mesh, 1, &quad_set, &error);
    TEST_ASSERT( iBase_SUCCESS == error );

    iMesh_addEntToSet(mesh, quad, quad_set, &error);
    TEST_ASSERT( iBase_SUCCESS == error );

    // Create the domain.
    Teuchos::RCP<FOOD::Domain> domain = Teuchos::rcp(
	new FOOD::Domain(mesh, quad_set) );

    // Tag the domain with degrees of freedom.
    iBase_EntityIterator domain_quad_iterator = 0; 
    error = domain.initEntArrIter(iBase_FACE, 
				  iMesh_QUADRILATERAL, 
				  1, 
				  &domain_quad_arr_iterator);
    TEST_ASSERT( iBase_SUCCESS == error );

    iBase_EntityHandle *domain_quads = 0;
    int domain_quads_allocated = 0;
    int domain_quads_size = 0;
    int domain_quads_have_data = 0;
    iMesh_getNextEntIter(mesh,
			 domain_quad_arr_iterator,
			 &domain_quads,
			 &domain_quads_allocated,
			 &domain_quads_size,
			 &domain_quads_have_data,
			 &error);
    TEST_ASSERT( iBase_SUCCESS == error );
    TEST_ASSERT( domain_quads_allocated == 1);
    TEST_ASSERT( domain_quads_size == 1);
    TEST_ASSERT( domain_quads_have_data != 0);

    iBase_TagHandle dof_tag;
    iMesh_createTag(mesh,
		    "DOF_TAG",
		    1,
		    iBase_DOUBLE,
		    &dof_tag,
		    &error,
		    7);
    TEST_ASSERT( iBase_SUCCESS == error );

    std::vector<double> dof_data(1, 2.94857);
    iMesh_setArrData(mesh,
		     domain_quads,
		     domain_quads_size,
		     dof_tag,
		     &dof_data[0],
		     (int) dof_data.size()*sizeof(double),
		     &error);
    TEST_ASSERT( iBase_SUCCESS == error );

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

    // Create the tensor template for this field.
    Teuchos::RCP<FOOD::TensorTemplate> tensor_template = Teuchos::rcp(
	new FOOD::TensorTemplate(0, 1, FOOD::REAL, quantity) );

    // Create the units for this field.
    Teuchos::RCP<FOOD::Unit> unit = Teuchos::rcp(
	new FOOD::Unit(quantity, 1.4, 4.3, "FOO_UNIT") );

    // Create the field and associate with a tag on the domain.
    FOOD::TensorField<double> field( domain, 
				     FOOD::CARTESIAN, 
				     tensor_template,
				     unit,
				     "QUAD_DATA_FIELD" );

    field.attachToTagData( dof_tag, 0.0 );

    // Check the data attachment.
    TEST_ASSERT( field.getTensorFieldDF()->size() == 1 );
    TEST_ASSERT( field.getTensorFieldDF()[0] == 2.94857 );
}

//---------------------------------------------------------------------------//
//                        end of tstTensorTemplate.cpp
//---------------------------------------------------------------------------//
