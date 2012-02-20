//----------------------------------*-C++-*----------------------------------//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
//
/*!
 * \file   mesh/test/tstDFuncKernel.cpp
 * \author Stuart Slattery
 * \brief  DFuncKernel class unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include <Types.hpp>
#include <DFuncKernel.hpp>

#include <iBase.h>
#include <iMesh.h>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Tuple.hpp>

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

TEUCHOS_UNIT_TEST( DFuncKernel, constructor_test )
{
    FOOD::DFuncKernel<double> hex_fem_div_1_kernel( iBase_REGION,
						    iMesh_HEXAHEDRON,
						    iBase_VERTEX,
						    iMesh_POINT,
						    FOOD::FOOD_CARTESIAN,
						    FOOD::FOOD_FEM,
						    FOOD::FOOD_HDIV,
						    FOOD::FOOD_SHARDSCN,
						    1 );
    TEST_ASSERT( hex_fem_div_1_kernel.getBasisCardinality() == 6 );
    TEST_ASSERT( hex_fem_div_1_kernel.getBasisDegree() == 1 );

    FOOD::DFuncKernel<double> tet_fem_curl_1_kernel( iBase_REGION,
						     iMesh_TETRAHEDRON,
						     iBase_VERTEX,
						     iMesh_POINT,
						     FOOD::FOOD_CARTESIAN,
						     FOOD::FOOD_FEM,
						     FOOD::FOOD_HCURL,
						     FOOD::FOOD_SHARDSCN,
						     1 );
    TEST_ASSERT( tet_fem_curl_1_kernel.getBasisCardinality() == 6 );
    TEST_ASSERT( tet_fem_curl_1_kernel.getBasisDegree() == 1 );

    FOOD::DFuncKernel<double> quad_fem_grad_2_kernel( iBase_FACE,
						      iMesh_QUADRILATERAL,
						      iBase_VERTEX,
						      iMesh_POINT,
						      FOOD::FOOD_CARTESIAN,
						      FOOD::FOOD_FEM,
						      FOOD::FOOD_HGRAD,
						      FOOD::FOOD_SHARDSCN,
						      2 );
    TEST_ASSERT( quad_fem_grad_2_kernel.getBasisCardinality() == 9 );
    TEST_ASSERT( quad_fem_grad_2_kernel.getBasisDegree() == 2 );
}

TEUCHOS_UNIT_TEST( DFuncKernel, hex_evaluation_test )
{
    typedef Intrepid::FieldContainer<double> MDArray;

    FOOD::DFuncKernel<double> hex_kernel( iBase_REGION,
					  iMesh_HEXAHEDRON,
					  iBase_VERTEX,
					  iMesh_POINT,
					  FOOD::FOOD_CARTESIAN,
					  FOOD::FOOD_FEM,
					  FOOD::FOOD_HGRAD,
					  FOOD::FOOD_SHARDSCN,
					  1 );
    
    MDArray coords(1,3);
    coords(0,0) = 0.5;
    coords(0,1) = 0.5;
    coords(0,2) = 0.5;

    MDArray values_at_coords( hex_kernel.getBasisCardinality(),
			      coords.dimension(0) );

    hex_kernel.evaluateValueBasis( values_at_coords, coords);

    TEST_ASSERT( values_at_coords(0,0) == 0.015625 );
    TEST_ASSERT( values_at_coords(1,0) == 0.046875 );
    TEST_ASSERT( values_at_coords(2,0) == 0.140625 );
    TEST_ASSERT( values_at_coords(3,0) == 0.046875 );
    TEST_ASSERT( values_at_coords(4,0) == 0.046875 );
    TEST_ASSERT( values_at_coords(5,0) == 0.140625 );
    TEST_ASSERT( values_at_coords(6,0) == 0.421875 );
    TEST_ASSERT( values_at_coords(7,0) == 0.140625 );
}

//---------------------------------------------------------------------------//
//                        end of tstDFuncKernel.cpp
//---------------------------------------------------------------------------//
