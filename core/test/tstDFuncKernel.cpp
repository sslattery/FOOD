//----------------------------------*-C++-*----------------------------------//
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
    FOOD::DFuncKernel<double> hex_fem_div_1_kernel( iMesh_HEXAHEDRON,
						    FOOD::FOOD_FEM,
						    FOOD::FOOD_DIVERGENCE,
						    1 );
    TEST_ASSERT( hex_fem_div_1_kernel.getBasis()->getCardinality() == 6 );
    TEST_ASSERT( hex_fem_div_1_kernel.getBasis()->getDegree() == 1 );

    FOOD::DFuncKernel<double> tet_fem_curl_1_kernel( iMesh_TETRAHEDRON,
						     FOOD::FOOD_FEM,
						     FOOD::FOOD_CURL,
						     1 );
    TEST_ASSERT( tet_fem_curl_1_kernel.getBasis()->getCardinality() == 6 );
    TEST_ASSERT( tet_fem_curl_1_kernel.getBasis()->getDegree() == 1 );

    FOOD::DFuncKernel<double> quad_fem_grad_2_kernel( iMesh_QUADRILATERAL,
						      FOOD::FOOD_FEM,
						      FOOD::FOOD_GRADIENT,
						      2 );
    TEST_ASSERT( quad_fem_grad_2_kernel.getBasis()->getCardinality() == 9 );
    TEST_ASSERT( quad_fem_grad_2_kernel.getBasis()->getDegree() == 2 );
}

TEUCHOS_UNIT_TEST( DFuncKernel, hex_evaluation_test )
{
    typedef Intrepid::FieldContainer<double> MDArray;

    FOOD::DFuncKernel<double> hex_kernel( iMesh_HEXAHEDRON,
					  FOOD::FOOD_FEM,
					  FOOD::FOOD_GRADIENT,
					  1 );
    
    MDArray coords(1,3);
    coords(0,0) = 0.5;
    coords(0,1) = 0.5;
    coords(0,2) = 0.5;

    MDArray values_at_coords( hex_kernel.getBasis()->getCardinality(),
			      coords.dimension(0) );

    hex_kernel.evaluateDF( values_at_coords, coords);

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
