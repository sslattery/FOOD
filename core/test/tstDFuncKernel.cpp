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
						    1,
						    1 );
    TEST_ASSERT( hex_fem_div_1_kernel.getBasis()->getCardinality() == 6 );
    TEST_ASSERT( hex_fem_div_1_kernel.getBasis()->getDegree() == 1 );
    TEST_ASSERT( hex_fem_div_1_kernel.getCubature()->getNumPoints() == 4 );
    TEST_ASSERT( hex_fem_div_1_kernel.getCubature()->getDimension() == 3 );
    TEST_ASSERT( hex_fem_div_1_kernel.getCubature()->getAccuracy() == 3 );

    FOOD::DFuncKernel<double> tet_fem_curl_1_kernel( iMesh_TETRAHEDRON,
						     FOOD::FOOD_FEM,
						     FOOD::FOOD_CURL,
						     1,
						     1 );
    TEST_ASSERT( tet_fem_curl_1_kernel.getBasis()->getCardinality() == 6 );
    TEST_ASSERT( tet_fem_curl_1_kernel.getBasis()->getDegree() == 1 );
    TEST_ASSERT( tet_fem_curl_1_kernel.getCubature()->getNumPoints() == 4 );
    TEST_ASSERT( tet_fem_curl_1_kernel.getCubature()->getDimension() == 3 );
    TEST_ASSERT( tet_fem_curl_1_kernel.getCubature()->getAccuracy() == 3 );

    FOOD::DFuncKernel<double> quad_fem_grad_2_kernel( iMesh_QUADRILATERAL,
						      FOOD::FOOD_FEM,
						      FOOD::FOOD_GRADIENT,
						      2,
						      1);
    TEST_ASSERT( quad_fem_grad_2_kernel.getBasis()->getCardinality() == 9 );
    TEST_ASSERT( quad_fem_grad_2_kernel.getBasis()->getDegree() == 2 );
    TEST_ASSERT( quad_fem_grad_2_kernel.getCubature()->getNumPoints() == 4 );
    TEST_ASSERT( quad_fem_grad_2_kernel.getCubature()->getDimension() == 2 );
    TEST_ASSERT( quad_fem_grad_2_kernel.getCubature()->getAccuracy() == 3 );
}

TEUCHOS_UNIT_TEST( DFuncKernel, evaluation_test )
{

}

//---------------------------------------------------------------------------//
//                        end of tstDFuncKernel.cpp
//---------------------------------------------------------------------------//
