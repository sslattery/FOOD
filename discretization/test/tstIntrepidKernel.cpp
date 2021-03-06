//----------------------------------*-C++-*----------------------------------//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
//
/*!
 * \file   mesh/test/tstIntrepidKernel.cpp
 * \author Stuart Slattery
 * \brief  IntrepidKernel class unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include <DiscretizationTypes.hpp>
#include <DFuncKernel.hpp>
#include <IntrepidKernel.hpp>

#include <iBase.h>
#include <iMesh.h>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Tuple.hpp>

#include <Intrepid_FieldContainer.hpp>
#include <Intrepid_HGRAD_HEX_C1_FEM.hpp>

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

TEUCHOS_UNIT_TEST( IntrepidKernel, intrepid_kernel_test)
{
    typedef Intrepid::FieldContainer<double> MDArray;
    Teuchos::RCP< FOOD::DFuncKernel<double> > kernel = 
	Teuchos::rcp( 
	    new FOOD::IntrepidKernel<double>( 
		Teuchos::rcp( new Intrepid::Basis_HGRAD_HEX_C1_FEM<double,MDArray>() ),
		iBase_REGION,
		iMesh_HEXAHEDRON,
		FOOD::FOOD_FEM,
		FOOD::FOOD_HGRAD ) );

    TEST_ASSERT( kernel->getCardinality()        == 8 );
    TEST_ASSERT( kernel->getDegree()             == 1 );
    TEST_ASSERT( kernel->getEntityType()         == iBase_REGION );
    TEST_ASSERT( kernel->getEntityTopology()     == iMesh_HEXAHEDRON );
    TEST_ASSERT( kernel->getCNType()             == FOOD::FOOD_SHARDSCN );
    TEST_ASSERT( kernel->getCoordType()          == FOOD::FOOD_CARTESIAN );
    TEST_ASSERT( kernel->getDiscretizationType() == FOOD::FOOD_FEM );
    TEST_ASSERT( kernel->getFunctionSpaceType()  == FOOD::FOOD_HGRAD );
}

//---------------------------------------------------------------------------//
//                        end of tstIntrepidKernel.cpp
//---------------------------------------------------------------------------//
