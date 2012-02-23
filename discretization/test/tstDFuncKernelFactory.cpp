//----------------------------------*-C++-*----------------------------------//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
//
/*!
 * \file   mesh/test/tstDFuncKernelFactory.cpp
 * \author Stuart Slattery
 * \brief  DFuncKernelFactory class unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include <DiscretizationTypes.hpp>
#include <CellTopologyFactory.hpp>
#include <DFuncKernelFactory.hpp>
#include <DFuncKernel.hpp>

#include <iBase.h>
#include <iMesh.h>

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

TEUCHOS_UNIT_TEST( DFuncKernelFactory, factory_test )
{
    FOOD::DFuncKernelFactory<double> kernel_factory;
    Teuchos::RCP< FOOD::DFuncKernel<double> > kernel = 
	kernel_factory.create( iMesh_HEXAHEDRON,
			       FOOD::FOOD_FEM,
			       FOOD::FOOD_HGRAD,
			       1 );

    TEST_ASSERT( kernel->getCardinality()        == 8 );
    TEST_ASSERT( kernel->getDegree()             == 1 );
    TEST_ASSERT( kernel->getTopology()           == iMesh_HEXAHEDRON );
    TEST_ASSERT( kernel->getCNType()             == FOOD::FOOD_SHARDSCN );
    TEST_ASSERT( kernel->getCoordType()          == FOOD::FOOD_CARTESIAN );
    TEST_ASSERT( kernel->getDiscretizationType() == FOOD::FOOD_FEM );
    TEST_ASSERT( kernel->getFunctionSpaceType()  == FOOD::FOOD_HGRAD );
}

//---------------------------------------------------------------------------//
//                        end of tstDFuncKernelFactory.cpp
//---------------------------------------------------------------------------//
