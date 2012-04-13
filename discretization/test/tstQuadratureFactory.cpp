//----------------------------------*-C++-*----------------------------------//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
//
/*!
 * \file   mesh/test/tstQuadratureFactory.cpp
 * \author Stuart Slattery
 * \brief  QuadratureFactory class unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include <QuadratureFactory.hpp>
#include <Quadrature.hpp>

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

TEUCHOS_UNIT_TEST( QuadratureFactory, factory_test )
{
    FOOD::QuadratureFactory<double> quadrature_factory;
    Teuchos::RCP< FOOD::Quadrature<double> > quadrature = 
	quadrature_factory.create( iBase_REGION,
				   iMesh_TETRAHEDRON,
				   10,
				   2 );

    TEST_ASSERT( quadrature->getNumPoints() == 4 );
    TEST_ASSERT( quadrature->getDimension() == 3 );
    TEST_ASSERT( quadrature->getEntityType() == iBase_REGION );
    TEST_ASSERT( quadrature->getEntityTopology() == iMesh_TETRAHEDRON );
}

//---------------------------------------------------------------------------//
//                        end of tstQuadratureFactory.cpp
//---------------------------------------------------------------------------//
