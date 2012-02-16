//----------------------------------*-C++-*----------------------------------//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
//
/*!
 * \file   mesh/test/tstUnit.cpp
 * \author Stuart Slattery
 * \brief  Unit class unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include <Quantity.hpp>
#include <Unit.hpp>

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

TEUCHOS_UNIT_TEST( Unit, constructor_test )
{
    Teuchos::Tuple<int,7> numerator;
    Teuchos::Tuple<int,7> denominator;
    for (int i = 0; i < 7; ++i)
    {
	numerator[i] = i;
	denominator[i] = 6 - i;
    }

    Teuchos::RCP<FOOD::Quantity> quantity = Teuchos::rcp(
	new FOOD::Quantity(numerator, denominator, "FOO_QUANTITY") );

    FOOD::Unit unit(quantity, 1.4, 4.3, "FOO_UNIT");

    TEST_ASSERT( unit.getQuantity() == quantity );
    TEST_ASSERT( unit.getScale() == 1.4 );
    TEST_ASSERT( unit.getName() == "FOO_UNIT" );
   
    unit.setName( "BAR_UNIT" );
    TEST_ASSERT( unit.getName() == "BAR_UNIT" );
    
    unit.setQuantityName( "BAR_QUANTITY" );
    TEST_ASSERT( unit.getQuantity()->getName() == "BAR_QUANTITY" );
}

//---------------------------------------------------------------------------//
//                        end of tstUnit.cpp
//---------------------------------------------------------------------------//
