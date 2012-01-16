//----------------------------------*-C++-*----------------------------------//
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

    TEST_ASSERT( unit.getUnitQuantity() == quantity );
    TEST_ASSERT( unit.getUnitScale() == 1.4 );
    TEST_ASSERT( unit.getUnitName() == "FOO_UNIT" );
   
    unit.setUnitName( "BAR_UNIT" );
    TEST_ASSERT( unit.getUnitName() == "BAR_UNIT" );
    
    unit.setUnitQuantityName( "BAR_QUANTITY" );
    TEST_ASSERT( unit.getUnitQuantity()->getQuantityName() == "BAR_QUANTITY" );
}

//---------------------------------------------------------------------------//
//                        end of tstUnit.cpp
//---------------------------------------------------------------------------//
