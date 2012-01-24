//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/test/tstQuantity.cpp
 * \author Stuart Slattery
 * \brief  Quantity class unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include <Quantity.hpp>

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

TEUCHOS_UNIT_TEST( Quantity, constructor_test )
{
    Teuchos::Tuple<int,7> numerator;
    Teuchos::Tuple<int,7> denominator;
    for (int i = 0; i < 7; ++i)
    {
	numerator[i] = i;
	denominator[i] = 6 - i;
    }

    FOOD::Quantity quantity(numerator, denominator, "FOO_QUANTITY");

    for (int i = 0; i < 7; ++i)
    {
	TEST_ASSERT( quantity.getNumerator()[i] == numerator[i] );
	TEST_ASSERT( quantity.getDenominator()[i] == denominator[i] );
    }
    TEST_ASSERT( quantity.getName() == "FOO_QUANTITY" );
}

TEUCHOS_UNIT_TEST( Quantity, derived_quantity_test )
{
    Teuchos::Tuple<int,7> first_numerator;
    Teuchos::Tuple<int,7> first_denominator;
    Teuchos::Tuple<int,7> second_numerator;
    Teuchos::Tuple<int,7> second_denominator;
    for (int i = 0; i < 7; ++i)
    {
	first_numerator[i] = i;
	first_denominator[i] = 6 - i;
	second_numerator[i] = i;
	second_denominator[i] = 2*i+1;
    }

    FOOD::Quantity first_quantity(first_numerator, 
				  first_denominator, 
				  "FIRST_QUANTITY");
    FOOD::Quantity second_quantity(second_numerator, 
				   second_denominator, 
				   "SECOND_QUANTITY");

    FOOD::Quantity derived_quantity = first_quantity*second_quantity;

    TEST_ASSERT( derived_quantity.getNumerator()[0]   == 0 );
    TEST_ASSERT( derived_quantity.getDenominator()[0] == 7 );
    TEST_ASSERT( derived_quantity.getNumerator()[1]   == 0 );
    TEST_ASSERT( derived_quantity.getDenominator()[1] == 6 );
    TEST_ASSERT( derived_quantity.getNumerator()[2]   == 0 );
    TEST_ASSERT( derived_quantity.getDenominator()[2] == 5 );
    TEST_ASSERT( derived_quantity.getNumerator()[3]   == 0 );
    TEST_ASSERT( derived_quantity.getDenominator()[3] == 4 );
    TEST_ASSERT( derived_quantity.getNumerator()[4]   == 0 );
    TEST_ASSERT( derived_quantity.getDenominator()[4] == 3 );
    TEST_ASSERT( derived_quantity.getNumerator()[5]   == 0 );
    TEST_ASSERT( derived_quantity.getDenominator()[5] == 2 );
    TEST_ASSERT( derived_quantity.getNumerator()[6]   == 0 );
    TEST_ASSERT( derived_quantity.getDenominator()[6] == 1 );
    TEST_ASSERT( derived_quantity.getName() == "DERIVED" );

    derived_quantity.setName("THIRD_QUANTITY");
    TEST_ASSERT( derived_quantity.getName() == "THIRD_QUANTITY" );
}

//---------------------------------------------------------------------------//
//                        end of tstQuantity.cpp
//---------------------------------------------------------------------------//
