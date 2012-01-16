//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/test/tstTensorTemplate.cpp
 * \author Stuart Slattery
 * \brief  TensorTemplate class unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include <Types.hpp>
#include <Quantity.hpp>
#include <TensorTemplate.hpp>

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

TEUCHOS_UNIT_TEST( TensorTemplate, constructor_test )
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

    FOOD::TensorTemplate tensor_template(1, 3, FOOD::REAL, quantity);

    TEST_ASSERT( tensor_template.getTensorTemplateOrder() == 1 );
    TEST_ASSERT( tensor_template.getTensorTemplateNumComponents() == 3 );
    TEST_ASSERT( tensor_template.getTensorTemplateAlgType() == FOOD::REAL );
    TEST_ASSERT( tensor_template.getTensorTemplateQuantity() == quantity );
}

//---------------------------------------------------------------------------//
//                        end of tstTensorTemplate.cpp
//---------------------------------------------------------------------------//
