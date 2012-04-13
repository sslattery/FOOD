//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file   Exception.cpp
 * \author Stuart Slattery
 * \brief  FOOD exception handling and error policy definition.
 */
//---------------------------------------------------------------------------//

#include <iostream>

#include "Exception.hpp"

#include <Teuchos_TestForException.hpp>

namespace FOOD
{
/*!
 * \brief Test for a precondition exception.
 */
void testPrecondition( bool throw_if_false, const std::string &msg )
{
    TEUCHOS_TEST_FOR_EXCEPTION( !throw_if_false, 
				PreconditionException,
				msg << std::endl );
}

/*!
 * \brief Test for a postcondition exception.
 */
void testPostcondition( bool throw_if_false, const std::string &msg )
{
    TEUCHOS_TEST_FOR_EXCEPTION( !throw_if_false, 
				PostconditionException,
				msg << std::endl );
}

/*!
 * \brief Test for a Invariant exception.
 */
void testInvariant( bool throw_if_false, const std::string &msg )
{
    TEUCHOS_TEST_FOR_EXCEPTION( !throw_if_false, 
				InvariantException,
				msg << std::endl );
}

} // end namespace FOOD

//---------------------------------------------------------------------------//
// end Exception.cpp
//---------------------------------------------------------------------------//
