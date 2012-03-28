//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file   Exception.hpp
 * \author Stuart Slattery
 * \brief  FOOD exception handling and error policy declaration.
 */
//---------------------------------------------------------------------------//

#ifndef FOOD_EXCEPTION_HPP
#define FOOD_EXCEPTION_HPP

#include <stdexcept>
#include <string>

namespace FOOD
{

/*!
 * \brief Exception class to be thrown when function preconditions are not
 * met.
 */
class PreconditionException : public std::runtime_error
{
  public:
    PreconditionException( const std::string &msg )
	: std::runtime_error( msg )
    { /* ... */ }
};

/*!
 * \brief Exception class to be thrown when function postconditions are not
 * met. 
 */
class PostconditionException : public std::runtime_error
{
  public:
    PostconditionException( const std::string &msg )
	: std::runtime_error( msg )
    { /* ... */ }
};

/*!
 * \brief Exception class to be thrown when a function alters an invariant.
 */
class InvariantException : public std::runtime_error
{
  public:
    InvariantException( const std::string &msg )
	: std::runtime_error( msg )
    { /* ... */ }
};

// Test for a precondition exception.
void testPrecondition( bool throw_if_false, const std::string &msg );

// Test for a postcondition exception.
void testPostcondition( bool throw_if_false, const std::string &msg );

// Test for a Invariant exception.
void testInvariant( bool throw_if_false, const std::string &msg );

} // end namespace FOOD

#endif // end FOOD_EXCEPTION_HPP

//---------------------------------------------------------------------------//
// end Exception.hpp
//---------------------------------------------------------------------------//
