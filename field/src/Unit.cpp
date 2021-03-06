//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file Unit.cpp
 * \author Stuart Slattery
 * \brief Unit definiton
 */
//---------------------------------------------------------------------------//

#include "Unit.hpp"

namespace FOOD
{

/*!
 * \brief Constructor.
 */
Unit::Unit( RCP_Quantity quantity, 
	    double scale, 
	    double offset,
	    const std::string &name)
    : d_quantity(quantity)
    , d_scale(scale)
    , d_offset(offset)
    , d_name(name)
{ /* ... */ }

/*!
 * \brief Destructor.
 */
Unit::~Unit()
{ /* ... */ }

/*!
 * \brief Set the name of this unit.
 */
void Unit::setName( const std::string &new_name )
{
    d_name = new_name;
}

/*!
 * \brief Set the name of the quantity this unit is a measure of.
 */
void Unit::setQuantityName( const std::string &new_name )
{
    d_quantity->setName(new_name);
}

} // end namespace FOOD

//---------------------------------------------------------------------------//
// end Unit.cpp
//---------------------------------------------------------------------------//
