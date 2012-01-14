//---------------------------------------------------------------------------//
// \file Unit.cpp
// \author Stuart Slattery
// \brief Unit definiton
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
 * \brief Multiplication operator.
 */
Unit Unit::operator*( const Unit &unit )
{
    
}

/*!
 * \brief Set the name of this unit.
 */
void Unit::setUnitName( const std::string &new_name )
{
    d_name = new_name;
}

/*!
 * \brief Set the name of the quantity this unit is a measure of.
 */
void Unit::setUnitQuantityName( const std::string &new_name )
{
    d_quantity->setQuantityName(new_name);
}

} // end namespace FOOD

//---------------------------------------------------------------------------//
// end Unit.cpp
//---------------------------------------------------------------------------//
