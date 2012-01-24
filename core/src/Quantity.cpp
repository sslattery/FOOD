//---------------------------------------------------------------------------//
// \file Quantity.cpp
// \author Stuart Slattery
// \brief Quantity defintion.
//---------------------------------------------------------------------------//

#include "Quantity.hpp"

namespace FOOD
{

/*!
 * \brief Constructor.
 * 
 * The seven basic quantities are defined as SI base units:
 *
 * http://physics.nist.gov/cuu/Units/units.html
 *
 * and are expected in the following order:
 *
 * [ length                    : meter ,
 *   mass                      : kilogram ,
 *   time                      : second ,
 *   electric current          : ampere ,
 *   thermodynamic temperature : kelvin ,
 *   amount of a substance     : mole ,
 *   luminous intensity        : candela ]
 */
Quantity::Quantity( const Teuchos::Tuple<int,7> &numerator,
		    const Teuchos::Tuple<int,7> &denominator,
		    const std::string &name )
    : d_numerator(numerator)
    , d_denominator(denominator)
    , d_name(name)
{ /* ... */ }

/*!
 * \brief Destructor.
 */
Quantity::~Quantity()
{ /* ... */ }

/*!
 * \brief Multiply operator.
 *
 * Derive a new quantity from two others.
 */
Quantity Quantity::operator*(const Quantity &quantity)
{
    int product = 0;
    Teuchos::Tuple<int,7> derived_numerator;
    Teuchos::Tuple<int,7> derived_denominator;

    for (int i = 0; i < 7; ++i)
    {
	product = quantity.getNumerator()[i] + d_numerator[i] -
		  quantity.getDenominator()[i] - d_denominator[i];

	if ( product > 0 )
	{
	    derived_numerator[i] = product;
	    derived_denominator[i] = 0;
	}
	else if ( product < 0 )
	{
	    derived_numerator[i] = 0;
	    derived_denominator[i] = -1*product;
	}
	else
	{
	    derived_numerator[i] = 0;
	    derived_denominator[i] = 0;
	}
    }

    Quantity derived_quantity(derived_numerator, 
			      derived_denominator, 
			      "DERIVED");
    return derived_quantity;
}

/*!
 * \brief Rename a quantity. Do this if derived.
 */
void Quantity::setName( const std::string &new_name )
{
    d_name = new_name;
}

} // end namepace FOOD

//---------------------------------------------------------------------------//
// end Quantity.cpp
//---------------------------------------------------------------------------//

