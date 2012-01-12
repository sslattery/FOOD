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
		    const std::string &name,
		    const std::string &abbr)
    : d_numerator(numerator)
    , d_denominator(denominator)
    , d_name(name)
    , d_abbr(abbr)
{ /* ... */ }

/*!
 * \brief Destructor.
 */
Quantity::~Quantity()
{ /* ... */ }

} // end namepace FOOD

//---------------------------------------------------------------------------//
// end Quantity.cpp
//---------------------------------------------------------------------------//

