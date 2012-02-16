//---------------------------------------------------------------------------//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
//
// \file Quantity.hpp
// \author Stuart Slattery
// \brief Quantity defintion.
//---------------------------------------------------------------------------//

#ifndef FOOD_QUANTITY_HPP
#define FOOD_QUANTITY_HPP

#include <string>

#include <Teuchos_Tuple.hpp>
#include <Teuchos_ArrayView.hpp>

namespace FOOD
{

class Quantity
{
  private:

    // Powers of the 7 basic SI quantities in the numerator.
    Teuchos::Tuple<int,7> d_numerator;

    // Powers of the 7 basic SI quantities in the denominator.
    Teuchos::Tuple<int,7> d_denominator;

    // Name of this quantity.
    std::string d_name;

  public:

    // Constructor.
    Quantity( const Teuchos::Tuple<int,7> &numerator,
	      const Teuchos::Tuple<int,7> &denominator,
	      const std::string &name );

    // Destructor.
    ~Quantity();

    //! Get the numerator powers for this quantity.
    const Teuchos::Tuple<int,7>& getNumerator() const
    { return d_numerator; }

    //! Get the denominator powers for this quantity.
    const Teuchos::Tuple<int,7>& getDenominator() const
    { return d_denominator; }

    //! Get the name for this quantity.
    const std::string& getName() const
    { return d_name; }

    // Multiplication operator. Used for creating derived quantities.
    Quantity operator*( const Quantity &quantity );

    // Rename a quantity. Do this if derived.
    void setName( const std::string &new_name );
};

} // end namespace FOOD

#endif // end FOOD_QUANTITY_HPP

//---------------------------------------------------------------------------//
// end Quantity.hpp
//---------------------------------------------------------------------------//


