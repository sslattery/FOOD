//---------------------------------------------------------------------------//
// \file Unit.hpp
// \author Stuart Slattery
// \brief Unit definiton
//---------------------------------------------------------------------------//

#ifndef FOOD_UNIT_HPP
#define FOOD_UNIT_HPP

#include <string>

#include "Quantity.hpp"

#include <Teuchos_RCP.hpp>

namespace FOOD
{

class Unit
{
 
  public:

    //@{
    //! Typedefs.
    typedef Teuchos::RCP<Quantity>          RCP_Quantity;
    //@}

  private:

    // Quantity this unit is a measure of.
    RCP_Quantity d_quantity;

    // Scale of this unit relative to the default unit of measure for the
    // associated quantity.
    double d_scale;

    // Offset of this unit relative to the default unit of measure for the
    // associated quantity.
    double d_offset;

    // Name of this unit.
    std::string d_name;

  public:

    // Constructor.
    Unit( RCP_Quantity quantity, 
	  double scale, 
	  double offset,
	  const std::string &name);

    // Destructor.
    ~Unit();

    //! Get the quantity associated with this unit.
    RCP_Quantity getQuantity() const
    { return d_quantity; }

    //! Get the scale of this unit with respect to the associated quantity.
    double getScale() const
    { return d_scale; }

    //! Get the offset of this unit with respect to the associated quantity.
    double getOffset() const
    { return d_offset; }

    //! Get the name of this unit.
    const std::string& getName() const
    { return d_name; }

    // Set the name of this unit.
    void setName( const std::string &new_name );

    // Set the name of the quantity this unit is a measure of.
    void setQuantityName( const std::string &new_name );
};

} // end namespace FOOD

#endif // end FOOD_UNIT_HPP

//---------------------------------------------------------------------------//
// end Unit.hpp
//---------------------------------------------------------------------------//

