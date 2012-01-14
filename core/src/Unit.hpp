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

namspace FOOD
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
    RCP_Quantity getUnitQuantity() const
    { return d_quantity; }

    //! Get the scale of this unit with respect to the associated quantity.
    double getUnitScale() const
    { return d_scale; }

    //! Get the offset of this unit with respect to the associated quantity.
    double getUnitOffset() const
    { return d_offset; }

    //! Get the name of this unit.
    const std::string& getUnitName() const
    { return d_name; }

    // Multiplication operator.
    Unit operator*( const Unit &unit );

    // Set the name of this unit.
    void setUnitName( const std::string &new_name );

    // Set the name of the quantity this unit is a measure of.
    void setUnitQuantityName( const std::string &new_name )
};

} // end namespace FOOD

#endif // end FOOD_UNIT_HPP

//---------------------------------------------------------------------------//
// end Unit.hpp
//---------------------------------------------------------------------------//

