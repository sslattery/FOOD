//---------------------------------------------------------------------------//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
//
// \file TensorTemplate.hpp
// \author Stuart Slattery
// \brief Tensor template definition.
//---------------------------------------------------------------------------//

#ifndef FOOD_TENSORTEMPLATE_HPP
#define FOOD_TENSORTEMPLATE_HPP

#include "Types.hpp"
#include "Quantity.hpp"

#include <Teuchos_RCP.hpp>

namespace FOOD
{

class TensorTemplate
{

  public:

    //@{
    //! Typedefs.
    typedef Teuchos::RCP<Quantity>                    RCP_Quantity;
    //@}

  public:

    // Tensor order.
    std::size_t d_order;
    
    // Number of tensor components.
    std::size_t d_num_comp;

    // Enumerated algebraic data type (real, complex, etc.).
    std::size_t d_alg_type;    

    // Physical quantity this tensor represents.
    RCP_Quantity d_quantity;
    
  public:

    // Constructor.
    TensorTemplate( const int order, 
		    const int num_comp, 
		    const int alg_type, 
		    RCP_Quantity quantity );

    // Destructor.
    ~TensorTemplate();

    // Get the order of this tensor template.
    std::size_t getOrder() const
    { return d_order; }

    // Get the number of components in this tensor template.
    std::size_t getNumComponents() const
    { return d_num_comp; }

    // Get the algebraic type of this tensor template.
    std::size_t getAlgType() const
    { return d_alg_type; }

    // Get the physical quantity this tensor template represents.
    RCP_Quantity getQuantity() const
    { return d_quantity; }
};

}

#endif // end FOOD_TENSORTEMPLATE_HPP

//---------------------------------------------------------------------------//
// end TensorTemplate.hpp
//---------------------------------------------------------------------------//

