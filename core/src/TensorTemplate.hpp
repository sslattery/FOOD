//---------------------------------------------------------------------------//
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
    std::size_t getTensorTemplateOrder() const
    { return d_order; }

    // Get the number of components in this tensor template.
    std::size_t getTensorTemplateNumComponents() const
    { return d_num_comp; }

    // Get the algebraic type of this tensor template.
    std::size_t getTensorTemplateAlgType() const
    { return d_alg_type; }

    // Get the physical quantity this tensor template represents.
    RCP_Quantity getTensorTemplateQuantity() const
    { return d_quantity; }
};

}

#endif // end FOOD_TENSORTEMPLATE_HPP

//---------------------------------------------------------------------------//
// end TensorTemplate.hpp
//---------------------------------------------------------------------------//

