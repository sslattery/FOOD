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
    int d_order;
    
    // Number of tensor components.
    int d_num_comp;

    // Enumerated algebraic data type (real, complex, etc.).
    int d_alg_type;    

    // Physical quantity this tensor represents.
    RCP_Quantity d_quantity;
    
  public:

    // Constructor.
    TensorTemplate( int order, int num_comp, 
		    int alg_type, RCP_Quantity quantity );

    // Destructor.
    ~TensorTemplate();

    // Get the order of this tensor template.
    int getTensorTemplateOrder() const
    { return d_order; }

    // Get the number of components in this tensor template.
    int getTensorTemplateNumComponents() const
    { return d_num_comp; }

    // Get the algebraic type of this tensor template.
    int getTensorTemplateAlgType() const
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

