//---------------------------------------------------------------------------//
// \file TensorTemplate.cpp
// \author Stuart Slattery
// \brief Tensor template definition.
//---------------------------------------------------------------------------//

#include "TensorTemplate.hpp"

namespace FOOD
{
/*!
 * \brief Constructor.
 */
TensorTemplate::TensorTemplate(int order, int num_comp, 
			       int alg_type, RCP_Quantity quantity )
    : d_order(order)
    , d_num_comp(num_comp)
    , d_alg_type(alg_type)
    , d_quantity(quantity)
{ /* ... */ }

/*!
 * \brief Destructor.
 */
TensorTemplate::~TensorTemplate()
{ /* ... */ }

}

//---------------------------------------------------------------------------//
// end TensorTemplate.cpp
//---------------------------------------------------------------------------//


