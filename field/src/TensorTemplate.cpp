//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file TensorTemplate.cpp
 * \author Stuart Slattery
 * \brief Tensor template definition.
 */
//---------------------------------------------------------------------------//

#include "TensorTemplate.hpp"

namespace FOOD
{
/*!
 * \brief Constructor.
 */
TensorTemplate::TensorTemplate(const int order, 
			       const int num_comp, 
			       const int alg_type )
    : d_order(order)
    , d_num_comp(num_comp)
    , d_alg_type(alg_type)
    , d_quantity(0)
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


