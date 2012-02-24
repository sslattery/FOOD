//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file DataTransferScheme.hpp
 * \author Stuart Slattery
 * \brief DataTransferScheme protocol definition.
 */
//---------------------------------------------------------------------------//

#ifndef FOOD_DATATRANSFERSCHEME_HPP
#define FOOD_DATATRANSFERSCHEME_HPP

#include "TensorField.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>

namespace FOOD
{

template<class Scalar>
class DataTransferScheme
{

  public:

    //@{
    //! Typedefs.
    typedef TensorField<Scalar>                      TensorField_t;
    typedef Teuchos::RCP<TensorField_t>              RCP_TensorField;
    //@}

  public:

    // Constructor.
    DataTransferScheme()
    { /* ... */ }

    // Destructor.
    virtual ~DataTransferScheme()
    { /* ... */ }

    // Setup for interpolation.
    virtual void setup() = 0;

    // Transfer the value of the degrees of freedom from the domain to the range. 
    virtual void transferValueDF() = 0;

    // Get the tensor field corresponding to the domain.
    virtual RCP_TensorField getDomainField() const = 0;

    // Get the tensor field corresponding to the range.
    virtual RCP_TensorField getRangeField() const = 0;
};

} // end namespace FOOD

#endif // end FOOD_DATATRANSFERSCHEME_HPP

//---------------------------------------------------------------------------//
// end DataTransferScheme.hpp
//---------------------------------------------------------------------------//
