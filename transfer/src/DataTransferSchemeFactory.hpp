//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file DataTransferSchemeFactory.hpp
 * \author Stuart Slattery
 * \brief DataTransferSchemeFactory declaration.
 */
//---------------------------------------------------------------------------//

#ifndef FOOD_DATATRANSFERSCHEME_HPP
#define FOOD_DATATRANSFERSCHEME_HPP

#include "TensorField.hpp"
#include "DataTransferScheme.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>

namespace FOOD
{

class DataTransferSchemeFactory
{

  public:

    // Constructor.
    DataTransferSchemeFactory();

    // Destructor.
    ~DataTransferSchemeFactory();

    // Factory method.
    template<class Scalar>
    Teuchos::RCP< DataTransferScheme<Scalar> >
    create( Teuchos::RCP< TensorField<Scalar> > domain,
	    Teuchos::RCP< TensorField<Scalar> > range );
};

} // end namespace FOOD

#endif // end FOOD_DATATRANSFERSCHEMEFACTORY_HPP

//---------------------------------------------------------------------------//
// end DataTransferSchemeFactory.hpp
//---------------------------------------------------------------------------//

