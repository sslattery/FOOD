//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file DataTransferSchemeFactory.cpp
 * \author Stuart Slattery
 * \brief DataTransferSchemeFactory definition.
 */
//---------------------------------------------------------------------------//

#include "DataTransferSchemeFactory.hpp"

namespace FOOD
{
/*!
 * \brief Constructor.
 */
DataTransferSchemeFactory::DataTransferSchemeFactory()
{ /* ... */ }

/*!
 * \brief Destructor.
 */
DataTransferSchemeFactory::~DataTransferSchemeFactory()
{ /* ... */ }

/*!
 * \brief Factory method.
 */
template<class Scalar>
Teuchos::RCP< DataTransferScheme<Scalar> >
DataTransferSchemeFactory::create( Teuchos::RCP< TensorField<Scalar> > domain,
				   Teuchos::RCP< TensorField<Scalar> > range )
{

}

} // end namespace FOOD

//---------------------------------------------------------------------------//
// end DataTransferSchemeFactory.cpp
//---------------------------------------------------------------------------//

