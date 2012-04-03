//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file QuadratureFactory.hpp
 * \author Stuart Slattery
 * \brief Factory method declaration for quadrature rule generation.
 */
//---------------------------------------------------------------------------//

#ifndef FOOD_QUADRATUREFACTORY_HPP
#define FOOD_QUADRATUREFACTORY_HPP

#include "Quadrature.hpp"

#include <Teuchos_RCP.hpp>

namespace FOOD
{

template<class Scalar>
class QuadratureFactory
{

  public:
    
    // Constructor.
    QuadratureFactory();

    // Destructor.
    ~QuadratureFactory();

    // Factory method.
    Teuchos::RCP< Quadrature<Scalar> > 
    create( const int entity_type,
	    const int entity_topology,
	    const int degree );
};

} // end namespace FOOD

#include "QuadratureFactory_Def.hpp"

#endif // end FOOD_QUADRATUREFACTORY_HPP

//---------------------------------------------------------------------------//
// end QuadratureFactory.hpp
//---------------------------------------------------------------------------//
