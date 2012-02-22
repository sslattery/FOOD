//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file DFuncKernelFactory.hpp
 * \author Stuart Slattery
 * \brief Factory method declaration for basis generation.
 */
//---------------------------------------------------------------------------//

#ifndef FOOD_DFUNCKERNELFACTORY_HPP
#define FOOD_DFUNCKERNELFACTORY_HPP

#include "DFuncKernel.hpp"

#include <Teuchos_RCP.hpp>

namespace FOOD
{

template<class Scalar>
class DFuncKernelFactory
{

  public:
    
    // Constructor.
    DFuncKernelFactory();

    // Destructor.
    ~DFuncKernelFactory();

    // Factory method.
    Teuchos::RCP< DFuncKernel<Scalar> > 
    create( const int entity_topology,
	    const int discretization_type,
	    const int function_space_type,
	    const int degree );
};

} // end namespace FOOD

#include "DFuncKernelFactory_Def.hpp"

#endif // end FOOD_DFUNCKERNELFACTORY_HPP

//---------------------------------------------------------------------------//
// end DFuncKernelFactory.hpp
//---------------------------------------------------------------------------//
