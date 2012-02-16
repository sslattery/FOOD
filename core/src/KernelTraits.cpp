//---------------------------------------------------------------------------//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
//
// /file KernelTraits.cpp
// /author Stuart Slattery
// /brief Specializations for trait debug strings.
//---------------------------------------------------------------------------//

#include "KernelTraits.hpp"

// Specialize the Evaluation and Data types for the TypeString object in the
// phalanx/src/Phalanx_TypeStrings.hpp file. 

// Scalar Types.
const std::string PHX::TypeString<double>::value = "double";
const std::string PHX::TypeString< Sacado::Fad::DFad<double> >::value = 
  "Sacado::Fad::DFad<double>";

// Evaluation types.
const std::string PHX::TypeString<FOOD::KernelTraits::Residual>::value = 
    "Residual";
const std::string PHX::TypeString<FOOD::KernelTraits::Jacobian>::value = 
    "Jacobian";

//---------------------------------------------------------------------------//
// end KernelTraits.cpp
//---------------------------------------------------------------------------//

