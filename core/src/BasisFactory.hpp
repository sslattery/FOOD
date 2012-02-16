//---------------------------------------------------------------------------//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
//
// \file BasisFactory.hpp
// \author Stuart Slattery
// \brief Factory method declaration for basis generation.
//---------------------------------------------------------------------------//

#ifndef FOOD_BASISFACTORY_HPP
#define FOOD_BASISFACTORY_HPP

#include "Types.hpp"

#include <iBase.h>
#include <iMesh.h>

#include <Teuchos_RCP.hpp>

#include <Shards_CellTopology.hpp>

#include <Intrepid_Basis.hpp>
#include <Intrepid_HCURL_HEX_I1_FEM.hpp>
#include <Intrepid_HCURL_QUAD_I1_FEM.hpp>
#include <Intrepid_HCURL_TET_I1_FEM.hpp>
#include <Intrepid_HCURL_TRI_I1_FEM.hpp>
#include <Intrepid_HDIV_HEX_I1_FEM.hpp>
#include <Intrepid_HDIV_QUAD_I1_FEM.hpp>
#include <Intrepid_HDIV_TET_I1_FEM.hpp>
#include <Intrepid_HDIV_TRI_I1_FEM.hpp>
#include <Intrepid_HGRAD_HEX_C1_FEM.hpp>
#include <Intrepid_HGRAD_HEX_C2_FEM.hpp>
#include <Intrepid_HGRAD_LINE_C1_FEM.hpp>
#include <Intrepid_HGRAD_QUAD_C1_FEM.hpp>
#include <Intrepid_HGRAD_QUAD_C2_FEM.hpp>
#include <Intrepid_HGRAD_TET_C1_FEM.hpp>
#include <Intrepid_HGRAD_TET_C2_FEM.hpp>
#include <Intrepid_HGRAD_TRI_C1_FEM.hpp>
#include <Intrepid_HGRAD_TRI_C2_FEM.hpp>
#include <Intrepid_HGRAD_WEDGE_C1_FEM.hpp>
#include <Intrepid_HGRAD_WEDGE_C2_FEM.hpp>

namespace FOOD
{

template<class Scalar, class ArrayScalar>
class BasisFactory
{

  public:
    
    // Constructor.
    BasisFactory();

    // Destructor.
    ~BasisFactory();

    // Factory method.
    Teuchos::RCP< Intrepid::Basis<Scalar,ArrayScalar> > 
    create( const int entity_topology,
	    const int discretization_type,
	    const int basis_function_space,
	    const int basis_degree );
};

} // end namespace FOOD

#include "BasisFactory_Def.hpp"

#endif // end FOOD_BASISFACTORY_HPP

//---------------------------------------------------------------------------//
// end BasisFactory.hpp
//---------------------------------------------------------------------------//
