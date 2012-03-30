//----------------------------------*-C++-*----------------------------------//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
//
/*!
 * \file   mesh/test/tstIntrepidQuadrature.cpp
 * \author Stuart Slattery
 * \brief  IntrepidQuadrature class unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "TopologyTools.hpp"
#include "CellTopologyFactory.hpp"
#include "Quadrature.hpp"
#include "IntrepidQuadrature.hpp"

#include <iBase.h>
#include <iMesh.h>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Tuple.hpp>

#include <Intrepid_DefaultCubatureFactory.hpp>

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//

template<class Ordinal>
Teuchos::RCP<const Teuchos::Comm<Ordinal> > getDefaultComm()
{
#ifdef HAVE_MPI
    return Teuchos::DefaultComm<Ordinal>::getComm();
#else
    return Teuchos::rcp(new Teuchos::SerialComm<Ordinal>() );
#endif
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEUCHOS_UNIT_TEST( IntrepidQuadrature, intrepid_quadrature_test)
{
    FOOD::CellTopologyFactory cell_topo_factory;
    Teuchos::RCP<shards::CellTopology> shards_topo = 
	cell_topo_factory.create( 
	    iMesh_TETRAHEDRON,
	    FOOD::TopologyTools::numLinearNodes( iMesh_TETRAHEDRON ) );

    Intrepid::DefaultCubatureFactory<double> cub_factory;
    Teuchos::RCP< Intrepid::Cubature<double> > intrepid_cub =
	cub_factory.create( *shards_topo, 2 );

    Teuchos::RCP< FOOD::Quadrature<double> > quadrature = Teuchos::rcp(
	new FOOD::IntrepidQuadrature<double>( intrepid_cub, 
					      iBase_REGION, 
					      iMesh_TETRAHEDRON ) );

    TEST_ASSERT( quadrature->getNumPoints() == 4 );
    TEST_ASSERT( quadrature->getDimension() == 3 );
    TEST_ASSERT( quadrature->getEntityType() == iBase_REGION );
    TEST_ASSERT( quadrature->getEntityTopology() == iMesh_TETRAHEDRON );
}

//---------------------------------------------------------------------------//
//                        end of tstIntrepidQuadrature.cpp
//---------------------------------------------------------------------------//
