//----------------------------------*-C++-*----------------------------------//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
//
/*!
 * \file   mesh/test/tstDFuncKernelFactory.cpp
 * \author Stuart Slattery
 * \brief  DFuncKernelFactory class unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include <DiscretizationTypes.hpp>
#include <CellTopologyFactory.hpp>
#include <DFuncKernelFactory.hpp>

#include <iBase.h>
#include <iMesh.h>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Tuple.hpp>

#include <Shards_Array.hpp>

#include <Intrepid_DFuncKernel.hpp>
#include <Intrepid_Cubature.hpp>
#include <Intrepid_DefaultCubatureFactory.hpp>
#include <Intrepid_FieldContainer.hpp>

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

TEUCHOS_UNIT_TEST( DFuncKernelFactory, factory_test )
{
    typedef Intrepid::FieldContainer<double> DoubleContainer;

    // Select a cell topology.
    FOOD::CellTopologyFactory cell_topo_factory;
    Teuchos::RCP<shards::CellTopology> cell_topo = 
	cell_topo_factory.create( iMesh_TETRAHEDRON, 0 );
    int dimension = cell_topo->getDimension();
    int node_count = cell_topo->getNodeCount();
    TEST_ASSERT( dimension == 3 );
    TEST_ASSERT( node_count == 4 );

    // Select a cubature rule.
    Intrepid::DefaultCubatureFactory<double> cubature_factory;
    int cubature_degree = 2;
    Teuchos::RCP< Intrepid::Cubature<double> > cubature = 
	cubature_factory.create(*cell_topo, cubature_degree);
    int num_cube_points = cubature->getNumPoints();
    TEST_ASSERT( num_cube_points == 4 );

    // Select a discrete dfunckernel.
    FOOD::DFuncKernelFactory<double,DoubleContainer> dfunckernel_factory;
    Teuchos::RCP< Intrepid::DFuncKernel<double,DoubleContainer> > dfunckernel = 
	dfunckernel_factory.create( iMesh_TETRAHEDRON,
			      FOOD::FOOD_FEM,
			      FOOD::FOOD_HGRAD,
			      1 );
    int num_fields = dfunckernel->getCardinality();
    TEST_ASSERT( num_fields == 4 );

    // Evaluate the dfunckernel.
}

//---------------------------------------------------------------------------//
//                        end of tstDFuncKernelFactory.cpp
//---------------------------------------------------------------------------//
