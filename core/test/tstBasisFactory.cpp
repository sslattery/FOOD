//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/test/tstBasisFactory.cpp
 * \author Stuart Slattery
 * \brief  BasisFactory class unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include <Types.hpp>
#include <CellTopologyFactory.hpp>
#include <BasisFactory.hpp>

#include <iBase.h>
#include <iMesh.h>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Tuple.hpp>

#include <Shards_Array.hpp>

#include <Intrepid_Basis.hpp>
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

TEUCHOS_UNIT_TEST( BasisFactory, factory_test )
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

    // Select a discrete basis.
    FOOD::BasisFactory<double,DoubleContainer> basis_factory;
    Teuchos::RCP< Intrepid::Basis<double,DoubleContainer> > basis = 
	basis_factory.create( iMesh_TETRAHEDRON,
			      FOOD::FOOD_FEM,
			      FOOD::FOOD_HGRAD,
			      1 );
    int num_fields = basis->getCardinality();
    TEST_ASSERT( num_fields == 4 );

    // Evaluate the basis.
}

//---------------------------------------------------------------------------//
//                        end of tstBasisFactory.cpp
//---------------------------------------------------------------------------//
