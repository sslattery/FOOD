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

TEUCHOS_UNIT_TEST( IntrepidQuadrature, integration_test )
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

    // Create a tetrahedron.
    int error;
    iMesh_Instance mesh;
    iMesh_newMesh( "", &mesh, &error, 0);
    TEST_ASSERT( iBase_SUCCESS == error );

    double vtx_coords[12] = { 0,0,0, 1,0,0, 1,1,0, 0,1,0 };
    int num_verts = 4;
    int new_coords_size = 3*num_verts;
    int new_vertex_handles_allocated = 0;
    int new_vertex_handles_size = 0;
    iBase_EntityHandle *vertex_handles = 0;
    iMesh_createVtxArr( mesh,
			num_verts,
			iBase_INTERLEAVED,
			vtx_coords,
			new_coords_size,
			&vertex_handles,
			&new_vertex_handles_allocated,
			&new_vertex_handles_size,
			&error );
    TEST_ASSERT( iBase_SUCCESS == error );

    int status = 0;
    iBase_EntityHandle tetrahedron;
    iMesh_createEnt( mesh,
		     iMesh_TETRAHEDRON,
		     vertex_handles,
		     new_vertex_handles_size,
		     &tetrahedron,
		     &status,
		     &error );  
    TEST_ASSERT( iBase_SUCCESS == error );

    // Integrate the tetrahedron.
    int cardinality = 4;
    int values_size = cardinality*quadrature->getNumPoints();
    Teuchos::ArrayRCP<double> integrated_values;
    Teuchos::ArrayRCP<double> values( values_size, 1.0 );
    quadrature->integrate( integrated_values,
			   values,
			   cardinality,
			   mesh,
			   tetrahedron );
    std::cout << "CELL INTEGRAL " << integrated_values[0] << std::endl;
    TEST_ASSERT( (int) integrated_values.size() == 1 );
    TEST_ASSERT( integrated_values[0] == 1.0 );
}

//---------------------------------------------------------------------------//
//                        end of tstIntrepidQuadrature.cpp
//---------------------------------------------------------------------------//
