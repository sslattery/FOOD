//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file QuadratureFactory_Def.hpp
 * \author Stuart Slattery
 * \brief Factory method definition for quadrature rule generation.
 */
//---------------------------------------------------------------------------//

#ifndef FOOD_QUADRATUREFACTORY_DEF_HPP
#define FOOD_QUADRATUREFACTORY_DEF_HPP

#include "Exception.hpp"

#include "TopologyTools.hpp"
#include "CellTopologyFactory.hpp"
#include "IntrepidQuadrature.hpp"

#include <iMesh.h>

#include <Teuchos_ENull.hpp>

#include <Shards_CellTopology.hpp>

#include <Intrepid_DefaultCubatureFactory.hpp>

namespace FOOD
{

/*!
 * \brief Constructor.
 */
template<class Scalar>
QuadratureFactory<Scalar>::QuadratureFactory()
{ /* ... */ }

/*
 * \brief Destructor.
 */
template<class Scalar>
QuadratureFactory<Scalar>::~QuadratureFactory()
{ /* ... */ }

/*! 
 * \brief Factory method.
 */
template<class Scalar>
Teuchos::RCP< Quadrature<Scalar> > 
QuadratureFactory<Scalar>::create( const int entity_type,
				   const int entity_topology,
				   const int degree )
{
    Teuchos::RCP< Quadrature<Scalar> > new_quadrature;

    switch ( entity_topology )
    {
	case iMesh_POINT:

	    testPrecondition( iMesh_POINT != entity_topology, 
			      "Topology not supported." );

	default:

	    CellTopologyFactory cell_topo_factory;
	    Teuchos::RCP<shards::CellTopology> shards_topo = 
		cell_topo_factory.create( 
		    entity_topology,
		    TopologyTools::numLinearNodes( entity_topology ) );

	    Intrepid::DefaultCubatureFactory<Scalar> cub_factory;
	    Teuchos::RCP< Intrepid::Cubature<Scalar> > intrepid_cub =
		cub_factory.create( *shards_topo,
				    degree );

	    new_quadrature = Teuchos::rcp( new IntrepidQuadrature<Scalar>(
					       intrepid_cub,
					       degree,
					       entity_type,
					       entity_topology ) );
    }
    
    testPostcondition( new_quadrature != Teuchos::null,
		       "Failure creating quadrature rule" );

    return new_quadrature;
}

} // end namespace FOOD

#endif // end FOOD_QUADRATUREFACTORY_DEF_HPP

//---------------------------------------------------------------------------//
// end QuadratureFactory_Def.hpp
//---------------------------------------------------------------------------//
