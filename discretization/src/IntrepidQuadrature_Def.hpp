//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file IntrepidQuadrature_Def.hpp
 * \author Stuart Slattery
 * \brief Quadrature definition for Intrepid cubature rule implementations.
 */
//---------------------------------------------------------------------------//

#ifndef FOOD_INTREPIDQUADRATURE_DEF_HPP
#define FOOD_INTREPIDQUADRATURE_DEF_HPP

#include "CellTopologyFactory.hpp"

#include <Teuchos_Tuple.hpp>

#include <Shards_CellTopology.hpp>

#include <Intrepid_CellTools.hpp>
#include <Intrepid_FunctionSpaceTools.hpp>

namespace FOOD
{
/*!
 * \brief Constructor.
 */
template<class Scalar>
IntrepidQuadrature<Scalar>::IntrepidQuadrature( 
    RCP_IntrepidCubature intrepid_cubature,
    const int degree,
    const int entity_type, 
    const int entity_topology )
    : d_intrepid_cubature( intrepid_cubature )
{
    this->b_num_points = d_intrepid_cubature->getNumPoints();
    this->b_dimension = d_intrepid_cubature->getDimension();
    this->b_degree = degree;
    this->b_type = entity_type;
    this->b_topology = entity_topology;
}

/*!
v * \brief Destructor.
 */
template<class Scalar>
IntrepidQuadrature<Scalar>::~IntrepidQuadrature()
{ /* ... */ }

/*!
 * \brief Get the quadrature rule.
 */
template<class Scalar>
void IntrepidQuadrature<Scalar>::getQuadratureRule( 
    Teuchos::ArrayRCP<Scalar> &coordinates,
    Teuchos::ArrayRCP<Scalar> &weights ) const
{
    MDArray points( this->b_num_points, this->b_dimension );
    MDArray point_weights( this->b_num_points );

    d_intrepid_cubature->getCubature( points, point_weights );

    coordinates = points.getData();
    weights = point_weights.getData();
}

} // end namespace FOOD

#endif // end FOOD_INTREPIDQUADRATURE_DEF_HPP

//---------------------------------------------------------------------------//
// end IntrepidQuadrature_Def.hpp
//---------------------------------------------------------------------------//

