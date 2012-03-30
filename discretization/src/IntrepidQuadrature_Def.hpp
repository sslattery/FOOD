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

namespace FOOD
{
/*!
 * \brief Constructor.
 */
template<class Scalar>
IntrepidQuadrature<Scalar>::IntrepidQuadrature( 
    RCP_IntrepidCubature intrepid_cubature,
    const int entity_type, 
    const int entity_topology )
    : d_intrepid_cubature( intrepid_cubature )
{
    this->b_num_points = d_intrepid_cubature->getNumPoints();
    this->b_dimension = d_intrepid_cubature->getDimension();
    this->b_type = entity_type;
    this->b_topology = entity_topology;
}

/*!
 * \brief Destructor.
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
    Teuchos::Tuple<int,2> points_dimensions;
    points_dimensions[0] = this->b_num_points;
    points_dimensions[1] = this->b_dimension;
    MDArray points( points_dimensions );
    
    Teuchos::Tuple<int,1> point_weights_dimensions;
    point_weights_dimensions[0] = this->b_num_points;
    MDArray point_weights( point_weights_dimensions );

    d_intrepid_cubature->getCubature( points, point_weights );

    coordinates = points.getData();
    weights = point_weights.getData();
}

} // end namespace FOOD

#endif // end FOOD_INTREPIDQUADRATURE_DEF_HPP

//---------------------------------------------------------------------------//
// end IntrepidQuadrature_Def.hpp
//---------------------------------------------------------------------------//

