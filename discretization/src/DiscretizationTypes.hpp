//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file   DiscretizationTypes.hpp
 * \author Stuart Slattery
 * \brief  Enumerated types for the FOOD discretization subpackage.
 */
//---------------------------------------------------------------------------//

#ifndef FOOD_DISCRETIZATIONTYPES_HPP
#define FOOD_DISCRETIZATIONTYPES_HPP

namespace FOOD
{

// Coordinate type.
enum FOOD_CoordType {
    FOOD_CoordType_MIN = 0,
    FOOD_CARTESIAN =  FOOD_CoordType_MIN,
    FOOD_CYLINDRICAL,
    FOOD_SPHERICAL,
    FOOD_CoordType_MAX =  FOOD_SPHERICAL
};

// Distribution function kernel discretization type.
enum FOOD_DiscretizationType {
    FOOD_DiscretizationType_MIN = 0,
    FOOD_FEM = FOOD_DiscretizationType_MIN,
    FOOD_FV,
    FOOD_FD,
    FOOD_DiscretizationType_MAX = FOOD_FD
};

// Basis function space types.
enum FOOD_FunctionSpaceType {
    FOOD_FunctionSpaceType_MIN = 0,
    FOOD_HGRAD = FOOD_FunctionSpaceType_MIN,
    FOOD_HDIV,
    FOOD_HCURL,
    FOOD_FunctionSpaceType_MAX = FOOD_HCURL
};

// Canonical numbering system types.
enum FOOD_CNType {
    FOOD_CNType_MIN = 0,
    FOOD_MBCN = FOOD_CNType_MIN,
    FOOD_SHARDSCN,
    FOOD_CNType_MAX = FOOD_SHARDSCN
};

} // end namespace FOOD

#endif // end FOOD_DISCRETIZATIONTYPES_HPP

//---------------------------------------------------------------------------//
// end DiscretizationTypes.hpp
//---------------------------------------------------------------------------//
