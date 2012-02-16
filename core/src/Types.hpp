//---------------------------------------------------------------------------//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
//
// \file   Types.hpp
// \author Stuart Slattery
// \brief  Enumerated types for FOOD.
//---------------------------------------------------------------------------//

#ifndef FOOD_TYPES_HPP
#define FOOD_TYPES_HPP

namespace FOOD
{

// Field precision.
enum FOOD_Precision {
    FOOD_Precision_MIN = 0,
    FOOD_BOOLEAN =  FOOD_Precision_MIN,
    FOOD_UCHAR,
    FOOD_INTEGER,
    FOOD_FLOAT,
    FOOD_DOUBLE,
    FOOD_QUAD,
    FOOD_Precision_MAX =  FOOD_QUAD
};

// Algorithm type.
enum FOOD_AlgType {
    FOOD_AlgType_MIN = 0,
    FOOD_LOGICAL =  FOOD_AlgType_MIN,
    FOOD_INTEGRAL,
    FOOD_REAL,
    FOOD_COMPLEX,
    FOOD_AlgType_MAX =  FOOD_COMPLEX
};

// Coordinate type.
enum FOOD_CoordType {
    FOOD_CoordType_MIN = 0,
    FOOD_CARTESIAN =  FOOD_CoordType_MIN,
    FOOD_CYLINDRICAL,
    FOOD_SPHERICAL,
    FOOD_CoordType_MAX =  FOOD_SPHERICAL
};

// Storage order.
enum FOOD_StorageHint {
    FOOD_StorageHint_MIN = 0,
    FOOD_BLOCKED =  FOOD_StorageHint_MIN,
    FOOD_INTERLEAVED,
    FOOD_MIXED,
    FOOD_PER_ENTITY,
    FOOD_StorageHint_MAX =  FOOD_PER_ENTITY
};

// Distribution function kernel discretization type.
enum FOOD_DiscretizationType {
    FOOD_DiscretizationType_MIN = 0,
    FOOD_FEM = FOOD_DiscretizationType_MIN,
    FOOD_FV,
    FOOD_FD,
    FOOD_DiscretizationType_MAX = FOOD_FD
};

// Basic operator types.
enum FOOD_BasisOperatorType {
    FOOD_BasisOperatorType_MIN = 0,
    FOOD_GRAD = FOOD_BasisOperatorType_MIN,
    FOOD_DIV,
    FOOD_CURL,
    FOOD_BasisOperatorType_MAX = FOOD_CURL
};

// Basis function space types.
enum FOOD_FunctionSpaceType {
    FOOD_FunctionSpaceType_MIN = 0,
    FOOD_HGRAD = FOOD_FunctionSpaceType_MIN,
    FOOD_HDIV,
    FOOD_HCURL,
    FOOD_FunctionSpaceType_MAX = FOOD_HCURL
};

} // end namespace FOOD

#endif // end FOOD_TYPES_HPP

//---------------------------------------------------------------------------//
// end Types.hpp
//---------------------------------------------------------------------------//
