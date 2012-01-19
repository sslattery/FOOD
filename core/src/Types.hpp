//---------------------------------------------------------------------------//
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

} // end namespace FOOD

#endif // end FOOD_TYPES_HPP

//---------------------------------------------------------------------------//
// end Types.hpp
//---------------------------------------------------------------------------//
