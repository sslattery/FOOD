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
enum Precision {
    Precision_MIN = 0,
    BOOLEAN =  Precision_MIN,
    UCHAR,
    INTEGER,
    FLOAT,
    DOUBLE,
    QUAD,
    Precision_MAX =  QUAD
};

// Algorithm type.
enum AlgType {
    AlgType_MIN = 0,
    LOGICAL =  AlgType_MIN,
    INTEGRAL,
    REAL,
    COMPLEX,
    AlgType_MAX =  COMPLEX
};

// Coordinate type.
enum CoordType {
    CoordType_MIN = 0,
    CARTESIAN =  CoordType_MIN,
    CYLINDRICAL,
    SPHERICAL,
    CoordType_MAX =  SPHERICAL
};

// Storage order.
enum StorageHint {
    StorageHint_MIN = 0,
    BLOCKED =  StorageHint_MIN,
    INTERLEAVED,
    MIXED,
    PER_ENTITY,
    StorageHint_MAX =  PER_ENTITY
};

} // end namespace FOOD

#endif // end FOOD_TYPES_HPP

//---------------------------------------------------------------------------//
// end Types.hpp
//---------------------------------------------------------------------------//
