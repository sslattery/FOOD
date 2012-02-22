//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file   FieldTypes.hpp
 * \author Stuart Slattery
 * \brief  Enumerated types for the FOOD field subpackage.
 */
//---------------------------------------------------------------------------//

#ifndef FOOD_FIELDTYPES_HPP
#define FOOD_FIELDTYPES_HPP

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

// Algebraic type.
enum FOOD_AlgType {
    FOOD_AlgType_MIN = 0,
    FOOD_LOGICAL =  FOOD_AlgType_MIN,
    FOOD_INTEGRAL,
    FOOD_REAL,
    FOOD_COMPLEX,
    FOOD_AlgType_MAX =  FOOD_COMPLEX
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

#endif // end FOOD_FIELDTYPES_HPP

//---------------------------------------------------------------------------//
// end FieldTypes.hpp
//---------------------------------------------------------------------------//
