//---------------------------------------------------------------------------//
// \file TypeTraits.hpp
// \author Stuart Slattery
// \brief TypeTraits defintion.
//---------------------------------------------------------------------------//

#ifndef FOOD_TYPETRAITS_HPP
#define FOOD_TYPETRAITS_HPP

#include "Types.hpp"

#include <iBase.h>

namespace FOOD
{

template<class T>
struct TypeTraits
{ /* ... */ };

template<>
struct TypeTraits<bool>
{
    static const int tag_type       = iBase_BYTES;
    static const int tag_size       = sizeof(bool);
    static const int precision      = FOOD_BOOLEAN;
};


template<> 
struct TypeTraits<char> 
{
    static const int tag_type       = iBase_BYTES;
    static const int tag_size       = sizeof(char);
    static const int precision      = FOOD_UCHAR;
};

template<> 
struct TypeTraits<unsigned char> 
{
    static const int tag_type       = iBase_BYTES;
    static const int tag_size       = sizeof(unsigned char);
    static const int precision      = FOOD_UCHAR;
};

template<> 
struct TypeTraits<short> 
{
    static const int tag_type       = iBase_INTEGER;
    static const int tag_size       = 1;
    static const int precision      = FOOD_INTEGER;
};

template<> 
struct TypeTraits<unsigned short> 
{
    static const int tag_type       = iBase_INTEGER;
    static const int tag_size       = 1;
    static const int precision      = FOOD_INTEGER;
};

template<> 
struct TypeTraits<int> 
{
    static const int tag_type       = iBase_INTEGER;
    static const int tag_size       = 1;
    static const int precision      = FOOD_INTEGER;
};

template<> 
struct TypeTraits<unsigned int> 
{
    static const int tag_type       = iBase_INTEGER;
    static const int tag_size       = 1;
    static const int precision      = FOOD_INTEGER;
};

template<> 
struct TypeTraits<long> 
{
    static const int tag_type       = iBase_INTEGER;
    static const int tag_size       = 1;
    static const int precision      = FOOD_INTEGER;
};

template<> 
struct TypeTraits<unsigned long> 
{
    static const int tag_type       = iBase_INTEGER;
    static const int tag_size       = 1;
    static const int precision      = FOOD_INTEGER;
};

template<> 
struct TypeTraits<float> 
{
    static const int tag_type       = iBase_DOUBLE;
    static const int tag_size       = 1;
    static const int precision      = FOOD_FLOAT;
};

template<> 
struct TypeTraits<double> 
{
    static const int tag_type       = iBase_DOUBLE;
    static const int tag_size       = 1;
    static const int precision      = FOOD_DOUBLE;
};

template<> 
struct TypeTraits<long double> 
{
    static const int tag_type       = iBase_DOUBLE;
    static const int tag_size       = 1;
    static const int precision      = FOOD_QUAD;
};

} // end namepsace FOOD

#endif // end FOOD_TYPETRAITS_HPP

//---------------------------------------------------------------------------//
// end TypeTraits.hpp
//---------------------------------------------------------------------------//
