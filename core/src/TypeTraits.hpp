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
    static const std::size_t tag_type  = iBase_BYTES;
    static const std::size_t tag_size  = sizeof(bool);
    static const std::size_t precision = FOOD_BOOLEAN;
};


template<> 
struct TypeTraits<char> 
{
    static const std::size_t tag_type  = iBase_BYTES;
    static const std::size_t tag_size  = sizeof(char);
    static const std::size_t precision = FOOD_UCHAR;
};

template<> 
struct TypeTraits<unsigned char> 
{
    static const std::size_t tag_type  = iBase_BYTES;
    static const std::size_t tag_size  = sizeof(unsigned char);
    static const std::size_t precision = FOOD_UCHAR;
};

template<> 
struct TypeTraits<short> 
{
    static const std::size_t tag_type  = iBase_INTEGER;
    static const std::size_t tag_size  = 1;
    static const std::size_t precision = FOOD_INTEGER;
};

template<> 
struct TypeTraits<unsigned short> 
{
    static const std::size_t tag_type  = iBase_INTEGER;
    static const std::size_t tag_size  = 1;
    static const std::size_t precision = FOOD_INTEGER;
};

template<> 
struct TypeTraits<int> 
{
    static const std::size_t tag_type  = iBase_INTEGER;
    static const std::size_t tag_size  = 1;
    static const std::size_t precision = FOOD_INTEGER;
};

template<> 
struct TypeTraits<unsigned int> 
{
    static const std::size_t tag_type  = iBase_INTEGER;
    static const std::size_t tag_size  = 1;
    static const std::size_t precision = FOOD_INTEGER;
};

template<> 
struct TypeTraits<long> 
{
    static const std::size_t tag_type  = iBase_INTEGER;
    static const std::size_t tag_size  = 1;
    static const std::size_t precision = FOOD_INTEGER;
};

template<> 
struct TypeTraits<unsigned long> 
{
    static const std::size_t tag_type  = iBase_INTEGER;
    static const std::size_t tag_size  = 1;
    static const std::size_t precision = FOOD_INTEGER;
};

template<> 
struct TypeTraits<float> 
{
    static const std::size_t tag_type  = iBase_DOUBLE;
    static const std::size_t tag_size  = 1;
    static const std::size_t precision = FOOD_FLOAT;
};

template<> 
struct TypeTraits<double> 
{
    static const std::size_t tag_type  = iBase_DOUBLE;
    static const std::size_t tag_size  = 1;
    static const std::size_t precision = FOOD_DOUBLE;
};

template<> 
struct TypeTraits<long double> 
{
    static const std::size_t tag_type  = iBase_DOUBLE;
    static const std::size_t tag_size  = 1;
    static const std::size_t precision = FOOD_QUAD;
};

} // end namepsace FOOD

#endif // end FOOD_TYPETRAITS_HPP

//---------------------------------------------------------------------------//
// end TypeTraits.hpp
//---------------------------------------------------------------------------//
