//---------------------------------------------------------------------------//
// \field KernelTraits.hpp
// \author Stuart Slattery
// \brief Traits declaration for local field evaluation kernels.
//---------------------------------------------------------------------------//

#ifndef FOOD_KERNELTRAITS_HPP
#define FOOD_KERNELTRAITS_HPP

#include <iBase.h>
#include <iMesh.h>

#include <Phalanx_ConfigDefs.hpp>
#include <Phalanx_Allocator_New.hpp>
#include <Phalanx_Traits_Base.hpp>
#include <Phalanx_TypeStrings.hpp>

#include <Sacado.hpp>
#include <Sacado_mpl_vector.hpp>
#include <Sacado_mpl_find.hpp>

#include <boost/mpl/map.hpp>
#include <boost/mpl/find.hpp>

namespace FOOD {

struct KernelTraits : public PHX::TraitsBase {
   
    //@{
    //! Scalar types.
    typedef double                            RealType;
    typedef Sacado::Fad::DFad<double>         FadType;
    //@}
    
    //@{
    //! Evaluation Types
    struct Residual { typedef RealType ScalarT; };
    struct Jacobian { typedef FadType ScalarT;  };
    typedef Sacado::mpl::vector<Residual, Jacobian> EvalTypes;
    //@}

    //! Residual
    typedef Sacado::mpl::vector<RealType> ResidualDataTypes;
  
    //! Jacobian
    typedef Sacado::mpl::vector<FadType> JacobianDataTypes;
    
    //! Map the key EvalType a vector of DataTypes
    typedef boost::mpl::map<
	boost::mpl::pair<Residual, ResidualDataTypes>,
	boost::mpl::pair<Jacobian, JacobianDataTypes>
	>::type EvalToDataMap;

    //! Allocator Type.
    typedef PHX::NewAllocator Allocator;

    //@{
    //! Evaluation Data.
    typedef void* SetupData;
    typedef const iBase_EntityHandle* EvalData;
    typedef void* PreEvalData;
    typedef void* PostEvalData;
    //@}
};
 

} // end namespace 

// Specialize the Evaluation and Data types for the TypeString object in
// phalanx/src/Phalanx_TypeString.hpp. 
namespace PHX
{

// Data Types
template<> 
struct TypeString<double> 
{ static const std::string value; };

template<> 
struct TypeString< Sacado::Fad::DFad<double> > 
{ static const std::string value; };

// Evaluation Types
template<> 
struct TypeString<FOOD::KernelTraits::Residual> 
{ static const std::string value; };

template<> 
struct TypeString<FOOD::KernelTraits::Jacobian> 
{ static const std::string value; };

}

#endif // end FOOD_KERNELTRAITS_HPP

//---------------------------------------------------------------------------//
// end FOOD::KernelTraits.hpp
//---------------------------------------------------------------------------//

