//---------------------------------------------------------------------------//
// \file DFuncKernel.hpp
// \author Stuart Slattery
// \brief Distribution function kernel declaration.
//---------------------------------------------------------------------------//

#ifndef FOOD_DFUNCKERNEL_HPP
#define FOOD_DFUNCKERNEL_HPP

#include <iMesh.h>

#include <Teuchos_RCP.hpp>

#include <Shards_Array.hpp>

#include <Intrepid_Basis.hpp>

namespace FOOD 
{

template<class ScalarType_T>
class DFuncKernel
{
  
  public:

    //@{
    //! Typedefs.
    typedef ScalarType_T                                  ScalarType;
    typedef Shards_Array<ScalarType>                      MDArray;
    typedef Intrepid::Basis<ScalarType,MDArray>           Basis_t;
    typedef Teuchos::RCP<Basis_t>                         RCP_Basis;

  private:

    // The basis for this kernel.
    RCP_Basis d_basis;

  public:

    // Constructor.
    DFuncKernel( const int entity_topology,
		 const int discretization_type,
	         const int basis_degree,
		 const int basis_operator_type );

    // Destructor.
    ~DFuncKernel();

    // Evaluate the degrees of freedom for this kernel at a specific
    // location. Coordinates are parametric.
    void evaluateDF( MDarray &dfunc_values, const MDarray &coords );

    // Evaluate the gradient of the degrees of freedom for this kernel at a
    // specific location. Coordinates are parametric.
    void evaluateGradDF( MDarray &dfunc_grad_values, const MDarray &coords );

    //! Get the basis for this kernel.
    RCP_Basis getDFuncKernelBasis() const
    { return d_basis; }
};

} // end namespace FOOD

#include "DFuncKernel_Def.hpp"

#endif // end FOOD_DFUNCKERNEL_HPP

//---------------------------------------------------------------------------//
// end DFuncKernel.hpp
//---------------------------------------------------------------------------//

