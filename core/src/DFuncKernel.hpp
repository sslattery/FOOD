//---------------------------------------------------------------------------//
// \file DFuncKernel.hpp
// \author Stuart Slattery
// \brief Distribution function kernel declaration.
//---------------------------------------------------------------------------//

#ifndef FOOD_DFUNCKERNEL_HPP
#define FOOD_DFUNCKERNEL_HPP

#include <iMesh.h>

#include <Teuchos_RCP.hpp>

#include <Intrepid_FieldContainer.hpp>
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
    typedef Intrepid::FieldContainer<ScalarType>          MDArray;
    typedef Intrepid::Basis<ScalarType,MDArray>           Basis_t;
    typedef Teuchos::RCP<Basis_t>                         RCP_Basis;

  private:

    // The basis for this kernel.
    RCP_Basis d_basis;

  public:

    // Constructor.
    DFuncKernel( const int entity_topology,
		 const int discretization_type,
		 const int basis_operator_type,
	         const int basis_degree );

    // Destructor.
    ~DFuncKernel();

    // Evaluate the degrees of freedom for this kernel at a specific
    // location. Coordinates are parametric.
    void evaluateDF( MDArray &dfunc_values, const MDArray &coords );

    // Evaluate the gradient of the degrees of freedom for this kernel at a
    // specific location. Coordinates are parametric.
    void evaluateGradDF( MDArray &dfunc_grad_values, const MDArray &coords );

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

