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
#include <Intrepid_Cubature.hpp>

namespace FOOD 
{

template<class Scalar>
class DFuncKernel
{
  
  public:

    //@{
    //! Typedefs.
    typedef Intrepid::FieldContainer<Scalar>          MDArray;
    typedef Intrepid::Basis<Scalar,MDArray>           Basis_t;
    typedef Teuchos::RCP<Basis_t>                     RCP_Basis;
    typedef Intrepid::Cubature<Scalar>                Cubature_t;
    typedef Teuchos::RCP<Cubature_t>                  RCP_Cubature;

  private:

    // The basis for this kernel.
    RCP_Basis d_basis;

    // The integration rule for this kernel.
    RCP_Cubature d_cubature;

  public:

    // Constructor.
    DFuncKernel( const int entity_topology,
		 const int discretization_type,
		 const int basis_operator_type,
	         const int basis_degree,
		 const int cubature_degree );

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

    //! Get the integration rule for this kernel.
    RCP_Cubature getDFuncKernelCubature() const
    { return d_cubature; }
};

} // end namespace FOOD

#include "DFuncKernel_Def.hpp"

#endif // end FOOD_DFUNCKERNEL_HPP

//---------------------------------------------------------------------------//
// end DFuncKernel.hpp
//---------------------------------------------------------------------------//

