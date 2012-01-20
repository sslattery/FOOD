//---------------------------------------------------------------------------//
// \file DFuncKernel.hpp
// \author
// \brief Protocol definition for distribution function kernels.
//---------------------------------------------------------------------------//

#ifndef FOOD_DFUNCKERNEL_HPP
#define FOOD_DFUNCKERNEL_HPP

#include <iMesh.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>

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
    typedef Teuchos::ArrayRCP<ScalarType>                 ScalarArray;
    typedef Intrepid::Basis<ScalarType,ScalarArray>       Basis_t;
    typedef Teuchos::RCP<Basis_t>                         RCP_Basis;

  private:

    // The basis for this kernel.
    RCP_Basis d_basis;

  public:

    // Constructor.
    DFuncKernel( int entity_topology,
	         int basis_degree,
		 int basis_operator_type );

    // Destructor.
    ~DFuncKernel();

    // Evaluate the degrees of freedom for this kernel at a specific
    // location. Coordinates are parametric.
    void evaluateDF( ScalarArray &dfunc_values,
		     const ScalarArray &coords );

    // Evaluate the gradient of the degrees of freedom for this kernel at a
    // specific location. Coordinates are parametric.
    void evaluateGradDF( ScalarArray &dfunc_grad_values,
			 const ScalarArray &coords );

    //! Get the basis for this kernel.
    RCP_Basis getDFuncKernelBasis() const
    { return d_basis; }
};

}

#endif // end FOOD_DFUNCKERNEL_HPP

//---------------------------------------------------------------------------//
// end DFuncKernel.hpp
//---------------------------------------------------------------------------//

