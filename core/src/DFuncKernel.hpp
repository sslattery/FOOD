//---------------------------------------------------------------------------//
// \file DFuncKernel.hpp
// \author Stuart Slattery
// \brief Distribution function kernel declaration.
//---------------------------------------------------------------------------//

#ifndef FOOD_DFUNCKERNEL_HPP
#define FOOD_DFUNCKERNEL_HPP

#include <iMesh.h>

#include <Teuchos_RCP.hpp>

#include <Shards_CellTopology.hpp>

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
    typedef Teuchos::RCP<shards::CellTopology>        RCP_CellTopology;
    typedef Intrepid::Cubature<Scalar>                Cubature_t;
    typedef Teuchos::RCP<Cubature_t>                  RCP_Cubature;

  private:

    // Entity topology for this kernel.
    std::size_t d_entity_topology;

    // Discretization type for this kernel.
    std::size_t d_discretization_type;

    // Basis operator space for this kernel.
    std::size_t d_basis_operator_type;

    // Basis degree for this kernel.
    std::size_t d_basis_degree;

    // Cubature degree for this kernel.
    std::size_t d_cubature_degree;

    // The basis for this kernel.
    RCP_Basis d_basis;

    // The cell topology for this kernel.
    RCP_CellTopology d_cell_topology;

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
    // location.
    void evaluateDF( MDArray &values_at_coords, 
		     const MDArray &coords );

    // Evaluate the gradient of the degrees of freedom for this kernel at a
    // specific location.
    void evaluateGradDF( MDArray &dfunc_grad_values, const MDArray &coords );

    //! Get the entity topology for the kernel.
    int getEntityTopology() const
    { return d_entity_topology; }

    //! Get the discretization type.
    int getDiscretizationType() const
    { return d_discretization_type; }

    //! Get the basis operator type for this kernel.
    int getBasisOperator() const
    { return d_basis_operator_type; }

    //! Get the basis degree for this kernel.
    int getBasisDegree() const
    { return d_basis_degree; }

    //! Get the cubature degree for this kernel.
    int getCubatureDegree() const
    { return d_cubature_degree; }

    //! Get the basis for this kernel.
    RCP_Basis getBasis() const
    { return d_basis; }

    //! Get the cell topology for this kernel.
    RCP_CellTopology getCellTopology() const
    { return d_cell_topology; }

    //! Get the integration rule for this kernel.
    RCP_Cubature getCubature() const
    { return d_cubature; }
};

} // end namespace FOOD

#include "DFuncKernel_Def.hpp"

#endif // end FOOD_DFUNCKERNEL_HPP

//---------------------------------------------------------------------------//
// end DFuncKernel.hpp
//---------------------------------------------------------------------------//

