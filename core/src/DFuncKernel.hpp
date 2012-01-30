//---------------------------------------------------------------------------//
// \file DFuncKernel.hpp
// \author Stuart Slattery
// \brief Distribution function kernel declaration.
//---------------------------------------------------------------------------//

#ifndef FOOD_DFUNCKERNEL_HPP
#define FOOD_DFUNCKERNEL_HPP

#include "Types.hpp"

#include <iBase.h>
#include <iMesh.h>

#include <Teuchos_RCP.hpp>

#include <Intrepid_FieldContainer.hpp>
#include <Intrepid_Basis.hpp>

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
    //@}

  private:

    // Entity type for this kernel.
    std::size_t d_entity_type;

    // Entity topology for this kernel.
    std::size_t d_entity_topology;

    // The coordinate type for this kernel.
    std::size_t d_coordinate_type;

    // Discretization type for this kernel.
    std::size_t d_discretization_type;

    // Basis function space for this kernel.
    std::size_t d_basis_function_space;

    // The basis for this kernel.
    RCP_Basis d_basis;

    // The cell topology for this kernel.
    RCP_CellTopology d_cell_topology;

  public:

    // Constructor.
    DFuncKernel( const int entity_type,
	         const int entity_topology,
		 const int coordinate_type,
		 const int discretization_type,
		 const int basis_function_space,
	         const int basis_degree );

    // Destructor.
    ~DFuncKernel();

    // Evaluate the basis value for this kernel at a specific parametric
    // location. 
    void evaluateValueBasis( MDArray &values_at_coords, 
			     const MDArray &coords );

    // Evaluate the gradient of the basis for this kernel at a specific
    // parametric location.
    void evaluateGradBasis( MDArray &dfunc_grad_values, const MDArray &coords );

    // Evaluate the divergence of the basis for this kernel at a specific
    // parametric location.
    void evaluateDivBasis( MDArray &dfunc_grad_values, const MDArray &coords );

    // Evaluate the curl of the basis for this kernel at a specific parametric
    // location. 
    void evaluateCurlBasis( MDArray &dfunc_grad_values, const MDArray &coords );

    //! Get the entity type for the kernel.
    int getEntityType() const
    { return d_entity_type; }

    //! Get the entity topology for the kernel.
    int getEntityTopology() const
    { return d_entity_topology; }

    //! Get the discretization type.
    int getDiscretizationType() const
    { return d_discretization_type; }

    //! Get the basis function space for this kernel.
    int getBasisFunctionSpace() const
    { return d_basis_function_space; }

    //! Get the basis for this kernel.
    RCP_Basis getBasis() const
    { return d_basis; }

    //! Get the basis cardinality.
    int getBasisCardinality() const
    { return d_basis->getCardinality(); }

    //! Get the basis degree.
    int getBasisDegree() const
    { return d_basis->getDegree(); }

    //! Get the reference cell coordinates for the degrees of freedom.
    //! MDArray(F,D) defined for the basis.
    void getBasisCoords( MDArray &coords )
    { return d_basis->getDoFCoords( coords ); }

    //! Get the cell topology for this kernel.
    RCP_CellTopology getCellTopology() const
    { return d_cell_topology; }
};

} // end namespace FOOD

#include "DFuncKernel_Def.hpp"

#endif // end FOOD_DFUNCKERNEL_HPP

//---------------------------------------------------------------------------//
// end DFuncKernel.hpp
//---------------------------------------------------------------------------//

