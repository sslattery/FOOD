//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file IntrepidKernel.hpp
 * \author Stuart Slattery
 * \brief Distribution function kernel declaration for Intrepid basis
 * implemenentations.
 */
//---------------------------------------------------------------------------//

#ifndef FOOD_INTREPIDKERNEL_HPP
#define FOOD_INTREPIDKERNEL_HPP

#include "DFuncKernel.hpp"

#include <iBase.h>
#include <iMesh.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include <Intrepid_Basis.hpp>
#include <Intrepid_FieldContainer.hpp>

namespace FOOD
{

template<class Scalar>
class IntrepidKernel : public DFuncKernel<Scalar>
{

  public:

    //@{
    //! Typedefs.
    typedef Intrepid::FieldContainer<Scalar>        MDArray;
    typedef Intrepid::Basis<Scalar,MDArray>         IntrepidBasis_t;
    typedef Teuchos::RCP<IntrepidBasis_t>           RCP_IntrepidBasis;
    //@}

  private:

    // The Intrepid basis implementation for this kernel.
    RCP_IntrepidBasis d_intrepid_basis;

  public:

    // Constructor.
    IntrepidKernel( RCP_IntrepidBasis intrepid_basis,
		    const int entity_type, 
		    const int entity_topology, 
		    const int discretization_type,
		    const int function_space_type );

    // Destructor.
    ~IntrepidKernel();

    // Evaluate the value of the distribution function kernel at a given set
    // of parametric coordinates.  
    void dfuncValue( Teuchos::ArrayRCP<Scalar> &values, 
		     const double param_coords[3] );

    // Evaluate the operator value of a distribution function kernel at a
    // given set of parametric coordinates. 
    void dfuncOperator( Teuchos::ArrayRCP<Scalar> &values, 
			const double param_coords[3] );

    // Transform a point from the physical frame in a physical cell to the
    // reference frame of the reference cell for the distribution function
    // kernel topology.
    void transformPoint( double param_coords[3],
			 const double physical_coords[3],
			 const iMesh_Instance mesh,
			 const iBase_EntityHandle physical_cell );

    // Transform the value of the distribution function kernel at a given set
    // of parametric coordinates back to the physical frame for the given
    // physical cell. 
    void transformValue( Teuchos::ArrayRCP<Scalar> &transformed_values, 
			 const Teuchos::ArrayRCP<Scalar> &values,
			 const double param_coords[3],
			 const iMesh_Instance mesh,
			 const iBase_EntityHandle physical_cell );

    // Transform the operator value of a distribution function kernel at a
    // given set of parametric coordinates back to the physical frame for the
    // given physical cell.
    void transformOperator( Teuchos::ArrayRCP<Scalar> &transformed_values,
			    const Teuchos::ArrayRCP<Scalar> &values,
			    const double param_coords[3],
			    const iMesh_Instance mesh,
			    const iBase_EntityHandle physical_cell );

    // Evaluate the distribution function using function coefficients and
    // physical frame distribution function kernel values. 
    void evaluate( Teuchos::ArrayRCP<Scalar> &function_values,
		   const Teuchos::ArrayRCP<Scalar> &coeffs,
		   const Teuchos::ArrayRCP<Scalar> &dfunc_values );
};

} // end namespace FOOD

#include "IntrepidKernel_Def.hpp"

#endif // end FOOD_INTREPIDKERNEL_HPP

//---------------------------------------------------------------------------//
// end IntrepidKernel.hpp
//---------------------------------------------------------------------------//
