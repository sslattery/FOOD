//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file InterpidKernel.hpp
 * \author Stuart Slattery
 * \brief Distribution function kernel declaration for Intrepid basis
 * implemenentations.
 */
//---------------------------------------------------------------------------//

#ifndef FOOD_INTERPIDKERNEL_HPP
#define FOOD_INTERPIDKERNEL_HPP

#include "DFuncKernel.hpp"

#include "Teuchos_RCP.hpp"

#include "Intrepid_Basis.hpp"

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
    IntrepidKernel( RCP_IntrepidBasis intrepid_basis );

    // Destructor.
    ~IntrepidKernel();

    // Evaluate the value of the distribution function kernel at a given set
    // of parametric coordinates.  
    void evaluateValue( Scalar **values, 
			int *num_values,
			const double param_coords[3] );

    // Evaluate the operator value of a distribution function kernel at a
    // given set of parametric coordinates. 
    void evaluateOperator( Scalar **values, 
			   int *num_values,
			   const double param_coords[3] );

    // Transform a point from the physical frame in a physical cell to the
    // reference frame of the reference cell for the distribution function
    // kernel topology.
    void transformPoint( double param_coords[3],
			 const double physical_coords[3],
			 const iBase_EntityHandle physical_cell );

    // Transform the value of the distribution function kernel at a given set
    // of parametric coordinates back to the physical frame for the given
    // physical cell. 
    void transformValue( Scalar **transformed_values, 
			 const Scalar *values,
			 const int num_values,
			 const double param_coords[3],
			 const iBase_EntityHandle physical_cell );

    // Transform the operator value of a distribution function kernel at a
    // given set of parametric coordinates back to the physical frame for the
    // given physical cell.
    void transformOperator( Scalar **transformed_values,
			    const Scalar *values,
			    const int num_values,
			    const double param_coords[3],
			    const iBase_EntityHandle physical_cell );
};

} // end namespace FOOD

#include "IntrepidKernel_Def.hpp"

#endif // end FOOD_INTERPIDKERNEL_HPP

//---------------------------------------------------------------------------//
// end IntrepidKernel.hpp
//---------------------------------------------------------------------------//
