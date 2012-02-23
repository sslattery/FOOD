//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file DFuncKernel.hpp
 * \author Stuart Slattery
 * \brief Interface definition for distribution function kernels.
 */
//---------------------------------------------------------------------------//

#ifndef FOOD_DFUNCKERNEL_HPP
#define FOOD_DFUNCKERNEL_HPP

#include <cstdlib>

#include "DiscretizationTypes.hpp"

#include <iBase.h>
#include <iMesh.h>

#include <Teuchos_ArrayRCP.hpp>

namespace FOOD
{

/*
 * \brief Interface definition for distribution function kernels associated
 * with an iMesh topology. This is heavily based on the Intrepid defintion as
 * the Intrepid data model is fairly complete. Templated on the degree of
 * freedom data type.
 */
template<class Scalar>
class DFuncKernel
{

  protected:

    // Distribution function kernel cardinality.
    int b_cardinality;

    // Distribution function kernel degree.
    int b_degree;

    // Distribution function kernel cell topology. (enum)
    std::size_t b_topology;

    // Reference cell canonical numbering. (enum)
    std::size_t b_cn_type;

    // Distribution function kernel coordinate type. (enum)
    std::size_t b_coord_type;

    // Discretization type. (enum)
    std::size_t b_discretization_type;

    // Function space type. (enum)
    std::size_t b_function_space_type;

  public:

    /*!
     * \brief Default constructor.
     */
    DFuncKernel()
    { /* ... */ }

    /*!
     * \brief Destructor.
     */
    virtual ~DFuncKernel()
    { /* ... */ }

    /*!
     * \brief Evaluate the value of the distribution function kernel at a
     * given set of parametric coordinates. 
     */
    virtual void dfuncValue( Teuchos::ArrayRCP<Scalar> values, 
			     const double param_coords[3] ) = 0;

    /*!
     * \brief Evaluate the operator value of a distribution function kernel at
     * a given set of parametric coordinates. The operator is defined by the
     * function space. (i.e. FOOD_HDIV space returns the divergence of the
     * distribution function kernel.) 
     */
    virtual void dfuncOperator( Teuchos::ArrayRCP<Scalar> values,
				const double param_coords[3] ) = 0;

    /*!
     * \brief Transform a point from the physical frame in a physical cell to
     * the reference frame of the reference cell for the distribution function
     * kernel topology. 
    */
    virtual void transformPoint( double param_coords[3],
				 const double physical_coords[3],
				 const iMesh_Instance mesh,
				 const iBase_EntityHandle physical_cell ) = 0;

    /*!
     * \brief Transform the value of the distribution function kernel at a
     * given set of parametric coordinates back to the physical frame for the
     * given physical cell. 
     */
    virtual void transformValue( Teuchos::ArrayRCP<Scalar> transformed_values, 
				 const Teuchos::ArrayRCP<Scalar> values,
				 const double param_coords[3],
				 const iMesh_Instance mesh,
				 const iBase_EntityHandle physical_cell ) = 0;

    /*!
     * \brief Transform the operator value of a distribution function kernel
     * at a given set of parametric coordinates back to the physical frame for
     * the given physical cell. The operator is defined by the function
     * space. (i.e. FOOD_HDIV space returns the divergence of the distribution
     * function kernel.)  
     */
    virtual void transformOperator( Teuchos::ArrayRCP<Scalar> transformed_values,
				    const Teuchos::ArrayRCP<Scalar> values,
				    const double param_coords[3],
				    const iMesh_Instance mesh,
				    const iBase_EntityHandle physical_cell ) = 0;

    /*!
     * \brief Evaluate the distribution function using function coefficients
     * and physical frame distribution function kernel values.
     */
    virtual void evaluate( Teuchos::ArrayRCP<Scalar> function_values,
			   const Teuchos::ArrayRCP<Scalar> coeffs,
			   const Teuchos::ArrayRCP<Scalar> dfunc_values ) = 0;

    //! Get the distribution function kernel cardinality.
    virtual int getCardinality() const
    { return b_cardinality; }

    //! Get the distribution function kernel degree.
    virtual int getDegree() const
    { return b_degree; }

    //! Get the distribution function kernel cell topology.
    virtual int getTopology() const
    { return b_topology; }

    //! Get the canonical numbering scheme for the reference cell.
    virtual int getCNType() const
    { return b_cn_type; }

    //! Get the distribution function kernel coordinate type.
    virtual int getCoordType() const
    { return b_coord_type; }

    //! Get the discretization type.
    virtual int getDiscretizationType() const
    { return b_discretization_type; }

    //! Get the function space type.
    virtual int getFunctionSpaceType() const
    { return b_function_space_type; }
};

} // end namespace FOOD

#endif // end FOOD_DFUNCKERNEL_HPP

//---------------------------------------------------------------------------//
// end DFuncKernel.hpp
//---------------------------------------------------------------------------//
