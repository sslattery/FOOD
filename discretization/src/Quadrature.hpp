//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file Quadrature.hpp
 * \author Stuart Slattery
 * \brief Interface definition for quadrature rules.
 */
//---------------------------------------------------------------------------//

#ifndef FOOD_QUADRATURE_HPP
#define FOOD_QUADRATURE_HPP

#include <cstdlib>

#include "DFuncKernel.hpp"

#include <iBase.h>
#include <iMesh.h>

#include <Teuchos_ArrayRCP.hpp>

namespace FOOD
{

/*
 * \brief Interface definition for quadrature rules associated with an iMesh
 * topology. This is heavily based on the Intrepid defintion as the Intrepid
 * data model is fairly complete. Templated on the degree of freedom data
 * type. 
 */
template<class Scalar>
class Quadrature
{

  public:

    //@{
    //! Typedefs.
    typedef DFuncKernel<Scalar>                     DFuncKernel_t;
    typedef Teuchos::RCP<DFuncKernel_t>             RCP_DFuncKernel;
    //@}

  protected:

    // Number of quadrature rule points.
    int b_num_points;

    // Quadrature rule dimension.
    int b_dimension;

    // Quadrature rule degree.
    int b_degree;

    // Quadrature cell type. (enum)
    std::size_t b_type;

    // Quadrature cell topology. (enum)
    std::size_t b_topology;

  public:

    /*!
     * \brief Default constructor.
     */
    Quadrature()
    { /* ... */ }

    /*!
     * \brief Destructor.
     */
    virtual ~Quadrature()
    { /* ... */ }

    /*!
     * \brief Get the quadrature rule.
     * \param coordinates Interleaved quadrature point coordinates.
     * \param weights Quadrature weights returned in the quadrature point
     * order. 
     */
    virtual void 
    getQuadratureRule( Teuchos::ArrayRCP<Scalar> &coordinates,
		       Teuchos::ArrayRCP<Scalar> &weights ) const = 0;

    /*!
     * \brief Integrate over a cell.
     * \param integrated_values Integration result.
     * \param values Function coefficients.
     * \param physical_cell The entity to integrate over.
     */
    virtual void 
    integrate( Teuchos::ArrayRCP<Scalar> &integrated_values,
	       const Teuchos::ArrayRCP<Scalar> &values,
	       const RCP_DFuncKernel &dfunckernel,
	       const iMesh_Instance mesh,
	       const iBase_EntityHandle physical_cell ) = 0;

    //! Get the number of quadrature points.
    virtual int getNumPoints() const
    { return b_num_points; }

    //! Get the quadrature rule dimension.
    virtual int getDimension() const
    { return b_dimension; }

    //! Get the quadrature rule degree.
    virtual int getDegree() const
    { return b_degree; }

    //! Get the quadrature rule cell type.
    virtual int getEntityType() const 
    { return b_type; }

    //! Get the quadrature rule cell topology.
    virtual int getEntityTopology() const
    { return b_topology; }
};

} // end namespace FOOD

#endif // end FOOD_QUADRATURE_HPP

//---------------------------------------------------------------------------//
// end Quadrature.hpp
//---------------------------------------------------------------------------//
