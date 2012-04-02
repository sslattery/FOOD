//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file IntrepidQuadrature.hpp
 * \author Stuart Slattery
 * \brief Quadrature declaration for Intrepid cubature rule
 * implemenentations. 
 */
//---------------------------------------------------------------------------//

#ifndef FOOD_INTREPIDQUADRATURE_HPP
#define FOOD_INTREPIDQUADRATURE_HPP

#include "Quadrature.hpp"

#include <iBase.h>
#include <iMesh.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include <Intrepid_Cubature.hpp>
#include <Intrepid_FieldContainer.hpp>

namespace FOOD
{

template<class Scalar>
class IntrepidQuadrature : public Quadrature<Scalar>
{
 
  public:

    //@{
    //! Typedefs.
    typedef Intrepid::FieldContainer<Scalar>        MDArray;
    typedef Intrepid::Cubature<Scalar>              IntrepidCubature_t;
    typedef Teuchos::RCP<IntrepidCubature_t>        RCP_IntrepidCubature;
    //@}

  private:

    // The Intrepid cubature implementation for this quadrature.
    RCP_IntrepidCubature d_intrepid_cubature;

  public:

    // Constructor.
    IntrepidQuadrature( RCP_IntrepidCubature intrepid_cubature,
			const int entity_type, 
			const int entity_topology );

    // Destructor.
    ~IntrepidQuadrature();

    // Get the quadrature rule.
    void getQuadratureRule( Teuchos::ArrayRCP<Scalar> &coordinates,
			    Teuchos::ArrayRCP<Scalar> &weights ) const;

    // Integrate over a cell.
    void integrate( Teuchos::ArrayRCP<Scalar> &integrated_values,
		    const Teuchos::ArrayRCP<Scalar> &values,
		    const int cardinality,
		    const iMesh_Instance mesh,
		    const iBase_EntityHandle physical_cell );
};

} // end namespace FOOD

#include "IntrepidQuadrature_Def.hpp"

#endif // end FOOD_INTREPIDQUADRATURE_HPP

//---------------------------------------------------------------------------//
// end IntrepidQuadrature.hpp
//---------------------------------------------------------------------------//

