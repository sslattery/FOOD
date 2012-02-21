//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file ConsistentScheme.hpp
 * \author Stuart Slattery
 * \brief Consistent finite element interpolation scheme declaration.
 */
//---------------------------------------------------------------------------//

#ifndef FOOD_CONSISTENTSCHEME_HPP
#define FOOD_CONSISTENTSCHEME_HPP

#include <map>

#include "TensorField.hpp"
#include "KDTree.hpp"

#include <iBase.h>
#include <iMesh.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include <Intrepid_FieldContainer.hpp>

namespace FOOD
{

template<class Scalar>
class ConsistentScheme
{

  public:

    //@{
    //! Typedefs.
    typedef TensorField<Scalar>                      TensorField_t;
    typedef Teuchos::RCP<TensorField_t>              RCP_TensorField;
    typedef Teuchos::RCP< KDTree<3> >                RCP_KDTree;
    typedef Intrepid::FieldContainer<Scalar>         MDArray;
    //@}

  private:

    // Degrees of freedom domain.
    RCP_TensorField d_dof_domain;

    // Degrees of freedom range.
    RCP_TensorField d_dof_range;

    // KDTree.
    RCP_KDTree d_kdtree;

    // Range to domain mapping.
    std::map<iBase_EntityHandle,iBase_EntityHandle> d_map;

  public:

    // Constructor.
    ConsistentScheme( RCP_TensorField dof_domain, RCP_TensorField dof_range );

    // Destructor.
    ~ConsistentScheme();

    // Setup for interpolation.
    void setup();

    // Perform value interpolation of the degrees of freedom from the domain
    // to the range.
    void interpolateValueDF();

    // Perform gradient interpolation of the degrees of freedom from the domain
    // to the range.
    void interpolateGradDF();
};

} // end namespace FOOD

#include "ConsistentScheme_Def.hpp"

#endif // end FOOD_CONSISTENTSCHEME_HPP

//---------------------------------------------------------------------------//
// end ConsistentScheme.hpp
//---------------------------------------------------------------------------//
