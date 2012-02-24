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

#include "DataTransferScheme.hpp"
#include "TensorField.hpp"
#include "KDTree.hpp"

#include <iBase.h>
#include <iMesh.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>

namespace FOOD
{

template<class Scalar>
class ConsistentScheme : public DataTransferScheme<Scalar>
{

  public:

    //@{
    //! Typedefs.
    typedef TensorField<Scalar>                      TensorField_t;
    typedef Teuchos::RCP<TensorField_t>              RCP_TensorField;
    typedef Teuchos::RCP< KDTree<3> >                RCP_KDTree;
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
    void transferValueDF();

    // Get the tensor field corresponding to the domain.
    RCP_TensorField getDomainField() const
    { return d_dof_domain; }

    // Get the tensor field corresponding to the range.
    RCP_TensorField getRangeField() const
    { return d_dof_range; }
};

} // end namespace FOOD

#include "ConsistentScheme_Def.hpp"

#endif // end FOOD_CONSISTENTSCHEME_HPP

//---------------------------------------------------------------------------//
// end ConsistentScheme.hpp
//---------------------------------------------------------------------------//
