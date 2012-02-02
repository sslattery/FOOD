//---------------------------------------------------------------------------//
// \file FEMInterpolate.hpp
// \author Stuart Slattery
// \brief Finite element interpolation declaration.
//---------------------------------------------------------------------------//

#ifndef FOOD_FEMINTERPOLATE_HPP
#define FOOD_FEMINTERPOLATE_HPP

#include <Map>

#include "TensorField.hpp"
#include "Octree.hpp"

#include <iBase.h>
#include <iMesh.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include <Intrepid_FieldContainer.hpp>

namespace FOOD
{

template<class Scalar>
class FEMInterpolate
{

  public:

    //@{
    //! Typedefs.
    typedef TensorField<Scalar>                      TensorField_t;
    typedef Teuchos::RCP<TensorField_t>              RCP_TensorField;
    typedef Teuchos::RCP<Octree>                     RCP_Octree;
    typedef Intrepid::FieldContainer<Scalar>         MDArray;
    //@}

  private:

    // Degrees of freedom domain.
    RCP_TensorField d_dof_domain;

    // Degrees of freedom range.
    RCP_TensorField d_dof_range;

    // Octree
    RCP_Octree d_octree;

    // Range to domain mapping.
    std::map<iBase_EntityHandle,iBase_EntityHandle> d_map;

    // Number of vertices in the range.
    int d_num_range_verts;

  public:

    // Constructor.
    FEMInterpolate( RCP_TensorField dof_domain, RCP_TensorField dof_range );

    // Destructor.
    ~FEMInterpolate();

    // Setup for interpolation.
    void setup();

    // Perform interpolation.
    void interpolate();
};

} // end namespace FOOD

#include "FEMInterpolation_Def.hpp"

#endif // end FOOD_FEMINTERPOLATE_HPP

//---------------------------------------------------------------------------//
// end FEMInterpolate.hpp
//---------------------------------------------------------------------------//
