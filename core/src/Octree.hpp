//---------------------------------------------------------------------------//
// \file Octree.hpp
// \author Stuart Slattery
// \breif Octree definition.
//---------------------------------------------------------------------------//

#ifndef FOOD_OCTREE_HPP
#define FOOD_OCTREE_HPP

#include "PointQuery.hpp"
#include "Domain.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>

namespace FOOD
{

class Octree
{

  public:

    //@{
    //! Typedefs.
    typedef Teuchos::RCP<Domain>     RCP_Domain;
    //@}

  private:

    // The domain this octree is generated for.
    RCP_Domain d_domain;

  public:

    // Constructor.
    Octree( RCP_Domain domain );

    // Destructor.
    ~Octree();
};

} // end namespace FOOD

#endif // end FOOD_OCTREE_HPP

//---------------------------------------------------------------------------//
// end Octree.hpp
//---------------------------------------------------------------------------//
