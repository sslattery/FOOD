//---------------------------------------------------------------------------//
// \file Octree.cpp
// \author Stuart Slattery
// \breif Octree declaration.
//---------------------------------------------------------------------------//

#include "Octree.hpp"

namespace FOOD
{

/*! 
 * \brief Constructor.
 */
Octree::Octree( RCP_Domain domain )
    : d_domain(domain)
    , d_root_node(0)
{ /* ... */ }

/*!
 * \brief Destructor.
 */
Octree::~Octree()
{ /* ... */ }

} // end namespace food

//---------------------------------------------------------------------------//
// end Octree.cpp
//---------------------------------------------------------------------------//
