//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file Domain.hpp
 * \author Stuart Slattery
 * \brief Domain definition
 */
//---------------------------------------------------------------------------//

#ifndef FOOD_DOMAIN_HPP
#define FOOD_DOMAIN_HPP

#include "Types.hpp"

#include <iBase.h>
#include <iMesh.h>

namespace FOOD
{

class Domain
{
  public:

    //@{
    //! Typedefs.
    typedef iMesh_Instance                             Mesh;
    typedef iBase_EntitySetHandle                      EntitySet;
    typedef iBase_EntityIterator*                      EntityIterator;
    typedef iBase_EntityArrIterator*                   EntityArrayIterator;
    typedef int                                        ErrorCode;
    //@}

  private:

    // Mesh instance to which the domain set belongs.
    Mesh d_mesh;

    // Mesh set over which the field is defined.
    EntitySet d_mesh_set;

  public:

    // Constructor.
    Domain(Mesh mesh_instance, EntitySet set_handle);

    // Destructor.
    ~Domain();

    //! Get the mesh instance.
    Mesh getMesh() const
    { return d_mesh; }

    //! Get the mesh set.
    EntitySet getMeshSet() const
    { return d_mesh_set; }

    // Initialize an iterator over the specified entity type and topology.
    ErrorCode initEntIter(const int entity_type,
			  const int entity_topology,
			  EntityIterator iterator);

    // Initalize an array iterator over the specified entity type and
    // topology. 
    ErrorCode initEntArrIter(const int entity_type,
			     const int entity_topology,
			     const int array_size,
			     EntityArrayIterator iterator);
};

}

#endif // end FOOD_DOMAIN_HPP

//---------------------------------------------------------------------------//
// end Domain.hpp
//---------------------------------------------------------------------------//
