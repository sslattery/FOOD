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

#include <cstddef>

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
    typedef iBase_EntityIterator*                      EntityIterator;
    typedef iBase_EntityArrIterator*                   EntityArrayIterator;
    //@}

  private:

    // iMesh_Instance instance to which the domain set belongs.
    iMesh_Instance d_mesh;

    // iMesh_Instance set over which the field is defined.
    iBase_EntitySetHandle d_mesh_set;

    // Canonical numbering system for this domain (enum).
    std::size_t d_cn;

  public:

    // Constructor.
    Domain( iMesh_Instance mesh_instance, 
	    iBase_EntitySetHandle set_handle,
	    const int cn );

    // Destructor.
    ~Domain();

    //! Get the mesh instance.
    iMesh_Instance getMesh() const
    { return d_mesh; }

    //! Get the mesh set.
    iBase_EntitySetHandle getMeshSet() const
    { return d_mesh_set; }

    //! Get the canonical numbering system.
    int getCN() const
    { return d_cn; }

    // Initialize an iterator over the specified entity type and topology.
    int initEntIter( const int entity_type,
		     const int entity_topology,
		     EntityIterator iterator );

    // Initalize an array iterator over the specified entity type and
    // topology. 
    int initEntArrIter( const int entity_type,
		        const int entity_topology,
		        const int array_size,
		        EntityArrayIterator iterator );
};

}

#endif // end FOOD_DOMAIN_HPP

//---------------------------------------------------------------------------//
// end Domain.hpp
//---------------------------------------------------------------------------//
