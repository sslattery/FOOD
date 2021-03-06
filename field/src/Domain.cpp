//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file Domain.cpp
 * \author Stuart Slattery
 * \brief Domain definition
 */
//---------------------------------------------------------------------------//

#include "Domain.hpp"

namespace FOOD
{
/*!
 * \brief Constructor. 
 *
 * The entities used as support for distribution functions are mesh entities
 * in the given mesh instance, and the entity set is a set of mesh entities
 * from that instance. The set may be the root set.  
 */
Domain::Domain( iMesh_Instance mesh_instance, 
		iBase_EntitySetHandle set_handle,
		const int cn )
    : d_mesh( mesh_instance )
    , d_mesh_set( set_handle )
    , d_cn( cn )
{ /* ... */ }

/*!
 * \brief Destructor.
 */
Domain::~Domain()
{ /* ... */ }

/*!
 * \brief Initialize an iterator over specified entity type and topology.
 *
 * Initialize an array iterator over specified entity type and topology.
 * Iterator returned can be used as input to functions returning entities 
 * for the iterator.  If all entities of a specified type and/or topology
 * are to be iterated, specify iBase_ALL_TYPES or iMesh_ALL_TOPOLOGIES, 
 * respectively.  Specified type or topology must be a value in the
 * iBase_EntityType or iMesh_EntityTopology enumerations, respectively.
 * 
 * \param entity_type Type of entity to iterate. Use enumerations in iMesh
 * definition. 
 * \param enity_toplogy Type of entity toplogy to iterate. Use enumerations in
 * iiMesh_Instance definition.
 * \return Iterator over specified entity type and topology.
 */
int Domain::initEntIter( const int entity_type, 
			 const int entity_topology,
			 EntityIterator iterator )
{
    int return_code;
    iMesh_initEntIter(d_mesh,
		      d_mesh_set,
		      entity_type,
		      entity_topology,
		      0,
		      iterator,
		      &return_code);
    return return_code;
}

/*!
 * \brief Initialize an array iterator over specified entity type, topology,
 *  and size.
 *
 * Initialize an array iterator over specified entity type, topology, and 
 * size, for a specified domain.  Iterator returned can be used 
 * as input to functions returning entities for the iterator.  If all 
 * entities of a specified type and/or topology are to be iterated, 
 * specify iBase_ALL_TYPES or iMesh_ALL_TOPOLOGIES, respectively.  
 * Specified type or topology must be a value in the iBase_EntityType or 
 * iMesh_EntityTopology enumerations, respectively.
 *
 * \param entity_type Type of entity to iterate. Use enumerations in iMesh
 * definition. 
 * \param enity_toplogy Type of entity toplogy to iterate. Use enumerations in
 * iMesh definition.
 * \param array_size Size of blocks for iterator to return.
 * \return Array iterator over specified entity type and topology.
 */
int Domain::initEntArrIter( const int entity_type,
			    const int entity_topology,
			    const int array_size,
			    EntityArrayIterator iterator )
{
    int return_code;
    iMesh_initEntArrIter(d_mesh,
			 d_mesh_set,
			 entity_type,
			 entity_topology,
			 array_size,
			 0,
			 iterator,
			 &return_code);
    return return_code;
}

}

//---------------------------------------------------------------------------//
// end Domain.cpp
//---------------------------------------------------------------------------//
