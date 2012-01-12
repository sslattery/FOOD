//---------------------------------------------------------------------------//
// \file Domain.cpp
// \author Stuart Slattery
// \brief Domain definition
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
Domain::Domain(iMesh_Instance mesh_instance, iBase_EntitySetHandle set_handle)
    : d_mesh( Teuchos::rcp(mesh_instance) )
    , d_mesh_set( Teuchos::rcp(set_handle) )
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
 * iMesh definition.
 * \return Iterator over specified entity type and topology.
 */
Domain::EntityIterator Domain::initEntIter(const int entity_type, 
					   const int entity_topology)
{
    ErrorCode return_code = 0;
    EntityIterator return_iterator = 0;
    iMesh_initEntIter(d_mesh->ptr(),
		      d_mesh_set->ptr(),
		      entity_type,
		      entity_topology,
		      return_iterator,
		      return_code);
    return return_iterator;
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
Domain::EntityIterator Domain::initEntArrIter(const int entity_type,
					      const int entity_topology,
					      const int array_size)
{
    ErrorCode return_code = 0;
    EntityIterator return_iterator = 0;
    iMesh_initEntArrIter(d_mesh->ptr(),
			 d_mesh_set->ptr(),
			 entity_type,
			 entity_topology,
			 array_size,
			 return_iterator,
			 return_code);
    return return_iterator;
}

}

//---------------------------------------------------------------------------//
// end Domain.cpp
//---------------------------------------------------------------------------//
