#ifndef _ITAPS_iGeom
#define _ITAPS_iGeom

/***************************************************************************//**
 * \ingroup VersionNumbers
 * \brief Compile time version number digits
 *
 * iMesh maintains a major, minor and patch digit in its version number.
 * Technically speaking, there is not much practical value in patch digit
 * for an interface specification. A patch release is typically only used
 * for bug fix releases. Although it is rare, sometimes a bug fix
 * necessitates an API change. So, we define a patch digit for iMesh.
 ******************************************************************************/
#define IGEOM_VERSION_MAJOR ITAPS_VERSION_MAJOR
#define IGEOM_VERSION_MINOR ITAPS_VERSION_MINOR
#define IGEOM_VERSION_PATCH ITAPS_VERSION_PATCH

/***************************************************************************//**
 * \ingroup VersionNumbers
 * \brief Version Comparison
 *
 * Evaluates to true at CPP time if the version of iMesh currently being
 * compiled is greater than or equal to the version specified.
 ******************************************************************************/
#define IGEOM_VERSION_GE(Maj,Min,Pat) ITAPS_VERSION_GE(Maj,Min,Pat)

/***************************************************************************//**
 * \ingroup VersionNumbers
 * \brief Compose string represention of the iMesh version number
 ******************************************************************************/
#define IGEOM_VERSION_STRING ITAPS_VERSION_STRING_(iGeom)

/***************************************************************************//**
 * \ingroup VersionNumbers
 * \brief Compose a symbol name derived from the current iMesh version number.
 ******************************************************************************/
#define IGEOM_VERSION_TAG ITAPS_VERSION_TAG_(iGeom)

/***************************************************************************//**
 * \ingroup VersionNumbers
 * \brief Define iMesh_newMesh symbol such that it depends on version number.
 *
 * Note: We ran into problems with this as it influences or is influenced by
 * fortran name mangling and so breaks fortran compilation. So, this is
 * currently disabled.
 ******************************************************************************/
#define IGEOM_NEW_GEOM_NAME__(A,B,C) A##_##B##_##C
#define IGEOM_NEW_GEOM_NAME_(A,B,C) IGEOM_NEW_GEOM_NAME__(A,B,C)
#define IGEOM_NEW_GEOM_NAME(A) IGEOM_NEW_GEOM_NAME_(A,IGEOM_VERSION_MAJOR,IGEOM_VERSION_MINOR)

#include "iBase.h"
#include "iGeom_protos.h"
#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************//**
 * \ingroup Datatypes
 * \brief iGeom instance
 ******************************************************************************/
typedef struct iGeom_Instance_Private* iGeom_Instance;

/***************************************************************************//**
 * \ingroup ErrorHandling
 * \brief Get a description of the error returned from the last iGeom function
 *
 * Get a description of the error returned from the last iGeom function
 ******************************************************************************/
void iGeom_getDescription(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    char* descr,
        /**< [in] Pointer to a character string to be filled with a
        description of the error from the last iGeom function */
    int descr_len
        /**< [in] Length of the character string pointed to by descr */
);

/***************************************************************************//**
 * \ingroup ErrorHandling
 * \brief Get the error type returned from the last iGeom function
 *
 * Get the error type returned from the last iGeom function.  Value
 * returned is a member of the iBase_ErrorType enumeration.
 ******************************************************************************/
void iGeom_getErrorType(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    int* error_type
        /**< [out] Error type returned from last iGeom function */
);

/***************************************************************************//**
 * \ingroup iGeomInitialization
 * \brief Construct a new iGeom instance
 *
 * Construct a new iGeom instance, using implementation-specific
 * options
 ******************************************************************************/
void iGeom_newGeom(
    const char* options,
        /**< [in] Pointer to implementation-specific options string */
    iGeom_Instance* instance_out,
        /**< [out] Pointer to iGeom instance handle returned from function */
    int* err,
        /**< [out] Returned Error status (see iBase_ErrorType) */
    int options_len
        /**< [in] Length of the character string pointed to by options */
);

/***************************************************************************//**
 * \ingroup iGeomInitialization
 * \brief Destroy an iGeom instance
 *
 * Destroy an iGeom instance
 ******************************************************************************/
void iGeom_dtor(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomInitialization
 * \brief Load a geom from a file
 *
 * Load a geom from a file.  If entity set is specified, loaded geom
 * is added to that set; specify zero if that is not desired.
 ******************************************************************************/
void iGeom_load(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const char* name,
        /**< [in] File name from which geom is to be loaded */
    const char* options,
        /**< [in] Pointer to implementation-specific options string */
    int* err,
        /**< [out] Returned Error status (see iBase_ErrorType) */
    int name_len,
        /**< [in] Length of the file name character string */
    int options_len
        /**< [in] Length of the options character string \note
        entity_set_handle Set to which loaded geom will be added, zero if not
        desired */
);

/***************************************************************************//**
 * \brief Save a geom to a file
 *
 * Save a geom to a file.  If entity set is specified, save only the
 * geom contained in that set.
 ******************************************************************************/
void iGeom_save(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const char* name,
        /**< [in] File name to which geom is to be saved */
    const char* options,
        /**< [in] Pointer to implementation-specific options string */
    int* err,
        /**< [out] Returned Error status (see iBase_ErrorType) */
    int name_len,
        /**< [in] Length of the file name character string */
    int options_len
        /**< [in] Length of the options character string \note
        entity_set_handle Entity set being saved */
);

/***************************************************************************//**
 * \ingroup iGeomInitialization
 * \brief Get handle of the root set for this instance
 *
 * Get handle of the root set for this instance.  All geom in
 * this instance can be accessed from this set.
 ******************************************************************************/
void iGeom_getRootSet(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle* root_set,
        /**< [out] Pointer to set handle returned from function */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomBounds
 * \brief Get the bounding box of the entire model
 *
 * Get the bounding box of the entire model
 ******************************************************************************/
void iGeom_getBoundBox(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    double* min_x,
        /**< [out] Minimum coordinate of bounding box */
    double* min_y,
        /**< [out] Minimum coordinate of bounding box */
    double* min_z,
        /**< [out] Minimum coordinate of bounding box */
    double* max_x,
        /**< [out] Maximum coordinate of bounding box */
    double* max_y,
        /**< [out] Maximum coordinate of bounding box */
    double* max_z,
        /**< [out] Maximum coordinate of bounding box */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get entities of specific type and/or topology in set or instance
 *
 * Get entities of specific type and/or topology in set or instance.  All 
 * entities of a given type or topology are requested by specifying
 * iBase_ALL_TOPOLOGIES or iBase_ALL_TYPES, respectively.  Specified type
 * or topology must be a value in the iBase_EntityType or iBase_EntityTopology
 * enumeration, respectively.
 ******************************************************************************/
void iGeom_getEntities(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle set_handle,
        /**< [in] Entity set being queried */
    int entity_type,
        /**< [in] Type of entities being requested */
    iBase_EntityHandle** entity_handles,
        /**< [in,out] Pointer to array of entity handles returned from function */
    int* entity_handles_allocated,
        /**< [in,out] Pointer to allocated size of entity_handles array */
    int* entity_handles_size,
        /**< [in,out] Pointer to occupied size of entity_handles array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the number of entities with the specified type in the instance or set
 *
 * Get the number of entities with the specified type in the instance 
 * or set.  If entity set handle is zero, return information for instance,
 * otherwise for set.  Value of entity type must be from the
 * iBase_EntityType enumeration.  If iBase_ALL_TYPES is specified,
 * total number of entities (excluding entity sets) is returned.
 ******************************************************************************/
void iGeom_getNumOfType(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle set_handle,
        /**< [in] Entity set being queried */
    int entity_type,
        /**< [in] Type of entity requested */
    int* num_type,
        /**< [out] Pointer to number of entities, returned from function */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the entity type for the specified entity
 *
 * Get the entity type for the specified entity.  Types
 * returned are values in the iBase_EntityType enumeration.
 ******************************************************************************/
void iGeom_getEntType(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [in] entity handle being queried */
    int* type,
        /**< [out] Pointer to location at which to store the returned type */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the entity type for the specified entities
 *
 * Get the entity type for the specified entities.  Types
 * returned are values in the iBase_EntityType enumeration.
 ******************************************************************************/
void iGeom_getArrType(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [in] Array of entity handles being queried */
    int entity_handles_size,
        /**< [in] Number of entities in entity_handles array */
    int** type,
        /**< [in,out] Pointer to array of types returned from function */
    int* type_allocated,
        /**< [in,out] Pointer to allocated size of type array */
    int* type_size,
        /**< [in,out] Pointer to occupied size of type array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomAdjacencies
 * \brief Get entities of specified type adjacent to an entity
 *
 * Get entities of specified type adjacent to an entity.  Specified type
 * must be value in the iBase_EntityType enumeration.
 ******************************************************************************/
void iGeom_getEntAdj(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [in] Entity handle being queried */
    int to_dimension,
        /**< [in] unknown */
    iBase_EntityHandle** adj_entities,
        /**< [in,out] Pointer to array of adjacent entities returned from function */
    int* adj_entities_allocated,
        /**< [in,out] Pointer to allocated size of adj_entities array */
    int* adj_entities_size,
        /**< [in,out] Pointer to occupied size of adj_entities array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomAdjacencies
 * \brief Get entities of specified type adjacent to entities
 *
 * Get entities of specified type adjacent to entities.  Specified type
 * must be value in the iBase_EntityType enumeration.  \em offset(i) is
 * index of first entity in adjacentEntityHandles array adjacent to 
 * entity_handles[i].
 ******************************************************************************/
void iGeom_getArrAdj(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [in] Array of entity handles being queried */
    int entity_handles_size,
        /**< [in] Number of entities in entity_handles array */
    int requested_entity_type,
        /**< [in] Type of adjacent entities requested */
    iBase_EntityHandle** adj_entity_handles,
        /**< [in,out] Pointer to array of adjacentEntityHandles returned from
        function */
    int* adj_entity_handles_allocated,
        /**< [in,out] Pointer to allocated size of adjacentEntityHandles array */
    int* adj_entity_handles_size,
        /**< [in,out] Pointer to occupied size of adjacentEntityHandles array */
    int** offset,
        /**< [in,out] Pointer to array of offsets returned from function */
    int* offset_allocated,
        /**< [in,out] Pointer to allocated size of offset array */
    int* offset_size,
        /**< [in,out] Pointer to occupied size of offset array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomAdjacencies
 * \brief Get "2nd order" adjacencies to an entity
 *
 * Get "2nd order" adjacencies to an entity, that is, from an entity, through
 * other entities of a specified "bridge" dimension, to other entities of another 
 * specified "to" dimension.
 ******************************************************************************/
void iGeom_getEnt2ndAdj(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [in] Entity from which adjacencies are requested */
    int bridge_dimension,
        /**< [in] Bridge dimension for 2nd order adjacencies */
    int to_dimension,
        /**< [in] Dimension of adjacent entities returned */
    iBase_EntityHandle** adjacent_entities,
        /**< [in,out] Adjacent entities */
    int* adjacent_entities_allocated,
        /**< [in,out] Allocated size of returned array */
    int* adjacent_entities_size,
        /**< [in,out] Occupied size of returned array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomAdjacencies
 * \brief Get "2nd order" adjacencies to an array of entities
 *
 * Get "2nd order" adjacencies to an array of entities, that is, from each entity, through
 * other entities of a specified "bridge" dimension, to other entities of another 
 * specified "to" dimension.
 ******************************************************************************/
void iGeom_getArr2ndAdj(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [in] Entities from which adjacencies are requested */
    int entity_handles_size,
        /**< [in] Number of entities whose adjacencies are requested */
    int order_adjacent_key,
        /**< [in] Bridge dimension for 2nd order adjacencies */
    int requested_entity_type,
        /**< [in] unknown */
    iBase_EntityHandle** adj_entity_handles,
        /**< [in,out] Adjacent entities */
    int* adj_entity_handles_allocated,
        /**< [in,out] Allocated size of returned array */
    int* adj_entity_handles_size,
        /**< [in,out] Occupied size of returned array */
    int** offset,
        /**< [in,out] Offset[i] is offset into adj_entity_handles of 2nd order
        adjacencies of ith entity in entity_handles */
    int* offset_allocated,
        /**< [in,out] Allocated size of offset array */
    int* offset_size,
        /**< [in,out] Occupied size of offset array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomAdjacencies
 * \brief Return whether two entities are adjacent
 *
 * Return whether two entities are adjacent.
 ******************************************************************************/
void iGeom_isEntAdj(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle1,
        /**< [in] First entity queried */
    iBase_EntityHandle entity_handle2,
        /**< [in] Second entity queried */
    int* are_adjacent,
        /**< [out] If returned non-zero, entities are adjacent, otherwise they
        are not */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomAdjacencies
 * \brief Return whether entity pairs are adjacent
 *
 * Return whether entity pairs are adjacent, i.e. if entity_handles_1[i] is
 * adjacent to entity_handles_2[i].  This function requires entity_handles_1_size
 * and entity_handles_2_size to be equal.
 ******************************************************************************/
void iGeom_isArrAdj(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles_1,
        /**< [in] First array of entities */
    int entity_handles_1_size,
        /**< [in] Number of entities in first array */
    const iBase_EntityHandle* entity_handles_2,
        /**< [in] Second array of entities */
    int entity_handles_2_size,
        /**< [in] Number of entities in second array */
    int** is_adjacent_info,
        /**< [in,out] Array of flags returned from function */
    int* is_adjacent_info_allocated,
        /**< [in,out] Allocated size of flags array */
    int* is_adjacent_info_size,
        /**< [in,out] Occupied size of flags array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief DON'T KNOW WHAT THIS FUNCTION IS
 *
 ******************************************************************************/
void iGeom_getTopoLevel(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    int* topo_level_out,
        /**< [out] */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get closest point to an entity
 *
 * Get closest point to a specified position on an entity
 ******************************************************************************/
void iGeom_getEntClosestPt(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [in] Entity being queried */
    double near_x,
        /**< [in] Coordinates of starting point */
    double near_y,
        /**< [in] Coordinates of starting point */
    double near_z,
        /**< [in] Coordinates of starting point */
    double* on_x,
        /**< [out] Closest point on entity */
    double* on_y,
        /**< [out] Closest point on entity */
    double* on_z,
        /**< [out] Closest point on entity */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get closest point for an array of entities and points
 *
 * Get closest point for an array of entities and points.  If either the number
 * of entities or number of coordinate triples is unity, then all points or 
 * entities are queried for that entity or point, respectively, otherwise each
 * point corresponds to each entity.  storage_order should be a value in
 * the iBase_StorageOrder enum; if input as iBase_UNKNOWN, order is in native order
 * with respect to implementation.
 ******************************************************************************/
void iGeom_getArrClosestPt(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [in] Entity(ies) being queried */
    int entity_handles_size,
        /**< [in] Number of entities being queried */
    int storage_order,
        /**< [in] Storage order of input points */
    const double* near_coordinates,
        /**< [in] Coordinates of starting point(s) */
    int near_coordinates_size,
        /**< [in] Number of values in near_coordinates array */
    double** on_coordinates,
        /**< [in,out] Coordinates of closest points */
    int* on_coordinates_allocated,
        /**< [in,out] Allocated size of closest point array */
    int* on_coordinates_size,
        /**< [in,out] Occupied size of closest point array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the normal vector on an entity at the given position
 *
 * Get the normal vector on an entity at the given position.
 ******************************************************************************/
void iGeom_getEntNrmlXYZ(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [in] Entity being queried */
    double x,
        /**< [in] Coordinates of starting point */
    double y,
        /**< [in] Coordinates of starting point */
    double z,
        /**< [in] Coordinates of starting point */
    double* nrml_i,
        /**< [out] Normal vector at starting point */
    double* nrml_j,
        /**< [out] Normal vector at starting point */
    double* nrml_k,
        /**< [out] Normal vector at starting point */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the normal vector on an entity(ies) at given position(s)
 *
 * Get the normal vector on an entity(ies) at given position(s).  If either the 
 * number of entities or number of coordinate triples is unity, then all points or 
 * entities are queried for that entity or point, respectively, otherwise each
 * point corresponds to each entity.  storage_order should be a value in
 * the iBase_StorageOrder enum; if input as iBase_UNKNOWN, order is in native order
 * with respect to implementation.
 ******************************************************************************/
void iGeom_getArrNrmlXYZ(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [in] Entity(ies) being queried */
    int entity_handles_size,
        /**< [in] Number of entities being queried */
    int storage_order,
        /**< [in] Storage order of coordinates */
    const double* coordinates,
        /**< [in] Starting coordinates */
    int coordinates_size,
        /**< [in] Number of values in coordinates array */
    double** normals,
        /**< [in,out] Normal coordinates */
    int* normals_allocated,
        /**< [in,out] Allocated size of normals array */
    int* normals_size,
        /**< [in,out] Occupied size of normals array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the normal vector AND closest point on an entity at given position
 *
 * Get the normal vector AND closest point on an entity at a given position.
 ******************************************************************************/
void iGeom_getEntNrmlPlXYZ(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [in] Entity being queried */
    double x,
        /**< [in] Starting coordinates */
    double y,
        /**< [in] Starting coordinates */
    double z,
        /**< [in] Starting coordinates */
    double* pt_x,
        /**< [out] Closest point */
    double* pt_y,
        /**< [out] Closest point */
    double* pt_z,
        /**< [out] Closest point */
    double* nrml_i,
        /**< [out] Normal at closest point */
    double* nrml_j,
        /**< [out] Normal at closest point */
    double* nrml_k,
        /**< [out] Normal at closest point */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the normal vector AND closest point on an entity(ies) at given position(s)
 *
 * Get the normal vector AND closest point on an entity(ies) at given position(s).  If either the 
 * number of entities or number of coordinate triples is unity, then all points or 
 * entities are queried for that entity or point, respectively, otherwise each
 * point corresponds to each entity.  storage_order should be a value in
 * the iBase_StorageOrder enum; if input as iBase_UNKNOWN, order is in native order
 * with respect to implementation.
 ******************************************************************************/
void iGeom_getArrNrmlPlXYZ(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [in] Entity(ies) being queried */
    int entity_handles_size,
        /**< [in] Number of entity(ies) being queried */
    int storage_order,
        /**< [in] Storage order in near_coordinates array */
    const double* near_coordinates,
        /**< [in] Starting coordinates */
    int near_coordinates_size,
        /**< [in] Number of values in near_coordinates array */
    double** on_coordinates,
        /**< [in,out] Closest point array */
    int* on_coordinates_allocated,
        /**< [in,out] Allocated size of closest point array */
    int* on_coordinates_size,
        /**< [in,out] Occupied size of closest point array */
    double** normals,
        /**< [in,out] Normal array */
    int* normals_allocated,
        /**< [in,out] Allocated size of normal array */
    int* normals_size,
        /**< [in,out] Occupied size of normal array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the tangent vector on an entity at given position
 *
 * Get the tangent vector on an entity at a given position.
 ******************************************************************************/
void iGeom_getEntTgntXYZ(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [in] Entity being queried */
    double x,
        /**< [in] Starting coordinates */
    double y,
        /**< [in] Starting coordinates */
    double z,
        /**< [in] Starting coordinates */
    double* tgnt_i,
        /**< [out] Tangent at closest point */
    double* tgnt_j,
        /**< [out] Tangent at closest point */
    double* tgnt_k,
        /**< [out] Tangent at closest point */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the tangent vector on an entity(ies) at given position(s)
 *
 * Get the tangent vector on an entity(ies) at given position(s).  If either the 
 * number of entities or number of coordinate triples is unity, then all points or 
 * entities are queried for that entity or point, respectively, otherwise each
 * point corresponds to each entity.  storage_order should be a value in
 * the iBase_StorageOrder enum; if input as iBase_UNKNOWN, order is in native order
 * with respect to implementation.
 ******************************************************************************/
void iGeom_getArrTgntXYZ(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [in] Entity(ies) being queried */
    int entity_handles_size,
        /**< [in] Number of entities being queried */
    int storage_order,
        /**< [in] Storage order of coordinates */
    const double* coordinates,
        /**< [in] Starting coordinates */
    int coordinates_size,
        /**< [in] Number of values in coordinates array */
    double** tangents,
        /**< [in,out] Tangent coordinates */
    int* tangents_allocated,
        /**< [in,out] Allocated size of tangents array */
    int* tangents_size,
        /**< [in,out] Occupied size of tangents array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomCurvature
 * \brief Get the two principle curvature vectors for a face at a point
 *
 * Get the two principle curvature vectors for a face at a point.  Magnitudes of
 * vectors are curvature, directions are directions of principal curvatures.
 ******************************************************************************/
void iGeom_getFcCvtrXYZ(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle face_handle,
        /**< [in] Face being queried */
    double x,
        /**< [in] Position being queried */
    double y,
        /**< [in] Position being queried */
    double z,
        /**< [in] Position being queried */
    double* cvtr1_i,
        /**< [out] Maximum curvature vector */
    double* cvtr1_j,
        /**< [out] Maximum curvature vector */
    double* cvtr1_k,
        /**< [out] Maximum curvature vector */
    double* cvtr2_i,
        /**< [out] Minimum curvature vector */
    double* cvtr2_j,
        /**< [out] Minimum curvature vector */
    double* cvtr2_k,
        /**< [out] Minimum curvature vector */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomCurvature
 * \brief Get the principle curvature vector for an edge at a point
 *
 * Get the principle curvature vector for an edge at a point.  Magnitude of
 * vector is the curvature, direction is direction of principal curvature.
 ******************************************************************************/
void iGeom_getEgCvtrXYZ(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle edge_handle,
        /**< [in] Edge being queried */
    double x,
        /**< [in] Position being queried */
    double y,
        /**< [in] Position being queried */
    double z,
        /**< [in] Position being queried */
    double* cvtr_i,
        /**< [out] Maximum curvature vector */
    double* cvtr_j,
        /**< [out] Maximum curvature vector */
    double* cvtr_k,
        /**< [out] Maximum curvature vector */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomCurvature
 * \brief Get the curvature(s) on an entity(ies) at given position(s)
 *
 * Get the curvature(s) on an entity(ies) at given position(s).  If either the 
 * number of entities or number of coordinate triples is unity, then all points or 
 * entities are queried for that entity or point, respectively, otherwise each
 * point corresponds to each entity.  storage_order should be a value in
 * the iBase_StorageOrder enum; if input as iBase_UNKNOWN, order is in native order
 * with respect to implementation.
 ******************************************************************************/
void iGeom_getEntArrCvtrXYZ(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [in] Entity(ies) being queried */
    int entity_handles_size,
        /**< [in] Number of entities being queried */
    int storage_order,
        /**< [in] Storage order of coordinates */
    const double* coords,
        /**< [in] Starting coordinates */
    int coords_size,
        /**< [in] Number of values in coordinates array */
    double** cvtr_1,
        /**< [in,out] First principal curvatures */
    int* cvtr_1_allocated,
        /**< [in,out] Allocated size of first curvature array */
    int* cvtr_1_size,
        /**< [in,out] Occupied size of first curvature array */
    double** cvtr_2,
        /**< [in,out] Second principal curvatures */
    int* cvtr_2_allocated,
        /**< [in,out] Allocated size of second curvature array */
    int* cvtr_2_size,
        /**< [in,out] Occupied size of second curvature array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get closest point, tangent, and curvature of edge
 *
 * Get closest point, tangent, and curvature of edge.
 ******************************************************************************/
void iGeom_getEgEvalXYZ(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle edge_handle,
        /**< [in] Edge being queried */
    double x,
        /**< [in] Point at which entity is being queried */
    double y,
        /**< [in] Point at which entity is being queried */
    double z,
        /**< [in] Point at which entity is being queried */
    double* on_x,
        /**< [out] Closest point at point being queried */
    double* on_y,
        /**< [out] Closest point at point being queried */
    double* on_z,
        /**< [out] Closest point at point being queried */
    double* tgnt_i,
        /**< [out] Tangent at point being queried */
    double* tgnt_j,
        /**< [out] Tangent at point being queried */
    double* tgnt_k,
        /**< [out] Tangent at point being queried */
    double* cvtr_i,
        /**< [out] Curvature at point being queried */
    double* cvtr_j,
        /**< [out] Curvature at point being queried */
    double* cvtr_k,
        /**< [out] Curvature at point being queried */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get closest point, tangent, and curvature of face
 *
 * Get closest point, tangent, and curvature of face.  If any of input
 * coordinate pointers are NULL, that value is not returned.
 ******************************************************************************/
void iGeom_getFcEvalXYZ(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle face_handle,
        /**< [in] Face being queried */
    double x,
        /**< [in] Point at which entity is being queried */
    double y,
        /**< [in] Point at which entity is being queried */
    double z,
        /**< [in] Point at which entity is being queried */
    double* on_x,
        /**< [out] Closest point at point being queried */
    double* on_y,
        /**< [out] Closest point at point being queried */
    double* on_z,
        /**< [out] Closest point at point being queried */
    double* nrml_i,
        /**< [out] Normal at point being queried */
    double* nrml_j,
        /**< [out] Normal at point being queried */
    double* nrml_k,
        /**< [out] Normal at point being queried */
    double* cvtr1_i,
        /**< [out] First principal curvature at point being queried */
    double* cvtr1_j,
        /**< [out] First principal curvature at point being queried */
    double* cvtr1_k,
        /**< [out] First principal curvature at point being queried */
    double* cvtr2_i,
        /**< [out] Second principal curvature at point being queried */
    double* cvtr2_j,
        /**< [out] Second principal curvature at point being queried */
    double* cvtr2_k,
        /**< [out] Second principal curvature at point being queried */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the closest point(s), tangent(s), and curvature(s) on an entity(ies) at given position(s)
 *
 * Get the closest point(s), tangent(s), and curvature(s) on an entity(ies) at given position(s).  
 * If either the number of entities or number of coordinate triples is unity, then all points or 
 * entities are queried for that entity or point, respectively, otherwise each
 * point corresponds to each entity.  storage_order should be a value in
 * the iBase_StorageOrder enum; if input as iBase_UNKNOWN, order is in native order
 * with respect to implementation.
 ******************************************************************************/
void iGeom_getArrEgEvalXYZ(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* edge_handles,
        /**< [???] Edge(s) being queried */
    int edge_handles_size,
        /**< [???] Number of edges being queried */
    int storage_order,
        /**< [???] Storage order of coordinates */
    const double* coords,
        /**< [???] Starting coordinates */
    int coords_size,
        /**< [???] Number of values in coordinates array */
    double** on_coords,
        /**< [???] Closest point array */
    int* on_coords_allocated,
        /**< [???] Allocated size of closest point array */
    int* on_coords_size,
        /**< [???] Occupied size of closest point array */
    double** tangent,
        /**< [???] Tangent array */
    int* tangent_allocated,
        /**< [???] Allocated size of tangent array */
    int* tangent_size,
        /**< [???] Occupied size of tangent array */
    double** cvtr,
        /**< [???] First principal curvatures */
    int* cvtr_allocated,
        /**< [???] Allocated size of first curvature array */
    int* cvtr_size,
        /**< [???] Occupied size of first curvature array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the closest point(s), tangent(s), and curvature(s) on an entity(ies) at given position(s)
 *
 * Get the closest point(s), tangent(s), and curvature(s) on an entity(ies) at given position(s).  
 * If either the number of entities or number of coordinate triples is unity, then all points or 
 * entities are queried for that entity or point, respectively, otherwise each
 * point corresponds to each entity.  storage_order should be a value in
 * the iBase_StorageOrder enum; if input as iBase_UNKNOWN, order is in native order
 * with respect to implementation.
 ******************************************************************************/
void iGeom_getArrFcEvalXYZ(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* face_handles,
        /**< [???] Face(s) being queried */
    int face_handles_size,
        /**< [???] Number of edges being queried */
    int storage_order,
        /**< [???] Storage order of coordinates */
    const double* coords,
        /**< [???] Starting coordinates */
    int coords_size,
        /**< [???] Number of values in coordinates array */
    double** on_coords,
        /**< [???] Closest point array */
    int* on_coords_allocated,
        /**< [???] Allocated size of closest point array */
    int* on_coords_size,
        /**< [???] Occupied size of closest point array */
    double** normal,
        /**< [???] Normal array */
    int* normal_allocated,
        /**< [???] Allocated size of normal array */
    int* normal_size,
        /**< [???] Occupied size of normal array */
    double** cvtr1,
        /**< [???] description unknown */
    int* cvtr1_allocated,
        /**< [???] Allocated size of first curvature array */
    int* cvtr1_size,
        /**< [???] Occupied size of first curvature array */
    double** cvtr2,
        /**< [???] Second principal curvatures */
    int* cvtr2_allocated,
        /**< [???] Allocated size of second curvature array */
    int* cvtr2_size,
        /**< [???] Occupied size of second curvature array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomBounds
 * \brief Get the bounding box of the specified entity
 *
 * Get the bounding box of the specified entity
 ******************************************************************************/
void iGeom_getEntBoundBox(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [???] Entity being queried */
    double* min_x,
        /**< [???] Minimum coordinate of bounding box */
    double* min_y,
        /**< [???] Minimum coordinate of bounding box */
    double* min_z,
        /**< [???] Minimum coordinate of bounding box */
    double* max_x,
        /**< [???] Maximum coordinate of bounding box */
    double* max_y,
        /**< [???] Maximum coordinate of bounding box */
    double* max_z,
        /**< [???] Maximum coordinate of bounding box */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomBounds
 * \brief Get the bounding box of the specified entities
 *
 * Get the bounding box of the specified entities.  Storage order passed back
 * will be member of iBase_StorageOrder enum.
 ******************************************************************************/
void iGeom_getArrBoundBox(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [???] Entity handles being queried */
    int entity_handles_size,
        /**< [???] Number of entities being queried */
    int storage_order,
        /**< [???] Storage order of coordinates passed back */
    double** min_corner,
        /**< [???] Minimum coordinates of bounding boxes */
    int* min_corner_allocated,
        /**< [???] Allocated size of minimum coordinates array */
    int* min_corner_size,
        /**< [???] Occupied size of minimum coordinates array */
    double** max_corner,
        /**< [???] Maximum coordinates of bounding boxes */
    int* max_corner_allocated,
        /**< [???] Allocated size of maximum coordinates array */
    int* max_corner_size,
        /**< [???] Occupied size of maximum coordinates array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get coordinates of specified vertex
 *
 * Get coordinates of specified vertex.
 ******************************************************************************/
void iGeom_getVtxCoord(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle vertex_handle,
        /**< [???] Geom vertex being queried */
    double* x,
        /**< [???] Pointer to x coordinate returned from function */
    double* y,
        /**< [???] Pointer to y coordinate returned from function */
    double* z,
        /**< [???] Pointer to z coordinate returned from function */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get coordinates of specified vertices
 *
 * Get coordinates of specified vertices.  If storage order is passed in
 * with a value other than iBase_UNDETERMINED, coordinates are returned
 * in the specified storage order, otherwise storage order is that native
 * to the implementation.  Storage order of returned coordinates is also
 * returned.
 ******************************************************************************/
void iGeom_getVtxArrCoords(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* vertex_handles,
        /**< [???] Array of geom vertex handles whose coordinates are being
        requested */
    int vertex_handles_size,
        /**< [???] Number of vertices in vertex_handles array */
    int storage_order,
        /**< [???] Storage order requested for coordinate data */
    double** coords,
        /**< [???] Pointer to array of coordinates returned from function */
    int* coords_allocated,
        /**< [???] Pointer to allocated size of coords array */
    int* coords_size,
        /**< [???] Pointer to occupied size of coords array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Intersect a ray with the model
 *
 * Intersect a ray with the model.  Storage orders passed back are members of the
 * iBase_StorageOrder enumeration; if output is iBase_UNKNOWN, order is in native order
 * with respect to implementation.
 ******************************************************************************/
void iGeom_getPntRayIntsct(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    double x,
        /**< [???] Point from which ray is fired */
    double y,
        /**< [???] Point from which ray is fired */
    double z,
        /**< [???] Point from which ray is fired */
    double dir_x,
        /**< [???] Direction in which ray is fired */
    double dir_y,
        /**< [???] Direction in which ray is fired */
    double dir_z,
        /**< [???] Direction in which ray is fired */
    iBase_EntityHandle** intersect_entity_handles,
        /**< [???] Entities intersected by ray */
    int* intersect_entity_handles_allocated,
        /**< [???] Allocated size of intersections array */
    int* intersect_entity_handles_size,
        /**< [???] Occupied size of intersections array */
    int storage_order,
        /**< [???] Storage order of coordinates passed back */
    double** intersect_coords,
        /**< [???] Coordinates of intersections */
    int* intersect_coords_allocated,
        /**< [???] Allocated size of coordinates array */
    int* intersect_coords_size,
        /**< [???] Occupied size of coordinates array */
    double** param_coords,
        /**< [???] Distances along ray of intersections */
    int* param_coords_allocated,
        /**< [???] Allocated size of param_coords array */
    int* param_coords_size,
        /**< [???] Occupied size of param_coords array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Intersect an array of rays with the model
 *
 * Intersect an array of rays with the model.  Storage orders passed back are members of the
 * iBase_StorageOrder enumeration; if input/output is iBase_UNKNOWN, order is in native order
 * with respect to implementation.
 ******************************************************************************/
void iGeom_getPntArrRayIntsct(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    int storage_order,
        /**< [???] Storage order of input coordinates */
    const double* coords,
        /**< [???] Points from which rays are fired */
    int coords_size,
        /**< [???] Number of points from which rays are fired */
    const double* directions,
        /**< [???] Directions in which rays are fired */
    int directions_size,
        /**< [???] Number of coordinates in directions array */
    iBase_EntityHandle** intersect_entity_handles,
        /**< [???] Entities intersected by ray */
    int* intersect_entity_handles_allocated,
        /**< [???] Allocated size of intersections array */
    int* intersect_entity_handles_size,
        /**< [???] Occupied size of intersections array */
    int** offset,
        /**< [???] Offset[i] is offset into intersect_entity_handles of ith ray */
    int* offset_allocated,
        /**< [???] Allocated size of offset array */
    int* offset_size,
        /**< [???] Occupied size of offset array */
    double** intersect_coords,
        /**< [???] Coordinates of intersections */
    int* intersect_coords_allocated,
        /**< [???] Allocated size of coordinates array */
    int* intersect_coords_size,
        /**< [???] Occupied size of coordinates array */
    double** param_coords,
        /**< [???] Distances along ray of intersections */
    int* param_coords_allocated,
        /**< [???] Allocated size of param_coords array */
    int* param_coords_size,
        /**< [???] Occupied size of param_coords array \note storage_order
        Storage order of coordinates passed back */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the entity on which a point is located
 *
 * Get the entity on which a point is located
 ******************************************************************************/
void iGeom_getPntClsf(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    double x,
        /**< [???] Point being queried */
    double y,
        /**< [???] Point being queried */
    double z,
        /**< [???] Point being queried */
    iBase_EntityHandle* entity_handle,
        /**< [???] Entity on which point is located */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the entities on which points are located
 *
 * Get the entities on which points are located.   Storage orders should be members of the
 * iBase_StorageOrder enumeration; if input is iBase_UNKNOWN, order is in native order
 * with respect to implementation.
 ******************************************************************************/
void iGeom_getPntArrClsf(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    int storage_order,
        /**< [???] Storage order of coordinates in coords */
    const double* coords,
        /**< [???] Points being queried */
    int coords_size,
        /**< [???] Number of entries in coords array */
    iBase_EntityHandle** entity_handles,
        /**< [???] Entities on which points are located */
    int* entity_handles_allocated,
        /**< [???] Allocated size of entity_handles array */
    int* entity_handles_size,
        /**< [???] Occupied size of entity_handles array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the sense of a face with respect to a region
 *
 * Get the sense of a face with respect to a region.  Sense returned is -1, 0, or 1,
 * representing "reversed", "both", or "forward".  "both" sense indicates that face bounds
 * the region once with each sense.
 ******************************************************************************/
void iGeom_getEntNrmlSense(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle face,
        /**< [???] Face being queried */
    iBase_EntityHandle region,
        /**< [???] Region being queried */
    int* sense_out,
        /**< [???] Sense of face with respect to region */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the senses of an array of faces with respect to an array of regions
 *
 * Get the senses of an array of faces with respect to an array of regions.  Sense returned 
 * is -1, 0, or 1, representing "reversed", "both", or "forward".  "both" sense indicates 
 * that face bounds the region once with each sense.
 ******************************************************************************/
void iGeom_getArrNrmlSense(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* face_handles,
        /**< [???] Faces being queried */
    int face_handles_size,
        /**< [???] Size of face handles array */
    const iBase_EntityHandle* region_handles,
        /**< [???] Regions being queried */
    int region_handles_size,
        /**< [???] Size of region handles array */
    int** sense,
        /**< [???] Senses of faces with respect to regions */
    int* sense_allocated,
        /**< [???] Allocated size of senses array */
    int* sense_size,
        /**< [???] Occupied size of senses array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the sense of an edge with respect to a face
 *
 * Get the sense of an edge with respect to a face.  Sense returned is -1, 0, or 1,
 * representing "reversed", "both", or "forward".  "both" sense indicates that edge bounds
 * the face once with each sense.
 ******************************************************************************/
void iGeom_getEgFcSense(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle edge,
        /**< [???] Edge being queried */
    iBase_EntityHandle face,
        /**< [???] Face being queried */
    int* sense_out,
        /**< [???] Sense of edge with respect to face */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the senses of an array of edges with respect to an array of faces
 *
 * Get the senses of an array of edges with respect to an array of faces.  Sense returned 
 * is -1, 0, or 1, representing "reversed", "both", or "forward".  "both" sense indicates 
 * that edge bounds the face once with each sense.
 ******************************************************************************/
void iGeom_getEgFcArrSense(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* edge_handles,
        /**< [???] Edges being queried */
    int edge_handles_size,
        /**< [???] Size of edge handles array */
    const iBase_EntityHandle* face_handles,
        /**< [???] Faces being queried */
    int face_handles_size,
        /**< [???] Size of face handles array */
    int** sense,
        /**< [???] Senses of faces with respect to regions */
    int* sense_allocated,
        /**< [???] Allocated size of senses array */
    int* sense_size,
        /**< [???] Occupied size of senses array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the sense of a vertex pair with respect to an edge
 *
 * Get the sense of a vertex pair with respect to an edge.  Sense returned is -1, 0, or 1,
 * representing "reversed", "both", or "forward".  "both" sense indicates that vertices
 * are identical and that vertex bounds both sides of the edge
 ******************************************************************************/
void iGeom_getEgVtxSense(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle edge,
        /**< [???] Edge being queried */
    iBase_EntityHandle vertex1,
        /**< [???] First vertex being queried */
    iBase_EntityHandle vertex2,
        /**< [???] Second vertex being queried */
    int* sense_out,
        /**< [???] Sense of vertex pair with respect to edge */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the senses of vertex pair with respect to a edges
 *
 * Get the senses of vertex pairs with respect to edges.  Sense returned is -1, 0, or 1,
 * representing "reversed", "both", or "forward".  "both" sense indicates that both vertices
 * in pair are identical and that vertex bounds both sides of the edge
 ******************************************************************************/
void iGeom_getEgVtxArrSense(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* edge_handles,
        /**< [???] Edges being queried */
    int edge_handles_size,
        /**< [???] Number of edges being queried */
    const iBase_EntityHandle* vertex_handles_1,
        /**< [???] First vertex being queried */
    int vertex_handles_1_size,
        /**< [???] Number of vertices in vertices array */
    const iBase_EntityHandle* vertex_handles_2,
        /**< [???] Second vertex being queried */
    int vertex_handles_2_size,
        /**< [???] Number of vertices in vertices array */
    int** sense,
        /**< [???] Sense of vertex pair with respect to edge */
    int* sense_allocated,
        /**< [???] Allocated size of sense array */
    int* sense_size,
        /**< [???] Occupied size of sense array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Return the measure (length, area, volume) of entities
 *
 * Return the measure (length, area, volume) of entities
 ******************************************************************************/
void iGeom_measure(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [???] Array of entities being queried */
    int entity_handles_size,
        /**< [???] Number of entities in entity array */
    double** measures,
        /**< [???] Measures of entities being queried */
    int* measures_allocated,
        /**< [???] Allocated size of measures array */
    int* measures_size,
        /**< [???] Occupied size of measures array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the geometric type of an entity
 *
 * Get the geometric type of an entity.  Specific types depend on implementation.
 ******************************************************************************/
void iGeom_getFaceType(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle face_handle,
        /**< [???] Face being queried */
    char* face_type,
        /**< [???] Face type */
    int* err,
        /**< [out] Returned Error status (see iBase_ErrorType) */
    int* face_type_length
        /**< [???] Length of face type string */
);

/***************************************************************************//**
 * \brief Return whether interface has information about parameterization
 *
 * Return whether an interface has information about parameterization (!=0) or not (0)
 ******************************************************************************/
void iGeom_getParametric(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    int* is_parametric,
        /**< [???] If non-zero, interface has information about parameterization */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Return whether an entity has a parameterization
 *
 * Return whether an entity has a parameterization (!= 0) or not (=0)
 ******************************************************************************/
void iGeom_isEntParametric(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [???] Entity being queried */
    int* is_parametric,
        /**< [???] Entity has a parameterization (!= 0) or not (=0) */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Return whether entities have parameterizations
 *
 * Return whether entities have parameterizations (!= 0) or not (=0)
 ******************************************************************************/
void iGeom_isArrParametric(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [???] Entities being queried */
    int entity_handles_size,
        /**< [???] Number of entities being queried */
    int** is_parametric,
        /**< [???] entity_handles[i] has a parameterization (!= 0) or not (=0) */
    int* is_parametric_allocated,
        /**< [???] Allocated size of is_parametric array */
    int* is_parametric_size,
        /**< [???] Occupied size of is_parametric array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomCoordEval 
 * \brief Return coordinate position at specified parametric position on entity
 *
 * Return coordinate position at specified parametric position on entity.
 ******************************************************************************/
void iGeom_getEntUVtoXYZ(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [???] Entity being queried */
    double u,
        /**< [???] Parametric coordinate being queried */
    double v,
        /**< [???] Parametric coordinate being queried */
    double* x,
        /**< [???] Spatial coordinate at parametric position being queried */
    double* y,
        /**< [???] Spatial coordinate at parametric position being queried */
    double* z,
        /**< [???] Spatial coordinate at parametric position being queried */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomCoordEval 
 * \brief Return coordinate positions at specified parametric position(s) on entity(ies)
 *
 * Return coordinate positions at specified parametric position(s) on entity(ies).
 * If either the number of entities or number of parametric coordinate pairs is unity, then 
 * all points or entities are queried for that entity or point, respectively, otherwise each
 * point corresponds to each entity.  storage_order should be a value in
 * the iBase_StorageOrder enum; if input as iBase_UNKNOWN, order is in native order
 * with respect to implementation.
 ******************************************************************************/
void iGeom_getArrUVtoXYZ(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [???] Entities being queried */
    int entity_handles_size,
        /**< [???] Number of entities being queried */
    int storage_order,
        /**< [???] Storage order of uv coordinates input and xyz coordinate output */
    const double* uv,
        /**< [???] Coordinates being queried */
    int uv_size,
        /**< [???] Number of coordinates in array */
    double** coordinates,
        /**< [???] Coordinates of parametric positions */
    int* coordinates_allocated,
        /**< [???] Allocated size of coordinates array */
    int* coordinates_size,
        /**< [???] Occupied size of coordinates array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomCoordEval 
 * \brief Return coordinate position at specified parametric position on entity
 *
 * Return coordinate position at specified parametric position on entity.
 ******************************************************************************/
void iGeom_getEntUtoXYZ(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [???] Entity being queried */
    double u,
        /**< [???] Parametric coordinate being queried */
    double* x,
        /**< [???] Spatial coordinate at parametric position being queried */
    double* y,
        /**< [???] Spatial coordinate at parametric position being queried */
    double* z,
        /**< [???] Spatial coordinate at parametric position being queried */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomCoordEval 
 * \brief Return coordinate positions at specified parametric position(s) on entity(ies)
 *
 * Return coordinate positions at specified parametric position(s) on entity(ies).
 * If either the number of entities or number of parametric coordinate pairs is unity, then 
 * all points or entities are queried for that entity or point, respectively, otherwise each
 * point corresponds to each entity.  storage_order should be a value in
 * the iBase_StorageOrder enum; if input as iBase_UNKNOWN, order is in native order
 * with respect to implementation.
 ******************************************************************************/
void iGeom_getArrUtoXYZ(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [???] Entities being queried */
    int entity_handles_size,
        /**< [???] Number of entities being queried */
    const double* u,
        /**< [???] Coordinates being queried */
    int u_size,
        /**< [???] Number of coordinates in array */
    int storage_order,
        /**< [???] Storage order of resulting coordinates */
    double** on_coords,
        /**< [???] Coordinates of parametric positions */
    int* on_coords_allocated,
        /**< [???] Allocated size of coordinates array */
    int* on_coords_size,
        /**< [???] Occupied size of coordinates array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomCoordEval 
 * \brief Return parametric position at specified spatial position on entity
 *
 * Return parametric position at specified spatial position on entity
 ******************************************************************************/
void iGeom_getEntXYZtoUV(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [???] Entity being queried */
    double x,
        /**< [???] Spatial coordinate being queried */
    double y,
        /**< [???] Spatial coordinate being queried */
    double z,
        /**< [???] Spatial coordinate being queried */
    double* u,
        /**< [???] Parametric coordinate at spatial position being queried */
    double* v,
        /**< [???] Parametric coordinate at spatial position being queried */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomCoordEval 
 * \brief Return parametric position at specified spatial position on entity
 *
 * Return parametric position at specified spatial position on entity
 ******************************************************************************/
void iGeom_getEntXYZtoU(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [???] Entity being queried */
    double x,
        /**< [???] Spatial coordinate being queried */
    double y,
        /**< [???] Spatial coordinate being queried */
    double z,
        /**< [???] Spatial coordinate being queried */
    double* u,
        /**< [???] Parametric coordinate at spatial position being queried */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomCoordEval 
 * \brief Return parametric positions at specified spatial position(s) on entity(ies)
 *
 * Return parametric positions at specified spatial position(s) on entity(ies).
 * If either the number of entities or number of spatial coordinate triples is unity, then 
 * all points or entities are queried for that entity or point, respectively, otherwise each
 * point corresponds to each entity.  storage_order should be a value in
 * the iBase_StorageOrder enum; if input as iBase_UNKNOWN, order is in native order
 * with respect to implementation.
 ******************************************************************************/
void iGeom_getArrXYZtoUV(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [???] Entities being queried */
    int entity_handles_size,
        /**< [???] Number of entities being queried */
    int storage_order,
        /**< [???] Storage order of spatial coordinates input */
    const double* coordinates,
        /**< [???] Coordinates being queried */
    int coordinates_size,
        /**< [???] Number of coordinates in array */
    double** uv,
        /**< [???] Coordinates of parametric positions */
    int* uv_allocated,
        /**< [???] Allocated size of coordinates array */
    int* uv_size,
        /**< [???] Occupied size of coordinates array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomCoordEval 
 * \brief Return spatial positions at specified parametric position(s) on entity(ies)
 *
 * Return spatial positions at specified parametric position(s) on entity(ies).
 * If either the number of entities or number of spatial coordinate triples is unity, then 
 * all points or entities are queried for that entity or point, respectively, otherwise each
 * point corresponds to each entity.  storage_order should be a value in
 * the iBase_StorageOrder enum; if input as iBase_UNKNOWN, order is in native order
 * with respect to implementation.
 ******************************************************************************/
void iGeom_getArrXYZtoU(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [???] Entities being queried */
    int entity_handles_size,
        /**< [???] Number of entities being queried */
    int storage_order,
        /**< [???] Storage order of spatial coordinates input */
    const double* coordinates,
        /**< [???] Coordinates being queried */
    int coordinates_size,
        /**< [???] Number of coordinates in array */
    double** u,
        /**< [???] Coordinates of parametric positions */
    int* u_allocated,
        /**< [???] Allocated size of coordinates array */
    int* u_size,
        /**< [???] Occupied size of coordinates array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomCoordEval 
 * \brief Return parametric position at specified spatial position on entity, based on parametric position hint
 *
 * Return parametric position at specified spatial position on entity, based on parametric
 * position hint.  For this function, u and v are input with parameters from which to start search.
 * Typically this will reduce the search time for new parametric coordinates.
 ******************************************************************************/
void iGeom_getEntXYZtoUVHint(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [???] Entity being queried */
    double x,
        /**< [???] Spatial coordinate being queried */
    double y,
        /**< [???] Spatial coordinate being queried */
    double z,
        /**< [???] Spatial coordinate being queried */
    double* u,
        /**< [???] Parametric coordinate at spatial position being queried */
    double* v,
        /**< [???] Parametric coordinate at spatial position being queried */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomCoordEval 
 * \brief Return parametric positions at specified spatial position(s) on entity(ies), based on parametric position hints
 *
 * Return parametric positions at specified spatial position(s) on entity(ies), based on 
 * parametric position hints.  If either the number of entities or number of spatial 
 * coordinate triples is unity, then all points or entities are queried for that entity 
 * or point, respectively, otherwise each point corresponds to each entity.  storage_order 
 * should be a value in the iBase_StorageOrder enum; if input as iBase_UNKNOWN, order is in 
 * native order with respect to implementation.
 ******************************************************************************/
void iGeom_getArrXYZtoUVHint(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [???] Entities being queried */
    int entity_handles_size,
        /**< [???] Number of entities being queried */
    int storage_order,
        /**< [???] Storage order of spatial coordinates input */
    const double* coords,
        /**< [???] Coordinates being queried */
    int coords_size,
        /**< [???] Number of coordinates in array */
    double** uv,
        /**< [???] Coordinates of parametric positions */
    int* uv_allocated,
        /**< [???] Allocated size of coordinates array */
    int* uv_size,
        /**< [???] Occupied size of coordinates array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomCoordEval 
 * \brief Get parametric range of entity
 *
 * Get parametric range of entity
 ******************************************************************************/
void iGeom_getEntUVRange(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [???] Entity being queried */
    double* u_min,
        /**< [???] Minimum parametric coordinate for entity */
    double* v_min,
        /**< [???] Minimum parametric coordinate for entity */
    double* u_max,
        /**< [???] Maximum parametric coordinate for entity */
    double* v_max,
        /**< [???] Maximum parametric coordinate for entity */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomCoordEval 
 * \brief Get parametric range of entity
 *
 * Get parametric range of entity
 ******************************************************************************/
void iGeom_getEntURange(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [???] Entity being queried */
    double* u_min,
        /**< [???] Minimum parametric coordinate for entity */
    double* u_max,
        /**< [???] Maximum parametric coordinate for entity */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomCoordEval 
 * \brief Get parametric range of entities
 *
 * Get parametric range of entities
 ******************************************************************************/
void iGeom_getArrUVRange(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [???] Entities being queried */
    int entity_handles_size,
        /**< [???] Number of entities being queried */
    int storage_order,
        /**< [???] Storage order of parametric coordinates being returned */
    double** uv_min,
        /**< [???] Minimum parametric coordinate for entities */
    int* uv_min_allocated,
        /**< [???] Allocated size of minimum parametric coordinate array */
    int* uv_min_size,
        /**< [???] Occupied size of minimum parametric coordinate array */
    double** uv_max,
        /**< [???] Maximum parametric coordinate for entities */
    int* uv_max_allocated,
        /**< [???] Allocated size of maximum parametric coordinate array */
    int* uv_max_size,
        /**< [???] Occupied size of maximum parametric coordinate array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomCoordEval 
 * \brief Get parametric range of entities
 *
 * Get parametric range of entities
 ******************************************************************************/
void iGeom_getArrURange(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [???] Entities being queried */
    int entity_handles_size,
        /**< [???] Number of entities being queried */
    double** u_min,
        /**< [???] Minimum parametric coordinate for entities */
    int* u_min_allocated,
        /**< [???] Allocated size of minimum parametric coordinate array */
    int* u_min_size,
        /**< [???] Occupied size of minimum parametric coordinate array */
    double** u_max,
        /**< [???] Maximum parametric coordinate for entities */
    int* u_max_allocated,
        /**< [???] Allocated size of maximum parametric coordinate array */
    int* u_max_size,
        /**< [???] Occupied size of maximum parametric coordinate array \note
        storage_order Storage order of parametric coordinates being returned */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomCoordEval 
 * \brief Return the face parametric coordinates for a parametric position on a bounding edge 
 *
 * Return the face parametric coordinates for a parametric position on a bounding edge
 ******************************************************************************/
void iGeom_getEntUtoUV(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle edge_handle,
        /**< [???] Edge being queried */
    iBase_EntityHandle face_handle,
        /**< [???] Face being queried */
    double in_u,
        /**< [???] Parametric position on edge */
    double* u,
        /**< [???] Corresponding parametric position on face */
    double* v,
        /**< [???] Corresponding parametric position on face */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomCoordEval 
 * \brief Return parametric coordinates on face of vertex position
 *
 * Return parametric coordinates on face of vertex position
 ******************************************************************************/
void iGeom_getVtxToUV(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle vertex_handle,
        /**< [???] Vertex being queried */
    iBase_EntityHandle face_handle,
        /**< [???] Face being queried */
    double* u,
        /**< [???] Corresponding parametric position on face */
    double* v,
        /**< [???] Corresponding parametric position on face */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomCoordEval 
 * \brief Return parametric coordinates on edge of vertex position
 *
 * Return parametric coordinates on edge of vertex position
 ******************************************************************************/
void iGeom_getVtxToU(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle vertex_handle,
        /**< [???] Vertex being queried */
    iBase_EntityHandle edge_handle,
        /**< [???] Edge being queried */
    double* u,
        /**< [???] Corresponding parametric position on face */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomCoordEval 
 * \brief Return the face parametric coordinates for a parametric position on bounding edges
 *
 * Return the face parametric coordinates for a parametric position on bounding edges
 ******************************************************************************/
void iGeom_getArrUtoUV(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* edge_handles,
        /**< [???] Edges being queried */
    int edge_handles_size,
        /**< [???] Number of edges being queried */
    const iBase_EntityHandle* face_handles,
        /**< [???] Faces being queried */
    int face_handles_size,
        /**< [???] Number of faces being queried */
    const double* u_in,
        /**< [???] Parametric positions on edges */
    int u_in_size,
        /**< [???] Number of parametric positions on edges */
    int storage_order,
        /**< [???] Storage order of coordinates returned */
    double** uv,
        /**< [???] Corresponding parametric positions on faces */
    int* uv_allocated,
        /**< [???] Allocated size of parameter array */
    int* uv_size,
        /**< [???] Occupied size of parameter array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomCoordEval 
 * \brief Return parametric coordinates on faces of vertex positions
 *
 * Return parametric coordinates on faces of vertex positions
 ******************************************************************************/
void iGeom_getVtxArrToUV(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* vertex_handles,
        /**< [???] Vertices being queried */
    int vertex_handles_size,
        /**< [???] Number of vertices being queried */
    const iBase_EntityHandle* face_handles,
        /**< [???] Faces being queried */
    int face_handles_size,
        /**< [???] Number of faces being queried */
    int storage_order,
        /**< [???] Storage order of coordinates returned */
    double** uv,
        /**< [???] Corresponding parametric positions on faces */
    int* uv_allocated,
        /**< [???] Allocated size of positions array */
    int* uv_size,
        /**< [???] Occupied size of positions array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomCoordEval 
 * \brief Return parametric coordinates on edges of vertex positions
 *
 * Return parametric coordinates on edges of vertex positions
 ******************************************************************************/
void iGeom_getVtxArrToU(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* vertex_handles,
        /**< [???] Vertices being queried */
    int vertex_handles_size,
        /**< [???] Number of vertices being queried */
    const iBase_EntityHandle* edge_handles,
        /**< [???] Edges being queried */
    int edge_handles_size,
        /**< [???] Number of edges being queried */
    double** u,
        /**< [???] Corresponding parametric positions on faces */
    int* u_allocated,
        /**< [???] Allocated size of positions array */
    int* u_size,
        /**< [???] Occupied size of positions array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Return the normal at a specified parametric position
 *
 * Return the normal at a specified parametric position
 ******************************************************************************/
void iGeom_getEntNrmlUV(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [???] Entity being queried */
    double u,
        /**< [???] Parametric position being queried */
    double v,
        /**< [???] Parametric position being queried */
    double* nrml_i,
        /**< [???] Normal at specified position */
    double* nrml_j,
        /**< [???] Normal at specified position */
    double* nrml_k,
        /**< [???] Normal at specified position */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Return the normals at a specified parametric positions
 *
 * Return the normals at a specified parametric positions.  If either the number of 
 * entities or number of spatial 
 * coordinate triples is unity, then all points or entities are queried for that entity 
 * or point, respectively, otherwise each point corresponds to each entity.  storage_order 
 * should be a value in the iBase_StorageOrder enum; if input as iBase_UNKNOWN, order is in 
 * native order with respect to implementation.
 ******************************************************************************/
void iGeom_getArrNrmlUV(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* face_handles,
        /**< [???] description unknown */
    int face_handles_size,
        /**< [???] description unknown */
    int storage_order,
        /**< [???] description unknown */
    const double* parameters,
        /**< [???] description unknown */
    int parameters_size,
        /**< [???] description unknown */
    double** normals,
        /**< [???] description unknown */
    int* normals_allocated,
        /**< [???] description unknown */
    int* normals_size,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief 
 *
 ******************************************************************************/
void iGeom_getEntTgntU(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [???] description unknown */
    double u,
        /**< [???] description unknown */
    double* tgnt_i,
        /**< [???] description unknown */
    double* tgnt_j,
        /**< [???] description unknown */
    double* tgnt_k,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief 
 *
 ******************************************************************************/
void iGeom_getArrTgntU(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* edge_handles,
        /**< [???] description unknown */
    int edge_handles_size,
        /**< [???] description unknown */
    int storage_order,
        /**< [???] description unknown */
    const double* parameters,
        /**< [???] description unknown */
    int parameters_size,
        /**< [???] description unknown */
    double** tangents,
        /**< [???] description unknown */
    int* tangents_allocated,
        /**< [???] description unknown */
    int* tangents_size,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief 
 *
 ******************************************************************************/
void iGeom_getEnt1stDrvt(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [???] description unknown */
    double u,
        /**< [???] description unknown */
    double v,
        /**< [???] description unknown */
    double** drvt_u,
        /**< [???] description unknown */
    int* drvt_u_allocated,
        /**< [???] description unknown */
    int* drvt_u_size,
        /**< [???] description unknown */
    double** drvt_v,
        /**< [???] description unknown */
    int* dvrt_v_allocated,
        /**< [???] description unknown */
    int* dvrt_v_size,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief 
 *
 ******************************************************************************/
void iGeom_getArr1stDrvt(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [???] description unknown */
    int entity_handles_size,
        /**< [???] description unknown */
    int storage_order,
        /**< [???] description unknown */
    const double* uv,
        /**< [???] description unknown */
    int uv_size,
        /**< [???] description unknown */
    double** dvtr_u,
        /**< [???] description unknown */
    int* dvrt_u_allocated,
        /**< [???] description unknown */
    int* dvrt_u_size,
        /**< [???] description unknown */
    int** u_offset,
        /**< [???] description unknown */
    int* u_offset_allocated,
        /**< [???] description unknown */
    int* u_offset_size,
        /**< [???] description unknown */
    double** dvrt_v,
        /**< [???] description unknown */
    int* dvrt_v_allocated,
        /**< [???] description unknown */
    int* dvrt_v_size,
        /**< [???] description unknown */
    int** v_offset,
        /**< [???] description unknown */
    int* v_offset_allocated,
        /**< [???] description unknown */
    int* v_offset_size,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief 
 *
 ******************************************************************************/
void iGeom_getEnt2ndDrvt(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [???] description unknown */
    double u,
        /**< [???] description unknown */
    double v,
        /**< [???] description unknown */
    double** drvt_uu,
        /**< [???] description unknown */
    int* drvt_uu_allocated,
        /**< [???] description unknown */
    int* drvt_uu_size,
        /**< [???] description unknown */
    double** drvt_vv,
        /**< [???] description unknown */
    int* dvrt_vv_allocated,
        /**< [???] description unknown */
    int* dvrt_vv_size,
        /**< [???] description unknown */
    double** drvt_uv,
        /**< [???] description unknown */
    int* dvrt_uv_allocated,
        /**< [???] description unknown */
    int* dvrt_uv_size,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief 
 *
 ******************************************************************************/
void iGeom_getArr2ndDrvt(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [???] description unknown */
    int entity_handles_size,
        /**< [???] description unknown */
    int storage_order,
        /**< [???] description unknown */
    const double* uv,
        /**< [???] description unknown */
    int uv_size,
        /**< [???] description unknown */
    double** dvtr_uu,
        /**< [???] description unknown */
    int* dvrt_uu_allocated,
        /**< [???] description unknown */
    int* dvrt_uu_size,
        /**< [???] description unknown */
    int** uu_offset,
        /**< [???] description unknown */
    int* uu_offset_allocated,
        /**< [???] description unknown */
    int* uu_offset_size,
        /**< [???] description unknown */
    double** dvtr_vv,
        /**< [???] description unknown */
    int* dvrt_vv_allocated,
        /**< [???] description unknown */
    int* dvrt_vv_size,
        /**< [???] description unknown */
    int** vv_offset,
        /**< [???] description unknown */
    int* vv_offset_allocated,
        /**< [???] description unknown */
    int* vv_offset_size,
        /**< [???] description unknown */
    double** dvrt_uv,
        /**< [???] description unknown */
    int* dvrt_uv_allocated,
        /**< [???] description unknown */
    int* dvrt_uv_size,
        /**< [???] description unknown */
    int** uv_offset,
        /**< [???] description unknown */
    int* uv_offset_allocated,
        /**< [???] description unknown */
    int* uv_offset_size,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomCurvature
 * \brief 
 *
 ******************************************************************************/
void iGeom_getFcCvtrUV(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [???] description unknown */
    double u,
        /**< [???] description unknown */
    double v,
        /**< [???] description unknown */
    double* cvtr1_i,
        /**< [???] description unknown */
    double* cvtr1_j,
        /**< [???] description unknown */
    double* cvtr1_k,
        /**< [???] description unknown */
    double* cvtr2_i,
        /**< [???] description unknown */
    double* cvtr2_j,
        /**< [???] description unknown */
    double* cvtr2_k,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomCurvature
 * \brief 
 *
 ******************************************************************************/
void iGeom_getFcArrCvtrUV(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* face_handles,
        /**< [???] description unknown */
    int face_handles_size,
        /**< [???] description unknown */
    int storage_order,
        /**< [???] description unknown */
    const double* uv,
        /**< [???] description unknown */
    int uv_size,
        /**< [???] description unknown */
    double** cvtr_1,
        /**< [???] description unknown */
    int* cvtr_1_allocated,
        /**< [???] description unknown */
    int* cvtr_1_size,
        /**< [???] description unknown */
    double** cvtr_2,
        /**< [???] description unknown */
    int* cvtr_2_allocated,
        /**< [???] description unknown */
    int* cvtr_2_size,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief 
 *
 ******************************************************************************/
void iGeom_isEntPeriodic(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [???] description unknown */
    int* in_u,
        /**< [???] description unknown */
    int* in_v,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief 
 *
 ******************************************************************************/
void iGeom_isArrPeriodic(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [???] description unknown */
    int entity_handles_size,
        /**< [???] description unknown */
    int** in_uv,
        /**< [???] description unknown */
    int* in_uv_allocated,
        /**< [???] description unknown */
    int* in_uv_size,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief 
 *
 ******************************************************************************/
void iGeom_isFcDegenerate(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle face_handle,
        /**< [???] description unknown */
    int* is_degenerate,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief 
 *
 ******************************************************************************/
void iGeom_isFcArrDegenerate(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* face_handles,
        /**< [???] description unknown */
    int face_handles_size,
        /**< [???] description unknown */
    int** degenerate,
        /**< [???] description unknown */
    int* degenerate_allocated,
        /**< [???] description unknown */
    int* degenerate_size,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief 
 *
 ******************************************************************************/
void iGeom_getTolerance(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    int* type,
        /**< [???] description unknown */
    double* tolerance,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief 
 *
 ******************************************************************************/
void iGeom_getEntTolerance(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [???] description unknown */
    double* tolerance,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief 
 *
 ******************************************************************************/
void iGeom_getArrTolerance(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [???] description unknown */
    int entity_handles_size,
        /**< [???] description unknown */
    double** tolerances,
        /**< [???] description unknown */
    int* tolerances_allocated,
        /**< [???] description unknown */
    int* tolerances_size,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomIterators
 * \brief Initialize an iterator over specified entity type, topology, and size
 *
 * Initialize an iterator over specified entity type, topology, and size,
 * for a specified set or instance.  Iterator returned can be used as input
 * to functions returning the entity for the iterator.  If all entities of 
 * a specified type and/or topology are to be iterated, specify 
 * iBase_ALL_TYPES or iGeom_ALL_TOPOLOGIES, respectively.  Specified type 
 * or topology must be a value in the iBase_EntityType or 
 * iGeom_EntityTopology enumerations, respectively.
 ******************************************************************************/
void iGeom_initEntIter(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle entity_set_handle,
        /**< [???] Entity set being iterated */
    int entity_dimension,
        /**< [???] dimension of entity to iterate */
    iBase_EntityIterator* entity_iterator,
        /**< [???] Pointer to iterator returned from function */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomIterators
 * \brief Initialize an array iterator over specified entity type, topology, and 
 *        size
 *
 * Initialize an array iterator over specified entity type, topology, and 
 * size, for a specified set or instance.  Iterator returned can be used 
 * as input to functions returning entities for the iterator.  If all 
 * entities of a specified type and/or topology are to be iterated, 
 * specify iBase_ALL_TYPES or iGeom_ALL_TOPOLOGIES, respectively.  
 * Specified type or topology must be a value in the iBase_EntityType or 
 * iGeom_EntityTopology enumerations, respectively.
 ******************************************************************************/
void iGeom_initEntArrIter(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle entity_set_handle,
        /**< [???] Entity set being iterated */
    int entity_dimension,
        /**< [???] dimension of entity to iterate */
    int requested_array_size,
        /**< [???] Size of chunks of handles returned for each value of the
        iterator */
    iBase_EntityArrIterator* entArr_iterator,
        /**< [???] Pointer to iterator returned from function */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomIterators
 * \brief Get entity corresponding to an iterator and increment iterator
 *
 * Get the entity corresponding to an array iterator, and increment the 
 * iterator.  Also return whether the next value of the iterator has
 * an entity (if non-zero, next iterator value is the end of the
 * iteration).
 ******************************************************************************/
void iGeom_getNextEntIter(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityIterator entity_iterator,
        /**< [???] Iterator being queried */
    iBase_EntityHandle* entity_handle,
        /**< [???] Pointer to an entity handle corresponding to the current
        value of iterator */
    int* has_data,
        /**< [???] Pointer to a flag indicating if the value(s) returned in
        entity_handles are valid. A non-zero value indicates the value(s) are
        valid. A zero value indicates the value(s) are NOT valid. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomIterators
 * \brief Get entities contained in array iterator and increment iterator
 *
 * Get the entities contained in an array iterator, and increment the 
 * iterator.  Also return whether the next value of the iterator has
 * any entities (if non-zero, next iterator value is the end of the
 * iteration).
 ******************************************************************************/
void iGeom_getNextEntArrIter(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityArrIterator entArr_iterator,
        /**< [???] Iterator being queried */
    iBase_EntityHandle** entity_handles,
        /**< [???] Pointer to array of entity handles contained in current
        value of iterator */
    int* entity_handles_allocated,
        /**< [???] Pointer to allocated size of entity_handles array
 */
    int* entity_handles_size,
        /**< [???] Pointer to occupied size of entity_handles array
 */
    int* has_data,
        /**< [???] Pointer to a flag indicating if the value(s) returned in
        entity_handles are valid. A non-zero value indicates the value(s) are
        valid. A zero value indicates the value(s) are NOT valid. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomIterators
 * \brief Reset the iterator
 *
 * Reset the iterator
 ******************************************************************************/
void iGeom_resetEntIter(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityIterator entity_iterator,
        /**< [???] Iterator to reset */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomIterators
 * \brief Reset the array iterator
 *
 * Reset the array iterator
 ******************************************************************************/
void iGeom_resetEntArrIter(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityArrIterator entArr_iterator,
        /**< [???] Iterator to reset */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomIterators
 * \brief Destroy the specified iterator
 *
 * Destroy the specified iterator
 ******************************************************************************/
void iGeom_endEntIter(
    iBase_EntityIterator entity_iterator,
        /**< [???] iGeom instance handle */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomIterators
 * \brief Destroy the specified array iterator
 *
 * Destroy the specified array iterator
 ******************************************************************************/
void iGeom_endEntArrIter(
    iBase_EntityArrIterator entArr_iterator,
        /**< [???] iGeom instance handle */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief 
 *
 ******************************************************************************/
void iGeom_copyEnt(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle source,
        /**< [???] description unknown */
    iBase_EntityHandle* copy,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief 
 *
 ******************************************************************************/
void iGeom_sweepEntAboutAxis(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle geom_entity,
        /**< [???] description unknown */
    double angle,
        /**< [???] description unknown */
    double axis_normal_x,
        /**< [???] description unknown */
    double axis_normal_y,
        /**< [???] description unknown */
    double axis_normal_z,
        /**< [???] description unknown */
    iBase_EntityHandle* geom_entity2,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Delete all entities and sets
 *
 ******************************************************************************/
void iGeom_deleteAll(
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Delete specified entity
 *
 * Delete specified entity
 ******************************************************************************/
void iGeom_deleteEnt(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [???] Entity to be deleted */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Create a sphere centered on the origin
 *
 ******************************************************************************/
void iGeom_createSphere(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    double radius,
        /**< [???] description unknown */
    iBase_EntityHandle* sphere_handle_out,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Create a prism centered on the origin
 *
 ******************************************************************************/
void iGeom_createPrism(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    double height,
        /**< [???] description unknown */
    int n_sides,
        /**< [???] description unknown */
    double major_rad,
        /**< [???] description unknown */
    double minor_rad,
        /**< [???] description unknown */
    iBase_EntityHandle* prism_handle_out,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief 
 *
 ******************************************************************************/
void iGeom_createBrick(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    double x,
        /**< [???] description unknown */
    double y,
        /**< [???] description unknown */
    double z,
        /**< [???] description unknown */
    iBase_EntityHandle* geom_entity,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief 
 *
 ******************************************************************************/
void iGeom_createCylinder(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    double height,
        /**< [???] description unknown */
    double major_rad,
        /**< [???] description unknown */
    double minor_rad,
        /**< [???] description unknown */
    iBase_EntityHandle* geom_entity,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief 
 *
 ******************************************************************************/
void iGeom_createTorus(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    double major_rad,
        /**< [???] description unknown */
    double minor_rad,
        /**< [???] description unknown */
    iBase_EntityHandle* geom_entity,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief 
 *
 ******************************************************************************/
void iGeom_moveEnt(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle geom_entity,
        /**< [???] description unknown */
    double x,
        /**< [???] description unknown */
    double y,
        /**< [???] description unknown */
    double z,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief 
 *
 ******************************************************************************/
void iGeom_rotateEnt(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle geom_entity,
        /**< [???] description unknown */
    double angle,
        /**< [???] description unknown */
    double axis_normal_x,
        /**< [???] description unknown */
    double axis_normal_y,
        /**< [???] description unknown */
    double axis_normal_z,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief 
 *
 ******************************************************************************/
void iGeom_reflectEnt(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle geom_entity,
        /**< [???] description unknown */
    double plane_normal_x,
        /**< [???] description unknown */
    double plane_normal_y,
        /**< [???] description unknown */
    double plane_normal_z,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief 
 *
 ******************************************************************************/
void iGeom_scaleEnt(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle geom_entity,
        /**< [???] description unknown */
    double scale_x,
        /**< [???] description unknown */
    double scale_y,
        /**< [???] description unknown */
    double scale_z,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief 
 *
 ******************************************************************************/
void iGeom_uniteEnts(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* geom_entities,
        /**< [???] description unknown */
    int geom_entities_size,
        /**< [???] description unknown */
    iBase_EntityHandle* geom_entity,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief 
 *
 ******************************************************************************/
void iGeom_subtractEnts(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle blank,
        /**< [???] description unknown */
    iBase_EntityHandle tool,
        /**< [???] description unknown */
    iBase_EntityHandle* geom_entity,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief 
 *
 ******************************************************************************/
void iGeom_intersectEnts(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity2,
        /**< [???] description unknown */
    iBase_EntityHandle entity1,
        /**< [???] description unknown */
    iBase_EntityHandle* geom_entity,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief 
 *
 ******************************************************************************/
void iGeom_sectionEnt(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle geom_entity,
        /**< [???] description unknown */
    double plane_normal_x,
        /**< [???] description unknown */
    double plane_normal_y,
        /**< [???] description unknown */
    double plane_normal_z,
        /**< [???] description unknown */
    double offset,
        /**< [???] description unknown */
    int reverse,
        /**< [???] description unknown */
    iBase_EntityHandle* geom_entity2,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief 
 *
 ******************************************************************************/
void iGeom_imprintEnts(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* geom_entities,
        /**< [???] description unknown */
    int geom_entities_size,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief 
 *
 ******************************************************************************/
void iGeom_mergeEnts(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* geom_entities,
        /**< [???] description unknown */
    int geom_entities_size,
        /**< [???] description unknown */
    double tolerance,
        /**< [???] description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iGeomEntitySets
 * \brief Create an entity set
 *
 * Create an entity set, either ordered (isList=1) or unordered 
 * (isList=0).  Unordered entity sets can contain a given entity or 
 * set only once.
 ******************************************************************************/
void iGeom_createEntSet(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    int isList,
        /**< [???] If non-zero, an ordered list is created, otherwise an
        unordered set is created. */
    iBase_EntitySetHandle* entity_set_created,
        /**< [???] Entity set created by function */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Destroy an entity set
 *
 * Destroy an entity set
 ******************************************************************************/
void iGeom_destroyEntSet(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle entity_set,
        /**< [???] Entity set to be destroyed */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Return whether a specified set is ordered or unordered
 *
 * Return whether a specified set is ordered (*is_list=1) or 
 * unordered (*is_list=0)
 ******************************************************************************/
void iGeom_isList(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle entity_set,
        /**< [???] Entity set being queried */
    int* is_list,
        /**< [???] Pointer to flag returned from function */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the number of entity sets contained in a set or interface
 *
 * Get the number of entity sets contained in a set or interface.  If
 * a set is input which is not the root set, num_hops indicates the 
 * maximum number of contained sets from entity_set_handle to one of the
 * contained sets, inclusive of the contained set.
 ******************************************************************************/
void iGeom_getNumEntSets(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle entity_set_handle,
        /**< [???] Entity set being queried */
    int num_hops,
        /**< [???] Maximum hops from entity_set_handle to contained set,
        inclusive of the contained set */
    int* num_sets,
        /**< [???] Pointer to the number of sets returned from function */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the entity sets contained in a set or interface
 *
 * Get the entity sets contained in a set or interface.  If
 * a set is input which is not the root set, num_hops indicates the 
 * maximum number of contained sets from entity_set_handle to one of the
 * contained sets, inclusive of the contained set.
 ******************************************************************************/
void iGeom_getEntSets(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle entity_set_handle,
        /**< [???] Entity set being queried */
    int num_hops,
        /**< [???] Maximum hops from entity_set_handle to contained set,
        inclusive of the contained set */
    iBase_EntitySetHandle** contained_set_handles,
        /**< [???] Pointer to array of set handles returned from function
 */
    int* contained_set_handles_allocated,
        /**< [???] Pointer to allocated length of contained_set_handles array
 */
    int* contained_set_handles_size,
        /**< [???] Pointer to occupied length of contained_set_handles array
 */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Add an entity to a set
 *
 * Add an entity to a set
 ******************************************************************************/
void iGeom_addEntToSet(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [???] The entity being added */
    iBase_EntitySetHandle entity_set,
        /**< [???] Pointer to the set being added to */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Remove an entity from a set
 *
 * Remove an entity from a set
 ******************************************************************************/
void iGeom_rmvEntFromSet(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [???] The entity being removed */
    iBase_EntitySetHandle entity_set,
        /**< [???] Pointer to the set being removed from */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Add an array of entities to a set
 *
 * Add an array of entities to a set
 ******************************************************************************/
void iGeom_addEntArrToSet(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [???] Array of entities being added */
    int entity_handles_size,
        /**< [???] Number of entities in entity_handles array */
    iBase_EntitySetHandle entity_set,
        /**< [???] Pointer to the set being added to */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Remove an array of entities from a set
 *
 * Remove an array of entities from a set
 ******************************************************************************/
void iGeom_rmvEntArrFromSet(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [???] Array of entities being remove */
    int entity_handles_size,
        /**< [???] Number of entities in entity_handles array */
    iBase_EntitySetHandle entity_set,
        /**< [???] Pointer to the set being removed from */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Add an entity set to a set
 *
 * Add an entity set to a set
 ******************************************************************************/
void iGeom_addEntSet(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle entity_set_to_add,
        /**< [???] The entity set being added */
    iBase_EntitySetHandle entity_set_handle,
        /**< [???] Pointer to the set being added to */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Remove an entity set from a set
 *
 * Remove an entity set from a set
 ******************************************************************************/
void iGeom_rmvEntSet(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle entity_set_to_remove,
        /**< [???] The entity set being removed */
    iBase_EntitySetHandle entity_set_handle,
        /**< [???] Pointer to the set being removed from */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Return whether an entity is contained in another set
 *
 * Return whether an entity is contained (*is_contained=1) or not 
 * contained (*is_contained=0) in another set
 ******************************************************************************/
void iGeom_isEntContained(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle containing_entity_set,
        /**< [???] Entity set being queried */
    iBase_EntityHandle contained_entity,
        /**< [???] Entity potentially contained in containing_entity_set
 */
    int* is_contained,
        /**< [???] Pointer to flag returned from function */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Return whether entities are contained in a set
 *
 * Return whether each entity is contained in the set.
 ******************************************************************************/
void iGeom_isEntArrContained(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle containing_set,
        /**< [???] Entity set being queried */
    const iBase_EntityHandle* entity_handles,
        /**< [???] List of entities for which to check containment. */
    int num_entity_handles,
        /**< [???] size of entity_handles array */
    int** is_contained,
        /**< [???] One value for each input entity, 1 if contained in set, zero
        otherwise. */
    int* is_contained_allocated,
        /**< [???] allocated size of is_contained array */
    int* is_contained_size,
        /**< [???] occupied size of is_contained array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Return whether an entity set is contained in another set
 *
 * Return whether a set is contained (*is_contained=1) or not contained
 * (*is_contained=0) in another set
 ******************************************************************************/
void iGeom_isEntSetContained(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle containing_entity_set,
        /**< [???] Entity set being queried */
    iBase_EntitySetHandle contained_entity_set,
        /**< [???] Entity set potentially contained in containing_entity_set */
    int* is_contained,
        /**< [???] Pointer to flag returned from function */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Add parent/child links between two sets
 *
 * Add parent/child links between two sets.  Makes parent point to child
 * and child point to parent.
 ******************************************************************************/
void iGeom_addPrntChld(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle parent_entity_set,
        /**< [???] Pointer to parent set */
    iBase_EntitySetHandle child_entity_set,
        /**< [???] Pointer to child set */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Remove parent/child links between two sets
 *
 * Remove parent/child links between two sets.
 ******************************************************************************/
void iGeom_rmvPrntChld(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle parent_entity_set,
        /**< [???] Pointer to parent set */
    iBase_EntitySetHandle child_entity_set,
        /**< [???] Pointer to child set */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Return whether two sets are related by parent/child links
 *
 * Return whether two sets are related (*is_child=1) or not (*is_child=0)
 * by parent/child links
 ******************************************************************************/
void iGeom_isChildOf(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle parent_entity_set,
        /**< [???] Pointer to parent set */
    iBase_EntitySetHandle child_entity_set,
        /**< [???] Pointer to child set */
    int* is_child,
        /**< [???] Pointer to flag returned from function */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the number of child sets linked from a specified set
 *
 * Get the number of child sets linked from a specified set.  If num_hops
 * is non-zero, this represents the maximum hops from entity_set to any
 * child in the count.
 ******************************************************************************/
void iGeom_getNumChld(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle entity_set,
        /**< [???] Entity set being queried */
    int num_hops,
        /**< [???] Maximum hops from entity_set_handle to child set, inclusive
        of the child set */
    int* num_child,
        /**< [???] Pointer to number of children returned from function */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the number of parent sets linked from a specified set
 *
 * Get the number of parent sets linked from a specified set.  If num_hops
 * is non-zero, this represents the maximum hops from entity_set to any
 * parent in the count.
 ******************************************************************************/
void iGeom_getNumPrnt(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle entity_set,
        /**< [???] Entity set being queried */
    int num_hops,
        /**< [???] Maximum hops from entity_set_handle to parent set, inclusive
        of the parent set */
    int* num_parent,
        /**< [???] Pointer to number of parents returned from function */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the child sets linked from a specified set
 *
 * Get the child sets linked from a specified set.  If num_hops
 * is non-zero, this represents the maximum hops from entity_set to any
 * child.
 ******************************************************************************/
void iGeom_getChldn(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle from_entity_set,
        /**< [???] Entity set being queried */
    int num_hops,
        /**< [???] Maximum hops from entity_set_handle to child set, inclusive
        of the child set */
    iBase_EntitySetHandle** entity_set_handles,
        /**< [???] Pointer to array of child sets returned from function */
    int* entity_set_handles_allocated,
        /**< [???] Pointer to allocated size of entity_set_handles array */
    int* entity_set_handles_size,
        /**< [???] Pointer to occupied size of entity_set_handles array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the parent sets linked from a specified set
 *
 * Get the parent sets linked from a specified set.  If num_hops
 * is non-zero, this represents the maximum hops from entity_set to any
 * parent.
 ******************************************************************************/
void iGeom_getPrnts(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle from_entity_set,
        /**< [???] Entity set being queried */
    int num_hops,
        /**< [???] Maximum hops from entity_set_handle to parent set, inclusive
        of the parent set */
    iBase_EntitySetHandle** entity_set_handles,
        /**< [???] Pointer to array of parent sets returned from function */
    int* entity_set_handles_allocated,
        /**< [???] Pointer to allocated size of entity_set_handles array */
    int* entity_set_handles_size,
        /**< [???] Pointer to occupied size of entity_set_handles array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Create a tag with specified name, size, and type
 *
 * Create a tag with specified name, size, and type.  Tag size is in
 * units of size of tag_type data types.  Value input for tag type must be 
 * value in iBase_TagType enumeration.
 ******************************************************************************/
void iGeom_createTag(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const char* tag_name,
        /**< [???] Character string indicating tag name */
    int tag_size,
        /**< [???] Size of each tag value, in units of number of tag_type entities */
    int tag_type,
        /**< [???] Data type for data stored in this tag */
    iBase_TagHandle* tag_handle,
        /**< [???] Pointer to tag handle returned from function */
    int* err,
        /**< [out] Returned Error status (see iBase_ErrorType) */
    int tag_name_len
        /**< [???] Length of tag name string */
);

/***************************************************************************//**
 * \brief Destroy a tag
 *
 * Destroy a tag.  If forced is non-zero and entities still have values
 * set for this tag, tag is deleted anyway and those values disappear,
 * otherwise tag is not deleted.
 ******************************************************************************/
void iGeom_destroyTag(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_TagHandle tag_handle,
        /**< [???] Handle of tag to be deleted */
    int forced,
        /**< [???] If non-zero, delete the tag even if entities have values set
        for that tag */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the name for a given tag handle
 *
 * Get the name for a given tag handle
 ******************************************************************************/
void iGeom_getTagName(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_TagHandle tag_handle,
        /**< [???] Tag handle being queried */
    char* name,
        /**< [???] Pointer to character string to store name returned from
        function */
    int* err,
        /**< [out] Returned Error status (see iBase_ErrorType) */
    int name_len
        /**< [???] Length of character string input to function */
);

/***************************************************************************//**
 * \brief Get size of a tag in units of numbers of tag data type
 *
 * Get size of a tag in units of numbers of tag data type
 ******************************************************************************/
void iGeom_getTagSizeValues(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_TagHandle tag_handle,
        /**< [???] Handle of tag being queried */
    int* tag_size,
        /**< [???] Pointer to tag size returned from function */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get size of a tag in units of bytes
 *
 * Get size of a tag in units of bytes
 ******************************************************************************/
void iGeom_getTagSizeBytes(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_TagHandle tag_handle,
        /**< [???] Handle of tag being queried */
    int* tag_size,
        /**< [???] Pointer to tag size returned from function */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get a the handle of an existing tag with the specified name
 *
 * Get a the handle of an existing tag with the specified name
 ******************************************************************************/
void iGeom_getTagHandle(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const char* tag_name,
        /**< [???] Name of tag being queried */
    iBase_TagHandle* tag_handle,
        /**< [???] Pointer to tag handle returned from function */
    int* err,
        /**< [out] Returned Error status (see iBase_ErrorType) */
    int tag_name_len
        /**< [???] Length of tag name string */
);

/***************************************************************************//**
 * \brief Get the data type of the specified tag handle
 *
 * Get the data type of the specified tag handle.  Tag type is a value in
 * the iBase_TagType enumeration.
 ******************************************************************************/
void iGeom_getTagType(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_TagHandle tag_handle,
        /**< [???] Handle for the tag being queried */
    int* tag_type,
        /**< [???] Pointer to tag type returned from function */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Set a tag value of arbitrary type on an entity set
 *
 * Set a tag value of arbitrary type on an entity set.  Tag data is 
 * passed as char* type,
 * but really represents pointer to arbitrary data.
 ******************************************************************************/
void iGeom_setEntSetData(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle entity_set_handle,
        /**< [???] Entity set on which tag is being set */
    iBase_TagHandle tag_handle,
        /**< [???] Tag being set on an entity set */
    const char* tag_value,
        /**< [???] Pointer to tag data being set on entity set */
    int tag_value_size,
        /**< [???] Size in bytes of tag data */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Set a tag value of integer type on an entity set
 *
 * Set a tag value of integer type on an entity set.
 ******************************************************************************/
void iGeom_setEntSetIntData(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle entity_set,
        /**< [???] Entity set on which tag is being set */
    iBase_TagHandle tag_handle,
        /**< [???] Tag being set on an entity set */
    int tag_value,
        /**< [???] Tag value being set on entity set */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Set a tag value of double type on an entity set
 *
 * Set a tag value of double type on an entity set.
 ******************************************************************************/
void iGeom_setEntSetDblData(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle entity_set,
        /**< [???] Entity set on which tag is being set */
    iBase_TagHandle tag_handle,
        /**< [???] Tag being set on an entity set */
    double tag_value,
        /**< [???] Tag value being set on entity set */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Set a tag value of entity handle type on an entity set
 *
 * Set a tag value of entity handle type on an entity set.
 ******************************************************************************/
void iGeom_setEntSetEHData(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle entity_set,
        /**< [???] Entity set on which tag is being set */
    iBase_TagHandle tag_handle,
        /**< [???] Tag being set on an entity set */
    iBase_EntityHandle tag_value,
        /**< [???] Tag value being set on entity set */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the value of a tag of arbitrary type on an entity set
 *
 * Get the value of a tag of arbitrary type on an entity set.  Tag data 
 * is passed back as void* type, but really represents arbitrary data.
 ******************************************************************************/
void iGeom_getEntSetData(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle entity_set_handle,
        /**< [???] Entity set on which tag is being set */
    iBase_TagHandle tag_handle,
        /**< [???] Tag being set on an entity set */
    void* tag_value,
        /**< [???] Pointer to tag data array being queried */
    int* tag_value_allocated,
        /**< [???] Pointer to tag data array allocated size */
    int* tag_value_size,
        /**< [???] Pointer to tag data array occupied size */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the value of a tag of integer type on an entity set
 *
 * Get the value of a tag of integer type on an entity set.
 ******************************************************************************/
void iGeom_getEntSetIntData(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle entity_set,
        /**< [???] Entity set on which tag is being set */
    iBase_TagHandle tag_handle,
        /**< [???] Tag being set on an entity set */
    int* out_data,
        /**< [???] Pointer to tag value returned from function */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the value of a tag of double type on an entity set
 *
 * Get the value of a tag of double type on an entity set.
 ******************************************************************************/
void iGeom_getEntSetDblData(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle entity_set,
        /**< [???] Entity set on which tag is being set */
    iBase_TagHandle tag_handle,
        /**< [???] Tag being set on an entity set */
    double* out_data,
        /**< [???] Pointer to tag value returned from function */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the value of a tag of entity handle type on an entity set
 *
 * Get the value of a tag of entity handle type on an entity set.
 ******************************************************************************/
void iGeom_getEntSetEHData(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle entity_set,
        /**< [???] Entity set on which tag is being set */
    iBase_TagHandle tag_handle,
        /**< [???] Tag being set on an entity set */
    iBase_EntityHandle* out_data,
        /**< [???] Pointer to tag value returned from function */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get all the tags associated with a specified entity set
 *
 * Get all the tags associated with a specified entity set
 ******************************************************************************/
void iGeom_getAllEntSetTags(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle entity_set_handle,
        /**< [???] Entity being queried */
    iBase_TagHandle** tag_handles,
        /**< [???] Pointer to array of tag_handles returned from function
 */
    int* tag_handles_allocated,
        /**< [???] Pointer to allocated size of tag_handles array */
    int* tag_handles_size,
        /**< [???] Pointer to occupied size of tag_handles array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Remove a tag value from an entity set
 *
 * Remove a tag value from an entity set
 ******************************************************************************/
void iGeom_rmvEntSetTag(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle entity_set_handle,
        /**< [???] Entity set from which tag is being removed */
    iBase_TagHandle tag_handle,
        /**< [???] Tag handle of tag being removed */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get tag values of arbitrary type for an array of entities
 *
 * Get tag values of arbitrary type for an array of entities.  Tag data 
 * is returned as void* type, but really represents arbitrary data.
 ******************************************************************************/
void iGeom_getArrData(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [???] Entity array on which tag is being set */
    int entity_handles_size,
        /**< [???] Number of entities in array */
    iBase_TagHandle tag_handle,
        /**< [???] Tag being set on an entity */
    void* tag_values,
        /**< [???] Pointer to tag data array being returned from function */
    int* tag_values_allocated,
        /**< [???] Pointer to allocated size of tag data array */
    int* tag_values_size,
        /**< [???] Pointer to occupied size of tag data array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get tag values of integer type for an array of entities
 *
 * Get tag values of integer type for an array of entities.
 ******************************************************************************/
void iGeom_getIntArrData(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [???] Entity array on which tag is being set */
    int entity_handles_size,
        /**< [???] Number of entities in array */
    iBase_TagHandle tag_handle,
        /**< [???] Tag being set on an entity */
    int** tag_values,
        /**< [???] Pointer to tag data array being returned from function */
    int* tag_values_allocated,
        /**< [???] Pointer to allocated size of tag data array */
    int* tag_values_size,
        /**< [???] Pointer to occupied size of tag data array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get tag values of double type for an array of entities
 *
 * Get tag values of double type for an array of entities.
 ******************************************************************************/
void iGeom_getDblArrData(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [???] Entity array on which tag is being set */
    int entity_handles_size,
        /**< [???] Number of entities in array */
    iBase_TagHandle tag_handle,
        /**< [???] Tag being set on an entity */
    double** tag_values,
        /**< [???] Pointer to tag data array being returned from function */
    int* tag_values_allocated,
        /**< [???] Pointer to allocated size of tag data array */
    int* tag_values_size,
        /**< [???] Pointer to occupied size of tag data array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get tag values of entity handle type for an array of entities
 *
 * Get tag values of entity handle type for an array of entities.
 ******************************************************************************/
void iGeom_getEHArrData(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [???] Entity array on which tag is being set */
    int entity_handles_size,
        /**< [???] Number of entities in array */
    iBase_TagHandle tag_handle,
        /**< [???] Tag being set on an entity */
    iBase_EntityHandle** tag_value,
        /**< [???] Pointer to tag data array being returned from function */
    int* tag_value_allocated,
        /**< [???] Pointer to allocated size of tag data array */
    int* tag_value_size,
        /**< [???] Pointer to occupied size of tag data array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Set tag values of arbitrary type on an array of entities
 *
 * Set tag values of arbitrary type on an array of entities.  Tag data is 
 * passed as char* type, but really represents pointer to arbitrary data.
 ******************************************************************************/
void iGeom_setArrData(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [???] Entity array on which tag is being set */
    int entity_handles_size,
        /**< [???] Number of entities in array */
    iBase_TagHandle tag_handle,
        /**< [???] Tag being set on an entity */
    const char* tag_values,
        /**< [???] Pointer to tag data being set on entity */
    int tag_values_size,
        /**< [???] Size in total bytes of tag data */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Set tag values of integer type on an array of entities
 *
 * Set tag values of integer type on an array of entities.
 ******************************************************************************/
void iGeom_setIntArrData(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [???] Entity array on which tag is being set */
    int entity_handles_size,
        /**< [???] Number of entities in array */
    iBase_TagHandle tag_handle,
        /**< [???] Tag being set on an entity */
    const int* tag_values,
        /**< [???] Pointer to tag data being set on entities */
    int tag_values_size,
        /**< [???] Size in total number of integers of tag data */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Set tag values of double type on an array of entities
 *
 * Set tag values of double type on an array of entities.
 ******************************************************************************/
void iGeom_setDblArrData(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [???] Entity array on which tag is being set */
    int entity_handles_size,
        /**< [???] Number of entities in array */
    iBase_TagHandle tag_handle,
        /**< [???] Tag being set on an entity */
    const double* tag_values,
        /**< [???] Pointer to tag data being set on entities */
    const int tag_values_size,
        /**< [???] Size in total number of doubles of tag data */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Set tag values of entity handle type on an array of entities
 *
 * Set tag values of entity handle type on an array of entities.
 ******************************************************************************/
void iGeom_setEHArrData(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [???] Entity array on which tag is being set */
    int entity_handles_size,
        /**< [???] Number of entities in array */
    iBase_TagHandle tag_handle,
        /**< [???] Tag being set on an entity */
    const iBase_EntityHandle* tag_values,
        /**< [???] Pointer to tag data being set on entities */
    int tag_values_size,
        /**< [???] Size in total number of entity handles of tag data */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Remove a tag value from an array of entities
 *
 * Remove a tag value from an array of entities
 ******************************************************************************/
void iGeom_rmvArrTag(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle* entity_handles,
        /**< [???] Entity from which tag is being removed */
    int entity_handles_size,
        /**< [???] Number of entities in entity array */
    iBase_TagHandle tag_handle,
        /**< [???] Tag handle of tag being removed */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the value of a tag of arbitrary type on an entity
 *
 * Get the value of a tag of arbitrary type on an entity.  Tag data 
 * is passed back as void* type, but really represents arbitrary data.
 ******************************************************************************/
void iGeom_getData(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [???] Entity on which tag is being set */
    iBase_TagHandle tag_handle,
        /**< [???] Tag being set on an entity */
    void* tag_value,
        /**< [???] Pointer to tag data array being queried */
    int* tag_value_allocated,
        /**< [???] Pointer to tag data array allocated size */
    int* tag_value_size,
        /**< [???] Pointer to tag data array occupied size */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the value of a tag of integer type on an entity
 *
 * Get the value of a tag of integer type on an entity.
 ******************************************************************************/
void iGeom_getIntData(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [???] Entity on which tag is being set */
    iBase_TagHandle tag_handle,
        /**< [???] Tag being set on an entity */
    int* out_data,
        /**< [???] Pointer to tag value returned from function */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the value of a tag of double type on an entity
 *
 * Get the value of a tag of double type on an entity.
 ******************************************************************************/
void iGeom_getDblData(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    const iBase_EntityHandle entity_handle,
        /**< [???] Entity on which tag is being set */
    const iBase_TagHandle tag_handle,
        /**< [???] Tag being set on an entity */
    double* out_data,
        /**< [???] Pointer to tag value returned from function */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get the value of a tag of entity handle type on an entity
 *
 * Get the value of a tag of entity handle type on an entity.
 ******************************************************************************/
void iGeom_getEHData(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [???] Entity on which tag is being set */
    iBase_TagHandle tag_handle,
        /**< [???] Tag being set on an entity */
    iBase_EntityHandle* out_data,
        /**< [???] Pointer to tag value returned from function */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Set a tag value of arbitrary type on an entity
 *
 * Set a tag value of arbitrary type on an entity.  Tag data is 
 * passed as char* type, but really represents pointer to arbitrary data.
 ******************************************************************************/
void iGeom_setData(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [???] Entity on which tag is being set */
    iBase_TagHandle tag_handle,
        /**< [???] Tag being set on an entity */
    const char* tag_value,
        /**< [???] Pointer to tag data being set on entity */
    int tag_value_size,
        /**< [???] Size in bytes of tag data */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Set a tag value of integer type on an entity
 *
 * Set a tag value of integer type on an entity.
 ******************************************************************************/
void iGeom_setIntData(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [???] Entity on which tag is being set */
    iBase_TagHandle tag_handle,
        /**< [???] Tag being set on an entity */
    int tag_value,
        /**< [???] Tag value being set on entity */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Set a tag value of double type on an entity
 *
 * Set a tag value of double type on an entity.
 ******************************************************************************/
void iGeom_setDblData(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [???] Entity on which tag is being set */
    iBase_TagHandle tag_handle,
        /**< [???] Tag being set on an entity */
    double tag_value,
        /**< [???] Tag value being set on entity */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Set a tag value of entity handle type on an entity
 *
 * Set a tag value of entity handle type on an entity.
 ******************************************************************************/
void iGeom_setEHData(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [???] Entity on which tag is being set */
    iBase_TagHandle tag_handle,
        /**< [???] Tag being set on an entity */
    iBase_EntityHandle tag_value,
        /**< [???] Tag value being set on entity */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Get all the tags associated with a specified entity handle
 *
 * Get all the tags associated with a specified entity handle
 ******************************************************************************/
void iGeom_getAllTags(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [???] Entity being queried */
    iBase_TagHandle** tag_handles,
        /**< [???] Pointer to array of tag_handles returned from function */
    int* tag_handles_allocated,
        /**< [???] Pointer to allocated size of tag_handles array */
    int* tag_handles_size,
        /**< [???] Pointer to occupied size of tag_handles array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Remove a tag value from an entity
 *
 * Remove a tag value from an entity
 ******************************************************************************/
void iGeom_rmvTag(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntityHandle entity_handle,
        /**< [???] Entity from which tag is being removed */
    iBase_TagHandle tag_handle,
        /**< [???] Tag handle of tag being removed */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Subtract contents of one entity set from another
 *
 * Subtract contents of one entity set from another
 ******************************************************************************/
void iGeom_subtract(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle entity_set_1,
        /**< [???] Entity set from which other set is being subtracted */
    iBase_EntitySetHandle entity_set_2,
        /**< [???] Entity set being subtracted from other set */
    iBase_EntitySetHandle* result_entity_set,
        /**< [???] Pointer to entity set returned from function */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Intersect contents of one entity set with another
 *
 * Intersect contents of one entity set with another
 ******************************************************************************/
void iGeom_intersect(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle entity_set_1,
        /**< [???] Entity set being intersected with another */
    iBase_EntitySetHandle entity_set_2,
        /**< [???] Entity set being intersected with another */
    iBase_EntitySetHandle* result_entity_set,
        /**< [???] Pointer to entity set returned from function */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \brief Unite contents of one entity set with another
 *
 * Unite contents of one entity set with another
 ******************************************************************************/
void iGeom_unite(
    iGeom_Instance instance,
        /**< [in] iGeom instance handle */
    iBase_EntitySetHandle entity_set_1,
        /**< [???] Entity set being united with another */
    iBase_EntitySetHandle entity_set_2,
        /**< [???] Entity set being united with another */
    iBase_EntitySetHandle* result_entity_set,
        /**< [???] Pointer to entity set returned from function */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \page igeom iGeom: ITAPS Geometry Interface 
 *
 * The ITAPS Geometry Interface iGeom provides a common interface for
 * accessing geometry and data associated with a mesh.  Applications written
 * to use this interface can use a variety of implementations, choosing
 * the one that best meets its needs.  They can also use tools written
 * to this interface.
 *
 * \section ITAPS Data Model
 *
 * The ITAPS interfaces use a data model composed of four basic data types:\n
 * \em Entity: basic topological entities in a model, e.g. vertices, 
 * edges, faces, regions. \n
 * \em Entity \em Set: arbitrary grouping of other entities and sets. 
 * Entity sets also support parent/child relations with other sets which
 * are distinct from entities contained in those sets.  Parent/child links
 * can be used to embed graph relationships between sets, e.g. to 
 * represent topological relationships between the sets. \n
 * \em Interface: the object with which model is associated and on which
 * functions in iGeom are called. \n
 * \em Tag: application data associated with objects of any of the other 
 * data types.  Each tag has a designated name, size, and data type.
 *
 * \section ITAPS Entity Type
 * Each entity has a specific Entity Type.  The Entity 
 * Type is one of VERTEX, EDGE, FACE, and REGION, and is synonymous with
 * the topological dimension of the entity.  Entity Type is an enumerated
 * type in the iBase_EntityType enumeration.
 *
 * \section ITAPS Entity-, Array-, and Iterator-Based Access
 *
 * The iGeom interface provides functions for accessing entities
 * individually, as arrays of entities, or using iterators.  These access
 * methods have different memory versus execution time tradeoffs, 
 * depending on the implementation.
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup iGeom iGeom
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup iGeomInitialization Initialization
 * \ingroup iGeom
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup iGeomEntities Entities 
 * \ingroup iGeom
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup iGeomEntitySets Entity Sets
 * \ingroup iGeom
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup iGeomEntitySetOperators Entity Set Operators
 * \ingroup EntitySets
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup iGeomAdjacencies Adjacencies
 * \ingroup iGeom
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup iGeomCurvature Curvature
 * \ingroup iGeom
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup iGeomIterators Entity Iterators
 * \ingroup iGeom
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup iGeomTags Tags
 * \ingroup iGeom
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup iGeomTagData Tag Data
 * \ingroup iGeom
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup iGeomTagsOnEnts Tag Data On Entities
 * \ingroup iGeomTagData
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup iGeomTagsOnSets Tag Data On Entity Sets
 * \ingroup iGeomTagData
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup TagsOnArr Tag Data On Arrays of Entities
 * \ingroup TagData
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup iGeomParentChildLinks Parent Child Links
 * \ingroup iGeom
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup iGeomCoordEval Parametric<-->Geometric Coordinate Evaluations
 * \ingroup iGeom
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup iGeomBounds Spatial Bounding Box
 * \ingroup iGeom
 ******************************************************************************/

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
