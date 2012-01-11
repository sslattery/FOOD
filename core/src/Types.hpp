//---------------------------------------------------------------------------//
// \file   Types.hpp
// \author Stuart Slattery
// \brief  Enumerated types for FOOD.
//---------------------------------------------------------------------------//

#ifndef FOOD_TYPES_HPP
#define FOOD_TYPES_HPP

namespace FOOD
{

// Error codes.
enum ErrorType {
    ErrorType_MIN = 0,
    SUCCESS =  ErrorType_MIN,
    MESH_ALREADY_LOADED,
    FILE_NOT_FOUND,
    FILE_WRITE_ERROR,
    NIL_ARRAY,
    BAD_ARRAY_SIZE,
    BAD_ARRAY_DIMENSION,
    INVALID_ENTITY_HANDLE,
    INVALID_ENTITY_COUNT,
    INVALID_ENTITY_TYPE,
    INVALID_ENTITY_TOPOLOGY,
    BAD_TYPE_AND_TOPO,
    ENTITY_CREATION_ERROR,
    INVALID_TAG_HANDLE,
    TAG_NOT_FOUND,
    TAG_ALREADY_EXISTS,
    TAG_IN_USE,
    INVALID_ENTITYSET_HANDLE,
    INVALID_ITERATOR_HANDLE,
    INVALID_ARGUMENT,
    MEMORY_ALLOCATION_FAILED,
    NOT_SUPPORTED,
    FAILURE,
    ErrorType_MAX =  FAILURE
};

// Field precision.
enum Precision {
    Precision_MIN = 0,
    BOOLEAN =  Precision_MIN,
    UCHAR,
    INTEGER,
    FLOAT,
    DOUBLE,
    QUAD,
    Precision_MAX =  QUAD
};

// Algorithm type.
enum AlgType {
    AlgType_MIN = 0,
    LOGICAL =  AlgType_MIN,
    INTEGRAL,
    REAL,
    COMPLEX,
    AlgType_MAX =  COMPLEX
};

// Coordinate type.
enum CoordType {
    CoordType_MIN = 0,
    CARTESIAN =  CoordType_MIN,
    CYLINDRICAL,
    SPHERICAL,
    CoordType_MAX =  SPHERICAL
};

// Storage order.
enum StorageHint {
    StorageHint_MIN = 0,
    BLOCKED =  StorageHint_MIN,
    INTERLEAVED,
    MIXED,
    PER_ENTITY,
    StorageHint_MAX =  PER_ENTITY
};

// Entity type.
enum EntityType {
    EntityType_MIN = 0,
    VERTEX =  EntityType_MIN,
    /**< A topological dimension 0 entity */
    EDGE,
    /**< A topological dimension 1 entity */
    FACE,
    /**< A topological dimension 2 entity */
    REGION,
    /**< A topological dimension 3 entity */
    ALL_TYPES,
    /**< used only in queires to request information about all types */
    EntityType_MAX =  ALL_TYPES
};

// Entity topologies.
enum iMesh_EntityTopology {
    iMesh_EntityTopology_MIN = 0,
    iMesh_POINT = iMesh_EntityTopology_MIN,
        /**< a 0D entity (e.g. a vertex) */
    iMesh_LINE_SEGMENT,
        /**< a 1D entity (e.g. an edge) */
    iMesh_POLYGON,
        /**< a general 2D entity */
    iMesh_TRIANGLE,
        /**< a specific 2D entity bounded by 3 edge entities */
    iMesh_QUADRILATERAL,
        /**< a specific 2D entity bounded by 4 edge entities */
    iMesh_POLYHEDRON,
        /**< a general 3D entity */
    iMesh_TETRAHEDRON,
        /**< a specific 3D entity bounded by 4 triangle entities */
    iMesh_HEXAHEDRON,
        /**< a specific 3D entity bounded by 6 quadrilateral entities */
    iMesh_PRISM,
        /**< a speicifc 3D entity bounded by a combination of 3 quadrilateral
            entities and 2 triangle entities */
    iMesh_PYRAMID,
        /**< a specific 3D entity bounded by a combination of 1 quadrilateral
             entity and 4 triangle entities */
    iMesh_SEPTAHEDRON,
        /**< a hexahedral entity with one collapsed edge */
    iMesh_ALL_TOPOLOGIES,
        /**< used only in queires to request information about all topologies */
    iMesh_EntityTopology_MAX = iMesh_ALL_TOPOLOGIES
};

// Tag types.
enum TagValueType {
    TagValueType_MIN = 0,
    BYTES =  TagValueType_MIN,
    /**< An opaque sequence of bytes, size always measured in bytes */
    INTEGER,
    /**< A value of type \c int */
    DOUBLE,
    /**< A value of type \c double */
    ENTITY_HANDLE,
    /**< A value of type \c  EntityHandle */
    ENTITY_SET_HANDLE,
    /**< A value of type \c  EntitySetHandle */
    TagValueType_MAX =  ENTITY_SET_HANDLE
};


} // end namespace FOOD

#endif // end FOOD_TYPES_HPP

//---------------------------------------------------------------------------//
// end Types.hpp
//---------------------------------------------------------------------------//
