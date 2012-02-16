# CMakeLists.txt for src

INCLUDE(TribitsLibraryMacros)

#
# A) Package-specific configuration options
#

#
# B) Define the header and source files (and include directories)
#

SET(HEADERS "")
SET(SOURCES "")

# Set headers (for IDEs and installing)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  BasisFactory.hpp
  BasisFactory_Def.hpp
  CellTopologyFactory.hpp
  DFuncKernel.hpp
  DFuncKernel_Def.hpp
  Domain.hpp
  FEMInterpolate.hpp
  FEMInterpolate_Def.hpp
  KDTree.hpp
  KDTree_Def.hpp
  KernelTraits.hpp
  Octree.hpp
  PointQuery.hpp
  Quantity.hpp
  Rendezvous.hpp
  Rendezvous_Def.hpp
  TensorField.hpp
  TensorField_Def.hpp
  TensorTemplate.hpp
  TopologyTools.hpp
  Types.hpp
  TypeTraits.hpp
  Unit.hpp
  )

# Set sources

APPEND_SET(SOURCES
  CellTopologyFactory.cpp
  Domain.cpp
  KernelTraits.cpp
  Octree.cpp
  PointQuery.cpp
  Quantity.cpp
  TensorTemplate.cpp
  TopologyTools.cpp
  Unit.cpp
  )

# Set addition headers and source based on package-specific
# configuration options

IF(${PACKAGE_NAME}_ENABLE_EXTENDED)
APPEND_SET(HEADERS
    )

  APPEND_SET(SOURCES
    )

ENDIF()

#
# C) Define the targets for install
#

TRIBITS_ADD_LIBRARY(
  food
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

INSTALL(FILES 
  ${HEADERS} 
  DESTINATION 
  ${CMAKE_INSTALL_PREFIX}/include/)