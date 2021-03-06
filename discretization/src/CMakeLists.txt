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
  CellTopologyFactory.hpp
  DFuncKernel.hpp
  DFuncKernelFactory.hpp
  DFuncKernelFactory_Def.hpp
  DiscretizationTypes.hpp
  IntrepidKernel.hpp
  IntrepidKernel_Def.hpp
  IntrepidQuadrature.hpp
  IntrepidQuadrature_Def.hpp
  Quadrature.hpp
  QuadratureFactory.hpp
  QuadratureFactory_Def.hpp
  TopologyTools.hpp
  )

# Set sources

APPEND_SET(SOURCES
  CellTopologyFactory.cpp
  TopologyTools.cpp
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
  food_discretization
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

INSTALL(FILES 
  ${HEADERS} 
  DESTINATION 
  ${CMAKE_INSTALL_PREFIX}/include/food)