INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ConsistentInterpolationExample1
  SOURCES example1.cpp
  COMM serial mpi
  DEPLIBS teuchos food_core food_transfer
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ConsistentInterpolationExample2
  SOURCES example2.cpp
  COMM serial mpi
  DEPLIBS teuchos food_core food_transfer
  )
