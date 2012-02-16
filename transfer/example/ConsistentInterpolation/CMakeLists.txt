INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ConsistentInterpolationExample1
  SOURCES small_box.cpp
  COMM serial mpi
  DEPLIBS teuchos food_core food_transfer
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ConsistentInterpolationExample2
  SOURCES linear_vector.cpp
  COMM serial mpi
  DEPLIBS teuchos food_core food_transfer
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ConsistentInterpolationExample3
  SOURCES quadratic_scalar.cpp
  COMM serial mpi
  DEPLIBS teuchos food_core food_transfer
  )