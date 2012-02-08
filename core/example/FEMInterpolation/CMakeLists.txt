INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  FEMInterpolationExample1
  SOURCES linear_vector.cpp
  COMM serial mpi
  DEPLIBS teuchos intrepid food
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  FEMInterpolationExample2
  SOURCES quadratic_scalar.cpp
  COMM serial mpi
  DEPLIBS teuchos intrepid food
  )
