INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TiedExample1
  SOURCES flat_surf.cpp
  COMM serial mpi
  DEPLIBS teuchos food_core food_transfer
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TiedExample2
  SOURCES curved_surf.cpp
  COMM serial mpi
  DEPLIBS teuchos food_core food_transfer
  )