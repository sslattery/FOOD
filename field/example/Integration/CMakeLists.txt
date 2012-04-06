INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  IntegrationExample1
  SOURCES example1.cpp
  COMM serial
  DEPLIBS teuchos food_discretization food_field
  )

