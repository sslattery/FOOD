INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TetPartInterpolateExample
  SOURCES cxx_main.cpp
  COMM serial
  DEPLIBS teuchos food
  )
