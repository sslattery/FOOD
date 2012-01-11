#
# See documentation in Trilinos preCopyrightTrilinos/ExtraExternalRepositories.cmake
#

INCLUDE(TribitsListHelpers)

SET( FOOD_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
  FOOD         .     SS
  )

PACKAGE_DISABLE_ON_PLATFORMS(FOOD Windows)
