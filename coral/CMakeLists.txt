project( coral )
cmake_minimum_required( VERSION 2.8.0 )

# These options only need to be set if building this project by itself
if ( NOT rhic_BINARY_DIR )
  set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${coral_SOURCE_DIR}/../CMake )

  set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib )
  set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib )
  set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin )

  ## Need GSL
  find_package( GSL REQUIRED )
endif()

set( coral_INCLUDE_DIRS
  ${coral_SOURCE_DIR}/../b3d/src

  ${coral_SOURCE_DIR}/include

  ${coral_SOURCE_DIR}/src/coral
  ${coral_SOURCE_DIR}/src/coral/CF2Source
  ${coral_SOURCE_DIR}/src/coral/CFCalc
  ${coral_SOURCE_DIR}/src/coral/Kernel
  ${coral_SOURCE_DIR}/src/coral/Source2CF
  ${coral_SOURCE_DIR}/src/coral/SourceCalc
  ${coral_SOURCE_DIR}/src/coral/WaveFunctions

  ${coral_SOURCE_DIR}/src/coralutils
  ${coral_SOURCE_DIR}/src/coralutils/Arrays
  ${coral_SOURCE_DIR}/src/coralutils/GSLMatrix
  ${coral_SOURCE_DIR}/src/coralutils/Misc
  ${coral_SOURCE_DIR}/src/coralutils/NumMath
  ${coral_SOURCE_DIR}/src/coralutils/ParameterMap
  ${coral_SOURCE_DIR}/src/coralutils/SpecialFunctions

  ${GSL_INCLUDE_DIRS}

  ${HDF5_C_INCLUDE_DIR}
)

set( coral_LIBRARY_DIRS
  ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
  ${GSL_LIBRARY_DIRS}
)
set( coral_LIBRARIES
  coral
  coralutils
  xgraph
  ${GSL_LIBRARIES}
)

if ( rhic_BINARY_DIR )
  # Export variables to the parent scope
  set( coral_INCLUDE_DIRS ${coral_INCLUDE_DIRS} PARENT_SCOPE )
  set( coral_LIBRARY_DIRS ${coral_LIBRARY_DIRS} PARENT_SCOPE )
  set( coral_LIBRARIES    ${coral_LIBRARIES}    PARENT_SCOPE )
endif()

add_subdirectory( src )

## Configure the package file
configure_file( coral-config.cmake.in ${PROJECT_BINARY_DIR}/coral-config.cmake @ONLY )
