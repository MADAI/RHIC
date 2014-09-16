project( isHydro3 )
cmake_minimum_required( VERSION 2.8.0 )

# These options only need to be set if building this project by itself
if ( NOT rhic_BINARY_DIR )
  set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib )
  set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib )
  set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin )

  ## Need coral library
  find_package( coral REQUIRED )
endif()

set( isHydro3_INCLUDE_DIRS
  ${PROJECT_SOURCE_DIR}/src
  ${coral_INCLUDE_DIRS}
)
set( isHydro3_LIBRARY_DIRS
  ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
  ${coral_LIBRARY_DIRS}
)
set( isHydro3_LIBRARIES
  ishydro3
  ${coral_LIBRARIES}
)

if ( rhic_BINARY_DIR )
  # Export variables to the parent scope
  set( isHydro3_INCLUDE_DIRS ${isHydro3_INCLUDE_DIRS} PARENT_SCOPE )
  set( isHydro3_LIBRARY_DIRS ${isHydro3_LIBRARY_DIRS} PARENT_SCOPE )
  set( isHydro3_LIBRARIES    ${isHydro3_LIBRARIES}    PARENT_SCOPE )
endif()

add_subdirectory( src )

## Configure the package file
configure_file( ishydro3-config.cmake.in ${PROJECT_BINARY_DIR}/ishydro3-config.cmake @ONLY )
