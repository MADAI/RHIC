project( b3d )
cmake_minimum_required( VERSION 2.8.0 )

# These options only need to be set if building this project by itself
if ( NOT rhic_BINARY_DIR )
  set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib )
  set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib )
  set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin )

  ## Need HDF5
  find_package( HDF5 REQUIRED
    COMPONENTS
      C
      CXX
    )

  ## Need coral library
  find_package( coral REQUIRED )
endif()

set( b3d_INCLUDE_DIRS
  ${b3d_SOURCE_DIR}/src
  ${coral_INCLUDE_DIRS}
  ${HDF5_C_INCLUDE_DIR}
)
set( b3d_LIBRARY_DIRS
  ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
  ${coral_LIBRARY_DIRS}
  ${HDF5_LIBRARY_DIRS}
)
set( b3d_LIBRARIES
  b3d
  ${coral_LIBRARIES}
  ${HDF5_LIBRARIES}
)

if ( rhic_BINARY_DIR )
  # Export variables to the parent scope
  set( b3d_INCLUDE_DIRS ${b3d_INCLUDE_DIRS} PARENT_SCOPE )
  set( b3d_LIBRARY_DIRS ${b3d_LIBRARY_DIRS} PARENT_SCOPE )
  set( b3d_LIBRARIES    ${b3d_LIBRARIES}    PARENT_SCOPE )
endif()

add_subdirectory( src )

## Configure the package file
configure_file( b3d-config.cmake.in ${PROJECT_BINARY_DIR}/b3d-config.cmake @ONLY )
