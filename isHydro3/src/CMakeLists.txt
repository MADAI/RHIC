set( isHydro3_SOURCE
  CCell.cpp
  CEos.cpp
  CHydro.cpp
  CMesh.cpp
  cornelius.cpp
  fullMesh.cpp
  octMesh.cpp
)

include_directories(
  ${isHydro3_INCLUDE_DIRS}
)

add_library( ishydro3
  ${isHydro3_SOURCE}
)

install( TARGETS ishydro3
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

install( FILES
    CCell.h
    CEos.h
    CHydro.h
    CMesh.h
    cornelius.h
    fullMesh.h
    hydroDef.h
    octMesh.h
  DESTINATION
    include/ishydro3
)
