set( b3d_SOURCE
  action.cc
  action_perform.cc
  action_perform_activate.cc
  action_perform_collide.cc
  action_perform_decay.cc
  action_perform_exitcell.cc
  action_perform_vizwrite.cc
  b3d.cc
  bjmaker.cc
  cell.cc
  collide.cc
  decay.cc
  findcollision.cc
  hydrotob3d.cc
  inelastic.cc
  osuhydrotob3d.cc
  part.cc
  resonances.cc
)

include_directories(
  ${b3d_INCLUDE_DIRS}
)

add_library( b3d ${b3d_SOURCE} )

install( TARGETS b3d
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

install( FILES
    b3d.h
    bjmaker.h
    hydrotob3d.h
    inelastic.h
    part.h
    resonances.h
  DESTINATION
    include/b3d
)
