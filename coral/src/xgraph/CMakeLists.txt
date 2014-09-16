set( xgraph_SOURCE
  xgraph.cc
)

find_package( X11 REQUIRED )

add_library( xgraph ${xgraph_SOURCE} )

install( TARGETS xgraph
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

install( TARGETS xgraph
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

install( FILES
    xgraph.h
  DESTINATION
    include/xgraph
)
