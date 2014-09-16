set( coralutils_SOURCE
  Arrays/arraycalc.cc
  Arrays/charray.cc
  Arrays/mclist.cc
  Arrays/mcprlist.cc
  Arrays/threed.cc
  Arrays/ylmarray.cc
  GSLMatrix/gslmatrix.cc
  Misc/misc.cc
  NumMath/random.cc
  ParameterMap/parametermap.cc
  SpecialFunctions/Bessel/bess.cc
  SpecialFunctions/CartesianHarmonics/chcalc.cc
  SpecialFunctions/CoulWave/coulwave.cc
  SpecialFunctions/Ylm/ylm.cc
)

include_directories(
  ${coral_INCLUDE_DIRS}
)

add_library( coralutils ${coralutils_SOURCE} )

install( TARGETS coralutils
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

install( FILES
    coralutils.h
    Arrays/arrays.h
    GSLMatrix/gslmatrix.h
    Misc/misc.h
    NumMath/constants.h
    NumMath/random.h
    ParameterMap/parametermap.h
    SpecialFunctions/Chebyshev/chebyshev.h
    SpecialFunctions/Chebyshev/chebyshev_fit.h
    SpecialFunctions/sf.h
  DESTINATION
    include/coralutils
)
