set( coral_SOURCE
  CF2Source/minimization.cc
  CF2Source/sfit.cc
  CF2Source/sfit_3dgaussian.cc
  CF2Source/sfit_blast.cc
  CF2Source/sfit_GX1d.cc
  CFCalc/cfcalc.cc
  Kernel/kernel.cc
  Source2CF/source2cf.cc
  SourceCalc/sourcecalc.cc
  SourceCalc/sourcecalc_blast.cc
  SourceCalc/sourcecalc_ellipticblast.cc
  SourceCalc/sourcecalc_gauss.cc
  SourceCalc/sourcecalc_getgausspars.cc
  SourceCalc/sourcecalc_GX1d.cc
  SourceCalc/sourcecalc_HDF5.cc
  SourceCalc/sourcecalc_HDF5_multibin.cc
  SourceCalc/sourcecalc_OSCAR.cc
  WaveFunctions/wf/reid93.cc
  WaveFunctions/wf/wf_generic.cc
  WaveFunctions/wf/wf_kpluspiminus_phaseshift.cc
  WaveFunctions/wf/wf_kpluspiminus_sqwell.cc
  WaveFunctions/wf/wf_kpluspiplus_sqwell.cc
  WaveFunctions/wf/wf_lambdalambda_phaseshift.cc
  WaveFunctions/wf/wf_lambdalambdaantiparspin_phaseshift.cc
  WaveFunctions/wf/wf_lambdalambdaparspin_phaseshift.cc
  WaveFunctions/wf/wf_nn_phaseshift.cc
  WaveFunctions/wf/wf_pipluspiminus_nostrong.cc
  WaveFunctions/wf/wf_pipluspiminus_phaseshift.cc
  WaveFunctions/wf/wf_pipluspiminus_sqwell.cc
  WaveFunctions/wf/wf_pipluspiplus_nostrong.cc
  WaveFunctions/wf/wf_pipluspiplus_phaseshift.cc
  WaveFunctions/wf/wf_pipluspiplus_sqwell.cc
  WaveFunctions/wf/wf_pkplus_phaseshift.cc
  WaveFunctions/wf/wf_pkplus_sqwell.cc
  WaveFunctions/wf/wf_plambda_phaseshift.cc
  WaveFunctions/wf/wf_pn_phaseshift.cc
  WaveFunctions/wf/wf_pp_phaseshift.cc
  WaveFunctions/wf/wf_pp_schrod.cc
  WaveFunctions/wf/wf_ppiminus_sqwell.cc
  WaveFunctions/wf/wf_ppiplus_phaseshift.cc
  WaveFunctions/wf/wf_ppiplus_sqwell.cc
  WaveFunctions/wf/wf_Xipi_phaseshift.cc
  WaveFunctions/wfcommon/kpi_phaseshifts.cc
  WaveFunctions/wfcommon/partwave.cc
  WaveFunctions/wfcommon/pipi_phaseshifts.cc
  WaveFunctions/wfcommon/planewave.cc
  WaveFunctions/wfcommon/reid93.cc
  WaveFunctions/wfcommon/wfcommon.cc
  WaveFunctions/wfcommon/wfsquarewell.cc
)

include_directories(
  ${coral_INCLUDE_DIRS}
)

add_library( coral ${coral_SOURCE} )

install( TARGETS coral
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

install( FILES
    coral.h
    CF2Source/minimization.h
    CF2Source/sfit.h
    CFCalc/cfcalc.h
    Kernel/kernel.h
    Source2CF/source2cf.h
    ${coral_SOURCE_DIR}/../b3d/src/part.h
    SourceCalc/sourcecalc.h
    WaveFunctions/wavefunction.h
  DESTINATION
    include/coral
)