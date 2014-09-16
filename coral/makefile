#enter make, then make install, make clean, make uninstall
#other commands: make coralutils, make coral, make install_coralutils, make install_coral
include ../../makefile_defs.mk
#the variables below were set in the include statement above
CPP=${MADAI_CPP} ${MADAI_COMPILERDEFS}
#e.g. /usr/bin/g++
OPT=${MADAI_CFLAGS}
#e.g. -O3
INSTALLDIR=${MADAI_INSTALLDIR}
#e.g. /User/scottepratt/local
GSLPATH=${MADAI_GSLPATH}
HDF5PATH=${MADAI_HDF5_HOME}
X11PATH=${MADAI_X11PATH}
#e.g. /usr/local
#################################################

INC=-I ${GSLPATH}/include -I ${HDF5PATH}/include -I ${X11PATH}/include -I include

all: coraldirs lib/libcoralutils.a lib/libcoral.a

coral : coraldirs lib/libcoralutils.a lib/libcoral.a

coralutils : lib/libcoralutils.a

coraldirs :
	mkdir -p build/coral build/coralutils lib include

CORAL_HFILES = include/minimization.h\
include/sfit.h\
include/cfcalc.h\
include/coral.h\
include/kernel.h\
include/source2cf.h\
include/sourcecalc.h\
include/wavefunction.h\
include/part.h

CORALUTILS_HFILES=include/arrays.h\
include/gslmatrix.h\
include/misc.h\
include/constants.h\
include/random.h\
include/parametermap.h\
include/chebyshev.h\
include/chebyshev_fit.h\
include/sf.h\
include/coralutils.h


CORAL_OBJFILES = build/coral/minimization.o\
build/coral/sfit.o\
build/coral/sfit_3dgaussian.o\
build/coral/sfit_blast.o\
build/coral/sfit_GX1d.o\
build/coral/cfcalc.o\
build/coral/kernel.o\
build/coral/source2cf.o\
build/coral/sourcecalc.o\
build/coral/sourcecalc_getgausspars.o\
build/coral/sourcecalc_ellipticblast.o\
build/coral/sourcecalc_gauss.o\
build/coral/sourcecalc_GX1d.o\
build/coral/sourcecalc_OSCAR.o\
build/coral/sourcecalc_HDF5.o\
build/coral/sourcecalc_HDF5_multibin.o\
build/coral/sourcecalc_OSCAR_multibin.o\
build/coral/sourcecalc_blast.o\
build/coral/wf_generic.o\
build/coral/wf_kpluspiminus_phaseshift.o\
build/coral/wf_kpluspiminus_sqwell.o\
build/coral/wf_kpluspiplus_sqwell.o\
build/coral/wf_lambdalambda_phaseshift.o\
build/coral/wf_lambdalambdaantiparspin_phaseshift.o\
build/coral/wf_lambdalambdaparspin_phaseshift.o\
build/coral/wf_nn_phaseshift.o\
build/coral/wf_pipluspiminus_nostrong.o\
build/coral/wf_pipluspiminus_phaseshift.o\
build/coral/wf_pipluspiminus_sqwell.o\
build/coral/wf_pipluspiplus_nostrong.o\
build/coral/wf_pipluspiplus_phaseshift.o\
build/coral/wf_pipluspiplus_sqwell.o\
build/coral/wf_pkplus_phaseshift.o\
build/coral/wf_pkplus_sqwell.o\
build/coral/wf_plambda_phaseshift.o\
build/coral/wf_pn_phaseshift.o\
build/coral/wf_pp_phaseshift.o\
build/coral/wf_pp_schrod.o\
build/coral/wf_ppiminus_sqwell.o\
build/coral/wf_ppiplus_phaseshift.o\
build/coral/wf_ppiplus_sqwell.o\
build/coral/wf_Xipi_phaseshift.o\
build/coral/kpi_phaseshifts.o\
build/coral/partwave.o\
build/coral/pipi_phaseshifts.o\
build/coral/planewave.o\
build/coral/reid93.o\
build/coral/wfcommon.o\
build/coral/wfsquarewell.o

build/coral/minimization.o : ${CORAL_HFILES} src/coral/CF2Source/minimization.cc
	${CPP} -c ${OPT} ${INC} src/coral/CF2Source/minimization.cc -o build/coral/minimization.o

build/coral/sfit.o : 	${CORAL_HFILES} src/coral/CF2Source/sfit.cc
	${CPP} -c ${OPT} ${INC} src/coral/CF2Source/sfit.cc -o build/coral/sfit.o

build/coral/sfit_3dgaussian.o : 	${CORAL_HFILES} src/coral/CF2Source/sfit_3dgaussian.cc
	${CPP} -c ${OPT} ${INC} src/coral/CF2Source/sfit_3dgaussian.cc -o build/coral/sfit_3dgaussian.o

build/coral/sfit_blast.o : 	${CORAL_HFILES} src/coral/CF2Source/sfit_blast.cc
	${CPP} -c ${OPT} ${INC} src/coral/CF2Source/sfit_blast.cc -o build/coral/sfit_blast.o

build/coral/sfit_GX1d.o : 	${CORAL_HFILES} src/coral/CF2Source/sfit_GX1d.cc
	${CPP} -c ${OPT} ${INC} src/coral/CF2Source/sfit_blast.cc -o build/coral/sfit_GX1d.o

build/coral/cfcalc.o : ${CORAL_HFILES} src/coral/CFCalc/cfcalc.cc
	${CPP} -c ${OPT} ${INC} src/coral/CFCalc/cfcalc.cc -o build/coral/cfcalc.o

build/coral/kernel.o : 	${CORAL_HFILES} src/coral/Kernel/kernel.cc
	${CPP} -c ${OPT} ${INC} src/coral/Kernel/kernel.cc -o build/coral/kernel.o

build/coral/source2cf.o : ${CORAL_HFILES} src/coral/Source2CF/source2cf.cc
	${CPP} -c ${OPT} ${INC} src/coral/Source2CF/source2cf.cc -o build/coral/source2cf.o

build/coral/sourcecalc.o : ${CORAL_HFILES} src/coral/SourceCalc/sourcecalc.cc
	${CPP} -c ${OPT} ${INC} src/coral/SourceCalc/sourcecalc.cc -o build/coral/sourcecalc.o

build/coral/sourcecalc_getgausspars.o : ${CORAL_HFILES} src/coral/SourceCalc/sourcecalc_getgausspars.cc
	${CPP} -c ${OPT} ${INC} src/coral/SourceCalc/sourcecalc_getgausspars.cc -o build/coral/sourcecalc_getgausspars.o

build/coral/sourcecalc_ellipticblast.o : ${CORAL_HFILES} src/coral/SourceCalc/sourcecalc_ellipticblast.cc
	${CPP} -c ${OPT} ${INC} src/coral/SourceCalc/sourcecalc_ellipticblast.cc -o build/coral/sourcecalc_ellipticblast.o

build/coral/sourcecalc_gauss.o : ${CORAL_HFILES} src/coral/SourceCalc/sourcecalc_gauss.cc
	${CPP} -c ${OPT} ${INC} src/coral/SourceCalc/sourcecalc_gauss.cc -o build/coral/sourcecalc_gauss.o

build/coral/sourcecalc_GX1d.o : ${CORAL_HFILES} src/coral/SourceCalc/sourcecalc_GX1d.cc
	${CPP} -c ${OPT} ${INC} src/coral/SourceCalc/sourcecalc_GX1d.cc -o build/coral/sourcecalc_GX1d.o

build/coral/sourcecalc_OSCAR.o : ${CORAL_HFILES} src/coral/SourceCalc/sourcecalc_OSCAR.cc
	${CPP} -c ${OPT} ${INC} src/coral/SourceCalc/sourcecalc_OSCAR.cc -o build/coral/sourcecalc_OSCAR.o

build/coral/sourcecalc_HDF5.o : ${CORAL_HFILES} src/coral/SourceCalc/sourcecalc_HDF5.cc
	${CPP} -c ${OPT} ${INC} src/coral/SourceCalc/sourcecalc_HDF5.cc -o build/coral/sourcecalc_HDF5.o

build/coral/sourcecalc_HDF5_multibin.o : ${CORAL_HFILES} src/coral/SourceCalc/sourcecalc_HDF5_multibin.cc
	${CPP} -c ${OPT} ${INC} src/coral/SourceCalc/sourcecalc_HDF5_multibin.cc -o build/coral/sourcecalc_HDF5_multibin.o
	
build/coral/sourcecalc_OSCAR_multibin.o : ${CORAL_HFILES} src/coral/SourceCalc/sourcecalc_OSCAR_multibin.cc
	${CPP} -c ${OPT} ${INC} src/coral/SourceCalc/sourcecalc_OSCAR_multibin.cc -o build/coral/sourcecalc_OSCAR_multibin.o

build/coral/sourcecalc_blast.o : ${CORAL_HFILES} src/coral/SourceCalc/sourcecalc_blast.cc
	${CPP} -c ${OPT} ${INC} src/coral/SourceCalc/sourcecalc_blast.cc -o build/coral/sourcecalc_blast.o

build/coral/wf_generic.o : ${CORAL_HFILES} src/coral/WaveFunctions/wf/wf_generic.cc
	${CPP} -c ${OPT} ${INC} src/coral/WaveFunctions/wf/wf_generic.cc -o build/coral/wf_generic.o

build/coral/wf_kpluspiminus_phaseshift.o : ${CORAL_HFILES} src/coral/WaveFunctions/wf/wf_kpluspiminus_phaseshift.cc
	${CPP} -c ${OPT} ${INC} src/coral/WaveFunctions/wf/wf_kpluspiminus_phaseshift.cc -o build/coral/wf_kpluspiminus_phaseshift.o

build/coral/wf_kpluspiminus_sqwell.o : ${CORAL_HFILES} src/coral/WaveFunctions/wf/wf_kpluspiminus_sqwell.cc
	${CPP} -c ${OPT} ${INC} src/coral/WaveFunctions/wf/wf_kpluspiminus_sqwell.cc -o build/coral/wf_kpluspiminus_sqwell.o

build/coral/wf_kpluspiplus_sqwell.o : ${CORAL_HFILES} src/coral/WaveFunctions/wf/wf_kpluspiplus_sqwell.cc
	${CPP} -c ${OPT} ${INC} src/coral/WaveFunctions/wf/wf_kpluspiplus_sqwell.cc -o build/coral/wf_kpluspiplus_sqwell.o

build/coral/wf_lambdalambda_phaseshift.o : ${CORAL_HFILES} src/coral/WaveFunctions/wf/wf_lambdalambda_phaseshift.cc
	${CPP} -c ${OPT} ${INC} src/coral/WaveFunctions/wf/wf_lambdalambda_phaseshift.cc -o build/coral/wf_lambdalambda_phaseshift.o

build/coral/wf_lambdalambdaantiparspin_phaseshift.o : ${CORAL_HFILES} src/coral/WaveFunctions/wf/wf_lambdalambdaantiparspin_phaseshift.cc
	${CPP} -c ${OPT} ${INC} src/coral/WaveFunctions/wf/wf_lambdalambdaantiparspin_phaseshift.cc -o build/coral/wf_lambdalambdaantiparspin_phaseshift.o

build/coral/wf_lambdalambdaparspin_phaseshift.o : ${CORAL_HFILES} src/coral/WaveFunctions/wf/wf_lambdalambdaparspin_phaseshift.cc
	${CPP} -c ${OPT} ${INC} src/coral/WaveFunctions/wf/wf_lambdalambdaparspin_phaseshift.cc -o build/coral/wf_lambdalambdaparspin_phaseshift.o

build/coral/wf_nn_phaseshift.o : ${CORAL_HFILES} src/coral/WaveFunctions/wf/wf_nn_phaseshift.cc
	${CPP} -c ${OPT} ${INC} src/coral/WaveFunctions/wf/wf_nn_phaseshift.cc -o build/coral/wf_nn_phaseshift.o

build/coral/wf_pipluspiminus_nostrong.o : ${CORAL_HFILES} src/coral/WaveFunctions/wf/wf_pipluspiminus_nostrong.cc
	${CPP} -c ${OPT} ${INC} src/coral/WaveFunctions/wf/wf_pipluspiminus_nostrong.cc -o build/coral/wf_pipluspiminus_nostrong.o

build/coral/wf_pipluspiminus_phaseshift.o : ${CORAL_HFILES} src/coral/WaveFunctions/wf/wf_pipluspiminus_phaseshift.cc
	${CPP} -c ${OPT} ${INC} src/coral/WaveFunctions/wf/reid93.cc -o build/coral/wf_pipluspiminus_phaseshift.o

build/coral/wf_pipluspiminus_sqwell.o : ${CORAL_HFILES} src/coral/WaveFunctions/wf/wf_pipluspiminus_sqwell.cc
	${CPP} -c ${OPT} ${INC} src/coral/WaveFunctions/wf/wf_pipluspiminus_sqwell.cc -o build/coral/wf_pipluspiminus_sqwell.o

build/coral/wf_pipluspiplus_nostrong.o : ${CORAL_HFILES} src/coral/WaveFunctions/wf/wf_pipluspiplus_nostrong.cc
	${CPP} -c ${OPT} ${INC} src/coral/WaveFunctions/wf/wf_pipluspiplus_nostrong.cc -o build/coral/wf_pipluspiplus_nostrong.o

build/coral/wf_pipluspiplus_phaseshift.o : ${CORAL_HFILES} src/coral/WaveFunctions/wf/wf_pipluspiplus_phaseshift.cc
	${CPP} -c ${OPT} ${INC} src/coral/WaveFunctions/wf/wf_pipluspiplus_phaseshift.cc -o build/coral/wf_pipluspiplus_phaseshift.o

build/coral/wf_pipluspiplus_sqwell.o : ${CORAL_HFILES} src/coral/WaveFunctions/wf/wf_pipluspiplus_sqwell.cc
	${CPP} -c ${OPT} ${INC} src/coral/WaveFunctions/wf/wf_pipluspiplus_sqwell.cc -o build/coral/wf_pipluspiplus_sqwell.o

build/coral/wf_pkplus_phaseshift.o : ${CORAL_HFILES} src/coral/WaveFunctions/wf/wf_pkplus_phaseshift.cc
	${CPP} -c ${OPT} ${INC} src/coral/WaveFunctions/wf/wf_pkplus_phaseshift.cc -o build/coral/wf_pkplus_phaseshift.o

build/coral/wf_pkplus_sqwell.o : ${CORAL_HFILES} src/coral/WaveFunctions/wf/wf_pkplus_sqwell.cc
	${CPP} -c ${OPT} ${INC} src/coral/WaveFunctions/wf/wf_pkplus_sqwell.cc -o build/coral/wf_pkplus_sqwell.o

build/coral/wf_plambda_phaseshift.o : ${CORAL_HFILES} src/coral/WaveFunctions/wf/wf_plambda_phaseshift.cc
	${CPP} -c ${OPT} ${INC} src/coral/WaveFunctions/wf/wf_plambda_phaseshift.cc -o build/coral/wf_plambda_phaseshift.o

build/coral/wf_pn_phaseshift.o : ${CORAL_HFILES} src/coral/WaveFunctions/wf/wf_pn_phaseshift.cc
	${CPP} -c ${OPT} ${INC} src/coral/WaveFunctions/wf/wf_pn_phaseshift.cc -o build/coral/wf_pn_phaseshift.o

build/coral/wf_pp_phaseshift.o : ${CORAL_HFILES} src/coral/WaveFunctions/wf/wf_pp_phaseshift.cc
	${CPP} -c ${OPT} ${INC} src/coral/WaveFunctions/wf/wf_pp_phaseshift.cc -o build/coral/wf_pp_phaseshift.o

build/coral/wf_pp_schrod.o : ${CORAL_HFILES} src/coral/WaveFunctions/wf/wf_pp_schrod.cc
	${CPP} -c ${OPT} ${INC} src/coral/WaveFunctions/wf/wf_pp_schrod.cc -o build/coral/wf_pp_schrod.o

build/coral/wf_ppiminus_sqwell.o : ${CORAL_HFILES} src/coral/WaveFunctions/wf/wf_ppiminus_sqwell.cc
	${CPP} -c ${OPT} ${INC} src/coral/WaveFunctions/wf/wf_ppiminus_sqwell.cc -o build/coral/wf_ppiminus_sqwell.o

build/coral/wf_ppiplus_phaseshift.o : ${CORAL_HFILES} src/coral/WaveFunctions/wf/wf_ppiplus_phaseshift.cc
	${CPP} -c ${OPT} ${INC} src/coral/WaveFunctions/wf/wf_ppiplus_phaseshift.cc -o build/coral/wf_ppiplus_phaseshift.o

build/coral/wf_ppiplus_sqwell.o : ${CORAL_HFILES} src/coral/WaveFunctions/wf/wf_ppiplus_sqwell.cc
	${CPP} -c ${OPT} ${INC} src/coral/WaveFunctions/wf/wf_ppiplus_sqwell.cc -o build/coral/wf_ppiplus_sqwell.o

build/coral/wf_Xipi_phaseshift.o : ${CORAL_HFILES} src/coral/WaveFunctions/wf/wf_Xipi_phaseshift.cc
	${CPP} -c ${OPT} ${INC} src/coral/WaveFunctions/wf/wf_Xipi_phaseshift.cc -o build/coral/wf_Xipi_phaseshift.o

build/coral/kpi_phaseshifts.o : ${CORAL_HFILES} src/coral/WaveFunctions/wfcommon/kpi_phaseshifts.cc
	${CPP} -c ${OPT} ${INC} src/coral/WaveFunctions/wfcommon/kpi_phaseshifts.cc -o build/coral/kpi_phaseshifts.o

build/coral/partwave.o : ${CORAL_HFILES} src/coral/WaveFunctions/wfcommon/partwave.cc
	${CPP} -c ${OPT} ${INC} src/coral/WaveFunctions/wfcommon/partwave.cc -o build/coral/partwave.o

build/coral/planewave.o : ${CORAL_HFILES} src/coral/WaveFunctions/wfcommon/planewave.cc
	${CPP} -c ${OPT} ${INC} src/coral/WaveFunctions/wfcommon/planewave.cc -o build/coral/planewave.o

build/coral/pipi_phaseshifts.o : ${CORAL_HFILES} src/coral/WaveFunctions/wfcommon/pipi_phaseshifts.cc
	${CPP} -c ${OPT} ${INC} src/coral/WaveFunctions/wfcommon/pipi_phaseshifts.cc -o build/coral/pipi_phaseshifts.o

build/coral/reid93.o : ${CORAL_HFILES} src/coral/WaveFunctions/wfcommon/reid93.cc
	${CPP} -c ${OPT} ${INC} src/coral/WaveFunctions/wfcommon/reid93.cc -o build/coral/reid93.o

build/coral/wfcommon.o : ${CORAL_HFILES} src/coral/WaveFunctions/wfcommon/wfcommon.cc
	${CPP} -c ${OPT} ${INC} src/coral/WaveFunctions/wfcommon/wfcommon.cc -o build/coral/wfcommon.o

build/coral/wfsquarewell.o : ${CORAL_HFILES} src/coral/WaveFunctions/wfcommon/wfsquarewell.cc
	${CPP} -c ${OPT} ${INC} src/coral/WaveFunctions/wfcommon/wfsquarewell.cc -o build/coral/wfsquarewell.o

CORALUTILS_OBJFILES = build/coralutils/arraycalc.o\
build/coralutils/charray.o\
build/coralutils/mclist.o\
build/coralutils/mcprlist.o\
build/coralutils/threed.o\
build/coralutils/ylmarray.o\
build/coralutils/gslmatrix.o\
build/coralutils/misc.o\
build/coralutils/random.o\
build/coralutils/parametermap.o\
build/coralutils/bess.o\
build/coralutils/chcalc.o\
build/coralutils/coulwave.o\
build/coralutils/ylm.o

build/coralutils/arraycalc.o : ${CORAL_HFILES} src/coralutils/Arrays/arraycalc.cc
	${CPP} -c ${OPT} ${INC} src/coralutils/Arrays/arraycalc.cc -o build/coralutils/arraycalc.o

build/coralutils/charray.o : ${CORAL_HFILES} src/coralutils/Arrays/charray.cc
	${CPP} -c ${OPT} ${INC} src/coralutils/Arrays/charray.cc -o build/coralutils/charray.o

build/coralutils/mclist.o : ${CORAL_HFILES} src/coralutils/Arrays/mclist.cc
	${CPP} -c ${OPT} ${INC} src/coralutils/Arrays/mclist.cc -o build/coralutils/mclist.o

build/coralutils/mcprlist.o : ${CORAL_HFILES} src/coralutils/Arrays/mcprlist.cc
	${CPP} -c ${OPT} ${INC} src/coralutils/Arrays/mcprlist.cc -o build/coralutils/mcprlist.o

build/coralutils/threed.o : ${CORAL_HFILES} src/coralutils/Arrays/threed.cc
	${CPP} -c ${OPT} ${INC} src/coralutils/Arrays/threed.cc -o build/coralutils/threed.o

build/coralutils/ylmarray.o : ${CORAL_HFILES} src/coralutils/Arrays/ylmarray.cc
	${CPP} -c ${OPT} ${INC} src/coralutils/Arrays/ylmarray.cc -o build/coralutils/ylmarray.o

build/coralutils/gslmatrix.o : ${CORAL_HFILES} src/coralutils/GSLMatrix/gslmatrix.cc
	${CPP} -c ${OPT} ${INC} src/coralutils/GSLMatrix/gslmatrix.cc -o build/coralutils/gslmatrix.o

build/coralutils/misc.o : ${CORAL_HFILES} src/coralutils/Misc/misc.cc
	${CPP} -c ${OPT} ${INC} src/coralutils/Misc/misc.cc -o build/coralutils/misc.o

build/coralutils/random.o : ${CORAL_HFILES} src/coralutils/NumMath/random.cc
	${CPP} -c ${OPT} ${INC} src/coralutils/NumMath/random.cc -o build/coralutils/random.o

build/coralutils/parametermap.o : ${CORAL_HFILES} src/coralutils/ParameterMap/parametermap.cc
	${CPP} -c ${OPT} ${INC} src/coralutils/ParameterMap/parametermap.cc -o build/coralutils/parametermap.o

build/coralutils/bess.o : ${CORAL_HFILES} src/coralutils/SpecialFunctions/Bessel/bess.cc
	${CPP} -c ${OPT} ${INC} src/coralutils/SpecialFunctions/Bessel/bess.cc -o build/coralutils/bess.o

build/coralutils/chcalc.o : ${CORAL_HFILES} src/coralutils/SpecialFunctions/CartesianHarmonics/chcalc.cc
	${CPP} -c ${OPT} ${INC} src/coralutils/SpecialFunctions/CartesianHarmonics/chcalc.cc -o build/coralutils/chcalc.o

build/coralutils/coulwave.o : ${CORAL_HFILES} src/coralutils/SpecialFunctions/CoulWave/coulwave.cc
	${CPP} -c ${OPT} ${INC} src/coralutils/SpecialFunctions/CoulWave/coulwave.cc -o build/coralutils/coulwave.o

build/coralutils/ylm.o : ${CORAL_HFILES} src/coralutils/SpecialFunctions/Ylm/ylm.cc
	${CPP} -c ${OPT} ${INC} src/coralutils/SpecialFunctions/Ylm/ylm.cc -o build/coralutils/ylm.o

include/minimization.h : src/coral/CF2Source/minimization.h
	cp -f src/coral/CF2Source/minimization.h include/

include/sfit.h : src/coral/CF2Source/sfit.h
	cp -f src/coral/CF2Source/sfit.h include/

include/cfcalc.h : src/coral/CFCalc/cfcalc.h
	cp -f src/coral/CFCalc/cfcalc.h include/

include/coral.h : src/coral/coral.h
	cp -f src/coral/coral.h include/

include/kernel.h : src/coral/Kernel/kernel.h
	cp -f src/coral/Kernel/kernel.h include/

include/source2cf.h : src/coral/Source2CF/source2cf.h
	cp -f src/coral/Source2CF/source2cf.h include/

include/sourcecalc.h : src/coral/SourceCalc/sourcecalc.h
	cp -r src/coral/SourceCalc/sourcecalc.h include/

include/wavefunction.h : src/coral/WaveFunctions/wavefunction.h
	cp -f src/coral/WaveFunctions/wavefunction.h include/

include/arrays.h : src/coralutils/Arrays/arrays.h
	cp -f src/coralutils/Arrays/arrays.h include/

include/gslmatrix.h : src/coralutils/GSLMatrix/gslmatrix.h
	cp -f src/coralutils/GSLMatrix/gslmatrix.h include/

include/misc.h : src/coralutils/Misc/misc.h
	cp -f src/coralutils/Misc/misc.h include/

include/constants.h : src/coralutils/NumMath/constants.h
	cp -f src/coralutils/NumMath/constants.h include/

include/random.h : src/coralutils/NumMath/random.h
	cp -f src/coralutils/NumMath/random.h include/

include/parametermap.h : src/coralutils/ParameterMap/parametermap.h
	cp -f src/coralutils/ParameterMap/parametermap.h include/

include/chebyshev.h : src/coralutils/SpecialFunctions/Chebyshev/chebyshev.h
	cp -f src/coralutils/SpecialFunctions/Chebyshev/chebyshev.h include/

include/chebyshev_fit.h : src/coralutils/SpecialFunctions/Chebyshev/chebyshev_fit.h
	cp -f src/coralutils/SpecialFunctions/Chebyshev/chebyshev_fit.h include/

include/sf.h : src/coralutils/SpecialFunctions/sf.h
	cp -f src/coralutils/SpecialFunctions/sf.h include/

include/coralutils.h : src/coralutils/coralutils.h
	cp -f src/coralutils/coralutils.h include/

include/part.h : ${MADAI_HOME}/rhic/b3d2/src/part.h
	cp -f ${MADAI_HOME}/rhic/b3d2/src/part.h include/

lib/libcoralutils.a : ${CORALUTILS_HFILES} ${CORALUTILS_OBJFILES}
	mkdir -p lib;\
	rm -f lib/libcoralutils.a;\
	ar -ru lib/libcoralutils.a ${CORALUTILS_OBJFILES};

lib/libcoral.a : ${CORAL_HFILES} ${CORALUTILS_HFILES} ${CORAL_OBJFILES}
	mkdir -p lib;\
	rm -f lib/libcoral.a;\
	ar -ru lib/libcoral.a ${CORAL_OBJFILES};

install : install_coral

install_coral : install_coralutils ${INSTALLDIR}/lib/madai/libcoral.a ${INSTALLDIR}/include/madai/minimization.h ${INSTALLDIR}/include/madai/sfit.h	${INSTALLDIR}/include/madai/cfcalc.h ${INSTALLDIR}/include/madai/coral.h	${INSTALLDIR}/include/madai/kernel.h 	${INSTALLDIR}/include/madai/source2cf.h 	${INSTALLDIR}/include/madai/kernel.h ${INSTALLDIR}/include/madai/part.h	${INSTALLDIR}/include/madai/source2cf.h 	${INSTALLDIR}/include/madai/sourcecalc.h	${INSTALLDIR}/include/madai/wavefunction.h

install_coralutils : installdirs ${INSTALLDIR}/lib/madai/libcoralutils.a 	${INSTALLDIR}/include/madai/coralutils.h ${INSTALLDIR}/include/madai/arrays.h 	${INSTALLDIR}/include/madai/gslmatrix.h ${INSTALLDIR}/include/madai/misc.h ${INSTALLDIR}/include/madai/constants.h	${INSTALLDIR}/include/madai/random.h 	${INSTALLDIR}/include/madai/parametermap.h 	${INSTALLDIR}/include/madai/chebyshev.h 	${INSTALLDIR}/include/madai/chebyshev_fit.h	${INSTALLDIR}/include/madai/sf.h 	${INSTALLDIR}/include/madai/coralutils.h

installdirs :
	mkdir -p ${INSTALLDIR}/lib/madai;\
	mkdir -p ${INSTALLDIR}/include/madai;\
	mkdir -p ${INSTALLDIR}/bin/madai;\
	mkdir -p ${INSTALLDIR}/progdata/madai;\

${INSTALLDIR}/lib/madai/libcoralutils.a : lib/libcoralutils.a
	cp -f lib/libcoralutils.a ${INSTALLDIR}/lib/madai/

${INSTALLDIR}/lib/madai/libcoral.a : lib/libcoral.a
	cp -f lib/libcoral.a ${INSTALLDIR}/lib/madai/

#########################

${INSTALLDIR}/include/madai/minimization.h : include/minimization.h
	cp -f include/minimization.h ${INSTALLDIR}/include/madai/

${INSTALLDIR}/include/madai/sfit.h : include/sfit.h
	cp -f include/sfit.h ${INSTALLDIR}/include/madai/

${INSTALLDIR}/include/madai/cfcalc.h : include/cfcalc.h
	cp -f include/cfcalc.h ${INSTALLDIR}/include/madai/

${INSTALLDIR}/include/madai/coral.h : include/coral.h
	cp -f include/coral.h ${INSTALLDIR}/include/madai/

${INSTALLDIR}/include/madai/kernel.h : include/kernel.h
	cp -f include/kernel.h ${INSTALLDIR}/include/madai/

${INSTALLDIR}/include/madai/source2cf.h : include/source2cf.h
	cp -f include/source2cf.h ${INSTALLDIR}/include/madai/

${INSTALLDIR}/include/madai/sourcecalc.h : include/sourcecalc.h
	cp -f include/sourcecalc.h ${INSTALLDIR}/include/madai/

${INSTALLDIR}/include/madai/wavefunction.h : include/wavefunction.h
	cp -f include/wavefunction.h ${INSTALLDIR}/include/madai/

${INSTALLDIR}/include/madai/arrays.h : include/arrays.h
	cp -f include/arrays.h ${INSTALLDIR}/include/madai/

${INSTALLDIR}/include/madai/gslmatrix.h : include/gslmatrix.h
	cp -f include/gslmatrix.h ${INSTALLDIR}/include/madai/

${INSTALLDIR}/include/madai/misc.h : include/misc.h
	cp -f include/misc.h ${INSTALLDIR}/include/madai/

${INSTALLDIR}/include/madai/constants.h : include/constants.h
	cp -f include/constants.h ${INSTALLDIR}/include/madai/

${INSTALLDIR}/include/madai/random.h : include/random.h
	cp -f include/random.h ${INSTALLDIR}/include/madai/

${INSTALLDIR}/include/madai/parametermap.h : include/parametermap.h
	cp -f include/parametermap.h ${INSTALLDIR}/include/madai/

${INSTALLDIR}/include/madai/chebyshev.h : include/chebyshev.h
	cp -f include/chebyshev.h ${INSTALLDIR}/include/madai/

${INSTALLDIR}/include/madai/chebyshev_fit.h : include/chebyshev_fit.h
	cp -f include/chebyshev_fit.h ${INSTALLDIR}/include/madai/

${INSTALLDIR}/include/madai/sf.h : include/sf.h
	cp -f include/sf.h ${INSTALLDIR}/include/madai/

${INSTALLDIR}/include/madai/coralutils.h : include/coralutils.h
	cp -f include/coralutils.h ${INSTALLDIR}/include/madai/

${INSTALLDIR}/include/madai/part.h : include/part.h
	cp -f include/part.h ${INSTALLDIR}/include/madai/

######################
uninstall:
	 rm -f ${INSTALLDIR}/lib/madai/libcoral.a ${INSTALLDIR}/lib/madai/libcoralutils.a ${INSTALLDIR}/include/madai/coralutils.h ${INSTALLDIR}/include/madai/coralutils.h ${INSTALLDIR}/include/madai/sf.h ${INSTALLDIR}/include/madai/chebyshev_fit.h ${INSTALLDIR}/include/madai/chebyshev.h ${INSTALLDIR}/include/madai/parametermap.h ${INSTALLDIR}/include/madai/random.h ${INSTALLDIR}/include/madai/constants.h ${INSTALLDIR}/include/madai/misc.h ${INSTALLDIR}/include/madai/gslmatrix.h ${INSTALLDIR}/include/madai/arrays.h ${INSTALLDIR}/include/madai/wavefunction.h ${INSTALLDIR}/include/madai/sourcecalc.h ${INSTALLDIR}/include/madai/source2cf.h ${INSTALLDIR}/include/madai/kernel.h ${INSTALLDIR}/include/madai/coral.h ${INSTALLDIR}/include/madai/cfcalc.h ${INSTALLDIR}/include/madai/sfit.h ${INSTALLDIR}/include/madai/minimization.h ${INSTALLDIR}/include/madai/part.h

clean :
	rm -f ${CORAL_HFILES} ${CORALUTILS_HFILES} ${CORALUTILS_OBJFILES} ${CORAL_OBJFILES} ${CORAL_X} lib/*.a
