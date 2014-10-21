include ../../makefile_defs.mk
# enter make all, make b3d, make install, make clean or make uninstall
### the following variables are set in the include statement above
CPP=${MADAI_CPP}
OPT=${MADAI_CFLAGS}
CORAL_INCLUDE=${INSTALLDIR}/include/madai
#e.g. /Users/scottepratt/git -- 
CORAL_HOME=${MADAI_HOME}/rhic/coral
HDF5_HOME=${MADAI_HDF5_HOME}
#e.g. /opt/local
INSTALLDIR=${MADAI_INSTALLDIR}
#e.g. /Users/scottepratt/local
GSLPATH=${MADAI_GSLPATH}
#e.g. /opt/local
HDF5_HOME=${MADAI_HDF5_HOME}
#e.g. /opt/local
#########################################################################
#HOPEFULLY, everything below need not be touched

all : b3d

#you must be sure coral is installed first
#coral :
#	make -C ${CORAL_HOME} -f ${CORAL_HOME}/makefile;

b3d: b3d_dirs lib/libb3d.a ${B3D_HFILES}

install : installdirs ${INSTALLDIR}/lib/madai/libb3d.a ${INSTALLDIR}/include/madai/b3d.h ${INSTALLDIR}/include/madai/sampler.h ${INSTALLDIR}/include/madai/resonances.h ${INSTALLDIR}/include/madai/resonances.h ${INSTALLDIR}/include/madai/part.h ${INSTALLDIR}/include/madai/inelastic.h ${INSTALLDIR}/include/madai/joshconvert.h ${INSTALLDIR}/include/madai/hydrotob3d.h ${INSTALLDIR}/include/madai/mutinfo.h ${INSTALLDIR}/include/madai/seinfo.h  ${INSTALLDIR}/include/madai/cell.h ${INSTALLDIR}/include/madai/regen.h ${INSTALLDIR}/progdata/madai/resinfo/resonances_pdg.dat ${INSTALLDIR}/progdata/madai/resinfo/resonances_pdg_weak.dat ${INSTALLDIR}/progdata/madai/resinfo/decays_pdg.dat ${INSTALLDIR}/progdata/madai/resinfo/decays_pdg_weak.dat ${INSTALLDIR}/progdata/madai/resinfo/resonances_nodecays.dat ${INSTALLDIR}/progdata/madai/resinfo/pionsonly.dat ${INSTALLDIR}/progdata/madai/resinfo/decays_none.dat

installdirs :
	mkdir -p ${INSTALLDIR}/include/madai;\
	mkdir -p ${INSTALLDIR}/lib/madai;\
	mkdir -p ${INSTALLDIR}/bin/madai;\
	mkdir -p ${INSTALLDIR}/progdata/madai/resinfo;\

clean :
	rm -f build/*.o include/*.h lib/*.a

INC=-I include -I ${GSLPATH}/include -I ${CORAL_INCLUDE} -I ${HDF5_HOME}/include

uninstall :
	rm -f ${INSTALLDIR}/lib/madai/libb3d.a ${INSTALLDIR}/include/madai/b3d.h ${INSTALLDIR}/include/madai/sampler.h ${INSTALLDIR}/include/madai/resonances.h ${INSTALLDIR}/include/madai/inelastic.h ${INSTALLDIR}/include/madai/joshconvert.h ${INSTALLDIR}/include/madai/hydrotob3d.h

B3D_HFILES = include/b3d.h\
include/sampler.h\
include/resonances.h\
include/inelastic.h\
include/part.h\
include/joshconvert.h\
include/hydrotob3d.h\
include/mutinfo.h\
include/seinfo.h\
include/regen.h\
include/cell.h

B3D_OBJFILES = build/b3d.o\
build/b3dsubs.o\
build/b3dio.o\
build/addaction.o\
build/cell.o\
build/collide.o\
build/scatter.o\
build/annihilate.o\
build/decay.o\
build/action.o\
build/action_perform.o\
build/action_perform_activate.o\
build/action_perform_collide.o\
build/action_perform_decay.o\
build/action_perform_regen.o\
build/action_perform_exitcell.o\
build/findcollision.o\
build/sampler.o\
build/mutinfo.o\
build/seinfo.o\
build/regen.o\
build/volume_element.o\
build/part.o\
build/inelastic.o\
build/resonances.o\
build/anal_spectra.o\
build/anal_balance.o\
build/anal_hbt.o\
build/anal_v2.o\
build/consolidate.o\
build/anal_realitydiff.o\
build/joshconvert.o\
build/hydrotob3d.o

#######################

b3d_dirs :
	mkdir -p build;\
	mkdir -p include;\
	mkdir -p lib;\
	mkdir -p ${INSTALLDIR}/include;\
	mkdir -p ${INSTALLDIR}/lib;\
	mkdir -p ${INSTALLDIR}/progdata
	mkdir -p ${INSTALLDIR}/progdata/madai/resinfo

#######################

include/b3d.h : src/b3d.h
	cp -f src/b3d.h include/

include/sampler.h : src/sampler.h
	cp src/sampler.h include/

include/resonances.h : src/resonances.h
	cp -f src/resonances.h include/

include/part.h : src/part.h
	cp -f src/part.h include/

include/inelastic.h : src/inelastic.h
	cp src/inelastic.h include/

include/joshconvert.h : src/joshconvert.h
	cp src/joshconvert.h include/
	
include/hydrotob3d.h : src/hydrotob3d.h
	cp src/hydrotob3d.h include/
	
include/mutinfo.h : src/mutinfo.h
	cp src/mutinfo.h include/

include/seinfo.h : src/seinfo.h
	cp src/seinfo.h include/
	
include/cell.h : src/cell.h
	cp src/cell.h include/

include/regen.h : src/regen.h
	cp src/regen.h include/

######################

lib/libb3d.a : ${B3D_HFILES} ${B3D_OBJFILES}
	rm -f lib/libb3d.a;\
	ar -ru lib/libb3d.a ${B3D_OBJFILES};

build/b3d.o : src/b3d.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/b3d.o src/b3d.cc

build/b3dsubs.o : src/b3dsubs.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/b3dsubs.o src/b3dsubs.cc
	
build/b3dio.o : src/b3dio.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/b3dio.o src/b3dio.cc

build/cell.o : src/cell.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/cell.o src/cell.cc

build/collide.o : src/collide.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/collide.o src/collide.cc

build/addaction.o : src/addaction.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/addaction.o src/addaction.cc

build/scatter.o : src/scatter.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/scatter.o src/scatter.cc

build/annihilate.o : src/annihilate.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/annihilate.o src/annihilate.cc

build/decay.o : src/decay.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/decay.o src/decay.cc

build/action.o : src/action.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/action.o src/action.cc

build/action_perform.o : src/action_perform.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/action_perform.o src/action_perform.cc

build/action_perform_activate.o : src/action_perform_activate.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/action_perform_activate.o src/action_perform_activate.cc

build/action_perform_collide.o : src/action_perform_collide.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/action_perform_collide.o src/action_perform_collide.cc

build/action_perform_decay.o : src/action_perform_decay.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/action_perform_decay.o src/action_perform_decay.cc

build/action_perform_regen.o : src/action_perform_regen.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/action_perform_regen.o src/action_perform_regen.cc

build/action_perform_exitcell.o : src/action_perform_exitcell.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/action_perform_exitcell.o src/action_perform_exitcell.cc

build/findcollision.o : src/findcollision.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/findcollision.o src/findcollision.cc

build/sampler.o : src/sampler.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/sampler.o src/sampler.cc

build/mutinfo.o : src/mutinfo.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/mutinfo.o src/mutinfo.cc

build/seinfo.o : src/seinfo.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/seinfo.o src/seinfo.cc
	
build/regen.o : src/regen.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/regen.o src/regen.cc

build/volume_element.o : src/volume_element.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/volume_element.o src/volume_element.cc

build/part.o : src/part.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/part.o src/part.cc

build/resonances.o : 	src/resonances.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/resonances.o src/resonances.cc

build/inelastic.o : src/inelastic.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/inelastic.o src/inelastic.cc
	
build/anal_spectra.o : src/anal_spectra.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/anal_spectra.o src/anal_spectra.cc

build/anal_balance.o : src/anal_balance.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/anal_balance.o src/anal_balance.cc
	
build/anal_hbt.o : src/anal_hbt.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/anal_hbt.o src/anal_hbt.cc

build/anal_v2.o : src/anal_v2.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/anal_v2.o src/anal_v2.cc

build/anal_realitydiff.o : src/anal_realitydiff.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/anal_realitydiff.o src/anal_realitydiff.cc
		
build/consolidate.o : src/consolidate.cc ${B3D_HFIES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/consolidate.o src/consolidate.cc
		
build/joshconvert.o : src/joshconvert.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/joshconvert.o src/joshconvert.cc

build/hydrotob3d.o : src/hydrotob3d.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} ${MADAI_COMPILERDEFS} -o build/hydrotob3d.o src/hydrotob3d.cc
	
############

${INSTALLDIR}/lib/madai/libb3d.a : lib/libb3d.a
	cp -f lib/libb3d.a ${INSTALLDIR}/lib/madai

${INSTALLDIR}/include/madai/b3d.h : include/b3d.h
	cp -f include/b3d.h ${INSTALLDIR}/include/madai/

${INSTALLDIR}/include/madai/sampler.h : include/sampler.h
	cp -f include/sampler.h ${INSTALLDIR}/include/madai/

${INSTALLDIR}/include/madai/resonances.h : include/resonances.h
	cp -f include/resonances.h ${INSTALLDIR}/include/madai/

${INSTALLDIR}/include/madai/part.h : include/part.h
	cp -f include/part.h ${INSTALLDIR}/include/madai/

${INSTALLDIR}/include/madai/inelastic.h : include/inelastic.h
	cp -f include/inelastic.h ${INSTALLDIR}/include/madai/

${INSTALLDIR}/include/madai/joshconvert.h : include/joshconvert.h
	cp -f include/joshconvert.h ${INSTALLDIR}/include/madai/
	
${INSTALLDIR}/include/madai/hydrotob3d.h : include/hydrotob3d.h
	cp -f include/hydrotob3d.h ${INSTALLDIR}/include/madai/
	
${INSTALLDIR}/include/madai/mutinfo.h : include/mutinfo.h
	cp -f include/mutinfo.h ${INSTALLDIR}/include/madai/

${INSTALLDIR}/include/madai/seinfo.h : include/seinfo.h
	cp -f include/seinfo.h ${INSTALLDIR}/include/madai/
	
${INSTALLDIR}/include/madai/regen.h : include/regen.h
	cp -f include/regen.h ${INSTALLDIR}/include/madai/
	
${INSTALLDIR}/include/madai/cell.h : include/cell.h
	cp -f include/cell.h ${INSTALLDIR}/include/madai/

${INSTALLDIR}/progdata/madai/resinfo/resonances_pdg_new.dat : resinfo/resonances_pdg_new.dat
	cp -f resinfo/resonances_pdg_new.dat ${INSTALLDIR}/progdata/madai/resinfo/

${INSTALLDIR}/progdata/madai/resinfo/resonances_pdg.dat : resinfo/resonances_pdg.dat
	cp -f resinfo/resonances_pdg.dat ${INSTALLDIR}/progdata/madai/resinfo/

${INSTALLDIR}/progdata/madai/resinfo/resonances_nodecays.dat : resinfo/resonances_nodecays.dat
	cp -f resinfo/resonances_nodecays.dat ${INSTALLDIR}/progdata/madai/resinfo/

${INSTALLDIR}/progdata/madai/resinfo/resonances_pdg_weak.dat : resinfo/resonances_pdg_weak.dat
	cp -f resinfo/resonances_pdg_weak.dat ${INSTALLDIR}/progdata/madai/resinfo/

${INSTALLDIR}/progdata/madai/resinfo/pionsonly.dat : resinfo/pionsonly.dat
	cp -f resinfo/pionsonly.dat ${INSTALLDIR}/progdata/madai/resinfo/

${INSTALLDIR}/progdata/madai/resinfo/decays_pdg.dat : resinfo/decays_pdg.dat
	cp -f resinfo/decays_pdg.dat ${INSTALLDIR}/progdata/madai/resinfo/

${INSTALLDIR}/progdata/madai/resinfo/decays_pdg_weak.dat : resinfo/decays_pdg_weak.dat
	cp -f resinfo/decays_pdg_weak.dat ${INSTALLDIR}/progdata/madai/resinfo/

${INSTALLDIR}/progdata/madai/resinfo/decays_none.dat : resinfo/decays_none.dat
	cp -f resinfo/decays_none.dat ${INSTALLDIR}/progdata/madai/resinfo/

joshtest : b3d joshtest.cc
	clang ${MADAI_CFLAGS} joshtest.cc -Iinclude -I${MADAI_GSLPATH}/lib/ -I/Users/scottepratt/git/rhic/coral/include -Llib -L${MADAI_GSLPATH}/lib -L/Users/scottepratt/git/rhic/coral/lib -lcoralutils -lb3d -lgsl -lgslcblas -o joshtest
