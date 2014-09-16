include ../../makefile_defs.mk
# enter make all, make b3d, make install, make clean or make uninstall
### the following variables are set in the include statement above
CPP=${MADAI_CPP}
OPT=${MADAI_CFLAGS}
CORAL_INCLUDE=${INSTALLDIR}/include/madai
#e.g. /Users/scottepratt/git -- 
CORAL_HOME=${MADAI_HOME}/rhic/CorAL/trunk
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

install : installdirs ${INSTALLDIR}/lib/madai/libb3d.a ${INSTALLDIR}/include/madai/b3d.h ${INSTALLDIR}/include/madai/hydrotob3d.h ${INSTALLDIR}/include/madai/resonances.h ${INSTALLDIR}/include/madai/resonances.h ${INSTALLDIR}/include/madai/part.h ${INSTALLDIR}/include/madai/inelastic.h ${INSTALLDIR}/include/madai/bjmaker.h ${INSTALLDIR}/progdata/madai/resinfo/resonances_pdg.dat ${INSTALLDIR}/progdata/madai/resinfo/resonances_pdg_weak.dat ${INSTALLDIR}/progdata/madai/resinfo/decays_pdg.dat ${INSTALLDIR}/progdata/madai/resinfo/decays_pdg_weak.dat ${INSTALLDIR}/progdata/madai/resinfo/pionsonly.dat ${INSTALLDIR}/progdata/madai/resinfo/nodecays.dat

installdirs :
	mkdir -p ${INSTALLDIR}/include/madai;\
	mkdir -p ${INSTALLDIR}/lib/madai;\
	mkdir -p ${INSTALLDIR}/bin/madai;\
	mkdir -p ${INSTALLDIR}/progdata/madai/resinfo;\

clean :
	rm -f lib/libb3d.a build/*.o include/*.h

INC=-I include -I ${GSLPATH}/include -I${CORAL_INCLUDE} -I ${HDF5_HOME}/include

uninstall :
	rm -f ${INSTALLDIR}/lib/madai/libb3d.a ${INSTALLDIR}/include/madai/b3d.h ${INSTALLDIR}/include/madai/hydrotob3d.h ${INSTALLDIR}/include/madai/resonances.h ${INSTALLDIR}/include/madai/bjmaker.h ${INSTALLDIR}/include/madai/inelastic.h \

B3D_HFILES = include/b3d.h\
include/hydrotob3d.h\
include/bjmaker.h\
include/resonances.h\
include/inelastic.h\
include/part.h

B3D_OBJFILES = build/b3d.o\
build/addaction.o\
build/bjmaker.o\
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
build/action_perform_exitcell.o\
build/findcollision.o\
build/hydrotob3d.o\
build/hydrotob3d_PR.o\
build/osuhydrotob3d.o\
build/part.o\
build/inelastic.o\
build/resonances.o

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

include/hydrotob3d.h : src/hydrotob3d.h
	cp src/hydrotob3d.h include/

include/bjmaker.h : src/bjmaker.h
	cp -f src/bjmaker.h include/

include/resonances.h : src/resonances.h
	cp -f src/resonances.h include/

include/part.h : src/part.h
	cp -f src/part.h include/

include/inelastic.h : src/inelastic.h
	cp src/inelastic.h include/

######################

lib/libb3d.a : ${B3D_HFILES} ${B3D_OBJFILES}
	rm -f lib/libb3d.a;\
	ar -ru lib/libb3d.a ${B3D_OBJFILES};\

build/b3d.o : src/b3d.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} -o build/b3d.o src/b3d.cc

build/bjmaker.o : src/bjmaker.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} -o build/bjmaker.o src/bjmaker.cc

build/cell.o : src/cell.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} -o build/cell.o src/cell.cc

build/collide.o : src/collide.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} -o build/collide.o src/collide.cc

build/addaction.o : src/addaction.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} -o build/addaction.o src/addaction.cc

build/scatter.o : src/scatter.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} -o build/scatter.o src/scatter.cc

build/annihilate.o : src/annihilate.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} -o build/annihilate.o src/annihilate.cc

build/decay.o : src/decay.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} -o build/decay.o src/decay.cc

build/action.o : src/action.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} -o build/action.o src/action.cc

build/action_perform.o : src/action_perform.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} -o build/action_perform.o src/action_perform.cc

build/action_perform_activate.o : src/action_perform_activate.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} -o build/action_perform_activate.o src/action_perform_activate.cc

build/action_perform_collide.o : src/action_perform_collide.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} -o build/action_perform_collide.o src/action_perform_collide.cc

build/action_perform_decay.o : src/action_perform_decay.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} -o build/action_perform_decay.o src/action_perform_decay.cc

build/action_perform_exitcell.o : src/action_perform_exitcell.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} -o build/action_perform_exitcell.o src/action_perform_exitcell.cc

build/findcollision.o : src/findcollision.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} -o build/findcollision.o src/findcollision.cc

build/hydrotob3d.o : src/hydrotob3d.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} -o build/hydrotob3d.o src/hydrotob3d.cc

build/hydrotob3d_PR.o : src/hydrotob3d_PR.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} -o build/hydrotob3d_PR.o src/hydrotob3d_PR.cc

build/osuhydrotob3d.o : src/osuhydrotob3d.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} -o build/osuhydrotob3d.o src/osuhydrotob3d.cc

build/part.o : src/part.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} -o build/part.o src/part.cc

build/resonances.o : 	src/resonances.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} -o build/resonances.o src/resonances.cc

build/inelastic.o : src/inelastic.cc ${B3D_HFILES}
	${CPP} -c ${OPT} ${INC} -o build/inelastic.o src/inelastic.cc

############

${INSTALLDIR}/lib/madai/libb3d.a : lib/libb3d.a
	cp -f lib/libb3d.a ${INSTALLDIR}/lib/madai

${INSTALLDIR}/include/madai/b3d.h : include/b3d.h
	cp -f include/b3d.h ${INSTALLDIR}/include/madai/

${INSTALLDIR}/include/madai/hydrotob3d.h : include/hydrotob3d.h
	cp -f include/hydrotob3d.h ${INSTALLDIR}/include/madai/

${INSTALLDIR}/include/madai/resonances.h : include/resonances.h
	cp -f include/resonances.h ${INSTALLDIR}/include/madai/

${INSTALLDIR}/include/madai/part.h : include/part.h
	cp -f include/part.h ${INSTALLDIR}/include/madai/

${INSTALLDIR}/include/madai/bjmaker.h : include/bjmaker.h
	cp -f include/bjmaker.h ${INSTALLDIR}/include/madai/

${INSTALLDIR}/include/madai/inelastic.h : include/inelastic.h
	cp -f include/inelastic.h ${INSTALLDIR}/include/madai/

${INSTALLDIR}/progdata/madai/resinfo/resonances_pdg_new.dat : resinfo/resonances_pdg_new.dat
	cp -f resinfo/resonances_pdg_new.dat ${INSTALLDIR}/progdata/madai/resinfo/

${INSTALLDIR}/progdata/madai/resinfo/resonances_pdg.dat : resinfo/resonances_pdg.dat
	cp -f resinfo/resonances_pdg.dat ${INSTALLDIR}/progdata/madai/resinfo/

${INSTALLDIR}/progdata/madai/resinfo/resonances_pdg_weak.dat : resinfo/resonances_pdg_weak.dat
	cp -f resinfo/resonances_pdg_weak.dat ${INSTALLDIR}/progdata/madai/resinfo/

${INSTALLDIR}/progdata/madai/resinfo/decays_pdg.dat : resinfo/decays_pdg.dat
	cp -f resinfo/decays_pdg.dat ${INSTALLDIR}/progdata/madai/resinfo/

${INSTALLDIR}/progdata/madai/resinfo/decays_pdg_weak.dat : resinfo/decays_pdg_weak.dat
	cp -f resinfo/decays_pdg_weak.dat ${INSTALLDIR}/progdata/madai/resinfo/

${INSTALLDIR}/progdata/madai/resinfo/pionsonly.dat : resinfo/pionsonly.dat
	cp -f resinfo/pionsonly.dat ${INSTALLDIR}/progdata/madai/resinfo/

${INSTALLDIR}/progdata/madai/resinfo/nodecays.dat : resinfo/nodecays.dat
	cp -f resinfo/nodecays.dat ${INSTALLDIR}/progdata/madai/resinfo/
