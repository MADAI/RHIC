#!/usr/bin/bash
include ../../makefile_defs.mk

# enter make, make parmaker, make install, make clean, make uninstall
# local executables can be made with make runHydro or make multiEos
include ../../makefile_defs.mk
#########################################################################
#HOPEFULLY, everything below need not be touched

INC=-I include -I ${MADAI_GSLPATH}/include -I ${MADAI_INSTALLDIR}/include/madai -I ${MADAI_HDF5_HOME}/include -I ${MADAI_GSLPATH}/include

LB=-L ${MADAI_INSTALLDIR}/lib/madai -L ${MADAI_GSLPATH}/lib -lcoralutils -lgsl -lgslcblas

all : utilsdirs lib/librhicutils.a bin/bmaker

utilsdirs :
	mkdir -p lib;\
	mkdir -p include;\
	mkdir -p bin;\
	mkdir -p build;\
	mkdir -p ${MADAI_INSTALLDIR}/lib/madai;\
	mkdir -p ${MADAI_INSTALLDIR}/include/madai;\
	mkdir -p ${MADAI_INSTALLDIR}/bin/madai

lib/librhicutils.a : build/qualifier.o
	rm -f lib/librhicutils.a;\
	ar -ru lib/librhicutils.a build/qualifier.o

build/qualifier.o : src/qualifier.cc include/qualifier.h
	${MADAI_CPP} -c src/qualifier.cc ${MADAI_CFLAGS} ${INC} -o build/qualifier.o

include/qualifier.h : src/qualifier.h
	cp src/qualifier.h include/qualifier.h

include/latin_center.h : src/latin_center.h
	cp src/latin_center.h include/latin_center.h

#########################

install : utilsdirs ${MADAI_INSTALLDIR}/lib/madai/librhicutils.a ${MADAI_INSTALLDIR}/include/madai/qualifier.h install_bmaker

${MADAI_INSTALLDIR}/lib/madai/librhicutils.a : lib/librhicutils.a
	rm -f ${MADAI_INSTALLDIR}/lib/madai/librhicutils.a;
	cp -f lib/librhicutils.a ${MADAI_INSTALLDIR}/lib/madai/

${MADAI_INSTALLDIR}/include/madai/qualifier.h : include/qualifier.h
	cp include/qualifier.h ${MADAI_INSTALLDIR}/include/madai/

clean :
	rm -f build/* lib/* bin/* include/*

uninstall :
	rm -f ${MADAI_INSTALLDIR}/lib/madai/librhicutils.a ${MADAI_INSTALLDIR}/include/madai/rhicutils.h ${MADAI_INSTALLDIR}/include/madai/qualifier.h ${MADAI_INSTALLDIR}/bin/madai/bmaker


CPP=${MADAI_CPP}
OPT=${MADAI_CFLAGS}

bin/bmaker: src/bmaker.cc
	${CPP} ${OPT} -o bin/bmaker src/bmaker.cc

install_bmaker: bin/bmaker
	cp -f bin/bmaker ${MADAI_INSTALLDIR}/bin/madai/
