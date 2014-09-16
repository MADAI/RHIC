# enter make, make install, make clean, make uninstall
# local executables can be made with make runHydro or make multiEos
include ../../makefile_defs.mk
#########################################################################
#HOPEFULLY, everything below need not be touched

INC=-I include -I${MADAI_GSLPATH}/include -I${MADAI_INSTALLDIR}/include/madai -I${MADAI_HDF5_HOME}/include -I${MADAI_GSLPATH}/include

LB=-L${MADAI_INSTALLDIR}/lib/madai -L${MADAI_GSLPATH}/lib -lcoralutils -lgsl -lgslcblas

all : lib/librhicstat.a

statdirs :
	mkdir -p lib;\
	mkdir -p include;\
	mkdir -p bin;\
	mkdir -p build;\
	mkdir -p ${MADAI_INSTALLDIR}/lib/madai;\
	mkdir -p ${MADAI_INSTALLDIR}/include/madai;\
	mkdir -p ${MADAI_INSTALLDIR}/bin/madai

lib/librhicstat.a : build/pca.o build/rhicstat.o build/rhicstat_pca.o build/quadfit.o build/mcmc.o build/mcmc_opt.o build/covariance.o build/gp.o build/locallinear.o build/zgetter.o build/runinfo.o
	rm -f lib/librhicstat.a;\
	ar -ru lib/librhicstat.a build/pca.o build/rhicstat.o build/rhicstat_pca.o build/quadfit.o build/mcmc.o build/mcmc_opt.o build/covariance.o build/gp.o build/locallinear.o build/zgetter.o build/runinfo.o

build/pca.o : src/pca.cc include/rhicstat.h
	${MADAI_CPP} -c src/pca.cc ${MADAI_CFLAGS} ${INC} -o build/pca.o
	
build/rhicstat.o : src/rhicstat.cc include/rhicstat.h
	${MADAI_CPP} -c src/rhicstat.cc ${MADAI_CFLAGS} ${INC} -o build/rhicstat.o

build/rhicstat_pca.o : src/rhicstat_pca.cc include/rhicstat.h
	${MADAI_CPP} -c src/rhicstat_pca.cc ${MADAI_CFLAGS} ${INC} -o build/rhicstat_pca.o

build/quadfit.o : src/quadfit.cc include/rhicstat.h
	${MADAI_CPP} -c src/quadfit.cc ${MADAI_CFLAGS} ${INC} -o build/quadfit.o

build/mcmc.o : src/mcmc.cc include/rhicstat.h
	${MADAI_CPP} -c src/mcmc.cc ${MADAI_CFLAGS} ${INC} -o build/mcmc.o

build/mcmc_opt.o : src/mcmc_opt.cc include/rhicstat.h
	${MADAI_CPP} -c src/mcmc_opt.cc ${MADAI_CFLAGS} ${INC} -o build/mcmc_opt.o

build/covariance.o : src/covariance.cc include/rhicstat.h
	${MADAI_CPP} -c src/covariance.cc ${MADAI_CFLAGS} ${INC} -o build/covariance.o

build/gp.o : src/gp.cc include/rhicstat.h
	${MADAI_CPP} -c src/gp.cc ${MADAI_CFLAGS} ${INC} -o build/gp.o

build/locallinear.o : src/locallinear.cc include/rhicstat.h
	${MADAI_CPP} -c src/locallinear.cc ${MADAI_CFLAGS} ${INC} -o build/locallinear.o

build/zgetter.o : src/zgetter.cc include/rhicstat.h
	${MADAI_CPP} -c src/zgetter.cc ${MADAI_CFLAGS} ${INC} -o build/zgetter.o

build/runinfo.o : src/runinfo.cc include/rhicstat.h
	${MADAI_CPP} -c src/runinfo.cc ${MADAI_CFLAGS} ${INC} -o build/runinfo.o

include/rhicstat.h : src/rhicstat.h
	cp src/rhicstat.h include/rhicstat.h

#########################

install : ${MADAI_INSTALLDIR}/lib/madai/librhicstat.a ${MADAI_INSTALLDIR}/include/madai/rhicstat.h

${MADAI_INSTALLDIR}/lib/madai/librhicstat.a : lib/librhicstat.a
	rm -f ${MADAI_INSTALLDIR}/lib/madai/librhicstat.a;
	cp -f lib/librhicstat.a ${MADAI_INSTALLDIR}/lib/madai/

${MADAI_INSTALLDIR}/include/madai/rhicstat.h : include/rhicstat.h
	cp include/rhicstat.h ${MADAI_INSTALLDIR}/include/madai/

clean :
	rm -f build/* lib/* bin/* include/*

uninstall :
	rm -f ${MADAI_INSTALLDIR}/lib/madai/librhicstat.a ${MADAI_INSTALLDIR}/include/madai/rhicstat.h
