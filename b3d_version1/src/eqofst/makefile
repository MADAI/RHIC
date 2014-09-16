include ../../../makefile_defs.mk
OPT=${MADAI_CFLAGS}
CH=${MADAI_H}
#
CORALINCLUDE=${MADAI_INSTALLDIR}/include
CORALLIB=${MADAI_INSTALLDIR}/lib
GSLLIB=${MADAI_GSLPATH}/lib
GSLINC=${MADAI_GSLPATH}/include

HFILES=eqofst.h species.h intrinsic.h

lib/libeqofst.a : eqofst.o eqofst_threephase.o intrinsic.o species.o
	ar r lib/libeqofst.a eqofst.o eqofst_threephase.o intrinsic.o species.o

eqofst.o : ${HFILES} eqofst.cc
	g++ ${OPT} -I ${CORALINCLUDE} -L${GSLLIB} -I${GSLINC} -c eqofst.cc -o eqofst.o
	
eqofst_threephase.o : ${HFILES} eqofst.cc
	g++ ${OPT} -I ${CORALINCLUDE} -L${GSLLIB} -I${GSLINC} -c eqofst_threephase.cc -o eqofst_threephase.o

intrinsic.o : ${HFILES} eqofst.cc
	g++ ${OPT} -I ${CORALINCLUDE} -L${GSLLIB} -I${GSLINC} -c intrinsic.cc -o intrinsic.o
	
species.o : ${HFILES} species.cc
	g++ ${OPT} -I ${CORALINCLUDE} -L${GSLLIB} -I${GSLINC} -c species.cc -o species.o
	
eqofst_test : eqofst_test.cc lib/libeqofst.a
	g++ -o eqofst_test ${OPT} eqofst_test.cc -I. -I${CORALINCLUDE} -I${GSLINC} -L${GSLLIB} -L./lib -L${CORALLIB} -leqofst -lcoralutils  -lgsl -lgslcblas

interpolator.o : interpolator.cc interpolator.h
	g++ -c ${OPT} interpolator.cc -I . -I ${CORALINCLUDE} -I ${MSU_HOME}/3dhydro -I ${GSLINC}  -o interpolator.o

interpolator_test : interpolator_test.cc interpolator.o ./lib/libeqofst.a
	g++ -o interpolator_test ${OPT} interpolator_test.cc  -I. -I${CORALINCLUDE} -I${GSLINC} -L./lib -L${CORALLIB} -L${GSLLIB} interpolator.o -leqofst -lcoralutils  -lgsl -lgslcblas
