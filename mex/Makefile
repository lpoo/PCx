
ROOTDIR   = ..

SRCDIR    = ${ROOTDIR}/SRC
MEXDIR    = ${ROOTDIR}/mex

TARGETDIR = ${MEXDIR}

NG_LIB    = ../Ng-Peyton/cholesky.a
PCX_LIB   = libPCx.a

# f77 under Solaris does not need the f2c libraries, but g77 does....
# if you have the f2c libraries on your system somewhere then just define
# F2C_LIB appropriately and you're done.
# let me know if you have a more elegant solution!

ifeq ($(FC),g77)
F2C_LIB   = -L../F2C -lf2c
F2C       = f2c
endif

CFLAGS = -O

MEX = mex
MEXDEBUG = -O
MEXOPTS = ${MEXDEBUG} -v -I${SRCDIR} 

mex: ${F2C} ${SRCDIR}/blkLVL.o ${SRCDIR}/Ng-Peyton.o ${TARGETDIR}/${PCX_LIB} 
	$(MEX) $(MEXOPTS) -output PCx PCx_mex.c matPCx_interface.c \
	${SRCDIR}/blkLVL.o ${SRCDIR}/Ng-Peyton.o $(PCX_LIB) \
	$(MISC_OBJS) $(NG_LIB) $(F2C_LIB)

${TARGETDIR}/${PCX_LIB}:
	cd ${SRCDIR}; touch *.c; make TARGETDIR=${TARGETDIR} MEX=MEX objs
	ar r $(PCX_LIB) *.o
	ranlib $(PCX_LIB)

f2c:
	cd ../F2C; make

${SRCDIR}/blkLVL.o:
	cd ${SRCDIR}; make TARGETDIR=./ blkLVL.o

${SRCDIR}/Ng-Peyton.o:
	cd ${SRCDIR}; make TARGETDIR=./ Ng-Peyton.o

clean:
	/bin/rm -f *.mex* *.o *.a

