
OBJS = presolve.o PCx2.o readmps.o lpmps.o memory.o hash.o \
       parameters.o split.o io.o timers.o PCx.o solve.o \
       dcolumns.o wrappers.o basics.o scale.o lp.o \
       jair.o cblas.o algorithm.o

PCX_LIB   =  libPCx.a
MYSOLVER_LIB = ../mysolver/libmysolver.a

#NG_LIB    is defined in the build script
#WSSMP_LIB is defined in the build script

FORTRAN = ${FC}

ifdef MEX
CFLAGS = -O -D$(PCx_ARCH) -DMEX $(MORECFLAGS)
#CFLAGS = -g -D$(PCx_ARCH) -DMEX
else
CFLAGS = -O -D$(PCx_ARCH) $(MORECFLAGS)
#CFLAGS = -g -D$(PCx_ARCH)
endif

FFLAGS = -O $(MOREFFLAGS)
LFLAGS = -O

include ../MAKEARCH/$(PCx_ARCH).make
# MISC_OBJS is defined in architecture-specific makefile 

# this is still the default:
PCx_NgPeyton: main.o $(PCX_LIB) blkLVL.o Ng-Peyton.o $(MISC_OBJS) $(NG_LIB)
	$(FORTRAN) $(LFLAGS) main.o blkLVL.o Ng-Peyton.o ${TARGETDIR}/$(PCX_LIB) \
	  $(MISC_OBJS) $(NG_LIB) -o ../PCx -lm

# might have to add the parameter -bmaxdata:44000000 for large applications.
# (or with en even larger number)
PCx_wssmp: main.o $(PCX_LIB) wssmp.o $(MISC_OBJS) $(WSSMP_LIB)
	$(FORTRAN) $(LFLAGS) main.o $(PCX_LIB) wssmp.o $(MISC_OBJS) \
	  $(WSSMP_LIB) -o ../PCx -lm

# template for other solvers
PCx_mysolver: main.o $(PCX_LIB) mysolver.o $(MISC_OBJS) $(MYSOLVER_LIB)
	$(CC) $(LFLAGS) main.o $(PCX_LIB) mysolver.o $(MISC_OBJS) \
	  $(MYSOLVER_LIB) -o ../PCx -lm

objs: $(OBJS)

clean:
	/bin/rm -f *.o $(PCX_LIB) ../PCx

$(PCX_LIB): $(OBJS)
	cd ${TARGETDIR}; ar r $(PCX_LIB) $(OBJS)
	$(RANLIB) $(PCX_LIB)

.c.o:
	${CC} ${CFLAGS} -c -o ${TARGETDIR}/$*.o $<

.f.o:
	${FC} ${FFLAGS} -c -o ${TARGETDIR}/$*.o $<

