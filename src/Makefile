#@(#) Makefile started by makeout(1) 2022-01-17 17:19:25
NULL = 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# platform-specific values you will probably change
F90=gfortran
F90FLAGS := -g -O -std=legacy
LDFLAGS :=  -g -std=legacy
LIBS := 

# You can either ignore it, or disable the "may be used uninitialized"
# warnings with -Wno-maybe-uninitialized or compile with optimizations
# (-O1 and more).
gfortran: F90=gfortran
gfortran: F90FLAGS := -g -O -Wall -std=f2018 -Wextra -Wno-maybe-uninitialized -Wno-uninitialized
gfortran: F90FLAGS := -g -O -std=legacy
gfortran: LDFLAGS :=  -g -std=legacy
gfortran: all
	@echo built with $(F90)

ifort: F90=ifort
ifort: F90FLAGS := -warn all -traceback
ifort: F90FLAGS := -g -O 
ifort: LDFLAGS :=  -g -O
ifort: all
	@echo built with $(F90)

nvfortran: F90=nvfortran
nvfortran: F90FLAGS := -Minform=inform -Mbackslash -traceback
nvfortran: F90FLAGS := -C all -Mbackslash
nvfortran: LDFLAGS :=  -g -O
nvfortran: all
	@echo built with $(F90)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

PROGFILES = \
    ../example/lsoda.f90\
    ../example/lsodar.f90\
    ../example/lsode.f90\
    ../example/lsodes.f90\
    ../example/lsodi.f90\
    ../example/lsodis.f90\
    ../example/lsodkr.f90\
    ../example/lsodpk.f90\
    ../example/lsoibt.f90\
$(NULL)

PROG = ${PROGFILES:.f90=}
PROGO= ${PROGFILES:.f90=.o}
SRCS =  M_odepack.f90
OBJS =  M_odepack.o
all: $(PROG) 
	@echo "That is all folks!"
$(PROG): $(OBJS)
	-$(F90) $(LDFLAGS) $@.f90 -o $@ $(OBJS) $(LIBS)||echo "ouch: $@.f90 "
clean:
	rm -f $(PROG) $(PROGO) $(OBJS) *.mod *.o
.SUFFIXES: $(SUFFIXES) .f90 .F90
.f90.o:
	$(F90) $(F90FLAGS) -c $<
.F90.o:
	$(F90) $(F90FLAGS) -c $<

.PHONY: help
help:
	: 'make [gfortran|nvfortran|ifort]                            '
	: '        -- build odepack and example programs              '
	: 'clean                                                      '
	: '        -- clean directory of object files and executables '
	: 'ford                                                       '
	: '        -- run ford(1) if you have it                      '
	: 'help                                                       '
	: '        -- display this text                               '
#	: 'doxygen -- run doxygen(1) if you have it                   '
#	: 'man     -- show all manpages as text                       '
#	: 'ship    -- rebuild archive files of manpages               '
#	: 'run     -- run manpage demo programs and test program      '
#	: 'test    -- run test program (requires additional modules)  '
#opkdmain.o : M_odepack.o
#opkda1.o   : M_odepack.o
