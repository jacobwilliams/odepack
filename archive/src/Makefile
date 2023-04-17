#@(#) Makefile started by makeout(1) 2022-01-17 17:19:25
NULL = 
# platform-specific values you will probably change
# 
LIBS := 
F90 := nvfortran
F90 := ifort
F90 := gfortran
#  -Wno-aliasing \
#  -Wno-all \
#  -Wno-array-bounds \
#  -Wno-c-binding-type \
#  -Wno-character-truncation \
#  -Wno-conversion \
#  -Wno-function-elimination \
#  -Wno-argument-mismatch \
#  -Wno-implicit-interface \
#  -Wno-implicit-procedure \

F90FLAGS := -O \
  -std=legacy \
$(NULL)

LDFLAGS := \

$(NULL)

PROGFILES = \
 demos/LSODA.f\
 demos/LSODAR.f\
 demos/LSODE.f\
 demos/LSODES.f\
 demos/LSODI.f\
 demos/LSODIS.f\
 demos/LSODKR.f\
 demos/LSODPK.f\
 demos/LSOIBT.f\
$(NULL)

PROG = ${PROGFILES:.f=}
PROGO= ${PROGFILES:.f=.o}
SRCS =  opkda1.f opkdmain.f opkda2.f
OBJS =  opkda1.o opkdmain.o opkda2.o
all: $(PROG) 
	@echo "That is all folks!"
$(PROG): $(OBJS)
	-$(F90) $(LDFLAGS) $@.f -o $@ $(OBJS) $(LIBS)||echo "ouch: $@.f "
clean:
	rm -f $(PROG) $(PROGO) $(OBJS) *.mod
.SUFFIXES: $(SUFFIXES) .f .F
.f.o:
	$(F90) $(F90FLAGS) -c $<
.F.o:
	$(F90) $(F90FLAGS) -c $<
