#	Fresco: for Fortran
#
#	Ian J. Thompson,     I-Thompson@llnl.gov
#	Dec 2012

FFLAGS = -O2 		# generic
#FFLAGS = -O2 -Wall -ftrace=full	# generic
#FFLAGS = -O2 -Wtabs	# gfortran 
#FFLAGS = -O2  -Wall	# generic
#FFLAGS = -O2 -m32	# 32-bit machines generic
#FFLAGS = -O2 -m64	# 64-bit machines generic

#FFLAGS = -O2 -dcfuns -w				# for nag f95
#FFLAGS = -C=all -C=undefined -dcfuns -w -O0 -gline     # for nag f95 debug
#FFLAGS = -O2 -qmaxmem=16384                            #   optimisation for xlf
#FFLAGS = -Wall -g -fbounds-check -ftrace=full          # debug for g95
#FFLAGS = -g
#FFLAGS = -g  -C  	 				# generic debug for f90
#FFLAGS = -g -check all -O0 -debug 			# debug for ifort
FC=g95
F90C=$(FC)

#               Select subroutine for CPU time:
#TIME = etime.o		# Comment out for g95 versions >= 0.91
#TIME = system_clock.o
TIME = cpu_time.o
#               Select subroutine for flush file buffers to disk
#FLUSH = flush.o
#FLUSH = 
 
COMMON= globx7.o frxx0.o frxx1.o frxx2.o frxx3.o frxx4.o frxx5.o frxx6.o frxx7a.o \
        frxx7b.o frxx8.o coulfg.o frxx9.o frxx10.o frxx11.o naglib.o \
        frxx13.o frxx16.o frxx17.o lapack.o cdc.o   $(TIME) $(FLUSH)

OBJ = $(COMMON) fresco.o 

target: fresco sfresco 
install: fresco sfresco 
	cp -fp fresco sfresco $(HOME)/bin/$(ARCH)

fresco: $(OBJ)
	$(FC) $(FFLAGS) -o fresco $(OBJ)

OBJS = $(COMMON) sfresco.o minuit.o

sfresco: $(OBJS)
	$(FC) $(FFLAGS) -o sfresco $(OBJS)

frxx1.o: globx7.o frxx1.f usescatwf.f veffpot.f bpmfus.f
frxx17.o: globx7.o frxx17.f usescatwf.f

.f.o:
	$(FC) $(FORM) -c $(FFLAGS)  $<

#########
clean:
	\rm -f *.o core sfresco fresco  *.mod
