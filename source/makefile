#	Fresco: for Fortran
#
#	Ian J. Thompson,     I-Thompson@llnl.gov
#	Dec 2012

FFLAGS = -O3 		# generic
#FFLAGS = -O2 -Wall -ftrace=full	# generic
#FFLAGS = -O2  -Wall	# generic
#FFLAGS = -O2 -m32	# 32-bit machines generic
#FFLAGS = -O2 -m64	# 64-bit machines generic

#FFLAGS = -O2 -dcfuns -w				# for nag f95
#FFLAGS = -C=all -C=undefined -dcfuns -w -O0 -gline     # for nag f95 debug
#FFLAGS = -O2 -qmaxmem=16384                            #   optimisation for xlf
#FFLAGS = -g
#FFLAGS = -g  -C  	 				# generic debug for f90
#FFLAGS = -g -check all -O0 -debug 			# debug for ifort
#FFLAGS = -g -fbacktrace -O0 -fcheck=all  		# debug for gfortran
FC=gfortran 
# F90C=$(FC)

#   Use NAG library if available for 2D interpolationm of external 3-body wave functions
# NAG = -lnag
NAG = nagstub.o

#               Select subroutine for CPU time:
#TIME = etime.o		# Comment out for g95 versions >= 0.91
#TIME = system_clock.o
TIME = cpu_time.o
#               Select subroutine for flush file buffers to disk
#FLUSH = flush.o
#FLUSH = 
 
COMMON= globx7.o frxx0.o frxx1.o frxx2.o frxx3.o frxx4.o frxx5.o frxx6.o frxx7a.o \
        frxx7b.o frxx8.o coulfg.o frxx9.o frxx10.o frxx11.o \
        frxx13.o frxx16.o frxx17.o lapack.o cdc.o   $(TIME) $(FLUSH)

OBJ = $(COMMON) fresco.o 

target: fresco sfresco 
install: fresco sfresco 
	cp -fp fresco sfresco $(FRESCOLIB)

#	cp -fp fresco sfresco $(HOME)/bin/$(ARCH)

fresco: $(OBJ)  nagstub.o
	$(FC) $(FFLAGS) -o fresco $(OBJ) $(NAG)

OBJS = $(COMMON) sfresco.o minuit.o

sfresco: $(OBJS)
	$(FC) $(FFLAGS) -o sfresco $(OBJS) $(NAG)

frxx0.o: globx7.o frxx0.f 
frxx1.o: globx7.o frxx1.f usescatwf.f veffpot.f bpmfus.f
frxx2.o: globx7.o frxx2.f 
frxx3.o: globx7.o frxx3.f 
frxx4.o: globx7.o frxx4.f 
frxx5.o: globx7.o frxx5.f 
frxx6.o: globx7.o frxx6.f 
frxx7a.o: globx7.o frxx7a.f 
frxx7b.o: globx7.o frxx7b.f 
frxx8.o: globx7.o frxx8.f 
frxx9.o: globx7.o frxx9.f 
frxx10.o: globx7.o frxx10.f 
frxx11.o: globx7.o frxx11.f 
frxx16.o: globx7.o frxx16.f 
frxx17.o: globx7.o frxx17.f usescatwf.f
cdc.o: globx7.o cdc.f 
coulfg.o: globx7.o coulfg.f 
lanecoup.o: globx7.o lanecoup.f 
sfresco.o: globx7.o sfresco.f 
fresco.o: globx7.o fresco.f 

.f.o:	globx7.o
	$(FC) $(FORM) -c $(FFLAGS)  $<

#########
clean:
	\rm -f *.o core sfresco fresco  *.mod
