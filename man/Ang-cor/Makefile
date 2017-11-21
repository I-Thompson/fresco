#
# For g77:
#   FC = g77
#   FFLAGS = -fno-f90 -fugly-complex
#
# For gfortran:
#   FC = gfortran
#   FFLAGS =
#

FC = gfortran
FFLAGS = -Wall -fbounds-check
#LDFLAGS = -L/usr/lib/cernlib/2005/lib/ -lkernlib -lpacklib -lmathlib
LDFLAGS = -L/usr/lib64/cernlib/2006/lib -lkernlib -lpacklib -lmathlib

ang-cor: ang-cor.f ang-cor.inc
	$(FC) $(FFLAGS) $(LDFLAGS) -o ang-cor ang-cor.f

.PHONY: clean
clean:
	rm ang-cor.dsdo
	rm ang-cor.out

