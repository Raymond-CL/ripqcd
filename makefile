# executable name
TARGET = ripqcd
# fortran compiler
FC = gfortran

# compiler options (must have)
opt = -fno-strict-overflow
# compiler options for debugging
optd = $(opt) \
	-Wall -Wextra -Wcharacter-truncation -Wunderflow -Wunused-parameter -Walign-commons\
	-Wno-conversion -Wno-intrinsic-shadow -Wno-unused-dummy-argument -Wno-compare-reals\
	-ffixed-form\
	-ffixed-line-length-72\
	-fmodule-private\
	-std=legacy\
	-fcheck=all\
	-fbacktrace
# compiler options for release
optr = $(opt) \
	-finit-local-zero\
	-ffast-math\
	-O3

FFLAGS = $(optd)

# directories
dirpwd = $(PWD)
dirmod = $(dirpwd)/mod
dirsrc = $(dirpwd)/src
dirobj = $(dirpwd)/obj

# fortran file names in order of dependency
mod = \
nrtype.f\
nr.f\
nrutil.f\
ran_state.f\
comvar.f\
sysio.f\
histogram.f\
gauss.f\
glauber.f\
qed.f\
qcd.f\
pdf.f\
ff.f\
sudakov.f

src = \
ran1.f\
vegas.f\
CT18Pdf.f\
EPPS16.f\
akk.f\
hard_dijet.f

prog = main.f

# compile routines
#$(target): $(mod) $(src) $(prog)
#	$(fc) $(opt) -c $^
#	$(fc) *.o -o $(target)
#	rm -rf *.o *.mod

#default: $(target)

debug: $(mod) $(src) $(prog)
	$(FC) $(optd) -c $^
	$(FC) *.o -o $(TARGET)
	rm -rf *.o *.mod
	./$(TARGET)

fast: $(mod) $(src) $(prog)
	$(FC) $(optr) -c $^
	$(FC) *.o -o $(TARGET)
	rm -rf *.o *.mod
#	./$(TARGET)

#.f.o:
#	$(FC) $(FFLAGS) -c $<

#include depend.mk

run:
	./$(TARGET)

clean:
	rm -rf *.o *.mod $(TARGET)

#depend depend.mk:
#	makedepf90 -W -Wconfused -u intrinsic -fixed -o $(TARGET) *.f > depend.mk

# ======================================================================
