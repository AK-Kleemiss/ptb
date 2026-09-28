PROG = ../build/ptb
WITH_DFTD4 = TRUE

SOURCES = accuracy.f90 timer.f90 la.f90 blas_level1.f90 blas_level2.f90 blas_level3.f90 lapack_eig.f90
$(foreach src, $(SOURCES), $(eval $(src) := $(patsubst %.f90, %.o, $(src))))

MODULES= bascom.o parcom.o com.o thresholds.o cbascom.o mocom.o dftd4.o aescom.o default_files.o \
			$(foreach src, $(SOURCES), $($(src)))

CFILES= symmetry_i.o

OBJS1 = matrix_ops.o metrics.o purification_settings.o norms.o  purification.o \
		rdcoord.o readl.o asym.o elem.o readl2.o readline.o prmat.o  timing.o echo.o solve.o rdstvint.o \
        ints.o intpack.o pop.o rdmat.o fermi.o reorder.o bas.o occ.o dtrf2.o lin.o energy.o \
        cma.o drsp.o mlpop.o calcrab.o pgtb.o eeq.o ncoord.o prdipole.o \
        sscal.o cbas.o gauss.o printmos.o qsort.o rdtm.o shellocc_ref.o wr_control_atoms.o \
        guess_qsh.o wrdenmat.o

MAINFILES= main.o 

OBJS = $(MODULES) $(CFILES) $(OBJS1) $(MAINFILES)
#--------------------------------------------------------------------------
ifeq ($(OS), Windows_NT)
  NAME := WINDOWS
  RM = rm
else
  UNAME_S := $(shell uname -s)
  ifeq ($(UNAME_S), Linux)
    NAME := LINUX
	RM = rm
	COMPILER ?= ifx
  endif
  ifeq ($(UNAME_S), Darwin)
    NAME := MAC
  endif
endif
#--------------------------------------------------------------------------
FC = ifx
CC  = icx
ifeq ($(NAME),WINDOWS)
	LINKER = $(FC) -static -I"$(MKLROOT)/include/intel64/lp64" -I"$(MKLROOT)/include" -Qiopenmp
	LIBS = \
	"$(MKLROOT)/lib/mkl_blas95_lp64.lib" \
	"$(MKLROOT)/lib/mkl_lapack95_lp64.lib" \
	"$(MKLROOT)/lib/mkl_intel_lp64.lib" \
	"$(MKLROOT)/lib/mkl_core.lib" \
	"$(MKLROOT)\lib\mkl_intel_thread.lib" 
#	"$(MKLROOT)\..\..\2024.0\lib\libiomp5md.lib"
	FFLAGS = -Od -traceback -Qiopenmp -nologo
	CCFLAGS = -O -DLINUX
endif
ifeq ($(NAME),LINUX)
	ifeq ($(COMPILER),ifx)
	    LINKER = $(FC) -static -I"$(MKLROOT)/include/intel64/lp64" -I"$(MKLROOT)/include" -O -fopenmp
	    LIBS = \
	    $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a \
	    $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a \
	    -Wl,--start-group \
	    $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
	    $(MKLROOT)/lib/intel64/libmkl_core.a \
	    $(MKLROOT)/lib/intel64/libmkl_intel_thread.a \
	    -Wl,--end-group -lpthread -lm
	    FFLAGS = -Og -g -align array64byte -fopenmp -qopenmp -traceback -nologo
	    CCFLAGS = -O -DLINUX
	else
	    FC = gfortran
	    CC = gcc
	    # MKLROOT is expected to be set in the environment (e.g. via
	    # `module load intel-oneapi-mkl` on the cluster); dynamic link
	    # against the single MKL dispatch lib (mkl_rt) so we don't need
	    # to hardcode interface/threading-layer static lib names, which
	    # differ across MKL layouts. Both candidate lib dirs are passed;
	    # -L on a nonexistent path is not fatal.
	    LINKER = $(FC) -I"$(MKLROOT)/include" -O3 -fopenmp
	    LIBS = \
		-L"$(MKLROOT)/lib" -L"$(MKLROOT)/lib/intel64" -Wl,-rpath,$(MKLROOT)/lib -Wl,-rpath,$(MKLROOT)/lib/intel64 \
		-lmkl_rt -lpthread -lm -ldl
	    # Production fallback when ifx is unavailable.  The fit launches many
	    # independent one-thread PTB processes, so runtime checks and an
	    # unoptimised build make the array jobs unnecessarily slow.
	    FFLAGS = -O3 -ffree-line-length-none -fopenmp
	    ifeq ($(DEBUG),1)
	      FFLAGS = -Og -g -fbacktrace -fcheck=all -ffree-line-length-none -fopenmp
	    endif
	    CCFLAGS = -O3 -DLINUX
	endif
endif
ifeq ($(NAME),MAC)
	@echo "Sorry, no working Intel Fortran Compiler for Mac available"
endif

all: $(PROG)
	@echo "I think this is $(NAME)"

.PHONY: all
.PHONY: clean
create_build_dir:
	@mkdir -p ../build


#--------------------------------------------------------------------------
# Setting the rules to make the object files:
#--------------------------------------------------------------------------
main.o: pgtb.o
purification.o: purification_settings.o symmetry_i.o norms.o metrics.o
pgtb.o: metrics.o purification_settings.o symmetry_i.o purification.o
metrics.o: matrix_ops.o

# implizite Regel zur Erzeugung von *.o aus *.F ausschalten
%.o: %.F
fortran.o : $(CUDA)/src/fortran.c
	gcc $(CCFLAGS) -I$(CUDA)/src/ -c  $(CUDA)/src/fortran.c

# aus *.f mache ein *.o
%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# aus *.c mache ein *.o
%.o: %.c
	$(CC) $(CCFLAGS) -c $< -o $@

# linken
$(PROG): $(OBJS) 
	$(LINKER) $(OBJS) $(LIBS) -o $(PROG)

# dependencies
$(blas_level1.f90): $(accuracy.f90)
$(blas_level2.f90): $(accuracy.f90)
$(blas_level3.f90): $(accuracy.f90)
$(lapack_eig.f90): $(accuracy.f90)
$(la.f90): $(blas_level1.f90) $(blas_level2.f90) $(blas_level3.f90) $(lapack_eig.f90)

#aufraeumen
clean:
	rm *.o 
	rm *.mod 
	rm *genmod* 
	rm $(PROG)
