# --- Variables ---
PROG = ../build/ptb_dev
ACCEL = /pc2/users/m/mpatbf1/libm/
OSTYPE = LINUXI
WITH_DFTD4 = FALSE

# --- Code-specific variables ---
SOURCES = accuracy.f90 cuda_context.f90 timer.f90 la.f90 blas_level1.f90 blas_level2.f90 blas_level3.f90 lapack_eig.f90
$(foreach src, $(SOURCES), $(eval $(src) := $(patsubst %.f90, %.o, $(src))))

MODULES= bascom.o parcom.o com.o thresholds.o cbascom.o mocom.o dftd4.o aescom.o \
			$(foreach src, $(SOURCES), $($(src)))

CFILES= symmetry_i.o

OBJS1 = matrix_ops.o metrics.o purification_settings.o norms.o  purification.o \
		 rdcoord.o readl.o asym.o elem.o readl2.o readline.o prmat.o  timing.o echo.o solve.o rdstvint.o \
       ints.o intpack.o pop.o rdmat.o fermi.o reorder.o bas.o occ.o dtrf2.o lin.o energy.o mrec.o \
       cma.o dtrf3.o drsp.o mlpop.o calcrab.o symmetry.o pgtb.o eeq.o ncoord.o prdipole.o itabrow.o \
       sscal.o cbas.o gauss.o printmos.o qsort.o rdtm.o shellocc_ref.o wr_control_atoms.o \
       egtbcorr.o guess_qsh.o eptb.o symtranslib.o spinden.o camm.o json_output.o

MAINFILES= main.o 

OBJS = $(MODULES) $(CFILES) $(OBJS2) $(MAINFILES)
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
	"$(MKLROOT)\lib\mkl_intel_thread.lib" \
	"$(MKLROOT)\..\..\2024.0\lib\libiomp5md.lib"
	FFLAGS = -O3 -Qiopenmp
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
	    FFLAGS = -O3 -align array64byte -fopenmp
	    CCFLAGS = -O -DLINUX
	else
	    FC = gfortran
	    CC = gcc
	    MKLROOT = /usr/lib/x86_64-linux-gnu
	    LINKER = $(FC) -static -I"/usr/include/mkl/intel64/lp64" -I"/usr/include" -O -g
	    LIBS = \
		-L"/usr/lib/x86_64-linux-gnu/lapack" -llapack -L"/usr/lib/x86_64-linux-gnu/blas" -lblas -lpthread -lm
	    #$(MKLROOT)/libmkl_blas95_ilp64.a \
	    #$(MKLROOT)/libmkl_lapack95_ilp64.a \
	    #-Wl,--start-group \
	    #$(MKLROOT)/libmkl_gf_ilp64.a \
	    #$(MKLROOT)/libmkl_core.a \
	    #$(MKLROOT)/libmkl_sequential.a \
	    #-Wl,--end-group -lpthread -lm
	    FFLAGS = -g -fbacktrace -fcheck=all
		#FFLAGS = -O3 -fopenmp
	    CCFLAGS = -O -DLINUX
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
# implizite Regel zur Erzeugung von *.o aus *.F ausschalten
%.o: %.F
fortran.o : $(CUDA)/src/fortran.c
	gcc $(CCFLAGS) -I$(CUDA)/src/ $(INC) -c  $(CUDA)/src/fortran.c

# aus *.f mache ein *.o
%.o: %.f
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

%.o: %.f90
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) -c $< -o $@

# aus *.c mache ein *.o
%.o: %.c
	@echo "making $@ from $<"
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
	rm -rf *.o *.mod *genmod* $(PROG) ../build
