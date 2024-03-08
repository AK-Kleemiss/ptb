PROG = ../build/ptb

SOURCES = accuracy.f90 la.f90 blas_level1.f90 blas_level2.f90 blas_level3.f90 lapack_eig.f90
$(foreach src, $(SOURCES), $(eval $(src) := $(patsubst %.f90, %.o, $(src))))

MODULES= bascom.o parcom.o com.o thresholds.o cbascom.o mocom.o dftd4.o aescom.o \
			$(foreach src, $(SOURCES), $($(src)))

CFILES= symmetry_i.o

OBJS2= rdcoord.o readl.o asym.o elem.o readl2.o readline.o prmat.o timing.o echo.o solve.o rdstvint.o \
       ints.o intpack.o pop.o rdmat.o fermi.o reorder.o bas.o occ.o dtrf2.o lin.o energy.o mrec.o \
       cma.o dtrf3.o drsp.o mlpop.o calcrab.o pgtb.o eeq.o ncoord.o prdipole.o itabrow.o \
       sscal.o cbas.o gauss.o printmos.o qsort.o rdtm.o shellocc_ref.o wr_control_atoms.o \
       egtbcorr.o symmetry.o guess_qsh.o eptb.o symmetry2.o symtranslib.o spinden.o camm.o json_output.o
       
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
  endif
  ifeq ($(UNAME_S), Darwin)
    NAME := MAC
  endif
endif

FC = ifx
CC  = icx
LINKER = $(FC) -static -I"$(MKLROOT)/include/intel64/lp64" -I"$(MKLROOT)/include"

all: $(PROG)
	@echo "I think this is $(NAME)"
ifeq ($(NAME),WINDOWS)
	@echo "Start making Windows executable"
    LINKER += -Qiopenmp
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
	@echo "Start making Linux executable"
    LINKER += -O -fopenmp
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
endif
ifeq ($(NAME),MAC)
	@echo "Sorry, no working Intel Fortran Compiler for Mac available"
endif

.PHONY: all
.PHONY: clean

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
	${RM} *.o
	${RM} *.mod
