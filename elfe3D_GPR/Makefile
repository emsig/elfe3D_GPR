# Makefile for elfe3D_GPR
# ==========================================
# Configuration variables (override via environment or command line):
#   FCC=<compiler>              - Fortran compiler (default: gfortran)
#   MUMPS_LIB_DIR=<path>        - Path to MUMPS library (relative or absolute)

# Compiler - allow override: make FCC=gfortran or FCC=gfortran-10
FCC ?= gfortran-10

# Path configuration - use environment variable or default to relative path
# Example: MUMPS_LIB_DIR=../path/to/lib make
MUMPS_LIB_DIR ?= ../../Setup_Aid/elfe3D/Dependencies/MUMPS_5.7.3/lib

# System information (for reference)
OS = $(shell uname)
ARCH = $(shell uname -m)
HOSTNAME = $(shell hostname)

# Object files created during compilation
OBJ =   mod_constant.o \
	mod_types_basic.o \
	mod_util.o \
	mod_vector_products.o \
	mod_define_model.o \
	mod_read_mesh.o \
	mod_calculate_matrices.o \
	mod_model_parameters.o \
	mod_interp_functions.o \
	mod_calculate_local_left.o \
	mod_calculate_global_source.o \
	mod_BC.o \
	mod_sparse_matrix_operations.o \
	mod_solvers.o \
	mod_calculate_tf.o \
	mod_error_estimates.o \
	mod_tetgen_operations.o \
	elfe3d_gpr.o

# Compiler flags
COMMON_FFLAGS = -O2 -march=native -mtune=native -std=f2008 -ffree-form -ffree-line-length-none -fopenmp

# Compiler-specific optimizations
ifeq ($(FCC),gfortran)
    FFLAGS = $(COMMON_FFLAGS) -mfma -malign-data=cacheline -ftree-vectorize -ftree-vectorizer-verbose=2 -fopt-info-vec -fall-intrinsics
else ifeq ($(FCC),gfortran-10)
    FFLAGS = $(COMMON_FFLAGS) -mfma -malign-data=cacheline -ftree-vectorize -ftree-vectorizer-verbose=2 -fopt-info-vec -fall-intrinsics
else
    FFLAGS = $(COMMON_FFLAGS)
endif

COMPFLAGS = -x f95-cpp-input

# MUMPS and related libraries
OBJUTIL = -L$(MUMPS_LIB_DIR) -L/usr/local/lib -L/usr/lib -lpthread -lm -lzmumps -lmumps_common -lesmumps -lpord -lmpiseq -lscotch -lmetis -lopenblas

##-------------------------------------------------------------------------------------
#
#  Default command.
#
.DEFAULT:
	@echo "Unknown target $@, try:  make help"

#
#  Local help file.
#
help:
	@echo "   make all    -  make all object, module and executable files"
	@echo "   make clean  -  remove all object, module and executable files"
	@echo ""
	@echo "Configuration (override via environment or command line):"
	@echo "   FCC=<compiler>              - Fortran compiler (default: $(FCC))"
	@echo "   MUMPS_LIB_DIR=<path>        - Path to MUMPS library (default: $(MUMPS_LIB_DIR))"
	@echo ""
	@echo "Examples:"
	@echo "   make all"
	@echo "   make FCC=gfortran-10 all"
	@echo "   make MUMPS_LIB_DIR=/opt/mumps/lib all"

#
# Dependencies of object files and their compilation rules
#
$(OBJ): %.o: %.f90 
	$(FCC) $(FFLAGS) $(COMPFLAGS) -c $< -o $@ 

#
# Linking stage
#
all: $(OBJ)
	$(FCC) $(FFLAGS) $(OBJ) $(OBJUTIL) -o elfe3d_gpr


# clean up
.PHONY: clean

clean:
	-rm *.o *.mod in/*.vol
	@echo "Clean complete"
