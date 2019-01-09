# Compiler options
HOSTNAME = $(firstword $(subst -, ,$(shell hostname)))
FC = gfortran
FCFLAGS = -O3 -ffast-math -funroll-loops -march=native 
DEBUG = -O0 -ffast-math -funroll-loops -march=native -fcheck=bounds -g -fbacktrace
DEBUGALL = -Wall -pedantic -fcheck=all -ffpe-trap=invalid,zero,overflow

# Required libraries: Lapack, FFTW3, HDF5
LIBS = -lm -L/usr/local/lib -L/usr/lib -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lfftw3 -llapack -lhdf5_fortran -lhdf5 -lhdf5_hl -lblas

# Source tree definitions
VPATH = src
BINDIR = bin
EXEC = scadyn

.SUFFIXES:
.SUFFIXES: .o .mod .f90 

# Hello world
yell = "Starting Make..."

# Includes and flag for putting .mod files to directory bin, home version
INCS = -I/usr/include -I/usr/local/include/ -I/usr/include/hdf5/serial/ -J${BINDIR}

# Dependency tree
OBJECTS = ${BINDIR}/common.o \
${BINDIR}/sfunctions.o \
${BINDIR}/h5io.o \
${BINDIR}/io.o \
${BINDIR}/integration_points.o \
${BINDIR}/translations.o \
${BINDIR}/mie.o \
${BINDIR}/possu.o \
${BINDIR}/sparse.o \
${BINDIR}/singularity_subtraction.o \
${BINDIR}/singularity_subtraction_N.o \
${BINDIR}/geometry.o \
${BINDIR}/sparse_mat.o \
${BINDIR}/precorrection.o \
${BINDIR}/projection.o \
${BINDIR}/build_G.o \
${BINDIR}/gmres_module.o \
${BINDIR}/setup.o \
${BINDIR}/T_matrix.o \
${BINDIR}/forces.o \
${BINDIR}/bessel.o \
${BINDIR}/shapebeam.o \
${BINDIR}/mueller.o \
${BINDIR}/integrator.o \
${BINDIR}/postprocessing.o \
${BINDIR}/main.o

###############################################################################
# Taito version for everything
FCTAITO = mpif90
INCSTAITO = -I/usr/include -I/usr/local/include/ -I${FFTW_ROOT}/include/ -I${H5ROOT}/include/ -m64 -I$(MKLROOT)/include/ -J${BINDIR}
LIBSTAITO = -L${FFTW_ROOT}/lib -lfftw3 -lfftw3_mpi -L${H5ROOT}/lib -lhdf5_fortran -lhdf5 -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_gf_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -ldl

ifeq ($(HOSTNAME),taito)
	yell = "Starting Make... When in Taito, remember to run 'module load gcc mkl fftw hdf5-serial' or Make will fail. Upon failure, load modules, clean and run Make again."
	FC := $(FCTAITO)
	INCS := $(INCSTAITO)
	LIBS := $(LIBSTAITO)
endif
###############################################################################

# No need to touch below, unless bad makefileing or messages need tweaking...
.PHONY: all clean
.SECONDARY: main-build

all: pre-build main-build post-build

pre-build:
	@echo $(yell)

main-build: ${EXEC} | $(BINDIR)

post-build:
	@echo "Target $(EXEC) compiled successfully"
	
debug: FCFLAGS = $(DEBUG)
debug: all

debugall: FCFLAGS = $(DEBUG) $(DEBUGALL)
debugall: all

$(BINDIR):
	@echo "Binary files are put into the directory $(BINDIR)"
	@mkdir -p ${BINDIR}

${BINDIR}/%.o: %.f90 |$(BINDIR)
	@echo "Compiling $^"
	@${FC} ${FCFLAGS} $(INCS) -c $^ -o $@ 

${EXEC}: ${OBJECTS}
	@echo "Linking the target"
	@${FC} ${FCFLAGS} ${OBJECTS} ${LIBS} -o $@

# Clean only objects
clean:
	@echo "Deleted all .o files"
	@rm -rf $(BINDIR)/*.o
	@rm -rf *.mod
	@rm -f *~

# Full clean
veryclean: clean
	@echo "Deleted all .mod files"
	@echo "Deleted executable and directory bin"
	@rm -f $(EXEC)
	@rm -rf $(BINDIR)

