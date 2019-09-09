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
BUILDDIR = build
EXEC = scadyn

.SUFFIXES:
.SUFFIXES: .o .mod .f90 

# Hello world
yell = "Starting Make..."

# Includes and flag for putting .mod files to directory bin, home version
INCS = -I/usr/include -I/usr/local/include/ -I/usr/include/hdf5/serial/ -J${BUILDDIR}

# Dependency tree
OBJECTS = ${BUILDDIR}/common.o \
${BUILDDIR}/sfunctions.o \
${BUILDDIR}/h5io.o \
${BUILDDIR}/io.o \
${BUILDDIR}/integration_points.o \
${BUILDDIR}/translations.o \
${BUILDDIR}/mie.o \
${BUILDDIR}/possu.o \
${BUILDDIR}/sparse.o \
${BUILDDIR}/singularity_subtraction.o \
${BUILDDIR}/singularity_subtraction_N.o \
${BUILDDIR}/geometry.o \
${BUILDDIR}/sparse_mat.o \
${BUILDDIR}/precorrection.o \
${BUILDDIR}/projection.o \
${BUILDDIR}/build_G.o \
${BUILDDIR}/gmres_module.o \
${BUILDDIR}/setup.o \
${BUILDDIR}/T_matrix.o \
${BUILDDIR}/forces.o \
${BUILDDIR}/bessel.o \
${BUILDDIR}/shapebeam.o \
${BUILDDIR}/mueller.o \
${BUILDDIR}/integrator.o \
${BUILDDIR}/postprocessing.o \
${BUILDDIR}/main.o

###############################################################################
# Taito version for everything
FCTAITO = mpif90
INCSTAITO = -I/usr/include -I/usr/local/include/ -I${FFTW_ROOT}/include/ -I${H5ROOT}/include/ -m64 -I$(MKLROOT)/include/ -J${BUILDDIR}
LIBSTAITO = -L${FFTW_ROOT}/lib -lfftw3 -lfftw3_mpi -L${H5ROOT}/lib -lhdf5_fortran -lhdf5 -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_gf_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -ldl

ifeq ($(HOSTNAME),taito)
	yell = "Starting Make... When in Taito, remember to run 'module load gcc mkl fftw hdf5-serial' or Make will fail. Upon failure, load modules, clean and run Make again."
	FC := $(FCTAITO)
	INCS := $(INCSTAITO)
	LIBS := $(LIBSTAITO)
endif
###############################################################################
# Puhti (test)version for everything
FCPUHTI = mpif90
INCSPUHTI = -I${FFTW_INSTALL_ROOT}/include/ -m64 -I$(MKLROOT)/include/ -I/${HOME}/include -J${BUILDDIR}
LIBSPUHTI = -L${FFTW_INSTALL_ROOT}/lib -lfftw3 -lfftw3_mpi -lhdf5_fortran -lhdf5 -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_gf_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -ldl

ifeq ($(HOSTNAME),puhti)
	yell = "Starting Make... When in Puhti, remember to run 'module load python-env gcc fftw/3.3.8-mpi hdf5', and remember to keep fingers crossed that all libraries are available (HDF5 is the problem child, who must be compiled locally at least for now, then loaded with 'module use <hdf5 install location>/include, and modify this Makefile accordingly'), or Make will fail. Upon failure, load modules, clean and run Make again."
	FC := $(FCPUHTI)
	INCS := $(INCSPUHTI)
	LIBS := $(LIBSPUHTI)
endif
###############################################################################

# No need to touch below, unless bad makefileing or messages need tweaking...
.PHONY: all clean
.SECONDARY: main-build

all: pre-build main-build post-build

pre-build:
	@echo $(yell)

main-build: ${EXEC} | $(BUILDDIR)

post-build:
	@echo "Target $(EXEC) compiled successfully"
	
debug: FCFLAGS = $(DEBUG)
debug: all

debugall: FCFLAGS = $(DEBUG) $(DEBUGALL)
debugall: all

$(BUILDDIR):
	@echo "Compiler files are put into the directory $(BUILDDIR)"
	@mkdir -p ${BUILDDIR}

${BUILDDIR}/%.o: %.f90 |$(BUILDDIR)
	@echo "Compiling $^"
	@${FC} ${FCFLAGS} $(INCS) -c $^ -o $@ 

${EXEC}: ${OBJECTS}
	@echo "Linking the target"
	@${FC} ${FCFLAGS} ${OBJECTS} ${LIBS} -o $@

# Clean only objects
clean:
	@echo "Deleted all .o files"
	@rm -rf $(BUILDDIR)/*.o
	@rm -rf *.mod
	@rm -f *~

# Full clean
veryclean: clean
	@echo "Deleted all .o and .mod files"
	@echo "Deleted executable and directory bin"
	@rm -f $(EXEC)
	@rm -rf $(BUILDDIR)

