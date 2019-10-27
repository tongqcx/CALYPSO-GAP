#FC = ifort -CA -CU -CB  -traceback -I /opt/intel/mkl/include
FC = ifort -I /opt/intel/mkl/include
LIB_MKL_ROOT=/opt/intel/mkl/lib/intel64
MKLLIBS= -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -liomp5 -lpthread
FLAGS =  -fopenmp  -c
BINPATH=./
EXEC = ngap.x
obj =  gap_constants.o math.o structure.o io.o linearalgebra.o gpr_base.o gpr_mb.o\
       blas.o exchng.o  read_module.o gap_module.o  setcellderv.o celldrv.o relax_module.o olinmin.o mlp_main.o main.o
${EXEC} : ${obj}
	${FC} -fopenmp -o ${BINPATH}${EXEC}  ${obj} ${LIB_MKL_PATH} ${MKLLIBS} 
.SUFFIXES: ${SUFFIXES} .f95 .f90 .f .F
%.o : %.f95
	$(FC) $(FLAGS) $< -o $*.o
%.o : %.f90
	$(FC) $(FLAGS) $< -o $*.o
%.o : %.F
	$(FC) -free $(FLAGS) $< -o $*.o
%.o : %.f
	$(FC) -free $(FLAGS) $< -o $*.o
clean:
	rm *.o  *.mod  
