# ------------------Macro-Defs---------------------
CULASPARSE_INC_PATH=#/opt/culasparse/include/
CULASPARSE_LIB_PATH_64=#/opt/culasparse/lib64
CUDA_LIB_PATH_64=/usr/local/cuda-8.0/lib64
CUDA_INCLUDE = /usr/local/cuda-8.0/include/ 
CUDA_COMMON_INCL = /usr/local/cuda-8.0/samples/common/inc/
CUSP_INCL = /opt/cusp/0.5.1/
F90=ifort
CU=nvcc
CU_COMPILER_OPTS = -Xptxas -O3
COMPILER_OPTS=-O3 -traceback -xHost
COMP_OPTS_EXT=-I${CUDA_INCLUDE} -I${CUDA_COMMON_INCL} -I${CUSP_INCL}
LIBSPATH=-L${CUDA_LIB_PATH_64}
LIBS= -lcudart -lstdc++ #-lcula_sparse -lcolamd -liomp5  #-lcublas -lcusparse
MODULES= #${CULASPARSE_INC_PATH}/cula_sparse_legacy.f90
# -------------------End-macro-Defs---------------------------

 OBJ =\
modcusp_library.o\
modcusp_library_intrf.o\
precision.o\
real_parameters.o\
parameters.o\
shared_data.o\
modfd_openclose_ops.o\
modfd_set_bc.o\
modfd_tecwrite.o\
mod_create_filenum.o\
modfd_create_geom.o\
modfd_solve_linearsys.o\
modfd_problem_setup.o\
modfd_calc_pre.o\
modfd_calc_temp.o\
modfd_calc_mom.o\
modfd_calc_integrals.o\
solver_sparsekit2.o\
FD2DC.o

FD2DC:$(OBJ)
	 $(F90) -o FD2DC $(OBJ) $(COMPILER_OPTS) $(LIBSPATH) $(LIBS)   

modcusp_library_intrf.o:./modcusp_library_intrf.f90 modcusp_library.o
	 $(F90) -c $(COMPILER_OPTS) ./modcusp_library_intrf.f90

modcusp_library.o:./modcusp_library.cu modcusp_library.h
	$(CU) -c $(COMP_OPTS_EXT) $(CU_COMPILER_OPTS) ./modcusp_library.cu

solver_sparsekit2.o:./solver_sparsekit2.f
	 $(F90) -c $(COMPILER_OPTS) ./solver_sparsekit2.f

precision.o:./precision.f90  
	 $(F90) -c $(COMPILER_OPTS) ./precision.f90 

real_parameters.o:./real_parameters.f90 precision.o 
	 $(F90) -c $(COMPILER_OPTS) ./real_parameters.f90 

parameters.o:./parameters.f90  
	 $(F90) -c $(COMPILER_OPTS) ./parameters.f90 

shared_data.o:./shared_data.f90 parameters.o precision.o 
	 $(F90) -c $(COMPILER_OPTS) $(COMP_OPTS_EXT) ./shared_data.f90 

modfd_openclose_ops.o:./modfd_openclose_ops.f90 parameters.o shared_data.o parameters.o shared_data.o 
	 $(F90) -c $(COMPILER_OPTS) ./modfd_openclose_ops.f90 

modfd_set_bc.o:./modfd_set_bc.f90 precision.o shared_data.o real_parameters.o  
	 $(F90) -c $(COMPILER_OPTS) ./modfd_set_bc.f90 

modfd_tecwrite.o:./modfd_tecwrite.f90 precision.o shared_data.o parameters.o precision.o  
	 $(F90) -c $(COMPILER_OPTS) ./modfd_tecwrite.f90 

mod_create_filenum.o:./mod_create_filenum.f90 precision.o 
	 $(F90) -c $(COMPILER_OPTS) ./mod_create_filenum.f90 

modfd_create_geom.o:./modfd_create_geom.f90 shared_data.o precision.o parameters.o real_parameters.o modfd_tecwrite.o   
	 $(F90) -c $(COMPILER_OPTS) ./modfd_create_geom.f90 

modfd_solve_linearsys.o:./modfd_solve_linearsys.f90 precision.o parameters.o shared_data.o real_parameters.o  
	 $(F90) -c $(COMPILER_OPTS) ./modfd_solve_linearsys.f90 

#modcu_BiCGSTAB.o:./modcu_BiCGSTAB.f90  
#	 $(F90) -c $(COMPILER_OPTS) ./modcu_BiCGSTAB.f90 

#modcu_BiCGSTAB.o
modfd_problem_setup.o:./modfd_problem_setup.f90 real_parameters.o parameters.o precision.o shared_data.o modfd_set_bc.o modfd_create_geom.o modfd_tecwrite.o modfd_solve_linearsys.o parameters.o modcusp_library_intrf.o
	 $(F90) -c $(COMPILER_OPTS) ./modfd_problem_setup.f90 

#modcu_BiCGSTAB.o
modfd_calc_pre.o:./modfd_calc_pre.f90 parameters.o real_parameters.o precision.o modfd_set_bc.o shared_data.o modfd_solve_linearsys.o solver_sparsekit2.o modcusp_library_intrf.o
	 $(F90) -c $(COMPILER_OPTS) ./modfd_calc_pre.f90 

#modcu_BiCGSTAB.o
modfd_calc_temp.o:./modfd_calc_temp.f90 parameters.o real_parameters.o precision.o modfd_set_bc.o shared_data.o modfd_solve_linearsys.o solver_sparsekit2.o modcusp_library_intrf.o
	 $(F90) -c $(COMPILER_OPTS) ./modfd_calc_temp.f90 

#modcu_BiCGSTAB.o 
modfd_calc_mom.o:./modfd_calc_mom.f90 parameters.o real_parameters.o shared_data.o modfd_set_bc.o precision.o modfd_solve_linearsys.o solver_sparsekit2.o modcusp_library_intrf.o
	 $(F90) -c $(COMPILER_OPTS) ./modfd_calc_mom.f90 

modfd_calc_integrals.o:./modfd_calc_integrals.f90 precision.o real_parameters.o shared_data.o parameters.o modfd_solve_linearsys.o  
	 $(F90) -c $(COMPILER_OPTS) ./modfd_calc_integrals.f90 

#modcu_BiCGSTAB.o
FD2DC.o:./FD2DC.f90 real_parameters.o precision.o modfd_openclose_ops.o modfd_problem_setup.o shared_data.o parameters.o modfd_set_bc.o modfd_calc_pre.o modfd_calc_temp.o modfd_calc_mom.o mod_create_filenum.o modfd_calc_integrals.o modfd_tecwrite.o modfd_create_geom.o solver_sparsekit2.o modcusp_library_intrf.o
	 $(F90) -c $(COMPILER_OPTS) FD2DC.f90
# This entry allows you to type " make clean " to get rid of
# all object and module files 
clean:
	rm -f -r f_{files,modd}* *.o *.mod *.M *.d V*.inc *.vo \
	V*.f *.dbg album F.err
 
