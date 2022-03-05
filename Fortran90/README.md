The folder src/ contains a Fortran90 version of the Very Fast Simulated Annealing global optimization algorithm. The files op_vfsa90 and param_vfsa90 are the input files of the algorithm

To compile:

gfortran -O3 -o vfsa \
	 precision.f90 \
	 variables.f90 \
	 count_file_lines_mod.f90 \
	 fun_cost_mod.f90 \
	 read_parameters.f90 \
	 sub_vfsa90.f90 \
	 vfsa90.f90 
