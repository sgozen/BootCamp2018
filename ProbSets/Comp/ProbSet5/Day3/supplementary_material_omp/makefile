
#######################################################################	
all:normalize_vec.exec dot_prod.exec est_pi.exec
#######################################################################	


## Execute normalize_vec.f90
##normalize_vec.exec: normalize_vec.f90
##	gfortran normalize_vec.f90 -fopenmp -o normalize_vec.exec

## Execute dot_prod.f90

normalize_vec.exec: normalize_vec.f90
	gfortran normalize_vec.f90 -fopenmp -o normalize_vec.exec


dot_prod.exec: dot_prod.f90
	gfortran dot_prod.f90 -fopenmp -o dot_prod.exec

est_pi.exec: est_pi.f90
	gfortran est_pi.f90 -fopenmp -o est_pi.exec 
