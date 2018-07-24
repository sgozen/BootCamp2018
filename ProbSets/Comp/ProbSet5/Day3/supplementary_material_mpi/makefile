#######################################################################	
all:broadcast.exec allreduce.exec scatter.exec est_pi.exec
#######################################################################	


## Execute normalize_vec.f90
##normalize_vec.exec: normalize_vec.f90
##	gfortran normalize_vec.f90 -fopenmp -o normalize_vec.exec

## Execute dot_prod.f90

broadcast.exec: broadcast.f90
	mpif90 -o3 broadcast.f90 -o broadcast.exec


allreduce.exec: allreduce.f90
	mpif90 -o3 allreduce.f90 -o allreduce.exec


scatter.exec: scatter.f90
	mpif90 -o3 scatter.f90 -o scatter.exec


est_pi.exec: est_pi.f90
	mpif90 -o3 est_pi.f90 -o est_pi.exec
