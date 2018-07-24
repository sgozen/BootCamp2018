program est_pi

  use omp_lib
  integer ::  N=100000000,i, seed=1,id,omp_get_num_threads
  real*8 :: TIME, sum, pi_para,pi_serial
  real*8, dimension(:), allocatable :: x, y
  allocate(x(N))
  allocate(y(N))
  
  sum = 0.0


call srand(seed)
do i =1,N
x(i) = rand()
y(i) = rand()
end do
TIME = -omp_get_wtime()

!$OMP PARALLEL DO REDUCTION(+:sum) private(i)
  do i = 1,N
        IF ((x(i)*x(i) + y(i)*y(i)) .LE. 1.0) THEN 
        sum = sum + 1.0
        END IF 
  enddo
!$OMP END PARALLEL DO
pi_para = (sum/N)*4.0
TIME = TIME + omp_get_wtime()
print*,'Estimated parallel Pi : ', pi_para

print*,'Time it took in parallel code : ', TIME


TIME = -omp_get_wtime()

  sum1 = 0.0
  sum = 0.0
  do i = 1,N
        IF ((x(i)*x(i) + y(i)*y(i)) .LE. 1.0) THEN 
        sum = sum + 1.0
        END IF 
  enddo


pi_serial = (sum/N)*4.0
TIME = TIME + omp_get_wtime()

print*,'Estimated serial Pi : ', pi_serial

print*,'Time it took in serial code : ', TIME
end program est_pi

