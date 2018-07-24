program est_pi


  use mpi
  integer ::  N,i
  integer :: ierror, my_rank, input,size,num,sum
  double precision :: result, sum1,pi
  real*8 ::  est,Nc,Ns,mu,TIME
  real*8, dimension(:), allocatable :: x, y

  Nc = 0.0
  Ns = 0.0
  N = 80000000
  sum1 = 0
  sum = 0
  est = 0.0
  pi = 0.0  

  call srand(1) !Fix the seed
  allocate(x(N))
  allocate(y(N))
  do i=1,N
   x(i) = rand()
   y(i) = rand()
  enddo
  
  TIME = -MPI_Wtime()
  call MPI_INIT(ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)
  
  ! size is number of tasks
   num = N/size
      

  do i = 1,num
        IF ((x(i+num*my_rank)*x(i+num*my_rank) + y(i+num*my_rank)*y(i+num*my_rank)) .LE. 1.0) THEN 

        sum1 = sum1 + 1.0
        END IF
  enddo

 call MPI_REDUCE(sum1,result,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierror) 


if (my_rank == 0) then
 pi = (result/(N*1.d0))*4.0
 
 print*,'Estimated parallel Pi : ', pi
end if
 
 call MPI_FINALIZE(ierror) 


TIME = TIME+MPI_Wtime()
if (my_rank == 0) then

print*, 'Time it took in parallel code was : ', TIME
end if
end program est_pi

