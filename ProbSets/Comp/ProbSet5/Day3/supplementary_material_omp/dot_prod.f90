!----------------------------------------------------------------------   
      program dot
!----------------------------------------------------------------------   

      use omp_lib

      implicit none

      integer      :: N = 100000000
      integer      :: i, num_threads,id
      real(8), dimension(:), allocatable  :: a, b
      real(8) :: time, sum, expected

!----------------------------------------------------------------------         
      allocate(a(N), b(N))

      num_threads = omp_get_max_threads()
      write(*,*) 'dot of vectors with length ', N, ' with ', num_threads, ' threads'


!.....set timer      
      time = -omp_get_wtime();

!$OMP PARALLEL DO
!.....init. arrays      
      do i=1,N
        a(i) = 1.0d0/2.0d0
        b(i) = i
       id = omp_get_num_threads()
 
      enddo

!$OMP END PARALLEL DO 

      sum = 0.0d0


!$OMP PARALLEL DO REDUCTION(+:sum)
      do i=1,n
        sum = sum + a(i) * b(i)
      end do
      
!$OMP END PARALLEL DO 
!.....get time      
      time = time + omp_get_wtime()

      expected = (N+1.0d0)*N/4.0d0;
      
      write(*,*) 'relative error ', abs(expected-sum)/expected
      write(*,*) 'took ', time, ' seconds'

      deallocate(a, b)

!----------------------------------------------------------------------   
      end program dot
!----------------------------------------------------------------------   

