!---------------------------------------------------------------------
!      program to demonstrate MPI_Scatter                                  
!---------------------------------------------------------------------  

      program scatter_number
   
      use mpi 
   
      implicit none
 
      integer :: i, rank, size, senddata(10), receivedata(2), ierror,root = 0
!----------------------------------------------------------------------      
      call MPI_Init(ierror)
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)
      call MPI_Comm_size(MPI_COMM_WORLD, size, ierror)
!----------------------------------------------------------------------   

      if (size .gt. 10) then
        if (rank == 0) then
          write(*,*) "do not use more than 10 processors"
          call MPI_Finalize(ierror)
        end if
      end if
      
!----------------------------------------------------------------------         
      if (rank  == 0) then
        do i = 1, 10, 1
       !   write(*,*) 'enter The first value for processor', i
       !   read(*,*) senddata(2*i-1)

          senddata(i) = i
       !   write(*,*) 'enter The second value for processor', i
       !   read(*,*) senddata(2*i)

        end do
      end if
!.....scatter the value of senddata of rank 0 to receivedata of all ranks
      call MPI_SCATTER(senddata, 2, MPI_INT, receivedata, &
                       2, MPI_INT, root, MPI_COMM_WORLD, ierror)
     
      if (rank<= size) then
       write (*,*) "I am rank", rank, "and the value is", receivedata
      end if
      call MPI_Finalize(ierror)
!----------------------------------------------------------------------   
      end program scatter_number
!---------------------------------------------------------------------- 
