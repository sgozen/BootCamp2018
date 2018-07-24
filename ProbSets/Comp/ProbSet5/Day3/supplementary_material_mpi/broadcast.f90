program broadcast
!----------------------------------------------------------------------   
!     program to demonstrate MPI_Bcast
!----------------------------------------------------------------------   

      use mpi 
 
      implicit none
      integer, parameter :: root = 0 !Define the root!
      integer :: rank, data, ierror

!----------------------------------------------------------------------         
      call MPI_Init(ierror)
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)
      
      data = 1
      if (rank == 0) then
        write(*,*) 'enter a value:'
        write(*,*) data
      end if

!.....broadcast the value of data of rank 0 to all ranks
      call MPI_BCAST(data, 1, MPI_INT, root, MPI_COMM_WORLD,ierror)
      write(*,*) "I am rank", rank, "and the value is", data

      call MPI_Finalize(ierror)
!----------------------------------------------------------------------   
      end program broadcast
!----------------------------------------------------------------------   
