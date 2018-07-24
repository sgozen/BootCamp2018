!----------------------------------------------------------------------                                                               
!  a grogram to demonstrate a reduce to all                              
!----------------------------------------------------------------------   

      program allreduce
      
      use mpi

      implicit none

      integer :: ierror, my_rank, input, result
      integer :: sum

!----------------------------------------------------------------------         
      call MPI_Init(ierror)
      call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierror)

!.....TBD calculate sum of all ranks
      input = my_rank + 1

      call MPI_Reduce(input,result,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD,ierror)
     
      if (my_rank .eq. 0) then
        write(*,*) 'Rank', my_rank, ": Sum =", result
      end if
      call MPI_Finalize(ierror)

!----------------------------------------------------------------------         
      end program allreduce
!----------------------------------------------------------------------
