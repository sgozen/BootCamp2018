program Exercise4
                  
    implicit none
    real*8 :: x, pi, sum, step
    integer :: N_step, i
   
   N_step = 10000
   sum = 0.0   
   step = 1.0/N_step

   do i = 1,N_step
        
        x = (i+0.5)*step
        sum = sum + 4.0/(1.0 + x*x) 
   
   enddo

    print *, sum*step
 
end program Exercise4

