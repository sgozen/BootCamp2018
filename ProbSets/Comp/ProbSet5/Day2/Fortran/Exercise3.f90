program Exercise3
                  
    implicit none
    real*8 :: a,b,c,x1,x2
   
    write(*,*) 'The value of a is : '    
    read(*,*) a
    write(*,*) 'The value of b is : '
    read(*,*) b
    write(*,*) 'The value of c is : '
    read(*,*) c    

    call Quadratic_sol(a,b,c,x1,x2)   

   
  print *, x1   
  print *, x2

end program Exercise3


subroutine Quadratic_sol(a,b,c,x1,x2)

    implicit none
    real*8,intent(in) :: a,b,c
    real*8,intent(inout) :: x1,x2
  
    x1 = (-b+((b**2.0 - 4.0*a*c))**(1.0/2.0))/(2.0*a)    
    x2 = (-b-((b**2.0 - 4.0*a*c))**(1.0/2.0))/(2.0*a)

end subroutine Quadratic_sol



