program Exercise1

   implicit none
   real*8, f,p,howspeed

   f = 0.4/100.0
   p = 100.0
   howspeed = 1.0/(f+(1.0-f)/p)
  
   
   write(*,*), 'Number of Cores? : ', p
   write(*,*), 'Fraction of Serial Code? : ', f
   write(*,*), 'Maximum Sppedup? : ', howspeed

end program Exercise1

