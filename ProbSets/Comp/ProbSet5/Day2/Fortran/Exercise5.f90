program Exercise5
  integer,parameter :: seed = 1
  integer ::  N,i
  real*8 :: sum1, est,x,y,Nc,Ns,mu 
  Nc = 0.0
  Ns = 0.0
  N = 10000000
  sum1 = 0.0
  est = 0.0


  call srand(seed)
  
  do i = 1,N
        x = rand()*2 - 1
        y = rand()*2 - 1
        IF ((x*x + y*y) .LE. 1.0) THEN 
            Nc = Nc +1.0 
        ELSE
            Ns = Ns +1.0
        END IF
        est =  4.0*Nc/(Nc+Ns) 
        sum1 = sum1 + est
  enddo
mu = sum1/(Nc+Ns)
print*, mu
end program Exercise5

