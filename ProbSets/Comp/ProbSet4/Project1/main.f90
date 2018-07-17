program BS

    implicit none
    real*8 :: S,K,r,V,T,call1,put1, payoff_sum

    integer :: num_sims

num_sims = 1000000
S = 100.0
K = 100.0
r = 0.05
V = 0.2
T = 1.0
call monte_carlo_call_price(num_sims,S,K,r,V,T,payoff_sum,call1)

call monte_carlo_put_price(num_sims,S,K,r,V,T,payoff_sum,put1)
print*,'The Call Price is : ', call1
print*,'The Put Price is : ', put1
end program BS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine monte_carlo_call_price(num_sims,S,K,r,V,T,payoff_sum,call1)
 implicit none
 real*8, intent(in) :: S,K,r,V,T
 real*8, intent(inout) :: call1,payoff_sum
 
 integer, intent(in) :: num_sims
 integer :: i
 real*8 :: gauss_bm, S_cur, S_adjust,result

S_adjust = S*exp(T*(r-0.5*V*V))
payoff_sum = 0.0
do i =1,num_sims

 call gaussian_box_muller(result)
 gauss_bm = result
 S_cur = S_adjust*exp(sqrt(V*V*T)*gauss_bm)
 payoff_sum = payoff_sum + max(S_cur-K,0.0)

enddo
call1 = payoff_sum/num_sims  * exp(-r*T)

end subroutine monte_carlo_call_price
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine monte_carlo_put_price(num_sims,S,K,r,V,T,payoff_sum,put1)
 implicit none
 real*8, intent(in) :: S,K,r,V,T
 real*8, intent(inout) :: put1, payoff_sum

 integer, intent(in) :: num_sims
 integer :: i
 real*8 :: gauss_bm, S_cur, S_adjust,result

S_adjust = S*exp(T*(r-0.5*V*V))
payoff_sum = 0.0

do i =1,num_sims

 call gaussian_box_muller(result)
 gauss_bm = result
 S_cur = S_adjust*exp(sqrt(V*V*T)*gauss_bm)
 payoff_sum = payoff_sum + max(K-S_cur,0.0)

enddo
put1 = payoff_sum/num_sims  * exp(-r*T)

end subroutine monte_carlo_put_price

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gaussian_box_muller(result)
 implicit none
 real*8, intent(inout) :: result
 real*8 :: x,y, euclid_sq

euclid_sq = 100.0


	do while  (euclid_sq >= 1.0)

	x = 2.0*rand() - 1.0
	y = 2.0*rand() - 1.0
	euclid_sq = x*x + y*y

	end do

result = x*sqrt(-2*log(euclid_sq)/euclid_sq)

end subroutine gaussian_box_muller
