real(8) function gen_gaussian(sigma)
! generate an Gaussian variable with variation=sigma, centered at zero
! test show that it take 4.57s to complete 1d8 loop
IMPLICIT NONE
	integer,save::flag=0
	real(8),save::v1,v2,s
	real(8) u1,u2,x,sigma
!	data flag/0/,s/0./
	real(8),external::rnd
!$OMP threadprivate(flag, v1,v2,s)
	if(flag==0)then
		s=0.
		do while(s>=1.or.s==0.)
            call random_number(u1)
            call random_number(u2)
			!u1=rnd(0d0,1d0)
			!u2=rnd(0d0,1d0)
			v1=2*u1-1
			v2=2*u2-1
			s=v1**2+v2**2
		end do
			x=v1*sqrt(-2*log(s)/s)
	!		print*, "x1=",x
        flag=1
	else
		x=v2*sqrt(-2*log(s)/s)
        flag=0
	!	print*, "x2=",x
	end if
	!flag=1-flag
	x=x*sigma
	gen_gaussian=x
	!print*, "x3=",x, sigma
end function

subroutine gen_gaussian_correlate(y1,y2, coeff)
	!Generate two random variables y1, and y2, with mean value <y1>=<y2>=0, 
	!dispersion <y1^2>=<y2^2>=1 and cross-correlation <y1y2>=coeff,where coeff<1
	implicit none
	real(8) y1p,y2p,y1,y2,coeff, gen_gaussian
	y1p=gen_gaussian(1d0); y2p=gen_gaussian(1d0)
	y1=y1p
	if(abs(coeff).le.1)then
		y2=y1p*coeff+y2p*sqrt(1-coeff**2)	
	else
		y2=coeff/abs(coeff)*y1
	end if
end subroutine
