real(8) function rnd(min,max)
		IMPLICIT NONE
		real(8) min,max,temp
		CALL RANDOM_NUMBER (temp)
		if(min<max)then
			rnd=temp*(max-min)+min
		end if
		if(min>max) then
			rnd=temp*(min-max)+max
		end if
end function
SUBROUTINE set_system_random_seed()
	implicit none
	INTEGER :: i, n, clock
	INTEGER, DIMENSION(:), ALLOCATABLE :: seed

	CALL RANDOM_SEED(size = n)
	ALLOCATE(seed(n))

	CALL SYSTEM_CLOCK(COUNT=clock)

	seed = clock + 37 * (/ (i - 1, i = 1, n) /)
	CALL RANDOM_SEED(PUT = seed)

	DEALLOCATE(seed)
END SUBROUTINE

subroutine same_random_seed(seed_value)
    implicit none
    integer seed_int, seed_value
    integer,allocatable::seed(:)
 
   	call random_seed(size=seed_int)
	allocate(seed(seed_int))
	seed(1:seed_int)=seed_value
    call random_seed(put=seed)
end subroutine
function rndI(min,max)
		IMPLICIT NONE
		Integer min,max, rndI
		real(8) temp
		CALL RANDOM_NUMBER (temp)
		if(min<max)then
			rndI=int(temp*(max-min+1)+min)
		end if
		if(min>max) then
			rndI=int(temp*(min-max+1)+max)
		end if
		if(min.eq.max) rndI=min
end function

real(8) function fPowerLaw(yta,xmin,xmax)
	implicit none
	real(8) yta,xmin,xmax
	real(8),external::rnd
	if (yta==-1d0) then
		fPowerLaw=exp(rnd(log(xmin),log(xmax)))
	else
		fPowerLaw=(rnd(xmin**(1+yta),xmax**(1+yta)))**(1d0/(1d0+yta))
	end if
end function