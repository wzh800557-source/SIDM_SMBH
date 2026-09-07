module fun_math
implicit none
contains	
	real(8) function farravg(arr,n)
		integer n
		real(8) arr(n)
		farravg=sum(arr)/dble(n)
	end function
	real(8) function farrsct(arr, n)
		integer n
		real(8) arr(n),mean
		mean=farravg(arr,n)
		farrsct=(sum((arr-mean)**2)/dble(n-1))**0.5d0
	end function
	real(8) function farrsct2(arr,mean, n)
		integer n
		real(8) arr(n), mean
		farrsct2=(sum((arr-mean)**2)/dble(n-1))**0.5d0
	end function
end module
