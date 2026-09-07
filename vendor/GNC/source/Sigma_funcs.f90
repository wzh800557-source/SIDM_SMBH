
subroutine get_sigma0(energyx,  funcs, sigma)
	use com_sts_type
	use constant
	use md_coeff
    use my_intgl
    use model_basic,only: fgx_g0
    use sidm_params, only: use_sidm, log_lambda_eff
	implicit none
	integer i,idid
	real(8) energyx, sigma, ecc,vel,rad,semi, int_out
    !real(8) intne, intpo
	external::funcs
!	integer, parameter::nintn=25
	int_out=0
	!print*, "inf=",inf
    idid=0
	!if(energyx.ne.0d0)then
        !intne=0d0; intpo=0d0
        !print*, emin_factor
    !print*, "energyx=",energyx
    if(energyx/emin_factor>1d0+1d-7)then
		call my_integral_none(emin_factor/energyx, 1d0, int_out, fcn,idid)
        if(idid.lt.0)then
            print*, "error! idid=",idid, energyx, emin_factor
            stop
        end if
    end if
        !if(emin_factor<energyx)then
        !    call my_integral(emin_factor/energyx, 1d0, intpo, fcn,idid)
        !end if
        !int_out=intne+intpo
	!else
!		call my_integral_none(inf, 1d0, int_out, fcn,idid)
!        stop
!	end if
	if (use_sidm) then
	    sigma = int_out + fgx_g0/energyx * log_lambda_eff(energyx)
	else
	    sigma = int_out + fgx_g0/energyx
	end if
	!print*, "sigma=", energyx, sigma
	!read(*,*)
	if(idid.eq.-2)print*, "sigma0=",sigma
contains
	subroutine fcn(n, x, y, f, par, ipar)
		implicit none
		integer n, ipar(100)
		real(8) x, y(n), f(n), par(100), Fx
		call funcs(x*energyx, fx)
		if (use_sidm) then
		    f(1) = fx * log_lambda_eff(x * energyx)
		else
		    f(1) = fx
		end if
        !print*, "x,energxy,E,fx=",x,energyx,x*energyx,fx
		!if(coeff_chattery.ge.3) 
        !print*, x, energyx,x*energyx, fx
		if(isnan(fx).or.isnan(x))then
			print*, emin_factor, x*energyx, fx
			stop
		end if

	end subroutine
end subroutine

subroutine get_sigma_funcs_cfs_grid(energyx,jum, funcs, cfs_grid, &
    nx, ny, jmin,jmax, sigma)
    use my_intgl
    use constant
    use sidm_params, only: use_sidm, log_lambda_eff
    implicit none
    integer idid, nx,ny, i,j
    real(8) energyx, jum, sigma, ecc
    real(8) yout, jmin, jmax, ss,smax,smin,dsmin, dsmax
    real(8) cfs_grid(nx,ny), cfs, sstep,fx,ds
    real(8) get_cfs
    real(8) lnL_local
    integer nstep
    external::funcs
    if(jum>1d0)then
        print*, "error, jum>1:", jum
        stop
    end if
    ecc=sqrt(1-jum**2)
    dsmin=-7d0; 
    smin=log10(10**dsmin+1d0)
    dsmax=log10(2/(1-ecc)-1d0)
    smax=log10(10**dsmax+1d0)

    nstep=nx*6
    sstep=(dsmax-dsmin)/real(nstep-1)
    yout=0
    do i=1, nstep
        ds=(dsmax-dsmin)*real(i-1)/real(nstep-1)+dsmin  
        ss=log10(10**ds+1d0)    
        cfs=get_cfs(jum,ds,cfs_grid,nx,ny, jmin, jmax, dsmin, dsmax)
        call funcs(10**ss*energyx, fx)

        if (use_sidm) then
            lnL_local = log_lambda_eff(10**ss * energyx)
        else
            lnL_local = 1.0d0
        end if
        yout=yout+fx*cfs*lnL_local*(10**ss-1)*sstep
    
    end do
    
    sigma=1d0/pi*yout*log(10d0)
    !if(sigma<0)then
    !    print*, "s, ecc, sstep, sigma=",10**ss, ecc, sstep, sigma
    !    stop
    !end if
end subroutine

real(8) function get_cfs(jum,y, cfs, nx, ny, jmin ,jmax,dsmin,dsmax)
    implicit none
    real(8) x, y
    integer nx, ny
    real(8) cfs(nx, ny)
    real(8) xmin, xmax, ymin, ymax, jmin, jmax,jum,ecc,dsmin,dsmax
    integer idx, idy
    integer i,j
    xmin=jmin; xmax=jmax
    x=log10(jum)
    ymin=dsmin
    ymax=dsmax

    if(x<xmin)then
        x=xmin
    endif
    if(x>xmax)then
        print*, "warnning, jmax, j=",jum, 10**x
        x=xmax
    end if
    if(y<ymin)then
        if(abs(y-ymin)>1d-6)then
            print*, "warnning, smin,s=",ymin, y
        end if
        y=ymin
    end if
    if(y>ymax)then
        if(abs(y-ymax)>1d-6)then
            print*, "warnning, smax,s=",ymax, y
        end if
        y=ymax
    end if

    call return_idxy(x,y,xmin,xmax,ymin,ymax,nx,ny,idx,idy,1)
    if((idx>nx.or.idx<-1).or.(idy>ny.or.idy<-1))then
        print*, "idx,idy=",idx,idy
        print*, "x, y, xmin, xmax, ymin, ymax=",  x, y, xmin, xmax, ymin, ymax
        stop
    end if
    !idx=min(max(idx,1), nx)
    !idy=min(max(idy,1), ny)
    !print*, "x,y=",x,y
    !print*, "idx,idy=",idx,idy
    !call linear_int(x,y,n, xev,yev)
    get_cfs=cfs(idx,idy)
    !print*, "get_cfs=",get_cfs
    !read(*,*)
end function


real(8) function get_cfs_blinear(jum,y, cfs, nx, ny, jmin ,jmax,dsmin,dsmax)
    implicit none
    real(8) x, y
    integer nx, ny
    real(8) cfs(nx, ny)
    real(8) xmin, xmax, ymin, ymax, jmin, jmax,jum,ecc,dsmin,dsmax
    real(8) xstep,ystep
    integer idx, idy
    integer i,j
    
    if(x>jmax)then
        print*, "warnning, jmax, j=",jum, 10**x
       ! x=xmax
    end if
    if(y>dsmax)then
        if(abs(y-dsmax)>1d-6)then
            print*, "warnning, dsmax,s=",dsmax, y
        end if
        y=ymax
    end if

    ystep=(dsmax-dsmin)/real(ny-1)
    xstep=(jmax-jmin)/real(nx-1)
    call linear_int_2d(jmin,dsmin,nx,ny,xstep, ystep,cfs,log10(jum),y,get_cfs_blinear)
end function

subroutine get_sigma_funcs_cfs_rk(energyx,jum, funcs, cfs_grid, &
    nx, ny, jmin,jmax, sigma)
    use my_intgl
    use constant
    use md_cfuns,only:dsmin_value
    implicit none
    integer idid, nx,ny
    real(8) energyx, jum, sigma, ecc
    real(8) yout, jmin, jmax,dsmin,dsmax
    real(8) cfs_grid(nx,ny)
    external::funcs

    !print*, jmin, jmax, jum
    dsmin=dsmin_value; 
    ecc=sqrt(1-jum**2)
    dsmax=log10(2/(1-ecc)-1d0)
    yout=0
    call my_integral_acc(log10(1d0+10**dsmin),log10(2d0/(1-ecc)), yout, 1d-13,1d-12, fcn, idid)
    sigma=1d0/pi*yout*log(10d0)
    !print*, "sigma=",sigma
contains
    subroutine fcn(n, x, y, f, par, ipar)
    use, intrinsic :: ieee_arithmetic
    implicit none
    integer n, ipar(100)
    real(8) x, y(n), f(n), par(100)
    real(8),external::get_cfs_blinear
    real(8) cfs, fx
    !cfs=get_cfs(jum,log10(10**x-1),cfs_grid, nx, ny, jmin, jmax,dsmin,dsmax)
    cfs=get_cfs_blinear(jum,log10(10**x-1),cfs_grid, nx, ny, jmin, jmax,dsmin,dsmax)
    !print*, "set dsmin,dsmax"
   ! stop   
    call funcs(10**x*energyx, fx)
    f(1)=fx*cfs*10**x
    !print*, "f1, fx,cfs, x, energyx=", f(1), fx, cfs, x, energyx
    end subroutine
end subroutine
