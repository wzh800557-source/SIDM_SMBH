

subroutine fx_g(x,fx)
	use model_basic
	implicit none
	real(8) x,fx, rpar(10), alpha
	real(8) v0, n0
	integer ipar(10)
	alpha=cct_share%alpha
	if(x<=0d0)then
		fx=fgx_g0*exp(x)
	else
		fx=fgx_g0*(x/ctl%x_boundary)**(alpha-1.5d0)
	end if
!    print*, alpha, ctl%x_boundary
!	fx=fx*(2*pi*v0**2)**(-3/2d0)*n0
	!print*, "x=",x
	!read(*,*)
end subroutine

subroutine fgx_mb_star(x, fx)
	use model_basic
	implicit none
	real(8) x, fx
	!real(8) yout(6)
	!print*, fc_share%xb(fc_share%nbin)
	!read(*,*)
	if(x>=10**fc_share%xb(1).and.x<=10**fc_share%xb(fc_share%nbin))then
		!call get_value_at_x_fc(fc_share, log10(x), yout, 1)
		call fc_share%get_value_l(log10(x),fx)
		!fx=yout(1)
	else if(x<=0) then
		fx=exp(x)*fgx_g0
	else if(x>10**fc_share%xb(fc_share%nbin).or.(x.ge.0.and.x<10**fc_share%xb(1)))then
		fx=0d0		
	end if
	!call fc_share%print()
	!print*, "x, fx=", x, fx
	!read(*,*)
end subroutine

subroutine fgx_mb(x, fx)
	use model_basic
	implicit none
	real(8) x, fx
	!real(8) yout(6)
    real(8) yout
	!if(fgx_g0.eq.0)then
    !    print*, "error! fgx_g0 =0"
    !    stop
    !end if
	if(x>ctl%x_boundary.and.x<=emax_factor)then
		!call get_value_at_x_fc(fc_share, log10(x), yout,1)
        call linear_int(fc_share%xb(1:fc_share%nbin),fc_share%fx(1:fc_share%nbin),&
            fc_share%nbin, log10(x),yout)
		fx=yout
	else if(x<=0) then
		fx=exp(x)*fgx_g0
	!else if(x>10**fc_share%xb(fc_share%nbin).or.(x<=10**fc_share%xb(1).and.x>0))then
    else if(x>emax_factor.or.(x<=ctl%x_boundary.and.x>0))then
		fx=0d0
	end if
	!call fc_share%print()
	!if(x>2)then
		!print*, "x, fx=", x, fx, ctl%x_boundary
		!read(*,*)
	!end if
end subroutine


subroutine fgx_mb_direct(x, fx)
	use model_basic
	implicit none
	real(8) x, fx
	real(8) yout(6)
	integer idx
	!if(fgx_g0.eq.0)then
    !    print*, "error! fgx_g0 =0"
    !    stop
    !end if
	if(x>ctl%x_boundary.and.x<=emax_factor)then
		call return_idx(log10(x),fc_share%xmin,fc_share%xmax,fc_share%nbin,idx,1)
		fx=fc_share%fx(idx)
		!fx=yout(1)
	else if(x<=0) then
		fx=exp(x)*fgx_g0
	!else if(x>10**fc_share%xb(fc_share%nbin).or.(x<=10**fc_share%xb(1).and.x>0))then
    else if(x>emax_factor.or.(x<=ctl%x_boundary.and.x>0))then
		fx=0d0
	end if
	!call fc_share%print()
	!print*, "x, fx=", x, fx
	!read(*,*)
end subroutine

subroutine output_coef_ba16(fn)
	use com_main_gw
	use md_coeff
	implicit none
	character*(*) fn
	integer i

	open(unit=999,file=trim(adjustl(fn))//"/coef_ba16.bin",access='stream', form='unformatted')
	do i=1, ctl%m_bins
		call df(i)%write_grid(999)
	end do
	call df_tot%write_grid(999)
	close(999)	
end subroutine

subroutine input_coef_ba16(fn)
	use com_main_gw
	use md_coeff
	implicit none
	character*(*) fn
	integer i

	open(unit=999,file=trim(adjustl(fn))//"/coef_ba16.bin",access='stream', form='unformatted',status='old')
	allocate(df(ctl%m_bins))
	do i=1, ctl%m_bins
		call df(i)%read_grid(999)
	end do
	call df_tot%read_grid(999)
	close(999)	
	
end subroutine


subroutine get_coeff_sigma_funcs_cfs_grid(energy,jum,barg, sigma110, &
    sigma111,sigma131,sigma130,sigma13_1,sigma330,sigma310)
use md_coeff
use md_cfuns
implicit none
real(8) energy,jum
real(8) sigma110, sigma130,sigma111,sigma310, sigma13_1,sigma330, sigma131
external::barg

call get_sigma_funcs_cfs_grid(energy,jum, barg, cfs%cfs_110,cfs%nj,cfs%ns,&
    cfs%jmin,cfs%jmax, sigma110)
!print*, "e,j,sigma110=", energy, jum, sigma110
call get_sigma_funcs_cfs_grid(energy,jum, barg, cfs%cfs_111,cfs%nj,cfs%ns,&
    cfs%jmin,cfs%jmax, sigma111)
!print*, "e,j,sigma111=", energy, jum, sigma111
call get_sigma_funcs_cfs_grid(energy,jum, barg, cfs%cfs_131,cfs%nj,cfs%ns,&
    cfs%jmin,cfs%jmax, sigma131)
!print*, "e,j,sigma131=", energy, jum, sigma131
call get_sigma_funcs_cfs_grid(energy,jum, barg, cfs%cfs_130,cfs%nj,cfs%ns,&
    cfs%jmin,cfs%jmax, sigma130)
!print*, "e,j,sigma130=", energy, jum, sigma130
call get_sigma_funcs_cfs_grid(energy,jum, barg, cfs%cfs_13_1,cfs%nj,cfs%ns,&
    cfs%jmin,cfs%jmax, sigma13_1)
!print*, "e,j,sigma13_1=", energy, jum, sigma13_1
!print*, cfs%cfs_330
!read(*,*)
call get_sigma_funcs_cfs_grid(energy,jum, barg, cfs%cfs_330,cfs%nj,cfs%ns,&
    cfs%jmin,cfs%jmax, sigma330)
!print*, "e,j,sigma330=", energy, jum, sigma330
call get_sigma_funcs_cfs_grid(energy,jum, barg, cfs%cfs_310,cfs%nj,cfs%ns,&
    cfs%jmin,cfs%jmax, sigma310)
!print*, "e,j,sigma330=", energy, jum, sigma310
end subroutine


subroutine get_coeff_xr(coejum,rj, sigma0,sigma110, &
    sigma111,sigma131,sigma130,sigma13_1,sigma330,sigma310)
	use md_coeff
    use md_cfuns
	implicit none
	type(coeff_type)::coe,coejum
	real(8) energy,jum,jc,rj
	real(8) sigma0, sigma110, sigma130,sigma111,sigma310, sigma13_1,sigma330, sigma131
	external::barg
    integer i

	coe%e_110=sigma110
	coe%e_0=-sigma0
	coe%ee=4d0/3d0*(sigma0+sigma13_1)
!	coe%j=1d0/jum*( (5-3*jum**2)/12d0*sigma0-jum**2*(m+ma)/(2*ma)*sigma111+sigma310-sigma330/3d0)
	coe%j_111=rj*(sigma110-sigma111)
	coe%j_rest=(5-10*rj)/3d0*sigma0+4*sigma310-4d0/3d0*sigma330+rj/2d0*sigma131&
        -3d0/2d0*rj*sigma111-4d0/3d0*rj*sigma130

	coe%jj= 10/3d0*(rj-rj**2)*sigma0+2*rj**2*sigma131-2*rj**2*sigma111+8*rj*sigma310-8d0/3d0*rj*sigma330&
        +rj**2*4d0/3d0*sigma13_1-8d0/3d0*rj**2*sigma130
	coe%ej= 4d0/3d0*rj*(sigma13_1-sigma130)
	coejum=coe

end subroutine

subroutine get_coeff_xj(coejum,j, sigma0,sigma110, &
    sigma111,sigma131,sigma130,sigma13_1,sigma330,sigma310)
	use md_coeff
    use md_cfuns
	implicit none
	type(coeff_type)::coe,coejum
	real(8) j,rj
	real(8) sigma0, sigma110, sigma130,sigma111,sigma310, sigma13_1,sigma330, sigma131
	rj=j**2
	coe%e_110=sigma110
	coe%e_0=-sigma0
	coe%ee=4d0/3d0*(sigma0+sigma13_1)
!	coe%j=1d0/jum*( (5-3*jum**2)/12d0*sigma0-jum**2*(m+ma)/(2*ma)*sigma111+sigma310-sigma330/3d0)
	coe%j_111=j*sigma111
	coe%j_rest=1d0/j*(5*(1-3*rj)/12d0*sigma0+sigma310-1d0/3d0*sigma330+rj/2d0*sigma111&
        -1d0/3d0*rj*sigma130-1d0/6d0*rj*sigma13_1)

	coe%jj=5/6d0*(1-rj)*sigma0+rj/2d0*sigma131-rj/2d0*sigma111+2*sigma310-2d0/3d0*sigma330&
        +rj/3d0*sigma13_1-2d0/3d0*rj*sigma130
	coe%ej=2d0/3d0*j*(sigma13_1-sigma130)
	coejum=coe

end subroutine


subroutine get_coeff_ej(coejum,jum, sigma0,sigma110, &
    sigma111,sigma131,sigma130,sigma13_1,sigma330,sigma310)
	use md_coeff
	implicit none
	type(coeff_type)::coe,coejum
	real(8) jum
    real(8) sigma0, sigma110, sigma130,sigma111,sigma310, sigma13_1,sigma330, sigma131
	
	coe%e_110=sigma110
	coe%e_0=-sigma0
	coe%ee=4d0/3d0*(sigma0+sigma13_1)
	coe%j_111=-jum*sigma111
	coe%j_rest=1d0/jum*((5-3*jum**2)/12d0*sigma0+sigma310-sigma330/3d0)
	coe%jj= (5-3*jum**2)/6d0*sigma0+jum**2/2d0*sigma131-jum**2/2d0*sigma111+2*sigma310-2d0/3d0*sigma330
	coe%ej= -2d0/3d0*jum*(sigma0+sigma130)
	coejum=coe

end subroutine

subroutine get_coenr(even, evjum, m, en, jc, coenr,idx,idy)
	use com_main_gw
	implicit none
	real(8) m
	type(coeff_type)::coenr
	real(8) even, evjum
	real(8) dj_111(20), de_110(20)
	integer idx, idy
	real(8) jc, de_0, dee,dj_rest,djj,dej, en, jmin, jmax, rdx,rdy
	integer i
	
	
	select case(ctl%jbin_type)
	case(jbin_type_lin)
		!
	case(jbin_type_log)
		!print*, "sample%jm=",sample%jm, evjum
		evjum=log10(evjum)
		!print*, "evjum=",evjum
	case(jbin_type_sqr)
		evjum=evjum**2
	end select
	select case(ctl%method_interpolate)
	case(method_int_nearst)
		if(ctl%grid_type.eq.sts_type_grid)then
			rdx=(even-dms%emin)/dc_grid_xstep
			rdy=(evjum-dms%jmin)/dc_grid_ystep
			idx=nint(rdx)+1
			idy=nint(rdy)+1
		end if
		
		if(idx<-1)then
			idx=1
			print*, "even=",even, dms%emin
		end if
		if(idx>dms%nbin_grid)then
			idx=dms%nbin_grid
			print*, "even=",even, dms%emax
		end if
		if(idy<-1)then
			idy=1
			print*, "evjum=",evjum, jmin
			stop
		end if
		if(idy>dms%nbin_grid)then
			idy=dms%nbin_grid
		end if
		do i=1, dms%n
			!evjum=dms%dc0%s2_de_110%ymin
			
			de_110(i)=dms%mb(i)%dc%s2_de_110%fxy(idx,idy)
			dj_111(i)=dms%mb(i)%dc%s2_dj_111%fxy(idx,idy)
		end do
		associate(dc0=>dms%dc0)
		!    dc0=>dms%dc0
			de_0=dc0%s2_de_0%fxy(idx,idy)
			dee=dc0%s2_dee%fxy(idx,idy)
			dj_rest=dc0%s2_dj_rest%fxy(idx,idy)
			djj=dc0%s2_djj%fxy(idx,idy)
			dej=dc0%s2_dej%fxy(idx,idy)
		end associate
	case(method_int_linear)
		do i=1, dms%n
			!print*, "dms%nbin_grid=",dms%nbin_grid
			call linear_int_2d(dms%emin,dms%jmin,dms%nbin_grid,dms%nbin_grid,dc_grid_xstep,dc_grid_ystep,&
				dms%mb(i)%dc%s2_de_110%fxy, even,evjum,de_110(i))
			call linear_int_2d(dms%emin,dms%jmin,dms%nbin_grid,dms%nbin_grid,dc_grid_xstep,dc_grid_ystep,&
				dms%mb(i)%dc%s2_dj_111%fxy, even,evjum,dj_111(i))
		end do
		associate(dc0=>dms%dc0)
			call linear_int_2d(dms%emin,dms%jmin,dms%nbin_grid,dms%nbin_grid,dc_grid_xstep,dc_grid_ystep,&
				dc0%s2_de_0%fxy, even,evjum,de_0)
			call linear_int_2d(dms%emin,dms%jmin,dms%nbin_grid,dms%nbin_grid,dc_grid_xstep,dc_grid_ystep,&
				dc0%s2_dee%fxy, even,evjum,dee)
			call linear_int_2d(dms%emin,dms%jmin,dms%nbin_grid,dms%nbin_grid,dc_grid_xstep,dc_grid_ystep,&
				dc0%s2_dj_rest%fxy, even,evjum,dj_rest)
			call linear_int_2d(dms%emin,dms%jmin,dms%nbin_grid,dms%nbin_grid,dc_grid_xstep,dc_grid_ystep,&
				dc0%s2_djj%fxy, even,evjum,djj)
			call linear_int_2d(dms%emin,dms%jmin,dms%nbin_grid,dms%nbin_grid,dc_grid_xstep,dc_grid_ystep,&
				dc0%s2_dej%fxy, even,evjum,dej)

			end associate
	end select


	coeNr%jj=djj*jc**2; 
	coeNr%e=de_0; coeNr%j=dj_rest
	do i=1, dms%n
		coeNr%e=coeNr%e+m/dms%mb(i)%mc*de_110(i)
		coeNr%j=coeNr%j+dj_111(i)*(m+dms%mb(i)%mc)/dms%mb(i)%mc/2d0; 
	end do
	coeNr%ee=dee*en*en;
	!coeNr%ee=dee*(10**even*ctl%energy0)**2

	coeNr%e=coeNr%e*en; 
	!coeNr%e=coeNr%e*(10**even*ctl%energy0)

	coeNR%j=  coeNR%j*jc
	coeNr%ej=  dej*en*jc

	if(ctl%chattery.ge.3)then
		print*, "=========get cej NR==============="
		print*, "enev, evjum, en, jc=", even, evjum, en, jc
		print*, "idx,idy=",idx,idy
		print*, "de_0, dee=", de_0, dee
		print*, "coenr%e, ee, ej=", coenr%e, coenr%ee,coenr%ej
		print*, "coenr%j, jj=", coenr%j, coenr%jj
		print*, "=========end of get cej NR========"
	end if
end subroutine



subroutine get_coeff_sigma_funcs_cfs_rk(energy,jum,barg, sigma110, &
	sigma111,sigma131,sigma130,sigma13_1,sigma330,sigma310)
use md_coeff
use md_cfuns
implicit none
real(8) energy,jum
real(8) sigma110, sigma130,sigma111,sigma310, sigma13_1,sigma330, sigma131
external::barg

call get_sigma_funcs_cfs_rk(energy,jum, barg, cfs%cfs_110,cfs%nj,cfs%ns,&
	cfs%jmin,cfs%jmax, sigma110)

call get_sigma_funcs_cfs_rk(energy,jum, barg, cfs%cfs_111,cfs%nj,cfs%ns,&
	cfs%jmin,cfs%jmax, sigma111)
call get_sigma_funcs_cfs_rk(energy,jum, barg, cfs%cfs_131,cfs%nj,cfs%ns,&
	cfs%jmin,cfs%jmax, sigma131)
call get_sigma_funcs_cfs_rk(energy,jum, barg, cfs%cfs_130,cfs%nj,cfs%ns,&
	cfs%jmin,cfs%jmax, sigma130)

call get_sigma_funcs_cfs_rk(energy,jum, barg, cfs%cfs_13_1,cfs%nj,cfs%ns,&
	cfs%jmin,cfs%jmax, sigma13_1)
!print*, cfs%cfs_330
!read(*,*)
call get_sigma_funcs_cfs_rk(energy,jum, barg, cfs%cfs_330,cfs%nj,cfs%ns,&
	cfs%jmin,cfs%jmax, sigma330)
call get_sigma_funcs_cfs_rk(energy,jum, barg, cfs%cfs_310,cfs%nj,cfs%ns,&
	cfs%jmin,cfs%jmax, sigma310)
end subroutine