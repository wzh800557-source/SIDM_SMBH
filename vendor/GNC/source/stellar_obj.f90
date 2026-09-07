module md_stellar_object
    !use com_sts_type
    use md_bk_species
    use md_coeff
	use md_chain_pointer
    !use md_sts_d21
    use com_sts_type
    use, intrinsic :: ieee_arithmetic
    type nejw_type
        real(8) e, j, w
        integer idx
    end type
    type dms_stellar_object
        integer n
        real(8) n_real    
        type(nejw_type),allocatable::nejw(:)
        type(s1d_type)::barge!, barge_norm
        !type(s1d_type)::gj
        type(s2d_hst_type)::nxj
        type(s2d_type)gxj
        type(s1d_type)::fden_simu, fden
        type(s1d_type)::fNa, fMa
        real(8) asymp
        contains
            procedure::init=>init_dms_stellar_object
            procedure::deallocation=>deallocate_dms_stellar_object            
    end type
	type dso_pointer
		type(dms_stellar_object),pointer::p
	end type
    private::init_dms_stellar_object,deallocate_dms_stellar_object
    integer,parameter::dj_n=8, dj_n2=8
contains
    subroutine init_dms_stellar_object(so, nbin_grid, nbin_gx,&
        emin, emax, tmin, tmax, rh, jb_type)
        implicit none
        class(dms_stellar_object)::so
        integer nbin_grid, nbin_gx, i
        real(8) emin, emax, tmin, tmax, rh
        real(8) rmin, rmax
		integer jb_type

        rmin=log10(0.5d0*rh/(10**emax));rmax=log10(0.5d0*rh/(10**emin))

        call so%barge%init(emin, emax, nbin_gx,sts_type_dstr)
        call so%barge%set_range()
		
        call so%fden%init(rmin, rmax, nbin_gx, sts_type_dstr)
        call so%fden%set_range()
        call so%fden_simu%init(rmin, rmax, nbin_gx, sts_type_dstr)
        call so%fden_simu%set_range()

        call so%fNa%init(rmin, rmax, nbin_gx, sts_type_dstr)
        call so%fNa%set_range()
        call so%fMa%init(rmin, rmax, nbin_gx, sts_type_dstr)
        call so%fMa%set_range()

        call so%gxj%init(nbin_gx, nbin_gx, emin, emax, tmin,tmax, sts_type_dstr)
        call so%gxj%set_range()
		call so%nxj%init(nbin_gx, nbin_gx, emin, emax, tmin,tmax, use_weight=.true.)
        call so%nxj%set_range()
		!call so%gxjcr%init(nbin_gx, nbin_gx, emin, emax, tmin,tmax, sts_type_dstr)
		!call so%gxjcr%set_range()

        so%n=0
        so%n_real=0
    end subroutine
    subroutine deallocate_dms_stellar_object(so)
        implicit none
        class(dms_stellar_object)::so
        call so%fden%deallocate()
        call so%fna%deallocate()
        call so%fma%deallocate()
        call so%barge%deallocate()
    end subroutine

    subroutine dms_so_get_nxj_from_nejw(so, jbtype)
        implicit none
        class(dms_stellar_object)::so
        integer ntasks, jbtype
        real(8) en(so%n), jm(so%n),we(so%n)
        if(so%n>0)then
            en(1:so%n)=so%nejw(1:so%n)%e
            jm(1:so%n)=so%nejw(1:so%n)%j
            we(1:so%n)=so%nejw(1:so%n)%w
			select case(jbtype)
			case(Jbin_type_lin)
				call so%nxj%get_stats_weight(en, jm, we, so%n)
			case(Jbin_type_log)
				call so%nxj%get_stats_weight(en, log10(jm), we, so%n)
			case(Jbin_type_sqr)
				call so%nxj%get_stats_weight(en, jm**2d0, we, so%n)
			case default
				print*, "dms_nxj_newj:error! define jbtype", jbtype
				stop
			end select
			
			so%n_real=sum(so%nejw(1:so%n)%w)
		else
			so%n_real=0
			return
		end if
    end subroutine

    subroutine get_asymp_norm_factor_one(dso,  x_boundary, norm)
        implicit none
        class(dms_stellar_object)::dso
		type(s1d_type)::bg
        real(8) cnorm, cn, norm,x_boundary

        if(dso%n_real>0)then
            
			call dso%barge%get_value_l(log10(x_boundary),cnorm)

            if( isnan(cnorm).or. cnorm.eq.0)then
                print*, "star:cnorm is nan or 0", cnorm
                call dso%barge%print("dso%barge")
                stop
            else
                norm=dso%asymp/cnorm
            end if
        else
            print*, "error! dso%n_real=", dso%n_real
            stop
        end if
    end subroutine
    subroutine normalize_barge_one(dso, norm)
        implicit none
        class(dms_stellar_object)::dso
        real(8) norm
        integer i
        
        dso%barge%fx=dso%barge%fx*norm
        
    end subroutine
	subroutine normalize_gxj_one(dso, norm)
        implicit none
        class(dms_stellar_object)::dso
        real(8) norm
        integer i
        
        dso%gxj%fxy=dso%gxj%fxy*norm
        
    end subroutine
    subroutine dms_so_get_fxj(so, n0, mbh, v0,jbtype)
        implicit none
        class(dms_stellar_object)::so
        integer i, j,jbtype
        real(8) jm ,x, n0, mbh, v0
		if(so%n.eq.0) return
		select case(jbtype)
		case(Jbin_type_lin)
			do i=1, so%nxj%nx
			x=10**so%nxj%xcenter(i)
				do j=1, so%nxj%ny
					jm=so%nxj%ycenter(j)
					so%gxj%fxy(i,j)=so%nxj%nxyw(i,j)/(x*log(10d0))&
					/so%nxj%xstep/so%nxj%ystep &
					*pi**(-1.5d0)*v0**6*x**2.5d0/jm/n0/mbh**3
					!print*, "xystep=",so%nxj%xstep, so%nxj%ystep
					!read(*,*)
				end do 
			end do
		case(Jbin_type_log)
			do i=1, so%nxj%nx
				x=10**so%nxj%xcenter(i)
				do j=1, so%nxj%ny
					jm=10**so%nxj%ycenter(j)
					so%gxj%fxy(i,j)=so%nxj%nxyw(i,j)/(x*log(10d0))&
					/so%nxj%xstep/so%nxj%ystep &
					*pi**(-1.5d0)*v0**6*x**2.5d0/(jm**2*log(10d0))/n0/mbh**3
					!print*, "xystep=",so%nxj%xstep, so%nxj%ystep
					!read(*,*)
				end do 
			end do
		case(Jbin_type_sqr)
			do i=1, so%nxj%nx
				x=10**so%nxj%xcenter(i)
				do j=1, so%nxj%ny
					jm=so%nxj%ycenter(j)**0.5d0
					so%gxj%fxy(i,j)=so%nxj%nxyw(i,j)/(x*log(10d0))&
					/so%nxj%xstep/so%nxj%ystep &
					*pi**(-1.5d0)*v0**6*x**2.5d0/2d0/n0/mbh**3
					!print*, "xystep=",so%nxj%xstep, so%nxj%ystep
					!read(*,*)
				end do 
			end do
		case default
			print*, "fxj error!"
			stop
		end select
    end subroutine
    subroutine get_barge_stellar(so,jbtype)
        implicit none
        class(dms_stellar_object)::so
        !real(8) sums
        integer i, j, idid, jbtype
        real(8) x,int_out
		if(so%n.eq.0)return
        do i=1, so%barge%nbin
            int_out=0
            !print*,"xmax=",mb%barge%xmax
            if(so%barge%nbin.ne.so%gxj%nx)then
                print*, "error! barge%nbin should = gxj%nx"
                stop
            end if
            so%barge%xb(i)=so%gxj%xcenter(i)
			select case(jbtype)
			case (Jbin_type_lin)
				do j=1, so%gxj%ny
					!==original==
					int_out=int_out+so%gxj%fxy(i,j)*so%gxj%ycenter(j)*so%gxj%ystep*2d0
					!===test=====
					!int_out=int_out+so%gxj%fxy(i,j)*so%gxj%ycenter(j)*so%gxj%ystep&
					!    /(1-so%gxj%ycenter(j)**2)**0.5d0
				end do
				so%barge%fx(i)=int_out
				!so%barge%nsam=so%n
				if(isnan(so%barge%fx(i)))then
					print*, "get_barge_stellar:fx is NaN:", int_out
					call so%gxj%print()
					stop "NaN detected in barge"
				end if
			case(Jbin_type_log)
				do j=1, so%gxj%ny
					!==original==
					int_out=int_out+so%gxj%fxy(i,j)*(10**so%gxj%ycenter(j))**2*so%gxj%ystep*2d0*log(10d0)
					!===test=====
					!int_out=int_out+so%gxj%fxy(i,j)*so%gxj%ycenter(j)*so%gxj%ystep&
					!    /(1-so%gxj%ycenter(j)**2)**0.5d0
				end do
				so%barge%fx(i)=int_out
				!so%barge%nsam=so%n
				if(isnan(so%barge%fx(i)))then
					print*, "get_barge_stellar:fx is NaN:", int_out
					call so%gxj%print()
					stop "NaN detected in barge"
				end if
			case(Jbin_type_sqr)
				do j=1, so%gxj%ny
					!==original==
					int_out=int_out+so%gxj%fxy(i,j)*so%gxj%ycenter(j)*so%gxj%ystep*2d0
					!===test=====
					!int_out=int_out+so%gxj%fxy(i,j)*so%gxj%ycenter(j)*so%gxj%ystep&
					!    /(1-so%gxj%ycenter(j)**2)**0.5d0
				end do
				so%barge%fx(i)=int_out
				!so%barge%nsam=so%n
				if(isnan(so%barge%fx(i)))then
					print*, "get_barge_stellar:fx is NaN:", int_out
					call so%gxj%print()
					stop "NaN detected in barge"
				end if
			end select	
        end do
    end subroutine
    subroutine get_dens(so, n0, v0, rh, emin,weight_asym)
        implicit none
        class(dms_stellar_object)::so
        real(8)  n0, v0, rh, emin, weight_asym
        integer i
        real(8) en(so%n), jm(so%n), we(so%n)

		if(so%n>0)then

			call get_fden(so%barge, so%fden, n0, &
			v0, rh, emin)
			!call so%barge%print("barge:after get_fden")
			!call so%fden%print("get_dens:after get_fden")
			!print*, "get_dens:finished get_fden!" 
			en(1:so%n)=so%nejw(1:so%n)%e
            jm(1:so%n)=so%nejw(1:so%n)%j
            we(1:so%n)=so%nejw(1:so%n)%w
			call get_fden_sample_particle(en(1:so%n), jm(1:so%n), we(1:so%n), so%n, so%fden_simu)
            !call so%fden_simu%print("fden_simu")
			call get_fna(so%fden_simu, so%fNa)
            !call so%fna%print("fna")
			!print*, "emin=",emin
			!print*, "r0=", log10(rh)
			!call so%fden%print("fden")
			!call so%fden_simu%print("fden_simu")
			!call so%barge%print("barge")
			!read(*,*)
		else
			so%fden_simu%fx=0
		end if
    end subroutine



end module