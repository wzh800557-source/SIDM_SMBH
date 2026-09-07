
module md_mass_bins 
	use md_stellar_object
	integer,parameter:: n_tot_comp=7
	type mass_bins
		type(dms_stellar_object):: bstar, star, sbh, bbh, ns, all, wd, bd
		type(dso_pointer)::dsp(n_tot_comp)

		real(8) mc, m1, m2 !m1<mass<m2
		integer nbin_grid, nbin_gx
		real(8) frac   
		real(8) emin, emax, jmin, jmax, mbh, v0, n0, barmin, rh
		type(diffuse_coeffient_type)::dc
		!type(s2d_type)::dt
		contains
		procedure::init=>init_mass_bins
		procedure::write_mb=>write_info_mass_bin
		procedure::read_mb=>read_info_mass_bin
		!procedure::get_nEJ0_fEJ0
		!procedure::get_barge
	end type
contains
	subroutine init_mass_bins(mb, nbin_grid, nbin_gx, emin, emax, jmin, jmax,&
			mbh, v0, n0, rh,jb_type)
		implicit none
		class(mass_bins),target::mb
		integer nbin_grid, nbin_gx
		real(8) emin, emax, jmin, jmax, mbh, n0, v0, rmin, rmax, rh
		integer jb_type, i

		mb%nbin_grid=nbin_grid; mb%nbin_gx=nbin_gx
		mb%emin=emin; mb%emax=emax; mb%jmin=jmin; mb%jmax=jmax
		mb%mbh=mbh; mb%v0=v0; mb%n0=n0; mb%rh=rh

		rmin=log10(0.5d0*rh/(10**emax));rmax=log10(0.5d0*rh/(10**emin))

		call mb%all%init(nbin_grid, nbin_gx,emin, emax, jmin, jmax, rh, jb_type)

		
		mb%dsp(1)%p=>mb%star
		mb%dsp(2)%p=>mb%sbh
		mb%dsp(3)%p=>mb%ns
		mb%dsp(4)%p=>mb%wd
		mb%dsp(5)%p=>mb%bd
		mb%dsp(6)%p=>mb%bstar
		mb%dsp(7)%p=>mb%bbh

		do i=1, n_tot_comp
			call mb%dsp(i)%p%init(nbin_grid, nbin_gx,emin, emax, jmin, jmax, rh, jb_type)
		end do

	end subroutine

    subroutine write_info_mass_bin(mb, funit)
        implicit none
        class(mass_bins)::mb
        integer funit
        call mb%dc%write_grid(funit)        
        write(funit) mb%mc, mb%m1, mb%m2, mb%nbin_grid, mb%nbin_gx
        write(funit) mb%emin, mb%emax, mb%jmin, mb%jmax, &
             mb%mbh, mb%v0, mb%n0, mb%rh
		write(funit) mb%all%n,mb%star%n, mb%sbh%n, mb%wd%n, mb%ns%n, mb%bd%n
		write(funit) mb%bstar%n, mb%bbh%n	 
        write(funit) mb%all%barge, mb%all%fden,mb%all%fden_simu, mb%all%asymp
		if(mb%star%n>0)then
			write(funit) mb%star%fden, mb%star%fden_simu, mb%star%barge,mb%star%asymp
		end if
		if(mb%sbh%n>0)then
			write(funit) mb%sbh%fden, mb%sbh%fden_simu, mb%sbh%barge,mb%sbh%asymp
		end if
		if(mb%ns%n>0)then
			write(funit) mb%ns%fden, mb%ns%fden_simu, mb%ns%barge,mb%ns%asymp
		end if
		if(mb%wd%n>0)then
			write(funit) mb%wd%fden, mb%wd%fden_simu, mb%wd%barge,mb%wd%asymp
		end if
		if(mb%bd%n>0)then
			write(funit) mb%bd%fden, mb%bd%fden_simu, mb%bd%barge,mb%bd%asymp
		end if
		if(mb%bstar%n>0)then
			write(funit) mb%bstar%fden, mb%bstar%fden_simu, mb%bstar%barge,mb%bstar%asymp
		end if
		if(mb%bbh%n>0)then
			write(funit) mb%bbh%fden, mb%bbh%fden_simu, mb%bbh%barge,mb%bbh%asymp
		end if
		

    end subroutine
    subroutine read_info_mass_bin(mb, funit)
        implicit none
        class(mass_bins)::mb
        integer funit
        call mb%dc%read_grid(funit)        
        read(funit) mb%mc, mb%m1, mb%m2, mb%nbin_grid, mb%nbin_gx
        read(funit) mb%emin, mb%emax, mb%jmin, mb%jmax, &
            mb%mbh, mb%v0, mb%n0, mb%rh
		read(funit) mb%all%n,mb%star%n, mb%sbh%n, mb%wd%n, mb%ns%n, mb%bd%n
		read(funit) mb%bstar%n, mb%bbh%n	 
		read(funit) mb%all%barge, mb%all%fden,mb%all%fden_simu, mb%all%asymp
		if(mb%star%n>0)then
			read(funit) mb%star%fden, mb%star%fden_simu, mb%star%barge,mb%star%asymp
		end if
		if(mb%sbh%n>0)then
			read(funit) mb%sbh%fden, mb%sbh%fden_simu, mb%sbh%barge,mb%sbh%asymp
		end if
		if(mb%ns%n>0)then
			read(funit) mb%ns%fden, mb%ns%fden_simu, mb%ns%barge,mb%ns%asymp
		end if
		if(mb%wd%n>0)then
			read(funit) mb%wd%fden, mb%wd%fden_simu, mb%wd%barge,mb%wd%asymp
		end if
		if(mb%bd%n>0)then
			read(funit) mb%bd%fden, mb%bd%fden_simu, mb%bd%barge,mb%bd%asymp
		end if
		if(mb%bstar%n>0)then
			read(funit) mb%bstar%fden, mb%bstar%fden_simu, mb%bstar%barge,mb%bstar%asymp
		end if
		if(mb%bbh%n>0)then
			read(funit) mb%bbh%fden, mb%bbh%fden_simu, mb%bbh%barge,mb%bbh%asymp
		end if

        !read(funit) mb%all%barge, mb%all%fden, mb%star%fden, mb%sbh%fden, &
        !    mb%bstar%fden, mb%bbh%fden, mb%ns%fden, mb%wd%fden
        !read(funit) mb%all%fden_simu, mb%star%fden_simu, mb%sbh%fden_simu, &
        !    mb%bstar%fden_simu, mb%bbh%fden_simu, mb%ns%fden_simu, mb%wd%fden_simu
        !read(funit) mb%star%barge, mb%sbh%barge, mb%bbh%barge, mb%ns%barge, &
		!	mb%wd%barge
        !read(funit) mb%all%asymp,mb%star%asymp,mb%bstar%asymp,mb%sbh%asymp, &
        !    mb%bbh%asymp, mb%ns%asymp, mb%wd%asymp

    end subroutine
	subroutine get_dens_mb(mb, n0, v0, rh, emin, weight_asym)
		implicit none
		type(mass_bins)::mb
		real(8) n0, v0, rh, emin, weight_asym
		integer i
		do i=1, n_tot_comp
			call get_dens(mb%dsp(i)%p, n0, v0, rh, emin,weight_asym)
		end do
	end subroutine	
	subroutine get_fMa_mb(mb)
		implicit none
		type(mass_bins)::mb
		integer i
		do i=1, n_tot_comp
			if(mb%dsp(i)%p%n>0)then
				call get_fMa(mb%dsp(i)%p%fNa, mb%dsp(i)%p%fMa, mb%mc)
			else
				mb%dsp(i)%p%fMa%fx=0
			end if
		end do
		
	end subroutine	
	subroutine get_barge_stellar_mb(mb, jbtype)
		implicit none
		type(mass_bins)::mb
		integer i, jbtype

		do i=1, n_tot_comp
			call get_barge_stellar(mb%dsp(i)%p,jbtype)
		end do
		call get_barge_stellar(mb%all,jbtype)
	end subroutine
end module
module md_dms
	use md_mass_bins
    
   !integer,parameter::n_massbin=2
   type diffuse_mspec
        integer n
       type(mass_bins),allocatable:: mb(:)
       type(mass_bins):: all
       real(8) weight_asym
       integer idx_ref
       !type(s1d_type)::barge0, barge_bh, barge_star
       type(diffuse_coeffient_type)::dc0
       integer nbin_grid, nbin_gx,jbin_type,grid_type
       real(8) emin, emax, jmin, jmax, mbh, v0, n0, rh, acmin, acmax
       real(8) x_boundary
       contains
       procedure:: set_diffuse_mspec
       procedure::init=>init_diffuse_mspec
       procedure::get_nxj
       procedure::get_fxj0
       procedure::get_dens0
       procedure::output_bin=>output_diffuse_mspec_bin
       procedure::input_bin=>input_diffuse_mspec_bin
       procedure::get_barge0
       procedure::normalize_barge
       procedure::get_asymp_norm_factor
       procedure::print_norm=>print_norm_dms
   end type
   private::init_diffuse_mspec, get_fxj0
   private::write_info_mass_bin,read_info_mass_bin,print_norm_dms
contains
   subroutine init_diffuse_mspec(dm, n)
        implicit none
        integer n, i
        class(diffuse_mspec)::dm

        if(allocated(dm%mb))then
            deallocate(dm%mb)
        end if
        allocate(dm%mb(n))
        dm%n=n

		call dm%dc0%init(dm%nbin_grid,  dm%emin, dm%emax, dm%jmin,dm%jmax,dm%grid_type)
        !call dm%dc0_bk%init(nbin_grid,  emin, emax, jmin,jmax)
        do i=1, dm%n
            call dm%mb(i)%dc%init(dm%nbin_grid,  dm%emin, dm%emax, dm%jmin,dm%jmax,dm%grid_type)
            !print*, "nbin_grid=", nbin_grid
            call dm%mb(i)%init(dm%nbin_grid, dm%nbin_gx, dm%emin, dm%emax, dm%jmin, dm%jmax, dm%mbh, &
			dm%v0, dm%n0, dm%rh, dm%jbin_type)
			!call dm%mb(i)%dt%init(nbin_gx,nbin_gx,emin, emax, jmin, jmax,sts_type_dstr)
        end do
		call dm%all%init(dm%nbin_grid, dm%nbin_gx, dm%emin, dm%emax, dm%jmin, dm%jmax, dm%mbh, &
		dm%v0, dm%n0, dm%rh, dm%jbin_type)      

    end subroutine


    subroutine set_diffuse_mspec(dm, nbin_grid, nbin_gx, emin, emax, jmin, jmax, mbh, &
         v0, n0, rh, xb, idx_ref,jb_type,grid_type)
        implicit none
        class(diffuse_mspec)::dm
        integer nbin_grid, nbin_gx, idx_ref,jb_type,grid_type
        real(8) emin, emax, jmin, jmax, mbh, v0, n0, rh
        real(8) acmin, acmax, xb
        dm%idx_ref=idx_ref
        dm%nbin_grid=nbin_grid; dm%nbin_gx=nbin_gx; dm%emin=emin; dm%emax=emax
        dm%mbh=mbh; 
        dm%v0=v0; dm%n0=n0; dm%rh=rh
        acmin=log10(0.5d0*rh/(10**emax)); acmax=log10(0.5d0*rh/(10**emin))
        dm%acmin=acmin;dm%acmax=acmax
        dm%x_boundary=xb
		dm%grid_type=grid_type
		
		call conv_j_space(jmin,dm%jmin,jb_type)
		call conv_j_space(jmax,dm%jmax,jb_type)

        dm%jbin_type=jb_type
		dm%grid_type=grid_type
        
    end subroutine

    subroutine set_mass_bin_mass_given(dm,masses, m1,m2, asym, n)
        implicit none
        class(diffuse_mspec) dm
        integer i, j, n
        real(8) masses(n), m1(n),m2(n), asym(n_tot_comp+1, n)
		!print*, "n=",n, n_tot_comp
        do i=1, n
			!print*, "i=",i
            dm%mb(i)%mc=masses(i)
            dm%mb(i)%m1=m1(i)
            dm%mb(i)%m2=m2(i)
            dm%mb(i)%all%asymp=asym(1, i)
			do j=1, n_tot_comp
				!print*, "j=",asym(j+1,i), asym(1,i)
				dm%mb(i)%dsp(j)%p%asymp=asym(j+1,i)*asym(1,i)
			end do
        end do
    end subroutine
	
    
    subroutine get_dc0(dm)
        implicit none
        integer n, nbin
        class(diffuse_mspec)::dm
    ! type(diffuse_coeffient_type):: dfs(n), dfs_tot
        integer i, j, k
        
        do i=1, dm%nbin_grid
            do j=1, dm%nbin_grid
                do k=1, dm%n
                    !dm%dc0%s2_de_110%fxy(i,j)=dm%dc0%s2_de_110%fxy(i,j)+dm%mb(k)%dc%s2_de_110%fxy(i,j)
                    dm%dc0%s2_de_0%fxy(i,j)=dm%dc0%s2_de_0%fxy(i,j)+dm%mb(k)%dc%s2_de_0%fxy(i,j)
                    dm%dc0%s2_dee%fxy(i,j)=dm%dc0%s2_dee%fxy(i,j)+dm%mb(k)%dc%s2_dee%fxy(i,j)
                    !dm%dc0%s2_dj_111%fxy(i,j)=dm%dc0%s2_dj_111%fxy(i,j)+dm%mb(k)%dc%s2_dj_111%fxy(i,j)
                    dm%dc0%s2_dj_rest%fxy(i,j)=dm%dc0%s2_dj_rest%fxy(i,j)+dm%mb(k)%dc%s2_dj_rest%fxy(i,j)
                    dm%dc0%s2_djj%fxy(i,j)=dm%dc0%s2_djj%fxy(i,j)+dm%mb(k)%dc%s2_djj%fxy(i,j)
                    dm%dc0%s2_dej%fxy(i,j)=dm%dc0%s2_dej%fxy(i,j)+dm%mb(k)%dc%s2_dej%fxy(i,j)
                end do
            end do
        end do        
    end subroutine

    subroutine get_fxj0(dm)
        implicit none
        class(diffuse_mspec)::dm
        integer i, maxidx(1)
        real(8) cnorm(6), norm

        do i=1, dm%n
            !print*, "dm%mb(i),i=",i
            call get_fxj(dm%mb(i),dm%jbin_type)
        end do
       ! print*, "get_fxj0:"
       ! print*, dm%mb(1)%star%nxj%nxyw(:,3)
       ! print*, dm%mb(1)%star%gxj%fxy(:,3)
        !print*, "dm all"
        call get_fxj(dm%all,dm%jbin_type)
       ! print*, dm%all%star%nxj%nxyw(:,3)
       ! print*, dm%all%star%gxj%fxy(:,3)
       ! stop
    end subroutine
    subroutine get_barge0(dm)
        implicit none
        class(diffuse_mspec)::dm
        integer i

        do i=1, dm%n
			call get_barge_stellar_mb(dm%mb(i),dm%jbin_type)
        end do
        call get_barge_stellar_mb(dm%all,dm%jbin_type)
        
    end subroutine
    
    subroutine get_asymp_norm_factor(dm)
        implicit none
        class(diffuse_mspec)::dm
        integer i
       
        call get_asymp_norm_factor_one(dm%mb(dm%idx_ref)%all,  &
			dm%x_boundary, dm%weight_asym)

    end subroutine
    subroutine normalize_barge(dm)
        implicit none
        class(diffuse_mspec)::dm
        real(8) starnorm, sbhnorm, bbhnorm, bstarnorm
        integer i,j
		do i=1, n_tot_comp
			dm%all%dsp(i)%p%barge%fx=0
			dm%all%all%barge%fx=0
		end do
        
        do i=1, dm%n
			call normalize_barge_one(dm%mb(i)%all,dm%weight_asym)
			do j=1,n_tot_comp
				call normalize_barge_one(dm%mb(i)%dsp(j)%p,dm%weight_asym)
				dm%all%dsp(j)%p%barge%fx=dm%all%dsp(j)%p%barge%fx+dm%mb(i)%dsp(j)%p%barge%fx
				dm%all%all%barge%fx=dm%all%all%barge%fx+dm%mb(i)%all%barge%fx
            end do
            
        end do
    
    end subroutine
	subroutine normalize_gxj(dm)
		implicit none
        class(diffuse_mspec)::dm
        real(8) starnorm, sbhnorm, bbhnorm, bstarnorm
        integer i,j
		do i=1, n_tot_comp
			dm%all%dsp(i)%p%gxj%fxy=0
			dm%all%all%gxj%fxy=0
		end do
        
        do i=1, dm%n
			call normalize_gxj_one(dm%mb(i)%all,dm%weight_asym)
			do j=1,n_tot_comp
				call normalize_gxj_one(dm%mb(i)%dsp(j)%p,dm%weight_asym)
				dm%all%dsp(j)%p%gxj%fxy=dm%all%dsp(j)%p%gxj%fxy+&
				dm%mb(i)%dsp(j)%p%gxj%fxy
				dm%all%all%gxj%fxy=dm%all%all%gxj%fxy+dm%mb(i)%all%gxj%fxy
            end do
            
        end do
	end subroutine
    subroutine print_norm_dms(dm,out_unit)
        implicit none 
        class(diffuse_mspec)::dm
        integer i,out_unit
        write(out_unit, fmt="(A20, F10.4)") "weight_asym=", dm%weight_asym
        write(out_unit, fmt="(2A10,10A13)") "N", "mc", "rNstar", "rNsbh", "rNbbh", "rNbstar","rNns", "rNwd", "rNbd"
        do i=1, dm%n
            write(out_unit, fmt="(I10, F10.2, 10F13.2)") i, dm%mb(i)%mc,&
                dm%mb(i)%star%n_real, dm%mb(i)%sbh%n_real, dm%mb(i)%bbh%n_real, dm%mb(i)%bstar%n_real, &
                dm%mb(i)%ns%n_real, dm%mb(i)%wd%n_real, dm%mb(i)%bd%n_real
        end do
		write(out_unit, fmt="(2A10,10A13)") "N", "mc", "sNstar", "sNsbh", "sNbbh", "sNbstar", "sNns", "sNwd", "sNbd"
        do i=1, dm%n
            write(out_unit, fmt="(I10, F10.2, 10I13)") i, dm%mb(i)%mc,&
                dm%mb(i)%star%n, dm%mb(i)%sbh%n, dm%mb(i)%bbh%n, dm%mb(i)%bstar%n, dm%mb(i)%ns%n, dm%mb(i)%wd%n, &
				 dm%mb(i)%bd%n
        end do
    end subroutine
    subroutine get_dens0(dm)
        implicit none
        class(diffuse_mspec)::dm
        integer i,j
        do i=1, dm%n
            !print*, "get_dns_mb"
			call get_dens_mb(dm%mb(i), dm%n0, dm%v0, dm%rh, dm%emin,dm%weight_asym)
			!call  dm%mb(i)%star%fden%print("dm%mb(i)%star%fden")
            !print*, "get_fma"
            call get_fma_mb(dm%mb(i))
            associate(mb=>dm%mb(i))
				mb%all%fden%fx=0; mb%all%fden_simu%fx=0
                mb%all%fma%fx=0
                do j=1, n_tot_comp
					mb%all%fden%fx=mb%all%fden%fx+mb%dsp(j)%p%fden%fx
					mb%all%fden_simu%fx=mb%all%fden_simu%fx+mb%dsp(j)%p%fden_simu%fx
                    mb%all%fma%fx=mb%all%fma%fx+mb%dsp(j)%p%fma%fx
				end do
                !print*, "fna"
                call get_fna(mb%all%fden_simu,mb%all%fna)
                !call get_fna(mb%all%fden_simu,mb%all%fna)
                !call mb%all%fden%print("all fden")
                !call mb%all%fden_simu%print("all fden_simu")
                !read(*,*)
            end associate
        end do


		associate(all=>dm%all)
			do j=1, n_tot_comp
				all%dsp(j)%p%fden%fx=0
				all%dsp(j)%p%fna%fx=0
				all%dsp(j)%p%fden_simu%fx=0
				all%dsp(j)%p%fma%fx=0

			end do    
			all%all%fden%fx=0
			all%all%fden_simu%fx=0
			all%all%fna%fx=0
			all%all%fma%fx=0
			do i=1, dm%n
				do j=1, n_tot_comp
					all%dsp(j)%p%fden%fx=all%dsp(j)%p%fden%fx+dm%mb(i)%dsp(j)%p%fden%fx
					all%dsp(j)%p%fden_simu%fx=all%dsp(j)%p%fden_simu%fx+dm%mb(i)%dsp(j)%p%fden_simu%fx
					all%dsp(j)%p%fma%fx=all%dsp(j)%p%fma%fx+dm%mb(i)%dsp(j)%p%fma%fx
					all%dsp(j)%p%fna%fx=all%dsp(j)%p%fna%fx+dm%mb(i)%dsp(j)%p%fna%fx
				end do
			end do
			do j=1, n_tot_comp
				all%all%fden%fx=all%all%fden%fx+all%dsp(j)%p%fden%fx
				all%all%fma%fx=all%all%fma%fx+all%dsp(j)%p%fma%fx
				all%all%fna%fx=all%all%fna%fx+all%dsp(j)%p%fna%fx
			end do
		end associate
    end subroutine
    
    subroutine get_nxj(dm)
        implicit none
        class(diffuse_mspec)::dm
        integer i, j, k,l
        do i=1, dm%n
            associate(mb=>dm%mb(i))
				do j=1, n_tot_comp
					call dms_so_get_nxj_from_nejw(mb%dsp(j)%p,dm%jbin_type)
				end do
                
                !if(mb%star%n.ne.0.or.mb%sbh%n.ne.0.or.mb%bbh%n.ne.0.or.mb%bstar%n.ne.0)then
                    !print*, "mc=", mb%mc
                    !if(mb%star%n.ne.0) print*, "num_star=", mb%star%n, mb%star%n_real
                    !if(mb%sbh%n.ne.0) print*, "num_sbh=", mb%sbh%n, mb%sbh%n_real
                    !if(mb%bbh%n.ne.0) print*, "num_bbh=", mb%bbh%n, mb%bbh%n_real
                    !if(mb%bstar%n.ne.0) print*, "num_bstar=", mb%bstar%n, mb%bstar%n_real
                !end if
                mb%all%nxj%nxyw=0
                do j=1, mb%nbin_gx
                    do k=1, mb%nbin_gx
						do l=1, n_tot_comp
							mb%all%nxj%nxyw(j,k)=mb%all%nxj%nxyw(j,k)+mb%dsp(l)%p%nxj%nxyw(j,k)
						end do
                    end do 
                end do
            end associate
        end do
       ! print*,"get_nxj:"
       ! print*, dm%mb(1)%star%nxj%nxyw(:,3)
        
        !associate(star0=>dm%all%star%nxj, bstar0=>dm%all%bstar%nxj, &
        !    sbh0=>dm%all%sbh%nxj, bbh0=>dm%all%bbh%nxj, all0=>dm%all%all%nxj, &
		!	wd0=>dm%all%wd%nxj, ns0=>dm%all%ns%nxj, bd0=>dm%all%bd%nxj)
		do j=1, n_tot_comp
			dm%all%dsp(j)%p%nxj%nxyw=0
			dm%all%dsp(j)%p%n=0
			dm%all%dsp(j)%p%n_real=0
		end do
		dm%all%all%nxj%nxyw=0
		
        do i=1, dm%n
			dm%mb(i)%all%n_real=0
			dm%mb(i)%all%n=0
			do j=1, n_tot_comp
				dm%all%dsp(j)%p%nxj%nxyw=dm%all%dsp(j)%p%nxj%nxyw+dm%mb(i)%dsp(j)%p%nxj%nxyw
				dm%mb(i)%all%n_real=dm%mb(i)%all%n_real+dm%mb(i)%dsp(j)%p%n_real
				dm%mb(i)%all%n=dm%mb(i)%all%n+dm%mb(i)%dsp(j)%p%n
				dm%all%dsp(j)%p%n=dm%all%dsp(j)%p%n+dm%mb(i)%dsp(j)%p%n
				dm%all%dsp(j)%p%n_real=dm%all%dsp(j)%p%n_real+dm%mb(i)%dsp(j)%p%n_real
			end do
			dm%all%all%nxj%nxyw=dm%all%all%nxj%nxyw+dm%mb(i)%all%nxj%nxyw
		end do    
       
    end subroutine
!    subroutine get_nejw(dm, estar, jstar, wstar, mstar, n_star)
!        implicit none
!        class(diffuse_mspec)::dm
!        integer n_star
!        real(8):: estar(n_star),jstar(n_star),wstar(n_star),mstar(n_star)
!        real(8) m1, m2
!        integer i
!
!        do i=1, dm%n
!            associate(mb=>dm%mb(i))
!                m1=mb%m1; m2=mb%m2
!                call get_ejw_from_particle(estar,jstar, wstar,mstar,n_star, &
!                    m1, m2, dm%mbh, dm%v0, mb%nejwstar,mb%n_star)
!             !   call get_ejw_from_particle(ps_bstar, m1, m2, dm%mbh, dm%v0, mb%nejwbstar,mb%n_bstar)
!             !   call get_ejw_from_particle(ps_bbh, m1, m2, dm%mbh, dm%v0, mb%nejwbbh,mb%n_bbh)
!             !   call get_ejw_from_particle(ps_sbh, m1, m2, dm%mbh, dm%v0, mb%nejwsbh,mb%n_sbh)
!              !  mb%n_tot=mb%n_star+mb%n_bstar+mb%n_bbh+mb%n_sbh
!              !  print*, "mb%n_tot=", mb%n_tot, m1, m2
!            end associate            
!        end do
!    end subroutine
    
    subroutine get_ejw_from_particle(estar, jstar, wstar, mstar, n_star,&
                 m1,  m2, mbh, v0, xb, nejw, nsam)
        implicit none
        integer n_star
        real(8):: estar(n_star),jstar(n_star),wstar(n_star),mstar(n_star), xstar(n_star)
        integer i, nsam, idx(n_star)
        type(nejw_type),allocatable::nejw(:)
        real(8) jmax, m1,m2, v0, mbh,xb
		!integer jb_type

        xstar=abs(estar(1:n_star))/v0**2
        !call get_n_from_particlex(mstar,xstar, n_star, m1, m2, xb, idx, nsam)
        call get_n_from_particle(mstar, n_star, m1, m2, idx, nsam)
        if(allocated(nejw)) deallocate(nejw)
        allocate(nejw(nsam))

        do i=1, nsam
            nejw(i)%e=log10(xstar(idx(i)))
            if(isnan(nejw(i)%e))then
                print*, "nejw=NaN:", estar(idx(i)), v0
                read(*,*)
            end if

			!select case(jb_type)
			!case(Jbin_type_lin)
				nejw(i)%j=jstar(idx(i))!/jmax
			!case(jbin_type_log)
			!	nejw(i)%j=log10(jstar(idx(i)))!/jmax
			!case(jbin_type_sqr)
			!	nejw(i)%j=jstar(idx(i))**2
			!end select
            
            nejw(i)%w=wstar(idx(i))
            nejw(i)%idx=idx(i)
        end do
!        nsam=0
!        do i=1, n_star
!            if(mstar(i).ge.m1.and.mstar(i).le.m2)then
!                nsam=nsam+1
!                nejw(nsam)%e=log10(abs(estar(i))/v0**2)
!                if(isnan(nejw(nsam)%e))then
!                    print*, "nejw=NaN:", estar(i), v0
!                    read(*,*)
!                end if
!            ! print*, abs(sp%en),v0**2
!            ! read(*,*)
!                !jmax=mbh/(2*sp%en)**0.5d0
!                nejw(nsam)%j=jstar(i)!/jmax
!                nejw(nsam)%w=wstar(i)
!                nejw(nsam)%idx=i
!            end if
!        end do
    end subroutine    


    subroutine output_diffuse_mspec_bin(dm, fl)
        implicit none
        class(diffuse_mspec)::dm
        character*(*) fl
        integer i 
        open(unit=2999,file=trim(adjustl(fl)), access='stream', form='unformatted')
        write(unit=2999) dm%n, dm%idx_ref, dm%nbin_grid, dm%nbin_gx, dm%emin, &
        dm%emax, dm%jmin,dm%jmax, dm%mbh, dm%v0, dm%n0, dm%rh, &
        dm%acmin,dm%acmax,dm%x_boundary,dm%jbin_type,dm%grid_type
    
        do i=1, dm%n
            call dm%mb(i)%write_mb(2999)
        end do
        call dm%dc0%write_grid(2999)
        !call dm%dc0_bk%write_grid(2999)
        write(2999) dm%weight_asym
		write(2999) dm%all%all%fna
		!print*, dm%nbin_grid, dm%nbin_gx, dm%emin, dm%emax, dm%jmin, dm%jmax, &
        !dm%mbh, dm%v0, dm%n0, dm%rh, dm%acmin, dm%acmax
		!print*, dm%all%all%fna%nbin
        !write(2999) dm%asymp_star, dm%asymp_sbh, dm%asymp_bstar, dm%asymp_bbh
       ! write(2999) dm%asymp_star_norm_factor, dm%asymp_sbh_norm_factor, &
        !    dm%asymp_bstar_norm_factor, dm%asymp_bbh_norm_factor    
        close(2999)

    end subroutine
    subroutine input_diffuse_mspec_bin(dm, fl)
        implicit none
        class(diffuse_mspec)::dm
        character*(*) fl
        integer i, n
        open(unit=2999,file=trim(adjustl(fl)), access='stream', form='unformatted', status='old')
        read(2999) n , dm%idx_ref, dm%nbin_grid, dm%nbin_gx, dm%emin, &
        dm%emax, dm%jmin,dm%jmax, dm%mbh, dm%v0, dm%n0, dm%rh, &
        dm%acmin,dm%acmax,dm%x_boundary,dm%jbin_type,dm%grid_type

        !call set_diffuse_mspec(dm, dm%nbin_grid, dm%nbin_gx, dm%emin, dm%emax, &
        !dm%jmin, dm%jmax, dm%mbh, dm%v0, dm%n0, dm%rh, dm%x_boundary, dm%idx_ref,dm%jbin_type)
		call dm%init(n)
        do i=1, dm%n
            call dm%mb(i)%read_mb(2999)
        end do
       ! print*, "read 1 finished"
        call dm%dc0%read_grid(2999)
        !call dm%dc0_bk%read_grid(999)
        !print*, "read 2 finished"
        read(2999) dm%weight_asym!, dm%idx_ref
		!print*, dm%weight_asym, dm%idx_ref
		!print*, "dm%all%all%fna%nbin=", dm%all%all%fna%nbin
        !read(2999) dm%nbin_grid, dm%nbin_gx, dm%emin, dm%emax, dm%jmin, dm%jmax, &
        !dm%mbh, dm%v0, dm%n0, dm%rh, dm%acmin, dm%acmax
		read(2999) dm%all%all%fna
		!print*, dm%nbin_grid, dm%nbin_gx, dm%emin, dm%emax, dm%jmin, dm%jmax, &
        !dm%mbh, dm%v0, dm%n0, dm%rh, dm%acmin, dm%acmax
		!print*, dm%all%all%fna%nbin
		!stop
        !read(2999) dm%asymp_star, dm%asymp_sbh, dm%asymp_bstar, dm%asymp_bbh
        !read(2999) dm%asymp_star_norm_factor, dm%asymp_sbh_norm_factor, &
        !dm%asymp_bstar_norm_factor, dm%asymp_bbh_norm_factor
        !call read_s2d(dm%fvm, 999)   
        close(2999)
       
    end subroutine
end module

subroutine get_n_from_particlex(mstar,xstar, n_star,&
    m1,  m2, x_b, idx, nsam)
    implicit none
    integer n_star
    real(8):: mstar(n_star), xstar(n_star), x_b
    integer i, nsam, idx(n_star)
    real(8) jmax, m1,m2
    nsam=0
    if(n_star.eq.0) return
    do i=1, n_star
        !print*, ps_arr%sp(i)%m, m1, m2
        if(mstar(i).ge.m1.and.mstar(i).le.m2 &
            .and.xstar(i).ge.x_b)then
            nsam=nsam+1
            idx(nsam)=i
        end if
    end do
    !print*, "nsam, n_star=",nsam, n_star
end subroutine

subroutine get_n_from_particle(mstar,n_star,&
    m1,  m2, idx, nsam)
    implicit none
    integer n_star
    real(8):: mstar(n_star)
    integer i, nsam, idx(n_star)
    real(8) jmax, m1,m2
    nsam=0
    if(n_star.eq.0) return
    do i=1, n_star
        !print*, ps_arr%sp(i)%m, m1, m2
        if(mstar(i).ge.m1.and.mstar(i).le.m2)then
            nsam=nsam+1
            idx(nsam)=i
        end if
    end do
    !print*, "nsam, n_star=",nsam, n_star
end subroutine