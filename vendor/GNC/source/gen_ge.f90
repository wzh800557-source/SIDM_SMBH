
subroutine get_ge_by_root(smsa, n,norm_in)
    use com_main_gw
    integer n
    type(particle_samples_arr_type)::smsa(n), smstot
    logical::norm_in

    if(rid.eq.mpi_master_id)then
        call smmerge_arr_single(smsa,n,smstot)
        print*, "smstot%n=", smstot%n
        call get_ge(dms, smstot,norm_in)
        !print*, "finished get ge"
        deallocate(smstot%sp)
    end if
end subroutine
subroutine get_ge(dm, sms_arr_single,norm_in)
    use com_main_gw
    implicit none
    type(diffuse_mspec)::dm
    type(particle_samples_arr_type)::sms_arr_single
    type(particle_samples_arr_type)::sms_arr_star, sms_arr_sbh
	type(particle_samples_arr_type)::sms_arr_wd, sms_arr_ns,sms_arr_bd
    logical::norm_in

    !print*, "get ge"
   ! print*, "sms_arr_by%n=", sms_arr_by%n
    call sams_get_weight_clone_single(sms_arr_single)
    if(norm_in)then
        dm%weight_asym=1d0
        !call set_asym_weight_arr(sms_arr_single, sms_arr_by)
    end if
    call set_real_weight_arr_single(sms_arr_single)

	call separate_to_species(sms_arr_single,  sms_arr_star, &
        sms_arr_sbh, sms_arr_wd, sms_arr_ns,sms_arr_bd)

	!print*, "start get_gx", sms_arr_bbh%n, rid
    !print*, sms_arr_bbh%n
	!print*, sms_arr_bd%n
    call gen_gx(dm, sms_arr_star,  sms_arr_sbh, sms_arr_wd, &
		sms_arr_ns,sms_arr_bd)
    if(norm_in)then
        !print*, "xbounday=",dms%x_boundary
        !call dm%mb(1)%sbh%barge%print("sbh")
        call dm%get_asymp_norm_factor()
        !print*, dm%weight_asym
        
        !print*, "gen_gx:finished!"
    ! dm%asymp_star_norm_factor=1d0
        !call dms%mb(1)%star%barge%print("bf")
        call dm%normalize_barge()
        !call dm%mb(5)%sbh%barge%print("sbh")
        !stop
        !call dms%mb(1)%star%barge%print("af")
        !call set_asym_weight_arr(sms_arr_single, sms_arr_by)

    end if
    print*, "start get_dens0"
    print*, "weight_asym=", dm%weight_asym
    call dm%get_dens0()
    print*, "finished get_dens0"
    call set_real_weight_arr_single(sms_arr_single)
    if(norm_in)then
        call separate_to_species(sms_arr_single,  sms_arr_star, &
         sms_arr_sbh, sms_arr_wd, sms_arr_ns,sms_arr_bd)
        call conv_dms_nejw(dm, sms_arr_star, &
        sms_arr_sbh, sms_arr_wd, sms_arr_ns, sms_arr_bd)
    end if    

    if(allocated(sms_arr_star%sp))then
        deallocate(sms_arr_star%sp)
    end if
    if(allocated(sms_arr_sbh%sp))then
        deallocate(sms_arr_sbh%sp)
    end if
	if(allocated(sms_arr_wd%sp))then
        deallocate(sms_arr_wd%sp)
    end if
	if(allocated(sms_arr_ns%sp))then
        deallocate(sms_arr_ns%sp)
    end if
	if(allocated(sms_arr_bd%sp))then
        deallocate(sms_arr_bd%sp)
    end if
    print*, "finished get_ge"

end subroutine
subroutine 	get_fma(fna, fma, mass)
    use com_sts_type
    implicit NONE
    type(s1d_type)::fna, fma
    real(8) mass
    integer i
    
    do i=1, fna%nbin
        fma%fx(i)=fna%fx(i)*mass
    end do	
end subroutine   

subroutine dm_get_dc_mpi(dm)
    use com_main_gw
    implicit none
    type(diffuse_mspec):: dm
    interface
        subroutine mb_get_dc_mpi(mb)
            use com_main_gw
            implicit none
            type(mass_bins),target::mb
        end subroutine
    end interface
    integer i
    do i=1, dm%n
     !   print*, "i, m=", i, dm%mb(i)%mc
        call mb_get_dc_mpi(dm%mb(i))
    end do
    !call get_dc0_bk(dm)
    call get_dc0(dm)
    !print*, "finished NR grid",rid
    !call dm%all%all%fna%print("all_fna")
    !call dm%favg_mass%print("favg_mass")
    !call dm%all%all%fma%print("all_fma")
    !call dm%fvm%print("fvm")
    !print*, "====================="
    !call dm%dc0%s2_dRRJ%print("s2_drrj")
    !call dm%dc0%s2_dRRJJ%print("s2_drrjj")
    !read(*,*)
end subroutine

subroutine mb_get_dc_mpi(mb)
    use com_main_gw
    use sidm_params, only: use_sidm, sidm_kappa0
    implicit none
    type(mass_bins),target::mb
    real(8),allocatable:: s0(:,:), s110(:,:),s111(:,:)
    real(8),allocatable:: s131(:,:), s130(:,:), s13_1(:,:)
    real(8),allocatable:: s330(:,:), s310(:,:),ycenter(:)
    integer n,i,j
    type(coeff_type)::cej
    real(8) kappa, sigma32, n0
    interface
        subroutine mb_get_dc_mpi_sigma(mc,nbin, xb, yb, &
            s0, s110, s111,s131,s130,s13_1,s330,s310, barge, asymp)
            use com_main_gw
            implicit none
            integer nbin
            real(8) mc, s0(nbin,nbin), s110(nbin,nbin),s111(nbin,nbin), s131(nbin,nbin)
            real(8) s130(nbin,nbin), s13_1(nbin,nbin), s330(nbin,nbin), s310(nbin,nbin)
            real(8) xb(nbin),yb(nbin),  asymp
            type(s1d_type),target::barge
        end subroutine
    end interface

    n=mb%nbin_grid
    allocate(ycenter(n))
    allocate(s0(n, n), s110(n,n), s111(n,n),s131(n,n), s130(n,n), s13_1(n,n), &
            s330(n,n), s310(n,n))
    select case(ctl%jbin_type)
    case(Jbin_type_lin)
        ycenter(1:n)=mb%dc%s2_de_110%ycenter(1:n)
    case(jbin_type_log)
        ycenter(1:n)=10**mb%dc%s2_de_110%ycenter(1:n)
	case(jbin_type_sqr)   
		ycenter(1:n)=mb%dc%s2_de_110%ycenter(1:n)**0.5d0
	case default
		print*, "error! define jbin_type_lin", ctl%jbin_type
		stop
    end select
	select case (ctl%model_intej)
	case(model_intej_fast)
        !print*, "mb%all%asymp=",mb%all%asymp
		call mb_get_dc_mpi_sigma(mb%mc,mb%nbin_grid, 10**mb%dc%s2_de_110%xcenter, &
			ycenter,s0, s110, s111,s131,s130,s13_1,s330,s310, mb%all%barge, mb%all%asymp)
	case default
		print*, "model_intej:error! define:",ctl%model_intej
		stop
	end select
    if (use_sidm) then
        kappa = sidm_kappa0()
    else
        kappa = (4*pi*mb%mc)**2*log(mb%mbh/mb%mc)
    end if

    !====================test=========================
    !kappa=(4*pi*mb%mc)**2*log(mb%mbh/mb%mc)!*factor
    !=================================================

    sigma32=(2*pi*mb%v0**2)**(-3/2d0)
    n0=mb%n0

    do i=1, n
        do j=1, n
            call get_coeff_ej(cej,ycenter(j), &
            s0(i,j), s110(i,j), &
            s111(i,j),s131(i,j),s130(i,j),s13_1(i,j),s330(i,j),s310(i,j))

            mb%dc%s2_de_110%fxy(i,j)=cej%e_110*sigma32*n0*kappa
            mb%dc%s2_de_0%fxy(i,j)=cej%e_0 *sigma32*n0*kappa
            if(cej%ee<0d0) cej%ee=abs(cej%ee)
            mb%dc%s2_dee%fxy(i,j)=cej%ee*sigma32*n0*kappa
            mb%dc%s2_dj_111%fxy(i,j)=cej%j_111 *sigma32*n0*kappa
            mb%dc%s2_dj_rest%fxy(i,j)=cej%j_rest *sigma32*n0*kappa
            if(cej%jj<0d0) cej%jj=abs(cej%jj)
            mb%dc%s2_djj%fxy(i,j)=cej%jj*sigma32*n0*kappa
            mb%dc%s2_dej%fxy(i,j)=cej%ej*sigma32*n0*kappa
        end do
    end do

end subroutine
subroutine get_steps_grid(dm,mi, dt)
	use com_main_gw
	implicit none
	type(diffuse_mspec)::dm
	type(s2d_type)::dt
	integer i, j,k, mi
	call dt%init(dm%nbin_grid,dm%nbin_grid,dm%emin,dm%emax,dm%jmin, dm%jmax,sts_type_dstr)
	do i=1, dm%n
		do j=1, dm%nbin_grid
			do k=1,dm%nbin_grid
				
			end do
		end do
	end do
end subroutine
subroutine set_dm_init(dm)
    use com_main_gw
    implicit none
    type(diffuse_mspec)::dm


	call dm%set_diffuse_mspec(ctl%grid_bins,ctl%gx_bins, &
        log10emin_factor, log10emax_factor, &
        jmin_value ,jmax_value, mbh,  ctl%v0, ctl%n0, rh,ctl%x_boundary, &
        ctl%idx_ref, ctl%jbin_type,ctl%grid_type)

    call dm%init(ctl%m_bins)
    !select case(ctl%bin_mass_mode) 
		
	call set_mass_bin_mass_given(dm, ctl%bin_mass, ctl%bin_mass_m1,&
		ctl%bin_mass_m2, ctl%asymptot, ctl%m_bins)
end subroutine


subroutine mb_get_dc_mpi_sigma(mc,nbin, xb, yb, &
        s0, s110, s111,s131,s130,s13_1,s330,s310, barge, asymp)
    use com_main_gw
    implicit none
    !type(mass_bins),target::mb
    integer nbin
    real(8) mc, s0(nbin,nbin), s110(nbin,nbin),s111(nbin,nbin), s131(nbin,nbin), s130(nbin,nbin), s13_1(nbin,nbin), &
        s330(nbin,nbin), s310(nbin,nbin)
    real(8) xb(nbin),yb(nbin)
    integer i, j, ierr
	integer nbg, ned, nblock,ntasks
    real(8) kappa, factor, energyx, jum, asymp
    type(s1d_type),target::barge
    type(coeff_type)::cej
    external::fgx_mb, fx_g, fgx_mb_star

    if(.not.allocated(ctl%cc))allocate(ctl%cc(1))
    cct_share=>ctl%cc(1)
    ctl%cc(1)%alpha=7d0/4d0
    
    fc_share=>barge!_norm
    fgx_g0=asymp
    if(rid.eq.0) print*, "mc, fgx_g0, bin=",mc, fgx_g0, nbin
   ! print*, "mc, fgx_g0=",mc, fgx_g0, nbin, asymp, rid
    !if(fgx_g0.eq.0)then
    !    call fc_share%print()
    !    read(*,*)
    !end if
    if(all(fc_share%fx.eq.0d0))then
        !call fc_share%print()
        !stop
        s0=0
        s110=0
        s111=0
        s131=0
        s130=0
        s13_1=0
        s330=0
        s310=0
        return
    end if
	nbg=ctl%nblock_mpi_bg
	ned=ctl%nblock_mpi_ed
	nblock=ctl%nblock_size
	ntasks=ctl%ntasks
    do i=nbg, ned
        energyx=xb(i)
        
        call get_sigma0(energyx, fgx_mb, s0(i,1))
        s0(i,2:nbin)=s0(i,1)
        !print*, "i,sigma0=",i,s0(i,1)
        do j=1, nbin
            jum=yb(j)
            !call get_coeff_sigma_funcs_cfs_rk(energyx,jum, fgx_mb, s110(i,j), &
            !    s111(i,j),s131(i,j),s130(i,j),s13_1(i,j),s330(i,j),s310(i,j))
            
            call get_coeff_sigma_funcs_cfs_grid(energyx,jum, fgx_mb, s110(i,j), &
                s111(i,j),s131(i,j),s130(i,j),s13_1(i,j),s330(i,j),s310(i,j))                
           
        end do
    end do
    call mpi_BARRIER(mpi_comm_world, ierr)
    call collect_data_mpi(s0, nbin,nbg, ned, nblock, ntasks)
    call mpi_BARRIER(mpi_comm_world, ierr)
    call collect_data_mpi(s110, nbin,nbg, ned, nblock, ntasks)
    call mpi_BARRIER(mpi_comm_world, ierr)
    call collect_data_mpi(s111, nbin,nbg, ned, nblock, ntasks)
    call mpi_BARRIER(mpi_comm_world, ierr)
    call collect_data_mpi(s131, nbin,nbg, ned, nblock, ntasks)
    call mpi_BARRIER(mpi_comm_world, ierr)
    call collect_data_mpi(s130, nbin,nbg, ned, nblock, ntasks)    
    call mpi_BARRIER(mpi_comm_world, ierr)
    call collect_data_mpi(s13_1, nbin,nbg, ned, nblock, ntasks)
    call mpi_BARRIER(mpi_comm_world, ierr)
    call collect_data_mpi(s330, nbin,nbg, ned, nblock, ntasks)
    call mpi_BARRIER(mpi_comm_world, ierr)
    call collect_data_mpi(s310, nbin,nbg, ned, nblock, ntasks)

end subroutine


subroutine conv_dms_nejw_obj(dm,en,jm,m,w_real,nobj, obj_type)
    use com_main_gw
    implicit none
    type(diffuse_mspec)::dm
    integer nobj
    real(8) en(nobj),jm(nobj), w_real(nobj), wobj(nobj), m(nobj)
    integer i, obj_type
    real(8) m1, m2
    integer,parameter::obj_type_star=1, obj_type_sbh=2, obj_type_bstar=3
	integer,parameter::obj_type_bbh=4, obj_type_wd=5, obj_type_ns=6, obj_type_bd=7

    
    if(nobj>0)then
        wobj(1:nobj)=w_real(1:nobj)/dble(ctl%ntasks)
	else
		return
    end if

    do i=1, dm%n
        associate(mb=>dm%mb(i))
            m1=mb%m1; m2=mb%m2
            select case(obj_type)
            case(obj_type_star)
                call get_ejw_from_particle(en(1:nobj),&
                jm(1:nobj),wobj(1:nobj),m(1:nobj), &
                nobj, m1, m2, dms%mbh, dms%v0, dms%x_boundary, mb%star%nejw,mb%star%n)
            case(obj_type_sbh)
                !call get_ejw_from_particle(bkobj%sp(1:nobj)%en,&
                !bkobj%sp(1:nobj)%jm,wobj(1:nobj),bkobj%sp(1:nobj)%m, &
                !nobj, m1, m2, dms%mbh, dms%v0, dms%x_boundary, mb%sbh%nejw,mb%sbh%n)
                call get_ejw_from_particle(en(1:nobj),&
                jm(1:nobj),wobj(1:nobj),m(1:nobj), &
                nobj, m1, m2, dms%mbh, dms%v0, dms%x_boundary, mb%sbh%nejw,mb%sbh%n)
            !case(obj_type_bstar)
            !    call get_ejw_from_particle(en(1:nobj),&
            !    jm(1:nobj),wobj(1:nobj),m(1:nobj), &
            !    nobj, m1, m2, dms%mbh, dms%v0, dms%x_boundary, mb%bstar%nejw,mb%bstar%n)             
            !case(obj_type_bbh)
            !    call get_ejw_from_particle(en(1:nobj),&
            !    jm(1:nobj),wobj(1:nobj),m(1:nobj), &
            !    nobj, m1, m2, dms%mbh, dms%v0, dms%x_boundary, mb%bbh%nejw,mb%bbh%n)             
			case(obj_type_wd)
				call get_ejw_from_particle(en(1:nobj),&
                jm(1:nobj),wobj(1:nobj),m(1:nobj), &
                nobj, m1, m2, dms%mbh, dms%v0, dms%x_boundary, mb%wd%nejw,mb%wd%n)
			case(obj_type_ns)
				call get_ejw_from_particle(en(1:nobj),&
                jm(1:nobj),wobj(1:nobj),m(1:nobj), &
                nobj, m1, m2, dms%mbh, dms%v0, dms%x_boundary, mb%ns%nejw,mb%ns%n)
			case(obj_type_bd)
				call get_ejw_from_particle(en(1:nobj),&
                jm(1:nobj),wobj(1:nobj),m(1:nobj), &
                nobj, m1, m2, dms%mbh, dms%v0, dms%x_boundary, mb%bd%nejw,mb%bd%n)
								
			case default
				print*, "error! define obj_type"
				stop
            end select
        end associate      
    end do
    
end subroutine

subroutine conv_dms_nejw(dm, bkstar, bksbh, bkwd, bkns,bkbd)
    use com_main_gw
    implicit none
    type(diffuse_mspec)::dm
    type(particle_samples_arr_type)::bkstar, bksbh, bkwd, bkns,bkbd
    integer i
    integer,parameter::obj_type_star=1, obj_type_sbh=2, obj_type_bstar=3
	integer,parameter::obj_type_bbh=4, obj_type_wd=5, obj_type_ns=6, obj_type_bd=7
    
    call conv_dms_newj_obj_one(dm, bkstar,obj_type_star)
    call conv_dms_newj_obj_one(dm, bksbh,obj_type_sbh)
    call conv_dms_newj_obj_one(dm, bkwd,obj_type_wd)
    call conv_dms_newj_obj_one(dm, bkns,obj_type_ns)
    call conv_dms_newj_obj_one(dm, bkbd,obj_type_bd)

    do i=1, dm%n
        associate(mb=>dm%mb(i))
            mb%all%n=mb%star%n+mb%sbh%n+mb%wd%n+mb%ns%n+&
				mb%bd%n            
        end associate            
    end do
end subroutine
subroutine conv_dms_newj_obj_one(dm, bk, i_obj)
    use com_main_gw
    implicit none
    type(particle_samples_arr_type)::bk
    type(diffuse_mspec)::dm
    integer i_obj
    real(8) en(bk%n), jm(bk%n), mass(bk%n), weight_real(bk%n)
    en=bk%sp(1:bk%n)%en
    jm=bk%sp(1:bk%n)%jm
    mass=bk%sp(1:bk%n)%m
    weight_real=bk%sp(1:bk%n)%weight_real

    call conv_dms_nejw_obj(dm,en, jm , mass, weight_real, bk%n, i_obj)
end subroutine

subroutine separate_to_species(bks, bkstar, bksbh, bkwd, bkns, bkbd)
    use com_main_gw
    implicit none
    type(particle_samples_arr_type)::bks, bkstar, bksbh, bkwd, bkns,bkbd
    
   call separate_to_species_single(bks, bkstar, bksbh, bkwd, bkns,bkbd)
end subroutine
subroutine separate_to_species_single(bks, bkstar, bksbh, bkwd, bkns, bkbd)
    use com_main_gw
    implicit none
    type(particle_samples_arr_type)::bks, bkstar, bksbh, bkwd, bkns, bkbd

    call sams_arr_select_type_single(bks, bkstar, star_type_MS)
    call sams_arr_select_type_single(bks, bksbh, star_type_BH)
	call sams_arr_select_type_single(bks, bkwd, star_type_WD)
	call sams_arr_select_type_single(bks, bkns, star_type_NS)
	call sams_arr_select_type_single(bks, bkbd, star_type_BD)

end subroutine

subroutine gen_gx(dm, bkstar,  bksbh, bkwd, bkns,bkbd)
    use com_main_gw
    use flux_diag
    implicit none
    integer, save :: iflux_snap = 0
    type(diffuse_mspec)::dm
    type(particle_samples_arr_type)::bkstar, bksbh, bkwd, bkns
	type(particle_samples_arr_type)::bkbd

    !print*, "gen_gx:bkbbh%n=", bkbbh%n
    call conv_dms_nejw(dm, bkstar, bksbh, bkwd,bkns,bkbd)
    !call print_num_boundary(dm)
    call dm%get_nxj()
   ! print*, "1"
   ! print*, "gen_gx:dm%all%bbh%n=", dm%all%bbh%n, rid
    call dm%get_fxj0()
   ! print*, "2"
    !print*, "gen_gx:3, finished!"
    call dm%get_barge0()
        call output_boundary_fluxes(dm, iflux_snap)
        iflux_snap = iflux_snap + 1
	!print*, "dm%get_barge0"
	!call dm%mb(1)%star%barge%print("dm%mb(1)%star%barge")
  !  print*, "3"
    !print*, "gen_gx:4, finished!"

end subroutine
subroutine print_num_all(dm)
    use com_main_gw
    implicit none
    type(diffuse_mspec)::dm
    integer i
    real(8) nb(ctl%m_bins,10),nbw(ctl%m_bins,10)
    
    nb=0
    write(*,fmt=*)  "print_num_all=========================="
    write(*,fmt="(A3, 10A13)") "i","mc", "star", "sbh", "bbh", "ns", "wd","bd"
    do i=1, dm%n
        call get_num_all(nb(i,1),  nbw(i,1),dm%mb(i)%star)
        call get_num_all(nb(i,2),  nbw(i,2),  dm%mb(i)%sbh)
        call get_num_all(nb(i,3),  nbw(i,3), dm%mb(i)%bbh)
		call get_num_all(nb(i,4),  nbw(i,4), dm%mb(i)%ns)
		call get_num_all(nb(i,5),  nbw(i,5), dm%mb(i)%wd)
		call get_num_all(nb(i,6),  nbw(i,6), dm%mb(i)%bd)
	end do	
	do i=1, dm%n
    	write(*,fmt="(I3,7F13.2)") i, dm%mb(i)%mc, nb(i,1:6)
	end do
	write(*,fmt="(A3, 8A13)") "i","mc", "starw", "sbhw", "bbhw", "nsw", "wdw","bdw"
	do i=1, dm%n
        write(*,fmt="(I3,7F13.2)") i, dm%mb(i)%mc, nbw(i,1:6)
    end do
	!write(*,fmt="(A3, 8A13)") "i","mc","j>0.5","","","j<0.5"
    !do i=1, dm%n
    !    write(*,fmt="(I3,7F13.2)") i, dm%mb(i)%mc, nbj1(i,1:3), nbj2(i,1:3)
    !end do
    write(*,fmt=*) "end of print_num_all==================="
end subroutine
subroutine check_boundary(str)
    use com_main_gw
    implicit none
    type(chain_pointer_type),pointer::pt
    integer n
    character*(*) str
    
    n=0
    pt=>bksams%head
    do while(associated(pt))
        if(pt%ob%exit_flag.eq.exit_normal)then
            if(pt%ob%en>ctl%energy_boundary)then
                n=n+1
            end if
        end if
        pt=>pt%next
    end do
    if(n.ne.200)then
        print*, str, " error, n=", n
        stop
    end if
end subroutine
subroutine print_num_boundary(dm)
    use com_main_gw
    implicit none
    type(diffuse_mspec)::dm
    integer i
    real(8) nb(ctl%m_bins,10),nbw(ctl%m_bins,10),nbj1(ctl%m_bins,10),nbj2(ctl%m_bins,10)
    
    nb=0
    write(*,fmt=*)  "print_num_boundary=========================="
    write(*,fmt="(A3, 10A13)") "i","mc", "star", "sbh", "bbh", "ns", "wd","bd"
    do i=1, dm%n
        call get_num_boundary(nb(i,1), nbj1(i,1), nbj2(i,1), nbw(i,1),dm%mb(i)%star)
        call get_num_boundary(nb(i,2), nbj1(i,2), nbj2(i,2), nbw(i,2),  dm%mb(i)%sbh)
        call get_num_boundary(nb(i,3), nbj1(i,3), nbj2(i,3), nbw(i,3), dm%mb(i)%bbh)
		call get_num_boundary(nb(i,4), nbj1(i,4), nbj2(i,4), nbw(i,4), dm%mb(i)%ns)
		call get_num_boundary(nb(i,5), nbj1(i,5), nbj2(i,5), nbw(i,5), dm%mb(i)%wd)
		call get_num_boundary(nb(i,6), nbj1(i,6), nbj2(i,6), nbw(i,6), dm%mb(i)%bd)
	end do	
	do i=1, dm%n
    	write(*,fmt="(I3,7F13.2)") i, dm%mb(i)%mc, nb(i,1:6)
	end do
	write(*,fmt="(A3, 8A13)") "i","mc", "starw", "sbhw", "bbhw", "nsw", "wdw","bdw"
	do i=1, dm%n
        write(*,fmt="(I3,7F13.2)") i, dm%mb(i)%mc, nbw(i,1:6)
    end do
	!write(*,fmt="(A3, 8A13)") "i","mc","j>0.5","","","j<0.5"
    !do i=1, dm%n
    !    write(*,fmt="(I3,7F13.2)") i, dm%mb(i)%mc, nbj1(i,1:3), nbj2(i,1:3)
    !end do
    write(*,fmt=*) "end of print_num_boundary==================="
end subroutine
subroutine get_num_all(nb,nbw, so)
    use com_main_gw
    implicit none
    real(8) nb, nbw
    type(dms_stellar_object)::so
    integer i
    nb=0;nbw=0; 
    do i=1, so%n
        if(so%nejw(i)%e>log10(ctl%x_boundary).and.so%nejw(i)%e<log10(emax_factor))then
            nbw=nbw+so%nejw(i)%w
            nb=nb+1
					
        end if
    end do
end subroutine
subroutine get_num_boundary(nb,nbj1,nbj2,nbw, so)
    use com_main_gw
    implicit none
    real(8) nb, nbw, nbj1,nbj2
    type(dms_stellar_object)::so
    integer i
    nb=0;nbw=0;nbj1=0; nbj2=0
    do i=1, so%n
        if(so%nejw(i)%e<log10(ctl%x_boundary))then
            nbw=nbw+so%nejw(i)%w
            nb=nb+1
			if(so%nejw(i)%j>0.5d0)then
				nbj1=nbj1+1
			else
				nbj2=nbj2+1
			end if	
        end if
		
        if(so%nejw(i)%e<log10emin_factor)then
            print*, "error!??:", so%nejw(i)%e, log10emin_factor
        end if
    end do
end subroutine

!subroutine get_fden_simu(dm, bkstar, bkbystar, bksbh, bkbbh)
!    use com_main_gw
!    implicit none
!    type(diffuse_mspec)::dm
!    type(particle_samples_arr_type)::bkstar, bksbh
!    type(samples_arr_type)::bkbystar, bkbbh
!    
!!    call gen_gx(dm, bkstar, bkbystar, bksbh, bkbbh)
!    
!    call get_fden_sample_particle(bkstar, dm%all%star%fden_simu)
!    call get_fden_sample_particle(bksbh, dm%all%sbh%fden_simu)
!    !call dm%all%star%fden_simu%print("star:fden_simu")
!   ! call dm%all%star%fden%print("star:fden")
!   ! print*, "dms%n0=",dms%n0
!   ! print*, dm%all%star%fden_simu%fx/dm%all%star%fden%fx
!   ! stop
!end subroutine
subroutine get_fden_sample_particle(en,jum,w,n,fden)
	use com_main_gw
	implicit none
	type(s1d_type)::fden
	integer i,n,j
	real(8) a(n), e(n), en(n),jum(n),w(n), ivr(n)
    real(8) ivrsum

	a(1:n)=mbh/(-10**(en(1:n))*ctl%energy0)/2d0
	e(1:n)=sqrt(1-jum(1:n)**2)

    do i=1, fden%nbin
        do j=1, n
            call get_vr(10**fden%xb(i),a(j),e(j),ivr(j))
        end do
        ivrsum=sum(ivr*w)
        fden%fx(i)=ivrsum/pi/(4*pi*10**(fden%xb(i)*2))
    end do
    
end subroutine

subroutine update_weights()
    use com_main_gw
    implicit none
   ! print*, "begin update weights"
    !print*, bksams_arr%n
    !call all_chain_to_arr_single(bksams, bksams_arr)
    !print*, "1"
    !print*, bysams_arr%n
	!call all_chain_to_arr_by(bysams, bysams_arr)
    !print*, "2"
	call bcast_dms_asym_weights(dms)
    !print*, "3"
    !print*, "rid, dms%weight_asym=",rid, dms%weight_asym
	!call set_weights_for_all_samples(bksams,bysams, bksams_arr, bysams_arr)

   !print*, "finished update weights"

    call set_clone_weight(bksams)
end subroutine
subroutine set_weights_for_all_samples(sms_single, sms_arr_single)
    use com_main_gw
    implicit none
    type(particle_samples_arr_type)::sms_arr_single
    type(chain_type)::sms_single,sms_by
    integer i
    real(8) wtot
    
    call sams_get_weight_clone_single(sms_arr_single)
    
    !call set_asym_weight_arr(sms_arr_single, sms_arr_by)
    call set_real_weight_arr_single(sms_arr_single)
    !print*,rid, dms%asymp_bbh,dms%asymp_bbh_norm_factor

    call set_clone_weight(sms_single)
    !call set_asym_weight(sms_single)
    call set_real_weight(sms_single)
    
end subroutine

subroutine update_arrays_single()
    use com_main_gw
	implicit none
    !print*, "1"
	call all_chain_to_arr_single(bksams,bksams_arr)
    !print*, "2"
	call set_sample_arr_indexs_rid_particle(bksams_arr,rid)
    !print*, "3"
	call convert_sams_pointer_arr(bksams, bksams_pointer_arr,type=1)
    !print*, "4"
	call bksams_arr%select(bksams_arr_norm, exit_normal, -1d0, -1d0)
    
    if(allocated(bksams_arr%sp))deallocate(bksams_arr%sp)
    !call check_bksams()
end subroutine	

subroutine update_density()
    use com_main_gw
	implicit none
    type(particle_samples_arr_type)::sms_arr_star, sms_arr_sbh
    type(particle_samples_arr_type)::sms_arr_single, sms_arr_wd, sms_arr_ns,sms_arr_bd

    call sams_get_weight_clone_single(bksams_arr_norm)
    call set_real_weight_arr_single(bksams_arr_norm)

    call separate_to_species(bksams_arr_norm, sms_arr_star, &
     sms_arr_sbh, sms_arr_wd, sms_arr_ns,sms_arr_bd)

    !print*, "start get_gx", sms_arr_bbh%n, rid
    !print*, sms_arr_bbh%n
	!print*, "from update density"
    call gen_gx(dms, sms_arr_star,  sms_arr_sbh, &
		sms_arr_wd, sms_arr_ns,sms_arr_bd)
    call dms%get_dens0()
    !print*, "finished get_vm"
end subroutine


subroutine 	get_fna(fden, fna)
	use constant
	use com_sts_type
	use my_intgl
	use, intrinsic :: ieee_arithmetic
	implicit NONE
	type(s1d_type):: fden, fna
	integer i, idid
	real(8) int_out
    !print*, "start get_fna"
	do i=1, fna%nbin
		int_out=0
		call my_integral_none(fden%xmin,fna%xb(i), int_out, fcn,idid)
		fna%fx(i)=4*pi*int_out 
	end do
	!fna%nsam=fden%nsam

	contains
	subroutine fcn(n, x, y, f, ipar, par)
		implicit none
		integer n, ipar(100)
		real(8) x, y(n), f(n), par(100), y_out
		call fden%get_value_l(x, y_out)
		f(1)=y_out *(10**x)**3*log(10d0)
		if(isnan(f(1)).or..not.(ieee_is_finite(f(1))))then
			print*, "in fna"
			print*, "x, f=", x, f(1), y_out 
			stop
		end if

	end subroutine	
end subroutine

subroutine get_fden(gx, fden, n0, v0, rh,xmin)
	use constant
	use com_sts_type
	use my_intgl
	use, intrinsic :: ieee_arithmetic
	implicit none
	type(s1d_type)::gx, fden
	real(8) n0, v0,  rh, int_out,xmin
	integer i, j, idid
	do i=1, fden%nbin
		int_out=0
		call my_integral_none(10**xmin,rh/10**fden%xb(i), int_out, fcn,idid)
		fden%fx(i)=2/pi**0.5*int_out*n0 ! in unit of AU^{-3}
		!print*, fden%fx(i), int_out, n0
		!read(*,*)
	end do
contains
	subroutine fcn(n, x, y, f, ipar, par)
		implicit none
		integer n, ipar(100)
		real(8) x, y(n), f(n), par(100), y_out
		call gx%get_value_l( log10(x), y_out)
		f(1)=y_out*sqrt(abs(rh/10**fden%xb(i)-x))
		if(rh/10**fden%xb(i)-x<-1e-1)then
			print*, "error, rh/r-x<0", rh/10**fden%xb(i), rh, fden%xb(i), x
			stop
		end if
		if(isnan(f(1)).or..not.(ieee_is_finite(f(1))))then
			print*, "in fden"
			print*, "x, f=", x, f(1), y_out, rh/10**fden%xb(i)-x, rh, fden%xb(i)
			call gx%print("fden_gx")
			stop
		end if
		!call gx%print("gx")
		
		!read(*,*)
	end subroutine
end subroutine

subroutine get_den_u(fdenu,rmin,rmax, nbin, r0, n0)
	use com_sts_type
	use constant
	implicit none
	type(s1d_type)::fdenu
	real(8) rmin, rmax
	integer i, nbin
	real(8) n0,r0,r,phi

	!rmin=log10(0.5d0*rh/(emax_factor));rmax=log10(0.5d0*rh/(emin_factor))
	call fdenu%init(rmin,rmax,nbin,sts_type_grid)
	call fdenu%set_range()
	do i=1, nbin
		r=10**fdenu%xb(i)
		phi=r0/r
		if(phi<20)then
			fdenu%fx(i)=n0*(2d0/pi**0.5*phi**0.5+exp(phi)*erfc(phi**0.5))
		else 
			fdenu%fx(i)=n0*(2d0/pi**0.5*phi**0.5+1d0/pi**0.5/phi**0.5*(1-1d0/phi/2d0))
		end if
	end do
end subroutine

subroutine get_fxj(mb, jbtype)
    use com_main_gw
    implicit none
    type(mass_bins)::mb
    integer i, j, jbtype
    real(8) jm ,x
    !print*, "mb%all"

    call dms_so_get_fxj(mb%all, mb%n0, mb%mbh, mb%v0,jbtype)
    do i=1, n_tot_comp
        !print*, "i=",i
        call dms_so_get_fxj(mb%dsp(i)%p, mb%n0, mb%mbh, mb%v0,jbtype)
    end do
end subroutine
