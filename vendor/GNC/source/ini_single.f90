subroutine get_init_samples(bksps_arr_ini)
	use com_main_gw
	implicit none
	type(particle_samples_arr_type)::bksps_arr_ini
	call get_init_samples_given(bksps_arr_ini)
end subroutine

subroutine init_particle_sample_one(sample,m,flag)
	use com_main_gw
	implicit none
	type(particle_sample_type)::sample
	real(8) m
	integer flag,imass

	sample%m=m
	!print*, "init_particle_sample_one type:", sample%obtype, sample%obidx
	call init_particle_sample_one_model_rnd(sample, flag)

end subroutine

subroutine get_numbers_each_bin(nstar_tot,nsbh_tot,nns_tot,nwd_tot,nbd_tot,&
	nstar,nsbh,nns,nwd,nbd,n)
	use com_main_gw
	implicit none
	integer n, nstar_tot, nsbh_tot,nns_tot,nwd_tot,nbd_tot
	integer i, j, nstar(n), nbstar(n), nbd(n)
	integer nsbh(n), nns(n), nwd(n)
	real(8) factor
	
	!print*, "idxwd=",idxwd
    !factor=mbh!/ctl%bin_mass(ctl%idx_ref)!*(ctl%x_boundary/emin_factor)**(3-7d0/4d0)
	
    do i=1, ctl%m_bins
		!print*, "i=",i
		!print*, "aymptot(:)=",ctl%asymptot(:,i)
        nstar(i)=ctl%asymptot(2,i)*ctl%bin_mass_particle_number(i)
        
        if(ctl%idxsbh.ne.-1)then
		    nsbh(i)=ctl%asymptot(3,i)*ctl%bin_mass_particle_number(i)
        else
            nsbh(i)=0
        end if
		if(ctl%idxns.ne.-1)then
			nns(i)=ctl%asymptot(4,i)*ctl%bin_mass_particle_number(i)
		else
			nns(i)=0
		end if
		if(ctl%idxwd.ne.-1)then
			nwd(i)=ctl%asymptot(5,i)*ctl%bin_mass_particle_number(i)
		else
			nwd(i)=0
		end if
		!print*, "ctl%cc(idxbd)%N_pre_in=", ctl%cc(idxbd)%N_pre_in,ctl%asymptot(6,i)
		if(ctl%idxbd.ne.-1)then
			nbd(i)=ctl%asymptot(6,i)*ctl%bin_mass_particle_number(i)
		else
			nbd(i)=0
		end if
		
    end do

    nstar_tot=sum(nstar)
    if(ctl%idxsbh.ne.-1)then
	    nsbh_tot=sum(nsbh)
    else
        nsbh_tot=0
    end if
	if(ctl%idxwd.ne.-1)then
	    nwd_tot=sum(nwd)
    else
        nwd_tot=0
    end if
	if(ctl%idxns.ne.-1)then
	    nns_tot=sum(nns)
    else
        nns_tot=0
    end if
	if(ctl%idxbd.ne.-1)then
	    nbd_tot=sum(nbd)
    else
        nbd_tot=0
    end if
	
	print*, "nstar_tot=", nstar_tot
	print*, "nsbh_tot=", nsbh_tot
	print*, "wd_tot,ns_tot,bd_tot=", nwd_tot, nns_tot, nbd_tot
end subroutine
subroutine get_init_samples_given(bksps_arr_ini)
	use com_main_gw
	implicit none
	type(particle_samples_arr_type)::bksps_arr_ini
	integer i, j, nstar(ctl%m_bins),  nbd(ctl%m_bins)
	integer nsbh(ctl%m_bins), nns(ctl%m_bins), nwd(ctl%m_bins)
	integer nsg0
	integer nstar_tot,   nsbh_tot, nwd_tot, nns_tot, nbd_tot

	call get_numbers_each_bin(nstar_tot, nsbh_tot,nns_tot,nwd_tot,nbd_tot,&
	nstar,nsbh,nns,nwd,nbd,ctl%m_bins)

	call bksps_arr_ini%init(nstar_tot+nsbh_tot+nwd_tot+nns_tot+nbd_tot)
	!print*, "bksps_arr_ini%n=", bksps_arr_ini%n, nstar(1)
    nsg0=0; 
	do i=1, ctl%m_bins
		do j=1, nstar(i)
			bksps_arr_ini%sp(j+nsg0)%obtype=star_type_ms
			bksps_arr_ini%sp(j+nsg0)%obidx=ctl%idxstar
			bksps_arr_ini%sp(j+nsg0)%m=ctl%bin_mass(i)
			call init_particle_sample_one_model_rnd(bksps_arr_ini%sp(j+nsg0), flag_ini_ini)
            !print*, "star%a_bin=", bksps_arr_ini%sp(j+nsg0)%byot%a_bin,rid
		end do
		nsg0=nsg0+nstar(i)
		!=============================================================================
		do j=1, nsbh(i)
			bksps_arr_ini%sp(j+nsg0)%obtype=star_type_bh
			bksps_arr_ini%sp(j+nsg0)%obidx=ctl%idxsbh
			!bksps_arr_ini%sp(j+nsg0)%m=ctl%bin_mass(i)
            bksps_arr_ini%sp(j+nsg0)%m=ctl%bin_mass(i)
			call init_particle_sample_one_model_rnd(bksps_arr_ini%sp(j+nsg0), flag_ini_ini)
           ! print*, j+nsg0, bksps_arr_ini%sp(j+nsg0)%exit_flag
		end do
		nsg0=nsg0+nsbh(i)
        !stop
		!=============================================================================
		do j=1, nwd(i)
			bksps_arr_ini%sp(j+nsg0)%obtype=star_type_wd
			bksps_arr_ini%sp(j+nsg0)%obidx=ctl%idxwd
            bksps_arr_ini%sp(j+nsg0)%m=ctl%bin_mass(i)
			call init_particle_sample_one_model_rnd(bksps_arr_ini%sp(j+nsg0), flag_ini_ini)
           ! print*, j+nsg0, bksps_arr_ini%sp(j+nsg0)%exit_flag
		end do
        !stop
		nsg0=nsg0+nwd(i)
		!=============================================================================
		do j=1, nns(i)
			bksps_arr_ini%sp(j+nsg0)%obtype=star_type_ns
			bksps_arr_ini%sp(j+nsg0)%obidx=ctl%idxns
            bksps_arr_ini%sp(j+nsg0)%m=ctl%bin_mass(i)
			call init_particle_sample_one_model_rnd(bksps_arr_ini%sp(j+nsg0), flag_ini_ini)
		end do
        !stop
		nsg0=nsg0+nns(i)
		!=============================================================================
		do j=1, nbd(i)
			bksps_arr_ini%sp(j+nsg0)%obtype=star_type_bd
			bksps_arr_ini%sp(j+nsg0)%obidx=ctl%idxbd
            bksps_arr_ini%sp(j+nsg0)%m=ctl%bin_mass(i)
			call init_particle_sample_one_model_rnd(bksps_arr_ini%sp(j+nsg0), flag_ini_ini)
		end do
        !stop
		nsg0=nsg0+nbd(i)
		!=================================================================
	end do
   ! print*, "2"
end subroutine

real(8) function star_Radius(mass)
use constant
IMPLICIT NONE
real(8) mass
if(mass.le.0.06d0)then
	star_Radius=0.1*rd_sun  ! Chabrier, G. \& Baraffe, I.\ 2000, \araa, 38, 337. Figure 3
elseif(mass.gt.0.06d0.and.mass.le.1)then
	star_Radius=mass**0.8d0*rd_sun  ! Stellar Structure and Evolution, Rudolf Kippenhahn, Alfred Weigert, Achim Weiss P253
else
	star_Radius=mass**0.56d0*rd_sun  ! Stellar Structure and Evolution, Rudolf Kippenhahn, Alfred Weigert, Achim Weiss P253
end if
end function

real(8) function white_dwarf_radius(mass)
	use constant
	implicit none
	real(8) mass
	if(mass<1.44d0)then
		white_dwarf_radius=0.01*rd_sun*mass**(-1/3d0)
	else
		print*, "error, white dwarf mass should be smaller than 1.44 solar mass"
	end if
end function

subroutine set_chain_samples(cbk, bksps_arr)
	use com_main_gw
	IMPLICIT NONE
	type(particle_samples_arr_type)::bksps_arr
	type(chain_type)::cbk

	call set_chain_samples_single(cbk, bksps_arr)
end subroutine

subroutine set_chain_samples_single(cbk, bksps_arr)
	use com_main_gw
	IMPLICIT NONE
	integer i,j, flag, obtype, typeidx,n
	type(chain_pointer_type),pointer::ptbk, ptby
	type(particle_samples_arr_type)::bksps_arr
	type(chain_type)::cbk
	!integer,parameter::flag_sg=1,flag_by=2

   ! bksps_arr%n=10
    !bysps_arr%n=10
	call cbk%init(bksps_arr%n)
    ptbk=>cbk%head
	do i=1, bksps_arr%n
		if(.not.allocated(ptbk%ob)) allocate(particle_sample_type::ptbk%ob)
        select type (ca=>ptbk%ob)
        type is(particle_sample_type)
            ca=bksps_arr%sp(i)
            ca%create_time=0d0
            ca%simu_bgtime=0d0
            !if(ctl%clone_scheme.ge.1)then
            !    ca%nhiar=get_lvl(ca%en)
            !else
            !    ca%nhiar=0
            !end if
            ca%en0=-mbh/2d0/ca%byot%a_bin
            ca%jm0=sqrt(1-ca%byot%e_bin**2)
			if(ca%jm0<jmin_value)then
				print*, "jm0=",ca%jm0
				stop
			end if
			!print*,"jm, ebin=", ca%jm, ca%byot%e_bin
            if(ctl%clone_scheme.ge.1)then
                !print*, "clone"
                call create_init_clone_particle(ptbk, ca%en0,0d0)
            end if
		end select
        ptbk=>ptbk%next
	end do
    
end subroutine

subroutine init_particle_sample_one_model_rnd(bkps, flag)
	use com_main_gw
	implicit none
	type(particle_sample_type)::bkps
	real(8) ecmax,logpd, beta, GET_T_GW
	real(8),external:: rnd, fpowerlaw,gen_ran_from_data
	type(core_comp_type),pointer::cc
	integer flag,midx
    
	if(bkps%obtype.eq.0.or.bkps%obidx.eq.0)then
		print*, "error:particle type not assigned", bkps%obtype, bkps%obidx
		stop
	end if

	if(bkps%m.le.0)then
		print*, "error:particle mass should be assigned", bkps%m
		stop
	end if
	!cc=>ctl%cc(bkps%obidx)!	  
	
100 select case (flag)
	case(flag_ini_or)
        if(bkps%byot%a_bin.eq.0d0) then
            print*, "error! flag_ini_or should assume abin first"
            stop
        end if

		bkps%byot%ms%m=bkps%m
		
		bkps%byot%ms%obtype=bkps%obtype
		bkps%byot%ms%obidx=bkps%obidx

		call set_star_radius(bkps%byot%ms)

		!select case(ctl%include_loss_cone)
		!case (1)
		if(bkps%byot%ms%radius.eq.0) then
			print*, "ini", bkps%byot%ms%radius
			stop
		end if
		call get_sample_r_td(bkps)

		bkps%rp=bkps%byot%a_bin*(1-bkps%byot%e_bin)
		
		
		call init_particle_sample_common(bkps)
		if(bkps%weight_real.eq.0)then
			print*, "error in ini run particle: bkps%weight_real=0,id=", &
			bkps%id, bkps%weight_real, bkps%weight_n, &
			bkps%weight_clone
			stop
		end if

	case(flag_ini_ini)
		
		call set_ini_byot_abin(bkps%byot%a_bin)		
		!bkps%m=m_input
		bkps%byot%ms%m=bkps%m
		
		bkps%byot%ms%obtype=bkps%obtype
		bkps%byot%ms%obidx=bkps%obidx

		call set_star_radius(bkps%byot%ms)

		!select case(ctl%include_loss_cone)
		!case (1)
		if(bkps%byot%ms%radius.eq.0) then
			print*, "ini", bkps%byot%ms%radius
			stop
		end if
		call get_sample_r_td(bkps)
		call set_jm_init(bkps)

		bkps%rp=bkps%byot%a_bin*(1-bkps%byot%e_bin)	
		
		call init_particle_sample_common(bkps)

		call get_mass_idx(bkps%m,midx)
		bkps%weight_n=ctl%Weight_n(midx)
        !bkps%weight_asym=1d0
    case default
        print*, "undefined flag in init_particle_sample_one_model_rnd"
        stop
	end select

end subroutine

