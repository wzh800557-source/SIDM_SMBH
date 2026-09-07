
subroutine set_ini_byot_abin(abin)
	use com_main_gw
	implicit none
	real(8) abin, fpowerlaw!,alpha_ini
	!real(8) m
	!integer midx
	!call get_mass_idx(m,midx)
	!alpha_ini=
	abin=fpowerlaw(ctl%alpha_ini, 0.01d0*rh/(ctl%x_boundary*2), rhmax)
end subroutine
subroutine init_particle_sample_common(bkps)
	use com_main_gw
	implicit none
	type(particle_sample_type)::bkps
	real(8),external:: rnd, fpowerlaw,gen_ran_from_data
	
	bkps%en=-MBH/(2*bkps%byot%a_bin)
	bkps%en0=bkps%en
	bkps%jm=sqrt(1-bkps%byot%e_bin**2)
    bkps%jm0=bkps%jm
	bkps%djp=0;bkps%elp=0
	bkps%id=int(rnd(0d0,100000000d0))
	bkps%byot%ms%id=bkps%id

	bkps%state_flag_last=state_ae_evl

	bkps%exit_flag=exit_normal
	bkps%rid=rid
	!print*, "2"
	call set_particle_sample_other(bkps)
	bkps%byot_ini=bkps%byot
	call track_init(bkps,0)
end subroutine

subroutine set_particle_sample_other(bkps)
	use com_main_gw
	implicit none
	type(particle_sample_type)::bkps
	real(8) rnd, period
	
!	print*, "sp%byin%a_bin=",sp%byin%a_bin,sp%byin%e_bin, sp%byin%a_bin, sp%byin%e_bin
!	read(*,*)

	bkps%byot%ms%m=bkps%m; bkps%byot%mm%m=mbh; bkps%byot%mtot=bkps%m+mbh
	bkps%byot%Inc=rnd(0d0,pi); bkps%byot%Om=rnd(0d0,2*pi); 
    bkps%byot%pe=rnd(0d0, 2*pi); 
	bkps%byot%bname='byot'
	bkps%byot_bf%bname='byot_bf'
	bkps%byot%ms%id=bkps%id
	bkps%byot%mm%obtype=star_type_bh
	bkps%byot%mm%radius=mbh/1d8
    bkps%byot%me=rnd(0d0,2*pi)
    bkps%byot%an_in_mode=an_in_mode_mean
    call by_em2st(bkps%byot)
    call by_split_from_rd(bkps%byot)
end subroutine

subroutine get_sample_r_td_single(sp)
	use com_main_gw
    implicit none
    type(particle_sample_type)::sp
    real(8) wsi

	select case(sp%obtype)
	case(star_type_ms)
		if(sp%byot%ms%radius.eq.0)then
			print*, "ms radius=0?? check"
			stop
		end if
		sp%r_td=(3*mbh/sp%m)**(1/3d0)*sp%byot%ms%radius
	case(star_type_bh,star_type_ns,star_type_wd,star_type_bd)
		sp%r_td=16*mbh_radius/(1+sp%byot%e_bin)
		
	case default
		print*, "error! star type not defined",sp%obtype
		stop
	end select
end subroutine


subroutine set_star_radius(pr)
    use com_main_gw
    implicit none
    type(particle)::pr
    real(8) star_Radius,white_dwarf_radius
    select case(pr%obtype)
    case(star_type_MS)
        pr%radius=star_Radius(pr%m)
	case(star_type_BD)
		pr%radius=pr%m/my_unit_vel_c**2
    case(star_type_BH)
        pr%radius=pr%m/my_unit_vel_c**2
    case(star_type_WD)
        pr%radius=white_dwarf_radius(pr%m)
    case(star_type_NS)
        pr%radius=10d3/AU_SI
    case default
        !print*, "error, "
		print*, "star type=",pr%obtype
        pr%radius=0d0
		stop
    end select
end subroutine