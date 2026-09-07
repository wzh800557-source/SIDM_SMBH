module model_basic
    use md_chain_pointer
    use md_bk_species
	use md_chain
	use md_events
	use md_dms
	use constant
	use com_sts_type
	real(8),save:: MBH, rh,log10rh, mbh_radius 
    real(8) rhmin,rhmax
	type core_comp_type
		integer bk_evl_mode
		real(8) blkmass, gamma
		real(8) n0, alpha, N_in_rh, Mass_in_rh!,r0_vanish
        real(8) N_pre_in
		integer N_pre
		integer nbd_tot     ! number of  particles in outer boundary
		!integer n_test_tot
		real(8) rb_min, rb_max, frac, rb, rb_ini
		real(8) alpha_ini, rb_min_ini, rb_max_ini
		integer clone_amplifier
		!integer 
		character*(4) str_bktypemodel, str_bkmassmodel, str_bkacmodel, str_bkacinimodel
		integer::bktypemodel, bkmassmodel, bkacmodel, bkacinimodel
	end type
	!type ini_par_type
	!end type
	type(core_comp_type),pointer::cct_share
	type(s1d_type),pointer::fc_share
    real(8) fgx_g0
    integer, parameter::NJbin_tot=8
	type control_type
		type(core_comp_type),allocatable::cc(:)

		real(8):: mass_ref
		real(8):: bin_mass_min, bin_mass_max
        real(8) n_basic
        !real(8) weight_asym
		real(8):: bin_mass(20)
		real(8):: bin_mass_m1(20), bin_mass_m2(20), asymptot(8,20)
		real(8):: ini_weight_n(20)
		!asymptot 1=tot, 2=star, 3=sbh,  4=ns, 5=wd, 6=brown dwarf, 7=bstar, 8=bbh
        !integer:: bin_mass_Nbalance(10,NJbin_tot,4)=0  ! 1=star, 2=sbh, 3= bbh, 4=bstar
		
		real(8):: Weight_n(20)
		real(8):: alpha_ini 
        !real(8):: bin_mass_Nbalance(20,4)=0  ! 1=star, 2=sbh, 3= bbh, 4=bstar
        !integer:: bin_mass_N(10,NJbin_tot,4)=0
        real(8):: bin_mass_N(20,4)=0
		real(8)  total_time
       ! real(8) burn_in_time        
		real(8) sigma, energy0, energy_min, energy_max!, energy_boudary_min
		real(8) clone_e0,  v0, n0!, rh_min, rh_max
		real(8) rbd
		real(8) x_boundary, energy_boundary 
        real(8) ts_spshot, tnr,ts_snap_input    !ts_spshot_tnr, the shapshottime in unit of 
                                                !two body relaxation timescale, in unit of Myr
        !real(8) mass_bk_avg,n2a,na            !n2a:N(<2a), na:N(<a) 
		real(8) update_dt

		character*(3) str_jbin_bd, str_fj_bd
		character*(4) time_unit, str_model_intej
        character*(4)  str_method_interpolate
        character*(200) cfs_file_dir, burn_in_dir

		integer:: bin_mass_particle_number(20)  ! the number of particles at each bin
		integer:: clone_factor(20)
		integer EJ_mode
		!integer test_flag 
		integer:: num_bk_comp
        integer num_clone_created, num_boundary_created
		integer num_boundary_elim
        integer boundary_method, boundary_fj
		integer num_update_per_snap!, num_step_per_update
		integer n_spshot, n_spshot_total
		integer include_loss_cone
		integer model_intej
		integer::same_rseed_evl, same_rseed_ini	
		integer::clone_scheme 
		integer::trace_all_sample  
		!logical:: enable_collision_correction  : depleted feature!
		integer::del_cross_clone
        integer burn_in_snap
		integer jbin_type
		integer idxstar, idxsbh,  idxns, idxwd,idxbd 
		
		integer method_interpolate
		!integer::npar_sam_tot
		integer::m_bins ! number of mass bin
		integer::grid_bins, gx_bins,grid_type
        integer::idx_ref
		integer:: debug=0
        integer:: bin_mass_flux_in(20,4)=0  ! 1=star, 2=sbh, 3= bbh, 4=bstar
        integer:: bin_mass_flux_out(20,4)=0  ! 1=star, 2=sbh, 3= bbh, 4=bstar
        integer:: bin_mass_emax_out(20,4)=0
		integer:: chattery ! <=2, normal; =3, add tdial details
        integer:: ntasks, ntask_total 
        integer:: seed_value, ntask_bg
		integer::nblock_mpi_bg, nblock_mpi_ed, nblock_size
		integer::output_track_td, output_track_emri
		integer::output_track_plunge
		!real(8) blkmass    ! blkmass:back ground mass
        logical burn_in_phase

	end type
	!type(ini_par_type)::ipt
	type(control_type),target::ctl
	type(chain_type)::Allsams

	type(particle_samples_arr_type)::bksams_arr_ini
	type(chain_type):: bkstars, bkbystars, bkbbhs, bksbhs
	type(chain_type):: bksams  ! includes only particle samples
	type(chain_type):: bksams_norm,bksams_merge

	type(particle_samples_arr_type):: bkstars_arr, bkbystars_arr
	type(particle_samples_arr_type):: bkbbhs_arr, bksbhs_arr
	type(particle_samples_arr_type):: bksbhs_arr_norm
	type(particle_samples_arr_type),target:: bksams_arr_norm   ! only normal samples
	type(particle_samples_arr_type),target:: bksams_arr_norm_sbh
	type(particle_samples_arr_type):: bksams_arr    ! include all samples
	!type(particle_samples_arr_type):: bksams_arr_merge  ! due to gw capture
	type(samples_type_pointer_arr),target::bksams_pointer_arr
	type(diffuse_mspec),target::dms

	!real(8),allocatable:: rcolld(:)
	!real(8),allocatable:: rc_weight(:)
	integer, parameter::ncd_maxsize=200000
    integer,parameter:: chattery_out_unit_0=1383829393
    integer chattery_out_unit

	real(8)::dc_grid_xstep, dc_grid_ystep
	!real(8) rcolld(ncd_maxsize)
	!real(8) rc_weight(ncd_maxsize)
    !real(8),allocatable::lambda_aux(:)
    real(8),parameter::my_unit_vel_c5=my_unit_vel_c**5
	real(8),parameter::my_unit_vel_c3=my_unit_vel_c**3
    
	real(8) clone_e0_factor     !the reference postion of energy
	real(8):: clone_emax

	real(8)::jmin_value, jmax_value
	real(8)::emin_value, emax_value
	integer,parameter::method_int_nearst=1, method_int_linear=2
		 
	integer,parameter::boundary_fj_iso=1, boundary_fj_ls=2
	integer,parameter::model_intej_fast=1 
	integer,parameter::task_mode_new=1, task_mode_append=2
	integer,parameter::snap_mode_new=1, snap_mode_append=2
    
    integer,parameter::outcome_ejection=15
	integer,parameter::ls_type_compact=1, ls_type_binary=2
	integer,parameter::MAX_LENGTH=100000
!    integer,parameter::default_track_size=10000
    integer(8),parameter::MAX_RUN_LENGTH=int(1d7) 
	integer,parameter::n_orb_track_block=100000
    integer,parameter::n_orb_track_block_max=n_orb_track_block*16
	real(8)::ini_bhb_mbh_r0, ini_by_3body_r0
	integer,parameter::collection_data_through_mpi=0
    integer nsize_chain_bk, nsize_chain_by, nsize_arr_bk, nsize_arr_by, &
    nsize_arr_bk_norm, nsize_arr_by_norm, nsize_arr_bk_pointer, &
    nsize_arr_by_pointer, nsize_tot_bk, nsize_tot_by, nsize_tot

    integer:: boundary_sts_emin_cros=0,boundary_sts_emin_dir=0
    integer:: boundary_sts_create=0,boundary_sts_replace=0
    integer::boundary_cross_up=0,boundary_cross_down=0
    integer:: boundary_cross_net=0,boundary_sts_emax_cros=0
contains
	real(8) function p(a)
	IMPLICIT NONE
		real(8) a
		P=2*PI*sqrt(a**3/Mbh)
	end function
	
	
	subroutine all_chain_to_arr_single(sps, sps_arr)
		implicit none
		type(chain_type),intent(in)::sps
		type(particle_samples_arr_type),intent(out)::sps_arr
	!	type(chain_type)::sps_tmp
		type(chain_pointer_type),pointer::sp
		integer i,n
		integer::nr
	!	print*, associated(sps_tmp%sp)
	!	call copy_sams(sps,sps_tmp)

		call sps%get_length(nr,type=1)

        !print*, "nr=", nr
        if(nr.eq.0)return
!		sp=>sps_tmp%sp(1)
!		print*, sps_tmp%sp(1)%ob%ac

		call sps_arr%init(nr)
		
        sp=>sps%head
		call sp%chain_to_arr_single(sps_arr%sp(1:nr),nr)

	end subroutine

	
	subroutine get_exit_flag_str(exit_flag, str_flag)
		implicit none
		integer exit_flag
		character*(100) str_flag
 !   integer,parameter::exit_normal=1, exit_tidal=2,exit_max_reach=3
!	integer,parameter::exit_boundary_min=4, exit_boundary_max=5 !, !exit_gw_capture=6
!!	integer,parameter::exit_invtransit=7, exit_tidal_empty=8, exit_tidal_full=9
!	integer,parameter::exit_ejection=11, exit_by_ionization=10, exit_gw_kl=12
!	integer,parameter::exit_by_exchange=13,exit_gw_bhb_mbh_enc=14
	!integer,parameter::exit_gw_3body_enc_byself=15,exit_gw_3body_enc_with_incoming=16

		select case(exit_flag)
		
		case (exit_normal)
			str_flag="NORMAL"
		case (exit_tidal)
			str_flag="TIDAL"
		case (exit_max_reach)
			str_flag="NMAX"
		case (exit_boundary_min)
			str_flag="BOUNDAY_MIN"
		case(exit_boundary_max)
			str_flag="BOUNDARY MAX"
		!case(exit_gw_capture)
		!	str_flag="GW CAPTURE"
		case(exit_invtransit)
			str_flag="INVERSE TRANSIT"
		case(exit_tidal_empty)
			str_flag="TD EMPTY"
		case(exit_tidal_full)
			str_flag="TD FULL"
		case(exit_ejection)
			str_flag="EJECT"
		case (exit_gw_iso)
			str_flag="GW_ISO"
		case default
			str_flag="Null"
		end select
	end subroutine

	integer function get_lvl(en)
        implicit none
        real(8) en
        get_lvl=int(log10(en/ctl%clone_e0))
    end function

	subroutine set_seed(same_seed,seed_value)
		implicit none
		integer same_seed
		integer seed_value
		if(same_seed>0)then
			call same_random_seed(seed_value)
		else
			call random_seed()
		!   call same_random_seed(seed_value+rid)
		end if
	end subroutine
end module


