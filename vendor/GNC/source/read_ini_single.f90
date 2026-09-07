
subroutine print_current_code_version()
	implicit none
	print*, "current code version: 1.0 (updated in 2023-8-13)"
end subroutine

subroutine readin_model_par(fl)
	use com_main_gw
	implicit none
	character*(*) fl
	character*(200) model, inacmodel, outacmodel, ecmodel, dejmodel, model_intej, str_method_int
    character*(200) time_unit, byctype, bytype ,taskmode
    real(8) mass(2)
    integer i,ier
	character*(200) strall, str(12), str_massbin_mode, &
		str_sg_data_mode,str_jbin_bd,str_fj_bd, str_plunge_cr
	integer nnum	
 !   real(8) ts_snap_input
    !print*, "0"
	open(unit=999,file=fl,status='old')
    call readpar_int_sp(ctl%ntasks, 999, "#","=",ier)
	!call readpar_int_sp(ctl%test_flag,999, "#", "=", ier)
	call readpar_str_sp(model_intej,999, "#","=")
	ctl%str_model_intej=trim(adjustl(model_intej))
	select case (trim(adjustl(ctl%str_model_intej)))
	case ("FAST")
		ctl%model_intej=model_intej_fast
	case default
		print*, "ERROR, define dejmodel:", trim(adjustl(ctl%str_model_intej))
		stop
	end select
	!print*, "1"
	call readpar_str_sp(str_method_int,999, "#","=")
	ctl%str_method_interpolate=trim(adjustl(str_method_int))
	select case(trim(adjustl(str_method_int)))
	case("near")
		ctl%method_interpolate=method_int_nearst
	case("li2d")
		ctl%method_interpolate=method_int_linear
	end select
	call READPAR_STR_AUTO_sp(strall,str,2,Nnum,999,"#",","," ")
	!print*, "str=",str(1),str(2)
	read(str(1),fmt=*) jmin_value
	read(str(2),fmt=*) jmax_value
	!print*, "jmax=",jmax_value
    call readpar_str_sp(ctl%cfs_file_dir, 999, "#", "=")
	call readpar_dbl_sp(mbh,999, "#","=",IER)
    mbh_radius=mbh/(my_unit_vel_c**2)
	
	call readpar_dbl_sp(emin_factor, 999, "#","=",IER)
    log10emin_factor=log10(emin_factor)
    call readpar_dbl_sp(emax_factor, 999, "#","=",IER)
    log10emax_factor=log10(emax_factor)
    call readpar_dbl_sp(ctl%x_boundary, 999, "#","=",IER)
	call readpar_str_sp(str_jbin_bd,999, "#","=")
	ctl%str_jbin_bd=trim(adjustl(str_jbin_bd))
	select case(trim(adjustl(ctl%str_jbin_bd)))
	case("LIN")
		ctl%jbin_type=Jbin_type_lin
	case("LOG")
		ctl%jbin_type=Jbin_type_log
	case("SQR")
		ctl%jbin_type=jbin_type_sqr
	case default
		print*, "error! define jbin type", ctl%str_jbin_bd
		stop
	end select
	
	call readpar_str_sp(str_fj_bd,999,"#","=")
	ctl%str_fj_bd=trim(adjustl(str_fj_bd))
	select case(trim(adjustl(ctl%str_fj_bd)))
	case("ISO")
		ctl%boundary_fj=boundary_fj_iso
	case("LS")
		ctl%boundary_fj=boundary_fj_ls
	case default
		print*, "error! define fj type", ctl%str_fj_bd
		stop
	end select
	!  
	
	call readpar_int_sp(ctl%seed_value, 999, "#","=",ier)
	!print*, "seed_value=", ctl%seed_value
	call readpar_int_sp(ctl%same_rseed_ini, 999, "#","=",ier)
	
	call READPAR_STR_AUTO_sp(strall,str,2,Nnum,999,"#",","," ")
	read(str(1),fmt=*) ctl%gx_bins
	read(str(2),fmt=*) ctl%grid_bins
	if(mod(ctl%grid_bins,ctl%ntasks).ne.0)then
		print*, "read_ini_single.f90: 101: grid_bins, ntasks=", ctl%grid_bins, ctl%ntasks
		print*, "read_ini_single.f90: 101: error! grid_bin number must be integer times of ntasks"
		stop
	end if

	call readpar_int_sp(ctl%same_rseed_evl, 999, "#","=",ier)
	!print*, "1??"
	!ctl%nbksam_tot=ctl%nbksam_task*(ctl%ntask_bg+ctl%ntasks)
	call readin_task_mode(999)
	call readin_snap_mode(999)
	call readpar_dbl_sp(ctl%alpha_ini, 999, "#", "=", ier)
	call readpar_int_sp(ctl%include_loss_cone, 999, "#","=",ier) 

	rtmin=1

	call readpar_int_sp(ctl%clone_scheme, 999, "#","=",ier)
!	print*, ctl%clone_scheme
	if(ctl%clone_scheme.ge.1)then
	!	call readpar_int_sp(clone_amplifier, 999, "#","=")
	    call readpar_dbl_sp(clone_e0_factor, 999, "#","=",IER)
	end if
	
	ctl%del_cross_clone=ctl%clone_scheme
	call readpar_int_sp(ctl%chattery, 999, "#","=",ier)
    call readpar_int_sp(ctl%trace_all_sample, 999, "#","=",ier)
	if(ctl%trace_all_sample.ge.1)then
		call readpar_int_sp(ctl%output_track_td, 999, "#","=",ier)
		call readpar_int_sp(ctl%output_track_plunge, 999, "#","=",ier)
		call readpar_int_sp(ctl%output_track_emri, 999, "#","=",ier)
	end if
	close(999)
    print*, "single finished"
	call readin_mass_bins("mfrac.in")
    
    call check_readin()
end subroutine
subroutine readin_mass_bins(fl)
	use com_main_gw
	implicit none
	character*(200) strall, str(12), str_massbin_mode
	integer nnum, funit, i,j, min_particle_n, num_comp
	logical have_comp(7)
	character*(*) fl
	real(8) mb(20),m1(20),m2(20),mstep, mass_imf, norm, max_weight_n
	real(8) wn_close

	funit=1999
	open(unit=funit,file=fl,status='old')
	!call readin_clone_factor(funit)
	call READPAR_STR_AUTO_sp(strall,str,12,Nnum,funit,"#",","," ")
	!print*, trim(adjustl(str(1)))
	!print*, trim(adjustl(str(2)))
	!print*, trim(adjustl(str(3)))
	!print*, trim(adjustl(str(4)))

	read(str(1),fmt=*) ctl%m_bins
	read(str(2),fmt=*) str_massbin_mode
	select case (trim(adjustl(str_massbin_mode)))
	case("GIVEN")
		do i=1, ctl%m_bins
			call skip_comments(funit,"#")
			read(funit,fmt=*) ctl%bin_mass_m1(i),ctl%bin_mass(i),ctl%bin_mass_m2(i),&
				ctl%asymptot(1,i), ctl%ini_weight_n(i),ctl%clone_factor(i) 
			read(funit,fmt=*) ctl%asymptot(2:6,i)
		end do
        do i=2,ctl%m_bins
            if(ctl%bin_mass_m1(i).eq.ctl%bin_mass_m2(i-1))then
                print*, "error! bin_mass_m1=bin_mass_m2", i
                print*, "m1=",ctl%bin_mass_m1(i), ctl%bin_mass_m2(i-1)
                stop
            end if
        end do
		!ctl%bin_mass_mode=bin_mass_mode_given
		have_comp=.false.
		num_comp=0
loopj:	do j=1, 7		
loopi:		do i=1, ctl%m_bins
				if(ctl%asymptot(j+1,i).eq.0)then
					cycle loopi
				else
					have_comp(j)=.true.
					num_comp=num_comp+1
					cycle loopj
				end if
			end do loopi
		end do loopj
		
		ctl%num_bk_comp=num_comp
		allocate(ctl%cc(num_comp))
		j=0
		do i=1, 7
			if(have_comp(i))then
				j=j+1
				select case(i)
				case(1)
					ctl%cc(j)%bktypemodel=star_type_MS
				case(2)
					ctl%cc(j)%bktypemodel=star_type_BH
				case(3)
					ctl%cc(j)%bktypemodel=star_type_NS
				case(4)
					ctl%cc(j)%bktypemodel=star_type_WD
				case(5)
					ctl%cc(j)%bktypemodel=star_type_BD										
				case default
				end select
			end if
		end do
	case default
		print*, "readin_mass_bins:error!"
		stop
	end select
	!ctl%bin_mass_Nbalance=0
    !ctl%bin_mass_Nbalance_in=0
    !ctl%bin_mass_Nbalance_ot=0
	!print*, ctl%grid_bins, ctl%gx_bins
	!read(*,*)

    ctl%idx_ref=1
    ctl%mass_ref=ctl%bin_mass(ctl%idx_ref)

	close(unit=funit)
	
end subroutine
subroutine find_close_number(nin, nout)
	implicit none
	real(8) nin, nout
	integer dg
	if(nin<=0)then
		print*, "nin should be larger than zero! stoped!"
		stop
	end if
	dg=int(log10(nin))+1
	if(nin>5*10**(dg-1)) then
		nout=5*10**(dg-1)
	else
		nout=10**(dg-1)
	end if
	print*, "nin, nout=",nin, nout
end subroutine
subroutine readin_snap_mode(funit)
	use com_main_gw
	implicit none
	character*(200) snapmode,time_unit
	integer funit
    integer ier

	!call readpar_str_sp(snapmode,funit,'#',"=")
	!ctl%str_snap_mode=trim(adjustl(snapmode))
	!select case (trim(adjustl(snapmode)))
	!case ("NEW")
		!ctl%insnapmode=snap_mode_new
        !print*, "start burn in reading"
        !call readpar_int_sp(ctl%num_step_per_update, funit, "#", "=", ier)
		call readpar_int_sp(ctl%num_update_per_snap, funit, "#", "=",ier)
		call readpar_dbl_sp(ctl%ts_snap_input, funit, "#","=",ier)
        call readpar_str_sp(time_unit,funit,"#","=")
		call readpar_int_sp(ctl%n_spshot,funit,'#',"=",ier)
		!ctl%n_spshot_bg=0
        ctl%time_unit=trim(adjustl(time_unit))
	!case ("APP")
	!	!ctl%insnapmode=snap_mode_append
    !    ! ts_spshot is in unit of Tnr
    !    !call readpar_int_sp(ctl%burn_in_snap, funit, "#", "=",ier)
    !    !call readpar_dbl_sp(ctl%burn_in_time, funit, "#", "=",ier)
    !    !call readpar_int_sp(ctl%num_step_per_update, funit, "#", "=", ier)
	!	call readpar_int_sp(ctl%num_update_per_snap, funit, "#", "=",ier)
	!	call readpar_dbl_sp(ctl%ts_snap_input, funit, "#","=",ier)
    !    call readpar_str_sp(time_unit,funit,"#","=")
	!	call readpar_int_sp(ctl%n_spshot,funit,'#',"=",ier)
    !    
	!	call readpar_int_sp(ctl%n_spshot_bg,funit,'#',"=",ier)
    !    ctl%time_unit=trim(adjustl(time_unit))
	!case default
	!	print*, "define snapmode ", trim(adjustl(snapmode))
	!	stop
	!end select
	
    !print*, "burn in =", ctl%burn_in_snap
end subroutine
subroutine set_simu_time()
	use com_main_gw
	implicit none

	!select case(trim(adjustl(ctl%time_unit)))
	!case ("TNR")
        call get_tnr_timescale_at_rh(ctl%tnr)
        ctl%ts_spshot=ctl%tnr*ctl%ts_snap_input
		!ctl%burn_in_time=ctl%tnr*ctl%burn_in_time
		ctl%update_dt=ctl%ts_spshot/dble(ctl%num_update_per_snap)
		ctl%n_spshot_total=ctl%n_spshot!+ctl%n_spshot_bg
        ctl%total_time=ctl%ts_spshot*ctl%n_spshot_total
		!print*,"tot,ts,n=",ctl%total_time,ctl%ts_spshot,ctl%n_spshot_total
	!case ("Myr")
	!	ctl%update_dt=ctl%ts_snap_input/dble(ctl%num_update_per_snap)
	!	ctl%n_spshot_total=ctl%n_spshot+ctl%n_spshot_bg
	!	ctl%ts_spshot=ctl%ts_snap_input
	!	ctl%total_time=ctl%ts_spshot*ctl%n_spshot_total
	!end select
end subroutine
subroutine readin_task_mode(funit)
	use com_main_gw
	implicit none
	character*(200) taskmode
	integer funit, ier

	!call readpar_str_sp(taskmode,999,'#',"=")
	!	print*, "taskmode=",taskmode
	!ctl%str_task_mode=trim(adjustl(taskmode))
	!select case (trim(adjustl(taskmode)))
	!case ("NEW")
		!ctl%intaskmode=task_mode_new
		ctl%ntask_bg=0
	!case ("APP")
	!	ctl%intaskmode=task_mode_append		
	!	call readpar_int_sp(ctl%ntask_bg,999,"#","=",ier)
	!case default
	!	print*, "define taskmode", trim(adjustl(taskmode))
	!	stop
	!end select
	ctl%ntask_total=ctl%ntask_bg+ctl%ntasks
end subroutine

subroutine print_model_par()
	use com_main_gw
	implicit none
	integer i,k
    !write(unit=*,fmt="(A25, L30)") "test j=", test_no_j_relaxation
	!write(unit=*,fmt="(A25, A30)") "DEJ_GRID_FILE=", adjustl(trim(DEJ_GRID_FILE_DIR))
	write(unit=*,fmt="(A25, I10)") "num of components=", ctl%num_bk_comp
	write(unit=*,fmt="(A25, 1P4E10.3)") "MBH,rh,rhmin,rhmax=",MBH, rh, rhmin,rhmax
	write(unit=*,fmt="(A25, 1P4E10.3)") "Jmin, Jmax=",jmin_value, jmax_value
	write(unit=*,fmt="(A25, 1P3E10.3)") "ExMIN,ExMAX, EB=", emin_factor, emax_factor, ctl%x_boundary
	write(unit=*,fmt="(A25, A25)") "Integration MODEL=", trim(adjustl(ctl%str_model_intej))
	write(unit=*,fmt="(A25, A25)") "Intepolate Method=", trim(adjustl(ctl%str_method_interpolate))
	if(ctl%include_loss_cone.ge.1)then
    	print*, "include_loss_cone=",ctl%include_loss_cone
	end if
	
	write(unit=*, fmt="(A25)") "mass bin mode: given"
    write(unit=*, fmt="(A25, 1P20E10.3)") "mass_bin=", (ctl%bin_mass(k),k=1, ctl%m_bins) 
	do i=1, ctl%m_bins
		write(unit=*,fmt="(A25,I8, 2F10.3, 2I10)") "I,IWN,WN,PN,CL=",&
			 i, ctl%ini_weight_n(i), ctl%weight_n(i), ctl%bin_mass_particle_number(i), ctl%clone_factor(i)
	end do
	write(unit=*,fmt="(A25,F10.3)") "ctl%N_basic=", ctl%N_basic
	write(unit=*,fmt="(A25, 3I10)") "NTASK, NTASK_TOT=", ctl%ntasks, ctl%ntask_total
	write(unit=*,fmt="(A25, 3I10)") "NSNAP, NSNAP_TOT=", ctl%n_spshot, ctl%n_spshot_total
	write(unit=*,fmt="(A25, A4,1F10.2)") "TIMEUNIT,TNR=", ctl%time_unit, ctl%Tnr
	write(unit=*,fmt="(A25, 1F10.2, 2A10)") "TOTAL_TIME=", ctl%total_time
	write(unit=*,fmt="(A25, 4L10)") "CLONE=", ctl%clone_scheme
	write(unit=*,fmt="(A25, 3L10)") "SEVSEED=",ctl%same_rseed_evl
	
end subroutine
subroutine check_readin()
    use com_main_gw
    implicit none
    integer idxsbh, idxstar

	call get_ccidx_from_type(idxstar,star_type_ms)
	call get_ccidx_from_type(idxsbh,star_type_bh)
end subroutine