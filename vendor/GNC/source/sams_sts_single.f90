
subroutine get_sams_events_sbh(sa, sams_sel_norm, sams_sel_norm_inbd,sams_emri_by,& 
		sams_emri_sg,sams_emri, sams_plunge,  eveset,isnap)
	use com_main_gw
	implicit none
	type(particle_samples_arr_type)::sa, sams_sel, sams_sel_norm, sams_sel_norm_inbd
	type(particle_samples_arr_type)::sams_emri, sams_emri_by, sams_emri_sg, sams_plunge
    type(event_sets_comp)::eveset
	integer i, isnap
    logical,external::selection_2gene_single,selection_within_bd
	logical,external::selection_emri_sgsource
	!logical,external::selection_fgw_gw, selection_fgw_kl, selection_td
	!print*, "sa%n=",sa%n

	call get_sams_sts_single(sa, sams_sel,  -1, eveset%etot, isnap)
	!print*, "sams_sel%n=", sams_sel%n
	call get_sams_sts_single(sa, sams_sel_norm,  exit_normal, eveset%enorm, isnap)	

    call get_sams_sts_condition_single(sams_sel_norm, sams_sel_norm_inbd,eveset%enorm_withinbd, &
        isnap, selection_within_bd)
    
	!call get_avgmass_single(sams_sel_norm,eveset%enorm%avg_bhmass(isnap))
	!print*, "3"
	call get_sams_event_rate(eveset%etot, eveset%enorm, eveset%etot, isnap)
	!print*, "4"
	call get_sams_event_rate(eveset%enorm, eveset%enorm, eveset%etot, isnap)
	
	call get_sams_event_rate(eveset%egw, eveset%enorm, eveset%etot, isnap)
    

    call get_sams_sts_single(sa, sams_plunge,  exit_plunge_single, eveset%egw_plunge, isnap)
	call get_sams_event_rate(eveset%egw_plunge, eveset%enorm, eveset%etot, isnap)
	print*, "sams_plunge%n=",sams_plunge%n
    !call get_sams_sts(sa, sams_sel,  exit_max_reach, eveset%emaxsteps,  isnap)
	!call get_sams_sts(sa, sams_sel,  exit_boundary_min, eveset%eemin, isnap)
	!print*, "start eemax"
	call get_sams_sts_single(sa, sams_sel,  exit_boundary_max, eveset%eemax, isnap)
	!call get_sams_event_rate(eveset%eemin, eveset%enorm, eveset%etot, isnap)
	call get_sams_event_rate(eveset%eemax, eveset%enorm, eveset%etot, isnap)
	!print*, "end eemax"
	!call get_sams_sts(sa, sams_sel,  exit_invtransit, eveset%einvtrans, isnap)
	!call get_sams_sts(sa, sams_sel, exit_other, eveset%eother, isnap)
	!call get_sams_sts(sa, sams_sel, exit_gw_bhb_mbh_enc_plunge, eveset%egw_plunge, isnap)
	
end subroutine

subroutine get_sams_events_star(sa, sams_sel_norm, sams_sel_norm_inbd,sams_sel_td, eveset,isnap)
	use com_main_gw
	implicit none
	type(particle_samples_arr_type)::sa, sams_sel, sams_sel_norm, sams_sel_norm_inbd
	type(particle_samples_arr_type)::sams_sel_td
    type(event_sets_star)::eveset
	integer i, isnap
	logical,external:: selection_within_bd, selection_td
	!logical,external::selection_fgw_gw, selection_fgw_kl, selection_td
	!print*, "sa%n=",sa%n
	call get_sams_sts_single(sa, sams_sel,  -1, eveset%etot, isnap)
	!print*, "sams_sel%n=", sams_sel%n
	call get_sams_sts_single(sa, sams_sel_norm,  exit_normal, eveset%enorm, isnap)	
    
    call get_sams_sts_condition_single(sams_sel_norm, sams_sel_norm_inbd,eveset%enorm_withinbd, &
        isnap, selection_within_bd)
	!call get_avgmass_single(sams_sel_norm,eveset%enorm%avg_bhmass(isnap))
	!print*, "3"
	call get_sams_event_rate(eveset%etot, eveset%enorm, eveset%etot, isnap)
	!print*, "4"
	call get_sams_event_rate(eveset%enorm, eveset%enorm, eveset%etot, isnap)
	if(ctl%include_loss_cone.ge.1)then
		call get_sams_sts_condition_single(sa, sams_sel_td,  eveset%etd, isnap,selection_td)
        !print*, "td full"
        !ctl%debug=2
		call get_sams_sts_single(sa, sams_sel,  exit_tidal_full , eveset%etdfull,  isnap)
        !print*, "td empty"
		call get_sams_sts_single(sa, sams_sel,  exit_tidal_empty, eveset%etdempty, isnap)
        !ctl%debug=0
		call get_sams_event_rate(eveset%etd, eveset%enorm, eveset%etot,  isnap)
		call get_sams_event_rate(eveset%etdfull, eveset%enorm, eveset%etot, isnap)
		call get_sams_event_rate(eveset%etdempty,eveset%enorm, eveset%etot, isnap)
	end if
	call get_sams_sts_single(sa, sams_sel,  exit_boundary_max, eveset%eemax, isnap)
	call get_sams_event_rate(eveset%eemax, eveset%enorm, eveset%etot, isnap)
	
end subroutine


logical function selection_within_bd(sp)
    use com_main_gw
    implicit none
    type(particle_sample_type)::sp
    if(sp%en<ctl%energy_boundary)then
        selection_within_bd=.true.
    else
        selection_within_bd=.false.
    end if
end function

logical function selection_emri_sgsource(sp)
	use com_main_gw
	implicit none
	type(particle_sample_type)::sp
	integer flag

	if(sp%exit_flag.eq.exit_boundary_max)then
		selection_emri_sgsource=.true.
	else
		selection_emri_sgsource=.false.
	end if
end function

subroutine get_sams_sts_single(sa, sams_sel, flag, eve, isnap)
	use com_main_gw
	implicit none
	integer flag, isnap, i
	type(particle_samples_arr_type)::sa, sams_sel
	type(event_rate_type):: eve
    real(8) wn

	!print*, "1"
	call sa%select(sams_sel,flag,-1d0,-1d0)
	!print*, "2:",sams_sel%n
	!call sams_get_weight(sams_sel,ctl%clone_e0)
    !print*, "sams_sel%n=",sams_sel%n
    !wn=0
	!do i=1, sams_sel%n
        !wn=wn+sams_sel%sp(i)%weight_real
        !print*,  i, wn, sams_sel%sp(i)%weight_real, sams_sel%sp(i)%weight_n, &
        !    sams_sel%sp(i)%weight_clone, sams_sel%sp(i)%obtype
        !if(mod(i,1000).eq.0)then
        !    read(*,*)
        !end if
    !end do
    !read(*,*)
	call sams_get_eve_num(sams_sel%sp(1:sams_sel%n)%weight_real,sams_sel%n, &
			 eve, isnap)
    !if(ctl%debug>1)then
    !    print*,dms%weight_asym
    !    do i=1, sams_sel%n
    !        print*, sams_sel%sp(i)%weight_real, sams_sel%sp(i)%weight_clone
    !    end do
    !    print*, eve%ntot_eve_simu_w(isnap), eve%ntot_eve_simu(isnap)
    !endif

end subroutine


subroutine get_sams_sts_condition_single(sa, sams_sel, eve, isnap,sub_condition)
	use com_main_gw
	implicit none
	integer flag, isnap
	type(particle_samples_arr_type)::sa, sams_sel
	type(event_rate_type):: eve
	logical, external::sub_condition

	call sams_arr_select_condition_single(sa, sams_sel,sub_condition)
	!call sams_get_weight(sams_sel,ctl%clone_e0)
	call sams_get_eve_num(sams_sel%sp(1:sams_sel%n)%weight_real,sams_sel%n, &
	eve, isnap)
end subroutine

subroutine get_avgmass_single(sa,avg_mass)
	use com_main_gw
	implicit none
	type(particle_samples_arr_type)::sa
	real(8) avg_mass
	integer i

	avg_mass=0
	do i=1,sa%n
		avg_mass=avg_mass+sa%sp(i)%m
	end do
	avg_mass=avg_mass/dble(sa%n)
end subroutine

logical function selection_td(sp)
	use com_main_gw
	implicit none
	type(particle_sample_type)::sp
	integer flag
	flag=sp%exit_flag

	selection_td=.false.
	if(flag.eq.exit_tidal_full.or.flag.eq.exit_tidal_empty.and.sp%m>0.1d0)then
		selection_td=.true.
	end if
end function