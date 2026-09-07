subroutine sams_get_eve_num(weights, n, eve, isnap)
	use com_main_gw
	implicit none
	integer n
	real(8) weights(n)
	!type(samples_arr_type)::sps
	type(event_rate_type)::eve
	integer i, isnap
	real(8) dt, neve, neve_w

	neve=0; neve_w=0

	do i=1, n
		neve_w=neve_w+weights(i)
		neve=neve+1
	end do
	eve%ntot_eve_simu_w(isnap)=neve_w/dble(ctl%ntask_total)
	eve%ntot_eve_simu(isnap)=neve/dble(ctl%ntask_total)
end subroutine


subroutine get_sams_event_rate(eve, evenorm, evetot, isnap)
	use com_main_gw
	implicit none
	type(event_rate_type)::eve, evenorm, evetot
	integer isnap
	real(8) norm_w, neve_w
	!norm_w=evenorm%ntot_eve_simu_w(isnap)
	neve_w=eve%ntot_eve_simu_w(isnap)
	eve%eve_rate(isnap)=real(neve_w)/eve%dt*1d3	  ! per Gyr
	!print*, "i,rate,w,dt=",isnap,eve%eve_rate(isnap),neve_w,eve%dt

	if(eve%eve_rate(isnap).ne.0d0)then
		eve%p_eve(isnap)=eve%eve_rate(isnap)/sqrt(neve_w)
	else
		eve%p_eve(isnap)=0d0
	end if
	
end subroutine
subroutine init_sams_events()
	use com_main_gw
	implicit none
	!call init_sams_eventsets(pteve, ctl%total_time, ctl%n_spshot_total)
	!call init_sams_byeventsets(byeve, ctl%total_time, ctl%n_spshot_total)
    call pteve_star%init(ctl%total_time, ctl%n_spshot_total)
    call pteve_sbh%init(ctl%total_time, ctl%n_spshot_total)
	call pteve_wd%init(ctl%total_time, ctl%n_spshot_total)
	call pteve_ns%init(ctl%total_time, ctl%n_spshot_total)
	call pteve_bd%init(ctl%total_time, ctl%n_spshot_total)
	
    !call pteve_bbh%init(ctl%total_time, ctl%n_spshot_total)
	
end subroutine

