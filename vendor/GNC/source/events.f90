module md_events
	implicit none
	type event_rate_type
		real(8),allocatable:: ntot_eve_simu_w(:)  ! weighted total number of eve events in (one snapshot of a) simulation
		!real(8),allocatable:: ntot_simu_w(:)    ! weighted total number of normal BBH in (one snapshot of a) simulation
		real(8),allocatable:: ntot_eve_simu(:)  ! total number of eve events in (one snapshot of a) simulation
		!real(8),allocatable:: ntot_simu(:)      ! total number of normal BBH in (one snapshot of a) simulation
		!integer nby_assume    ! total number of binary BH assumed at a given time in reality
		real(8),allocatable:: eve_rate(:)       
		integer nsp
		real(8),allocatable:: t_snap(:)         ! the time of the snapshot
		real(8),allocatable:: p_eve(:)
		real(8) dt
		real(8),allocatable:: numavg_flyby(:),numavg_kl(:), numavg_encnt(:) ! average number of flybys or KLs experienced per Myr
		contains
		 	procedure write_sample
		  	procedure read_sample
		  	generic :: write(unformatted) => write_sample
		  	generic :: read(unformatted) => read_sample
	end type
    type event_sets_star
        integer nsp
        real(8) tot_time
        type(event_rate_type)::etot,enorm, enorm_withinbd
        type(event_rate_type)::etd, etdfull,etdempty
        type(event_rate_type)::eemin,eemax, emaxsteps
        contains
        procedure::init=>init_event_sets_star
    end type
	type,extends(event_sets_star)::event_sets_comp
		type(event_rate_type)::egw
		type(event_rate_type)::egw_emri, eother, egw_plunge
        contains
        procedure::init=>init_event_sets_comp
	end type
	
	real(8) :: frac_nbh, frac_mbh, frac_rate_supply
	!type(event_sets_bbh)::pteve_bbh
	type(event_sets_star)::pteve_star
    type(event_sets_comp)::pteve_sbh,pteve_wd, pteve_ns, pteve_bd
    private::init_event_sets_star,init_event_sets_comp
contains
	subroutine init_event(eve, tspan, nsp)
		implicit none
		integer nsp,i
		type(event_rate_type):: eve
		real(8) tspan
		eve%nsp=nsp
		!eve%nby_assume=num0
		if(allocated(eve%ntot_eve_simu))then
			deallocate(eve%ntot_eve_simu_w, eve%t_snap, eve%p_eve,eve%ntot_eve_simu, &
                       eve%eve_rate, eve%numavg_flyby, eve%numavg_kl, eve%numavg_encnt)
		end if
		allocate(eve%ntot_eve_simu_w(nsp),eve%t_snap(nsp), &
				 eve%p_eve(nsp), eve%ntot_eve_simu(nsp), &
                 eve%eve_rate(nsp), eve%numavg_flyby(nsp), eve%numavg_kl(nsp), &
				  eve%numavg_encnt(nsp) )
		do i=1, nsp
			eve%t_snap(i)=tspan/dble(nsp)*dble(i)
		end do
!		eve%eve_rate=0; 
		eve%ntot_eve_simu_w=0; eve%dt=dble(tspan)/dble(nsp)
		eve%ntot_eve_simu=0; 
	end subroutine

	subroutine print_events(et,enm,j)
		implicit none
		integer j
		type(event_rate_type)::et,enm
		print*, "isnap=",j
		print*, "Ntot, Ntotw=", et%ntot_eve_simu(j), et%ntot_eve_simu_w(j)
		print*, "Nnorm, Nnormw=", enm%ntot_eve_simu(j),enm%ntot_eve_simu_w(j)
	end subroutine
	subroutine output_eveset_avgscatter_rate(eveset, fn, frac)
		use fun_math
		implicit none
		class(event_sets_star)::eveset
		character*(*) fn
		integer i
		real(8) frac
		integer nbg, ned, nar

		nbg=dble(eveset%nsp)*frac
		ned=eveset%nsp
		nar=ned-nbg+1
		open(unit=999,file=trim(adjustl(fn))//"_avgrates.txt")
		select type (eveset)
		type is(event_sets_star)
			write(unit=999,fmt="(20A25)") "Rate_td(Gyr-1)", "S_td", "R_tdfull", "S_tdfull",&
            "R_tdempty","s_tmpfull", "R_emax", "s_emax"
			
			write(unit=999,fmt="(10F25.10)")  farravg(eveset%etd%eve_rate(nbg:ned),nar), &
				farrsct(eveset%etd%eve_rate(nbg:ned),nar), farravg(eveset%etdfull%eve_rate(nbg:ned),nar),&
				farrsct(eveset%etdfull%eve_rate(nbg:ned),nar), farravg(eveset%etdempty%eve_rate(nbg:ned),nar),&
				farrsct(eveset%etdempty%eve_rate(nbg:ned),nar), farravg(eveset%eemax%eve_rate(nbg:ned),nar),&
				farrsct(eveset%eemax%eve_rate(nbg:ned),nar)

		type is(event_sets_comp)
			write(unit=999,fmt="(20A25)") "Rate_plunge(Gyr-1)", "S_plunge"
			
			write(unit=999,fmt="(10F25.10)")  farravg(eveset%egw_plunge%eve_rate(nbg:ned), nar),&
			farrsct(eveset%egw_plunge%eve_rate(nbg:ned), nar)
		end select
		close(999)
	end subroutine
	subroutine output_eveset_rate(eveset,fn)
		implicit none
		class(event_sets_star)::eveset
		character*(*) fn
		integer i

		open(unit=999,file=trim(adjustl(fn))//"_rates.txt")
		select type (eveset)
		type is(event_sets_star)
			write(unit=999,fmt="(20A25)") "Tsnap(Myr)","Rate_td(Gyr-1)", "R_tdfull", &
            "R_tmpfull", "R_emax"
			do i=1, eveset%nsp
				write(unit=999,fmt="(10F20.10)")  eveset%etot%t_snap(i), eveset%etd%eve_rate(i), &
					eveset%etdfull%eve_rate(i), eveset%etdempty%eve_rate(i), eveset%eemax%eve_rate(i)
			end do
		type is(event_sets_comp)
			write(unit=999,fmt="(20A20)") "Tsnap(Myr)","Rate_td(Gyr-1)", "R_tdfull", &
			"R_tmpfull", "R_emax", "R_plunge"
			do i=1, eveset%nsp
				write(unit=999,fmt="(10F20.10)")  eveset%etot%t_snap(i), eveset%etd%eve_rate(i), &
					eveset%etdfull%eve_rate(i), eveset%etdempty%eve_rate(i), eveset%eemax%eve_rate(i), &
					eveset%egw_plunge%eve_rate(i)
			end do
		
		end select
		close(999)
	end subroutine

     ! Unformatted writing for the sample derived type
     subroutine write_sample(eve, unit, iostat, iomsg)
       class(event_rate_type), intent(in) :: eve
       integer, intent(in) :: unit
       integer, intent(out) :: iostat
       character(*), intent(inout) :: iomsg

       integer i

       ! Write a record giving sizes for the allocation
       write(unit, iostat=iostat, iomsg=iomsg) eve%ntot_eve_simu_w, eve%t_snap, eve%ntot_eve_simu, &
                       eve%eve_rate,eve%p_eve, eve%numavg_flyby, eve%numavg_kl, eve%numavg_encnt

     end subroutine write_sample

     ! Unformatted reading for the sample derived type
     subroutine read_sample(eve, unit, iostat, iomsg)
       class(event_rate_type), intent(inout) :: eve
       integer, intent(in) :: unit
       integer, intent(out) :: iostat
       character(*), intent(inout) :: iomsg
       integer i
       integer sizeb, sizec

       ! We first have a record telling us the sizes of components
       read(unit, iostat=iostat, iomsg=iomsg) eve%ntot_eve_simu_w, eve%t_snap, eve%ntot_eve_simu, &
	   eve%eve_rate,eve%p_eve, eve%numavg_flyby, eve%numavg_kl, eve%numavg_encnt

     end subroutine read_sample
     subroutine init_event_sets_star(eveset, tot_time, nsap)
		implicit none
		class(event_sets_star)::eveset
		real(8) tot_time
		integer nsap
		call init_event(eveset%etot, tot_time, nsap)
		call init_event(eveset%enorm, tot_time, nsap)
        call init_event(eveset%enorm_withinbd, tot_time, nsap)        
		call init_event(eveset%etd, tot_time, nsap)
		call init_event(eveset%etdfull, tot_time, nsap)
		call init_event(eveset%etdempty, tot_time, nsap)
		call init_event(eveset%emaxsteps, tot_time, nsap)
		call init_event(eveset%eemin, tot_time, nsap)
		call init_event(eveset%eemax, tot_time, nsap)
		eveset%tot_time=tot_time
		eveset%nsp=nsap
	end subroutine
	 subroutine init_event_sets_comp(eveset, tot_time, nsap)
		implicit none
		class(event_sets_comp)::eveset
		real(8) tot_time
		integer nsap
        call eveset%event_sets_star%init(tot_time, nsap)
		call init_event(eveset%egw, tot_time, nsap)
		call init_event(eveset%egw_emri, tot_time, nsap)
        call init_event(eveset%egw_plunge, tot_time, nsap)
		call init_event(eveset%eother, tot_time, nsap)
		call init_event(eveset%egw_plunge, tot_time, nsap)

	end subroutine
	
    subroutine output_eveset_N(eveset, fout, flag)
        implicit none
        class(event_sets_star)::eveset
        character*(*) fout
        integer i, flag

        open(unit=999,file=trim(adjustl(fout))//"_event_Nweight.txt")
        select type (eveset)
        type is(event_sets_star)
			select case(flag)
			case(1)
				write(unit=999,fmt="(20A25)") "Tsnap","N_all","N_norm", "N_norm_bd", "N_emax"
				do i=1, eveset%nsp
					write(unit=999,fmt="(20F25.10)") eveset%etot%t_snap(i), eveset%etot%ntot_eve_simu_w(i), &
					eveset%enorm%ntot_eve_simu_w(i), eveset%enorm_withinbd%ntot_eve_simu_w(i), &
					eveset%eemax%ntot_eve_simu_w(i)
				end do
			case(2)
				write(unit=999,fmt="(20A25)") "Tsnap","N_all","N_norm", "N_norm_bd", "N_td", "N_tdfull", &
				"N_tmpfull", "N_emax"
				do i=1, eveset%nsp
					write(unit=999,fmt="(20F25.10)") eveset%etot%t_snap(i), eveset%etot%ntot_eve_simu_w(i), &
					eveset%enorm%ntot_eve_simu_w(i), eveset%enorm_withinbd%ntot_eve_simu_w(i), &
					eveset%etd%ntot_eve_simu_w(i), eveset%etdfull%ntot_eve_simu_w(i), &
					eveset%etdempty%ntot_eve_simu_w(i), eveset%eemax%ntot_eve_simu_w(i)
				end do
			end select
        type is(event_sets_comp)
			select case(flag)
			case(1)
				write(unit=999,fmt="(20A25)") "Tsnap","N_all","N_norm", "N_norm_bd", "N_emax", "N_plunge"
				do i=1, eveset%nsp
					write(unit=999,fmt="(20F25.10)") eveset%etot%t_snap(i), eveset%etot%ntot_eve_simu_w(i), &
					eveset%enorm%ntot_eve_simu_w(i), eveset%enorm_withinbd%ntot_eve_simu_w(i), &
					eveset%eemax%ntot_eve_simu_w(i), &
					eveset%egw_plunge%ntot_eve_simu_w(i)
					
				end do
			case(2)
					write(unit=999,fmt="(20A25)") "Tsnap","N_all","N_norm", "N_norm_bd", "N_td", "N_tdfull", &
				"N_tmpfull", "N_emax", "N_plunge"
				do i=1, eveset%nsp
					write(unit=999,fmt="(20F25.10)") eveset%etot%t_snap(i), eveset%etot%ntot_eve_simu_w(i), &
					eveset%enorm%ntot_eve_simu_w(i), eveset%enorm_withinbd%ntot_eve_simu_w(i), &
					eveset%etd%ntot_eve_simu_w(i), eveset%etdfull%ntot_eve_simu_w(i), &
					eveset%etdempty%ntot_eve_simu_w(i), eveset%eemax%ntot_eve_simu_w(i), &
					eveset%egw_plunge%ntot_eve_simu_w(i)
					
				end do
			
			end select
         end select
        close(999)
        

        open(unit=999,file=trim(adjustl(fout))//"_event_N.txt")
        select type (eveset)
        type is(event_sets_star)
		 	select case(flag)
			case(1)
            write(unit=999,fmt="(20A25)") "Tsnap","N_all","N_norm", "N_norm_bd", "N_emax"
            do i=1, eveset%nsp
                write(unit=999,fmt="(20F25.10)") eveset%etot%t_snap(i), eveset%etot%ntot_eve_simu(i), &
                eveset%enorm%ntot_eve_simu(i), eveset%enorm_withinbd%ntot_eve_simu(i), &
                eveset%eemax%ntot_eve_simu(i)
            end do
			case(2)
				write(unit=999,fmt="(20A25)") "Tsnap","N_all","N_norm", "N_norm_bd", "N_td", "N_tdfull", &
				"N_tdempty", "N_emax"
				do i=1, eveset%nsp
					write(unit=999,fmt="(20F25.10)") eveset%etot%t_snap(i), eveset%etot%ntot_eve_simu(i), &
					eveset%enorm%ntot_eve_simu(i), eveset%enorm_withinbd%ntot_eve_simu(i), &
					eveset%etd%ntot_eve_simu(i), eveset%etdfull%ntot_eve_simu(i), &
					eveset%etdempty%ntot_eve_simu(i), eveset%eemax%ntot_eve_simu(i)
				end do
			end select
        type is(event_sets_comp)
			select case(flag)
			case(1)
            write(unit=999,fmt="(20A25)") "Tsnap","N_all","N_norm", "N_norm_bd", "N_emax", "N_plunge"
            do i=1, eveset%nsp
                write(unit=999,fmt="(20F25.10)") eveset%etot%t_snap(i), eveset%etot%ntot_eve_simu(i), &
                eveset%enorm%ntot_eve_simu(i), eveset%enorm_withinbd%ntot_eve_simu(i), &
                eveset%eemax%ntot_eve_simu(i), eveset%egw_plunge%ntot_eve_simu(i)
            end do
			case(2)
				write(unit=999,fmt="(20A25)") "Tsnap","N_all","N_norm", "N_norm_bd", "N_td", "N_tdfull", &
				"N_tdempty", "N_emax", "N_plunge"
				do i=1, eveset%nsp
					write(unit=999,fmt="(20F25.10)") eveset%etot%t_snap(i), eveset%etot%ntot_eve_simu(i), &
					eveset%enorm%ntot_eve_simu(i), eveset%enorm_withinbd%ntot_eve_simu(i), &
					eveset%etd%ntot_eve_simu(i), eveset%etdfull%ntot_eve_simu(i), &
					eveset%etdempty%ntot_eve_simu(i), eveset%eemax%ntot_eve_simu(i), &
					eveset%egw_plunge%ntot_eve_simu(i)
				end do
			end select
         end select
        close(999)
    end subroutine
    subroutine output_eveset_txt(eveset, fout,flag)
        implicit none
        class(event_sets_star)::eveset
        character*(*) fout
		integer flag

!        call output_events_single_bin(pteve, trim(adjustl(fout))//"events")
        call output_eveset_N(eveset, fout,flag)
		if(flag>1)then
        	call output_eveset_rate(eveset, fout)
			call output_eveset_avgscatter_rate(eveset,fout, 0.5d0)
		end if
    end subroutine
!	subroutine input_events_single_bin(eveset,fl)
!		implicit none
!		character*(*) fl
!		integer i,n
!		type(event_sets)::eveset
!	
!		open(unit=999,file=trim(adjustl(fl))//".eve", form='unformatted',access='stream',status='old')
!		read(unit=999) eveset%tot_time,  eveset%nsp
!		call init_sams_eventsets(eveset, eveset%tot_time,  eveset%nsp)
!		read(unit=999) eveset%etot, eveset%enorm, eveset%etd, eveset%etdfull, eveset%etdempty,&
!		eveset%egw, eveset%emaxsteps, eveset%eemin, eveset%eemax, eveset%einvtrans
!						!print*, egw_bhb_mbh%norm_rate
!	!					print*, enorm%norm_rate
!		close(unit=999)
!	end subroutine
!	subroutine  output_events_single_bin(eveset,fl)
!		implicit none
!		character*(*) fl
!		integer i,n
!		type(event_sets)::eveset
!	
!		open(unit=999,file=trim(adjustl(fl))//".eve", form='unformatted',access='stream')
!		write(unit=999) eveset%tot_time,  eveset%nsp
!		write(unit=999) eveset%etot, eveset%enorm, eveset%etd, eveset%etdfull, eveset%etdempty,&
!		eveset%egw, eveset%emaxsteps, eveset%eemin, eveset%eemax, eveset%einvtrans
!		close(unit=999)
!	end subroutine

!	subroutine input_events_bin(eveset,fl)
!		implicit none
!		character*(*) fl
!		integer i,n
!		type(by_event_sets)::eveset
!	
!		open(unit=999,file=trim(adjustl(fl))//".eve", form='unformatted',access='stream',status='old')
!		read(unit=999) eveset%tot_time,  eveset%nsp
!		call init_sams_byeventsets(eveset, eveset%tot_time,  eveset%nsp)
!		read(unit=999) eveset%etot, eveset%enorm, eveset%etd, eveset%etdfull, eveset%etdempty,&
!		eveset%egw, eveset%egw_kl, eveset%egw_iso, eveset%ebyion, eveset%emaxsteps, &
!		eveset%eemin, eveset%eemax, eveset%einvtrans, eveset%ebyexchange, eveset%egw_fgw_gw,&
!		eveset%egw_fgw_kl, eveset%egw_bhb_mbh,eveset%egw_embhb, eveset%egw_bkenc
!						!print*, egw_bhb_mbh%norm_rate
!	!					print*, enorm%norm_rate
!		close(unit=999)
!	end subroutine
!	subroutine  output_events_bin(eveset,fl)
!		implicit none
!		character*(*) fl
!		integer i,n
!		type(by_event_sets)::eveset
!	
!		open(unit=999,file=trim(adjustl(fl))//".eve", form='unformatted',access='stream')
!		write(unit=999) eveset%tot_time,  eveset%nsp
!		
!		write(unit=999) eveset%etot, eveset%enorm, eveset%etd, eveset%etdfull, eveset%etdempty,&
!		eveset%egw, eveset%egw_kl, eveset%egw_iso, eveset%ebyion, eveset%emaxsteps, &
!		eveset%eemin, eveset%eemax, eveset%einvtrans, eveset%ebyexchange, eveset%egw_fgw_gw,&
!		eveset%egw_fgw_kl, eveset%egw_bhb_mbh,eveset%egw_embhb, eveset%egw_bkenc
!		close(unit=999)
!	end subroutine
end module



