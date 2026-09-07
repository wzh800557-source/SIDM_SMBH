subroutine RR_mpi(ch,total_time)
use com_main_gw
! timeunit=0.158yr
IMPLICIT NONE
	integer n
	integer steps,i
	external::fl_rcnum
    integer,save::try=1
	!integer,parameter:: chunk=4
	logical :: flag_normal
	type(chain_type)::ch
	type(chain_pointer_type),pointer::pt,ps
	real(8) t1,t2,total_time
	interface 
	subroutine run_one_sample_particle(pt, run_time)
		use com_main_gw
		implicit none
		type(chain_pointer_type),target::pt
		real(8) run_time
	end subroutine
	end interface
	!integer,parameter::imax_trace=10000
	!integer i_trace

    if(ctl%chattery.ge.1)then
		write(chattery_out_unit,fmt=*) 'proc',rid,"starting..."
    end if
	if(rid.eq.mpi_master_id)then
		write(chattery_out_unit,*) "simu begin"
		if(ctl%chattery.ge.1) write(chattery_out_unit,fmt=*) "total number of samples:", bksams%n
		if(ctl%chattery.ge.1) write(chattery_out_unit,fmt=*) "total time:", total_time,"Myr"
		if(ctl%chattery.ge.1) write(chattery_out_unit,fmt=*) "total number of procs:", ctl%ntasks
	end if
    if(ctl%chattery.ge.1)then
		write(chattery_out_unit,fmt="(A25, A5, 10A7)") "-------result","type","idx", "nhiar","cpuid","eid", "ngene"
	endif
	

	dc_grid_xstep=dms%dc0%s2_dee%xstep
	dc_grid_ystep=dms%dc0%s2_dee%ystep
	
	!ctl%bin_mass_Nbalance=0
    !ctl%bin_mass_Nbalance_ot=0
    pt=>ch%head
    !if(rid.eq.7)then
    !	print*, "rid,i, ac,ec=",rid,i, sp(449)%ob%ac, sp(449)%ob%ec
    !end if
    !call pt%bg%output()
    !call bksams%output_screen(2)

loop1:	do while(associated(pt))        
		if(pt%ob%exit_flag.eq.exit_normal)then
			select case (ctl%trace_all_sample)
			case(-1)
				if(pt%ob%obtype.eq.star_type_MS)then
					pt%ob%write_down_track=record_track_detail
				else
					pt%ob%write_down_track=0
				end if
			case(-2)
				if(pt%ob%obtype.eq.star_type_BH)then
					pt%ob%write_down_track=record_track_detail
				else
					pt%ob%write_down_track=0
				end if
				
			end select
			
			call run_one_sample(pt,total_time)
			
		end if
		
        ps=>pt
        pt=>pt%next
        !print*, "1"
        if(ctl%clone_scheme.ge.1)then
            if(ps%ob%exit_flag.eq.exit_invtransit &
                .and.ctl%del_cross_clone.ge.1.or.(ps%ob%exit_flag.eq.exit_boundary_min))then
				!print*, "associated(ps)?",associated(ps), ps%ob%exit_flag
                !call delete_clone_particle(pt%prev)
                if(associated(ps%prev))then
                    if(allocated(ps%ob))then
                        deallocate(ps%ob)
                    end if
					!call bksams%output_screen()
                    call chain_pointer_delete_item_chain_type(ps)
					!stop
                else
                    if(.not. associated(pt)) exit loop1
                    !print*, "delete head node cycle"
                    !cycle 
                    !call ch%output_screen(2,10)
                    call destroy_attach_pointer_chain_type(ch%head)
                    ch%head=>pt
                    call pt%set_head()                    
                    !print*, "1"
                    pt%prev=>null()
                    
					
                    
                end if
                cycle
            end if
        end if
        !call check("after delete in RRMPI")
    end do loop1
end subroutine

