program main
    use com_main_gw
!	use md_coeff
	implicit none
	integer i,j,ierror
	real(8) t1,t2
	character*(4) tmprid, tmpssnapid, tmpssnapid_prev, tmpsupdate_max
	logical ex
    INTEGER*4 getpid, pid
	character*(20) pid_root_str, pid_cld_str

    call readin_model_par("model.in")
	
	!allocate(smsa(ctl%ntask_total))
	call init_mpi()
    call init_model()
    if(rid.eq.0)then
	    call print_model_par()
        call print_current_code_version()
    end if
    !rid=2
	if(rid.eq.0)then
		write(*,*) 'root process id = ', proc_id
	end if
	if(rid.eq.1)then
		write(*,*) 'cld process id = ', proc_id
	end if
    call cfs%input_bin(trim(adjustl(ctl%cfs_file_dir)))
    print*, rid, "readin cfs finished", cfs%nj, cfs%ns
        
	mpi_master_id=0
	!print*, "use collision method:", use_collision_method
		
	write(unit=tmprid,fmt="(I4)") rid+1+ctl%ntask_bg

	call prepare_ini_data(tmprid)
	
	if(rid.eq.0)then
		call cpu_time(t1)
	end if
	call set_seed(ctl%same_rseed_evl, ctl%seed_value+rid)
    call update_arrays_single()
    call set_clone_weight(bksams)
    call set_real_weight(bksams)

	do i=1, ctl%n_spshot_total
		write(unit=tmpssnapid,fmt="(I4)") i
		call run_snapshots(tmprid, tmpssnapid, pid_root_str, pid_cld_str, i)
		!if(ctl%include_binary.ge.1)then
		!	call run_binary(tmprid, tmpssnapid, i)	
		!end if
	end do

	if(rid.eq.0)then
		call cpu_time(t2)
		print*, "total time:", t2-t1
	end if
	print*, "finished main", rid
!	call mpi_barrier(mpi_comm_world,ierror)
	call stop_mpi()
!	print*, i_get_nu_m
	call deallocate_chains_arrs()
end

subroutine run_snapshots(tmprid, tmpssnapid, pid_root_str, pid_cld_str,i)
	use com_main_gw
	implicit none
	integer i, j, k
	character*(*) tmprid, tmpssnapid, pid_cld_str,pid_root_str
	character*(4) tmpj
	character*(100) str_, str_update
	type(particle_samples_arr_type)::smstot,smstot_sel
	real(8) cur_time_f,cur_time_i
	type(sts_fc_type) ::fc
	type(particle_samples_arr_type),allocatable::smsa(:)
	logical update_dms
	
	allocate(smsa(ctl%ntasks))

    ctl%burn_in_phase=.false.
	str_=trim(adjustl(tmprid))//"_"//trim(adjustl(tmpssnapid))

	do j=1, ctl%num_update_per_snap
		write(unit=tmpj, fmt="(I4)") j
		!print*, sms%sp(1)%ob%write_down_track
		!read(*,*)
        cur_time_i=ctl%ts_spshot*dble(i-1)+ctl%update_dt*dble(j-1)
		cur_time_f=ctl%update_dt+cur_time_i
        write(chattery_out_unit,fmt=*) "start.. snap, i, j, cur_time_i, f=", i, j, cur_time_i, cur_time_f
		if(ctl%trace_all_sample.ge.1)then
			update_dms=.false.
		else
			update_dms=.true.
		end if
        call run_one_snap(cur_time_i, cur_time_f, smsa, ctl%ntasks,update_dms)
		!call check("end of snapshot")
        write(chattery_out_unit,fmt=*) "finished snapshot i, rid, update_j=", i, rid, j
		if(rid.eq.mpi_master_id)then
			print*, "start output", rid
			str_update=trim(adjustl(str_))//"_"//trim(adjustl(tmpj))
			call dms%output_bin("output/ecev/dms/dms_"//trim(adjustl(tmpssnapid))//"_"//trim(adjustl(tmpj)))
			!select case(ctl%intaskmode)
			!case(task_mode_new)
#ifdef HDF5				
				call output_dms_hdf5_pdf(dms, "output/ecev/dms/dms_"//trim(adjustl(tmpssnapid))//"_"//trim(adjustl(tmpj)))
#else				
!				call output_all_barge_txt(dms,"output/ecev/dms/dms_"//trim(adjustl(tmpssnapid))//"_"//trim(adjustl(tmpj)))
#endif				
!			case(task_mode_append)
!#ifdef HDF5				
!				call output_dms_hdf5_pdf(dms, "output/ecev/dms/dms_"//trim(adjustl(tmpssnapid))//"_"//trim(adjustl(tmpj))//"_append")
!#else
!				call output_all_barge_txt(dms,"output/ecev/dms/dms_"//trim(adjustl(tmpssnapid))//"_"//trim(adjustl(tmpj))//"_append")
!#endif					
!			end select
			print*, "output finished", rid
            !call dms%dc0%s2_de_0%output_txt("output/ecev/txt/s2de0_"//trim(adjustl(tmpssnapid))//"_"//trim(adjustl(tmpj)))
            !call dms%dc0%s2_de_110%output_txt("output/ecev/txt/s2de110_"//trim(adjustl(tmpssnapid))//"_"//trim(adjustl(tmpj)))
		endif
	end do
   
	print*, "start output bins"
	!call update_weights()
    !print*, "finished update weights"
!===================================
	call bksams%output_bin("output/ecev/bin/single/samchn"//trim(adjustl(str_)))	
	if(ctl%trace_all_sample.ge.1)then
		call all_chain_to_arr_single(bksams,bksams_arr)
    	call output_sams_sg_track_txt(Bksams_arr, "output/indvd/")
		
	end if
!	call bksams_arr%output_txt("output/ecev/txt/single/sample"//trim(adjustl(str_)))
!	call bysams_arr%output_txt("output/ecev/txt/by/sample"//trim(adjustl(str_)))
!===================================	
	if (.not.(i.eq.ctl%n_spshot_total)) then
		print*, "create_sams"
        !call check("before convert_sams")
        call convert_sams(bksams)	
	end if

end subroutine
