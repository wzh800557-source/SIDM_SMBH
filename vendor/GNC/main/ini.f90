program ini
	use com_main_gw
	implicit none
	integer i
	character*(4) tmpi
	logical ex
	type(particle_samples_arr_type),allocatable::smsa(:)
	type(sts_fc_type) ::fc
	integer ierr
    real(8) t1,t2
	!type(particle_samples_arr_type)::sps_arr

    call readin_model_par("model.in")
   
	!if(ctl%include_binary.ge.1)then
    !!	call readin_binary_par("by.in")
	!end if
	!call read_output_in("output.in")
	call init_mpi()
	call init_model()	

    if(rid.eq.mpi_master_id)then
        print*, "mpi_master_id=", mpi_master_id
    end if

	call cfs%input_bin(trim(adjustl(ctl%cfs_file_dir)))
    print*,  "readin cfs finished", cfs%nj, cfs%ns, rid
	
    if(rid.eq.0)then
	    call print_model_par()
        call print_current_code_version()
        call cpu_time(t1)
    end if

	allocate(smsa(ctl%ntasks))
	
	!print*, ctl%ntasks
	!if(ctl%single_bk_evl_mode.ge.1)then
	print*, "proc rid start", rid
	call set_seed(ctl%same_rseed_ini, ctl%seed_value+rid)

	write(unit=tmpi,fmt="(I4)") rid+ctl%ntask_bg+1
	!print*, "start set_bkchain_sample"
	!print*, "bksams:", associated(bksams%sp)
	!print*, "bksams arr:", associated(bksams_arr%sp)
	!print*, "sps arr:" , associated(sps_arr%sp),size(sps_arr%sp)
    print*, "start dms init"
	call set_dm_init(dms)
      call load_df_from_table(dms, "f_evolved.txt")
	!print*, "get_ge_by_root",rid

	call get_init_samples(bksams_arr_ini)
    do i=1, bksams_arr_ini%n
        if(bksams_arr_ini%sp(i)%exit_flag.ne.exit_normal)then
            print*, "i=",i,bksams_arr_ini%sp(i)%exit_flag
            stop
        end if
    end do

	call set_chain_samples(bksams, bksams_arr_ini)
    !stop
	!call all_chain_to_arr_single(Bksams,smsa(i))
	!print*, "start collect"
	dms%weight_asym=1d0
	call collection_data_single_bympi(smsa, ctl%ntasks)
	call mpi_barrier(mpi_comm_world,ierr)
	
    call get_ge_by_root(smsa, ctl%ntasks,.true.)
   ! print*, "update_weights"
    call update_weights()
    !print*, "get_dms"
	call get_dms(dms)
    if(rid.eq.0)then
        call dms%print_norm(6)
        call print_num_boundary(dms)
		call print_num_all(dms)
    end if
    !print*, "rid=", rid, 10**dms%dc0%s2_de_0%xcenter
    !call mpi_barrier(mpi_comm_world,ierr)
    !print*, "rid=", rid, dms%dc0%s2_de_0%ycenter
    !call get_num_at_boundary(bksams, bysams)
    !call print_bin_mass_N()
	call MPI_barrier(MPI_COMM_WORLD, ierr)
	!print*, "start set weights"
	!call update_weights()
	
	call bksams%output_bin("output/ini/bin/single/samchn"//trim(adjustl(tmpi)))

	if(rid.eq.0)then
		!select case(ctl%intaskmode)
		!case(task_mode_new)
#ifdef HDF5			
			call output_dms_hdf5_pdf(dms, "output/ini/hdf5/dms_0_0")
#else
!			call output_all_barge_txt(dms,"output/ini/hdf5/dms_0_0")
#endif
			!call output_coef_ba16_txt(dms%dc0, "output/ini/txt/dc")
			call dms%output_bin("output/ini/bin/dms.bin")
		!case(task_mode_append)
!ifdef HDF5			
!			call output_dms_hdf5_pdf(dms, "output/ini/hdf5/dms_0_0_append")
!else
!			call output_all_barge_txt(dms,"output/ini/hdf5/dms_0_0_append")
!endif	
!			!call output_coef_ba16_txt(dms%dc0, "output/ini/txt/dc")
!			call dms%output_bin("output/ini/bin/dms_append.bin")
!		case default
!			print*, "error, define taks mode"
!			stop
!		end select
		print*, "output control..."
		!call output_control("output/")
		!print*, "output encounter rate"
		!call output_encounter_rate("output/time")
		write(*,*) "init finished"	
        call cpu_time(t2)
        print*, "total_time=",t2-t1
	end if
	!call print_minj()
	call stop_mpi()
	
end