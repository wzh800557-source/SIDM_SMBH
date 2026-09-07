
subroutine pro_single(isnap)
    use com_main_gw
    implicit none
    !type(diffuse_mspec)::dms_gene 
    type(particle_samples_arr_type)::bksma_arr_sel
    type(chain_type)::bksma_sel
    integer jsp, isnap

    print*, "pro star"
    call pro_single_star(isnap)
	if(ctl%idxsbh.ne.-1)then
		print*, "pro sbh"
		call pro_single_compact(isnap, "BH", star_type_BH, pteve_sbh)
	end if
	if(ctl%idxns.ne.-1)then
		print*, "pro ns"
		call pro_single_compact(isnap, "NS", star_type_NS, pteve_ns)
	end if
	if(ctl%idxwd.ne.-1)then
		print*, "pro wd"
		call pro_single_compact(isnap, "WD", star_type_WD, pteve_wd)
	end if
	if(ctl%idxbd.ne.-1)then
		print*, "pro bd"
		call pro_single_compact(isnap, "BD", star_type_BD, pteve_bd)
	end if
	
end subroutine

subroutine pro_single_compact(isnap,str_comp, comp_type, eve_sets)
    use com_main_gw
    implicit none
    !type(diffuse_mspec)::dms_gene 
    type(particle_samples_arr_type)::bksma_arr_sel, bksma_arr_sel_norm
    type(particle_samples_arr_type)::bksma_arr_sel_norm_bd
	type(particle_samples_arr_type)::bksma_by_source_emri, bksma_emri
	type(particle_samples_arr_type)::bksma_sg_source_emri, bksma_plunge
	type(event_sets_comp)::eve_sets
	character*(*) str_comp
	integer comp_type
    !type(chain_type)::bksma_sel
    integer jsp, isnap
    character*(4) tmpspid

    call sams_get_weight_clone_single(bksams_arr)
    call set_real_weight_arr_single(bksams_arr)

    call sams_arr_select_type_single(bksams_arr, bksma_arr_sel,comp_type) 
    print*, "arr selection finished"
	if(bksma_arr_sel%n.eq.0) return
    !call chain_select_type_single(bksams, bksma_sel, star_type_BH)
    !print*, "chain selection finished"

    !call get_sts_one_species_single(bksma_sel,bksma_arr_sel, ctl%cc(1), isnap)
    !print*, "get_sts_one_species_single finished"
    call get_sams_events_sbh(bksma_arr_sel, bksma_arr_sel_norm, bksma_arr_sel_norm_bd, &
        bksma_by_source_emri, bksma_sg_source_emri, bksma_emri, bksma_plunge, eve_sets, isnap)
	print*, "get_sams_events_single finshed!"
	call print_events(eve_sets%etot,eve_sets%enorm,isnap)
	print*, "print finished"
	!call get_ge_profile_of_sample_single(sma_arr,isnap, exit_normal, ge_profile(isnap))
	!print*, "sma:chain:"
	!call sma%sp(1)%output()
    
	write(*,fmt="(8A10)") "tot", "tot_sbh", "norm", "norm_inbd"
	write(*,fmt="(8I10)") Bksams_arr%n, bksma_arr_sel%n ,bksma_arr_sel_norm%n, &
		bksma_arr_sel_norm_bd%n
!#ifdef HDF5		
!	call output_sts_star_samples_hdf5(bksma_plunge, "output/pro/"//trim(adjustl(str_comp))//&
!		"/emri_plunge_"//trim(adjustl(tmpspid)))
!#endif    
end subroutine
subroutine pro_single_star(isnap)
    use com_main_gw
    implicit none
    !type(diffuse_mspec)::dms_gene 
    type(particle_samples_arr_type)::bksma_arr_sel, bksma_arr_sel_norm
    type(particle_samples_arr_type)::bksma_arr_sel_norm_bd, bksma_arr_sel_td
    type(chain_type)::bksma_sel
    integer jsp, isnap
	character*(4) tmpspid

    call sams_get_weight_clone_single(bksams_arr)
    call set_real_weight_arr_single(bksams_arr)

   ! print*, "pro star"
    call sams_arr_select_type_single(bksams_arr, bksma_arr_sel,star_type_MS) 
    print*, "arr selection finished"
    !call chain_select_type_single(bksams, bksma_sel, star_type_MS)
    !print*, "chain selection finished"
    call get_sams_events_star(bksma_arr_sel,bksma_arr_sel_norm,&
	 bksma_arr_sel_norm_bd, bksma_arr_sel_td,  pteve_star, isnap)
	print*, "get_sams_events_single finshed!"
	call print_events(pteve_star%etot,pteve_star%enorm,isnap)
	print*, "print finished"
	!call chain_select(bksma_sel,bksma_arr_sel_norm,exit_normal,obj_type=1)
    !call save_one_species_single(ctl%cc(1), isnap, &
    !    "output/pro/MS/")

    write(*,fmt="(8A10)") "tot", "tot_star", "norm", "norm_inbd"
	write(*,fmt="(8I10)") Bksams_arr%n, bksma_arr_sel%n ,bksma_arr_sel_norm%n, &
        bksma_arr_sel_norm_bd%n
	write(unit=tmpspid,fmt="(I4)") isnap
!#ifdef HDF5	
	!call output_sts_star_samples_hdf5(bksma_arr_sel_norm, "output/pro/MS/"//trim(adjustl(tmpspid)))
	!call output_sts_star_samples_hdf5(bksma_arr_sel_td, "output/pro/MS/td_"//trim(adjustl(tmpspid)))
!#endif	
    !write(*,fmt="(8I10)") Bksams_arr%n, bksma_arr_sel%n ,bksma_arr_sel_norm%n, &
    !    bksma_arr_sel_norm_bd%n

end subroutine

subroutine get_sts_one_species_single(sma,sma_arr,pteve,isnap)
	use com_main_gw
	implicit none
	type (particle_samples_arr_type)::sma_arr
	type(chain_type)::sma
	integer isnap
	!type(core_comp_type)::cc
    type(event_sets_star)::pteve
	
	
end subroutine

subroutine init_pro()
	use com_main_gw
	implicit none

	aomax=rh/2d0/emin_factor
	aomin=rh/2d0/emax_factor !in unit of rh
!	print*, emin_factor,rh,aomin,aomax
!	stop
	
end subroutine