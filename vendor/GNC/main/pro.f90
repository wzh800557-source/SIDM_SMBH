program pro
  use flux_diag
	use com_main_gw
	implicit none
	integer i,isnap,jsp, n
	!type(particle_samples_arr_type)::bksma_sel, bksma_arr_sel
    type(diffuse_mspec)::dms_gene 
    character*(200) tmpspid, tmpssnapid,tmpj, fdir
    logical ex

    call readin_model_par("model.in")
	call init_model()
 	call init_pro()
  flux_enabled = .true.

	call init_sams_events()
	do isnap=1, ctl%n_spshot_total
		print*, "isnap=",isnap
        write(unit=tmpssnapid,fmt="(I4)") isnap
        write(unit=tmpj, fmt="(I4)") ctl%num_update_per_snap
        print*, trim(adjustl(tmpssnapid))//"_"//trim(adjustl(tmpj))
        fdir="output/ecev/dms/dms_"//trim(adjustl(tmpssnapid))//"_"//trim(adjustl(tmpj))
        inquire(file=fdir, exist=ex)
        if(ex)then
            call dms%input_bin(trim(adjustl(fdir)))
        else
            print*, isnap, ctl%num_update_per_snap, "not exist"
        end if
    
		call gethering_samples_single("output/ecev/",isnap, bksams, bksams_arr,ex)
		!call set_weights(bksams_arr,bysams_arr)
		print*,"gathering finished"
        if(ex)then
            call pro_single(isnap)
        end if

	end do	
    
	call output_all_sts("output/")
    
end 

