

subroutine print_flux()
    use com_main_gw
    implicit none
    integer i
    write(chattery_out_unit,fmt=*) "print_flux"
    do i=1, ctl%m_bins
        write(chattery_out_unit,fmt="(2I4, 10I9)") rid, i, ctl%bin_mass_flux_in(i,1:ctl%num_bk_comp), &
            ctl%bin_mass_flux_out(i,1:ctl%num_bk_comp), &
            ctl%bin_mass_emax_out(i,1:ctl%num_bk_comp)
    end do
    ctl%bin_mass_flux_in=0
    ctl%bin_mass_flux_out=0
    ctl%bin_mass_emax_out=0
end subroutine
subroutine reset_create_time_zero(sps)
	use com_main_gw
	implicit none
	type(chain_type)::sps
	type(chain_pointer_type),pointer::ps
	real(8) e0
	integer i

	ps=>sps%head
    do while(associated(ps))
        select type(ca=>ps%ob)
        type is(particle_sample_type)
            if(ps%ob%exit_flag.eq.exit_normal)then
                ps%ob%create_time=0
               ! print*, "exit_time=",sp%ob%exit_time, sp%ob%m, sp%ob%source, sp%ob%exit_flag
            end if

        end select
        ps=>ps%next
	end do
end subroutine

subroutine reset_create_time(sps)
	use com_main_gw
	implicit none
	type(chain_type)::sps
	type(chain_pointer_type),pointer::sp

    sp=>sps%head
    do while(associated(sp))
        if(.not.allocated(sp%ob)) then
            print*, "warnning, why it is not allocated?"
            print*, allocated(sp%ob)
            call sps%output_screen()
            stop
        end if
        select type(ca=>sp%ob)
        type is(particle_sample_type)
            if(sp%ob%exit_flag.eq.exit_normal)then
                sp%ob%simu_bgtime=sp%ob%exit_time
                !print*, "exit_time=",sp%ob%exit_time, sp%ob%m, sp%ob%source, sp%ob%exit_flag
            end if
        !type is(sample_type)
        !    if(sp%ob%exit_flag.eq.exit_normal)then
        !        sp%ob%simu_bgtime=sp%ob%exit_time
        !    end if
            !print*, "exit_time=",sp%ob%exit_time, sp%ob%m, sp%ob%source, sp%ob%exit_flag
        end select
        sp=>sp%next
    end do
    !read(*,*)
end subroutine

subroutine reset_j_for_boundary(sps)
	use com_main_gw
	implicit none
	type(chain_type)::sps
	type(chain_pointer_type),pointer::sp

    sp=>sps%head
    do while(associated(sp))
        if(.not.allocated(sp%ob)) then
            print*, "warnning, why it is not allocated?"
            print*, allocated(sp%ob)
            call sps%output_screen()
            stop
        end if
        select type(ca=>sp%ob)
        type is(particle_sample_type)
            if(sp%ob%exit_flag.eq.exit_normal)then
				if(sp%ob%en>ctl%energy_boundary)then
					call set_jm_init(sp%ob)
					sp%ob%jm0=sp%ob%jm
                end if
            end if
        !type is(sample_type)
        !    if(sp%ob%exit_flag.eq.exit_normal)then
        !        if(sp%ob%en>ctl%energy_boundary)then
		!			call set_jm_init(sp%ob)
		!			sp%ob%jm0=sp%ob%jm
        !        end if
        !    end if
            !print*, "exit_time=",sp%ob%exit_time, sp%ob%m, sp%ob%source, sp%ob%exit_flag
        end select
		sp%ob%rp=sp%ob%byot%a_bin*(1-sp%ob%byot%e_bin)
        sp=>sp%next
    end do
    !read(*,*)
end subroutine

subroutine set_real_weight_arr_single(sms_single)
	use com_main_gw
	implicit none
	integer i
	type(particle_samples_arr_type)::sms_single

	do i=1, sms_single%n
        !sms_single%sp(i)%weight_asym=dms%weight_asym
        call get_sample_weight_real(sms_single%sp(i))
	end do
end subroutine

subroutine set_real_weight(sms)
	use com_main_gw
	implicit none
	integer i
	type(chain_type)::sms
	type(chain_pointer_type),pointer::pt

    pt=>sms%head
    do while(associated(pt).and.allocated(pt%ob))
        !pt%ob%weight_asym=dms%weight_asym
        call get_sample_weight_real(pt%ob)
        pt=>pt%next
    end do
	
end subroutine

subroutine convert_sams_convert()
	use com_main_gw
	implicit none

    call arr_to_chain_single(bksams_arr_norm, bksams)
    call reset_create_time(bksams)
    call convert_sams_pointer_arr(bksams, bksams_pointer_arr,type=1)
end subroutine


subroutine arr_to_chain_single(bkarr, chain)
    use com_main_gw
    implicit none
    type(chain_type)::chain   
    type(chain_pointer_type),pointer::pt
    integer i
    type(particle_samples_arr_type) ::bkarr

    call chain%init(bkarr%n)
    pt=>chain%head
    do i=1, bkarr%n
        allocate(particle_sample_type::pt%ob)
        select type(ca=>pt%ob)
        type is(particle_sample_type)
            ca=bkarr%sp(i)
        end select
        pt%idx=i
        pt=>pt%next
    end do
end subroutine

subroutine convert_sams()
	use com_main_gw
	implicit none

    call refine_chain(bksams)
    call convert_sams_pointer_arr(bksams, bksams_pointer_arr,type=1)
    call reset_create_time(bksams)
end subroutine

subroutine convert_sams_no_refine()
	use com_main_gw
	implicit none

    call convert_sams_pointer_arr(bksams, bksams_pointer_arr,type=1)
    call reset_create_time(bksams)
end subroutine

subroutine convert_sams_copy()
	use com_main_gw
	implicit none
	!type(chain_type) ::bysams_out
	type(chain_type) ::bksams_out
        !call sms%head%output()
        !print*, "1"
        !call sms%output_screen(2)
    call chain_select(bksams,bksams_out, exit_normal)
    call bksams_out%copy(bksams)
    call bksams_out%destory()
    call reset_create_time(bksams)
    call convert_sams_pointer_arr(bksams, bksams_pointer_arr,type=1)
    
end subroutine

subroutine run_one_snap(cur_time_i, cur_time_f,  smsa, n,update_dms)
    use com_main_gw
    implicit none
    real(8) cur_time_f, cur_time_i, cur_time_f_istep
    integer n, ierr
	type(particle_samples_arr_type) smsa(n)
    integer,save::norm=5
    real(8) t1, t2
	logical::update_dms

    if(rid.eq.0)then
		call cpu_time(t1)
	end if
    
    call update_arrays_single()
        
    write(unit=chattery_out_unit,fmt=*)  "sg:rid, cur_time_f=", rid, cur_time_f
    call RR_mpi(bksams, cur_time_f)
    call mpi_BARRIER(mpi_comm_world, ierr)
    call collection_data_single_bympi(smsa, ctl%ntasks)
    if(rid.eq.0)then
        print*, "reset time"
    end if
    call reset_create_time(bksams)

    if(rid.eq.0)then
		call cpu_time(t2)
        write(unit=chattery_out_unit,fmt=*)  "one snap running time:", t2-t1
	end if
    
    !stop
	if(update_dms)then
		call set_dm_init(dms)
		!call get_ge_by_root(smsa,bysmsa, ctl%ntasks,.false.)
		if(norm>0)then
			call get_ge_by_root(smsa,ctl%ntasks,.true.)
			call update_weights()
			norm=norm-1
		else
			call get_ge_by_root(smsa, ctl%ntasks,.false.)
		end if
        
    	call get_dms(dms)

		if(rid.eq.0)then
			call dms%print_norm(chatterY_out_unit)
			call print_num_boundary(dms)
			call print_num_all(dms)
		end if
	end if
    
    !call print_balance()
    call mpi_barrier(mpi_comm_world, ierr)
    call reset_create_time(bksams)
	call reset_j_for_boundary(bksams)
    !call check("end of snapshot")
end subroutine

subroutine print_minj()
	use com_main_gw
	implicit none
	type(chain_pointer_type),pointer::pt
	real(8) jmin, jmin0
	integer sid
	pt=>bksams%head
	jmin=1; jmin0=1
	do while(associated(pt))
		if(jmin>pt%ob%jm.and.pt%ob%en>ctl%energy_boundary)then
			jmin=pt%ob%jm
			jmin0=pt%ob%jm0
			sid=pt%ob%id
		end if
		pt=>pt%next
	end do
	print*, "rid, jmin=",rid, jmin,jmin0, sid
end subroutine

subroutine get_vr(r, a,e,vr)
    implicit none
    real(8) r,a,e, ra,rp,vr
    ra=a*(1+e)
    rp=a*(1-e)
    if(r>rp.and.r<ra)then
        vr=r/a*((r-rp)*(ra-r))**(-0.5)
    else
        vr=0
    end if
end subroutine
subroutine get_collection_memory_usage(smsa,n)
    use com_main_gw
    implicit none
    integer n,i, nsize_smsa
    type(particle_samples_arr_type)::smsa(n)
    nsize_smsa=0; 
    do i=1, n
        nsize_smsa=nsize_smsa+sizeof(smsa(i)%sp)/1024
    end do
    nsize_tot=nsize_tot+nsize_smsa
    write(chattery_out_unit,fmt=*) "nsize_smsa=", nsize_smsa
end subroutine
subroutine show_memory_usage()
		use com_main_gw
		implicit none
		integer i, nsize_collision_aux, nsize_other, nsize_cfs

		call get_memo_usage(proc_id)
        nsize_chain_bk=bksams%head%get_sizeof()
		nsize_chain_bk=nsize_chain_bk/1024
		nsize_arr_bk=sizeof(bksams_arr%sp)/1024
		nsize_arr_bk_norm=sizeof(bksams_arr_norm%sp)/1024
		nsize_arr_bk_pointer=sizeof(bksams_pointer_arr%pt)/1024
		nsize_tot_bk=nsize_arr_bk+nsize_chain_bk+nsize_arr_bk_norm+nsize_arr_bk_pointer
		write(chattery_out_unit,fmt=*) "rid:nsize_bk, bksam, bkarr, bkarrnorm, bkp:", &
                nsize_tot_bk, nsize_chain_bk, nsize_arr_bk, &
                nsize_arr_bk_norm , nsize_arr_bk_pointer
        
        !nsize_collision_aux=0
        nsize_chain_by=0
        nsize_arr_by=0
        nsize_arr_by_norm=0
        nsize_arr_by_pointer=0
        nsize_tot_by=0
        
        call cfs%get_size(nsize_cfs)
        !nsize_other=sizeof(rcolld)/1024+sizeof(rc_weight)/1024&
        !    +sizeof(dms)/1024+nsize_cfs
        nsize_tot=nsize_tot_bk+nsize_tot_by!+nsize_other
        write(chattery_out_unit,fmt=*) "rid:cfs:", nsize_cfs!, nsize_collision_aux
    
end subroutine
subroutine show_total_memory_usage()
        use com_main_gw
        implicit none
        integer ierr
        integer nsizetot_all(ctl%ntasks)
        
        call mpi_gather(nsize_tot, 1, MPI_INTEGER,nsizetot_all, 1,MPI_INteger,0,&
             MPI_comm_world,ierr )
        if(rid.eq.0)then
            write(*,fmt=*) "tot memory usage:", sum(nsizetot_all), " kb"
        end if
end subroutine


subroutine print_single_data_size()
    use com_main_gw
    implicit none
    integer n

    call bksams%get_length(n)
    write(chattery_out_unit,fmt="(A20, 20I7)") "rid, bksams, bkarr, bks_norm=",&
        rid, n, bksams_arr%n, bksams_arr_norm%n

end subroutine


subroutine refine_chain(chain)
    use com_main_gw
    implicit none
    type(chain_type)::chain
    type(chain_pointer_type),pointer::pt,ps
    ps=>chain%head
    do while(associated(ps))
        pt=>ps%next
        if(ps%ob%exit_flag.ne.exit_normal)then
            if(allocated(ps%ob))then
                deallocate(ps%ob)
            end if
            if(.not.associated(ps%prev))then
                !ps point to the head of the chain
                print*, "ps point to the head of the chain"
                if(associated(pt))then
                !call chain%output_screen()
                    chain%head=>pt
                    if((.not.associated(pt%prev))) then
                        print*, "??? not associted pt%prev", associated(pt%prev)
                        stop
                    end if
                    pt%prev=>null()
                    call destroy_attach_pointer_chain_type(ps)
                    call pt%set_head()
                else
                    chain%head=>null()
                end if
            else
                call chain_pointer_delete_item_chain_type(ps)
            end if
        end if
        ps=>pt
    end do
    !pt=>null()
    !ps=>null()
end subroutine


subroutine particle_sample_get_weight_clone(en, clone,  amplifier, e0, weight_clone)
    use com_main_gw
	implicit none
	real(8) e0, en
	integer nlvl, clone
	integer amplifier
    real(8) weight_clone

	if(clone.ge.1)then
		if(en/e0<0)then
			!sp%weight=1d0
			!obidx=sp%obidx
			weight_clone=1d0 !sp%weight_asym
		else
            nlvl=int(log10(en/e0))			    
            !if(nlvl>log10clone_emax) nlvl=int(log10clone_emax)-1
            !print*, log10clone_emax
            !read(*,*)
            if(nlvl<0) nlvl=0
			!obidx=sp%obidx
			!sp%weight_clone=dble()**(-dble(nlvl))
			weight_clone=dble(amplifier)**(-dble(nlvl))
			!sp%weight_real=sp%weight_clone*sp%weight_asym
			!print*, "w=", sp%weight, sp%weight0
            
            !if(.not.ieee_is_finite(weight_clone))then
                !print*, "weight_clone=",weight_clone, amplifier, nlvl, en, e0
                !read(*,*)
            !end if
            if(isnan(weight_clone).or..not.ieee_is_finite(weight_clone))then
                print*, "error!weight_clone=", weight_clone
                print*, "amplifier, nlvl, sp%en, e0=", amplifier, nlvl, en, e0
                stop
            end if
		end if
	else
		!obidx=sp%obidx
		weight_clone=1d0
		!sp%weight_real=sp%weight_asym
	end if
end subroutine

subroutine prepare_ini_data(tmprid)
	use com_main_gw
	implicit none
	character*(4) tmprid


    call bksams%input_bin("output/ini/bin/single/samchn"//trim(adjustl(tmprid)))
    call dms%input_bin("output/ini/bin/dms.bin")
    print*, "readin dms finished!"

    call update_arrays_single()
    !if(rid.eq.0) print*, "size of one pksample, bysample (B)", sizeof(bksams_arr%sp(1)), sizeof(bysams_arr%sp(1))
    
end subroutine
