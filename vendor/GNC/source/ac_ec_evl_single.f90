subroutine run_one_sample(pt, run_time)
	use com_main_gw
	use md_coeff
	implicit none
	type(chain_pointer_type),target::pt
	real(8) elp, en0,en1,  steps, run_time, GET_T_GW
    real(8) tgw, period, npi,npf,fpowerlaw
	real(8) rp_out, rtmax, total_time, time, time_next,time_dt, time_create
	integer out_flag_boundary
	real(8) ipdi,ipdf, wsi
    integer flag, flag_dedj,out_flag_clone
    integer,parameter::flag_sg=1,flag_by=2
	integer(8),parameter:: ouput_freq=1,pause_freq=1
	type(coeff_type)::coeNr, coeRR, coeGW
    logical:: init_out, particle_cloned
	integer,save::num_out=0
	interface 
		!subroutine run_one_sample_particle_inside_cluster(pt, time, total_time)
		!	use com_main_gw
		!!	implicit none
		!	type(chain_pointer_type),target::pt
		!	real(8) time, total_time
		!end subroutine
		subroutine run_one_sample_inside_cluster(pt, time, total_time)
			use com_main_gw
			implicit none
			type(chain_pointer_type),target::pt
			real(8) time, total_time
		end subroutine
	end interface

	associate(sample=>pt%ob)

		time=sample%simu_bgtime*1d6*2*pi
		total_time=run_time*1d6*2*pi

		if(sample%weight_real.eq.0d0.or.sample%weight_n.eq.0d0)then
			print*, "ac_ec_evl_single:error, sample%weight_real or n=0"
			call sample%print("ac_ec_evl single")
			stop
		end if
		call reset_sample_init(sample, total_time, time)

		!call check_boundary("1")
		out_flag_boundary=0

		!print*,"source=", sample%source
100  	if(sample%en>ctl%energy_boundary)then
			!if(sample%en0.ne.sample%en.or.sample%jm0.ne.sample%jm)then
			!	print*, "stop", sample%en0, sample%en,sample%jm0, sample%jm
			!	stop
			!end if
			!print*, "run_bd: time, total_time=",time/2d6/pi, total_time/2d6/pi, sample%id
			!select case(ctl%boundary_method)
			!case(boundary_method_fix)
				call run_boundary_state(sample, total_time, time, out_flag_boundary)
			!case(boundary_method_reflesh)
			!	call run_boundary_state_reflesh(sample, total_time, time, out_flag_boundary)
			!end select
			if(out_flag_boundary.eq.100)then
				!if(sample%jm0<0.02d0)then
				!	num_out=num_out+1
				!	print*, "================================="
				!	print*, "num_out=", num_out
				!end if
				time_create=time
				sample%create_time=time_create/2d6/pi

				if(ctl%chattery.ge.3)then
					print*, "cross at", time/2d6/pi
					print*, "en0, jm0=",sample%en0/ctl%energy0, sample%jm0
					print*, "en, jm, time_create=",sample%en/ctl%energy0, sample%jm
				end if
				select type(ca=>sample)
				type is(particle_sample_type)
					call run_one_sample_particle_inside_cluster(pt,time, total_time)
				!type is (sample_type)
				!	call run_one_sample_inside_cluster(pt,time, total_time)
				end select
				!print*, "finished out",sample%exit_flag
				!read(*,*)
				if(sample%exit_flag.eq.exit_boundary_min)then
					sample%jm=sample%jm0; 
					sample%en=sample%en0	
					sample%byot%a_bin=-mbh/2d0/sample%en
					sample%byot%e_bin=(1-sample%jm**2)**0.5
					time=time_create
					!sample%exit_flag=exit_normal
					call set_star_radius(sample%byot%ms)
					!print*, "emin_relesh:sample%radius=",sample%byot%ms%radius
					call get_sample_r_td(sample)
					select type(ca=>sample)
					type is(particle_sample_type)
						call init_particle_sample_common(ca)
					!type is(sample_type)
					!	call init_sample_common(ca)
					end select
					goto 100
				end if
				call update_samples(sample, pt, time_create/1d6/2d0/pi, flag_ini_or)
				if(ctl%chattery>=3)then
					print*, "create at:", time_create/1d6/2d0/pi, bksams%head%ed%ob%en/ctl%energy0, &
				sample%en/ctl%energy0
				end if
			else
				sample%exit_time=time/2d6/pi
			end if
		else
			select type(ca=>sample)
			type is(particle_sample_type)
				call run_one_sample_particle_inside_cluster(pt,time, total_time)
			!type is (sample_type)
			!	call run_one_sample_inside_cluster(pt,time, total_time)
			end select
			if(pt%ob%exit_flag.eq.exit_boundary_min)then
				ctl%num_boundary_elim=ctl%num_boundary_elim+1
			end if
		end if

		!print*, "sample%en,jm=", sample%en/ctl%energy0, sample%jm
		if(ctl%chattery.ge.1)then
			if(ctl%chattery.eq.1)then
				if(sample%exit_flag.ne.exit_boundary_min)then
					select type(sample)
					type is(particle_sample_type)
						call print_results_single(pt,  pt%idx, pt%ed%idx, sample)
					!type is(sample_type)
					!	call print_results(pt,  pt%idx, pt%ed%idx, sample)
					end select
				end if
			else
				select type(sample)
				type is(particle_sample_type)
					call print_results_single(pt,  pt%idx, pt%ed%idx, sample)
				!type is(sample_type)
				!	call print_results(pt,  pt%idx, pt%ed%idx, sample)
				end select
			end if
			if(ctl%chattery.ge.4.or.ctl%debug.ge.1)then
				read(*,*)
			end if
		end if
	end associate

end subroutine

subroutine run_one_sample_particle_inside_cluster(pt, time, total_time)
	use com_main_gw
	use md_coeff
	implicit none
	type(chain_pointer_type)::pt
	!type(particle_sample_type),pointer::sample
	real(8) elp, en0,en1,jm0, jm1,  steps, run_time, GET_T_GW
    real(8) tgw, period,fpowerlaw
	real(8) total_time, time, time_next,time_dt, time_create
	integer(8) j
	real(8) rnd
    integer flag, out_flag_clone,flag_pass_rp
    integer,parameter::flag_sg=1,flag_by=2
	integer(8),parameter:: ouput_freq=1,pause_freq=1
	type(coeff_type)::coeNr
	integer,save::num_out=0
	real(8) dt_block, tp, jlc,ratio

    associate(sample=>pt%ob)
		!if(sample%en>ctl%energy_boundary.and.sample%jm<0.04d0) then
		!	print*, sample%jm,sample%jm0,sample%en
		!	ctl%chattery=4
		!else
		!	ctl%chattery=0
		!end if
		j=0
		
		time_next=time  ! initial time_next
		
		call get_mass_idx(sample%m,sample_mass_idx)

		loop1:do while(time<total_time)
			call update_sample_ej(sample)
			call get_coeff(sample,coeNr)
			!print*, "sample%r_td=",sample%r_td
			call get_sample_r_td(sample)
			!print*, "sample%r_td=",sample%r_td
			call if_sample_within_lc(sample)
			call get_step(sample,coeNr,steps, total_time, time)
		!	print*, "--step finished--"
			en0=sample%en
			jm0=sample%jm
			!print*, "steps, r=", steps, mbh/(-2d0*sample%en)
			!read(*,*)
			if(steps>1d99)then
				print*, "af get_steps steps=",steps,ieee_is_finite(steps)
				call sample%print("sample")
				stop
			end if

			period=P(sample%byot%a_bin)
			time_dt=steps*period
			
			!if(ctl%chattery.ge.4) then
			!	write(*,*) "steps,source=",steps,sample%source
			!	write(*,*) "ac,ec=", sample%byot%a_bin, sample%byot%e_bin
			!	write(*,*) "en, x, en0=", sample%en, sample%en/ctl%energy0, ctl%energy_max
			!	write(*,*) "ceoNr%e, ceoNr%ee=", coenr%e, coenr%ee
			!	!read(*,*)
			!end if

			!=====change da dj according to the current period
			!if(sample%jm>1)then
			!    Print*, "?bf, jm=", sample%jm
			!    stop
			!end if
			!call check("before dedj")

			call get_de_dj(sample, coeNR, time,time_dt,steps, period)
			call get_move_result(sample,sample%den,sample%djp,steps,&
				sample_enf, sample_jmf, sample_mef, sample_af)

			if(ctl%include_loss_cone.ge.1 .and.&
				sample%en<ctl%energy_boundary)then
				
				if(sample%within_jt.eq.1)then	
						call if_sample_pass_rp(sample, steps,flag_pass_rp)
						if(flag_pass_rp.ge.1)then
							select case(sample%obtype)
							case(star_type_ms)
								if(abs(sample%djp0)<sqrt(sample%r_td*mbh*2))then
									sample%exit_flag=exit_tidal_empty
								else
									sample%exit_flag=exit_tidal_full
								end if
								if(ctl%trace_all_sample.ge.record_track_nes.or.&
								sample%write_down_track.ge.record_track_detail)then
									call add_track(time/1d6/(2*pi), sample,state_td)
								end if
								exit loop1
							case(star_type_bh, star_type_ns, star_type_wd, star_type_bd)
								sample%exit_flag=exit_plunge_single

								if(ctl%trace_all_sample.ge.record_track_nes.or.&
								sample%write_down_track.ge.record_track_detail)then
									call add_track(time/2d6/pi,sample,state_plunge)
								end if
								exit loop1
							case default
								print*, "define star type:", sample%obtype
								stop
							end select
						end if				
				end if
				
			end if

			call update_track(sample, j)
			!print*, ctl%trace_all_sample.ge.record_track_detail
			if(sample%write_down_track.ge.record_track_detail&
				.or.ctl%trace_all_sample.ge.record_track_detail)then
				call add_track(time/1d6/(2*pi), sample,state_ae_evl)
			end if

			j=j+1
	!!$			if(mod(j,1000).eq.0) print*, "rid,j=",rid,j,MAX_RUN_LENGTH
			if(j>MAX_RUN_LENGTH)then
			sample%exit_flag=exit_max_reach
			print*, "j=",j
			exit loop1
			end if
			if(j.eq.MAX_RUN_LENGTH/20.or.j.eq.MAX_RUN_LENGTH/2)then
				print*, "single:warning, j, rid=", j, rid
				print*, "step,ao,eo, rp=", steps, sample%byot%a_bin, sample%byot%e_bin, &
						sample%byot%a_bin*(1-sample%byot%e_bin),sample%r_td
				print*, "within_jt=",sample%within_jt
				print*, "sample%obtype=",sample%obtype
				print*, time, period
			end if

			!if(en1>ctl%energy_boundary)then
			!    print*, "trans:", -mbh/2d0/en0/rh, -mbh/2d0/en1/rh, &
			!    -mbh/2d0/ctl%energy_boundary/rh
			!end if
			!call check("before clone")

			if(steps>1d99.or.isnan(steps).or.steps<0)then
				print*, "bf get_dedj steps=",steps,ieee_is_finite(steps)
				call sample%print("ac_ec_evl_single")
				print*, "time, rid=", time, rid
				stop
			end if

			call move_de_dj_one(sample,sample_enf, sample_jmf, sample_mef,sample_af)
			
			!if(sample%jm>1)then
			!    Print*, "?af, jm=", sample%jm
			!    stop
			!end if
			en1=sample%en
			!print*, "time, en1=", time/2d0/pi/1d6, en1
			time_next=time+time_dt

			
			if(en1<ctl%energy_max)then
				sample%exit_flag=exit_boundary_max
				boundary_sts_emax_cros=boundary_sts_emax_cros+1
				!block
				!real(8) wsi
				!call get_wsi(mbh_spin, sample%byot%Inc,wsi)
				!print*,"emax", sample%byot%inc,wsi
				!read(*,*)
				!end block
				
				if(ctl%trace_all_sample.ge.record_track_nes)then
					call add_track(time_next/2d6/pi,sample,state_emax)
				end if
				exit loop1
			end if
			
			if(ctl%clone_scheme.ge.1)then
				!if(ctl%chattery.ge.3) then
				!    print*, "clone",ctl%del_cross_clone, en0, en1
				!end if
				call clone_scheme(pt, en0, en1, ctl%clone_factor(sample_mass_idx),&
					time_next/1d6/2d0/pi, out_flag_clone)
					!if(out_flag_clone.ge.1)then
					!    particle_cloned=.true.
					!end if
					if(out_flag_clone.eq.100)then
						sample%exit_flag=exit_invtransit
						exit loop1
					end if
			end if
			!call check("after clone")

			!======================important notice====================

			if(en1>ctl%energy_boundary)then
				sample%exit_flag=exit_boundary_min
				exit loop1
			end if

			time=time_next
			!print*, "time, total_time=",time/2d6/pi,total_time/2d6/pi
			!print*, "sample%nhiar=",sample%nhiar
			!read(*,*)
			!print*,"run3:",allocated(pt%ob%track), sizeof(pt%ob%track)
		end do loop1
		if(ctl%chattery.ge.4) print*, "exit time=", time/1d6/2d0/pi, sample%exit_flag
		sample%exit_time=time_next/1d6/2d0/pi
		
		if(ctl%chattery.ge.5) then
				write(*,*) "finished, rid, flag, ac=",rid, sample%byot%a_bin
				write(*,fmt="(A25, 1PE11.4)") "-------time exit:",sample%exit_time
				write(*,fmt="(A25, I7)") "-------exit flag:",sample%exit_flag
				print*, "Enter to go"	
				read(*,*)
		end if
	end associate
end subroutine

subroutine print_results_single(pt,  id, eid,sample)
	use com_main_gw
	implicit none
	type(chain_pointer_type):: pt
	type(particle_sample_type)::sample
	integer id, eid,  spid
	character*(3) str_type

	!call get_binary_types(sample%bytype,str_type)
	spid=sample%id
	call get_star_type(sample%obtype,str_type)
	select case(sample%exit_flag)
		case(exit_normal)
			write(chattery_out_unit,fmt="(A25, A5, 5I10, 2I12)") "-------time out",&
                str_type,id,rid,eid, spid, pt%ed%ob%id
                !print*, sample%source
		case(exit_tidal_empty)
			write(chattery_out_unit,fmt="(A25, A5, 5I10, 2I12)") "-------td empty",&
                str_type,id,rid,eid,  spid, pt%ed%ob%id
		case(exit_tidal_full)
			write(chattery_out_unit,fmt="(A25, A5, 5I10, 2I12)") "-------td full",&
                str_type,id,rid,eid,  spid, pt%ed%ob%id
		case(exit_max_reach)
			write(chattery_out_unit,fmt="(A25, A5, 5I10, 2I12)") "-------MORE STEPS NEEDED",str_type,id,rid,eid,&
                      spid, pt%ed%ob%id
			write(chattery_out_unit,fmt="(A25, 1P4E10.3)") "--ac,ec,ain,ein=",&
                sample%byot%a_bin,sample%byot%e_bin
		case (exit_plunge_single)
			write(chattery_out_unit,fmt="(A25, A5, 5I10, 2I12)") "----plunge",&
                str_type,id,rid,eid, spid, pt%ed%ob%id
		case(exit_boundary_min)
			write(chattery_out_unit,fmt="(A25, A5, 5I10, 2I12)") "-------exit_to_emin",str_type,id, &
                rid,eid, spid, pt%ed%ob%id
		case(exit_boundary_max)
			write(chattery_out_unit,fmt="(A25, A5, 5I10, 2I12)") "-------exit_to_emax",str_type,id, &
                rid,eid, spid, pt%ed%ob%id
            !
            if(isnan(sample%byot%a_bin).or.sample%byot%a_bin.eq.0d0)then
                print*, "aout,eout=",sample%byot%a_bin, sample%byot%e_bin,sample%id
                call sample%print("print_results_single")
                print*, "????", rid
                read(*,*)
            end if
		case(exit_ejection)
			write(chattery_out_unit,fmt="(A25, A5, 5I10, 2I12)") "-------escape",str_type,id, &
                rid,eid, spid, pt%ed%ob%id
		case(exit_invtransit)
		case(exit_other)
			write(chattery_out_unit,fmt="(A25, A5, 5I10, 2I12)") "-------other",str_type,id, &
                rid,eid, spid, pt%ed%ob%id
		case default
			write(chattery_out_unit,*) "define state:", sample%exit_flag
			!stop
	end select
end subroutine
