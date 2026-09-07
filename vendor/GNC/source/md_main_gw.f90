module com_main_gw
	use model_basic
	use md_events
    use md_coeff
    use MPI_comu
    use md_cfuns
	implicit none
!	real(8),parameter::aomax=2, aomin=1d-5 !in unit of rh
	real(8)::aomax, aomin !in unit of rh
	real(8),parameter::einmin=0d0,einmax=0.99d0
	real(8) ::rtmin
	real(8),parameter::rtmin_kl=10
	real(8),parameter::pmin_collision=0.01
	real(8) pmax_collision
	real(8),parameter::period_min=0.1d0*2*pi/365d0
	real(8),parameter::period_max=1d6*2*pi/365d0

    integer::evaluate_dej_method  
    integer,parameter::evaluate_dej_method_rk=1, evaluate_dej_method_grid=2,&
         evaluate_dej_method_cfs=3
	real(8)::min_mp, max_mp
	real(8)::min_ms, max_ms
	real(8)::sample_jc, sample_jlc, sample_jm, sample_jlc_dimless
	real(8)::sample_enf, sample_jmf, sample_mef, sample_af
	integer::sample_mass_idx

	character*(200) DEJ_GRID_FILE_DIR
!	character*(200) record_file_three_body_encounter
!    character*(200) record_file_three_body_encounter_detail
!	character*(200) record_file_kozai_Lidov		
	type(sts_fc_type),allocatable:: ge_profile(:)

contains
	subroutine update_sample_ej(sample)
		implicit none
		class(particle_sample_type)::sample
		real(8) ratio
		
       
		if((abs(sample%en+mbh/(2d0*sample%byot%a_bin)))>1d-5.or.&
			abs(sample%jm-sqrt(1-sample%byot%e_bin**2))>1d-5)then
            print*, "warnning: sample%en, aout=",sample%en, sample%byot%a_bin, &
            -mbh/(2d0*sample%byot%a_bin)
			stop
        end if
		sample%en=-mbh/(2d0*sample%byot%a_bin)
        sample%jm=sqrt(1-sample%byot%e_bin**2)

		sample_jc=sqrt(mbh*sample%byot%a_bin)
		sample_jm=sample%jm*sample_jc

		if(ctl%include_loss_cone.ge.1)then
			sample%rp=sample%byot%a_bin*(1-sample%byot%e_bin)
			ratio=min(2d0, sample%r_td/sample%byot%a_bin)
			sample_jlc_dimless=sqrt(ratio)*sqrt(2d0-ratio)
			sample_jlc=sample_jlc_dimless*sample_jc
		end if
	end subroutine
    subroutine get_coeff(sample, coeNR)
		!use com_main_gw
		!use md_coeff
		implicit none

		real(8) drrjj,drrj, ac, ec, sample_mass
		integer i,idx,idy
		real(8) evjum, even, den, fpowerlaw, dedtmax,dldtmax, period
		type(coeff_type)::coenr
        
		class(particle_sample_type)::sample
	
		if(ctl%chattery.ge.3)then
			print*, "start:x, jm, ac, ec=", sample%en/ctl%energy0, sample%jm, sample%byot%a_bin, &
				sample%byot%e_bin
		end if
		evjum=sample%jm	
		if(sample%jm<jmin_value)then
			sample%jm=2*jmin_value-sample%jm
			sample%byot%e_bin=(1-sample%jm**2)**0.5d0
			evjum=jmin_value
		end if

		if(sample%jm>jmax_value)then
			sample%jm=jmax_value
			sample%byot%e_bin=(1-sample%jm**2)**0.5d0
			evjum=jmax_value
			!stop
		end if
		!		
		even=log10(sample%en/ctl%energy0)
		if(even>log10emax_factor) then
			even=log10emax_factor
		end if
		if(even<log10emin_factor)then
			even=log10emin_factor
		end if
		!print*, "even,evjum=", even, evjum
		call get_coenr(even, evjum, sample%m, sample%en, sample_jc, coenr,idx,idy)
		
        
		
	end subroutine
	
    subroutine create_clone_particle(pt,lvl,amplifier,time)
        !use com_main_gw
        implicit none
        type(chain_pointer_type)::pt
    !	integer,parameter::number_of_clone=9
        type(chain_pointer_type),pointer::ps, pe
        integer lvl, i
        real(8) time
        integer amplifier
        pe=>pt%ed
        call pt%ed%create_chain(amplifier-1)
		!print*, "nodes created:cl"
        ctl%num_clone_created=ctl%num_clone_created+1
        !print*, "create_chain clone"
        !call bysams%output_screen()
        !read(*,*)
        ps=>pe%next
        do i=1, amplifier-1
           ! print*, allocated(pe%ob)
            call pt%copy(ps)
            ps%ob%create_time=time
            ps%ob%simu_bgtime=time

            !ps%ob%source=lvl
            !call ps%ob%get_weight_clone(ctl%clone_scheme, &
            !            amplifier, ctl%clone_e0)
      call particle_sample_get_weight_clone(ps%ob%en, ctl%clone_scheme, &
                        amplifier,ctl%clone_e0,ps%ob%weight_clone)                        
            ps=>ps%next
        end do
    end subroutine
    subroutine reset_sample_init(sample, total_time, time)
        !use com_main_gw
        implicit none
        class(particle_sample_type)::sample
        real(8) total_time, time
        integer ntrack_esti

            if(ctl%chattery.ge.3) then
                if(ctl%ntasks.gt.1)then
                    write(unit=chattery_out_unit,fmt="(A25)") "init", star_type(sample%obtype), rid
                    write(unit=chattery_out_unit,fmt="(A25, 1PE10.3, I10)") "time=",time/(2*pi)/1d6
                    write(unit=chattery_out_unit,fmt="(A25, 1P3E10.3)") "ac,ec,en=", sample%byot%a_bin,&
                        sample%byot%e_bin,sample%en/ctl%energy0
                else
                    write(*,fmt="(A25)") "init", star_type(sample%obtype), rid
                    write(*,fmt="(A25, 1PE10.3, I10)") "time=",time/(2*pi)/1d6
                    write(*,fmt="(A25, 1P3E10.3)") "ac,ec,en=", sample%byot%a_bin,sample%byot%e_bin,sample%en/ctl%energy0
                    !write(*,fmt="(A25, 1P3E10.3)") "blkmass_avg=", ctl%mass_bk_avg
                end if
                
            end if
            sample%exit_flag=exit_normal
    end subroutine

    subroutine create_one_to_chain(pt_chain,sample, time, obj_type, obtype, typeidx,flag)
        implicit none
        integer obtype, typeidx,obj_type, flag
        logical attach
        type(chain_pointer_type)::pt_chain
        type(chain_pointer_type),pointer::pt
        class(particle_sample_type)::sample
        real(8) time, jmi,eni,fpowerlaw
        integer,parameter::flag_sg=1,flag_by=2

        pt=>pt_chain%ed

        !call insert_after_item_chain_pointer(pt) 
        call pt_chain%ed%create_chain(1)
		!print*, "nodes created:bd"

        ctl%num_boundary_created=ctl%num_boundary_created+1
        pt=>pt%next
        select case(obj_type)
        case(flag_sg)
            allocate(particle_sample_type::pt%ob)
        end select
    
        call pt%ob%init()
        !call set_init_by_type(ps%ob%obtype, ps%ob%obidx)
        pt%ob%obtype=obtype; pt%ob%obidx=typeidx
        pt%ob%weight_real=sample%weight_real
        !pt%ob%weight_n=sample%weight_n
        !print*, sample%byot_bf%e_bin, sample%jm
        
        select case(flag)
		case(flag_ini_or)
			pt%ob%en=sample%en0
			pt%ob%jm=sample%jm0
			pt%ob%byot%e_bin=(1-pt%ob%jm**2)**0.5d0
			pt%ob%byot%a_bin=-mbh/2d0/pt%ob%en
            !print*, "en0, jm0", pt%ob%en, pt%ob%jm
        case default
            !print*, "error! define flag"
            !stop
        end select

        select case(obj_type)
        case(flag_sg)
            !print*, "sg type:", pt%ob%obtype, pt%ob%obidx
            select type(ca=>pt%ob)
            type is(particle_sample_type)
                call init_particle_sample_one(ca, sample%m, flag)
                !print*, ca%byot%a_bin
            class default
                print*, "error! why?"
                stop
            end select
        end select
        !

        pt%ob%create_time=time
        pt%ob%simu_bgtime=time
        !if(ctl%clone_scheme.ge.1)then
        !    pt%ob%nhiar=get_lvl(pt%ob%en)
        !else
        !    pt%ob%nhiar=0
        !end if
        !print*, "pt%ob%nhair=", pt%ob%nhiar, pt%ob%en, pt%ob%en0
        !read(*,*)
        pt%ob%weight_N=sample%weight_N
        !pt%ob%weight_asym=sample%weight_asym
        !pt%ob%weight_clone=1d0
        !pt%ob%weight_real=pt%ob%weight_clone*pt%ob%weight_asym*pt%ob%weight_n
        
        if(ctl%clone_scheme.ge.1)then
            call create_init_clone_particle(pt, pt%ob%en0,time)
            !print*, "create?"
            !stop
        end if
        !print*, "en0, jm0", pt%ob%en0, pt%ob%jm0
        !read(*,*)
        
    end subroutine
	subroutine if_sample_within_lc(sample)
		implicit none
		class(particle_sample_type)::sample
		sample%rp=sample%byot%a_bin*(1-sample%byot%e_bin)
		if(sample%rp<sample%r_td)then
			sample%within_jt=1
		else
			sample%within_jt=0
		end if
	end subroutine
	subroutine if_sample_pass_rp(sample,steps, flag)
		implicit none
		class(particle_sample_type)::sample
		integer flag
		real(8) ipdi,ipdf
		real(8) steps, npi,npf

		npi=sample%byot%me
		npf=npi+steps*pi*2
		!print*, "loss cone: sample%me=",npi
		ipdi=npi/2d0/pi-int(npi/2d0/pi)
		ipdf=npf/2d0/pi-int(npf/2d0/pi)
		!print*, "ipdi,ipdf=",ipdi,ipdf
		if(ipdi<0.5.and.ipdf>0.5.or.steps.ge.1d0)then
			flag=1
		else
			flag=0
		end if
	end subroutine
    subroutine get_step(sample,coeNr, steps, ttot, tnow)
        implicit none
        integer(8) nstemp
        integer num
        real(8) tgw_otby,get_t_gw, r_mean
        real(8)  period, steps
        real(8),external:: rnd
        real(8) NminNr, NminRR, ttot, tnow, ac, na,mass_bk_avg
        real(8) e2,e1, jnr2, jrr2, jgwr, j1,time_dt_nr
        real(8) nmjm, nmjr, nmjl, nmjmr, nmjrr, jm, nmen, nmendrift
        real(8) nmende2, nmendee2,deb2, njdrift
        real(8) ntmax, nmingwe, nmingwj!,lambda(ctl%m_bins)
        real(8),parameter::avvr=0.31
        real(8) tvrr_step, ratio,wsi
        class(particle_sample_type)::sample
        type(coeff_type)::coenr
        !print*, "1"
        !print*, sample%byot%a_bin

        period=P(sample%byot%a_bin)
        !print*, "1.5"
		!if(sample%en/ctl%energy0<0.5)then
		!	steps=5e5/Period
		!	return
		!end if
        !Jmax=sqrt(mbh*sample%byot%a_bin)
        !print*, "2"
        
        if(.not.isnan(sample%byot%a_bin))then
            sample%en=-mbh/(2*sample%byot%a_bin)
        else
            print*, "get_step:aot = NaN", sample%byot%a_bin, sample%byot%e_bin
            stop
        end if
        if(sample%byot%e_bin<1d0)then
            sample%jm=sqrt(1-sample%byot%e_bin**2)
        else
            select type(sample)
            type is(particle_sample_type)
                print*, "particle:get_step:eot >1 ???", sample%byot%a_bin, sample%byot%e_bin
                print*, "en,jm=",sample%en, sample%jm
            !type is(sample_type)
            !    print*, "by:get_step:eot >1 ???", sample%byot%a_bin, sample%byot%e_bin
            end select
            print*, "id,type=",sample%id,sample%obtype, sample%obidx
            print*, "m=",sample%m
            read(*,*)
        end if
        jm=sample%jm*sample_jc
		
        call get_steps_nr_EJ(sample%en, sample%jm, coenr, sample_jc,  time_dt_nr)
		
		steps=time_dt_nr/period
		ntmax=(ttot-tnow)/Period
        steps=min(steps,ntmax)

        if(steps.eq.1d6) then
            print*, "warnning steps=", steps
            print*, "nmen, nmjm, nmjr=", nmen, nmjm, nmjr, ntmax
            print*, "nmendrift,steps=", nmendrift,steps
            print*, "e1,e2=",e1,e2, coeNr%ee, coeNr%e
            print*, "P,ac=", Period, sample%byot%a_bin, sample%byot%e_bin, sample%en
            print*, "id, rid=",sample%id, rid
            stop
        end if
        
        !print*, "steps=", steps
        if(isnan(steps))then
            print*,"1:steps is nan:steps,  nmen=", &
                        steps,  nmjm, nmjr, nmjm, nmen
            print*, "jnr2, e1, e2, jm=", jnr2, e1, e2, jm
            print*, "sample%id=", sample%id
            print*, "sample%byot%a_bin,ec, period=",sample%byot%a_bin,sample%byot%e_bin, period
            read(*,*)
        end if
        
        if(sample%r_td<0)then
            print*, "get_step:error, sample%r_td<0", sample%r_td
            print*, "sample%id=",sample%id
            call print_binary(sample%byot)
            select type(sample)
            type is(particle_sample_type)
                print*, "particle"
            end select
            stop
        end if

        if(ctl%include_loss_cone.ge.1 .and.&
		sample%en<ctl%energy_boundary)then
            
				jnr2=coenr%jj*period
				nmjl=(max(0.1*sample_jlc, 0.25*abs(jm-sample_jlc)))**2/jnr2
				steps=min(steps, nmjl)
			!end if
        end if
        if(steps.le.0)then
            print*, "ao,eo=",sample%byot%a_bin, sample%byot%e_bin
            print*, "nmjl,  nmjm, nmjr, nmjm, nmen=",nmjl,  nmjm, nmjr, nmjm, nmen, ntmax
            print*, "ttot, tnow=", ttot, tnow
            stop
        end if
       
        if(isnan(steps).or.steps<0)then
            print*,"2:steps is nan:steps,  nmjl, nmjr, nmjm, nmen=", &
                        steps,  nmjl, nmjr, nmjm, nmen
            print*, "jnr2,jmin,jmax,rtd/ac=", jnr2, sample_jlc,sample_jc, sample%r_td/sample%byot%a_bin
            !print*, "coeNr%jj,Period=", coeNr%jj,Period
            
            print*, "sample%byot%a_bin,ec, period=",sample%byot%a_bin,sample%byot%e_bin, period
           ! print*, "lambda_aux=",lambda_aux
            read(*,*)
        end if
        
    end subroutine
	subroutine get_steps_nr_EJ(en, jm,coenr, jc,time_dt_nr)
		implicit none
		real(8) steps,period, jc, jm,en, enev
		type(coeff_type)::coenr
		real(8) time_dt_e, time_dt_j, time_dt_nr
		if(coenr%ee.ne.0)then
			time_dt_e=min((en*0.15)**2/coenr%ee, abs(en*0.15)/abs(coenr%e))
		else
			time_dt_e=1d6!*period
		end if
		if(coenr%jj.ne.0)then
			time_dt_j=min((jc*0.1)**2/coenr%jj, &
			(0.4d0*(1.0075-jm)*jc)**2/coenr%jj)
		else
			time_dt_j=1d6!*period
		end if
		time_dt_nr=min(time_dt_e,time_dt_j)

		if(ctl%chattery.ge.3)then
			print*, "=========get steps NR==============="
			print*, "sample%x, jm=", en/ctl%energy0, jm
			print*, "time_dt_nr=", time_dt_nr
			print*, "=========end of get steps NR========"
        end if
	end subroutine

	subroutine get_steps_nr_xj(en, jm,coenr, time_dt_nr)
		implicit none
		real(8) steps,period,  jm,en, enev
		type(coeff_type)::coenr
		real(8) time_dt_e, time_dt_j, time_dt_nr
	
		if(coenr%ee.ne.0)then
			enev=en/ctl%energy0
			time_dt_e=min((enev*0.15)**2/coenr%ee, abs(enev*0.15)/abs(coenr%e))
		else
			time_dt_e=1d6!*period
		end if
		if(coenr%jj.ne.0)then
			time_dt_j=min(0.1d0**2/coenr%jj, &
			(0.4d0*(1.0075-jm))**2/coenr%jj, &
			(0.25*abs(jm))**2/coenr%jj)
		else
			time_dt_j=1d6!*period
		end if
		
		time_dt_nr=min(time_dt_e,time_dt_j)

		if(ctl%chattery.ge.3)then
			print*, "=========get steps NR==============="
			print*, "sample%x, jm=", en/ctl%energy0, jm
			print*, "time_dt_nr=", time_dt_nr
			print*, "=========end of get steps NR========"
        end if
	end subroutine
	subroutine get_sample_r_td(sp)
		implicit none
		class(particle_sample_type)::sp
		real(8) wsi

		select type(sp)
		type is(particle_sample_type)
			call get_sample_r_td_single(sp)
		class default
			print*, "?? which type??"
			stop
		end select
		
	end subroutine

    subroutine set_jm_init(bkps)
		implicit none
		class(particle_sample_type)::bkps
		real(8) rtd, jlc,tmp
		real(8),external::fpowerlaw, rnd
		
		select case(ctl%boundary_fj)
		case(boundary_fj_iso)
            bkps%jm=fpowerlaw(1d0,0.0044d0,0.99999d0)
            bkps%byot%e_bin=(1-bkps%jm**2)**0.5

		case(boundary_fj_ls)
200		 	bkps%jm=fpowerlaw(1d0,0.01d0,0.99999d0)
			select type (bkps)
				type is(particle_sample_type)
					select case(bkps%obtype)
					case(star_type_ms)
						rtd=(3*mbh/bkps%m)**(1/3d0)*bkps%byot%ms%radius
						jlc=(1-(1-rtd/ctl%rbd)**2)**0.5
					case(star_type_bh,star_type_ns,star_type_wd,star_type_bd)
						jlc=4*(mbh_radius/bkps%byot%a_bin)**0.5d0
					case default
						print*, "star type??", bkps%obtype
						stop
					end select
			end select
			tmp=rnd(0d0,1d0)
			if(bkps%jm<jlc) goto 200
			if(tmp>log(bkps%jm/jlc)/log(1d0/jlc)) goto 200
			!end if
			bkps%byot%e_bin=(1-bkps%jm**2)**0.5
		case default
			print*, "define flag INI", ctl%boundary_fj
			stop
		end select
	end subroutine
    subroutine update_samples(sample, pt,time,flag_bd)
        implicit none
        class(particle_sample_type)::sample
        type(chain_pointer_type)::pt
        integer  obtype, typeidx, flag_bd
        real(8) time
        integer,parameter::flag_sg=1,flag_by=2
        !note that time is in unit of Myr!

        if(ctl%chattery.ge.4)then
            print*, "new particle created", pt%idx,rid
        end if
        !select case(ctl%ini_sample_mode)
        !case(ini_sample_mode_given)
            !call select_reini_types(sample%m, flag, obtype, typeidx)
            obtype=sample%obtype
            typeidx=sample%obidx

            if(ctl%chattery.ge.4)then
                print*, "created type:", obtype!, flag
            end if
            select type (sample)
            type is (particle_sample_type)
                call create_one_to_chain(bksams%head%ed, sample, time, flag_sg, &
                    obtype, typeidx,flag_bd)
                !call create_one_to_chain(pt, sample, time, flag_sg, &
                !    obtype, typeidx,flag_bd)
                !    print*, "created bksams ed"
                
                !   print*, "=>sg"
            end select
            if(ctl%chattery.ge.4)then
                print*, "particle created"
            end if
        !case(ini_sample_mode_mobse)
        !    print*, "write the code here!"
        !end select
    end subroutine
    subroutine get_type_idx(sample, idxob)
        implicit none
        class(particle_sample_type)::sample
        integer idxob
        select type(sample)
            type is(particle_sample_type)
                select case(sample%obtype)
                case(star_type_ms)
                    idxob=1
                case(star_type_bh)
                    idxob=2
				case(star_type_ns)
					idxob=3
				case(star_type_wd)
					idxob=4
				case(star_type_bd)
					idxob=5
                case default
                    print*, "error in get_type_idx:sg",sample%id
                    stop
                end select
        end select
    end subroutine
    subroutine update_track(sp,j)
        implicit none
        class(particle_sample_type)::sp
        integer(8) j
        
        if((sp%write_down_track.ge.record_track_nes&
        .or.ctl%trace_all_sample.ge.record_track_all))then
            !print*, "sp%track_step=",sp%track_step
            !read(*,*)
            if(mod(j,sp%track_step).eq.0)then
                if(sp%length_to_expand>MAX_LENGTH)then
                    sp%track_step=sp%track_step*10 
                    call track_compress(sp, 10)
                    print*,"track compressed"
                end if
            end if
    end if
    end subroutine
    subroutine track_compress(sp, ns)
        use model_basic
        implicit none
        class(particle_sample_type)::sp
        type(track_type),allocatable::tk(:)
        integer ns,j,i
    
        allocate(tk(sp%length))
        tk=sp%track
        j=0
        do i=1, sp%length, ns
            j=j+1
            sp%track(j)=tk(i)
        end do
        sp%length=j
    end subroutine

    subroutine get_de_dj(sample,coeNr,  time, dt,steps, period)        !use com_main_gw
        !use md_coeff
        implicit none
        class(particle_sample_type)::sample
        real(8),intent(in) ::dt, steps,time
        real(8) den, gen_gaussian
        real(8) y1,y2,y3,y4,rho, deb2, jc, jum, n2j
        real(8) dpe,ai, ei, Eni,Enf, Ji,Jf, af, ef
        integer jb, ju
        real(8) ipdi, ipdf
        real(8) period, npi, npf
        real(8),save::t0=0
        type(coeff_type)::coeNr 
        
        if(isnan(dt).or.(.not.ieee_is_finite(dt)).or.(.not.ieee_is_finite(steps)))then
            print*, "dt=!, dt, steps=", dt, steps
            stop
        end if

        if(isnan(sample%jm).or.isnan(sample%en))then
            print*
            print*, "ac,ec=",sample%byot%a_bin,sample%byot%e_bin
            print*, "en,jm=",sample%en, sample%jm,dt
            read(*,*)
        end if
        !if(sample%en>ctl%energy_min.or.sample%en<ctl%energy_max)then
        !    flag=1
        !    return
        !end if
        !---for debug only, can be removed later
        !ai=sample%byot%a_bin
        !ei=sample%byot%e_bin
        !---

		call get_de_dj_nr(coenr, dt, steps, sample%den, sample%djp, &
			sample%djp0)
        
        if(isnan(sample%djp))then
            print*, "djp=",sample%djp, coeNR%jj, dt, coeNR%j
        end if
        
        
        
    end subroutine
	
	subroutine get_de_dj_nr(coenr, dt, steps, den, djp,djp0)
		implicit none
		real(8),intent(out):: den, djp, djp0
		type(coeff_type)::coenr
		real(8) jum, rho, y1, y2, dt,n2j,steps

		!jum=sample%jm*sqrt(mbh*sample%byot%a_bin)
        rho=coeNR%ej/sqrt(abs(coeNR%ee*coeNR%jj))  
		!print*, "jm,rho=", sample%jm, rho
        !if(rho>=1)then
        !    print*, "rho,ej,ee,jj=", rho, coeNR%ej, coeNR%ee, coeNR%jj
        !    read(*,*)
        !endif
        call gen_gaussian_correlate(y1,y2,rho)
        y1=max(min(y1,6d0),-6d0)
        y2=max(min(y2,6d0),-6d0)
        den=coeNR%e*dt+y1*sqrt(coeNR%ee*dt)
        !sample%den=coeNR%e*dt+y1*sqrt((coeNR%ee+coeNr%e**2)*dt)
        n2j=sqrt(coeNR%jj*dt)
        !if(n2j<jum/4d0)then
            djp=coeNR%j*dt+y2*n2j
		!	print*, "jum=",jum
		!	print*, "1:djp=",sample%djp, y2
        !else
        !    y3=gen_gaussian(1d0)
        !    y4=gen_gaussian(1d0)
        !    sample%djp=sqrt((jum+n2j*y3)**2+(n2j*y4)**2)-jum
		!	print*, "2:djp=",sample%djp, y3,y4
		!	read(*,*)
        !    !print*, sample%djp, jum, sample%jm
        !    !read(*,*)
        !end if
		
		if(ctl%chattery.ge.3)then
			print*, "======start dedj========"
			print*, "dedrift, scatter=", coenr%e*dt, y1*sqrt(coeNR%ee*dt)
			print*, "djdrift, scatter=", coenr%j*dt, y2*n2j
			print*, "======end dedj========"
		end if
        djp0=coeNR%j*dt/steps+y2*sqrt(coeNR%jj*dt/steps)
	end subroutine
    subroutine get_sample_weight_real(sp)
        implicit none
        class(particle_sample_type)::sp
        sp%weight_real=sp%weight_clone*dms%weight_asym*sp%weight_n*ctl%n_basic
        !if(sp%weight_real.eq.0d0)then
        !    print*, "get_sample_weight_real:sp%wreal=0"
        !    print*, "sp%weight_clone,dms%weight_asym,sp%weight_n,ctl%n_basic"
        !    print*, sp%weight_clone,dms%weight_asym,sp%weight_n,ctl%n_basic
        !    call sp%print("get_sample_weight_real")
        !    stop
        !end if
    end subroutine

    subroutine output_sample_track_txt(sp,fl)
        class(particle_sample_type)::sp
        character*(*) fl
        select type(ca=>sp)
        type is (particle_sample_type)
            call output_sg_sample_track_txt(ca,fl)
       ! type is (sample_type)
        !    call output_by_sample_track_txt(ca,fl)
        end select
    end subroutine
    subroutine add_track(t,sp, state_flag)
        !use model_basic
        implicit none
        class(particle_sample_type) ::sp
        type(track_type),allocatable:: tk(:)
        real(8) t
        integer i, state_flag
        i=sp%length
       ! print*, "i=",i, sp%length_to_expand
        if(i.eq.sp%length_to_expand)then
          if(ctl%chattery.ge.4) then
              print*, "track expanding"
          end if
          sp%length_to_expand=sp%length_to_expand+track_length_expand_block
          allocate(tk(i))
          tk(1:i)=sp%track(1:i)
          if(allocated(sp%track)) deallocate(sp%track)
          allocate(sp%track(sp%length_to_expand))
          sp%track(1:i)=tk(1:i)
        end if
        i=i+1
        sp%length=i
        sp%track(i)%time=t
        sp%track(i)%ac=sp%byot%a_bin
        sp%track(i)%ec=sp%byot%e_bin
        sp%track(i)%incout=sp%byot%inc
        sp%track(i)%omout=sp%byot%om
        sp%track(i)%state_flag=state_flag
		
     !	print*, sp%length
     end subroutine
     subroutine run_boundary_state(sample, total_time, time_now, flag)
        implicit none
        class(particle_sample_type)::sample
        real(8) total_time, time_now, time_dt, even, evjum, period
        real(8) deb, jum, rho, y1, y2, n2j,y3,y4,gen_gaussian, djp, den
        real(8) en, jc,  time_dt_nr,time_dt_t,steps, djp0
        integer jb, ju
        integer flag,idx,idy
        type(coeff_type)::coenr_step, coenr
        flag=0
        if(time_now.ge.total_time)then
            return
        end if

        !even=log10(ctl%x_boundary); evjum=jmin_value
        !deb=ctl%energy_min-ctl%energy_boundary
        !en=ctl%energy_boundary
        !jc=mbh/(-2d0*en)**0.5d0
        !call get_coenr(even, evjum,sample%m, en, jc, coenr_step)
        !time_dt=(0.5*deb)**2/coenr_step%ee
        !print*, "deb, coenr%ee=", deb, coenr%ee
        !print*, "time, time_dt,time_tot=", time_now/2d0/pi/1d6, time_dt/2d0/pi/1d6, &
        !    total_time/2d0/pi/1d6
        !print*, "en, eb=" ,sample%en/ctl%energy0, ctl%x_boundary
         
!=============================================       
        evjum=max(min(sample%jm, jmax_value), jmin_value)
        even=log10(sample%en/ctl%energy0)
        en=sample%en
        jc=mbh/(-2d0*en)**0.5d0

        call get_coenr(even, evjum,sample%m, en, jc, coeNr,idx,idy)
		period=P(sample%byot%a_bin)
100     do while (time_now<total_time) 
            !time_dt=min((0.15*deb)**2/coenr%ee, total_time-time_now)
            time_dt_t=total_time-time_now
            !time_dt_e=(sample%en*0.15)**2/coenr%ee
			
			call get_steps_nr_EJ(en, sample%jm, coenr, jc,  time_dt_nr)
			
            time_dt=min(time_dt_nr, time_dt_t)
            time_now=time_now+time_dt
			steps=time_dt/period


            !jum=sample%jm*sqrt(mbh*sample%byot%a_bin)
			call get_de_dj_nr(coenr,time_dt,steps, den, djp, djp0)
            
    !=============================================
			if(ctl%chattery.ge.4)then
				print*, "bd:sample%en, jm, den, djm=", sample%en, &
					sample%jm, den, djp
			end if

            !print*, "time, den, sample%en+den=",time_now/2d0/1d6/pi,&
            !    den/ctl%energy0, (sample%en+den)/ctl%energy0
            if(sample%en+den<ctl%energy_boundary)then
                
                sample%den=den
                sample%djp=djp
                !call move_de_dj_one(sample,den,djp,steps)
				call get_move_result(sample,den,djp,steps,&
				sample_enf, sample_jmf, sample_mef, sample_af)
				call move_de_dj_one(sample,sample_enf, sample_jmf, sample_mef,sample_af)

                flag=100
                !print*, "exit time_now=", time_now/2d0/pi/1d6
                !print*, "sample%en, sample%byot%a_bin=", sample%en, sample%byot%a_bin
                !read(*,*)
                return
            end if
        
        end do
     end subroutine
	 
	 subroutine get_move_result(sample,den,djp,steps, &
			enf, jmf, mef, af)
		implicit none
		class(particle_sample_type)::sample
		real(8) den, djp,steps
		real(8) ai,ei,Eni, Enf, Ji, Jf, af
		real(8) jmf, mef

		sample%byot%a_bin=-mbh/(2*sample%en)
		ai=sample%byot%a_bin
		ei=sample%byot%e_bin		
        Eni=sample%en
        Ji=sample%jm*sqrt(mbh*ai)
        Enf=Eni+den
        af=-mbh/(2*Enf)
        jf=Ji+djp
        jmf=jf/sqrt(mbh*af)
        
		mef=sample%byot%me+(steps-int(steps))*2*pi
		!print*, "me=", sample%byot%me
10      if(jmf<jmin_value) then
            jmf=2*jmin_value-jmf
        end if
        if(jmf>jmax_value)then
			jmf=jmax_value
        end if
        if(jmf<jmin_value) goto 10
        
	 end subroutine
	 subroutine move_de_dj_one(sample, enf, jmf, mef,af)
		implicit none
		class(particle_sample_type)::sample
		real(8) den, djp,steps
		real(8) enf, jmf, mef,af
		real(8) eni, ai, ei,ji

		sample%byot%a_bin=-mbh/(2*sample%en)
		ai=sample%byot%a_bin
		ei=sample%byot%e_bin
		ji=sample%jm
		sample%byot_bf%a_bin=sample%byot%a_bin
		sample%byot_bf%e_bin=sample%byot%e_bin
	
		sample%byot%a_bin=af
		sample%en=enf
		sample%jm=jmf
        sample%byot%e_bin=sqrt(1-sample%jm**2)

		if(mef>2*pi)then
			sample%byot%me=mef-2*pi
		else
			sample%byot%me=mef
		end if
		
                
        if(sample%byot%e_bin.eq.1d0)then
        !	print*, "error: sample%byot%e_bin=1d0, coeNr%j, coeNr%jj, sample%jm, sample%djp=", &
        !			coeNR%j, coeNR%jj, sample%jm, sample%djp
        !	stop
            sample%byot%e_bin=0.9999999999d0
        end if

		if(isnan(sample%jm).and.sample%en<0d0)then
			print*, "error! sample%jm=NaN", sample%jm
			print*, "sample%m,id=",sample%m,sample%id
			select type(sample)
			type is(particle_sample_type)
				print*, "type is particle"
			end select

            print*, "enmin,enmax=",ctl%energy_min,ctl%energy_max
            print*, "jm, en=", sample%jm, sample%en, sample%den
            print*, "sample%byot%a_bin,djp=",sample%byot%a_bin,sample%djp
			stop
        end if

		if(ctl%chattery.ge.3)then
			print*, "========start de_dj======"
			print*, "befor: sample%en,jm,ac,ec=",mbh/2./ai/ctl%energy0,&
			Ji, ai, ei
			print*, "after: sample%en,jm,ac,ec=",sample%en/ctl%energy0, sample%jm, &
				sample%byot%a_bin, sample%byot%e_bin
			print*, "========end of de_dj======"
		end if
	 end subroutine
	 
end module



subroutine get_rh_vh_nh(mbh, rh, nh, vh)
        use constant
        use sidm_params, only: use_sidm, n0_phys_pc3, use_gravothermal_n0, &
                               use_spike_n0, spike_gamma, spike_mass_factor
        implicit none
        real(8) mbh, rh, nh, vh
        real(8),parameter::rgc=8.32d3
        real(8),parameter::r0=3.1
        real(8),parameter::n0=2d4
        real(8),parameter::PI_loc=3.14159265358979d0
        real(8) :: n0_use

        n0_use = n0
        if (use_sidm .and. use_gravothermal_n0) n0_use = n0_phys_pc3
        rh=r0*(Mbh/4d6)**0.55*pc
        vh=sqrt(mbh/rh)
        if (use_sidm .and. use_spike_n0) then
            nh = spike_mass_factor*mbh*(3d0-spike_gamma)/(4d0*PI_loc*rh**3)
        else
            nh=n0_use/pc**3*(mbh/4e6)**(-0.65)
        end if
end subroutine
subroutine init_model_ctl()
	use com_main_gw
	implicit none
	integer i
    real(8) ra

    !real(8),parameter::massfactor=0.36
	!mbh=4d6
	!rh=2.31*(Mbh/4d6)**0.5*206264	
    
	ctl%grid_type=sts_type_grid
	!method_int=method_int_linear

	emin_value=(1-jmax_value**2)**0.5
	emax_value=(1-jmin_value**2)**0.5

	call get_rh_vh_nh(mbh, rh, ctl%n0, ctl%v0)
    log10rh=log10(rh)
    rhmin=rh/emax_factor/2d0; 
    rhmax=rh/emin_factor/2d0
    !===================test=============
    !rh=2.31*(Mbh/4d6)**0.5*206264	
    !====================================

!	ctl%sigma=(mbh/3.16d8)**(1d0/4.42d0)*200/29.79    ! K & Ho 2013
!	print*, ctl%sigma
	ctl%sigma=sqrt(mbh/rh)
!	print*, ctl%sigma
!	stop
	ctl%energy0=-mbh/rh
	ctl%energy_min=ctl%energy0*emin_factor; ctl%energy_max=ctl%energy0*emax_factor
    
	!ctl%n0=MBH*(3d0-7d0/4d0)/(4d0*pi*rh**3);
    !ctl%n0=MBH/(pi*rh**3)*massfactor;
    !print*, ctl%n0*pc**3*massfactor
    !ra=rh/pc/rgc*180/pi*3600.
    !ctl%n0=1.2e6/pc**3*(ra/10)**(-2.0)*(mbh/4d6)**(-0.65)
    
    !print*, "ctl%n0=",ctl%n0
    !===========test================
    !ctl%n0=MBH*(3d0-7d0/4d0)/(4d0*pi*rh**3)
    !================================
    !print*, "ctl%n0=", ctl%n0*pc**3
    !stop
    !ctl%n0=3.7d4/pc**3*(mbh/4e6)**0.5
    !print*, "ctl%n0=",  ctl%n0*pc**3

	!ctl%energy_boudary_min=energy_min
	!ctl%r0bg=rh*a0_rh

	call set_simu_time()
	
	ctl%rbd=rh*0.5d0/ctl%x_boundary
	if(ctl%clone_scheme.ge.1)then
		ctl%clone_e0=clone_e0_factor/10*ctl%energy0
		clone_emax=ctl%energy_max/ctl%clone_e0
        log10clone_emax=log10(clone_emax)
	end if
	do i=1, ctl%num_bk_comp
		associate(ci=>ctl%cc(i))
			!ci%r0_vanish=18.91*(mbh/4e6)**(-1d0/6d0)*(ci%blkmass)**(2d0/3d0)
			ci%n0=ci%n_in_rh/(4d0/(3d0-ci%alpha)*pi*rh**3);
		end associate
	end do
	
    ctl%energy_boundary=ctl%x_boundary*ctl%energy0
    !ctl%bd_up=ctl%energy_boundary*10**(-ctl%bd_thickness)
    !ctl%bd_dn=ctl%energy_boundary*10**(ctl%bd_thickness)
    !if(rid.eq.0)then
    !    print*, "bd_thickness=", ctl%bd_thickness
    !    print*, "bd_up, bd, bd_dn, bd_clone=", ctl%bd_up/ctl%energy0, ctl%x_boundary, &
    !        ctl%bd_dn/ctl%energy0, ctl%clone_e0/ctl%energy0
    !end if
    !if(ctl%bd_up>ctl%energy_min)then
    !    print*, "error!,ctl%bd_up>ctl%energy_min",ctl%bd_up,ctl%energy_min
    !    stop
    !end if
    !if(ctl%bd_dn<ctl%clone_e0)then
    !    print*, "error!,ctl%bd_dn<ctl%clone_e0",ctl%bd_dn,ctl%clone_e0
    !    stop
    !end if
    !ctl%birth_position=3.5*rh
    !ctl%x_birth=0.5d0/3.5d0
    !ctl%rh_min=1d0/emax_factor/2d0; 
    !ctl%rh_max=1d0/emin_factor/2d0
    
    do i=1, ctl%num_bk_comp
        ctl%cc(i)%alpha_ini=7d0/4d0
    end do
	ctl%bin_mass_particle_number(1:ctl%m_bins)=mbh/ctl%ini_weight_n(1:ctl%m_bins)*ctl%asymptot(1,1:ctl%m_bins)
    ctl%n_basic=minval(ctl%ini_weight_n(1:ctl%m_bins))
!    do i=1, ctl%num_bk_comp
!        ctl%n_basic=min(ctl%n_basic, ctl%cc(i)%N_pre_in)
!    end do
!    do i=1, ctl%by_sp_num
!        ctl%n_basic=min(ctl%n_basic, ctl%bs(i)%N_pre_in)
!    end do
!    do i=1, ctl%num_bk_comp
!        ctl%cc(I)%N_pre=int(ctl%cc(i)%N_pre_in/ctl%n_basic)
!        if(abs(ctl%cc(I)%N_pre-anint(ctl%cc(i)%N_pre_in/ctl%n_basic))>1d-6)then
!            print*, "error! N_pre should be integer times:", i, ctl%cc(i)%N_pre_in
!        end if
!    end do
!    do i=1, ctl%by_sp_num
!        ctl%bs(I)%N_pre=int(ctl%bs(i)%N_pre_in/ctl%n_basic)
!        if(abs(ctl%cc(I)%N_pre-anint(ctl%cc(i)%N_pre_in/ctl%n_basic))>1d-6)then
!            print*, "error! N_pre should be integer times:", i, ctl%cc(i)%N_pre_in
!        end if
!    end do
	do i=1, ctl%m_bins
		if(abs(int(ctl%ini_weight_n(i)/ctl%n_basic)-(ctl%ini_weight_n(i)/ctl%n_basic))>1d-1)then
			print*, i, int(ctl%ini_weight_n(i)/ctl%n_basic),(ctl%ini_weight_n(i)/ctl%n_basic)
			print*, "error! set particle weighting an integer multiple of the minimium one"
			stop
		end if
	end do
    ctl%weight_n(1:ctl%m_bins)=int(ctl%ini_weight_n(1:ctl%m_bins)/ctl%n_basic)
	!do i=1, ctl%m_bins
	!	ctl%ini_Np(i)=mbh/ctl%ini_pre(i)*ctl%asymptot(1,i)
	!	print*,"i=", i, ctl%ini_Np(i)
	!end do
	!read(*,*)
    ctl%num_boundary_created=0
	ctl%num_boundary_elim=0
    ctl%num_clone_created=0
    
	call get_ccidx_from_type(ctl%idxstar,star_type_ms)
	call get_ccidx_from_type(ctl%idxsbh,star_type_bh)
	call get_ccidx_from_type(ctl%idxns, star_type_ns)
	call get_ccidx_from_type(ctl%idxwd, star_type_wd)
	call get_ccidx_from_type(ctl%idxbd, star_type_bd)


	ctl%nblock_size=int(ctl%grid_bins/ctl%ntasks)
    ctl%nblock_mpi_bg=ctl%nblock_size*rid+1
    ctl%nblock_mpi_ed=ctl%nblock_size*(rid+1)
end subroutine
subroutine init_model()
	!use com_main_gw
	implicit none
	call init_model_ctl()
    call init_chattery()
	call init_pro()
end subroutine
subroutine init_chattery()
    use com_main_gw
    implicit none
    character*(5) tmprid
    logical,save::first=.true.
    if(ctl%chattery.ge.1)then
        if(first)then
            first=.false.
        else
            close(unit=chattery_out_unit)
        end if
        chattery_out_unit=chattery_out_unit_0+rid
        write(unit=tmprid,fmt="(I5)") rid
        open(unit=chattery_out_unit,file="output/chattery_"//trim(adjustl(tmprid)), &
            status='replace')

    end if
    chattery_out_unit=6
end subroutine
