!subroutine sams_get_weight_single(sps,e0)
!	use com_main_gw
!	implicit none
!	type(particle_samples_arr_type)::sps
!	real(8) e0
!	integer i
!	do i=1, sps%n
!		call particle_sample_get_weight(sps%sp(i), e0)
!	end do
!	!print*, "3"
!end subroutine

subroutine sams_get_weight_clone_single(sps)
	use com_main_gw
	implicit none
	type(particle_samples_arr_type)::sps
	integer i, amplifier, obidx
    real(8) en
	integer mass_idx
	!print*, "sps%n=",sps%n
	do i=1, sps%n
        !print*, "i=",i
		obidx=sps%sp(i)%obidx
		call get_mass_idx(sps%sp(i)%m, mass_idx)
		amplifier=ctl%clone_factor(mass_idx)
		!call sps%sp(i)%get_weight_clone(ctl%clone_scheme, amplifier, ctl%clone_e0)
        
        if(sps%sp(i)%exit_flag.eq.exit_boundary_max)then
            en=-mbh/2d0/sps%sp(i)%byot_bf%a_bin
            !print*, sps%sp(i)%en/ctl%energy0, en/ctl%energy0
            !en=ctl%energy_max*0.999d0
        else
            en=SPS%sp(i)%en
        end if
        call particle_sample_get_weight_clone(en, ctl%clone_scheme, &
            amplifier,ctl%clone_e0,sps%sp(i)%weight_clone)
            !print*, "en, weight=", en, sps%sp(i)%weight_clone, sps%sp(i)%exit_flag, &
            !    sps%sp(i)%byot_bf%a_bin, log10(sps%sp(i)%en/ctl%energy0)
		if(isnan(sps%sp(i)%weight_clone))then
			print*, "sams_get_weight_clone_single:NaN", sps%sp(i)%weight_clone, &
				 sps%sp(i)%obtype, sps%sp(i)%obidx, sps%sp(i)%en, &
				sps%sp(i)%id
			stop
		end if
	end do
end subroutine

subroutine set_clone_weight(sms)
	use com_main_gw
	implicit none
	type(chain_type)::sms
	type(chain_pointer_type),pointer::pt
	integer i, obidx, amplifier
    real(8) en
	integer mass_idx

    pt=>sms%head
    do while(associated(pt).and.allocated(pt%ob))
        obidx=pt%ob%obidx
        select type(ca=>pt%ob)
        type is(particle_sample_type)
			call get_mass_idx(ca%m, mass_idx)
            amplifier=ctl%clone_factor(mass_idx)
        !type is(sample_type)
		!	call get_mass_idx(ca%m, mass_idx)
		!	amplifier=ctl%clone_factor(mass_idx)
        end select
       ! call pt%ob%get_weight_clone(ctl%clone_scheme, amplifier, ctl%clone_e0)
        !if(pt%ob%exit_flag.eq.exit_normal)then
            en=pt%ob%en
        !else
        !    en=-mbh/2d0/pt%ob%byot_bf%a_bin
        !end if
        call particle_sample_get_weight_clone(en, ctl%clone_scheme, &
        amplifier,ctl%clone_e0,pt%ob%weight_clone)
        pt=>pt%next
    end do		
    
end subroutine

subroutine sams_select_merge_single(sps, sps_out)
	use com_main_gw
	implicit none	
	type(chain_type)::sps, sps_out
	type(chain_pointer_type),pointer::ps,psout
	integer nsel,i, exitflag
	!print*, "sps%n=",sps%n
	!do i=1, sps%n
    call chain_select_by_condition(sps,sps_out, selection)
	!end do
contains
	logical function selection(pt)
		implicit none
		type(chain_pointer_type)::pt
		!print*,"sams_select_merge_single:",pt%ob%exit_flag
		!if(pt%ob%N_gene.ge.2)then
		!	print*, pt%ob%exit_flag, associated(pt%Append_left),  associated(pt%Append_right)
		!end if
		if(associated(pt%Append_left))then
			!.and.pt%ob%exit_flag.eq.exit_gw_capture_particle)then
			selection=.True.
		else
			selection=.False.
		end if
	end function
end subroutine

subroutine chain_select_type_single(sps, sps_out,obtype)
	use com_main_gw
	implicit none	
	type(chain_type)::sps, sps_out
	type(chain_pointer_type),pointer::ps,psout
	integer nsel,i, exitflag,obtype

    call chain_select_by_condition(sps,sps_out, selection)

contains
	logical function selection(pt)
		implicit none
		type(chain_pointer_type)::pt
		if(pt%ob%obtype.eq.obtype)then
			selection=.True.
		else
			selection=.False.
		end if
	end function
end subroutine

subroutine output_sams_merge_hierarchy(pt, fl)
	use com_main_gw
	implicit none
	character*(*) fl
	integer i
	type(chain_pointer_type),target::pt
	type(chain_pointer_type),pointer::psleft, psright, ps

	open(unit=999, file=trim(adjustl(fl))//"_hierarch.txt")
	ps=>pt
	write(999,fmt="(A15, A6, 20A20)" ) "", "nhiar", "m", "weight_real", "weight_clone", "weight_asym", &
		"weight_N", "ctime", "abin"
	do while(associated(ps))
		call write_txt("-", ps)
		ps=>ps%next
	end do
	close(999)
contains
	recursive subroutine write_txt(marker, pt_now)
		implicit none
		character*(*) marker
		character*(100) marker_here
		type(chain_pointer_type),target::pt_now

		write(unit=999,fmt="(A15,I6, 20F20.6)") trim(adjustl(marker)), pt_now%ob%m,  &
			pt_now%ob%weight_real, pt_now%ob%weight_clone, &
			pt_now%ob%weight_N, pt_now%ob%create_time, pt_now%ob%byot%a_bin
		marker_here=trim(adjustl(marker))//"-"
		if(associated(pt_now%append_left))then			
			call write_txt( marker_here, pt_now%append_left)
		end if
		if(associated(pt_now%append_right))then
			call write_txt( marker_here, pt_now%append_right)
		end if
	end subroutine
end subroutine

subroutine sams_arr_select_type_single(sps, sps_out, obtype)
	use com_main_gw
	implicit none
	type(particle_samples_arr_type)::sps,sps_out
	integer nsel,i, obtype
	nsel=0
	do i=1, sps%n
		if(selection()) nsel=nsel+1
	end do
	!print*, "nsel=",nsel
	call sps_out%init(nsel)
	nsel=0
	do i=1, sps%n
		if(selection()) then
			nsel=nsel+1
			sps_out%sp(nsel)=sps%sp(i)
		end if		
	end do
contains
	function selection()
		implicit none
		logical selection
		selection=.false.
		
		if(sps%sp(i)%obtype.eq.obtype)then
				selection=.true.
		end if
	end function
end subroutine	



subroutine output_all_sts(fout)
	use com_main_gw
	implicit none
	character*(*) fout
	integer i,flag

	!call output_events_single_bin(pteve, trim(adjustl(fout))//"events")
    if(ctl%include_loss_cone.ge.1)then
        flag=2
    else
        flag=1
    end if
    call output_eveset_txt(pteve_star, fout//"/pro/MS/",flag)
	if(ctl%idxsbh.ne.-1)then
    	call output_eveset_txt(pteve_sbh, fout//"/pro/BH/",flag)
	end if
	if(ctl%idxbd.ne.-1)then
		call output_eveset_txt(pteve_bd, fout//"/pro/BD/",flag)
	end if
	if(ctl%idxns.ne.-1)then
		call output_eveset_txt(pteve_ns, fout//"/pro/NS/",flag)
	end if
	if(ctl%idxwd.ne.-1)then
		call output_eveset_txt(pteve_wd, fout//"/pro/WD/",flag)
	end if
!
    
!	close(999)
end subroutine

subroutine get_memo_usage(pid)
	implicit none
	integer pid
	character*(200) system_command
	character*(10) pid_str
	write(unit=pid_str,fmt="(I10)") pid
	system_command="ps axo pid,rss | grep "//trim(adjustl(pid_str))
    print*, "pid, memory occupy(kB)"
	call system(trim(adjustl(system_command)))
end subroutine



subroutine get_ccidx_from_type(idx,bktype)
	use com_main_gw
	implicit none
	integer i
	integer idx, bktype
    idx=-1
    
	do i=1, ctl%num_bk_comp
        !print*, i, ctl%cc(i)%bktypemodel, bktype
		if(ctl%cc(i)%bktypemodel .eq.bktype)then	
			idx=i
			return
		end if
	end do

end subroutine

subroutine get_mass_idx(m,idx)
    use com_main_gw
    implicit none
    integer i, idx
    real(8) m
    idx=-1
    do i=1, ctl%m_bins
        if(m.ge.ctl%bin_mass_m1(i).and.m.le.ctl%bin_mass_m2(i))then
            idx=i
			!print*, m, idx
            return
        end if
    end do
    print*, "get_mass_idx idx=-1: m=", m
    !print*, "error! m=",m
    !stop
    
end subroutine


subroutine get_javg_coef(dm)
    use com_main_gw
    integer i, j, k,l
    real(8) ss
    type(diffuse_mspec)::dm
 !   integer n
 !   type(sts_fc_type)::s1_dee
 !   type(sts_fc_type)::s1_de
    
    do i=1, dm%n
        associate(s1_dee=>dm%mb(i)%dc%s1_dee, s1_de=>dm%mb(i)%dc%s1_de)
            call s1_dee%init(emin_factor,emax_factor, dm%mb(i)%dc%s2_de_0%nx, fc_spacing_linear)
            s1_dee%xb=dm%mb(i)%dc%s2_de_0%xcenter
            do j=1, dm%mb(i)%dc%s2_de_0%nx
                ss=0
                do k=1, dm%mb(i)%dc%s2_de_0%ny
                    do l=1, dm%n
                        associate(s2d=>dm%mb(l)%dc%s2_dee)
                            s2d%ystep=s2d%ycenter(2)-s2d%ycenter(1)
                            select case(dm%jbin_type)
                            case(jbin_type_lin)
                                ss=ss+s2d%fxy(j,k)*s2d%ycenter(k)*s2d%ystep
                            case(jbin_type_log)                                
                                ss=ss+s2d%fxy(j,k)*(10**s2d%ycenter(k))**2*s2d%ystep*log(10d0)
							case(jbin_type_sqr)
								ss=ss+s2d%fxy(j,k)/2d0*s2d%ystep
                            end select
                        end associate
                    end do
                    s1_dee%fx(j)=ss*2d0
                end do
            end do
			
            call s1_de%init(emin_factor,emax_factor, dm%mb(i)%dc%s2_de_0%nx, fc_spacing_linear)
            s1_de%xb=dm%mb(i)%dc%s2_de_0%xcenter
            do j=1, dm%mb(i)%dc%s2_de_0%nx
                ss=0
                do k=1, dm%mb(i)%dc%s2_de_0%ny
                    do l=1, dm%n
                        associate(s2dl1=>dm%mb(l)%dc%s2_de_0,s2dl2=>dm%mb(l)%dc%s2_de_110)
                            s2dl1%ystep=s2dl1%ycenter(2)-s2dl1%ycenter(1)
                            
                            select case(dm%jbin_type)
                            case(jbin_type_lin)
                                ss=ss+(s2dl1%fxy(j,k)+dm%mb(i)%mc/dm%mb(l)%mc*s2dl2%fxy(j,k))&
                                *s2dl1%ycenter(k)*s2dl1%ystep
                            case(jbin_type_log)      
                                ss=ss+(s2dl1%fxy(j,k)+dm%mb(i)%mc/dm%mb(l)%mc*s2dl2%fxy(j,k))&
                                *(10**s2dl1%ycenter(k))**2*s2dl1%ystep*log(10d0)
							case(jbin_type_sqr)
								ss=ss+(s2dl1%fxy(j,k)+dm%mb(i)%mc/dm%mb(l)%mc*s2dl2%fxy(j,k))&
                                /2d0*s2dl1%ystep
                            end select

                        end associate
                    end do
                    s1_de%fx(j)=ss*2d0
                end do
            end do   
        end associate     
    end do
end subroutine

subroutine print_javg_coef_theory_from_pow(xb,n, mc,m, b0,xmin,xmax,alpha,m1, res,res2)
    use com_main_gw
    integer i, n, m, j
    real(8) xb(n), b0(m), xmin, alpha(m), m1, aj,xmax
    real(8) sigma32, n0,kappa, de, mc(m), ss, const, res(n), res2(n)

    sigma32=(2*pi*ctl%v0**2)**(-3/2d0)
    n0=ctl%n0
!    kappa=(4*pi*mc)**2*log(mbh)
    const=16*pi**2*log(mbh)*sigma32*n0
    
    do i=1, n
        ss=0
        do j=1, m
            aj=b0(j)/xmin**alpha(j)
            ss=ss+(m1*mc(j)*aj/(alpha(j)-1.5d0)*(xmax**(alpha(j)-1.5d0)-xb(i)**(alpha(j)-1.5d0))-&
              mc(j)**2*(aj/(alpha(j)+1)/xb(i)**2.5d0*(xb(i)**(alpha(j)+1)-xmin**(alpha(j)+1))+b0(j)/xb(i)**2.5d0))
        end do
        res(i)=const*ss*xb(i)**1.5d0

        ss=0
        do j=1, m
            aj=b0(j)/xmin**alpha(j)
            !print*, alpha(j).eq.0.5d0, alpha(j)
            if(alpha(j).ne.0.5d0)then
            ss=ss+mc(j)**2*(aj/(alpha(j)-0.5d0)*(xmax**(alpha(j)-0.5d0)-xb(i)**(alpha(j)-0.5d0))+&
              (aj/(alpha(j)+1)/(xb(i)**1.5d0)*(xb(i)**(alpha(j)+1)-xmin**(alpha(j)+1))+b0(j)/xb(i)**1.5d0))
            else
            ss=ss+mc(j)**2*(aj*log(xmax/xb(i))+&
                (aj/(alpha(j)+1)/(xb(i)**1.5d0)*(xb(i)**(alpha(j)+1)-xmin**(alpha(j)+1))+b0(j)/xb(i)**1.5d0))    
            end if
        end do
        res2(i)=4d0/3d0*const*ss*xb(i)**0.5d0

        write(unit=2,fmt=*) xb(i), res(i), res2(i)
    end do
end subroutine


subroutine get_boundary_flux(en0, en1,sample)
    use com_main_gw
    implicit none
    real(8) en0, en1
    type(particle_sample_type)::sample
    integer idxmass, idxob
    if(en0>ctl%energy_boundary.and.en1<ctl%energy_boundary)then
        call get_mass_idx(sample%m, idxmass)
        call get_type_idx(sample,idxob)
        if(idxob.ne.-1)then
            ctl%bin_mass_flux_in(idxmass,idxob)=ctl%bin_mass_flux_in(idxmass,idxob)+sample%weight_n
        end if
        return
    end if
   
    if(en1>ctl%energy_boundary.and.en0<ctl%energy_boundary.and.en1<ctl%energy_min)then
        call get_mass_idx(sample%m, idxmass)
        call get_type_idx(sample,idxob)
        if(idxob.ne.-1)then
            ctl%bin_mass_flux_out(idxmass,idxob)=ctl%bin_mass_flux_out(idxmass,idxob)+sample%weight_n
        end if
        !print*, "|==<-|"
        return
    end if
    if(en1<ctl%energy_max)then
        call get_mass_idx(sample%m, idxmass)
        call get_type_idx(sample,idxob)
        if(idxob.ne.-1)then
            ctl%bin_mass_emax_out(idxmass,idxob)=ctl%bin_mass_emax_out(idxmass,idxob)+sample%weight_n
        end if
    end if
end subroutine

subroutine get_num_at_boundary(bksam,bysam)
    use com_main_gw
    implicit none
    type(chain_type)::bksam, bysam
    type(chain_pointer_type),pointer::pt
    integer obmidx,obtidx,jmidx,typeidx

    ctl%bin_mass_N=0
    pt=>bksam%head
    do while(associated(pt))
        if(pt%ob%en>ctl%energy_boundary.and.pt%ob%en<ctl%energy_min&
            .and.pt%ob%exit_flag.eq.exit_normal)then
            call get_mass_idx(pt%ob%m, obmidx)
            call get_type_idx(pt%ob,obtidx)
            if(obtidx.ne.-1)then
                ctl%bin_mass_N(obmidx,obtidx)=ctl%bin_mass_N(obmidx,obtidx)+pt%ob%weight_n
            end if
        end if
        pt=>pt%next
    end do
    pt=>bysam%head
    do while(associated(pt))
        if(allocated(pt%ob).and.pt%ob%en>ctl%energy_boundary.and.pt%ob%en<ctl%energy_min&
            .and.pt%ob%exit_flag.eq.exit_normal)then
            call get_mass_idx(pt%ob%m, obmidx)
            call get_type_idx(pt%ob,obtidx)
            if(obtidx.ne.-1)then
                ctl%bin_mass_N(obmidx,obtidx)=ctl%bin_mass_N(obmidx,obtidx)+pt%ob%weight_n
            end if
        end if
        pt=>pt%next
    end do

    
end subroutine
subroutine print_bin_mass_N()
    use com_main_gw
    implicit none
    integer i,j
    write(chattery_out_unit,fmt=*) "BIN_MASS_N:"
    write(chattery_out_unit,fmt="(30A12)") "id", "mI", "star", "sbh", "bbh"
    do i=1, ctl%m_bins
        do j=1, ctl%num_bk_comp
            if(.not.ctl%bin_mass_N(i,j).eq.0)then
           ! write(*, fmt="(A8, I8)") "COMP=", j
            write(chattery_out_unit,fmt="(A8, 2A12, 30I12.3)") "BM", rid, i, ctl%bin_mass_N(i,1:3)!, &
     !       ctl%bin_mass_Nbalance_ot(i,1:ctl%num_bk_comp+ctl%by_sp_num)
            end if
        end do
    end do
end subroutine

subroutine print_boundary_sts()
    use com_main_gw
    implicit none
    boundary_cross_net=boundary_cross_up-boundary_cross_down
    write(chattery_out_unit,fmt="(8A10)") "emax_cros", "emin_dir", "emin_cros", "replace", "create", &
        "up", "down","net"
    write(chattery_out_unit,fmt="(8I10)") boundary_sts_emax_cros, &
        boundary_sts_emin_dir, boundary_sts_emin_cros, &
        boundary_sts_replace, boundary_sts_create, boundary_cross_up, &
        boundary_cross_down, boundary_cross_net

    boundary_sts_emin_cros=0; boundary_sts_replace=0; 
    boundary_sts_create=0; boundary_sts_emin_dir=0;
    boundary_cross_up=0; boundary_cross_down=0; boundary_sts_emax_cros=0
end subroutine


subroutine check(str_pos)
    use com_main_gw
    implicit none
    type(chain_pointer_type),pointer::pt
    character*(*) str_pos
    pt=>bksams%head
    do while(associated(pt))
        if(pt%ob%en<ctl%energy_min&
        .and.pt%ob%exit_flag.eq.exit_normal)then
            print*, str_pos
            print*, "error!", pt%ob%en/ctl%clone_e0, ctl%clone_e0, pt%ob%id
            stop
        end if
        pt=>pt%next
    end do

end subroutine

subroutine get_star_type(type_int, str)
    use com_main_gw
    implicit none
    integer type_int
    character*(*) str
    select case(type_int)
    case (star_type_BH)
        str="BH"
    case (star_type_MS)
        str="MS"
    case (star_type_NS)
        str="NS"
    case (star_type_WD)
        str="WD"
    case (star_type_BD)
        str="BD"		
    case default
        str='UNKNOWN'
      !  print*, "unknown type type_int=",type_int
     !   stop
     !   print*, "define star_type:",type_int
     !   stop        
    end select
end subroutine

subroutine deallocate_chains_arrs()
    use com_main_gw
    implicit none
    integer i,j
    !print*, "1"
    call bksams%destory()
    if(allocated(bksams_arr%sp))deallocate(bksams_arr%sp)
    if(allocated(bksams_arr_norm%sp))deallocate(bksams_arr_norm%sp)
    if(allocated(bksams_pointer_arr%pt))deallocate(bksams_pointer_arr%pt)
    !print*, "4"
    do i=1, dms%n
        if(allocated(dms%mb(i)%bstar%nejw)) deallocate(dms%mb(i)%bstar%nejw)
        if(allocated(dms%mb(i)%bbh%nejw)) deallocate(dms%mb(i)%bbh%nejw)
        if(allocated(dms%mb(i)%sbh%nejw)) deallocate(dms%mb(i)%sbh%nejw)
        if(allocated(dms%mb(i)%bbh%nejw)) deallocate(dms%mb(i)%bbh%nejw)
        if(allocated(dms%mb(i)%all%nejw)) deallocate(dms%mb(i)%all%nejw)
		if(allocated(dms%mb(i)%bd%nejw)) deallocate(dms%mb(i)%bd%nejw)
		if(allocated(dms%mb(i)%wd%nejw)) deallocate(dms%mb(i)%wd%nejw)
		if(allocated(dms%mb(i)%ns%nejw)) deallocate(dms%mb(i)%ns%nejw)
        call dms%mb(i)%bstar%deallocation()
        call dms%mb(i)%star%deallocation()
        call dms%mb(i)%bbh%deallocation()
        call dms%mb(i)%sbh%deallocation()
        call dms%mb(i)%all%deallocation()
		call dms%mb(i)%ns%deallocation()
		call dms%mb(i)%wd%deallocation()
		call dms%mb(i)%bd%deallocation()
		do j=1, n_tot_comp
			Nullify(dms%mb(i)%dsp(j)%p)
		end do	
    end do
    deallocate(dms%mb)
   ! print*, "5"
    if(allocated(cfs%cfs_110))then
        deallocate(cfs%cfs_110,cfs%cfs_111,cfs%cfs_131,cfs%cfs_13_1, &
        cfs%cfs_130,cfs%cfs_330,cfs%cfs_310, cfs%jum,cfs%s)
    !    this%cfs_010,this%cfs_0_11)
    end if
    if(allocated(df)) deallocate(df)
end subroutine