module md_particle_sample
    use,intrinsic::ieee_arithmetic
    use md_binary
	use md_triple
    real(8)::log10clone_emax
    type track_type
        real(8) time, ac, ec, ain, ein
        real(8) Incin, incout, omin, omout
        integer state_flag, ms_star_type, mm_star_type
    end type
    type particle_sample_type
		integer id                    ! a unque id for the particle
		integer rid, idx  !rid=the proc id, idx=the index in array
		integer obtype, obidx
		integer state_flag_last         ! flag of the last dynamic, use to determine the type of exit_flag
		integer exit_flag
		integer length, length_to_expand
		integer track_step
		integer write_down_track, within_jt
		
		real(8) r_td, m, en0, jm0
		!type(particle)::p
		real(8) pd, rp, tgw, simu_bgtime
		real(8) En, Jm, create_time,exit_time
		real(8) djp, elp, den, djp0
		real(8) weight_clone  ! the weight factor due to clone
		real(8) weight_N      ! the weight factor due to number of particles set by simulation
							  ! weight_N=1 by default, and can be changed due to collisions.
		!real(8) weight_asym   ! the weight factor due to asymptotic boundary condition 
		real(8) weight_real   ! the real weight used for calculation
		type(binary) byot, byot_ini, byot_bf

		type(track_type),allocatable:: track(:)
		contains
		procedure::track_init
		procedure::init=>particle_sample_init
		procedure::read_info=>read_sample_info
		procedure::write_info=>write_sample_info
		procedure::print=>print_particle_sample
		!procedure::reset_to_ini=>reset_to_ini_particle
	end type
    integer,parameter::track_length_expand_block=100
	integer,parameter::nint_particle=10, nreal_particle=19

	integer,parameter::flag_ini_ini=1
    integer,parameter::flag_ini_or=5
    integer,parameter::record_track_nes=1, record_track_detail=2, record_track_all=3

	integer,parameter:: star_type_BH=1,star_type_NS=2,star_type_MS=3, star_type_WD=4
	integer,parameter:: star_type_BD=5  ! brown dwarf
	integer,parameter:: star_type_other=6

    integer,parameter::exit_normal=1
    integer,parameter::exit_other=100, exit_ejection=11

	integer,parameter::exit_by_gw_embhb=17, exit_emri_single=180, exit_plunge_single=181
    
    integer,parameter::exit_tidal=2
    integer,parameter::exit_merge_eject=200
    integer,parameter::exit_max_reach=3
	integer,parameter::exit_boundary_min=4, exit_boundary_max=5
	integer,parameter::exit_invtransit=7, exit_tidal_empty=8, exit_tidal_full=9,  exit_gw_iso=99

	integer,parameter:: state_ae_evl=1
	integer,parameter::state_emax=19, state_plunge=199
	integer,parameter::state_td=71

	private::particle_sample_init !, write_sample_info, read_sample_info
	private::print_particle_sample!, particle_sample_get_weight_clone
	!private::reset_to_ini_particle

contains
subroutine print_particle_sample(this, str)
	implicit none
	class(particle_sample_type)::this
	character*(*) str
	print*, str
	write(*,fmt="(A20,4I15)") "rid,idx,length=", & 
			this%rid, this%idx, this%length
	write(*,fmt="(A20,4I15)") "id, exit_flag=", & 
			this%id,this%EXIT_FLAG
    write(*,fmt="(A20,4I15)") "state_last=", this%state_flag_last
						
	write(*,fmt="(A20,4F15.6)") "otby:a,e, I, Om=", & 
			this%byot%a_bin, this%byot%e_bin, this%byot%Inc, this%byot%Om
	write(*,fmt="(A30,6E15.6)") "ctime, w(real,  clone, n)=", & 
			this%create_time, this%weight_real, this%weight_clone, &
            this%weight_n
end subroutine
character*(4) function star_type(type_int)
	implicit none
	integer type_int

	select case(type_int)
	case (star_type_BH)
		star_type="BH"
	case (star_type_MS)
		star_type="MS"
	case (star_type_NS)
		star_type="NS"
	case (star_type_WD)
		star_type="WD"
	case(star_type_BD)
		star_type="BD"            
	case default
		star_type='UNKNOWN'
		!   print*, "define star_type:",type_int
		!   stop        
	end select
end function

	subroutine track_init(sp,n)
		implicit none
		class(particle_sample_type)::sp
		integer n
		if(n>0)then
			if(allocated(sp%track))then
			deallocate(sp%track)
			endif
			allocate(sp%track(n))
		end if
		sp%length=0
		sp%length_to_expand=n
		!print*, "track init", n
	end subroutine
	subroutine particle_sample_init(sp)
		implicit none
		class(particle_sample_type) ::sp
		sp%write_down_track=0
		sp%within_jt=0
		sp%byot%a_bin=0; sp%byot%e_bin=0;
		sp%track_step=1; sp%pd=0d0
		sp%exit_flag=0;sp%m=0; sp%weight_real=-1d99
		sp%weight_clone=-1d99;sp%weight_N=-1d99
        sp%exit_time=0d0
		!if(allocated(sp%track)) deallocate(sp%track)
		!print*,"init:", sizeof(sp%track), allocated(sp%track)
		call track_init(sp,0)
	end subroutine

subroutine read_sample_info(sp,funit)
	integer funit
	integer mypos
	class(particle_sample_type)::sp
    logical,save::first=.true.

	call sp%init()
!	print*, sp%byot%a_bin
	read(unit=funit) sp%exit_time,sp%r_td,sp%m ,sp%en0,sp%jm0, &
                   sp%pd,sp%rp,sp%tgw, sp%djp0, sp%state_flag_last
	
	read(unit=funit) sp%obtype, sp%obidx, sp%rid, sp%idx, sp%id
!	print*, sp%obtype, sp%obidx
	read(unit=funit) sp%byot,sp%byot_ini, sp%byot_bf
	!sp%byot%a_bin=actmp; sp%byot%e_bin=ectmp
	!---------------------------------------------------------
	!print*, sp%byot%a_bin, sp%byin%a_bin
   ! inquire(unit=funit, pos=mypos)
	!call print_binary(sp%byin)
	!read(*,*)
   ! print*, "pos2=",mypos
	read(unit=funit) sp%length_to_expand,sp%exit_flag,sp%create_time,sp%Jm,sp%En, sp%den, &
	               sp%write_down_track,sp%track_step,sp%weight_clone, sp%djp, sp%elp, &
                    sp%within_jt, sp%length, sp%weight_real, &
					sp%weight_N
					if(isnan(sp%weight_real))then
						print*, "error! sp%weight_real = NaN", &
							sp%weight_real, sp%weight_N
						stop
					end if
					!print*, "read",sp%byot%a_bin
 	if(sp%length>0)Then
		!print*, "sp%length=",sp%length
		if(allocated(sp%track))deallocate(sp%track)
		allocate(sp%track(sp%length))
		read(unit=funit) sp%track(1:sp%length)
		!print*, "read: sizeoftrack:", sizeof(sp%track), sp%length
		!stop
	end if
	if(sp%byot%a_bin.eq.0d0.and.first)then
		write(*,*) "warnning!-------"
		write(*,*) sp%exit_flag
		print*,  sp%byot%a_bin,sp%byot%e_bin,sp%exit_time,sp%r_td,sp%m , sp%en0,sp%jm0,&
        sp%pd,sp%rp,sp%tgw, sp%djp0
        first=.false.
	end if
end subroutine
subroutine write_sample_info(sp,funit)
	implicit none
	integer funit
	class(particle_sample_type)::sp
    logical,save::first=.true.

	write(unit=funit) sp%exit_time,sp%r_td,sp%m ,sp%en0,sp%jm0, &
                   sp%pd,sp%rp,sp%tgw, sp%djp0, sp%state_flag_last
	
	write(unit=funit) sp%obtype, sp%obidx, sp%rid, sp%idx, sp%id

	write(unit=funit) sp%byot,sp%byot_ini, sp%byot_bf
	write(unit=funit) sp%length_to_expand,sp%exit_flag,sp%create_time,sp%Jm,sp%En, sp%den, &
	               sp%write_down_track,sp%track_step,sp%weight_clone, sp%djp, sp%elp, &
                    sp%within_jt, sp%length,  sp%weight_real, &
					sp%weight_N
					if(isnan(sp%weight_real))then
						print*, "error! sp%weight_real = NaN",  &
							sp%weight_real, sp%weight_N
						stop
					end if
					!print*, "write",sp%byot%a_bin
 	if(sp%length>0)Then
		!print*, "write sp: rid, idx, m",&
			 !sp%rid, sp%idx, sp%m, sp%exit_flag,sp%length
		!	 call sp%print("write")
		!stop
		print*, "sp%length=",sp%length
		write(unit=funit) sp%track(1:sp%length)
	end if
	!print*, "write: sizeoftrack:", sizeof(sp%track), sp%length
	!read(*,*)
	if(sp%byot%a_bin.eq.0d0.and.first)then
		write(*,*) "sg:warnning!-------"
		write(*,*) sp%exit_flag
		print*,  sp%byot%a_bin,sp%byot%e_bin,sp%exit_time,sp%r_td,sp%m, sp%id, &
        sp%en, sp%jm, sp%en0,sp%jm0,sp%pd,sp%rp,sp%tgw, sp%djp0
        first=.false.
	end if
end subroutine




end module