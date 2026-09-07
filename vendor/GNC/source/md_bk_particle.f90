module md_bk_species
    use md_coeff
    use com_sts_type
	use md_particle_sample
	use md_chain
    implicit none
    
    type particle_samples_arr_type
        integer n 
        type(particle_sample_type),allocatable::sp(:)
        contains
			procedure::init=>init_particle_sample_arr
			!procedure::add_member=>add_member_particle_sample_arr
			procedure::output_bin=>output_particle_sams_arr_bin
			procedure::input_bin=>input_particle_sams_arr_bin
			procedure::select=>sams_arr_select_single
			procedure::output_txt=>output_particle_sams_txt
			!procedure::get_int_total
			!procedure::get_real_total
			!procedure::convert_to_int_arrays
			!procedure::convert_to_real_arrays
    end type	
	type samples_type_pointer
		type(chain_pointer_type),pointer::sp=>null()
		integer rid, index
	end type
	type samples_type_pointer_arr
		integer n
		type(samples_type_pointer),allocatable::pt(:)
		contains
			procedure::init=>init_pointer_arr
	end type


   private::init_particle_sample_arr
   private::input_particle_sams_arr_bin, output_particle_sams_arr_bin
   private::sams_arr_select_single,sams_selection_function
   private::output_particle_sams_txt
contains
subroutine init_pointer_arr(this,n)
	implicit none
	class(samples_type_pointer_arr)::this
	integer n
	if(allocated(this%pt))deallocate(this%pt)        
	allocate(this%pt(n))
	this%n=n
end subroutine

subroutine convert_sams_pointer_arr(sps, sps_pointer,type)
	!use com_main_gw
	implicit none	
	type(chain_type)::sps
	type(samples_type_pointer_arr)::sps_pointer
	type(chain_pointer_type),pointer::ps
	integer nsel,i, n
	integer,optional::type
	integer typeI

	nsel=0
	typeI=0
	if(present(type))then
		typeI=type
	end if
	call sps%get_length(nsel,type=typeI)
	call sps_pointer%init(nsel)
	nsel=0
	ps=>sps%head
	do while (associated(ps))
		select case(typeI)
		case(0)
			nsel=nsel+1
			sps_pointer%pt(nsel)%sp=>ps        
		case(1)
			select type(ca=>ps%ob)
			type is(particle_sample_type)
				nsel=nsel+1
				sps_pointer%pt(nsel)%sp=>ps        
			end select
		end select
		ps=>ps%next
	end do

end subroutine

subroutine init_particle_sample_arr(bksps,n)
    implicit none
    integer n,i
    class(particle_samples_arr_type)::bksps
    !print*, "start arr ini", n
    !print*, associated(bksps%sp)
    if(allocated(bksps%sp)) then
    !	print*, "deallocating arr...", bksps%n
        deallocate(bksps%sp)
    !	print*, "deallocate arr success!"
    end if
    !print*, "start allocation"
    allocate(bksps%sp(n))
    bksps%n=n
    !print*, "allocate success"
    do i=1, n
        call bksps%sp(i)%init()
    end do
    !print*, "finished init"
end subroutine
    subroutine sams_arr_select_condition_single(sps, sps_out,selection_func,ipar,rpar)
		implicit none	
		type(particle_samples_arr_type)::sps,sps_out
		integer nsel,i, exitflag, nhiar
		real(8) timebg,timeed
		logical::selection_func
		logical::selected
		integer,optional:: ipar(10)
		real(8),optional:: rpar(10)
		nsel=0
		do i=1, sps%n
			if(present(ipar))then
				selected=selection_func(sps%sp(i),ipar,rpar)
			else
				selected=selection_func(sps%sp(i))
			end if
			if(selected)	 nsel=nsel+1
		end do
		call sps_out%init( nsel)
		nsel=0
		do i=1, sps%n
			if(present(ipar))then
				selected=selection_func(sps%sp(i),ipar,rpar)
			else
				selected=selection_func(sps%sp(i))
			end if
			if(selected) then
				nsel=nsel+1
	!			PRINT*, "NSEL=",NSEL
				sps_out%sp(nsel)=sps%sp(i)
	!			PRINT*, SPS_OUT%SP(NSEL)%EXIT_FLAG, sps_out%sp(nsel)%agw, sps%sp(nsel)%agw
			end if		
		end do
	end subroutine
	subroutine copy_particle_sample_arr(scopy, sp)
		!scopy=>sp
		implicit none
		type(particle_samples_arr_type)::scopy, sp
		integer i
		call sp%init(scopy%n)
		do i=1, scopy%n
			sp%sp(i)=scopy%sp(i)
		end do
	end subroutine
    subroutine smmerge_arr_single(sma,n,smam)
		implicit none
		integer n,i,j, nsam
		type(particle_samples_arr_type)::sma(n), smam
		nsam=0
		do i=1, n
			nsam=nsam+sma(i)%n
		end do
		!print*, "nsam=",nsam
		call smam%init(nsam)
		nsam=0
		do i=1, n
			do j=1, sma(i)%n
				nsam=nsam+1
				!print*, nsam
				smam%sp(nsam)=sma(i)%sp(j)
			end do
		end do
	end subroutine

	subroutine set_sample_arr_indexs_rid_particle(sams_arr,rid)
		!use com_main_gw
		implicit none
		type(particle_samples_arr_type)::sams_arr
		integer i,rid
		do i=1, sams_arr%n
			sams_arr%sp(i)%idx=i
			sams_arr%sp(i)%rid=rid
		end do
	end subroutine
    subroutine output_particle_sams_arr_bin(bksps, fl)
        implicit none
        character*(*) fl
        integer i,n
        class(particle_samples_arr_type)::bksps
        open(unit=999,file=trim(adjustl(fl))//".bin", form='unformatted',access='stream')
        write(unit=999) bksps%n
        !print*, "output_sams_arr_bin:", sps%n
        do i=1, bksps%n
            call bksps%sp(i)%write_info(999)
        end do
        close(unit=999)
    end subroutine
    
    subroutine input_particle_sams_arr_bin(bksps, fl)
        implicit none
        character*(*) fl
        integer i,n
        class(particle_samples_arr_type)::bksps
        open(unit=999,file=trim(adjustl(fl))//".bin", form='unformatted',access='stream', status='old')
        read(unit=999) bksps%n
        !print*, "output_sams_arr_bin:", sps%n
        call bksps%init(bksps%n)
        do i=1, bksps%n
            call bksps%sp(i)%read_info(999)
        end do
        close(unit=999)
    end subroutine
   
	subroutine sams_arr_select_single(sps, sps_out, exitflag, timebg, timeed)
		implicit none
		class(particle_samples_arr_type)::sps, sps_out
		integer nsel,i, exitflag
		real(8) timebg,timeed
		!logical sams_selection_function
		nsel=0
		do i=1, sps%n
			if(sams_selection_function(sps,i,exitflag,timebg,timeed)) nsel=nsel+1
		end do
		!print*, "nsel=",nsel
		call sps_out%init(nsel)
		nsel=0
		do i=1, sps%n
			if(sams_selection_function(sps,i,exitflag,timebg,timeed)) then
				nsel=nsel+1
				sps_out%sp(nsel)=sps%sp(i)
			end if		
		end do
	!contains

	end subroutine
	logical function sams_selection_function(sps,i,exitflag,timebg,timeed)
		implicit none
		class(particle_samples_arr_type)::sps
		integer i,exitflag
		real(8) timebg,timeed
		sams_selection_function=.false.
		if(isnan(sps%sp(i)%weight_clone).or.isnan(sps%sp(i)%weight_N))then
			print*, "selection error:", sps%sp(i)%weight_clone, sps%sp(i)%weight_N, &
			sps%sp(i)%id
			stop
		end if
		!if(sps%sp(i)%id.eq.98172504)then
		!	print*, "98172504%exit_flag=",sps%sp(i)%exit_flag
		!end if
		!if(abs(sps%sp(i)%weight_real-sps%sp(i)%weight_asym*sps%sp(i)%weight_n*&
		!	sps%sp(i)%weight_clone)>0.1)then
		!	print*, "error!, id=", sps%sp(i)%id, sps%sp(i)%weight_real, &
		!	sps%sp(i)%weight_asym, sps%sp(i)%weight_n,  &
		!	sps%sp(i)%weight_clone
		!	stop
		!end if
		if(sps%sp(i)%exit_flag.eq.exitflag.or.(exitflag.eq.-1) &
			.or.(exitflag.eq.-2.and.sps%sp(i)%exit_flag.ne.exit_invtransit))then
			
			if(sps%sp(i)%exit_time>timebg.or.timebg<0)then
				if(sps%sp(i)%exit_time<timeed.or.timeed<0)then
					sams_selection_function=.true.
				end if
			end if
		
		end if
	end function

	subroutine output_particle_sams_txt(sps, fl)
		implicit none
		character*(*) fl
		integer i,n
		class(particle_samples_arr_type)::sps
		real(8) x,y,z
        real(8),external::rnd
        
		!type(particle_sample_type),pointer::sp
		open(unit=999,file=trim(adjustl(fl))//".txt")
		write(unit=999,fmt="(13A20, 3A10)") "m", "aout", "eout", "inc","om","pe", "rp", "en", "jm","weight", "x","y","z", &
			"exit_flag", "obtype"
		do i=1, sps%n
			associate(sp=>sps%sp(i))		
				sp%byot%an_in_mode=an_in_mode_mean
                sp%byot%me=rnd(0d0,2*pi)
				call by_em2st(sp%byot)
				call by_split_from_rd(sp%byot)
				write(unit=999,fmt="(13E20.10, 3I10)") sp%m, sp%byot%a_bin,  sp%byot%e_bin, &
			sp%byot%inc, sp%byot%om, sp%byot%pe ,sp%rp, sp%en, sp%jm,sp%weight_real,sp%byot%rd%x,&
                 sp%exit_flag, sp%obtype
			end associate
		end do
		close(unit=999)
	end subroutine

end module