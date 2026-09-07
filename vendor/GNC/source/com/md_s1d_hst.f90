
module md_s1d_hst_basic_type
	!use my_intgl
	!use com_sts_type
	use md_s1d_basic_type
	integer,parameter:: f_log=0, f_linear=1
	integer,parameter::dct_x=1,dct_y=2
	
	type,extends(s1d_basic_type):: s1d_hst_basic_type
		integer, allocatable:: nb(:)
		real(8), allocatable:: fxw(:), nbw(:)
		real(8) nsw
        integer ns
		logical use_weight
        ! xb: the central position of each bin
        ! nb: unweighted number of samples per bin
		!nbw:   weighted number of samples per bin
		! fx:  unweighted distribution function of samples per pin
        ! fxw:   weighted distribution function of samples per pin
        ! nsw:      weighted total number of samples in all bins
        ! ns:      unweighted total number of samples in all bins

        contains
        procedure::print=>print_s1_hst
        procedure::get_s1d_hst
        procedure::get_s1d_hst_weight        
        generic::get_hst=>get_s1d_hst, get_s1d_hst_weight
	end type
	!interface get_hst
	!	module procedure get_hst_no_weight
	!	module procedure get_hst_weight
	!end interface
    private::print_s1_hst,get_s1d_hst,get_s1d_hst_weight
    
contains
	subroutine init_s1_hst_basic(this, xmin,xmax, n, use_weight)
		implicit none
		class(s1d_hst_basic_type)::this
		integer n
		real(8) xmin,xmax
		logical,optional:: use_weight
		if(n==0)then
			print*, "init_s1d:warnning s1d%nbin=0"
		end if
        !print*, "1"
		call init_s1d_basic(this, xmin,xmax,n, sts_type_dstr)
        !print*, "2"
		if(present(use_weight))then
			this%use_weight=use_weight
		else
			this%use_weight=.false.
		end if
        !print*, "3"
		if(allocated(this%nb))then
			deallocate(this%nb)
			if(this%use_weight)then
				deallocate(this%nbw, this%fxw)
			end if
		end if
		!print*, "4"
		allocate(this%nb(n))
		this%nb=0;this%ns=0
		if(this%use_weight)then
			allocate(this%nbw(n), this%fxw(n))
			this%nbw=0; this%fxw=0; this%nsw=0
		end if
		
	end subroutine

    subroutine print_s1_hst(this,str_)
        implicit none
        class(s1d_hst_basic_type)::this
        character*(*),optional ::str_
        character*(15) str_bin_type		
        integer i
        !print*, "?",str_,trim(adjustl(str_))
        !print*, present(str_)
		if(present(str_))then
        !print*, "0"
            write(*,*) "for sts s1d=", trim(adjustl(str_))
        end if
        !print*, "???"
		select case (this%bin_type)
        case (sts_type_grid)
            str_bin_type="GRID"
        case (sts_type_dstr)
            str_bin_type="DSTR"
        end select
        !print*, "1"
        
        !print*, "3"
		if(this%use_weight)then
			write(*, fmt="(10A15)") "bin_type", "nbin", "xmin", "xmax", "xstep", "NS", "NSW"
			write(*, fmt="(A15, I15, 3E15.5, I15, E15.5 )") str_bin_type, this%nbin, &
            this%xmin, this%xmax, this%xstep, this%ns, this%nsw
			write(unit=*, fmt="(20A20)") "X", "FX", "FXW", "NB", "NBW"
			do i=1, this%nbin
				write(unit=*, fmt="(3E20.10, I20, E20.10)") this%xb(i), this%fx(i), this%fxw(i), this%nb(i), this%nbw(i)
			end do
		else
			write(*, fmt="(10A15)") "bin_type", "nbin", "xmin", "xmax", "xstep", "NS"
			!print*, "2"
			write(*, fmt="(A15, I15, 3E15.5, I15 )") str_bin_type, this%nbin, &
            this%xmin, this%xmax, this%xstep, this%ns
			write(unit=*, fmt="(20A20)") "X", "FX", "NB"
			do i=1, this%nbin
				write(unit=*, fmt="(2E20.10, I20)") this%xb(i), this%fx(i), this%nb(i)
			end do
		end if
        !print*, "4"
    end subroutine
	
	
	subroutine get_s1d_hst_weight(s1_hst,x,w,n)
		implicit none
		class(s1d_hst_basic_type)::s1_hst
		!type(sts_fc_type)::fc
		integer n
		real(8) x(n),w(n)
		integer i, indx
		!call init_fc(fc,s1d%xmin,s1d%xmax,s1d%nbin,fc_spacing_linear)
		if(.not.s1_hst%use_weight)then
			print*, "error! s1_hst%use_weight=FALSE"
			stop
		end if
		call get_dstr_num_in_each_bin_weight(x,w,n,s1_hst%xmin,s1_hst%xstep, s1_hst%nbin, &
			s1_hst%nbw, s1_hst%nsw)
		s1_hst%fxw=s1_hst%nbw/s1_hst%xstep

		call get_dstr_num_in_each_bin(x(1:n),n,s1_hst%xmin,s1_hst%xstep, s1_hst%nbin, &
			s1_hst%nb, s1_hst%ns)
		s1_hst%fx=dble(s1_hst%nb)/s1_hst%xstep
	end subroutine
	subroutine get_s1d_hst(s1_hst,x,n)
		implicit none
		class(s1d_hst_basic_type)::s1_hst
		!type(sts_fc_type)::fc
		integer n
		real(8) x(n)!,w(n)
		integer i, indx
		!call init_fc(fc,s1d%xmin,s1d%xmax,s1d%nbin,fc_spacing_linear)
		if(s1_hst%use_weight)then
			print*, "warnning: s1_hst%use_weight=True, but there is no input weighting data "
		endif
		!w=1d0
        !print*, "start"

		call get_dstr_num_in_each_bin(x(1:n),n,s1_hst%xmin,s1_hst%xstep, s1_hst%nbin, &
			s1_hst%nb, s1_hst%ns)
		s1_hst%fx=dble(s1_hst%nb)/s1_hst%xstep
        !print*, "end"
	end subroutine	

	!subroutine output_sts_1d(s1d,fn)
	!	implicit none
	!	type(s1d_hst_basic_type)::s1d
	!	character*(*) fn
	!	integer i
	!	character*(8) tmp
!
	!	open(unit=999,file=trim(adjustl(fn)))
	!	do i=1, s1d%nbin
	!		write(unit=999,fmt="(3E20.10E4)") s1d%xb(i), s1d%fx(i), s1d%rnx(i)
	!	end do
	!	close(unit=999)
!
	!end subroutine
	subroutine output_s1d_hst(s1d_hst, fn)
        implicit none
        type(s1d_hst_basic_type)::s1d_hst
        character*(*) fn
        integer i
        open(unit=999,file=trim(adjustl(fn))//"_s1d_hst.txt")
		if(s1d_hst%use_weight)then
			write(unit=999,fmt="(10A27)") "xb", "fx", "fxw", "nb", "nbw"
			do i=1, s1d_hst%nbin
				!if(isnan(fc%fxw(i)).or.fc%ns.eq.0)then
				!	print*, "warnning: i, fx, fxw, ns=", i,fc%fx(i), fc%fxw(i), fc%ns
				!end if
				write(unit=999,fmt="(1P3E27.12E4, I27, 1PE27.12E4)") s1d_hst%xb(i), s1d_hst%fx(i), &
					s1d_hst%fxw(i), s1d_hst%nb(i), s1d_hst%nbw(i)
			end do
		else
			write(unit=999,fmt="(10A27)") "xb", "fx", "nb"
			do i=1, s1d_hst%nbin
				write(unit=999,fmt="(1P2E27.12E4, I27)") s1d_hst%xb(i), s1d_hst%fx(i), s1d_hst%nb(i)
			end do
		end if
		close(999)
    end subroutine
	
end module

module md_s1d_hst_type
	use md_s1d_hst_basic_type
	type,extends(s1d_hst_basic_type):: s1d_hst_type
	contains
		procedure::init=>init_s1d_hst
		procedure::read_s1d_hst
		procedure::write_s1d_hst
		generic :: read(unformatted) => read_s1d_hst
		generic :: write(unformatted) => write_s1d_hst
	end type
	private init_s1d_hst
	private::read_s1d_hst
	private::write_s1d_hst
contains 
	subroutine init_s1d_hst(this, xmin,xmax, n, use_weight)
		implicit none
		class(s1d_hst_type)::this
		integer n
		real(8) xmin,xmax
		logical,optional:: use_weight
		if(present(use_weight))then
			call init_s1_hst_basic(this,  xmin,xmax,n, use_weight)
		else
			call init_s1_hst_basic(this,  xmin,xmax,n, .false.)
		end if
	end subroutine

	subroutine read_s1d_hst(s1d,file_unit, iostat, iomsg)
		implicit none
		class(s1d_hst_type), intent(inout)::s1d
		integer,intent(in):: file_unit
		integer,intent(out)::iostat
		character(*), intent(inout) :: iomsg
		integer n
		read(unit=file_unit) s1d%nbin, s1d%xmin,s1d%xmax
		read(unit=file_unit) s1d%use_weight
		n=s1d%nbin
		call s1d%init(s1d%xmin,s1d%xmax,n,s1d%use_weight)
		read(unit=file_unit) s1d%xb(1:n), s1d%fx(1:n), s1d%nb(1:n), s1d%ns
		if(s1d%use_weight)then
			read(unit=file_unit)  s1d%nbw(1:n), s1d%fxw(1:n), s1d%nsw
		endif
	end subroutine
	subroutine write_s1d_hst(s1d,file_unit, iostat, iomsg)
		implicit none
		class(s1d_hst_type), intent(in)::s1d
		integer,intent(in):: file_unit
		integer,intent(out)::iostat
		character(*), intent(inout) :: iomsg
		integer n

		write(unit=file_unit) s1d%nbin, s1d%xmin,s1d%xmax
		write(unit=file_unit) s1d%use_weight
		n=s1d%nbin
		write(unit=file_unit) s1d%xb(1:n), s1d%fx(1:n), s1d%nb(1:n), s1d%ns
		if(s1d%use_weight)then
			write(unit=file_unit)  s1d%nbw(1:n), s1d%fxw(1:n), s1d%nsw
		end if
	end subroutine

end module
