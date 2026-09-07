
module md_s1d_type
	use md_s1d_basic_type
	type,extends(s1d_basic_type)::s1d_type
		integer type_int_size
		integer type_real_size
		integer type_log_size
		real(8),allocatable::pn(:)
	contains
		procedure::init=>init_s1d
		procedure,private::read_s1d
		procedure,private::write_s1d
		!procedure::conv_to_array=>conv_s1d_int_real_arrays
		!procedure::conv_from_array=>conv_int_real_arrays_s1d
#ifdef HDF5		
		procedure::save_hdf5=>save_s1d_hdf5
		procedure::read_hdf5=>read_s1d_hdf5
#endif		
		generic::read(unformatted)=>read_s1d
		generic::write(unformatted)=>write_s1d
	end type
	private init_s1d
	private::read_s1d
	private::write_s1d
	!private::conv_s1d_int_real_arrays
	!private::conv_int_real_arrays_s1d
#ifdef HDF5
	private::save_s1d_hdf5
	private::read_s1d_hdf5
#endif

contains
	subroutine init_s1d(this,  xmin,xmax,n, bin_type)
		implicit none
		class(s1d_type)::this
		integer n,bin_type
		real(8) xmin,xmax
		call init_s1d_basic(this,  xmin,xmax,n, bin_type)
		if(allocated(this%pn)) deallocate(this%pn)
		allocate(this%pn(this%nbin))

		this%type_int_size=2
		this%type_real_size=this%nbin*3+3
		this%type_log_size=1
	end subroutine
	subroutine read_s1d(s1d,file_unit, iostat, iomsg)
		implicit none
		class(s1d_type), intent(inout)::s1d
		integer,intent(in):: file_unit
		integer,intent(out)::iostat
		character(*), intent(inout) :: iomsg
		integer n, bin_type
		read(unit=file_unit) s1d%nbin, s1d%xmin,s1d%xmax, bin_type
		n=s1d%nbin
		call s1d%init(s1d%xmin,s1d%xmax,n,bin_type)
		read(unit=file_unit) s1d%xb(1:n), s1d%fx(1:n)
	end subroutine
	subroutine write_s1d(s1d,file_unit, iostat, iomsg)
		implicit none
		class(s1d_type), intent(in)::s1d
		integer,intent(in):: file_unit
		integer,intent(out)::iostat
		character(*), intent(inout) :: iomsg
		integer n

		write(unit=file_unit) s1d%nbin, s1d%xmin,s1d%xmax,s1d%bin_type
		n=s1d%nbin
		write(unit=file_unit) s1d%xb(1:n), s1d%fx(1:n)
	end subroutine
	subroutine conv_s1d_int_real_arrays(s1d, intarr, realarr,logarr)
		!use com_main_gw
		implicit none
		class(s1d_type)::s1d
		integer intarr(s1d%type_int_size)
		logical logarr(s1d%type_log_size)
		real(8) realarr(s1d%type_real_size)

		intarr(1:2)=(/s1d%nbin,s1d%bin_type/)
		logarr=(/s1d%is_spline_prepared/)
		realarr(1:s1d%nbin)=s1d%xb(1:s1d%nbin)
		realarr(s1d%nbin+1:2*s1d%nbin)=s1d%fx(1:s1d%nbin)
		if(s1d%is_spline_prepared)then
			realarr(s1d%nbin*2+1:3*s1d%nbin)=s1d%y2(1:s1d%nbin)
		end if
		realarr(s1d%nbin*3+1:3*s1d%nbin+3)=&
				(/s1d%xmin,s1d%xmax,s1d%xstep/)

	end subroutine
	subroutine conv_int_real_arrays_s1d(s1d, intarr,  realarr,logarr)
		!use com_main_gw
		implicit none
		class(s1d_type)::s1d
		integer nint, nreal
		integer intarr(s1d%type_int_size)
		real(8) realarr(s1d%type_real_size)
		logical logarr(s1d%type_log_size)
		
		s1d%nbin=intarr(1); s1d%bin_type=intarr(2)
		s1d%is_spline_prepared=logarr(1)
		s1d%xb(1:s1d%nbin)=realarr(1:s1d%nbin)
		s1d%fx(1:s1d%nbin)=realarr(1+s1d%nbin:s1d%nbin*2)
		if(s1d%is_spline_prepared)then
			if(.not.allocated(s1d%y2)) allocate(s1d%y2(s1d%nbin))
			s1d%y2(1:s1d%nbin)=realarr(1+s1d%nbin*2:s1d%nbin*3)
		end if
		s1d%xmin=realarr(s1d%nbin*3+1)
		s1d%xmax=realarr(s1d%nbin*3+2)
		s1d%xstep=realarr(s1d%nbin*3+3)
	end subroutine
#ifdef HDF5	
	subroutine save_s1d_hdf5(s1d, group_id,  tablename)
		use md_hdf5_table
		implicit none
		class(s1d_type)::s1d
		character*(*) tablename
		INTEGER(HID_T) :: group_id      ! Group identifier
        !INTEGER(HSIZE_T), PARAMETER :: nfields  = 11            ! nfields
        type(hdf5_table_type)::htable

		if(s1d%nbin.eq.0) return

		call htable%init_table(2, s1d%nbin,  tablename)        
		htable%field_names=(/"   X","  FX"/)	

		htable%field_types(1:2)=H5T_NATIVE_DOUBLE

		call htable%prepare_write_table(group_id)

		CALL htable%write_column_real(s1d%xb)
		CALL htable%write_column_real(s1d%fx)
	end subroutine

	subroutine read_s1d_hdf5(s1d, group_id, tablename)
		use md_hdf5_table
		implicit none
		class(s1d_type)::s1d
		character*(*) tablename
		INTEGER(HID_T) :: group_id      ! Group identifier
        !INTEGER(HSIZE_T), PARAMETER :: nfields  = 11            ! nfields
        type(hdf5_table_type)::htable

		if(s1d%nbin.eq.0) return
        call htable%init_table(2, s1d%nbin,  tablename)
        
        !nrecords=s1d%nbin;
        htable%field_names=(/"   X","  FX"/)	

		htable%field_types(1:2)=H5T_NATIVE_DOUBLE 

		call htable%prepare_read_table(group_id)

		CALL htable%read_column_real(s1d%xb)
        CALL htable%read_column_real(s1d%fx) 
		
	end subroutine
#endif
!==========================================================================================

end module
