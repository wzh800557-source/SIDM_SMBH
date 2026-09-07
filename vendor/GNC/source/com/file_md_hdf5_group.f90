module md_hdf5_group
	use hdf5
	use h5lt
	use h5tb
    type hdf5_group_type
    
		INTEGER(HID_T) :: group_id_upper, group_id      ! Group identifier		
		character*(200) group_name
	contains
		procedure::create=>hdf5_group_create
		procedure::close=>hdf5_group_close
	end type
	private::hdf5_group_create, hdf5_group_close
contains
	subroutine hdf5_group_create(this, group_id_upper, gname)
		implicit none
        class(hdf5_group_type)::this
        INTEGER(HID_T) :: group_id_upper
        character*(*) gname
        integer error
		call h5gcreate_f(group_id_upper, gname, this%group_id, error)
		this%group_name=trim(adjustl(gname))
	end subroutine
	subroutine hdf5_group_close(this)
		implicit none
        class(hdf5_group_type)::this
        integer error
		call h5gclose_f(this%group_id, error)
	end subroutine
end module
