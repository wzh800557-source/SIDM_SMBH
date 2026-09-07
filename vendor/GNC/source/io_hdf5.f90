#ifdef HDF5

subroutine save_stsfcweight_hdf5(x,w,n,group_id,tname,xmin,xmax,rxn,fc_flag)
	use hdf5
	use h5lt
	use com_main_gw
	implicit none
	character*(*) tname
    integer rxn, n, fc_flag
	INTEGER(HID_T) :: group_id      ! Group identifier
	real(8) x(n), w(n),xmin,xmax
	integer error
	type(sts_fc_type)::fc

	call fc%init(xmin,xmax, rxn, fc_flag,use_weight=.true.)
	call get_fc_weight(x,w,n, fc)
	call save_fc_hdf5( fc,group_id,tname)	
end subroutine
subroutine save_stsfcweight_auto_hdf5(x,w,n,group_id,groupname,rxn,fc_flag)
	use hdf5
	use h5lt
	use com_main_gw
	implicit none
	character*(*) groupname
    integer rxn, n, fc_flag
	INTEGER(HID_T) :: group_id      ! Group identifier
	real(8) x(n), w(n),xmin,xmax
	integer error
	type(sts_fc_type)::fc

	xmin=minval(x(1:n)); xmax=maxval(x(1:n))	

!	CALL h5gcreate_f(group_id, trim(adjustl(groupname)), subgroup_id, error)
	call fc%init(xmin,xmax, rxn, fc_flag,use_weight=.true.)
	call get_fc_weight(x,w,n, fc)
	call save_fc_hdf5(fc,group_id, trim(adjustl(groupname)))
end subroutine

subroutine input_dms_hdf5_pdf(dm,fl)
	use md_hdf5	
	use com_main_gw
	implicit none
	type(diffuse_mspec)::dm
	character*(*) fl
	integer i, error
	character*(4) tmpi
	type(hdf5_file_type)::hdf5_file

	type(hdf5_table_type)::hdf5_table
	INTEGER(HID_T) :: group_id,sub_group_id       ! File identifier	
    !print*, "000"
	call h5eset_auto_f(0,error)
	call hdf5_file%open(trim(adjustl(fl))//".hdf5",flag=hdf5_file_flag_read)
	if(hdf5_file%error.eq.0)then
		!print*, "dm%n=", dm%n
		do i=1, dm%n
			write(unit=tmpi, fmt="(I4)") i
			CALL h5gopen_f(hdf5_file%file_id, trim(adjustl(tmpi)), group_id, error)
			!print*, "star to read"
			call read_dms_hdf5_pdf(dm%mb(i)%all, group_id, "all")
			call read_dms_hdf5_pdf(dm%mb(i)%star, group_id, "star")
			call read_dms_hdf5_pdf(dm%mb(i)%bd, group_id, "bd")
			call read_dms_hdf5_pdf(dm%mb(i)%sbh, group_id, "sbh")
			call read_dms_hdf5_pdf(dm%mb(i)%wd, group_id, "wd")
			call read_dms_hdf5_pdf(dm%mb(i)%ns, group_id, "ns")
			
			CALL h5gclose_f(group_id, error)
		end do
		call read_dms_hdf5_pdf(dm%all%all, hdf5_file%file_id, "all")
		call read_dms_hdf5_pdf(dm%all%star, hdf5_file%file_id, "star")
		call read_dms_hdf5_pdf(dm%all%sbh, hdf5_file%file_id, "sbh")
		call read_dms_hdf5_pdf(dm%all%wd, hdf5_file%file_id, "wd")
		call read_dms_hdf5_pdf(dm%all%ns, hdf5_file%file_id, "ns")
		call read_dms_hdf5_pdf(dm%all%bd, hdf5_file%file_id, "bd")
		
	else
		print*, "error! file may not exist!"
		stop
	end if
	CALL hdf5_file%close()

end subroutine

subroutine output_dms_hdf5_pdf(dm, fl)
	use md_hdf5	
	use com_main_gw
	implicit none
	type(diffuse_mspec)::dm
	character*(*) fl
	integer i, error
	character*(4) tmpi
	type(hdf5_file_type)::hdf5_file
	type(hdf5_table_type)::hdf5_table
	INTEGER(HID_T) :: group_id,sub_group_id       ! File identifier	
    !print*, "000"
	call hdf5_file%open(trim(adjustl(fl))//".hdf5")
	!print*, "001"
	do i=1, dm%n
		write(unit=tmpi, fmt="(I4)") i
		CALL h5gcreate_f(hdf5_file%file_id, trim(adjustl(tmpi)), group_id, error)
		call save_dms_hdf5_pdf(dm%mb(i)%all, group_id, "all")
		call save_dms_hdf5_pdf(dm%mb(i)%star, group_id, "star")
		call save_dms_hdf5_pdf(dm%mb(i)%bd, group_id, "bd")
		call save_dms_hdf5_pdf(dm%mb(i)%sbh, group_id, "sbh")
		call save_dms_hdf5_pdf(dm%mb(i)%wd, group_id, "wd")
		call save_dms_hdf5_pdf(dm%mb(i)%ns, group_id, "ns")
		
		CALL h5gclose_f(group_id, error)
	end do
	call save_dms_hdf5_pdf(dm%all%all, hdf5_file%file_id, "all")
	call save_dms_hdf5_pdf(dm%all%star, hdf5_file%file_id, "star")
	call save_dms_hdf5_pdf(dm%all%sbh, hdf5_file%file_id, "sbh")
	call save_dms_hdf5_pdf(dm%all%wd, hdf5_file%file_id, "wd")
	call save_dms_hdf5_pdf(dm%all%ns, hdf5_file%file_id, "ns")
	call save_dms_hdf5_pdf(dm%all%bd, hdf5_file%file_id, "bd")

    call h5gcreate_f(hdf5_file%file_id, "dej", sub_group_id, error)
    call output_de_hdf5(dm,sub_group_id)
    call h5gclose_f(sub_group_id, error)
    !print*, "002"
	CALL hdf5_file%close()
    !print*, "003"
end subroutine
subroutine output_de_hdf5(dm,group_id)
    use md_hdf5
    use com_main_gw
    implicit none
    type(diffuse_mspec)::dm
    integer(HID_T)::  group_id
    type(s2d_type)::s2d
    integer i
    
    call dm%mb(1)%dc%s2_de_0%save_hdf5(group_id, "de_0_1")
        
	call s2d%init(dm%nbin_grid, dm%nbin_grid, dm%emin,dm%emax,dm%jmin,dm%jmax, sts_type_dstr)
	!call s2d%set_range()
	s2d%xcenter=dm%mb(1)%dc%s2_de_110%xcenter
	s2d%ycenter=dm%mb(1)%dc%s2_de_110%ycenter
	s2d%fxy=0
	s2d%fxy=dm%dc0%s2_de_0%fxy
	do i=1, dm%n
		s2d%fxy=s2d%fxy+dm%mb(1)%mc/dm%mb(i)%mc*dm%mb(i)%dc%s2_de_110%fxy
	end do
	call s2d%save_hdf5(group_id, "de1")
	call dm%mb(1)%dc%s2_de_110%save_hdf5(group_id, "de_110_1")

    if(dm%n>1)then

		call s2d%init(dm%nbin_grid, dm%nbin_grid, dm%emin,dm%emax,dm%jmin,dm%jmax, sts_type_dstr)
		!call s2d%set_range()
		s2d%xcenter=dm%mb(1)%dc%s2_de_110%xcenter
		s2d%ycenter=dm%mb(1)%dc%s2_de_110%ycenter
		s2d%fxy=0
		s2d%fxy=dm%dc0%s2_de_0%fxy
		do i=1, dm%n
			s2d%fxy=s2d%fxy+dm%mb(2)%mc/dm%mb(i)%mc*dm%mb(i)%dc%s2_de_110%fxy
		end do
		call s2d%save_hdf5(group_id, "de2")
        
		call dm%mb(2)%dc%s2_de_110%save_hdf5(group_id, "de_110_2")

    end if
    
	call dm%dc0%s2_dee%save_hdf5(group_id, "dee")
	call dm%dc0%s2_dej%save_hdf5(group_id, "dej")
	call dm%dc0%s2_djj%save_hdf5(group_id, "djj")

        s2d%fxy=dm%dc0%s2_dj_rest%fxy
        do i=1, dm%n
            s2d%fxy=s2d%fxy+(dm%mb(1)%mc+dm%mb(i)%mc)/dm%mb(i)%mc/2d0 &
                *dm%mb(i)%dc%s2_dj_111%fxy
        end do
		   
	call s2d%save_hdf5(group_id, "dj1")
     
    if(dm%n>1)then
        
            s2d%fxy=dm%dc0%s2_dj_rest%fxy
            do i=1, dm%n
                s2d%fxy=s2d%fxy+(dm%mb(2)%mc+dm%mb(i)%mc)/dm%mb(i)%mc/2d0 &
                    *dm%mb(i)%dc%s2_dj_111%fxy
            end do
            call s2d%save_hdf5(group_id, "dj2")             
    end if
end subroutine

subroutine read_dms_hdf5_pdf(so, group_id, str_)
	use md_hdf5	
	use com_main_gw
	implicit none
	type(dms_stellar_object)::so
	character*(*) str_
	integer i, error
	character*(4) tmpi
	!type(hdf5_file_type)::hdf5_file
	!type(hdf5_table_type)::hdf5_table
	INTEGER(HID_T) :: group_id,sub_group_id,s2d_group_id,lapl_id       ! File identifier	
	logical:: exist
	!if(so%n_real>0)then
		
		CALL h5gopen_f(group_id, trim(adjustL(str_)), sub_group_id, error)
		if(error.eq.0)then
			call so%fden%read_hdf5(sub_group_id, "fden")
			!call so%fden%print("fden")
			!read(*,*)
			call so%fna%read_hdf5(sub_group_id, "fNa")
			call so%barge%read_hdf5(sub_group_id,  "fgx")
			call so%fMa%read_hdf5(sub_group_id,  "fMa")
			!call read_fc_hdf5(sub_group_id, so%faniso, "fanosi")
			call so%fden_simu%read_hdf5(sub_group_id, "fden_simu")	
			
			!call h5gopen_f(sub_group_id, "gxj", s2d_group_id, error)    
			call so%gxj%read_hdf5(sub_group_id, "gxj")
			!call read_s2d_hdf5(s2d_group_id, so%gxj)
			!call so%gxj%print("gxj")
			!call h5gclose_f(s2d_group_id, error)
			so%n_real=1
			so%n=1
		else
			so%n_real=0
			so%n=0
			!print*, trim(adjustL(str_)), " not exist"
			!stop
		end if
		CALL h5gclose_f(sub_group_id, error)
	!end if
end subroutine

subroutine save_dms_hdf5_pdf(so, group_id, str_)
	use md_hdf5	
	use com_main_gw
	implicit none
	type(dms_stellar_object)::so
	character*(*) str_
	integer i, error
	character*(4) tmpi
	!type(hdf5_file_type)::hdf5_file
	!type(hdf5_table_type)::hdf5_table
	INTEGER(HID_T) :: group_id,sub_group_id,s2d_group_id       ! File identifier	
	if(so%n_real>0)then
		CALL h5gcreate_f(group_id, trim(adjustL(str_)), sub_group_id, error)
		call so%fden%save_hdf5(sub_group_id, "fden")
		call so%fNa%save_hdf5(sub_group_id, "fNa")
		call so%barge%save_hdf5(sub_group_id, "fgx")
		call so%fMa%save_hdf5(sub_group_id, "fMa")
		!call save_fc_hdf5(sub_group_id, so%faniso, "fanosi")
		!call so%faniso%print("faniso")
		call so%fden_simu%save_hdf5(sub_group_id, "fden_simu")	
		call so%gxj%save_hdf5(sub_group_id,"gxj")
		CALL h5gclose_f(sub_group_id, error)
	end if
end subroutine
 

#endif