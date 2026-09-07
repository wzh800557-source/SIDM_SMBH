!module md_sts
!	type(stats_1d_type)
!		real(8),allocatable::bin(:),fprb(:),fcmp(:),fncp(:)
!		integer,allocatable::ffra(:),fcfr(:),fncf(:)
!		real(8) mean,std
!	end type
!end module

module md_s2_hst_type
	!use md_dm2
	use md_s1d_hst_type
	use md_s2d_basic_type
	use md_fc_type
	type,extends(s2d_basic_type):: s2d_hst_basic_type
        !type (s1d_type),allocatable::cy(:),cx(:)
		!type (s1d_type):: hst_x,hst_y   !not yet implimented
        !type(dm2),allocatable::da_x_proj(:),da_y_proj(:)
        !type(sts_fc_type),allocatable::fc_x_proj(:), fc_y_proj(:)
		!real(8),allocatable:: xcenter(:),ycenter(:),fxy(:,:)
	!	real(8),allocatable::fxbg(:,:), fybg(:,:)
		!  x, y: the value of the coordinates, f(x,y): the joint probability
		real(8),allocatable:: nxyw(:,:)
		integer,allocatable:: nxy(:,:)
		real(8),allocatable:: fxyw(:,:)
		! n(x,y) the number in each bins
		real(8),allocatable::xmean(:),ymean(:)
		real(8),allocatable::xsct(:),ysct(:)
	    !real(8) xmin,xmax,ymin,ymax,xstep,ystep
        !integer nx,ny
	!	integer intrp_scale
		logical use_weight
		contains
		procedure::get_stats_weight=>get_stats_2d_weight
		procedure::print=>print_s2d
		!procedure::set_range=>set_range_s2d
        procedure::output_txt=>output_sts_2d
#ifdef HDF5
		!procedure::save_hdf5=>save_hdf5_s2d
		!procedure::read_hdf5=>read_hdf5_s2d
#endif
	end type
    private::output_sts_2d
!#ifdef HDF5
    !private::save_hdf5_s2d,read_hdf5_s2d
!#endif
contains
	subroutine init_s2d_hst_basic(s2d,nx,ny,xmin,xmax,ymin,ymax,use_weight)
		implicit none
		class(s2d_hst_basic_type)::s2d
		integer nx,ny,bin_type
        integer i
		real(8) xmin,xmax,ymin,ymax
		logical,optional::use_weight

		if(nx.eq.0.or.ny.eq.0)then
			print*, "init_s2d: warnning s2d%nx=0 or s2d%ny=0"
		end if
		call init_s2d_basic(s2d,nx,ny,xmin,xmax,ymin,ymax,sts_type_dstr)
		if(present(use_weight))then
			s2d%use_weight=use_weight
		else
			s2d%use_weight=.false.
		end if
		if (allocated(s2d%nxy))then
			!print*, "deallocating..."
			deallocate(s2d%nxy)
			deallocate(s2d%xmean,s2d%ymean)
			!deallocate(s2d%da_x_proj,s2d%da_y_proj)
			!deallocate(s2d%fc_x_proj,s2d%fc_y_proj)
			deallocate(s2d%xsct,s2d%ysct )
			if(s2d%use_weight)then
				deallocate(s2d%fxyw,s2d%nxyw)
			end if
			!deallocate(s2d%cy, s2d%cx, s2d%fx_proj, s2d%fy_proj)
		end if
		
		allocate(s2d%nxy(s2d%nx,s2d%ny))
		!allocate(s2d%da_y_proj(s2d%nx),s2d%da_x_proj(s2d%ny))
		!allocate(s2d%fc_y_proj(s2d%nx),s2d%fc_x_proj(s2d%ny))
		allocate(s2d%ymean(s2d%nx),s2d%xmean(s2d%ny))
        !allocate(s2d%fx_proj(s2d%nx), s2d%fy_proj(s2d%ny))
		allocate(s2d%xsct(s2d%ny),s2d%ysct(s2d%nx))

		s2d%nxy=0

		if(s2d%use_weight)then
			allocate(s2d%nxyw(s2d%nx,s2d%ny), s2d%fxyw(s2d%nx,s2d%ny) )
			s2d%fxyw=0;s2d%nxyw=0
		end if
		!allocate(s2d%cy(s2d%ny),s2d%cx(s2d%nx))
		
        !do i=1, s2d%nx
        !    call s2d%fx_proj(i)%init(s2d%xmin,s2d%xmax,s2d%nx,fc_spacing_linear)
        !end do
        !do i=1, s2d%ny
        !    call s2d%fy_proj(i)%init(s2d%ymin,s2d%ymax,s2d%ny,fc_spacing_linear)
        !end do
        
		
	end subroutine

	subroutine print_s2d(this,str_)
		implicit none
        class(s2d_hst_basic_type)::this
        character*(15) str_bin_type
		character*(*) , optional::str_
        integer i
		if(present(str_))then
            write(*,*) "s2d=", trim(adjustl(str_))
        end if

        write(unit=*, fmt="(10A15)")  "nx", "xmin", "xmax", "xstep", &
			 "ny", "ymin", "ymax", "ystep"
        write(unit=*, fmt="(A15, 2(I15, 3E15.5))")  this%nx,&
            this%xmin, this%xmax, this%xstep,  this%ny, this%ymin,this%ymax, this%ystep
        
		do i=1, this%nx
			write(*, fmt="(100E12.3)") this%fxy(i,:)
		end do
    end subroutine




	subroutine norm_s2d(s2d)
		implicit none
		type(s2d_hst_basic_type)::s2d
		integer i,j
		real(8) sums
		sums=0
		do i=1, s2d%nx
			do j=1, s2d%ny
				sums=sums+s2d%fxy(i,j)
			end do
		end do
		if(sums.ne.0)then
			s2d%fxy=s2d%fxy/sums
		end if
	end subroutine

	subroutine get_stats_2d(x,y,n,s2d)

		implicit none
		type(s2d_hst_basic_type)::s2d
		integer n,rxn,i
		real(8) x(n),y(n)
		real(8) xmin,xmax,ymin,ymax
		!call d2to1(x,y,n,s2d%xmin,s2d%xmax,s2d%nx, s2d%ymin,s2d%ymax,s2d%da_y_proj)
		!do i=1, s2d%nx
		!	call arravg(s2d%da_y_proj(i)%y,s2d%da_y_proj(i)%n,s2d%ymean(i))
		!end do
		call cal_bin2_arr(x,y,n,s2d%xmin,s2d%xmax,s2d%nx,s2d%ymin,s2d%ymax,s2d%ny,&
					s2d%xcenter,s2d%ycenter,s2d%fxy,sts_type_dstr,1)
		call cal_bin2_arr(x,y,n,s2d%xmin,s2d%xmax,s2d%nx,s2d%ymin,s2d%ymax,s2d%ny,&
					s2d%xcenter,s2d%ycenter,s2d%fxyw,sts_type_dstr,0)

	end subroutine

	subroutine get_stats_2d_weight(s2d, x,y,w,n)
		implicit none
		class(s2d_hst_basic_type)::s2d
		integer n,rxn,i, ns
		real(8) x(n),y(n),w(n)
		!logical,optional::get_data_in_bins
		!real(8) xmin,xmax,ymin,ymax

		
		!print*,"xcenter=",s2d%xcenter
		!print*, "ycenter=", s2d%ycenter
		
		call cal_bin2_arr_weight(x,y,w,n,s2d%xmin,s2d%xmax,s2d%nx,s2d%ymin,s2d%ymax,s2d%ny,&
					s2d%xcenter,s2d%ycenter,s2d%fxyw,sts_type_dstr,2)
		!print*, "s2d%fxy(1,1)=",isnan(s2d%fxy(1,1)),s2d%fxy(1,1)
		call cal_bin2_arr_weight(x,y,w,n,s2d%xmin,s2d%xmax,s2d%nx,s2d%ymin,s2d%ymax,s2d%ny,&
					s2d%xcenter,s2d%ycenter,s2d%nxyw,sts_type_dstr,0)
		!print*, "333"
	end subroutine

	subroutine get_stats_2d_xyzweight(x,y,z,w,n,s2d)
		implicit none
		class(s2d_hst_basic_type)::s2d
		integer n,i,j,idx,idy
		real(8) x(n),y(n),z(n),w(n)
		s2d%xstep=(s2d%xmax-s2d%xmin)/real(s2d%nx)
		s2d%ystep=(s2d%ymax-s2d%ymin)/real(s2d%ny)
		do i=1, s2d%nx
			s2d%xcenter(i)=s2d%xmin+s2d%xstep*(i-0.5d0)
		end do
		do i=1, s2d%ny
			s2d%ycenter(i)=s2d%ymin+s2d%ystep*(i-0.5d0)
		end do
		do i=1, n
			if(x(i)>s2d%xmin.and.x(i)<s2d%xmax.and.y(i)>s2d%ymin.and.y(i)<s2d%ymax)then
				select case (s2d%bin_type)
				case (sts_type_grid)
					idx=(x(i)-s2d%xmin)/s2d%xstep+2
					idy=(y(i)-s2d%ymin)/s2d%ystep+2
				case (sts_type_dstr)
					idx=int((x(i)-s2d%xmin)/s2d%xstep+1)
					idy=int((y(i)-s2d%ymin)/s2d%ystep+1)
				end select
				if(idx>=0.and.idy>=0.and.idx<=s2d%nx.and.idy<=s2d%ny)then
					s2d%fxy(idx,idy)=s2d%fxy(idx,idy)+w(i)*z(i)
					s2d%fxyw(idx,idy)=s2d%fxyw(idx,idy)+w(i)
					s2d%nxy(idx,idy)=s2d%nxy(idx,idy)+1
				end if
			end if
		end do
		do i=1, s2d%nx
			do j=1,s2d%ny
				if(s2d%nxy(i,j).ne.0)then
					s2d%fxy(i,j)=s2d%fxy(i,j)/s2d%fxyw(i,j)	
				end if
			end do
		end do
	end subroutine

	subroutine output_sts_2d(s2d,fn)
		implicit none
		class(s2d_hst_basic_type)::s2d
		character*(*) fn
		integer i
		character*(8) tmp
        if(.not.allocated(s2d%xcenter))return
		open(unit=999,file=trim(adjustl(fn))//"_fxy.txt")
		write(unit=tmp,fmt="(I4)") s2d%ny
		do i=1, s2d%nx
			write(unit=999,fmt="("//trim(adjustl(tmp))//"E20.10E4)") s2d%fxy(:,i)
		end do
		close(unit=999)

		!open(unit=999,file=trim(adjustl(fn))//"_rnxy.txt")
		!write(unit=tmp,fmt="(I4)") s2d%ny
		!do i=1, s2d%nx
		!	write(unit=999,fmt="("//trim(adjustl(tmp))//"E20.10E4)") s2d%fxyw(:,i)
		!end do
		!close(unit=999)

		open(unit=999,file=trim(adjustl(fn))//"_fxy_x.txt")
		do i=1, s2d%nx
			write(unit=999,fmt="(E20.10)") s2d%xcenter(i)
		end do
		close(unit=999)

		open(unit=999,file=trim(adjustl(fn))//"_fxy_y.txt")
		do i=1, s2d%ny
			write(unit=999,fmt="(E20.10)") s2d%ycenter(i)
		end do
		close(unit=999)

		!open(unit=999,file=trim(adjustl(fn))//"_2dymean.txt")
		!do i=1, s2d%nx
		!	write(unit=999,fmt="(4E20.10)") s2d%da_y_proj(i)%center, s2d%ymean(i)
		!end do
		!close(unit=999)

	end subroutine


end module



module md_s2d_type
    use md_s2_hst_type
	use bin_constants
	type,extends(s2d_basic_type) ::s2d_type
	contains
			procedure::init=>init_s2d
			procedure::read_s2d
			procedure::write_s2d
#ifdef HDF5		
			procedure::save_hdf5=>save_s2d_hdf5
			procedure::read_hdf5=>read_s2d_hdf5
			procedure::output_hdf5=>output_s2d_hdf5
#endif								
			generic ::read(unformatted)=>read_s2d
			generic ::write(unformatted)=>write_s2d
	end type
	private::init_s2d
#ifdef HDF5
	private::save_s2d_hdf5, write_s2d
	private::read_s2d_hdf5, read_s2d
#endif	
!=============================================================
	type,extends(s2d_hst_basic_type)::s2d_hst_type
	contains
			procedure::init=>init_s2d_hst
			procedure::read_s2d_hst
			procedure::write_s2d_hst
#ifdef HDF5		
			procedure::save_hdf5=>save_s2d_hst_hdf5
			procedure::read_hdf5=>read_s2d_hst_hdf5
#endif								
			generic ::read(unformatted)=>read_s2d_hst
			generic ::write(unformatted)=>write_s2d_hst
	end type

	private::init_s2d_hst
	private::read_s2d_hst,write_s2d_hst
#ifdef HDF5
	private::save_s2d_Hst_hdf5
	private::read_s2d_Hst_hdf5
#endif	
contains
	subroutine init_s2d(s2d,nx,ny,xmin,xmax,ymin,ymax,bin_type)
		implicit none
		class(s2d_type)::s2d
		integer nx,ny,bin_type
		integer i
		real(8) xmin,xmax,ymin,ymax
		
		call init_s2d_basic(s2d,nx,ny,xmin,xmax,ymin,ymax,bin_type)
	end subroutine

	subroutine read_s2d(s2d,file_unit, iostat, iomsg)
		implicit none
		class(s2d_type),intent(inout)::s2d
		integer nx,ny,i, bin_type
		integer,intent(in):: file_unit
		integer, intent(out) :: iostat
		character(*), intent(inout) :: iomsg

		read(unit=file_unit) nx, ny,s2d%xmin,s2d%xmax,s2d%ymin,s2d%ymax, bin_type
		call s2d%init(nx,ny,s2d%xmin,s2d%xmax,s2d%ymin,s2d%ymax, bin_type)
		read(unit=file_unit) s2d%xcenter(1:nx), s2d%ycenter(1:ny), s2d%fxy(1:nx,1:ny)
		read(unit=file_unit) s2d%xstep, s2d%ystep
	end subroutine

	subroutine write_s2d(s2d,file_unit, iostat, iomsg)
		implicit none
		class(s2d_type),intent(in)::s2d
		integer nx,ny,i
		integer,intent(in):: file_unit
		integer, intent(out) :: iostat
		character(*), intent(inout) :: iomsg

		write(unit=file_unit) s2d%nx,s2d%ny,s2d%xmin,s2d%xmax,s2d%ymin,s2d%ymax, s2d%bin_type
		nx=s2d%nx;	ny=s2d%ny;
		write(unit=file_unit) s2d%xcenter(1:nx), s2d%ycenter(1:ny), s2d%fxy(1:nx,1:ny)
		write(unit=file_unit) s2d%xstep, s2d%ystep
		
	end subroutine

#ifdef HDF5
	subroutine save_s2d_hdf5(s2d,group_id, s2dname)
		use hdf5
		use h5lt	
		implicit none
		character*(*) s2dname
		INTEGER(HID_T) :: group_id,sub_group_id      ! Group identifier
		INTEGER(HSIZE_T), DIMENSION(2) :: dims2  ! Dataset dimensions
		INTEGER(HSIZE_T), DIMENSION(1) :: dimsx,dimsy  ! Dataset dimensions
		integer error
		class(s2d_type)::s2d
	
		call h5gcreate_f(group_id, trim(adjustl(s2dname)), sub_group_id, error)

		dims2=(/s2d%nx,s2d%ny/)
		dimsx=s2d%nx;dimsy=s2d%ny
		!print*, "1"
		!print*, "allocated?", allocated(s2d%xcenter), allocated(s2d%ycenter), allocated(s2d%fxy)
		call h5ltmake_dataset_double_f(sub_group_id, "X",  1, dimsx,   s2d%xcenter, error)
		call h5ltmake_dataset_double_f(sub_group_id, "Y",  1, dimsy,   s2d%ycenter, error)
		!print*, "2"
		call h5ltmake_dataset_double_f(sub_group_id, "FXY", 2, dims2,  s2d%fxy(:,:), error)
		call h5gclose_f(sub_group_id, error) 
	end subroutine
	
	subroutine read_s2d_hdf5(s2d,group_id, s2dname)
		use hdf5
		use h5lt	
		implicit none
		character*(*) s2dname
		INTEGER(HID_T) :: group_id,sub_group_id      ! Group identifier
		INTEGER(HSIZE_T), DIMENSION(2) :: dims2  ! Dataset dimensions
		INTEGER(HSIZE_T), DIMENSION(1) :: dimsx,dimsy  ! Dataset dimensions
		integer error
		class(s2d_type)::s2d	

		call h5gcreate_f(group_id, trim(adjustl(s2dname)), sub_group_id, error)

		dims2=(/s2d%nx,s2d%ny/)
		dimsx=s2d%nx;dimsy=s2d%ny
		
		call h5ltread_dataset_f(sub_group_id, "X",H5T_NATIVE_DOUBLE,  s2d%xcenter, dimsx,  error)
		call H5LTread_dataset_f(sub_group_id, "Y",H5T_NATIVE_DOUBLE,   s2d%ycenter, dimsy,  error)
		call H5LTread_dataset_f(sub_group_id, "FXY",H5T_NATIVE_DOUBLE,  s2d%fxy(:,:), dims2, error)
		call h5gclose_f(sub_group_id, error) 
	end subroutine
	subroutine output_s2d_hdf5(s2d,s2dname, fl)
		use md_hdf5
		implicit none
		type(hdf5_file_type)::hf
		class(s2d_type)::s2d
		character*(*) s2dname, fl
		call hf%open(trim(adjustl(fl))//".hdf5")
		call s2d%save_hdf5(hf%file_id,s2dname)
		call hf%close()
	end subroutine
#endif

!=====================================================================================================
	subroutine init_s2d_hst(s2d,nx,ny,xmin,xmax,ymin,ymax,use_weight)
		implicit none
		class(s2d_hst_type)::s2d
		integer nx,ny
		integer i
		real(8) xmin,xmax,ymin,ymax
		logical,optional::use_weight
		
		call init_s2d_hst_basic(s2d,nx,ny,xmin,xmax,ymin,ymax,use_weight)
	end subroutine
	subroutine read_s2d_hst(s2d,file_unit, iostat, iomsg)
		implicit none
		class(s2d_hst_type),intent(inout)::s2d
		integer nx,ny,i
		logical proj, use_weight
		integer,intent(in):: file_unit
		integer, intent(out) :: iostat
		character(*), intent(inout) :: iomsg

		read(unit=file_unit) s2d%nx,s2d%ny,s2d%xmin,s2d%xmax,s2d%ymin,s2d%ymax, use_weight
		nx=s2d%nx;	ny=s2d%ny;
		!print*, "nx, ny=", nx, ny
		call s2d%init(nx,ny,s2d%xmin,s2d%xmax,s2d%ymin,s2d%ymax, use_weight)
		read(unit=file_unit) s2d%xcenter(1:nx), s2d%ycenter(1:ny), s2d%fxy(1:nx,1:ny)
		read(unit=file_unit) s2d%nxy(1:nx,1:ny), s2d%xmean(1:ny),s2d%ymean(1:nx),s2d%xsct(1:ny),s2d%ysct(1:nx)
		read(unit=file_unit) s2d%xstep, s2d%ystep
	end subroutine

	subroutine write_s2d_hst(s2d,file_unit, iostat, iomsg)
		implicit none
		class(s2d_hst_type),intent(in)::s2d
		integer nx,ny,i
		integer,intent(in):: file_unit
		logical proj
		integer, intent(out) :: iostat
		character(*), intent(inout) :: iomsg

		write(unit=file_unit) s2d%nx,s2d%ny,s2d%xmin,s2d%xmax,s2d%ymin,s2d%ymax, s2d%use_weight
		nx=s2d%nx;	ny=s2d%ny;
		write(unit=file_unit) s2d%xcenter(1:nx), s2d%ycenter(1:ny), s2d%fxy(1:nx,1:ny)
		write(unit=file_unit) s2d%nxy(1:nx,1:ny), s2d%xmean(1:ny),s2d%ymean(1:nx),s2d%xsct(1:ny),s2d%ysct(1:nx)
		!proj=allocated(s2d%da_x_proj(1)%x)
		write(unit=file_unit) s2d%xstep, s2d%ystep

	end subroutine
#ifdef HDF5
    subroutine save_s2d_hst_hdf5(s2d, group_id, s2dname)
        use md_hdf5
        implicit none
        class(s2d_hst_type)::s2d
        integer(HID_T):: sub_group_id, group_id
        character*(*) s2dname
        integer error
        INTEGER(HSIZE_T), DIMENSION(2) :: dims2  ! Dataset dimensions
        INTEGER(HSIZE_T), DIMENSION(1) :: dimsx,dimsy  ! Dataset dimensions


        call h5gcreate_f(group_id, trim(adjustl(s2dname)), sub_group_id, error)
        
    
        dims2=(/s2d%nx,s2d%ny/)
        dimsx=s2d%nx;dimsy=s2d%ny
        
        call h5ltmake_dataset_double_f(sub_group_id, "X",  1, dimsx,   s2d%xcenter, error)
        call h5ltmake_dataset_double_f(sub_group_id, "Y",  1, dimsy,   s2d%ycenter, error)
        call h5ltmake_dataset_double_f(sub_group_id, "FXY", 2, dims2,  s2d%fxy(:,:), error)
        call h5ltmake_dataset_double_f(sub_group_id, "FXYW", 2, dims2, s2d%fxyw(:,:), error)
    !  	call h5ltmake_dataset_double_f(group_id, "NXY", 2, dims2, s2d%nxy(:,:), error)
        call h5ltmake_dataset_int_f(sub_group_id, "NXY", 2, dims2,  s2d%nxy(:,:), error)
        call h5gclose_f(sub_group_id, error)   
    end subroutine
	subroutine read_s2d_hst_hdf5(s2d, group_id, s2dname)
        use md_hdf5
		use h5lt		
        implicit none
        class(s2d_hst_type)::s2d
        integer(HID_T):: sub_group_id, group_id
        character*(*) s2dname
        integer error
        INTEGER(HSIZE_T), DIMENSION(2) :: dims2  ! Dataset dimensions
        INTEGER(HSIZE_T), DIMENSION(1) :: dimsx,dimsy  ! Dataset dimensions


        call h5gopen_f(group_id, trim(adjustl(s2dname)), sub_group_id, error)
            
        dims2=(/s2d%nx,s2d%ny/)
        dimsx=s2d%nx;dimsy=s2d%ny
        
        call h5ltread_dataset_f(sub_group_id, "X",H5T_NATIVE_DOUBLE,  s2d%xcenter, dimsx,  error)
        call H5LTread_dataset_f(sub_group_id, "Y",H5T_NATIVE_DOUBLE,   s2d%ycenter, dimsy,  error)
        call H5LTread_dataset_f(sub_group_id, "FXY",H5T_NATIVE_DOUBLE,  s2d%fxy(:,:), dims2, error)
        call H5LTread_dataset_f(sub_group_id, "FXYW",H5T_NATIVE_DOUBLE, s2d%fxyw(:,:), dims2, error)
    !  	call H5LTread_dataset_double_double_f(group_id, "NXY", 2, dims2, s2d%nxy(:,:), error)
        call H5LTread_dataset_f(sub_group_id, "NXY",H5T_NATIVE_INTEGER,  s2d%nxy(:,:), dims2, error)

        call h5gclose_f(sub_group_id, error)   
    end subroutine
#endif
end module
