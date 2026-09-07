module md_cfuns
    type ::cfuns_type
        integer:: nj, ns
        real(8) jmin, jmax
        real(8),allocatable:: cfs_110(:,:)
        real(8),allocatable:: cfs_111(:,:)
        real(8),allocatable:: cfs_131(:,:)
        real(8),allocatable:: cfs_13_1(:,:)
        real(8),allocatable:: cfs_130(:,:)
        real(8),allocatable:: cfs_330(:,:)
        real(8),allocatable:: cfs_310(:,:)
        !real(8),allocatable:: cfs_010(:,:)
        !real(8),allocatable:: cfs_0_11(:,:)
        real(8),allocatable:: jum(:)
        real(8),allocatable:: s(:,:)
    contains
        procedure::get_size=>get_size_of_cfs
        procedure::input_bin=>input_bin_cfs
        procedure::output_bin=>output_bin_cfs
        procedure::init=>init_cfs
        procedure::get_cfs_s
        procedure::get_cfs_s_mpi
    end type
    private::input_bin_cfs, output_bin_cfs,get_size_of_cfs
    type(cfuns_type)::cfs
    real(8),parameter::dsmin_value=-8
contains
    subroutine init_cfs(this, nj, ns, jmin, jmax)
        implicit none
        class(cfuns_type)::this
        integer nj, ns
        real(8) jmin, jmax
        this%nj=nj; this%ns=ns; this%jmin=jmin; this%jmax=jmax
        if(allocated(this%cfs_110))then
            deallocate(this%cfs_110,this%cfs_111,this%cfs_131,this%cfs_13_1, &
            this%cfs_130,this%cfs_330,this%cfs_310, this%jum,this%s)
        !    this%cfs_010,this%cfs_0_11)
        end if
        allocate(this%cfs_110(nj,ns),this%cfs_111(nj,ns),this%cfs_131(nj,ns),&
        this%cfs_13_1(nj,ns),this%cfs_130(nj,ns),this%cfs_330(nj,ns),&
        this%cfs_310(nj,ns), this%s(nj, ns), this%jum(nj))
        !this%cfs_010(nj,ns),this%cfs_0_11(ns,ns))
    end subroutine
    subroutine get_size_of_cfs(this,n)
        implicit none
        class(cfuns_type)::this
        integer n
        n=sizeof(this%cfs_110)/1024*8
        n=n+sizeof(this%jum)/1024
    end subroutine
    subroutine input_bin_cfs(this, fl)
        implicit none
        integer nj,ns
        class(cfuns_type)::this
        real(8) jmin, jmax
        character*(*) fl
        integer,parameter::funit=19923923
        open(unit=funit, file=trim(adjustl(fl))//".bin", &
            status='old',access='stream',form='unformatted')
        read(funit) nj, ns,jmin, jmax
        call this%init(nj, ns,jmin, jmax)
        read(funit) this%cfs_110,this%cfs_111,this%cfs_130,this%cfs_13_1,&
            this%cfs_131,this%cfs_310,this%cfs_330, this%s, this%jum
            
        close(unit=funit)
    end subroutine
    subroutine output_bin_cfs(this, fl)
        implicit none
        integer nj, ns
        class(cfuns_type)::this
        character*(*) fl
        integer,parameter::funit=19923923
        open(unit=funit, file=trim(adjustl(fl))//".bin", &
            access='stream',form='unformatted')
        write(funit) this%nj, this%ns,this%jmin, this%jmax
        write(funit) this%cfs_110,this%cfs_111,this%cfs_130,this%cfs_13_1,&
            this%cfs_131,this%cfs_310,this%cfs_330, this%s, this%jum
           
        close(unit=funit)
    end subroutine
    subroutine get_cfs_s(this)
        implicit none
        class(cfuns_type)::this
        call set_grid_cfs(this%jum, this%s, this%nj, this%ns, this%jmin, this%jmax)
        call get_cfunction_grid_s(1,1,1, this%cfs_111, this%nj, this%ns, this%jmin, this%jmax)
        call get_cfunction_grid_s(1,1,0, this%cfs_110, this%nj, this%ns, this%jmin, this%jmax)
        call get_cfunction_grid_s(1,3,0, this%cfs_130, this%nj, this%ns, this%jmin, this%jmax)
        call get_cfunction_grid_s(1,3,-1, this%cfs_13_1, this%nj, this%ns, this%jmin, this%jmax)
        call get_cfunction_grid_s(1,3,1, this%cfs_131, this%nj, this%ns, this%jmin, this%jmax)
        call get_cfunction_grid_s(3,1,0, this%cfs_310, this%nj, this%ns, this%jmin, this%jmax)
        !print*, "start 330"
        call get_cfunction_grid_s(3,3,0, this%cfs_330, this%nj, this%ns, this%jmin, this%jmax)
        !print*, "cfs 0 1 0"
        !call get_cfunction_grid_s(0,1,0, this%cfs_010, this%nj, this%ns, this%jmin, this%jmax)
        !print*, "cfs 0 -1 1"
        !call get_cfunction_grid_s(0,-1,1, this%cfs_0_11, this%nj, this%ns, this%jmin, this%jmax)
        !print*, this%cfs_330
        !read(*,*)
    end subroutine
    subroutine get_cfs_s_mpi(this)
        implicit none
        class(cfuns_type)::this
        call set_grid_cfs(this%jum, this%s, this%nj, this%ns, this%jmin, this%jmax)
        print*, "start 111"
        call get_cfunction_grid_s_mpi(1,1,1, this%cfs_111, this%nj, this%ns, this%jmin, this%jmax)
        print*, "start 110"
        call get_cfunction_grid_s_mpi(1,1,0, this%cfs_110, this%nj, this%ns, this%jmin, this%jmax)
        print*, "start 130"
        call get_cfunction_grid_s_mpi(1,3,0, this%cfs_130, this%nj, this%ns, this%jmin, this%jmax)
        print*, "start 13-1"
        call get_cfunction_grid_s_mpi(1,3,-1, this%cfs_13_1, this%nj, this%ns, this%jmin, this%jmax)
        print*, "start 131"
        call get_cfunction_grid_s_mpi(1,3,1, this%cfs_131, this%nj, this%ns, this%jmin, this%jmax)
        print*, "start 310"
        call get_cfunction_grid_s_mpi(3,1,0, this%cfs_310, this%nj, this%ns, this%jmin, this%jmax)
        print*, "start 330"
        call get_cfunction_grid_s_mpi(3,3,0, this%cfs_330, this%nj, this%ns, this%jmin, this%jmax)
        !print*, "cfs 0 1 0"
        !call get_cfunction_grid_s(0,1,0, this%cfs_010, this%nj, this%ns, this%jmin, this%jmax)
        !print*, "cfs 0 -1 1"
        !call get_cfunction_grid_s(0,-1,1, this%cfs_0_11, this%nj, this%ns, this%jmin, this%jmax)
        !print*, this%cfs_330
        !read(*,*)
    end subroutine

end module
subroutine set_grid_cfs(jum,s, nj,ns, jmin,jmax)
    implicit none
    integer nj, ns
    real(8) jum(nj),s(nj, ns)
    real(8) jmin, jmax,smax,smin,e,tmin,tmax,ds,dsmin,dsmax
    integer i,j
    dsmin=-7d0; 
    smin=log10(10**dsmin+1d0)
    do i=1, nj
        jum(i)=(jmax-jmin)*real(i-1)/real(nj-1)+jmin
        e=sqrt(1-(10d0**jum(i))**2)
        dsmax=log10(2/(1-e)-1d0)
        smax=log10(10**dsmax+1d0)
        do j=1, ns
            ds=(dsmax-dsmin)*real(j-1)/real(ns-1)+dsmin
            s(i,j)=log10(10**ds+1d0)
        enddo
    enddo
end subroutine
subroutine get_cfunction_grid_s(l,m,n,cfs,xbin, ybin, jmin, jmax)
    implicit none
    integer l,m,n
    integer xbin, ybin
    real(8) cfs(xbin, ybin)
    real(8) s, e, t, smin,smax, jmin,jmax, jum
    real(8) dsmin,dsmax,ds
    integer i,j
    dsmin=-7d0; 
    smin=log10(10**dsmin+1d0)
    do i=1, xbin
        ! e value
        jum=(jmax-jmin)*real(i-1)/real(xbin-1)+jmin
        e=sqrt(1-(10d0**jum)**2)
        dsmax=log10(2/(1-e)-1d0)
        smax=log10(10**dsmax+1d0)
        !print*, "e,smax=",e,smax
        do j=1, ybin
            ! s value
            ds=(dsmax-dsmin)*real(j-1)/real(ybin-1)+dsmin
            s=log10(10**ds+1d0)
            !s=1d0;e=emax
            !print*, "l,m,n=",l,m,n
            call get_cfunction(10**s,e,cfs(i,j),l,m,n)
            !print*, "s,cfs=", s, cfs(i,j)
            !read(*,*)
        end do
    end do
end subroutine

subroutine get_cfunction_grid_s_mpi(l,m,n,cfs,xbin, ybin, jmin, jmax)
    use MPI_comu
    implicit none
    integer l,m,n
    integer xbin, ybin
    real(8) cfs(xbin, ybin)
    real(8) s, e, t, smin,smax, jmin,jmax, jum
    real(8) dsmin,dsmax,ds
    integer i,j,ierr
	integer nbg, ned, nblock, ntasks

    dsmin=-7d0; 
    smin=log10(10**dsmin+1d0)

    ctl%nblock_size=int(xbin/ctl%ntasks)
    ctl%nblock_mpi_bg=ctl%nblock_size*rid+1
    ctl%nblock_mpi_ed=ctl%nblock_size*(rid+1)
    
	nbg=ctl%nblock_mpi_bg
	ned=ctl%nblock_mpi_ed
	nblock=ctl%nblock_size
	ntasks=ctl%ntasks

    do i=nbg, ned
        ! e value
        jum=(jmax-jmin)*real(i-1)/real(xbin-1)+jmin
        e=sqrt(1-(10d0**jum)**2)
        dsmax=log10(2/(1-e)-1d0)
        smax=log10(10**dsmax+1d0)
        !print*, "e,smax=",e,smax
        do j=1, ybin
            ! s value
            ds=(dsmax-dsmin)*real(j-1)/real(ybin-1)+dsmin
            s=log10(10**ds+1d0)
            !s=1d0;e=emax
            !print*, "l,m,n=",l,m,n
            call get_cfunction(10**s,e,cfs(i,j),l,m,n)
            !print*, "s,cfs=", s, cfs(i,j)
            !read(*,*)
        end do
    end do
    call mpi_barrier(mpi_comm_world,ierr)
    call collect_data_mpi(cfs, xbin,nbg, ned, nblock, ntasks)
end subroutine

subroutine get_cfunction(s,e,cfs,ld,md,nd)
    use constant
    use my_intgl
    implicit none
    integer ld,md,nd,idid
    real(8) s,e, cfs
    real(8) ymin,ymax, y_out
    ymin=-pi/2d0
    ymax=asin(max(min(1d0, (2/s-1)/e ),-1d0))
    y_out=0d0

    !print*, "s,e=",s,e
    !print*, "ymin,ymax=",ymin,ymax
    if(ymin.ne.ymax)then
        call my_integral_none(ymin, ymax, y_out, fcn, idid)
    else
        y_out=0
    endif
    cfs=2d0**(1d0-ld)*y_out
contains
    subroutine fcn(n, x, y, f, par, ipar)
        use, intrinsic :: ieee_arithmetic
        implicit none
        integer n, ipar(100)
        real(8) x, y(n), f(n), par(100)

        f(1)=(1+e*sin(x))**(ld+nd)/(1-e*sin(x))**(nd+md/2d0) &
        *(abs(2-s-s*e*sin(x)))**(md/2d0)
        !if(2-s-s*e*sin(x)<0d0)then
        !    print*, "error!", 2-s-s*e*sin(x)
        !end if
        !print*, "x, f(1)=",x, f(1), (1+e*sin(x))**(ld+nd), (1-e*sin(x))**(nd+md/2d0), &
        !(2-s-s*e*sin(x))**(md/2d0)
        
    end subroutine
end subroutine