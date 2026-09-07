! ============================================================================
!  load_df_from_table.f90  --  GNC two-way coupling: re-seed the initial DF
!                              from an evolved fluid profile.
!
!  Integration:
!    * this file is in ~/GNC/source and listed in the makefile OBJS.
!    * in ~/GNC/main/ini.f90 , after `call set_dm_init(dms)`:
!          call load_df_from_table(dms, "f_evolved.txt")
!    * f_evolved.txt is produced by emit_gnc_df.py on GNC's OWN (x,j) grid
!      (same emin/emax, jmin/jmax, log spacing, 24x24), so we copy directly
!      into gxj%fxy with a size check -- no axis-array reads needed.
!
!  DF storage (confirmed in dms.f90 + compiler): DF = dm%mb(i)%star%gxj%fxy(:,:)
!    of type s2d_basic_type with size members %nx (energy), %ny (j).
!    Aggregate dm%all rebuilt by normalize_gxj (dms.f90).
!
!  NOTE: direct copy assumes emit_gnc_df.py's x,j grid matches GNC's gxj grid
!  ordering. Verify orientation once by comparing GNC's recovered density to the
!  input profile (dms_0_0.hdf5 /star/fden); if flipped, reverse a row/col index.
! ============================================================================
subroutine load_df_from_table(dm, fname)
    use com_main_gw
    implicit none
    type(diffuse_mspec)      :: dm
    character(*), intent(in) :: fname
    integer :: nx, nj, i, a, b, u
    real(8), allocatable :: xin(:), jin(:), fin(:,:)

    open(newunit=u, file=trim(fname), status='old', action='read')
    read(u,*)                       ! comment header
    read(u,*) nx, nj
    allocate(xin(nx), jin(nj), fin(nx,nj))
    read(u,*) xin
    read(u,*) jin
    do a = 1, nx
        read(u,*) fin(a,:)
    end do
    close(u)

    do i = 1, dm%n
        if (dm%mb(i)%star%gxj%nx /= nx .or. dm%mb(i)%star%gxj%ny /= nj) then
            print*, "load_df_from_table: GRID MISMATCH gnc=", &
                    dm%mb(i)%star%gxj%nx, dm%mb(i)%star%gxj%ny, " table=", nx, nj
            stop
        end if
        do a = 1, nx
            do b = 1, nj
                dm%mb(i)%star%gxj%fxy(a,b) = max(fin(a,b), 0.d0)   ! clip inner f<0
            end do
        end do
    end do

    call normalize_gxj(dm)          ! re-aggregate dm%all from dm%mb
    print*, "load_df_from_table: injected ", trim(fname), " into gxj%fxy (", nx, "x", nj, ")"
    deallocate(xin, jin, fin)
end subroutine load_df_from_table
