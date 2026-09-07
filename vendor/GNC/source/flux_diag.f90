module flux_diag
  use com_main_gw
  implicit none
  logical, save :: flux_enabled = .false.
contains

  subroutine output_boundary_fluxes(dm, isnap)
    type(diffuse_mspec), intent(in) :: dm
    integer, intent(in) :: isnap
    integer :: i, j, ib, nE, nJ, uu, nden
    real(8) :: f, dfdx, DE, DEE, dxE, dJ
    real(8) :: FM, FE, Fp, denb
    real(8), allocatable :: xE(:), yJ(:)
    real(8) :: lemin, lemax, ljmin, ljmax
    character(16) :: sid

    if (.not. flux_enabled) return
    nE = dm%nbin_grid
    nJ = dm%nbin_gx
    if (nE < 2 .or. nJ < 2) return

    allocate(xE(nE), yJ(nJ))
    lemin = log(dm%emin);  lemax = log(dm%emax)
    do i = 1, nE
       xE(i) = exp(lemin + (lemax-lemin)*dble(i-1)/dble(nE-1))
    end do
    ljmin = log(max(dm%jmin,1d-30));  ljmax = log(dm%jmax)
    do j = 1, nJ
       yJ(j) = exp(ljmin + (ljmax-ljmin)*dble(j-1)/dble(nJ-1))
    end do
    nden = size(dm%mb(1)%all%fden%fx)

    write(sid,'(I0)') isnap
    open(newunit=uu, file='output/pro/fluxes_'//trim(adjustl(sid))//'.txt', &
         status='replace', action='write')
    write(uu,'(A)')     '# FP boundary fluxes (mb(1))'
    write(uu,'(A,2I6)') '# nE nJ = ', nE, nJ
    write(uu,'(A)') '# ib   E_b            F_M            F_E(L)         F_p            rho_b'
    do ib = 1, nE
       FM = 0d0; FE = 0d0; Fp = 0d0
       do i = 1, ib
          if (i < nE) then
             dxE = xE(i+1) - xE(i)
          else
             dxE = xE(i) - xE(i-1)
          end if
          do j = 1, nJ
             if (j < nJ) then
                dJ = yJ(j+1) - yJ(j)
             else
                dJ = yJ(j) - yJ(j-1)
             end if
             f = dm%mb(1)%all%gxj%fxy(i,j)
             if (i > 1 .and. i < nE) then
                dfdx = ( dm%mb(1)%all%gxj%fxy(i+1,j) &
                       - dm%mb(1)%all%gxj%fxy(i-1,j) ) / ( xE(i+1) - xE(i-1) )
             else
                dfdx = 0d0
             end if
             DE  = dm%mb(1)%dc%s2_de_110%fxy(i,j) + dm%mb(1)%dc%s2_de_0%fxy(i,j)
             DEE = dm%mb(1)%dc%s2_dee%fxy(i,j)
             FE = FE + ( DE*f - DEE*dfdx ) * xE(i) * dxE * dJ
             FM = FM + ( DE*f - DEE*dfdx )          * dxE * dJ
          end do
       end do
       denb = 0d0
       if (nden > 0) denb = dm%mb(1)%all%fden%fx(min(ib,nden))
       Fp = denb * (2d0*xE(ib))
       write(uu,'(I6,5ES15.6)') ib, xE(ib), FM, FE, Fp, denb
    end do
    close(uu)
    deallocate(xE, yJ)
  end subroutine output_boundary_fluxes

end module flux_diag
