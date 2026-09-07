module sidm_params
  use unitless_value, only: PI
  implicit none
  real(8), save :: sigma_over_m = 13.6d0
  real(8), save :: w_damp = 3.0d0
  real(8), save :: n_power = 4.0d0
  real(8), save :: v0_code = 1.0d0
  real(8), save :: m_particle = 1.0d0
  logical, save :: use_sidm = .true.

  real(8), parameter :: AU_CM = 1.495979d13
  real(8), parameter :: MSUN_G = 1.98892d33
  real(8), parameter :: VUNIT_KMS = 29.7846d0
  real(8), parameter :: CM2G_TO_CODE = MSUN_G/AU_CM**2

  real(8), save :: n0_phys_pc3 = 2860.0d0
  logical, save :: use_gravothermal_n0 = .true.
  real(8), save :: spike_gamma = 2.3333333d0
  real(8), save :: spike_mass_factor = 2.0d0
  logical, save :: use_spike_n0 = .false.

  ! sigma0 is the zero-velocity total cross section in the t-channel Born
  ! kernel.  The run-local sidm_kernel.in file contains sigma0/m [cm^2/g]
  ! and w [km/s].
  real(8), save :: sigma0_cgs = 30.0d0
  real(8), save :: omega_kms = 80.0d0

contains
  subroutine initialize_sidm_kernel(v0_in)
    implicit none
    real(8), intent(in) :: v0_in
    integer :: iu, ios
    logical :: have_kernel

    if (v0_in <= 0.0d0 .or. v0_in /= v0_in) &
      error stop "SIDM_KERNEL: invalid GNC velocity unit"
    v0_code = v0_in
    inquire(file="sidm_kernel.in", exist=have_kernel)
    if (have_kernel) then
      open(newunit=iu,file="sidm_kernel.in",status="old",action="read",iostat=ios)
      if (ios /= 0) error stop "SIDM_KERNEL: cannot open sidm_kernel.in"
      read(iu,*,iostat=ios) sigma0_cgs, omega_kms
      close(iu)
      if (ios /= 0 .or. sigma0_cgs <= 0.0d0 .or. omega_kms <= 0.0d0) &
        error stop "SIDM_KERNEL: malformed physical kernel parameters"
      sigma_over_m = sigma0_cgs * CM2G_TO_CODE
      w_damp = omega_kms / VUNIT_KMS
    end if
    if (w_damp <= 0.0d0 .or. sigma_over_m <= 0.0d0) &
      error stop "SIDM_KERNEL: non-positive code parameters"
    print *, "SIDM_KERNEL_INITIALIZED", sigma0_cgs, omega_kms, &
             sigma_over_m, w_damp, v0_code
  end subroutine initialize_sidm_kernel

  real(8) function log_lambda_eff(energy_dimless)
    implicit none
    real(8), intent(in) :: energy_dimless
    real(8) :: v_over_wdamp, abs_e, aa

    if (.not. use_sidm) then
      log_lambda_eff = 1.0d0
      return
    end if
    abs_e = abs(energy_dimless)
    v_over_wdamp = sqrt(2.0d0*abs_e)*v0_code/w_damp
    aa = v_over_wdamp**2
    if (v_over_wdamp <= 0.0d0) then
      log_lambda_eff = tiny(1.0d0)
    else if (aa < 1.0d-3) then
      ! Stable series for [(a+2)log(1+a)-2a]/a.
      log_lambda_eff = aa**2/6.0d0 - aa**3/6.0d0 + &
                       3.0d0*aa**4/20.0d0 - 2.0d0*aa**5/15.0d0
    else
      log_lambda_eff = ((aa+2.0d0)*log(1.0d0+aa)-2.0d0*aa)/aa
    end if
    if (log_lambda_eff /= log_lambda_eff .or. log_lambda_eff < 0.0d0) &
      error stop "SIDM_KERNEL: non-finite effective logarithm"
    log_lambda_eff = max(log_lambda_eff, tiny(1.0d0))
  end function log_lambda_eff

  real(8) function log_lambda_eff_avg()
    implicit none
    log_lambda_eff_avg = log_lambda_eff(1.0d0)
  end function log_lambda_eff_avg

  real(8) function sidm_kappa0()
    implicit none
    ! Combined with log_lambda_eff, this is exactly the viscosity moment of
    ! the t-channel Born kernel mapped onto GNC's Landau normalization.
    sidm_kappa0 = 2.0d0*PI*sigma_over_m*m_particle*w_damp**4
  end function sidm_kappa0

  real(8) function orbital_avg_validity(energy_dimless, v0, n0_local)
    implicit none
    real(8), intent(in) :: energy_dimless, v0, n0_local
    real(8) :: t_relax, period, lnL, kap
    lnL = log_lambda_eff(energy_dimless)
    kap = sidm_kappa0()
    t_relax = 0.34d0*v0**3/(kap*lnL*n0_local*m_particle**2)
    period = 1.0d0/(2.0d0*max(energy_dimless,1.0d-30))**1.5d0
    if (period > 0.0d0) then
      orbital_avg_validity = t_relax/period
    else
      orbital_avg_validity = 1.0d30
    end if
  end function orbital_avg_validity
end module sidm_params
