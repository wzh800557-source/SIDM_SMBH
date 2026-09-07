module sidm_params
    use unitless_value, only: PI
    implicit none
    real(8), save :: sigma_over_m  = 13.6d0   ! FP-valid validation value
    real(8), save :: w_damp        = 3.0d0
    real(8), save :: n_power       = 4.0d0
    real(8), save :: m_particle    = 1.0d0
    logical, save :: use_sidm      = .true.
    ! ---- physical<->code unit conversion (code units: length=AU, mass=Msun, G=1) ----
    real(8), parameter :: AU_CM       = 1.495979d13     ! cm per AU
    real(8), parameter :: MSUN_G      = 1.98892d33      ! g  per Msun
    real(8), parameter :: VUNIT_KMS   = 29.7846d0       ! sqrt(G Msun/AU) [km/s]
    real(8), parameter :: CM2G_TO_CODE= MSUN_G/AU_CM**2 ! (cm^2/g) -> (AU^2/Msun)
    ! gravothermal-derived density at handoff radius [Msun/pc^3, m=1 Msun]
    real(8), save :: n0_phys_pc3 = 2860.0d0  ! halo-continuity spike-matched n(r_h)
    logical, save :: use_gravothermal_n0 = .true.
    ! Gondolo-Silk DM spike: M_DM(<r_h)=spike_mass_factor*M_BH, rho~r^-spike_gamma
    real(8), save :: spike_gamma       = 2.3333333d0  ! (9-2g)/(4-g): 7/3 for initial NFW g=1
    real(8), save :: spike_mass_factor = 2.0d0
    logical, save :: use_spike_n0      = .false. ! use matched n0 via gravothermal path

    ! --- PHYSICAL CALIBRATION (Jiang+25 LRD values) ---
    ! call set_sidm_physical(30.0d0, 80.0d0) => sigma0/m=30 cm^2/g, w=80 km/s
    !   => sigma_over_m=2.666d8, w_damp=2.686d0 (code units)
    ! WARNING: at this physical cross-section t_relax ~ 1e4 yr << orbital period
    !   => Fokker-Planck (orbit-averaging) BREAKS DOWN near the BH; the inner
    !      cusp is in the collisional FLUID regime. Do NOT run GNC here; use the
    !      gravothermal fluid code for r < r_h. (The code default 13.6 = 1.5e-6
    !      cm^2/g is effectively collisionless.)
    real(8), save :: sigma0_cgs = 30.0d0   ! reference physical value
    real(8), save :: omega_kms  = 80.0d0   ! reference physical value
    real(8), parameter :: LNL_MIN = 0.1d0
    real(8), parameter :: LNL_MAX = 30.0d0
contains
    real(8) function log_lambda_eff(energy_dimless)
        implicit none
        real(8), intent(in) :: energy_dimless
        real(8) :: v_over_wdamp, abs_e
        if (.not. use_sidm) then
            log_lambda_eff = 1.0d0; return
        end if
        abs_e = abs(energy_dimless)
        if (abs_e < 1.0d-30) abs_e = 1.0d-30
        v_over_wdamp = sqrt(2.0d0 * abs_e) / w_damp
        if (v_over_wdamp < 1.0d-10) then
            log_lambda_eff = 0.5d0 * n_power * log(1.0d0 / v_over_wdamp)
        else
            log_lambda_eff = 0.5d0 * log(1.0d0 + (1.0d0 / v_over_wdamp)**n_power)
        end if
        if (log_lambda_eff /= log_lambda_eff) log_lambda_eff = LNL_MIN
        if (log_lambda_eff < LNL_MIN) log_lambda_eff = LNL_MIN
        if (log_lambda_eff > LNL_MAX) log_lambda_eff = LNL_MAX
    end function
    real(8) function log_lambda_eff_avg()
        implicit none
        log_lambda_eff_avg = log_lambda_eff(1.0d0)
    end function
    real(8) function sidm_kappa0()
        implicit none
        sidm_kappa0 = 2.0d0 * PI * sigma_over_m * m_particle * w_damp**4
    end function
    real(8) function orbital_avg_validity(energy_dimless, v0, n0_local)
        implicit none
        real(8), intent(in) :: energy_dimless, v0, n0_local
        real(8) :: t_relax, period, lnL, kap
        lnL = log_lambda_eff(energy_dimless)
        kap = sidm_kappa0()
        t_relax = 0.34d0 * v0**3 / (kap * lnL * n0_local * m_particle**2)
        period  = 1.0d0 / (2.0d0 * max(energy_dimless, 1.0d-30))**1.5d0
        if (period > 0.0d0) then
            orbital_avg_validity = t_relax / period
        else
            orbital_avg_validity = 1.0d30
        end if
    end function
end module sidm_params
