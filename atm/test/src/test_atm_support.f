module test_atm_support

  use const_def
  use math_lib
  use atm_def
  use atm_lib
  use chem_lib, only: composition_info
  use chem_def
  use eos_def
  use eos_lib
  use kap_def
  use kap_lib

  implicit none

  logical, parameter :: SKIP_PARTIALS = .TRUE.

  logical :: test_verbosely
  integer :: eos_handle, kap_handle
  real(dp) :: cgrav, Pextra_factor
  ! composition info
  integer, parameter :: species = 7
  integer, pointer :: chem_id(:), net_iso(:)
  real(dp) :: X, Y, Z, XC, XN, XO, xa(species), abar, zbar, z53bar

  integer :: ierr, max_iters, max_steps, iters
  real(dp) :: logg, Teff, M, R, L, kap_guess, tau, tau_base, tau_phot, &
       lnT, dlnT_dL, dlnT_dlnR, dlnT_dlnM, & 
       lnP, dlnP_dL, dlnP_dlnR, dlnP_dlnM, &
       atol, rtol, kap, T, err, P, Prad, &
       dlnT_dlnkap, dlnP_dlnkap, &
       logTeff, log_gsurf, log_Rsurf, log_M, g, &
       T_eq, kap_v, gamma, T_int


contains

  subroutine do_test_atm( &
       test_verbosely_in, cgrav_in, eos_handle_in, kap_handle_in)
    logical, intent(in) :: test_verbosely_in
    real(dp), intent(in) :: cgrav_in
    integer, intent(in) :: eos_handle_in, kap_handle_in

    test_verbosely = test_verbosely_in
    cgrav = cgrav_in
    eos_handle = eos_handle_in
    kap_handle = kap_handle_in

    call test_table(ATM_TABLE_TAU_1M1, 'tau=0.1')
    call test_table(ATM_TABLE_TAU_1, 'tau=1')
    call test_table(ATM_TABLE_TAU_10, 'tau=10')              
    call test_table(ATM_TABLE_TAU_100, 'tau=100')             

    call test_table(ATM_TABLE_WD_TAU_25, 'WD tau=25')
    call test_table(ATM_TABLE_DB_WD_TAU_25, 'DB WD tau=25')
    call test_table(ATM_TABLE_PHOTOSPHERE, 'photosphere')

    call test_T_tau_varying(ATM_T_TAU_EDDINGTON, 'Eddington', -1._dp)
    call test_T_tau_varying(ATM_T_TAU_KRISHNA_SWAMY, 'Krishna-Swamy', -1._dp)
    call test_T_tau_varying(ATM_T_TAU_SOLAR_HOPF, 'solar Hopf', -1._dp)
    call test_T_tau_varying(ATM_T_TAU_TRAMPEDACH_SOLAR, 'Trampedach solar', -1._dp)
    call test_T_tau_varying(ATM_T_TAU_EDDINGTON, 'Eddington', 100._dp)

    call test_T_tau_uniform('fixed', 100._dp)
    call test_T_tau_uniform('iterated', 150._dp)

    call test_irradiated()

  end subroutine do_test_atm

  !****

  subroutine test_table(table_id, label)

    integer, intent(in)      :: table_id
    character(*), intent(in) :: label

    include 'formats'

    ierr = 0

    if (test_verbosely) then
       write(*,*)
       write(*,*) 'test_table: ', label
    endif

    if (table_id == ATM_TABLE_WD_TAU_25) then
       logg = 7.1d0
       Teff = 5000
       M = 0.8*Msun
       R = sqrt(cgrav*M / exp10(logg))
       !write(*,*) 'R/Rsun', R/Rsun
       L = pi*crad*clight*R*R*Teff*Teff*Teff*Teff
       !write(*,*) 'L/Lsun', L/Lsun
    elseif (table_id == ATM_TABLE_DB_WD_TAU_25) then
       logg = 8.0d0
       Teff = 25000
       M = 0.8*Msun
       R = sqrt(cgrav*M / exp10(logg))
       !write(*,*) 'R/Rsun', R/Rsun
       L = pi*crad*clight*R*R*Teff*Teff*Teff*Teff
       !write(*,*) 'L/Lsun', L/Lsun    else
    else
       M = Msun
       R = Rsun
       L = Lsun
       Teff = atm_Teff(L, R)
    end if

    Z = 0.02d0

    call atm_eval_table( &
       L, R, M, cgrav, table_id, Z, SKIP_PARTIALS, &
       Teff, &
       lnT, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
       lnP, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
       ierr)
    if (ierr /= 0) then
       if (test_verbosely) write(*,*) 'failed in atm_eval_table'
       call mesa_error(__FILE__,__LINE__)
    end if

    T = exp(lnT)
    P = exp(lnP)

    Prad = radiation_pressure(T)

    if (test_verbosely) write(*,*)
    if (test_verbosely) write(*,'(99a16)') &
         'T', 'log_T', 'log P', 'M/Msun', 'L/Lsun', 'R/Rsun', 'logPgas'
    if (test_verbosely) write(*,'(i16,99f16.8)') &
         floor(0.5d0 + T), lnT/ln10, log10(P), M/Msun, L/Lsun, R/Rsun, log10(P-Prad)
    if (test_verbosely) write(*,*)

  end subroutine test_table

  !****

  subroutine test_T_tau_varying(T_tau_id, label, tau_base_in)

    integer, intent(in)      :: T_tau_id
    character(*), intent(in) :: label
    real(dp), intent(in)     :: tau_base_in

    real(dp) :: errtol
    integer  :: max_steps

    include 'formats'

    ierr = 0

    if (tau_base_in < 0._dp) then
       call atm_get_T_tau_base(T_tau_id, tau_base, ierr)
       if (ierr /= 0) then
          if (test_verbosely) write(*,*) 'failed in atm_eval_T_tau_varying'
          return
       end if
    else
       tau_base = tau_base_in
    endif

    if (test_verbosely) then
       write(*,*)
       write(*,*) 'test_T_tau_varying: ', label
    endif

    ! TEST SOLAR VALUES
    logTeff = log10(5776d0)
    log_gsurf = log10(cgrav*Msun/(Rsun*Rsun))
    log_Rsurf = log10(Rsun)
    log_M = log10(Msun)

    Z = 0.02d0
    X = 0.70d0
    XC = 3.2724592105263235d-03
    XN = 9.5023842105263292d-04
    XO = 8.8218000000000601d-03
    call set_composition()

    R = exp10(log_Rsurf)
    M = exp10(log_M)

    g = exp10(log_gsurf)
    Teff = exp10(logTeff)
    L = pi*crad*clight*R*R*Teff*Teff*Teff*Teff

    errtol = 1.E-9_dp
    max_steps = 500

    call atm_eval_T_tau_varying( &
       tau_base, L, R, M, cgrav, &
       T_tau_id, eos_proc, kap_proc, &
       errtol, max_steps, SKIP_PARTIALS, &
       Teff, &
       lnT, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
       lnP, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
       ierr)
    if (ierr /= 0) then
       if (test_verbosely) write(*,*) 'failed in atm_eval_T_tau_varying'
       call mesa_error(__FILE__,__LINE__)
    end if

    T = exp(lnT)
    P = exp(lnP)

    if (test_verbosely) write(*,*)
    if (test_verbosely) write(*,'(99a16)') 'T', 'log T', 'log P', 'M/Msun', 'L/Lsun', 'X', 'Z'
    if (test_verbosely) write(*,'(i16,99f16.8)') &
         floor(0.5d0 + T), log10(T), log10(P), M/Msun, L/Lsun, X, Z
    if (test_verbosely) write(*,*)

    deallocate(chem_id, net_iso)

  end subroutine test_T_tau_varying

  !****

  subroutine test_T_tau_uniform (T_tau_opacity, tau_base_in)

    character(*), intent(in) :: T_tau_opacity
    real(dp), intent(in)     :: tau_base_in

    real(dp) :: tau_base
    real(dp) :: errtol

    include 'formats'

    if (test_verbosely) then
       write(*,*)
       write(*,*) 'test_T_tau_uniform: ', TRIM(T_tau_opacity)
    endif

    tau_base = tau_base_in
    if (tau_base <= 0._dp) then
       tau_base = 2._dp/3._dp
    end if

    Pextra_factor = 1._dp

    select case (T_tau_opacity)

    case ('fixed')

       M = 1.9892000000000002D+32
       R = 6.3556231577545586D+10
       L = 2.4015399190199118D+32
       Teff = atm_Teff(L, R)

       kap_guess = 5.8850802481174469D-02

       max_iters = 0

    case ('iterated')

       logTeff = log10(5776._dp)
       log_gsurf = log10(cgrav*Msun/(Rsun*Rsun))
       log_Rsurf = log10(Rsun)
       log_M = log10(Msun)

       R = exp10(log_Rsurf)
       M = exp10(log_M)
       Teff = exp10(logTeff)
       L = pi*crad*clight*R*R*Teff*Teff*Teff*Teff

       kap_guess = 0.5d0

       max_iters = 30

    case default

       write(*,*) 'unknown value for T_tau_opacity ' // trim(T_tau_opacity)
       
       stop

    end select

    Z = 0.02d0
    X = 0.70d0
    XC = 3.2724592105263235d-03
    XN = 9.5023842105263292d-04
    XO = 8.8218000000000601d-03
    call set_composition()

    errtol = 1.E-6_dp

    call atm_eval_T_tau_uniform( &
       tau_base, L, R, M, cgrav, kap_guess, Pextra_factor, &
       ATM_T_TAU_EDDINGTON, eos_proc, kap_proc, errtol, max_iters, SKIP_PARTIALS, &
       Teff, kap, &
       lnT, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
       lnP, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
       ierr)    
    if (ierr /= 0) then
       if (test_verbosely) write(*,*) 'failed in atm_eval_T_tau_uniform'
       call mesa_error(__FILE__,__LINE__)
    end if
    if (test_verbosely) write(*,1) 'atm_get_grey logP surf', lnP/ln10
    if (test_verbosely) write(*,1) 'atm_get_grey logT surf', lnT/ln10
    if (test_verbosely) write(*,1) 'atm_get_grey tau_phot', tau_base
    if (test_verbosely) write(*,*)

    deallocate(chem_id, net_iso)

  end subroutine test_T_tau_uniform

  !****

  subroutine test_irradiated()

    real(dp) :: errtol
    type (Kap_General_Info), pointer :: rq

    include 'formats'

    if (test_verbosely) then
       write(*,*)
       write(*,*) 'test_irradiated'
    endif

    ierr = 0

    X = 0.70d0
    Z = 1d-2
    XC = 3.2724592105263235d-03
    XN = 9.5023842105263292d-04
    XO = 8.8218000000000601d-03

    call set_composition()

    ! at these conditions, the appropriate lowT opacity table is Freedman11
    ! this is around logR = 3.5, way off the default tables (max logR = 1)
    call kap_ptr(kap_handle,rq,ierr)
    if (ierr /= 0) return
    rq% kap_lowT_option = kap_lowT_Freedman11

    ! must set up tables again after changing options
    call kap_setup_tables(kap_handle, ierr)

    errtol = 1.E-6_dp
    max_iters = 30

    T_eq = 1000.
    kap_v = 4.E-3_dp
    kap_guess = 1.5d-2
    gamma = 0._dp
    P = 1d6
    M = 1.5d0*M_jupiter
    tau = 10 ! just a guess for use in getting R
    R = sqrt(cgrav*M*tau/(P*kap_guess)) ! g = P*kap/tau = G*M/R^2
    T_int = 900
    L = pi*crad*clight*R*R*T_int*T_int*T_int*T_int
    
    call atm_eval_irradiated( &
       L, R, M, cgrav, T_eq, P, kap_guess, kap_v, gamma, &
       eos_proc, kap_proc, errtol, max_iters, SKIP_PARTIALS, &
       Teff, kap, tau, &
       lnT, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
       ierr)
    if (ierr /= 0) then
       if (test_verbosely) write(*,*) 'bad return from atm_eval_irradiated'
       write(*,*) 'ierr',ierr
       call mesa_error(__FILE__,__LINE__)
    end if


    if (test_verbosely) write(*,*) 
    if (test_verbosely) write(*,1) 'test_grey_irradiated: kap', kap
    if (test_verbosely) write(*,1) 'M/M_jupiter', M/M_jupiter
    if (test_verbosely) write(*,1) 'R/R_jupiter', R/R_jupiter
    if (test_verbosely) write(*,1) 'R/Rsun', R/Rsun
    if (test_verbosely) write(*,1) 'kap_v', kap_v
    if (test_verbosely) write(*,1) 'P', P
    if (test_verbosely) write(*,*) 

    if (test_verbosely) write(*,1) 'tau', tau
    if (test_verbosely) write(*,1) 'Teff', Teff
    if (test_verbosely) write(*,1) 'T_eq', T_eq
    if (test_verbosely) write(*,1) 'T_int', T_int
    if (test_verbosely) write(*,1) 'T', exp(lnT)
    if (test_verbosely) write(*,1) 'logT', lnT/ln10
    if (test_verbosely) write(*,*)         

    deallocate(chem_id, net_iso)

  end subroutine test_irradiated

  !****

  subroutine set_composition()
    real(dp) :: xz, z2bar, ye, mass_correction, sumx, &
         dabar_dx(species), dzbar_dx(species), dmc_dx(species)
    integer :: i
    real(dp) :: norm
    type (Kap_General_Info), pointer :: rq
    call kap_ptr(kap_handle,rq,ierr)
    rq% Zbase = Z
    Y = 1-(X+Z)
    allocate(chem_id(species), net_iso(num_chem_isos), stat=ierr)
    if (ierr /= 0) then
       write(*,*) 'allocate failed'
       call mesa_error(__FILE__,__LINE__)
    end if
    chem_id(:) = (/ ih1, ihe4, ic12, in14, io16, ine20, img24 /)
    net_iso(:) = 0
    forall (i=1:species) net_iso(chem_id(i)) = i
    xa(:) = (/ X, Y, xc, xn, xo, 0d0, 0d0 /)
    xa(species) = 1 - sum(xa(:))
 
    norm = 0d0
    do i=1,species
      xa(i) = max(0d0, xa(i))
      norm = norm + xa(i)
    end do
    do i=1,species
      xa(i) = xa(i) / norm
    end do

    call composition_info( &
         species, chem_id, xa, X, Y, xz, abar, zbar, z2bar, z53bar, ye, &
         mass_correction, sumx, dabar_dx, dzbar_dx, dmc_dx)
  end subroutine set_composition

  !****

  subroutine eos_proc( &
       lnP, lnT, &
       lnRho, res, dres_dlnRho, dres_dlnT, &
       ierr)
    
    use eos_def, only: num_eos_basic_results
    use eos_lib, only: eosPT_get, radiation_pressure

    real(dp), intent(in)  :: lnP
    real(dp), intent(in)  :: lnT
    real(dp), intent(out) :: lnRho
    real(dp), intent(out) :: res(:)
    real(dp), intent(out) :: dres_dlnRho(:)
    real(dp), intent(out) :: dres_dlnT(:)
    integer, intent(out)  :: ierr

    real(dp) :: T, P, Prad, Pgas, logPgas, rho
    real(dp) :: logRho, dlnRho_dlnPgas, dlnRho_dlnT
    real(dp), dimension(num_eos_d_dxa_results, species) :: dres_dxa

    T = exp(lnT)
    P = exp(lnP)

    Prad = radiation_pressure(T)
    Pgas = max(1E-99_dp, P - Prad)
    logPgas = log10(Pgas)

    call eosPT_get( &
         eos_handle, &
         species, chem_id, net_iso, xa, &
         Pgas, logPgas, T, lnT/ln10, &
         Rho, logRho, dlnRho_dlnPgas, dlnRho_dlnT, &
         res, dres_dlnRho, dres_dlnT, dres_dxa, ierr)

    lnRho = logRho*ln10

  end subroutine eos_proc

  !****

  subroutine kap_proc( &
       lnRho, lnT, res, dres_dlnRho, dres_dlnT, &
       kap, dlnkap_dlnRho, dlnkap_dlnT, &
       ierr)

    use kap_def, only: num_kap_fracs
    use kap_lib, only: kap_get
    use eos_def, only: i_lnfree_e, i_eta

    real(dp), intent(in)  :: lnRho
    real(dp), intent(in)  :: lnT
    real(dp), intent(in)  :: res(:)
    real(dp), intent(in)  :: dres_dlnRho(:)
    real(dp), intent(in)  :: dres_dlnT(:)
    real(dp), intent(out) :: kap
    real(dp), intent(out) :: dlnkap_dlnRho
    real(dp), intent(out) :: dlnkap_dlnT
    integer, intent(out)  :: ierr

    real(dp) :: kap_fracs(num_kap_fracs), dlnkap_dxa(species)

    call kap_get( &
         kap_handle, species, chem_id, net_iso, xa, &
         lnRho/ln10, lnT/ln10, res(i_lnfree_e), dres_dlnRho(i_lnfree_e), dres_dlnT(i_lnfree_e), &
         res(i_eta), dres_dlnRho(i_eta), dres_dlnT(i_eta), &
         kap_fracs, kap, dlnkap_dlnRho, dlnkap_dlnT, dlnkap_dxa, ierr)
    
  end subroutine kap_proc

end module test_atm_support




