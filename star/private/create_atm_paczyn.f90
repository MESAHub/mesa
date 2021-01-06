! ***********************************************************************
!
!   Copyright (C) 2011-2019  Bill Paxton & The MESA Team
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
!
! ***********************************************************************

module create_atm_paczyn

  ! Uses

  use star_private_def
  use const_def
  use chem_def, only: num_chem_isos
  use atm_lib
  use atm_def
  use num_lib
  use mtx_lib
  use star_utils, only: get_XYZ

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: r_X = 1
  integer, parameter :: r_Z = 2
  integer, parameter :: r_rho = 3
  integer, parameter :: r_lntau = 4
  integer, parameter :: r_L = 5
  integer, parameter :: r_M = 6
  integer, parameter :: r_cgrav = 7
  integer, parameter :: r_rsurf = 8
  integer, parameter :: r_lnTeff = 9
  integer, parameter :: atm_lrpar = 9
  
  integer, parameter :: i_id = 1
  integer, parameter :: i_save_structure_info = 2
  integer, parameter :: atm_lipar = 2

  integer, parameter :: var_lnT = 1
  integer, parameter :: var_lnP = 2
  integer, parameter :: var_lnr = 3
  integer, parameter :: var_m = 4
  
  integer, parameter :: num_vars = 4

  real(dp), parameter :: lntau_Teff = ln2 - ln3

  ! Access specifiers

  private

  public :: get_Paczynski_atm_surf_PT
  public :: do_create_Paczynski_atm
  public :: get_mlt_basics ! someday move this out to mlt where it belongs

  ! Routines

contains

  subroutine get_Paczynski_atm_surf_PT( &
       s, R_surf, L_surf, skip_partials, Teff, &
       lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
       lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
       ierr)

    type (star_info), pointer :: s
    real(dp), intent(in) :: R_surf, L_surf
    logical, intent(in) :: skip_partials
    real(dp), intent(out) :: Teff, &
         lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
         lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap
    integer, intent(out) :: ierr

    logical, parameter :: save_structure_info = .false.
    real(dp) :: delta_L_surf, delta_R_surf, &
         lnP_surf1, lnT_surf1, Teff1, &
         lnP_surf2, lnT_surf2, Teff2, &
         lnP_surf3, lnT_surf3, Teff3, &
         lnP_surf4, lnT_surf4, Teff4, &
         R0_div_R, R0_div_R1, R0_div_R2, R0_div_R3, R0_div_R4
    integer :: j, op_err
    logical :: trace
    logical, parameter :: dbg = .false.

    include 'formats'

    ierr = 0
    dlnT_dlnM = 0
    dlnT_dlnkap = 0
    dlnP_dlnM = 0
    dlnP_dlnkap = 0
    trace = s% trace_atm_Paczynski_grey

    if (skip_partials) then
       call do_create_Paczynski_atm( &
            s, save_structure_info, trace, &
            L_surf, R_surf, lnP_surf, lnT_surf, Teff, R0_div_R, ierr)
       if (ierr == 0) &
            s% prev_create_atm_R0_div_R = R0_div_R ! save hint for next time
       dlnT_dL = 0
       dlnT_dlnR = 0
       dlnP_dL = 0
       dlnP_dlnR = 0
       return
    end if

    delta_L_surf = 1d-7*L_surf
    delta_R_surf = 1d-7*R_surf

    !$OMP PARALLEL DO PRIVATE(j,op_err) SCHEDULE(dynamic,2)
    do j=0,4
       op_err = 0
       select case(j)
       case (0)
          call do_create_Paczynski_atm( &
               s, save_structure_info, trace, L_surf, R_surf, &
               lnP_surf, lnT_surf, Teff, R0_div_R, op_err)
       case (1)
          call do_create_Paczynski_atm( &
               s, .false., .false., L_surf + delta_L_surf, R_surf, &
               lnP_surf1, lnT_surf1, Teff1, R0_div_R1, op_err)
       case (2)
          call do_create_Paczynski_atm( &
               s, .false., .false., L_surf - delta_L_surf, R_surf, &
               lnP_surf2, lnT_surf2, Teff2, R0_div_R2, op_err)
       case (3)
          call do_create_Paczynski_atm( &
               s, .false., .false., L_surf, R_surf + delta_R_surf, &
               lnP_surf3, lnT_surf3, Teff3, R0_div_R3, op_err)
       case (4)
          call do_create_Paczynski_atm( &
               s, .false., .false., L_surf, R_surf - delta_R_surf, &
               lnP_surf4, lnT_surf4, Teff4, R0_div_R4, op_err)
       end select
       if (op_err /= 0) ierr = op_err
    end do
    !$OMP END PARALLEL DO

    if (ierr /= 0) then
       if (dbg) write(*,2) 'get_Paczynski_atm_surf_PT ierr', ierr
       return
    end if

    s% prev_create_atm_R0_div_R = R0_div_R ! save hint for next time

    dlnP_dL = (lnP_surf1 - lnP_surf2)/(2*delta_L_surf)
    dlnT_dL = (lnP_surf1 - lnP_surf2)/(2*delta_L_surf)

    dlnP_dlnR = (lnP_surf3 - lnP_surf4)*R_surf/(2*delta_R_surf)
    dlnT_dlnR = (lnT_surf3 - lnT_surf4)*R_surf/(2*delta_R_surf)

  end subroutine get_Paczynski_atm_surf_PT

  !****

  ! create an atmosphere for given base conditions.
  ! inspired by B. Paczynski, 1969, Acta Astr., vol. 19, 1.
  ! takes into account dilution when tau < 2/3,
  ! and calls mlt to get gradT allowing for convection.

  subroutine do_create_Paczynski_atm ( &
       s, save_structure_info, trace, &
       L_surf, R_surf, lnP_surf, lnT_surf, Teff, R0_div_R, ierr)

    type (star_info), pointer :: s
    logical, intent(in) :: save_structure_info, trace
    real(dp), intent(in) :: L_surf, R_surf
    real(dp), intent(out) :: lnP_surf, lnT_surf, Teff, R0_div_R
    integer, intent(out) :: ierr

    real(dp) :: delr, R0, errtol, integration_R_surf, prev_create_atm_R0_div_R
    integer :: i
    integer, parameter :: max_tries = 10
    logical :: okay
    logical, parameter :: dbg = .false.

    include 'formats'

    ierr = 0
    prev_create_atm_R0_div_R = s% prev_create_atm_R0_div_R
    if (prev_create_atm_R0_div_R <= 0) prev_create_atm_R0_div_R = 1.01d0
    R0 = R_surf*prev_create_atm_R0_div_R
    R0_div_R = prev_create_atm_R0_div_R
    errtol = s% Paczynski_atm_R_surf_errtol
    okay = .false.
    integration_R_surf = 0 ! to keep gfortran quiet

    do i=1,max_tries
       call do1_create_Paczynski_atm( &
            s, save_structure_info, R0, L_surf, R_surf, &
            integration_R_surf, lnP_surf, lnT_surf, Teff, ierr)
       if (ierr /= 0) then
          write(*,*) 'failed in do1_create_Paczynski_atm'
          return
          stop 'do_create_Paczynski_atm'
       end if
       delr = integration_R_surf - R_surf
       if (dbg .or. trace) then
          write(*,*)
          write(*,1) 'R0/Rsun', R0/Rsun
          write(*,1) 'integration_R_surf/Rsun', (R_surf+delr)/Rsun
          write(*,1) 'R_surf/Rsun', R_surf/Rsun
          write(*,1) 'delr/Rsun', delr/Rsun
          write(*,1) 'errtol', errtol
          write(*,1) 'abs(delr)/R_surf', abs(delr)/R_surf
          write(*,*)
       end if
       ! want to adjust R0 to get abs(delr) < errtol*R_surf
       R0 = R0 - delr
       if (abs(delr) < errtol*R_surf) then
          R0_div_R = R0/R_surf ! save hint for next time
          okay = .true.
          exit
       end if
       if (dbg) write(*,2) 'try again', i
    end do

    if (.not. okay) then
       ierr = -1
       write(*,*) 'do_create_Paczynski_atm failed to accept'
       return

       stop 'do_create_Paczynski_atm'
    end if

  end subroutine do_create_Paczynski_atm

  !****

  subroutine do1_create_Paczynski_atm( &
       s, save_structure_info, R0, L_surf, R_surf, &
       integration_R_surf, lnP_surf, lnT_surf, Teff, ierr)

    type (star_info), pointer :: s
    logical, intent(in) :: save_structure_info
    real(dp), intent(in) :: R0, L_surf, R_surf
    real(dp), intent(out) :: integration_R_surf, lnP_surf, lnT_surf, Teff
    integer, intent(out) :: ierr

    integer, parameter :: nrdens = num_vars
    real(dp) :: P0, rho0, T0, X, Y, Z, lntau_stop, lnTeff, x_min, x_max, &
         tau_surface, lntau_start, errtol, init_step_size, max_step_size

    real(dp) :: rtol(1) ! relative error tolerance(s)
    real(dp) :: atol(1) ! absolute error tolerance(s)
    integer, parameter :: itol = 0 ! switch for rtol and atol
    integer, parameter :: lout = 0 ! set to 6 for debugging

    integer :: iout
    ! iout = 0 don't call solout
    ! iout = 1 call solout without info for interpolating
    ! iout = 2 call solout with info for interpolating

    integer, pointer :: iwork(:), ipar_decsol(:), ipar(:)
    real(dp), pointer :: work(:), rpar_decsol(:), rpar(:)

    integer :: liwork, lwork, max_steps, idid, lipar, lrpar, lrd, lid
    integer :: ijac, mljac, mujac, imas, mlmas, mumas, nzmax, which_solver

    logical, parameter :: use_dopri5 = .false.
    logical, parameter :: use_dopri853 = .true.

    real(dp), target :: vars_ary(num_vars)
    real(dp), pointer :: vars(:)

    integer :: caller_id, nvar_blk, nz_blk
    real(dp), dimension(:), pointer :: lblk, dblk, ublk ! =(nvar,nvar,nz)

    logical, parameter :: dbg = .false.

    include 'formats'

    vars => vars_ary

    ierr = 0

    lipar = atm_lipar
    lrpar = atm_lrpar
    lid = 0
    lrd = 0

    if (use_dopri5) then
       call dopri5_work_sizes(num_vars, nrdens, liwork, lwork)
    else if (use_dopri853) then
       call dop853_work_sizes(num_vars, nrdens, liwork, lwork)
    else
       stop 'debug: need to pick which solver for create_atm'
    end if
    allocate(ipar(lipar), rpar(lrpar), iwork(liwork), work(lwork), &
         ipar_decsol(lid), rpar_decsol(lrd), stat=ierr)
    if (ierr /= 0) then
       write(*,*) 'allocate failed in do1_create_Paczynski_atm'
       return
    end if

    lntau_stop = log(s% tau_factor*s% tau_base)
    lntau_start = -4d0*ln10 ! log(1d-4)

    if (lntau_stop > lntau_Teff) then
       iout = 2 ! interpolate lnTeff
    else
       iout = 1 ! don't need to interpolate
    end if

    call get_XYZ(s, s% xa(:,1), X, Y, Z)
    rpar(r_X) = X
    rpar(r_Z) = Z
    rpar(r_L) = L_surf
    rpar(r_cgrav) = s% cgrav(1)
    rpar(r_M) = s% m_grav(1)
    rpar(r_rsurf) = R_surf
    rpar(r_lnTeff) = missing_value

    rho0 = 1d-12
    if (s% lnT(1)/ln10 > 5) rho0 = exp10((2*s% lnT(1)/ln10 - 22))
    ! eos requires logQ = logRho - 2*logT + 12 > -10
    T0 = pow(L_surf/(8*pi*boltz_sigma*R0*R0),0.25d0)
    P0 = get_P0(ierr)
    if (ierr /= 0) then
       call dealloc
       return
    end if

    max_step_size = s% atm_build_dlogtau*ln10
    init_step_size = 0.5d0*max_step_size

    ipar(i_id) = s% id
    if (save_structure_info) then
       ipar(i_save_structure_info) = 1
    else
       ipar(i_save_structure_info) = 0
    end if

    max_steps = 1000

    iwork(:) = 0
    work(:) = 0
    iwork(5) = nrdens

    errtol = 1d-7
    rtol(:) = errtol
    atol(:) = errtol

    x_min = -1d99
    x_max = 1d99

    vars(var_lnT) = log(T0)
    vars(var_lnP) = log(P0)
    vars(var_lnr) = log(R0)
    vars(var_m) = 0d0

    s% atm_structure_num_pts = 0
    if (associated(s% atm_structure)) deallocate(s% atm_structure)
    nullify(s% atm_structure, lblk, dblk, ublk)
    caller_id = 0
    nvar_blk = 0
    nz_blk = 0

    if (use_dopri5) then
       call dopri5( &
            num_vars, create_atm_fcn, lntau_start, vars, lntau_stop, &
            init_step_size, max_step_size, max_steps, &
            rtol, atol, itol, &
            create_atm_solout, iout, &
            work, lwork, iwork, liwork, &
            lrpar, rpar, lipar, ipar, &
            lout, idid)
       if (idid < 0) write(*,*) 'do1_create_Paczynski_atm failed in dopri5: idid', idid
    else if (use_dopri853) then
       call dop853( &
            num_vars, create_atm_fcn, lntau_start, vars, lntau_stop, &
            init_step_size, max_step_size, max_steps, &
            rtol, atol, itol, &
            create_atm_solout, iout, &
            work, lwork, iwork, liwork, &
            lrpar, rpar, lipar, ipar, &
            lout, idid)
       if (idid < 0) write(*,*) 'do1_create_Paczynski_atm failed in dopri853: idid', idid
    else
       stop 'debug: need to pick which solver for do1_create_Paczynski_atm'
    end if

    if (idid < 0) then
       ierr = -1
       call dealloc
       return
    end if

    integration_R_surf = exp(vars(var_lnr))
    lnP_surf = vars(var_lnP)
    lnT_surf = vars(var_lnT)
    lnTeff = rpar(r_lnTeff)
    if (lnTeff <= 0d0) lnTeff = lnT_surf
    Teff = exp(lnTeff)
    
    !write(*,2) 'atm dM/Msun, R0/Rsun, R(1)/Rsun', s% model_number, &
    !  vars(var_m)/Msun, R0/Rsun, s% r(1)/Rsun

    if (dbg) then
       write(*,*)
       write(*,1) 'results for do1_create_Paczynski_atm'
       write(*,2) 'atm_structure_num_pts', s% atm_structure_num_pts
       write(*,1) 'at bottom of atmosphere: tau', exp(rpar(r_lntau))
       write(*,1) 'logRho', log10(rpar(r_rho))
       write(*,1) 'logT_surf', lnT_surf/ln10
       write(*,1) 'logP_surf', lnP_surf/ln10
       write(*,1) 'integration_R_surf', integration_R_surf
       write(*,1) 'Teff', Teff
       write(*,*)
    end if

    call dealloc()

  contains

    subroutine dealloc ()
      deallocate(rpar, ipar, iwork, work, ipar_decsol, rpar_decsol)
    end subroutine dealloc

         real(dp) function get_P0(ierr)
            use eos_support, only: get_eos
            use eos_def
            integer, intent(out) :: ierr
            real(dp) :: Pgas, Prad, &
               res(num_eos_basic_results), &
               d_eos_dlnd(num_eos_basic_results), &
               d_eos_dlnT(num_eos_basic_results), &
               d_eos_dabar(num_eos_basic_results), &
               d_eos_dzbar(num_eos_basic_results), &
               dres_dxa(num_eos_basic_results,s% species)
            call get_eos( &
               s, 0, s% xa(:,1), rho0, log10(rho0), T0, log10(T0), &
               res, d_eos_dlnd, d_eos_dlnT, &
               dres_dxa, ierr)
            if (ierr /= 0) then
               write(*,*) 'get_P0 failed in get_eos'
               stop
               return
            end if
            Pgas = exp(res(i_lnPgas))
            Prad = crad * T0*T0*T0*T0 / 3
            get_P0 = Prad + Pgas
         end function get_P0

  end subroutine do1_create_Paczynski_atm

  !****

  subroutine create_atm_solout( &
       nr, lntau_prev, lntau, n, vars, rwork_y, iwork_y, interp_y, &
       lrpar, rpar, lipar, ipar, irtrn)
    ! nr is the step number.
    ! lntau is the current value; lntau_prev is the previous lntau value.
    ! vars has the current variables.
    ! irtrn negative means terminate integration.
    ! rwork_y and iwork_y hold info for interp_y
    integer, intent(in) :: nr, n, lrpar, lipar
    real(dp), intent(in) :: lntau_prev, lntau
    real(dp), intent(inout) :: vars(:) ! (n)
    ! y can be modified if necessary to keep it in valid range of possible solutions.
    real(dp), intent(inout), target :: rwork_y(*)
    integer, intent(inout), target :: iwork_y(*)
    integer, intent(inout), pointer :: ipar(:) ! (lipar)
    real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
    interface
       ! this subroutine can be called from your solout routine.
       ! it computes interpolated values for y components during the just completed step.
       real(dp) function interp_y(i, s, rwork_y, iwork_y, ierr)
         use const_def, only: dp
         integer, intent(in) :: i ! result is interpolated approximation of y(i) at x=s.
         real(dp), intent(in) :: s ! interpolation x value (between xold and x).
         real(dp), intent(inout), target :: rwork_y(*)
         integer, intent(inout), target :: iwork_y(*)
         integer, intent(out) :: ierr
       end function interp_y
    end interface
    integer, intent(out) :: irtrn ! < 0 causes solver to return to calling program.
    integer :: ierr
    logical, parameter :: dbg = .false.
    logical :: save_structure_info
    real(dp) :: derivs(n)
    include 'formats'
    irtrn = 0
    rpar(r_lntau) = lntau
    if ( .false. .and. abs(lntau_Teff - lntau) < 0.1d0*abs(lntau_Teff)) then
       write(*,2) 'create_atm_solout: step, log tau, T, log P', &
            nr, lntau/ln10, exp(vars(var_lnT)), vars(var_lnP)/ln10
    end if
    ierr = 0
    save_structure_info = (ipar(i_save_structure_info) /= 0)
    call get_atm_fcn_info( &
         n, lntau, vars, derivs, lrpar, rpar, lipar, ipar, save_structure_info, ierr)
    if (ierr /= 0) then
       write(*,*) 'failed in create_atm_solout'
       stop 'create_atm_solout'
       irtrn = -1
    end if
    if (dbg) then
       write(*,2) 'lntau_prev', nr, lntau_prev
       write(*,2) 'lntau_Teff', nr, lntau_Teff
       write(*,2) 'lntau', nr, lntau
       write(*,*) 'lntau_prev <= lntau_Teff', lntau_prev <= lntau_Teff
       write(*,*) 'lntau_Teff < lntau', lntau_Teff < lntau
    end if
    if (lntau_prev <= lntau_Teff .and. lntau_Teff < lntau) then ! set lnTeff
       rpar(r_lnTeff) = interp_y(var_lnT, lntau_Teff, rwork_y, iwork_y, ierr)
       if (dbg) write(*,2) 'Teff T', nr, exp(rpar(r_lnTeff)), exp(vars(var_lnT))
       if (ierr /= 0) then
          write(*,*) 'interp_y failed in create_atm_solout'
          stop 'create_atm_solout'
          irtrn = -1
       end if
    end if
  end subroutine create_atm_solout

  !****
  
  subroutine create_atm_fcn(n, lntau, h, vars, derivs, lr, rpar, li, ipar, ierr)
    integer, intent(in) :: n, lr, li
    real(dp), intent(in) :: lntau, h
    real(dp), intent(inout) :: vars(:) ! (n)
    real(dp), intent(inout) :: derivs(:) ! (n)
    integer, intent(inout), pointer :: ipar(:) ! (li)
    real(dp), intent(inout), pointer :: rpar(:) ! (lr)
    integer, intent(out) :: ierr ! nonzero means retry with smaller step.
    logical, parameter :: save_structure_info = .false.
    call get_atm_fcn_info( &
         n, lntau, vars, derivs, lr, rpar, li, ipar, save_structure_info, ierr)
  end subroutine create_atm_fcn

  !****
  
      subroutine get_atm_fcn_info( &
            n, lntau, vars, derivs, lr, rpar, li, ipar, save_structure_info, ierr)
         use utils_lib, only: realloc_double2
         use eos_def
         use eos_support, only: solve_eos_given_PgasT_auto
         use kap_def, only: num_kap_fracs
         use kap_support, only: get_kap

         integer, intent(in) :: n, lr, li
         real(dp), intent(in) :: lntau
         real(dp), intent(inout) :: vars(:) ! (n)
         real(dp) :: derivs(:) ! (n) ! dvars/dlntau
         integer, intent(inout), pointer :: ipar(:) ! (li)
         real(dp), intent(inout), pointer :: rpar(:) ! (lr)
         logical, intent(in) :: save_structure_info
         integer, intent(out) :: ierr ! nonzero means retry with smaller step.

         real(dp), parameter :: LOGPGAS_TOL = 1E-6_dp
         real(dp), parameter :: LOGRHO_TOL = 1E-6_dp

         type (star_info), pointer :: s
         integer :: sz, k, id
         real(dp) :: &
            Z, X, tau, logRho, logT, logP, logPgas, lnm, m, mstar, lnr, r_surf, delr, cgrav, &
            r, L, rho, T, Pgas, Prad, P, opacity, chiRho, chiT, Cp, Cv, grada, gradr_factor, &
            sfactor, f, dm_dlntau, dlnr_dlntau, &
            dlnP_dm, dlnP_dlntau, gradT, dlnT_dlntau, max_conv_vel, dt, &
            kap_fracs(num_kap_fracs), dlnkap_dlnd, dlnkap_dlnT, res(num_eos_basic_results), &
            d_eos_dlnd(num_eos_basic_results), d_eos_dlnT(num_eos_basic_results), &
            d_eos_dabar(num_eos_basic_results), d_eos_dzbar(num_eos_basic_results), &
            mlt_basics(num_mlt_results)
         real(dp), pointer :: xa(:)

         logical, parameter :: dbg = .false.

         include 'formats'

         ierr = 0
         id = ipar(1)
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'bad id for get_atm_fcn_info', id
            stop
            return
         end if

         Z = rpar(r_Z)
         X = rpar(r_X)
         L = rpar(r_L)
         cgrav = rpar(r_cgrav)
         m = rpar(r_M)
         mstar = rpar(r_M)
         r_surf = rpar(r_rsurf)

         tau = exp(lntau)
         lnr = vars(var_lnr)

         logT = vars(var_lnT)/ln10
         logP = vars(var_lnP)/ln10

         P = exp10(logT)
         T = exp10(logP)
         r = exp(lnr)

         Prad = crad * T*T*T*T / 3
         Pgas = P - Prad
         if (Pgas <= 0) then
            ierr = -1
            if (dbg) then
               write(*,*) 'Pgas <= 0'
               stop
            end if
            return
         end if
         logPgas = log10(Pgas)

         xa(1:s% species) => s% xa(1:s% species,1)
         call solve_eos_given_PgasT_auto( &
              s, 0, Z, X, s% abar(1), s% zbar(1), xa, &
              logT, logPgas, LOGRHO_TOL, LOGPGAS_TOL, &
              logRho, res, d_eos_dlnd, d_eos_dlnT, d_eos_dabar, d_eos_dzbar, &
              ierr)
         if (ierr /= 0) then
            if (dbg) then
               write(*,*) 'failed in solve_eos_given_PgasT_auto'
               stop
            end if
            return
         end if
         rho = exp10(logRho)
         rpar(r_rho) = rho

         call get_kap( &
            s, 1, s% zbar(1), xa, logRho, logT, &
            s% lnfree_e(1), s% d_eos_dlnd(i_lnfree_e,1), s% d_eos_dlnT(i_lnfree_e,1), &
            s% eta(1), s% d_eos_dlnd(i_eta,1), s% d_eos_dlnT(i_eta,1), &
            kap_fracs, opacity, dlnkap_dlnd, dlnkap_dlnT, ierr)
         if (ierr /= 0) then
            if (dbg) then
               write(*,*) 'failed in get_kap'
               stop
            end if
            return
         end if

         chiRho = res(i_chiRho)
         chiT = res(i_chiT)
         Cp = res(i_cp)
         Cv = res(i_Cv)
         grada = res(i_grad_ad)
            
         call get_mlt_basics(s, &
            cgrav, m, mstar, r, L, X, T, rho, P, chiRho, chiT, &
            Cp, opacity, grada, tau, mlt_basics, ierr)
            
         if (ierr /= 0) then
            if (dbg) then
               write(*,*) 'failed in get_mlt_basics'
               stop
            end if
            return
         end if

         if (tau < 2d0/3d0) then
            sfactor = (2d0*crad*T*T*T*sqrt(r))/&
                        (3d0*cgrav*m*rho)*pow(L/(8d0*pi*boltz_sigma),0.25d0)
            f = 1d0 - 1.5d0*tau
         else
            sfactor = 0
            f = 0
         end if

         dm_dlntau = -4*pi*r*r*tau/opacity

         dlnP_dm = -(cgrav*m/(4*pi*r*r*r*r*P))*(1 + f*sfactor) ! Paczynski, 1969, eqn 33
         dlnP_dlntau = dlnP_dm*dm_dlntau

         gradT = mlt_basics(mlt_gradT)

         dlnT_dlntau = gradT*dlnP_dlntau
         dlnr_dlntau = -tau/(r*rho*opacity)

         derivs(var_lnT) = dlnT_dlntau
         derivs(var_lnP) = dlnP_dlntau
         derivs(var_lnr) = dlnr_dlntau
         derivs(var_m) = -dm_dlntau

         if (.not. save_structure_info) return

         k = s% atm_structure_num_pts + 1
         if (.not. associated(s% atm_structure)) then
            sz = 100
            allocate(s% atm_structure(num_results_for_build_atm,sz))
         else
            sz = size(s% atm_structure,dim=2)
            if (k >= sz) then
               sz = 2*sz + 100
               call realloc_double2( &
                  s% atm_structure,num_results_for_build_atm,sz,ierr)
            end if
         end if

         s% atm_structure_num_pts = k

         s% atm_structure(atm_delta_r,k) = r - r_surf
         s% atm_structure(atm_lnP,k) = logP*ln10
         s% atm_structure(atm_lnd,k) = logRho*ln10
         s% atm_structure(atm_lnT,k) = logT*ln10
         s% atm_structure(atm_gradT,k) = gradT
         s% atm_structure(atm_kap,k) = opacity
         s% atm_structure(atm_gamma1,k) = res(i_gamma1)
         s% atm_structure(atm_grada,k) = res(i_grad_ad)
         s% atm_structure(atm_chiT,k) = chiT
         s% atm_structure(atm_chiRho,k) = chiRho
         s% atm_structure(atm_cp,k) = cp
         s% atm_structure(atm_cv,k) = cv
         s% atm_structure(atm_tau,k) = tau
         s% atm_structure(atm_lnfree_e,k) = res(i_lnfree_e)
         s% atm_structure(atm_dlnkap_dlnT,k) = dlnkap_dlnT
         s% atm_structure(atm_dlnkap_dlnd,k) = dlnkap_dlnd
         s% atm_structure(atm_lnPgas,k) = log(Pgas)
         s% atm_structure(atm_gradr,k) = mlt_basics(mlt_gradr)

         if (k == 1) return
         if (s% atm_structure(atm_lnP,k) > s% atm_structure(atm_lnP,k-1)) return

         write(*,2) 'pressure inversion k', k
         write(*,2) 's% atm_structure(atm_lnP,k)', k, s% atm_structure(atm_lnP,k)
         write(*,2) 's% atm_structure(atm_lnP,k-1)', k-1, s% atm_structure(atm_lnP,k-1)
         write(*,1) 'lnP(k) - lnP(k-1)', s% atm_structure(atm_lnP,k) - s% atm_structure(atm_lnP,k-1)
         write(*,2) 'logT', k, logT
         write(*,2) 'logRho', k, logRho
         write(*,2) 'delta_r/Rsun', k, s% atm_structure(atm_delta_r,k)/Rsun
         write(*,*)

      end subroutine get_atm_fcn_info
      
      
      subroutine get_mlt_basics(s, &
            cgrav, m, mstar, r, L, X, T, rho, P, chiRho, chiT, &
            Cp, opacity, grada, tau, mlt_basics, ierr)
         use mlt_info, only: do1_mlt_eval
         type (star_info), pointer :: s
         real(dp), intent(in) :: cgrav, m, mstar, r, L, X, T, &
            rho, P, chiRho, chiT, Cp, opacity, grada, tau
         real(dp), intent(out) :: mlt_basics(num_mlt_results)
         integer, intent(out) :: ierr
         
         real(dp) :: &
            gradr_factor, d_gradr_factor_dw, gradL_composition_term, normal_mlt_gradT_factor, &
            max_conv_vel, dt, csound, &
            alfa, beta, &
            T_00, T_m1, rho_00, rho_m1, P_00, P_m1, &
            chiRho_for_partials_00, chiT_for_partials_00, &
            chiRho_for_partials_m1, chiT_for_partials_m1, &
            chiRho_00, d_chiRho_00_dlnd, d_chiRho_00_dlnT, &
            chiRho_m1, d_chiRho_m1_dlnd, d_chiRho_m1_dlnT, &
            chiT_00, d_chiT_00_dlnd, d_chiT_00_dlnT, &
            chiT_m1, d_chiT_m1_dlnd, d_chiT_m1_dlnT, &
            Cp_00, d_Cp_00_dlnd, d_Cp_00_dlnT, &
            Cp_m1, d_Cp_m1_dlnd, d_Cp_m1_dlnT, &
            opacity_00, d_opacity_00_dlnd, d_opacity_00_dlnT, &
            opacity_m1, d_opacity_m1_dlnd, d_opacity_m1_dlnT, &
            grada_00, d_grada_00_dlnd, d_grada_00_dlnT, &
            grada_m1, d_grada_m1_dlnd, d_grada_m1_dlnT         
         real(dp), target :: mlt_partials1_ary(num_mlt_partials*num_mlt_results)
         real(dp), pointer :: mlt_partials1(:), mlt_partials(:,:)
         integer :: mixing_type
         
         mlt_partials1 => mlt_partials1_ary
         mlt_partials(1:num_mlt_partials,1:num_mlt_results) => &
            mlt_partials1(1:num_mlt_partials*num_mlt_results)

         if (s% rotation_flag) then
            gradr_factor = s% ft_rot(1)/s% fp_rot(1)
         else
            gradr_factor = 1
         end if
         gradL_composition_term = 0
         d_gradr_factor_dw = 0d0

         normal_mlt_gradT_factor = 1d0 ! when .not. s% conv_vel_flag
         
         max_conv_vel = 1d99
         dt = -1
         
         csound = 0 ! not used when dt <= 0
         ! blend info not used when k == 0
         alfa=0d0; beta=0d0
         T_00=0d0; T_m1=0d0; rho_00=0d0; rho_m1=0d0; P_00=0d0; P_m1=0d0
         chiRho_for_partials_00=0d0; chiT_for_partials_00=0d0
         chiRho_for_partials_m1=0d0; chiT_for_partials_m1=0d0
         chiRho_00=0d0; d_chiRho_00_dlnd=0d0; d_chiRho_00_dlnT=0d0
         chiRho_m1=0d0; d_chiRho_m1_dlnd=0d0; d_chiRho_m1_dlnT=0d0
         chiT_00=0d0; d_chiT_00_dlnd=0d0; d_chiT_00_dlnT=0d0
         chiT_m1=0d0; d_chiT_m1_dlnd=0d0; d_chiT_m1_dlnT=0d0
         Cp_00=0d0; d_Cp_00_dlnd=0d0; d_Cp_00_dlnT=0d0
         Cp_m1=0d0; d_Cp_m1_dlnd=0d0; d_Cp_m1_dlnT=0d0
         opacity_00=0d0; d_opacity_00_dlnd=0d0; d_opacity_00_dlnT=0d0
         opacity_m1=0d0; d_opacity_m1_dlnd=0d0; d_opacity_m1_dlnT=0d0
         grada_00=0d0; d_grada_00_dlnd=0d0; d_grada_00_dlnT=0d0
         grada_m1=0d0; d_grada_m1_dlnd=0d0; d_grada_m1_dlnT=0d0            

         call do1_mlt_eval( &
            s, 0, cgrav, m, mstar, r, L, X, T, rho, P, & 
            chiRho, chiT, Cp, opacity, grada, &
            
            ! not used
               alfa, beta, &
               T_00, T_m1, rho_00, rho_m1, P_00, P_m1, &
               chiRho_for_partials_00, chiT_for_partials_00, &
               chiRho_for_partials_m1, chiT_for_partials_m1, &
               chiRho_00, d_chiRho_00_dlnd, d_chiRho_00_dlnT, &
               chiRho_m1, d_chiRho_m1_dlnd, d_chiRho_m1_dlnT, &
               chiT_00, d_chiT_00_dlnd, d_chiT_00_dlnT, &
               chiT_m1, d_chiT_m1_dlnd, d_chiT_m1_dlnT, &
               Cp_00, d_Cp_00_dlnd, d_Cp_00_dlnT, &
               Cp_m1, d_Cp_m1_dlnd, d_Cp_m1_dlnT, &
               opacity_00, d_opacity_00_dlnd, d_opacity_00_dlnT, &
               opacity_m1, d_opacity_m1_dlnd, d_opacity_m1_dlnT, &
               grada_00, d_grada_00_dlnd, d_grada_00_dlnT, &
               grada_m1, d_grada_m1_dlnd, d_grada_m1_dlnT, &            
            
            
            gradr_factor, d_gradr_factor_dw, gradL_composition_term, &
            s% alpha_semiconvection, s% semiconvection_option, &
            s% thermohaline_coeff, s% thermohaline_option, &
            s% dominant_iso_for_thermohaline(1), &
            s% mixing_length_alpha, s% alt_scale_height_flag, s% remove_small_D_limit, &
            s% MLT_option, s% Henyey_MLT_y_param, s% Henyey_MLT_nu_param, &
            normal_mlt_gradT_factor, &
            max_conv_vel, dt, tau, .false., &
            mixing_type, mlt_basics, mlt_partials1, ierr)
      
      end subroutine get_mlt_basics
      

end module create_atm_paczyn

