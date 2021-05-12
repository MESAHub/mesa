! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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

module atm_T_tau_uniform

  ! Uses

  use const_def
  use math_lib
  use utils_lib, only: mesa_error
  use utils_lib, only: is_bad

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: eval_T_tau_uniform
  public :: build_T_tau_uniform

  ! Procedures

contains

  ! Evaluate atmosphere data from T-tau relation with uniform opacity

  subroutine eval_T_tau_uniform( &
       tau_surf, L, R, M, cgrav, kap_guess, Pextra_factor, &
       T_tau_id, eos_proc, kap_proc, errtol, max_iters, skip_partials, &
       Teff, kap, &
       lnT, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
       lnP, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
       ierr)

    use atm_def, only: atm_eos_iface, atm_kap_iface
    use atm_utils, only: eval_Teff_g
    use eos_def, only: num_eos_basic_results, i_chiRho, i_chiT

    real(dp), intent(in)       :: tau_surf
    real(dp), intent(in)       :: L
    real(dp), intent(in)       :: R
    real(dp), intent(in)       :: M
    real(dp), intent(in)       :: cgrav
    real(dp), intent(in)       :: kap_guess
    real(dp), intent(in)       :: Pextra_factor
    integer, intent(in)        :: T_tau_id
    procedure(atm_eos_iface)   :: eos_proc
    procedure(atm_kap_iface)   :: kap_proc
    real(dp), intent(in)       :: errtol
    integer, intent(in)        :: max_iters
    logical, intent(in)        :: skip_partials
    real(dp), intent(out)      :: Teff
    real(dp), intent(out)      :: kap
    real(dp), intent(out)      :: lnT
    real(dp), intent(out)      :: dlnT_dL
    real(dp), intent(out)      :: dlnT_dlnR
    real(dp), intent(out)      :: dlnT_dlnM
    real(dp), intent(out)      :: dlnT_dlnkap
    real(dp), intent(out)      :: lnP
    real(dp), intent(out)      :: dlnP_dL
    real(dp), intent(out)      :: dlnP_dlnR
    real(dp), intent(out)      :: dlnP_dlnM
    real(dp), intent(out)      :: dlnP_dlnkap
    integer, intent(out)       :: ierr

    real(dp) :: g
    integer  :: iters
    real(dp) :: lnRho
    real(dp) :: res(num_eos_basic_results)
    real(dp) :: dres_dlnRho(num_eos_basic_results)
    real(dp) :: dres_dlnT(num_eos_basic_results)
    real(dp) :: dlnkap_dlnRho
    real(dp) :: dlnkap_dlnT
    real(dp) :: kap_prev
    real(dp) :: err
    real(dp) :: chiRho
    real(dp) :: chiT
    real(dp) :: dlnkap_dlnP_T
    real(dp) :: dlnkap_dlnT_P
    real(dp) :: dlnP_dL_
    real(dp) :: dlnT_dL_
    real(dp) :: dlnP_dlnR_
    real(dp) :: dlnT_dlnR_
    real(dp) :: dlnP_dlnM_
    real(dp) :: dlnT_dlnM_
    
    include 'formats'

    ierr = 0

    ! Sanity checks

    if (L <= 0._dp .OR. R <= 0._dp .OR. M <= 0._dp) then
       ierr = -1
       write(*,*) 'atm: eval_T_tau_uniform: L, R, or M bad', L, R, M
       return
    end if

    ! Evaluate the effective temperature & gravity

    call eval_Teff_g(L, R, M, cgrav, Teff, g)

    ! Evaluate atmosphere data at optical depth tau_surf,
    ! using kap_guess as the opacity

    kap = kap_guess

    call eval_data( &
         tau_surf, Teff, g, L, M, cgrav, &
         kap, Pextra_factor, T_tau_id, skip_partials, &
         lnT, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
         lnP, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
         ierr)
    if (ierr /= 0) then
       write(*,*) 'atm: Call to eval_data failed in eval_T_tau_uniform'
       return
    end if
    if (is_bad(lnT)) then
       ierr = -1
       write(*,*) 'atm: eval_T_tau_uniform logT from eval_data', lnT/ln10
       return
    end if
    if (is_bad(lnP)) then
       ierr = -1
       write(*,*) 'atm: eval_T_tau_uniform logP from eval_data', lnP/ln10
       return
    end if

    ! Iterate to find a consistent opacity

    iterate_loop : do iters = 1, max_iters

       ! Calculate the density & eos results

       call eos_proc( &
            lnP, lnT, &
            lnRho, res, dres_dlnRho, dres_dlnT, &
            ierr)
       if (ierr /= 0) then
          write(*,*) 'atm: call to eos_proc failed in eval_T_tau_uniform logP logT logPrad', &
             lnP/ln10, lnT/ln10, log10(crad*exp(4d0*lnT)/3d0)
          return
       end if

       ! Update the opacity

       kap_prev = kap

       call kap_proc( &
            lnRho, lnT, res, dres_dlnRho, dres_dlnT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, &
            ierr)
       if (ierr /= 0) then
          write(*,*) 'atm: Call to kap_proc failed in eval_T_tau_uniform logRho logT', lnRho/ln10, lnT/ln10
          return
       end if

       ! Check for convergence

       err = abs(kap_prev - kap)/(errtol + errtol*kap)

       if (err < 1._dp) exit iterate_loop

       kap = kap_prev + 0.5_dp*(kap - kap_prev) ! under correct

       ! Re-evaluate atmosphere data

       call eval_data( &
            tau_surf, Teff, g, L, M, cgrav, &
            kap, Pextra_factor, T_tau_id, skip_partials, &
            lnT, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
            ierr)
       if (ierr /= 0) then
          write(*,*) 'Call to eval_data failed in eval_T_tau_uniform'
          return
       end if

    end do iterate_loop

    if (max_iters > 0 .AND. iters > max_iters) then
       write(*,*) 'atm: Exceeded max_iters iterations in eval_T_tau_uniform'
       ierr = -1
       return
    end if

    ! If necessary, fix up the partials to account for the implicit
    ! dependence of the opacity on the final P and T

    if (max_iters > 0 .AND. .NOT. skip_partials ) then

       chiRho = res(i_chiRho)
       chiT = res(i_chiT)

       dlnkap_dlnP_T = dlnkap_dlnRho/chiRho
       dlnkap_dlnT_P = dlnkap_dlnT - dlnkap_dlnRho*chiT/chiRho

       dlnP_dL_ = (dlnP_dL + dlnkap_dlnT_P*(dlnP_dlnkap*dlnT_dL - dlnP_dL*dlnT_dlnkap))/&
            (1._dp - dlnkap_dlnP_T*dlnP_dlnkap - dlnkap_dlnT_P*dlnT_dlnkap)
       dlnT_dL_ = (dlnT_dL + dlnkap_dlnP_T*(dlnT_dlnkap*dlnP_dL - dlnT_dL*dlnP_dlnkap))/&
            (1._dp - dlnkap_dlnP_T*dlnP_dlnkap - dlnkap_dlnT_P*dlnT_dlnkap)

       dlnP_dlnR_ = (dlnP_dlnR + dlnkap_dlnT_P*(dlnP_dlnkap*dlnT_dlnR - dlnP_dlnR*dlnT_dlnkap))/ &
            (1._dp - dlnkap_dlnP_T*dlnP_dlnkap - dlnkap_dlnT_P*dlnT_dlnkap)
       dlnT_dlnR_ = (dlnT_dlnR + dlnkap_dlnP_T*(dlnT_dlnkap*dlnP_dlnR - dlnT_dlnR*dlnP_dlnkap))/ &
            (1._dp - dlnkap_dlnP_T*dlnP_dlnkap - dlnkap_dlnT_P*dlnT_dlnkap)

       dlnP_dlnM_ = (dlnP_dlnM + dlnkap_dlnT_P*(dlnP_dlnkap*dlnT_dlnM - dlnP_dlnM*dlnT_dlnkap))/ &
            (1._dp - dlnkap_dlnP_T*dlnP_dlnkap - dlnkap_dlnT_P*dlnT_dlnkap)
       dlnT_dlnM_ = (dlnT_dlnM + dlnkap_dlnP_T*(dlnT_dlnkap*dlnP_dlnM - dlnT_dlnM*dlnP_dlnkap))/ &
            (1._dp - dlnkap_dlnP_T*dlnP_dlnkap - dlnkap_dlnT_P*dlnT_dlnkap)

       dlnP_dL = dlnP_dL_
       dlnT_dL = dlnT_dL_

       dlnP_dlnR = dlnP_dlnR_
       dlnT_dlnR = dlnT_dlnR_

       dlnP_dlnM = dlnP_dlnM_
       dlnT_dlnM = dlnT_dlnM_

       dlnP_dlnkap = 0._dp
       dlnT_dlnkap = 0._dp

    endif

    ! Finish

    return

  end subroutine eval_T_tau_uniform

  !****

  ! Build atmosphere structure data from T-tau relation with uniform
  ! opacity

  subroutine build_T_tau_uniform( &
       tau_surf, L, R, M, cgrav, kap, Pextra_factor, tau_outer, &
       T_tau_id, eos_proc, kap_proc, errtol, dlogtau, &
       atm_structure_num_pts, atm_structure, &
       ierr)

    use atm_def, only: atm_eos_iface, atm_kap_iface, num_results_for_build_atm
    use atm_utils, only: eval_Teff_g
    use num_lib, only: dopri5_work_sizes, dopri5

    real(dp), intent(in)      :: tau_surf
    real(dp), intent(in)      :: L
    real(dp), intent(in)      :: R
    real(dp), intent(in)      :: M
    real(dp), intent(in)      :: cgrav
    real(dp), intent(in)      :: kap
    real(dp), intent(in)      :: Pextra_factor
    real(dp), intent(in)      :: tau_outer
    integer, intent(in)       :: T_tau_id
    procedure(atm_eos_iface)  :: eos_proc
    procedure(atm_kap_iface)  :: kap_proc
    real(dp), intent(in)      :: errtol
    real(dp), intent(in)      :: dlogtau
    integer, intent(out)      :: atm_structure_num_pts
    real(dp), pointer         :: atm_structure(:,:)
    integer, intent(out)      :: ierr

    integer, parameter :: INIT_NUM_PTS = 100
    integer, parameter :: NUM_VARS = 1
    integer, parameter :: NRDENS = 0
    integer, parameter :: MAX_STEPS = 0
    integer, parameter :: LRPAR = 0
    integer, parameter :: LIPAR = 0
    integer, parameter :: IOUT = 1
    integer, parameter :: LOUT = 0

    real(dp)          :: Teff
    real(dp)          :: g
    integer           :: liwork
    integer           :: lwork
    real(dp), pointer :: work(:)
    integer, pointer  :: iwork(:)
    real(dp), target  :: rpar_ary(LRPAR)
    integer, target   :: ipar_ary(LIPAR)
    real(dp), target  :: y_ary(NUM_VARS)
    real(dp), pointer :: rpar(:)
    integer, pointer  :: ipar(:)
    real(dp), pointer :: y(:)
    real(dp)          :: lntau_surf
    real(dp)          :: lntau_outer
    real(dp)          :: rtol(NUM_VARS)
    real(dp)          :: atol(NUM_VARS)
    real(dp)          :: dlntau
    real(dp)          :: dlntau_max
    integer           :: idid

    ierr = 0

    ! Sanity check

    if (dlogtau <= 0.) then
       write(*,*) 'atm: Invalid dlogtau in build_T_tau_uniform:', dlogtau
       call mesa_error(__FILE__,__LINE__)
    end if

    ! Evaluate the effective temperature & gravity

    call eval_Teff_g(L, R, M, cgrav, Teff, g)

    ! Allocte atm_structure at its initial size

    allocate(atm_structure(num_results_for_build_atm,INIT_NUM_PTS))

    atm_structure_num_pts = 0

    if (tau_outer > tau_surf) return

    ! Allocate work rrays for the integrator

    call dopri5_work_sizes(NUM_VARS, NRDENS, liwork, lwork)
    allocate(work(lwork), iwork(liwork), stat=ierr)
    if (ierr /= 0) then
       write(*,*) 'atm: allocate failed in build_T_tau_uniform'
       deallocate(atm_structure)
       return
    end if

    work = 0._dp
    iwork = 0

    ! Set pointers (simply because dopri5 wants pointer args)

    rpar => rpar_ary
    ipar => ipar_ary

    y => y_ary

    ! Set starting values for the integrator (y(1) = delta_r)

    y(1) = 0._dp

    ! Integrate from the atmosphere base outward

    lntau_outer = log(tau_outer)
    lntau_surf = log(tau_surf)

    dlntau_max = -dlogtau*ln10
    dlntau = 0.5_dp*dlntau_max

    rtol = errtol
    atol = errtol

    call dopri5( &
         NUM_VARS, build_fcn, lntau_surf, y, lntau_outer, &
         dlntau, dlntau_max, MAX_STEPS, & 
         rtol, atol, 1, & 
         build_solout, IOUT, & 
         work, lwork, iwork, liwork, & 
         LRPAR, rpar, LIPAR, ipar, & 
         LOUT, idid)
    if (idid < 0) then
       write(*,*) 'atm: Call to dopri5 failed in build_T_tau_uniform: idid=', idid
       ierr = -1
    end if

    ! Reverse the data ordering (since by convention data in atm_structure
    ! should be ordered outward-in)

    atm_structure(:,:atm_structure_num_pts) = &
         atm_structure(:,atm_structure_num_pts:1:-1)

    ! Deallocate arrays

    deallocate(work, iwork)

  contains

    subroutine build_fcn(n, x, h, y, f, lr, rpar, li, ipar, ierr)

      use atm_def

      integer, intent(in) :: n, lr, li
      real(dp), intent(in) :: x, h
      real(dp), intent(inout) :: y(:)
      real(dp), intent(inout) :: f(:)
      integer, intent(inout), pointer :: ipar(:)
      real(dp), intent(inout), pointer :: rpar(:)
      integer, intent(out) :: ierr

      real(dp) :: atm_structure_sgl(num_results_for_build_atm)
      real(dp) :: tau
      real(dp) :: rho

      ierr = 0

      ! Evaluate structure data

      call build_data(x, y(1), atm_structure_sgl, ierr)
      if (ierr /= 0) then
         ! Non-zero ierr signals we must use a smaller stepsize;
         ! whereas in fact we want to terminate the integration.
         ! Needs fixing!
         return
      end if

      ! Set up the rhs for the optical depth equation
      ! dr/dlntau = -tau/(kappa*rho)

      tau = exp(x)
      rho = exp(atm_structure_sgl(atm_lnd))

      f(1) = -tau/(kap*rho)

      ! Finish

    end subroutine build_fcn

    !****

    subroutine build_solout( &
         nr, xold, x, n, y, rwork_y, iwork_y, interp_y, lrpar, rpar, lipar, ipar, irtrn)

      use utils_lib, only: realloc_double2

      integer, intent(in) :: nr, n, lrpar, lipar
      real(dp), intent(in) :: xold, x
      real(dp), intent(inout) :: y(:)
      real(dp), intent(inout), target :: rwork_y(*)
      integer, intent(inout), target :: iwork_y(*)
      integer, intent(inout), pointer :: ipar(:)
      real(dp), intent(inout), pointer :: rpar(:)
      interface
         real(dp) function interp_y(i, s, rwork_y, iwork_y, ierr)
           use const_def, only: dp
           integer, intent(in) :: i
           real(dp), intent(in) :: s
           real(dp), intent(inout), target :: rwork_y(*)
           integer, intent(inout), target :: iwork_y(*)
           integer, intent(out) :: ierr
         end function interp_y
      end interface
      integer, intent(out) :: irtrn

      integer :: ierr
      integer :: sz
      real(dp) :: atm_structure_sgl(num_results_for_build_atm)

      ierr = 0

      ! Evaluate structure data

      call build_data(x, y(1), atm_structure_sgl, ierr)
      if (ierr /= 0) then
         irtrn = -1
         return
      end if

      ! If necessary, expand arrays

      atm_structure_num_pts = atm_structure_num_pts + 1
      sz = size(atm_structure, dim=2)

      if (atm_structure_num_pts > sz) then
         sz = 2*sz + 100
         call realloc_double2( &
              atm_structure,num_results_for_build_atm,sz,ierr)
      end if

      ! Store data

      atm_structure(:,atm_structure_num_pts) = atm_structure_sgl

      ! Finish

      return

    end subroutine build_solout

    !***

    subroutine build_data(lntau, delta_r, atm_structure_sgl, ierr)

      use atm_def
      use atm_utils, only: eval_Paczynski_gradr
      use eos_def

      real(dp), intent(in)  :: lntau
      real(dp), intent(in)  :: delta_r
      real(dp), intent(out) :: atm_structure_sgl(:)
      integer, intent(out)  :: ierr

      real(dp) :: tau
      real(dp) :: lnT
      real(dp) :: dlnT_dL
      real(dp) :: dlnT_dlnR
      real(dp) :: dlnT_dlnM
      real(dp) :: dlnT_dlnkap
      real(dp) :: lnP
      real(dp) :: dlnP_dL
      real(dp) :: dlnP_dlnR
      real(dp) :: dlnP_dlnM
      real(dp) :: dlnP_dlnkap
      real(dp) :: lnRho
      real(dp) :: res(num_eos_basic_results)
      real(dp) :: dres_dlnRho(num_eos_basic_results)
      real(dp) :: dres_dlnT(num_eos_basic_results)
      real(dp) :: gradr

      ierr = 0

      ! Evaluate temperature and pressure at optical depth tau

      tau = exp(lntau)

      call eval_data( &
           tau, Teff, g, L, M, cgrav, &
           kap, Pextra_factor, T_tau_id, .TRUE., &
           lnT, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
           lnP, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
           ierr)
      if (ierr /= 0) then
         write(*,*) 'atm: Call to eval_data failed in build_data'
         return
      end if

      ! Calculate the density & eos results

      call eos_proc( &
           lnP, lnT, &
           lnRho, res, dres_dlnRho, dres_dlnT, &
           ierr)
      if (ierr /= 0) then
         write(*,*) 'atm: Call to eos_proc failed in build_data'
         return
      end if

      ! Evaluate radiative temperature gradient

      gradr = eval_Paczynski_gradr(exp(lnT), exp(lnP), exp(lnRho), tau, kap, L, M, R, cgrav)

      ! Store data

      atm_structure_sgl(atm_xm) = 0._dp ! We assume negligible mass in the atmosphere
      atm_structure_sgl(atm_delta_r) = delta_r
      atm_structure_sgl(atm_lnP) = lnP
      atm_structure_sgl(atm_lnd) = lnRho
      atm_structure_sgl(atm_lnT) = lnT
      atm_structure_sgl(atm_gradT) = gradr ! by assumption, atm is radiative
      atm_structure_sgl(atm_kap) = kap
      atm_structure_sgl(atm_gamma1) = res(i_gamma1)
      atm_structure_sgl(atm_grada) = res(i_grad_ad)
      atm_structure_sgl(atm_chiT) = res(i_chiT)
      atm_structure_sgl(atm_chiRho) = res(i_chiRho)
      atm_structure_sgl(atm_cp) = res(i_Cp)
      atm_structure_sgl(atm_cv) = res(i_Cv)
      atm_structure_sgl(atm_tau) = tau
      atm_structure_sgl(atm_lnfree_e) = res(i_lnfree_e)
      atm_structure_sgl(atm_dlnkap_dlnT) = 0._dp
      atm_structure_sgl(atm_dlnkap_dlnd) = 0._dp
      atm_structure_sgl(atm_lnPgas) = res(i_lnPgas)
      atm_structure_sgl(atm_gradr) = gradr

    end subroutine build_data

  end subroutine build_T_tau_uniform

  !****

  ! Evaluate atmosphere data

  subroutine eval_data( &
       tau, Teff, g, L, M, cgrav, &
       kap, Pextra_factor, T_tau_id, skip_partials, &
       lnT, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
       lnP, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
       ierr)

    use atm_t_tau_relations, only: eval_T_tau

    real(dp), intent(in)       :: tau
    real(dp), intent(in)       :: Teff
    real(dp), intent(in)       :: g
    real(dp), intent(in)       :: L
    real(dp), intent(in)       :: M
    real(dp), intent(in)       :: cgrav
    real(dp), intent(in)       :: kap
    real(dp), intent(in)       :: Pextra_factor
    integer, intent(in)        :: T_tau_id
    logical, intent(in)        :: skip_partials
    real(dp), intent(out)      :: lnT
    real(dp), intent(out)      :: dlnT_dL
    real(dp), intent(out)      :: dlnT_dlnR
    real(dp), intent(out)      :: dlnT_dlnM
    real(dp), intent(out)      :: dlnT_dlnkap
    real(dp), intent(out)      :: lnP
    real(dp), intent(out)      :: dlnP_dL
    real(dp), intent(out)      :: dlnP_dlnR
    real(dp), intent(out)      :: dlnP_dlnM
    real(dp), intent(out)      :: dlnP_dlnkap
    integer, intent(out)       :: ierr

    real(dp) :: P0
    real(dp) :: Pextra
    real(dp) :: Pfactor
    real(dp) :: P
    real(dp) :: dlogg_dlnR
    real(dp) :: dlogg_dlnM
    real(dp) :: dP0_dlnR
    real(dp) :: dP0_dlnkap
    real(dp) :: dP0_dL
    real(dp) :: dP0_dlnM
    real(dp) :: dPfactor_dlnR
    real(dp) :: dPfactor_dlnkap
    real(dp) :: dPfactor_dL
    real(dp) :: dPfactor_dlnM
    real(dp) :: dlnTeff_dL
    real(dp) :: dlnTeff_dlnR
    real(dp) :: dlnT_dlnTeff
    
    include 'formats'

    ierr = 0

    ! The analytic P(tau) relation is
    !    P = (tau*g/kap)*[1 + (kap/tau)*(L/M)/(6*pi*c*G)]
    ! The factor in square brackets comes from including nonzero Prad at tau=0
    ! see, e.g., Cox & Giuli, Section 20.1

    P0 = tau*g/kap

    Pextra = Pextra_factor*(kap/tau)*(L/M)/(6._dp*pi*clight*cgrav)

    Pfactor = 1._dp + Pextra
    P = P0*Pfactor
    lnP = log(P)
    if (is_bad(lnP)) then
       ierr = -1
       write(*,1) 'bad logP in atm_t_tau_uniform eval_data', lnP/ln10
       write(*,1) 'tau', tau
       write(*,1) 'g', g
       write(*,1) 'kap', kap
       write(*,1) 'P0', P0
       write(*,1) 'Pextra', Pextra
       write(*,1) 'P', P
       !stop 'atm_t_tau_uniform'
       return
    end if

    call eval_T_tau(T_tau_id, tau, Teff, lnT, ierr)
    if (ierr /= 0) then
       write(*,*) 'atm: Call to eval_T_tau failed in eval_data' 
       return
    end if

    ! Set up partials

    if (.NOT. skip_partials) then

       dlogg_dlnR = -2._dp
       dlogg_dlnM = 1._dp

       dP0_dL = 0._dp
       dP0_dlnR = dlogg_dlnR*P0
       dP0_dlnM = dlogg_dlnM*P0
       dP0_dlnkap = -P0

       dPfactor_dL = Pextra/L
       dPfactor_dlnR = 0._dp
       dPfactor_dlnM = -Pextra
       dPfactor_dlnkap = Pextra

       dlnP_dL = (dP0_dL*Pfactor + P0*dPfactor_dL)/P
       dlnP_dlnR = (dP0_dlnR*Pfactor + P0*dPfactor_dlnR)/P
       dlnP_dlnM = (dP0_dlnM*Pfactor + P0*dPfactor_dlnM)/P
       dlnP_dlnkap = (dP0_dlnkap*Pfactor + P0*dPfactor_dlnkap)/P

       dlnTeff_dL = 0.25_dp/L
       dlnTeff_dlnR = -0.5_dp

       dlnT_dlnTeff = 1._dp

       dlnT_dL = dlnTeff_dL*dlnT_dlnTeff
       dlnT_dlnR = dlnTeff_dlnR*dlnT_dlnTeff
       dlnT_dlnM = 0._dp
       dlnT_dlnkap = 0._dp

    else

       dlnP_dL = 0._dp
       dlnP_dlnR  = 0._dp
       dlnP_dlnM = 0._dp
       dlnP_dlnkap = 0._dp

       dlnT_dL = 0._dp
       dlnT_dlnR = 0._dp
       dlnT_dlnM = 0._dp
       dlnT_dlnkap = 0._dp

    endif

    ! Finish

    return

  end subroutine eval_data

end module atm_T_tau_uniform
