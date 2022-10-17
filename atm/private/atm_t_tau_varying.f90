! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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

module atm_T_tau_varying
   
   ! Uses
   
   use const_def
   use math_lib
   use utils_lib, only : mesa_error
   
   ! No implicit typing
   
   implicit none
   
   ! Access specifiers
   
   private
   
   public :: eval_T_tau_varying
   public :: build_T_tau_varying
   
   ! Procedures

contains
   
   ! Evaluate atmosphere data from T-tau relation with varying opacity
   
   subroutine eval_T_tau_varying(&
      tau_surf, L, R, M, cgrav, &
      T_tau_id, eos_proc, kap_proc, &
      errtol, max_steps, skip_partials, &
      Teff, &
      lnT, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
      lnP, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
      ierr)
      
      use atm_def, only : atm_eos_iface, atm_kap_iface
      
      real(dp), intent(in) :: tau_surf
      real(dp), intent(in) :: L
      real(dp), intent(in) :: R
      real(dp), intent(in) :: M
      real(dp), intent(in) :: cgrav
      integer, intent(in) :: T_tau_id
      procedure(atm_eos_iface) :: eos_proc
      procedure(atm_kap_iface) :: kap_proc
      real(dp), intent(in) :: errtol
      integer, intent(in) :: max_steps
      logical, intent(in) :: skip_partials
      real(dp), intent(in) :: Teff
      real(dp), intent(out) :: lnT
      real(dp), intent(out) :: dlnT_dL
      real(dp), intent(out) :: dlnT_dlnR
      real(dp), intent(out) :: dlnT_dlnM
      real(dp), intent(out) :: dlnT_dlnkap
      real(dp), intent(out) :: lnP
      real(dp), intent(out) :: dlnP_dL
      real(dp), intent(out) :: dlnP_dlnR
      real(dp), intent(out) :: dlnP_dlnM
      real(dp), intent(out) :: dlnP_dlnkap
      integer, intent(out) :: ierr
      
      real(dp), parameter :: DLNTEFF = 1.E-4_dp
      
      real(dp) :: g
      real(dp) :: lnTeff
      real(dp) :: lnTeff_p
      real(dp) :: lnTeff_m
      real(dp) :: lnT_p
      real(dp) :: lnT_m
      real(dp) :: lnP_p
      real(dp) :: lnP_m
      integer :: ierr_p
      integer :: ierr_m
      real(dp) :: dlnTeff_dlnR
      real(dp) :: dlnTeff_dL
      real(dp) :: dlnT_dlnTeff
      real(dp) :: dlnP_dlnTeff
      
      ierr = 0
      
      ! Sanity checks
      
      if (L <= 0._dp .OR. R <= 0._dp .OR. M <= 0._dp) then
         ierr = -1
         return
      end if
      
      ! Evaluate the gravity
      
      g = cgrav * M / (R * R)
      
      ! Perform the necessary evaluations
      
      if (skip_partials) then
         
         ! No partial evaluation required
         
         call eval_data(&
            tau_surf, Teff, g, &
            T_tau_id, eos_proc, kap_proc, errtol, max_steps, &
            lnT, lnP, ierr)
         
         if (ierr /= 0) then
            write(*, *) 'atm: Call to eval_data failed in atm_t_tau_varying'
            return
         endif
         
         dlnT_dL = 0._dp
         dlnT_dlnR = 0._dp
         dlnT_dlnM = 0._dp
         dlnT_dlnkap = 0._dp
         
         dlnP_dL = 0._dp
         dlnP_dlnR = 0._dp
         dlnP_dlnM = 0._dp
         dlnP_dlnkap = 0._dp
      
      else
         
         ! Partials required, use finite differencing in Teff (we could
         ! in principle get dlnT_dlnTeff from the T-tau relation, but
         ! for consistency with dln_dlnTeff we use the same finite
         ! differencing)
         
         lnTeff = log(Teff)
         
         lnTeff_p = lnTeff + DLNTEFF
         lnTeff_m = lnTeff - DLNTEFF
         
         !$OMP SECTIONS
         
         !$OMP SECTION
         
         call eval_data(&
            tau_surf, exp(lnTeff), g, &
            T_tau_id, eos_proc, kap_proc, errtol, max_steps, &
            lnT, lnP, ierr)
         
         !$OMP SECTION
         
         call eval_data(&
            tau_surf, exp(lnTeff_p), g, &
            T_tau_id, eos_proc, kap_proc, errtol, max_steps, &
            lnT_p, lnP_p, ierr_p)
         
         !$OMP SECTION
         
         call eval_data(&
            tau_surf, exp(lnTeff_m), g, &
            T_tau_id, eos_proc, kap_proc, errtol, max_steps, &
            lnT_m, lnP_m, ierr_m)
         
         !$OMP END SECTIONS
         
         if (ierr_p /= 0) ierr = ierr_p
         if (ierr_m /= 0) ierr = ierr_m
         if (ierr /= 0) then
            write(*, *) 'Call to eval_data failed in atm_t_tau_varying_opacity'
            return
         endif
         
         ! Set up the partials
         
         dlnTeff_dlnR = -0.5_dp
         dlnTeff_dL = 0.25_dp / L
         
         dlnT_dlnTeff = (lnT_p - lnT_m) / (lnTeff_p - lnTeff_m)
         dlnT_dL = dlnT_dlnTeff * dlnTeff_dL
         dlnT_dlnR = dlnT_dlnTeff * dlnTeff_dlnR
         dlnT_dlnM = 0._dp
         dlnT_dlnkap = 0._dp
         
         dlnP_dlnTeff = (lnP_p - lnP_m) / (lnTeff_p - lnTeff_m)
         dlnP_dL = dlnP_dlnTeff * dlnTeff_dL
         dlnP_dlnR = dlnP_dlnTeff * dlnTeff_dlnR
         dlnP_dlnM = 0._dp
         dlnP_dlnkap = 0._dp
      
      endif
      
      ! Finish
      
      return
   
   end subroutine eval_T_tau_varying
   
   !****
   
   ! Evaluate atmosphere data from T-tau relation with varying opacity
   
   subroutine eval_data(&
      tau_surf, Teff, g, &
      T_tau_id, eos_proc, kap_proc, errtol, max_steps, &
      lnT, lnP, ierr)
      
      use atm_def, only : atm_eos_iface, atm_kap_iface
      
      real(dp), intent(in) :: tau_surf
      real(dp), intent(in) :: Teff
      real(dp), intent(in) :: g
      integer, intent(in) :: T_tau_id
      procedure(atm_eos_iface) :: eos_proc
      procedure(atm_kap_iface) :: kap_proc
      real(dp), intent(in) :: errtol
      integer, intent(in) :: max_steps
      real(dp), intent(out) :: lnT
      real(dp), intent(out) :: lnP
      integer, intent(out) :: ierr
      
      real(dp), parameter :: TAU_OUTER_FACTOR = 1.E-5_dp
      integer, parameter :: MAX_TRIES = 3
      real(dp), parameter :: TRY_SCALE = 10._dp
      
      real(dp) :: tau_outer_curr
      real(dp) :: errtol_curr
      integer :: try
      
      ! Try doing the integration, relaxing tau_outer and errtol if failures occur
      
      tau_outer_curr = tau_surf * TAU_OUTER_FACTOR
      errtol_curr = errtol
      
      try_loop : do try = 1, MAX_TRIES
         
         call eval_data_try(&
            tau_surf, Teff, g, tau_outer_curr, &
            T_tau_id, eos_proc, kap_proc, errtol_curr, max_steps, &
            lnT, lnP, ierr)
         if (ierr == 0) exit try_loop
         
         tau_outer_curr = tau_outer_curr * TRY_SCALE
         errtol_curr = errtol * TRY_SCALE
      
      end do try_loop
      
      ! Finish
      
      return
   
   end subroutine eval_data
   
   !****
   
   subroutine eval_data_try(&
      tau_surf, Teff, g, tau_outer, &
      T_tau_id, eos_proc, kap_proc, errtol, max_steps, &
      lnT, lnP, ierr)
      
      use eos_lib, only : radiation_pressure
      use atm_def, only : atm_eos_iface, atm_kap_iface
      use atm_T_tau_relations, only : eval_T_tau
      use num_lib, only : dopri5_work_sizes, dopri5
      
      real(dp), intent(in) :: tau_surf
      real(dp), intent(in) :: Teff
      real(dp), intent(in) :: g
      real(dp), intent(in) :: tau_outer
      integer, intent(in) :: T_tau_id
      procedure(atm_eos_iface) :: eos_proc
      procedure(atm_kap_iface) :: kap_proc
      real(dp), intent(in) :: errtol
      integer, intent(in) :: max_steps
      real(dp), intent(out) :: lnT
      real(dp), intent(out) :: lnP
      integer, intent(out) :: ierr
      
      integer, parameter :: NUM_VARS = 1
      integer, parameter :: NRDENS = 0
      integer, parameter :: LRPAR = 0
      integer, parameter :: LIPAR = 0
      integer, parameter :: IOUT = 0
      integer, parameter :: LOUT = 0
      
      real(dp), parameter :: RHO_OUTER = 1.E-10_dp
      real(dp), parameter :: DLNTAU_SURF = 1.E-3_dp
      real(dp), parameter :: DLNTAU_MAX = 0._dp
      
      integer :: liwork
      integer :: lwork
      real(dp), pointer :: work(:)
      integer, pointer :: iwork(:)
      real(dp), target :: rpar_ary(LRPAR)
      integer, target :: ipar_ary(LIPAR)
      real(dp), target :: y_ary(NUM_VARS)
      real(dp), pointer :: rpar(:)
      integer, pointer :: ipar(:)
      real(dp), pointer :: y(:)
      real(dp) :: lnTeff
      real(dp) :: T_outer
      real(dp) :: Pgas_outer
      real(dp) :: P_outer
      real(dp) :: lntau_outer
      real(dp) :: lntau_surf
      real(dp) :: dlntau
      real(dp) :: rtol(NUM_VARS)
      real(dp) :: atol(NUM_VARS)
      integer :: idid
      
      ierr = 0
      
      ! Allocate work arrays for the integrator
      
      call dopri5_work_sizes(NUM_VARS, NRDENS, liwork, lwork)
      allocate(work(lwork), iwork(liwork), stat = ierr)
      if (ierr /= 0) then
         write(*, *) 'allocate failed in eval_data_try'
         return
      end if
      
      work = 0._dp
      iwork = 0
      
      ! Set pointers (simply because dopri5 wants pointer args)
      
      rpar => rpar_ary
      ipar => ipar_ary
      
      y => y_ary
      
      ! Set values at the outer boundary (y(1) = lnP; estimate
      ! Pgas from low-density ideal gas law)
      
      lnTeff = log(Teff)
      
      call eval_T_tau(T_tau_id, tau_outer, Teff, lnT, ierr)
      if (ierr /= 0) then
         write(*, *) 'atm: Call to eval_T_tau failed in eval_data_try'
         return
      end if
      
      T_outer = exp(lnT)
      Pgas_outer = cgas * RHO_OUTER * T_outer
      P_outer = Pgas_outer + radiation_pressure(T_outer)
      lnP = log(P_outer)
      
      y(1) = lnP
      
      ! Integrate inward from tau_outer to tau_surf
      
      lntau_outer = log(tau_outer)
      lntau_surf = log(tau_surf)
      
      dlntau = DLNTAU_SURF
      
      rtol = errtol
      atol = errtol
      
      call dopri5(&
         NUM_VARS, eval_fcn, lntau_outer, y, lntau_surf, &
         dlntau, DLNTAU_MAX, max_steps, &
         rtol, atol, 1, &
         eval_solout, IOUT, &
         work, lwork, iwork, liwork, &
         LRPAR, rpar, LIPAR, ipar, &
         LOUT, idid)
      if (idid < 0) then
         write(*, *) 'Call to dopri5 failed in eval_data_try: idid=', idid
         ierr = -1
         return
      end if
      
      ! Store the final pressure and temperature
      
      lnP = y(1)
      
      call eval_T_tau(T_tau_id, tau_surf, Teff, lnT, ierr)
      if (ierr /= 0) then
         write(*, *) 'atm: Call to eval_T_tau failed in eval_data_try'
         return
      end if
      
      ! Deallocate arrays
      
      deallocate(work, iwork)
      
      ! Finish
      
      return
   
   contains
      
      subroutine eval_fcn(n, x, h, y, f, lr, rpar, li, ipar, ierr)
         
         use eos_def, only : num_eos_basic_results
         use atm_T_tau_relations, only : eval_T_tau
         
         integer, intent(in) :: n, lr, li
         real(dp), intent(in) :: x, h
         real(dp), intent(inout) :: y(:)
         real(dp), intent(inout) :: f(:)
         integer, intent(inout), pointer :: ipar(:)
         real(dp), intent(inout), pointer :: rpar(:)
         integer, intent(out) :: ierr
         
         real(dp) :: tau
         real(dp) :: lnP
         real(dp) :: lnT
         real(dp) :: lnRho
         real(dp) :: res(num_eos_basic_results)
         real(dp) :: dres_dlnRho(num_eos_basic_results)
         real(dp) :: dres_dlnT(num_eos_basic_results)
         real(dp) :: kap
         real(dp) :: dlnkap_dlnRho
         real(dp) :: dlnkap_dlnT
         
         ierr = 0
         
         ! Calculate the temperature
         
         tau = exp(x)
         
         lnP = y(1)
         
         call eval_T_tau(T_tau_id, tau, Teff, lnT, ierr)
         if (ierr /= 0) then
            write(*, *) 'atm: Call to eval_T_tau failed in eval_fcn'
            return
         end if
         
         ! Calculate the density and eos results
         
         call eos_proc(&
            lnP, lnT, &
            lnRho, res, dres_dlnRho, dres_dlnT, &
            ierr)
         if (ierr /= 0) then
            write(*, *) 'atm: Call to eos_proc failed in eval_fcn'
            return
         end if
         
         ! Calculate the opacity
         
         call kap_proc(&
            lnRho, lnT, res, dres_dlnRho, dres_dlnT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, &
            ierr)
         if (ierr /= 0) then
            write(*, *) 'atm: Call to kap_proc failed in eval_fcn'
            return
         end if
         
         ! Set up the rhs for the hydrostatic eqm equation
         ! dlnP/dlntau = tau*g/P*kappa
         
         f(1) = tau * g / (kap * exp(lnP))
         
         ! Finish
         
         return
      
      end subroutine eval_fcn
      
      !****
      
      subroutine eval_solout(&
         nr, xold, x, n, y, rwork_y, iwork_y, interp_y, lrpar, rpar, lipar, ipar, irtrn)
         
         integer, intent(in) :: nr, n, lrpar, lipar
         real(dp), intent(in) :: xold, x
         real(dp), intent(inout) :: y(:)
         real(dp), intent(inout), target :: rwork_y(*)
         integer, intent(inout), target :: iwork_y(*)
         integer, intent(inout), pointer :: ipar(:)
         real(dp), intent(inout), pointer :: rpar(:)
         interface
            real(dp) function interp_y(i, s, rwork_y, iwork_y, ierr)
               use const_def, only : dp
               integer, intent(in) :: i
               real(dp), intent(in) :: s
               real(dp), intent(inout), target :: rwork_y(*)
               integer, intent(inout), target :: iwork_y(*)
               integer, intent(out) :: ierr
            end function interp_y
         end interface
         integer, intent(out) :: irtrn ! < 0 causes solver to return to calling program.
         
         irtrn = 0 ! for ifort
         
         ! Dummy routine that's never called
         
         call mesa_error(__FILE__, __LINE__, 'Bogus call to eval_solout')
      
      end subroutine eval_solout
   
   end subroutine eval_data_try
   
   !****
   
   ! Build atmosphere structure data from T-tau relation with varying
   ! opacity
   
   subroutine build_T_tau_varying(&
      tau_surf, L, R, Teff, M, cgrav, lnP_surf, tau_outer, &
      T_tau_id, eos_proc, kap_proc, errtol, dlogtau, &
      atm_structure_num_pts, atm_structure, &
      ierr)
      
      use atm_def
      use atm_utils, only : eval_Teff_g
      use num_lib, only : dopri5_work_sizes, dopri5
      
      real(dp), intent(in) :: tau_surf
      real(dp), intent(in) :: L
      real(dp), intent(in) :: R
      real(dp), intent(in) :: Teff
      real(dp), intent(in) :: M
      real(dp), intent(in) :: cgrav
      real(dp), intent(in) :: lnP_surf
      real(dp), intent(in) :: tau_outer
      integer, intent(in) :: T_tau_id
      procedure(atm_eos_iface) :: eos_proc
      procedure(atm_kap_iface) :: kap_proc
      real(dp), intent(in) :: errtol
      real(dp), intent(in) :: dlogtau
      integer, intent(out) :: atm_structure_num_pts
      real(dp), pointer :: atm_structure(:, :)
      integer, intent(out) :: ierr
      
      integer, parameter :: INIT_NUM_PTS = 100
      integer, parameter :: NUM_VARS = 2
      integer, parameter :: NRDENS = 0
      integer, parameter :: MAX_STEPS = 0
      integer, parameter :: LRPAR = 0
      integer, parameter :: LIPAR = 0
      integer, parameter :: IOUT = 1
      integer, parameter :: LOUT = 0
      
      real(dp) :: g
      integer :: liwork
      integer :: lwork
      real(dp), pointer :: work(:)
      integer, pointer :: iwork(:)
      real(dp), target :: rpar_ary(LRPAR)
      integer, target :: ipar_ary(LIPAR)
      real(dp), target :: y_ary(NUM_VARS)
      real(dp), pointer :: rpar(:)
      integer, pointer :: ipar(:)
      real(dp), pointer :: y(:)
      real(dp) :: lntau_surf
      real(dp) :: lntau_outer
      real(dp) :: rtol(NUM_VARS)
      real(dp) :: atol(NUM_VARS)
      real(dp) :: dlntau
      real(dp) :: dlntau_max
      integer :: idid
      
      ierr = 0
      
      ! Sanity checks
      
      if (L <= 0._dp .OR. R <= 0._dp .OR. M <= 0._dp) then
         ierr = -1
         return
      end if
      
      if (dlntau <= 0.) then
         write(*, *) 'Invalid dlntau in build_atm_uniform:', dlntau
         call mesa_error(__FILE__, __LINE__)
      end if
      
      ! Evaluate the gravity
      
      g = cgrav * M / (R * R)
      
      ! Allocte atm_structure at its initial size
      
      allocate(atm_structure(num_results_for_build_atm, INIT_NUM_PTS))
      
      atm_structure_num_pts = 0
      
      if (tau_outer > tau_surf) return
      
      ! Allocate work rrays for the integrator
      
      call dopri5_work_sizes(NUM_VARS, NRDENS, liwork, lwork)
      allocate(work(lwork), iwork(liwork), stat = ierr)
      if (ierr /= 0) then
         write(*, *) 'atm: allocate failed in build_T_tau_varying'
         deallocate(atm_structure)
         return
      end if
      
      work = 0._dp
      iwork = 0
      
      ! Set pointers (simply because dopri5 wants pointer args)
      
      rpar => rpar_ary
      ipar => ipar_ary
      
      y => y_ary
      
      ! Set starting values for the integrator (y(1) = delta_r; y(2) = lnP)
      
      y(1) = 0._dp
      y(2) = lnP_surf
      
      ! Integrate from the atmosphere base outward
      
      lntau_outer = log(tau_outer)
      lntau_surf = log(tau_surf)
      
      dlntau_max = -dlogtau * ln10
      dlntau = 0.5_dp * dlntau_max
      
      rtol = errtol
      atol = errtol
      
      call dopri5(&
         NUM_VARS, build_fcn, lntau_surf, y, lntau_outer, &
         dlntau, dlntau_max, MAX_STEPS, &
         rtol, atol, 1, &
         build_solout, IOUT, &
         work, lwork, iwork, liwork, &
         LRPAR, rpar, LIPAR, ipar, &
         LOUT, idid)
      if (idid < 0) then
         write(*, *) 'atm: Call to dopri5 failed in build_T_tau_varying: idid=', idid
         ierr = -1
      end if
      
      ! Reverse the data ordering (since by convention data in atm_structure
      ! should be ordered outward-in)
      
      atm_structure(:, :atm_structure_num_pts) = &
         atm_structure(:, atm_structure_num_pts:1:-1)
      
      ! Deallocate arrays
      
      deallocate(work, iwork)
   
   contains
      
      subroutine build_fcn(n, x, h, y, f, lr, rpar, li, ipar, ierr)
         
         integer, intent(in) :: n, lr, li
         real(dp), intent(in) :: x, h
         real(dp), intent(inout) :: y(:)
         real(dp), intent(inout) :: f(:)
         integer, intent(inout), pointer :: ipar(:)
         real(dp), intent(inout), pointer :: rpar(:)
         integer, intent(out) :: ierr
         
         real(dp) :: tau
         real(dp) :: P
         real(dp) :: rho
         real(dp) :: kap
         
         real(dp) :: atm_structure_sgl(num_results_for_build_atm)
         
         ierr = 0
         
         ! Evaluate structure data
         
         call build_data(x, y(1), y(2), atm_structure_sgl, ierr)
         if (ierr /= 0) then
            ! Non-zero ierr signals we must use a smaller stepsize;
            ! whereas in fact we want to terminate the integration.
            ! Needs fixing!
            return
         end if
         
         ! Set up the rhs for the optical depth and hydrostatic
         ! equilibrium equations
         ! dr/dlntau = -tau/(kappa*rho)
         !
         
         tau = exp(x)
         P = exp(y(2))
         rho = exp(atm_structure_sgl(atm_lnd))
         kap = atm_structure_sgl(atm_kap)
         
         f(1) = -tau / (kap * rho)
         f(2) = tau * g / (kap * P)
         
         ! Finish
      
      end subroutine build_fcn
      
      !****
      
      subroutine build_solout(&
         nr, xold, x, n, y, rwork_y, iwork_y, interp_y, lrpar, rpar, lipar, ipar, irtrn)
         
         use utils_lib, only : realloc_double2
         
         integer, intent(in) :: nr, n, lrpar, lipar
         real(dp), intent(in) :: xold, x
         real(dp), intent(inout) :: y(:)
         real(dp), intent(inout), target :: rwork_y(*)
         integer, intent(inout), target :: iwork_y(*)
         integer, intent(inout), pointer :: ipar(:)
         real(dp), intent(inout), pointer :: rpar(:)
         interface
            real(dp) function interp_y(i, s, rwork_y, iwork_y, ierr)
               use const_def, only : dp
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
         
         call build_data(x, y(1), y(2), atm_structure_sgl, ierr)
         if (ierr /= 0) then
            irtrn = -1
            return
         end if
         
         ! If necessary, expand arrays
         
         atm_structure_num_pts = atm_structure_num_pts + 1
         sz = size(atm_structure, dim = 2)
         
         if (atm_structure_num_pts > sz) then
            sz = 2 * sz + 100
            call realloc_double2(&
               atm_structure, num_results_for_build_atm, sz, ierr)
         end if
         
         ! Store data
         
         atm_structure(:, atm_structure_num_pts) = atm_structure_sgl
         
         ! Finish
         
         return
      
      end subroutine build_solout
      
      !***
      
      subroutine build_data(lntau, delta_r, lnP, atm_structure_sgl, ierr)
         
         use eos_def
         use eos_lib, only : radiation_pressure
         use atm_T_tau_relations, only : eval_T_tau
         use atm_utils, only : eval_Paczynski_gradr
         
         real(dp), intent(in) :: lntau
         real(dp), intent(in) :: delta_r
         real(dp), intent(in) :: lnP
         real(dp), intent(out) :: atm_structure_sgl(:)
         integer, intent(out) :: ierr
         
         real(dp) :: tau
         real(dp) :: lnT
         real(dp) :: lnRho
         real(dp) :: res(num_eos_basic_results)
         real(dp) :: dres_dlnRho(num_eos_basic_results)
         real(dp) :: dres_dlnT(num_eos_basic_results)
         real(dp) :: kap
         real(dp) :: dlnkap_dlnRho
         real(dp) :: dlnkap_dlnT
         real(dp) :: gradr
         
         ierr = 0
         
         ! Calculate the temperature
         
         tau = exp(lntau)
         
         call eval_T_tau(T_tau_id, tau, Teff, lnT, ierr)
         if (ierr /= 0) then
            write(*, *) 'atm: Call to eval_T_tau failed in build_data'
            return
         end if
         
         ! Calculate the density and eos results
         
         call eos_proc(&
            lnP, lnT, &
            lnRho, res, dres_dlnRho, dres_dlnT, &
            ierr)
         if (ierr /= 0) then
            write(*, *) 'atm: Call to eos_proc failed in build_data'
            return
         end if
         
         ! Calculate the opacity
         
         call kap_proc(&
            lnRho, lnT, res, dres_dlnRho, dres_dlnT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, &
            ierr)
         if (ierr /= 0) then
            write(*, *) 'atm: Call to kap_proc failed in build_data'
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
         atm_structure_sgl(atm_dlnkap_dlnT) = dlnkap_dlnT
         atm_structure_sgl(atm_dlnkap_dlnd) = dlnkap_dlnRho
         atm_structure_sgl(atm_lnPgas) = res(i_lnPgas)
         atm_structure_sgl(atm_gradr) = gradr
      
      end subroutine build_data
   
   end subroutine build_T_tau_varying

end module atm_T_tau_varying
