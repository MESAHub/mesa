! ***********************************************************************
!
!   Copyright (C) 2010-2024  The MESA Team
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
! ***********************************************************************

module parasite_model

   use const_def, only: dp, pi
   use num_lib
   use math_lib
   use fingering_modes
   use parasite_model_matrices

   use f95_lapack

   implicit none

   real(dp), parameter :: C2 = 0.33_dp

   private
   public :: eval_parasite_saturation

   !real(dp), parameter :: CH = 1.66_dp

contains

   subroutine eval_parasite_saturation(Pr, tau, R_0, H_B, D_B, &
      lam_hat, l2_hat, k_z, N, safety, sigma_max, k_z_max, w, ierr)

      ! Evaluate the parasitic-mode saturation defined by eqn.(28) of Fraser, Reifenstein, & Garaud (2024, FRG24).
      ! This is implemented as a root finding problem with w the search variable. safety controls the level of
      ! performance-enhancing optimizations/approximations:
      !
      ! safety = 0: Low-tau limit, low-Rm limit, ditch even-parity block a (which has the ordinary KH mode), use HG19 away from large R_0
      ! safety = 1: low-tau limit, uses low-Rm limit for odd-parity block, still calculates even-parity block
      ! (I wonder if there's a middle ground here that's basically safety=2 but we use HG19 away from large R_0?
      !  would require knowing a priori what constitutes "large R_0")
      ! safety = 2: low-tau limit
      ! safety = 3 (max): use full model (I don't think this option will ever be needed since tau is always tiny)
      ! safety = 4 (actual max for now while debugging): use the original implementation; gives same answers as safety=3 but needlessly slower

      real(dp), intent(in)  :: Pr
      real(dp), intent(in)  :: tau
      real(dp), intent(in)  :: R_0
      real(dp), intent(in)  :: H_B
      real(dp), intent(in)  :: D_B
      real(dp), intent(in)  :: lam_hat
      real(dp), intent(in)  :: l2_hat
      real(dp), intent(in)  :: k_z(:)
      integer, intent(in)   :: N ! +/- bounds of mode expansion in eqn. (40) of FRG24
      integer, intent(in)   :: safety
      real(dp), intent(out) :: sigma_max
      real(dp), intent(out) :: k_z_max
      real(dp), intent(out) :: w
      integer, intent(out)  :: ierr

      integer, parameter  :: IMAX = 50
      real(dp), parameter :: EPSX = 2e-12_dp ! may want to put more thought into optimal x and y tolerances
      real(dp), parameter :: EPSY = 0._dp

      real(dp) :: l_hat
      real(dp) :: w_1
      real(dp) :: w_2
      real(dp) :: f_1
      real(dp) :: f_2
      real(dp) :: dfdx
      integer  :: i

      real(dp), pointer :: rpar(:) => null() ! not used, but needed to pass to safe_root_with_brackets
      integer, pointer  :: ipar(:) => null() ! not used, but needed to pass to safe_root_with_brackets
      integer           :: lrpar = 0         ! not used, but needed to pass to safe_root_with_brackets
      integer           :: lipar = 0         ! not used, but needed to pass to safe_root_with_brackets

      ierr = 0

      l_hat = sqrt(l2_hat)

      ! Set initial bracketing interval

      w_1 = 2.0_dp*pi*lam_hat/l_hat ! this is the BGS13 solution
      w_2 = w_1 + sqrt(2.0_dp*H_B)  ! this is the HG19 solution for large HB. I wonder if we should use the actual HG19 solution

      f_1 = root_func_(w_1, dfdx, lrpar, rpar, lipar, ipar, ierr)
      if (ierr /= 0) return

      f_2 = root_func_(w_2, dfdx, lrpar, rpar, lipar, ipar, ierr)
      if (ierr /= 0) return

      ! Check whether f_1 < 0, and if necessary adjust until w_1 is
      ! valid lower bound

      i = 0

      do while (f_1 > 0._dp)

         write(*,*) 'f_1 > 0, resetting brackets:', w_1, '-->', w_1/10d0

         if (i > 20) then ! we tried decreasing w_1 by 20 orders of magnitude, so need to break
            write(*, *) 'Can''t find valid lower bracket for FRG w search'
            ierr = -1
            return
         end if

         i = i + 1

         ! Set new lower bound

         w_2 = w_1
         f_2 = f_1

         w_1 = w_1/10._dp
         f_1 = root_func_(w_1, dfdx, lrpar, rpar, lipar, ipar, ierr)
         if (ierr /= 0) return

      end do

      ! Check whether f_2 > 0, and if necessary adjust until w_2 is
      ! valid upper bound

      i = 0

      do while (f_2 < 0._dp)

         !write(*,*) 'f_2 < 0, resetting brackets:', w_2, '-->', 10*w_2

         if (i > 20) then ! we tried increasing w_2 by 20 orders of magnitude, so need to break
            write(*, *) 'Can''t find valid upper bracket for FRG w search'
            ierr = -1
         end if

         i = i + 1

         ! Set new lower bound

         w_1 = w_2
         f_1 = f_2

         w_2 = w_2*10._dp
         f_2 = root_func_(w_2, dfdx, lrpar, rpar, lipar, ipar, ierr)
         if (ierr /= 0) return
         
      end do

      ! Call bracketed root finder

      w = safe_root_with_brackets(root_func_, &
         w_1, w_2, f_1, f_2, imax, epsx, epsy, lrpar, rpar, lipar, ipar, ierr)
      if (ierr /= 0) return

      ! Set sigma_max and k_max

      call find_fastest_parasite(w, k_z, Pr, tau, R_0, H_B, D_B, lam_hat, l_hat, N, safety, sigma_max, k_z_max, ierr)
      
   contains

      function root_func_(x, dfdx, lrpar, rpar, lipar, ipar, ierr) result(f)

         real(dp), intent(in)             :: x
         real(dp), intent(out)            :: dfdx
         integer, intent(in)              :: lrpar
         real(dp), intent(inout), pointer :: rpar(:)
         integer, intent(in)              :: lipar
         integer, intent(inout), pointer  :: ipar(:)
         integer, intent(out)             :: ierr
         real(dp)                         :: f

         real(dp) :: sigma_max, k_z_max

         ! Evaluate the root func, defined as rhs - lhs of eqn. (28) of FRG24

         associate (w => x)
            call find_fastest_parasite(w, k_z, Pr, tau, R_0, H_B, D_B, lam_hat, l_hat, N, safety, sigma_max, k_z_max, ierr)
         end associate

         f = sigma_max*C2 - lam_hat

         dfdx = 0._dp

      end function root_func_

   end subroutine eval_parasite_saturation
         
   !****

   subroutine find_fastest_parasite(w, k_z, Pr, tau, R_0, H_B, D_B, lam_hat, l_hat, N, safety, sigma_max, k_z_max, ierr)

      real(dp), intent(in)  :: w
      real(dp), intent(in)  :: k_z(:)
      real(dp), intent(in)  :: Pr
      real(dp), intent(in)  :: tau
      real(dp), intent(in)  :: R_0
      real(dp), intent(in)  :: H_B
      real(dp), intent(in)  :: D_B
      real(dp), intent(in)  :: lam_hat
      real(dp), intent(in)  :: l_hat
      integer, intent(in)   :: N
      integer, intent(in)   :: safety
      real(dp), intent(out) :: sigma_max
      real(dp), intent(out) :: k_z_max
      integer, intent(out)  :: ierr

      real(dp) :: sigma_max_i
      integer  :: i
      integer  :: i_max
      integer  :: n_k_z

      ! Over the range of wavenumbers k_z(:), find the fastest-growing parasite mode

      ierr = 0

      sigma_max = -HUGE(0._dp)
      k_z_max = -HUGE(0._dp)

      n_k_z = SIZE(k_z)

      k_loop : do i = 1, n_k_z

         call find_fastest_parasite_k_z(w, k_z(i), Pr, tau, R_0, H_B, D_B, lam_hat, l_hat, N, safety, sigma_max_i, ierr)
         if (ierr /= 0) return

         if (sigma_max_i > sigma_max) then 
            sigma_max = sigma_max_i
            k_z_max = k_z(i)
            i_max = i
         end if

      end do k_loop

      ! Check for marginal cases (commented out for now)

      ! if (sigma_max > 0._dp .AND. &
      !     (i_max == 1 .OR. i_max == n_k_z)) then
      !    write(*,*) 'warning: most unstable growth at edge of k range:', k_z(1), k_z(n_k_z), k_z(i_max)
      !    ierr = -1
      ! end if

   end subroutine find_fastest_parasite

   !****

   subroutine find_fastest_parasite_k_z(w, k_z, Pr, tau, R_0, H_B, D_B, lam_hat, l_hat, N, safety, sigma_max, ierr)

      real(dp), intent(in)  :: w
      real(dp), intent(in)  :: k_z
      real(dp), intent(in)  :: Pr
      real(dp), intent(in)  :: tau
      real(dp), intent(in)  :: R_0
      real(dp), intent(in)  :: H_B
      real(dp), intent(in)  :: D_B
      real(dp), intent(in)  :: lam_hat
      real(dp), intent(in)  :: l_hat
      integer, intent(in)   :: N
      integer, intent(in)   :: safety
      real(dp), intent(out) :: sigma_max
      integer, intent(out)  :: ierr

      real(dp) :: sigma_max_odd
      real(dp) :: sigma_max_even

      ! For the given flow speed w and vertical wavenumber k_z, find
      ! the fastest-growing parasite mode

      ierr = 0

      select case(safety)
      case (0)
         call eval_max_eigval_(build_parasite_matrix_LPN_QS, sigma_max_odd, parity='ODD')
         sigma_max = sigma_max_odd
      case (1)
         call eval_max_eigval_(build_parasite_matrix_LPN, sigma_max_even, parity='EVEN')
         call eval_max_eigval_(build_parasite_matrix_LPN_QS, sigma_max_odd, parity='ODD')
         sigma_max = MAX(sigma_max_even, sigma_max_odd)
      case (2)
         call eval_max_eigval_(build_parasite_matrix_LPN, sigma_max_even, parity='EVEN')
         call eval_max_eigval_(build_parasite_matrix_LPN, sigma_max_odd, parity='ODD')
         sigma_max = MAX(sigma_max_even, sigma_max_odd)
      case (3)
         call eval_max_eigval_(build_parasite_matrix, sigma_max_even, parity='EVEN')
         call eval_max_eigval_(build_parasite_matrix, sigma_max_odd, parity='ODD')
         sigma_max = MAX(sigma_max_even, sigma_max_odd)
      case (4)
         call eval_max_eigval_(build_parasite_matrix, sigma_max)
      case default
         write(*, *) '** invalid safety option in find_fastest_parasite_k'
         ierr = -1
      end select

   contains

      subroutine eval_max_eigval_(build_matrix, sigma_max, parity)

         interface
            subroutine build_matrix(w, k_z, Pr, tau, R_0, H_B, D_B, lam_hat, l_hat, N, L, parity)
               use const_def, only: dp
               real(dp), intent(in)               :: w
               real(dp), intent(in)               :: k_z
               real(dp), intent(in)               :: Pr
               real(dp), intent(in)               :: tau
               real(dp), intent(in)               :: R_0
               real(dp), intent(in)               :: H_B
               real(dp), intent(in)               :: D_B
               real(dp), intent(in)               :: lam_hat
               real(dp), intent(in)               :: l_hat
               integer, intent(in)                :: N
               real(dp), allocatable, intent(out) :: L(:,:)
               character(*), intent(in), optional :: parity
            end subroutine build_matrix
         end interface
         real(dp), intent(out)              :: sigma_max
         character(*), intent(in), optional :: parity

         real(dp), allocatable :: L(:,:)
         integer               :: s
         real(dp), allocatable :: sigma_r(:)
         real(dp), allocatable :: sigma_i(:)

         ! Build the matrix

         call build_matrix(w, k_z, Pr, tau, R_0, H_B, D_B, lam_hat, l_hat, N, L, parity)

         ! Reorder elements to improve numerical stability

         s = SIZE(L, 1)

         L = L(s:1:-1,s:1:-1)

         ! Calculate its eigenvalues

         allocate(sigma_r(s))
         allocate(sigma_i(s))

         call LA_GEEV(L, sigma_r, sigma_i)

         ! Extract the maximal real part

         sigma_max = MAXVAL(sigma_r)

      end subroutine eval_max_eigval_

   end subroutine find_fastest_parasite_k_z

end module parasite_model
