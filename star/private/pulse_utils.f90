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
! ***********************************************************************

module pulse_utils

  ! Uses

  use star_private_def
  use const_def
  use num_lib
  use star_utils
  
  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: set_segment_indices
  public :: eval_face
  public :: eval_face_X
  public :: eval_face_A_ast
  public :: eval_face_rho
  public :: eval_center
  public :: eval_center_X
  public :: eval_center_rho
  public :: eval_center_d2

  ! Procedures

contains

  subroutine set_segment_indices (s, k_a, k_b, include_last_face)

    type(star_info), intent(in)       :: s
    integer, allocatable, intent(out) :: k_a(:)
    integer, allocatable, intent(out) :: k_b(:)
    logical, intent(in)               :: include_last_face

    logical, parameter :: DEBUG = .FALSE.

    real(dp) :: grad_mu(s%nz)
    logical  :: mask(s%nz)
    integer  :: k
    integer  :: i(s%nz)
    integer  :: n_mk
    integer  :: n_sg
    integer  :: sg

    ! Set up index ranges for segments which are delineated by
    ! composition jumps

    if (s%add_double_points_to_pulse_data) then

       ! Calculate grad_mu

       grad_mu(1) = 0d0

       do k = 2, s%nz-1
          grad_mu(k) = log(s%mu(k)/s%mu(k-1))/log(s%Peos(k)/s%Peos(k-1))
       end do

       if (include_last_face) then
          grad_mu(k) = log(s%mu(s%nz)/s%mu(s%nz-1))/log(s%Peos(s%nz)/s%Peos(s%nz-1))
       else
          grad_mu(k) = 0d0
       endif

       ! Set up the mask marking faces which will have a double point

       mask = ABS(grad_mu) > s%threshold_grad_mu_for_double_point

       if (s%max_number_of_double_points > 0) then

          ! Limit the number of marked faces

          n_mk = MIN(s%max_number_of_double_points, s%nz)

          call qsort(i, s%nz, -ABS(grad_mu))

          mask(i(n_mk+1:)) = .FALSE.

       endif

    else
       
       mask = .FALSE.

    endif

    ! Use the mask to set up the index ranges

    n_sg = COUNT(mask) + 1

    allocate(k_a(n_sg))
    allocate(k_b(n_sg))

    sg = 1

    k_a(sg) = 1

    do k = 1, s%nz

       if (mask(k)) then

          k_b(sg) = k - 1
          sg = sg + 1
          k_a(sg) = k

          if (DEBUG) then
             write(*, 100) k, grad_mu(k)
100          format('placing double point at k=', I6, 1X, '(grad_mu=', 1PE10.3, ')')
          endif

       end if

    end do

    k_b(sg) = s%nz

    ! Finish

    return

  end subroutine set_segment_indices

  !****
  
  real(dp) function eval_face (dq, v, k, k_a, k_b, v_lo, v_hi) result (v_face)

    real(dp), intent(in)           :: dq(:)
    real(dp), intent(in)           :: v(:)
    integer, intent(in)            :: k
    integer, intent(in)            :: k_a
    integer, intent(in)            :: k_b
    real(dp), intent(in), optional :: v_lo
    real(dp), intent(in), optional :: v_hi

    ! Evaluate v at face k, by interpolating (or extrapolating, if
    ! k==k_a or k==k_b+1) from cells k_a:k_b

    if (k_b < k_a) stop 'eval_face: invalid segment indices'
    if (k < k_a .OR. k > k_b+1) stop 'eval_face: out-of-bounds interpolation'

    if (k_b == k_a) then
            
       ! Using a single cell

       v_face = v(k_a)

    else

       ! Using multiple cells
       
       if (k == k_a) then
          v_face = v(k_a) - dq(k_a)*(v(k_a+1) - v(k_a))/(dq(k_a+1) + dq(k_a))
       elseif (k == k_a+1) then
          v_face = v(k_a) + dq(k_a)*(v(k_a+1) - v(k_a))/(dq(k_a+1) + dq(k_a))
       elseif (k == k_b) then
          v_face = v(k_b) - dq(k_b)*(v(k_b) - v(k_b-1))/(dq(k_b) + dq(k_b-1))
       elseif (k == k_b+1) then
          v_face = v(k_b) + dq(k_b)*(v(k_b) - v(k_b-1))/(dq(k_b) + dq(k_b-1))
       else
          v_face = interp_val_to_pt(v(k_a:k_b), k-k_a+1, k_b-k_a+1, dq(k_a:k_b), 'pulse_utils : eval_face')
       endif

    end if

    ! Apply limits

    if (PRESENT(v_lo)) then
       v_face = MAX(v_face, v_lo)
    endif

    if (PRESENT(v_hi)) then
       v_face = MIN(v_face, v_hi)
    end if

    ! Finish

    return
    
  end function eval_face

  !****
  
  real(dp) function eval_face_X (s, i, k, k_a, k_b) result (X_face)

    type(star_info), intent(in)    :: s
    integer, intent(in)            :: i
    integer, intent(in)            :: k
    integer, intent(in)            :: k_a
    integer, intent(in)            :: k_b

    ! Evaluate the abundance for species i at face k, by interpolating
    ! (or extrapolating, if k==k_a or k==k_b+1) from cells k_a:k_b

    if (k_b < k_a) stop 'eval_face_X: invalid segment indices'
    if (k < k_a .OR. k > k_b+1) stop 'eval_face_X: out-of-bounds interpolation'

    if (i >= 1) then
    
       if (k_b == k_a) then
            
          ! Using a single cell

          X_face = s%xa(i,k_a)

       else

          ! Using multiple cells
       
          if (k == k_a) then
             X_face = s%xa(i,k_a) - s%dq(k_a)*(s%xa(i,k_a+1) - s%xa(i,k_a))/(s%dq(k_a+1) + s%dq(k_a))
          elseif (k == k_a+1) then
             X_face = s%xa(i,k_a) + s%dq(k_a)*(s%xa(i,k_a+1) - s%xa(i,k_a))/(s%dq(k_a+1) + s%dq(k_a))
          elseif (k == k_b) then
             X_face = s%xa(i,k_b) - s%dq(k_b)*(s%xa(i,k_b) - s%xa(i,k_b-1))/(s%dq(k_b) + s%dq(k_b-1))
          elseif (k == k_b+1) then
             X_face = s%xa(i,k_b) + s%dq(k_b)*(s%xa(i,k_b) - s%xa(i,k_b-1))/(s%dq(k_b) + s%dq(k_b-1))
          else
             X_face = interp_val_to_pt(s%xa(i,k_a:k_b), k-k_a+1, k_b-k_a+1, s%dq(k_a:k_b), 'pulse_utils : eval_face_X')
          endif

       end if

       ! Apply limits

       X_face = MIN(1d0, MAX(0d0, X_face))

    else

       X_face = 0d0

    endif

    ! Finish

    return
    
  end function eval_face_X

  !****

  real(dp) function eval_face_A_ast (s, k, k_a, k_b) result (A_ast_face)

    type(star_info), intent(in) :: s
    integer, intent(in)         :: k
    integer, intent(in)         :: k_a
    integer, intent(in)         :: k_b

    real(dp) :: A_ast_1
    real(dp) :: A_ast_2

    ! Evaluate A*=N2*r/g (A_ast) at face k, using data from faces
    ! k_a:k_b+1. If k==k_a or k==k_b+1, then use extrapolation from
    ! neighboring faces

    if (k_b < k_a) stop 'eval_face_A_ast: invalid segment indices'
    if (k < k_a .OR. k > k_b+1) stop 'eval_face_A_ast: out-of-bounds interpolation'

    if (.not. s% calculate_Brunt_N2) stop 'eval_face_A_ast: must have calculate_Brunt_N2 = .true.'

    if (k_b == k_a) then

       A_ast_face = s%brunt_N2(k)*s%r(k)/s%grav(k)

    else

       if (k == k_a) then

          A_ast_1 = s%brunt_N2(k_a+1)*s%r(k_a+1)/s%grav(k_a+1)
          A_ast_2 = s%brunt_N2(k_a+2)*s%r(k_a+2)/s%grav(k_a+2)

          A_ast_face = A_ast_1 - s%dq(k_a)*(A_ast_2 - A_ast_1)/s%dq(k_a+1)

       elseif (k == k_b+1) then

          A_ast_1 = s%brunt_N2(k_b-1)*s%r(k_b-1)/s%grav(k_b)
          A_ast_2 = s%brunt_N2(k_b)*s%r(k_b)/s%grav(k_b)
          
          A_ast_face = A_ast_2 + s%dq(k_b)*(A_ast_2 - A_ast_1)/s%dq(k_b-1)

       else

          A_ast_face = s%brunt_N2(k)*s%r(k)/s%grav(k)

       endif

    end if
    
    ! Finish

    return

  end function eval_face_A_ast

  !****

  real(dp) function eval_face_rho (s, k, k_a, k_b) result (rho_face)


    type(star_info), intent(in) :: s
    integer, intent(in)         :: k
    integer, intent(in)         :: k_a
    integer, intent(in)         :: k_b

    real(dp) :: r
    real(dp) :: dm
    real(dp) :: dlnr

    ! Evaluate rho at face k, using data from cells k_a:k_b

    if (k_b < k_a) stop 'eval_face_rho: invalid segment indices'
    if (k < k_a .OR. k > k_b+1) stop 'eval_face_rho: out-of-bounds interpolation'

    if (k_b == k_a) then

       rho_face = s%rho(k)

    else

       if (k == k_a) then

          r = s%r(k_a)

          dm = 0.5d0*s%dm(k_a)
          dlnr = 1d0 - s%rmid(k_a)/r

       elseif (k == k_b+1) then

          r = s%r(k_b+1)

          dm = 0.5d0*s%dm(k_b)
          dlnr = s%rmid(k_b)/r - 1d0

       else

          r = s%r(k)

          dm = 0.5d0*(s%dm(k) + s%dm(k-1))
          dlnr = (s%rmid(k-1) - s%rmid(k))/r

       endif

       rho_face = dm/(pi4*r*r*r*dlnr)

    endif

    ! Finish

    return

  end function eval_face_rho

  !****

  real(dp) function eval_center (r, v, k_a, k_b, v_lo, v_hi) result (v_center)
 
    real(dp), intent(in)           :: r(:)
    real(dp), intent(in)           :: v(:)
    integer, intent(in)            :: k_a
    integer, intent(in)            :: k_b
    real(dp), intent(in), optional :: v_lo
    real(dp), intent(in), optional :: v_hi

    real(dp) :: r_1
    real(dp) :: r_2
    real(dp) :: v_1
    real(dp) :: v_2

    ! Evaluate v at the center, by extrapolating from cells/faces
    ! k_a:k_b

    if (k_b < k_a) stop 'eval_center: invalid segment indices'

    if (k_a == k_b) then

       ! Using a single cell/face

       v_center = v(k_a)

    else

       ! Using the innermost two cells/faces in k_a:k_b; fit a parabola,
       ! with dv/dr = 0 at the center

       r_1 = r(k_b)
       r_2 = r(k_b-1)

       v_1 = v(k_b)
       v_2 = v(k_b-1)

       v_center = (v_1*r_2*r_2 - v_2*r_1*r_1)/(r_2*r_2 - r_1*r_1)

    endif

    ! Apply limits

    if (PRESENT(v_lo)) then
       v_center = MAX(v_center, v_lo)
    endif

    if (PRESENT(v_hi)) then
       v_center = MIN(v_center, v_hi)
    end if

    ! Finish

    return

  end function eval_center

  !****

  real(dp) function eval_center_X (s, i, k_a, k_b) result (X_center)

    type(star_info), intent(in) :: s
    integer, intent(in)         :: i
    integer, intent(in)         :: k_a
    integer, intent(in)         :: k_b

    real(dp) :: r_1
    real(dp) :: r_2
    real(dp) :: X_1
    real(dp) :: X_2

    ! Evaluate the abundance for species i at the center, by
    ! extrapolating from cells k_a:k_b

    if (i >= 1) then

       if (k_b < k_a) stop 'eval_center: invalid segment indices'

       if (k_a == k_b) then

          ! Using a single cell

          X_center = s%xa(i,k_a)

       else

          ! Using the innermost two cells/faces in k_a:k_b; fit a parabola,
          ! with dv/dr = 0 at the center
          
          r_1 = s%rmid(k_b)
          r_2 = s%rmid(k_b-1)

          X_1 = s%xa(i,k_b)
          X_2 = s%xa(i,k_b-1)

          X_center = (X_1*r_2*r_2 - X_2*r_1*r_1)/(r_2*r_2 - r_1*r_1)

       endif

       ! Apply limits

       X_center = MIN(1d0, MAX(0d0, X_center))

    else

       X_center = 0d0

    endif

    ! Finish

    return

  end function eval_center_X

  !****

  real(dp) function eval_center_rho (s, k_b) result (rho_center)
 
    type(star_info), intent(in) :: s
    integer, intent(in)         :: k_b

    real(dp) :: r_1
    real(dp) :: rho_1
    real(dp) :: M_1

    ! Evaluate rho at the center, by extrapolating from cell k_b

    ! Fit a parabola with drho/dr = 0 at the center, which conserves mass

    r_1 = s%rmid(k_b)
    rho_1 = s%rho(k_b)
    M_1 = s%m(k_b) - 0.5d0*s%dm(k_b)

    rho_center = 3d0*(5d0*M_1/(pi*r_1*r_1*r_1) - 4d0*rho_1)/8d0

    ! Finish

    return

  end function eval_center_rho

  !****

  real(dp) function eval_center_d2 (r, v, k_a, k_b) result (d2v_center)

    real(dp), intent(in) :: r(:)
    real(dp), intent(in) :: v(:)
    integer, intent(in)  :: k_a
    integer, intent(in)  :: k_b

    real(dp) :: r_1
    real(dp) :: r_2
    real(dp) :: v_1
    real(dp) :: v_2

    ! Evaluate d2v/dr2 at the center, by extrapolating from
    ! cells/faces k_a:k_b

    if (k_b < k_a) stop 'eval_center_d2: invalid segment indices'

    if (k_a == k_b) then

       ! Using a single cell/face

       d2v_center = 0d0

    else

       ! Using the innermost two cells/faces; fit a parabola, with
       ! dv/dq = 0 at the center

       r_1 = r(k_b)
       r_2 = r(k_b-1)

       v_1 = v(k_b)
       v_2 = v(k_b-1)

       d2v_center = 2d0*(v_2 - v_1)/(r_2*r_2 - r_1*r_1)

    endif

    ! Finish

    return

  end function eval_center_d2

end module pulse_utils
