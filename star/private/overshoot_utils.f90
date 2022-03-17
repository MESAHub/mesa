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

module overshoot_utils

  ! Uses

  use num_lib
  use star_private_def

  ! No implicit typing

  implicit none

  ! Access specifiers

  private
  
  public :: eval_conv_bdy_k
  public :: eval_conv_bdy_r
  public :: eval_conv_bdy_Hp
  public :: eval_over_bdy_params

  ! Procedures

contains

  subroutine eval_conv_bdy_k (s, i, k, ierr)

    type(star_info), pointer :: s
    integer, intent(in)      :: i
    integer, intent(out)     :: k
    integer, intent(out)     :: ierr

    ! Evaluate the index k of the cell containing the i'th convective
    ! boundary

    ierr = 0
    
    if (s%top_conv_bdy(i)) then
       k = s%conv_bdy_loc(i)
    else
       k = s%conv_bdy_loc(i) - 1
    endif

    if (k >= s%nz .OR. k < 1) then
       write(*,*) 'Invalid cell for convective boundary: i, k, nz=', i, k, s%nz
       ierr = -1
       return
    endif

    ! Finish

    return

  end subroutine eval_conv_bdy_k

  !****

  subroutine eval_conv_bdy_r (s, i, r, ierr)

    type(star_info), pointer :: s
    integer, intent(in)      :: i
    real(dp), intent(out)    :: r
    integer, intent(out)     :: ierr

    integer  :: k
    real(dp) :: w

    ! Evaluate the radius r at the i'th convective boundary

    ! Find the convective boundary cell

    ierr = 0

    call eval_conv_bdy_k(s, i, k, ierr)
    if (ierr /= 0) return

    ! Interpolate r based on the fact that r^3 varies linearly with q
    ! across the (constant-density) cell

    w = s%cz_bdy_dq(k)/s%dq(k)

    if (w < 0._dp .OR. w > 1._dp) then
       write(*,*) 'Invalid weight for convective boundary: i, w=', i, w
       ierr = -1
       return
    end if

    associate (k_o => k, &
               k_i => k+1)

      r = pow((      w)*s%r(k_i)*s%r(k_i)*s%r(k_i) + &
                 (1._dp-w)*s%r(k_o)*s%r(k_o)*s%r(k_o), 1._dp/3._dp)

    end associate

    ! Finish

    return

  end subroutine eval_conv_bdy_r

  !****

  subroutine eval_conv_bdy_Hp (s, i, Hp, ierr)

    type(star_info), pointer :: s
    integer, intent(in)      :: i
    real(dp), intent(out)    :: Hp
    integer, intent(out)     :: ierr

    integer  :: k
    real(dp) :: r
    real(dp) :: w
    real(dp) :: x0
    real(dp) :: x1
    real(dp) :: x2
    real(dp) :: x
    real(dp) :: a0
    real(dp) :: a1
    real(dp) :: a2
    real(dp) :: P
    real(dp) :: rho
    real(dp) :: r_top
    real(dp) :: r_bot
    real(dp) :: dr

    ! Evaluate the pressure scale height Hp at the i'th convective boundary

    ! Find the convective boundary cell

    ierr = 0

    call eval_conv_bdy_k(s, i, k, ierr)
    if (ierr /= 0) return

    ! Evaluate the radius at the convective boundary

    call eval_conv_bdy_r(s, i, r, ierr)
    if (ierr /= 0) return

    ! Interpolate the pressure and density at the boundary, using a
    ! quadratic fit across the boundary cell and its neighbors (the
    ! x's are fractional mass distances from the outer edge of cell
    ! k-1); then, evaluate the pressure scale height

    associate (k_o => k-1, &
               k_m => k, &
               k_i => k+1)

      x0 = s%dq(k_o)/2._dp
      x1 = s%dq(k_o) + s%dq(k_m)/2._dp
      x2 = s%dq(k_o) + s%dq(k_m) + s%dq(k_i)/2._dp

      x = s%dq(k_o) + s%cz_bdy_dq(k)

      call two_piece_linear_coeffs(x, x0, x1, x2, a0, a1, a2, ierr)
      if (ierr /= 0) return

      P = exp(a0*s%lnPeos(k_o) + a1*s%lnPeos(k_m) + a2*s%lnPeos(k_i))
      rho = exp(a0*s%lnd(k_o) + a1*s%lnd(k_m) + a2*s%lnd(k_i))

      ! Evaluate the pressure scale height

      Hp = P/(rho*s%cgrav(k_m)* &
           (s%M_center + s%xmstar*s%conv_bdy_q(i))/(r*r))

    end associate

    ! (Possibly) limit the scale height using the size of the
    ! convection zone

    if (s% ctrl% limit_overshoot_Hp_using_size_of_convection_zone) then

       ! Determine the radial extent of the convection zone (note that
       ! r_top/r_bot don't coincide exactly with the r calculated
       ! above)

       if (s%top_conv_bdy(i)) then

          if (i == 1) then
             r_bot = s%R_center
          else
             if (s%top_conv_bdy(i-1)) then
                write(*,*) 'Double top boundary in overshoot; i=', i
                ierr = -1
                return
             end if
             r_bot = s%r(s%conv_bdy_loc(i-1))
          endif

          r_top = s%r(k)

       else

          r_bot = s%r(k+1)

          if (i == s%num_conv_boundaries) then
             r_top = s%r(1)
          else
             if (.NOT. s%top_conv_bdy(i+1)) then
                write(*,*) 'Double bottom boundary in overshoot; i=', i
                ierr = -1
                return
             endif
             r_top = s%r(s%conv_bdy_loc(i+1))
          endif

       endif
          
       dr = r_top - r_bot

       ! Apply the limit

       if (s% ctrl% overshoot_alpha > 0d0) then
          if (s% ctrl% overshoot_alpha*Hp > dr) Hp = dr/s% ctrl% overshoot_alpha
       else
          if (s%alpha_mlt(k)*Hp > dr) Hp = dr/s% ctrl% mixing_length_alpha
       end if

    end if

    ! Finish

    return

  end subroutine eval_conv_bdy_Hp

  !****

  subroutine eval_over_bdy_params (s, i, f0, k, r, D, vc, ierr)

    type(star_info), pointer :: s
    integer, intent(in)      :: i
    real(dp), intent(in)     :: f0
    integer, intent(out)     :: k
    real(dp), intent(out)    :: r
    real(dp), intent(out)    :: D
    real(dp), intent(out)    :: vc
    integer, intent(out)     :: ierr

    integer  :: k_cb
    real(dp) :: r_cb
    real(dp) :: Hp_cb
    real(dp) :: w
    real(dp) :: lambda

    ! Evaluate parameters (cell index k, radius r, diffusion
    ! coefficients D and cdc) for the overshoot boundary associated
    ! with the i'th convective boundary

    ! Find the convective boundary cell

    ierr = 0

    call eval_conv_bdy_k(s, i, k_cb, ierr)
    if (ierr /= 0) return

    ! Evaluate the radius at the convective boundary

    call eval_conv_bdy_r(s, i, r_cb, ierr)
    if (ierr /= 0) return

    ! Evaluate the pressure scale height at the convective boundary

    call eval_conv_bdy_Hp(s, i, Hp_cb, ierr)
    if (ierr /= 0) return

    ! Search for the overshoot boundary cell

    ierr = 0

    if (s%top_conv_bdy(i)) then

       ! Overshooting outward -- search inward

       r = r_cb - f0*Hp_cb

       if (r <= s%r(s%nz)) then

          r = s%r(s%nz)
          k = s%nz - 1

       else

          search_in_loop: do k = k_cb, s%nz-1
             if (s%r(k+1) <= r) exit search_in_loop
          end do search_in_loop

       endif

    else

       ! Overshooting inward -- search outward

       r = r_cb + f0*Hp_cb

       if (r >=  s%r(1)) then

          r = s%r(1)
          k = 1

       else

          search_out_loop : do k = k_cb, 1, -1
             if (s%r(k) > r) exit search_out_loop
          end do search_out_loop

       endif

    endif

    if (.NOT. (s%r(k+1) <= r .AND. s%r(k) >= r)) then
       write(*,*) 'r_ob not correctly bracketed: r(k+1), r, r(k)=', s%r(k+1), r, s%r(k)
       ierr = -1
       return
    end if

    ! Interpolate mixing parameters
 
    w = (s%r(k)*s%r(k)*s%r(k) - r*r*r)/ &
        (s%r(k)*s%r(k)*s%r(k) - s%r(k+1)*s%r(k+1)*s%r(k+1))

    lambda = (1._dp-w)*s%mlt_mixing_length(k) + w*s%mlt_mixing_length(k+1)

    if (s%conv_vel(k) /= 0._dp .AND. s%conv_vel(k+1) /= 0._dp) then

       ! Both faces of cell have non-zero mixing; interpolate vc between faces

       vc = (1._dp-w)*s%conv_vel(k) + w*s%conv_vel(k+1)

    elseif (s%conv_vel(k) /= 0._dp .AND. s%conv_vel(k+1) == 0._dp) then

       ! Outer face of cell has non-zero mixing; interpolate vc
       ! between this face and r_cb, assuming vc = 0 at the latter

        if(s%r(k) /= r_cb) then
          w = (s%r(k)*s%r(k)*s%r(k) - r*r*r)/ &
           (s%r(k)*s%r(k)*s%r(k) - r_cb*r_cb*r_cb)
        else
          w = 0d0
        end if

       vc = (1._dp-w)*s%conv_vel(k)

    elseif (s%conv_vel(k) == 0._dp .AND. s%conv_vel(k+1) /= 0._dp) then

       ! Inner face of cell has non-zero mixing; interpolate vc
       ! between this face and r_cb, assuming vc = 0 at the latter

       if(s%r(k+1) /= r_cb) then
          w = (r_cb*r_cb*r_cb - r*r*r)/ &
           (r_cb*r_cb*r_cb - s%r(k+1)*s%r(k+1)*s%r(k+1))
       else
          w = 0d0
       end if

       vc = w*s%conv_vel(k+1)

    else

       ! Neither face of cell has non-zero mixing; return

       vc = 0._dp

    endif

    ! Evaluate the diffusion coefficient

    D = vc*lambda/3._dp

    ! Finish

    ierr = 0

    return

  end subroutine eval_over_bdy_params

end module overshoot_utils
