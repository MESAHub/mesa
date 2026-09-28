! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! ***********************************************************************

module overshoot_utils

  use num_lib
  use star_private_def

  implicit none

  private
  public :: eval_conv_bdy_k
  public :: eval_conv_bdy_r
  public :: eval_conv_bdy_Hp
  public :: eval_over_bdy_params

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
    end if

    if (k >= s%nz .OR. k < 1) then
       write(*,*) 'Invalid cell for convective boundary: i, k, nz=', i, k, s%nz
       ierr = -1
       return
    end if

    return

  end subroutine eval_conv_bdy_k


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

    if (is_bad_num(w) .or. w < 0d0 .or. w > 1d0) then
       write(*,'(a,i0,1x,es26.16e3)') &
          'Invalid weight for convective boundary: i, w=', i, w
       ierr = -1
       return
    end if

    associate (k_o => k, &
               k_i => k+1)

      ! Preserve a boundary that lies exactly on a mesh face.
      if (w == 0d0) then
         r = s%r(k_o)
      else if (w == 1d0) then
         r = s%r(k_i)
      else
         r = pow(w*s%r(k_i)*s%r(k_i)*s%r(k_i) + &
              (1d0-w)*s%r(k_o)*s%r(k_o)*s%r(k_o), 1d0/3d0)
      end if

    end associate

    return

  end subroutine eval_conv_bdy_r


  subroutine eval_conv_bdy_Hp (s, i, Hp, ierr)

    type(star_info), pointer :: s
    integer, intent(in)      :: i
    real(dp), intent(out)    :: Hp
    integer, intent(out)     :: ierr

    integer  :: k
    real(dp) :: r
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

    if (s%limit_overshoot_Hp_using_size_of_convection_zone) then

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
          end if

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
             end if
             r_top = s%r(s%conv_bdy_loc(i+1))
          end if

       end if

       dr = r_top - r_bot

       ! Apply the limit

       if (s%overshoot_alpha > 0d0) then
          if (s%overshoot_alpha*Hp > dr) Hp = dr/s%overshoot_alpha
       else
          if (s%alpha_mlt(k)*Hp > dr) Hp = dr/s%alpha_mlt(k)
       end if

    end if

    return

  end subroutine eval_conv_bdy_Hp


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
    real(dp) :: ri, ro

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

       end if

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

       end if

    end if

    ri = s%r(k+1)
    ro = s%r(k)

    if (is_bad_num(ri) .or. is_bad_num(ro) .or. &
        is_bad_num(r) .or. is_bad_num(r_cb)) then
       write(*,'(a,i0)') 'Nonfinite overshoot interpolation radius; k=', k
       ierr = -1
       return
    end if

    if (.not. (ri >= 0d0 .and. ro > ri .and. ri <= r .and. r <= ro)) then
       write(*,'(a,i0,3(1x,es26.16e3))') &
          'Invalid overshoot interpolation cell: k, ri, r, ro=', k, ri, r, ro
       ierr = -1
       return
    end if

    ! Interpolate the mixing length on the original mesh faces.

    w = ((ro-r)/(ro-ri))* &
        (1d0 + r/ro + pow2(r/ro))/ &
        (1d0 + ri/ro + pow2(ri/ro))
    w = min(1d0, max(0d0, w))
    lambda = (1d0-w)*s%mlt_mixing_length(k) + &
             w*s%mlt_mixing_length(k+1)

    ! Use vc = 0 at r_cb only if it brackets r with the other face.
    if (s%conv_vel(k) /= 0d0 .and. s%conv_vel(k+1) == 0d0) then
       if (ri <= r_cb .and. r_cb < ro .and. r >= r_cb) ri = r_cb
    else if (s%conv_vel(k) == 0d0 .and. s%conv_vel(k+1) /= 0d0) then
       if (ri < r_cb .and. r_cb <= ro .and. r <= r_cb) ro = r_cb
    end if

    w = ((ro-r)/(ro-ri))* &
        (1d0 + r/ro + pow2(r/ro))/ &
        (1d0 + ri/ro + pow2(ri/ro))
    w = min(1d0, max(0d0, w))
    vc = (1d0-w)*s%conv_vel(k) + w*s%conv_vel(k+1)

    ! Evaluate the diffusion coefficient

    D = vc*lambda/3._dp

    ierr = 0

    return

  end subroutine eval_over_bdy_params

end module overshoot_utils
