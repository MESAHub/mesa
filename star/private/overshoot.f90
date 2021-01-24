! ***********************************************************************
!
!   Copyright (C) 2010-2017 Rich Townsend, Bill Paxton
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

module overshoot

  ! Uses

  use const_def
  use num_lib
  use star_private_def

  use overshoot_utils
  use overshoot_exp
  use overshoot_step
  
  ! No implicit typing

  implicit none

  ! Access specifers

  private

  public :: add_overshooting

  ! Procedures

contains

  subroutine add_overshooting (s, ierr)

    type(star_info), pointer :: s
    integer, intent(out)     :: ierr

    logical, parameter :: DEBUG = .FALSE.

    integer  :: i
    integer  :: j
    logical  :: is_core
    logical  :: match_zone_type
    logical  :: match_zone_loc
    logical  :: match_bdy_loc
    integer  :: k_cb
    integer  :: k_a
    integer  :: k_b
    real(dp) :: D(s%nz)
    real(dp) :: vc(s%nz)
    integer  :: k
    integer  :: dk
    real(dp) :: rho
    real(dp) :: cdc
    real(dp) :: r_cb

    include 'formats'

    ! Initialize

    ierr = 0

    if (DEBUG) then
       write(*, 3) 'add_overshooting; model, n_conv_bdy=', s%model_number, s%num_conv_boundaries
    end if

    ! Loop over convective boundaries, from center to surface

    conv_bdy_loop : do i = 1, s%num_conv_boundaries

       ! Skip this boundary if it's too close to the center

       if (s%conv_bdy_q(i) < s%min_overshoot_q) then
          if (DEBUG) then
             write(*,*) 'skip since s%conv_bdy_q(i) < min_overshoot_q', i
          endif
          cycle conv_bdy_loop
       endif

       ! Skip this boundary if it's at the surface, since we don't
       ! overshoot there

       if (s%conv_bdy_loc(i) == 1) then
          if (DEBUG) then
             write(*,*) 'skip since s%conv_bdy_loc(i) == 1', i
          endif
          cycle conv_bdy_loop
       end if

       ! Loop over overshoot criteria

       criteria_loop : do j = 1, NUM_OVERSHOOT_PARAM_SETS

          if (s%overshoot_scheme(j) == '') cycle criteria_loop

          ! Check if the criteria match the current boundary

          select case (s%overshoot_zone_type(j))
          case ('burn_H')
             match_zone_type = s%burn_h_conv_region(i)
          case ('burn_He')
             match_zone_type = s%burn_he_conv_region(i)
          case ('burn_Z')
             match_zone_type = s%burn_z_conv_region(i)
          case ('nonburn')
             match_zone_type = .NOT. ( &
                  s%burn_h_conv_region(i) .OR. &
                  s%burn_he_conv_region(i) .OR. &
                  s%burn_z_conv_region(i) )              
          case ('any')
             match_zone_type = .TRUE.
          case default
             write(*,*) 'Invalid overshoot_zone_type: j, s%overshoot_zone_type(j)=', j, s%overshoot_zone_type(j)
             ierr = -1
             return
          end select

          is_core = (i == 1 .AND. s%R_center == 0d0 .AND. s%top_conv_bdy(i))

          select case (s%overshoot_zone_loc(j))
          case ('core')
             match_zone_loc = is_core
          case ('shell')
             match_zone_loc = .NOT. is_core
          case ('any')
             match_zone_loc = .TRUE.
          case default
             write(*,*) 'Invalid overshoot_zone_loc: j, s%overshoot_zone_loc(j)=', j, s%overshoot_zone_loc(j)
             ierr = -1
             return
          end select

          select case (s%overshoot_bdy_loc(j))
          case ('bottom')
             match_bdy_loc = .NOT. s%top_conv_bdy(i)
          case ('top')
             match_bdy_loc = s%top_conv_bdy(i)
          case ('any')
             match_bdy_loc = .TRUE.
          case default
             write(*,*) 'Invalid overshoot_bdy_loc: j, s%overshoot_bdy_loc(j)=', j, s%overshoot_bdy_loc(j)
             ierr = -1
             return
          end select

          if (.NOT. (match_zone_type .AND. match_zone_loc .AND. match_bdy_loc)) cycle criteria_loop

          if (DEBUG) then
             write(*,*) 'Overshooting at convective boundary: i, j=', i, j
             write(*,*) '  s%overshoot_scheme=', TRIM(s%overshoot_scheme(j))
             write(*,*) '  s%overshoot_zone_type=', TRIM(s%overshoot_zone_type(j))
             write(*,*) '  s%overshoot_zone_loc=', TRIM(s%overshoot_zone_loc(j))
             write(*,*) '  s%overshoot_bdy_loc=', TRIM(s%overshoot_bdy_loc(j))
          endif

          ! Special-case check for an overshoot scheme of 'none' (this can be used
          ! to turn *off* overshoot for specific boundary configurations)

          if (s%overshoot_scheme(j) == 'none') exit criteria_loop

          ! Evaluate convective boundary (_cb) parameters

          call eval_conv_bdy_k(s, i, k_cb, ierr)
          if (ierr /= 0) return

          ! Evaluate the overshoot diffusion coefficient and velocity
          ! using the appropriate scheme-dependent routine

          select case (s%overshoot_scheme(j))
          case ('exponential')
             call eval_overshoot_exp(s, i, j, k_a, k_b, D, vc, ierr)
          case ('step')
             call eval_overshoot_step(s, i, j, k_a, k_b, D, vc, ierr)
          case ('other')
             call s% other_overshooting_scheme(s% id, i, j, k_a, k_b, D, vc, ierr)
          case default
             write(*,*) 'Invalid overshoot_scheme:', s%overshoot_scheme(j)
             ierr = -1
          end select

          if (ierr /= 0) return

          ! Update the model

          if (s%top_conv_bdy(i)) then
             dk = -1
          else
             dk = 1
          endif

          face_loop : do k = k_a, k_b, dk

             ! Check if the overshoot will be stabilized by the stratification

             if (s%overshoot_brunt_B_max > 0._dp .and. s% calculate_Brunt_B) then
                
                if (.not. s% calculate_Brunt_N2) &
                   stop 'add_overshooting: when overshoot_brunt_B_max > 0, must have calculate_Brunt_N2 = .true.'

                ! (NB: we examine B(k+dk) rather than B(k), as the latter
                ! would allow the overshoot region to eat into a composition
                ! gradient). Special case for k == 1 or k == nz

                if (k > 1 .AND. k < s%nz) then
                   if (s%unsmoothed_brunt_B(k+dk) > s%overshoot_brunt_B_max) exit face_loop
                else
                   if (s%unsmoothed_brunt_B(k) > s%overshoot_brunt_B_max) exit face_loop
                endif

             endif

             ! Check whether D has dropped below the minimum

             if (D(k) < s%overshoot_D_min) then

                ! Update conv_bdy_dq to reflect where D drops below the minimum

                if (.not. s% conv_vel_flag) then
                   if (s%top_conv_bdy(i)) then
                      s%cz_bdy_dq(k) = find0(0._dp, D(k)-s%overshoot_D_min, s%dq(k), s%D_mix(k+1)-s%overshoot_D_min)
                      if (s%cz_bdy_dq(k) < 0._dp .OR. s%cz_bdy_dq(k) > s%dq(k)) then
                         write(*,*) 'k, k_a, k_b', k, k_a, k_b
                         write(*,*) 's%top_conv_bdy(i)=', s%top_conv_bdy(i)
                         write(*,*) 'D(k)', D(k)
                         write(*,*) 's%D_mix(k+1)', s%D_mix(k+1)
                         write(*,*) 's%overshoot_D_min', s%overshoot_D_min
                         write(*,*) 'Invalid location for overshoot boundary: cz_bdy_dq, dq=', s%cz_bdy_dq(k), s%dq(k)
                         ierr = -1
                         return
                      end if
                   else
                      s%cz_bdy_dq(k-1) = find0(0._dp, s%D_mix(k-1)-s%overshoot_D_min, s%dq(k-1), D(k)-s%overshoot_D_min)
                      if (s%cz_bdy_dq(k-1) < 0._dp .OR. s%cz_bdy_dq(k-1) > s%dq(k-1)) then
                         write(*,*) 'k, k_a, k_b', k, k_a, k_b
                         write(*,*) 's%top_conv_bdy(i)=', s%top_conv_bdy(i)
                         write(*,*) 'D(k)', D(k)
                         write(*,*) 's%D_mix(k-1)', s%D_mix(k-1)
                         write(*,*) 's%overshoot_D_min', s%overshoot_D_min
                         write(*,*) 'Invalid location for overshoot boundary: cz_bdy_dq, dq=', s%cz_bdy_dq(k-1), s%dq(k)
                         ierr = -1
                         return
                      end if
                   endif
                end if

                exit face_loop

             endif

             ! Revise mixing coefficients

             if (k > 1) then
                rho = (s%dq(k-1)*s%rho(k) + s%dq(k)*s%rho(k-1))/ &
                      (s%dq(k-1) + s%dq(k))
             else
                rho = s%rho(k)
             endif
       
             cdc = (pi4*s%r(k)*s%r(k)*rho)*(pi4*s%r(k)*s%r(k)*rho)*D(k) ! gm^2/sec

             call eval_conv_bdy_r(s, i, r_cb, ierr)
             if (ierr /= 0) then
                write(*,*) 'error calling eval_conv_bdy_r in overshoot:add_overshooting'
                return
             end if

             if (s% r(k)*dk >= r_cb*dk) then

                s%cdc(k) = cdc
                s%D_mix(k) = D(k)
                if (.not. s% conv_vel_flag) s%conv_vel(k) = vc(k)

             elseif (D(k) > s%D_mix(k)) then

                s%cdc(k) = cdc
                s%D_mix(k) = D(k)
                if (.not. s% conv_vel_flag) s%conv_vel(k) = vc(k)
                s%mixing_type(k) = overshoot_mixing

             end if

          end do face_loop

          ! If we're still in the convection zone, set the rest of the
          ! zone to be non-convective

          s%cdc(k:k_cb:dk) = 0._dp
          s%D_mix(k:k_cb:dk) = 0._dp
          if (.not. s% conv_vel_flag) s%conv_vel(k:k_cb:dk) = 0._dp
          s%mixing_type(k:k_cb:dk) = no_mixing
          
          ! Finish (we apply at most a single overshoot scheme to each boundary)

          exit criteria_loop

       end do criteria_loop

    end do conv_bdy_loop

    ! Perform a sanity check on D_mix

    check_loop : do k = 1, s%nz

       if (is_bad_num(s%D_mix(k))) then

          ierr = -1

          if (s%report_ierr) then
             write(*,3) 'mixing_type, D_mix', k, s%mixing_type(k), s%D_mix(k)
          end if

          if (s% stop_for_bad_nums) stop 'add_overshooting'

       end if

    end do check_loop

    ! Finish

    return

  end subroutine add_overshooting

end module overshoot
