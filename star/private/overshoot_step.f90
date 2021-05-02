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
! ***********************************************************************

module overshoot_step

  ! Uses

  use star_private_def

  use overshoot_utils

  ! No implicit typing

  implicit none

  ! Access specifers

  private

  public :: eval_overshoot_step

  ! Procedures

contains

  subroutine eval_overshoot_step (s, i, j, k_a, k_b, D, vc, ierr)

    type(star_info), pointer :: s
    integer, intent(in)      :: i
    integer, intent(in)      :: j
    integer, intent(out)     :: k_a
    integer, intent(out)     :: k_b
    real(dp), intent(out)    :: D(:)
    real(dp), intent(out)    :: vc(:)
    integer, intent(out)     :: ierr

    real(dp) :: f
    real(dp) :: f0
    real(dp) :: D0
    real(dp) :: Delta0
    real(dp) :: w
    real(dp) :: factor
    real(dp) :: Hp_cb
    integer  :: k_ob
    real(dp) :: r_ob
    real(dp) :: D_ob
    real(dp) :: vc_ob
    logical  :: outward
    integer  :: dk
    integer  :: k
    real(dp) :: r
    real(dp) :: dr

    ! Evaluate the overshoot diffusion coefficient D(k_a:k_b) and
    ! mixing velocity vc(k_a:k_b) at the i'th convective boundary,
    ! using the j'th set of overshoot parameters. The overshoot
    ! follows a simple step scheme

    ierr = 0

    ! Extract parameters

    f = s%overshoot_f(j)
    f0 = s%overshoot_f0(j)

    D0 = s%overshoot_D0(j)
    Delta0 = s%overshoot_Delta0(j)

    if (f <= 0._dp .OR. f0 <= 0._dp) then
       write(*,*) 'ERROR: for step overshooting, must set f and f0 > 0'
       write(*,*) 'see description of overshooting in star/defaults/control.defaults'
       ierr = -1
       return
    end if

    ! Apply mass limits

    if (s%star_mass < s%overshoot_mass_full_on(j)) then
       if (s%star_mass > s%overshoot_mass_full_off(j)) then
          w = (s%star_mass - s%overshoot_mass_full_off(j)) / &
              (s%overshoot_mass_full_on(j) - s%overshoot_mass_full_off(j))
          factor = 0.5_dp*(1._dp - cospi(w))
          f = f*factor
          f0 = f0*factor
       else
          f = 0._dp
          f0 = 0._dp
       endif
    endif

    ! Evaluate convective boundary (_cb) parameters

    call eval_conv_bdy_Hp(s, i, Hp_cb, ierr)
    if (ierr /= 0) return

    ! Evaluate overshoot boundary (_ob) parameters

    call eval_over_bdy_params(s, i, f0, k_ob, r_ob, D_ob, vc_ob, ierr)
    if (ierr /= 0) return

    ! Loop over cell faces, adding overshoot until D <= overshoot_D_min

    outward = s%top_conv_bdy(i)

    if (outward) then
       k_a = k_ob
       k_b = 1
       dk = -1
    else
       k_a = k_ob+1
       k_b = s%nz
       dk = 1
    endif

    face_loop : do k = k_a, k_b, dk

       ! Evaluate the step factor

       r = s%r(k)

       if (outward) then
          dr = r - r_ob
       else
          dr = r_ob - r
       endif

       if (dr < f*Hp_cb) then
          factor = 1._dp
       else
          factor = 0._dp
       endif

       ! Store the diffusion coefficient and velocity

       D(k) = (D0 + Delta0*D_ob)*factor
       if(D_ob /= 0d0) then
          vc(k) = (D0/D_ob + Delta0)*vc_ob*factor
       else
          vc(k) = 0d0
       end if
       ! Check for early overshoot completion

       if (D(k) < s%overshoot_D_min) then
          k_b = k
          exit face_loop
       endif
          
    end do face_loop

    ! Finish

    ierr = 0

    return

  end subroutine eval_overshoot_step

end module overshoot_step
