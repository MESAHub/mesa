! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton and Pablo Marchant
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

module mod_other_binary_wind_transfer
   use binary_def, only : binary_info, binary_ptr
   
   ! NOTE: remember to set true:
   ! use_other_binary_wind_transfer = .true.
   
   ! you can add your own routine for use instead of the default ones
   
   ! here's how to do it.
   
   ! Before doing anything, let's make sure your working copy of run_binary_extras works.
   ! edit the extras_binary_controls routine
   !      subroutine extras_binary_controls(binary_id, ierr)
   !         integer :: binary_id
   !         integer, intent(out) :: ierr
   !         type (binary_info), pointer :: b
   !         ierr = 0
   !         call binary_ptr(binary_id, b, ierr)
   !         if (ierr /= 0) then
   !            write(*,*) 'failed in binary_ptr'
   !            return
   !         end if
   !        write(*,*) 'hello from extra_binary_controls'
   !      end subroutine extras_binary_controls
   
   ! then, in your work directory, do ./mk and ./rn to check that it is okay.
   ! assuming that worked, now edit extra_binary_controls to set the procedure pointer to other_wind_transfer
   
   !      subroutine extras_binary_controls(binary_id, ierr)
   !         integer :: binary_id
   !         integer, intent(out) :: ierr
   !         type (binary_info), pointer :: b
   !         ierr = 0
   !         call binary_ptr(binary_id, b, ierr)
   !         if (ierr /= 0) then
   !            write(*,*) 'failed in binary_ptr'
   !            return
   !         end if
   !        b% other_binary_wind_transfer => wind_transfer_routine
   !      end subroutine extras_controls
   
   !      subroutine wind_transfer_routine(binary_id, s_i, ierr)
   !         integer, intent(in) :: binary_id, s_i
   !         integer, intent(out) :: ierr
   !         type (binary_info), pointer :: b
   !         ierr = 0
   !         call binary_ptr(binary_id, b, ierr)
   !         if (ierr /= 0) then
   !            write(*,*) 'failed in binary_ptr'
   !            return
   !         end if
   !         ! for example, a fixed fraction of the wind
   !         b% wind_xfer_fraction(s_i) = 0.5
   !      end subroutine wind_transfer_routine
   
   ! A more elaborate example for wind transfer is the Wind Roche Lobe Overflow mechanism.
   ! The implementation of Abate et al. 2013 is:
   
   !       subroutine WRLOF_wind_transfer(binary_id, s_i, ierr)
   !         integer, intent(in) :: binary_id, s_i
   !         integer, intent(out) :: ierr
   !         real(dp) :: r_dust, x, q2, c1, c2, c3
   !         type (binary_info), pointer :: b
   !         ! Wind roche lobe overflow as implemented by:
   !         ! Abate et al. 2013, A&A, 552, A26
   !
   !         ierr = 0
   !         call binary_ptr(binary_id, b, ierr)
   !         if (ierr /= 0) then
   !            write(*,*) 'failed in binary_ptr'
   !            return
   !         end if
   !
   !         ! Dust radius based on Hofner 2007, ASPC, 378, 145
   !         r_dust = 0.5d0 * b% r(s_i) * (b% s_donor % Teff / 1500d0)**2.5d0
   !         x = r_dust / b% rl(s_i)
   !
   !         ! constants from Abate et al. eq. 5
   !         q2 = ( b% m(3-s_i)/b% m(s_i) )**2d0
   !         c1 = -0.284d0
   !         c2 = 0.918d0
   !         c3 = -0.234d0
   !
   !         ! WRLOF transfer fraction from Abate et al. eq. 9
   !         b% wind_xfer_fraction(s_i) = 25d0 / 9d0 * q2 * (c1*x**2 + c2*x + c3)
   !
   !         b% wind_xfer_fraction(s_i) = min(b% wind_xfer_fraction(s_i), 5d-1)
   !
   !      end subroutine WRLOF_wind_transfer
   
   ! NOTE: if you'd like to have some inlist controls for your routine,
   ! you can use the x_ctrl array of real(dp) variables that is in &controls
   ! e.g., in the &controls inlist, you can set
   !     x_ctrl(1) = my_special_param
   ! then in your routine, you can access that by
   !     s% x_ctrl(1)
   ! of course before you can use s, you need to get it using the id argument.
   ! here's an example of how to do that -- add these lines at the start of your routine:
   !
   !         use star_lib, only: star_ptr
   !         type (star_info), pointer :: s
   !         call star_ptr(id, s, ierr)
   !         if (ierr /= 0) then ! OOPS
   !            return
   !         end if
   !
   ! To get the binary pointer using the provided binary_id, add these lines.
   !
   !      type (binary_info), pointer :: b
   !      call binary_ptr(binary_id, b, ierr)
   !      if (ierr /= 0) then ! failure in  binary_ptr
   !         return
   !      end if
   !
   ! for integer control values, you can use x_integer_ctrl
   ! for logical control values, you can use x_logical_ctrl
   
   implicit none


contains
   
   subroutine null_other_binary_wind_transfer(binary_id, s_i, ierr)
      integer, intent(in) :: binary_id, s_i
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then
         write(*, *) 'failed in binary_ptr'
         return
      end if
      b% wind_xfer_fraction(s_i) = 0
   end subroutine null_other_binary_wind_transfer

end module mod_other_binary_wind_transfer

