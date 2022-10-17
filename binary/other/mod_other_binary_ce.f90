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

module mod_other_binary_ce
   use star_def
   use binary_def, only : binary_info, binary_ptr
   
   ! NOTE: remember to set true:
   ! use_other_CE_rlo_mdot = .true.
   
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
   !        b% other_CE_rlo_mdot => my_ce_rlo_mdot
   !      end subroutine extras_controls
   
   !      subroutine my_ce_rlo_mdot(binary_id, rlo_mdot, ierr)
   !         use const_def, only: dp
   !         integer, intent(in) :: binary_id
   !         real(dp), intent(out) :: rlo_mdot
   !         integer, intent(out) :: ierr
   !         type (binary_info), pointer :: b
   !         ierr = 0
   !         call binary_ptr(binary_id, b, ierr)
   !         if (ierr /= 0) then
   !            write(*,*) 'failed in binary_ptr'
   !            return
   !         end if
   !         rlo_mdot = 0d0
   !      end subroutine my_ce_rlo_mdot
   
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
   
   subroutine null_other_CE_init(binary_id, restart, ierr)
      use const_def, only : dp
      integer, intent(in) :: binary_id
      logical, intent(in) :: restart
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      write(*, *) "WARNING: using null_other_CE_init"
      ierr = 0
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then
         write(*, *) 'failed in binary_ptr'
         return
      end if
   end subroutine null_other_CE_init
   
   subroutine null_other_CE_rlo_mdot(binary_id, mdot, ierr)
      use const_def, only : dp
      integer, intent(in) :: binary_id
      real(dp), intent(out) :: mdot
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      write(*, *) "WARNING: using null_other_CE_rlo_mdot"
      ierr = 0
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then
         write(*, *) 'failed in binary_ptr'
         return
      end if
      mdot = -1d-99
   end subroutine null_other_CE_rlo_mdot
   
   integer function null_other_CE_binary_evolve_step(binary_id)
      use const_def, only : dp
      integer, intent(in) :: binary_id
      integer :: ierr
      type (binary_info), pointer :: b
      write(*, *) "WARNING: using null_other_CE_binary_evolve_step"
      ierr = 0
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then
         write(*, *) 'failed in binary_ptr'
         return
      end if
      null_other_CE_binary_evolve_step = keep_going
   end function null_other_CE_binary_evolve_step
   
   integer function null_other_CE_binary_finish_step(binary_id)
      integer, intent(in) :: binary_id
      integer :: ierr
      type (binary_info), pointer :: b
      write(*, *) "WARNING: using null_other_CE_binary_check_model"
      ierr = 0
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then
         write(*, *) 'failed in binary_ptr'
         return
      end if
      null_other_CE_binary_finish_step = keep_going
   end function null_other_CE_binary_finish_step

end module mod_other_binary_ce

