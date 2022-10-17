! ***********************************************************************
!
!   Copyright (C) 2015  Bill Paxton and Pablo Marchant
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

module mod_other_binary_extras
   use binary_def
   use const_def
   
   implicit none
   
   private warn_run_star_extras
   public

contains
   
   integer function null_extras_binary_startup(binary_id, restart, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: binary_id
      integer, intent(out) :: ierr
      logical, intent(in) :: restart
      
      null_extras_binary_startup = keep_going
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then ! failure in  binary_ptr
         return
      end if
   
   end function  null_extras_binary_startup
   
   integer function null_extras_binary_start_step(binary_id, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: binary_id
      integer, intent(out) :: ierr
      
      null_extras_binary_start_step = keep_going
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then ! failure in  binary_ptr
         return
      end if
   
   end function  null_extras_binary_start_step
   
   !Return either keep_going, retry or terminate
   integer function null_extras_binary_check_model(binary_id)
      type (binary_info), pointer :: b
      integer, intent(in) :: binary_id
      integer :: ierr
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then ! failure in  binary_ptr
         return
      end if
      null_extras_binary_check_model = keep_going
   
   end function null_extras_binary_check_model
   
   
   ! returns either keep_going or terminate.
   ! note: cannot request retry; extras_binary_check_model can do that.
   integer function null_extras_binary_finish_step(binary_id)
      type (binary_info), pointer :: b
      integer, intent(in) :: binary_id
      integer :: ierr
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then ! failure in  binary_ptr
         return
      end if
      null_extras_binary_finish_step = keep_going
   
   end function null_extras_binary_finish_step
   
   
   subroutine null_extras_binary_after_evolve(binary_id, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: binary_id
      integer, intent(out) :: ierr
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then ! failure in  binary_ptr
         return
      end if
   
   end subroutine null_extras_binary_after_evolve
   
   integer function null_how_many_extra_binary_history_columns(binary_id)
      integer, intent(in) :: binary_id
      type (binary_info), pointer :: b
      integer :: ierr
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then ! failure in  binary_ptr
         return
      end if
      null_how_many_extra_binary_history_columns = 0
      
      call warn_run_star_extras(b% warn_binary_extra)
   
   end function null_how_many_extra_binary_history_columns
   
   subroutine null_data_for_extra_binary_history_columns(binary_id, n, names, vals, ierr)
      integer, intent(in) :: binary_id
      integer, intent(in) :: n
      character (len = maxlen_binary_history_column_name) :: names(n)
      real(dp) :: vals(n)
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then ! failure in  binary_ptr
         return
      end if
      
      call warn_run_star_extras(b% warn_binary_extra)
   
   end subroutine null_data_for_extra_binary_history_columns
   
   
   integer function null_how_many_extra_binary_history_header_items(binary_id)
      use binary_def, only : binary_info
      integer, intent(in) :: binary_id
      null_how_many_extra_binary_history_header_items = 0
   end function null_how_many_extra_binary_history_header_items
   
   subroutine null_data_for_extra_binary_history_header_items(&
      binary_id, n, names, vals, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: binary_id, n
      character (len = maxlen_binary_history_column_name) :: names(n)
      real(dp) :: vals(n)
      integer, intent(out) :: ierr
      ierr = 0
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then
         write(*, *) 'failed in binary_ptr'
         return
      end if
   end subroutine null_data_for_extra_binary_history_header_items
   
   
   subroutine warn_run_star_extras(warn)
      logical, intent(in) :: warn
      
      if(warn) then
         write(*, *) "WARNING: run_binary_extras has changed"
         write(*, *) "and you are calling a null version of this routine"
         write(*, *) "If you had customized your run_binary_extras.f file those functions are not being called"
         write(*, *) "See $MESA_DIR/binary/work/src/run_binary_extras.f for what you need to do now"
         write(*, *) "To disable this warning set warn_run_binary_extra=.false. in your binary_job inlist"
         write(*, *) "run_star_extras has also changed, so set warn_run_star_extra=.false. in star1's star_job inlist"
         write(*, *) "MESA exited due to run_binary_extras warning."
         stop
      end if
   
   end subroutine warn_run_star_extras


end module mod_other_binary_extras
      
