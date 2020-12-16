! ***********************************************************************
!
!   Copyright (C) 2015  Bill Paxton
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

   module other_extras
      use star_def
      use star_lib
      use const_def
      use utils_lib

      implicit none

      private warn_run_star_extras
      public

      contains

      ! Each null routine should have a call to warn_run_star_extras
      ! This adds a warning message to let people know they are
      ! calling the null_ version and not the version in their
      ! run_star_extras.f file

      subroutine null_extras_startup(id,restart,ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         logical, intent(in) :: restart
         call star_ptr(id, s, ierr)
         if (ierr /= 0) then ! failure in  ptr
            return
         end if

         call warn_run_star_extras(s%job%warn_run_star_extras, "extras_startup")

      end subroutine null_extras_startup


      ! return either retry, keep_going or terminate
      integer function null_extras_check_model(id)
         type (star_info), pointer :: s
         integer, intent(in) :: id
         integer :: ierr
         call star_ptr(id, s, ierr)
         if (ierr /= 0) then ! failure in  ptr
            return
         end if
         null_extras_check_model = keep_going

         call warn_run_star_extras(s%job%warn_run_star_extras, "extras_check_model")

      end function null_extras_check_model


      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function null_extras_start_step(id)
         type (star_info), pointer :: s
         integer, intent(in) :: id
         integer :: ierr
         call star_ptr(id, s, ierr)
         if (ierr /= 0) then ! failure in  ptr
            return
         end if
         null_extras_start_step = keep_going

         !call warn_run_star_extras(s%job%warn_run_star_extras, "extras_start_step")

      end function null_extras_start_step


      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function null_extras_finish_step(id)
         type (star_info), pointer :: s
         integer, intent(in) :: id
         integer :: ierr
         call star_ptr(id, s, ierr)
         if (ierr /= 0) then ! failure in  ptr
            return
         end if
         null_extras_finish_step = keep_going

         call warn_run_star_extras(s%job%warn_run_star_extras, "extras_finish_step")

      end function null_extras_finish_step


      subroutine null_extras_after_evolve(id, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         call star_ptr(id, s, ierr)
         if (ierr /= 0) then ! failure in  ptr
            return
         end if

         call warn_run_star_extras(s%job%warn_run_star_extras, "extras_after_evolve")

      end subroutine null_extras_after_evolve


      integer function null_how_many_extra_history_columns(id)
         integer, intent(in) :: id
         type (star_info), pointer :: s
         integer :: ierr
         call star_ptr(id, s, ierr)
         if (ierr /= 0) then ! failure in  ptr
            return
         end if
         null_how_many_extra_history_columns=0

         call warn_run_star_extras(s%job%warn_run_star_extras, "how_many_extra_history_columns")

      end function null_how_many_extra_history_columns


      subroutine null_data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id
         integer, intent(in) :: n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         call star_ptr(id, s, ierr)
         if (ierr /= 0) then ! failure in  ptr
            return
         end if

         call warn_run_star_extras(s%job%warn_run_star_extras, "data_for_extra_history_columns")

      end subroutine null_data_for_extra_history_columns


      integer function null_how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         type (star_info), pointer :: s
         integer :: ierr
         call star_ptr(id, s, ierr)
         if (ierr /= 0) then ! failure in  ptr
            return
         end if
         null_how_many_extra_profile_columns=0

         call warn_run_star_extras(s%job%warn_run_star_extras, "how_many_extra_profile_columns")

      end function null_how_many_extra_profile_columns


      subroutine null_data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id
         integer, intent(in) :: n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         call star_ptr(id, s, ierr)
         if (ierr /= 0) then ! failure in  ptr
            return
         end if

         call warn_run_star_extras(s%job%warn_run_star_extras, "data_for_extra_profile_columns")

      end subroutine null_data_for_extra_profile_columns


      subroutine warn_run_star_extras(warn, routine)
         logical, intent(in) :: warn
         character(len=*), intent(in) :: routine

         if (warn) then
            write(*,*) "WARNING: you are calling a null version of " // trim(routine)
            write(*,*)
            write(*,*) "If you have customized your run_star_extras.f file, your routine is not being called."
            write(*,*) "Check that your extras_controls has a line like"
            write(*,*) "    s% " // trim(routine) // " => " // trim(routine)
            write(*,*) "See $MESA_DIR/star/job/standard_run_star_extras.inc for the standard example."
            write(*,*)
            write(*,*) "This error can also occur if you switched MESA versions without recompiling."
            write(*,*) "Do a ./clean and ./mk in your work directory to get recompiled."
            write(*,*)
            write(*,*) "To disable this warning set"
            write(*,*) "    warn_run_star_extras = .false."
            write(*,*) "in your star_job inlist."
            write(*,*)
            write(*,*) "MESA exited due to run_star_extras warning."
            stop
         end if

      end subroutine warn_run_star_extras


   end module other_extras
