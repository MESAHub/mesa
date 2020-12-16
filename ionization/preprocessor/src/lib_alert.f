!  
! alert
! francois hebert, 09/09/02   
!
! this subroutine offers a quick and dirty way to print a message to stdout
! and either stop, or not, the program
!
! default behavior:
! flag = 1 is an error  (usually stop)
! flag = 0 is a warning (usually no stop)
!

      module lib_alert

      implicit none

      contains


      subroutine alert(flag, message)

         integer, intent(in) :: flag
         character (*), intent(in) :: message

         !
         ! whether to stop on warning/error
         !
         logical, parameter :: stop_on_warning = .false.
         logical, parameter :: stop_on_error = .true.

         if (flag == 0) then
            write (*,*) 'WARNING: ', trim(message)
            if (stop_on_warning) call mesa_error(__FILE__,__LINE__)
         else if (flag == 1) then
            write (*,*) 'ERROR:   ', trim(message)
            if (stop_on_error) call mesa_error(__FILE__,__LINE__)
         end if

      end subroutine alert


      end module lib_alert
