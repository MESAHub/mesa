!!****** stl/rad_photo_rate [1.0] *
!!
!!   NAME
!!     rad_photo_rate
!!
!!   FUNCTION
!!     Compute the rates of atomic processes
!!
!!   AUTHOR
!!     bakl
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
module nag_sub
 use kinds,                           only: dp, sp
  use termination,            only:   stop_program, stop_memory
#include "cp_common_uses.h"

  implicit none
  private
    
  character(len=*), parameter, private :: mdl_name = 'lnag_sub' 
  
  public :: nag_d01gaf 
  
contains

!!******
!!
!!   NAME
!!     nag_d01gaf
!!
!!   FUNCTION
!!     This subroutine integrates a function (y) specified
!!     numerically at n points (x), where n is at least 4,
!!     over the range x(1) to x(n).  the points need not be
!!     equally spaced, but should be distinct and in ascending
!!     or descending order.  an error estimate is returned.
!!     the method is due to gill and miller.
!!     see SUBROUTINE D01GAF in NAG
!!
!!   INPUTS
!!      see SUBROUTINE D01GAF in NAG
!!
!!   OUTPUT
!!       - 
!!   SOURCE
!!
!!*** **********************************************************************
  subroutine nag_d01gaf (X,y,N,ans,err)
    real (kind=dp), dimension(:), intent(in)   :: X 
    real (kind=dp), dimension(:), intent(in)   ::  y 
    integer, intent(in) :: N 
    real (kind=dp), intent(out)  :: ans, err
    character(len=*), parameter ::  subrtn_name = 'd01gaf' &
                                       , fullPathSubrtn = mdl_name//'.'//subrtn_name
    integer  ::  ifail
 

    call d01gaf(X,y,N,ans,err,ifail)
    if(ifail > 0)  then
        if (ifail.eq.1) then
            call stop_program(subrtn_name,mdl_name,__LINE__, 'less than 4 points supplied');
        else if (ifail.eq.2) then
            call stop_program(subrtn_name,mdl_name,__LINE__, 'points not in increasing or decreasing order');
        else if (ifail.eq.3) then
            call stop_program(subrtn_name,mdl_name,__LINE__, 'points not all distinct');
        end if
    endif
  end subroutine nag_d01gaf

end module nag_sub