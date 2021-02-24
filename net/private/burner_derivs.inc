      subroutine derivs(x,h,y,dydt,nvar,ierr)
         use const_def, only: dp
         integer, intent(in) :: nvar
         real(dp) :: x,h,y(:),dydt(:)
         integer, intent(out) :: ierr
      end subroutine derivs
