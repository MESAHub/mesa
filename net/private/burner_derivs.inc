      subroutine derivs(x,y,dydt,nvar,ierr,k)
         use const_def, only: dp
         integer, intent(in) :: nvar
         real(dp) :: x,y(:),dydt(:)
         integer, intent(out) :: ierr
         integer, intent(in) :: k
      end subroutine derivs
