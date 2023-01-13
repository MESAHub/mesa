      subroutine jakob(x,y,dfdy,nvar,ierr)
         use const_def, only: dp
         integer, intent(in) :: nvar
         real(dp) :: x,y(:)
         real(dp) :: dfdy(:,:)
         integer, intent(out) :: ierr
      end subroutine jakob
