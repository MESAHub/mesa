      subroutine jakob(x,y,dfdy,nvar,ierr,k)
         use const_def, only: dp
         integer, intent(in) :: nvar
         real(dp) :: x,y(:)
         real(dp) :: dfdy(:,:)
         integer, intent(out) :: ierr
         integer, intent(in) :: k
      end subroutine jakob
