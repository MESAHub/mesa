      subroutine jakob(x,h,y,dfdy,nvar,ierr)
         use const_def, only: dp
         integer, intent(in) :: nvar
         real(dp) :: x,h,y(:)
         real(dp), pointer :: dfdy(:,:)
         integer, intent(out) :: ierr
      end subroutine jakob
