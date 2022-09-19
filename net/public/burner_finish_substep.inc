      subroutine burner_finish_substep(nstp, time, y, ierr)
         use const_def, only: dp
         integer,intent(in) :: nstp
         real(dp), intent(in) :: time, y(:)
         integer, intent(out) :: ierr
      end subroutine burner_finish_substep
