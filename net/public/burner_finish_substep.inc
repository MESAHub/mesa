      subroutine burner_finish_substep(nstp, time, y, ierr)
         integer,intent(in) :: nstp
         double precision, intent(in) :: time, y(:)
         integer, intent(out) :: ierr
      end subroutine burner_finish_substep
