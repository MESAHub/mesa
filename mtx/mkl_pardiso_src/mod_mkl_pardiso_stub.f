
      module mod_mkl_pardiso
      
      implicit none

      integer, parameter :: num_pardiso_ipar_decsol = 0
      integer, parameter :: num_pardiso_rpar_decsol = 0

      
      contains
      
      
      logical function use_mkl_pardiso()
         use_mkl_pardiso = .false.
      end function use_mkl_pardiso

      
      subroutine do_mkl_pardiso_work_sizes(n,nzmax,lrd,lid)
         integer, intent(in) :: n,nzmax
         integer, intent(out) :: lrd,lid
         lid = 0
         lrd = 0
      end subroutine do_mkl_pardiso_work_sizes

      
      subroutine do_mkl_pardiso(iop,n,nzmax,ia,ja,sa,b,lrd,rpar_decsol,lid,ipar_decsol,ierr)
         integer, intent(in) :: iop, n, nzmax, lrd, lid
         integer, intent(inout) :: ia(n+1), ja(nzmax)
         double precision, intent(inout) :: sa(nzmax)
         double precision, intent(inout) :: b(n)
         double precision, intent(inout), target :: rpar_decsol(lrd)
         integer, intent(inout), target :: ipar_decsol(lid)
         integer, intent(out) :: ierr
         ierr = -1
         write(*,*) 'mkl pardiso not loaded'
         call mesa_error(__FILE__,__LINE__)
      end subroutine do_mkl_pardiso

      
      subroutine do_mkl_pardiso_nrhs(iop,nrhs,n,nzmax,ia,ja,sa,b,lrd,rpar_decsol,lid,ipar_decsol,ierr)
         integer, intent(in) :: iop, nrhs, n, nzmax, lrd, lid
         integer, intent(inout) :: ia(n+1), ja(nzmax)
         double precision, intent(inout) :: sa(nzmax)
         double precision, intent(inout) :: b(n)
         double precision, intent(inout), target :: rpar_decsol(lrd)
         integer, intent(inout), target :: ipar_decsol(lid)
         integer, intent(out) :: ierr
         ierr = -1
         write(*,*) 'mkl pardiso not loaded'
         call mesa_error(__FILE__,__LINE__)
      end subroutine do_mkl_pardiso_nrhs


      end module mod_mkl_pardiso
