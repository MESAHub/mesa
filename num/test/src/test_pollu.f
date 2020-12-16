      module test_pollu
      use num_def
      use num_lib
      use test_int_support,only:i_nfcn,i_njac
      implicit none

      contains


      subroutine pollu_derivs(n, x, h, vars, dvars_dx, lrpar,rpar,lipar,ipar, ierr)
         integer, intent(in) :: n, lrpar, lipar
         double precision, intent(in) :: x, h
         double precision, intent(inout) :: vars(n)
         double precision, intent(out) :: dvars_dx(n)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         double precision :: yprime(n)
         integer, intent(out) :: ierr
         ierr = 0
         ipar(i_nfcn) = ipar(i_nfcn) + 1
         call pollu_feval(n,x,h,vars,yprime,dvars_dx,ierr,rpar,ipar)
      end subroutine pollu_derivs


      subroutine pollu_jacob(n,x,h,y,f,dfdy,ld_dfdy,lrpar,rpar,lipar,ipar,ierr)
         integer, intent(in) :: n, ld_dfdy, lrpar, lipar
         double precision, intent(in) :: x,h
         double precision, intent(inout) :: y(n)
         double precision, intent(out) :: f(n), dfdy(ld_dfdy,n)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         integer :: nzo, i, j
         double precision :: yprime(n)
         ierr = 0
         ipar(i_njac) = ipar(i_njac) + 1
         call pollu_jeval(ld_dfdy,n,x,y,yprime,dfdy,ierr,rpar,ipar) 
         if (ierr == 0) call pollu_derivs(n, x, y, f, lrpar,rpar,lipar,ipar, ierr)      
      end subroutine pollu_jacob


      subroutine pollu_sjac(n,x,h,y,f,nzmax,ia,ja,values,lrpar,rpar,lipar,ipar,ierr)  
         ! sparse jacobian. format either compressed row or compressed column.
         use mtx_lib,only:dense_to_row_sparse_with_diag,dense_to_col_sparse_with_diag
         use test_int_support,only:ipar_sparse_format
         integer, intent(in) :: n, nzmax, lrpar, lipar
         double precision, intent(in) :: x, h
         double precision, intent(inout) :: y(n)
         integer, intent(out) :: ia(n+1), ja(nzmax)
         double precision, intent(out) :: f(n), values(nzmax)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr ! nonzero means terminate integration
         double precision :: dfdy(n,n)
         integer :: ld_dfdy, nz
         ld_dfdy = n
         ierr = 0
         call pollu_jacob(n,x,h,y,f,dfdy,ld_dfdy,lrpar,rpar,lipar,ipar,ierr)
         if (ierr /= 0) return
         if (ipar(ipar_sparse_format) == 0) then
            call dense_to_row_sparse_with_diag(n,n,dfdy,nzmax,nz,ia,ja,values,ierr)
         else
            call dense_to_col_sparse_with_diag(n,n,dfdy,nzmax,nz,ia,ja,values,ierr)
         end if
      end subroutine pollu_sjac


      subroutine pollu_solout(nr,xold,x,n,y,rwork,iwork,interp_y,lrpar,rpar,lipar,ipar,irtrn)
         ! nr is the step number.
         ! x is the current x value; xold is the previous x value.
         ! y is the current y value.
         ! irtrn negative means terminate integration.
         ! rwork and iwork hold info for 
         integer, intent(in) :: nr, n, lrpar, lipar
         double precision, intent(in) :: xold, x
         double precision, intent(inout) :: y(n)
         real(dp), intent(inout), target :: rwork(*)
         integer, intent(inout), target :: iwork(*)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         interface
            ! this subroutine can be called from your solout routine.
            ! it computes interpolated values for y components during the just completed step.
            double precision function interp_y(i,s,rwork,iwork,ierr)
               integer, intent(in) :: i ! result is interpolated approximation of y(i) at x=s.
               double precision, intent(in) :: s ! interpolation x value (between xold and x).
               double precision, intent(inout), target :: rwork(*)
               integer, intent(inout), target :: iwork(*)
               integer, intent(out) :: ierr
            end function interp_y
         end interface
         integer, intent(out) :: irtrn
         irtrn = 0
         !write(*,*) nr
      end subroutine pollu_solout
      
      
      subroutine do_test_pollu(which_solver,which_decsol,numerical_jacobian,show_all,quiet)
         use test_support,only:show_results,show_statistics,check_results
         use test_int_support,only:do_test_stiff_int
         integer, intent(in) :: which_solver,which_decsol
         logical, intent(in) :: numerical_jacobian,show_all,quiet

         integer, parameter :: n = 20 ! the number of variables in the "pollu" system of ODEs
         integer, parameter :: ndisc = 0, nzo = 86
         double precision :: y(n), yprime(n), yexact(n)
         integer, parameter :: lrpar = 1, lipar = 3, iout=1
         logical :: consis
         integer, parameter :: licn = 5*nzo, lirn = 3*nzo
         double precision :: h0, t(0:ndisc+1), atol(1), rtol(1)
         integer :: i, mujac, mljac, matrix_type_spec, ierr, imas, mlmas, mumas, m1, m2, itol, nstep
         integer :: ivect(nzo), jvect(nzo)
         real(dp), target :: rpar_ary(lrpar) 
         integer, target :: ipar_ary(lipar)
         real(dp), pointer :: rpar(:)
         integer, pointer :: ipar(:)
         integer :: caller_id, nvar, nz
         real(dp), dimension(:), pointer :: lblk, dblk, ublk ! =(nvar,nvar,nz)
         
         nullify(lblk, dblk, ublk)
         caller_id = 0
         nvar = 0
         nz = 0
                  
         rpar => rpar_ary
         ipar => ipar_ary

         if (.not. quiet) write(*,*) 'pollu'

         t(0)   = 0
         t(1)   = 60d0
         
         itol = 0 ! scalar tolerances
         rtol(1) = 1d-5
         atol(1) = 1d-5
         h0 = 1d-7 ! initial step size
         
         call pollu_init(n,y,yprime,consis)
         nstep=0   
         mljac = n ! square matrix
         mujac = n
         matrix_type_spec = square_matrix_type

         imas = 0
         mlmas = 0
         mumas = 0        
         
         m1 = 0
         m2 = 0     

         call do_test_stiff_int(which_solver,which_decsol,numerical_jacobian, &
            pollu_derivs,pollu_jacob,pollu_sjac,pollu_solout,iout, &
            null_fcn_blk_dble,null_jac_blk_dble, &
            caller_id,nvar,nz,lblk,dblk,ublk, &
            n,ndisc,mljac,mujac,matrix_type_spec,null_mas,imas,mlmas,mumas,m1,m2, &
            t,rtol,atol,itol,h0,y,nstep,lrpar,rpar,lipar,ipar,quiet,ierr)
         if (ierr /= 0) then
            write(*,*) 'test_pollu ierr', ierr
            call mesa_error(__FILE__,__LINE__)
         end if
      
         call pollu_solut(n,0d0,yexact)
         call check_results(n,y,yexact,rtol(1)*2,ierr)
         if (ierr /= 0) then
            write(*,*) 'check results ierr', ierr
            call mesa_error(__FILE__,__LINE__) ! do_test_vdpol
         end if
         
         if (quiet) return
         
         call show_results(n,y,yexact,show_all)
         call show_statistics(ipar(i_nfcn),ipar(i_njac),nstep,show_all)

      end subroutine do_test_pollu
      
      
      end module test_pollu
