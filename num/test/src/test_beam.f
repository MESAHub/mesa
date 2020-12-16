      module test_beam
      use num_def
      use num_lib
      use test_int_support, only: i_nfcn, i_njac
      use utils_lib, only: mesa_error

      implicit none

      contains


      subroutine beam_derivs(n, x, h, vars, dvars_dx, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: n, lrpar, lipar
         real(dp), intent(in) :: x, h
         real(dp), intent(inout) :: vars(:) ! (n)
         real(dp), intent(inout) :: dvars_dx(:) ! (n)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         ierr = 0
         ipar(i_nfcn) = ipar(i_nfcn) + 1
         call beam_feval(n,x,vars,dvars_dx,ierr,rpar,ipar)
      end subroutine beam_derivs


      subroutine beam_jacob(n,x,h,y,f,dfdy,ld_dfdy,lrpar,rpar,lipar,ipar,ierr)
         integer, intent(in) :: n, ld_dfdy, lrpar, lipar
         real(dp), intent(in) :: x, h
         real(dp), intent(inout) :: y(:) ! (n)
         real(dp), intent(inout) :: dfdy(:,:) !(ld_dfdy,n)
         real(dp), intent(inout) :: f(:) !(n)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         real(dp) :: yprime(n)
         integer, intent(out) :: ierr
         integer :: nz, i, j
         ierr = 0
         ipar(i_njac) = ipar(i_njac) + 1
         call beam_jeval(ld_dfdy,n,x,y,yprime,dfdy,ierr,rpar,ipar)
         if (ierr == 0) call beam_derivs(n, x, h, y, f, lrpar, rpar, lipar, ipar, ierr)
      end subroutine beam_jacob


      subroutine beam_sjac(n,x,h,y,f,nzmax,ia,ja,values,lrpar,rpar,lipar,ipar,ierr)  
         ! sparse jacobian. format either compressed row or compressed column.
         use mtx_lib,only:dense_to_row_sparse_with_diag,dense_to_col_sparse_with_diag
         use test_int_support,only:ipar_sparse_format
         integer, intent(in) :: n, nzmax, lrpar, lipar
         real(dp), intent(in) :: x, h
         real(dp), intent(inout) :: y(:) ! (n)
         real(dp), intent(inout) :: f(:) ! (n) ! dy/dx
         integer, intent(inout) :: ia(:) ! (n+1)
         integer, intent(inout) :: ja(:) ! (nzmax)
         real(dp), intent(inout) :: values(:) ! (nzmax)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr ! nonzero means terminate integration
         real(dp) :: dfdy(n,n)
         integer :: ld_dfdy, nz
         ld_dfdy = n
         ierr = 0
         call beam_jacob(n,x,h,y,f,dfdy,ld_dfdy,lrpar,rpar,lipar,ipar,ierr)
         if (ierr /= 0) return
         if (ipar(ipar_sparse_format) == 0) then
            call dense_to_row_sparse_with_diag(n,n,dfdy,nzmax,nz,ia,ja,values,ierr)
         else
            call dense_to_col_sparse_with_diag(n,n,dfdy,nzmax,nz,ia,ja,values,ierr)
         end if
      end subroutine beam_sjac


      subroutine beam_solout(nr,xold,x,n,y,rwork,iwork,interp_y,lrpar,rpar,lipar,ipar,irtrn)
         ! nr is the step number.
         ! x is the current x value; xold is the previous x value.
         ! y is the current y value.
         ! irtrn negative means terminate integration.
         ! rwork and iwork hold info for 
         integer, intent(in) :: nr, n, lrpar, lipar
         real(dp), intent(in) :: xold, x
         real(dp), intent(inout) :: y(:) ! (n)
         real(dp), intent(inout), target :: rwork(*)
         integer, intent(inout), target :: iwork(*)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         interface
            ! this subroutine can be called from your solout routine.
            ! it computes interpolated values for y components during the just completed step.
            real(dp) function interp_y(i,s,rwork,iwork,ierr)
               use const_def, only: dp
               integer, intent(in) :: i ! result is interpolated approximation of y(i) at x=s.
               real(dp), intent(in) :: s ! interpolation x value (between xold and x).
               real(dp), intent(inout), target :: rwork(*)
               integer, intent(inout), target :: iwork(*)
               integer, intent(out) :: ierr
            end function interp_y
         end interface
         integer, intent(out) :: irtrn
         irtrn = 0
      end subroutine beam_solout
      
      
      subroutine do_test_beam(which_solver,which_decsol,numerical_jacobian,show_all,quiet)
         use test_support,only:show_results,show_statistics,check_results
         use test_int_support,only:do_test_stiff_int
         integer, intent(in) :: which_solver,which_decsol
         logical, intent(in) :: numerical_jacobian,show_all,quiet

         integer, parameter :: n = 80 ! the number of variables in the "beam" system of ODEs
         real(dp), target :: y_ary(n), yprime(n), yexact(n)
         real(dp), pointer :: y(:)
         integer, parameter :: lrpar = 1, lipar = 3, iout=1
         logical :: consis
         integer, parameter :: ndisc = 0, n_soln=9
         real(dp) :: result(n_soln), soln(n_soln)
         real(dp) :: h0, t(0:ndisc+1), atol(1), rtol(1)
         integer :: i, mujac, mljac, matrix_type_spec, ierr, indsol(n_soln), imas, mlmas, mumas, m1, m2, itol, nstep
         real(dp), target :: rpar_ary(lrpar) 
         integer, target :: ipar_ary(lipar)
         real(dp), pointer :: rpar(:)
         integer, pointer :: ipar(:)
         integer :: caller_id, nvar, nz
         real(dp), dimension(:), pointer :: lblk, dblk, ublk ! =(nvar,nvar,nz)
         real(dp), dimension(:), pointer :: uf_lblk, uf_dblk, uf_ublk ! =(nvar,nvar,nz)
         
         nullify(lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk)
         caller_id = 0
         nvar = 0
         nz = 0
         
         rpar => rpar_ary
         ipar => ipar_ary
         y => y_ary
            
         if (.not. quiet) write(*,*) 'beam'
         
         t(0)   = 0
         t(1)   = 5d0
         
         itol = 0 ! scalar tolerances
         rtol(1) = 1d-3
         atol(1) = 1d-3
         h0 = 1d-4 ! initial step size
         
         m1 = n/2
         m2 = 0     
         
         mljac = n
         mujac = n
         matrix_type_spec = square_matrix_type

         imas = 0
         mlmas = 0
         mumas = 0
         
         if (.not. numerical_jacobian) then
            write(*,*) 'beam test only supports numerical jacobian'
            return
         end if
         
         call beam_init(n,y,yprime,consis)
         nstep=0   
         call do_test_stiff_int(which_solver,which_decsol,numerical_jacobian, &
               beam_derivs,beam_jacob,beam_sjac,beam_solout,iout, &
               null_fcn_blk_dble,null_jac_blk_dble, &
               caller_id,nvar,nz,lblk,dblk,ublk,uf_lblk,uf_dblk,uf_ublk, &
               n,ndisc,mljac,mujac,matrix_type_spec,null_mas,imas,mlmas,mumas,m1,m2, &
               t,rtol,atol,itol,h0,y,nstep,lrpar,rpar,lipar,ipar,quiet,ierr)
         if (ierr /= 0) then
            write(*,*) 'test_beam ierr', ierr
            call mesa_error(__FILE__,__LINE__)
         end if
         
         call beam_solut(n,0d0,yexact)
         indsol(1:n_soln) = (/ 1, 10, 20, 30, 40, 50, 60, 70, 80 /)
         do i=1,n_soln
            result(i) = y(indsol(i))
            soln(i) = yexact(indsol(i))
         end do
         
         call check_results(n,y,yexact,rtol(1)*1d2,ierr)
         if (ierr /= 0) then
            write(*,*) 'check results ierr', ierr
            call mesa_error(__FILE__,__LINE__) ! do_test_vdpol
         end if
         
         if (quiet) return
         
         call show_results(n_soln,result,soln,show_all)
         call show_statistics(ipar(i_nfcn),ipar(i_njac),nstep,show_all)

      end subroutine do_test_beam
            
      
      end module test_beam
