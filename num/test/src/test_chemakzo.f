      module test_chemakzo
      use num_def
      use num_lib
      use test_int_support, only: i_nfcn, i_njac
      use utils_lib, only: mesa_error

      implicit none

      contains


      subroutine chemakzo_derivs(n, x, h, vars, dvars_dx, lrpar,rpar,lipar,ipar, ierr)
         integer, intent(in) :: n, lrpar, lipar
         real(dp), intent(in) :: x, h
         real(dp), intent(inout) :: vars(:) ! (n)
         real(dp), intent(inout) :: dvars_dx(:) ! (n)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         real(dp) :: yprime(n)
         integer, intent(out) :: ierr
         ierr = 0
         ipar(i_nfcn) = ipar(i_nfcn) + 1
         call chemakzo_feval(n,x,vars,yprime,dvars_dx,ierr,rpar,ipar)
      end subroutine chemakzo_derivs


      subroutine chemakzo_jacob(n,x,h,y,f,dfdy,ld_dfdy,lrpar,rpar,lipar,ipar,ierr)
         integer, intent(in) :: n, ld_dfdy, lrpar, lipar
         real(dp), intent(in) :: x, h
         real(dp), intent(inout) :: y(:) ! (n)
         real(dp), intent(inout) :: f(:), dfdy(:,:) ! (ld_dfdy,n)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         real(dp) :: yprime(n)
         integer, intent(out) :: ierr
         integer :: nz, i, j
         ierr = 0
         ipar(i_njac) = ipar(i_njac) + 1
         call chemakzo_jeval(ld_dfdy,n,x,y,yprime,dfdy,ierr,rpar,ipar)
         if (ierr == 0) call chemakzo_derivs(n, x, h, y, f, lrpar,rpar,lipar,ipar, ierr)
      end subroutine chemakzo_jacob


      subroutine chemakzo_sjac(n,x,h,y,f,nzmax,ia,ja,values,lrpar,rpar,lipar,ipar,ierr)  
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
         call chemakzo_jacob(n,x,h,y,f,dfdy,ld_dfdy,lrpar,rpar,lipar,ipar,ierr)
         if (ierr /= 0) return
         if (ipar(ipar_sparse_format) == 0) then
            call dense_to_row_sparse_with_diag(n,n,dfdy,nzmax,nz,ia,ja,values,ierr)
         else
            call dense_to_col_sparse_with_diag(n,n,dfdy,nzmax,nz,ia,ja,values,ierr)
         end if
      end subroutine chemakzo_sjac


      subroutine chemakzo_solout(nr,xold,x,n,y,rwork,iwork,interp_y,lrpar,rpar,lipar,ipar,irtrn)
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
         real(dp) :: xout, val
         integer :: i, ierr
         include 'formats'
         irtrn = 0
         !write(*,2) 'x', nr, x
         xout = 100
         if (x >= xout .and. xout > xold) then
            write(*,*)
            write(*,1) 'dense output for x=', xout
            ierr = 0
            do i=1,n
               val = interp_y(i,xout,rwork,iwork,ierr)
               if (ierr /= 0) then
                  write(*,2) 'interp_y failed', i, val
               else
                  write(*,2) 'val', i, val
               end if
            end do
            write(*,*)
         end if
         
      end subroutine chemakzo_solout


      subroutine chemakzo_mas_band(n,am,lmas,lrpar,rpar,lipar,ipar)
         integer, intent(in) :: n, lmas, lrpar, lipar
         real(dp), intent(inout) :: am(:,:) ! (lmas,n)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         real(dp) :: yprime(n), t
         integer :: ierr
         ierr = 0
         call chemakzo_meval(lmas,n,t,yprime,am,ierr,rpar,ipar)
      end subroutine chemakzo_mas_band


      subroutine chemakzo_mas_full(n,am,lmas,lrpar,rpar,lipar,ipar)
         integer, intent(in) :: n, lmas, lrpar, lipar
         real(dp), intent(inout) :: am(:,:) ! (lmas,n)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         real(dp) :: yprime(n), t
         integer :: ierr, i
         ierr = 0
         call chemakzo_meval(lmas,n,t,yprime,am,ierr,rpar,ipar)
         do i=2,n
            am(i,i) = am(1,i)
            am(1,i) = 0
         end do
      end subroutine chemakzo_mas_full
      
      
      subroutine do_test_chemakzo(which_solver,which_decsol,m_band,numerical_jacobian,show_all,quiet)
         use test_support,only:show_results,show_statistics,check_results
         use test_int_support,only:do_test_stiff_int
         integer, intent(in) :: which_solver,which_decsol
         logical, intent(in) :: m_band,numerical_jacobian,show_all,quiet

         integer, parameter :: n = 6 ! the number of variables in the "chemakzo" system of ODEs
         real(dp), target :: y_ary(n), yprime(n), yexact(n)
         real(dp), pointer :: y(:)
         integer, parameter :: lrpar = 1, lipar = 3, iout=2 ! test dense output
         logical :: consis
         integer, parameter :: ndisc = 0
         real(dp) :: h0, t(0:ndisc+1), atol(1), rtol(1)
         integer :: i, mujac, mljac, matrix_type_spec, ierr, imas, mlmas, mumas, m1, m2, itol, nstep
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
            
         if (.not. quiet) write(*,*) 'chemakzo'
         
         t(0)   = 0
         t(1)   = 180d0
         
         itol = 0 ! scalar tolerances
         rtol(1) = 1d-8
         atol(1) = 1d-8
         h0 = 1d-10 ! initial step size
         
         mljac = n ! square matrix
         mujac = n
         matrix_type_spec = square_matrix_type

         imas = 1
         m1 = 0
         m2 = 0     
         
         call chemakzo_init(n,y,yprime,consis)
         nstep=0   
         if (m_band) then
            write(*,*) 'M banded'
            ! mass matrix is diagonal
            mlmas = 0 
            mumas = 0
            call do_test_stiff_int(which_solver,which_decsol,numerical_jacobian, &
               chemakzo_derivs,chemakzo_jacob,chemakzo_sjac,chemakzo_solout,iout, &
               null_fcn_blk_dble,null_jac_blk_dble, &
               caller_id,nvar,nz,lblk,dblk,ublk,uf_lblk,uf_dblk,uf_ublk, &
               n,ndisc,mljac,mujac,matrix_type_spec,chemakzo_mas_band,imas,mlmas,mumas,m1,m2, &
               t,rtol,atol,itol,h0,y,nstep,lrpar,rpar,lipar,ipar,quiet,ierr)
         else
            write(*,*) 'M full'
            mlmas = n 
            mumas = n 
            call do_test_stiff_int(which_solver,which_decsol,numerical_jacobian, &
               chemakzo_derivs,chemakzo_jacob,chemakzo_sjac,chemakzo_solout,iout, &
               null_fcn_blk_dble,null_jac_blk_dble, &
               caller_id,nvar,nz,lblk,dblk,ublk,uf_lblk,uf_dblk,uf_ublk, &
               n,ndisc,mljac,mujac,matrix_type_spec,chemakzo_mas_full,imas,mlmas,mumas,m1,m2, &
               t,rtol,atol,itol,h0,y,nstep,lrpar,rpar,lipar,ipar,quiet,ierr)
         end if
         if (ierr /= 0) then
            write(*,*) 'chemakzo ierr', ierr
            call mesa_error(__FILE__,__LINE__)
         end if
         
         call chemakzo_solut(n,0d0,yexact)
         call check_results(n,y,yexact,rtol(1)*10,ierr)
         if (ierr /= 0) then
            write(*,*) 'check results ierr', ierr
            call mesa_error(__FILE__,__LINE__) ! do_test_vdpol
         end if
         
         if (quiet) return

         call show_results(n,y,yexact,show_all)
         call show_statistics(ipar(i_nfcn),ipar(i_njac),nstep,show_all)

      end subroutine do_test_chemakzo
     
      end module test_chemakzo
