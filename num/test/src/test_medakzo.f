      module test_medakzo
      use num_def
      use num_lib
      use mtx_lib
      use mtx_def
      use test_int_support, only: i_nfcn, i_njac
      use utils_lib, only: mesa_error

      implicit none
      
      integer :: mljac, mujac

      contains
      
      
      subroutine medakzo_feval_for_blk_dble(neqn,t,y,yprime,f,ierr,rpar,ipar)
      integer neqn,ierr,ipar(:)
      double precision t,y(:),yprime(:),f(:),rpar(:)

      integer N,i,j
      double precision zeta,dzeta,dzeta2,k,c,phi,alpha,beta,gama,dum
      parameter(k=100d0,c=4d0)
      
      include 'formats'

      N      = neqn/2
      dzeta  = 1d0/dble(N)
      dzeta2 = dzeta*dzeta
      dum    = (dzeta-1d0)*(dzeta-1d0)/c
      alpha  = 2d0*(dzeta-1d0)*dum/c
      beta   = dum*dum

      if (t.le.5d0) then
         phi = 2d0
      else
         phi = 0d0
      endif

      f(1) = (phi-2d0*y(1)+y(3))*beta/dzeta2+alpha*(y(3)-phi)/(2d0*dzeta)-k*y(1)*y(2)
      f(2) = -k*y(1)*y(2)

      do 10 j=2,N-1
         i     = 2*j-1
         zeta  = j*dzeta
         dum   = (zeta-1d0)*(zeta-1d0)/c
         alpha = 2d0*(zeta-1d0)*dum/c
         beta  = dum*dum
         gama = (y(i-2)-2d0*y(i)+y(i+2))*beta/dzeta2+alpha*(y(i+2)-y(i-2))/(2d0*dzeta)
         f(i) = gama-k*y(i)*y(i+1)
         i     = 2*j
         f(i) = -k*y(i)*y(i-1)         
   10 continue

      f(2*N-1) = -k*y(2*N-1)*y(2*N)
      f(2*N)   = -k*y(2*N-1)*y(2*N)

      return
      end subroutine medakzo_feval_for_blk_dble


      subroutine medakzo_jeval_for_blk_dble(ldim,neqn,t,y,yprime,dfdy,ierr,rpar,ipar)
      integer ldim,neqn,ierr,ipar(:)
      double precision t,y(:),yprime(:),dfdy(:,:),rpar(:)

      integer N,i,j
      double precision zeta,dzeta,dzeta2,alpha,beta,k,c,dum,bz
      parameter(k=100d0,c=4d0)

      do 20 j=1,neqn
         do 10 i=1,5
            dfdy(i,j) = 0d0
   10    continue
   20 continue

      N      = neqn/2
      dzeta  = 1d0/dble(N)
      dzeta2 = dzeta*dzeta
      dum    = (dzeta-1d0)*(dzeta-1d0)/c
      alpha  = 2d0*(dzeta-1d0)*dum/c
      beta   = dum*dum

      dfdy(3,1) = -beta*2d0/dzeta2-k*y(2)
      dfdy(1,3) = beta/dzeta2+alpha/(2d0*dzeta)
      dfdy(2,2) = -k*y(1)
      dfdy(4,1) = -k*y(2)
      dfdy(3,2) = -k*y(1)

      do 30 j=2,N-1
         i          = 2*j-1
         zeta       = j*dzeta
         dum        = (zeta-1d0)*(zeta-1d0)/c
         alpha      = 2d0*(zeta-1d0)*dum/c
         beta       = dum*dum
         bz         = beta/dzeta2
         dfdy(5,i-2) = bz-alpha/(2d0*dzeta)
         dfdy(3,i)   = -2d0*bz-k*y(i+1)
         dfdy(1,i+2) = bz+alpha/(2d0*dzeta)
         dfdy(2,i+1) = -k*y(i)
         i          = 2*j
         dfdy(4,i-1) = -k*y(i)
         dfdy(3,i)   = -k*y(i-1)
   30 continue

      dfdy(3,2*N-1) = -k*y(2*N)
      dfdy(2,2*N)   = -k*y(2*N-1)
      dfdy(4,2*N-1) = -k*y(2*N)
      dfdy(3,2*N)   = -k*y(2*N-1)

      return
      end subroutine medakzo_jeval_for_blk_dble


      subroutine medakzo_fcn_blk_dble(n,caller_id,nvar,nz,x,h,y,f,lrpar,rpar,lipar,ipar,ierr)
         use const_def, only: dp
         integer, intent(in) :: n, caller_id, nvar, nz, lrpar, lipar
         real(dp), intent(in) :: x,h
         real(dp), intent(inout), pointer :: y(:) ! (n)
         real(dp), intent(inout), pointer :: f(:) ! (n) ! dy/dx
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         real(dp), target :: yprime_ary(n)
         real(dp), pointer :: yprime(:)
         ierr = 0
         ipar(i_nfcn) = ipar(i_nfcn) + 1
         yprime => yprime_ary
         call medakzo_feval_for_blk_dble(n,x,y,yprime,f,ierr,rpar,ipar)
      end subroutine medakzo_fcn_blk_dble


      subroutine medakzo_derivs(n, x, h, vars, dvars_dx, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: n, lrpar, lipar
         real(dp), intent(in) :: x,h
         real(dp), intent(inout) :: vars(:) ! (n)
         real(dp), intent(inout) :: dvars_dx(:) ! (n)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         real(dp) :: yprime(n)
         integer, intent(out) :: ierr
         include 'formats'
         ierr = 0
         ipar(i_nfcn) = ipar(i_nfcn) + 1
         call medakzo_feval(n,x,vars,yprime,dvars_dx,ierr,rpar,ipar)
      end subroutine medakzo_derivs


      subroutine medakzo_jac_blk_dble(n,caller_id,nvar,nz,x,h,y,f,lblk1,dblk1,ublk1,lrpar,rpar,lipar,ipar,ierr)
         use const_def,only: dp
         integer,intent(in) :: n,caller_id,nvar,nz,lrpar,lipar
         real(dp),intent(in) :: x,h
         real(dp),intent(inout), pointer :: y(:) ! (n)
         real(dp),intent(inout), pointer :: f(:) ! (n) ! dy/dx
         real(dp),dimension(:),pointer,intent(inout) :: lblk1,dblk1,ublk1 ! =(nvar,nvar,nz)
         integer,intent(inout),pointer :: ipar(:) ! (lipar)
         real(dp),intent(inout),pointer :: rpar(:) ! (lrpar)
         integer,intent(out) :: ierr
         
         real(dp),dimension(:,:,:),pointer :: lblk,dblk,ublk ! =(nvar,nvar,nz)
         integer, parameter :: ld_dfdy = 5 ! for medakzo
         real(dp), target :: dfdy1(ld_dfdy*n)
         real(dp), pointer :: dfdy(:,:)
         !real(dp), pointer :: banded(:,:,:)
         integer :: i, k         
         ierr = 0
         dfdy(1:ld_dfdy,1:n) => dfdy1(1:ld_dfdy*n)

         ierr = 0
         ipar(i_njac) = ipar(i_njac) + 1
         call medakzo_jeval_for_blk_dble(ld_dfdy,n,x,y,f,dfdy,ierr,rpar,ipar)
         if (ierr == 0) call medakzo_fcn_blk_dble(n,caller_id,nvar,nz,x,h,y,f,lrpar,rpar,lipar,ipar,ierr)         
         
         !banded(1:ld_dfdy,1:nvar,1:nz) => dfdy1(1:ld_dfdy*n)
         lblk(1:nvar,1:nvar,1:nz) => lblk1(1:nvar*nvar*nz)
         dblk(1:nvar,1:nvar,1:nz) => dblk1(1:nvar*nvar*nz)
         ublk(1:nvar,1:nvar,1:nz) => ublk1(1:nvar*nvar*nz)
         
         ! convert from banded to block tridiagonal
         ! lblk(:,:,1) is not used; ublk(:,:,nz) is not used.
         k = 1         
         dblk(1,1,k) = dfdy(3,1) ! partial of f(1,k) wrt var(1,k)    dfdy(3,i)
         dblk(1,2,k) = dfdy(2,2) ! partial of f(1,k) wrt var(2,k)    dfdy(2,i+1)
         dblk(2,1,k) = dfdy(4,1) ! partial of f(2,k) wrt var(1,k)    dfdy(4,i)
         dblk(2,2,k) = dfdy(3,2) ! partial of f(2,k) wrt var(2,k)    dfdy(3,i+1)
         ublk(1,1,k) = dfdy(1,3) ! partial of f(1,k) wrt var(1,k+1)  dfdy(1,i+2)
         
!dfdy(1,i+2) partial of f(1,k) wrt var(1,k+1)
!dfdy(2,i+1) partial of f(1,k) wrt var(2,k)
!dfdy(3,i)   partial of f(1,k) wrt var(1,k)
!dfdy(3,i+1) partial of f(2,k) wrt var(2,k)
!dfdy(4,i)   partial of f(2,k) wrt var(1,k)
!dfdy(5,i-2) partial of f(1,k) wrt var(1,k-1)         

         do k=2,nz-1
            i = 2*k-1            
            ! set lblk
            lblk(1,1,k) = dfdy(5,i-2) ! partial of f(1,k) wrt var(1,k-1)
            lblk(1,2,k) = 0 ! partial of f(1,k) wrt var(2,k-1)
            lblk(2,1,k) = 0 ! partial of f(2,k) wrt var(1,k-1)
            lblk(2,2,k) = 0 ! partial of f(2,k) wrt var(2,k-1)
            ! set dblk
            dblk(1,1,k) = dfdy(3,i)   ! partial of f(1,k) wrt var(1,k)  dfdy(3,i)
            dblk(1,2,k) = dfdy(2,i+1) ! partial of f(1,k) wrt var(2,k)  dfdy(2,i+1)
            dblk(2,1,k) = dfdy(4,i)   ! partial of f(2,k) wrt var(1,k)  dfdy(4,i)
            dblk(2,2,k) = dfdy(3,i+1) ! partial of f(2,k) wrt var(2,k)  dfdy(3,i+1)
            ! set ublk
            ublk(1,1,k) = dfdy(1,i+2) ! partial of f(1,k) wrt var(1,k+1)   dfdy(1,i+2)
            ublk(2,1,k) = 0 ! partial of f(2,k) wrt var(1,k+1)
            ublk(1,2,k) = 0 ! partial of f(1,k) wrt var(2,k+1)
            ublk(2,2,k) = 0 ! partial of f(2,k) wrt var(2,k+1)
         end do   
         
         k = nz
         i = 2*k-1
         dblk(1,1,k) = dfdy(3,i)   ! partial of f(1,k) wrt var(1,k)
         dblk(1,2,k) = dfdy(2,i+1) ! partial of f(1,k) wrt var(2,k)
         dblk(2,1,k) = dfdy(4,i)   ! partial of f(2,k) wrt var(1,k)
         dblk(2,2,k) = dfdy(3,i+1) ! partial of f(2,k) wrt var(2,k)
         
      end subroutine medakzo_jac_blk_dble


      subroutine medakzo_jacob(n,x,h,y,f,dfdy,ld_dfdy,lrpar,rpar,lipar,ipar,ierr)
         integer, intent(in) :: n, ld_dfdy, lrpar, lipar
         real(dp), intent(in) :: x, h
         real(dp), intent(inout) :: y(:)
         real(dp), intent(inout) :: f(:), dfdy(:,:)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         real(dp) :: yprime(n)
         integer, intent(out) :: ierr
         integer :: nz, i, j
         include 'formats'
         ierr = 0
         ipar(i_njac) = ipar(i_njac) + 1
         call medakzo_jeval(ld_dfdy,n,x,y,yprime,dfdy,ierr,rpar,ipar)
         if (ierr == 0) call medakzo_derivs(n, x, h, y, f, lrpar,rpar,lipar,ipar, ierr)
         !write(*,*)
         !write(*,2)'medakzo_jacob', ipar(i_njac), x, y(1), dfdy(3,1:2)
      end subroutine medakzo_jacob


      subroutine medakzo_sjac(n,x,h,y,f,nzmax,ia,ja,values,lrpar,rpar,lipar,ipar,ierr)  
         ! sparse jacobian. format either compressed row or compressed column.
         use mtx_lib,only:band_to_row_sparse_with_diag,band_to_col_sparse_with_diag,mtx_rcond_banded
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
         call medakzo_jacob(n,x,h,y,f,dfdy,ld_dfdy,lrpar,rpar,lipar,ipar,ierr)
         if (ierr /= 0) return
         if (ipar(ipar_sparse_format) == 0) then
            call band_to_row_sparse_with_diag(n,mljac,mujac,dfdy,ld_dfdy,nzmax,nz,ia,ja,values,ierr)
         else
            call band_to_col_sparse_with_diag(n,mljac,mujac,dfdy,ld_dfdy,nzmax,nz,ia,ja,values,ierr)
         end if
      end subroutine medakzo_sjac
      
      
      subroutine do_test_medakzo(which_solver,which_decsol,numerical_jacobian,show_all,quiet)
         use test_support,only:show_results,show_statistics,check_results
         use test_int_support,only:do_test_stiff_int
         integer, intent(in) :: which_solver,which_decsol
         logical, intent(in) :: numerical_jacobian,show_all,quiet

         integer, parameter :: nvar = 2, nz = 200
         integer, parameter :: n = nvar*nz ! the number of variables in the "medakzo" system of ODEs
         real(dp), target :: y_ary(n), yprime(n), yexact(n)
         real(dp), pointer :: y(:)
         integer, parameter :: lrpar = 1, lipar = 3, iout=1
         logical :: consis
         integer, parameter :: ndisc = 1, n_soln=11
         real(dp) :: result(n_soln), soln(n_soln), h0, t(0:ndisc+1), atol(1), rtol(1)
         integer :: i, j, k, matrix_type_spec, ierr, imas, mlmas, mumas, m1, m2, itol, nstep
         real(dp), target :: rpar_ary(lrpar) 
         integer, target :: ipar_ary(lipar)
         real(dp), pointer :: rpar(:)
         integer, pointer :: ipar(:)
         integer :: caller_id, nvar_blk_dble, nz_blk_dble
         real(dp), dimension(:), pointer :: lblk, dblk, ublk ! =(nvar,nvar,nz)
         real(dp), dimension(:), pointer :: uf_lblk, uf_dblk, uf_ublk ! =(nvar,nvar,nz)
         logical, parameter :: dbg = .false.

         include 'formats'
                  
         rpar => rpar_ary
         ipar => ipar_ary
         y => y_ary
            
         if (.not. quiet) write(*,*) 'medakzo'

         nullify(lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk)
         caller_id = 0
         nvar_blk_dble = 0
         nz_blk_dble = 0

         t(0)   = 0d0
         if (dbg) then
            t(1)   = 0.05d0
            t(2)   = 0.20d0
         else
            t(1)   = 5d0
            t(2)   = 20d0
         end if
         
         itol = 0 ! scalar tolerances
         rtol = 1d-6
         atol = 1d-6
         h0 = 1d-9 ! initial step size
         
         matrix_type_spec = banded_matrix_type
         mljac = 2
         mujac = 2

         imas = 0
         mlmas = 0
         mumas = 0        
         
         m1 = 0
         m2 = 0     
         
         call medakzo_init(n,y,yprime,consis)
         nstep=0        

         if (nvar_blk_dble == 0) then
            call do_test_stiff_int(which_solver,which_decsol,numerical_jacobian, &
               medakzo_derivs,medakzo_jacob,medakzo_sjac,medakzo_solout,iout, &
               null_fcn_blk_dble,null_jac_blk_dble, &
               caller_id,nvar_blk_dble,nz_blk_dble,lblk,dblk,ublk,uf_lblk,uf_dblk,uf_ublk, &
               n,ndisc,mljac,mujac,matrix_type_spec,null_mas,imas,mlmas,mumas,m1,m2, &
               t,rtol,atol,itol,h0,y,nstep,lrpar,rpar,lipar,ipar,quiet,ierr)
         else
            call do_test_stiff_int(which_solver,which_decsol,numerical_jacobian, &
               null_fcn,null_jac,null_sjac,medakzo_solout,iout, &
               medakzo_fcn_blk_dble,medakzo_jac_blk_dble, &
               caller_id,nvar_blk_dble,nz_blk_dble,lblk,dblk,ublk,uf_lblk,uf_dblk,uf_ublk, &
               n,ndisc,mljac,mujac,matrix_type_spec,null_mas,imas,mlmas,mumas,m1,m2, &
               t,rtol,atol,itol,h0,y,nstep,lrpar,rpar,lipar,ipar,quiet,ierr)
         end if
         if (ierr /= 0) then
            write(*,*) 'test_medakzo ierr', ierr
            call mesa_error(__FILE__,__LINE__)
         end if
         
         call medakzo_solut(n,0d0,yexact)
         j = 1
         do k = 1, n/2, max(1,(n/2-1)/11)
            if (j > n_soln) exit
            result(j) = y(1+2*(k-1))
            soln(j) = yexact(1+2*(k-1))
            j = j+1
         end do

         if (.not. dbg) then
            call check_results(n,y,yexact,rtol(1)*50,ierr)
            if (ierr /= 0) then
               write(*,*) 'check results ierr', ierr
               !call mesa_error(__FILE__,__LINE__) ! do_test_medakzo
            end if
         end if
         
         if (quiet) return
         
         call show_results(n_soln,result,soln,show_all)
         call show_statistics(ipar(i_nfcn),ipar(i_njac),nstep,show_all)

      end subroutine do_test_medakzo


      subroutine medakzo_solout(nr,xold,x,n,y,rwork,iwork,interp_y,lrpar,rpar,lipar,ipar,irtrn)
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
         real(dp) :: xout, y1, y2
         integer :: ierr
         include 'formats'
         !if (mod(nr,10) == 0) write(*,2) 'step', nr, x, y(1:2)
         !if (nr >= 100) stop
         ierr = 0
         irtrn = 0
      end subroutine medakzo_solout
            
      
      end module test_medakzo
