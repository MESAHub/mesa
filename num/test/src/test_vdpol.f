      module test_vdpol
      use num_def
      use num_lib
      use mtx_lib
      use mtx_def
      use test_int_support, only: i_nfcn, i_njac
      use utils_lib, only: mesa_error

      implicit none
      
      logical, parameter :: dbg = .false.
      
      integer :: cnt = 0
      

      contains


      subroutine vdpol_fcn_blk_dble(n,caller_id,nvar,nz,x,h,y,f,lrpar,rpar,lipar,ipar,ierr)
         use const_def, only: dp
         integer, intent(in) :: n, caller_id, nvar, nz, lrpar, lipar
         real(dp), intent(in) :: x, h
         real(dp), intent(inout), pointer :: y(:) ! (n)
         real(dp), intent(inout), pointer :: f(:) ! (n) ! dy/dx
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         include 'formats'
         call vdpol_derivs(n, x, h, y, f, lrpar, rpar, lipar, ipar, ierr)
         !write(*,2) 'vdpol_fcn_blk_dble', ipar(i_nfcn), x
      end subroutine vdpol_fcn_blk_dble


      subroutine vdpol_derivs(n, x, h, vars, dvars_dx, lrpar,rpar,lipar,ipar, ierr)
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
         call vdpol_feval(n,x,vars,yprime,dvars_dx,ierr,rpar,ipar)

         if (dbg) then
            cnt = cnt+1
            write(*,*) 'func cnt', cnt
            write(*,*) 'func n', n
            write(*,*) 'func x', x
            write(*,*) 'func y(1)', vars(1)
            write(*,*) 'func y(2)', vars(2)
            write(*,*) 'func f(1)', dvars_dx(1)
            write(*,*) 'func f(2)', dvars_dx(2)
            write(*,*)
            if (cnt == 5) stop 'vdpol_derivs'
         end if

      end subroutine vdpol_derivs


      subroutine vdpol_jac_blk_dble(n,caller_id,nvar,nz,x,h,y,f,lblk1,dblk1,ublk1,lrpar,rpar,lipar,ipar,ierr)
         use const_def,only: dp
         integer,intent(in) :: n,caller_id,nvar,nz,lrpar,lipar
         real(dp),intent(in) :: x, h
         real(dp),intent(inout), pointer :: y(:) ! (n)
         real(dp),intent(inout), pointer :: f(:) ! (n) ! dy/dx
         real(dp),dimension(:),pointer,intent(inout) :: lblk1,dblk1,ublk1 ! =(nvar,nvar,nz)
         integer,intent(inout),pointer :: ipar(:) ! (lipar)
         real(dp),intent(inout),pointer :: rpar(:) ! (lrpar)
         integer,intent(out) :: ierr
         
         real(dp),dimension(:,:,:),pointer :: lblk,dblk,ublk ! =(nvar,nvar,nz)
         integer, parameter :: ld_dfdy = 2 ! for vdpol
         real(dp), target :: dfdy1(ld_dfdy*n)
         real(dp), pointer :: dfdy(:,:)
         integer :: i, k 
         include 'formats'        
         ierr = 0
         dfdy(1:ld_dfdy,1:n) => dfdy1(1:ld_dfdy*n)
         call vdpol_jacob(n,x,h,y,f,dfdy,ld_dfdy,lrpar,rpar,lipar,ipar,ierr)
         !write(*,2) 'vdpol_jac_blk_dble', ipar(i_nfcn), x
         if (ierr /= 0) return
         lblk(1:nvar,1:nvar,1:nz) => lblk1(1:nvar*nvar*nz)
         dblk(1:nvar,1:nvar,1:nz) => dblk1(1:nvar*nvar*nz)
         ublk(1:nvar,1:nvar,1:nz) => ublk1(1:nvar*nvar*nz)
         ! convert from square to block tridiagonal
         dblk(:,:,1) = dfdy(:,:)
         ublk(:,:,1) = 0
         lblk(:,:,1) = 0
      end subroutine vdpol_jac_blk_dble


      subroutine vdpol_jacob(n,x,h,y,f,dfdy,ld_dfdy,lrpar,rpar,lipar,ipar,ierr)
         integer, intent(in) :: n, ld_dfdy, lrpar, lipar
         real(dp), intent(in) :: x, h
         real(dp), intent(inout) :: y(:)
         real(dp), intent(inout) :: f(:), dfdy(:,:)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         real(dp) :: yprime(n)
         integer, intent(out) :: ierr
         integer :: nz, i, j
         ierr = 0
         ipar(i_njac) = ipar(i_njac) + 1
         call vdpol_jeval(ld_dfdy,n,x,y,yprime,dfdy,ierr,rpar,ipar)
         if (ierr == 0) call vdpol_derivs(n, x, h, y, f, lrpar,rpar,lipar,ipar, ierr)

      
         if (.false.) then
            write(*,*) 'jac_fcn y(1)', y(1)
            write(*,*) 'jac_fcn y(2)', y(2)
            write(*,*) 'jac_fcn dfdy(1,1)', dfdy(1,1)
            write(*,*) 'jac_fcn dfdy(1,2)', dfdy(1,2)
            write(*,*) 'jac_fcn dfdy(2,1)', dfdy(2,1)
            write(*,*) 'jac_fcn dfdy(2,2)', dfdy(2,2)
            write(*,*)
         end if

      end subroutine vdpol_jacob


      subroutine vdpol_sjac(n,x,h,y,f,nzmax,ia,ja,values,lrpar,rpar,lipar,ipar,ierr)  
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
         call vdpol_jacob(n,x,h,y,f,dfdy,ld_dfdy,lrpar,rpar,lipar,ipar,ierr)
         if (ierr /= 0) return
         if (ipar(ipar_sparse_format) == 0) then
            call dense_to_row_sparse_with_diag(n,n,dfdy,nzmax,nz,ia,ja,values,ierr)
         else
            call dense_to_col_sparse_with_diag(n,n,dfdy,nzmax,nz,ia,ja,values,ierr)
         end if
      end subroutine vdpol_sjac


      subroutine vdpol_solout(nr,xold,x,n,y,rwork,iwork,interp_y,lrpar,rpar,lipar,ipar,irtrn)
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
         
         if (dbg .and. nr > 450) stop
         
         
         ierr = 0
         irtrn = 0
         xout = rpar(1)
         if (nr.eq.1) then
            write (6,99) x,y(1),y(2),nr-1
            xout=0.2d0
         else
            do
               if (x.ge.xout) then
                  y1 = interp_y(1,xout,rwork,iwork,ierr)
                  if (ierr /= 0) exit
                  y2 = interp_y(2,xout,rwork,iwork,ierr)
                  write (6,99) xout,y1,y2,nr-1
                  if (ierr /= 0) exit
                  xout=xout+0.2d0
                  cycle
               end if
               exit
            end do
         end if
         if (ierr /= 0) then
            write(*,*) 'problem with interp_y in vdpol_solout'
            irtrn = -1
         end if
         rpar(1) = xout
  99     format(1x,'x =',f5.2,'    y =',2e18.10,'    nstep =',i8)
      end subroutine vdpol_solout
      
      
      subroutine do_test_vdpol(which_solver,which_decsol,numerical_jacobian,show_all,quiet)
         use test_support,only:show_results,show_statistics,check_results
         use test_int_support,only:do_test_stiff_int
         integer, intent(in) :: which_solver,which_decsol
         logical, intent(in) :: numerical_jacobian,show_all,quiet

         integer, parameter :: nvar = 2, nz = 1
         integer, parameter :: n = nvar*nz ! the number of variables in the "vdpol" system of ODEs
         real(dp), target :: y_ary(n), yprime(n), yexact(n)
         real(dp), pointer :: y(:)
         integer, parameter :: lrpar = 1, lipar = 3
         logical :: consis
         integer, parameter :: ndisc = 0
         real(dp) :: h0, t(0:ndisc+1), atol(1), rtol(1)
         integer :: i, mujac, mljac, matrix_type_spec, ierr, imas, mlmas, mumas, m1, m2, itol, iout, nstep
         real(dp), target :: rpar_ary(lrpar) 
         integer, target :: ipar_ary(lipar)
         real(dp), pointer :: rpar(:)
         integer, pointer :: ipar(:)
         integer :: caller_id, nvar_blk_dble, nz_blk_dble
         real(dp), dimension(:), pointer :: lblk, dblk, ublk ! =(nvar,nvar,nz)
         real(dp), dimension(:), pointer :: uf_lblk, uf_dblk, uf_ublk ! =(nvar,nvar,nz)
                  
         rpar => rpar_ary
         ipar => ipar_ary
         y => y_ary
            
         if (.not. quiet) write(*,*) 'vdpol'

         nullify(lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk)
         caller_id = 0
         nvar_blk_dble = 0
         nz_blk_dble = 0
         
         t(0)   = 0
         t(1)   = 2d0
         
         itol = 0 ! scalar tolerances
         rtol(1) = 1d-8
         atol(1) = 1d-8
         h0 = 1d-10 ! initial step size
         
         rtol(1) = 1d-6
         atol(1) = 1d-6
         h0 = 1d-8 ! initial step size
         
         rtol(1) = 1d-4
         atol(1) = 1d-4
         h0 = 1d-4 ! initial step size
         
         mljac = n ! square matrix
         mujac = n
         matrix_type_spec = square_matrix_type

         imas = 0
         mlmas = 0
         mumas = 0        
         
         m1 = 0
         m2 = 0     
         
         if (show_all) then
            iout = 1
         else
            iout = 0
         end if
         
         ipar = 0
         
         call vdpol_init(n,y,yprime,consis)
         nstep=0  

         if (nvar_blk_dble == 0) then
            call do_test_stiff_int(which_solver,which_decsol,numerical_jacobian, &
               vdpol_derivs,vdpol_jacob,vdpol_sjac,vdpol_solout,iout, &
               null_fcn_blk_dble,null_jac_blk_dble, &
               caller_id,nvar_blk_dble,nz_blk_dble,lblk,dblk,ublk,uf_lblk,uf_dblk,uf_ublk, &
               n,ndisc,mljac,mujac,matrix_type_spec,null_mas,imas,mlmas,mumas,m1,m2, &
               t,rtol,atol,itol,h0,y,nstep,lrpar,rpar,lipar,ipar,quiet,ierr)
         else
            call do_test_stiff_int(which_solver,which_decsol,numerical_jacobian, &
               null_fcn,null_jac,null_sjac,vdpol_solout,iout, &
               vdpol_fcn_blk_dble,vdpol_jac_blk_dble, &
               caller_id,nvar_blk_dble,nz_blk_dble,lblk,dblk,ublk,uf_lblk,uf_dblk,uf_ublk, &
               n,ndisc,mljac,mujac,matrix_type_spec,null_mas,imas,mlmas,mumas,m1,m2, &
               t,rtol,atol,itol,h0,y,nstep,lrpar,rpar,lipar,ipar,quiet,ierr)
         end if
         if (ierr /= 0) then
            write(*,*) 'test_vdpol ierr', ierr
            call mesa_error(__FILE__,__LINE__)
         end if
         
         call vdpol_solut(n,0d0,yexact)
         !call check_results(n,y,yexact,rtol(1)*2,ierr)
         if (ierr /= 0) then
            write(*,*) 'check results ierr', ierr
            call mesa_error(__FILE__,__LINE__) ! do_test_vdpol
         end if
         
         if (quiet) return
         
         call show_results(n,y,yexact,show_all)
         call show_statistics(ipar(i_nfcn),ipar(i_njac),nstep,show_all)

      end subroutine do_test_vdpol
            
      
      end module test_vdpol
