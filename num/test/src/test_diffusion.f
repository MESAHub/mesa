      module test_diffusion
      use num_def
      use num_lib
      use math_lib
      use test_int_support, only: i_nfcn, i_njac
      use utils_lib, only: mesa_error

      implicit none
      
      integer :: mljac, mujac, nstep
      
      integer, parameter :: nz=48
      integer, parameter :: diff_mujac=1, diff_mljac=1, diff_ldjac=diff_mujac+diff_mljac+1

      real(dp) :: y0(nz), yexact(nz), dfdy2(4,nz), rcond, work(3*nz)
      integer :: ipiv(nz), iwork(nz), info

      real(dp) :: sig_dm(nz)

      contains


      subroutine diffusion_op(n,x,h,y,f,dfdy,ld_dfdy,lrpar,rpar,lipar,ipar,ierr)
         integer, intent(in) :: n, ld_dfdy, lrpar, lipar
         real(dp), intent(in) :: x, h
         real(dp), intent(inout) :: y(n)
         real(dp), intent(inout) :: f(n), dfdy(ld_dfdy,n)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         
         real(dp) :: sig1, sig2
         integer :: i, j
         ierr = 0

         f = 0; dfdy=0
         
         sig2 = 0
         do i=1,n
            sig1 = sig2
            if (i < n) then
               sig2 = sig_dm(i+1)
            else
               sig2 = 0
            end if
            if (ld_dfdy > 0) then
               if (i > 1) then
                  dfdy(3,i-1) = sig1
               end if
               dfdy(2,i) = -(sig1+sig2);
               if (i < n) then
                  dfdy(1,i+1) = sig2;
               end if
            end if
            if (i == n) then
               f(i) = sig1*(y(i-1)-y(i));
            else if (i == 1) then
               f(i) = -sig2*(y(i)-y(i+1));
            else
               f(i) = sig1*(y(i-1)-y(i)) - sig2*(y(i)-y(i+1));
            end if
         end do

      end subroutine diffusion_op
      

      subroutine diffusion_derivs(n, x, h, y, f, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: n, lrpar, lipar
         real(dp), intent(in) :: x, h
         real(dp), intent(inout) :: y(:) ! (n)
         real(dp), intent(inout) :: f(:) ! (n)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         
         real(dp) :: dfdy(0,n)
         
         ierr = 0
         ipar(i_nfcn) = ipar(i_nfcn) + 1
         call diffusion_op(n,x,h,y,f,dfdy,0,lrpar,rpar,lipar,ipar,ierr)
      end subroutine diffusion_derivs


      subroutine diffusion_jacob(n,x,h,y,f,dfdy,ld_dfdy,lrpar,rpar,lipar,ipar,ierr)
         use mtx_lib,only:mtx_rcond_banded
         integer, intent(in) :: n, ld_dfdy, lrpar, lipar
         real(dp), intent(in) :: x, h
         real(dp), intent(inout) :: y(:)
         real(dp), intent(inout) :: f(:), dfdy(:,:)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         real(dp) :: yprime(n)
         integer, intent(out) :: ierr
         include 'formats'
         ierr = 0
         ipar(i_njac) = ipar(i_njac) + 1
         call diffusion_op(n,x,h,y,f,dfdy,ld_dfdy,lrpar,rpar,lipar,ipar,ierr)
         
         return
         
         dfdy2(2,1:n) = dfdy(1,1:n)
         dfdy2(3,1:n) = dfdy(2,1:n) - 1
         dfdy2(4,1:n) = dfdy(3,1:n)
         
         call mtx_rcond_banded('N', n, n, 1, 1, dfdy2, 4, ipiv, rcond, work, iwork, info)
         write(*,2) 'diffusion_jacob rcond', info, x, safe_log10(rcond)
         
      end subroutine diffusion_jacob


      subroutine diffusion_sjac(n,x,h,y,f,nzmax,ia,ja,values,lrpar,rpar,lipar,ipar,ierr)  
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
         call diffusion_jacob(n,x,h,y,f,dfdy,ld_dfdy,lrpar,rpar,lipar,ipar,ierr)
         if (ierr /= 0) return
         if (ipar(ipar_sparse_format) == 0) then
            call band_to_row_sparse_with_diag(n,mljac,mujac,dfdy,ld_dfdy,nzmax,nz,ia,ja,values,ierr)
         else
            call band_to_col_sparse_with_diag(n,mljac,mujac,dfdy,ld_dfdy,nzmax,nz,ia,ja,values,ierr)
         end if
      end subroutine diffusion_sjac


      subroutine diffusion_solout(nr,xold,x,n,y,rwork,iwork,interp_y,lrpar,rpar,lipar,ipar,irtrn)
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
      end subroutine diffusion_solout
      
      
      subroutine do_test_diffusion(which_solver,which_decsol,numerical_jacobian,show_all,quiet)
         use test_support,only:show_results,show_statistics,check_results
         use test_int_support,only:do_test_stiff_int
         integer, intent(in) :: which_solver,which_decsol
         logical, intent(in) :: numerical_jacobian,show_all,quiet

         real(dp), parameter :: sig = 1d-2

         real(dp), parameter :: ystart = 1d0, tend = 1d4
         integer, parameter :: n = nz
         integer, parameter :: lrpar = 1, lipar = 3, iout=1
         integer, parameter :: ndisc = 0, n_soln=10
         real(dp), target :: y_ary(n)
         real(dp), pointer :: y(:)
         real(dp) :: result(n_soln), soln(n_soln), h0, atol(1), rtol(1), t(0:ndisc+1)
         integer :: i, j, k, matrix_type_spec, ierr, imas, mlmas, mumas, m1, m2, itol
         real(dp), target :: rpar_ary(lrpar) 
         integer, target :: ipar_ary(lipar)
         real(dp), pointer :: rpar(:)
         integer, pointer :: ipar(:)
         integer :: caller_id, nvar_blk_dble, nz_blk_dble
         real(dp), dimension(:), pointer :: lblk, dblk, ublk ! =(nvar,nvar,nz)
         real(dp), dimension(:), pointer :: uf_lblk, uf_dblk, uf_ublk ! =(nvar,nvar,nz)
         
         nullify(lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk)
         caller_id = 0
         nvar_blk_dble = 0
         nz_blk_dble = 0
                  
         rpar => rpar_ary
         ipar => ipar_ary
         y => y_ary
            
         if (.not. quiet) write(*,*) 'diffusion'

         t(0)   = 0d0
         t(1)   = tend
         
         itol = 0 ! scalar tolerances
         rtol = 1d-6
         atol = 1d-6
         h0 = atol(1)*1d-1 ! initial step size
         
         matrix_type_spec = banded_matrix_type
         mljac = diff_mljac
         mujac = diff_mujac

         imas = 0
         mlmas = 0
         mumas = 0        
         
         m1 = 0
         m2 = 0     
         
         k=nz/2
         y(1:k) = 0
         y(k+1:nz) = ystart
         
         y0 = y
         
         do k=1,nz
            if (k == 1) then
               sig_dm(k) = 0;
            else if (5*k <= n) then
               sig_dm(k) = 0;
            else if (5*k >= 4*n) then
               sig_dm(k) = 0;
            else
               sig_dm(k) = sig;
            end if
         end do
         
         nstep=0        

         call do_test_stiff_int(which_solver,which_decsol,numerical_jacobian, &
               diffusion_derivs,diffusion_jacob,diffusion_sjac,diffusion_solout,iout, &
               null_fcn_blk_dble,null_jac_blk_dble, &
               caller_id,nvar_blk_dble,nz_blk_dble,lblk,dblk,ublk,uf_lblk,uf_dblk,uf_ublk, &
               n,ndisc,mljac,mujac,matrix_type_spec,null_mas,imas,mlmas,mumas,m1,m2, &
               t,rtol,atol,itol,h0,y,nstep,lrpar,rpar,lipar,ipar,quiet,ierr)
         if (ierr /= 0) then
            write(*,*) 'test_diffusion ierr', ierr
            call mesa_error(__FILE__,__LINE__)
         end if
         
         !call write_diffusion_results
         
         call set_yexact
         i=0
         do k=10,nz-10,(nz-10)/12
            i=i+1
            if (i > n_soln) exit
            soln(i) = yexact(k)
            result(i) = y(k)
         end do
         
         call show_results(n_soln,result,soln,show_all)
         call show_statistics(ipar(i_nfcn),ipar(i_njac),nstep,show_all)
         
         contains
         
         subroutine set_yexact
            ! for nz=48, sig = 1d-2, ystart = 1d0, tend = 1d4
            yexact( 1)=    0.0000000000000000D+00
            yexact( 2)=    0.0000000000000000D+00
            yexact( 3)=    0.0000000000000000D+00
            yexact( 4)=    0.0000000000000000D+00
            yexact( 5)=    0.0000000000000000D+00
            yexact( 6)=    0.0000000000000000D+00
            yexact( 7)=    0.0000000000000000D+00
            yexact( 8)=    0.0000000000000000D+00
            yexact( 9)=    2.5602881146875722D-01
            yexact(10)=    2.5830832258765279D-01
            yexact(11)=    2.6284365888436573D-01
            yexact(12)=    2.6958764648712957D-01
            yexact(13)=    2.7847002063795223D-01
            yexact(14)=    2.8939802361837114D-01
            yexact(15)=    3.0225720584222654D-01
            yexact(16)=    3.1691243229910382D-01
            yexact(17)=    3.3320909577496399D-01
            yexact(18)=    3.5097453679942375D-01
            yexact(19)=    3.7001966806159553D-01
            yexact(20)=    3.9014079814263225D-01
            yexact(21)=    4.1112164592868455D-01
            yexact(22)=    4.3273553313243934D-01
            yexact(23)=    4.5474773813802244D-01
            yexact(24)=    4.7691799008736224D-01
            yexact(25)=    4.9900307794850141D-01
            yexact(26)=    5.2075954544461744D-01
            yexact(27)=    5.4194643935567555D-01
            yexact(28)=    5.6232807598367041D-01
            yexact(29)=    5.8167678861286254D-01
            yexact(30)=    5.9977561767405640D-01
            yexact(31)=    6.1642090507187908D-01
            yexact(32)=    6.3142475475257254D-01
            yexact(33)=    6.4461732303943897D-01
            yexact(34)=    6.5584890447885325D-01
            yexact(35)=    6.6499178183713092D-01
            yexact(36)=    6.7194181237134964D-01
            yexact(37)=    6.7661972646501933D-01
            yexact(38)=    6.7897211907369082D-01
            yexact(39)=    1.0000000000000000D+00
            yexact(40)=    1.0000000000000000D+00
            yexact(41)=    1.0000000000000000D+00
            yexact(42)=    1.0000000000000000D+00
            yexact(43)=    1.0000000000000000D+00
            yexact(44)=    1.0000000000000000D+00
            yexact(45)=    1.0000000000000000D+00
            yexact(46)=    1.0000000000000000D+00
            yexact(47)=    1.0000000000000000D+00
            yexact(48)=    1.0000000000000000D+00
         end subroutine set_yexact
         
         subroutine write_diffusion_results
            use utils_lib, only: mkdir
            use const_def
            character (len=100) :: filename,dir
            integer :: k, ierr, iounit
            dir='plot_data'
            call mkdir(dir)
            filename = trim(dir)//'/'//'diffusion.data'
            ierr = 0
            open(newunit=iounit, file=trim(filename), action='write', status='replace', iostat=ierr)
            if (ierr == 0) then
               write(*,*) 'write burn results to ' // trim(filename)
               write(iounit,'(a)') 'burn log'
               write(iounit,'(99(a,1x))') 'y', 'y0', 'sigma'
               do k=1,nz
                  write(iounit,'(99e24.10)') y(k), y0(k), sig_dm(k)
               end do
               close(iounit)
            else
               write(*,*) 'failed to open internals file ' // trim(filename)
            end if
         end subroutine write_diffusion_results

      end subroutine do_test_diffusion
      
      
            
      
      end module test_diffusion
