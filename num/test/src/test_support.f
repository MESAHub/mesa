      module test_support
      use num_def
      use num_lib
      use math_lib
      use const_def, only: dp, arg_not_provided
      use utils_lib, only: mesa_error
      
      implicit none

      contains

      
      subroutine show_results(nv,y,expect,show_all)
         integer, intent(in) :: nv
         real(dp), dimension(nv), intent(in) :: y, expect
         logical, intent(in) :: show_all
         integer :: i
         include 'formats'
         if (show_all) then
            do i=1,nv
               write(*,2) 'calculated', i, y(i)
            end do
            do i=1,nv
               write(*,2) 'reference', i, expect(i)
               write(*,2) 'lg(abs rel diff)', i, &
                  log10(abs(y(i)-expect(i))/max(1d-299,abs(expect(i))))
            end do
         else
            do i=1,nv
               write(*,2) 'calculated', i, y(i)
            end do
            do i=1,nv
               write(*,2) 'reference', i, expect(i)
            end do
         end if
         write(*,*)
      end subroutine show_results
      
      
      subroutine show_statistics(nfcn,njac,nstep,show_all)
         integer, intent(in) :: nfcn,njac,nstep
         logical, intent(in) :: show_all
         if (.not. show_all) return
         write(*,*) 'number of steps         ', nstep
         write(*,*) 'number of function evals', nfcn
         write(*,*) 'number of jacobians     ', njac
         write(*,*) 'functions + jacobians   ', nfcn+njac
         write(*,*)
      end subroutine show_statistics

      
      subroutine check_results(nv,y,expect,tol,ierr)
         integer, intent(in) :: nv
         real(dp), dimension(nv), intent(in) :: y, expect
         real(dp), intent(in) :: tol
         integer, intent(out) :: ierr
         integer :: i
         logical :: okay
         include 'formats.dek'
         okay = .true.
         ierr = 0
         do i=1,nv
            if (abs(expect(i)) < tol) cycle
            if (abs(y(i)-expect(i)) > tol) then
               write(*,2) 'check results result#', i
               write(*,1) 'log10 abs diff', log10(abs(y(i)-expect(i)))
               write(*,1) 'y(i)', y(i)
               write(*,1) 'expect(i)', expect(i)
               write(*,1) 'log10 tol', log10(tol)
               write(*,*)
               ierr = -1
               okay = .false.
            end if
         end do
         if (okay) return
         write(*,*)
         do i=1,nv
            write(*,2) 'y expected', i, y(i), expect(i)
         end do
         write(*,*)
         call mesa_error(__FILE__,__LINE__)
         return
      end subroutine check_results

      real(dp) function f(x,dfdx,lrpar,rpar,lipar,ipar,ierr)
         integer, intent( in ) :: lrpar, lipar
         real(dp), intent( in ) :: x
         real(dp), intent( out ) :: dfdx
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent( out ) :: ierr
         ierr = 0
         f = x-3*sin(1-x)
         dfdx = 1+3*cos(1-x)
      end function f
      
      subroutine test_root_with_brackets
         integer, parameter :: lrpar=0, lipar=0
         real(dp) :: x, dfdx, y
         real(dp) :: x1, x3 ! bounds for x
            ! values of f at x1 and x3 must have opposite sign
            ! return value for safe_root will be bracketed by x1 and x3
         real(dp) :: y1, y3 ! f(x1) and f(x3)
         integer :: imax ! max number of iterations for search
         real(dp) :: epsx, epsy 
         ! stop seaching when x is determined to within epsx
         ! or when abs(f(x)) is less than epsy
         integer :: i, ierr
         real(dp) :: expected_root = 0.74800611d0
         real(dp), target :: rpar_ary(lrpar)
         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:) ! (lipar)
         real(dp), pointer :: rpar(:) ! (lrpar)
         include 'formats'
         ipar => ipar_ary
         rpar => rpar_ary
         ierr = 0
         imax = 100
         x1 = 0
         x3 = 10
         y1 = f(x1,dfdx,lrpar,rpar,lipar,ipar,ierr)
         y3 = f(x3,dfdx,lrpar,rpar,lipar,ipar,ierr)
         epsx = 1d-6
         epsy = 1d-6
         x = safe_root_with_brackets( &
            f,x1,x3,y1,y3,imax,epsx,epsy,lrpar,rpar,lipar,ipar,ierr)
         if (abs(x-expected_root) > 1d-6) call mesa_error(__FILE__,__LINE__)
         write(*,1) 'root', x
      end subroutine test_root_with_brackets
      

      real(dp) function test_f(x,dfdx,lrpar,rpar,lipar,ipar,ierr)
         ! returns with ierr = 0 if was able to evaluate f and df/dx at x
         real(dp), intent(in) :: x
         real(dp), intent(out) :: dfdx
         integer, intent(in) :: lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         test_f = tanh(x) - 0.4621171572600098d0
         dfdx = 1/cosh(x)**2
         ierr = 0       
      end function test_f
      
      
      subroutine test_root2
         real(dp) :: x ! provide starting guess on input
         real(dp) :: x1,x3 ! bounds for x
         real(dp) :: y1,y3 ! f(x1) and f(x3)
         integer, parameter :: imax = 50, lipar = 0, lrpar = 0
         real(dp) :: dx
         real(dp), parameter :: epsx = 1d-10, epsy = 1d-10
         integer :: ierr      
         real(dp), target :: rpar_ary(lrpar)
         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:) ! (lipar)
         real(dp), pointer :: rpar(:) ! (lrpar)
         include 'formats'
         ipar => ipar_ary
         rpar => rpar_ary
         x = -1d0
         dx = 0.1d0
         ierr = 0      
         write(*,*) 'test_root2'         
         call look_for_brackets(x,dx,x1,x3,test_f,y1,y3,imax,lrpar,rpar,lipar,ipar,ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         write(*,1) 'x1', x1
         write(*,1) 'x3', x3         
         write(*,1) 'y1', y1
         write(*,1) 'y3', y3
         x = safe_root_with_brackets( &
            test_f,x1,x3,y1,y3,imax,epsx,epsy,lrpar,rpar,lipar,ipar,ierr)   
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         write(*,1) 'safe_root', x
         write(*,*)         
      end subroutine test_root2
      
      
      subroutine test_root3
         real(dp) :: x ! provide starting guess on input
         real(dp) :: x1, x3 ! bounds for x
         real(dp) :: y1, y3 ! f(x1) and f(x3)
         integer, parameter :: newt_imax = 10, imax = 50, lipar = 0, lrpar = 0
         real(dp) :: dx
         real(dp), parameter :: epsx = 1d-10, epsy = 1d-10
         integer :: ierr      
         real(dp), target :: rpar_ary(lrpar)
         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:) ! (lipar)
         real(dp), pointer :: rpar(:) ! (lrpar)
         include 'formats'
         ipar => ipar_ary
         rpar => rpar_ary
         dx = 0.1d0
         x1 = arg_not_provided
         x3 = arg_not_provided
         y1 = arg_not_provided
         y3 = arg_not_provided
         ierr = 0      
         write(*,*) 'test_root3'         
         x = 0.1d0 ! not too bad initial guess.  newton should find it okay.
         x = safe_root_with_guess( &
            test_f, x, dx, x1, x3, y1, y3, &
            newt_imax, imax, epsx, epsy, lrpar, rpar, lipar, ipar, ierr)   
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         write(*,1) 'first safe_root_with_guess', x
         x = -1d0 ! really bad guess will make it give up on newton
         x = safe_root_with_guess( &
            test_f, x, dx, x1, x3, y1, y3, &
            newt_imax, imax, epsx, epsy, lrpar, rpar, lipar, ipar, ierr)   
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         write(*,1) 'second safe_root_with_guess', x
         write(*,*)         
      end subroutine test_root3


      subroutine van_der_Pol_derivs(n,x,h,y,f,lrpar,rpar,lipar,ipar,ierr)
         integer, intent(in) :: n, lrpar, lipar
         real(dp), intent(in) :: x,h
         real(dp), intent(inout) :: y(:) ! (n)
         real(dp), intent(inout) :: f(:) ! (n)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr ! nonzero means retry with smaller timestep.
         include 'formats'
         ierr = 0
         f(1) = y(2)
         f(2) = ((1 - y(1)*y(1))*y(2) - y(1))/rpar(1)
         ! the derivatives do not depend on x         
      end subroutine van_der_Pol_derivs


      subroutine solout(nr,xold,x,n,y,rwork,iwork,interp_y,lrpar,rpar,lipar,ipar,irtrn)
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
         ! --- prints solution at equidistant output-points
         ! --- by using "contd8", the continuous collocation solution
         real(dp) :: xout, y1, y2
         integer, parameter :: iprint = 6
         integer :: ierr
         xout = rpar(2)
         irtrn = 1
         if (ipar(1) /= 1) return ! no output
        
         if (nr.eq.1) then
            write (6,99) x,y(1),y(2),nr-1
            xout=x+0.2d0
         else
            do 
               if (x >= xout-1d-10) then
                  ierr = 0
                  y1 = interp_y(1,xout,rwork,iwork,ierr)
                  if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
                  y2 = interp_y(2,xout,rwork,iwork,ierr)
                  if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
                  write (6,99) xout,y1,y2,nr-1
                  xout=xout+0.2d0
                  cycle
               end if
               exit
            end do
        end if
 99     format(1x,'x =',f5.2,'    y =',2e18.10,'    nstep =',i6)
        rpar(2) = xout
        end subroutine solout


      subroutine test_dopri(do_853,show_all)
         logical, intent(in) :: do_853,show_all
         integer, parameter :: nv = 2  ! the number of variables in the van der Pol system of ODEs
         real(dp), parameter :: eps = 1d-3 ! stiffness parameter for van der Pol
         real(dp) :: rtol(1) ! relative error tolerance(s)
         real(dp) :: atol(1) ! absolute error tolerance(s)
         real(dp) :: x ! starting value for the interval of integration
         real(dp) :: xend ! ending value for the interval of integration
         real(dp) :: expect(nv), yprime(nv)
         character (len=64) :: str
         character (len=256) :: dir, fname
         integer, parameter :: lrpar = 2, lipar = 1, nrdens = nv
         integer, parameter :: liwork = nrdens+100, lwork = 11*nv+8*nrdens+100
         real(dp) :: max_abs_yp2, h, max_step_size
         integer :: io_unit, i, lout, iout, idid, itol, j
         integer :: check_liwork, check_lwork, max_steps, ierr
         real(dp), target :: y_ary(nv)
         real(dp), pointer :: y(:)
         real(dp), target :: rpar_ary(lrpar), work_ary(lwork)
         integer, target :: ipar_ary(lipar), iwork_ary(liwork)
         integer, pointer :: ipar(:) ! (lipar)
         real(dp), pointer :: rpar(:) ! (lrpar)
         integer, pointer :: iwork(:)
         real(dp), pointer :: work(:)
         
         ipar => ipar_ary
         rpar => rpar_ary
         work => work_ary
         iwork => iwork_ary
         y => y_ary
         
         write(*,*)        
         write(*,*) 'vdpol'
         if (do_853) then
            write(*,*) 'dop853'
         else
            write(*,*) 'dopri5'
         end if

         x = 0
         xend = 2.0
         y(1) = 2
         y(2) = 0
         
         lout = 0
         max_steps = 0
         max_step_size = 9

         itol = 0 ! scalar tolerances
         iout = 2 ! want dense output
         
         rtol(1) = 1d-4
         atol(1) = 1d-4
         h = 1d-6
         
         rpar(1) = eps
         rpar(2) = 0
         if (show_all) then
            ipar(1) = 1
         else
            ipar(1) = 0
         end if
         
         iwork = 0
         work = 0

         iwork(5)=nrdens ! want dense output for all components
         iwork(4)=1 ! test for stiffness at each step
         
         if (do_853) then
            call dopri5_work_sizes(nv,nrdens,check_liwork,check_lwork)
         else
            call dop853_work_sizes(nv,nrdens,check_liwork,check_lwork)
         end if
         
         if (check_liwork > liwork .or. check_lwork > lwork) then
            write(*,*) 'need to enlarge work arrays for dopri5'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         ierr = 0
         if (do_853) then
            call dop853( &
                  nv,van_der_Pol_derivs,x,y,xend, &
                  h,max_step_size,max_steps, &
                  rtol,atol,itol, &
                  solout,iout,work,lwork,iwork,liwork, &
                  lrpar,rpar,lipar,ipar,lout,idid)
         else
            call dopri5( &
                  nv,van_der_Pol_derivs,x,y,xend, &
                  h,max_step_size,max_steps, &
                  rtol,atol,itol, &
                  solout,iout,work,lwork,iwork,liwork, &
                  lrpar,rpar,lipar,ipar,lout,idid)
         end if

         if (idid /= 1) then ! trouble
            write(*,*) 'idid', idid
            call mesa_error(__FILE__,__LINE__)
         end if
         
         expect(1:2) = (/ 1.7632345401889102d+00, -8.3568868191466206d-01 /)
         
         call show_results(nv,y,expect,show_all)
         if (.not. show_all) return
         
         ! typical: fcn=   21530     step=  1468     accpt=  1345     rejct=  122
         write (6,91) (iwork(j),j=17,20)
 91      format(' fcn=',i8,'     step=',i6,'     accpt=',i6,'     rejct=',i5)
 
         write(*,*)
      
      end subroutine test_dopri



      subroutine test_cash_karp(show_all)
         logical, intent(in) :: show_all
         
         integer, parameter :: nv = 2  ! the number of variables in the van der Pol system of ODEs
         real(dp), parameter :: eps = 1d-3 ! stiffness parameter for van der Pol
         real(dp) :: rtol(1) ! relative error tolerance(s)
         real(dp) :: atol(1) ! absolute error tolerance(s)
         real(dp) :: x ! starting value for the interval of integration
         real(dp) :: xend ! ending value for the interval of integration
         real(dp) :: expect(nv), yprime(nv)
         character (len=64) :: str
         character (len=256) :: dir, fname
         integer, parameter :: lrpar = 2, lipar = 1
         real(dp) :: max_abs_yp2, h, max_step_size
         integer :: io_unit, i, lout, iout, idid, itol, j
         integer :: liwork, lwork, max_steps, ierr
         real(dp), pointer :: work(:)
         integer, pointer :: iwork(:)
         real(dp), target :: y_ary(nv)
         real(dp), pointer :: y(:)
         real(dp), target :: rpar_ary(lrpar)
         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:) ! (lipar)
         real(dp), pointer :: rpar(:) ! (lrpar)
         ipar => ipar_ary
         rpar => rpar_ary

         write(*,*)        
         write(*,*) 'vdpol'
         write(*,*) 'cash_karp'
         
         y => y_ary

         x = 0
         xend = 2.0
         y(1) = 2
         y(2) = 0
         
         lout = 6
         max_steps = 10000
         max_step_size = 9

         itol = 0 ! scalar tolerances
         iout = 0 ! no intermediate output
         
         rtol(1) = 1d-4
         atol(1) = 1d-4
         h = 1d-6
         
         rpar(1) = eps
         rpar(2) = 0
         ipar(1) = 0

         call cash_karp_work_sizes(nv,liwork,lwork)
         allocate(work(lwork), iwork(liwork))
         
         iwork = 0
         work = 0
         
         ierr = 0
         call cash_karp( &
               nv,van_der_Pol_derivs,x,y,xend, &
               h,max_step_size,max_steps, &
               rtol,atol,itol, &
               solout,iout,work,lwork,iwork,liwork, &
               lrpar,rpar,lipar,ipar,lout,idid)

         if (idid /= 1) then ! trouble
            write(*,*) 'idid', idid
            call mesa_error(__FILE__,__LINE__)
         end if

         expect(1:2) = (/ 1.7632345401889102d+00, -8.3568868191466206d-01 /)
         
         call show_results(nv,y,expect,show_all)
         
         if (.not. show_all) return
         
         write (6,91) (iwork(j),j=1,4)
 91      format(' fcn=',i8,'     step=',i6,'     accpt=',i6,'     rejct=',i5)
 
         write(*,*)
         
         deallocate(work, iwork)
      
      end subroutine test_cash_karp

      
      subroutine test_binary_search
         integer, parameter :: n = 100
         integer :: k, loc(3)
         
         real(dp) :: val(3)
         real(dp), target :: vec_ary(n)
         real(dp), pointer :: vec(:)
         include 'formats'
         vec => vec_ary
         
         do k=1,n
            vec(k) = dble(k)*dble(k)
         end do

         write(*,*) 
         write(*,*) 'binary_search, increasing values'
         
         loc = -1
         val = [0d0, dble(n/3)**2 +2, vec(n)+1d0]
         do k=1,3
            loc(k) = binary_search(n, vec, 0, val(k))
            if(loc(k) == 0 .and. val(k) < vec(1))then
               write(*,2) 'val is less than vec(1):', k, val(k), vec(1)
            else if (loc(k) == n .and. val(k) > vec(n))then
               write(*,1) 'val is greater than vec(n):', val(k), vec(n)
            else if (vec(loc(k)) <= val(k) .and. val(k) < vec(loc(k)+1)) then
               write(*,1) 'vec(result)', vec(loc(k))
               write(*,1) 'val', val(k)
               write(*,1) 'vec(result+1)', vec(loc(k)+1)
            else
               write(*,*) 'binary_search failed for increasing-value array'
               call mesa_error(__FILE__,__LINE__)
            end if
            write(*,*) 'okay'
         enddo
         
         ! test decreasing values
         loc = -1
         where(vec /= 0d0) vec = -vec
         where(val /= 0d0) val = -val

         write(*,*)
         write(*,*) 'binary_search, decreasing values'
         do k=1,3
            loc(k) = binary_search(n, vec, 0, val(k))
            if(loc(k) == 0 .and. val(k) > vec(1))then
               write(*,2) 'val is greater than vec(n):', k, val(k), vec(1)
            else if (loc(k) == n .and. val(k) < vec(n))then
               write(*,1) 'val is less than vec(1):', val(k), vec(n)
            else if (vec(loc(k)) >= val(k) .and. val(k) > vec(loc(k)+1)) then
               write(*,1) 'vec(result+1)', vec(loc(k)+1)
               write(*,1) 'val', val(k)
               write(*,1) 'vec(result)', vec(loc(k))
            else
               write(*,*) 'binary_search failed for decreasing-value array'
               call mesa_error(__FILE__,__LINE__)
            end if
            write(*,*) 'okay'
         enddo         
         write(*,*)
      
      end subroutine test_binary_search
      
      
      subroutine test_qsort
         use const_def
         integer, parameter :: n = 100
         integer :: ord(n), i
         real*8 :: a(n)
         include 'formats'
         write(*,*)
         write(*,*) 'qsort into increasing order'
         do i=1,n
            a(i) = sinpi(2.1*dble(i)/dble(n))
         end do
         call qsort(ord, n, a)
         do i=1,n
            write(*,3) 'ord, a(ord)', i, ord(i), a(ord(i))
         end do
         write(*,*)
      end subroutine test_qsort
      
      
      real*8 function g(x) result(y)
         real*8, intent(in) :: x
         y = (x-3)*(x-8)
      end function g
      
      
      subroutine test_find0_quadratic
         real*8 :: xx1, yy1, xx2, yy2, xx3, yy3, x, y
         integer :: ierr
         include 'formats.dek'
         write(*,*) 'test_find0_quadratic'
         xx1 = 1
         yy1 = g(xx1)
         xx2 = 2
         yy2 = g(xx2)
         xx3 = 4
         yy3 = g(xx3)
         x = find0_quadratic(xx1, yy1, xx2, yy2, xx3, yy3, ierr)
         if (ierr /= 0) then
            call mesa_error(__FILE__,__LINE__)
         end if
         y = g(x)
         write(*,1) 'x', x
         write(*,1) 'y', y

         xx1 = 6
         yy1 = g(xx1)
         xx2 = 7
         yy2 = g(xx2)
         xx3 = 8
         yy3 = g(xx3)
         x = find0_quadratic(xx1, yy1, xx2, yy2, xx3, yy3, ierr)
         if (ierr /= 0) then
            call mesa_error(__FILE__,__LINE__)
         end if
         y = g(x)
         write(*,1) 'x', x
         write(*,1) 'y', y
         write(*,*)         
      end subroutine test_find0_quadratic
      
      
      subroutine test_find_max_quadratic
         real*8 :: x1, y1, x2, y2, x3, y3, dx1, dx2, xmax, ymax
         integer :: ierr
         include 'formats.dek'
         write(*,*) 'test_find_max_quadratic'
         dx1 = 1; dx2 = 2; y1 = 1; y2 = 10; y3 = 8
         x1 = 0; x2 = x1 + dx1; x3 = x2 + dx2
         ierr = 0
         call find_max_quadratic(x1, y1, x2, y2, x3, y3, xmax, ymax, ierr)
         if (ierr /= 0) then
            call mesa_error(__FILE__,__LINE__)
         end if
         write(*,1) 'xmax', xmax, 1.85d0
         write(*,1) 'ymax', ymax, 12.4083d0
         write(*,*)         
      end subroutine test_find_max_quadratic

      
            
      end module test_support
