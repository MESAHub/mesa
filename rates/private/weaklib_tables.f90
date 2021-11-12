! ***********************************************************************
!
!   Copyright (C) 2017-2019 Josiah Schwab & The MESA Team
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************


module weaklib_tables

  use const_def
  use math_lib
  use utils_lib, only: mesa_error, is_bad
  use rates_def

  implicit none

  type, extends(weak_rate_table) :: weaklib_rate_table
     integer :: i_ldecay = 1
     integer :: i_lcapture = 2
     integer :: i_lneutrino = 3

   contains

     procedure :: setup => setup_weaklib_table
     procedure :: interpolate => interpolate_weaklib_table

     final :: deallocate_weaklib_rate_table, deallocate_weaklib_rate_table_array
  end type weaklib_rate_table

  interface weaklib_rate_table
     module procedure new_weaklib_rate_table
  end interface weaklib_rate_table

contains


  function new_weaklib_rate_table(T9s, lYeRhos)
    real(dp), intent(in), dimension(:) :: T9s, lYeRhos
    type(weaklib_rate_table) :: new_weaklib_rate_table
    integer :: N

    new_weaklib_rate_table% num_T9 = size(T9s)
    allocate(new_weaklib_rate_table% T9s(new_weaklib_rate_table% num_T9))
    new_weaklib_rate_table% T9s = T9s

    new_weaklib_rate_table% num_lYeRho = size(lYeRhos)
    allocate(new_weaklib_rate_table% lYeRhos(new_weaklib_rate_table% num_lYeRho))
    new_weaklib_rate_table% lYeRhos = lYeRhos

    allocate(new_weaklib_rate_table% data(4, &
         new_weaklib_rate_table% num_T9, &
         new_weaklib_rate_table% num_lYeRho, 3))
  end function new_weaklib_rate_table


  subroutine setup_weaklib_table(table, ierr)
    class(weaklib_rate_table), intent(inout) :: table
    integer, intent(out) :: ierr

    integer :: ii
    do ii = 1, size(table% data, dim=4)
       if (weak_bicubic) then
          call create_bicubic_interpolant(table,ii,ierr)
       else
          call create_pm_T9_interpolant(table,ii,ierr)
       end if
    end do

  contains

    subroutine create_pm_T9_interpolant(table,ii,ierr)
      ! piecewise monotonic interpolation in T9 for each lYeRho in table
      use interp_1d_def
      use interp_1d_lib
      class(weaklib_rate_table), intent(inout) :: table
      integer, intent(in) :: ii
      integer, intent(out) :: ierr
      integer :: j, m, n

      integer :: nx       ! length of x vector (>= 2)
      real(dp), pointer :: x(:) ! (nx)    ! junction points, strictly monotonic
      real(dp), pointer :: f1(:), f(:,:) ! (4,nx)  ! data & interpolation coefficients
      integer, parameter :: nwork = pm_work_size
      real(dp), pointer :: work(:) ! =(nx,nwork)

      ierr = 0

      nx = table % num_T9
      allocate(x(nx), f1(4*nx), work(nx*nwork), stat=ierr)
      if (ierr /= 0) return

      f(1:4,1:nx) => f1(1:4*nx)

      x = table % T9s

      do j=1,table % num_lYeRho
         do m=1,nx
            f(1,m) = table % data(1,m,j,ii)
         end do
         call interp_pm(x, nx, f1, nwork, work, 'create_pm_T9_interpolant', ierr)
         if (ierr /= 0) return
         do n=1,nx
            do m=1,4
               table % data(m,n,j,ii) = real(f(m,n),kind=dp)
            end do
         end do
      end do

      deallocate(x, f1, work)

    end subroutine create_pm_T9_interpolant


    subroutine create_bicubic_interpolant(table,ii,ierr)
      use interp_2d_lib_db
      class(weaklib_rate_table), intent(inout) :: table
      integer, intent(in) :: ii
      integer, intent(out) :: ierr
      integer :: ibcxmin                   ! bc flag for x=xmin
      real(dp), allocatable :: bcxmin(:)    ! bc data vs. y at x=xmin
      integer :: ibcxmax                   ! bc flag for x=xmax
      real(dp), allocatable :: bcxmax(:)     ! bc data vs. y at x=xmax
      integer :: ibcymin                   ! bc flag for y=ymin
      real(dp), allocatable :: bcymin(:)   ! bc data vs. x at y=ymin
      integer :: ibcymax                   ! bc flag for y=ymax
      real(dp), allocatable :: bcymax(:)   ! bc data vs. x at y=ymax
      integer :: il1, il2, j, k, m

      integer :: num_T9, num_lYeRho
      real(dp), pointer :: f1(:), f3(:,:,:)

      num_T9 = table % num_T9
      num_lYeRho = table % num_lYeRho

      allocate(bcxmin(num_lYeRho), bcxmax(num_lYeRho))
      allocate(bcymin(num_T9), bcymax(num_T9))

      ! just use "not a knot" bc's at edges of tables
      ibcxmin = 0; bcxmin(:) = 0
      ibcxmax = 0; bcxmax(:) = 0
      ibcymin = 0; bcymin(:) = 0
      ibcymax = 0; bcymax(:) = 0

      allocate(f1(4*num_T9*num_lYeRho))

      f3(1:4,1:num_T9,1:num_lYeRho) => f1(1:4*num_T9*num_lYeRho)
      do k = 1,num_T9
         do j = 1,4
            do m = 1,num_lYeRho
               f3(j,k,m) = table % data(j,k,m,ii)
            end do
         end do
      end do
      call interp_mkbicub_db( &
           table % T9s, num_T9, &
           table % lYeRhos, num_lYeRho, &
           f1, num_T9, &
           ibcxmin, bcxmin, ibcxmax, bcxmax, &
           ibcymin, bcymin, ibcymax, bcymax, &
           il1, il2, ierr)
      do k = 1,num_T9
         do j = 1,4
            do m = 1,num_lYeRho
               table % data(j,k,m,ii) = f3(j,k,m)
            end do
         end do
      end do
      deallocate(f1, bcxmin, bcxmax, bcymin, bcymax)
    end subroutine create_bicubic_interpolant

  end subroutine setup_weaklib_table

  subroutine interpolate_weaklib_table(table, T9, lYeRho, &
       lambda, dlambda_dlnT, dlambda_dlnRho, &
       Qneu, dQneu_dlnT, dQneu_dlnRho, ierr)
    use const_def, only : dp
    class(weaklib_rate_table), intent(inout) :: table
    real(dp), intent(in) :: T9, lYeRho
    real(dp), intent(out) :: lambda, dlambda_dlnT, dlambda_dlnRho
    real(dp), intent(out) :: Qneu, dQneu_dlnT, dQneu_dlnRho
    integer, intent(out) :: ierr

    integer :: ix, jy          ! target cell in the spline data
    real(dp) :: x0, xget, x1      ! x0 <= xget <= x1;  x0 = xs(ix), x1 = xs(ix+1)
    real(dp) :: y0, yget, y1      ! y0 <= yget <= y1;  y0 = ys(jy), y1 = ys(jy+1)
    real(dp) :: xp, xpi, xp2, xpi2, cx, cxi, hx2, cxd, cxdi, hx, hxi
    real(dp) :: yp, ypi, yp2, ypi2, cy, cyi, hy2, cyd, cydi, hy, hyi

    real(dp) :: delta_T9, dT9, dlYeRho, delta_lYeRho, y_alfa, y_beta, x_alfa, x_beta
    integer :: iT9, ilYeRho

    real(dp) :: ldecay, d_ldecay_dT9, d_ldecay_dlYeRho, &
         lcapture, d_lcapture_dT9, d_lcapture_dlYeRho, &
         lneutrino, d_lneutrino_dT9, d_lneutrino_dlYeRho

    real(dp) :: decay, capture

    logical :: dbg = .false.

    xget = T9
    yget = lYeRho
    
    if (weak_bicubic) then
       call setup_for_bicubic_interpolations
    else
       call setup_for_linear_interp
    endif

    if (weak_bicubic) then

       call do_bicubic_interpolations( &
            table % data(1:4,1:table%num_T9,1:table%num_lYeRho,table%i_ldecay), &
            ldecay, d_ldecay_dT9, d_ldecay_dlYeRho, ierr)

       call do_bicubic_interpolations( &
            table % data(1:4,1:table%num_T9,1:table%num_lYeRho,table%i_lcapture), &
            lcapture, d_lcapture_dT9, d_lcapture_dlYeRho, ierr)

       call do_bicubic_interpolations( &
            table % data(1:4,1:table%num_T9,1:table%num_lYeRho,table%i_lneutrino), &
            lneutrino, d_lneutrino_dT9, d_lneutrino_dlYeRho, ierr)

    else

       call do_linear_interp( &
            table % data(1:4,1:table%num_T9,1:table%num_lYeRho,table%i_ldecay), &
            ldecay, d_ldecay_dT9, d_ldecay_dlYeRho, ierr)

       call do_linear_interp( &
            table % data(1:4,1:table%num_T9,1:table%num_lYeRho,table%i_lcapture), &
            lcapture, d_lcapture_dT9, d_lcapture_dlYeRho, ierr)

       call do_linear_interp( &
            table % data(1:4,1:table%num_T9,1:table%num_lYeRho,table%i_lneutrino), &
            lneutrino, d_lneutrino_dT9, d_lneutrino_dlYeRho, ierr)

    end if

    decay = exp10(ldecay)
    capture = exp10(lcapture)
    lambda = decay + capture

    ! lrates are log10
    ! T9 is linear
    ! lYeRho is log10

    ! drate_dlnT = T9*drate_dT9 = ln10 * rate * d_lrate_dT9
    dlambda_dlnT = ln10*T9*(decay*d_ldecay_dT9 + capture*d_lcapture_dT9)
    ! drate_dlnRho = drate_dlnYeRho (assuming Ye held fixed)
    !              = rate * dlnrate_dlnYeRho = rate * dlrate_dlYeRho
    dlambda_dlnRho = decay*d_ldecay_dlYeRho + capture*d_lcapture_dlYeRho
    Qneu = exp10(lneutrino)/lambda
    dQneu_dlnT = ln10*T9*Qneu*d_lneutrino_dT9 - dlambda_dlnT*Qneu/lambda
    dQneu_dlnRho = Qneu*d_lneutrino_dlYeRho - dlambda_dlnRho*Qneu/lambda
        
  contains

    subroutine find_location ! set ix, jy; x is T9; y is lYeRho
      integer i, j
      real(dp) :: del
      include 'formats'
      ! x0 <= T9 <= x1
      ix = table % num_T9-1 ! since weak_num_T9 is small, just do a linear search
      do i = 2, table % num_T9-1
         if (T9 > table% T9s(i)) cycle
         ix = i-1
         exit
      end do

      ! y0 <= lYeRho <= y1
      jy = table % num_lYeRho-1 ! since weak_num_lYeRho is small, just do a linear search
      do j = 2, table % num_lYeRho-1
         if (lYeRho > table % lYeRhos(j)) cycle
         jy = j-1
         exit
      end do

      x0 = table % T9s(ix)
      x1 = table % T9s(ix+1)
      y0 = table % lYeRhos(jy)
      y1 = table % lYeRhos(jy+1)

    end subroutine find_location

    subroutine setup_for_bicubic_interpolations
      integer i, j
      real(dp) :: del

      include 'formats'

      call find_location

      ! set factors for interpolation

      hx=x1-x0
      hxi=1.0d0/hx
      hx2=hx*hx

      xp=(xget-x0)*hxi
      xpi=1.0d0-xp
      xp2=xp*xp
      xpi2=xpi*xpi

      cx=xp*(xp2-1.0d0)
      cxi=xpi*(xpi2-1.0d0)
      cxd=3.0d0*xp2-1.0d0
      cxdi=-3.0d0*xpi2+1.0d0

      hy=y1-y0
      hyi=1.0d0/hy
      hy2=hy*hy

      yp=(yget-y0)*hyi
      ypi=1.0d0-yp
      yp2=yp*yp
      ypi2=ypi*ypi

      cy=yp*(yp2-1.0d0)
      cyi=ypi*(ypi2-1.0d0)
      cyd=3.0d0*yp2-1.0d0
      cydi=-3.0d0*ypi2+1.0d0

      if (dbg) then
         write(*,2) 'T9', ix, x0, T9, x1
         write(*,2) 'lYeRho', jy, y0, lYeRho, y1
         write(*,1) 'xpi', xpi
         write(*,1) 'ypi', ypi
         write(*,'(A)')
      end if

    end subroutine setup_for_bicubic_interpolations

    subroutine do_bicubic_interpolations(fin, fval, df_dx, df_dy, ierr)
      ! derived from routines in the PSPLINE package written by Doug McCune 
      real(dp), dimension(:,:,:) :: fin ! the spline data array, dimensions (4, nx, ny)
      real(dp), intent(out) :: fval, df_dx, df_dy
      integer, intent(out) :: ierr

      real(dp), parameter :: one_sixth = 1d0/6d0
      real(dp), parameter :: z36th = 1d0/36d0

      ierr = 0

      ! bicubic spline interpolation
      fval = &
           xpi*( &
           ypi*fin(1,ix,jy)  +yp*fin(1,ix,jy+1)) &
           +xp*(ypi*fin(1,ix+1,jy)+yp*fin(1,ix+1,jy+1)) &
           +one_sixth*hx2*( &
           cxi*(ypi*fin(2,ix,jy) +yp*fin(2,ix,jy+1))+ &
           cx*(ypi*fin(2,ix+1,jy)+yp*fin(2,ix+1,jy+1))) &
           +one_sixth*hy2*( &
           xpi*(cyi*fin(3,ix,jy) +cy*fin(3,ix,jy+1))+ &
           xp*(cyi*fin(3,ix+1,jy)+cy*fin(3,ix+1,jy+1))) &
           +z36th*hx2*hy2*( &
           cxi*(cyi*fin(4,ix,jy) +cy*fin(4,ix,jy+1))+ &
           cx*(cyi*fin(4,ix+1,jy)+cy*fin(4,ix+1,jy+1)))

      include 'formats'
      if (dbg) then
         write(*,1) 'fin(1,ix,jy)', fin(1,ix,jy)
         write(*,1) 'fin(1,ix,jy+1)', fin(1,ix,jy+1)
         write(*,1) 'fin(1,ix+1,jy)', fin(1,ix+1,jy)
         write(*,1) 'fin(1,ix+1,jy+1)', fin(1,ix+1,jy+1)
         write(*,1) 'fval', fval

         write(*,'(A)')
         call mesa_error(__FILE__,__LINE__,'debug: do_bicubic_interpolations')
      end if

      ! derivatives of bicubic splines
      df_dx = &
           hxi*( &
           -(ypi*fin(1,ix,jy)  +yp*fin(1,ix,jy+1)) &
           +(ypi*fin(1,ix+1,jy)+yp*fin(1,ix+1,jy+1))) &
           +one_sixth*hx*( &
           cxdi*(ypi*fin(2,ix,jy) +yp*fin(2,ix,jy+1))+ &
           cxd*(ypi*fin(2,ix+1,jy)+yp*fin(2,ix+1,jy+1))) &
           +one_sixth*hxi*hy2*( &
           -(cyi*fin(3,ix,jy)  +cy*fin(3,ix,jy+1)) &
           +(cyi*fin(3,ix+1,jy)+cy*fin(3,ix+1,jy+1))) &
           +z36th*hx*hy2*( &
           cxdi*(cyi*fin(4,ix,jy) +cy*fin(4,ix,jy+1))+ &
           cxd*(cyi*fin(4,ix+1,jy)+cy*fin(4,ix+1,jy+1)))

      df_dy = &
           hyi*( &
           xpi*(-fin(1,ix,jy) +fin(1,ix,jy+1))+ &
           xp*(-fin(1,ix+1,jy)+fin(1,ix+1,jy+1))) &
           +one_sixth*hx2*hyi*( &
           cxi*(-fin(2,ix,jy) +fin(2,ix,jy+1))+ &
           cx*(-fin(2,ix+1,jy)+fin(2,ix+1,jy+1))) &
           +one_sixth*hy*( &
           xpi*(cydi*fin(3,ix,jy) +cyd*fin(3,ix,jy+1))+ &
           xp*(cydi*fin(3,ix+1,jy)+cyd*fin(3,ix+1,jy+1))) &
           +z36th*hx2*hy*( &
           cxi*(cydi*fin(4,ix,jy) +cyd*fin(4,ix,jy+1))+ &
           cx*(cydi*fin(4,ix+1,jy)+cyd*fin(4,ix+1,jy+1)))

    end subroutine do_bicubic_interpolations

    subroutine setup_for_linear_interp
      include 'formats'

      call find_location

      dT9 = T9 - x0
      delta_T9 = x1 - x0
      x_beta = dT9 / delta_T9 ! fraction of x1 result
      x_alfa = 1.0d0 - x_beta ! fraction of x0 result
      if (x_alfa < 0 .or. x_alfa > 1) then
         write(*,1) 'weaklib: x_alfa', x_alfa
         write(*,1) 'T9', T9
         write(*,1) 'x0', x0
         write(*,1) 'x1', x1
         call mesa_error(__FILE__,__LINE__)
      end if

      dlYeRho = lYeRho - y0
      delta_lYeRho = y1 - y0
      y_beta = dlYeRho / delta_lYeRho ! fraction of y1 result
      y_alfa = 1 - y_beta ! fraction of y0 result     
      if (is_bad(y_alfa) .or. y_alfa < 0 .or. y_alfa > 1) then
         write(*,1) 'weaklib: y_alfa', y_alfa
         write(*,1) 'T9', T9
         write(*,1) 'x0', x0
         write(*,1) 'dT9', dT9
         write(*,1) 'delta_T9', delta_T9
         write(*,1) 'lYeRho', lYeRho
         write(*,1) 'y0', y0
         write(*,1) 'dlYeRho', dlYeRho
         write(*,1) 'y1', y1
         write(*,1) 'delta_lYeRho', delta_lYeRho
         write(*,1) 'y_beta', y_beta
         !call mesa_error(__FILE__,__LINE__,'weak setup_for_linear_interp')
      end if

      if (dbg) then
         write(*,2) 'T9', ix, x0, T9, x1
         write(*,2) 'lYeRho', jy, y0, lYeRho, y1
         write(*,1) 'x_alfa, x_beta', x_alfa, x_beta
         write(*,1) 'y_alfa, y_beta', y_alfa, y_beta
         write(*,'(A)')
      end if
      
    end subroutine setup_for_linear_interp

    subroutine do_linear_interp(f, fval, df_dx, df_dy, ierr)
      use interp_1d_lib
      use utils_lib, only: is_bad
      real(dp), dimension(:,:,:) :: f ! (4, nx, ny)         
      real(dp), intent(out) :: fval, df_dx, df_dy
      integer, intent(out) :: ierr

      real(dp) :: fx0, fx1, fy0, fy1

      include 'formats'

      ierr = 0

      fx0 = y_alfa*f(1,ix,jy) + y_beta*f(1,ix,jy+1)
      fx1 = y_alfa*f(1,ix+1,jy) + y_beta*f(1,ix+1,jy+1)

      fy0 = x_alfa*f(1,ix,jy) + x_beta*f(1,ix+1,jy)
      fy1 = x_alfa*f(1,ix,jy+1) + x_beta*f(1,ix+1,jy+1)

      fval = x_alfa*fx0 + x_beta*fx1
      df_dx = (fx1 - fx0)/(x1 - x0)
      df_dy = (fy1 - fy0)/(y1 - y0)    

      if (is_bad(fval)) then
         ierr = -1
         return

         write(*,1) 'x_alfa', x_alfa
         write(*,1) 'x_beta', x_beta
         write(*,1) 'fx0', fx0
         write(*,1) 'fx1', fx1
         write(*,1) 'y_alfa', y_alfa
         write(*,1) 'y_beta', y_beta
         write(*,1) 'f(1,ix,jy)', f(1,ix,jy)
         write(*,1) 'f(1,ix,jy+1)', f(1,ix,jy+1)
         !call mesa_error(__FILE__,__LINE__,'weak do_linear_interp')
      end if

    end subroutine do_linear_interp

  end subroutine interpolate_weaklib_table


  subroutine deallocate_weaklib_rate_table(self)
    type(weaklib_rate_table), intent(inout) :: self
    if(allocated(self % T9s)) deallocate(self % T9s)
    if(allocated(self % lYeRhos)) deallocate(self % lYeRhos)
    if(allocated(self % data)) deallocate(self % data)
  end subroutine deallocate_weaklib_rate_table

  subroutine deallocate_weaklib_rate_table_array(self)
    type(weaklib_rate_table), intent(inout) :: self(:)
    integer :: i
    do i = 1, size(self)
       call deallocate_weaklib_rate_table(self(i))
    end do
  end subroutine deallocate_weaklib_rate_table_array


end module weaklib_tables
