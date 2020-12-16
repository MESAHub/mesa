      module interp_2D_support
      use const_def
      use interp_2d_lib_db
      use interp_2d_lib_sg
      use math_lib
      use utils_lib, only: mesa_error

      implicit none


      integer, parameter :: num_xpts = 31
      integer, parameter :: num_ypts = 35
      integer, parameter :: sz_per_pt = 4
      
      real, parameter :: pi_sg = 3.14159
      
      real(dp), target :: f_db_ary(sz_per_pt*num_xpts*num_ypts)
      real(dp), pointer :: f_db(:,:,:), f_db1(:)
      real, target :: f_sg_ary(sz_per_pt*num_xpts*num_ypts)
      real, pointer :: f_sg(:,:,:), f_sg1(:)
      integer :: ilinx                           ! =1: x grid is "nearly" equally spaced
      integer :: iliny                           ! =1: y grid is "nearly" equally spaced
      real, target :: xpts_sg_ary(num_xpts), ypts_sg_ary(num_ypts)
      real, pointer :: xpts_sg(:), ypts_sg(:)
      real(dp), target :: xpts_db_ary(num_xpts), ypts_db_ary(num_ypts)
      real(dp), pointer :: xpts_db(:), ypts_db(:)

      contains


      real(dp) function test_fcn_db(x,y)
         real(dp), intent(IN) :: x, y
         real(dp) :: r, r0
         r = DSQRT(x*x+y*y)
         r0 = PI/2
         test_fcn_db = exp(-r/r0) * sin(x) * cos(y)
      end function test_fcn_db


      real function test_fcn_sg(x,y)
         real, intent(IN) :: x, y
         real :: r, r0
         r = SQRT(x*x+y*y)
         r0 = pi_sg/2
         test_fcn_sg = real(test_fcn_db(dble(x),dble(y)))
      end function test_fcn_sg

      
      subroutine get_2D_test_values_db
         
         integer :: i, j
         real(dp) :: x, xmin, xmax, dx
         real(dp) :: y, ymin, ymax, dy
         
         xmin = -PI; xmax = PI; dx = (xmax - xmin) / (num_xpts-1)
         ymin = -PI; ymax = PI; dy = (ymax - ymin) / (num_ypts-1)
         
         do i = 1, num_xpts
            x = xmin + (i-1)*dx
            xpts_db(i) = x
            do j = 1, num_ypts
               y = ymin + (j-1)*dy
               if (i == 1) ypts_db(j) = y
               f_db(1,i,j) = test_fcn_db(x,y)
            end do
         end do

      end subroutine get_2D_test_values_db

      
      subroutine get_2D_test_values_sg
         
         integer :: i, j
         real :: x, xmin, xmax, dx
         real :: y, ymin, ymax, dy
         
         xmin = -pi_sg; xmax = pi_sg; dx = (xmax - xmin) / (num_xpts-1)
         ymin = -pi_sg; ymax = pi_sg; dy = (ymax - ymin) / (num_ypts-1)
         
         do i = 1, num_xpts
            x = xmin + (i-1)*dx
            xpts_sg(i) = x
            do j = 1, num_ypts
               y = ymin + (j-1)*dy
               if (i == 1) ypts_sg(j) = y
               f_sg(1,i,j) = test_fcn_sg(x,y)
            end do
         end do

      end subroutine get_2D_test_values_sg
      
      
      subroutine setup_to_interp_2D_db(bicub_flag)
         logical, intent(in) :: bicub_flag
         integer :: nf2                             ! 2nd dimension of f, nf2.ge.nx
         integer :: ibcxmin                         ! bc flag for x=xmin
         real(dp) :: bcxmin(num_ypts)       ! bc data vs. y at x=xmin
         integer :: ibcxmax                         ! bc flag for x=xmax
         real(dp) :: bcxmax(num_ypts)       ! bc data vs. y at x=xmax
         integer :: ibcymin                         ! bc flag for y=ymin
         real(dp) :: bcymin(num_xpts)       ! bc data vs. x at y=ymin
         integer :: ibcymax                         ! bc flag for y=ymax
         real(dp) :: bcymax(num_xpts)       ! bc data vs. x at y=ymax
         integer :: ier                             ! =0 on exit if there is no error.

         nf2 = num_xpts
         
         if (bicub_flag) then
            !..just use "not a knot" bc's
            ibcxmin = 0; bcxmin = 0
            ibcxmax = 0; bcxmax = 0
            ibcymin = 0; bcymin = 0
            ibcymax = 0; bcymax = 0
            call interp_mkbicub_db(xpts_db,num_xpts,ypts_db,num_ypts,f_db1,nf2,
     >         ibcxmin,bcxmin,ibcxmax,bcxmax,
     >         ibcymin,bcymin,ibcymax,bcymax,
     >         ilinx,iliny,ier)
         else
            call interp_mkbipm_db(xpts_db,num_xpts,ypts_db,num_ypts,f_db1,nf2,ier)
         end if

         if (ier /= 0) then
            write(*,*) 'error'
            call mesa_error(__FILE__,__LINE__)
         end if

      end subroutine setup_to_interp_2D_db
      
      
      subroutine setup_to_interp_2D_sg(bicub_flag)
         logical, intent(in) :: bicub_flag
         integer :: nf2                             ! 2nd dimension of f, nf2.ge.nx
         integer :: ibcxmin                         ! bc flag for x=xmin
         real :: bcxmin(num_ypts)       ! bc data vs. y at x=xmin
         integer :: ibcxmax                         ! bc flag for x=xmax
         real :: bcxmax(num_ypts)       ! bc data vs. y at x=xmax
         integer :: ibcymin                         ! bc flag for y=ymin
         real :: bcymin(num_xpts)       ! bc data vs. x at y=ymin
         integer :: ibcymax                         ! bc flag for y=ymax
         real :: bcymax(num_xpts)       ! bc data vs. x at y=ymax
         integer :: ier                             ! =0 on exit if there is no error.

         nf2 = num_xpts
         
         if (bicub_flag) then
            !..just use "not a knot" bc's
            ibcxmin = 0; bcxmin = 0
            ibcxmax = 0; bcxmax = 0
            ibcymin = 0; bcymin = 0
            ibcymax = 0; bcymax = 0
            call interp_mkbicub_sg(xpts_sg,num_xpts,ypts_sg,num_ypts,f_sg1,nf2,
     >         ibcxmin,bcxmin,ibcxmax,bcxmax,
     >         ibcymin,bcymin,ibcymax,bcymax,
     >         ilinx,iliny,ier)
         else
            call interp_mkbipm_sg(xpts_sg,num_xpts,ypts_sg,num_ypts,f_sg1,nf2,ier)
         end if

         if (ier /= 0) then
            write(*,*) 'error'
            call mesa_error(__FILE__,__LINE__)
         end if

      end subroutine setup_to_interp_2D_sg
      
      
      subroutine eval_2D_interp_db(bicub_flag,x,y,z,dz_dx,dz_dy)
         logical, intent(in) :: bicub_flag
         real(dp), intent(IN) :: x, y
         real(dp), intent(OUT) :: z,dz_dx,dz_dy
         
         integer :: ict(6) ! code specifying output desired
         real(dp) :: fval(6) ! results
         integer :: nf2, ix, jy
         integer :: ier ! error code =0 ==> no error
         
         nf2 = num_xpts
         ier = 0
         fval = 0
         
         if (bicub_flag) then
            ict(1) = 1 ! want f at (x,y)
            ict(2) = 1 ! want df_dx at (x,y)
            ict(3) = 1 ! want df_dy at (x,y)
            ict(4) = 0 ! skip d2f_dx2
            ict(5) = 0 ! skip d2f_dy2
            ict(6) = 0 ! skip d2f_dx_dy
            call interp_evbicub_db(x,y,xpts_db,num_xpts,ypts_db,num_ypts,ilinx,iliny,
     >         f_db1,nf2,ict,fval,ier)
            z = fval(1)
            dz_dx = fval(2)
            dz_dy = fval(3)
         else
            call interp_evbipm_db(
     >         x,y,xpts_db,num_xpts,ypts_db,num_ypts,f_db1,nf2,z,ier)
            dz_dx = 0
            dz_dy = 0
         end if

         if (ier /= 0) then
            write(*,*) 'error'
            call mesa_error(__FILE__,__LINE__)
         end if
      
      end subroutine eval_2D_interp_db
      
      
      subroutine eval_2D_interp_sg(bicub_flag,x,y,z,dz_dx,dz_dy)
         logical, intent(in) :: bicub_flag
         real, intent(IN) :: x, y
         real, intent(OUT) :: z,dz_dx,dz_dy
         
         integer :: ict(6) ! code specifying output desired
         real :: fval(6) ! results
         integer :: nf2, ix, jy
         integer :: ier ! error code =0 ==> no error
         

         nf2 = num_xpts
         
         ier = 0
         fval = 0
         if (bicub_flag) then
            ict(1) = 1 ! want f at (x,y)
            ict(2) = 1 ! want df_dx at (x,y)
            ict(3) = 1 ! want df_dy at (x,y)
            ict(4) = 0 ! skip d2f_dx2
            ict(5) = 0 ! skip d2f_dy2
            ict(6) = 0 ! skip d2f_dx_dy
            call interp_evbicub_sg(x,y,xpts_sg,num_xpts,ypts_sg,num_ypts,ilinx,iliny,
     >         f_sg1,nf2,ict,fval,ier)
            z = fval(1)
            dz_dx = fval(2)
            dz_dy = fval(3)
         else
            call interp_evbipm_sg(
     >            x,y,xpts_sg,num_xpts,ypts_sg,num_ypts,f_sg1,nf2,z,ier)
            dz_dx = 0
            dz_dy = 0
         end if
         
         if (ier /= 0) then
            write(*,*) 'error'
            call mesa_error(__FILE__,__LINE__)
         end if
         
      end subroutine eval_2D_interp_sg
      

      end module interp_2D_support




