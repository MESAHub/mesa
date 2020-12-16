      program test_interp
      use const_lib
      use interp_2d_lib_db
      use interp_2d_lib_sg
      use interp_2d_support
      use utils_lib, only: mesa_error

      implicit none
      
      character (len=32) :: my_mesa_dir
      integer :: ierr

      my_mesa_dir = '../..'         
      call const_init(my_mesa_dir,ierr)     
      if (ierr /= 0) then
         write(*,*) 'const_init failed'
         call mesa_error(__FILE__,__LINE__)
      end if
      
      call math_init()
      
      call test
      
      contains
      

      subroutine test
      
         xpts_sg => xpts_sg_ary
         ypts_sg => ypts_sg_ary
         xpts_db => xpts_db_ary
         ypts_db => ypts_db_ary
      
         f_db1 => f_db_ary
         f_db(1:sz_per_pt,1:num_xpts,1:num_ypts) => f_db_ary(1:sz_per_pt*num_xpts*num_ypts)
      
         f_sg1 => f_sg_ary
         f_sg(1:sz_per_pt,1:num_xpts,1:num_ypts) => f_sg_ary(1:sz_per_pt*num_xpts*num_ypts)

         !write(*,*)      
         !call TEST_RENKA790_DB
      
         !write(*,*)      
         !call TEST_RENKA790_SG

         !write(*,*)      
         !call TEST_AKIMA_DB
      
         !write(*,*)      
         !call TEST_AKIMA_SG
      
         write(*,*)      
         call test2D_db(.true.)
      
         write(*,*)      
         call test2D_db(.false.)
      
         write(*,*)      
         call test2D_sg(.true.)
      
         write(*,*)      
         call test2D_sg(.false.)

         write(*,*)
      
      end subroutine test


      subroutine test2D_db(bicub_flag)
         logical, intent(in) :: bicub_flag
         
         integer :: x_points, y_points, i, j
         real(dp) :: x_max, x_min, y_max, y_min, dx, dy, x, y, z, dz_dx, dz_dy, tmp(4)
         
         include 'formats'
         
         write(*,*) 'bicub_flag', bicub_flag
         
         x_points = 2
         y_points = 3
         x_max = 0.8d0*PI; x_min = 0.1d0
         y_max = 0.6d0*PI; y_min = 0.2d0
         dx = (x_max - x_min) / (x_points - 1)
         dy = (y_max - y_min) / (y_points - 1)

         call get_2D_test_values_db
         call setup_to_interp_2D_db(bicub_flag)         

         write(*,*)
         write(*,*) 'interpolant coefficients at midpoint'
         tmp(1:4) = f_db(1:4,num_xpts/2,num_ypts/2)
         write(*,1) 'tmp', tmp
         
         do j = 1, y_points
            y = y_min + (j-1) * dy
            do i = 1, x_points
               x = x_min + (i-1) * dx
               call eval_2D_interp_db(bicub_flag,x,y,z,dz_dx,dz_dy)
               if (bicub_flag) then
                  write(*,3) 'test2D_db', i, j, x,y,z,dz_dx,dz_dy
               else
                  write(*,3) 'test2D_db', i, j, x,y,z
               end if
            end do
         end do
         
         write(*,*)

      end subroutine test2D_db
      

      subroutine test2D_sg(bicub_flag)
         logical, intent(in) :: bicub_flag
         
         integer :: x_points, y_points, i, j
         real :: x_max, x_min, y_max, y_min, dx, dy, x, y, z, dz_dx, dz_dy
         
         include 'formats'
         
         x_points = 2
         y_points = 3
         x_max = 0.8*pi_sg; x_min = 0.1
         y_max = 0.6*pi_sg; y_min = 0.2
         dx = (x_max - x_min) / (x_points - 1)
         dy = (y_max - y_min) / (y_points - 1)

         call get_2D_test_values_sg
         call setup_to_interp_2D_sg(bicub_flag)
         
         do j = 1, y_points
            y = y_min + (j-1) * dy
            do i = 1, x_points
               x = x_min + (i-1) * dx
               call eval_2D_interp_sg(bicub_flag,x,y,z,dz_dx,dz_dy)
               if (bicub_flag) then
                  write(*,3) 'test2D_sg', i, j, x,y,z,dz_dx,dz_dy
               else
                  write(*,3) 'test2D_sg', i, j, x,y,z
               end if
            end do
         end do
         
         write(*,*)

      end subroutine test2D_sg


      end program




