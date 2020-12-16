      program test_interp
      use const_lib
      use interp_1d_support
      use interp_1d_support_sg
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
         write(*,*) 'test db'
         call do_test
         write(*,*) 'test sg'
         call do_test_sg
      end subroutine test
      
      end program




