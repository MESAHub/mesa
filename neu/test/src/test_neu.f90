program test_neu
   use neu_support
   use const_lib
   use utils_lib, only : mesa_error
   implicit none
   
   character (len = 32) :: my_mesa_dir
   integer :: ierr
   
   my_mesa_dir = '../..'
   call const_init(my_mesa_dir, ierr)
   if (ierr /= 0) then
      write(*, *) 'const_init failed'
      call mesa_error(__FILE__, __LINE__)
   end if
   call math_init()
   call do_test_neutrinos()

end program
