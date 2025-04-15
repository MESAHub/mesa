program test_colors
   use colors_def
   use colors_lib
   use const_def
   implicit none

   integer :: ierr

   write(*,*) 'Testing colors module initialization...'
   
   ! Initialize colors module
   call colors_init(ierr)
   if (ierr /= 0) then
      write(*,*) 'Error: colors_init failed with code', ierr
      stop 1
   end if

   write(*,*) 'Colors module initialized successfully.'
   write(*,*) 'Test passed!'

   ! Clean up
   call colors_shutdown()

end program test_colors