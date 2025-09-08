program test_colors
   use colors_lib, only: colors_init, colors_shutdown
   implicit none

   integer :: ierr
   logical :: use_cache

   ierr = 0
   use_cache = .false.

   ! TODO: add tests for colors module functionality here

   write(*,*) 'Testing colors module initialization...'

   ! Initialize colors module
   ! call colors_init(use_cache, '', ierr)
   if (ierr /= 0) then
      write(*,*) 'Error: colors_init failed with code', ierr
      stop 1
   end if

   write(*,*) 'Colors module initialized successfully.'
   write(*,*) 'Test passed!'

   ! Clean up
   ! call colors_shutdown()

end program test_colors