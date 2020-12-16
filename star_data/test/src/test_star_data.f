      program test_star_data
      use star_data_def
      use star_data_lib
      implicit none

      integer :: ierr
      character(len=256) :: my_mesa_dir

         my_mesa_dir = '../..'         
         call star_data_init(my_mesa_dir,ierr)     
         if (ierr /= 0) then
            write(*,*) 'star_data_init failed'
            error stop 1
         end if        

         call do_test_star_data

      contains


      subroutine do_test_star_data
      
         write(*,*) 'done testing star_data'
                        
      end subroutine do_test_star_data 


      end program




