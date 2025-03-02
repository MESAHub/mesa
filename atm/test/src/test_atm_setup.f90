      module test_atm_setup
      use const_def
      use math_lib
      use atm_lib
      use chem_lib, only: chem_init
      use eos_def
      use eos_lib
      use kap_def
      use kap_lib
      use test_atm_support

      implicit none


      contains


      subroutine setup
         use const_lib
         
         logical, parameter :: use_cache = .true.
         character(len=256) :: my_mesa_dir

         ierr = 0

         my_mesa_dir = '../..'         
         call const_init(my_mesa_dir,ierr)     
         if (ierr /= 0) then
            write(*,*) 'const_init failed'
            call mesa_error(__FILE__,__LINE__)
         end if        
         cgrav = standard_cgrav
         
         call math_init()

         call chem_init('isotopes.data', ierr)
         if (ierr /= 0) then
            write(*,*) 'chem_init failed'
            call mesa_error(__FILE__,__LINE__)
         end if        
         call eos_init(' ', use_cache, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)         
         eos_handle = alloc_eos_handle(ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)        
         call kap_init(use_cache, ' ', ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)         
         kap_handle = alloc_kap_handle(ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)  
         call atm_init(.false., ierr)
         if (ierr /= 0) then
            if (test_verbosely) write(*,*) 'bad return from atm_init'
            call mesa_error(__FILE__,__LINE__)
         end if   
         
         Pextra_factor = 1
               
      end subroutine setup

      end module test_atm_setup




