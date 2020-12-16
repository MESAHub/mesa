      program test_const
      use const_def
      use const_lib
      implicit none

      integer :: ierr
      character(len=256) :: my_mesa_dir

         my_mesa_dir = '../..'         
         call const_init(my_mesa_dir,ierr)     
         if (ierr /= 0) then
            write(*,*) 'const_init failed'
            error stop 1
         end if        

         call do_test_const

      contains


      subroutine do_test_const
         real(dp) :: val
         integer :: ierr
 
 1       format(a40,3x,1pe24.16)
      
         write(*,1) 'pi', pi
         write(*,1) 'ln10', ln10
         write(*,1) 'boltz_sigma', boltz_sigma
         write(*,1) 'boltz_sigma*4', boltz_sigma*4
         write(*,1) 'boltz_sigma*4/clight', boltz_sigma*4/clight
         write(*,1) 'boltz_sigma/clight', boltz_sigma/clight
         write(*,1) 'crad', crad

         write(*,1) 'secyer', secyer
         write(*,1) 'Msun', Msun
         write(*,1) 'Rsun', Rsun
         write(*,1) 'Lsun', Lsun
         write(*,1) 'ly', ly
         write(*,1) 'm_earth', m_earth
         write(*,1) 'au', au
         write(*,1) 'amu', amu
         write(*,1) 'mn', mn
         write(*,1) 'mp', mp
         write(*,1) 'me', me
         write(*,1) 'planck_h', planck_h
         write(*,1) 'qe', qe
         write(*,1) 'avo', avo
         write(*,1) 'clight', clight
         write(*,1) 'kerg', kerg
                        
      end subroutine do_test_const 


      end program




