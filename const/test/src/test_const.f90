program test_const

   use const_def
   use const_lib, only: const_init

   implicit none

   integer :: ierr
   character(len=256) :: my_mesa_dir

   my_mesa_dir = '../..'
   call const_init(my_mesa_dir, ierr)
   if (ierr /= 0) then
      write (*, *) 'const_init failed'
      error stop 1
   end if

   call do_test_const

contains

   subroutine do_test_const

      character(*), parameter :: fmt1 = "(a40, 3x, 1pe24.16)"

      write (*, fmt=fmt1) 'pi', pi
      write (*, fmt=fmt1) 'ln10', ln10
      write (*, fmt=fmt1) 'boltz_sigma', boltz_sigma
      write (*, fmt=fmt1) 'boltz_sigma*4', boltz_sigma*4
      write (*, fmt=fmt1) 'boltz_sigma*4/clight', boltz_sigma*4/clight
      write (*, fmt=fmt1) 'boltz_sigma/clight', boltz_sigma/clight
      write (*, fmt=fmt1) 'crad', crad
      write (*, fmt=fmt1) 'secyer', secyer
      write (*, fmt=fmt1) 'Msun', Msun
      write (*, fmt=fmt1) 'Rsun', Rsun
      write (*, fmt=fmt1) 'Lsun', Lsun
      write (*, fmt=fmt1) 'ly', ly
      write (*, fmt=fmt1) 'm_earth', m_earth
      write (*, fmt=fmt1) 'au', au
      write (*, fmt=fmt1) 'amu', amu
      write (*, fmt=fmt1) 'mn', mn
      write (*, fmt=fmt1) 'mp', mp
      write (*, fmt=fmt1) 'me', me
      write (*, fmt=fmt1) 'planck_h', planck_h
      write (*, fmt=fmt1) 'qe', qe
      write (*, fmt=fmt1) 'avo', avo
      write (*, fmt=fmt1) 'clight', clight
      write (*, fmt=fmt1) 'kerg', kerg

   end subroutine do_test_const

end program test_const
