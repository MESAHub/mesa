!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!  \brief Defines the basic variable types
!  \note
!       Data type definitions; tested on:
!           - IBM AIX xlf90
!           - SGI IRIX  f90
!           - CRAY T3E  f90
!           - DEC ALPHA f90
!           - NAG_F90
!           - SUN
!           - HITACHI
!  \par History
!       Adapted for CP2K by JGH
!  \author Matthias Krack
! *****************************************************************************
module kinds

   implicit none

   private
   public :: sp, dp, print_kind_info, dp_size, sp_size, int_size, int_8, int_4
   public :: default_string_length, default_path_length

#if __SGL
   integer, parameter :: sp = selected_real_kind(6, 30)
   integer, parameter :: dp = selected_real_kind(6, 30)
   ! we rely on this (libraries) but do not check this
   integer, parameter :: dp_size = 4, &
                         int_size = bit_size(0)/8, &
                         sp_size = 4
#else
   integer, parameter :: sp = selected_real_kind(6, 30)
   integer, parameter :: dp = selected_real_kind(14, 200)
   ! we rely on this (libraries) but do not check this
   integer, parameter :: dp_size = 8, &
                         int_size = bit_size(0)/8, &
                         sp_size = 4
#endif

   ! this int holds more than the normal 4byte ints
   ! on standard machines it ought to be an 8byte int but this is not guaranteed
   ! this should also be different from the default integer size (which we require to be 4 bytes)
   integer, parameter :: int_8 = selected_int_kind(10)
   integer, parameter :: int_4 = selected_int_kind(5)
   integer, parameter :: default_string_length = 80
   integer, parameter :: default_path_length = 250
   character(len=1), parameter, public :: default_blank_character(2) = [" ", char(9)]

contains

! *****************************************************************************
!  \brief Print information about the used data types.
!  \par History
!    Adapted by JGH for Cp2k
!   \author Matthias Krack
! *****************************************************************************
   subroutine print_kind_info(iw)

      integer, intent(in) :: iw

      write (iw, '( /, T2, A )') 'DATA TYPE INFORMATION:'

      write (iw, '( /,T2,A,T79,A,2(/,T2,A,T75,I6),3(/,T2,A,T67,E14.8) )') &
         'REAL: Data type name:', 'dp', '      Kind value:', kind(0.0_dp), &
         '      Precision:', precision(0.0_dp), &
         '      Smallest non-negligible quantity relative to 1:', &
         epsilon(0.0_dp), &
         '      Smallest positive number:', tiny(0.0_dp), &
         '      Largest representable number:', huge(0.0_dp)
      write (iw, '( /,T2,A,T79,A,2(/,T2,A,T75,I6),3(/,T2,A,T67,E14.8) )') &
         '      Data type name:', 'sp', '      Kind value:', kind(0.0_sp), &
         '      Precision:', precision(0.0_sp), &
         '      Smallest non-negligible quantity relative to 1:', &
         epsilon(0.0_sp), &
         '      Smallest positive number:', tiny(0.0_sp), &
         '      Largest representable number:', huge(0.0_sp)
      write (iw, '( /,T2,A,T72,A,4(/,T2,A,T61,I20) )') &
         'integer: Data type name:', '(default)', '         Kind value:', &
         kind(0), &
         '         Bit size:', bit_size(0), &
         '         Largest representable number:', huge(0)
      write (iw, '( /,T2,A,T72,A,/,T2,A,T75,I6,/ )') &
         'LOGICAL: Data type name:', '(default)', &
         '         Kind value:', kind(.true.)
      write (iw, '( /,T2,A,T72,A,/,T2,A,T75,I6,/ )') &
         'character: Data type name:', '(default)', &
         '           Kind value:', kind('C')

   end subroutine print_kind_info

end module kinds
