!
! type
! francois hebert, dec 26 2009
!
! define some numerical types
!
   
      module def_type

      implicit none


      ! a few integer and real types

      integer, parameter :: int4 = selected_int_kind(9)
      integer, parameter :: int2 = selected_int_kind(4)

      integer, parameter :: sp = selected_real_kind(6)
      integer, parameter :: dp = selected_real_kind(12)
      integer, parameter :: qp = selected_real_kind(24)


      end module def_type