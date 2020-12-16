!
! args
! francois hebert, jul 18 2010
! 
! for storing and passing variables through the integrators
!

      module def_args

      use def_type, only:dp

      implicit none


      type datastruct

         ! initial conditions
         real (dp) :: ztr    = 0.0_dp
         real (dp) :: neinf  = 0.0_dp
         real (dp) :: chi    = 0.0_dp
         real (dp) :: kt     = 0.0_dp
         real (dp) :: tau    = 0.0_dp
         real (dp) :: cut    = 0.0_dp

         ! background composition information
         integer :: species  = 0
         real (dp), dimension(:), allocatable :: as, zs, xs, niinf


         ! output and diagnostic values
         real (dp) :: znet   = 0.0_dp
         real (dp) :: rws    = 0.0_dp
         real (dp) :: rfinal = 0.0_dp
         real (dp) :: nerws  = 0.0_dp !n_e at r_ws
         real (dp), dimension(:), allocatable :: rex
         !real (dp), dimension(:), allocatable :: nienc !N_i,enc at r_ws


         ! arrays for profile outputs
         integer :: output_length
         real (dp), dimension(:), allocatable :: r, v, z, zb, ne, nb, qitot, qi

         ! intermediate values used in calculation
         real (dp) :: mu     = 0.0_dp
         real (dp) :: dv0    = 0.0_dp
         real (dp) :: xi     = 0.0_dp
         real (dp) :: neloc  = 0.0_dp
         real (dp) :: nbloc  = 0.0_dp
         real (dp) :: qitotloc = 0.0_dp
         real (dp), dimension(:), allocatable :: qiloc

      end type datastruct


      end module def_args
   
