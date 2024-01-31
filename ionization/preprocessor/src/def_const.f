!
! const
! francois hebert, 2009/12/26
! (some of these are from bill paxton's mesa code)
!
! physical & mathematical constants
!

      module def_const

      use def_type, only:dp

      implicit none

      ! math constants

      real (dp), parameter :: sqrt2      = 1.4142135623730950488_dp
      real (dp), parameter :: sqrt3      = 1.73205081_dp

      real (dp), parameter :: pi         = 3.1415926535897932384_dp
      real (dp), parameter :: twopi      = 2.0_dp*pi
      real (dp), parameter :: pi_by2     = pi/2.0_dp
      real (dp), parameter :: fthpi      = 4.0_dp*pi/3.0_dp

      real (dp), parameter :: zero       = 0.0_dp
      real (dp), parameter :: one        = 1.0_dp
      real (dp), parameter :: two        = 2.0_dp
      real (dp), parameter :: one_third  = 1/3.0_dp
      real (dp), parameter :: two_thirds = 2/3.0_dp


      ! physical constants in cgs

      real (dp), parameter :: clight     = 2.99792458e10_dp
      real (dp), parameter :: ggrav      = 6.6742e-8_dp 
      real (dp), parameter :: hplanck    = 6.6260693e-27_dp 
      real (dp), parameter :: hbar       = hplanck / twopi
      real (dp), parameter :: kb         = 1.3806505e-16_dp 
      real (dp), parameter :: mn         = 1.6749286e-24_dp ! neutron mass (g)
      real (dp), parameter :: mp         = 1.6726231e-24_dp ! proton mass (g)
      real (dp), parameter :: me         = 9.1093826e-28_dp ! electron mass (g)
      real (dp), parameter :: qe         = 4.80320440e-10_dp 

      real (dp), parameter :: rbohr      = 5.2917721e-9_dp        ! Bohr radius (cm)
      !real (dp), parameter :: rbohr      = hbar*hbar / (me*qe*qe) ! Bohr radius (cm)
      real (dp), parameter :: fine       = qe*qe / (hbar*clight)  !fine structure constant
      real (dp), parameter :: hion       = 13.605698140_dp        ! hydrogen ionization energy (eV)
      real (dp), parameter :: avogadro   = 6.0221367e23_dp

      real (dp), parameter :: sbsigma    = 5.670400e-5_dp   ! Stefan-Boltzmann constant
      real (dp), parameter :: crad       = sbsigma*4/clight ! radiation constant


      ! astronomical constants
      ! solar age, L, and R values from Bahcall et al, ApJ 618 (2005) 1049-1056.

      real (dp), parameter :: msun       = 1.9892e33_dp  ! solar mass (g)
      real (dp), parameter :: rsun       = 6.9598e10_dp  ! solar radius (cm)
      real (dp), parameter :: lsun       = 3.8418e33_dp  ! solar luminosity (erg s^-1)

      real (dp), parameter :: mearth     = 5.9764e27_dp  ! earth mass (g)
      real (dp), parameter :: rearth     = 6.37e8_dp     ! earth radius (cm)

      real (dp), parameter :: ly         = 9.460528e17_dp      ! light year (cm)
      real (dp), parameter :: pc         = 3.261633_dp * ly    ! parsec (cm)
      real (dp), parameter :: au         = 1.495978921e13_dp   ! astronomical unit (cm)


      ! unit conversions
      ! format: a_b = a per b

      real (dp), parameter :: deg_rad    = pi/180 ! radians per degree
      real (dp), parameter :: rad_deg    = 180/pi ! degrees per radian

      real (dp), parameter :: ev_erg     = 6.241509e11_dp      ! eVs per erg
      reaL (dp), parameter :: erg_ev     = 1.60217733e-12_dp   ! ergs per eV

      real (dp), parameter :: sec_yr     = 3.1558149984e7_dp   ! seconds per year
      real (dp), parameter :: sec_gyr    = 3.1558149984e16_dp  ! seconds per gigayear


      end module def_const

