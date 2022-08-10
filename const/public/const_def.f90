! ***********************************************************************
!
!   Copyright (C) 2010-2020  The MESA Team
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
!
! ***********************************************************************

      module const_def
      implicit none


      ! real number precision options: single, double, quad
      integer, parameter :: sp = selected_real_kind(p=5)
      integer, parameter :: dp = selected_real_kind(p=15)
      integer, parameter :: qp = selected_real_kind(p=30)

      ! integer precision options
      integer, parameter :: i4 = selected_int_kind(9)
      integer, parameter :: i8 = selected_int_kind(14)


      integer, parameter :: strlen = 256 ! for character (len=strlen)


!
! mathematical and physical constants (in cgs)
!

! math constants
      real(dp), parameter :: pi = 3.1415926535897932384626433832795028841971693993751d0
      real(dp), parameter :: pi2 = pi * pi
      real(dp), parameter :: pi4 = 4*pi
      real(dp), parameter :: eulercon = 0.577215664901532861d0
      real(dp), parameter :: eulernum = 2.71828182845904523536028747135266249d0
      real(dp), parameter :: ln2 = 6.9314718055994529D-01 ! = log(2d0)
      real(dp), parameter :: ln3 = 1.0986122886681096D+00 ! = log(3d0)
      real(dp), parameter :: lnPi = 1.14472988584940017414343_dp ! = log(pi)
      real(dp), parameter :: ln10 = 2.3025850929940455_dp ! = log(10d0)
      real(dp), parameter :: iln10 = 0.43429448190325187_dp ! = 1d0/log(10d0)
      real(dp), parameter :: a2rad = pi/180.0d0 ! angle to radians
      real(dp), parameter :: rad2a = 180.0d0/pi ! radians to angle
      real(dp), parameter :: one_third = 1d0/3d0
      real(dp), parameter :: two_thirds = 2d0/3d0
      real(dp), parameter :: four_thirds = 4d0/3d0
      real(dp), parameter :: five_thirds = 5d0/3d0
      real(dp), parameter :: one_sixth = 1d0/6d0
      real(dp), parameter :: four_thirds_pi = four_thirds*pi
      real(dp), parameter :: ln4pi3 = 1.4324119583011810d0 ! = log(4*pi/3)
      real(dp), parameter :: two_13 = 1.2599210498948730d0 ! = pow(2d0,1d0/3d0)
      real(dp), parameter :: four_13 = 1.5874010519681994d0 ! = pow(4d0,1d0/3d0)
      real(dp), parameter :: sqrt2 = 1.414213562373095d0 ! = sqrt(2)
      real(dp), parameter :: sqrt_2_div_3 = 0.816496580927726d0 ! = sqrt(2/3)

! exact physical constants

      ! CODATA 2018
      real(dp), parameter :: avo = 6.02214076d23 ! Avogadro constant (mole^-1)
      real(dp), parameter :: amu = 1d0 / avo ! atomic mass unit (g)
      real(dp), parameter :: clight = 2.99792458d10 ! speed of light in vacuum (cm s^-1)
      real(dp), parameter :: qe = (clight/10d0) * 1.602176634d-19 ! elementary charge (esu == (g cm^3 s^-2)^(1/2))
      real(dp), parameter :: kerg = 1.380649d-16
      real(dp), parameter :: boltzm = kerg ! Boltzmann constant (erg K^-1)
      real(dp), parameter :: planck_h = 6.62607015d-27 ! Planck constant (erg s)
      real(dp), parameter :: hbar = planck_h / (2*pi)
      real(dp), parameter :: cgas = boltzm*avo ! ideal gas constant (erg K^-1)
      real(dp), parameter :: ev2erg = 1.602176634d-12 ! electron volt (erg)
      real(dp), parameter :: mev_to_ergs = 1d6*ev2erg
      real(dp), parameter :: mev_amu = mev_to_ergs/amu
      real(dp), parameter :: Qconv = mev_to_ergs*avo
      real(dp), parameter :: kev = kerg / ev2erg ! converts temp to ev (ev K^-1)
      real(dp), parameter :: boltz_sigma = (pi*pi * boltzm*boltzm*boltzm*boltzm) / (60 * hbar*hbar*hbar * clight*clight) ! Stefan-Boltzmann constant (erg cm^-2 K^-4 s^-1)
      real(dp), parameter :: crad = boltz_sigma*4/clight ! radiation density constant, AKA "a" (erg cm^-3 K^-4); Prad = crad * T^4 / 3

      ! IAU
      real(dp), parameter :: au = 1.49597870700D13 ! (cm) - exact value defined by IAU 2009, 2012
      real(dp), parameter :: pc = (3.600D3 * rad2a) * au ! (cm) parsec, by definition
      real(dp), parameter :: dayyer = 365.25d0 ! days per (Julian) year
      real(dp), parameter :: secyer = secday*dayyer ! seconds per year
      real(dp), parameter :: secday = 24*60*60  ! seconds in a day
      real(dp), parameter :: secyer = secday*dayyer ! seconds per year
      real(dp), parameter :: ly = clight*secyer ! light year (cm)

! inexact but very well measured physical constants

      real(dp), parameter :: mn = 1.67492749804d-24 ! neutron mass (g)
      real(dp), parameter :: mp = 1.67262192369d-24 ! proton mass (g)
      real(dp), parameter :: me = 9.1093837015d-28  ! electron mass (g)

      real(dp), parameter :: rbohr = 5.29177210903d-9 ! Bohr radius (cm)
      real(dp), parameter :: fine = 7.2973525693d-3   ! fine-structure constant
      real(dp), parameter :: hion = 13.605693122994d0 ! Rydberg constant (eV)

      real(dp), parameter :: sige = 6.6524587321d-25 ! Thomson cross section (cm^2)

      real(dp), parameter :: weinberg_theta  = 0.22290d0 ! sin**2(theta_weinberg)
      real(dp), parameter :: num_neu_fam = 3.0d0 ! number of neutrino flavors = 3.02 plus/minus 0.005 (1998)

! the following quantities are not exact

      real(dp), parameter :: standard_cgrav = 6.67430d-8 ! gravitational constant (g^-1 cm^3 s^-2)

      ! IAU 2015 Resolution B3
      ! standard gravitational parameters = G*M, units cm^3 s^-2
      real(dp), parameter :: mu_sun = 1.3271244d26
      real(dp), parameter :: mu_earth = 3.986004d20
      real(dp), parameter :: mu_jupiter = 1.2668653d23

   ! astronomical constants
      real(dp), parameter :: agesun = 4.57d9  ! solar age (years) from Bahcall et al, ApJ 618 (2005) 1049-1056.
      real(dp), parameter :: Msun = mu_sun / standard_cgrav ! solar mass (g); gravitational mass, not baryonic
      real(dp), parameter :: Rsun = 6.957d10 ! solar radius (cm), IAU 2015 Resolution B3
      real(dp), parameter :: Lsun = 3.828d33 ! solar luminosity (erg s^-1), IAU 2015 Resolution B3
      real(dp), parameter :: Teffsun = 5772.0d0! solar effective temperature (K), IAU 2015 Resolution B3
      real(dp), parameter :: loggsun = 4.4380676273031332_dp ! log10(mu_sun/(Rsun*Rsun)), can't call log10 because we don't have math_lib at this point
      real(dp), parameter :: mbolsun = 4.74d0  ! Bolometric magnitude of the Sun, IAU 2015 Resolution B2

      real(dp), parameter :: m_earth = mu_earth/standard_cgrav! earth mass (g)
      real(dp), parameter :: r_earth = 6.3781d8 ! earth equatorial radius (cm)
      real(dp), parameter :: r_earth_polar = 6.3568d8 ! earth polar radius (cm)

      real(dp), parameter :: m_jupiter = mu_jupiter/standard_cgrav ! jupiter mass (g)
      real(dp), parameter :: r_jupiter = 7.1492d9 ! jupiter equatorial radius (cm)
      real(dp), parameter :: r_jupiter_polar = 6.6854d9 ! jupiter polar radius (cm)
      real(dp), parameter :: semimajor_axis_jupiter = 7.7857d13 ! jupiter semimajor axis (cm)


      ! many routines allow either a value, a log value, or both as args
      ! omitted args are indicated by passing 'arg_not_provided'

      real(dp), parameter :: arg_not_provided = -9d99
      real(dp), parameter :: missing_value = arg_not_provided

      character (len=strlen) :: mesa_dir
      character (len=strlen) :: mesa_data_dir ! = trim(mesa_dir) // '/data'
      character (len=strlen) :: mesa_caches_dir
      character (len=strlen) :: mesa_temp_caches_dir !Temp storage must be local to run

      logical :: use_mesa_temp_cache ! If false then dont use mesa_temp_caches_dir

      ! mixing types
      ! NOTE: some packages may depend on the order
      integer, parameter :: crystallized = -1
      integer, parameter :: no_mixing = 0
      integer, parameter :: convective_mixing = 1
      integer, parameter :: overshoot_mixing = 2
      integer, parameter :: semiconvective_mixing = 3
      integer, parameter :: thermohaline_mixing = 4
      integer, parameter :: rotation_mixing = 5
      integer, parameter :: rayleigh_taylor_mixing = 6
      integer, parameter :: minimum_mixing = 7
      integer, parameter :: anonymous_mixing = 8  ! AKA "WTF_mixing"
      integer, parameter :: leftover_convective_mixing = 9  ! for regions with non-zero conv_vel that are not unstable to convection
                                                            ! used for time dependent convection
      integer, parameter :: phase_separation_mixing = 10
      
      integer, parameter :: number_of_mixing_types =  phase_separation_mixing+1


      contains

      subroutine do_const_init(mesa_dir_init, ierr)
         character (len=*), intent(in) :: mesa_dir_init
         integer, intent(out) :: ierr

         integer :: i, iounit
         character (len=40) :: version_number
         character (len=strlen) :: filename, temp_caches_disable

         ierr = 0

         call get_environment_variable("MESA_CACHES_DIR", mesa_caches_dir)
         !write(*,*) 'MESA_CACHES_DIR "' // trim(mesa_caches_dir) // '"'

         mesa_dir = mesa_dir_init
         if (len_trim(mesa_dir) == 0) then
            call get_environment_variable("MESA_DIR", mesa_dir)
         end if

         if (len_trim(mesa_dir) > 0) then
            mesa_data_dir = trim(mesa_dir) // '/data'
         else
            write(*,*) 'ERROR: you must provide the path to your mesa directory,'
            write(*,*) 'either in your inlist or by setting the MESA_DIR environment variable.'
            ierr = -1
            return
         end if

         !write(*,*) 'mesa_data_dir ' // trim(mesa_data_dir)

         call get_environment_variable("MESA_TEMP_CACHES_DIR", mesa_temp_caches_dir)

         if (len_trim(mesa_temp_caches_dir) == 0) then
            mesa_temp_caches_dir = './.mesa_temp_cache'
         end if

         !Opt out for temp caches
         use_mesa_temp_cache=.true.
         call get_environment_variable("MESA_TEMP_CACHES_DISABLE", temp_caches_disable)

         if (len_trim(temp_caches_disable) > 0) then
            use_mesa_temp_cache = .false.
         end if

      end subroutine do_const_init


      end module const_def
