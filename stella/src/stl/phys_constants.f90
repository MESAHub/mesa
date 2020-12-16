module phys_constants

  use kinds,                 only: dp
  use math_constants,        only: p_pi          

  implicit none
  private
  
!  public ::  p_h_bar, p_hplanc, p_c_light, p_boltzk, p_avogar, p_atom_m, p_m_e, p_echarg, &
!              p_cg, p_cms, p_rsol, p_Rgas, &
!              p_carad, p_csigm, p_ergev, p_gradev, p_radc, p_ctomp, p_ccaps, p_ccapz, p_a_fine
              
  public :: write_phys_constants              
       

 !   *** Planck constant [J*s] ***
  real (kind=dp), parameter, public :: p_h_bar = 1.0545715960e-27_dp
  real (kind=dp), parameter, public :: p_h= p_h_bar * (2.0_dp*p_pi) ! p_h = 6.62606875735612676959e-27
!   *** Speed of light in vacuum [cm/s] ***
  real (kind=dp), parameter, public :: p_c= 2.9979245800e+10_dp
!   *** Boltzmann constant [erg/K] ***
  real (kind=dp), parameter, public :: p_boltzk= 1.3806503000e-16_dp
  real (kind=dp), parameter, public :: p_k= p_boltzk ! alias
!   *** Avogadro's Constant  [1/mol] ***
  real (kind=dp), parameter, public :: p_avogar= 6.0221419900e+23_dp ! p_avogar= 6.0221419900e+17_dp p_avogar= 6.0221419900e+23_dp
!   *** Atomic mass unit [g]       ***
  real (kind=dp), parameter, public :: p_atom_m= 1.6605387280e-24_dp
 !   *** Electron mass [g] ***
  real (kind=dp), parameter, public :: p_m_e   = 9.1093818800e-28_dp
 !   *** Elementary charge [CGS] ***
  real (kind=dp), parameter, public :: p_echarg= 4.8032042000e-10_dp
 !   ***  Newtonian constant of gravitation G   [CGS] ***
  real (kind=dp), parameter, public :: p_cg    = 6.6730000000e-08_dp
 !   *** Solar Mass [g]  ***
  real (kind=dp), parameter, public :: p_cms   = 1.9890000000e+33_dp
 !   ***  Solar radius [cm] ***
  real (kind=dp), parameter, public :: p_rsol  = 6.9600000000e+10_dp
 !   ***  Radiation constant  [erg/cm^3 K^4] ***
  real (kind=dp), parameter, public :: p_carad = 7.5657665159e-15_dp
 !   ***  Stefan-Boltzman const [erg/s cm^2 K^4] ***
  real (kind=dp), parameter, public :: p_csigm = 5.6703993512e-05_dp
 !   ***  Electron Volt (energy) [erg] ***
  real (kind=dp), parameter, public :: p_ergev = 1.6021764630e-12_dp
!  !   ***   ***
!   real (kind=dp), parameter, public :: p_gradev= 1.1604505957e+04_dp
 !   ***  Radiation constant  [erg cm^3 K^4]  {stated above}***
  real (kind=dp), parameter, public :: p_radc  = 7.5657665159e-15_dp
 !   ***   ***
  real (kind=dp), parameter, public :: p_ctomp = 4.0062049945e-01_dp
 !   ***   ***
  real (kind=dp), parameter, public :: p_ccaps = 2.6901219544e+01_dp
 !   ***   ***
  real (kind=dp), parameter, public :: p_ccapz = 9.8964056126e+00_dp
 !   *** Molar gas constant [??] { [erg/K mol] => 8.3144721451e+07_dp }  ***
  real (kind=dp), parameter, public :: p_Rgas = p_boltzk * p_avogar ! stella => 8.3144721451e-01_dp
!   *** Fine-structure constant ***
  real (kind=dp), parameter, public :: p_a_fine= 7.2973525330e-03_dp

!!   transformation coefficients
!   *** eV to erg ***
  real (kind=dp), parameter, public :: p_ev2erg = p_ergev  ! p_ev2erg = 1.6021764630e-12_dp
!   *** erg to ev ***
  real (kind=dp), parameter, public :: p_erg2ev =  1._dp /  p_ev2erg !  p_erg2ev = 6.24150974061588121083e11
!   *** freq to eV ***
  real (kind=dp), parameter, public :: p_freq2ev = p_h * p_erg2ev !     p_freq2ev = 4.1356672691e-15
!   *** eV to freq ***
  real (kind=dp), parameter, public :: p_ev2freq = 1._dp / p_freq2ev ! p_ev2freq = 2.4179894922e14
!   *** T [K] to eV ***
  real (kind=dp), parameter, public :: p_T2ev = p_k * p_erg2ev  !  p_T2ev = 8.61734229583e-5
!   *** eV to T [K] ***
  real (kind=dp), parameter, public :: p_ev2T = 1._dp / p_T2ev  ! p_ev2T = 11604.5059563 
  
  contains

  SUBROUTINE write_phys_constants(output_unit)

    INTEGER, INTENT(IN)                      :: output_unit

!   ---------------------------------------------------------------------------

    WRITE (UNIT=output_unit,FMT="(T2,A,/,/,(T2,A))")&
      "*** Fundamental physical constants (SI units) ***",&
      "*** Literature: B. J. Mohr and B. N. Taylor,",&
      "***             CODATA recommended values of the fundamental physical",&
      "***             constants: 1998, Rev. Mod. Phys. 72(2), 351 (2000)."

    WRITE (UNIT=output_unit,FMT="(/,T2,A,T61,ES20.14)")&
      "Speed of light in vacuum [cm/s]",p_c
    WRITE (UNIT=output_unit,FMT="(T2,A,T61,ES20.14)")&
      "Planck constant (h) [J*s]",p_h
    WRITE (UNIT=output_unit,FMT="(T2,A,T61,ES20.14)")&
      "Planck constant (h-bar) [J*s]",p_h_bar
    WRITE (UNIT=output_unit,FMT="(T2,A,T61,ES20.14)")&
      "Elementary charge [cgs]",p_echarg
    WRITE (UNIT=output_unit,FMT="(T2,A,T61,ES20.14)")&
      "Electron mass [g]",p_m_e
    WRITE (UNIT=output_unit,FMT="(T2,A,T61,ES20.14)")&
      "Fine-structure constant",p_a_fine
    WRITE (UNIT=output_unit,FMT="(T2,A,T61,ES20.14)")&
      "Avogrado constant [1/mol]",p_avogar
    WRITE (UNIT=output_unit,FMT="(T2,A,T61,ES20.14)")&
      "Boltzmann constant [erg/K]",p_boltzk
    WRITE (UNIT=output_unit,FMT="(T2,A,T61,ES20.14)")&
      "Atomic mass unit [g]",p_atom_m
    WRITE (UNIT=output_unit,FMT="(T2,A,T61,ES20.14)")&
      "Sun radius [cm]",p_rsol

  END SUBROUTINE write_phys_constants

! *****************************************************************************
  
end module phys_constants

!******************************************************************************
