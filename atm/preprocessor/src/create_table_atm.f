! ***********************************************************************
!
!   Copyright (C) 2010-2019  Aaron Dotter & The MESA Team
!
!   This file is part of MESA.
!
!   MESA is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

!     this program creates a base layer of atmosphere boundary conditions
!     this can be combined with tables created from model atmosphere 
!     structures. the main purposes of the files generated from this
!     program is to provide a smooth base that covers all values of Teff
!     and log(g). the model atmosphere grids have limited range of log(g)
!     at fixed Teff and may also have holes.

program create_table_atm

  use atm_def
  use atm_lib
  use chem_def
  use chem_lib, only: chem_init, basic_composition_info
  use const_def
  use const_lib
  use eos_def
  use eos_lib
  use kap_lib
  use utils_lib
  use math_lib

  implicit none

  logical, parameter :: debug = .false.

  integer :: eos_handle, kap_handle, io_out
  integer, parameter :: num_isos = 7, num_logg = 15, num_Teff = 85

  integer, pointer :: chem_id(:), net_iso(:)
  real(dp), pointer :: xa(:)

  integer :: ierr, i_Teff, i_logg
  character(len=256) :: clogZ, output_file, ctau_base
  logical, parameter :: use_cache = .true.
  real(dp) :: M, R, L, X, Y, Z, XC, XN, XO, XNe, XMg, abar, zbar, z2bar, z53bar, kap, err
  real(dp) :: Pextra_factor, Teff, lnP, lnT, tau_base, Teff_out, ye
  real(dp) :: dabar_dx(num_isos), dzbar_dx(num_isos), Xsun, Ysun, Zsun, logZ
  real(dp) :: dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap
  real(dp) :: dlnP_dL, dlnP_dlnR, dlnP_dlnm, dlnP_dlnkap, Xbbn, Ybbn
  real(dp) :: logg_array(num_logg), Teff_array(num_Teff)
  real(dp) :: Pgas(num_logg, num_Teff), T(num_logg, num_Teff)
  integer, parameter :: max_iters = 100
  real(dp), parameter :: errtol = 1d-6
  integer :: iters

  !process command line args
  if( COMMAND_ARGUMENT_COUNT() /= 3 ) then
     stop 'usage: ./create_table_atm [logZ/Zsolar] [output file] [tau_base]'
  endif

  call GET_COMMAND_ARGUMENT(1,clogZ)
  read(clogZ,*) logZ

  call GET_COMMAND_ARGUMENT(2,output_file)

  !mesa initialization
  ierr = 0
  Pextra_factor = 1

  call mesa_init()

  call set_table_composition()

  call GET_COMMAND_ARGUMENT(3,ctau_base)
  read(ctau_base,*) tau_base

  if ( tau_base <= 0 ) stop ' tau > 0 '

  !for table creation:
  M = Msun

  logg_array(:) = (/ -1d0, -0.5d0, 0d0, 0.5d0, 1d0, 1.5d0, 2d0, 2.5d0, 3d0, 3.5d0, 4d0, 4.5d0, 5d0, 5.5d0, 6d0 /)

  Teff_array(1) = 2.0d3
  do i_Teff = 2, num_Teff
     if( Teff_array(i_Teff-1) < 1.3d4) Teff_array(i_Teff) = Teff_array(i_Teff-1) + 2.5d2
     if( Teff_array(i_Teff-1) >=1.3d4) Teff_array(i_Teff) = Teff_array(i_Teff-1) + 1d3
     if( Teff_array(i_Teff-1) >=5d4) Teff_array(i_Teff) = Teff_array(i_Teff-1) + 5d4
  enddo

  if(debug) write(*,*) trim(output_file)

  do i_Teff = 1, num_Teff
     do i_logg = 1, num_logg
        R = sqrt ( standard_cgrav*M / 10d0**logg_array(i_logg) )               
        L = pi*crad*clight * R**2 * Teff_array(i_Teff)**4
        call atm_eval_T_tau_uniform( &
             tau_base, L, R, M, standard_cgrav, 0.2d0*(1 + X), Pextra_factor, &
             ATM_T_TAU_EDDINGTON, eos_proc, kap_proc, errtol, max_iters, .TRUE., &
             Teff, kap, &
             lnT, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
             lnP, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
             ierr)    
        if (ierr /= 0) then
           Pgas(i_logg, i_Teff) = -1
           T(i_logg, i_Teff) = -1
           write(*,*) 'failed in atm_get_gray_and_kap', logg_array(i_logg), Teff_array(i_Teff)
        else           
           T(i_logg, i_Teff) = exp(lnT)
           Pgas(i_logg, i_Teff) = max( 0d0, exp(lnP) - Radiation_Pressure(T(i_logg, i_Teff)) )
        end if
     enddo
  enddo


  open(newunit=io_out, file=output_file)
  !write Pgas
  write(io_out,'("#Teff(K)| Pgas@",15("  log g =",f5.2,1x))') logg_array(1:num_logg)
  do i_Teff = 1, num_Teff
     write(io_out,'(1p,99e15.7)') Teff_array(i_Teff), Pgas(1:num_logg, i_Teff)
  enddo
  !write T
  write(io_out,'("#Teff(K)|    T@",15("  log g =",f5.2,1x))') logg_array(1:num_logg)           
  do i_Teff = 1, num_Teff
     write(io_out,'(1p,99e15.7)') Teff_array(i_Teff), T(1:num_logg, i_Teff)
  enddo
  close (io_out)

  call mesa_shutdown

contains

  subroutine mesa_init

    call const_init(' ',ierr)
    if(ierr/=0) call mesa_error(__FILE__,__LINE__)

    call chem_init('isotopes.data',ierr)
    if(ierr/=0) call mesa_error(__FILE__,__LINE__)

    call eos_init(' ', ' ', ' ', use_cache, ierr)
    if(ierr/=0) call mesa_error(__FILE__,__LINE__)

    eos_handle = alloc_eos_handle(ierr)
    if (ierr /= 0) call mesa_error(__FILE__,__LINE__)

    call kap_init(use_cache, ' ', ierr) 
    if(ierr/=0) call mesa_error(__FILE__,__LINE__)

    kap_handle = alloc_kap_handle(ierr)
    if (ierr /= 0) call mesa_error(__FILE__,__LINE__)

    call atm_init(.false., ierr)
    if (ierr /= 0) call mesa_error(__FILE__,__LINE__)

  end subroutine mesa_init

  !****

  subroutine mesa_shutdown
    call eos_shutdown
    call free_eos_handle(eos_handle)
    call kap_shutdown
    call free_kap_handle(kap_handle)
    call atm_shutdown
  end subroutine mesa_shutdown

  !****

  subroutine set_table_composition()

    integer :: i
    real(dp) :: ye, mass_correction, sumx

    Xbbn   = 0.75d0
    Ybbn   = 0.25d0

    Xsun   = 7.348d-1
    Ysun   = 2.348d-1
    Zsun   = 1.685d-2

    Z = exp10(logZ)*Zsun
    Y = 1.5d0*Z + Ybbn
    X = 1d0 - (Y+Z)

    if(debug)then
       write(*,*) 'X = ', X
       write(*,*) 'Y = ', Y
       write(*,*) 'Z = ', Z
    endif

    XC  = 0.1718d0*Z
    XN  = 0.0503d0*Z
    XO  = 0.4674d0*Z
    XNe = 0.1048d0*Z
    XMg = 0.2057d0*Z ! 1-(XC+XN+XO+XNe)

    allocate(xa(num_isos), chem_id(num_isos), net_iso(num_chem_isos))

    chem_id(:) = (/ ih1, ihe4, ic12, in14, io16, ine20, img24 /)
    net_iso(:) = 0
    do i=1,num_isos
       net_iso(chem_id(i)) = i
    end do
    xa(:) = (/ X, Y, xc, xn, xo, xne, xmg /)
    
    call basic_composition_info( &
       num_isos, chem_id, xa, X, Y, Z, &
       abar, zbar, z2bar, z53bar, ye, mass_correction, sumx)

    if(debug)then
       write(*,*) ' xa = ', xa
       write(*,*) ' sum(xa) = ', sum(xa)
    endif

  end subroutine set_table_composition

  !****

  subroutine eos_proc( &
       lnP, lnT, &
       lnRho, res, dres_dlnRho, dres_dlnT, &
       ierr)
    
    use eos_def, only: num_eos_basic_results
    use eos_lib, only: eosPT_get, radiation_pressure

    real(dp), intent(in)  :: lnP
    real(dp), intent(in)  :: lnT
    real(dp), intent(out) :: lnRho
    real(dp), intent(out) :: res(:)
    real(dp), intent(out) :: dres_dlnRho(:)
    real(dp), intent(out) :: dres_dlnT(:)
    integer, intent(out)  :: ierr

    real(dp) :: T, P, Prad, Pgas, logPgas, rho
    real(dp) :: logRho, dlnRho_dlnPgas, dlnRho_dlnT
    real(dp), dimension(num_eos_basic_results) :: &
         dres_dabar, dres_dzbar

    T = exp(lnT)
    P = exp(lnP)

    Prad = radiation_pressure(T)
    Pgas = max(1E-99_dp, P - Prad)
    logPgas = log10(Pgas)

    call eosPT_get( &
         eos_handle, Z, X, abar, zbar, &
         num_isos, chem_id, net_iso, xa, &
         Pgas, logPgas, T, lnT/ln10, &
         Rho, logRho, dlnRho_dlnPgas, dlnRho_dlnT, &
         res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, ierr)

    lnRho = logRho*ln10

  end subroutine eos_proc

  !****

  subroutine kap_proc( &
       lnRho, lnT, res, dres_dlnRho, dres_dlnT, &
       kap, dlnkap_dlnRho, dlnkap_dlnT, &
       ierr)

    use kap_def, only: num_kap_fracs
    use kap_lib, only: kap_get
    use eos_def, only: i_lnfree_e, i_eta

    real(dp), intent(in)  :: lnRho
    real(dp), intent(in)  :: lnT
    real(dp), intent(in)  :: res(:)
    real(dp), intent(in)  :: dres_dlnRho(:)
    real(dp), intent(in)  :: dres_dlnT(:)
    real(dp), intent(out) :: kap
    real(dp), intent(out) :: dlnkap_dlnRho
    real(dp), intent(out) :: dlnkap_dlnT
    integer, intent(out)  :: ierr

    real(dp) :: kap_fracs(num_kap_fracs), dlnkap_dxa(num_isos)

    call kap_get( &
         kap_handle, num_isos, chem_id, net_iso, xa, &
         lnRho/ln10, lnT/ln10, res(i_lnfree_e), dres_dlnRho(i_lnfree_e), dres_dlnT(i_lnfree_e), &
         res(i_eta), dres_dlnRho(i_eta), dres_dlnT(i_eta), &
         kap_fracs, kap, dlnkap_dlnRho, dlnkap_dlnT, dlnkap_dxa, ierr)
    
  end subroutine kap_proc

end program create_table_atm
