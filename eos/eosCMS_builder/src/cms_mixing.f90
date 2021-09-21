! ***********************************************************************
!
!   Copyright (C) 2020  Aaron Dotter
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
! ***********************************************************************

! program does the additive volume rule between the pure-H
! and pure-He tables provided by CMS; it produces tables with
! any value 0 <= X <= 1 with Z==0.
!
! for some reason this program does not work with -fopenmp
! even though it does not explicitly use any OpenMP. instead
! it dies with a segfault. i have no idea why!

module cms

   use const_def, only: dp, ln10
   
   implicit none

   logical, parameter :: DBG = .false.
   integer, parameter :: NT = 121
   integer, parameter :: NP = 441

   !unit conversion from CMS units to cgs
   real(dp), parameter :: GPa_to_dyncm2 = 1.0E10_dp
   real(dp), parameter :: MJ_to_erg = 1.0E13_dp
   real(dp), parameter :: kg_to_g   = 1.0E3_dp
   real(dp), parameter :: MJ_kg_to_erg_g = MJ_to_erg/kg_to_g

   real(dp), parameter :: atomic_mass_H   = 1.00784_dp
   real(dp), parameter :: atomic_mass_He  = 4.002603_dp

   type eos_vars
      real(dp) :: logT
      real(dp) :: logP
      real(dp) :: logRho
      real(dp) :: logU
      real(dp) :: logS
      real(dp) :: dlnRho_dlnT_constP
      real(dp) :: dlnRho_dlnP_constT
      real(dp) :: dlnS_dlnT_constP
      real(dp) :: dlnS_dlnP_constT
      real(dp) :: grad_ad
   end type eos_vars

   type table
      character(len=256) :: filename
      type(eos_vars) :: data(NP,NT)
      real(dp) :: H_mass_frac, He_mass_frac
   end type table

contains


   subroutine read_one(t)
      type(table), intent(out) :: t
      integer :: io, i, j
      open(newunit=io, file=trim(t% filename))
      read(io,*) !header
      do i=1,NT
         read(io,*) ! sub-header
         do j = 1, NP
            read(io,'(0p10e15.6)') t% data(j,i)% logT, t% data(j,i)% logP, t% data(j,i)% logRho, &
               t% data(j,i)% logU, t% data(j,i)% logS, t% data(j,i)% dlnRho_dlnT_constP, &
               t% data(j,i)% dlnRho_dlnP_constT, t% data(j,i)% dlnS_dlnT_constP, &
               t% data(j,i)% dlnS_dlnP_constT, t% data(j,i)% grad_ad

            !convert to cgs units: need to do P, U, S
            t% data(j,i)% logP = t% data(j,i)% logP + log10(GPa_to_dyncm2)
            t% data(j,i)% logS = t% data(j,i)% logS + log10(MJ_kg_to_erg_g)
            t% data(j,i)% logU = t% data(j,i)% logU + log10(MJ_kg_to_erg_g)
         enddo
      enddo
      close(io)
   end subroutine read_one

   subroutine write_one(tab)
      type(table), intent(in) :: tab
      integer :: io, i, j
      real(dp) :: chiRho, chiT !U, P, T, S, rho, dE_dP_constT, dE_dT_constP, cp, cv, lnfree_e, mu

      open(newunit=io,file=trim(tab% filename), action='write')
      !write(io,'(a)') '#log T [K]        log P [GPa]   log rho [g/cc]  log U [MJ/kg] log S [MJ/kg/K]' &
      !     // '  dlrho/dlT_P   dlrho/dlP_T    dlS/dlT_P      dlS/dlP_T           grad_ad        dE/dP_T        dE/dT_P'

      write(io,'(a2,f5.3)') 'X=', tab% H_mass_frac
      write(io,'(a)') '#        logT            logP          logRho           logU           logS ' &
         // '   dlrho/dlT_P    dlrho/dlP_T     dlS/dlT_P       dlS/dlP_T        grad_ad         chiRho           chiT'
      do i=1,NT
         do j=1,NP
            !from CMS eqn 5
            chiT = tab% data(j,i)% dlnRho_dlnT_constP / tab% data(j,i)% dlnRho_dlnP_constT
            chiRho = 1._dp / tab% data(j,i)% dlnRho_dlnP_constT

            write(io,'(1p99e15.6)') tab% data(j,i)% logT, tab% data(j,i)% logP, tab% data(j,i)% logRho, &
               tab% data(j,i)% logU, tab% data(j,i)% logS, tab% data(j,i)% dlnRho_dlnT_constP, &
               tab% data(j,i)% dlnRho_dlnP_constT, tab% data(j,i)% dlnS_dlnT_constP, &
               tab% data(j,i)% dlnS_dlnP_constT, tab% data(j,i)% grad_ad, chiRho, chiT
         enddo
      enddo
      close(io)
   end subroutine write_one

   function exp10(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      y=exp(ln10*x)
   end function exp10
   

   subroutine blend_tables(X,Y,Z,ierr)
      type(table), intent(inout) :: X, Y, Z
      integer, intent(out) :: ierr
      integer :: i, j

      ierr = 0

      Z% H_mass_Frac = X% H_mass_frac
      Z% He_mass_frac = Y% He_mass_frac

      do i=1,NT
         do j=1,NP
            call additive_volume(X% data(j,i), Y% data(j,i), Z% data(j,i), X% H_mass_frac, Y% He_mass_frac, ierr)
            if(ierr/=0) then
               write(*,*) 'additive volume failed!'
               return
            endif
         enddo
      enddo
   end subroutine blend_tables


   subroutine additive_volume(eosX,eosY,eosXY,mass_frac_X,mass_frac_Y,ierr)
      type(eos_vars), intent(in) :: eosX, eosY
      real(dp), intent(in) :: mass_frac_X, mass_frac_Y
      type(eos_vars), intent(out) :: eosXY
      integer, intent(out) :: ierr
      real(dp), parameter :: tol = 1.0E-5_dp
      real(dp), parameter :: kerg = 1.380649D-16
      real(dp), parameter :: avo =  6.02214076d23
      real(dp), parameter :: amu = 1.0_dp/avo 

      real(dp) :: Nx, Ny, Ntot, Abar
      real(dp) :: rhoXY, rhoX, rhoY
      real(dp) :: dlnRho_dlnP_XY, dlnRho_dlnT_XY
      real(dp) :: UX, UY, UXY
      real(dp) :: SX, SY, SXY, Smix
      real(dp) :: dlnS_dlnP_XY, dlnS_dlnT_XY
      real(dp) :: grad_ad_num, grad_ad_denom, grad_ad_XY

      ierr = 0

      !safety checks: X+Y == 1
      if ( abs(1.0_dp - mass_frac_X - mass_frac_Y) > tol )then
         write(*,*) 'additive_volume: X + Y != 1'
         ierr = -1
      endif

      !logT1 == logT2
      if (abs(eosX% logT - eosY% logT) > tol ) then
         write(*,*) 'additive_volume: logT_X != logT_Y'
         ierr = -1
      endif

      !logP1 == log
      if (abs(eosX% logP - eosY% logP) > tol ) then
         write(*,*) 'additive_volume: logP_X != logP_Y'
         ierr = -1
      endif

      if(mass_frac_X < 1.0E-12_dp)then
         eosXY = eosY
         return
      elseif(mass_frac_Y < 1.0E-12)then
         eosXY = eosX
         return
      endif

      Nx = mass_frac_X / atomic_mass_H
      Ny = mass_frac_Y / atomic_mass_He
      Ntot = Nx + Ny
      Nx = Nx / Ntot
      Ny = Ny / Ntot

      Abar = Nx*atomic_mass_H + Ny*atomic_mass_He


      ! *** density ***
      rhoX = exp10(eosX% logRho)
      rhoY = exp10(eosY% logRho)
      !CMS19 equation 9
      rhoXY = (rhoX*rhoY)/(rhoY*mass_frac_X + rhoX*mass_frac_Y)

      !derivatives
      dlnRho_dlnP_XY = mass_frac_X*(rhoXY/rhoX)*eosX% dlnRho_dlnP_constT + mass_frac_Y*(rhoXY/rhoY)*eosY% dlnRho_dlnP_constT
      dlnRho_dlnT_XY = mass_frac_X*(rhoXY/rhoX)*eosX% dlnRho_dlnT_constP + mass_frac_Y*(rhoXY/rhoY)*eosY% dlnRho_dlnT_constP


      ! *** energy ***
      UX = exp10( max( min(eosX% logU, 18.0_dp), 9.0_dp))
      UY = exp10( max( min(eosY% logU, 18.0_dp), 9.0_dp))
      UXY = mass_frac_X*UX + mass_frac_Y*UY


      ! *** entropy ***
      SX = exp10( max( min(eosX% logS, 16._dp), 0._dp))
      SY = exp10( max( min(eosY% logS, 16._dp), 0._dp))

      !CMS equation 11
      Smix = - kerg * (Nx*log(Nx) + Ny*log(Ny))/(Abar*amu)

      if(DBG)then
         write(*,*) ' X, Y = ', mass_frac_X, mass_frac_Y
         write(*,*) ' A_X, A_Y = ', atomic_mass_H , atomic_mass_He
         write(*,*) ' rhoX, rhoY = ', rhoX, rhoY
         write(*,*) ' SX, SY = ', SX, SY
         stop
      endif

      !CMS equation 10
      SXY = mass_frac_X*SX + mass_frac_Y*SY + Smix

      !derivative
      dlnS_dlnP_XY = mass_frac_X*(SX/SXY)*eosX% dlnS_dlnP_constT + mass_frac_Y*(SY/SXY)*eosY% dlnS_dlnP_constT
      dlnS_dlnT_XY = mass_frac_X*(SX/SXY)*eosX% dlnS_dlnT_constP + mass_frac_Y*(SY/SXY)*eosY% dlnS_dlnT_constP

      !grad_ad
      grad_ad_num = (mass_frac_X*SX*eosX% dlnS_dlnP_constT + mass_frac_Y*SY*eosY% dlnS_dlnP_constT)
      grad_ad_denom = (mass_frac_X*SX*eosX% dlnS_dlnT_constP + mass_frac_Y*SY*eosY% dlnS_dlnT_constP)
      grad_ad_XY = - grad_ad_num/grad_ad_denom

      grad_ad_XY = min(max(grad_ad_XY,0.1_dp),0.5_dp)

      !final results in eosXY
      eosXY% logT   = eosX% logT
      eosXY% logP   = eosX% logP
      eosXY% logRho = log10(rhoXY)
      eosXY% logU   = log10(UXY)
      eosXY% logS   = log10(SXY)
      eosXY% dlnRho_dlnT_constP = dlnRho_dlnT_XY
      eosXY% dlnRho_dlnP_constT = dlnRho_dlnP_XY
      eosXY% dlnS_dlnT_constP = dlnS_dlnT_XY
      eosXY% dlnS_dlnP_constT = dlnS_dlnP_XY
      eosXY% grad_ad = grad_ad_XY

   end subroutine additive_volume

end module cms



program cms_mixing

   use cms

   implicit none

   integer :: ierr
   real(dp) :: X, Y
   type(table) :: H, J, K
   character(len=256) :: arg

   if(command_argument_count()<2) then
      write(*,*) './cms_mixing [X] [filename]'
      stop
   endif

   call get_command_argument(1,arg)
   read(arg,*) X
   Y = 1.0_dp - X

   call get_command_argument(2, H% filename)

   call get_command_argument(3, J% filename)
   
   call get_command_argument(4, K% filename)

   write(*,*) 'H  input:  ', trim(H% filename)
   write(*,*) 'He input:  ', trim(J% filename)
   write(*,*) 'output to: ', trim(K% filename)
   write(*,*) ' X = ', X
   write(*,*) ' Y = ', Y
   write(*,*) ' X+Y=', X + Y

   H% H_mass_frac  = X
   H% He_mass_frac = 0.0_dp
   call read_one(H)

   J% H_mass_frac  = 0.0_dp
   J% He_mass_frac = Y   
   call read_one(J)

   call blend_tables(H,J,K,ierr)

   if(ierr==0) call write_one(K)

end program cms_mixing
