! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
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

      module neu_def
      use const_def
      implicit none

      real(dp), parameter :: log10Tmin_neu = 7d0, Tmin_neu = 1d7
         ! for T less than this, neu results are all 0

      ! results returned
      
      integer, parameter :: ineu        = 1     ! loss rate (nonnegative) in units of ergs / gram / second
      integer, parameter :: idneu_dT    = 2     ! partial of rate wrt temperature
      integer, parameter :: idneu_dRho  = 3     ! partial of rate wrt density
      integer, parameter :: idneu_dabar = 4     ! partial of rate wrt mean atomic weight
      integer, parameter :: idneu_dzbar = 5     ! partial of rate wrt mean charge
   
      integer, parameter :: num_neu_rvs = 5     ! number of result values per rate
      
      integer, parameter :: pair_neu_type = 1   ! pair production (for reactions like e+ + e- => nu_e + nubar_e)
      integer, parameter :: plas_neu_type = 2   ! plasmon neutrinos (for collective reactions like gamma_plasmon => nu_e + nubar_e)
      integer, parameter :: phot_neu_type = 3   ! photon neutrinos (for reactions like e- + gamma => e- + nu_e + nubar_e)
      integer, parameter :: brem_neu_type = 4   ! bremsstrahlung (for reactions like e- + (z,a) => e- + (z,a) + nu + nubar)
      integer, parameter :: reco_neu_type = 5   ! recombination (for reactions like e- (continuum) => e- (bound) + nu_e + nubar_e)
      
      integer, parameter :: num_neu_types = 5

      end module neu_def

