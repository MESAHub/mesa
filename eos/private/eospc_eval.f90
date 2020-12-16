! ***********************************************************************
!
!   Copyright (C) 2017-2019  Bill Paxton & The MESA Team
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

      module eospc_eval
      use eos_def
      use auto_diff
      use const_def, only: avo, crad, ln10, arg_not_provided, amu, kerg, dp, qp
      use utils_lib, only: is_bad, mesa_error
      use math_lib

      implicit none
         
      
      contains

      
      subroutine Get_PC_alfa( & 
            rq, logRho, logT, Z, abar, zbar, &
            alfa, d_alfa_dlogT, d_alfa_dlogRho, &
            ierr)
         use const_def
         type (EoS_General_Info), pointer :: rq
         real(dp), intent(in) :: logRho, logT, Z, abar, zbar
         real(dp), intent(out) :: alfa, d_alfa_dlogT, d_alfa_dlogRho
         integer, intent(out) :: ierr
         real(dp) :: logGe0, logGe, logGe_lo, logGe_hi, &
            A, B, dA_dlnT, dA_dlnRho, dB_dlnT, dB_dlnRho, dlogGe_dlogT, dlogGe_dlogRho, &
            logT_lo, logT_hi, logRho_lo, logRho_hi
         
         include 'formats.dek'

         ierr = 0
         
         d_alfa_dlogT = 0d0
         d_alfa_dlogRho = 0d0

         logRho_lo = rq% logRho2_PC_limit ! don't use PC for logRho < this
         logRho_hi = rq% logRho1_PC_limit ! okay for pure PC for logRho > this
         
         if (rq% PC_use_Gamma_limit_instead_of_T) then
            !gamma_e = (qe**2)*((4.0d0/3.0d0)*pi*avo*Rho*zbar/abar)**one_third/(kerg*T)
            !logGe = logGe0 + logRho/3 - logT
            ! where Ge0 = (qe**2)*((4.0d0/3.0d0)*pi*avo*zbar/abar)**one_third/kerg
            logGe0 = log10( & 
                 qe*qe*pow((4.0d0/3.0d0)*pi*avo*zbar/abar,(1d0/3d0))/kerg)
            logGe = logGe0 + logRho/3 - logT
            logGe_lo = rq% log_Gamma_e_all_HELM ! HELM for logGe <= this
            logGe_hi = rq% log_Gamma_e_all_PC ! PC for logGe >= this
            logT_lo = logGe0 + logRho_lo/3 - logGe_lo
            logT_hi = logGe0 + logRho_hi/3 - logGe_hi
            if (logRho <= logRho_lo .or. logGe <= logGe_lo) then
               alfa = 1d0 ! no PC
            else if (logRho >= logRho_hi .and. logGe >= logGe_hi) then
               alfa = 0d0 ! pure PC
            else if (logT >= logT_hi) then ! blend in logGe
               alfa = (logGe_hi - logGe)/(logGe_hi - logGe_lo)
               dlogGe_dlogT = -1d0
               dlogGe_dlogRho = 1d0/3d0
               d_alfa_dlogT = -dlogGe_dlogT/(logGe_hi - logGe_lo)
               d_alfa_dlogRho = -dlogGe_dlogRho/(logGe_hi - logGe_lo)
            else ! blend in logRho
               if (logT >= logT_lo) logRho_lo = (logT_lo + logGe_lo - logGe0)*3
               alfa = (logRho_hi - logRho)/(logRho_hi - logRho_lo)
               d_alfa_dlogRho = -1d0/(logRho_hi - logRho_lo)
            end if
         else
            logT_lo = rq% logT1_PC_limit ! okay for pure PC for logT < this (like logT_all_OPAL)
            logT_hi = rq% logT2_PC_limit ! don't use PC for logT > this (like logT_all_HELM)
            if (logRho <= logRho_lo .or. logT >= logT_hi) then
               alfa = 1d0 ! no PC
            else if (logRho >= logRho_hi .and. logT <= logT_lo) then
               alfa = 0d0 ! pure PC
            else if (logT >= logT_lo) then ! blend in logT
               alfa = (logT - logT_lo)/(logT_hi - logT_lo)
               d_alfa_dlogT = 1/(logT_hi - logT_lo)
            else ! blend in logRho
               alfa = (logRho_hi - logRho)/(logRho_hi - logRho_lo)
               d_alfa_dlogRho = -1d0/(logRho_hi - logRho_lo)
            end if
         end if

      end subroutine Get_PC_alfa


      subroutine get_PC_for_eosdt( &
            handle, dbg, Z, X, abar, zbar, &
            species, chem_id, net_iso, xa, &
            rho, logRho, T, logT, remaining_fraction, &
            res, d_dlnd, d_dlnT, d_dxa, &
            skip, ierr)
         integer, intent(in) :: handle
         logical, intent(in) :: dbg
         real(dp), intent(in) :: &
            Z, X, abar, zbar, remaining_fraction
         integer, intent(in) :: species
         integer, pointer :: chem_id(:), net_iso(:)
         real(dp), intent(in) :: xa(:)
         real(dp), intent(in) :: rho, logRho, T, logT
         real(dp), intent(inout), dimension(nv) :: &
            res, d_dlnd, d_dlnT
         real(dp), intent(inout), dimension(nv, species) :: d_dxa
         logical, intent(out) :: skip
         integer, intent(out) :: ierr
         type (EoS_General_Info), pointer :: rq
         rq => eos_handles(handle)
         call Get_PC_Results(rq, &
            Z, X, abar, zbar, &
            species, chem_id, net_iso, xa, &
            rho, logRho, T, logT, &
            res, d_dlnd, d_dlnT, d_dxa, ierr)
         skip = .false.

         ! zero phase information
         res(i_phase:i_latent_ddlnRho) = 0d0
         d_dlnT(i_phase:i_latent_ddlnRho) = 0d0
         d_dlnd(i_phase:i_latent_ddlnRho) = 0d0

         ! zero all components
         res(i_frac:i_frac+num_eos_frac_results-1) = 0.0
         d_dlnd(i_frac:i_frac+num_eos_frac_results-1) = 0.0
         d_dlnT(i_frac:i_frac+num_eos_frac_results-1) = 0.0

         ! mark this one
         res(i_frac_PC) = 1.0

      end subroutine get_PC_for_eosdt


      subroutine Get_PC_Results( &
            rq, Z, X, abar, zbar, &
            species, chem_id, net_iso, xa, &
            Rho, logRho, T, logT, &
            res, d_dlnRho_c_T, d_dlnT_c_Rho, d_dxa, ierr)
         use pc_eos
         use chem_def, only: chem_isos
         type (EoS_General_Info), pointer :: rq
         real(dp), intent(in) :: Z, X, abar, zbar
         integer, intent(in) :: species
         integer, pointer :: chem_id(:), net_iso(:)
         real(dp), intent(in) :: xa(:)
         real(dp), intent(in) :: Rho, logRho, T, logT
         real(dp), intent(inout) :: res(:) ! (nv)
         real(dp), intent(inout) :: d_dlnRho_c_T(:) ! (nv)
         real(dp), intent(inout) :: d_dlnT_c_Rho(:) ! (nv)
         real(dp), intent(inout) :: d_dxa(:,:) ! (nv, species)
         integer, intent(out) :: ierr
         
         real(dp) :: start_crystal, full_crystal
         real(dp), dimension(species) :: AY, AZion, ACMI
         integer :: i, j
         
         include 'formats.dek'
         
         ierr = 0
         AZion(1:species) = chem_isos% Z(chem_id(1:species))
         ACMI(1:species) = chem_isos% W(chem_id(1:species)) ! this really is atomic weight.
         do j=1,species
            if (xa(j) < rq% mass_fraction_limit_for_PC) then
               AY(j) = 0
            else
               AY(j) = xa(j)/ACMI(j)
            end if
         end do

         start_crystal = rq% PC_Gamma_start_crystal
         full_crystal = rq% PC_Gamma_full_crystal


         call do1(.false.,Rho,T,start_crystal,full_crystal,res,d_dlnT_c_Rho,d_dlnRho_c_T,ierr)
         if (ierr /= 0) return

         ! composition derivatives not provided
         d_dxa = 0
         
         if (is_bad(res(i_lnS))) then
            ierr = -1
            write(*,1) 'res(i_lnS), logRho, logT', res(i_lnS), logRho, logT
            stop 'Get_PC_Results'
         end if

         contains
         
         subroutine do1(show,RHO_real,T_real,start_crystal,full_crystal,res,d_dlnT_c_Rho,d_dlnRho_c_T,ierr)
            logical, intent(in) :: show
            real(dp), intent(in) :: start_crystal, full_crystal
            real(dp), intent(in) :: RHO_real, T_real
            real(dp), intent(out) :: res(:), d_dlnT_c_Rho(:), d_dlnRho_c_T(:) ! (nv)
            integer, intent(out) :: ierr

            integer :: j, LIQSOL
            type(auto_diff_real_2var_order1) :: dsp, dse, dpe, &
               DENS,GAMI,CHI,TPT,GAMImean,TEMP, &
               PnkT,UNkT,SNk,CVNkt,CV,CHIR,CHIT,Tnk,P,Cp,gamma1,gamma3,grad_ad, N, &
               Prad,Pgas,PRADnkT,dE_dRho,dS_dT,dS_dRho,mu,lnfree_e, lnPgas, lnE, lnS, RHO, T
            real(dp) :: Zmean,CMImean,Z2mean
            real(dp), parameter :: UN_T6=0.3157746d0

            include 'formats.dek'
            
            ierr = 0
            
            T = T_real
            T%d1val1 = 1d0

            RHO = RHO_real
            RHO%d1val2 = 1d0
         
            TEMP=T*1d-6/UN_T6 ! T [au]

            if (show) then
               write(*,1) 'RHO', RHO
               write(*,1) 'TEMP', TEMP
               write(*,1) 'start_crystal', start_crystal
               write(*,1) 'full_crystal', full_crystal
               write(*,1) 'AZion(1:species)', AZion(1:species)
               write(*,1) 'ACMI(1:species)', ACMI(1:species)
               write(*,1) 'AY(1:species)', AY(1:species)
               write(*,1) 'xa(1:species)', xa(1:species)
            end if

            call MELANGE9( &
               species,AY,AZion,ACMI,RHO,TEMP, &
               start_crystal,full_crystal,PRADnkT, &
               DENS,Zmean,CMImean,Z2mean,GAMImean,CHI,TPT,LIQSOL, &
               PnkT,UNkT,SNk,CVNkt,CHIR,CHIT,ierr)
            
            if (ierr /= 0) then
               return
               write(*,1) 'RHO', RHO
               write(*,1) 'T', T
               write(*,1) 'logRho', log10(RHO)
               write(*,1) 'logT', log10(T)
               write(*,*) 'ierr from MELANGE9'
               stop 'debug eos'
            end if
            
            if (show) then
               write(*,1) 'PRADnkT', PRADnkT
               write(*,1) 'DENS', DENS
               write(*,1) 'Zmean', Zmean
               write(*,1) 'CMImean', CMImean
               write(*,1) 'Z2mean', Z2mean
               write(*,1) 'GAMI', GAMImean
               write(*,1) 'CHI', CHI
               write(*,1) 'TPT', TPT
               write(*,1) 'PnkT', PnkT
               write(*,1) 'UNkT', UNkT
               write(*,1) 'SNk', SNk
               write(*,1) 'CV', CVNkt
               write(*,1) 'CHIR', CHIR
               write(*,1) 'CHIT', CHIT
               write(*,*)
            end if
            
            Tnk=8.31447d7/CMImean*RHO*T ! n_i kT [erg/cc]
            Pgas = PnkT*Tnk
            if (rq% include_radiation) then
               ! results from MELANGE9 do not include radiation.  add it now.
               CHIT=CHIT*(PnkT/(PnkT+PRADnkT))  + 4.d0*PRADnkT/(PnkT+PRADnkT)
               CHIR=CHIR*(PnkT/(PnkT+PRADnkT))
               PnkT=PnkT+PRADnkT
               UNkT=UNkT+3.d0*PRADnkT
               SNk=SNk+4.d0*PRADnkT
               CVNkt=CVNkt+12.d0*PRADnkT
            end if
            P = PnkT*Tnk
            N = 1/(abar*amu)
            CV = CVNkt*N*kerg
            gamma3 = 1d0 + P/rho * chit/(T*CV)
            gamma1 = chit*(gamma3-1d0) + chir
            grad_ad = (gamma3-1d0)/gamma1
            Cp = CV * gamma1/chir
            dE_dRho = (1d0-chiT)*P/(rho*rho)
            dS_dT = CV/T
            dS_dRho = -P*chiT/(rho*rho * T)
            mu = abar / (1d0 + zbar)
            lnfree_e = log(zbar/abar) ! for complete ionization
            lnPgas = log(Pgas)
            lnE = safe_log(UNkt*N*kerg*T)
            lnS = safe_log(SNk*N*kerg)


            res(i_lnPgas) = lnPgas % val
            res(i_lnE) = lnE % val
            res(i_lnS) = lnS % val
            res(i_grad_ad) = grad_ad % val
            res(i_chiRho) = CHIR % val
            res(i_chiT) = CHIT % val
            res(i_Cp) = Cp % val
            res(i_Cv) = CV % val
            res(i_dE_dRho) = dE_dRho % val
            res(i_dS_dT) = dS_dT % val
            res(i_dS_dRho) = dS_dRho % val
            res(i_mu) = mu % val
            res(i_lnfree_e) = lnfree_e % val
            res(i_gamma1) = gamma1 % val
            res(i_gamma3) = gamma3 % val
            res(i_eta) = CHI % val


            d_dlnT_c_Rho(i_lnPgas) = lnPgas % d1val1 * T % val
            d_dlnT_c_Rho(i_lnE) = lnE % d1val1 * T % val
            d_dlnT_c_Rho(i_lnS) = lnS % d1val1 * T % val
            d_dlnT_c_Rho(i_grad_ad) = grad_ad % d1val1 * T % val
            d_dlnT_c_Rho(i_chiRho) = CHIR % d1val1 * T % val
            d_dlnT_c_Rho(i_chiT) = CHIT % d1val1 * T % val
            d_dlnT_c_Rho(i_Cp) = Cp % d1val1 * T % val
            d_dlnT_c_Rho(i_Cv) = CV % d1val1 * T % val
            d_dlnT_c_Rho(i_dE_dRho) = dE_dRho % d1val1 * T % val
            d_dlnT_c_Rho(i_dS_dT) = dS_dT % d1val1 * T % val
            d_dlnT_c_Rho(i_dS_dRho) = dS_dRho % d1val1 * T % val
            d_dlnT_c_Rho(i_mu) = mu % d1val1 * T % val
            d_dlnT_c_Rho(i_lnfree_e) = lnfree_e % d1val1 * T % val
            d_dlnT_c_Rho(i_gamma1) = gamma1 % d1val1 * T % val
            d_dlnT_c_Rho(i_gamma3) = gamma3 % d1val1 * T % val
            d_dlnT_c_Rho(i_eta) = CHI % d1val1 * T % val


            d_dlnRho_c_T(i_lnPgas) = lnPgas % d1val2 * RHO % val
            d_dlnRho_c_T(i_lnE) = lnE % d1val2 * RHO % val
            d_dlnRho_c_T(i_lnS) = lnS % d1val2 * RHO % val
            d_dlnRho_c_T(i_grad_ad) = grad_ad % d1val2 * RHO % val
            d_dlnRho_c_T(i_chiRho) = CHIR % d1val2 * RHO % val
            d_dlnRho_c_T(i_chiT) = CHIT % d1val2 * RHO % val
            d_dlnRho_c_T(i_Cp) = Cp % d1val2 * RHO % val
            d_dlnRho_c_T(i_Cv) = CV % d1val2 * RHO % val
            d_dlnRho_c_T(i_dE_dRho) = dE_dRho % d1val2 * RHO % val
            d_dlnRho_c_T(i_dS_dT) = dS_dT % d1val2 * RHO % val
            d_dlnRho_c_T(i_dS_dRho) = dS_dRho % d1val2 * RHO % val
            d_dlnRho_c_T(i_mu) = mu % d1val2 * RHO % val
            d_dlnRho_c_T(i_lnfree_e) = lnfree_e % d1val2 * RHO % val
            d_dlnRho_c_T(i_gamma1) = gamma1 % d1val2 * RHO % val
            d_dlnRho_c_T(i_gamma3) = gamma3 % d1val2 * RHO % val
            d_dlnRho_c_T(i_eta) = CHI % d1val2 * RHO % val

         end subroutine do1
     
      end subroutine Get_PC_Results
         

      end module eospc_eval
      
