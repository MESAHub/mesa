! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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

      module kap_eval
      use utils_lib, only: is_bad, mesa_error
      use kap_def
      use const_def, only: dp, ln10, sige
      use math_lib
      
      implicit none
      
            
      contains
      
      
      subroutine Get_kap_Results( &
           rq, zbar, X, Z, XC, XN, XO, XNe, logRho, logT, &
           lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
           eta, d_eta_dlnRho, d_eta_dlnT, &
           kap_fracs, kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         use const_def
         use utils_lib, only: is_bad
         use auto_diff

         type (Kap_General_Info), pointer :: rq
         real(dp), intent(in) :: X, XC, XN, XO, XNe, Z, zbar
         real(dp), intent(in) :: logRho, logT
         real(dp), intent(in) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
         real(dp), intent(in) :: eta, d_eta_dlnRho, d_eta_dlnT
         real(dp), intent(out) :: kap_fracs(num_kap_fracs)
         real(dp), intent(out) :: kap ! opacity
         real(dp), intent(out) :: dlnkap_dlnRho ! partial derivative at constant T
         real(dp), intent(out) :: dlnkap_dlnT   ! partial derivative at constant Rho
         integer, intent(out) :: ierr ! 0 means AOK.

         real(dp) :: kap_rad, dlnkap_rad_dlnRho, dlnkap_rad_dlnT
         real(dp) :: kap_compton, dlnkap_compton_dlnRho, dlnkap_compton_dlnT

         real(dp) :: logT_Compton_blend_lo, logT_Compton_blend_hi
         real(dp) :: logR_Compton_blend_lo, logR_Compton_blend_hi

         real(dp) :: Rho, T
         type(auto_diff_real_2var_order1) :: logT_auto, logRho_auto, logR_auto
         type(auto_diff_real_2var_order1) :: logkap_rad, logkap_compton
         type(auto_diff_real_2var_order1) :: blend_logT, blend_logR, blend

         real(dp) :: frac_Type2, frac_highT, frac_lowT
         
         logical :: dbg

         include 'formats'

         dbg = rq% dbg
         if (dbg) dbg = & ! check limits
            logT >= rq% logT_lo .and. logT <= rq% logT_hi .and. &
            logRho >= rq% logRho_lo .and. logRho <= rq% logRho_hi .and. &
            X >= rq% X_lo .and. X <= rq% X_hi .and. &
            Z >= rq% Z_lo .and. Z <= rq% Z_hi

         ! auto diff
         ! var1: logRho
         ! var2: logT

         Rho = exp10(logRho)
         logRho_auto% val = logRho
         logRho_auto% d1val1 = 1d0
         logRho_auto% d1val2 = 0d0

         T = exp10(logT)
         logT_auto% val = logT
         logT_auto% d1val1 = 0d0
         logT_auto% d1val2 = 1d0

         logR_auto = logRho_auto - 3d0*logT_auto + 18d0

         if (dbg) write(*,1) 'logRho', logRho
         if (dbg) write(*,1) 'logT', logT
         if (dbg) write(*,1) 'logR', logR_auto% val

         ! initialize fracs
         kap_fracs = 0
         frac_Type2 = 0d0
         frac_lowT = 0d0
         frac_highT = 0d0

         ! blend to Compton at high T
         logT_Compton_blend_hi = rq% logT_Compton_blend_hi
         logT_Compton_blend_lo = logT_Compton_blend_hi - 0.50d0

         ! blend in logT
         if (logT_auto >= logT_Compton_blend_hi) then
            blend_logT = 1d0
         else if (logT_auto > logT_Compton_blend_lo .and. logT_auto <= logT_Compton_blend_hi) then
            blend_logT = (logT_auto - logT_Compton_blend_lo) / (logT_Compton_blend_hi - logT_Compton_blend_lo)
         else
            blend_logT = 0d0
         end if
         ! quintic smoothing
         blend_logT = -blend_logT*blend_logT*blend_logT*(-10d0 + blend_logT*(15d0 - 6d0*blend_logT))


         ! blend to Compton at low R
         logR_Compton_blend_lo = rq% logR_Compton_blend_lo
         logR_Compton_blend_hi = logR_Compton_blend_lo + 0.50d0

         ! blend in logR
         if (logR_auto <= logR_Compton_blend_lo) then
            blend_logR = 1d0
         else if (logR_auto > logR_Compton_blend_lo .and. logR_auto <= logR_Compton_blend_hi) then
            blend_logR = (logR_Compton_blend_hi - logR_auto) / (logR_Compton_blend_hi - logR_Compton_blend_lo)
         else
            blend_logR = 0d0
         endif
         ! quintic smoothing
         blend_logR = -blend_logR*blend_logR*blend_logR*(-10d0 + blend_logR*(15d0 - 6d0*blend_logR))


         ! smoothly combine blends
         blend = 1d0 - (1d0-blend_logT)*(1d0-blend_logR)
         kap_fracs(i_frac_Compton) = blend% val


         if (blend .gt. 0) then ! at least some compton

            if (rq % use_other_compton_opacity) then
               call rq% other_compton_opacity(rq% handle, Rho, T, &
                  lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
                  eta, d_eta_dlnRho, d_eta_dlnT, &
                  kap_compton, dlnkap_compton_dlnRho, dlnkap_compton_dlnT, ierr)
            else
               call Compton_Opacity(rq, Rho, T, &
                  lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
                  eta, d_eta_dlnRho, d_eta_dlnT, &
                  kap_compton, dlnkap_compton_dlnRho, dlnkap_compton_dlnT, ierr)
            end if

            if (ierr /= 0) return

         else ! no Compton

            kap_compton = 1d-30
            dlnkap_compton_dlnRho = 0d0
            dlnkap_compton_dlnT = 0d0

         end if

         ! pack into auto_diff type
         logkap_compton = log10(kap_compton)
         logkap_compton% d1val1 = dlnkap_compton_dlnRho
         logkap_compton% d1val2 = dlnkap_compton_dlnT


         if (blend .lt. 1) then ! at least some tables

            if (rq% use_other_radiative_opacity) then
               call rq% other_radiative_opacity( &
                  rq% handle, X, Z, XC, XN, XO, XNe, logRho, logT, &
                  frac_lowT, frac_highT, frac_Type2, &
                  kap_rad, dlnkap_rad_dlnRho, dlnkap_rad_dlnT, ierr)
            else
               call Get_kap_Results_blend_T( &
                  rq, X, Z, XC, XN, XO, XNe, logRho, logT, &
                  frac_lowT, frac_highT, frac_Type2, &
                  kap_rad, dlnkap_rad_dlnRho, dlnkap_rad_dlnT, ierr)
            end if
            if (ierr /= 0) return

            ! revise reported fractions based on Compton blend
            frac_lowT = (1d0 - blend% val) * frac_lowT
            frac_highT = (1d0 - blend% val) * frac_highT
            ! the value of frac_Type2 doesn't need revised if it represents
            ! the fraction of highT opacities provided by the Type2 tables
            ! frac_Type2 = (1d0 - blend% val) * frac_Type2

         else ! no tables

            kap_rad = 1d-30
            dlnkap_rad_dlnRho = 0d0
            dlnkap_rad_dlnT = 0d0

         end if

         ! pack into auto_diff type
         logkap_rad = log10(kap_rad)
         logkap_rad% d1val1 = dlnkap_rad_dlnRho
         logkap_rad% d1val2 = dlnkap_rad_dlnT

         ! pack kap_fracs from tables
         kap_fracs(i_frac_lowT) = frac_lowT
         kap_fracs(i_frac_highT) = frac_highT
         kap_fracs(i_frac_Type2) = frac_Type2

         ! do blend
         logkap_rad = blend * logkap_compton + (1d0 - blend) * logkap_rad

         ! unpack auto_diff
         kap_rad = exp10(logkap_rad% val)
         dlnkap_rad_dlnRho = logkap_rad% d1val1
         dlnkap_rad_dlnT = logkap_rad% d1val2

         call combine_rad_with_conduction( &
            rq, Rho, logRho, T, logT, zbar, &
            kap_rad, dlnkap_rad_dlnRho, dlnkap_rad_dlnT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         
      end subroutine Get_kap_Results
      
      
      subroutine Get_kap_Results_blend_T( &
           rq, X, Z, XC, XN, XO, XNe, logRho_in, logT_in, &
           frac_lowT, frac_highT, frac_Type2, kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)

        use kap_def
        use utils_lib, only: is_bad

        ! INPUT
        type (Kap_General_Info), pointer :: rq
        real(dp), intent(in) :: X, Z, XC, XN, XO, XNe ! composition    
        real(dp), intent(in) :: logRho_in ! density
        real(dp), intent(in) :: logT_in ! temperature

        ! OUTPUT
        real(dp), intent(out) :: frac_lowT, frac_highT, frac_Type2
        real(dp), intent(out) :: kap ! opacity
        real(dp), intent(out) :: dlnkap_dlnRho ! partial derivative at constant T
        real(dp), intent(out) :: dlnkap_dlnT   ! partial derivative at constant Rho
        integer, intent(out) :: ierr ! 0 means AOK.

        real(dp) :: logRho, Rho, logT, T
        real(dp) :: frac_Type1, dZ, frac_X, frac_dZ, &
             dXC, dXO, fC, fN, fO, fNe, fHeavy, ZHeavy, dXsum

        real(dp) :: lowT_logT_max, logT_min

        logical :: clipped_Rho
        
        logical :: dbg

        real(dp) :: alfa, beta, &
             logKap, logKap1, logKap2, &
             kap1, dlnkap1_dlnRho, dlnkap1_dlnT, &
             kap2, dlnkap2_dlnRho, dlnkap2_dlnT, &
             lower_bdy, upper_bdy, alfa0, d_alfa0_dlnT, &
             d_alfa_dlnT, d_beta_dlnT, dkap_dlnT


        include 'formats'

        dbg = .false.

        logRho = logRho_in; Rho = exp10(logRho)
        logT = logT_in; T = exp10(logT)

        if (dbg) write(*,1) 'Get_kap_blend_logT logT', logT

        clipped_Rho = .false.
        if (logRho < kap_min_logRho) then
           logRho = kap_min_logRho
           rho = exp10(kap_min_logRho)
           clipped_Rho = .true.
        end if

        frac_lowT = 0d0
        frac_highT = 0d0
        frac_Type2 = 0d0

        if (rq% use_Type2_opacities .and. &
            associated(kap_co_z_tables(rq% kap_CO_option)% ar)) then
           logT_min = kap_co_z_tables(rq% kap_CO_option)% ar(1)% x_tables(1)% logT_min
        else
           logT_min = kap_z_tables(rq% kap_option)% ar(1)% x_tables(1)% logT_min
        end if

        lower_bdy = max(rq% kap_blend_logT_lower_bdy, logT_min)
        if (dbg) write(*,1) 'Get_kap_blend_logT lower_bdy', lower_bdy, &
             rq% kap_blend_logT_lower_bdy, logT_min
        if (logT <= lower_bdy) then ! all lowT
           if (dbg) write(*,1) 'all lowT'
           call Get_kap_lowT_Results(rq, &
                X, Z, XC, XN, XO, XNe, logRho, logT, &
                frac_Type2, kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
           if (clipped_Rho) dlnkap_dlnRho = 0d0
           frac_lowT = 1d0
           return
        end if


        select case (rq% kap_lowT_option)
        case (kap_lowT_kapCN)
           lowT_logT_max = kapCN_max_logT
        case (kap_lowT_AESOPUS)
           lowT_logT_max = kA % max_logT
        case default
           lowT_logT_max = kap_lowT_z_tables(rq% kap_lowT_option)% ar(1)% x_tables(1)% logT_max           
        end select


        upper_bdy = min(rq% kap_blend_logT_upper_bdy, lowT_logT_max)
        if (dbg) write(*,1) 'Get_kap_blend_logT upper_bdy', upper_bdy, &
             rq% kap_blend_logT_upper_bdy, lowT_logT_max
        if (logT >= upper_bdy) then ! no lowT
           if (dbg) write(*,1) 'no lowT'
           call Get_kap_highT_Results(rq, &
                X, Z, XC, XN, XO, XNe, logRho, logT, &
                frac_Type2, kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
           if (clipped_Rho) dlnkap_dlnRho = 0d0
           frac_highT = 1d0
           return
        end if

        ! in blend region

        call Get_kap_lowT_Results(rq, &
             X, Z, XC, XN, XO, XNe, logRho, logT, &
             frac_Type2, kap2, dlnkap2_dlnRho, dlnkap2_dlnT, ierr)
        if (clipped_Rho) dlnkap2_dlnRho = 0d0
        if (ierr /= 0) return

        call Get_kap_highT_Results(rq, &
             X, Z, XC, XN, XO, XNe, logRho, logT, &
             frac_Type2, kap1, dlnkap1_dlnRho, dlnkap1_dlnT, ierr)
        if (clipped_Rho) dlnkap1_dlnRho = 0d0
        if (ierr /= 0) return

        alfa0 = (logT - lower_bdy) / (upper_bdy - lower_bdy)
        d_alfa0_dlnT = 1d0/(upper_bdy - lower_bdy)/ln10

        ! must smooth the transitions near alfa = 0.0 and 1.0  
        ! Rich Townsend's smoothing function for this
        alfa = -alfa0*alfa0*alfa0*(-10d0 + alfa0*(15d0 - 6d0*alfa0))
        d_alfa_dlnT = 30d0*(alfa0 - 1d0)*(alfa0 - 1d0)*alfa0*alfa0*d_alfa0_dlnT

        beta = 1d0 - alfa
        d_beta_dlnT = -d_alfa_dlnT


        logKap1 = log10(kap1); logKap2 = log10(kap2)

        logKap = alfa*logKap1 + beta*logKap2
        dlnkap_dlnRho = alfa*dlnkap1_dlnRho + beta*dlnkap2_dlnRho
        dlnkap_dlnT = &
             alfa*dlnkap1_dlnT + beta*dlnkap2_dlnT + &
             d_alfa_dlnT*logKap1*ln10 + d_beta_dlnT*logKap2*ln10
        kap = exp10(logKap)

        frac_lowT = beta
        frac_highT = alfa

        if (dbg) then
           write(*,1) 'alfa', alfa
           write(*,1) 'kap1 (lowT)', kap1
           write(*,1) 'beta', beta
           write(*,1) 'kap2 (highT)', kap2
           write(*,1) 'kap', kap
           write(*,*)
        end if

        if (dbg) write(*,1) 'Get_kap_blend_logT kap', kap

      end subroutine Get_kap_Results_blend_T


      subroutine Get_kap_lowT_Results( &
           rq, &
           X, Z, XC, XN, XO, XNe, logRho_in, logT_in, &
           frac_Type2, kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)

        use kap_def

        use kap_eval_fixed
        use kap_eval_co
        use kapcn
        use kap_aesopus

        use utils_lib, only: is_bad

        ! INPUT
        type (Kap_General_Info), pointer :: rq
        real(dp), intent(in) :: X, Z, XC, XN, XO, XNe ! composition    
        real(dp), intent(in) :: logRho_in ! density
        real(dp), intent(in) :: logT_in ! temperature
        ! free_e := total combined number per nucleon of free electrons and positrons

        ! OUTPUT
        real(dp), intent(out) :: frac_Type2
        real(dp), intent(out) :: kap ! opacity
        real(dp), intent(out) :: dlnkap_dlnRho ! partial derivative at constant T
        real(dp), intent(out) :: dlnkap_dlnT   ! partial derivative at constant Rho
        integer, intent(out) :: ierr ! 0 means AOK.

        real(dp) :: logRho, Rho, logT, T
        real(dp) :: frac_Type1, Zbase, dZ, frac_X, frac_dZ, &
             dXC, dXO, fCO, fC, fN, fO, fNe, fHeavy, ZHeavy, dXsum

        real(dp) :: logKap

        logical :: dbg = .false.

        ierr = -1 ! should be set by each case, otherwise something is wrong
        frac_Type2 = 0d0

        logRho = logRho_in; Rho = exp10(logRho)
        logT = logT_in; T = exp10(logT)

        Zbase = rq% Zbase

        select case (rq% kap_lowT_option)

        case (kap_lowT_kapCN)
           if (dbg) write(*,*) 'Calling kapCN for lowT'

           if (Zbase < 0d0) then
              write(*,*) 'must supply Zbase for kapCN opacities', Zbase
              call mesa_error(__FILE__,__LINE__)
              ierr = -1
              return
           end if

           fC = XC / (kapCN_ZC * Zbase)
           fN = XN / (kapCN_ZN * Zbase)
           call kapCN_get(Zbase, X, fC, fN, logRho, logT, &
                kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)

        case (kap_lowT_AESOPUS)
           if (dbg) write(*,*) 'Calling AESOPUS for lowT'

           if (Zbase < 0d0) then
              write(*,*) 'must supply Zbase for AESOPUS opacities', Zbase
              ierr = -1
              return
           end if

           call AESOPUS_get(Zbase, X, XC, XN, XO, logRho, logT, &
                kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)

           
        case default
           if (rq% kap_lowT_option == kap_lowT_test) then
              if (dbg) write(*,*) 'Calling test for lowT'
           else if (rq% kap_lowT_option == kap_lowT_Freedman11) then
              if (dbg) write(*,*) 'Calling Freedman for lowT'
           else
              if (dbg) write(*,*) 'Calling Ferg for lowT'
           end if

           call Get1_kap_fixed_metal_Results( &
                kap_lowT_z_tables(rq% kap_lowT_option)% ar, num_kap_lowT_Zs(rq% kap_lowT_option), rq, Z, X, Rho, logRho, T, logT, &
                logKap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
           kap = exp10(logKap)
        end select

        return

      end subroutine Get_kap_lowT_Results


      subroutine Get_kap_highT_Results(rq, &
           X, Z, XC_in, XN_in, XO_in, XNe_in, logRho_in, logT_in, &
           frac_Type2, kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)

        use kap_def

        use kap_eval_fixed
        use kap_eval_co

        use utils_lib, only: is_bad

        ! INPUT
        type (Kap_General_Info), pointer :: rq
        real(dp), intent(in) :: X, Z, XC_in, XN_in, XO_in, XNe_in ! composition    
        real(dp), intent(in) :: logRho_in ! density
        real(dp), intent(in) :: logT_in ! temperature

        ! OUTPUT
        real(dp), intent(out) :: frac_Type2
        real(dp), intent(out) :: kap ! opacity
        real(dp), intent(out) :: dlnkap_dlnRho ! partial derivative at constant T
        real(dp), intent(out) :: dlnkap_dlnT   ! partial derivative at constant Rho
        integer, intent(out) :: ierr ! 0 means AOK.

        real(dp) :: logRho, Rho, logT, T
        real(dp) :: XC, XN, XO, XNe
        real(dp) :: frac_Type1, Zbase, dZ, frac_X, frac_dZ, &
             dXC, dXO, fC, fN, fO, fNe, fHeavy, ZHeavy, dXsum

        real(dp) :: max_frac_Type2

        real(dp) :: beta, kap_beta, logkap_beta, dlnkap_beta_dlnRho, dlnkap_beta_dlnT, &
             alfa, kap_alfa, logkap_alfa, dlnkap_alfa_dlnRho, dlnkap_alfa_dlnT

        logical, parameter :: dbg = .false.

        include 'formats'

        logRho = logRho_in; Rho = exp10(logRho)
        logT = logT_in; T = exp10(logT)

        Zbase = rq% Zbase

        if (rq% use_Type2_opacities) then
           if (Zbase < 0d0) then
              write(*,1) 'must supply Zbase for Type2 opacities', Zbase
              call mesa_error(__FILE__,__LINE__)
              ierr = -1
              return
           end if
        end if

        if (rq% use_Type2_opacities) then
           max_frac_Type2 = 1d0
           XC = XC_in
           XN = XN_in
           XO = XO_in
           XNe = XNe_in
        else
           max_frac_Type2 = 0d0
           XC = 0d0
           XN = 0d0
           XO = 0d0
           XNe = 0d0
        endif


        dXC = 0d0
        dxO = 0d0
        frac_Type2 = max_frac_Type2


        ! Type2 opacities use a base metallicity and enhancements to C and O
        ! from the 3 definitions
        ! Z = Zbase + dXC + dXO ! total metals
        ! XC = dXC + base_fC*Zbase ! total mass fraction of carbon
        ! XO = dXO + base_fO*Zbase ! total mass fraction of oxygen
        ! we get expressions for the 3 parameters, Zbase, dXC, and dXO
        ! using the base values for fractions of C and O

        if (frac_Type2 > 0d0 .and. &
            associated(kap_co_z_tables(rq% kap_CO_option)% ar)) then

           fC = kap_co_z_tables(rq% kap_CO_option)% ar(1)% Zfrac_C
           fN = kap_co_z_tables(rq% kap_CO_option)% ar(1)% Zfrac_N
           fO = kap_co_z_tables(rq% kap_CO_option)% ar(1)% Zfrac_O
           fNe = kap_co_z_tables(rq% kap_CO_option)% ar(1)% Zfrac_Ne

           dXC = max(0d0, xC - fC*Zbase)
           dXO = max(0d0, xO - fO*Zbase)

           if (X >= rq% kap_Type2_full_off_X) then
              frac_X = 0d0
           else if (X <= rq% kap_Type2_full_on_X) then
              frac_X = 1d0
           else ! blend
              frac_X = (rq% kap_Type2_full_off_X - X) / &
                   (rq% kap_Type2_full_off_X - rq% kap_Type2_full_on_X)
           end if

           dZ = Z - Zbase
           if (dZ <= rq% kap_Type2_full_off_dZ) then
              frac_dZ = 0d0
           else if (dZ >= rq% kap_Type2_full_on_dZ) then
              frac_dZ = 1d0
           else ! blend
              frac_dZ = (rq% kap_Type2_full_off_dZ - dZ) / &
                   (rq% kap_Type2_full_off_dZ - rq% kap_Type2_full_on_dZ)
           end if

           frac_Type2 = frac_Type2*frac_X*frac_dZ
           if (frac_Type2 == 0d0) then
              dXC = 0d0
              dXO = 0d0
           end if

           if (is_bad(dXC) .or. is_bad(dXO)) then
              write(*,1) 'X', X
              write(*,1) 'Z', Z
              write(*,1) 'Zbase', Zbase
              write(*,1) 'XC', XC
              write(*,1) 'XN', XN
              write(*,1) 'XO', XO
              write(*,1) 'XNe', XNe
              write(*,1) 'logRho', logRho
              write(*,1) 'logT', logT
              write(*,1) 'dXC', dXC
              write(*,1) 'dXO', dXO
              write(*,1) 'frac_X', frac_X
              write(*,1) 'frac_dZ', frac_dZ
              write(*,1) 'frac_Type2', frac_Type2
              ierr = -1
              return
           end if

        end if

        frac_Type1 = 1d0 - frac_Type2
        if (dbg) then
           write(*,1) 'max_frac_Type2', max_frac_Type2
           write(*,1) 'frac_Type2', frac_Type2
           write(*,1) 'frac_Type1', frac_Type1
           write(*,*)
           write(*,1) 'X', X
           write(*,1) 'dXC', dXC
           write(*,1) 'dXO', dXO
           write(*,1) 'Zbase', Zbase
        end if

        ! do blend
        beta = frac_Type1
        alfa = frac_Type2

        if (beta > 0) then ! get value from fixed metal tables
           if (rq% use_Type2_opacities .and. rq% use_Zbase_for_Type1) then
              ! when Z > Zbase, use Zbase instead of Z.
              call Get1_kap_fixed_metal_Results( &
                   kap_z_tables(rq% kap_option)% ar, num_kap_Zs(rq% kap_option), rq, min(Z,Zbase), X, Rho, logRho, T, logT, &
                   logkap_beta, dlnkap_beta_dlnRho, dlnkap_beta_dlnT, ierr)
           else
              call Get1_kap_fixed_metal_Results( &
                   kap_z_tables(rq% kap_option)% ar, num_kap_Zs(rq% kap_option), rq, Z, X, Rho, logRho, T, logT, &
                   logkap_beta, dlnkap_beta_dlnRho, dlnkap_beta_dlnT, ierr)
           end if
           kap_beta = exp10(logkap_beta)
           if (dbg) then
              write(*,1) 'Get_kap_fixed_metal_Results'
              write(*,1) 'Z', Z
              write(*,1) 'X', X
              write(*,1) 'Rho', Rho
              write(*,1) 'logRho', logRho
              write(*,1) 'T', T
              write(*,1) 'logT', logT
              write(*,1) 'kap_beta', kap_beta
              write(*,1) 'beta', beta
              write(*,*)
           end if
        else
           kap_beta = 0d0; dlnkap_beta_dlnRho = 0d0; dlnkap_beta_dlnT = 0d0
        end if

        if (alfa > 0d0) then ! get value from C/O enhanced tables
           call Get1_kap_CO_Results( &
                rq, Zbase, X, dXC, dXO, Rho, logRho, T, logT, &
                logkap_alfa, dlnkap_alfa_dlnRho, dlnkap_alfa_dlnT, ierr)
           kap_alfa = exp10(logkap_alfa)
           if (dbg) then
              write(*,1) 'Get_kap_CO_Results'
              write(*,1) 'Z', Z
              write(*,1) 'X', X
              write(*,1) 'Rho', Rho
              write(*,1) 'logRho', logRho
              write(*,1) 'T', T
              write(*,1) 'logT', logT
              write(*,1) 'kap_alfa', kap_alfa
              write(*,1) 'alfa', alfa
              write(*,*)
           end if
        else
           kap_alfa = 0d0; dlnkap_alfa_dlnRho = 0d0; dlnkap_alfa_dlnT = 0d0
        end if

        kap = alfa*kap_alfa + beta*kap_beta
        if (dbg) then
           write(*,1) 'alfa', alfa
           write(*,1) 'kap_alfa (Type2)', kap_alfa
           write(*,1) 'beta', beta
           write(*,1) 'kap_beta (Type1)', kap_beta
           write(*,1) 'kap', kap
        end if

        if (kap < 1d-30) then
           kap = 1d-30
           dlnkap_dlnRho = 0d0
           dlnkap_dlnT = 0d0
        else
           dlnkap_dlnRho = (alfa*kap_alfa*dlnkap_alfa_dlnRho + beta*kap_beta*dlnkap_beta_dlnRho) / kap
           dlnkap_dlnT = (alfa*kap_alfa*dlnkap_alfa_dlnT + beta*kap_beta*dlnkap_beta_dlnT) / kap
        end if

      end subroutine Get_kap_highT_Results


      subroutine combine_rad_with_conduction( &
            rq, Rho, logRho, T, logT, zbar, &
            kap_rad, dlnkap_rad_dlnRho, dlnkap_rad_dlnT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)

         use condint, only: do_electron_conduction
         type (Kap_General_Info), pointer :: rq
         real(dp), intent(in) :: Rho, logRho, T, logT, zbar
         real(dp), intent(inout) :: kap_rad, dlnkap_rad_dlnRho, dlnkap_rad_dlnT
         real(dp), intent(out) :: kap, dlnkap_dlnRho, dlnkap_dlnT
         integer, intent(out) :: ierr ! 0 means AOK.
      
         real(dp) :: kap_ec, dlnkap_ec_dlnRho, dlnkap_ec_dlnT
         logical, parameter :: dbg = .false.
         
         include 'formats'
         
         ierr = 0
         
         if (.not. rq% include_electron_conduction) then
            kap = kap_rad
            dlnkap_dlnRho = dlnkap_rad_dlnRho
            dlnkap_dlnT = dlnkap_rad_dlnT
            return
         end if

         if (rq% use_other_elect_cond_opacity) then
            call rq% other_elect_cond_opacity( &
               rq% handle, zbar, logRho, logT, &
               kap_ec, dlnkap_ec_dlnRho, dlnkap_ec_dlnT, ierr)
         else
            call do_electron_conduction( &
               rq, zbar, logRho, logT, &
               kap_ec, dlnkap_ec_dlnRho, dlnkap_ec_dlnT, ierr)
         end if
         if (ierr /= 0) return

         if (is_bad(kap_ec)) then
            write(*,*) 'kap_ec', kap_ec
            call mesa_error(__FILE__,__LINE__,'combine_rad_with_conduction')
         end if
         if (dbg) write(*,1) 'kap_ec', kap_ec
         if (dbg) write(*,1) 'dlnkap_ec_dlnRho', dlnkap_ec_dlnRho
         if (dbg) write(*,1) 'dlnkap_ec_dlnT', dlnkap_ec_dlnT
         if (dbg) write(*,*)
         
         kap = 1d0 / (1d0/kap_rad + 1d0/kap_ec)
         if (dbg) write(*,1) 'kap_rad', kap_rad
         if (dbg) write(*,1) 'kap', kap
         if (dbg) write(*,1) 'log10(kap)', log10(kap)
         
         if (is_bad(kap)) then
            ierr = -1; return
            write(*,1) 'kap', kap
            call mesa_error(__FILE__,__LINE__,'Get_kap_Results')
         end if

         dlnkap_dlnRho = (kap/kap_rad) * dlnkap_rad_dlnRho + (kap/kap_ec) * dlnkap_ec_dlnRho

         if (is_bad(dlnkap_dlnRho)) then
            ierr = -1; return
            write(*,1) 'dlnkap_dlnRho', dlnkap_dlnRho
            write(*,1) 'kap', kap
            write(*,1) 'dkap_dlnRho', kap * dlnkap_dlnRho
            write(*,1) 'dkap_ec_dlnRho', kap_ec * dlnkap_ec_dlnRho
            write(*,1) 'dkap_rad_dlnRho', kap_rad * dlnkap_rad_dlnRho
            write(*,1) 'kap_rad', kap_rad
            write(*,1) 'kap_ec', kap_ec
            call mesa_error(__FILE__,__LINE__,'combine_rad_with_conduction')
         end if
         
         dlnkap_dlnT = (kap/kap_rad) * dlnkap_rad_dlnT + (kap/kap_ec) * dlnkap_ec_dlnT
         
         if (is_bad(dlnkap_dlnT)) then
            ierr = -1; return
            write(*,1) 'dlnkap_dlnT', dlnkap_dlnT
            write(*,1) 'kap', kap
            write(*,1) 'dkap_dlnT', kap * dlnkap_dlnT
            write(*,1) 'dkap_ec_dlnT', kap_ec * dlnkap_ec_dlnT
            write(*,1) 'dkap_rad_dlnT', kap_rad * dlnkap_rad_dlnT
            write(*,1) 'kap_rad', kap_rad
            write(*,1) 'kap_ec', kap_ec
            call mesa_error(__FILE__,__LINE__,'combine_rad_with_conduction')
         end if

         if (dbg) write(*,1) 'dlnkap_dlnRho', dlnkap_dlnRho
         if (dbg) write(*,1) 'dlnkap_dlnT', dlnkap_dlnT
         if (dbg) call mesa_error(__FILE__,__LINE__,'combine_rad_with_conduction')
      
      
      end subroutine combine_rad_with_conduction


      subroutine Compton_Opacity(rq, &
            Rho_in, T_in, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            eta_in, d_eta_dlnRho, d_eta_dlnT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         use eos_lib
         use eos_def
         use const_def
         use auto_diff

         type (Kap_General_Info), pointer :: rq

         ! evaluates the poutanen 2017, apj 835, 119 fitting formula for the compton opacity.
         ! coefficients from table 1 between 2 and 300 kev

         real(dp), intent(in) :: Rho_in, T_in, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
         real(dp), intent(in) :: eta_in, d_eta_dlnRho, d_eta_dlnT
         real(dp), intent(out) :: kap, dlnkap_dlnRho, dlnkap_dlnT
         integer, intent(out) :: ierr
         
         type(auto_diff_real_2var_order1) :: T, rho, free_e, eta, kap_auto
         type(auto_diff_real_2var_order1) :: zeta, f1, f2, f3, alpha, tbr, theta, tkev, mfp

         real(dp) :: free_e_val

         ! poutanen table 1
         real(dp), parameter :: &
            t0     =  43.3d0, &
            alpha0 =  0.885d0, &
            c01    =  0.682d0, &
            c02    = -0.0454d0, &
            c11    =  0.240d0, &
            c12    =  0.0043d0, &
            c21    =  0.050d0, &
            c22    = -0.0067d0, &
            c31    = -0.037d0, &
            c32    =  0.0031d0
         
         include 'formats'
         
         ierr = 0

         ! set up auto diff
         ! var1: Rho
         ! var2: T
         
         Rho = Rho_in
         Rho% d1val1 = 1d0
         Rho% d1val2 = 0d0

         T = T_in
         T% d1val1 = 0d0
         T% d1val2 = 1d0

         free_e_val = exp(lnfree_e)
         free_e = free_e_val
         free_e % d1val1 = free_e_val*d_lnfree_e_dlnRho/Rho_in
         free_e % d1val2 = free_e_val*d_lnfree_e_dlnT/T_in

         eta = eta_in
         eta % d1val1 = d_eta_dlnRho/Rho_in
         eta % d1val2 = d_eta_dlnT/T_in

         theta = T * kerg / (me * clight * clight)
         tkev  = T * kev * 1.0d-3
         zeta  = exp(c01 * eta + c02 * eta*eta)
         f1    = 1.0d0 + c11 * zeta + c12 * zeta * zeta
         f2    = 1.0d0 + c21 * zeta + c22 * zeta * zeta
         f3    = 1.0d0 + c31 * zeta + c32 * zeta * zeta
         alpha = alpha0 * f3
         tbr   = t0 * f2
         mfp   = f1 * (1.0d0 + pow(tkev/tbr,alpha))

         ! equation 31
         kap_auto = sige*(free_e/amu)/mfp

         ! unpack auto_diff
         kap = kap_auto% val
         dlnkap_dlnRho = Rho% val * kap_auto% d1val1 / kap
         dlnkap_dlnT = T% val * kap_auto% d1val2 / kap
         
      end subroutine Compton_Opacity


      end module kap_eval
      
