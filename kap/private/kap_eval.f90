! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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
           rq, zbar, X, Z, XC, XN, XO, XNe, logRho_in, logT_in, &
           lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
           frac_Type2, kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         use const_def
         use utils_lib, only: is_bad

         type (Kap_General_Info), pointer :: rq
         real(dp), intent(in) :: X, XC, XN, XO, XNe, Z, zbar
         real(dp), intent(in) :: logRho_in, logT_in
         real(dp), intent(in) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
         real(dp), intent(out) :: frac_Type2
         real(dp), intent(out) :: kap ! opacity
         real(dp), intent(out) :: dlnkap_dlnRho ! partial derivative at constant T
         real(dp), intent(out) :: dlnkap_dlnT   ! partial derivative at constant Rho
         integer, intent(out) :: ierr ! 0 means AOK.
         
         real(dp) :: logR, logT, T, logRho, Rho, logKap
         real(dp) :: logT_Compton_blend_hi, logR_Compton_blend_lo, &
            kap_rad, dlnkap_rad_dlnRho, dlnkap_rad_dlnT, &
            beta, kap_beta, dlnkap_beta_dlnRho, dlnkap_beta_dlnT, &
            alfa, kap_alfa, dlnkap_alfa_dlnRho, dlnkap_alfa_dlnT
         
         logical :: dbg

         include 'formats.dek'

         dbg = rq% dbg
         if (dbg) dbg = & ! check limits
            logT >= rq% logT_lo .and. logT <= rq% logT_hi .and. &
            logRho >= rq% logRho_lo .and. logRho <= rq% logRho_hi .and. &
            X >= rq% X_lo .and. X <= rq% X_hi .and. &
            Z >= rq% Z_lo .and. Z <= rq% Z_hi

         logRho = logRho_in; Rho = exp10(logRho)
         logT = logT_in; T = exp10(logT)
         logR = logRho - 3*logT + 18
         if (dbg) write(*,1) 'logRho', logRho
         if (dbg) write(*,1) 'logT', logT
         if (dbg) write(*,1) 'logR', logR
         
         frac_Type2 = 0d0


         logT_Compton_blend_hi = rq% logT_Compton_blend_hi
         logR_Compton_blend_lo = rq% logR_Compton_blend_lo
         if (logT >= logT_Compton_blend_hi .or. logR <= logR_Compton_blend_lo) then ! just use compton

            kap_rad = 1d-30
            dlnkap_rad_dlnRho = 0d0
            dlnkap_rad_dlnT = 0d0

         else

            call Get_kap_Results_blend_T( &
                 rq, X, Z, XC, XN, XO, XNe, logRho, logT, &
                 frac_Type2, kap_rad, dlnkap_rad_dlnRho, dlnkap_rad_dlnT, ierr)
                 
         end if
         
         call combine_rad_with_compton_and_conduction( &
            rq, Rho, logRho, T, logT, zbar, &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            kap_rad, dlnkap_rad_dlnRho, dlnkap_rad_dlnT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         
      end subroutine Get_kap_Results
      
      
      subroutine Get_kap_Results_blend_T( &
           rq, X, Z, XC, XN, XO, XNe, logRho_in, logT_in, &
           frac_Type2, kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)

        use kap_def
        use utils_lib, only: is_bad

        ! INPUT
        type (Kap_General_Info), pointer :: rq
        real(dp), intent(in) :: X, Z, XC, XN, XO, XNe ! composition    
        real(dp), intent(in) :: logRho_in ! density
        real(dp), intent(in) :: logT_in ! temperature

        ! OUTPUT
        real(dp), intent(out) :: frac_Type2
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


      subroutine combine_rad_with_compton_and_conduction( &
            rq, Rho, logRho, T, logT, zbar, &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            kap_rad, dlnkap_rad_dlnRho, dlnkap_rad_dlnT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)

         use condint, only: do_electron_conduction
         type (Kap_General_Info), pointer :: rq
         real(dp), intent(in) :: Rho, logRho, T, logT, zbar
         real(dp), intent(in) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
         real(dp), intent(inout) :: kap_rad, dlnkap_rad_dlnRho, dlnkap_rad_dlnT
         real(dp), intent(out) :: kap, dlnkap_dlnRho, dlnkap_dlnT
         integer, intent(out) :: ierr ! 0 means AOK.
      
         real(dp) :: &
            kap_ec, dlnkap_ec_dlnRho, dlnkap_ec_dlnT, &
            kap_compton, dlnkap_compton_dlnRho, dlnkap_compton_dlnT, kap_highT, &
            dkap_compton_dlnRho, dkap_compton_dlnT, &
            alfa, d_alfa_dlnT, d_alfa_dlnRho, beta, d_beta_dlnT, d_beta_dlnRho,&
            alfa0, d_alfa0_dlnT, d_alfa0_dlnRho, logKap, logKap_rad, logKap_compton
         
         real(dp) :: logT_Compton_blend_lo, logT_Compton_blend_hi
         real(dp) :: logR, logR_Compton_blend_lo, logR_Compton_blend_hi
         
         logical, parameter :: dbg = .false.
         
         include 'formats'
         
         ierr = 0
         
         if (dbg) then
            write(*,1) 'logT', logT
            write(*,1) 'logRho', logRho
            write(*,1) 'lnfree_e', lnfree_e
            write(*,1) 'd_lnfree_e_dlnRho', d_lnfree_e_dlnRho
            write(*,1) 'd_lnfree_e_dlnT', d_lnfree_e_dlnT
            call test_Compton_Opacity( &
               Rho, T, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, ierr)
            stop 'test_Compton_Opacity'
         end if

         logT_Compton_blend_hi = rq% logT_Compton_blend_hi
         logT_Compton_blend_lo = logT_Compton_blend_hi - 0.50d0

         !also blend to Compton at low R
         logR = logRho - 3d0*logT + 18d0
         logR_Compton_blend_lo = rq% logR_Compton_blend_lo
         logR_Compton_blend_hi = logR_Compton_blend_lo + 0.50d0

         if (logT_Compton_blend_hi < 7.5d0) then
            write(*,1) 'Get_kap_Results: logT_Compton_blend_hi', logT_Compton_blend_hi
            write(*,1) 'bogus??? expect 7.99 for OP and 8.69 for OPAL'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         if (logT > logT_Compton_blend_lo .or. logR < logR_Compton_blend_hi) then ! combine kap_rad with rad_compton
         
            if (dbg) write(*,*) 'combine kap_rad with rad_compton'
            
            call Compton_Opacity( &
               Rho, T, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
               kap_compton, dlnkap_compton_dlnRho, dlnkap_compton_dlnT, ierr)
            if (ierr /= 0) return
            if (dbg) write(*,1) 'kap_compton', kap_compton
            if (dbg) write(*,1) 'dlnkap_compton_dlnRho', dlnkap_compton_dlnRho
            if (dbg) write(*,1) 'dlnkap_compton_dlnT', dlnkap_compton_dlnT
            if (dbg) write(*,*) 'logT >= logT_Compton_blend_hi', logT >= logT_Compton_blend_hi, &
               logT, logT_Compton_blend_hi
            if (logT >= logT_Compton_blend_hi .or. logR <= logR_Compton_blend_lo) then
               kap_rad = kap_compton
               dlnkap_rad_dlnRho = dlnkap_compton_dlnRho
               dlnkap_rad_dlnT = dlnkap_compton_dlnT
            else
               !blend in temperature
               if (logT > logT_Compton_blend_lo .and. logT <= logT_Compton_blend_hi) then
                  alfa0 = (logT - logT_Compton_blend_lo) &
                     / (logT_Compton_blend_hi - logT_Compton_blend_lo)
                  d_alfa0_dlnT = 1d0/(logT_Compton_blend_hi - logT_Compton_blend_lo)/ln10
                  
                  ! must smooth the transitions near alfa = 0.0 and 1.0  
                  ! use quintic smoothing function for this with 1st and 2nd derivs = 0 at ends
                  alfa = -alfa0*alfa0*alfa0*(-10d0 + alfa0*(15d0 - 6d0*alfa0))
                  d_alfa_dlnT = 30d0*(alfa0 - 1d0)*(alfa0 - 1d0)*alfa0*alfa0*d_alfa0_dlnT

                  beta = 1d0 - alfa
                  d_beta_dlnT = -d_alfa_dlnT
                  
                  if (dbg) then
                     write(*,1) 'alfa', alfa
                     write(*,1) 'beta', beta
                  end if

                  logKap_rad = log10(kap_rad)
                  logKap_compton = log10(kap_compton)
                  
                  logKap = alfa*logKap_compton + beta*logKap_rad
                  dlnkap_rad_dlnRho = &
                     alfa*dlnkap_compton_dlnRho + beta*dlnkap_rad_dlnRho
                  dlnkap_rad_dlnT = &
                     alfa*dlnkap_compton_dlnT + beta*dlnkap_rad_dlnT + &
                     d_alfa_dlnT*logKap_compton*ln10 + d_beta_dlnT*logKap_rad*ln10
                  kap_rad = exp10(logKap)
               end if
               !blend in R
               if (logR > logR_Compton_blend_lo .and. logR <= logR_Compton_blend_hi) then
                  alfa0 = (logR - logR_Compton_blend_lo) &
                     / (logR_Compton_blend_hi - logR_Compton_blend_lo)
                  d_alfa0_dlnT = -3d0/(logR_Compton_blend_hi - logR_Compton_blend_lo)/ln10
                  d_alfa0_dlnRho = 1d0/(logR_Compton_blend_hi - logR_Compton_blend_lo)/ln10
                  
                  ! must smooth the transitions near alfa = 0.0 and 1.0  
                  ! use quintic smoothing function for this with 1st and 2nd derivs = 0 at ends
                  alfa = -alfa0*alfa0*alfa0*(-10d0 + alfa0*(15d0 - 6d0*alfa0))
                  d_alfa_dlnT = 30d0*(alfa0 - 1d0)*(alfa0 - 1d0)*alfa0*alfa0*d_alfa0_dlnT
                  d_alfa_dlnRho = 30d0*(alfa0 - 1d0)*(alfa0 - 1d0)*alfa0*alfa0*d_alfa0_dlnRho

                  beta = 1d0 - alfa
                  d_beta_dlnT = -d_alfa_dlnT
                  d_beta_dlnRho = -d_alfa_dlnRho
                  
                  if (dbg) then
                     write(*,1) 'alfa', alfa
                     write(*,1) 'beta', beta
                  end if

                  logKap_rad = log10(kap_rad)
                  logKap_compton = log10(kap_compton)
                  
                  logKap = beta*logKap_compton + alfa*logKap_rad
                  dlnkap_rad_dlnRho = &
                     beta*dlnkap_compton_dlnRho + alfa*dlnkap_rad_dlnRho + &
                     d_beta_dlnRho*logKap_compton*ln10 + d_alfa_dlnRho*logKap_rad*ln10
                  dlnkap_rad_dlnT = &
                     beta*dlnkap_compton_dlnT + alfa*dlnkap_rad_dlnT + &
                     d_beta_dlnT*logKap_compton*ln10 + d_alfa_dlnT*logKap_rad*ln10
                  kap_rad = exp10(logKap)
               end if
            end if
            if (dbg) write(*,1) 'kap_rad', kap_rad
            if (dbg) write(*,1) 'dlnkap_rad_dlnRho', dlnkap_rad_dlnRho
            if (dbg) write(*,1) 'dlnkap_rad_dlnT', dlnkap_rad_dlnT
            if (dbg) write(*,*)
         end if
         
         if (.not. rq% include_electron_conduction) then
            kap = kap_rad
            dlnkap_dlnRho = dlnkap_rad_dlnRho
            dlnkap_dlnT = dlnkap_rad_dlnT
            return
         end if
         
         call do_electron_conduction( &
            zbar, logRho, logT, &
            kap_ec, dlnkap_ec_dlnRho, dlnkap_ec_dlnT, ierr)
         if (ierr /= 0) return

         if (is_bad(kap_ec)) then
            write(*,*) 'kap_ec', kap_ec
            stop 'Get_kap_Results'
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
            stop 'Get_kap_Results'
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
            stop 'Get_kap_Results'
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
            stop 'Get_kap_Results'
         end if

         if (dbg) write(*,1) 'dlnkap_dlnRho', dlnkap_dlnRho
         if (dbg) write(*,1) 'dlnkap_dlnT', dlnkap_dlnT
         if (dbg) stop 'combine_rad_with_compton_and_conduction'
      
      
      end subroutine combine_rad_with_compton_and_conduction

      subroutine test_Compton_Opacity( &
            Rho, T, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, ierr)
         real(dp), intent(in) :: Rho, T, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
         integer, intent(out) :: ierr
         
         logical :: doing_d_dlnd
         real(dp) :: kap, dlnkap_dlnRho, dlnkap_dlnT, &
            dbg_var, d_dbg_var_dlnd, d_dbg_var_dlnT, &
            lnd, lnT, dx_0, dvardx_0, err, dvardx, xdum
         include 'formats'
         ierr = 0
         
         call Compton_Opacity_test( &
            Rho, T, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, &
            dbg_var, d_dbg_var_dlnd, d_dbg_var_dlnT, ierr)
         if (ierr /= 0) stop 'failed in test_Compton_Opacity'
         
         !doing_d_dlnd = .true.
         doing_d_dlnd = .false.
      
         lnd = log(Rho)
         lnT = log(T)

         if (doing_d_dlnd) then
            dx_0 = lnd*1d-4
            dvardx_0 = d_dbg_var_dlnd ! analytic value of partial
         else
            dx_0 = lnT*1d-4
            dvardx_0 = d_dbg_var_dlnT ! analytic value of partial
         end if
         err = 0d0
         dvardx = dfridr(dx_0,err)
         xdum = (dvardx - dvardx_0)/max(abs(dvardx_0),1d-3)         
         if (doing_d_dlnd) then
            write(*,1) 'compton d_dlnRho analytic, numeric, diff, rel diff', &
                  dvardx_0, dvardx, err, xdum
         else ! doing d_dlnT
            write(*,1) 'compton d_dlnT analytic, numeric, diff, rel diff', &
                  dvardx_0, dvardx, err, xdum
         end if
         write(*,*)
         
      
         contains
         
         
         real(dp) function dfridr_func(delta_x) result(val)
            real(dp), intent(in) :: delta_x
            integer :: ierr
            real(dp) :: log_var
            include 'formats'
            ierr = 0
            if (doing_d_dlnd) then
               log_var = (lnd + delta_x)/ln10
               call Compton_Opacity_test( &
                  exp10(log_var), T, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
                  kap, dlnkap_dlnRho, dlnkap_dlnT, &
                  dbg_var, d_dbg_var_dlnd, d_dbg_var_dlnT, ierr)
            else
               log_var = (lnT + delta_x)/ln10
               call Compton_Opacity_test( &
                  Rho, exp10(log_var), lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
                  kap, dlnkap_dlnRho, dlnkap_dlnT, &
                  dbg_var, d_dbg_var_dlnd, d_dbg_var_dlnT, ierr)
            end if
            val = dbg_var
         end function dfridr_func

         real(dp) function dfridr(hx,err) ! from Frank
            real(dp), intent(in) :: hx
            real(dp), intent(out) :: err
            !  this routine returns the first derivative of a function func(x)
            !  at the point x, by ridders method of polynomial extrapolation.
            !  value hx is the initial step size;
            !  it should be an increment for which func changes substantially.
            !  an estimate of the error in the first derivative is returned in err.
            integer, parameter :: ntab = 20
            integer :: i,j
            real(dp) :: x,errt,fac,hh,a(ntab,ntab),xdum,ydum
            real(dp), parameter :: con2=2d0, con=sqrt(con2), big=1d50, safe=2d0
            include 'formats'
            dfridr = 0d0
            hh = hx
            ! 2nd order central difference
            a(1,1) = (dfridr_func(hh) - dfridr_func(-hh))/(2d0*hh)
            write(*,2) 'dfdx hh', 1, a(1,1), hh
            err = big
            ! succesive columns in the neville tableu will go to smaller stepsizes
            ! and higher orders of extrapolation
            do i=2,ntab
               hh = hh/con
               a(1,i) = (dfridr_func(hh) - dfridr_func(-hh))/(2d0*hh)
               !write(*,2) 'dfdx hh', i, a(1,i), hh
               ! compute extrapolations of various orders; the error stratagy is to compare
               ! each new extrapolation to one order lower but both at the same stepsize
               ! and at the previous stepsize
               fac = con2
               do j=2,i
                  a(j,i) = (a(j-1,i)*fac - a(j-1,i-1))/(fac-1d0)
                  fac = con2*fac
                  errt = max(abs(a(j,i)-a(j-1,i)),abs(a(j,i)-a(j-1,i-1)))
                  if (errt <= err) then
                     err = errt
                     dfridr = a(j,i)
                     write(*,3) 'dfridr err', i, j, dfridr, err
                  end if
               end do
               ! if higher order is worse by a significant factor safe, then bail
               if (abs(a(i,i) - a(i-1,i-1)) >= safe*err) then
                  write(*,1) 'higher order is worse', err, a(i,i), a(i-1,i-1)
                  return
               end if
            end do
         end function dfridr
         
      end subroutine test_Compton_Opacity
      

      subroutine Compton_Opacity( &
            Rho_in, T_in, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         real(dp), intent(in) :: Rho_in, T_in, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
         real(dp), intent(out) :: kap, dlnkap_dlnRho, dlnkap_dlnT
         integer, intent(out) :: ierr
         real(dp) :: dbg_var, d_dbg_var_dlnd, d_dbg_var_dlnT
         call Compton_Opacity_test( &
            Rho_in, T_in, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, &
            dbg_var, d_dbg_var_dlnd, d_dbg_var_dlnT, ierr)
      end subroutine Compton_Opacity
      

      subroutine Compton_Opacity_test( &
            Rho_in, T_in, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, &
            dbg_var, d_dbg_var_dlnd, d_dbg_var_dlnT, ierr)
         use eos_lib
         use eos_def
         use const_def

         ! same as Thompson at low T and low RHO.
         ! high T and high degeneracy both lengthen mean free path and thus decrease opacity.
         ! formula for approximating this effect from Buchler & Yueh, 1976, Apj, 210:440-446.

         ! this approximation breaks down for high degeneracy (eta > 4), but
         ! by then conduction should be dominant anyway, so things should be okay.

         real(dp), intent(in) :: Rho_in, T_in, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
         real(dp), intent(out) :: kap, dlnkap_dlnRho, dlnkap_dlnT, &
            dbg_var, d_dbg_var_dlnd, d_dbg_var_dlnT
         integer, intent(out) :: ierr
         
         real(dp) :: T, rho, edensity, efermi, eta, psi, psi2, Gbar_inverse, kT_mc2, free_e, tmp, &
            d_free_e_dlnRho, d_free_e_dlnT, d_edensity_dlnRho, d_edensity_dlnT, &
            edensity_one_third, d_efermi_dlnRho, d_efermi_dlnT, d_kT_mc2_dlnT, &
            d_eta_dlnRho, d_eta_dlnT, &
            d_psi_dlnRho, d_psi_dlnT, d_psi2_dlnRho, d_psi2_dlnT, &
            d_Gbar_inverse_dlnRho, d_Gbar_inverse_dlnT
         real(dp), parameter :: rho_limit = 5d7
         real(dp), parameter :: T_limit = 2d10
         real(dp), parameter :: pow_3o8pi = 0.242430689510993d0 ! pow(3/(8*pi),2d0/3d0)
         
         logical :: hit_rho_limit, hit_T_limit
         
         include 'formats.dek'
         
         
         
         
         
         ! NOTE: to get correct numerical estimates need to call eos to get new free_e when change lnd/lnT
         
         
         
         
         
         ierr = 0
         free_e = exp(lnfree_e)
         d_free_e_dlnRho = free_e*d_lnfree_e_dlnRho
         d_free_e_dlnT = free_e*d_lnfree_e_dlnT
         
         T = T_in
         hit_T_limit = (T > T_limit)
         if (hit_T_limit) T = T_limit
         rho = rho_in
         hit_rho_limit = (rho > rho_limit)
         if (hit_rho_limit) rho = rho_limit

         edensity = free_e*rho/amu
         if (edensity <= 0) then
            ierr = -1; return
         end if
         d_edensity_dlnRho = d_free_e_dlnRho*rho/amu + edensity
         d_edensity_dlnT = d_free_e_dlnT*rho/amu
         
         tmp = planck_h*planck_h/(2d0*me)*pow_3o8pi
         edensity_one_third = pow(edensity,1d0/3d0)
         efermi = tmp*edensity_one_third*edensity_one_third
         d_efermi_dlnRho = (2d0/3d0)*tmp*d_edensity_dlnRho/edensity_one_third
         d_efermi_dlnT = (2d0/3d0)*tmp*d_edensity_dlnT/edensity_one_third
         
         dbg_var = efermi
         d_dbg_var_dlnd = d_efermi_dlnRho
         d_dbg_var_dlnT = d_efermi_dlnT
         
         kT_mc2 = T * (kerg / (me * clight*clight))
         d_kT_mc2_dlnT = kT_mc2

         eta = (efermi - me*clight*clight)/(kerg*T)
         d_eta_dlnRho = d_efermi_dlnRho/(kerg*T)
         d_eta_dlnT = d_efermi_dlnT/(kerg*T) - eta

         if (eta > 4d0) then
            eta = 4d0
            d_eta_dlnRho = 0d0
            d_eta_dlnT = 0d0
         else if (eta < -10d0) then  
            eta = -10
            d_eta_dlnRho = 0d0
            d_eta_dlnT = 0d0
         end if

         psi = exp(eta*(0.8168d0 - eta*0.05522d0)) ! coefficients from B&Y
         d_psi_dlnRho = psi*(0.8168d0 - 2*0.05522d0*eta)*d_eta_dlnRho
         d_psi_dlnT = psi*(0.8168d0 - 2*0.05522d0*eta)*d_eta_dlnT
         
         psi2 = psi*psi
         d_psi2_dlnRho = 2d0*psi*d_psi_dlnRho
         d_psi2_dlnT = 2d0*psi*d_psi_dlnT

         ! formula for Gbar_inverse from B&Y (incorrectly given as Gbar in their paper!)
         Gbar_inverse = &
            1.129d0+0.2965d0*psi-0.005594d0*psi2 &
            +(11.47d0+0.3570d0*psi+0.1078d0*psi2)*kT_mc2 &
            +(0.1678d0*psi-3.249d0-0.04706d0*psi2)*kT_mc2*kT_mc2
            
         d_Gbar_inverse_dlnRho = &
            0.2965d0*d_psi_dlnRho-0.005594d0*d_psi2_dlnRho &
            +(0.3570d0*d_psi_dlnRho+0.1078d0*d_psi2_dlnRho)*kT_mc2 &
            +(0.1678d0*d_psi_dlnRho-0.04706d0*d_psi2_dlnRho)*kT_mc2*kT_mc2
            
         d_Gbar_inverse_dlnT = &
            0.2965d0*d_psi_dlnT-0.005594d0*d_psi2_dlnT &
            +(0.3570d0*d_psi_dlnT+0.1078d0*d_psi2_dlnT)*kT_mc2 &
            +(0.1678d0*d_psi_dlnT-0.04706d0*d_psi2_dlnT)*kT_mc2*kT_mc2 &
            +(11.47d0+0.3570d0*psi+0.1078d0*psi2)*d_kT_mc2_dlnT &
            +(0.1678d0*psi-3.249d0-0.04706d0*psi2)*2d0*kT_mc2*d_kT_mc2_dlnT
         
         !dbg_var = Gbar_inverse
         !d_dbg_var_dlnd = d_Gbar_inverse_dlnRho
         !d_dbg_var_dlnT = d_Gbar_inverse_dlnT
         
         kap = sige*edensity/(rho*Gbar_inverse)
         
         if (kap <= 0.0d0) then
            write(*,*) 'lnfree_e', lnfree_e
            write(*,*) 'Rho', Rho
            write(*,*) 'T', T
            write(*,*) 'edensity', edensity
            write(*,*) 'efermi', efermi
            write(*,*) 'sige', sige
            write(*,*) 'Gbar_inverse', Gbar_inverse
            write(*,*) 'eta',eta
            write(*,*) 'psi', psi
            write(*,*) 'bad compton kap', kap
            call mesa_error(__FILE__,__LINE__)
         end if
         
         dlnkap_dlnRho = &
            d_edensity_dlnRho/edensity - d_Gbar_inverse_dlnRho/Gbar_inverse - 1d0
            
         dlnkap_dlnT = &
            d_edensity_dlnT/edensity - d_Gbar_inverse_dlnT/Gbar_inverse
            
         if (hit_rho_limit) dlnkap_dlnRho = 0d0
         if (hit_T_limit) dlnkap_dlnT = 0d0
         
         !dbg_var = log(kap)
         !d_dbg_var_dlnd = dlnkap_dlnRho
         !d_dbg_var_dlnT = dlnkap_dlnT

      end subroutine Compton_Opacity_test


      end module kap_eval
      
