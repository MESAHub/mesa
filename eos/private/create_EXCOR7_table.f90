! ***********************************************************************
!
!   Copyright (C) 2020  The MESA Team
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

module create_EXCOR7_table
   
   use const_def
   use chem_def
   use utils_lib, only : is_bad
   use math_lib
   
   implicit none
   
   public :: do_create_EXCOR7_table
   private


contains
   
   subroutine do_create_EXCOR7_table(fname)
      character (len = *), intent(in) :: fname
      real(dp) :: logRS, logRS_min, logRS_max, dlogRS
      real(dp) :: logGAME, logGAME_min, logGAME_max, dlogGAME
      integer :: nlogRS, nlogGAME, io_unit, i, j
      real(dp) :: RS, GAME, FXC, UXC, PXC, CVXC, SXC, PDTXC, PDRXC
      include 'formats'
      
      logRS_min = -3.5d0
      logRS_max = 0.0d0
      dlogRS = 1d-2
      nlogRS = (logRS_max - logRS_min) / dlogRS + 1
      
      logGAME_min = -2d0
      logGAME_max = 4d0
      dlogGAME = 1d-2
      nlogGAME = (logGAME_max - logGAME_min) / dlogGAME + 1
      
      !write(*,'(a)') 'create ' // trim(fname)
      
      open(newunit = io_unit, file = trim(fname))
      
      write(io_unit, '(99(a14))') 'num logRS', 'logRS min', 'logRS max', 'del logRS', &
         'num logGAME', 'logGAME min', 'logGAME max', 'del logGAME'
      
      write(io_unit, '(2(i10,4x,3(f14.4)),i10)') &
         nlogRS, logRS_min, logRS_max, dlogRS, &
         nlogGAME, logGAME_min, logGAME_max, dlogGAME, 0
      
      do i = 1, nlogRS
         logRS = logRS_min + (i - 1) * dlogRS
         RS = exp10(logRS)
         write(io_unit, '(/,7x,a)') 'logRS'
         write(io_unit, '(2x,f14.6/)') logRS
         write(io_unit, '(99(a26))') &
            'logGAME', 'FXC', 'UXC', 'PXC', 'CVXC', 'SXC', 'PDTXC', 'PDRXC', 'RS', 'GAME'
         do j = 1, nlogGAME
            logGAME = logGAME_min + (j - 1) * dlogGAME
            GAME = exp10(logGAME)
            if (is_bad(GAME)) then
               write(*, 3) 'GAME', i, j, GAME, logGAME, logGAME_min, dlogGAME
            end if
            call EXCOR7(RS, GAME, FXC, UXC, PXC, CVXC, SXC, PDTXC, PDRXC)
            write(io_unit, '(99(1pd26.16))') &
               logGAME, FXC, UXC, PXC, CVXC, SXC, PDTXC, PDRXC, RS, GAME
         end do
         write(io_unit, *)
      end do
      
      write(io_unit, *)
      write(io_unit, *)
      close(io_unit)
   
   end subroutine do_create_EXCOR7_table
   
   
   ! ==============  ELECTRON EXCHANGE AND CORRELATION   ================ *
   subroutine EXCOR7(RS, GAME, FXC, UXC, PXC, CVXC, SXC, PDTXC, PDRXC)
      !                                                       Version 09.06.07
      ! Exchange-correlation contribution for the electron gas
      ! Stems from TANAKA1 v.03.03.96. Added derivatives.
      ! Input: RS - electron density parameter =electron-sphere radius in a.u.
      !        GAME - electron Coulomb coupling parameter
      ! Output: FXC - excess free energy of e-liquid per kT per one electron
      !             according to Tanaka & Ichimaru 85-87 and Ichimaru 93
      !         UXC - internal energy contr.[per 1 electron, kT]
      !         PXC - pressure contribution divided by (n_e kT)
      !         CVXC - heat capacity divided by N_e k
      !         SXC - entropy divided by N_e k
      !         PDTXC,PDRXC = PXC + d ln PXC / d ln(T,\rho)
      real(dp), intent(in) :: RS, GAME
      real(dp), intent(out) :: FXC, UXC, PXC, CVXC, SXC, PDTXC, PDRXC
      
      real(dp) :: THETA, SQTH, THETA2, THETA3, THETA4, EXP1TH
      real(dp) :: CHT1, SHT1, CHT2, SHT2
      real(dp) :: T1, T2, T1DH, T1DHH, T2DH, T2DHH
      real(dp) :: A0, A0DH, A0DHH, A1, A1DH, A1DHH, A, AH, ADH, ADHH
      real(dp) :: B0, B0DH, B0DHH, B1, B1DH, B1DHH, B, BH, BDH, BDHH
      real(dp) :: C, CDH, CDHH, C3, C3DH, C3DHH
      real(dp) :: D0, D0DH, D0DHH, D1, D1DH, D1DHH, D, DH, DDH, DDHH
      real(dp) :: E0, E0DH, E0DHH, E1, E1DH, E1DHH, E, EH, EDH, EDHH
      real(dp) :: DISCR, DIDH, DIDHH
      real(dp) :: S1, S1H, S1DH, S1DHH, S1DG, S1DHG
      real(dp) :: B2, B2DH, B2DHH, SQGE, B3, B3DH, B3DHH
      real(dp) :: S2, S2H, S2DH, S2DHH, S2DG, S2DGG, S2DHG
      real(dp) :: R3, R3H, R3DH, R3DHH, R3DG, R3DGG, R3DHG
      real(dp) :: S3, S3DH, S3DHH, S3DG, S3DGG, S3DHG
      real(dp) :: B4, B4DH, B4DHH
      real(dp) :: C4, C4DH, C4DHH, C4DG, C4DGG, C4DHG
      real(dp) :: S4A, S4AH, S4ADH, S4ADHH, S4ADG, S4ADGG, S4ADHG
      real(dp) :: S4B, S4BDH, S4BDHH, UP1, DN1, UP2, DN2
      real(dp) :: S4C, S4CDH, S4CDHH, S4CDG, S4CDGG, S4CDHG
      real(dp) :: S4, S4DH, S4DHH, S4DG, S4DGG, S4DHG
      real(dp) :: FXCDH, FXCDHH, FXCDG, FXCDGG, FXCDHG
      real(dp) :: PDLH, PDLG
      include 'formats'
      THETA = 0.543d0 * RS / GAME ! non-relativistic degeneracy parameter
      SQTH = sqrt(THETA)
      THETA2 = THETA * THETA
      THETA3 = THETA2 * THETA
      THETA4 = THETA3 * THETA
      if (THETA.gt..005d0) then
         CHT1 = cosh(1.d0 / THETA)
         SHT1 = sinh(1.d0 / THETA)
         CHT2 = cosh(1.d0 / SQTH)
         SHT2 = sinh(1.d0 / SQTH)
         T1 = SHT1 / CHT1 ! dtanh(1.d0/THETA)
         T2 = SHT2 / CHT2 ! dtanh(1./sqrt(THETA))
         T1DH = -1.d0 / ((THETA * CHT1) * (THETA * CHT1)) ! d T1 / d\theta
         T1DHH = 2.d0 / pow3(THETA * CHT1) * (CHT1 - SHT1 / THETA)
         T2DH = -0.5d0 * SQTH / ((THETA * CHT2) * (THETA * CHT2))
         T2DHH = (0.75d0 * SQTH * CHT2 - 0.5d0 * SHT2) / pow3(THETA * CHT2)
      else
         T1 = 1.d0
         T2 = 1.d0
         T1DH = 0.d0
         T2DH = 0.d0
         T1DHH = 0.d0
         T2DHH = 0.d0
      endif
      A0 = 0.75d0 + 3.04363d0 * THETA2 - 0.09227d0 * THETA3 + 1.7035d0 * THETA4
      A0DH = 6.08726d0 * THETA - 0.27681d0 * THETA2 + 6.814d0 * THETA3
      A0DHH = 6.08726d0 - 0.55362d0 * THETA + 20.442d0 * THETA2
      A1 = 1d0 + 8.31051d0 * THETA2 + 5.1105d0 * THETA4
      A1DH = 16.62102d0 * THETA + 20.442d0 * THETA3
      A1DHH = 16.62102d0 + 61.326d0 * THETA2
      A = 0.610887d0 * A0 / A1 * T1 ! HF fit of Perrot and Dharma-wardana
      AH = A0DH / A0 - A1DH / A1 + T1DH / T1
      ADH = A * AH
      ADHH = ADH * AH + A * (A0DHH / A0 - pow2(A0DH / A0) - A1DHH / A1 + pow2(A1DH / A1) &
         + T1DHH / T1 - pow2(T1DH / T1))
      B0 = 0.341308d0 + 12.070873d0 * THETA2 + 1.148889d0 * THETA4
      B0DH = 24.141746d0 * THETA + 4.595556d0 * THETA3
      B0DHH = 24.141746d0 + 13.786668d0 * THETA2
      B1 = 1d0 + 10.495346d0 * THETA2 + 1.326623d0 * THETA4
      B1DH = 20.990692d0 * THETA + 5.306492d0 * THETA3
      B1DHH = 20.990692d0 + 15.919476d0 * THETA2
      B = SQTH * T2 * B0 / B1
      BH = 0.5d0 / THETA + T2DH / T2 + B0DH / B0 - B1DH / B1
      BDH = B * BH
      BDHH = BDH * BH + B * (-0.5d0 / THETA2 + T2DHH / T2 - pow2(T2DH / T2) &
         + B0DHH / B0 - pow2(B0DH / B0) - B1DHH / B1 + pow2(B1DH / B1))
      D0 = 0.614925d0 + 16.996055d0 * THETA2 + 1.489056d0 * THETA4
      D0DH = 33.99211d0 * THETA + 5.956224d0 * THETA3
      D0DHH = 33.99211d0 + 17.868672d0 * THETA2
      D1 = 1d0 + 10.10935d0 * THETA2 + 1.22184d0 * THETA4
      D1DH = 20.2187d0 * THETA + 4.88736d0 * THETA3
      D1DHH = 20.2187d0 + 14.66208d0 * THETA2
      D = SQTH * T2 * D0 / D1
      DH = 0.5d0 / THETA + T2DH / T2 + D0DH / D0 - D1DH / D1
      DDH = D * DH
      DDHH = DDH * DH + D * (-0.5d0 / THETA2 + T2DHH / T2 - pow2(T2DH / T2) &
         + D0DHH / D0 - pow2(D0DH / D0) - D1DHH / D1 + pow2(D1DH / D1))
      E0 = 0.539409d0 + 2.522206d0 * THETA2 + 0.178484d0 * THETA4
      E0DH = 5.044412d0 * THETA + 0.713936d0 * THETA3
      E0DHH = 5.044412d0 + 2.141808d0 * THETA2
      E1 = 1d0 + 2.555501d0 * THETA2 + 0.146319d0 * THETA4
      E1DH = 5.111002d0 * THETA + 0.585276d0 * THETA3
      E1DHH = 5.111002d0 + 1.755828d0 * THETA2
      E = THETA * T1 * E0 / E1
      EH = 1.d0 / THETA + T1DH / T1 + E0DH / E0 - E1DH / E1
      EDH = E * EH
      EDHH = EDH * EH + E * (T1DHH / T1 - pow2(T1DH / T1) + E0DHH / E0 - pow2(E0DH / E0) &
         - E1DHH / E1 + pow2(E1DH / E1) - 1.0d0 / THETA2)
      EXP1TH = exp(-1.d0 / THETA)
      C = (0.872496d0 + 0.025248d0 * EXP1TH) * E
      CDH = 0.025248d0 * EXP1TH / THETA2 * E + C * EDH / E
      CDHH = 0.025248d0 * EXP1TH / THETA2 * (EDH + (1.0d0 - 2.0d0 * THETA) / THETA2 * E) &
         + CDH * EDH / E + C * EDHH / E - C * pow2(EDH / E)
      DISCR = SQRT(4.0d0 * E - D * D)
      DIDH = 0.5d0 / DISCR * (4.0d0 * EDH - 2.0d0 * D * DDH)
      DIDHH = (-pow2((2.0d0 * EDH - D * DDH) / DISCR) + 2.0d0 * EDHH - DDH * DDH - D * DDHH) / DISCR
      S1 = -C / E * GAME
      S1H = CDH / C - EDH / E
      S1DH = S1 * S1H
      S1DHH = S1DH * S1H + S1 * (CDHH / C - pow2(CDH / C) - EDHH / E + pow2(EDH / E))
      S1DG = -C / E ! => S1DGG=0
      S1DHG = S1DG * (CDH / C - EDH / E)
      B2 = B - C * D / E
      B2DH = BDH - (CDH * D + C * DDH) / E + C * D * EDH / (E * E)
      B2DHH = BDHH - (CDHH * D + 2d0 * CDH * DDH + C * DDHH) / E + (2d0 * (CDH * D + C * DDH - C * D * EDH / E) * EDH + C * D * EDHH) / (E * E)
      SQGE = SQRT(GAME)
      S2 = -2.d0 / E * B2 * SQGE
      S2H = B2DH / B2 - EDH / E
      S2DH = S2 * S2H
      S2DHH = S2DH * S2H + S2 * (B2DHH / B2 - pow2(B2DH / B2) - EDHH / E + pow2(EDH / E))
      S2DG = 0.5d0 * S2 / GAME
      S2DGG = -0.5d0 * S2DG / GAME
      S2DHG = 0.5d0 * S2DH / GAME
      R3 = E * GAME + D * SQGE + 1.0d0
      R3DH = EDH * GAME + DDH * SQGE
      R3DHH = EDHH * GAME + DDHH * SQGE
      R3DG = E + 0.5d0 * D / SQGE
      R3DGG = -0.25d0 * D / (GAME * SQGE)
      R3DHG = EDH + 0.5d0 * DDH / SQGE
      B3 = A - C / E
      B3DH = ADH - CDH / E + C * EDH / (E * E)
      B3DHH = ADHH - CDHH / E + (2.0d0 * CDH * EDH + C * EDHH) / (E * E) - 2d0 * C * EDH * EDH / pow3(E)
      C3 = (D / E * B2 - B3) / E
      C3DH = (DDH * B2 + D * B2DH + B3 * EDH) / (E * E) - 2d0 * D * B2 * EDH / pow3(E) - B3DH / E
      C3DHH = (-B3DHH &
         + (DDHH * B2 + 2d0 * DDH * B2DH + D * B2DHH + B3DH * EDH + B3 * EDHH + B3DH * EDH) / E &
         - 2.0d0 * ((DDH * B2 + D * B2DH + B3 * EDH + DDH * B2 + D * B2DH) * EDH + D * B2 * EDHH) / (E * E) &
         + 6.0d0 * D * B2 * EDH * EDH / pow3(E)) / E
      S3 = C3 * log(R3)
      S3DH = S3 * C3DH / C3 + C3 * R3DH / R3
      S3DHH = (S3DH * C3DH + S3 * C3DHH) / C3 - S3 * pow2(C3DH / C3) &
         + (C3DH * R3DH + C3 * R3DHH) / R3 - C3 * pow2(R3DH / R3)
      S3DG = C3 * R3DG / R3
      S3DGG = C3 * (R3DGG / R3 - pow2(R3DG / R3))
      S3DHG = (C3DH * R3DG + C3 * R3DHG) / R3 - C3 * R3DG * R3DH / (R3 * R3)
      B4 = 2.d0 - D * D / E
      B4DH = EDH * (D / E) * (D / E) - 2d0 * D * DDH / E
      B4DHH = EDHH * (D / E) * (D / E) + 2d0 * EDH * (D / E) * (D / E) * (DDH / D - EDH / E) &
         - 2d0 * (DDH * DDH + D * DDHH) / E + 2d0 * D * DDH * EDH / (E * E)
      C4 = 2d0 * E * SQGE + D
      C4DH = 2d0 * EDH * SQGE + DDH
      C4DHH = 2d0 * EDHH * SQGE + DDHH
      C4DG = E / SQGE
      C4DGG = -0.5d0 * E / (GAME * SQGE)
      C4DHG = EDH / SQGE
      S4A = 2.0d0 / E / DISCR
      S4AH = EDH / E + DIDH / DISCR
      S4ADH = -S4A * S4AH
      S4ADHH = -S4ADH * S4AH - S4A * (EDHH / E - (EDH / E) * (EDH / E) + DIDHH / DISCR - pow2(DIDH / DISCR))
      S4B = D * B3 + B4 * B2
      S4BDH = DDH * B3 + D * B3DH + B4DH * B2 + B4 * B2DH
      S4BDHH = DDHH * B3 + 2d0 * DDH * B3DH + D * B3DHH + B4DHH * B2 + 2d0 * B4DH * B2DH + B4 * B2DHH
      S4C = atan(C4 / DISCR) - atan(D / DISCR)
      UP1 = C4DH * DISCR - C4 * DIDH
      DN1 = DISCR * DISCR + C4 * C4
      UP2 = DDH * DISCR - D * DIDH
      DN2 = DISCR * DISCR + D * D
      S4CDH = UP1 / DN1 - UP2 / DN2
      S4CDHH = (C4DHH * DISCR - C4 * DIDHH) / DN1 &
         - UP1 * 2d0 * (DISCR * DIDH + C4 * C4DH) / (DN1 * DN1) &
         - (DDHH * DISCR - D * DIDHH) / DN2 + UP2 * 2d0 * (DISCR * DIDH + D * DDH) / (DN2 * DN2)
      S4CDG = C4DG * DISCR / DN1
      S4CDGG = C4DGG * DISCR / DN1 - 2d0 * C4 * DISCR * pow2(C4DG / DN1)
      S4CDHG = (C4DHG * DISCR + C4DG * DIDH - C4DG * DISCR / DN1 * 2d0 * (DISCR * DIDH + C4 * C4DH)) / DN1
      S4 = S4A * S4B * S4C
      S4DH = S4ADH * S4B * S4C + S4A * S4BDH * S4C + S4A * S4B * S4CDH
      S4DHH = S4ADHH * S4B * S4C + S4A * S4BDHH * S4C + S4A * S4B * S4CDHH &
         + 2d0 * (S4ADH * S4BDH * S4C + S4ADH * S4B * S4CDH + S4A * S4BDH * S4CDH)
      S4DG = S4A * S4B * S4CDG
      S4DGG = S4A * S4B * S4CDGG
      S4DHG = S4A * S4B * S4CDHG + S4CDG * (S4ADH * S4B + S4A * S4BDH)
      FXC = S1 + S2 + S3 + S4
      FXCDH = S1DH + S2DH + S3DH + S4DH
      FXCDG = S1DG + S2DG + S3DG + S4DG
      FXCDHH = S1DHH + S2DHH + S3DHH + S4DHH
      FXCDGG = S2DGG + S3DGG + S4DGG
      FXCDHG = S1DHG + S2DHG + S3DHG + S4DHG
      PXC = (GAME * FXCDG - 2d0 * THETA * FXCDH) / 3.d0
      UXC = GAME * FXCDG - THETA * FXCDH
      SXC = (GAME * S2DG - S2 + GAME * S3DG - S3 + S4A * S4B * (GAME * S4CDG - S4C)) - THETA * FXCDH
      if (abs(SXC).lt.1.d-9 * abs(THETA * FXCDH)) SXC = 0.d0 ! accuracy loss
      CVXC = 2d0 * THETA * (GAME * FXCDHG - FXCDH) - THETA * THETA * FXCDHH - GAME * GAME * FXCDGG
      if (abs(CVXC).lt.1.d-9 * abs(GAME * GAME * FXCDGG)) CVXC = 0.d0 ! accuracy
      PDLH = THETA * (GAME * FXCDHG - 2d0 * FXCDH - 2d0 * THETA * FXCDHH) / 3.d0
      PDLG = GAME * (FXCDG + GAME * FXCDGG - 2d0 * THETA * FXCDHG) / 3.d0
      PDRXC = PXC + (PDLG - 2d0 * PDLH) / 3.d0
      PDTXC = GAME * (THETA * FXCDHG - GAME * FXCDGG / 3.d0) - THETA * (FXCDH / 0.75d0 + THETA * FXCDHH / 1.5d0)
      return
   end subroutine EXCOR7


end module create_EXCOR7_table
