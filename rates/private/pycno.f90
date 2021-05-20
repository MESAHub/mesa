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





      module pycno
      use rates_def
      use utils_lib
      use const_def, only: dp
      use math_lib
      
      implicit none
      
      
      contains
      
      subroutine FL_epsnuc_3alf(T, Rho, Y, UE, r, drdT, drdRho)       
         real(dp), intent(in) :: T ! temperature
         real(dp), intent(in) :: Rho ! density
         real(dp), intent(in) :: Y ! helium mass fraction
         real(dp), intent(in) :: UE ! electron molecular weight
         real(dp), intent(out) :: r ! rate in ergs/g/sec
         real(dp), intent(out) :: drdT ! partial wrt temperature
         real(dp), intent(out) :: drdRho ! partial wrt density


         real(dp) :: T6, R6, R6T, R6T13, R6T16, T62, T612, T623, T653, T632, T613, U, AF
         real(dp) :: G1, dG1dRho, dG1dT, G2, dG2dRho, dG2dT
         real(dp) :: B1, dB1dRho, dB1dT, B2, dB2dRho, dB2dT
         real(dp) :: dUdT, dUdRho, U32, U52, dAFdT, dAFdRho
         real(dp) :: E1, dE1dT, dE1dRho, E2, dE2dT, dE2dRho
         real(dp) :: F1, dF1dT, dF1dRho, F2, dF2dT, dF2dRho
         real(dp) :: dR6dRho, dR6TdRho, dR6T13dRho, dR6T16dRho
         real(dp) :: dT6dT, dT612dT, dT62dT, dT613dT, dT623dT, dT632dT, dT653dT
         
         ! DEBUG
         real(dp), parameter :: AF_0  =  1.9005324047511074D+00
         real(dp), parameter :: B1_denom_0  =  2.9602238143383192D-01
         real(dp), parameter :: B1_0  =  1.2227955158250397D-08
         real(dp), parameter :: B2_denom_0  =  1.7563773044362474D+00
         real(dp), parameter :: B2_0  =  1.0173166567483392D-14
         real(dp), parameter :: E1_0  =  -2.2308014220480969D+00
         real(dp), parameter :: F1_0  =  1.5176626709750911D-04
         real(dp), parameter :: E2_0  =  -2.2350904778008243D+01
         real(dp), parameter :: F2_0  =  2.7741209323605414D-13
         real(dp), parameter :: T_0  =  7.9432823472428218D+07
         real(dp), parameter :: RHO_0  =  3.1622776917911558D+09
         real(dp), parameter :: r_0  =  2.2348420508311778D+20
         real(dp), parameter :: G1_0  =  1.5177849505266735D-04
         real(dp), parameter :: G2_0  =  2.8758525980353755D-13
         real(dp), parameter :: U_0  =  1.0723431204522564D+00

         real(dp) :: tmp, &
               B1_denom, dB1_denom_dRho, dB1_denom_dT, &
               B2_denom, dB2_denom_dRho, dB2_denom_dT
         real(dp) :: A1, dA1dT, B1_numerator, dB1_numerator_dT
         real(dp) :: A2, dA2dT, B2_numerator, dB2_numerator_dT

         include 'formats'
         
         R6=RHO*1d-6
         dR6dRho = 1d-6
         R6T=2d0*R6/UE
         dR6TdRho = 2d0*dR6dRho/UE

         R6T16=pow(R6T,1d0/6d0)
         dR6T16dRho = (1d0/6d0)*dR6TdRho*R6T16/R6T
         R6T13=R6T16*R6T16       
         dR6T13dRho = 2*R6T16*dR6T16dRho
         T6=T*1d-6
         dT6dT=1d-6
         dT62dT=2d0*T6*dT6dT
         T613=pow(T6,1d0/3d0)
         dT613dT=(1d0/3d0)*dT6dT*T613/T6
         T623=T613*T613
         dT623dT=2*T613*dT613dT
         T653=T623*T6
         dT653dT = dT623dT*T6 + T623*dT6dT
         
         T62=T6*T6
         T612=sqrt(T6)
         dT612dT=0.5d0*dT6dT/T612
         T632=T6*T612
         dT632dT=1.5d0*T612*dT6dT        
         
         U=1.35D0*R6T13/T623
         dUdT = -U * dT623dT / T623
         dUdRho = U * dR6T13dRho / R6T13
         
         U32 = U*sqrt(U)
         U52 = U*U32

         if (U < 1) then ! strong screening regime, eqn 4.8a in F&L
         
            A1 = pow2(1d0-4.222D-2*T623) + 2.643D-5*T653
            dA1dT = -2d0*4.222d-2*dT623dT*(1d0 - 4.222D-2*T623) + 2.643D-5*dT653dT
            
            B1_denom=A1*T623
            dB1_denom_dT = dA1dT*T623 + A1*dT623dT
            
            B1_numerator = 16.16D0*exp(-134.92d0/T613)
            dB1_numerator_dT = B1_numerator*134.92d0*dT613dT/(T613*T613)
            
            B1=B1_numerator/B1_denom
            dB1dT = dB1_numerator_dT/B1_denom - B1*dB1_denom_dT/B1_denom
            
            A2 = pow2(1d0-2.807D-2*T623) + 2.704D-6*T653
            dA2dT = -2*2.807D-2*dT623dT*(1d0-2.807D-2*T623) + 2.704D-6*dT653dT
            
            B2_denom=A2*T623
            dB2_denom_dT = dA2dT*T623 + A2*dT623dT
            
            B2_numerator = 244.6D0*pow5(1D0+3.528D-3*T623) * exp(-235.72D0/T613)
            dB2_numerator_dT = B2_numerator* &
                  (5D0*3.528D-3*dT623dT/(1D0+3.528D-3*T623) + 235.72D0*dT613dT/T623)

            B2=B2_numerator/B2_denom
            dB2dT = dB2_numerator_dT/B2_denom - B2*dB2_denom_dT/B2_denom
            
            if (5.458D3 > R6T) then
            
               E1 = -1065.1D0/T6
               dE1dT = -E1 * dT6dT / T6
               
               F1 = exp(E1)/T632
               dF1dT = F1 * (dE1dT - dT632dT/T632)
               
               B1=B1+F1
               dB1dT = dB1dT + dF1dT
               
            endif
            
            if (1.836D4 > R6T) then
            
               E2 = -3336.4D0/T6
               dE2dT = -E2 * dT6dT / T6
               
               F2 = exp(E2)/T632
               dF2dT = F2 * (dE2dT - dT632dT/T632)
            
               B2=B2+F2
               dB2dT = dB2dT + dF2dT
               
            endif
            
            G1=B1*exp(60.492D0*R6T13/T6)
            dG1dT = G1*(dB1dT/B1 - 60.492D0*R6T13*dT6dT/(T6*T6))
            dG1dRho=0
            
            G2=B2*exp(106.35D0*R6T13/T6)
            dG2dT = G2*(dB2dT/B2 - 106.35D0*R6T13*dT6dT/(T6*T6))
            dG2dRho=0            

         else ! pycnonuclear regime, eqn 4.8b in F&L
         
            AF=1d0/U32 + 1d0
            dAFdT = -1.5d0 * dUdT/U52
            dAFdRho = -1.5d0 * dUdRho/U52
         
            B1_denom=T612*(pow2(1.0d0-5.680D-2*R6T13)+8.815D-7*T62)
            dB1_denom_dT = B1_denom*dT612dT/T612 + T612*8.815D-7*dT62dT
            dB1_denom_dRho = -2*5.680D-2*(1.0d0-5.680D-2*R6T13)*T612*dR6T13dRho
            
            B1=1.178D0*AF*exp(-77.554d0/R6T16)/B1_denom
            dB1dT = B1 * (dAFdT/AF - dB1_denom_dT/B1_denom)
            dB1dRho = B1 * (dAFdRho/AF + 77.554d0*dR6T16dRho/(R6T16*R6T16) - dB1_denom_dRho/B1_denom)
         
            B2_denom=T612*(pow2(1.0d0-3.791D-2*R6T13)+5.162D-8*T62)
            dB2_denom_dT = B2_denom*dT612dT/T612 + T612*5.162D-8*dT62dT
            
            tmp = pow(Rho/UE,1d0/3d0)
            dB2_denom_dRho = T612*(-0.000252733d0 + 9.58112d-8*tmp)*tmp/Rho
            
            B2=13.48D0*AF*pow5(1.0d0+5.070D-3*R6T13)*exp(-135.08D0/R6T16)/B2_denom
            dB2dT = B2 * (dAFdT/AF - dB2_denom_dT/B2_denom)
            dB2dRho = B2 * (dAFdRho/AF + 135.08D0*dR6T16dRho/(R6T16*R6T16) - dB2_denom_dRho/B2_denom)            
            
            if (5.458D3 > R6T) then
            
               E1 = (60.492d0*R6T13 - 1065.1D0)/T6
               dE1dT = -E1 * dT6dT / T6
               dE1dRho = 60.492d0*dR6T13dRho/T6
               
               F1 = exp(E1)/T632
               dF1dT = F1 * (dE1dT - dT632dT/T632)
               dF1dRho = F1 * dE1dRho
               
               !write(*,1) 'E1', E1
               !write(*,1) 'F1', F1

               G1=B1+F1
               dG1dT = dB1dT + dF1dT
               dG1dRho = dB1dRho + dF1dRho
               
            else
            
               G1=B1; dG1dRho = dB1dRho; dG1dT = dB1dT
               
            endif
            
            if (1.836D4 > R6T) then
            
               E2 = (106.35D0*R6T13 - 3336.4D0)/T6
               dE2dT = -E2 * dT6dT / T6
               dE2dRho = 106.35D0*dR6T13dRho/T6
               
               F2 = exp(E2)/T632
               dF2dT = F2 * (dE2dT - dT632dT/T632)
               dF2dRho = F2 * dE2dRho
               
               !write(*,1) 'E2', E2
               !write(*,1) 'F2', F2
            
               G2=B2+F2
               dG2dT = dB2dT + dF2dT
               dG2dRho = dB2dRho + dF2dRho
               
            else
            
               G2=B2; dG2dRho = dB2dRho; dG2dT = dB2dT
               
            endif

         endif
      
         r=5.120D29*G1*G2*Y*Y*Y*R6*R6 ! ergs/g/sec, eqn 4.7 in F&L

         if (r < 1d-99 .or. G1 < 1d-99 .or. G2 < 1d-99) then
            drdT = 0
            drdRho = 0
         else
            drdT = r * (dG1dT/G1 + dG2dT/G2)
            drdRho = r * (dG1dRho/G1 + dG2dRho/G2 + 2*dR6dRho/R6)

            return
         
            write(*,1) 'T', T
            write(*,1) 'RHO', RHO
            write(*,1) 'r', r
            write(*,1) 'G1', G1
            write(*,1) 'G2', G2
            write(*,1) 'U', U
            write(*,*)
            
            write(*,1) 'abs(Rho_0 - Rho)', abs(Rho_0 - Rho)
            
            if (.true. .and. abs(Rho_0 - Rho) > 1d-2) then
               write(*,*)
               write(*,1) 'analytic drdRho', drdRho
               write(*,1) 'numeric drdRho', (r_0 - r) / (Rho_0 - Rho)
               write(*,*)
               write(*,1) 'analytic dG1dRho', dG1dRho
               write(*,1) 'numeric dG1dRho', (G1_0 - G1) / (Rho_0 - Rho)
               write(*,*)
               write(*,1) 'analytic dG2dRho', dG2dRho
               write(*,1) 'numeric dG2dRho', (G2_0 - G2) / (Rho_0 - Rho)
               write(*,*)
               write(*,1) 'analytic dUdRho', dUdRho
               write(*,1) 'numeric dUdRho', (U_0 - U) / (Rho_0 - Rho)
               write(*,*)
               write(*,1) 'analytic AF', dAFdRho 
               write(*,1) 'numeric AF', (AF_0 - AF) / (Rho_0 - Rho)
               write(*,*)
               write(*,1) 'analytic B1_denom', dB1_denom_dRho
               write(*,1) 'numeric B1_denom', (B1_denom_0 - B1_denom) / (Rho_0 - Rho)
               write(*,*)
               write(*,1) 'analytic B1', dB1dRho
               write(*,1) 'numeric B1', (B1_0 - B1) / (Rho_0 - Rho)
               write(*,*)
               write(*,1) 'analytic B2_denom', dB2_denom_dRho 
               write(*,1) 'numeric B2_denom', (B2_denom_0 - B2_denom) / (Rho_0 - Rho)
               write(*,*)
               write(*,1) 'analytic B2', dB2dRho 
               write(*,1) 'numeric B2', (B2_0 - B2) / (Rho_0 - Rho)
               write(*,*)
               write(*,1) 'analytic E1', dE1dRho 
               write(*,1) 'numeric E1', (E1_0 - E1) / (Rho_0 - Rho)
               write(*,*)
               write(*,1) 'analytic F1', dF1dRho 
               write(*,1) 'numeric F1', (F1_0 - F1) / (Rho_0 - Rho)
               write(*,*)
               write(*,1) 'analytic E2', dE2dRho 
               write(*,1) 'numeric E2', (E2_0 - E2) / (Rho_0 - Rho)
               write(*,*)
               write(*,1) 'analytic F2', dF2dRho 
               write(*,1) 'numeric F2', (F2_0 - F2) / (Rho_0 - Rho)
               write(*,*)
               stop 'FL_epsnuc_3alf' 
            end if
                        
         end if

      end subroutine FL_epsnuc_3alf


      end module pycno
      
      


