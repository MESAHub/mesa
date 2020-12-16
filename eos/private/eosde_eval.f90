! ***********************************************************************
!
!   Copyright (C) 2014-2019  Bill Paxton & The MESA Team
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

      module eosDE_eval
      use eos_def
      use const_def
      use math_lib

      implicit none


      contains


      subroutine Get_eosDE_Results(rq, &
               Z_in, X_in, abar, zbar, &
               species, chem_id, net_iso, xa, &
               aE, alogE, aRho, alogRho, logT_guess, &
               T, logT, res, d_dlnRho_c_T, d_dlnT_c_Rho,  &
               dlnT_dlnE_c_Rho, dlnT_dlnd_c_E, &
               dlnPgas_dlnE_c_Rho, dlnPgas_dlnd_c_E, &
               ierr)

         use utils_lib, only: is_bad
         
         ! INPUT
         
         type (EoS_General_Info), pointer :: rq ! general information about the request

         real(dp), intent(in) :: Z_in ! the desired Z
         real(dp), intent(in) :: X_in ! the desired X
            
         real(dp), intent(in) :: abar, zbar 
         
         integer, intent(in) :: species
         integer, pointer :: chem_id(:) ! maps species to chem id
            ! index from 1 to species
            ! value is between 1 and num_chem_isos         
         integer, pointer :: net_iso(:) ! maps chem id to species number
            ! index from 1 to num_chem_isos (defined in chem_def)
            ! value is 0 if the iso is not in the current net
            ! else is value between 1 and number of species in current net
         real(dp), intent(in) :: xa(:) ! mass fractions
         
         real(dp), intent(in) :: aE, alogE
            
         real(dp), intent(in) :: aRho, alogRho
            
         real(dp), intent(in) :: logT_guess

         
         ! OUTPUT    
              
         real(dp), intent(out) :: T, logT
         real(dp), intent(inout) :: res(:) ! (nv)
         real(dp), intent(inout) :: d_dlnRho_c_T(:) ! (nv)
         real(dp), intent(inout) :: d_dlnT_c_Rho(:) ! (nv)
         real(dp), intent(out) :: &
            dlnT_dlnE_c_Rho, dlnT_dlnd_c_E, &
            dlnPgas_dlnE_c_Rho, dlnPgas_dlnd_c_E
         integer, intent(out) :: ierr
         
         ! LOCALS
         
         real(dp), dimension(nv) ::  &
            res1, d_dlnRho_c_T1, d_dlnT_c_Rho1
         real(dp), dimension(nv) ::  &
            res2, d_dlnRho_c_T2, d_dlnT_c_Rho2
         real(dp) :: X, Z, &
               dlnT_dlnE_c_Rho1, dlnT_dlnd_c_E1, &
               dlnPgas_dlnE_c_Rho1, dlnPgas_dlnd_c_E1, &
               dlnT_dlnE_c_Rho2, dlnT_dlnd_c_E2, &
               dlnPgas_dlnE_c_Rho2, dlnPgas_dlnd_c_E2
         integer :: iregion
         character (len=256) :: message
   
         integer, parameter :: pure_helm = 1
         integer, parameter :: pure_opal_scvh = 2
         integer, parameter :: blend_in_x = 3
         integer, parameter :: blend_in_y = 4
         integer, parameter :: blend_corner_out = 5
         
         real(dp) :: Rho, logRho, energy, logE, logV, alfa, beta, &
               T1, logT1, T2, logT2, tiny
         
         logical, parameter :: dbg = .false.

         include 'formats.dek'
         
         ierr = 0
         tiny = rq% tiny_fuzz

         if (is_bad(X_in) .or. is_bad(Z_in)) then
            ierr = -1
            if (dbg) write(*,*) 'error from is_bad'
            return
         end if
         
         X = X_in; Z = Z_in
         if (X < tiny) X = 0
         if (Z < tiny) Z = 0
         
         call get_DE_args( &
            aE, alogE, aRho, alogRho, energy, logE, Rho, logRho, ierr)
         if (ierr /= 0) then
            if (dbg) write(*,*) 'error from get_DE_args'
            return
         end if

         call Get_DE_Results_using_DT( &
            rq, Z, X, abar, zbar, &
            species, chem_id, net_iso, xa, &
            energy, logE, Rho, logRho, logT_guess, &
            T, logT, res, d_dlnRho_c_T, d_dlnT_c_Rho,  &
            dlnT_dlnE_c_Rho, dlnT_dlnd_c_E, &
            dlnPgas_dlnE_c_Rho, dlnPgas_dlnd_c_E, &
            ierr)

         ! zero blend fractions; not supported for eosDE
         res(i_frac:i_frac+num_eos_frac_results-1) = 0.0
         d_dlnRho_c_T(i_frac:i_frac+num_eos_frac_results-1) = 0.0
         d_dlnT_c_Rho(i_frac:i_frac+num_eos_frac_results-1) = 0.0

         ! zero phase information
         res(i_phase:i_latent_ddlnRho) = 0d0
         d_dlnT_c_Rho(i_phase:i_latent_ddlnRho) = 0d0
         d_dlnRho_c_T(i_phase:i_latent_ddlnRho) = 0d0

         
      end subroutine Get_eosDE_Results
      
      
      subroutine get_DE_args( &
            aE, alogE, aRho, alogRho, energy, logE, Rho, logRho, ierr)       
         real(dp), intent(in) :: aE, alogE, aRho, alogRho
         real(dp), intent(out) :: energy, logE, Rho, logRho
         integer, intent(out) :: ierr
         include 'formats.dek'
         ierr = 0
         energy = aE; logE = alogE
         if (aE == arg_not_provided .and. alogE == arg_not_provided) then
            ierr = -2; return
         end if
         if (alogE == arg_not_provided) logE = log10(energy)
         if (aE == arg_not_provided) energy = exp10(logE)
         if (energy <= 0) then
            ierr = -1
            return
         end if
         Rho = aRho; logRho = alogRho
         if (Rho == arg_not_provided .and. logRho == arg_not_provided) then
            ierr = -3; return
         end if
         if (logRho == arg_not_provided) logRho = log10(Rho)
         if (Rho == arg_not_provided) Rho = exp10(logRho)
         if (Rho <= 0) then
            ierr = -1
            return
         end if
      end subroutine get_DE_args


      subroutine Get_DE_Results_using_DT( &
               rq, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               energy, logE, Rho, logRho, logT_guess, &
               T, logT, res, d_dlnRho_c_T, d_dlnT_c_Rho,  &
               dlnT_dlnE_c_Rho, dlnT_dlnd_c_E, &
               dlnPgas_dlnE_c_Rho, dlnPgas_dlnd_c_E, &
               ierr)
         use eosDT_eval, only: get_T
         use utils_lib, only: is_bad
         
         type (EoS_General_Info), pointer :: rq ! general information about the request
         real(dp), intent(in) :: Z, X, abar, zbar
         integer, intent(in) :: species
         integer, pointer :: chem_id(:)    
         integer, pointer :: net_iso(:)
         real(dp), intent(in) :: xa(:)
         real(dp), intent(in) :: energy, logE, Rho, logRho, logT_guess
         real(dp), intent(out) :: T, logT
         real(dp), intent(inout) :: res(:) ! (nv)
         real(dp), intent(inout) :: d_dlnRho_c_T(:) ! (nv)
         real(dp), intent(inout) :: d_dlnT_c_Rho(:) ! (nv)
         real(dp), intent(out) :: &
            dlnT_dlnE_c_Rho, dlnT_dlnd_c_E, &
            dlnPgas_dlnE_c_Rho, dlnPgas_dlnd_c_E
         integer, intent(out) :: ierr
         
         logical, parameter :: basic_flag = .false.
         integer:: i, eos_calls, max_iter, which_other
         real(dp) ::  &
            other, other_tol, logT_tol, &
            logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, logT_result
         
         logical, parameter :: dbg = .false.
         
         include 'formats.dek'
         
         ierr = 0

         which_other = i_lnE
         other = logE*ln10
         other_tol = 1d-8
         logT_tol = 1d-8
               
         logT_bnd1 = arg_not_provided
         logT_bnd2 = arg_not_provided
         other_at_bnd1 = arg_not_provided
         other_at_bnd2 = arg_not_provided

         max_iter = 20
         eos_calls = 0
         
         if (dbg) write(*,1) 'logRho', logRho
         if (dbg) write(*,1) 'logE', logE
         if (dbg) write(*,1) 'logT_guess', logT_guess
         if (dbg) write(*,1) 'T_guess', exp10(logT_guess)
         if (dbg) write(*,1) 'T_guess_gas', 2*energy*abar*mp/(3*kerg*(1+zbar)) ! ideal gas
         if (dbg) write(*,1) 'T_guess_rad', pow(energy/crad,0.25d0)
         
         call get_T( &
               rq% handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               logRho, which_other, other, &
               logT_tol, other_tol, max_iter, logT_guess,  &
               logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
               logT_result, res, d_dlnRho_c_T, d_dlnT_c_Rho,  &
               eos_calls, ierr)
         if (ierr /= 0 .or. (dbg .and. abs(logT_result) > 20)) then
            if (dbg) then
               write(*,*) 'failed in get_T for Get_DE_Results_using_DT'
               write(*,1) 'Z = ', Z
               write(*,1) 'X = ', X
               write(*,1) 'abar = ', abar
               write(*,1) 'zbar = ', zbar
               write(*,1) 'logRho = ', logRho
               write(*,1) 'logE = ', logE
               write(*,1) 'logT_guess = ', logT_guess
               write(*,*)
               stop 'Get_DE_Results_using_DT'
            end if
            return
         end if
         
         logT = logT_result
         T = exp10(logT)

         dlnT_dlnd_c_E = -rho*res(i_dE_dRho)/(T*res(i_Cv))
         dlnPgas_dlnd_c_E = &
            d_dlnRho_c_T(i_lnPgas) + d_dlnT_c_Rho(i_lnPgas)*dlnT_dlnd_c_E

         dlnT_dlnE_c_Rho = energy/(T*res(i_Cv))
         dlnPgas_dlnE_c_Rho = d_dlnT_c_Rho(i_lnPgas)*dlnT_dlnE_c_Rho

         if (dbg) write(*,1) 'logT', logT
         if (dbg) write(*,1) 'T', T

      end subroutine Get_DE_Results_using_DT
      
      
      subroutine Get_eos_gamma_DE_Results( &
            rq, abar, energy, log10E, rho, log10Rho, gamma, &
            T, log10T, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
            dlnT_dlnE_c_Rho, dlnT_dlnd_c_E, &
            dlnPgas_dlnE_c_Rho, dlnPgas_dlnd_c_E, ierr)
         use utils_lib, only: is_bad
         type (EoS_General_Info), pointer :: rq
         real(dp), intent(in) :: abar, energy, log10E, Rho, log10Rho, gamma
         real(dp), intent(out) :: T, log10T
         real(dp), intent(inout), dimension(:) :: &
            res, d_dlnRho_const_T, d_dlnT_const_Rho
         real(dp), intent(out) :: & 
            dlnT_dlnE_c_Rho, dlnT_dlnd_c_E, &
            dlnPgas_dlnE_c_Rho, dlnPgas_dlnd_c_E
         integer, intent(out) :: ierr
         
         real(dp) :: avo_k_div_abar, P, entropy
         include 'formats'
         
         ierr = 0
         
         res(1:nv) = 0d0
         d_dlnRho_const_T(1:nv) = 0d0
         d_dlnT_const_Rho(1:nv) = 0d0

         avo_k_div_abar = avo*kerg/abar
         
         P = (gamma - 1d0)*energy*rho
         T = (gamma - 1d0)*energy/avo_k_div_abar
         log10T = log10(T)

         res(i_Cv) = avo_k_div_abar/(gamma - 1) ! energy/T
         
         entropy = res(i_Cv)*log(P/pow(rho,gamma)) + 1d9 ! offset to keep it > 0
         
         if (is_bad(entropy) .or. entropy <= 0d0) then
            if (.false.) then
!$OMP critical (eosde_eval_crit1)
               write(*,1) 'entropy', entropy
               write(*,1) 'energy', energy
               write(*,1) 'T', T
               write(*,1) 'P', P
               write(*,1) 'rho', rho
               write(*,1) 'pow(rho,gamma)', pow(rho,gamma)
               write(*,1) 'res(i_Cv)', res(i_Cv)
               write(*,1) 'P/pow(rho,gamma)', P/pow(rho,gamma)
               write(*,1) 'log(P/pow(rho,gamma)', log(P/pow(rho,gamma))
               write(*,1) 'gamma', gamma
               stop 'bad gamma law entropy'
!$OMP end critical (eosde_eval_crit1)
            end if
            entropy = 1d-99
         end if
         
         res(i_lnPgas) = log(P) ! treat P as Pgas
         res(i_lnE) = log10E*ln10
         res(i_lnS) = log(entropy)
         res(i_Cp) = gamma*res(i_Cv)
         res(i_grad_ad) = (gamma - 1d0)/gamma ! dlnT_dlnP|S
         res(i_chiRho) = 1d0 ! dlnP_dlnRho|T
         res(i_chiT) = 1d0 ! dlnP_dlnT|Rho
         res(i_dE_dRho) = 0d0
         res(i_dS_dT) = 0d0
         res(i_dS_dRho) = 0d0
         res(i_mu) = 1d0/abar
         res(i_gamma1) = gamma
         res(i_gamma3) = gamma
         res(i_lnfree_e) = -1d99
         res(i_eta) = 0d0
         
         d_dlnRho_const_T(i_lnPgas) = 1
         
         d_dlnT_const_Rho(i_lnPgas) = 1
         d_dlnT_const_Rho(i_lnE) = 1
         
         dlnT_dlnd_c_E = 0
         dlnPgas_dlnd_c_E = 1

         dlnT_dlnE_c_Rho = 1
         dlnPgas_dlnE_c_Rho = 1
         
      end subroutine Get_eos_gamma_DE_Results
      
      
      subroutine Do_eosDE_Interpolations( &
                   nv, nx, x, ny, y, fin1, i, j, &
                   x0, xget, x1, y0, yget, y1, &
                   fval, df_dx, df_dy, ierr)
         integer, intent(in) :: nv, nx, ny
         real(dp), intent(in) :: x(:) ! (nx)
         real(dp), intent(in) :: y(:) ! (ny)
         real(dp), intent(in), pointer :: fin1(:) ! =(4,nv,nx,ny)
         integer, intent(in) :: i, j           ! target cell in f
         real(dp), intent(in) :: x0, xget, x1      ! x0 <= xget <= x1;  x0 = xs(i), x1 = xs(i+1)
         real(dp), intent(in) :: y0, yget, y1      ! y0 <= yget <= y1;  y0 = ys(j), y1 = ys(j+1)
         real(dp), intent(inout), dimension(nv) :: &
            fval, df_dx, df_dy
         integer, intent(out) :: ierr
   
         real(dp), parameter :: z36th = 1d0/36d0
         real(dp) :: xp, xpi, xp2, xpi2, ax, axbar, bx, bxbar, cx, cxi, hx2, cxd, cxdi, hx, hxi
         real(dp) :: yp, ypi, yp2, ypi2, ay, aybar, by, bybar, cy, cyi, hy2, cyd, cydi, hy, hyi
         real(dp) :: sixth_hx2, sixth_hy2, z36th_hx2_hy2
         real(dp) :: sixth_hx, sixth_hxi_hy2, z36th_hx_hy2
         real(dp) :: sixth_hx2_hyi, sixth_hy, z36th_hx2_hy
         integer :: k, ip1, jp1
         real(dp), pointer :: fin(:,:,:,:)
         
         include 'formats.dek'
         
         ierr = 0
         
         fin(1:4,1:nv,1:nx,1:ny) => fin1(1:4*nv*nx*ny)
         
         hx=x1-x0
         hxi=1.0d0/hx
         hx2=hx*hx
   
         xp=(xget-x0)*hxi
         xpi=1.0d0-xp
         xp2=xp*xp
         xpi2=xpi*xpi

         ax=xp2*(3.d0-2.d0*xp)
         axbar=1.0d0-ax
         
         bx=-xp2*xpi
         bxbar=xpi2*xp
   
         cx=xp*(xp2-1.0d0)
         cxi=xpi*(xpi2-1.0d0)
         cxd=3.0d0*xp2-1.0d0
         cxdi=-3.0d0*xpi2+1.0d0
   
         hy=y1-y0
         hyi=1.0d0/hy
         hy2=hy*hy
   
         yp=(yget-y0)*hyi
         ypi=1.0d0-yp
         yp2=yp*yp
         ypi2=ypi*ypi

         ay=yp2*(3.0d0-2.0d0*yp)
         aybar=1.d0-ay
         
         by=-yp2*ypi
         bybar=ypi2*yp
   
         cy=yp*(yp2-1.0d0)
         cyi=ypi*(ypi2-1.0d0)
         cyd=3.0d0*yp2-1.0d0
         cydi=-3.0d0*ypi2+1.0d0
                  
         sixth_hx2 = one_sixth*hx2
         sixth_hy2 = one_sixth*hy2
         z36th_hx2_hy2 = z36th*hx2*hy2
         
         sixth_hx = one_sixth*hx
         sixth_hxi_hy2 = one_sixth*hxi*hy2
         z36th_hx_hy2 = z36th*hx*hy2
         
         sixth_hx2_hyi = one_sixth*hx2*hyi
         sixth_hy = one_sixth*hy
         z36th_hx2_hy = z36th*hx2*hy
         
         ip1 = i+1
         jp1 = j+1
         
         do k=1,nv
            ! bicubic spline interpolation
            fval(k) = &
                  xpi*( &
                     ypi*fin(1,k,i,j)  +yp*fin(1,k,i,jp1)) &
                     +xp*(ypi*fin(1,k,ip1,j)+yp*fin(1,k,ip1,jp1)) &
                  +sixth_hx2*( &
                     cxi*(ypi*fin(2,k,i,j) +yp*fin(2,k,i,jp1))+ &
                     cx*(ypi*fin(2,k,ip1,j)+yp*fin(2,k,ip1,jp1))) &
                  +sixth_hy2*( &
                     xpi*(cyi*fin(3,k,i,j) +cy*fin(3,k,i,jp1))+ &
                     xp*(cyi*fin(3,k,ip1,j)+cy*fin(3,k,ip1,jp1))) &
                  +z36th_hx2_hy2*( &
                     cxi*(cyi*fin(4,k,i,j) +cy*fin(4,k,i,jp1))+ &
                     cx*(cyi*fin(4,k,ip1,j)+cy*fin(4,k,ip1,jp1)))
            
            ! derivatives of bicubic splines
            df_dx(k) = &
                  hxi*( &
                     -(ypi*fin(1,k,i,j)  +yp*fin(1,k,i,jp1)) &
                     +(ypi*fin(1,k,ip1,j)+yp*fin(1,k,ip1,jp1))) &
                  +sixth_hx*( &
                     cxdi*(ypi*fin(2,k,i,j) +yp*fin(2,k,i,jp1))+ &
                     cxd*(ypi*fin(2,k,ip1,j)+yp*fin(2,k,ip1,jp1))) &
                  +sixth_hxi_hy2*( &
                     -(cyi*fin(3,k,i,j)  +cy*fin(3,k,i,jp1)) &
                     +(cyi*fin(3,k,ip1,j)+cy*fin(3,k,ip1,jp1))) &
                  +z36th_hx_hy2*( &
                     cxdi*(cyi*fin(4,k,i,j) +cy*fin(4,k,i,jp1))+ &
                     cxd*(cyi*fin(4,k,ip1,j)+cy*fin(4,k,ip1,jp1)))

            df_dy(k) = &
                  hyi*( &
                     xpi*(-fin(1,k,i,j) +fin(1,k,i,jp1))+ &
                     xp*(-fin(1,k,ip1,j)+fin(1,k,ip1,jp1))) &
                  +sixth_hx2_hyi*( &
                     cxi*(-fin(2,k,i,j) +fin(2,k,i,jp1))+ &
                     cx*(-fin(2,k,ip1,j)+fin(2,k,ip1,jp1))) &
                  +sixth_hy*( &
                     xpi*(cydi*fin(3,k,i,j) +cyd*fin(3,k,i,jp1))+ &
                     xp*(cydi*fin(3,k,ip1,j)+cyd*fin(3,k,ip1,jp1))) &
                  +z36th_hx2_hy*( &
                     cxi*(cydi*fin(4,k,i,j) +cyd*fin(4,k,i,jp1))+ &
                     cx*(cydi*fin(4,k,ip1,j)+cyd*fin(4,k,ip1,jp1)))
         end do
         
      end subroutine Do_eosDE_Interpolations


      end module eosDE_eval
      
