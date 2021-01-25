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
!
! ***********************************************************************

      module eosPT_eval
      use eos_def
      use const_def
      use math_lib
      use utils_lib, only: mesa_error

      implicit none

         
      integer, parameter :: i_doing_which = 1
      integer, parameter :: i_which_other = 2
      integer, parameter :: i_handle = 3
      integer, parameter :: i_count = 4
      integer, parameter :: i_species = 5
      
      integer, parameter :: eos_lipar = 5

      integer, parameter :: r_other_value = 1
      integer, parameter :: r_Z = 2
      integer, parameter :: r_X = 3
      integer, parameter :: r_abar = 4
      integer, parameter :: r_zbar = 5
      integer, parameter :: r_Pgas = 6
      integer, parameter :: r_T = 7
      integer, parameter :: r_the_other_log = 8

      integer, parameter :: eos_lrpar = 8
      
      integer, parameter :: doing_get_T = 1
      integer, parameter :: doing_get_Pgas = 2
      integer, parameter :: doing_get_Pgas_for_Rho = 3

      

      
      contains


      subroutine Get_eosPT_Results(rq, &
               Z_in, X_in, abar, zbar, &
               species, chem_id, net_iso, xa, &
               aPgas, alogPgas, atemp, alogtemp, &
               Rho, logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas, &
               res, d_dlnRho_c_T, d_dlnT_c_Rho, &
               ierr)
         use utils_lib, only: is_bad
         type (EoS_General_Info), pointer :: rq
         real(dp), intent(in) :: Z_in, X_in, abar, zbar
         integer, intent(in) :: species
         integer, pointer :: chem_id(:), net_iso(:)
         real(dp), intent(in) :: xa(:)
         real(dp), intent(in) :: aPgas, alogPgas, atemp, alogtemp
         real(dp), intent(out) :: Rho, logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas
         real(dp), intent(inout) :: res(:), d_dlnRho_c_T(:), d_dlnT_c_Rho(:) ! (nv)
         integer, intent(out) :: ierr
         
         real(dp) :: X, Z, T, logT
         real(dp) :: Pgas, logPgas, Prad, tiny
         logical, parameter :: dbg = .false.
         
         logical :: skip

         include 'formats'
         
         ierr = 0
         tiny = rq% tiny_fuzz
         
         if (is_bad(X_in) .or. is_bad(Z_in)) then
            ierr = -1
            return
         end if
         
         X = X_in; Z = Z_in
         if (X < tiny) X = 0d0
         if (Z < tiny) Z = 0d0
         
         if (X > 1d0) then
            if (X > 1.0001D0) then
               write(*,1) 'Get_eosPT_Results: X bad', X
               ierr = -1
               return
               stop 'eosPT'
            end if
            X = 1d0
         end if
         
         call get_PT_args( &
            aPgas, alogPgas, atemp, alogtemp, Pgas, logPgas, T, logT, ierr)
         if (ierr /= 0) then
            if (dbg) write(*,*) 'error from get_PT_args'
            return
         end if
         
         if (Pgas <= 0) then
            ierr = -1
            return
         end if
         
         if (is_bad(Pgas) .or. is_bad(T)) then
            ierr = -1
            return
         end if
         
         call Get1_eosPT_Results(rq, &
            Z, X, abar, zbar, &
            species, chem_id, net_iso, xa, &
            Pgas, logPgas, T, logT, &
            Rho, logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas, &
            res, d_dlnRho_c_T, d_dlnT_c_Rho, ierr)

         ! zero blend fractions; not supported for eosPT
         res(i_frac:i_frac+num_eos_frac_results-1) = 0.0
         d_dlnRho_c_T(i_frac:i_frac+num_eos_frac_results-1) = 0.0
         d_dlnT_c_Rho(i_frac:i_frac+num_eos_frac_results-1) = 0.0

         ! zero phase information
         res(i_phase:i_latent_ddlnRho) = 0d0
         d_dlnT_c_Rho(i_phase:i_latent_ddlnRho) = 0d0
         d_dlnRho_c_T(i_phase:i_latent_ddlnRho) = 0d0

                   
      end subroutine Get_eosPT_Results


      subroutine Get1_eosPT_Results(rq, &
               Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               Pgas, logPgas, T, logT, &
               Rho, logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas, &
               res, d_dlnRho_c_T, d_dlnT_c_Rho, ierr)

         use utils_lib, only: is_bad
         use eosPT_load_tables, only: load_single_eosPT_table
                  
         type (EoS_General_Info), pointer :: rq
         real(dp), intent(in) :: Z, X, abar, zbar
         integer, intent(in) :: species
         integer, pointer :: chem_id(:), net_iso(:)
         real(dp), intent(in) :: xa(:)
         real(dp), intent(in) :: Pgas, logPgas, T, logT
         real(dp), intent(out) :: Rho, logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas
         real(dp), intent(inout) :: res(:), d_dlnRho_c_T(:), d_dlnT_c_Rho(:) ! (nv)
         integer, intent(out) :: ierr
         
         real(dp), dimension(nv) :: &
            res1, d_dlnRho_c_T1, d_dlnT_c_Rho1
         real(dp), dimension(nv) :: &
            res2, d_dlnRho_c_T2, d_dlnT_c_Rho2
         integer :: iregion, ix, iz
         character (len=256) :: message
   
         integer, parameter :: pure_helm = 1
         integer, parameter :: pure_opal_scvh = 2
         integer, parameter :: blend_in_x = 3
         integer, parameter :: blend_in_y = 4
         integer, parameter :: blend_corner_out = 5
         
         real(dp), parameter :: dZ_transition = 0.001d0 ! must be > 0
         real(dp) :: Pg, logPg, temp, logtemp, Prad, logW, alfa, beta, &
               Rho1, logRho1, dlnRho_dlnPgas_c_T1, dlnRho_dlnT_c_Pgas1, &
               Rho2, logRho2, dlnRho_dlnPgas_c_T2, dlnRho_dlnT_c_Pgas2, tiny
         type (EosPT_XZ_Info), pointer :: ep
         
         logical, parameter :: dbg = .false.

         include 'formats'
         
         ierr = 0
         tiny = rq% tiny_fuzz
         
         call load_single_eosPT_table(rq,ep,x,z,ix,iz,ierr)
         if (ierr /= 0) then
            if (dbg) write(*,*) 'Load_single_eosPE_Table ierr', ierr
            return
         end if
         
         logW = logPgas - 4d0*logT
         
         call get_alfa_beta(ierr)
         if (ierr /= 0) then
            if (dbg) write(*,*) 'error from get_alfa_beta'
            return
         end if
         
         if (beta /= 0D0) then
            ! before calling OPAL_SCVH, check that Z isn't too large. 
            if (Z > eosPT_Zs(num_eosPT_Zs)) then
               beta = 0d0
            else if (Z > eosPT_Zs(num_eosPT_Zs) - dZ_transition) then
               beta = beta*(eosPT_Zs(num_eosPT_Zs) - Z)/dZ_transition
            end if
            alfa = 1d0 - beta
         end if
         
         if (beta > 0.9999d0) then
            beta = 1d0; alfa = 0d0
         else if (alfa > 0.9999d0) then
            alfa = 1d0; beta = 0d0
         end if

         if (dbg) write(*,1) 'alfa', alfa
         if (dbg) write(*,1) 'beta', beta
         
         if (beta == 1d0) then
            call Get_PT_OPAL_SCVH_Results( &
               rq, Z, X, abar, zbar, Pgas, logPgas, T, logT, &
               Rho, logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas, &
               res, d_dlnRho_c_T, d_dlnT_c_Rho, ierr)
            return
         end if
         
         if (alfa == 1d0) then
            Pg = Pgas
            logPg = logPgas
            temp = T
            logtemp = logT
            call Get_PT_Results_using_DT( &
               rq, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               Pg, logPg, temp, logtemp, &
               Rho, logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas, &
               res, d_dlnRho_c_T, d_dlnT_c_Rho, ierr)
            return
         end if
               
         call Get_PT_OPAL_SCVH_Results( &
                  rq, Z, X, abar, zbar, Pgas, logPgas, T, logT, &
                  Rho2, logRho2, dlnRho_dlnPgas_c_T2, dlnRho_dlnT_c_Pgas2, &
                  res2, d_dlnRho_c_T2, d_dlnT_c_Rho2, ierr)
         if (dbg .and. logRho2 > rq% logRho2_OPAL_SCVH_limit) &
            write(*,*) 'logRho2 > rq% logRho2_OPAL_SCVH_limit', logRho2
         if (ierr /= 0 .or. logRho2 > rq% logRho2_OPAL_SCVH_limit) then
            ierr = 0
            beta = 0d0
            alfa = 1d0
         end if
      
         Pg = Pgas
         logPg = logPgas
         temp = T
         logtemp = logT
         call Get_PT_Results_using_DT( &
            rq, Z, X, abar, zbar, &
            species, chem_id, net_iso, xa, &
            Pg, logPg, temp, logtemp, &
            Rho1, logRho1, dlnRho_dlnPgas_c_T1, dlnRho_dlnT_c_Pgas1, &
            res1, d_dlnRho_c_T1, d_dlnT_c_Rho1, &
            ierr)
         if (ierr /= 0) then
            if (dbg) then
               write(*,*) 'error from Get_PT_Results_using_DT'
               stop
            end if
            return
         end if
      
         beta = 1d0 - alfa
         res = alfa*res1 + beta*res2
         logRho = alfa*logRho1 + beta*logRho2
         Rho = exp10(logRho)
         dlnRho_dlnPgas_c_T = alfa*dlnRho_dlnPgas_c_T1 + beta*dlnRho_dlnPgas_c_T2
         dlnRho_dlnT_c_Pgas = alfa*dlnRho_dlnT_c_Pgas1 + beta*dlnRho_dlnT_c_Pgas2
         d_dlnRho_c_T = alfa*d_dlnRho_c_T1 + beta*d_dlnRho_c_T2
         d_dlnT_c_Rho = alfa*d_dlnT_c_Rho1 + beta*d_dlnT_c_Rho2


         contains
         
         subroutine get_alfa_beta(ierr)
            integer, intent(out) :: ierr
            
            real(dp) :: logW_margin, logT_margin, &
               logW1, logW2, logW3, logW4, logT1, logT2, logT3, logT4, c_dx, c_dy
            
            include 'formats'
            
            c_dx = 0
            c_dy = 0

            logW_margin = 0.3d0
            logT_margin = 0.1d0
         
            logW1 = ep% logW_max
            logW2 = ep% logW_max - logW_margin
            logW3 = ep% logW_min + logW_margin
            logW4 = ep% logW_min
         
            logT1 = rq% logT_all_HELM
            logT2 = rq% logT_all_OPAL
            logT3 = ep% logT_min + logT_margin
            logT4 = ep% logT_min
                        
            if (dbg) write(*,1) 'logW', logW
            if (dbg) write(*,1) 'logW1', logW1
            if (dbg) write(*,1) 'logW2', logW2
            if (dbg) write(*,1) 'logW3', logW3
            if (dbg) write(*,1) 'logW4', logW4
            if (dbg) write(*,1) 'logT', logT
            if (dbg) write(*,1) 'logT1', logT1
            if (dbg) write(*,1) 'logT2', logT2
            if (dbg) write(*,1) 'logT3', logT3
            if (dbg) write(*,1) 'logT4', logT4
         
            if (logW <= logW4 .or. logW >= logW1 &
               .or. logT <= logT4 .or. logT >= logT1) then
               iregion = pure_helm
               if (dbg) write(*,*) 'case 1'
               if (dbg) write(*,*) 'logW <= logW4', logW <= logW4, logW, logW4
               if (dbg) write(*,*) 'logW >= logW1', logW >= logW1
               if (dbg) write(*,*) 'logT <= logT4', logT <= logT4
               if (dbg) write(*,*) 'logT >= logT1', logT >= logT1
            else if (logT > logT2) then
               c_dy = (logT - logT2) / (logT1 - logT2)
               if (logW > logW2) then
                  c_dx = (logW - logW2) / (logW1 - logW2)
                  iregion = blend_corner_out
                  if (dbg) write(*,*) 'case 2'
               else if (logW > logW3) then
                  iregion = blend_in_y
                  if (dbg) write(*,*) 'case 3'
               else ! logW > logW4
                  c_dx = (logW - logW3) / (logW4 - logW3)
                  iregion = blend_corner_out
                  if (dbg) write(*,*) 'case 4'
               end if
            else if (logT > logT3) then
               if (logW > logW2) then
                  c_dx = (logW - logW2) / (logW1 - logW2)
                  iregion = blend_in_x
                  if (dbg) write(*,*) 'case 5'
               else if (logW > logW3) then
                  iregion = pure_opal_scvh
                  if (dbg) write(*,*) 'case 6'
               else ! logW > logW4
                  c_dx = (logW - logW3) / (logW4 - logW3)
                  iregion = blend_in_x
                  if (dbg) write(*,*) 'case 7'
               end if
            else ! logT > logT4
               c_dy = (logT - logT3) / (logT4 - logT3)
               if (logW > logW2) then
                  c_dx = (logW - logW2) / (logW1 - logW2)
                  iregion = blend_corner_out
                  if (dbg) write(*,*) 'case 8'
               else if (logW > logW3) then
                  iregion = blend_in_y
                  if (dbg) write(*,*) 'case 9'
               else ! logW > logW4
                  c_dx = (logW - logW3) / (logW4 - logW3)
                  iregion = blend_corner_out
                  if (dbg) write(*,*) 'case 10'
               end if
            end if
         
            if (iregion == pure_helm) then
               alfa = 1d0
               beta = 0d0
            else if (iregion == pure_opal_scvh) then
               alfa = 0d0
               beta = 1d0
            else if (iregion == blend_in_y) then
               alfa = c_dy
               beta = 1d0 - alfa
            else if (iregion == blend_in_x) then
               alfa = c_dx
               beta = 1d0 - alfa
            else if (iregion == blend_corner_out) then
               alfa = min(1d0, sqrt(c_dx*c_dx + c_dy*c_dy))
               beta = 1d0 - alfa
            !else if (iregion == blend_corner_in) then
            !   beta = min(1d0, sqrt(c_dx*c_dx + c_dy*c_dy))
            !   alfa = 1 - beta
            else
               ierr = -1
               return
            end if
         end subroutine get_alfa_beta
         
      end subroutine Get1_eosPT_Results
      
      
      subroutine get_PT_args( &
            aPg, alogPg, atemp, alogtemp, Pgas, logPgas, T, logT, ierr)       
         real(dp), intent(in) :: aPg, alogPg
         real(dp), intent(in) :: atemp, alogtemp
         real(dp), intent(out) :: Pgas, logPgas, T, logT
         integer, intent(out) :: ierr
         include 'formats'
         ierr = 0
         T = atemp; logT = alogtemp
         if (atemp == arg_not_provided .and. alogtemp == arg_not_provided) then
            ierr = -2; return
         end if
         if (alogtemp == arg_not_provided) logT = log10(T)
         if (atemp == arg_not_provided) T = exp10(logT)
         if (T <= 0) then
            ierr = -1
            return
         end if
         Pgas = aPg; logPgas = alogPg
         if (Pgas == arg_not_provided .and. logPgas == arg_not_provided) then
            ierr = -3; return
         end if
         if (logPgas == arg_not_provided) logPgas = log10(Pgas)
         if (Pgas == arg_not_provided) Pgas = exp10(logPgas)
         if (Pgas <= 0) then
            ierr = -1
            return
         end if
      end subroutine get_PT_args

      
      subroutine Get_PT_OPAL_SCVH_Results( &
               rq, Z, X, abar, zbar, Pgas, logPgas, T, logT, &
               Rho, logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas, &
               res, d_dlnRho_c_T, d_dlnT_c_Rho, ierr)
         use chem_def
         type (EoS_General_Info), pointer :: rq ! general information about the request
         real(dp), intent(in) :: Z, X, abar, zbar, Pgas, logPgas, T, logT
         real(dp), intent(out) :: Rho, logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas
         real(dp), intent(inout) :: res(:) ! (nv)
         real(dp), intent(inout) :: d_dlnRho_c_T(:) ! (nv)
         real(dp), intent(inout) :: d_dlnT_c_Rho(:) ! (nv)
         integer, intent(out) :: ierr

         real(dp), dimension(nv, num_eosPT_Zs) :: &
            res_zx, d_dlnRho_c_T_zx, d_dlnT_c_Rho_zx, d_res_dX
         real(dp), dimension(nv) :: d_dX, d_dZ
         real(dp), dimension(3) :: c, &
            logRho_z, dlnRho_dlnPgas_c_T_z, dlnRho_dlnT_c_Pgas_z    
         real(dp) :: d_res_dZ(nv), dZ, denom, tiny

         character (len=256) :: message
         integer :: iz, j, ci
         logical, parameter :: OPAL_SCVH_dbg = .false.
         
         include 'formats'
         tiny = rq% tiny_fuzz

         ierr = 0
         
         if (num_eosPT_Zs < 3) then
            write(*, *) 'error: Get_PT_OPAL_SCVH_Results assumes num_eosPT_Zs >= 3'
            call mesa_error(__FILE__,__LINE__)
         end if
      
         if (eosPT_Zs(1) /= 0) then
            write(*, *) 'error: Get_PT_OPAL_SCVH_Results assumes eosPT_Zs(1) == 0'
            call mesa_error(__FILE__,__LINE__)
         end if
      
         if (abs(eosPT_Zs(1) - 2d0*eosPT_Zs(2) + eosPT_Zs(3)) > tiny) then
            write(*, *) 'error: Get_PT_OPAL_SCVH_Results assumes equal spaced Zs(1:3)'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         if (Z < 1d-20) then
            call Get_PT_OPAL_SCVH_for_X( &
               rq, 1, X, abar, zbar, Pgas, logPgas, T, logT, &
               logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas, &
               res, d_dlnRho_c_T, d_dlnT_c_Rho, ierr)
            Rho = exp10(logRho)
            return
         end if
                  
         do iz = 1, 3
            call Get_PT_OPAL_SCVH_for_X( &
               rq, iz, X, abar, zbar, Pgas, logPgas, T, logT, &
               logRho_z(iz), dlnRho_dlnPgas_c_T_z(iz), dlnRho_dlnT_c_Pgas_z(iz), &
               res_zx(:, iz), d_dlnRho_c_T_zx(:, iz), d_dlnT_c_Rho_zx(:, iz), ierr)
            if (ierr /= 0) return
         end do

         dZ = eosPT_Zs(2) - eosPT_Zs(1)
         denom = 2d0*dZ*dZ
         c(1) = (2d0*dZ*dZ - 3d0*dZ*Z + Z*Z)/denom
         c(2) = 2d0*(2d0*dZ-Z)*Z/denom
         c(3) = Z*(Z-dZ)/denom
            
         ! zero these for now
         res_zx(i_phase:i_latent_ddlnRho,:) = 0d0
         res_zx(i_frac:i_frac+num_eos_frac_results-1,:) = 0d0

         d_dlnT_c_Rho_zx(i_phase:i_latent_ddlnRho,:) = 0d0
         d_dlnT_c_Rho_zx(i_frac:i_frac+num_eos_frac_results-1,:) = 0d0

         d_dlnRho_c_T_zx(i_phase:i_latent_ddlnRho,:) = 0d0
         d_dlnRho_c_T_zx(i_frac:i_frac+num_eos_frac_results-1,:) = 0d0


         res(:) = c(1)*res_zx(:, 1) + c(2)*res_zx(:, 2) + c(3)*res_zx(:, 3)
         
         logRho = dot_product(c(:), logRho_z(:))
         dlnRho_dlnPgas_c_T = dot_product(c(:), dlnRho_dlnPgas_c_T_z(:))
         dlnRho_dlnT_c_Pgas = dot_product(c(:), dlnRho_dlnT_c_Pgas_z(:))
         Rho = exp10(logRho)
         
         d_dlnRho_c_T(:) = &
           c(1)*d_dlnRho_c_T_zx(:, 1) + &
           c(2)*d_dlnRho_c_T_zx(:, 2) + &
           c(3)*d_dlnRho_c_T_zx(:, 3)
     
         d_dlnT_c_Rho(:) = &
           c(1)*d_dlnT_c_Rho_zx(:, 1) + &
           c(2)*d_dlnT_c_Rho_zx(:, 2) + &
           c(3)*d_dlnT_c_Rho_zx(:, 3)
              
      end subroutine Get_PT_OPAL_SCVH_Results

      
      subroutine Get_PT_OPAL_SCVH_for_X( &
               rq, iz, X, abar, zbar, Pgas, logPgas, T, logT, &
               logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas, &
               res, d_dlnRho_c_T, d_dlnT_c_Rho, ierr)
         type (EoS_General_Info), pointer :: rq ! general information about the request
         integer, intent(in) :: iz ! the index in eosPT_Zs
         real(dp), intent(in) :: X, abar, zbar, Pgas, logPgas, T, logT
         real(dp), intent(out) :: logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas
         real(dp), dimension(:,:), pointer :: d_dx ! (nv, species)
         real(dp), intent(inout) :: res(:) ! (nv)
         real(dp), intent(inout) :: d_dlnRho_c_T(:) ! (nv)
         real(dp), intent(inout) :: d_dlnT_c_Rho(:) ! (nv)
         integer, intent(out) :: ierr

         real(dp), dimension(nv, 4) :: &
               res_zx, d_dlnRho_c_T_zx, d_dlnT_c_Rho_zx
         real(dp), dimension(4) :: logRho_zx, dlnRho_dlnPgas_c_T_zx, dlnRho_dlnT_c_Pgas_zx, c
         real(dp) :: dX, denom, delX, coef, tiny
         character (len=256) :: message
         integer :: ix, ix_lo, ix_hi, j, num_Xs
         logical, parameter :: dbg_for_X = .false.
         
         include 'formats'

         ierr = 0
         tiny = rq% tiny_fuzz
         
         if (num_eosPT_Xs <= 5) then
            write(*, *) 'error: Get_PT_OPAL_SCVH_for_X assumes num_eosPT_Xs > 5'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         if (eosPT_Xs(1) /= 0) then
            write(*, *) 'error: Get_PT_OPAL_SCVH_for_X assumes eosPT_Xs(1) == 0'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         num_Xs = num_eosPT_Xs_for_Z(iz)
         
         if (X < tiny .or. num_Xs == 1) then
            call Get_eosPT_XTable_Results(rq, 1, iz, Pgas, logPgas, T, logT, &
               logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas, res, d_dlnRho_c_T, d_dlnT_c_Rho, &
               ierr)
            if (dbg_for_X) write(*,1) 'dbg_for_X logRho', logRho
            return
         end if
         
         dX = eosPT_Xs(2)-eosPT_Xs(1)
         
         do ix = 3, num_Xs
            if (abs(dX - (eosPT_Xs(ix) - eosPT_Xs(ix-1))) > tiny) then
               write(*, *) 'error: Get_PT_OPAL_SCVH_for_X assumes equal spaced Xs'
               call mesa_error(__FILE__,__LINE__)
            end if
         end do
         
         ix_hi = -1
         ix_lo = -1
         if (X <= eosPT_Xs(2)) then
            ix_lo = 1; ix_hi = 3
         else if (X >= eosPT_Xs(num_Xs-1)) then
            ix_lo = num_Xs-2; ix_hi = num_Xs
         else
            do ix = 3, num_Xs-1
               if (X <= eosPT_Xs(ix)) then
                  ix_lo = ix-2; ix_hi = ix+1; exit
               end if
            end do
         end if
         
         if (ix_hi < 0) then
            write(*, *) 'error: Get_PT_OPAL_SCVH_for_X logic bug'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         if (dbg_for_X) then
            write(*, *) 'X', X
            write(*, *) 'ix_lo', ix_lo
            write(*, *) 'ix_hi', ix_hi
         end if
         
         do ix=ix_lo, ix_hi
            j = ix-ix_lo+1
            call Get_eosPT_XTable_Results(rq, ix, iz, Pgas, logPgas, T, logT, &
               logRho_zx(j), dlnRho_dlnPgas_c_T_zx(j), dlnRho_dlnT_c_Pgas_zx(j), &
               res_zx(:, j), d_dlnRho_c_T_zx(:, j), d_dlnT_c_Rho_zx(:, j), ierr)
            if (ierr /= 0) return
         end do

         ! zero these for now
         res_zx(i_phase:i_latent_ddlnRho,:) = 0d0
         res_zx(i_frac:i_frac+num_eos_frac_results-1,:) = 0d0

         d_dlnRho_c_T_zx(i_phase:i_latent_ddlnRho,:) = 0d0
         d_dlnRho_c_T_zx(i_frac:i_frac+num_eos_frac_results-1,:) = 0d0

         d_dlnT_c_Rho_zx(i_phase:i_latent_ddlnRho,:) = 0d0
         d_dlnT_c_Rho_zx(i_frac:i_frac+num_eos_frac_results-1,:) = 0d0
         
         delX = X - eosPT_Xs(ix_lo)
         if (ix_hi-ix_lo==2) then
         
            denom = 2*dX*dX
            c(1) = (2*dX*dX - 3*dX*delX + delX*delX)/denom
            c(2) = 2*(2*dX-delX)*delX/denom
            c(3) = delX*(delX-dX)/denom
            res(:) = c(1)*res_zx(:, 1) + c(2)*res_zx(:, 2) + c(3)*res_zx(:, 3)
            logRho = dot_product(c(1:3), logRho_zx(1:3))
            dlnRho_dlnPgas_c_T = dot_product(c(1:3), dlnRho_dlnPgas_c_T_zx(1:3))
            dlnRho_dlnT_c_Pgas = dot_product(c(1:3), dlnRho_dlnT_c_Pgas_zx(1:3))
            
            d_dlnRho_c_T(:) = &
               c(1)*d_dlnRho_c_T_zx(:, 1) + &
               c(2)*d_dlnRho_c_T_zx(:, 2) + &
               c(3)*d_dlnRho_c_T_zx(:, 3)
            d_dlnT_c_Rho(:) = &
               c(1)*d_dlnT_c_Rho_zx(:, 1) + &
               c(2)*d_dlnT_c_Rho_zx(:, 2) + &
               c(3)*d_dlnT_c_Rho_zx(:, 3)
     
            c(1) = (2*delX-3*dX)/denom
            c(2) = 4*(dX-delX)/denom
            c(3) = (2*delX-dX)/denom
            
         else
         
            coef = (X-eosPT_Xs(ix_lo+1))/dX 
            ! coef = fractional location of X between 2nd and 3rd X's for fit.
            ! coef is weight for the quadratic based on points 2, 3, 4 of fit.
            ! (1-coef) is weight for quadratic based on points 1, 2, 3 of fit.
            if (coef < 0 .or. coef > 1) then
               ierr = -1
               return
            end if
            
            c(1) = -coef*(coef-1)*(coef-1)/2
            c(2) = (2 - coef*coef*(5 - 3*coef))/2
            c(3) = coef*(1 + coef*(4 - 3*coef))/2
            c(4) = coef*coef*(coef-1)/2
            res(:) = c(1)*res_zx(:, 1) + &
                        (c(2)*res_zx(:, 2) + &
                           (c(3)*res_zx(:, 3) + &
                              c(4)*res_zx(:, 4)))
     
            logRho = dot_product(c(1:4), logRho_zx(1:4))
            dlnRho_dlnPgas_c_T = dot_product(c(1:4), dlnRho_dlnPgas_c_T_zx(1:4))
            dlnRho_dlnT_c_Pgas = dot_product(c(1:4), dlnRho_dlnT_c_Pgas_zx(1:4))
            
            d_dlnRho_c_T(:) = &
               c(1)*d_dlnRho_c_T_zx(:, 1) + &
                  (c(2)*d_dlnRho_c_T_zx(:, 2) + &
                     (c(3)*d_dlnRho_c_T_zx(:, 3) + &
                           c(4)*d_dlnRho_c_T_zx(:, 4)))
            d_dlnT_c_Rho(:) = &
               c(1)*d_dlnT_c_Rho_zx(:, 1) + &
                  (c(2)*d_dlnT_c_Rho_zx(:, 2) + &
                     (c(3)*d_dlnT_c_Rho_zx(:, 3) + &
                           c(4)*d_dlnT_c_Rho_zx(:, 4)))
     
            denom = 2*dX
            c(1) = (-1 + coef*(3 - 2*coef))/denom
            c(2) = coef*(-7 + 6*coef)/denom
            c(3) = (1 + coef*(5 - 6*coef))/denom
            c(4) = coef*(-1 + 2*coef)/denom
     
         end if
         
      end subroutine Get_PT_OPAL_SCVH_for_X
         
      
      subroutine Get_eosPT_XTable_Results( &
               rq, ix, iz, Pgas, logPgas_in, T, logT_in, &
               logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas, &
               res, d_dlnRho_c_T, d_dlnT_c_Rho, ierr)
         use eosPT_load_tables, only : load_single_eosPT_table_by_id
         type (EoS_General_Info), pointer :: rq
         integer, intent(in) :: ix, iz
         real(dp), intent(in) :: Pgas, logPgas_in, T, logT_in
         real(dp), intent(out) :: logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas
         real(dp), intent(inout) :: res(:) ! (nv)
         real(dp), intent(inout) :: d_dlnRho_c_T(:) ! (nv)
         real(dp), intent(inout) :: d_dlnT_c_Rho(:) ! (nv)
            ! only does partials wrt lnRho and lnT
         integer, intent(out) :: ierr
         
         real(dp) :: f(nv), df_dlogW_c_T(nv), df_dlogT_c_W(nv), &
            df_dlogRho_c_T(nv), df_dlogT_c_Rho(nv)
         real(dp) :: logW0, logW1, logT0, logT1
         real(dp) :: logPgas, logT, logW, Prad, P, dlnPgas_dlnRho_c_T, dlnPgas_dlnT_c_Rho, &
               dlnW_dlnd, dlnW_dlnT, dlogW_dlogRho_c_T, dlogW_dlogT_c_Rho 
         integer :: iW, jtemp
         type (EosPT_XZ_Info), pointer :: ep
         
         character (len=256) :: message
         
         include 'formats'

         logPgas = logPgas_in
         logT = logT_in
         logW = logPgas - 4*logT

         ierr = 0
         
         call load_single_eosPT_table_by_id(rq,ep,ix,iz,ierr)
         if (ierr /= 0) return
         
         call Locate_logW(rq, ep, logW, iW, logW0, logW1, ierr)
         if (ierr /= 0) return
         
         call Locate_logT(rq, ep, logT, jtemp, logT0, logT1, ierr)
         if (ierr /= 0) return

         call Do_EoS_Interpolations( &
                  nv, ep% num_logWs, ep% logWs, ep% num_logTs, ep% logTs, &
                  ep% tbl1, iW, jtemp, logW0, logW, logW1, logT0, logT, logT1, &
                  f, df_dlogW_c_T, df_dlogT_c_W, ierr)
         if (ierr /= 0) return
            
         logRho = f(eosPT_ilnRho)/ln10
         
         ! logRho = logRho(logW(logPgas,logT),logT)
         
         ! dlnRho/dlnPgas|T = dlogRho/dlogPgas|T
         !     = dlogRho/dlogW|T * dlogW/dlogPgas|T
         !     = dlnRho/dlogW|T * dlogRho/dlnRho * dlogW/dlogPgas|T
         ! dlnRho/dlogW|T = df_dlogW_c_T(eosPT_ilnRho)
         ! dlogRho/dlnRho = 1/ln10
         ! dlogW/dlogPgas|T = 1
         dlnRho_dlnPgas_c_T = df_dlogW_c_T(eosPT_ilnRho)/ln10
         
         ! dlnRho_dlnT|Pgas = dlogRho/dlogT|Pgas
         !     = dlogRho/dlogW|Pgas * dlogW/dlogT|Pgas + dlogRho/dlogT|W
         !     = (dlnRho/dlogW|Pgas * dlogW/dlogT|Pgas + dlnRho/dlogT|W)/ln10
         ! dlnRho/dlogW|Pgas = df_dlogW_c_T(eosPT_ilnRho)
         ! dlogW/dlogT|Pgas = -4
         ! dlnRho/dlogT|W = df_dlogT_c_W(eosPT_ilnRho)
         dlnRho_dlnT_c_Pgas = &
               (-4*df_dlogW_c_T(eosPT_ilnRho) + df_dlogT_c_W(eosPT_ilnRho))/ln10
         
         Prad = crad*T*T*T*T/3
         P = Pgas + Prad
         
         ! dlnPgas_dlnT|Rho = (T/Pgas)*dPgas_dT|Rho
         !     = (T/Pgas)*(dP_dT|Rho - dPrad_dT)
         !     = (T/Pgas)*((P/T)dlnP_dlnT|Rho - dPrad_dT)
         !     = (P/Pgas)*dlnP_dlnT|Rho - (T/Pgas)*dPrad_dT
         !     = (P*dlnP_dlnT|Rho - T*dPrad_dT)/Pgas
         dlnPgas_dlnT_c_Rho = (P*f(i_chiT) - 4*Prad/3)/Pgas         
         dlnPgas_dlnRho_c_T = f(i_chiRho)
         
         ! logW = logPgas(logRho,logT) - 4*logT
         dlogW_dlogRho_c_T = dlnPgas_dlnRho_c_T
         dlogW_dlogT_c_Rho = dlnPgas_dlnT_c_Rho - 4
         
         ! now let f = f(logW(logRho,logT),logT)
         df_dlogRho_c_T = df_dlogW_c_T*dlogW_dlogRho_c_T
         df_dlogT_c_Rho = df_dlogW_c_T*dlogW_dlogT_c_Rho + df_dlogT_c_W

         res(i_lnPgas) = logPgas*ln10
         res(i_lnE) = f(i_lnE)
         res(i_lnS) = f(i_lnS)
         res(i_grad_ad) = f(i_grad_ad)
         res(i_chiRho) = f(i_chiRho)
         res(i_chiT) = f(i_chiT)
         res(i_Cp) = f(i_Cp)
         res(i_Cv) = f(i_Cv)
         res(i_dE_dRho) = f(i_dE_dRho)
         res(i_dS_dT) = f(i_dS_dT)
         res(i_dS_dRho) = f(i_dS_dRho)
         res(i_mu) = f(i_mu)
         res(i_lnfree_e) = f(i_lnfree_e)
         res(i_gamma1) = f(i_gamma1)
         res(i_gamma3) = f(i_gamma3)
         res(i_eta) = f(i_eta)
         
         d_dlnRho_c_T(i_lnPgas) = dlnPgas_dlnRho_c_T
         d_dlnRho_c_T(i_lnE) = df_dlogRho_c_T(i_lnE)/ln10
         d_dlnRho_c_T(i_lnS) = df_dlogRho_c_T(i_lnS)/ln10
         d_dlnRho_c_T(i_grad_ad) = df_dlogRho_c_T(i_grad_ad)/ln10
         d_dlnRho_c_T(i_chiRho) = df_dlogRho_c_T(i_chiRho)/ln10
         d_dlnRho_c_T(i_chiT) = df_dlogRho_c_T(i_chiT)/ln10
         d_dlnRho_c_T(i_Cp) = df_dlogRho_c_T(i_Cp)/ln10
         d_dlnRho_c_T(i_Cv) = df_dlogRho_c_T(i_Cv)/ln10
         d_dlnRho_c_T(i_dE_dRho) = df_dlogRho_c_T(i_dE_dRho)/ln10
         d_dlnRho_c_T(i_dS_dT) = df_dlogRho_c_T(i_dS_dT)/ln10
         d_dlnRho_c_T(i_dS_dRho) = df_dlogRho_c_T(i_dS_dRho)/ln10
         d_dlnRho_c_T(i_mu) = df_dlogRho_c_T(i_mu)/ln10
         d_dlnRho_c_T(i_lnfree_e) = df_dlogRho_c_T(i_lnfree_e)/ln10
         d_dlnRho_c_T(i_gamma1) = df_dlogRho_c_T(i_gamma1)/ln10
         d_dlnRho_c_T(i_gamma3) = df_dlogRho_c_T(i_gamma3)/ln10
         d_dlnRho_c_T(i_eta) = df_dlogRho_c_T(i_eta)/ln10
      
         d_dlnT_c_Rho(i_lnPgas) = dlnPgas_dlnT_c_Rho
         d_dlnT_c_Rho(i_lnE) = df_dlogT_c_Rho(i_lnE)/ln10
         d_dlnT_c_Rho(i_lnS) = df_dlogT_c_Rho(i_lnS)/ln10
         d_dlnT_c_Rho(i_grad_ad) = df_dlogT_c_Rho(i_grad_ad)/ln10
         d_dlnT_c_Rho(i_chiRho) = df_dlogT_c_Rho(i_chiRho)/ln10
         d_dlnT_c_Rho(i_chiT) = df_dlogT_c_Rho(i_chiT)/ln10
         d_dlnT_c_Rho(i_Cp) = df_dlogT_c_Rho(i_Cp)/ln10
         d_dlnT_c_Rho(i_Cv) = df_dlogT_c_Rho(i_Cv)/ln10
         d_dlnT_c_Rho(i_dE_dRho) = df_dlogT_c_Rho(i_dE_dRho)/ln10
         d_dlnT_c_Rho(i_dS_dT) = df_dlogT_c_Rho(i_dS_dT)/ln10
         d_dlnT_c_Rho(i_dS_dRho) = df_dlogT_c_Rho(i_dS_dRho)/ln10
         d_dlnT_c_Rho(i_mu) = df_dlogT_c_Rho(i_mu)/ln10
         d_dlnT_c_Rho(i_lnfree_e) = df_dlogT_c_Rho(i_lnfree_e)/ln10
         d_dlnT_c_Rho(i_gamma1) = df_dlogT_c_Rho(i_gamma1)/ln10
         d_dlnT_c_Rho(i_gamma3) = df_dlogT_c_Rho(i_gamma3)/ln10
         d_dlnT_c_Rho(i_eta) = df_dlogT_c_Rho(i_eta)/ln10

      end subroutine Get_eosPT_XTable_Results
      
      
      subroutine Locate_logW( &
            rq, ep, logW, iW, logW0, logW1, ierr)
         type (EoS_General_Info), pointer :: rq
         type (EosPT_XZ_Info), pointer :: ep
         real(dp), intent(inout) :: logW
         integer, intent(out) :: iW
         real(dp), intent(out) :: logW0, logW1
         integer, intent(out) :: ierr
         
         logical, parameter :: dbg = .false.

         include 'formats'
      
         ierr = 0
         iW = int((logW - ep% logW_min) / ep% del_logW + 1d-4) + 1
         
         if (iW < 1 .or. iW >= ep% num_logWs) then
            
            if (iW < 1) then
               iW = 1
               logW0 = ep% logW_min
               logW1 = logW0 + ep% del_logW
               logW = logW0
            else
               iW = ep% num_logWs-1
               logW0 = ep% logW_min + (iW-1) * ep% del_logW
               logW1 = logW0 + ep% del_logW
               logW = logW1
            end if
            
         else
         
            logW0 = ep% logW_min + (iW-1) * ep% del_logW
            logW1 = logW0 + ep% del_logW

         end if
         
         if (dbg) then
            write(*,*) 'Locate_logW iW', iW
            write(*,1) 'logW0', logW0
            write(*,1) 'logW', logW
            write(*,1) 'logW1', logW1
            write(*,1) 'ep% logW_min', ep% logW_min
            write(*,1) 'ep% del_logW', ep% del_logW
            write(*,*)
         end if

      end subroutine Locate_logW
      
      
      subroutine Locate_logT( &
            rq, ep, logT, iT, logT0, logT1, ierr)
         type (EoS_General_Info), pointer :: rq
         type (EosPT_XZ_Info), pointer :: ep
         real(dp), intent(inout) :: logT
         integer, intent(out) :: iT
         real(dp), intent(out) :: logT0, logT1
         integer, intent(out) :: ierr
         
         logical, parameter :: dbg = .false.

         include 'formats'
      
         ierr = 0
         iT = int((logT - ep% logT_min) / ep% del_logT + 1d-4) + 1
         
         if (iT < 1 .or. iT >= ep% num_logTs) then
            
            if (iT < 1) then
               iT = 1
               logT0 = ep% logT_min
               logT1 = logT0 + ep% del_logT
               logT = logT0
            else
               iT = ep% num_logTs-1
               logT0 = ep% logT_min + (iT-1) * ep% del_logT
               logT1 = logT0 + ep% del_logT
               logT = logT1
            end if
            
         else
         
            logT0 = ep% logT_min + (iT-1) * ep% del_logT
            logT1 = logT0 + ep% del_logT

         end if
         
         if (dbg) then
            write(*,*) 'Locate_logT iT', iT
            write(*,1) 'logT0', logT0
            write(*,1) 'logT', logT
            write(*,1) 'logT1', logT1
            write(*,1) 'ep% logT_min', ep% logT_min
            write(*,1) 'ep% del_logT', ep% del_logT
            write(*,*)
         end if

      end subroutine Locate_logT


      subroutine Do_EoS_Interpolations( &
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
         
         include 'formats'
         
         ierr = 0
         
         fin(1:4,1:nv,1:nx,1:ny) => fin1(1:4*nv*nx*ny)
         
         hx=x1-x0
         hxi=1.0d0/hx
         hx2=hx*hx
   
         xp=(xget-x0)*hxi
         xpi=1.0d0-xp
         xp2=xp*xp
         xpi2=xpi*xpi

         ax=xp2*(3.0d0-2.0d0*xp)
         axbar=1.d0-ax
         
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
         aybar=1.0d0-ay
         
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
         
      end subroutine Do_EoS_Interpolations
         

      subroutine Get_PT_Results_using_DT( &
               rq, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               Pgas, logPgas, T, logT, &
               Rho, logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas, &
               res, d_dlnRho_c_T, d_dlnT_c_Rho, ierr)
         use eosDT_eval, only: get_Rho
         use utils_lib, only: is_bad
         
         type (EoS_General_Info), pointer :: rq ! general information about the request
         real(dp), intent(in) :: Z, X, abar, zbar
         integer, intent(in) :: species
         integer, pointer :: chem_id(:)    
         integer, pointer :: net_iso(:)
         real(dp), intent(in) :: xa(:)
         real(dp), intent(inout) :: Pgas, logPgas, T, logT
         real(dp), intent(out) :: Rho, logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas
         real(dp), intent(inout) :: res(:) ! (nv)
         real(dp), intent(inout) :: d_dlnRho_c_T(:) ! (nv)
         real(dp), intent(inout) :: d_dlnT_c_Rho(:) ! (nv)
         integer, intent(out) :: ierr
         
         logical, parameter :: basic_flag = .false.
         integer:: i, eos_calls, max_iter, which_other
         real(dp) :: &
            logRho_guess, rho_guess, other, other_tol, logRho_tol, Prad, f, dfdx, &
            logRho_bnd1, logRho_bnd2, other_at_bnd1, other_at_bnd2, logRho_result
         
         logical, parameter :: dbg = .false.
         
         include 'formats'
         
         ierr = 0

         which_other = i_lnPgas
         other = logPgas*ln10
         other_tol = 1d-8
         logRho_tol = 1d-8
         
         ! guess based on fully ionized, ideal gas of ions and electrons
         rho_guess = Pgas*abar*mp/(kerg*T*(1+zbar))
         logRho_guess = log10(rho_guess)
      
         logRho_bnd1 = arg_not_provided
         logRho_bnd2 = arg_not_provided
         other_at_bnd1 = arg_not_provided
         other_at_bnd2 = arg_not_provided

         max_iter = 20
         eos_calls = 0
         
         if (dbg) write(*,1) 'rho_guess', rho_guess
         if (dbg) write(*,1) 'logRho_guess', logRho_guess
         
         call get_Rho( &
               rq% handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               logT, which_other, other, &
               logRho_tol, other_tol, max_iter, logRho_guess, &
               logRho_bnd1, logRho_bnd2, other_at_bnd1, other_at_bnd2, &
               logRho_result, res, d_dlnRho_c_T, d_dlnT_c_Rho, &
               eos_calls, ierr)
         if (ierr /= 0) then
            if (.false.) then
               write(*,*) 'failed in get_Rho for Get_PT_Results_using_DT'
               write(*,1) 'Z = ', Z
               write(*,1) 'X = ', X
               write(*,1) 'abar = ', abar
               write(*,1) 'zbar = ', zbar
               write(*,1) 'logT = ', logT
               write(*,1) 'Pgas = ', Pgas
               write(*,1) 'logPgas = ', logPgas
               write(*,1) 'logRho_tol = ', logRho_tol
               write(*,1) 'other_tol = ', other_tol
               write(*,1) 'logRho_guess = ', logRho_guess
               write(*,*)
            end if
            return
         end if
         
         logRho = logRho_result
         Rho = exp10(logRho)
         
         if (dbg) write(*,1) 'Rho', Rho
         if (dbg) write(*,1) 'logRho', logRho
         if (dbg) write(*,*)
         if (dbg) write(*,1) 'Pgas input', Pgas
         if (dbg) write(*,1) 'logPgas input', logPgas
         if (dbg) write(*,1) 'Pgas match', exp(res(i_lnPgas))
         if (dbg) write(*,1) 'logPgas match', res(i_lnPgas)/ln10
         if (dbg) write(*,*)
         if (dbg) write(*,1) 'get_Rho: grad_ad', res(i_grad_ad)
         if (dbg) write(*,*)
         
         call do_partials
         
         contains
         
         subroutine do_partials ! dlnRho_dlnPgas_c_T and dlnRho_dlnT_c_Pgas
            real(dp) :: Prad, P, dP_dRho, dPgas_dRho, &
                  dP_dT, dPrad_dT, dPgas_dT, dRho_dPgas, dRho_dT
            include 'formats'
            
            Prad = crad*T*T*T*T/3
            P = Pgas + Prad
            dP_dRho = res(i_chiRho)*P/Rho
            dPgas_dRho = dP_dRho ! const T, so dP_dRho = dPgas_dRho
            dRho_dPgas = 1/dPgas_dRho ! const T
            dlnRho_dlnPgas_c_T = dRho_dPgas*Pgas/Rho ! const T
            
            dPrad_dT = 4*crad*T*T*T/3
            dP_dT = res(i_chiT)*P/T
            dPgas_dT = dP_dT - dPrad_dT ! const Rho
            dRho_dT = -dPgas_dT/dPgas_dRho ! const Pgas
            dlnRho_dlnT_c_Pgas = dRho_dT*T/Rho
            
         end subroutine do_partials

      end subroutine Get_PT_Results_using_DT


      subroutine get_T( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               logPgas, which_other, other_value, &
               logT_tol, other_tol, max_iter, logT_guess, &
               logT_bnd1, logT_bnd2,  other_at_bnd1, other_at_bnd2, &
               logT_result, Rho, logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas, &
               res, d_dlnRho_c_T, d_dlnT_c_Rho, &
               eos_calls, ierr)
         
         integer, intent(in) :: handle

         real(dp), intent(in) :: Z ! the metals mass fraction
         real(dp), intent(in) :: X ! the hydrogen mass fraction
            
         real(dp), intent(in) :: abar, zbar
         
         integer, intent(in) :: species
         integer, pointer :: chem_id(:)    
         integer, pointer :: net_iso(:)
         real(dp), intent(in) :: xa(:)
         
         real(dp), intent(in) :: logPgas ! log10 of density
         integer, intent(in) :: which_other
         real(dp), intent(in) :: other_value ! desired value for the other variable
         real(dp), intent(in) :: other_tol
         
         real(dp), intent(in) :: logT_tol
         integer, intent(in) :: max_iter ! max number of iterations        

         real(dp), intent(in) :: logT_guess
         real(dp), intent(in) :: logT_bnd1, logT_bnd2 ! bounds for logT
            ! set to arg_not_provided if do not know bounds
         real(dp), intent(in) :: other_at_bnd1, other_at_bnd2 ! values at bounds
            ! if don't know these values, just set to arg_not_provided (defined in c_def)
         
         real(dp), intent(out) :: logT_result
         real(dp), intent(out) :: Rho, logRho ! density
         real(dp), intent(out) :: dlnRho_dlnPgas_c_T
         real(dp), intent(out) :: dlnRho_dlnT_c_Pgas

         real(dp), intent(inout) :: res(:) ! (nv)
         real(dp), intent(inout) :: d_dlnRho_c_T(:) ! (nv)
         real(dp), intent(inout) :: d_dlnT_c_Rho(:) ! (nv)
         
         integer, intent(out) :: eos_calls
         integer, intent(out) :: ierr ! 0 means AOK.
         
         call do_safe_get_Pgas_T( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               logPgas, which_other, other_value, doing_get_T, &
               logT_guess, logT_result, logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
               logT_tol, other_tol, max_iter, &
               Rho, logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas, &
               res, d_dlnRho_c_T, d_dlnT_c_Rho, &
               eos_calls, ierr)
      
      end subroutine get_T
      

      subroutine get_Pgas( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               logT, which_other, other_value, &
               logPgas_tol, other_tol, max_iter, logPgas_guess, &
               logPgas_bnd1, logPgas_bnd2, other_at_bnd1, other_at_bnd2, &
               logPgas_result, Rho, logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas, &
               res, d_dlnRho_c_T, d_dlnT_c_Rho, &
               eos_calls, ierr)
     
         use const_def
         
         integer, intent(in) :: handle

         real(dp), intent(in) :: Z ! the metals mass fraction
         real(dp), intent(in) :: X ! the hydrogen mass fraction
            
         real(dp), intent(in) :: abar, zbar
         
         integer, intent(in) :: species
         integer, pointer :: chem_id(:)    
         integer, pointer :: net_iso(:)
         real(dp), intent(in) :: xa(:)
         
         real(dp), intent(in) :: logT ! log10 of temperature

         integer, intent(in) :: which_other
         real(dp), intent(in) :: other_value ! desired value for the other variable
         real(dp), intent(in) :: other_tol
         
         real(dp), intent(in) :: logPgas_tol

         integer, intent(in) :: max_iter ! max number of Newton iterations        

         real(dp), intent(in) :: logPgas_guess
         real(dp), intent(in) :: logPgas_bnd1, logPgas_bnd2 ! bounds for logPgas
            ! set to arg_not_provided if do not know bounds
         real(dp), intent(in) :: other_at_bnd1, other_at_bnd2 ! values at bounds
            ! if don't know these values, just set to arg_not_provided (defined in c_def)
            
         real(dp), intent(out) :: logPgas_result
         real(dp), intent(out) :: Rho, logRho ! density
         real(dp), intent(out) :: dlnRho_dlnPgas_c_T
         real(dp), intent(out) :: dlnRho_dlnT_c_Pgas

         real(dp), intent(inout) :: res(:) ! (nv)
         real(dp), intent(inout) :: d_dlnRho_c_T(:) ! (nv)
         real(dp), intent(inout) :: d_dlnT_c_Rho(:) ! (nv)

         integer, intent(out) :: eos_calls
         integer, intent(out) :: ierr ! 0 means AOK.
         
         call do_safe_get_Pgas_T( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               logT, which_other, other_value, doing_get_Pgas, &
               logPgas_guess, logPgas_result, logPgas_bnd1, logPgas_bnd2, other_at_bnd1, other_at_bnd2, &
               logPgas_tol, other_tol, max_iter, &
               Rho, logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas, &
               res, d_dlnRho_c_T, d_dlnT_c_Rho, &
               eos_calls, ierr)

      end subroutine get_Pgas
      

      subroutine get_Pgas_for_Rho( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               logT, logRho_want, &
               logPgas_tol, logRho_tol, max_iter, logPgas_guess_in, &
               logPgas_bnd1, logPgas_bnd2, logRho_at_bnd1, logRho_at_bnd2, &
               logPgas_result, Rho, logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas, &
               res, d_dlnRho_c_T, d_dlnT_c_Rho, &
               eos_calls, ierr)

         integer, intent(in) :: handle

         real(dp), intent(in) :: Z ! the metals mass fraction
         real(dp), intent(in) :: X ! the hydrogen mass fraction
            
         real(dp), intent(in) :: abar, zbar
         
         integer, intent(in) :: species
         integer, pointer :: chem_id(:)    
         integer, pointer :: net_iso(:)
         real(dp), intent(in) :: xa(:)
         
         real(dp), intent(in) :: logT ! log10 of temperature

         real(dp), intent(in) :: logRho_want ! log10 of desired density
         real(dp), intent(in) :: logRho_tol
         
         real(dp), intent(in) :: logPgas_tol

         integer, intent(in) :: max_iter ! max number of Newton iterations        

         real(dp), intent(in) :: logPgas_guess_in ! log10 of gas pressure
            ! if = arg_not_provided, then will use ideal gas for guess
         real(dp), intent(in) :: logPgas_bnd1, logPgas_bnd2 ! bounds for logPgas
            ! if don't know bounds, just set to arg_not_provided (defined in const_def)
         real(dp), intent(in) :: logRho_at_bnd1, logRho_at_bnd2 ! values at bounds
            ! if don't know these values, just set to arg_not_provided (defined in const_def)

         real(dp), intent(out) :: logPgas_result ! log10 of gas pressure
         real(dp), intent(out) :: Rho, logRho ! density corresponding to logPgas_result
         real(dp), intent(out) :: dlnRho_dlnPgas_c_T
         real(dp), intent(out) :: dlnRho_dlnT_c_Pgas

         real(dp), intent(inout) :: res(:) ! (nv)
         real(dp), intent(inout) :: d_dlnRho_c_T(:) ! (nv)
         real(dp), intent(inout) :: d_dlnT_c_Rho(:) ! (nv)

         integer, intent(out) :: eos_calls
         integer, intent(out) :: ierr ! 0 means AOK.

         real(dp) :: logPgas_guess
         
         logPgas_guess = logPgas_guess_in
         if (logPgas_guess == arg_not_provided) then
            ! Pgas = rho*kerg*T*(1+zbar)/(abar*mp)
            logPgas_guess = logRho_want + logT + log10(kerg*(1+zbar)/(abar*mp))
         end if

         call do_safe_get_Pgas_T( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               logT, 0, logRho_want, doing_get_Pgas_for_Rho, &
               logPgas_guess, logPgas_result, logPgas_bnd1, logPgas_bnd2, &
               logRho_at_bnd1, logRho_at_bnd2, logPgas_tol, logRho_tol, max_iter, &
               Rho, logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas, &
               res, d_dlnRho_c_T, d_dlnT_c_Rho, &
               eos_calls, ierr)
      
      end subroutine get_Pgas_for_Rho

      
      subroutine do_safe_get_Pgas_T( &
               handle, Z, XH1, abar, zbar, &
               species, chem_id, net_iso, xa, &
               the_other_log, which_other, other_value, doing_which, &
               initial_guess, x, xbnd1, xbnd2, other_at_bnd1, other_at_bnd2, &
               xacc, yacc, ntry, &
               Rho, logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas, &
               res, d_dlnRho_c_T, d_dlnT_c_Rho, &
               eos_calls, ierr)
         use const_def
         use utils_lib, only: is_bad
         use num_lib, only: brent_safe_zero, look_for_brackets
         use chem_def, only: num_chem_isos
         use eosPT_load_tables, only: load_single_eosPT_table
         integer, intent(in) :: handle
         real(dp), intent(in) :: Z, XH1, abar, zbar
         integer, intent(in) :: species
         integer, pointer :: chem_id(:)    
         integer, pointer :: net_iso(:)
         real(dp), intent(in) :: xa(:)
         integer, intent(in) :: which_other
         real(dp), intent(in) :: other_value
         integer, intent(in) :: doing_which
         real(dp), intent(in) :: initial_guess ! for x
         real(dp), intent(out) :: x ! if doing_Pgas, then logPgas, else logT
         real(dp), intent(in) :: the_other_log
         real(dp), intent(in) :: xbnd1, xbnd2, other_at_bnd1, other_at_bnd2
         real(dp), intent(in) :: xacc, yacc ! tolerances
         integer, intent(in) :: ntry ! max number of iterations        
         real(dp), intent(out) :: Rho, logRho ! density
         real(dp), intent(out) :: dlnRho_dlnPgas_c_T
         real(dp), intent(out) :: dlnRho_dlnT_c_Pgas
         real(dp), intent(inout) :: res(:) ! (nv)
         real(dp), intent(inout) :: d_dlnRho_c_T(:) ! (nv)
         real(dp), intent(inout) :: d_dlnT_c_Rho(:) ! (nv)
         integer, intent(out) :: eos_calls, ierr
         
         integer :: i, j, lrpar, lipar, max_iter, irho, ix, iz
         integer, parameter :: lrextras=4
         real(dp), parameter :: dx = 0.1d0
         integer, pointer :: ipar(:)
         real(dp), pointer :: rpar(:)
         real(dp) :: Pgas, T, xb1, xb3, y1, y3, dfdx, f, logPgas, logT
         logical, parameter :: basic_flag = .false.
         type (EoS_General_Info), pointer :: rq
         type (EosPT_XZ_Info), pointer :: ep
         
         logical, parameter :: dbg = .false.
         
         include 'formats'
         
         ierr = 0

         call get_eos_ptr(handle, rq, ierr)
         if (ierr /= 0) then
            write(*, *) 'get_eos_ptr returned ierr', ierr
            return
         end if

         call load_single_eosPT_table(rq,ep,XH1,z,ix,iz,ierr)
         if (ierr /= 0) return
         
         eos_calls = 0
         x = initial_guess
            
         if (doing_which /= doing_get_T) then
            Pgas = arg_not_provided
            T = exp10(the_other_log)
         else
            T = arg_not_provided
            Pgas = exp10(the_other_log)
         end if

         lipar = eos_lipar + species + num_chem_isos
         lrpar = eos_lrpar + lrextras + nv*3 + species
        
         allocate(rpar(lrpar),ipar(lipar),stat=ierr)
         if (ierr /= 0) then
            write(*, *) 'allocate ierr', ierr
            return
         end if

         ipar(i_doing_which) = doing_which
         ipar(i_which_other) = which_other
         ipar(i_handle) = handle
         ipar(i_count) = eos_calls
         ipar(i_species) = species
         i = eos_lipar
         do j=1,species 
            ipar(i+j) = chem_id(j)
         end do
         i = i+species
         do j=1,num_chem_isos 
            ipar(i+j) = net_iso(j)
         end do
         i = i+num_chem_isos

         rpar(r_other_value) = other_value
         rpar(r_Z) = Z
         rpar(r_X) = XH1
         rpar(r_abar) = abar
         rpar(r_zbar) = zbar
         rpar(r_Pgas) = Pgas
         rpar(r_T) = T
         rpar(r_the_other_log) = the_other_log
         i = eos_lrpar
         i = i+nv ! res
         i = i+nv ! d_dlnRho_c_T
         i = i+nv ! d_dlnT_c_Rho
         i = i+1 ! Rho
         i = i+1 ! logRho
         i = i+1 ! dlnRho_dlnPgas_c_T
         i = i+1 ! dlnRho_dlnT_c_Pgas
         rpar(i+1:i+species) = xa(1:species); i = i+species
         if (i /= lrpar) then
            write(*,3) 'i /= lrpar', i, lrpar
            stop 'bad value for lrpar in do_safe_get_Pgas_T'
         end if
         
         
         xb1 = xbnd1; xb3 = xbnd2
         if (xb1 == arg_not_provided .or. xb3 == arg_not_provided .or. xb1 == xb3) then
         
            if (dbg) then
               write(*,*)
               write(*,*) 'call look_for_brackets'
               write(*,2) 'ntry', ntry
               write(*,1) 'x', x
               write(*,1) 'dx', dx
               write(*,1) 'Z', Z
               write(*,1) 'XH1', XH1
               write(*,1) 'abar', abar
               write(*,1) 'zbar', zbar
               write(*,1) 'Pgas', Pgas
               write(*,1) 'T', T
               write(*,1) 'the_other_log', the_other_log
               write(*,*)
            end if

            call look_for_brackets(x, dx, xb1, xb3, get_f_df, y1, y3, &
                     ntry, lrpar, rpar, lipar, ipar, ierr)
            if (ierr /= 0) then
               if (dbg) then
                  write(*, *) 'look_for_brackets returned ierr', ierr
                  write(*,1) 'x', x
                  write(*,1) 'dx', dx
                  write(*,1) 'xb1', xb1
                  write(*,1) 'xb3', xb3
                  write(*,*) 'ntry', ntry
                  write(*,*) 'lrpar', lrpar
                  write(*,*) 'lipar', lipar
               end if
               call dealloc
               return
            end if
            !write(*,*) 'done look_for_brackets'
         else
            if (other_at_bnd1 == arg_not_provided) then
               y1 = get_f_df(xb1, dfdx, lrpar, rpar, lipar, ipar, ierr)
               if (ierr /= 0) then
                  call dealloc
                  return
               end if
            else
               y1 = other_at_bnd1 - other_value
            end if
            if (other_at_bnd2 == arg_not_provided) then
               y3 = get_f_df(xb3, dfdx, lrpar, rpar, lipar, ipar, ierr)
               if (ierr /= 0) then
                  call dealloc
                  return
               end if
            else
               y3 = other_at_bnd2 - other_value
            end if
         end if
         
         if (dbg) then
            write(*,*)
            write(*,*) 'call brent_safe_zero'
            write(*,1) 'xb1', xb1
            write(*,1) 'xb3', xb3
            write(*,1) 'y1', y1
            write(*,1) 'y3', y3
         end if

         x = brent_safe_zero( &
            xb1, xb3, 1d-14, 0.5d0*xacc, 0.5d0*yacc, get_f_df, y1, y3, &
            lrpar, rpar, lipar, ipar, ierr)
         if (ierr /= 0) then
            call dealloc
            return
         end if
      
         i = eos_lrpar
         res = rpar(i+1:i+nv); i = i+nv
         d_dlnRho_c_T = rpar(i+1:i+nv); i = i+nv
         d_dlnT_c_Rho = rpar(i+1:i+nv); i = i+nv
         Rho = rpar(i+1); i = i+1
         logRho = rpar(i+1); i = i+1
         dlnRho_dlnPgas_c_T = rpar(i+1); i = i+1
         dlnRho_dlnT_c_Pgas = rpar(i+1); i = i+1
         i = i + species ! xa
         if (i /= lrpar) then
            write(*,3) 'i /= lrpar', i, lrpar
            stop 'bad value for lrpar at end of do_safe_get_Pgas_T'
         end if
         
         eos_calls = ipar(4)
         
         call dealloc
         
         !write(*,*) 'do_safe_get_Pgas_T eos_calls', eos_calls
         
         contains
         
         subroutine dealloc
            deallocate(rpar,ipar)
         end subroutine dealloc
         
      end subroutine do_safe_get_Pgas_T


      real(dp) function get_f_df(x, dfdx, lrpar, rpar, lipar, ipar, ierr)
         use eos_def, only:EoS_General_Info, get_eos_ptr
         use chem_def, only: num_chem_isos
         integer, intent(in) :: lrpar, lipar
         real(dp), intent(in) :: x
         real(dp), intent(out) :: dfdx
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr

         real(dp) :: new, logT, logPgas
         real(dp) :: Z, XH1, abar, zbar, Pgas, T, other, the_other_log
         
         real(dp), dimension(nv) :: d_res_d_abar, d_res_d_zbar
         
         integer :: i, which_other, handle, doing_which, species, rpar_irho
         real(dp), dimension(:), pointer :: &
            res, d_dlnRho_c_T, d_dlnT_c_Rho, xa
         integer, dimension(:), pointer :: chem_id, net_iso
         real(dp), pointer :: Rho, logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas
         real(dp) :: dfdx_alt
         logical, parameter :: basic_flag = .false.
         type (EoS_General_Info), pointer :: rq

         include 'formats'
         
         ierr = 0
         get_f_df = 0
         
         doing_which = ipar(i_doing_which)
         which_other = ipar(i_which_other)
         handle = ipar(i_handle)
         species = ipar(i_species)
         i = eos_lipar
         chem_id => ipar(i+1:i+species); i = i+species
         net_iso => ipar(i+1:i+num_chem_isos); i = i+num_chem_isos
         if (i /= lipar) then
         end if
         
         call get_eos_ptr(handle, rq, ierr)
         if (ierr /= 0) then
            write(*, *) 'get_eos_ptr returned ierr', ierr
            return
         end if
         dfdx = 0

         other = rpar(r_other_value)
         Z = rpar(r_Z)
         XH1 = rpar(r_X)
         abar = rpar(r_abar)
         zbar = rpar(r_zbar)
         Pgas = rpar(r_Pgas)
         T = rpar(r_T)
         the_other_log = rpar(r_the_other_log)
         
         i = eos_lrpar
         res => rpar(i+1:i+nv); i = i+nv
         d_dlnRho_c_T => rpar(i+1:i+nv); i = i+nv
         d_dlnT_c_Rho => rpar(i+1:i+nv); i = i+nv
         rpar_irho = i+1
         Rho => rpar(i+1); i = i+1
         logRho => rpar(i+1); i = i+1
         dlnRho_dlnPgas_c_T => rpar(i+1); i = i+1
         dlnRho_dlnT_c_Pgas => rpar(i+1); i = i+1
         xa => rpar(i+1:i+species); i = i+species
         if (i /= lrpar) stop 'bad value for lrpar in eosPT get_f_df'
         
         if (doing_which /= doing_get_T) then
            logPgas = x
            Pgas = exp10(logPgas)
            logT = the_other_log
            T = arg_not_provided
         else
            logT = x
            T = exp10(logT)
            logPgas = the_other_log
            Pgas = arg_not_provided
         end if
         
         ierr = 0
         call Get_eosPT_Results(rq, &
               Z, XH1, abar, zbar, &
               species, chem_id, net_iso, xa, &
               Pgas, logPgas, T, logT, &
               Rho, logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas, &
               res, d_dlnRho_c_T, d_dlnT_c_Rho, &
               ierr)
         if (ierr /= 0) then
 22          format(a30, e26.16)
            if (.true.) then
               write(*, *) 'Get_eosPT_Results returned ierr', ierr
               write(*, 22) 'Z', Z
               write(*, 22) 'XH1', XH1
               write(*, 22) 'abar', abar
               write(*, 22) 'zbar', zbar
               write(*, 22) 'Pgas', Pgas
               write(*, 22) 'logPgas', logPgas
               write(*, 22) 'T', T
               write(*, 22) 'logT', logT
               write(*,*)
            end if
            return
         end if
         
         ipar(4) = ipar(4)+1 ! count eos calls
         
         if (doing_which == doing_get_Pgas_for_Rho) then
            new = logRho
         else
            new = res(which_other)
         end if
         get_f_df = new - other
         
         ! f = f(lnRho(lnPgas,lnT),lnT)
         if (doing_which == doing_get_T) then
            dfdx = (d_dlnT_c_Rho(which_other) &
                  + dlnRho_dlnT_c_Pgas*d_dlnRho_c_T(which_other))*ln10
         else if (doing_which == doing_get_Pgas) then
            dfdx = dlnRho_dlnPgas_c_T*d_dlnRho_c_T(which_other)*ln10
         else if (doing_which == doing_get_Pgas_for_Rho) then
            dfdx = dlnRho_dlnPgas_c_T
         else
            stop 'bad value for doing_which in eosPT_eval'
         end if
         
         if (.false. .and. abs(other - 3.5034294596213336d+01) < 1d-14) then
         !if (.true.) then
            if (doing_which /= doing_get_T) then
               write(*,2) 'logPgas, f, dfdx, f/dfdx', ipar(4), logPgas, get_f_df, dfdx, get_f_df/dfdx
            else
               write(*,2) 'logT, f, dfdx, f/dfdx', ipar(4), logT, get_f_df, dfdx, get_f_df/dfdx
            end if
            !write(*,1) 'new', new
            !write(*,1) 'other', other
            !write(*,1) 'get_f_df', get_f_df
            !write(*,1) 'dfdx', dfdx
            !write(*,*)
            !if (ipar(4) > 25) stop
         end if
         
         !write(*,2) 'get_f_df Rho', rpar_irho, rpar(rpar_irho), Rho
         
      end function get_f_df


      end module eosPT_eval
      
