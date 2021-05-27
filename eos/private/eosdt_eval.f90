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
!
! ***********************************************************************

      module eosDT_eval
      use eos_def
      use const_def, only: avo, crad, ln10, arg_not_provided, mp, kerg, dp, qp
      use utils_lib, only: is_bad, mesa_error
      use math_lib
      use eosDT_load_tables, only: load_single_eosDT_table_by_id
      use eos_HELM_eval
      use eoscms_eval, only: Get_CMS_alfa, get_CMS_for_eosdt
      use skye, only: get_Skye_for_eosdt, get_Skye_alfa, get_Skye_alfa_simple


      implicit none

      logical, parameter :: return_ierr_beyond_table_bounds = .true.

      integer, parameter :: i_doing_Rho = 1
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
      integer, parameter :: r_rho = 6
      integer, parameter :: r_T = 7
      integer, parameter :: r_the_other_log = 8
      integer, parameter :: eos_lrpar = 8      

      integer, parameter :: use_none = 1
      integer, parameter :: use_all = 2
      integer, parameter :: blend_in_x = 3
      integer, parameter :: blend_in_y = 4
      integer, parameter :: blend_corner_out = 5
      integer, parameter :: blend_corner_in = 6
      integer, parameter :: blend_diagonal = 7
      integer, parameter :: blend_in_Z = 8
      
      abstract interface
         subroutine get_values_for_eosdt_interface( &
               handle, dbg, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               rho, logRho, T, logT, remaining_fraction, &
               res, d_dlnd, d_dlnT, d_dxa, &
               skip, ierr)
            use const_def, only: dp
            use eos_def, only: nv         
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
         end subroutine get_values_for_eosdt_interface
      end interface
      
         
      contains


      subroutine Test_one_eosDT_component(rq, which_eos, &
            Z, X, abar, zbar, &
            species, chem_id, net_iso, xa, &
            arho, alogrho, atemp, alogtemp, &
            res, d_dlnd, d_dlnT, d_dxa, ierr)
         type (EoS_General_Info), pointer :: rq
         real(dp), intent(in) :: Z, X, abar, zbar
         integer, intent(in) :: which_eos, species
         integer, pointer :: chem_id(:), net_iso(:)
         real(dp), intent(in) :: xa(:)
         real(dp), intent(in) :: arho, alogrho, atemp, alogtemp
         real(dp), intent(inout), dimension(nv) :: res, d_dlnd, d_dlnT
         real(dp), intent(inout), dimension(nv, species) :: d_dxa
         integer, intent(out) :: ierr
         
         real(dp) :: rho, logRho, T, logT
         logical :: skip

         include 'formats'
         
         T = atemp; logT = alogtemp
         if (atemp == arg_not_provided .and. alogtemp == arg_not_provided) then
            ierr = -1; return
         end if
         if (alogtemp == arg_not_provided) logT = log10(T)
         if (atemp == arg_not_provided) T = exp10(logT)
         
         if (T <= 0) then
            ierr = -1
            return
         end if
         
         Rho = arho; logrho = alogrho
         if (arho == arg_not_provided .and. alogrho == arg_not_provided) then
            ierr = -1; return
         end if
         if (alogrho == arg_not_provided) logRho = log10(Rho)
         if (arho == arg_not_provided) Rho = exp10(logRho)
         
         if (Rho <= 0) then
            ierr = -1
            return
         end if
         if (is_bad(Rho) .or. is_bad(T)) then
            ierr = -1
            return
         end if

         select case(which_eos)
         case(i_eos_HELM)
            call get_helm_for_eosdt( &
               rq% handle, dbg, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               rho, logRho, T, logT, 1d0, &
               res, d_dlnd, d_dlnT, d_dxa, &
               skip, ierr)
         case(i_eos_OPAL_SCVH)
            call get_opal_scvh_for_eosdt( &
               rq% handle, dbg, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               rho, logRho, T, logT, 1d0, &
               res, d_dlnd, d_dlnT, d_dxa, &
               skip, ierr)
         case(i_eos_FreeEOS)
            call get_FreeEOS_for_eosdt( &
               rq% handle, dbg, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               rho, logRho, T, logT, 1d0, &
               res, d_dlnd, d_dlnT, d_dxa, &
               skip, ierr)
         case(i_eos_PC)
            call get_PC_for_eosdt( &
               rq% handle, dbg, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               rho, logRho, T, logT, 1d0, &
               res, d_dlnd, d_dlnT, d_dxa, &
               skip, ierr)
         case(i_eos_Skye)
            call get_Skye_for_eosdt( &
               rq% handle, dbg, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               rho, logRho, T, logT, 1d0, &
               res, d_dlnd, d_dlnT, d_dxa, &
               skip, ierr)
         case(i_eos_CMS)
            call get_CMS_for_eosdt( &
               rq% handle, dbg, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               rho, logRho, T, logT, 1d0, &
               res, d_dlnd, d_dlnT, d_dxa, &
               skip, ierr)
         case default
            ierr = -1
         end select
         
         if (ierr /= 0) then
            write(*,*) 'failed in Test_one_eosDT_component', which_eos
            return
         end if
         
         if (skip) then
            write(*,*) 'skipped - no results Test_one_eosDT_component', which_eos
            return
         end if
         
      end subroutine Test_one_eosDT_component


      subroutine Get_eosDT_Results(rq, &
            Z, X, abar, zbar, &
            species, chem_id, net_iso, xa, &
            arho, alogrho, atemp, alogtemp, &
            res, d_dlnd, d_dlnT, d_dxa, ierr)
         type (EoS_General_Info), pointer :: rq
         real(dp), intent(in) :: Z, X, abar, zbar          
         integer, intent(in) :: species
         integer, pointer :: chem_id(:), net_iso(:)
         real(dp), intent(in) :: xa(:)
         real(dp), intent(in) :: arho, alogrho, atemp, alogtemp
         real(dp), intent(inout), dimension(nv) :: res, d_dlnd, d_dlnT
         real(dp), intent(inout), dimension(nv, species) :: d_dxa
         integer, intent(out) :: ierr
         
         real(dp) :: rho, logRho, T, logT
         logical :: skip, dbg

         include 'formats'
         

         T = atemp; logT = alogtemp
         if (atemp == arg_not_provided .and. alogtemp == arg_not_provided) then
            ierr = -1; return
         end if
         if (alogtemp == arg_not_provided) logT = log10(T)
         if (atemp == arg_not_provided) T = exp10(logT)
         
         if (T <= 0) then
            ierr = -1
            return
         end if
         
         Rho = arho; logrho = alogrho
         if (arho == arg_not_provided .and. alogrho == arg_not_provided) then
            ierr = -1; return
         end if
         if (alogrho == arg_not_provided) logRho = log10(Rho)
         if (arho == arg_not_provided) Rho = exp10(logRho)
         
         if (Rho <= 0) then
            ierr = -1
            return
         end if
         if (is_bad(Rho) .or. is_bad(T)) then
            ierr = -1
            return
         end if

         dbg = rq% dbg
         if (dbg) dbg = & ! check limits
            logT >= rq% logT_lo .and. logT <= rq% logT_hi .and. &
            logRho >= rq% logRho_lo .and. logRho <= rq% logRho_hi .and. &
            X >= rq% X_lo .and. X <= rq% X_hi .and. &
            Z >= rq% Z_lo .and. Z <= rq% Z_hi
         
         call get_level1_for_eosdt( &
            rq% handle, dbg, Z, X, abar, zbar, &
            species, chem_id, net_iso, xa, &
            rho, logRho, T, logT, 1d0, &
            res, d_dlnd, d_dlnT, d_dxa, &
            skip, ierr)
         if (skip) ierr = -1
         if (ierr /= 0) return
         
         if (eos_test_partials) then   
            eos_test_partials_val = abar
            eos_test_partials_dval_dx = 0
            write(*,*) 'eos_test_partials'
         end if
         
      end subroutine Get_eosDT_Results
      
      
      subroutine get_level1_for_eosdt( & ! CMS
            handle, dbg, Z, X, abar, zbar, &
            species, chem_id, net_iso, xa, &
            rho, logRho, T, logT, remaining_fraction, &
            res, d_dlnd, d_dlnT, d_dxa, &
            skip, ierr)
         use eoscms_eval, only: Get_CMS_alfa, get_CMS_for_eosdt
         integer, intent(in) :: handle
         logical, intent(in) :: dbg
         real(dp), intent(in) :: Z, X, abar, zbar, remaining_fraction
         integer, intent(in) :: species
         integer, pointer :: chem_id(:), net_iso(:)
         real(dp), intent(in) :: xa(:)
         real(dp), intent(in) :: rho, logRho, T, logT
         real(dp), intent(inout), dimension(nv) :: &
            res, d_dlnd, d_dlnT
         real(dp), intent(inout), dimension(nv, species) :: d_dxa
         logical, intent(out) :: skip
         integer, intent(out) :: ierr
         
         real(dp) :: alfa, d_alfa_dlogT, d_alfa_dlogRho
         type (EoS_General_Info), pointer :: rq
         procedure (get_values_for_eosdt_interface), pointer :: get_1st, get_2nd

         include 'formats'
         
         ierr = 0
         rq => eos_handles(handle)

         if (rq% use_CMS) then
            call Get_CMS_alfa( &
               rq, logRho, logT, Z, abar, zbar, &
               alfa, d_alfa_dlogT, d_alfa_dlogRho, &
               ierr)
            if (ierr /= 0) return
         else
            alfa = 1d0 ! no CMS
            d_alfa_dlogT = 0d0
            d_alfa_dlogRho = 0d0
         end if
         
         if (dbg) write(*,1) 'CMS', (1d0 - alfa)*remaining_fraction
         
         get_1st => get_CMS_for_eosdt
         get_2nd => get_level2_for_eosdt
         call combine_for_eosdt( &
            get_1st, get_2nd, alfa*remaining_fraction, &
            alfa, d_alfa_dlogT, d_alfa_dlogRho, &
            rq, dbg, Z, X, abar, zbar, &
            species, chem_id, net_iso, xa, &
            rho, logRho, T, logT, &
            res, d_dlnd, d_dlnT, d_dxa, &
            skip, ierr)
         if (ierr /= 0 .and. rq% okay_to_convert_ierr_to_skip) then
            skip = .true.
            ierr = 0
         end if
            
      end subroutine get_level1_for_eosdt
      
      
      subroutine get_level2_for_eosdt( & ! Skye
            handle, dbg, Z, X, abar, zbar, &
            species, chem_id, net_iso, xa, &
            rho, logRho, T, logT, remaining_fraction, &
            res, d_dlnd, d_dlnT, d_dxa, &
            skip, ierr)
         integer, intent(in) :: handle
         logical, intent(in) :: dbg
         real(dp), intent(in) :: Z, X, abar, zbar, remaining_fraction
         integer, intent(in) :: species
         integer, pointer :: chem_id(:), net_iso(:)
         real(dp), intent(in) :: xa(:)
         real(dp), intent(in) :: rho, logRho, T, logT
         real(dp), intent(inout), dimension(nv) :: &
            res, d_dlnd, d_dlnT
         real(dp), intent(inout), dimension(nv, species) :: d_dxa
         logical, intent(out) :: skip
         integer, intent(out) :: ierr
         
         real(dp) :: alfa, d_alfa_dlogT, d_alfa_dlogRho
         type (EoS_General_Info), pointer :: rq
         procedure (get_values_for_eosdt_interface), pointer :: get_1st, get_2nd

         include 'formats'
         
         ierr = 0
         rq => eos_handles(handle)

         if (rq% use_Skye) then
            if (rq% use_simple_Skye_blends) then
               call Get_Skye_alfa_simple( &
                  rq, logRho, logT, Z, abar, zbar, &
                  alfa, d_alfa_dlogT, d_alfa_dlogRho, &
                  ierr)
            else
               call Get_Skye_alfa( &
                  rq, logRho, logT, Z, abar, zbar, &
                  alfa, d_alfa_dlogT, d_alfa_dlogRho, &
                  ierr)
            end if
            if (ierr /= 0) return               
         else
            alfa = 1d0 ! no Skye
            d_alfa_dlogT = 0d0
            d_alfa_dlogRho = 0d0
         end if 
            
         if (dbg) write(*,1) 'Skye', (1d0 - alfa)*remaining_fraction
         get_1st => get_Skye_for_eosdt

         get_2nd => get_level3_for_eosdt
         call combine_for_eosdt( &
            get_1st, get_2nd, alfa*remaining_fraction, &
            alfa, d_alfa_dlogT, d_alfa_dlogRho, &
            rq, dbg, Z, X, abar, zbar, &
            species, chem_id, net_iso, xa, &
            rho, logRho, T, logT, &
            res, d_dlnd, d_dlnT, d_dxa, &
            skip, ierr)
         if (ierr /= 0 .and. rq% okay_to_convert_ierr_to_skip) then
            skip = .true.
            ierr = 0
         end if
            
      end subroutine get_level2_for_eosdt
      
      
      subroutine get_level3_for_eosdt( & ! PC
            handle, dbg, Z, X, abar, zbar, &
            species, chem_id, net_iso, xa, &
            rho, logRho, T, logT, remaining_fraction, &
            res, d_dlnd, d_dlnT, d_dxa, &
            skip, ierr)
         use eospc_eval, only: Get_PC_alfa, get_PC_for_eosdt
         integer, intent(in) :: handle
         logical, intent(in) :: dbg
         real(dp), intent(in) :: Z, X, abar, zbar, remaining_fraction
         integer, intent(in) :: species
         integer, pointer :: chem_id(:), net_iso(:)
         real(dp), intent(in) :: xa(:)
         real(dp), intent(in) :: rho, logRho, T, logT
         real(dp), intent(inout), dimension(nv) :: &
            res, d_dlnd, d_dlnT
         real(dp), intent(inout), dimension(nv, species) :: d_dxa
         logical, intent(out) :: skip
         integer, intent(out) :: ierr
         
         real(dp) :: alfa, d_alfa_dlogT, d_alfa_dlogRho
         type (EoS_General_Info), pointer :: rq
         procedure (get_values_for_eosdt_interface), pointer :: get_1st, get_2nd

         include 'formats'
         
         ierr = 0
         rq => eos_handles(handle)

         if (rq% use_PC) then
            call Get_PC_alfa( & 
               rq, logRho, logT, Z, abar, zbar, &
               alfa, d_alfa_dlogT, d_alfa_dlogRho, &
               ierr)
            if (ierr /= 0) return
         else
            alfa = 1d0 ! no PC
            d_alfa_dlogT = 0d0
            d_alfa_dlogRho = 0d0
         end if 
            
         if (dbg) write(*,1) 'PC', (1d0 - alfa)*remaining_fraction
         get_1st => get_PC_for_eosdt

         get_2nd => get_level4_for_eosdt
         call combine_for_eosdt( &
            get_1st, get_2nd, alfa*remaining_fraction, &
            alfa, d_alfa_dlogT, d_alfa_dlogRho, &
            rq, dbg, Z, X, abar, zbar, &
            species, chem_id, net_iso, xa, &
            rho, logRho, T, logT, &
            res, d_dlnd, d_dlnT, d_dxa, &
            skip, ierr)
         if (ierr /= 0 .and. rq% okay_to_convert_ierr_to_skip) then
            skip = .true.
            ierr = 0
         end if
            
      end subroutine get_level3_for_eosdt
      
      
      subroutine get_level4_for_eosdt( & ! FreeEOS
            handle, dbg, Z, X, abar, zbar, &
            species, chem_id, net_iso, xa, &
            rho, logRho, T, logT, remaining_fraction, &
            res, d_dlnd, d_dlnT, d_dxa, &
            skip, ierr)
         use skye, only: get_Skye_for_eosdt, get_Skye_alfa
         integer, intent(in) :: handle
         logical, intent(in) :: dbg
         real(dp), intent(in) :: Z, X, abar, zbar, remaining_fraction
         integer, intent(in) :: species
         integer, pointer :: chem_id(:), net_iso(:)
         real(dp), intent(in) :: xa(:)
         real(dp), intent(in) :: rho, logRho, T, logT
         real(dp), intent(inout), dimension(nv) :: &
            res, d_dlnd, d_dlnT
         real(dp), intent(inout), dimension(nv, species) :: d_dxa
         logical, intent(out) :: skip
         integer, intent(out) :: ierr
         
         real(dp) :: alfa, d_alfa_dlogT, d_alfa_dlogRho
         type (EoS_General_Info), pointer :: rq
         procedure (get_values_for_eosdt_interface), pointer :: get_1st, get_2nd

         include 'formats'
         
         ierr = 0
         rq => eos_handles(handle)

         if (rq% use_FreeEOS) then
            call Get_FreeEOS_alfa( & 
               rq, dbg, logRho, logT, Z, abar, zbar, &
               alfa, d_alfa_dlogT, d_alfa_dlogRho, &
               ierr)
            if (ierr /= 0) return
         else
            alfa = 1d0 ! no FreeEOS
            d_alfa_dlogT = 0d0
            d_alfa_dlogRho = 0d0
         end if

         if (dbg) write(*,1) 'FreeEOS', (1d0 - alfa)*remaining_fraction
         get_1st => get_FreeEOS_for_eosdt
         
         get_2nd => get_level5_for_eosdt
         call combine_for_eosdt( &
            get_1st, get_2nd, alfa*remaining_fraction, &
            alfa, d_alfa_dlogT, d_alfa_dlogRho, &
            rq, dbg, Z, X, abar, zbar, &
            species, chem_id, net_iso, xa, &
            rho, logRho, T, logT, &
            res, d_dlnd, d_dlnT, d_dxa, &
            skip, ierr)
         if (ierr /= 0 .and. rq% okay_to_convert_ierr_to_skip) then
            skip = .true.
            ierr = 0
         end if
            
      end subroutine get_level4_for_eosdt
      
      
      subroutine get_level5_for_eosdt( &  ! OPAL/SCVH
            handle, dbg, Z, X, abar, zbar, &
            species, chem_id, net_iso, xa, &
            rho, logRho, T_in, logT_in, remaining_fraction, &
            res, d_dlnd, d_dlnT, d_dxa, &
            skip, ierr)
         integer, intent(in) :: handle
         logical, intent(in) :: dbg
         real(dp), intent(in) :: &
            Z, X, abar, zbar, remaining_fraction
         integer, intent(in) :: species
         integer, pointer :: chem_id(:), net_iso(:)
         real(dp), intent(in) :: xa(:)
         real(dp), intent(in) :: rho, logRho, T_in, logT_in
         real(dp), intent(inout), dimension(nv) :: &
            res, d_dlnd, d_dlnT
         real(dp), intent(inout), dimension(nv, species) :: d_dxa
         logical, intent(out) :: skip
         integer, intent(out) :: ierr
         
         real(dp) :: alfa, d_alfa_dlogT, d_alfa_dlogRho, &
            logT_HELM, T_HELM, logQ, logQ2, T, logT
         type (EoS_General_Info), pointer :: rq
         procedure (get_values_for_eosdt_interface), pointer :: get_1st, get_2nd

         include 'formats'
         
         ierr = 0
         rq => eos_handles(handle)
         
         T = T_in
         logT = logT_in

         if (rq% use_OPAL_SCVH) then 
            call get_opal_scvh_alfa_and_partials( &
               rq, logT, logRho, Z, &
               alfa, d_alfa_dlogRho, d_alfa_dlogT, ierr)
            if (ierr /= 0) return
         else
            alfa = 1d0 ! no OPAL_SCVH
            d_alfa_dlogT = 0d0
            d_alfa_dlogRho = 0d0
         end if
         
         if (dbg) write(*,1) 'OPAL/SCVH', (1d0 - alfa)*remaining_fraction
         if (dbg) write(*,1) 'HELM', alfa*remaining_fraction
         
         get_1st => get_opal_scvh_for_eosdt
         get_2nd => get_helm_for_eosdt
         call combine_for_eosdt( &
            get_1st, get_2nd, alfa*remaining_fraction, &
            alfa, d_alfa_dlogT, d_alfa_dlogRho, &
            rq, dbg, Z, X, abar, zbar, &
            species, chem_id, net_iso, xa, &
            rho, logRho, T, logT, &
            res, d_dlnd, d_dlnT, d_dxa, &
            skip, ierr)
         if (ierr /= 0 .and. rq% okay_to_convert_ierr_to_skip) then
            skip = .true.
            ierr = 0
         end if
            
      end subroutine get_level5_for_eosdt


      subroutine Get_FreeEOS_alfa( &
            rq, dbg, logRho, logT, Z, abar, zbar, &
            alfa, d_alfa_dlogT, d_alfa_dlogRho, &
            ierr)
         use const_def
         use auto_diff
         type (EoS_General_Info), pointer :: rq
         logical, intent(in) :: dbg
         real(dp), intent(in) :: logRho, logT, Z, abar, zbar
         real(dp), intent(out) :: alfa, d_alfa_dlogT, d_alfa_dlogRho
         integer, intent(out) :: ierr
         real(dp) :: logQ_cut_lo, logQ_cut_hi
         type(auto_diff_real_2var_order1) :: logT_auto, logRho_auto, logQ_auto
         type(auto_diff_real_2var_order1) :: blend, blend_logT, blend_logRho, blend_logQ
         type(auto_diff_real_2var_order1) :: blend_cut, blend_logT_cut, blend_logRho_cut, blend_logQ_cut

         include 'formats'

         ierr = 0

         ! auto diff
         ! var1: logRho
         ! var2: logT

         logRho_auto% val = logRho
         logRho_auto% d1val1 = 1
         logRho_auto% d1val2 = 0

         logT_auto% val = logT
         logT_auto% d1val1 = 0
         logT_auto% d1val2 = 1

         logQ_auto = logRho_auto - 2d0*logT_auto + 12d0

         ! The FreeEOS region usually looks like
         !
         !                 ________C________
         !                /                 |
         !               /                  D
         !              /              __E__|
         !             B              /
         !            /              /
         !           /              F
         !          /              /
         !         /              /
         !        /              /
         !       /              |
         !      /               G
         !     /________A_______|
         !
         ! where blend A and C come from blend_logT,
         ! blend B comes from blend_logQ,
         ! blend D usually comes from PC/Skye,
         ! blend E comes from blend_logT_cut,
         ! blend F comes from blend_logQ_cut,
         ! blend G comes from blend_logRho_cut.


         ! logT blend
         if (logT_auto < rq% logT_min_FreeEOS_lo) then
            blend_logT = 0d0
         else if (logT_auto < rq% logT_min_FreeEOS_hi) then
            blend_logT = (logT_auto - rQ% logT_min_FreeEOS_lo) / (rq% logT_min_FreeEOS_hi - rq% logT_min_FreeEOS_lo)
         else if (logT_auto < rq% logT_max_FreeEOS_lo) then
            blend_logT = 1d0
         else if (logT_auto < rq% logT_max_FreeEOS_hi) then
            blend_logT = (logT_auto - rQ% logT_max_FreeEOS_hi) / (rq% logT_max_FreeEOS_lo - rq% logT_max_FreeEOS_hi)
         else
            blend_logT = 0
         end if


         ! logRho blend
         if (logRho_auto < rq% logRho_min_FreeEOS_lo) then
            blend_logRho = 0d0
         else if (logRho_auto < rq% logRho_min_FreeEOS_lo) then
            blend_logRho = (logRho_auto - rQ% logRho_min_FreeEOS_lo) / (rq% logRho_min_FreeEOS_hi - rq% logRho_min_FreeEOS_lo)
         else if (logRho_auto < rq% logRho_max_FreeEOS_lo) then
            blend_logRho = 1d0
         else if (logRho_auto < rq% logRho_max_FreeEOS_hi) then
            blend_logRho = (logRho_auto - rQ% logRho_max_FreeEOS_hi) / (rq% logRho_max_FreeEOS_lo - rq% logRho_max_FreeEOS_hi)
         else
            blend_logRho = 0
         end if


         ! logQ blend
         if (logQ_auto < rq% logQ_min_FreeEOS_lo) then
            blend_logQ = 0d0
         else if (logQ_auto < rq% logQ_min_FreeEOS_hi) then
            blend_logQ = (logQ_auto - rQ% logQ_min_FreeEOS_lo) / (rq% logQ_min_FreeEOS_hi - rq% logQ_min_FreeEOS_lo)
         else if (logQ_auto < rq% logQ_max_FreeEOS_lo) then
            blend_logQ = 1d0
         else if (logQ_auto < rq% logQ_max_FreeEOS_hi) then
            blend_logQ = (logQ_auto - rQ% logQ_max_FreeEOS_hi) / (rq% logQ_max_FreeEOS_lo - rq% logQ_max_FreeEOS_hi)
         else
            blend_logQ = 0
         end if


         ! cut blend logRho
         if (logRho_auto < rq% logRho_cut_FreeEOS_lo) then
            blend_logRho_cut = 1d0
         else if (logRho_auto < rq% logRho_cut_FreeEOS_hi) then
            blend_logRho_cut = (logRho_auto - rQ% logRho_cut_FreeEOS_hi) / (rq% logRho_cut_FreeEOS_lo - rq% logRho_cut_FreeEOS_hi)
         else
            blend_logRho_cut = 0d0
         end if

         ! cut blend logT
         if (logT_auto < rq% logT_cut_FreeEOS_lo) then
            blend_logT_cut = 0d0
         else if (logT_auto < rq% logT_cut_FreeEOS_hi) then
            blend_logT_cut = (logT_auto - rQ% logT_cut_FreeEOS_lo) / (rq% logT_cut_FreeEOS_hi - rq% logT_cut_FreeEOS_lo)
         else
            blend_logT_cut = 1d0
         end if

         ! cut blend logQ
         if (Z <= rq% logQ_cut_FreeEOS_lo_Z_max) then
            logQ_cut_hi = rq% logQ_cut_lo_Z_FreeEOS_hi
            logQ_cut_lo = rq% logQ_cut_lo_Z_FreeEOS_lo
         else
            logQ_cut_hi = rq% logQ_cut_hi_Z_FreeEOS_hi
            logQ_cut_lo = rq% logQ_cut_hi_Z_FreeEOS_lo
         end if

         if (logQ_auto < logQ_cut_lo) then
            blend_logQ_cut = 1d0
         else if (logQ_auto < logQ_cut_hi) then
            blend_logQ_cut = (logQ_auto - logQ_cut_hi) / (logQ_cut_lo - logQ_cut_hi)
         else
            blend_logQ_cut = 0d0
         end if

         ! combine cut blends
         blend_cut = 0
         if (blend_logT_cut == 0) then
            if (blend_logQ_cut == 0) then
               blend_cut = blend_logRho_cut
            else if (blend_logQ_cut == 1) then
               blend_cut = 1d0
            else
               blend_cut = max(blend_logRho_cut, blend_logQ_cut)
            endif
         elseif (blend_logT_cut == 1) then
            blend_cut = 1d0
         else
            if (blend_logQ_cut == 0) then
               blend_cut = blend_logT_cut
            else if (blend_logQ_cut == 1) then
               blend_cut = 1d0
            else
               blend_cut = max(blend_logT_cut, blend_logQ_cut)
            endif
         end if

         ! combine all blends
         blend = blend_cut * blend_logT * blend_logRho * blend_logQ

         ! unpack auto_diff
         alfa = 1d0 - blend% val
         d_alfa_dlogRho = -blend% d1val1
         d_alfa_dlogT = -blend% d1val2

      end subroutine Get_FreeEOS_alfa
      

      subroutine get_opal_scvh_for_eosdt( &
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
         call get1_for_eosdt( &
            handle, eosdt_OPAL_SCVH, dbg, &
            Z, X, abar, zbar, &
            species, chem_id, net_iso, xa, &
            rho, logRho, T, logT, remaining_fraction, &
            res, d_dlnd, d_dlnT, d_dxa, &
            skip, ierr)

         ! zero phase information
         res(i_phase:i_latent_ddlnRho) = 0d0
         d_dlnT(i_phase:i_latent_ddlnRho) = 0d0
         d_dlnd(i_phase:i_latent_ddlnRho) = 0d0

         ! zero all components
         res(i_frac:i_frac+num_eos_frac_results-1) = 0.0
         d_dlnd(i_frac:i_frac+num_eos_frac_results-1) = 0.0
         d_dlnT(i_frac:i_frac+num_eos_frac_results-1) = 0.0

         ! mark this one
         res(i_frac_OPAL_SCVH) = 1.0

      end subroutine get_opal_scvh_for_eosdt      


      subroutine get_FreeEOS_for_eosdt( &
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
         call get1_for_eosdt( &
            handle, eosdt_max_FreeEOS, dbg, Z, X, abar, zbar, &
            species, chem_id, net_iso, xa, &
            rho, logRho, T, logT, remaining_fraction, &
            res, d_dlnd, d_dlnT, d_dxa, &
            skip, ierr)

         ! zero phase information
         res(i_phase:i_latent_ddlnRho) = 0d0
         d_dlnT(i_phase:i_latent_ddlnRho) = 0d0
         d_dlnd(i_phase:i_latent_ddlnRho) = 0d0

         ! zero all components
         res(i_frac:i_frac+num_eos_frac_results-1) = 0.0
         d_dlnd(i_frac:i_frac+num_eos_frac_results-1) = 0.0
         d_dlnT(i_frac:i_frac+num_eos_frac_results-1) = 0.0

         ! mark this one
         res(i_frac_FreeEOS) = 1.0

      end subroutine get_FreeEOS_for_eosdt      


      subroutine get_opal_scvh_alfa_and_partials( &
         rq, logT, logRho, Z, alfa, d_alfa_dlogRho, d_alfa_dlogT, ierr)
         type (EoS_General_Info), pointer :: rq
         real(dp), intent(in) :: logT, logRho, Z
         real(dp), intent(out) :: alfa, d_alfa_dlogRho, d_alfa_dlogT
         integer, intent(out) :: ierr
         
         integer :: iregion
         real(dp) :: logRho1_max, logRho1, logRho2, logRho5, logRho6, logRho7, &
            logRho8, logT5, logT6, logT3, logT4
         real(dp) :: logQ1, logQ2, logQ3, logQ4, logQmax, Z_all_HELM, Z_no_HELM
         real(dp) :: beta, logRho_lo, logRho_hi, &
            logT1, logT2, logT7, logT8, logRho3, logRho4
         real(dp) :: logQ, A, B, dA_dlnT, dA_dlnRho, dB_dlnT, dB_dlnRho
         real(dp) :: c_dx, c_dy, d_dx_dlogT, d_dx_dlogRho, d_dy_dlogT, d_dy_dlogRho
         real(dp), parameter :: tiny = 1d-20
         
         logical :: debug
         
         include 'formats'

         logRho1_max = 3.71d0

         logQ1 = 5.69d0 ! SCVH full off for logQ > this
         logQ2 = 5.68d0 ! must have logQ2 < logQ1
         logQmax = rq% logQ_max_OPAL_SCVH ! 5.3
         logQ3 = rq% logQ_min_OPAL_SCVH ! -8.0
         logQ4 = rq% logQ_min_OPAL_SCVH ! -8.0
      
         logRho5 = rq% logRho_min_OPAL_SCVH_limit ! -14.299
         logRho6 = logRho5 - 1d-3 ! -14.3
         logRho7 = -14.90d0
         logRho8 = -14.99d0

         logRho1 = rq% logRho1_OPAL_SCVH_limit ! 3.50
         if (logRho1 > logRho1_max) then
            write(*,*) 'sorry: value for logRho1_OPAL_SCVH_limit is too large.  max allowed is ', &
               logRho1_max
            ierr = -1
            return
         end if
         logRho2 = rq% logRho2_OPAL_SCVH_limit ! 3.48

         logT8 = rq% logT_low_all_HELM ! 2.2
         logT7 = rq% logT_low_all_SCVH ! 2.3
         logT6 = 4.890d0   
         logT5 = 4.899d0   ! problems with blend here so just jump
         logT2 = rq% logT_all_OPAL ! 7.6
         logT1 = rq% logT_all_HELM ! 7.7
         
         Z_all_HELM = rq% Z_all_HELM
         Z_no_HELM = rq% Z_all_OPAL
         
         if (logT >= logT1) then ! just use other

            alfa = 1d0
            beta = 0d0
         
         else
         
            logT3 = (logRho1 - logQ1 + 12d0)/2d0
            logT4 = (logRho2 - logQ2 + 12d0)/2d0

            logRho3 = logQ1 + 2*logT7 - 12d0
            logRho4 = logQ2 + 2*logT7 - 12d0
            
            if (.false.) then
               write(*,*)
               write(*,1) 'logRho1', logRho1
               write(*,1) 'logRho2', logRho2
               write(*,1) 'logRho3', logRho3
               write(*,1) 'logRho4', logRho4
               write(*,1) 'logRho5', logRho5
               write(*,1) 'logRho6', logRho6
               write(*,1) 'logRho7', logRho7
               write(*,1) 'logRho8', logRho8
               write(*,*)
               write(*,1) 'logT1', logT1
               write(*,1) 'logT2', logT2
               write(*,1) 'logT3', logT3
               write(*,1) 'logT4', logT4
               write(*,1) 'logT5', logT5
               write(*,1) 'logT6', logT6
               write(*,1) 'logT7', logT7
               write(*,1) 'logT8', logT8
               write(*,*)
               write(*,1) 'logQ1', logQ1
               write(*,1) 'logQ2', logQ2
               write(*,*)
               stop 'eosdt_eval'
            end if

            ! check validity of Rho's and T's for region boundaries
            if (logRho1 <= logRho2 .or. logRho2 <= logRho3 .or. &
                logRho3 <= logRho4 .or. logRho4 <= logRho5 .or. &
                logRho5 <= logRho6 .or. logRho6 <= logRho7 .or. &
                logRho7 <= logRho8 .or. &
                logT1 <= logT2 .or. logT2 <= logT3 .or. logT3 <= logT4 .or. &
                logT7 <= logT8) then
               write(*,*)
               write(*,*)
               write(*,*)
               write(*,'(a)') 'must have strictly decreasing values for eos logT + logRho region boundaries'
               write(*,*)
               write(*,1) 'logRho1', logRho1
               write(*,1) 'logRho2', logRho2
               write(*,1) 'logRho3', logRho3
               write(*,1) 'logRho4', logRho4
               write(*,1) 'logRho5', logRho5
               write(*,1) 'logRho6', logRho6
               write(*,1) 'logRho7', logRho7
               write(*,1) 'logRho8', logRho8
               write(*,*)
               write(*,1) 'logT1', logT1
               write(*,1) 'logT2', logT2
               write(*,1) 'logT3', logT3
               write(*,1) 'logT4', logT4
               write(*,1) 'logT5', logT5
               write(*,1) 'logT6', logT6
               write(*,1) 'logT7', logT7
               write(*,1) 'logT8', logT8
               write(*,*)
               write(*,1) 'logQ1', logQ1
               write(*,1) 'logQ2', logQ2
               write(*,*)
               write(*,*)
               if (logT3 <= logT4) then
                  write(*,'(a)') 'must have logRho1 > logRho2 + logQ1 - logQ2'
                  write(*,1) 'must have logRho1 > ', logRho2 + logQ1 - logQ2
                  write(*,1) 'logRho1', logRho1
                  write(*,'(a)') 'logRho1_OPAL_SCVH_limit sets logRho1'
                  write(*,'(a)') 'logRho2_OPAL_SCVH_limit sets logRho2'
                  write(*,1) 'max allowed logRho1 is', logRho1_max
               end if
               write(*,*)
               ierr = -1
               return
            end if
         
            call determine_region_opal_scvh
         
            call set_alfa_and_partials
            if (ierr /= 0) return
            
         end if

         
         contains
         
         
         subroutine determine_region_opal_scvh
            logical, parameter :: dbg = .false.
            real(dp) :: logRho_hi, logRho_lo, d_logRho_dlogT, &
               d_alfa_dlogQ, dlogQ_dlogRho, dlogQ_dlogT, Z_all_HELM
            
            include 'formats'

            logQ = logRho - 2d0*logT + 12d0
            d_dx_dlogRho=0d0
            d_dy_dlogT=0d0

            ! in high-Z fall back to HELM
            if (Z >= rq% Z_all_HELM) then
               iregion = use_none
               if (dbg) then
                  write(*,1) 'iregion = use_none'
                  write(*,1) 'Z Z_all_HELM', Z, rq% Z_all_HELM
                  stop
               end if
               return
            end if
            
            ! blends in T/Rho

            if (logT >= logT1 .or. logT <= logT8 .or. logRho >= logRho1 .or. &
                  logQ <= logQ4 .or. logQ >= logQmax .or. &
                  (logRho <= logRho6 .and. logT <= logT6)) then
               if (dbg) then
                  write(*,*) 'logT >= logT1', logT >= logT1
                  write(*,*) 'logT <= logT8', logT <= logT8
                  write(*,*) 'logRho >= logRho1', logRho >= logRho1
                  write(*,*) 'logQ <= logQ4', logQ <= logQ4, logQ, logQ4
                  write(*,*) 'logQ >= logQmax', logQ >= logQmax
                  write(*,*) 'logRho <= logRho6 .and. logT <= logT6', logRho <= logRho6 .and. logT <= logT6
                  write(*,1) 'iregion = use_none 1 logT logT5 logT6', logT, logT5, logT6
               end if
               iregion = use_none
               
            else if (logQ <= logQ3 .and. logT >= logT5) then ! blend in Q
               d_alfa_dlogQ = 1d0/(logQ4 - logQ3)
               alfa = (logQ - logQ3)*d_alfa_dlogQ
               dlogQ_dlogRho = 1d0
               dlogQ_dlogT = -2d0
               d_dx_dlogRho = d_alfa_dlogQ*dlogQ_dlogRho
               d_dy_dlogT = d_alfa_dlogQ*dlogQ_dlogT
               c_dx = alfa
               if (dbg) then
                  write(*,*) 'iregion = blend_diagonal'
                  write(*,1) 'logQ3', logQ3
                  write(*,1) 'logQ', logQ
                  write(*,1) 'logQ4', logQ4
                  write(*,1) 'd_alfa_dlogQ', d_alfa_dlogQ
                  write(*,1) 'c_dx', c_dx
                  write(*,1) 'd_dx_dlogRho', d_dx_dlogRho
                  write(*,1) 'd_dy_dlogT', d_dy_dlogT
               end if
               iregion = blend_diagonal
               
            else if (logT >= logT2) then         
               if (dbg) write(*,*) 'logT >= logT2', logT, logT2
               if (logT1 - logT2 < 0.01d0) then
                  d_dy_dlogT = 0d0
                  ! bad blend partials cause problems for 150M_z1m4_pre_ms_to_collapse
                  ! have tried to fix, but failed.  hence this awful workaround.
               else
                  d_dy_dlogT = 1/(logT1 - logT2) 
               end if
               c_dy = (logT - logT2)*d_dy_dlogT
               if (logRho > logRho2) then
                  if (dbg) write(*,*) 'logRho > logRho2', logRho, logRho2
                  d_dx_dlogRho = 1/(logRho1 - logRho2)
                  c_dx = (logRho - logRho2)*d_dx_dlogRho 
                  if (dbg) write(*,*) 'iregion = blend_corner_out'
                  iregion = blend_corner_out
               else ! logRho <= logRho2
                  if (dbg) write(*,*) 'logRho <= logRho2', logRho, logRho2
                  if (dbg) write(*,*) 'iregion = blend_in_y'
                  iregion = blend_in_y
               end if        
                
            else if (logT >= logT3) then  ! NOTE: this assumes logT3 > logT4
               if (dbg) write(*,*) 'logT >= logT3', logT, logT3
               if (logRho > logRho2) then
                  if (dbg) write(*,*) 'logRho > logRho2', logRho, logRho2
                  d_dx_dlogRho = 1/(logRho1 - logRho2)
                  c_dx = (logRho - logRho2)*d_dx_dlogRho
                  if (dbg) write(*,*) 'iregion = blend_in_x'
                  iregion = blend_in_x
               else
                  if (dbg) write(*,*) 'logRho <= logRho2', logRho, logRho2
                  if (dbg) write(*,*) 'iregion = use_all'
                  iregion = use_all
               end if    
               
            else if (logT >= logT4) then         
               if (dbg) write(*,*) 'logT >= logT4', logT, logT4
               logRho_hi = logQ1 + 2*logT - 12
               if (logRho >= logRho_hi) then
                  if (dbg) write(*,*) 'logRho >= logRho_hi', logRho, logRho_hi
                  if (dbg) write(*,*) 'iregion = use_none 2'
                  iregion = use_none
               else if (logRho > logRho2) then
                  if (dbg) write(*,*) 'logRho > logRho2', logRho, logRho2
                  d_dx_dlogRho = 1/(logRho_hi - logRho2)
                  c_dx = (logRho - logRho2)*d_dx_dlogRho
                  if (dbg) write(*,*) 'iregion = blend_in_x'
                  iregion = blend_in_x
               else ! logRho <= logRho2
                  if (dbg) write(*,*) 'logRho <= logRho2', logRho, logRho2
                  if (dbg) write(*,*) 'iregion = use_all'
                  iregion = use_all
               end if           
                       
            else if (logRho > logRho4) then         
               if (dbg) write(*,*) 'logRho > logRho4', logRho, logRho4
               if (logT > logT7) then
                  A = ((logQ1+2*logT4-12) - logRho3)/(logT4-logT7)
                  logRho_hi = logRho3 + (logT-logT7)*A
                  B = (logRho2-logRho4)/(logT4-logT7)
                  logRho_lo = logRho4 + (logT-logT7)*B
                  if (logRho >= logRho_hi) then
                     if (dbg) write(*,*) 'logRho >= logRho_hi', logRho, logRho_hi
                     if (dbg) write(*,*) 'iregion = use_none 3'
                     iregion = use_none
                  else if (logRho >= logRho_lo) then
                     if (dbg) write(*,*) 'logRho >= logRho_lo', logRho, logRho_lo
                     c_dx = (logRho - logRho_lo)/(logRho_hi - logRho_lo)
                     d_dx_dlogRho = 1/(logRho3 - logRho4 + (A - B)*(logT - logT7))
                     if (dbg) write(*,*) 'iregion = blend_in_x'
                     iregion = blend_in_x               
                  else ! logRho < logRho_lo
                     if (dbg) write(*,*) 'logRho < logRho_lo', logRho, logRho_lo
                     if (dbg) write(*,*) 'iregion = use_all'
                     iregion = use_all               
                  end if            
               else ! logT is > logT8            
                  if (dbg) write(*,*) 'logT > logT8', logT, logT8
                  if (logRho > logRho3) then
                     if (dbg) write(*,*) 'logRho > logRho3', logRho, logRho3
                     if (dbg) write(*,*) 'iregion = use_none 4'
                     iregion = use_none
                  else ! logRho is > logRho4
                     if (dbg) write(*,*) 'logRho <= logRho3', logRho, logRho3
                     d_dx_dlogRho = 1/(logRho3 - logRho4)
                     c_dx = (logRho - logRho4)*d_dx_dlogRho
                     d_dy_dlogT = 1/(logT8 - logT7)
                     c_dy = (logT - logT7)*d_dy_dlogT
                     if (dbg) write(*,*) 'iregion = blend_corner_out'
                     iregion = blend_corner_out
                  end if            
               end if
               
            else if (logRho >= logRho5 .or. logT > logT5) then
               if (dbg) write(*,*) 'iregion = use_all'
               iregion = use_all
                              
            else if (logT >= logT6) then
               if (logRho <= logRho6) then
                  d_dy_dlogT = 1/(logT6 - logT5)
                  c_dy = (logT - logT5)*d_dy_dlogT
                  if (dbg) write(*,*) 'iregion = blend_in_y'
                  iregion = blend_in_y
               else
                  d_dx_dlogRho = 1/(logRho5 - logRho6)
                  c_dx = (logRho - logRho6)*d_dx_dlogRho
                  d_dy_dlogT = 1/(logT5 - logT6)
                  c_dy = (logT - logT6)*d_dy_dlogT
                  if (dbg) write(*,*) 'iregion = blend_corner_in'
                  iregion = blend_corner_in
               end if
               
            else 
               if (dbg) write(*,*) 'logRho > logRho6', logRho, logRho6
               d_dx_dlogRho = 1/(logRho6 - logRho5)
               c_dx = (logRho - logRho5)*d_dx_dlogRho
               if (dbg) write(*,*) 'iregion = blend_in_x'
               iregion = blend_in_x
            end if

            if (dbg) stop 'determine_region'
            
         end subroutine determine_region_opal_scvh


         subroutine set_alfa_and_partials ! alfa = fraction other
            logical, parameter :: dbg = .false.
            
            real(dp) :: zfactor
            
            include 'formats'
            
            d_alfa_dlogT = 0d0
            d_alfa_dlogRho = 0
            
            if (iregion == use_none .or. Z >= Z_all_HELM) then
               if (dbg) write(*,*) 'iregion == use_none'
               alfa = 1
            else if (iregion == use_all) then
               if (dbg) write(*,*) 'iregion == use_all'
               alfa = 0
            else if (iregion == blend_in_y) then
               if (dbg) write(*,*) 'iregion == blend_in_y'
               alfa = c_dy
               d_alfa_dlogT = d_dy_dlogT
            else if (iregion == blend_in_x) then
               if (dbg) write(*,*) 'iregion == blend_in_x'
               alfa = c_dx
               d_alfa_dlogRho = d_dx_dlogRho
            else if (iregion == blend_diagonal) then
               if (dbg) write(*,*) 'iregion == blend_diagonal'
               alfa = c_dx
               d_alfa_dlogRho = d_dx_dlogRho
               d_alfa_dlogT = d_dy_dlogT
            else if (iregion == blend_corner_out) then
               if (dbg) write(*,*) 'iregion == blend_corner_out'
               alfa = sqrt(c_dx*c_dx + c_dy*c_dy)
               if (alfa >= 1d0) then
                  alfa = 1
               else if (alfa < 1d-10) then
                  alfa = 0
               else
                  d_alfa_dlogT = c_dy*d_dy_dlogT/alfa                  
                  d_alfa_dlogRho = c_dx*d_dx_dlogRho/alfa                                   
               end if
            else if (iregion == blend_corner_in) then
               if (dbg) write(*,*) 'iregion == blend_corner_in'
               beta = sqrt(c_dx*c_dx + c_dy*c_dy)
               alfa = 1d0 - beta
               if (alfa >= 1d0) then
                  alfa = 1d0
               else if (alfa < 1d-10) then
                  alfa = 0
               else
                  d_alfa_dlogT = -c_dy*d_dy_dlogT/beta                  
                  d_alfa_dlogRho = -c_dx*d_dx_dlogRho/beta                  
               end if
            else
               ierr = -1
               return
            end if

            if (Z > Z_no_HELM .and. Z < Z_all_HELM .and. alfa < 1d0) then
               ! reduce alfa to reduce the HELM fraction
               zfactor = (Z - Z_no_HELM)/(Z_all_HELM - Z_no_HELM)
               alfa = alfa*zfactor
               d_alfa_dlogRho = d_alfa_dlogRho*zfactor
               d_alfa_dlogT = d_alfa_dlogT*zfactor
            end if
            
         end subroutine set_alfa_and_partials


      end subroutine get_opal_scvh_alfa_and_partials
      
      
      subroutine combine_for_eosdt( &
            get_1st, get_2nd, remaining_fraction, &
            alfa_in, d_alfa_dlogT_in, d_alfa_dlogRho_in, &
            rq, dbg, Z, X, abar, zbar, &
            species, chem_id, net_iso, xa, &
            rho, logRho, T, logT, &
            res, d_dlnd, d_dlnT, d_dxa, &
            skip, ierr)
         use eosdt_support, only : Do_Blend
         procedure (get_values_for_eosdt_interface), pointer :: get_1st, get_2nd
         type (EoS_General_Info), pointer :: rq
         logical, intent(in) :: dbg
         real(dp), intent(in) :: Z, X, abar, zbar, &
            alfa_in, d_alfa_dlogT_in, d_alfa_dlogRho_in, remaining_fraction
         integer, intent(in) :: species
         integer, pointer :: chem_id(:), net_iso(:)
         real(dp), intent(in) :: xa(:)
         real(dp), intent(in) :: rho, logRho, T, logT
         real(dp), intent(inout), dimension(nv) :: res, d_dlnd, d_dlnT
         real(dp), intent(inout), dimension(nv, species) :: d_dxa
         logical, intent(out) :: skip
         integer, intent(out) :: ierr
         
         real(dp), dimension(nv) :: &
            res_1, d_dlnd_1, d_dlnT_1, res_2, d_dlnd_2, d_dlnT_2
         real(dp), dimension(:,:), allocatable :: d_dxa_1, d_dxa_2
         real(dp) :: alfa, d_alfa_dlogT, d_alfa_dlogRho
         logical :: skip_1st, skip_2nd
         logical, parameter :: linear_blend = .false.
         
         include 'formats'
         
         ierr = 0
         skip = .false.
         
         allocate(d_dxa_1(nv, species), d_dxa_2(nv, species))

         alfa = alfa_in
         d_alfa_dlogT = d_alfa_dlogT_in
         d_alfa_dlogRho = d_alfa_dlogRho_in

         if (alfa == 0d0) then ! pure 1st
            call get_1st(rq% handle, dbg, &
               Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               rho, logRho, T, logT, remaining_fraction, &
               res, d_dlnd, d_dlnT, d_dxa, &
               skip_1st, ierr)
            if (ierr /= 0) then
               if (.not. rq% okay_to_convert_ierr_to_skip) return
               if (dbg) write(*,*) 'ierr => skip 1st in combine_for_eosdt'
               ierr = 0; skip_1st = .true.
            end if
            if (skip_1st) then ! switch to pure 2nd
               alfa = 1d0; d_alfa_dlogT = 0d0; d_alfa_dlogRho = 0d0
            else
               return
            end if
         end if
         
         if (alfa < 1d0) then ! some of 1st
            call get_1st(rq% handle, dbg, &
               Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               rho, logRho, T, logT, remaining_fraction, &
               res_2, d_dlnd_2, d_dlnT_2, d_dxa_2, &
               skip_1st, ierr)
            if (ierr /= 0) then
               if (.not. rq% okay_to_convert_ierr_to_skip) return
               if (dbg) write(*,*) 'ierr => skip 1st in combine_for_eosdt'
               ierr = 0; skip_1st = .true.
            end if
            if (skip_1st) then ! switch to pure 2nd
               alfa = 1d0; d_alfa_dlogT = 0d0; d_alfa_dlogRho = 0d0
            end if         
         end if
         
         if (alfa == 1d0) then ! no 1st
            call get_2nd(rq% handle, dbg, &
               Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               Rho, logRho, T, logT, remaining_fraction, &
               res, d_dlnd, d_dlnT, d_dxa, &
               skip_2nd, ierr)
            if (ierr /= 0) then
               if (.not. rq% okay_to_convert_ierr_to_skip) return
               if (dbg) write(*,*) 'ierr => skip 2nd in combine_for_eosdt'
               ierr = 0; skip_2nd = .true.
            end if
            if (skip_2nd) skip = .true.
            return
         end if
         
         ! blend 1st and 2nd
         
         call get_2nd( &
            rq% handle, dbg, Z, X, abar, zbar, &
            species, chem_id, net_iso, xa, &
            Rho, logRho, T, logT, remaining_fraction, &
            res_1, d_dlnd_1, d_dlnT_1, d_dxa_1, &
            skip_2nd, ierr)
         if (ierr /= 0) then
            if (.not. rq% okay_to_convert_ierr_to_skip) return
            if (dbg) write(*,*) 'ierr => skip 2nd in combine_for_eosdt'
            ierr = 0; skip_2nd = .true.
         end if
         if (skip_2nd) then
            skip = .true.
            return
         end if
         call Do_Blend( &
            rq, species, Rho, logRho, T, logT, &
            alfa, d_alfa_dlogT, d_alfa_dlogRho, linear_blend, &
            res_1, d_dlnd_1, d_dlnT_1, d_dxa_1, &
            res_2, d_dlnd_2, d_dlnT_2, d_dxa_2, &
            res, d_dlnd, d_dlnT, d_dxa)
                   
      end subroutine combine_for_eosdt


      subroutine get1_for_eosdt( &
            handle, which_eosdt, dbg, Z, X, abar, zbar, &
            species, chem_id, net_iso, xa, &
            rho, logRho, T, logT, remaining_fraction, &
            res, d_dlnd, d_dlnT, d_dxa, &
            skip, ierr)
         use chem_def, only: chem_isos
         integer, intent(in) :: handle
         integer, intent(in) :: which_eosdt
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
         real(dp), dimension(nv) :: d_dX, d_dZ
         logical, intent(out) :: skip
         integer, intent(out) :: ierr
         type (EoS_General_Info), pointer :: rq
         type (DT_xz_Info), pointer :: xz
         integer :: i
         rq => eos_handles(handle)
         if (which_eosdt == eosdt_max_FreeEOS) then
            xz => FreeEOS_xz_struct
         else
            xz => eosDT_xz_struct
         end if
         call Get1_eosdt_Results( &
            rq, which_eosdt, xz, Z, X, Rho, logRho, T, logT, &
            res, d_dlnd, d_dlnT, d_dX, d_dZ, ierr)
         do i=1,species
            select case(chem_isos% Z(chem_id(i))) ! charge
            case (1) ! X
               d_dxa(:,i) = d_dX
            case (2) ! Y
               d_dxa(:,i) = 0
            case default ! Z
               d_dxa(:,i) = d_dZ
            end select
         end do
         skip = .false.
      end subroutine get1_for_eosdt


      subroutine Get1_eosdt_Results( & ! blend in Z
               rq, which_eosdt, xz, Z, X, Rho, logRho, T, logT, &
               res, dlnd, dlnT, dX, dZ, ierr)
         use chem_def
         type (EoS_General_Info), pointer :: rq
         integer, intent(in) :: which_eosdt
         type (DT_xz_Info), pointer :: xz
         real(dp), intent(in) :: Z, X, Rho, logRho, T, logT
         real(dp), intent(inout), dimension(nv) :: res, dlnd, dlnT, dX, dZ
         integer, intent(out) :: ierr

         real(dp), dimension(nv, 2) :: res_zx, dlnd_zx, dlnT_zx, dX_zx
         real(dp) :: denom, c(2), dcdZ(2), tiny

         integer :: iz, j, ci
         
         include 'formats'

         ierr = 0
         tiny = rq% tiny_fuzz
         
         if (xz% nZs < 3) then
            write(*, *) 'error: Get1_eosdt_Results assumes nZs >= 3'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         if (xz% Zs(1) /= 0) then
            write(*, *) 'error: Get1_eosdt_Results assumes eos_Zs(1) == 0'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         if (abs(xz% Zs(1) - 2*xz% Zs(2) + xz% Zs(3)) > tiny) then
            write(*, *) 'error: Get1_eosdt_Results assumes equal spaced Zs(1:3)'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         if (Z <= max(1d-20,xz% Zs(1))) then
            call Get1_eosdt_for_X( &
                  rq, which_eosdt, xz, 1, X, &
                  Rho, logRho, T, logT, &
                  res, dlnd, dlnT, dX, ierr)
            dZ = 0
            return
         end if

         ! zero these for now
         res_zx(i_phase:i_latent_ddlnRho,:) = 0d0
         res_zx(i_frac:i_frac+num_eos_frac_results-1,:) = 0d0

         dlnd_zx(i_phase:i_latent_ddlnRho,:) = 0d0
         dlnd_zx(i_frac:i_frac+num_eos_frac_results-1,:) = 0d0

         dlnT_zx(i_phase:i_latent_ddlnRho,:) = 0d0
         dlnT_zx(i_frac:i_frac+num_eos_frac_results-1,:) = 0d0        
         
         if (Z >= xz% Zs(xz% nZs)) then
            call Get1_eosdt_for_X( &
                  rq, which_eosdt, xz, xz% nZs, X, &
                  Rho, logRho, T, logT, &
                  res, dlnd, dlnT, dX, ierr)
            dZ = 0
            return
         end if

         do iz = 2, xz% nZs
            if (Z < xz% Zs(iz)) then
               call do_interp2(iz-1,iz,ierr)
               if (ierr /= 0) return
               exit
            end if
         end do
         
         do j=1,nv
         
            res(j) = c(1)*res_zx(j,1) + c(2)*res_zx(j,2)
            
            dlnd(j) = &
               c(1)*dlnd_zx(j,1) + c(2)*dlnd_zx(j,2)
               
            dlnT(j) = &
               c(1)*dlnT_zx(j,1) + c(2)*dlnT_zx(j,2)

            dX(j) = &
               c(1)*dX_zx(j,1) + c(2)*dX_zx(j,2)

            dZ(j) = &
               dcdZ(1)*res_zx(j,1) + dcdZ(2)*res_zx(j,2)

         end do
            
         contains
         
         subroutine do_interp2(iz1, iz2, ierr)
            integer, intent(in) :: iz1, iz2
            integer, intent(out) :: ierr
            real(dp) :: Z1, Z2
            include 'formats'
            ierr = 0
            Z1 = xz% Zs(iz1)
            Z2 = xz% Zs(iz2)
            c(2) = (Z - Z1) / (Z2 - Z1)
            c(1) = 1d0 - c(2)
            dcdZ(2) = 1d0/(Z2 - Z1)
            dcdZ(1) = -dcdZ(2)
            call Get1_eosdt_for_X( &
                  rq, which_eosdt, xz, iz1, X, Rho, logRho, T, logT, &
                  res_zx(:,1), dlnd_zx(:,1), dlnT_zx(:,1), dX_zx(:,1), &
                  ierr)
            if (ierr /= 0) return
            call Get1_eosdt_for_X( &
                  rq, which_eosdt, xz, iz2, X, Rho, logRho, T, logT, &
                  res_zx(:,2), dlnd_zx(:,2), dlnT_zx(:,2), dX_zx(:,2), &
                  ierr)
            if (ierr /= 0) return
         end subroutine do_interp2
     
      end subroutine Get1_eosdt_Results

      
      subroutine Get1_eosdt_for_X( &
               rq, which_eosdt, xz, iz, X, Rho, logRho, T, logT, &
               res, dlnd, dlnT, d_dX, ierr)
         type (EoS_General_Info), pointer :: rq
         integer, intent(in) :: which_eosdt
         type (DT_xz_Info), pointer :: xz
         integer, intent(in) :: iz ! the index in eos_Zs
         real(dp), intent(in) :: X, Rho, logRho, T, logT
         real(dp), intent(inout), dimension(nv) :: &
            res, dlnd, dlnT, d_dX
         integer, intent(out) :: ierr

         real(dp), dimension(nv, 4) :: &
               res_zx, dlnd_zx, dlnT_zx
         real(dp) :: dX, dX1, dX2, dX3, c(4), dcdX(4), denom, delX, coef, dcoef_dX, alfa, beta, dalfa_dX, dbeta_dX, tiny
         character (len=256) :: message
         integer :: ix, ix_lo, ix_hi, j, num_Xs
         logical, parameter :: dbg_for_X = dbg ! .or. .true.
         logical :: what_we_use_is_equal_spaced
         
         include 'formats'
         
         ierr = 0
         tiny = rq% tiny_fuzz
         
         num_Xs = xz% nXs_for_Z(iz)
         
         if (xz% Xs_for_Z(1,iz) /= 0d0) then
            write(*, *) 'error: Get1_eosdt_for_X assumes xz% nXs_for_Z(1) == 0'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         if (X < tiny .or. num_Xs == 1) then
            call Get1_eosdt_XTable_Results( &
               rq, which_eosdt, 1, iz, Rho, logRho, T, logT, &
               res, dlnd, dlnT, ierr)
            d_dX = 0
            return
         end if
         
         if (X >= xz% Xs_for_Z(num_Xs,iz)) then

            call Get1_eosdt_XTable_Results( &
               rq, which_eosdt, num_Xs, iz, Rho, logRho, T, logT, &
               res, dlnd, dlnT, ierr)
            d_dX = 0

            if (is_bad(res(i_lnS))) then
               ierr = -1
               if (.not. stop_for_is_bad) return
               write(*,1) 'res(i_lnS), logRho, logT', res(i_lnS), logRho, logT
               stop 'Get1_eosdt_for_X num_Xs'
            end if
            
            return
         end if

         if (rq% eosDT_use_linear_interp_for_X .or. num_Xs == 2) then
            call do_linear
            return
         end if
         
         ix_hi = -1
         if (X <= xz% Xs_for_Z(2,iz)) then
            ix_lo = 1; ix_hi = 3
         else if (X >= xz% Xs_for_Z(num_Xs-1,iz)) then
            ix_lo = num_Xs-2; ix_hi = num_Xs
         else
            do ix = 3, num_Xs-1
               if (X <= xz% Xs_for_Z(ix,iz)) then
                  ix_lo = ix-2; ix_hi = ix+1; exit
               end if
            end do
         end if
         
         if (ix_hi < 0) then
            write(*, *) 'X', X
            write(*, *) 'ix_lo', ix_lo
            write(*, *) 'ix_hi', ix_hi
            write(*, *) 'error: Get1_eosdt_for_X logic bug'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         if (dbg_for_X) then
            write(*, *) 'X', X
            write(*, *) 'ix_lo', ix_lo
            write(*, *) 'ix_hi', ix_hi
         end if

         what_we_use_is_equal_spaced = .true.
         dX1 = xz% Xs_for_Z(ix_lo+1,iz)-xz% Xs_for_Z(ix_lo,iz)
         dX2 = xz% Xs_for_Z(ix_lo+2,iz)-xz% Xs_for_Z(ix_lo+1,iz)
         if (ix_hi-ix_lo==2) then ! check that the 3 table X's are equal spaced
            if (abs(dX1 - dX2) > tiny) what_we_use_is_equal_spaced = .false.
         else ! check that the 4 table X's are equal spaced 
            dX3 = xz% Xs_for_Z(ix_hi,iz)-xz% Xs_for_Z(ix_lo+2,iz)
            if (abs(dX1 - dX2) > tiny .or. abs(dX2 - dX3) > tiny) &
               what_we_use_is_equal_spaced = .false.
         end if
         
         if (.not. what_we_use_is_equal_spaced) then
            call do_linear
            if (is_bad(d_dX(1))) then
               stop 'Get1_eosdt_for_X bad d_dX; linear'
            end if
            return
         end if
         
         do ix=ix_lo, ix_hi
            j = ix-ix_lo+1
            call Get1_eosdt_XTable_Results( &
               rq, which_eosdt, ix, iz, Rho, logRho, T, logT, &
               res_zx(:, j), dlnd_zx(:, j), dlnT_zx(:, j), &
               ierr)
            if (ierr /= 0) return
         end do

         ! zero these for now
         res_zx(i_phase:i_latent_ddlnRho,:) = 0d0
         res_zx(i_frac:i_frac+num_eos_frac_results-1,:) = 0d0

         dlnd_zx(i_phase:i_latent_ddlnRho,:) = 0d0
         dlnd_zx(i_frac:i_frac+num_eos_frac_results-1,:) = 0d0

         dlnT_zx(i_phase:i_latent_ddlnRho,:) = 0d0
         dlnT_zx(i_frac:i_frac+num_eos_frac_results-1,:) = 0d0

         
         delX = X - xz% Xs_for_Z(ix_lo,iz)
         dX = dX1
         
         if (ix_hi-ix_lo==2) then
         
            denom = 2*dX*dX
            c(1) = (2*dX*dX - 3*dX*delX + delX*delX)/denom
            c(2) = 2*(2*dX-delX)*delX/denom
            c(3) = delX*(delX-dX)/denom
            res(:) = c(1)*res_zx(:, 1) + c(2)*res_zx(:, 2) + c(3)*res_zx(:, 3)
            
            dlnd(:) = &
               c(1)*dlnd_zx(:,1) + &
               c(2)*dlnd_zx(:,2) + &
               c(3)*dlnd_zx(:,3)
            dlnT(:) = &
               c(1)*dlnT_zx(:,1) + &
               c(2)*dlnT_zx(:,2) + &
               c(3)*dlnT_zx(:,3)

            dcdx(1) = (-3*dX + 2*delX)/denom
            dcdx(2) = 2*(2*dX-2*delX)/denom
            dcdx(3) = (2*delX-dX)/denom

            d_dX(:) = &
               dcdX(1)*res_zx(:,1) + &
               dcdX(2)*res_zx(:,2) + &
               dcdX(3)*res_zx(:,3)

            if (is_bad(d_dX(1))) then
               stop 'Get1_eosdt_for_X bad d_dX; 3'
            end if
            
         else
         
            coef = (X-xz% Xs_for_Z(ix_lo+1,iz))/dX 
            ! coef = fractional location of X between 2nd and 3rd X's for fit.
            ! coef is weight for the quadratic based on points 2, 3, 4 of fit.
            ! (1-coef) is weight for quadratic based on points 1, 2, 3 of fit.
            coef = min(1d0,max(0d0,coef))
            c(1) = -coef*(coef-1)*(coef-1)/2
            c(2) = (2 - coef*coef*(5 - 3*coef))/2
            c(3) = coef*(1 + coef*(4 - 3*coef))/2
            c(4) = coef*coef*(coef-1)/2
            res(:) = c(1)*res_zx(:, 1) + &
                        (c(2)*res_zx(:, 2) + &
                           (c(3)*res_zx(:, 3) + &
                              c(4)*res_zx(:, 4)))
            
            dlnd(:) = &
               c(1)*dlnd_zx(:, 1) + &
                  (c(2)*dlnd_zx(:, 2) + &
                     (c(3)*dlnd_zx(:, 3) + &
                           c(4)*dlnd_zx(:, 4)))
            dlnT(:) = &
               c(1)*dlnT_zx(:, 1) + &
                  (c(2)*dlnT_zx(:, 2) + &
                     (c(3)*dlnT_zx(:, 3) + &
                           c(4)*dlnT_zx(:, 4)))

            dcoef_dX = 1d0/dX
            dcdX = 0
            dcdX(1) = -(3*coef*coef-4*coef+1)/2*dcoef_dX
            dcdX(2) = (9*coef*coef-10*coef)/2*dcoef_dX
            dcdX(3) = -(9*coef*coef-8*coef-1)/2*dcoef_dX
            dcdX(4) = coef*(3*coef-2)/2*dcoef_dX

            d_dX(:) = &
               dcdX(1)*res_zx(:,1) + &
               dcdX(2)*res_zx(:,2) + &
               dcdX(3)*res_zx(:,3) + &
               dcdX(4)*res_zx(:,4)

            if (is_bad(d_dX(1))) then
               stop 'Get1_eosdt_for_X bad d_dX; 4'
            end if

         end if
         
         contains
         
         subroutine do_linear
         
            do ix = 2, num_Xs
               if (xz% Xs_for_Z(ix,iz) >= X) exit
            end do
         
            j = 1
            call Get1_eosdt_XTable_Results( &
               rq, which_eosdt, ix-1, iz, Rho, logRho, T, logT, &
               res_zx(:,j), dlnd_zx(:,j), dlnT_zx(:,j), &
               ierr)
            if (ierr /= 0) then
               if (.not. stop_for_is_bad) return
               stop 'Get1_eosdt_for_X'
            end if
         
            j = 2
            call Get1_eosdt_XTable_Results( &
               rq, which_eosdt, ix, iz, Rho, logRho, T, logT, &
               res_zx(:,j), dlnd_zx(:,j), dlnT_zx(:,j), &
               ierr)
            if (ierr /= 0) then
               if (.not. stop_for_is_bad) return
               stop 'Get1_eosdt_for_X'
            end if

            ! zero these for now
            res_zx(i_phase:i_latent_ddlnRho,:) = 0d0
            res_zx(i_frac:i_frac+num_eos_frac_results-1,:) = 0d0

            dlnd_zx(i_phase:i_latent_ddlnRho,:) = 0d0
            dlnd_zx(i_frac:i_frac+num_eos_frac_results-1,:) = 0d0

            dlnT_zx(i_phase:i_latent_ddlnRho,:) = 0d0
            dlnT_zx(i_frac:i_frac+num_eos_frac_results-1,:) = 0d0

         
            alfa = (X - xz% Xs_for_Z(ix,iz))/(xz% Xs_for_Z(ix-1,iz) - xz% Xs_for_Z(ix,iz))
            beta = 1d0 - alfa

            dalfa_dX = 1d0 / (xz% Xs_for_Z(ix-1,iz) - xz% Xs_for_Z(ix,iz))
            dbeta_dX = -dalfa_dX
         
            do j=1,nv
            
               res(j) = alfa*res_zx(j,1) + beta*res_zx(j,2)
               
               dlnd(j) = &
                  alfa*dlnd_zx(j,1) + beta*dlnd_zx(j,2)
                  
               dlnT(j) = &
                  alfa*dlnT_zx(j,1) + beta*dlnT_zx(j,2)

               d_dX(j) = &
                  dalfa_dX*res_zx(j,1) + dbeta_dX*res_zx(j,2)

            end do
         
         end subroutine do_linear
         
      end subroutine Get1_eosdt_for_X

      
      subroutine Locate_logQ(rq, ep, logQ, iQ, logQ0, logQ1, ierr)
         type (EoS_General_Info), pointer :: rq
         type (EosDT_xz_Info), pointer :: ep
         real(dp), intent(inout) :: logQ
         integer, intent(out) :: iQ
         real(dp), intent(out) :: logQ0, logQ1
         integer, intent(out) :: ierr      
         ierr = 0
         iQ = int((logQ - ep% logQ_min)/ep% del_logQ + 1d-4) + 1         
         if (iQ < 1 .or. iQ >= ep% num_logQs) then            
            if (iQ < 1) then
               iQ = 1
               logQ0 = ep% logQ_min
               logQ1 = logQ0 + ep% del_logQ
               logQ = logQ0
               if (return_ierr_beyond_table_bounds) ierr = -1
            else
               iQ = ep% num_logQs-1
               logQ0 = ep% logQ_min + (iQ-1)*ep% del_logQ
               logQ1 = logQ0 + ep% del_logQ
               logQ = logQ1
               if (return_ierr_beyond_table_bounds) ierr = -1
            end if            
         else         
            logQ0 = ep% logQ_min + (iQ-1)*ep% del_logQ
            logQ1 = logQ0 + ep% del_logQ
         end if
      end subroutine Locate_logQ
      
      
      subroutine Locate_logT(rq, ep, logT, iT, logT0, logT1, ierr)
         type (EoS_General_Info), pointer :: rq
         type (EosDT_xz_Info), pointer :: ep
         real(dp), intent(inout) :: logT
         integer, intent(out) :: iT
         real(dp), intent(out) :: logT0, logT1
         integer, intent(out) :: ierr      
         ierr = 0
         iT = int((logT - ep% logT_min)/ep% del_logT + 1d-4) + 1        
         if (iT < 1 .or. iT >= ep% num_logTs) then           
            if (iT < 1) then
               iT = 1
               logT0 = ep% logT_min
               logT1 = logT0 + ep% del_logT
               logT = logT0
               if (return_ierr_beyond_table_bounds) ierr = -1
            else
               iT = ep% num_logTs-1
               logT0 = ep% logT_min + (iT-1)*ep% del_logT
               logT1 = logT0 + ep% del_logT
               logT = logT1
               if (return_ierr_beyond_table_bounds) ierr = -1
            end if            
         else         
            logT0 = ep% logT_min + (iT-1)*ep% del_logT
            logT1 = logT0 + ep% del_logT
         end if
      end subroutine Locate_logT
      
      
      subroutine Get1_eosdt_XTable_Results( &
            rq, which_eosdt, ix, iz, Rho, logRho_in, T, logT_in, &
            res, d_dlnd, d_dlnT, ierr)
         use eosdt_support, only: Do_EoS_Interpolations
         type (EoS_General_Info), pointer :: rq
         integer, intent(in) :: which_eosdt
         integer, intent(in) :: ix, iz
         real(dp), intent(in) :: Rho, logRho_in, T, logT_in
         real(dp), intent(inout), dimension(nv) :: res, d_dlnd, d_dlnT
         integer, intent(out) :: ierr
         
         real(dp), parameter :: ln10sq = ln10*ln10
         real(dp) :: &
            fval(nv), df_dx(nv), df_dy(nv), &
            df_dlnd(nv), df_dlnT(nv), &
            energy, entropy, P, Pgas, Prad, x, y, &
            dx_dlnd, dx_dlnT, dy_dlnd, dy_dlnT, &
            chiT, chiRho, Cv, gamma1, numeric, &
            dS_dlnd, dS_dlnT, dE_dlnd, dE_dlnT, &
            dPgas_dlnd, dPgas_dlnT, dP_dlnd, dP_dlnT
         real(dp) :: logQ0, logQ1, logT0, logT1, logRho0, logRho1
         integer :: iQ, jtemp, k, j, irho
         type (EosDT_xz_Info), pointer :: ep
         logical, parameter :: show = .false.
         real(dp) :: logRho, logT, logQ
         
         include 'formats'

         logRho = logRho_in
         logT = logT_in
         logQ = logRho - 2*logT + 12

         ierr = 0 
         call load_single_eosDT_table_by_id(rq, which_eosdt, ep, ix, iz, ierr)
         if (ierr /= 0) return

         call Locate_logQ(rq, ep, logQ, iQ, logQ0, logQ1, ierr)
         if (ierr /= 0) then
            write(*,1) 'eosDT failed in Locate_logQ', logQ
            return
         end if
      
         call Locate_logT(rq, ep, logT, jtemp, logT0, logT1, ierr)
         if (ierr /= 0) then
            write(*,1) 'eosDT failed in Locate_logT', logT
            return
         end if
         
         call Do_EoS_Interpolations( &
            1, nv, nv, ep% num_logQs, ep% logQs, ep% num_logTs, ep% logTs, &
            ep% tbl1, iQ, jtemp, logQ0, logQ, logQ1, logT0, logT, logT1, &
            fval, df_dx, df_dy, ierr)               
         if (ierr /= 0) then
            write(*,1) 'failed in Do_EoS_Interpolations'
            return
         end if
         
         if (is_bad(fval(i_lnS))) then
            ierr = -1
            if (.not. stop_for_is_bad) return
            write(*,1) 'fval(i_lnS), logRho, logT', fval(i_lnS), logRho, logT
            stop 'after Do_Interp_with_2nd_derivs'
         end if
         
         res(i_lnPgas) = fval(i_lnPgas)
         res(i_lnE) = fval(i_lnE)
         res(i_lnS) = fval(i_lnS)
         
         if (is_bad(res(i_lnS))) then
            ierr = -1
            if (.not. stop_for_is_bad) return
            write(*,1) 'res(i_lnS), logRho, logT', res(i_lnS), logRho, logT
            stop 'after interpolation'
         end if
         
         if (is_bad(res(i_lnS)) .or. res(i_lnS) > ln10*100) then
            ierr = -1
            if (.not. stop_for_is_bad) return
            write(*,1) 'res(i_lnS), logRho, logT', res(i_lnS), logRho, logT
            stop 'after interpolation'
         end if
         
         res(i_grad_ad) = fval(i_grad_ad)
         res(i_chiRho) = fval(i_chiRho)
         res(i_chiT) = fval(i_chiT)
      
         res(i_Cp) = fval(i_Cp)
         res(i_Cv) = fval(i_Cv)
      
         res(i_dE_dRho) = fval(i_dE_dRho)
         res(i_dS_dT) = fval(i_dS_dT)
         res(i_dS_dRho) = fval(i_dS_dRho)
      
         res(i_mu) = fval(i_mu)
         res(i_lnfree_e) = fval(i_lnfree_e)
         res(i_gamma1) = fval(i_gamma1)
         res(i_gamma3) = fval(i_gamma3)
         res(i_eta) = fval(i_eta)

         ! convert df_dx and df_dy to df_dlogRho_c_T and df_dlogT_c_Rho
         
         ! df_dx is df_dlogQ at const T
         ! df_dy is df_dlogT_c_Rho at const Q
         ! logQ = logRho - 2*logT + 12
         
         ! f = f(logQ(logRho,logT),logT)
         ! df/dlogRho|T = df/dlogQ|T * dlogQ/dlogRho|T = df_dx
         ! df/dlogT|Rho = df/dlogT|Q + df/dlogQ|T * dlogQ/dlogT|Rho = df_dy - 2*df_dx
            
         do k=1,nv
            df_dlnd(k) = df_dx(k)/ln10
            df_dlnT(k) = df_dy(k)/ln10 - 2d0*df_dlnd(k)
         end do

         d_dlnd(i_lnPgas) = df_dlnd(i_lnPgas)
         d_dlnd(i_lnE) = df_dlnd(i_lnE)
         d_dlnd(i_lnS) = df_dlnd(i_lnS)
         d_dlnd(i_grad_ad) = df_dlnd(i_grad_ad)
         d_dlnd(i_chiRho) = df_dlnd(i_chiRho)
         d_dlnd(i_chiT) = df_dlnd(i_chiT)
      
         d_dlnd(i_Cp) = df_dlnd(i_Cp)
         d_dlnd(i_Cv) = df_dlnd(i_Cv)
         d_dlnd(i_dE_dRho) = df_dlnd(i_dE_dRho)
         d_dlnd(i_dS_dT) = df_dlnd(i_dS_dT)
         d_dlnd(i_dS_dRho) = df_dlnd(i_dS_dRho)
         d_dlnd(i_mu) = df_dlnd(i_mu)
         d_dlnd(i_lnfree_e) = df_dlnd(i_lnfree_e)
         d_dlnd(i_gamma1) = df_dlnd(i_gamma1)
         d_dlnd(i_gamma3) = df_dlnd(i_gamma3)
         d_dlnd(i_eta) = df_dlnd(i_eta)
   
         d_dlnT(i_lnPgas) = df_dlnT(i_lnPgas)
         d_dlnT(i_lnE) = df_dlnT(i_lnE)
         d_dlnT(i_lnS) = df_dlnT(i_lnS)
         d_dlnT(i_grad_ad) = df_dlnT(i_grad_ad)
         d_dlnT(i_chiRho) = df_dlnT(i_chiRho)
         d_dlnT(i_chiT) = df_dlnT(i_chiT)
         d_dlnT(i_Cp) = df_dlnT(i_Cp)
         d_dlnT(i_Cv) = df_dlnT(i_Cv)
         d_dlnT(i_dE_dRho) = df_dlnT(i_dE_dRho)
         d_dlnT(i_dS_dT) = df_dlnT(i_dS_dT)
         d_dlnT(i_dS_dRho) = df_dlnT(i_dS_dRho)
         d_dlnT(i_mu) = df_dlnT(i_mu)
         d_dlnT(i_lnfree_e) = df_dlnT(i_lnfree_e)
         d_dlnT(i_gamma1) = df_dlnT(i_gamma1)
         d_dlnT(i_gamma3) = df_dlnT(i_gamma3)
         d_dlnT(i_eta) = df_dlnT(i_eta)
         
         if (is_bad(d_dlnd(i_lnS)) .or. is_bad(d_dlnT(i_lnS))) then
            ierr = -1
            if (.not. stop_for_is_bad) return
            write(*,1) 'fval(i_lnS)', fval(i_lnS)
            write(*,1) 'd_dlnd(i_lnS)', d_dlnd(i_lnS)
            write(*,1) 'd_dlnT(i_lnS)', d_dlnT(i_lnS)
            stop 'Get1_eosdt_XTable_Results'
         end if

      end subroutine Get1_eosdt_XTable_Results


      subroutine get_T( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               logRho, which_other, other_value, &
               logT_tol, other_tol, max_iter, logT_guess, &
               logT_bnd1, logT_bnd2,  other_at_bnd1, other_at_bnd2, &
               logT_result, res, d_dlnRho_c_T, d_dlnT_c_Rho, d_dxa_c_TRho, &
               eos_calls, ierr)
         
         integer, intent(in) :: handle

         real(dp), intent(in) :: Z ! the metals mass fraction
         real(dp), intent(in) :: X ! the hydrogen mass fraction
            
         real(dp), intent(in) :: abar, zbar
         
         integer, intent(in) :: species
         integer, pointer :: chem_id(:)    
         integer, pointer :: net_iso(:)
         real(dp), intent(in) :: xa(:)
         
         real(dp), intent(in) :: logRho ! log10 of density
         integer, intent(in) :: which_other ! from eos_def.  e.g., i_P for pressure
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
         real(dp), intent(inout), dimension(nv) :: res, d_dlnRho_c_T, d_dlnT_c_Rho
         real(dp), intent(inout), dimension(:,:) :: d_dxa_c_TRho

         integer, intent(out) :: eos_calls
         integer, intent(out) :: ierr ! 0 means AOK.

         logical, parameter :: doing_Rho = .false.
         
         call do_safe_get_Rho_T( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               logRho, which_other, other_value, doing_Rho, &
               logT_guess, logT_result, logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
               logT_tol, other_tol, max_iter, res, d_dlnRho_c_T, d_dlnT_c_Rho, d_dxa_c_TRho, &
               eos_calls, ierr)
      
      end subroutine get_T
      

      subroutine get_Rho( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               logT, which_other, other_value, &
               logRho_tol, other_tol, max_iter, logRho_guess, &
               logRho_bnd1, logRho_bnd2, other_at_bnd1, other_at_bnd2, &
               logRho_result, res, d_dlnRho_c_T, d_dlnT_c_Rho, d_dxa_c_TRho, &
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

         integer, intent(in) :: which_other ! from eos_def.
         real(dp), intent(in) :: other_value ! desired value for the other variable
         real(dp), intent(in) :: other_tol
         
         real(dp), intent(in) :: logRho_tol

         integer, intent(in) :: max_iter ! max number of Newton iterations        

         real(dp), intent(in) :: logRho_guess
         real(dp), intent(in) :: logRho_bnd1, logRho_bnd2 ! bounds for logrho
            ! set to arg_not_provided if do not know bounds
         real(dp), intent(in) :: other_at_bnd1, other_at_bnd2 ! values at bounds
            ! if don't know these values, just set to arg_not_provided (defined in c_def)
            
         real(dp), intent(out) :: logRho_result
         real(dp), intent(inout), dimension(nv) :: res, d_dlnRho_c_T, d_dlnT_c_Rho
         real(dp), intent(inout), dimension(:,:) :: d_dxa_c_TRho

         integer, intent(out) :: eos_calls
         integer, intent(out) :: ierr ! 0 means AOK.

         logical, parameter :: doing_Rho = .true.
         real(dp) :: Prad
         
         call do_safe_get_Rho_T( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               logT, which_other, other_value, doing_Rho, &
               logRho_guess, logRho_result, logRho_bnd1, logRho_bnd2, other_at_bnd1, other_at_bnd2, &
               logRho_tol, other_tol, max_iter, res, d_dlnRho_c_T, d_dlnT_c_Rho, d_dxa_c_TRho, &
               eos_calls, ierr)

      end subroutine get_Rho
      

      subroutine do_safe_get_Rho_T( &
               handle, Z, XH1, abar, zbar, &
               species, chem_id, net_iso, xa, &
               the_other_log, which_other, other_value, doing_Rho, &
               initial_guess, x, xbnd1, xbnd2, other_at_bnd1, other_at_bnd2, &
               xacc, yacc, ntry, &
               res, d_dlnRho_c_T, d_dlnT_c_Rho, d_dxa_c_TRho, &
               eos_calls, ierr)
         use const_def
         use chem_def, only: num_chem_isos
         use num_lib, only: safe_root_with_guess
         integer, intent(in) :: handle
         real(dp), intent(in) :: Z, XH1, abar, zbar
         integer, intent(in) :: species
         integer, pointer :: chem_id(:)    
         integer, pointer :: net_iso(:)
         real(dp), intent(in) :: xa(:)
         integer, intent(in) :: which_other ! 0 means total P
         real(dp), intent(in) :: other_value
         logical, intent(in) :: doing_Rho
         real(dp), intent(in) :: initial_guess ! for x
         real(dp), intent(out) :: x ! if doing_Rho, then logRho, else logT
         real(dp), intent(in) :: the_other_log
         real(dp), intent(in) :: xbnd1, xbnd2, other_at_bnd1, other_at_bnd2
         real(dp), intent(in) :: xacc, yacc ! tolerances
         integer, intent(in) :: ntry ! max number of iterations        
         real(dp), intent(inout), dimension(nv) :: res, d_dlnRho_c_T, d_dlnT_c_Rho
         real(dp), dimension(:,:) :: d_dxa_c_TRho
         integer, intent(out) :: eos_calls, ierr
         
         integer :: i, j, ix, iz
         integer, parameter :: lrpar = 0, lipar = 0, newt_imax = 6
         real(dp), parameter :: dx = 0.1d0
         integer, pointer :: ipar(:)
         real(dp), pointer :: rpar(:)
         real(dp) :: the_other_val, logRho, logT, rho, T, x1, x3, y1, y3
         type (EoS_General_Info), pointer :: rq

         include 'formats'
         
         ierr = 0
         
         call get_eos_ptr(handle, rq, ierr)
         if (ierr /= 0) then
            write(*, *) 'get_eos_ptr returned ierr', ierr
            return
         end if

         x1 = arg_not_provided
         x3 = arg_not_provided
         y1 = arg_not_provided
         y3 = arg_not_provided
         
         eos_calls = 0
         the_other_val = exp10(the_other_log)
         nullify(ipar, rpar)

         x = safe_root_with_guess( &
            f, initial_guess, dx, x1, x3, y1, y3, &
            min(ntry,newt_imax), ntry, xacc, yacc, &
            lrpar, rpar, lipar, ipar, ierr)
         
         contains
         
         real(dp) function f(x, dfdx, lrpar, rpar, lipar, ipar, ierr)
            ! returns with ierr = 0 if was able to evaluate f and df/dx at x
            ! if df/dx not available, it is okay to set it to 0
            use const_def, only: dp
            integer, intent(in) :: lrpar, lipar
            real(dp), intent(in) :: x
            real(dp), intent(out) :: dfdx
            integer, intent(inout), pointer :: ipar(:) ! (lipar)
            real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
            integer, intent(out) :: ierr
            
            real(dp) :: Pgas, Prad, energy, entropy, dPgas_dlnT, dPrad_dlnT, &
               dPgas_dlnRho, erad, egas, derad_dlnT, degas_dlnT, derad_dlnRho
            
            include 'formats'
            ierr = 0
            eos_calls = eos_calls + 1
            f = 0; dfdx = 0
            
            if (x > 50d0) then
               ierr = -1
               return
            end if
            
            if (doing_Rho) then
               logRho = x
               rho = exp10(logRho)
               logT = the_other_log
               T = the_other_val
            else
               logT = x
               T = exp10(logT)
               logRho = the_other_log
               rho = the_other_val
            end if
            
            call Get_eosDT_Results(rq, Z, XH1, abar, zbar, &
                  species, chem_id, net_iso, xa, &
                  rho, logRho, T, logT, &
                  res, d_dlnRho_c_T, d_dlnT_c_Rho, d_dxa_c_TRho, ierr)

            Pgas = exp(res(i_lnPgas))
            Prad = crad*T*T*T*T/3d0
            energy = exp(res(i_lnE))
            entropy = exp(res(i_lnS))

            if (ierr /= 0) then
               if (.false.) then
                  write(*,*) 'Get_eosDT_Results returned ierr', ierr
                  write(*,1) 'Z', Z
                  write(*,1) 'XH1', XH1
                  write(*,1) 'abar', abar
                  write(*,1) 'zbar', zbar
                  write(*,1) 'rho', rho
                  write(*,1) 'logRho', logRho
                  write(*,1) 'T', T
                  write(*,1) 'logT', logT
                  write(*,*)
                  stop 'do_safe_get_Rho_T'
               end if
               return
            end if
            
            if (is_bad(res(i_Cv))) then
               ierr = -1
               if (.not. stop_for_is_bad) return
               write(*,1) 'res(i_Cv)', res(i_Cv)
               write(*,1) 'Z', Z
               write(*,1) 'XH1', XH1
               write(*,1) 'abar', abar
               write(*,1) 'zbar', zbar
               write(*,1) 'rho', rho
               write(*,1) 'logRho', logRho
               write(*,1) 'T', T
               write(*,1) 'logT', logT
               write(*,*)
               stop 'do_safe_get_Rho_T'
            end if
            
            if (which_other == -1) then ! other_value is egas
               erad = crad*pow4(T)/rho
               egas = energy - erad
               f = egas - other_value
               if (doing_Rho) then 
                  derad_dlnRho = -erad
                  dfdx = energy*d_dlnRho_c_T(i_lnE)*ln10 - derad_dlnRho
               else
                  derad_dlnT = 4d0*erad
                  degas_dlnT = energy*d_dlnT_c_Rho(i_lnE) - derad_dlnT
                  dfdx = degas_dlnT*ln10
               end if
            else if (which_other == 0) then ! other_value is log10P
               f = log10(Pgas + Prad) - other_value
               if (doing_Rho) then 
                  dPgas_dlnRho = Pgas*d_dlnRho_c_T(i_lnPgas)
                  dfdx = dPgas_dlnRho/(Pgas + Prad)*ln10
               else
                  dPgas_dlnT = Pgas*d_dlnT_c_Rho(i_lnPgas)
                  dPrad_dlnT = 4d0*Prad
                  dfdx = (dPgas_dlnT + dPrad_dlnT)/(Pgas + Prad)*ln10
               end if
            else
               f = res(which_other) - other_value
               if (doing_Rho) then
                  dfdx = d_dlnRho_c_T(which_other)*ln10
               else
                  dfdx = d_dlnT_c_Rho(which_other)*ln10
               end if
            end if
            
         end function f
         
      end subroutine do_safe_get_Rho_T


      end module eosDT_eval
