! ***********************************************************************
!
!   Copyright (C) 2010-2019 Bill Paxton & The MESA Team
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

      module eval_weak
      
      use const_def, only: dp
      use math_lib
      use utils_lib, only: mesa_error
      use rates_def
      use suzuki_tables
      implicit none


      contains
      
      subroutine do_eval_weak_reaction_info( &
            n, ids, reaction_ids, T9, YeRho, &
            eta, d_eta_dlnT, d_eta_dlnRho, &
            lambda, dlambda_dlnT, dlambda_dlnRho, &
            Q, dQ_dlnT, dQ_dlnRho, &
            Qneu, dQneu_dlnT, dQneu_dlnRho, &
            ierr)
         use const_def, only: ln10, kerg, mev_to_ergs
         use chem_def
         use interp_1d_def
         use utils_lib, only: is_bad
         integer, intent(in) :: n, ids(:), reaction_ids(:)
         real(dp), intent(in) :: T9, YeRho, eta, d_eta_dlnT, d_eta_dlnRho
         real(dp), dimension(:), intent(inout), pointer :: &
            lambda, dlambda_dlnT, dlambda_dlnRho, &
            Q, dQ_dlnT, dQ_dlnRho, &
            Qneu, dQneu_dlnT, dQneu_dlnRho
         integer, intent(out) :: ierr
         
         real(dp) :: alfa, beta, d_alfa_dlnT, alfa_hi_Z, beta_hi_Z, d_alfa_hi_Z_dlnT
         integer :: i, ir, cid

         include 'formats'
         
         call do_eval_weaklib_reaction_info( &
            n, ids, T9, YeRho, &
            eta, d_eta_dlnT, d_eta_dlnRho, &
            lambda, dlambda_dlnT, dlambda_dlnRho, &
            Q, dQ_dlnT, dQ_dlnRho, &
            Qneu, dQneu_dlnT, dQneu_dlnRho, &
            ierr)
         if (ierr /= 0) then
            !write(*,*) 'failed in do_eval_weaklib_reaction_info'
            return
         end if

         !if (T9 >= max(T9_weaklib_full_on, T9_weaklib_full_on_hi_Z)) return
         
         ! revise lambda using rate for low T         
         ! alfa is fraction from weaklib
         if (T9 >= T9_weaklib_full_on) then
            alfa = 1d0
            d_alfa_dlnT = 0d0
         else if (T9 > T9_weaklib_full_off) then
            alfa = (T9 - T9_weaklib_full_off)/ &
                 (T9_weaklib_full_on - T9_weaklib_full_off)
            d_alfa_dlnT = T9 / &
                 (T9_weaklib_full_on - T9_weaklib_full_off)
         else
            alfa = 0d0
            d_alfa_dlnT = 0d0
         end if
         beta = 1d0 - alfa ! beta is fraction for low T

         if (T9 >= T9_weaklib_full_on_hi_Z) then
            alfa_hi_Z = 1d0
            d_alfa_hi_Z_dlnT = 0d0
         else if (T9 > T9_weaklib_full_off_hi_Z) then
            alfa_hi_Z = (T9 - T9_weaklib_full_off_hi_Z)/ &
                 (T9_weaklib_full_on_hi_Z - T9_weaklib_full_off_hi_Z)
            d_alfa_hi_Z_dlnT = T9 / &
                 (T9_weaklib_full_on_hi_Z - T9_weaklib_full_off_hi_Z)
         else
            alfa_hi_Z = 0d0
            d_alfa_hi_Z_dlnT = 0d0
         end if
         beta_hi_Z = 1d0 - alfa_hi_Z

         do i = 1, n
            ir = reaction_ids(i)
            if (ir == 0) cycle
            cid = weak_reaction_info(1,ir)
            if (weak_lowT_rate(ir) <= 0d0) cycle
            if (cid <= 0) cycle  
            if (ids(i) <= 0) then
               lambda(i) = weak_lowT_rate(ir)
               dlambda_dlnT(i) = 0d0
               dlambda_dlnRho(i) = 0d0
            else if (chem_isos% Z(cid) >= weaklib_blend_hi_Z) then
               lambda(i) = alfa_hi_Z*lambda(i) + beta_hi_Z*weak_lowT_rate(ir)
               dlambda_dlnT(i) = alfa_hi_Z*dlambda_dlnT(i) + &
                    d_alfa_hi_Z_dlnT * (lambda(i) - weak_lowT_rate(ir))
               dlambda_dlnRho(i) = alfa_hi_Z*dlambda_dlnRho(i)
            else
               lambda(i) = alfa*lambda(i) + beta*weak_lowT_rate(ir)
               dlambda_dlnT(i) = alfa*dlambda_dlnT(i) + &
                    d_alfa_dlnT * (lambda(i) - weak_lowT_rate(ir))
               dlambda_dlnRho(i) = alfa*dlambda_dlnRho(i)
            end if
         end do
         
      end subroutine do_eval_weak_reaction_info

      
      subroutine do_eval_weaklib_reaction_info( &
            n, ids, T9_in, YeRho_in, &
            eta, d_eta_dlnT, d_eta_dlnRho, &
            lambda, dlambda_dlnT, dlambda_dlnRho, &
            Q, dQ_dlnT, dQ_dlnRho, &
            Qneu, dQneu_dlnT, dQneu_dlnRho, &
            ierr)
         use const_def, only: ln10, kerg, mev_to_ergs
         use chem_def
         use interp_1d_def
         use utils_lib, only: is_bad, integer_dict_lookup
         integer, intent(in) :: n, ids(:)
         real(dp), intent(in) :: T9_in, YeRho_in, eta, d_eta_dlnT, d_eta_dlnRho
         real(dp), dimension(:), intent(inout), pointer :: &
            lambda, dlambda_dlnT, dlambda_dlnRho, &
            Q, dQ_dlnT, dQ_dlnRho, &
            Qneu, dQneu_dlnT, dQneu_dlnRho
         integer, intent(out) :: ierr
         
         logical, parameter :: dbg = .false.
         
         real(dp) :: T, T9, YeRho, lYeRho
         integer :: i, ir, in, out, rxn_idx
         logical :: neg
         real(dp) :: decay, capture, Qx, Qn, conv, mue, d_mue_dlnRho, d_mue_dlnT
         character(len=iso_name_length) :: weak_lhs, weak_rhs
         real(dp) :: delta_Q, Vs
         integer, parameter :: nwork = pm_work_size

         real(dp) :: &
              s_lambda, s_dlambda_dlnT, s_dlambda_dlnRho, &
              s_Qneu, s_dQneu_dlnT, s_dQneu_dlnRho,&
              s_delta_Q, s_Vs


         character(len=2*iso_name_length+1) :: key

         class(weak_rate_table), pointer :: table

         include 'formats'

         ierr = 0
         
         T9 = T9_in
         YeRho = YeRho_in
         lYeRho = log10(YeRho_in)
         if (is_bad(lYeRho)) then
            ierr = -1
            return
            
            write(*,1) 'lYeRho', lYeRho
            write(*,1) 'YeRho_in', YeRho_in
            write(*,1) 'log10(YeRho_in)', log10(YeRho_in)
            !stop 'weak lYeRho'
         end if
         
         if (n == 0) then
            write(*,*) 'problem in eval_weak_reaction_info: n == 0'
            write(*,2) 'n', n
            write(*,1) 'T9', T9
           return
         end if
         
       do i = 1, n

            lambda(i) = 0d0
            dlambda_dlnT(i) = 0d0
            dlambda_dlnRho(i) = 0d0
            Q(i) = 0d0
            dQ_dlnT(i) = 0d0
            dQ_dlnRho(i) = 0d0
            Qneu(i) = 0d0
            dQneu_dlnT(i) = 0d0
            dQneu_dlnRho(i) = 0d0

            ir = ids(i)
            if (ir <= 0) cycle
            if (ir > size(weak_reactions_tables)) then
               call show_stuff
               write(*,*) 'bad ir', ir, i
               call mesa_error(__FILE__,__LINE__)
               ierr = -1
               return
            end if

            table => weak_reactions_tables(ir) % t

            T9 = T9_in
            YeRho = YeRho_in
            lYeRho = log10(YeRho_in)

            ! convert to MeV
            conv = kerg/mev_to_ergs
            T = T9*1d9
            mue = eta*conv*T
            d_mue_dlnRho = d_eta_dlnRho*conv*T
            if (d_eta_dlnT == 0) then
               d_mue_dlnT = 0
            else
               d_mue_dlnT = d_eta_dlnT*conv*T + mue
            end if

         ! clip small values to edge of table
         if (T9 < table % T9s(1)) &
            T9 = table % T9s(1)
         if (lYeRho < table % lYeRhos(1)) &
            lYeRho = table % lYeRhos(1)

         ! clip large values to edge of table
         if (T9 > table % T9s(table % num_T9)) &
            T9 = table % T9s(table % num_T9)
         if (lYeRho > table % lYeRhos(table % num_lYeRho)) &
            lYeRho = table % lYeRhos(table % num_lYeRho)

            call table% interpolate(T9, lYeRho, &
               lambda(i), dlambda_dlnT(i), dlambda_dlnRho(i), & 
               Qneu(i), dQneu_dlnT(i), dQneu_dlnRho(i), &
               delta_Q, Vs, ierr)

            in = weak_lhs_nuclide_id(ir)
            out = weak_rhs_nuclide_id(ir)
            Qx = chem_isos% mass_excess(in) - chem_isos% mass_excess(out)

            if (use_suzuki_tables) then
               ! now, if there's a suzuki reaction, use that one instead
               call create_weak_dict_key(weak_lhs_nuclide_name(ir), weak_rhs_nuclide_name(ir), key)
               call integer_dict_lookup(suzuki_reactions_dict, key, rxn_idx, ierr)
               if (ierr /=0) then
                  if (dbg) write(*,*) key, "is not a reaction included in the Suzuki tables"
                  ierr = 0
               else
                  if (dbg) write(*,*) key, "is a reaction included in the Suzuki tables"
                  table => suzuki_reactions_tables(rxn_idx) %t
                  call table % interpolate(T9, lYeRho, &
                       s_lambda, s_dlambda_dlnT, s_dlambda_dlnRho, &
                       s_Qneu, s_dQneu_dlnT, s_dQneu_dlnRho, &
                       s_delta_Q, s_Vs, ierr)
                  if (ierr == 0) then
                     lambda(i) = s_lambda
                     dlambda_dlnT(i) = s_dlambda_dlnT
                     dlambda_dlnRho(i) = s_dlambda_dlnRho
                     Qneu(i) = s_Qneu
                     dQneu_dlnT(i) = s_dQneu_dlnT
                     dQneu_dlnRho(i) = s_dQneu_dlnRho
                     delta_Q = s_delta_Q
                     Vs = s_Vs
                     !write(*,*) lYeRho, T9, s_lambda, s_Qneu
                  end if
                  ierr = 0
               end if
            end if

          ! neg is true for electron capture and positron emission
          neg = ((chem_isos% Z(in) - chem_isos% Z(out)) == 1.0d0)

          ! in the past, these Q values used to include terms
          ! associated with the electron and ion chemical potentials
          ! these terms are now handled elsewhere, so Q is just the change in rest mass.
          ! since Qx is made from atomic mass excesses, it includes the electron rest mass.

          if (neg) then ! electron capture and positron emission
             Q(i) = Qx
             dQ_dlnT(i) = 0
             dQ_dlnRho(i) = 0
          else ! positron capture and electron emission
             Q(i) = Qx
             dQ_dlnT(i) = 0
             dQ_dlnRho(i) = 0
          end if

          if (lambda(i) < 1d-30) then
             Qneu(i) = 0
             dQneu_dlnT(i) = 0
             dQneu_dlnRho(i) = 0
          end if

            if (is_bad(Qneu(i))) then
               ierr = -1
               return

               write(*,2) 'lambda', i, lambda(i)
               write(*,2) 'Qneu', i, Qneu(i)
               call show_stuff
            end if
            
         end do
                  
         if (is_bad(lYeRho)) then
            ierr = -1
            return
            call show_stuff
         end if
      
      
         contains
         
         subroutine show_stuff
            include 'formats'
            write(*,1) 'T9', T9
            write(*,1) 'lYeRho', lYeRho
            write(*,1) 'eta', eta
            write(*,1) 'mue', mue
            write(*,*)
            do i = 1, n
               ir = ids(i)
               in = weak_lhs_nuclide_id(i)
               out = weak_rhs_nuclide_id(i)
               if (.true.) then
                  write(*,4) 'ir, i, n', ir, i, n
                  if (ir <= 0 .or. ir > ubound(weak_lhs_nuclide_id,dim=1)) cycle
                  weak_lhs = chem_isos% name(weak_lhs_nuclide_id(ir))
                  weak_rhs = chem_isos% name(weak_rhs_nuclide_id(ir))
                  write(*,'(a30,3i5)') weak_lhs // weak_rhs, ir, i, n
                  !write(*,1) 'chem_isos% mass_excess(in)', chem_isos% mass_excess(in)
                  !write(*,1) 'chem_isos% mass_excess(out)', chem_isos% mass_excess(out)
                  write(*,2) 'Qx', i, chem_isos% mass_excess(in) - chem_isos% mass_excess(out)
                  write(*,2) 'Q', i, Q(i)
                  write(*,2) 'Qneu', i, Qneu(i)
                  write(*,*)
               end if
            end do
            call mesa_error(__FILE__,__LINE__)
         end subroutine show_stuff         
         
      end subroutine do_eval_weaklib_reaction_info
      



      end module eval_weak

