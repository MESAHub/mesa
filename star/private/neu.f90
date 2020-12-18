! ***********************************************************************
!
!   Copyright (C) 2013  Bill Paxton
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

      module neu

      use star_private_def
      use const_def
      use utils_lib

      implicit none

      private
      public :: do_neu_for_cell, do_clear_neu_for_cell


      contains




      subroutine do_clear_neu_for_cell(s,k,ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         ierr = 0
         s% non_nuc_neu(k) = 0
         s% d_nonnucneu_dlnd(k) = 0
         s% d_nonnucneu_dlnT(k) = 0
         s% nonnucneu_plas(k) = 0
         s% nonnucneu_brem(k) = 0
         s% nonnucneu_phot(k) = 0
         s% nonnucneu_pair(k) = 0
         s% nonnucneu_reco(k) = 0
      end subroutine do_clear_neu_for_cell


      subroutine do_neu_for_cell(s,k,ierr)
         use neu_def
         use neu_lib
         use chem_def, only: chem_isos
         use const_def,only:ln10
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr

         real(dp) :: loss(num_neu_rvs) ! total from all sources
         real(dp) :: sources(num_neu_types, num_neu_rvs)
         real(dp) :: log10_rho, log10_T
         real(dp), parameter :: log10_Tlim = 7.5d0
         logical :: flags(num_neu_types) ! true if should include the type of loss
         integer :: j

         include 'formats'

         ierr = 0

         if (s% non_nuc_neu_factor <= 0d0) then
            call do_clear_neu_for_cell(s,k,ierr)
            return
         end if

         flags = .true.
         !flags(reco_neu_type) = .false.

         log10_rho = s% lnd(k)/ln10
         log10_T = s% lnT(k)/ln10

         if (s% use_other_neu) then
            call s% other_neu( &
               s% id, k, s% T(k), log10_T, s% rho(k), log10_rho, &
               s% abar(k), s% zbar(k), &
               log10_Tlim, flags, loss, sources, ierr)
         else
            call neu_get( &
               s% T(k), log10_T, s% rho(k), log10_rho, &
               s% abar(k), s% zbar(k), &
               log10_Tlim, flags, loss, sources, ierr)
         end if

         if (ierr /= 0) then
            if (s% report_ierr) then
               write(*,3) 'do_neu_for_cell: neu_get ierr', ierr, k
               write(*,1) 'T=', s% T(k)
               write(*,1) 'log10_T=', log10_T
               write(*,1) 'rho=', s% rho(k)
               write(*,1) 'log10_rho=', log10_rho
               write(*,1) 'abar', s% abar(k)
               write(*,1) 'zbar', s% zbar(k)
               write(*,*)
               return
               stop
            end if
            return
         end if

         if (s% non_nuc_neu_factor /= 1d0) loss(:) = loss(:)*s% non_nuc_neu_factor
         s% non_nuc_neu(k) = loss(ineu)
         s% d_nonnucneu_dlnd(k) = loss(idneu_dRho)*s% rho(k)
         s% d_nonnucneu_dlnT(k) = loss(idneu_dT)*s% T(k)

         s% nonnucneu_plas(k) = sources(plas_neu_type,ineu)
         s% nonnucneu_brem(k) = sources(brem_neu_type,ineu)
         s% nonnucneu_phot(k) = sources(phot_neu_type,ineu)
         s% nonnucneu_pair(k) = sources(pair_neu_type,ineu)
         s% nonnucneu_reco(k) = sources(reco_neu_type,ineu)

         if (is_bad(s% non_nuc_neu(k))) then
            ierr = -1
            if (s% report_ierr) write(*,*) 'do_neu_for_cell ierr for cell', k
            if (s% stop_for_bad_nums) then
               write(*,2) 'bad s% non_nuc_neu(k)', k, s% non_nuc_neu(k)
               stop 'do_neu_for_cell'
            end if
            return
         end if

         if (k == s% trace_k) then
            write(*,5) 'non_nuc_neu', k, &
               s% solver_iter, s% model_number, s% solver_adjust_iter, &
                        s% non_nuc_neu(k)
         end if

      end subroutine do_neu_for_cell



      end module neu

