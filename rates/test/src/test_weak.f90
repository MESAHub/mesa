! ***********************************************************************
!
!   Copyright (C) 2011-2020  Bill Paxton & The MESA Team
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

module test_weak
   use rates_lib
   use rates_def
   use math_lib
   use const_def
   use utils_lib, only: mesa_error
   use num_lib, only: dfridr
   
   implicit none
   
   contains
   
   subroutine do_test_weak
      use const_lib
      use chem_lib
      use chem_def, only: iso_name_length, chem_isos
      use rates_def, only: Coulomb_Info
      integer :: ierr, i, ir, nr
      integer, pointer :: ids(:), reaction_ids(:)
      type(Coulomb_Info), pointer :: cc
      real(dp), dimension(:), pointer :: &
         lambda, dlambda_dlnT, dlambda_dlnRho, &
         Q, dQ_dlnT, dQ_dlnRho, &
         Qneu, dQneu_dlnT, dQneu_dlnRho
      real(dp) :: logT, T, T9, dT9, dlnT, YeRho, &
         ye, rho, logRho, dlogRho, eta, d_eta_dlnT, d_eta_dlnRho
      character(len=iso_name_length) :: weak_lhs, weak_rhs
      character(len=2*iso_name_length+1) :: key

      real(dp) :: abar, zbar, z2bar

      real(dp) :: dvardx, dvardx_0, dx_0, err, var_0, xdum
      logical :: doing_d_dlnd
      
      include 'formats'

      ierr = 0
               
      write(*,*) 'check weak_info_list'
      weak_lhs = 'o14'
      weak_rhs = 'n14'
      i = get_weak_info_list_id(weak_lhs, weak_rhs)
      if (i <= 0 .or. i > num_weak_info_list_reactions) then
         write(*,*) 'get_weak_info_list_id failed'
         call mesa_error(__FILE__,__LINE__)
      end if
      call create_weak_dict_key(weak_lhs, weak_rhs, key)
      write(*,1) trim(key)
      write(*,1) 'halflife', weak_info_list_halflife(i)
      write(*,1) 'Qneu', weak_info_list_Qneu(i)
      write(*,*)
      
      d_eta_dlnT = 0
      d_eta_dlnRho = 0
      
      if (.false.) then ! TESTING
         logT =     7.5904236599874348D+00
         logRho =     1.0657946486820271D+00
         ye =     8.2724691280321605D-01
         eta =    -5.3262903257381922D+00
         d_eta_dlnT =    -1.5299344982339016D+00
         d_eta_dlnRho =     9.9482489248846617D-01
      else if (.true.) then ! TESTING
         logT =  log10(9.0d8)
         ye =    0.5d0
         logRho =  log10(4.5d5)
         eta =    10d0

         ! call for cell          565
         ! logT =    8.3534130765231005
         ! logRho =    2.4507395003327828
         ! T =    225638433.79026267
         ! Rho =    282.31860559981624
         ! abar =    4.0424973056746829
         ! zbar =    2.0238390731055702
         ! z2bar =    4.2987387813744071
         ! ye =   0.50064079703023978
         ! eta =   -5.3287260155378711
         ! d_eta_dlnT=   -1.5713400060794886
         ! d_eta_dlnRho =    1.0016532307086357

      else
         logT = 7.5d0
         logRho = 4d0
         Ye = 0.5d0
         eta = 0d0
      end if
      
      T = exp10(logT)
      T9 = T*1d-9
      rho = exp10(logRho)
      YeRho = Ye*rho


      write(*,1) 'logT', logT
      write(*,1) 'logRho', logRho
      write(*,1) 'ye', ye
      write(*,1) 'eta', eta
      write(*,1) 'T9', T9
      write(*,1) 'lYeRho', log10(YeRho)
      write(*,*)
      
      nr = num_weak_reactions
      allocate( &
         ids(nr), reaction_ids(nr), &
         lambda(nr), dlambda_dlnT(nr), dlambda_dlnRho(nr), &
         Q(nr), dQ_dlnT(nr), dQ_dlnRho(nr), &
         Qneu(nr), dQneu_dlnT(nr), dQneu_dlnRho(nr), &
         stat=ierr)
      if (ierr /= 0) return
      do i = 1, nr
         ids(i) = i
         reaction_ids(i) = i
      end do
      
      write(*,*)
      write(*,2) 'nr', nr
      write(*,*)
      
      call eval_weak_reaction_info( &
         nr, ids, reaction_ids, cc, T9, YeRho, &
         eta, d_eta_dlnT, d_eta_dlnRho, &
         lambda, dlambda_dlnT, dlambda_dlnRho, &
         Q, dQ_dlnT, dQ_dlnRho, &
         Qneu, dQneu_dlnT, dQneu_dlnRho, &
         ierr)
      if (ierr /= 0) then
         write(*,*) 'failed in eval_weak_reaction_info'
         call mesa_error(__FILE__,__LINE__)
      end if
      
      if (.true.) then
         write(*,'(30x,99a16)') &
            'halflife', 'Qneu', 'Qtotal'
         do i = 1, nr
            ir = ids(i)
            if (ir <= 0 .or. ir > size(weak_lhs_nuclide_id)) then
               write(*,*) 'ir', ir
               call mesa_error(__FILE__,__LINE__)
            end if
            if (Qneu(i) < 1d-20) cycle
            if (lambda(i) < 1d-20) cycle
            weak_lhs = chem_isos% name(weak_lhs_nuclide_id(ir))
            weak_rhs = chem_isos% name(weak_rhs_nuclide_id(ir))
            write(*,'(a30,99(1pe16.6))') weak_lhs // weak_rhs, &
               ln2/lambda(i), Qneu(i), Q(i)
         end do
         write(*,*)
      else
         write(*,'(30x,5a12,a20)') 'Q', 'Qneu', 'lambda'
         do i = 1, nr
            ir = ids(i)
            if (ir <= 0 .or. ir > size(weak_lhs_nuclide_id)) then
               write(*,*) 'ir', ir
               call mesa_error(__FILE__,__LINE__)
            end if
            weak_lhs = chem_isos% name(weak_lhs_nuclide_id(ir))
            weak_rhs = chem_isos% name(weak_rhs_nuclide_id(ir))
            write(*,'(a30,5f12.6,e20.12)') weak_lhs // weak_rhs, &
               Q(i), Qneu(i), lambda(i)
         end do
         write(*,*)
         
         if (.false.) then
         write(*,'(a30,5a12,a20)') 'd_dT9', 'Q', 'Qneu', 'lambda'
         do i = 1, nr
            ir = ids(i)
            weak_lhs = chem_isos% name(weak_lhs_nuclide_id(ir))
            weak_rhs = chem_isos% name(weak_rhs_nuclide_id(ir))
            write(*,'(a30,5f12.6,e20.12)') weak_lhs // weak_rhs, &
               dQ_dlnT(i), dQneu_dlnT(i), dlambda_dlnT(i)
         end do
         write(*,*)
         end if
      
         if (.false.) then
         write(*,'(a30,5a12,a20)') 'd_d_rho', 'Q', 'Qneu', 'lambda'
         do i = 1, nr
            ir = ids(i)
            weak_lhs = chem_isos% name(weak_lhs_nuclide_id(ir))
            weak_rhs = chem_isos% name(weak_rhs_nuclide_id(ir))
            write(*,'(a30,5f12.6,e20.12)') weak_lhs // weak_rhs, &
               dQ_dlnRho(i), dQneu_dlnRho(i), dlambda_dlnRho(i)
         end do
         write(*,*)
         end if

      end if
      
      write(*,*) 'done'

      if (.false.) then ! dfridr tests for partials

         do i=1, num_weak_reactions

            ir = ids(i)
            weak_lhs = chem_isos% name(weak_lhs_nuclide_id(ir))
            weak_rhs = chem_isos% name(weak_rhs_nuclide_id(ir))
            write(*,'(a30)') weak_lhs // weak_rhs

         doing_d_dlnd = .true.
         doing_d_dlnd = .false.

         var_0 = lambda(i)

         if (doing_d_dlnd) then
            dx_0 = max(1d-14, abs(logRho*ln10*1d-6))
            dvardx_0 = dlambda_dlnRho(i)
         else
            dx_0 = max(1d-14, abs(logT*ln10*1d-6))
            dvardx_0 = dlambda_dlnT(i)
         end if
         err = 0d0
         dvardx = dfridr(dx_0,dfridr_weak_reaction_info,err)
         xdum = (dvardx - dvardx_0)/max(abs(dvardx_0),1d-50)
         write(*,1) 'analytic, numeric, est err in numeric, rel diff', &
               dvardx_0, dvardx, err, xdum
         if (doing_d_dlnd) then
            write(*,*) 'doing dlnd'
         else ! doing d_dlnT
            write(*,*) 'doing dlnT'
         end if
         write(*,*) 'test net'
         write(*,*)

         end do

         call mesa_error(__FILE__,__LINE__,'test rate')

      end if

      
      deallocate( &
         ids, reaction_ids, &
         lambda, dlambda_dlnT, dlambda_dlnRho, &
         Q, dQ_dlnT, dQ_dlnRho, &
         Qneu, dQneu_dlnT, dQneu_dlnRho)


      contains

      real(dp) function dfridr_weak_reaction_info(delta_x) result(val)
         real(dp), intent(in) :: delta_x
         integer :: ierr
         real(dp) :: logYeRho, logT, var, log_var
         include 'formats'
         ierr = 0

         if (doing_d_dlnd) then

            logYeRho = log10(YeRho)
            log_var = logYeRho + delta_x/ln10
            var = exp10(log_var)

            call eval_weak_reaction_info( &
                  nr, ids, reaction_ids, cc, T9, var, &
                  eta, d_eta_dlnT, d_eta_dlnRho, &
                  lambda, dlambda_dlnT, dlambda_dlnRho, &
                  Q, dQ_dlnT, dQ_dlnRho, &
                  Qneu, dQneu_dlnT, dQneu_dlnRho, &
                  ierr)

         else

            logT = log10(T9) + 9d0
            log_var = logT + delta_x/ln10
            var = exp10(log_var - 9d0)

            call eval_weak_reaction_info( &
                  nr, ids, reaction_ids, cc, var, YeRho, &
                  eta, d_eta_dlnT, d_eta_dlnRho, &
                  lambda, dlambda_dlnT, dlambda_dlnRho, &
                  Q, dQ_dlnT, dQ_dlnRho, &
                  Qneu, dQneu_dlnT, dQneu_dlnRho, &
                  ierr)

         end if

         if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'failed in call eval_weak_reaction_info')
         val = lambda(i)

         end function dfridr_weak_reaction_info

   end subroutine do_test_weak


end module test_weak

