! ***********************************************************************
!
!   Copyright (C) 2010-2021  The MESA Team
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

         
      integer, parameter :: doing_get_T = 1
      integer, parameter :: doing_get_Pgas = 2
      
      contains


      subroutine Get_eosPT_Results(rq, &
               Z_in, X_in, abar, zbar, &
               species, chem_id, net_iso, xa, &
               aPgas, alogPgas, atemp, alogtemp, &
               Rho, logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas, &
               res, d_dlnRho_c_T, d_dlnT_c_Rho, d_dxa_c_TRho, &
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
         real(dp), intent(inout) :: d_dxa_c_TRho(:,:) ! (nv, species)
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
               call mesa_error(__FILE__,__LINE__,'eosPT')
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
         
         call Get_PT_Results_using_DT( &
            rq, Z, X, abar, zbar, &
            species, chem_id, net_iso, xa, &
            Pgas, logPgas, T, logT, &
            Rho, logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas, &
            res, d_dlnRho_c_T, d_dlnT_c_Rho, d_dxa_c_TRho, ierr)

         if (ierr /= 0) then
            if (dbg) write(*,*) 'error from Get_PT_Results_using_DT'
            return
         end if

      end subroutine Get_eosPT_Results


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
         

      subroutine Get_PT_Results_using_DT( &
               rq, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               Pgas, logPgas, T, logT, &
               Rho, logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas, &
               res, d_dlnRho_c_T, d_dlnT_c_Rho, d_dxa_c_TRho, ierr)
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
         real(dp), intent(inout) :: d_dxa_c_TRho(:,:) ! (nv, species)
         integer, intent(out) :: ierr
         
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
               logRho_result, res, d_dlnRho_c_T, d_dlnT_c_Rho, d_dxa_c_TRho, &
               eos_calls, ierr)
         if (ierr /= 0) then
            if (dbg) then
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
               write(*,'(A)')
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
               res, d_dlnRho_c_T, d_dlnT_c_Rho, d_dxa_c_TRho, &
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
         real(dp), intent(inout) :: d_dxa_c_TRho(:,:) ! (nv, species)
         
         integer, intent(out) :: eos_calls
         integer, intent(out) :: ierr ! 0 means AOK.
         
         call do_safe_get_Pgas_T( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               logPgas, which_other, other_value, doing_get_T, &
               logT_guess, logT_result, logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
               logT_tol, other_tol, max_iter, &
               Rho, logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas, &
               res, d_dlnRho_c_T, d_dlnT_c_Rho, d_dxa_c_TRho, &
               eos_calls, ierr)
      
      end subroutine get_T
      

      subroutine get_Pgas( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               logT, which_other, other_value, &
               logPgas_tol, other_tol, max_iter, logPgas_guess, &
               logPgas_bnd1, logPgas_bnd2, other_at_bnd1, other_at_bnd2, &
               logPgas_result, Rho, logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas, &
               res, d_dlnRho_c_T, d_dlnT_c_Rho, d_dxa_c_TRho, &
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
         real(dp), intent(inout) :: d_dxa_c_TRho(:,:) ! (nv, species)

         integer, intent(out) :: eos_calls
         integer, intent(out) :: ierr ! 0 means AOK.
         
         call do_safe_get_Pgas_T( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               logT, which_other, other_value, doing_get_Pgas, &
               logPgas_guess, logPgas_result, logPgas_bnd1, logPgas_bnd2, other_at_bnd1, other_at_bnd2, &
               logPgas_tol, other_tol, max_iter, &
               Rho, logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas, &
               res, d_dlnRho_c_T, d_dlnT_c_Rho, d_dxa_c_TRho, &
               eos_calls, ierr)

      end subroutine get_Pgas


      subroutine do_safe_get_Pgas_T( &
               handle, Z, XH1, abar, zbar, &
               species, chem_id, net_iso, xa, &
               the_other_log, which_other, other_value, doing_which, &
               initial_guess, x, xbnd1, xbnd2, other_at_bnd1, other_at_bnd2, &
               xacc, yacc, ntry, &
               Rho, logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas, &
               res, d_dlnRho_c_T, d_dlnT_c_Rho, d_dxa_c_TRho, &
               eos_calls, ierr)
         use const_def
         use utils_lib, only: is_bad
         use num_lib, only: brent_safe_zero, look_for_brackets
         use chem_def, only: num_chem_isos
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
         real(dp), intent(inout) :: d_dxa_c_TRho(:,:) ! (nv, species)
         integer, intent(out) :: eos_calls, ierr

         integer :: i, j, lrpar, lipar, max_iter, irho, ix, iz
         real(dp), parameter :: dx = 0.1d0
         integer, pointer :: ipar(:)
         real(dp), pointer :: rpar(:)
         real(dp) :: Pgas, T, xb1, xb3, y1, y3, dfdx, f, logPgas, logT
         type (EoS_General_Info), pointer :: rq

         logical, parameter :: dbg = .false.

         include 'formats'

         ierr = 0

         call get_eos_ptr(handle, rq, ierr)
         if (ierr /= 0) then
            write(*, *) 'get_eos_ptr returned ierr', ierr
            return
         end if

         eos_calls = 0
         x = initial_guess

         if (doing_which /= doing_get_T) then
            Pgas = arg_not_provided
            T = exp10(the_other_log)
         else
            T = arg_not_provided
            Pgas = exp10(the_other_log)
         end if

         lipar = 0; nullify(ipar)
         lrpar = 0; nullify(rpar)

         xb1 = xbnd1; xb3 = xbnd2
         if (xb1 == arg_not_provided .or. xb3 == arg_not_provided .or. xb1 == xb3) then

            if (dbg) then
               write(*,'(A)')
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
               write(*,'(A)')
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
               return
            end if
            !write(*,*) 'done look_for_brackets'
         else
            if (other_at_bnd1 == arg_not_provided) then
               y1 = get_f_df(xb1, dfdx, lrpar, rpar, lipar, ipar, ierr)
               if (ierr /= 0) then
                  return
               end if
            else
               y1 = other_at_bnd1 - other_value
            end if
            if (other_at_bnd2 == arg_not_provided) then
               y3 = get_f_df(xb3, dfdx, lrpar, rpar, lipar, ipar, ierr)
               if (ierr /= 0) then
                  return
               end if
            else
               y3 = other_at_bnd2 - other_value
            end if
         end if

         if (dbg) then
            write(*,'(A)')
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
            return
         end if

         contains

            real(dp) function get_f_df(x, dfdx, lrpar, rpar, lipar, ipar, ierr)
               integer, intent(in) :: lrpar, lipar
               real(dp), intent(in) :: x
               real(dp), intent(out) :: dfdx
               integer, intent(inout), pointer :: ipar(:) ! (lipar)
               real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
               integer, intent(out) :: ierr

               real(dp) :: new, logT, logPgas

               include 'formats'

               ierr = 0
               get_f_df = 0

               dfdx = 0

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
                  res, d_dlnRho_c_T, d_dlnT_c_Rho, d_dxa_c_TRho, &
                  ierr)
               if (ierr /= 0) then
22                format(a30, e26.16)
                  if (dbg) then
                     write(*, *) 'Get_eosPT_Results returned ierr', ierr
                     write(*, 22) 'Z', Z
                     write(*, 22) 'XH1', XH1
                     write(*, 22) 'abar', abar
                     write(*, 22) 'zbar', zbar
                     write(*, 22) 'Pgas', Pgas
                     write(*, 22) 'logPgas', logPgas
                     write(*, 22) 'T', T
                     write(*, 22) 'logT', logT
                     write(*,'(A)')
                  end if
                  return
               end if

               eos_calls = eos_calls+1 ! count eos calls

               new = res(which_other)
               get_f_df = new - other_value

               ! f = f(lnRho(lnPgas,lnT),lnT)
               if (doing_which == doing_get_T) then
                  dfdx = (d_dlnT_c_Rho(which_other) &
                     + dlnRho_dlnT_c_Pgas*d_dlnRho_c_T(which_other))*ln10
               else if (doing_which == doing_get_Pgas) then
                  dfdx = dlnRho_dlnPgas_c_T*d_dlnRho_c_T(which_other)*ln10
               else
                  call mesa_error(__FILE__,__LINE__,'bad value for doing_which in eosPT_eval')
               end if

            end function get_f_df

      end subroutine do_safe_get_Pgas_T

      end module eosPT_eval

