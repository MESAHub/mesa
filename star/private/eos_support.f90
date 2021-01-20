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

module eos_support

  ! Uses

  use const_def
  use star_private_def
  use utils_lib, only : is_bad

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: MAX_ITER_FOR_SOLVE = 100

  ! Access specifiers

  private

  public :: get_eos
  public :: solve_eos_given_DE
  public :: solve_eos_given_DEgas
  public :: solve_eos_given_DP
  public :: solve_eos_given_DS
  public :: solve_eos_given_PT
  public :: solve_eos_given_PgasT
  public :: solve_eos_given_PgasT_auto
  
  ! Procedures

contains

  ! Get eos results data given density & temperature

  subroutine get_eos( &
       s, k, xa, &
       Rho, logRho, T, logT, &
       res, dres_dlnRho, dres_dlnT, &
       dres_dxa, ierr)

    use eos_lib, only: eosDT_get
    use eos_def, only: num_eos_basic_results, num_eos_d_dxa_results, num_helm_results, i_lnE

    type (star_info), pointer :: s
    integer, intent(in) :: k ! 0 means not being called for a particular cell
    real(dp), intent(in) :: xa(:), Rho, logRho, T, logT
    real(dp), dimension(num_eos_basic_results), intent(out) :: &
         res, dres_dlnRho, dres_dlnT
    real(dp), intent(out) :: dres_dxa(num_eos_d_dxa_results,s% species)
    integer, intent(out) :: ierr

    real(dp), dimension(num_eos_basic_results) :: dres_dabar, dres_dzbar
    integer :: j
    logical :: off_table

    include 'formats'

    ierr = 0

    if (s% doing_timing) &
       s% timing_num_get_eos_calls = s% timing_num_get_eos_calls + 1

    if (s% use_other_eos) then
       call s% other_eosDT_get( &
            s% id, k, &
            s% eos_handle, s% species, s% chem_id, s% net_iso, xa, &
            Rho, logRho, T, logT, &
            res, dres_dlnRho, dres_dlnT, dres_dxa, ierr)
    else
       call eosDT_get( &
            s% eos_handle, s% species, s% chem_id, s% net_iso, xa, &
            Rho, logRho, T, logT, &
            res, dres_dlnRho, dres_dlnT, dres_dxa, ierr)
    end if

    if (ierr /= 0) then
       s% retry_message = 'get_eos failed'
       if (s% report_ierr) then
          !$OMP critical (get_eos_critical)
          write(*,*) 'get_eos ierr', ierr
          write(*,2) 'k', k
          do j=1,s% species
             write(*,2) 'xa(j) ' // trim(s% nameofequ(j+s% nvar_hydro)), j, xa(j)
          end do
          write(*,1) 'log10Rho', logRho
          write(*,1) 'log10T', logT
          if (s% stop_for_bad_nums .and. &
               is_bad(logRho+logT)) stop 'do_eos_for_cell'
          !$OMP end critical (get_eos_critical)
       end if
       return
    end if

  end subroutine get_eos

  !****

  ! Solve for temperature & eos results data given density & energy

  subroutine solve_eos_given_DE( &
       s, k, z, xh, abar, zbar, xa, &
       logRho, logE, logT_guess, logT_tol, logE_tol, &
       logT, res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, &
       ierr)

    use eos_def
    use eos_lib, only: eosDT_get_T

    type (star_info), pointer :: s
    integer, intent(in) :: k ! 0 indicates not for a particular cell.
    real(dp), intent(in) :: &
         z, xh, abar, zbar, xa(:), logRho, logE, &
         logT_guess, logT_tol, logE_tol
    real(dp), intent(out) :: logT
    real(dp), dimension(num_eos_basic_results), intent(out) :: &
         res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar
    integer, intent(out) :: ierr

    integer :: eos_calls
    real(dp) :: eos_x, eos_z

    include 'formats'

    ierr = 0

    if (s% doing_timing) s% timing_num_solve_eos_calls = s% timing_num_solve_eos_calls + 1

    if (s% use_fixed_XZ_for_eos) then
       eos_x = s% fixed_X_for_eos
       eos_z = s% fixed_Z_for_eos
    else
       eos_x = xh
       eos_z = z
    end if

    if (s% use_other_eos) then
       call s% other_eosDT_get_T( &
            s% id, k, s% eos_handle, eos_z, eos_x, abar, zbar, &
            s% species, s% chem_id, s% net_iso, xa, &
            logRho, i_lnE, logE*ln10, &
            logT_tol, logE_tol*ln10, MAX_ITER_FOR_SOLVE, logT_guess, &
            arg_not_provided, arg_not_provided, arg_not_provided, arg_not_provided, &
            logT, res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, &
            eos_calls, ierr)
    else
       call eosDT_get_T( &
            s% eos_handle, eos_z, eos_x, abar, zbar, &
            s% species, s% chem_id, s% net_iso, xa, &
            logRho, i_lnE, logE*ln10, &
            logT_tol, logE_tol*ln10, MAX_ITER_FOR_SOLVE, logT_guess,  &
            arg_not_provided, arg_not_provided, arg_not_provided, arg_not_provided, &
            logT, res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, &
            eos_calls, ierr)
    end if

  end subroutine solve_eos_given_DE
  
  !****

  ! Solve for temperature & eos results data given density & gas energy

  subroutine solve_eos_given_DEgas( &
       s, k, z, xh, abar, zbar, xa, &
       logRho, egas, logT_guess, logT_tol, egas_tol, &
       logT, res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, &
       ierr)

    use eos_def
    use eos_lib, only: eosDT_get_T_given_egas

    type (star_info), pointer :: s
    integer, intent(in) :: k ! 0 indicates not for a particular cell.
    real(dp), intent(in) :: &
         z, xh, abar, zbar, xa(:), logRho, egas, &
         logT_guess, logT_tol, egas_tol
    real(dp), intent(out) :: logT
    real(dp), dimension(num_eos_basic_results), intent(out) :: &
         res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar
    integer, intent(out) :: ierr

    integer :: eos_calls
    real(dp) :: eos_x, eos_z

    include 'formats'

    ierr = 0

    if (s% doing_timing) s% timing_num_solve_eos_calls = s% timing_num_solve_eos_calls + 1

    if (s% use_fixed_XZ_for_eos) then
       eos_x = s% fixed_X_for_eos
       eos_z = s% fixed_Z_for_eos
    else
       eos_x = xh
       eos_z = z
    end if

    if (s% use_other_eos) then
       write(*,*) 'cannot call solve_eos_given_DEgas with use_other_eos'
       ierr = -1
       return
    else
      call eosDT_get_T_given_egas( &
         s% eos_handle, eos_z, eos_x, abar, zbar, &
         s% species, s% chem_id, s% net_iso, xa, &            
         logRho, egas, logT_tol, egas_tol, MAX_ITER_FOR_SOLVE, logT_guess, &
         arg_not_provided, arg_not_provided, arg_not_provided, arg_not_provided, &
         logT, res, dres_dlnRho, dres_dlnT, &
         dres_dabar, dres_dzbar, eos_calls, ierr)
    end if

  end subroutine solve_eos_given_DEgas

  !****

  ! Solve for temperature & eos results data given density & pressure

  subroutine solve_eos_given_DP( &
       s, k, z, xh, abar, zbar, xa, &
       logRho, logP, logT_guess, logT_tol, logP_tol, &
       logT, res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, &
       ierr)

    use eos_def
    use eos_lib, only: eosDT_get_T_given_Ptotal

    type (star_info), pointer :: s
    integer, intent(in) :: k ! 0 indicates not for a particular cell.
    real(dp), intent(in) :: &
         z, xh, abar, zbar, xa(:), logRho, logP, &
         logT_guess, logT_tol, logP_tol
    real(dp), intent(out) :: logT
    real(dp), dimension(num_eos_basic_results), intent(out) :: &
         res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar
    integer, intent(out) :: ierr

    integer :: eos_calls
    real(dp) :: eos_x, eos_z

    include 'formats'

    ierr = 0

    if (s% doing_timing) s% timing_num_solve_eos_calls = s% timing_num_solve_eos_calls + 1

    if (s% use_fixed_XZ_for_eos) then
       eos_x = s% fixed_X_for_eos
       eos_z = s% fixed_Z_for_eos
    else
       eos_x = xh
       eos_z = z
    end if

    if (s% use_other_eos) then
       write(*,*) 'cannot call solve_eos_given_DP with use_other_eos set'
       ierr = -1
       return
    else
       call eosDT_get_T_given_Ptotal( &
            s% eos_handle, eos_z, eos_x, abar, zbar, &
            s% species, s% chem_id, s% net_iso, xa, &
            logRho, logP, logT_tol, logP_tol, MAX_ITER_FOR_SOLVE, logT_guess, &
            arg_not_provided, arg_not_provided, arg_not_provided, arg_not_provided, &
            logT, res, dres_dlnRho, dres_dlnT, &
            dres_dabar, dres_dzbar, eos_calls, ierr)
    end if
          
  end subroutine solve_eos_given_DP

  !****

  ! Solve for temperature & eos results data for a given density &
  ! entropy

  subroutine solve_eos_given_DS( &
       s, k, z, xh, abar, zbar, xa, &
       logRho, logS, logT_guess, logT_tol, logS_tol, &
       logT, res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, &
       ierr)

    use eos_def
    use eos_lib, only: eosDT_get_T

    type (star_info), pointer :: s
    integer, intent(in) :: k ! 0 indicates not for a particular cell.
    real(dp), intent(in) :: &
         z, xh, abar, zbar, xa(:), logRho, logS, &
         logT_guess, logT_tol, logS_tol
    real(dp), intent(out) :: logT
    real(dp), dimension(num_eos_basic_results), intent(out) :: &
         res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar
    integer, intent(out) :: ierr

    integer :: eos_calls
    real(dp) :: eos_x, eos_z

    include 'formats'

    ierr = 0
    
    if (s% use_fixed_XZ_for_eos) then
       eos_x = s% fixed_X_for_eos
       eos_z = s% fixed_Z_for_eos
    else
       eos_x = xh
       eos_z = z
    end if

    if (s% use_other_eos) then
       call s% other_eosDT_get_T( &
            s% id, k, s% eos_handle, eos_z, eos_x, abar, zbar, &
            s% species, s% chem_id, s% net_iso, xa, &
            logRho, i_lnS, logS*ln10, &
            logT_tol, logS_tol*ln10, MAX_ITER_FOR_SOLVE, logT_guess, &
            arg_not_provided, arg_not_provided, arg_not_provided, arg_not_provided, &
            logT, res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, &
            eos_calls, ierr)
    else
       call eosDT_get_T( &
            s% eos_handle, eos_z, eos_x, abar, zbar, &
            s% species, s% chem_id, s% net_iso, xa, &
            logRho, i_lnS, logS*ln10, &
            logT_tol, logS_tol*ln10, MAX_ITER_FOR_SOLVE, logT_guess,  &
            arg_not_provided, arg_not_provided, arg_not_provided, arg_not_provided, &
            logT, res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, &
            eos_calls, ierr)
    end if

  end subroutine solve_eos_given_DS

  !****

  ! Solve for density & eos results data given pressure & temperature

  subroutine solve_eos_given_PT( &
       s, k, z, xh, abar, zbar, xa, &
       logT, logP, logRho_guess, logRho_tol, logP_tol, &
       logRho, res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, &
       ierr)

    use eos_def
    use eos_lib, only: eosDT_get_Rho_given_Ptotal

    type (star_info), pointer :: s
    integer, intent(in) :: k ! 0 indicates not for a particular cell.
    real(dp), intent(in) :: &
         z, xh, abar, zbar, xa(:), logT, logP, &
         logRho_guess, logRho_tol, logP_tol
    real(dp), intent(out) :: logRho
    real(dp), dimension(num_eos_basic_results), intent(out) :: &
         res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar
    integer, intent(out) :: ierr

    integer :: eos_calls
    real(dp) :: eos_x, eos_z

    include 'formats'

    ierr = 0

    if (s% doing_timing) s% timing_num_solve_eos_calls = s% timing_num_solve_eos_calls + 1

    if (s% use_fixed_XZ_for_eos) then
       eos_x = s% fixed_X_for_eos
       eos_z = s% fixed_Z_for_eos
    else
       eos_x = xh
       eos_z = z
    end if

    if (s% use_other_eos) then
       write(*,*) 'cannot call solve_eos_given_PT with use_other_eos set'
       ierr = -1
       return
    else
       call eosDT_get_Rho_given_Ptotal( &
            s% eos_handle, eos_z, eos_x, abar, zbar, &
            s% species, s% chem_id, s% net_iso, xa, &
            logT, logP, logRho_tol, logP_tol, MAX_ITER_FOR_SOLVE, logRho_guess, &
            arg_not_provided, arg_not_provided, arg_not_provided, arg_not_provided, &
            logRho, res, dres_dlnRho, dres_dlnT, &
            dres_dabar, dres_dzbar, eos_calls, ierr)
    end if

  end subroutine solve_eos_given_PT

  !****

  ! Solve for density & eos results data given gas pressure &
  ! temperature

  subroutine solve_eos_given_PgasT( &
       s, k, z, xh, abar, zbar, xa, &
       logT, logPgas, logRho_guess, logRho_tol, logPgas_tol, &
       logRho, res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, &
       ierr)

    use eos_def
    use eos_lib, only: eosDT_get_Rho

    type (star_info), pointer :: s
    integer, intent(in) :: k ! 0 indicates not for a particular cell.
    real(dp), intent(in) :: &
         z, xh, abar, zbar, xa(:), logT, logPgas, &
         logRho_guess, logRho_tol, logPgas_tol
    real(dp), intent(out) :: logRho
    real(dp), dimension(num_eos_basic_results), intent(out) :: &
         res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar
    integer, intent(out) :: ierr

    integer :: eos_calls
    real(dp) :: eos_x, eos_z

    include 'formats'

    ierr = 0

    if (s% doing_timing) s% timing_num_solve_eos_calls = s% timing_num_solve_eos_calls + 1

    if (s% use_fixed_XZ_for_eos) then
       eos_x = s% fixed_X_for_eos
       eos_z = s% fixed_Z_for_eos
    else
       eos_x = xh
       eos_z = z
    end if

    if (s% use_other_eos) then
       call s% other_eosDT_get_Rho( &
            s% id, k, s% eos_handle, eos_z, eos_x, abar, zbar, &
            s% species, s% chem_id, s% net_iso, xa, &
            logT, i_lnPgas, logPgas*ln10, &
            logRho_tol, logPgas_tol*ln10, MAX_ITER_FOR_SOLVE, logRho_guess, &
            arg_not_provided, arg_not_provided, arg_not_provided, arg_not_provided, &
            logRho, res, dres_dlnRho, dres_dlnT, &
            dres_dabar, dres_dzbar, eos_calls, ierr)
    else if (is_bad(logPgas)) then
       write(*,*) 'bad logPgas for solve_eos_given_PgasT', logPgas
       ierr = -1
       return
    else
       call eosDT_get_Rho( &
            s% eos_handle, eos_z, eos_x, abar, zbar, &
            s% species, s% chem_id, s% net_iso, xa, &
            logT, i_lnPgas, logPgas*ln10, &
            logRho_tol, logPgas_tol*ln10, MAX_ITER_FOR_SOLVE, logRho_guess, &
            arg_not_provided, arg_not_provided, arg_not_provided, arg_not_provided, &
            logRho, res, dres_dlnRho, dres_dlnT, &
            dres_dabar, dres_dzbar, eos_calls, ierr)
       if (ierr /= 0 .and. s% report_ierr) then
          write(*,*) 'Call to eosDT_get_Rho failed in solve_eos_given_PgasT'
          write(*,2) 'logPgas', k, logPgas
          write(*,2) 'logT', k, logT
          write(*,2) 'logRho_guess', k, logRho_guess
       end if
    end if

  end subroutine solve_eos_given_PgasT

  !****

  ! Solve for density & eos results data given gas pressure &
  ! temperature, with logRho_guess calculated automatically via an
  ! initial call to eos_gamma_PT_get

  subroutine solve_eos_given_PgasT_auto( &
       s, k, z, xh, abar, zbar, xa, &
       logT, logPgas, logRho_tol, logPgas_tol, &
       logRho, res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, &
       ierr)

    use eos_def
    use eos_lib, only: eos_gamma_PT_get

    type (star_info), pointer :: s
    integer, intent(in) :: k ! 0 indicates not for a particular cell.
    real(dp), intent(in) :: &
         z, xh, abar, zbar, xa(:), logT, logPgas, &
         logRho_tol, logPgas_tol
    real(dp), intent(out) :: logRho
    real(dp), dimension(num_eos_basic_results), intent(out) :: &
         res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar
    integer, intent(out) :: ierr

    real(dp) :: dres_dxa(num_eos_d_dxa_results,s% species)
    
    real(dp) :: rho_guess, logRho_guess, gamma, &
         dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas
    
    gamma = 5d0/3d0
    call eos_gamma_PT_get( &
       s% eos_handle, abar, exp10(logPgas), logPgas, exp10(logT), logT, gamma, &
       rho_guess, logRho_guess, res, dres_dlnRho, dres_dlnT, &
       ierr)
    if (ierr /= 0) then
       ierr = 0
       logRho_guess = arg_not_provided
    end if

    call solve_eos_given_PgasT( &
       s, k, z, xh, abar, zbar, xa, &
       logT, logPgas, logRho_guess, logRho_tol, logPgas_tol, &
       logRho, res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, &
       ierr)

  end subroutine solve_eos_given_PgasT_auto
         
  !****

end module eos_support
