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

module atm_support

  ! Uses

  use star_private_def
  use const_def
  use utils_lib, only: mesa_error, is_bad

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: get_atm_PT
  public :: get_atm_PT_legacy_grey_and_kap
  public :: get_atm_tau_base
  public :: get_T_tau_id
  public :: build_atm

  ! Procedures

contains

  subroutine get_atm_PT( &
       s, tau_surf, L, R, M, cgrav, skip_partials, &
       Teff, &
       lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
       lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
       ierr)

    type (star_info), pointer :: s
    real(dp), intent(in)      :: tau_surf, L, R, M, cgrav
    logical, intent(in)       :: skip_partials
    real(dp), intent(in)      :: Teff
    real(dp), intent(out)     :: lnT_surf
    real(dp), intent(out)     :: dlnT_dL
    real(dp), intent(out)     :: dlnT_dlnR
    real(dp), intent(out)     :: dlnT_dlnM
    real(dp), intent(out)     :: dlnT_dlnkap
    real(dp), intent(out)     :: lnP_surf
    real(dp), intent(out)     :: dlnP_dL
    real(dp), intent(out)     :: dlnP_dlnR
    real(dp), intent(out)     :: dlnP_dlnM
    real(dp), intent(out)     :: dlnP_dlnkap
    integer, intent(out)      :: ierr

    real(dp) :: kap
    real(dp) :: tau
    real(dp) :: r_surf
    real(dp) :: L_surf
    
    if (is_bad(tau_surf)) then
       write(*,*) 'tau_surf', tau_surf
       ierr = -1
       return
       call mesa_error(__FILE__,__LINE__,'bad tau_surf arg for get_atm_PT')
    end if

    ! Evaluate surface temperature and pressure by dispatching to the
    ! appropriate internal routine

    select case (s% atm_option)

    case ('T_tau')

       call get_T_tau( &
            s, tau_surf, L, R, M, cgrav, &
            s% atm_T_tau_relation, s% atm_T_tau_opacity, skip_partials, &
            Teff, kap, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
            ierr)

    case ('table')

       call get_table( &
            s, skip_partials, L, R, M, cgrav, &
            Teff, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
            ierr)

    case ('irradiated_grey')
    
       call get_irradiated( &
            s, s% atm_irradiated_opacity, skip_partials, L, R, M, cgrav, &
            Teff, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
            ierr)

    case default

       call get_legacy(s, ierr)

    end select

    ! Finish

  end subroutine get_atm_PT

  !****

  subroutine get_atm_PT_legacy_grey_and_kap( &
       s, tau_surf, L, R, M, cgrav, skip_partials, &
       Teff, kap, &
       lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
       lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
       ierr)

    type(star_info), pointer :: s
    real(dp), intent(in)     :: tau_surf, L, R, M, cgrav
    logical, intent(in)      :: skip_partials
    real(dp), intent(in)     :: Teff
    real(dp), intent(out)    :: kap
    real(dp), intent(out)    :: lnT_surf
    real(dp), intent(out)    :: dlnT_dL
    real(dp), intent(out)    :: dlnT_dlnR
    real(dp), intent(out)    :: dlnT_dlnM
    real(dp), intent(out)    :: dlnT_dlnkap
    real(dp), intent(out)    :: lnP_surf
    real(dp), intent(out)    :: dlnP_dL
    real(dp), intent(out)    :: dlnP_dlnR
    real(dp), intent(out)    :: dlnP_dlnM
    real(dp), intent(out)    :: dlnP_dlnkap
    integer, intent(out)     :: ierr

    ! This routine is for legacy callers who require the old
    ! grey_and_kap behavior

    call get_T_tau( &
       s, tau_surf, L, R, M, cgrav, 'Eddington', 'iterated', skip_partials, &
       Teff, kap, &
       lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
       lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
       ierr)

    ! Finish

    return

  end subroutine get_atm_PT_legacy_grey_and_kap

  !****

  subroutine get_atm_tau_base(s, tau_base, ierr)

    use atm_lib, only: &
         atm_get_T_tau_base, &
         atm_get_table_base

    type(star_info), pointer :: s
    real(dp), intent(out)    :: tau_base
    integer, intent(out)     :: ierr

    integer :: T_tau_id
    integer :: table_id

    ierr = 0

    ! Get the base optical depth

    select case (s% atm_option)

    case ('T_tau')

       call get_T_tau_id(s% atm_T_tau_relation, T_tau_id, ierr)
       if (ierr /= 0) then
         s% retry_message = 'Call to get_T_tau_id failed in get_atm_tau_base'
         if (s% report_ierr) write(*, *) s% retry_message
          return
       endif

       call atm_get_T_tau_base(T_tau_id, tau_base, ierr)
       if (ierr /= 0) then
         s% retry_message = 'Call to atm_get_T_tau_base failed in get_atm_tau_base'
         if (s% report_ierr) write(*, *) s% retry_message
          return
       end if

    case ('table')

       call get_table_id(s% atm_table, table_id, ierr)
       if (ierr /= 0) then
         s% retry_message = 'Call to get_table_id failed in get_atm_tau_base'
         if (s% report_ierr) write(*, *) s% retry_message
          return
       endif

       call atm_get_table_base(table_id, tau_base, ierr)
       if (ierr /= 0) then
         s% retry_message = 'Call to atm_get_table_base failed in get_atm_tau_base'
         if (s% report_ierr) write(*, *) s% retry_message
          return
       end if

    case ('irradiated_grey')

       tau_base = 2._dp/3._dp

    case default

       ! All other options use this

       tau_base = 2._dp/3._dp
       
    end select

    ! Finish

    return

  end subroutine get_atm_tau_base

  !****

  subroutine get_T_tau( &
       s, tau_surf, L, R, M, cgrav, T_tau_relation, T_tau_opacity, skip_partials, &
       Teff, kap, &
       lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
       lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
       ierr)

    use atm_lib, only: &
         atm_eval_T_tau_uniform, &
         atm_eval_T_tau_varying
    use kap_support, only : prepare_kap

    type (star_info), pointer :: s
    real(dp), intent(in)      :: tau_surf, L, R, M, cgrav
    character(*), intent(in)  :: T_tau_relation
    character(*), intent(in)  :: T_tau_opacity
    logical, intent(in)       :: skip_partials
    real(dp), intent(in)      :: Teff
    real(dp), intent(out)     :: kap
    real(dp), intent(out)     :: lnT_surf
    real(dp), intent(out)     :: dlnT_dL
    real(dp), intent(out)     :: dlnT_dlnR
    real(dp), intent(out)     :: dlnT_dlnM
    real(dp), intent(out)     :: dlnT_dlnkap
    real(dp), intent(out)     :: lnP_surf
    real(dp), intent(out)     :: dlnP_dL
    real(dp), intent(out)     :: dlnP_dlnR
    real(dp), intent(out)     :: dlnP_dlnM
    real(dp), intent(out)     :: dlnP_dlnkap
    integer, intent(out)      :: ierr

    integer  :: T_tau_id
    real(dp) :: kap_guess

    include 'formats'

    ! Sanity check on L

    if (L < 0._dp) then
       s% retry_message = 'get_T_tau -- L <= 0'
       if (s% report_ierr) then
          write(*,2) 'get_T_tau: L <= 0', s% model_number, L
          !call mesa_error(__FILE__,__LINE__)
       end if
       ierr = -1
       return
    end if

    ! Get the T-tau id

    call get_T_tau_id(T_tau_relation, T_tau_id, ierr)
    if (ierr /= 0) then
       s% retry_message = 'Call to get_T_tau_id failed in get_T_tau'
       if (s% report_ierr) write(*, *) s% retry_message
       return
    end if

    ! Evaluate temperature and pressure by dispatching to the
    ! appropriate atm_lib routine, with the supplied T-tau relation

    select case (T_tau_opacity)

    case ('fixed')
       
       ! ok to use s% opacity(1) for fixed
       call atm_eval_T_tau_uniform( &
            tau_surf, L, R, M, cgrav, s% opacity(1), s% Pextra_factor, &
            T_tau_id, eos_proc_for_get_T_tau, kap_proc_for_get_T_tau, &
            s%atm_T_tau_errtol, 0, skip_partials, &
            Teff, kap, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
            ierr)

    case ('iterated')

       call prepare_kap(s, ierr)
       if (ierr /= 0) then
         s% retry_message = 'Call to prepare_kap failed in get_T_tau'
         if (s% report_ierr) write(*, *) s% retry_message
          return
       endif
       
       ! need to start iterations from same kap each time, so use opacity_start
       if (s% solver_iter > 0) then
          kap_guess = s% opacity_start(1)
       else
          kap_guess = s% opacity(1)
       end if
       call atm_eval_T_tau_uniform( &
            tau_surf, L, R, M, cgrav, kap_guess, s% Pextra_factor, &
            T_tau_id, eos_proc_for_get_T_tau, kap_proc_for_get_T_tau, &
            s%atm_T_tau_errtol, s%atm_T_tau_max_iters, skip_partials, &
            Teff, kap, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
            ierr)

    case ('varying')

       call prepare_kap(s, ierr)
       if (ierr /= 0) then
         s% retry_message = 'Call to prepare_kap failed in get_T_tau'
         if (s% report_ierr) write(*, *) s% retry_message
          return
       endif

       call atm_eval_T_tau_varying( &
            tau_surf, L, R, M, cgrav, &
            T_tau_id, eos_proc_for_get_T_tau, kap_proc_for_get_T_tau, &
            s% atm_T_tau_errtol, s% atm_T_tau_max_steps, skip_partials, &
            Teff, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
            ierr)

       kap = 0._dp ! This value is not used

    case default
       
       write(*,*) 'Unknown value for atm_T_tau_opacity: ' // trim(T_tau_opacity)
       call mesa_error(__FILE__,__LINE__,'Please amend your inlist file to correct this problem')

    end select

    ! Finish

    return

  contains

    subroutine eos_proc_for_get_T_tau( &
         lnP, lnT, &
         lnRho, res, dres_dlnRho, dres_dlnT, &
         ierr)

      real(dp), intent(in)  :: lnP
      real(dp), intent(in)  :: lnT
      real(dp), intent(out) :: lnRho
      real(dp), intent(out) :: res(:)
      real(dp), intent(out) :: dres_dlnRho(:)
      real(dp), intent(out) :: dres_dlnT(:)
      integer, intent(out)  :: ierr
      include 'formats'
      
      call eos_proc( &
           s, lnP, lnT, &
           lnRho, res, dres_dlnRho, dres_dlnT, &
           ierr)

    end subroutine eos_proc_for_get_T_tau

    !****

    subroutine kap_proc_for_get_T_tau( &
         lnRho, lnT, res, dres_dlnRho, dres_dlnT, &
         kap, dlnkap_dlnRho, dlnkap_dlnT, &
         ierr)

      real(dp), intent(in)  :: lnRho
      real(dp), intent(in)  :: lnT
      real(dp), intent(in)  :: res(:)
      real(dp), intent(in)  :: dres_dlnRho(:)
      real(dp), intent(in)  :: dres_dlnT(:)
      real(dp), intent(out) :: kap
      real(dp), intent(out) :: dlnkap_dlnRho
      real(dp), intent(out) :: dlnkap_dlnT
      integer, intent(out)  :: ierr
      include 'formats'
      
      call kap_proc( &
           s, lnRho, lnT, res, dres_dlnRho, dres_dlnT, &
           kap, dlnkap_dlnRho, dlnkap_dlnT, &
           ierr)

    end subroutine kap_proc_for_get_T_tau
    
  end subroutine get_T_tau

  !****

  subroutine get_T_tau_id (T_tau_relation, T_tau_id, ierr)

    use atm_def, only: &
         ATM_T_TAU_EDDINGTON, &
         ATM_T_TAU_SOLAR_HOPF, &
         ATM_T_TAU_KRISHNA_SWAMY, &
         ATM_T_TAU_TRAMPEDACH_SOLAR
    
    character(*), intent(in) :: T_tau_relation
    integer, intent(out)     :: T_tau_id
    integer, intent(out)     :: ierr

    ierr = 0

    ! Get the T-tau id

    select case (T_tau_relation)
    case ('Eddington')
       T_tau_id = ATM_T_TAU_EDDINGTON
    case ('solar_Hopf')
       T_tau_id = ATM_T_TAU_SOLAR_HOPF
    case ('Krishna_Swamy')
       T_tau_id = ATM_T_TAU_KRISHNA_SWAMY
    case ('Trampedach_solar')
       T_tau_id = ATM_T_TAU_TRAMPEDACH_SOLAR
    case default
       write(*,*) 'Unknown value for atm_T_tau_relation: ' // trim(T_tau_relation)
       call mesa_error(__FILE__,__LINE__,'Please amend your inlist file to correct this problem')
    end select

    ! Finish

    return

  end subroutine get_T_tau_id

  !****

  subroutine get_table( &
       s, skip_partials, L, R, M, cgrav, &
       Teff, &
       lnT, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
       lnP, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
       ierr)

    use atm_lib, only: &
         atm_eval_table, &
         atm_get_table_alfa_beta, &
         atm_get_T_tau_base, &
         atm_get_table_base

    type (star_info), pointer :: s
    logical, intent(in)       :: skip_partials
    real(dp), intent(in)      :: L, R, M, cgrav
    real(dp), intent(in)      :: Teff
    real(dp), intent(out)     :: lnT
    real(dp), intent(out)     :: dlnT_dL
    real(dp), intent(out)     :: dlnT_dlnR
    real(dp), intent(out)     :: dlnT_dlnM
    real(dp), intent(out)     :: dlnT_dlnkap
    real(dp), intent(out)     :: lnP
    real(dp), intent(out)     :: dlnP_dL
    real(dp), intent(out)     :: dlnP_dlnR
    real(dp), intent(out)     :: dlnP_dlnM
    real(dp), intent(out)     :: dlnP_dlnkap
    integer, intent(out)      :: ierr

    real(dp) :: Z
    real(dp) :: tau_base
    integer  :: T_tau_id
    integer  :: table_id
    real(dp) :: alfa
    real(dp) :: beta
    real(dp) :: lnT_a
    real(dp) :: dlnT_dL_a
    real(dp) :: dlnT_dlnR_a
    real(dp) :: dlnT_dlnM_a
    real(dp) :: dlnT_dlnkap_a
    real(dp) :: lnP_a
    real(dp) :: dlnP_dL_a
    real(dp) :: dlnP_dlnR_a
    real(dp) :: dlnP_dlnM_a
    real(dp) :: dlnP_dlnkap_a
    real(dp) :: lnT_b
    real(dp) :: dlnT_dL_b
    real(dp) :: dlnT_dlnR_b
    real(dp) :: dlnT_dlnM_b
    real(dp) :: dlnT_dlnkap_b
    real(dp) :: lnP_b
    real(dp) :: dlnP_dL_b
    real(dp) :: dlnP_dlnR_b
    real(dp) :: dlnP_dlnM_b
    real(dp) :: dlnP_dlnkap_b
    real(dp) :: kap_a
    real(dp) :: tau_b

    include 'formats'

    ! Check that tau_factor is correct

    if (s% tau_factor /= 1._dp) then
       write(*,*) 'Cannot use atm_option == ''table'' with tau_factor /= 1'
       call mesa_error(__FILE__,__LINE__,'Please amend your inlist file to correct this problem')
    end if

    Z = s% Z(1)

    ! Sanity check on L

    if (L < 0._dp) then
       s% retry_message = 'atm get_table: L < 0'
       if (s% report_ierr) then
          write(*,2) 'atm get_table: L < 0', s% model_number, L
          !call mesa_error(__FILE__,__LINE__)
       end if
       ierr = -1
       return
    end if

    ! Get the table id

    call get_table_id(s% atm_table, table_id, ierr)
    if (ierr /= 0) then
       s% retry_message = 'get_table_id failed in get_table'
       if (s% report_ierr) write(*, *) s% retry_message
       return
    end if

    ! Set up alfa and beta for doing table blending with off-table
    ! option

    call atm_get_table_alfa_beta( &
         L, Teff, R, M, cgrav, table_id, alfa, beta, ierr)
    if (ierr /= 0) then
       s% retry_message = 'Call to atm_get_table_alfa_beta failed in get_table'
       if (s% report_ierr) write(*, *) s% retry_message
       return
    endif
    
    ! If completely off the table, may need to reset tau_base to the
    ! T_Tau value to get the expected off-table behavior.
    if(beta == 0._dp) then
       call get_T_tau_id(s% atm_T_tau_relation, T_tau_id, ierr)
       if (ierr /= 0) then
          if (s% report_ierr) then
             write(*,*) 'Call to get_T_tau_id failed in get_table'
          endif
          return
       endif

       call atm_get_T_tau_base(T_tau_id, tau_base, ierr)
       if (ierr /= 0) then
          if (s% report_ierr) then
             write(*,*) 'Call to atm_get_T_tau_base failed in get_table'
          end if
          return
       end if

       if(s% tau_base /= tau_base) then
          write(*,*) "outside range for atm_table,"
          write(*,*) "setting T_tau value for tau_base =", tau_base
          s% tau_base = tau_base
       end if
    else ! check to see if need to switch back to table value for tau_base
       call atm_get_table_base(table_id, tau_base, ierr)
       if (ierr /= 0) then
          if (s% report_ierr) then
             write(*,*) 'Call to atm_get_table_base failed in get_table'
          end if
          return
       end if

       if(s% tau_base /= tau_base) then
          write(*,*) "back on atm_table,"
          write(*,*) "resetting to atm_table tau_base =", tau_base
          s% tau_base = tau_base
       end if
    end if
    
    ! Evaluate temperature and pressure from the table

    if (beta /= 0._dp) then

       call atm_eval_table( &
            L, R, M, cgrav, table_id, Z, skip_partials, &
            Teff, &
            lnT_b, dlnT_dL_b, dlnT_dlnR_b, dlnT_dlnM_b, dlnT_dlnkap_b, &
            lnP_b, dlnP_dL_b, dlnP_dlnR_b, dlnP_dlnM_b, dlnP_dlnkap_b, &
            ierr)
       if (ierr /= 0) then
          s% retry_message = 'Call to atm_eval_table failed in get_table'
          if (s% report_ierr) write(*, *) s% retry_message
          return
       end if

    else

       lnT_b = 0._dp
       dlnT_dL_b = 0._dp
       dlnT_dlnR_b = 0._dp
       dlnT_dlnM_b = 0._dp
       dlnT_dlnkap_b = 0._dp

       lnP_b = 0._dp
       dlnP_dL_b = 0._dp
       dlnP_dlnR_b = 0._dp
       dlnP_dlnM_b = 0._dp
       dlnP_dlnkap_b = 0._dp
       
    endif

    ! Evaluate temperature and pressure from the backup atmosphere
    ! option

    if (alfa /= 0._dp) then

       select case (s% atm_off_table_option)

       case ('T_tau')

          call get_T_tau( &
               s, s% tau_base, L, R, M, cgrav, &
               s% atm_T_tau_relation, s% atm_T_tau_opacity, skip_partials, &
               Teff, kap_a, &
               lnT_a, dlnT_dL_a, dlnT_dlnR_a, dlnT_dlnM_a, dlnT_dlnkap_a, &
               lnP_a, dlnP_dL_a, dlnP_dlnR_a, dlnP_dlnM_a, dlnP_dlnkap_a, &
               ierr)
          if (ierr /= 0) then
             s% retry_message = 'Call to get_T_tau failed in get_table'
             if (s% report_ierr) write(*, *) s% retry_message
             return
          end if

       case ('')

          write(*,*) 'Attempt to interpolate outside atmosphere table'
          call mesa_error(__FILE__,__LINE__,'Try setting the which_off_table_option control in your inlist file')
          stop

       case default

          write(*,*) 'Unknown value for atm_off_table_option: ' // trim(s% atm_off_table_option)
          call mesa_error(__FILE__,__LINE__,'Please amend your inlist file to correct this problem')
          
       end select

    else

       lnT_a = 0._dp
       dlnT_dL_a = 0._dp
       dlnT_dlnR_a = 0._dp
       dlnT_dlnM_a = 0._dp
       dlnT_dlnkap_a = 0._dp

       lnP_a = 0._dp
       dlnP_dL_a = 0._dp
       dlnP_dlnR_a = 0._dp
       dlnP_dlnM_a = 0._dp
       dlnP_dlnkap_a = 0._dp
       
    end if

    ! Blend the results together

    lnT = alfa*lnT_a + beta*lnT_b
    lnP = alfa*lnP_a + beta*lnP_b

    if (.not. skip_partials) then
       
       dlnT_dL = alfa*dlnT_dL_a + beta*dlnT_dL_b
       dlnT_dlnR = alfa*dlnT_dlnR_a + beta*dlnT_dlnR_b
       dlnT_dlnM = alfa*dlnT_dlnM_a + beta*dlnT_dlnM_b
       dlnT_dlnkap = alfa*dlnT_dlnkap_a + beta*dlnT_dlnkap_b
    
       dlnP_dL = alfa*dlnP_dL_a + beta*dlnP_dL_b
       dlnP_dlnR = alfa*dlnP_dlnR_a + beta*dlnP_dlnR_b
       dlnP_dlnM = alfa*dlnP_dlnM_a + beta*dlnP_dlnM_b
       dlnP_dlnkap = alfa*dlnP_dlnkap_a + beta*dlnP_dlnkap_b

    endif

    ! Set the effective temperature

    ! if (alfa /= 0._dp) then
    !    if (beta /= 0._dp) then
    !       if (Teff_a /= Teff_b) then
    !          ierr = -1
    !          s% retry_message = 'Mismatch between Teff values in get_tables'
    !          write(*,*) 'Mismatch between Teff values in get_tables: ', Teff_a, Teff_b
    !          !call mesa_error(__FILE__,__LINE__)
    !       end if
    !    endif
    !    Teff = Teff_a
    ! else
    !    Teff = Teff_b
    ! endif

    ! Finish

    return

  end subroutine get_table

  !****

  subroutine get_table_id (table_name, table_id, ierr)

    use atm_def, only: &
         ATM_TABLE_TAU_100, &
         ATM_TABLE_TAU_10, &
         ATM_TABLE_TAU_1, &
         ATM_TABLE_TAU_1M1, &
         ATM_TABLE_PHOTOSPHERE, &
         ATM_TABLE_WD_TAU_25, &
         ATM_TABLE_DB_WD_TAU_25

    character(*), intent(in) :: table_name
    integer, intent(out)     :: table_id
    integer, intent(out)     :: ierr

    ierr = 0

    ! Get the table id

    select case (table_name)
    case ('tau_100')
       table_id = ATM_TABLE_TAU_100
    case ('tau_10')
       table_id = ATM_TABLE_TAU_10
    case ('tau_1')
       table_id = ATM_TABLE_TAU_1
    case ('tau_1m1')
       table_id = ATM_TABLE_TAU_1M1
    case ('photosphere')
       table_id = ATM_TABLE_PHOTOSPHERE
    case ('WD_tau_25')
       table_id = ATM_TABLE_WD_TAU_25
    case ('DB_WD_tau_25')
       table_id = ATM_TABLE_DB_WD_TAU_25
    case default
       write(*,*) 'Unknown value for atm_table: ' // trim(table_name)
       call mesa_error(__FILE__,__LINE__,'Please amend your inlist file to correct this problem')
    end select

    ! Finish

    return

  end subroutine get_table_id

  !****
    
  subroutine get_irradiated( &
       s, irradiated_opacity, skip_partials, L, R, M, cgrav, &
       Teff, &
       lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
       lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
       ierr)

    use atm_lib, only: atm_eval_irradiated
    use kap_support, only: prepare_kap

    type (star_info), pointer :: s
    character(*), intent(in)  :: irradiated_opacity
    logical, intent(in)       :: skip_partials
    real(dp), intent(in)      :: L, R, M, cgrav
    real(dp), intent(in)      :: Teff
    real(dp), intent(out)     :: lnT_surf
    real(dp), intent(out)     :: dlnT_dL
    real(dp), intent(out)     :: dlnT_dlnR
    real(dp), intent(out)     :: dlnT_dlnM
    real(dp), intent(out)     :: dlnT_dlnkap
    real(dp), intent(out)     :: lnP_surf
    real(dp), intent(out)     :: dlnP_dL
    real(dp), intent(out)     :: dlnP_dlnR
    real(dp), intent(out)     :: dlnP_dlnM
    real(dp), intent(out)     :: dlnP_dlnkap
    integer, intent(out)      :: ierr

    real(dp) :: kap_for_atm
    real(dp) :: kap
    real(dp) :: tau_surf

    include 'formats'
    
    if (s% solver_iter > 0) then
       kap_for_atm = s% opacity_start(1)
    else
       kap_for_atm = s% opacity(1)
    end if

    ! Sanity check on L

    if (L < 0._dp) then
       s% retry_message = 'get_irradiated: L <= 0'
       if (s% report_ierr) then
          write(*,2) 'get_irradiated: L <= 0', s% model_number, L
          !call mesa_error(__FILE__,__LINE__)
       end if
       ierr = -1
       return
    end if

    ! Evaluate temperature and pressure by dispatching to the
    ! appropriate atm_lib routine

    select case (irradiated_opacity)

    case ('fixed')

       call atm_eval_irradiated( &
            L, R, M, cgrav, s% atm_irradiated_T_eq, s% atm_irradiated_P_surf, &
            kap_for_atm, s% atm_irradiated_kap_v, s% atm_irradiated_kap_v_div_kap_th, &
            eos_proc_for_get_irradiated, kap_proc_for_get_irradiated, &
            s% atm_irradiated_errtol, 0, skip_partials, &
            Teff, kap, tau_surf, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            ierr)

    case ('iterated')

       call prepare_kap(s, ierr)
       if (ierr /= 0) then
          s% retry_message = 'Failed in call to prepare_kap'
          if (s% report_ierr) write(*, *) s% retry_message
          return
       endif

       call atm_eval_irradiated( &
            L, R, M, cgrav, s% atm_irradiated_T_eq, s% atm_irradiated_P_surf, &
            kap_for_atm, s% atm_irradiated_kap_v, s% atm_irradiated_kap_v_div_kap_th, &
            eos_proc_for_get_irradiated, kap_proc_for_get_irradiated, &
            s% atm_irradiated_errtol, s% atm_irradiated_max_iters, skip_partials, &
            Teff, kap, tau_surf, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            ierr)
       
    case default

       write(*,*) 'Unknown value for atm_irradiated_opacity: ' // trim(irradiated_opacity)
       call mesa_error(__FILE__,__LINE__,'Please amend your inlist file to correct this problem')

    end select

    ! Set up remaining values

    lnP_surf = log(s% atm_irradiated_P_surf)

    dlnP_dL = 0._dp
    dlnP_dlnR = 0._dp
    dlnP_dlnM = 0._dp
    dlnP_dlnkap = 0._dp

    ! Update tau_factor

    s% tau_factor = tau_surf/s% tau_base
    
    ! Finish

    return

  contains

    subroutine eos_proc_for_get_irradiated( &
         lnP, lnT, &
         lnRho, res, dres_dlnRho, dres_dlnT, &
         ierr)

      real(dp), intent(in)  :: lnP
      real(dp), intent(in)  :: lnT
      real(dp), intent(out) :: lnRho
      real(dp), intent(out) :: res(:)
      real(dp), intent(out) :: dres_dlnRho(:)
      real(dp), intent(out) :: dres_dlnT(:)
      integer, intent(out)  :: ierr
      
      call eos_proc( &
           s, lnP, lnT, &
           lnRho, res, dres_dlnRho, dres_dlnT, &
           ierr)

    end subroutine eos_proc_for_get_irradiated

    !****

    subroutine kap_proc_for_get_irradiated( &
         lnRho, lnT, res, dres_dlnRho, dres_dlnT, &
         kap, dlnkap_dlnRho, dlnkap_dlnT, &
         ierr)

      real(dp), intent(in)  :: lnRho
      real(dp), intent(in)  :: lnT
      real(dp), intent(in)  :: res(:)
      real(dp), intent(in)  :: dres_dlnRho(:)
      real(dp), intent(in)  :: dres_dlnT(:)
      real(dp), intent(out) :: kap
      real(dp), intent(out) :: dlnkap_dlnRho
      real(dp), intent(out) :: dlnkap_dlnT
      integer, intent(out)  :: ierr
      
      call kap_proc( &
           s, lnRho, lnT, res, dres_dlnRho, dres_dlnT, &
           kap, dlnkap_dlnRho, dlnkap_dlnT, &
           ierr)

    end subroutine kap_proc_for_get_irradiated
    
  end subroutine get_irradiated

  !****

  subroutine get_legacy (s, ierr)

    type(star_info), pointer :: s
    integer, intent(out)     :: ierr

    ! Handles legacy / invalid choices for atm_option. This is
    ! to guide users toward the newer atmosphere controls; in the long
    ! term, it can be removed

    select case (s% atm_option)

    case ('simple_photosphere')

       write(*,*) 'The ''simple_photosphere'' choice for atm_option is deprecated'
       write(*,*) 'To achieve the same results, set:'
       write(*,*) '   atm_option = ''T_tau'''
       write(*,*) '   atm_T_tau_relation = ''Eddington'''
       write(*,*) '   atm_T_tau_opacity = ''fixed'''

    case ('grey_and_kap')

       write(*,*) 'The ''grey_and_kap'' choice for atm_option is deprecated'
       write(*,*) 'To achieve the same results, set:'
       write(*,*) '   atm_option = ''T_tau'''
       write(*,*) '   atm_T_tau_relation = ''Eddington'''
       write(*,*) '   atm_T_tau_opacity = ''iterated'''

    case ('Eddington_grey')

       write(*,*) 'The ''Eddington_grey'' choice for atm_option is deprecated'
       write(*,*) 'To achieve the same results, set:'
       write(*,*) '   atm_option = ''T_tau'''
       write(*,*) '   atm_T_tau_relation = ''Eddington'''
       write(*,*) '   atm_T_tau_opacity = ''varying'''

    case ('solar_Hopf_grey')

       write(*,*) 'The ''solar_Hopf_grey'' choice for atm_option is deprecated'
       write(*,*) 'To achieve the same results, set:'
       write(*,*) '   atm_option = ''T_tau'''
       write(*,*) '   atm_T_tau_relation = ''solar_Hopf'''
       write(*,*) '   atm_T_tau_opacity = ''varying'''

    case ('Krishna_Swamy')

       write(*,*) 'The ''Krishna_Swamy'' choice for atm_option is deprecated'
       write(*,*) 'To achieve the same results, set:'
       write(*,*) '   atm_option = ''T_tau'''
       write(*,*) '   atm_T_tau_relation = ''Krishna_Swamy'''
       write(*,*) '   atm_T_tau_opacity = ''varying'''

    case ('tau_100_tables')

       write(*,*) 'The ''tau_100_tables'' choice for atm_option is deprecated'
       write(*,*) 'To achieve the same results, set:'
       write(*,*) '   atm_option = ''table'''
       write(*,*) '   atm_table = ''tau_100'''

    case ('tau_10_tables')

       write(*,*) 'The ''tau_10_tables'' choice for atm_option is deprecated'
       write(*,*) 'To achieve the same results, set:'
       write(*,*) '   atm_option = ''table'''
       write(*,*) '   atm_table = ''tau_10'''

    case ('tau_1_tables')

       write(*,*) 'The ''tau_1_tables'' choice for atm_option is deprecated'
       write(*,*) 'To achieve the same results, set:'
       write(*,*) '   atm_option = ''table'''
       write(*,*) '   atm_table = ''tau_1'''

    case ('tau_1m1_tables')

       write(*,*) 'The ''tau_1m1_tables'' choice for atm_option is deprecated'
       write(*,*) 'To achieve the same results, set:'
       write(*,*) '   atm_option = ''table'''
       write(*,*) '   atm_table = ''tau_1m1'''

    case ('photosphere_tables')

       write(*,*) 'The ''tau_10_tables'' choice for atm_option is deprecated'
       write(*,*) 'To achieve the same results, set:'
       write(*,*) '   atm_option = ''table'''
       write(*,*) '   atm_table = ''photosphere'''

    case ('grey_irradiated')

       write(*,*) 'The ''grey_irradiated'' choice for atm_option is deprecated'
       write(*,*) 'To achieve the same results, set:'
       write(*,*) '   atm_option = ''irradiated_grey'''
       write(*,*) '   atm_irradiated_opacity = ''fixed'' (if atm_grey_irradiated_simple_kap_th = .true.)'
       write(*,*) '   atm_irradiated_opacity = ''varying'' (if atm_grey_irradiated_simple_kap_th = .false.)'

    case default
 
      write(*,*) 'Unknown value for atm_option: ' // trim(s% atm_option)

    end select

    ! Stop because we can't continue

    ierr = -1 ! ifort complains if this isn't set
    call mesa_error(__FILE__,__LINE__,'Please amend your inlist file to correct this problem')

    ! Finish

    return

  end subroutine get_legacy

  !****

  subroutine build_atm( &
       s, L, R, M, cgrav, ierr)

    type(star_info), pointer :: s
    real(dp), intent(in)     :: L, R, M, cgrav
    integer, intent(out)     :: ierr
       
    ! Create the atmosphere structure by dispatching to the
    ! appropriate internal routine

    select case (s% atm_option)

    case ('T_tau')

       call build_T_tau( &
            s, s% tau_factor*s% tau_base, L, R, M, cgrav, &
            s% atm_T_tau_relation, s% atm_T_tau_opacity, &
            s% atm_structure_num_pts, s% atm_structure, &
            ierr)

    case default

       write(*,*) 'Cannot create atm structure for atm_option: ' // TRIM(s% atm_option)
       call mesa_error(__FILE__,__LINE__,'Please amend your inlist file to correct this problem')

    end select

    ! Finish

  end subroutine build_atm

  !****

  subroutine build_T_tau( &
       s, tau_surf, L, R, M, cgrav, T_tau_relation, T_tau_opacity, &
       atm_structure_num_pts, atm_structure, &
       ierr)

    use atm_lib, only: &
         atm_build_T_tau_uniform, &
         atm_build_T_tau_varying
    
    type(star_info), pointer :: s
    real(dp), intent(in)     :: tau_surf, L, R, M, cgrav
    character(*), intent(in) :: T_tau_relation
    character(*), intent(in) :: T_tau_opacity
    integer, intent(out)     :: atm_structure_num_pts
    real(dp), pointer        :: atm_structure(:,:)
    integer, intent(out)     :: ierr

    real(dp) :: Teff
    real(dp) :: kap
    real(dp) :: lnT_surf
    real(dp) :: dlnT_dL
    real(dp) :: dlnT_dlnR
    real(dp) :: dlnT_dlnM
    real(dp) :: dlnT_dlnkap
    real(dp) :: lnP_surf
    real(dp) :: dlnP_dL
    real(dp) :: dlnP_dlnR
    real(dp) :: dlnP_dlnM
    real(dp) :: dlnP_dlnkap
    integer  :: T_tau_id

    ierr = 0

    ! First evaluate the temperature, pressure and opacity

    call get_T_tau( &
         s, tau_surf, L, R, M, cgrav, T_tau_relation, T_tau_opacity, .TRUE., &
         Teff, kap, &
         lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
         lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
         ierr)
    if (ierr /= 0) then
       s% retry_message = 'Call to get_T_tau failed in build_T_tau'
       if (s% report_ierr) write(*, *) s% retry_message
    end if

    ! Get the T-tau id

    call get_T_tau_id(T_tau_relation, T_tau_id, ierr)
    if (ierr /= 0) then
       s% retry_message = 'Call to get_T_tau_id failed in build_T_tau'
       if (s% report_ierr) write(*, *) s% retry_message
       return
    end if

    ! Build the atmosphere structure by dispatching to the
    ! appropriate atm_lib routine, with the supplied T-tau relation

    select case (T_tau_opacity)

    case ('fixed', 'iterated')

       call atm_build_T_tau_uniform( &
            tau_surf, L, R, M, cgrav, kap, s% Pextra_factor, s% atm_build_tau_outer, &
            T_tau_id, eos_proc_for_build_T_tau, kap_proc_for_build_T_tau, &
            s% atm_build_errtol, s% atm_build_dlogtau, &
            atm_structure_num_pts, atm_structure, &
            ierr)
       if (ierr /= 0) then
          s% retry_message = 'Call to atm_build_T_tau_uniform failed in build_T_tau'
          if (s% report_ierr) write(*, *) s% retry_message
          return
       end if

    case ('varying')

       call atm_build_T_tau_varying( &
            tau_surf, L, R, M, cgrav, lnP_surf, s% atm_build_tau_outer, &
            T_tau_id, eos_proc_for_build_T_tau, kap_proc_for_build_T_tau, &
            s% atm_build_errtol, s% atm_build_dlogtau, &
            atm_structure_num_pts, atm_structure, &
            ierr)
       if (ierr /= 0) then
          s% retry_message = 'Call to atm_build_T_tau_varying failed in build_T_tau'
          if (s% report_ierr) write(*, *) s% retry_message
          return
       end if

    case default
       
       write(*,*) 'Unknown value for atm_T_tau_opacity: ' // trim(T_tau_opacity)
       call mesa_error(__FILE__,__LINE__,'Please amend your inlist file to correct this problem')

    end select

    ! Finish

    return

  contains

    subroutine eos_proc_for_build_T_tau( &
         lnP, lnT, &
         lnRho, res, dres_dlnRho, dres_dlnT, &
         ierr)

      real(dp), intent(in)  :: lnP
      real(dp), intent(in)  :: lnT
      real(dp), intent(out) :: lnRho
      real(dp), intent(out) :: res(:)
      real(dp), intent(out) :: dres_dlnRho(:)
      real(dp), intent(out) :: dres_dlnT(:)
      integer, intent(out)  :: ierr
      
      call eos_proc( &
           s, lnP, lnT, &
           lnRho, res, dres_dlnRho, dres_dlnT, &
           ierr)

    end subroutine eos_proc_for_build_T_tau

    !****

    subroutine kap_proc_for_build_T_tau( &
         lnRho, lnT, res, dres_dlnRho, dres_dlnT, &
         kap, dlnkap_dlnRho, dlnkap_dlnT, &
         ierr)

      real(dp), intent(in)  :: lnRho
      real(dp), intent(in)  :: lnT
      real(dp), intent(in)  :: res(:)
      real(dp), intent(in)  :: dres_dlnRho(:)
      real(dp), intent(in)  :: dres_dlnT(:)
      real(dp), intent(out) :: kap
      real(dp), intent(out) :: dlnkap_dlnRho
      real(dp), intent(out) :: dlnkap_dlnT
      integer, intent(out)  :: ierr
      
      call kap_proc( &
           s, lnRho, lnT, res, dres_dlnRho, dres_dlnT, &
           kap, dlnkap_dlnRho, dlnkap_dlnT, &
           ierr)

    end subroutine kap_proc_for_build_T_tau
    
  end subroutine build_T_tau
    
  !****

  subroutine eos_proc( &
       s, lnP, lnT, &
       lnRho, res, dres_dlnRho, dres_dlnT, &
       ierr)

    use eos_def, only: num_eos_basic_results
    use eos_lib, only: radiation_pressure, eos_gamma_PT_get_rho_energy
    use eos_support, only: solve_eos_given_PgasT

    type(star_info), pointer :: s
    real(dp), intent(in)     :: lnP
    real(dp), intent(in)     :: lnT
    real(dp), intent(out)    :: lnRho
    real(dp), intent(out)    :: res(num_eos_basic_results)
    real(dp), intent(out)    :: dres_dlnRho(num_eos_basic_results)
    real(dp), intent(out)    :: dres_dlnT(num_eos_basic_results)
    integer, intent(out)     :: ierr

    real(dp), parameter :: LOGRHO_TOL = 1E-11_dp
    real(dp), parameter :: LOGPGAS_TOL = 1E-11_dp

    real(dp) :: T
    real(dp) :: P
    real(dp) :: Prad
    real(dp) :: Pgas
    real(dp) :: rho, gamma, energy
    real(dp) :: logRho, logRho_guess
    real(dp) :: dlnRho_dlnPgas
    real(dp) :: dlnRho_dlnT
    real(dp) :: dres_dxa(num_eos_basic_results,s% species)

    T = exp(lnT)
    P = exp(lnP)

    Prad = radiation_pressure(T)
    Pgas = MAX(1.E-99_dp, P - Prad)
       
    gamma = 5d0/3d0
    call eos_gamma_PT_get_rho_energy( &
       s% abar(1), P, T, gamma, rho, energy, ierr)
    if (ierr /= 0) then
       s% retry_message = 'Call to eos_gamma_PT_get_rho_energy failed in eos_proc'
       if (s% report_ierr) write(*, *) trim(s% retry_message)
    end if
    
    logRho_guess = log10(rho)

    call solve_eos_given_PgasT( &
         s, 0, s% xa(:,1), &
         lnT/ln10, log10(Pgas), logRho_guess, LOGRHO_TOL, LOGPGAS_TOL, &
         logRho, res, dres_dlnRho, dres_dlnT, dres_dxa, &
         ierr)
    if (ierr /= 0) then
       s% retry_message = 'Call to solve_eos_given_PgasT failed in eos_proc'
       if (s% report_ierr) write(*, *) trim(s% retry_message)
    end if

    lnRho = logRho*ln10

    ! Finish
    

    return

  end subroutine eos_proc

  !****

  subroutine kap_proc( &
       s, lnRho, lnT, res, dres_dlnRho, dres_dlnT, &
       kap, dlnkap_dlnRho, dlnkap_dlnT, &
       ierr)

    use eos_def, only: i_lnfree_e, i_eta
    use kap_def, only: num_kap_fracs
    use kap_support, only: get_kap

    type(star_info), pointer :: s
    real(dp), intent(in)     :: lnRho
    real(dp), intent(in)     :: lnT
    real(dp), intent(in)     :: res(:)
    real(dp), intent(in)     :: dres_dlnRho(:)
    real(dp), intent(in)     :: dres_dlnT(:)
    real(dp), intent(out)    :: kap
    real(dp), intent(out)    :: dlnkap_dlnRho
    real(dp), intent(out)    :: dlnkap_dlnT
    integer, intent(out)     :: ierr

    real(dp) :: kap_fracs(num_kap_fracs)
    
    include 'formats'

    call get_kap( &
         s, 0, s% zbar(1), s% xa(:,1), lnRho/ln10, lnT/ln10, &
         res(i_lnfree_e), dres_dlnRho(i_lnfree_e), dres_dlnT(i_lnfree_e), &
         res(i_eta), dres_dlnRho(i_eta), dres_dlnT(i_eta), &
         kap_fracs, kap, dlnkap_dlnRho, dlnkap_dlnT, &
         ierr)
    if (ierr /= 0) then
       s% retry_message = 'Call to get_kap failed in kap_proc'
       if (s% report_ierr) write(*, *) s% retry_message
    end if

    ! Finish
    
    if (kap <= 0d0 .or. is_bad(kap)) then
       write(*,1) 'bad kap', kap
       write(*,1) 's% zbar(1)', s% zbar(1)
       write(*,1) 'lnRho/ln10', lnRho/ln10
       write(*,1) 'lnT/ln10', lnT/ln10
       write(*,1) 'res(i_eta)', res(i_eta)
       ierr = -1
       return
       call mesa_error(__FILE__,__LINE__,'atm kap_proc')
    end if

    return

  end subroutine kap_proc

end module atm_support
