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

module micro

  ! Uses

  use star_private_def
  use const_def
  use star_utils, only: foreach_cell
  use utils_lib, only: is_bad

  ! No implicit typing

  implicit none

  ! Parameters

  logical, parameter :: dbg = .false.
  real(dp), parameter :: fmax_fit = 3.47189d0

  ! Access specifiers

  private

  public :: set_micro_vars
  public :: set_eos_with_mask
  public :: do_eos_for_cell
  public :: store_eos_for_cell
  public :: do_kap_for_cell
  public :: shutdown_microphys

contains

  subroutine set_micro_vars( &
       s, nzlo, nzhi, skip_eos, skip_net, skip_neu, skip_kap, ierr)

    use kap_support, only: prepare_kap
    use star_utils, only: start_time, update_time
    use net_lib, only: net_work_size
    use net, only: do_net
    use chem_def, only: icno, ipp
    use rates_def, only: i_rate

    type (star_info), pointer :: s
    integer, intent(in) :: nzlo, nzhi
    logical, intent(in) :: skip_eos, skip_net, skip_neu, skip_kap
    integer, intent(out) :: ierr

    integer :: j, k, op_err, k_bad, res
    integer(8) :: time0
    real(dp) :: total, alfa, beta, X_lo, X_lo_fac, X_hi, X_hi_fac

    include 'formats'

    ierr = 0
    if (dbg) then
       write(*,*)
       write(*,*) 'set_micro_vars'
       write(*,*) 'nzlo', nzlo
       write(*,*) 'nzhi', nzhi
       write(*,*) 'skip_net', skip_net
       write(*,*) 'skip_kap', skip_kap
    end if

    if (.not. skip_eos) then
       call set_eos(ierr)
       if (ierr /= 0) return
    end if

    do k=nzlo, nzhi
       if (s% csound_start(k) < 0) then
          s% csound_start(k) = s% csound(k)
       end if
    end do

    do k=nzlo, nzhi
       if (k == 1) then
          s% rho_face(k) = s% rho(k)
          if (.not. s% u_flag) s% P_face(k) = s% P(k)
          s% csound_face(1) = s% csound(1)
       else
          alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
          beta = 1 - alfa
          s% rho_face(k) = alfa*s% rho(k) + beta*s% rho(k-1)
          if (.not. s% u_flag) s% P_face(k) = alfa*s% P(k) + beta*s% P(k-1)
          s% csound_face(k) = alfa*s% csound(k) + beta*s% csound(k-1)
       end if
    end do

    if (.not. (skip_kap .and. skip_neu)) then

       if (.not. skip_kap) then
          call prepare_kap(s, ierr)
          if (ierr /= 0) return
          if (s% use_other_opacity_factor) then
             call s% other_opacity_factor(s% id, ierr)
             if (ierr /= 0) return
          else
             s% extra_opacity_factor(1:s% nz) = s% opacity_factor
          end if
       endif
       if (s% use_other_opacity_factor) then
          call s% other_opacity_factor(s% id, ierr)
          if (ierr /= 0) return
       else
          s% extra_opacity_factor(1:s% nz) = s% opacity_factor
       end if

       if (s% doing_timing) call start_time(s, time0, total)
       !!$OMP PARALLEL DO PRIVATE(k,op_err) SCHEDULE(dynamic,2)
       !$OMP PARALLEL DO PRIVATE(k,op_err) SCHEDULE(dynamic,2)
       do k = nzlo, nzhi
          op_err = 0
          call do1_neu_kap(s,k,op_err)
          if (op_err /= 0) ierr = op_err
       end do
       !$OMP END PARALLEL DO
       if (s% doing_timing) call update_time(s, time0, total, s% time_neu_kap)
       if (ierr /= 0) then
          if (s% report_ierr) write(*,*) 'do1_neu_kap returned ierr', ierr
          return
       end if

    end if

    if (ierr /= 0) return

    if (.not. skip_net) then

       if (dbg) write(*,*) 'micro: call do_net'
       call do_net(s, nzlo, nzhi, .false., ierr)
       if (dbg) write(*,*) 'micro: done do_net'
       if (ierr /= 0) then
          if (s% report_ierr) write(*,*) 'do_net returned ierr', ierr
          return
       end if

    end if

  contains

    subroutine set_eos(ierr)
      integer, intent(out) :: ierr
      integer :: k
      real(dp) :: alfa
      include 'formats'
      ierr = 0
      if (dbg) write(*,*) 'call do_eos'
      if (s% doing_timing) call start_time(s, time0, total)
      call do_eos(s,nzlo,nzhi,ierr)
      if (s% doing_timing) call update_time(s, time0, total, s% time_eos)
      if (ierr /= 0) then
         if (s% report_ierr) write(*,*) 'do_eos returned ierr', ierr
         return
      end if
    end subroutine set_eos

    subroutine debug(str)
      use chem_def
      character (len=*), intent(in) :: str
      integer :: k, j
      include 'formats'
      k = 1469
      do j=1,1 !s% species
         if (.true. .or. s% xa(j,k) > 1d-9) &
              write(*,1) trim(str) // ' xin(net_iso(i' // &
              trim(chem_isos% name(s% chem_id(j))) // '))= ', s% xa(j,k)
      end do
    end subroutine debug

    subroutine do1_neu_kap(s,k,ierr)
      use neu_def,only:Tmin_neu
      use neu, only: do_neu_for_cell, do_clear_neu_for_cell
      type (star_info), pointer :: s
      integer, intent(in) :: k
      integer, intent(out) :: ierr
      if (s% T(k) >= Tmin_neu .and. s% non_nuc_neu_factor > 0d0) then
         call do_neu_for_cell(s,k,ierr)
      else
         call do_clear_neu_for_cell(s,k,ierr)
      end if
      if (ierr /= 0) return
      if (skip_kap) return
      call do_kap_for_cell(s,k,ierr)
      if (ierr /= 0) return
    end subroutine do1_neu_kap

  end subroutine set_micro_vars

  !****

  subroutine do_eos(s,nzlo,nzhi,ierr)

    type (star_info), pointer :: s
    integer, intent(in) :: nzlo, nzhi
    integer, intent(out) :: ierr
    logical, parameter :: use_omp = .true.
    include 'formats'
    integer :: k

    ierr = 0

    if (dbg) write(*,*) 'do_eos call foreach_cell', nzlo, nzhi

    call foreach_cell(s,nzlo,nzhi,use_omp,do_eos_for_cell,ierr)

  end subroutine do_eos

  !****

  subroutine set_eos_with_mask (s,nzlo,nzhi,mask,ierr)

    type (star_info), pointer :: s
    integer, intent(in) :: nzlo, nzhi
    logical, intent(in) :: mask(:)
    integer, intent(out) :: ierr
    include 'formats'
    integer :: k
    logical :: lerr(s%nz)

    ierr = 0

    if (dbg) write(*,*) 'set_eos_with_mask'

    !$OMP PARALLEL DO PRIVATE(ierr) SCHEDULE(dynamic,2)
    do k = nzlo, nzhi
       if (mask(k)) then
          call do_eos_for_cell(s, k, ierr)
          lerr(k) = ierr /= 0
       else
          lerr(k) = .FALSE.
       endif
    end do
    !$OMP END PARALLEL DO

    if (ANY(lerr(nzlo:nzhi))) then
       ierr = 1
    else
       ierr = 0
    endif

  end subroutine set_eos_with_mask

  !****
  
  subroutine do_eos_for_cell(s,k,ierr)

    use const_def
    use chem_def
    use chem_lib
    use eos_def
    use eos_support
    use net_def, only: net_general_info
    use star_utils, only: write_eos_call_info, lookup_nameofvar

    type (star_info), pointer :: s
    integer, intent(in) :: k
    integer, intent(out) :: ierr

    real(dp), dimension(num_eos_basic_results) :: &
         res, res_a, res_b, d_dlnd, d_dlnT, d_eos_dabar, d_eos_dzbar
    real(dp) :: &
         sumx, dx, dxh_a, dxh_b, &
         Rho, logRho, lnd, lnE, logT, T, energy, logQ, frac
    integer, pointer :: net_iso(:)
    integer :: j, species, i_var, i_var_sink
    real(dp), parameter :: epsder = 1d-4, Z_limit = 0.5d0

    logical, parameter :: testing = .false.

    include 'formats'

    ierr = 0

    net_iso => s% net_iso
    species = s% species
    call basic_composition_info( &
         species, s% chem_id, s% xa(1:species,k), s% X(k), s% Y(k), s% Z(k), &
         s% abar(k), s% zbar(k), s% z2bar(k), s% z53bar(k), &
         s% ye(k), s% mass_correction(k), sumx)

    logT = s% lnT(k)/ln10

    logRho = s% lnd(k)/ln10

    if (s% solver_test_eos_partials) then
       eos_test_partials = (k == s% solver_test_partials_k .and. &
          s% solver_call_number == s% solver_test_partials_call_number .and. &
          s% solver_iter == s% solver_test_partials_iter_number )
    end if

    call get_eos( &
         s, k, s% xa(:,k), &
         s% rho(k), logRho, s% T(k), logT, &
         res, s% d_eos_dlnd(:,k), s% d_eos_dlnT(:,k), &
         s% d_eos_dxa(:,:,k), ierr)
    if (ierr /= 0) then
       if (s% report_ierr) then
          write(*,*) 'do_eos_for_cell: get_eos ierr', ierr
          !stop 'do_eos_for_cell'
       end if
       return
    end if
    
    if (s% solver_test_eos_partials .and. eos_test_partials) then
       s% solver_test_partials_val = eos_test_partials_val
       s% solver_test_partials_dval_dx = eos_test_partials_dval_dx
    end if

    s% lnPgas(k) = res(i_lnPgas)
    s% Pgas(k) = exp(s% lnPgas(k))

    call store_stuff(ierr)
    if (ierr /= 0) return

    if (s% solver_test_partials_write_eos_call_info .and. &
       k == s% solver_test_partials_k .and. &
       s% solver_call_number == s% solver_test_partials_call_number .and. &
       s% solver_iter == s% solver_test_partials_iter_number) then
       call write_eos_call_info(s,k)
       stop 'do_eos_for_cell: write_eos_call_info'
    end if

  contains

    subroutine store_stuff(ierr)

      integer, intent(out) :: ierr

      include 'formats'

      call store_eos_for_cell(s, k, &
         res, s% d_eos_dlnd(:,k), s% d_eos_dlnT(:,k), &
         s% d_eos_dxa(:,:,k), ierr)

      if (ierr /= 0) then
         if (s% report_ierr) then
            call write_eos_call_info(s,k)
            write(*,2) 'store_eos_for_cell failed', k
            if (s% stop_for_bad_nums) stop 'do_eos_for_cell'
            return
         end if
         return
      end if

      if (k == s% trace_k) then
         write(*,5) 'grada', k, s% solver_iter, s% solver_adjust_iter, &
              s% model_number, s% grada(k)
      end if
      if (s% model_number == -1) then
         write(*,4) 'grada', k, s% solver_iter, s% model_number, s% grada(k)
      end if

    end subroutine store_stuff

  end subroutine do_eos_for_cell

  !****

  subroutine store_eos_for_cell(s, k, res, d_dlnd, d_dlnT, d_dxa, ierr)

    use eos_def
    use chem_def
    use star_utils, only: eval_csound, write_eos_call_info

    type (star_info), pointer :: s
    integer, intent(in) :: k
    real(dp), intent(in), dimension(num_eos_basic_results) :: res, d_dlnd, d_dlnT
    real(dp), intent(in) :: d_dxa(num_eos_d_dxa_results,s% species)
    integer, intent(out) :: ierr
    integer :: i, j

    real(dp) :: T, rho, &
         dlnd_dV, Pgas, Prad, &
         P, PV, PT, E, EV, ET, CP, dCp_dV, dCp_dT, &
         chiT, dchiT_dlnd, dchiT_dlnT, &
         chiRho, dchiRho_dlnd, dchiRho_dlnT, &
         Q, dQ_dlnd, dQ_dlnT, QV, QT

    include 'formats'

    ierr = 0

    ! check values
    do i = 1, num_eos_basic_results
       if (is_bad(res(i) + d_dlnd(i) + d_dlnT(i))) then
          ierr = -1
          if (s% report_ierr) then
             !$OMP critical (micro_crit0)
             write(*,2) trim(eosDT_result_names(i)), k, res(i)
             write(*,2) 'd_dlnd ' // trim(eosDT_result_names(i)), k, d_dlnd(i)
             write(*,2) 'd_dlnT ' // trim(eosDT_result_names(i)), k, d_dlnT(i)
             write(*,*)
             call write_eos_call_info(s,k)
             if (s% stop_for_bad_nums) stop 'store_eos_for_cell'
             !$OMP end critical (micro_crit0)
          end if
          return
       end if
    end do

    T = s% T(k)
    rho = s% rho(k)
    if (T > 1d15 .or. rho > 1d15) then
       ierr = -1
       if (s% report_ierr) then
          write(*,4) 'bad T or rho for eos', k, s% solver_iter, s% model_number
          write(*,2) 'T', k, T
          write(*,2) 'rho', k, rho
       end if
       return
    end if
    s% Prad(k) = crad * T*T*T*T / 3
    s% P(k) = s% Prad(k) + s% Pgas(k)
    s% lnP(k) = log(s% P(k))
    s% lnS(k) = res(i_lnS)
    s% lnE(k) = res(i_lnE)
    s% energy(k) = exp(s% lnE(k))
    s% entropy(k) = exp(s% lnS(k))
    s% grada(k) = res(i_grad_ad)
    s% dE_dRho(k) = res(i_dE_drho)
    s% Cv(k) = res(i_Cv)
    s% cp(k) = res(i_cp)
    s% chiRho(k) = res(i_chiRho)
    s% chiT(k) = res(i_chiT)
    s% gamma1(k) = res(i_gamma1)
    s% gamma3(k) = res(i_gamma3)
    s% eta(k) = res(i_eta)
    s% gam(k) = s% z53bar(k)*qe*qe * &
         pow((4d0/3d0)*pi*avo*rho*s% zbar(k)/s% abar(k),one_third) / (kerg*T)
    s% mu(k) = res(i_mu)
    s% lnfree_e(k) = res(i_lnfree_e)

    s% phase(k) = res(i_phase)
    s% latent_ddlnT(k) = res(i_latent_ddlnT)
    s% latent_ddlnRho(k) = res(i_latent_ddlnRho)

    s% eos_frac_OPAL_SCVH(k) = res(i_frac_OPAL_SCVH)
    s% eos_frac_HELM(k) = res(i_frac_HELM)
    s% eos_frac_Skye(k) = res(i_frac_Skye)
    s% eos_frac_PC(k) = res(i_frac_PC)
    s% eos_frac_FreeEOS(k) = res(i_frac_FreeEOS)
    s% eos_frac_CMS(k) = res(i_frac_CMS)

    s% chiRho_for_partials(k) = s% Pgas(k)*d_dlnd(i_lnPgas)/s% P(k)
    s% chiT_for_partials(k) = (s% Pgas(k)*d_dlnT(i_lnPgas) + 4d0*s% Prad(k))/s% P(k)
    s% dE_drho_for_partials(k) = d_dlnd(i_lnE)*s% energy(k)/s% rho(k)
    s% Cv_for_partials(k) = d_dlnT(i_lnE)*s% energy(k)/s% T(k)
    s% dS_drho_for_partials(k) = d_dlnd(i_lnS)*s% entropy(k)/s% rho(k)
    s% dS_dT_for_partials(k) = d_dlnT(i_lnS)*s% entropy(k)/s% T(k)
    do j=1, s% species
       s% dlnE_dxa_for_partials(j,k) = d_dxa(i_lnE,j)
       s% dlnP_dxa_for_partials(j,k) = s% Pgas(k)*d_dxa(i_lnPgas,j)/s% P(k)
    end do
    
    s% QQ(k) = s% chiT(k)/(s% rho(k)*s% T(k)*s% chiRho(k)) ! thermal expansion coefficient
    s% d_QQ_dlnd(k) = s% QQ(k)*(d_dlnd(i_chiT)/s% chiT(k) - d_dlnd(i_chiRho)/s% chiRho(k) - 1d0)
    s% d_QQ_dlnT(k) = s% QQ(k)*(d_dlnT(i_chiT)/s% chiT(k) - d_dlnT(i_chiRho)/s% chiRho(k) - 1d0)
    s% csound(k) = eval_csound(s,k,ierr)

    ! check some key values
    if (s% gamma1(k) <= 0 .or. &
       s% grada(k) <= 0 .or. &
       s% csound(k) <= 0 .or. is_bad(s% csound(k)) .or. &
       s% chiT(k) <= 0 .or. &
       s% chiRho(k) <= 0 .or. &
       s% Cv(k) <= 0 .or. &
       s% Cp(k) <= 0) then
       ierr = -1
    end if

    if (ierr /= 0) then
       if (s% report_ierr) then
          !$OMP critical (micro_crit1)
          write(*,2) 's% cp(k)', k, s% cp(k)
          write(*,2) 's% csound(k)', k, s% csound(k)
          write(*,2) 's% lnP(k)', k, s% lnP(k)
          write(*,2) 's% gam(k)', k, s% gam(k)
          write(*,2) 's% P(k)', k, s% P(k)
          write(*,2) 's% Pgas(k)', k, s% Pgas(k)
          write(*,2) 's% rho(k)', k, s% rho(k)
          write(*,2) 's% T(k)', k, s% T(k)
          write(*,2) 'logRho', k, s% lnd(k)/ln10
          write(*,2) 'logT', k, s% lnT(k)/ln10
          write(*,2) 'abar', k, s% abar(k)
          write(*,2) 'zbar', k, s% zbar(k)
          write(*,*) 
          call write_eos_call_info(s,k)
          stop 'store_eos_for_cell'
          !$OMP end critical (micro_crit1)
       end if
       ierr = -1
    end if

  end subroutine store_eos_for_cell

  !****

  subroutine do_kap_for_cell(s,k,ierr)

    use const_def,only:ln10
    use net_def,only:net_general_info
    use rates_def, only:i_rate
    use chem_def
    use kap_def
    use kap_lib
    use kap_support, only: fraction_of_op_mono, get_kap
    use eos_def, only : i_lnfree_e, i_eta
    use eos_lib, only: Radiation_Pressure
    use star_utils, only: get_XYZ

    type (star_info), pointer :: s
    integer, intent(in) :: k
    integer, intent(out) :: ierr

    integer, pointer :: net_iso(:)
    integer :: i, iz, kh
    real(dp) :: &
         log10_rho, log10_T, dlnkap_dlnd, dlnkap_dlnT, &
         opacity_max, opacity_max0, opacity_max1, zbar, &
         lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
         eta, d_eta_dlnRho, d_eta_dlnT, &
         log_r, log_r_in, log_r_out, log_r_frac, frac, min_cno_for_kap_limit, &
         P, Prad, Pgas, Ledd_factor, Ledd_kap, Ledd_log, &
         a, b, da_dlnd, da_dlnT, db_dlnd, db_dlnT, opacity_factor
    real(dp), dimension(num_kap_fracs) :: kap_fracs
    character (len=100) :: message
    real(dp), pointer :: xa(:)
    logical :: test_partials

    include 'formats'

    ierr = 0

    log10_rho = s% lnd(k)/ln10
    log10_T = s% lnT(k)/ln10

    lnfree_e = s% lnfree_e(k)
    d_lnfree_e_dlnRho = s% d_eos_dlnd(i_lnfree_e,k)
    d_lnfree_e_dlnT = s% d_eos_dlnT(i_lnfree_e,k)

    eta = s% eta(k)
    d_eta_dlnRho = s% d_eos_dlnd(i_eta,k)
    d_eta_dlnT = s% d_eos_dlnT(i_eta,k)

    if (s% use_starting_composition_for_kap) then
       xa(1:s% species) => s% xa_start(1:s% species,k)
       zbar = s% zbar_start(k)
    else
       xa(1:s% species) => s% xa(1:s% species,k)
       zbar = s% zbar(k)
    end if
    
    call get_kap( &
         s, k, zbar, xa, log10_rho, log10_T, &
         lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
         eta, d_eta_dlnRho, d_eta_dlnT, &
         kap_fracs, s% opacity(k), dlnkap_dlnd, dlnkap_dlnT, ierr)

    ! unpack fractions
    s% kap_frac_lowT(k) = kap_fracs(i_frac_lowT)
    s% kap_frac_highT(k) = kap_fracs(i_frac_highT)
    s% kap_frac_Type2(k) = kap_fracs(i_frac_Type2)
    s% kap_frac_Compton(k) = kap_fracs(i_frac_Compton)

    if (is_bad_num(s% opacity(k)) .or. ierr /= 0) then
       if (s% report_ierr) then
          !$omp critical (star_kap_get)
          call show_stuff()
          stop 'debug1: do_kap_for_cell'
          !$omp end critical (star_kap_get)
       end if
       ierr = -1
       return
    end if

    if (s% opacity(k) < 1d-99) then
       s% opacity(k) = 1d-99
       dlnkap_dlnd = 0
       dlnkap_dlnT = 0
    end if

    opacity_factor = s% extra_opacity_factor(k)
    if (s% min_logT_for_opacity_factor_off > 0) then
       if (log10_T >= s% max_logT_for_opacity_factor_off .or. &
            log10_T <= s% min_logT_for_opacity_factor_off) then
          opacity_factor = 1
       else if (log10_T > s% max_logT_for_opacity_factor_on) then
          opacity_factor = 1 + (opacity_factor-1)* &
               (log10_T - s% max_logT_for_opacity_factor_off)/ &
               (s% max_logT_for_opacity_factor_on - s% max_logT_for_opacity_factor_off)
       else if (log10_T < s% min_logT_for_opacity_factor_on) then
          opacity_factor = 1 + (opacity_factor-1)* &
               (log10_T - s% min_logT_for_opacity_factor_off)/ &
               (s% min_logT_for_opacity_factor_on - s% min_logT_for_opacity_factor_off)
       end if
    end if

    s% opacity(k) = s% opacity(k)*opacity_factor
    s% d_opacity_dlnd(k) = s% opacity(k)*dlnkap_dlnd
    s% d_opacity_dlnT(k) = s% opacity(k)*dlnkap_dlnT

    if (s% opacity(k) > s% opacity_max .and. s% opacity_max > 0) then
       s% opacity(k) = s% opacity_max
       s% d_opacity_dlnd(k) = 0
       s% d_opacity_dlnT(k) = 0
    end if

    if (ierr /= 0 .or. is_bad_num(s% opacity(k))) then
       if (s% stop_for_bad_nums) then
          !$omp critical (star_kap_get_bad_num)
          write(*,*)
          write(*,2) 's% opacity(k)', k, s% opacity(k)
          write(*,2) 's% kap_frac_Type2(k)', k, s% kap_frac_Type2(k)
          write(*,*)
          call show_stuff()
          stop 'do_kap_for_cell'
          !$omp end critical (star_kap_get_bad_num)
       end if
       if (s% report_ierr) then
          return
          !$omp critical (star_kap_get_bad_num2)
          write(*,*) 'do_kap_for_cell: kap_get failure for cell ', k
          call show_stuff()
          stop 'debug: do_kap_for_cell'
          !$omp end critical (star_kap_get_bad_num2)
       end if
       ierr = -1
       return
    end if


      !test_partials = (k == s% solver_test_partials_k)
      test_partials = .false.

      if (test_partials) then
         s% solver_test_partials_val = s% opacity(k)
         s% solver_test_partials_var = s% i_lnd
         s% solver_test_partials_dval_dx = s% d_opacity_dlnd(k)
         write(*,*) 'do_kap_for_cell', s% solver_test_partials_var
      end if



  contains

    subroutine show_stuff()
      use star_utils, only: write_eos_call_info

      include 'formats'

      real(dp) :: xc, xo, xh, xhe, Z
      integer :: i, iz
      
      write(*,*) 'eos info'
      call write_eos_call_info(s,k)
      
      write(*,*) 'kap info'
      call get_XYZ(s, s% xa(:,k), xh, xhe, Z)
      xc = 0; xo = 0
      do i=1, s% species
         iz = chem_isos% Z(s% chem_id(i))
         select case(iz)
         case (6)
            xc = xc + s% xa(i,k)
         case (8)
            xo = xo + s% xa(i,k)
         end select
      end do
      write(*,2) 'show opacity info', k
      write(*,*) 'kap_option ' // trim(kap_option_str(s% kap_rq% kap_option))
      write(*,*) 'kap_CO_option ' // trim(kap_CO_option_str(s% kap_rq% kap_CO_option))
      write(*,*) 'kap_lowT_option ' // trim(kap_lowT_option_str(s% kap_rq% kap_lowT_option))
      write(*,1) 'logT = ', log10_T
      write(*,1) 'logRho = ', log10_rho
      write(*,1) 'z = ', z
      write(*,1) 'Zbase = ', s% kap_rq% Zbase
      write(*,1) 'xh = ', xh
      write(*,1) 'xc = ', xc
      write(*,1) 'xo = ', xo
      write(*,1) 'lnfree_e = ', lnfree_e
      write(*,1) 'd_lnfree_e_dlnRho = ', d_lnfree_e_dlnRho
      write(*,1) 'd_lnfree_e_dlnT = ', d_lnfree_e_dlnT
      write(*,1) 'abar = ', s% abar(k)
      write(*,1) 'zbar = ', s% zbar(k)
      write(*,*)
      write(*,*) 'cubic_interpolation_in_X = ', s% kap_rq% cubic_interpolation_in_X
      write(*,*) 'cubic_interpolation_in_Z = ', s% kap_rq% cubic_interpolation_in_Z
      write(*,*) 'include_electron_conduction = ', s% kap_rq% include_electron_conduction
      write(*,*) 'use_Zbase_for_Type1 = ', s% kap_rq% use_Zbase_for_Type1
      write(*,*) 'use_Type2_opacities = ', s% kap_rq% use_Type2_opacities
      write(*,1) 'kap_Type2_full_off_X = ', s% kap_rq% kap_Type2_full_off_X
      write(*,1) 'kap_Type2_full_on_X = ', s% kap_rq% kap_Type2_full_on_X
      write(*,1) 'kap_Type2_full_off_dZ = ', s% kap_rq% kap_Type2_full_off_dZ
      write(*,1) 'kap_Type2_full_on_dZ = ', s% kap_rq% kap_Type2_full_on_dZ
      write(*,*)
      write(*,*)
      write(*,1) 'rho = ', s% rho(k)
      write(*,1) 'lnrho = ', s% lnd(k)
      write(*,1) 'T = ', s% T(k)
      write(*,1) 'lnT = ', s% lnT(k)
      write(*,1) 'logKap = ', log10(s% opacity(k))
      write(*,1) 'opacity = ', s% opacity(k)
      write(*,1) 'dlnkap_dlnd = ', dlnkap_dlnd
      write(*,1) 'dlnkap_dlnT = ', dlnkap_dlnT
      write(*,1) 'd_opacity_dlnd = ', s% d_opacity_dlnd(k)
      write(*,1) 'd_opacity_dlnT = ', s% d_opacity_dlnT(k)
      write(*,*)
      write(*,1) 'logQ = ', s% lnd(k)/ln10 - 2*s% lnT(k)/ln10 + 12
      write(*,*)
      write(*,*)
      write(*,1) 'kap_frac_Type2', s% kap_frac_Type2(k)
      write(*,1) 'extra_opacity_factor', s% extra_opacity_factor(k)
      write(*,*)
      write(*,*)
    end subroutine show_stuff

  end subroutine do_kap_for_cell

  !****

  subroutine shutdown_microphys
    use eos_lib, only: eos_shutdown
    use kap_lib, only: kap_shutdown
    use net_lib
    call eos_shutdown
    call kap_shutdown
    call net_shutdown
  end subroutine shutdown_microphys


end module micro

