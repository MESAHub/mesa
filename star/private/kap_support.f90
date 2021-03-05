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

module kap_support

  ! Uses

  use star_private_def
  use const_def
  use auto_diff

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: prepare_kap
  public :: setup_for_op_mono
  public :: fraction_of_op_mono
  public :: frac_op_mono
  public :: get_kap

  ! Procedures
  
contains

  subroutine prepare_kap (s, ierr)

    type(star_info), pointer :: s
    integer, intent(out)     :: ierr

    ! Set opacity parameters and monochromatic stuff

    call setup_for_op_mono(s, .true., ierr)
    if (ierr /= 0) return

  end subroutine prepare_kap

  !****

  subroutine setup_for_op_mono(s, check_if_need, ierr)

    use kap_lib, only: load_op_mono_data, get_op_mono_params
    use utils_lib, only: utils_OMP_GET_MAX_THREADS

    type (star_info), pointer :: s
    logical, intent(in)       :: check_if_need
    integer, intent(out)      :: ierr

    integer  :: k, nptot, ipe, nrad, n
    logical  :: need_op
    type(auto_diff_real_2var_order1) :: beta

    include 'formats'

    ierr = 0
    if (s% op_mono_n > 0) return ! already setup

    if (check_if_need) then
       if (s% high_logT_op_mono_full_off < 0d0 .or. &
            s% low_logT_op_mono_full_off > 99d0) return
       need_op = .false.
       do k=1,s% nz
          beta = fraction_of_op_mono(s, k)
          if (beta > 0d0) then
             need_op = .true.
             exit
          end if
       end do
       if (.not. need_op) return
    end if

    call load_op_mono_data( &
         s% op_mono_data_path, s% op_mono_data_cache_filename, ierr)
    if (ierr /= 0) then
       write(*,*) 'error while loading OP data, ierr = ',ierr
       return
    end if

    call get_op_mono_params(nptot, ipe, nrad)

    n = utils_OMP_GET_MAX_THREADS()

    if (n /= s% op_mono_n .or. &
         nptot /= s% op_mono_nptot .or. &
         ipe /= s% op_mono_ipe .or. &
         nrad /= s% op_mono_nrad) then
       if (associated(s% op_mono_umesh1)) deallocate(s% op_mono_umesh1)
       if (associated(s% op_mono_semesh1)) deallocate(s% op_mono_semesh1)
       if (associated(s% op_mono_ff1)) deallocate(s% op_mono_ff1)
       if (associated(s% op_mono_rs1)) deallocate(s% op_mono_rs1)
       if (s% use_op_mono_alt_get_kap) then
          allocate( &
             s% op_mono_umesh1(nptot*n), s% op_mono_semesh1(nptot*n), s% op_mono_ff1(nptot*ipe*6*6*n), &
             s% op_mono_rs1(nptot*6*6*n), stat=ierr)
       else
          allocate( &
             s% op_mono_umesh1(nptot*n), s% op_mono_semesh1(nptot*n), s% op_mono_ff1(nptot*ipe*4*4*n), &
             s% op_mono_rs1(nptot*4*4*n), stat=ierr)
       end if
       if (ierr /= 0) return
       s% op_mono_n = n
       s% op_mono_nptot = nptot
       s% op_mono_ipe = ipe
       s% op_mono_nrad = nrad
    end if

  end subroutine setup_for_op_mono

  !****

  real(dp) function fraction_of_op_mono(s, k) result(beta)
    ! returns the real(dp) value of the blend function for cell k
    type (star_info), pointer :: s
    integer, intent(in)       :: k

    type(auto_diff_real_2var_order1) :: frac

    frac = frac_op_mono(s, s% lnd(k)/ln10, s% lnT(k)/ln10)
    beta = frac% val

  end function fraction_of_op_mono

  !****

  type(auto_diff_real_2var_order1) function frac_op_mono(s, logRho, logT) result(beta)
    ! returns an auto_diff type: var1 = lnd, var2 = lnT (derivs w.r.t. ln *not* log)
    use utils_lib, only: mesa_error
    type (star_info), pointer :: s
    real(dp), intent(in)       :: logRho, logT

    real(dp) :: high_full_off, high_full_on, low_full_off, low_full_on
    type(auto_diff_real_2var_order1) :: alfa, log10_T

    include 'formats'

    ! default to 0
    beta = 0
    
    high_full_off = s% high_logT_op_mono_full_off
    high_full_on = s% high_logT_op_mono_full_on
    low_full_off = s% low_logT_op_mono_full_off
    low_full_on = s% low_logT_op_mono_full_on

    if (high_full_off < high_full_on &
          .or. high_full_on < low_full_on &
          .or. low_full_on < low_full_off) then
       write(*,*) "ERROR: OP mono controls must satisfy"
       write(*,*) "high_logT_op_mono_full_off >= high_logT_op_mono_full_on"
       write(*,*) "high_logT_op_mono_full_on >= low_logT_op_mono_full_on"
       write(*,*) "low_logT_op_mono_full_on >= high_logT_op_mono_full_off"
       call mesa_error(__FILE__,__LINE__,'fraction_of_op_mono')
    end if

    if (high_full_off <= 0d0 .or. high_full_on <= 0d0) then
       beta = 0._dp
       return
    end if

    ! auto_diff version of inputs
    log10_T = logT
    log10_T % d1val2 = 1d0/ln10

    ! alfa is fraction standard opacity

    if (log10_T >= high_full_off .or. log10_T <= low_full_off) then
       alfa = 1d0
    else if (log10_T <= high_full_on .and. log10_T >= low_full_on) then
       alfa = 0d0
    else if (log10_T > high_full_on) then ! between high_on and high_off
       if (high_full_off - high_full_on > 1d-10) then
          alfa = (log10_T - high_full_on) / (high_full_off - high_full_on)
       else
          alfa = 1d0
       end if
    else ! between low_off and low_on
       if (low_full_on - low_full_off > 1d-10) then
          alfa = (log10_T - low_full_on) / (low_full_off - low_full_on)
       else
          alfa = 1d0
       end if
    end if

    ! beta is fraction of op mono
    ! transform linear blend to smooth quintic blend
    alfa = -alfa*alfa*alfa*(-10d0 + alfa*(15d0 - 6d0*alfa))
    beta = 1d0 - alfa 

  end function frac_op_mono

  !****

  subroutine get_kap( &
       s, k, zbar, xa, logRho, logT, &
       lnfree_e, dlnfree_e_dlnRho, dlnfree_e_dlnT, &
       eta, d_eta_dlnRho, d_eta_dlnT, &
       kap_fracs, kap, dlnkap_dlnd, dlnkap_dlnT, ierr)

    use utils_lib
    use num_lib
    use kap_def, only: kap_test_partials, &
       kap_test_partials_val, kap_test_partials_dval_dx, &
       num_kap_fracs
    use kap_lib, only: &
         load_op_mono_data, get_op_mono_params, &
         get_op_mono_args, kap_get_op_mono, kap_get
    use chem_def, only: ih1, ihe3, ihe4, chem_isos
    use star_utils, only: get_XYZ, lookup_nameofvar

    type (star_info), pointer :: s
    integer, intent(in) :: k
    real(dp), intent(in) :: &
       zbar, logRho, logT
    real(dp), intent(in) :: lnfree_e, dlnfree_e_dlnRho, dlnfree_e_dlnT
    real(dp), intent(in) :: eta, d_eta_dlnRho, d_eta_dlnT
    real(dp), intent(in) :: xa(:)
    real(dp), intent(out) :: kap_fracs(num_kap_fracs)
    real(dp), intent(out) :: kap, dlnkap_dlnd, dlnkap_dlnT
    integer, intent(out) :: ierr

    real(dp) :: dlnkap_dxa(s% species)

    integer :: i, iz, nptot, ipe, nrad, thread_num, sz, offset
    type(auto_diff_real_2var_order1) :: beta, lnkap, lnkap_op
    real(dp) :: kap_op, dlnkap_op_dlnRho, dlnkap_op_dlnT, &
       Z, xh, xhe, xc, xn, xo, xne
    integer, pointer :: net_iso(:)
    real, pointer :: &
         umesh(:), semesh(:), ff(:,:,:,:), rs(:,:,:)
    integer :: nel, izzp(s% species)
    real(dp) :: fap(s% species), fac(s% species), gp1(s% species)
    logical :: screening
    real(dp), parameter :: &
         eps = 1d-6, minus_eps = -eps, one_plus_eps = 1d0 + eps

    integer :: i_var

    include 'formats'

    if (s% doing_timing) s% timing_num_get_kap_calls = s% timing_num_get_kap_calls + 1

    kap = -1d99
    if (k > 0 .and. k <= s% nz) s% kap_frac_op_mono(k) = 0

    net_iso => s% net_iso
    xc = 0; xn = 0; xo = 0; xne = 0
    do i=1, s% species
       iz = chem_isos% Z(s% chem_id(i))
       select case(iz)
       case (6)
          xc = xc + xa(i)
       case (7)
          xn = xn + xa(i)
       case (8)
          xo = xo + xa(i)
       case (10)
          xne = xne + xa(i)
       end select
    end do

    if (xc < minus_eps .or. xn < minus_eps .or. &
         xo < minus_eps .or. xne < minus_eps .or. &
         xc > one_plus_eps .or. xn > one_plus_eps .or. &
         xo > one_plus_eps .or. xne > one_plus_eps) then
       ierr = -1
       if (xc < 0d0 .or. xc > 1d0) write(s% retry_message,2) 'get_kap: xc', k, xc
       if (xn < 0d0 .or. xn > 1d0) write(s% retry_message,2) 'get_kap: xn', k, xn
       if (xo < 0d0 .or. xo > 1d0) write(s% retry_message,2) 'get_kap: xo', k, xo
       if (xne < 0d0 .or. xne > 1d0) write(s% retry_message,2) 'get_kap: xne', k, xne
       if (s% report_ierr) write(*, *) s% retry_message
       return
       stop 'bad mass fraction: get_kap'
    end if

    call get_XYZ(s, xa, xh, xhe, Z)

    if (xh < 0d0 .or. xhe < 0d0 .or. Z < 0d0) then
       ierr = -1
       if (xh < 0d0) write(s% retry_message,2) 'xh', k, xh
       if (xhe < 0d0) write(s% retry_message,2) 'xhe', k, xhe
       if (Z < 0d0) write(s% retry_message,2) 'Z', k, Z
       if (s% report_ierr) write(*, *) s% retry_message
       return
       stop 'negative mass fraction: get_kap'
    end if

    if (s% use_simple_es_for_kap) then
       kap = 0.2d0*(1 + xh)
       dlnkap_dlnd = 0
       dlnkap_dlnT = 0
       return
    end if

    beta = frac_op_mono(s, logRho, logT)
    if (k > 0 .and. k <= s% nz) s% kap_frac_op_mono(k) = beta % val

    if (beta > 0d0) then

       call get_op_mono_args( &
            s% species, xa, s% op_mono_min_X_to_include, s% chem_id, &
            s% op_mono_factors, nel, izzp, fap, fac, ierr)
       if (ierr /= 0) then
          write(*,*) 'error in get_op_mono_args, ierr = ',ierr
          return
       end if

       if (associated(s% op_mono_umesh1)) then

          thread_num = utils_OMP_GET_THREAD_NUM() ! in range 0 to op_mono_n-1
          if (thread_num < 0) then
             write(*,3) 'thread_num < 0', thread_num, s% op_mono_n
             ierr = -1
             return
          end if
          if (thread_num >= s% op_mono_n) then
             write(*,3) 'thread_num >= s% op_mono_n', thread_num, s% op_mono_n
             ierr = -1
             return
          end if
          nptot = s% op_mono_nptot
          ipe = s% op_mono_ipe
          nrad = s% op_mono_nrad

          sz = nptot; offset = thread_num*sz
          umesh(1:nptot) => s% op_mono_umesh1(offset+1:offset+sz)
          semesh(1:nptot) => s% op_mono_semesh1(offset+1:offset+sz)
          if (s% use_op_mono_alt_get_kap) then
             sz = nptot*ipe*6*6; offset = thread_num*sz
             ff(1:nptot,1:ipe,1:6,1:6) => s% op_mono_ff1(offset+1:offset+sz)
             sz = nptot*6*6; offset = thread_num*sz
             rs(1:nptot,1:6,1:6) => s% op_mono_rs1(offset+1:offset+sz)
             sz = nptot*nrad*6*6; offset = thread_num*sz
          else
             sz = nptot*ipe*4*4; offset = thread_num*sz
             ff(1:nptot,1:ipe,1:4,1:4) => s% op_mono_ff1(offset+1:offset+sz)
             sz = nptot*4*4; offset = thread_num*sz
             rs(1:nptot,1:4,1:4) => s% op_mono_rs1(offset+1:offset+sz)
             sz = nptot*nrad*4*4; offset = thread_num*sz
          end if

       else

          call load_op_mono_data( &
               s% op_mono_data_path, s% op_mono_data_cache_filename, ierr)
          if (ierr /= 0) then
             write(*,*) 'error while loading OP data, ierr = ',ierr
             return
          end if

          call get_op_mono_params(nptot, ipe, nrad)
          if (s% use_op_mono_alt_get_kap) then
             allocate( &
                umesh(nptot), semesh(nptot), ff(nptot,ipe,6,6), &
                rs(nptot,6,6), stat=ierr)
          else
             allocate( &
                umesh(nptot), semesh(nptot), ff(nptot,ipe,4,4), &
                rs(nptot,4,4), stat=ierr)
          end if
          if (ierr /= 0) return

       end if

       if (s% solver_test_kap_partials) then
          kap_test_partials = (k == s% solver_test_partials_k .and. &
             s% solver_call_number == s% solver_test_partials_call_number .and. &
             s% solver_iter == s% solver_test_partials_iter_number )
       end if

       screening = .true.
       if (s% use_other_kap) then
          call s% other_kap_get_op_mono( &
               s% kap_handle, zbar, logRho, logT, &
               s% use_op_mono_alt_get_kap, &
               nel, izzp, fap, fac, screening, umesh, semesh, ff, rs, &
               kap_op, dlnkap_op_dlnRho, dlnkap_op_dlnT, ierr)
       else
          call kap_get_op_mono( &
               s% kap_handle, zbar, logRho, logT, &
               s% use_op_mono_alt_get_kap, &
               nel, izzp, fap, fac, screening, umesh, semesh, ff, rs, &
               kap_op, dlnkap_op_dlnRho, dlnkap_op_dlnT, ierr)
       end if

       if (s% solver_test_kap_partials .and. kap_test_partials) then
          s% solver_test_partials_val = kap_test_partials_val
          s% solver_test_partials_dval_dx = kap_test_partials_dval_dx
       end if

       if (.not. associated(s% op_mono_umesh1)) deallocate(umesh, semesh, ff, rs)

       if (ierr /= 0) then
          s% retry_message = 'error in op_mono kap'
          if (s% report_ierr) write(*, *) s% retry_message
          beta = 0
          if (k > 0 .and. k <= s% nz) s% kap_frac_op_mono(k) = 0
          ierr = 0
       else if (beta == 1d0) then
          kap = kap_op
          dlnkap_dlnT = dlnkap_op_dlnT
          dlnkap_dlnd = dlnkap_op_dlnRho
          return
       end if

       lnkap_op = log(kap_op)
       lnkap_op% d1val1 = dlnkap_op_dlnRho
       lnkap_op% d1val2 = dlnkap_op_dlnT

    end if

    if (s% solver_test_kap_partials) then
       kap_test_partials = (k == s% solver_test_partials_k .and. &
          s% solver_call_number == s% solver_test_partials_call_number .and. &
          s% solver_iter == s% solver_test_partials_iter_number )
    end if

    if (s% use_other_kap) then
       call s% other_kap_get( &
            s% id, k, s% kap_handle, s% species, s% chem_id, s% net_iso, xa, &
            logRho, logT, &
            lnfree_e, dlnfree_e_dlnRho, dlnfree_e_dlnT, &
            eta, d_eta_dlnRho, d_eta_dlnT, &
            kap_fracs, kap, dlnkap_dlnd, dlnkap_dlnT, dlnkap_dxa, ierr)
    else
       call kap_get( &
            s% kap_handle, s% species, s% chem_id, s% net_iso, xa, &
            logRho, logT, &
            lnfree_e, dlnfree_e_dlnRho, dlnfree_e_dlnT, &
            eta, d_eta_dlnRho, d_eta_dlnT, &
            kap_fracs, kap, dlnkap_dlnd, dlnkap_dlnT, dlnkap_dxa, ierr)
    end if

    if (ierr /= 0) then
       return
    end if

    if (s% solver_test_kap_partials .and. kap_test_partials) then
       s% solver_test_partials_val = kap_test_partials_val
       s% solver_test_partials_dval_dx = kap_test_partials_dval_dx
    end if

    if (beta == 0d0) then
       return
    end if

    ! OP mono already packed in lnkap_op
    ! pack lnkap with standard opacities
    lnkap = log(kap)
    lnkap% d1val1 = dlnkap_dlnd
    lnkap% d1val2 = dlnkap_dlnT

    ! Blend opacities
    lnkap = (1d0-beta)*lnkap + beta*lnkap_op

    ! unpack autodiff type
    kap = exp(lnkap% val)
    dlnkap_dlnd = lnkap% d1val1
    dlnkap_dlnT = lnkap% d1val2

  end subroutine get_kap

end module kap_support
