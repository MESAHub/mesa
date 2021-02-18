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

module test_ecapture

  use const_def
  use eos_def
  use eos_lib
  use rates_def
  use rates_lib
  use num_lib
  use math_lib
  use utils_lib, only: mesa_error

  implicit none

  real(dp) :: X, Z, Y, abar, zbar, z2bar, z53bar, ye
  integer, parameter :: species = 7
  integer, parameter :: h1=1, he4=2, c12=3, n14=4, o16=5, ne20=6, mg24=7
  integer, pointer, dimension(:) :: net_iso, chem_id
  real(dp), dimension(species) :: xa, ya, za, aa, za52
  integer :: handle, i, ierr

  real(dp) :: log10T, T, log10Rho, Rho

contains


  subroutine do_test_ecapture

    call Init_Composition
    call Setup_eos(handle)

    log10Rho = 9.5d0
    Rho = exp10(log10Rho)
    log10T = 8.5d0
    T = exp10(log10T)

    ! first test the phase space integral functions

    write(*,*)
    write(*,*) 'do_test_psi_derivs'

    call do_test_psi_derivs

    write(*,*) 'done'
    write(*,*)

    ! check that the coulomb corrections are behaving

    write(*,*)
    write(*,*) 'do_test_coulomb'

    call do_test_coulomb

    write(*,*) 'done'
    write(*,*)


    ! check that the special weak reactions are working

    write(*,*)
    write(*,*) 'do_test_special_weak'

    call do_test_special_weak(.false.)
    call do_test_special_weak(.true.)

    write(*,*) 'done'
    write(*,*)

    ! deallocate the eos tables
    call free_eos_handle(handle)
    call eos_shutdown

  end subroutine do_test_ecapture


  subroutine do_test_psi_derivs

    real(dp) :: q
    logical :: doing_d_dlnd, doing_I

    real(dp) :: dvardx, dvardx_0, dx_0, err, var_0, xdum
    real(dp) :: IJ, dIJ_dlnT, dIJ_dlnRho

    include 'formats'

    ! you only need this for debugging
    write(*,*) '   intentionally skipping test'
    return

    q = -7d0

    doing_d_dlnd = .true.
    doing_I = .true.

    call I_or_J_wrapper(log10T, log10Rho, IJ, dIJ_dlnT, dIJ_dlnRho)
    var_0 = IJ
    write(*,*) -q, var_0

    if (doing_d_dlnd) then
       dx_0 = max(1d-14, abs(log10Rho*ln10*1d-6))
       dvardx_0 = dIJ_dlnRho
    else
       dx_0 = max(1d-14, abs(log10T*ln10*1d-6))
       dvardx_0 = dIJ_dlnT
    end if
    err = 0d0
    dvardx = dfridr(dx_0,dfridr_IJ,err)
    xdum = (dvardx - dvardx_0)/max(abs(dvardx_0),1d-50)
    write(*,1) 'analytic, numeric, est err in numeric, rel diff', &
         dvardx_0, dvardx, err, xdum
    if (doing_d_dlnd) then
       write(*,*) 'doing dlnd'
    else ! doing d_dlnT
       write(*,*) 'doing dlnT'
    end if
    write(*,*) 'test psi derivs'
    write(*,*)

  contains

    subroutine I_or_J_wrapper(logT, logRho, IJ, dIJ_dlnT, dIJ_dlnRho)

      real(dp), intent(in) :: logT, logRho
      real(dp), intent(out) :: IJ, dIJ_dlnT, dIJ_dlnRho

      real(dp) :: D9, T8

      real(dp) :: beta   ! mec2 / kT
      real(dp) :: zeta   ! Q_n / kT
      real(dp) :: mu     ! chemical potential
      real(dp) :: eta    ! chemical potential / kT
      real(dp) :: deta_dlnT, deta_dlnRho ! and derivs

      real(dp) :: I, J    ! phase space integral
      real(dp) :: dI_dlnT, dI_dlnRho ! and derivatives
      real(dp) :: dJ_dlnT, dJ_dlnRho ! and derivatives

      T8 = exp10(logT-8d0)
      D9 = exp10(logRho-9d0)

      ! effectively assumes Ye = 0.5
      mu = 10.1d0 * pow(D9,1d0/3d0)
      beta = 59.3d0 / T8

      zeta = q * beta
      eta = mu * beta
      deta_dlnRho = eta / 3.0d0
      deta_dlnT = -eta

      if (q < 0) then
         call psi_Iec_and_Jec(beta, zeta, eta, deta_dlnT, deta_dlnRho, &
              I, dI_dlnT, dI_dlnRho, J, dJ_dlnT, dJ_dlnRho)
      else
         call psi_Iee_and_Jee(beta, zeta, eta, deta_dlnT, deta_dlnRho, &
              I, dI_dlnT, dI_dlnRho, J, dJ_dlnT, dJ_dlnRho)
      end if

      if (doing_I) then
         IJ = I
         dIJ_dlnT = dI_dlnT
         dIJ_dlnRho = dI_dlnRho
      else
         IJ = J
         dIJ_dlnT = dJ_dlnT
         dIJ_dlnRho = dJ_dlnRho
      end if

    end subroutine I_or_J_wrapper


    real(dp) function dfridr_IJ(delta_x) result(val)
      real(dp), intent(in) :: delta_x
      real(dp) :: log_var

      include 'formats'
      ierr = 0

      if (doing_d_dlnd) then

         log_var = log10Rho + delta_x/ln10
         call I_or_J_wrapper(log10T, log_var, IJ, dIJ_dlnT, dIJ_dlnRho)

      else

         log_var = log10T + delta_x/ln10
         call I_or_J_wrapper(log_var, log10Rho, IJ, dIJ_dlnT, dIJ_dlnRho)

      end if

      val = IJ

    end function dfridr_IJ

  end subroutine do_test_psi_derivs



  subroutine do_test_coulomb

    use eval_coulomb
    use rates_def, only: Coulomb_Info, which_mui_coulomb, which_vs_coulomb

    real(dp) :: mu1, mu2, vs
    real(dp) :: z1, z2
    type(Coulomb_Info), pointer :: cc

    include 'formats'

    which_mui_coulomb = PCR2009
    which_vs_coulomb = Itoh2002

    z1 = 12
    z2 = 11

    allocate(cc)

    call coulomb_set_context(cc, T, Rho, log10T, log10Rho, &
         zbar, abar, z2bar, species, ya, za52)

    mu1 = do_mui_coulomb(cc, z1)
    mu2 = do_mui_coulomb(cc, z2)

    vs = do_vs_coulomb(cc, z1)

    deallocate(cc)

    write(*,'(6X, A4, 3F26.16)') 'mu', mu1, mu2, abs(mu1-mu2) * kev * T
    write(*,'(6X, A4, F26.16)') 'vs', vs

  end subroutine do_test_coulomb


  subroutine do_test_special_weak(use_special)

    use const_lib
    use utils_lib
    use eos_def
    use eos_lib
    use chem_def
    use chem_lib
    use rates_def, only: Coulomb_Info
    use eval_coulomb

    logical, intent(in) :: use_special

    real(dp) :: Rho, T, Pgas, log10Rho, log10T
    real(dp) :: dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, d_dlnRho_const_T, d_dlnT_const_Rho
    real(dp), dimension(num_eos_basic_results) :: res, d_dlnd, d_dlnT, d_dabar, d_dzbar
    integer :: ierr

    integer :: i, ir, nr
    integer, pointer :: ids(:), reaction_ids(:)
    type(Coulomb_Info), pointer :: cc

    real(dp), dimension(:), pointer :: &
         lambda, dlambda_dlnT, dlambda_dlnRho, &
         Q, dQ_dlnT, dQ_dlnRho, &
         Qneu, dQneu_dlnT, dQneu_dlnRho
    real(dp) :: lntwo, logT, T9, dT9, dlnT, YeRho, &
         ye, logRho, dlogRho, eta, d_eta_dlnT, d_eta_dlnRho, &
         Prad, energy, entropy
    character(len=iso_name_length), dimension(2) :: weak_lhs, weak_rhs

    allocate(cc)

    nr = 2

    allocate( &
         ids(nr), reaction_ids(nr), &
         lambda(nr), dlambda_dlnT(nr), dlambda_dlnRho(nr), &
         Q(nr), dQ_dlnT(nr), dQ_dlnRho(nr), &
         Qneu(nr), dQneu_dlnT(nr), dQneu_dlnRho(nr), &
         stat=ierr)


    ! pick reactions
    weak_lhs(1) = 'mg24'
    weak_rhs(1) = 'na24'

    weak_lhs(2) = 'na24'
    weak_rhs(2) = 'mg24'

    do i = 1, nr
       ids(i) = get_weak_rate_id(weak_lhs(i), weak_rhs(i))
       reaction_ids(i) = 0
    enddo

    logT = 8.3d0
    logRho = 9.8d0
    Ye = 0.5d0

    T = exp10(logT)
    T9 = T*1d-9
    rho = exp10(logRho)
    YeRho = Ye*rho

    ! get a set of results for given temperature and density
    call eosDT_get_legacy( &
         handle, Z, X, abar, zbar, &
         species, chem_id, net_iso, xa, &
         Rho, logRho, T, logT,   &
         !res, d_dlnd, d_dlnT, Pgas, Prad, energy, entropy, ierr)
         res, d_dlnd, d_dlnT, d_dabar, d_dzbar, ierr)

    call coulomb_set_context(cc, T, Rho, log10T, log10Rho, &
         zbar, abar, z2bar, species, ya, za)

    eta = res(i_eta)
    d_eta_dlnT = d_dlnT(i_eta)
    d_eta_dlnRho = d_dlnd(i_eta)


    if (use_special) then

       do_ecapture = .true.
       ecapture_states_file = 'test_special.states'
       ecapture_transitions_file = 'test_special.transitions'
       which_mui_coulomb = PCR2009
       which_vs_coulomb = Itoh2002

       call eval_ecapture_reaction_info( &
            nr, ids, cc, T9, YeRho, &
            eta, d_eta_dlnT, d_eta_dlnRho, &
            lambda, dlambda_dlnT, dlambda_dlnRho, &
            Q, dQ_dlnT, dQ_dlnRho, &
            Qneu, dQneu_dlnT, dQneu_dlnRho, &
            ierr)

       write(*,*) "special weak rates"
    else
       do_ecapture = .false.
       call eval_weak_reaction_info( &
            nr, ids, reaction_ids, cc, T9, YeRho, &
            eta, d_eta_dlnT, d_eta_dlnRho, &
            lambda, dlambda_dlnT, dlambda_dlnRho, &
            Q, dQ_dlnT, dQ_dlnRho, &
            Qneu, dQneu_dlnT, dQneu_dlnRho, &
            ierr)

       write(*,*) "weaklib weak rates"

    end if

    do i = 1, nr
       write(*,'(6X, 2A6, ES26.16)') weak_lhs(i), weak_rhs(i), lambda(i)
    enddo
    write(*,*)

    deallocate( &
         ids, reaction_ids, &
         lambda, dlambda_dlnT, dlambda_dlnRho, &
         Q, dQ_dlnT, dQ_dlnRho, &
         Qneu, dQneu_dlnT, dQneu_dlnRho, cc)


  end subroutine do_test_special_weak

  subroutine Init_Composition
    use chem_def
    use chem_lib

    real(dp) :: frac, dabar_dx(species), dzbar_dx(species),   &
         sumx, xh, xhe, xz, mass_correction, dmc_dx(species)

    allocate(net_iso(num_chem_isos), chem_id(species), stat=ierr)
    if (ierr /= 0) stop 'allocate failed'
    X = 0.70d0
    Z = 0.02d0

    net_iso(:) = 0

    chem_id(h1) = ih1; net_iso(ih1) = h1
    chem_id(he4) = ihe4; net_iso(ihe4) = he4
    chem_id(c12) = ic12; net_iso(ic12) = c12
    chem_id(n14) = in14; net_iso(in14) = n14
    chem_id(o16) = io16; net_iso(io16) = o16
    chem_id(ne20) = ine20; net_iso(ine20) = ne20
    chem_id(mg24) = img24; net_iso(img24) = mg24

    xa(h1) = 0.0d0
    xa(he4) = 0.0d0
    xa(c12) = 0.0d0
    xa(n14) = 0.0d0
    xa(o16) = 0.0d0
    xa(ne20) = 0.5d0
    xa(mg24) = 0.5d0

    do i = 1, species
       za(i) = chem_isos% Z(chem_id(i))
       aa(i) = chem_isos% W(chem_id(i))
       ya(i) = xa(i) / aa(i)
       za52(i) = pow(real(chem_isos% Z(chem_id(i)),kind=dp),5.0d0/2.d0)
    enddo

    call composition_info( &
         species, chem_id, xa, xh, xhe, xz, abar, zbar, z2bar, z53bar, ye,   &
         mass_correction, sumx, dabar_dx, dzbar_dx, dmc_dx)

  end subroutine Init_Composition


  subroutine Setup_eos(handle)
    ! allocate and load the eos tables
    use eos_def
    use eos_lib
    integer, intent(out) :: handle

    integer :: ierr
    logical, parameter :: use_cache = .true.

    call eos_init('', '', '', use_cache, ierr)
    if (ierr /= 0) then
       write(*,*) 'eos_init failed in Setup_eos'
       call mesa_error(__FILE__,__LINE__)
    end if

    handle = alloc_eos_handle(ierr)
    if (ierr /= 0) then
       write(*,*) 'failed trying to allocate eos handle'
       call mesa_error(__FILE__,__LINE__)
    end if

  end subroutine Setup_eos



  

  
end module test_ecapture
