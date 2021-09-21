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

module pulse_osc

  ! Uses

  use star_private_def
  use const_def
  use utils_lib
  use chem_def
  use atm_def
  use atm_support

  use pulse_utils

  ! No implicit typing

  implicit none

  ! Parameter definitions (these values are for OSC format version 2K,
  ! with 14 abundances)

  integer, parameter :: ICONST = 15
  integer, parameter :: IVAR = 22
  integer, parameter :: IABUND = 14

  ! Access specifiers

  private

  public :: get_osc_data
  public :: write_osc_data

contains

  subroutine get_osc_data (id, &
       add_center_point, keep_surface_point, add_atmosphere, global_data, point_data, ierr)

    integer, intent(in)                :: id
    logical, intent(in)                :: add_center_point
    logical, intent(in)                :: keep_surface_point
    logical, intent(in)                :: add_atmosphere
    real(dp), allocatable, intent(out) :: global_data(:)
    real(dp), allocatable, intent(out) :: point_data(:,:)
    integer, intent(out)               :: ierr

    type(star_info), pointer :: s
    integer, allocatable     :: k_a(:)
    integer, allocatable     :: k_b(:)
    integer                  :: n_sg
    integer                  :: h1
    integer                  :: h2
    integer                  :: he3
    integer                  :: he4
    integer                  :: li7
    integer                  :: be7
    integer                  :: be9
    integer                  :: c12
    integer                  :: c13
    integer                  :: n14
    integer                  :: n15
    integer                  :: o16
    integer                  :: o17
    integer                  :: o18
    integer                  :: si28
    integer                  :: n_atm
    integer                  :: nn_atm
    integer                  :: n_env
    integer                  :: nn_env
    integer                  :: nn
    real(dp)                 :: r_outer
    real(dp)                 :: m_outer
    real(dp)                 :: P_c
    real(dp)                 :: rho_c
    real(dp)                 :: d2P_dr2_c
    integer                  :: j
    integer                  :: k
    integer                  :: sg

    ! Get model data for OSC output

    call get_star_ptr(id, s, ierr)
    if (ierr /= 0) then
       write(*,*) 'bad star id for get_osc_data'
       return
    end if

    ! Set up segment indices

    call set_segment_indices(s, k_a, k_b, add_center_point)

    n_sg = SIZE(k_a)

    ! Set up specicies indices

    h1 = s%net_iso(ih1)
    h2 = s%net_iso(ih2)
    he3 = s%net_iso(ihe3)
    he4 = s%net_iso(ihe4)
    li7 = s%net_iso(ili7)
    be7 = s%net_iso(ibe7)
    be9 = s%net_iso(ibe9)
    c12 = s%net_iso(ic12)
    c13 = s%net_iso(ic13)
    n14 = s%net_iso(in14)
    n15 = s%net_iso(in15)
    o16 = s%net_iso(io16)
    o17 = s%net_iso(io17)
    o18 = s%net_iso(io18)
    si28 = s%net_iso(isi28)

    ! Determine data dimensiones

    if (add_atmosphere) then
       call build_atm(s, s%L(1), s%r(1), s%m_grav(1), s%cgrav(1), ierr)
       if (ierr /= 0) then
          write(*,*) 'failed in build_atm'
          return
       end if
       n_atm = s%atm_structure_num_pts
       nn_atm = n_atm - 1
    else
       n_atm = 0
       nn_atm = 0
    end if

    n_env = s%nz

    if (keep_surface_point) then
       nn_env = n_env + n_sg - 1
    else
       nn_env = n_env - 1 + n_sg - 1
    endif

    if (add_center_point) then
       nn = nn_env + nn_atm + 1
    else
       nn = nn_env + nn_atm
    endif

    ! Store global data

    allocate(global_data(ICONST))
    global_data = 0d0

    r_outer = Rsun*s%photosphere_r
    m_outer = s%m_grav(1)

    global_data(1) = m_outer
    global_data(2) = r_outer
    global_data(3) = s%L(1)
    global_data(4) = s%initial_z
    global_data(5) = 1d0 - (s%initial_y + s%initial_z)
    global_data(6) = s%mixing_length_alpha
    if (s%largest_conv_mixing_region /= 0) then
       k = s%mixing_region_bottom(s%largest_conv_mixing_region)
       global_data(7) = s%X(k)
       global_data(8) = s%Y(k)
    end if

    if (s%interpolate_rho_for_pulse_data) then
       rho_c = eval_center_rho(s, k_b(n_sg))
    else
       rho_c = eval_center(s%rmid, s%rho, k_a(n_sg), k_b(n_sg))
    endif

    ! at the centre d²P/dr² = -4πGρ²/3
    d2P_dr2_c = -four_thirds*pi*s% cgrav(s% nz)*rho_c**2
    P_c = s%Peos(s% nz) - 0.5*d2P_dr2_c*s% rmid(s% nz)**2
    global_data(9) = r_outer**2*d2P_dr2_c/P_c
    global_data(10) = r_outer**2*eval_center_d2(s%rmid, s%rho, k_a(n_sg), k_b(n_sg)) / rho_c

    global_data(11) = s%star_age
    if (s%rotation_flag) then
       global_data(12) = s%omega(1)
    else
       global_data(12) = 0d0
    end if

    ! global_data(13) should be the initial rotation rate, but we lack
    ! that datum
    
    ! Store point data

    allocate(point_data(IVAR+IABUND,nn))
    point_data = 0d0
    
    j = 1

    ! Atmosphere (we skip the point at the base of the atm to
    ! smooth the transition)

    atm_loop : do k = 1, n_atm-1
       call store_point_data_atm(j, k, k_a(1), k_b(1))
       j = j + 1
    end do atm_loop
    
    ! Envelope

    sg = 1

    env_loop : do k = 1, n_env

       if (k == 1 .AND. .NOT. keep_surface_point) cycle env_loop

       call store_point_data_env(j, k, k_a(sg), k_b(sg))
       j = j + 1

       if (k == k_b(sg)+1) then

          sg = sg + 1

          call store_point_data_env(j, k, k_a(sg), k_b(sg))
          j = j + 1
             
       endif

    end do env_loop

    ! Center

    if (add_center_point) then
       call store_point_data_ctr(j, k_a(n_sg), k_b(n_sg))
       j = j + 1
    end if

    ! Tidy up

    if (ASSOCIATED(s%atm_structure)) then
       deallocate(s%atm_structure)
    end if

    ! Finish

    return

  contains

    subroutine store_point_data_atm (j, k, k_a, k_b)

      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: k_a
      integer, intent(in) :: k_b

      real(dp) :: grav
      real(dp) :: N2

      ! Store data associated with atmosphere point k into the point_data array at
      ! position j

      associate ( &
           r => point_data(1,j), &
           lnq => point_data(2,j), &
           T => point_data(3,j), &
           P => point_data(4,j), &
           rho => point_data(5,j), &
           nabla => point_data(6,j), &
           L => point_data(7,j), &
           kap => point_data(8,j), &
           eps => point_data(9,j), &
           Gamma_1 => point_data(10,j), &
           nabla_ad => point_data(11,j), &
           delta => point_data(12,j), &
           c_P => point_data(13, j), &
           rec_mu_e => point_data(14,j), &
           A_ast => point_data(15,j), &
           omega => point_data(16,j), &
           kap_T => point_data(17,j), &
           kap_rho => point_data(18,j), &
           eps_T => point_data(19,j), &
           eps_rho => point_data(20,j), &
           P_tot_gas => point_data(21,j), &
           nabla_rad => point_data(22,j), &
           X_H1 => point_data(23,j), &
           X_H2 => point_data(24,j), &
           X_He3 => point_data(25,j), &
           X_He4 => point_data(26,j), &
           X_Li7 => point_data(27,j), &
           X_Be7 => point_data(28,j), &
           X_C12 => point_data(29,j), &
           X_C13 => point_data(30,j), &
           X_N14=> point_data(31,j), &
           X_N15 => point_data(32,j), &
           X_O16 => point_data(33,j), &
           X_O17 => point_data(34,j), &
           X_Be9 => point_data(35,j), &
           X_Si28 => point_data(36,j))
      
        r = s%r(1) + s%atm_structure(atm_delta_r,k)
        lnq = log(s%m_grav(1)/m_outer)
        T = exp(s%atm_structure(atm_lnT,k))
        P = exp(s%atm_structure(atm_lnP,k))
        rho = exp(s%atm_structure(atm_lnd,k))
        nabla = s%atm_structure(atm_gradT,k)
        L = s%L(1)
        kap = s%atm_structure(atm_kap,k)
        eps = 0d0
        Gamma_1 = s%atm_structure(atm_gamma1,k)
        nabla_ad = s%atm_structure(atm_grada,k)
        delta = s%atm_structure(atm_chiT,k)/s%atm_structure(atm_chiRho,k)
        c_P = s%atm_structure(atm_cp,k)
        rec_mu_e = exp(s%atm_structure(atm_lnfree_e,k)) ! check

        grav = s%cgrav(1)*s%m_grav(1)/(r*r)
        N2 = grav*grav*(rho/P)*delta*(nabla_ad - nabla)
        A_ast = N2*r/grav

        if (s%rotation_flag) then
           omega = s%omega(1)
        else
           omega = 0d0
        endif
        kap_T = s%atm_structure(atm_dlnkap_dlnT,k)
        kap_rho = s%atm_structure(atm_dlnkap_dlnd,k)
        eps_T = 0d0
        eps_rho = 0d0
        P_tot_gas = P/exp(s%atm_structure(atm_lnPgas,k))
        nabla_rad = s%atm_structure(atm_gradr,k)
        X_H1 = eval_face_X(s, h1, 1, k_a, k_b)
        X_H2 = eval_face_X(s, h2, 1, k_a, k_b)
        X_He3 = eval_face_X(s, he3, 1, k_a, k_b)
        X_He4 = eval_face_X(s, he4, 1, k_a, k_b)
        X_Li7 = eval_face_X(s, li7, 1, k_a, k_b)
        X_Be7 = eval_face_X(s, be7, 1, k_a, k_b)
        X_C12 = eval_face_X(s, c12, 1, k_a, k_b)
        X_C13 = eval_face_X(s, c13, 1, k_a, k_b)
        X_N14 = eval_face_X(s, n14, 1, k_a, k_b)
        X_N15 = eval_face_X(s, n15, 1, k_a, k_b)
        X_O16 = eval_face_X(s, o16, 1, k_a, k_b)
        X_O17 = eval_face_X(s, o17, 1, k_a, k_b)
        X_Be9 = eval_face_X(s, be9, 1, k_a, k_b)
        X_Si28 = eval_face_X(s, si28, 1, k_a, k_b)

      end associate

      ! Finish

      return

    end subroutine store_point_data_atm

    !****

    subroutine store_point_data_env (j, k, k_a, k_b)

      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: k_a
      integer, intent(in) :: k_b

      ! Store data associated with envelope face k into the point_data array at
      ! position j

      associate ( &
           r => point_data(1,j), &
           lnq => point_data(2,j), &
           T => point_data(3,j), &
           P => point_data(4,j), &
           rho => point_data(5,j), &
           nabla => point_data(6,j), &
           L => point_data(7,j), &
           kap => point_data(8,j), &
           eps => point_data(9,j), &
           Gamma_1 => point_data(10,j), &
           nabla_ad => point_data(11,j), &
           delta => point_data(12,j), &
           c_P => point_data(13, j), &
           rec_mu_e => point_data(14,j), &
           A_ast => point_data(15,j), &
           omega => point_data(16,j), &
           kap_T => point_data(17,j), &
           kap_rho => point_data(18,j), &
           eps_T => point_data(19,j), &
           eps_rho => point_data(20,j), &
           P_tot_gas => point_data(21,j), &
           nabla_rad => point_data(22,j), &
           X_H1 => point_data(23,j), &
           X_H2 => point_data(24,j), &
           X_He3 => point_data(25,j), &
           X_He4 => point_data(26,j), &
           X_Li7 => point_data(27,j), &
           X_Be7 => point_data(28,j), &
           X_C12 => point_data(29,j), &
           X_C13 => point_data(30,j), &
           X_N14=> point_data(31,j), &
           X_N15 => point_data(32,j), &
           X_O16 => point_data(33,j), &
           X_O17 => point_data(34,j), &
           X_Be9 => point_data(35,j), &
           X_Si28 => point_data(36,j))
      
        r = s%r(k)
        lnq = log(s%m_grav(k)/m_outer)
        T = eval_face(s%dq, s%T, k, 1, s%nz)
        P = eval_face(s%dq, s%Peos, k, 1, s%nz)
        if (s%interpolate_rho_for_pulse_data) then
           rho = eval_face(s%dq, s%rho, k, k_a, k_b)
        else
           rho = eval_face_rho(s, k, k_a, k_b)
        endif
        nabla = s%gradT(k) ! Not quite right; gradT can be discontinuous
        L = s%L(k)
        kap = eval_face(s%dq, s%opacity, k, k_a, k_b)
        eps = eval_face(s%dq, s%eps_nuc, k, k_a, k_b) + eval_face(s%dq, s%eps_grav_ad%val, k, k_a, k_b)
        Gamma_1 = eval_face(s%dq, s%gamma1, k, k_a, k_b)
        nabla_ad = eval_face(s%dq, s%grada, k, k_a, k_b)
        delta = eval_face(s%dq, s%chiT, k, k_a, k_b)/eval_face(s%dq, s%chiRho, k, k_a, k_b)
        c_P = eval_face(s%dq, s%cp, k, k_a, k_b)
        rec_mu_e = exp(eval_face(s%dq, s%lnfree_e, k, k_a, k_b)) ! check
        A_ast = eval_face_A_ast(s, k, k_a, k_b)
        if (s%rotation_flag) then
           omega = s%omega(k) ! Not quite right; omega can be discontinuous
        else
           omega = 0d0
        endif
        kap_T = eval_face(s%dq, s%d_opacity_dlnT, k, k_a, k_b)/kap
        kap_rho = eval_face(s%dq, s%d_opacity_dlnd, k, k_a, k_b)/kap
        eps_T = eval_face(s%dq, s%d_epsnuc_dlnT, k, k_a, k_b)
        eps_rho = eval_face(s%dq, s%d_epsnuc_dlnd, k, k_a, k_b)
        P_tot_gas = eval_face(s%dq, s%Peos, k, 1, s%nz)/eval_face(s%dq, s%Pgas, k, k_a, k_b)
        nabla_rad = s%gradr(k)  ! Not quite right; gradr can be discontinuous
        X_H1 = eval_face_X(s, h1, k, k_a, k_b)
        X_H2 = eval_face_X(s, h2, k, k_a, k_b)
        X_He3 = eval_face_X(s, he3, k, k_a, k_b)
        X_He4 = eval_face_X(s, he4, k, k_a, k_b)
        X_Li7 = eval_face_X(s, li7, k, k_a, k_b)
        X_Be7 = eval_face_X(s, be7, k, k_a, k_b)
        X_C12 = eval_face_X(s, c12, k, k_a, k_b)
        X_C13 = eval_face_X(s, c13, k, k_a, k_b)
        X_N14 = eval_face_X(s, n14, k, k_a, k_b)
        X_N15 = eval_face_X(s, n15, k, k_a, k_b)
        X_O16 = eval_face_X(s, o16, k, k_a, k_b)
        X_O17 = eval_face_X(s, o17, k, k_a, k_b)
        X_Be9 = eval_face_X(s, be9, k, k_a, k_b)
        X_Si28 = eval_face_X(s, si28, k, k_a, k_b)

      end associate

      ! Finish

      return

    end subroutine store_point_data_env

    !****

    subroutine store_point_data_ctr (j, k_a, k_b)

      integer, intent(in) :: j
      integer, intent(in) :: k_a
      integer, intent(in) :: k_b

      ! Store data for the center into the point_data array at position j

      associate ( &
           r => point_data(1,j), &
           lnq => point_data(2,j), &
           T => point_data(3,j), &
           P => point_data(4,j), &
           rho => point_data(5,j), &
           nabla => point_data(6,j), &
           L => point_data(7,j), &
           kap => point_data(8,j), &
           eps => point_data(9,j), &
           Gamma_1 => point_data(10,j), &
           nabla_ad => point_data(11,j), &
           delta => point_data(12,j), &
           c_P => point_data(13, j), &
           rec_mu_e => point_data(14,j), &
           A_ast => point_data(15,j), &
           omega => point_data(16,j), &
           kap_T => point_data(17,j), &
           kap_rho => point_data(18,j), &
           eps_T => point_data(19,j), &
           eps_rho => point_data(20,j), &
           P_tot_gas => point_data(21,j), &
           nabla_rad => point_data(22,j), &
           X_H1 => point_data(23,j), &
           X_H2 => point_data(24,j), &
           X_He3 => point_data(25,j), &
           X_He4 => point_data(26,j), &
           X_Li7 => point_data(27,j), &
           X_Be7 => point_data(28,j), &
           X_C12 => point_data(29,j), &
           X_C13 => point_data(30,j), &
           X_N14=> point_data(31,j), &
           X_N15 => point_data(32,j), &
           X_O16 => point_data(33,j), &
           X_O17 => point_data(34,j), &
           X_Be9 => point_data(35,j), &
           X_Si28 => point_data(36,j))

        r = 0d0
        lnq = log(TINY(0d0))
        T = eval_center(s%rmid, s%T, 1, s%nz)
        P = P_c
        rho = rho_c
        nabla = eval_center(s%r, s%gradT, k_a, k_b)
        L = 0d0
        kap = eval_center(s%rmid, s%opacity, k_a, k_b)
        eps = eval_center(s%rmid, s%eps_nuc, k_a, k_b) + eval_center(s%rmid, s%eps_grav_ad%val, k_a, k_b)
        Gamma_1 = eval_center(s%rmid, s%gamma1, k_a, k_b)
        nabla_ad = eval_center(s%rmid, s%grada, k_a, k_b)
        delta = eval_center(s%rmid, s%chiT, k_a, k_b)/eval_center(s%rmid, s%chiRho, k_a, k_b)
        c_P = eval_center(s%rmid, s%cp, k_a, k_b)
        rec_mu_e = exp(eval_center(s%rmid, s%lnfree_e, k_a, k_b)) ! check
        A_ast = point_data(15,j)
        if (s%rotation_flag) then
           omega = eval_center(s%r, s%omega, k_a, k_b)
        else
           omega = 0d0
        endif
        kap_T = eval_center(s%rmid, s%d_opacity_dlnT, k_a, k_b)/kap
        kap_rho = eval_center(s%rmid, s%d_opacity_dlnd, k_a, k_b)/kap
        eps_T = eval_center(s%rmid, s%d_epsnuc_dlnT, k_a, k_b)
        eps_rho = eval_center(s%rmid, s%d_epsnuc_dlnd, k_a, k_b)
        P_tot_gas = eval_center(s%rmid, s%Peos, 1, s%nz)/eval_center(s%rmid, s%Pgas, k_a, k_b)
        nabla_rad = eval_center(s%r, s%gradr, k_a, k_b)
        X_H1 = eval_center_X(s, h1, k_a, k_b)
        X_H2 = eval_center_X(s, h2, k_a, k_b)
        X_He3 = eval_center_X(s, he3, k_a, k_b)
        X_He4 = eval_center_X(s, he4, k_a, k_b)
        X_Li7 = eval_center_X(s, li7, k_a, k_b)
        X_Be7 = eval_center_X(s, be7, k_a, k_b)
        X_C12 = eval_center_X(s, c12, k_a, k_b)
        X_C13 = eval_center_X(s, c13, k_a, k_b)
        X_N14 = eval_center_X(s, n14, k_a, k_b)
        X_N15 = eval_center_X(s, n15, k_a, k_b)
        X_O16 = eval_center_X(s, o16, k_a, k_b)
        X_O17 = eval_center_X(s, o17, k_a, k_b)
        X_Be9 = eval_center_X(s, be9, k_a, k_b)
        X_Si28 = eval_center_X(s, si28, k_a, k_b)

      end associate

      ! Finish

      return

    end subroutine store_point_data_ctr

  end subroutine get_osc_data

  !****

  subroutine write_osc_data (id, filename, global_data, point_data, ierr)

    integer, intent(in)      :: id
    character(*), intent(in) :: filename
    real(dp), intent(in)     :: global_data(:)
    real(dp), intent(in)     :: point_data(:,:)
    integer, intent(out)     :: ierr

    type(star_info), pointer :: s
    integer                  :: iounit
    integer                  :: n_global
    integer                  :: n_point
    integer                  :: nn
    integer                  :: i
    integer                  :: j
    integer                  :: k

    ! Write OSC data to file

    call get_star_ptr(id, s, ierr)
    if (ierr /= 0) then
       write(*,*) 'bad star id for write_osc_data'
       return
    end if

    ! Open the file

    open(newunit=iounit, file=TRIM(filename), status='REPLACE', iostat=ierr)
    if (ierr /= 0) then
       write(*,*) 'failed to open '//TRIM(filename)
       return
    end if

    ! Write the data

    n_global = SIZE(global_data)
    n_point = SIZE(point_data, 1)

    nn = SIZE(point_data, 2)

    write(iounit, *) 'OSC file'
    write(iounit, *) 'Created by MESAstar version', version_number
    write(iounit, *)
    write(iounit, *)

    write(iounit, '(a)') '14 H1 H2 He3 He4 Li7 Be7 C12 C13 N14 N15 O16 O17 Be9 Si28'
    write(iounit, '(5I10)') nn, ICONST, IVAR, IABUND, 2000

    write(iounit, s%format_for_osc_data) (global_data(i), i=1,n_global)

    do k = 1, nn
       write(iounit, s%format_for_osc_data) (point_data(j,k), j=1,n_point)
    end do

    ! Close the file

    close(iounit)

    ! Finish

    return

  end subroutine write_osc_data

end module pulse_osc
