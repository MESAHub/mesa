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

module pulse_fgong

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

  ! Parameter definitions (these values are for FGONG format verions
  ! 300 & 1300)

  integer, parameter :: ICONST = 15
  integer, parameter :: IVAR = 40

  ! Access specifiers

  private

  public :: get_fgong_data
  public :: write_fgong_data

contains

  subroutine get_fgong_data (id, &
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
    integer                  :: c12
    integer                  :: c13
    integer                  :: n14
    integer                  :: n15
    integer                  :: o16
    integer                  :: o17
    integer                  :: o18
    integer                  :: ne20
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

    ! Get FGONG data

    call get_star_ptr(id, s, ierr)
    if (ierr /= 0) then
       write(*,*) 'bad star id for get_fgong_data'
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
    c12 = s%net_iso(ic12)
    c13 = s%net_iso(ic13)
    n14 = s%net_iso(in14)
    n15 = s%net_iso(in15)
    o16 = s%net_iso(io16)
    o17 = s%net_iso(io17)
    o18 = s%net_iso(io18)
    ne20 = s%net_iso(ine20)

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

    ! global_data(7)-global_data(10) are not used

    if (s%interpolate_rho_for_pulse_data) then
       rho_c = eval_center_rho(s, k_b(n_sg))
    else
       rho_c = eval_center(s%rmid, s%rho, k_a(n_sg), k_b(n_sg))
    endif

    ! at the centre d²P/dr² = -4πGρ²/3
    d2P_dr2_c = -four_thirds*pi*s% cgrav(s% nz)*rho_c**2
    P_c = s%Peos(s% nz) - 0.5*d2P_dr2_c*s% rmid(s% nz)**2
    global_data(11) = r_outer**2*d2P_dr2_c/P_c
    global_data(12) = r_outer**2*eval_center_d2(s%rmid, s%rho, k_a(n_sg), k_b(n_sg)) / rho_c
    global_data(13) = s%star_age
    global_data(14) = s%Teff
    global_data(15) = standard_cgrav

    ! Store point data

    allocate(point_data(IVAR,nn))
    point_data = 0d0
    
    j = 1

    ! Atmosphere (we skip the point at the base of the atm to smooth
    ! the transition)

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

      ! Store data associated with atmosphere point k into the
      ! point_data array at position j

      associate ( &
           r => point_data(1,j), &
           lnq => point_data(2,j), &
           T => point_data(3,j), &
           P => point_data(4,j), &
           rho => point_data(5,j), &
           X => point_data(6,j), &
           L => point_data(7,j), &
           kap => point_data(8,j), &
           eps => point_data(9,j), &
           Gamma_1 => point_data(10,j), &
           nabla_ad => point_data(11,j), &
           delta => point_data(12,j), &
           c_P => point_data(13,j), &
           rec_mu_e => point_data(14,j), &
           A_ast => point_data(15,j), &
           r_X => point_data(16,j), &
           Z => point_data(17,j), &
           R_r => point_data(18,j), &
           eps_grav => point_data(19,j), &
           X_He3 => point_data(21,j), &
           X_C12 => point_data(22,j), &
           X_C13 => point_data(23,j), &
           X_N14 => point_data(24,j), &
           X_O16 => point_data(25,j), &
           X_H2 => point_data(29,j), &
           X_He4 => point_data(30,j), &
           X_Li7 => point_data(31,j), &
           X_Be7 => point_data(32,j), &
           X_N15 => point_data(33,j), &
           X_O17 => point_data(34,j), &
           X_O18 => point_data(35,j), &
           X_Ne20 => point_data(36,j))

        r = s%r(1) + s%atm_structure(atm_delta_r,k)
        lnq = log(s%m_grav(1)/m_outer)
        T = exp(s%atm_structure(atm_lnT,k))
        P = exp(s%atm_structure(atm_lnP,k))
        rho = exp(s%atm_structure(atm_lnd,k))
        X = eval_face(s%dq, s%X, 1, k_a, k_b, v_lo=0d0, v_hi=1d0)
        L = s%L(1)
        kap = s%atm_structure(atm_kap,k)
        eps = 0d0
        Gamma_1 = s%atm_structure(atm_gamma1,k)
        nabla_ad = s%atm_structure(atm_grada,k)
        delta = s%atm_structure(atm_chiT,k)/s%atm_structure(atm_chiRho,k)
        c_P = s%atm_structure(atm_cp,k)
        rec_mu_e = exp(s%atm_structure(atm_lnfree_e,k)) ! check

        grav = s%cgrav(1)*s%m_grav(1)/(r*r)
        N2 = grav*grav*(rho/P)*delta*(nabla_ad - s%atm_structure(atm_gradT,k))
        A_ast = N2*r/grav

        r_X = 0d0 ! dxdt_nuc_h1
        Z = MIN(1d0, MAX(0d0, 1d0 - &
             eval_face(s%dq, s%X, 1, k_a, k_b) - &
             eval_face(s%dq, s%Y, 1, k_a, k_b)))
        R_r = r_outer - r
        eps_grav = 0d0

        X_He3 = eval_face_X(s, he3, 1, k_a, k_b)
        X_C12 = eval_face_X(s, c12, 1, k_a, k_b)
        X_C13 = eval_face_X(s, c13, 1, k_a, k_b)
        X_N14 = eval_face_X(s, n14, 1, k_a, k_b)
        X_O16 = eval_face_X(s, o16, 1, k_a, k_b)
        X_H2 = eval_face_X(s, h2, 1, k_a, k_b)
        X_He4 = eval_face_X(s, he4, 1, k_a, k_b)
        X_Li7 = eval_face_X(s, li7, 1, k_a, k_b)
        X_Be7 = eval_face_X(s, be7, 1, k_a, k_b)
        X_N15 = eval_face_X(s, n15, 1, k_a, k_b)
        X_O17 = eval_face_X(s, o17, 1, k_a, k_b)
        X_O18 = eval_face_X(s, o18, 1, k_a, k_b)
        X_Ne20 = eval_face_X(s, ne20, 1, k_a, k_b)

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

      ! Store data for envelope face k into the point_data array at
      ! position j

      associate ( &
           r => point_data(1,j), &
           lnq => point_data(2,j), &
           T => point_data(3,j), &
           P => point_data(4,j), &
           rho => point_data(5,j), &
           X => point_data(6,j), &
           L => point_data(7,j), &
           kap => point_data(8,j), &
           eps => point_data(9,j), &
           Gamma_1 => point_data(10,j), &
           nabla_ad => point_data(11,j), &
           delta => point_data(12,j), &
           c_P => point_data(13,j), &
           rec_mu_e => point_data(14,j), &
           A_ast => point_data(15,j), &
           r_X => point_data(16,j), &
           Z => point_data(17,j), &
           R_r => point_data(18,j), &
           eps_grav => point_data(19,j), &
           X_He3 => point_data(21,j), &
           X_C12 => point_data(22,j), &
           X_C13 => point_data(23,j), &
           X_N14 => point_data(24,j), &
           X_O16 => point_data(25,j), &
           X_H2 => point_data(29,j), &
           X_He4 => point_data(30,j), &
           X_Li7 => point_data(31,j), &
           X_Be7 => point_data(32,j), &
           X_N15 => point_data(33,j), &
           X_O17 => point_data(34,j), &
           X_O18 => point_data(35,j), &
           X_Ne20 => point_data(36,j))
           
        r = s%r(k)
        lnq = log(s%m_grav(k)/m_outer)
        T = eval_face(s%dq, s%T, k, 1, s%nz)
        P = eval_face(s%dq, s%Peos, k, 1, s%nz)
        if (s%interpolate_rho_for_pulse_data) then
           rho = eval_face(s%dq, s%rho, k, k_a, k_b)
        else
           rho = eval_face_rho(s, k, k_a, k_b)
        end if
        X = eval_face(s%dq, s%X, k, k_a, k_b, v_lo=0d0, v_hi=1d0)
        L = s%L(k)
        kap = eval_face(s%dq, s%opacity, k, k_a, k_b)
        eps = eval_face(s%dq, s%eps_nuc, k, k_a, k_b) + eval_face(s%dq, s%eps_grav_ad%val, k, k_a, k_b)
        Gamma_1 = eval_face(s%dq, s%gamma1, k, k_a, k_b)
        nabla_ad = eval_face(s%dq, s%grada, k, k_a, k_b)
        delta = eval_face(s%dq, s%chiT, k, k_a, k_b)/eval_face(s%dq, s%chiRho, k, k_a, k_b)
        c_P = eval_face(s%dq, s%cp, k, k_a, k_b)
        rec_mu_e = exp(eval_face(s%dq, s%lnfree_e, k, k_a, k_b))
        if (r <= s%fgong_zero_A_inside_r*Rsun .AND. s%mixing_type(k) /= no_mixing) then
           A_ast = 0d0
        else
           A_ast = eval_face_A_ast(s, k, k_a, k_b)
        end if
        r_X = eval_face(s%dq, s%dxdt_nuc(h1,:), k, k_a, k_b)
        Z = MIN(1d0, MAX(0d0, 1d0 - &
             eval_face(s%dq, s%X, k, k_a, k_b) - &
             eval_face(s%dq, s%Y, k, k_a, k_b)))
        R_r = r_outer - s%r(k)
        eps_grav = eval_face(s%dq, s%eps_grav_ad%val, k, k_a, k_b)
        X_He3 = eval_face_X(s, he3, k, k_a, k_b)
        X_C12 = eval_face_X(s, c12, k, k_a, k_b)
        X_C13 = eval_face_X(s, c13, k, k_a, k_b)
        X_N14 = eval_face_X(s, n14, k, k_a, k_b)
        X_O16 = eval_face_X(s, o16, k, k_a, k_b)
        X_H2 = eval_face_X(s, h2, k, k_a, k_b)
        X_He4 = eval_face_X(s, he4, k, k_a, k_b)
        X_Li7 = eval_face_X(s, li7, k, k_a, k_b)
        X_Be7 = eval_face_X(s, be7, k, k_a, k_b)
        X_N15 = eval_face_X(s, n15, k, k_a, k_b)
        X_O17 = eval_face_X(s, o17, k, k_a, k_b)
        X_O18 = eval_face_X(s, o18, k, k_a, k_b)
        X_Ne20 = eval_face_X(s, ne20, k, k_a, k_b)

      end associate

      ! Finish

      return

    end subroutine store_point_data_env

    !****

    subroutine store_point_data_ctr (j, k_a, k_b)

      integer, intent(in) :: j
      integer, intent(in) :: k_a
      integer, intent(in) :: k_b

      real(dp) :: dq_center
      
      ! Store data for the center into the point_data array at position j

      associate ( &
           r => point_data(1,j), &
           lnq => point_data(2,j), &
           T => point_data(3,j), &
           P => point_data(4,j), &
           rho => point_data(5,j), &
           X => point_data(6,j), &
           L => point_data(7,j), &
           kap => point_data(8,j), &
           eps => point_data(9,j), &
           Gamma_1 => point_data(10,j), &
           nabla_ad => point_data(11,j), &
           delta => point_data(12,j), &
           c_P => point_data(13,j), &
           rec_mu_e => point_data(14,j), &
           A_ast => point_data(15,j), &
           r_X => point_data(16,j), &
           Z => point_data(17,j), &
           R_r => point_data(18,j), &
           eps_grav => point_data(19,j), &
           X_He3 => point_data(21,j), &
           X_C12 => point_data(22,j), &
           X_C13 => point_data(23,j), &
           X_N14 => point_data(24,j), &
           X_O16 => point_data(25,j), &
           X_H2 => point_data(29,j), &
           X_He4 => point_data(30,j), &
           X_Li7 => point_data(31,j), &
           X_Be7 => point_data(32,j), &
           X_N15 => point_data(33,j), &
           X_O17 => point_data(34,j), &
           X_O18 => point_data(35,j), &
           X_Ne20 => point_data(36,j))
           
        r = 0d0
        lnq = log(TINY(0d0))
        T = eval_center(s%rmid, s%T, 1, s%nz)
        P = P_c
        rho = rho_c
        X = eval_center(s%rmid, s%X, k_a, k_b, v_lo=0d0, v_hi=1d0)
        L = 0d0
        kap = eval_center(s%rmid, s%opacity, k_a, k_b)
        eps = eval_center(s%rmid, s%eps_nuc, k_a, k_b) + eval_center(s%rmid, s%eps_grav_ad%val, k_a, k_b)
        Gamma_1 = eval_center(s%rmid, s%gamma1, k_a, k_b)
        nabla_ad = eval_center(s%rmid, s%grada, k_a, k_b)
        delta = eval_center(s%rmid, s%chiT, k_a, k_b)/eval_center(s%rmid, s%chiRho, k_a, k_b)
        c_P = eval_center(s%rmid, s%cp, k_a, k_b)
        rec_mu_e = exp(eval_center(s%rmid, s%lnfree_e, k_a, k_b))
        A_ast = 0d0
        r_X = eval_center(s%rmid, s%dxdt_nuc(h1,:), k_a, k_b)
        Z = MIN(1d0, MAX(0d0, 1d0 - &
             eval_center(s%rmid, s%X, k_a, k_b) - &
             eval_center(s%rmid, s%Y, k_a, k_b)))
        R_r = r_outer
        eps_grav = eval_center(s%rmid, s%eps_grav_ad%val, k_a, k_b)
        X_He3 = eval_center_X(s, he3, k_a, k_b)
        X_C12 = eval_center_X(s, c12, k_a, k_b)
        X_C13 = eval_center_X(s, c13, k_a, k_b)
        X_N14 = eval_center_X(s, n14, k_a, k_b)
        X_O16 = eval_center_X(s, o16, k_a, k_b)
        X_H2 = eval_center_X(s, h2, k_a, k_b)
        X_He4 = eval_center_X(s, he4, k_a, k_b)
        X_Li7 = eval_center_X(s, li7, k_a, k_b)
        X_Be7 = eval_center_X(s, be7, k_a, k_b)
        X_N15 = eval_center_X(s, n15, k_a, k_b)
        X_O17 = eval_center_X(s, o17, k_a, k_b)
        X_O18 = eval_center_X(s, o18, k_a, k_b)
        X_Ne20 = eval_center_X(s, ne20, k_a, k_b)

      end associate

      ! Finish

      return

    end subroutine store_point_data_ctr

  end subroutine get_fgong_data

  !****

  subroutine write_fgong_data (id, filename, global_data, point_data, ierr)

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
    character(len=strlen)    :: format_for_fgong_data

    ! Write FGONG data to file

    call get_star_ptr(id, s, ierr)
    if (ierr /= 0) then
       write(*,*) 'bad star id for write_fgong_info'
       return
    end if

    ! Set float format from version number (ivers)

    if (s% fgong_ivers == 300) then
       format_for_fgong_data = '(1P5E16.9,x)'
    else if (s% fgong_ivers == 1300) then
       format_for_fgong_data = '(1P,5(X,E26.18E3))'
    else
       write(*,*) ''
       write(*,'(a,i4)') 'bad fgong_ivers: must be 300 or 1300, not ', s% fgong_ivers
       ierr = 1
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

    do k = 1, 4
       write(iounit, *) trim(s% fgong_header(k))
    end do

    write(iounit,'(4I10)') nn, ICONST, IVAR, s% fgong_ivers

    write(iounit, format_for_fgong_data) (global_data(i), i=1,n_global)

    do k = 1, nn
       write(iounit, format_for_fgong_data) (point_data(j,k), j=1,n_point)
    end do

    ! Close the file

    close(iounit)

    ! Finish

    return

  end subroutine write_fgong_data

end module pulse_fgong
