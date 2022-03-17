! ***********************************************************************
!
!   Copyright (C) 2010-2018  The MESA Team
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

module pulse_gyre

  ! Uses

  use star_private_def
  use const_def
  use utils_lib
  use atm_def
  use atm_support

  use pulse_utils

  ! No implicit typing

  implicit none

  ! Parameters

  integer, parameter :: NCOL = 18

  ! Access specifiers

  private

  public :: get_gyre_data
  public :: write_gyre_data

contains

  subroutine get_gyre_data (id, &
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
    integer                  :: n_atm
    integer                  :: nn_atm
    integer                  :: n_env
    integer                  :: nn_env
    integer                  :: nn
    real(dp)                 :: r_outer
    real(dp)                 :: m_outer
    integer                  :: j
    integer                  :: k
    integer                  :: sg

    ! Get model data for GYRE output

    call get_star_ptr(id, s, ierr)
    if (ierr /= 0) then
       write(*,*) 'bad star id for get_gyre_data'
       return
    end if

    ! Set up segment indices

    call set_segment_indices(s, k_a, k_b, add_center_point)

    n_sg = SIZE(k_a)

    ! Determine data dimensiones

    if (add_atmosphere) then
       call build_atm(s, s%L(1), s%r(1), s%Teff, s%m_grav(1), s%cgrav(1), ierr)
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

    allocate(global_data(3))

    r_outer = Rsun*s%photosphere_r
    m_outer = s%m_grav(1)

    global_data(1) = m_outer
    global_data(2) = r_outer
    global_data(3) = s%L(1)

    ! Store point data

    allocate(point_data(NCOL,nn))

    j = 1

    ! Atmosphere (we skip the point at the base of the atm to smooth
    ! the transition)

    atm_loop : do k = 1, n_atm-1
       call store_point_data_atm(j, k)
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

    ! Check that all point data has correctly been calculated

    if (j /= nn+1) call mesa_error(__FILE__,__LINE__,'Invalid cell index in get_gyre_data')

    ! Reverse the data ordering (GYRE format is center-to-surface)

    point_data = point_data(:,nn:1:-1)

    ! Tidy up

    if (ASSOCIATED(s%atm_structure)) then
       deallocate(s%atm_structure)
    end if

    ! Finish

    return

  contains

    subroutine store_point_data_atm (j, k)

      integer, intent(in) :: j
      integer, intent(in) :: k

      real(dp) :: grav

      ! Store data associated with atmosphere point k into the
      ! point_data array at position j

      associate ( &
           r => point_data(1,j), &
           m => point_data(2,j), &
           L => point_data(3,j), &
           P => point_data(4,j), &
           T => point_data(5,j), &
           rho => point_data(6,j), &
           nabla => point_data(7,j), &
           N2 => point_data(8,j), &
           Gamma_1 => point_data(9,j), &
           nabla_ad => point_data(10,j), &
           delta => point_data(11,j), &
           kap => point_data(12,j), &
           kap_T => point_data(13,j), &
           kap_rho => point_data(14,j), &
           eps => point_data(15,j), &
           eps_T => point_data(16,j), &
           eps_rho => point_data(17,j), &
           omega => point_data(18,j))

        r = s%r(1) + s%atm_structure(atm_delta_r,k)
        m = s%m_grav(1) !+ s%atm_structure(atm_delta_m,k)
        L = s%L(1)
        P = exp(s%atm_structure(atm_lnP,k))
        rho = exp(s%atm_structure(atm_lnd,k))
        T = exp(s%atm_structure(atm_lnT,k))
        Gamma_1 = s%atm_structure(atm_gamma1,k)
        nabla_ad = s%atm_structure(atm_grada,k)
        delta = s%atm_structure(atm_chiT,k)/s%atm_structure(atm_chiRho,k)
        nabla = s%atm_structure(atm_gradT,k)

        grav = s%cgrav(1)*m/(r*r)
        N2 = grav*grav*(rho/P)*delta*(nabla_ad - nabla)

        kap = s%atm_structure(atm_kap,k)
        kap_rho = kap*s%atm_structure(atm_dlnkap_dlnd,k)*kap
        kap_T = kap*s%atm_structure(atm_dlnkap_dlnT,k)*kap
        eps = 0d0
        eps_rho = 0d0
        eps_T = 0d0
        if (s%rotation_flag) then
           omega = s%omega(1)
        else
           omega = 0d0
        end if

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

      ! Store data associated with envelope face k into the point_data
      ! array at position j

      associate ( &
           r => point_data(1,j), &
           m => point_data(2,j), &
           L => point_data(3,j), &
           P => point_data(4,j), &
           T => point_data(5,j), &
           rho => point_data(6,j), &
           nabla => point_data(7,j), &
           N2 => point_data(8,j), &
           Gamma_1 => point_data(9,j), &
           nabla_ad => point_data(10,j), &
           delta => point_data(11,j), &
           kap => point_data(12,j), &
           kap_T => point_data(13,j), &
           kap_rho => point_data(14,j), &
           eps => point_data(15,j), &
           eps_T => point_data(16,j), &
           eps_rho => point_data(17,j), &
           omega => point_data(18,j))

        r = s%r(k)
        m = s%m_grav(k)
        L = s%L(k)
        P = eval_face(s%dq, s%Peos, k, 1, s%nz)
        if (s% ctrl% interpolate_rho_for_pulse_data) then
           rho = eval_face(s%dq, s%rho, k, k_a, k_b)
        else
           rho = eval_face_rho(s, k, k_a, k_b)
        end if
        T = eval_face(s%dq, s%T, k, 1, s%nz)
        N2 = eval_face_A_ast(s, k, k_a, k_b)*s%grav(k)/s%r(k)
        Gamma_1 = eval_face(s%dq, s%gamma1, k, k_a, k_b)
        nabla_ad = eval_face(s%dq, s%grada, k, k_a, k_b)
        delta = eval_face(s%dq, s%chiT, k, k_a, k_b)/eval_face(s%dq, s%chiRho, k, k_a, k_b)
        nabla = s%gradT(k) ! Not quite right; gradT can be discontinuous
        kap = eval_face(s%dq, s%opacity, k, k_a, k_b)
        kap_rho = eval_face(s%dq, s%d_opacity_dlnd, k, k_a, k_b)
        kap_T = eval_face(s%dq, s%d_opacity_dlnT, k, k_a, k_b)
        eps = eval_face(s%dq, s%eps_nuc, k, k_a, k_b)
        eps_rho = eval_face(s%dq, s%d_epsnuc_dlnd, k, k_a, k_b)
        eps_T = eval_face(s%dq, s%d_epsnuc_dlnT, k, k_a, k_b)
        if (s%rotation_flag) then
           omega = s%omega(k) ! Not quite right; omega can be discontinuous
        else
           omega = 0d0
        end if

      end associate

      ! Finish

      return

    end subroutine store_point_data_env

    !****

    subroutine store_point_data_ctr (j, k_a, k_b)

      integer, intent(in) :: j
      integer, intent(in) :: k_a
      integer, intent(in) :: k_b

      real(dp) :: d2P_dr2_c

      ! Store data for the center into the point_data array at position j

      associate ( &
           r => point_data(1,j), &
           m => point_data(2,j), &
           L => point_data(3,j), &
           P => point_data(4,j), &
           T => point_data(5,j), &
           rho => point_data(6,j), &
           nabla => point_data(7,j), &
           N2 => point_data(8,j), &
           Gamma_1 => point_data(9,j), &
           nabla_ad => point_data(10,j), &
           delta => point_data(11,j), &
           kap => point_data(12,j), &
           kap_T => point_data(13,j), &
           kap_rho => point_data(14,j), &
           eps => point_data(15,j), &
           eps_T => point_data(16,j), &
           eps_rho => point_data(17,j), &
           omega => point_data(18,j))

        r = 0d0
        m = 0d0
        L = 0d0
        P = eval_center(s%rmid, s%Peos, 1, s%nz)
        if (s% ctrl% interpolate_rho_for_pulse_data) then
           rho = eval_center(s%rmid, s%rho, k_a, k_b)
        else
           rho = eval_center_rho(s, k_b)
        end if

        ! at the centre d²P/dr² = -4πGρ²/3
        d2P_dr2_c = -four_thirds*pi*s% cgrav(s% nz)*rho**2
        P = s%Peos(s% nz) - 0.5*d2P_dr2_c*s% rmid(s% nz)**2
        T = eval_center(s%rmid, s%T, 1, s%nz)
        N2 = 0d0
        Gamma_1 = eval_center(s%rmid, s%gamma1, k_a, k_b)
        nabla_ad = eval_center(s%rmid, s%grada, k_a, k_b)
        delta = eval_center(s%rmid, s%chiT, k_a, k_b)/eval_center(s%rmid, s%chiRho, k_a, k_b)
        nabla = eval_center(s%r, s%gradT, k_a, k_b)
        kap = eval_center(s%rmid, s%opacity, k_a, k_b)
        kap_rho = eval_center(s%rmid, s%d_opacity_dlnd, k_a, k_b)
        kap_T = eval_center(s%rmid, s%d_opacity_dlnT, k_a, k_b)
        eps = eval_center(s%rmid, s%eps_nuc, k_a, k_b)
        eps_rho = eval_center(s%rmid, s%d_epsnuc_dlnd, k_a, k_b)
        eps_T = eval_center(s%rmid, s%d_epsnuc_dlnT, k_a, k_b)
        if (s%rotation_flag) then
           omega = eval_center(s%r, s%omega, k_a, k_b)
        else
           omega = 0d0
        end if

      end associate

      ! Finish

      return

    end subroutine store_point_data_ctr
    
  end subroutine get_gyre_data

  !****

  subroutine write_gyre_data (id, filename, global_data, point_data, ierr)

    integer, intent(in)      :: id
    character(*), intent(in) :: filename
    real(dp), intent(in)     :: global_data(:)
    real(dp), intent(in)     :: point_data(:,:)
    integer, intent(out)     :: ierr

    type(star_info), pointer :: s
    integer                  :: iounit
    integer                  :: nn
    integer                  :: j

    ! Write GYRE data to file

    call get_star_ptr(id, s, ierr)
    if (ierr /= 0) then
       write(*,*) 'bad star id for write_gyre_data'
       return
    end if

    ! Open the file

    open(newunit=iounit, file=TRIM(filename), status='REPLACE', iostat=ierr)
    if (ierr /= 0) then
       write(*,*) 'failed to open '//TRIM(filename)
       return
    end if

    ! Write the data

    nn = SIZE(point_data, 2)

    write(iounit, 100) nn, global_data, GYRE_MODEL_VERSION
100 format(I6, 3(1X,1PE26.16), 1X, I6)

    do j = 1, nn
       write(iounit, 110) j, point_data(:,j)
110    format(I6, 99(1X,1PE26.16))
    end do

    ! Close the file
    
    close(iounit)

    ! Finish

    return

  end subroutine write_gyre_data

end module pulse_gyre
