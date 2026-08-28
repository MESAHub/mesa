! ***********************************************************************
!
!   Copyright (C) 2010-2018  The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! ***********************************************************************

module pulse_gyre

  use star_private_def
  use const_def, only: dp, pi, four_thirds, rsun, no_mixing, sqrt_2_div_3
  use utils_lib
  use atm_def
  use atm_support
  use eps_grav
  use pulse_utils
  use star_utils, only: get_rho_face_val, get_scale_height_face_val

  implicit none

  private
  public :: get_gyre_data
  public :: write_gyre_data

  integer, parameter :: gyre_schema_101_ncols = 18
  integer, parameter :: gyre_schema_120_ncols = 19
  integer, parameter :: gyre_schema_130_ncols = 38

  integer, parameter :: i_tdc_l_conv0 = 20
  integer, parameter :: i_tdc_a0 = 21
  integer, parameter :: i_tdc_d_conv_h = 22
  integer, parameter :: i_tdc_hp_face = 23
  integer, parameter :: i_tdc_alpha_mlt_face = 24
  integer, parameter :: i_tdc_cp_face = 25
  integer, parameter :: i_tdc_chiT_face = 26
  integer, parameter :: i_tdc_chiRho_face = 27
  integer, parameter :: i_tdc_gradL_face = 28
  integer, parameter :: i_tdc_gradT_face = 29
  integer, parameter :: i_tdc_alpha_C = 30
  integer, parameter :: i_tdc_alpha_S = 31
  integer, parameter :: i_tdc_alpha_D = 32
  integer, parameter :: i_tdc_alpha_R = 33
  integer, parameter :: i_tdc_alpha_Pt = 34
  integer, parameter :: i_tdc_alpha_M = 35
  integer, parameter :: i_tdc_include_mlt_corr = 36
  integer, parameter :: i_tdc_mlt_Pturb_factor = 37
  integer, parameter :: i_tdc_mlt_Pturb0 = 38

contains

  subroutine get_gyre_data (id, &
       add_center_point, keep_surface_point, add_atmosphere, global_data, point_data, ierr)

    integer, intent(in)                        :: id
    logical, intent(in)                        :: add_center_point
    logical, intent(in)                        :: keep_surface_point
    logical, intent(in)                        :: add_atmosphere
    real(dp), allocatable, intent(out)         :: global_data(:)
    real(dp), allocatable, target, intent(out) :: point_data(:,:)
    integer, intent(out)                       :: ierr

    type(star_info), pointer :: s
    integer, allocatable     :: k_a(:)
    integer, allocatable     :: k_b(:)
    integer                  :: n_sg
    integer                  :: n_atm
    integer                  :: nn_atm
    integer                  :: n_env
    integer                  :: nn_env
    integer                  :: nn
    real(dp), pointer        :: r(:) => NULL()
    real(dp), pointer        :: m(:) => NULL()
    real(dp), pointer        :: L(:) => NULL()
    real(dp), pointer        :: P(:) => NULL()
    real(dp), pointer        :: T(:) => NULL()
    real(dp), pointer        :: rho(:) => NULL()
    real(dp), pointer        :: nabla(:) => NULL()
    real(dp), pointer        :: N2(:) => NULL()
    real(dp), pointer        :: Gamma_1(:) => NULL()
    real(dp), pointer        :: nabla_ad(:) => NULL()
    real(dp), pointer        :: delta(:) => NULL()
    real(dp), pointer        :: kap(:) => NULL()
    real(dp), pointer        :: kap_kap_T(:) => NULL()
    real(dp), pointer        :: kap_kap_rho(:) => NULL()
    real(dp), pointer        :: eps_nuc(:) => NULL()
    real(dp), pointer        :: eps_eps_T(:) => NULL()
    real(dp), pointer        :: eps_eps_rho(:) => NULL()
    real(dp), pointer        :: eps_grav(:) => NULL()
    real(dp), pointer        :: Omega_rot(:) => NULL()
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

    if (s%gyre_data_schema == 130 .and. .not. s%gyre_write_tdc_lna_background) then
       write(*,*) 'gyre_data_schema = 130 requires gyre_write_tdc_lna_background'
       ierr = -1
       return
    end if
    if (s%gyre_write_tdc_lna_background) then
       if (s%gyre_data_schema /= 130) then
          write(*,*) 'gyre_write_tdc_lna_background requires gyre_data_schema = 130'
          ierr = -1
          return
       end if
       if (trim(s%MLT_option) /= 'TDC') then
          write(*,*) 'gyre_write_tdc_lna_background requires MLT_option = ''TDC'''
          ierr = -1
          return
       end if
       if (s%harmonic_dissipation_length_beta > 0d0) then
          write(*,'(a)') 'WARNING: GYRE schema 130 does not support the harmonic TDC mixing length.'
          write(*,'(a)') 'GYRE will use Lambda = mixing_length_alpha*Hp_face from the exported data.'
       end if
    end if

    ! Set up segment indices

    call set_segment_indices(s, k_a, k_b, add_center_point)

    n_sg = SIZE(k_a)

    ! Determine data dimensions

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
    end if

    if (add_center_point) then
       nn = nn_env + nn_atm + 1
    else
       nn = nn_env + nn_atm
    end if

    ! Allocate arrays & set up data pointers

    allocate(global_data(3))

    ! Set up data pointers

    select case(s%gyre_data_schema)

    case(101,110)

       allocate(point_data(gyre_schema_101_ncols,nn))

       r => point_data(1,:)
       m => point_data(2,:)
       L => point_data(3,:)
       P => point_data(4,:)
       T => point_data(5,:)
       rho => point_data(6,:)
       nabla => point_data(7,:)
       N2 => point_data(8,:)
       Gamma_1 => point_data(9,:)
       nabla_ad => point_data(10,:)
       delta => point_data(11,:)
       kap => point_data(12,:)
       kap_kap_T => point_data(13,:)
       kap_kap_rho => point_data(14,:)
       eps_nuc => point_data(15,:)
       eps_eps_T => point_data(16,:)
       eps_eps_rho => point_data(17,:)
       Omega_rot => point_data(18,:)

    case(120)

       allocate(point_data(gyre_schema_120_ncols,nn))

       r => point_data(1,:)
       m => point_data(2,:)
       L => point_data(3,:)
       P => point_data(4,:)
       T => point_data(5,:)
       rho => point_data(6,:)
       nabla => point_data(7,:)
       N2 => point_data(8,:)
       Gamma_1 => point_data(9,:)
       nabla_ad => point_data(10,:)
       delta => point_data(11,:)
       kap => point_data(12,:)
       kap_kap_T => point_data(13,:)
       kap_kap_rho => point_data(14,:)
       eps_nuc => point_data(15,:)
       eps_eps_T => point_data(16,:)
       eps_eps_rho => point_data(17,:)
       eps_grav => point_data(18,:)
       Omega_rot => point_data(19,:)

    case(130)

       allocate(point_data(gyre_schema_130_ncols,nn))

       r => point_data(1,:)
       m => point_data(2,:)
       L => point_data(3,:)
       P => point_data(4,:)
       T => point_data(5,:)
       rho => point_data(6,:)
       nabla => point_data(7,:)
       N2 => point_data(8,:)
       Gamma_1 => point_data(9,:)
       nabla_ad => point_data(10,:)
       delta => point_data(11,:)
       kap => point_data(12,:)
       kap_kap_T => point_data(13,:)
       kap_kap_rho => point_data(14,:)
       eps_nuc => point_data(15,:)
       eps_eps_T => point_data(16,:)
       eps_eps_rho => point_data(17,:)
       eps_grav => point_data(18,:)
       Omega_rot => point_data(19,:)

    case default

       write(*,*) 'invalid gyre_data_schema'
       ierr = -1
       return

    end select

    ! If necessary, update the eps_grav data in the model

    if (ASSOCIATED(eps_grav)) then

       do k = 1, s%nz
          call eval_eps_grav_and_partials(s, k, ierr)
          if (ierr /= 0) then
             write(*,*) 'failed in call to eval_eps_grav_and_partials'
             return
          end if
       end do

    end if

    ! Store global data

    r_outer = Rsun*s%photosphere_r
    m_outer = s%m_grav(1)

    global_data(1) = m_outer
    global_data(2) = r_outer
    global_data(3) = s%L(1)

    ! Store point data

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

       end if

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

    return

  contains

    subroutine store_point_data_atm (j, k)

      integer, intent(in) :: j
      integer, intent(in) :: k

      real(dp) :: grav

      ! Store data associated with atmosphere point k into the
      ! point_data array at position j

      r(j) = s%r(1) + s%atm_structure(atm_delta_r,k)
      m(j) = s%m_grav(1)  !+ s%atm_structure(atm_delta_m,k)
      L(j) = s%L(1)

      P(j) = exp(s%atm_structure(atm_lnP,k))
      rho(j) = exp(s%atm_structure(atm_lnd,k))
      T(j) = exp(s%atm_structure(atm_lnT,k))

      Gamma_1(j) = s%atm_structure(atm_gamma1,k)
      nabla_ad(j) = s%atm_structure(atm_grada,k)
      delta(j) = s%atm_structure(atm_chiT,k)/s%atm_structure(atm_chiRho,k)
      nabla(j) = s%atm_structure(atm_gradT,k)

      grav = s%cgrav(1)*m(j)/(r(j)*r(j))
      N2(j) = grav*grav*(rho(j)/P(j))*delta(j)*(nabla_ad(j) - nabla(j))

      kap(j) = s%atm_structure(atm_kap,k)
      kap_kap_T(j) = kap(j)*s%atm_structure(atm_dlnkap_dlnT,k)
      kap_kap_rho(j) = kap(j)*s%atm_structure(atm_dlnkap_dlnd,k)

      eps_nuc(j) = 0d0
      eps_eps_T(j) = 0d0
      eps_eps_rho(j) = 0d0
      if (ASSOCIATED(eps_grav)) eps_grav(j) = 0d0

      if (s%rotation_flag) then
         Omega_rot(j) = s%omega(1)
      else
         Omega_rot(j) = 0d0
      end if

      call store_tdc_lna_point_defaults(j)

      return

    end subroutine store_point_data_atm


    subroutine store_point_data_env (j, k, k_a, k_b)

      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: k_a
      integer, intent(in) :: k_b

      ! Store data associated with envelope face k into the point_data
      ! array at position j

      r(j) = s%r(k)
      m(j) = s%m_grav(k)
      L(j) = s%L(k)

      P(j) = eval_face(s%dq, s%Peos, k, 1, s%nz)
      if (s%interpolate_rho_for_pulse_data) then
         rho(j) = eval_face(s%dq, s%rho, k, k_a, k_b)
      else
         rho(j) = eval_face_rho(s, k, k_a, k_b)
      end if
      T(j) = eval_face(s%dq, s%T, k, 1, s%nz)

      N2(j) = eval_face_A_ast(s, k, k_a, k_b)*s%grav(k)/s%r(k)
      Gamma_1(j) = eval_face(s%dq, s%gamma1, k, k_a, k_b)
      nabla_ad(j) = eval_face(s%dq, s%grada, k, k_a, k_b)
      delta(j) = eval_face(s%dq, s%chiT, k, k_a, k_b)/eval_face(s%dq, s%chiRho, k, k_a, k_b)
      nabla(j) = s%gradT(k)  ! Not quite right; gradT can be discontinuous

      kap(j) = eval_face(s%dq, s%opacity, k, k_a, k_b)
      kap_kap_T(j) = eval_face(s%dq, s%d_opacity_dlnT, k, k_a, k_b)
      kap_kap_rho(j) = eval_face(s%dq, s%d_opacity_dlnd, k, k_a, k_b)

      eps_nuc(j) = eval_face(s%dq, s%eps_nuc, k, k_a, k_b)
      eps_eps_T(j) = eval_face(s%dq, s%d_epsnuc_dlnT, k, k_a, k_b)
      eps_eps_rho(j) = eval_face(s%dq, s%d_epsnuc_dlnd, k, k_a, k_b)
      if (ASSOCIATED(eps_grav)) eps_grav(j) = eval_face(s%dq, s%eps_grav_ad%val, k, k_a, k_b)

      if (s%rotation_flag) then
         Omega_rot(j) = s%omega(k)  ! Not quite right; omega can be discontinuous
      else
         Omega_rot = 0d0
      end if

      call store_tdc_lna_point_env(j, k, k_a, k_b)

      return

    end subroutine store_point_data_env


    subroutine store_point_data_ctr (j, k_a, k_b)

      integer, intent(in) :: j
      integer, intent(in) :: k_a
      integer, intent(in) :: k_b

      real(dp) :: d2P_dr2_c

      ! Store data for the center into the point_data array at position j

      r(j) = 0d0
      m(j) = 0d0
      L(j) = 0d0

      P(j) = eval_center(s%rmid, s%Peos, 1, s%nz)
      if (s%interpolate_rho_for_pulse_data) then
         rho(j) = eval_center(s%rmid, s%rho, k_a, k_b)
      else
         rho(j) = eval_center_rho(s, k_b)
      end if
      ! at the centre d²P/dr² = -4πGρ²/3
      d2P_dr2_c = -four_thirds*pi*s% cgrav(s% nz)*rho(j)**2
      P(j) = s%Peos(s% nz) - 0.5d0*d2P_dr2_c*s% rmid(s% nz)**2
      T(j) = eval_center(s%rmid, s%T, 1, s%nz)

      N2(j) = 0d0
      Gamma_1(j) = eval_center(s%rmid, s%gamma1, k_a, k_b)
      nabla_ad(j) = eval_center(s%rmid, s%grada, k_a, k_b)
      delta(j) = eval_center(s%rmid, s%chiT, k_a, k_b)/eval_center(s%rmid, s%chiRho, k_a, k_b)
      nabla(j) = eval_center(s%r, s%gradT, k_a, k_b)

      kap(j) = eval_center(s%rmid, s%opacity, k_a, k_b)
      kap_kap_T(j) = eval_center(s%rmid, s%d_opacity_dlnT, k_a, k_b)
      kap_kap_rho(j) = eval_center(s%rmid, s%d_opacity_dlnd, k_a, k_b)

      eps_nuc(j) = eval_center(s%rmid, s%eps_nuc, k_a, k_b)
      eps_eps_T(j) = eval_center(s%rmid, s%d_epsnuc_dlnT, k_a, k_b)
      eps_eps_rho(j) = eval_center(s%rmid, s%d_epsnuc_dlnd, k_a, k_b)
      if (ASSOCIATED(eps_grav)) eps_grav(j) = eval_center(s%rmid, s%eps_grav_ad%val, k_a, k_b)

      if (s%rotation_flag) then
         Omega_rot(j) = eval_center(s%r, s%omega, k_a, k_b)
      else
         Omega_rot(j) = 0d0
      end if

      call store_tdc_lna_point_defaults(j)

      return

    end subroutine store_point_data_ctr


    subroutine store_tdc_lna_point_defaults(j)

      integer, intent(in) :: j

      if (SIZE(point_data, 1) < gyre_schema_130_ncols) return

      point_data(i_tdc_l_conv0:gyre_schema_130_ncols,j) = 0d0
      if (.not. s%gyre_write_tdc_lna_background) return

      point_data(i_tdc_alpha_C,j) = s%TDC_alpha_C
      point_data(i_tdc_alpha_S,j) = s%TDC_alpha_S
      point_data(i_tdc_alpha_D,j) = s%TDC_alpha_D
      point_data(i_tdc_alpha_R,j) = s%TDC_alpha_R
      point_data(i_tdc_alpha_Pt,j) = s%TDC_alpha_Pt
      point_data(i_tdc_alpha_M,j) = s%TDC_alpha_M
      point_data(i_tdc_include_mlt_corr,j) = merge(1d0, 0d0, s%include_mlt_corr_to_TDC)
      point_data(i_tdc_mlt_Pturb_factor,j) = s%mlt_Pturb_factor

    end subroutine store_tdc_lna_point_defaults


    subroutine store_tdc_lna_point_env(j, k, k_a, k_b)

      integer, intent(in) :: j
      integer, intent(in) :: k
      integer, intent(in) :: k_a
      integer, intent(in) :: k_b

      real(dp) :: Cp_face
      real(dp) :: Hp_face
      real(dp) :: Y_face
      real(dp) :: denom

      call store_tdc_lna_point_defaults(j)

      if (.not. s%gyre_write_tdc_lna_background) return

      Hp_face = get_scale_height_face_val(s, k)
      Cp_face = eval_face(s%dq, s%Cp, k, k_a, k_b)
      Y_face = s%gradT(k) - s%gradL(k)

      ! Schema 130 exports the legacy Lambda = mixing_length_alpha*Hp_face.
      point_data(i_tdc_hp_face,j) = Hp_face
      point_data(i_tdc_alpha_mlt_face,j) = s%mixing_length_alpha
      point_data(i_tdc_cp_face,j) = Cp_face
      point_data(i_tdc_chiT_face,j) = eval_face(s%dq, s%chiT, k, k_a, k_b)
      point_data(i_tdc_chiRho_face,j) = eval_face(s%dq, s%chiRho, k, k_a, k_b)
      point_data(i_tdc_gradL_face,j) = s%gradL(k)
      point_data(i_tdc_gradT_face,j) = s%gradT(k)

      if (.not. tdc_lna_has_convective_velocity(k)) return

      point_data(i_tdc_l_conv0,j) = s%L_conv(k)
      point_data(i_tdc_a0,j) = s%mlt_vc(k)/sqrt_2_div_3
      if (s%mlt_Pturb_factor > 0d0) &
         point_data(i_tdc_mlt_Pturb0,j) = &
            s%mlt_Pturb_factor*get_rho_face_val(s,k)*s%mlt_vc(k)*s%mlt_vc(k)/3d0
      ! Infer D from F_conv = rho*T*Cp*D*(gradT - gradL)/Hp.
      denom = 4d0*pi*r(j)*r(j)*rho(j)*T(j)*Cp_face*Y_face
      if (denom > 0d0 .and. s%L_conv(k) > 0d0) &
         point_data(i_tdc_d_conv_h,j) = s%L_conv(k)*Hp_face/denom

    end subroutine store_tdc_lna_point_env


    logical function tdc_lna_has_convective_velocity(k) result(has_velocity)

      integer, intent(in) :: k

      has_velocity = s%mixing_length_alpha > 0d0 .and. &
         trim(s%MLT_option) == 'TDC' .and. &
         k > s%TDC_num_outermost_cells_forced_nonturbulent .and. &
         k <= s%nz - s%TDC_num_innermost_cells_forced_nonturbulent .and. &
         s%mlt_mixing_type(k) /= no_mixing .and. &
         s%mlt_vc(k) > 0d0

    end function tdc_lna_has_convective_velocity

  end subroutine get_gyre_data


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

    select case(s%gyre_data_schema)
    case(101,120,130)
    case default
       write(*,*) 'invalid gyre_data_schema'
       ierr = -1
       return
    end select

    ! Open the file

    open(newunit=iounit, file=TRIM(filename), status='REPLACE', iostat=ierr)
    if (ierr /= 0) then
       write(*,*) 'failed to open '//TRIM(filename)
       return
    end if

    ! Write the data

    nn = SIZE(point_data, 2)

    write(iounit, 100) nn, global_data, s%gyre_data_schema
100 format(I6, 3(1X,1PE26.16), 1X, I6)

    do j = 1, nn
       write(iounit, 110) j, point_data(:,j)
110    format(I6, 99(1X,1PE26.16))
    end do

    ! Close the file

    close(iounit)

    return

  end subroutine write_gyre_data

end module pulse_gyre
