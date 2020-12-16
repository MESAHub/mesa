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

module pulse_cafein

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

  integer, parameter :: NCOL = 35

  ! Access specifiers

  private

  public :: get_cafein_data
  public :: write_cafein_data

contains

  subroutine get_cafein_data (id, &
       add_center_point, keep_surface_point, add_atmosphere, global_data, point_data, ierr)

    integer, intent(in)                :: id
    logical, intent(in)                :: add_center_point
    logical, intent(in)                :: keep_surface_point
    logical, intent(in)                :: add_atmosphere
    real(dp), allocatable, intent(out) :: global_data(:)
    real(dp), allocatable, intent(out) :: point_data(:,:)
    integer, intent(out)               :: ierr

    type(star_info), pointer :: s
    integer                  :: n_atm
    integer                  :: nn_atm
    integer                  :: n_env
    integer                  :: nn_env
    integer                  :: nn
    real(dp)                 :: R_star
    real(dp)                 :: M_star
    real(dp)                 :: P_c
    real(dp)                 :: T_c
    real(dp)                 :: rho_c
    real(dp)                 :: s_units
    real(dp)                 :: g_units
    real(dp)                 :: cm_units
    real(dp)                 :: erg_units
    real(dp)                 :: K_units
    real(dp), allocatable    :: l_rad(:)
    integer                  :: j
    integer                  :: k

    ! Get model data for CAFein output

    call get_star_ptr(id, s, ierr)
    if (ierr /= 0) then
       write(*,*) 'bad star id for get_cafein_data'
       return
    end if

    ! Determine data dimensiones

    if (add_atmosphere) then
       call build_atm(s, ierr)
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
       nn_env = n_env
    else
       nn_env = n_env - 1
    endif

    if (add_center_point) then
       nn = nn_env + nn_atm + 1
    else
       nn = nn_env + nn_atm
    endif

    ! Store global data

    allocate(global_data(13))

    M_star = s%m_grav(1)
    R_star = Rsun*s%photosphere_r

    P_c = eval_center(s%rmid, s%P, 1, s%nz)
    T_c = eval_center(s%rmid, s%T, 1, s%nz)
    rho_c = eval_center(s%rmid, s%rho, 1, s%nz)

    s_units = SQRT(pow3(R_star)/(s%cgrav(1)*M_star))
    g_units = M_star
    cm_units = R_star
    erg_units = s%cgrav(1)*M_star**2/R_star
    K_units = s%cgrav(1)*M_star**2/(R_star*cgas)

    global_data(1) = M_star/Msun
    global_data(2) = R_star/Rsun
    global_data(3) = s%L(1)
    global_data(4) = P_c
    global_data(5) = T_c
    global_data(6) = rho_c
    global_data(7) = s%Teff
    global_data(8) = s_units
    global_data(9) = g_units
    global_data(10) = cm_units
    global_data(11) = erg_units
    global_data(12) = K_units
    global_data(13) = 0.d0

    ! Store point data. IMPORTANT NOTE: CAFein is ambiguous about
    ! which luminosity (total l or radiative l_rad) should be stored
    ! in output files (both as the Lr field, and for calculation of
    ! the dlnLr/dlnr, c_3 and c_4 fields). From context (i.e., reading
    ! the CAFein code), both dlnL/dlnr and the c_3/c_4 coefficients
    ! should be calculated using l_rad; however, it seems more useful
    ! for the Lr field to contain l. This is what the following code
    ! does.

    allocate(point_data(NCOL,nn))

    allocate(l_rad(nn))

    j = 1

    ! Atmosphere (we skip the point at the base of the atm to smooth
    ! the transition)

    atm_loop : do k = 1, n_atm-1
       call store_point_data_atm(j, k)
       j = j + 1
    end do atm_loop
    
    ! Envelope

    env_loop : do k = 1, n_env

       if (k == 1 .AND. .NOT. keep_surface_point) cycle env_loop

       call store_point_data_env(j, k)
       j = j + 1

    end do env_loop

    ! Center

    if (add_center_point) then
       call store_point_data_ctr(j)
       j = j + 1
    end if

    ! Check that all point data has correctly been calculated

    if (j /= nn+1) stop 'Invalid cell index in get_cafein_data'

    ! Reverse the data ordering (CAFein format is center-to-surface)

    point_data = point_data(:,nn:1:-1)

    l_rad = l_rad(nn:1:-1)

   ! Set remaining point data

    associate ( &
         r => point_data(5,:), &
         U => point_data(8,:), &
         N2 => point_data(10,:), &
         V => point_data(14,:), &
         nabla_ad => point_data(15,:), &
         c_2_fit => point_data(22,:), &
         dlnLr_dlnr => point_data(25,:), & 
         surf_r_rad_beg => point_data(27,:), &
         surf_r_rad_end => point_data(28,:), &
         U_U_surf => point_data(34,:), &
         c_2 => point_data(35,:))

      c_2 = c_2 + nabla_ad*(log_deriv(r, nabla_ad, dy_a=0d0) + V)
      c_2_fit = c_2

      dlnLr_dlnr = log_deriv(r, l_rad, dy_a=3d0)

      U_U_surf = U/U(nn)

      surf_r_rad_beg = R_star

      do j = 1, nn
         surf_r_rad_end = r(j)
         if (N2(j) < 0d0) exit
      end do

    end associate

    ! Convert dimensioned data to dimensionless

    associate ( &
         m => point_data(1,:), &
         logrho => point_data(2,:) , &
         logP => point_data(3,:), &
         logT => point_data(4,:), &
         r => point_data(5,:), &
         N2 => point_data(10,:), &
         L2_ll1 => point_data(11,:), &
         g => point_data(12,:), &
         l => point_data(13,:), &
         C_P => point_data(17,:), &
         P_scale => point_data(26,:), &
         surf_r_rad_beg => point_data(27,:), &
         surf_r_rad_end => point_data(28,:), &
         entropy => point_data(30,:), &
         kap => point_data(31,:))

      m = m/g_units

      r = r/cm_units
      l = l/(erg_units/s_units)

      logrho = logrho - log10(g_units/cm_units**3)
      logP = logP - log10(erg_units/cm_units**3)
      logT = logT - log10(K_units)

      N2 = N2*s_units**2
      L2_ll1 = L2_ll1*s_units**2

      g = g/(cm_units/s_units**2)

      C_p = C_p/(erg_units/(K_units*g_units))

      P_scale = P_scale/cm_units

      surf_r_rad_beg = surf_r_rad_beg/cm_units
      surf_r_rad_end = surf_r_rad_end/cm_units

      entropy = entropy/(erg_units/(K_units*g_units))

      kap = kap/(cm_units**2/g_units)

    end associate

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

      real(dp) :: rho
      real(dp) :: P
      real(dp) :: T
      real(dp) :: g
      real(dp) :: kap_rho
      real(dp) :: kap_T
      real(dp) :: kap_ad
      real(dp) :: eps

      ! Store data associated with atmosphere point k into the
      ! point_data array at position j

      associate ( &
           m => point_data(1,j), &
           logrho => point_data(2,j) , &
           logP => point_data(3,j), &
           logT => point_data(4,j), &
           r => point_data(5,j), &
           V_g => point_data(6,j), &
           As => point_data(7,j), &
           U => point_data(8,j), &
           c_1 => point_data(9,j), &
           N2 => point_data(10,j), &
           L2_ll1 => point_data(11,j), &
           g => point_data(12,j), &
           l => point_data(13,j), &
           V => point_data(14,j), &
           nabla_ad => point_data(15,j), &
           nabla => point_data(16,j), &
           C_P => point_data(17,j), &
           delta => point_data(18,j), &
           kap_S => point_data(19,j), &
           eps_ad => point_data(20,j), &
           eps_S => point_data(21,j), &
           c_3 => point_data(23,j), &
           c_4 => point_data(24,j), &
           P_scale => point_data(26,j), &
           Gamma_1 => point_data(29,j), &
           entropy => point_data(30,j), &
           kap => point_data(31,j), &
           chi_rho => point_data(32,j), &
           chi_T => point_data(33,j), &
           c_2 => point_data(35,j))

        m = s%m_grav(1) !+ s%atm_structure(atm_delta_m,k)
        r = s%r(1) + s%atm_structure(atm_delta_r,k)
        l = s%L(1)

        logrho = s%atm_structure(atm_lnd,k)/log(10d0)
        logP = s%atm_structure(atm_lnP,k)/log(10d0)
        logT = s%atm_structure(atm_lnT,k)/log(10d0)

        rho = pow(10d0, logrho)
        P = pow(10d0, logP)
        T = pow(10d0, logT)

        nabla_ad = s%atm_structure(atm_grada,k)
        nabla = s%atm_structure(atm_gradT,k)
        C_P = s%atm_structure(atm_cp,k)
        delta = s%atm_structure(atm_chiT,k)/s%atm_structure(atm_chiRho,k)
        Gamma_1 = s%atm_structure(atm_gamma1,k)
        entropy = 0d0 ! No entropy data in surface layers
        chi_rho = s%atm_structure(atm_chiRho,k)
        chi_T = s%atm_structure(atm_chiT,k)

        kap = s%atm_structure(atm_kap,k)
        kap_rho = s%atm_structure(atm_dlnkap_dlnd,k)
        kap_T = s%atm_structure(atm_dlnkap_dlnT,k)
        kap_ad = nabla_ad*kap_T + kap_rho/Gamma_1
        kap_S = kap_T - delta*kap_rho
        
        eps = 0d0
        eps_ad = 0d0
        eps_S = 0d0

        g = s%cgrav(1)*m/(r*r)
        N2 = g*g*(rho/P)*delta*(nabla_ad - nabla)
        L2_ll1 = Gamma_1*P/(rho*r**2)
        P_scale = P/(rho*g)

        V = rho*g*r/P
        V_g = V/Gamma_1
        As = N2*r/g
        U = 4d0*pi*rho*r**3/m

        l_rad(j) = l

        c_1 = (r/R_star)**3*(M_star/m)
        c_2 = (kap_ad - 4d0*nabla_ad)*V*nabla ! Note -- we omit the nabla_ad*(dnabla_ad + V) term for now
        c_3 = 0d0
        c_4 = 4d0*pi*r**3*rho*T*c_P/l_rad(j)*SQRT(s%cgrav(1)*M_star/R_star**3)

      end associate

      ! Finish

      return

    end subroutine store_point_data_atm

    !****

    subroutine store_point_data_env (j, k)

      integer, intent(in) :: j
      integer, intent(in) :: k

      real(dp) :: rho
      real(dp) :: P
      real(dp) :: T
      real(dp) :: kap_rho
      real(dp) :: kap_T
      real(dp) :: kap_ad
      real(dp) :: eps
      real(dp) :: eps_rho
      real(dp) :: eps_T

      ! Store data associated with envelope face k into the point_data
      ! array at position j

      ! Store data associated with atmosphere point k into the
      ! point_data array at position j

      associate ( &
           m => point_data(1,j), &
           logrho => point_data(2,j) , &
           logP => point_data(3,j), &
           logT => point_data(4,j), &
           r => point_data(5,j), &
           V_g => point_data(6,j), &
           As => point_data(7,j), &
           U => point_data(8,j), &
           c_1 => point_data(9,j), &
           N2 => point_data(10,j), &
           L2_ll1 => point_data(11,j), &
           g => point_data(12,j), &
           l => point_data(13,j), &
           V => point_data(14,j), &
           nabla_ad => point_data(15,j), &
           nabla => point_data(16,j), &
           C_P => point_data(17,j), &
           delta => point_data(18,j), &
           kap_S => point_data(19,j), &
           eps_ad => point_data(20,j), &
           eps_S => point_data(21,j), &
           c_3 => point_data(23,j), &
           c_4 => point_data(24,j), &
           P_scale => point_data(26,j), &
           Gamma_1 => point_data(29,j), &
           entropy => point_data(30,j), &
           kap => point_data(31,j), &
           chi_rho => point_data(32,j), &
           chi_T => point_data(33,j), &
           c_2 => point_data(35,j))

        m = s%m_grav(k)
        r = s%r(k)
        l = s%L(k)
        
        if (s%interpolate_rho_for_pulse_data) then
           rho = eval_face(s%dq, s%rho, k, 1, s%nz)
        else
           rho = eval_face_rho(s, k, 1, s%nz)
        end if
        P = eval_face(s%dq, s%P, k, 1, s%nz)
        T = eval_face(s%dq, s%T, k, 1, s%nz)

        logrho = log10(rho)
        logP = log10(P)
        logT = log10(T)

        nabla_ad = eval_face(s%dq, s%grada, k, 1, s%nz)
        nabla = s%gradT(k)
        C_P = eval_face(s%dq, s%cp, k, 1, s%nz)
        delta = eval_face(s%dq, s%chiT, k, 1, s%nz)/eval_face(s%dq, s%chiRho, k, 1, s%nz)
        Gamma_1 = eval_face(s%dq, s%gamma1, k, 1, s%nz)
        entropy = exp(eval_face(s%dq, s%lnS, k, 1, s%nz))
        chi_rho = eval_face(s%dq, s%chiRho, k, 1, s%nz)
        chi_T = eval_face(s%dq, s%chiT, k, 1, s%nz)

        kap = eval_face(s%dq, s%opacity, k, 1, s%nz)
        kap_rho = eval_face(s%dq, s%d_opacity_dlnd, k, 1, s%nz)/kap
        kap_T = eval_face(s%dq, s%d_opacity_dlnT, k, 1, s%nz)/kap
        kap_ad = nabla_ad*kap_T + kap_rho/Gamma_1
        kap_S = kap_T - delta*kap_rho
        
        eps = eval_face(s%dq, s%eps_nuc, k, 1, s%nz)
         if (ABS(eps) > 1D-99) then
           eps_rho = eval_face(s%dq, s%d_epsnuc_dlnd, k, 1, s%nz)/eps
           eps_T = eval_face(s%dq, s%d_epsnuc_dlnT, k, 1, s%nz)/eps
        else
           eps_rho = 0d0
           eps_T = 0d0
        endif
        eps_ad = nabla_ad*eps_T + eps_rho/Gamma_1
        eps_S = eps_T - delta*eps_rho

        g = s%grav(k)
        N2 = eval_face_A_ast(s, k, 1, s%nz)*g/r
        L2_ll1 = Gamma_1*P/(rho*r**2)
        P_scale = s%scale_height(k)
        
        V = rho*g*r/P
        V_g = V/Gamma_1
        As = N2*r/g
        U = 4d0*pi*rho*r**3/m

        l_rad(j) = 16d0*pi*r*crad*clight*T**4*nabla*V/(3d0*kap*rho)

        c_1 = (r/R_star)**3*(M_star/m)
        c_2 = (kap_ad - 4d0*nabla_ad)*V*nabla ! Note -- we omit the nabla_ad*(dnabla_ad + V) term for now
        c_3 = 4d0*pi*r**3*rho*eps/l_rad(j)
        c_4 = 4d0*pi*r**3*rho*T*c_P/l_rad(j)*SQRT(s%cgrav(1)*M_star/R_star**3)
        
      end associate

      ! Finish

      return

    end subroutine store_point_data_env

    !****

    subroutine store_point_data_ctr (j)

      integer, intent(in) :: j

      real(dp) :: rho
      real(dp) :: P
      real(dp) :: T
      real(dp) :: kap_rho
      real(dp) :: kap_T
      real(dp) :: kap_ad
      real(dp) :: eps
      real(dp) :: eps_rho
      real(dp) :: eps_T
      real(dp) :: cgrav

      ! Store data for the center into the point_data array at position j

      associate ( &
           m => point_data(1,j), &
           logrho => point_data(2,j) , &
           logP => point_data(3,j), &
           logT => point_data(4,j), &
           r => point_data(5,j), &
           V_g => point_data(6,j), &
           As => point_data(7,j), &
           U => point_data(8,j), &
           c_1 => point_data(9,j), &
           N2 => point_data(10,j), &
           L2_ll1 => point_data(11,j), &
           g => point_data(12,j), &
           l => point_data(13,j), &
           V => point_data(14,j), &
           nabla_ad => point_data(15,j), &
           nabla => point_data(16,j), &
           C_P => point_data(17,j), &
           delta => point_data(18,j), &
           kap_S => point_data(19,j), &
           eps_ad => point_data(20,j), &
           eps_S => point_data(21,j), &
           c_3 => point_data(23,j), &
           c_4 => point_data(24,j), &
           P_scale => point_data(26,j), &
           Gamma_1 => point_data(29,j), &
           entropy => point_data(30,j), &
           kap => point_data(31,j), &
           chi_rho => point_data(32,j), &
           chi_T => point_data(33,j), &
           c_2 => point_data(35,j))

        m = 0d0
        r = 0d0
        l = 0d0

        if (s%interpolate_rho_for_pulse_data) then
           rho = eval_center(s%rmid, s%rho, 1, s%nz)
        else
           rho = eval_center_rho(s, s%nz)
        end if
        P = eval_center(s%rmid, s%P, 1, s%nz)
        T = eval_center(s%rmid, s%T, 1, s%nz)

        logrho = log10(rho)
        logP = log10(P)
        logT = log10(T)

        nabla_ad = eval_center(s%rmid, s%grada, 1, s%nz)
        nabla = eval_center(s%r, s%gradT, 1, s%nz)
        C_P = eval_center(s%rmid, s%cp, 1, s%nz)
        delta = eval_center(s%rmid, s%chiT, 1, s%nz)/eval_center(s%rmid, s%chiRho, 1, s%nz)
        Gamma_1 = eval_center(s%rmid, s%gamma1, 1, s%nz)
        entropy = exp(eval_center(s%rmid, s%lnS, 1, s%nz))
        chi_rho = eval_center(s%rmid, s%chiRho, 1, s%nz)
        chi_T = eval_center(s%rmid, s%chiT, 1, s%nz)
        
        kap = eval_center(s%rmid, s%opacity, 1, s%nz)
        kap_rho = eval_center(s%rmid, s%d_opacity_dlnd, 1, s%nz)/kap
        kap_T = eval_center(s%rmid, s%d_opacity_dlnT, 1, s%nz)/kap
        kap_ad = nabla_ad*kap_T + kap_rho/Gamma_1
        kap_S = kap_T - delta*kap_rho

        eps = eval_center(s%rmid, s%eps_nuc, 1, s%nz)
        if (ABS(eps) > 1D-99) then
           eps_rho = eval_center(s%rmid, s%d_epsnuc_dlnd, 1, s%nz)/eps
           eps_T = eval_center(s%rmid, s%d_epsnuc_dlnT, 1, s%nz)/eps
        else
           eps_rho = 0d0
           eps_T = 0d0
        endif
        eps_ad = nabla_ad*eps_T + eps_rho/Gamma_1
        eps_S = eps_T - delta*eps_rho
 
        g = 0d0
        N2 = 0d0
        L2_ll1 = HUGE(0d0)
        P_scale = eval_center(s%r, s%scale_height, 1, s%nz)
        
        V = 0d0
        V_g = 0d0
        As = 0d0
        U = 3d0

        l_rad(j) = 0d0

        cgrav = eval_center(s%r, s%cgrav, 1, s%nz)

        c_1 = 3d0*(M_star/R_star**3)/(4d0*pi*rho)
        c_2 = 0d0
        c_3 = 9d0*eps*kap*P/(16*pi*crad*clight*T**4*nabla*cgrav)
        c_4 = 9d0*T*c_P*kap*P/(16*pi*crad*clight*T**4*nabla*cgrav)*SQRT(s%cgrav(s%nz)*M_star/R_star**3)

      end associate

      ! Finish

      return

    end subroutine store_point_data_ctr

    !****

    function log_deriv (x, y, dy_a, dy_b) result (dy)

      real(dp), intent(in)           :: x(:)
      real(dp), intent(in)           :: y(:)
      real(dp), intent(in), optional :: dy_a
      real(dp), intent(in), optional :: dy_b
      real(dp)                       :: dy(SIZE(x))

      integer :: n
      integer :: j

      ! Evaluate the logarithmic derivative of y wrt x, using simple
      ! finite differences

      n = SIZE(x)

      if (PRESENT(dy_a)) then
         dy(1) = dy_a
      else
         dy(1) = x(1)/y(1) * (y(2) - y(1))/(x(2) - x(1))
      endif

      do j = 2, n-1
         dy(j) = x(j)/y(j) * (y(j+1) - y(j-1))/(x(j+1) - x(j-1))
      end do

      if (PRESENT(dy_b)) then
         dy(n) = dy_b
      else
         dy(n) = x(n)/y(n) * (y(n) - y(n-1))/(x(n) - x(n-1))
      endif

      ! Finish

      return

    end function log_deriv
    
  end subroutine get_cafein_data

  !****

  subroutine write_cafein_data (id, filename, global_data, point_data, ierr)

    integer, intent(in)      :: id
    character(*), intent(in) :: filename
    real(dp), intent(in)     :: global_data(:)
    real(dp), intent(in)     :: point_data(:,:)
    integer, intent(out)     :: ierr

    integer :: iounit
    integer :: nn
    integer :: j

    ! Write CAFein data to file

    ! Open the file

    open(newunit=iounit, file=TRIM(filename), status='REPLACE', iostat=ierr)
    if (ierr /= 0) then
       write(*,*) 'failed to open '//TRIM(filename)
       return
    end if

    ! Write the data

    nn = SIZE(point_data, 2)

    write(iounit, 100) &
         '..................M(Msun)' // '..................R(Rsun)' // '..................L(Lsun)' // &
         '..............Pc(erg/cm3)' // '....................Tc(K)' // '..............rhoC(g/cm3)' // &
         '..................Teff(K)' // '................sec_units' // '..............grams_units' // &
         '.................cm_units' // '................erg_units' // '..................K_units' // &
         '...................unused'
100 format(A)

    write(iounit, 110) global_data
110 format(13(1X,E24.14E3))

    write(iounit, *)
    write(iounit, 100) &
         '........................m' // '...................logrho' // '.....................logP' // &
         '.....................logT' // '........................r' // '.......................Vg' // &
         '....................Astar' // '........................U' // '.......................c1' // &
         '.......................N2' // '..............L2/(l(l+1))' // '........................g' // &
         '.......................Lr' // '........................V' // '............del_ad(noDim)' // &
         '...............del(noDim)' // '.......................Cp' // '......................v_t' // &
         '...........kappa_s(noDim)' // '............eps_ad(noDim)' // '.............eps_s(noDim)' // &
         '................c2_fitted' // '.......................c3' // '.......................c4' // &
         '...............dlnLr/dlnr' // '..................P_scale' // '........StartRadSurfLayer' // &
         '..........EndRadSurfLayer' // '...................Gamma1' // '..................entropy' // &
         '....................kappa' // '...................chiRho' // '.....................chiT' // &
         '..................U/Usurf' // '.......................c2'

    do j = 1, nn
       write(iounit, 120) point_data(:,j)
120    format(35(1X,E24.14E3))
    end do

    ! Finish

    ! Close the file
    
    close(iounit)

    ! Finish

    return

  end subroutine write_cafein_data

end module pulse_cafein
