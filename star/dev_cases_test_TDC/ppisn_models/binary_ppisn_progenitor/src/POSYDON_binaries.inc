! Contains the custom MESA routines for POSYDON integration
! author: Mathieu Renzo - copy paste of POSYDON repo and adapted to MESA r15140

subroutine POSYDON_tsync(id, sync_type, Ftid, qratio, m, r_phot, osep, t_sync, ierr)
  integer, intent(in) :: id
  character (len=strlen), intent(in) :: sync_type !synchronization timescale
  real(dp), intent(in) :: Ftid ! efficiency of tidal synchronization. (time scale / Ftid ).
  real(dp), intent(in) :: qratio !mass_other_star/mass_this_star
  real(dp), intent(in) :: m
  real(dp), intent(in) :: r_phot
  real(dp), intent(in) :: osep ! orbital separation (cm)
  real(dp), intent(out) :: t_sync
  integer, intent(out) :: ierr
  real(dp) :: rGyr_squared , moment_of_inertia
  real(dp) :: one_div_t_sync_conv, one_div_t_sync_rad, one_div_t_sync ! t_sync_rad, t_sync_conv
  type (binary_info), pointer :: b
  type (star_info), pointer :: s

  ierr = 0
  call star_ptr(id, s, ierr)
  if (ierr /= 0) then
     write(*,*) 'failed in star_ptr'
     return
  end if

  call binary_ptr(s% binary_id, b, ierr)
  if (ierr /= 0) then
     write(*,*) 'failed in binary_ptr'
     return
  end if
  moment_of_inertia = dot_product(s% i_rot(:s% nz), s% dm_bar(:s%nz))
  rGyr_squared = (moment_of_inertia/(m*r_phot*r_phot))

  ! Implemented the option for both equilibrium and dynamical tides
  if (sync_type == "Hut_conv") then
     !sync_type .eq. "Hut_conv"!Convective envelope + Radiative core
     ! eq. (11) of Hut, P. 1981, A&A, 99, 126
     t_sync = 3.0d0*k_div_T(b, s,.true.)*(qratio*qratio/rGyr_squared)*pow6(r_phot/osep)
     ! invert it.
     !write(*,*) 'star id', s% id
     if (b% point_mass_i /= 1 .and. b% s1% id == s% id) then
        b% s1% xtra(2) = 1d0/t_sync
        !write(*,*) 'two timescales ', b% s1% xtra1, b% s1% xtra2
     else if (b% point_mass_i /= 2 .and. b% s2% id == s% id) then
        b% s2% xtra(2) = 1d0/t_sync
        !write(*,*) 'two timescales ', b% s2% xtra1, b% s2% xtra2
     else
        write(*,*) 'something is not going well with the stars IDs '
     end if
     t_sync = 1d0/t_sync
     !write(*,*) 'Hut_conv ', t_sync
  else if (sync_type == "Hut_rad") then
     !sync_type .eq. "Hut_rad"! Radiative envelope + convective core
     ! eq. (11) of Hut, P. 1981, A&A, 99, 126
     t_sync = 3.0*k_div_T(b, s,.false.)*(qratio*qratio/rGyr_squared)*pow6(r_phot/osep)
     ! invert it.
     !write(*,*) 'star id', s% id
     if (b% point_mass_i /= 1 .and. b% s1% id == s% id) then
        b% s1% xtra(3) = 1d0/t_sync
        !write(*,*) 'two timescales ', b% s1% xtra1, b% s1% xtra2
     else if (b% point_mass_i /= 2 .and. b% s2% id == s% id) then
        b% s2% xtra(3) = 1d0/t_sync
        !write(*,*) 'two timescales ', b% s2% xtra1, b% s2% xtra2
     else
        write(*,*) 'something is not going well with the stars IDs '
     end if
     t_sync = 1d0/t_sync
     !write(*,*) 'Hut_rad ', t_sync
  else if (sync_type == "structure_dependent") then !  Checks if the core is radiative or not and uses equation from Hut_con or Hut_rad respectively (Hut word refers to the envelope status)
     !sync_type .eq. "Hut_conv"!Convective envelope + Radiative core
     ! eq. (11) of Hut, P. 1981, A&A, 99, 126
     one_div_t_sync_conv = 3.0d0*k_div_T_posydon(b, s, .true.)*(qratio*qratio/rGyr_squared)*pow6(r_phot/osep)
     !one_div_t_sync_conv2 = 3.0d0*k_div_T_posydon(b, s,2)*(qratio*qratio/rGyr_squared)*pow6(r_phot/osep)
     one_div_t_sync_rad = 3.0d0*k_div_T_posydon(b, s, .false.)*(qratio*qratio/rGyr_squared)*pow6(r_phot/osep)
     !write(*,*) 'star id', s% id
     if (b% point_mass_i /= 1 .and. b% s1% id == s% id) then
        b% s1% xtra(3) = 1d0/one_div_t_sync_rad
        b% s1% xtra(2) = 1d0/one_div_t_sync_conv
        !write(*,*) 'two timescales ', b% s1% xtra1, b% s1% xtra2
     else if (b% point_mass_i /= 2 .and. b% s2% id == s% id) then
        b% s2% xtra(3) = 1d0/one_div_t_sync_rad
        b% s2% xtra(2) = 1d0/one_div_t_sync_conv
        !write(*,*) 'two timescales ', b% s2% xtra1, b% s2% xtra2
     else
        write(*,*) 'something is not going well with the stars IDs '
     end if
     !write(*,*) 'two 1/timescales ', one_div_t_sync_conv , one_div_t_sync_rad
     !write(*,*) 'two timescales ', b% s1% ixtra1, b% s1% ixtra2
     one_div_t_sync = MAX(one_div_t_sync_conv,one_div_t_sync_rad)
     !one_div_t_sync = one_div_t_sync_conv1 + one_div_t_sync_conv2 + one_div_t_sync_rad
     t_sync = 1d0/one_div_t_sync
     !write(*,*) 't_tides in years', t_sync / secyer
  else if (sync_type == "Orb_period") then ! sync on timescale of orbital period
     t_sync = b% period ! synchronize on timescale of orbital period
  else
     ierr = -1
     write(*,*) 'unrecognized sync_type', sync_type
     return
  end if
  t_sync = t_sync / Ftid
end subroutine POSYDON_tsync

subroutine get_tsync(id, sync_type, Ftid, qratio, m, r_phot, osep, t_sync, ierr)
  integer, intent(in) :: id
  character (len=strlen), intent(in) :: sync_type ! synchronization timescale
  real(dp), intent(in) :: Ftid ! efficiency of tidal synchronization. (time scale / Ftid ).
  real(dp), intent(in) :: qratio ! mass_other_star/mass_this_star
  real(dp), intent(in) :: m
  real(dp), intent(in) :: r_phot
  real(dp), intent(in) :: osep ! orbital separation (cm)
  real(dp), intent(out) :: t_sync
  integer, intent(out) :: ierr
  real(dp) :: rGyr_squared, moment_of_inertia
  type (binary_info), pointer :: b
  type (star_info), pointer :: s

  include 'formats'

  ierr = 0

  call star_ptr(id, s, ierr)
  if (ierr /= 0) then
     write(*,*) 'failed in star_ptr'
     return
  end if

  call binary_ptr(s% binary_id, b, ierr)
  if (ierr /= 0) then
     write(*,*) 'failed in binary_ptr'
     return
  end if
  ! calculate the gyration radius squared
  moment_of_inertia = dot_product(s% i_rot(:s% nz), s% dm_bar(:s% nz))
  rGyr_squared = (moment_of_inertia/(m*r_phot*r_phot))
  if (sync_type == "Hut_conv") then
     ! eq. (11) of Hut, P. 1981, A&A, 99, 126
     t_sync = 3.0*k_div_T(b, s, .true.)*(qratio*qratio/rGyr_squared)*pow6(r_phot/osep)
     ! invert it.
     t_sync = 1d0/t_sync
  else if (sync_type == "Hut_rad") then
     ! eq. (11) of Hut, P. 1981, A&A, 99, 126
     t_sync = 3.0*k_div_T(b, s,.false.)*(qratio*qratio/rGyr_squared)*pow6(r_phot/osep)
     ! invert it.
     t_sync = 1d0/t_sync
  else if (sync_type == "Orb_period") then ! sync on timescale of orbital period
     t_sync = b% period ! synchronize on timescale of orbital period
  else
     ierr = -1
     write(*,*) 'unrecognized sync_type', sync_type
     return
  end if
  t_sync = t_sync / Ftid
end subroutine get_tsync

subroutine POSYDON_sync_spin_to_orbit(id, nz, osep, qratio, rl, dt_next, Ftid,sync_type, sync_mode, ierr)
  use const_def, only: dp, strlen
  integer, intent(in) :: id
  integer, intent(in) :: nz
  real(dp), intent(in) :: osep ! orbital separation (cm)
  real(dp), intent(in) :: qratio ! mass_other_star/mass_this_star
  real(dp), intent(in) :: rl ! roche lobe radius (cm)
  real(dp), intent(in) :: dt_next ! next timestep
  real(dp), intent(in) :: Ftid ! efficiency of tidal synchronization. (time scale / Ftid ).

  character (len=strlen), intent(in) :: sync_type ! synchronization timescale
  character (len=strlen), intent(in) :: sync_mode ! where to put/take angular momentum
  integer, intent(out) :: ierr
  type (star_info), pointer :: s
  type (binary_info), pointer :: b

  integer :: k
  real(dp), dimension(nz) :: j_sync, delta_j
  real(dp) :: t_sync, m, r_phot, omega_orb
  real(dp) :: a1,a2

  include 'formats'
  ierr = 0

  t_sync = 0
  call star_ptr(id, s, ierr)
  if (ierr /= 0) then
     write(*,*) 'failed in star_ptr'
     return
  end if

  call binary_ptr(s% binary_id, b, ierr)
  if (ierr /= 0) then
     write(*,*) 'failed in binary_ptr'
     return
  end if

  if (is_donor(b, s)) then
     m = b% m(b% d_i)
     r_phot = b% r(b% d_i)
  else
     m = b% m(b% a_i)
     r_phot = b% r(b% a_i)
  end if

  omega_orb = 2d0*pi/b% period
  do k=1,nz
     j_sync(k) = omega_orb*s% i_rot(k)
  end do

  if (.not. b% use_other_tsync) then !Default tidal synchronization timescale calculation
     call get_tsync(s% id, sync_type, Ftid, qratio, m, r_phot, osep, t_sync, ierr)
     if (ierr/=0) return
  else
     call b% other_tsync(s% id, sync_type, Ftid, qratio, m, r_phot, osep, t_sync, ierr)
     if (ierr/=0) return
  end if
  a1 = f2(b% eccentricity)
  a2 = pow(1-pow2(b% eccentricity), 1.5d0)*f5(b% eccentricity)

  ! Tides apply only to the envelope. (Qin et al. 2018 implementation)
  !if (.not. b% have_radiative_core(id)) then ! convective core
  !    !write(*,*) 'applying tides only in radiative envelope'
  !    do k=1,nz
  !       if (s% mixing_type(k) /= convective_mixing) then
  !           delta_j(k) = (1d0 - exp(-a2*dt_next/t_sync))*(s% j_rot(k) - a1/a2*j_sync(k))
  !       else
  !           delta_j(k) = 0.0
  !       end if
  !    end do
  !else
  !    !write(*,*) 'applying tides only in convective regions'
  !    do k=1,nz
  !       if (s% mixing_type(k) == convective_mixing) then
  !           delta_j(k) = (1d0 - exp(-a2*dt_next/t_sync))*(s% j_rot(k) - a1/a2*j_sync(k))
  !       else
  !           delta_j(k) = 0.0
  !       end if
  !    end do
  !end if

  ! Tides apply in all layers
  ! write(*,*) 'applying tides in all layers'
  do k=1,nz
     delta_j(k) = (1d0 - exp(-a2*dt_next/t_sync))*(s% j_rot(k) - a1/a2*j_sync(k))
  end do


  if (b% point_mass_i /= 1 .and. b% s1% id == s% id) then
     b% t_sync_1 = t_sync
  else
     b% t_sync_2 = t_sync
  end if

  if (.not. b% doing_first_model_of_run) then
     do k=1,nz
        s% extra_jdot(k) = s% extra_jdot(k) - delta_j(k)/dt_next
     end do
  end if
end subroutine POSYDON_sync_spin_to_orbit

real(dp) function f2(e)
  real(dp), intent(in) :: e

  f2 = 1d0

  ! Hut 1981, A&A, 99, 126, definition of f2 after eq. 11
  if (e > 0d0) then
     f2 = 1d0 + 15d0/2d0 * pow2(e) + 45d0/8d0 * pow4(e) + 5d0/16d0 * pow6(e)
  end if

end function f2

real(dp) function f3(e)
  real(dp), intent(in) :: e

  f3 = 1d0

  ! Hut 1981, A&A, 99, 126, definition of f3 after eq. 11
  if (e > 0d0) then
     f3 = 1d0 + 15d0/4d0*pow2(e) + 15d0/8d0 * pow4(e) + 5d0/64d0 * pow6(e)
  end if

end function f3


real(dp) function f4(e)
  real(dp), intent(in) :: e

  f4 = 1d0

  ! Hut 1981, A&A, 99, 126, definition of f4 after eq. 11
  if (e > 0d0) then
     f4 = 1d0 + 3d0/2d0 * pow2(e) + 1d0/8d0 * pow4(e)
  end if

end function f4


real(dp) function f5(e)
  real(dp), intent(in) :: e

  f5 = 1d0

  ! Hut 1981, A&A, 99, 126, definition of f5 after eq. 11
  if (e > 0d0) then
     f5 = 1d0 + 3d0*pow2(e) + 3d0/8d0 * pow4(e)
  end if

end function f5


real(dp) function k_div_T(b, s, has_convective_envelope)
  type(binary_info), pointer :: b
  type(star_info), pointer :: s
  logical, intent(in) :: has_convective_envelope

  integer :: k,i, h1
  real(dp) osep, qratio, m, r_phot,porb, m_env, r_env, tau_conv, P_tid, f_conv,E2, Xs

  ! k/T computed as in Hurley, J., Tout, C., Pols, O. 2002, MNRAS, 329, 897
  ! Kudos to Francesca Valsecchi for help implementing and testing this

  k_div_T = 0d0

  osep = b% separation
  qratio = b% m(b% a_i) / b% m(b% d_i)
  if (is_donor(b, s)) then
     m = b% m(b% d_i)
     r_phot = b% r(b% d_i)
  else
     qratio = 1.0d0/qratio
     m = b% m(b% a_i)
     r_phot = b% r(b% a_i)
  end if
  porb = b% period

  if (has_convective_envelope) then
     m_env = 0d0
     r_env = 0d0
     do k=1, s% nz
        if (s% mixing_type(k) /= convective_mixing .and. &
             s% rho(k) > 1d5*s% rho(1)) then
           r_env = (r_phot - s% r(k))/Rsun
           m_env = (s% m(1) - s% m(k))/Msun
           exit
        end if
     end do
     tau_conv = 0.431d0*pow(m_env*r_env* &
          (r_phot/Rsun-r_env/2d0)/3d0/s% L_phot,one_third) * secyer
     P_tid = 1d0/abs(1d0/porb-s% omega_avg_surf/(2d0*pi))
     f_conv = min(1.0d0, pow(P_tid/(2d0*tau_conv),b% tidal_reduction))

     k_div_T = 2d0/21d0*f_conv/tau_conv*m_env/(m/Msun)
  else ! radiative envelope
     ! New fitting E2 (Qin et al. 2018)
     do i = s% nz, 1, -1
        if (s% brunt_N2(i) >= 0) exit
     end do
     !write(*,*) i
     h1 = s% net_iso(ih1)
     Xs = s% xa(h1,1)
     ! E2 is different for H-rich and He stars (Qin et al. 2018)
     if (Xs < 0.4d0) then ! HeStar
        E2 = exp10(-0.93_dp)*pow(s% r(i)/r_phot,6.7_dp)! HeStars
     else
        E2 = exp10(-0.42_dp)*pow(s% r(i)/r_phot,7.5_dp)! H-rich stars
        !write(*,*) E2, s% r(i)
     end if
     if (isnan(E2)) then  !maybe this won't be used.
        k_div_T = 1d-20
     else
        k_div_T = sqrt(standard_cgrav*m*r_phot*r_phot/pow5(osep)/(Msun/pow3(Rsun)))
        k_div_T = k_div_T*pow(1d0+qratio,5d0/6d0)
        k_div_T = k_div_T * E2
     end if
  end if

end function k_div_T


real(dp) function k_div_T_posydon(b, s, conv_layer_calculation)
  type(binary_info), pointer :: b
  type(star_info), pointer :: s
  !logical, intent(in) :: has_convective_envelope
  logical, intent(in) :: conv_layer_calculation

  integer :: k,i, h1, top_bound_zone, bot_bound_zone
  real(dp) :: osep, qratio, m, r_phot,porb, m_env, Dr_env, Renv_middle, tau_conv, P_tid, f_conv,E2, Xs, m_conv_core
  real(dp) :: k_div_T_posydon_new, conv_mx_top, conv_mx_bot, conv_mx_top_r, conv_mx_bot_r ,omega_conv_region,r_top, r_bottom
  integer :: n_conv_regions_posydon
  integer,  dimension (max_num_mixing_regions) :: n_zones_of_region, bot_bdy, top_bdy
  real(dp), dimension (max_num_mixing_regions) :: cz_bot_mass_posydon
  real(dp) :: cz_bot_radius_posydon(max_num_mixing_regions)
  real(dp), dimension (max_num_mixing_regions) :: cz_top_mass_posydon, cz_top_radius_posydon

  ! k/T computed as in Hurley, J., Tout, C., Pols, O. 2002, MNRAS, 329, 897
  ! Kudos to Francesca Valsecchi for help implementing and testing this

  k_div_T_posydon = 0d0

  osep = b% separation
  qratio = b% m(b% a_i) / b% m(b% d_i)
  if (is_donor(b, s)) then
     m = b% m(b% d_i)
     r_phot = b% r(b% d_i)
  else
     qratio = 1.0d0/qratio
     m = b% m(b% a_i)
     r_phot = b% r(b% a_i)
  end if
  porb = b% period

  if (conv_layer_calculation) then
     m_conv_core = mass_conv_core(s)
     !write(*,'(g0)') 'm_conv_core', m_conv_core
     !write(*,'(g0)') "s% nz, s% n_conv_regions, s% num_conv_boundaries", s% nz, s% n_conv_regions, s% num_conv_boundaries
     !do k_boundary = 1, s% num_conv_boundaries
     !   write(*,'(g0)') "k_boundary, s% conv_bdy_loc(k_boundary), s% top_conv_bdy(k_boundary)", &
     !   k_boundary, s% conv_bdy_loc(k_boundary), s% top_conv_bdy(k_boundary)
     !end do
     !do k=1, s% n_conv_regions
     !   top_k = 2 * s% n_conv_regions -2*(k-1)
     !   bot_k = 2 * s% n_conv_regions -2*(k-1) -1
     !   write(*,'(g0)') "trying with k in n_conv_regions", k, top_k, &
     !   s% top_conv_bdy(top_k),s% r(s% conv_bdy_loc(top_k))/Rsun,&
     !   bot_k, &
     !   s% top_conv_bdy(bot_k), s% r(s% conv_bdy_loc(bot_k))/ Rsun
     !end do
     n_zones_of_region(:)=0
     bot_bdy(:)=0
     top_bdy(:)=0
     cz_bot_mass_posydon(:)=0.0d0
     cz_bot_radius_posydon(:)=0.0d0
     cz_top_mass_posydon(:)=0.0d0
     cz_top_radius_posydon(:)=0.0d0
     n_conv_regions_posydon = 0

     call loop_conv_layers(s,n_conv_regions_posydon, n_zones_of_region, bot_bdy, top_bdy, &
          cz_bot_mass_posydon, cz_bot_radius_posydon, cz_top_mass_posydon, cz_top_radius_posydon)

     if (n_conv_regions_posydon > 0) then
        do k=1, n_conv_regions_posydon ! from inside out
           m_env = 0.0d0
           Dr_env = 0.0d0
           Renv_middle = 0.0d0
           if ((cz_bot_mass_posydon(k) / Msun) >=  m_conv_core) then ! if the conv. region is not inside the conv. core
              m_env = (cz_top_mass_posydon(k) - cz_bot_mass_posydon(k)) / Msun
              Dr_env = cz_top_radius_posydon(k) - cz_bot_radius_posydon(k)  !depth of the convective layer, length of the eddie
              ! Corresponding to the Renv term in eq.31 of Hurley et al. 2002
              ! and to (R-Renv) term in eq. 4 of Rasio et al. 1996  (different notation)
              Renv_middle = (cz_top_radius_posydon(k) + cz_bot_radius_posydon(k) )*0.5d0  !middle of the convective layer
              ! Corresponding to the (R-0.5d0*Renv) in eq.31 of Hurley et al 2002
              ! and to the Renv in eq. 4 of Rasio et al. 1996
              ! where it represented the base of the convective layer (different notation)
              tau_conv = 0.431_dp*pow(m_env*Dr_env* &
                   Renv_middle/3d0/s% L_phot,1.0d0/3.0d0) * secyer
              P_tid = 1d0/abs(1d0/porb-s% omega(top_bdy(k))/(2d0*pi))
              f_conv = min(1.0d0, pow(P_tid/(2d0*tau_conv), b% tidal_reduction))
              !write(*,'(g0)') 'porb, p_from_omega, f_conv = ', porb, &
              !1                / (s% omega(top_bdy(k))/(2d0*pi)), &
              !1                /(s% omega_avg_surf/(2d0*pi)), f_conv
              k_div_T_posydon_new = 2d0/21d0*f_conv/tau_conv*m_env/(m/Msun)
              !write(*,'(g0)') 'tau_conv, K/T = ', tau_conv, k_div_T_posydon_new, m_env, (m/Msun)
              if (k_div_T_posydon_new >= k_div_T_posydon) then
                 k_div_T_posydon = k_div_T_posydon_new
                 !write(*,'(g0)') 'M_env, DR_env, Renv_middle, omega_conv_region in conv region ', k ,' is ', &
                 ! m_env, Dr_env, Renv_middle, s% omega(top_bdy(k)), 'spanning number of zones = ', &
                 !top_bdy(k) , bot_bdy(k), &
                 !n_zones_of_region(k)
              end if
           end if
        end do
     end if
  else ! assuming a radiative star
     ! New fitting E2 (Qin et al. 2018)
     do i = s% nz, 1, -1
        if (s% brunt_N2(i) >= 0d0) exit
     end do
     !write(*,*) i
     if (i == 0) then ! expected in a fully convective star
        E2 = 1d-99
     else
        h1 = s% net_iso(ih1)
        Xs = s% xa(h1,1)
        ! E2 is different for H-rich and He stars (Qin et al. 2018)
        if (Xs < 0.4d0) then ! HeStar
           E2 = exp10(-0.93_dp)*pow(s% r(i)/r_phot, 6.7_dp)! HeStars
        else
           E2 = exp10(-0.42_dp)*pow(s% r(i)/r_phot, 7.5_dp)! H-rich stars
           !write(*,*) E2, s% r(i)
        end if
     end if

     if (isnan(E2)) then  !maybe this won't be used.
        k_div_T_posydon = 1d-99
     else
        k_div_T_posydon = sqrt(standard_cgrav*m*r_phot*r_phot/pow5(osep)/(Msun/pow3(Rsun)))
        k_div_T_posydon = k_div_T_posydon*pow(1d0+qratio,5d0/6d0)
        k_div_T_posydon = k_div_T_posydon * E2
     end if
  end if

end function k_div_T_posydon




real(dp) function acc_radius(b, m_acc) !Calculates Sch. radius of compact object (or surface radius in case of NS) in cm
  type(binary_info), pointer :: b
  real(dp) :: m_acc, a
  real(dp) :: r_isco, Z1, Z2, eq_initial_bh_mass

  if (m_acc/Msun <= 2.50d0) then ! NS
     !Radius refernces for NS:
     ! 1) Miller, M. C., Lamb, F. K., Dittmann, A. J., et al. 2019, ApJL, 887, L2
     ! 2) Riley, T. E., Watts, A. L., Bogdanov, S., et al., 2019, ApJL, 887, L21
     ! 3) Landry, P., Essick, R., & Chatziioannou, K. 2020
     ! 4) E.R. Most, L.R. Weih, L. Rezzolla and J. Schaffner-Bielich, 2018, Phys. Rev. Lett. 120, 261103
     ! 5) Abbott, B. P., Abbott, R., Abbott, T. D., et al. 2020, ApJL, 892, L3
     acc_radius = 12.5E5_dp !* 10 ** 5 !in cm
  else ! Event horizon for Kerr-BH
     ! this part is only relevant for BH accretors
     if (b% initial_bh_spin < 0d0) then
        b% initial_bh_spin = 0d0
        write(*,*) "initial_bh_spin is smaller than zero. It has been set to zero."
     else if (b% initial_bh_spin > 1d0) then
        b% initial_bh_spin = 1d0
        write(*,*) "initial_bh_spin is larger than one. It has been set to one."
     end if
     ! compute isco radius from eq. 2.21 of Bardeen et al. (1972), ApJ, 178, 347
     Z1 = 1d0 + pow(1d0 - pow2(b% initial_bh_spin),one_third) &
          * (pow(1d0 + b% initial_bh_spin,one_third) + pow(1d0 - b% initial_bh_spin,one_third))
     Z2 = sqrt(3d0*pow2(b% initial_bh_spin) + pow2(Z1))
     r_isco = 3d0 + Z2 - sqrt((3d0 - Z1)*(3d0 + Z1 + 2d0*Z2))
     ! compute equivalent mass at zero spin from eq. (3+1/2) (ie. the equation between (3) and (4))
     ! of Bardeen (1970), Nature, 226, 65, taking values with subscript zero to correspond to
     ! zero spin (r_isco = sqrt(6)).

     if (initial_mass(2) > 2.5_dp) then ! If it was already a BH then take the initial mass m2
        eq_initial_bh_mass = b% eq_initial_bh_mass
     else if (initial_mass(2) <= 2.5_dp) then! If it was initially a NS then take 2.5Msun as eq_initial_mass
        eq_initial_bh_mass = 2.5_dp * Msun * sqrt(r_isco/6d0)
     end if

     a = sqrt(two_thirds) &
          *(eq_initial_bh_mass/min(b% m(b% point_mass_i),sqrt(6d0)* eq_initial_bh_mass)) &
          *(4._dp - sqrt(18._dp*pow2(eq_initial_bh_mass/ &
          min(b% m(b% point_mass_i),sqrt(6d0)* eq_initial_bh_mass)) - 2._dp))
     !Podsiadlowski et al. (2003) assuming a initially non-rotating BH
     acc_radius = (1.0_dp + sqrt(1.0_dp - a*a)) * b% s_donor% cgrav(1) * m_acc/pow2(clight)
  end if
end function acc_radius



! Note that the function below is unused and untested by MAthieu

!! Eddington accreton limits for NS and BH
subroutine POSYDON_mdot_edd(binary_id, mdot_edd, mdot_edd_eta, ierr)
  use const_def, only: dp
  integer, intent(in) :: binary_id
  real(dp), intent(out) :: mdot_edd
  real(dp), intent(out) :: mdot_edd_eta
  integer, intent(out) :: ierr
  real(dp) :: r_isco, Z1, Z2, eq_initial_bh_mass
  type (binary_info), pointer :: b
  ierr = 0
  call binary_ptr(binary_id, b, ierr)
  if (ierr /= 0) then
     write(*,*) 'failed in binary_ptr'
     return
  end if
  if (b% m(2)/Msun > 2.50_dp) then ! M2 > 2.5 Msol for BHs
     ! this part is only relevant for BH accretors
     if (b% initial_bh_spin < 0d0) then
        b% initial_bh_spin = 0d0
        write(*,*) "initial_bh_spin is smaller than zero. It has been set to zero."
     else if (b% initial_bh_spin > 1d0) then
        b% initial_bh_spin = 1d0
        write(*,*) "initial_bh_spin is larger than one. It has been set to one."
     end if
     ! compute isco radius from eq. 2.21 of Bardeen et al. (1972), ApJ, 178, 347
     Z1 = 1d0 + pow(1d0 - pow2(b% initial_bh_spin),one_third) &
          * (pow(1d0 + b% initial_bh_spin,one_third) + pow(1d0 - b% initial_bh_spin,one_third))
     Z2 = sqrt(3d0*pow2(b% initial_bh_spin) + pow2(Z1))
     r_isco = 3d0 + Z2 - sqrt((3d0 - Z1)*(3d0 + Z1 + 2d0*Z2))
     ! compute equivalent mass at zero spin from eq. (3+1/2) (ie. the equation between (3) and (4))
     ! of Bardeen (1970), Nature, 226, 65, taking values with subscript zero to correspond to
     ! zero spin (r_isco = sqrt(6)).

     if (initial_mass(2) > 2.5_dp) then ! If it was already a BH then take the initial mass m2
        eq_initial_bh_mass = b% eq_initial_bh_mass
     else if (initial_mass(2) <= 2.5_dp) then! If it was initially a NS then take 2.5 as eq_initial_mass
        eq_initial_bh_mass = 2.5_dp * Msun * sqrt(r_isco/6d0)
     end if

     !! mdot_edd_eta for BH following Podsiadlowski, Rappaport & Han (2003), MNRAS, 341, 385
     mdot_edd_eta = 1d0 - sqrt(1d0 - &
          pow2(min(b% m(b% a_i),sqrt(6d0)*eq_initial_bh_mass)/(3d0*eq_initial_bh_mass)))
  else ! NS
     !! mdot_edd_eta for NS accretors
     mdot_edd_eta = b% s_donor% cgrav(1) * b% m(2) / (pow2(clight) * acc_radius(b, b% m(2)))
  end if
  mdot_edd = 4d0*pi*b% s_donor% cgrav(1)*b% m(b% a_i) &
       /(clight*0.2d0*(1d0+b% s_donor% surface_h1)*mdot_edd_eta)
  !b% s1% x_ctrl(1) used to adjust the Eddington limit in inlist1
  mdot_edd = mdot_edd * b% s1% x_ctrl(1)

end subroutine POSYDON_mdot_edd

! Include functions for single stars
include 'POSYDON_single_stars.inc'


! real(dp) function eval_rlobe(m1, m2, a) result(rlobe)
!   real(dp), intent(in) :: m1, m2, a
!   real(dp) :: q
!   q = pow(m1/m2,one_third)
!   ! Roche lobe size for star of mass m1 with a
!   ! companion of mass m2 at separation a, according to
!   ! the approximation of Eggleton 1983, apj 268:368-369
!   rlobe = a*0.49d0*q*q/(0.6d0*q*q + log10(q+1))
! end function eval_rlobe
