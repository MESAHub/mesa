! Contains the custom MESA routines for POSYDON integration for single stars
! author: Mathieu Renzo - copy paste of POSYDON repo and adapted to MESA r15140

function f_ov_fcn_of_mass(m) result(f_ov)
  use const_def, only: pi
  real(dp), intent(in) :: m
  ! real(dp), parameter :: pi = 3.1415926535897932384626433832795028841971693993751d0
  real(dp) :: f_ov, frac
  real(dp), parameter :: f1 = 1.6d-2, f2=4.15d-2
  if(m < 4.0d0) then
     frac = 0.0d0
  elseif(m >= 8.0d0) then
     frac = 1.0d0
  else
     frac = 0.5d0 * (1.0d0 - cos(pi*0.25d0*(m-4.0d0)))
  endif
  f_ov = f1 + (f2-f1)*frac
end function f_ov_fcn_of_mass

real(dp) function mass_conv_core(s)
  type (star_info), pointer :: s
  integer :: j, nz, k
  real(dp) :: dm_limit
  include 'formats'
  mass_conv_core = 0.0d0
  dm_limit = s% conv_core_gap_dq_limit*s% xmstar
  nz = s% nz
  do j = 1, s% n_conv_regions
     ! ignore possible small gap at center
     if (s% cz_bot_mass(j) <= s% m(nz) + dm_limit) then
        mass_conv_core = s% cz_top_mass(j)/Msun
        ! jump over small gaps
        do k = j+1, s% n_conv_regions
           if (s% cz_bot_mass(k) - s% cz_top_mass(k-1) >= dm_limit) exit
           mass_conv_core = s% cz_top_mass(k)/Msun
        end do
        exit
     end if
  end do
end function mass_conv_core


real(dp) function lambda_CE(s, adjusted_energy, star_core_mass_CE)
  type (star_info), pointer :: s
  integer :: k
  real(dp) :: E_bind, E_bind_shell, star_core_mass_CE
  real(dp) :: adjusted_energy(:)

  if (s% m(1) <= (star_core_mass_CE)) then
     lambda_CE = 1d99 ! no envelope, so immediately have a "succesfull envelope ejection"
  else
     E_bind = 0.0d0
     E_bind_shell = 0.0d0
     do k=1, s% nz
        if (s% m(k) > (star_core_mass_CE)) then !envelope is defined as everything above star_core_mass_CE.
           E_bind_shell = s% dm(k) * adjusted_energy(k) - (s% cgrav(1) * s% m(k) * s% dm_bar(k))/(s% r(k))
           E_bind = E_bind+ E_bind_shell
        end if
     end do
     lambda_CE = - s% cgrav(1) * (s% m(1)) * ((s% m(1)) - star_core_mass_CE)/(E_bind * s% r(1))
  end if
end function lambda_CE


real(dp) function get_ion_info(s,id,k)
  use ionization_def, only: num_ion_vals
  use ionization_lib, only: eval_ionization
  integer, intent(in) :: id, k
  integer :: ierr
  real(dp) :: ionization_res(num_ion_vals)
  type (star_info), pointer :: s
  ierr = 0
  call eval_ionization( &
       1d0 - (s% X(k) + s% Y(k)), s% X(k), s% Rho(k), s% lnd(k)/ln10, &
       s% T(k), s% lnT(k)/ln10, ionization_res, ierr)
  if (ierr /= 0) ionization_res = 0
  get_ion_info = ionization_res(id)
end function get_ion_info


! simpler version of the same function at star/private/report.f90
subroutine set_core_info(s, k, &
     core_k, core_m, core_r)
  type (star_info), pointer :: s
  integer, intent(in) :: k
  integer, intent(out) :: core_k
  real(dp), intent(out) :: &
       core_m, core_r

  integer :: j, jm1, j00
  real(dp) :: dm1, d00, qm1, q00, core_q, &
       core_lgP, core_g, core_X, core_Y, core_edv_H, core_edv_He, &
       core_scale_height, core_dlnX_dr, core_dlnY_dr, core_dlnRho_dr

  include 'formats'

  if (k == 1) then
     core_q = 1d0
  else
     jm1 = maxloc(s% xa(:,k-1), dim=1)
     j00 = maxloc(s% xa(:,k), dim=1)
     qm1 = s% q(k-1) - 0.5d0*s% dq(k-1) ! center of k-1
     q00 = s% q(k) - 0.5d0*s% dq(k) ! center of k
     dm1 = s% xa(j00,k-1) - s% xa(jm1,k-1)
     d00 = s% xa(j00,k) - s% xa(jm1,k)
     if (dm1*d00 > 0d0) then
        write(*,2) 'bad args for set_core_info', k, dm1, d00
        ! call mesa_error(__FILE__,__LINE__)
        core_q = 0.5d0*(qm1 + q00)
     else if (dm1 == 0d0 .and. d00 == 0d0) then
        core_q = 0.5d0*(qm1 + q00)
     else if (dm1 == 0d0) then
        core_q = qm1
     else if (d00 == 0d0) then
        core_q = q00
     else
        core_q = find0(qm1, dm1, q00, d00)
     end if
  end if

  call get_info_at_q(s, core_q, &
       core_k, core_m, core_r)

end subroutine set_core_info


subroutine get_info_at_q(s, bdy_q, &
     kbdy, bdy_m, bdy_r)

  type (star_info), pointer :: s
  real(dp), intent(in) :: bdy_q
  integer, intent(out) :: kbdy
  real(dp), intent(out) :: &
       bdy_m, bdy_r

  real(dp) :: x, x0, x1, x2, alfa, beta, bdy_omega_crit
  integer :: k, ii, klo, khi

  include 'formats'

  bdy_m=0.0d0; bdy_r=0.0d0;
  kbdy = 0

  if (bdy_q <= 0.0d0) return
  k = k_for_q(s,bdy_q)
  if (k >= s% nz) then
     kbdy = s% nz
     return
  end if
  if (k <= 1) then
     bdy_m = s% star_mass
     bdy_r = s% r(1)/Rsun
     kbdy = 1
     return
  end if

  kbdy = k+1

  bdy_m = (s% M_center + s% xmstar*bdy_q)/Msun

  x = s% q(k-1) - bdy_q
  x0 = s% dq(k-1)/2.0d0
  x1 = s% dq(k)/2 + s% dq(k-1)
  x2 = s% dq(k+1)/2.0d0 + s% dq(k) + s% dq(k-1)

  alfa = max(0d0, min(1d0, (bdy_q - s% q(k+1))/s% dq(k)))

  bdy_r = pow( &
       interp2(s% r(k)*s% r(k)*s% r(k), s% r(k+1)*s% r(k+1)*s% r(k+1)),1d0/3d0)/Rsun

contains

  real(dp) function interp2(f0, f1)
    real(dp), intent(in) :: f0, f1
    interp2 = alfa*f0 + (1.0d0-alfa)*f1
  end function interp2

end subroutine get_info_at_q




integer function k_for_q(s, q)
  ! return k s.t. q(k) >= q > q(k)-dq(k)
  type (star_info), pointer :: s
  real(dp), intent(in) :: q
  integer :: k, nz
  nz = s% nz
  if (q >= 1.0d0) then
     k_for_q = 1; return
  else if (q <= s% q(nz)) then
     k_for_q = nz; return
  end if
  do k = 1, nz-1
     if (q > s% q(k+1)) then
        k_for_q = k; return
     end if
  end do
  k_for_q = nz
end function k_for_q



subroutine loop_conv_layers(s,n_conv_regions_posydon, n_zones_of_region, bot_bdy, top_bdy, &
     cz_bot_mass_posydon, cz_bot_radius_posydon, cz_top_mass_posydon, cz_top_radius_posydon)
  type (star_info), pointer :: s
  ! integer, intent(out) :: ierr

  logical :: in_convective_region
  integer :: k, j, nz
  logical, parameter :: dbg = .false.
  integer, intent(out) :: n_conv_regions_posydon
  !integer :: max_num_mixing_regions
  !max_num_mixing_regions = 100
  !integer, intent(out), dimension (:), allocatable :: n_zones_of_region, bot_bdy, top_bdy
  !real(dp),intent(out), dimension (:), allocatable :: cz_bot_mass_posydon, cz_bot_radius_posydon
  !real(dp),intent(out), dimension (:), allocatable :: cz_top_mass_posydon, cz_top_radius_posydon
  integer :: min_zones_for_convective_tides
  integer ::  pot_n_zones_of_region, pot_bot_bdy, pot_top_bdy
  real(dp) :: pot_cz_bot_mass_posydon, pot_cz_bot_radius_posydon
  integer, intent(out), dimension (max_num_mixing_regions) :: n_zones_of_region, bot_bdy, top_bdy
  real(dp),intent(out), dimension (max_num_mixing_regions) :: cz_bot_mass_posydon
  real(dp),intent(out) :: cz_bot_radius_posydon(max_num_mixing_regions)
  real(dp),intent(out), dimension (max_num_mixing_regions) :: cz_top_mass_posydon, cz_top_radius_posydon

  include 'formats'
  !ierr = 0
  min_zones_for_convective_tides = 10
  nz = s% nz
  n_zones_of_region(:)=0
  bot_bdy(:)=0
  top_bdy(:)=0
  cz_bot_mass_posydon(:)=0.0d0
  cz_bot_radius_posydon(:)=0.0d0
  cz_top_mass_posydon(:)=0.0d0
  cz_top_radius_posydon(:)=0.0d0
  n_conv_regions_posydon = 0
  pot_cz_bot_mass_posydon = 0.0d0
  pot_cz_bot_radius_posydon = 0.0d0
  pot_bot_bdy = 0
  pot_n_zones_of_region = 0

  in_convective_region = (s% mixing_type(nz) == convective_mixing)
  if (in_convective_region) then
     pot_cz_bot_mass_posydon = s% M_center
     pot_cz_bot_radius_posydon = 0.0d0
     pot_bot_bdy = nz
  end if

  !write(*,*) 'initial in_convective_region', in_convective_region

  do k=nz-1, 2, -1
     if (in_convective_region) then
        if (s% mixing_type(k) /= convective_mixing) then ! top of convective region
           pot_top_bdy = k
           pot_n_zones_of_region = pot_bot_bdy - pot_top_bdy
           if (pot_n_zones_of_region >= min_zones_for_convective_tides) then
              if (n_conv_regions_posydon < max_num_mixing_regions) then
                 n_conv_regions_posydon = n_conv_regions_posydon + 1
              end if
              cz_top_mass_posydon(n_conv_regions_posydon) = &
                   s% M_center + (s% q(k) - s% cz_bdy_dq(k))*s% xmstar
              cz_bot_mass_posydon(n_conv_regions_posydon) = pot_cz_bot_mass_posydon
              cz_top_radius_posydon(n_conv_regions_posydon) = s% r(k)/Rsun
              cz_bot_radius_posydon(n_conv_regions_posydon) = pot_cz_bot_radius_posydon
              top_bdy(n_conv_regions_posydon) = pot_top_bdy
              bot_bdy(n_conv_regions_posydon) = pot_bot_bdy
              n_zones_of_region(n_conv_regions_posydon) = pot_n_zones_of_region
           end if
           in_convective_region = .false.
        end if
     else
        if (s% mixing_type(k) == convective_mixing) then ! bottom of convective region
           pot_cz_bot_mass_posydon = &
                s% M_center + (s% q(k) - s% cz_bdy_dq(k))*s% xmstar
           pot_cz_bot_radius_posydon = s% r(k)/Rsun
           pot_bot_bdy = k
           in_convective_region = .true.
        end if
     end if
  end do
  if (in_convective_region) then
     pot_top_bdy = 1
     pot_n_zones_of_region = pot_bot_bdy - pot_top_bdy
     if (pot_n_zones_of_region >= min_zones_for_convective_tides) then
        if (n_conv_regions_posydon < max_num_mixing_regions) then
           n_conv_regions_posydon = n_conv_regions_posydon + 1
        end if
        cz_top_mass_posydon(n_conv_regions_posydon) = s% mstar
        cz_top_radius_posydon(n_conv_regions_posydon) = s% r(1)/Rsun
        top_bdy(n_conv_regions_posydon) = 1
        cz_bot_mass_posydon(n_conv_regions_posydon) = pot_cz_bot_mass_posydon
        cz_bot_radius_posydon(n_conv_regions_posydon) = pot_cz_bot_radius_posydon
        bot_bdy(n_conv_regions_posydon) = pot_bot_bdy
        n_zones_of_region(n_conv_regions_posydon) = pot_n_zones_of_region
     end if
  end if

  !write(*,*)
  !write(*,2) 'set_mixing_info n_conv_regions_posydon', n_conv_regions_posydon
  !do j = 1, n_conv_regions_posydon
  !   write(*,2) 'conv region', j, cz_bot_mass_posydon(j)/Msun, cz_top_mass_posydon(j)/Msun
  !   write(*,2) 'conv region', j, cz_bot_radius_posydon(j), cz_top_radius_posydon(j)
  !end do
  !write(*,*)
end subroutine loop_conv_layers
