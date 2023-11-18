module velocity_limiter

    use const_def
    use num_lib
    use math_lib
    use utils_lib

    implicit none
  
    private
  
    public :: calculate_drag_limited_velocity
  
  contains

  ! Calling the limiter requires knowledge of the radius of the falling element
  ! In the Rayleigh-Taylor limit of thermohaline mixing the horizontal wavenumber of the 
  ! finger instability sets the relevant scale. We do calculate this in fingering_modes.f90 
  ! 
  ! Gravity can be calculated as   
  ! G = s% cgrav(k)
  ! grav = G*s% m_grav(k)/pow2(s% r(k))
  ! 
  ! grad_mu is s% gradL_composition_term(k) 
  ! 
  ! phi = 1 for ideal gas
  !
  ! Q parameter = s% chiT(k) / s% chiRho(k) = dlnRho/dlnT = \delta (see C&G 14.24)
  ! Required to go from drho/rho = \phi dmu/mu -> Grad_mu 
  !  
  ! Returns v_max (cgs) 

  ! Example of call in run_star_extras: 
  
    ! names(1) = 'log_v_max'
    !   v_max = 0d0 

    !   do k=1,s% nz
    !      grav = s% cgrav(k)*s% m_grav(k)/pow2(s% r(k))
    !      radius = 1d4 ! Set it to typical horizontal scale of thermohaline fingers in FRG23 for now (100m is typical in RGB, see HG19 table 1)
    !      grad_mu = s% gradL_composition_term(k) 
    !      Q = s% chiT(k) / s% chiRho(k) 
    !      call calculate_drag_limited_velocity(radius, grav, grad_mu, Q ,v_max)
    !      vals(k,1) = safe_log10(v_max)
    !   end do

  
  subroutine calculate_drag_limited_velocity(radius, gravity, grad_mu, Q, v_max)
  
    real(dp), intent(in) :: radius
    real(dp), intent(in) :: gravity 
    real(dp), intent(in) :: grad_mu
    real(dp), intent(in) :: Q 
    real(dp), intent(out) :: v_max

    ! Balance buoyancy force with drag 
    ! Assume isothermal blob (limit where thermohaline mixing -> Rayleigh-Taylor)
    ! Only calculate in thermohaline unstable regions. Otherwise make it 0 

    if (Q*grad_mu < 0d0) then 
       v_max = sqrt(radius * gravity * abs(Q*grad_mu)) 
    else 
       v_max = 0d0  
    endif 

  end subroutine calculate_drag_limited_velocity
  
    
  end module velocity_limiter