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
  ! B is s% gradL_composition_term(k) 
  !
  ! B = (phi/delta)*grad_mu
  ! phi = delta = 1 for ideal gas, so this often reduces to B = grad_mu
  !
  ! delta parameter = s% chiT(k) / s% chiRho(k) = -dlnRho/dlnT
  ! so delta*B = phi*grad_mu ~ drho/rho
  !  
  ! Returns v_max (cgs) 

  ! Example of call in run_star_extras: 
  
    ! names(1) = 'log_v_max'
    !   v_max = 0d0 

    !   do k=1,s% nz
    !      grav = s% cgrav(k)*s% m_grav(k)/pow2(s% r(k))
    !      radius = 1d4 ! Set it to typical horizontal scale of thermohaline fingers in FRG23 for now (100m is typical in RGB, see HG19 table 1)
    !      B = s% gradL_composition_term(k) 
    !      delta = s% chiT(k) / s% chiRho(k) 
    !      call calculate_drag_limited_velocity(radius, grav, B, delta ,v_max)
    !      vals(k,1) = safe_log10(v_max)
    !   end do

  
  subroutine calculate_drag_limited_velocity(radius, gravity, B, delta, v_max)
  
    real(dp), intent(in) :: radius
    real(dp), intent(in) :: gravity 
    real(dp), intent(in) :: B
    real(dp), intent(in) :: delta 
    real(dp), intent(out) :: v_max

    ! Balance buoyancy force with drag 
    ! Assume isothermal blob (limit where thermohaline mixing -> Rayleigh-Taylor)
    ! Only calculate in thermohaline unstable regions. Otherwise make it 0 

    ! B ~ grad_mu
    if (delta*B < 0d0) then 
       v_max = sqrt(radius * gravity * abs(delta*B)) 
    else 
       v_max = 0d0  
    endif 

  end subroutine calculate_drag_limited_velocity
  
    
  end module velocity_limiter
