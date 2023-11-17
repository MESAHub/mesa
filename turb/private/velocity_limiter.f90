module velocity_limiter

    implicit none
  
    private
  
    public :: calculate_drag_limited_velocity
  
  contains
  
    subroutine calculate_drag_limited_velocity(radius, gravity, grad_mu, phi ,v_max)
  
      real(dp), intent(in) :: radius
      real(dp), intent(in) :: gravity 
      real(dp), intent(in) :: grad_mu
      real(dp), intent(in) :: phi
      real(dp), intent(out) :: v_max
  
      real(dp) :: scale_height
      real(dp) :: delta_rho
      real(dp) :: buoyancy_force, drag_force
  
      ! Balance buoyancy force with drag 
      ! Assume isothermal blob (limit where thermohaline mixing -> Rayleigh-Taylor)

  
      v_max = sqrt(radius * gravity * phi * grad_mu) 
  
    end subroutine calculate_drag_limited_velocity
    
  end module velocity_limiter