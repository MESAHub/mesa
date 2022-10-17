!!% This module contains code to interpolate gravity darkening coefficients.
!!% The coefficients scale surface-averaged Teff and Luminosity as a function
!!% of the orientation of the observer w.r.t. the rotation axis of the star
!!% and the ratio of the surface angular velocity to the Keplerian angular
!!% velocity. The orientation angle is denoted "inclination" and the angular
!!% velocity is denoted "omega". Permissible ranges are
!!%
!!%         0 (no rotation) <= omega <= 1 (critical rotation)
!!%         0 (equator) <= inclination in radians <= pi/2 (pole)
!!% 
!!%
!!% The coefficients are obtained via 2D interpolation using the interp_2d
!!% module in tables. The tables were computed using the code at
!!%
!!% https://github.com/aarondotter/GDit
!!%
!!% where the reader can also find documentation on how the calculations are
!!% performed. The table included in mesa/data/star_data consists of a square
!!% array with 21 values each of omega and inclination, equally spaced between
!!% the limits quoted above. The interpolated results fall with roughly 0.001
!!% of exact values calculated using the code referenced above over the full
!!% range.

module gravity_darkening
   
   use interp_2d_lib_db
   use utils_lib
   
   implicit none
   
   logical :: GD_initialized = .false.
   
   integer, parameter :: ndim = 21
   integer, parameter :: size_per_point = 4
   
   real(dp) :: omega_grid(ndim), inclination_grid(ndim)
   real(dp), target :: C_T_ary(size_per_point * ndim * ndim)
   real(dp), target :: C_L_ary(size_per_point * ndim * ndim)
   real(dp), pointer :: C_T(:, :, :), C_T1(:), C_L(:, :, :), C_L1(:)
   character(len = 256) :: coefficient_filename
   
   integer :: ibcxmin, ibcxmax, ibcymin, ibcymax
   real(dp) :: bcxmin(ndim), bcxmax(ndim), bcymin(ndim), bcymax(ndim)
   
   private
   public :: gravity_darkening_Teff_coeff, gravity_darkening_L_coeff

contains
   
   subroutine GD_init(ierr)
      use const_def, only : mesa_data_dir
      integer, intent(out) :: ierr
      real(dp) :: dummy(4)
      integer :: io, i, j, ilinx, iliny
      
      C_T1 => C_T_ary
      C_T(1:size_per_point, 1:ndim, 1:ndim) => C_T_ary(1:size_per_point * ndim * ndim)
      
      C_L1 => C_L_ary
      C_L(1:size_per_point, 1:ndim, 1:ndim) => C_L_ary(1:size_per_point * ndim * ndim)
      
      ibcxmin = 0; bcxmin = 0
      ibcxmax = 0; bcxmax = 0
      ibcymin = 0; bcymin = 0
      ibcymax = 0; bcymax = 0
      
      coefficient_filename = trim(mesa_data_dir) // '/star_data/gravity_darkening_coefficients.data'
      
      open(newunit = io, file = trim(coefficient_filename), status = 'old', action = 'read')
      read(io, *) !skip header
      do i = 1, ndim
         do j = 1, ndim
            read(io, *) dummy(1:4)
            if(i==1) inclination_grid(j) = dummy(2)
            if(j==1) omega_grid(i) = dummy(1)
            C_T(1, i, j) = dummy(3)
            C_L(1, i, j) = dummy(4)
         enddo
      enddo
      close(io)
      
      ! construct interpolant for C_T
      call interp_mkbicub_db(omega_grid, ndim, inclination_grid, ndim, C_T1, ndim, &
         ibcxmin, bcxmin, ibcxmax, bcxmax, ibcymin, bcymin, ibcymax, bcymax, &
         ilinx, iliny, ierr)
      
      !construct interpolant for C_L
      call interp_mkbicub_db(omega_grid, ndim, inclination_grid, ndim, C_L1, ndim, &
         ibcxmin, bcxmin, ibcxmax, bcxmax, ibcymin, bcymin, ibcymax, bcymax, &
         ilinx, iliny, ierr)
      
      if(ierr==0) GD_initialized = .true.
   end subroutine GD_init
   
   
   function GD_coeff(omega, inclination, C1) result(coeff)
      use const_def, only : pi
      real(dp), intent(in) :: omega, inclination
      real(dp), pointer :: C1(:)
      real(dp) :: coeff, coeff_eval(6), safe_omega, safe_incl
      integer :: ict(6) = [1, 0, 0, 0, 0, 0], ierr
      if(.not.GD_initialized) call GD_init(ierr)
      !ensure that omega and inclination are within table bounds
      safe_omega = min(max(omega, 0.0d0), 1.0d0)
      safe_incl = min(max(inclination, 0.0d0), 0.5d0 * pi)
      call interp_evbicub_db(safe_omega, safe_incl, omega_grid, ndim, inclination_grid, ndim, &
         1, 1, C1, ndim, ict, coeff_eval, ierr)
      if(ierr==0)then
         coeff = coeff_eval(1)
      else
         coeff = 1.0d0
      endif
   end function GD_coeff
   
   
   function gravity_darkening_Teff_coeff(omega, inclination) result(Teff_coeff)
      real(dp), intent(in) :: omega, inclination
      real(dp) :: Teff_coeff
      Teff_coeff = GD_coeff(omega, inclination, C_T1)
   end function gravity_darkening_Teff_coeff
   
   
   function gravity_darkening_L_coeff(omega, inclination) result(L_coeff)
      real(dp), intent(in) :: omega, inclination
      real(dp) :: L_coeff
      L_coeff = GD_coeff(omega, inclination, C_L1)
   end function gravity_darkening_L_coeff


end module gravity_darkening


