! ***********************************************************************
!
!   Copyright (C) 2026  The MESA Team
!
! ***********************************************************************

module pgstar_sed

   use star_private_def
   use const_def, only: dp, strlen
   use pgstar_support
   use star_pgstar

   implicit none

contains

   subroutine do_SED_Plot1(s, id, device_id, xleft, xright, ybot, ytop, subplot, title, txt_scale, ierr)
      type(star_info), pointer :: s; integer :: id, device_id, ierr; real :: xleft, xright, ybot, ytop, txt_scale; logical :: subplot; character(len=*) :: title
      call render_sed_core(s, id, device_id, xleft, xright, ybot, ytop, subplot, title, txt_scale, &
         s% pg% SED1_lambda_min, s% pg% SED1_lambda_max, s% pg% SED1_flux_min, s% pg% SED1_flux_max, ierr)
   end subroutine do_SED_Plot1

   subroutine do_SED_Plot2(s, id, device_id, xleft, xright, ybot, ytop, subplot, title, txt_scale, ierr)
      type(star_info), pointer :: s; integer :: id, device_id, ierr; real :: xleft, xright, ybot, ytop, txt_scale; logical :: subplot; character(len=*) :: title
      call render_sed_core(s, id, device_id, xleft, xright, ybot, ytop, subplot, title, txt_scale, &
         s% pg% SED2_lambda_min, s% pg% SED2_lambda_max, s% pg% SED2_flux_min, s% pg% SED2_flux_max, ierr)
   end subroutine do_SED_Plot2

   subroutine do_SED_Plot3(s, id, device_id, xleft, xright, ybot, ytop, subplot, title, txt_scale, ierr)
      type(star_info), pointer :: s; integer :: id, device_id, ierr; real :: xleft, xright, ybot, ytop, txt_scale; logical :: subplot; character(len=*) :: title
      call render_sed_core(s, id, device_id, xleft, xright, ybot, ytop, subplot, title, txt_scale, &
         s% pg% SED3_lambda_min, s% pg% SED3_lambda_max, s% pg% SED3_flux_min, s% pg% SED3_flux_max, ierr)
   end subroutine do_SED_Plot3

   subroutine do_SED_Plot4(s, id, device_id, xleft, xright, ybot, ytop, subplot, title, txt_scale, ierr)
      type(star_info), pointer :: s; integer :: id, device_id, ierr; real :: xleft, xright, ybot, ytop, txt_scale; logical :: subplot; character(len=*) :: title
      call render_sed_core(s, id, device_id, xleft, xright, ybot, ytop, subplot, title, txt_scale, &
         s% pg% SED4_lambda_min, s% pg% SED4_lambda_max, s% pg% SED4_flux_min, s% pg% SED4_flux_max, ierr)
   end subroutine do_SED_Plot4

   subroutine do_SED_Plot5(s, id, device_id, xleft, xright, ybot, ytop, subplot, title, txt_scale, ierr)
      type(star_info), pointer :: s; integer :: id, device_id, ierr; real :: xleft, xright, ybot, ytop, txt_scale; logical :: subplot; character(len=*) :: title
      call render_sed_core(s, id, device_id, xleft, xright, ybot, ytop, subplot, title, txt_scale, &
         s% pg% SED5_lambda_min, s% pg% SED5_lambda_max, s% pg% SED5_flux_min, s% pg% SED5_flux_max, ierr)
   end subroutine do_SED_Plot5

   subroutine do_SED_Plot6(s, id, device_id, xleft, xright, ybot, ytop, subplot, title, txt_scale, ierr)
      type(star_info), pointer :: s; integer :: id, device_id, ierr; real :: xleft, xright, ybot, ytop, txt_scale; logical :: subplot; character(len=*) :: title
      call render_sed_core(s, id, device_id, xleft, xright, ybot, ytop, subplot, title, txt_scale, &
         s% pg% SED6_lambda_min, s% pg% SED6_lambda_max, s% pg% SED6_flux_min, s% pg% SED6_flux_max, ierr)
   end subroutine do_SED_Plot6

   subroutine do_SED_Plot7(s, id, device_id, xleft, xright, ybot, ytop, subplot, title, txt_scale, ierr)
      type(star_info), pointer :: s; integer :: id, device_id, ierr; real :: xleft, xright, ybot, ytop, txt_scale; logical :: subplot; character(len=*) :: title
      call render_sed_core(s, id, device_id, xleft, xright, ybot, ytop, subplot, title, txt_scale, &
         s% pg% SED7_lambda_min, s% pg% SED7_lambda_max, s% pg% SED7_flux_min, s% pg% SED7_flux_max, ierr)
   end subroutine do_SED_Plot7

   subroutine do_SED_Plot8(s, id, device_id, xleft, xright, ybot, ytop, subplot, title, txt_scale, ierr)
      type(star_info), pointer :: s; integer :: id, device_id, ierr; real :: xleft, xright, ybot, ytop, txt_scale; logical :: subplot; character(len=*) :: title
      call render_sed_core(s, id, device_id, xleft, xright, ybot, ytop, subplot, title, txt_scale, &
         s% pg% SED8_lambda_min, s% pg% SED8_lambda_max, s% pg% SED8_flux_min, s% pg% SED8_flux_max, ierr)
   end subroutine do_SED_Plot8

   subroutine do_SED_Plot9(s, id, device_id, xleft, xright, ybot, ytop, subplot, title, txt_scale, ierr)
      type(star_info), pointer :: s; integer :: id, device_id, ierr; real :: xleft, xright, ybot, ytop, txt_scale; logical :: subplot; character(len=*) :: title
      call render_sed_core(s, id, device_id, xleft, xright, ybot, ytop, subplot, title, txt_scale, &
         s% pg% SED9_lambda_min, s% pg% SED9_lambda_max, s% pg% SED9_flux_min, s% pg% SED9_flux_max, ierr)
   end subroutine do_SED_Plot9

   subroutine render_sed_core(s, id, device_id, vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
                              w_min, w_max, f_min, f_max, ierr)
      use colors_def, only: Colors_General_Info
      use colors_lib, only: colors_ptr
      use pgstar_colors

      type(star_info), pointer :: s
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character(len=*), intent(in) :: title
      real, intent(in) :: w_min, w_max, f_min, f_max
      integer, intent(out) :: ierr
      
      type(Colors_General_Info), pointer :: cs
      real, allocatable :: x_wave(:), y_flux(:), filter_x(:), filter_y(:)
      real :: wave_min, wave_max, flux_min, flux_max
      integer :: i, n_lam, k
      character(len=strlen) :: final_title

      ierr = 0
      if (s% colors_handle <= 0) return

      call colors_ptr(s% colors_handle, cs, ierr)
      if (ierr /= 0 .or. .not. cs% pgstar_sed_valid) then
         call pgsvp(vp_xleft, vp_xright, vp_ybot, vp_ytop)
         call pgsch(txt_scale)
         call pgbox('BCNST', 0.0, 0, 'BCNST', 0.0, 0)
         call pgptxt(0.5, 0.5, 0.0, 0.5, 'Waiting for active Colors spectrum cache step...')
         return
      end if

      if (.not. allocated(cs% current_sed_wavelengths)) return
      n_lam = size(cs% current_sed_wavelengths)
      allocate(x_wave(n_lam), y_flux(n_lam), stat=ierr)
      if (ierr /= 0) return

      x_wave = real(cs% current_sed_wavelengths)
      
      where (cs% current_sed_fluxes > 0.0_dp)
         y_flux = real(log10(cs% current_sed_fluxes))
      elsewhere
         y_flux = -50.0 
      end where

      wave_min = merge(3000.0e0, w_min, w_min == -101.0e0)
      wave_max = merge(12000.0e0, w_max, w_max == -101.0e0)
      
      flux_min = f_min
      flux_max = f_max

      ! SAFE MASKING: Ignore the -50.0 padding bounds when auto-scaling limits
      if (flux_min == -101.0e0) then
         if (any(y_flux > -49.0e0)) then
            flux_min = minval(y_flux, mask=(y_flux > -49.0e0)) - 0.2e0
         else
            flux_min = -20.0e0
         end if
      end if

      if (flux_max == -101.0e0) then
         if (any(y_flux > -49.0e0)) then
            flux_max = maxval(y_flux, mask=(y_flux > -49.0e0)) + 0.2e0
         else
            flux_max = -10.0e0
         end if
      end if

      call pgsvp(vp_xleft, vp_xright, vp_ybot, vp_ytop)
      call pgswin(wave_min, wave_max, flux_min, flux_max)
      
      call pgsch(txt_scale)
      call pgsci(clr_Foreground)
      call pgbox('BCNST', 0.0, 0, 'BCNST', 0.0, 0)
      
      final_title = title
      if (len_trim(final_title) == 0) final_title = 'Synthetic Stellar SED Spectrum'
      
      call show_title_pgstar(s, final_title)
      call show_xaxis_label_pgstar(s, 'Wavelength (\A)')
      call show_left_yaxis_label_pgstar(s, 'log Flux (erg s\u-1\d cm\u-2\d \A\u-1\d)')

      ! Render base continuous spectrum spectrum
      call pgsci(clr_Foreground)
      call pgslw(s% pg% pgstar_lw)
      call pgline(n_lam, x_wave, y_flux)
      call pgslw(1)

      ! Overlay filter curves bounded nicely into the top 20% vertical region
      if (cs% filters_loaded) then
         do i = 1, size(cs% filters)
            if (.not. allocated(cs% filters(i)% wavelengths)) cycle
            
            k = size(cs% filters(i)% wavelengths)
            allocate(filter_x(k), filter_y(k))
            
            filter_x = real(cs% filters(i)% wavelengths)
            filter_y = flux_max - ((1.0e0 - real(cs% filters(i)% transmission)) * (flux_max - flux_min) * 0.20e0)
            
            call pgsci(clr_Cyan + mod(i, 4)) 
            call pgsls(2) 
            call pgline(k, filter_x, filter_y)
            
            call pgsch(txt_scale * 0.7e0)
            call pgptxt(filter_x(k/2), flux_max - 0.04e0 * (flux_max - flux_min), 0.0e0, 0.5e0, trim(cs% filters(i)% name))
            call pgsch(txt_scale)
            
            deallocate(filter_x, filter_y)
         end do
      end if

      call pgsls(1) 
      call pgsci(clr_Foreground)
      deallocate(x_wave, y_flux)

   end subroutine render_sed_core

end module pgstar_sed