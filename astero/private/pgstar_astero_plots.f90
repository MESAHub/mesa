! ***********************************************************************
!
!   Copyright (C) 2012  The MESA Team
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
 
      module pgstar_astero_plots
      use star_lib
      use star_def
      use astero_support

      implicit none
         

      contains
      
      
      subroutine astero_pgstar_plots_info(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         
         integer :: i, plot_id
         type (pgstar_win_file_data), pointer :: p
         type (star_info), pointer :: s

         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         i = 1
         plot_id = i_Other + i - 1
         p => s% pgstar_win_file_ptr(plot_id)
         p_echelle => p
         p% plot => echelle_plot
         p% okay_to_call_do_plot_in_grid = .true.
         p% do_plot_in_grid => do_echelle_plot_in_grid
         p% id = plot_id
         p% name = 'Echelle'
         p% win_flag = echelle_win_flag
         p% win_width = echelle_win_width
         p% win_aspect_ratio = echelle_win_aspect_ratio
         p% file_flag = echelle_file_flag
         p% file_dir = echelle_file_dir
         p% file_prefix = echelle_file_prefix
         p% file_interval = echelle_file_interval
         p% file_width = echelle_file_width
         p% file_aspect_ratio = echelle_file_aspect_ratio
         
         if (nl(1) > 0) then
            i = i+1
            plot_id = i_Other + i - 1
            p => s% pgstar_win_file_ptr(plot_id)
            p_ratios => p
            p% plot => ratios_plot
            p% okay_to_call_do_plot_in_grid = .true.
            p% do_plot_in_grid => do_ratios_plot_in_grid
            p% id = plot_id
            p% name = 'Ratios'
            p% win_flag = ratios_win_flag
            p% win_width = ratios_win_width
            p% win_aspect_ratio = ratios_win_aspect_ratio
            p% file_flag = ratios_file_flag
            p% file_dir = ratios_file_dir
            p% file_prefix = ratios_file_prefix
            p% file_interval = ratios_file_interval
            p% file_width = ratios_file_width
            p% file_aspect_ratio = ratios_file_aspect_ratio
         end if
         
      end subroutine astero_pgstar_plots_info
            

      subroutine echelle_plot(id, device_id, ierr)
         integer, intent(in) :: id, device_id
         integer, intent(out) :: ierr

         type (star_info), pointer :: s

         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         call pgslct(device_id)
         call pgbbuf()
         call pgeras()
         
         call do_echelle_plot(id, device_id, &
            echelle_xleft, echelle_xright, &
            echelle_ybot, echelle_ytop, &
            .false., echelle_title, echelle_txt_scale, ierr)

         call pgebuf()
      
      end subroutine echelle_plot


      subroutine do_echelle_plot_in_grid( &
            id, device_id, xleft, xright, ybot, ytop, txt_scale, ierr)         
         integer, intent(in) :: id, device_id
         real, intent(in) :: xleft, xright, ybot, ytop, txt_scale
         integer, intent(out) :: ierr
         call do_echelle_plot( &
            id, device_id, xleft, xright, ybot, ytop, .true., echelle_title, txt_scale, ierr)
      end subroutine do_echelle_plot_in_grid


      subroutine do_echelle_plot( &
            id, device_id, xleft, xright, ybot, ytop, subplot, title, txt_scale, ierr)
            
         use utils_lib
         use const_def
         
         integer, intent(in) :: id, device_id
         real, intent(in) :: xleft, xright, ybot, ytop, txt_scale
         logical, intent(in) :: subplot
         character (len=*), intent(in) :: title
         integer, intent(out) :: ierr

         type (star_info), pointer :: s
         real :: xmin, xmax, ymin, ymax, dx, dy, plot_delta_nu, freq, marker_scale, &
            x_obs, y_obs, x_model, y_model, y_txt, xpt_min, xpt_max, xmargin
         character (len=256) :: str
         integer :: i, &
            l0_color, l0_shape, l1_color, l1_shape, &
            l2_color, l2_shape, l3_color, l3_shape, &
            model_color, model_shape
                     
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         plot_delta_nu = echelle_delta_nu
         if (plot_delta_nu <= 0) plot_delta_nu = delta_nu
         if (plot_delta_nu <= 0) then
            write(*,*) 'must supply value for echelle_delta_nu'
            ierr = -1
            return
         end if
         
         xpt_min = 1e9
         xpt_max = -1e9
         ymin = 1e9
         ymax = 0         
         if (nl(0) > 0) then
            ymin = min(ymin,minval(freq_target(0,1:nl(0))))
            ymax = max(ymax,maxval(freq_target(0,1:nl(0))))
            xpt_min = min(xpt_min,minval(mod(freq_target(0,1:nl(0)),plot_delta_nu)))
            xpt_max = max(xpt_max,maxval(mod(freq_target(0,1:nl(0)),plot_delta_nu)))
         end if
         if (nl(1) > 0) then
            ymin = min(ymin,minval(freq_target(1,1:nl(1))))
            ymax = max(ymax,maxval(freq_target(1,1:nl(1))))
            xpt_min = min(xpt_min,minval(mod(freq_target(1,1:nl(1)),plot_delta_nu)))
            xpt_max = max(xpt_max,maxval(mod(freq_target(1,1:nl(1)),plot_delta_nu)))
         end if
         if (nl(2) > 0) then
            ymin = min(ymin,minval(freq_target(2,1:nl(2))))
            ymax = max(ymax,maxval(freq_target(2,1:nl(2))))
            xpt_min = min(xpt_min,minval(mod(freq_target(2,1:nl(2)),plot_delta_nu)))
            xpt_max = max(xpt_max,maxval(mod(freq_target(2,1:nl(2)),plot_delta_nu)))
         end if
         if (nl(3) > 0) then
            ymin = min(ymin,minval(freq_target(3,1:nl(3))))
            ymax = max(ymax,maxval(freq_target(3,1:nl(3))))
            xpt_min = min(xpt_min,minval(mod(freq_target(3,1:nl(3)),plot_delta_nu)))
            xpt_max = max(xpt_max,maxval(mod(freq_target(3,1:nl(3)),plot_delta_nu)))
         end if
         dy = ymax - ymin
         dy = max(dy, 1.0)
         ymin = ymin - dy*0.25
         ymax = ymax + dy*0.15

         xmargin = max(plot_delta_nu/5, (plot_delta_nu - (xpt_max - xpt_min))/2)
         xmin = -xmargin
         xmax = 2*plot_delta_nu + xmargin
         
         call pgsave
         
         call pgsch(txt_scale)
         call pgsvp(xleft, xright, ybot, ytop)
         call pgswin(xmin, xmax, ymin, ymax)
         call pgscf(1)
         call pgsci(1)
         call pgstar_show_box(s,'BCNST1','BCNSTV1')
         call pgstar_show_xaxis_label(s, &
            "Frequency mod \(0530)\d\(0639)\u (\(0638)Hz) (duplicated at x+\(0530)\d\(0639)\u)")
         call pgstar_show_left_yaxis_label(s,"Frequency (\(0638)Hz)")
         if (.not. subplot) then
            call pgstar_show_model_number(s)
            call pgstar_show_age(s)
         end if
         call pgstar_show_title(s, title)
         
         call pgslw(1)
         
         ! label
         y_obs = ymin + dy*0.12
         y_txt = ymin + dy*0.17
         if (nl(3) > 0) then
            dx = (xmax-xmin)/4d0
         else
            dx = (xmax-xmin)/3d0
         end if
         
         l0_color = clr_Teal
         l0_shape = 0840 ! circle
         
         l1_color = clr_Crimson
         l1_shape = 0842 ! triangle
         
         l2_color = clr_BrightBlue
         l2_shape = 0841 ! square
         
         l3_color = clr_Coral
         l3_shape = 0843 ! diamond
            
         model_color = clr_Silver
         model_shape = 0828 ! bullet
         
         
         x_obs = xmin + dx/2
         call pgsci(l0_color)
         call pgsch(1.6*txt_scale)
         call pgpt1(x_obs, y_obs, l0_shape)
         call pgsci(1)
         call pgsch(txt_scale)
         call pgptxt(x_obs, y_txt, 0.0, 0.5, 'l=0')
         
         x_obs = x_obs+dx
         call pgsci(l1_color)
         call pgsch(1.6*txt_scale)
         call pgpt1(x_obs, y_obs, l1_shape)
         call pgsci(1)
         call pgsch(1.0*txt_scale)
         call pgptxt(x_obs, y_txt, 0.0, 0.5, 'l=1')
         
         x_obs = x_obs+dx
         call pgsci(l2_color)
         call pgsch(1.6*txt_scale)
         call pgpt1(x_obs, y_obs, l2_shape)
         call pgsci(1)
         call pgsch(1.0*txt_scale)
         call pgptxt(x_obs, y_txt, 0.0, 0.5, 'l=2')
         
         if (nl(3) > 0) then
            x_obs = x_obs+dx
            call pgsci(l3_color)
            call pgsch(1.6*txt_scale)
            call pgpt1(x_obs, y_obs, l3_shape)
            call pgsci(1)
            call pgsch(1.0*txt_scale)
            call pgptxt(x_obs, y_txt, 0.0, 0.5, 'l=3')
         end if
         
         marker_scale = 2.4*txt_scale
         call pgsch(marker_scale)
         if (nl(0) > 0) then
            do i=1,nl(0)
               call show_obs(freq_target(0,i), l0_color, l0_shape)
               !if (have_radial) call show_model( &
               if (model_freq_corr(0,i) > 0) call show_model( &
                  freq_target(0,i), model_freq_corr(0,i), 0d0, 0d0, &
                  0d0, 0d0, 0d0, l0_color)
            end do
         end if
         
         if (nl(1) > 0) then
            do i=1,nl(1)
               call show_obs(freq_target(1,i), l1_color, l1_shape)
               if (have_nonradial) then
                  call show_model( &
                     freq_target(1,i), model_freq_corr(1,i), model_freq_corr_alt_up(1,i), &
                     model_freq_corr_alt_down(1,i), &
                     model_inertia(1,i), model_inertia_alt_up(1,i), model_inertia_alt_down(1,i), &
                     l1_color)
               end if
            end do
         end if
         
         if (nl(2) > 0) then
            do i=1,nl(2)
               call show_obs(freq_target(2,i), l2_color, l2_shape)
               if (have_nonradial) then
                  call show_model( &
                     freq_target(2,i), model_freq_corr(2,i), model_freq_corr_alt_up(2,i), &
                     model_freq_corr_alt_down(2,i), &
                     model_inertia(2,i), model_inertia_alt_up(2,i), model_inertia_alt_down(2,i), &
                     l2_color)
               end if
            end do
         end if
         
         if (nl(3) > 0) then
            do i=1,nl(3)
               call show_obs(freq_target(3,i), l3_color, l3_shape)
               if (have_nonradial) then
                  call show_model( &
                     freq_target(3,i), model_freq_corr(3,i), model_freq_corr_alt_up(3,i), &
                     model_freq_corr_alt_down(3,i), &
                     model_inertia(3,i), model_inertia_alt_up(3,i), model_inertia_alt_down(3,i), &
                     l3_color)
               end if
            end do
         end if
         
         call pgsci(clr_SlateGray)
         call pgsls(1)
         call pgslw(8)
         
         call pgmove(0., ymax - dy*0.08)
         call pgdraw(0., ymax)
         call pgmove(plot_delta_nu, ymax - dy*0.08)
         call pgdraw(plot_delta_nu, ymax)
         call pgmove(2*plot_delta_nu, ymax - dy*0.08)
         call pgdraw(2*plot_delta_nu, ymax)
         
         call pgmove(0., ymin + dy*0.08)
         call pgdraw(0., ymin)
         call pgmove(plot_delta_nu, ymin + dy*0.08)
         call pgdraw(plot_delta_nu, ymin)
         call pgmove(2*plot_delta_nu, ymin + dy*0.08)
         call pgdraw(2*plot_delta_nu, ymin)
         

         call pgunsa
      
         call show_pgstar_annotations(s, &
            show_echelle_annotation1, &
            show_echelle_annotation2, &
            show_echelle_annotation3)


         contains
         
         
         subroutine show_obs(freq, color, shape)
            real(dp), intent(in) :: freq
            integer, intent(in) :: color, shape
            y_obs = freq
            x_obs = mod(freq,plot_delta_nu)               
            call pgsci(color)
            call pgpt1(x_obs, y_obs, shape)
            call pgpt1(x_obs + plot_delta_nu, y_obs, shape)
         end subroutine show_obs
         
         subroutine show_model( &
               freq_obs, freq, freq_alt_up, freq_alt_down, &
               inertia, inertia_alt_up, inertia_alt_down, color)
            real(dp), intent(in) :: freq_obs, freq, freq_alt_up, freq_alt_down, &
               inertia, inertia_alt_up, inertia_alt_down
            integer, intent(in) :: color
            real :: x_obs, x_model_alt_up, x_model_alt_down
            real :: y_model_alt_up, y_model_alt_down
            real :: y_model_alt_shift
            include 'formats'
            y_model_alt_shift = echelle_model_alt_y_shift
            call pgsci(color)
            y_model = freq
            x_obs = mod(freq_obs, plot_delta_nu)
            x_model = (freq - freq_obs) + x_obs
            call pgpt1(x_model, y_model, model_shape)
            call pgmove(x_obs, y_obs)
            call pgdraw(x_model, y_model)
            call pgpt1(x_model + plot_delta_nu, y_model, model_shape)
            call pgmove(x_obs + plot_delta_nu, y_obs)
            call pgdraw(x_model + plot_delta_nu, y_model)            
            if (freq_alt_up > 0d0 .and. show_echelle_next_best_at_higher_frequency) then
               y_model_alt_up = freq_alt_up + y_model_alt_shift
               x_model_alt_up = (freq_alt_up - freq_obs) + x_obs
               call pgsch(marker_scale*real(sqrt(inertia/inertia_alt_up)))
               call pgpt1(x_model_alt_up, y_model_alt_up, model_shape)
               call pgpt1(x_model_alt_up + plot_delta_nu, y_model_alt_up, model_shape)
               call pgsch(marker_scale)
            end if            
            if (freq_alt_down > 0d0 .and. show_echelle_next_best_at_lower_frequency) then
               y_model_alt_down = freq_alt_down - y_model_alt_shift
               x_model_alt_down = (freq_alt_down - freq_obs) + x_obs
               call pgsch(marker_scale*real(sqrt(inertia/inertia_alt_down)))
               call pgpt1(x_model_alt_down, y_model_alt_down, model_shape)
               call pgpt1(x_model_alt_down + plot_delta_nu, y_model_alt_down, model_shape)
               call pgsch(marker_scale)
            end if           
         end subroutine show_model
         
      end subroutine do_echelle_plot
      

      subroutine ratios_plot(id, device_id, ierr)
         integer, intent(in) :: id, device_id
         integer, intent(out) :: ierr

         type (star_info), pointer :: s

         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         call pgslct(device_id)
         call pgbbuf()
         call pgeras()
         
         call do_ratios_plot(id, device_id, &
            ratios_xleft, ratios_xright, &
            ratios_ybot, ratios_ytop, &
            .false., ratios_title, ratios_txt_scale, ierr)

         call pgebuf()
      
      end subroutine ratios_plot


      subroutine do_ratios_plot_in_grid( &
            id, device_id, xleft, xright, ybot, ytop, txt_scale, ierr)         
         integer, intent(in) :: id, device_id
         real, intent(in) :: xleft, xright, ybot, ytop, txt_scale
         integer, intent(out) :: ierr
         call do_ratios_plot( &
            id, device_id, xleft, xright, ybot, ytop, .true., ratios_title, txt_scale, ierr)
      end subroutine do_ratios_plot_in_grid


      subroutine do_ratios_plot( &
            id, device_id, xleft, xright, ybot, ytop, subplot, title, txt_scale, ierr)
            
         use utils_lib
         use const_def
         
         integer, intent(in) :: id, device_id
         real, intent(in) :: xleft, xright, ybot, ytop, txt_scale
         logical, intent(in) :: subplot
         character (len=*), intent(in) :: title
         integer, intent(out) :: ierr

         type (star_info), pointer :: s
         real :: xmin, xmax, ymin, ymax, dx, dy, freq, &
            x_obs, y_obs, x_model, y_model, y_txt, sig_max
         character (len=256) :: str
         logical :: show_model
         integer :: i, n, i0, i1, l0_first, l1_first, &
            r01_color, r01_shape, r10_color, r10_shape, &
            r02_color, r02_shape, model_color, model_shape
                     
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         if (chi2_seismo_r_010_fraction <= 0d0 .and. &
             chi2_seismo_r_02_fraction <= 0d0) then
            return
         end if
         
         if (nl(1) <= 0 .or. ratios_n <= 0) then
            return
         end if
         
         n = ratios_n
         l0_first = ratios_l0_first
         l1_first = ratios_l1_first

         xmax = -HUGE(xmax)
         xmin = HUGE(xmin)
         do i=1,n
            i0 = i + l0_first
            i1 = i + l1_first            
            if (ratios_r01(i) > xmax) xmax = ratios_r01(i)
            if (ratios_r01(i) < xmin) xmin = ratios_r01(i)
            if (ratios_r10(i) > xmax) xmax = ratios_r10(i)
            if (ratios_r10(i) < xmin) xmin = ratios_r10(i)
         end do
         do i=1,nl(0)
            if (sigmas_r02(i) == 0) cycle
            if (ratios_r02(i) > xmax) xmax = ratios_r02(i)
            if (ratios_r02(i) < xmin) xmin = ratios_r02(i)
         end do
         sig_max = max(maxval( &
            sigmas_r01(1:n)), maxval(sigmas_r10(1:n)), maxval(sigmas_r02(1:n)))
         xmin = xmin - ratios_margin_sig_factor*sig_max
         xmax = xmax + ratios_margin_sig_factor*sig_max
         dx = xmax - xmin
         dx = max(dx, 0.02)
         xmin = xmin - dx*0.1
         xmax = xmax + dx*0.1
         
         ymin = freq_target(0,1 + l0_first)
         ymax = freq_target(1,n + l1_first)
         do i=2,nl(0)
            if (sigmas_r02(i) == 0d0) cycle
            if (freq_target(0,i) > ymax) ymax = freq_target(0,i)
            if (freq_target(0,i) < ymin) ymin = freq_target(0,i)
         end do
         dy = ymax - ymin
         dy = max(dy, 1.0)
         ymin = ymin - dy*0.25
         ymax = ymax + dy*0.12
         
         call pgsave

         call pgsvp(xleft, xright, ybot, ytop)
         call pgswin(xmin, xmax, ymin, ymax)
         call pgscf(1)
         call pgsci(1)
         call pgstar_show_box(s,'BCNST1','BCNSTV1')
         call pgstar_show_xaxis_label(s,"Ratio")
         call pgstar_show_left_yaxis_label(s,"Frequency (\(0638)Hz)")
         if (.not. subplot) then
            call pgstar_show_model_number(s)
            call pgstar_show_age(s)
         end if
         call pgstar_show_title(s, title)
         
         call pgslw(1)
         
         r01_color = clr_Teal
         r01_shape = 0840 ! circle
         
         r10_color = clr_Crimson
         r10_shape = 0842 ! triangle
         
         r02_color = clr_BrightBlue
         r02_shape = 0841 ! square
         
         model_color = clr_Silver
         model_shape = 0828 ! bullet
         
         ! label
         y_obs = ymin + dy*0.06
         y_txt = ymin + dy*0.10
         dx = (xmax-xmin)/4d0
         
         x_obs = xmin+dx
         call pgsci(r01_color)
         call pgsch(1.6*txt_scale)
         call pgpt1(x_obs, y_obs, r01_shape)
         call pgsci(1)
         call pgsch(1.0*txt_scale)
         call pgptxt(x_obs, y_txt, 0.0, 0.5, 'r01')
         
         x_obs = x_obs+dx
         call pgsci(r10_color)
         call pgsch(1.6*txt_scale)
         call pgpt1(x_obs, y_obs, r10_shape)
         call pgsci(1)
         call pgsch(1.0*txt_scale)
         call pgptxt(x_obs, y_txt, 0.0, 0.5, 'r10')
         
         x_obs = x_obs+dx
         call pgsci(r02_color)
         call pgsch(1.6*txt_scale)
         call pgpt1(x_obs, y_obs, r02_shape)
         call pgsci(1)
         call pgsch(1.0*txt_scale)
         call pgptxt(x_obs, y_txt, 0.0, 0.5, 'r02')
                  
         show_model = &
            (model_ratios_n == ratios_n .and. &
               model_ratios_l0_first == ratios_l0_first .and. &
                  model_ratios_l1_first == ratios_l1_first)
         
         call pgsch(2.4*txt_scale)
         do i=1,n
            call show_r01(i)
            call show_r10(i)
         end do
         
         do i=1,nl(0)
            call show_r02(i)
         end do
         
         call pgunsa
      
         call show_pgstar_annotations(s, &
            show_ratios_annotation1, &
            show_ratios_annotation2, &
            show_ratios_annotation3)

         contains
         
         subroutine show_r01(i)
            integer, intent(in) :: i
            real :: x_obs, y_obs, sig_obs, x_model, y_model
            include 'formats'
            y_obs = freq_target(0,i + l0_first)
            x_obs = ratios_r01(i)
            sig_obs = sigmas_r01(i)
            call pgsci(r01_color)
            call pgmove(x_obs-sig_obs, y_obs)
            call pgdraw(x_obs+sig_obs, y_obs)
            if (show_model) then
               y_model = model_freq_corr(0,i + l0_first)
               x_model = model_ratios_r01(i)
               call pgmove(x_obs, y_obs)
               call pgdraw(x_model, y_model)
               call pgpt1(x_model, y_model, model_shape)
            end if
            call pgpt1(x_obs, y_obs, r01_shape)
         end subroutine show_r01
         
         subroutine show_r10(i)
            integer, intent(in) :: i
            real :: x_obs, y_obs, sig_obs, x_model, y_model
            include 'formats'
            y_obs = freq_target(1,i + l1_first)
            x_obs = ratios_r10(i)
            sig_obs = sigmas_r10(i)
            call pgsci(r10_color)
            call pgmove(x_obs-sig_obs, y_obs)
            call pgdraw(x_obs+sig_obs, y_obs)
            if (show_model) then
               y_model = model_freq_corr(1,i + l1_first)
               x_model = model_ratios_r10(i)
               call pgmove(x_obs, y_obs)
               call pgdraw(x_model, y_model)
               call pgpt1(x_model, y_model, model_shape)
            end if
            call pgpt1(x_obs, y_obs, r10_shape)
         end subroutine show_r10
         
         subroutine show_r02(i)
            integer, intent(in) :: i
            real :: x_obs, y_obs, sig_obs, x_model, y_model
            include 'formats'
            y_obs = freq_target(0,i + l0_first)
            x_obs = ratios_r02(i)
            sig_obs = sigmas_r02(i)
            if (sig_obs == 0d0) return
            call pgsci(r02_color)
            call pgmove(x_obs-sig_obs, y_obs)
            call pgdraw(x_obs+sig_obs, y_obs)
            if (show_model) then
               y_model = model_freq_corr(0,i + l0_first)
               x_model = model_ratios_r02(i)
               call pgmove(x_obs, y_obs)
               call pgdraw(x_model, y_model)
               call pgpt1(x_model, y_model, model_shape)
            end if
            call pgpt1(x_obs, y_obs, r02_shape)
         end subroutine show_r02
         
      end subroutine do_ratios_plot
      
      
      subroutine write_plot_to_file(s, p, file_prefix, number, ierr)
         use star_lib, only: pgstar_write_plot_to_file
         type (star_info), pointer :: s
         type (pgstar_win_file_data), pointer :: p
         character (len=*), intent(in) :: file_prefix
         integer, intent(in) :: number
         integer, intent(out) :: ierr

         character (len=256) :: format_string, num_str, name, extension
         integer :: len
         
         ierr = 0
         
         if (len_trim(file_prefix) == 0 .or. .not. associated(p)) return
         
         write(format_string, '( "(i",i2.2,".",i2.2,")" )') num_digits, num_digits
         write(num_str, format_string) number
         
         if (len_trim(p% file_dir) > 0) then
            name = trim(p% file_dir) // '/' // trim(file_prefix)
         else
            name = file_prefix
         end if
         
         extension = 'png' ! s% file_extension
         name = trim(name) // '_sample' // trim(num_str) // '.' // trim(extension)
         
         write(*,'(a)') 'write plot to file ' // trim(name)
         call pgstar_write_plot_to_file(s, p, name, ierr)
      
      end subroutine write_plot_to_file
      

      end module pgstar_astero_plots
      
      
      
      
