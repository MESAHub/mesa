! ***********************************************************************
!
!   Copyright (C) 2015-2019  The MESA Team
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

module pgstar_rti
   
   use star_private_def
   use const_def
   use pgstar_support
   use star_pgstar
   
   implicit none


contains
   
   
   subroutine rti_Plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      
      ierr = 0
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) return
      
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      
      call do_rti_Plot(s, id, device_id, &
         s% pg% rti_xleft, s% pg% rti_xright, &
         s% pg% rti_ybot, s% pg% rti_ytop, .false., &
         s% pg% rti_title, s% pg% rti_txt_scale, ierr)
      if (ierr /= 0) return
      
      call pgebuf()
   
   end subroutine rti_Plot
   
   
   subroutine do_rti_Plot(s, id, device_id, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      use chem_def
      use net_def
      use utils_lib
      
      type (star_info), pointer :: s
      integer, intent(in) :: id, device_id
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      
      integer :: i, ii, n, step_min, step_max
      real :: xmin, xmax, ymin_L_axis, ymax_L_axis, &
         ymin_mass_axis, ymax_mass_axis, dx
      real, allocatable, dimension(:) :: xvec, star_mass, star_M_center, log_xmstar, &
         he_core_mass, &
         c_core_mass, &
         o_core_mass, &
         si_core_mass, &
         fe_core_mass, &
         shock_mass
      logical :: have_star_mass, have_log_xmstar, &
         have_he_core_mass, &
         have_c_core_mass, &
         have_o_core_mass, &
         have_si_core_mass, &
         have_fe_core_mass, &
         have_shock_mass
      
      integer :: ix, k
      real :: xleft, xright, now
      real :: dxmin = -1.d0
      
      include 'formats'
      
      ierr = 0
      
      step_min = s% pg% rti_step_xmin
      if (step_min <= 0) step_min = 1
      step_max = s% pg% rti_step_xmax
      if (step_max <= 0 .or. step_max > s% model_number) step_max = s% model_number
      
      if (step_min >= s% model_number) step_min = 1
      
      if (s% pg% rti_max_width > 0) &
         step_min = max(step_min, step_max - s% pg% rti_max_width)
      
      n = count_hist_points(s, step_min, step_max)
      if (n <= 1) return
      step_min = max(step_min, step_max - n + 1)
      
      call integer_dict_lookup(s% history_names_dict, s% pg% rti_xaxis_name, ix, ierr)
      if (ierr /= 0) ix = -1
      if (ix <= 0) then
         write(*, '(A)')
         write(*, *) 'ERROR: failed to find ' // &
            trim(s% pg% rti_xaxis_name) // ' in rti data'
         write(*, '(A)')
         ierr = -1
      end if
      
      allocate(xvec(n), &
         star_mass(n), &
         log_xmstar(n), &
         star_M_center(n), &
         he_core_mass(n), &
         c_core_mass(n), &
         o_core_mass(n), &
         si_core_mass(n), &
         fe_core_mass(n), &
         shock_mass(n), &
         stat = ierr)
      if (ierr /= 0) then
         write(*, *) 'allocate failed for PGSTAR rti'
         return
      end if
      
      call get_hist_points(s, step_min, step_max, n, ix, xvec, ierr)
      if (ierr /= 0) then
         write(*, *) 'pgstar get_hist_points failed ' // trim(s% pg% rti_xaxis_name)
         call dealloc
         ierr = 0
         return
      end if
      
      if (s% pg% rti_xaxis_in_seconds .and. s% pg% rti_xaxis_name=='star_age') then
         do k = 1, n
            xvec(k) = xvec(k) * secyer
         end do
      else if (s% pg% rti_xaxis_in_Myr .and. s% pg% rti_xaxis_name=='star_age') then
         do k = 1, n
            xvec(k) = xvec(k) * 1d-6
         end do
      end if
      
      now = xvec(n)
      if (s% pg% rti_xaxis_time_from_present .and. s% pg% rti_xaxis_name=='star_age') then
         do k = 1, n
            xvec(k) = xvec(k) - now
         end do
      end if
      
      if (s% pg% rti_xaxis_log) then
         do k = 1, n
            xvec(k) = log10(max(1e-50, abs(xvec(k))))
         end do
      end if
      
      if(s% pg% rti_xmin<-100d0) s% pg% rti_xmin = xvec(1)
      if(s% pg% rti_xmax<-100d0) s% pg% rti_xmax = xvec(n)
      
      xmin = max(s% pg% rti_xmin, xvec(1))
      xmax = min(s% pg% rti_xmax, xvec(n))
      
      call set_xleft_xright(&
         n, xvec, xmin, xmax, s% pg% rti_xmargin, &
         s% pg% rti_xaxis_reversed, dxmin, xleft, xright)
      
      have_star_mass = get1_yvec('star_mass', star_mass)
      if (.not. have_star_mass) then
         write(*, *) 'failed to find star_mass in history data'
         ierr = -1
      end if
      if (ierr /= 0) return
      
      have_log_xmstar = get1_yvec('log_xmstar', log_xmstar)
      if (have_log_xmstar) then
         do i = 1, n
            star_M_center(i) = &
               star_mass(i) - real(exp10(dble(log_xmstar(i))) / Msun)
         end do
      else
         star_M_center(:) = 0
      end if
      
      have_he_core_mass = get1_yvec('he_core_mass', he_core_mass)
      have_c_core_mass = get1_yvec('c_core_mass', c_core_mass)
      have_o_core_mass = get1_yvec('o_core_mass', o_core_mass)
      have_si_core_mass = get1_yvec('si_core_mass', si_core_mass)
      have_fe_core_mass = get1_yvec('fe_core_mass', fe_core_mass)
      have_shock_mass = get1_yvec('shock_mass', shock_mass)
      
      call pgsave
      call pgsch(txt_scale)
      
      dx = (xmax - xmin) / 250.0
      
      call init_rti_plot
      call setup_mass_yaxis
      call plot_total_mass_line
      call plot_rti_data
      call plot_mass_lines
      
      call show_annotations(s, &
         s% pg% show_rti_annotation1, &
         s% pg% show_rti_annotation2, &
         s% pg% show_rti_annotation3)
      
      call pgsch(txt_scale)
      call finish_rti_plot
      
      call pgunsa
      
      call dealloc
   
   
   contains
      
      subroutine dealloc
         deallocate(xvec, star_mass, star_M_center, log_xmstar, &
            he_core_mass, &
            c_core_mass, &
            o_core_mass, &
            si_core_mass, &
            fe_core_mass, &
            shock_mass)
      end subroutine dealloc
      
      subroutine plot_total_mass_line
         call pgsave
         call pgsch(txt_scale * 0.8)
         call pgsci(clr_Gray)
         call pgslw(2)
         call show_left_yaxis_label_pgmtxt_pgstar(s, 1.0, 1.0, 'M\dtotal\u', -0.8)
         call pgslw(s% pg% pgstar_lw)
         call plot_mass_line(star_mass)
         call pgunsa
      end subroutine plot_total_mass_line
      
      
      subroutine plot_mass_line(yvec)
         real, intent(in) :: yvec(:)
         if (any(yvec > 1e-2 * ymax_mass_axis)) then
            call pgsave
            call pgslw(s% pg% rti_line_weight)
            call pgline(n, xvec, yvec)
            call pgunsa
         end if
      end subroutine plot_mass_line
      
      
      subroutine setup_mass_yaxis
         real :: dy, ymin, ymax
         include 'formats'
         ymax = s% pg% rti_mass_max
         if (ymax <= 0) ymax = maxval(star_mass)
         ymin = s% pg% rti_mass_min
         if (ymin < 0) ymin = minval(star_M_center)
         if (ymax <= ymin) ymax = ymin + 1
         dy = ymax - ymin
         if (s% pg% rti_mass_min /= -101d0) ymin = ymin - s% pg% rti_mass_margin * dy
         if (s% pg% rti_mass_max /= -101d0) ymax = ymax + s% pg% rti_mass_margin * dy
         call pgswin(xleft, xright, ymin, ymax)
         call pgscf(1)
         call pgsci(1)
         ymin_mass_axis = ymin
         ymax_mass_axis = ymax
      end subroutine setup_mass_yaxis
      
      
      subroutine init_rti_plot
         call pgsvp(vp_xleft, vp_xright, vp_ybot, vp_ytop)
      end subroutine init_rti_plot
      
      
      subroutine finish_rti_plot
         character(len = 256) :: xlabel
         call show_box_pgstar(s, 'BCNST1', 'BCNMSTV1')
         
         xlabel = ''
         if (s% pg% rti_xaxis_log) then
            xlabel = 'log ' // s% pg% rti_xaxis_name
         else
            xlabel = s% pg% rti_xaxis_name
         end if
         
         if (s% pg% rti_xaxis_name =='star_age') then
            if (s% pg% rti_xaxis_in_seconds) then
               xlabel = trim(xlabel) // ' (s)'
            else if (s% pg% rti_xaxis_in_Myr) then
               xlabel = trim(xlabel) // ' (Myr)'
            else
               xlabel = trim(xlabel) // ' (yr)'
            end if
         end if
         
         call show_xaxis_label_pgstar(s, trim(xlabel), 1.0)
         
         call show_left_yaxis_label_pgstar(s, 'M/M\d\(2281)')
         if (.not. subplot) then
            call show_model_number_pgstar(s)
            call show_age_pgstar(s)
         end if
         call show_title_pgstar(s, title)
         
         call show_pgstar_decorator(s% pg%id, s% pg% rti_use_decorator, &
            s% pg% rti_pgstar_decorator, 0, ierr)
      
      end subroutine finish_rti_plot
      
      
      logical function get1_yvec(name, vec)
         character (len = *) :: name
         real, dimension(:), allocatable :: vec
         get1_yvec = get1_hist_yvec(s, step_min, step_max, n, name, vec)
      end function get1_yvec
      
      
      subroutine plot_rti_data
         use history_specs, only : rti_offset
         type (pgstar_hist_node), pointer :: pg
         integer :: i_rti_type_first, i_rti_type_last
         integer :: k, cnt, num_specs, step
         
         include 'formats'
         
         i_rti_type_first = 0
         num_specs = size(s% history_column_spec, dim = 1)
         do k = 1, num_specs
            if (s% history_column_spec(k) == rti_offset + 1) then
               i_rti_type_first = k
               exit
            end if
         end do
         if (i_rti_type_first == 0) then
            write(*, *) 'i_rti_type_first == 0'
            return
         end if
         
         i_rti_type_last = 0
         cnt = 1
         do k = i_rti_type_first + 1, num_specs
            i_rti_type_last = k - 1
            cnt = cnt + 1
            if (s% history_column_spec(k) /= rti_offset + cnt) exit
         end do
         
         call pgsave
         call pgslw(s% pg% rti_line_weight)
         if (.not. associated(s% pg% pgstar_hist)) then
            write(*, *) '.not. associated(s% pg% pgstar_hist)'
         end if
         pg => s% pg% pgstar_hist
         
         do
            if (.not. associated(pg)) exit
            step = pg% step
            if (step < step_min) exit
            if (step <= step_max .and. mod(step, s% pg% rti_interval) == 0) then
               call draw_rti_for_step(&
                  pg, step, i_rti_type_first, i_rti_type_last, real(xvec(step - step_min + 1)), &
                  star_mass(step - step_min + 1), star_M_center(step - step_min + 1))
            end if
            pg => pg% next
         end do
         call pgunsa
      
      end subroutine plot_rti_data
      
      
      subroutine draw_rti_for_step(&
         pg, step, i_rti_type_first, i_rti_type_last, xval, mass, mass_center)
         type (pgstar_hist_node), pointer :: pg
         integer, intent(in) :: step, i_rti_type_first, i_rti_type_last
         real, intent(in) :: xval, mass, mass_center
         real :: qbot, qtop, mbot, mtop, xmass
         integer :: k, rti_type
         include 'formats'
         qbot = 0
         xmass = mass - mass_center
         do k = i_rti_type_first, i_rti_type_last, 2
            rti_type = int(pg% vals(k))
            if (rti_type < 0) exit
            qtop = pg% vals(k + 1)
            if (rti_type > 0) then
               mbot = mass_center + xmass * qbot
               mtop = mass_center + xmass * qtop
               call draw1(xval, mbot, mtop, clr_Blue)
            end if
            qbot = qtop
         end do
      end subroutine draw_rti_for_step
      
      
      subroutine draw1(xval, y1, y2, clr)
         real, intent(in) :: xval, y1, y2
         integer, intent(in) :: clr
         real :: top, bot
         if (y1 < y2) then
            bot = y1; top = y2
         else
            bot = y2; top = y1
         end if
         if (top < 1d-50) return
         call pgsci(clr)
         call pgmove(xval, bot)
         call pgdraw(xval, top)
      end subroutine draw1
      
      
      subroutine plot_mass_lines
         integer :: i
         
         include 'formats'
         
         call pgsave
         call pgsch(txt_scale * 0.8)
         call pgslw(2)
         
         if (have_he_core_mass) then
            call pgsci(clr_Teal)
            call show_left_yaxis_label_pgmtxt_pgstar(s, 0.77, 0.5, 'He', -0.8)
            call plot_mass_line(he_core_mass)
         else
            write(*, *) "please add 'he_core_mass' to your history_columns.list"
         end if
         
         if (have_c_core_mass) then
            call pgsci(clr_LightOliveGreen)
            call show_left_yaxis_label_pgmtxt_pgstar(s, 0.59, 0.5, 'C', -0.8)
            call plot_mass_line(c_core_mass)
         else
            write(*, *) "please add 'c_core_mass' to your history_columns.list"
         end if
         
         if (have_o_core_mass) then
            call pgsci(clr_SeaGreen)
            call show_left_yaxis_label_pgmtxt_pgstar(s, 0.41, 0.5, 'O', -0.8)
            call plot_mass_line(o_core_mass)
         else
            write(*, *) "please add 'o_core_mass' to your history_columns.list"
         end if
         
         if (have_si_core_mass) then
            call pgsci(clr_Lilac)
            call show_left_yaxis_label_pgmtxt_pgstar(s, 0.23, 0.5, 'Si', -0.8)
            call plot_mass_line(si_core_mass)
         else
            write(*, *) "please add 'si_core_mass' to your history_columns.list"
         end if
         
         if (have_fe_core_mass) then
            call pgsci(clr_Crimson)
            call show_left_yaxis_label_pgmtxt_pgstar(s, 0.00, 0.0, 'Iron', -0.8)
            call plot_mass_line(fe_core_mass)
         else
            write(*, *) "please add 'fe_core_mass' to your history_columns.list"
         end if
         
         if (have_shock_mass) then
            call pgsci(clr_IndianRed)
            call show_right_yaxis_label_pgmtxt_pgstar(s, 0.5, 0.5, 'Outward Shock', -0.8)
            call plot_mass_line(shock_mass)
         else
            write(*, *) "please add 'shock_mass' to your history_columns.list"
         end if
         
         call pgunsa
      
      end subroutine plot_mass_lines
   
   
   end subroutine do_rti_Plot


end module pgstar_rti

