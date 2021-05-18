! ***********************************************************************
!
!   Copyright (C) 2010-2019  THe MESA Team
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

      module pgstar_kipp

      use star_private_def
      use const_def
      use pgstar_support

      implicit none


      contains


      subroutine Kipp_Plot(id, device_id, ierr)
         integer, intent(in) :: id, device_id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         call pgslct(device_id)
         call pgbbuf()
         call pgeras()

         call do_Kipp_Plot(s, id, device_id, &
            s% Kipp_xleft, s% Kipp_xright, &
            s% Kipp_ybot, s% Kipp_ytop, .false., &
            s% Kipp_title, s% Kipp_txt_scale, ierr)
         if (ierr /= 0) return

         call pgebuf()

      end subroutine Kipp_Plot


      subroutine do_Kipp_Plot(s, id, device_id, &
            vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
         use chem_def
         use net_def
         use utils_lib

         type (star_info), pointer :: s
         integer, intent(in) :: id, device_id
         real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
         logical, intent(in) :: subplot
         character (len=*), intent(in) :: title
         integer, intent(out) :: ierr

         integer :: i, ii, n, step_min, step_max
         real :: xmin, xmax, ymin_L_axis, ymax_L_axis, &
            ymin_mass_axis, ymax_mass_axis, dx, burn_type_cutoff
         real, allocatable, dimension(:) :: xvec, &
            log_L, &
            log_Lneu, &
            log_LH, &
            log_LHe, &
            star_mass, &
            log_xmstar, &
            star_M_center, &
            he_core_mass, &
            c_core_mass, &
            o_core_mass, &
            si_core_mass, &
            fe_core_mass
         logical :: &
            have_log_L, &
            have_log_Lneu, &
            have_log_LH, &
            have_log_LHe, &
            have_star_mass, &
            have_log_xmstar, &
            have_he_core_mass, &
            have_c_core_mass, &
            have_o_core_mass, &
            have_si_core_mass, &
            have_fe_core_mass
         integer, parameter :: max_mix_type = leftover_convective_mixing
         integer :: mix_clr(max_mix_type), max_mix_type_to_show
         logical :: showed_this_mix_type(max_mix_type)

         integer :: ix,k
         real :: xleft,xright,now
         real :: dxmin=-1.d0

         include 'formats'

         ierr = 0

         mix_clr(convective_mixing) = clr_convection
         mix_clr(overshoot_mixing) = clr_overshoot
         mix_clr(semiconvective_mixing) = clr_semiconvection
         mix_clr(thermohaline_mixing) = clr_thermohaline
         mix_clr(rotation_mixing) = clr_rotation
         mix_clr(rayleigh_taylor_mixing) = clr_rayleigh_taylor
         mix_clr(minimum_mixing) = clr_minimum
         mix_clr(anonymous_mixing) = clr_anonymous
         mix_clr(leftover_convective_mixing) = clr_leftover_convection

         max_mix_type_to_show = thermohaline_mixing
         showed_this_mix_type = .false.

         step_min = s% Kipp_step_xmin
         if (step_min <= 0) step_min = 1
         step_max = s% Kipp_step_xmax
         if (step_max <= 0 .or. step_max > s% model_number) step_max = s% model_number

         if (step_min >= s% model_number) step_min = 1

         if (s% Kipp_max_width > 0) &
            step_min = max(step_min, step_max - s% Kipp_max_width)

         n = count_hist_points(s, step_min, step_max)
         if (n <= 1) return
         step_min = max(step_min, step_max-n+1)

         call integer_dict_lookup(s% history_names_dict, s% kipp_xaxis_name, ix, ierr)
         if (ierr /= 0) ix = -1
         if (ix <= 0) then
            write(*,*)
            write(*,*) 'ERROR: failed to find ' // &
               trim(s% kipp_xaxis_name) // ' in kipp data'
            write(*,*)
            ierr = -1
         end if

         allocate(xvec(n), &
            log_L(n), &
            log_Lneu(n), &
            log_LH(n), &
            log_LHe(n), &
            star_mass(n), &
            log_xmstar(n), &
            star_M_center(n), &
            he_core_mass(n), &
            c_core_mass(n), &
            o_core_mass(n), &
            si_core_mass(n), &
            fe_core_mass(n), &
            stat=ierr)
         if (ierr /= 0) then
            write(*,*) 'allocate failed for PGSTAR Kipp'
            return
         end if

         call get_hist_points(s, step_min, step_max, n, ix, xvec, ierr)
         if (ierr /= 0) then
            write(*,*) 'pgstar get_hist_points failed ' // trim(s% kipp_xaxis_name)
            call dealloc
            ierr = 0
            return
         end if

         if (s% kipp_xaxis_in_seconds .and. s% kipp_xaxis_name=='star_age')THEN
            do k=1,n
               xvec(k) = xvec(k)*secyer
            end do
         else if (s% kipp_xaxis_in_Myr .and. s% kipp_xaxis_name=='star_age')THEN
            do k=1,n
               xvec(k) = xvec(k)*1d-6
            end do
         end if

         now=xvec(n)
         if (s% kipp_xaxis_time_from_present .and. s% kipp_xaxis_name=='star_age') then
            do k=1,n
               xvec(k) = xvec(k)-now
            end do
         end if

         if (s% kipp_xaxis_log) then
            do k=1,n
               xvec(k) = log10(max(tiny(xvec(k)),abs(xvec(k))))
            end do
         end if

         if(s% kipp_xmin<-100d0) s% kipp_xmin=xvec(1)
         if(s% kipp_xmax<-100d0) s% kipp_xmax=xvec(n)

         xmin=max(s% kipp_xmin,xvec(1))
         xmax=min(s% kipp_xmax,xvec(n))
         
         burn_type_cutoff = s% Kipp_burn_type_cutoff

         call set_xleft_xright( &
            n, xvec, xmin, xmax, s% kipp_xmargin, &
            s% kipp_xaxis_reversed, dxmin, xleft, xright)

         have_star_mass = get1_yvec('star_mass', star_mass)
         if (.not. have_star_mass) then
            write(*,*) 'PGSTAR Kipp failed to find star_mass in history data'
            ierr = -1
         end if
         if (ierr /= 0) return
         have_log_xmstar = get1_yvec('log_xmstar', log_xmstar)
         if (have_log_xmstar) then
            do i = 1, n
               star_M_center(i) = &
                  star_mass(i) - real(exp10(dble(log_xmstar(i)))/Msun)
            end do
         else
            star_M_center(:) = 0
         end if

         if (s% Kipp_show_luminosities) then
            have_log_L = get1_yvec('log_L', log_L)
            have_log_Lneu = get1_yvec('log_Lneu', log_Lneu)
            have_log_LH = get1_yvec('log_LH', log_LH)
            have_log_LHe = get1_yvec('log_LHe', log_LHe)
         else
            have_log_L = .false.
            have_log_Lneu = .false.
            have_log_LH = .false.
            have_log_LHe = .false.
         end if

         if (s% Kipp_show_mass_boundaries) then
            have_he_core_mass = get1_yvec('he_core_mass', he_core_mass)
            have_c_core_mass = get1_yvec('c_core_mass', c_core_mass)
            have_o_core_mass = get1_yvec('o_core_mass', o_core_mass)
            have_si_core_mass = get1_yvec('si_core_mass', si_core_mass)
            have_fe_core_mass = get1_yvec('fe_core_mass', fe_core_mass)
         else
            have_he_core_mass = .false.
            have_c_core_mass = .false.
            have_o_core_mass = .false.
            have_si_core_mass = .false.
            have_fe_core_mass = .false.
         end if

         call pgsave
         call pgsch(txt_scale)

         dx = (xmax - xmin)/250.0
         call init_Kipp_plot
         call setup_mass_yaxis
         call plot_total_mass_line
         if (s% Kipp_show_burn) then
            call plot_burn_data(dx)
            if (s% Kipp_show_mixing) call plot_mix_data
         else if (s% Kipp_show_mixing) then
            call plot_mix_data
         end if
         if (s% Kipp_show_mass_boundaries) call plot_mass_lines
         if (s% Kipp_show_luminosities) call plot_L_lines

         call show_annotations(s, &
            s% show_Kipp_annotation1, &
            s% show_Kipp_annotation2, &
            s% show_Kipp_annotation3)

         call pgsch(txt_scale)
         call finish_Kipp_plot

         call pgunsa
         
         call dealloc


         contains
         
         subroutine dealloc
            deallocate(xvec, &
               log_L, &
               log_Lneu, &
               log_LH, &
               log_LHe, &
               star_mass, &
               log_xmstar, &
               star_M_center, &
               he_core_mass, &
               c_core_mass, &
               o_core_mass, &
               si_core_mass, &
               fe_core_mass)
         end subroutine dealloc

         subroutine init_Kipp_plot
            call pgsvp(vp_xleft, vp_xright, vp_ybot, vp_ytop)
         end subroutine init_Kipp_plot


         subroutine finish_Kipp_plot
            character(len=256) :: xlabel
            if (s% Kipp_show_luminosities) then
               call setup_L_yaxis
               call show_box_pgstar(s,'','CMSTV')
               call show_right_yaxis_label_pgstar(s,'log (L\d\(2281)\u)')
            end if
            call setup_mass_yaxis
            if (s% Kipp_show_luminosities) then
               call show_box_pgstar(s,'BCNST1','BNSTV1')
            else
               call show_box_pgstar(s,'BCNST1','BCNMSTV1')
            end if

            xlabel=''
            if (s% kipp_xaxis_log) then
               xlabel='log '// s% kipp_xaxis_name
            else
               xlabel=s% kipp_xaxis_name
               end if

            if (s% kipp_xaxis_name =='star_age') then
               if (s% kipp_xaxis_in_seconds) then
                  xlabel=trim(xlabel)//' (s)'
               else if (s% Kipp_xaxis_in_Myr) then
                  xlabel=trim(xlabel)//' (Myr)'
               else
                  xlabel=trim(xlabel)//' (yr)'
               end if
            end if

            call show_xaxis_label_pgstar(s,trim(xlabel),1.0)

            call show_left_yaxis_label_pgstar(s,'M/M\d\(2281)')
            if (.not. subplot) then
               call show_model_number_pgstar(s)
               call show_age_pgstar(s)
            end if
            call show_title_pgstar(s, title)
            call show_mix_legend
            call show_burn_legend
            
            call show_pgstar_decorator(s%id, s% kipp_use_decorator,s% kipp_pgstar_decorator, 0, ierr)

            
         end subroutine finish_Kipp_plot


         logical function get1_yvec(name, vec)
            character (len=*) :: name
            real, dimension(:), allocatable :: vec
            get1_yvec = get1_hist_yvec(s, step_min, step_max, n, name, vec)
         end function get1_yvec


         subroutine setup_L_yaxis
            real :: dy, ymin, ymax
            ymin = -2
            ymax = ymin
            if (have_log_L) ymax = maxval(log_L)
            if (have_log_LH) ymax = max(ymax, maxval(log_LH))
            if (have_log_LHe) ymax = max(ymax, maxval(log_LHe))
            if (ymax <= ymin) ymax = ymin+1
            if (s% Kipp_lgL_min /= -101d0) ymin = s% Kipp_lgL_min
            if (s% Kipp_lgL_max /= -101d0) ymax = s% Kipp_lgL_max
            dy = ymax - ymin
            if (s% Kipp_lgL_min /= -101d0) ymin = ymin - s% Kipp_lgL_margin*dy
            if (s% Kipp_lgL_max /= -101d0) ymax = ymax + s% Kipp_lgL_margin*dy
            call pgswin(xleft, xright, ymin, ymax)
            call pgscf(1)
            call pgsci(1)
            ymin_L_axis = ymin
            ymax_L_axis = ymax
         end subroutine setup_L_yaxis


         subroutine setup_mass_yaxis
            real :: dy, ymin, ymax
            ymax = s% Kipp_mass_max
            if (ymax <= 0) ymax = maxval(star_mass)
            ymin = s% Kipp_mass_min
            if (ymin < 0) ymin = minval(star_M_center)
            if (ymax <= ymin) ymax = ymin+1
            dy = ymax - ymin
            if (s% Kipp_mass_min /= -101d0) ymin = ymin - s% Kipp_mass_margin*dy
            if (s% Kipp_mass_max /= -101d0) ymax = ymax + s% Kipp_mass_margin*dy
            call pgswin(xleft, xright, ymin, ymax)
            call pgscf(1)
            call pgsci(1)
            ymin_mass_axis = ymin
            ymax_mass_axis = ymax
         end subroutine setup_mass_yaxis


         subroutine plot_burn_data(dx)
            use history_specs, only: burning_offset
            real, intent(in) :: dx
            type (pgstar_hist_node), pointer :: pg
            real :: burn_max, burn_min, burn_type, x, xnext
            integer :: i_burn_type_first, i_burn_type_last
            integer :: k, cnt, num_specs, step

            include 'formats'

            if (.not. s% Kipp_show_burn) return
            i_burn_type_first = 0
            num_specs = size(s% history_column_spec, dim=1)
            do k = 1, num_specs
               if (s% history_column_spec(k) == burning_offset + 1) then
                  i_burn_type_first = k
                  exit
               end if
            end do
            if (i_burn_type_first == 0) return

            i_burn_type_last = 0
            cnt = 1
            do k=i_burn_type_first+1, num_specs
               i_burn_type_last = k-1
               cnt = cnt+1
               if (s% history_column_spec(k) /= burning_offset + cnt) exit
            end do

            burn_max = 0.1
            burn_min = -0.1
            pg => s% pgstar_hist
            do
               if (.not. associated(pg)) exit
               step = pg% step
               if (step < step_min) exit
               if (step <= step_max) then
               do k = i_burn_type_first, i_burn_type_last, 2
                  burn_type = pg% vals(k)
                  if (burn_type < -9990) exit
                  if (burn_type /= 0) then
                     if (burn_type > 0) then
                        burn_max = max(burn_max, burn_type)
                     else if (burn_type < -1) then
                        burn_min = min(burn_min, burn_type)
                     end if
                  end if
                  if (pg% vals(k+1) >= 1d0) exit
               end do

               end if
               pg => pg% next
            end do

            call pgsave
            call pgslw(s% Kipp_burn_line_weight)
            pg => s% pgstar_hist
            xnext = xmax
            do
               if (.not. associated(pg)) exit
               step = pg% step
               if (step < step_min) exit
               x = real(xvec(step-step_min+1))
               if (step <= step_max .and. x <= xnext) then
                  call draw_burn_for_step( &
                     pg, i_burn_type_first, i_burn_type_last, &
                     burn_max, burn_min, burn_type_cutoff, x, &
                     star_mass(step-step_min+1), star_M_center(step-step_min+1))
                  xnext = max(x - dx, xmin)
               end if
               pg => pg% next
            end do
            call pgunsa
         end subroutine plot_burn_data


         subroutine draw_burn_for_step( &
               pg, i_burn_type_first, i_burn_type_last, &
               bmax, bmin, burn_type_cutoff, xval, mass, mass_center)

            type (pgstar_hist_node), pointer :: pg
            integer, intent(in) :: i_burn_type_first, i_burn_type_last
            real, intent(in) :: bmax, bmin, burn_type_cutoff, xval, mass, mass_center

            real :: burn_qbot, burn_qtop, mbot, mtop, xmass
            real :: burn_type, color_frac
            integer :: k, colormap_index, clr, mid_map
            include 'formats'
            xmass = mass - mass_center
            burn_qbot = 0
            mid_map = colormap_size/2
            do k = i_burn_type_first, i_burn_type_last, 2
               burn_type = pg% vals(k)
               burn_qtop = pg% vals(k+1)
               if (burn_type < -9990) exit
               if (abs(burn_type) > burn_type_cutoff) then
                  mbot = mass_center + xmass*burn_qbot
                  mtop = mass_center + xmass*burn_qtop
                  if (burn_type > 0.0) then
                     color_frac = 1.0 - max(0.0, min(1.0, burn_type/bmax))
                     colormap_index = &
                        colormap_size - int(0.6*color_frac*(colormap_size - mid_map))
                  else ! burn_type < 0.0
                     color_frac = 1.0 - max(0.0, min(1.0, burn_type/bmin))
                     colormap_index = 1 + int(0.6*color_frac*mid_map)
                  end if
                  clr = colormap_offset + colormap_index
                  call pgsci(clr)
                  call draw1(xval, mbot, mtop, clr)
               end if
               burn_qbot = burn_qtop
               if (burn_qbot >= 1d0) exit
            end do
         end subroutine draw_burn_for_step


         subroutine plot_L_lines
            integer :: i, cnt, n
            real :: coords(4), fjusts(4)

            logical, parameter :: dbg = .false.

            include 'formats'

            cnt = 0
            if (have_log_L) cnt = cnt + 1
            if (have_log_Lneu) cnt = cnt + 1
            if (have_log_LH) cnt = cnt + 1
            if (have_log_LHe) cnt = cnt + 1
            select case(cnt)
            case (1)
               coords(1) = 0.5; fjusts(1) = 0.5
            case (2)
               coords(1) = 0.80; fjusts(1) = 1.0
               coords(2) = 0.20; fjusts(2) = 0.0
            case (3)
               coords(1) = 0.90; fjusts(1) = 1.0
               coords(2) = 0.50; fjusts(2) = 0.5
               coords(3) = 0.10; fjusts(3) = 0.0
            case (4)
               coords(1) = 0.95; fjusts(1) = 1.0
               coords(2) = 0.75; fjusts(2) = 0.7
               coords(3) = 0.25; fjusts(3) = 0.5
               coords(4) = 0.05; fjusts(4) = 0.0
            case default
               return
            end select

            call pgsave

            call setup_L_yaxis

            call pgsch(txt_scale*0.8)
            call pgslw(2)

            n = 0

            if (have_log_L) then
               n = n+1
               call pgsci(clr_Crimson)
               call show_right_yaxis_label_pgmtxt_pgstar(s,coords(n),fjusts(n),'logL',-1.2)
               call plot_L_line(log_L)
            end if

            if (have_log_Lneu) then
               n = n+1
               call pgsci(clr_Tan)
               call show_right_yaxis_label_pgmtxt_pgstar(s,coords(n),fjusts(n),'logL\d\gn',-1.2)
               call plot_L_line(log_Lneu)
            end if

            if (have_log_LH) then
               n = n+1
               call pgsci(clr_Goldenrod)
               call show_right_yaxis_label_pgmtxt_pgstar(s,coords(n),fjusts(n),'logL\dHe',-1.2)
               call plot_L_line(log_LHe)
            end if

            if (have_log_LHe) then
               n = n+1
               call pgsci(clr_Silver)
               call show_right_yaxis_label_pgmtxt_pgstar(s,coords(n),fjusts(n),'logL\dH',-1.2)
               call plot_L_line(log_LH)
            end if

            call pgunsa

         end subroutine plot_L_lines


         subroutine plot_total_mass_line

            call pgsave
            call pgsch(txt_scale*0.8)

            call pgsci(clr_Gray)
            call pgslw(2)
            call show_left_yaxis_label_pgmtxt_pgstar(s,1.0,1.0,'M\dtotal\u',-0.8)
            call pgslw(s% pgstar_lw)
            call plot_mass_line(star_mass)

            call pgunsa

         end subroutine plot_total_mass_line


         subroutine plot_mass_lines
            integer :: i

            include 'formats'

            call pgsave
            call pgsch(txt_scale*0.8)
            call pgslw(2)

            if (have_he_core_mass) then
               call pgsci(clr_Teal)
               call show_left_yaxis_label_pgmtxt_pgstar(s,0.77,0.5,'He',-0.8)
               call plot_mass_line(he_core_mass)
            end if

            if (have_c_core_mass) then
               call pgsci(clr_LightOliveGreen)
               call show_left_yaxis_label_pgmtxt_pgstar(s,0.59,0.5,'C',-0.8)
               call plot_mass_line(c_core_mass)
            end if

            if (have_o_core_mass) then
               call pgsci(clr_SeaGreen)
               call show_left_yaxis_label_pgmtxt_pgstar(s,0.41,0.5,'O',-0.8)
               call plot_mass_line(o_core_mass)
            end if

            if (have_si_core_mass) then
               call pgsci(clr_Lilac)
               call show_left_yaxis_label_pgmtxt_pgstar(s,0.23,0.5,'Si',-0.8)
               call plot_mass_line(si_core_mass)
            end if

            if (have_fe_core_mass) then
               call pgsci(clr_Crimson)
               call show_left_yaxis_label_pgmtxt_pgstar(s,0.00,0.0,'Iron',-0.8)
               call plot_mass_line(fe_core_mass)
            end if

            call pgunsa

         end subroutine plot_mass_lines


         subroutine show_mix_legend
            integer :: mix_type

            call pgsave
            call pgslw(2)

            mix_type = convective_mixing
            call pgsci(mix_clr(mix_type))
            call show_xaxis_label_pgmtxt_pgstar(s, 0.05, 0.0, 'conv', 0.0)

            mix_type = leftover_convective_mixing
            call pgsci(mix_clr(mix_type))
            call show_xaxis_label_pgmtxt_pgstar(s, 0.275, 0.5, 'left', 0.0)

            mix_type = overshoot_mixing
            call pgsci(mix_clr(mix_type))
            call show_xaxis_label_pgmtxt_pgstar(s, 0.5, 0.5, 'over', 0.0)

            mix_type = semiconvective_mixing
            call pgsci(mix_clr(mix_type))
            call show_xaxis_label_pgmtxt_pgstar(s, 0.725, 0.5, 'semi', 0.0)

            mix_type = thermohaline_mixing
            call pgsci(mix_clr(mix_type))
            call show_xaxis_label_pgmtxt_pgstar(s, 0.95, 1.0, 'thrm', 0.0)

            mix_type = rayleigh_taylor_mixing
            call pgsci(mix_clr(mix_type))
            call show_xaxis_label_pgmtxt_pgstar(s, 0.85, 0.75, 'RTI', 0.0)

            call pgunsa

         end subroutine show_mix_legend


         subroutine show_burn_legend
            integer :: colormap_index

            call pgsave
            call pgslw(2)

            colormap_index = int(colormap_size*0.85)
            call pgsci(colormap_offset + colormap_index)
            call show_xaxis_label_pgmtxt_pgstar(s, 0.17, 0.5, 'burning', 1.0)

            colormap_index = int(colormap_size*0.15)
            call pgsci(colormap_offset + colormap_index)
            call show_xaxis_label_pgmtxt_pgstar(s, 0.82, 0.5, 'cooling', 1.0)

            call pgunsa

         end subroutine show_burn_legend


         subroutine plot_mass_line(yvec)
            real, intent(in) :: yvec(:)
            if (any(yvec > 1e-2*ymax_mass_axis)) then
               call pgsave
               call pgslw(s% Kipp_masses_line_weight)
               call pgline(n, xvec, yvec)
               call pgunsa
            end if
         end subroutine plot_mass_line


         subroutine plot_L_line(yvec)
            real, intent(in) :: yvec(:)
            if (any(yvec > ymin_L_axis)) then
               call pgsave
               call pgslw(s% Kipp_luminosities_line_weight)
               call pgline(n, xvec, yvec)
               call pgunsa
            end if
         end subroutine plot_L_line


         subroutine plot_mix_data
            use history_specs, only: mixing_offset
            type (pgstar_hist_node), pointer :: pg
            integer :: i_mix_type_first, i_mix_type_last
            integer :: k, cnt, num_specs, step

            include 'formats'

            i_mix_type_first = 0
            num_specs = size(s% history_column_spec, dim=1)
            do k = 1, num_specs
               if (s% history_column_spec(k) == mixing_offset + 1) then
                  i_mix_type_first = k
                  exit
               end if
            end do
            if (i_mix_type_first == 0) then
               !write(*,*) 'i_mix_type_first == 0'
               return
            end if

            i_mix_type_last = 0
            cnt = 1
            do k=i_mix_type_first+1, num_specs
               i_mix_type_last = k-1
               cnt = cnt+1
               if (s% history_column_spec(k) /= mixing_offset + cnt) exit
            end do

            call pgsave
            call pgslw(s% Kipp_mix_line_weight)
            if (.not. associated(s% pgstar_hist)) then
               write(*,*) '.not. associated(s% pgstar_hist)'
            end if
            pg => s% pgstar_hist
            do
               if (.not. associated(pg)) exit
               step = pg% step
               if (step < step_min) exit
               if (step <= step_max .and. mod(step, s% Kipp_mix_interval) == 0) then
                  call draw_mix_for_step( &
                     pg, step, i_mix_type_first, i_mix_type_last, real(xvec(step-step_min+1)), &
                     star_mass(step-step_min+1), star_M_center(step-step_min+1))
               end if
               pg => pg% next
            end do
            call pgunsa
         end subroutine plot_mix_data


         subroutine draw_mix_for_step( &
               pg, step, i_mix_type_first, i_mix_type_last, xval, mass, mass_center)
            type (pgstar_hist_node), pointer :: pg
            integer, intent(in) :: step, i_mix_type_first, i_mix_type_last
            real, intent(in) :: xval, mass, mass_center
            real :: mix_qbot, mix_qtop, mbot, mtop, xmass
            integer :: k, mix_type
            include 'formats'
            mix_qbot = 0
            xmass = mass - mass_center
            do k = i_mix_type_first, i_mix_type_last, 2
               mix_type = int(pg% vals(k))
               if (mix_type < 0) exit
               mix_qtop = pg% vals(k+1)
               if (mix_type > 0 .and. mix_type <= max_mix_type_to_show) then
                  mbot = mass_center + xmass*mix_qbot
                  mtop = mass_center + xmass*mix_qtop
                  call draw1(xval, mbot, mtop, mix_clr(mix_type))
                  showed_this_mix_type(mix_type) = .true.
               end if
               mix_qbot = mix_qtop
            end do
         end subroutine draw_mix_for_step


         subroutine draw1(xval,y1,y2,clr)
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


      end subroutine do_Kipp_Plot



      end module pgstar_kipp

