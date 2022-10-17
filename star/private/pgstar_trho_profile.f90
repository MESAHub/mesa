! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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

module pgstar_trho_profile
   
   use star_private_def
   use const_def
   use pgstar_support
   use star_pgstar
   
   implicit none


contains
   
   
   subroutine TRho_Profile_plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      
      ierr = 0
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) return
      
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      
      call do_TRho_Profile_plot(s, id, device_id, &
         s% pg% TRho_Profile_xleft, s% pg% TRho_Profile_xright, &
         s% pg% TRho_Profile_ybot, s% pg% TRho_Profile_ytop, .false., &
         s% pg% TRho_Profile_title, s% pg% TRho_Profile_txt_scale, ierr)
      
      call pgebuf()
   
   end subroutine TRho_Profile_plot
   
   
   subroutine do_TRho_Profile_plot(s, id, device_id, &
      xleft, xright, ybot, ytop, subplot, title, txt_scale_in, ierr)
      use utils_lib
      
      type (star_info), pointer :: s
      integer, intent(in) :: id, device_id
      real, intent(in) :: xleft, xright, ybot, ytop, txt_scale_in
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      
      integer :: nz, k
      real :: xmin, xmax, ymin, ymax, xpos, ypos, dx, dy, &
         txt_scale, vpxmin, vpxmax, vpymin, vpymax, vpymargin, vpwinheight, lgT1, lgT2
      real, allocatable, dimension(:) :: xvec, yvec
      character (len = 128) :: str
      real, parameter :: lgrho1 = -8, lgrho2 = 5
      
      include 'formats'
      
      ierr = 0
      nz = s% nz
      allocate (xvec(nz), yvec(nz))
      
      txt_scale = txt_scale_in
      
      if (s% pg% TRho_switch_to_Column_Depth) then
         do k = 1, nz
            xvec(k) = safe_log10(s% xmstar * sum(s% dq(1:k - 1)) / (pi4 * s% r(k) * s% r(k)))
         end do
      else ! log rho
         do k = 1, nz
            xvec(k) = s% lnd(k) / ln10
         end do
      end if
      if (s% pg% TRho_switch_to_mass) then
         do k = 1, nz
            xvec(k) = safe_log10((s% xmstar - s% m(k)) / Msun)
         end do
      end if
      xmin = s% pg% TRho_Profile_xmin
      xmax = s% pg% TRho_Profile_xmax
      dx = xmax - xmin
      
      call pgsave
      call pgsch(txt_scale)
      
      ! log T
      do k = 1, nz
         yvec(k) = s% lnT(k) / ln10
      end do
      ymin = s% pg% TRho_Profile_ymin
      ymax = s% pg% TRho_Profile_ymax
      dy = ymax - ymin
      
      call pgsvp(xleft, xright, ybot, ytop)
      call pgswin(xmin, xmax, ymin, ymax)
      call pgscf(1)
      call pgsci(1)
      call show_box_pgstar(s, 'BCNST1', 'BCMNSTV1')
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (s% pg% show_TRho_accretion_mesh_borders) then
         if(s% pg% TRho_switch_to_mass) then
            call do_accretion_mesh_borders(safe_log10((s% xmstar&
               - s% m(s% k_const_mass)) / Msun), &
               safe_log10((s% xmstar&
                  - s% m(s% k_below_const_q)) / Msun), &
               safe_log10((s% xmstar&
                  - s% m(s% k_below_just_added)) / Msun), &
               ymin, ymax)
         end if
         if(s% pg% TRho_switch_to_Column_Depth) then
            call do_accretion_mesh_borders(safe_log10(s% xmstar * sum(s% &
               dq(1:s% k_const_mass - 1)) / (pi4 * s% r(s% k_const_mass)&
               * s% r(s% k_const_mass))), &
               safe_log10(s% xmstar * sum(s% &
                  dq(1:s% k_below_const_q - 1)) / (pi4 * s% r(s% k_below_const_q)&
                  * s% r(s% k_below_const_q))), &
               safe_log10(s% xmstar * sum(s% &
                  dq(1:s% k_below_just_added - 1)) / (pi4 * s% r(s% k_below_just_added)&
                  * s% r(s% k_below_just_added))), &
               ymin, ymax)
         end if
         
         if(.not. s% pg% TRho_switch_to_Column_Depth .and. .not. s% pg% TRho_switch_to_mass) then
            call do_accretion_mesh_borders(s% lnd(s% k_const_mass) / ln10, &
               s% lnd(s% k_below_const_q) / ln10, &
               s% lnd(s% k_below_just_added) / ln10, &
               ymin, ymax)
         end if
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      end if
      if (s% pg% TRho_switch_to_Column_Depth) then
         call show_xaxis_label_pgstar(s, 'log column depth (g cm\u-2\d)')
      end if
      if(.not. s% pg% TRho_switch_to_Column_Depth .and. .not. s% pg% &
         TRho_switch_to_mass) then
         call show_xaxis_label_pgstar(s, 'log Density (g cm\u-3\d)')
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(s% pg% TRho_switch_to_mass) then
         call show_xaxis_label_pgstar(s, 'log M - m (Msun)')
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call show_left_yaxis_label_pgstar(s, 'log Temperature (K)')
      
      if (.not. subplot) then
         call show_model_number_pgstar(s)
         call show_age_pgstar(s)
      end if
      call show_title_pgstar(s, title)
      
      if (.not. s% pg% TRho_switch_to_Column_Depth .and. .not. s% pg% TRho_switch_to_mass) then
         if (s% pg% show_TRho_Profile_kap_regions) call do_kap_regions
         if (s% pg% show_TRho_Profile_eos_regions) call do_eos_regions
         ! for now, show eos regions will imply showing gamma1 4/3 also
         if (s% pg% show_TRho_Profile_gamma1_4_3rd .or. s% pg% show_TRho_Profile_eos_regions) call do_gamma1_4_3rd
         if (s% pg% show_TRho_Profile_degeneracy_line) call do_degeneracy_line
         if (s% pg% show_TRho_Profile_Pgas_Prad_line) call do_Pgas_Prad_line
         if (s% pg% show_TRho_Profile_burn_lines) call do_burn_lines
      end if
      
      if (len_trim(s% pg% TRho_Profile_fname) > 0) then
         
         call mesa_error(__FILE__, __LINE__, 'NEED TO ADD ABILITY TO SHOW EXTRA PROFILE FOR COMPARISON')
      
      end if
      
      call show_profile_line(s, xvec, yvec, txt_scale, xmin, xmax, ymin, ymax, &
         s% pg% show_TRho_Profile_legend, s% pg% TRho_Profile_legend_coord, &
         s% pg% TRho_Profile_legend_disp1, s% pg% TRho_Profile_legend_del_disp, &
         s% pg% TRho_Profile_legend_fjust, &
         s% pg% show_TRho_Profile_mass_locs)
      
      if (s% pg% show_TRho_Profile_text_info) &
         call do_show_Profile_text_info(&
            s, txt_scale, xmin, xmax, ymin, ymax, &
            s% pg% TRho_Profile_text_info_xfac, s% pg% TRho_Profile_text_info_dxfac, &
            s% pg% TRho_Profile_text_info_yfac, s% pg% TRho_Profile_text_info_dyfac, &
            .false., .false.)
      
      call show_annotations(s, &
         s% pg% show_TRho_Profile_annotation1, &
         s% pg% show_TRho_Profile_annotation2, &
         s% pg% show_TRho_Profile_annotation3)
      
      deallocate(xvec, yvec)
      
      call show_pgstar_decorator(s%id, s% pg% TRho_Profile_use_decorator, &
         s% pg% TRho_Profile_pgstar_decorator, 0, ierr)
      
      call pgunsa
   
   
   contains
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine do_accretion_mesh_borders(x_Lagrange, x_Eulerian, x_just_added, min_T, max_T)
         real(dp), intent(in) :: x_Lagrange, x_Eulerian, x_just_added
         real, intent(in) :: min_T, max_T
         call pgsci(clr_RoyalPurple)
         call stroke_line(real(x_Lagrange), min_T, real(x_Lagrange), max_T)
         call pgsci(clr_RoyalBlue)
         call stroke_line(real(x_Eulerian), min_T, real(x_Eulerian), max_T)
         call pgsci(clr_Tan)
         call stroke_line(real(x_just_added), min_T, real(x_just_added), max_T)
         call pgsci(clr_Gray)
      end subroutine do_accretion_mesh_borders
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      
      subroutine do_degeneracy_line
         call pgsave
         call pgsch(txt_scale * 0.9)
         call pgsci(clr_Gray)
         call pgsls(Line_Type_Dash)
         call pgline(size(psi4_logT), psi4_logRho, psi4_logT)
         call pgsls(Line_Type_Solid)
         xpos = -0.2 ! 1.9 ! psi4_logRho(1)
         ypos = 4 ! 5.9 ! psi4_logT(1)-dy*0.04
         if (inside(xpos, ypos)) call pgptxt(xpos, ypos, 0.0, 0.5, '\ge\dF\u/kT\(0248)4')
         call pgunsa
      end subroutine do_degeneracy_line
      
      
      subroutine add_TR_line(logR1, logT1, logR2, logT2)
         real, intent(in) :: logR1, logT1, logR2, logT2
         real :: logRho1, logRho2
         logRho1 = logR1 + 3 * logT1 - 18
         logRho2 = logR2 + 3 * logT2 - 18
         call pgmove(logRho1, logT1)
         call pgdraw(logRho2, logT2)
      end subroutine add_TR_line
      
      
      subroutine show_label(xpos, ypos, angle, justification, txt)
         real, intent(in) :: xpos, ypos, angle, justification
         character (len = *), intent(in) :: txt
         if (inside(xpos, ypos)) call pgptxt(xpos, ypos, angle, justification, txt)
      end subroutine show_label
      
      
      subroutine do_kap_regions
         real :: logT_lo, logT_hi, logT_max, logR, logRho_lo, logRho_hi, logRho_min
         real, parameter :: min_logR_for_freedman = 1
         real, parameter :: freg_blend_logT2 = 4.10
         real, parameter :: freg_blend_logT1 = 3.93
         
         call pgsave
         
         call pgsci(clr_Coral)
         call pgsls(Line_Type_Solid)
         logT_lo = 2.7; logT_hi = 8.7; logT_max = 10.3
         call add_TR_line(-8.0, logT_lo, -8.0, logT_hi)
         call add_TR_line(1.0, logT_lo, 1.0, logT_hi)
         call add_TR_line(1.0, logT_lo, -8.0, logT_lo)
         call add_TR_line(1.0, logT_hi, -8.0, logT_hi)
         call add_TR_line(1.0, 2.7, -8.0, 2.7)
         call add_TR_line(1.0, freg_blend_logT1, -8.0, freg_blend_logT1)
         call add_TR_line(1.0, freg_blend_logT2, -8.0, freg_blend_logT2)
         call add_TR_line(1.0, 8.2, -8.0, 8.2)
         
         call pgsci(1)
         call add_TR_line(-8.0, logT_hi, -8.0, logT_max)
         call add_TR_line(-8.0, logT_max, 8.0, logT_max)
         !call add_TR_line(8.0, logT_lo, 8.0, logT_hi)
         !call add_TR_line(1.0, logT_lo, 8.0, logT_lo)
         
         ! Freedman
         call pgsci(clr_Tan)
         call pgmove(-8.8, 1.88)
         call pgdraw(-3.36, 1.88)
         call pgdraw(-1.5, 2.5)
         call pgdraw(-2.6, 3.6)
         call pgdraw(-11.3, 3.6)
         call pgdraw(-9.5, 1.88)
         call pgdraw(-8.8, 1.88)
         
         call pgsci(1)
         call show_label(-4.9, 2.47, 0.0, 0.5, 'FREEDMAN')
         call show_label(-8.5, 3.3, 0.0, 0.5, 'FERGUSON')
         call show_label(-7.5, 5.1, 0.0, 0.5, 'OPAL/OP')
         call show_label(5.5, 9.0, 0.0, 0.5, 'COMPTON')
         call show_label(1.8, 8.35, 0.0, 0.5, 'BLEND')
         call show_label(-8.5, (freg_blend_logT1 + freg_blend_logT2) / 2, 0.0, 0.5, 'BLEND')
         call show_label(0.2, 3.9, 0.0, 1.0, '\(0636)\drad\u = \(0636)\dcond\u')
         call pgsci(clr_Crimson)
         call show_label(3.8, 9.4, 0.0, 0.5, 'e\u-\de\u+\d')
         call pgsci(1)
         
         call show_label(-6.8, 6.9, 0.0, 0.5, 'logR = -8')
         call show_label(5.0, 6.9, 0.0, 0.5, 'logR = 1')
         call show_label(2.8, 3.8, 0.0, 0.5, 'logR = 8')
         
         ! show where electron to baryon ratio is twice that expected
         call pgsci(clr_Crimson)
         call pgsls(Line_Type_Dash)
         call pgline(size(elect_data_logT), elect_data_logRho, elect_data_logT)
         ! show where kap_rad == kap_cond
         call pgsci(clr_LightSkyBlue)
         call pgsls(Line_Type_Dot)
         call pgline(size(kap_rad_cond_eq_logT), kap_rad_cond_eq_logRho, kap_rad_cond_eq_logT)
         call pgunsa
      end subroutine do_kap_regions
      
      
      subroutine do_eos_regions
         integer :: ierr
         real :: logRho0, logRho1, logRho2, logRho3, logRho4, logRho5, logRho6
         real :: logT1, logT2, logT3, logT4, logT5, logT6
         
         call pgsave
         
         ! blend from table to non-table
         call pgsci(clr_LightSkyGreen)
         call pgsls(Line_Type_Dash)
         
         logT1 = s% eos_rq% logT_min_for_all_Skye
         logT2 = s% eos_rq% logT_min_for_any_Skye
         logT3 = 0 ! s% eos_rq% logT_min_FreeEOS_lo
         logT4 = 0 ! s% eos_rq% logT_min_FreeEOS_lo
         
         logRho1 = s% eos_rq% logRho_min_for_all_Skye
         logRho2 = s% eos_rq% logRho_min_for_any_Skye
         logRho3 = s% eos_rq% logQ_min_FreeEOS_lo + 2 * logT1 - 12
         logRho4 = s% eos_rq% logQ_min_FreeEOS_hi + 2 * logT2 - 12
         logRho5 = s% eos_rq% logQ_min_FreeEOS_lo + 2 * logT3 - 12
         logRho6 = s% eos_rq% logQ_min_FreeEOS_hi + 2 * logT4 - 12
         
         call stroke_line(logRho1, logT1, logRho3, logT1)
         call stroke_line(logRho2, logT2, logRho4, logT2)
         call stroke_line(logRho3, logT1, logRho5, logT3)
         call stroke_line(logRho4, logT2, logRho6, logT4)
         
         call stroke_line(logRho1, logT1, logRho1, logT3)
         call stroke_line(logRho2, logT2, logRho2, logT4)
         
         ! blend from OPAL to SCVH
         call pgsci(clr_LightSkyBlue)
         call pgsls(Line_Type_Dot)
         
         logRho0 = logRho1
         
         logT1 = s% eos_rq% logT_cut_FreeEOS_hi
         logT2 = s% eos_rq% logT_cut_FreeEOS_lo
         logT3 = s% eos_rq% logT_min_FreeEOS_hi
         logT4 = s% eos_rq% logT_min_FreeEOS_lo
         logT5 = 0.5 * (logRho0 - s% eos_rq% logQ_max_OPAL_SCVH + 12)
         logT6 = s% eos_rq% logT_low_all_HELM
         
         logRho1 = s% eos_rq% logQ_cut_lo_Z_FreeEOS_hi + 2 * logT1 - 12
         logRho2 = s% eos_rq% logQ_cut_lo_Z_FreeEOS_lo + 2 * logT2 - 12
         logRho3 = s% eos_rq% logQ_cut_lo_Z_FreeEOS_hi + 2 * logT3 - 12
         logRho4 = s% eos_rq% logQ_cut_lo_Z_FreeEOS_lo + 2 * logT4 - 12
         logRho5 = s% eos_rq% logRho_min_OPAL_SCVH_limit
         logRho6 = s% eos_rq% logQ_max_OPAL_SCVH + 2 * logT6 - 12
         
         call stroke_line(logRho0, logT1, logRho2, logT1)
         call stroke_line(logRho2, logT1, logRho4, logT3)
         call stroke_line(logRho4, logT3, logRho5, logT3)
         
         call stroke_line(logRho0, logT2, logRho1, logT2)
         call stroke_line(logRho1, logT2, logRho3, logT4)
         call stroke_line(logRho3, logT4, logRho5, logT4)
         
         call stroke_line(logRho0, logT5, logRho6, logT6)
         call stroke_line(logRho5, logT6, logRho6, logT6)
         
         call pgsci(1)
         call show_label(1.0, 3.2, 0.0, 0.5, 'HELM')
         call show_label(-7.2, 5.8, 0.0, 0.5, 'FreeEOS')
         call show_label(-1.5, 3.7, 0.0, 0.5, 'OPAL/SCVH')
         call show_label(-1.5, 9.7, 0.0, 0.5, 'HELM/Skye EOS')
         call show_label(6.0, 4.5, 0.0, 0.5, 'Skye EOS')
         
         call pgunsa
      end subroutine do_eos_regions
      
      
      subroutine do_gamma1_4_3rd
         call pgsave
         ! show where gamma1 = 4/3
         call pgsci(clr_Gold)
         call pgsls(Line_Type_Solid)
         call pgslw(3)
         call show_label(3.0, 9.3, 0.0, 0.5, '\(0529)\d1\u < 4/3')
         call pgslw(4)
         call pgline(size(gamma_4_thirds_logT), gamma_4_thirds_logRho, gamma_4_thirds_logT)
         call pgunsa
      end subroutine do_gamma1_4_3rd
      
      
      subroutine stroke_line(x1, y1, x2, y2)
         real, intent(in) :: x1, y1, x2, y2
         call pgmove(x1, y1)
         call pgdraw(x2, y2)
      end subroutine stroke_line
      
      
      subroutine do_Pgas_Prad_line
         lgT1 = log10(3.2d7) + (lgRho1 - log10(0.7d0)) / 3.0
         lgT2 = log10(3.2d7) + (lgRho2 - log10(0.7d0)) / 3.0
         call pgsave
         call pgsch(txt_scale * 0.9)
         call pgsci(clr_Gray)
         call pgsls(Line_Type_Dash)
         call pgmove(lgRho1, lgT1)
         call pgdraw(lgRho2, lgT2)
         call pgsls(Line_Type_Solid)
         xpos = -4 ! lgRho1-dx*0.065
         ypos = 6.5 ! lgT1-dy*0.025
         if (inside(xpos, ypos)) call pgptxt(xpos, ypos, 0.0, 0.0, 'P\drad\u\(0248)P\dgas\u')
         call pgunsa
      end subroutine do_Pgas_Prad_line
      
      
      subroutine do_burn_lines
         call pgsave
         call pgsch(txt_scale * 0.9)
         call pgsci(clr_Gray)
         call pgsls(Line_Type_Dash)
         call write_burn_line(hydrogen_burn_logRho, hydrogen_burn_logT, 'H burn')
         call write_burn_line(helium_burn_logRho, helium_burn_logT, 'He burn')
         call write_burn_line(carbon_burn_logRho, carbon_burn_logT, 'C burn')
         call write_burn_line(oxygen_burn_logRho, oxygen_burn_logT, 'O burn')
         call pgsls(Line_Type_Solid)
         call pgunsa
      end subroutine do_burn_lines
      
      
      logical function inside(xpos, ypos)
         real, intent(in) :: xpos, ypos
         inside = .false.
         if (xpos <= s% pg% TRho_Profile_xmin .or. xpos >= s% pg% TRho_Profile_xmax) return
         if (ypos <= s% pg% TRho_Profile_ymin .or. ypos >= s% pg% TRho_Profile_ymax) return
         inside = .true.
      end function inside
      
      
      subroutine write_burn_line(logRho, logT, label)
         real, dimension(:), allocatable :: logRho, logT
         character (len = *), intent(in) :: label
         integer :: sz
         real :: xpos, ypos
         character (len = 128) :: str
         sz = size(logRho)
         call pgline(sz, logRho, logT)
         if (.not. s% pg% show_TRho_Profile_burn_labels) return
         xpos = logRho(sz)
         ypos = logT(sz)
         if (.not. inside(xpos, ypos)) return
         write(str, '(a)') trim(label)
         call pgptxt(xpos, ypos, 0.0, 1.0, trim(adjustl(str)))
      end subroutine write_burn_line
   
   
   end subroutine do_TRho_Profile_plot
   
   
   subroutine do_show_Profile_text_info(&
      s, txt_scale, xmin, xmax, ymin, ymax, xfac, dxfac, yfac, dyfac, &
      xaxis_reversed, yaxis_reversed)
      type (star_info), pointer :: s
      real, intent(in) :: txt_scale, xmin, xmax, ymin, ymax, xfac, dxfac, yfac, dyfac
      logical, intent(in) :: xaxis_reversed, yaxis_reversed
      
      real :: dxpos, xpos0, dxval, ypos, dypos
      real(dp) :: age
      integer :: cnt
      
      include 'formats'
      
      call pgsave
      call pgsch(0.7 * txt_scale)
      call pgsci(1)
      dxpos = 0
      xpos0 = xmin + xfac * (xmax - xmin)
      dxval = dxfac * (xmax - xmin)
      if (xaxis_reversed) dxval = -dxval
      ypos = ymin + yfac * (ymax - ymin)
      dypos = dyfac * (ymax - ymin)
      if (yaxis_reversed) dypos = -dypos
      
      cnt = 0; ypos = ypos + dypos
      cnt = write_info_line_flt(cnt, ypos, xpos0, dxpos, dxval, &
         'mass', s% star_mass)
      
      cnt = 0; ypos = ypos + dypos
      cnt = write_info_line_flt(cnt, ypos, xpos0, dxpos, dxval, &
         'H rich', s% star_mass - s% he_core_mass)
      
      cnt = 0; ypos = ypos + dypos
      cnt = write_info_line_flt(cnt, ypos, xpos0, dxpos, dxval, &
         'He core', s% he_core_mass)
      
      cnt = 0; ypos = ypos + dypos
      cnt = write_info_line_flt(cnt, ypos, xpos0, dxpos, dxval, &
         'CO core', s% co_core_mass)
      
      cnt = 0; ypos = ypos + dypos
      cnt = write_info_line_flt(cnt, ypos, xpos0, dxpos, dxval, &
         'lg mdot', safe_log10(abs(s% star_mdot)))
      
      cnt = 0; ypos = ypos + dypos
      cnt = write_info_line_flt2(cnt, ypos, xpos0, dxpos, dxval, &
         'Teff', s% Teff)
      
      cnt = 0; ypos = ypos + dypos
      cnt = write_info_line_flt(cnt, ypos, xpos0, dxpos, dxval, &
         'lg L', s% log_surface_luminosity)
      
      cnt = 0; ypos = ypos + dypos
      cnt = write_info_line_flt(cnt, ypos, xpos0, dxpos, dxval, &
         'lg LH', safe_log10(s% power_h_burn))
      
      cnt = 0; ypos = ypos + dypos
      cnt = write_info_line_flt(cnt, ypos, xpos0, dxpos, dxval, &
         'lg LHe', safe_log10(s% power_he_burn))
      
      cnt = 0; ypos = ypos + dypos
      cnt = write_info_line_flt(cnt, ypos, xpos0, dxpos, dxval, &
         'lg R', s% log_surface_radius)
      
      cnt = 0; ypos = ypos + dypos
      cnt = write_info_line_flt(cnt, ypos, xpos0, dxpos, dxval, &
         'max lg T', s% log_max_temperature)
      
      cnt = 0; ypos = ypos + dypos
      cnt = write_info_line_exp(cnt, ypos, xpos0, dxpos, dxval, &
         'lg dt yr', log10(s% time_step))
      
      cnt = 0; ypos = ypos + dypos
      age = s% star_age
      if (s% pg% pgstar_show_age_in_seconds) then
         cnt = write_info_line_exp(cnt, ypos, xpos0, dxpos, dxval, &
            'age sec', age * secyer)
      else
         cnt = write_info_line_exp(cnt, ypos, xpos0, dxpos, dxval, &
            'age yr', age)
      end if
      
      call pgunsa
   
   end subroutine do_show_Profile_text_info


end module pgstar_trho_profile

