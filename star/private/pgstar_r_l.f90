! ***********************************************************************
!
!   Copyright (C) 2014  Bill Paxton
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

      module pgstar_r_l

      use star_private_def
      use const_def
      use pgstar_support

      implicit none


      contains


      subroutine R_L_Plot(id, device_id, ierr)
         integer, intent(in) :: id, device_id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         call pgslct(device_id)
         call pgbbuf()
         call pgeras()

         call do_R_L_Plot(s, id, device_id, &
            s% R_L_xleft, s% R_L_xright, &
            s% R_L_ybot, s% R_L_ytop, .false., &
            s% R_L_title, s% R_L_txt_scale, ierr)
         if (ierr /= 0) return

         call pgebuf()

      end subroutine R_L_Plot


      subroutine do_R_L_Plot(s, id, device_id, &
            xleft, xright, ybot, ytop, subplot, &
            title, txt_scale, ierr)
         use pgstar_hist_track, only: null_decorate, do_Hist_Track
         type (star_info), pointer :: s
         integer, intent(in) :: id, device_id
         real, intent(in) :: xleft, xright, ybot, ytop, txt_scale
         logical, intent(in) :: subplot
         character (len=*), intent(in) :: title
         integer, intent(out) :: ierr
         logical, parameter :: &
            reverse_xaxis = .false., reverse_yaxis = .false.
         ierr = 0
         call do_Hist_Track(s, id, device_id, &
            xleft, xright, ybot, ytop, subplot, title, txt_scale, &
            'luminosity', 'radius', &
            'L/L\d\(2281)', 'R/R\d\(2281)', &
            s% R_L_L_min, s% R_L_L_max, &
            s% R_L_L_margin, s% R_L_dL_min, &
            s% R_L_R_min, s% R_L_R_max, &
            s% R_L_R_margin, s% R_L_dR_min, &
            s% R_L_step_min, s% R_L_step_max, &
            reverse_xaxis, reverse_yaxis, .false., .false., &
            s% show_R_L_target_box, s% R_L_target_n_sigma, &
            s% R_L_target_L, s% R_L_target_R, &
            s% R_L_target_L_sigma, s% R_L_target_R_sigma, &
            s% show_R_L_annotation1, &
            s% show_R_L_annotation2, &
            s% show_R_L_annotation3, &
            s% R_L_fname, &
            s% R_L_use_decorator, &
            s% R_L_pgstar_decorator, &
            null_decorate, ierr)
      end subroutine do_R_L_Plot


      end module pgstar_r_l

