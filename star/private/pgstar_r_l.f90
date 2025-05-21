! ***********************************************************************
!
!   Copyright (C) 2014  The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! ***********************************************************************

      module pgstar_r_l

      use star_private_def
      use const_def, only: dp
      use pgstar_support
      use star_pgstar

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
            s% pg% R_L_xleft, s% pg% R_L_xright, &
            s% pg% R_L_ybot, s% pg% R_L_ytop, .false., &
            s% pg% R_L_title, s% pg% R_L_txt_scale, ierr)
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
            s% pg% R_L_L_min, s% pg% R_L_L_max, &
            s% pg% R_L_L_margin, s% pg% R_L_dL_min, &
            s% pg% R_L_R_min, s% pg% R_L_R_max, &
            s% pg% R_L_R_margin, s% pg% R_L_dR_min, &
            s% pg% R_L_step_min, s% pg% R_L_step_max, &
            reverse_xaxis, reverse_yaxis, .false., .false., &
            s% pg% show_R_L_target_box, s% pg% R_L_target_n_sigma, &
            s% pg% R_L_target_L, s% pg% R_L_target_R, &
            s% pg% R_L_target_L_sigma, s% pg% R_L_target_R_sigma, &
            s% pg% show_R_L_annotation1, &
            s% pg% show_R_L_annotation2, &
            s% pg% show_R_L_annotation3, &
            s% pg% R_L_fname, &
            s% pg% R_L_use_decorator, &
            s% pg% R_L_pgstar_decorator, &
            null_decorate, ierr)
      end subroutine do_R_L_Plot

      end module pgstar_r_l
