! ***********************************************************************
!
!   Copyright (C) 2013  The MESA Team
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

      module pgstar_dPg_dnu

      use star_private_def
      use const_def, only: dp
      use pgstar_support
      use star_pgstar

      implicit none

      contains

      subroutine dPg_dnu_Plot(id, device_id, ierr)
         integer, intent(in) :: id, device_id
         integer, intent(out) :: ierr

         type (star_info), pointer :: s

         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         call pgslct(device_id)
         call pgbbuf()
         call pgeras()

         call do_dPg_dnu_Plot(s, id, device_id, &
            s% pg% dPg_dnu_xleft, s% pg% dPg_dnu_xright, &
            s% pg% dPg_dnu_ybot, s% pg% dPg_dnu_ytop, .false., &
            s% pg% dPg_dnu_title, s% pg% dPg_dnu_txt_scale, ierr)

         call pgebuf()

      end subroutine dPg_dnu_Plot


      subroutine do_dPg_dnu_Plot(s, id, device_id, &
            xleft, xright, ybot, ytop, subplot, title, txt_scale, ierr)
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
            'delta_nu', 'delta_Pg', &
            'delta nu', 'delta Pg', &
            s% pg% dPg_dnu_delta_nu_min, s% pg% dPg_dnu_delta_nu_max, &
            s% pg% dPg_dnu_delta_nu_margin, s% pg% dPg_dnu_d_delta_nu_min, &
            s% pg% dPg_dnu_delta_Pg_min, s% pg% dPg_dnu_delta_Pg_max, &
            s% pg% dPg_dnu_delta_Pg_margin, s% pg% dPg_dnu_d_delta_Pg_min, &
            s% pg% dPg_dnu_step_min, s% pg% dPg_dnu_step_max, &
            reverse_xaxis, reverse_yaxis, .false., .false., &
            .false., 0, 0.0, 0.0, 0.0, 0.0, &
            s% pg% show_dPg_dnu_annotation1, &
            s% pg% show_dPg_dnu_annotation2, &
            s% pg% show_dPg_dnu_annotation3, &
            s% pg% dPg_dnu_fname, &
            s% pg% dPg_dnu_use_decorator, s% pg% dPg_dnu_pgstar_decorator, &
            null_decorate, ierr)
      end subroutine do_dPg_dnu_Plot

      end module pgstar_dPg_dnu
