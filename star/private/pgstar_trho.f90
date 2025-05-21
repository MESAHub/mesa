! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
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

      module pgstar_trho

      use star_private_def
      use const_def, only: dp
      use pgstar_support
      use star_pgstar

      implicit none


      contains


      subroutine TRho_Plot(id, device_id, ierr)
         integer, intent(in) :: id, device_id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         call pgslct(device_id)
         call pgbbuf()
         call pgeras()

         call do_TRho_Plot(s, id, device_id, &
            s% pg% TRho_xleft, s% pg% TRho_xright, &
            s% pg% TRho_ybot, s% pg% TRho_ytop, .false., &
            s% pg% TRho_title, s% pg% TRho_txt_scale, ierr)
         if (ierr /= 0) return

         call pgebuf()

      end subroutine TRho_Plot


      subroutine do_TRho_Plot(s, id, device_id, &
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
            'log_center_Rho', 'log_center_T', &
            'log Rho\dc\u (g cm\u-3\d)', 'log T\dc\u (K)', &
            s% pg% TRho_logRho_min, s% pg% TRho_logRho_max, &
            s% pg% TRho_logRho_margin, s% pg% TRho_logRho_dlogRho_min, &
            s% pg% TRho_logT_min, s% pg% TRho_logT_max, &
            s% pg% TRho_logT_margin, s% pg% TRho_logT_dlogT_min, &
            s% pg% TRho_step_min, s% pg% TRho_step_max, &
            reverse_xaxis, reverse_yaxis, .false., .false., &
            .false., 0, 0.0, 0.0, 0.0, 0.0, &
            s% pg% show_TRho_annotation1, &
            s% pg% show_TRho_annotation2, &
            s% pg% show_TRho_annotation3, &
            s% pg% TRho_fname, &
            s% pg% TRho_use_decorator, &
            s% pg% TRho_pgstar_decorator, &
            do_degeneracy_line, ierr)
      end subroutine do_TRho_Plot


      subroutine do_degeneracy_line(id, ierr)
         use pgstar_colors
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (.not. s% pg% show_TRho_degeneracy_line) return
         call pgsave
         call pgsci(clr_Gray)
         call pgsls(Line_Type_Dash)
         call pgline(size(psi4_logT), psi4_logRho, psi4_logT)
         call pgsls(Line_Type_Solid)
         call pgunsa
      end subroutine do_degeneracy_line


      end module pgstar_trho

