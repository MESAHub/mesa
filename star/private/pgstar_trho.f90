! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
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

      module pgstar_trho

      use star_private_def
      use const_def
      use pgstar_support

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
            s% TRho_xleft, s% TRho_xright, &
            s% TRho_ybot, s% TRho_ytop, .false., &
            s% TRho_title, s% TRho_txt_scale, ierr)
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
            s% TRho_logRho_min, s% TRho_logRho_max, &
            s% TRho_logRho_margin, s% TRho_logRho_dlogRho_min, &
            s% TRho_logT_min, s% TRho_logT_max, &
            s% TRho_logT_margin, s% TRho_logT_dlogT_min, &
            s% TRho_step_min, s% TRho_step_max, &
            reverse_xaxis, reverse_yaxis, .false., .false., &
            .false., 0, 0.0, 0.0, 0.0, 0.0, &
            s% show_TRho_annotation1, &
            s% show_TRho_annotation2, &
            s% show_TRho_annotation3, &
            s% TRho_fname, &
            s% TRho_use_decorator, &
            s% TRho_pgstar_decorator, &
            do_degeneracy_line, ierr)
      end subroutine do_TRho_Plot


      subroutine do_degeneracy_line(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (.not. s% show_TRho_degeneracy_line) return
         call pgsave
         call pgsci(clr_Gray)
         call pgsls(Line_Type_Dash)
         call pgline(size(psi4_logT), psi4_logRho, psi4_logT)
         call pgsls(Line_Type_Solid)
         call pgunsa
      end subroutine do_degeneracy_line


      end module pgstar_trho

