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

      module pgstar_tmaxrho

      use star_private_def
      use const_def
      use pgstar_support
      use star_pgstar

      implicit none


      contains


      subroutine TmaxRho_Plot(id, device_id, array_ix, ierr)
         integer, intent(in) :: id, device_id, array_ix
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         call pgslct(device_id)
         call pgbbuf()
         call pgeras()

         call do_TmaxRho_Plot(s, id, device_id, &
            s% pg% TmaxRho_xleft, s% pg% TmaxRho_xright, &
            s% pg% TmaxRho_ybot, s% pg% TmaxRho_ytop, .false., &
            s% pg% TmaxRho_title, s% pg% TmaxRho_txt_scale, ierr)
         if (ierr /= 0) return

         call pgebuf()

      end subroutine TmaxRho_Plot


      subroutine do_TmaxRho_Plot(s, id, device_id, &
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
            'log_center_Rho', 'log_max_T', &
            'log Rho\dc\u (g cm\u-3\d)', 'log T\dmax\u (K)', &
            s% pg% TmaxRho_logRho_min, s% pg% TmaxRho_logRho_max, &
            s% pg% TmaxRho_logRho_margin, s% pg% TmaxRho_logRho_dlogRho_min, &
            s% pg% TmaxRho_logT_min, s% pg% TmaxRho_logT_max, &
            s% pg% TmaxRho_logT_margin, s% pg% TmaxRho_logT_dlogT_min, &
            s% pg% TmaxRho_step_min, s% pg% TmaxRho_step_max, &
            reverse_xaxis, reverse_yaxis, .false., .false., &
            .false., 0, 0.0, 0.0, 0.0, 0.0, &
            s% pg% show_TmaxRho_annotation1, &
            s% pg% show_TmaxRho_annotation2, &
            s% pg% show_TmaxRho_annotation3, &
            s% pg% TmaxRho_fname, &
            s% pg% TmaxRho_use_decorator, &
            s% pg% TmaxRho_pgstar_decorator, &
            do_degeneracy_line, ierr)
      end subroutine do_TmaxRho_Plot


      subroutine do_degeneracy_line(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (.not. s% pg% show_TmaxRho_degeneracy_line) return
         call pgsave
         call pgsci(clr_Gray)
         call pgsls(Line_Type_Dash)
         call pgline(size(psi4_logT), psi4_logRho, psi4_logT)
         call pgsls(Line_Type_Solid)
         call pgunsa
      end subroutine do_degeneracy_line


      end module pgstar_tmaxrho

