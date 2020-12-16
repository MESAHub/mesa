! ***********************************************************************
!
!   Copyright (C) 2013  Bill Paxton
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

      module pgstar_dPg_dnu

      use star_private_def
      use const_def
      use pgstar_support

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
            s% dPg_dnu_xleft, s% dPg_dnu_xright, &
            s% dPg_dnu_ybot, s% dPg_dnu_ytop, .false., &
            s% dPg_dnu_title, s% dPg_dnu_txt_scale, ierr)

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
            s% dPg_dnu_delta_nu_min, s% dPg_dnu_delta_nu_max, &
            s% dPg_dnu_delta_nu_margin, s% dPg_dnu_d_delta_nu_min, &
            s% dPg_dnu_delta_Pg_min, s% dPg_dnu_delta_Pg_max, &
            s% dPg_dnu_delta_Pg_margin, s% dPg_dnu_d_delta_Pg_min, &
            s% dPg_dnu_step_min, s% dPg_dnu_step_max, &
            reverse_xaxis, reverse_yaxis, .false., .false., &
            .false., 0, 0.0, 0.0, 0.0, 0.0, &
            s% show_dPg_dnu_annotation1, &
            s% show_dPg_dnu_annotation2, &
            s% show_dPg_dnu_annotation3, &
            s% dPg_dnu_fname, &
            s% dPg_dnu_use_decorator, s% dPg_dnu_pgstar_decorator, &
            null_decorate, ierr)
      end subroutine do_dPg_dnu_Plot


      end module pgstar_dPg_dnu






