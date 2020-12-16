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

      module pgstar_l_teff

      use star_private_def
      use const_def
      use pgstar_support

      implicit none


      contains


      subroutine L_Teff_Plot(id, device_id, ierr)
         integer, intent(in) :: id, device_id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         call pgslct(device_id)
         call pgbbuf()
         call pgeras()

         call do_L_Teff_Plot(s, id, device_id, &
            s% L_Teff_xleft, s% L_Teff_xright, &
            s% L_Teff_ybot, s% L_Teff_ytop, .false., &
            s% L_Teff_title, s% L_Teff_txt_scale, ierr)
         if (ierr /= 0) return

         call pgebuf()

      end subroutine L_Teff_Plot


      subroutine do_L_Teff_Plot(s, id, device_id, &
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
            reverse_xaxis = .true., reverse_yaxis = .false.
         ierr = 0
         call do_Hist_Track(s, id, device_id, &
            xleft, xright, ybot, ytop, subplot, title, txt_scale, &
            'effective_T', 'luminosity', &
            'Teff', 'L/L\d\(2281)', &
            s% L_Teff_Teff_min, s% L_Teff_Teff_max, &
            s% L_Teff_Teff_margin, s% L_Teff_dTeff_min, &
            s% L_Teff_L_min, s% L_Teff_L_max, &
            s% L_Teff_L_margin, s% L_Teff_dL_min, &
            s% L_Teff_step_min, s% L_Teff_step_max, &
            reverse_xaxis, reverse_yaxis, .false., .false., &
            s% show_L_Teff_target_box, s% L_Teff_target_n_sigma, &
            s% L_Teff_target_Teff, s% L_Teff_target_L, &
            s% L_Teff_target_Teff_sigma, s% L_Teff_target_L_sigma, &
            s% show_L_Teff_annotation1, &
            s% show_L_Teff_annotation2, &
            s% show_L_Teff_annotation3, &
            s% L_Teff_fname, &
            s% L_Teff_use_decorator, &
            s% L_Teff_pgstar_decorator, &
            null_decorate, ierr)
      end subroutine do_L_Teff_Plot


      end module pgstar_l_teff

