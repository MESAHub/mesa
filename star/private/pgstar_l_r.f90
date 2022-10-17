! ***********************************************************************
!
!   Copyright (C) 2014  The MESA Team
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either rersion 2 of the License,
!   or (at your option) any later rersion.
!
!   You should hare receired a copy of the MESA MANIFESTO along with
!   this software; if not, it is arailable at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without eren the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should hare receired a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

module pgstar_l_r
   
   use star_private_def
   use const_def
   use pgstar_support
   use star_pgstar
   
   implicit none


contains
   
   
   subroutine L_r_Plot(id, device_id, ierr)
      integer, intent(in) :: id, device_id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) return
      
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      
      call do_L_r_Plot(s, id, device_id, &
         s% pg% L_r_xleft, s% pg% L_r_xright, &
         s% pg% L_r_ybot, s% pg% L_r_ytop, .false., &
         s% pg% L_r_title, s% pg% L_r_txt_scale, ierr)
      if (ierr /= 0) return
      
      call pgebuf()
   
   end subroutine L_r_Plot
   
   
   subroutine do_L_r_Plot(s, id, device_id, &
      xleft, xright, ybot, ytop, subplot, &
      title, txt_scale, ierr)
      use pgstar_hist_track, only : null_decorate, do_Hist_Track
      type (star_info), pointer :: s
      integer, intent(in) :: id, device_id
      real, intent(in) :: xleft, xright, ybot, ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      logical, parameter :: &
         reverse_xaxis = .false., reverse_yaxis = .false.
      ierr = 0
      call do_Hist_Track(s, id, device_id, &
         xleft, xright, ybot, ytop, subplot, title, txt_scale, &
         'radius', 'luminosity', &
         'R/R\d\(2281)', 'L/L\d\(2281)', &
         s% pg% L_r_r_min, s% pg% L_r_r_max, &
         s% pg% L_r_r_margin, s% pg% L_r_dr_min, &
         s% pg% L_r_L_min, s% pg% L_r_L_max, &
         s% pg% L_r_L_margin, s% pg% L_r_dL_min, &
         s% pg% L_r_step_min, s% pg% L_r_step_max, &
         reverse_xaxis, reverse_yaxis, .false., .false., &
         s% pg% show_L_r_target_box, s% pg% L_r_target_n_sigma, &
         s% pg% L_r_target_r, s% pg% L_r_target_L, &
         s% pg% L_r_target_r_sigma, s% pg% L_r_target_L_sigma, &
         s% pg% show_L_r_annotation1, &
         s% pg% show_L_r_annotation2, &
         s% pg% show_L_r_annotation3, &
         s% pg% L_r_fname, &
         s% pg% L_r_use_decorator, &
         s% pg% L_r_pgstar_decorator, &
         null_decorate, ierr)
   end subroutine do_L_r_Plot


end module pgstar_l_r

