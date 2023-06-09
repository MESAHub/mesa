! ***********************************************************************
!
!   Copyright (C) 2013  The MESA Team
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

module pgstar_hist_track

   use star_private_def
   use const_def
   use pgstar_support
   use star_pgstar

   implicit none


contains


   subroutine History_Track_plot(id, device_id, array_ix, ierr)
      integer, intent(in) :: id, device_id, array_ix
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call pgslct(device_id)
      call pgbbuf()
      call pgeras()
      call do_History_Track_plot(s, id, device_id, array_ix, &
         s% pg% History_Track_xleft(array_ix), s% pg% History_Track_xright(array_ix), &
         s% pg% History_Track_ybot(array_ix), s% pg% History_Track_ytop(array_ix), .false., &
         s% pg% History_Track_title(array_ix), s% pg% History_Track_txt_scale(array_ix), ierr)
      if (ierr /= 0) return
      call pgebuf()
   end subroutine History_Track_plot


   subroutine do_History_Track_plot(s, id, device_id, array_ix, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, ierr)
      type (star_info), pointer :: s
      integer, intent(in) :: id, device_id, array_ix
      real, intent(in) :: vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale
      logical, intent(in) :: subplot
      character (len = *), intent(in) :: title
      integer, intent(out) :: ierr
      call do_Hist_Track(s, id, device_id, array_ix, &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
         s% pg% History_Track_xname(array_ix), &
         s% pg% History_Track_yname(array_ix), &
         s% pg% History_Track_xaxis_label(array_ix), &
         s% pg% History_Track_yaxis_label(array_ix), &
         s% pg% History_Track_xmin(array_ix), &
         s% pg% History_Track_xmax(array_ix), &
         s% pg% History_Track_xmargin(array_ix), &
         s% pg% History_Track_dxmin(array_ix), &
         s% pg% History_Track_ymin(array_ix), &
         s% pg% History_Track_ymax(array_ix), &
         s% pg% History_Track_ymargin(array_ix), &
         s% pg% History_Track_dymin(array_ix), &
         s% pg% History_Track_step_min(array_ix), &
         s% pg% History_Track_step_max(array_ix), &
         s% pg% History_Track_reverse_xaxis(array_ix), &
         s% pg% History_Track_reverse_yaxis(array_ix), &
         s% pg% History_Track_log_xaxis(array_ix), &
         s% pg% History_Track_log_yaxis(array_ix), &
         s% pg% show_History_Track_target_box(array_ix), &
         s% pg% History_Track_n_sigma(array_ix), &
         s% pg% History_Track_xtarget(array_ix), &
         s% pg% History_Track_ytarget(array_ix), &
         s% pg% History_Track_xsigma(array_ix), &
         s% pg% History_Track_ysigma(array_ix), &
         s% pg% show_History_Track_annotation1(array_ix), &
         s% pg% show_History_Track_annotation2(array_ix), &
         s% pg% show_History_Track_annotation3(array_ix), &
         s% pg% History_Track_fname(array_ix), &
         s% pg% History_Track_use_decorator(array_ix), &
         s% pg% History_Track_pgstar_decorator(array_ix), &
         null_decorate, ierr)
   end subroutine do_History_Track_plot


   subroutine null_decorate(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      ierr = 0
   end subroutine null_decorate


   subroutine do_Hist_Track(s, id, device_id, array_ix, &
      vp_xleft, vp_xright, vp_ybot, vp_ytop, subplot, title, txt_scale, &
      xname, yname, xaxis_label, yaxis_label, &
      given_xmin, given_xmax, xmargin, dxmin, &
      given_ymin, given_ymax, ymargin, dymin, &
      step_min_in, step_max_in, &
      reverse_xaxis, reverse_yaxis, &
      log_xaxis, log_yaxis, &
      show_target_box, n_sigma, &
      xtarget, ytarget, xsigma, ysigma, &
      show_annotation1, show_annotation2, show_annotation3, &
      fname, &
      use_decorator, pgstar_decorator, &
      decorate, ierr)
      use utils_lib

      type (star_info), pointer :: s
      integer, intent(in) :: &
         id, device_id, step_min_in, step_max_in, n_sigma, array_ix
      real, intent(in) :: &
         vp_xleft, vp_xright, vp_ybot, vp_ytop, txt_scale, &
         xtarget, ytarget, xsigma, ysigma, &
         given_xmin, given_xmax, xmargin, dxmin, &
         given_ymin, given_ymax, ymargin, dymin
      character (len = *), intent(in) :: &
         title, xname, yname, xaxis_label, yaxis_label, fname
      logical, intent(in) :: subplot, &
         reverse_xaxis, reverse_yaxis, log_xaxis, log_yaxis, show_target_box, &
         show_annotation1, show_annotation2, show_annotation3, use_decorator
      interface
         subroutine decorate(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
         end subroutine decorate
      end interface
      integer, intent(out) :: ierr
      procedure(pgstar_decorator_interface), pointer :: pgstar_decorator

      real :: xmin, xmax, ymin, ymax, xleft, xright, ybot, ytop
      integer :: i, j, j_min, j_max, step_min, step_max
      real :: dx, dy, xplus, xminus, yplus, yminus
      real, dimension(:), allocatable :: xvec, yvec
      character (len = strlen) :: str
      integer :: k, n
      integer :: ix, iy
      integer :: file_data_len
      real, allocatable, dimension(:) :: file_data_xvec, file_data_yvec

      logical, parameter :: dbg = .false.

      include 'formats'

      ierr = 0

      call integer_dict_lookup(s% history_names_dict, xname, ix, ierr)
      if (ierr /= 0) ix = -1
      if (ix <= 0) then
         write(*, '(A)')
         write(*, *) 'ERROR: failed to find ' // &
            trim(xname) // ' in history data'
         write(*, '(A)')
         ierr = -1
      end if

      call integer_dict_lookup(s% history_names_dict, yname, iy, ierr)
      if (ierr /= 0) iy = -1
      if (iy <= 0) then
         write(*, '(A)')
         write(*, *) 'ERROR: failed to find ' // &
            trim(yname) // ' in history data'
         write(*, '(A)')
         ierr = -1
      end if
      if (ierr /= 0) return
      step_min = max(step_min_in, 1)
      if (step_max_in >= 0) then
         step_max = min(step_max_in, s% model_number)
      else
         step_max = s% model_number
      end if
      n = count_hist_points(s, step_min, step_max)
      allocate(xvec(n), yvec(n), stat = ierr)
      if (ierr /= 0) then
         write(*, *) 'allocate failed for PGSTAR'
         return
      end if

      call get_hist_points(s, step_min, step_max, n, ix, xvec, ierr)
      if (ierr /= 0) then
         write(*, *) 'pgstar do_Hist_Track get_hist_points failed ' // trim(xname)
         ierr = 0
         !return
      end if
      call get_hist_points(s, step_min, step_max, n, iy, yvec, ierr)
      if (ierr /= 0) then
         write(*, *) 'pgstar do_Hist_Track get_hist_points failed ' // trim(yname)
         ierr = 0
         !return
      end if

      if (log_xaxis) then
         do k = 1, n
            xvec(k) = log10(max(tiny(xvec(k)), abs(xvec(k))))
         end do
      end if

      if (log_yaxis) then
         do k = 1, n
            yvec(k) = log10(max(tiny(yvec(k)), abs(yvec(k))))
         end do
      end if

      call set_xleft_xright(&
         n, xvec, given_xmin, given_xmax, xmargin, &
         reverse_xaxis, dxmin, xleft, xright)

      call set_ytop_ybot(&
         n, yvec, given_ymin, given_ymax, -101.0, ymargin, &
         reverse_yaxis, dymin, ybot, ytop)

      call pgsave
      call pgsch(txt_scale)
      call pgsvp(vp_xleft, vp_xright, vp_ybot, vp_ytop)
      call pgswin(xleft, xright, ybot, ytop)
      call pgscf(1)
      call pgsci(1)
      call show_box_pgstar(s, 'BCNST1', 'BCNSTV1')

      if (log_xaxis) then
         call show_xaxis_label_pgstar(s, 'log ' // xaxis_label)
      else
         call show_xaxis_label_pgstar(s, xaxis_label)
      end if

      if (log_yaxis) then
         call show_left_yaxis_label_pgstar(s, 'log ' // yaxis_label)
      else
         call show_left_yaxis_label_pgstar(s, yaxis_label)
      end if

      if (.not. subplot) then
         call show_model_number_pgstar(s)
         call show_age_pgstar(s)
      end if
      call show_title_pgstar(s, title)

      call pgslw(s% pg% pgstar_lw)

      call show_file_track

      if (show_target_box) then
         call pgsci(clr_Silver)
         if (n_sigma >= 0) then
            j_min = n_sigma
            j_max = n_sigma
         else
            j_min = 1
            j_max = -n_sigma
         end if
         do j = j_min, j_max
            dx = xsigma * j
            xplus = xtarget + dx
            xminus = xtarget - dx
            if (log_xaxis) then
               xplus = log10(max(tiny(xplus), xplus))
               xminus = log10(max(tiny(xminus), xminus))
            end if
            dy = ysigma * j
            yplus = ytarget + dy
            yminus = ytarget - dy
            if (log_yaxis) then
               yplus = log10(max(tiny(yplus), yplus))
               yminus = log10(max(tiny(yminus), yminus))
            end if
            call pgmove(xminus, yminus)
            call pgdraw(xplus, yminus)
            call pgdraw(xplus, yplus)
            call pgdraw(xminus, yplus)
            call pgdraw(xminus, yminus)
         end do
      end if

      call pgsci(clr_Teal)
      call pgline(n, xvec, yvec)
      call pgsci(clr_Crimson)
      call pgsch(2.8 * txt_scale)
      call pgpt1(xvec(n), yvec(n), 0902)

      call show_annotations(s, &
         show_annotation1, &
         show_annotation2, &
         show_annotation3)

      call decorate(s% id, ierr)

      call pgunsa

      call show_pgstar_decorator(s%id, use_decorator, pgstar_decorator, 0, ierr)

      deallocate(xvec, yvec)

   contains


      subroutine show_file_track
         integer :: k
         if (len_trim(fname) == 0) return
         if (.not. read_values_from_file(fname, &
            file_data_xvec, file_data_yvec, file_data_len)) then
            write(*, *) &
               'bad filename for History tracks plot ' // trim(fname)
            return
         end if
         if (log_xaxis) then
            do k = 1, file_data_len
               file_data_xvec(k) = log10(max(tiny(file_data_xvec(k)), abs(file_data_xvec(k))))
            end do
         end if
         if (log_yaxis) then
            do k = 1, file_data_len
               file_data_yvec(k) = log10(max(tiny(file_data_yvec(k)), abs(file_data_yvec(k))))
            end do
         end if
         call pgsci(clr_Goldenrod)
         call pgline(file_data_len, file_data_xvec, file_data_yvec)
         deallocate(file_data_xvec, file_data_yvec)
      end subroutine show_file_track


   end subroutine do_Hist_Track


end module pgstar_hist_track

