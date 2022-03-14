! ***********************************************************************
!
!   Copyright (C) 2013-2019  The MESA Team
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

      module pgstar_summary_burn

      use star_private_def
      use const_def
      use pgstar_support
      use star_pgstar

      implicit none


      contains


      subroutine summary_burn_plot(id, device_id, ierr)
         implicit none
         integer, intent(in) :: id, device_id
         integer, intent(out) :: ierr

         type (star_info), pointer :: s
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         call pgslct(device_id)
         call pgbbuf()
         call pgeras()

         call do_summary_burn_plot(s, id, device_id, &
            s% pg% Summary_Burn_xleft, s% pg% Summary_Burn_xright, &
            s% pg% Summary_Burn_ybot, s% pg% Summary_Burn_ytop, .false., &
            s% pg% Summary_Burn_title, s% pg% Summary_Burn_txt_scale, ierr)

         call pgebuf()

      end subroutine summary_burn_plot


      subroutine do_summary_burn_plot(s, id, device_id, &
            winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, ierr)
         use utils_lib
         use chem_def
         use net_def
         use const_def, only: Msun, Rsun
         type (star_info), pointer :: s
         integer, intent(in) :: id, device_id
         real, intent(in) :: winxmin, winxmax, winymin, winymax
         logical, intent(in) :: subplot
         character (len=*), intent(in) :: title
         real, intent(in) :: txt_scale
         integer, intent(out) :: ierr

         character (len=strlen) :: yname, xaxis_name, str
         logical :: xaxis_reversed
         real, allocatable, dimension(:) :: xvec, yvec, yvec2, yvec3
         real :: xmin, xmax, xleft, xright, dx, windy, dy, &
            ymin, ymax, xaxis_min, xaxis_max, xmargin, &
            legend_xmin, legend_xmax, legend_ymin, legend_ymax
         integer :: lw, lw_sav, grid_min, grid_max, npts, nz
         integer :: docat(num_categories), eps_k(num_categories), num_cat
         real :: xpos_nuc(num_categories)
         real :: xnuc_cat(num_categories), coord
         real(dp) :: eps_max, eps, maxv

         include 'formats'

         ierr = 0
         xaxis_name = s% pg% Summary_Burn_xaxis_name
         xaxis_min = s% pg% Summary_Burn_xmin
         xaxis_max = s% pg% Summary_Burn_xmax
         xaxis_reversed = s% pg% Summary_Burn_xaxis_reversed

         nz = s% nz
         windy = winymax - winymin

         legend_xmin = winxmax - 0.01
         legend_xmax = 0.99
         legend_ymin = winymin
         legend_ymax = winymax

         allocate (xvec(nz), yvec(nz), yvec2(nz), yvec3(nz))

         xmargin = 0
         call set_xaxis_bounds( &
            s, xaxis_name, xaxis_min, xaxis_max, xaxis_reversed, xmargin, &
            xvec, xmin, xmax, xleft, xright, dx, &
            grid_min, grid_max, npts, ierr)

         if (ierr == 0) then
            call pgsave
            call pgsch(txt_scale)
            call plot(ierr)
            call pgunsa
         end if

         deallocate(xvec, yvec, yvec2, yvec3)


         contains


         subroutine plot(ierr)
            use rates_def
            integer, intent(out) :: ierr

            integer :: j, ii, jj, i, cnt, k
            logical, parameter :: dbg = .false.
            real :: ybot, yvec_min, yvec_max

            include 'formats'

            ymax = 1.02
            ymin = 0.0
            ybot = -0.02

            lw = s% pg% pgstar_lw
            call pgqlw(lw_sav)

            call pgsvp(winxmin, winxmax, winymin, winymax)
            if (.not. subplot) then
               call show_model_number_pgstar(s)
               call show_age_pgstar(s)
            end if
            call show_title_pgstar(s, title, s% pg% Summary_Burn_title_shift)

            ! logT
            do k=grid_min,grid_max
               yvec(k) = s% lnT(k)/ln10
            end do
            ymax = maxval(yvec(grid_min:grid_max))
            ymin = minval(yvec(grid_min:grid_max))
            dy = ymax-ymin
            ymax = ymax + 0.1*dy
            ymin = ymin - 0.1*dy

            call pgswin(xleft, xright, ymin+ybot, ymax)

            call pgscf(1)
            call pgsci(1)
            call pgsch(txt_scale)
            call pgbox('BCNST',0.0,0,'CMSTV',0.0,0)
            call pgsci(clr_Goldenrod)
            call pgmtxt('R',5.0,0.5,0.5,'log T')
            call pgslw(lw)
            call pgline(npts, xvec(grid_min:grid_max), yvec(grid_min:grid_max))
            call pgslw(lw_sav)

            do k=grid_min,grid_max
               yvec(k) = safe_log10(s% non_nuc_neu(k) + s% eps_nuc_neu_total(k))
               yvec2(k) = s% lnd(k)/ln10
               yvec3(k) = safe_log10(abs(s% eps_nuc(k)))
            end do

            ymax = maxval(yvec(grid_min:grid_max))
            ymax = max(ymax,maxval(yvec2(grid_min:grid_max)))
            ymax = max(ymax,maxval(yvec3(grid_min:grid_max)))
            if (ymax <= 10) then
               ymax = 11.1
            else
               ymax = ymax*1.1
            end if

            ymin = minval(yvec(grid_min:grid_min))
            ymin = min(ymin,minval(yvec2(grid_min:grid_min)))
            ymin = min(ymin,minval(yvec3(grid_min:grid_min)))
            ymin = max(-6.6, ymin)

            ymax = max(ymax, ymin+1)
            dy = ymax-ymin
            ymin = ymin-dy*0.1

            call pgswin(xleft, xright, ymin+ybot, ymax)
            call pgscf(1)
            call pgsci(1)
            call pgsch(txt_scale)
            call pgbox('',0.0,0,'BNSTV',0.0,0)

            call pgsci(clr_Lilac)
            call pgmtxt('L',4.9,0.5,0.5,'log \gr')
            call pgslw(lw)
            call pgline(npts, xvec(grid_min:grid_max), yvec2(grid_min:grid_max))
            call pgslw(lw_sav)

            call pgsci(clr_LightSteelBlue)
            call pgsch(txt_scale)
            call pgmtxt('L',3.7,0.5,0.5,'log \ge\dnuc\u')
            call pgqlw(lw)
            call pgslw(8)
            call pgline(npts, xvec(grid_min:grid_max), yvec3(grid_min:grid_max))
            call pgslw(lw_sav)

            call pgsci(clr_LightSkyGreen)
            call pgsch(txt_scale)
            call pgmtxt('L',2.5,0.5,0.5,'log \ge\d\gn\u burn+thermal')
            call pgqlw(lw)
            call pgslw(8)
            call pgline(npts, xvec(grid_min:grid_max), yvec(grid_min:grid_max))
            call pgslw(lw_sav)

            ! label peaks of eps_nuc
            ! (stack the labels if they are too close in the x-direction)
            num_cat = 0
            do i = 1, num_categories
               if (i == iphoto) cycle
               maxv = safe_log10(maxval(s% eps_nuc_categories(i,1:nz)))
               if (maxv > -1) then
                  num_cat = num_cat + 1
                  docat(num_cat) = i
                  xnuc_cat(num_cat) = maxv
               end if
            end do
            call pgsci(1)
            call pgsch(txt_scale*0.8)
            do ii = 1, num_cat
               eps_max = -100; i = 0
               do jj = 1, num_cat
                  if (xnuc_cat(jj) < -1) cycle
                  if (xnuc_cat(jj) > eps_max) then
                     eps_max = xnuc_cat(jj)
                     i = jj
                  end if
               end do
               if (i == 0) exit
               if (dbg) write(*,2) 'place ' // category_name(docat(i)), i, eps_max
               xnuc_cat(i) = -1e10 ! mark as done
               eps_max = -100
               eps_k(i) = 0
               do k = 1, nz ! if limit this to grid_min:grid_max, locations jump around too much
                  eps = s% eps_nuc_categories(docat(i),k)
                  if(eps > eps_max) then
                     eps_max = eps
                     eps_k(i) = k
                  end if
               end do
               if(eps_k(i) > 0) then
                  k = eps_k(i)
                  xpos_nuc(i) = xvec(k)
                  if (xpos_nuc(i) < xmin .or. xpos_nuc(i) > xmax) cycle
                  cnt = 0
                  do j = 1, num_cat ! compare location to ones already placed
                     if (j == i) cycle
                     if (xnuc_cat(j) > -1) cycle ! haven't done this one yet
                     if (abs(xpos_nuc(i) - xpos_nuc(j)) < 0.1*dx) then
                        cnt = cnt + 1
                        if (dbg) write(*,*) 'conflicts with ' // category_name(docat(j))
                     end if
                  end do
                  if (cnt < 3) then ! only show 3 max
                     str = category_name(docat(i))
                     if (str(1:5) == 'burn_') then
                        str = str(6:len_trim(str))
                     else if (str == 'tri_alfa') then
                        str = '3a'
                     else if (str == 'c12_c12') then
                        str = 'c+c'
                     else if (str == 'c12_o16') then
                        str = 'c+o'
                     else if (str == 'o16_o16') then
                        str = 'o+o'
                     end if

                     coord = (xpos_nuc(i) - xmin)/(xmax - xmin)
                     if (xaxis_reversed) coord = 1.0 - coord
                     call show_title_label_pgmtxt_pgstar( &
                        s, coord, 0.5, trim(str), 1.03*(cnt+1))
                     !call pgptxt(xpos_nuc(i), ylb, 0.0, 0.5, trim(str))
                  end if
               end if
            end do

            call pgsci(1)
            ierr = 0
            call show_xaxis_name(s,xaxis_name,ierr)
            if (ierr == 0) then ! show mix regions at bottom of plot
               call pgslw(10)
               call show_mix_regions_on_xaxis( &
                  s,ymin+ybot,ymax,grid_min,grid_max,xvec)
            end if

         call show_pgstar_decorator(s%id, s% pg% summary_burn_use_decorator, &
            s% pg% summary_burn_pgstar_decorator, 0, ierr)


         end subroutine plot


      end subroutine do_summary_burn_plot


      end module pgstar_summary_burn

