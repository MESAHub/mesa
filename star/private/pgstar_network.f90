   ! ***********************************************************************
   !
   !   Copyright (C) 2015-2019  The MESA Team
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
      module pgstar_network

      use star_private_def
      use const_def, only: dp
      use pgstar_support
      use star_pgstar
      use pgstar_colors

      implicit none
      private

      public :: network_plot, do_network_plot

      contains

      subroutine network_plot(id, device_id, ierr)
         integer, intent(in) :: id, device_id
         integer, intent(out) :: ierr

         type (star_info), pointer :: s

         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         call pgslct(device_id)
         call pgbbuf()
         call pgeras()

         call do_Network_plot(s, id, device_id, &
            s% pg% Network_xleft, s% pg% Network_xright, &
            s% pg% Network_ybot, s% pg% Network_ytop, .false., &
            s% pg% Network_title, s% pg% Network_txt_scale, ierr)

         call pgebuf()

      end subroutine network_plot


      subroutine do_network_plot(s, id, device_id, &
            winxmin, winxmax, winymin, winymax, subplot, &
            title, txt_scale, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, device_id
         real, intent(in) :: winxmin, winxmax, winymin, winymax
         logical, intent(in) :: subplot
         character (len=*), intent(in) :: title
         real, intent(in) :: txt_scale
         integer, intent(out) :: ierr
         call do_network_panel(s, id, device_id, &
            winxmin, winxmax, winymin, winymax, subplot, &
            title, txt_scale, ierr)
      end subroutine do_network_plot


      subroutine do_network_panel(s, id, device_id, &
            winxmin, winxmax, winymin, winymax, subplot, title, txt_scale, &
            ierr)
         use utils_lib
         use chem_def
         use net_def
         use const_def, only: Msun

         type (star_info), pointer :: s
         integer, intent(in) :: id, device_id
         real, intent(in) :: &
            winxmin, winxmax, winymin, winymax
         character (len=*), intent(in) :: title
         real, intent(in) :: txt_scale
         logical, intent(in) :: subplot
         integer, intent(out) :: ierr

         real :: xleft, xright, chScale, xmargin

         include 'formats'
         ierr = 0

         chScale = txt_scale

         xmargin = 0
         call plot(ierr)


         contains

         subroutine plot(ierr)
            use chem_def
            integer, intent(out) :: ierr

            integer :: i, j

            integer :: z,n,zmax,zmin,nmin,nmax
            integer :: clr,mid_map
            real :: abun,xhigh,xlow
            real :: ymin,ymax,log10_min_abun,log10_max_abun
            real,parameter :: pad=2.5,step=0.5

            include 'formats'
            ierr = 0

            call pgsave
            call pgsch(txt_scale)
            call pgsvp(winxmin, winxmax, winymin, winymax)

            zmax=0
            nmax=0
            zmin=HUGE(zmin)
            nmin=HUGE(nmin)


            log10_min_abun=s% pg% Network_log_mass_frac_min
            log10_max_abun=s% pg% Network_log_mass_frac_max

            do i=1,s%species

               Z=chem_isos%Z(s%chem_id(i))
               N=chem_isos%N(s%chem_id(i))

               zmax=max(Z,zmax)
               nmax=max(n,nmax)

               zmin=min(Z,zmin)
               nmin=min(n,nmin)

            end do

            if (s% pg% network_zmax > -100) then
               ymax = s% pg% network_zmax
            else
               ymax = zmax
            end if

            if (s% pg% network_zmin > -100) then
               ymin = s% pg% network_zmin
            else
               ymin = zmin
            end if

            if (s% pg% network_nmax > -100) then
               xright = s% pg% network_nmax
            else
               xright= nmax
            end if

            if (s% pg% network_nmin > -100) then
               xleft = s% pg% network_nmin
            else
               xleft = nmin
            end if

            !Set xaxis and yaxis bounds
            call pgswin(xleft-5,xright+pad,ymin-pad,ymax+pad)
            !Create a box with ticks
            call show_box_pgstar(s,'BCNST','BCNSTV')
            !Labels
            call show_xaxis_name(s,'N',ierr)
            call show_left_yaxis_label_pgstar(s,'Z',-1.5)

            if (.not. subplot) then
               call show_model_number_pgstar(s)
               call show_age_pgstar(s)
            end if
            call show_title_pgstar(s, title)

            mid_map = colormap_length/2
            do i=1,s%species

               Z=chem_isos%Z(s%chem_id(i))
               N=chem_isos%N(s%chem_id(i))
               abun=(dot_product(s%xa(i,1:s%nz),s%dm(1:s%nz))/msun)/&
                     ((s%star_mass)-(s%m_center/msun))

               abun=safe_log10(dble(abun))

               if(z<ymin .or. z>ymax .or. n<xleft .or.n>xright)CYCLE

               if (s% pg% Network_show_element_names) THEN
                  call pgsci(clr_Foreground)
                  call pgtext(xleft-3.5,z*1.0-0.25,el_name(Z))
               end if

               !Plot colored dots for mass fractions
               if(s% pg% Network_show_mass_fraction) then
                  if(abun>log10_min_abun .and. abun < log10_max_abun)THEN
                     do j=mid_map,colormap_length
                        xlow=log10_min_abun+(j-mid_map)*(log10_max_abun-log10_min_abun)/(colormap_length-mid_map)
                        xhigh=log10_min_abun+(j-mid_map+1)*(log10_max_abun-log10_min_abun)/(colormap_length-mid_map)
                        if(abun>=xlow .and. abun<xhigh)THEN
                           clr = colormap_offset + (colormap_length-(j-mid_map))
                           call pgsci(clr)
                        end if
                     end do

                     call pgrect(n-step,n+step,z-step,z+step)
                  end if
               end if

               !Plot box centered on the (N,Z)
               call pgsci(clr_Foreground)
               call pgline(5,[n-step,n+step,n+step,n-step,n-step],[z-step,z-step,z+step,z+step,z-step])
            end do

            call pgunsa

            if(s% pg% network_show_colorbar)then
               call network_colorbar_legend(winxmin, winxmax, winymin, winymax,log10_min_abun,log10_max_abun)
            end if

            call show_pgstar_decorator(s%id,s% pg% network_use_decorator,s% pg% network_pgstar_decorator, 0, ierr)


         end subroutine plot


      end subroutine do_network_panel


      subroutine network_colorbar_legend(winxmin, winxmax, winymin, winymax,abun_min,abun_max)
         real,intent(in) :: winxmin, winxmax, winymin, winymax,abun_min,abun_max
         real :: legend_xmin,legend_xmax,legend_ymin,legend_ymax
         real :: xmin,xmax,ymin,ymax
         real :: dx, dyline, xpts(2),yt,yb,text
         character(len=16) :: str

         integer :: i,j,clr,mid_map,num_cms

         call PGQWIN(xmin, xmax, ymin, ymax)

         legend_xmin = winxmax - 0.01
         legend_xmax = 0.99
         legend_ymin = winymin
         legend_ymax = winymax

         mid_map = colormap_length/2
         num_cms=colormap_length-mid_map
         dyline = (ymax-ymin)/num_cms
         dx = 0.1

         xpts(1) = 2.0*dx
         xpts(2) = xpts(1) + 2.0*dx

         call pgsave
         call pgsvp(legend_xmin, legend_xmax, legend_ymin, legend_ymax)
         call pgswin(0.0, 1.0, ymin, ymax)
         do j=mid_map,colormap_length
            i=j-mid_map
            clr = colormap_offset + (colormap_length-i+1)
            call pgsci(clr)
            yt = ymin + (i)*dyline
            yb = ymin + (i-1)*dyline

            call pgrect(xpts(1),xpts(2),yb,yt)
         end do

         call pgsci(clr_Foreground)
         do j=1,5
            text=abun_min+(j-1)*(abun_max-abun_min)/4.0
            write(str,'(F8.3)') text
            call pgptxt(xpts(2) + 0.025, ymin+(j-1)*(ymax-ymin)/4.0, 0.0, 0.0, trim(str))
         end do

         call pgunsa

      end subroutine network_colorbar_legend

      end module pgstar_network
