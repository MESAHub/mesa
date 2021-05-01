   ! ***********************************************************************
   !
   !   Copyright (C) 2015-2019  Bill Paxton & The MESA Team
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

      module pgstar_network

      use star_private_def
      use const_def
      use pgstar_support

      implicit none


      contains

      subroutine network_plot(id, device_id, ierr)
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

         call do_Network_plot(s, id, device_id, &
            s% Network_xleft, s% Network_xright, &
            s% Network_ybot, s% Network_ytop, .false., &
            s% Network_title, s% Network_txt_scale, ierr)

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
         implicit none

         type (star_info), pointer :: s
         integer, intent(in) :: id, device_id
         real, intent(in) :: &
            winxmin, winxmax, winymin, winymax
         character (len=*), intent(in) :: title
         real, intent(in) :: txt_scale
         logical, intent(in) :: subplot
         integer, intent(out) :: ierr

         character (len=strlen) :: str
         real :: xmin, xmax, xleft, xright, dx, dylbl, chScale, windy, xmargin, &
            ymin, ymax
         integer :: lw, lw_sav, grid_min, grid_max, npts, i, nz

         include 'formats'
         ierr = 0

         chScale = txt_scale

         xmargin = 0
         call plot(ierr)


         contains

         subroutine plot(ierr)
            use chem_def
            integer, intent(out) :: ierr

            integer :: lw, lw_sav, k,i,j
            real :: ybot, eps

            integer :: z,n,zmax,zmin,nmin,nmax
            integer :: base_z,base_n,clr,mid_map
            real :: abun,xhigh,xlow
            real :: ymin,ymax,r,g,b,log10_min_abun,log10_max_abun
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


            log10_min_abun=s%Network_log_mass_frac_min
            log10_max_abun=s%Network_log_mass_frac_max

            do i=1,s%species

               Z=chem_isos%Z(s%chem_id(i))
               N=chem_isos%N(s%chem_id(i))

               zmax=max(Z,zmax)
               nmax=max(n,nmax)

               zmin=min(Z,zmin)
               nmin=min(n,nmin)

            end do

            if (s% network_zmax > -100) then
               ymax = s% network_zmax
            else
               ymax = zmax
            end if

            if (s% network_zmin > -100) then
               ymin = s% network_zmin
            else
               ymin = zmin
            end if

            if (s% network_nmax > -100) then
               xright = s% network_nmax
            else
               xright= nmax
            end if

            if (s% network_nmin > -100) then
               xleft = s% network_nmin
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

            mid_map = colormap_size/2
            do i=1,s%species

               Z=chem_isos%Z(s%chem_id(i))
               N=chem_isos%N(s%chem_id(i))
               abun=(dot_product(s%xa(i,1:s%nz),s%dm(1:s%nz))/msun)/&
                     ((s%mstar/Msun)-(s%m_center/msun))

               abun=safe_log10(dble(abun))

               if(z.lt.ymin .or. z.gt.ymax .or. n.lt.xleft .or.n.gt.xright)CYCLE

               if (s% Network_show_element_names) THEN
                  call pgsci(1)
                  call pgtext(xleft-3.5,z*1.0-0.25,el_name(Z))
               end if

               !Plot colored dots for mass fractions
               if(s% Network_show_mass_fraction) then
                  if(abun>log10_min_abun .and. abun < log10_max_abun)THEN
                     do j=mid_map,colormap_size
                        xlow=log10_min_abun+(j-mid_map)*(log10_max_abun-log10_min_abun)/(colormap_size-mid_map)
                        xhigh=log10_min_abun+(j-mid_map+1)*(log10_max_abun-log10_min_abun)/(colormap_size-mid_map)
                        if(abun>=xlow .and. abun<xhigh)THEN
                           clr = colormap_offset + (colormap_size-(j-mid_map))
                           call pgsci(clr)
                        end if
                     end do

                     call pgrect(n-step,n+step,z-step,z+step)
                  end if
               end if

               !Plot box centered on the (N,Z)
               call pgsci(1)
               call pgline(5,(/n-step,n+step,n+step,n-step,n-step/),(/z-step,z-step,z+step,z+step,z-step/))
            end do

            call pgunsa

            if(s% network_show_colorbar)then
               call network_colorbar_legend(winxmin, winxmax, winymin, winymax,log10_min_abun,log10_max_abun)
            end if
            
            call show_pgstar_decorator(s%id,s% network_use_decorator,s% network_pgstar_decorator, 0, ierr)


         end subroutine plot


      end subroutine do_network_panel
      
      
      subroutine network_colorbar_legend(winxmin, winxmax, winymin, winymax,abun_min,abun_max)
         real,intent(in) :: winxmin, winxmax, winymin, winymax,abun_min,abun_max
         real :: legend_xmin,legend_xmax,legend_ymin,legend_ymax
         real :: xmin,xmax,ymin,ymax
         real :: ymx, dx, dyline, ypos, xpts(2),yt,yb,text
         character(len=16) :: str
         
         integer :: i,j,clr,mid_map,num_cms

         call PGQWIN(xmin, xmax, ymin, ymax)

         legend_xmin = winxmax - 0.01
         legend_xmax = 0.99
         legend_ymin = winymin
         legend_ymax = winymax

         mid_map = colormap_size/2
         num_cms=colormap_size-mid_map
         dyline = (ymax-ymin)/num_cms
         dx = 0.1
         
         xpts(1) = 2.0*dx
         xpts(2) = xpts(1) + 2.0*dx

         call pgsave
         call pgsvp(legend_xmin, legend_xmax, legend_ymin, legend_ymax)
         call pgswin(0.0, 1.0, ymin, ymax)
         do j=mid_map,colormap_size
            i=j-mid_map
            clr = colormap_offset + (colormap_size-i+1)
            call pgsci(clr)
            yt = ymin + (i)*dyline
            yb = ymin + (i-1)*dyline
         
            call pgrect(xpts(1),xpts(2),yb,yt)
         end do

         call pgsci(1)
         do j=1,5
            text=abun_min+(j-1)*(abun_max-abun_min)/4.0
            write(str,'(F8.3)') text
            call pgptxt(xpts(2) + 0.025, ymin+(j-1)*(ymax-ymin)/4.0, 0.0, 0.0, trim(str))
         end do
         
         call pgunsa

      end subroutine network_colorbar_legend

      end module pgstar_network
