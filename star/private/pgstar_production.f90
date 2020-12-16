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

      module pgstar_production

      use star_private_def
      use const_def
      use pgstar_support

      implicit none


      contains


      subroutine production_plot(id, device_id, ierr)
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

         call do_Production_plot(s, id, device_id, &
            s% Production_xleft, s% Production_xright, &
            s% Production_ybot, s% Production_ytop, .false., &
            s% Production_title, s% Production_txt_scale, ierr)

         call pgebuf()

      end subroutine production_plot


      subroutine do_production_plot(s, id, device_id, &
            winxmin, winxmax, winymin, winymax, subplot, &
            title, txt_scale, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, device_id
         real, intent(in) :: winxmin, winxmax, winymin, winymax
         logical, intent(in) :: subplot
         character (len=*), intent(in) :: title
         real, intent(in) :: txt_scale
         integer, intent(out) :: ierr
         call do_production_panel(s, id, device_id, &
            winxmin, winxmax, winymin, winymax, subplot, &
            title, txt_scale, ierr)
      end subroutine do_production_plot


      subroutine do_production_panel(s, id, device_id, &
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
         integer, parameter :: num_colors = 14
         integer :: colors(num_colors)

         include 'formats'
         ierr = 0

         colors(:) = (/ &
               clr_Gold, clr_LightSkyBlue, clr_Crimson, clr_Goldenrod, clr_MediumSlateBlue, &
               clr_Coral, clr_LightSkyGreen, clr_DarkGray, clr_Lilac, &
               clr_Tan, clr_IndianRed, clr_Teal, clr_Silver, clr_BrightBlue /)

         chScale = txt_scale

         xmargin = 0
         call plot(ierr)


         contains

         subroutine plot(ierr)
            use chem_def
            use chem_lib
            use adjust_xyz, only: get_xa_for_standard_metals
            integer, intent(out) :: ierr

            integer :: lw, lw_sav, k,i,j
            real :: ybot, eps

            integer :: amin,amax,z,n,a,plot_a,zmin,zmax
            integer :: min_zone,max_zone,alternate,skip_cnt
            real :: xhigh,xlow,extra_pad
            real :: min_mass,max_mass,yloc
            real,parameter :: point_size=0.1
            real :: ymin,ymax,r,g,b,log10_min_abun,log10_max_abun
            real,parameter :: pad=1.0
            real :: last_x,last_y,log_sa
            logical :: z_in_use
            real(dp),dimension(1:solsiz) :: scaled_abun,scaled_abun_init
            real(dp),dimension(:),allocatable :: init_comp,abun

            real(dp) :: initial_z,initial_y,initial_h1,initial_h2
            real(dp) :: initial_he3,initial_he4,xsol_he3,xsol_he4,la,lac

            include 'formats.inc'
            ierr = 0

            call pgsave
            call pgsch(txt_scale)
            call pgsvp(winxmin, winxmax, winymin, winymax)
            
            amax=0
            amin=0.0
            zmax=0
            zmin=HUGE(zmin)
            ymax=-HUGE(ymax)
            ymin=HUGE(ymin)


            if (s% production_min_mass > 0.0) then
               min_mass = s% production_min_mass
            else
               min_mass = 0.d0
            end if

            if (s% production_max_mass > -100) then
               max_mass = s% production_max_mass
            else
               max_mass = s%mstar/msun
            end if

            min_zone=1
            max_zone=s%nz

            do i=1,s%nz
               if(s%m(i).le.max_mass) then
                  min_zone=i-1
                  exit
               end if
            end do

            do i=min_zone,s%nz
               if(s%m(i).le.min_mass) then
                  max_zone=i-1
                  exit
               end if
            end do

            allocate(abun(1:s%species),init_comp(1:s%species))
            abun=0.d0

            !Get initial composotion

            !Stolen from star/private/create_initial_model
            initial_z = s% initial_z
            initial_y = s% initial_y
            if (initial_y < 0) initial_y = max(0d0, min(1d0, 0.24d0 + 2*initial_z))
            initial_h1 = max(0d0, min(1d0, 1d0 - (initial_z + initial_y)))
            initial_h2 = chem_Xsol('h2')
            xsol_he3 = chem_Xsol('he3')
            xsol_he4 = chem_Xsol('he4')
            initial_he3 = initial_y*xsol_he3/(xsol_he3 + xsol_he4)
            initial_he4 = initial_y*xsol_he4/(xsol_he3 + xsol_he4)
            !

            call get_xa_for_standard_metals( &
               s, s%species, s%chem_id, s%net_iso,&
               initial_h1,initial_h2,initial_he3,initial_he4, &
               s% job% initial_zfracs, s% job% dump_missing_metals_into_heaviest, init_comp, ierr)

            !compute abundances
            do i=1,s%species

               Z=chem_isos%Z(s%chem_id(i))
               N=chem_isos%N(s%chem_id(i))

               zmin=min(Z,zmin)
               zmax=max(Z,zmax)

               abun(i)=(dot_product(s%xa(i,min_zone:max_zone),s%dm(min_zone:max_zone))/msun)/&
                     ((max_mass-min_mass)-(s%m_center/msun))

            end do

            !Get stable isotope abundances
            call get_stable_mass_frac(s%chem_id,s%species,dble(abun),scaled_abun) 
            call get_stable_mass_frac(s%chem_id,s%species,init_comp,scaled_abun_init)           

            do i=1,solsiz
               
               la=safe_log10(scaled_abun(i))
               lac=safe_log10(scaled_abun_init(i))
   
               !Remove low abundance isotopes, low in star and low in solar can lead to large production factor
               if(la .lt.s%production_min_mass_frac .or. lac .lt.s%production_min_mass_frac) then
                  scaled_abun(i)=-HUGE(ymin)
               else
                  scaled_abun(i)=real(la-lac)
                  ymax=max(ymax,real(scaled_abun(i),kind=kind(ymax)))
                  ymin=min(ymin,real(scaled_abun(i),kind=kind(ymin)))
               end if
      
               if(zmax==izsol(i)) amax=int(iasol(i))

            end do

            if (s% production_ymax > -100) then
               ymax = s% production_ymax
            else
               ymax = ymax+0.01
            end if

            if (s% production_ymin > -100) then
               ymin = s% production_ymin
            else
               ymin = ymin
            end if

            if (s% production_amax > -100) then
               xright = s% production_amax
            else
               xright= amax
            end if

            if (s% production_amin > -100) then
               xleft = s% production_amin
            else
               xleft = amin
            end if

            if (s% Production_show_element_names) THEN
               extra_pad=2.5
            else
               extra_pad=1.0
            end if

            !Set xaxis and yaxis bounds
            call pgswin(xleft-1,xright+1,ymin-pad,ymax+extra_pad*pad)
            !Create a box with ticks
            call show_box_pgstar(s,'BCNST','BCNSTV')
            !Labels
            call show_xaxis_name(s,'A',ierr)
            call show_left_yaxis_label_pgstar(s,'Log Production Factor',0.0)

            if (.not. subplot) then
               call show_model_number_pgstar(s)
               call show_age_pgstar(s)
            end if
            call show_title_pgstar(s, title)           

            alternate=-1
            i=1
            outer: do 
               if(i.gt.solsiz) exit outer
      
               ! Z is greater than zmax
               if(izsol(i).gt.zmax) exit outer

               !Sets color
               call set_line_style(izsol(i))

               !Shows element name, alternates between two levels to spread them out
               if (s% Production_show_element_names &
                  .and.(iasol(i).le.xright).and.(iasol(i).ge.xleft)) THEN
                     yloc=(ymax*1.0)+abs(0.75+(alternate)/2.0)
                     call pgtext(iasol(i)*1.0,yloc,el_name(izsol(i)))
                     alternate=alternate*(-1.0)
               end if

               !When we have more than one isotope per element we draw a line, this marks the last
               ! x cordinate we "saw" for an element
               last_x=-HUGE(last_x)
               last_y=-HUGE(last_y)

               inner: do j=i,solsiz
               
                  if(izsol(j).eq.izsol(i))then
                     if((scaled_abun(j).ge. ymin) .and. (scaled_abun(j) .le. ymax)&
                        .and.(iasol(j).le.xright).and.(iasol(j).ge.xleft)) then
                        a=iasol(j)

                        !Draw point at values
                        call PGCIRC(A*1.0,real(scaled_abun(j)),point_size)
                        !Not the first isotope of an element
                        if(last_x>0)then
                           !Then draw a line between isotopes of same element
                           call pgline(2,(/last_x,A*1.0/),(/last_y,real(scaled_abun(j)*1.0)/))
                        end if

                        !Save last x,y pair we saw
                        last_x=A*1.0
                        last_y=scaled_abun(j)               
                     end if
                  else
                     exit inner
                  end if

               end do inner
               !Jump to next element
               i=j
            end do outer

            call pgunsa
            deallocate(abun,init_comp)
            
            call show_pgstar_decorator(s%id,s% production_use_decorator,&
                  s% production_pgstar_decorator, 0, ierr)

            
         end subroutine plot

         subroutine set_line_style(cnt)
            integer, intent(in) :: cnt
            integer :: iclr, itype
            iclr = cnt - num_colors*(cnt/num_colors) + 1
            call pgsci(colors(iclr))
            if (cnt >= num_colors) then
               itype = Line_Type_Dot
            else
               itype = Line_Type_Solid
            end if
            call pgsls(itype)
         end subroutine set_line_style

      end subroutine do_production_panel

      end module pgstar_production
