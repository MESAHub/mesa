! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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

      module sample_pgstar_plot
      
      ! you can add your own pgstar plots in addition to the standard ones.
      ! don't edit this file
      ! instead copy the contents to your src/run_star_extras file and edit that.

      ! edit the extras_controls routine to set s% other_pgstar_plots_info
         !    s% other_pgstar_plots_info => my_pgstar_plots_info
      

      use star_lib
      use star_def
      use math_lib

      implicit none
      
      
      ! basics needed for every pgstar plot
      logical :: my_win_flag, my_file_flag
      integer :: my_file_interval
      character (len=256) :: my_file_dir, my_file_prefix
      real :: &
         my_win_width, my_win_aspect_ratio, &
         my_file_width, my_file_aspect_ratio
         
      ! optional extra controls
      character (len=256) :: my_xaxis_by
      real :: &
         my_xmin, my_xmax, &
         my_ymin_left, my_ymax_left, my_dymin_left, &
         my_ymin_right, my_ymax_right, my_dymin_right
      
         
      namelist /my_pgstar/ &
         my_win_flag, my_file_flag, &
         my_file_interval, &
         my_file_dir, my_file_prefix, &
         my_win_width, my_win_aspect_ratio, &
         my_file_width, my_file_aspect_ratio, &
         my_xaxis_by, my_xmin, my_xmax, &
         my_ymin_left, my_ymax_left, my_dymin_left, &
         my_ymin_right, my_ymax_right, my_dymin_right
      

      contains

      
      subroutine my_pgstar_plots_info(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         
         integer, parameter :: num_Other_plots = 1 ! can have up to max_num_Other_plots
         integer :: i, plot_id
         type (pgstar_win_file_data), pointer :: p
         type (star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         call set_my_namelist_defaults
         call read_my_pgstar_namelist('inlist_for_my_pgstar_plots', ierr)
         if (ierr /= 0) return
         
         do i = 1, num_Other_plots
            plot_id = i_Other + i - 1
            p => s% pgstar_win_file_ptr(plot_id)
            p% plot => my_plot
            p% id = plot_id
            p% name = 'My_Plot'
            p% win_flag = my_win_flag
            p% win_width = my_win_width
            p% win_aspect_ratio = my_win_aspect_ratio
            p% file_flag = my_file_flag
            p% file_dir = my_file_dir
            p% file_prefix = my_file_prefix
            p% file_interval = my_file_interval
            p% file_width = my_file_width
            p% file_aspect_ratio = my_file_aspect_ratio
         end do
         
      end subroutine my_pgstar_plots_info
      
      
      subroutine set_my_namelist_defaults
      
         my_win_flag = .false.

         my_win_width = 7
         my_win_aspect_ratio = 0.62 ! aspect_ratio = height/width
         
         my_xaxis_by = 'mass' ! same choices as for main window xaxis_by
         my_xmin = -101 ! only used if > -100
         my_xmax = -101 ! only used if > -100
         
         my_ymin_left = -101 ! only used if > -100
         my_ymax_left = -101 ! only used if > -100        
         my_dymin_left = -101 ! only used if > -100
         
         my_ymin_right = -101 ! only used if > -100
         my_ymax_right = -101 ! only used if > -100        
         my_dymin_right = -101 ! only used if > -100 
         
         ! file output
         my_file_flag = .false.
         my_file_dir = 'pgstar_out'
         my_file_prefix = 'profile'
         my_file_interval = 5 ! output when mod(model_number,my_file_interval)==0
         my_file_width = -1 ! negative means use same value as for window
         my_file_aspect_ratio = -1 ! negative means use same value as for window
         
      end subroutine set_my_namelist_defaults
      
      
      subroutine read_my_pgstar_namelist(filename, ierr)
         use utils_lib
         character(*), intent(in) :: filename
         integer, intent(out) :: ierr

         integer :: unit 
         
         ierr = 0
         
         open(newunit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
         if (ierr /= 0) then
            write(*,'(a)') 'Failed to open control namelist file '// trim(filename)
            return
         end if
         read(unit, nml=my_pgstar, iostat=ierr)  
         close(unit)
         
         if (ierr /= 0) then
            write(*, *) 
            write(*, *) 
            write(*, *) 
            write(*, *) 
            write(*, '(a)') &
               'Failed while trying to read control namelist file: ' // trim(filename)
            write(*, '(a)') &
               'Perhaps the following runtime error message will help you find the problem.'
            write(*, *) 
            open(newunit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
            read(unit, nml=my_pgstar)
            close(unit)
            return
         end if
      
      end subroutine read_my_pgstar_namelist


      subroutine my_plot(id, device_id, ierr)
         integer, intent(in) :: id, device_id
         integer, intent(out) :: ierr
         
         real :: winxmin, winxmax, winymin, winymax, label_scale

         type (star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         call pgslct(device_id)
         call pgbbuf()
         call pgeras()

         winxmin = 0.14
         winxmax = 0.85
         winymin = 0.13
         winymax = 0.92
         label_scale = 1.2
         
         call do_my_plot(s, device_id, &
            winxmin, winxmax, winymin, winymax, &
            label_scale, ierr)

         call pgebuf()
      
      end subroutine my_plot


      subroutine do_my_plot(s, device_id, &
            winxmin, winxmax, winymin, winymax, label_scale, ierr)
            
         use utils_lib
         use const_def

         type (star_info), pointer :: s
         integer, intent(in) :: device_id
         real, intent(in) :: winxmin, winxmax, winymin, winymax, label_scale
         integer, intent(out) :: ierr
         
         real :: windy, xmargin
         real :: xmin, xmax, xleft, xright, dx, tmp, ymin, ymax, ymin2, ymax2, dy
         integer :: grid_min, grid_max, npts, nz
         real, pointer, dimension(:) :: xvec, yvec, yvec2, yvec3
         
         logical :: dbg = .false.
         
         include 'formats.inc'
         ierr = 0
         xmargin = 0

         nz = s% nz
         allocate (xvec(nz), yvec(nz), yvec2(nz), yvec3(nz))
         
         call set_pgstar_xaxis_bounds( &
            s, my_xaxis_by, my_xmin, my_xmax, xmargin, &
            xvec, xmin, xmax, xleft, xright, dx, &
            grid_min, grid_max, npts, ierr)
         if (ierr /= 0) return
         
         if (dbg) then
            write(*,1) 'my_xaxis_by ' // trim(my_xaxis_by)
            write(*,1) 'my_xmin', my_xmin
            write(*,1) 'my_xmax', my_xmax
            write(*,2) 'grid_min', grid_min
            write(*,2) 'grid_max', grid_max
            write(*,1) 'xmin', xmin
            write(*,1) 'xmax', xmax
            write(*,1) 'xleft', xleft
            write(*,1) 'xright', xright
            write(*,1) 'dx', dx
         end if

         call plot(ierr)
         if (ierr /= 0) return

         deallocate(xvec, yvec, yvec2, yvec3)
         
         
         contains
         
         
         subroutine plot(ierr)
            use rates_def, only: i_rate
            use chem_def, only: ipp, icno
            integer, intent(out) :: ierr
            
            integer :: lw, lw_sav, k
            real :: ybot, eps, &
               default_ymax_left, default_ymin_left, &
               default_ymax_right, default_ymin_right
            character (len=128) :: str
         
            include 'formats.inc'
            ierr = 0
            
            call pgsave
                       
            lw = 6
            call pgqlw(lw_sav)
            
            call pgsvp(winxmin, winxmax, winymin, winymax)
            
            ! title
            call pgmtxt('T',1.5,0.5,0.5,'My Plot')
         
            call pgsch(label_scale)
            write(str,'(i9)') s% model_number
            call pgmtxt('T',1.8,0.9,0.5,str)
            
            ! xlabel
            call pgsci(1)
            call pgsch(label_scale)
            call show_pgstar_xaxis_by(s,my_xaxis_by,ierr)
            if (ierr /= 0) return

            
            ! left axis
            
            default_ymax_left = 10
            default_ymin_left = -2
         
            yvec = s% lnd(1:nz)/ln10
            yvec2 = s% lnT(1:nz)/ln10

            if (my_ymax_left > -100) then
               ymax = my_ymax_left
            else
               ymax = max(default_ymax_left,maxval(yvec(grid_min:grid_max)))
               ymax2 = max(default_ymax_left,maxval(yvec2(grid_min:grid_max)))
               if (ymax2 > ymax) ymax = ymax2
            end if
            
            if (my_ymin_left > -100) then
               ymin = my_ymin_left
            else
               ymin = max(default_ymin_left,minval(yvec(grid_min:grid_max)))
               ymin2 = max(default_ymin_left,minval(yvec2(grid_min:grid_max)))
               if (ymin2 < ymin) ymin = ymin2
            end if
            
            dy = ymax-ymin
            if (dy == 0) dy = 1
            if (my_dymin_left > -100) dy = my_dymin_left
            
            ymax = ymax + 0.1*dy
            ymin = ymin - 0.1*dy
            
            if (dbg) then
               write(*,1) 'left axis xleft, xright', xleft, xright
               write(*,1) 'left axis ymin, ymax, dy', ymin, ymax, dy
            end if

            call pgswin(xleft, xright, ymin, ymax)
            call pgscf(1)
            call pgsci(1)
            call pgsch(label_scale)
            call pgbox('',0.0,0,'BNSTV',0.0,0)

            call pgsci(clr_Teal)
            call pgsch(label_scale)
            call pgmtxt('L',3.6,0.5,0.5,'log density (g cm\u-3\d)')
            call pgslw(lw)
            call pgline(npts, xvec(grid_min:grid_max), yvec(grid_min:grid_max))
            call pgslw(lw_sav)

            call pgsci(clr_Coral)
            call pgmtxt('L',5.3,0.5,0.5,'log T (K)')
            call pgslw(lw)
            call pgline(npts, xvec(grid_min:grid_max), yvec2(grid_min:grid_max))
            call pgslw(lw_sav)
            
            
            ! right axis

            lw = 8
            
            default_ymax_right = 0
            default_ymin_right = -10
            
            do k=1,nz
               yvec(k) = safe_log10(s% eps_nuc_categories(ipp,k))
               yvec2(k) = safe_log10(s% eps_nuc_categories(icno,k))
            end do
            
            if (my_ymax_right > -100) then
               ymax = my_ymax_right
            else
               ymax = max(default_ymax_right,maxval(yvec(grid_min:grid_max)))
               ymax2 = max(default_ymax_right,maxval(yvec2(grid_min:grid_max)))
               if (ymax2 > ymax) ymax = ymax2
            end if
            if (my_ymin_right > -100) then
               ymin = my_ymin_right
            else
               ymin = max(default_ymin_right,minval(yvec(grid_min:grid_max)))
               ymin2 = max(default_ymin_right,minval(yvec2(grid_min:grid_max)))
               if (ymin2 < ymin) ymin = ymin2
            end if
            
            dy = ymax-ymin
            if (dy == 0) dy = 1
            if (my_dymin_right > -100) dy = my_dymin_right
            
            ymax = ymax + 0.1*dy
            ymin = ymin - 0.1*dy

            if (dbg) then
               write(*,1) 'right axis xleft, xright', xleft, xright
               write(*,1) 'right axis ymin, ymax, dy', ymin, ymax, dy
            end if
            
            call pgswin(xleft, xright, ymin, ymax)

            call pgscf(1)
            call pgsci(1)
            call pgsch(label_scale)
            call pgbox('BCNST',0.0,0,'CMSTV',0.0,0)
            
            call pgsci(clr_FireBrick)
            call pgmtxt('R',5.6,0.5,0.5,'log eps PP (erg g\u-1\d s\u-1\d)')
            call pgslw(lw)
            call pgline(npts, xvec(grid_min:grid_max), yvec2(grid_min:grid_max))
            call pgslw(lw_sav)
            
            call pgsci(clr_RoyalBlue)
            call pgmtxt('R',3.9,0.5,0.5,'log eps CNO (erg g\u-1\d s\u-1\d)')
            call pgslw(lw)
            call pgline(npts, xvec(grid_min:grid_max), yvec(grid_min:grid_max))
            call pgslw(lw_sav)

            call pgunsa
            
         end subroutine plot
      
         
      end subroutine do_my_plot
      
      

      
      

      end module sample_pgstar_plot
      
      

! inlist_for_my_pgstar_plots

! remove the leading !'s to use this


!&my_pgstar
!
!         my_win_flag = .true.
!
!         my_win_width = 7
!         my_win_aspect_ratio = 0.62 ! aspect_ratio = height/width
!         
!         my_xaxis_by = 'mass' ! same choices as for main window xaxis_by
!         my_xmin = -101 ! only used if > -100
!         my_xmax = -101 ! only used if > -100
!         
!         ! file output
!         my_file_flag = .false.
!         my_file_dir = 'pgstar_out'
!         my_file_prefix = 'profile'
!         my_file_interval = 5 ! output when mod(model_number,my_file_interval)==0
!         my_file_width = -1 ! negative means use same value as for window
!         my_file_aspect_ratio = -1 ! negative means use same value as for window
!
!
!/ ! end of my_pgstar namelist
      
