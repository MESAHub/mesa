! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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
!
! ***********************************************************************

module colors_lib
   use math_lib
   use colors_def
   use const_def
   ! library for calculating theoretical estimates of magnitudes and colors
   ! from Teff, L, M, and [M/H].
   
   ! Color-magnitude data shipped with MESA is from:
   ! Lejeune, Cuisinier, Buser (1998) A&AS 130, 65-75.
   ! However, you add your own bolometric corrections files for mesa to use
   
   ! The data interface for the library is defined in colors_def
   ! Th easiest way to get output is to add the columns to your history_columns.list file
   
   ! The prefered way for users (in a run_star_extras routine) for accessing the colors data is to
   ! call either get_by_by_name, get_abs_mag_by_name or get_abs_bolometric_mag. Other routines are there
   ! to hook into the rest of MESA.
   
   ! Routines get_bc will return the coefficents from interpolating over log Teff, log g, [M/H]
   ! even though the tables are defined as Teff, log g, [M/H]. get_abs_mag routines return
   ! data thats been turned into an absolute magnitude. A color can be computed by taking the difference between
   ! two get_bc or two get_abs_mag calls.
   
   ! Names for the filters should be unique accross all data files (left to the user to enforce this).
   ! Name matching is perfomed in a case sensitive manner.
   ! The names themselves are not important as far as MESA is concerned, you can name each filter (including the
   ! ones MESA ships by defaults) by what ever name you want by editing the data file(s) and changing the names in the header.
   ! MESA does not rely on any particlaur band exisiting.
   
   implicit none


contains ! the procedure interface for the library
   ! client programs should only call these routines.
   
   
   subroutine colors_init(num_files, fnames, num_colors, ierr)
      use mod_colors, only : do_colors_init
      integer, intent(in) :: num_files
      integer, dimension(:), intent(in) :: num_colors
      character(len = *), dimension(:), intent(in) :: fnames
      integer, intent(out) :: ierr
      
      ierr = 0
      
      !$OMP critical (color_init)
      if (.not. color_is_initialized) then
         call do_colors_init(num_files, fnames, num_colors, ierr)
      endif
      !$OMP end critical (color_init)
      
      if(ierr/=0)THEN
         ierr = -1
         write(*, *) "colors_init failed"
         return
      endif
   
   end subroutine colors_init
   
   
   subroutine colors_shutdown ()
      
      use mod_colors, only : free_colors_all
      
      if (.not. color_is_initialized) return
      
      call free_colors_all()
      
      color_is_initialized = .FALSE.
   
   end subroutine colors_shutdown
   
   
   subroutine get_bcs_one(log_Teff, log_g, M_div_h, results, thead, n_colors, ierr)
      use mod_colors, only : Eval_Colors
      ! input
      real(dp), intent(in) :: log_Teff ! log10 of surface temp
      real(dp), intent(in) :: M_div_H ! [M/H]
      ! output
      real(dp), dimension(:), intent(out) :: results
      real(dp), intent(in) :: log_g
      integer, intent(in) :: n_colors
      integer, intent(inout) :: ierr
      type (lgt_list), intent(inout), pointer :: thead
      
      results(:) = -99.9d0
      ierr = 0
      
      if (.not. color_is_initialized) then
         ierr = -1
         return
      endif
      
      call Eval_Colors(log_Teff, log_g, M_div_h, results, thead, n_colors, ierr)
   
   end subroutine get_bcs_one
   
   real(dp) function get_bc_by_name(name, log_Teff, log_g, M_div_h, ierr)
      ! input
      character(len = *), intent(in) :: name
      real(dp), intent(in) :: log_Teff ! log10 of surface temp
      real(dp), intent(in) :: log_g ! log_10 of surface gravity
      real(dp), intent(in) :: M_div_h ! [M/H]
      real(dp), dimension(max_num_bcs_per_file) :: results
      type (lgt_list), pointer :: thead => null()
      integer, intent(inout) :: ierr
      integer :: i, j, n_colors
      
      get_bc_by_name = -99.9d0
      ierr = 0
      
      if (.not. color_is_initialized) then
         ierr = -1
         return
      endif
      
      do i = 1, num_thead
         thead => thead_all(i)%thead
         n_colors = thead_all(i)%n_colors
         do j = 1, n_colors
            if(trim(name)==trim(thead_all(i)%color_names(j)).or. &
               trim(name)=='bc_' // trim(thead_all(i)%color_names(j)).or. &
               trim(name)=='abs_mag_' // trim(thead_all(i)%color_names(j)).or. &
               trim(name)=='lum_band_' // trim(thead_all(i)%color_names(j)).or. &
               trim(name)=='log_lum_band_' // trim(thead_all(i)%color_names(j))&
               ) then
               
               call get_bcs_one(log_Teff, log_g, M_div_h, results, thead, n_colors, ierr)
               if(ierr/=0) return
               
               get_bc_by_name = results(j)
               
               return
            end if
         end do
      
      end do
   
   end function get_bc_by_name
   
   real(dp) function get_bc_by_id(id, log_Teff, log_g, M_div_h, ierr)
      ! input
      integer, intent(in) :: id
      real(dp), intent(in) :: log_Teff ! log10 of surface temp
      real(dp), intent(in) :: log_g ! log_10 of surface gravity
      real(dp), intent(in) :: M_div_h ! [M/H]
      integer, intent(inout) :: ierr
      character(len = strlen) :: name
      
      get_bc_by_id = -99.9d0
      ierr = 0
      
      if (.not. color_is_initialized) then
         ierr = -1
         return
      endif
      
      name = get_bc_name_by_id(id, ierr)
      if(ierr/=0) return
      
      get_bc_by_id = get_bc_by_name(name, log_Teff, log_g, M_div_h, ierr)
   
   end function get_bc_by_id
   
   integer function get_bc_id_by_name(name, ierr)
      ! input
      character(len = *), intent(in) :: name
      integer, intent(inout) :: ierr
      integer :: i, j, k, n_colors
      
      get_bc_id_by_name = -1
      ierr = 0
      
      if (.not. color_is_initialized) then
         ierr = -1
         return
      endif
      
      k = 0
      do i = 1, num_thead
         do j = 1, thead_all(i)%n_colors
            k = k + 1
            if(trim(name)==trim(thead_all(i)%color_names(j)).or. &
               trim(name)=='bc_' // trim(thead_all(i)%color_names(j)).or. &
               trim(name)=='abs_mag_' // trim(thead_all(i)%color_names(j))) then
               get_bc_id_by_name = k
               return
            end if
         end do
      end do
   
   end function get_bc_id_by_name
   
   character(len = strlen) function get_bc_name_by_id(id, ierr)
      ! input
      integer, intent(in) :: id
      integer, intent(inout) :: ierr
      integer :: i, j, k, n_colors
      
      get_bc_name_by_id = ''
      ierr = 0
      
      if (.not. color_is_initialized) then
         ierr = -1
         return
      endif
      
      k = 1
      do i = 1, num_thead
         do j = 1, thead_all(i)%n_colors
            if(k==id) then
               get_bc_name_by_id = thead_all(i)%color_names(j)
               return
            end if
            k = k + 1
         end do
      end do
   
   end function get_bc_name_by_id
   
   real(dp) function get_abs_bolometric_mag(lum)
      use const_def
      real(dp), intent(in) :: lum ! Luminsoity in lsun units
      
      get_abs_bolometric_mag = mbolsun - 2.5d0 * log10(lum)
   
   end function get_abs_bolometric_mag
   
   real(dp) function get_abs_mag_by_name(name, log_Teff, log_g, M_div_h, lum, ierr)
      ! input
      character(len = *) :: name
      real(dp), intent(in) :: log_Teff ! log10 of surface temp
      real(dp), intent(in) :: M_div_h ! [M/H]
      real(dp), intent(in) :: log_g ! log_10 of surface gravity
      real(dp), intent(in) :: lum ! Luminsoity in lsun units
      integer, intent(inout) :: ierr
      
      ierr = 0
      get_abs_mag_by_name = -99.9d0
      
      if (.not. color_is_initialized) then
         ierr = -1
         return
      endif
      
      get_abs_mag_by_name = get_abs_bolometric_mag(lum) - &
         get_bc_by_name(name, log_Teff, log_g, M_div_h, ierr)
   
   end function get_abs_mag_by_name
   
   real(dp) function get_abs_mag_by_id(id, log_Teff, log_g, M_div_h, lum, ierr)
      ! input
      integer, intent(in) :: id
      real(dp), intent(in) :: log_Teff ! log10 of surface temp
      real(dp), intent(in) :: log_g ! log_10 of surface gravity
      real(dp), intent(in) :: M_div_h ! [M/H]
      real(dp), intent(in) :: lum ! Luminsoity in lsun units
      integer, intent(inout) :: ierr
      character(len = strlen) :: name
      
      ierr = 0
      get_abs_mag_by_id = -99.9d0
      
      if (.not. color_is_initialized) then
         ierr = -1
         return
      endif
      
      name = get_bc_name_by_id(id, ierr)
      if(ierr/=0) return
      
      get_abs_mag_by_id = get_abs_mag_by_name(name, log_Teff, log_g, M_div_h, lum, ierr)
   
   end function get_abs_mag_by_id
   
   subroutine get_all_bc_names(names, ierr)
      character(len = strlen), dimension(:) :: names
      integer, intent(inout) :: ierr
      integer :: i, j, cnt
      
      names(:) = ''
      
      if (.not. color_is_initialized) then
         ierr = -1
         return
      endif
      
      cnt = 1
      do i = 1, num_thead
         do j = 1, thead_all(i)%n_colors
            names(cnt) = trim(thead_all(i)%color_names(j))
            cnt = cnt + 1
         end do
      end do
   
   end subroutine get_all_bc_names
   
   subroutine get_bcs_all(log_Teff, log_g, M_div_h, results, ierr)
      ! input
      real(dp), intent(in) :: log_Teff ! log10 of surface temp
      real(dp), intent(in) :: M_div_h ! [M/H]
      ! output
      real(dp), dimension(:), intent(out) :: results
      real(dp), intent(in) :: log_g
      integer, intent(inout) :: ierr
      type (lgt_list), pointer :: thead => null()
      integer :: i, iStart, iEnd
      
      ierr = 0
      results(:) = -99.d0
      
      if (.not. color_is_initialized) then
         ierr = -1
         return
      endif
      
      do i = 1, num_thead
         thead => thead_all(i)%thead
         iStart = (i - 1) * thead_all(i)%n_colors + 1
         iEnd = i * thead_all(i)%n_colors
         call get_bcs_one(log_Teff, log_g, M_div_h, results(iStart:iEnd), thead, thead_all(i)%n_colors, ierr)
         if(ierr/=0) return
      end do
   
   end subroutine get_bcs_all
   
   !Returns in lsun units
   real(dp) function get_lum_band_by_name(name, log_Teff, log_g, M_div_h, lum, ierr)
      ! input
      character(len = *) :: name
      real(dp), intent(in) :: log_Teff ! log10 of surface temp
      real(dp), intent(in) :: M_div_h ! [M/H]
      real(dp), intent(in) :: log_g ! log_10 of surface gravity
      real(dp), intent(in) :: lum ! Total luminsoity in lsun units
      real(dp) :: solar_abs_mag, star_abs_mag
      integer, intent(inout) :: ierr
      
      ierr = 0
      get_lum_band_by_name = -99.d0
      
      if (.not. color_is_initialized) then
         ierr = -1
         return
      endif
      
      ! Filter dependent terms
      solar_abs_mag = get_abs_mag_by_name(name, safe_log10(Teffsun), loggsun, 0.d0, 1.d0, ierr)
      if(ierr/=0) return
      
      star_abs_mag = get_abs_mag_by_name(name, log_Teff, log_g, M_div_h, lum, ierr)
      if(ierr/=0) return
      
      get_lum_band_by_name = exp10((star_abs_mag - solar_abs_mag) / (-2.5d0))
   
   end function get_lum_band_by_name
   
   !Returns in lsun units
   real(dp) function get_lum_band_by_id(id, log_Teff, log_g, M_div_h, lum, ierr)
      ! input
      integer, intent(in) :: id
      real(dp), intent(in) :: log_Teff ! log10 of surface temp
      real(dp), intent(in) :: log_g ! log_10 of surface gravity
      real(dp), intent(in) :: M_div_h ! [M/H]
      real(dp), intent(in) :: lum ! Total luminsoity in lsun units
      real(dp) :: solar_abs_mag, star_abs_mag
      integer, intent(inout) :: ierr
      
      ierr = 0
      get_lum_band_by_id = -99.d0
      
      if (.not. color_is_initialized) then
         ierr = -1
         return
      endif
      
      ! Filter dependent terms
      solar_abs_mag = get_abs_mag_by_id(id, safe_log10(Teffsun), loggsun, 0.d0, 1.d0, ierr)
      if(ierr/=0) return
      
      star_abs_mag = get_abs_mag_by_id(id, log_Teff, log_g, M_div_h, lum, ierr)
      if(ierr/=0) return
      
      get_lum_band_by_id = exp10((star_abs_mag - solar_abs_mag) / (-2.5d0))
   
   end function get_lum_band_by_id

end module colors_lib

