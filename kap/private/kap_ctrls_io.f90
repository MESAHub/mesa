! ***********************************************************************
!
! Copyright (C) 2020 Bill Paxton and The MESA Team
!
! MESA is free software; you can use it and/or modify
! it under the combined terms and restrictions of the MESA MANIFESTO
! and the GNU General Library Public License as published
! by the Free Software Foundation; either version 2 of the License,
! or (at your option) any later version.
!
! You should have received a copy of the MESA MANIFESTO along with
! this software; if not, it is available at the mesa website:
! http://mesa.sourceforge.net/
!
! MESA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU Library General Public License for more details.
!
! You should have received a copy of the GNU Library General Public License
! along with this software; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

   module kap_ctrls_io

   use const_def
   use kap_def
   use math_lib

   implicit none

   public :: read_namelist, write_namelist
   private

   real(dp) :: Zbase

   integer :: kap_option, kap_CO_option, kap_lowT_option

   character(len=strlen) :: kap_file_prefix, kap_CO_prefix, kap_lowT_prefix

   ! user table info
   integer :: user_num_kap_Xs = 0
   real(dp), dimension(kap_max_dim) :: user_kap_Xs = -1d0
   integer :: user_num_kap_Zs = 0
   real(dp), dimension(kap_max_dim) :: user_kap_Zs= -1d0
   integer, dimension(kap_max_dim) :: user_num_kap_Xs_for_this_Z = 0

   integer :: user_num_kap_CO_Xs = 0
   real(dp), dimension(kap_max_dim) :: user_kap_CO_Xs = -1d0
   integer :: user_num_kap_CO_Zs = 0
   real(dp), dimension(kap_max_dim) :: user_kap_CO_Zs= -1d0
   integer, dimension(kap_max_dim) :: user_num_kap_CO_Xs_for_this_Z = 0

   integer :: user_num_kap_lowT_Xs = 0
   real(dp), dimension(kap_max_dim) :: user_kap_lowT_Xs = -1d0
   integer :: user_num_kap_lowT_Zs = 0
   real(dp), dimension(kap_max_dim) :: user_kap_lowT_Zs= -1d0
   integer, dimension(kap_max_dim) :: user_num_kap_lowT_Xs_for_this_Z = 0


   real(dp) :: kap_blend_logT_upper_bdy, kap_blend_logT_lower_bdy

   logical :: cubic_interpolation_in_X
   logical :: cubic_interpolation_in_Z
   logical :: include_electron_conduction

   logical :: use_Zbase_for_Type1
   logical :: use_Type2_opacities

   real(dp) :: kap_Type2_full_off_X, kap_Type2_full_on_X
   real(dp) :: kap_Type2_full_off_dZ, kap_Type2_full_on_dZ

   logical :: show_info

   ! hooks
   logical :: use_other_elect_cond_opacity, use_other_compton_opacity

   ! debugging
   logical :: dbg
   real(dp) :: logT_lo, logT_hi
   real(dp) :: logRho_lo, logRho_hi
   real(dp) :: X_lo, X_hi
   real(dp) :: Z_lo, Z_hi

   logical :: read_extra_kap_inlist1
   character (len=128) :: extra_kap_inlist1_name

   logical :: read_extra_kap_inlist2
   character (len=128) :: extra_kap_inlist2_name

   logical :: read_extra_kap_inlist3
   character (len=128) :: extra_kap_inlist3_name

   logical :: read_extra_kap_inlist4
   character (len=128) :: extra_kap_inlist4_name

   logical :: read_extra_kap_inlist5
   character (len=128) :: extra_kap_inlist5_name


   namelist /kap/ &

      Zbase, & 

      kap_file_prefix, kap_CO_prefix, kap_lowT_prefix, aesopus_filename, &

      user_num_kap_Xs, user_kap_Xs, &
      user_num_kap_Zs, user_kap_Zs, user_num_kap_Xs_for_this_Z, &

      user_num_kap_CO_Xs, user_kap_CO_Xs, &
      user_num_kap_CO_Zs, user_kap_CO_Zs, user_num_kap_CO_Xs_for_this_Z, &

      user_num_kap_lowT_Xs, user_kap_lowT_Xs, &
      user_num_kap_lowT_Zs, user_kap_lowT_Zs, user_num_kap_lowT_Xs_for_this_Z, &

      kap_blend_logT_upper_bdy, kap_blend_logT_lower_bdy, &

      cubic_interpolation_in_X, cubic_interpolation_in_Z, &

      include_electron_conduction, &

      use_Zbase_for_Type1, use_Type2_opacities, &

      kap_Type2_full_off_X, kap_Type2_full_on_X, &
      kap_Type2_full_off_dZ, kap_Type2_full_on_dZ, &

      show_info, &

      use_other_elect_cond_opacity, use_other_compton_opacity, &

      read_extra_kap_inlist1, extra_kap_inlist1_name, &
      read_extra_kap_inlist2, extra_kap_inlist2_name, &
      read_extra_kap_inlist3, extra_kap_inlist3_name, &
      read_extra_kap_inlist4, extra_kap_inlist4_name, &
      read_extra_kap_inlist5, extra_kap_inlist5_name


   contains


   ! read a "namelist" file and set parameters
   subroutine read_namelist(handle, inlist, ierr)
      integer, intent(in) :: handle
      character (len=*), intent(in) :: inlist
      integer, intent(out) :: ierr ! 0 means AOK.
      type (Kap_General_Info), pointer :: rq
      integer :: iz, j
      include 'formats'
      call get_kap_ptr(handle,rq,ierr)
      if (ierr /= 0) return
      call set_default_controls
      call read_controls_file(rq, inlist, 1, ierr)
      if (ierr /= 0) return
   end subroutine read_namelist


   recursive subroutine read_controls_file(rq, filename, level, ierr)
      use ISO_FORTRAN_ENV, only: IOSTAT_END
      character(*), intent(in) :: filename
      type (Kap_General_Info), pointer :: rq
      integer, intent(in) :: level
      integer, intent(out) :: ierr
      logical :: read_extra1, read_extra2, read_extra3, read_extra4, read_extra5
      character (len=128) :: message, extra1, extra2, extra3, extra4, extra5
      integer :: unit

      ierr = 0
      if (level >= 10) then
         write(*,*) 'ERROR: too many levels of nested extra controls inlist files'
         ierr = -1
         return
      end if

      if (len_trim(filename) > 0) then
         open(newunit=unit, file=trim(filename), &
            action='read', delim='quote', status='old', iostat=ierr)
         if (ierr /= 0) then
            if (level == 1) then
               ierr = 0 ! no inlist file so just use defaults
               call store_controls(rq, ierr)
            else
               write(*, *) 'Failed to open kap namelist file ', trim(filename)
            end if
            return
         end if
         read(unit, nml=kap, iostat=ierr)
         close(unit)
         if (ierr == IOSTAT_END) then ! end-of-file means didn't find an &kap namelist
            ierr = 0
            write(*, *) 'WARNING: Failed to find kap namelist in file: ', trim(filename)
            call store_controls(rq, ierr)
            close(unit)
            return
         else if (ierr /= 0) then
            write(*, *)
            write(*, *)
            write(*, *)
            write(*, *)
            write(*, '(a)') 'Failed while trying to read kap namelist file: ' // trim(filename)
            write(*, '(a)') 'Perhaps the following runtime error message will help you find the problem.'
            write(*, *)
            open(newunit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
            read(unit, nml=kap)
            close(unit)
            return
         end if
      end if

      call store_controls(rq, ierr)

      if (len_trim(filename) == 0) return

      ! recursive calls to read other inlists

      read_extra1 = read_extra_kap_inlist1
      read_extra_kap_inlist1 = .false.
      extra1 = extra_kap_inlist1_name
      extra_kap_inlist1_name = 'undefined'

      read_extra2 = read_extra_kap_inlist2
      read_extra_kap_inlist2 = .false.
      extra2 = extra_kap_inlist2_name
      extra_kap_inlist2_name = 'undefined'

      read_extra3 = read_extra_kap_inlist3
      read_extra_kap_inlist3 = .false.
      extra3 = extra_kap_inlist3_name
      extra_kap_inlist3_name = 'undefined'

      read_extra4 = read_extra_kap_inlist4
      read_extra_kap_inlist4 = .false.
      extra4 = extra_kap_inlist4_name
      extra_kap_inlist4_name = 'undefined'

      read_extra5 = read_extra_kap_inlist5
      read_extra_kap_inlist5 = .false.
      extra5 = extra_kap_inlist5_name
      extra_kap_inlist5_name = 'undefined'

      if (read_extra1) then
         call read_controls_file(rq, extra1, level+1, ierr)
         if (ierr /= 0) return
      end if

      if (read_extra2) then
         call read_controls_file(rq, extra2, level+1, ierr)
         if (ierr /= 0) return
      end if

      if (read_extra3) then
         call read_controls_file(rq, extra3, level+1, ierr)
         if (ierr /= 0) return
      end if

      if (read_extra4) then
         call read_controls_file(rq, extra4, level+1, ierr)
         if (ierr /= 0) return
      end if

      if (read_extra5) then
         call read_controls_file(rq, extra5, level+1, ierr)
         if (ierr /= 0) return
      end if

   end subroutine read_controls_file


   subroutine set_default_controls
      include 'kap.defaults'
   end subroutine set_default_controls


   subroutine store_controls(rq, ierr)
      type (Kap_General_Info), pointer :: rq

      integer :: i, ierr

      rq% Zbase = Zbase

      rq% cubic_interpolation_in_X = cubic_interpolation_in_X
      rq% cubic_interpolation_in_Z = cubic_interpolation_in_Z
      rq% include_electron_conduction = include_electron_conduction
      rq% use_Zbase_for_Type1 = use_Zbase_for_Type1
      rq% use_Type2_opacities = use_Type2_opacities

      ! check for limits on full_off/on options
      if (kap_Type2_full_off_X > 0.71d0) then
         write(*,*) "kap_Type2_full_off_X must be smaller than 0.71"
         ierr = -1
         return
      end if
      if (kap_Type2_full_on_X > 0.71d0) then
         write(*,*) "kap_Type2_full_on_X must be smaller than 0.71"
         ierr = -1
         return
      end if
      if (kap_Type2_full_off_X < kap_Type2_full_on_X) then
         write(*,*) "kap_Type2_full_off_X has to be bigger than kap_Type2_full_on_X"
         ierr = -1
         return
      end if
      if (kap_Type2_full_off_dZ > kap_Type2_full_on_dZ) then
         write(*,*) "kap_Type2_full_off_dZ has to be smaller than kap_Type2_full_on_dZ"
         ierr = -1
         return
      end if

      rq% kap_Type2_full_off_X = kap_Type2_full_off_X
      rq% kap_Type2_full_on_X = kap_Type2_full_on_X
      rq% kap_Type2_full_off_dZ = kap_Type2_full_off_dZ
      rq% kap_Type2_full_on_dZ = kap_Type2_full_on_dZ


      if (kap_blend_logT_upper_bdy > 0) rq% kap_blend_logT_upper_bdy = kap_blend_logT_upper_bdy
      if (kap_blend_logT_lower_bdy > 0) rq% kap_blend_logT_lower_bdy = kap_blend_logT_lower_bdy


      kap_option = 0
      do i=1,kap_options_max
         if (kap_file_prefix == kap_option_str(i)) then
            kap_option = i
            exit
         end if
      end do
      if (kap_option == 0) then
         write(*,*) 'WARNING: unknown kap_file_prefix (assuming user table): ' // trim(kap_file_prefix)
         kap_option = kap_user
         kap_option_str(kap_user) = trim(kap_file_prefix)

         if (user_num_kap_Xs == 0 .or. user_num_kap_Zs == 0) then
            write(*,*) 'ERROR: must set user_num_kap_Xs, user_num_kap_Zs, and related variables'
            ierr = -1
            return
         end if

         if (user_num_kap_Xs > kap_max_dim .or. user_num_kap_Zs > kap_max_dim) then
            write(0,*) ' failed in kap_read_config_file: maximum X or Z dimensions exceeded'
            write(0,*) ' maximum dimension is ', kap_max_dim
            write(0,*) ' num_kap_Xs = ', num_kap_Xs
            write(0,*) ' num_kap_Zs = ', num_kap_Zs
            ierr = -1
            return
         endif

         num_kap_Xs(kap_user) = user_num_kap_Xs
         kap_Xs(:, kap_user) = user_kap_Xs

         num_kap_Zs(kap_user) = user_num_kap_Zs
         kap_Zs(:, kap_user) = user_kap_Zs

         num_kap_Xs_for_this_Z(:, kap_user) = user_num_kap_Xs_for_this_Z

      end if
      rq% kap_option = kap_option


      kap_CO_option = 0
      do i=1,kap_CO_options_max
         if (kap_CO_prefix == kap_CO_option_str(i)) then
            kap_CO_option = i
            exit
         end if
      end do
      if (kap_CO_option == 0) then
         write(*,*) 'WARNING: unknown kap_CO_prefix (assuming user table): ' // trim(kap_CO_prefix)
         kap_CO_option = kap_CO_user
         kap_CO_option_str(kap_CO_user) = trim(kap_CO_prefix)

         if (user_num_kap_CO_Xs == 0 .or. user_num_kap_CO_Zs == 0) then
            write(*,*) 'ERROR: must set user_num_kap_CO_Xs, user_num_kap_CO_Zs, and related variables'
            ierr = -1
            return
         end if

         num_kap_CO_Xs(kap_CO_user) = user_num_kap_CO_Xs
         kap_CO_Xs(:, kap_CO_user) = user_kap_CO_Xs

         num_kap_CO_Zs(kap_CO_user) = user_num_kap_CO_Zs
         kap_CO_Zs(:, kap_CO_user) = user_kap_CO_Zs

         num_kap_CO_Xs_for_this_Z(:, kap_CO_user) = user_num_kap_CO_Xs_for_this_Z

      end if
      rq% kap_CO_option = kap_CO_option


      kap_lowT_option = 0
      do i=1,kap_lowT_options_max
         if (kap_lowT_prefix == kap_lowT_option_str(i)) then
            kap_lowT_option = i
            exit
         end if
      end do
      if (kap_lowT_option == 0) then
         write(*,*) 'WARNING: unknown kap_lowT_prefix (assuming user table): ' // trim(kap_lowT_prefix)
         kap_lowT_option = kap_lowT_user
         kap_lowT_option_str(kap_lowT_user) = trim(kap_lowT_prefix)

         if (user_num_kap_lowT_Xs == 0 .or. user_num_kap_lowT_Zs == 0) then
            write(*,*) 'ERROR: must set user_num_kap_lowT_Xs, user_num_kap_lowT_Zs, and related variables'
            ierr = -1
            return
         end if

         if (user_num_kap_lowT_Xs > kap_max_dim .or. user_num_kap_lowT_Zs > kap_max_dim) then
            write(0,*) ' failed in kap_read_config_file: maximum X or Z dimensions exceeded'
            write(0,*) ' maximum dimension is ', kap_max_dim
            write(0,*) ' num_kap_lowT_Xs = ', num_kap_lowT_Xs
            write(0,*) ' num_kap_lowT_Zs = ', num_kap_lowT_Zs
            ierr = -1
            return
         endif

         num_kap_lowT_Xs(kap_lowT_user) = user_num_kap_lowT_Xs
         kap_lowT_Xs(:, kap_lowT_user) = user_kap_lowT_Xs

         num_kap_lowT_Zs(kap_lowT_user) = user_num_kap_lowT_Zs
         kap_lowT_Zs(:, kap_lowT_user) = user_kap_lowT_Zs

         num_kap_lowT_Xs_for_this_Z(:, kap_lowT_user) = user_num_kap_lowT_Xs_for_this_Z

      end if
      rq% kap_lowT_option = kap_lowT_option


      ! this parameter needs to be removed
      rq% min_logT_for_logR_gt_1 = 3.3d0

      rq% show_info = show_info

      rq% use_other_elect_cond_opacity = use_other_elect_cond_opacity
      rq% use_other_compton_opacity = use_other_compton_opacity

   end subroutine store_controls


   subroutine write_namelist(handle, filename, ierr)
      integer, intent(in) :: handle
      character(*), intent(in) :: filename
      integer, intent(out) :: ierr
      type (Kap_General_Info), pointer :: rq
      integer :: iounit
      open(newunit=iounit, file=trim(filename), &
         action='write', status='replace', iostat=ierr)
      if (ierr /= 0) then
         write(*,*) 'failed to open ' // trim(filename)
         return
      endif
      call get_kap_ptr(handle,rq,ierr)
      if (ierr /= 0) then
         close(iounit)
         return
      end if
      call set_controls_for_writing(rq)
      write(iounit, nml=kap, iostat=ierr)
      close(iounit)
   end subroutine write_namelist


   subroutine set_controls_for_writing(rq)
      type (Kap_General_Info), pointer :: rq
   end subroutine set_controls_for_writing


   end module kap_ctrls_io
