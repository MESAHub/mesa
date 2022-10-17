! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************


! create zams basically just does a series of create pre-MS,
! run until it reaches zams, and save it to a file for future use.
! that can be efficient if you are going to do many runs
! using the same composition but different masses.

! however, making a zams file is a big job and probably only
! worth the effort if you are going to run a very large number
! of cases.   For most applications it is much easier just to
! create a single model using create_pre_main_sequence_model.



module run_star_extras
   
   use star_lib
   use star_def
   use const_def
   use math_lib
   use run_star_support
   
   implicit none
   
   include "test_suite_extras_def.inc"
   
   
   ! controls for create zams
   character (len = 256) :: zams_name
   real(dp) :: create_z, create_y, mlo, mhi, dmass
   namelist /create_zams_job/ &
      zams_name, create_z, create_y, mlo, mhi, dmass


contains

include "test_suite_extras.inc"
   
   
   subroutine do_run
      
      integer :: id, ierr
      type (star_info), pointer :: s
      character (len = 128) :: zams_inlist
      real(dp) :: dt
      
      write(*, *) 'do create zams'
      
      ierr = 0
      call test_suite_startup(s, .false., ierr)
      
      call do_read_star_job('inlist', ierr)
      if (failed('do_read_star_job')) return
      
      id = id_from_read_star_job
      id_from_read_star_job = 0
      
      call star_ptr(id, s, ierr)
      if (failed('star_ptr')) return
      
      call starlib_init(s, ierr)
      if (failed('star_init')) return
      
      s% inlist_fname = 'inlist'
      
      call star_set_kap_and_eos_handles(id, ierr)
      if (failed('set_star_kap_and_eos_handles')) return
      
      call star_setup(id, 'inlist', ierr)
      if (failed('star_setup')) return
      
      call extras_controls(s% id, ierr)
      if (failed('extras_controls')) return
      
      call do_star_job_controls_before(id, s, .false., ierr)
      if (failed('do_star_job_controls_before')) return
      
      zams_inlist = 'inlist_zams_specification'
      call do_create_zams(&
         s, zams_inlist, s% job% history_columns_file, &
         s% job% profile_columns_file, ierr)
      
      call test_suite_after_evolve(s, ierr)
   
   
   contains
      
      
      logical function failed(str)
         character (len = *), intent(in) :: str
         failed = (ierr /= 0)
         if (failed) then
            write(*, *) trim(str) // ' ierr', ierr
         end if
      end function failed
   
   
   end subroutine do_run
   
   
   subroutine do_create_zams(&
      s, zams_inlist, history_columns_file_in, profile_columns_file_in, ierr)
      use mtx_lib, only : lapack_decsol
      use num_lib
      type (star_info), pointer :: s
      character (len = *) :: zams_inlist, history_columns_file_in, profile_columns_file_in
      integer, intent(out) :: ierr
      
      integer :: io_ms_mod, io_ms_index
      real(dp) :: init_m
      integer :: i, j, k, n, id, result, result_reason, pre_ms_relax_num_steps
      character (len = 256) :: ms_file
      character (len = 1024) :: line
      logical :: okay
      
      1 format(a40, 1pe26.16)
      2 format(a40, i6, 1pe26.16)
      3 format(a15, 2x, f15.6)
      14 format(a40, e24.14)
      
      ierr = 0
      id = s% id
      
      pre_ms_relax_num_steps = 1
      
      call read_zams_controls(s, zams_inlist, ierr)
      if (failed('read_zams_controls')) return
      
      okay = .true.
      
      ierr = 0
      
      do j = 1, 10
         write(*, *)
      end do
      
      ms_file = trim(zams_name) // '.data'
      write(*, *) 'creating ' // trim(ms_file)
      open(newunit = io_ms_index, file = trim(ms_file), action = 'write', status = 'replace')
      
      ms_file = trim(zams_name) // '_mod.data'
      open(newunit = io_ms_mod, file = trim(ms_file), action = 'write', status = 'replace')
      n = (mhi - mlo) / dmass + 1
      
      write(*, 1) 'mlo', mlo
      write(*, 1) 'mhi', mhi
      
      mass_loop : do i = 1, n
         
         init_m = exp10(mlo + (i - 1) * dmass)
         
         s% mesh_delta_coeff = 0.5d0
         if (init_m > 1) s% mesh_delta_coeff = 0.8d0
         if (init_m > 80) s% mesh_delta_coeff = 1
         
         if (init_m > exp10(mhi)) exit
         do j = 1, 10
            write(*, *)
         end do
         
         s% initial_z = create_z
         s% initial_y = create_y
         s% initial_mass = init_m
         write(*, 14) 'do ' // trim(zams_name), s% initial_mass
         
         if (i==1) call write_index_head
         
         call star_create_pre_ms_model(&
            id, s% job% pre_ms_T_c, s% job% pre_ms_guess_rho_c, &
            s% job% pre_ms_d_log10_P, &
            s% job% pre_ms_logT_surf_limit, s% job% pre_ms_logP_surf_limit, &
            s% job% initial_zfracs, &
            s% job% dump_missing_metals_into_heaviest, &
            .false., '', s% job% pre_ms_relax_num_steps, ierr)
         if (failed('star_create_pre_ms_model')) exit
         
         call evolve_to_zams(s, id, ierr)
         if (failed('evolve_to_zams')) exit
         
         call write_model(id, io_ms_mod, io_ms_index, ierr)
         if (failed('write_model')) exit
      
      end do mass_loop
      
      11 format(3x, f15.8, i15)
      write(io_ms_index, 11) -1d0, -1 ! marks end of index
      write(io_ms_index, *) ! blank line at end of index
      
      close (io_ms_mod)
      open(newunit = io_ms_mod, file = trim(ms_file), action = 'read', status = 'old', iostat = ierr)
      if (failed('open mods to read')) return
      
      do
         read(io_ms_mod, fmt = '(a)', iostat = ierr) line
         if (ierr /= 0) then
            ierr = 0; exit
         end if
         write(io_ms_index, fmt = '(a)') trim(line)
      end do
      
      close (io_ms_mod)
      close (io_ms_index)
      
      call free_star(id, ierr)
      if (failed('free_star')) return
      
      call starlib_shutdown
      
      write(*, *)
      if (okay) then
         write(*, '(a)') 'finished create main sequence'
      else
         write(*, '(a)') 'failed during attempt to create main sequence'
      end if
      write(*, *)
   
   contains
      
      subroutine write_index_head
         use chem_def
         use net_def
         integer :: i, time_vals(8)
         character (len = 10) :: date_str, time_str, zone_str
         type (star_info), pointer :: s
         1 format(a32, 2x, 1pe26.14)
         2 format(a32, 2x, i9)
         3 format(a32, 3x, a8)
         4 format(a32, 3x, a)
         call star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*, *) 'write_model: star_ptr failed'
            return
         end if
         write(io_ms_index, '(a,/)') '          1 -- mesa/star zams'
         call date_and_time(date_str, time_str, zone_str, time_vals)
         ! write property list
         write(io_ms_index, 3) 'year_month_day_when_created', date_str(1:8)
         write(io_ms_index, 4) 'version_number', version_number
         write(io_ms_index, 4) 'compiler', compiler_name
         write(io_ms_index, 4) 'build', compiler_version_name
         write(io_ms_index, 4) 'MESA_SDK_version', mesasdk_version_name
         write(io_ms_index, 4) 'math_backend', math_backend
         write(io_ms_index, 4) 'net_name', "'basic.net'"
         write(io_ms_index, 2) 'species', 8
         write(io_ms_index, 1) 'initial_z', create_z
         write(io_ms_index, 1) 'initial_y', create_y
         write(io_ms_index, *) ! blank line for end of property list
         write(io_ms_index, '(a)') '          M/Msun           n_shells'
      end subroutine write_index_head
      
      logical function failed(str)
         character (len = *), intent(in) :: str
         failed = (ierr /= 0)
         if (failed) then
            write(*, *)
            write(*, *) trim(str) // ' ierr', ierr
            okay = .false.
            !call mesa_error(__FILE__,__LINE__)
         end if
      end function failed
   
   end subroutine do_create_zams
   
   
   subroutine write_model(id, io_ms_mod, io_ms_index, ierr)
      use chem_def
      integer, intent(in) :: id, io_ms_mod, io_ms_index
      integer, intent(out) :: ierr
      
      integer :: k, j, species, nz
      type (star_info), pointer :: s
      1 format(a32, 2x, 1pe26.16)
      2 format(a32, 2x, i9)
      11 format(3x, f15.8, i15)
      
      call star_ptr(id, s, ierr)
      if (ierr /= 0) then
         write(*, *) 'write_model: star_ptr failed'
         return
      end if
      
      species = s% species
      nz = s% nz
      
      write(io_ms_index, 11) s% star_mass, nz
      
      ! write property list
      write(io_ms_mod, 1) 'M/Msun', s% star_mass
      write(io_ms_mod, 2) 'n_shells', nz
      write(io_ms_mod, *) ! blank line for end of property list
      
      write(io_ms_mod, fmt = '(7x, a9, 1x, 99(a24, 1x))', advance = 'no') &
         'lnd', 'lnT', 'lnR', 'L', 'dq'
      do j = 1, species
         write(io_ms_mod, fmt = '(a24, 1x)', advance = 'no') trim(chem_isos% name(s% chem_id(j)))
      end do
      write(io_ms_mod, *)
      do k = 1, nz
         write(io_ms_mod, fmt = '(i5, 1x)', advance = 'no') k
         write(io_ms_mod, fmt = '(99(1pe24.16, 1x))', advance = 'no') &
            s% lnd(k), s% lnT(k), s% lnR(k), s% L(k), s% dq(k)
         do j = 1, species
            write(io_ms_mod, fmt = '(1pe24.16, 1x)', advance = 'no') s% xa(j, k)
         end do
         write(io_ms_mod, *)
      end do
      write(io_ms_mod, *)
   
   end subroutine write_model
   
   
   subroutine evolve_to_zams(s, id, ierr)
      type (star_info), pointer :: s
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      integer, parameter :: lipar = 0, lrpar = 0
      logical, parameter :: restore_at_end = .true.
      integer, target :: ipar_ary(lipar)
      real(dp), target :: rpar_ary(lrpar)
      integer, pointer :: ipar(:)
      real(dp), pointer :: rpar(:)
      ipar => ipar_ary
      rpar => rpar_ary
      call star_evolve_to_check_point(&
         id, before_evolve_to_zams, evolve_to_zams_adjust_model, evolve_to_zams_check_model, &
         evolve_to_zams_finish_step, restore_at_end, &
         lipar, ipar, lrpar, rpar, ierr)
   end subroutine evolve_to_zams
   
   
   subroutine before_evolve_to_zams(s, id, lipar, ipar, lrpar, rpar, ierr)
      use star_def, only : star_info
      type (star_info), pointer :: s
      integer, intent(in) :: id, lipar, lrpar
      integer, intent(inout), pointer :: ipar(:) ! (lipar)
      real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
      integer, intent(out) :: ierr
      ierr = 0
   end subroutine before_evolve_to_zams
   
   integer function evolve_to_zams_adjust_model(s, id, lipar, ipar, lrpar, rpar)
      use star_def, only : star_info
      type (star_info), pointer :: s
      integer, intent(in) :: id, lipar, lrpar
      integer, intent(inout), pointer :: ipar(:) ! (lipar)
      real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
      evolve_to_zams_adjust_model = keep_going
   end function evolve_to_zams_adjust_model
   
   integer function evolve_to_zams_check_model(s, id, lipar, ipar, lrpar, rpar)
      use star_def, only : star_info
      type (star_info), pointer :: s
      integer, intent(in) :: id, lipar, lrpar
      integer, intent(inout), pointer :: ipar(:) ! (lipar)
      real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
      evolve_to_zams_check_model = bare_bones_check_model(id)
      if (evolve_to_zams_check_model /= keep_going) return
      if (s% L_nuc_burn_total >= s% L(1) / Lsun) then
         evolve_to_zams_check_model = terminate
         s% termination_code = t_extras_check_model
      end if
   end function evolve_to_zams_check_model
   
   
   integer function evolve_to_zams_finish_step(s)
      type (star_info), pointer :: s
      evolve_to_zams_finish_step = keep_going
   end function evolve_to_zams_finish_step
   
   
   subroutine read_zams_controls(s, zams_inlist, ierr)
      use utils_lib
      type (star_info), pointer :: s
      character (len = *), intent(in) :: zams_inlist
      integer, intent(out) :: ierr
      
      character (len = 256) :: filename, message
      integer :: unit
      
      11 format(a30, f16.6)
      
      ierr = 0
      
      ! set defaults
      create_z = 2d-2
      zams_name = 'z2m2'
      dmass = 0.1d0
      mlo = 0
      mhi = mlo + 2 * dmass
      
      filename = zams_inlist
      open(newunit = unit, file = trim(filename), action = 'read', delim = 'quote', iostat = ierr)
      if (ierr /= 0) then
         write(*, *) 'Failed to open control namelist file ', trim(filename)
      else
         read(unit, nml = create_zams_job, iostat = ierr)
         close(unit)
         if (ierr /= 0) then
            write(*, *) 'Failed while trying to read control namelist file ', trim(filename)
            write(*, '(a)') &
               'The following runtime error message might help you find the problem'
            write(*, *)
            open(newunit = unit, file = trim(filename), action = 'read', &
               delim = 'quote', status = 'old', iostat = ierr)
            read(unit, nml = create_zams_job)
            close(unit)
         end if
      end if
   
   end subroutine read_zams_controls


include 'standard_run_star_extras.inc'

end module run_star_extras
      
