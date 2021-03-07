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
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      
      implicit none

      include "test_suite_extras_def.inc"
      
      ! declarations for xtra_coeff_os
      real(dp) :: &
         xtra_coef_os_full_on, &
         xtra_coef_os_full_off, &
         xtra_coef_os_above_nonburn, &
         xtra_coef_os_below_nonburn, &
         xtra_coef_os_above_burn_h, &
         xtra_coef_os_below_burn_h, &
         xtra_coef_os_above_burn_he, &
         xtra_coef_os_below_burn_he, &
         xtra_coef_os_above_burn_z, &
         xtra_coef_os_below_burn_z, &
         xtra_dist_os_above_nonburn, &
         xtra_dist_os_below_nonburn, &
         xtra_dist_os_above_burn_h, &
         xtra_dist_os_below_burn_h, &
         xtra_dist_os_above_burn_he, &
         xtra_dist_os_below_burn_he, &
         xtra_dist_os_above_burn_z, &
         xtra_dist_os_below_burn_z      
      namelist /xtra_coeff_os/ &
         xtra_coef_os_full_on, &
         xtra_coef_os_full_off, &
         xtra_coef_os_above_nonburn, &
         xtra_coef_os_below_nonburn, &
         xtra_coef_os_above_burn_h, &
         xtra_coef_os_below_burn_h, &
         xtra_coef_os_above_burn_he, &
         xtra_coef_os_below_burn_he, &
         xtra_coef_os_above_burn_z, &
         xtra_coef_os_below_burn_z, &
         xtra_dist_os_above_nonburn, &
         xtra_dist_os_below_nonburn, &
         xtra_dist_os_above_burn_h, &
         xtra_dist_os_below_burn_h, &
         xtra_dist_os_above_burn_he, &
         xtra_dist_os_below_burn_he, &
         xtra_dist_os_above_burn_z, &
         xtra_dist_os_below_burn_z
      ! end of declarations for xtra_coeff_os
      
      ! these routines are called by the standard run_star check_model
      contains

      include "test_suite_extras.inc"
      
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         if (s% use_other_mesh_delta_coeff_factor) then ! setup for xtra_coeff_os
            call read_inlist_xtra_coeff_os(ierr)
            if (ierr /= 0) return
            s% other_mesh_delta_coeff_factor => other_mesh_delta_coeff_factor
         end if
         
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  
      end subroutine extras_controls


  !****     code for xtra_coeff_os
      
      subroutine read_inlist_xtra_coeff_os(ierr)
         use utils_lib
         integer, intent(out) :: ierr
         character (len=256) :: filename, message
         integer :: unit
         
         filename = 'inlist_xtra_coeff_os'
         
         write(*,*) 'read_inlist_xtra_coeff_os'
         
         ! set defaults
         xtra_coef_os_full_on = 1d-4
         xtra_coef_os_full_off = 0.1d0
         xtra_coef_os_above_nonburn = 1d0
         xtra_coef_os_below_nonburn = 1d0
         xtra_coef_os_above_burn_h = 1d0
         xtra_coef_os_below_burn_h = 1d0
         xtra_coef_os_above_burn_he = 1d0
         xtra_coef_os_below_burn_he = 1d0
         xtra_coef_os_above_burn_z = 1d0
         xtra_coef_os_below_burn_z = 1d0
         xtra_dist_os_above_nonburn = 0.2d0
         xtra_dist_os_below_nonburn = 0.2d0   
         xtra_dist_os_above_burn_h = 0.2d0
         xtra_dist_os_below_burn_h = 0.2d0   
         xtra_dist_os_above_burn_he = 0.2d0
         xtra_dist_os_below_burn_he = 0.2d0   
         xtra_dist_os_above_burn_z = 0.2d0
         xtra_dist_os_below_burn_z = 0.2d0

         open(newunit=unit, file=trim(filename), action='read', delim='quote', iostat=ierr)
         if (ierr /= 0) then
            write(*, *) 'Failed to open control namelist file ', trim(filename)
         else
            read(unit, nml=xtra_coeff_os, iostat=ierr)  
            close(unit)
            if (ierr /= 0) then
               write(*, *) 'Failed while trying to read control namelist file ', trim(filename)
               write(*, '(a)') &
                  'The following runtime error message might help you find the problem'
               write(*, *) 
               open(newunit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
               read(unit, nml=xtra_coeff_os)
               close(unit)
            end if  
         end if

      end subroutine read_inlist_xtra_coeff_os


      subroutine other_mesh_delta_coeff_factor(id, ierr)
         use const_def
         use math_lib
         use chem_def
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: he_cntr, full_off, full_on, alfa_os
         integer :: k, kk, nz, max_eps_loc
         real(dp) :: xtra_coef, xtra_dist, coef, Hp, r_extra, max_eps
         real(dp) :: eps, eps_h, eps_he, eps_z
         logical :: in_convective_region
         logical, parameter :: dbg = .false.
         
         include 'formats'
         
         !write(*,*) 'enter other_mesh_delta_coeff_factor'
         ierr = 0
         if (xtra_coef_os_above_nonburn == 1d0 .and. &
             xtra_coef_os_below_nonburn == 1d0 .and. &
             xtra_coef_os_above_burn_h == 1d0 .and. &
             xtra_coef_os_below_burn_h == 1d0 .and. &
             xtra_coef_os_above_burn_he == 1d0 .and. &
             xtra_coef_os_below_burn_he == 1d0 .and. &
             xtra_coef_os_above_burn_z == 1d0 .and. &
             xtra_coef_os_below_burn_z == 1d0) return
             
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         nz = s% nz
         he_cntr = s% xa(s% net_iso(ihe4),nz)
         full_off = xtra_coef_os_full_off
         full_on = xtra_coef_os_full_on
         if (he_cntr >= full_off) then
            alfa_os = 0
         else if (he_cntr <= full_on) then
            alfa_os = 1
         else
            alfa_os = (full_off - he_cntr)/(full_off - full_on)
         end if
         !write(*,1) 'alfa_os', alfa_os
         if (alfa_os == 0) return

         ! first go from surface to center doing below convective boundaries
         in_convective_region = (s% mixing_type(1) == convective_mixing)
         k = 2
         max_eps = -1d99
         max_eps_loc = -1

         ! these variables store the component reaction rates at
         ! max_eps_loc
         eps_h = max_eps
         eps_he = max_eps
         eps_z = max_eps

         do while (k <= nz)
            eps = s% eps_nuc(k)
            if (in_convective_region) then
               if (s% mixing_type(k) == convective_mixing) then
                  if (eps > max_eps) then
                     max_eps = eps
                     eps_h = s% eps_nuc_categories(ipp,k) + &
                             s% eps_nuc_categories(icno,k)
                     eps_he = s% eps_nuc_categories(i3alf,k)
                     eps_z = eps - (eps_h + eps_he)
                     max_eps_loc = k
                  end if
               else
                  in_convective_region = .false.
                  if (max_eps < 1d0) then
                     xtra_coef = xtra_coef_os_below_nonburn
                     xtra_dist = xtra_dist_os_below_nonburn
                  else if (eps_h > 0.5d0*max_eps) then
                     xtra_coef = xtra_coef_os_below_burn_h
                     xtra_dist = xtra_dist_os_below_burn_h
                  else if (eps_he > 0.5d0*max_eps) then
                     xtra_coef = xtra_coef_os_below_burn_he
                     xtra_dist = xtra_dist_os_below_burn_he
                  else
                     xtra_coef = xtra_coef_os_below_burn_z
                     xtra_dist = xtra_dist_os_below_burn_z
                  end if
                  xtra_coef = xtra_coef*alfa_os + (1-alfa_os)
                  if (xtra_coef > 0 .and. xtra_coef /= 1) then
                     coef = xtra_coef
                     do
                        if (s% mixing_type(k) /= overshoot_mixing) exit
                        if (coef < s% mesh_delta_coeff_factor(k)) then
                           s% mesh_delta_coeff_factor(k) = coef
                           !write(*,2) 'below mesh_delta_coeff_factor(k)', &
                           !   k, s% mesh_delta_coeff_factor(k)
                        end if
                        if (k == nz) exit
                        k = k+1
                     end do
                     if (xtra_dist > 0) then
                        Hp = s% Peos(k)/(s% rho(k)*s% grav(k))
                        r_extra = max(0d0, s% r(k) - xtra_dist*Hp)
                        if (dbg) write(*,2) 'extra below overshoot region', &
                           k, s% r(k)/Rsun, Hp/Rsun, r_extra/Rsun
                        do
                           if (s% r(k) < r_extra) exit
                           if (coef < s% mesh_delta_coeff_factor(k)) then
                              s% mesh_delta_coeff_factor(k) = coef
                              !write(*,2) 'extra below mesh_delta_coeff_factor(k)', &
                              !   k, s% mesh_delta_coeff_factor(k)
                           end if
                           if (k == nz) exit
                           k = k+1
                        end do
                     end if
                  end if
                  if (dbg) write(*,2) 'done with extra below overshoot region', k
                  if (dbg) write(*,*)
               end if
            else if (s% mixing_type(k) == convective_mixing) then
               in_convective_region = .true.
               max_eps = eps
               max_eps_loc = k
            end if
            k = k+1
         end do

         ! now go from center to surface doing above convective boundaries
         in_convective_region = (s% mixing_type(nz) == convective_mixing)
         k = nz-1
         max_eps = -1d99
         max_eps_loc = -1

         eps_h = max_eps
         eps_he = max_eps
         eps_z = max_eps

         do while (k >= 1)
            eps = s% eps_nuc(k)
            if (in_convective_region) then
               if (s% mixing_type(k) == convective_mixing) then
                  if (eps > max_eps) then
                     max_eps = eps
                     eps_h = s% eps_nuc_categories(ipp,k) + &
                             s% eps_nuc_categories(icno,k)
                     eps_he = s% eps_nuc_categories(i3alf,k)
                     eps_z = eps - (eps_h + eps_he)
                     max_eps_loc = k
                  end if
               else
                  in_convective_region = .false.
                  if (max_eps < 1d0) then
                     xtra_coef = xtra_coef_os_above_nonburn
                     xtra_dist = xtra_dist_os_above_nonburn
                  else if (eps_h > 0.5d0*max_eps) then
                     xtra_coef = xtra_coef_os_above_burn_h
                     xtra_dist = xtra_dist_os_above_burn_h
                  else if (eps_he > 0.5d0*max_eps) then
                     xtra_coef = xtra_coef_os_above_burn_he
                     xtra_dist = xtra_dist_os_above_burn_he
                  else
                     xtra_coef = xtra_coef_os_above_burn_z
                     xtra_dist = xtra_dist_os_above_burn_z
                  end if
                  xtra_coef = xtra_coef*alfa_os + (1-alfa_os)
                  if (dbg) write(*,2) 'xtra_coeff to surf', s% model_number, xtra_coef

                  if (xtra_coef > 0 .and. xtra_coef /= 1) then
                     coef = xtra_coef
                     do
                        if (s% mixing_type(k) /= overshoot_mixing) exit
                        if (coef < s% mesh_delta_coeff_factor(k)) then
                           s% mesh_delta_coeff_factor(k) = coef
                           !write(*,2) 'above mesh_delta_coeff_factor(k)', &
                           !   k, s% mesh_delta_coeff_factor(k)
                        end if
                        if (k == 1) exit
                        k = k-1
                     end do
                     if (xtra_dist > 0) then
                        Hp = s% Peos(k)/(s% rho(k)*s% grav(k))
                        r_extra = min(s% r(1), s% r(k) + xtra_dist*Hp)
                        if (dbg) write(*,2) 'extra above overshoot region', &
                           k, s% r(k)/Rsun, Hp/Rsun, r_extra/Rsun
                        do
                           if (s% r(k) > r_extra) exit
                           if (coef < s% mesh_delta_coeff_factor(k)) then
                              s% mesh_delta_coeff_factor(k) = coef
                              !write(*,2) 'extra above mesh_delta_coeff_factor(k)', &
                              !   k, s% mesh_delta_coeff_factor(k)
                           end if
                           if (k == 1) exit
                           k = k-1
                        end do
                     end if
                  end if
                  if (dbg) write(*,2) 'done with extra above overshoot region', k
                  if (dbg) write(*,*)
               end if
            else if (s% mixing_type(k) == convective_mixing) then
               in_convective_region = .true.
               max_eps = eps
               max_eps_loc = k
            end if
            k = k-1
         end do

      end subroutine other_mesh_delta_coeff_factor

  !****     end of code for xtra_coeff_os
        
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)
      end subroutine extras_startup
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_after_evolve(s, ierr)
      end subroutine extras_after_evolve
      

      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_profile_columns
      

      ! returns either keep_going, retry, or terminate.
      integer function extras_finish_step(id)
         use chem_def
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
      end function extras_finish_step
      
      

      end module run_star_extras
      
