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
      integer :: h1_index
      real(dp) :: initial_h1_mass, h1_acc_abund, h1_init_surf_abund
      logical :: all_ok = .true.

      
      contains

      include "test_suite_extras.inc"
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns

         s% other_photo_read => extras_photo_read
         s% other_photo_write => extras_photo_write

      end subroutine extras_controls


      subroutine extras_photo_read(id, iounit, ierr)
         integer, intent(in) :: id, iounit
         integer, intent(out) :: ierr
         ierr = 0
         read(iounit,iostat=ierr) initial_h1_mass, h1_init_surf_abund
      end subroutine extras_photo_read


      subroutine extras_photo_write(id, iounit)
         integer, intent(in) :: id, iounit
         write(iounit) initial_h1_mass, h1_init_surf_abund
      end subroutine extras_photo_write

      
      subroutine extras_startup(id, restart, ierr)
         use chem_def, only: ih1
!         use adjust_xyz, only: get_xa_for_accretion
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)

         ! get iso index for hydrogen
         h1_index = s% net_iso(ih1)

         ! get abundance of accreted material (for testing profile)
         h1_acc_abund = s% accretion_h1

         if (.not. restart) then

            ! total initial h1 mass in star
            initial_h1_mass = 0
            do k=1,s%nz
               initial_h1_mass = initial_h1_mass + s%xa(h1_index,k) * s%dm(k)
            end do

            ! initial surface abundance (for testing profile)
            h1_init_surf_abund = s%xa( h1_index, 1)
            
         end if

      end subroutine extras_startup
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         logical :: ok
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_after_evolve(s, ierr)

         ok =  check_accretion_profile_ok( s%nz, s%species, s%dm, s%xa, s% star_age * s% star_mdot)
         if (.not. ok) write(*,*) 'check_accretion_profile_ok failed'
         all_ok = all_ok .and. ok

         if (.not. all_ok) then
            write(*,*) '!!!!!!!!!  there are FAILED accretion tests  !!!!!!!!!!!!'
         else
            write(*,*) 'All accretion tests passed successfully'
         end if

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
         how_many_extra_history_columns = 1
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         integer :: k
         real (dp) :: h1_mass, expected_h1_mass_added, h1_mass_added

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         h1_mass = 0.0d0
         do k = 1, s% nz
            h1_mass = h1_mass + s% dm(k) * s% xa(h1_index,k)
         end do
         h1_mass = h1_mass / Msun
         h1_mass_added = h1_mass - initial_h1_mass / Msun

         names(1) = 'h1_mass_added'
         vals(1) = h1_mass_added

      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 1
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
         names(1) = 'zbar_div_abar'
         do k=1,s% nz
            vals(k,1) = s% zbar(k)/s% abar(k)
         end do
      end subroutine data_for_extra_profile_columns
      

      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         logical :: ok

         ierr = 0
         call star_ptr(id, s, ierr)
         ! default return is terminate
         if (ierr /= 0) return

         ok = check_initial_accretion_ok( s%nz, s%species, s%dm, s%xa, s% star_age * s% star_mdot)
         if (.not. ok) write(*,*) 'check_initial_accretion_ok failed for step', s% model_number
         all_ok = all_ok .and. ok

         extras_finish_step = keep_going

      end function extras_finish_step
 

      logical function check_initial_accretion_ok( nz, species, dm, xa, accreted_mass )
         integer, intent(in) :: nz, species
         real (dp), intent(in) :: dm(nz), xa(species,nz), accreted_mass
         integer :: k
         real (dp) :: h1_mass, expected_h1_mass_added, h1_mass_added

         ! get current h1 mass
         h1_mass = 0.0d0
         do k = 1,nz
            h1_mass = h1_mass + dm(k) * xa(h1_index,k)
         end do
         h1_mass = h1_mass / Msun
         h1_mass_added = h1_mass - initial_h1_mass / Msun

         ! compare
         expected_h1_mass_added = accreted_mass * h1_acc_abund

         if ( abs(h1_mass_added-expected_h1_mass_added)/expected_h1_mass_added > 0.005d0 ) then
            write(*,*) 'test of buildup of h1 mass FAILED'
            write(*,*) 'expected h1 mass added: ', expected_h1_mass_added, accreted_mass, h1_acc_abund
            write(*,*) 'current h1 mass added: ', h1_mass_added, h1_mass, initial_h1_mass / Msun
            check_initial_accretion_ok = .false.
            stop 'check_initial_accretion_ok'
         else
            write(*,*) 'added h1 mass increasing at expected rate; expect;', expected_h1_mass_added, &
             'got:', h1_mass_added, 'difference(%):', abs(h1_mass_added-expected_h1_mass_added)/expected_h1_mass_added*100
            write(*,*) 'expected h1 mass added: ', expected_h1_mass_added, accreted_mass, h1_acc_abund
            write(*,*) 'current h1 mass added: ', h1_mass_added, h1_mass, initial_h1_mass / Msun
            check_initial_accretion_ok = .true.
         end if

         
         return
      end function check_initial_accretion_ok

      logical function check_accretion_profile_ok( nz, species, dm, xa, accreted_mass )
         integer, intent(in) :: nz, species
         real (dp), intent(in) :: dm(nz), xa(species,nz), accreted_mass
         real (dp) :: xm, accreted_mass_g
         integer :: k
         logical :: prof_all_ok

         ! need in units of grams for comparing to profile
         accreted_mass_g = accreted_mass * Msun

         ! will reset to false if a problem is discovered while scanning through profile
         prof_all_ok = .true.

         ! check that abundance matches accreted material in 80% of total accreted material
         ! and that it matches the initial surface abundance between 1.2 and 2 times accreted material
         k = 1
         xm = 0
         do while ( xm < 0.8d0 * accreted_mass_g )
            if ( abs( xa(h1_index,k) - h1_acc_abund ) > 0.005d0 ) prof_all_ok = .false.
            xm = xm + dm(k)
            k = k+1
         end do
         if (.not. all_ok) then
            write (*,*) 'abundance of accreted material in profile FAILED to match specified abundance'
         end if
         ! skip through transition region near base of accreted material
         do while ( xm < 1.2d0 * accreted_mass_g )
            xm = xm + dm(k)
            k = k+1
         end do
         ! now in non-accreted material
         do while ( xm < 2 * accreted_mass_g )
            if ( abs( xa(h1_index,k) - h1_init_surf_abund ) > 0.001d0 ) prof_all_ok = .false.
            xm = xm + dm(k)
            k = k+1
         end do
         if (.not. all_ok) then
            write (*,*) 'abundance of material in profile below accreted layer FAILED to match initial abundance'
         end if

         if (prof_all_ok) then
            write(*,*) 'with accreted mass of ', accreted_mass, ' Msun'
            write(*,*) 'abundance profile has clean transition from accreted to underlying abundances'
         end if
            
         check_accretion_profile_ok = prof_all_ok
         
         return
      end function check_accretion_profile_ok

      end module run_star_extras
      
