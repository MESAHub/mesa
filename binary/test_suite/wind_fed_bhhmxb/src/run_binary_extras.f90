! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
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
      module run_binary_extras 

      use star_lib
      use star_def
      use const_def
      use math_lib
      use binary_def
      
      implicit none

      integer, parameter :: lx_tested_super_edd = 1
      integer, parameter :: lx_tested_sub_edd = 2
      
      include "binary_test_suite_extras_def.inc"

      contains

      include "binary_test_suite_extras.inc"
      
      subroutine extras_binary_controls(binary_id, ierr)
         integer :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         ! Set these function pointers to point to the functions you wish to use in
         ! your run_binary_extras. Any which are not set, default to a null_ version
         ! which does nothing.
         b% how_many_extra_binary_history_header_items => how_many_extra_binary_history_header_items
         b% data_for_extra_binary_history_header_items => data_for_extra_binary_history_header_items
         b% how_many_extra_binary_history_columns => how_many_extra_binary_history_columns
         b% data_for_extra_binary_history_columns => data_for_extra_binary_history_columns

         b% extras_binary_startup=> extras_binary_startup
         b% extras_binary_start_step=> extras_binary_start_step
         b% extras_binary_check_model=> extras_binary_check_model
         b% extras_binary_finish_step => extras_binary_finish_step
         b% extras_binary_after_evolve=> extras_binary_after_evolve

         ! Once you have set the function pointers you want, then uncomment this (or set it in your star_job inlist)
         ! to disable the printed warning message,
          b% warn_binary_extra =.false.
         
      end subroutine extras_binary_controls


      integer function how_many_extra_binary_history_header_items(binary_id)
         use binary_def, only: binary_info
         integer, intent(in) :: binary_id
         how_many_extra_binary_history_header_items = 0
      end function how_many_extra_binary_history_header_items

      subroutine data_for_extra_binary_history_header_items( &
           binary_id, n, names, vals, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id, n
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
      end subroutine data_for_extra_binary_history_header_items


      integer function how_many_extra_binary_history_columns(binary_id)
         use binary_def, only: binary_info
         integer, intent(in) :: binary_id
         how_many_extra_binary_history_columns = 0
      end function how_many_extra_binary_history_columns
      
      subroutine data_for_extra_binary_history_columns(binary_id, n, names, vals, ierr)
         use const_def, only: dp
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(in) :: n
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         real(dp) :: beta
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         
      end subroutine data_for_extra_binary_history_columns
      
      
      integer function extras_binary_startup(binary_id,restart,ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         logical, intent(in) :: restart
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if
         if (.not. restart) then
            b% lxtra(lx_tested_super_edd) = .false. ! this is used to verify that the system has a super edd phase
            b% lxtra(lx_tested_sub_edd) = .false. ! this is used to verify that the system has a sub edd phase
         end if

         call test_suite_startup(b, restart, ierr)

         extras_binary_startup = keep_going
      end function extras_binary_startup
      
      integer function extras_binary_start_step(binary_id,ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr

         extras_binary_start_step = keep_going
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if
      
      end function  extras_binary_start_step
      
      !Return either keep_going, retry or terminate
      integer function extras_binary_check_model(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if  
         extras_binary_check_model = keep_going
        
      end function extras_binary_check_model
      
      
      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_binary_finish_step(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
         real(dp) :: frac_error
         real(dp) :: actual_mtransfer_rate !actual rate of mass falling into the BH
         real(dp) :: expected_mdot_system !rate of mass loss from vicinity of BH (includes rest mass energy from radiation)
         real(dp) :: expected_mdot_wind !rate of mass loss from the system as a wind from donor
         real(dp) :: donor_wind_mdot
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if  
         extras_binary_finish_step = keep_going

         ! Test that we are not going above Eddington
         if (b% component_mdot(2) > b% mdot_edd*(1-b% mdot_edd_eta)*(1+1e-10)) then
            write(*,*) "ERROR: BH is accreting above Eddington"
            stop
         end if

         ! Test that mass loss from the vicinity of the accretor is correctly set
         ! this also serves to check the accretion luminosity is consistent
         actual_mtransfer_rate = b% mdot_wind_transfer(b% d_i)+b% step_mtransfer_rate
         if (-actual_mtransfer_rate > b% mdot_edd) then
            ! only mdot_edd falls into the black hole, a fraction of that is also removed as radiation
            expected_mdot_system = actual_mtransfer_rate + b% mdot_edd
            expected_mdot_system = expected_mdot_system - b% accretion_luminosity/clight**2
            b% lxtra(lx_tested_super_edd) = .true.
         else
            expected_mdot_system = -b% accretion_luminosity/clight**2
            b% lxtra(lx_tested_sub_edd) = .true.
         end if
         frac_error = abs((b% mdot_system_transfer(2) - expected_mdot_system)/expected_mdot_system)
         if (frac_error > 1e-15) then
            if (b% model_number > 1) then
               write(*,*) "ERROR: mass lost from the vicinity of the BH is incorrect", &
                  b% mdot_system_transfer(2), expected_mdot_system, frac_error
               stop
            end if
         end if

         ! verify that wind mass loss and rlof match mdot of the donor
         ! as well as checking that wind accretion is set right
         donor_wind_mdot = b% s_donor% mstar_dot - b% mtransfer_rate
         frac_error = abs((b% mdot_wind_transfer(b% d_i)-donor_wind_mdot*b% wind_xfer_fraction(b% d_i)) &
            /b% mdot_wind_transfer(b% d_i))
         if (frac_error > 1e-10) then ! there's a loss of precision when adding wind and mtransfer_rate, so test is not that strigent
               write(*,*) "ERROR: unexpected difference in expected wind mdot", &
                  b% mdot_wind_transfer(b% d_i), donor_wind_mdot*b% wind_xfer_fraction(b% d_i), frac_error
               stop
         end if
         expected_mdot_wind = donor_wind_mdot*(1-b% wind_xfer_fraction(b% d_i))
         frac_error = abs((b% mdot_system_wind(1)-expected_mdot_wind)/expected_mdot_wind)
         if (frac_error > 1e-10) then
            write(*,*) "ERROR: calculation of mass leaving system from vicinity of donor failed", &
                  b% mdot_system_wind(1), expected_mdot_wind, frac_error
            stop
         end if
         
      end function extras_binary_finish_step
      
      subroutine extras_binary_after_evolve(binary_id, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if      

         if (b% lxtra(lx_tested_super_edd) .and. b% lxtra(lx_tested_sub_edd)) then
            write(*,*) "Properly tested sub and super Eddington phases"
         else
            write(*,*) "System did not had both sub and super Eddington phases", &
               b% lxtra(lx_tested_super_edd), b% lxtra(lx_tested_sub_edd)
         end if 
 

         call test_suite_after_evolve(b, ierr)

      end subroutine extras_binary_after_evolve     
      
      end module run_binary_extras
