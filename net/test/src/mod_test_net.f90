! ***********************************************************************
!
!   Copyright (C) 2009  Bill Paxton, Frank Timmes
!
!   This file is part of MESA.
!
!   MESA is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   MESA is distributed in the hope that it will be useful, 
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

      module mod_test_net
      
      implicit none
         

      contains

      
      subroutine test(quiet)
         use net_def
         use net_lib
         use mod_one_zone_burn, only: do_one_burn
         use test_net_support
         use utils_lib, only: mesa_error
         
         logical, intent(in) :: quiet
         
         integer :: ierr
         
         qt = quiet
               
         call load_libs
      
         test_logT = 7.833d0
         test_logRho = 2d0
         screening_mode = extended_screening
      
         if (.false.) then
            call Do_One_Test('approx21_cr60_plus_co56.net',.false.)
            stop
         end if

         if (.false.) then
            write(*,*) 'inlist_one_zone_burn'
            call do_one_burn('inlist_one_zone_burn',qt)      
            stop
         end if

         if (.false.) then
            write(*,*) 'test_one_zone_burn_const_density'
            call do_one_burn('inlist_one_zone_burn_const_density',qt)      
            stop
         end if

         if (.false.) then
            write(*,*) 'test_one_zone_burn_const_P'
            call do_one_burn('inlist_test_one_zone_burn_const_P',qt)      
            stop
         end if
      
         if (.false.) then
            call Do_One_Test('pp_extras_alternate.net',.false.)
            stop
         end if
      
         if (.false.) then
            test_logRho = 6d0
            test_logT = 9.6d0
            call Do_One_Test('approx21.net',.false.)
            stop
         end if
      
         if (.false.) then
            call Do_One_Test('approx21.net',.false.); stop
         end if
      
         if (.false.) then
            call Do_One_Test_and_show_Qs('pp_and_cno_extras.net',.false.)
            stop
         end if
      
         if (.false.) then
            call Do_One_Test('wd_o_ne_ignite.net',.false.)
            stop
         end if
         
         if (.not. qt) write(*,*) ' **************** basic **************** '
      
         ! 1st one calls test_net_setup -- after than call change_net
         call test_net_setup('basic.net')
         call do_test_net(.false.,.false.)

         if (.not. qt) write(*,*)
         if (.not. qt) write(*,*)
         if (.not. qt) write(*,*) ' **************** o18_and_ne22 **************** '
      
         call change_net('o18_and_ne22.net')
         call do_test_net(.false.,.false.)    
      
         if (.not. qt) write(*,*)
         if (.not. qt) write(*,*)
         if (.not. qt) write(*,*) ' **************** pp_extras **************** '
      
         call change_net('pp_extras.net')
         call do_test_net(.false.,.false.)      
      
         if (.not. qt) write(*,*)
         if (.not. qt) write(*,*)
         if (.not. qt) write(*,*) ' **************** cno_extras **************** '
      
         call change_net('cno_extras.net')
         call do_test_net(.false.,.false.)   

         if (.not. qt) write(*,*)
         if (.not. qt) write(*,*)
         if (.not. qt) write(*,*) ' **************** approx21 **************** '
      
         call change_net('approx21.net')
         test_logRho = 6d0
         test_logT = 9.6d0
         call do_test_net(.false.,.false.)            

         if (.not. qt) write(*,*)
         if (.not. qt) write(*,*)
         if (.not. qt) write(*,*) ' **************** approx21_plus_co56 **************** '
      
         call change_net('approx21_plus_co56.net')
         test_logRho = 6d0
         test_logT = 9.6d0
         call do_test_net(.false.,.false.)            

         if (.not. qt) write(*,*)
         if (.not. qt) write(*,*)
         if (.not. qt) write(*,*) ' **************** approx21_cr60_plus_co56 **************** '
      
         call change_net('approx21_cr60_plus_co56.net')
         call do_test_net(.false.,.false.)            
         
         if (.not. qt) write(*,*)
         if (.not. qt) write(*,*)
      
         if (.not. qt) write(*,*) 'test_one_zone_burn_small_net'
         call do_one_burn('inlist_test_one_zone_burn_small_net',qt)
!      
!         if (.not. qt) write(*,*) 'test_one_zone_burn_big_net'
!         call do_one_burn('inlist_test_one_zone_burn_big_net',qt)
!      
         if (.not. qt) write(*,*) 'test_one_zone_burn_const_P'
         call do_one_burn('inlist_test_one_zone_burn_const_P',qt)      

         call test_net_cleanup
         call net_shutdown
         if (.not. qt) write(*,*)
      
      
      end subroutine test

      
      end module mod_test_net




