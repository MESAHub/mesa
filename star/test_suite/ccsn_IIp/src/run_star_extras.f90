! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the teerr of the gnu general library public license as published
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
      use auto_diff
      use utils_lib, only: mesa_error, is_bad
            
      implicit none

      integer, parameter :: X_VEL_FRAC_C = 1 ! fraction of c to limit v_center to 
      integer, parameter :: X_STOP_M = 2 ! stop_m

      integer, parameter :: X_STOP_M_FRAC_HE = 7 !for setting stop_m in part2, fraction of He layer.

      integer, parameter :: X_NI_MASS = 12 ! Amount of Ni mass to add Msun

      integer, parameter :: X_RTI_DAYS_OFF = 14 ! days afterwhich to turn off rti if >0
      integer, parameter :: X_MASS_BELOW_SURF = 16 ! mass below surface for stop_m Msun

      integer, parameter :: X_CSM_MDOT = 18 ! mass of csm to add if >0 - part 5 only
      integer, parameter :: X_CSM_MASS = 19 ! mass of csm to add if >0 - part 5 only
      
      integer, parameter :: X_MLT_ALPHA = 21 ! use this mlta_alpha when h1>x_ctrl(X_MLT_H_LIM)
      integer, parameter :: X_MLT_OTHER = 22 ! else use this mlt alpha
      integer, parameter :: X_MLT_H_LIM = 23 ! h limit to switch mlt alpha's

      integer, parameter :: X_CSM_RHO0 = 28
      integer, parameter :: X_CSM_T0 = 29

      integer, parameter :: X_DELTA_LGL_AGE = 32 ! when to change delta_lgL options
      integer, parameter :: X_DELTA_LGL_LIM = 33 ! s% delta_lgL_limit
      integer, parameter :: X_DELTA_LGL_HARD_LIM = 34 ! s% delta_lgL_hard_limit
      
      integer, parameter :: X_NI_MASS_START = 35 ! where to put Ni above M_center
      integer, parameter :: X_NI_MASS_END = 36 ! where to stop putting Ni56 above he core 

      integer, parameter :: X_STELLA_MIN_CNTR_U = 37 !min center velocity for stella.
      
      integer, parameter :: X_SMOOTH_XA_1_START = 45 ! boxcar smooth start mass above M_center
      integer, parameter :: X_SMOOTH_XA_1_END = 46 ! boxcar smooth end mass above he core
      integer, parameter :: X_SMOOTH_XA_1_BOXCAR_MASS = 47 ! boxcar smooth boxcar size

      integer, parameter :: X_SMOOTH_XA_2_START = 48 ! boxcar smooth start mass above M_center
      integer, parameter :: X_SMOOTH_XA_2_END = 49 ! boxcar smooth end mass above he core
      integer, parameter :: X_SMOOTH_XA_2_BOXCAR_MASS = 50 ! boxcar smooth boxcar size

      integer, parameter :: X_MAGNETAR_L_CNTR = 55 ! L_center - Magnetar is only enabled if this is greater than 0
      integer, parameter :: X_MAGNETAR_START_UP = 56 ! start ramping up magnetar at this time in days 
      integer, parameter :: X_MAGNETAR_END_UP = 57 ! stop ramping up magnetar at this time  in days 
      integer, parameter :: X_MAGNETAR_START_DOWN = 58 ! start ramping down magnetar at this time  in days 
      integer, parameter :: X_MAGNETAR_END_DOWN = 59 ! stop ramping down magnetar at this time  in days 

      integer, parameter :: X_FORCE_STOP_M = 98
      integer, parameter :: X_DEFAULT_STOP_M = 99

      integer, parameter :: I_INLIST_PART = 1 ! inlist part number
      integer, parameter :: I_SMOOTH_XA_1_NUM_ITERS = 3 ! boxcar smooth num iters
      integer, parameter :: I_SMOOTH_XA_2_NUM_ITERS = 4 ! boxcar smooth num iters

      integer, parameter :: L_V_CNTR = 1 ! Whether to allow v_center to change

      integer, parameter :: INLIST_INFALL=-1, INLIST_END_INFALL=-2, INLIST_EDEP=-3
      integer, parameter :: INLIST_SHOCK_PART1=1, INLIST_SHOCK_PART2=2, INLIST_SHOCK_PART3=3
      integer, parameter :: INLIST_SHOCK_PART4=4, INLIST_SHOCK_PART5=5


      include "test_suite_extras_def.inc"
      include 'stella/stella_def.inc'

      
      real(dp) :: &
         initial_time, &
         initial_nico, initial_M_center, initial_he_core_mass, initial_mass, &
         start_m, stop_m
         
      real(dp), parameter :: h1_limit = 0.1 ! We use this to check that RTI mixing worked
      real(dp) :: max_mass_h ! Mass cororidnate where h1< h1_limit

      contains

      include "test_suite_extras.inc"


      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return         
         
         include 'stella/stella_controls.inc'
         if (ierr /= 0) return
         
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_start_step => extras_start_step
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns
         s% other_wind => low_density_wind_routine
         s% other_alpha_mlt => alpha_mlt_routine

         s% other_photo_read => extras_photo_read
         s% other_photo_write => extras_photo_write

      end subroutine extras_controls


      include 'stella/stella.inc'


      subroutine alpha_mlt_routine(id, ierr)
         use chem_def, only: ih1
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, h1
         real(dp) :: alpha_H, alpha_other, H_limit
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         alpha_H = s% x_ctrl(X_MLT_ALPHA)
         alpha_other = s% x_ctrl(X_MLT_OTHER)
         H_limit = s% x_ctrl(X_MLT_H_LIM)
         h1 = s% net_iso(ih1)
         if (alpha_H <= 0 .or. alpha_other <= 0 .or. h1 <= 0) return
         do k=1,s% nz
            if (s% xa(h1,k) >= H_limit) then
               s% alpha_mlt(k) = alpha_H
            else
               s% alpha_mlt(k) = alpha_other
            end if
         end do
      end subroutine alpha_mlt_routine


      subroutine low_density_wind_routine(id, Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z, w, ierr)
         use star_def
         integer, intent(in) :: id
         real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z ! surface values (cgs)
         real(dp), intent(out) :: w ! wind in units of Msun/year (value is >= 0)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: msum, lgrho_limit, lnd_limit
         integer :: k, i_lnd
         include 'formats'
         w = 0
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         i_lnd = s% i_lnd
         if (i_lnd <= 0) return
         lgrho_limit = -14.5d0 ! remove material from surface with density below this
         lnd_limit = ln10*lgrho_limit
         msum = 0d0
         do k=1, s% nz
            if (s% xh(i_lnd,k) >= lnd_limit) exit
            msum = msum + s% dm(k)
         end do
         if (msum == 0d0) return
         w = (msum/Msun)/(s% dt/secyer)
         write(*,1) 'low_density_wind_routine lg(Mdot) msum/Msun', safe_log10(w), msum/Msun
      end subroutine low_density_wind_routine
      
      
      subroutine set_nico_mass(s, i_ni56, i_co56, new_ni, final_call, mass_ni56, ierr)
         use chem_def, only: io16
         type (star_info), pointer :: s
         integer, intent(in) :: i_ni56, i_co56
         logical, intent(in) :: final_call
         real(dp), intent(in) :: new_ni
         real(dp), intent(out) :: mass_ni56
         integer, intent(out) :: ierr
         
         integer :: i_o16, k, j, n, nz, species, jmax, kcut
         real(dp) :: old_nico, nico_change, &
            old_o16, new_o16, alfa_o16, alfa_nico, sum_dm, &
            xo16, dm, m00, mp1, mass_co56, mcut, &
            old_xnico, new_xnico, dxnico, x, f, old, new, xsum, &
            check_o16, check_ni56, check_co56, min_m, max_m, new_ni_frac
         include 'formats'
         ierr = 0
         i_o16 = s% net_iso(io16)
         if (i_o16 <= 0) then
            call mesa_error(__FILE__,__LINE__,'need to have o16 in net for set_nico_mass')
         end if
         nz = s% nz
         species = s% species
         do k=1,nz ! fixup abundances before revise
            do j=1,species
               s% xa(j,k) = max(0d0, min(1d0, s% xa(j,k)))
            end do
            xsum = sum(s% xa(1:species,k))
            if (abs(xsum - 1d0) > 1d-12) then
               do j=1,species
                  s% xa(j,k) = s% xa(j,k)/xsum
               end do
            end if
         end do
         min_m = s% x_ctrl(X_NI_MASS_START)*Msun + s% M_center
         max_m = s% x_ctrl(X_NI_MASS_END)*Msun + s% he_core_mass*Msun
         
         if (s% u_flag) then
            do k=nz,1,-1
               if (s% u(k) > s% x_ctrl(X_STELLA_MIN_CNTR_U)) then 
                  ! prepare for removal before give to Stella
                  if (s% m(k) > min_m) then
                     min_m = s% m(k)
                     write(*,1) 'increase min_m to ', min_m/Msun
                     exit
                  end if
               end if
            end do 
         end if        
         
         write(*,1) 'max_m min_m he_core_mass new_ni', &
            max_m/Msun, min_m/Msun, s% he_core_mass, new_ni
         if (max_m > 0d0) then
            do k=1,nz ! replace ni56 + co56 by o16
               s% xa(i_o16,k) = s% xa(i_o16,k) + s% xa(i_ni56,k) + s% xa(i_co56,k)
               s% xa(i_co56,k) = 0d0
               s% xa(i_ni56,k) = 0d0
            end do
            max_m = min(max_m, s% m(1))
            if (min_m < 0) then
               min_m = max_m - new_ni*Msun
               mp1 = s% m(1)
               do k=1,nz ! add ni56 from max_m to min_m
                  m00 = mp1
                  mp1 = m00 - s% dm(k)
                  if (mp1 >= max_m) cycle
                  if (m00 <= max_m .and. mp1 >= min_m) then
                     s% xa(1:species,k) = 0d0
                     s% xa(i_ni56,k) = 1d0
                  else
                     dm = min(m00, max_m) - max(mp1, min_m)
                     f = max(0d0, min(1d0, 1d0 - dm/s% dm(k)))
                     do j=1,species
                        s% xa(j,k) = f*s% xa(j,k)
                     end do
                     s% xa(i_ni56,k) = 0d0
                     s% xa(i_ni56,k) = 1d0 - sum(s% xa(1:species,k))
                  end if
               end do
            else ! min_m >= 0
               min_m = max(min_m, s% M_center)
               if (min_m >= max_m) then
                  write(*,1) 'min_m >= max_m', min_m/Msun, max_m/Msun, s% M_center/Msun
                  call mesa_error(__FILE__,__LINE__,'set_nico_mass')
               end if
               new_ni_frac = new_ni*Msun/(max_m - min_m)
               write(*,1) 'new_ni_frac min_m/Msun max_m/Msun', new_ni_frac, min_m/Msun, max_m/Msun
               mp1 = s% m(1)
               do k=1,nz ! add ni56 from max_m to min_m
                  m00 = mp1
                  mp1 = m00 - s% dm(k)
                  if (m00 <= min_m) exit
                  if (mp1 >= max_m) cycle
                  if (m00 <= max_m .and. mp1 >= min_m) then
                     dm = s% dm(k)
                  else
                     dm = min(m00, max_m) - max(mp1, min_m)
                  end if
                  s% xa(i_ni56,k) = 0d0
                  jmax = maxloc(s% xa(1:species,k),dim=1)
                  if (s% xa(jmax,k) < new_ni_frac) then
                     write(*,2) 'too much ni for cell', k
                     call mesa_error(__FILE__,__LINE__,'set_nico_mass')
                  end if
                  s% xa(i_ni56,k) = new_ni_frac*dm/s% dm(k)
                  s% xa(jmax,k) = s% xa(jmax,k) - s% xa(i_ni56,k)
                  if (is_bad(s% xa(i_ni56,k))) then
                     write(*,2) 's% xa(i_ni56,k)', k, s% xa(i_ni56,k)
                  end if
               end do
            end if
            mass_ni56 = dot_product(s% dm(1:nz), s% xa(i_ni56,1:nz))/Msun
            if (abs(mass_ni56 - new_ni) > 1d-6*new_ni) then
               write(*,1) 'bad new ni', new_ni, mass_ni56
               call mesa_error(__FILE__,__LINE__,'set_nico_mass')
            end if
            write(*,1) 'new ni56 mass', mass_ni56
            !call mesa_error(__FILE__,__LINE__,'set_nico_mass')
            return
         end if
         
         write(*,*) 'rescale Ni profile'
         do k=1,nz
            do j=1,species
               s% xa(j,k) = max(0d0, min(1d0, s% xa(j,k)))
            end do
            xsum = sum(s% xa(1:species,k))
            if (abs(xsum - 1d0) > 1d-12) then
               do j=1,species
                  s% xa(j,k) = s% xa(j,k)/xsum
               end do
            end if
         end do
         
         kcut = nz
         if (final_call .and. stella_skip_inner_dm > 0d0) then
            mcut = s% M_center + stella_skip_inner_dm*Msun
            do k = nz, 1, -1
               if (s% m(k) >= mcut .and. &
                   (stella_skip_inner_v_limit < 0d0 .or. &
                    s% u(k) >= stella_skip_inner_v_limit)) then
                  kcut = k
                  exit
               end if
            end do
         end if

         do k=1,nz ! replace co56 by ni56
            s% xa(i_ni56,k) = s% xa(i_ni56,k) + s% xa(i_co56,k)
            s% xa(i_co56,k) = 0d0
         end do
         ! only count cells 1:kcut
         old_o16 = dot_product(s% dm(1:kcut), s% xa(i_o16,1:kcut))/Msun
         mass_ni56 = dot_product(s% dm(1:kcut), s% xa(i_ni56,1:kcut))/Msun
         mass_co56 = dot_product(s% dm(1:kcut), s% xa(i_co56,1:kcut))/Msun
         old_nico = mass_ni56 + mass_co56
         nico_change = new_ni - old_nico
         if (new_ni < 0d0) then
            write(*,1) 'new_ni', new_ni
            call mesa_error(__FILE__,__LINE__,'bad mass set_nico_mass set_nico_mass')
         end if
         alfa_nico = new_ni/old_nico
         do k=1,kcut
            xo16 = s% xa(i_o16,k)
            old_xnico = s% xa(i_ni56,k) + s% xa(i_co56,k)
            new_xnico = alfa_nico*old_xnico
            dxnico = new_xnico - old_xnico
            s% xa(i_co56,k) = 0d0
            s% xa(i_ni56,k) = new_xnico
            if (dxnico <= xo16) then
               s% xa(i_o16,k) = xo16 - dxnico
            else
               j = maxloc(s% xa(1:species,k),dim=1)
               s% xa(j,k) = s% xa(j,k) - dxnico
               if (s% xa(j,k) <= 0d0) then
                  write(*,3) 'failed in set_nico_mass', j, k, s% xa(j,k), dxnico
                  call mesa_error(__FILE__,__LINE__,'set_nico_mass')
               end if
            end if
         end do
         check_ni56 = dot_product(s% dm(1:kcut), s% xa(i_ni56,1:kcut))/Msun
         if (abs(check_ni56 - new_ni) > 1d-10*max(1d-10,new_ni)) then
            write(*,1) 'check_ni56 - new_ni', check_ni56 - new_ni, check_ni56, new_ni
            write(*,1) 'alfa_nico', alfa_nico
            call mesa_error(__FILE__,__LINE__,'bad check mass in set_nico_mass')
         end if
         mass_ni56 = new_ni
         write(*,1) 'revised mass Ni56', check_ni56
      end subroutine set_nico_mass
            
      subroutine extras_startup(id, restart, ierr)
         use chem_def, only: ini56, ico56, ih1, ihe4, io16
         use interp_2d_lib_db, only: interp_mkbicub_db
         use eos_lib, only: eosDT_get_T
         use eos_def
         use atm_lib, only: atm_L
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: xni56, xmax, &
            min_mass, max_mass, boxcar_mass, tot_h1, tot_he4
         integer :: i, j, k, kk, k1, k_max_v, ni56, co56, o16, he4, h1, &
            max_iter, eos_calls
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)

         xni56 = 0d0
         ni56 = s% net_iso(ini56)
         co56 = s% net_iso(ico56)
         o16 = s% net_iso(io16)
         he4 = s% net_iso(ihe4)
         h1 = s% net_iso(ih1)
         if (o16 <= 0 .or. he4 <= 0 .or. h1 <= 0) call mesa_error(__FILE__,__LINE__,'missing o16, he4, or h1')
         
         if (s% eos_rq% logRho_min_OPAL_SCVH_limit > -12d0) then
            write(*,'(A)')
            write(*,*)'FIX: have set_logRho_OPAL_SCVH_limits too large'
            write(*,*) 'causes errors in HELM/OPAL blend partials'
            write(*,'(A)')
            call mesa_error(__FILE__,__LINE__,'extras_startup')
         end if
         
         if (.not. restart) then
            initial_nico = 0
            stop_m = 0
            initial_M_center = s% M_center
            initial_mass = s% star_mass
            initial_time = s% time
            initial_he_core_mass = s% he_core_mass
            max_mass_h = -1
         
            if (s% x_integer_ctrl(I_INLIST_PART) == INLIST_EDEP) then
               if (s% total_mass_for_inject_extra_ergs_sec > 0) then ! doing edep
                  if (s% v_flag) then
                     do k=1,s% nz
                        s% xh(s% i_v,k) = 0d0
                        s% v(k) = 0d0
                     end do
                  end if
               end if
               return
            end if

            if (s% x_integer_ctrl(I_INLIST_PART) == INLIST_SHOCK_PART4) then
               ! This is to check that RTI worked
               ! Find outer most location where H1<h1_limit then later we will check if this has changed
               do k=1, s%nz
                  if(s% xa(h1, k) < h1_limit) then
                     max_mass_h = s% m(k)
                     exit
                  end if
               end do
               if(max_mass_h < 0) call mesa_error(__FILE__,__LINE__,'Error: Could not find h1 limit')
            end if

            if (s% x_integer_ctrl(I_INLIST_PART) == INLIST_SHOCK_PART5 .and. &
                s% x_ctrl(X_CSM_MASS) > 0) call add_csm ! part 5, add csm
            
            if (s% x_ctrl(X_SMOOTH_XA_1_BOXCAR_MASS) > 0d0 .and. s% x_integer_ctrl(I_SMOOTH_XA_1_NUM_ITERS) > 0) then
               min_mass = s% x_ctrl(X_SMOOTH_XA_1_START) + s% M_center/Msun
               max_mass = s% x_ctrl(X_SMOOTH_XA_1_END) + s% he_core_mass
               boxcar_mass = s% x_ctrl(X_SMOOTH_XA_1_BOXCAR_MASS)
               call smooth_xa_by_boxcar_mass( &
                  s% id, min_mass, max_mass, boxcar_mass, s% x_integer_ctrl(I_SMOOTH_XA_1_NUM_ITERS), ierr)
               if (ierr /= 0) return
            end if

            if (ni56 > 0 .and. co56 > 0) then
               if (s% x_ctrl(X_NI_MASS) > 0) then
                  call set_nico_mass(s, ni56, co56, s% x_ctrl(X_NI_MASS), .false., xni56, ierr)
                  if (ierr /= 0) return
               end if
               initial_nico = xni56
            end if
            
            if (s% x_ctrl(X_SMOOTH_XA_2_BOXCAR_MASS) > 0d0 .and. s% x_integer_ctrl(I_SMOOTH_XA_2_NUM_ITERS) > 0) then
               min_mass = s% x_ctrl(X_SMOOTH_XA_2_START) + s% M_center/Msun
               max_mass = s% x_ctrl(X_SMOOTH_XA_2_END) + s% he_core_mass
               boxcar_mass = s% x_ctrl(X_SMOOTH_XA_2_BOXCAR_MASS)
               call smooth_xa_by_boxcar_mass( &
                  s% id, min_mass, max_mass, boxcar_mass, s% x_integer_ctrl(I_SMOOTH_XA_2_NUM_ITERS), ierr)
               if (ierr /= 0) return
            end if
         
            if (s% x_integer_ctrl(I_INLIST_PART) == INLIST_SHOCK_PART1) &
               s% cumulative_energy_error = 0d0 ! set to 0 at start of part1
               
            start_m = s% shock_mass
            if (start_m == 0d0) then ! use max v
               k_max_v = maxloc(s% u(1:s% nz),dim=1)
               start_m = s% m(k_max_v)/Msun
               if (start_m < 2d0*s% M_center/Msun) then
                  write(*,2) 'no shock; use max_v for start_m', k_max_v, start_m
               else
                  start_m = s% M_center/Msun
                  write(*,2) 'no shock; use M_center for start_m', s% nz, start_m
               end if
            else
               write(*,2) 'use shock for start_m', s% shock_k, s% shock_mass
            end if
            write(*,2) 'M_center', s% nz, s% M_center/Msun
            
            write(*,1) 's% x_ctrl(X_FORCE_STOP_M)', s% x_ctrl(X_FORCE_STOP_M)
            write(*,2) 's% x_integer_ctrl(I_INLIST_PART)', s% x_integer_ctrl(I_INLIST_PART)
         
            if (s% x_ctrl(X_FORCE_STOP_M) > 0d0) then
               stop_m = s% x_ctrl(X_FORCE_STOP_M)
            else
               stop_m = 0d0
               select case(s% x_integer_ctrl(I_INLIST_PART))
               case (INLIST_SHOCK_PART1)
                  call find_inlist_part1_stop_m()
               case (INLIST_SHOCK_PART2)
                  call find_inlist_part2_stop_m()
               case (INLIST_SHOCK_PART3)
                  call find_inlist_part3_stop_m()
               case (INLIST_SHOCK_PART4)
                  stop_m = s% star_mass - s% x_ctrl(X_MASS_BELOW_SURF)
               case (INLIST_SHOCK_PART5)
                  if (s% x_ctrl(X_MASS_BELOW_SURF) > 0d0) &
                     stop_m = s% star_mass - s% x_ctrl(X_MASS_BELOW_SURF)
               end select
            end if
            
         end if ! not restart
         
         write(*,1) 's% x_ctrl(X_MASS_BELOW_SURF)', s% x_ctrl(X_MASS_BELOW_SURF)
         write(*,1) 's% star_mass', s% star_mass
         write(*,1) 'start_m', start_m
         write(*,1) 'stop_m', stop_m
                  
         if (s% x_ctrl(X_MASS_BELOW_SURF) > 0d0) &
               stop_m = min(stop_m, s% star_mass - s% x_ctrl(X_MASS_BELOW_SURF))
         
         if (start_m > stop_m .and. stop_m > 0d0) then
            write(*,1) 'start_m > stop_m', start_m, stop_m
            call mesa_error(__FILE__,__LINE__,'extras_startup')
         end if

         if (stop_m > 0d0) then
            s% x_ctrl(X_STOP_M) = stop_m
            write(*,1) 'stop when shock reaches', stop_m
            write(*,1) 's% he_core_mass', s% he_core_mass
            write(*,1) 's% co_core_mass', s% co_core_mass
            write(*,1) 's% M_center/Msun', s% M_center/Msun
            write(*,1) 's% star_mass', s% star_mass
            write(*,'(A)')
            !stop
         end if
                  
         contains

         subroutine find_inlist_part1_stop_m()
            include 'formats'
            do k=1,s% nz
               if (s% xa(o16,k) > s% xa(he4,k)) then
                  stop_m = (0.15d0*s% M_center + 0.85d0*s% m(k))/Msun
                     if (stop_m <= start_m) then
                        write(*,2) '1st choice for stop_m too small', k, stop_m
                        stop_m = (0.05d0*s% M_center + 0.95d0*s% m(k))/Msun
                        if (stop_m <= start_m) then
                           write(*,2) '2nd choice for stop_m too small', k, stop_m
                           stop_m = (0.01d0*s% M_center + 0.99d0*s% m(k))/Msun
                           if (stop_m <= start_m) then
                              write(*,2) '3rd choice for stop_m too small', k, stop_m
                              stop_m = 1.001d0*start_m
                           end if
                        end if
                     end if
                  exit
               end if
            end do
            if (stop_m == 0d0) then
               if (s% x_ctrl(X_DEFAULT_STOP_M) <= 0d0) call mesa_error(__FILE__,__LINE__,'failed to find stop_m')
               stop_m = s% x_ctrl(X_DEFAULT_STOP_M)
            end if

         end subroutine find_inlist_part1_stop_m

         subroutine find_inlist_part2_stop_m()
            real(dp) :: f
            include 'formats'
            k1 = 0
            do k=1,s% nz
               if (s% xa(he4,k) > s% xa(h1,k)) then
                  k1 = k
                  exit
               end if
            end do
            if (k1 > 0) then
               do k=k1+1,s% nz
                  if (s% xa(o16,k) > s% xa(he4,k)) then
                     f = s% x_ctrl(X_STOP_M_FRAC_HE)
                     stop_m = ((1d0 - f)*s% m(k) + f*s% m(k1))/Msun
                     exit
                  end if
               end do
            end if
            if (stop_m == 0d0) then
               if (s% x_ctrl(X_DEFAULT_STOP_M) <= 0d0) call mesa_error(__FILE__,__LINE__,'failed to find stop_m')
               stop_m = s% x_ctrl(X_DEFAULT_STOP_M)
            end if


         end subroutine find_inlist_part2_stop_m


         subroutine find_inlist_part3_stop_m()

            do k=1,s% nz
               if (s% xa(he4,k) > s% xa(h1,k)) then
                  stop_m = 0.95d0*s% m(k)/Msun + 0.05d0*(s% star_mass - 0.1d0)
                  exit
               end if
            end do
            if (stop_m == 0d0) call mesa_error(__FILE__,__LINE__,'failed to find stop_m')

         end subroutine find_inlist_part3_stop_m
         
         subroutine add_csm
            real(dp) :: xni56, xmax, &
            logT, P_hse, Z, &
            logT_bnd1, logT_bnd2, logP_at_bnd1, logP_at_bnd2, &
            logT_guess, logT_result, logT_tol, logP_tol, &
            csm_mass, csm_mdot, windv, r0, rho0, T0, L0, r, f, rho, dm, dV, dq
            real(dp), dimension(num_eos_basic_results) :: &
               res, d_dlnd, d_dlnT
            real(dp), allocatable :: d_dxa(:,:)
            include 'formats'
            ierr = 0

            allocate(d_dxa(num_eos_d_dxa_results, s% species))

            csm_mass = s% x_ctrl(X_CSM_MASS)*Msun
            csm_mdot = s% x_ctrl(X_CSM_MDOT)*Msun/secyer
            do kk=2,s% nz-2
               if (s% m(1) - s% m(kk) < csm_mass) cycle
               r0 = s% r(kk)
               if (s% x_ctrl(X_CSM_RHO0) > 0) then
                  rho0 = s% x_ctrl(X_CSM_RHO0)
               else
                  rho0 = s% rho(kk)
               end if
               if (s% x_ctrl(X_CSM_T0) > 0) then
                  T0 = s% x_ctrl(X_CSM_T0)
               else
                  T0 = s% T(kk)
               end if
               L0 = s% L(1)
               windv = csm_mdot/(4*pi*r0*r0*rho0) ! mdot velocity
               !windv = sqrt(2*s% cgrav(kk)*s% m(kk)/r0) ! escape velocity
               write(*,1) 'old log(r(1)/Rsun), R/Rsun', log10(s% r(1)/Rsun), s% r(1)/Rsun
               write(*,2) 'rho0, T0, r0, csm mass, csm v', &
                  kk, rho0, T0, r0/Rsun, (s% m(1) - s% m(kk))/Msun, windv
               dm = sum(s% dm(1:kk-1))/(kk-1) 
               dq = dm/s% xmstar
               do k = 1, kk-1
                  s% dq(k) = dq
                  s% q(k+1) = s% q(k) - dq
                  s% dm(k) = dm
                  s% m(k+1) = s% m(k) - dm
               end do
               logT_bnd1 = arg_not_provided
               logT_bnd2 = arg_not_provided
               logP_at_bnd1 = arg_not_provided
               logP_at_bnd2 = arg_not_provided
               max_iter = 100
               logT_tol = 1d-5
               logP_tol = 1d-8
               k = kk
               write(*,'(A)')
               write(*,2) 'T rho P logR u', k, s% T(k), s% rho(k), s% Peos(k), &
                  log10(s% r(k)/Rsun), s% u(k)
               do k = kk-1, 1, -1
                  r = s% r(k+1) + 0.5d0*(s% r(k+1) - s% r(k+2))
                  f = r0/r; f = f*f
                  rho = rho0*f ! rho proportional to 1/r^2
                  dV = dm/rho
                  r = s% r(k+1)
                  r = pow(dV/(4d0*pi/3d0) + r*r*r, 1d0/3d0)
                  s% r(k) = r
                  s% lnR(k) = log(s% r(k))
                  s% xh(s% i_lnR,k) = s% lnR(k)
                  s% rho(k) = rho
                  s% lnd(k) = log(s% rho(k))
                  s% xh(s% i_lnd,k) = s% lnd(k)
         
         
                  s% u(k) = windv !* r/r0
         
         
                  s% xh(s% i_u,k) = s% u(k)
         
         
                  if (.true.) then ! set T to give P for HSE
                     r = s% r(k+1)
                     P_hse = s% Peos(k+1) - &
                        s% cgrav(k+1)*s% m(k+1)*(s% dm(k+1)+s% dm(k))/(8*pi*r*r*r*r)
                     Z = max(0d0, min(1d0, 1d0 - (s% X(k) + s% Y(k))))
                     logT_guess = s% lnT(k+1)/ln10
                     call eosDT_get_T( &
                        s% eos_handle, &
                        s% species, s% chem_id, s% net_iso, s% xa(:,k), &
                        s% lnd(k)/ln10, i_logPtot, log10(P_hse), &
                        logT_tol, logP_tol, max_iter, logT_guess, &
                        logT_bnd1, logT_bnd2, logP_at_bnd1, logP_at_bnd2, &
                        logT_result, res, d_dlnd, d_dlnT, &
                        d_dxa, eos_calls, ierr)
                     if (ierr /= 0) return
                     s% lnT(k) = logT_result*ln10
                     s% T(k) = exp(s% lnT(k))
                     if (.false. .and. s% T(k) > s% T(k+1)) then ! can happen at base
                        s% T(k) = s% T(k+1)
                        s% lnT(k) = s% lnT(k+1)
                     end if
                  else
                     s% T(k) = T0*f
                     s% lnT(k) = log(s% T(k))
                  end if
                  s% xh(s% i_lnT,k) = s% lnT(k)
         
                  if (.true.) then ! set to black body L
                     s% L(k) = atm_L(s% T(k), s% r(k))
                  else
                     s% L(k) = L0
                  end if
                  s% xh(s% i_lum,k) = s% L(k)

                  write(*,2) 'T rho P logR u', k, s% T(k), s% rho(k), s% Peos(k), &
                     log10(s% r(k)/Rsun), s% u(k)
               end do
               write(*,'(A)')
               write(*,1) 'new log(r(1)/Rsun), R/Rsun', log10(s% r(1)/Rsun), s% r(1)/Rsun
               write(*,'(A)')
               exit
            end do
         end subroutine add_csm

      end subroutine extras_startup
      
      
      subroutine extras_after_evolve(id, ierr)
         use chem_def, only: ini56, ico56, ih1
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: dt, h1_mass
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (s% x_integer_ctrl(I_INLIST_PART) == INLIST_SHOCK_PART5 .and. &
             save_stella_data_when_terminate) then
            call write_stella_data(s, ierr)
            if (ierr /= 0) return
         end if
         
         ! Check that RTI worked
         if (s% x_integer_ctrl(I_INLIST_PART) == INLIST_SHOCK_PART4 .and. s% rti_flag) then
            ! This is only needed for the test suite and can be removed if doing science
            call check_rti()
         end if

         call test_suite_after_evolve(s, ierr)

         contains

         subroutine check_rti()
            logical,parameter :: dbg=.false.

            ! This is to check that RTI worked
            ! Find outer most location where H1<0.1 then later we will check if this has changed
            h1_mass = -1

            if(dbg) write(*,*) "max_mass_h", max_mass_h/msun

            do k=1, s%nz-1
               !write(*,*) k, s% m(k), max_mass_h, s% m(k) .ge. max_mas_h, s% m(k+1).le. max_mass_h,2.3377068540442852E+034 >=  6.6430170767931355E+033
               if(s% m(k)>= max_mass_h .and. s% m(k+1)<= max_mass_h) then
                  h1_mass = s% xa(s% net_iso(s% net_iso(ih1)),k)
                  exit
               end if
            end do
            if(h1_mass < 0) call mesa_error(__FILE__,__LINE__,'Error: Could not find h1 limit')

            if(h1_mass < h1_limit*2) then
               write(*,*) "h1_limit ",h1_limit
               write(*,*) "max_mass_h", max_mass_h/msun
               write(*,*) "h1_mass", h1_mass
               ierr = -1
               call mesa_error(__FILE__,__LINE__,'Error: RTI mxiing did not seem to work')
            end if

         end subroutine check_rti


      end subroutine extras_after_evolve
      

      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         use chem_def, only: ini56
         integer, intent(in) :: id
         integer :: ierr, k
         type (star_info), pointer :: s
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         include 'formats'
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
         integer :: iounit
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
      end subroutine data_for_extra_history_columns
      
      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
         if (s% x_integer_ctrl(I_INLIST_PART) == INLIST_SHOCK_PART3) &
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
         real(dp) :: dtau, kap
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (n == 0) return
         if (s% x_integer_ctrl(I_INLIST_PART) == INLIST_SHOCK_PART3) then
            if (n /= 1) call mesa_error(__FILE__,__LINE__,'bad num cols for data_for_extra_profile_columns')
            names(1) = 'du'
            vals(1,1) = 0
            do k=2,s% nz
               vals(k,1) = 0 ! s% xtra1_array(k)
            end do
            return
         end if
      end subroutine data_for_extra_profile_columns
      

      ! returns either keep_going or terminate.
      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         real(dp) :: alfa, beta, age_days, L_center, next_R_center, &
            start_ramp_up, end_ramp_up, &
            start_ramp_down, end_ramp_down, factor_end, factor_start

         include 'formats'
         extras_start_step = keep_going
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         if (s% x_integer_ctrl(I_INLIST_PART) == INLIST_SHOCK_PART1 .and. s% model_number >= 1000) &
            s% max_timestep = 0 ! turn off limit
            
         age_days = s% star_age*365.25d0
         
         if (s% x_ctrl(X_RTI_DAYS_OFF) > 0 .and. age_days >= s% x_ctrl(X_RTI_DAYS_OFF) .and. &
               s% RTI_C > 0d0) then
            call turn_off_rti()
         end if
         
         call enable_magnetar(s)
                           
         if (age_days >= s% x_ctrl(X_DELTA_LGL_AGE)) then
            call adjust_delta_lgL()
         end if
         
         if (s% x_logical_ctrl(L_V_CNTR) .and. s% dt > 0d0) then
            call adjust_v_center()
         end if


         contains

         subroutine turn_off_rti()
            include 'formats'
            s% RTI_C = 0d0
            s% RTI_log_max_boost = 0d0 
            s% RTI_m_full_boost = -1d0
            s% RTI_m_no_boost = 0d0
            s% dedt_RTI_diffusion_factor = 1d0
            s% composition_RTI_diffusion_factor = 1d0
            write(*,'(A)')
            write(*,2) 'turn off RTI', s% model_number, age_days
            write(*,'(A)')
         end subroutine turn_off_rti


         subroutine adjust_v_center()

            s% v_center = &
               s% v_center - s% cgrav(s% nz)*s% M_center/(s% R_center*s% R_center)
            s% v_center = max(s% v_center, -s% x_ctrl(X_VEL_FRAC_C)*clight)
            next_R_center = s% R_center + s% v_center*s% dt
            if (next_R_center < s% center_R_lower_limit) then
               s% v_center = (s% center_R_lower_limit - s% R_center)/s% dt
            end if

         end subroutine adjust_v_center


         subroutine adjust_delta_lgL()

            s% delta_lgL_limit = s% x_ctrl(X_DELTA_LGL_LIM)
            s% delta_lgL_hard_limit = s% x_ctrl(X_DELTA_LGL_HARD_LIM)

         end subroutine adjust_delta_lgL

      end function extras_start_step
         

      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         use chem_def, only: ih1, ico56, ini56
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         real(dp) :: age_days, log_L, shock_mass, log_Lnuc_burn, xmax, xni56
         integer :: k
         include 'formats'
         extras_finish_step = keep_going
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         if (s% x_ctrl(X_STOP_M) <= 0) return
         shock_mass = s% shock_mass
         if (shock_mass >= s% x_ctrl(X_STOP_M)) then
            write(*,'(a,2f12.5)') 'shock has reached target location', &
               shock_mass, s% x_ctrl(X_STOP_M)
            extras_finish_step = terminate
            s% termination_code = t_extras_finish_step
            call restore_nico_mass()
         else if (shock_mass >= 0.9995d0*s% x_ctrl(X_STOP_M)) then
               write(*,1) 'shock has reached this fraction of target', &
                  shock_mass/s% x_ctrl(X_STOP_M)
         end if
         
         if (s% x_integer_ctrl(I_INLIST_PART) == INLIST_SHOCK_PART5 .and. &
             s% model_number == save_stella_data_for_model_number) then
            call write_stella_data(s, ierr)
            if (ierr /= 0) return
         end if

         contains

         subroutine restore_nico_mass()

            ! restore Ni+Co mass
            if (s% net_iso(ini56) > 0 .and. s% net_iso(ico56) > 0) then
               if (s% x_ctrl(X_NI_MASS) > 0) then
                  call set_nico_mass( &
                     s, s% net_iso(ini56), s% net_iso(ico56), s% x_ctrl(X_NI_MASS), .true., xni56, ierr)
                  if (ierr /= 0) return
               end if
            end if

         end subroutine restore_nico_mass


         
      end function extras_finish_step
           
      subroutine extras_photo_read(id, iounit, ierr)
         integer, intent(in) :: id, iounit
         integer, intent(out) :: ierr
         integer :: inlist_part
         type (star_info), pointer :: s
         ierr = 0
   
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
   
         read(iounit,iostat=ierr) initial_nico, initial_M_center, initial_mass, initial_time, initial_he_core_mass
         read(iounit,iostat=ierr) start_m, stop_m, inlist_part, max_mass_h

         if(inlist_part/= s% x_integer_ctrl(I_INLIST_PART)) then
            call mesa_error(__FILE__,__LINE__,'Error: Photo was saved for different inlist')
            ierr=-1
            return
         end if

         end subroutine extras_photo_read
   
         subroutine extras_photo_write(id, iounit)
         integer, intent(in) :: id, iounit
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
   
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
   
         write(iounit) initial_nico, initial_M_center, initial_mass, initial_time, initial_he_core_mass
         write(iounit) start_m, stop_m, s% x_integer_ctrl(I_INLIST_PART), max_mass_h
   
         end subroutine extras_photo_write


         subroutine enable_magnetar(s)
            type(star_info), pointer :: s
            real(dp) :: age_days, L_center, &
                        start_ramp_up, end_ramp_up, &
                        start_ramp_down, end_ramp_down

            L_center = s% x_ctrl(X_MAGNETAR_L_CNTR)
            if (L_center > 0) then  ! magnetar
               start_ramp_up = s% x_ctrl(X_MAGNETAR_START_UP)
               end_ramp_up = s% x_ctrl(X_MAGNETAR_END_UP)
               start_ramp_down = s% x_ctrl(X_MAGNETAR_START_DOWN)
               end_ramp_down = s% x_ctrl(X_MAGNETAR_END_DOWN)
               if (age_days <= start_ramp_up) then
                  s% L_center = 0d0
               else if (age_days < end_ramp_up) then
                  s% L_center = &
                     L_center*(age_days - start_ramp_up)/(end_ramp_up - start_ramp_up)
               else if (age_days < start_ramp_down) then
                  s% L_center = L_center
               else if (age_days <= end_ramp_down) then
                  s% L_center = &
                     L_center*(end_ramp_down - age_days)/(end_ramp_down - start_ramp_down)
               else
                  s% L_center = 0d0
               end if
            end if

         end subroutine enable_magnetar

      end module run_star_extras
      