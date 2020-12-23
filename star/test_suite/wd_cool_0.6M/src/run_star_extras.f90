! ***********************************************************************
!
!   Copyright (C) 2011  Bill Paxton
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

      
      ! these routines are called by the standard run_star check_model
      contains

      include "test_suite_extras.inc"

      subroutine blouin_elect_cond_opacity( &
         zbar, logRho, logT, &
         kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         use const_def, only: dp
         use auto_diff
         use kap_lib, only: kap_get_elect_cond_opacity
         real(dp), intent(in) :: zbar ! average ionic charge (for electron conduction)
         real(dp), intent(in) :: logRho ! the density
         real(dp), intent(in) :: logT ! the temperature
         real(dp), intent(out) :: kap ! electron conduction opacity
         real(dp), intent(out) :: dlnkap_dlnRho ! partial derivative at constant T
         real(dp), intent(out) :: dlnkap_dlnT   ! partial derivative at constant Rho
         integer, intent(out) :: ierr ! 0 means AOK.

         ! this implements the correction formulae from Blouin et al. (2020)
         ! https://ui.adsabs.harvard.edu/abs/2020ApJ...899...46B/abstract

         real(dp), parameter :: alpha_H = -0.52, alpha_He = -0.46
         real(dp), parameter :: a_H = 2.0, a_He = 1.25
         real(dp), parameter :: b_H = 10.0, b_He = 2.5
         real(dp), parameter :: logRho0_H = 5.45, logRho0_He = 6.50
         real(dp), parameter :: logT0_H = 8.40, logT0_He = 8.57
         real(dp), parameter :: sigRho_H = 5.14, sigRho_He = 6.20
         real(dp), parameter :: sigT_H = 0.45, sigT_He = 0.55

         type(auto_diff_real_2var_order1) :: logRho_auto, logT_auto, lnkap_auto
         type(auto_diff_real_2var_order1) :: Rhostar, Tstar
         type(auto_diff_real_2var_order1) :: g_H, g_He, H_H, H_He
         type(auto_diff_real_2var_order1) :: log_correction, log_correction_H, log_correction_He

         real(dp) :: alfa, frac_H, frac_He

         ! auto_diff
         ! var1: lnRho
         ! var2: lnT

         logRho_auto = logRho
         logRho_auto% d1val1 = iln10
         logT_auto = logT
         logT_auto% d1val2 = iln10

         ! correction for H
         Rhostar = logRho_auto - logRho0_H
         Tstar = logT_auto - logT0_H

         g_H = a_H * exp(&
            -pow2(Tstar*cos(alpha_H) + Rhostar*sin(alpha_H))/pow2(sigT_H)  &
            -pow2(Tstar*sin(alpha_H) - Rhostar*cos(alpha_H))/pow2(sigRho_H))

         H_H = 0.5d0 * tanh(b_H*(g_H-0.5d0)) + 0.5d0

         log_correction_H = -log(1d0 + g_H * H_H)

         ! correction for He
         Rhostar = logRho_auto - logRho0_He
         Tstar = logT_auto - logT0_He

         g_He = a_He * exp(&
            -pow2(Tstar*cos(alpha_He) + Rhostar*sin(alpha_He))/pow2(sigT_He)  &
            -pow2(Tstar*sin(alpha_He) - Rhostar*cos(alpha_He))/pow2(sigRho_He))

         H_He = 0.5d0 * tanh(b_He*(g_He-0.5d0)) + 0.5d0

         log_correction_He = -log(1d0 + g_He * H_He)


         ! combined correction
         !
         ! The thermal conductivity is tabulated at Zbar = {1,2,3,4,6,...}
         ! and linear interpolation in logK vs logZbar is applied.
         !
         ! Therefore, we apply the Blouin+ 2020 corrections in a manner
         ! equivalent to individually correcting the Zbar = {1,2} tables.

         if (Zbar .le. 1d0) then ! all H
            frac_H = 1d0
            frac_He = 0d0
         else if (Zbar .le. 2d0) then ! mix H and He
            alfa = (log10(Zbar) - log10(1d0)) / (log10(2d0) - log10(1d0))
            frac_H = 1d0 - alfa
            frac_He = alfa
         else if (Zbar .le. 3d0) then ! mix He and no correction
            alfa = (log10(Zbar) - log10(2d0)) / (log10(3d0) - log10(2d0))
            frac_H = 0d0
            frac_He = 1d0 - alfa
         else ! no correction
            frac_H = 0d0
            frac_He = 0d0
         end if

         ! blend H correction, He correction, and other correction (none)
         log_correction = frac_H * log_correction_H + frac_He * log_correction_He


         ! call standard routines (Cassisi/Potekhin)
         call kap_get_elect_cond_opacity( &
            zbar, logRho, logT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)

         ! pack results in auto_diff
         lnkap_auto = log(kap)
         lnkap_auto% d1val1 = dlnkap_dlnRho
         lnkap_auto% d1val2 = dlnkap_dlnT

         ! apply correction factor
         lnkap_auto = lnkap_auto + log_correction

         ! unpack auto_diff
         kap = exp(lnkap_auto% val)
         dlnkap_dlnRho = lnkap_auto% d1val1
         dlnkap_dlnT = lnkap_auto% d1val2

      end subroutine blouin_elect_cond_opacity


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

         ! internal kap hook
         s% kap_rq% other_elect_cond_opacity => blouin_elect_cond_opacity

      end subroutine extras_controls
      
      
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
         how_many_extra_profile_columns = 2
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, kmax
         real(dp) :: max_frac, frac
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         names(1) = 'zbar_div_abar'
         names(2) = 'eosPC_fraction'
         max_frac = 0d0
         kmax = 0
         do k=1,s% nz
            vals(k,1) = s% zbar(k)/s% abar(k)
            frac = s% eos_frac_PC(k)
            vals(k,2) = frac
            if (frac > max_frac) then
               max_frac = frac
               kmax = k
            end if
         end do
         if (.false. .and. kmax > 0) write(*,*) 'max PC fraction', &
            s% model_number, kmax, max_frac, &
            s% lnd(kmax)/ln10, s% lnT(kmax)/ln10, &
            1d0 - (s% x(kmax) + s% y(kmax)), s% abar(kmax), s% zbar(kmax)
      end subroutine data_for_extra_profile_columns
            

      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
      end function extras_finish_step
      
      

      end module run_star_extras
      
