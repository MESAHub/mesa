! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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
      ! these routines are called by the standard run_star check_model
      contains

      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)
         ! the extras functions in this file will not be called
         ! unless you set their function pointers as done below.
         ! otherwise we use a null_ version which does nothing (except warn).
         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns
         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items
      end subroutine extras_controls

      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_startup

      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0
      end function extras_start_step

      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         use chem_def
         use eos_def
         integer, intent(in) :: id
         integer :: ierr, k, i_accr, iXC, iXO, iXne20, iXNe22, iXNa, iXMg, iXHe
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going
         if (.false. .and. s% star_mass_h1 < 0.35d0) then
            ! stop when star hydrogen mass drops to specified level
            extras_check_model = terminate
            write(*, *) 'have reached desired hydrogen mass'
            return
         end if
         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'
         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
         ! Set iX to the index of o16 and iY to ne20
         iXC = -1
         iXO = -1
         iXNe20 = -1
         iXNe22 = -1
         iXNa = -1
         iXMg = -1
         iXHe = -1
         do k = 1,s% species
             if (chem_isos% name(s% chem_id(k)) == 'c12') then
                iXC = k
             end if
             if (chem_isos% name(s% chem_id(k)) == 'o16') then
                iXO = k
             end if
             if (chem_isos% name(s% chem_id(k)) == 'ne20') then
                iXNe20 = k
             end if
             if (chem_isos% name(s% chem_id(k)) == 'ne22') then
                iXNe22 = k
             end if
             if (chem_isos% name(s% chem_id(k)) == 'na23') then
                iXNa = k
             end if
             if (chem_isos% name(s% chem_id(k)) == 'mg24') then
                iXMg = k
             end if
             if (chem_isos% name(s% chem_id(k)) == 'he4') then
                iXHe = k
             end if
         end do
         if (iXC == -1 .or. iXO == -1 .or. iXNe20 == -1 .or. iXNe22 == -1 .or. iXNa == -1 .or. iXMg == -1 .or. iXHe == -1) then
            write (*,*) 'Could not find elements specified!!'
         end if
         ! Find the base of the oxygen layer
         !do k = 1, s% nz
          !  if (s% xa(iX, k) .le. 0.1) then
                !write(*,*) log10(s% rho(k)), s% xa(iX, k), s% gam(k), s% chiRho(k), s% chiT(k), s% d_eos_dxa(i_lnPgas,iX,k), s% d_eos_dxa(i_lnPgas,iY,k)
           !     exit
            !end if
         !end do
         ! terminate the model when the base density reaches a particular value
         if (.false.) then
         write(*,*) 'base density=', log10(s% rho(s% nz))
         if (s% rho(s% nz) > 1d10) then
            extras_check_model = terminate
            termination_code_str(t_xtra1) = 'base density'
            s% termination_code = t_xtra1
         end if
         end if
      end function extras_check_model

      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 3
      end function how_many_extra_history_columns

      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         use chem_def, only: i3alf
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n), rcore
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: i_max,kc,k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         ! note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in history column options.
         ! it must not include the new column names you are adding here.
         ! find the maximum in the helium burning rate
         i_max = maxloc(s% eps_nuc_categories(i3alf,:), dim=1)
         ! to find the maximum in the total epsilon:
         !i_max = maxloc(s% eps_nuc, dim=1)
         names(1) = "max_eps_he_lgT"
         vals(1) = log10(s% T(i_max))
         names(2) = "mass_core_cryst"
         vals(2) = s% crystal_core_boundary_mass
         do k = s%nz,1,-1
               if(s% m(k) >= s% crystal_core_boundary_mass) then
                  kc = k
                  exit
               end if
         end do
         names(3) = "r_core_cryst"
         vals(3) = exp(s% lnR(kc))
      end subroutine data_for_extra_history_columns

      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 14
      end function how_many_extra_profile_columns

      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use chem_def
         use eos_def
         use eos_lib
         integer :: iXC, iXO, iXne20, iXNe22, iXNa, iXMg, iXHe
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         real(dp), allocatable :: xa1_c12(:), xa1_o16(:), xa1_ne20(:), xa1_ne22(:), xa1_na23(:), xa1_mg24(:), xa1_he4(:)
         real(dp), dimension(num_eos_basic_results) :: res, d_dlnd, d_dlnT
         integer, intent(out) :: ierr
         real(dp), allocatable, dimension(:,:) :: d_dxa
         integer :: k
         real(dp) :: P1, S1, my_chiT, eps, chiX_C12, chiX_O16, chiX_Ne20, chiX_Ne22, chiX_Na23, chiX_Mg24,chiX_He4, bs_C12, bs_O16, bs_Ne20, bs_Ne22, bs_Na23, bs_Mg24, bs_He4, mu1  !!!
         !real(dp) :: chimu_c12,chimu_o16,chimu_ne20,chimu_mg24
         real(dp) :: plnxc_plnye, plnxo_plnye, plnxne20_plnye, plnxne22_plnye, plnxmg_plnye, ln_ye
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         allocate(d_dxa(num_eos_d_dxa_results,s% species))
         allocate(xa1_c12(s% species))
         allocate(xa1_o16(s% species))
         allocate(xa1_ne20(s% species))
         allocate(xa1_ne22(s% species))
         allocate(xa1_na23(s% species))
         allocate(xa1_mg24(s% species))
         allocate(xa1_he4(s% species))
         names(1) = 'chiX_C12'
         names(2) = 'chiX_O16'
         names(3) = 'chiX_Ne20'
         names(4) = 'chiX_Ne22'
         names(5) = 'chiX_Na23'
         names(6) = 'chiX_Mg24'
         names(7) = 'chiX_He4'
         names(8) = 'bs_C12'
         names(9) = 'bs_O16'
         names(10) = 'bs_Ne20'
         names(11) = 'bs_Ne22'
         names(12) = 'bs_N23'
         names(13) = 'bs_Mg24'
         names(14) = 'bs_He4'
         iXC = -1
         iXO = -1
         iXNe20 = -1
         iXNe22 = -1
         iXNa = -1
         iXMg = -1
         iXHe = -1
         do k = 1,s% species
             if (chem_isos% name(s% chem_id(k)) == 'c12') then
                iXC = k
             end if
             if (chem_isos% name(s% chem_id(k)) == 'o16') then
                iXO = k
             end if
             if (chem_isos% name(s% chem_id(k)) == 'ne20') then
                iXNe20 = k
             end if
             if (chem_isos% name(s% chem_id(k)) == 'ne22') then
                iXNe22 = k
             end if
              if (chem_isos% name(s% chem_id(k)) == 'na23') then
                iXNa = k
             end if
             if (chem_isos% name(s% chem_id(k)) == 'mg24') then
                iXMg = k
             end if
             if (chem_isos% name(s% chem_id(k)) == 'he4') then
                iXHe = k
             end if
         end do
         if (iXC == -1 .or. iXO == -1 .or. iXNe20 == -1 .or. iXNe22 == -1 .or. iXMg == -1 .or. iXNa == -1 .or. iXHe == -1) then
            write (*,*) 'Could not find elements specified!!'
         end if
         do k = 1, nz
            eps = 1d-4
            call eosDT_get( &
               s% eos_handle, s% species, s% chem_id, s% net_iso, s% xa(:, k), &
               s% rho(k), log10(s% rho(k)), s% T(k), log10(s% T(k)), &
               res, d_dlnd, d_dlnT, d_dxa, ierr)
            P1 = res(i_lnPgas)
            S1 = exp(res(i_lnS))
            mu1= res(i_mu)!s% mu(k)
            ln_ye=res(i_lnfree_e)
            xa1_c12 = s% xa(:, k)
            xa1_c12(iXC) = xa1_c12(iXC)*(1d0+eps)
            xa1_o16 = s% xa(:, k)
            xa1_o16(iXO) = xa1_o16(iXO)*(1d0+eps)
            xa1_ne20 = s% xa(:, k)
            xa1_ne20(iXNe20) = xa1_ne20(iXNe20)*(1d0+eps)
            xa1_ne22 = s% xa(:, k)
            xa1_ne22(iXNe22) = xa1_ne22(iXNe22)*(1d0+eps)
            xa1_na23 = s% xa(:, k)
            xa1_na23(iXNa) = xa1_na23(iXNa)*(1d0+eps)
            xa1_mg24 = s% xa(:, k)
            xa1_mg24(iXMg) = xa1_mg24(iXMg)*(1d0+eps)
            xa1_he4 = s% xa(:, k)
            xa1_he4(iXHe) = xa1_he4(iXHe)*(1d0+eps)
            call eosDT_get( &
               s% eos_handle, s% species, s% chem_id, s% net_iso, xa1_c12, &
               s% rho(k), log10(s% rho(k)), s% T(k), log10(s% T(k)), &
               res, d_dlnd, d_dlnT, d_dxa, ierr)
            chiX_C12 = (res(i_lnPgas)-P1)/eps
            bs_C12 = - s% xa(iXC,k)*(exp(res(i_lnS))-S1)/(eps* s% xa(iXC,k)) !!-X ds/dX as in Medin & Cumming (2015)
            plnxc_plnye= (log(s% xa(iXC, k))-log(xa1_c12(iXC)))/(ln_ye-res(i_lnfree_e))
            call eosDT_get( &
               s% eos_handle, s% species, s% chem_id, s% net_iso, xa1_o16, &
               s% rho(k), log10(s% rho(k)), s% T(k), log10(s% T(k)), &
               res, d_dlnd, d_dlnT, d_dxa, ierr)
            chiX_O16 = (res(i_lnPgas)-P1)/eps
            bs_O16 = - s% xa(iXO,k)*(exp(res(i_lnS))-S1)/(eps* s% xa(iXO,k))
            plnxo_plnye= (log(s% xa(iXO, k))-log(xa1_o16(iXO)))/(ln_ye-res(i_lnfree_e))
            call eosDT_get( &
               s% eos_handle, s% species, s% chem_id, s% net_iso, xa1_ne20, &
               s% rho(k), log10(s% rho(k)), s% T(k), log10(s% T(k)), &
               res, d_dlnd, d_dlnT, d_dxa, ierr)
            chiX_Ne20 = (res(i_lnPgas)-P1)/eps
            bs_Ne20 = - s% xa(iXNe20,k)*(exp(res(i_lnS))-S1)/(eps* s% xa(iXNe20,k))
            plnxne20_plnye= (log(s% xa(iXNe20, k))-log(xa1_ne20(iXNe20)))/(ln_ye-res(i_lnfree_e))
            call eosDT_get( &
               s% eos_handle, s% species, s% chem_id, s% net_iso, xa1_ne22, &
               s% rho(k), log10(s% rho(k)), s% T(k), log10(s% T(k)), &
               res, d_dlnd, d_dlnT, d_dxa, ierr)
            chiX_Ne22 = (res(i_lnPgas)-P1)/eps
            bs_Ne22 = - s% xa(iXNe22,k)*(exp(res(i_lnS))-S1)/(eps* s% xa(iXNe22,k))
            plnxne22_plnye= (log(s% xa(iXNe22, k))-log(xa1_ne20(iXNe22)))/(ln_ye-res(i_lnfree_e))
            call eosDT_get( &
               s% eos_handle, s% species, s% chem_id, s% net_iso, xa1_na23, &
               s% rho(k), log10(s% rho(k)), s% T(k), log10(s% T(k)), &
               res, d_dlnd, d_dlnT, d_dxa, ierr)
            chiX_Na23 = (res(i_lnPgas)-P1)/eps
            bs_Na23 = - s% xa(iXNa,k)*(exp(res(i_lnS))-S1)/(eps* s% xa(iXNa,k))
            call eosDT_get( &
               s% eos_handle, s% species, s% chem_id, s% net_iso, xa1_mg24, &
               s% rho(k), log10(s% rho(k)), s% T(k), log10(s% T(k)), &
               res, d_dlnd, d_dlnT, d_dxa, ierr)
            chiX_Mg24 = (res(i_lnPgas)-P1)/eps
            bs_Mg24 = - s% xa(iXMg,k)*(exp(res(i_lnS))-S1)/(eps* s% xa(iXMg,k))
            plnxmg_plnye= (log(s% xa(iXMg, k))-log(xa1_mg24(iXMg)))/(ln_ye-res(i_lnfree_e))
            call eosDT_get( &
               s% eos_handle, s% species, s% chem_id, s% net_iso, xa1_he4, &
               s% rho(k), log10(s% rho(k)), s% T(k), log10(s% T(k)), &
               res, d_dlnd, d_dlnT, d_dxa, ierr)
            chiX_He4 = (res(i_lnPgas)-P1)/eps
            bs_He4 = - s% xa(iXHe,k)*(exp(res(i_lnS))-S1)/(eps* s% xa(iXHe,k))
            vals(k,1) = chiX_C12
            vals(k,2) = chiX_O16
            vals(k,3) = chiX_Ne20
            vals(k,4) = chiX_Ne22
            vals(k,5) = chiX_Na23
            vals(k,6) = chiX_Mg24
            vals(k,7) = chiX_He4
            vals(k,8) = bs_C12
            vals(k,9) = bs_O16
            vals(k,10) = bs_Ne20
            vals(k,11) = bs_Ne22
            vals(k,12) = bs_Na23
            vals(k,13) = bs_Mg24
            vals(k,14) = bs_He4
          end do
      end subroutine data_for_extra_profile_columns

      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items

      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return
         ! here is an example for adding an extra history header item
         ! also set how_many_extra_history_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha
      end subroutine data_for_extra_history_header_items

      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items

      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return
         ! here is an example for adding an extra profile header item
         ! also set how_many_extra_profile_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha
      end subroutine data_for_extra_profile_header_items

      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         ! to save a profile,
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.
         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step

      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve
      end module run_star_extras

