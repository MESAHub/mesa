! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
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
      use auto_diff

      implicit none

      include "test_suite_extras_def.inc"

! here are the x controls used below

!alpha_mlt_routine
         !alpha_H = s% x_ctrl(21)
         !alpha_other = s% x_ctrl(22)
         !H_limit = s% x_ctrl(23)

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
         s% other_alpha_mlt => alpha_mlt_routine
!         s% other_adjust_mdot => other_adjust_mdot
s% other_wind => my_erupt_other_wind

      end subroutine extras_controls


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
         alpha_H = s% x_ctrl(21)
         alpha_other = s% x_ctrl(22)
         H_limit = s% x_ctrl(23)
         h1 = s% net_iso(ih1)
         !write(*,1) 'alpha_H', alpha_H
         !write(*,1) 'alpha_other', alpha_other
         !write(*,1) 'H_limit', H_limit
         !write(*,2) 'h1', h1
         !write(*,2) 's% nz', s% nz
         if (alpha_H <= 0 .or. alpha_other <= 0 .or. h1 <= 0) return
         do k=1,s% nz
            if (s% xa(h1,k) >= H_limit) then
               s% alpha_mlt(k) = alpha_H
            else
               s% alpha_mlt(k) = alpha_other
            end if
            !write(*,2) 'alpha_mlt', k, s% alpha_mlt(k),
         end do
         !stop
      end subroutine alpha_mlt_routine


!
!subroutine other_adjust_mdot(id, ierr)
!use star_def
!integer, intent(in) :: id
!integer, intent(out) :: ierr
!integer :: k
!real(dp) :: E_excess, eruptive_ml, cell_energy, E_diff, E_bind
!type(star_info), pointer :: s
!ierr = 0
!call star_ptr(id, s, ierr) ! retrieve star id
!
!E_excess = 0d0
!eruptive_ml = 0d0
!cell_energy = 0d0 ! kinetic + rot + turb + pot , no internal.
!E_diff = 0d0
!E_bind = 0d0
!!  s% xtra1_array(k) , cell Energy excess
!! cell_energy = eval_cell_section_total_energy(s,k,k), cell total energy w/o internal
!
!k = 1
!do while (s% T(k) < 1d7 .and. k <= s% nz)
!cell_energy = eval_cell_section_total_energy(s,k,k)!/s% dm(k)
!E_excess = E_excess + s% xtra1_array(k) ! g/s
! k = k+1
!end do
!
!
!E_diff = E_excess - E_bind
!
!delta_M = s% m(1) - s% m(kinside)
!
!write(*,*) '- eruptive_ml', - eruptive_ml
!write(*,*) 's% delta_mdot', s% mstar_dot
!
!s% mstar_dot = s% mstar_dot - eruptive_ml
!end subroutine other_adjust_mdot

subroutine get_mloss_forwind(id, ierr, mdot_erupt)
   !! ------------------------------------------------------------------------
   !! Implements a streamlined "Shelley-style" eruptive mass-loss check:
   !!   1) Identify base of region where tau(k) < c/csound(k).
   !!   2) Accumulate local-excess energy vs. binding for each zone above that.
   !!   3) Compare mass-averaged E_excess to mass-averaged E_bind.
   !!   4) If E_excess > |E_bind| => remove all mass above that base.
   !!   5) Use an average local dynamical time for the timescale => Mdot.
   !!
   !! No partial sub-steps are done here; MESA’s adaptive timestep
   !! will typically subdivide the actual ejection if it's large.
   !!
   !! On exit, mdot_erupt (Msun/yr) is >= 0.
   !! ------------------------------------------------------------------------
   use star_def
   use const_def, only: clight => clight, Gconst => standard_cgrav
   implicit none

   integer, intent(in)    :: id
   integer, intent(out)   :: ierr
   real(dp), intent(out)  :: mdot_erupt

   type(star_info), pointer :: s

   ! Indices, loop control
   integer :: k, nz, k_massloss

   ! Summation variables
   real(dp) :: sumExcess, sumBind, sumDM, sumDt
   real(dp) :: avgExcess, avgBind, avgDt

   ! Local zone stuff
   real(dp) :: dtLoc, L_excess_loc

   ! Final mass to remove, final dt
   real(dp) :: DeltaM, eruptive_ml

   ! Ejection efficiency factor from Shelley’s simple approach
   ! (feel free to parametrize or pass as an argument instead)
   real(dp), parameter :: xi = 0.3d0

   ! Handy references
   real(dp), parameter :: Msun_cgs = 1.98840987d33
   real(dp), parameter :: sec_to_yr = 3.17098d-8

   ! --------------------------------------------------------------------
   ! attach star pointer
   ! --------------------------------------------------------------------
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) then
      mdot_erupt = 0d0
      return
   end if

   nz = s% nz
   if (nz < 2) then
      mdot_erupt = 0d0
      return
   end if

   ! --------------------------------------------------------------------
   ! 1) find the base of region with tau(k) < c / csound(k)
   !    i.e. from k=1 (surface) inward, stop if tau(k) >= c/csound(k).
   ! --------------------------------------------------------------------
   k_massloss = -1

   ! Initialize sums
   sumExcess = 0d0
   sumBind   = 0d0
   sumDM     = 0d0
   sumDt     = 0d0

   do k = 1, nz

      if ( s% tau(k) >= clight / s% csound(k) ) then
         ! found base => store, then exit
         k_massloss = k
         exit
      end if

      ! local dynamical time
      if (s% m(k) > 0d0 .and. s% r(k) > 0d0) then
         dtLoc = sqrt( (s% r(k)**3) / (standard_cgrav * s% m(k)) )
      else
         dtLoc = 0d0
      end if

      ! local "excess" luminosity (above Eddington)
      ! L_Edd(k) ~ 4 pi c G Mgrav(k)/kappa(k).
      L_excess_loc = s% L(k) - (4d0*acos(-1.d0)*clight*standard_cgrav*s% m_grav(k)/s% opacity(k))
      if (L_excess_loc < 0d0) L_excess_loc = 0d0

      ! accumulate sums:
      ! sumExcess is effectively integral(L_excess_loc * dtLoc) over dm(k)
      sumExcess = sumExcess + L_excess_loc * dtLoc * s% dm(k)

      ! sumBind ~ negative of gravitational potential
      ! snippet uses " - (G m(k)/r(k))*dm(k)"
      sumBind   = sumBind   - ( standard_cgrav * s% m(k) / s% r(k) ) * s% dm(k)

      sumDM     = sumDM     + s% dm(k)

      ! sumDt is for average timescale => we do dtLoc * dm(k),
      ! then we’ll divide by sumDM
      sumDt     = sumDt     + dtLoc * s% dm(k)

   end do

   ! If we never found k_massloss, means even the center is tau < c/csound => big star?
   ! That is unusual, but to be safe:
   if (k_massloss < 0) then
      k_massloss = nz
   end if

   ! --------------------------------------------------------------------
   ! 2) find mass-averaged E_excess vs. E_bind
   ! --------------------------------------------------------------------
   if (sumDM > 0d0) then
      avgExcess = sumExcess / sumDM
      avgBind   = sumBind   / sumDM
      avgDt     = sumDt     / sumDM
   else
      avgExcess = 0d0
      avgBind   = 0d0
      avgDt     = 0d0
   end if

   ! --------------------------------------------------------------------
   ! 3) check if overlaying layers are unbound => E_excess > |E_bind|
   ! --------------------------------------------------------------------
   if (avgExcess > abs(avgBind)) then
      ! star is unbinding all mass above k_massloss
      DeltaM = s% m(1) - s% m(k_massloss)
      if (DeltaM < 0d0) DeltaM = 0d0
   else
      DeltaM = 0d0
   end if

   ! --------------------------------------------------------------------
   ! 4) define the eruptive mass-loss rate [Msun/yr]
   !    negative sign => mass is lost
   ! --------------------------------------------------------------------
   if (avgDt > 0d0 .and. DeltaM > 0d0) then
      eruptive_ml = - xi * (DeltaM / Msun_cgs) / (avgDt * sec_to_yr)
   else
      eruptive_ml = 0d0
   end if

   ! MESA expects a non-negative s% mstar_dot, so we return positive.
   mdot_erupt = -eruptive_ml
   if (mdot_erupt < 0d0) mdot_erupt = 0d0

end subroutine get_mloss_forwind






subroutine my_erupt_other_wind(id, Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z, wind, ierr)
   ! This subroutine sets "wind" = Mdot_total (in Msun/yr).
   use star_def
   implicit none

   integer, intent(in) :: id
   real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z
   real(dp), intent(out) :: wind   ! Msun/yr
   integer, intent(out) :: ierr

   type(star_info), pointer :: s
   real(dp) :: mdot_dutch, mdot_erupt_local
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return

   ! -- Evaluate standard "Dutch" (or other) wind
   !call some_standard_wind(id, Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z, mdot_dutch, ierr)

   ! -- Evaluate Shelley's eruptive mass loss
   call get_mloss_forwind(id, ierr, mdot_erupt_local)
   if (ierr /= 0) return

   ! -- Combine them
   wind = wind + mdot_erupt_local! +mdot_dutch
   write(*,*) 'mdot_erupt =' , mdot_erupt_local
   ! If you want to store each piece in s%xtra(...) for reference:
   !s% xtra(1) = mdot_dutch
  ! s% xtra(2) = mdot_erupt_local
   !s% xtra(3) = wind

end subroutine my_erupt_other_wind

!subroutine other_adjust_mdot(id, ierr)
!  use star_def
!  implicit none
!  integer, intent(in)    :: id
!  integer, intent(out)   :: ierr
!  integer                :: k, nz, k_massloss
!  real(dp)               :: dt, dt_sum, avg_dt
!  real(dp)               :: L_excess
!  real(dp)               :: cum_excess, cum_bind, mass_accum
!  real(dp)               :: avg_excess, avg_bind, Delta_M, eruptive_ml
!  ! xi is our free efficiency parameter (you can later include an extra eta factor)
!  real(dp), parameter  :: xi = 0.3d0
!  type(star_info), pointer :: s
!
!  ierr = 0
!  call star_ptr(id, s, ierr)
!  if (ierr /= 0) return
!
!  ! In your model the mass coordinate is such that:
!  ! • The surface is at index 1,
!  ! • τ starts small at the surface and grows inward,
!  ! • The center is at index nz.
!  nz = s%nz
!
!  cum_excess = 0d0
!  cum_bind   = 0d0
!  mass_accum = 0d0
!  dt_sum     = 0d0
!  k_massloss = -1
!
!  ! Loop from the surface (k=1) inward until the local optical depth exceeds tau_crit.
!  ! (tau_crit = clight/s%csound(k))
!  do k = 1, nz
!     dt = sqrt( (s%r(k)**3) / (standard_cgrav * s%m(k)) )
!     ! Use xtra1_array as a proxy for (L_rad - L_Edd) at zone k.
!     L_excess = s%xtra1_array(k)
!     if (L_excess < 0d0) L_excess = 0d0
!
!     ! In your models τ is smallest at the surface (k=1) and increases inward.
!     ! We integrate only over zones with τ < clight/s%csound(k)
!     if ( s%tau(k) < clight / s%csound(k) ) then
!        cum_excess = cum_excess + L_excess * dt * s%dm(k)
!        cum_bind   = cum_bind   - (standard_cgrav * s%m(k) / s%r(k)) * s%dm(k)
!        mass_accum = mass_accum + s%dm(k)
!        dt_sum = dt_sum + dt * s%dm(k)
!     else
!        ! Once we reach a zone where τ ≥ clight/s%csound, we define that as the base.
!        k_massloss = k
!        exit
!     end if
!  end do
!
!  if (mass_accum > 0d0) then
!    avg_excess = cum_excess / mass_accum
!    avg_bind   = cum_bind   / mass_accum
!    avg_dt = dt_sum / mass_accum
!  else
!    avg_excess = 0d0
!    avg_bind   = 0d0
!    avg_dt = 0d0
!  end if
!
!  ! According to the paper, if the mass-averaged excess energy exceeds the absolute
!  ! value of the mass-averaged binding energy, then material above that point can be ejected.
!  if (avg_excess > abs(avg_bind)) then
!     ! ΔM is the envelope mass above the base of the ejection region.
!     ! In your ordering, the surface is at index 1.
!     Delta_M = s%m(1) - s%m(k_massloss)
!  else
!     Delta_M = 0d0
!  end if
!
!  if (avg_dt > 0d0) then
!     eruptive_ml = - xi * Delta_M / avg_dt
!  else
!     eruptive_ml = 0d0
!  end if
!
!  write(*,*) '- eruptive_ml = ', eruptive_ml
!  write(*,*) 'original mstar_dot = ', s%mstar_dot
!  ! Add eruptive_ml (which is negative) to the star's mass loss rate.
!  s%mstar_dot = s%mstar_dot + eruptive_ml
!  write(*,*) 'adjusted mstar_dot = ', s%mstar_dot
!
!end subroutine other_adjust_mdot




!
!
!real(dp) function eval_cell_section_total_energy( &
!      s, klo, khi) result(sum_total)
!   type (star_info), pointer :: s
!   integer, intent(in) :: klo, khi  ! sum from klo to khi
!   real(dp) :: &
!      total_internal_energy, total_gravitational_energy, &
!      total_radial_kinetic_energy, total_rotational_kinetic_energy, &
!      total_turbulent_energy
!   real(dp), allocatable, dimension(:) :: total_energy_profile
!   allocate(total_energy_profile(1:s% nz))
!   call eval_deltaM_total_energy_integrals( &
!      s, klo, khi, s% mstar, .false., &
!      total_energy_profile, &
!      total_internal_energy, total_gravitational_energy, &
!      total_radial_kinetic_energy, total_rotational_kinetic_energy, &
!      total_turbulent_energy, sum_total)
!end function eval_cell_section_total_energy
!!
!subroutine eval_deltaM_total_energy_integrals( &
!      s, klo, khi, deltaM, save_profiles, &
!      total_energy_profile, &
!      total_internal_energy, total_gravitational_energy, &
!      total_radial_kinetic_energy, total_rotational_kinetic_energy, &
!      total_turbulent_energy, sum_total)
!   type (star_info), pointer :: s
!   integer, intent(in) :: klo, khi  ! sum from klo to khi
!   real(dp), intent(in) :: deltaM
!   logical, intent(in) :: save_profiles
!   real(dp), intent(out), dimension(:) :: total_energy_profile
!   real(dp), intent(out) :: &
!      total_internal_energy, total_gravitational_energy, &
!      total_radial_kinetic_energy, total_rotational_kinetic_energy, &
!      total_turbulent_energy, sum_total
!   integer :: k
!   real(dp) :: dm, sum_dm, cell_total, cell1, d_dv00, d_dvp1, d_dlnR00, d_dlnRp1
!   include 'formats'
!
!   total_internal_energy = 0d0
!   total_gravitational_energy = 0d0
!   total_radial_kinetic_energy = 0d0
!   total_rotational_kinetic_energy = 0d0
!   total_turbulent_energy = 0d0
!   sum_total = 0d0
!
!   if (klo < 1 .or. khi > s% nz .or. klo > khi) return
!
!   sum_dm = 0
!   do k=klo,khi
!      if (sum_dm >= deltaM) exit
!      cell_total = 0
!      dm = s% dm(k)
!      if (sum_dm + dm > deltaM) dm = deltaM - sum_dm
!      cell1 = dm*s% energy(k)
!      cell_total = cell_total + cell1
!      total_internal_energy = total_internal_energy + cell1
!      if (s% v_flag .or. s% u_flag) then
!         cell1 = dm*cell_specific_KE(s,k,d_dv00,d_dvp1)
!         cell_total = cell_total + cell1
!         total_radial_kinetic_energy = total_radial_kinetic_energy + cell1
!      end if
!      cell1 = dm*cell_specific_PE(s,k,d_dlnR00,d_dlnRp1)
!      cell_total = cell_total + cell1
!      total_gravitational_energy = total_gravitational_energy + cell1
!      if (s% rotation_flag) then
!         cell1 = dm*cell_specific_rotational_energy(s,k)
!         total_rotational_kinetic_energy = total_rotational_kinetic_energy + cell1
!         if (s% include_rotation_in_total_energy) &
!            cell_total = cell_total + cell1
!      end if
!      if (s% RSP2_flag) then
!         cell1 = dm*pow2(s% w(k))
!         cell_total = cell_total + cell1
!         total_turbulent_energy = total_turbulent_energy + cell1
!      end if
!      if (s% rsp_flag) then
!         cell1 = dm*s% RSP_Et(k)
!         cell_total = cell_total + cell1
!         total_turbulent_energy = total_turbulent_energy + cell1
!      end if
!      if (save_profiles) then
!         total_energy_profile(k) = cell_total
!      end if
!   end do
!
!   sum_total = total_gravitational_energy + &
!      total_radial_kinetic_energy + total_turbulent_energy ! + total_internal_energy
!
!   if (s% include_rotation_in_total_energy) &
!      sum_total = sum_total + total_rotational_kinetic_energy
!
!end subroutine eval_deltaM_total_energy_integrals
!
!real(dp) function cell_specific_PE(s,k,d_dlnR00,d_dlnRp1)
!   type (star_info), pointer :: s
!   integer, intent(in) :: k
!   real(dp), intent(out) :: d_dlnR00,d_dlnRp1
!   cell_specific_PE = cell_specific_PE_qp(s,k,d_dlnR00,d_dlnRp1)
!end function cell_specific_PE



      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         integer :: k
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)
         do k = 1,s% nz
         s% xtra1_array(k) = 0d0
         end do
         if (.not. s% x_logical_ctrl(37)) return

      end subroutine extras_startup


      subroutine extras_after_evolve(id, ierr)
         use num_lib, only: find0
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt, m
         integer :: k, nz
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         nz = s% nz
         write(*,'(A)')
         select case (s% x_integer_ctrl(5))
         case (7)
            ! put target info in TestHub output
            testhub_extras_names(1) = 'fe_core_mass'; testhub_extras_vals(1) = s% fe_core_mass

            if(s% fe_core_mass < 1d0) then
               write(*,1) "Bad fe_core_mass", s%fe_core_mass
            else
               if(s% fe_core_infall > s% fe_core_infall_limit) then
                  write(*,'(a)') 'all values are within tolerance'
               else
                  write(*,'(a)') "Bad fe core infall"
               end if
            end if
         end select
         call test_suite_after_evolve(s, ierr)
         if (.not. s% x_logical_ctrl(37)) return
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
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         if (.not. s% x_logical_ctrl(37)) return
         if (extras_finish_step == terminate) &
             s% termination_code = t_extras_finish_step
      end function extras_finish_step


!
!real(dp) function cell_start_specific_KE(s,k)
!   type (star_info), pointer :: s
!   integer, intent(in) :: k
!   cell_start_specific_KE = cell_start_specific_KE_qp(s,k)
!end function cell_start_specific_KE
!
!
!real(qp) function cell_start_specific_KE_qp(s,k)
!   ! for consistency with dual cells at faces, use <v**2> instead of <v>**2
!   type (star_info), pointer :: s
!   integer, intent(in) :: k
!   real(qp) :: qhalf, v0, v1, v2, Mbar
!   qhalf = 0.5d0
!   if (s% use_mass_corrections) then
!      Mbar = s% mass_correction(k)
!   else
!      Mbar = 1.0_qp
!   end if
!   if (s% u_flag) then
!      v0 = s% u_start(k)
!      cell_start_specific_KE_qp = qhalf*Mbar*v0**2
!   else if (s% v_flag) then
!      v0 = s% v_start(k)
!      if (k < s% nz) then
!         v1 = s% v_start(k+1)
!      else
!         v1 = s% v_center
!      end if
!      v2 = qhalf*(v0**2 + v1**2)
!      cell_start_specific_KE_qp = qhalf*Mbar*v2
!   else  ! ignore kinetic energy if no velocity variables
!      cell_start_specific_KE_qp = 0d0
!   end if
!end function cell_start_specific_KE_qp
!
!
!real(dp) function cell_specific_KE(s,k,d_dv00,d_dvp1)
!   type (star_info), pointer :: s
!   integer, intent(in) :: k
!   real(dp), intent(out) :: d_dv00, d_dvp1
!   cell_specific_KE = cell_specific_KE_qp(s,k,d_dv00,d_dvp1)
!end function cell_specific_KE
!
!
!real(qp) function cell_specific_KE_qp(s,k,d_dv00,d_dvp1)
!   ! for consistency with dual cells at faces, use <v**2> instead of <v>**2
!   type (star_info), pointer :: s
!   integer, intent(in) :: k
!   real(dp), intent(out) :: d_dv00, d_dvp1
!   real(dp) :: dv2_dv00, dv2_dvp1
!   real(qp) :: qhalf, v0, v1, v2, Mbar
!   qhalf = 0.5d0
!   if (s% use_mass_corrections) then
!      Mbar = s% mass_correction(k)
!   else
!      Mbar = 1.0_qp
!   end if
!   if (s% u_flag) then
!      v0 = s% u(k)
!      cell_specific_KE_qp = qhalf*Mbar*v0**2
!      d_dv00 = s% u(k)
!      d_dvp1 = 0d0
!   else if (s% v_flag) then
!      v0 = s% v(k)
!      if (k < s% nz) then
!         v1 = s% v(k+1)
!         dv2_dvp1 = s% v(k+1)
!      else
!         v1 = s% v_center
!         dv2_dvp1 = 0d0
!      end if
!      v2 = qhalf*(v0**2 + v1**2)
!      dv2_dv00 = s% v(k)
!      cell_specific_KE_qp = qhalf*Mbar*v2
!      d_dv00 = qhalf*Mbar*dv2_dv00
!      d_dvp1 = qhalf*Mbar*dv2_dvp1
!   else  ! ignore kinetic energy if no velocity variables
!      cell_specific_KE_qp = 0d0
!      d_dv00 = 0d0
!      d_dvp1 = 0d0
!   end if
!end function cell_specific_KE_qp
!
!
!real(dp) function cell_specific_PE(s,k,d_dlnR00,d_dlnRp1)
!   type (star_info), pointer :: s
!   integer, intent(in) :: k
!   real(dp), intent(out) :: d_dlnR00,d_dlnRp1
!   cell_specific_PE = cell_specific_PE_qp(s,k,d_dlnR00,d_dlnRp1)
!end function cell_specific_PE
!
!
!real(qp) function cell_specific_PE_qp(s,k,d_dlnR00,d_dlnRp1)
!   ! for consistency with dual cells at faces, <m/r>_cntr => (m(k)/r(k) + m(k+1)/r(k+1))/2 /= m_cntr/r_cntr
!   ! i.e., use avg of m/r at faces of cell rather than ratio of cell center mass over cell center r.
!   type (star_info), pointer :: s
!   integer, intent(in) :: k
!   real(dp), intent(out) :: d_dlnR00,d_dlnRp1
!   real(qp) :: qhalf, rp1, r00, mp1, m00, Gp1, G00, gravp1, grav00, Mbar
!   real(dp) :: d_grav00_dlnR00, d_gravp1_dlnRp1
!   include 'formats'
!   qhalf = 0.5d0
!   if (s% use_mass_corrections) then
!      Mbar = s% mass_correction(k)
!   else
!      Mbar = 1.0_qp
!   end if
!   if (k == s% nz) then
!      rp1 = s% R_center
!      mp1 = s% m_center
!      Gp1 = s% cgrav(s% nz)
!   else
!      rp1 = s% r(k+1)
!      mp1 = s% m_grav(k+1)
!      Gp1 = s% cgrav(k+1)
!   end if
!   if (rp1 <= 0d0) then
!      gravp1 = 0d0
!      d_gravp1_dlnRp1 = 0d0
!   else
!      gravp1 = -Gp1*mp1/rp1
!      d_gravp1_dlnRp1 = -gravp1
!   end if
!   r00 = s% r(k)
!   m00 = s% m_grav(k)
!   G00 = s% cgrav(k)
!   grav00 = -G00*m00/r00
!   d_grav00_dlnR00 = -grav00
!   cell_specific_PE_qp = qhalf*Mbar*(gravp1 + grav00)
!   d_dlnR00 = qhalf*Mbar*d_grav00_dlnR00
!   d_dlnRp1 = qhalf*Mbar*d_gravp1_dlnRp1
!   if (is_bad(cell_specific_PE_qp)) then
!      write(*,2) 'cell_specific_PE_qp', k, cell_specific_PE_qp
!      write(*,2) 'gravp1', k, gravp1
!      write(*,2) 'grav00', k, grav00
!      call mesa_error(__FILE__,__LINE__,'cell_specific_PE')
!   end if
!end function cell_specific_PE_qp
!
!
!real(dp) function cell_start_specific_PE(s,k)
!   type (star_info), pointer :: s
!   integer, intent(in) :: k
!   cell_start_specific_PE = cell_start_specific_PE_qp(s,k)
!end function cell_start_specific_PE
!
!
!real(dp) function cell_start_specific_PE_qp(s,k)
!   ! for consistency with dual cells at faces, <m/r>_cntr => (m(k)/r(k) + m(k+1)/r(k+1))/2 /= m_cntr/r_cntr
!   ! i.e., use avg of m/r at faces of cell rather than ratio of cell center mass over cell center r.
!   type (star_info), pointer :: s
!   integer, intent(in) :: k
!   real(qp) :: qhalf, rp1, r00, mp1, m00, Gp1, G00, gravp1, grav00, Mbar
!   include 'formats'
!   qhalf = 0.5d0
!   if (s% use_mass_corrections) then
!      Mbar = s% mass_correction_start(k)
!   else
!      Mbar = 1.0_qp
!   end if
!   if (k == s% nz) then
!      rp1 = s% R_center
!      mp1 = s% m_center
!      Gp1 = s% cgrav(s% nz)
!   else
!      rp1 = s% r_start(k+1)
!      mp1 = s% m_grav_start(k+1)
!      Gp1 = s% cgrav(k+1)
!   end if
!   if (rp1 <= 0d0) then
!      gravp1 = 0d0
!   else
!      gravp1 = -Gp1*mp1/rp1
!   end if
!   r00 = s% r_start(k)
!   m00 = s% m_grav_start(k)
!   G00 = s% cgrav(k)
!   grav00 = -G00*m00/r00
!   cell_start_specific_PE_qp = qhalf*Mbar*(gravp1 + grav00)
!   if (is_bad(cell_start_specific_PE_qp)) then
!      write(*,2) 'cell_start_specific_PE_qp', k, cell_start_specific_PE_qp
!      write(*,2) 'gravp1', k, gravp1
!      write(*,2) 'grav00', k, grav00
!      call mesa_error(__FILE__,__LINE__,'cell_start_specific_PE_qp')
!   end if
!end function cell_start_specific_PE_qp
!
!
!real(dp) function cell_specific_rotational_energy(s,k)
!   type (star_info), pointer :: s
!   integer, intent(in) :: k
!   real(dp) :: e_00, e_p1
!   e_00 = s% i_rot(k)% val*s% omega(k)*s% omega(k)
!   if (k < s% nz) then
!      e_p1 = s% i_rot(k+1)% val*s% omega(k+1)*s% omega(k+1)
!   else
!      e_p1 = 0
!   end if
!   cell_specific_rotational_energy = 0.5d0*(e_p1 + e_00)
!end function cell_specific_rotational_energy
!
!
!subroutine get_dke_dt_dpe_dt(s, k, dt, &
!      dke_dt, d_dkedt_dv00, d_dkedt_dvp1, &
!      dpe_dt, d_dpedt_dlnR00, d_dpedt_dlnRp1, ierr)
!   type (star_info), pointer :: s
!   integer, intent(in) :: k
!   real(dp), intent(in) :: dt
!   real(dp), intent(out) :: &
!      dke_dt, d_dkedt_dv00, d_dkedt_dvp1, &
!      dpe_dt, d_dpedt_dlnR00, d_dpedt_dlnRp1
!   integer, intent(out) :: ierr
!   real(dp) :: PE_start, PE_new, KE_start, KE_new, q1
!   real(dp) :: dpe_dlnR00, dpe_dlnRp1, dke_dv00, dke_dvp1
!   include 'formats'
!   ierr = 0
!   ! rate of change in specific PE (erg/g/s)
!   PE_start = cell_start_specific_PE_qp(s,k)
!   PE_new = cell_specific_PE_qp(s,k,dpe_dlnR00,dpe_dlnRp1)
!   q1 = PE_new - PE_start
!   dpe_dt = q1/dt  ! erg/g/s
!   d_dpedt_dlnR00 = dpe_dlnR00/dt
!   d_dpedt_dlnRp1 = dpe_dlnRp1/dt
!   ! rate of change in specific KE (erg/g/s)
!   KE_start = cell_start_specific_KE_qp(s,k)
!   KE_new = cell_specific_KE_qp(s,k,dke_dv00,dke_dvp1)
!   q1 = KE_new - KE_start
!   dke_dt = q1/dt  ! erg/g/s
!   d_dkedt_dv00 = dke_dv00/dt
!   d_dkedt_dvp1 = dke_dvp1/dt
!end subroutine get_dke_dt_dpe_dt
!

      end module run_star_extras

