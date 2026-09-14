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

module run_star_extras

  use star_lib
  use star_def
  use const_def
  use math_lib
  use chem_def
  use num_lib
  use binary_def
  use ionization_def

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
    s% other_wind => my_wind

  end subroutine extras_controls


  subroutine extras_startup(id, restart, ierr)
    integer, intent(in) :: id
    logical, intent(in) :: restart
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    ! use 11-16 to avoid overlap with lxtra(1) and lxtra(2) used in run_binary_extras.f90
    if (.not. restart) then
       s% lxtra(11) = .false.
       s% lxtra(12) = .false.
       s% lxtra(13) = .false.
       s% lxtra(14) = .false.
       s% lxtra(15) = .false.
       s% lxtra(16) = .false.
       s% lxtra(17) = .false.
       s% lxtra(18) = .false.
       s% lxtra(19) = .false.
       s% lxtra(20) = .false.
       s% lxtra(21) = .false.
    end if
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


  ! returns either keep_going, retry, backup, or terminate.
  integer function extras_check_model(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    real(dp) :: error, atol, rtol
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_check_model = keep_going

    ! check we don't overshoot too much the radii of interest
    if ((s% lxtra(11) .eqv. .false.) .and. (s% r(1)/Rsun >= 20.0d0)) then

        error = 1d0*(s%r(1)/Rsun - 20.0d0)

        if (error > 0.1d0) then
            extras_check_model = retry
        end if

    end if

    !check we don't overshoot too much the radii of interest
    if ((s% lxtra(12) .eqv. .false.) .and. (s% r(1)/Rsun >= 25.0d0)) then

        error = 1d0*(s%r(1)/Rsun - 25.0d0)

        if (error > 1d0) then
            print*, "error", error
            extras_check_model = retry
        end if

    end if

    !check we don't overshoot too much the radii of interest
    if ((s% lxtra(13) .eqv. .false.) .and. (s% r(1)/Rsun >= 30.0d0)) then

        error = 1d0*(s%r(1)/Rsun - 30.0d0)

        if (error > 1d0) then
            print*, "error", error
            extras_check_model = retry

        end if

    end if

    !check we don't overshoot too much the radii of interest
    if ((s% lxtra(14) .eqv. .false.) .and. (s% r(1)/Rsun >= 50.0d0)) then

        error = 1d0*(s%r(1)/Rsun - 50.0d0)

        if (error > 1d0) then
            print*, "error", error
            extras_check_model = retry

        end if

    end if

    !check we don't overshoot too much the radii of interest
    if ((s% lxtra(15) .eqv. .false.) .and. (s% r(1)/Rsun >= 100.0d0)) then

        error = 1d0*(s%r(1)/Rsun - 100.0d0)

        if (error > 5d0) then
            print*, "error", error
            extras_check_model = retry

        end if

    end if

    !check we don't overshoot too much the radii of interest
    if ((s% lxtra(16) .eqv. .false.) .and. (s% r(1)/Rsun >= 500.0d0)) then

        error = 1d0*(s%r(1)/Rsun - 500.0d0)

        if (error > 10d0) then
            print*, "error", error
            extras_check_model = retry

        end if

    end if

    !check we don't overshoot too much the radii of interest
    if ((s% lxtra(17) .eqv. .false.) .and. (s% r(1)/Rsun >= 1000.0d0)) then

        error = 1d0*(s%r(1)/Rsun - 1000.0d0)

        if (error > 20d0) then
            print*, "error", error
            extras_check_model = retry

        end if

    end if

    !check we don't overshoot too much the radii of interest
    if ((s% lxtra(18) .eqv. .false.) .and. (s% r(1)/Rsun >= 1200.0d0)) then

        error = 1d0*(s%r(1)/Rsun - 1200.0d0)

        if (error > 20d0) then
            print*, "error", error
            extras_check_model = retry

        end if

    end if

    !check we don't overshoot too much the radii of interest
    if ((s% lxtra(19) .eqv. .false.) .and. (s% r(1)/Rsun >= 1300.0d0)) then

        error = 1d0*(s%r(1)/Rsun - 1300.0d0)

        if (error > 20d0) then
            print*, "error", error
            extras_check_model = retry

        end if

    end if

    !check we don't overshoot too much the radii of interest
    if ((s% lxtra(20) .eqv. .false.) .and. (s% r(1)/Rsun >= 1400.0d0)) then

        error = 1d0*(s%r(1)/Rsun - 1400.0d0)

        if (error > 20d0) then
            print*, "error", error
            extras_check_model = retry

        end if

    end if

    !check we don't overshoot too much the radii of interest
    if ((s% lxtra(21) .eqv. .false.) .and. (s% r(1)/Rsun >= 1500.0d0)) then


        error = 1d0*(s%r(1)/Rsun - 1500.0d0)

        if (error > 20d0) then
            print*, "error", error
            extras_check_model = retry

        end if

    end if

    ! by default, indicate where (in the code) MESA terminated
    if (extras_check_model == terminate) s% termination_code = t_extras_check_model
  end function extras_check_model


  integer function how_many_extra_history_columns(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    how_many_extra_history_columns =  13
  end function how_many_extra_history_columns


  subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
    use chem_def, only: chem_isos
    integer, intent(in) :: id, n
    character (len=maxlen_history_column_name) :: names(n)
    real(dp) :: vals(n)
    integer, intent(out) :: ierr
    ! POSYDON output
    real(dp) :: env_binding_E, total_env_binding_E
    integer :: i, k, nz
    integer :: i1, k1, k2, j
    integer :: h1, he4, c12, o16
    real(dp) :: he_core_mass_1cent,  he_core_mass_10cent, he_core_mass_30cent
    real(dp) :: he_core_radius_1cent, he_core_radius_10cent, he_core_radius_30cent
    real(dp) ::  lambda_CE_1cent, lambda_CE_10cent, lambda_CE_30cent
    real(dp),  dimension (:), allocatable ::  adjusted_energy, energy
    real(dp) :: rec_energy_HII_to_HI, &
         rec_energy_HeII_to_HeI, &
         rec_energy_HeIII_to_HeII, &
         diss_energy_H2, &
         frac_HI, frac_HII, &
         frac_HeI, frac_HeII, frac_HeIII, &
         avg_charge_He, energy_comp
    logical :: have_30_value, have_10_value, have_1_value, have_co_value
    logical :: sticking_to_energy_without_recombination_corr
    ! -------------------------------------
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    nz = s% nz

    ! output info about the ENV.: binding energy
    total_env_binding_E = 0.0d0
    do k=1,s% nz-1
       if (s% m(k) > (s% he_core_mass * Msun)) then !envelope is defined to be H-rich
          env_binding_E = s% dm(k) * (s% energy(k) - (s% cgrav(k) * s% m(k+1))/s% r(k+1))
          total_env_binding_E = total_env_binding_E + env_binding_E
       end if
    end do

    names(1) = 'envelope_binding_energy'
    vals(1) = total_env_binding_E

    ! lambda_CE calculation for different core definitions.
    ! get energy from the EOS and adjust the different contributions from recombination/dissociation to internal energy
    ! allocate(adjusted_energy(s% nz))
    allocate(energy(s% nz))
    ! adjusted_energy=0.0d0
    energy=0.0d0

    ! do k=1, s% nz
    !   ! the following lines compute the fractions of HI, HII, HeI, HeII and HeIII
    !   ! things like ion_ifneut_H are defined in $MESA_DIR/ionization/public/ionization.def
    !   ! this file can be checked for additional ionization output available
    !   frac_HI = get_ion_info(s,ion_ifneut_H,k)
    !   frac_HII = 1.0d0 - frac_HI

    !   ! ionization module provides neutral fraction and average charge of He.
    !   ! use these two to compute the mass fractions of HeI and HeII
    !   frac_HeI = get_ion_info(s,ion_ifneut_He,k)
    !   avg_charge_He = get_ion_info(s,ion_iZ_He,k)
    !   ! the following is the solution to the equations
    !   !   avg_charge_He = 2*fracHeIII + 1*fracHeII
    !   !               1 = fracHeI + fracHeII + fracHeIII
    !   frac_HeII = 2d0 - 2d0*frac_HeI - avg_charge_He
    !   frac_HeIII = 1d0 - frac_HeII - frac_HeI

    !   ! recombination energies from https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
    !   rec_energy_HII_to_HI = avo*13.59843449d0*frac_HII*ev2erg*s% X(k)
    !   diss_energy_H2 = avo*4.52d0/2d0*ev2erg*s% X(k)
    !   rec_energy_HeII_to_HeI = avo*24.58738880d0*(frac_HeII+frac_HeIII)*ev2erg*s% Y(k)/4d0
    !   rec_energy_HeIII_to_HeII = avo*54.4177650d0*frac_HeIII*ev2erg*s% Y(k)/4d0

    !   adjusted_energy(k) = s% energy(k) &
    !         - rec_energy_HII_to_HI &
    !         - rec_energy_HeII_to_HeI &
    !         - rec_energy_HeIII_to_HeII &
    !         - diss_energy_H2

    !     !write(*,*) "s% energy(k):", s% energy(k), " adjusted_energy, ", adjusted_energy(k)
    !     !write(*,*) "frac HII", frac_HII
    !   if (adjusted_energy(k) < 0d0 .or. adjusted_energy(k) > s% energy(k)) then
    !       write(*,*) "Error when computing adjusted energy in CE, ", &
    !           "s% energy(k):", s% energy(k), " adjusted_energy, ", adjusted_energy(k)
    !       !sticking_to_energy_without_recombination_corr = .true.
    !   end if

    !   if(.false.) then
    !       ! for debug, check the mismatch between the EOS energy and that of a gas+radiation
    !       energy_comp = 3.0d0*avo*boltzm*s% T(k)/(2*s% mu(k)) + crad*pow4(s% T(k))/s% rho(k) &
    !           + rec_energy_HII_to_HI &
    !           + rec_energy_HeII_to_HeI &
    !           + rec_energy_HeIII_to_HeII &
    !           + diss_energy_H2

    !       write(*,*) "compare energies", k, s%m(k)/Msun, s% energy(k), energy_comp, &
    !           (s% energy(k)-energy_comp)/s% energy(k)
    !   end if

    ! end do

   do k=1, s% nz
      energy(k) = s% energy(k)
   end do

    ! to do. lambda_CEE calculation for He star envelope too.s
    he_core_mass_1cent = 0.0d0
    he_core_mass_10cent = 0.0d0
    he_core_mass_30cent = 0.0d0
    !for MS stars = he core is 0, lambda calculated for whole star
    !for He stars = he core is whole star, lambda calculated for whole star

    he_core_radius_1cent = 0.0d0
    he_core_radius_10cent = 0.0d0
    he_core_radius_30cent = 0.0d0

    h1 = s% net_iso(ih1)
    he4 = s% net_iso(ihe4)
    have_30_value = .false.
    have_10_value = .false.
    have_1_value = .false.
    if (h1 /= 0 .and. he4 /= 0) then
       do k=1, s% nz
          if (.not. have_30_value) then
             if (s% xa(h1,k) <=  0.3d0 .and. &
                  s% xa(he4,k) >= 0.1d0) then
                he_core_mass_30cent = s% m(k)
                he_core_radius_30cent = s% r(k)
                have_30_value = .true.
             end if
          end if
          if (.not. have_10_value) then
             if (s% xa(h1,k) <= 0.1d0 .and. &
                  s% xa(he4,k) >= 0.1d0) then
                he_core_mass_10cent = s% m(k)
                he_core_radius_10cent = s% r(k)
                have_10_value = .true.
             end if
          end if
          if (.not. have_1_value) then
             if (s% xa(h1,k) <= 0.01d0 .and. &
                  s% xa(he4,k) >= 0.1d0) then
                he_core_mass_1cent = s% m(k)
                he_core_radius_1cent = s% r(k)
                have_1_value = .true.
             end if
          end if
       end do
    end if

    lambda_CE_1cent = lambda_CE(s,energy, he_core_mass_1cent)
    lambda_CE_10cent = lambda_CE(s,energy, he_core_mass_10cent)
    lambda_CE_30cent = lambda_CE(s,energy, he_core_mass_30cent)

    names(2) = 'lambda_CE_1cent'
    vals(2) = lambda_CE_1cent
    names(3) = 'lambda_CE_10cent'
    vals(3) = lambda_CE_10cent
    names(4) = 'lambda_CE_30cent'
    vals(4) = lambda_CE_30cent

    names(5) = 'he_core_mass_1cent'
    vals(5) = he_core_mass_1cent
    names(6) = 'he_core_mass_10cent'
    vals(6) = he_core_mass_10cent
    names(7) = 'he_core_mass_30cent'
    vals(7) = he_core_mass_30cent
    names(8) = 'he_core_radius_1cent'
    vals(8) = he_core_radius_1cent
    names(9) = 'he_core_radius_10cent'
    vals(9) = he_core_radius_10cent
    names(10) = 'he_core_radius_30cent'
    vals(10) = he_core_radius_30cent

    do k=1, s% nz
      energy(k) = 0
    end do

    lambda_CE_1cent = lambda_CE(s,energy, he_core_mass_1cent)
    lambda_CE_10cent = lambda_CE(s,energy, he_core_mass_10cent)
    lambda_CE_30cent = lambda_CE(s,energy, he_core_mass_30cent)

    names(11) = 'lambda_CE_1cent_grav'
    vals(11) = lambda_CE_1cent
    names(12) = 'lambda_CE_10cent_grav'
    vals(12) = lambda_CE_10cent
    names(13) = 'lambda_CE_30cent_grav'
    vals(13) = lambda_CE_30cent

  end subroutine data_for_extra_history_columns


  integer function how_many_extra_profile_columns(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    how_many_extra_profile_columns = 0
  end function how_many_extra_profile_columns


  subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
    integer, intent(in) :: id, n, nz
    character (len=maxlen_profile_column_name) :: names(n)
    real(dp) :: vals(nz,n)
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    integer :: k

    ! real(dp) :: rec_energy_HII_to_HI, &
    !      rec_energy_HeII_to_HeI, &
    !      rec_energy_HeIII_to_HeII, &
    !      diss_energy_H2, &
    !      frac_HI, frac_HII, &
    !      frac_HeI, frac_HeII, frac_HeIII, &
    !      avg_charge_He, energy_comp

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    ! note: do NOT add the extra names to profile_columns.list
    ! the profile_columns.list is only for the built-in profile column options.
    ! it must not include the new column names you are adding here.

    ! here is an example for adding a profile column
    !if (n /= 1) stop 'data_for_extra_profile_columns'
    !names(1) = 'beta'
    !do k = 1, nz
    !   vals(k,1) = s% Pgas(k)/s% P(k)
    !end do

    ! names(1) = 'u_no_rec'
    ! names(2) = 'u_yes_rec'

    ! do k=1, s% nz
    !   ! the following lines compute the fractions of HI, HII, HeI, HeII and HeIII
    !   ! things like ion_ifneut_H are defined in $MESA_DIR/ionization/public/ionization.def
    !   ! this file can be checked for additional ionization output available
    !   frac_HI = get_ion_info(s,ion_ifneut_H,k)
    !   frac_HII = 1.0d0 - frac_HI

    !   ! ionization module provides neutral fraction and average charge of He.
    !   ! use these two to compute the mass fractions of HeI and HeII
    !   frac_HeI = get_ion_info(s,ion_ifneut_He,k)
    !   avg_charge_He = get_ion_info(s,ion_iZ_He,k)
    !   ! the following is the solution to the equations
    !   !   avg_charge_He = 2*fracHeIII + 1*fracHeII
    !   !               1 = fracHeI + fracHeII + fracHeIII
    !   frac_HeII = 2d0 - 2d0*frac_HeI - avg_charge_He
    !   frac_HeIII = 1d0 - frac_HeII - frac_HeI

    !   ! recombination energies from https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
    !   rec_energy_HII_to_HI = avo*13.59843449d0*frac_HII*ev2erg*s% X(k)
    !   diss_energy_H2 = avo*4.52d0/2d0*ev2erg*s% X(k)
    !   rec_energy_HeII_to_HeI = avo*24.58738880d0*(frac_HeII+frac_HeIII)*ev2erg*s% Y(k)/4d0
    !   rec_energy_HeIII_to_HeII = avo*54.4177650d0*frac_HeIII*ev2erg*s% Y(k)/4d0

    !     vals(k,1) = s% energy(k) &
    !         - rec_energy_HII_to_HI &
    !         - rec_energy_HeII_to_HeI &
    !         - rec_energy_HeIII_to_HeII &
    !         - diss_energy_H2

    !     vals(k,2) = s% energy(k)
    ! end do

  end subroutine data_for_extra_profile_columns


  integer function how_many_extra_history_header_items(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    how_many_extra_history_header_items = 3
  end function how_many_extra_history_header_items


  subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
    integer, intent(in) :: id, n
    character (len=maxlen_history_column_name) :: names(n)
    real(dp) :: vals(n)
    type(star_info), pointer :: s
    integer, intent(out) :: ierr
    integer :: i
    real(dp) :: Initial_X, Initial_Y, Initial_Z, initial_m

    ierr = 0
    call star_ptr(id,s,ierr)
    if(ierr/=0) return

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
  ! note: cannot request retry or backup; extras_check_model can do that.
  integer function extras_finish_step(id)
    integer, intent(in) :: id
    integer :: ierr, i
    real(dp) :: envelope_mass_fraction, L_He, L_tot, min_center_h1_for_diff, &
         critmass, feh, rot_full_off, rot_full_on, frac2, TAMS_h1_treshold
    real(dp), parameter :: huge_dt_limit = 3.15d16 ! ~1 Gyr
    real(dp), parameter :: new_varcontrol_target = 1d-3
    real(dp), parameter :: Zsol = 0.0142_dp
    logical :: diff_test1, diff_test2, diff_test3, is_ne_biggest
    !character (len=strlen) :: photoname, fname
    character (len=30) :: photoname, fname
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_finish_step = keep_going

    ! write(fname, fmt="(a11)") 'recombination.data'
    ! print*, "saving profile "// fname
    ! call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
    ! extras_finish_step = terminate

    if (s%lxtra(11) .eqv. .false.) then
       ! save profile for R=20Rsun
       if (s% r(1)/Rsun >= 20) then
          s% lxtra(11) = .true.
          write(fname, fmt="(a11)") '20Rsun.data'
          print*, "saving profile "// fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          write(fname, fmt="(a10)") '20Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          write(fname, fmt="(a6)") '20Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if

    if (s%lxtra(12) .eqv. .false.) then
       ! save profile for R=25Rsun
       if (s% r(1)/Rsun >= 25) then
          s% lxtra(12) = .true.
          fname = '25Rsun.data'
          print*, "saving profile "// fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          fname = '25Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          fname = '25Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if

    if (s%lxtra(13) .eqv. .false.) then
       ! save profile for R=30Rsun
       if (s% r(1)/Rsun >= 30) then
          s% lxtra(13) = .true.
          fname = '30Rsun.data'
          print*, "saving profile "// fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          fname = '30Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          fname = '30Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if

    if (s%lxtra(14) .eqv. .false.) then
       ! save profile for R=50Rsun
       if (s% r(1)/Rsun >= 50) then
          s% lxtra(14) = .true.
          fname = '50Rsun.data'
          print*, "saving profile "//fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          fname = '50Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          fname = '50Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if

    if (s%lxtra(15) .eqv. .false.) then
       ! save profile for R=100Rsun
       if (s% r(1)/Rsun >= 100) then
          s% lxtra(15) = .true.
          write(fname, fmt="(a12)") '100Rsun.data'
          print*, "saving profile "//fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          write(fname, fmt="(a11)") '100Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          write(fname, fmt="(a7)") '100Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if

    if (s%lxtra(16) .eqv. .false.) then
       ! save profile for R=500Rsun
       if (s% r(1)/Rsun >= 500) then
          s% lxtra(16) = .true.
          fname = '500Rsun.data'
          print*, "saving profile "//fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          fname = '500Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          fname = '500Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if

    if (s%lxtra(17) .eqv. .false.) then
       ! save profile for R=1000Rsun
       if (s% r(1)/Rsun >= 1000) then
          s% lxtra(17) = .true.
          write(fname, fmt="(a13)") '1000Rsun.data'
          print*, "saving profile "//fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          write(fname, fmt="(a12)") '1000Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          write(fname, fmt="(a8)") '1000Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if

    if (s%lxtra(18) .eqv. .false.) then
       ! save profile for R=1200Rsun
       if (s% r(1)/Rsun >= 1200) then
          s% lxtra(18) = .true.
          write(fname, fmt="(a13)") '1200Rsun.data'
          print*, "saving profile "//fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          write(fname, fmt="(a12)") '1200Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          write(fname, fmt="(a8)") '1200Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if

    if (s%lxtra(19) .eqv. .false.) then
       ! save profile for R=1300Rsun
       if (s% r(1)/Rsun >= 1300) then
          s% lxtra(19) = .true.
          write(fname, fmt="(a13)") '1300Rsun.data'
          print*, "saving profile "//fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          write(fname, fmt="(a12)") '1300Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          write(fname, fmt="(a8)") '1300Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if

    if (s%lxtra(20) .eqv. .false.) then
       ! save profile for R=1400Rsun
       if (s% r(1)/Rsun >= 1400) then
          s% lxtra(20) = .true.
          write(fname, fmt="(a13)") '1400Rsun.data'
          print*, "saving profile "//fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          write(fname, fmt="(a12)") '1400Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          write(fname, fmt="(a8)") '1400Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if

    if (s%lxtra(21) .eqv. .false.) then
       ! save profile for R=1500Rsun
       if (s% r(1)/Rsun >= 1500) then
          s% lxtra(21) = .true.
          write(fname, fmt="(a13)") '1500Rsun.data'
          print*, "saving profile "//fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          write(fname, fmt="(a12)") '1500Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          write(fname, fmt="(a8)") '1500Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if

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


    subroutine my_wind(id, Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z, w, ierr)
   use star_def
   type (star_info), pointer :: s
   integer, intent(in) :: id
   real(dp), parameter :: Zbjork = 0.014d0, T_high = 11000, T_low = 10000
   real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z ! surface values (cgs)
   real(dp), intent(out) :: w ! wind in units of Msun/year (value is >= 0)
   integer, intent(out) :: ierr
   real(dp) :: M1, L1, T1, R1, w1, w2, alfa
   w = 0
   ierr = 0
   call star_ptr(id, s, ierr)

   M1 = Msurf
   L1 = Lsurf
   T1 = Tsurf

   if (T1 <= T_low) then
      call eval_lowT(w)
      !print*, "Decin"

   else if (T1 >= T_high) then
        call eval_highT(w)

   else ! transition
      call eval_lowT(w1)
      call eval_highT(w2)
      alfa = (T1 - T_low) / (T_high - T_low)
      w = (1-alfa)*w1 + alfa*w2
   end if

contains

    subroutine eval_de_Jager_wind(w4)
       real(dp), intent(out) :: w4
       real(dp) :: log10w

       log10w = 1.769d0*log10(L1/Lsun) - 1.676d0*log10(T1) - 8.158d0
       w4 = exp10(log10w)
   end subroutine eval_de_Jager_wind

   subroutine eval_Decin_wind(w)
      ! Nieuwenhuijzen, H.; de Jager, C. 1990, A&A, 231, 134 (eqn 2)
      real(dp), intent(out) :: w
      real(dp) :: log10w
      include 'formats'
      !log10w = -17.44 &
      !         +2.80*log10(L1/Lsun) &
      !         -0.13*s% initial_mass! &
      log10w = -20.71 &
               +3.50*log10(L1/Lsun) &
               -0.17*s% initial_mass+0.49! &
      !write(*,*) "Check Decin wind", log10w
      w = 10**(log10w)
   end subroutine eval_Decin_wind

  subroutine eval_lowT(w)
        real(dp), intent(out) :: w
        if (s% x_character_ctrl(2) == 'Decin') then
                call eval_Decin_wind(w)
        end if

        if (s% x_character_ctrl(2) == 'nowind') then
                call eval_Decin_wind(w1)
                w = min(w1, 1d-6)
        end if

   end subroutine eval_lowT

   subroutine eval_highT(w5)
      real(dp), intent(out) :: w5

      if (X < 0.4d0) then
         w5 = 1d-11 * pow(L1/Lsun, 1.29d0) * pow(Y, 1.7d0) * sqrt(Z)
         !print*, "Nugis & Lamers"
      else if (s% x_character_ctrl(1) == 'Bjorklund') then
         call eval_Bjorklund_wind(w5)
      else if (s% x_character_ctrl(1) == 'Vink') then
         call eval_Vink_wind(w5)
      else if (s% x_character_ctrl(1) == 'Sabhahit') then
        call eval_Sabhahit_wind(w5)
      end if
   end subroutine eval_highT

   subroutine eval_Bjorklund_wind(w6)
      real(dp), intent(out) :: w6
      real(dp) :: logw, Meff

      Meff = M1 * (1d0 - 0.34d0 * L1 / (pi4 * clight * s% cgrav(1) * M1)) ! effective mass

      logw = -5.52d0 &
           + 2.39d0 * log10(L1 / (1d6 * Lsun)) &
           - 1.48d0 * log10(Meff / (4.5d1 * Msun)) &
           + 2.12d0 * log10(T1 / 4.5d4) &
           + (0.75d0 - 1.87d0 * log10(T1 / 4.5d4)) * log10(Z / Zbjork)
      w6 = exp10(logw)
   end subroutine eval_Bjorklund_wind

    subroutine eval_Vink_wind(w)
    real(dp), intent(inout) :: w
    real(dp) :: alfa, w1, w2, Teff_jump, logMdot, dT, vinf_div_vesc
    real(dp), parameter :: Zsolar = 0.019d0 ! for Vink et al formula

    ! alfa = 1 for hot side, = 0 for cool side
    if (T1 > 27500d0) then
       alfa = 1
    else if (T1 < 22500d0) then
       alfa = 0
    else ! use Vink et al 2001, eqns 14 and 15 to set "jump" temperature
       Teff_jump = 1d3*(61.2d0 + 2.59d0*(-13.636d0 + 0.889d0*log10(Z/Zsolar)))
       dT = 100d0
       if (T1 > Teff_jump + dT) then
          alfa = 1
       else if (T1 < Teff_jump - dT) then
          alfa = 0
       else
          alfa = 0.5d0*(T1 - (Teff_jump - dT)) / dT
       end if
    end if

    if (alfa > 0) then ! eval hot side wind (eqn 24)
       vinf_div_vesc = 2.6d0 ! this is the hot side galactic value
       vinf_div_vesc = vinf_div_vesc*pow(Z/Zsolar,0.13d0) ! corrected for Z
       logMdot = &
          - 6.697d0 &
          + 2.194d0*log10(L1/Lsun/1d5) &
          - 1.313d0*log10(M1/Msun/30d0) &
          - 1.226d0*log10(vinf_div_vesc/2d0) &
          + 0.933d0*log10(T1/4d4) &
          - 10.92d0*pow2(log10(T1/4d4)) &
          + 0.85d0*log10(Z/Zsolar)
       w1 = exp10(logMdot)
    else
       w1 = 0
    end if

    if (alfa < 1) then ! eval cool side wind (eqn 25)
       vinf_div_vesc = 1.3d0 ! this is the cool side galactic value
       vinf_div_vesc = vinf_div_vesc*pow(Z/Zsolar,0.13d0) ! corrected for Z
       logMdot = &
          - 6.688d0 &
          + 2.210d0*log10(L1/Lsun/1d5) &
          - 1.339d0*log10(M1/Msun/30d0) &
          - 1.601d0*log10(vinf_div_vesc/2d0) &
          + 1.07d0*log10(T1/2d4) &
          + 0.85d0*log10(Z/Zsolar)
       w2 = exp10(logMdot)
    else
       w2 = 0
    end if

    w = alfa*w1 + (1 - alfa)*w2

   ! if (dbg) write(*,*) 'vink wind', w

 end subroutine eval_Vink_wind

subroutine eval_Sabhahit_wind(w)
            real(dp), intent(inout) :: w
            real(dp) :: alfa, w1, w2, Teff_jump, logMdot, dT, vinf_div_vesc
            real(dp), parameter :: Zsolar = 0.019d0 ! for Vink et al formula
            ! alfa = 1 for hot side, = 0 for cool side
            if (T1 > 27500d0) then
               alfa = 1
            else if (T1 < 22500d0) then
               alfa = 0
            else ! use Vink et al 2001, eqns 14 and 15 to set "jump" temperature
               Teff_jump = 1d3*(61.2d0 + 2.59d0*(-13.636d0 + 0.889d0*log10(Z/Zsolar)))
               dT = 100d0
               if (T1 > Teff_jump + dT) then
                  alfa = 1
               else if (T1 < Teff_jump - dT) then
                  alfa = 0
               else
                  alfa = 0.5d0*(T1 - (Teff_jump - dT)) / dT
               end if
            end if
            if (alfa > 0) then ! eval hot side wind (eqn 10 in Sabhahit)
               vinf_div_vesc = 2.6d0 ! this is the hot side galactic value
               vinf_div_vesc = vinf_div_vesc*pow(Z/Zsolar,0.13d0) ! corrected for Z
               logMdot = &
                  - 6.697d0 &
                  + 2.194d0*log10(L1/Lsun/1d5) &
                  - 1.313d0*log10(M1/Msun/30d0) &
                  - 1.226d0*log10(vinf_div_vesc/2d0) &
                  + 0.933d0*log10(T1/4d4) &
                  - 10.92d0*pow2(log10(T1/4d4)) &
                  + 0.85d0*log10(Z/Zsolar)
               w1 = exp10(logMdot)
            else
               w1 = 0
            end if
            if (alfa < 1) then ! eval cool side wind (eqn 25)
               vinf_div_vesc = 1.3d0 ! this is the cool side galactic value
               vinf_div_vesc = vinf_div_vesc*pow(Z/Zsolar,0.13d0) ! corrected for Z
               logMdot = &
                  - 9.552 &
                  + 4.77d0*log10(1+X) &
                  + 4.77d0*log10(L1/Lsun/1d5) &
                  - 3.99d0*log10(M1/Msun/30d0) &
                  - 1.226d0*log10(vinf_div_vesc/2d0) &
                  + 0.5d0*log10(Z/Zsolar)
               w2 = exp10(logMdot)
            else
               w2 = 0
            end if
            w = alfa*w1 + (1 - alfa)*w2
      !      if (dbg) write(*,*) 'Sabhahit wind', w
         end subroutine eval_Sabhahit_wind

end subroutine my_wind

real(dp) function lambda_CE(s, adjusted_energy, star_core_mass_CE)
  type (star_info), pointer :: s
  integer :: k
  real(dp) :: E_bind, E_bind_shell, star_core_mass_CE, E_bind_grav_shell, E_bind_thermal_shell, E_bind_grav, E_bind_thermal
  real(dp) :: adjusted_energy(:)

  if (s% m(1) <= (star_core_mass_CE)) then
     lambda_CE = 1d99 ! no envelope, so immediately have a "succesfull envelope ejection"
  else
     E_bind = 0.0d0
     E_bind_shell = 0.0d0
     E_bind_grav = 0.0d0
     E_bind_grav_shell = 0.0d0
     E_bind_thermal = 0.0d0
     E_bind_thermal_shell = 0.0d0
     do k=1, s% nz
        if (s% m(k) > (star_core_mass_CE)) then !envelope is defined as everything above star_core_mass_CE.

           E_bind_grav_shell = (s% cgrav(1) * s% m(k) * s% dm_bar(k))/(s% r(k))
           E_bind_grav = E_bind_grav + E_bind_grav_shell

           E_bind_thermal_shell = s% dm(k) * adjusted_energy(k)
           E_bind_thermal = E_bind_thermal + E_bind_thermal_shell

        end if
     end do

    !  write(*,*) "E_bind_grav", E_bind_grav, "E_bind_thermal", E_bind_thermal

     E_bind = E_bind_thermal - E_bind_grav
     lambda_CE = - s% cgrav(1) * (s% m(1)) * ((s% m(1)) - star_core_mass_CE)/(E_bind * s% r(1))
  end if

end function lambda_CE


real(dp) function get_ion_info(s,id,k)
  use ionization_def, only: num_ion_vals
  use ionization_lib, only: eval_ionization
  integer, intent(in) :: id, k
  integer :: ierr
  real(dp) :: ionization_res(num_ion_vals)
  type (star_info), pointer :: s
  ierr = 0
  call eval_ionization( &
       1d0 - (s% X(k) + s% Y(k)), s% X(k), s% Rho(k), s% lnd(k)/ln10, &
       s% T(k), s% lnT(k)/ln10, ionization_res, ierr)
  if (ierr /= 0) ionization_res = 0
  get_ion_info = ionization_res(id)
end function get_ion_info




  !include 'POSYDON_single_stars.inc'

end module run_star_extras
