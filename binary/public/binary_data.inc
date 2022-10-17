! state saved and restored at each step (only for donor if doing binary)

integer :: model_number, model_number_old

real(dp) :: &
      ! age in yrs
      binary_age, binary_age_old, &
      mtransfer_rate, mtransfer_rate_old, &
      angular_momentum_j, angular_momentum_j_old, &
      ! separation in cms
      separation, separation_old, &
      eccentricity, eccentricity_old, &
      ! (r-rl)/rl for each star
      rl_relative_gap(2), rl_relative_gap_old(2), &
      ! photosphere radius for each star (as computed by the evolution code) (cm)
      r(2), r_old(2), &
      ! roche lobe radius for each star (as derived from masses and separation) (cm)
      rl(2), rl_old(2), &
      ! star masses
      m(2), m_old(2), &
      dt, dt_old, &
      env(2), env_old(2), &
      period, period_old, &
      max_timestep, max_timestep_old, &
      change_factor, change_factor_old, &
      min_binary_separation

logical :: &
      using_jdot_mb(2), using_jdot_mb_old(2), &
      CE_flag, CE_flag_old, & ! true during a CE phase
      CE_init, CE_init_old    ! true if variables for CE have been initiated

integer :: &
      CE_num1, CE_num1_old, & !count the number of times each star has undergone a CE
      CE_num2, CE_num2_old

real(dp) :: &
      CE_lambda1, CE_lambda1_old, & ! lambda value, updated through a CE ejection, zero before CE.
      CE_lambda2, CE_lambda2_old, &
      CE_Ebind1, CE_Ebind1_old, & !same as lambda for binding energy of ejected layers
      CE_Ebind2, CE_Ebind2_old

! during CE, this counts how much time the star has been detached
real(dp) :: &
      CE_years_detached, CE_years_detached_old

integer :: generations ! 1 means no old values, 2 means have old

! array indexes for the donor and the accretor
integer :: d_i = 1, d_i_old, a_i = 2, a_i_old, &
      ! index for star that is taken as point mass
      point_mass_i, point_mass_i_old

! flag indicating if mass transfer due to RLOF is ignored or not
logical :: ignore_rlof_flag, ignore_rlof_flag_old

! flag indicating if star 2 is modeled as twin of 1
logical :: model_twins_flag, model_twins_flag_old

integer :: dt_why_reason, dt_why_reason_old

! logical values indicating if we have star pointers for each star.
! There are saved in photos, but since they're set once at the beginning of the
! run, there's no need for old values.
! These are used in runs where point_mass_i switches from zero to another
! value, to save photos of both stars in a consistent way.
logical :: have_star_1, have_star_2

! *************** end of info that is saved and restored for restarts  ***************

! work variables -- not saved or restored. must be set at each step.

real(dp) :: time_step ! in yrs

real(dp) :: jdot, jdot_mb, jdot_gr, jdot_ml, jdot_ls, jdot_missing_wind, extra_jdot

real(dp) :: fixed_delta_mdot

!rate of mass loss from the system due to inefficient mass accretion by each star
real(dp) :: mdot_system_transfer(2)
!rate of mass loss from the system due to winds
real(dp) :: mdot_system_wind(2)
!rate of mass loss from the system from circumbinary coplanar toroid
real(dp) :: mdot_system_cct

!rate of mass transfer from one star to the other due to winds
real(dp) :: mdot_wind_transfer(2)
!mass transfer fraction from winds
real(dp) :: wind_xfer_fraction(2)

real(dp) :: fixed_xfer_fraction ! defined as 1-alpha-beta-delta, for the fixed alpha,beta,delta set in controls

!rate of mass change for each component, after accounting for wind mass loss,
!wind accretion and roche lobe overflow
real(dp) :: component_mdot(2)

! Eddington mass-accretion rate, in units of Msun/secyer
real(dp) :: mdot_edd
! Efficiency for Eddignton limited mass-accretion rate
real(dp) :: mdot_edd_eta
! Accretion luminosity from a point mass compact object
real(dp) :: accretion_luminosity
! 'Equivalent' BH initial mass that would have been required to evolve to the
! initial mass and spin parameter of the BH.
real(dp) :: eq_initial_bh_mass

! When the donor reaches critical rotation during the steps of the implicit MT calculation,
! things can start to go very slowly, as for each iteration the implicit wind calculation
! will do many steps. This flag checks if during the steps of calculation of the transfer
! rate this happens, to simplify the calculation of the implicit wind.
logical :: donor_started_implicit_wind

integer :: star_ids(2), star_extra_ids(2)
! pointers that are adjusted to donor and accretor
type (star_info), pointer :: s_donor, s_accretor
! pointers to fixed stars
type (star_info), pointer :: s1, s2

integer :: donor_id, accretor_id
character (len = strlen) :: last_photo_filename

logical :: evolve_both_stars

! Warning from when run_binary_extras changed to a hook
logical :: warn_binary_extra

! variable used to store mtransfer_rate used in a step, as it is rewritten
! by binary_check_model
real(dp) :: step_mtransfer_rate

! if true, then the following step will use a timestep reduced by dt_reduction_factor_for_j
logical :: have_to_reduce_timestep_due_to_j

! info for circularisation and eccentricity enhancement
real(dp) :: edot, edot_tidal, edot_enhance, extra_edot

! info for implicit rlo
real(dp) :: mdot_lo, mdot_hi
real(dp) :: implicit_function_lo, implicit_function_hi
logical :: have_mdot_lo = .false., have_mdot_hi = .false.
integer :: num_tries = 0

! info for Ritter and Kolb & Ritter prescriptions
real(dp) :: ritter_h, ritter_exponent, mdot_thin, mdot_thin0, mdot_thick

! type of accretion (none=0, ballistic=1, Keplerian=2)
integer :: accretion_mode = 0
! ratio of accreted specific angular momentum to that of a Keplerian orbit at R star
real(dp) :: acc_am_div_kep_am

! sync timescales for both stars
real(dp) :: t_sync_1, t_sync_2

! parameters to integrate over 1 orbit (anomaly and time coordinates)
real(dp), pointer :: theta_co(:), time_co(:), mdot_donor_theta(:)
! these are used to calculate edot_enhance
real(dp), pointer :: edot_theta(:), e1(:), e2(:), e3(:)

! If this is true, then the binary timestep routines will not apply
! hard limits. Set to false at the end of each step.
logical :: ignore_hard_limits_this_step

! Parameters for common envelope evolution
integer :: CE_nz ! number of points in model at CE onset
real(dp) :: CE_initial_radius ! radius (in cm) at the onset of CE
real(dp) :: CE_initial_separation ! separation (in cm) at the onset of CE
real(dp) :: CE_initial_Mdonor ! donor mass (in grams) at the onset of CE
real(dp) :: CE_initial_Maccretor ! accretor mass (in grams) at the onset of CE
real(dp) :: CE_initial_age ! age (in seconds) at the onset of CE
real(dp) :: CE_initial_model_number ! binary model_number at the onset of CE
! values for binary are stored as well, they need not match the ones of the star
! undergoing CE
real(dp) :: CE_b_initial_age ! binary age (in seconds) at the onset of CE
real(dp) :: CE_b_initial_model_number ! binary model_number at the onset of CE
! arrays used to interpolate properties from the star at the onset of CE
real(dp), pointer :: CE_m(:), CE_entropy(:)
real(dp), pointer :: CE_U_in(:), CE_U_out(:), CE_Omega_in(:), CE_Omega_out(:)
!! values computed at each step during CE
!real(dp) :: CE_Ebind ! binding energy of removed layers (in ergs)
!real(dp) :: CE_delta_Ecore ! change in binding energy of internal layers (in ergs)
!real(dp) :: CE_Uout ! thermal energy of removed layers (includes recombination energy)
!real(dp) :: CE_lambda ! lambda value if star is cut at this point: -standard_cgrav*Mi*Mf/Ri/(Ebind)

! more extras -- for run_binary_extras applications with state.
! if your extra state info is all held in the following variables,
! restarts and retries will work without doing anything special.
integer, dimension(binary_num_xtra_vals) :: ixtra, ixtra_old
real(dp), dimension(binary_num_xtra_vals) :: xtra, xtra_old
logical, dimension(binary_num_xtra_vals) :: lxtra, lxtra_old

! "other" procedures
procedure(other_rlo_mdot_interface), pointer, nopass :: &
      other_rlo_mdot => null()

procedure(other_check_implicit_rlo_interface), pointer, nopass :: &
      other_check_implicit_rlo => null()

procedure(other_implicit_function_to_solve_interface), pointer, nopass :: &
      other_implicit_function_to_solve => null()

procedure(other_tsync_interface), pointer, nopass :: &
      other_tsync => null()

procedure(other_sync_spin_to_orbit_interface), pointer, nopass :: &
      other_sync_spin_to_orbit => null()

procedure(other_mdot_edd_interface), pointer, nopass :: &
      other_mdot_edd => null()

procedure(other_adjust_mdots_interface), pointer, nopass :: &
      other_adjust_mdots => null()

procedure(other_accreted_material_j_interface), pointer, nopass :: &
      other_accreted_material_j => null()

procedure(other_jdot_interface), pointer, nopass :: &
      other_jdot_mb => null()

procedure(other_jdot_interface), pointer, nopass :: &
      other_jdot_gr => null()

procedure(other_jdot_interface), pointer, nopass :: &
      other_jdot_ml => null()

procedure(other_jdot_interface), pointer, nopass :: &
      other_extra_jdot => null()

procedure(other_jdot_interface), pointer, nopass :: &
      other_jdot_ls => null()

procedure(other_jdot_interface), pointer, nopass :: &
      other_jdot_missing_wind => null()

procedure(other_binary_wind_transfer_interface), pointer, nopass :: &
      other_binary_wind_transfer => null()

procedure(other_edot_interface), pointer, nopass :: &
      other_edot_tidal => null()

procedure(other_edot_interface), pointer, nopass :: &
      other_edot_enhance => null()

procedure(other_edot_interface), pointer, nopass :: &
      other_extra_edot => null()

procedure(other_e2_interface), pointer, nopass :: &
      other_e2 => null()

procedure(other_CE_init_interface), pointer, nopass :: &
      other_CE_init => null()

procedure(other_CE_rlo_mdot_interface), pointer, nopass :: &
      other_CE_rlo_mdot => null()

procedure(other_CE_binary_evolve_step_interface), pointer, nopass :: &
      other_CE_binary_evolve_step => null()

procedure(other_CE_binary_finish_step_interface), pointer, nopass :: &
      other_CE_binary_finish_step => null()

procedure(extras_binary_startup_interface), pointer, nopass :: &
      extras_binary_startup => null()

procedure(extras_binary_start_step_interface), pointer, nopass :: &
      extras_binary_start_step => null()

procedure(extras_binary_check_model_interface), pointer, nopass :: &
      extras_binary_check_model => null()

procedure(extras_binary_finish_step_interface), pointer, nopass :: &
      extras_binary_finish_step => null()

procedure(extras_binary_after_evolve_interface), pointer, nopass :: &
      extras_binary_after_evolve => null()

procedure(how_many_extra_binary_history_columns_interface), pointer, nopass :: &
      how_many_extra_binary_history_columns => null()

procedure(data_for_extra_binary_history_columns_interface), pointer, nopass :: &
      data_for_extra_binary_history_columns => null()

procedure(how_many_extra_binary_history_header_items_interface), pointer, nopass :: &
      how_many_extra_binary_history_header_items => null()

procedure(data_for_extra_binary_history_header_items_interface), pointer, nopass :: &
      data_for_extra_binary_history_header_items => null()

procedure (other_binary_photo_write_interface), pointer, nopass :: &
      other_binary_photo_write => null()
procedure (other_binary_photo_read_interface), pointer, nopass :: &
      other_binary_photo_read => null()


! info for binary_history.data
integer, pointer :: history_column_spec(:)
logical, pointer :: history_column_is_integer(:) ! true if integer value; false if double
real(dp), pointer :: history_column_values(:)
logical :: need_to_update_binary_history_now

logical :: doing_first_model_of_run = .false.
! open new history file is used to determine when the binary_history.data file
! needs to be open to append to it (when doing a restart), or if it needs to
! be rewritten (when starting a fresh run)
logical :: open_new_history_file = .false.

! info on terminal output
integer :: recent_binary_log_header
logical :: just_wrote_binary_terminal_header
