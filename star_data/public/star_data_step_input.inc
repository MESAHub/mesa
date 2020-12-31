! star_data_step_input.inc

   ! input for taking a step.
   ! modified during the step so that at end holds the new model.
       
      integer :: model_number
      
      real(dp) :: time ! elapsed time (seconds)
      real(dp) :: dt, dt_next ! timestep (seconds)
   
      ! total baryonic mass (grams)
      real(dp) :: mstar
      real(dp) :: xmstar ! baryonic mass external to M_center
         ! thus mstar = M_center + xmstar

      ! the number of zones in the model
      integer :: nz

      ! dq(k) is the fraction of xmstar between points k and k+1
      ! NOTE: near the surface, dq's can be extremely small
      ! the code should only use q(k)-q(k+1) in cases where 0 is an okay value.
      real(dp), pointer, dimension(:) :: dq ! (nz)
        
      ! vectors of structure variables: xh(i,k) holds hydro structure variable i for point k.
         ! note: abundances are in xa, not in xh.
         ! also note that mass and rotation are not included in xh
      real(dp), pointer, dimension(:,:) :: xh ! (nvar,nz)
   
      ! composition variables: xa(j,k) holds mass fraction of species j for location k
      ! j is an integer from 1 to s% species
      ! chem_id(j) gives the iso identifier from chem_def
      ! net_iso(chem_id(j)) = j
      real(dp), pointer, dimension(:,:) :: xa ! (species,nz)
      
      ! rotation
      real(dp), pointer, dimension(:) :: j_rot ! (nz)
        ! j_rot(k) is specific angular momentum at outer edge of cell k; = i_rot*omega
      real(dp), pointer, dimension(:) :: omega ! (nz)
      
      
      ! center boundary
      
      ! typically, mesa/star models an entire star and the center boundary has M, R, v, and L all = 0.
      ! however, in some situations it is useful to model only the envelope
      ! and use a nonzero center boundary.
      real(dp) :: M_center ! total baryonic mass (g) interior to innermost cell (usually 0)
      real(dp) :: v_center ! velocity (cm/s) of inner edge of innermost cell (usually 0)
      real(dp) :: R_center ! radius (cm) of inner edge of innermost cell (usually 0)
         ! R_center is automatically updated at each timestep if velocity_center is nonzero.
      real(dp) :: L_center ! luminosity (erg/s) at inner edge of innermost cell (usually 0)
      
      ! one reason for using a nonzero center boundary is to resolve small features in the envelope.
      ! the full model (with 0 center boundary) includes the envelope of course, but because of the
      ! limits of real(dp) we are limited to resolving envelope cell masses of about 10^-14*mstar.
      ! If that is adequate for your problem, then there is no reason to use a nonzero center boundary.
      ! But if you need to resolve smaller cell masses and you can ignore the center details,
      ! then a nonzero center can be the answer.  The code is written so that it can resolve cell masses
      ! down to about 10^-14*xmstar, where xmstar is the mass of the envelope exterior to M_center.
      ! So, for example, if xmstar = 10^-16*mstar,
      ! you can now resolve down to 10^-30*mstar in the envelope.
      ! The nonzero center boundaries are a big help with mass resolution,
      ! but they do not improve radial resolution -- that is still limited to about 10^-14*R.
      ! The remeshing algorithm knows that it must not split cells that are near that limit.
      ! In practice, this limit on radial resolution has not been a issue.
      ! For example, with R of 10 km for a neutron star, the ability to resolve 10^-8 cm is quite enough!
      ! However for mass resolution with a 1.6 Msun neutron star, we need to resolve small mass cells
      ! down to 10^10 g or less while the smallest amount of the total mass we can resolve is ~ 10^19 g.
      ! So for the NS envelope case, we need to use a nonzero center boundary.

      ! flags indicating extra variables and equations in addition to the minimal set
      logical :: v_flag ! true for velocity variables at cell boundaries     
      logical :: u_flag ! true for velocity variables at cell centers     
      logical :: rotation_flag ! true for rotation variables     
      logical :: RTI_flag ! true for Rayleigh-Taylor instabilities     
      logical :: conv_vel_flag ! true for conv_vel as a solver variable     
      logical :: w_div_wc_flag ! true for w_div_wc as a solver variable     
      logical :: j_rot_flag ! true for j_rot as a solver variable     
      logical :: RSP_flag ! true if we using RSP      
      logical :: Eturb_flag ! true for turbulent energy as a solver variable    

      ! index definitions for the variables (= 0 if variable not in use)
         integer :: i_lnd ! ln(cell density average by mass)
         integer :: i_lnT ! ln cell temperature average by mass
         integer :: i_lnR ! ln radius at outer face of cell
         integer :: i_lum ! luminosity at outer face of cell
         integer :: i_v ! Lagrangian velocity at outer face of cell
         integer :: i_u ! Lagrangian velocity at center of cell
         integer :: i_alpha_RTI ! alpha for Rayleigh-Taylor instabilities
         integer :: i_eturb_RSP, i_erad_RSP, i_Fr_RSP ! only when RSP_flag is true
         integer :: i_ln_cvpv0 ! ln(conv_vel + v0) at the outer face of cell.
                               ! v0 is an arbitrary user defined value.
         integer :: i_w_div_wc ! ratio of omega to critical value without the Eddington factor
         integer :: i_j_rot ! specific angular momentum associated to cell face
         integer :: i_eturb ! cell specific turbulent energy. only when Eturb_flag is true.

      ! index definitions for the equations (= 0 if equation not in use)
         integer :: i_dlnR_dt ! dlnR/dt = v/r
         integer :: i_dv_dt ! momentum conservation using v
         integer :: i_du_dt ! momentum conservation using u
         integer :: i_dlnd_dt ! mass conservation
         integer :: i_dlnE_dt ! energy conservation
         integer :: i_equL ! luminosity
         integer :: i_dalpha_RTI_dt ! alpha for Rayleigh-Taylor instabilities
         integer :: i_deturb_RSP_dt, i_derad_RSP_dt, i_dFr_RSP_dt ! only when RSP_flag is true.
         integer :: i_dln_cvpv0_dt ! equation for convection velocities
         integer :: i_equ_w_div_wc ! equation for w_div_wc
         integer :: i_dj_rot_dt ! equation for specific angular momentum
         integer :: i_deturb_dt ! turbulent energy equation. only when Eturb_flag is true.

      ! names of variables and equations
         character (len=name_len), dimension(:), pointer :: nameofvar, nameofequ ! (nvar)

         
      ! chem variables, if any, follow the structure variables
         integer :: i_chem1 ! index of the first chem variable (0 means no chem variables)
            ! if i_chem1 /= 0, then i_chem1 == nvar_hydro + 1
         integer :: i_dchem1 ! index of the first delta chem variable (0 means no delta chem variables)
         
      ! chem equations follow the structure equations
         integer :: equchem1 ! index of the first chem equation (0 means no chem equations)
            ! if equchem1 /= 0, then equchem1 == nvar_hydro + 1
         integer :: equdchem1 ! index of the first delta chem equation (0 means no delta chem equations)
        

      integer :: nvar_hydro ! the number of structure variables per zone      
      integer :: species ! species = number of isotopes in the current net
      integer :: nvar_chem ! the number of chemical abundance variables per zone
      ! = 0 for operator splitting with separate evolution of abundances
      ! = species to solve for all abundances along with structure
      integer :: nvar ! = nvar_hydro + nvar_chem

      ! net_name identifies which net we are currently using
      character (len=net_name_len) :: net_name
      ! num_reactions = number of reactions in the current net
      integer :: num_reactions
      
      ! chem info
      integer, pointer :: chem_id(:) ! maps species to chem id
         ! index from 1 to species in current net
         ! value is between 1 and num_chem_isos         
      integer, pointer :: net_iso(:) ! maps chem id to species number
         ! index from 1 to num_chem_isos
         ! value is 0 if the iso is not in the current net
         ! else is value between 1 and number of species in current net
         
        
      logical :: D_omega_flag ! D_omega_flag is true if we are using D_omega for mixing        
      logical :: am_nu_rot_flag! am_nu_rot_flag is true if we are using smoothed am_nu_rot variables

      real(dp) :: dt_limit_ratio ! for timestep low-pass limiter

           
         
      ! history
         integer, pointer :: profile_column_spec(:)
         integer :: number_of_history_columns
            ! number_of_history_columns = num_standard_history_columns + num extra history columns
         real(dp), pointer :: history_values(:) ! (number_of_history_columns)
         logical, pointer :: history_value_is_integer(:) ! (number_of_history_columns)
         character (len=maxlen_history_column_name), pointer :: &
            history_names(:) ! (number_of_history_columns)
         type (integer_dict), pointer :: history_names_dict
            ! maps names to index in history_values
         integer :: model_number_of_history_values
         logical :: need_to_set_history_names_etc
         ! the following does not include "extra" columns
         integer, pointer :: history_column_spec(:) ! (num_standard_history_columns)
      
      ! profile
         character (len=strlen) :: &
            model_profile_filename, model_controls_filename, model_data_filename
         character (len=strlen) :: &
            most_recent_profile_filename, most_recent_controls_filename, &
            most_recent_model_data_filename
         
      ! pgstar
         type (pgstar_hist_node), pointer :: pgstar_hist => null()
      
      ! random number generator (marsaglia and zaman). period is 2**43 with
      ! 900 million different sequences. the state of the generator (for restarts)
         integer :: rand_i97, rand_j97
         real(dp) :: rand_u(rand_u_len), rand_c, rand_cd, rand_cm
      
      real(dp) :: cumulative_energy_error, cumulative_extra_heating
      integer :: total_num_solver_iterations, total_num_solver_relax_iterations
      integer :: total_num_solver_calls_made, total_num_solver_relax_calls_made
      integer :: total_num_solver_calls_converged, total_num_solver_relax_calls_converged
      integer :: num_retries

      logical :: have_initial_energy_integrals
      real(dp) :: total_energy_initial          
         
      integer :: mesh_call_number
      integer :: solver_call_number
      integer :: diffusion_call_number

      integer :: recent_log_header ! model_number when last output header info to terminal
      
      integer :: phase_of_evolution
