! star_data_work.inc
   ! working storage for data derived from inputs while taking a step.
   ! can be undefined when start to take a step (e.g. by calling fill_arrays_with_NaNs).
   ! not saved in models or in photos.
      
      real(dp) :: mstar_dot ! (gm/second)
      
      real(dp), pointer, dimension(:) :: D_omega
      
      ! when using the implicit wind calculation, which keeps omega_div_omega_crit
      ! below a certain threshold, store the explicit value of mstar_dot
      real(dp) :: explicit_mstar_dot
      real(dp) :: adjust_J_q ! q down to which angular momentum is adjusted when using do_adjust_J

      ! q is defined in terms of dq:
      !  q(1) = 1.0,
      !  for nz > k > 1, q(k) = q(k-1) - dq(k-1),
      !  q(nz) = dq(nz)
      real(dp), pointer, dimension(:) :: q

      real(dp), pointer :: m(:) ! baryonic mass coord
         ! m(k) is enclosed baryonic mass at outer edge of cell k
         ! m(k) = s% M_center + s% q(k)*s% xmstar
      real(dp), pointer :: dm(:)
         ! dm(k) is baryonic mass of cell k
         ! dm(k) = s% dq(k)*s% xmstar
      
      real(dp), pointer :: dm_bar(:) 
         ! baryonic mass associated with cell edge
         ! for k=2,nz-1 dm_bar(k) = (dm(k-1)+dm(k))/2
         ! dm_bar(1) = dm(1)/2
         ! dm_bar(nz) = dm(nz-1)/2 + dm(nz)
      
      real(dp), pointer :: T(:) ! Temperature
      real(dp), pointer :: lnT(:) ! log(T)
      real(dp), pointer :: rho(:) ! density
      real(dp), pointer :: lnd(:) ! log(rho)
      real(dp), pointer :: lnPgas(:) ! log(Pgas)
      real(dp), pointer :: L(:) ! L(k) is luminosity at outer edge of cell k
      real(dp), pointer :: r(:) ! r(k) is radius at outer edge of cell k
      real(dp), pointer :: lnR(:) ! lnR(k) is log(r(k))
      real(dp), pointer :: w(:) ! cell turbulent velocity.  only used when RSP2_flag is true
      real(dp), pointer :: v(:) ! only used when v_flag is true
         ! v(k) is velocity at outer edge of cell k
      real(dp), pointer :: u(:) ! only used when u_flag is true
         ! u(k) is velocity at center of cell k
      real(dp), pointer :: lnE(:) ! log(specific internal energy)
      real(dp), pointer :: alpha_RTI(:) ! only used when RTI_flag is true
         ! alpha_RTI(k) is measure of Rayleigh-Taylor instability at center of cell k

      ! composition information -- cell average values
      
      real(dp), pointer :: X(:) ! mass fraction hydrogen
      real(dp), pointer :: Y(:) ! mass fraction helium
      real(dp), pointer :: Z(:) ! mass fraction metals
      real(dp), pointer :: abar(:) ! average atomic weight
      real(dp), pointer :: zbar(:) ! average charge per nucleus
      real(dp), pointer :: z2bar(:) ! average charge squared per nucleus
      real(dp), pointer :: z53bar(:) ! average charge^5/3 per nucleus
      real(dp), pointer :: ye(:)
         ! ye is mean number free electrons per nucleon, assuming complete ionization
         
      ! gravitational vs baryonic mass
      real(dp), pointer :: mass_correction(:), mass_correction_start(:) ! = dm_gravitational/dm_baryonic
         ! calculated by chem module
      real(dp), pointer :: m_grav(:), m_grav_start(:)
         ! enclosed gravitational mass at cell outer edge
         ! the variable called "m" is the enclosed baryonic mass
            ! = number of enclosed baryons * atomic mass unit
         ! NOTE: if use_mass_corrections = .false., 
         ! then m_grav = baryonic mass

      
      ! eos variables -- cell average values
      
      real(dp), pointer :: energy(:) ! specific internal energy (erg/g)
      real(dp), pointer :: egas(:) ! specific gas energy (erg/g)
      real(dp), pointer :: erad(:) ! specific radiation energy (erg/g)
      real(dp), pointer :: Pgas(:) ! gas pressure (Pgas = Pion + Pelec)
      real(dp), pointer :: Prad(:) ! radiation pressure
      real(dp), pointer :: Peos(:) ! eos pressure (Peos = Prad + Pgas)
      real(dp), pointer :: lnPeos(:) ! log(Peos)
      real(dp), pointer :: grada(:)   ! d ln T / d ln P | adiabatic
      real(dp), pointer :: dE_dRho(:) ! dE_dRho at constant T&abundances
      real(dp), pointer :: Cv(:)   ! dE/dT at constant density&abundances
      real(dp), pointer :: Cp(:)         ! heat capacity at constant Peos&abundances
         ! Cp := dh_dT at constant Peos, specific heat at constant pressure
         ! where h is enthalpy, h = E + Peos/Rho
      real(dp), pointer :: lnS(:)   ! log(specific entropy)
      real(dp), pointer :: gamma1(:)      ! adiabatic coefficient.
         ! (d_ln(p)/d_ln(rho)) at constant entropy.
      real(dp), pointer :: gamma3(:)      ! gamma3 := 1 + dlnT_dlnRho at constant S, 
      real(dp), pointer :: eta(:)        ! electron degeneracy parameter
      real(dp), pointer :: gam(:)        ! plasma interaction parameter
      real(dp), pointer :: mu(:)         
         ! mean molecular weight per particle (ions + free electrons)
      real(dp), pointer :: lnfree_e(:) ! log(number of free electrons per nucleon)
      real(dp), pointer :: phase(:) ! phase (0 for liquid, 1 for solid, in-between for blend)
      real(dp), pointer :: latent_ddlnT(:) ! latent heat of phase transition per unit lnT (erg/g) because we blur the transition
      real(dp), pointer :: latent_ddlnRho(:) ! latent heat of phase transition per unit lnRho (erg/g) because we blur the transition
      real(dp), pointer :: chiRho(:) ! dlnPeos_dlnRho at constant T
      real(dp), pointer :: chiT(:) ! dlnPeos_dlnT at constant Rho
      real(dp), pointer :: QQ(:) ! thermal expansion coefficient

      ! eos blend information
      real(dp), pointer :: eos_frac_OPAL_SCVH(:)
      real(dp), pointer :: eos_frac_HELM(:)
      real(dp), pointer :: eos_frac_Skye(:)
      real(dp), pointer :: eos_frac_PC(:)
      real(dp), pointer :: eos_frac_FreeEOS(:)
      real(dp), pointer :: eos_frac_CMS(:)
      real(dp), pointer :: eos_frac_ideal(:)

      real(dp), pointer :: d_eos_dlnd(:,:) ! (num_eos_basic_results,nz)
      real(dp), pointer :: d_eos_dlnT(:,:) ! (num_eos_basic_results,nz)
      real(dp), pointer :: d_eos_dxa(:,:,:) ! (num_eos_d_dxa_results,species,nz)

      ! eos partials for use in calculating equation partials for Jacobian matrix
      ! calculated using partials of interpolating polynomials of lnPgas, lnE, lnS where appropriate
      real(dp), pointer :: chiRho_for_partials(:) ! dlnP_dlnRho 
      real(dp), pointer :: chiT_for_partials(:) ! dlnP_dlnT
      real(dp), pointer :: dE_dRho_for_partials(:)
      real(dp), pointer :: Cv_for_partials(:) ! Cv
      real(dp), pointer :: dS_dRho_for_partials(:)
      real(dp), pointer :: dS_dT_for_partials(:)
      real(dp), pointer :: dlnE_dxa_for_partials(:,:) ! (species,nz)
      real(dp), pointer :: dlnPeos_dxa_for_partials(:,:) ! (species,nz)
      
      real(dp), pointer :: csound(:) ! at cell center -- set when call eos
      real(dp), pointer :: csound_face(:) ! interpolated 
      real(dp), pointer :: rho_face(:) ! interpolated rho at cell outer boundary      
      real(dp), pointer :: scale_height(:) ! cm
      real(dp), pointer :: v_div_csound(:) ! velocity/csound_face
      real(dp), pointer :: entropy(:) ! NOTE: kerg/baryon units
      real(dp), pointer :: grav(:) ! uses gravitational mass
      real(dp), pointer :: tau(:) ! optical depth at cell outer boundary
      real(dp), pointer :: dr_div_csound(:) ! sound crossing time for cell

      real(dp), pointer :: ergs_error(:)

      real(dp), pointer :: gradr_factor(:)

      real(dp), pointer :: rmid(:) ! radius at center by mass (cm)
      
      ! rotation variables only calculated if rotatation_flag is true
         ! angular velocity at cell boundary
      type(auto_diff_real_star_order1), pointer, dimension(:) :: i_rot(:), fp_rot(:), ft_rot(:)

      real(dp), pointer :: w_div_w_crit_roche(:) ! fraction of critical rotation at the equator,
                                                 ! without the contribution of the Eddington factor
      real(dp), pointer :: r_polar(:) ! radius in polar direction
      real(dp), pointer :: r_equatorial(:) ! radius in equatorial direction    
      real(dp) :: rotational_mdot_boost

      ! angular momentum fluxes
      type(auto_diff_real_star_order1), pointer, dimension(:) :: j_flux(:)
           
      real(dp), pointer :: am_nu_rot(:)
      real(dp), pointer :: am_nu_non_rot(:) ! diffusion driven by other sources, e.g. convection
      
      real(dp), pointer :: am_nu_omega(:) ! for diffusion of rotation angular velocity omega
      real(dp), pointer :: am_nu_j(:) ! for diffusion of specific angular momentum j
      
      ! sig = (4 pi r^2 rho)^2 nu / dm
      real(dp), pointer :: am_sig_omega(:)
      real(dp), pointer :: am_sig_j(:)
      
      real(dp), pointer :: domega_dlnR(:) ! defined at cell boundaries
      real(dp), pointer :: richardson_number(:) ! defined at cell boundaries
      
      ! diffusion coefs defined at cell boundaries (cm^2/s)

      real(dp), pointer :: D_mix_non_rotation(:) ! from non-rotation sources
      real(dp), pointer :: D_mix_rotation(:) ! from non-rotation sources
         ! convection, semiconvection, thermohaline, overshooting
         
      real(dp), pointer :: D_visc(:) ! kinematic shear viscosity
      real(dp), pointer :: D_DSI(:) ! dynamical shear instability 
      real(dp), pointer :: D_SH(:) ! Solberg-Hoiland instability
      real(dp), pointer :: D_SSI(:) ! secular shear instability 
      real(dp), pointer :: D_ES(:) ! Eddington-Sweet circulation 
      real(dp), pointer :: D_GSF(:) ! Goldreich-Schubert-Fricke instability
       
      real(dp), pointer :: D_ST(:) ! Spruit dynamo mixing diffusivity
      real(dp), pointer :: nu_ST(:) ! Spruit dynamo effective viscosity 
      real(dp), pointer :: D_ST_start(:), D_ST_start_old(:) ! For time smoothing
      real(dp), pointer :: nu_ST_start(:), nu_ST_start_old(:) 
      real(dp), pointer :: omega_shear(:) ! max(1d-30,min(1d30,abs(domega_dlnr(k)/omega(k))))

      real(dp), pointer :: dynamo_B_r(:) ! magnetic field (Gauss)
      real(dp), pointer :: dynamo_B_phi(:) ! magnetic field (Gauss)


      ! kap
      real(dp), pointer :: opacity(:) ! cell average
      real(dp), pointer :: d_opacity_dlnd(:)
         ! this is partial of opacity wrt natural log of density at constant T
      real(dp), pointer :: d_opacity_dlnT(:)
         ! this is partial of opacity wrt natural log of T at constant density
      real(dp), pointer :: kap_frac_lowT(:) ! fraction of opacity from lowT tables
      real(dp), pointer :: kap_frac_highT(:) ! fraction of opacity from highT tables
      real(dp), pointer :: kap_frac_Type2(:) ! fraction of opacity from Type2 tables
      real(dp), pointer :: kap_frac_Compton(:) ! fraction of opacity from Compton_Opacity
      real(dp), pointer :: kap_frac_op_mono(:) ! fraction of opacity from OP mono
      real, pointer, dimension(:) :: & ! working storage
         op_mono_umesh1, op_mono_semesh1, op_mono_ff1, op_mono_rs1
      integer :: op_mono_nptot, op_mono_ipe, op_mono_nrad, op_mono_n
      real(dp), pointer :: op_mono_factors(:) ! (species)
      procedure (), pointer, nopass :: set_op_mono_factors => null()
         ! called whenever the net changes
   
      ! net
      integer :: screening_mode_value ! values defined in screen_def
         ! screening_mode matches the screening_mode control parameter
      real(dp), pointer :: rate_factors(:) ! (num_reactions)
         ! when reaction rates are calculated, they are multiplied by the
         ! corresponding values in this array.  
      procedure (), pointer, nopass :: set_rate_factors => null()
         ! called whenever the net changes

      integer, pointer :: which_rates(:) ! (rates_reaction_id_max)  controls for rates
      procedure (), pointer, nopass :: set_which_rates => null()

      real(dp), pointer :: eps_nuc(:) ! sum of reaction_eps_nuc for all reactions in net
         ! thermal ergs per gram per second from nuclear reactions
         ! (including losses to neutrinos produced in the reactions)
      real(dp), pointer :: eps_nuc_categories(:,:) ! (num_categories, nz)
      real(dp), pointer :: d_epsnuc_dlnd(:) ! partial wrt density
      real(dp), pointer :: d_epsnuc_dlnT(:) ! partial wrt temperature
      real(dp), pointer :: d_epsnuc_dx(:,:)  ! (species,nz)
         ! d_ex_dx(j, k) is partial of eps_nuc(k) wrt species(j)

      real(dp), pointer :: eps_nuc_neu_total(:) ! erg/gm/sec, (positive) energy lost to neutrinos in nuclear reactions

      real(dp), pointer :: dxdt_nuc(:,:) ! (species,nz)
         ! rate of change of mass fractions from nuclear reactions
         ! dxdt_nuc(j,k) is time derivative from nuclear burning of species(j) at point(k)
      real(dp), pointer :: d_dxdt_nuc_dRho(:,:) ! (species,nz)
         ! d_dxdt_nuc_dRho(j, k) is partial of dxdt_nuc(k) wrt density
      real(dp), pointer :: d_dxdt_nuc_dT(:,:) ! (species,nz)
         ! d_dxdt_nuc_dRho(j, k) is partial of dxdt_nuc(k) wrt temperature
      real(dp), pointer :: d_dxdt_nuc_dx(:,:,:) ! (species,species,nz)
         ! d_dxdt_nuc_dx(i,j,k) is partial of dxdt_nuc for species i
         ! wrt abundance of species j at location k

      real(dp), pointer :: dxdt_mix(:,:) ! (species,nz)
         ! rate of change of mass fractions from mixing
         ! dxdt_mix(j,k) is time derivative from mixing of species(j) at point(k)
      real(dp), pointer :: d_dxdt_mix_dxm1(:) ! (nz)
      real(dp), pointer :: d_dxdt_mix_dx00(:) ! (nz)
      real(dp), pointer :: d_dxdt_mix_dxp1(:) ! (nz)
         ! partial of dxdt_mix(k) wrt to abundance in cell k and neighboring cells

      ! integrated eps_nuc_categories (ergs/sec)
      real(dp), pointer :: luminosity_by_category(:,:) ! (num_categories, nz)

      ! non-nuclear-reaction neutrino losses
      real(dp), pointer :: non_nuc_neu(:) ! positive
      real(dp), pointer :: d_nonnucneu_dlnd(:)
      real(dp), pointer :: d_nonnucneu_dlnT(:)
   
      real(dp), pointer :: nonnucneu_plas(:)
      real(dp), pointer :: nonnucneu_brem(:)
      real(dp), pointer :: nonnucneu_phot(:)
      real(dp), pointer :: nonnucneu_pair(:) 
      real(dp), pointer :: nonnucneu_reco(:) 

      ! irradiation heat profile
      real(dp), pointer :: irradiation_heat(:) ! (ergs/gm/sec)      
      
      ! eps_heat = eps_nuc - non_nuc_neu + extra_heat + irradiation_heat
      real(dp), pointer :: eps_heat(:) ! (ergs/gm/sec)
      
      ! extra factor for opacities (can be set by user)
      real(dp), pointer, dimension(:) :: extra_opacity_factor
      
      ! extra pressure profile (can be set by user)
      real(dp), pointer :: extra_pressure(:) ! (ergs/gm/sec)
      
      ! extra heat profile (can be set by user)
      type(auto_diff_real_star_order1), pointer, dimension(:) :: extra_heat
      
      ! extra gravity (can be set by user)  added to -G*m/r^2 in momentum equation
      type(auto_diff_real_star_order1), pointer, dimension(:) :: extra_grav
      
      ! extra torque (can be set by user)
      real(dp), pointer, dimension(:) :: &
         extra_jdot, extra_omegadot, &
         d_extra_jdot_domega_m1, d_extra_omegadot_domega_m1, &
         d_extra_jdot_domega_00, d_extra_omegadot_domega_00, &
         d_extra_jdot_domega_p1, d_extra_omegadot_domega_p1

      ! velocity time centering
      real(dp), pointer, dimension(:) :: vc
      real(dp) :: d_vc_dv
   
      ! gravitational constant (can be set by user and can be vary within model)
      real(dp), pointer :: cgrav(:)      
      
      ! local (per cell) factor for mesh_delta_coeff.
      real(dp), pointer :: mesh_delta_coeff_factor(:)      

      ! used by amr_split_merge to track which cells have undergone a remesh
      logical, pointer :: amr_split_merge_has_undergone_remesh(:)

      ! -T*Ds/Dt energy generation rate related to change in entropy
      type(auto_diff_real_star_order1), pointer, dimension(:) :: eps_grav_ad(:) ! '_ad' for auto diff
      ! partials
      real(dp), pointer :: d_eps_grav_dx(:,:)  ! (species,nz)
         ! d_eps_grav_dx(j, k) is partial of eps_grav(k) wrt species(j)
         ! these derivatives generally occur as part of the composition term
         
      real(dp), pointer :: eps_grav_composition_term(:) 
         ! only set when include_composition_in_eps_grav is true
         
      real(dp) :: mdot_acoustic_surface ! infalling material does work at surface
         ! this is ergs/sec added to the star by accretion work at the surface.
      real(dp) :: mdot_adiabatic_surface ! outgoing material can carry more energy
         ! with it than the surface entropy suggests when the thermal time is long.
         ! this is ergs/sec lost from the star.

      ! these store mesh information before remeshing
      integer :: prev_mesh_nz
      real(dp), pointer, dimension(:,:) :: prev_mesh_xh ! (nvar_hydro,prev_mesh_nz)
      real(dp), pointer, dimension(:,:) :: prev_mesh_xa ! (species,prev_mesh_nz)
      real(dp), pointer, dimension(:) :: & ! (prev_mesh_nz)
         prev_mesh_j_rot, prev_mesh_omega, prev_mesh_dq, prev_mesh_mlt_vc
      logical :: prev_mesh_species_or_nvar_hydro_changed ! specifies if either species or number of hydro variables
                                                         ! have been altered since prev_mesh info was stored

      real(dp) :: gradT_excess_alpha ! 1 for full efficiency boost; 0 for no boost
      real(dp) :: gradT_excess_min_beta
      real(dp) :: gradT_excess_max_lambda

      real(dp), pointer :: adjust_mlt_gradT_fraction(:)
         ! this vector defaults to -1.
         ! if f(k) k is >= 0 and <= 1, then
         ! gradT(k) from mlt is replaced by f(k)*grada_face(k) + (1-f(k))*gradr(k)
         
      real(dp), pointer :: gradT_excess_effect(:)
         ! this reports effect of MLT++
         ! it is the fraction of grada used to modify gradT, so 0 means no effect.
         
      real(dp), pointer :: superad_reduction_factor(:) ! only set when use_superad_reduction = .true.
         ! reports effect of superad_reduction.  value = 1d0 for no effect.

      ! mixing length alpha for MLT (can be set by user and can vary within model)
      real(dp), pointer :: alpha_mlt(:)

      ! NOTE: the following values are from the most recent call on mlt
         ! values for these at start of step have "_start" at end of name
         real(dp), pointer :: mlt_mixing_length(:) ! mixing_length_alpha*scale_height as computed by mlt
         integer, pointer :: mlt_mixing_type(:) ! as defined in mesa/const. before overshooting has been added.
         real(dp), pointer :: mlt_D(:) ! mlt diffusion coeff. (cm^2/sec)
         real(dp), pointer :: mlt_Gamma(:) ! local mlt convection efficiency.
         real(dp), pointer :: mlt_gradT(:) ! prior to adjustments
         real(dp), pointer :: mlt_cdc(:) 
         real(dp), pointer :: L_conv(:)
         
         integer, pointer :: tdc_num_iters(:)
         
      ! mlt results related to gradT are evaluated at each iteration and at end of step
         real(dp), pointer :: gradT_sub_grada(:) ! difference at face
         real(dp), pointer :: grada_face(:) ! grada evaluated at face
         real(dp), pointer :: gradT(:) ! desired dlnT/dlnP at face
         real(dp), pointer :: gradr(:) ! grad_rad at face
         type(auto_diff_real_star_order1), pointer, dimension(:) :: &
            gradT_ad, Y_face_ad, gradr_ad, mlt_vc_ad, grada_face_ad, &
            gradL_ad, scale_height_ad, Lambda_ad, mlt_D_ad, mlt_Gamma_ad
         
         
         logical, pointer :: fixed_gradr_for_rest_of_solver_iters(:)
         
      ! Convective turbulent velocity w -- when RSP2_flag is true
         
         type(auto_diff_real_star_order1), pointer, dimension(:) :: &
            PII_ad, Chi_ad, Eq_ad, COUPL_ad, Lr_ad, Lc_ad, Lt_ad
         real(dp), pointer, dimension(:) :: &
            Ptrb, Pvsc, Pvsc_start, Hp_face, Y_face, Y_face_start, PII, Chi, Eq, Uq, &
            SOURCE, DAMP, DAMPR, COUPL, Lr, Lc, Lt, Lt_start
         integer :: RSP2_num_periods
         real(dp) :: RSP2_period, &
            RSP2_DeltaR, RSP2_DeltaMag, RSP2_GRPDV, RSP2_GREKM, RSP2_GREKM_avg_abs
         
      ! RSP: Radial Stellar Pulsations. RSP derived from Radec Smolec's Program

         real(dp), pointer, dimension(:) :: &
            Vol, Vol_start, RSP_w, RSP_w_start, RSP_Et, COUPL_start, &
            Fr, Fr_start, Lc_start, Ptrb_start, Chi_start, f_Edd
         integer :: RSP_num_periods
         real(dp) :: RSP_dt, RSP_period, g_Edd, &
            rsp_DeltaR, rsp_DeltaMag, rsp_GRPDV, rsp_GREKM, rsp_GREKM_avg_abs
         real(dp), dimension(max_num_rsp_LINA_modes) :: &
            rsp_LINA_periods, rsp_LINA_growth_rates
      

      ! adjusted (non-local) mixing info -- evaluated at start of step.
         
         real(dp), pointer :: cdc(:) ! convective diffusion coeff. (gm^2/sec)
         real(dp), pointer :: D_mix(:) ! diffusion coeff. (cm^2/sec)
         integer, pointer :: mixing_type(:) ! mixing types are defined in mesa/mlt  
         real(dp), pointer :: cz_bdy_dq(:) ! sub-cell location of convection boundary
            ! gives dq distance inward from outer edge of cell(k)
            ! where the convection zone boundary happens.
            ! i.e., q location of boundary is q(k) - cz_bdy_dq(k)
            ! this is where gradr == gradL by linear interpolation in the cell.
      
         integer :: num_conv_boundaries ! boundaries of regions with mixing_type = convective_mixing
         ! boundaries are ordered from center to surface
         real(dp), pointer :: conv_bdy_q(:) ! (num_conv_boundaries)
            ! subcell location of boundary
         logical, pointer :: top_conv_bdy(:) ! (num_conv_boundaries)
         integer, pointer :: conv_bdy_loc(:) ! (num_conv_boundaries)
            ! if top_conv_bdy, top of region is between loc and loc+1
            ! else bottom of region is between loc and loc-1
         logical, pointer :: burn_h_conv_region(:) ! (num_conv_boundaries)
         logical, pointer :: burn_he_conv_region(:) ! (num_conv_boundaries)
         logical, pointer :: burn_z_conv_region(:) ! (num_conv_boundaries)
      
         integer :: num_mix_boundaries 
         ! boundaries of regions where D_mix > mix_region_min_D_mix
         ! mix_region_min_D_mix is a parameter in &controls.
         ! boundaries given in order from center to surface
         ! no bottom boundary at loc=nz included if center is mixed
         ! however, do include top boundary at loc=1 if surface is mixed
         integer :: num_mix_regions ! number of mix boundaries with top_mix_bdy = true
         real(dp), pointer :: mix_bdy_q(:) ! (num_mix_boundaries)
            ! sub-cell location where estimate D_mix = mix_region_min_D_mix
         integer, pointer :: mix_bdy_loc(:) ! (num_mix_boundaries)
         logical, pointer :: top_mix_bdy(:) ! (num_conv_boundaries)
         logical, pointer :: burn_h_mix_region(:) ! (num_conv_boundaries)
         logical, pointer :: burn_he_mix_region(:) ! (num_conv_boundaries)
         logical, pointer :: burn_z_mix_region(:) ! (num_conv_boundaries)

         real(dp), pointer :: sig(:) ! sig(k) is mixing flow across face k in (gm sec^1)
            ! sig(k) is D_mix*(4*pi*r(k)**2*rho_face)**2/dmavg
         real(dp), pointer :: sig_raw(:) ! before applying limits

         
      real(dp), pointer :: grad_density(:) ! smoothed dlnRho/dlnP at cell boundary
      real(dp), pointer :: grad_temperature(:) ! smoothed dlnT/dlnP at cell boundary
      real(dp), pointer :: gradL(:) ! Ledoux gradient at cell boundary
         ! gradL = grada + brunt_B
      real(dp), pointer :: gradL_composition_term(:)
      integer, pointer :: dominant_iso_for_thermohaline(:)
         ! chem id for dominant isotope in setting chemical gradient term.
         ! i.e. iso with largest local dX*(1+Z)/A
      
   
      ! element diffusion
      real(dp), pointer :: diffusion_D_self(:,:) ! self coeff for chem diffusion. [cm^2/sec]
      real(dp), pointer :: extra_diffusion_factor(:,:) ! diffusion scaling factor (species,nz)
      real(dp), pointer :: edv(:,:) ! velocity [cm/s] (species,nz)
      real(dp), pointer :: v_rad(:,:) ! velocity from radiative levitation [cm/s] (species,nz)
      real(dp), pointer :: g_rad(:,:) ! acceleration from radiative levitation [cm/s^2] (species,nz)
      real(dp), pointer :: typical_charge(:,:) ! (species,nz)
      real(dp), pointer :: diffusion_dX(:,:) ! (species,nz)
      real(dp), pointer :: E_field(:) ! (nz)
      real(dp), pointer :: eps_WD_sedimentation(:) ! (nz)
      real(dp), pointer :: eps_diffusion(:) ! (nz)
      real(dp), pointer :: g_field_element_diffusion(:) ! (nz)
      
      
      ! for Riemann solver      
      type(auto_diff_real_star_order1), pointer, dimension(:) :: &
         u_face_ad, P_face_ad ! 'ad' for auto_diff
      real(dp), pointer, dimension(:) :: &
         u_face_val, &
         d_uface_domega, &
         RTI_du_diffusion_kick, &
         u_face_start
      real(dp), pointer, dimension(:) :: P_face_start
      real(dp), pointer, dimension(:) :: &
         abs_du_div_cs, &
         abs_du_plus_cs


      ! power from nuclear reactions integrated from center out (ergs/second)
      real(dp), pointer :: L_nuc_burn(:) ! excluding photodisintegrations
      real(dp), pointer :: L_nuc_by_category(:,:) 

      
      ! Brunt-Vaisala information
      real(dp), pointer :: brunt_N2(:) ! at outer boundary of cell
      real(dp), pointer :: brunt_N2_composition_term(:)
      real(dp), pointer :: brunt_B(:) ! (at points)
         ! this is the Brassard et al B := -(chiY/chiT)*(dlnY/dlnP)
      real(dp), pointer :: unsmoothed_brunt_B(:) ! pre-smoothing
      real(dp), pointer :: smoothed_brunt_B(:) ! with smoothing applied for gradL
         
         
      ! asteroseismology info
      real(dp) :: delta_Pg ! g-mode period spacing for l=1 (seconds)
         ! sqrt(2)*pi^2/(integral of brunt_N/r dr)
      real(dp) :: nu_max ! estimate from scaling relation (microHz)
         ! nu_max = nu_max_sun * M/Msun / ((R/Rsun)^2 (Teff/Teff_sun)^0.5)
         ! with nu_max_sun = 3100 microHz, Teff_sun = 5777
      real(dp) :: acoustic_cutoff ! 1e6*0.5*g*sqrt(gamma1*rho/P) at surface (micro Hz)
      ! asteroseismology
      real(dp) :: delta_nu ! in micro Hz.  = 10^6/(2*photosphere_acoustic_r)
      
      
      ! atmosphere structure information (only created on demand)
      integer :: atm_structure_num_pts
      real(dp), pointer :: atm_structure(:,:) ! (num_results_for_create_atm,atm_structure_num_pts)
         ! defined at points in atmosphere
         ! atm_structure(:,1) is base of atmosphere
         ! atm_structure(:,atm_structure_num_pts) is top of atmosphere
        
      ! changes in values for cells during solver iterations 
      
         real(dp), pointer :: dxh_lnR(:) ! lnR(k) = dxh_lnR(k) + lnR_start(k); dxh_lnR from solver.
         real(dp), pointer :: dxh_lnd(:) ! lnd(k) = dxh_lnd(k) + lnd_start(k); dxh_lnd from solver.
         real(dp), pointer :: dxh_lnT(:) ! lnT(k) = dxh_lnT(k) + lnT_start(k); dxh_lnT from solver.
         real(dp), pointer :: dxh_v(:) ! v(k) = dxh_v(k) + v_start(k); dxh_v from solver.
         real(dp), pointer :: dxh_u(:) ! u(k) = dxh_u(k) + u_start(k); dxh_u from solver.
         real(dp), pointer :: dxh_alpha_RTI(:) ! alpha_RTI(k) = dxh_alpha_RTI(k) + alpha_RTI_start(k); dxh_alpha_RTI from solver.
        
      ! Rayleigh Taylor Instabilities  --  only used when RTI_flag is true

         real(dp), pointer :: dPdr_dRhodr_info(:) ! only used when RTI_flag is true
         real(dp), pointer :: dPdr_info(:) ! only used when RTI_flag is true
         real(dp), pointer :: dRhodr_info(:) ! only used when RTI_flag is true

         real(dp), pointer :: source_plus_alpha_RTI(:) ! only used when RTI_flag is true
         real(dp), pointer :: source_minus_alpha_RTI(:) ! only used when RTI_flag is true

         real(dp), pointer :: dudt_RTI(:) ! for cell
         real(dp), pointer :: dedt_RTI(:) ! for cell
      
         real(dp), pointer :: eta_RTI(:) 
            ! eta_RTI(k) is diffusion coefficient at cell boundary
         real(dp), pointer :: etamid_RTI(:)
         real(dp), pointer :: boost_for_eta_RTI(:)
         
         real(dp), pointer :: sig_RTI(:) ! sig_RTI(k) is mixing flow across face k in (gm sec^1)
            ! sig_RTI(k) is eta_RTI*(4*pi*r(k)**2*rho_face)**2/dmavg
         real(dp), pointer :: sigmid_RTI(:) ! sigmid_RTI(k) is mixing flow at center k in (gm sec^1)
            ! sigmid_RTI(k) is etamid_RTI*(4*pi*rmid(k)**2*rho)**2/dm

         ! the number of models that are currently stored
            integer :: generations ! either 1 or 2
               ! 1 means have current model only
               ! 2 means have current and old

         real(dp), pointer, dimension(:,:) :: xh_old, xa_old
         real(dp), pointer, dimension(:) :: &
            dq_old, q_old, j_rot_old, omega_old, mlt_vc_old
         real(dp) :: time_old, dt_old, mstar_old, xmstar_old, mstar_dot_old, &
            M_center_old, v_center_old, R_center_old, L_center_old
         integer :: nz_old, model_number_old
   
      
   ! other info
      ! NOTE: much of this is computed at the end of the step.
      ! so it is available for "extras_check_model" for example,
      ! but not for hook routines called during the step.
      ! if you aren't sure, check the sources by grep'ing in star/private.
         
      
      ! mixing regions (from surface inward)
      ! a "mixing region" is a run of points with the same value of mixing_type
      ! i.e., a region of convective_mixing, overshoot_mixing,
      ! semiconvective_mixing, or thermohaline_mixing
      
      integer :: num_mixing_regions ! the number of different regions in the current model
         ! from 0 to max_num_mixing_regions
      integer :: mixing_region_type(max_num_mixing_regions) ! values from mlt_def
      integer :: mixing_region_top(max_num_mixing_regions) ! index of point at top of region
      integer :: mixing_region_bottom(max_num_mixing_regions) ! index of point at bottom of region
      integer :: largest_conv_mixing_region ! index into the mixing_region arrays
      
      ! convection regions
      integer :: n_conv_regions
      real(dp), dimension(max_num_mixing_regions) :: cz_bot_mass, cz_top_mass
      
      ! top & bottom q of 2 largest mixing regions (mx1 the largest; mx2 the 2nd largest)
      ! conv_mx for convective mixing region; mx for any kind of mixing
      real(dp) :: & ! in Mstar units --- i.e., q values, m/Mstar
         conv_mx1_top, conv_mx1_bot, conv_mx2_top, conv_mx2_bot, &
         mx1_top, mx1_bot, mx2_top, mx2_bot
      ! top & bottom r (in Rsun units) of 2 largest mixing regions
      real(dp) :: &
         conv_mx1_top_r, conv_mx1_bot_r, conv_mx2_top_r, conv_mx2_bot_r, &
         mx1_top_r, mx1_bot_r, mx2_top_r, mx2_bot_r
      
      ! burn zones
         ! up to 3 zones where eps_nuc > burn_min1 erg/g/s
         ! for each zone have 4 numbers: start1, start2, end2, end1
         ! start1 is mass of inner edge where first goes > burn_min1 (or -20 if none such)
         ! start2 is mass of inner edge where first zone reaches burn_min2 erg/g/sec (or -20 if none such)
         ! end2 is mass of outer edge where first zone drops back below burn_min2 erg/g/s
         ! end1 is mass of outer edge where first zone ends (i.e. eps_nuc < burn_min1)
         ! similar for second and third zones
      real(dp) :: burn_zone_mass(4,3)
         
      ! power emitted from nuclear reactions after subtracting neutrinos (Lsun units)         
      real(dp) :: center_eps_burn(num_categories)
      real(dp) :: L_nuc_burn_total, L_by_category(num_categories)
      real(dp) :: power_nuc_neutrinos
      real(dp) :: power_nonnuc_neutrinos
      real(dp) :: power_neutrinos
   
      real(dp) :: power_nuc_burn_old
      real(dp) :: power_h_burn_old
      real(dp) :: power_he_burn_old
      real(dp) :: power_z_burn_old
      real(dp) :: power_photo_old
      
      
      ! surface       
      real(dp) :: L_surf ! (Lsun) -- depends on control tau_for_L_surf
         ! if tau_for_L_surf > 0, then is optical depth to set L_surf; else use L_surf = L(1)
      real(dp) :: v_surf ! (cm/second)
      real(dp) :: Lrad_div_Ledd_avg_surf
      real(dp) :: w_div_w_crit_avg_surf
      real(dp) :: omega_avg_surf
      real(dp) :: omega_crit_avg_surf
      real(dp) :: j_rot_avg_surf
      real(dp) :: v_rot_avg_surf ! cm/s
      real(dp) :: v_crit_avg_surf
      real(dp) :: v_div_v_crit_avg_surf      
      real(dp) :: opacity_avg_surf
      real(dp) :: logT_avg_surf
      real(dp) :: logRho_avg_surf
      real(dp) :: T_surf, P_surf ! results from atm
      real(dp) :: log_surface_luminosity ! log10(L_surf)
      real(dp) :: log_surface_radius ! log10(total radius in Rsun units)
      
      ! fractional abundances by mass near surface
         ! averages over a fraction = 'surface_avg_abundance_dq' of star mass
      real(dp) :: surface_h1
      real(dp) :: surface_he3
      real(dp) :: surface_he4
      real(dp) :: surface_c12
      real(dp) :: surface_n14
      real(dp) :: surface_o16
      real(dp) :: surface_ne20
   
      ! photosphere
      real(dp) :: L_phot ! photosphere luminosity (Lsun units)
      real(dp) :: photosphere_L ! Lsun units (= L_phot)
      real(dp) :: photosphere_r ! Rsun units
      real(dp) :: photosphere_m ! Msun units
      real(dp) :: photosphere_black_body_T ! black body with R and L at photosphere
      real(dp) :: photosphere_v ! cm/s
      real(dp) :: photosphere_csound ! cm/s
      real(dp) :: photosphere_opacity
      real(dp) :: photosphere_column_density ! g/cm^2
      real(dp) :: photosphere_logg
      real(dp) :: photosphere_acoustic_r ! sound travel time from center to photosphere (seconds)
      real(dp) :: photosphere_T
      integer :: photosphere_cell_k

      ! center region -- averages over center_avg_value_dq
      real(dp) :: log_max_temperature ! log10(max temperature in kelvin)
      real(dp) :: log_center_temperature ! log10(center temperature in kelvin)
      real(dp) :: log_center_density ! log10(center density in g/cm^3)
      real(dp) :: log_center_pressure ! log10(center pressure in dynes/cm^2)
      real(dp) :: center_degeneracy 
         ! the electron chemical potential in units of k*T
      real(dp) :: center_abar 
      real(dp) :: center_zbar 
      real(dp) :: center_mu 
      real(dp) :: center_ye 
      real(dp) :: center_eps_nuc, d_center_eps_nuc_dlnT, d_center_eps_nuc_dlnd
      
      real(dp) :: center_non_nuc_neu 
      real(dp) :: center_entropy ! in kerg per baryon 
      real(dp) :: fe_core_infall, non_fe_core_infall, non_fe_core_rebound
      real(dp) :: max_infall_speed_mass ! Msun
      real(dp) :: center_omega
      real(dp) :: center_omega_div_omega_crit 
         
      real(dp) :: center_h1  ! fractional abundance by mass
      real(dp) :: center_he3
      real(dp) :: center_he4
      real(dp) :: center_c12
      real(dp) :: center_n14
      real(dp) :: center_o16
      real(dp) :: center_ne20
      real(dp) :: center_si28

      real(dp) :: max_entropy ! in kerg per baryon 

      ! total baryonic mass
      real(dp) :: star_mass ! total stellar baryonic mass (Msun)
      real(dp) :: star_mdot ! d(star_mass)/dt (Msun/year)
      real(dp) :: star_mass_h1
      real(dp) :: star_mass_he3
      real(dp) :: star_mass_he4
      real(dp) :: star_mass_c12
      real(dp) :: star_mass_n14
      real(dp) :: star_mass_o16
      real(dp) :: star_mass_ne20
            
      real(dp) :: mass_conv_core ! (Msun) 0 if core is non-convective
      real(dp) :: mass_semiconv_core ! (Msun)
         ! = mass_conv_core if there is no semiconvective region above convective core

      real(dp) :: mass_bdy_core_overshooting ! (Msun) 0 if no core overshooting
      real(dp) :: radius_bdy_core_overshooting ! (Rsun) 0 if no core overshooting

      ! abundance boundaries
      
      real(dp) :: he_core_mass ! baryonic (Msun)
      real(dp) :: he_core_radius ! Rsun
      real(dp) :: he_core_lgT
      real(dp) :: he_core_lgRho
      real(dp) :: he_core_L ! Lsun
      real(dp) :: he_core_v
      real(dp) :: he_core_omega ! (s^-1)
      real(dp) :: he_core_omega_div_omega_crit
      integer :: he_core_k ! boundary is in this cell
      
      real(dp) :: co_core_mass ! baryonic (Msun)
      real(dp) :: co_core_radius ! Rsun
      real(dp) :: co_core_lgT
      real(dp) :: co_core_lgRho
      real(dp) :: co_core_L ! Lsun
      real(dp) :: co_core_v
      real(dp) :: co_core_omega ! (s^-1)
      real(dp) :: co_core_omega_div_omega_crit
      integer :: co_core_k ! boundary is in this cell
      
      real(dp) :: one_core_mass ! baryonic (Msun)
      real(dp) :: one_core_radius ! Rsun
      real(dp) :: one_core_lgT
      real(dp) :: one_core_lgRho
      real(dp) :: one_core_L ! Lsun
      real(dp) :: one_core_v
      real(dp) :: one_core_omega ! (s^-1)
      real(dp) :: one_core_omega_div_omega_crit
      integer :: one_core_k ! boundary is in this cell

      real(dp) :: fe_core_mass ! baryonic (Msun)
      real(dp) :: fe_core_radius ! Rsun
      real(dp) :: fe_core_lgT
      real(dp) :: fe_core_lgRho
      real(dp) :: fe_core_L ! Lsun
      real(dp) :: fe_core_v
      real(dp) :: fe_core_omega ! (s^-1)
      real(dp) :: fe_core_omega_div_omega_crit     
      integer :: fe_core_k ! boundary is in this cell

      real(dp) :: neutron_rich_core_mass ! baryonic (Msun)
      real(dp) :: neutron_rich_core_radius ! Rsun
      real(dp) :: neutron_rich_core_lgT
      real(dp) :: neutron_rich_core_lgRho
      real(dp) :: neutron_rich_core_L ! Lsun
      real(dp) :: neutron_rich_core_v
      real(dp) :: neutron_rich_core_omega ! (s^-1)
      real(dp) :: neutron_rich_core_omega_div_omega_crit     
      integer :: neutron_rich_core_k ! boundary is in this cell
      
      ! info about outermost outward moving shock
      ! excluding locations with q > max_q_for_outer_mach1_location
      real(dp) :: shock_mass ! baryonic (Msun)
      real(dp) :: shock_q
      real(dp) :: shock_radius ! (Rsun)
      real(dp) :: shock_velocity
      real(dp) :: shock_csound
      real(dp) :: shock_lgT
      real(dp) :: shock_lgRho
      real(dp) :: shock_lgP
      real(dp) :: shock_gamma1
      real(dp) :: shock_entropy
      real(dp) :: shock_tau
      real(dp) :: shock_pre_lgRho
      integer :: shock_k

      real(dp) :: min_gamma1

      ! timescales
      real(dp) :: dynamic_timescale 
         ! dynamic timescale (seconds) -- estimated by 2*pi*sqrt(r^3/(G*m))
      real(dp) :: kh_timescale 
         ! kelvin-helmholtz timescale (years) -- proportional to gravitational energy divided by luminosity
      real(dp) :: nuc_timescale 
         ! nuclear timescale (years) -- proportional to mass divided by luminosity

      ! times in other units (just for convenience)
      real(dp) :: star_age, time_days, time_years ! time is in seconds
      real(dp) :: time_step, dt_days, dt_years ! dt is in seconds
      
      real(dp) :: dt_start ! dt during first try for step
      
      ! a few things from start of step that will need later
      real(dp) :: prev_Ledd, starting_T_center, &
         surf_r_equatorial, surf_csound, surf_rho
   
      ! for solvers
      
         integer, pointer :: solver_iwork(:)
         real(dp), pointer :: solver_work(:)
         
         real(dp), pointer :: xh_start(:,:) ! (nvar_hydro,nz)
         real(dp), pointer :: xa_start(:,:) ! (species,nz)         
         real(dp), pointer :: xa_sub_xa_start(:,:) ! (species,nz)          
         real(dp), pointer :: lnd_start(:) ! (nz)
         real(dp), pointer :: Peos_start(:) ! (nz)
         real(dp), pointer :: lnPeos_start(:) ! (nz)
         real(dp), pointer :: lnPgas_start(:) ! (nz)
         real(dp), pointer :: lnT_start(:) ! (nz)
         real(dp), pointer :: energy_start(:) ! (nz)
         real(dp), pointer :: erad_start(:) ! (nz)
         real(dp), pointer :: egas_start(:) ! (nz)
         real(dp), pointer :: Prad_start(:) ! (nz)
         real(dp), pointer :: Pgas_start(:) ! (nz)
         real(dp), pointer :: lnS_start(:) ! (nz)
         real(dp), pointer :: lnR_start(:) ! (nz)
         real(dp), pointer :: w_start(:) ! (nz)
         real(dp), pointer :: Hp_face_start(:) ! (nz)
         real(dp), pointer :: v_start(:) ! (nz)
         real(dp), pointer :: u_start(:) ! (nz)
         real(dp), pointer :: L_start(:) ! (nz)
         real(dp), pointer :: r_start(:) ! (nz)
         real(dp), pointer :: rmid_start(:) ! (nz)
         real(dp), pointer :: rho_start(:) ! (nz)
         real(dp), pointer :: omega_start(:) ! (nz)
         real(dp), pointer :: ye_start(:) ! (nz)
         real(dp), pointer :: j_rot_start(:) ! (nz)
         real(dp), pointer :: i_rot_start(:) ! (nz)
         real(dp), pointer :: dxdt_nuc_start(:,:)
         real(dp), pointer :: luminosity_by_category_start(:,:) ! (num_categories,nz)
         real(dp), pointer :: grada_start(:) ! (nz)
         real(dp), pointer :: chiT_start(:) ! (nz)
         real(dp), pointer :: chiRho_start(:) ! (nz)
         real(dp), pointer :: Cp_start(:) ! (nz)
         real(dp), pointer :: Cv_start(:) ! (nz)
         real(dp), pointer :: dE_dRho_start(:) ! (nz)
         real(dp), pointer :: gam_start(:) ! (nz)
         real(dp), pointer :: T_start(:) ! (nz)
         real(dp), pointer :: mu_start(:) ! (nz)
         real(dp), pointer :: phase_start(:) ! (nz)
         real(dp), pointer :: latent_ddlnT_start(:) ! (nz)
         real(dp), pointer :: latent_ddlnRho_start(:) ! (nz)
         real(dp), pointer :: zbar_start(:) ! (nz)
         real(dp), pointer :: eps_nuc_start(:) ! (nz)
         real(dp), pointer :: non_nuc_neu_start(:) ! (nz)
         real(dp), pointer :: alpha_RTI_start(:) ! (nz)
         real(dp), pointer :: csound_start(:) ! (nz)
         real(dp), pointer :: opacity_start(:) ! (nz)

         real(dp) :: shock_mass_start ! Msun units like shock_mass
         
         ! for Lagrangian time derivatives where k < k_below_just_added
         real(dp) :: surf_lnS
         integer :: num_surf_revisions
            
         logical :: have_burner_storage
         integer :: burner_num_threads
         real(dp), pointer :: burner_storage(:)
            
         real(dp), pointer, dimension(:) :: & ! (nz)
            max_burn_correction, avg_burn_correction

         integer :: k_below_const_q
         integer :: k_const_mass
         real(dp) :: adjust_mass_outer_frac_sub1 ! for k 1:k_below_const_q-1
         real(dp) :: adjust_mass_mid_frac_sub1 ! for k k_below_const_q:k_const_mass-1
         real(dp) :: adjust_mass_inner_frac_sub1 ! for k k_const_mass:nz

         integer :: k_for_test_CpT_absMdot_div_L 
         integer :: k_below_just_added 
            ! "just added" == added during the current step
            ! i.e., largest k s.t. for all k' < k, cell k' has material just added this step.
            ! equivalently, this is the smallest k st. for all k' > k,
            !      material below k' was already in the star at the start of this step.
            ! if no mass was added, then this == 1.
               
      ! for remesh
         type (star_info), pointer :: other_star_info 
            ! holds pointers to "other" allocated arrays
            ! these are used during remesh
      
      ! misc other stuff        

         logical :: doing_finish_load_model, doing_relax, doing_struct_burn_mix, &
            RSP_just_set_velocities, &
            have_new_cz_bdy_info, using_gold_tolerances, need_to_setvars, have_new_generation, &
            using_velocity_time_centering, okay_to_set_mixing_info, &
            okay_to_set_mlt_vc, have_mlt_vc
         
         real(dp), pointer :: xa_removed(:) ! mass fractions for removed mass (1:species)
         real(dp) :: h1_czb_mass
         real(dp) :: total_angular_momentum
         real(dp) :: total_abs_angular_momentum
         real(dp) :: angular_momentum_removed
         real(dp) :: angular_momentum_added
         real(dp) :: residual_norm
         real(dp) :: max_residual
         real(dp) :: initial_timestep    
         real(dp) :: initial_L_center
         real(dp) :: initial_R_center
         real(dp) :: initial_v_center
         real(dp) :: tau_center  
         real(dp) :: super_eddington_wind_mdot
         real(dp) :: min_conv_time_scale, max_conv_time_scale, max_QHSE_time_scale

         integer :: num_solver_iterations ! for this step
         integer :: num_diffusion_solver_iters
         integer :: num_diffusion_solver_steps
         integer :: num_rotation_solver_steps
         integer :: result_reason

         logical :: just_wrote_terminal_header
         logical :: doing_first_model_of_run
         logical :: doing_first_model_after_restart
         logical :: was_in_implicit_wind_limit_old
      
         integer :: Tlim_dXnuc_species
         integer :: Tlim_dXnuc_cell
         
         integer :: Tlim_dXnuc_drop_species
         integer :: Tlim_dXnuc_drop_cell
         
         integer :: Tlim_dX_species
         integer :: Tlim_dX_cell
         
         integer :: Tlim_dt_div_min_dr_div_cs_cell
         
         integer :: Tlim_dX_div_X_species
         integer :: Tlim_dX_div_X_cell
         
         integer :: Tlim_dlgL_nuc_category
         integer :: Tlim_dlgL_nuc_cell
   
         logical :: need_to_update_history_now
         logical :: need_to_save_profiles_now
         integer :: save_profiles_model_priority
         
         ! extras -- for saving extra info along with model data for restarts
         integer :: len_extra_iwork
         integer :: len_extra_work
         integer, pointer, dimension(:) :: extra_iwork, extra_iwork_old
         real(dp), pointer, dimension(:) :: extra_work, extra_work_old
         
         ! more extras -- for run_star_extras applications with state.
         ! if your extra state info is all held in the following variables,
         ! restarts and retries will work without doing anything special.

         integer, dimension(star_num_xtra_vals) :: ixtra, ixtra_old
         real(dp), dimension(star_num_xtra_vals) :: xtra, xtra_old
         logical, dimension(star_num_xtra_vals) :: lxtra, lxtra_old


      ! some extra arrays for use in run_star_extra.f when you need an
      ! array which is guaranteed to have a length of s% nz
      ! these arrays are NOT automatically managed during restarts

         integer, pointer, dimension(:) :: &
            ixtra1_array, ixtra2_array, ixtra3_array, &
            ixtra4_array, ixtra5_array, ixtra6_array
         real(dp), pointer, dimension(:) :: &
            xtra1_array, xtra2_array, xtra3_array, &
            xtra4_array, xtra5_array, xtra6_array


      ! equation residuals, etc
      
         ! equ(i,k) is residual for equation i of cell k            
         real(dp), dimension(:,:), pointer :: equ ! (nvar,nz);  equ => equ1
         real(dp), dimension(:), pointer :: equ1 ! (nvar*nz); data for equ
         real(dp), dimension(:,:,:), pointer :: ublk, dblk, lblk ! (nvar,nvar,nz) 
         ! dblk(i,j,k) = dequ(i,k)/dx(j,k)
         ! lblk(i,j,k) = dequ(i,k)/dx(j,k-1)
         ! ublk(i,j,k) = dequ(i,k)/dx(j,k+1)
         
         real(dp), pointer, dimension(:) :: AF1 ! for factored jacobian

         real(dp), pointer :: solver_dx(:,:), solver_dx1(:)
         integer :: solver_iter, solver_adjust_iter
         real(dp), pointer :: x_scale(:,:), x_scale1(:)

         ! residual_weight and correction_weight are used by the solver in deciding if have a good enough solution.
         ! these are = 1d0 in most case, but are sometimes reduced to tell the solver "don't worry, be happy."
         ! weight = 0d0 means just ignore this item in making decisions.
         real(dp), dimension(:,:), pointer :: residual_weight, correction_weight ! (nvar,nz)
         real(dp), dimension(:), pointer :: residual_weight1, correction_weight1 ! (nvar,nz)

         
      ! sizeB info
         integer :: max_zone1, max_var1, max_zone2, max_var2
         real(dp) :: abs_max_corr1, abs_max_corr2
      
      ! bcyclic solver
         type(ulstore), pointer :: bcyclic_odd_storage(:)
         type(ulstore_qp), pointer :: bcyclic_odd_storage_qp(:)
      
      ! misc
         
         character (len=strlen) :: retry_message
         integer :: retry_message_k

         real(dp), pointer :: eps_pre_mix(:) ! (nz)
         real(dp), pointer, dimension(:) :: burn_avg_epsnuc
         integer, pointer, dimension(:) :: burn_num_iters

         integer :: termination_code
         integer :: dX_nuc_drop_max_k
         integer :: dX_nuc_drop_max_j
         integer :: solver_test_partials_var
         integer :: solver_test_partials_dx_sink
         
         logical :: doing_check_partials
         logical :: have_mixing_info
         logical :: doing_solver_iterations
         logical :: doing_flash_wind
         logical :: doing_rlo_wind
         
         real(dp) :: save1, save2, save3, save4, save5, save6, save7
         
         real(dp) :: dt_next_unclipped
         real(dp) :: L_for_BB_outer_BC
         real(dp) :: solver_test_partials_val, solver_test_partials_dval_dx

         ! variation control value for timestep controller
         real(dp) :: dt_limit_ratio_old
      
         ! surface info
         real(dp) :: L_surf_old ! surface luminosity (Lsun units)
         real(dp) :: L_phot_old ! photosphere luminosity (Lsun units)
         
         ! convection regions
         integer :: n_conv_regions_old
         real(dp), dimension(max_num_mixing_regions) :: &
            cz_bot_mass_old, cz_top_mass_old
      
         real(dp), pointer, dimension(:) :: max_abs_xa_corr ! from most recent solver iteration

         real(dp) :: total_angular_momentum_old
         real(dp) :: Teff_old

         ! For eps_mdot  (source term only used when mdot /= 0)
         real(dp), pointer, dimension(:) :: eps_mdot
         real(dp), pointer, dimension(:)  :: dm_before_adjust_mass
         real(dp), pointer, dimension(:)  :: total_energy_profile_before_adjust_mass
         real(dp), pointer, dimension(:)  :: total_energy_profile_after_adjust_mass
         ! End of eps_mdot

         logical :: eps_grav_form_for_energy_eqn
         real(dp), pointer, dimension(:)  :: &
            dL_dm, energy_sources, energy_others, dwork_dm, dkedt, dpedt, dedt, detrbdt
            
         real(dp) :: &
            total_internal_energy_old, &
            total_gravitational_energy_old, &
            total_radial_kinetic_energy_old, &
            total_turbulent_energy_old, &
            total_rotational_kinetic_energy_old, &
            total_energy_old, &
         
            total_internal_energy_before_adjust_mass, &
            total_gravitational_energy_before_adjust_mass, &
            total_radial_kinetic_energy_before_adjust_mass, &
            total_turbulent_energy_before_adjust_mass, &
            total_rotational_kinetic_energy_before_adjust_mass, &
            total_energy_before_adjust_mass, &      
         
            total_internal_energy_after_adjust_mass, &
            total_gravitational_energy_after_adjust_mass, &
            total_radial_kinetic_energy_after_adjust_mass, &
            total_turbulent_energy_after_adjust_mass, &
            total_rotational_kinetic_energy_after_adjust_mass, &
            total_energy_after_adjust_mass, &      
            
            total_internal_energy_start, &
            total_gravitational_energy_start, &
            total_radial_kinetic_energy_start, &
            total_turbulent_energy_start, &
            total_rotational_kinetic_energy_start, &
            total_energy_start, &      
         
            total_internal_energy_end, &
            total_gravitational_energy_end, &
            total_radial_kinetic_energy_end, &
            total_turbulent_energy_end, &
            total_rotational_kinetic_energy_end, &
            total_energy_end, &      
             
            total_internal_energy, &
            total_gravitational_energy, &
            total_radial_kinetic_energy, &
            total_turbulent_energy, &
            total_rotational_kinetic_energy, &
            total_energy, &
            
            surface_cell_specific_total_energy_old, &

            total_energy_change_from_mdot, & 
            energy_change_from_do_adjust_mass_and_calculate_eps_mdot, &    
            
            total_energy_sources_and_sinks, &
            virial_thm_P_avg, &
            total_eps_grav, &
            total_eps_mdot, &
            work_outward_at_surface, &
            work_inward_at_center, &
            non_epsnuc_energy_change_from_split_burn, &
            total_nuclear_heating, &
            total_non_nuc_neu_cooling, &
            total_irradiation_heating, &
            total_WD_sedimentation_heating, &
            total_energy_from_diffusion, &
            total_extra_heating, &
            error_in_energy_conservation
         
         real(dp) :: &
            cumulative_energy_error_old, &
            cumulative_extra_heating_old
            
         real(dp) :: mesh_adjust_IE_conservation
         real(dp) :: mesh_adjust_PE_conservation
         real(dp) :: mesh_adjust_KE_conservation
         real(dp) :: mesh_adjust_Eturb_conservation
         
         character (len=strlen) :: most_recent_photo_name
         
         real(dp) :: &
            total_internal_energy_initial, &
            total_gravitational_energy_initial, &
            total_radial_kinetic_energy_initial, &
            total_turbulent_energy_initial, &
            total_rotational_kinetic_energy_initial
         
         logical :: using_revised_net_name
         character (len=net_name_len) :: revised_net_name
                  
         integer :: num_extra_profile_cols
         character (len=80), pointer :: extra_profile_col_names(:)
         real(dp), pointer :: extra_profile_col_vals(:,:) ! (nz,num_extra)

      ! system clock 
         integer(8) :: starting_system_clock_time
         integer(8) :: current_system_clock_time
         integer(8) :: system_clock_rate
         real(dp) :: total_elapsed_time ! (current - starting)/rate
         integer(8) :: system_clock_at_start_of_step
         
      ! timing information
         logical :: doing_timing

!         use star_utils, only: start_time, update_time
!         integer(8) :: time0
!         real(dp) :: total
!         if (s% doing_timing) call start_time(s, time0, total)
!         <work>
!         if (s% doing_timing) call update_time(s, time0, total, s% time_<work>)
         
         ! if change this list, also must change
                  ! total_times function in star_utils
                  ! and init to set the time to 0.0 at start.
         real(dp) :: &
            time_evolve_step, &
            time_remesh, &
            time_adjust_mass, &
            time_conv_premix, &
            time_element_diffusion, &
            time_struct_burn_mix, &
            time_solver_matrix, &
            time_solve_mix, &
            time_solve_burn, &
            time_solve_omega_mix, &
            time_eos, &
            time_neu_kap, &
            time_nonburn_net, &
            time_mlt, &
            time_set_hydro_vars, &
            time_set_mixing_info, &
            time_total
            ! change total_times function in star_utils if change this list
         integer :: timing_num_get_eos_calls, timing_num_solve_eos_calls, timing_num_get_kap_calls
         
         
         
           
