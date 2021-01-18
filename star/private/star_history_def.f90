! ***********************************************************************
!
!   Copyright (C) 2014-2019  Bill Paxton & The MESA Team
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

      module star_history_def

      use star_def

      implicit none
   ! history column options

      integer, parameter :: h_model_number = 1
      integer, parameter :: h_log_star_age = h_model_number + 1
      integer, parameter :: h_star_age = h_log_star_age + 1
      integer, parameter :: h_log_star_mass = h_star_age + 1
      integer, parameter :: h_star_mass = h_log_star_mass + 1
      integer, parameter :: h_delta_mass = h_star_mass + 1
      integer, parameter :: h_log_xmstar = h_delta_mass + 1
      integer, parameter :: h_star_mdot = h_log_xmstar + 1
      integer, parameter :: h_log_abs_mdot = h_star_mdot + 1
      integer, parameter :: h_time_step = h_log_abs_mdot + 1
      integer, parameter :: h_e_thermal = h_time_step + 1
      integer, parameter :: h_species = h_e_thermal + 1
      integer, parameter :: h_Tsurf_factor = h_species + 1
      integer, parameter :: h_tau_factor = h_Tsurf_factor + 1
      integer, parameter :: h_log_tau_center = h_tau_factor + 1
      integer, parameter :: h_tau_surface = h_log_tau_center + 1
      integer, parameter :: h_num_zones = h_tau_surface + 1

      integer, parameter :: h_star_age_sec = h_num_zones + 1
      integer, parameter :: h_star_age_min = h_star_age_sec + 1
      integer, parameter :: h_star_age_hr = h_star_age_min + 1
      integer, parameter :: h_day = h_star_age_hr + 1
      integer, parameter :: h_star_age_day = h_day + 1
      integer, parameter :: h_log_star_age_sec = h_star_age_day + 1
      integer, parameter :: h_time_step_sec = h_log_star_age_sec + 1
      integer, parameter :: h_log_dt_sec = h_time_step_sec + 1
      integer, parameter :: h_time_step_days = h_log_dt_sec + 1
      integer, parameter :: h_log_dt_days = h_time_step_days + 1

      integer, parameter :: h_m_center_gm = h_log_dt_days + 1
      integer, parameter :: h_r_center_km = h_m_center_gm + 1
      integer, parameter :: h_r_center_cm = h_r_center_km + 1
      integer, parameter :: h_m_center = h_r_center_cm + 1
      integer, parameter :: h_r_center = h_m_center + 1
      integer, parameter :: h_L_center = h_r_center + 1
      integer, parameter :: h_log_L_center_ergs_s = h_L_center + 1
      integer, parameter :: h_log_L_center = h_log_L_center_ergs_s + 1
      integer, parameter :: h_infall_div_cs = h_log_L_center + 1
      integer, parameter :: h_v_center_kms = h_infall_div_cs + 1
      integer, parameter :: h_v_center = h_v_center_kms + 1

      integer, parameter :: h_dlnR_dlnM = h_v_center + 1
      integer, parameter :: h_mdot_timescale = h_dlnR_dlnM + 1
      integer, parameter :: h_kh_div_mdot_timescales = h_mdot_timescale + 1

      integer, parameter :: h_star_gravitational_mass = h_kh_div_mdot_timescales + 1
      integer, parameter :: h_star_mass_grav_div_mass = h_star_gravitational_mass + 1

      integer, parameter :: h_conv_mx1_top = h_star_mass_grav_div_mass + 1
      integer, parameter :: h_conv_mx1_bot = h_conv_mx1_top + 1
      integer, parameter :: h_conv_mx2_top = h_conv_mx1_bot + 1
      integer, parameter :: h_conv_mx2_bot = h_conv_mx2_top + 1
      integer, parameter :: h_mx1_top = h_conv_mx2_bot + 1
      integer, parameter :: h_mx1_bot = h_mx1_top + 1
      integer, parameter :: h_mx2_top = h_mx1_bot + 1
      integer, parameter :: h_mx2_bot = h_mx2_top + 1

      integer, parameter :: h_conv_mx1_top_r = h_mx2_bot + 1
      integer, parameter :: h_conv_mx1_bot_r = h_conv_mx1_top_r + 1
      integer, parameter :: h_conv_mx2_top_r = h_conv_mx1_bot_r + 1
      integer, parameter :: h_conv_mx2_bot_r = h_conv_mx2_top_r + 1
      integer, parameter :: h_mx1_top_r = h_conv_mx2_bot_r + 1
      integer, parameter :: h_mx1_bot_r = h_mx1_top_r + 1
      integer, parameter :: h_mx2_top_r = h_mx1_bot_r + 1
      integer, parameter :: h_mx2_bot_r = h_mx2_top_r + 1

      integer, parameter :: h_mix_relr_regions = h_mx2_bot_r + 1
      integer, parameter :: h_mixing_regions = h_mix_relr_regions + 1
      integer, parameter :: h_epsnuc_M_1 = h_mixing_regions + 1
      integer, parameter :: h_epsnuc_M_2 = h_epsnuc_M_1 + 1
      integer, parameter :: h_epsnuc_M_3 = h_epsnuc_M_2 + 1
      integer, parameter :: h_epsnuc_M_4 = h_epsnuc_M_3 + 1
      integer, parameter :: h_epsnuc_M_5 = h_epsnuc_M_4 + 1
      integer, parameter :: h_epsnuc_M_6 = h_epsnuc_M_5 + 1
      integer, parameter :: h_epsnuc_M_7 = h_epsnuc_M_6 + 1
      integer, parameter :: h_epsnuc_M_8 = h_epsnuc_M_7 + 1
      integer, parameter :: h_burning_regions = h_epsnuc_M_8 + 1
      integer, parameter :: h_burn_relr_regions = h_burning_regions + 1

      integer, parameter :: h_power_h_burn = h_burn_relr_regions + 1
      integer, parameter :: h_power_he_burn = h_power_h_burn + 1

      integer, parameter :: h_h_rich_layer_mass = h_power_he_burn + 1
      integer, parameter :: h_he_rich_layer_mass = h_h_rich_layer_mass + 1
      integer, parameter :: h_c_rich_layer_mass = h_he_rich_layer_mass + 1
      integer, parameter :: h_o_rich_layer_mass = h_c_rich_layer_mass + 1
      integer, parameter :: h_si_rich_layer_mass = h_o_rich_layer_mass + 1

      integer, parameter :: h_he_core_mass = h_si_rich_layer_mass + 1
      integer, parameter :: h_he_core_radius = h_he_core_mass + 1
      integer, parameter :: h_he_core_lgT = h_he_core_radius + 1
      integer, parameter :: h_he_core_lgRho = h_he_core_lgT + 1
      integer, parameter :: h_he_core_L = h_he_core_lgRho + 1
      integer, parameter :: h_he_core_v = h_he_core_L + 1
      integer, parameter :: h_he_core_omega = h_he_core_v + 1
      integer, parameter :: h_he_core_omega_div_omega_crit = h_he_core_omega + 1
      integer, parameter :: h_he_core_k = h_he_core_omega_div_omega_crit + 1

      integer, parameter :: h_c_core_mass = h_he_core_k + 1
      integer, parameter :: h_c_core_radius = h_c_core_mass + 1
      integer, parameter :: h_c_core_lgT = h_c_core_radius + 1
      integer, parameter :: h_c_core_lgRho = h_c_core_lgT + 1
      integer, parameter :: h_c_core_L = h_c_core_lgRho + 1
      integer, parameter :: h_c_core_v = h_c_core_L + 1
      integer, parameter :: h_c_core_omega = h_c_core_v + 1
      integer, parameter :: h_c_core_omega_div_omega_crit = h_c_core_omega + 1
      integer, parameter :: h_c_core_k = h_c_core_omega_div_omega_crit + 1

      integer, parameter :: h_o_core_mass = h_c_core_k + 1
      integer, parameter :: h_o_core_radius = h_o_core_mass + 1
      integer, parameter :: h_o_core_lgT = h_o_core_radius + 1
      integer, parameter :: h_o_core_lgRho = h_o_core_lgT + 1
      integer, parameter :: h_o_core_L = h_o_core_lgRho + 1
      integer, parameter :: h_o_core_v = h_o_core_L + 1
      integer, parameter :: h_o_core_omega = h_o_core_v + 1
      integer, parameter :: h_o_core_omega_div_omega_crit = h_o_core_omega + 1
      integer, parameter :: h_o_core_k = h_o_core_omega_div_omega_crit + 1

      integer, parameter :: h_si_core_mass = h_o_core_k + 1
      integer, parameter :: h_si_core_radius = h_si_core_mass + 1
      integer, parameter :: h_si_core_lgT = h_si_core_radius + 1
      integer, parameter :: h_si_core_lgRho = h_si_core_lgT + 1
      integer, parameter :: h_si_core_L = h_si_core_lgRho + 1
      integer, parameter :: h_si_core_v = h_si_core_L + 1
      integer, parameter :: h_si_core_omega = h_si_core_v + 1
      integer, parameter :: h_si_core_omega_div_omega_crit = h_si_core_omega + 1
      integer, parameter :: h_si_core_k = h_si_core_omega_div_omega_crit + 1

      integer, parameter :: h_fe_core_mass = h_si_core_k + 1
      integer, parameter :: h_fe_core_radius = h_fe_core_mass + 1
      integer, parameter :: h_fe_core_lgT = h_fe_core_radius + 1
      integer, parameter :: h_fe_core_lgRho = h_fe_core_lgT + 1
      integer, parameter :: h_fe_core_L = h_fe_core_lgRho + 1
      integer, parameter :: h_fe_core_v = h_fe_core_L + 1
      integer, parameter :: h_fe_core_omega = h_fe_core_v + 1
      integer, parameter :: h_fe_core_omega_div_omega_crit = h_fe_core_omega + 1
      integer, parameter :: h_fe_core_k = h_fe_core_omega_div_omega_crit + 1

      integer, parameter :: h_neutron_rich_core_mass = h_fe_core_k + 1
      integer, parameter :: h_neutron_rich_core_radius = h_neutron_rich_core_mass + 1
      integer, parameter :: h_neutron_rich_core_lgT = h_neutron_rich_core_radius + 1
      integer, parameter :: h_neutron_rich_core_lgRho = h_neutron_rich_core_lgT + 1
      integer, parameter :: h_neutron_rich_core_L = h_neutron_rich_core_lgRho + 1
      integer, parameter :: h_neutron_rich_core_v = h_neutron_rich_core_L + 1
      integer, parameter :: h_neutron_rich_core_omega = h_neutron_rich_core_v + 1
      integer, parameter :: h_neutron_rich_core_omega_div_omega_crit = h_neutron_rich_core_omega + 1
      integer, parameter :: h_neutron_rich_core_k = h_neutron_rich_core_omega_div_omega_crit + 1

      integer, parameter :: h_log_max_T = h_neutron_rich_core_k + 1
      integer, parameter :: h_log_cntr_dr_cm = h_log_max_T + 1
      integer, parameter :: h_log_cntr_T = h_log_cntr_dr_cm + 1
      integer, parameter :: h_log_center_T = h_log_cntr_T + 1
      integer, parameter :: h_log_cntr_Rho = h_log_center_T + 1
      integer, parameter :: h_log_center_Rho = h_log_cntr_Rho + 1
      integer, parameter :: h_log_cntr_P = h_log_center_Rho + 1
      integer, parameter :: h_log_center_P = h_log_cntr_P + 1

      integer, parameter :: h_max_T = h_log_center_P + 1
      integer, parameter :: h_center_T = h_max_T + 1
      integer, parameter :: h_center_Rho = h_center_T + 1
      integer, parameter :: h_center_P = h_center_Rho + 1

      integer, parameter :: h_log_mesh_adjust_IE_conservation = h_center_P + 1
      integer, parameter :: h_log_mesh_adjust_PE_conservation = h_log_mesh_adjust_IE_conservation + 1
      integer, parameter :: h_log_mesh_adjust_KE_conservation = h_log_mesh_adjust_PE_conservation + 1

      integer, parameter :: h_avg_abs_v_div_cs = h_log_mesh_adjust_KE_conservation + 1
      integer, parameter :: h_log_avg_abs_v_div_cs = h_avg_abs_v_div_cs + 1
      integer, parameter :: h_max_abs_v_div_cs = h_log_avg_abs_v_div_cs + 1
      integer, parameter :: h_log_max_abs_v_div_cs = h_max_abs_v_div_cs + 1

      integer, parameter :: h_avg_abs_v = h_log_max_abs_v_div_cs + 1
      integer, parameter :: h_log_avg_abs_v = h_avg_abs_v + 1
      integer, parameter :: h_max_abs_v = h_log_avg_abs_v + 1
      integer, parameter :: h_log_max_abs_v = h_max_abs_v + 1

      integer, parameter :: h_total_IE_plus_KE_start = h_log_max_abs_v + 1
      integer, parameter :: h_log_total_IE_plus_KE = h_total_IE_plus_KE_start + 1
      integer, parameter :: h_total_IE_plus_KE = h_log_total_IE_plus_KE + 1

      integer, parameter :: h_total_energy_plus_L_surf = h_total_IE_plus_KE + 1

      integer, parameter :: h_total_internal_energy_after_adjust_mass = h_total_energy_plus_L_surf + 1
      integer, parameter :: h_total_gravitational_energy_after_adjust_mass = &
         h_total_internal_energy_after_adjust_mass + 1
      integer, parameter :: h_total_turbulent_energy_after_adjust_mass = &
         h_total_gravitational_energy_after_adjust_mass + 1
      integer, parameter :: h_total_radial_kinetic_energy_after_adjust_mass = &
         h_total_turbulent_energy_after_adjust_mass + 1
      integer, parameter :: h_total_rotational_kinetic_energy_after_adjust_mass = &
         h_total_radial_kinetic_energy_after_adjust_mass + 1
      integer, parameter :: h_total_energy_after_adjust_mass = &
         h_total_rotational_kinetic_energy_after_adjust_mass + 1

      integer, parameter :: h_total_internal_energy = h_total_energy_after_adjust_mass + 1
      integer, parameter :: h_total_gravitational_energy = h_total_internal_energy + 1
      integer, parameter :: h_total_turbulent_energy = h_total_gravitational_energy + 1
      integer, parameter :: h_total_radial_kinetic_energy = h_total_turbulent_energy + 1
      integer, parameter :: h_total_rotational_kinetic_energy = h_total_radial_kinetic_energy + 1
      integer, parameter :: h_total_energy_foe = h_total_rotational_kinetic_energy + 1
      integer, parameter :: h_total_energy = h_total_energy_foe + 1

      integer, parameter :: h_log_total_internal_energy = h_total_energy + 1
      integer, parameter :: h_log_total_gravitational_energy = h_log_total_internal_energy + 1
      integer, parameter :: h_log_total_turbulent_energy = h_log_total_gravitational_energy + 1
      integer, parameter :: h_log_total_radial_kinetic_energy = h_log_total_turbulent_energy + 1
      integer, parameter :: h_log_total_rotational_kinetic_energy = h_log_total_radial_kinetic_energy + 1
      integer, parameter :: h_log_total_energy = h_log_total_rotational_kinetic_energy + 1
      
      integer, parameter :: h_total_IE_div_IE_plus_KE = h_log_total_energy + 1

      integer, parameter :: h_total_entropy = h_total_IE_div_IE_plus_KE + 1

      integer, parameter :: h_rms_dvdt_div_v = h_total_entropy + 1

      integer, parameter :: h_virial_thm_P_avg = h_rms_dvdt_div_v + 1
      integer, parameter :: h_virial_thm_rel_err = h_virial_thm_P_avg + 1
      integer, parameter :: h_total_eps_grav = h_virial_thm_rel_err + 1
      integer, parameter :: h_work_outward_at_surface = h_total_eps_grav + 1
      integer, parameter :: h_work_inward_at_center = h_work_outward_at_surface + 1
      integer, parameter :: h_total_nuclear_heating = h_work_inward_at_center + 1
      integer, parameter :: h_total_non_nuc_neu_cooling = h_total_nuclear_heating + 1
      integer, parameter :: h_total_WD_sedimentation_heating = h_total_non_nuc_neu_cooling + 1
      integer, parameter :: h_total_irradiation_heating = h_total_WD_sedimentation_heating + 1
      integer, parameter :: h_total_extra_heating = h_total_irradiation_heating + 1

      integer, parameter :: h_total_energy_sources_and_sinks = h_total_extra_heating + 1
      integer, parameter :: h_log_rel_error_in_energy_conservation = h_total_energy_sources_and_sinks + 1
      integer, parameter :: h_rel_error_in_energy_conservation = h_log_rel_error_in_energy_conservation + 1
      integer, parameter :: h_error_in_energy_conservation = h_rel_error_in_energy_conservation + 1

      integer, parameter :: h_cumulative_energy_error = h_error_in_energy_conservation + 1
      integer, parameter :: h_rel_cumulative_energy_error = h_cumulative_energy_error + 1
      integer, parameter :: h_abs_rel_E_err = h_rel_cumulative_energy_error + 1
      integer, parameter :: h_log_rel_E_err = h_abs_rel_E_err + 1
      integer, parameter :: h_tot_E_equ_err = h_log_rel_E_err + 1
      integer, parameter :: h_tot_E_err = h_tot_E_equ_err + 1
      integer, parameter :: h_rel_E_err = h_tot_E_err + 1
      integer, parameter :: h_rel_run_E_err = h_rel_E_err + 1
      integer, parameter :: h_log_rel_run_E_err = h_rel_run_E_err + 1
      integer, parameter :: h_log_rel_cumulative_energy_error = h_log_rel_run_E_err + 1

      integer, parameter :: h_log_residual_norm = h_log_rel_cumulative_energy_error + 1
      integer, parameter :: h_log_max_residual = h_log_residual_norm + 1

      integer, parameter :: h_log_max_dvdt_residual = h_log_max_residual + 1
      integer, parameter :: h_log_max_drdt_residual = h_log_max_dvdt_residual + 1
      integer, parameter :: h_log_max_lnd_residual = h_log_max_drdt_residual + 1
      integer, parameter :: h_log_max_dEdt_residual = h_log_max_lnd_residual + 1

      integer, parameter :: h_max_abs_v_residual = h_log_max_dEdt_residual + 1
      integer, parameter :: h_log_max_abs_v_residual = h_max_abs_v_residual + 1
      integer, parameter :: h_avg_v_residual = h_log_max_abs_v_residual + 1
      integer, parameter :: h_log_avg_v_residual = h_avg_v_residual + 1

      integer, parameter :: h_log_max_abs_E_residual = h_log_avg_v_residual + 1
      integer, parameter :: h_log_avg_E_residual = h_log_max_abs_E_residual + 1
      integer, parameter :: h_max_abs_E_residual = h_log_avg_E_residual + 1
      integer, parameter :: h_avg_E_residual = h_max_abs_E_residual + 1

      integer, parameter :: h_u_surf_km_s = h_avg_E_residual + 1
      integer, parameter :: h_u_surf = h_u_surf_km_s + 1
      integer, parameter :: h_u_div_csound_max = h_u_surf + 1
      integer, parameter :: h_u_div_csound_surf = h_u_div_csound_max + 1

      integer, parameter :: h_center_zbar = h_u_div_csound_surf + 1
      integer, parameter :: h_center_abar = h_center_zbar + 1
      integer, parameter :: h_center_mu = h_center_abar + 1
      integer, parameter :: h_center_ye = h_center_mu + 1
      integer, parameter :: h_max_entropy = h_center_ye + 1
      integer, parameter :: h_center_entropy = h_max_entropy + 1
      integer, parameter :: h_v_surf_div_escape_v = h_center_entropy + 1
      integer, parameter :: h_v_surf_km_s = h_v_surf_div_escape_v + 1
      integer, parameter :: h_v_surf = h_v_surf_km_s + 1
      integer, parameter :: h_v_surf_div_v_kh = h_v_surf + 1
      integer, parameter :: h_v_div_csound_max = h_v_surf_div_v_kh + 1
      integer, parameter :: h_v_div_csound_surf = h_v_div_csound_max + 1
      integer, parameter :: h_log_dt = h_v_div_csound_surf + 1
      integer, parameter :: h_log_LH = h_log_dt + 1
      integer, parameter :: h_log_LHe = h_log_LH + 1
      integer, parameter :: h_power_photo = h_log_LHe + 1
      integer, parameter :: h_Lnuc_photo = h_power_photo + 1
      
      integer, parameter :: h_Lsurf_m = h_Lnuc_photo + 1
      integer, parameter :: h_luminosity_ergs_s = h_Lsurf_m + 1
      integer, parameter :: h_log_L_ergs_s = h_luminosity_ergs_s + 1
      integer, parameter :: h_luminosity = h_log_L_ergs_s + 1
      integer, parameter :: h_log_L = h_luminosity + 1
      
      integer, parameter :: h_power_z_burn = h_log_L + 1
      integer, parameter :: h_log_LZ = h_power_z_burn + 1

      integer, parameter :: h_log_Lneu_nuc = h_log_LZ + 1
      integer, parameter :: h_log_Lneu_nonnuc = h_log_Lneu_nuc + 1
      integer, parameter :: h_log_Lneu = h_log_Lneu_nonnuc + 1
      integer, parameter :: h_log_R_cm = h_log_Lneu + 1
      integer, parameter :: h_radius_cm = h_log_R_cm + 1
      integer, parameter :: h_radius = h_radius_cm + 1
      integer, parameter :: h_log_R = h_radius + 1
      integer, parameter :: h_log_Teff = h_log_R + 1
      integer, parameter :: h_Teff = h_log_Teff + 1
      integer, parameter :: h_effective_T = h_Teff + 1
      integer, parameter :: h_gravity = h_effective_T + 1
      integer, parameter :: h_log_g = h_gravity + 1
      integer, parameter :: h_log_L_div_Ledd = h_log_g + 1
      integer, parameter :: h_lum_div_Ledd = h_log_L_div_Ledd + 1
      integer, parameter :: h_max_L_rad_div_Ledd_div_phi_Joss = h_lum_div_Ledd + 1
      integer, parameter :: h_max_L_rad_div_Ledd = h_max_L_rad_div_Ledd_div_phi_Joss + 1

      integer, parameter :: h_gamma1_min = h_max_L_rad_div_Ledd + 1
      integer, parameter :: h_logT_max = h_gamma1_min + 1
      integer, parameter :: h_logQ_max = h_logT_max + 1
      integer, parameter :: h_logQ_min = h_logQ_max + 1

      integer, parameter :: h_avg_skipped_setvars_per_step = h_logQ_min + 1
      integer, parameter :: h_avg_setvars_per_step = h_avg_skipped_setvars_per_step + 1
      integer, parameter :: h_avg_solver_setvars_per_step = h_avg_setvars_per_step + 1
      
      integer, parameter :: h_num_steps_skipped_1st_setvars = h_avg_solver_setvars_per_step + 1
      integer, parameter :: h_fraction_of_steps_skipped_1st_setvars = h_num_steps_skipped_1st_setvars + 1
      integer, parameter :: h_num_retries = h_fraction_of_steps_skipped_1st_setvars + 1
      integer, parameter :: h_h1_czb_mass = h_num_retries + 1
      integer, parameter :: h_surf_c12_minus_o16 = h_h1_czb_mass + 1
      integer, parameter :: h_surf_num_c12_div_num_o16 = h_surf_c12_minus_o16 + 1

      integer, parameter :: h_min_Pgas_div_P = h_surf_num_c12_div_num_o16 + 1
      integer, parameter :: h_log_center_eps_nuc = h_min_Pgas_div_P + 1
      integer, parameter :: h_d_center_eps_nuc_dlnT = h_log_center_eps_nuc + 1
      integer, parameter :: h_d_center_eps_nuc_dlnd = h_d_center_eps_nuc_dlnT + 1

      integer, parameter :: h_center_eps_nuc = h_d_center_eps_nuc_dlnd + 1
      integer, parameter :: h_center_non_nuc_neu = h_center_eps_nuc + 1

      integer, parameter :: h_center_dL_dm = h_center_non_nuc_neu + 1
      integer, parameter :: h_center_eps_grav = h_center_dL_dm + 1
      integer, parameter :: h_center_degeneracy = h_center_eps_grav + 1
      integer, parameter :: h_center_gamma = h_center_degeneracy + 1

      integer, parameter :: h_center_dlogT = h_center_gamma + 1
      integer, parameter :: h_center_dlogRho = h_center_dlogT + 1

      integer, parameter :: h_center_dlnT_dt = h_center_dlogRho + 1
      integer, parameter :: h_center_dlnd_dt = h_center_dlnT_dt + 1

      integer, parameter :: h_envelope_mass = h_center_dlnd_dt + 1
      integer, parameter :: h_envelope_fraction_left = h_envelope_mass + 1

      integer, parameter :: h_tau10_mass = h_envelope_fraction_left + 1
      integer, parameter :: h_tau10_radius = h_tau10_mass + 1
      integer, parameter :: h_tau10_lgP = h_tau10_radius + 1
      integer, parameter :: h_tau10_T = h_tau10_lgP + 1
      integer, parameter :: h_tau10_lgT = h_tau10_T + 1
      integer, parameter :: h_tau10_lgRho = h_tau10_lgT + 1
      integer, parameter :: h_tau10_L = h_tau10_lgRho + 1
      integer, parameter :: h_tau100_mass = h_tau10_L + 1
      integer, parameter :: h_tau100_radius = h_tau100_mass + 1
      integer, parameter :: h_tau100_lgP = h_tau100_radius + 1
      integer, parameter :: h_tau100_T = h_tau100_lgP + 1
      integer, parameter :: h_tau100_lgT = h_tau100_T + 1
      integer, parameter :: h_tau100_lgRho = h_tau100_lgT + 1
      integer, parameter :: h_tau100_L = h_tau100_lgRho + 1
      integer, parameter :: h_dynamic_timescale = h_tau100_L + 1
      integer, parameter :: h_kh_timescale = h_dynamic_timescale + 1
      integer, parameter :: h_nuc_timescale = h_kh_timescale + 1
      integer, parameter :: h_log_abs_Lgrav = h_nuc_timescale + 1
      integer, parameter :: h_eps_grav_integral = h_log_abs_Lgrav + 1
      integer, parameter :: h_log_extra_L = h_eps_grav_integral + 1
      integer, parameter :: h_extra_L = h_log_extra_L + 1
      integer, parameter :: h_log_power_nuc_burn = h_extra_L + 1
      integer, parameter :: h_log_Lnuc_ergs_s = h_log_power_nuc_burn + 1
      integer, parameter :: h_log_Lnuc = h_log_Lnuc_ergs_s + 1
      integer, parameter :: h_mass_ext_to_max_eps_nuc = h_log_Lnuc + 1
      integer, parameter :: h_mass_loc_of_max_eps_nuc = h_mass_ext_to_max_eps_nuc + 1

      integer, parameter :: h_diffusion_time_H_He_bdy = h_mass_loc_of_max_eps_nuc + 1
      integer, parameter :: h_temperature_H_He_bdy = h_diffusion_time_H_He_bdy + 1

      integer, parameter :: h_max_abs_v_velocity = h_temperature_H_He_bdy + 1
      integer, parameter :: h_max_abs_v_csound = h_max_abs_v_velocity + 1
      integer, parameter :: h_max_abs_v_v_div_cs = h_max_abs_v_csound + 1
      integer, parameter :: h_max_abs_v_lgT = h_max_abs_v_v_div_cs + 1
      integer, parameter :: h_max_abs_v_lgRho = h_max_abs_v_lgT + 1
      integer, parameter :: h_max_abs_v_lgP = h_max_abs_v_lgRho + 1
      integer, parameter :: h_max_abs_v_mass = h_max_abs_v_lgP + 1
      integer, parameter :: h_max_abs_v_radius = h_max_abs_v_mass + 1
      integer, parameter :: h_max_abs_v_radius_cm = h_max_abs_v_radius + 1
      integer, parameter :: h_max_abs_v_lgR = h_max_abs_v_radius_cm + 1
      integer, parameter :: h_max_abs_v_lgR_cm = h_max_abs_v_lgR + 1
      integer, parameter :: h_max_abs_v_L = h_max_abs_v_lgR_cm + 1
      integer, parameter :: h_max_abs_v_gamma1 = h_max_abs_v_L + 1
      integer, parameter :: h_max_abs_v_entropy = h_max_abs_v_gamma1 + 1
      integer, parameter :: h_max_abs_v_E0 = h_max_abs_v_entropy + 1
      integer, parameter :: h_max_abs_v_eps_nuc = h_max_abs_v_E0 + 1

      integer, parameter :: h_total_ni_co_56 = h_max_abs_v_eps_nuc + 1
      
      integer, parameter :: h_inner_mach1_velocity = h_total_ni_co_56 + 1
      integer, parameter :: h_inner_mach1_csound = h_inner_mach1_velocity + 1
      integer, parameter :: h_inner_mach1_v_div_cs = h_inner_mach1_csound + 1
      integer, parameter :: h_inner_mach1_lgT = h_inner_mach1_v_div_cs + 1
      integer, parameter :: h_inner_mach1_lgRho = h_inner_mach1_lgT + 1
      integer, parameter :: h_inner_mach1_lgP = h_inner_mach1_lgRho + 1
      integer, parameter :: h_inner_mach1_q = h_inner_mach1_lgP + 1
      integer, parameter :: h_inner_mach1_tau = h_inner_mach1_q + 1
      integer, parameter :: h_inner_mach1_mass = h_inner_mach1_tau + 1
      integer, parameter :: h_inner_mach1_radius = h_inner_mach1_mass + 1
      integer, parameter :: h_inner_mach1_gamma1 = h_inner_mach1_radius + 1
      integer, parameter :: h_inner_mach1_entropy = h_inner_mach1_gamma1 + 1
      integer, parameter :: h_inner_mach1_k = h_inner_mach1_entropy + 1

      integer, parameter :: h_outer_mach1_velocity = h_inner_mach1_k + 1
      integer, parameter :: h_outer_mach1_csound = h_outer_mach1_velocity + 1
      integer, parameter :: h_outer_mach1_v_div_cs = h_outer_mach1_csound + 1
      integer, parameter :: h_outer_mach1_lgT = h_outer_mach1_v_div_cs + 1
      integer, parameter :: h_outer_mach1_lgRho = h_outer_mach1_lgT + 1
      integer, parameter :: h_outer_mach1_lgP = h_outer_mach1_lgRho + 1
      integer, parameter :: h_outer_mach1_q = h_outer_mach1_lgP + 1
      integer, parameter :: h_outer_mach1_tau = h_outer_mach1_q + 1
      integer, parameter :: h_outer_mach1_mass = h_outer_mach1_tau + 1
      integer, parameter :: h_outer_mach1_radius = h_outer_mach1_mass + 1
      integer, parameter :: h_outer_mach1_gamma1 = h_outer_mach1_radius + 1
      integer, parameter :: h_outer_mach1_entropy = h_outer_mach1_gamma1 + 1
      integer, parameter :: h_outer_mach1_k = h_outer_mach1_entropy + 1

      integer, parameter :: h_shock_velocity = h_outer_mach1_k + 1
      integer, parameter :: h_shock_csound = h_shock_velocity + 1
      integer, parameter :: h_shock_v_div_cs = h_shock_csound + 1
      integer, parameter :: h_shock_lgT = h_shock_v_div_cs + 1
      integer, parameter :: h_shock_lgRho = h_shock_lgT + 1
      integer, parameter :: h_shock_lgP = h_shock_lgRho + 1
      integer, parameter :: h_shock_q = h_shock_lgP + 1
      integer, parameter :: h_shock_tau = h_shock_q + 1
      integer, parameter :: h_shock_mass_gm = h_shock_tau + 1
      integer, parameter :: h_shock_radius_cm = h_shock_mass_gm + 1
      integer, parameter :: h_shock_mass = h_shock_radius_cm + 1
      integer, parameter :: h_shock_radius = h_shock_mass + 1
      integer, parameter :: h_shock_gamma1 = h_shock_radius + 1
      integer, parameter :: h_shock_entropy = h_shock_gamma1 + 1
      integer, parameter :: h_shock_pre_lgRho = h_shock_entropy + 1
      integer, parameter :: h_shock_k = h_shock_pre_lgRho + 1

      integer, parameter :: h_trace_mass_location = h_shock_k + 1
      integer, parameter :: h_trace_mass_radius = h_trace_mass_location + 1
      integer, parameter :: h_trace_mass_lgT = h_trace_mass_radius + 1
      integer, parameter :: h_trace_mass_lgRho = h_trace_mass_lgT + 1
      integer, parameter :: h_trace_mass_L = h_trace_mass_lgRho + 1
      integer, parameter :: h_trace_mass_v = h_trace_mass_L + 1

      integer, parameter :: h_max_T_shell_binding_energy = h_trace_mass_v + 1
      integer, parameter :: h_max_T_lgP_thin_shell = h_max_T_shell_binding_energy + 1
      integer, parameter :: h_max_T_lgP = h_max_T_lgP_thin_shell + 1
      integer, parameter :: h_max_T_entropy = h_max_T_lgP + 1
      integer, parameter :: h_max_T_mass = h_max_T_entropy + 1
      integer, parameter :: h_max_T_radius = h_max_T_mass + 1
      integer, parameter :: h_max_T_lgT = h_max_T_radius + 1
      integer, parameter :: h_max_T_lgRho = h_max_T_lgT + 1
      integer, parameter :: h_max_T_L = h_max_T_lgRho + 1
      integer, parameter :: h_max_T_eps_nuc = h_max_T_L + 1

      integer, parameter :: h_max_conv_vel_div_csound = h_max_T_eps_nuc + 1
      integer, parameter :: h_max_gradT_div_grada = h_max_conv_vel_div_csound + 1
      integer, parameter :: h_max_gradT_sub_grada = h_max_gradT_div_grada + 1
      integer, parameter :: h_min_log_mlt_Gamma = h_max_gradT_sub_grada + 1

      integer, parameter :: h_dt_cell_collapse = h_min_log_mlt_Gamma + 1
      integer, parameter :: h_dt_div_dt_cell_collapse = h_dt_cell_collapse + 1

      integer, parameter :: h_min_dr_div_cs_k = h_dt_div_dt_cell_collapse + 1
      integer, parameter :: h_min_dr_div_cs = h_min_dr_div_cs_k + 1
      integer, parameter :: h_log_min_dr_div_cs = h_min_dr_div_cs + 1
      integer, parameter :: h_min_dr_div_cs_yr = h_log_min_dr_div_cs + 1
      integer, parameter :: h_log_min_dr_div_cs_yr = h_min_dr_div_cs_yr + 1
      integer, parameter :: h_dt_div_min_dr_div_cs = h_log_min_dr_div_cs_yr + 1
      integer, parameter :: h_log_dt_div_min_dr_div_cs = h_dt_div_min_dr_div_cs + 1

      integer, parameter :: h_surface_optical_depth = h_log_dt_div_min_dr_div_cs + 1
      integer, parameter :: h_log_surf_optical_depth = h_surface_optical_depth + 1

      integer, parameter :: h_log_surf_cell_opacity = h_log_surf_optical_depth + 1
      integer, parameter :: h_log_surf_cell_density = h_log_surf_cell_opacity + 1
      integer, parameter :: h_surface_cell_temperature = h_log_surf_cell_density + 1
      integer, parameter :: h_log_surf_cell_temperature = h_surface_cell_temperature + 1
      integer, parameter :: h_log_surf_cell_P = h_log_surf_cell_temperature + 1
      integer, parameter :: h_log_surf_cell_pressure = h_log_surf_cell_P + 1
      integer, parameter :: h_log_surf_cell_z = h_log_surf_cell_pressure + 1
      integer, parameter :: h_surface_cell_entropy = h_log_surf_cell_z + 1

      integer, parameter :: h_gradT_excess_min_beta = h_surface_cell_entropy + 1
      integer, parameter :: h_gradT_excess_max_lambda = h_gradT_excess_min_beta + 1
      integer, parameter :: h_gradT_excess_alpha = h_gradT_excess_max_lambda + 1

      integer, parameter :: h_log_Ledd = h_gradT_excess_alpha + 1
      integer, parameter :: h_compactness = h_log_Ledd + 1
      integer, parameter :: h_compactness_parameter = h_compactness + 1
      integer, parameter :: h_max_infall_speed = h_compactness_parameter + 1
      integer, parameter :: h_non_fe_core_rebound = h_max_infall_speed + 1
      integer, parameter :: h_non_fe_core_infall = h_non_fe_core_rebound + 1
      integer, parameter :: h_fe_core_infall = h_non_fe_core_infall + 1

      integer, parameter :: h_cz_bot_mass = h_fe_core_infall + 1
      integer, parameter :: h_cz_mass = h_cz_bot_mass + 1
      integer, parameter :: h_cz_log_xmsun = h_cz_mass + 1
      integer, parameter :: h_cz_xm = h_cz_log_xmsun + 1
      integer, parameter :: h_cz_log_xmass = h_cz_xm + 1
      integer, parameter :: h_cz_logT = h_cz_log_xmass + 1
      integer, parameter :: h_cz_logRho = h_cz_logT + 1
      integer, parameter :: h_cz_logP = h_cz_logRho + 1
      integer, parameter :: h_cz_log_column_depth = h_cz_logP + 1
      integer, parameter :: h_cz_log_radial_depth = h_cz_log_column_depth + 1
      integer, parameter :: h_cz_bot_radius = h_cz_log_radial_depth + 1
      integer, parameter :: h_cz_csound = h_cz_bot_radius + 1
      integer, parameter :: h_cz_scale_height = h_cz_csound + 1
      integer, parameter :: h_cz_grav = h_cz_scale_height + 1
      integer, parameter :: h_cz_log_eps_nuc = h_cz_grav + 1
      integer, parameter :: h_cz_t_heat = h_cz_log_eps_nuc + 1
      integer, parameter :: h_cz_eta = h_cz_t_heat + 1

      integer, parameter :: h_cz_log_tau = h_cz_eta + 1
      integer, parameter :: h_cz_opacity = h_cz_log_tau + 1
      integer, parameter :: h_cz_luminosity = h_cz_opacity + 1
      integer, parameter :: h_cz_zone = h_cz_luminosity + 1
      integer, parameter :: h_cz_omega = h_cz_zone + 1
      integer, parameter :: h_cz_omega_div_omega_crit = h_cz_omega + 1

      integer, parameter :: h_cz_top_mass = h_cz_omega_div_omega_crit + 1
      integer, parameter :: h_cz_top_log_xmsun = h_cz_top_mass + 1
      integer, parameter :: h_cz_top_xm = h_cz_top_log_xmsun + 1
      integer, parameter :: h_cz_top_log_xmass = h_cz_top_xm + 1
      integer, parameter :: h_cz_top_logT = h_cz_top_log_xmass + 1
      integer, parameter :: h_cz_top_logRho = h_cz_top_logT + 1
      integer, parameter :: h_cz_top_logP = h_cz_top_logRho + 1
      integer, parameter :: h_cz_top_log_column_depth = h_cz_top_logP + 1
      integer, parameter :: h_cz_top_log_radial_depth = h_cz_top_log_column_depth + 1
      integer, parameter :: h_cz_top_radius = h_cz_top_log_radial_depth + 1
      integer, parameter :: h_cz_top_csound = h_cz_top_radius + 1
      integer, parameter :: h_cz_top_scale_height = h_cz_top_csound + 1
      integer, parameter :: h_cz_top_grav = h_cz_top_scale_height + 1
      integer, parameter :: h_cz_top_log_eps_nuc = h_cz_top_grav + 1
      integer, parameter :: h_cz_top_t_heat = h_cz_top_log_eps_nuc + 1
      integer, parameter :: h_cz_top_eta = h_cz_top_t_heat + 1

      integer, parameter :: h_cz_top_log_tau = h_cz_top_eta + 1
      integer, parameter :: h_cz_top_opacity = h_cz_top_log_tau + 1
      integer, parameter :: h_cz_top_luminosity = h_cz_top_opacity + 1
      integer, parameter :: h_cz_top_zone = h_cz_top_luminosity + 1
      integer, parameter :: h_cz_top_omega = h_cz_top_zone + 1
      integer, parameter :: h_cz_top_omega_div_omega_crit = h_cz_top_omega + 1

      integer, parameter :: h_trace_mass_omega = h_cz_top_omega_div_omega_crit + 1
      integer, parameter :: h_trace_mass_omega_div_omega_crit = h_trace_mass_omega + 1

      integer, parameter :: h_kh_mdot_limit = h_trace_mass_omega_div_omega_crit + 1
      integer, parameter :: h_rotational_mdot_boost = h_kh_mdot_limit + 1
      integer, parameter :: h_log_rotational_mdot_boost = h_rotational_mdot_boost + 1

      integer, parameter :: h_logL_for_BB_outer_BC = h_log_rotational_mdot_boost + 1
      integer, parameter :: h_luminosity_for_BB_outer_BC = h_logL_for_BB_outer_BC + 1
      integer, parameter :: h_i_rot_total = h_luminosity_for_BB_outer_BC + 1
      integer, parameter :: h_surf_avg_j_rot = h_i_rot_total + 1
      integer, parameter :: h_surf_avg_omega = h_surf_avg_j_rot + 1
      integer, parameter :: h_surf_avg_omega_crit = h_surf_avg_omega + 1
      integer, parameter :: h_surf_avg_omega_div_omega_crit = h_surf_avg_omega_crit + 1

      integer, parameter :: h_surf_avg_v_rot = h_surf_avg_omega_div_omega_crit + 1
      integer, parameter :: h_surf_avg_v_crit = h_surf_avg_v_rot + 1
      integer, parameter :: h_surf_avg_v_div_v_crit = h_surf_avg_v_crit + 1

      integer, parameter :: h_v_wind_Km_per_s = h_surf_avg_v_div_v_crit + 1
      integer, parameter :: h_surf_escape_v = h_v_wind_Km_per_s + 1

      integer, parameter :: h_surf_avg_logT = h_surf_escape_v + 1
      integer, parameter :: h_surf_avg_logRho = h_surf_avg_logT + 1
      integer, parameter :: h_surf_avg_opacity = h_surf_avg_logRho + 1
      integer, parameter :: h_surf_avg_Lrad_div_Ledd = h_surf_avg_opacity + 1

      integer, parameter :: h_center_omega = h_surf_avg_Lrad_div_Ledd + 1
      integer, parameter :: h_surf_r_equatorial_div_r_polar = h_center_omega + 1
      integer, parameter :: h_surf_r_equatorial_div_r = h_surf_r_equatorial_div_r_polar + 1
      integer, parameter :: h_surf_r_polar_div_r = h_surf_r_equatorial_div_r + 1
      integer, parameter :: h_center_omega_div_omega_crit = h_surf_r_polar_div_r + 1

      integer, parameter :: h_total_angular_momentum = h_center_omega_div_omega_crit + 1
      integer, parameter :: h_log_total_angular_momentum = h_total_angular_momentum + 1

      integer, parameter :: h_min_t_eddy = h_log_total_angular_momentum + 1
      integer, parameter :: h_elapsed_time = h_min_t_eddy + 1
      
      integer, parameter :: h_num_hydro_merges = h_elapsed_time + 1
      integer, parameter :: h_num_hydro_splits = h_num_hydro_merges + 1
      
      integer, parameter :: h_RSP_DeltaR = h_num_hydro_splits + 1
      integer, parameter :: h_RSP_DeltaMag = h_RSP_DeltaR + 1
      integer, parameter :: h_RSP_GRPDV = h_RSP_DeltaMag + 1
      integer, parameter :: h_RSP_GREKM_avg_abs = h_RSP_GRPDV + 1
      integer, parameter :: h_RSP_GREKM = h_RSP_GREKM_avg_abs + 1
      
      integer, parameter :: h_rsp_phase = h_RSP_GREKM + 1
      integer, parameter :: h_rsp_period_in_days = h_rsp_phase + 1
      integer, parameter :: h_rsp_num_periods = h_rsp_period_in_days + 1
      
      integer, parameter :: h_RSP_LINA_period_F_days = h_rsp_num_periods + 1
      integer, parameter :: h_RSP_LINA_period_O1_days = h_RSP_LINA_period_F_days + 1
      integer, parameter :: h_RSP_LINA_period_O2_days = h_RSP_LINA_period_O1_days + 1
      integer, parameter :: h_RSP_LINA_growth_rate_F = h_RSP_LINA_period_O2_days + 1
      integer, parameter :: h_RSP_LINA_growth_rate_O1 = h_RSP_LINA_growth_rate_F + 1
      integer, parameter :: h_RSP_LINA_growth_rate_O2 = h_RSP_LINA_growth_rate_O1 + 1

      integer, parameter :: h_total_num_solver_iterations = h_RSP_LINA_growth_rate_O2 + 1
      integer, parameter :: h_total_num_solver_calls_made = h_total_num_solver_iterations + 1
      integer, parameter :: h_total_num_solver_calls_converged = h_total_num_solver_calls_made + 1
      integer, parameter :: h_total_num_solver_calls_failed = h_total_num_solver_calls_converged + 1
      
      integer, parameter :: h_total_num_solver_relax_iterations = h_total_num_solver_calls_failed + 1
      integer, parameter :: h_total_num_solver_relax_calls_made = h_total_num_solver_relax_iterations + 1
      integer, parameter :: h_total_num_solver_relax_calls_converged = h_total_num_solver_relax_calls_made + 1
      integer, parameter :: h_total_num_solver_relax_calls_failed = h_total_num_solver_relax_calls_converged + 1

      integer, parameter :: h_total_step_attempts = h_total_num_solver_relax_calls_failed + 1
      integer, parameter :: h_total_step_retries = h_total_step_attempts + 1
      integer, parameter :: h_total_step_redos = h_total_step_retries + 1
      integer, parameter :: h_total_steps_taken = h_total_step_redos + 1
      integer, parameter :: h_total_steps_finished = h_total_steps_taken + 1
      
      integer, parameter :: h_total_relax_step_attempts = h_total_steps_finished + 1
      integer, parameter :: h_total_relax_step_retries = h_total_relax_step_attempts + 1
      integer, parameter :: h_total_relax_step_redos = h_total_relax_step_retries + 1
      integer, parameter :: h_total_relax_steps_taken = h_total_relax_step_redos + 1
      integer, parameter :: h_total_relax_steps_finished = h_total_relax_steps_taken + 1
      
      integer, parameter :: h_avg_num_solver_iters = h_total_relax_steps_finished + 1
      integer, parameter :: h_num_solver_iterations = h_avg_num_solver_iters + 1
      integer, parameter :: h_num_iters = h_num_solver_iterations + 1

      integer, parameter :: h_photosphere_cell_density = h_num_iters + 1
      integer, parameter :: h_photosphere_cell_log_density = h_photosphere_cell_density + 1
      
      integer, parameter :: h_photosphere_cell_log_opacity = h_photosphere_cell_log_density + 1
      integer, parameter :: h_photosphere_cell_opacity = h_photosphere_cell_log_opacity + 1

      integer, parameter :: h_photosphere_cell_log_free_e = h_photosphere_cell_opacity + 1
      integer, parameter :: h_photosphere_cell_free_e = h_photosphere_cell_log_free_e + 1

      integer, parameter :: h_photosphere_cell_k = h_photosphere_cell_free_e + 1
      integer, parameter :: h_photosphere_cell_log_T = h_photosphere_cell_k + 1
      integer, parameter :: h_photosphere_cell_T = h_photosphere_cell_log_T + 1
      integer, parameter :: h_photosphere_black_body_T = h_photosphere_cell_T + 1
      integer, parameter :: h_photosphere_xm = h_photosphere_black_body_T + 1
      integer, parameter :: h_photosphere_logg = h_photosphere_xm + 1
      integer, parameter :: h_photosphere_T = h_photosphere_logg + 1
      integer, parameter :: h_photosphere_m = h_photosphere_T + 1
      integer, parameter :: h_photosphere_log_L = h_photosphere_m + 1
      integer, parameter :: h_photosphere_L = h_photosphere_log_L + 1
      integer, parameter :: h_log_one_div_yphot = h_photosphere_L + 1
      integer, parameter :: h_one_div_yphot = h_log_one_div_yphot + 1
      integer, parameter :: h_photosphere_column_density = h_one_div_yphot + 1
      integer, parameter :: h_photosphere_log_column_density = h_photosphere_column_density + 1
      integer, parameter :: h_photosphere_opacity = h_photosphere_log_column_density + 1
      integer, parameter :: h_photosphere_csound = h_photosphere_opacity + 1
      integer, parameter :: h_photosphere_v_div_cs = h_photosphere_csound + 1
      integer, parameter :: h_v_phot_km_s = h_photosphere_v_div_cs + 1
      integer, parameter :: h_photosphere_v_km_s = h_v_phot_km_s + 1
      integer, parameter :: h_photosphere_log_r = h_photosphere_v_km_s + 1
      integer, parameter :: h_photosphere_r = h_photosphere_log_r + 1

      integer, parameter :: h_min_opacity = h_photosphere_r + 1
      integer, parameter :: h_log_min_opacity = h_min_opacity + 1
      
      integer, parameter :: h_delta_nu = h_log_min_opacity + 1
      integer, parameter :: h_delta_Pg = h_delta_nu + 1
      integer, parameter :: h_nu_max = h_delta_Pg + 1
      integer, parameter :: h_acoustic_radius = h_nu_max + 1
      integer, parameter :: h_acoustic_cutoff = h_acoustic_radius + 1
      integer, parameter :: h_gs_per_delta_nu = h_acoustic_cutoff + 1
      integer, parameter :: h_ng_for_nu_max = h_gs_per_delta_nu + 1
      integer, parameter :: h_log_delta_Pg = h_ng_for_nu_max + 1
      integer, parameter :: h_nu_max_3_4th_div_delta_nu = h_log_delta_Pg + 1

      integer, parameter :: h_int_k_r_dr_nu_max_Sl1 = h_nu_max_3_4th_div_delta_nu + 1
      integer, parameter :: h_int_k_r_dr_2pt0_nu_max_Sl1 = h_int_k_r_dr_nu_max_Sl1 + 1
      integer, parameter :: h_int_k_r_dr_0pt5_nu_max_Sl1 = h_int_k_r_dr_2pt0_nu_max_Sl1 + 1

      integer, parameter :: h_int_k_r_dr_nu_max_Sl2 = h_int_k_r_dr_0pt5_nu_max_Sl1 + 1
      integer, parameter :: h_int_k_r_dr_2pt0_nu_max_Sl2 = h_int_k_r_dr_nu_max_Sl2 + 1
      integer, parameter :: h_int_k_r_dr_0pt5_nu_max_Sl2 = h_int_k_r_dr_2pt0_nu_max_Sl2 + 1

      integer, parameter :: h_int_k_r_dr_nu_max_Sl3 = h_int_k_r_dr_0pt5_nu_max_Sl2 + 1
      integer, parameter :: h_int_k_r_dr_2pt0_nu_max_Sl3 = h_int_k_r_dr_nu_max_Sl3 + 1
      integer, parameter :: h_int_k_r_dr_0pt5_nu_max_Sl3 = h_int_k_r_dr_2pt0_nu_max_Sl3 + 1

      integer, parameter :: h_log_Lnuc_sub_log_L = h_int_k_r_dr_0pt5_nu_max_Sl3 + 1
      integer, parameter :: h_cz_top_zone_logdq = h_log_Lnuc_sub_log_L + 1

      integer, parameter :: h_mass_semiconv_core = h_cz_top_zone_logdq + 1
      integer, parameter :: h_mass_conv_core = h_mass_semiconv_core + 1

      integer, parameter :: h_trace_mass_lgP = h_mass_conv_core + 1
      integer, parameter :: h_trace_mass_g = h_trace_mass_lgP + 1
      integer, parameter :: h_trace_mass_X = h_trace_mass_g + 1
      integer, parameter :: h_trace_mass_Y = h_trace_mass_X + 1
      integer, parameter :: h_trace_mass_edv_H = h_trace_mass_Y + 1
      integer, parameter :: h_trace_mass_edv_He = h_trace_mass_edv_H + 1
      integer, parameter :: h_trace_mass_scale_height = h_trace_mass_edv_He + 1
      integer, parameter :: h_trace_mass_dlnX_dr = h_trace_mass_scale_height + 1
      integer, parameter :: h_trace_mass_dlnY_dr = h_trace_mass_dlnX_dr + 1
      integer, parameter :: h_trace_mass_dlnRho_dr = h_trace_mass_dlnY_dr + 1

      integer, parameter :: h_k_below_const_q = h_trace_mass_dlnRho_dr + 1
      integer, parameter :: h_q_below_const_q = h_k_below_const_q + 1
      integer, parameter :: h_logxq_below_const_q = h_q_below_const_q + 1

      integer, parameter :: h_k_const_mass = h_logxq_below_const_q + 1
      integer, parameter :: h_q_const_mass = h_k_const_mass + 1
      integer, parameter :: h_logxq_const_mass = h_q_const_mass + 1

      integer, parameter :: h_k_below_just_added = h_logxq_const_mass + 1
      integer, parameter :: h_q_below_just_added = h_k_below_just_added + 1
      integer, parameter :: h_logxq_below_just_added = h_q_below_just_added + 1

      integer, parameter :: h_k_for_test_CpT_absMdot_div_L = h_logxq_below_just_added + 1
      integer, parameter :: h_q_for_test_CpT_absMdot_div_L = h_k_for_test_CpT_absMdot_div_L + 1
      integer, parameter :: h_logxq_for_test_CpT_absMdot_div_L = h_q_for_test_CpT_absMdot_div_L + 1

      integer, parameter :: h_k_CpTMdot_lt_L = h_logxq_for_test_CpT_absMdot_div_L + 1
      integer, parameter :: h_q_CpTMdot_lt_L = h_k_CpTMdot_lt_L + 1
      integer, parameter :: h_logxq_CpTMdot_lt_L = h_q_CpTMdot_lt_L + 1

      integer, parameter :: h_tot_E = h_logxq_CpTMdot_lt_L + 1
      integer, parameter :: h_log_tot_E = h_tot_E + 1
      integer, parameter :: h_tot_KE = h_log_tot_E + 1
      integer, parameter :: h_log_tot_KE = h_tot_KE + 1
      integer, parameter :: h_tot_IE = h_log_tot_KE + 1
      integer, parameter :: h_log_tot_IE = h_tot_IE + 1
      integer, parameter :: h_tot_PE = h_log_tot_IE + 1
      integer, parameter :: h_log_tot_PE = h_tot_PE + 1

      integer, parameter :: h_tot_Et = h_log_tot_PE + 1
      integer, parameter :: h_log_tot_Et = h_tot_Et + 1

      integer, parameter :: h_tot_IE_div_IE_plus_KE = h_log_tot_Et + 1      

      integer, parameter :: h_surface_extra_Pgas = h_tot_IE_div_IE_plus_KE + 1
      integer, parameter :: h_min_L = h_surface_extra_Pgas + 1
      integer, parameter :: h_min_dL_dm_m = h_min_L + 1
      integer, parameter :: h_min_dL_dm = h_min_dL_dm_m + 1
      integer, parameter :: h_burn_solver_maxsteps = h_min_dL_dm + 1
      integer, parameter :: h_rotation_solver_steps = h_burn_solver_maxsteps + 1
      integer, parameter :: h_diffusion_solver_steps = h_rotation_solver_steps + 1
      integer, parameter :: h_diffusion_solver_iters = h_diffusion_solver_steps + 1
      integer, parameter :: h_total_radiation = h_diffusion_solver_iters + 1
      integer, parameter :: h_total_energy_plus_total_radiation = h_total_radiation + 1
      integer, parameter :: h_grav_dark_L_polar = h_total_energy_plus_total_radiation + 1
      integer, parameter :: h_grav_dark_Teff_polar = h_grav_dark_L_polar + 1
      integer, parameter :: h_grav_dark_L_equatorial = h_grav_dark_Teff_polar + 1
      integer, parameter :: h_grav_dark_Teff_equatorial = h_grav_dark_L_equatorial + 1

      integer, parameter :: h_apsidal_constant_k2 = h_grav_dark_Teff_equatorial + 1

      integer, parameter :: h_lg_Lnuc = h_apsidal_constant_k2 + 1
      integer, parameter :: h_H_rich = h_lg_Lnuc + 1
      integer, parameter :: h_N_cntr = h_H_rich + 1
      integer, parameter :: h_lg_Lneu = h_N_cntr + 1
      integer, parameter :: h_He_core = h_lg_Lneu + 1
      integer, parameter :: h_O_cntr = h_He_core + 1
      integer, parameter :: h_lg_Lphoto = h_O_cntr + 1
      integer, parameter :: h_C_core = h_lg_Lphoto + 1
      integer, parameter :: h_O_core = h_C_core + 1
      integer, parameter :: h_Si_core = h_O_core + 1
      integer, parameter :: h_Fe_core = h_Si_core + 1
      integer, parameter :: h_Ne_cntr = h_Fe_core + 1
      integer, parameter :: h_Mass = h_Ne_cntr + 1
      integer, parameter :: h_H_cntr = h_Mass + 1
      integer, parameter :: h_Si_cntr = h_H_cntr + 1
      integer, parameter :: h_lg_Mdot = h_Si_cntr + 1
      integer, parameter :: h_He_cntr = h_lg_Mdot + 1
      integer, parameter :: h_eta_cntr = h_He_cntr + 1
      integer, parameter :: h_gam_cntr = h_eta_cntr + 1
      integer, parameter :: h_v_div_cs = h_gam_cntr + 1
      integer, parameter :: h_zones = h_v_div_cs + 1
      integer, parameter :: h_lg_Dsurf = h_zones + 1
      integer, parameter :: h_C_cntr = h_lg_Dsurf + 1
      integer, parameter :: h_retries = h_C_cntr + 1

      integer, parameter :: h_col_id_max = h_retries

      character (len=maxlen_history_column_name) :: history_column_name(h_col_id_max)
      type (integer_dict), pointer :: history_column_names_dict


      contains


      subroutine history_column_names_init(ierr)
         use utils_lib, only: integer_dict_define
         integer, intent(out) :: ierr

         integer :: i, cnt
         ierr = 0
         cnt = 0
         history_column_name(:) = ''

         history_column_name(h_model_number) = 'model_number'
         history_column_name(h_log_star_age) = 'log_star_age'
         history_column_name(h_star_age) = 'star_age'
         history_column_name(h_log_star_age_sec) = 'log_star_age_sec'
         history_column_name(h_star_age_sec) = 'star_age_sec'
         history_column_name(h_star_age_min) = 'star_age_min'
         history_column_name(h_star_age_hr) = 'star_age_hr'
         history_column_name(h_star_age_day) = 'star_age_day'
         history_column_name(h_day) = 'day'

         history_column_name(h_star_mass) = 'star_mass'
         history_column_name(h_log_star_mass) = 'log_star_mass'
         history_column_name(h_log_xmstar) = 'log_xmstar'
         history_column_name(h_delta_mass) = 'delta_mass'
         history_column_name(h_star_mdot) = 'star_mdot'
         history_column_name(h_log_abs_mdot) = 'log_abs_mdot'
         history_column_name(h_time_step) = 'time_step'
         history_column_name(h_time_step_sec) = 'time_step_sec'
         history_column_name(h_e_thermal) = 'e_thermal'
         history_column_name(h_num_zones) = 'num_zones'
         history_column_name(h_Tsurf_factor) = 'Tsurf_factor'
         history_column_name(h_log_tau_center) = 'log_tau_center'
         history_column_name(h_tau_factor) = 'tau_factor'
         history_column_name(h_tau_surface) = 'tau_surface'
         history_column_name(h_species) = 'species'

         history_column_name(h_m_center_gm) = 'm_center_gm'
         history_column_name(h_r_center_cm) = 'r_center_cm'
         history_column_name(h_r_center_km) = 'r_center_km'
         history_column_name(h_m_center) = 'm_center'
         history_column_name(h_r_center) = 'r_center'
         history_column_name(h_L_center) = 'L_center'
         history_column_name(h_log_L_center_ergs_s) = 'log_L_center_ergs_s'
         history_column_name(h_log_L_center) = 'log_L_center'
         history_column_name(h_infall_div_cs) = 'infall_div_cs'
         history_column_name(h_v_center) = 'v_center'
         history_column_name(h_v_center_kms) = 'v_center_kms'

         history_column_name(h_mdot_timescale) = 'mdot_timescale'
         history_column_name(h_kh_div_mdot_timescales) = 'kh_div_mdot_timescales'

         history_column_name(h_star_gravitational_mass) = 'star_gravitational_mass'
         history_column_name(h_star_mass_grav_div_mass) = 'star_mass_grav_div_mass'

         history_column_name(h_conv_mx1_top) = 'conv_mx1_top'
         history_column_name(h_conv_mx1_bot) = 'conv_mx1_bot'
         history_column_name(h_conv_mx2_top) = 'conv_mx2_top'
         history_column_name(h_conv_mx2_bot) = 'conv_mx2_bot'
         history_column_name(h_mx1_top) = 'mx1_top'
         history_column_name(h_mx1_bot) = 'mx1_bot'
         history_column_name(h_mx2_top) = 'mx2_top'
         history_column_name(h_mx2_bot) = 'mx2_bot'

         history_column_name(h_conv_mx1_top_r) = 'conv_mx1_top_r'
         history_column_name(h_conv_mx1_bot_r) = 'conv_mx1_bot_r'
         history_column_name(h_conv_mx2_top_r) = 'conv_mx2_top_r'
         history_column_name(h_conv_mx2_bot_r) = 'conv_mx2_bot_r'
         history_column_name(h_mx1_top_r) = 'mx1_top_r'
         history_column_name(h_mx1_bot_r) = 'mx1_bot_r'
         history_column_name(h_mx2_top_r) = 'mx2_top_r'
         history_column_name(h_mx2_bot_r) = 'mx2_bot_r'

         history_column_name(h_mix_relr_regions) = 'mix_relr_regions'
         history_column_name(h_mixing_regions) = 'mixing_regions'

         history_column_name(h_epsnuc_M_1) = 'epsnuc_M_1'
         history_column_name(h_epsnuc_M_2) = 'epsnuc_M_2'
         history_column_name(h_epsnuc_M_3) = 'epsnuc_M_3'
         history_column_name(h_epsnuc_M_4) = 'epsnuc_M_4'
         history_column_name(h_epsnuc_M_5) = 'epsnuc_M_5'
         history_column_name(h_epsnuc_M_6) = 'epsnuc_M_6'
         history_column_name(h_epsnuc_M_7) = 'epsnuc_M_7'
         history_column_name(h_epsnuc_M_8) = 'epsnuc_M_8'
         history_column_name(h_burning_regions) = 'burning_regions'
         history_column_name(h_burn_relr_regions) = 'burn_relr_regions'

         history_column_name(h_log_cntr_dr_cm) = 'log_cntr_dr_cm'
         history_column_name(h_log_max_T) = 'log_max_T'
         history_column_name(h_log_cntr_T) = 'log_cntr_T'
         history_column_name(h_log_center_T) = 'log_center_T'
         history_column_name(h_log_cntr_Rho) = 'log_cntr_Rho'
         history_column_name(h_log_center_Rho) = 'log_center_Rho'
         history_column_name(h_log_cntr_P) = 'log_cntr_P'
         history_column_name(h_log_center_P) = 'log_center_P'

         history_column_name(h_max_T) = 'max_T'
         history_column_name(h_center_T) = 'center_T'
         history_column_name(h_center_Rho) = 'center_Rho'
         history_column_name(h_center_P) = 'center_P'

         history_column_name(h_center_zbar) = 'center_zbar'
         history_column_name(h_center_abar) = 'center_abar'
         history_column_name(h_center_mu) = 'center_mu'
         history_column_name(h_center_ye) = 'center_ye'
         history_column_name(h_center_entropy) = 'center_entropy'
         history_column_name(h_max_entropy) = 'max_entropy'
         history_column_name(h_max_infall_speed) = 'max_infall_speed'
         history_column_name(h_fe_core_infall) = 'fe_core_infall'
         history_column_name(h_non_fe_core_infall) = 'non_fe_core_infall'
         history_column_name(h_non_fe_core_rebound) = 'non_fe_core_rebound'
         history_column_name(h_compactness) = 'compactness'
         history_column_name(h_compactness_parameter) = 'compactness_parameter'
         history_column_name(h_v_surf_div_escape_v) = 'v_surf_div_escape_v'
         history_column_name(h_v_surf_km_s) = 'v_surf_km_s'
         history_column_name(h_v_surf) = 'v_surf'
         history_column_name(h_v_surf_div_v_kh) = 'v_surf_div_v_kh'
         history_column_name(h_v_div_csound_max) = 'v_div_csound_max'
         history_column_name(h_v_div_csound_surf) = 'v_div_csound_surf'
         history_column_name(h_log_dt) = 'log_dt'
         history_column_name(h_log_dt_sec) = 'log_dt_sec'
         history_column_name(h_time_step_days) = 'time_step_days'
         history_column_name(h_log_dt_days) = 'log_dt_days'
         

         history_column_name(h_power_h_burn) = 'power_h_burn'
         history_column_name(h_log_LH) = 'log_LH'
         
         history_column_name(h_power_he_burn) = 'power_he_burn'
         history_column_name(h_log_LHe) = 'log_LHe'
         
         history_column_name(h_power_photo) = 'power_photo'
         history_column_name(h_Lnuc_photo) = 'Lnuc_photo'
         
         history_column_name(h_power_z_burn) = 'power_z_burn'
         history_column_name(h_log_LZ) = 'log_LZ'
         
         history_column_name(h_luminosity) = 'luminosity'
         history_column_name(h_Lsurf_m) = 'Lsurf_m'
         history_column_name(h_log_L) = 'log_L'
         history_column_name(h_luminosity_ergs_s) = 'luminosity_ergs_s'
         history_column_name(h_log_L_ergs_s) = 'log_L_ergs_s'

         history_column_name(h_log_mesh_adjust_IE_conservation) = 'log_mesh_adjust_IE_conservation'
         history_column_name(h_log_mesh_adjust_PE_conservation) = 'log_mesh_adjust_PE_conservation'
         history_column_name(h_log_mesh_adjust_KE_conservation) = 'log_mesh_adjust_KE_conservation'

         history_column_name(h_total_IE_plus_KE_start) = 'total_IE_plus_KE_start'
         history_column_name(h_log_total_IE_plus_KE) = 'log_total_IE_plus_KE'
         history_column_name(h_total_IE_plus_KE) = 'total_IE_plus_KE'

         history_column_name(h_avg_abs_v_div_cs) = 'avg_abs_v_div_cs'
         history_column_name(h_log_avg_abs_v_div_cs) = 'log_avg_abs_v_div_cs'
         history_column_name(h_max_abs_v_div_cs) = 'max_abs_v_div_cs'
         history_column_name(h_log_max_abs_v_div_cs) = 'log_max_abs_v_div_cs'

         history_column_name(h_avg_abs_v) = 'avg_abs_v'
         history_column_name(h_log_avg_abs_v) = 'log_avg_abs_v'
         history_column_name(h_max_abs_v) = 'max_abs_v'
         history_column_name(h_log_max_abs_v) = 'log_max_abs_v'

         history_column_name(h_total_energy_plus_L_surf) = 'total_energy_plus_L_surf'

         history_column_name(h_total_internal_energy_after_adjust_mass) = 'total_internal_energy_after_adjust_mass'
         history_column_name(h_total_gravitational_energy_after_adjust_mass) = &
            'total_gravitational_energy_after_adjust_mass'
         history_column_name(h_total_radial_kinetic_energy_after_adjust_mass) = &
            'total_radial_kinetic_energy_after_adjust_mass'
         history_column_name(h_total_rotational_kinetic_energy_after_adjust_mass) = &
            'total_rotational_kinetic_energy_after_adjust_mass'
         history_column_name(h_total_turbulent_energy_after_adjust_mass) = 'total_turbulent_energy_after_adjust_mass'
         history_column_name(h_total_energy_after_adjust_mass) = 'total_energy_after_adjust_mass'

         history_column_name(h_total_internal_energy) = 'total_internal_energy'
         history_column_name(h_total_gravitational_energy) = 'total_gravitational_energy'
         history_column_name(h_total_radial_kinetic_energy) = 'total_radial_kinetic_energy'
         history_column_name(h_total_rotational_kinetic_energy) = 'total_rotational_kinetic_energy'
         history_column_name(h_total_turbulent_energy) = 'total_turbulent_energy'
         history_column_name(h_total_energy_foe) = 'total_energy_foe'
         history_column_name(h_total_energy) = 'total_energy'

         history_column_name(h_log_total_internal_energy) = 'log_total_internal_energy'
         history_column_name(h_log_total_gravitational_energy) = 'log_total_gravitational_energy'
         history_column_name(h_log_total_radial_kinetic_energy) = 'log_total_radial_kinetic_energy'
         history_column_name(h_log_total_rotational_kinetic_energy) = 'log_total_rotational_kinetic_energy'
         history_column_name(h_log_total_turbulent_energy) = 'log_total_turbulent_energy'
         history_column_name(h_log_total_energy) = 'log_total_energy'

         history_column_name(h_rms_dvdt_div_v) = 'rms_dvdt_div_v'

         history_column_name(h_total_entropy) = 'total_entropy'
         history_column_name(h_total_IE_div_IE_plus_KE) = 'total_IE_div_IE_plus_KE'

         history_column_name(h_virial_thm_P_avg) = 'virial_thm_P_avg'
         history_column_name(h_virial_thm_rel_err) = 'virial_thm_rel_err'
         history_column_name(h_total_eps_grav) = 'total_eps_grav'
         history_column_name(h_work_outward_at_surface) = 'work_outward_at_surface'
         history_column_name(h_work_inward_at_center) = 'work_inward_at_center'
         history_column_name(h_total_nuclear_heating) = 'total_nuclear_heating'
         history_column_name(h_total_non_nuc_neu_cooling) = 'total_non_nuc_neu_cooling'
         history_column_name(h_total_irradiation_heating) = 'total_irradiation_heating'
         history_column_name(h_total_WD_sedimentation_heating) = 'total_WD_sedimentation_heating'
         history_column_name(h_total_extra_heating) = 'total_extra_heating'

         history_column_name(h_total_energy_sources_and_sinks) = 'total_energy_sources_and_sinks'
         history_column_name(h_error_in_energy_conservation) = 'error_in_energy_conservation'
         history_column_name(h_rel_error_in_energy_conservation) = 'rel_error_in_energy_conservation'
         history_column_name(h_log_rel_error_in_energy_conservation) = 'log_rel_error_in_energy_conservation'

         history_column_name(h_cumulative_energy_error) = 'cumulative_energy_error'
         history_column_name(h_rel_cumulative_energy_error) = 'rel_cumulative_energy_error'
         
         history_column_name(h_abs_rel_E_err) = 'abs_rel_E_err'
         history_column_name(h_log_rel_E_err) = 'log_rel_E_err'
         
         history_column_name(h_tot_E_equ_err) = 'tot_E_equ_err'
         history_column_name(h_tot_E_err) = 'tot_E_err'
         history_column_name(h_rel_E_err) = 'rel_E_err'
         
         history_column_name(h_rel_run_E_err) = 'rel_run_E_err'
         history_column_name(h_log_rel_run_E_err) = 'log_rel_run_E_err'
         history_column_name(h_log_rel_cumulative_energy_error) = 'log_rel_cumulative_energy_error'

         history_column_name(h_log_residual_norm) = 'log_residual_norm'
         history_column_name(h_log_max_residual) = 'log_max_residual'

         history_column_name(h_log_max_dvdt_residual) = 'log_max_dvdt_residual'
         history_column_name(h_log_max_drdt_residual) = 'log_max_drdt_residual'
         history_column_name(h_log_max_lnd_residual) = 'log_max_lnd_residual'
         history_column_name(h_log_max_dEdt_residual) = 'log_max_dEdt_residual'

         history_column_name(h_log_max_abs_E_residual) = 'log_max_abs_E_residual'
         history_column_name(h_log_avg_E_residual) = 'log_avg_E_residual'
         history_column_name(h_max_abs_E_residual) = 'max_abs_E_residual'
         history_column_name(h_avg_E_residual) = 'avg_E_residual'

         history_column_name(h_u_surf_km_s) = 'u_surf_km_s'
         history_column_name(h_u_surf) = 'u_surf'
         history_column_name(h_u_div_csound_max) = 'u_div_csound_max'
         history_column_name(h_u_div_csound_surf) = 'u_div_csound_surf'

         history_column_name(h_max_abs_v_residual) = 'max_abs_v_residual'
         history_column_name(h_log_max_abs_v_residual) = 'log_max_abs_v_residual'
         history_column_name(h_avg_v_residual) = 'avg_v_residual'
         history_column_name(h_log_avg_v_residual) = 'log_avg_v_residual'

         history_column_name(h_log_Lneu_nuc) = 'log_Lneu_nuc'
         history_column_name(h_log_Lneu_nonnuc) = 'log_Lneu_nonnuc'
         history_column_name(h_log_Lneu) = 'log_Lneu'
         history_column_name(h_log_R) = 'log_R'
         history_column_name(h_radius) = 'radius'
         history_column_name(h_log_R_cm) = 'log_R_cm'
         history_column_name(h_radius_cm) = 'radius_cm'
         history_column_name(h_Teff) = 'Teff'
         history_column_name(h_log_Teff) = 'log_Teff'
         history_column_name(h_effective_T) = 'effective_T'

         history_column_name(h_gravity) = 'gravity'
         history_column_name(h_log_g) = 'log_g'
         history_column_name(h_log_L_div_Ledd) = 'log_L_div_Ledd'
         history_column_name(h_lum_div_Ledd) = 'lum_div_Ledd'
         history_column_name(h_max_L_rad_div_Ledd_div_phi_Joss) = 'max_L_rad_div_Ledd_div_phi_Joss'
         history_column_name(h_max_L_rad_div_Ledd) = 'max_L_rad_div_Ledd'
         history_column_name(h_log_Ledd) = 'log_Ledd'

         history_column_name(h_gradT_excess_alpha) = 'gradT_excess_alpha'
         history_column_name(h_gradT_excess_min_beta) = 'gradT_excess_min_beta'
         history_column_name(h_gradT_excess_max_lambda) = 'gradT_excess_max_lambda'

         history_column_name(h_gamma1_min) = 'gamma1_min'
         history_column_name(h_logT_max) = 'logT_max'
         history_column_name(h_logQ_max) = 'logQ_max'
         history_column_name(h_logQ_min) = 'logQ_min'

         history_column_name(h_num_hydro_merges) = 'num_hydro_merges'
         history_column_name(h_num_hydro_splits) = 'num_hydro_splits'

         history_column_name(h_RSP_DeltaR) = 'rsp_DeltaR'
         history_column_name(h_RSP_DeltaMag) = 'rsp_DeltaMag'
         history_column_name(h_RSP_GRPDV) = 'rsp_GRPDV'
         history_column_name(h_RSP_GREKM) = 'rsp_GREKM'
         history_column_name(h_RSP_GREKM_avg_abs) = 'rsp_GREKM_avg_abs'
         
         history_column_name(h_rsp_phase) = 'rsp_phase'
         history_column_name(h_rsp_period_in_days) = 'rsp_period_in_days'
         history_column_name(h_rsp_num_periods) = 'rsp_num_periods'

         history_column_name(h_RSP_LINA_period_F_days) = 'rsp_LINA_period_F_days'
         history_column_name(h_RSP_LINA_period_O1_days) = 'rsp_LINA_period_O1_days'
         history_column_name(h_RSP_LINA_period_O2_days) = 'rsp_LINA_period_O2_days'
         history_column_name(h_RSP_LINA_growth_rate_F) = 'rsp_LINA_growth_rate_F'
         history_column_name(h_RSP_LINA_growth_rate_O1) = 'rsp_LINA_growth_rate_O1'
         history_column_name(h_RSP_LINA_growth_rate_O2) = 'rsp_LINA_growth_rate_O2'

         history_column_name(h_total_num_solver_iterations) = 'total_num_solver_iterations'
         history_column_name(h_total_num_solver_calls_made) = 'total_num_solver_calls_made'
         history_column_name(h_total_num_solver_calls_converged) = 'total_num_solver_calls_converged'
         history_column_name(h_total_num_solver_calls_failed) = 'total_num_solver_calls_failed'

         history_column_name(h_total_num_solver_relax_iterations) = 'total_num_solver_relax_iterations'
         history_column_name(h_total_num_solver_relax_calls_made) = 'total_num_solver_relax_calls_made'
         history_column_name(h_total_num_solver_relax_calls_converged) = 'total_num_solver_relax_calls_converged'
         history_column_name(h_total_num_solver_relax_calls_failed) = 'total_num_solver_relax_calls_failed'


         history_column_name(h_total_step_attempts) = 'total_step_attempts'
         history_column_name(h_total_step_retries) = 'total_step_retries'
         history_column_name(h_total_step_redos) = 'total_step_redos'
         history_column_name(h_total_steps_taken) = 'total_steps_taken'
         history_column_name(h_total_steps_finished) = 'total_steps_finished'
      
         history_column_name(h_total_relax_step_attempts) = 'total_relax_step_attempts'
         history_column_name(h_total_relax_step_retries) = 'total_relax_step_retries'
         history_column_name(h_total_relax_step_redos) = 'total_relax_step_redos'
         history_column_name(h_total_relax_steps_taken) = 'total_relax_steps_taken'
         history_column_name(h_total_relax_steps_finished) = 'total_relax_steps_finished'

         history_column_name(h_avg_num_solver_iters) = 'avg_num_solver_iters'
         history_column_name(h_num_solver_iterations) = 'num_solver_iterations'
         history_column_name(h_num_iters) = 'num_iters'

         history_column_name(h_avg_skipped_setvars_per_step) = 'avg_skipped_setvars_per_step'
         history_column_name(h_avg_setvars_per_step) = 'avg_setvars_per_step'
         history_column_name(h_avg_solver_setvars_per_step) = 'avg_solver_setvars_per_step'
         
         history_column_name(h_num_steps_skipped_1st_setvars) = 'num_steps_skipped_1st_setvars'
         history_column_name(h_fraction_of_steps_skipped_1st_setvars) = 'fraction_of_steps_skipped_1st_setvars'
         history_column_name(h_num_retries) = 'num_retries'

         history_column_name(h_total_num_solver_iterations) = 'total_num_solver_iterations'
         
         
         
         history_column_name(h_h1_czb_mass) = 'h1_czb_mass'
         history_column_name(h_surf_c12_minus_o16) = 'surf_c12_minus_o16'
         history_column_name(h_surf_num_c12_div_num_o16) = 'surf_num_c12_div_num_o16'

         history_column_name(h_min_Pgas_div_P) = 'min_Pgas_div_P'
         history_column_name(h_center_degeneracy) = 'center_degeneracy'
         history_column_name(h_center_eps_grav) = 'center_eps_grav'
         history_column_name(h_center_non_nuc_neu) = 'center_non_nuc_neu'
         history_column_name(h_center_dL_dm) = 'center_dL_dm'
         history_column_name(h_log_center_eps_nuc) = 'log_center_eps_nuc'
         history_column_name(h_d_center_eps_nuc_dlnT) = 'd_center_eps_nuc_dlnT'
         history_column_name(h_d_center_eps_nuc_dlnd) = 'd_center_eps_nuc_dlnd'
         history_column_name(h_center_eps_nuc) = 'center_eps_nuc'
         history_column_name(h_center_gamma) = 'center_gamma'

         history_column_name(h_center_dlogT) = 'center_dlogT'
         history_column_name(h_center_dlogRho) = 'center_dlogRho'

         history_column_name(h_center_dlnT_dt) = 'center_dlnT_dt'
         history_column_name(h_center_dlnd_dt) = 'center_dlnd_dt'

         history_column_name(h_h_rich_layer_mass) = 'h_rich_layer_mass'
         history_column_name(h_he_rich_layer_mass) = 'he_rich_layer_mass'
         history_column_name(h_c_rich_layer_mass) = 'c_rich_layer_mass'
         history_column_name(h_o_rich_layer_mass) = 'o_rich_layer_mass'
         history_column_name(h_si_rich_layer_mass) = 'si_rich_layer_mass'

         history_column_name(h_he_core_mass) = 'he_core_mass'
         history_column_name(h_he_core_radius) = 'he_core_radius'
         history_column_name(h_he_core_lgT) = 'he_core_lgT'
         history_column_name(h_he_core_lgRho) = 'he_core_lgRho'
         history_column_name(h_he_core_L) = 'he_core_L'
         history_column_name(h_he_core_v) = 'he_core_v'
         history_column_name(h_he_core_omega) = 'he_core_omega'
         history_column_name(h_he_core_omega_div_omega_crit) = 'he_core_omega_div_omega_crit'
         history_column_name(h_he_core_k) = 'he_core_k'

         history_column_name(h_c_core_mass) = 'c_core_mass'
         history_column_name(h_c_core_radius) = 'c_core_radius'
         history_column_name(h_c_core_lgT) = 'c_core_lgT'
         history_column_name(h_c_core_lgRho) = 'c_core_lgRho'
         history_column_name(h_c_core_L) = 'c_core_L'
         history_column_name(h_c_core_v) = 'c_core_v'
         history_column_name(h_c_core_omega) = 'c_core_omega'
         history_column_name(h_c_core_omega_div_omega_crit) = 'c_core_omega_div_omega_crit'
         history_column_name(h_c_core_k) = 'c_core_k'

         history_column_name(h_o_core_mass) = 'o_core_mass'
         history_column_name(h_o_core_radius) = 'o_core_radius'
         history_column_name(h_o_core_lgT) = 'o_core_lgT'
         history_column_name(h_o_core_lgRho) = 'o_core_lgRho'
         history_column_name(h_o_core_L) = 'o_core_L'
         history_column_name(h_o_core_v) = 'o_core_v'
         history_column_name(h_o_core_omega) = 'o_core_omega'
         history_column_name(h_o_core_omega_div_omega_crit) = 'o_core_omega_div_omega_crit'
         history_column_name(h_o_core_k) = 'o_core_k'

         history_column_name(h_si_core_mass) = 'si_core_mass'
         history_column_name(h_si_core_radius) = 'si_core_radius'
         history_column_name(h_si_core_lgT) = 'si_core_lgT'
         history_column_name(h_si_core_lgRho) = 'si_core_lgRho'
         history_column_name(h_si_core_L) = 'si_core_L'
         history_column_name(h_si_core_v) = 'si_core_v'
         history_column_name(h_si_core_omega) = 'si_core_omega'
         history_column_name(h_si_core_omega_div_omega_crit) = 'si_core_omega_div_omega_crit'
         history_column_name(h_si_core_k) = 'si_core_k'

         history_column_name(h_fe_core_mass) = 'fe_core_mass'
         history_column_name(h_fe_core_radius) = 'fe_core_radius'
         history_column_name(h_fe_core_lgT) = 'fe_core_lgT'
         history_column_name(h_fe_core_lgRho) = 'fe_core_lgRho'
         history_column_name(h_fe_core_L) = 'fe_core_L'
         history_column_name(h_fe_core_v) = 'fe_core_v'
         history_column_name(h_fe_core_omega) = 'fe_core_omega'
         history_column_name(h_fe_core_omega_div_omega_crit) = 'fe_core_omega_div_omega_crit'
         history_column_name(h_fe_core_k) = 'fe_core_k'

         history_column_name(h_neutron_rich_core_mass) = 'neutron_rich_core_mass'
         history_column_name(h_neutron_rich_core_radius) = 'neutron_rich_core_radius'
         history_column_name(h_neutron_rich_core_lgT) = 'neutron_rich_core_lgT'
         history_column_name(h_neutron_rich_core_lgRho) = 'neutron_rich_core_lgRho'
         history_column_name(h_neutron_rich_core_L) = 'neutron_rich_core_L'
         history_column_name(h_neutron_rich_core_v) = 'neutron_rich_core_v'
         history_column_name(h_neutron_rich_core_omega) = 'neutron_rich_core_omega'
         history_column_name(h_neutron_rich_core_omega_div_omega_crit) = 'neutron_rich_core_omega_div_omega_crit'
         history_column_name(h_neutron_rich_core_k) = 'neutron_rich_core_k'

         history_column_name(h_envelope_mass) = 'envelope_mass'
         history_column_name(h_envelope_fraction_left) = 'envelope_fraction_left'

         history_column_name(h_tau10_mass) = &
            'tau10_mass' ! mass in solar units where optical depth = 10
         history_column_name(h_tau10_radius) = &
            'tau10_radius' ! radius in solar units where optical depth = 10
         history_column_name(h_tau10_lgP) = 'tau10_lgP' ! estimate for log10(P) at tau = 10
         history_column_name(h_tau10_T) = 'tau10_T' ! estimate for T at tau = 10
         history_column_name(h_tau10_lgT) = 'tau10_lgT' ! estimate for log10(T) at tau = 10
         history_column_name(h_tau10_lgRho) = 'tau10_lgRho' ! estimate for log10(density) at tau = 10
         history_column_name(h_tau10_L) = 'tau10_L' ! estimate for L/Lsun at tau = 10
         history_column_name(h_tau100_mass) = &
            'tau100_mass' ! location in solar units where optical depth = 100
         history_column_name(h_tau100_radius) = &
            'tau100_radius' ! location in solar units where optical depth = 100
         history_column_name(h_tau100_lgP) = 'tau100_lgP' ! estimates for values at tau = 100
         history_column_name(h_tau100_T) = 'tau100_T'
         history_column_name(h_tau100_lgT) = 'tau100_lgT'
         history_column_name(h_tau100_lgRho) = 'tau100_lgRho'
         history_column_name(h_tau100_L) = 'tau100_L'
         history_column_name(h_dynamic_timescale) = 'dynamic_timescale'
         history_column_name(h_dlnR_dlnM) = 'dlnR_dlnM'
         history_column_name(h_kh_timescale) = 'kh_timescale'
         history_column_name(h_nuc_timescale) = 'nuc_timescale'
         history_column_name(h_log_abs_Lgrav) = 'log_abs_Lgrav'
         history_column_name(h_eps_grav_integral) = 'eps_grav_integral'
         history_column_name(h_extra_L) = 'extra_L'
         history_column_name(h_log_extra_L) = 'log_extra_L'
         history_column_name(h_log_power_nuc_burn) = 'log_power_nuc_burn'
         history_column_name(h_log_Lnuc_ergs_s) = 'log_Lnuc_ergs_s'
         history_column_name(h_log_Lnuc) = 'log_Lnuc'
         history_column_name(h_log_Lneu) = 'log_Lneu'
         history_column_name(h_mass_loc_of_max_eps_nuc) = 'mass_loc_of_max_eps_nuc'
         history_column_name(h_mass_ext_to_max_eps_nuc) = 'mass_ext_to_max_eps_nuc'

         history_column_name(h_diffusion_time_H_He_bdy) = 'diffusion_time_H_He_bdy'
         history_column_name(h_temperature_H_He_bdy) = 'temperature_H_He_bdy'

         history_column_name(h_max_abs_v_velocity) = 'max_abs_v_velocity'
         history_column_name(h_max_abs_v_csound)  = 'max_abs_v_csound'
         history_column_name(h_max_abs_v_v_div_cs)  = 'max_abs_v_v_div_cs'
         history_column_name(h_max_abs_v_lgT)  = 'max_abs_v_lgT'
         history_column_name(h_max_abs_v_lgRho)  = 'max_abs_v_lgRho'
         history_column_name(h_max_abs_v_lgP)  = 'max_abs_v_lgP'
         history_column_name(h_max_abs_v_mass)  = 'max_abs_v_mass'
         history_column_name(h_max_abs_v_L)  = 'max_abs_v_L'
         history_column_name(h_max_abs_v_gamma1)  = 'max_abs_v_gamma1'
         history_column_name(h_max_abs_v_entropy)  = 'max_abs_v_entropy'
         history_column_name(h_max_abs_v_eps_nuc)  = 'max_abs_v_eps_nuc'
         history_column_name(h_max_abs_v_E0)  = 'max_abs_v_E0'

         history_column_name(h_total_ni_co_56)  = 'total_ni_co_56'

         history_column_name(h_max_abs_v_radius)  = 'max_abs_v_radius'
         history_column_name(h_max_abs_v_radius_cm)  = 'max_abs_v_radius_cm'
         history_column_name(h_max_abs_v_lgR)  = 'max_abs_v_lgR'
         history_column_name(h_max_abs_v_lgR_cm)  = 'max_abs_v_lgR_cm'

         history_column_name(h_inner_mach1_velocity) = 'inner_mach1_velocity'
         history_column_name(h_inner_mach1_csound)  = 'inner_mach1_csound'
         history_column_name(h_inner_mach1_v_div_cs)  = 'inner_mach1_v_div_cs'
         history_column_name(h_inner_mach1_lgT)  = 'inner_mach1_lgT'
         history_column_name(h_inner_mach1_lgRho)  = 'inner_mach1_lgRho'
         history_column_name(h_inner_mach1_lgP)  = 'inner_mach1_lgP'
         history_column_name(h_inner_mach1_q)  = 'inner_mach1_q'
         history_column_name(h_inner_mach1_tau)  = 'inner_mach1_tau'
         history_column_name(h_inner_mach1_mass)  = 'inner_mach1_mass'
         history_column_name(h_inner_mach1_radius)  = 'inner_mach1_radius'
         history_column_name(h_inner_mach1_gamma1)  = 'inner_mach1_gamma1'
         history_column_name(h_inner_mach1_entropy)  = 'inner_mach1_entropy'
         history_column_name(h_inner_mach1_k)  = 'inner_mach1_k'

         history_column_name(h_outer_mach1_velocity) = 'outer_mach1_velocity'
         history_column_name(h_outer_mach1_csound)  = 'outer_mach1_csound'
         history_column_name(h_outer_mach1_v_div_cs)  = 'outer_mach1_v_div_cs'
         history_column_name(h_outer_mach1_lgT)  = 'outer_mach1_lgT'
         history_column_name(h_outer_mach1_lgRho)  = 'outer_mach1_lgRho'
         history_column_name(h_outer_mach1_lgP)  = 'outer_mach1_lgP'
         history_column_name(h_outer_mach1_mass)  = 'outer_mach1_mass'
         history_column_name(h_outer_mach1_q)  = 'outer_mach1_q'
         history_column_name(h_outer_mach1_tau)  = 'outer_mach1_tau'
         history_column_name(h_outer_mach1_radius)  = 'outer_mach1_radius'
         history_column_name(h_outer_mach1_gamma1)  = 'outer_mach1_gamma1'
         history_column_name(h_outer_mach1_entropy)  = 'outer_mach1_entropy'
         history_column_name(h_outer_mach1_k)  = 'outer_mach1_k'

         history_column_name(h_shock_velocity) = 'shock_velocity'
         history_column_name(h_shock_csound)  = 'shock_csound'
         history_column_name(h_shock_v_div_cs)  = 'shock_v_div_cs'
         history_column_name(h_shock_lgT)  = 'shock_lgT'
         history_column_name(h_shock_lgRho)  = 'shock_lgRho'
         history_column_name(h_shock_lgP)  = 'shock_lgP'
         history_column_name(h_shock_q)  = 'shock_q'
         history_column_name(h_shock_tau)  = 'shock_tau'
         history_column_name(h_shock_mass)  = 'shock_mass'
         history_column_name(h_shock_radius)  = 'shock_radius'
         history_column_name(h_shock_mass_gm)  = 'shock_mass_gm'
         history_column_name(h_shock_radius_cm)  = 'shock_radius_cm'
         history_column_name(h_shock_gamma1)  = 'shock_gamma1'
         history_column_name(h_shock_entropy)  = 'shock_entropy'
         history_column_name(h_shock_pre_lgRho)  = 'shock_pre_lgRho'
         history_column_name(h_shock_k)  = 'shock_k'

         history_column_name(h_trace_mass_location) = 'trace_mass_location'
         history_column_name(h_trace_mass_radius) = 'trace_mass_radius'
         history_column_name(h_trace_mass_lgT) = 'trace_mass_lgT'
         history_column_name(h_trace_mass_lgRho) = 'trace_mass_lgRho'
         history_column_name(h_trace_mass_L) = 'trace_mass_L'
         history_column_name(h_trace_mass_v) = 'trace_mass_v'
         history_column_name(h_trace_mass_lgP) = 'trace_mass_lgP'
         history_column_name(h_trace_mass_g) = 'trace_mass_g'
         history_column_name(h_trace_mass_X) = 'trace_mass_X'
         history_column_name(h_trace_mass_Y) = 'trace_mass_Y'
         history_column_name(h_trace_mass_edv_H) = 'trace_mass_edv_H'
         history_column_name(h_trace_mass_edv_He) = 'trace_mass_edv_He'
         history_column_name(h_trace_mass_scale_height) = 'trace_mass_scale_height'
         history_column_name(h_trace_mass_dlnX_dr) = 'trace_mass_dlnX_dr'
         history_column_name(h_trace_mass_dlnY_dr) = 'trace_mass_dlnY_dr'
         history_column_name(h_trace_mass_dlnRho_dr) = 'trace_mass_dlnRho_dr'

         history_column_name(h_max_T_shell_binding_energy) = 'max_T_shell_binding_energy'
         history_column_name(h_max_T_lgP_thin_shell) = 'max_T_lgP_thin_shell'
         history_column_name(h_max_T_lgT) = 'max_T_lgT'
         history_column_name(h_max_T_lgP) = 'max_T_lgP'
         history_column_name(h_max_T_entropy) = 'max_T_entropy'
         history_column_name(h_max_T_mass) = 'max_T_mass'
         history_column_name(h_max_T_radius) = 'max_T_radius'
         history_column_name(h_max_T_lgRho) = 'max_T_lgRho'
         history_column_name(h_max_T_L) = 'max_T_L'
         history_column_name(h_max_T_eps_nuc) = 'max_T_eps_nuc'

         history_column_name(h_log_surf_optical_depth) = 'log_surf_optical_depth'
         history_column_name(h_surface_optical_depth) = 'surface_optical_depth'

         history_column_name(h_log_surf_cell_opacity) = 'log_surf_cell_opacity'
         history_column_name(h_log_surf_cell_density) = 'log_surf_cell_density'
         history_column_name(h_surface_cell_temperature) = 'surface_cell_temperature'
         history_column_name(h_log_surf_cell_temperature) = 'log_surf_cell_temperature'
         history_column_name(h_log_surf_cell_pressure) = 'log_surf_cell_pressure'
         history_column_name(h_log_surf_cell_P) = 'log_surf_cell_P'
         history_column_name(h_log_surf_cell_z) = 'log_surf_cell_z'
         history_column_name(h_surface_cell_entropy) = 'surface_cell_entropy'

         history_column_name(h_max_conv_vel_div_csound) = 'max_conv_vel_div_csound'
         history_column_name(h_max_gradT_div_grada) = 'max_gradT_div_grada'
         history_column_name(h_max_gradT_sub_grada) = 'max_gradT_sub_grada'
         history_column_name(h_min_log_mlt_Gamma) = 'min_log_mlt_Gamma'

         history_column_name(h_dt_cell_collapse) = 'dt_cell_collapse'
         history_column_name(h_dt_div_dt_cell_collapse) = 'dt_div_dt_cell_collapse'

         history_column_name(h_min_dr_div_cs_k) = 'min_dr_div_cs_k'
         history_column_name(h_min_dr_div_cs) = 'min_dr_div_cs'
         history_column_name(h_log_min_dr_div_cs) = 'log_min_dr_div_cs'
         history_column_name(h_min_dr_div_cs_yr) = 'min_dr_div_cs_yr'
         history_column_name(h_log_min_dr_div_cs_yr) = 'log_min_dr_div_cs_yr'
         history_column_name(h_dt_div_min_dr_div_cs) = 'dt_div_min_dr_div_cs'
         history_column_name(h_log_dt_div_min_dr_div_cs) = 'log_dt_div_min_dr_div_cs'

         history_column_name(h_cz_bot_mass) = 'cz_bot_mass'
         history_column_name(h_cz_mass) = 'cz_mass'
         history_column_name(h_cz_log_xmsun) = 'cz_log_xmsun'
         history_column_name(h_cz_log_xmass) = 'cz_log_xmass'
         history_column_name(h_cz_xm) = 'cz_xm'
         history_column_name(h_cz_logT) = 'cz_logT'
         history_column_name(h_cz_logRho) = 'cz_logRho'
         history_column_name(h_cz_logP) = 'cz_logP'
         history_column_name(h_cz_log_column_depth) = 'cz_log_column_depth'
         history_column_name(h_cz_log_radial_depth) = 'cz_log_radial_depth'

         history_column_name(h_cz_log_tau) = 'cz_log_tau'
         history_column_name(h_cz_opacity) = 'cz_opacity'
         history_column_name(h_cz_eta) = 'cz_eta'
         history_column_name(h_cz_log_eps_nuc) = 'cz_log_eps_nuc'
         history_column_name(h_cz_t_heat) = 'cz_t_heat'
         history_column_name(h_cz_csound) = 'cz_csound'
         history_column_name(h_cz_scale_height) = 'cz_scale_height'
         history_column_name(h_cz_grav) = 'cz_grav'
         history_column_name(h_cz_luminosity) = 'cz_luminosity'
         history_column_name(h_cz_bot_radius) = 'cz_bot_radius'
         history_column_name(h_cz_zone) = 'cz_zone'
         history_column_name(h_cz_omega) = 'cz_omega'
         history_column_name(h_cz_omega_div_omega_crit) = 'cz_omega_div_omega_crit'

         history_column_name(h_cz_top_mass) = 'cz_top_mass'
         history_column_name(h_cz_top_log_xmsun) = 'cz_top_log_xmsun'
         history_column_name(h_cz_top_log_xmass) = 'cz_top_log_xmass'
         history_column_name(h_cz_top_xm) = 'cz_top_xm'
         history_column_name(h_cz_top_logT) = 'cz_top_logT'
         history_column_name(h_cz_top_logRho) = 'cz_top_logRho'
         history_column_name(h_cz_top_logP) = 'cz_top_logP'
         history_column_name(h_cz_top_log_column_depth) = 'cz_top_log_column_depth'
         history_column_name(h_cz_top_log_radial_depth) = 'cz_top_log_radial_depth'

         history_column_name(h_cz_top_log_tau) = 'cz_top_log_tau'
         history_column_name(h_cz_top_opacity) = 'cz_top_opacity'
         history_column_name(h_cz_top_eta) = 'cz_top_eta'
         history_column_name(h_cz_top_log_eps_nuc) = 'cz_top_log_eps_nuc'
         history_column_name(h_cz_top_t_heat) = 'cz_top_t_heat'
         history_column_name(h_cz_top_csound) = 'cz_top_csound'
         history_column_name(h_cz_top_scale_height) = 'cz_top_scale_height'
         history_column_name(h_cz_top_grav) = 'cz_top_grav'
         history_column_name(h_cz_top_luminosity) = 'cz_top_luminosity'
         history_column_name(h_cz_top_radius) = 'cz_top_radius'
         history_column_name(h_cz_top_zone) = 'cz_top_zone'
         history_column_name(h_cz_top_omega) = 'cz_top_omega'
         history_column_name(h_cz_top_omega_div_omega_crit) = 'cz_top_omega_div_omega_crit'

         history_column_name(h_trace_mass_omega) = 'trace_mass_omega'
         history_column_name(h_trace_mass_omega_div_omega_crit) = 'trace_mass_omega_div_omega_crit'
         history_column_name(h_kh_mdot_limit) = 'kh_mdot_limit'
         history_column_name(h_rotational_mdot_boost) = 'rotational_mdot_boost'
         history_column_name(h_log_rotational_mdot_boost) = 'log_rotational_mdot_boost'

         history_column_name(h_logL_for_BB_outer_BC) = 'logL_for_BB_outer_BC'
         history_column_name(h_luminosity_for_BB_outer_BC) = 'luminosity_for_BB_outer_BC'

         history_column_name(h_i_rot_total) = 'i_rot_total'
         history_column_name(h_surf_avg_j_rot) = 'surf_avg_j_rot'
         history_column_name(h_surf_avg_omega) = 'surf_avg_omega'
         history_column_name(h_surf_avg_omega_crit) = 'surf_avg_omega_crit'
         history_column_name(h_surf_avg_omega_div_omega_crit) = 'surf_avg_omega_div_omega_crit'

         history_column_name(h_surf_avg_v_rot) = 'surf_avg_v_rot'
         history_column_name(h_surf_avg_v_crit) = 'surf_avg_v_crit'
         history_column_name(h_surf_avg_v_div_v_crit) = 'surf_avg_v_div_v_crit'

         history_column_name(h_surf_avg_logT) = 'surf_avg_logT'
         history_column_name(h_surf_avg_logRho) = 'surf_avg_logRho'
         history_column_name(h_surf_avg_opacity) = 'surf_avg_opacity'
         history_column_name(h_surf_avg_Lrad_div_Ledd) = 'surf_avg_Lrad_div_Ledd'

         history_column_name(h_surf_escape_v) = 'surf_escape_v'
         history_column_name(h_v_wind_Km_per_s) = 'v_wind_Km_per_s'

         history_column_name(h_center_omega) = 'center_omega'
         history_column_name(h_center_omega_div_omega_crit) = 'center_omega_div_omega_crit'

         history_column_name(h_surf_r_equatorial_div_r_polar) = 'surf_r_equatorial_div_r_polar'
         history_column_name(h_surf_r_equatorial_div_r) = 'surf_r_equatorial_div_r'
         history_column_name(h_surf_r_polar_div_r) = 'surf_r_polar_div_r'

         history_column_name(h_total_angular_momentum) = 'total_angular_momentum'
         history_column_name(h_log_total_angular_momentum) = 'log_total_angular_momentum'

         history_column_name(h_min_t_eddy) = 'min_t_eddy'
         history_column_name(h_elapsed_time) = 'elapsed_time'

         history_column_name(h_photosphere_cell_log_density) = 'photosphere_cell_log_density'
         history_column_name(h_photosphere_cell_density) = 'photosphere_cell_density'
         
         history_column_name(h_photosphere_cell_log_opacity) = 'photosphere_cell_log_opacity'
         history_column_name(h_photosphere_cell_opacity) = 'photosphere_cell_opacity'

         history_column_name(h_photosphere_cell_log_free_e) = 'photosphere_cell_log_free_e'
         history_column_name(h_photosphere_cell_free_e) = 'photosphere_cell_free_e'

         history_column_name(h_photosphere_cell_k) = 'photosphere_cell_k'
         history_column_name(h_photosphere_cell_log_T) = 'photosphere_cell_log_T'
         history_column_name(h_photosphere_cell_T) = 'photosphere_cell_T'
         history_column_name(h_photosphere_black_body_T) = 'photosphere_black_body_T'
         history_column_name(h_photosphere_T) = 'photosphere_T'
         history_column_name(h_photosphere_logg) = 'photosphere_logg'
         history_column_name(h_photosphere_m) = 'photosphere_m'
         history_column_name(h_photosphere_xm) = 'photosphere_xm'
         history_column_name(h_photosphere_L) = 'photosphere_L'
         history_column_name(h_photosphere_r) = 'photosphere_r'
         history_column_name(h_photosphere_log_L) = 'photosphere_log_L'
         history_column_name(h_photosphere_log_r) = 'photosphere_log_r'
         history_column_name(h_photosphere_v_km_s) = 'photosphere_v_km_s'
         history_column_name(h_v_phot_km_s) = 'v_phot_km_s'
         history_column_name(h_one_div_yphot) = 'one_div_yphot'
         history_column_name(h_log_one_div_yphot) = 'log_one_div_yphot'
         history_column_name(h_photosphere_column_density) = 'photosphere_column_density'
         history_column_name(h_photosphere_log_column_density) = 'photosphere_log_column_density'
         history_column_name(h_photosphere_opacity) = 'photosphere_opacity'
         history_column_name(h_photosphere_csound) = 'photosphere_csound'
         history_column_name(h_photosphere_v_div_cs) = 'photosphere_v_div_cs'

         history_column_name(h_min_opacity) = 'min_opacity'
         history_column_name(h_log_min_opacity) = 'log_min_opacity'

         history_column_name(h_delta_nu) = 'delta_nu'
         history_column_name(h_delta_Pg) = 'delta_Pg'
         history_column_name(h_log_delta_Pg) = 'log_delta_Pg'
         history_column_name(h_nu_max_3_4th_div_delta_nu) = 'nu_max_3_4th_div_delta_nu'

         history_column_name(h_nu_max) = 'nu_max'
         history_column_name(h_acoustic_cutoff) = 'acoustic_cutoff'
         history_column_name(h_acoustic_radius) = 'acoustic_radius'
         history_column_name(h_gs_per_delta_nu) = 'gs_per_delta_nu'
         history_column_name(h_ng_for_nu_max) = 'ng_for_nu_max'

         history_column_name(h_int_k_r_dr_nu_max_Sl1) = 'int_k_r_dr_nu_max_Sl1'
         history_column_name(h_int_k_r_dr_2pt0_nu_max_Sl1) = 'int_k_r_dr_2pt0_nu_max_Sl1'
         history_column_name(h_int_k_r_dr_0pt5_nu_max_Sl1) = 'int_k_r_dr_0pt5_nu_max_Sl1'

         history_column_name(h_int_k_r_dr_nu_max_Sl2) = 'int_k_r_dr_nu_max_Sl2'
         history_column_name(h_int_k_r_dr_2pt0_nu_max_Sl2) = 'int_k_r_dr_2pt0_nu_max_Sl2'
         history_column_name(h_int_k_r_dr_0pt5_nu_max_Sl2) = 'int_k_r_dr_0pt5_nu_max_Sl2'

         history_column_name(h_int_k_r_dr_nu_max_Sl3) = 'int_k_r_dr_nu_max_Sl3'
         history_column_name(h_int_k_r_dr_2pt0_nu_max_Sl3) = 'int_k_r_dr_2pt0_nu_max_Sl3'
         history_column_name(h_int_k_r_dr_0pt5_nu_max_Sl3) = 'int_k_r_dr_0pt5_nu_max_Sl3'

         history_column_name(h_mass_conv_core) = 'mass_conv_core'
         history_column_name(h_mass_semiconv_core) = 'mass_semiconv_core'

         history_column_name(h_cz_top_zone_logdq) = 'cz_top_zone_logdq'
         history_column_name(h_log_Lnuc_sub_log_L) = 'log_Lnuc_sub_log_L'

         history_column_name(h_surface_extra_Pgas) = 'surface_extra_Pgas'

         history_column_name(h_min_L) = 'min_L'
         history_column_name(h_min_dL_dm_m) = 'min_dL_dm_m'
         history_column_name(h_min_dL_dm) = 'min_dL_dm'

         history_column_name(h_k_below_const_q) = 'k_below_const_q'
         history_column_name(h_q_below_const_q) = 'q_below_const_q'
         history_column_name(h_logxq_below_const_q) = 'logxq_below_const_q'

         history_column_name(h_k_const_mass) = 'k_const_mass'
         history_column_name(h_q_const_mass) = 'q_const_mass'
         history_column_name(h_logxq_const_mass) = 'logxq_const_mass'

         history_column_name(h_k_below_just_added) = 'k_below_just_added'
         history_column_name(h_q_below_just_added) = 'q_below_just_added'
         history_column_name(h_logxq_below_just_added) = 'logxq_below_just_added'

         history_column_name(h_k_for_test_CpT_absMdot_div_L) = 'k_for_test_CpT_absMdot_div_L'
         history_column_name(h_q_for_test_CpT_absMdot_div_L) = 'q_for_test_CpT_absMdot_div_L'
         history_column_name(h_logxq_for_test_CpT_absMdot_div_L) = 'logxq_for_test_CpT_absMdot_div_L'

         history_column_name(h_k_CpTMdot_lt_L) = 'k_CpTMdot_lt_L'
         history_column_name(h_q_CpTMdot_lt_L) = 'q_CpTMdot_lt_L'
         history_column_name(h_logxq_CpTMdot_lt_L) = 'logxq_CpTMdot_lt_L'

         history_column_name(h_burn_solver_maxsteps) = 'burn_solver_maxsteps'
         history_column_name(h_rotation_solver_steps) = 'rotation_solver_steps'
         history_column_name(h_diffusion_solver_steps) = 'diffusion_solver_steps'
         history_column_name(h_diffusion_solver_iters) = 'diffusion_solver_iters'
         history_column_name(h_total_radiation) = 'total_radiation'
         history_column_name(h_total_energy_plus_total_radiation) = &
            'total_energy_plus_total_radiation'

         history_column_name(h_tot_E) = 'tot_E'
         history_column_name(h_log_tot_E) = 'log_tot_E'

         history_column_name(h_tot_KE) = 'tot_KE'
         history_column_name(h_log_tot_KE) = 'log_tot_KE'

         history_column_name(h_tot_PE) = 'tot_PE'
         history_column_name(h_log_tot_PE) = 'log_tot_PE'

         history_column_name(h_tot_IE) = 'tot_IE'
         history_column_name(h_log_tot_IE) = 'log_tot_IE'

         history_column_name(h_tot_Et) = 'tot_Et'
         history_column_name(h_log_tot_Et) = 'log_tot_Et'

         history_column_name(h_tot_IE_div_IE_plus_KE) = 'tot_IE_div_IE_plus_KE'

         history_column_name(h_grav_dark_L_polar) = 'grav_dark_L_polar'
         history_column_name(h_grav_dark_Teff_polar) = 'grav_dark_Teff_polar'
         
         history_column_name(h_grav_dark_L_equatorial) = 'grav_dark_L_equatorial'
         history_column_name(h_grav_dark_Teff_equatorial) = 'grav_dark_Teff_equatorial'

         history_column_name(h_apsidal_constant_k2) = 'apsidal_constant_k2'
         
! items corresponding to names on terminal output lines
         history_column_name(h_lg_Lnuc) = 'lg_Lnuc'
         history_column_name(h_H_rich) = 'H_rich'
         history_column_name(h_N_cntr) = 'N_cntr'
         history_column_name(h_lg_Lneu) = 'lg_Lneu'
         history_column_name(h_He_core) = 'He_core'
         history_column_name(h_O_cntr) = 'O_cntr'
         history_column_name(h_lg_Lphoto) = 'lg_Lphoto'
         history_column_name(h_C_core) = 'C_core'
         history_column_name(h_O_core) = 'O_core'
         history_column_name(h_Si_core) = 'Si_core'
         history_column_name(h_Fe_core) = 'Fe_core'
         history_column_name(h_Ne_cntr) = 'Ne_cntr'
         history_column_name(h_Mass) = 'Mass'
         history_column_name(h_H_cntr) = 'H_cntr'
         history_column_name(h_Si_cntr) = 'Si_cntr'
         history_column_name(h_lg_Mdot) = 'lg_Mdot'
         history_column_name(h_He_cntr) = 'He_cntr'
         history_column_name(h_eta_cntr) = 'eta_cntr'
         history_column_name(h_gam_cntr) = 'gam_cntr'
         history_column_name(h_v_div_cs) = 'v_div_cs'
         history_column_name(h_zones) = 'zones'
         history_column_name(h_lg_Dsurf) = 'lg_Dsurf'
         history_column_name(h_C_cntr) = 'C_cntr'
         history_column_name(h_retries) = 'retries'

         cnt = 0
         do i=1,h_col_id_max
            if (len_trim(history_column_name(i)) == 0) then
               write(*,*) 'missing name for log column id', i
               if (i > 1) write(*,*) 'following ' // trim(history_column_name(max(1,i-1))) ! bp: get rid of bogus compiler warning
               write(*,*)
               cnt = cnt+1
            end if
            !write(*,'(a)') '    ! ' // trim(history_column_name(i))
         end do
         !stop

         if (cnt > 0) then
            ierr = -1
            return
         end if

         nullify(history_column_names_dict)
         do i=1,h_col_id_max
            call integer_dict_define(history_column_names_dict, history_column_name(i), i, ierr)
            if (ierr /= 0) then
               write(*,*) 'FATAL ERROR: history_column_names_init failed in integer_dict_define'
               return
            end if
         end do

      end subroutine history_column_names_init


      subroutine history_column_names_shutdown()
        use utils_lib, only: integer_dict_free
        if (ASSOCIATED(history_column_names_dict)) call integer_dict_free(history_column_names_dict)
      end subroutine history_column_names_shutdown


      integer function do_get_history_id(cname)
         use utils_lib
         character (len=*), intent(in)  :: cname
         ! returns id for the log column if there is a matching name
         ! returns 0 otherwise.
         integer :: ierr, value
         call integer_dict_lookup(history_column_names_dict, cname, value, ierr)
         if (ierr /= 0) value = 0
         do_get_history_id = value
      end function do_get_history_id




      end module star_history_def

