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

      module star_profile_def

      use star_def

      implicit none

   ! profile column options

      ! abundances -- names from chem_Name array in chem_def
      ! nuclear reaction eps - names from reaction_Name array in rates_def
      ! nuclear reaction category eps - names from category_name array in net_def
      ! items specified in another file - include 'filename'
      ! star_data items from the following list

      integer, parameter :: p_zone = 1
      integer, parameter :: p_k = p_zone + 1
      integer, parameter :: p_log_Lrad_div_Ledd = p_k + 1
      integer, parameter :: p_log_Lconv_div_L = p_log_Lrad_div_Ledd + 1
      integer, parameter :: p_log_Lrad_div_L = p_log_Lconv_div_L + 1
      integer, parameter :: p_log_Lconv = p_log_Lrad_div_L + 1
      integer, parameter :: p_log_Lrad = p_log_Lconv + 1
      integer, parameter :: p_lum_conv_MLT = p_log_Lrad + 1
      integer, parameter :: p_lum_conv = p_lum_conv_MLT + 1
      integer, parameter :: p_lum_conv_div_lum_Edd = p_lum_conv + 1
      integer, parameter :: p_lum_rad_div_L_Edd = p_lum_conv_div_lum_Edd + 1
      integer, parameter :: p_lum_rad_div_L = p_lum_rad_div_L_Edd + 1

      integer, parameter :: p_lum_conv_div_L = p_lum_rad_div_L + 1
      integer, parameter :: p_lum_conv_div_lum_rad = p_lum_conv_div_L + 1
      integer, parameter :: p_lum_plus_lum_adv = p_lum_conv_div_lum_rad + 1
      integer, parameter :: p_lum_adv = p_lum_plus_lum_adv + 1
      integer, parameter :: p_lum_rad = p_lum_adv + 1
      integer, parameter :: p_log_abs_lum_erg_s = p_lum_rad + 1
      integer, parameter :: p_lum_erg_s = p_log_abs_lum_erg_s + 1
      integer, parameter :: p_luminosity = p_lum_erg_s + 1
      integer, parameter :: p_log_g = p_luminosity + 1
      integer, parameter :: p_grav = p_log_g + 1
      integer, parameter :: p_r_div_g = p_grav + 1
      integer, parameter :: p_g_div_r = p_r_div_g + 1
      integer, parameter :: p_eps_nuc_plus_nuc_neu = p_g_div_r + 1

      integer, parameter :: p_eps_nuc_minus_non_nuc_neu = p_eps_nuc_plus_nuc_neu + 1
      integer, parameter :: p_net_nuclear_energy = p_eps_nuc_minus_non_nuc_neu + 1
      integer, parameter :: p_net_energy = p_net_nuclear_energy + 1
      integer, parameter :: p_logL = p_net_energy + 1
      integer, parameter :: p_log_Ledd = p_logL + 1
      integer, parameter :: p_log_L_div_Ledd = p_log_Ledd + 1
      integer, parameter :: p_lum_div_Ledd = p_log_L_div_Ledd + 1
      integer, parameter :: p_signed_log_power = p_lum_div_Ledd + 1
      integer, parameter :: p_vel_km_per_s = p_signed_log_power + 1
      integer, parameter :: p_log_abs_dvdt_div_v = p_vel_km_per_s + 1
      integer, parameter :: p_log_abs_v = p_log_abs_dvdt_div_v + 1
      integer, parameter :: p_superad_reduction_factor = p_log_abs_v + 1
      integer, parameter :: p_gradT_excess_effect = p_superad_reduction_factor + 1
      integer, parameter :: p_log_diff_grads = p_gradT_excess_effect + 1
      integer, parameter :: p_diff_grads = p_log_diff_grads + 1

      integer, parameter :: p_v = p_diff_grads + 1
      integer, parameter :: p_velocity = p_v + 1
      integer, parameter :: p_rmid = p_velocity + 1
      integer, parameter :: p_logR_cm = p_rmid + 1
      integer, parameter :: p_radius_km = p_logR_cm + 1
      integer, parameter :: p_radius_cm = p_radius_km + 1
      integer, parameter :: p_radius = p_radius_cm + 1
      integer, parameter :: p_logR = p_radius + 1
      integer, parameter :: p_log_q = p_logR + 1
      integer, parameter :: p_q = p_log_q + 1
      integer, parameter :: p_log_dq = p_q + 1
      integer, parameter :: p_dq = p_log_dq + 1

      integer, parameter :: p_avg_charge_H = p_dq + 1
      integer, parameter :: p_avg_charge_He = p_avg_charge_H + 1
      integer, parameter :: p_avg_charge_C = p_avg_charge_He + 1
      integer, parameter :: p_avg_charge_N = p_avg_charge_C + 1
      integer, parameter :: p_avg_charge_O = p_avg_charge_N + 1
      integer, parameter :: p_avg_charge_Ne = p_avg_charge_O + 1
      integer, parameter :: p_avg_charge_Mg = p_avg_charge_Ne + 1
      integer, parameter :: p_avg_charge_Si = p_avg_charge_Mg + 1
      integer, parameter :: p_avg_charge_Fe = p_avg_charge_Si + 1

      integer, parameter :: p_neutral_fraction_H = p_avg_charge_Fe + 1
      integer, parameter :: p_neutral_fraction_He = p_neutral_fraction_H + 1
      integer, parameter :: p_neutral_fraction_C = p_neutral_fraction_He + 1
      integer, parameter :: p_neutral_fraction_N = p_neutral_fraction_C + 1
      integer, parameter :: p_neutral_fraction_O = p_neutral_fraction_N + 1
      integer, parameter :: p_neutral_fraction_Ne = p_neutral_fraction_O + 1
      integer, parameter :: p_neutral_fraction_Mg = p_neutral_fraction_Ne + 1
      integer, parameter :: p_neutral_fraction_Si = p_neutral_fraction_Mg + 1
      integer, parameter :: p_neutral_fraction_Fe = p_neutral_fraction_Si + 1

      integer, parameter :: p_pgas_div_p = p_neutral_fraction_Fe + 1
      integer, parameter :: p_pgas_div_ptotal = p_pgas_div_p + 1

      integer, parameter :: p_prad_div_pgas = p_pgas_div_ptotal + 1
      integer, parameter :: p_prad_div_pgas_div_L_div_Ledd = p_prad_div_pgas + 1

      integer, parameter :: p_log_dmbar_m_div_r = p_prad_div_pgas_div_L_div_Ledd + 1
      integer, parameter :: p_dmbar_m_div_r = p_log_dmbar_m_div_r + 1
      integer, parameter :: p_m_div_r = p_dmbar_m_div_r + 1
      integer, parameter :: p_mass_grams = p_m_div_r + 1
      integer, parameter :: p_log_mass = p_mass_grams + 1
      integer, parameter :: p_mass = p_log_mass + 1
      integer, parameter :: p_mmid = p_mass + 1

      integer, parameter :: p_dm = p_mmid + 1
      integer, parameter :: p_dm_bar = p_dm + 1

      integer, parameter :: p_m_grav = p_dm_bar + 1
      integer, parameter :: p_mass_correction_factor = p_m_grav + 1
      integer, parameter :: p_m_grav_div_m_baryonic = p_mass_correction_factor + 1

      integer, parameter :: p_xr = p_m_grav_div_m_baryonic + 1
      integer, parameter :: p_xr_cm = p_xr + 1
      integer, parameter :: p_xr_div_R = p_xr_cm + 1
      integer, parameter :: p_log_xr = p_xr_div_R + 1
      integer, parameter :: p_log_xr_cm = p_log_xr + 1
      integer, parameter :: p_log_xr_div_R = p_log_xr_cm + 1

      integer, parameter :: p_xm = p_log_xr_div_R + 1
      integer, parameter :: p_xq = p_xm + 1
      integer, parameter :: p_logxq = p_xq + 1
      integer, parameter :: p_logxm = p_logxq + 1

      integer, parameter :: p_log_radial_depth = p_logxm + 1
      integer, parameter :: p_log_column_depth = p_log_radial_depth + 1
      integer, parameter :: p_dRstar_div_dr =p_log_column_depth + 1
      integer, parameter :: p_dr_div_R =p_dRstar_div_dr + 1

      integer, parameter :: p_log_dr_div_R = p_dr_div_R + 1
      integer, parameter :: p_r_div_R = p_log_dr_div_R + 1
      integer, parameter :: p_log_dr_div_rmid = p_r_div_R + 1
      integer, parameter :: p_dr_div_rmid = p_log_dr_div_rmid + 1
      integer, parameter :: p_dlogR = p_dr_div_rmid + 1
      integer, parameter :: p_log_dr = p_dlogR + 1

      integer, parameter :: p_lnR_residual = p_log_dr + 1
      integer, parameter :: p_lnd_residual = p_lnR_residual + 1
      integer, parameter :: p_equL_residual = p_lnd_residual + 1
      integer, parameter :: p_cell_internal_energy_fraction_start = p_equL_residual + 1
      integer, parameter :: p_cell_internal_energy_fraction = p_cell_internal_energy_fraction_start + 1
      
      integer, parameter :: p_log_rel_E_err = p_cell_internal_energy_fraction + 1
      integer, parameter :: p_ergs_error_integral = p_log_rel_E_err + 1
      integer, parameter :: p_ergs_rel_error_integral = p_ergs_error_integral + 1
      integer, parameter :: p_ergs_error = p_ergs_rel_error_integral + 1
      integer, parameter :: p_E_residual = p_ergs_error + 1
      integer, parameter :: p_Et_residual = p_E_residual + 1
      integer, parameter :: p_log_Et_residual = p_Et_residual + 1
      integer, parameter :: p_dvdt_residual = p_log_Et_residual + 1
      integer, parameter :: p_v_residual = p_dvdt_residual + 1
      integer, parameter :: p_log_E_residual = p_v_residual + 1
      integer, parameter :: p_log_lnR_residual = p_log_E_residual + 1
      integer, parameter :: p_log_lnd_residual = p_log_lnR_residual + 1
      integer, parameter :: p_log_equL_residual = p_log_lnd_residual + 1
      integer, parameter :: p_log_dvdt_residual = p_log_equL_residual + 1
      integer, parameter :: p_log_v_residual = p_log_dvdt_residual + 1

      integer, parameter :: p_t_rad = p_log_v_residual + 1
      integer, parameter :: p_log_t_rad = p_t_rad + 1
      integer, parameter :: p_log_dt_cs_div_dr = p_log_t_rad + 1
      integer, parameter :: p_dt_cs_div_dr = p_log_dt_cs_div_dr + 1
      integer, parameter :: p_dr_div_cs = p_dt_cs_div_dr + 1
      integer, parameter :: p_log_dr_div_cs = p_dr_div_cs + 1
      integer, parameter :: p_dr_div_cs_yr = p_log_dr_div_cs + 1
      integer, parameter :: p_log_dr_div_cs_yr = p_dr_div_cs_yr + 1

      integer, parameter :: p_cell_collapse_time = p_log_dr_div_cs_yr + 1
      integer, parameter :: p_log_cell_collapse_time = p_cell_collapse_time + 1

      integer, parameter :: p_log_acoustic_depth = p_log_cell_collapse_time + 1
      integer, parameter :: p_log_acoustic_radius = p_log_acoustic_depth + 1
      integer, parameter :: p_acoustic_depth = p_log_acoustic_radius + 1
      integer, parameter :: p_acoustic_radius = p_acoustic_depth + 1
      integer, parameter :: p_acoustic_r_div_R_phot = p_acoustic_radius + 1

      integer, parameter :: p_x = p_acoustic_r_div_R_phot + 1
      integer, parameter :: p_log_x = p_x + 1
      integer, parameter :: p_y = p_log_x + 1
      integer, parameter :: p_log_y = p_y + 1
      integer, parameter :: p_z = p_log_y + 1
      integer, parameter :: p_log_z = p_z + 1

      integer, parameter :: p_compression_gradient = p_log_z + 1
      integer, parameter :: p_logdq = p_compression_gradient + 1
      integer, parameter :: p_dq_ratio = p_logdq + 1
      integer, parameter :: p_tau_eff = p_dq_ratio + 1
      integer, parameter :: p_tau_eff_div_tau = p_tau_eff + 1
      integer, parameter :: p_xtau = p_tau_eff_div_tau + 1
      integer, parameter :: p_tau = p_xtau + 1
      integer, parameter :: p_kap_frac_Type2 = p_tau + 1
      integer, parameter :: p_kap_frac_op_mono = p_kap_frac_Type2 + 1
      integer, parameter :: p_extra_opacity_factor = p_kap_frac_op_mono + 1
      integer, parameter :: p_log_kap_times_factor = p_extra_opacity_factor + 1
      integer, parameter :: p_log_opacity = p_log_kap_times_factor + 1
      integer, parameter :: p_energy = p_log_opacity + 1
      integer, parameter :: p_logM = p_energy + 1
      integer, parameter :: p_logtau_sub_xlogtau = p_logM + 1
      integer, parameter :: p_xlogtau = p_logtau_sub_xlogtau + 1
      integer, parameter :: p_logtau = p_xlogtau + 1
      integer, parameter :: p_temperature = p_logtau + 1
      integer, parameter :: p_logT_face = p_temperature + 1
      integer, parameter :: p_logT_bb = p_logT_face + 1
      integer, parameter :: p_logT_face_div_logT_bb = p_logT_bb + 1
      integer, parameter :: p_logT = p_logT_face_div_logT_bb + 1
      integer, parameter :: p_density = p_logT + 1
      integer, parameter :: p_rho = p_density + 1
      integer, parameter :: p_logRho = p_rho + 1
      integer, parameter :: p_pgas = p_logRho+ 1
      integer, parameter :: p_logPgas = p_pgas + 1
      integer, parameter :: p_prad = p_logPgas + 1
      integer, parameter :: p_pressure = p_prad + 1
      integer, parameter :: p_logP = p_pressure + 1
      integer, parameter :: p_logE = p_logP + 1
      integer, parameter :: p_grada = p_logE + 1
      integer, parameter :: p_dE_dRho = p_grada + 1
      integer, parameter :: p_Cv = p_dE_dRho + 1
      integer, parameter :: p_cp = p_Cv + 1
      integer, parameter :: p_thermal_time_to_surface = p_cp + 1
      integer, parameter :: p_log_thermal_time_to_surface = p_thermal_time_to_surface + 1
      integer, parameter :: p_log_CpT = p_log_thermal_time_to_surface + 1
      integer, parameter :: p_log_CpT_absMdot_div_L = p_log_CpT + 1
      integer, parameter :: p_logS = p_log_CpT_absMdot_div_L + 1
      integer, parameter :: p_logS_per_baryon = p_logS + 1
      integer, parameter :: p_gamma1 = p_logS_per_baryon + 1
      integer, parameter :: p_gamma3 = p_gamma1 + 1
      integer, parameter :: p_eta = p_gamma3 + 1
      integer, parameter :: p_theta_e = p_eta + 1
      integer, parameter :: p_gam = p_theta_e + 1
      integer, parameter :: p_mu = p_gam + 1

      integer, parameter :: p_eos_frac_OPAL_SCVH = p_mu + 1
      integer, parameter :: p_eos_frac_HELM = p_eos_frac_OPAL_SCVH + 1
      integer, parameter :: p_eos_frac_Skye = p_eos_frac_HELM + 1
      integer, parameter :: p_eos_frac_PC = p_eos_frac_Skye + 1
      integer, parameter :: p_eos_frac_FreeEOS = p_eos_frac_PC + 1
      integer, parameter :: p_eos_frac_CMS = p_eos_frac_FreeEOS + 1

      integer, parameter :: p_log_rho_times_r3 = p_eos_frac_CMS + 1
      integer, parameter :: p_rho_times_r3 = p_log_rho_times_r3 + 1
      integer, parameter :: p_v_times_t_div_r = p_rho_times_r3 + 1
      integer, parameter :: p_d_v_div_r = p_v_times_t_div_r + 1
      integer, parameter :: p_v_div_r = p_d_v_div_r + 1

      integer, parameter :: p_log_c_div_tau = p_v_div_r + 1
      integer, parameter :: p_log_v_escape = p_log_c_div_tau + 1
      integer, parameter :: p_v_escape = p_log_v_escape + 1
      integer, parameter :: p_v_div_vesc = p_v_escape + 1
      integer, parameter :: p_v_div_v_escape = p_v_div_vesc + 1
      integer, parameter :: p_v_div_cs = p_v_div_v_escape + 1
      integer, parameter :: p_v_div_csound = p_v_div_cs + 1
      integer, parameter :: p_csound_face = p_v_div_csound + 1
      integer, parameter :: p_log_csound = p_csound_face + 1
      integer, parameter :: p_csound = p_log_csound + 1
      integer, parameter :: p_scale_height = p_csound + 1
      integer, parameter :: p_omega = p_scale_height + 1

      integer, parameter :: p_log_omega = p_omega + 1
      integer, parameter :: p_log_j_rot = p_log_omega + 1

      integer, parameter :: p_log_J_inside = p_log_j_rot + 1
      integer, parameter :: p_log_J_div_M53 = p_log_J_inside + 1
      integer, parameter :: p_log_abs_dlnR_domega = p_log_J_div_M53 + 1
      integer, parameter :: p_log_abs_shear = p_log_abs_dlnR_domega + 1
      integer, parameter :: p_shear = p_log_abs_shear + 1
      integer, parameter :: p_i_rot = p_shear + 1
      integer, parameter :: p_j_rot = p_i_rot + 1
      integer, parameter :: p_v_rot = p_j_rot + 1
      integer, parameter :: p_fp_rot = p_v_rot + 1


      integer, parameter :: p_ft_rot = p_fp_rot + 1
      integer, parameter :: p_ft_rot_div_fp_rot = p_ft_rot + 1
      integer, parameter :: p_w_div_w_crit_roche = p_ft_rot_div_fp_rot + 1
      integer, parameter :: p_w_div_w_crit_roche2 = p_w_div_w_crit_roche + 1
      integer, parameter :: p_log_am_nu_rot = p_w_div_w_crit_roche2 + 1
      integer, parameter :: p_log_am_nu_non_rot = p_log_am_nu_rot + 1
      integer, parameter :: p_log_am_nu = p_log_am_nu_non_rot + 1

      integer, parameter :: p_eps_WD_sedimentation = p_log_am_nu + 1
      integer, parameter :: p_log_eps_WD_sedimentation = p_eps_WD_sedimentation + 1

      integer, parameter :: p_eps_diffusion = p_log_eps_WD_sedimentation + 1
      integer, parameter :: p_log_eps_diffusion = p_eps_diffusion + 1
      
      integer, parameter :: p_log_e_field = p_log_eps_diffusion + 1
      integer, parameter :: p_e_field = p_log_e_field + 1
      integer, parameter :: p_log_g_field_element_diffusion = p_e_field + 1
      integer, parameter :: p_g_field_element_diffusion = p_log_g_field_element_diffusion + 1

      integer, parameter :: p_log_eE_div_mg_element_diffusion = p_g_field_element_diffusion + 1
      integer, parameter :: p_eE_div_mg_element_diffusion = p_log_eE_div_mg_element_diffusion + 1

      integer, parameter :: p_r_polar = p_eE_div_mg_element_diffusion + 1
      integer, parameter :: p_log_r_polar = p_r_polar + 1
      integer, parameter :: p_r_equatorial = p_log_r_polar + 1
      integer, parameter :: p_log_r_equatorial = p_r_equatorial + 1
      integer, parameter :: p_r_e_div_r_p = p_log_r_equatorial + 1

      integer, parameter :: p_omega_crit = p_r_e_div_r_p + 1
      integer, parameter :: p_omega_div_omega_crit = p_omega_crit + 1

      integer, parameter :: p_am_log_sig = p_omega_div_omega_crit + 1
      integer, parameter :: p_am_log_sig_omega = p_am_log_sig + 1
      integer, parameter :: p_am_log_sig_j = p_am_log_sig_omega + 1

      integer, parameter :: p_am_log_nu_omega = p_am_log_sig_j + 1
      integer, parameter :: p_am_log_nu_j = p_am_log_nu_omega + 1
      integer, parameter :: p_am_log_nu_rot = p_am_log_nu_j + 1
      integer, parameter :: p_am_log_nu_non_rot = p_am_log_nu_rot + 1

      integer, parameter :: p_richardson_number = p_am_log_nu_non_rot + 1
      integer, parameter :: p_am_domega_dlnR = p_richardson_number + 1
      integer, parameter :: p_am_log_D_visc = p_am_domega_dlnR + 1
      integer, parameter :: p_am_log_D_DSI = p_am_log_D_visc + 1
      integer, parameter :: p_am_log_D_SH = p_am_log_D_DSI + 1
      integer, parameter :: p_am_log_D_SSI = p_am_log_D_SH + 1
      integer, parameter :: p_am_log_D_ES = p_am_log_D_SSI + 1
      integer, parameter :: p_am_log_D_GSF = p_am_log_D_ES + 1
      integer, parameter :: p_am_log_D_ST = p_am_log_D_GSF + 1
      integer, parameter :: p_am_log_nu_ST = p_am_log_D_ST + 1

      integer, parameter :: p_dynamo_log_B_r = p_am_log_nu_ST + 1
      integer, parameter :: p_dynamo_log_B_phi = p_dynamo_log_B_r + 1
      integer, parameter :: p_grada_face = p_dynamo_log_B_phi + 1

      integer, parameter :: p_gradr_div_grada = p_grada_face + 1
      integer, parameter :: p_gradr_sub_grada = p_gradr_div_grada + 1
      integer, parameter :: p_entropy = p_gradr_sub_grada + 1
      integer, parameter :: p_free_e = p_entropy + 1
      integer, parameter :: p_logfree_e = p_free_e + 1
      integer, parameter :: p_chiRho = p_logfree_e + 1
      integer, parameter :: p_chiT = p_chiRho + 1
      integer, parameter :: p_chiRho_for_partials = p_chiT + 1
      integer, parameter :: p_chiT_for_partials = p_chiRho_for_partials + 1
      integer, parameter :: p_rel_diff_chiRho_for_partials = p_chiT_for_partials + 1
      integer, parameter :: p_rel_diff_chiT_for_partials = p_rel_diff_chiRho_for_partials + 1
      integer, parameter :: p_QQ = p_rel_diff_chiT_for_partials + 1

      integer, parameter :: p_phase = p_QQ + 1
      integer, parameter :: p_latent_ddlnT = p_phase + 1
      integer, parameter :: p_latent_ddlnRho = p_latent_ddlnT + 1

      integer, parameter :: p_x_mass_fraction_H = p_latent_ddlnRho + 1
      integer, parameter :: p_y_mass_fraction_He = p_x_mass_fraction_H + 1
      integer, parameter :: p_z_mass_fraction_metals = p_y_mass_fraction_He + 1

      integer, parameter :: p_abar = p_z_mass_fraction_metals + 1
      integer, parameter :: p_zbar = p_abar + 1
      integer, parameter :: p_z2bar = p_zbar + 1
      integer, parameter :: p_ye = p_z2bar + 1

      integer, parameter :: p_dkap_dlnrho_face = p_ye + 1
      integer, parameter :: p_dkap_dlnT_face = p_dkap_dlnrho_face + 1
      integer, parameter :: p_opacity = p_dkap_dlnt_face + 1
      
      integer, parameter :: p_deps_dlnd_face = p_opacity + 1
      integer, parameter :: p_deps_dlnT_face = p_deps_dlnd_face + 1
      integer, parameter :: p_d_epsnuc_dlnd = p_deps_dlnT_face + 1
      integer, parameter :: p_d_lnepsnuc_dlnd = p_d_epsnuc_dlnd + 1
      integer, parameter :: p_d_lnepsnuc_dlnT = p_d_lnepsnuc_dlnd + 1
      integer, parameter :: p_d_epsnuc_dlnT = p_d_lnepsnuc_dlnT + 1
      integer, parameter :: p_eps_nuc_start = p_d_epsnuc_dlnT + 1
      integer, parameter :: p_signed_log_eps_nuc = p_eps_nuc_start + 1 


      integer, parameter :: p_log_abs_eps_nuc = p_signed_log_eps_nuc + 1
      integer, parameter :: p_eps_nuc = p_log_abs_eps_nuc + 1

      integer, parameter :: p_burn_sum_xa_err_before_fix = p_eps_nuc + 1
      integer, parameter :: p_burn_log_abs_sum_xa_err_before_fix = p_burn_sum_xa_err_before_fix + 1

      integer, parameter :: p_eps_nuc_neu_total = p_burn_log_abs_sum_xa_err_before_fix + 1
      integer, parameter :: p_non_nuc_neu = p_eps_nuc_neu_total + 1

      integer, parameter :: p_nonnucneu_plas = p_non_nuc_neu + 1
      integer, parameter :: p_nonnucneu_brem = p_nonnucneu_plas + 1
      integer, parameter :: p_nonnucneu_phot = p_nonnucneu_brem + 1
      integer, parameter :: p_nonnucneu_pair = p_nonnucneu_phot + 1
      integer, parameter :: p_nonnucneu_reco = p_nonnucneu_pair + 1
      integer, parameter :: p_log_irradiation_heat = p_nonnucneu_reco + 1
      integer, parameter :: p_alpha_mlt = p_log_irradiation_heat + 1
      integer, parameter :: p_cgrav_factor = p_alpha_mlt + 1
      integer, parameter :: p_log_extra_L = p_cgrav_factor + 1
      integer, parameter :: p_extra_L = p_log_extra_L + 1
      integer, parameter :: p_extra_jdot = p_extra_L + 1
      integer, parameter :: p_extra_omegadot = p_extra_jdot + 1
      integer, parameter :: p_extra_grav = p_extra_omegadot + 1
      integer, parameter :: p_extra_heat = p_extra_grav + 1
      integer, parameter :: p_div_v = p_extra_heat + 1
      integer, parameter :: p_d_v_div_r_dm = p_div_v + 1 
      integer, parameter :: p_d_v_div_r_dr = p_d_v_div_r_dm + 1

      integer, parameter :: p_dvdt_grav = p_d_v_div_r_dr + 1
      integer, parameter :: p_dvdt_dPdm = p_dvdt_grav + 1
      integer, parameter :: p_du = p_dvdt_dPdm + 1
      integer, parameter :: p_P_face = p_du + 1

      integer, parameter :: p_log_P_face = p_P_face + 1

      integer, parameter :: p_hse_ratio = p_log_P_face + 1
      integer, parameter :: p_hse_ratio_gyre = p_hse_ratio + 1

      integer, parameter :: p_dvdt_RTI_diffusion = p_hse_ratio_gyre + 1
      integer, parameter :: p_dlnddt_RTI_diffusion = p_dvdt_RTI_diffusion + 1

      integer, parameter :: p_dlnP_dlnR = p_dlnddt_RTI_diffusion + 1
      integer, parameter :: p_dlnRho_dlnR = p_dlnP_dlnR + 1

      integer, parameter :: p_gradP_div_rho = p_dlnRho_dlnR + 1
      integer, parameter :: p_dPdr_div_grav = p_gradP_div_rho + 1
      integer, parameter :: p_env_eps_grav = p_dPdr_div_grav + 1

      integer, parameter :: p_log_abs_eps_grav_dm_div_L = p_env_eps_grav + 1
      integer, parameter :: p_eps_grav_composition_term = p_log_abs_eps_grav_dm_div_L + 1

      integer, parameter :: p_eps_grav_plus_eps_mdot = p_eps_grav_composition_term + 1
      integer, parameter :: p_ergs_eps_grav_plus_eps_mdot = p_eps_grav_plus_eps_mdot + 1

      integer, parameter :: p_ergs_mdot = p_ergs_eps_grav_plus_eps_mdot + 1
      integer, parameter :: p_eps_mdot = p_ergs_mdot + 1

      integer, parameter :: p_dm_eps_grav = p_eps_mdot + 1
      
      integer, parameter :: p_log_xm_div_delta_m = p_dm_eps_grav + 1
      integer, parameter :: p_xm_div_delta_m = p_log_xm_div_delta_m + 1
      integer, parameter :: p_eps_grav = p_xm_div_delta_m + 1
      integer, parameter :: p_signed_log_eps_grav = p_eps_grav + 1
      integer, parameter :: p_log_conv_L_div_L = p_signed_log_eps_grav + 1
      integer, parameter :: p_conv_L_div_L = p_log_conv_L_div_L + 1
      integer, parameter :: p_mlt_Zeta = p_conv_L_div_L + 1
      integer, parameter :: p_mlt_Gamma = p_mlt_Zeta + 1
      integer, parameter :: p_mlt_mixing_length = p_mlt_Gamma + 1
      integer, parameter :: p_mlt_Pturb = p_mlt_mixing_length + 1
      integer, parameter :: p_mlt_mixing_type = p_mlt_Pturb + 1

      integer, parameter :: p_grada_sub_gradT = p_mlt_mixing_type + 1
      integer, parameter :: p_gradT_sub_a = p_grada_sub_gradT + 1
      integer, parameter :: p_gradT_div_grada = p_gradT_sub_a + 1

      integer, parameter :: p_gradT_rel_err = p_gradT_div_grada + 1
      integer, parameter :: p_gradr_sub_gradT = p_gradT_rel_err + 1
      integer, parameter :: p_gradT_sub_gradr = p_gradr_sub_gradT + 1
      integer, parameter :: p_gradT_div_gradr = p_gradT_sub_gradr + 1
      integer, parameter :: p_log_gradT_div_gradr = p_gradT_div_gradr + 1

      integer, parameter :: p_log_mlt_Gamma = p_log_gradT_div_gradr + 1
      integer, parameter :: p_log_mlt_vc = p_log_mlt_Gamma + 1
      integer, parameter :: p_conv_vel_div_mlt_vc = p_log_mlt_vc + 1
      integer, parameter :: p_mlt_vc = p_conv_vel_div_mlt_vc + 1
      integer, parameter :: p_super_ad = p_mlt_vc + 1

      integer, parameter :: p_delta_r = p_super_ad + 1
      integer, parameter :: p_delta_L = p_delta_r + 1
      integer, parameter :: p_delta_cell_vol = p_delta_L + 1
      integer, parameter :: p_delta_entropy = p_delta_cell_vol + 1
      integer, parameter :: p_delta_T = p_delta_entropy + 1
      integer, parameter :: p_delta_rho = p_delta_T + 1
      integer, parameter :: p_delta_eps_nuc = p_delta_rho + 1
      integer, parameter :: p_delta_mu = p_delta_eps_nuc + 1
      integer, parameter :: p_log_D_conv = p_delta_mu + 1

      integer, parameter :: p_log_D_smooth = p_log_D_conv + 1
      integer, parameter :: p_log_D_leftover = p_log_D_smooth + 1
      integer, parameter :: p_log_D_semi = p_log_D_leftover + 1
      integer, parameter :: p_log_D_anon = p_log_D_semi + 1
      integer, parameter :: p_log_D_ovr = p_log_D_anon + 1
      integer, parameter :: p_log_D_thrm = p_log_D_ovr + 1
      integer, parameter :: p_log_D_rayleigh_taylor = p_log_D_thrm + 1
      integer, parameter :: p_log_D_minimum = p_log_D_rayleigh_taylor + 1
      integer, parameter :: p_log_D_omega = p_log_D_minimum + 1

      integer, parameter :: p_log_D_mix_rotation = p_log_D_omega + 1
      integer, parameter :: p_log_D_mix_non_rotation = p_log_D_mix_rotation + 1
      integer, parameter :: p_log_D_mix = p_log_D_mix_non_rotation + 1
      integer, parameter :: p_conv_vel = p_log_D_mix + 1
      integer, parameter :: p_dt_times_conv_vel_div_mixing_length = p_conv_vel + 1
      integer, parameter :: p_log_dt_times_conv_vel_div_mixing_length = p_dt_times_conv_vel_div_mixing_length + 1
      
      integer, parameter :: p_log_lambda_RTI_div_Hrho = p_log_dt_times_conv_vel_div_mixing_length + 1
      integer, parameter :: p_lambda_RTI = p_log_lambda_RTI_div_Hrho + 1
      integer, parameter :: p_dPdr_info = p_lambda_RTI + 1
      integer, parameter :: p_dRhodr_info = p_dPdr_info + 1
      
      integer, parameter :: p_source_plus_alpha_RTI = p_dRhodr_info + 1
      integer, parameter :: p_log_source_RTI = p_source_plus_alpha_RTI + 1
      integer, parameter :: p_log_source_plus_alpha_RTI = p_log_source_RTI + 1
      integer, parameter :: p_source_minus_alpha_RTI = p_log_source_plus_alpha_RTI + 1
      integer, parameter :: p_log_source_minus_alpha_RTI = p_source_minus_alpha_RTI + 1

      integer, parameter :: p_dudt_RTI = p_log_source_minus_alpha_RTI + 1
      integer, parameter :: p_dedt_RTI = p_dudt_RTI + 1
      integer, parameter :: p_eta_RTI = p_dedt_RTI + 1
      integer, parameter :: p_log_eta_RTI = p_eta_RTI + 1

      integer, parameter :: p_boost_for_eta_RTI = p_log_eta_RTI + 1
      integer, parameter :: p_log_boost_for_eta_RTI = p_boost_for_eta_RTI + 1

      integer, parameter :: p_alpha_RTI = p_log_boost_for_eta_RTI + 1
      integer, parameter :: p_log_alpha_RTI = p_alpha_RTI + 1
      integer, parameter :: p_log_etamid_RTI = p_log_alpha_RTI + 1
      integer, parameter :: p_log_sig_RTI = p_log_etamid_RTI + 1
      integer, parameter :: p_log_sigmid_RTI = p_log_sig_RTI + 1

      integer, parameter :: p_burn_avg_epsnuc = p_log_sigmid_RTI + 1
      integer, parameter :: p_log_burn_avg_epsnuc = p_burn_avg_epsnuc + 1
      integer, parameter :: p_burn_num_iters = p_log_burn_avg_epsnuc + 1

      integer, parameter :: p_log_sig_raw_mix = p_burn_num_iters + 1
      integer, parameter :: p_log_sig_mix = p_log_sig_raw_mix + 1
      integer, parameter :: p_log_conv_vel = p_log_sig_mix + 1
      integer, parameter :: p_conv_vel_div_L_vel = p_log_conv_vel + 1
      integer, parameter :: p_conv_vel_div_csound = p_conv_vel_div_L_vel + 1
      integer, parameter :: p_log_tau_conv_yrs = p_conv_vel_div_csound + 1
      integer, parameter :: p_mixing_type = p_log_tau_conv_yrs + 1
      integer, parameter :: p_conv_mixing_type = p_mixing_type + 1

      integer, parameter :: p_log_mlt_D_mix = p_conv_mixing_type + 1
      integer, parameter :: p_log_cp_T_div_t_sound = p_log_mlt_D_mix + 1
      integer, parameter :: p_log_t_thermal = p_log_cp_T_div_t_sound + 1
      integer, parameter :: p_log_t_sound = p_log_t_thermal + 1
      integer, parameter :: p_pressure_scale_height_cm = p_log_t_sound + 1
      integer, parameter :: p_pressure_scale_height = p_pressure_scale_height_cm + 1

      integer, parameter :: p_grada_sub_actual_gradT = p_pressure_scale_height + 1
      integer, parameter :: p_gradT_sub_actual_gradT = p_grada_sub_actual_gradT + 1
      integer, parameter :: p_actual_gradT = p_gradT_sub_actual_gradT + 1

      integer, parameter :: p_grad_superad_actual = p_actual_gradT + 1
      integer, parameter :: p_grad_superad = p_grad_superad_actual + 1
      integer, parameter :: p_gradT_sub_grada = p_grad_superad + 1
      integer, parameter :: p_gradT = p_gradT_sub_grada + 1
      integer, parameter :: p_gradr = p_gradT + 1

      integer, parameter :: p_d_gradT_dlnd00 = p_gradr + 1
      integer, parameter :: p_d_gradT_dlnT00 = p_d_gradT_dlnd00 + 1
      integer, parameter :: p_d_gradT_dlndm1 = p_d_gradT_dlnT00 + 1
      integer, parameter :: p_d_gradT_dlnTm1 = p_d_gradT_dlndm1 + 1
      integer, parameter :: p_d_gradT_dlnR = p_d_gradT_dlnTm1 + 1
      integer, parameter :: p_d_gradT_dln_cvpv0 = p_d_gradT_dlnR + 1
      integer, parameter :: p_d_gradT_dL = p_d_gradT_dln_cvpv0 + 1

      integer, parameter :: p_accel_div_grav = p_d_gradT_dL + 1

      integer, parameter :: p_dlnd_dt_const_q = p_accel_div_grav + 1
      integer, parameter :: p_dlnT_dt_const_q = p_dlnd_dt_const_q + 1

      integer, parameter :: p_dlnd = p_dlnT_dt_const_q + 1
      integer, parameter :: p_dlnT = p_dlnd + 1
      integer, parameter :: p_dlnR = p_dlnT + 1

      integer, parameter :: p_dlnd_dt = p_dlnR + 1
      integer, parameter :: p_dlnT_dt = p_dlnd_dt + 1
      integer, parameter :: p_dlnR_dt = p_dlnT_dt + 1
      integer, parameter :: p_dr_dt = p_dlnR_dt + 1
      integer, parameter :: p_du_dt = p_dr_dt + 1
      integer, parameter :: p_dv_dt = p_du_dt + 1

      integer, parameter :: p_signed_dlnd = p_dv_dt + 1
      integer, parameter :: p_signed_dlnT = p_signed_dlnd + 1

      integer, parameter :: p_dt_dm_eps_grav = p_signed_dlnT + 1

      integer, parameter :: p_dm_de = p_dt_dm_eps_grav + 1
      integer, parameter :: p_dt_dL = p_dm_de + 1

      integer, parameter :: p_ds_from_eps_grav = p_dt_dL + 1
      integer, parameter :: p_del_entropy = p_ds_from_eps_grav + 1
      integer, parameter :: p_cno_div_z = p_del_entropy + 1

      integer, parameter :: p_dE = p_cno_div_z + 1
      integer, parameter :: p_dr = p_dE + 1
      integer, parameter :: p_dv = p_dr + 1

      integer, parameter :: p_dr_ratio = p_dv + 1
      integer, parameter :: p_dt_dv_div_dr = p_dr_ratio + 1

      integer, parameter :: p_dlog_h1_dlogP = p_dt_dv_div_dr + 1
      integer, parameter :: p_dlog_he3_dlogP = p_dlog_h1_dlogP + 1
      integer, parameter :: p_dlog_he4_dlogP = p_dlog_he3_dlogP + 1
      integer, parameter :: p_dlog_c12_dlogP = p_dlog_he4_dlogP + 1
      integer, parameter :: p_dlog_c13_dlogP = p_dlog_c12_dlogP + 1
      integer, parameter :: p_dlog_n14_dlogP = p_dlog_c13_dlogP + 1
      integer, parameter :: p_dlog_o16_dlogP = p_dlog_n14_dlogP + 1
      integer, parameter :: p_dlog_ne20_dlogP = p_dlog_o16_dlogP + 1
      integer, parameter :: p_dlog_mg24_dlogP = p_dlog_ne20_dlogP + 1
      integer, parameter :: p_dlog_si28_dlogP = p_dlog_mg24_dlogP + 1

      integer, parameter :: p_dlog_pp_dlogP = p_dlog_si28_dlogP + 1
      integer, parameter :: p_dlog_cno_dlogP = p_dlog_pp_dlogP + 1
      integer, parameter :: p_dlog_3alf_dlogP = p_dlog_cno_dlogP + 1

      integer, parameter :: p_dlog_burn_c_dlogP = p_dlog_3alf_dlogP + 1
      integer, parameter :: p_dlog_burn_n_dlogP = p_dlog_burn_c_dlogP + 1
      integer, parameter :: p_dlog_burn_o_dlogP = p_dlog_burn_n_dlogP + 1

      integer, parameter :: p_dlog_burn_ne_dlogP = p_dlog_burn_o_dlogP + 1
      integer, parameter :: p_dlog_burn_na_dlogP = p_dlog_burn_ne_dlogP + 1
      integer, parameter :: p_dlog_burn_mg_dlogP = p_dlog_burn_na_dlogP + 1

      integer, parameter :: p_dlog_cc_dlogP = p_dlog_burn_mg_dlogP + 1
      integer, parameter :: p_dlog_co_dlogP = p_dlog_cc_dlogP + 1
      integer, parameter :: p_dlog_oo_dlogP = p_dlog_co_dlogP + 1

      integer, parameter :: p_dlog_burn_si_dlogP = p_dlog_oo_dlogP + 1
      integer, parameter :: p_dlog_burn_s_dlogP = p_dlog_burn_si_dlogP + 1
      integer, parameter :: p_dlog_burn_ar_dlogP = p_dlog_burn_s_dlogP + 1
      integer, parameter :: p_dlog_burn_ca_dlogP = p_dlog_burn_ar_dlogP + 1
      integer, parameter :: p_dlog_burn_ti_dlogP = p_dlog_burn_ca_dlogP + 1
      integer, parameter :: p_dlog_burn_cr_dlogP = p_dlog_burn_ti_dlogP + 1
      integer, parameter :: p_dlog_burn_fe_dlogP = p_dlog_burn_cr_dlogP + 1

      integer, parameter :: p_dlog_pnhe4_dlogP = p_dlog_burn_fe_dlogP + 1
      integer, parameter :: p_dlog_photo_dlogP = p_dlog_pnhe4_dlogP + 1
      integer, parameter :: p_dlog_other_dlogP = p_dlog_photo_dlogP + 1

      integer, parameter :: p_dlnX_dr = p_dlog_other_dlogP + 1
      integer, parameter :: p_dlnY_dr = p_dlnX_dr + 1
      integer, parameter :: p_dlnRho_dr = p_dlnY_dr + 1

      integer, parameter :: p_logR_kap = p_dlnRho_dr + 1
      integer, parameter :: p_logW = p_logR_kap + 1
      integer, parameter :: p_logV = p_logW + 1
      integer, parameter :: p_logQ = p_logV + 1
      integer, parameter :: p_log_mdot_cs = p_logQ + 1
      integer, parameter :: p_log_mdot_v = p_log_mdot_cs + 1
      integer, parameter :: p_log_L_div_CpTMdot = p_log_mdot_v + 1
      integer, parameter :: p_cs_at_cell_bdy = p_log_L_div_CpTMdot + 1

      integer, parameter :: p_total_energy_integral_outward = p_cs_at_cell_bdy + 1
      integer, parameter :: p_total_energy_integral = p_total_energy_integral_outward + 1
      integer, parameter :: p_total_energy_sign = p_total_energy_integral + 1
      integer, parameter :: p_total_energy = p_total_energy_sign + 1
      
      integer, parameter :: p_Pturb = p_total_energy + 1
      integer, parameter :: p_log_Pturb = p_Pturb + 1
      integer, parameter :: p_Eturb = p_log_Pturb + 1
      integer, parameter :: p_log_Eturb = p_Eturb + 1
      integer, parameter :: p_avQ = p_log_Eturb + 1
      integer, parameter :: p_Hp_face = p_avQ + 1
      integer, parameter :: p_Y_face = p_Hp_face + 1
      integer, parameter :: p_PII_face = p_Y_face + 1
      integer, parameter :: p_Chi = p_PII_face + 1
      integer, parameter :: p_COUPL = p_Chi + 1
      integer, parameter :: p_SOURCE = p_COUPL + 1
      integer, parameter :: p_DAMP = p_SOURCE + 1
      integer, parameter :: p_DAMPR = p_DAMP + 1
      integer, parameter :: p_Eq = p_DAMPR + 1
      integer, parameter :: p_Uq = p_Eq + 1
      integer, parameter :: p_Lr = p_Uq + 1
      integer, parameter :: p_Lr_div_L = p_Lr + 1
      integer, parameter :: p_Lc = p_Lr_div_L + 1
      integer, parameter :: p_Lc_div_L = p_Lc + 1
      integer, parameter :: p_Lt = p_Lc_div_L + 1
      integer, parameter :: p_Lt_div_L = p_Lt + 1

      integer, parameter :: p_rsp_vt_div_cs = p_Lt_div_L + 1
      integer, parameter :: p_rsp_vt = p_rsp_vt_div_cs + 1
      integer, parameter :: p_rsp_log_erad = p_rsp_vt + 1
      integer, parameter :: p_rsp_erad = p_rsp_log_erad + 1
      integer, parameter :: p_rsp_logEt = p_rsp_erad + 1
      integer, parameter :: p_rsp_Et = p_rsp_logEt + 1
      integer, parameter :: p_rsp_Pt = p_rsp_Et + 1
      integer, parameter :: p_rsp_Lr = p_rsp_Pt + 1
      integer, parameter :: p_rsp_Lc = p_rsp_Lr + 1
      integer, parameter :: p_rsp_Lt = p_rsp_Lc + 1
      integer, parameter :: p_rsp_Eq = p_rsp_Lt + 1
      integer, parameter :: p_rsp_Uq = p_rsp_Eq + 1
      integer, parameter :: p_rsp_PII_face = p_rsp_Uq + 1
      integer, parameter :: p_rsp_src_snk = p_rsp_PII_face + 1
      integer, parameter :: p_rsp_src = p_rsp_src_snk + 1
      integer, parameter :: p_rsp_sink = p_rsp_src + 1
      integer, parameter :: p_rsp_damp = p_rsp_sink + 1
      integer, parameter :: p_rsp_dampR = p_rsp_damp + 1
      integer, parameter :: p_rsp_Hp_face = p_rsp_dampR + 1
      integer, parameter :: p_rsp_Y_face = p_rsp_Hp_face + 1
      integer, parameter :: p_rsp_log_heat_exchange_timescale = p_rsp_Y_face + 1
      integer, parameter :: p_rsp_log_dt_div_heat_exchange_timescale = p_rsp_log_heat_exchange_timescale + 1
      integer, parameter :: p_rsp_heat_exchange_timescale = p_rsp_log_dt_div_heat_exchange_timescale + 1
      integer, parameter :: p_rsp_avQ = p_rsp_heat_exchange_timescale + 1
      integer, parameter :: p_rsp_Chi = p_rsp_avQ + 1
      integer, parameter :: p_rsp_gradT = p_rsp_Chi + 1
      integer, parameter :: p_rsp_Lr_div_L = p_rsp_gradT + 1
      integer, parameter :: p_rsp_Lc_div_L = p_rsp_Lr_div_L + 1
      integer, parameter :: p_rsp_Lt_div_L = p_rsp_Lc_div_L + 1

      integer, parameter :: p_rsp_WORK = p_rsp_Lt_div_L + 1
      integer, parameter :: p_rsp_WORKQ = p_rsp_WORK + 1
      integer, parameter :: p_rsp_WORKT = p_rsp_WORKQ + 1
      integer, parameter :: p_rsp_WORKC = p_rsp_WORKT + 1

      integer, parameter :: p_d_u_div_rmid_start = p_rsp_WORKC + 1
      integer, parameter :: p_d_u_div_rmid = p_d_u_div_rmid_start + 1

      integer, parameter :: p_conv_vel_residual = p_d_u_div_rmid + 1
      integer, parameter :: p_log_conv_vel_residual = p_conv_vel_residual + 1
      integer, parameter :: p_dconv_vel_dt = p_log_conv_vel_residual + 1

      integer, parameter :: p_cell_ie_div_star_ie = p_dconv_vel_dt + 1
      integer, parameter :: p_log_cell_specific_IE = p_cell_ie_div_star_ie + 1
      integer, parameter :: p_log_cell_ie_div_star_ie = p_log_cell_specific_IE + 1

      integer, parameter :: p_cell_specific_PE = p_log_cell_ie_div_star_ie + 1
      integer, parameter :: p_cell_specific_IE = p_cell_specific_PE + 1
      integer, parameter :: p_cell_specific_KE = p_cell_specific_IE + 1
      integer, parameter :: p_cell_IE_div_IE_plus_KE = p_cell_specific_KE + 1
      integer, parameter :: p_cell_KE_div_IE_plus_KE = p_cell_IE_div_IE_plus_KE + 1

      integer, parameter :: p_grada_sub_gradr = p_cell_KE_div_IE_plus_KE + 1
      integer, parameter :: p_gradL_sub_gradr = p_grada_sub_gradr + 1
      integer, parameter :: p_gradL = p_gradL_sub_gradr + 1
      integer, parameter :: p_sch_stable = p_gradL + 1
      integer, parameter :: p_ledoux_stable = p_sch_stable + 1
      integer, parameter :: p_grad_density = p_ledoux_stable + 1
      integer, parameter :: p_grad_temperature = p_grad_density + 1

      integer, parameter :: p_dominant_isoA_for_thermohaline = p_grad_temperature + 1

      integer, parameter :: p_dominant_isoZ_for_thermohaline = p_dominant_isoA_for_thermohaline + 1
      integer, parameter :: p_gradL_composition_term = p_dominant_isoZ_for_thermohaline + 1
      integer, parameter :: p_log_brunt_nonB = p_gradL_composition_term + 1

      integer, parameter :: p_log_brunt_B = p_log_brunt_nonB + 1
      integer, parameter :: p_brunt_nonB = p_log_brunt_B + 1
      integer, parameter :: p_brunt_B = p_brunt_nonB + 1
      integer, parameter :: p_lamb_S2 = p_brunt_B + 1
      integer, parameter :: p_lamb_S = p_lamb_S2 + 1
      integer, parameter :: p_lamb_Sl1 = p_lamb_S + 1
      integer, parameter :: p_lamb_Sl2 = p_lamb_Sl1 + 1
      integer, parameter :: p_lamb_Sl3 = p_lamb_Sl2 + 1
      integer, parameter :: p_lamb_Sl10 = p_lamb_Sl3 + 1
      integer, parameter :: p_sign_brunt_N2 = p_lamb_Sl10 + 1
      integer, parameter :: p_brunt_A_div_x2 = p_sign_brunt_N2 + 1
      integer, parameter :: p_brunt_A = p_brunt_A_div_x2 + 1

      integer, parameter :: p_brunt_N2 = p_brunt_A + 1
      integer, parameter :: p_brunt_N2_structure_term = p_brunt_N2 + 1
      integer, parameter :: p_brunt_N2_composition_term = p_brunt_N2_structure_term + 1
      integer, parameter :: p_log_brunt_N2_structure_term = p_brunt_N2_composition_term + 1
      integer, parameter :: p_log_brunt_N2_composition_term = p_log_brunt_N2_structure_term + 1
      integer, parameter :: p_log_brunt_N2_dimensionless = p_log_brunt_N2_composition_term + 1
      integer, parameter :: p_brunt_N2_dimensionless = p_log_brunt_N2_dimensionless + 1
      integer, parameter :: p_brunt_N_dimensionless = p_brunt_N2_dimensionless + 1
      integer, parameter :: p_brunt_N = p_brunt_N_dimensionless + 1
      integer, parameter :: p_brunt_frequency = p_brunt_N + 1
      integer, parameter :: p_log_brunt_N2 = p_brunt_frequency + 1
      integer, parameter :: p_log_brunt_N = p_log_brunt_N2 + 1
      integer, parameter :: p_brunt_N_div_r_integral = p_log_brunt_N + 1
      integer, parameter :: p_brunt_nu = p_brunt_N_div_r_integral + 1
      integer, parameter :: p_log_brunt_nu = p_brunt_nu + 1
      integer, parameter :: p_log_lamb_Sl1 = p_log_brunt_nu + 1
      integer, parameter :: p_log_lamb_Sl2 = p_log_lamb_Sl1 + 1
      integer, parameter :: p_log_lamb_Sl3 = p_log_lamb_Sl2 + 1
      integer, parameter :: p_log_lamb_Sl10 = p_log_lamb_Sl3 + 1
      integer, parameter :: p_brunt_N2_sub_omega2 = p_log_lamb_Sl10 + 1
      integer, parameter :: p_sl2_sub_omega2 = p_brunt_N2_sub_omega2 + 1
      integer, parameter :: p_k_r_integral = p_sl2_sub_omega2 + 1

      integer, parameter :: p_num_steps = p_k_r_integral + 1
      integer, parameter :: p_mtx_solve = p_num_steps + 1
      integer, parameter :: p_mtx_factor = p_mtx_solve + 1
      
      integer, parameter :: p_max_abs_xa_corr = p_mtx_factor + 1
      integer, parameter :: p_log_zFe = p_max_abs_xa_corr + 1
      integer, parameter :: p_zFe = p_log_zFe + 1
      integer, parameter :: p_log_u_residual = p_zFe + 1
      integer, parameter :: p_u_residual = p_log_u_residual + 1
      integer, parameter :: p_u = p_u_residual + 1
      integer, parameter :: p_u_face = p_u + 1
      integer, parameter :: p_dPdr_dRhodr_info = p_u_face + 1
      integer, parameter :: p_signed_log_ergs_err = p_dPdr_dRhodr_info + 1
      integer, parameter :: p_RTI_du_diffusion_kick = p_signed_log_ergs_err + 1
      integer, parameter :: p_log_du_kick_div_du = p_RTI_du_diffusion_kick + 1
      
      integer, parameter :: p_col_id_max = p_log_du_kick_div_du

      character (len=maxlen_profile_column_name) :: profile_column_name(p_col_id_max)
      type (integer_dict), pointer :: profile_column_names_dict


      contains


      subroutine profile_column_names_init(ierr)
         use utils_lib, only: integer_dict_define
         integer, intent(out) :: ierr
         integer :: i, cnt

         ierr = 0

         cnt = 0
         profile_column_name(:) = ''

         profile_column_name(p_zone) = 'zone'
         profile_column_name(p_k) = 'k'
         profile_column_name(p_luminosity) = 'luminosity'
         profile_column_name(p_lum_erg_s) = 'lum_erg_s'
         profile_column_name(p_log_abs_lum_erg_s) = 'log_abs_lum_erg_s'
         profile_column_name(p_log_Lrad_div_Ledd) = 'log_Lrad_div_Ledd'
         profile_column_name(p_log_Lrad_div_L) = 'log_Lrad_div_L'
         profile_column_name(p_lum_plus_lum_adv) = 'lum_plus_lum_adv'
         profile_column_name(p_lum_adv) = 'lum_adv'
         profile_column_name(p_lum_rad) = 'lum_rad'
         profile_column_name(p_lum_conv) = 'lum_conv'
         profile_column_name(p_lum_conv_MLT) = 'lum_conv_MLT'

         profile_column_name(p_lum_conv_div_lum_Edd) = 'p_lum_conv_MLT'
         profile_column_name(p_lum_rad_div_L_Edd) = 'lum_rad_div_L_Edd'
         profile_column_name(p_lum_rad_div_L) = 'lum_rad_div_L'
         profile_column_name(p_lum_conv_div_L) = 'lum_conv_div_L'
         profile_column_name(p_lum_conv_div_lum_rad) = 'lum_conv_div_lum_rad'
         profile_column_name(p_log_Lrad) = 'log_Lrad'
         profile_column_name(p_log_Lconv) = 'log_Lconv'
         profile_column_name(p_log_Lconv_div_L) = 'log_Lconv_div_L'
         profile_column_name(p_lum_rad_div_L) = 'lum_rad_div_L'

         profile_column_name(p_grav) = 'grav'
         profile_column_name(p_log_g) = 'log_g'
         profile_column_name(p_r_div_g) = 'r_div_g'
         profile_column_name(p_g_div_r) = 'g_div_r'
         profile_column_name(p_net_nuclear_energy) = 'net_nuclear_energy'
         
         profile_column_name(p_eps_nuc_plus_nuc_neu) = 'eps_nuc_plus_nuc_neu'
         profile_column_name(p_eps_nuc_minus_non_nuc_neu) = 'eps_nuc_minus_non_nuc_neu'
         profile_column_name(p_net_energy) = 'net_energy'
         profile_column_name(p_logL) = 'logL'
         profile_column_name(p_log_Ledd) = 'log_Ledd'
         profile_column_name(p_lum_div_Ledd) = 'lum_div_Ledd'
         profile_column_name(p_log_L_div_Ledd) = 'log_L_div_Ledd'
         profile_column_name(p_signed_log_power) = 'signed_log_power'
         profile_column_name(p_log_abs_dvdt_div_v) = 'log_abs_dvdt_div_v'
         profile_column_name(p_log_abs_v) = 'log_abs_v'
         profile_column_name(p_log_diff_grads) = 'log_diff_grads'
         profile_column_name(p_diff_grads) = 'diff_grads'
         profile_column_name(p_gradT_excess_effect) = 'gradT_excess_effect'
         profile_column_name(p_superad_reduction_factor) = 'superad_red'

         profile_column_name(p_v) = 'v'
         profile_column_name(p_velocity) = 'velocity'
         profile_column_name(p_vel_km_per_s) = 'vel_km_per_s'
         profile_column_name(p_radius_km) = 'radius_km'
         profile_column_name(p_radius_cm) = 'radius_cm'
         profile_column_name(p_radius) = 'radius'
         profile_column_name(p_rmid) = 'rmid'
         profile_column_name(p_logR_cm) = 'logR_cm'
         profile_column_name(p_logR) = 'logR'
         profile_column_name(p_log_q) = 'log_q'
         profile_column_name(p_q) = 'q'
         profile_column_name(p_dq) = 'dq'
         profile_column_name(p_log_dq) = 'log_dq'
         profile_column_name(p_logtau_sub_xlogtau) = 'logtau_sub_xlogtau'
         profile_column_name(p_xlogtau) = 'xlogtau'
         profile_column_name(p_logtau) = 'logtau'
         profile_column_name(p_pgas_div_p) = 'pgas_div_p'
         profile_column_name(p_pgas_div_ptotal) = 'pgas_div_ptotal'

         profile_column_name(p_prad_div_pgas) = 'prad_div_pgas'
         profile_column_name(p_prad_div_pgas_div_L_div_Ledd) = 'prad_div_pgas_div_L_div_Ledd'

         profile_column_name(p_m_div_r) = 'm_div_r'
         profile_column_name(p_dmbar_m_div_r) = 'dmbar_m_div_r'
         profile_column_name(p_log_dmbar_m_div_r) = 'log_dmbar_m_div_r'

         profile_column_name(p_log_mass) = 'log_mass'
         profile_column_name(p_mass) = 'mass'
         profile_column_name(p_mass_grams) = 'mass_grams'
         profile_column_name(p_mmid) = 'mmid'

         profile_column_name(p_dm) = 'dm'
         profile_column_name(p_dm_bar) = 'dm_bar'

         profile_column_name(p_xr) = 'xr'
         profile_column_name(p_xr_cm) = 'xr_cm'
         profile_column_name(p_xr_div_R) = 'xr_div_R'
         profile_column_name(p_log_xr) = 'log_xr'
         profile_column_name(p_log_xr_cm) = 'log_xr_cm'
         profile_column_name(p_log_xr_div_R) = 'log_xr_div_R'

         profile_column_name(p_m_grav) = 'm_grav'
         profile_column_name(p_mass_correction_factor) = 'mass_correction_factor'
         profile_column_name(p_m_grav_div_m_baryonic) = 'm_grav_div_m_baryonic'

         profile_column_name(p_avg_charge_H) = 'avg_charge_H'
         profile_column_name(p_avg_charge_He) = 'avg_charge_He'
         profile_column_name(p_avg_charge_C) = 'avg_charge_C'
         profile_column_name(p_avg_charge_N) = 'avg_charge_N'
         profile_column_name(p_avg_charge_O) = 'avg_charge_O'
         profile_column_name(p_avg_charge_Ne) = 'avg_charge_Ne'
         profile_column_name(p_avg_charge_Mg) = 'avg_charge_Mg'
         profile_column_name(p_avg_charge_Si) = 'avg_charge_Si'
         profile_column_name(p_avg_charge_Fe) = 'avg_charge_Fe'
         profile_column_name(p_neutral_fraction_H) = 'neutral_fraction_H'
         profile_column_name(p_neutral_fraction_He) = 'neutral_fraction_He'
         profile_column_name(p_neutral_fraction_C) = 'neutral_fraction_C'
         profile_column_name(p_neutral_fraction_N) = 'neutral_fraction_N'
         profile_column_name(p_neutral_fraction_O) = 'neutral_fraction_O'
         profile_column_name(p_neutral_fraction_Ne) = 'neutral_fraction_Ne'
         profile_column_name(p_neutral_fraction_Mg) = 'neutral_fraction_Mg'
         profile_column_name(p_neutral_fraction_Si) = 'neutral_fraction_Si'
         profile_column_name(p_neutral_fraction_Fe) = 'neutral_fraction_Fe'

         profile_column_name(p_log_x) = 'log_x'
         profile_column_name(p_x) = 'x'
         profile_column_name(p_log_y) = 'log_y'
         profile_column_name(p_y) = 'y'
         profile_column_name(p_log_z) = 'log_z'
         profile_column_name(p_z) = 'z'

         profile_column_name(p_xm) = 'xm'
         profile_column_name(p_xq) = 'xq'
         profile_column_name(p_logxm) = 'logxm'
         profile_column_name(p_logxq) = 'logxq'
         profile_column_name(p_logdq) = 'logdq'

         profile_column_name(p_log_radial_depth) = 'log_radial_depth'
         profile_column_name(p_log_column_depth) = 'log_column_depth'
         profile_column_name(p_dRstar_div_dr) = 'dRstar_div_dr'
         profile_column_name(p_dr_div_R) = 'dr_div_R'
         profile_column_name(p_log_dr_div_R) = 'log_dr_div_R'
         profile_column_name(p_r_div_R) = 'r_div_R'
         profile_column_name(p_dr_div_rmid) = 'dr_div_rmid'
         profile_column_name(p_log_dr_div_rmid) = 'log_dr_div_rmid'
         profile_column_name(p_log_dr) = 'log_dr'
         profile_column_name(p_dlogR) = 'dlogR'

         profile_column_name(p_lnR_residual) = 'lnR_residual'
         profile_column_name(p_lnd_residual) = 'lnd_residual'

         profile_column_name(p_cell_internal_energy_fraction_start) = 'cell_internal_energy_fraction_start'
         profile_column_name(p_cell_internal_energy_fraction) = 'cell_internal_energy_fraction'
         profile_column_name(p_ergs_error) = 'ergs_error'
         profile_column_name(p_log_rel_E_err) = 'log_rel_E_err'
         profile_column_name(p_ergs_error_integral) = 'ergs_error_integral'
         profile_column_name(p_ergs_rel_error_integral) = 'ergs_rel_error_integral'
         profile_column_name(p_E_residual) = 'E_residual'
         profile_column_name(p_equL_residual) = 'equL_residual'
         profile_column_name(p_dvdt_residual) = 'dvdt_residual'
         profile_column_name(p_v_residual) = 'v_residual'

         profile_column_name(p_Et_residual) = 'Et_residual'
         profile_column_name(p_log_Et_residual) = 'log_Et_residual'
         profile_column_name(p_log_E_residual) = 'log_E_residual'
         profile_column_name(p_log_lnR_residual) = 'log_lnR_residual'
         profile_column_name(p_log_lnd_residual) = 'log_lnd_residual'
         profile_column_name(p_log_equL_residual) = 'log_equL_residual'
         profile_column_name(p_log_dvdt_residual) = 'log_dvdt_residual'
         profile_column_name(p_log_v_residual) = 'log_v_residual'

         profile_column_name(p_t_rad) = 't_rad'
         profile_column_name(p_log_t_rad) = 'log_t_rad'

         profile_column_name(p_log_dt_cs_div_dr) = 'log_dt_cs_div_dr'
         profile_column_name(p_dt_cs_div_dr) = 'dt_cs_div_dr'
         profile_column_name(p_dr_div_cs) = 'dr_div_cs'
         profile_column_name(p_log_dr_div_cs) = 'log_dr_div_cs'
         profile_column_name(p_dr_div_cs_yr) = 'dr_div_cs_yr'
         profile_column_name(p_log_dr_div_cs_yr) = 'log_dr_div_cs_yr'

         profile_column_name(p_cell_collapse_time) = 'cell_collapse_time'
         profile_column_name(p_log_cell_collapse_time) = 'log_cell_collapse_time'

         profile_column_name(p_log_acoustic_depth) = 'log_acoustic_depth'
         profile_column_name(p_log_acoustic_radius) = 'log_acoustic_radius'
         profile_column_name(p_acoustic_depth) = 'acoustic_depth'
         profile_column_name(p_acoustic_radius) = 'acoustic_radius'
         profile_column_name(p_acoustic_r_div_R_phot) = 'acoustic_r_div_R_phot'

         profile_column_name(p_compression_gradient) = 'compression_gradient'
         profile_column_name(p_dq_ratio) = 'dq_ratio'
         profile_column_name(p_tau_eff) = 'tau_eff'
         profile_column_name(p_tau_eff_div_tau) = 'tau_eff_div_tau'
         profile_column_name(p_xtau) = 'xtau'
         profile_column_name(p_tau) = 'tau'
         profile_column_name(p_extra_opacity_factor) = 'extra_opacity_factor'
         profile_column_name(p_log_kap_times_factor) = 'log_kap_times_factor'
         profile_column_name(p_log_opacity) = 'log_opacity'
         profile_column_name(p_kap_frac_Type2) = 'kap_frac_Type2'
         profile_column_name(p_kap_frac_op_mono) = 'kap_frac_op_mono'
         profile_column_name(p_energy) = 'energy'
         profile_column_name(p_logM) = 'logM'
         profile_column_name(p_temperature) = 'temperature'
         profile_column_name(p_logT_face) = 'logT_face'
         profile_column_name(p_logT_bb) = 'logT_bb'
         profile_column_name(p_logT_face_div_logT_bb) = 'logT_face_div_logT_bb'
         profile_column_name(p_logT) = 'logT'


         profile_column_name(p_rho) = 'rho'
         profile_column_name(p_density) = 'density'
         profile_column_name(p_logRho) = 'logRho'
         profile_column_name(p_pgas) = 'pgas'
         profile_column_name(p_logPgas) = 'logPgas'
         profile_column_name(p_prad) = 'prad'
         profile_column_name(p_pressure) = 'pressure'
         profile_column_name(p_logP) = 'logP'
         profile_column_name(p_logE) = 'logE'
         profile_column_name(p_grada) = 'grada'
         profile_column_name(p_dE_dRho) = 'dE_dRho'
         profile_column_name(p_cv) = 'cv'
         profile_column_name(p_thermal_time_to_surface) = 'thermal_time_to_surface'
         profile_column_name(p_log_thermal_time_to_surface) = 'log_thermal_time_to_surface'
         profile_column_name(p_cp) = 'cp'
         profile_column_name(p_log_CpT) = 'log_CpT'
         profile_column_name(p_log_CpT_absMdot_div_L) = 'log_CpT_absMdot_div_L'
         profile_column_name(p_logS) = 'logS'
         profile_column_name(p_logS_per_baryon) = 'logS_per_baryon'
         profile_column_name(p_gamma1) = 'gamma1'
         profile_column_name(p_gamma3) = 'gamma3'
         profile_column_name(p_eta) = 'eta'
         profile_column_name(p_theta_e) = 'theta_e'
         profile_column_name(p_gam) = 'gam'
         profile_column_name(p_mu) = 'mu'

         profile_column_name(p_eos_frac_OPAL_SCVH) = 'eos_frac_OPAL_SCVH'
         profile_column_name(p_eos_frac_HELM) = 'eos_frac_HELM'
         profile_column_name(p_eos_frac_Skye) = 'eos_frac_Skye'
         profile_column_name(p_eos_frac_PC) = 'eos_frac_PC'
         profile_column_name(p_eos_frac_FreeEOS) = 'eos_frac_FreeEOS'
         profile_column_name(p_eos_frac_CMS) = 'eos_frac_CMS'

         profile_column_name(p_log_rho_times_r3) = 'log_rho_times_r3'
         profile_column_name(p_rho_times_r3) = 'rho_times_r3'
         profile_column_name(p_v_times_t_div_r) = 'v_times_t_div_r'
         profile_column_name(p_v_div_r) = 'v_div_r'
         profile_column_name(p_d_v_div_r) = 'd_v_div_r'
         
         profile_column_name(p_log_c_div_tau) = 'log_c_div_tau'
         profile_column_name(p_log_v_escape) = 'log_v_escape'
         profile_column_name(p_v_escape) = 'v_escape'
         profile_column_name(p_v_div_v_escape) = 'v_div_v_escape'
         profile_column_name(p_v_div_vesc) = 'v_div_vesc'
         profile_column_name(p_v_div_cs) = 'v_div_cs'
         profile_column_name(p_v_div_csound) = 'v_div_csound'
         profile_column_name(p_log_csound) = 'log_csound'
         profile_column_name(p_csound) = 'csound'
         profile_column_name(p_csound_face) = 'csound_face'

         profile_column_name(p_omega) = 'omega'
         profile_column_name(p_log_omega) = 'log_omega'
         profile_column_name(p_log_j_rot) = 'log_j_rot'
         profile_column_name(p_log_J_div_M53) = 'log_J_div_M53'
         profile_column_name(p_log_J_inside) = 'log_J_inside'
         profile_column_name(p_shear) = 'shear'
         profile_column_name(p_log_abs_shear) = 'log_abs_shear'
         profile_column_name(p_log_abs_dlnR_domega) = 'log_abs_dlnR_domega'
         profile_column_name(p_i_rot) = 'i_rot'
         profile_column_name(p_j_rot) = 'j_rot'
         profile_column_name(p_v_rot) = 'v_rot'
         profile_column_name(p_fp_rot) = 'fp_rot'
         profile_column_name(p_ft_rot) = 'ft_rot'
         profile_column_name(p_ft_rot_div_fp_rot) = 'ft_rot_div_fp_rot'
         profile_column_name(p_w_div_w_crit_roche) = 'w_div_w_crit_roche'
         profile_column_name(p_w_div_w_crit_roche2) = 'w_div_w_crit_roche2'
         profile_column_name(p_log_am_nu_rot) = 'log_am_nu_rot'
         profile_column_name(p_log_am_nu_non_rot) = 'log_am_nu_non_rot'
         profile_column_name(p_log_am_nu) = 'log_am_nu'

         profile_column_name(p_eps_WD_sedimentation) = 'eps_WD_sedimentation'
         profile_column_name(p_log_eps_WD_sedimentation) = 'log_eps_WD_sedimentation'

         profile_column_name(p_eps_diffusion) = 'eps_diffusion'
         profile_column_name(p_log_eps_diffusion) = 'log_eps_diffusion'
         
         profile_column_name(p_log_e_field) = 'log_e_field'
         profile_column_name(p_e_field) = 'e_field'
         profile_column_name(p_log_g_field_element_diffusion) = 'log_g_field_element_diffusion'
         profile_column_name(p_g_field_element_diffusion) = 'g_field_element_diffusion'
         profile_column_name(p_log_eE_div_mg_element_diffusion) = 'log_eE_div_mg_element_diffusion'
         profile_column_name(p_eE_div_mg_element_diffusion) = 'eE_div_mg_element_diffusion'

         profile_column_name(p_r_polar) = 'r_polar'
         profile_column_name(p_log_r_polar) = 'log_r_polar'
         profile_column_name(p_r_equatorial) = 'r_equatorial'
         profile_column_name(p_log_r_equatorial) = 'log_r_equatorial'
         profile_column_name(p_r_e_div_r_p) = 'r_e_div_r_p'
         profile_column_name(p_omega_crit) = 'omega_crit'
         profile_column_name(p_omega_div_omega_crit) = 'omega_div_omega_crit'

         profile_column_name(p_am_log_sig) = 'am_log_sig'
         profile_column_name(p_am_log_sig_omega) = 'am_log_sig_omega'
         profile_column_name(p_am_log_sig_j) = 'am_log_sig_j'

         profile_column_name(p_am_log_nu_omega) = 'am_log_nu_omega'
         profile_column_name(p_am_log_nu_j) = 'am_log_nu_j'
         profile_column_name(p_am_log_nu_rot) = 'am_log_nu_rot'
         profile_column_name(p_am_log_nu_non_rot) = 'am_log_nu_non_rot'

         profile_column_name(p_am_domega_dlnR) = 'am_domega_dlnR'
         profile_column_name(p_richardson_number) = 'richardson_number'
         profile_column_name(p_am_log_D_visc) = 'am_log_D_visc'
         profile_column_name(p_am_log_D_DSI) = 'am_log_D_DSI'
         profile_column_name(p_am_log_D_SH) = 'am_log_D_SH'
         profile_column_name(p_am_log_D_SSI) = 'am_log_D_SSI'
         profile_column_name(p_am_log_D_ES) = 'am_log_D_ES'
         profile_column_name(p_am_log_D_GSF) = 'am_log_D_GSF'
         profile_column_name(p_am_log_D_ST) = 'am_log_D_ST'
         profile_column_name(p_am_log_nu_ST) = 'am_log_nu_ST'

         profile_column_name(p_dynamo_log_B_r) = 'dynamo_log_B_r'
         profile_column_name(p_dynamo_log_B_phi) = 'dynamo_log_B_phi'

         profile_column_name(p_grada_face) = 'grada_face'
         profile_column_name(p_gradr_div_grada) = 'gradr_div_grada'
         profile_column_name(p_gradr_sub_grada) = 'gradr_sub_grada'
         profile_column_name(p_scale_height) = 'scale_height'

         profile_column_name(p_entropy) = 'entropy'
         profile_column_name(p_free_e) = 'free_e'
         profile_column_name(p_logfree_e) = 'logfree_e'
         profile_column_name(p_chiRho) = 'chiRho'
         profile_column_name(p_chiT) = 'chiT'
         profile_column_name(p_QQ) = 'QQ'

         profile_column_name(p_phase) = 'eos_phase'
         profile_column_name(p_latent_ddlnT) = 'latent_ddlnT'
         profile_column_name(p_latent_ddlnRho) = 'latent_ddlnRho'

         profile_column_name(p_chiRho_for_partials) = 'chiRho_for_partials'
         profile_column_name(p_chiT_for_partials) = 'chiT_for_partials'
         profile_column_name(p_rel_diff_chiRho_for_partials) = 'rel_diff_chiRho_for_partials'
         profile_column_name(p_rel_diff_chiT_for_partials) = 'rel_diff_chiT_for_partials'

         profile_column_name(p_x_mass_fraction_H) = 'x_mass_fraction_H'
         profile_column_name(p_y_mass_fraction_He) = 'y_mass_fraction_He'
         profile_column_name(p_z_mass_fraction_metals) = 'z_mass_fraction_metals'

         profile_column_name(p_abar) = 'abar'
         profile_column_name(p_zbar) = 'zbar'
         profile_column_name(p_z2bar) = 'z2bar'
         profile_column_name(p_ye) = 'ye'
         profile_column_name(p_opacity) = 'opacity'
         profile_column_name(p_dkap_dlnrho_face) = 'dkap_dlnrho_face'
         profile_column_name(p_dkap_dlnT_face) = 'dkap_dlnT_face'

         profile_column_name(p_eps_nuc_start) = 'eps_nuc_start'
         profile_column_name(p_eps_nuc) = 'eps_nuc'
         profile_column_name(p_signed_log_eps_nuc) = 'signed_log_eps_nuc'
         profile_column_name(p_log_abs_eps_nuc) = 'log_abs_eps_nuc'
         profile_column_name(p_burn_sum_xa_err_before_fix) = 'burn_sum_xa_err_before_fix'
         profile_column_name(p_burn_log_abs_sum_xa_err_before_fix) = 'burn_log_abs_sum_xa_err_before_fix'

         profile_column_name(p_d_epsnuc_dlnd) = 'd_epsnuc_dlnd'
         profile_column_name(p_d_epsnuc_dlnT) = 'd_epsnuc_dlnT'
         profile_column_name(p_d_lnepsnuc_dlnd) = 'd_lnepsnuc_dlnd'
         profile_column_name(p_d_lnepsnuc_dlnT) = 'd_lnepsnuc_dlnT'
         profile_column_name(p_deps_dlnd_face) = 'deps_dlnd_face'
         profile_column_name(p_deps_dlnT_face) = 'deps_dlnT_face'
         profile_column_name(p_eps_nuc_neu_total) = 'eps_nuc_neu_total'

         profile_column_name(p_non_nuc_neu) = 'non_nuc_neu'
         profile_column_name(p_nonnucneu_plas) = 'nonnucneu_plas'
         profile_column_name(p_nonnucneu_brem) = 'nonnucneu_brem'
         profile_column_name(p_nonnucneu_phot) = 'nonnucneu_phot'
         profile_column_name(p_nonnucneu_pair) = 'nonnucneu_pair'
         profile_column_name(p_nonnucneu_reco) = 'nonnucneu_reco'

         profile_column_name(p_log_irradiation_heat) = 'log_irradiation_heat'
         profile_column_name(p_extra_L) = 'extra_L'
         profile_column_name(p_log_extra_L) = 'log_extra_L'
         profile_column_name(p_extra_jdot) = 'extra_jdot'
         profile_column_name(p_extra_omegadot) = 'extra_omegadot'
         profile_column_name(p_extra_grav) = 'extra_grav'
         profile_column_name(p_extra_heat) = 'extra_heat'
         profile_column_name(p_alpha_mlt) = 'alpha_mlt'
         profile_column_name(p_cgrav_factor) = 'cgrav_factor'
         profile_column_name(p_div_v) = 'div_v'

         profile_column_name(p_d_v_div_r_dm) = 'd_v_div_r_dm'
         profile_column_name(p_d_v_div_r_dr) = 'd_v_div_r_dr'
         profile_column_name(p_dvdt_grav) = 'dvdt_grav'
         profile_column_name(p_dvdt_dPdm) = 'dvdt_dPdm'
         profile_column_name(p_du) = 'du'
         profile_column_name(p_P_face) = 'P_face'
         profile_column_name(p_log_P_face) = 'log_P_face'

         profile_column_name(p_hse_ratio) = 'hse_ratio'
         profile_column_name(p_hse_ratio_gyre) = 'hse_ratio_gyre'

         profile_column_name(p_dlnP_dlnR) = 'dlnP_dlnR'
         profile_column_name(p_dlnRho_dlnR) = 'dlnRho_dlnR'
         profile_column_name(p_gradP_div_rho) = 'gradP_div_rho'
         profile_column_name(p_dPdr_div_grav) = 'dPdr_div_grav'

         profile_column_name(p_dvdt_RTI_diffusion) = 'dvdt_RTI_diffusion'
         profile_column_name(p_dlnddt_RTI_diffusion) = 'dlnddt_RTI_diffusion'
         profile_column_name(p_log_abs_eps_grav_dm_div_L) = 'log_abs_eps_grav_dm_div_L'
         profile_column_name(p_eps_grav_composition_term) = 'eps_grav_composition_term'

         profile_column_name(p_dt_dm_eps_grav) = 'dt_dm_eps_grav'
         profile_column_name(p_dm_eps_grav) = 'dm_eps_grav'
         profile_column_name(p_dm_de) = 'dm_de'
         profile_column_name(p_dt_dL) = 'dt_dL'

         profile_column_name(p_eps_grav_plus_eps_mdot) = 'eps_grav_plus_eps_mdot'
         profile_column_name(p_ergs_eps_grav_plus_eps_mdot) = 'ergs_eps_grav_plus_eps_mdot'
         profile_column_name(p_ergs_mdot) = 'ergs_mdot'
         profile_column_name(p_eps_mdot) = 'eps_mdot'
         
         profile_column_name(p_log_xm_div_delta_m) = 'log_xm_div_delta_m'
         profile_column_name(p_xm_div_delta_m) = 'xm_div_delta_m'
         profile_column_name(p_eps_grav) = 'eps_grav'
         profile_column_name(p_env_eps_grav) = 'env_eps_grav'
         profile_column_name(p_signed_log_eps_grav) = 'signed_log_eps_grav'
         profile_column_name(p_mlt_mixing_length) = 'mlt_mixing_length'
         profile_column_name(p_log_conv_L_div_L) = 'log_conv_L_div_L'
         profile_column_name(p_conv_L_div_L) = 'conv_L_div_L'
         profile_column_name(p_mlt_Zeta) = 'mlt_Zeta'
         profile_column_name(p_mlt_Gamma) = 'mlt_Gamma'
         profile_column_name(p_mlt_Pturb) = 'mlt_Pturb'
         profile_column_name(p_mlt_mixing_type) = 'mlt_mixing_type'

         profile_column_name(p_grada_sub_gradT) = 'grada_sub_gradT'
         profile_column_name(p_gradT_sub_grada) = 'gradT_sub_grada'
         profile_column_name(p_grad_superad) = 'grad_superad'
         profile_column_name(p_grad_superad_actual) = 'grad_superad_actual'

         profile_column_name(p_gradT_sub_a) = 'gradT_sub_a'
         profile_column_name(p_gradT_div_grada) = 'gradT_div_grada'

         profile_column_name(p_gradT_rel_err) = 'gradT_rel_err'
         profile_column_name(p_gradr_sub_gradT) = 'gradr_sub_gradT'
         profile_column_name(p_gradT_sub_gradr) = 'gradT_sub_gradr'
         profile_column_name(p_gradT_div_gradr) = 'gradT_div_gradr'
         profile_column_name(p_log_gradT_div_gradr) = 'log_gradT_div_gradr'

         profile_column_name(p_log_mlt_Gamma) = 'log_mlt_Gamma'
         profile_column_name(p_log_mlt_vc) = 'log_mlt_vc'
         profile_column_name(p_conv_vel_div_mlt_vc) = 'conv_vel_div_mlt_vc'
         profile_column_name(p_mlt_vc) = 'mlt_vc'
         profile_column_name(p_super_ad) = 'super_ad'

         profile_column_name(p_delta_r) = 'delta_r'
         profile_column_name(p_delta_L) = 'delta_L'
         profile_column_name(p_delta_cell_vol) = 'delta_cell_vol'
         profile_column_name(p_delta_entropy) = 'delta_entropy'
         profile_column_name(p_delta_T) = 'delta_T'
         profile_column_name(p_delta_rho) = 'delta_rho'
         profile_column_name(p_delta_eps_nuc) = 'delta_eps_nuc'
         profile_column_name(p_delta_mu) = 'delta_mu'
         profile_column_name(p_log_D_conv) = 'log_D_conv'
         profile_column_name(p_log_D_smooth) = 'log_D_smooth'
         profile_column_name(p_log_D_leftover) = 'log_D_leftover'
         profile_column_name(p_log_D_semi) = 'log_D_semi'
         profile_column_name(p_log_D_ovr) = 'log_D_ovr'
         profile_column_name(p_log_D_anon) = 'log_D_anon'
         profile_column_name(p_log_D_thrm) = 'log_D_thrm'
         profile_column_name(p_log_D_rayleigh_taylor) = 'log_D_rayleigh_taylor'
         profile_column_name(p_log_D_minimum) = 'log_D_minimum'
         profile_column_name(p_log_D_mix_non_rotation) = 'log_D_mix_non_rotation'
         profile_column_name(p_log_D_mix_rotation) = 'log_D_mix_rotation'
         profile_column_name(p_log_D_omega) = 'log_D_omega'

         profile_column_name(p_log_lambda_RTI_div_Hrho) = 'log_lambda_RTI_div_Hrho'
         profile_column_name(p_lambda_RTI) = 'lambda_RTI'
         profile_column_name(p_dPdr_info) = 'dPdr_info'
         profile_column_name(p_dRhodr_info) = 'dRhodr_info'
         
         profile_column_name(p_source_plus_alpha_RTI) = 'source_plus_alpha_RTI'
         profile_column_name(p_source_minus_alpha_RTI) = 'source_minus_alpha_RTI'
         profile_column_name(p_log_source_RTI) = 'log_source_RTI'
         profile_column_name(p_log_source_plus_alpha_RTI) = 'log_source_plus_alpha_RTI'
         profile_column_name(p_log_source_minus_alpha_RTI) = 'log_source_minus_alpha_RTI'

         profile_column_name(p_dudt_RTI) = 'dudt_RTI'
         profile_column_name(p_dedt_RTI) = 'dedt_RTI'

         profile_column_name(p_eta_RTI) = 'eta_RTI'
         profile_column_name(p_log_eta_RTI) = 'log_eta_RTI'
         profile_column_name(p_boost_for_eta_RTI) = 'boost_for_eta_RTI'
         profile_column_name(p_log_boost_for_eta_RTI) = 'log_boost_for_eta_RTI'

         profile_column_name(p_alpha_RTI) = 'alpha_RTI'
         profile_column_name(p_log_alpha_RTI) = 'log_alpha_RTI'
         profile_column_name(p_log_etamid_RTI) = 'log_etamid_RTI'
         profile_column_name(p_log_sig_RTI) = 'log_sig_RTI'
         profile_column_name(p_log_sigmid_RTI) = 'log_sigmid_RTI'

         profile_column_name(p_log_D_mix) = 'log_D_mix'
         profile_column_name(p_log_sig_raw_mix) = 'log_sig_raw_mix'
         profile_column_name(p_log_sig_mix) = 'log_sig_mix'

         profile_column_name(p_burn_avg_epsnuc) = 'burn_avg_epsnuc'
         profile_column_name(p_log_burn_avg_epsnuc) = 'log_burn_avg_epsnuc'
         profile_column_name(p_burn_num_iters) = 'burn_num_iters'

         profile_column_name(p_dt_times_conv_vel_div_mixing_length) = 'dt_times_conv_vel_div_mixing_length'
         profile_column_name(p_log_dt_times_conv_vel_div_mixing_length) = 'log_dt_times_conv_vel_div_mixing_length'
         profile_column_name(p_conv_vel) = 'conv_vel'
         profile_column_name(p_log_conv_vel) = 'log_conv_vel'
         profile_column_name(p_conv_vel_div_csound) = 'conv_vel_div_csound'
         profile_column_name(p_conv_vel_div_L_vel) = 'conv_vel_div_L_vel'
         profile_column_name(p_conv_mixing_type) = 'conv_mixing_type'
         profile_column_name(p_log_tau_conv_yrs) = 'log_tau_conv_yrs'
         profile_column_name(p_mixing_type) = 'mixing_type'
         profile_column_name(p_log_mlt_D_mix) = 'log_mlt_D_mix'
         profile_column_name(p_log_cp_T_div_t_sound) = 'log_cp_T_div_t_sound'
         profile_column_name(p_log_t_thermal) = 'log_t_thermal'
         profile_column_name(p_log_t_sound) = 'log_t_sound'
         profile_column_name(p_pressure_scale_height_cm) = 'pressure_scale_height_cm'
         profile_column_name(p_pressure_scale_height) = 'pressure_scale_height'
         profile_column_name(p_actual_gradT) = 'actual_gradT'
         profile_column_name(p_gradT_sub_actual_gradT) = 'gradT_sub_actual_gradT'
         profile_column_name(p_grada_sub_actual_gradT) = 'grada_sub_actual_gradT'

         profile_column_name(p_d_gradT_dlnd00) = 'd_gradT_dlnd00'
         profile_column_name(p_d_gradT_dlnT00) = 'd_gradT_dlnT00'
         profile_column_name(p_d_gradT_dlndm1) = 'd_gradT_dlndm1'
         profile_column_name(p_d_gradT_dlnTm1) = 'd_gradT_dlnTm1'
         profile_column_name(p_d_gradT_dlnR) = 'd_gradT_dlnR'
         profile_column_name(p_d_gradT_dL) = 'd_gradT_dL'
         profile_column_name(p_d_gradT_dln_cvpv0) = 'd_gradT_dln_cvpv0'

         profile_column_name(p_gradT) = 'gradT'
         profile_column_name(p_gradr) = 'gradr'

         profile_column_name(p_accel_div_grav) = 'accel_div_grav'

         profile_column_name(p_dlnd_dt_const_q) = 'dlnd_dt_const_q'
         profile_column_name(p_dlnT_dt_const_q) = 'dlnT_dt_const_q'

         profile_column_name(p_dlnd_dt) = 'dlnd_dt'
         profile_column_name(p_dlnT_dt) = 'dlnT_dt'
         profile_column_name(p_dlnR_dt) = 'dlnR_dt'
         profile_column_name(p_dr_dt) = 'dr_dt'
         profile_column_name(p_dv_dt) = 'dv_dt'
         profile_column_name(p_du_dt) = 'du_dt'

         profile_column_name(p_dlnd) = 'dlnd'
         profile_column_name(p_dlnT) = 'dlnT'
         profile_column_name(p_dlnR) = 'dlnR'


         profile_column_name(p_del_entropy) = 'del_entropy'
         profile_column_name(p_ds_from_eps_grav) = 'ds_from_eps_grav'

         profile_column_name(p_cno_div_z) = 'cno_div_z'

         profile_column_name(p_signed_dlnd) = 'signed_dlnd'
         profile_column_name(p_signed_dlnT) = 'signed_dlnT'

         profile_column_name(p_dE) = 'dE'
         profile_column_name(p_dr) = 'dr'
         profile_column_name(p_dv) = 'dv'
         profile_column_name(p_dr_ratio) = 'dr_ratio'
         profile_column_name(p_dt_dv_div_dr) = 'dt_dv_div_dr'

         profile_column_name(p_dlog_h1_dlogP) = 'dlog_h1_dlogP'
         profile_column_name(p_dlog_he3_dlogP) = 'dlog_he3_dlogP'
         profile_column_name(p_dlog_he4_dlogP) = 'dlog_he4_dlogP'
         profile_column_name(p_dlog_c12_dlogP) = 'dlog_c12_dlogP'
         profile_column_name(p_dlog_c13_dlogP) = 'dlog_c13_dlogP'
         profile_column_name(p_dlog_n14_dlogP) = 'dlog_n14_dlogP'
         profile_column_name(p_dlog_o16_dlogP) = 'dlog_o16_dlogP'
         profile_column_name(p_dlog_ne20_dlogP) = 'dlog_ne20_dlogP'
         profile_column_name(p_dlog_mg24_dlogP) = 'dlog_mg24_dlogP'
         profile_column_name(p_dlog_si28_dlogP) = 'dlog_si28_dlogP'

         profile_column_name(p_dlog_pp_dlogP) = 'dlog_pp_dlogP'
         profile_column_name(p_dlog_cno_dlogP) = 'dlog_cno_dlogP'
         profile_column_name(p_dlog_3alf_dlogP) = 'dlog_3alf_dlogP'

         profile_column_name(p_dlog_burn_c_dlogP) = 'dlog_burn_c_dlogP'
         profile_column_name(p_dlog_burn_n_dlogP) = 'dlog_burn_n_dlogP'
         profile_column_name(p_dlog_burn_o_dlogP) = 'dlog_burn_o_dlogP'

         profile_column_name(p_dlog_burn_ne_dlogP) = 'dlog_burn_ne_dlogP'
         profile_column_name(p_dlog_burn_na_dlogP) = 'dlog_burn_na_dlogP'
         profile_column_name(p_dlog_burn_mg_dlogP) = 'dlog_burn_mg_dlogP'

         profile_column_name(p_dlog_cc_dlogP) = 'dlog_cc_dlogP'
         profile_column_name(p_dlog_co_dlogP) = 'dlog_co_dlogP'
         profile_column_name(p_dlog_oo_dlogP) = 'dlog_oo_dlogP'

         profile_column_name(p_dlog_burn_si_dlogP) = 'dlog_burn_si_dlogP'
         profile_column_name(p_dlog_burn_s_dlogP) = 'dlog_burn_s_dlogP'
         profile_column_name(p_dlog_burn_ar_dlogP) = 'dlog_burn_ar_dlogP'
         profile_column_name(p_dlog_burn_ca_dlogP) = 'dlog_burn_ca_dlogP'
         profile_column_name(p_dlog_burn_ti_dlogP) = 'dlog_burn_ti_dlogP'
         profile_column_name(p_dlog_burn_cr_dlogP) = 'dlog_burn_cr_dlogP'
         profile_column_name(p_dlog_burn_fe_dlogP) = 'dlog_burn_fe_dlogP'

         profile_column_name(p_dlog_pnhe4_dlogP) = 'dlog_pnhe4_dlogP'
         profile_column_name(p_dlog_photo_dlogP) = 'dlog_photo_dlogP'
         profile_column_name(p_dlog_other_dlogP) = 'dlog_other_dlogP'

         profile_column_name(p_total_energy) = 'total_energy'
         profile_column_name(p_total_energy_sign) = 'total_energy_sign'
         profile_column_name(p_total_energy_integral_outward) = 'total_energy_integral_outward'
         profile_column_name(p_total_energy_integral) = 'total_energy_integral'

         profile_column_name(p_Pturb) = 'Pturb'
         profile_column_name(p_log_Pturb) = 'log_Pturb'
         profile_column_name(p_Eturb) = 'Eturb'
         profile_column_name(p_log_Eturb) = 'log_Eturb'
         profile_column_name(p_avQ) = 'avQ'
         profile_column_name(p_Hp_face) = 'Hp_face'
         profile_column_name(p_Y_face) = 'Y_face'
         profile_column_name(p_PII_face) = 'PII_face'
         profile_column_name(p_Chi) = 'Chi'
         profile_column_name(p_COUPL) = 'COUPL'
         profile_column_name(p_SOURCE) = 'SOURCE'
         profile_column_name(p_DAMP) = 'DAMP'
         profile_column_name(p_DAMPR) = 'DAMPR'
         profile_column_name(p_Eq) = 'Eq'
         profile_column_name(p_Uq) = 'Uq'
         profile_column_name(p_Lr) = 'Lr'
         profile_column_name(p_Lr_div_L) = 'Lr_div_L'
         profile_column_name(p_Lc) = 'Lc'
         profile_column_name(p_Lc_div_L) = 'Lc_div_L'
         profile_column_name(p_Lt) = 'Lt'
         profile_column_name(p_Lt_div_L) = 'Lt_div_L'

         profile_column_name(p_rsp_Et) = 'rsp_Et'
         profile_column_name(p_rsp_logEt) = 'rsp_logEt'
         profile_column_name(p_rsp_erad) = 'rsp_erad'
         profile_column_name(p_rsp_log_erad) = 'rsp_log_erad'
         profile_column_name(p_rsp_vt_div_cs) = 'rsp_vt_div_cs'
         profile_column_name(p_rsp_vt) = 'rsp_vt'
         profile_column_name(p_rsp_Pt) = 'rsp_Pt'
         profile_column_name(p_rsp_Lr) = 'rsp_Lr'
         profile_column_name(p_rsp_Lc) = 'rsp_Lc'
         profile_column_name(p_rsp_Lt) = 'rsp_Lt'
         profile_column_name(p_rsp_Eq) = 'rsp_Eq'
         profile_column_name(p_rsp_Uq) = 'rsp_Uq'
         profile_column_name(p_rsp_PII_face) = 'rsp_PII_face'
         profile_column_name(p_rsp_src) = 'rsp_src'
         profile_column_name(p_rsp_sink) = 'rsp_sink'
         profile_column_name(p_rsp_damp) = 'rsp_damp'
         profile_column_name(p_rsp_src_snk) = 'rsp_src_snk'
         profile_column_name(p_rsp_dampR) = 'rsp_dampR'
         profile_column_name(p_rsp_Y_face) = 'rsp_Y_face'
         profile_column_name(p_rsp_Hp_face) = 'rsp_Hp_face'
         profile_column_name(p_rsp_Chi) = 'rsp_Chi'
         profile_column_name(p_rsp_heat_exchange_timescale) = 'rsp_heat_exchange_timescale'
         profile_column_name(p_rsp_log_dt_div_heat_exchange_timescale) = 'rsp_log_dt_div_heat_exchange_timescale'
         profile_column_name(p_rsp_log_heat_exchange_timescale) = 'rsp_log_heat_exchange_timescale'
         profile_column_name(p_rsp_avQ) = 'rsp_avQ'
         profile_column_name(p_rsp_gradT) = 'rsp_gradT'
         profile_column_name(p_rsp_Lr_div_L) = 'rsp_Lr_div_L'
         profile_column_name(p_rsp_Lc_div_L) = 'rsp_Lc_div_L'
         profile_column_name(p_rsp_Lt_div_L) = 'rsp_Lt_div_L'

         profile_column_name(p_rsp_WORK) = 'rsp_WORK'
         profile_column_name(p_rsp_WORKQ) = 'rsp_WORKQ'
         profile_column_name(p_rsp_WORKT) = 'rsp_WORKT'
         profile_column_name(p_rsp_WORKC) = 'rsp_WORKC'

         profile_column_name(p_d_u_div_rmid) = 'd_u_div_rmid'
         profile_column_name(p_d_u_div_rmid_start) = 'd_u_div_rmid_start'

         profile_column_name(p_conv_vel_residual) = 'conv_vel_residual'
         profile_column_name(p_log_conv_vel_residual) = 'log_conv_vel_residual'
         profile_column_name(p_dconv_vel_dt) = 'dconv_vel_dt'

         profile_column_name(p_cell_specific_IE) = 'cell_specific_IE'
         profile_column_name(p_cell_ie_div_star_ie) = 'cell_ie_div_star_ie'
         profile_column_name(p_log_cell_specific_IE) = 'log_cell_specific_IE'
         profile_column_name(p_log_cell_ie_div_star_ie) = 'log_cell_ie_div_star_ie'
         
         profile_column_name(p_cell_specific_PE) = 'cell_specific_PE'
         profile_column_name(p_cell_specific_KE) = 'cell_specific_KE'
         profile_column_name(p_cell_IE_div_IE_plus_KE) = 'cell_IE_div_IE_plus_KE'
         profile_column_name(p_cell_KE_div_IE_plus_KE) = 'cell_KE_div_IE_plus_KE'

         profile_column_name(p_log_L_div_CpTMdot) = 'log_L_div_CpTMdot'
         profile_column_name(p_log_mdot_cs) = 'log_mdot_cs'
         profile_column_name(p_log_mdot_v) = 'log_mdot_v'
         profile_column_name(p_cs_at_cell_bdy) = 'cs_at_cell_bdy'
         profile_column_name(p_grad_density) = 'grad_density'
         profile_column_name(p_grad_temperature) = 'grad_temperature'

         profile_column_name(p_gradL) = 'gradL'
         profile_column_name(p_grada_sub_gradr) = 'grada_sub_gradr'
         profile_column_name(p_gradL_sub_gradr) = 'gradL_sub_gradr'
         profile_column_name(p_sch_stable) = 'sch_stable'
         profile_column_name(p_ledoux_stable) = 'ledoux_stable'

         profile_column_name(p_dominant_isoA_for_thermohaline) = 'dominant_isoA_for_thermohaline'
         profile_column_name(p_dominant_isoZ_for_thermohaline) = 'dominant_isoZ_for_thermohaline'
         profile_column_name(p_gradL_composition_term) = 'gradL_composition_term'

         profile_column_name(p_log_brunt_B) = 'log_brunt_B'
         profile_column_name(p_log_brunt_nonB) = 'log_brunt_nonB'
         profile_column_name(p_brunt_B) = 'brunt_B'
         profile_column_name(p_brunt_nonB) = 'brunt_nonB'

         profile_column_name(p_brunt_N2) = 'brunt_N2'
         profile_column_name(p_brunt_N2_structure_term) = 'brunt_N2_structure_term'
         profile_column_name(p_brunt_N2_composition_term) = 'brunt_N2_composition_term'
         profile_column_name(p_log_brunt_N2_structure_term) = 'log_brunt_N2_structure_term'
         profile_column_name(p_log_brunt_N2_composition_term) = 'log_brunt_N2_composition_term'

         profile_column_name(p_brunt_A_div_x2) = 'brunt_A_div_x2'
         profile_column_name(p_brunt_A) = 'brunt_A'
         profile_column_name(p_log_brunt_N2_dimensionless) = 'log_brunt_N2_dimensionless'
         profile_column_name(p_brunt_N2_dimensionless) = 'brunt_N2_dimensionless'
         profile_column_name(p_brunt_N_dimensionless) = 'brunt_N_dimensionless'
         profile_column_name(p_brunt_N) = 'brunt_N'
         profile_column_name(p_brunt_frequency) = 'brunt_frequency'
         profile_column_name(p_brunt_nu) = 'brunt_nu'
         profile_column_name(p_log_brunt_N) = 'log_brunt_N'
         profile_column_name(p_log_brunt_N2) = 'log_brunt_N2'
         profile_column_name(p_sign_brunt_N2) = 'sign_brunt_N2'
         profile_column_name(p_lamb_S2) = 'lamb_S2'
         profile_column_name(p_lamb_S) = 'lamb_S'
         profile_column_name(p_lamb_Sl1) = 'lamb_Sl1'
         profile_column_name(p_lamb_Sl2) = 'lamb_Sl2'
         profile_column_name(p_lamb_Sl3) = 'lamb_Sl3'
         profile_column_name(p_lamb_Sl10) = 'lamb_Sl10'
         profile_column_name(p_brunt_N_div_r_integral) = 'brunt_N_div_r_integral'
         profile_column_name(p_brunt_N2_sub_omega2) = 'brunt_N2_sub_omega2'
         profile_column_name(p_sl2_sub_omega2) = 'sl2_sub_omega2'
         profile_column_name(p_k_r_integral) = 'k_r_integral'

         profile_column_name(p_num_steps) = 'num_steps'
         profile_column_name(p_mtx_solve) = 'mtx_solve'
         profile_column_name(p_mtx_factor) = 'mtx_factor'

         profile_column_name(p_log_brunt_nu) = 'log_brunt_nu'
         profile_column_name(p_log_lamb_Sl1) = 'log_lamb_Sl1'
         profile_column_name(p_log_lamb_Sl2) = 'log_lamb_Sl2'
         profile_column_name(p_log_lamb_Sl3) = 'log_lamb_Sl3'
         profile_column_name(p_log_lamb_Sl10) = 'log_lamb_Sl10'

         profile_column_name(p_logR_kap) = 'logR_kap'
         profile_column_name(p_logW) = 'logW'
         profile_column_name(p_logQ) = 'logQ'
         profile_column_name(p_logV) = 'logV'

         profile_column_name(p_dlnX_dr) = 'dlnX_dr'
         profile_column_name(p_dlnY_dr) = 'dlnY_dr'
         profile_column_name(p_dlnRho_dr) = 'dlnRho_dr'

         profile_column_name(p_log_zFe) = 'log_zFe'
         profile_column_name(p_zFe) = 'zFe'
         profile_column_name(p_log_u_residual) = 'log_u_residual'
         profile_column_name(p_u_residual) = 'u_residual'
         profile_column_name(p_u) = 'u'
         profile_column_name(p_u_face) = 'u_face'
         profile_column_name(p_dPdr_dRhodr_info) = 'dPdr_dRhodr_info'
         profile_column_name(p_signed_log_ergs_err) = 'signed_log_ergs_err'
         profile_column_name(p_RTI_du_diffusion_kick) = 'RTI_du_diffusion_kick'
         profile_column_name(p_log_du_kick_div_du) = 'log_du_kick_div_du'
         profile_column_name(p_max_abs_xa_corr) = 'max_abs_xa_corr'

         cnt = 0
         do i=1,p_col_id_max
            if (len_trim(profile_column_name(i)) == 0) then
               write(*,*) 'missing name for profile column id', i
               if (i > 1) write(*,*) 'following ' // trim(profile_column_name(max(1,i-1))) ! bp: get rid of bogus compiler warning
               write(*,*)
               cnt = cnt+1
            end if
         end do

         if (cnt > 0) then
            ierr = -1
            return
         end if

         nullify(profile_column_names_dict)
         do i=1,p_col_id_max
            call integer_dict_define(profile_column_names_dict, profile_column_name(i), i, ierr)
            if (ierr /= 0) then
               write(*,*) 'FATAL ERROR: profile_column_names_init failed in integer_dict_define'
               return
            end if
         end do

      end subroutine profile_column_names_init


      subroutine profile_column_names_shutdown()
        use utils_lib, only: integer_dict_free
        if (ASSOCIATED(profile_column_names_dict)) call integer_dict_free(profile_column_names_dict)
      end subroutine profile_column_names_shutdown


      integer function do_get_profile_id(cname)
         use utils_lib
         character (len=*), intent(in)  :: cname
         ! returns id for the profile column if there is a matching name
         ! returns 0 otherwise.
         integer :: ierr, value
         call integer_dict_lookup(profile_column_names_dict, cname, value, ierr)
         if (ierr /= 0) value = 0
         do_get_profile_id = value
      end function do_get_profile_id



      end module star_profile_def

