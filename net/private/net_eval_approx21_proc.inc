            ierr = 0
            
            call approx21_special_reactions( &
               btemp, bden, abar, zbar, n% y, &
               g% which_rates(ir_he4_he4_he4_to_c12) == use_rate_3a_FL87, &
               Qconv*reaction_Qs(ir_he4_he4_he4_to_c12), &
               rate_screened, rate_screened_dT, rate_screened_dRho, &
               dratdumdy1, dratdumdy2, ierr)
            if (ierr /= 0) return            
            
            call approx21_dydt( &
               n% y, rate_screened, rate_screened, &
               dydt1, .false., g% fe56ec_fake_factor, g% min_T_for_fe56ec_fake_factor, &
               g% fe56ec_n_neut, btemp, ierr)
            if (ierr /= 0) return
               
            fII = approx21_eval_PPII_fraction(n% y, rate_screened)
            
            if (g% add_co56_to_approx21) then
               call get_approx21_plus_co56_eps_info( &
                     dydt1, rate_screened, .true., eps_total, eps_neu_total, ierr)
            else
               call get_approx21_eps_info( &
                     dydt1, rate_screened, .true., eps_total, eps_neu_total, ierr)
            end if
                  
            if (ierr /= 0) return               
            eps_nuc = eps_total - eps_neu_total
            
            do i=1,num_isos
               dxdt(i) = chem_isos% Z_plus_N(chem_id(i))*dydt1(i)
            end do

            if (just_dxdt) return

            call approx21_dfdy( &
               n% y, dfdy, &
               g% fe56ec_fake_factor, g% min_T_for_fe56ec_fake_factor, g% fe56ec_n_neut, &
               rate_screened, rate_screened_dT, rate_screened_dRho, &
               dratdumdy1, dratdumdy2, btemp, ierr)
            if (ierr /= 0) return

            call approx21_dfdT_dfdRho( & 
               
               ! NOTE: currently this gives d_eps_total_dy -- should fix to account for neutrinos too
               
               n% y, mion, dfdy, rate_screened, rate_screened_dT, rate_screened_dRho, &
               g% fe56ec_fake_factor, g% min_T_for_fe56ec_fake_factor, &
               g% fe56ec_n_neut, btemp, dfdT, dfdRho, d_epsnuc_dy, ierr)
            if (ierr /= 0) return

            if (g% add_co56_to_approx21) then
               call get_approx21_plus_co56_eps_info( &
                  dfdT, rate_screened_dT, .false., deps_total_dT, deps_neu_dT, ierr)
            else
               call get_approx21_eps_info( &
                  dfdT, rate_screened_dT, .false., deps_total_dT, deps_neu_dT, ierr)
            end if
            if (ierr /= 0) return
            d_eps_nuc_dT = deps_total_dT - deps_neu_dT
                              
            if (g% add_co56_to_approx21) then
               call get_approx21_plus_co56_eps_info( &
                  dfdRho, rate_screened_dRho, .false., deps_total_dRho, deps_neu_dRho, ierr)
            else
               call get_approx21_eps_info( &
                  dfdRho, rate_screened_dRho, .false., deps_total_dRho, deps_neu_dRho, ierr)
            end if
            if (ierr /= 0) return             
            d_eps_nuc_dRho = deps_total_dRho - deps_neu_dRho
            
            call approx21_d_epsneu_dy( &
               n% y, rate_screened, &
               reaction_neuQs(irpp_to_he3), &   
               reaction_neuQs(ir34_pp2), &  
               reaction_neuQs(ir34_pp3), &  
               reaction_neuQs(irc12_to_n14), &  
               reaction_neuQs(irn14_to_c12), &  
               reaction_neuQs(iro16_to_n14), &  
               d_epsneu_dy, ierr)
            if (ierr /= 0) return

            do i=1,num_isos
               ci = chem_id(i)
               Z_plus_N = dble(chem_isos% Z_plus_N(ci))
               d_eps_nuc_dx(i) = (d_epsnuc_dy(i) - d_epsneu_dy(i))/Z_plus_N 
               d_dxdt_dRho(i) = Z_plus_N*dfdRho(i)
               d_dxdt_dT(i) = Z_plus_N*dfdT(i)
               do j=1, num_isos
                  d_dxdt_dx(i,j) = &
                     dfdy(i,j)*Z_plus_N/chem_isos% Z_plus_N(chem_id(j))
               end do
            end do
