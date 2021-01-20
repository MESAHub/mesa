      module test_eos_support

      use eos_def
      use eos_lib
      use chem_lib
      use chem_def
      use const_def
      use eos_support
      use math_lib
      
      implicit none

      real(dp), dimension(:), pointer :: rho_vector, T_vector
      real(dp), dimension(:,:), pointer :: results, results_saved
      
      contains
      
      
      subroutine get_logP20_contour
         real(dp) :: Z, X, logPgas, logT, logRho, logP, T, Prad, Pgas, P
         real(dp) :: logT_min, logT_max, dlogT
         integer :: iounit
         include 'formats'
         call Setup_eos
         Z =  0.02d0
         X =  0.0d0
         logT_min = 6d0
         logT_max = 8.59d0
         dlogT = 0.01d0
         
         write(*,'(99(a20,4x))') 'logRho', 'logT', 'logP'
         
         logT = logT_min
         do while (logT < logT_max)
            T = exp10(logT)
            logP = 20d0
            P = exp10(logP)
            Prad = crad*T*T*T*T/3
            Pgas = P - Prad
            logPgas = log10(Pgas)
            call test1_eosPT(Z, X, logPgas, logT, .true., .true., logRho, logP)
            write(*,'(99(f20.10,4x))') logRho, logT, logP
            logT = logT + dlogT
         end do
         
         write(*,*)
         
      end subroutine get_logP20_contour
      
      
      subroutine get_logP23_contour
         real(dp) :: Z, X, logPgas, logT, logRho, logP, T, Prad, Pgas, P
         real(dp) :: logT_min, logT_max, dlogT
         integer :: iounit
         include 'formats'
         call Setup_eos
         Z =  0.02d0
         X =  0.0d0
         logT_min = 7d0
         logT_max = 8.9d0
         dlogT = 0.01d0
         
         write(*,'(99(a20,4x))') 'logRho', 'logT', 'logP'
         
         logT = logT_min
         do while (logT < logT_max)
            T = exp10(logT)
            logP = 23d0
            P = exp10(logP)
            Prad = crad*T*T*T*T/3
            Pgas = P - Prad
            logPgas = log10(Pgas)
            call test1_eosPT(Z, X, logPgas, logT, .true., .true., logRho, logP)
            write(*,'(99(f20.10,4x))') logRho, logT, logP
            logT = logT + dlogT
         end do
         
      end subroutine get_logP23_contour
      
      
      subroutine test1_eosPT_for_ck(quietly)
         logical, intent(in) :: quietly
         real(dp) :: Z, X, logPgas, logT, logRho, logP
         ! pick args for test so that use both HELM and tables to get results.
         Z =  0.02d0
         X =  0.7d0
         logT = 3.0d0
         logPgas = 4.2d0
         call test1_eosPT(Z, X, logPgas, logT, .false., quietly, logRho, logP)
      end subroutine test1_eosPT_for_ck

      
      subroutine test_eosPT(which)
         integer, intent(in) :: which
         real(dp) :: Z, X, logPgas, logT, logRho, logP
         logical :: save_use_max_SCVH_for_PT
         Z =  0.02d0
         X =  0.6d0
         logT = 4d0
         logPgas = 5d0
         save_use_max_SCVH_for_PT = rq% use_max_SCVH_for_PT
         if (which == 1) then
            rq% use_max_SCVH_for_PT = .true.
         else
            rq% use_max_SCVH_for_PT = .false.
         end if
         write(*,*) 'test_eosPT rq% use_max_SCVH_for_PT', rq% use_max_SCVH_for_PT
         call test1_eosPT(Z, X, logPgas, logT, .false., .false., logRho, logP)
         rq% use_max_SCVH_for_PT = save_use_max_SCVH_for_PT
      end subroutine test_eosPT
            
      
      subroutine test1_eosPT(Z, X, logPgas, logT, do_compare, quietly, logRho, logP)
         logical, intent(in) :: quietly
         real(dp) :: Z, X, logPgas, logT
         real(dp), intent(out) :: logRho, logP
         logical, intent(in) :: do_compare
         real(dp) :: &
               P, Pgas, Prad, T, Rho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas, &
               res(num_eos_basic_results), d_dlnd(num_eos_basic_results), &
               d_dlnT(num_eos_basic_results), &
               d_dabar(num_eos_basic_results), d_dzbar(num_eos_basic_results), &
               res2(num_eos_basic_results), d_dlnd2(num_eos_basic_results), &
               d_dabar2(num_eos_basic_results), d_dzbar2(num_eos_basic_results), &
               d_dlnT2(num_eos_basic_results)
         integer:: ierr, i
         logical, parameter :: use_log10_for_other = .false.
         character (len=eos_name_length) :: names(num_eos_basic_results)
         
         
         include 'formats'
 
         ierr = 0

         call Init_Composition(X, Z, 0d0, 0d0) ! sets abar and zbar
         
         if (.false.) then ! TESTING
            z =     0d0
            x =     0.72d0
            call Init_Composition(X, Z, 0d0, 0d0) ! sets abar and zbar
            abar =     1.2966413082679851D+00
            zbar =     1.1021433867453336D+00
            logPgas  = 4.8066181993619859D+00
            logT  = 3.7569035961895620D+00
         end if
         
         T = exp10(logT)
         Pgas = exp10(logPgas)

         if (.not. quietly) then
            write(*,*) 'test1_eosPT'
            write(*,1) 'logT', logT
            write(*,1) 'logPgas', logPgas
            write(*,1) 'T', T
            write(*,1) 'Pgas', Pgas
            write(*,1) 'Z', Z
            write(*,1) 'X', X
            write(*,1) 'abar', abar
            write(*,1) 'zbar', zbar
            write(*,*)
         end if

         call eosPT_get( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               Pgas, logPgas, T, logT, &
               Rho, logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas, &
               res, d_dlnd, d_dlnT, d_dabar, d_dzbar, ierr)
         if (ierr /= 0) then
            write(*,*) 'ierr in eosPT_get for test1_eosPT'
            stop 1
         end if
         
         Prad = crad*T*T*T*T/3
         P = Pgas + Prad
         logP = log10(P)
         
         if (quietly) return
         
         write(*,*)
         write(*,1) 'rho', rho
         write(*,1) 'logRho', logRho
         write(*,1) 'logT', logT
         write(*,1) 'logPgas', logPgas
         write(*,*)
     
         call eos_result_names(names)
         
         if (.not. do_compare) then ! simple form of output
            write(*,1) 'dlnRho_dlnPgas_c_T', dlnRho_dlnPgas_c_T
            write(*,1) 'dlnRho_dlnT_c_Pgas', dlnRho_dlnT_c_Pgas
            write(*,*)
            write(*,1) 'P', P
            write(*,1) 'E', exp(res(i_lnE))
            write(*,1) 'S', exp(res(i_lnS))
            write(*,*)
            do i=4, num_eos_basic_results
               write(*,1) trim(names(i)), res(i)
            end do
            write(*,*)
            return
         end if

         call eosDT_get_legacy( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               rho, logRho, T, logT, &
               res2, d_dlnd2, d_dlnT2, d_dabar2, d_dzbar2, ierr)
         if (ierr /= 0) then
            write(*,*) 'ierr in eosDT_get for test1_eosPT'
            stop 1
         end if
     
         write(*,*)
         
         write(*,1) 'dlnRho_dlnPgas_c_T', dlnRho_dlnPgas_c_T
         write(*,1) 'dlnRho_dlnT_c_Pgas', dlnRho_dlnT_c_Pgas
         do i=1, num_eos_basic_results
            write(*,1) trim(names(i)), res(i), res2(i), &
                  (res(i)-res2(i)) / max(1d0, abs(res(i)), abs(res2(i)))
         end do
         write(*,*)
         
         do i=1, num_eos_basic_results
            write(*,1) 'd_dlnd ' // trim(names(i)), d_dlnd(i), d_dlnd2(i), &
                  (d_dlnd(i)-d_dlnd2(i)) / max(1d0, abs(d_dlnd(i)), abs(d_dlnd2(i)))
         end do
         write(*,*)
         
         do i=1, num_eos_basic_results
            write(*,1) 'd_dlnT ' // trim(names(i)), d_dlnT(i), &
                  d_dlnT2(i), (d_dlnT(i)-d_dlnT2(i)) / max(1d0, abs(d_dlnT(i)), abs(d_dlnT2(i)))
         end do
         write(*,*)
         
      end subroutine test1_eosPT
      
      
      subroutine test1_eosDE_for_ck(quietly)
         logical, intent(in) :: quietly
         real(dp) :: Z, X, logRho, logE, logT_guess, logT, logPgas
         logical :: do_compare
         
         !write(*,*) 'test in table'
         Z =  0.02d0
         X =  0.70d0
         logRho = -4.011d0
         logE = 14.011d0
         logT_guess = 4d0
         do_compare = .false.
         call test1_eosDE( &
            Z, X, logRho, logE, logT_guess, do_compare, quietly, logT, logPgas)
     
         !write(*,*) 'test off table'
         Z =  0.02d0
         X =  0.70d0
         logRho = 6d0
         logE = 17.2d0
         logT_guess = 8d0
         call test1_eosDE( &
            Z, X, logRho, logE, logT_guess, do_compare, quietly, logT, logPgas)
     
      end subroutine test1_eosDE_for_ck

      
      subroutine test_eosDE
         real(dp) :: Z, X, logRho, logE, logT_guess, logT, logPgas
         logical :: do_compare, quietly

         do_compare = .false.
         quietly = .false.

         Z =  9.9961454947386075D-01
         X =  4.5202251013642486D-77
         
         ! off table
         logRho = 6d0
         logE = 17.2d0
         logT_guess = 8.0d0
         
         ! in table
         logRho = -4.011d0
         logE = 14.011d0
         logT_guess = 4d0
         
         ! test1
         logRho = -7.7958201868593138D+00
         logE = 1.2067605390816265D+01
         logT_guess = 4.4670369613818757D+00
         
         call test1_eosDE( &
            Z, X, logRho, logE, logT_guess, do_compare, quietly, logT, logPgas)
     
         return
         
         
         
         ! test2
         logRho = -1.0793437870766159D+01
         logE = 1.5705471856124511D+01
         logT_guess = 4.4670369613818757D+00
         call test1_eosDE( &
            Z, X, logRho, logE, logT_guess, do_compare, quietly, logT, logPgas)
         
         write(*,*)
         write(*,*) 'dlogT', 4.4670383913163896D+00 - 4.4670369624727995D+00
         write(*,*) 'dlogRho', -7.7958123910391270D+00 - (-7.7958201868593138D+00)
         write(*,*) 'dlogT/dlogRho', &
            (4.4670383913163896D+00 - 4.4670369624727995D+00)/&
                 (-7.7958123910391270D+00 - (-7.7958201868593138D+00))
         write(*,*)
     
      end subroutine test_eosDE
      
      
      subroutine test1_eosDE( &
            Z, X, logRho, logE, logT_guess, do_compare, quietly, logT, logPgas)
         logical, intent(in) :: quietly
         real(dp) :: Z, X, logRho, logE, logT_guess
         real(dp), intent(out) :: logT, logPgas
         logical, intent(in) :: do_compare
         real(dp) :: &
               energy, rho, T, &
               dlnT_dlnE_c_Rho, dlnT_dlnd_c_E, &
               dlnPgas_dlnE_c_Rho, dlnPgas_dlnd_c_E, &
               res(num_eos_basic_results), d_dlnd(num_eos_basic_results), &
               d_dlnT(num_eos_basic_results), &
               d_dabar(num_eos_basic_results), d_dzbar(num_eos_basic_results), &
               res2(num_eos_basic_results), d_dlnd2(num_eos_basic_results), &
               d_dabar2(num_eos_basic_results), d_dzbar2(num_eos_basic_results), &
               d_dlnT2(num_eos_basic_results)    
         integer:: ierr, i
         logical, parameter :: use_log10_for_other = .false.
         character (len=eos_name_length) :: names(num_eos_basic_results)
         integer :: handle1, species1
         real(dp) :: abar1, zbar1, z53bar1
         integer, pointer, dimension(:) :: chem_id1, net_iso1
         real(dp) :: xa1(species)
                  
         include 'formats'
 
         ierr = 0

         call Init_Composition(X, Z, 0d0, 0d0) ! sets abar and zbar
         
         handle1 = handle
         species1 = species
         abar1 = abar
         zbar1 = zbar
         z53bar1 = z53bar
         chem_id1 => chem_id
         net_iso1 => net_iso
         xa1 = xa
         
         if (.false.) then ! TESTING
            z = 8.4114889753945277D-01
            x = 5.8622173195316740D-76
            call Init_Composition(X, Z, 0d0, 0d0) ! sets abar and zbar
            abar =  1.8102836455287726D+01
            zbar =  8.9391852734580954D+00
            logRho  = -1.0005343880781572D+01
            logE  = 1.5782735151459207D+01
         end if
         
         energy = exp10(logE)
         rho = exp10(logRho)

         if (.not. quietly) then
            write(*,*) 'test1_eosDE'
            write(*,1) 'logRho', logRho
            write(*,1) 'logE', logE
            write(*,1) 'logT_guess', logT_guess
            write(*,1) 'Z', Z
            write(*,1) 'X', X
            write(*,1) 'abar', abar
            write(*,1) 'zbar', zbar
            write(*,*)
         end if
         
         call eosDE_get( &
               handle1, Z, X, abar1, zbar1, &
               species1, chem_id1, net_iso1, xa1, &
               energy, logE, rho, logRho, logT_guess, &
               T, logT, res, d_dlnd, d_dlnT, d_dabar, d_dzbar, &
               dlnT_dlnE_c_Rho, dlnT_dlnd_c_E, &
               dlnPgas_dlnE_c_Rho, dlnPgas_dlnd_c_E, &
               ierr)
         if (ierr /= 0) then
            write(*,*) 'ierr in eosDE_get for test1_eosDE'
            stop 1
         end if
         
         if (quietly) return
         
         write(*,*)
         write(*,1) 'T', T
         write(*,1) 'logT', logT
         write(*,*)
     
         call eos_result_names(names)
         logPgas = res(i_lnPgas)/ln10
         if (.not. do_compare) then ! simple form of output
            write(*,1) 'logPgas', res(i_lnPgas)/ln10
            write(*,1) 'logS', res(i_lnS)/ln10
            do i=4, num_eos_basic_results
               write(*,1) trim(names(i)), res(i)
            end do
            write(*,1) 'dlnT_dlnE_c_Rho', dlnT_dlnE_c_Rho
            write(*,1) 'dlnT_dlnd_c_E', dlnT_dlnd_c_E
            write(*,1) 'dlnPgas_dlnE_c_Rho', dlnPgas_dlnE_c_Rho
            write(*,1) 'dlnPgas_dlnd_c_E', dlnPgas_dlnd_c_E
            write(*,*)
            return
         end if

         call eosDT_get_legacy( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               rho, logRho, T, logT, &
               res2, d_dlnd2, d_dlnT2, d_dabar2, d_dzbar2, ierr)
         if (ierr /= 0) then
            write(*,*) 'ierr in eosDT_get for test1_eosDE'
            stop 1
         end if
     
         write(*,*)
         
         do i=1, num_eos_basic_results
            write(*,1) trim(names(i)), res(i), res2(i), &
                  (res(i)-res2(i)) / max(1d0, abs(res(i)), abs(res2(i)))
         end do
         write(*,*)
         
         do i=1, num_eos_basic_results
            write(*,1) 'd_dlnd ' // trim(names(i)), d_dlnd(i), d_dlnd2(i), &
                  (d_dlnd(i)-d_dlnd2(i)) / max(1d0, abs(d_dlnd(i)), abs(d_dlnd2(i)))
         end do
         write(*,*)
         
         do i=1, num_eos_basic_results
            write(*,1) 'd_dlnT ' // trim(names(i)), d_dlnT(i), &
                  d_dlnT2(i), (d_dlnT(i)-d_dlnT2(i)) / max(1d0, abs(d_dlnT(i)), abs(d_dlnT2(i)))
         end do
         write(*,*)
           
      end subroutine test1_eosDE
      
      
      subroutine test_eosDT(which)
         integer, intent(in) :: which
         logical :: save_use_FreeEOS

         include 'formats'
         
         save_use_FreeEOS = rq% use_FreeEOS
         rq% use_FreeEOS = .false.
         if (which == 1 .or. which == 2) then
            rq% use_FreeEOS = .true.
         end if
         
         write(*,*) 'test_eosDT rq% use_FreeEOS', rq% use_FreeEOS
         if (which == 2) then
            call do_test_eosDT_new
            stop
         end if
         call do_test_eosDT
         rq% use_FreeEOS = save_use_FreeEOS
         
      end subroutine test_eosDT
      
      
      subroutine do_test_eosDT
         real(dp) :: logT, logRho, T, Rho, X, Z, logPgas, Pgas, Prad, P, logP, XC, XO, &
               logRho1_OPAL_SCVH_limit, logRho2_OPAL_SCVH_limit, Z_all_HELM, Z_all_OPAL
         real(dp), dimension(num_eos_basic_results) :: &
               res, d_dlnd, d_dlnT, d_dabar, d_dzbar
         integer :: info, i_var, i
         real(dp) :: dx_0, var_0, dvardx_0, dvardx, xdum, err, lnd, lnT, energy, entropy, &
            var, dvar_dlnd, dvar_dlnT, dse, dpe, dsp, dS_dT, dE_dT, dP_dT, dPrad_dT, &
            dE_dRho, dS_dRho, dlnE_dlnd, dlnE_dlnT, dlnPgas_dlnT, dlnS_dlnd, dlnS_dlnT
         logical :: doing_d_dlnd

         include 'formats'
         
         
         if (.false.) then
                                                rho = 952987205766.07312d0
                                                T =  3162.2778229206556d0
                                              logRho = log10(rho)
                                                logT = log10(T)
         else
                                              logRho = -3.7440634305717460D+00
                                                logT = 3.9830887848306653D+00
                                                rho = exp10(logRho)
                                                T = exp10(logT)
         end if
                                                   
                 
         XC = 0
         XO = 0
                                                   z =  2.0745049657628911D-02
                                                   x =  1.5302086314738124D-30
         call Init_Composition(X, Z, XC, XO)
         
         write(*,*) 'change abar and zbar for test'
         abar = 1.2198718692960793d0
         zbar = 1.0733094612841689d0
         
         write(*,1) 'abar = ', abar
         write(*,1) 'zbar = ', zbar
         
         if (.true.) then
            rq% Z_all_HELM = 0.04d0
            rq% Z_all_OPAL = 0.039d0
            write(*,1) 'Z_all_HELM Z_all_OPAL', Z_all_HELM, Z_all_OPAL
         end if

         write(*,*) 'use_PC', rq% use_PC
         write(*,1) 'log_Gamma_e_all_HELM', rq% log_Gamma_e_all_HELM
         write(*,1) 'log_Gamma_e_all_PC', rq% log_Gamma_e_all_PC
         write(*,1) 'logRho2_PC_limit', rq% logRho2_PC_limit
         write(*,1) 'logRho1_PC_limit', rq% logRho1_PC_limit

         write(*,*) 'call eosDT_get'
         info = 0
         call eosDT_get_legacy( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               Rho, logRho, T, logT, &
               res, d_dlnd, d_dlnT, d_dabar, d_dzbar, info)
         if (info /= 0) then
            write(*,*) 'failed in test1_eos'
            stop 1
         end if
         
         write(*,1) 'T', T
         write(*,1) 'logT', logT
         write(*,1) 'Rho', Rho
         write(*,1) 'logRho', logRho
         write(*,1) 'abar', abar
         write(*,1) 'zbar', zbar
         write(*,1) 'X', X
         write(*,1) 'Z', Z
         write(*,1) 'rq% logRho_min_OPAL_SCVH_limit', rq% logRho_min_OPAL_SCVH_limit
         write(*,1) 'logQ', logRho - 2d0*logT + 12d0
         write(*,*)
         
         if (.false.) then ! dfridr
         
            !stop

            i_var = i_gamma1
         
            !doing_d_dlnd = .true.
            doing_d_dlnd = .false.
         
            lnd = logRho*ln10
            lnT = logT*ln10

            if (doing_d_dlnd) then
               dx_0 = max(1d-6, abs(lnd*1d-6))
               dvardx_0 = d_dlnd(i_var)
            else
               dx_0 = max(1d-6, abs(lnT*1d-6))
               dvardx_0 = d_dlnT(i_var)
            end if
            err = 0d0
            dvardx = dfridr(dx_0,err)
            xdum = (dvardx - dvardx_0)/max(abs(dvardx_0),1d-4)
            if (doing_d_dlnd) then
               write(*,1) 'd_dlnRho analytic, numeric, diff, rel diff', &
                     dvardx_0, dvardx, err, xdum
            else ! doing d_dlnT
               write(*,1) 'd_dlnT analytic, numeric, diff, rel diff', &
                     dvardx_0, dvardx, err, xdum
            end if
            write(*,*)
            
            return
         
         end if
         
         energy = exp(res(i_lnE))
         entropy = exp(res(i_lnS))
         Pgas = exp(res(i_lnPgas))
         Prad = Radiation_Pressure(T)
         P = Pgas + Prad
         dPrad_dT = 4d0*Prad/T
         
         
         if (.true.) then ! simple form of output
            write(*,*)
            write(*,1) 'rho', rho
            write(*,1) 'logRho', logRho
            write(*,1) 'logT', logT
            write(*,1) 'logPgas', res(i_lnPgas)/ln10
            write(*,*)
            write(*,1) 'P', P
            write(*,1) 'E', exp(res(i_lnE))
            write(*,1) 'S', exp(res(i_lnS))
            write(*,*)
            do i=4, num_eos_basic_results
               write(*,1) trim(eos_names(i)), res(i)
            end do
            write(*,*)
            return
         end if

         
         
         
         dlnS_dlnT = d_dlnT(i_lnS)
         dlnS_dlnd = d_dlnd(i_lnS)
         dlnE_dlnT = d_dlnT(i_lnE)
         dlnE_dlnd = d_dlnd(i_lnE)
         dlnPgas_dlnT = d_dlnT(i_lnPgas)
         
         dS_dT = dlnS_dlnT*entropy/T
         dS_dRho = dlnS_dlnd*entropy/rho
         dE_dT = dlnE_dlnT*energy/T
         dE_dRho = dlnE_dlnd*energy/rho
         dP_dT = dlnPgas_dlnT*Pgas/T + dPrad_dT
         
         write(*,1) 'chiRho', res(i_chiRho)
         write(*,1) 'd_dlnd(i_chiRho)', d_dlnd(i_chiRho)
         write(*,1) 'd_dlnT(i_chiRho)', d_dlnT(i_chiRho)
         !stop
         
         write(*,1) 'grad_ad', res(i_grad_ad)
         write(*,1) 'd_dlnd(i_grad_ad)', d_dlnd(i_grad_ad)
         write(*,1) 'd_dlnT(i_grad_ad)', d_dlnT(i_grad_ad)
         write(*,1) 'chiT', res(i_chiT)
         write(*,1) 'd_dlnd(i_chiT)', d_dlnd(i_chiT)
         write(*,1) 'd_dlnT(i_chiT)', d_dlnT(i_chiT)
         write(*,1) 'Pgas', Pgas
         write(*,1) 'logPgas', res(i_lnPgas)/ln10
         !return
         
         write(*,1) 'P', P
         write(*,1) 'logP', log10(P)
         write(*,1) 'Pgas', Pgas
         write(*,1) 'logPgas', res(i_lnPgas)/ln10
         write(*,1) 'chiRho', res(i_chiRho)
         write(*,1) 'Cp', res(i_Cp)
         write(*,1) 'Cv', res(i_Cv)
         write(*,1) 'dE_dRho', res(i_dE_dRho)
         write(*,1) 'dS_dT', res(i_dS_dT)
         write(*,1) 'dS_dRho', res(i_dS_dRho)
         write(*,1) 'mu', res(i_mu)
         write(*,1) 'gamma1', res(i_gamma1)
         write(*,1) 'gamma3', res(i_gamma3)
         write(*,1) 'eta', res(i_eta)
         write(*,1) 'lnfree_e', res(i_lnfree_e)
         write(*,*)
         
         ! dse = T ∂S/∂T|_⍴ / ∂E/∂T|_⍴  - 1.0d0
         ! dpe = (⍴^2 ∂E/∂⍴|_T + T ∂P/∂T_⍴) / P - 1.0d0
         ! dsp = -(∂S/∂⍴|_T * ⍴^2) / ∂P/∂T|_⍴ - 1.0d0
         dse = T*dS_dT/dE_dT - 1d0
         dpe = (rho*rho*dE_dRho + T*dP_dT)/P - 1d0
         dsp = -rho*rho*dS_dRho/dP_dT - 1d0

         write(*,1) 'dS_dT', dS_dT
         write(*,1) 'dE_dT', dE_dT
         write(*,1) 'dP_dT', dP_dT
         write(*,1) 'dE_dRho', dE_dRho
         write(*,1) 'dS_dRho', dS_dRho
         write(*,*)
         
         write(*,1) 'T', T, 1.00000000000000D+08
         write(*,1) 'rho', rho, 1.00000000000000D+02
         write(*,1) 'Pgas', Pgas, 1.23652998034023D+18
         write(*,1) 'logT', logT, 8.00000000000000D+00
         write(*,1) 'logRho', logRho, 2.00000000000000D+00
         write(*,1) 'logPgas', res(i_lnPgas)/ln10, 1.80922046505303D+01
         write(*,1) 'logE', res(i_lnE)/ln10, 1.64204279937121D+01
         write(*,1) 'logS', res(i_lnS)/ln10, 9.21518793962630D+00
         write(*,1) 'dlnPgas_dlnT', d_dlnd(i_chiT), 9.99407523743355E-01
         write(*,1) 'dlnPgas_dlnd', d_dlnT(i_chiT), 1.00063526337044D+00
         write(*,1) 'd2lnPgas_dlnd_dlnT', 0d0, -1.00356261967249E-03
         write(*,1) 'dlnE_dlnT', d_dlnd(i_lnE), 1.86915994288048D+00
         write(*,1) 'dlnE_dlnd', d_dlnT(i_lnE), -2.87077214678670E-01
         write(*,1) 'd2lnE_dlnd_dlnT', 0d0, -1.43403816995508E-01
         write(*,1) 'dlnS_dlnT', d_dlnd(i_lnS), 2.99837896996580E-01
         write(*,1) 'dlnS_dlnd', d_dlnT(i_lnS), -1.36754620123785E-01
         write(*,1) 'd2lnS_dlnd_dlnT', 0d0, -1.43403816995508E-01
         write(*,1) 'dse', dse, 8.73492878028514E-09
         write(*,1) 'dpe', dpe, -6.04495641759601E-01
         write(*,1) 'dsp', dsp, 8.01573628057714E-01
         write(*,*)
         
         stop
         
         contains
         
         
         real(dp) function dfridr_func(delta_x) result(val)
            real(dp), intent(in) :: delta_x
            integer :: ierr
            real(dp) :: log_var
            include 'formats'
            ierr = 0
            if (doing_d_dlnd) then
               log_var = (lnd + delta_x)/ln10
               var = exp10(log_var)
               call eosDT_get_legacy( &
                  handle, Z, X, abar, zbar, &
                  species, chem_id, net_iso, xa, &
                  var, log_var, T, logT, &
                  res, d_dlnd, d_dlnT, d_dabar, d_dzbar, info)
            else
               log_var = (lnT + delta_x)/ln10
               var = exp10(log_var)
               call eosDT_get_legacy( &
                  handle, Z, X, abar, zbar, &
                  species, chem_id, net_iso, xa, &
                  Rho, logRho, var, log_var, &
                  res, d_dlnd, d_dlnT, d_dabar, d_dzbar, info)
            end if
            val = res(i_var)
         end function dfridr_func

         real(dp) function dfridr(hx,err) ! from Frank
            real(dp), intent(in) :: hx
            real(dp), intent(out) :: err
            !  this routine returns the first derivative of a function func(x)
            !  at the point x, by ridders method of polynomial extrapolation.
            !  value hx is the initial step size;
            !  it should be an increment for which func changes substantially.
            !  an estimate of the error in the first derivative is returned in err.
            integer, parameter :: ntab = 20
            integer :: i,j
            real(dp) :: x,errt,fac,hh,a(ntab,ntab),xdum,ydum
            real(dp), parameter :: con2=2d0, con=sqrt(con2), big=1d50, safe=2d0
            include 'formats'
            dfridr = 0d0
            hh = hx
            ! 2nd order central difference
            a(1,1) = (dfridr_func(hh) - dfridr_func(-hh))/(2d0*hh)
            write(*,2) 'dfdx hh', 1, a(1,1), hh
            err = big
            ! succesive columns in the neville tableu will go to smaller stepsizes
            ! and higher orders of extrapolation
            do i=2,ntab
               hh = hh/con
               a(1,i) = (dfridr_func(hh) - dfridr_func(-hh))/(2d0*hh)
               !write(*,2) 'dfdx hh', i, a(1,i), hh
               ! compute extrapolations of various orders; the error stratagy is to compare
               ! each new extrapolation to one order lower but both at the same stepsize
               ! and at the previous stepsize
               fac = con2
               do j=2,i
                  a(j,i) = (a(j-1,i)*fac - a(j-1,i-1))/(fac-1d0)
                  fac = con2*fac
                  errt = max(abs(a(j,i)-a(j-1,i)),abs(a(j,i)-a(j-1,i-1)))
                  if (errt <= err) then
                     err = errt
                     dfridr = a(j,i)
                     write(*,3) 'dfridr err', i, j, dfridr, err
                  end if
               end do
               ! if higher order is worse by a significant factor safe, then bail
               if (abs(a(i,i) - a(i-1,i-1)) >= safe*err) then
                  write(*,1) 'higher order is worse', err, a(i,i), a(i-1,i-1)
                  return
               end if
            end do
         end function dfridr

      end subroutine do_test_eosDT
      
         
      
      subroutine do_test_eosDT_new
         use num_lib, only : dfridr
         real(dp) :: logT, logRho, T, Rho, X, Z, logPgas, Pgas, Prad, P, logP, XC, XO
         real(dp), dimension(num_eos_basic_results) :: &
            res, d_dlnd, d_dlnT, d_dabar, d_dzbar
         real(dp), dimension(num_eos_d_dxa_results,species) :: d_dxa
         integer :: info, i_var, i
         real(dp) :: dx_0, var_0, dvardx_0, dvardx, xdum, err, lnd, lnT
         logical :: doing_d_dlnd

         real(dp) :: xh, xhe, zz, abar_ci, zbar_ci, z2bar_ci, z53bar_ci, ye_ci, mass_correction, sumx
         real(dp), dimension(species) :: dabar_dx, dzbar_dx, dmc_dx


         include 'formats'


         logRho = -1.23d0
         logT = 8.1d0
         rho = exp10(logRho)
         T = exp10(logT)

         XC = 0
         XO = 0
         z =  0.014
         x =  0.314
         call Init_Composition(X, Z, XC, XO)

         ! need this call to get dabar_dx, dzbar_dx
         ! might want to pass these in eventually
         call composition_info( &
            species, chem_id, xa, xh, xhe, zz, &
            abar_ci, zbar_ci, z2bar_ci, z53bar_ci, ye_ci, mass_correction, &
            sumx, dabar_dx, dzbar_dx, dmc_dx)

         write(*,*) 'call eosDT_get'
         info = 0
         call eosDT_get( &
            handle, species, chem_id, net_iso, xa, &
            Rho, logRho, T, logT, &
            res, d_dlnd, d_dlnT, d_dxa, info)
         if (info /= 0) then
            write(*,*) 'failed in test1_eos'
            stop 1
         end if

         write(*,1) 'T', T
         write(*,1) 'logT', logT
         write(*,1) 'Rho', Rho
         write(*,1) 'logRho', logRho
         write(*,1) 'abar', abar
         write(*,1) 'zbar', zbar
         write(*,1) 'X', X
         write(*,1) 'Z', Z
         write(*,1) 'rq% logRho_min_OPAL_SCVH_limit', rq% logRho_min_OPAL_SCVH_limit
         write(*,1) 'logQ', logRho - 2d0*logT + 12d0
         write(*,*)
         write(*,1) 'frac_OPAL_SCVH', res(i_frac_OPAL_SCVH)
         write(*,1) 'frac_HELM', res(i_frac_HELM)
         write(*,1) 'frac_Skye', res(i_frac_Skye)
         write(*,1) 'frac_PC', res(i_frac_PC)
         write(*,1) 'frac_FreeEOS', res(i_frac_FreeEOS)
         write(*,1) 'frac_CMS', res(i_frac_CMS)
         write(*,*)

         if (.true.) then ! dfridr

            i_var = i_lnPgas

            dx_0 = 1d-4
            !dvardx_0 = dzbar_dx(1) !d_dxa(i_var, 1)
            dvardx_0 = d_dxa(i_var, 1)

            err = 0d0
            dvardx = dfridr(dx_0,dfridr_func,err)
            xdum = (dvardx - dvardx_0)/max(abs(dvardx_0),1d-4)
            write(*,1) 'd_dxa analytic, numeric, diff, rel diff', &
               dvardx_0, dvardx, err, xdum
            write(*,*)

         end if


         contains


         real(dp) function dfridr_func(delta_x) result(val)
            real(dp), intent(in) :: delta_x
            integer :: ierr
            real(dp) :: log_var
            real(dp) ::  xa_var(species)
            include 'formats'
            ierr = 0
            xa_var = xa
            xa_var(1) = xa_var(1) + delta_x ! h1
            !xa_var(2) = xa_var(2) - delta_x ! he4

            ! need this call to get dabar_dx, dzbar_dx
            ! might want to pass these in eventually
            call composition_info( &
               species, chem_id, xa_var, xh, xhe, zz, &
               abar_ci, zbar_ci, z2bar_ci, z53bar_ci, ye_ci, mass_correction, &
               sumx, dabar_dx, dzbar_dx, dmc_dx)
            !val = zbar_ci
            !write(*,*) abar_ci, sumx

            call eosDT_get( &
               handle, species, chem_id, net_iso, xa_var, &
               Rho, logRho, T, logT, &
               res, d_dlnd, d_dlnT, d_dxa, info)
            val = res(i_var)
         end function dfridr_func

      end subroutine do_test_eosDT_new

      
      subroutine test_HELM
         
         include 'helm_def.dek'
         real(dp) :: logT, logRho, T, Rho, X, Z, logPgas, Pgas, Prad, P, XC, XO, &
               logRho1_OPAL_SCVH_limit, logRho2_OPAL_SCVH_limit, Z_all_HELM
         real(dp) :: hres(num_helm_results), mhd_res(50)
         real(dp), dimension(num_eos_basic_results) :: &
               res, d_dlnd, d_dlnT, d_dabar, d_dzbar
         integer :: info, i_var
         real(dp) :: dx_0, var_0, dvardx_0, dvardx, xdum, err, lnd, lnT, energy, entropy, &
            var, dvar_dlnd, dvar_dlnT, dse, dpe, dsp, dS_dT, dE_dT, dP_dT, dPrad_dT, &
            dE_dRho, dS_dRho, dlnE_dlnd, dlnE_dlnT, dlnPgas_dlnT, dlnS_dlnd, dlnS_dlnT, &
            dlnPgas_dlnRho, d2P_dT_drho, d2E_dT_drho, d2S_dT_drho, dPgas_dT
         logical :: doing_d_dlnd
      
         include 'formats'
                                                   
                                              logRho = 3.2d0
                                              rho = exp10(logRho)
                                                logT = 6.7d0
                                                T = exp10(logT)
         
         open(20,file='mesa-eosDT_test.data')
         read(20,*)
         read(20,*) mhd_res(1:39)
         close(20)
         
         XC = 0
         XO = 0
         
         z = 1d0
         x = 0d0
         call Init_Composition(X, Z, XC, XO)
         write(*,*) 'change abar and zbar for test'
         abar = 16.301805415246282d0
         zbar = 8.1509027076231408d0

         write(*,*) 'call eos_get_helm_results'         
         info = 0
         call eos_get_helm_results( &
               X, abar, zbar, Rho, logRho, T, logT, &
               rq% coulomb_temp_cut_HELM, rq% coulomb_den_cut_HELM, &
               .true., .false., .true., &
               rq% logT_ion_HELM, rq% logT_neutral_HELM, hres, info)
         if (info /= 0) then
            write(*,*) 'failed in test1_eos'
            stop 1
         end if
         
         write(*,1) 'T', T
         write(*,1) 'logT', logT
         write(*,1) 'Rho', Rho
         write(*,1) 'logRho', logRho
         write(*,1) 'abar', abar
         write(*,1) 'zbar', zbar
         write(*,1) 'X', X
         write(*,1) 'Z', Z
         write(*,1) 'rq% logRho_min_OPAL_SCVH_limit', rq% logRho_min_OPAL_SCVH_limit
         write(*,*)
         
         energy = hres(h_etot)
         entropy = hres(h_stot)
         P = hres(h_ptot)
         Prad = hres(h_prad)
         Pgas = hres(h_pgas)
         dPrad_dT = hres(h_dpradt)
         
         dS_dT = hres(h_dst)
         dS_dRho = hres(h_dsd)
         dE_dT = hres(h_det)
         dE_dRho = hres(h_ded)
         dP_dT = hres(h_dpt)
         
         d2P_dT_drho = hres(h_dpdt)
         d2E_dT_drho = hres(h_dedt)
         d2S_dT_drho = hres(h_dsdt)
         
         dlnS_dlnT = dS_dT*T/entropy
         dlnS_dlnd = dS_dRho*rho/entropy
         dlnE_dlnT = dE_dT*T/energy
         dlnE_dlnd = dE_dRho*rho/energy
         dPgas_dT = dP_dT - dPrad_dT
         dlnPgas_dlnT = dPgas_dT*T/Pgas
         dlnPgas_dlnRho = hres(h_dpd)*rho/Pgas
         
         !                      value         d_dRho        d_dT            d2_dRho2         d2_dRho_dT       d2_dT2  
         write(*,1) 'H Pg', hres(h_pgas), hres(h_dpgasd), hres(h_dpgast), hres(h_dpgasdd), hres(h_dpgasdt), hres(h_dpgastt)
         write(*,1) 'H  E', hres(h_etot), hres(h_ded), hres(h_det), hres(h_dedd), hres(h_dedt), hres(h_dett)
         write(*,1) 'H  S', hres(h_stot), hres(h_dsd), hres(h_dst), hres(h_dsdd), hres(h_dsdt), hres(h_dstt)
         stop
         
         dse = T*dS_dT/dE_dT - 1.0d0
         dpe = (dE_dRho*rho*rho + T*dP_dT)/P - 1.0d0
         dsp = -dS_dRho*rho*rho/dP_dT - 1.0d0

         write(*,*)
         write(*,1) 'logT', logT, mhd_res(1)
         write(*,1) 'logRho', logRho, mhd_res(2)
         write(*,*)
         write(*,'(45x,2a26)') 'helm', 'mhd'
         write(*,1) 'logPgas', log10(Pgas), mhd_res(3)
         write(*,1) 'logE', log10(energy), mhd_res(4)
         write(*,1) 'logS', log10(entropy), mhd_res(5)
         write(*,1) 'dlnPgas_dlnT', dlnPgas_dlnT, mhd_res(6)
         write(*,1) 'dlnPgas_dlnd', dlnPgas_dlnRho, mhd_res(7)
         write(*,1) 'd2lnPgas_dlnd_dlnT', 0d0, mhd_res(8)
         write(*,1) 'dlnE_dlnT', dlnE_dlnT, mhd_res(9)
         write(*,1) 'dlnE_dlnd', dlnE_dlnd, mhd_res(10)
         write(*,1) 'd2lnE_dlnd_dlnT', 0d0, mhd_res(11)
         write(*,1) 'dlnS_dlnT', dlnS_dlnT, mhd_res(12)
         write(*,1) 'dlnS_dlnd', dlnS_dlnd, mhd_res(13)
         write(*,1) 'd2lnS_dlnd_dlnT', 0d0, mhd_res(14)
         write(*,1) 'dse', dse, mhd_res(22)
         write(*,1) 'dpe', dpe, mhd_res(23)
         write(*,1) 'dsp', dsp, mhd_res(24)
         write(*,1) 'dS_dT', dS_dT, mhd_res(25)
         write(*,1) 'dS_dRho', dS_dRho, mhd_res(26)
         write(*,1) 'dE_dT', dE_dT, mhd_res(27)
         write(*,1) 'dE_dRho', dE_dRho, mhd_res(28)
         write(*,1) 'dP_dT', dP_dT, mhd_res(29)
         write(*,1) 'dPrad_dT', dPrad_dT, mhd_res(30)
         write(*,1) 'dPgas_dT', dPgas_dT, mhd_res(31)
         write(*,1) 'Prad', Prad, mhd_res(39)
         write(*,*)
         
         stop
         
         contains
         
         real(dp) function dfridr_func(delta_x) result(val)
            real(dp), intent(in) :: delta_x
            integer :: ierr
            real(dp) :: log_var
            include 'formats'
            ierr = 0
            if (doing_d_dlnd) then
               log_var = (lnd + delta_x)/ln10
               var = exp10(log_var)
               call eos_get_helm_results( &
                  X, abar, zbar, var, log_var, T, logT, &
                  rq% coulomb_temp_cut_HELM, rq% coulomb_den_cut_HELM, &
                  .true., .false., .true., &
                  rq% logT_ion_HELM, rq% logT_neutral_HELM, hres, info)
            else
               log_var = (lnT + delta_x)/ln10
               var = exp10(log_var)
               call eos_get_helm_results( &
                  X, abar, zbar, Rho, logRho, var, log_var, &
                  rq% coulomb_temp_cut_HELM, rq% coulomb_den_cut_HELM, &
                  .true., .false., .true., &
                  rq% logT_ion_HELM, rq% logT_neutral_HELM, hres, info)
            end if
            val = res(i_var)
         end function dfridr_func

         real(dp) function dfridr(hx,err) ! from Frank
            real(dp), intent(in) :: hx
            real(dp), intent(out) :: err
            !  this routine returns the first derivative of a function func(x)
            !  at the point x, by ridders method of polynomial extrapolation.
            !  value hx is the initial step size;
            !  it should be an increment for which func changes substantially.
            !  an estimate of the error in the first derivative is returned in err.
            integer, parameter :: ntab = 20
            integer :: i,j
            real(dp) :: x,errt,fac,hh,a(ntab,ntab),xdum,ydum
            real(dp), parameter :: con2=2d0, con=sqrt(con2), big=1d50, safe=2d0
            include 'formats'
            dfridr = 0d0
            hh = hx
            ! 2nd order central difference
            a(1,1) = (dfridr_func(hh) - dfridr_func(-hh))/(2d0*hh)
            write(*,2) 'dfdx hh', 1, a(1,1), hh
            err = big
            ! succesive columns in the neville tableu will go to smaller stepsizes
            ! and higher orders of extrapolation
            do i=2,ntab
               hh = hh/con
               a(1,i) = (dfridr_func(hh) - dfridr_func(-hh))/(2d0*hh)
               !write(*,2) 'dfdx hh', i, a(1,i), hh
               ! compute extrapolations of various orders; the error stratagy is to compare
               ! each new extrapolation to one order lower but both at the same stepsize
               ! and at the previous stepsize
               fac = con2
               do j=2,i
                  a(j,i) = (a(j-1,i)*fac - a(j-1,i-1))/(fac-1d0)
                  fac = con2*fac
                  errt = max(abs(a(j,i)-a(j-1,i)),abs(a(j,i)-a(j-1,i-1)))
                  if (errt <= err) then
                     err = errt
                     dfridr = a(j,i)
                     write(*,3) 'dfridr err', i, j, dfridr, err
                  end if
               end do
               ! if higher order is worse by a significant factor safe, then bail
               if (abs(a(i,i) - a(i-1,i-1)) >= safe*err) then
                  write(*,1) 'higher order is worse', err, a(i,i), a(i-1,i-1)
                  return
               end if
            end do
         end function dfridr

      end subroutine test_HELM
      
      
      subroutine test_eosDT_partials
         real(dp) :: logT, logRho, T, Rho, X, Z, logPgas, &
               Pgas, Prad, P, lnd, dlnd, Rho2, logRho2, &
               dlnd_analytic, dlnd_numerical
         real(dp), dimension(num_eos_basic_results) :: &
               res, d_dlnd, d_dlnT, d_dabar, d_dzbar, &
               res2, d_dlnd2, d_dlnT2, d_dabar2, d_dzbar2
         integer :: info, i
      
         include 'formats'

         info = 0

         call Setup_eos
         X =    0.8d0
         Z =    0.02d0
         abar = 1.1817162154266392d0
         zbar = 1.0635445938839752d0

         logT = 6.5d0
         T = exp10(logT)
         logRho = -1d0
         Rho = exp10(logRho)
         
         call eosDT_get_legacy( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               Rho, logRho, T, logT, &
               res, d_dlnd, d_dlnT, d_dabar, d_dzbar, info)
         if (info /= 0) then
            write(*,*) 'failed in test1_eos'
            stop 1
         end if
         
         lnd = logRho*ln10
         dlnd = lnd*1d-4
         logRho2 = (lnd + dlnd)/ln10
         Rho2 = exp10(logRho2)

         call eosDT_get_legacy( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               Rho2, logRho2, T, logT, &
               res2, d_dlnd2, d_dlnT2, d_dabar2, d_dzbar2, info)
         if (info /= 0) then
            write(*,*) 'failed in test1_eos'
            stop 1
         end if
         
         d_dlnd2(:) = (res2(:) - res(:)) / dlnd
         
         do i=1,num_eos_basic_results
            dlnd_analytic = d_dlnd(i)
            dlnd_numerical = d_dlnd2(i)
            write(*,2) 'dlnd_analytic, dlnd_numerical, diff, rel diff', i, &
               dlnd_analytic, dlnd_numerical, dlnd_analytic - dlnd_numerical, &
               (dlnd_analytic - dlnd_numerical) / max(1d-99, abs(dlnd_analytic), abs(dlnd_numerical))
            write(*,*)
         end do
         
         write(*,1) 'logRho', logRho
         write(*,1) 'logRho2', logRho2
         write(*,1) 'logT', logT
         write(*,*)
         stop

      end subroutine test_eosDT_partials
      
      
      subroutine test_HELM_OPAL_transition_T
         real(dp), dimension(num_eos_basic_results) :: &
               res, d_dlnd, d_dlnT, d_dabar, d_dzbar
         real(dp) :: X, Zinit, dXO, dXC, Z, Y, T, rho, &
               Prad, Pgas, P, logT_all_HELM, logT_all_OPAL
         logical, parameter :: quietly = .true.
         integer :: ierr
         
         include 'formats'

         call Setup_eos

         ! pure Helium
         X = 0.00d0
         Zinit = 0.00d0
         dXO = 0.00d0
         dXC = 0.00d0
         
         
         Z = Zinit + dXC + dXO
         Y = 1 - (X+Z)
      
         T = 1d7; rho = 1d-3

         write(*,1) 'T', T
         write(*,1) 'rho', rho
         write(*,*)
         
         logT_all_HELM = 6d0
         logT_all_OPAL = logT_all_HELM - 0.1d0
         call do1
         
         logT_all_HELM = 7d0
         logT_all_OPAL = logT_all_HELM - 0.1d0
         call do1
         
         logT_all_HELM = 8d0
         logT_all_OPAL = logT_all_HELM - 0.1d0
         call do1
      
         stop 'test_HELM_OPAL_transition_T'
         
         contains
         
         subroutine do1
            include 'formats'
            rq% logT_all_HELM = logT_all_HELM
            rq% logT_all_OPAL = logT_all_OPAL
            write(*,1) 'logT_all_HELM', logT_all_HELM
            call Do_One_TRho(quietly,T,Rho,X,Zinit,dXC,dXO,Y,Z,res) ! opal
            Prad = crad*T*T*T*T/3
            Pgas = exp(res(i_lnPgas))
            P = Pgas + Prad
            write(*,1) 'P', P
            write(*,*)
         end subroutine do1
         
      end subroutine test_HELM_OPAL_transition_T
      
      
      
      
      
      
      subroutine Do_One(quietly)
         logical, intent(in) :: quietly
         real(dp) :: T, rho, log10_rho, log10_T
         real(dp), dimension(num_eos_basic_results) :: res, d_dlnd, d_dlnT
         integer :: info, i
         real(dp) :: XC, XO, Ah, Zh, Yh, Ahe, Zhe, Yhe, Az, Zz, Yz, ddabar, ddzbar, &
               helm_ddXh, helm_ddXz, opal_ddXh, opal_ddXz, XC0, XO0, &
               helm_P, helm_PX, helm_PZ, opal_P, opal_PX, opal_PZ, X1, X2, Z1, Z2, &
               abar1, zbar1, dXC, dlnP_dabar, dlnP_dzbar, dlnP_dXC, &
               dabar_dZ, dzbar_dZ, dlnP_dZ, P, logRhoguess
         
         if (.true.) then
            ! pure Helium
            X = 0.00d0
            Zinit = 0.00d0
            dXO = 0.00d0
            dXC = 0.00d0
            call doit('pure Helium')
            ! pure Hydrogen
            X = 1.00d0
            Zinit = 0.00d0
            dXO = 0.00d0
            dXC = 0.00d0
            call doit('pure Hydrogen')
            ! mixed Z with H&He 3:1 ratio
            Zinit = 0.03d0
            X = 0.75d0*(1 - Zinit)
            dXO = 0.00d0
            dXC = 0.00d0
            call doit('mixed Z with H&He 3:1 ratio')
            ! solar
            X = 0.70d0
            Zinit = 0.02d0
            dXO = 0.00d0
            dXC = 0.00d0
            call doit('solar')
         end if
         
         if (.true.) then ! do get_Rho and get_T
            X = 0.70d+00
            Zinit = 0.02d0
            dXO = 0.00d0
            dXC = 0.00d0
            T = exp10(4.8d0)
            rho = 1d-7
            Z = Zinit + dXC + dXO
            Y = 1 - (X+Z)
            call Init_Composition(X, Zinit, dXC, dXO)
            res(i_lnS) = log(2.9680645120000000d+09)
            call test_get_Rho_T
            if (.not. quietly) write(*,*)
         end if
         
         contains

         
         subroutine doit(str)
            character (len=*), intent(in) :: str
            
            if (.false.) then
               T = 2d8; rho = 100
               call Do_One_TRho(quietly,T,Rho,X,Zinit,dXC,dXO,Y,Z,res) ! scvh
               stop
            end if
            
            if (.not. quietly) write(*,*) trim(str)

            
            T = 1d6; rho = 1d-2
            call Do_One_TRho(quietly,T,Rho,X,Zinit,dXC,dXO,Y,Z,res) ! opal
            T = 1d4; rho = 1d-1
            call Do_One_TRho(quietly,T,Rho,X,Zinit,dXC,dXO,Y,Z,res) ! scvh
            T = 1d5; rho = 1d-1
            call Do_One_TRho(quietly,T,Rho,X,Zinit,dXC,dXO,Y,Z,res) ! opal-scvh overlap
            T = 2d8; rho = 1d2
            call Do_One_TRho(quietly,T,Rho,X,Zinit,dXC,dXO,Y,Z,res)  ! helm
            
            if (.not. quietly) write(*,*)
            
         end subroutine


         subroutine test_get_Rho_T ! using most recent values from subroutine Do_One_TRho
            real(dp) :: tol, othertol, &
               result, result_log10, log10_T, log10_rho, lnS, Prad, Pgas, logP, &
               clipped_log10rho, clipped_log10temp, &
               logRho_guess, logRho_bnd1, logRho_bnd2, other_at_bnd1, other_at_bnd2, &
               logT_guess, logT_bnd1, logT_bnd2
            integer :: i, which_other, max_iter, eos_calls, ierr
            real(dp), dimension(num_eos_basic_results) :: &
                  d_dlnd, d_dlnT, d_dabar, d_dzbar
            logical, parameter :: use_log10_for_other = .false.
            
            if (.not. quietly) write(*,*)
                        
            log10_rho = log10(rho)
            log10_T = log10(T)
            lnS = res(i_lnS)

            ierr = 0
            max_iter = 100
            tol = 1d-5
            othertol = 1d-12

 1          format(a30,1pe24.16)
            
            if (.not. quietly) then
               write(*,*)
               write(*,1) ' tolerance', tol
            end if
            if (.not. quietly) write(*,*)
            result = rho*0.5d0 ! initial guess
            result_log10 = log10(result)
            res = 0
            logRho_guess = result_log10
            logRho_bnd1 = arg_not_provided
            logRho_bnd2 = arg_not_provided
            other_at_bnd1 = arg_not_provided
            other_at_bnd2 = arg_not_provided
            call eosDT_get_Rho_legacy( &
                  handle, Z, X, abar, zbar, &
                  species, chem_id, net_iso, xa, &
                  log10_T, i_lnS, lnS, &
                  tol, othertol, max_iter, logRho_guess, &
                  logRho_bnd1, logRho_bnd2, other_at_bnd1, other_at_bnd2, &
                  result_log10, res, d_dlnd, d_dlnT, &
                  d_dabar, d_dzbar, eos_calls, ierr)
            result = exp10(result_log10)
            if (ierr /= 0) then
               write(*,*) 'ierr in test_get_Rho_T'
               stop 2
            end if
            if (.not. quietly) then
               write(*,1) 'actual logRho', log10_rho
               write(*,1) ' guess logRho', logRho_guess
               write(*,1) ' found logRho', result_log10
               write(*,1) '  wanted logS', lnS/ln10
               write(*,1) '     got logS', res(i_lnS)/ln10
               write(*,*)
            end if
         
            result = T*2d0 ! initial guess
            result_log10 = log10(result)
            res = 0
            logT_guess = result_log10
            logT_bnd1 = 3
            logT_bnd2 = 9
            call eosDT_get_T_legacy( &
                  handle, Z, X, abar, zbar, &
                  species, chem_id, net_iso, xa, &
                  log10_rho, i_lnS, lnS, &
                  tol, othertol, max_iter, logT_guess, &
                  logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
                  result_log10, res, d_dlnd, d_dlnT, &
                  d_dabar, d_dzbar, eos_calls, ierr)
            result = exp10(result_log10)
            if (ierr /= 0) then
               write(*,*) 'ierr in test_get_Rho_T'
               stop 1
            end if
            if (.not. quietly) then
               write(*,*)
               write(*,1) 'actual logT', log10_T
               write(*,1) ' guess logT', logT_guess
               write(*,1) ' found logT', result_log10
               write(*,1) '  wanted logS', lnS/ln10
               write(*,1) '     got logS', res(i_lnS)/ln10
               write(*,*)
            end if

            T = exp10(log10_T)
            Prad = crad*T*T*T*T/3
            Pgas = exp(res(i_lnPgas))
            logP = log10(Prad + Pgas)
            result = T*2d0 ! initial guess
            result_log10 = log10(result)
            res = 0
            logT_guess = result_log10
            logT_bnd1 = 3
            logT_bnd2 = 9
            call eosDT_get_T_given_Ptotal( &
                  handle, Z, X, abar, zbar, &
                  species, chem_id, net_iso, xa, &
                  log10_rho, logP, &
                  tol, othertol, max_iter, logT_guess, &
                  logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
                  result_log10, res, d_dlnd, d_dlnT, &
                  d_dabar, d_dzbar, eos_calls, ierr)
            result = exp10(result_log10)
            if (ierr /= 0) then
               write(*,*) 'ierr in test_get_Rho_T (eosDT_get_T_given_Ptotal)'
               stop 1
            end if
            if (.not. quietly) then
               write(*,*)
               write(*,1) 'actual logT', log10_T
               write(*,1) ' guess logT', logT_guess
               write(*,1) ' found logT', result_log10
               write(*,1) '  wanted logP', logP
               T = exp10(result_log10)
               Prad = crad*T*T*T*T/3
               Pgas = exp(res(i_lnPgas))
               logP = log10(Prad + Pgas)
               write(*,1) '     got logP', logP
               write(*,*)
            end if

         end subroutine test_get_Rho_T

      end subroutine Do_One
      
      
      
      
      
      subroutine Do_One_Special_FreeEOS_test(quietly)
         logical, intent(in) :: quietly
         real(dp) :: T, rho, log10_rho, log10_T
         real(dp), dimension(num_eos_basic_results) :: res, d_dlnd, d_dlnT
         integer :: info, i
         real(dp) :: XC, XO, Ah, Zh, Yh, Ahe, Zhe, Yhe, Az, Zz, Yz, ddabar, ddzbar, &
               helm_ddXh, helm_ddXz, opal_ddXh, opal_ddXz, XC0, XO0, &
               helm_P, helm_PX, helm_PZ, opal_P, opal_PX, opal_PZ, X1, X2, Z1, Z2, &
               abar1, zbar1, dXC, dlnP_dabar, dlnP_dzbar, dlnP_dXC, &
               dabar_dZ, dzbar_dZ, dlnP_dZ, P, logRhoguess
         
         if (.true.) then
            ! solar
            X = 0.70d0
            Zinit = 0.02d0
            dXO = 0.00d0
            dXC = 0.00d0
            call doit('solar')
         end if
         
         contains

         
         subroutine doit(str)
            character (len=*), intent(in) :: str
                        
            if (.not. quietly) write(*,*) trim(str)

            
            T = 1d6; rho = 1d-2
            call Do_One_TRho(quietly,T,Rho,X,Zinit,dXC,dXO,Y,Z,res) ! opal region
            T = 1d5; rho = 1d-1
            call Do_One_TRho(quietly,T,Rho,X,Zinit,dXC,dXO,Y,Z,res) ! opal-scvh overlap region
            T = 1d4; rho = 1d-1
            call Do_One_TRho(quietly,T,Rho,X,Zinit,dXC,dXO,Y,Z,res) ! scvh region
            
            if (.not. quietly) write(*,*)
            
         end subroutine


         subroutine test_get_Rho_T ! using most recent values from subroutine Do_One_TRho
            real(dp) :: tol, othertol, &
               result, result_log10, log10_T, log10_rho, lnS, Prad, Pgas, logP, &
               clipped_log10rho, clipped_log10temp, &
               logRho_guess, logRho_bnd1, logRho_bnd2, other_at_bnd1, other_at_bnd2, &
               logT_guess, logT_bnd1, logT_bnd2
            integer :: i, which_other, max_iter, eos_calls, ierr
            real(dp), dimension(num_eos_basic_results) :: &
                  d_dlnd, d_dlnT, d_dabar, d_dzbar
            logical, parameter :: use_log10_for_other = .false.
            
            if (.not. quietly) write(*,*)
                        
            log10_rho = log10(rho)
            log10_T = log10(T)
            lnS = res(i_lnS)

            ierr = 0
            max_iter = 100
            tol = 1d-5
            othertol = 1d-12

 1          format(a30,1pe24.16)
            
            if (.not. quietly) then
               write(*,*)
               write(*,1) ' tolerance', tol
            end if
            if (.not. quietly) write(*,*)
            result = rho*0.5d0 ! initial guess
            result_log10 = log10(result)
            res = 0
            logRho_guess = result_log10
            logRho_bnd1 = arg_not_provided
            logRho_bnd2 = arg_not_provided
            other_at_bnd1 = arg_not_provided
            other_at_bnd2 = arg_not_provided
            call eosDT_get_Rho_legacy( &
                  handle, Z, X, abar, zbar, &
                  species, chem_id, net_iso, xa, &
                  log10_T, i_lnS, lnS, &
                  tol, othertol, max_iter, logRho_guess, &
                  logRho_bnd1, logRho_bnd2, other_at_bnd1, other_at_bnd2, &
                  result_log10, res, d_dlnd, d_dlnT, &
                  d_dabar, d_dzbar, eos_calls, ierr)
            result = exp10(result_log10)
            if (ierr /= 0) then
               write(*,*) 'ierr in test_get_Rho_T'
               stop 2
            end if
            if (.not. quietly) then
               write(*,1) 'actual logRho', log10_rho
               write(*,1) ' guess logRho', logRho_guess
               write(*,1) ' found logRho', result_log10
               write(*,1) '  wanted logS', lnS/ln10
               write(*,1) '     got logS', res(i_lnS)/ln10
               write(*,*)
            end if
         
            result = T*2d0 ! initial guess
            result_log10 = log10(result)
            res = 0
            logT_guess = result_log10
            logT_bnd1 = 3
            logT_bnd2 = 9
            call eosDT_get_T_legacy( &
                  handle, Z, X, abar, zbar, &
                  species, chem_id, net_iso, xa, &
                  log10_rho, i_lnS, lnS, &
                  tol, othertol, max_iter, logT_guess, &
                  logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
                  result_log10, res, d_dlnd, d_dlnT, &
                  d_dabar, d_dzbar, eos_calls, ierr)
            result = exp10(result_log10)
            if (ierr /= 0) then
               write(*,*) 'ierr in test_get_Rho_T'
               stop 1
            end if
            if (.not. quietly) then
               write(*,*)
               write(*,1) 'actual logT', log10_T
               write(*,1) ' guess logT', logT_guess
               write(*,1) ' found logT', result_log10
               write(*,1) '  wanted logS', lnS/ln10
               write(*,1) '     got logS', res(i_lnS)/ln10
               write(*,*)
            end if

            T = exp10(log10_T)
            Prad = crad*T*T*T*T/3
            Pgas = exp(res(i_lnPgas))
            logP = log10(Prad + Pgas)
            result = T*2d0 ! initial guess
            result_log10 = log10(result)
            res = 0
            logT_guess = result_log10
            logT_bnd1 = 3
            logT_bnd2 = 9
            call eosDT_get_T_given_Ptotal( &
                  handle, Z, X, abar, zbar, &
                  species, chem_id, net_iso, xa, &
                  log10_rho, logP, &
                  tol, othertol, max_iter, logT_guess, &
                  logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
                  result_log10, res, d_dlnd, d_dlnT, &
                  d_dabar, d_dzbar, eos_calls, ierr)
            result = exp10(result_log10)
            if (ierr /= 0) then
               write(*,*) 'ierr in test_get_Rho_T (eosDT_get_T_given_Ptotal)'
               stop 1
            end if
            if (.not. quietly) then
               write(*,*)
               write(*,1) 'actual logT', log10_T
               write(*,1) ' guess logT', logT_guess
               write(*,1) ' found logT', result_log10
               write(*,1) '  wanted logP', logP
               T = exp10(result_log10)
               Prad = crad*T*T*T*T/3
               Pgas = exp(res(i_lnPgas))
               logP = log10(Prad + Pgas)
               write(*,1) '     got logP', logP
               write(*,*)
            end if

         end subroutine test_get_Rho_T

      end subroutine Do_One_Special_FreeEOS_test
      
      
      
      subroutine test1_eosDT_get_T
         
         real(dp) :: &
               energy, abar, zbar, X, Z, logRho, logT_tol, other_tol, other, &
               logT_guess, logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
               logT_result, new_energy, &
               res(num_eos_basic_results), d_dlnd(num_eos_basic_results), &
               d_dabar(num_eos_basic_results), d_dzbar(num_eos_basic_results), &
               d_dlnT(num_eos_basic_results)      
         integer:: ierr, which_other, eos_calls, max_iter
         logical, parameter :: use_log10_for_other = .false.
         
 1       format(a40,1pe26.16)
 
         call Setup_eos
 
         ierr = 0

                                     Z =    2.0057910422533909D-02
                                     X =    6.8067344312014655D-01
                                  abar =    1.3213531256488331D+00
                                  zbar =    1.1105746510172319D+00
                                logRho =   -7.2923785725386976D+00
                                energy =    1.9116604229463992D+12
                              logT_tol =    9.9999999999999995D-07
                             other_tol =    2.3025850929940459D-07
                            logT_guess =    3.6800855464963247D+00
                             logT_bnd1 =   -8.9999999999999999D+99
                             logT_bnd2 =   -8.9999999999999999D+99
                         other_at_bnd1 =   -8.9999999999999999D+99
                         other_at_bnd2 =   -8.9999999999999999D+99


         which_other = i_lnE
         other = log(energy)
         max_iter = 100


         write(*,1) 'logRho', logRho
         write(*,1) 'logT_guess', logT_guess
         write(*,1) 'logT_bnd1', logT_bnd1
         write(*,1) 'logT_bnd2', logT_bnd2
         write(*,1) 'energy', energy
         write(*,1) 'other', other
         write(*,1) 'Z', Z
         write(*,1) 'X', X
         write(*,1) 'abar', abar
         write(*,1) 'zbar', zbar
         write(*,1) 'logT_tol', logT_tol
         write(*,1) 'other_tol', other_tol
         write(*,*)

         call eosDT_get_T_legacy( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               logRho, which_other, other, &
               logT_tol, other_tol, max_iter, logT_guess, &
               logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
               logT_result, res, d_dlnd, d_dlnT, d_dabar, d_dzbar, &
               eos_calls, ierr)
         if (ierr /= 0) then
            write(*,*) 'ierr in test1_eosDT_get_T'
            stop 1
         end if
         write(*,*)
         !write(*,1) 'actual logT', log10_T
         write(*,1) 'guess logT', logT_guess
         write(*,1) 'found logT', logT_result
         write(*,1) 'wanted logE', other/ln10
         write(*,1) 'got logE', res(i_lnE)/ln10
         write(*,*)
         write(*,*) 'eos_calls', eos_calls
         write(*,*)
         
      end subroutine test1_eosDT_get_T
      
      
      subroutine test1_eosDT_get_rho
         
         real(dp) :: &
               P, Prad, Pgas, T, rho_guess, abar, zbar, X, Z, logT, logRho_tol, other_tol, other, &
               logRho_guess, logRho_bnd1, logRho_bnd2, other_at_bnd1, other_at_bnd2, logRho_result, &
               res(num_eos_basic_results), d_dlnd(num_eos_basic_results), Rho, logRho, &
               d_dabar(num_eos_basic_results), d_dzbar(num_eos_basic_results), &
               d_dlnT(num_eos_basic_results), &
               dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas
         integer:: ierr, which_other, eos_calls, max_iter
         logical, parameter :: use_log10_for_other = .false.
         
 1       format(a40,1pe26.16)
 
         call Setup_eos
 
         ierr = 0

         Z =     2.0062018516311897D-02
         X =     6.7967739184154219D-01
         abar =     1.3226610221256363D+00
         zbar =     1.1110240673800316D+00
         logT =     3.5999013770057040D+00
         other =     7.7701772367967239D+00
         logRho_tol =     9.9999999999999995D-07
         other_tol =     2.3025850929940460D-06
         logRho_guess =    -8.0208039694507622D+00

         logRho_bnd1 = arg_not_provided
         logRho_bnd2 = arg_not_provided
         other_at_bnd1 = arg_not_provided
         other_at_bnd2 = arg_not_provided

         which_other = i_lnPgas
         max_iter = 100
         
         write(*,1) 'logT', logT
         write(*,1) 'logRho_guess', logRho_guess
         write(*,1) 'Z', Z
         write(*,1) 'X', X
         write(*,1) 'abar', abar
         write(*,1) 'zbar', zbar
         write(*,1) 'logRho_tol', logRho_tol
         write(*,1) 'other_tol', other_tol
         write(*,*)

         call eosDT_get_Rho_legacy( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               logT, which_other, other, &
               logRho_tol, other_tol, max_iter, logRho_guess, &
               logRho_bnd1, logRho_bnd2, other_at_bnd1, other_at_bnd2, &
               logRho_result, res, d_dlnd, d_dlnT, d_dabar, d_dzbar, &
               eos_calls, ierr)
         if (ierr /= 0) then
            write(*,*) 'ierr in test1_eosDT_get_rho'
            stop 1
         end if
         write(*,*)
         write(*,1) 'guess logRho', logRho_guess
         write(*,1) 'found logRho', logRho_result
         write(*,1) 'wanted logPgas', other/ln10
         write(*,1) 'got logPgas', res(i_lnPgas)/ln10
         write(*,*)
         write(*,*) 'eos_calls', eos_calls
         write(*,*)
         
      end subroutine test1_eosDT_get_rho


      
      
      
      subroutine test1_eosPT_get_T
         
         real(dp) :: &
               energy, abar, zbar, X, Z, logPgas, logT_tol, other_tol, other, &
               logT_guess, logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
               logT_result, new_energy, &
               res(num_eos_basic_results), d_dlnd(num_eos_basic_results), &
               d_dabar(num_eos_basic_results), d_dzbar(num_eos_basic_results), &
               d_dlnT(num_eos_basic_results), &
               Rho, log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas         
         integer:: ierr, which_other, eos_calls, max_iter
         logical, parameter :: use_log10_for_other = .false.
         
 1       format(a40,1pe26.16)
 
         call Setup_eos
 
         ierr = 0
         
         write(*,*) 'test1_eosPT_get_T'

                                     Z =    0.02d0
                                     X =    0.70d0
                                  abar =    1.2966353559153956d0
                                  zbar =    1.1021400447995373d0
                                logPgas =   15d0
                                energy =    exp(3.5034294596213336d+01)
                              logT_tol =    1d-6
                             other_tol =    ln10*1d-6
                            logT_guess =    7d0
                             logT_bnd1 =   arg_not_provided
                             logT_bnd2 =   arg_not_provided
                         other_at_bnd1 =   arg_not_provided
                         other_at_bnd2 =   arg_not_provided


         which_other = i_lnE
         other = log(energy)
         max_iter = 100


         write(*,1) 'logPgas', logPgas
         write(*,1) 'logT_guess', logT_guess
         write(*,1) 'logT_bnd1', logT_bnd1
         write(*,1) 'logT_bnd2', logT_bnd2
         write(*,1) 'energy', energy
         write(*,1) 'other', other
         write(*,1) 'Z', Z
         write(*,1) 'X', X
         write(*,1) 'abar', abar
         write(*,1) 'zbar', zbar
         write(*,1) 'logT_tol', logT_tol
         write(*,1) 'other_tol', other_tol
         write(*,*)

         call eosPT_get_T( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               logPgas, which_other, other, &
               logT_tol, other_tol, max_iter, logT_guess, &
               logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
               logT_result, Rho, log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
               res, d_dlnd, d_dlnT, d_dabar, d_dzbar, &
               eos_calls, ierr)
         if (ierr /= 0) then
            write(*,*) 'ierr in test1_eosPT_get_T'
            stop 1
         end if
         write(*,*)
         write(*,1) 'guess logT', logT_guess
         write(*,1) 'found logT', logT_result
         write(*,1) 'wanted logE', other/ln10
         write(*,1) 'got logE', res(i_lnE)/ln10
         write(*,*)
         write(*,*) 'eos_calls', eos_calls
         write(*,*)
         
      end subroutine test1_eosPT_get_T
      
      
      subroutine test1_eosPT_get_Pgas
         
         real(dp) :: &
               P, Prad, T, Pgas_guess, abar, zbar, X, Z, energy, logT, logPgas_tol, other_tol, other, &
               logPgas_guess, logPgas_bnd1, logPgas_bnd2, other_at_bnd1, other_at_bnd2, logPgas_result, &
               res(num_eos_basic_results), d_dlnd(num_eos_basic_results), Pgas, logPgas, &
               d_dabar(num_eos_basic_results), d_dzbar(num_eos_basic_results), &
               d_dlnT(num_eos_basic_results), &
               Rho, log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
               dlnPgas_dlnPgas_const_T, dlnPgas_dlnT_const_Pgas
         integer:: ierr, which_other, eos_calls, max_iter
         logical, parameter :: use_log10_for_other = .false.
         
 1       format(a40,1pe26.16)
 
         call Setup_eos
 
         ierr = 0

         call Init_Composition(X, Z, 0d0, 0d0) ! sets abar and zbar
         
         write(*,*) 'test1_eosPT_get_Pgas'


                                     Z =    0.02d0
                                     X =    0.70d0
                                  abar =    1.2966353559153956d0
                                  zbar =    1.1021400447995373d0
                                  logT =    6.9d0
                                energy = 1.6413485676831915d+15
                              logPgas_tol =    1d-6
                             other_tol =    ln10*1d-6
                            logPgas_guess = 14
                             logPgas_bnd1 =   arg_not_provided
                             logPgas_bnd2 =   arg_not_provided
                         other_at_bnd1 =   arg_not_provided
                         other_at_bnd2 =   arg_not_provided


         which_other = i_lnE
         other = log(energy)
         max_iter = 100


         write(*,1) 'logT', logT
         write(*,1) 'logPgas_guess', logPgas_guess
         write(*,1) 'logPgas_bnd1', logPgas_bnd1
         write(*,1) 'logPgas_bnd2', logPgas_bnd2
         write(*,1) 'energy', energy
         write(*,1) 'other', other
         write(*,1) 'Z', Z
         write(*,1) 'X', X
         write(*,1) 'abar', abar
         write(*,1) 'zbar', zbar
         write(*,1) 'logPgas_tol', logPgas_tol
         write(*,1) 'other_tol', other_tol
         write(*,*)

         call eosPT_get_Pgas( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               logT, which_other, other, &
               logPgas_tol, other_tol, max_iter, logPgas_guess, &
               logPgas_bnd1, logPgas_bnd2, other_at_bnd1, other_at_bnd2, &
               logPgas_result, Rho, log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
               res, d_dlnd, d_dlnT, d_dabar, d_dzbar, &
               eos_calls, ierr)
         if (ierr /= 0) then
            write(*,*) 'ierr in test1_eosPT_get_Pgas'
            stop 1
         end if
         write(*,*)
         write(*,1) 'guess logPgas', logPgas_guess
         write(*,1) 'found logPgas', logPgas_result
         write(*,1) 'wanted logE', other/ln10
         write(*,1) 'got logE', res(i_lnE)/ln10
         write(*,1) 'got log10Rho', log10Rho
         write(*,*)
         write(*,*) 'eos_calls', eos_calls
         write(*,*)
         
         
      end subroutine test1_eosPT_get_Pgas


      
      
      subroutine test1_eosPT_get_Pgas_for_Rho
         
         real(dp) :: &
               P, Prad, T, Pgas_guess, abar, zbar, X, Z, energy, logT, logPgas_tol, logRho_tol, logRho_want, &
               logPgas_guess, logPgas_bnd1, logPgas_bnd2, logRho_at_bnd1, logRho_at_bnd2, logPgas_result, &
               res(num_eos_basic_results), d_dlnd(num_eos_basic_results), Pgas, logPgas, &
               d_dabar(num_eos_basic_results), d_dzbar(num_eos_basic_results), &
               d_dlnT(num_eos_basic_results), &
               Rho, logRho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
               dlnPgas_dlnPgas_const_T, dlnPgas_dlnT_const_Pgas
         integer:: ierr, which_other, eos_calls, max_iter
         logical, parameter :: use_log10_for_other = .false.
         
 1       format(a40,1pe26.16)
 
         call Setup_eos
 
         ierr = 0

         call Init_Composition(X, Z, 0d0, 0d0) ! sets abar and zbar
         
         write(*,*) 'test1_eosPT_get_Pgas_for_Rho'


                                     Z =    1.4331866202912824D-02
                                     X =    7.4801683389855433D-01
                                  abar =    1.2372753933790150D+00
                                  zbar =    1.0813484340002424D+00
                                  logT =    7.6013756122152900D+00
                              logPgas_tol =    1d-6
                              logRho_tol = 1d-6
                            logPgas_guess = 1.7508981610296928D+01
                             logPgas_bnd1 =   arg_not_provided
                             logPgas_bnd2 =   arg_not_provided
                         logRho_at_bnd1 =   arg_not_provided
                         logRho_at_bnd2 =   arg_not_provided

         max_iter = 100
         logRho_want = 1.7645028058855445D+00

         write(*,1) 'logT', logT
         write(*,1) 'logPgas_guess', logPgas_guess
         write(*,1) 'logPgas_bnd1', logPgas_bnd1
         write(*,1) 'logPgas_bnd2', logPgas_bnd2
         write(*,1) 'logRho_want', logRho_want
         write(*,1) 'Z', Z
         write(*,1) 'X', X
         write(*,1) 'abar', abar
         write(*,1) 'zbar', zbar
         write(*,1) 'logPgas_tol', logPgas_tol
         write(*,1) 'logRho_tol', logRho_tol
         write(*,*)

         call eosPT_get_Pgas_for_Rho( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               logT, logRho_want, &
               logPgas_tol, logRho_tol, max_iter, logPgas_guess, &
               logPgas_bnd1, logPgas_bnd2, logRho_at_bnd1, logRho_at_bnd2, &
               logPgas_result, Rho, logRho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
               res, d_dlnd, d_dlnT, d_dabar, d_dzbar, &
               eos_calls, ierr)
         if (ierr /= 0) then
            write(*,*) 'ierr in test1_eosPT_get_Pgas_for_Rho'
            stop 1
         end if
         write(*,*)
         write(*,1) 'guess logPgas', logPgas_guess
         write(*,1) 'found logPgas', logPgas_result
         write(*,1) 'wanted logRho', logRho_want
         write(*,1) 'got logRho', logRho
         write(*,*)
         write(*,*) 'eos_calls', eos_calls
         write(*,*)
         
         
      end subroutine test1_eosPT_get_Pgas_for_Rho


      
      
      subroutine test_components
         real(dp) :: &
            X_test, Z, XC_test, XO_test, &
            logT, logRho, Pgas, Prad, energy, entropy
         real(dp), dimension(num_eos_basic_results) :: &
            res, d_dlnd, d_dlnT
         integer:: ierr
         include 'formats'
         
         write(*,*) 'test_components'
         
         call Setup_eos
 
         ierr = 0
         
         X_test = 0.12d0
         Z = 0.03d0
         XC_test = 0d0
         XO_test = 0d0

         call Init_Composition(X_test, Z, XC_test, XO_test)

         logT = 6.0d0
         logRho = 3.5d0
         
         write(*,1) 'logT', logT
         write(*,1) 'logRho', logRho
         write(*,*)
         write(*,1) 'xa h1', xa(h1)
         write(*,1) 'xa he4', xa(he4)
         write(*,1) 'xa c12', xa(c12)
         write(*,1) 'xa n14', xa(n14)
         write(*,1) 'xa o16', xa(o16)
         write(*,1) 'xa ne20', xa(ne20)
         write(*,1) 'xa mg24', xa(mg24)
         write(*,*)
         write(*,1) 'X', X
         write(*,1) 'Y', 1d0 - (X + Z)
         write(*,1) 'Z', Z
         write(*,*)
         write(*,1) 'abar', abar
         write(*,1) 'zbar', zbar
         write(*,*)
         write(*,'(a55,9a26)') ' ', 'Pgas', 'energy', 'entropy'
         call test1(1, 'opal_scvh')
         call test1(2, 'helm')
         call test1(4, 'pc')
         write(*,*)

         contains
         
         subroutine test1(which_eos, str)
            integer, intent(in) :: which_eos
            character (len=*), intent(in) :: str
            include 'formats'
            call eosDT_test_component( &
               handle, which_eos, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               exp10(logRho), logRho, exp10(logT), logT, &
               res, d_dlnd, d_dlnT, &
               Pgas, Prad, energy, entropy, ierr)       
            if (ierr /= 0) then
               write(*,1) trim(str) // ' no results'
            else
               write(*,1) trim(str), Pgas, energy, entropy
            end if  
         end subroutine test1
         
      end subroutine test_components

      
      
      subroutine test1_eosDT_get_T_given_egas
         
         real(dp) :: &
               X, Z, abar, zbar, logRho, egas_want, egas_tol, logT_tol, logT_guess, &
               logT_bnd1, logT_bnd2, egas_at_bnd1, egas_at_bnd2, logT_result, erad, egas, energy, &
               res(num_eos_basic_results), d_dlnd(num_eos_basic_results), Pgas, logPgas, &
               d_dabar(num_eos_basic_results), d_dzbar(num_eos_basic_results), &
               d_dlnT(num_eos_basic_results)
         integer:: ierr, eos_calls, max_iter
         
 1       format(a40,1pe26.16)
 
         call Setup_eos
 
         ierr = 0
         
                                     Z =    0.7D-02
                                     X =    7.3D-01

         call Init_Composition(X, Z, 0d0, 0d0) ! sets abar and zbar
         
         write(*,*) 'test1_eosDT_get_T_given_egas'


                                  abar =    1.2559567378472252D+00
                                  zbar =    1.0864043570945732D+00
                                  logRho =   -9.4201625429594529D+00
                                  
                              egas_want = 2.0596457989663662D+12
                              egas_tol = egas_want*1d-11
                              logT_tol = 1d-11
                            logT_guess = 3.6962155439999007D+00
                            
                             logT_bnd1 =   arg_not_provided
                             logT_bnd2 =   arg_not_provided
                         egas_at_bnd1 =   arg_not_provided
                         egas_at_bnd2 =   arg_not_provided
                         

         max_iter = 100

         write(*,1) 'logRho', logRho
         write(*,1) 'logT_guess', logT_guess
         write(*,1) 'egas_want', egas_want
         write(*,1) 'Z', Z
         write(*,1) 'X', X
         write(*,1) 'abar', abar
         write(*,1) 'zbar', zbar
         write(*,1) 'logT_tol', logT_tol
         write(*,1) 'egas_tol', egas_tol
         write(*,*)

         call eosDT_get_T_given_egas( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               logRho, egas_want, &
               logT_tol, egas_tol, max_iter, logT_guess, &
               logT_bnd1, logT_bnd2, egas_at_bnd1, egas_at_bnd2, &
               logT_result, res, d_dlnd, d_dlnT, d_dabar, d_dzbar, &
               eos_calls, ierr)
         if (ierr /= 0) then
            write(*,*) 'ierr in eosDT_get_T_given_egas'
            stop 1
         end if
         energy = exp(res(i_lnE))
         erad = crad*exp10(logT_result)**4/exp10(logRho)
         egas = energy - erad
         write(*,*)
         write(*,1) 'guess logT', logT_guess
         write(*,1) 'found logT', logT_result
         write(*,1) 'wanted egas', egas_want
         write(*,1) 'got egas', egas
         write(*,1) '(want - got)/got', (egas_want - egas)/egas
         write(*,*)
         write(*,*) 'eos_calls', eos_calls
         write(*,*)
         
         
      end subroutine test1_eosDT_get_T_given_egas




         
      
      subroutine Do_One_TRho(quietly,T,Rho,X,Zinit,dXC,dXO,Y,Z,res)
         logical, intent(in) :: quietly
         real(dp), intent(in) :: T, Rho, X, Zinit, dXC, dXO
         real(dp), intent(out) :: Y, Z
         real(dp), intent(out), dimension(num_eos_basic_results) :: res

         real(dp), dimension(num_eos_basic_results) :: &
               d_dlnd, d_dlnT, d_dabar, d_dzbar
         integer :: info, i
         real(dp) :: dlnT, dlnRho, lnRho_2, Prad, Pgas, P

  101    format(a30,4x,1pe24.16)
  102    format(a30,3x,1pe24.16)
         
         
         Z = Zinit + dXC + dXO
         Y = 1 - (X+Z)
                        
         call Init_Composition(X, Zinit, dXC, dXO)
         
         if (.not. quietly) then
            write(*,*)
            write(*,*)
            write(*,102) 'X', X
            write(*,102) 'Y', Y
            write(*,102) 'Z', Z
            write(*,102) 'abar', abar
            write(*,102) 'zbar', zbar
            write(*,102) 'logRho', log10(Rho)
            write(*,102) 'logT', log10(T)
            write(*,102) 'T6', T * 1d-6
            write(*,*)
         end if
         
         call eosDT_get_legacy( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               Rho, arg_not_provided, T, arg_not_provided, &
               res, d_dlnd, d_dlnT, d_dabar, d_dzbar, info)
         if (info /= 0) then
            write(*,*) 'info', info, 'Rho', Rho, 'T', T
            write(*,*) 'failed in Do_One_TRho'
            stop 1
         end if
         
         if (.not. quietly) then
         
            write(*,*) 'eosDT_get'
            Prad = crad*T*T*T*T/3
            Pgas = exp(res(i_lnPgas))
            P = Pgas + Prad
            write(*,101) 'P', P
            write(*,101) 'E', exp(res(i_lnE))
            write(*,101) 'S', exp(res(i_lnS))
            do i = 4, 9
               write(*,101) trim(eos_names(i)), res(i)
            end do
            write(*,101) trim(eos_names(i_gamma1)), res(i_gamma1)
            write(*,101) trim(eos_names(i_gamma3)), res(i_gamma3)
            write(*,101) trim(eos_names(i_eta)), res(i_eta)
            
            if (.false.) then ! debugging
               do i = 1, num_eos_basic_results
                  write(*,101) 'd_dlnd ' // trim(eos_names(i)), d_dlnd(i)
               end do
               write(*,*)
               do i = 1, num_eos_basic_results
                  write(*,101) 'd_dlnT ' // trim(eos_names(i)), d_dlnT(i)
               end do
               write(*,*)
            end if
            
         end if

      end subroutine Do_One_TRho

      
      subroutine test_dirac_integrals
         real(dp) :: dk, T, eta, theta, fdph, fdmh, fdeta, fdtheta, theta_e
 1       format(a40,1pe26.16)
         eta = 1.46722890948893d0
         T = 11327678.5183021d0
         theta = (kerg*T)/(me*clight*clight)
         !zt   1.55623520289424d0
         theta_e = 0.929542529701454d0
         call eos_fermi_dirac_integral(-0.5d0, eta, theta, fdmh, fdeta, fdtheta)
         call eos_fermi_dirac_integral(0.5d0, eta, theta, fdph, fdeta, fdtheta)
         write(*,*)
         write(*,*) 'test_dirac_integrals'
         write(*,1) 'calculated theta_e', fdmh/fdph
         write(*,*)
         stop
      end subroutine test_dirac_integrals


      end module test_eos_support  
