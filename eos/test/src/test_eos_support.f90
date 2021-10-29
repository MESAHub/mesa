      module test_eos_support

      use eos_def
      use eos_lib
      use chem_lib
      use chem_def
      use const_def
      use eos_support
      use math_lib
      
      implicit none

      contains
      
      
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
         Z =  0.02d0
         X =  0.6d0
         logT = 4d0
         logPgas = 5d0
         call test1_eosPT(Z, X, logPgas, logT, .false., .false., logRho, logP)
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
               d_dxa2(num_eos_d_dxa_results, species), &
               d_dlnT2(num_eos_basic_results)
         integer:: ierr, i
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
               handle, &
               species, chem_id, net_iso, xa, &
               Pgas, logPgas, T, logT, &
               Rho, logRho, dlnRho_dlnPgas_c_T, dlnRho_dlnT_c_Pgas, &
               res, d_dlnd, d_dlnT, d_dxa, ierr)
         if (ierr /= 0) then
            write(*,*) 'ierr in eosPT_get for test1_eosPT'
            call mesa_error(__FILE__,__LINE__)
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

         names = eosDT_result_names
         
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

         call eosDT_get( &
               handle, &
               species, chem_id, net_iso, xa, &
               rho, logRho, T, logT, &
               res2, d_dlnd2, d_dlnT2, d_dxa2, ierr)
         if (ierr /= 0) then
            write(*,*) 'ierr in eosDT_get for test1_eosPT'
            call mesa_error(__FILE__,__LINE__)
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
                  d_dlnd, d_dlnT
            real(dp), dimension(num_eos_d_dxa_results, species) :: &
                  d_dxa
            
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
            call eosDT_get_Rho( &
                  handle, &
                  species, chem_id, net_iso, xa, &
                  log10_T, i_lnS, lnS, &
                  tol, othertol, max_iter, logRho_guess, &
                  logRho_bnd1, logRho_bnd2, other_at_bnd1, other_at_bnd2, &
                  result_log10, res, d_dlnd, d_dlnT, &
                  d_dxa, eos_calls, ierr)
            result = exp10(result_log10)
            if (ierr /= 0) then
               write(*,*) 'ierr in test_get_Rho_T'
               call mesa_error(__FILE__,__LINE__)
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
            call eosDT_get_T( &
                  handle, &
                  species, chem_id, net_iso, xa, &
                  log10_rho, i_lnS, lnS, &
                  tol, othertol, max_iter, logT_guess, &
                  logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
                  result_log10, res, d_dlnd, d_dlnT, &
                  d_dxa, eos_calls, ierr)
            result = exp10(result_log10)
            if (ierr /= 0) then
               write(*,*) 'ierr in test_get_Rho_T'
               call mesa_error(__FILE__,__LINE__)
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
            call eosDT_get_T( &
                  handle, &
                  species, chem_id, net_iso, xa, &
                  log10_rho, i_logPtot, logP, &
                  tol, othertol, max_iter, logT_guess, &
                  logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
                  result_log10, res, d_dlnd, d_dlnT, &
                  d_dxa, eos_calls, ierr)
            result = exp10(result_log10)
            if (ierr /= 0) then
               write(*,*) 'ierr in test_get_Rho_T (eosDT_get_T_given_Ptotal)'
               call mesa_error(__FILE__,__LINE__)
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
      
      
      
      subroutine test1_eosPT_get_T
         
         real(dp) :: &
               energy, abar, zbar, X, Z, logPgas, logT_tol, other_tol, other, &
               logT_guess, logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
               logT_result, new_energy, &
               res(num_eos_basic_results), d_dlnd(num_eos_basic_results), &
               d_dxa(num_eos_d_dxa_results, species), &
               d_dlnT(num_eos_basic_results), &
               Rho, log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas         
         integer:: ierr, which_other, eos_calls, max_iter
         
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
               handle, &
               species, chem_id, net_iso, xa, &
               logPgas, which_other, other, &
               logT_tol, other_tol, max_iter, logT_guess, &
               logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
               logT_result, Rho, log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
               res, d_dlnd, d_dlnT, d_dxa, &
               eos_calls, ierr)
         if (ierr /= 0) then
            write(*,*) 'ierr in test1_eosPT_get_T'
            call mesa_error(__FILE__,__LINE__)
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
            call eosDT_get_component( &
               handle, which_eos, &
               species, chem_id, net_iso, xa, &
               exp10(logRho), logRho, exp10(logT), logT, &
               res, d_dlnd, d_dlnT, d_dxa, ierr)
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
               d_dxa(num_eos_d_dxa_results, species), &
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

         call eosDT_get_T( &
               handle, &
               species, chem_id, net_iso, xa, &
               logRho, i_egas, egas_want, &
               logT_tol, egas_tol, max_iter, logT_guess, &
               logT_bnd1, logT_bnd2, egas_at_bnd1, egas_at_bnd2, &
               logT_result, res, d_dlnd, d_dlnT, d_dxa, &
               eos_calls, ierr)
         if (ierr /= 0) then
            write(*,*) 'ierr in eosDT_get_T_given_egas'
            call mesa_error(__FILE__,__LINE__)
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
               d_dlnd, d_dlnT
         real(dp), dimension(num_eos_d_dxa_results, species) :: &
               d_dxa
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
         
         call eosDT_get( &
               handle, &
               species, chem_id, net_iso, xa, &
               Rho, arg_not_provided, T, arg_not_provided, &
               res, d_dlnd, d_dlnT, d_dxa, info)
         if (info /= 0) then
            write(*,*) 'info', info, 'Rho', Rho, 'T', T
            write(*,*) 'failed in Do_One_TRho'
            call mesa_error(__FILE__,__LINE__)
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
