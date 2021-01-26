c ***********************************************************************
!
!   Copyright (C) 2006  Bill Paxton, Frank Timmes
!
!   This file is part of MESA.
!
!   MESA is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
c ***********************************************************************

      module opal_scvh_driver
      use gauss_fermi
      implicit none

      contains
      
      
      subroutine interpolate_opal_scvh(
     >               opal_only, scvh_only, include_radiation, search_for_SCVH,
     >               logT_in,logRho_in,temp_in,den_in,abar_in,zbar_in,X_in,Z_in,
     >               logPgas, logE, logS, chiRho, chiT, 
     >               Cp, Cv, dE_dRho, dS_dT, dS_dRho, 
     >               mu, logNe, gamma1, gamma3, grad_ad, eta, dlnPgas_dlnY,
     >               data_dir,info)
      use scvh_core
      use const_def
      implicit none
      logical, intent(in) :: opal_only, scvh_only,search_for_SCVH, include_radiation
      double precision, intent(in) :: logT_in,logRho_in,temp_in,den_in,abar_in,zbar_in,X_in,Z_in
      double precision, intent(out) ::  
     >               logPgas, logE, logS, chiRho, chiT, 
     >               Cp, Cv, dE_dRho, dS_dT, dS_dRho, 
     >               mu, logNe, gamma1, gamma3, grad_ad, eta, dlnPgas_dlnY
      character (len=*), intent(in) ::  data_dir
      integer, intent(out) :: info ! returned = 0 if AOK

      double precision :: X, Y, Z, lnY, dlnY, Y0, Y1, X0, X1, logPgas0, logPgas1

      info = 0
      X = X_in
      Z = Z_in

      ! NOT DOING dlnPgas_dlnY ANY MORE.   remove it.

      call do_opal_scvh(
     >               opal_only, scvh_only,include_radiation, search_for_SCVH,
     >               logT_in,logRho_in,temp_in,den_in,abar_in,zbar_in,X,Z,
     >               logPgas, logE, logS, chiRho, chiT, 
     >               Cp, Cv, dE_dRho, dS_dT, dS_dRho, 
     >               mu, logNe, gamma1, gamma3, grad_ad, eta,
     >               data_dir,info)
      if (info /= 0) return
      
      dlnPgas_dlnY = 0

      end subroutine


      subroutine do_opal_scvh(
     >               opal_only, scvh_only, include_radiation, search_for_SCVH,
     >               logT_in,logRho_in,temp_in,den_in,abar_in,zbar_in,X_in,Z_in,
     >               logPgas, logE, logS, chiRho, chiT, 
     >               Cp, Cv, dE_dRho, dS_dT, dS_dRho, 
     >               mu, logNe, gamma1, gamma3, grad_ad, eta,
     >               data_dir,info)
      use scvh_core
      use const_def
      implicit none
      logical, intent(in) :: opal_only, scvh_only, search_for_SCVH, include_radiation
      double precision, intent(in) :: logT_in,logRho_in,temp_in,den_in,abar_in,zbar_in,X_in,Z_in
      double precision, intent(out) ::  
     >               logPgas, logE, logS, chiRho, chiT, 
     >               Cp, Cv, dE_dRho, dS_dT, dS_dRho, 
     >               mu, logNe, gamma1, gamma3, grad_ad, eta
      character (len=*), intent(in) ::  data_dir
      integer, intent(out) :: info ! returned = 0 if AOK


!..logT and logRho are log10's of temp and den


      save  ! Why do we need this?
      

!..mixes the opal and scvh equation of state
!..input:  temp,den,abar,zbar,xmassh1
!..output: the rest


      double precision :: logT,logRho,temp,den,abar,zbar,X,Z
!..for the opal portion
      double precision 
     >               logPgas_opal, logE_opal, logS_opal, chiRho_opal, chiT_opal, 
     >               Cp_opal, Cv_opal, dE_dRho_opal, dS_dT_opal, dS_dRho_opal, 
     >               mu_opal, gamma1_opal, gamma3_opal, grad_ad_opal, logNe_opal

!..for the scvh portion
      double precision 
     >               logPgas_scvh, logE_scvh, logS_scvh, chiRho_scvh, chiT_scvh, 
     >               Cp_scvh, Cv_scvh, dE_dRho_scvh, dS_dT_scvh, dS_dRho_scvh, 
     >               mu_scvh, gamma1_scvh, gamma3_scvh, grad_ad_scvh, logNe_scvh

!..mixed region
      double precision eta_ele, alfa, beta, kt, theta, fd, fdeta, fdtheta, 
     >      new_eta, ne, coef, f, eta_min, eta_max
      double precision, parameter :: eostol = 1.0d-4
      double precision, parameter :: fpmin = 1.0d-14

!..locals
      integer i
      double precision c_x, c_y, c_dx, c_dy, logRho_min, logRho_max, xnhp, xnhepp
      double precision log_free_e, log_free_e0, log_free_e1, Prad

      double precision, parameter :: tiny_eta_ele = -20d0
      double precision :: a
      
!..some constants
      double precision, parameter ::  CLN = ln10
      double precision ::  logAvo != log10(avo)
      double precision ::  mecc != me * clight * clight

!..loading the scvh tables
      integer          ifirst
      data             ifirst/0/


!..for the opal results

      integer        num_opal_results
      parameter      (num_opal_results = 11)
      double precision opal_results(num_opal_results)
      double precision t6,r,logRho_lo,logRho_hi,
     >      logRho_scvh,den_scvh,logT_scvh,temp_scvh,logQ_scvh
      logical :: scvh_okay
      integer :: irad
      double precision :: pout_opal, dpoutdd_opal, dpoutdt_opal, eout_opal, deoutdd_opal,
     >      sout_opal, dsoutdd_opal, dsoutdt_opal, pgas_out, prad_out

!..define the interior boundary between opal and scvh

      double precision, parameter :: logT1 = 6.60d0
      double precision, parameter :: logT2 = 6.50d0
      double precision, parameter :: logT3 = 4.0d0
      double precision, parameter :: logT4 = 3.4d0
      double precision, parameter :: logT5 = 3.3d0
      
      double precision, parameter :: logRho1 = 2.2d0
      double precision, parameter :: logRho2 = 1.2d0
      
      double precision, parameter :: logRho3 = -2.0d0 !-3.0d0
      double precision, parameter :: logRho4 = -3.8d0 !-3.5d0
      double precision, parameter :: logRho5 = -5.80d0
      double precision, parameter :: logRho6 = -6.8d0 !-6.3d0
      
      include 'formats'
      
      logAvo = log10(avo)
      mecc = me*clight*clight

      logT = logT_in
      logRho = logRho_in
      temp = temp_in
      den = den_in
      abar = abar_in
      zbar = zbar_in
      X = X_in
      Z = Z_in
      Prad = (crad*temp**4)/3

      if (ifirst .eq. 0) then
         call setup_scvh(data_dir)
         ifirst = 1
      end if


      if (scvh_only) then
         alfa = 0
      else if (opal_only) then
         alfa = 1
      else if (logT >= logT1) then
         alfa = 1
      else if (logT >= logT2) then
         if (logRho >= logRho1) then
            alfa = (logT-logT2)/(logT1-logT2)
         else if (logRho >= logRho2) then
            c_dy = (logT-logT2)/(logT1-logT2)
            c_dx = (logRho-logRho1)/(logRho2-logRho1)
            alfa = min(1d0,sqrt(c_dx**2 + c_dy**2))
         else ! logRho < logRho2
            alfa = 1
         end if
      else if (logT >= logT3) then
         if (logRho >= logRho1) then
            alfa = 0
         else if (logRho >= logRho4) then
            logRho_lo = logRho4 + (logT-logT3)*(logRho2-logRho4)/(logT2-logT3)
            logRho_hi = logRho3 + (logT-logT3)*(logRho1-logRho3)/(logT2-logT3)
            if (logRho >= logRho_hi) then
               alfa = 0
            else if (logRho >= logRho_lo) then
               alfa = (logRho-logRho_hi)/(logRho_lo-logRho_hi)
            else ! logRho < logRho_lo
               alfa = 1
            end if
         else ! logRho < logRho4
            alfa = 1
         end if
      else if (logT >= logT4) then
         if (logRho >= logRho3) then
            alfa = 0
         else if (logRho >= logRho6) then
            logRho_lo = logRho6 + (logRho4-logRho6)*(logT-logT4)/(logT3-logT4)
            logRho_hi = logRho5 + (logRho3-logRho5)*(logT-logT4)/(logT3-logT4)
            if (logRho >= logRho_hi) then
               alfa = 0
            else if (logRho >= logRho_lo) then
               alfa = (logRho-logRho_hi)/(logRho_lo-logRho_hi)
            else ! logRho < logRho_lo
               alfa = 1
            end if
         else ! logRho < logRho6
            alfa = 1
         end if
      else if (logT >= logT5) then
         if (logRho >= logRho5) then
            alfa = 0
         else if (logRho >= logRho6) then
            c_dy = (logT-logT4)/(logT5-logT4)
            c_dx = (logRho-logRho6)/(logRho5-logRho6)
            alfa = 1d0 - min(1d0,sqrt(c_dx**2 + c_dy**2))
         else
            alfa = (logT-logT5)/(logT4-logT5)
         end if
      else ! logT < logT5
         alfa = 0
      end if

      if (alfa < 0 .or. alfa > 1) then
         write(*,*) 'error in region code of opal_scvh_driver'
         stop 1
      end if
      alfa = 0.5d0 * (1.0d0 - cos(pi * alfa)) ! smooth the transitions
      beta = 1d0 - alfa
      
      if (alfa .ne. 0D0 .or. opal_only) then
         t6 = temp * 1d-6
         r = den
         if (include_radiation) then
            irad = 1
         else
            irad = 0
         end if
         call opal_eos (X, Z, t6, r, num_opal_results, 
     >            irad, opal_results, pgas_out, prad_out, data_dir, info)
         if (info /= 0) then
            return

            write(*,*) 'opal_eos returned info', info
            write(*,*) ' opal_only', opal_only
            write(*,*) ' alfa', alfa
            write(*,*) ' beta', beta
            write(*,*) ' X', X
            write(*,*) ' Z', Z
            write(*,*) 't6', t6
            write(*,*) ' r', r
            stop 1
         end if
         
         pout_opal    = opal_results(1) * 1d12
         dpoutdd_opal = opal_results(6) * pout_opal / den
         dpoutdt_opal = opal_results(7) * pout_opal / temp
         
         if (.false.) then
            write(*,*) 'log10(pout_opal)', log10(pout_opal)
            write(*,*) 'pout_opal', pout_opal
            write(*,*) 't6*r', t6*r
            write(*,*) 'opal_results(1)/(t6*r)', opal_results(1)/(t6*r)
            write(*,*)
         end if
   
         eout_opal    = opal_results(2) * 1d12
         ! NOTE: rather than use the opal tabulated value for cv,
         ! it is better to calculate cv from other things.
         deoutdd_opal = opal_results(4) * 1d12
   
         sout_opal    = opal_results(3) * 1d6
         dsoutdd_opal = (deoutdd_opal - pout_opal / den**2) / temp
         dsoutdt_opal = opal_results(5) * 1d6 / temp

         mu_opal = opal_results(10)
         logNe_opal = opal_results(11)

         if (include_radiation .and. Prad >= pout_opal) then
            write(*,*) 'opal_scvh_driver: Prad >= pout_opal'
            stop 1
         end if
         
         logPgas_opal = log10(pgas_out)
         logE_opal = log10(eout_opal)
         logS_opal = log10(sout_opal)
         chiRho_opal = dpoutdd_opal*den/pout_opal
         chiT_opal = dpoutdt_opal*temp/pout_opal
         gamma1_opal = opal_results(8)
         ! gam1 / (gam3 - 1) = gam2 / (gam2 - 1)
         ! gam2 / (gam2 - 1) = opal_results(9)
         gamma3_opal = 1 + opal_results(8)/opal_results(9)
         Cv_opal = chiT_opal * pout_opal / (den * temp * opal_results(8)/opal_results(9)) ! C&G 9.93
         Cp_opal = Cv_opal*gamma1_opal/chiRho_opal
         grad_ad_opal = 1/opal_results(9)
         dE_dRho_opal = deoutdd_opal
         dS_dRho_opal = dsoutdd_opal
         dS_dT_opal = dsoutdt_opal

         !write(*,1) 'opal_scvh_driver gamma1_opal', gamma1_opal


      end if

      if (beta .ne. 0D0 .or. scvh_only) then
         logRho_scvh = logRho; den_scvh = den
         logT_scvh = logT; temp_scvh= temp
         scvh_okay = .false.
         do i=1,10 ! reduce logQ if necessary
            call interpolate_scvh(
     >               include_radiation, search_for_SCVH,
     >               logT_scvh, logRho_scvh, temp_scvh, den_scvh, X,
     >               logPgas_scvh, logE_scvh, logS_scvh, chiRho_scvh, chiT_scvh, 
     >               Cp_scvh, Cv_scvh, dE_dRho_scvh, dS_dT_scvh, dS_dRho_scvh, 
     >               mu_scvh, gamma1_scvh, gamma3_scvh, grad_ad_scvh, logNe_scvh,
     >               info)
            if (info == 0) then
               scvh_okay = .true.; exit
            end if
            if (info == -2 .and. search_for_SCVH) then ! have hit edge of table.  reject it.
               ! move perpendicular to constant logQ
               logQ_scvh = logRho_scvh - 2*logT_scvh + 12
               logQ_scvh = logQ_scvh - 0.01d0 ! try this next
               logRho_scvh = logRho_scvh - 0.01d0
               den_scvh = 10**logRho_scvh
               logT_scvh = (logRho_scvh - logQ_scvh + 12)/2
               temp_scvh = 10**logT_scvh
               info = 0
               cycle
            end if
            if (info < 0) exit
         end do
         if (.not. scvh_okay) then 
            info = -1
            return

            write(*,*) 'opal_scvh_driver: failed for scvh'
            write(*,1) 'logT', logT
            write(*,1) 'logRho', logRho
            write(*,1) 'logRho_scvh', logRho_scvh
            write(*,*)

            stop 1
         end if
      end if
      
      if (beta .eq. 0d0) then ! just opal

         !write(*,*) 'just opal'
         logPgas = logPgas_opal
         logE = logE_opal
         logS = logS_opal
         chiRho = chiRho_opal
         chiT = chiT_opal
         Cp = Cp_opal
         Cv = Cv_opal
         dE_dRho = dE_dRho_opal
         dS_dT = dS_dT_opal
         dS_dRho = dS_dRho_opal
         mu = mu_opal
         gamma1 = gamma1_opal
         gamma3 = gamma3_opal
         grad_ad = grad_ad_opal
         logNe = logNe_opal

      else if (alfa .eq. 0d0) then ! just scvh
         
         !write(*,*) 'just scvh'
         logPgas = logPgas_scvh
         logE = logE_scvh
         logS = logS_scvh
         chiRho = chiRho_scvh
         chiT = chiT_scvh
         Cp = Cp_scvh
         Cv = Cv_scvh
         dE_dRho = dE_dRho_scvh
         dS_dT = dS_dT_scvh
         dS_dRho = dS_dRho_scvh
         mu = mu_scvh
         gamma1 = gamma1_scvh
         gamma3 = gamma3_scvh
         grad_ad = grad_ad_scvh
         logNe = logNe_scvh

      else ! combine alfa * opal + beta * scvh
         
         !write(*,*) 'combine alfa * opal + beta * scvh', alfa, beta
         
         if (.false.) then
            logPgas = alfa*logPgas_opal + beta*logPgas_scvh
            logE = alfa*logE_opal + beta*logE_scvh
            logS = alfa*logS_opal + beta*logS_scvh
            chiRho = alfa*chiRho_opal + beta*chiRho_scvh
            chiT = alfa*chiT_opal + beta*chiT_scvh
            Cp = alfa*Cp_opal + beta*Cp_scvh
            Cv = alfa*Cv_opal + beta*Cv_scvh
            dE_dRho = alfa*dE_dRho_opal + beta*dE_dRho_scvh
            dS_dT = alfa*dS_dT_opal + beta*dS_dT_scvh
            dS_dRho = alfa*dS_dRho_opal + beta*dS_dRho_scvh
            mu = alfa*mu_opal + beta*mu_scvh
            gamma1 = alfa*gamma1_opal + beta*gamma1_scvh
            gamma3 = alfa*gamma3_opal + beta*gamma3_scvh
            grad_ad = alfa*grad_ad_opal + beta*grad_ad_scvh
            logNe = alfa*logNe_opal + beta*logNe_scvh
         else
         
            ! new way
            call blend(
     >         alfa, beta, den, temp, Prad, 
     >         logPgas_opal, logPgas_scvh, 
     >         logS_opal, logS_scvh, dS_dT_opal, dS_dT_scvh, dS_dRho_opal, dS_dRho_scvh,
     >         chiT_opal, chiT_scvh, chiRho_opal, chiRho_scvh, mu_opal, mu_scvh, logNe_opal, logNe_scvh,
     >         logE_opal, logE_scvh, Cv_opal, Cv_scvh, dE_dRho_opal, dE_dRho_scvh,
     >         gamma1_opal, gamma1_scvh, gamma3_opal, gamma3_scvh, grad_ad_opal, grad_ad_scvh,
     >         logPgas, logE, logS, chiRho, chiT, Cp, Cv, dE_dRho, dS_dT, dS_dRho,
     >         mu, gamma1, gamma3, grad_ad, logNe)
         
            
         end if

      end if


      ! calculate eta_ele as function of Ne and T
      ! get guess to start the process

      log_free_e = logNe - logAvo - logRho
      log_free_e0 = -2d0
      log_free_e1 = -4d0
      ! when number of free electrons per nucleon gets too low, eta becomes meaningless
      ! so cut it off


      kt      = kerg * temp
      theta = kt/mecc

      if (log_free_e < log_free_e1) then

         eta_ele = tiny_eta_ele
         call dfermi(0.5d0, eta_ele, theta, fd, fdeta, fdtheta)

      else

         eta_min = -1200d0
         eta_max = 1200d0
         eta_ele = 0d0

         coef     = 4 * pi * (2 * me * kt) ** 1.5d0 / planck_h**3
         ne      = 10**logNe

         do i = 1, 100

            ! get fermi dirac integral
            call dfermi(0.5d0, eta_ele, theta, fd, fdeta, fdtheta)

            if (fdeta < 0) then
               write(*,*) fd, fdeta
               write(*,*) 'expected fdeta > 0'
               stop 1
            end if

            f = coef * fd - ne
            if (f > 0) then ! fd too large, make eta smaller to reduce it
               eta_max = eta_ele
            else  ! fd too small, make eta larger to reduce it
               eta_min = eta_ele
            end if
            new_eta = 0.5d0 * (eta_min + eta_max)
            if (abs(eta_ele - new_eta) < eostol) exit
            if (i == 100) then
               write(*,*) 'logT', logT
               write(*,*) 'logRho', logRho
               write(*,*) 'failed to find eta_ele by root solve'
               stop 1
            end if

            eta_ele = new_eta
            if (eta_ele > 1000d0) then
                eta_ele = 1000d0; exit
            end if
            if (eta_ele < -1000d0) then
                eta_ele = -1000d0; exit
            end if

         end do

         if (log_free_e < log_free_e0) then
            alfa = (log_free_e - log_free_e1) / (log_free_e0 - log_free_e1)
            alfa = 0.5d0 * (1d0 - cos(pi*alfa))
            beta = 1 - alfa
            eta_ele = alfa * eta_ele + beta * tiny_eta_ele
         end if

      end if

      eta = eta_ele

      return

      end subroutine
      
      
      subroutine blend(
     >         alfa, beta, den, temp, Prad, 
     >         logPgas_1, logPgas_2, 
     >         logS_1, logS_2, dS_dT_1, dS_dT_2, dS_dRho_1, dS_dRho_2,
     >         chiT_1, chiT_2, chiRho_1, chiRho_2, mu_1, mu_2, logNe_1, logNe_2,
     >         logE_1, logE_2, Cv_1, Cv_2, dE_dRho_1, dE_dRho_2,
     >         gamma1_1, gamma1_2, gamma3_1, gamma3_2, grad_ad_1, grad_ad_2,
     >         logPgas, logE, logS, chiRho, chiT, Cp, Cv, dE_dRho, dS_dT, dS_dRho,
     >         mu, gamma1, gamma3, grad_ad, logNe)
         double precision, intent(in) :: 
     >         alfa, beta, den, temp, Prad, 
     >         logPgas_1, logPgas_2, 
     >         logS_1, logS_2, dS_dT_1, dS_dT_2, dS_dRho_1, dS_dRho_2,
     >         chiT_1, chiT_2, chiRho_1, chiRho_2, mu_1, mu_2, logNe_1, logNe_2,
     >         logE_1, logE_2, Cv_1, Cv_2, dE_dRho_1, dE_dRho_2,
     >         gamma1_1, gamma1_2, gamma3_1, gamma3_2, grad_ad_1, grad_ad_2
         double precision, intent(out) :: 
     >         logPgas, logE, logS, chiRho, chiT, Cp, Cv, dE_dRho, dS_dT, dS_dRho,
     >         mu, gamma1, gamma3, grad_ad, logNe
         double precision :: 
     >      Pgas_1, Pgas_2, P_1, P_2, Pgas, P, s_1, s_2, s, dP_dT, dP_dRho,
     >      e_1, e_2, e, dE_dT, dP_dT_1, dP_dT_2, dP_dRho_1, dP_dRho_2, x   
         
         include 'formats'
         
         Pgas_1 = 10**logPgas_1
         Pgas_2 = 10**logPgas_2
         P_1 = Pgas_1 + Prad
         P_2 = Pgas_2 + Prad
         Pgas = alfa*Pgas_1 + beta*Pgas_2
         P = Pgas + Prad
         s_1 = 10**logS_1
         s_2 = 10**logS_2
         s = alfa*s_1 + beta*s_2
         ds_dT = alfa*dS_dT_1 + beta*dS_dT_2
         ds_dRho = alfa*dS_dRho_1 + beta*dS_dRho_2
         e_1 = 10**logE_1
         e_2 = 10**logE_2
         e = alfa*e_1 + beta*e_2
         dE_dT = alfa*Cv_1 + beta*Cv_2
         dE_dRho = alfa*dE_dRho_1 + beta*dE_dRho_2
      
         dP_dT_1 = chiT_1*P_1/temp
         dP_dT_2 = chiT_2*P_2/temp
         dP_dRho_1 = chiRho_1*P_1/den
         dP_dRho_2 = chiRho_2*P_2/den
      
         dP_dT = alfa*dP_dT_1 + beta*dP_dT_2
         dP_dRho = alfa*dP_dRho_1 + beta*dP_dRho_2
               
         logPgas = log10(Pgas)
         logE = log10(e)
         logS = log10(s)
         chiRho = dP_dRho*den/P
         chiT = dP_dT*temp/P
         Cv = dE_dT
         
         x = P*chiT/(den*temp*cv)
         
         !x = dP_dT/(dE_dT*den)
         gamma3 = x + 1
         gamma1 = chiT*x + chiRho
         grad_ad = x/gamma1
         Cp = Cv*gamma1/chiRho
         mu = alfa*mu_1 + beta*mu_2
         logNe = alfa*logNe_1 + beta*logNe_2
         
         return
         
         write(*,1) 'alfa', alfa
         write(*,1) 'beta', beta
         write(*,*)
         write(*,1) 'P', P
         write(*,1) 'P_2', P_2
         write(*,1) 'P_1', P_1
         write(*,*)
         write(*,1) 'dP_dRho', dP_dRho
         write(*,1) 'dP_dRho_2', dP_dRho_2
         write(*,1) 'dP_dRho_1', dP_dRho_1
         write(*,*)
         write(*,1) 'dP_dT', dP_dT
         write(*,1) 'dP_dT_2', dP_dT_2
         write(*,1) 'dP_dT_1', dP_dT_1
         write(*,*)
         write(*,1) 'E', E
         write(*,1) 'E_2', E_2
         write(*,1) 'E_1', E_1
         write(*,*)
         write(*,1) 'chiRho', chiRho
         write(*,1) 'chiRho_2', chiRho_2
         write(*,1) 'chiRho_1', chiRho_1
         write(*,*)
         write(*,1) 'chiT', chiT
         write(*,1) 'chiT_2', chiT_2
         write(*,1) 'chiT_1', chiT_1
         write(*,*)
         write(*,1) 'Cv', Cv
         write(*,1) 'Cv_2', Cv_2
         write(*,1) 'Cv_1', Cv_1
         write(*,*)
         write(*,1) 'gamma1', gamma1
         write(*,1) 'gamma1_2', gamma1_2
         write(*,1) 'gamma1_1', gamma1_1
         write(*,*)
         write(*,1) 'gamma3', gamma3
         write(*,1) 'gamma3_2', gamma3_2
         write(*,1) 'gamma3_1', gamma3_1
         write(*,*)
         write(*,1) 'grad_ad', grad_ad
         write(*,1) 'grad_ad_2', grad_ad_2
         write(*,1) 'grad_ad_1', grad_ad_1
         write(*,*) 'opal_scvh_driver'
         stop 1
         
      end subroutine blend
      

      end module opal_scvh_driver
