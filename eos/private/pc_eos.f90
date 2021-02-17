      module pc_eos
      use utils_lib
      use math_lib
      use auto_diff
      use const_def, only: dp, PI

      implicit none
      
      contains
      
!* Equation of state for fully ionized electron-ion plasmas (EOS EIP)
! A.Y.Potekhin & G.Chabrier, Contrib. Plasma Phys., 50 (2010) 82, 
!       and references therein
! Please communicate comments/suggestions to Alexander Potekhin:
!                                            palex@astro.ioffe.ru
! Previously distributed versions (obsolete):
!   eos2000, eos2002, eos2004, eos2006, eos2007, eos2009, eos10, eos11,
!     and eos13.
! Last update: 28.05.15. All updates since 2008 are listed below.
!*   L I S T   O F   S U B R O U T I N E S :
!  MAIN (normally commented-out) - example driving routine.
!  MELANGE9 - for arbitrary ionic mixture, renders total (ion+electron)
!          pressure, internal energy, entropy, heat capacity (all
!          normalized to the ionic ideal-gas values), logarithmic
!          derivatives of pressure over temperature and density.
!  EOSFI8 - nonideal (ion-ion + ion-electron + electron-electron)
!          contributions to the free and internal energies, pressure,
!          entropy, heat capacity, derivatives of pressure over
!          logarithm of temperature and over logarithm of density (all
!          normalized to the ionic ideal-gas values) for one ionic
!          component in a mixture.
!  FITION9 - ion-ion interaction contributions to the free and internal
!          energies, pressure, entropy, heat capacity, derivatives of
!          pressure over logarithms of temperature and density.
!  FSCRliq8 - ion-electron (screening) contributions to the free and
!          internal energies, pressure, entropy, heat capacity,
!          derivatives of pressure over logarithms of temperature and
!          density in the liquid phase for one ionic component in a
!          mixture.
!  FSCRsol8 - ion-electron (screening) contributions to the free and
!          internal energies, pressure, entropy, heat capacity,
!          derivatives of pressure over logarithms of temperature and
!          density for monoionic solid.
!  FHARM12 - harmonic (including static-lattice and zero-point)
!          contributions to the free and internal energies, pressure,
!          entropy, heat capacity, derivatives of pressure over
!          logarithms of temperature and density for solid OCP.
!  HLfit12 - the same as FHARM12, but only for thermal contributions
!  ANHARM8 - anharmonic contributions to the free and internal energies,
!          pressure, entropy, heat capacity, derivatives of pressure
!          over logarithms of temperature and density for solid OCP.
!  CORMIX - correction to the linear mixing rule for the Coulomb
!          contributions to the thermodynamic functions in the liquid.
!  ELECT11 - for an ideal electron gas of arbitrary degeneracy and
!          relativity at given temperature and electron chemical
!          potential, renders number density (in atomic units), free
!          energy, pressure, internal energy, entropy, heat capacity 
!          (normalized to the electron ideal-gas values), logarithmic
!          derivatives of pressure over temperature and density.
!  EXCOR7 - electron-electron (exchange-correlation) contributions to
!          the free and internal energies, pressure, entropy, heat
!          capacity, derivatives of pressure over logarithm of
!          temperature and over logarithm of density (all normalized
!          to the classical electron ideal-gas values).
!  FERINV7 - inverse non-relativistic Fermi integrals of orders -1/2,
!          1/2, 3/2, 5/2, and their first and second derivatives.
!  BLIN9 - relativistic Fermi-Dirac integrals of orders 1/2, 3/2, 5/2,
!          and their first, second, and some third derivatives.
!  CHEMFIT7 - electron chemical potential at given density and
!          temperature, and its first derivatives over density and
!          temperature and the second derivative over temperature.
!*   I M P R O V E M E N T S   O F   2 0 0 8  -  2 0 0 9 :
!  FHARM8 uses a fit HLfit8 to the thermal free energy of the harmonic
!   Coulomb lattice, which is more accurate than its predecessor FHARM7.
!   Resulting corrections amount up to 20% for the ion heat capacity.
!   Accordingly, S/R D3fit and FthCHA7 deleted (not used anymore).
!  BLIN7 upgraded to BLIN8:
!      - cleaned (a never-reached if-else branch deleted);
!      - Sommerfeld (high-\chi) expansion improved;
!      - some third derivatives added.
!  CORMIX added (and MELANGE7 upgraded to MELANGE8 accordingly).
!  ANHARM7 upgraded to ANHARM8, more consistent with Wigner-Kirkwood.
!  Since the T- and rho-dependences of individual Z values in a mixture
!    are not considered, the corresponding inputs (AYLR, AYLT) are
!    excluded from MELANGE8 (and EOSFI7 changed to EOSFI8 accordingly).
!  ELECT7 upgraded to ELECT9 (high-degeneracy behaviour is improved)
!*   P O S T - P U B L I C A T I O N    (2 0 1 0 +)   IMPROVEMENTS :
!  ELECT9 upgraded (smooth match of two fits at chi >> 1)
!  BLIN8 replaced by BLIN9 - smooth fit interfaces at chi=0.6 and 14.
!  MELANGE8 replaced by MELANGE9 - slightly modified input/output
! 08.08.11 - corrected mistake (unlikely to have an effect) in CHEMFIT7
! 16.11.11 - ELECT9 upgraded to ELECT11 (additional output)
! 20.04.12 - FHARM8 and HLfit8 upgraded to FHARM12 and HLfit12:
!   output of HLfit12 does not include zero-point vibr., but provides U1
! 22.12.12 - MELANGE9 now includes a correction to the linear mixing
!   rule (LMR) for the Madelung energy in the random bcc multi-ion
!   lattice.
! 14.05.13 - an accidental error in programming the newly introduced
!   correction to the LMR is fixed.
! 20.05.13 - calculation of the Wigner-Kirkwood quantum diffraction term
!   for the liquid plasma is moved from EOSFI8 into MELANGE9.
! 10.12.14 - slight cleaning of the text (no effect on the results)
! 28.05.15 - an accidental error in Wigner-Kirkwood entropy correction
!   is fixed (it was in the line "Stot=Stot+FWK*DENSI" since 20.05.13).
!***********************************************************************
!                           MAIN program:               Version 02.06.09
! This driving routine allows one to compile and run this code "as is".
! In practice, however, one usually needs to link subroutines from this
! file to another (external) code, therefore the MAIN program is
! normally commented-out.
! C%C      implicit double precision (A-H), double precision (O-Z)
! C%C      parameter(MAXY=10,UN_T6=.3157746,EPS=1.d-7)
! C%C      dimension AY(MAXY),AZion(MAXY),ACMI(MAXY)
! C%C      write(*,'('' Introduce the chemical composition (up to'',I3,
! C%C     *  '' ion species):''/
! C%C     *  '' charge number Z_i, atomic weight A_i,'',
! C%C     *  '' partial number density x_i, derivatives d x_i / d ln T'',
! C%C     *  '' and d x_i / d ln rho''/
! C%C     /   '' (non-positive Z, A, or x=1 terminates the input)'')') MAXY
! C%C      NMIX=0
! C%C    3 continue
! C%C      XSUM=0.
! C%C      do IX=1,MAXY
! C%C         write(*,'(''Z, A ('',I2,''): ''$)') IX
! C%C         read*,AZion(IX),ACMI(IX)
! C%C        if (AZion(IX).le.0..or.ACMI(IX).le.0.) goto 2
! C%C         write(*,'(''x ('',I2,''): ''$)') IX
! C%C         read*,AY(IX)
! C%C         XSUM=XSUM+AY(IX)
! C%C        if (AY(IX).le.0.) goto 2
! C%C         NMIX=IX
! C%C        if (dabs(XSUM-1.d0).lt.EPS) goto 2
! C%C      enddo
! C%C    2 continue
! C%C      if (NMIX.eq.0) then
! C%C         print*,'There must be at least one set of positive (x,Z,A).'
! C%C        goto 3
! C%C      endif
! C%C      write(*,114)
! C%C      do IX=1,NMIX
! C%C         write(*,113) IX,AZion(IX),ACMI(IX),AY(IX)
! C%C      enddo
! C%C    9 continue
! C%C      write(*,'('' Input T (K) (<0 to stop): ''$)')
! C%C      read*,T
! C%C      if (T.le.0.) stop
! C%C   10 continue
! C%C      write(*,'('' Input RHO [g/cc] (<0 to new T): ''$)')
! C%C      read*,RHO
! C%C      if (RHO.le.0.) goto 9
! C%C      RHOlg=dlog10(RHO)
! C%C      Tlg=dlog10(T)
! C%C      T6=10.d0**(Tlg-6.d0)
! C%C      RHO=10.d0**RHOlg
! C%C      write(*,112)
! C%C    1 continue
! C%C      TEMP=T6/UN_T6 ! T [au]
! C%C      call MELANGE9(NMIX,AY,AZion,ACMI,RHO,TEMP,150.0,175.0, ! input
! C%C     *   PRADnkT, ! additional output - radiative pressure
! C%C     *   DENS,Zmean,CMImean,Z2mean,GAMI,CHI,TPT,LIQSOL, ! output param.
! C%C     *   PnkT,UNkT,SNk,CV,CHIR,CHIT) ! output dimensionless TD functions
! C%C      Tnk=8.31447d13/CMImean*RHO*T6 ! n_i kT [erg/cc]
! C%C      P=PnkT*Tnk/1.d12 ! P [Mbar]
! C%C      TEGRAD=CHIT/(CHIT**2+CHIR*CV/PnkT) ! from Maxwell relat.
!   --------------------   OUTPUT   --------------------------------   *
! Here in the output we have:
! RHO - mass density in g/cc
! P - total pressure in Mbar (i.e. in 1.e12 dyn/cm^2)
! PnkT=P/nkT, where n is the number density of ions, T temperature
! CV - heat capacity at constant volume, divided by number of ions, /k
! CHIT - logarithmic derivative of pressure \chi_T
! CHIR - logarithmic derivative of pressure \chi_\rho
! UNkT - internal energy divided by NkT, N being the number of ions
! SNk - entropy divided by number of ions, /k
! GAMI - ionic Coulomb coupling parameter
! TPT=T_p/T, where T_p is the ion plasma temperature
! CHI - electron chemical potential, divided by kT
! LIQSOL = 0 in the liquid state, = 1 in the solid state
! C%C      write(*,111) RHO,T6,P,PnkT,CV,CHIT,CHIR,UNkT,SNk,GAMI,TPT,CHI,
! C%C     *  LIQSOL
! C%C      goto 10
! C%C  112 format(/
! C%C     *  ' rho [g/cc]     T6 [K]      P [Mbar]   P/(n_i kT)  Cv/(N k)',
! C%C     *  '     chi_T       chi_r      U/(N k T)    S/(N k)    Gamma_i',
! C%C     *  '      T_p/T    chi_e liq/sol')
! C%C  111 format(1P,12E12.3,I2)
! C%C  113 format(I3,2F8.3,1PE12.4)
! C%C  114 format('       Z      CMI        x_j')
! C%C      end
      
      subroutine MELANGE9(NMIX,AY,AZion,ACMI,RHO,TEMP, &
         GAMIlo,GAMIhi,PRADnkT, &
         DENS,Zmean,CMImean,Z2mean,GAMImean,CHI,TPT,LIQSOL, &
         PnkT,UNkT,SNk,CV,CHIR,CHIT,ierr)
!                                                       Version 20.05.13
!                                           slight optimization 10.12.14
!          Wigner-Kirkwood correction for the entropy corrected 28.05.15
! Stems from MELANGE8 v.26.12.09.
! Difference: output PRADnkT instead of input KRAD
! + EOS of fully ionized electron-ion plasma mixture.     
! Limitations:
! (a) inapplicable in the regimes of
!      (1) bound-state formation,
!      (2) quantum liquid,
!      (3) presence of positrons;
! (b) for the case of a composition gradually depending on RHO or TEMP,
!  second-order functions (CV,CHIR,CHIT in output) should not be trusted
! Choice of the liquid or solid regime - criterion GAMI [because the
!     choice based on comparison of total (non-OCP) free energies can be
!     sometimes dangerous because of the fit uncertainties ("Local field
!     correction" in solid and quantum effects in liquid are unknown)].
! Input: NMIX - number of different elements;
!        AY - their partial number densities,
!        AZion and ACMI - their charge and mass numbers,
!        RHO - total mass density [g/cc]
!        TEMP - temperature [in a.u.=2Ryd=3.1577e5 K].
!        GAMIlo - begin mixing liquid and solid solutions
!        GAMIhi - phase transition complete into fully solid
! NB: instead of RHO, a true input is CHI, defined below
!     Hence, disagreement between RHO and DENS is the fit error (<0.4%)
! Output:
!         AY - rescaled so that to sum up to 1
!         DENS - electron number density [in a.u.=6.7483346e24 cm^{-3}]
!         Zmean=<Z>, CMImean=<A> - mean ion charge and mass numbers,
!         Z2mean=<Z^2> - mean-square ion charge number
!         GAMImean - effective ion-ion Coulomb coupling constant
!         CHI = mu_e/kT, where mu_e is the electron chem.potential
!         TPT - effective ionic quantum parameter (T_p/T)
!         LIQSOL=0/1 for liquid/solid
!         SNk - dimensionless entropy per 1 ion
!         UNkT - internal energy per kT per ion
!         PnkT - pressure / n_i kT, where n_i is the ion number density
!         PRADnkT - radiative pressure / n_i kT
!         CV - heat capacity per ion, div. by Boltzmann const.
!         CHIR - inverse compressibility -(d ln P / d ln V)_T ("\chi_r")
!         CHIT = (d ln P / d ln T)_V ("\chi_T")
      integer, intent(in) :: NMIX
      integer, intent(out) :: LIQSOL, ierr
      real(dp), intent(in) :: AZion(:), ACMI(:), GAMIlo, GAMIhi
      real(dp), intent(inout) :: AY(:)
      type(auto_diff_real_2var_order1), intent(in) :: RHO, TEMP
      type(auto_diff_real_2var_order1), intent(out) :: PRADnkT, DENS, GAMImean
      real(dp), intent(out) :: Zmean, CMImean, Z2mean
      type(auto_diff_real_2var_order1), intent(out) :: CHI, TPT, PnkT, UNkT, SNk, CV, CHIR, CHIT
      
      integer :: IX, I, J
      real(dp) :: Y, Z13, Z52, Z53, Z73, Z321, Zion, CMI
      type(auto_diff_real_2var_order1) :: UINTRAD, PRESSRAD
      type(auto_diff_real_2var_order1) :: FEid,PEid,UEid,SEid,CVE,CHITE,CHIRE
      type(auto_diff_real_2var_order1) :: DlnDH,DlnDT,DlnDHH,DlnDTT,DlnDHT
      type(auto_diff_real_2var_order1) :: DTE, PRESSE, UINTE, RS, RSI
      type(auto_diff_real_2var_order1) :: GAME, alfa, beta, TPT2
      type(auto_diff_real_2var_order1) :: UINT, PRESS, CVtot, Stot, PDLT, PDLR
      type(auto_diff_real_2var_order1) :: DENSI, PRESSI, GAMI, DNI, PRI
      type(auto_diff_real_2var_order1) :: FC1,UC1,PC1,SC1,CV1,PDT1,PDR1
      type(auto_diff_real_2var_order1) :: FC2,UC2,PC2,SC2,CV2,PDT2,PDR2
      type(auto_diff_real_2var_order1) :: FC1_0,UC1_0,PC1_0,SC1_0,CV1_0,PDT1_0,PDR1_0
      type(auto_diff_real_2var_order1) :: FC2_0,UC2_0,PC2_0,SC2_0,CV2_0,PDT2_0,PDR2_0
      type(auto_diff_real_2var_order1) :: FC1_1,UC1_1,PC1_1,SC1_1,CV1_1,PDT1_1,PDR1_1
      type(auto_diff_real_2var_order1) :: FC2_1,UC2_1,PC2_1,SC2_1,CV2_1,PDT2_1,PDR2_1
      type(auto_diff_real_2var_order1) :: FMIX,UMIX,PMIX,CVMIX,PDTMIX,PDRMIX


      type(auto_diff_real_2var_order1) :: CTP, FWK, UWK
      type(auto_diff_real_2var_order1) :: X, X1, X2, RZ, DeltaG

      real(dp), parameter :: TINY=1.d-7
      real(dp), parameter :: AUM=1822.888d0 ! a.m.u./m_e
      real(dp), parameter :: GAMIMELT=175.d0 ! OCP value of Gamma_i for melting (not used, replaced by GAMIlo/hi)
      real(dp), parameter :: RSIMELT=140.d0 ! ion density parameter of quantum melting
      real(dp), parameter :: RAD=2.5568570411948021d-07 ! Radiation constant (=4\sigma/c) (in a.u.)


      ierr = 0
      Y=0.d0
      do IX=1,NMIX
         Y=Y+AY(IX)
      enddo
      if (abs(Y-1.d0).gt.TINY) then
        do IX=1,NMIX
           AY(IX)=AY(IX)/Y
        enddo
!         print*,'MELANGE9: partial densities (and derivatives)',
!     *    ' are rescaled by factor',1./Y
      endif
      Zmean=0d0
      Z2mean=0d0
      Z52=0d0
      Z53=0d0
      Z73=0d0
      Z321=0d0 ! corr.26.12.09
      CMImean=0d0
      do IX=1,NMIX
         if (AY(IX) < TINY) cycle
         Zmean=Zmean+AY(IX)*AZion(IX)
         Z2mean=Z2mean+AY(IX)*AZion(IX)*AZion(IX)
         Z13=pow(AZion(IX),1d0/3d0)
         Z53=Z53+AY(IX)*pow5(Z13)
         Z73=Z73+AY(IX)*pow7(Z13)
         Z52=Z52+AY(IX)*pow5(sqrt(AZion(IX)))
         Z321=Z321+AY(IX)*AZion(IX)*pow3(sqrt(AZion(IX)+1.d0)) ! 26.12.09
         CMImean=CMImean+AY(IX)*ACMI(IX)
      enddo
! (0) Photons:
      UINTRAD=RAD*TEMP*TEMP*TEMP*TEMP
      PRESSRAD=UINTRAD/3d0
!      CVRAD=4.*UINTRAD/TEMP
! (1) ideal electron gas (including relativity and degeneracy)  -----  *
      DENS=RHO/11.20587d0*Zmean/CMImean ! number density of electrons [au]
      call CHEMFIT(DENS,TEMP,CHI)

! NB: CHI can be used as true input instead of RHO or DENS
      call ELECT11(TEMP,CHI, &
        DENS,FEid,PEid,UEid,SEid,CVE,CHITE,CHIRE, &
        DlnDH,DlnDT,DlnDHH,DlnDTT,DlnDHT)
! NB: at this point DENS is redefined (the difference can be ~0.1%)
      DTE=DENS*TEMP
      PRESSE=PEid*DTE ! P_e [a.u.]
      UINTE=UEid*DTE ! U_e / V [a.u.]
! (2) non-ideal Coulomb EIP  ----------------------------------------  *
      RS=pow(0.75d0/PI/DENS,1d0/3d0) ! r_s - electron density parameter
      RSI=RS*CMImean*Z73*AUM ! R_S - ion density parameter
      GAME=1d0/RS/TEMP ! electron Coulomb parameter Gamma_e
      GAMImean=Z53*GAME   ! effective Gamma_i - ion Coulomb parameter
      if (RSI.lt.RSIMELT) then ! doesn't happen in "typical" wd cases
         LIQSOL=0 ! liquid regime
      else if (GAMImean.lt.GAMIlo) then
         LIQSOL=0 ! liquid regime
      else if (GAMImean.gt.GAMIhi) then
         LIQSOL=1 ! solid regime
      else ! blend of liquid and solid
         LIQSOL=-1
         alfa = (GAMImean - GAMIlo)/(GAMIhi - GAMIlo) ! 1 for solid, 0 for liquid
         beta = 1d0 - alfa
      endif
! Calculate partial thermodynamic quantities and combine them together:
      UINT=UINTE
      PRESS=PRESSE
      CVtot=CVE*DENS
      Stot=SEid*DENS
      PDLT=PRESSE*CHITE ! d P_e[a.u.] / d ln T
      PDLR=PRESSE*CHIRE ! d P_e[a.u.] / d ln\rho
      DENSI=DENS/Zmean ! number density of all ions
      PRESSI=DENSI*TEMP ! ideal-ions total pressure (normalization)
      TPT2=0d0
      CTP=4.d0*PI/AUM/(TEMP*TEMP) ! common coefficient for TPT2.10.12.14
! Add Coulomb+xc nonideal contributions, and ideal free energy:
      do IX=1,NMIX
         if (AY(IX).lt.TINY) cycle ! skip this species
         Zion=AZion(IX)
         if (Zion.eq.0d0) cycle ! skip neutrons
         CMI=ACMI(IX)
         GAMI=pow(Zion,5d0/3d0)*GAME ! Gamma_i for given ion species
         DNI=DENSI*AY(IX) ! number density of ions of given type
         PRI=DNI*TEMP ! = ideal-ions partial pressure (normalization)         
         if (LIQSOL == 0 .or. LIQSOL == 1) then
            call EOSFI8(LIQSOL,CMI,Zion,RS,GAMI, &
               FC1,UC1,PC1,SC1,CV1,PDT1,PDR1, &
               FC2,UC2,PC2,SC2,CV2,PDT2,PDR2,ierr)
            if (ierr /= 0) return
         else
            call EOSFI8(0,CMI,Zion,RS,GAMI, &
               FC1_0,UC1_0,PC1_0,SC1_0,CV1_0,PDT1_0,PDR1_0, &
               FC2_0,UC2_0,PC2_0,SC2_0,CV2_0,PDT2_0,PDR2_0,ierr)
            if (ierr /= 0) return
            call EOSFI8(1,CMI,Zion,RS,GAMI, &
               FC1_1,UC1_1,PC1_1,SC1_1,CV1_1,PDT1_1,PDR1_1, &
               FC2_1,UC2_1,PC2_1,SC2_1,CV2_1,PDT2_1,PDR2_1,ierr)
            if (ierr /= 0) return
            FC1 = alfa*FC1_1 + beta*FC1_0
            UC1 = alfa*UC1_1 + beta*UC1_0
            PC1 = alfa*PC1_1 + beta*PC1_0
            SC1 = alfa*SC1_1 + beta*SC1_0
            CV1 = alfa*CV1_1 + beta*CV1_0
            PDT1 = alfa*PDT1_1 + beta*PDT1_0
            PDR1 = alfa*PDR1_1 + beta*PDR1_0
            FC2 = alfa*FC2_1 + beta*FC2_0
            UC2 = alfa*UC2_1 + beta*UC2_0
            PC2 = alfa*PC2_1 + beta*PC2_0
            SC2 = alfa*SC2_1 + beta*SC2_0
            CV2 = alfa*CV2_1 + beta*CV2_0
            PDT2 = alfa*PDT2_1 + beta*PDT2_0
            PDR2 = alfa*PDR2_1 + beta*PDR2_0
         end if         
! First-order TD functions:
         UINT=UINT+UC2*PRI ! internal energy density (e+i+Coul.)
         Stot=Stot+DNI*(SC2-log(AY(IX))) !entropy per unit volume[a.u.]
         PRESS=PRESS+PC2*PRI ! pressure (e+i+Coul.) [a.u.]
! Second-order functions (they take into account compositional changes):
         CVtot=CVtot+DNI*CV2 ! C_V (e+i+Coul.)/ V (optim.10.12.14)
         PDLT=PDLT+PRI*PDT2 ! d P / d ln T
         PDLR=PDLR+PRI*PDR2 ! d P / d ln\rho
         TPT2=TPT2+CTP*DNI/ACMI(IX)*AZion(IX)**2 ! opt.10.12.14
      enddo ! next IX
! Wigner-Kirkwood perturbative correction for liquid:
      TPT=sqrt(TPT2) ! effective T_p/T - ion quantum parameter
! (in the case of a mixture, this estimate is crude)
      if (LIQSOL.eq.0) then
         FWK=TPT2/24.d0 ! Wigner-Kirkwood (quantum diffr.) term
         ! MESA doesn't warn/error when this term gets large because
         ! it is not clear that we are better off falling back to HELM
         !
         ! if (FWK.gt..7d0) then
         !    print*,'MELANGE9: strong quantum effects in liquid!'
         !    ierr = -1
         !    return
         ! endif
         UWK=2.d0*FWK
         UINT=UINT+UWK*PRESSI
         Stot=Stot+FWK*DENSI ! corrected 28.05.15
         PRESS=PRESS+FWK*PRESSI
         CVtot=CVtot-UWK*DENSI ! corrected by JWS 17.04.20
         PDLT=PDLT-FWK*PRESSI
         PDLR=PDLR+UWK*PRESSI
      endif
! Corrections to the linear mixing rule:
      if (LIQSOL.eq.0) then ! liquid phase
         call CORMIX(RS,GAME,Zmean,Z2mean,Z52,Z53,Z321, &
           FMIX,UMIX,PMIX,CVMIX,PDTMIX,PDRMIX)
      else ! solid phase (only Madelung contribution) [22.12.12]
        FMIX=0d0
        do I=1,NMIX
           if (AY(I).lt.TINY) cycle
          do J=I+1,NMIX
             if (AY(J).lt.TINY) cycle
             RZ=AZion(J)/AZion(I)
             X2=AY(J)/(AY(I)+AY(J))
             X1=dim(1.d0,X2)
             if (X2.lt.TINY) cycle
             X=X2/RZ+(1.d0-1.d0/RZ)*pow(X2, RZ)
             GAMI=pow(AZion(I),5d0/3d0)*GAME ! Gamma_i corrected 14.05.13
             DeltaG=.012d0*(1.d0-1.d0/pow2(RZ))*(X1+X2*pow(RZ,5d0/3d0))
             DeltaG=DeltaG*X/X2*dim(1.d0,X)/X1
             FMIX=FMIX+AY(I)*AY(J)*GAMI*DeltaG
          enddo
        enddo
         UMIX=FMIX
         PMIX=FMIX/3.d0
         CVMIX=0d0
         PDTMIX=0d0
         PDRMIX=FMIX/2.25d0
      endif
      UINT=UINT+UMIX*PRESSI
      Stot=Stot+DENSI*(UMIX-FMIX)
      PRESS=PRESS+PMIX*PRESSI
      CVtot=CVtot+DENSI*CVMIX
      PDLT=PDLT+PRESSI*PDTMIX
      PDLR=PDLR+PRESSI*PDRMIX
! First-order:
      PRADnkT=PRESSRAD/PRESSI ! radiative pressure / n_i k T
!      CVtot=CVtot+CVRAD
!      Stot=Stot+CVRAD/3.
      PnkT=PRESS/PRESSI ! P / n_i k T
      UNkT=UINT/PRESSI ! U / N_i k T
!      UNkT=UNkT+UINTRAD/PRESSI
      SNk=Stot/DENSI ! S / N_i k
! Second-order:
      CV=CVtot/DENSI ! C_V per ion
      CHIR=PDLR/PRESS ! d ln P / d ln\rho
      CHIT=PDLT/PRESS ! d ln P / d ln T      
!      CHIT=CHIT+4.*PRESSRAD/PRESS ! d ln P / d ln T
      return
      end subroutine MELANGE9
      
      subroutine EOSFI8(LIQSOL,CMI,Zion,RS,GAMI, &
        FC1,UC1,PC1,SC1,CV1,PDT1,PDR1, &
        FC2,UC2,PC2,SC2,CV2,PDT2,PDR2,ierr)
!                                                       Version 16.09.08
!                 call FHARM8 has been replaced by call FHARM12 27.04.12
!                           Wigner-Kirkwood correction excluded 20.05.13
!                                               slight cleaning 10.12.14
! Non-ideal parts of thermodynamic functions in the fully ionized plasma
! Stems from EOSFI5 and EOSFI05 v.04.10.05
! Input: LIQSOL=0/1(liquid/solid), 
!        Zion,CMI - ion charge and mass numbers,
!        RS=r_s (electronic density parameter),
!        GAMI=Gamma_i (ion coupling),
! Output: FC1 and UC1 - non-ideal "ii+ie+ee" contribution to the 
!         free and internal energies (per ion per kT),
!         PC1 - analogous contribution to pressure divided by (n_i kT),
!         CV1 - "ii+ie+ee" heat capacity per ion [units of k]
!         PDT1=(1/n_i kT)*(d P_C/d ln T)_V
!         PDR1=(1/n_i kT)*(d P_C/d ln\rho)_T
! FC2,UC2,PC2,SC2,CV2 - analogous to FC1,UC1,PC1,SC1,CV1, but including
! the part corresponding to the ideal ion gas. This is useful for 
! preventing accuracy loss in some cases (e.g., when SC2 << SC1).
! FC2 does not take into account the entropy of mixing S_{mix}: in a
! mixture, S_{mix}/(N_i k) has to be added externally (see MELANGE9).
! FC2 does not take into account the ion spin degeneracy either.
! When needed, the spin term must be added to the entropy externally.
      integer, intent(in) :: LIQSOL
      real(dp), intent(in) :: CMI, Zion
      type(auto_diff_real_2var_order1), intent(in) :: RS, GAMI
      type(auto_diff_real_2var_order1), intent(out) :: FC1,UC1,PC1,SC1,CV1,PDT1,PDR1
      type(auto_diff_real_2var_order1), intent(out) :: FC2,UC2,PC2,SC2,CV2,PDT2,PDR2
      integer, intent(out) :: ierr

      type(auto_diff_real_2var_order1) :: GAME,FXC,UXC,PXC,CVXC,SXC,PDTXC,PDRXC
      type(auto_diff_real_2var_order1) :: TPT, COTPT, FidION
      type(auto_diff_real_2var_order1) :: FION, UION, PION, CVii, PDTii, PDRii
      type(auto_diff_real_2var_order1) :: FItot, UItot, PItot, CVItot, SCItot
      type(auto_diff_real_2var_order1) :: PDTi, PDRi
      type(auto_diff_real_2var_order1) :: Fharm,Uharm,Pharm,CVharm,Sharm,PDTharm,PDRharm
      type(auto_diff_real_2var_order1) :: Fah,Uah,Pah,CVah,PDTah,PDRah
      type(auto_diff_real_2var_order1) :: FSCR,USCR,PSCR,S_SCR,CVSCR,PDTSCR,PDRSCR
      type(auto_diff_real_2var_order1) :: FC0, UC0, PC0, SC0, CV0, PDT0, PDR0
      
      real(dp), parameter :: TINY=1.d-20 
      real(dp), parameter :: AUM=1822.888d0 ! a.m.u/m_e

      ierr = 0
      if (LIQSOL.ne.1.and.LIQSOL.ne.0) then
         ierr = -1
         return
         !stop 'EOSFI8: invalid LIQSOL'
      end if
      if (CMI.le..1d0)  then
         ierr = -1
         return
         !stop 'EOSFI8: too small CMI'
      end if
      if (Zion.le..1d0)  then
         ierr = -1
         return
         !stop 'EOSFI8: too small Zion'
      end if
      if (RS.le..0d0)  then
         ierr = -1
         return
         !stop 'EOSFI8: invalid RS'
      end if
      if (GAMI.le..0d0)  then
         ierr = -1
         return
         !stop 'EOSFI8: invalid GAMI'
      end if
      GAME=GAMI/pow(Zion,5d0/3d0)
      call EXCOR7(RS,GAME,FXC,UXC,PXC,CVXC,SXC,PDTXC,PDRXC) ! "ee"("xc")
! Calculate "ii" part:
      COTPT=sqrt(3d0/AUM/CMI)/pow(Zion,7d0/6d0) ! auxiliary coefficient
      TPT=GAMI/sqrt(RS)*COTPT              ! T_p/T
      FidION=1.5d0*log(TPT*TPT/GAMI)-1.323515d0
! 1.3235=1+0.5*ln(6/pi); FidION = F_{id.ion gas}/(N_i kT), but without
! the term x_i ln x_i = -S_{mix}/(N_i k).
      if (LIQSOL.eq.0) then                 ! liquid
         call FITION9(GAMI, &
           FION,UION,PION,CVii,PDTii,PDRii)
         FItot=FION+FidION
         UItot=UION+1.5d0
         PItot=PION+1.0d0
         CVItot=CVii+1.5d0
         SCItot=UItot-FItot
         PDTi=PDTii+1.d0
         PDRi=PDRii+1.d0
      else                                  ! solid
         call FHARM12(GAMI,TPT, &
           Fharm,Uharm,Pharm,CVharm,Sharm,PDTharm,PDRharm) ! harm."ii"
         call ANHARM8(GAMI,TPT,Fah,Uah,Pah,CVah,PDTah,PDRah) ! anharm.
         FItot=Fharm+Fah
         FION=FItot-FidION
         UItot=Uharm+Uah
         UION=UItot-1.5d0 ! minus 1.5=ideal-gas, in order to get "ii"
         PItot=Pharm+Pah
         PION=PItot-1.d0 ! minus 1=ideal-gas
         PDTi=PDTharm+PDTah
         PDRi=PDRharm+PDRah
         PDTii=PDTi-1.d0 ! minus 1=ideal-gas
         PDRii=PDRi-1.d0 ! minus 1=ideal-gas
         CVItot=CVharm+CVah
         SCItot=Sharm+Uah-Fah
         CVii=CVItot-1.5d0 ! minus 1.5=ideal-gas
      endif
! Calculate "ie" part:
      if (LIQSOL.eq.1) then
         call FSCRsol8(RS,GAMI,Zion,TPT, &
           FSCR,USCR,PSCR,S_SCR,CVSCR,PDTSCR,PDRSCR,ierr)
         if (ierr /= 0) return
      else
         call FSCRliq8(RS,GAME,Zion, &
           FSCR,USCR,PSCR,CVSCR,PDTSCR,PDRSCR,ierr)
         if (ierr /= 0) return
         S_SCR=USCR-FSCR
      endif
! Total excess quantities ("ii"+"ie"+"ee", per ion):
      FC0=FSCR+Zion*FXC
      UC0=USCR+Zion*UXC
      PC0=PSCR+Zion*PXC
      SC0=S_SCR+Zion*SXC
      CV0=CVSCR+Zion*CVXC
      PDT0=PDTSCR+Zion*PDTXC
      PDR0=PDRSCR+Zion*PDRXC
      FC1=FION+FC0
      UC1=UION+UC0
      PC1=PION+PC0
      SC1=(UION-FION)+SC0
      CV1=CVii+CV0
      PDT1=PDTii+PDT0
      PDR1=PDRii+PDR0
! Total excess + ideal-ion quantities
      FC2=FItot+FC0
      UC2=UItot+UC0
      PC2=PItot+PC0
      SC2=SCItot+SC0
      CV2=CVItot+CV0
      PDT2=PDTi+PDT0      
      PDR2=PDRi+PDR0
      return
      end subroutine EOSFI8

! ==================  ELECTRON-ION COULOMB LIQUID  =================== *
      subroutine FITION9(GAMI, &
           FION,UION,PION,CVii,PDTii,PDRii)
!                                                       Version 11.09.08
! Non-ideal contributions to thermodynamic functions of classical OCP,
!       corrected at small density for a mixture.
!   Stems from FITION00 v.24.05.00.
! Input: GAMI - ion coupling parameter
! Output: FION - ii free energy / N_i kT
!         UION - ii internal energy / N_i kT
!         PION - ii pressure / n_i kT
!         CVii - ii heat capacity / N_i k
!         PDTii = PION + d(PION)/d ln T = (1/N_i kT)*(d P_{ii}/d ln T)
!         PDRii = PION + d(PION)/d ln\rho
!   Parameters adjusted to Caillol (1999).
      use pc_support
      type(auto_diff_real_2var_order1), intent(in) :: GAMI
      type(auto_diff_real_2var_order1), intent(out) :: FION, UION, PION, CVii, PDTii, PDRii
      type(auto_diff_real_2var_order1) :: F0, U0

      real(dp), parameter :: A1=-.907347d0
      real(dp), parameter :: A2=.62849d0
      real(dp) :: A3
      real(dp), parameter :: C1=.004500d0
      real(dp), parameter :: G1=170.0d0
      real(dp), parameter :: C2=-8.4d-5
      real(dp), parameter :: G2=.0037d0
      real(dp), parameter :: SQ32=.8660254038d0 ! SQ32=sqrt(3)/2
      
      real(dp) :: &
         xFION, dFION_dlnGAMI, &
         xUION, dUION_dlnGAMI, &
         xPION, dPION_dlnGAMI, &
         xCVii, dCVii_dlnGAMI, &
         xPDTii, dPDTii_dlnGAMI, &
         xPDRii, dPDRii_dlnGAMI
      real(dp) :: dlnGAMI_dT, dlnGAMI_dRho
      integer :: ierr
      logical :: skip
      logical, parameter :: use_FITION9_table = .false.
      logical, parameter :: debug_FITION9_table = .false.


      if (.not. use_FITION9_table .or. debug_FITION9_table) then

         A3=-SQ32-A1/sqrt(A2)
         F0=A1*(sqrt(GAMI*(A2+GAMI)) &
            - A2*log(sqrt(GAMI/A2)+sqrt(1d0+GAMI/A2))) &
            + 2d0*A3*(sqrt(GAMI)-atan(sqrt(GAMI)))
         U0=pow3(sqrt(GAMI))*(A1/sqrt(A2+GAMI)+A3/(1.d0+GAMI))
         !   This is the zeroth approximation. Correction:
         UION=U0+C1*GAMI*GAMI/(G1+GAMI)+C2*GAMI*GAMI/(G2+GAMI*GAMI)
         FION=F0+C1*(GAMI-G1*log(1.d0+GAMI/G1)) &
            + C2/2d0*log(1.d0+GAMI*GAMI/G2)
         CVii=-0.5d0*pow3(sqrt(GAMI))*(A1*A2/pow3(sqrt(A2+GAMI)) &
            + A3*(1.d0-GAMI)/pow2(1.d0+GAMI)) &
            - GAMI*GAMI*(C1*G1/pow2(G1+GAMI)+C2*(G2-GAMI*GAMI)/pow2(G2+GAMI*GAMI))
         PION=UION/3.0d0
         PDRii=(4.0d0*UION-CVii)/9.0d0 ! p_{ii} + d p_{ii} / d ln\rho
         PDTii=CVii/3.0d0 ! p_{ii} + d p_{ii} / d ln T

      endif
      
      if (use_FITION9_table .or. debug_FITION9_table) then
         ierr = 0
         call get_FITION9(GAMI%val, &
            xFION, dFION_dlnGAMI, &
            xUION, dUION_dlnGAMI, &
            xPION, dPION_dlnGAMI, &
            xCVii, dCVii_dlnGAMI, &
            xPDTii, dPDTii_dlnGAMI, &
            xPDRii, dPDRii_dlnGAMI, &
            skip, ierr)
         if (ierr /= 0) stop 'failed in call get_FITION9'
      else
         skip = .true.
      end if


      if (debug_FITION9_table .and. .not. skip) then

         if (.not. check1(FION, xFION, 'FION')) return
         if (.not. check1(UION, xUION, 'UION')) return
         if (.not. check1(PION, xPION, 'PION')) return
         if (.not. check1(CVii, xCVii, 'CVii')) return
         if (.not. check1(PDTii, xPDTii, 'PDTii')) return
         if (.not. check1(PDRii, xPDRii, 'PDRii')) return

      endif


      if (use_FITION9_table .and. .not. skip) then

         ! values
         FION% val = xFION
         UION% val = xUION
         PION% val = xPION
         CVii% val = xCVii
         PDTii% val = xPDTii
         PDRii% val = xPDRii

         ! dT (val1) derivatives
         ! dlnRs_dT = RS% d1val1 / RS% val = 0
         dlnGAMI_dT = GAMI% d1val1 / GAMI% val

         FION% d1val1 = dFION_dlnGAMI * dlnGAMI_dT
         UION% d1val1 = dUION_dlnGAMI * dlnGAMI_dT
         PION% d1val1 = dPION_dlnGAMI * dlnGAMI_dT
         CVii% d1val1 = dCVii_dlnGAMI * dlnGAMI_dT
         PDTii% d1val1 = dPDTii_dlnGAMI * dlnGAMI_dT
         PDRii% d1val1 = dPDRii_dlnGAMI * dlnGAMI_dT

         ! dRho (val2) derivatives
         dlnGAMI_dRho = GAMI% d1val2 / GAMI% val

         FION% d1val2 = dFION_dlnGAMI * dlnGAMI_dRho
         UION% d1val2 = dUION_dlnGAMI * dlnGAMI_dRho
         PION% d1val2 = dPION_dlnGAMI * dlnGAMI_dRho
         CVii% d1val2 = dCVii_dlnGAMI * dlnGAMI_dRho
         PDTii% d1val2 = dPDTii_dlnGAMI * dlnGAMI_dRho
         PDRii% d1val2 = dPDRii_dlnGAMI * dlnGAMI_dRho

      end if


      contains
      
      logical function check1(v, xv, str)
         type(auto_diff_real_2var_order1), intent(in) :: v
         real(dp), intent(in) :: xv
         character (len=*), intent(in) :: str
         real(dp) :: val
         real(dp), parameter :: atol = 1d-8, rtol = 1d-6
         include 'formats'
         val = v%val
         check1 = .false.
         if (is_bad(xv)) then
            write(*,*) 'is_bad ' // trim(str), xv, val, GAMI
            return
         end if
         if (abs(val - xv) > atol + rtol*max(abs(val),abs(xv))) then
            write(*,*) 'rel mismatch ' // trim(str), &
               (val - xv)/max(abs(val),abs(xv),1d-99), &
               xv, val, GAMI%val
            stop 'FITION9'
            return
         end if
         check1 = .true.
      end function check1

      end subroutine FITION9

      subroutine FSCRliq8(RS,GAME,Zion, &
           FSCR,USCR,PSCR,CVSCR,PDTSCR,PDRSCR,ierr) ! fit to the el.-ion scr.
!                                                       Version 11.09.08
!                                                       cleaned 16.06.09
! Stems from FSCRliq7 v. 09.06.07. Included a check for RS=0.
!   INPUT: RS - density parameter, GAME - electron Coulomb parameter,
!          Zion - ion charge number,
!   OUTPUT: FSCR - screening (e-i) free energy per kT per 1 ion,
!           USCR - internal energy per kT per 1 ion (screen.contrib.)
!           PSCR - pressure divided by (n_i kT) (screen.contrib.)
!           CVSCR - heat capacity per 1 ion (screen.contrib.)
!           PDTSCR,PDRSCR = PSCR + d PSCR / d ln(T,\rho)
      use pc_support
      type(auto_diff_real_2var_order1), intent(in) :: RS,GAME
      real(dp), intent(in) :: Zion
      type(auto_diff_real_2var_order1), intent(out) :: FSCR,USCR,PSCR,CVSCR,PDTSCR,PDRSCR
      integer, intent(out) :: ierr

      type(auto_diff_real_2var_order1) :: SQG, SQR, SQZ1, SQZ, CDH0, CDH, ZLN, Z13
      type(auto_diff_real_2var_order1) :: X, CTF, P01, P03, PTX
      type(auto_diff_real_2var_order1) :: TX, TXDG, TXDGG, TY1, TY1DG, TY1DGG
      type(auto_diff_real_2var_order1) :: TY2, TY2DX, TY2DXX, TY, TYX, TYDX, TYDG, P1
      type(auto_diff_real_2var_order1) :: COR1, COR1DX, COR1DG, COR1DXX, COR1DGG, COR1DXG
      type(auto_diff_real_2var_order1) :: U0, U0DX, U0DG, U0DXX, U0DGG, U0DXG
      type(auto_diff_real_2var_order1) :: D0DG, D0, D0DX, D0DXX
      type(auto_diff_real_2var_order1) :: COR0, COR0DX, COR0DG, COR0DXX, COR0DGG, COR0DXG
      type(auto_diff_real_2var_order1) :: RELE, Q1, Q2, H1U, H1D, H1, H1X, H1DX, H1DXX
      type(auto_diff_real_2var_order1) :: UP, UPDX, UPDG, UPDXX, UPDGG, UPDXG
      type(auto_diff_real_2var_order1) :: DN1, DN1DX, DN1DG, DN1DXX, DN1DGG, DN1DXG
      type(auto_diff_real_2var_order1) :: DN, DNDX, DNDG, DNDXX, DNDGG, DNDXG
      type(auto_diff_real_2var_order1) :: FX, FXDG, FDX, FG, FDG, FDGDH, FDXX, FDGG, FDXG
      
      real(dp), parameter :: XRS=.0140047d0
      real(dp), parameter :: TINY=1.d-19
      
      real(dp) :: &
         xFSCR, dFSCR_dlnRS, dFSCR_dlnGAME, &
         xUSCR, dUSCR_dlnRS, dUSCR_dlnGAME, &
         xPSCR, dPSCR_dlnRS, dPSCR_dlnGAME, &
         xCVSCR, dCVSCR_dlnRS, dCVSCR_dlnGAME, &
         xPDTSCR, dPDTSCR_dlnRS, dPDTSCR_dlnGAME, &
         xPDRSCR, dPDRSCR_dlnRS, dPDRSCR_dlnGAME
      real(dp) :: dlnRs_dT, dlnRs_dRho, dlnGAME_dT, dlnGAME_dRho
      logical :: skip
      logical, parameter :: use_FSCRliq8_table = .false.
      logical, parameter :: debug_FSCRliq8_table = .false.

      ierr = 0
      if (RS.lt.0d0) then
         ierr = -1
         return
         !stop 'FSCRliq8: RS < 0'
      end if
      if (RS.lt.TINY) then
         FSCR=0.d0
         USCR=0.d0
         PSCR=0.d0
         CVSCR=0.d0
         PDTSCR=0.d0
         PDRSCR=0.d0
         return
      endif
      
      if (use_FSCRliq8_table .or. debug_FSCRliq8_table) then
         ierr = 0
         call get_FSCRliq8(int(Zion), RS%val, GAME%val, &
            xFSCR, dFSCR_dlnRS, dFSCR_dlnGAME, &
            xUSCR, dUSCR_dlnRS, dUSCR_dlnGAME, &
            xPSCR, dPSCR_dlnRS, dPSCR_dlnGAME, &
            xCVSCR, dCVSCR_dlnRS, dCVSCR_dlnGAME, &
            xPDTSCR, dPDTSCR_dlnRS, dPDTSCR_dlnGAME, &
            xPDRSCR, dPDRSCR_dlnRS, dPDRSCR_dlnGAME, skip, ierr)
         if (ierr /= 0) return      
      else
         skip = .true.
      endif

      if (.not. use_FSCRliq8_table .or. debug_FSCRliq8_table .or. skip) then

         SQG=sqrt(GAME)
         SQR=sqrt(RS)
         SQZ1=sqrt(1d0+Zion)
         SQZ=sqrt(Zion)
         CDH0=Zion/1.73205d0 ! 1.73205=sqrt(3.)
         CDH=CDH0*(SQZ1*SQZ1*SQZ1-SQZ*SQZ*SQZ-1d0)
         SQG=sqrt(GAME)
         ZLN=log(Zion)
         Z13=exp(ZLN/3.d0) ! Zion**(1./3.)
         X=XRS/RS ! relativity parameter
         CTF=Zion*Zion*.2513d0*(Z13-1d0+.2d0/sqrt(Z13))
         ! Thomas-Fermi constant; .2513=(18/175)(12/\pi)^{2/3}
         P01=1.11d0*exp(0.475d0*ZLN)
         P03=0.2d0+0.078d0*ZLN*ZLN
         PTX=1.16d0+0.08d0*ZLN
         TX=pow(GAME,PTX)
         TXDG=PTX*TX/GAME
         TXDGG=(PTX-1.d0)*TXDG/GAME
         TY1=1d0/(1.d-3*Zion*Zion+2d0*GAME)
         TY1DG=-2d0*TY1*TY1
         TY1DGG=-4d0*TY1*TY1DG
         TY2=1d0+6d0*RS*RS
         TY2DX=-12d0*RS*RS/X
         TY2DXX=-3d0*TY2DX/X
         TY=RS*RS*RS/TY2*(1d0+TY1)
         TYX=3d0/X+TY2DX/TY2
         TYDX=-TY*TYX
         TYDG=RS*RS*RS*TY1DG/TY2
         P1=(Zion-1d0)/9d0
         COR1=1d0+P1*TY
         COR1DX=P1*TYDX
         COR1DG=P1*TYDG
         COR1DXX=P1*(TY*(3d0/(X*X)+(TY2DX/TY2)*(TY2DX/TY2)-TY2DXX/TY2)-TYDX*TYX)
         COR1DGG=P1*RS*RS*RS*TY1DGG/TY2
         COR1DXG=-P1*TYDG*TYX
         U0=0.78d0*sqrt(GAME/Zion)*RS*RS*RS
         U0DX=-3d0*U0/X
         U0DG=.5d0*U0/GAME
         U0DXX=-4.d0*U0DX/X
         U0DGG=-.5d0*U0DG/GAME
         U0DXG=-3.d0*U0DG/X
         D0DG=Zion*Zion*Zion
         D0=GAME*D0DG+21d0*RS*RS*RS
         D0DX=-63d0*RS*RS*RS/X
         D0DXX=252d0*RS*RS*RS/(X*X)
         COR0=1d0+U0/D0
         COR0DX=(U0DX-U0*D0DX/D0)/D0
         COR0DG=(U0DG-U0*D0DG/D0)/D0
         COR0DXX=(U0DXX-(2d0*U0DX*D0DX+U0*D0DXX)/D0+2d0*(D0DX/D0)*(D0DX/D0))/D0
         COR0DGG=(U0DGG-2d0*U0DG*D0DG/D0+2d0*U0*(D0DG/D0)*(D0DG/D0))/D0
         COR0DXG=(U0DXG-(U0DX*D0DG+U0DG*D0DX)/D0+2d0*U0*D0DX*D0DG/(D0*D0))/D0
         ! Relativism:
         RELE=sqrt(1.d0+X*X)
         Q1=0.18d0/sqrt(sqrt(Zion))
         Q2=0.2d0+0.37d0/sqrt(Zion)
         H1U=1d0+X*X/5.d0
         H1D=1d0+Q1*X+Q2*X*X
         H1=H1U/H1D
         H1X=0.4d0*X/H1U-(Q1+2d0*Q2*X)/H1D
         H1DX=H1*H1X
         H1DXX=H1DX*H1X &
            + H1*(0.4d0/H1U-(0.4d0*X/H1U)*(0.4d0*X/H1U)-2d0*Q2/H1D+pow2((Q1+2d0*Q2*X)/H1D))
         UP=CDH*SQG+P01*CTF*TX*COR0*H1
         UPDX=P01*CTF*TX*(COR0DX*H1+COR0*H1DX)
         UPDG=.5d0*CDH/SQG+P01*CTF*(TXDG*COR0+TX*COR0DG)*H1
         UPDXX=P01*CTF*TX*(COR0DXX*H1+2d0*COR0DX*H1DX+COR0*H1DXX)
         UPDGG=-.25d0*CDH/(SQG*GAME) &
            + P01*CTF*(TXDGG*COR0+2d0*TXDG*COR0DG+TX*COR0DGG)*H1
         UPDXG=P01*CTF*(TXDG*(COR0DX*H1+COR0*H1DX) &
            + TX*(COR0DXG*H1+COR0DG*H1DX))
         DN1=P03*SQG+P01/RS*TX*COR1
         DN1DX=P01*TX*(COR1/XRS+COR1DX/RS)
         DN1DG=.5d0*P03/SQG+P01/RS*(TXDG*COR1+TX*COR1DG)
         DN1DXX=P01*TX/XRS*(2d0*COR1DX+X*COR1DXX)
         DN1DGG=-.25d0*P03/(GAME*SQG) &
            + P01/RS*(TXDGG*COR1+2d0*TXDG*COR1DG+TX*COR1DGG)
         DN1DXG=P01*(TXDG*(COR1/XRS+COR1DX/RS)+TX*(COR1DG/XRS+COR1DXG/RS))
         DN=1d0+DN1/RELE
         DNDX=DN1DX/RELE-X*DN1/(RELE*RELE*RELE)
         DNDXX=(DN1DXX-((2d0*X*DN1DX+DN1)-3.d0*X*X*DN1/(RELE*RELE))/(RELE*RELE))/RELE
         DNDG=DN1DG/RELE
         DNDGG=DN1DGG/RELE
         DNDXG=DN1DXG/RELE-X*DN1DG/(RELE*RELE*RELE)
         FSCR=-UP/DN*GAME
         FX=(UP*DNDX/DN-UPDX)/DN
         FXDG=((UPDG*DNDX+UPDX*DNDG+UP*DNDXG-2d0*UP*DNDX*DNDG/DN)/DN-UPDXG)/DN
         FDX=FX*GAME
         FG=(UP*DNDG/DN-UPDG)/DN
         FDG=FG*GAME-UP/DN
         FDGDH=SQG*DNDG/(DN*DN) ! d FDG / d CDH
         FDXX=((UP*DNDXX+2d0*(UPDX*DNDX-UP*DNDX*DNDX/DN))/DN-UPDXX)/DN*GAME
         FDGG=2d0*FG+GAME*((2d0*DNDG*(UPDG-UP*DNDG/DN)+UP*DNDGG)/DN-UPDGG)/DN
         FDXG=FX+GAME*FXDG
         USCR=GAME*FDG
         CVSCR=-GAME*GAME*FDGG
         PSCR=(X*FDX+GAME*FDG)/3.d0
         PDTSCR=-GAME*GAME*(X*FXDG+FDGG)/3.d0
         PDRSCR=(12d0*PSCR+X*X*FDXX+2d0*X*GAME*FDXG+GAME*GAME*FDGG)/9d0

      endif

      if (debug_FSCRliq8_table .and. .not. skip) then

         if (.not. check1(FSCR, xFSCR, 'FSCR')) return
         if (.not. check1(USCR, xUSCR, 'USCR')) return
         if (.not. check1(PSCR, xPSCR, 'PSCR')) return
         if (.not. check1(CVSCR, xCVSCR, 'CVSCR')) return
         if (.not. check1(PDTSCR, xPDTSCR, 'PDTSCR')) return
         if (.not. check1(PDRSCR, xPDRSCR, 'PDRSCR')) return

      endif

      if (use_FSCRliq8_table .and. .not. skip) then

         ! values
         FSCR% val = xFSCR
         USCR% val = xUSCR
         PSCR% val = xPSCR
         CVSCR% val = xCVSCR
         PDTSCR% val = xPDTSCR
         PDRSCR% val = xPDRSCR

         ! dT (val1) derivatives
         ! dlnRs_dT = RS% d1val1 / RS% val = 0
         dlnGAME_dT = GAME% d1val1 / GAME% val

         FSCR% d1val1 = dFSCR_dlnGAME * dlnGAME_dT
         USCR% d1val1 = dUSCR_dlnGAME * dlnGAME_dT
         PSCR% d1val1 = dPSCR_dlnGAME * dlnGAME_dT
         CVSCR% d1val1 = dCVSCR_dlnGAME * dlnGAME_dT
         PDTSCR% d1val1 = dPDTSCR_dlnGAME * dlnGAME_dT
         PDRSCR% d1val1 = dPDRSCR_dlnGAME * dlnGAME_dT

         ! dRho (val2) derivatives
         dlnRs_dRho = RS% d1val2 / RS% val
         dlnGAME_dRho = GAME% d1val2 / GAME% val

         FSCR% d1val2 = dFSCR_dlnRS * dlnRS_dRho + dFSCR_dlnGAME * dlnGAME_dRho
         USCR% d1val2 = dUSCR_dlnRS * dlnRS_dRho + dUSCR_dlnGAME * dlnGAME_dRho
         PSCR% d1val2 = dPSCR_dlnRS * dlnRS_dRho + dPSCR_dlnGAME * dlnGAME_dRho
         CVSCR% d1val2 = dCVSCR_dlnRS * dlnRS_dRho + dCVSCR_dlnGAME * dlnGAME_dRho
         PDTSCR% d1val2 = dPDTSCR_dlnRS * dlnRS_dRho + dPDTSCR_dlnGAME * dlnGAME_dRho
         PDRSCR% d1val2 = dPDRSCR_dlnRS * dlnRS_dRho + dPDRSCR_dlnGAME * dlnGAME_dRho
         
      end if

      contains

      logical function check1(v, xv, str)
         type(auto_diff_real_2var_order1), intent(in) :: v
         real(dp), intent(in) :: xv
         character (len=*), intent(in) :: str
         real(dp) :: val
         real(dp), parameter :: atol = 1d-8, rtol = 1d-6
         include 'formats'
         val = v%val
         check1 = .false.
         if (is_bad(xv)) then
            write(*,*) 'is_bad ' // trim(str), xv, val, int(Zion), RS, GAME
            return
         end if
         if (abs(val - xv) > atol + rtol*max(abs(val),abs(xv))) then
            write(*,*) 'rel mismatch ' // trim(str), &
               (val - xv)/max(abs(val),abs(xv),1d-99), &
               xv, val, int(Zion), RS%val, GAME%val
            stop 'FSCRliq8'
            return
         end if
         check1 = .true.
      end function check1

      end subroutine FSCRliq8

! ==============   SUBROUTINES FOR THE SOLID STATE   ================= *
      subroutine FSCRsol8(RS,GAMI,Zion,TPT, &
           FSCR,USCR,PSCR,S_SCR,CVSCR,PDTSCR,PDRSCR,ierr)
!                                                       Version 28.05.08
!                    undefined zero variable Q1DXG is wiped out 21.06.10
!                                               cosmetic change 16.05.13
! Fit to the el.-ion screening in bcc or fcc Coulomb solid
! Stems from FSCRsol8 v.09.06.07. Included a check for RS=0.
!   INPUT: RS - el. density parameter, GAMI - ion coupling parameter,
!          Zion - ion charge, TPT=T_p/T - ion quantum parameter
!   OUTPUT: FSCR - screening (e-i) free energy per kT per 1 ion,
!           USCR - internal energy per kT per 1 ion (screen.contrib.)
!           PSCR - pressure divided by (n_i kT) (screen.contrib.)
!           S_SCR - screening entropy contribution / (N_i k)
!           CVSCR - heat capacity per 1 ion (screen.contrib.)
!           PDTSCR,PDRSCR = PSCR + d PSCR / d ln(T,\rho)
      type(auto_diff_real_2var_order1), intent(in) :: RS,GAMI
      real(dp), intent(in) :: Zion
      type(auto_diff_real_2var_order1) :: TPT
      type(auto_diff_real_2var_order1), intent(out) :: FSCR,USCR,PSCR,S_SCR,CVSCR,PDTSCR,PDRSCR
      integer, intent(out) :: ierr

      type(auto_diff_real_2var_order1) :: XSR, P1, P2, Finf, FinfX, FinfDX, FinfDXX
      real(dp) :: Z13, ZLN, R1, R2, R3
      type(auto_diff_real_2var_order1) :: Q1U, Q1D, Q1, Q1X, Q1XDX, Q1DX, Q1DXX
      type(auto_diff_real_2var_order1) :: Y0, Y0DX, Y0DG, Y0DXX, Y0DGG, Y0DXG
      type(auto_diff_real_2var_order1) :: Y1, Y1DX, Y1DG, Y1DXX, Y1DGG, Y1DXG
      type(auto_diff_real_2var_order1) :: SA, SUPA, SUPADX, SUPADG, SUPADXX, SUPADGG, SUPADXG
      type(auto_diff_real_2var_order1) :: EM2, EM2Y1
      type(auto_diff_real_2var_order1) :: SB, SUPB, SUPBDX, SUPBDG, SUPBDXX, SUPBDGG, SUPBDXG
      type(auto_diff_real_2var_order1) :: SUP, SUPX, SUPDX, SUPG, SUPDG, SUPDXX, SUPDGG, SUPDXG
      type(auto_diff_real_2var_order1) :: GR3, GR3X, GR3DX, GR3G, GR3DG, GR3DXX, GR3DGG, GR3DXG
      type(auto_diff_real_2var_order1) :: W, WDX, WDG, WDXX, WDGG, WDXG
      type(auto_diff_real_2var_order1) :: FDX, FDG, FDXX, FDGG, FDXG
      
      real(dp) :: AP(4)
      real(dp) :: ENAT
      real(dp) :: TINY
      real(dp) :: PX


      AP(1) = 1.1866d0
      AP(2) = 0.684d0
      AP(3) = 17.9d0
      AP(4) = 41.5d0

      ENAT=2.7182818285d0
      TINY=1.d-19
      PX=0.205d0

!      data AP/1.1857,.663,17.1,40./,PX/.212/ ! for fcc lattice
      ierr = 0
      if (RS.lt.0d0) then
         ierr = -1
         return
         !stop 'FSCRsol8: RS < 0'
      end if
      if (RS.lt.TINY) then
         FSCR=0.d0
         USCR=0.d0
         PSCR=0.d0
         S_SCR=0.d0
         CVSCR=0.d0
         PDTSCR=0.d0
         PDRSCR=0.d0
         return
      endif
      XSR=0.0140047d0/RS ! relativity parameter
      Z13=pow(Zion,1d0/3d0)
      P1=0.00352d0*(1d0-AP(1)/pow(Zion,0.267d0)+0.27d0/Zion)
      P2=1d0+2.25d0/Z13* &
           (1d0+AP(2)*pow5(Zion)+0.222d0*pow6(Zion))/(1d0+0.222d0*pow6(Zion))
      ZLN=log(Zion)
      Finf=sqrt(P2/(XSR*XSR)+1d0)*Z13*Z13*P1 ! The TF limit
      FinfX=-P2/((P2+XSR*XSR)*XSR)
      FinfDX=Finf*FinfX
      FinfDXX=FinfDX*FinfX-FinfDX*(P2+3d0*XSR*XSR)/((P2+XSR*XSR)*XSR)
      R1=AP(4)/(1d0+ZLN)
      R2=0.395d0*ZLN+0.347d0/Zion/sqrt(Zion)
      R3=1d0/(1d0+ZLN*sqrt(ZLN)*0.01d0+0.097d0/(Zion*Zion))
      Q1U=R1+AP(3)*XSR*XSR
      Q1D=1d0+R2*XSR*XSR
      Q1=Q1U/Q1D
      Q1X=2d0*XSR*(AP(3)/Q1U-R2/Q1D)
      Q1XDX=Q1X/XSR+4d0*XSR*XSR*((R2/Q1D)*(R2/Q1D)-(AP(3)/Q1U)*(AP(3)/Q1U))
      Q1DX=Q1*Q1X
      Q1DXX=Q1DX*Q1X+Q1*Q1XDX
! New quantum factor, in order to suppress CVSCR at TPT >> 1
      if (TPT.lt.6d0/PX) then
         Y0=(PX*TPT)*(PX*TPT)
         Y0DX=Y0/XSR
         Y0DG=2d0*Y0/GAMI
         Y0DXX=0d0
         Y0DGG=Y0DG/GAMI
         Y0DXG=Y0DG/XSR
         Y1=exp(Y0)
         Y1DX=Y1*Y0DX
         Y1DG=Y1*Y0DG
         Y1DXX=Y1*(Y0DX*Y0DX+Y0DXX)
         Y1DGG=Y1*(Y0DG*Y0DG+Y0DGG)
         Y1DXG=Y1*(Y0DX*Y0DG+Y0DXG)
         SA=1.d0+Y1
         SUPA=log(SA)
         SUPADX=Y1DX/SA
         SUPADG=Y1DG/SA
         SUPADXX=(Y1DXX-Y1DX*Y1DX/SA)/SA
         SUPADGG=(Y1DGG-Y1DG*Y1DG/SA)/SA
         SUPADXG=(Y1DXG-Y1DX*Y1DG/SA)/SA
         EM2=ENAT-2.d0
         SB=ENAT-EM2/Y1
         SUPB=log(SB)
         EM2Y1=EM2/(Y1*Y1*SB)
         SUPBDX=EM2Y1*Y1DX
         SUPBDG=EM2Y1*Y1DG
         SUPBDXX=EM2Y1*(Y1DXX-2d0*Y1DX*Y1DX/Y1-Y1DX*SUPBDX)
         SUPBDGG=EM2Y1*(Y1DGG-2d0*Y1DG*Y1DG/Y1-Y1DG*SUPBDG)
         SUPBDXG=EM2Y1*(Y1DXG-2d0*Y1DX*Y1DG/Y1-Y1DG*SUPBDX)
         SUP=sqrt(SUPA/SUPB)
         SUPX=0.5d0*(SUPADX/SUPA-SUPBDX/SUPB)
         SUPDX=SUP*SUPX
         SUPG=0.5d0*(SUPADG/SUPA-SUPBDG/SUPB)
         SUPDG=SUP*SUPG
         SUPDXX=SUPDX*SUPX &
            + SUP*0.5d0*(SUPADXX/SUPA-(SUPADX/SUPA)*(SUPADX/SUPA) &
                  - SUPBDXX/SUPB+(SUPBDX/SUPB)*(SUPBDX/SUPB))
         SUPDGG=SUPDG*SUPG &
            + SUP*0.5d0*(SUPADGG/SUPA-(SUPADG/SUPA)*(SUPADG/SUPA) &
                    - SUPBDGG/SUPB+(SUPBDG/SUPB)*(SUPBDG/SUPB))
         SUPDXG=SUPDX*SUPG &
            + SUP*0.5d0*((SUPADXG-SUPADX*SUPADG/SUPA)/SUPA &
                    - (SUPBDXG-SUPBDX*SUPBDG/SUPB)/SUPB)
      else
         SUP=PX*TPT
         SUPDX=0.5d0*PX*TPT/XSR
         SUPDG=PX*TPT/GAMI
         SUPDXX=-0.5d0*SUPDX/XSR
         SUPDGG=0d0
         SUPDXG=SUPDX/GAMI
      endif
      GR3=pow(GAMI/SUP,R3)
      GR3X=-R3*SUPDX/SUP
      GR3DX=GR3*GR3X
      GR3DXX=GR3DX*GR3X-R3*GR3*(SUPDXX/SUP-(SUPDX/SUP)*(SUPDX/SUP))
      GR3G=R3*(1d0/GAMI-SUPDG/SUP)
      GR3DG=GR3*GR3G
      GR3DGG=GR3DG*GR3G+GR3*R3*((SUPDG/SUP)*(SUPDG/SUP)-SUPDGG/SUP-1d0/(GAMI*GAMI))
      GR3DXG=GR3DG*GR3X+GR3*R3*(SUPDX*SUPDG/(SUP*SUP)-SUPDXG/SUP)
      W=1d0+Q1/GR3
      WDX=Q1DX/GR3-Q1*GR3DX/(GR3*GR3*GR3)
      WDG=-Q1*GR3DG/(GR3*GR3)
      WDXX=Q1DXX/GR3-(2d0*Q1DX*GR3DX+Q1*(GR3DXX-2d0*GR3DX*GR3DX/GR3))/(GR3*GR3)
      WDGG=Q1*(2d0*GR3DG*GR3DG/GR3-GR3DGG)/(GR3*GR3)
      WDXG=-(Q1DX*GR3DG+Q1*(GR3DXG-2d0*GR3DX*GR3DG/GR3))/(GR3*GR3)
      FSCR=-GAMI*Finf*W
      FDX=-GAMI*(FinfDX*W+Finf*WDX)
      FDXX=-GAMI*(FinfDXX*W+2d0*FinfDX*WDX+Finf*WDXX)
      FDG=-Finf*W-GAMI*Finf*WDG
      FDGG=-2d0*Finf*WDG-GAMI*Finf*WDGG
      FDXG=-FinfDX*W-Finf*WDX-GAMI*(FinfDX*WDG+Finf*WDXG)
      S_SCR=-GAMI*GAMI*Finf*WDG
      USCR=S_SCR+FSCR
      CVSCR=-GAMI*GAMI*FDGG
      PSCR=(XSR*FDX+GAMI*FDG)/3d0
      PDTSCR=GAMI*GAMI*(XSR*Finf*(FinfX*WDG+WDXG)-FDGG)/3d0
      PDRSCR=(12d0*PSCR+XSR*XSR*FDXX+2d0*XSR*GAMI*FDXG &
         + GAMI*GAMI*FDGG)/9d0
      return
      end subroutine FSCRsol8

      subroutine ANHARM8(GAMI,TPT,Fah,Uah,Pah,CVah,PDTah,PDRah)
! ANHARMONIC free energy                                Version 27.07.07
!                                                       cleaned 16.06.09
! Stems from ANHARM8b. Difference: AC=0., B1=.12 (.1217 - over accuracy)
! Input: GAMI - ionic Gamma, TPT=Tp/T - ionic quantum parameter
! Output: anharm.free en. Fah=F_{AH}/(N_i kT), internal energy Uah,
!   pressure Pah=P_{AH}/(n_i kT), specific heat CVah = C_{V,AH}/(N_i k),
!   PDTah = Pah + d Pah / d ln T, PDRah = Pah + d Pah / d ln\rho
      type(auto_diff_real_2var_order1), intent(in) :: GAMI,TPT
      type(auto_diff_real_2var_order1), intent(out) :: Fah,Uah,Pah,CVah,PDTah,PDRah
      integer, parameter :: NM=3
      
      integer :: N
      real(dp) :: AA(3)
      type(auto_diff_real_2var_order1) :: CK, TPT2, TPT4, TQ, TK2, SUP
      type(auto_diff_real_2var_order1) :: CN, SUPGN, ACN, PN
      type(auto_diff_real_2var_order1) :: B1

      B1 = 0.12d0 ! coeff.at \eta^2/\Gamma at T=0

      ! Farouki & Hamaguchi'93
      AA(1) = 10.9d0
      AA(2) = 247.0d0
      AA(3) = 1.765d5

      
      CK=B1/AA(1) ! fit coefficient
      TPT2=TPT*TPT
      TPT4=TPT2*TPT2
      TQ=B1*TPT2/GAMI ! quantum dependence
      TK2=CK*TPT2
      SUP=exp(-TK2) ! suppress.factor of class.anharmonicity
      Fah=0.d0
      Uah=0.d0
      Pah=0.d0
      CVah=0.d0
      PDTah=0.d0
      PDRah=0.d0
      SUPGN=SUP
      do N=1,NM
         CN=1d0*N
         SUPGN=SUPGN/GAMI ! SUP/Gamma^n
         ACN=AA(N)
         Fah=Fah-ACN/CN*SUPGN
         Uah=Uah+(ACN*(1d0+2d0*TK2/CN))*SUPGN
         PN=AA(N)/3d0+TK2*AA(N)/CN
         Pah=Pah+PN*SUPGN
         CVah=CVah+((CN+1d0)*AA(N)+(4d0-2d0/CN)*AA(N)*TK2 &
            + 4d0*AA(N)*CK*CK/CN*TPT4)*SUPGN
         PDTah=PDTah+(PN*(1d0+CN+2d0*TK2)-2d0/CN*AA(N)*TK2)*SUPGN
         PDRah=PDRah+(PN*(1d0-CN/3d0-TK2)+AA(N)/CN*TK2)*SUPGN
      enddo
      Fah=Fah-TQ
      Uah=Uah-TQ
      Pah=Pah-TQ/1.5d0
      PDRah=PDRah-TQ/4.5d0
      return
      end subroutine ANHARM8

      subroutine FHARM12(GAMI,TPT, &
         Fharm,Uharm,Pharm,CVth,Sharm,PDTharm,PDRharm)
! Thermodynamic functions of a harmonic crystal, incl.stat.Coul.lattice
! 
!                                                       Version 27.04.12
! Stems from FHARM8 v.15.02.08
! Replaced HLfit8 with HLfit12: rearranged output.
! Input: GAMI - ionic Gamma, TPT=T_{p,i}/T
! Output: Fharm=F/(N_i T), Uharm=U/(N_i T), Pharm=P/(n_i T),
! CVth=C_V/N_i, Sharm=S/N_i
! PDTharm = Pharm + d Pharm / d ln T, PDRharm = Pharm + d Pharm/d ln\rho
      type(auto_diff_real_2var_order1), intent(in) :: GAMI,TPT
      type(auto_diff_real_2var_order1), intent(out) :: Fharm,Uharm,Pharm,CVth,Sharm,PDTharm,PDRharm

      type(auto_diff_real_2var_order1) :: Fth,Uth,Sth,U0,E0
      type(auto_diff_real_2var_order1) :: F,U,U1
      
      real(dp), parameter :: CM = .895929256d0 ! Madelung
      
      call HLfit12(TPT,F,U,CVth,Sth,U1,1)
      U0=-CM*GAMI ! perfect lattice
      E0=1.5d0*U1*TPT ! zero-point energy
      Uth=U+E0
      Fth=F+E0
      Uharm=U0+Uth
      Fharm=U0+Fth
      Sharm=Sth
      Pharm=U0/3d0+Uth/2d0
      PDTharm=0.5d0*CVth
      PDRharm=U0/2.25d0+0.75d0*Uth-0.25d0*CVth
      return
      end subroutine FHARM12

      subroutine HLfit12(eta,F,U,CV,S,U1,LATTICE)
!                                                       Version 24.04.12
! Stems from HLfit8 v.03.12.08;
!   differences: E0 excluded from  U and F;
!   U1 and d(CV)/d\ln(T) are added on the output.
! Fit to thermal part of the thermodynamic functions.
! Baiko, Potekhin, & Yakovlev (2001).
! Zero-point lattice quantum energy 1.5u_1\eta EXCLUDED (unlike HLfit8).
! Input: eta=Tp/T, LATTICE=1 for bcc, 2 for fcc
! Output: F and U (normalized to NkT) - due to phonon excitations,
!   CV and S (normalized to Nk) in the HL model,
!   U1 - the 1st phonon moment,

      type(auto_diff_real_2var_order1) :: eta ! can be modified, not sure if this is an intended side-effect
      type(auto_diff_real_2var_order1), intent(out) :: F, U, CV, S, U1
      integer, intent(in) :: LATTICE

      real(dp) :: CLM, ALPHA, BETA, GAMMA
      real(dp) :: A1, A2, A3, A4, A6, A8
      real(dp) :: B0, B2, B4, B5, B6, B7, B9, B11
      real(dp) :: C9, C11
      type(auto_diff_real_2var_order1) :: UP, DN, EA, EB, EG, UP1, UP2, DN1, DN2, E0
      
      real(dp) :: EPS
      real(dp) :: TINY
      
      EPS=1.d-5
      TINY=1.d-99

      if (LATTICE.eq.1) then ! bcc lattice
         CLM=-2.49389d0 ! 3*ln<\omega/\omega_p>
         U1=0.5113875d0
         ALPHA=0.265764d0
         BETA=0.334547d0
         GAMMA=0.932446d0
         A1=0.1839d0
         A2=0.593586d0
         A3=0.0054814d0
         A4=5.01813d-4
         A6=3.9247d-7
         A8=5.8356d-11
         B0=261.66d0
         B2=7.07997d0
         B4=0.0409484d0
         B5=0.000397355d0
         B6=5.11148d-5
         B7=2.19749d-6
         C9=0.004757014d0
         C11=0.0047770935d0
      elseif (LATTICE.eq.2) then ! fcc lattice
         CLM=-2.45373d0
         U1=0.513194d0
         ALPHA=0.257591d0
         BETA=0.365284d0
         GAMMA=0.9167070d0
         A1=0.0d0
         A2=0.532535d0
         A3=0.0d0
         A4=3.76545d-4
         A6=2.63013d-7
         A8=6.6318d-11
         B0=303.20d0
         B2=7.7255d0
         B4=0.0439597d0
         B5=0.000114295d0
         B6=5.63434d-5
         B7=1.36488d-6
         C9=0.00492387d0
         C11=0.00437506d0
      else
         stop 'HLfit: unknown lattice type'
      endif
      if (eta.gt.1d0/EPS) then ! asymptote of Eq.(13) of BPY'01
         U=3d0/(C11*eta*eta*eta)
         F=-U/3d0
         CV=4d0*U
      else if (eta.lt.EPS) then ! Eq.(17) of BPY'01
         if (eta.lt.TINY) eta = TINY !stop 'HLfit8: eta is too small'
         F=3d0*log(eta)+CLM-1.5d0*U1*eta+eta*eta/24.d0
         U=3d0-1.5d0*U1*eta+eta*eta/12d0
         CV=3d0-eta*eta/12d0
      else
         B9=A6*C9
         B11=A8*C11

         !UP=1d0+A1*eta+A2*eta**2+A3*eta**3+A4*eta**4+A6*eta**6+A8*eta**8
         UP=1d0+eta*(A1+eta*(A2+eta*(A3+eta*(A4+eta*eta*(A6+eta*eta*A8)))))

         !DN=B0+B2*eta**2+B4*eta**4+B5*eta**5+B6*eta**6+B7*eta**7+B9*eta**9+B11*eta**11
         DN=B0+eta*eta*(B2+eta*eta*(B4+eta*(B5+eta*(B6+eta*(B7+eta*eta*(B9+eta*eta*B11))))))

         EA=exp(-ALPHA*eta)
         EB=exp(-BETA*eta)
         EG=exp(-GAMMA*eta)
         F=log(1.d0-EA)+log(1.d0-EB)+log(1.d0-EG)-UP/DN ! F_{thermal}/NT

         !UP1=A1+2d0*A2*eta+3.*A3*eta**2+4.*A4*eta**3+6d0*A6*eta**5+8.*A8*eta**7
         UP1=A1+eta*(2d0*A2+eta*(3.d0*A3+eta*(4.d0*A4+eta*eta*(6d0*A6+eta*eta*8d0*A8))))

         !UP2=2d0*A2+6d0A3*eta+12d0*A4*eta**2+30.*A6*eta**4+56d0A8*eta**6
         UP2=2d0*A2+eta*(6d0*A3+eta*(12d0*A4+eta*eta*(30d0*A6+eta*eta*56d0*A8)))

         !DN1=2d0*B2*eta+4.*B4*eta**3+5.*B5*eta**4+6d0*B6*eta**5+7.*B7*eta**6+9.*B9*eta**8+11d0*B11*eta**10.
         DN1=eta*(2d0*B2+eta*eta*(4d0*B4+eta*(5d0*B5+eta*(6d0*B6+eta*(7d0*B7+eta*eta*(9d0*B9+eta*eta*11d0*B11))))))

         !DN2=2d0*B2+12d0*B4*eta**2+20.*B5*eta**3+30.*B6*eta**4+42d0*B7*eta**5+72d0*B9*eta**7+110.*B11*eta**9
         DN2=2d0*B2+eta*eta*(12d0*B4+eta*(20d0*B5+eta*(30d0*B6+eta*(42d0*B7+eta*eta*(72d0*B9+eta*eta*110d0*B11)))))

         U=ALPHA*EA/(1.d0-EA)+BETA*EB/(1.d0-EB)+GAMMA*EG/(1.d0-EG) &
            - (UP1*DN-DN1*UP)/(DN*DN) ! int.en./NT/eta
         CV=ALPHA*ALPHA*EA/((1.d0-EA)*(1.d0-EA))+BETA*BETA*EB/((1.d0-EB)*(1.d0-EB)) &
            + GAMMA*GAMMA*EG/((1.d0-EG)*(1.d0-EG)) &
            + ((UP2*DN-DN2*UP)*DN-2d0*(UP1*DN-DN1*UP)*DN1)/(DN*DN*DN) ! cV/eta^2
         U=U*eta
         CV=CV*eta*eta
      end if
      S=U-F
      return
      end subroutine HLfit12

      subroutine CORMIX(RS,GAME,Zmean,Z2mean,Z52,Z53,Z321, &
        FMIX,UMIX,PMIX,CVMIX,PDTMIX,PDRMIX)
!                                                       Version 02.07.09
! Correction to the linear mixing rule for moderate to small Gamma
! Input: RS=r_s (if RS=0, then OCP, otherwise EIP)
!        GAME=\Gamma_e
!        Zmean=<Z> (average Z of all ions, without electrons)
!        Z2mean=<Z^2>, Z52=<Z^2.5>, Z53=<Z^{5/3}>, Z321=<Z(Z+1)^1.5>
! Output: FMIX=\Delta f - corr.to the reduced free energy f=F/N_{ion}kT
!         UMIX=\Delta u - corr.to the reduced internal energy u
!         PMIX=\Delta u - corr.to the reduced pressure P=P/n_{ion}kT
!         CVMIX=\Delta c - corr.to the reduced heat capacity c_V
!         PDTMIX=(1/n_{ion}kT)d\Delta P / d ln T
!               = \Delta p +  d \Delta p / d ln T
!         PDRMIX=(1/n_{ion}kT)d\Delta P / d ln n_e
! (composition is assumed fixed: Zmean,Z2mean,Z52,Z53=constant)
      type(auto_diff_real_2var_order1), intent(in) :: RS,GAME
      real(dp), intent(in) :: Zmean,Z2mean,Z52,Z53,Z321
      type(auto_diff_real_2var_order1), intent(out) :: FMIX,UMIX,PMIX,CVMIX,PDTMIX,PDRMIX

      type(auto_diff_real_2var_order1) :: GAMImean, Dif0, DifR, DifFDH, D
      type(auto_diff_real_2var_order1) :: P3, D0, GP, FMIX0, Q, R, GQ, G, GDG, UDG
      
      real(dp) :: TINY

      TINY=1.d-9
      
      GAMImean=GAME*Z53
      if (RS.lt.TINY) then ! OCP
         Dif0=Z52-sqrt(Z2mean*Z2mean*Z2mean/Zmean)
      else
         Dif0=Z321-sqrt((Z2mean+Zmean)*(Z2mean+Zmean)*(Z2mean+Zmean)/Zmean)
      endif
      DifR=Dif0/Z52
      DifFDH=Dif0*GAME*sqrt(GAME/3d0) ! F_DH - F_LM(DH)
      D=Z2mean/(Zmean*Zmean)
      if (abs(D-1.d0).lt.TINY) then ! no correction
         FMIX=0d0
         UMIX=0d0
         PMIX=0d0
         CVMIX=0d0
         PDTMIX=0d0
         PDRMIX=0d0
         return
      endif
      P3=pow(D,-0.2d0)
      D0=(2.6d0*DifR+14d0*DifR*DifR*DifR)/(1.d0-P3)
      GP=D0*pow(GAMImean,P3)
      FMIX0=DifFDH/(1d0+GP)
      Q=D*D*0.0117d0
      R=1.5d0/P3-1d0
      GQ=Q*GP
      FMIX=FMIX0/pow(1d0+GQ,R)
      G=1.5d0-P3*GP/(1d0+GP)-R*P3*GQ/(1d0+GQ)
      UMIX=FMIX*G
      PMIX=UMIX/3.d0
      GDG=-P3*P3*(GP/((1.d0+GP)*(1.d0+GP))+R*GQ/((1.d0+GQ)*(1.d0+GQ))) ! d G /d ln Gamma
      UDG=UMIX*G+FMIX*GDG ! d u_mix /d ln Gamma
      CVMIX=UMIX-UDG
      PDTMIX=PMIX-UDG/3d0
      PDRMIX=PMIX+UDG/9d0
      return
      end subroutine CORMIX

! ===================  IDEAL ELECTRON GAS  =========================== *
      subroutine ELECT11(TEMP,CHI, &
        DENS,FEid,PEid,UEid,SEid,CVE,CHITE,CHIRE, &
        DlnDH,DlnDT,DlnDHH,DlnDTT,DlnDHT)
!                                                       Version 17.11.11
! ELECT9 v.04.03.09 + smooth match of two fits at chi >> 1 + add.outputs
! Compared to ELECTRON v.06.07.00, this S/R is completely rewritten: 
!        numerical differentiation is avoided now.
! Compared to ELECT7 v.06.06.07,
!    - call BLIN7 is changed to call BLIN9,
!    - Sommerfeld expansion is used at chi >~ 28 i.o. 1.e4
!    - Sommerfeld expansion is corrected: introduced DeltaEF, D1 and D2.
! Ideal electron-gas EOS.
! Input: TEMP - T [a.u.], CHI=\mu/T
! Output: DENS - electron number density n_e [a.u.],
!         FEid - free energy / N_e kT, UEid - internal energy / N_e kT,
!         PEid - pressure (P_e) / n_e kT, SEid - entropy / N_e k,
!         CVE - heat capacity / N_e k,
!         CHITE=(d ln P_e/d ln T)_V, CHIRE=(d ln P_e/d ln n_e)_T
!         DlnDH=(d ln n_e/d CHI)_T = (T/n_e) (d n_e/d\mu)_T
!         DlnDT=(d ln n_e/d ln T)_CHI
!         DlnDHH=(d^2 ln n_e/d CHI^2)_T
!         DlnDTT=(d^2 ln n_e/d (ln T)^2)_CHI
!         DlnDHT=d^2 ln n_e/d (ln T) d CHI
      type(auto_diff_real_2var_order1), intent(in) :: TEMP,CHI
      type(auto_diff_real_2var_order1), intent(out) :: DENS,FEid,PEid,UEid,SEid,CVE,CHITE,CHIRE,DlnDH,DlnDT,DlnDHH,DlnDTT,DlnDHT

      type(auto_diff_real_2var_order1) :: X2, FP, FM
      type(auto_diff_real_2var_order1) :: DENSa,FEida,PEida,UEida,SEida,CVEa,CHITEa,CHIREa,DlnDHa,DlnDTa,DlnDHHa,DlnDTTa,DlnDHTa
      type(auto_diff_real_2var_order1) :: DENSb,FEidb,PEidb,UEidb,SEidb,CVEb,CHITEb,CHIREb,DlnDHb,DlnDTb,DlnDHHb,DlnDTTb,DlnDHTb
      
      type(auto_diff_real_2var_order1) :: CHI1
      type(auto_diff_real_2var_order1) :: CHI2
      type(auto_diff_real_2var_order1) :: XMAX
      type(auto_diff_real_2var_order1) :: DCHI1
      type(auto_diff_real_2var_order1) :: DCHI2
      type(auto_diff_real_2var_order1) :: XSCAL2

      CHI1=0.6d0
      CHI2=28.d0
      XMAX=20.d0
      DCHI1=0.1d0
      DCHI2=CHI2-CHI1-DCHI1
      XSCAL2=XMAX/DCHI2
      
      X2=(CHI-CHI2)*XSCAL2
      if (X2.lt.-XMAX) then
         call ELECT11a(TEMP,CHI, &
           DENS,FEid,PEid,UEid,SEid,CVE,CHITE,CHIRE,&
           DlnDH,DlnDT,DlnDHH,DlnDTT,DlnDHT)
      elseif (X2.gt.XMAX) then
         call ELECT11b(TEMP,CHI, &
           DENS,FEid,PEid,UEid,SEid,CVE,CHITE,CHIRE,&
           DlnDH,DlnDT,DlnDHH,DlnDTT,DlnDHT)
      else
         call FERMI10(X2,XMAX,FP,FM)
         call ELECT11a(TEMP,CHI, &
           DENSa,FEida,PEida,UEida,SEida,CVEa,CHITEa,CHIREa, &
           DlnDHa,DlnDTa,DlnDHHa,DlnDTTa,DlnDHTa)
         call ELECT11b(TEMP,CHI, &
           DENSb,FEidb,PEidb,UEidb,SEidb,CVEb,CHITEb,CHIREb,&
           DlnDHb,DlnDTb,DlnDHHb,DlnDTTb,DlnDHTb)
         DENS=DENSa*FP+DENSb*FM
         FEid=FEida*FP+FEidb*FM
         PEid=PEida*FP+PEidb*FM
         UEid=UEida*FP+UEidb*FM
         SEid=SEida*FP+SEidb*FM
         CVE=CVEa*FP+CVEb*FM
         CHITE=CHITEa*FP+CHITEb*FM
         CHIRE=CHIREa*FP+CHIREb*FM
         DlnDH=DlnDHa*FP+DlnDHb*FM
         DlnDT=DlnDTa*FP+DlnDTb*FM
         DlnDHH=DlnDHHa*FP+DlnDHHb*FM
         DlnDHT=DlnDHTa*FP+DlnDHTb*FM
         DlnDTT=DlnDTTa*FP+DlnDTTb*FM
      endif
      return
      end subroutine ELECT11

      subroutine ELECT11a(TEMP,CHI, &
        DENS,FEid,PEid,UEid,SEid,CVE,CHITE,CHIRE, &
        DlnDH,DlnDT,DlnDHH,DlnDTT,DlnDHT)
!                                                       Version 16.11.11
! This is THE FIRST PART of ELECT9 v.04.03.09.
      type(auto_diff_real_2var_order1), intent(in) :: TEMP,CHI
      type(auto_diff_real_2var_order1), intent(out) :: DENS,FEid,PEid,UEid,SEid,CVE,CHITE,CHIRE,DlnDH,DlnDT,DlnDHH,DlnDTT,DlnDHT

      type(auto_diff_real_2var_order1) :: TEMR
      type(auto_diff_real_2var_order1) :: W0,W0DX,W0DT,W0DXX,W0DTT,W0DXT
      type(auto_diff_real_2var_order1) :: W1,W1DX,W1DT,W1DXX,W1DTT,W1DXT
      type(auto_diff_real_2var_order1) :: W2,W2DX,W2DT,W2DXX,W2DTT,W2DXT
      type(auto_diff_real_2var_order1) :: W0XXX,W0XTT,W0XXT
      type(auto_diff_real_2var_order1) :: TPI, DENR, PR, U
      type(auto_diff_real_2var_order1) :: dndT, dPdT, dUdT, dndH, dPdH, dUdH      
      type(auto_diff_real_2var_order1) :: dndHH, dndHT, dndTT

      type(auto_diff_real_2var_order1) :: BOHR
      type(auto_diff_real_2var_order1) :: PI2
      type(auto_diff_real_2var_order1) :: BOHR2
      type(auto_diff_real_2var_order1) :: BOHR3
      
      BOHR=137.036d0
      PI2=PI*PI
      BOHR2=BOHR*BOHR
      BOHR3=BOHR2*BOHR !cleaned 15/6

      TEMR=TEMP/BOHR2 ! T in rel.units (=T/mc^2)
      call BLIN9(TEMR,CHI, &
        W0,W0DX,W0DT,W0DXX,W0DTT,W0DXT, &
        W1,W1DX,W1DT,W1DXX,W1DTT,W1DXT, &
        W2,W2DX,W2DT,W2DXX,W2DTT,W2DXT, &
        W0XXX,W0XTT,W0XXT)
      TPI=TEMR*sqrt(2.d0*TEMR)/PI2 ! common pre-factor
      DENR=TPI*(W1*TEMR+W0)
      PR=TEMR*TPI/3.d0*(W2*TEMR+2d0*W1)
      U=TEMR*TPI*(W2*TEMR+W1)
! (these are density, pressure, and internal energy in the rel.units)
      PEid=PR/(DENR*TEMR)
      UEid=U/(DENR*TEMR)
      FEid=CHI-PEid
      DENS=DENR*BOHR3 ! converts from rel.units to a.u.
      SEid=UEid-FEid
! derivatives over T at constant chi:
      dndT=TPI*(1.5d0*W0/TEMR+2.5d0*W1+W0DT+TEMR*W1DT) ! (d n_e/dT)_\chi
      dPdT=TPI/3.d0*(5.d0*W1+2d0*TEMR*W1DT+3.5d0*TEMR*W2+TEMR*TEMR*W2DT)!dP/dT
      dUdT=TPI*(2.5d0*W1+TEMR*W1DT+3.5d0*TEMR*W2+TEMR*TEMR*W2DT)!dU/dT_\chi
! derivatives over chi at constant T:
      dndH=TPI*(W0DX+TEMR*W1DX) ! (d n_e/d\chi)_T
      dndHH=TPI*(W0DXX+TEMR*W1DXX) ! (d^2 n_e/d\chi)_T
      dndTT=TPI*(0.75d0*W0/pow2(TEMR)+3d0*W0DT/TEMR+W0DTT+3.75d0*W1/TEMR+5d0*W1DT+TEMR*W1DTT)
      dndHT=TPI*(1.5d0*W0DX/TEMR+W0DXT+2.5d0*W1DX+TEMR*W1DXT)
      DlnDH=dndH/DENR ! (d ln n_e/d\chi)_T
      DlnDT=dndT*TEMR/DENR ! (d ln n_e/d ln T)_\chi
      DlnDHH=dndHH/DENR-pow2(DlnDH) ! (d^2 ln n_e/d\chi^2)_T
      DlnDTT=pow2(TEMR)/DENR*dndTT+DlnDT-pow2(DlnDT) ! d^2 ln n_e/d ln T^2
      DlnDHT=TEMR/DENR*(dndHT-dndT*DlnDH) ! d^2 ln n_e/d\chi d ln T
      dPdH=TPI/3d0*TEMR*(2d0*W1DX+TEMR*W2DX) ! (d P_e/d\chi)_T
      dUdH=TPI*TEMR*(W1DX+TEMR*W2DX) ! (d U_e/d\chi)_T
      CVE=(dUdT-dUdH*dndT/dndH)/DENR
      CHITE=TEMR/PR*(dPdT-dPdH*dndT/dndH)
      CHIRE=DENR/PR*dPdH/dndH ! (dndH*TEMR*PEid) ! DENS/PRE*dPdH/dndH
      return
      end subroutine ELECT11a

      subroutine ELECT11b(TEMP,CHI, &
        DENS,FEid,PEid,UEid,SEid,CVE,CHITE,CHIRE, &
        DlnDH,DlnDT,DlnDHH,DlnDTT,DlnDHT)
!                                                       Version 17.11.11
! Stems from ELECT9b v.19.01.10, Diff. - additional output.
! Sommerfeld expansion at very large CHI.
      type(auto_diff_real_2var_order1), intent(in) :: TEMP,CHI
      type(auto_diff_real_2var_order1), intent(out) :: DENS,FEid,PEid,UEid,SEid,CVE,CHITE,CHIRE,DlnDH,DlnDT,DlnDHH,DlnDTT,DlnDHT

      type(auto_diff_real_2var_order1) :: TEMR, EF, DeltaEF, G, PF, F, DF, P, DelP, S, U
      type(auto_diff_real_2var_order1) :: DENR, DT, D1, D2
      type(auto_diff_real_2var_order1) :: TPI, dndH, dndT, dndHH, dndHT, dndTT
      type(auto_diff_real_2var_order1) :: W0,W0DX,W0DT,W0DXX,W0DTT,W0DXT,W1,W1DX,W1DT,W1DXX,W1DTT,W1DXT,W2,W2DX,W2DT,W2DXX,W2DTT,W2DXT,W0XXX,W0XTT,W0XXT      

      type(auto_diff_real_2var_order1) :: BOHR
      type(auto_diff_real_2var_order1) :: PI2
      type(auto_diff_real_2var_order1) :: BOHR2
      type(auto_diff_real_2var_order1) :: BOHR3

      BOHR=137.036d0
      PI2=PI*PI
      BOHR2=BOHR*BOHR
      BOHR3=BOHR2*BOHR !cleaned 15/6

      TEMR=TEMP/BOHR2 ! T in rel.units (=T/mc^2)
      EF=CHI*TEMR ! Fermi energy in mc^2 - zeroth aprox. = CMU1
      DeltaEF=PI2*TEMR*TEMR/6.d0*(1.d0+2.d0*EF*(2.d0+EF)) &
             /(EF*(1.d0+EF)*(2.d0+EF)) ! corr. [page 125, equiv.Eq.(6) of PC'10]]
      EF=EF+DeltaEF ! corrected Fermi energy (14.02.09)
      G=1.d0+EF ! electron Lorentz-factor
      if (EF.gt.1.d-5) then ! relativistic expansion (Yak.&Shal.'89)
        PF=sqrt(G*G-1.d0) ! Fermi momentum [rel.un.=mc]
        F=(PF*(1d0+2.d0*PF*PF)*G-PF*PF*PF/0.375d0-log(PF+G))/8.d0/PI2!F/V
        DF=-TEMR*TEMR*PF*G/6.d0 ! thermal correction to F/V
        P=(PF*G*(PF*PF/1.5d0-1.d0)+log(PF+G))/8.d0/PI2 ! P(T=0)
        DelP=TEMR*TEMR*PF*(PF*PF+2.d0)/G/18.d0 ! thermal correction to P
        CVE=PI2*TEMR*G/(PF*PF)
      else ! nonrelativistic limit
        PF=sqrt(2.d0*EF)
        F=pow5(PF)*0.1d0/PI2
        DF=-TEMR*TEMR*PF/6.d0
        P=F/1.5d0
        DelP=TEMR*TEMR*PF/9.d0
        CVE=PI2*TEMR/EF/2.d0
      endif
      F=F+DF
      P=P+DelP
      S=-2.d0*DF ! entropy per unit volume [rel.un.]
      U=F+S
      CHIRE=pow5(PF)/(9.d0*PI2*P*G)
      CHITE=2.d0*DelP/P
      DENR=PF*PF*PF/3.d0/PI2 ! n_e [rel.un.=\Compton^{-3}]
      DENS=DENR*BOHR3 ! conversion to a.u.(=\Bohr_radius^{-3})
! derivatives over chi at constant T and T at constant chi:
      TPI=TEMR*sqrt(2.d0*TEMR)/PI2 ! common pre-factor
      call SOMMERF(TEMR,CHI, &
        W0,W0DX,W0DT,W0DXX,W0DTT,W0DXT, &
        W1,W1DX,W1DT,W1DXX,W1DTT,W1DXT, &
        W2,W2DX,W2DT,W2DXX,W2DTT,W2DXT,&
        W0XXX,W0XTT,W0XXT)
      dndH=TPI*(W0DX+TEMR*W1DX) ! (d n_e/d\chi)_T
      dndT=TPI*(1.5d0*W0/TEMR+2.5d0*W1+W0DT+TEMR*W1DT) ! (d n_e/dT)_\chi
      dndHH=TPI*(W0DXX+TEMR*W1DXX) ! (d^2 n_e/d\chi)_T
      dndTT=TPI*(0.75d0*W0/pow2(TEMR)+3d0*W0DT/TEMR+W0DTT+3.75d0*W1/TEMR+5d0*W1DT+TEMR*W1DTT)
      dndHT=TPI*(1.5d0*W0DX/TEMR+W0DXT+2.5d0*W1DX+TEMR*W1DXT)
      DlnDH=dndH/DENR ! (d ln n_e/d\chi)_T
      DlnDT=dndT*TEMR/DENR ! (d ln n_e/d ln T)_\chi
      DlnDHH=dndHH/DENR-pow2(DlnDH) ! (d^2 ln n_e/d\chi^2)_T
      DlnDTT=pow2(TEMR)/DENR*dndTT+DlnDT-pow2(DlnDT) ! d^2 ln n_e/d ln T^2
      DlnDHT=TEMR/DENR*(dndHT-dndT*DlnDH) ! d^2 ln n_e/d\chi d ln T
      DT=DENR*TEMR
      PEid=P/DT
      UEid=U/DT
      FEid=F/DT
      SEid=S/DT
! Empirical corrections of 16.02.09:
      D1=DeltaEF/EF
      D2=D1*(4.d0-2.d0*(PF/G))
      CVE=CVE/(1.d0+D2)
      SEid=SEid/(1.d0+D1)
      CHITE=CHITE/(1.d0+D2)
      return
      end subroutine ELECT11b

      subroutine SOMMERF(TEMR,CHI, &
        W0,W0DX,W0DT,W0DXX,W0DTT,W0DXT, &
        W1,W1DX,W1DT,W1DXX,W1DTT,W1DXT, &
        W2,W2DX,W2DT,W2DXX,W2DTT,W2DXT, &
        W0XXX,W0XTT,W0XXT)
!                                                       Version 17.11.11
! Sommerfeld expansion for the Fermi-Dirac integrals
! Input: TEMR=T/mc^2; CHI=(\mu-mc^2)/T
! Output: Wk - Fermi-Dirac integral of the order k+1/2
!         WkDX=dWk/dCHI, WkDT = dWk/dT, WkDXX=d^2 Wk / d CHI^2,
!         WkDTT=d^2 Wk / d T^2, WkDXT=d^2 Wk /dCHIdT,
!         W0XXX=d^3 W0 / d CHI^3, W0XTT=d^3 W0 /(d CHI d^2 T),
!         W0XXT=d^3 W0 /dCHI^2 dT
! [Draft source: yellow book pages 124-127]

      type(auto_diff_real_2var_order1), intent(in) :: TEMR,CHI
      type(auto_diff_real_2var_order1), intent(out) :: W0,W0DX,W0DT,W0DXX,W0DTT,W0DXT
      type(auto_diff_real_2var_order1), intent(out) :: W1,W1DX,W1DT,W1DXX,W1DTT,W1DXT
      type(auto_diff_real_2var_order1), intent(out) :: W2,W2DX,W2DT,W2DXX,W2DTT,W2DXT
      type(auto_diff_real_2var_order1), intent(out) :: W0XXX,W0XTT,W0XXT

      type(auto_diff_real_2var_order1) :: CMU, CMU1, PIT26, CN0, CN1, CN2
      type(auto_diff_real_2var_order1) :: CJ00,CJ10,CJ20,CJ01,CJ11,CJ21,CJ02,CJ12,CJ22,CJ03,CJ13,CJ23,CJ04,CJ14,CJ24,CJ05

      if (CHI.lt..5d0) stop 'SOMMERF: non-degenerate (small CHI)'
      if (TEMR.le.0.d0) stop 'SOMMERF: T < 0'
      CMU1=CHI*TEMR ! chemical potential in rel.units
      CMU=1.d0+CMU1
      call SUBFERMJ(CMU1, &
        CJ00,CJ10,CJ20, &
        CJ01,CJ11,CJ21, &
        CJ02,CJ12,CJ22, &
        CJ03,CJ13,CJ23, &
        CJ04,CJ14,CJ24,CJ05)
      PIT26=pow2(PI*TEMR)/6.d0
!!!     PITAU4=pow2(PIT26)*0.7d0
      CN0=sqrt(.5d0/TEMR)/TEMR
      CN1=CN0/TEMR
      CN2=CN1/TEMR
      W0=CN0*(CJ00+PIT26*CJ02) ! +CN0*PITAU4*CJ04
      W1=CN1*(CJ10+PIT26*CJ12) ! +CN1*PITAU4*CJ14
      W2=CN2*(CJ20+PIT26*CJ22) ! +CN2*PITAU4*CJ24
      W0DX=CN0*TEMR*(CJ01+PIT26*CJ03) ! +CN0*PITAU4*CJ05
      W1DX=CN0*(CJ11+PIT26*CJ13)
      W2DX=CN1*(CJ21+PIT26*CJ23)
      W0DT=CN1*(CMU1*CJ01-1.5d0*CJ00+PIT26*(CMU1*CJ03+.5d0*CJ02))
!!!     +  CN1*PITAU4*(CMU1*CJ05+2.5d0*CJ04)
      W1DT=CN2*(CMU1*CJ11-2.5d0*CJ10+PIT26*(CMU1*CJ13-.5d0*CJ12))
      W2DT=CN2/TEMR*(CMU1*CJ21-3.5d0*CJ20+PIT26*(CMU1*CJ23-1.5d0*CJ22))
      W0DXX=CN0*pow2(TEMR)*(CJ02+PIT26*CJ04)
      W1DXX=CN0*TEMR*(CJ12+PIT26*CJ14)
      W2DXX=CN0*(CJ22+PIT26*CJ24)
      W0DXT=CN0*(CMU1*CJ02-.5d0*CJ01+PIT26*(CMU1*CJ04+1.5d0*CJ03))
      W1DXT=CN1*(CMU1*CJ12-1.5d0*CJ11+PIT26*(CMU1*CJ14+.5d0*CJ13))
      W2DXT=CN2*(CMU1*CJ22-2.5d0*CJ21+PIT26*(CMU1*CJ24-.5d0*CJ23))
      W0DTT=CN2*(3.75d0*CJ00-3.d0*CMU1*CJ01+pow2(CMU1)*CJ02+PIT26*(-.25d0*CJ02+CMU1*CJ03+pow2(CMU1)*CJ04))
      W1DTT=CN2/TEMR*(8.75d0*CJ10-5.d0*CMU1*CJ11+pow2(CMU1)*CJ12+PIT26*(.75d0*CJ12-CMU1*CJ13+pow2(CMU1)*CJ14))
      W2DTT=CN2/pow2(TEMR)*(15.75d0*CJ20-7.d0*CMU1*CJ21+pow2(CMU1)*CJ22+PIT26*(3.75d0*CJ22-3.d0*CMU1*CJ23+pow2(CMU1)*CJ24))
      W0XXX=CN0*pow3(TEMR)*(CJ03+PIT26*CJ05)
      W0XXT=CN0*TEMR*(CMU1*CJ03+.5d0*CJ02+PIT26*(CMU1*CJ05+2.5d0*CJ04))
      W0XTT=CN1*(.75d0*CJ01-CMU1*CJ02+pow2(CMU1)*CJ03+PIT26*(.75d0*CJ03+3.d0*CMU1*CJ04+pow2(CMU1)*CJ05))
      return
      end subroutine SOMMERF

      subroutine SUBFERMJ(CMU1, &
        CJ00,CJ10,CJ20, &
        CJ01,CJ11,CJ21, &
        CJ02,CJ12,CJ22, &
        CJ03,CJ13,CJ23, &
        CJ04,CJ14,CJ24,CJ05)
!                                                       Version 17.11.11
! Supplement to SOMMERF

      type(auto_diff_real_2var_order1), intent(in) :: CMU1
      type(auto_diff_real_2var_order1), intent(out) :: CJ00,CJ10,CJ20,CJ01,CJ11,CJ21,CJ02,CJ12,CJ22,CJ03,CJ13,CJ23,CJ04,CJ14,CJ24,CJ05

      type(auto_diff_real_2var_order1) :: X0, X3, X5, CL, CMU

      if (CMU1.le.0.d0) stop 'SUBFERMJ: small CHI'
      CMU=1.d0+CMU1
      X0=sqrt(CMU1*(2.d0+CMU1))
      X3=pow3(X0)
      X5=pow5(X0)
      if (X0.lt.1d-4) then
         CJ00=X3/3.d0
         CJ10=0.1d0*X5
         CJ20=pow7(X0)/28.d0
      else
         CL=log(X0+CMU)
         CJ00=0.5d0*(X0*CMU-CL) ! J_{1/2}^0
         CJ10=X3/3d0-CJ00 ! J_{3/2}^0
         CJ20=(0.75d0*CMU-2d0)/3d0*X3+1.25d0*CJ00 ! J_{5/2}^0
      endif
      CJ01=X0 ! J_{1/2}^1
      CJ11=CJ01*CMU1 ! J_{3/2}^1
      CJ21=CJ11*CMU1 ! J_{5/2}^1
      CJ02=CMU/X0 ! J_{1/2}^2
      CJ12=CMU1/X0*(3.d0+2.d0*CMU1) ! J_{3/2}^2
      CJ22=pow2(CMU1)/X0*(5.d0+3.d0*CMU1) ! J_{5/2}^2
      CJ03=-1.d0/X3 ! J_{1/2}^3
      CJ13=CMU1/X3*(2.d0*pow2(CMU1)+6.d0*CMU1+3.d0)
      CJ23=pow2(CMU1)/X3*(6.d0*pow2(CMU1)+2.d1*CMU1+1.5d1)
      CJ04=3.d0*CMU/X5
      CJ14=-3.d0*CMU1/X5
      CJ24=pow2(CMU1)/X5*(6.d0*pow3(CMU1)+3.d1*pow2(CMU1)+45.d0*CMU1+15.d0)
      CJ05=(-12.d0*pow2(CMU1)-24.d0*CMU1-15.d0)/(X5*pow2(X0))
      return
      end subroutine SUBFERMJ

      subroutine FERMI10(X,XMAX,FP,FM)
!                                                       Version 20.01.10
! Fermi distribution function and its 3 derivatives
! Input: X - argument f(x)
!        XMAX - max|X| where it is assumed that 0 < f(x) < 1.
! Output: FP = f(x)
!         FM = 1-f(x)
      type(auto_diff_real_2var_order1), intent(in) :: X
      type(auto_diff_real_2var_order1) :: XMAX  ! not sure if this side-effect is desired
      type(auto_diff_real_2var_order1), intent(out) :: FP, FM
      
      if (XMAX.lt.3.d0) XMAX = 3d0 !stop 'FERMI: XMAX'
      if (X.gt.XMAX) then
         FP=0.d0
         FM=1.d0
      elseif (X.lt.-XMAX) then
         FP=1.d0
         FM=0.d0
      else
         FP=1.d0/(exp(X)+1.d0)
         FM=1.d0-FP
      endif
      return
      end subroutine FERMI10

! ==============  ELECTRON EXCHANGE AND CORRELATION   ================ *
      subroutine EXCOR7(RS,GAME,FXC,UXC,PXC,CVXC,SXC,PDTXC,PDRXC)
!                                                       Version 09.06.07
! Exchange-correlation contribution for the electron gas
! Stems from TANAKA1 v.03.03.96. Added derivatives.
! Input: RS - electron density parameter =electron-sphere radius in a.u.
!        GAME - electron Coulomb coupling parameter
! Output: FXC - excess free energy of e-liquid per kT per one electron
!             according to Tanaka & Ichimaru 85-87 and Ichimaru 93
!         UXC - internal energy contr.[per 1 electron, kT]
!         PXC - pressure contribution divided by (n_e kT)
!         CVXC - heat capacity divided by N_e k
!         SXC - entropy divided by N_e k
!         PDTXC,PDRXC = PXC + d ln PXC / d ln(T,\rho)
      use pc_support
      type(auto_diff_real_2var_order1), intent(in) :: RS, GAME
      type(auto_diff_real_2var_order1), intent(out) :: FXC,UXC,PXC,CVXC,SXC,PDTXC,PDRXC

      type(auto_diff_real_2var_order1) :: THETA, SQTH, THETA2, THETA3, THETA4, EXP1TH
      type(auto_diff_real_2var_order1) :: CHT1, SHT1, CHT2, SHT2
      type(auto_diff_real_2var_order1) :: T1, T2, T1DH, T1DHH, T2DH, T2DHH
      type(auto_diff_real_2var_order1) :: A0, A0DH, A0DHH, A1, A1DH, A1DHH, A, AH, ADH, ADHH
      type(auto_diff_real_2var_order1) :: B0, B0DH, B0DHH, B1, B1DH, B1DHH, B, BH, BDH, BDHH
      type(auto_diff_real_2var_order1) :: C, CDH, CDHH, C3, C3DH, C3DHH
      type(auto_diff_real_2var_order1) :: D0, D0DH, D0DHH, D1, D1DH, D1DHH, D, DH, DDH, DDHH
      type(auto_diff_real_2var_order1) :: E0, E0DH, E0DHH, E1, E1DH, E1DHH, E, EH, EDH, EDHH
      type(auto_diff_real_2var_order1) :: DISCR, DIDH, DIDHH
      type(auto_diff_real_2var_order1) :: S1, S1H, S1DH, S1DHH, S1DG, S1DHG
      type(auto_diff_real_2var_order1) :: B2, B2DH, B2DHH, SQGE, B3, B3DH, B3DHH
      type(auto_diff_real_2var_order1) :: S2, S2H, S2DH, S2DHH, S2DG, S2DGG, S2DHG
      type(auto_diff_real_2var_order1) :: R3, R3H, R3DH, R3DHH, R3DG, R3DGG, R3DHG
      type(auto_diff_real_2var_order1) :: S3, S3DH, S3DHH, S3DG, S3DGG, S3DHG
      type(auto_diff_real_2var_order1) :: B4, B4DH, B4DHH
      type(auto_diff_real_2var_order1) :: C4, C4DH, C4DHH, C4DG, C4DGG, C4DHG
      type(auto_diff_real_2var_order1) :: S4A, S4AH, S4ADH, S4ADHH, S4ADG, S4ADGG, S4ADHG
      type(auto_diff_real_2var_order1) :: S4B, S4BDH, S4BDHH, UP1, DN1, UP2, DN2
      type(auto_diff_real_2var_order1) :: S4C, S4CDH, S4CDHH, S4CDG, S4CDGG, S4CDHG
      type(auto_diff_real_2var_order1) :: S4, S4DH, S4DHH, S4DG, S4DGG, S4DHG
      type(auto_diff_real_2var_order1) :: FXCDH, FXCDHH, FXCDG, FXCDGG, FXCDHG
      type(auto_diff_real_2var_order1) :: PDLH, PDLG

      real(dp) :: &
         xFXC, dFXC_dlnRS, dFXC_dlnGAME, &
         xUXC, dUXC_dlnRS, dUXC_dlnGAME, &
         xPXC, dPXC_dlnRS, dPXC_dlnGAME, &
         xCVXC, dCVXC_dlnRS, dCVXC_dlnGAME, &
         xSXC, dSXC_dlnRS, dSXC_dlnGAME, &
         xPDTXC, dPDTXC_dlnRS, dPDTXC_dlnGAME, &
         xPDRXC, dPDRXC_dlnRS, dPDRXC_dlnGAME
      real(dp) :: dlnRs_dT, dlnRs_dRho, dlnGAME_dT, dlnGAME_dRho
      integer :: ierr
      logical :: skip
      logical, parameter :: use_EXCOR7_table = .true.
      logical, parameter :: debug_EXCOR7_table = .false.

      if (use_EXCOR7_table .or. debug_EXCOR7_table) then
         ierr = 0
         call get_EXCOR7(RS%val, GAME%val, &
            xFXC, dFXC_dlnRS, dFXC_dlnGAME, &
            xUXC, dUXC_dlnRS, dUXC_dlnGAME, &
            xPXC, dPXC_dlnRS, dPXC_dlnGAME, &
            xCVXC, dCVXC_dlnRS, dCVXC_dlnGAME, &
            xSXC, dSXC_dlnRS, dSXC_dlnGAME, &
            xPDTXC, dPDTXC_dlnRS, dPDTXC_dlnGAME, &
            xPDRXC, dPDRXC_dlnRS, dPDRXC_dlnGAME, skip, ierr)
         if (ierr /= 0) return
      else
         skip = .true.
      endif

      if (.not. use_EXCOR7_table .or. debug_EXCOR7_table .or. skip) then
      
         THETA=0.543d0*RS/GAME ! non-relativistic degeneracy parameter
         SQTH=sqrt(THETA)
         THETA2=THETA*THETA
         THETA3=THETA2*THETA
         THETA4=THETA3*THETA
         if (THETA.gt..007d0) then
            CHT1=cosh(1.d0/THETA)
            SHT1=sinh(1.d0/THETA)
            CHT2=cosh(1.d0/SQTH)
            SHT2=sinh(1.d0/SQTH)
            T1=SHT1/CHT1 ! dtanh(1.d0/THETA)
            T2=SHT2/CHT2 ! dtanh(1./sqrt(THETA))
            T1DH=-1.d0/((THETA*CHT1)*(THETA*CHT1)) ! d T1 / d\theta
            T1DHH=2.d0/pow3(THETA*CHT1)*(CHT1-SHT1/THETA)
            T2DH=-0.5d0*SQTH/((THETA*CHT2)*(THETA*CHT2))
            T2DHH=(0.75d0*SQTH*CHT2-0.5d0*SHT2)/pow3(THETA*CHT2)
         else
            T1=1.d0
            T2=1.d0
            T1DH=0.d0
            T2DH=0.d0
            T1DHH=0.d0
            T2DHH=0.d0
         endif
         A0=0.75d0+3.04363d0*THETA2-0.09227d0*THETA3+1.7035d0*THETA4
         A0DH=6.08726d0*THETA-0.27681d0*THETA2+6.814d0*THETA3
         A0DHH=6.08726d0-0.55362d0*THETA+20.442d0*THETA2
         A1=1d0+8.31051d0*THETA2+5.1105d0*THETA4
         A1DH=16.62102d0*THETA+20.442d0*THETA3
         A1DHH=16.62102d0+61.326d0*THETA2
         A=0.610887d0*A0/A1*T1 ! HF fit of Perrot and Dharma-wardana
         AH=A0DH/A0-A1DH/A1+T1DH/T1
         ADH=A*AH
         ADHH=ADH*AH+A*(A0DHH/A0-pow2(A0DH/A0)-A1DHH/A1+pow2(A1DH/A1) &
            + T1DHH/T1-pow2(T1DH/T1))
         B0=0.341308d0+12.070873d0*THETA2+1.148889d0*THETA4
         B0DH=24.141746d0*THETA+4.595556d0*THETA3
         B0DHH=24.141746d0+13.786668d0*THETA2
         B1=1d0+10.495346d0*THETA2+1.326623d0*THETA4
         B1DH=20.990692d0*THETA+5.306492d0*THETA3
         B1DHH=20.990692d0+15.919476d0*THETA2
         B=SQTH*T2*B0/B1
         BH=0.5d0/THETA+T2DH/T2+B0DH/B0-B1DH/B1
         BDH=B*BH
         BDHH=BDH*BH+B*(-0.5d0/THETA2+T2DHH/T2-pow2(T2DH/T2) &
            + B0DHH/B0-pow2(B0DH/B0)-B1DHH/B1+pow2(B1DH/B1))
         D0=0.614925d0+16.996055d0*THETA2+1.489056d0*THETA4
         D0DH=33.99211d0*THETA+5.956224d0*THETA3
         D0DHH=33.99211d0+17.868672d0*THETA2
         D1=1d0+10.10935d0*THETA2+1.22184d0*THETA4
         D1DH=20.2187d0*THETA+4.88736d0*THETA3
         D1DHH=20.2187d0+14.66208d0*THETA2
         D=SQTH*T2*D0/D1
         DH=0.5d0/THETA+T2DH/T2+D0DH/D0-D1DH/D1
         DDH=D*DH
         DDHH=DDH*DH+D*(-0.5d0/THETA2+T2DHH/T2-pow2(T2DH/T2) &
            + D0DHH/D0-pow2(D0DH/D0)-D1DHH/D1+pow2(D1DH/D1))
         E0=0.539409d0+2.522206d0*THETA2+0.178484d0*THETA4
         E0DH=5.044412d0*THETA+0.713936d0*THETA3
         E0DHH=5.044412d0+2.141808d0*THETA2
         E1=1d0+2.555501d0*THETA2+0.146319d0*THETA4
         E1DH=5.111002d0*THETA+0.585276d0*THETA3
         E1DHH=5.111002d0+1.755828d0*THETA2
         E=THETA*T1*E0/E1
         EH=1.d0/THETA+T1DH/T1+E0DH/E0-E1DH/E1
         EDH=E*EH
         EDHH=EDH*EH+E*(T1DHH/T1-pow2(T1DH/T1)+E0DHH/E0-pow2(E0DH/E0) &
            - E1DHH/E1+pow2(E1DH/E1)-1.0d0/THETA2)
         EXP1TH=exp(-1.d0/THETA)
         C=(0.872496d0+0.025248d0*EXP1TH)*E
         CDH=0.025248d0*EXP1TH/THETA2*E+C*EDH/E
         CDHH=0.025248d0*EXP1TH/THETA2*(EDH+(1.0d0-2.0d0*THETA)/THETA2*E) &
            + CDH*EDH/E+C*EDHH/E-C*pow2(EDH/E)
         DISCR=SQRT(4.0d0*E-D*D)
         DIDH=0.5d0/DISCR*(4.0d0*EDH-2.0d0*D*DDH)
         DIDHH=(-pow2((2.0d0*EDH-D*DDH)/DISCR)+2.0d0*EDHH-DDH*DDH-D*DDHH)/DISCR
         S1=-C/E*GAME
         S1H=CDH/C-EDH/E
         S1DH=S1*S1H
         S1DHH=S1DH*S1H+S1*(CDHH/C-pow2(CDH/C)-EDHH/E+pow2(EDH/E))
         S1DG=-C/E ! => S1DGG=0
         S1DHG=S1DG*(CDH/C-EDH/E)
         B2=B-C*D/E
         B2DH=BDH-(CDH*D+C*DDH)/E+C*D*EDH/(E*E)
         B2DHH=BDHH-(CDHH*D+2d0*CDH*DDH+C*DDHH)/E+(2d0*(CDH*D+C*DDH-C*D*EDH/E)*EDH+C*D*EDHH)/(E*E)
         SQGE=SQRT(GAME)
         S2=-2.d0/E*B2*SQGE
         S2H=B2DH/B2-EDH/E
         S2DH=S2*S2H
         S2DHH=S2DH*S2H+S2*(B2DHH/B2-pow2(B2DH/B2)-EDHH/E+pow2(EDH/E))
         S2DG=0.5d0*S2/GAME
         S2DGG=-0.5d0*S2DG/GAME
         S2DHG=0.5d0*S2DH/GAME
         R3=E*GAME+D*SQGE+1.0d0
         R3DH=EDH*GAME+DDH*SQGE
         R3DHH=EDHH*GAME+DDHH*SQGE
         R3DG=E+0.5d0*D/SQGE
         R3DGG=-0.25d0*D/(GAME*SQGE)
         R3DHG=EDH+0.5d0*DDH/SQGE
         B3=A-C/E
         B3DH=ADH-CDH/E+C*EDH/(E*E)
         B3DHH=ADHH-CDHH/E+(2.0d0*CDH*EDH+C*EDHH)/(E*E)-2d0*C*EDH*EDH/pow3(E)
         C3=(D/E*B2-B3)/E
         C3DH=(DDH*B2+D*B2DH+B3*EDH)/(E*E)-2d0*D*B2*EDH/pow3(E)-B3DH/E
         C3DHH=(-B3DHH &
            + (DDHH*B2+2d0*DDH*B2DH+D*B2DHH+B3DH*EDH+B3*EDHH+B3DH*EDH)/E &
            - 2.0d0*((DDH*B2+D*B2DH+B3*EDH+DDH*B2+D*B2DH)*EDH+D*B2*EDHH)/(E*E) &
            + 6.0d0*D*B2*EDH*EDH/pow3(E))/E
         S3=C3*log(R3)
         S3DH=S3*C3DH/C3+C3*R3DH/R3
         S3DHH=(S3DH*C3DH+S3*C3DHH)/C3-S3*pow2(C3DH/C3) &
            + (C3DH*R3DH+C3*R3DHH)/R3-C3*pow2(R3DH/R3)
         S3DG=C3*R3DG/R3
         S3DGG=C3*(R3DGG/R3-pow2(R3DG/R3))
         S3DHG=(C3DH*R3DG+C3*R3DHG)/R3-C3*R3DG*R3DH/(R3*R3)
         B4=2.d0-D*D/E
         B4DH=EDH*(D/E)*(D/E)-2d0*D*DDH/E
         B4DHH=EDHH*(D/E)*(D/E)+2d0*EDH*(D/E)*(D/E)*(DDH/D-EDH/E) &
            - 2d0*(DDH*DDH+D*DDHH)/E+2d0*D*DDH*EDH/(E*E)
         C4=2d0*E*SQGE+D
         C4DH=2d0*EDH*SQGE+DDH
         C4DHH=2d0*EDHH*SQGE+DDHH
         C4DG=E/SQGE
         C4DGG=-0.5d0*E/(GAME*SQGE)
         C4DHG=EDH/SQGE
         S4A=2.0d0/E/DISCR
         S4AH=EDH/E+DIDH/DISCR
         S4ADH=-S4A*S4AH
         S4ADHH=-S4ADH*S4AH - S4A*(EDHH/E-(EDH/E)*(EDH/E)+DIDHH/DISCR-pow2(DIDH/DISCR))
         S4B=D*B3+B4*B2
         S4BDH=DDH*B3+D*B3DH+B4DH*B2+B4*B2DH
         S4BDHH=DDHH*B3+2d0*DDH*B3DH+D*B3DHH+B4DHH*B2+2d0*B4DH*B2DH+B4*B2DHH
         S4C=atan(C4/DISCR)-atan(D/DISCR)
         UP1=C4DH*DISCR-C4*DIDH
         DN1=DISCR*DISCR+C4*C4
         UP2=DDH*DISCR-D*DIDH
         DN2=DISCR*DISCR+D*D
         S4CDH=UP1/DN1-UP2/DN2
         S4CDHH=(C4DHH*DISCR-C4*DIDHH)/DN1 &
            - UP1*2d0*(DISCR*DIDH+C4*C4DH)/(DN1*DN1) &
            - (DDHH*DISCR-D*DIDHH)/DN2+UP2*2d0*(DISCR*DIDH+D*DDH)/(DN2*DN2)
         S4CDG=C4DG*DISCR/DN1
         S4CDGG=C4DGG*DISCR/DN1-2d0*C4*DISCR*pow2(C4DG/DN1)
         S4CDHG=(C4DHG*DISCR+C4DG*DIDH-C4DG*DISCR/DN1*2d0*(DISCR*DIDH+C4*C4DH))/DN1
         S4=S4A*S4B*S4C
         S4DH=S4ADH*S4B*S4C+S4A*S4BDH*S4C+S4A*S4B*S4CDH
         S4DHH=S4ADHH*S4B*S4C+S4A*S4BDHH*S4C+S4A*S4B*S4CDHH &
            + 2d0*(S4ADH*S4BDH*S4C+S4ADH*S4B*S4CDH+S4A*S4BDH*S4CDH)
         S4DG=S4A*S4B*S4CDG
         S4DGG=S4A*S4B*S4CDGG
         S4DHG=S4A*S4B*S4CDHG+S4CDG*(S4ADH*S4B+S4A*S4BDH)
         FXC=S1+S2+S3+S4
         FXCDH=S1DH+S2DH+S3DH+S4DH
         FXCDG=S1DG+S2DG+S3DG+S4DG
         FXCDHH=S1DHH+S2DHH+S3DHH+S4DHH
         FXCDGG=S2DGG+S3DGG+S4DGG
         FXCDHG=S1DHG+S2DHG+S3DHG+S4DHG
         PXC=(GAME*FXCDG-2d0*THETA*FXCDH)/3.d0
         UXC=GAME*FXCDG-THETA*FXCDH
         SXC=(GAME*S2DG-S2+GAME*S3DG-S3+S4A*S4B*(GAME*S4CDG-S4C))-THETA*FXCDH
         if (abs(SXC).lt.1.d-9*abs(THETA*FXCDH)) SXC=0.d0 ! accuracy loss
         CVXC=2d0*THETA*(GAME*FXCDHG-FXCDH)-THETA*THETA*FXCDHH-GAME*GAME*FXCDGG
         if (abs(CVXC).lt.1.d-9*abs(GAME*GAME*FXCDGG)) CVXC=0.d0 ! accuracy
         PDLH=THETA*(GAME*FXCDHG-2d0*FXCDH-2d0*THETA*FXCDHH)/3.d0
         PDLG=GAME*(FXCDG+GAME*FXCDGG-2d0*THETA*FXCDHG)/3.d0
         PDRXC=PXC+(PDLG-2d0*PDLH)/3.d0
         PDTXC=GAME*(THETA*FXCDHG-GAME*FXCDGG/3.d0)-THETA*(FXCDH/0.75d0+THETA*FXCDHH/1.5d0)

      endif


      if (debug_EXCOR7_table .and. .not. skip) then

         if (.not. check1(FXC, xFXC, 'FXC')) return
         if (.not. check1(UXC, xUXC, 'UXC')) return
         if (.not. check1(PXC, xPXC, 'PXC')) return
         if (.not. check1(CVXC, xCVXC, 'CVXC')) return
         if (.not. check1(SXC, xSXC, 'SXC')) return
         if (.not. check1(PDTXC, xPDTXC, 'PDTXC')) return
         if (.not. check1(PDRXC, xPDRXC, 'PDRXC')) return

      endif

      if (use_EXCOR7_table .and. .not. skip) then

         ! values
         FXC = xFXC
         UXC = xUXC
         PXC = xPXC
         CVXC = xCVXC
         SXC = xSXC
         PDTXC = xPDTXC
         PDRXC = xPDRXC

         ! dT (val1) derivatives
         ! numerically, RS% d1val1 / RS% val is not quite 0
         dlnRs_dT = RS% d1val1 / RS% val
         dlnGAME_dT = GAME% d1val1 / GAME% val

         FXC% d1val1 = dFXC_dlnRS * dlnRS_dT + dFXC_dlnGAME * dlnGAME_dT
         UXC% d1val1 = dUXC_dlnRS * dlnRS_dT + dUXC_dlnGAME * dlnGAME_dT
         PXC% d1val1 = dPXC_dlnRS * dlnRS_dT + dPXC_dlnGAME * dlnGAME_dT
         CVXC% d1val1 = dCVXC_dlnRS * dlnRS_dT + dCVXC_dlnGAME * dlnGAME_dT
         SXC% d1val1 = dSXC_dlnRS * dlnRS_dT + dSXC_dlnGAME * dlnGAME_dT
         PDTXC% d1val1 = dPDTXC_dlnRS * dlnRS_dT + dPDTXC_dlnGAME * dlnGAME_dT
         PDRXC% d1val1 = dPDRXC_dlnRS * dlnRS_dT + dPDRXC_dlnGAME * dlnGAME_dT

         ! dRho (val2) derivatives
         dlnRs_dRho = RS% d1val2 / RS% val
         dlnGAME_dRho = GAME% d1val2 / GAME% val

         FXC% d1val2 = dFXC_dlnRS * dlnRS_dRho + dFXC_dlnGAME * dlnGAME_dRho
         UXC% d1val2 = dUXC_dlnRS * dlnRS_dRho + dUXC_dlnGAME * dlnGAME_dRho
         PXC% d1val2 = dPXC_dlnRS * dlnRS_dRho + dPXC_dlnGAME * dlnGAME_dRho
         CVXC% d1val2 = dCVXC_dlnRS * dlnRS_dRho + dCVXC_dlnGAME * dlnGAME_dRho
         SXC% d1val2 = dSXC_dlnRS * dlnRS_dRho + dSXC_dlnGAME * dlnGAME_dRho
         PDTXC% d1val2 = dPDTXC_dlnRS * dlnRS_dRho + dPDTXC_dlnGAME * dlnGAME_dRho
         PDRXC% d1val2 = dPDRXC_dlnRS * dlnRS_dRho + dPDRXC_dlnGAME * dlnGAME_dRho

      end if
      
      contains
      
      logical function check1(v, xv, str)
         type(auto_diff_real_2var_order1), intent(in) :: v
         real(dp), intent(in) :: xv
         character (len=*), intent(in) :: str
         real(dp) :: val
         real(dp), parameter :: atol = 1d-8, rtol = 1d-6
         include 'formats'
         val = v%val
         check1 = .false.
         if (is_bad(xv)) then
            write(*,*) 'is_bad ' // trim(str), xv, val, RS, GAME
            return
         end if
         if (abs(val - xv) > atol + rtol*max(abs(val),abs(xv))) then
            write(*,*) 'rel mismatch ' // trim(str), &
               (val - xv)/max(abs(val),abs(xv),1d-99), &
               xv, val, RS%val, GAME%val
            stop 'EXCOR7'
            return
         end if
         check1 = .true.
      end function check1


      end subroutine EXCOR7

! ======================  AUXILIARY SUBROUTINES   ==================== *
      subroutine FERINV7(F,N,X,XDF,XDFF) ! Inverse Fermi intergals
!                                                       Version 24.05.07
! X_q(f)=F^{-1}_q(f) : H.M.Antia 93 ApJS 84, 101
! q=N-1/2=-1/2,1/2,3/2,5/2 (N=0,1,2,3)
! Input: F - argument, N=q+1/2
! Output: X=X_q, XDF=dX/df, XDFF=d^2 X / df^2
! Relative error: N = 0     1      2      3
!        for X:    3.e-9, 4.2e-9, 2.3e-9, 6.2e-9
! jump at f=4:
!         for XDF: 6.e-7, 5.4e-7, 9.6e-8, 3.1e-7
!       for XDFF: 4.7e-5, 4.8e-5, 2.3e-6, 1.5e-6
      type(auto_diff_real_2var_order1) :: F ! can be modified, not sure if this is an intended side-effect
      integer, intent(in) :: N
      type(auto_diff_real_2var_order1), intent(out) :: X, XDF, XDFF
      integer :: I
      type(auto_diff_real_2var_order1) :: P, T, T1, T2, UP, UP1, UP2, DOWN, DOWN1, DOWN2
      type(auto_diff_real_2var_order1) :: R, R1, R2, RT
      
      ! The next four are really parameters but there isn't a clean way to initialize them
      ! at declaration time. - Adam Jermyn 4/2/2020
      real(dp) :: A(0:5,0:3) ! read only after initialization
      real(dp) :: B(0:6,0:3) ! read only after initialization
      real(dp) :: C(0:6,0:3) ! read only after initialization
      real(dp) :: D(0:6,0:3) ! read only after initialization

      integer, parameter :: LA(0:3) = (/5,4,3,2/)
      integer, parameter :: LB(0:3) = (/6,3,4,3/)
      integer, parameter :: LD(0:3) = (/6,5,5,6/)

         A(0:5,0) = (/ &
            -1.570044577033d4,1.001958278442d4,-2.805343454951d3, &
                  4.121170498099d2,-3.174780572961d1,1.d0/) ! X_{-1/2}
         A(0:5,1) = (/ &
            1.999266880833d4,5.702479099336d3,6.610132843877d2, &
                  3.818838129486d1,1.d0,0d0/) ! X_{1/2}
         A(0:5,2) = (/ &
            1.715627994191d2,1.125926232897d2,2.056296753055d1, &
                  1.d0,0d0,0d0/)
         A(0:5,3) = (/ &
            2.138969250409d2,3.539903493971d1,1.d0,0d0,0d0,0d0/) ! X_{5/2}
         B(0:6,0) = (/ &
            -2.782831558471d4,2.886114034012d4,-1.274243093149d4, &
                  3.063252215963d3,-4.225615045074d2,3.168918168284d1, &
                  -1.008561571363d0/) ! X_{-1/2}
         B(0:6,1) = (/ &
            1.771804140488d4,-2.014785161019d3,9.130355392717d1, &
                  -1.670718177489d0,0d0,0d0,0d0/) ! X_{1/2}
         B(0:6,2) = (/ &
            2.280653583157d2,1.193456203021d2,1.16774311354d1, &
                  -3.226808804038d-1,3.519268762788d-3,0d0,0d0/) ! X_{3/2}
         B(0:6,3) = (/ &
            7.10854551271d2,9.873746988121d1,1.067755522895d0, &
                  -1.182798726503d-2,0d0,0d0,0d0/) ! X_{5/2}
         C(0:6,0) = (/ &
            2.206779160034d-8,-1.437701234283d-6,6.103116850636d-5, &
            -1.169411057416d-3,1.814141021608d-2,-9.588603457639d-2, &
            1.d0/)
         C(0:6,1) = (/ &
             -1.277060388085d-2,7.187946804945d-2,-4.262314235106d-1, &
            4.997559426872d-1,-1.285579118012d0,-3.930805454272d-1, &
            1.d0/)
         C(0:6,2) = (/ &
             -6.321828169799d-3,-2.183147266896d-2,-1.05756279932d-1, &
             -4.657944387545d-1,-5.951932864088d-1,3.6844711771d-1, &
            1.d0/)
         C(0:6,3) = (/ &
           -3.312041011227d-2,1.315763372315d-1,-4.820942898296d-1, &
             5.099038074944d-1,5.49561349863d-1,-1.498867562255d0, &
            1.d0/)
         D(0:6,0) = (/ &
             8.827116613576d-8,-5.750804196059d-6,2.429627688357d-4, &
             -4.601959491394d-3,6.932122275919d-2,-3.217372489776d-1, &
             3.124344749296d0/) ! X_{-1/2}
         D(0:6,1) = (/ &
          -9.745794806288d-3,5.485432756838d-2,-3.29946624326d-1, &
             4.077841975923d-1,-1.145531476975d0,-6.067091689181d-2, &
            0d0/)
         D(0:6,2) = (/ &
          -4.381942605018d-3,-1.5132365041d-2,-7.850001283886d-2, &
           -3.407561772612d-1,-5.074812565486d-1,-1.387107009074d-1, &
            0d0/)
         D(0:6,3) = (/ &
          -2.315515517515d-2,9.198776585252d-2,-3.835879295548d-1, &
             5.415026856351d-1,-3.847241692193d-1,3.739781456585d-2, &
             -3.008504449098d-2/) ! X_{5/2}
      
      if (N.lt.0d0 .or.N.gt.3d0) stop 'FERINV7: Invalid subscript'
      if (F.le.0.d0) F = 1d-99 !stop 'FERINV7: Non-positive argument'
      if (F.lt.4.d0) then
         T=F
         UP=0.d0
         UP1=0.d0
         UP2=0.d0
         DOWN=0.d0
         DOWN1=0.d0
         DOWN2=0.d0
         do I=LA(N),0,-1
            UP=UP*T+A(I,N)
           if (I.ge.1) UP1=UP1*T+A(I,N)*I
           if (I.ge.2) UP2=UP2*T+A(I,N)*I*(I-1)
         enddo
         do I=LB(N),0,-1
            DOWN=DOWN*T+B(I,N)
           if (I.ge.1) DOWN1=DOWN1*T+B(I,N)*I
           if (I.ge.2) DOWN2=DOWN2*T+B(I,N)*I*(I-1)
         enddo
         X=log(T*UP/DOWN)
         XDF=1.d0/T+UP1/UP-DOWN1/DOWN
         XDFF=-1.d0/(T*T)+UP2/UP-pow2(UP1/UP)-DOWN2/DOWN+pow2(DOWN1/DOWN)
      else
         P=-1.d0/(.5d0+N) ! = -1/(1+\nu) = power index
         T=pow(F,P) ! t - argument of the rational fraction
         T1=P*T/F ! dt/df
         T2=P*(P-1.d0)*T/(F*F) ! d^2 t / df^2
         UP=0.d0
         UP1=0.d0
         UP2=0.d0
         DOWN=0.d0
         DOWN1=0.d0
         DOWN2=0.d0
         do I=6,0,-1
            UP=UP*T+C(I,N)
           if (I.ge.1) UP1=UP1*T+C(I,N)*I
           if (I.ge.2) UP2=UP2*T+C(I,N)*I*(I-1)
         enddo
         do I=LD(N),0,-1
            DOWN=DOWN*T+D(I,N)
           if (I.ge.1) DOWN1=DOWN1*T+D(I,N)*I
           if (I.ge.2) DOWN2=DOWN2*T+D(I,N)*I*(I-1)
         enddo
         R=UP/DOWN
         R1=(UP1-UP*DOWN1/DOWN)/DOWN ! dR/dt
         R2=(UP2-(2.0d0*UP1*DOWN1+UP*DOWN2)/DOWN+2.d0*UP*pow2(DOWN1/DOWN))/DOWN
         X=R/T
         RT=(R1-R/T)/T
         XDF=T1*RT
         XDFF=T2*RT+T1*T1*(R2-2d0*RT)/T
      endif
      return
      end subroutine FERINV7

      subroutine BLIN9(TEMP,CHI, &
        W0,W0DX,W0DT,W0DXX,W0DTT,W0DXT, &
        W1,W1DX,W1DT,W1DXX,W1DTT,W1DXT, &
        W2,W2DX,W2DT,W2DXX,W2DTT,W2DXT, &
        W0XXX,W0XTT,W0XXT)
!                                                       Version 21.01.10
! Stems from BLIN8 v.24.12.08
! Difference - smooth matching of different CHI ranges
! Input: TEMP=T/mc^2; CHI=(\mu-mc^2)/T
! Output: Wk - Fermi-Dirac integral of the order k+1/2
!         WkDX=dWk/dCHI, WkDT = dWk/dT, WkDXX=d^2 Wk / d CHI^2,
!         WkDTT=d^2 Wk / d T^2, WkDXT=d^2 Wk /dCHIdT,
!         W0XXX=d^3 W0 / d CHI^3, W0XTT=d^3 W0 /(d CHI d^2 T),
!         W0XXT=d^3 W0 /dCHI^2 dT
      type(auto_diff_real_2var_order1), intent(in) :: TEMP,CHI
      type(auto_diff_real_2var_order1), intent(out) :: W0,W0DX,W0DT,W0DXX,W0DTT,W0DXT
      type(auto_diff_real_2var_order1), intent(out) :: W1,W1DX,W1DT,W1DXX,W1DTT,W1DXT
      type(auto_diff_real_2var_order1), intent(out) :: W2,W2DX,W2DT,W2DXX,W2DTT,W2DXT
      type(auto_diff_real_2var_order1), intent(out) :: W0XXX,W0XTT,W0XXT

      type(auto_diff_real_2var_order1) :: X1, X2, FP, FM
      type(auto_diff_real_2var_order1) :: W0a,W0DXa,W0DTa,W0DXXa,W0DTTa,W0DXTa
      type(auto_diff_real_2var_order1) :: W1a,W1DXa,W1DTa,W1DXXa,W1DTTa,W1DXTa
      type(auto_diff_real_2var_order1) :: W2a,W2DXa,W2DTa,W2DXXa,W2DTTa,W2DXTa
      type(auto_diff_real_2var_order1) :: W0XXXa,W0XTTa,W0XXTa
      type(auto_diff_real_2var_order1) :: W0b,W0DXb,W0DTb,W0DXXb,W0DTTb,W0DXTb
      type(auto_diff_real_2var_order1) :: W1b,W1DXb,W1DTb,W1DXXb,W1DTTb,W1DXTb
      type(auto_diff_real_2var_order1) :: W2b,W2DXb,W2DTb,W2DXXb,W2DTTb,W2DXTb
      type(auto_diff_real_2var_order1) :: W0XXXb,W0XTTb,W0XXTb
      
      type(auto_diff_real_2var_order1) :: CHI1
      type(auto_diff_real_2var_order1) :: CHI2
      type(auto_diff_real_2var_order1) :: XMAX
      type(auto_diff_real_2var_order1) :: DCHI1
      type(auto_diff_real_2var_order1) :: DCHI2
      type(auto_diff_real_2var_order1) :: XSCAL1
      type(auto_diff_real_2var_order1) :: XSCAL2

      CHI1=0.6d0
      CHI2=14.d0
      XMAX=30.d0
      DCHI1=0.1d0
      DCHI2=CHI2-CHI1-DCHI1
      XSCAL1=XMAX/DCHI1
      XSCAL2=XMAX/DCHI2
      
      X1=(CHI-CHI1)*XSCAL1
      X2=(CHI-CHI2)*XSCAL2
      if (X1.lt.-XMAX) then
         call BLIN9a(TEMP,CHI, &
           W0,W0DX,W0DT,W0DXX,W0DTT,W0DXT, &
           W1,W1DX,W1DT,W1DXX,W1DTT,W1DXT, &
           W2,W2DX,W2DT,W2DXX,W2DTT,W2DXT, &
           W0XXX,W0XTT,W0XXT)
      elseif (X2.lt.XMAX) then ! match two fits
        if (X1.lt.XMAX) then ! match fits "a" and "b"
           call FERMI10(X1,XMAX,FP,FM)
           call BLIN9a(TEMP,CHI, &
             W0a,W0DXa,W0DTa,W0DXXa,W0DTTa,W0DXTa, &
             W1a,W1DXa,W1DTa,W1DXXa,W1DTTa,W1DXTa, &
             W2a,W2DXa,W2DTa,W2DXXa,W2DTTa,W2DXTa, &
             W0XXXa,W0XTTa,W0XXTa)
           call BLIN9b(TEMP,CHI, &
             W0b,W0DXb,W0DTb,W0DXXb,W0DTTb,W0DXTb, &
             W1b,W1DXb,W1DTb,W1DXXb,W1DTTb,W1DXTb, &
             W2b,W2DXb,W2DTb,W2DXXb,W2DTTb,W2DXTb, &
             W0XXXb,W0XTTb,W0XXTb)
        else ! match fits "b" and "c"
           call FERMI10(X2,XMAX,FP,FM)
           call BLIN9b(TEMP,CHI, &
             W0a,W0DXa,W0DTa,W0DXXa,W0DTTa,W0DXTa, &
             W1a,W1DXa,W1DTa,W1DXXa,W1DTTa,W1DXTa, &
             W2a,W2DXa,W2DTa,W2DXXa,W2DTTa,W2DXTa, &
             W0XXXa,W0XTTa,W0XXTa)
           call BLIN9c(TEMP,CHI, &
             W0b,W0DXb,W0DTb,W0DXXb,W0DTTb,W0DXTb, &
             W1b,W1DXb,W1DTb,W1DXXb,W1DTTb,W1DXTb, &
             W2b,W2DXb,W2DTb,W2DXXb,W2DTTb,W2DXTb, &
             W0XXXb,W0XTTb,W0XXTb)
        endif
         W0=W0a*FP+W0b*FM
         W0DX=W0DXa*FP+W0DXb*FM !! +(W0a-W0b)*F1
         W0DT=W0DTa*FP+W0DTb*FM
         W0DXX=W0DXXa*FP+W0DXXb*FM !! +2.d0*(W0DXa-W0DXb)*F1+(W0a-W0b)*F2
         W0DTT=W0DTTa*FP+W0DTTb*FM
         W0DXT=W0DXTa*FP+W0DXTb*FM !! +(W0DTa-W0DTb)*F1
         W0XXX=W0XXXa*FP+W0XXXb*FM !! +3.d0*(W0DXXa-W0DXXb)*F1+3.d0*(W0DXa-W0DXb)*F2+(W0a-W0b)*F3
         W0XTT=W0XTTa*FP+W0XTTb*FM !! +(W0DTTa-W0DTTb)*F1
         W0XXT=W0XXTa*FP+W0XXTb*FM !! +2.d0*(W0DXTa-W0DXTb)*F1+(W0DTa-W0DTb)*F2
         W1=W1a*FP+W1b*FM
         W1DX=W1DXa*FP+W1DXb*FM !! +(W1a-W1b)*F1
         W1DT=W1DTa*FP+W1DTb*FM
         W1DXX=W1DXXa*FP+W1DXXb*FM !! +2.d0*(W1DXa-W1DXb)*F1+(W1a-W1b)*F2
         W1DTT=W1DTTa*FP+W1DTTb*FM
         W1DXT=W1DXTa*FP+W1DXTb*FM !! +(W1DTa-W1DTb)*F1
         W2=W2a*FP+W2b*FM
         W2DX=W2DXa*FP+W2DXb*FM !! +(W2a-W2b)*F1
         W2DT=W2DTa*FP+W2DTb*FM
         W2DXX=W2DXXa*FP+W2DXXb*FM !! +2.d0*(W2DXa-W2DXb)*F1+(W2a-W2b)*F2
         W2DTT=W2DTTa*FP+W2DTTb*FM
         W2DXT=W2DXTa*FP+W2DXTb*FM !! 
      else
         call BLIN9c(TEMP,CHI, &
           W0,W0DX,W0DT,W0DXX,W0DTT,W0DXT, &
           W1,W1DX,W1DT,W1DXX,W1DTT,W1DXT, &
           W2,W2DX,W2DT,W2DXX,W2DTT,W2DXT, &
           W0XXX,W0XTT,W0XXT)
      endif
      return
      end subroutine BLIN9

      subroutine BLIN9a(TEMP,CHI, &
        W0,W0DX,W0DT,W0DXX,W0DTT,W0DXT, &
        W1,W1DX,W1DT,W1DXX,W1DTT,W1DXT, &
        W2,W2DX,W2DT,W2DXX,W2DTT,W2DXT, &
        W0XXX,W0XTT,W0XXT)
!                                                       Version 19.01.10
! First part of BILN9: small CHI. Stems from BLIN9 v.24.12.08
      type(auto_diff_real_2var_order1), intent(in) :: TEMP,CHI
      type(auto_diff_real_2var_order1), intent(out) :: W0,W0DX,W0DT,W0DXX,W0DTT,W0DXT
      type(auto_diff_real_2var_order1), intent(out) :: W1,W1DX,W1DT,W1DXX,W1DTT,W1DXT
      type(auto_diff_real_2var_order1), intent(out) :: W2,W2DX,W2DT,W2DXX,W2DTT,W2DXT
      type(auto_diff_real_2var_order1), intent(out) :: W0XXX,W0XTT,W0XXT

      type(auto_diff_real_2var_order1) :: W,WDX,WDT,WDXX,WDTT,WDXT,WDXXX,WDXTT,WDXXT
      type(auto_diff_real_2var_order1) :: ECHI, SQ, DN

      integer :: I, J, K

      ! The next three are really parameters but there isn't a clean way to initialize them
      ! at declaration time. - Adam Jermyn 4/2/2020
      real(dp) :: AC(5,0:2) ! read only after initialization
      real(dp) :: AU(5,0:2) ! read only after initialization
      real(dp) :: AA(5,0:2) ! read only after initialization

      AC(1:5,0) = (/ &
         0.37045057d0, .41258437d0, &
           9.777982d-2, 5.3734153d-3, 3.8746281d-5/) ! c_i^0
      AC(1:5,1) = (/ &
           .39603109d0, .69468795d0,  &
           .22322760d0, 1.5262934d-2, 1.3081939d-4/) ! c_i^1
      AC(1:5,2) = (/ &
           .76934619d0, 1.7891437d0,  &
           .70754974d0, 5.6755672d-2, 5.5571480d-4/) ! c_i^2
      AU(1:5,0) = (/ &
         0.43139881d0, 1.7597537d0,  &
           4.1044654d0, 7.7467038d0, 13.457678d0/) ! \chi_i^0
      AU(1:5,1) = (/ &
           .81763176d0, 2.4723339d0,  &
           5.1160061d0, 9.0441465d0, 15.049882d0/) ! \chi_i^1
      AU(1:5,2) = (/ &
           1.2558461d0, 3.2070406d0,  &
           6.1239082d0, 10.316126d0, 16.597079d0/) ! \chi_i^2

     do J=0,2
        do I=1,5
           AA(I,J)=exp(-AU(I,J))
        enddo
     enddo

        do K=0,2
           W=0.d0
           WDX=0.d0
           WDT=0.d0
           WDXX=0.d0
           WDTT=0.d0
           WDXT=0.d0
           WDXXX=0.d0
           WDXTT=0.d0
           WDXXT=0.d0
             ECHI=exp(-CHI)
            do I=1,5
               SQ=sqrt(1.d0+AU(I,K)*TEMP/2.d0)
               DN=AA(I,K)+ECHI ! e^{-\chi_i}+e^{-\chi})
               W=W+AC(I,K)*SQ/DN
               WDX=WDX+AC(I,K)*SQ/(DN*DN)
               WDT=WDT+AC(I,K)*AU(I,K)/(SQ*DN)
               WDXX=WDXX+AC(I,K)*SQ*(ECHI-AA(I,K))/pow3(DN)
               WDTT=WDTT-AC(I,K)*AU(I,K)*AU(I,K)/(DN*SQ*SQ*SQ)
               WDXT=WDXT+AC(I,K)*AU(I,K)/(SQ*DN*DN)
               WDXXX=WDXXX+AC(I,K)*SQ* &
                 (ECHI*ECHI-4.d0*ECHI*AA(I,K)+AA(I,K)*AA(I,K))/(DN*DN*DN*DN)
               WDXTT=WDXTT-AC(I,K)*AU(I,K)*AU(I,K)/(DN*DN*SQ*SQ*SQ)
               WDXXT=WDXXT+AC(I,K)*AU(I,K)*(ECHI-AA(I,K))/(SQ*DN*DN*DN)
            enddo
             WDX=WDX*ECHI
             WDT=0.25d0*WDT
             WDXX=WDXX*ECHI
             WDTT=0.0625d0*WDTT
             WDXT=0.25d0*WDXT*ECHI
             WDXXX=WDXXX*ECHI
             WDXTT=0.0625d0*WDXTT*ECHI
             WDXXT=0.25d0*WDXXT*ECHI
          if (K.eq.0) then
             W0=W
             W0DX=WDX
             W0DT=WDT
             W0DXX=WDXX
             W0DTT=WDTT
             W0DXT=WDXT
             W0XXX=WDXXX
             W0XTT=WDXTT
             W0XXT=WDXXT
          elseif (K.eq.1) then
             W1=W
             W1DX=WDX
             W1DT=WDT
             W1DXX=WDXX
             W1DTT=WDTT
             W1DXT=WDXT
          else
             W2=W
             W2DX=WDX
             W2DT=WDT
             W2DXX=WDXX
             W2DTT=WDTT
             W2DXT=WDXT
          endif
        enddo ! next K
      return
      end subroutine BLIN9a

      subroutine BLIN9b(TEMP,CHI, &
        W0,W0DX,W0DT,W0DXX,W0DTT,W0DXT, &
        W1,W1DX,W1DT,W1DXX,W1DTT,W1DXT, &
        W2,W2DX,W2DT,W2DXX,W2DTT,W2DXT, &
        W0XXX,W0XTT,W0XXT)
!                                                       Version 19.01.10
! Second part of BILN9: intermediate CHI. Stems from BLIN8 v.24.12.08
      type(auto_diff_real_2var_order1), intent(in) :: TEMP
      type(auto_diff_real_2var_order1) :: CHI ! can be modified, not sure if this is an intended side-effect
      type(auto_diff_real_2var_order1), intent(out) :: W0,W0DX,W0DT,W0DXX,W0DTT,W0DXT
      type(auto_diff_real_2var_order1), intent(out) :: W1,W1DX,W1DT,W1DXX,W1DTT,W1DXT
      type(auto_diff_real_2var_order1), intent(out) :: W2,W2DX,W2DT,W2DXX,W2DTT,W2DXT
      type(auto_diff_real_2var_order1), intent(out) :: W0XXX,W0XTT,W0XXT

      type(auto_diff_real_2var_order1) :: W,WDX,WDT,WDXX,WDTT,WDXT,WDXXX,WDXTT,WDXXT
      type(auto_diff_real_2var_order1) :: SQCHI, CE, ECHI, DE, D, XICHI, DXI
      type(auto_diff_real_2var_order1) :: H, HX, HDX, HXX, HDXX, HT, HDT, HDTT
      type(auto_diff_real_2var_order1) :: HTX, HDXT, HDXXT, HDXTT, HXXX, HDXXX
      type(auto_diff_real_2var_order1) :: V, VX, VDX, VT, VDT, VXX, VDXX, VDXXX
      type(auto_diff_real_2var_order1) :: VXXT, VDTT, VXT, VDXT, VDXXT, VDXTT

      integer :: I, J, K
      real(dp), parameter :: AX(5) = &
         (/7.265351d-2, .2694608d0,  &
              .533122d0, .7868801d0, .9569313d0/) ! x_i
      real(dp), parameter :: AXI(5) = &
         (/.26356032d0, 1.4134031d0,  &
               3.5964258d0, 7.0858100d0, 12.640801d0/) ! \xi_i
      real(dp), parameter :: AH(5) = &
         (/3.818735d-2, .1256732d0,  &
              .1986308d0, .1976334d0, .1065420d0/) ! H_i
      real(dp), parameter :: AV(5) = &
         (/.29505869d0, .32064856d0, 7.3915570d-2,  &
              3.6087389d-3, 2.3369894d-5/) ! \bar{V}_i
      real(dp), parameter :: EPS=1.d-3
      
      if (CHI.lt.EPS) CHI = EPS !stop 'BLIN9b: CHI is too small'
        do K=0,2
           W=0.d0
           WDX=0.d0
           WDT=0.d0
           WDXX=0.d0
           WDTT=0.d0
           WDXT=0.d0
           WDXXX=0.d0
           WDXTT=0.d0
           WDXXT=0.d0
             SQCHI=sqrt(CHI)
            do I=1,5
               CE=AX(I)-1.d0
               ECHI=exp(CE*CHI)
               DE=1.d0+ECHI
               D=1.d0+AX(I)*CHI*TEMP/2.d0
               H=pow(CHI,K+1d0)*SQCHI*sqrt(D)/DE
               HX=(K+1.5d0)/CHI+.25d0*AX(I)*TEMP/D-ECHI*CE/DE
               HDX=H*HX
               HXX=(K+1.5d0)/(CHI*CHI)+.125d0*pow2(AX(I)*TEMP/D)+ECHI*pow2(CE/DE)
               HDXX=HDX*HX-H*HXX
               HT=.25d0*AX(I)*CHI/D
               HDT=H*HT
               HDTT=-H*HT*HT
               HTX=1.d0/CHI-.5d0*AX(I)*TEMP/D
               HDXT=HDX*HT+HDT*HTX
               HDXXT=HDXX*HT+HDX*HT*HTX+HDXT*HTX &
                  + HDT*(.25d0*pow2(AX(I)*TEMP/D)-1.d0/(CHI*CHI))
               HDXTT=HDXT*HT-HDX*.125d0*pow2(AX(I)*CHI/D)+HDTT*HTX &
                  + HDT*.5d0*AX(I)*(TEMP*.5d0*AX(I)*CHI/(D*D)-1.d0/D)
               HXXX=(2d0*K+3d0)/pow3(CHI)+.125d0*pow3(AX(I)*TEMP/D) &
                  - ECHI*(1.d0-ECHI)*pow3(CE/DE)
               HDXXX=HDXX*HX-2.d0*HDX*HXX+H*HXXX
               XICHI=AXI(I)+CHI
               DXI=1.d0+XICHI*TEMP/2.d0
               V=pow(XICHI,1d0*K)*sqrt(XICHI*DXI)
               VX=(K+.5d0)/XICHI+.25d0*TEMP/DXI
               VDX=V*VX
               VT=.25d0*XICHI/DXI
               VDT=V*VT
               VXX=(K+0.5d0)/(XICHI*XICHI)+.125d0*pow2(TEMP/DXI)
               VDXX=VDX*VX-V*VXX
               VDXXX=VDXX*VX-2d0*VDX*VXX &
                    + V*((2d0*K+1d0)/pow3(XICHI)+.125d0*pow3(TEMP/DXI))
               VXXT=(1.d0-0.5d0*TEMP*XICHI/DXI)/DXI
               VDTT=-V*VT*VT
               VXT=1.d0/XICHI-0.5d0*TEMP/DXI
               VDXT=VDT*VXT+VDX*VT
               VDXXT=VDXT*VX+VDX*.25d0*VXXT-VDT*VXX-V*.25d0*TEMP/DXI*VXXT
               VDXTT=VDTT*VXT-VDT*0.5d0*VXXT+VDXT*VT &
                  - VDX*.125d0*pow2(XICHI/DXI)
               W=W+AH(I)*pow(AX(I),K)*H+AV(I)*V
               WDX=WDX+AH(I)*pow(AX(I),K)*HDX+AV(I)*VDX
               WDT=WDT+AH(I)*pow(AX(I),K)*HDT+AV(I)*VDT
               WDXX=WDXX+AH(I)*pow(AX(I),K)*HDXX+AV(I)*VDXX
               WDTT=WDTT+AH(I)*pow(AX(I),K)*HDTT+AV(I)*VDTT
               WDXT=WDXT+AH(I)*pow(AX(I),K)*HDXT+AV(I)*VDXT
               WDXXX=WDXXX+AH(I)*pow(AX(I),K)*HDXXX+AV(I)*VDXXX
               WDXTT=WDXTT+AH(I)*pow(AX(I),K)*HDXTT+AV(I)*VDXTT
               WDXXT=WDXXT+AH(I)*pow(AX(I),K)*HDXXT+AV(I)*VDXXT
            enddo
          if (K.eq.0) then
             W0=W
             W0DX=WDX
             W0DT=WDT
             W0DXX=WDXX
             W0DTT=WDTT
             W0DXT=WDXT
             W0XXX=WDXXX
             W0XTT=WDXTT
             W0XXT=WDXXT
          elseif (K.eq.1) then
             W1=W
             W1DX=WDX
             W1DT=WDT
             W1DXX=WDXX
             W1DTT=WDTT
             W1DXT=WDXT
          else
             W2=W
             W2DX=WDX
             W2DT=WDT
             W2DXX=WDXX
             W2DTT=WDTT
             W2DXT=WDXT
          endif
        enddo ! next K
      return
      end subroutine BLIN9b

      subroutine BLIN9c(TEMP,CHI, &
        W0,W0DX,W0DT,W0DXX,W0DTT,W0DXT, &
        W1,W1DX,W1DT,W1DXX,W1DTT,W1DXT, &
        W2,W2DX,W2DT,W2DXX,W2DTT,W2DXT, &
        W0XXX,W0XTT,W0XXT)
!                                                       Version 19.01.10
! Third part of BILN9: large CHI. Stems from BLIN8 v.24.12.08
      type(auto_diff_real_2var_order1), intent(in) :: CHI, TEMP
      ! type(auto_diff_real_2var_order1) :: CHI ! can be modified, not sure if this is an intended side-effect
      type(auto_diff_real_2var_order1), intent(out) :: W0,W0DX,W0DT,W0DXX,W0DTT,W0DXT
      type(auto_diff_real_2var_order1), intent(out) :: W1,W1DX,W1DT,W1DXX,W1DTT,W1DXT
      type(auto_diff_real_2var_order1), intent(out) :: W2,W2DX,W2DT,W2DXX,W2DTT,W2DXT
      type(auto_diff_real_2var_order1), intent(out) :: W0XXX,W0XTT,W0XXT

      type(auto_diff_real_2var_order1), dimension(0:2) :: AM,AMDX,AMDT,AMDXX,AMDTT,AMDXT
      type(auto_diff_real_2var_order1) :: CNU, CHINU, F, FDX, FDXX, FDXXX, C, D, DM
      type(auto_diff_real_2var_order1) :: W, WDX, WDT, WDXX, WDXT, WDTT, WDXXX, WDXXT, WDXTT
      type(auto_diff_real_2var_order1) :: R, RX, RDX, RDT, RXX, RDXX, RDTT, RXT, RDXT
      type(auto_diff_real_2var_order1) :: RXXX, RDXXX, RXTT, RDXTT, RXXT, RDXXT
      type(auto_diff_real_2var_order1) :: FMX1, FMX2, FMX, CKM
      type(auto_diff_real_2var_order1) :: FMT1, FMT2, FMT, FMXX, FMXXX, FMTT
      type(auto_diff_real_2var_order1) :: AMDXXX, FMT1DX, FMT2DX, FMXT, FMTTX
      type(auto_diff_real_2var_order1) :: AMDXTT, FMX1DT, FMX2DT, FMXXT, AMDXXT
      type(auto_diff_real_2var_order1) :: SQ2T, SQ2T3
      type(auto_diff_real_2var_order1) :: A, ADX, ADT, ADXX, ADTT, ADXT, ADXTT, ADXXT
      type(auto_diff_real_2var_order1) :: XT1, Aln, FJ0, ASQ3, ASQ3DX, BXT, BXXT
      type(auto_diff_real_2var_order1) :: FJ0DX, FJ0DT, FJ0DXX, FJ0DTT, FJ0DXT
      type(auto_diff_real_2var_order1) :: FJ0XXX, FJ0XXT, FJ0XTT
      type(auto_diff_real_2var_order1) :: FJ1, FJ1DX, FJ1DT, FJ1DXX, FJ1DXT, FJ1DTT
      type(auto_diff_real_2var_order1) :: FJ2, FJ2DX, FJ2DT, FJ2DXX, FJ2DXT, FJ2DTT
      
      integer :: J, K
      real(dp), parameter :: PI26=PI*PI/6.

      if (CHI*TEMP.lt..1d0) then
        do K=0,2
           W=0.d0
           WDX=0.d0
           WDT=0.d0
           WDXX=0.d0
           WDTT=0.d0
           WDXT=0.d0
           WDXXX=0.d0
           WDXTT=0.d0
           WDXXT=0.d0
            do J=0,4 ! for nonrel.Fermi integrals from k+1/2 to k+4.5
               CNU=K+J+0.5d0 ! nonrelativistic Fermi integral index \nu
               CHINU=pow(CHI,1d0*(K+J))*sqrt(CHI) ! \chi^\nu
               F=CHINU*(CHI/(CNU+1.d0)+PI26*CNU/CHI & ! nonrel.Fermi
                  + .7d0*PI26*PI26*CNU*(CNU-1.d0)*(CNU-2.d0)/pow3(CHI))
               FDX=CHINU*(1d0+PI26*CNU*(CNU-1.d0)/pow2(CHI) &
                  +.7d0*PI26*PI26*CNU*(CNU-1.d0)*(CNU-2.d0)*(CNU-3.d0)/pow4(CHI))
               FDXX=CHINU/CHI*CNU*(1.d0+PI26*(CNU-1.d0)*(CNU-2.d0)/pow2(CHI) &
                  +.7d0*PI26*PI26*(CNU-1.d0)*(CNU-2.d0)*(CNU-3.d0)*(CNU-4.d0)/pow4(CHI))
               FDXXX=CHINU/pow2(CHI)*CNU*(CNU-1.d0)* &
                  (1.d0+PI26*(CNU-2.d0)*(CNU-3.d0)/pow2(CHI) &
                  +.7d0*PI26*PI26*(CNU-2.d0)*(CNU-3.d0)*(CNU-4.d0)*(CNU-5.d0)/pow4(CHI))
              if (J.eq.0) then
                 W=F
                 WDX=FDX
                 WDXX=FDXX
                 WDXXX=FDXXX
              elseif (J.eq.1) then
                 C=.25d0*TEMP
                 W=W+C*F ! Fermi-Dirac, expressed through Fermi
                 WDX=WDX+C*FDX
                 WDXX=WDXX+C*FDXX
                 WDT=F/4.d0
                 WDXT=FDX/4.d0
                 WDTT=0.d0
                 WDXXX=WDXXX+C*FDXXX
                 WDXXT=FDXX/4.d0
                 WDXTT=0.d0
              else
                 C=-C/(1d0*J)*(2d0*J-3d0)/4d0*TEMP
                 W=W+C*F
                 WDX=WDX+C*FDX
                 WDT=WDT+C*(1d0*J)/TEMP*F
                 WDXX=WDXX+C*FDXX
                 WDTT=WDTT+C*(1d0*J)*(1d0*(J-1))/pow2(TEMP)*F
                 WDXT=WDXT+C*(1d0*J)/TEMP*FDX
                 WDXXX=WDXXX+C*FDXXX
                 WDXTT=WDXTT+C*(1d0*J)*(1d0*(J-1))/pow2(TEMP)*FDX
                 WDXXT=WDXXT+C*(1d0*J)/TEMP*FDXX
              endif
            enddo ! next J
          if (K.eq.0) then
             W0=W
             W0DX=WDX
             W0DT=WDT
             W0DXX=WDXX
             W0DTT=WDTT
             W0DXT=WDXT
             W0XXX=WDXXX
             W0XTT=WDXTT
             W0XXT=WDXXT
          elseif (K.eq.1) then
             W1=W
             W1DX=WDX
             W1DT=WDT
             W1DXX=WDXX
             W1DTT=WDTT
             W1DXT=WDXT
          else
             W2=W
             W2DX=WDX
             W2DT=WDT
             W2DXX=WDXX
             W2DTT=WDTT
             W2DXT=WDXT
          endif
        enddo ! next K
!   ----------------------------------------------------------------   *
      else ! CHI > 14, CHI*TEMP > 0.1: general high-\chi expansion
         D=1.d0+CHI*TEMP/2.d0
         R=sqrt(CHI*D)
         RX=.5d0/CHI+.25d0*TEMP/D
         RDX=R*RX
         RDT=.25d0*pow2(CHI)/R
         RXX=-.5d0/pow2(CHI)-.125d0*pow2(TEMP/D)
         RDXX=RDX*RX+R*RXX
         RDTT=-.25d0*RDT*CHI/D
         RXT=.25d0/D-.125d0*CHI*TEMP/(D*D)
         RDXT=RDT*RX+R*RXT
         RXXX=1.d0/pow3(CHI)+.125d0*pow3(TEMP/D)
         RDXXX=RDXX*RX+2.d0*RDX*RXX+R*RXXX
         RXTT=-.25d0/(D*D)*CHI+.125d0*pow2(CHI)*TEMP/pow3(D)
         RDXTT=RDTT*RX+2.d0*RDT*RXT+R*RXTT
         RXXT=-RXT*TEMP/D
         RDXXT=RDXT*RX+RDX*RXT+RDT*RXX+R*RXXT
        do K=0,2
           DM=K+.5d0+(K+1.d0)*CHI*TEMP/2.d0
           AM(K)=pow(CHI,1d0*K)*DM/R
           FMX1=.5d0*(K+1.d0)*TEMP/DM
           FMX2=.25d0*TEMP/D
           FMX=(K-.5d0)/CHI+FMX1-FMX2
           AMDX(K)=AM(K)*FMX
           CKM=.5d0*(K+1.d0)/DM
           FMT1=CKM*CHI
           FMT2=.25d0*CHI/D
           FMT=FMT1-FMT2
           AMDT(K)=AM(K)*FMT
           FMXX=-(K-.5d0)/pow2(CHI)-FMX1*FMX1+2.d0*FMX2*FMX2
           AMDXX(K)=AMDX(K)*FMX+AM(K)*FMXX
           FMTT=2.d0*FMT2*FMT2-FMT1*FMT1
           AMDTT(K)=AMDT(K)*FMT+AM(K)*FMTT
           AMDXT(K)=AMDX(K)*FMT+AM(K)*(CKM*(1.d0-CKM*CHI*TEMP) &
              - .25d0/D+.125d0*CHI*TEMP/(D*D))
          if (K.eq.0) then
             FMXXX=(2d0*K-1d0)/pow3(CHI)+2.d0*pow3(FMX1)-8.d0*pow3(FMX2)
             AMDXXX=AMDXX(K)*FMX+2.d0*AMDX(K)*FMXX+AM(K)*FMXXX
             FMT1DX=CKM-TEMP*CHI*CKM*CKM
             FMT2DX=(.25d0-CHI*TEMP*.125d0/D)/D
             FMXT=FMT1DX-FMT2DX
             FMTTX=4.d0*FMT2*FMT2DX-2.d0*FMT1*FMT1DX
             AMDXTT=AMDXT(K)*FMT+AMDT(K)*FMXT+AMDX(K)*FMTT+AM(K)*FMTTX
             FMX1DT=CKM-CHI*TEMP*CKM*CKM
             FMX2DT=.25d0/D*(1.d0-.5d0*CHI*TEMP/D)
             FMXXT=4.d0*FMX2*FMX2DT-2.d0*FMX1*FMX1DT
             AMDXXT=AMDXT(K)*FMX+AMDX(K)*FMXT+AMDT(K)*FMXX+AM(K)*FMXXT
          endif
        enddo
           SQ2T=sqrt(2.d0*TEMP)
           SQ2T3=SQ2T*SQ2T*SQ2T
           A=1.d0+CHI*TEMP+SQ2T*R
           ADX=TEMP+SQ2T*RDX
           ADT=CHI+R/SQ2T+SQ2T*RDT
           ADXX=SQ2T*RDXX
           ADTT=-R/SQ2T3+2.d0/SQ2T*RDT+SQ2T*RDTT
           ADXT=1.d0+RDX/SQ2T+SQ2T*RDXT
           ADXTT=-RDX/SQ2T3+2.d0/SQ2T*RDXT+SQ2T*RDXTT
           ADXXT=RDXX/SQ2T+SQ2T*RDXXT
           XT1=CHI+1.d0/TEMP
           Aln=log(A)
           FJ0=.5d0*XT1*R-Aln/SQ2T3
           ASQ3=A*SQ2T3
           ASQ3DX=ADX*SQ2T3
           FJ0DX=.5d0*(R+XT1*RDX)-ADX/ASQ3
           FJ0DT=.5d0*(XT1*RDT-R/pow2(TEMP))-ADT/ASQ3 &
              + .75d0/(TEMP*TEMP*SQ2T)*Aln
           FJ0DXX=RDX+.5d0*XT1*RDXX+pow2(ADX/A)/SQ2T3-ADXX/ASQ3
           FJ0DTT=R/pow3(TEMP)-RDT/pow2(TEMP)+.5d0*XT1*RDTT &
              + 3.d0/(ASQ3*TEMP)*ADT &
              + pow2(ADT/A)/SQ2T3-ADTT/ASQ3-1.875d0/(TEMP*TEMP*TEMP*SQ2T)*Aln
           BXT=1.5d0/TEMP*ADX+ADX*ADT/A-ADXT
           BXXT=1.5d0/TEMP*ADXX+(ADXX*ADT+ADX*ADXT)/A &
              - pow2(ADX/A)*ADT-ADXXT
           FJ0DXT=.5d0*(RDT-RDX/pow2(TEMP)+XT1*RDXT)+BXT/ASQ3
           FJ0XXX=RDXX*1.5d0+.5d0*XT1*RDXXX &
              +(2.d0*ADX*(ADXX/A-pow2(ADX/A)) &
              - SQ2T*RDXXX+ADXX/ASQ3*ASQ3DX)/ASQ3
           FJ0XTT=RDX/pow3(TEMP)-RDXT/pow2(TEMP)+.5d0*(RDTT+XT1*RDXTT) &
              + 3.d0/TEMP*(ADXT-ADT/ASQ3*ASQ3DX)/ASQ3 &
              + (2.d0*ADT*(ADXT/A-ADT*ADX/(A*A)) &
              - ADXTT+ADTT*ASQ3DX/ASQ3)/ASQ3-1.875d0/(TEMP*TEMP*TEMP*SQ2T)*ADX/A
           FJ0XXT=.5d0*(RDXT-RDXX/pow2(TEMP)+RDXT+XT1*RDXXT) &
              +(BXXT-BXT*ASQ3DX/ASQ3)/ASQ3
         W0=FJ0+PI26*AM(0)
         W0DX=FJ0DX+PI26*AMDX(0)
         W0DT=FJ0DT+PI26*AMDT(0)
         W0DXX=FJ0DXX+PI26*AMDXX(0)
         W0DTT=FJ0DTT+PI26*AMDTT(0)
         W0DXT=FJ0DXT+PI26*AMDXT(0)
         W0XXX=FJ0XXX+PI26*AMDXXX
         W0XTT=FJ0XTT+PI26*AMDXTT
         W0XXT=FJ0XXT+PI26*AMDXXT
           FJ1=(R*R*R/1.5d0-FJ0)/TEMP
           FJ1DX=(2.d0*R*R*RDX-FJ0DX)/TEMP
           FJ1DT=(2.d0*R*R*RDT-FJ0DT-FJ1)/TEMP
           FJ1DXX=(4.d0*R*RDX*RDX+2.d0*R*R*RDXX-FJ0DXX)/TEMP
           FJ1DTT=(4.d0*R*RDT*RDT+2.d0*R*R*RDTT-FJ0DTT-2.d0*FJ1DT)/TEMP
           FJ1DXT=(4.d0*R*RDX*RDT+2.d0*R*R*RDXT-FJ0DXT-FJ1DX)/TEMP
         W1=FJ1+PI26*AM(1)
         W1DX=FJ1DX+PI26*AMDX(1)
         W1DT=FJ1DT+PI26*AMDT(1)
         W1DXX=FJ1DXX+PI26*AMDXX(1)
         W1DTT=FJ1DTT+PI26*AMDTT(1)
         W1DXT=FJ1DXT+PI26*AMDXT(1)
           FJ2=(.5d0*CHI*R*R*R-1.25d0*FJ1)/TEMP
           FJ2DX=(.5d0*R*R*R+1.5d0*CHI*R*R*RDX-1.25d0*FJ1DX)/TEMP
           FJ2DT=(1.5d0*CHI*R*R*RDT-1.25d0*FJ1DT-FJ2)/TEMP
           FJ2DXX=(3.d0*R*RDX*(R+CHI*RDX)+1.5d0*CHI*R*R*RDXX-1.25d0*FJ1DXX)/TEMP
           FJ2DTT=(3.d0*CHI*R*(RDT*RDT+.5d0*R*RDTT)-1.25d0*FJ1DTT-2.d0*FJ2DT)/TEMP
           FJ2DXT=(1.5d0*R*RDT*(R+2.d0*CHI*RDX)+1.5d0*CHI*R*R*RDXT-1.25d0*FJ1DXT-FJ2DX)/TEMP
         W2=FJ2+PI26*AM(2)
         W2DX=FJ2DX+PI26*AMDX(2)
         W2DT=FJ2DT+PI26*AMDT(2)
         W2DXX=FJ2DXX+PI26*AMDXX(2)
         W2DTT=FJ2DTT+PI26*AMDTT(2)
         W2DXT=FJ2DXT+PI26*AMDXT(2)
      endif
      return
      end subroutine BLIN9c

      subroutine CHEMFIT(DENS,TEMP,CHI)
!                                                       Version 07.06.07
! This is merely an interface to CHEMFIT7 for compatibility purposes.
! Input:  DENS - electron density [a.u.=6.7483346e24 cm^{-3}],
! TEMP - temperature [a.u.=2Ryd=3.1577e5 K]
! Output: CHI=\mu/TEMP, where \mu - electron chem.pot.w/o rest-energy
      type(auto_diff_real_2var_order1), intent(in) :: DENS, TEMP
      type(auto_diff_real_2var_order1), intent(out) :: CHI

      type(auto_diff_real_2var_order1) :: DENR,TEMR,CMU1,CMUDENR,CMUDT,CMUDTT
      
      DENR=DENS/2.5733806d6 ! n_e in rel.un.=\lambda_{Compton}^{-3}
      TEMR=TEMP/1.8778865d4 ! T in rel.un.=(mc^2/k)=5.93e9 K
      call CHEMFIT7(DENR,TEMR,CHI,CMU1,0,CMUDENR,CMUDT,CMUDTT)
      return
      end subroutine CHEMFIT

      subroutine CHEMFIT7(DENR,TEMR,CHI,CMU1,KDERIV, &
        CMUDENR,CMUDT,CMUDTT)
!                                                       Version 28.05.07
!                                                     corrected 08.08.11
!                                               cosmetic change 16.05.13
! Fit to the chemical potential of free electron gas described in:
!     G.Chabrier & A.Y.Potekhin, Phys.Rev.E 58, 4941 (1998)
! Stems from CHEMFIT v.10.10.96. The main difference - derivatives.
!  All quantities are by default in relativistic units
! Input:  DENR - electron density, TEMR - temperature
!         KDERIV=0 if the derivatives are not required
! Output: CHI=CMU1/TEMR, where CMU1 = \mu-1 - chem.pot.w/o rest-energy
!         CMUDENR = (d\mu/d n_e)_T
!         CMUDT = (d\mu/dT)_V
!         CMUDTT = (d^2\mu/dT^2)_V
! CMUDENR,CMUDT, and CMUDTT =0 on output, if KREDIV=0
      type(auto_diff_real_2var_order1), intent(in) :: DENR,TEMR
      integer, intent(in) :: KDERIV
      type(auto_diff_real_2var_order1), intent(out) :: CHI,CMU1,CMUDENR,CMUDT,CMUDTT

      type(auto_diff_real_2var_order1) :: PF0, TF, THETA, THETA32, Q2, T1, U3, THETAC, THETAG
      type(auto_diff_real_2var_order1) :: D3, Q3, Q1, SQT, G, H, CT, F, X, XDF, XDFF
      type(auto_diff_real_2var_order1) :: THETA52, CHIDY, CHIDYY
      type(auto_diff_real_2var_order1) :: Q1D, Q1DD, Q2D, Q2DD, U3D, D3D, D3DD, Q3D, Q3DD
      type(auto_diff_real_2var_order1) :: GDY, GDT, GDYY, GDTT, GDYT, GH
      type(auto_diff_real_2var_order1) :: HDY, HDT, HDYY, HDTT, HDYT
      type(auto_diff_real_2var_order1) :: CTY, CTT, CTDY, CTDT, CTDYY, CTDTT, CTDYT
      type(auto_diff_real_2var_order1) :: CHIDT, CHIDTT, CHIDYT
      
      real(dp), parameter :: PARA=1.612d0
      real(dp), parameter :: PARB=6.192d0
      real(dp), parameter :: PARC=0.0944d0
      real(dp), parameter :: PARF=5.535d0
      real(dp), parameter :: PARG=0.698d0
      
      PF0=pow(29.6088132d0*DENR,1d0/3d0) ! Classical Fermi momentum
      if (PF0.gt.1.d-4) then
         TF=sqrt(1.d0+PF0*PF0)-1.d0 ! Fermi temperature
      else
         TF=.5d0*PF0*PF0
      endif
      THETA=TEMR/TF
      THETA32=THETA*sqrt(THETA)
      Q2=12.d0+8.d0/THETA32
      T1=exp(-THETA) ! former ('96) 1/T
      U3=T1*T1+PARA
      THETAC=pow(THETA,PARC)
      THETAG=pow(THETA,PARG)
      D3=PARB*THETAC*T1*T1+PARF*THETAG
      Q3=1.365568127d0-U3/D3 ! 1.365...=2/\pi^{1/3}
      if (THETA.gt.1.d-5) then 
         Q1=1.5d0*T1/(1.d0-T1)
      else
         Q1=1.5d0/THETA
      endif
      SQT=sqrt(TEMR)
      G=(1.d0+Q2*TEMR*Q3+Q1*SQT)*TEMR
      H=(1.d0+.5d0*TEMR/THETA)*(1.d0+Q2*TEMR)
      CT=1.d0+G/H
      F=(2d0/3d0)/THETA32
      call FERINV7(F,1,X,XDF,XDFF)
      CHI=X & ! non-relativistic result
         - 1.5d0*log(CT) ! Relativistic fit
      CMU1=TEMR*CHI ! Fit to chemical potential w/o mc^2
      if (KDERIV.eq.0) then ! DISMISS DERIVATIVES
         CMUDENR=0.d0
         CMUDT=0.d0
         CMUDTT=0.d0
         return
      endif
! CALCULATE DERIVATIVES:
! 1: derivatives of CHI over THETA and T
! (a): Non-relativistic result:
      THETA52=THETA32*THETA
      CHIDY=-XDF/THETA52 ! d\chi/d\theta
      CHIDYY=(XDFF/pow4(THETA)-2.5d0*CHIDY)/THETA ! d^2\chi/d\theta^2
! (b): Relativistic corrections:
      if (THETA.gt.1.d-5) then 
         Q1D=-Q1/(1.d0-T1)
         Q1DD=-Q1D*(1.d0+T1)/(1.d0-T1)
      else
         Q1D=-1.5d0/pow2(THETA)
         Q1DD=-2.d0*Q1D/THETA
      endif
      Q2D=-12.d0/THETA52 ! d q_2 / d \theta
      Q2DD=30.d0/(THETA52*THETA) ! d^2 q_2 / d \theta^2
      U3D=-2.d0*T1*T1
      D3D=PARF*PARG*THETAG/THETA+PARB*T1*T1*THETAC*(PARC/THETA-2.d0)
      D3DD=PARF*PARG*(PARG-1.d0)*THETAG/pow2(THETA) &
         + PARB*T1*T1*THETAC*(PARC*(PARC-1.d0)/pow2(THETA)-4.d0*PARC/THETA+4.d0)
      Q3D=(D3D*U3/D3-U3D)/D3
      Q3DD=(2.d0*U3D+(2.d0*U3D*D3D+U3*D3DD)/D3-2.d0*U3*pow2(D3D/D3))/D3
      GDY=TEMR*(Q1D*SQT+(Q2D*Q3+Q2*Q3D)*TEMR) ! dG/d\theta
      GDT=1.d0+1.5d0*Q1*SQT+2.d0*Q2*Q3*TEMR
      GDYY=TEMR*(Q1DD*SQT+(Q2DD*Q3+2.d0*Q2D*Q3D+Q2*Q3DD)*TEMR)
      GDTT=.75d0*Q1/SQT+2.d0*Q2*Q3
      GDYT=1.5d0*Q1D*SQT+2.d0*(Q2D*Q3+Q2*Q3D)*TEMR
      HDY=(-.5d0/pow2(THETA)+Q2D+.5d0*(Q2D-Q2/THETA)/THETA*TEMR)*TEMR
      HDT=(.5d0+Q2*TEMR)/THETA+Q2
      HDYY=TEMR/pow3(THETA)+Q2DD*TEMR &
         + TEMR*TEMR*(.5d0*Q2DD-Q2D/THETA+Q2/pow2(THETA))/THETA
      HDTT=Q2/THETA
      HDYT=Q2D*(1.d0+TEMR/THETA)-(.5d0+Q2*TEMR)/pow2(THETA)
      CTY=GDY/G-HDY/H
      CTT=GDT/G-HDT/H
      GH=G/H
      CTDY=GH*CTY
      CTDT=GH*CTT
      CTDYY=CTDY*CTY+GH*(GDYY/G-pow2(GDY/G)-HDYY/H+pow2(HDY/H))
      CTDTT=CTDT*CTT+GH*(GDTT/G-pow2(GDT/G)-HDTT/H+pow2(HDT/H))
      CTDYT=CTDT*CTY+GH*(GDYT/G-GDY*GDT/pow2(G)-HDYT/H+HDY*HDT/pow2(H))
      CHIDY=CHIDY-1.5d0*CTDY/CT
      CHIDT=-1.5d0*CTDT/CT
      CHIDYY=CHIDYY+1.5d0*(pow2(CTDY/CT)-CTDYY/CT)
      CHIDTT=1.5d0*(pow2(CTDT/CT)-CTDTT/CT)
      CHIDYT=1.5d0*(CTDY*CTDT/pow2(CT)-CTDYT/CT)
      CMUDENR=-pow2(THETA*PF0)/(3.d0*DENR*(1.d0+TF))*CHIDY
      CMUDT=CHI+THETA*CHIDY+TEMR*CHIDT
      CMUDTT=2.d0*(CHIDY/TF+CHIDT+THETA*CHIDYT)+THETA/TF*CHIDYY+TEMR*CHIDTT
      return
      end subroutine CHEMFIT7
      
      end module pc_eos
