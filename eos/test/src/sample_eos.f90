! ***********************************************************************
!
!   Copyright (C) 2008-2019  Bill Paxton & The MESA Team
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
! ***********************************************************************

      program sample_eos
      use eos_def
      use eos_lib
      use chem_def
      use chem_lib
      use const_lib
      use math_lib

      implicit none

! this program shows how to setup and use the mesa eos in an interactive manner.


      real(dp) :: X, Z, Y, abar, zbar, z2bar, z53bar, ye
      integer, parameter :: species = 7
      integer, parameter :: h1=1, he4=2, c12=3, n14=4, o16=5, ne20=6, mg24=7
      integer, pointer, dimension(:) :: net_iso, chem_id
      real(dp) :: xa(species)


      call Sample
      
      contains
      
      subroutine Sample
         
         integer :: handle
         real(dp) :: Rho, T, Pgas, log10Rho
         real(dp) :: dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, d_dlnRho_const_T, d_dlnT_const_Rho
         real(dp), dimension(num_eos_basic_results) :: res, d_dlnd, d_dlnT, d_dabar, d_dzbar
         real(dp) :: d_dxa(num_eos_basic_results,species) 
         real(dp) :: xz, frac, dabar_dx(species), dzbar_dx(species), sumx,mass_correction, dmc_dx(species)
         integer :: ierr
         character (len=32) :: my_mesa_dir
         character (len=80) :: string


! explicitly set my_mesa_dir to your $MESA_DIR, or use a blank string, in which case your $MESA_DIR is automagically used

         my_mesa_dir = '../..'
!         my_mesa_dir = ''         


! initialize 
         ierr = 0

         call const_init(my_mesa_dir,ierr)     
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         
         call math_init()
         
         call chem_init('isotopes.data', ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)

         call Setup_eos(handle)
         allocate(net_iso(num_chem_isos), chem_id(species), stat=ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)

         call Init_Composition



! keep returning to here
100   xa(:)  = 0.0d0

! get the temperature, density and composition

      write(6,*)  
      write(6,*) 'give the temperature, density, and mass fractions (h1, he4, c12, n14, o16, ne20, mg24) =>'
      write(6,*) 'hit return for T = 1e9 K, Rho = 1e4 g/cc, x(c12) = 1 ; enter -1 to stop'
      write(6,*)  
      read(5,'(a)') string

! stop
      if (string(1:2) .eq. '-1') then
       call Shutdown_eos(handle)
       deallocate(net_iso, chem_id)
       if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
       stop 'normal termination'

! read the conditions
      else
       if (string(1:6) .ne. '      ') then
        read(string,*) T,Rho, xa(h1), xa(he4), xa(c12), xa(n14), xa(o16), xa(ne20), xa(mg24)

! or set some defaults
       else
        T = 1.0d9 ; Rho = 1.0d4 ; xa(c12) = 1.0d0
       end if
      end if


! get some composition variables
         call composition_info( &
               species, chem_id, xa, X, Y, xz, abar, zbar, z2bar, z53bar, ye, mass_correction, &
               sumx, dabar_dx, dzbar_dx, dmc_dx)
         

! call the density-temperature based eos

! composition derivatives in terms or abar and zbar
!         call eosDT_get_legacy( &
!               handle, Z, X, abar, zbar, &
!               species, chem_id, net_iso, xa, &
!               Rho, log10(Rho), T, log10(T), &
!               res, d_dlnd, d_dlnT, d_dabar, d_dzbar, ierr)


! composition derivatives by isotope

         call eosDT_get( &
               handle, species, chem_id, net_iso, xa, &
               Rho, log10(Rho), T, log10(T), &
               res, d_dlnd, d_dlnT, d_dxa, ierr)


! report the results

         call PrettyOut(Rho, T, abar, zbar, res, d_dlnd, d_dlnT)


! example calling the pressure-temperature based eos
!         Pgas = exp(res(i_lnPgas))
!         call eosPT_get( &
!               handle, Z, X, abar, zbar, &
!               species, chem_id, net_iso, xa, &
!               Pgas, log10(Pgas), T, log10(T), &
!               Rho, log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
!               res, d_dlnd, d_dlnT, d_dabar, d_dzbar, ierr)
      

! back for another one
      goto 100
         
      end subroutine Sample
      



      subroutine Setup_eos(handle)
         ! allocate and load the eos tables
         use eos_def
         use eos_lib
         integer, intent(out) :: handle

         integer :: ierr
         logical, parameter :: use_cache = .true.

         call eos_init(' ', use_cache, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)

         write(*,*) 'loading eos tables'
         
         handle = alloc_eos_handle(ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)

      end subroutine Setup_eos
      

      
      subroutine Shutdown_eos(handle)
         use eos_def
         use eos_lib
         integer, intent(in) :: handle
         call free_eos_handle(handle)
         call eos_shutdown
      end subroutine Shutdown_eos



      subroutine Init_Composition
         use chem_lib
       
         net_iso(:) = 0
         
         chem_id(h1) = ih1; net_iso(ih1) = h1
         chem_id(he4) = ihe4; net_iso(ihe4) = he4
         chem_id(c12) = ic12; net_iso(ic12) = c12
         chem_id(n14) = in14; net_iso(in14) = n14
         chem_id(o16) = io16; net_iso(io16) = o16
         chem_id(ne20) = ine20; net_iso(ine20) = ne20
         chem_id(mg24) = img24; net_iso(img24) = mg24
         
      end subroutine Init_Composition





      subroutine PrettyOut(Rho, T, abar, zbar, res, d_dlnd, d_dlnT)
         use eos_def
         use eos_lib
         use chem_lib

! declare the pass
         real(dp) :: Rho, T, abar, zbar
         real(dp), dimension(num_eos_basic_results) :: res, d_dlnd, d_dlnT, d_dabar, d_dzbar


! local variables
         real(dp) :: my_ptot, my_dptotdd, my_dptotdt, my_dptotddd, my_dptotddt, my_dptotdtt, &
                     my_etot, my_detotdd, my_detotdt, my_detotddd, my_detotddt, my_detotdtt, &
                     my_stot, my_dstotdd, my_dstotdt, my_dstotddd, my_dstotddt, my_dstotdtt, &
                     my_pgas, my_dpgasdd, my_dpgasdt, my_dpgasddd, my_dpgasddt, my_dpgasdtt, &
                     my_egas, my_degasdd, my_degasdt, my_degasddd, my_degasddt, my_degasdtt, &
                     my_sgas, my_dsgasdd, my_dsgasdt, my_dsgasddd, my_dsgasddt, my_dsgasdtt, &
                     my_prad, my_dpraddd, my_dpraddt, my_dpradddd, my_dpradddt, my_dpraddtt, &
                     my_erad, my_deraddd, my_deraddt, my_deradddd, my_deradddt, my_deraddtt, &
                     my_srad, my_dsraddd, my_dsraddt, my_dsradddd, my_dsradddt, my_dsraddtt, &
                     my_xni, my_dxnidd, my_dxnidt, my_xne, my_dxnedd, my_dxnedt, my_eta, my_detadd, my_detadt, &
                     my_cv, my_dcvdd, my_dcvdt, my_cp, my_dcpdd, my_dcpdt, & 
                     my_g1, my_dg1dd, my_dg1dt, my_g2, my_dg2dd, my_dg2dt, my_g3, my_dg3dd, my_dg3dt, & 
                     my_grad, my_dgraddd, my_dgraddt, my_chit, my_dchitdd, my_dchitdt, my_chir, my_dchirdd, my_dchirdt, &
                     my_cs, my_dcsdd, my_dcsdt, x, z, xx, yy, ww, dfk, dse, dpe, dsp

! popular format statements
01    format(1x,t2,a,t16,a,t31,a,t46,a,t59,a,t74,a,t89,a)
02    format(1x,t2,a,1p7e15.6)
03    format(1x,t2,a7,1pe14.6,t24,a7,1pe14.6,t46,a7,1pe14.6,t68,a7,1pe14.6)


! indices for the res, d_dlnd, and d_dlnt arrays are defined in $MESA_DIR/eos/public/eos_def.f90 

! radiation
         my_prad     = crad * one_third * T * T * T * T
         my_dpraddd  = 0.0d0
         my_dpraddt  = 4.0d0 * my_prad/T
         my_dpradddd = 0.0d0
         my_dpradddt = 0.0d0
         my_dpraddtt = 4.0d0 * crad * T * T 

         my_erad    = 3.0d0 * my_prad / Rho
         my_deraddd = -my_erad / Rho 
         my_deraddt = 3.0d0 * my_dpraddt / Rho
         my_deradddd = -2.0d0 * my_deraddd / Rho
         my_deradddt = -my_deraddt / Rho
         my_deraddtt = 3.0d0 * my_dpraddtt / Rho

         my_srad    = (my_prad/Rho + my_erad) / T
         my_dsraddd = (my_dpraddd/Rho - my_prad/(Rho*Rho) + my_deraddd) / T
         my_dsraddt = (my_dpraddt/Rho + my_deraddt - my_srad) / T
         my_dsradddd = ((my_dpradddd &
                    - (2.0d0*my_dpraddd - 2.0d0*my_prad/Rho)/Rho)/Rho &
                     + my_deradddd) / T
         my_dsradddt = ((my_dpradddt - my_dpraddt/Rho)/Rho &
                    + my_deradddt - my_dsraddd) / T
         my_dsraddtt = ((my_dpraddtt/Rho + my_deraddtt - my_dsraddt) - my_dsraddt)/ T


! gas and totals
! pressure
         my_pgas     = exp(res(i_lnPgas))
         my_dpgasdd  = my_pgas/Rho * d_dlnd(i_lnPgas)  
         my_dpgasdt  = my_pgas/T   * d_dlnt(i_lnPgas)  

         my_ptot    = my_pgas    + my_prad  
         my_dptotdd = my_dpgasdd + my_dpraddd
         my_dptotdt = my_dpgasdt + my_dpraddt

! energy
         my_etot     = exp(res(i_lnE))
         my_detotdd  = my_etot/Rho * d_dlnd(i_lnE)  
         my_detotdt  = my_etot/T   * d_dlnt(i_lnE)  

         my_egas    = my_etot    - my_erad
         my_degasdd = my_detotdd - my_deraddd
         my_degasdt = my_detotdt - my_deraddt

! entropy
         my_stot    = exp(res(i_lnS))
         my_dstotdd = my_stot/Rho * d_dlnd(i_lnS)  
         my_dstotdt = my_stot/T   * d_dlnt(i_lnS)  

         my_sgas    = my_stot    - my_srad
         my_dsgasdd = my_dstotdd - my_dsraddd
         my_dsgasdt = my_dstotdt - my_dsraddt


! number densities and electron degeneracy 
         my_xni    = avo * Rho/abar
         my_dxnidd = avo / abar
         my_dxnidt = 0.0d0

         my_xne    = abar * my_xni * exp(res(i_lnfree_e))
         my_dxnedd = abar * my_xni * exp(res(i_lnfree_e))/Rho * d_dlnd(i_lnfree_e) + abar * my_dxnidd * exp(res(i_lnfree_e))
         my_dxnedt = abar * my_xni * exp(res(i_lnfree_e))/T   * d_dlnt(i_lnfree_e) + abar * my_dxnidt * exp(res(i_lnfree_e))

         my_eta    = res(i_eta)
         my_detadd = d_dlnd(i_eta) / Rho
         my_detadt = d_dlnt(i_eta) / T

! total specific heats
         my_cv     = res(i_Cv)
         my_dcvdd  = d_dlnd(i_Cv) / Rho
         my_dcvdt  = d_dlnt(i_Cv) / T

         my_cp     = res(i_Cp)
         my_dcpdd  = d_dlnd(i_Cp) / Rho
         my_dcpdt  = d_dlnt(i_Cp) / T

! total gammas
         my_grad    = res(i_grad_ad)
         my_dgraddd = d_dlnd(i_grad_ad) / Rho
         my_dgraddt = d_dlnt(i_grad_ad) / T

         my_g1    = res(i_gamma1)
         my_dg1dd = d_dlnd(i_gamma1) / Rho
         my_dg1dt = d_dlnt(i_gamma1) / T

         my_g2    = 1.0d0/(1.0d0 - my_grad)
         my_dg2dd = my_g2*my_g2 * my_dgraddd
         my_dg2dt = my_g2*my_g2 * my_dgraddt

         my_g3    = res(i_gamma3)
         my_dg3dd = d_dlnd(i_gamma3) / Rho
         my_dg3dt = d_dlnt(i_gamma3) / T

! total adiabatic exponents
         my_chit    = res(i_chiT)
         my_dchitdd = d_dlnd(i_chiT) / Rho
         my_dchitdt = d_dlnt(i_chiT) / T

         my_chir    = res(i_chiRho)
         my_dchirdd = d_dlnd(i_chiRho) / Rho
         my_dchirdt = d_dlnt(i_chiRho) / T

! total sound speed
         x         = my_etot + clight*clight
         z         = 1.0d0 + x*Rho/my_ptot
         xx        = 1.0d0/z
         ww        = x*Rho/my_ptot
         dfk       = my_g1*xx
         my_cs     = clight*sqrt(dfk)
         yy        = 0.5d0*my_cs/dfk
         my_dcsdd = yy*((my_dg1dd - dfk)*xx * (my_detotdd*Rho - ww*my_dptotdd + x)/my_ptot)
         my_dcsdt = yy*((my_dg1dt - dfk)*xx * (my_detotdt*Rho - ww*my_dptotdt)/my_ptot)


! maxwell relations; each is at flaoting point if the consistency is perfect
         dse = T*my_dstotdt/my_detotdt - 1.0d0
         dpe = (my_detotdd*Rho*Rho + T*my_dptotdt)/my_ptot - 1.0d0
         dsp = -my_dstotdd*(Rho*Rho/my_dptotdt) - 1.0d0



! and finally some second derivatives

! pressure
         my_dpgasddd = my_pgas/Rho * my_dchirdd  
         my_dpgasddt = my_pgas/Rho * my_dchirdt 
         my_dpgasdtt = my_pgas/T * my_dchitdt 

         my_dptotddd = my_dpgasddd + my_dpradddd
         my_dptotddt = my_dpgasddt + my_dpradddt
         my_dptotdtt = my_dpgasdtt + my_dpraddtt

! energy
         my_detotddd = d_dlnd(i_dE_dRho)/Rho
         my_detotddt = d_dlnT(i_dE_dRho)/T
         my_detotdtt = my_dcvdt 

         my_degasddd = my_detotddd - my_deradddd
         my_degasddt = my_detotddt - my_deradddt
         my_degasdtt = my_detotdtt - my_deraddtt

! entropy
         my_dstotddd = d_dlnd(i_dS_dRho)/Rho
         my_dstotddt = d_dlnT(i_dS_dRho)/T
         my_dstotdtt = d_dlnT(i_dS_dT)/T

         my_dsgasddd = my_dstotddd - my_dsradddd
         my_dsgasddt = my_dstotddt - my_dsradddt
         my_dsgasdtt = my_dstotdtt - my_dsraddtt



! write 'em

         write(6,03) 'T     =',T,       'Rho   =',Rho,      'abar  =',abar,   'zbar  =',zbar
         write(6,03) 'h1    =',xa(h1),  'he4   =',xa(he4),  'c12   =',xa(c12), 'n14   =',xa(n14)
         write(6,03) 'o16   =',xa(o16), 'ne20  =',xa(ne20), 'mg24  =',xa(mg24)

         write(6,*)  ' '
         write(6,01)  'quantity','value','d/d(Rho)','d/d(T)','d^2/d(Rho)^2','d^2/d(Rho)d(T)','d^2/d(T)^2'

! pressure, energy, entropy
         write(6,02) 'p tot   =', my_ptot, my_dptotdd, my_dptotdt, my_dptotddd, my_dptotddt, my_dptotdtt
         write(6,02) 'p gas   =', my_pgas, my_dpgasdd, my_dpgasdt, my_dpgasddd, my_dpgasddt, my_dpgasdtt
         write(6,02) 'p rad   =', my_prad, my_dpraddd, my_dpraddt, my_dpradddd, my_dpradddt, my_dpraddtt

         write(6,*)
         write(6,02) 'e tot   =', my_etot, my_detotdd, my_detotdt, my_detotddd, my_detotddt, my_detotdtt
         write(6,02) 'e gas   =', my_egas, my_degasdd, my_degasdt, my_degasddd, my_degasddt, my_degasdtt
         write(6,02) 'e rad   =', my_erad, my_deraddd, my_deraddt, my_deradddd, my_deradddt, my_deraddtt

         write(6,*)
         write(6,02) 's tot   =', my_stot, my_dstotdd, my_dstotdt, my_dstotddd, my_dstotddt, my_dstotdtt
         write(6,02) 's gas   =', my_sgas, my_dsgasdd, my_dsgasdt, my_dsgasddd, my_dsgasddt, my_dsgasdtt
         write(6,02) 's rad   =', my_srad, my_dsraddd, my_dsraddt, my_dsradddd, my_dsradddt, my_dsraddtt


! ion and free electron matter number density, electron degeneracy parameter
         write(6,*)
         write(6,02) 'n_ion   =',my_xni, my_dxnidd, my_dxnidt
         write(6,02) 'n_ele   =',my_xne, my_dxnedd, my_dxnedt
         write(6,02) 'eta_e   =',my_eta, my_detadd, my_detadt

! specific heats, 3 gammas, sound speed, chit and chid for the gas and the total
         write(6,02) 'cv      =',my_cv, my_dcvdd, my_dcvdt
         write(6,02) 'cp      =',my_cp, my_dcpdd, my_dcpdt

         write(6,02) 'gamma_1 =',my_g1, my_dg1dd, my_dg1dt
         write(6,02) 'gamma_2 =',my_g2, my_dg2dd, my_dg2dt
         write(6,02) 'gamma_3 =',my_g3, my_dg3dd, my_dg3dt
         write(6,02) 'grad_ad =',my_grad, my_dgraddd, my_dgraddt

         write(6,02) 'chi_t   =',my_chit, my_dchitdd, my_dchitdt
         write(6,02) 'chi_d   =',my_chir, my_dchirdd, my_dchirdt

         write(6,02) 'c_sound =',my_cs, my_dcsdd, my_dcsdt

         write(6,*)
         write(6,03) 'dsp   =',dse,'dpe   =',dpe,'dsp   =',dsp

      end subroutine PrettyOut


      end program sample_eos

