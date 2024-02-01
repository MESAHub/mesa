!..here is the tabular helmholtz free energy eos:
!..
!..routine helmeos computes the pressure, energy and entropy via tables

      module helm
      use const_def, only: dp
      use math_lib

      implicit none
      
      
      logical, parameter :: dbg = .false.
      !logical, parameter :: dbg = .true.
      
      
      private :: dbg

      contains


      subroutine helmeos2( &
         T, logT, Rho, logRho, abar_in, zbar_in, &
         coulomb_temp_cut, coulomb_den_cut, helm_res, &
         clip_to_table_boundaries, include_radiation,  &
         include_elec_pos, off_table, ierr)
      use eos_def
      use const_def, only: pi, avo
      use utils_lib, only: is_bad
      implicit none
      real(dp), intent(in) :: T, logT, Rho, logRho
      real(dp), intent(in) :: abar_in, zbar_in
      real(dp), intent(in) :: coulomb_temp_cut, coulomb_den_cut
      real(dp), intent(inout) :: helm_res(num_helm_results)
      logical, intent(in) ::  &
         clip_to_table_boundaries, include_radiation,  &
         include_elec_pos
      logical, intent(out) :: off_table
      integer, intent(out) :: ierr
      
      logical :: skip_elec_pos

      include 'formats'
      
      ierr = 0
      off_table = .false.
      
      skip_elec_pos = .not. include_elec_pos
      call helmeos2aux( &
         T, logT, Rho, logRho, abar_in, zbar_in,   &
         coulomb_temp_cut, coulomb_den_cut, helm_res,  &
         clip_to_table_boundaries, include_radiation, (.not. include_elec_pos), off_table, ierr)
      if (off_table) return
     
      end subroutine helmeos2


      subroutine helmeos2aux( &
            temp_in, logtemp_in, den_in, logden_in, abar_in, zbar_in,  &
            coulomb_temp_cut, coulomb_den_cut, helm_res,  &
            clip_to_table_boundaries, include_radiation, must_skip_elec_pos, off_table, ierr)

      use eos_def
      use helm_polynomials
      use const_def, asol => crad
      use utils_lib, only: is_bad
      
      implicit none

      real(dp), intent(in) :: temp_in, logtemp_in, den_in, logden_in
      real(dp), intent(in) :: abar_in, zbar_in, coulomb_temp_cut, coulomb_den_cut
      real(dp), intent(inout) :: helm_res(num_helm_results)
      logical, intent(in) :: clip_to_table_boundaries, include_radiation, must_skip_elec_pos
      logical, intent(out) :: off_table
      integer, intent(out) :: ierr

      real(dp) :: h ! = planck_h
      type (Helm_Table), pointer :: ht
      

!..declare local variables
      include 'helm_declare_local_variables.dek'
      

!..given a temperature temp [K], density den [g/cm**3], and a composition 
!..characterized by abar and zbar, this routine returns most of the other 
!..thermodynamic quantities. of prime interest is the pressure [erg/cm**3], 
!..specific thermal energy [erg/gr], the entropy [erg/g/K], along with 
!..their derivatives with respect to temperature, density, abar, and zbar.
!..other quantites such the normalized chemical potential eta (plus its
!..derivatives), number density of electrons and positron pair (along 
!..with their derivatives), adiabatic indices, specific heats, and 
!..relativistically correct sound speed are also returned.
!..
!..this routine assumes planckian photons, an ideal gas of ions, 
!..and an electron-positron gas with an arbitrary degree of relativity
!..and degeneracy. interpolation in a table of the helmholtz free energy
!..is used to return the electron-positron thermodynamic quantities.
!..all other derivatives are analytic.
!..
!..references: cox & giuli chapter 24 ; timmes & swesty apj 1999

!..this routine assumes a call to subroutine read_helm_table has
!..been performed prior to calling this routine.


!..declare

      real(dp) abar, zbar, temp, logtemp, den, logden
      logical skip_elec_pos
      
!..for the interpolations
      integer          iat, jat
      real(dp) dth, dt2, dti, dt2i, dt3i, dd, dd2, ddi, dd2i, dd3i, &
                       xt, xd, mxt, mxd, fi(36), &
                       din, dindd, dinda, dindz, dindda, dinddz, dindaa, &
                       dindaz, dindzz, dinddaa, dinddaz, &
                       w0t, w1t, w2t, w0mt, w1mt, w2mt, &
                       w0d, w1d, w2d, w0md, w1md, w2md, &
                       dpepdd_in, dpepddd_in, dpepddt_in

      ! real(dp) psi0, dpsi0, ddpsi0, dddpsi0, &
      !                  psi1, dpsi1, ddpsi1, dddpsi1, &
      !                  psi2, dpsi2, ddpsi2, dddpsi2, &
      !                  h5

      ! real(dp) xpsi0, xdpsi0, xddpsi0, &
      !                  xpsi1, xdpsi1, xddpsi1, h3

      real(dp) si0t, si1t, si2t, si0mt, si1mt, si2mt, &
                       si0d, si1d, si2d, si0md, si1md, si2md, &
                       dsi0t, dsi1t, dsi2t, dsi0mt, dsi1mt, dsi2mt, &
                       dsi0d, dsi1d, dsi2d, dsi0md, dsi1md, dsi2md, &
                       ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, &
                       ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md, &
                       dddsi0t, dddsi1t, dddsi2t, &
                       dddsi0mt, dddsi1mt, dddsi2mt, &
                       dddsi0d, dddsi1d, dddsi2d, &
                       dddsi0md, dddsi1md, dddsi2md

      real(dp) free, df_d, df_t, df_dd, df_tt, df_dt, &
                       df_ttt, df_dtt, df_ddt, df_ddd

         ht => eos_ht

         ierr = 0
         off_table = .false.
         
         h = planck_h
         third  = 1.0d0/3.0d0
         sioncon = (2.0d0 * pi * amu * kerg)/(h*h)
         sifac  = 8.6322745944370191d-45
         kergavo = kerg * avo
         asoli3  = asol/3.0d0
         clight2 = clight*clight
         eostol = 1.0d-13
         fpmin  = 1.0d-14
         !..note: sifac = h**3/(2.0d0*pi*amu)**1.5d0
         forth   = 4.0d0/3.0d0
         fiveth  = 5.0d0/3.0d0
         teninth = 10.0d0/9.0d0
         esqu    = qe*qe
         forthpi = forth * pi
         
         abar = abar_in
         zbar = zbar_in
         temp = temp_in
         logtemp = logtemp_in
         den = den_in
         logden = logden_in

!..for very low T, convert all H to H2.  adjust abar and zbar accordingly.
         
         ! NOTE: table lookup uses din rather than den
         ytot1 = 1.0d0/abar
         ye    = ytot1 * zbar
         din     = ye*den
         
         skip_elec_pos = must_skip_elec_pos
         if (.not. skip_elec_pos) then ! see if need to set it true
            if (temp < ht% templo) then
               ierr = 1
               return
            end if
         
            if (din < ht% denlo) then
               ierr = 1
               return
            end if

            if (temp > ht% temphi) then            
               ierr = 1
               return
            end if
            
            if (din > ht% denhi) then
               ierr = 1
               return
            end if
         end if


         if (abar < 0d0) then
            ierr = 1
            return
         end if

!..very neutron rich compositions may need to be bounded, 
!..avoid that extrema for now in order to increase efficiency
!       ye    = max(1.0d-16, ye)


!..initialize local variables
       include 'helm_initialize_local_variables.dek'
       if (ierr /= 0) then
         if (dbg) write(*,*) 'failed in helm_initialize_local_variables'
         return
       end if

!..radiation section:
       if (include_radiation) then
         include 'helm_radiation.dek'
          if (ierr /= 0) then
            if (dbg) write(*,*) 'failed in helm_radiation'
            return
          end if
       end if

!..ion section:
       include 'helm_ideal_ions.dek'
       if (ierr /= 0) then
         if (dbg) write(*,*) 'failed in helm_ideal_ions'
         return
       end if

!..electron-positron section:
       if (.not. skip_elec_pos) then
         include 'helm_electron_positron.dek'
          if (ierr /= 0) then
            if (dbg) write(*,*) 'failed in helm_electron_positron'
            return
          end if
       end if

!..coulomb section:
       if ((ht% with_coulomb_corrections) .and. (.not. skip_elec_pos)) then
         include 'helm_coulomb2.dek'
          if (ierr /= 0) then
            if (dbg) write(*,*) 'failed in helm_coulomb2'
            return
          end if
       end if

!..sum the gas and total (gas + radiation) components
       include 'helm_sum_totals.dek'
       if (ierr /= 0) then
         if (dbg) write(*,*) 'failed in helm_sum_totals'
         return
       end if

       ! error out for very low pgas,egas,sgas (previously just set hard floor at this value)
       if(pgas < 1d-20 .or. egas < 1d-20 .or. sgas < 1d-20) then
          ierr = 1
          if (dbg) write(*,*) 'failed in helm, sums too small'
          return
       end if

!..compute the derivative quantities (cv, gamma1 ...etc)
       include 'helm_gammas.dek'
       if (ierr /= 0) then
         if (dbg) write(*,*) 'failed in helm_gammas'
         return
       end if

!..maxwell relations; each is zero if the consistency is perfect
!..if you don't need this, save three divides and comment this out
       x   = den * den
       dse = temp*dentrdt/denerdt - 1.0d0
       dpe = (denerdd*x + temp*dpresdt)/pres - 1.0d0
       dsp = -dentrdd*x/dpresdt - 1.0d0

!..store results
      include 'helm_store_results.dek'
      helm_res(h_crp:h_valid) = 0

!..debugging printout      
      if (.false.) then
         include 'helm_print_results.dek'
      end if

      end subroutine helmeos2aux

      subroutine show_h5(fi, ia, ja, w0t, w1t, w2t, w0mt, w1mt, w2mt, w0d, w1d, w2d, w0md, w1md, w2md)
            integer, intent(in) :: ia, ja
            real(dp), intent(in) :: fi(36), w0t, w1t, w2t, w0mt, w1mt, w2mt, w0d, w1d, w2d, w0md, w1md, w2md
            write(*,'(99a15)') 'w0t', 'w1t', 'w2t', 'w0mt', 'w1mt', 'w2mt', 'w0d', 'w1d', 'w2d', 'w0md', 'w1md', 'w2md'
            write(*,'(99e15.6)') w0t, w1t, w2t, w0mt, w1mt, w2mt, w0d, w1d, w2d, w0md, w1md, w2md
            write(*,'(a30,99e26.16)') 'fi(1)*w0d*w0t', fi(1)*w0d*w0t, fi(1),w0d,w0t
            write(*,'(a30,99e26.16)') 'fi(2)*w0md*w0t', fi(2)*w0md*w0t, fi(2),w0md,w0t
            write(*,'(a30,99e26.16)') 'fi(3)*w0d*w0mt', fi(3)*w0d*w0mt, fi(3),w0d,w0mt 
            write(*,'(a30,99e26.16)') 'fi(4)*w0md*w0mt', fi(4)*w0md*w0mt, fi(4),w0md,w0mt
            write(*,'(a30,99e26.16)') '1 + 2 + 3 + 4', fi(1)*w0d*w0t + fi(2)*w0md*w0t + fi(3)*w0d*w0mt + fi(4)*w0md*w0mt
            write(*,'(A)')
            write(*,'(a30,99e26.16)') 'fi(5)*w0d*w1t', fi(5)*w0d*w1t, fi(5),w0d,w1t  
            write(*,'(a30,99e26.16)') 'fi(6)*w0md*w1t', fi(6)*w0md*w1t, fi(6),w0md,w1t
            write(*,'(a30,99e26.16)') 'fi(7)*w0d*w1mt', fi(7)*w0d*w1mt, fi(7),w0d,w1mt 
            write(*,'(a30,99e26.16)') 'fi(8)*w0md*w1mt', fi(8)*w0md*w1mt, fi(8),w0md,w1mt
            write(*,'(a30,99e26.16)') 'fi(9)*w0d*w2t', fi(9)*w0d*w2t, fi(9),w0d,w2t
            write(*,'(a30,99e26.16)') 'fi(10)*w0md*w2t', fi(10)*w0md*w2t, fi(10),w0md,w2t
            write(*,'(a30,99e26.16)') 'fi(11)*w0d*w2mt', fi(11)*w0d*w2mt,  fi(11),w0d,w2mt 
            write(*,'(a30,99e26.16)') 'fi(12)*w0md*w2mt', fi(12)*w0md*w2mt, fi(12),w0md,w2mt 
            write(*,'(a30,99e26.16)') 'fi(13)*w1d*w0t', fi(13)*w1d*w0t, fi(13),w1d,w0t   
            write(*,'(a30,99e26.16)') 'fi(14)*w1md*w0t', fi(14)*w1md*w0t, fi(14),w1md,w0t
            write(*,'(a30,99e26.16)') 'fi(15)*w1d*w0mt', fi(15)*w1d*w0mt, fi(15),w1d,w0mt   
            write(*,'(a30,99e26.16)') 'fi(16)*w1md*w0mt', fi(16)*w1md*w0mt, fi(16),w1md,w0mt
            write(*,'(a30,99e26.16)') 'fi(17)*w2d*w0t', fi(17)*w2d*w0t,  fi(17),w2d,w0t  
            write(*,'(a30,99e26.16)') 'fi(18)*w2md*w0t', fi(18)*w2md*w0t, fi(18),w2md,w0t
            write(*,'(a30,99e26.16)') 'fi(19)*w2d*w0mt', fi(19)*w2d*w0mt, fi(19),w2d,w0mt  
            write(*,'(a30,99e26.16)') 'fi(20)*w2md*w0mt', fi(20)*w2md*w0mt, fi(20),w2md,w0mt 
            write(*,'(a30,99e26.16)') 'fi(21)*w1d*w1t', fi(21)*w1d*w1t, fi(21),w1d,w1t   
            write(*,'(a30,99e26.16)') 'fi(22)*w1md*w1t', fi(22)*w1md*w1t, fi(22),w1md,w1t
            write(*,'(a30,99e26.16)') 'fi(23)*w1d*w1mt', fi(23)*w1d*w1mt, fi(23),w1d,w1mt  
            write(*,'(a30,99e26.16)') 'fi(24)*w1md*w1mt', fi(24)*w1md*w1mt, fi(24),w1md,w1mt
            write(*,'(a30,99e26.16)') 'fi(25)*w2d*w1t', fi(25)*w2d*w1t, fi(25),w2d,w1t   
            write(*,'(a30,99e26.16)') 'fi(26)*w2md*w1t', fi(26)*w2md*w1t, fi(26),w2md,w1t
            write(*,'(a30,99e26.16)') 'fi(27)*w2d*w1mt', fi(27)*w2d*w1mt, fi(27),w2d,w1mt  
            write(*,'(a30,99e26.16)') 'fi(28)*w2md*w1mt', fi(28)*w2md*w1mt, fi(28),w2md,w1mt
            write(*,'(a30,99e26.16)') 'fi(29)*w1d*w2t', fi(29)*w1d*w2t, fi(29),w1d,w2t   
            write(*,'(a30,99e26.16)') 'fi(30)*w1md*w2t', fi(30)*w1md*w2t, fi(30),w1md,w2t
            write(*,'(a30,99e26.16)') 'fi(31)*w1d*w2mt', fi(31)*w1d*w2mt, fi(31),w1d,w2mt  
            write(*,'(a30,99e26.16)') 'fi(32)*w1md*w2mt', fi(32)*w1md*w2mt, fi(32),w1md,w2mt
            write(*,'(a30,99e26.16)') 'fi(33)*w2d*w2t', fi(33)*w2d*w2t, fi(33),w2d,w2t   
            write(*,'(a30,99e26.16)') 'fi(34)*w2md*w2t', fi(34)*w2md*w2t, fi(34),w2md,w2t
            write(*,'(a30,99e26.16)') 'fi(35)*w2d*w2mt', fi(35)*w2d*w2mt, fi(35),w2d,w2mt  
            write(*,'(a30,99e26.16)') 'fi(36)*w2md*w2mt', fi(36)*w2md*w2mt, fi(36),w2md,w2mt
      end subroutine show_h5

      
      end module
      


