! ***********************************************************************
!
!   Copyright (C) 2021  The MESA Team
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


      module phase_separation

      use star_private_def
      use const_def

      implicit none

      logical, parameter :: dbg = .false.

      ! offset to higher phase than 0.5 to avoid interference
      ! between phase separation mixing and latent heat for Skye.
      real(dp), parameter :: eos_phase_boundary = 0.9d0
      
      ! The saved_data type saves composition, energy,
      ! and crystallization data from star_info
      type saved_data
         real(dp), allocatable :: xa(:,:)
         real(dp), allocatable :: energy(:)
         real(dp) :: crystal_core_boundary_mass
      end type saved_data

      private
      public :: do_phase_separation

      contains

      subroutine do_phase_separation(s, dt, ierr)
         type (star_info), pointer :: s
         real(dp), intent(in) :: dt
         integer, intent(out) :: ierr

         if(s% phase_separation_option == 'CO') then
            call do_2component_phase_separation(s, dt, 'CO', ierr)
            ! -1 argument means iterate through all freezing material
         else if(s% phase_separation_option == 'ONe') then
            call do_2component_phase_separation(s, dt, 'ONe', ierr)
            ! -1 argument means iterate through all freezing material
         else if(s% phase_separation_option == 'distillation') then
            call do_distillation(s, dt, ierr)
            call smooth_eps_phase_sep(s,dt,ierr)
         else
            write(*,*) 'invalid phase_separation_option'
            stop
         end if
      end subroutine do_phase_separation

      subroutine do_2component_phase_separation(s, dt, components, ierr)
         use chem_def, only: chem_isos, ic12, io16, ine20
         use chem_lib, only: chem_get_iso_id
         type (star_info), pointer :: s
         real(dp), intent(in) :: dt
         character (len=*), intent(in) :: components
         integer, intent(out) :: ierr
         
         real(dp) :: XNe, XO, XC, pad
         integer :: k, k_bound, kstart, net_ic12, net_io16, net_ine20
         logical :: save_Skye_use_ion_offsets

         ! Set phase separation mixing mass negative at beginning of phase separation
         s% phase_sep_mixing_mass = -1d0
         s% eps_phase_separation(1:s%nz) = 0d0
         
         if(s% phase(s% nz) < eos_phase_boundary) then
            s% crystal_core_boundary_mass = 0d0
            return
         end if
            
         net_ic12 = s% net_iso(ic12)
         net_io16 = s% net_iso(io16)
         net_ine20 = s% net_iso(ine20)
         
         ! Find zone of phase transition from liquid to solid
         k_bound = -1
         do k = s%nz,1,-1
            if(s% phase(k-1) <= eos_phase_boundary .and. s% phase(k) > eos_phase_boundary) then
               k_bound = k
               exit
            end if
         end do
         
         XC = s% xa(net_ic12,k_bound)
         XO = s% xa(net_io16,k_bound)
         XNe = s% xa(net_ine20,k_bound)
         ! Check that we're still in C/O or O/Ne dominated material as appropriate,
         ! otherwise skip phase separation
         if(components == 'CO'.and. XO + XC < 0.9d0) return
         if(components == 'ONe'.and. XNe + XO < 0.8d0) return ! O/Ne mixtures tend to have more byproducts of burning mixed in
         
         ! If there is a phase transition, reset the composition at the boundary
         if(k_bound > 0) then

            ! core boundary needs to be padded by a minimal amount (less than a zone worth of mass)
            ! to account for loss of precision during remeshing.
            pad = s% min_dq * s% m(1) * 0.5d0
            do k = s%nz,1,-1
               if(s% m(k) > s% crystal_core_boundary_mass + pad) then
                  kstart = k
                  exit
               end if
            end do

            ! calculate energy associated with phase separation, ignoring the ionization
            ! energy term that Skye sometimes calculates
            save_Skye_use_ion_offsets = s% eos_rq% Skye_use_ion_offsets
            s% eos_rq% Skye_use_ion_offsets = .false.
            call update_model_(s,1,s%nz,.false.)
            do k=1,s% nz
               s% eps_phase_separation(k) = s% energy(k)
            end do
            
            ! loop runs outward starting at previous crystallization boundary
            do k = kstart,1,-1
               ! Start by checking if this material should be crystallizing
               if(s% phase(k) <= eos_phase_boundary) then
                  s% crystal_core_boundary_mass = s% m(k+1)
                  exit
               end if

               call move_one_zone(s,k,components)
               ! crystallized out to k now, liquid starts at k-1.
               ! now mix the liquid material outward until stably stratified
               call mix_outward(s, k-1, 0)
               
            end do

            call update_model_(s,1,s%nz,.false.)

            ! phase separation heating term for use by energy equation
            do k=1,s% nz
               s% eps_phase_separation(k) = (s% eps_phase_separation(k) - s% energy(k)) / dt
            end do
            s% eos_rq% Skye_use_ion_offsets = save_Skye_use_ion_offsets
            s% need_to_setvars = .true.
         end if

         ierr = 0
      end subroutine do_2component_phase_separation

      subroutine do_distillation(s, dt, ierr)
         use chem_def, only: chem_isos, ine22
         use chem_lib, only: chem_get_iso_id
         type (star_info), pointer :: s
         real(dp), intent(in) :: dt
         integer, intent(out) :: ierr
         
         type(saved_data) :: sd
         real(dp) :: L_distill, L_max, &
              retry_scale_factor, scale_factor_low, scale_factor_high, tol_low, tol_high
         integer :: k, iter, net_ine22
         logical :: save_Skye_use_ion_offsets, distilled

         call alloc_saved_data_(s, sd)
         s% eps_phase_separation(1:s%nz) = 0d0

         net_ine22 = s% net_iso(ine22)
         if(net_ine22 < 0) then
            write(*,*) 'distillation requires ne22 be included in net'
            stop
         end if

         ! Set phase separation mixing mass negative at beginning of distillation
         s% phase_sep_mixing_mass = -1d0
         
         ! calculate energy associated with phase separation, ignoring the ionization
         ! energy term that Skye sometimes calculates
         save_Skye_use_ion_offsets = s% eos_rq% Skye_use_ion_offsets
         s% eos_rq% Skye_use_ion_offsets = .false.
         call update_model_(s,1,s%nz,.false.)
         call save_model_(s,sd) ! sets saved data sd to starting composition, energy, and crystal core location

         call distill_loop(s,sd,1d0,distilled)
         
         ! checking luminosity
         L_max = s% L_phot * Lsun
         call calc_L_distill(s,sd,L_distill)

         if(dbg) print *, "distill iter, scale_factor, Lum_ratio", 0, 1d0, L_distill/L_max

         ! try to converge amount of distillation toward not producing too much luminosity
         if(L_distill > L_max .and. distilled) then
            ! reset model and try again with scale_factor for changes
            retry_scale_factor = L_max/L_distill
            scale_factor_low = 0d0
            scale_factor_high = 1d0

            tol_low = 0.9d0
            tol_high = 1.0d0 ! initially try for something that strictly does not exceed WD luminosity
            
            iter = 0
            do while((L_distill < tol_low*L_max .or. L_distill > tol_high*L_max) .and. iter < 12)
               iter = iter + 1
               
               ! reset model
               s% phase_sep_mixing_mass = -1d0
               call restore_composition_(s,sd)
               call update_model_(s,1,s%nz,.true.)

               call distill_loop(s,sd,retry_scale_factor,distilled)

               call calc_L_distill(s,sd,L_distill)
               if(dbg) print *, "retry iter, scale_factor, Lum_ratio", iter, retry_scale_factor, L_distill/L_max

               ! simple bisection seems to work
               if(L_distill < tol_low*L_max) then
                  ! need to try a larger scale factor
                  scale_factor_low = retry_scale_factor
                  retry_scale_factor = 0.5d0*(retry_scale_factor + scale_factor_high)
               else
                  scale_factor_high = retry_scale_factor
                  retry_scale_factor = 0.5d0*(retry_scale_factor + scale_factor_low)
               end if

               ! relax tolerances after 4 and 8 iters, prefer L_distill < L
               if(iter > 8) then
                  tol_low = 0.3d0
                  tol_high = 1.1d0
               else if(iter > 4) then
                  tol_low = 0.7d0
                  tol_high = 1.02d0
               end if
            end do
            if(dbg) print *, "converged after iter, L_dist, L_max, ratio", iter, L_distill/Lsun, L_max/Lsun, L_distill/L_max
         end if
         
         call update_model_(s,1,s%nz,.false.)
         
         ! phase separation heating term for use by energy equation
         do k=1,s% nz
            s% eps_phase_separation(k) = (sd% energy(k) - s% energy(k)) / dt
         end do
         s% eos_rq% Skye_use_ion_offsets = save_Skye_use_ion_offsets
         s% need_to_setvars = .true.
         
         ierr = 0
      end subroutine do_distillation

      subroutine distill_loop(s,sd,scale_factor,distilling)
         use chem_def, only: chem_isos, ic12, io16, ine20, ine22
         use chem_lib, only: chem_get_iso_id
         
         type (star_info), pointer    :: s
         type(saved_data), intent(in) :: sd
         real(dp), intent(in)         :: scale_factor
         logical, intent(out)         :: distilling
         
         real(dp) :: Gamma_melt, GammaC_melt, GammaC_melt_CO, GammaC, L_distill, &
              XNe, XNe_out, XO, XC, pad, distill_final_XNe, XNe_crit, xne_num, xo_num
         integer :: k, kstart, net_ic12, net_io16, net_ine20, net_ine22

         distill_final_XNe = 0.3143d0
         
         net_ic12 = s% net_iso(ic12)
         net_io16 = s% net_iso(io16)
         net_ine20 = s% net_iso(ine20)
         net_ine22 = s% net_iso(ine22)
         
         pad = s% min_dq * s% m(1) * 0.5d0
         do k = s%nz,1,-1
            if(s% m(k) > s% crystal_core_boundary_mass + pad) then
               kstart = k
               exit
            end if
         end do

         distilling = .false.
         do k = kstart, 1, -1
            XC = s% xa(net_ic12,k)
            XO = s% xa(net_io16,k)
            XNe = s% xa(net_ine20,k) + s% xa(net_ine22,k)
            XNe_out = s% xa(net_ine20,k-1) + s% xa(net_ine22,k-1)
            if(XC + XO + XNe < 0.9d0) then
               exit
            end if

            ! some number fractions for expressions from Simon
            xne_num = XNe*s% abar(k)/22d0 ! assumes 22Ne is dominant isotope for this conversion
            xo_num = XO*s% abar(k)/16d0
            
            ! Check whether we are in a regime where distillation should occur
            GammaC = s% gam(k) * pow(6d0,5d0/3d0) / s% z53bar(k) ! <Gamma> * 6^(5/3) / <Z^(5/3)>
            Gamma_melt = blouin_Gamma_melt_CO(xo_num)

            ! no correction for presence of Ne yet
            ! conversion assumes binary C/O mixture, so xc = 1 - xo
            GammaC_melt_CO = Gamma_melt * pow(6d0,5d0/3d0) / ((1-xo_num)*pow(6d0,5d0/3d0) + xo_num*pow(8d0,5d0/3d0))
            ! this way is incorrect:
            ! GammaC_melt_CO = Gamma_melt * pow(6d0,5d0/3d0) / s% z53bar(k) ! <Gamma>_m * 6^(5/3) / <Z^(5/3)>

            ! now include correction for Ne
            GammaC_melt = GammaC_melt_CO & ! <Gamma>_m * 6^(5/3) / <Z^(5/3)> (assuming pure C/O)
                 + 1096.69d0*xne_num*xo_num & ! corrections to phase diagram accounting for presence of Ne
                 - 3410.33d0*xne_num*xo_num*xo_num & ! fits from Simon Blouin (priv comm)
                 + 2408.44d0*xne_num*xo_num*xo_num*xo_num

            XNe_crit = blouin_XNe_crit(GammaC_melt)
            
            if(XNe > max(XNe_crit,1d-7) .and. &
                 XNe < distill_final_XNe .and. &
                 GammaC > GammaC_melt) then
            
               ! should be distilling in this zone
               distilling = .true. ! so we know not to do 2 component separation later
               
               ! must distill to xNe = 0.2 (XNe = 0.3143) and xC = 0.8 (no O).
               ! Once distillation starts in a zone, it stays liquid and continues
               ! distilling until reaching xNe = 0.2.
               call distill_at_boundary(s,k,XNe_crit,distill_final_XNe,GammaC,GammaC_melt_CO,scale_factor)

               ! mix from zone k-1 outward
               call mix_outward(s, k-1, 0) ! TODO: add inlist option for how many cells to include in last argument
               XNe = s% xa(net_ine20,k) + s% xa(net_ine22,k)
               XNe_out = s% xa(net_ine20,k-1) + s% xa(net_ine22,k-1)
               
               if(XNe >= distill_final_XNe .or. XNe_out < 1d-7) then
                  if( k == s% nz .or. s% crystal_core_boundary_mass + pad > s% m(min(k+1,s%nz)) ) then
                     ! done distilling this zone and everything inward from it, mark crystallized
                     s% crystal_core_boundary_mass = s% m(k)
                  end if
               end if
               
               ! break out of loop once phase sep luminosity approaches that of the star
               call calc_L_distill(s,sd,L_distill)
               if(L_distill >  0.9d0 * s% L_phot * Lsun .and. scale_factor < 1d0) exit
            else if (GammaC > GammaC_melt .and. (.not. distilling)) then
               ! also check that we're done with everything inward from this point
               if( k == s% nz .or. s% crystal_core_boundary_mass + pad > s% m(min(k+1,s%nz)) ) then
                  ! zone won't distill, but is ready to phase separate C/O
                  call move_one_zone_for_distill(s,k)
                  ! crystallized out to k now, liquid starts at k-1.
                  ! now mix the liquid material outward until stably stratified
                  call mix_outward(s, k-1, 0)
                  s% crystal_core_boundary_mass = s% m(k)

                  ! break out of loop once phase sep luminosity approaches that of the star
                  call calc_L_distill(s,sd,L_distill)
                  if(L_distill >  0.9d0 * s% L_phot * Lsun) exit
               end if
            else if (GammaC < GammaC_melt) then
               exit ! zones outward from here are not cool enough to freeze yet
            end if
         end do
      end subroutine distill_loop
      
      subroutine distill_at_boundary(s,k,XNe_crit,distill_final_XNe,GammaC,GammaC_melt,scale_factor)
        use chem_def, only: chem_isos, ic12, io16, ine20, ine22
        use chem_lib, only: chem_get_iso_id
        type(star_info), pointer :: s
        integer, intent(in) :: k
        real(dp), intent(in) :: XNe_crit, distill_final_XNe, GammaC, GammaC_melt, scale_factor
        
        real(dp) :: XC, XO, XNe, XC_out, XO_out, XNe_out, Xout_sum, target_XNe
        real(dp) :: Delta_XC, Delta_XO, Delta_XNe, max_Delta_XC, max_Delta_XO, max_Delta_XNe
        real(dp) :: dXC, dXO, dXNe_tot, dXNe20, dXNe22, dq_ratio, scale
        real(dp) :: delta_binding_energy, Lmax
        integer :: net_ic12, net_io16, net_ine20, net_ine22

        dq_ratio = s% dq(k) / s% dq(k-1)
        
        net_ic12 = s% net_iso(ic12)
        net_io16 = s% net_iso(io16)
        net_ine20 = s% net_iso(ine20)
        net_ine22 = s% net_iso(ine22)
        
        XC = s% xa(net_ic12,k)
        XO = s% xa(net_io16,k)
        XNe = s% xa(net_ine20,k) + s% xa(net_ine22,k)
        
        XC_out = s% xa(net_ic12,k-1)
        XO_out = s% xa(net_io16,k-1)
        XNe_out = s% xa(net_ine20,k-1) + s% xa(net_ine22,k-1)

        ! distillation not complete until GammaC = 208,
        ! so draw a line between start of distillation and GammaC = 208, XNe = distill_final_XNe,
        ! only go as far along that line as current GammaC allows.
        target_XNe = XNe_crit + (distill_final_XNe - XNe_crit + 1d-5)*(GammaC - GammaC_melt)/(208d0 - GammaC_melt)
        if(XNe >= target_XNe) then
           ! not cool enough to continue distillation yet
           return
        end if

        ! Net effect of distillation is that crystals enriched in oxygen float upward.
        ! Need to limit toward xNe = 0.2, xC = 0.8. Start by pushing O outward in exchange
        ! for C/Ne mixture until O is depleted. Then exchange C/Ne until reaching critical
        ! Ne value.

        ! start by calculating difference between zone composition (XC,XO,XNe)
        ! and target compostion ~(1-distill_final_XNe,0,distill_final_XNe) (not accounting for trace impurities)
        Delta_XO = -XO ! get rid of all O
        Delta_XNe = distill_final_XNe - XNe + 1d-5 ! reach target Ne fraction (with some pad to alleviate roundoff error)
        Delta_XC = XO - Delta_XNe ! All O that doesn't become Ne must become C
        
        ! Which element will limit the size of composition step.
        ! C and Ne need to increase, so check how much is available in next zone out
        max_Delta_XC = min(Delta_XC,XC_out/dq_ratio)
        max_Delta_XNe = min(Delta_XNe,XNe_out/dq_ratio,target_XNe-XNe)

        ! O needs to go to zero by getting pushed into next zone out,
        ! so check how much next zone out can accept
        max_Delta_XO = min(XO,(1d0-XO_out)/dq_ratio)

        ! check which gives the smallest scale factor, and possibly scale downward further by input scale_factor
        scale = scale_factor*min(max_Delta_XC/Delta_XC, max_Delta_XO/XO, max_Delta_XNe/Delta_XNe)

        ! Now rescale all changes to the same scale factor
        dXC = scale*Delta_XC
        dXO = scale*Delta_XO
        dXNe_tot = scale*Delta_XNe
        dXNe20 = dXNe_tot*s% xa(net_ine20,k-1)/XNe_out ! proportional to composition of outer zone because that's what distills inward
        dXNe22 = dXNe_tot*s% xa(net_ine22,k-1)/XNe_out

        ! for debugging
        if(dbg) then
           print *, "k, nz, dq_ratio", k, s%nz, dq_ratio
           print *, "before distill at boundary; XC, XO, XNe, XC_out, XO_out, XNe_out", &
                XC, XO, XNe, XC_out, XO_out, XNe_out
           
           print *, "max_Delta_XC, max_Delta_XO, max_Delta_XNe, dXC, dXO, dXNe", &
                max_Delta_XC, max_Delta_XO, max_Delta_XNe, dXC, dXO, dXNe_tot
        end if
        
        s% xa(net_ic12,k) = s% xa(net_ic12,k) + dXC
        s% xa(net_io16,k) = s% xa(net_io16,k) + dXO ! dXO should be negative, so we're pushing O out of this zone
        s% xa(net_ine20,k) = s% xa(net_ine20,k) + dXNe20
        s% xa(net_ine22,k) = s% xa(net_ine22,k) + dXNe22
        
        ! use dq_ratio to conserve total mass of each element,
        ! accounting for the fact that zone k-1 has different mass than zone k.
        s% xa(net_ic12,k-1) = s% xa(net_ic12,k-1) - dXC*dq_ratio
        s% xa(net_io16,k-1) = s% xa(net_io16,k-1) - dXO*dq_ratio
        s% xa(net_ine20,k-1) = s% xa(net_ine20,k-1) - dXNe20*dq_ratio
        s% xa(net_ine22,k-1) = s% xa(net_ine22,k-1) - dXNe22*dq_ratio
        
        ! for debugging
        if(dbg) then
           XC = s% xa(net_ic12,k)
           XO = s% xa(net_io16,k)
           XNe = s% xa(net_ine20,k) + s% xa(net_ine22,k)
           
           XC_out = s% xa(net_ic12,k-1)
           XO_out = s% xa(net_io16,k-1)
           XNe_out = s% xa(net_ine20,k-1) + s% xa(net_ine22,k-1)
           print *, "after distill at boundary;  XC, XO, XNe, XC_out, XO_out, XNe_out", &
                XC, XO, XNe, XC_out, XO_out, XNe_out
        end if
        
        call update_model_(s,k-1,k,.false.)
        
      end subroutine distill_at_boundary

      subroutine move_one_zone(s,k,components)
        use chem_def, only: chem_isos, ic12, io16, ine20
        use chem_lib, only: chem_get_iso_id
        type(star_info), pointer :: s
        integer, intent(in) :: k
        character (len=*), intent(in) :: components
        
        real(dp) :: XC, XO, XNe, XC1, XO1, XNe1, dXO, dXNe, Xfac
        integer :: net_ic12, net_io16, net_ine20
        
        net_ic12 = s% net_iso(ic12)
        net_io16 = s% net_iso(io16)
        net_ine20 = s% net_iso(ine20)
        
        if(components == 'CO') then
           XO = s% xa(net_io16,k)
           XC = s% xa(net_ic12,k)
        
           ! Call Blouin phase diagram.
           ! Need to rescale temporarily because phase diagram assumes XO + XC = 1
           Xfac = XO + XC
           XO = XO/Xfac
           XC = XC/Xfac
           
           dXO = blouin_delta_xo(XO)
           
           s% xa(net_io16,k) = Xfac*(XO + dXO)
           s% xa(net_ic12,k) = Xfac*(XC - dXO)
           
           ! Redistribute change in C,O into zone k-1,
           ! conserving total mass of C,O
           XC1 = s% xa(net_ic12,k-1)
           XO1 = s% xa(net_io16,k-1)
           s% xa(net_ic12,k-1) = XC1 + Xfac*dXO * s% dq(k) / s% dq(k-1)
           s% xa(net_io16,k-1) = XO1 - Xfac*dXO * s% dq(k) / s% dq(k-1)
        else if(components == 'ONe') then
           XNe = s% xa(net_ine20,k)
           XO = s% xa(net_io16,k)
        
           ! Call Blouin phase diagram.
           ! Need to rescale temporarily because phase diagram assumes XO + XNe = 1
           Xfac = XO + XNe
           XO = XO/Xfac
           XNe = XNe/Xfac
           
           dXNe = blouin_delta_xne(XNe)
           
           s% xa(net_ine20,k) = Xfac*(XNe + dXNe)
           s% xa(net_io16,k) = Xfac*(XO - dXNe)
           
           ! Redistribute change in Ne,O into zone k-1,
           ! conserving total mass of Ne,O
           XO1 = s% xa(net_io16,k-1)
           XNe1 = s% xa(net_ine20,k-1)
           s% xa(net_io16,k-1) = XO1 + Xfac*dXNe * s% dq(k) / s% dq(k-1)
           s% xa(net_ine20,k-1) = XNe1 - Xfac*dXNe * s% dq(k) / s% dq(k-1)
        else
           write(*,*) 'invalid components option in phase separation'
           stop
        end if

        call update_model_(s,k-1,s%nz,.true.)
        
      end subroutine move_one_zone

      subroutine move_one_zone_for_distill(s,k)
        use chem_def, only: chem_isos, ic12, io16, ine20, ine22
        use chem_lib, only: chem_get_iso_id
        type(star_info), pointer :: s
        integer, intent(in) :: k
        
        real(dp) :: XC, XO, XNe, dXC, dXO, dXNe, dXNe20, dXNe22, Xfac, Xnew
        real(dp) :: xo_num, xne_num, dxo_num, dxne_num, abar_new
        integer :: net_ic12, net_io16, net_ine20, net_ine22
        
        net_ic12 = s% net_iso(ic12)
        net_io16 = s% net_iso(io16)
        net_ine20 = s% net_iso(ine20)
        net_ine22 = s% net_iso(ine22)
        
        XC = s% xa(net_ic12,k)
        XO = s% xa(net_io16,k)
        XNe = s% xa(net_ine20,k) + s% xa(net_ine22,k)

        ! rescale to C + O + Ne = 1
        Xfac = XO + XC + XNe
        XC = XC/Xfac
        XO = XO/Xfac
        XNe = XNe/Xfac
           
        xo_num = XO*s% abar(k)/16d0
        xne_num = XNe*s% abar(k)/22d0 ! assumes 22Ne is dominant isotope for this conversion

        ! change in number fractions
        dxo_num = blouin_delta_xo_3component(xo_num,xne_num)
        dxne_num = blouin_delta_xne_3component(xo_num,xne_num) ! this is <=0 by construction

        ! convert changes to mass fractions
        abar_new = &
             16d0*(xo_num + dxo_num) + &
             22d0*(xne_num + dxne_num) + &
             12d0*(1d0 - xo_num - dxo_num - xne_num - dxne_num)
        Xnew = 16d0*(xo_num + dxo_num)/abar_new ! placeholder for oxygen
        dXO = Xnew - XO
        Xnew = 22d0*(xne_num + dxne_num)/abar_new ! placeholder for neon
        dXNe = XNew - XNe
        dXC = -dXO -dXNe

        ! allot change in isotopes proportionally
        dXNe20 = dXNe*s% xa(net_ine20,k)/(XNe*Xfac)
        dXNe22 = dXNe*s% xa(net_ine22,k)/(XNe*Xfac)

        ! adjust mass fractions in the MESA model
        s% xa(net_ic12,k) = Xfac*(XC + dXC)
        s% xa(net_io16,k) = Xfac*(XO + dXO)
        s% xa(net_ine20,k) = Xfac*(s% xa(net_ine20,k) + dXNe20)
        s% xa(net_ine22,k) = Xfac*(s% xa(net_ine22,k) + dXNe22)
        
        ! Redistribute change in C,O,Ne into zone k-1, conserving total mass of C,O,Ne
        s% xa(net_ic12,k-1) = s% xa(net_ic12,k-1) - Xfac*dXC * s% dq(k) / s% dq(k-1)
        s% xa(net_io16,k-1) = s% xa(net_io16,k-1) - Xfac*dXO * s% dq(k) / s% dq(k-1)
        s% xa(net_ine20,k-1) = s% xa(net_ine20,k-1) - Xfac*dXNe20 * s% dq(k) / s% dq(k-1)
        s% xa(net_ine22,k-1) = s% xa(net_ine22,k-1) - Xfac*dXNe22 * s% dq(k) / s% dq(k-1)

        if(dbg) then
           print *, "move one zone for distill"
           print *, "abundances in zone", k, s% xa(net_ic12,k), s% xa(net_io16,k), s% xa(net_ine20,k), s% xa(net_ine22,k)
           print *, "abundances in zone", k-1, s% xa(net_ic12,k-1), s% xa(net_io16,k-1), s% xa(net_ine20,k-1), s% xa(net_ine22,k-1)
           print *, "dXC, dXO, dXNe", dXC, dXO, dXNe
           print *, "dq ratio (k-1/k)", s%dq(k-1)/s%dq(k)
        end if
        
        call update_model_(s,k-1,s%nz,.true.)
        
      end subroutine move_one_zone_for_distill

      
      ! mix composition outward until reaching stable composition profile
      subroutine mix_outward(s,kbot,min_mix_zones)
        type(star_info), pointer :: s
        integer, intent(in)      :: kbot, min_mix_zones
        
        real(dp) :: avg_xa(s%species)
        real(dp) :: mass, B_term, grada, gradr
        integer :: k, l, ktop
        logical :: use_brunt
        
        use_brunt = s% phase_separation_mixing_use_brunt

        do k=kbot-min_mix_zones,1,-1
           ktop = k

           if (s% m(ktop) > s% phase_sep_mixing_mass) then
              s% phase_sep_mixing_mass = s% m(ktop)
           end if

           mass = SUM(s%dm(ktop:kbot))
           do l = 1, s%species
              avg_xa(l) = SUM(s%dm(ktop:kbot)*s%xa(l,ktop:kbot))/mass
           end do

           ! some potential safeguards from conv_premix
           ! avg_xa = MAX(MIN(avg_xa, 1._dp), 0._dp)
           ! avg_xa = avg_xa/SUM(avg_xa)

           do l = 1, s%species
              s%xa(l,ktop:kbot) = avg_xa(l)
           end do

           ! updates, eos, opacities, mu, etc now that abundances have changed,
           ! but only in the cells near the boundary where we need to check here.
           ! Will call full update over mixed region after exiting loop.
           call update_model_(s, ktop-1, ktop+1, use_brunt)

           if(use_brunt) then
              B_term = s% unsmoothed_brunt_B(ktop)
              grada = s% grada_face(ktop)
              gradr = s% gradr(ktop)
              if(B_term + grada - gradr > 0d0) then
                 ! stable against further mixing, so exit loop
                 exit
              end if
           else ! simpler calculation based on mu gradient
              if(s% mu(ktop) >= s% mu(ktop-1)) then
                 ! stable against further mixing, so exit loop
                 exit
              end if
           end if

        end do

        ! Call a final update over all mixed cells now.
        call update_model_(s, ktop, kbot+1, .true.)
       
      end subroutine mix_outward

      real(dp) function blouin_delta_xo(Xin)
        real(dp), intent(in) :: Xin ! mass fraction
        real(dp) :: Xnew ! mass fraction
        real(dp) :: xo, dxo ! number fractions
        real(dp) :: a0, a1, a2, a3, a4, a5

        ! Convert input mass fraction to number fraction, assuming C/O mixture
        xo = (Xin/16d0)/(Xin/16d0 + (1d0 - Xin)/12d0)
        
        a0 = 0d0
        a1 = -0.311540d0
        a2 = 2.114743d0
        a3 = -1.661095d0
        a4 = -1.406005d0
        a5 = 1.263897d0

        dxo = &
             a0 + &
             a1*xo + &
             a2*xo*xo + &
             a3*xo*xo*xo + &
             a4*xo*xo*xo*xo + &
             a5*xo*xo*xo*xo*xo

        xo = xo + dxo

        ! Convert back to mass fraction
        Xnew = 16d0*xo/(16d0*xo + 12d0*(1d0-xo))
        
        blouin_delta_xo = Xnew - Xin
      end function blouin_delta_xo

      real(dp) function blouin_Gamma_melt_CO(xo)
        real(dp), intent(in) :: xo ! number fraction of oxygen
        real(dp) :: a0, a1, a2, a3, a4, a5, Gamma

        a0 = 178d0
        a1 = 167.178104d0
        a2 = -3.973461d0
        a3 = -741.863826d0
        a4 = 876.516929d0
        a5 = -297.857813d0
        
        Gamma = &
             a0 + &
             a1*xo + &
             a2*xo*xo + &
             a3*xo*xo*xo + &
             a4*xo*xo*xo*xo + &
             a5*xo*xo*xo*xo*xo

        blouin_Gamma_melt_CO = Gamma
      end function blouin_Gamma_melt_CO
      
      real(dp) function blouin_delta_xne(Xin)
        real(dp), intent(in) :: Xin ! mass fraction
        real(dp) :: Xnew ! mass fraction
        real(dp) :: xne, dxne ! number fractions
        real(dp) :: a0, a1, a2, a3, a4, a5

        ! Convert input mass fraction to number fraction, assuming O/Ne mixture
        xne = (Xin/20d0)/(Xin/20d0 + (1d0 - Xin)/16d0)
        
        a0 = 0d0
        a1 = -0.120299d0
        a2 = 1.304399d0
        a3 = -1.722625d0
        a4 = 0.393996d0
        a5 = 0.144529d0

        dxne = &
             a0 + &
             a1*xne + &
             a2*xne*xne + &
             a3*xne*xne*xne + &
             a4*xne*xne*xne*xne + &
             a5*xne*xne*xne*xne*xne

        xne = xne + dxne

        ! Convert back to mass fraction
        Xnew = 20d0*xne/(20d0*xne + 16d0*(1d0-xne))
        
        blouin_delta_xne = Xnew - Xin
      end function blouin_delta_xne

      real(dp) function blouin_delta_xo_3component(xo,xne)
        real(dp), intent(in) :: xo, xne ! number fractions
        real(dp) :: a0, a1, a2, a3, a4, a5

        a0 = 0d0
        a1 = -0.311540d0
        a2 = 2.114743d0
        a3 = -1.661095d0
        a4 = -1.406005d0
        a5 = 1.263897d0

        blouin_delta_xo_3component = &
             (a0 + 0.640125d0*xne) + &
             (a1 + 2.218484d0*xne)*xo + &
             (a2 - 4.599227d0*xne)*xo*xo + &
             a3*xo*xo*xo + &
             a4*xo*xo*xo*xo + &
             a5*xo*xo*xo*xo*xo
        
      end function blouin_delta_xo_3component

      real(dp) function blouin_delta_xne_3component(xo,xne)
        real(dp), intent(in) :: xo, xne ! number fractions
        blouin_delta_xne_3component = min(0d0,-0.611587d0*xne + 0.782489d0*xne*xo)
      end function blouin_delta_xne_3component
      
      real(dp) function blouin_XNe_crit(GammaC)
        real(dp), intent(in) :: GammaC ! Carbon coupling
        real(dp) :: XNe1, XNe2, GammaC1, GammaC2, XNe_crit

        ! 2 points from the fit line defining the boundary
        ! (figure 3; Blouin, Daligault, & Saumon 2021 ApJL 911:L5)
        ! using right y-axis for mass fractions rather than number fractions
        XNe1 = 0.042d0
        XNe2 = 0.005d0
        GammaC1 = 130d0
        GammaC2 = 180d0

        XNe_crit = XNe1 + (GammaC - GammaC1)*(XNe2 - XNe1)/(GammaC2-GammaC1)

        ! set a floor of zero to avoid returning negative values for large GammaC
        blouin_XNe_crit = max(XNe_crit,0d0)
      end function blouin_XNe_crit

      subroutine calc_L_distill(s,sd,L_distill)
        type(star_info), pointer     :: s
        type(saved_data), intent(in) :: sd
        real(dp), intent(out)        :: L_distill
        
        real(dp) :: eps_tmp
        integer  :: k
        
        L_distill = 0d0
        do k = 1, s% nz
           eps_tmp = (sd% energy(k) - s% energy(k)) / s% dt
           L_distill = L_distill + eps_tmp*s% dm(k)
        end do
      end subroutine calc_L_distill
      
      subroutine update_model_ (s, kc_t, kc_b, do_brunt)

        use turb_info, only: set_mlt_vars
        use brunt, only: do_brunt_B
        use micro
        
        type(star_info), pointer :: s
        integer, intent(in)      :: kc_t
        integer, intent(in)      :: kc_b
        logical, intent(in)      :: do_brunt
        
        integer  :: ierr
        integer  :: kf_t
        integer  :: kf_b

        logical :: mask(s%nz)

        mask(:) = .true.

        ! Update the model to reflect changes in the abundances across
        ! cells kc_t:kc_b (the mask part of this call is unused, mask=true for all zones).
        ! Do updates at constant (P,T) rather than constant (rho,T).
        s%fix_Pgas = .true.
        call set_eos_with_mask(s, kc_t, kc_b, mask, ierr)
        if (ierr /= 0) then
           write(*,*) 'phase_separation: error from call to set_eos_with_mask'
           stop
        end if
        s%fix_Pgas = .false.
        
        ! Update opacities across cells kc_t:kc_b (this also sets rho_face
        ! and related quantities on faces kc_t:kc_b)        
        call set_micro_vars(s, kc_t, kc_b, &
             skip_eos=.TRUE., skip_net=.TRUE., skip_neu=.TRUE., skip_kap=.FALSE., ierr=ierr)
        if (ierr /= 0) then
           write(*,*) 'phase_separation: error from call to set_micro_vars'
           stop
        end if

        ! This is expensive, so only do it if we really need to.
        if(do_brunt) then
           ! Need to make sure we can set brunt for mix_outward calculation.
           if(.not. s% calculate_Brunt_B) then
              stop "phase separation requires s% calculate_Brunt_B = .true."
           end if
           call do_brunt_B(s, kc_t, kc_b, ierr) ! for unsmoothed_brunt_B
           if (ierr /= 0) then
              write(*,*) 'phase_separation: error from call to do_brunt_B'
              stop
           end if
        end if

        ! Finally update MLT for interior faces
        
        kf_t = kc_t
        kf_b = kc_b + 1
        
        call set_mlt_vars(s, kf_t+1, kf_b-1, ierr)
        if (ierr /= 0) then
           write(*,*) 'phase_separation: failed in call to set_mlt_vars during update_model_'
           stop
        endif
        
        ! Finish
        
        return
        
      end subroutine update_model_

      subroutine save_model_ (s, sd)
        type(star_info), pointer        :: s
        type(saved_data), intent(inout) :: sd

        sd% crystal_core_boundary_mass = s% crystal_core_boundary_mass
        sd% xa(:,1:s%nz) = s% xa(:,1:s%nz)
        sd% energy(1:s%nz) = s% energy(1:s%nz)
      end subroutine save_model_

      subroutine restore_composition_ (s, sd)
        type(star_info), pointer     :: s
        type(saved_data), intent(in) :: sd

        s% crystal_core_boundary_mass = sd% crystal_core_boundary_mass
        s% xa(:,1:s%nz) = sd% xa(:,1:s%nz)
      end subroutine restore_composition_

      subroutine alloc_saved_data_ (s, sd)
        type(star_info), pointer :: s
        type(saved_data)         :: sd

        ! Allocate cell data arrays
        allocate(sd%xa(s%species,s%nz))
        allocate(sd%energy(s%nz))
      end subroutine alloc_saved_data_

      subroutine smooth_eps_phase_sep(s,dt,ierr)
         type (star_info), pointer :: s
         real(dp), intent(in) :: dt
         integer, intent(out) :: ierr

         real(dp) :: integrated_luminosity
         integer :: k, kmid

         integrated_luminosity = dot_product(s% dm(1:s%nz), s% eps_phase_separation(1:s%nz))

         ! redistribute evenly through the inner half of the star
         do k = 1,s%nz
            if(s% q(k) < 0.5d0) then
               kmid = k
               exit
            end if
         end do

         s% eps_phase_separation(:) = 0d0
         s% eps_phase_separation(kmid:s%nz) = integrated_luminosity/s% m(kmid)

      end subroutine smooth_eps_phase_sep
      
      end module phase_separation



