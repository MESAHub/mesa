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

      private
      public :: do_phase_separation

      logical, parameter :: dbg = .false.

      ! offset to higher phase than 0.5 to avoid interference
      ! between phase separation mixing and latent heat for Skye.
      real(dp), parameter :: eos_phase_boundary = 0.9d0
      
      contains

      subroutine do_phase_separation(s, dt, ierr)
         type (star_info), pointer :: s
         real(dp), intent(in) :: dt
         integer, intent(out) :: ierr

         if(s% phase_separation_option == 'CO') then
            call do_2component_phase_separation(s, dt, 'CO', ierr)
         else if(s% phase_separation_option == 'ONe') then
            call do_2component_phase_separation(s, dt, 'ONe', ierr)
         else if(s% phase_separation_option == 'distillation') then
            call do_distillation(s, dt, ierr)
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
         
         real(dp) :: dq_crystal, XNe, XO, XC, pad
         integer :: k, k_bound, kstart, net_ic12, net_io16, net_ine20
         logical :: save_Skye_use_ion_offsets

         s% eps_phase_separation(1:s%nz) = 0d0
         
         if(s% phase(s% nz) < eos_phase_boundary) then
            s% crystal_core_boundary_mass = 0d0
            return
         end if

         ! Set phase separation mixing mass negative at beginning of phase separation
         s% phase_sep_mixing_mass = -1d0

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
            dq_crystal = 0d0

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

               call move_one_zone(s,k,components,dq_crystal)
               ! crystallized out to k now, liquid starts at k-1.
               ! now mix the liquid material outward until stably stratified
               if(dq_crystal > 0d0) then
                  call mix_outward(s, k-1)
               end if
               
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
         use chem_def, only: chem_isos, ic12, io16, ine20, ine22
         use chem_lib, only: chem_get_iso_id
         type (star_info), pointer :: s
         real(dp), intent(in) :: dt
         integer, intent(out) :: ierr
         
         real(dp) :: XNe, XNe_out, XO, XC, pad, distill_critical_XNe
         integer :: k, k_bound, kstart, net_ic12, net_io16, net_ine20, net_ine22
         logical :: save_Skye_use_ion_offsets

         distill_critical_XNe = 0.3143d0
         
         s% eps_phase_separation(1:s%nz) = 0d0
         
         if(s% phase(s% nz) < eos_phase_boundary) then
            s% crystal_core_boundary_mass = 0d0
            return
         end if

         ! Set phase separation mixing mass negative at beginning of phase separation
         s% phase_sep_mixing_mass = -1d0

         net_ic12 = s% net_iso(ic12)
         net_io16 = s% net_iso(io16)
         net_ine20 = s% net_iso(ine20)
         net_ine22 = s% net_iso(ine22)

         if(net_ine22 < 0) then
            write(*,*) 'distillation requires ne22 be included in net'
            stop
         end if
         
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
         XNe = s% xa(net_ine20,k_bound) + s% xa(net_ine22,k_bound)
         
         if(XC + XO + XNe < 0.9d0) return

         ! check if there's enough Ne for distillation,
         ! or if we should just do 2component C/O phase separation
         if(XNe < 0.015d0) then ! TODO: make this composition check depend on GammaC
            call do_2component_phase_separation(s, dt, 'CO', ierr)
            return
         end if
         
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

               XC = s% xa(net_ic12,k)
               XO = s% xa(net_io16,k)
               XNe = s% xa(net_ine20,k) + s% xa(net_ine22,k)
               XNe_out = s% xa(net_ine20,k-1) + s% xa(net_ine22,k-1)
               ! must distill to xNe = 0.2 (XNe = 0.3143) and xC = 0.8 (no O),
               ! as long as there is enough Ne in the next zone out to proceed
               ! TODO: make check on XNe_out more physical
               do while(XNe < distill_critical_XNe .and. XNe_out > 0.01d0)
                  call distill_at_boundary(s,k,distill_critical_XNe)
                  ! mix from zone k-1 outward
                  call mix_outward(s, k-1)
                  XO = s% xa(net_io16,k)
                  XNe = s% xa(net_ine20,k) + s% xa(net_ine22,k)
                  XNe_out = s% xa(net_ine20,k-1) + s% xa(net_ine22,k-1)
               end do
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
      end subroutine do_distillation

      subroutine distill_at_boundary(s,k,distill_critical_XNe)
        use chem_def, only: chem_isos, ic12, io16, ine20, ine22
        use chem_lib, only: chem_get_iso_id
        type(star_info), pointer :: s
        integer, intent(in) :: k
        real(dp), intent(in) :: distill_critical_XNe
        
        real(dp) :: XC, XO, XNe, XC_out, XO_out, XNe_out, Xout_sum, Xfac
        real(dp) :: dXC, dXO, dXNe_tot, dXNe20, dXNe22, dq_ratio
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

        print *, "distill at boundary; k, XC, XO, XNe, XC_out, XO_out, XNe_out", &
             k, XC, XO, XNe, XC_out, XO_out, XNe_out
        
        ! Need to rescale temporarily because phase diagram assumes XC + XO + XNe = 1?
        Xfac = XO + XC + XNe
        ! TODO: think about whether we want to do this rescaling
        !XC = XC/Xfac
        !XO = XO/Xfac
        !XNe = XNe/Xfac

        ! Net effect of distillation is that crystals enriched in oxygen float upward.
        ! Need to limit toward xNe = 0.2, xC = 0.8. Start by pushing O outward in exchange
        ! for C/Ne mixture until O is depleted. Then exchange C/Ne until reaching critical
        ! Ne value.

        if(XO > 0d0 .and. XNe < distill_critical_XNe) then
           ! exchange O for a mixture of C and Ne
           Xout_sum = XC_out + XNe_out
           if(XO <= Xout_sum) then
              ! Can exchange all O in zone k for C and Ne
              dXO = s% xa(net_io16,k)
           else
              ! can only exchange as much O as we have C and Ne in external zone,
              ! but next iteration should get rid of the rest of the O
              dXO = Xout_sum
           end if
           ! make relative fractions of exchanged C12, Ne20, and Ne22 proportional to fractions in zone k-1,
           ! with dXC + dXNe20 + dXNe22 = dXO
           dXC    = dXO*s% xa(net_ic12,k-1) /Xout_sum
           dXNe20 = dXO*s% xa(net_ine20,k-1)/Xout_sum
           dXNe22 = dXO*s% xa(net_ine22,k-1)/Xout_sum
           
           s% xa(net_ic12,k) = s% xa(net_ic12,k) + dXC
           s% xa(net_io16,k) = s% xa(net_io16,k) - dXO
           s% xa(net_ine20,k) = s% xa(net_ine20,k) + dXNe20
           s% xa(net_ine22,k) = s% xa(net_ine22,k) + dXNe22

           ! use dq_ratio to conserve total mass of each element,
           ! accounting for the fact that zone k-1 has different mass than zone k.
           s% xa(net_ic12,k-1) = s% xa(net_ic12,k-1) - dXC*dq_ratio
           s% xa(net_io16,k-1) = s% xa(net_io16,k-1) + dXO*dq_ratio
           s% xa(net_ine20,k-1) = s% xa(net_ine20,k-1) - dXNe20*dq_ratio
           s% xa(net_ine22,k-1) = s% xa(net_ine22,k-1) - dXNe22*dq_ratio
        else if(XNe < distill_critical_XNe) then
           ! continue increasing XNe until reaching critical value
           if(XNe_out <= (distill_critical_XNe - XNe)) then
              ! pull all Ne from zone k-1 into zone k, exchange for C
              dXNe20 = s% xa(net_ine20,k-1)
              dXNe22 = s% xa(net_ine22,k-1)
              dXC = dXNe20 + dXNe22
           else
              ! Only exchange as much as needed to reach critical Ne
              dXNe_tot = distill_critical_XNe - XNe
              dXNe20 = dXNe_tot*s% xa(net_ine20,k-1)/XNe_out
              dXNe22 = dXNe_tot*s% xa(net_ine22,k-1)/XNe_out
              dXC = dXNe_tot
           end if
           
           s% xa(net_ine20,k) = s% xa(net_ine20,k) + dXNe20
           s% xa(net_ine22,k) = s% xa(net_ine22,k) + dXNe22
           s% xa(net_ic12,k) = s% xa(net_ic12,k) - dXC
           
           ! use dq_ratio to conserve total mass of each element,
           ! accounting for the fact that zone k-1 has different mass than zone k.
           s% xa(net_ine20,k-1) = s% xa(net_ine20,k-1) - dXNe20*dq_ratio
           s% xa(net_ine22,k-1) = s% xa(net_ine22,k-1) - dXNe22*dq_ratio
           s% xa(net_ic12,k-1) = s% xa(net_ic12,k-1) + dXC*dq_ratio
        else ! for debugging
           write(*,*) 'should not have ended up here...'
           write(*,*) 'XC, XO, XNe', XC, XO, XNe
           stop
        end if

        call update_model_(s,k-1,s%nz,.true.)
        
      end subroutine distill_at_boundary

      subroutine move_one_zone(s,k,components,dq_crystal)
        use chem_def, only: chem_isos, ic12, io16, ine20
        use chem_lib, only: chem_get_iso_id
        type(star_info), pointer :: s
        integer, intent(in) :: k
        character (len=*), intent(in) :: components
        real(dp), intent(inout) :: dq_crystal
        
        real(dp) :: XC, XO, XNe, XC1, XO1, XNe1, dXO, dXNe, Xfac
        integer :: net_ic12, net_io16, net_ine20
        
        dq_crystal = dq_crystal + s% dq(k)
        
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
      
      ! mix composition outward until reaching stable composition profile
      subroutine mix_outward(s,kbot)
        type(star_info), pointer :: s
        integer, intent(in)      :: kbot
        
        real(dp) :: avg_xa(s%species)
        real(dp) :: mass, B_term, grada, gradr
        integer :: k, l, ktop
        logical :: use_brunt
        
        use_brunt = s% phase_separation_mixing_use_brunt

        do k=kbot,1,-1
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
        call update_model_(s, ktop, kbot, .true.)
       
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
      
      
      end module phase_separation



