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

      integer, parameter :: FIXED_PT_MODE = 5
      integer, parameter :: FIXED_DT_MODE = 6      

      real(dp), parameter :: eos_phase_boundary = 0.9d0
      
      contains


      subroutine do_phase_separation(s, ierr)
         use chem_def, only: chem_isos, ic12, io16
         use chem_lib, only: chem_get_iso_id
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         
         real(dp) :: dq_crystal, XO, XC, pad
         integer :: k, k_bound, k_new, kstart, net_ic12, net_io16
         logical :: do_premix

         do_premix = .true.

         if(s% phase(s% nz) < eos_phase_boundary) then
            s% crystal_core_boundary_mass = 0d0
            return
         end if

         net_ic12 = s% net_iso(ic12)
         net_io16 = s% net_iso(io16)
         
         ! Find zone of phase transition from liquid to solid
         k_bound = -1
         do k = s%nz,1,-1
            if(s% phase(k-1) <= eos_phase_boundary .and. s% phase(k) > eos_phase_boundary) then
               k_bound = k
               exit
            end if
         end do

         ! Check that we're still in C/O dominated material, otherwise skip phase separation
         XO = s% xa(net_io16,k_bound)
         XC = s% xa(net_ic12,k_bound)
         if (XO + XC < 0.9d0) return
         
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
            
            ! print *, "kstart, k_bound, phase(k_bound), phase(k_bound - 1)", kstart, k_bound, s%phase(k_bound), s%phase(k_bound - 1)
            ! print *, "mass(kstart), crystal_core_boundary_mass, mass(kstart+1)", s% m(kstart)/msun, s% crystal_core_boundary_mass/msun, s% m(kstart+1)/msun
         
            k_new = k_bound
            ! loop runs outward starting at previous crystallization boundary
            do k = kstart,1,-1
               ! Start by checking if this material should be crystallizing
               if(s% phase(k) <= eos_phase_boundary) then
                  k_new = k+1 ! one zone inward from where material becomes liquid
                  s% crystal_core_boundary_mass = s% m(k+1)
                  exit
               end if

               call move_one_zone(s,k,dq_crystal)
               ! crystallized out to k now, liquid starts at k-1.
               ! now mix the liquid material outward until stably stratified
               if(do_premix .and. dq_crystal > 0d0) then
                  call mix_outward(s, k-1)
               end if
               
            end do
            ! print *, "new k after starting phase sep", k_new
            ! print *, "phase(k+1), phase(k), phase(k-1)", s% phase(k_new+1), s% phase(k_new), s% phase(k_new-1)

            s% need_to_setvars = .true.
         end if

         ierr = 0
      end subroutine do_phase_separation

      subroutine move_one_zone(s,k,dq_crystal)
        use chem_def, only: chem_isos, ic12, io16
        use chem_lib, only: chem_get_iso_id
        type(star_info), pointer :: s
        integer, intent(in) :: k
        real(dp), intent(inout) :: dq_crystal
        
        real(dp) :: XC, XO, XC1, XO1, dXO, Xfac, dqsum
        integer :: net_ic12, net_io16
        integer :: update_mode(s%nz)
        
        update_mode(:) = FIXED_DT_MODE
        
        dq_crystal = dq_crystal + s% dq(k)
        
        net_ic12 = s% net_iso(ic12)
        net_io16 = s% net_iso(io16)
        
        XO = s% xa(net_io16,k)
        XC = s% xa(net_ic12,k)
        
        ! Try a real phase diagram (Blouin 2021)
        ! Need to rescale temporarily because phase diagram assumes XO + XC = 1
        Xfac = XO + XC
        XO = XO/Xfac
        XC = XC/Xfac

        dXO = blouin_delta_xo(XO)

        ! print *, "k, XO, dXO", k, XO, dXO
        s% xa(net_io16,k) = Xfac*(XO + dXO)
        s% xa(net_ic12,k) = Xfac*(XC - dXO)
        
        ! Redistribute change in X,O into zone k-1,
        ! conserving total mass of X,O
        XC1 = s% xa(net_ic12,k-1)
        XO1 = s% xa(net_io16,k-1)
        s% xa(net_ic12,k-1) = XC1 + Xfac*dXO * s% dq(k) / s% dq(k-1)
        s% xa(net_io16,k-1) = XO1 - Xfac*dXO * s% dq(k) / s% dq(k-1)

        call update_model_(s,update_mode,k-1,s%nz,.true.)
        
      end subroutine move_one_zone
      
      ! mix composition outward until no more negative molecular weight gradient
      subroutine mix_outward(s,kbot)
        use chem_def, only: chem_isos, ihe4, ic12, io16

        type(star_info), pointer :: s
        integer, intent(in)      :: kbot
        
        real(dp) :: avg_xa(s%species)
        real(dp) :: mass, XHe_out, dXC_top, dXC_bot, dXO_top, dXO_bot, B_term, grada, gradr
        integer :: k, l, ktop, net_ihe4, net_ic12, net_io16
        integer :: update_mode(s%nz)
        logical :: use_brunt
        
        update_mode(:) = FIXED_DT_MODE
        use_brunt = s% phase_separation_mixing_use_brunt
        net_ihe4 = s% net_iso(ihe4)
        net_ic12 = s% net_iso(ic12)
        net_io16 = s% net_iso(io16)

        do k=kbot,1,-1
           ktop = k

           mass = SUM(s%dm(ktop:kbot))
           do l = 1, s%species
              avg_xa(l) = SUM(s%dm(ktop:kbot)*s%xa(l,ktop:kbot))/mass
           end do

           ! some potential safeguards from conv_premix
           ! avg_xa = MAX(MIN(avg_xa, 1._dp), 0._dp)
           ! avg_xa = avg_xa/SUM(avg_xa)

           XHe_out = max(s%xa(net_ihe4,ktop),s%xa(net_ihe4,ktop-1))
           if(XHe_out < s% eos_rq% mass_fraction_limit_for_Skye) then
              ! ok to mix all species
              do l = 1, s%species
                 s%xa(l,ktop:kbot) = avg_xa(l)
              end do
           else
              ! Mixing He can cause energy problems for eps_phase_separation
              ! when using Skye, so once we encounter enough He that it is
              ! included in Skye energy calculation, stop mixing all species.
              ! Instead, just flatten out the O16 profile, and mix in exchange
              ! for C12.
              dXO_top = avg_xa(net_io16) - s%xa(net_io16,ktop)
              dXO_bot = avg_xa(net_io16) - s%xa(net_io16,kbot)
              s%xa(net_io16,ktop:kbot) = avg_xa(net_io16)
              dXC_top = -dXO_top
              dXC_bot = -dXO_bot
              s%xa(net_ic12,ktop) = s%xa(net_ic12,ktop) + dXC_top
              s%xa(net_ic12,ktop+1:kbot) = s%xa(net_ic12,kbot) + dXC_bot
           end if


           ! updates, eos, opacities, mu, etc now that abundances have changed,
           ! but only in the cells near the boundary where we need to check here.
           ! Will call full update over mixed region after exiting loop.
           call update_model_(s, update_mode, ktop-1, ktop+1, use_brunt)

           if(use_brunt) then
              B_term = s% unsmoothed_brunt_B(k)
              grada = s% grada_face(k)
              gradr = s% gradr(k)
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
        call update_model_(s, update_mode, ktop, kbot, .true.)
        ! print *, "kbot, ktop, m(ktop), mu_in, mu_out", kbot, ktop, s% m(ktop)/msun, s% mu(ktop), s% mu(ktop-1)
       
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
      
      subroutine update_model_ (s, update_mode, kc_t, kc_b, do_brunt)

        use turb_info, only: set_mlt_vars
        use brunt, only: do_brunt_B
        use micro
        
        type(star_info), pointer :: s
        integer, intent(in)      :: update_mode(:)
        integer, intent(in)      :: kc_t
        integer, intent(in)      :: kc_b
        logical, intent(in)      :: do_brunt
        
        integer  :: ierr
        integer  :: kf_t
        integer  :: kf_b
        
        ! Update the model to reflect changes in the abundances across
        ! cells kc_t:kc_b
        
        call set_eos_with_mask(s, kc_t, kc_b, update_mode==FIXED_DT_MODE, ierr)
        if (ierr /= 0) then
           write(*,*) 'phase_separation: error from call to set_eos_with_mask'
           stop
        end if
        
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



