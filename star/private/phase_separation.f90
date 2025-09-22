      module phase_separation

      use star_private_def
      use const_def

      implicit none

      logical, parameter :: dbg = .false.

      ! offset to higher phase than 0.5 to avoid interference
      ! between phase separation mixing and latent heat for Skye.
      real(dp), parameter :: eos_phase_boundary = 0.9d0
      private
      public :: do_phase_separation

      contains

      subroutine do_phase_separation(s, dt, ierr)
         type (star_info), pointer :: s
         real(dp), intent(in) :: dt
         integer, intent(out) :: ierr

         ! 'CO' or 'ONe' will implement 2-species phase separation, for 'ONe' 22Ne is included
         if(s% phase_separation_option == 'CO') then
            call do_2component_phase_separation(s, dt, 'CO', ierr)
         else if(s% phase_separation_option == 'ONe') then
            call do_2component_phase_separation(s, dt, 'ONe', ierr)
         else if(s% phase_separation_option == '3c') then
            call do_2component_phase_separation(s, dt, '3c', ierr)
         else
            write(*,*) 'invalid phase_separation_option'
            stop
         end if
      end subroutine do_phase_separation

      subroutine do_2component_phase_separation(s, dt, components, ierr)
         use chem_def, only: chem_isos, ic12, io16, ine20, ine22, ina23, img24
         use chem_lib, only: chem_get_iso_id
         type (star_info), pointer :: s
         real(dp), intent(in) :: dt
         character (len=*), intent(in) :: components
         integer, intent(out) :: ierr
         
         real(dp) :: XNe20, XNe22, XO, XC, XNa, XMg , pad
         integer :: k, k_bound, kstart, net_ic12, net_io16, net_ine20, net_ine22, net_ina23, net_img24
         logical :: save_Skye_use_ion_offsets

         ! Set phase separation mixing mass negative at beginning of phase separation
         s% phase_sep_mixing_mass = -1d0
         s% eps_phase_separation(1:s%nz) = 0d0
         
         if(s% phase(s% nz) < eos_phase_boundary) then !!! prevent to move the core size inwards if the core is suddently "melted" leaving everything liquid under phi<0.9
            if (s% crystal_core_boundary_mass>0d0)then
               s% crystal_core_boundary_mass=s% crystal_core_boundary_mass
               return
            else 
               s% crystal_core_boundary_mass = 0d0
               return
            end if
         end if
         
         net_ic12 = s% net_iso(ic12)
         net_io16 = s% net_iso(io16)
         net_ine20 = s% net_iso(ine20)
         net_ine22 = s% net_iso(ine22)
         net_ina23 = s% net_iso(ina23)
         net_img24 = s% net_iso(img24)

         
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
         XNe20 = s% xa(net_ine20,k_bound)
         XNe22 = s% xa(net_ine22,k_bound)
         XNa = s% xa(net_ina23,k_bound)
         XMg = s% xa(net_img24,k_bound)
         
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
                  if (s% crystal_core_boundary_mass>s% m(k+1)) then
                     s% crystal_core_boundary_mass=s% crystal_core_boundary_mass
                      exit
                  else
                     s% crystal_core_boundary_mass = s% m(k+1)
                     exit
                  end if                 
               end if

               call move_one_zone(s,k,components)
               ! crystallized out to k now, liquid starts at k-1.
               ! now mix the liquid material outward until stably stratified
               call mix_outward(s, k-1, 0)
               
            end do

            call update_model_(s,1,s%nz,.false.)

            do k=1,s% nz
               s% eps_phase_separation(k) = (s% eps_phase_separation(k) - s% energy(k)) / dt
            end do

            s% eos_rq% Skye_use_ion_offsets = save_Skye_use_ion_offsets
            s% need_to_setvars = .true.
         end if
         ierr = 0
      end subroutine do_2component_phase_separation


      subroutine move_one_zone(s,k,components)
        use chem_def, only: chem_isos, ic12, io16, ine20, ine22, ina23, img24
        use chem_lib, only: chem_get_iso_id
        type(star_info), pointer :: s
        integer, intent(in) :: k
        real(dp), dimension(2) :: dXNe
        real(dp), dimension(4) :: Dd 
        real(dp) :: dx1_
        character (len=*), intent(in) :: components
        
        real(dp) :: XC, XO, XNe20, XNe22, XNa, XMg, XC1, XO1, XNe120, XNe122, XNa1, XMg1, dXO, Xfac
        integer :: net_ic12, net_io16, net_ine20, net_ine22, net_ina23, net_img24

        net_ic12 = s% net_iso(ic12)
        net_io16 = s% net_iso(io16)
        net_ine20 = s% net_iso(ine20)
        net_ine22 = s% net_iso(ine22)
        net_ina23 = s% net_iso(ina23)
        net_img24 = s% net_iso(img24)
        
        XO = s% xa(net_io16,k)
        XC = s% xa(net_ic12,k)
        XNe20 = s% xa(net_ine20,k)
        XNe22 = s% xa(net_ine22,k)
        XNa = s% xa(net_ina23,k)
        XMg = s% xa(net_img24,k)  
  
        if (components .ne. '3c') then   
        if(XO + XC > 0.7d0 .and. XC > XNe20 + XNe22) then
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
           !write(*,*) 'phase CO',XO,XC,XNe20+XNe22
        else if(XO + XNe20 + XNe22> 0.7d0 .and. XNe20 + XNe22 > XC) then
        
           ! Call Blouin phase diagram.
           ! Need to rescale temporarily because phase diagram assumes XO + XNe = 1
           Xfac = XO + XNe20 + XNe22
           XO = XO/Xfac
           XNe20 = XNe20/Xfac
           XNe22 = XNe22/Xfac
           
           dXNe = blouin_delta_xne(XNe20,XNe22)
          ! write(*,*) 'dXNe', dXNe
           
           s% xa(net_ine20,k) = Xfac*(XNe20 + dXNe(1))
           s% xa(net_ine22,k) = Xfac*(XNe22 + dXNe(2))
           s% xa(net_io16,k) = Xfac*(XO - sum(dXNe))
           
           ! Redistribute change in Ne,O into zone k-1,
           ! conserving total mass of Ne,O
           XO1 = s% xa(net_io16,k-1)
           XNe120 = s% xa(net_ine20,k-1)
           XNe122 = s% xa(net_ine22,k-1)
           s% xa(net_io16,k-1) = XO1 + Xfac*sum(dXNe) * s% dq(k) / s% dq(k-1)
           s% xa(net_ine20,k-1) = XNe120 - Xfac*dXNe(1) * s% dq(k) / s% dq(k-1)
           s% xa(net_ine22,k-1) = XNe122 - Xfac*dXNe(2) * s% dq(k) / s% dq(k-1)
           !write(*,*) 'phase ONe',XO,XC,XNe20+XNe22
        end if

        else if (components == '3c') then
          ! check the abundances to decide which table use for interpolation
          if (XO + XC + XNe20 + XNe22 > 0.7d0 .and. XC > XMg .and. XC > XNa) then
           Xfac = XO + XC + XNe20 + XNe22
           XO = XO/Xfac
           XC = XC/Xfac
           XNe20 = XNe20/Xfac
           XNe22 = XNe22/Xfac
           
           ! call the deltas resulting from interpolation (in mass fraction)
           call medin_cumming_3p_d_cone(XC,XO,XNe20,XNe22,Dd)

           ! apply fractionation as given by the deltas from interpolation
           s% xa(net_ic12,k) = Xfac*(XC + Dd(1))
           s% xa(net_io16,k) = Xfac*(XO + Dd(2))
           s% xa(net_ine20,k) = Xfac*(XNe20 + Dd(3))
           s% xa(net_ine22,k) = Xfac*(XNe22 + Dd(4))
           
           XC1 = s% xa(net_ic12,k-1)
           XO1 = s% xa(net_io16,k-1)
           XNe120 = s% xa(net_ine20,k-1)
           XNe122 = s% xa(net_ine22,k-1)

           s% xa(net_ic12,k-1) = XC1 - Xfac*Dd(1) * s% dq(k) / s% dq(k-1)
           s% xa(net_io16,k-1) = XO1 - Xfac*Dd(2) * s% dq(k) / s% dq(k-1)
           s% xa(net_ine20,k-1) = XNe120 - Xfac*(Dd(3)) * s% dq(k) / s% dq(k-1)
           s% xa(net_ine22,k-1) = XNe122 - Xfac*(Dd(4)) * s% dq(k) / s% dq(k-1)
          ! write(*,*) 'phase 3 CONe abundances',XC,XO,XNe20+XNe22 

         else if (XO  + XNe20 + XNe22 + XMg > 0.7d0 .and. XMg > XC .and. XMg > XNa) then
           Xfac = XO + XNe20 + XNe22 + XMg
           XMg = XMg/Xfac
           XO = XO/Xfac
           XNe20 = XNe20/Xfac
           XNe22 = XNe22/Xfac

          ! call the deltas resulting from interpolation (in mass fraction)
           call medin_cumming_3p_d_neomg(XMg,XO,XNe20,XNe22,Dd)

           ! apply fractionation as given by the deltas from interpolation
           
           s% xa(net_img24,k) = Xfac*(XMg + Dd(1))
           s% xa(net_io16,k) = Xfac*(XO + Dd(2))
           s% xa(net_ine20,k) = Xfac*(XNe20 + Dd(3))
           s% xa(net_ine22,k) = Xfac*(XNe22 + Dd(4))
           
           XMg1 = s% xa(net_img24,k-1)
           XO1 = s% xa(net_io16,k-1)
           XNe120 = s% xa(net_ine20,k-1)
           XNe122 = s% xa(net_ine22,k-1)

           s% xa(net_img24,k-1) = XMg1 - Xfac*Dd(1) * s% dq(k) / s% dq(k-1)
           s% xa(net_io16,k-1) = XO1 - Xfac*Dd(2) * s% dq(k) / s% dq(k-1)
           s% xa(net_ine20,k-1) = XNe120 - Xfac*Dd(3) * s% dq(k) / s% dq(k-1)
           s% xa(net_ine22,k-1) = XNe122 - Xfac*Dd(4) * s% dq(k) / s% dq(k-1)
          ! write(*,*) 'phase 3 ONeMg abundances',XO,XNe20+XNe22,XMg 

         else if (XO  + XNe20 + XNe22 + XNa > 0.7d0 .and. XNa > XC .and. XNa > XMg) then
           Xfac = XO + XNe20 + XNe22 + XNa
           XNa = XNa/Xfac
           XO = XO/Xfac
           XNe20 = XNe20/Xfac
           XNe22 = XNe22/Xfac

           ! call the deltas resulting from interpolation (in mass fraction)
           call medin_cumming_3p_d_onena(XNa,XO,XNe20,XNe22,Dd)

           ! apply fractionation as given by the deltas from interpolation
         
           s% xa(net_ina23,k) = Xfac*(XNa + Dd(1))
           s% xa(net_io16,k) = Xfac*(XO + Dd(2))
           s% xa(net_ine20,k) = Xfac*(XNe20 + Dd(3))
           s% xa(net_ine22,k) = Xfac*(XNe22 + Dd(4))
           
           XNa1 = s% xa(net_ina23,k-1)
           XO1 = s% xa(net_io16,k-1)
           XNe120 = s% xa(net_ine20,k-1)
           XNe122 = s% xa(net_ine22,k-1)

           s% xa(net_ina23,k-1) = XNa1 - Xfac*Dd(1) * s% dq(k) / s% dq(k-1)
           s% xa(net_io16,k-1) = XO1 - Xfac*Dd(2) * s% dq(k) / s% dq(k-1)
           s% xa(net_ine20,k-1) = XNe120 - Xfac*Dd(3) * s% dq(k) / s% dq(k-1)
           s% xa(net_ine22,k-1) = XNe122 - Xfac*Dd(4) * s% dq(k) / s% dq(k-1)
          ! write(*,*) 'phase 3 ONeNa abundances',XO,XNe20+XNe22,XNa 
         
         else if (XC  + XO + XMg > 0.7d0 .and. XMg > XNa .and. XMg > XNe20+XNe22) then
           Xfac = XC + XO + XMg
           XC = XC/Xfac
           XO = XO/Xfac
           XMg = XMg/Xfac

           ! call the deltas resulting from interpolation (in mass fraction)
           call medin_cumming_3p_d_comg(XC,XMg,XO,Dd)

           ! apply fractionation as given by the deltas from interpolation
      
           s% xa(net_ic12,k) = Xfac*(XC + Dd(1))
           s% xa(net_img24,k) = Xfac*(XMg + Dd(2))
           s% xa(net_io16,k) = Xfac*(XO - (Dd(1) + Dd(2)))
           
           XC1 = s% xa(net_ic12,k-1)
           XO1 = s% xa(net_io16,k-1)
           XMg1 = s% xa(net_img24,k-1)

           s% xa(net_ic12,k-1) = XC1 - Xfac*Dd(1) * s% dq(k) / s% dq(k-1)
           s% xa(net_img24,k-1) = XMg1 - Xfac*Dd(2) * s% dq(k) / s% dq(k-1)
           s% xa(net_io16,k-1) = XO1 + Xfac*(Dd(1)+Dd(2)) * s% dq(k) / s% dq(k-1)
          ! write(*,*) 'phase 3 COMg abundances',XC,XO,XMg 

         end if
        end if

        call update_model_(s,k-1,s%nz,.true.)
        
      end subroutine move_one_zone

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

      function blouin_delta_xne(Xin20,Xin22)
        real(dp), intent(in) :: Xin20, Xin22! mass fraction
        real(dp) :: Xnew1, Xnew2 ! mass fraction
        real(dp) :: xne, dxne, xne1, xne2 ! number fractions
        real(dp) :: a0, a1, a2, a3, a4, a5

        real(dp), dimension(2) :: blouin_delta_xne

        ! Convert input mass fraction to number fraction, assuming O/Ne mixture
        xne1 =(Xin20/20d0)/(Xin20/20d0 + Xin22/22d0 + (1d0 - Xin20 - Xin22)/16d0)
        xne2 =(Xin22/22d0)/(Xin20/20d0 + Xin22/22d0 + (1d0 - Xin20 - Xin22)/16d0)

        ! isotope 22Ne is added to the Ne separation along with 20Ne
        xne =((Xin22/22d0)+(Xin20/20d0))/(Xin20/20d0 + Xin22/22d0 + (1d0 - Xin20 - Xin22)/16d0)
        
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

        xne1 = xne1 + dxne*xne1/xne
        xne2 = xne2 + dxne*xne2/xne
        xne = xne1 + xne2
        ! Convert back to mass fraction
        Xnew1 = (20d0*xne1)/(20d0*xne1 + 22d0*xne2 + 16d0*(1d0-xne))
        Xnew2 = (22d0*xne2)/(20d0*xne1 + 22d0*xne2 + 16d0*(1d0-xne))
        
        blouin_delta_xne(1) = Xnew1 - (Xin20)
        blouin_delta_xne(2) = Xnew2 - (Xin22)
      end function blouin_delta_xne
 

     subroutine tab_interp_medin_cumming_dx1(x1_,x2_,components,dx1_)
            use interp_2D_lib_db, only: interp_mkbicub_db, interp_evbicub_db
            use const_def, only: mesa_data_dir
            use utils_lib, only: mesa_error, mkdir, is_bad
            implicit none
            integer, parameter :: num_x1 = 998, num_x2 = 998
            integer :: ilinx,iliny,ibcxmin,ibcxmax,ibcymin,ibcymax,iounit,ict(6),ierr,i,j,k
            real(dp) :: bcxmin(num_x1), bcxmax(num_x1)
            real(dp) :: bcymin(num_x2), bcymax(num_x2)
            real(dp), pointer, dimension(:) :: x1_l, x2_l, deltax1_sob_f1, deltax1_sob_values
            real(dp), pointer :: deltax1_sob_f(:,:,:)
            real(dp) :: deltax1,x1l,x2l
            real(dp), intent(in) :: x1_,x2_        ! target of this interpolation
            character (len=*), intent(in) :: components
            real(dp) :: fval(6)         ! output data
            real(dp), intent(out) :: dx1_
            integer :: ier

            ict = 0
            ict(1) = 1      

            iounit=999
      ! setup interpolation table for x1 x2 dx1 
            if (components=='CONe') then
               open(unit=iounit, file='CONe_deltaC.dat', action='read',status='old')
            else if  (components=='NeOMg') then
               open(unit=iounit, file='NeOMg_deltaMg.dat', action='read',status='old')
            else if  (components=='ONeNa') then
               open(unit=iounit, file='ONeNa_deltaNa.dat', action='read',status='old')
            else if  (components=='COMg') then
               open(unit=iounit, file='COMg_deltaC.dat', action='read',status='old')
            end if
            allocate(x1_l(num_x1), x2_l(num_x2), &
            deltax1_sob_f1(4*num_x1*num_x2))
            deltax1_sob_f(1:4,1:num_x1,1:num_x2) => &
            deltax1_sob_f1(1:4*num_x1*num_x2)
   do j=1,num_x1
      do i=1,num_x2
         read(iounit,*) x1l, x2l, deltax1
         x1_l(j)=x1l
         if (j == 1) then
            x2_l(i) =x2l 
         end if
         deltax1_sob_f(1,j,i) = deltax1
      end do
   end do
   close(iounit)
   ! just use "not a knot" bc's at edges of tables
   ibcxmin = 0; bcxmin(1:num_x1) = 0
   ibcxmax = 0; bcxmax(1:num_x1) = 0
   ibcymin = 0; bcymin(1:num_x2) = 0
   ibcymax = 0; bcymax(1:num_x2) = 0
   call interp_mkbicub_db( &
      x1_l, num_x1, x2_l, num_x2, deltax1_sob_f1, num_x1, &
      ibcxmin,bcxmin,ibcxmax,bcxmax, &
      ibcymin,bcymin,ibcymax,bcymax, &
      ilinx,iliny,ierr)
   if (ierr /= 0) then
      write(*,*) 'interp_mkbicub_db error'
      ierr = -1
      call mesa_error(__FILE__,__LINE__)
   end if

   do j=1,num_x1
      do i=1,num_x2
         do k=1,4
            if (is_bad(deltax1_sob_f(k,j,i))) then
            write(*,*) 'deltax1_sob_f', i, j, k, deltax1_sob_f(k,j,i)
            end if
         end do
      end do
   end do

   call interp_evbicub_db( &
            x1_, x2_, x1_l, num_x1, x2_l, num_x2, &
            ilinx, iliny, deltax1_sob_f1, num_x1, ict, fval, ier)

    dx1_=fval(1)  ! delta_x1 from 2d interpolation

      end subroutine tab_interp_medin_cumming_dx1


      subroutine tab_interp_medin_cumming_dx2(x1_,x2_,components,dx2_)
            !use utils_lib
            use interp_2D_lib_db, only: interp_mkbicub_db, interp_evbicub_db
            use const_def, only: mesa_data_dir
            use utils_lib, only: mesa_error, mkdir, is_bad
            implicit none
            integer, parameter :: num_x1 = 998, num_x2 = 998
            integer :: ilinx,iliny,ibcxmin,ibcxmax,ibcymin,ibcymax,iounit,ict(6),ierr,i,j,k
            real(dp) :: bcxmin(num_x1), bcxmax(num_x1)
            real(dp) :: bcymin(num_x2), bcymax(num_x2)
            real(dp), pointer, dimension(:) :: x1_l, x2_l, deltax1_sob_f1, deltax1_sob_values
            real(dp), pointer :: deltax1_sob_f(:,:,:)
            real(dp) :: deltax1,x1l,x2l
            real(dp), intent(in) :: x1_,x2_        ! target of this interpolation
            character (len=*), intent(in) :: components
            real(dp) :: fval(6)         ! output data
            real(dp), intent(out) :: dx2_
            integer :: ier

            ict = 0
            ict(1) = 1      

            iounit=998
      ! setup interpolation table for tau sob eta
            if (components=='CONe') then
               open(unit=iounit, file='CONe_deltaO.dat', action='read',status='old')
            else if  (components=='NeOMg') then
               open(unit=iounit, file='NeOMg_deltaO.dat', action='read',status='old')
            else if  (components=='ONeNa') then
               open(unit=iounit, file='ONeNa_deltaO.dat', action='read',status='old')
            else if  (components=='COMg') then
               open(unit=iounit, file='COMg_deltaMg.dat', action='read',status='old')
            end if
            allocate(x1_l(num_x1), x2_l(num_x2), &
            deltax1_sob_f1(4*num_x1*num_x2))
            deltax1_sob_f(1:4,1:num_x1,1:num_x2) => &
            deltax1_sob_f1(1:4*num_x1*num_x2)
   do j=1,num_x1
      do i=1,num_x2
         read(iounit,*) x1l, x2l, deltax1
         x1_l(j)=x1l
         if (j == 1) then
            x2_l(i) =x2l 
         end if
         deltax1_sob_f(1,j,i) = deltax1
      end do
   end do
   close(iounit)
   ! just use "not a knot" bc's at edges of tables
   ibcxmin = 0; bcxmin(1:num_x1) = 0
   ibcxmax = 0; bcxmax(1:num_x1) = 0
   ibcymin = 0; bcymin(1:num_x2) = 0
   ibcymax = 0; bcymax(1:num_x2) = 0

   call interp_mkbicub_db( &
      x1_l, num_x1, x2_l, num_x2, deltax1_sob_f1, num_x1, &
      ibcxmin,bcxmin,ibcxmax,bcxmax, &
      ibcymin,bcymin,ibcymax,bcymax, &
      ilinx,iliny,ierr)
   if (ierr /= 0) then
      write(*,*) 'interp_mkbicub_db error'
      ierr = -1
      call mesa_error(__FILE__,__LINE__)
   end if

   do j=1,num_x1
      do i=1,num_x2
         do k=1,4
            if (is_bad(deltax1_sob_f(k,j,i))) then
            write(*,*) 'deltax1_sob_f', i, j, k, deltax1_sob_f(k,j,i)
            end if
         end do
      end do
   end do
 
   call interp_evbicub_db( &
            x1_, x2_, x1_l, num_x1, x2_l, num_x2, &
            ilinx, iliny, deltax1_sob_f1, num_x1, ict, fval, ier)

    dx2_=fval(1)  ! delta_x2 from 2d interpolation

      end subroutine tab_interp_medin_cumming_dx2

       
      subroutine medin_cumming_3p_d_cone(X1,X2,X3_1,X3_2,Dd)
        real(dp), intent(in) :: X1, X2, X3_1, X3_2 ! mass fraction
        real(dp), dimension(4),intent(out) :: Dd
        real(dp) :: Xnew1, Xnew2, Xnew3_1, Xnew3_2, Xfac ! mass fraction
        real(dp) :: xc, dxc, xo, dxo, xne1, xne2 ! number fractions
        real(dp) :: dx1_,dx2_
        integer :: i,j

        Xfac = X1 + X2 + X3_1 + X3_2

        xc = (X1/12)/(X1/12 + X2/16 + X3_1/20 + X3_2/22)
        xo = (X2/16)/(X1/12 + X2/16 + X3_1/20 + X3_2/22)

        xne1 = (X3_1/20)/(X1/12 + X2/16 + X3_1/20 + X3_2/22)
        xne2 = (X3_2/22)/(X1/12 + X2/16 + X3_1/20 + X3_2/22)

        call tab_interp_medin_cumming_dx1(xc,xo,'CONe',dx1_)
        call tab_interp_medin_cumming_dx2(xc,xo,'CONe',dx2_)
        dxc=dx1_
        dxo=dx2_

        !write(*,*) 'delta_xc: ',dxc,' delta_xo: ', dxo      
        
        xc = xc + dxc
        xo = xo + dxo
        
        ! convert deltas in number fraction to mass fraction
        Xnew1 = 12*xc/(12*xc + 16*xo + 20*(1-xc-xo)*(xne1)/(xne1+xne2)+22*(1-xc-xo)*(xne2)/(xne1+xne2))
        Xnew2 = 16*xo/(12*xc + 16*xo + 20*(1-xc-xo)*(xne1)/(xne1+xne2)+22*(1-xc-xo)*(xne2)/(xne1+xne2))

        Xnew3_1 = (20*(1-xc-xo)*(xne1)/(xne1+xne2))/(12*xc + 16*xo + 20*(1-xc-xo)*(xne1)/(xne1+xne2)+22*(1-xc-xo)*(xne2)/(xne1+xne2))
        Xnew3_2 = (22*(1-xc-xo)*(xne2)/(xne1+xne2))/(12*xc + 16*xo + 20*(1-xc-xo)*(xne1)/(xne1+xne2)+22*(1-xc-xo)*(xne2)/(xne1+xne2))
         
         Dd=[0,0,0,0]
         Dd(1)= Xnew1 - X1
         Dd(2)= Xnew2 - X2
         Dd(3)= Xnew3_1 - X3_1
         Dd(4)= Xnew3_2 - X3_2

         !write(*,*) 'delta_XC: ',Dd(1),' delta_XO: ', Dd(2), 'delta_XNe:', Dd(3)+Dd(4) 

      end subroutine medin_cumming_3p_d_cone

      subroutine medin_cumming_3p_d_neomg(X1,X2,X3_1,X3_2,Dd)
        real(dp), intent(in) :: X1, X2, X3_1, X3_2 ! mass fraction
        real(dp), dimension(4),intent(out) :: Dd
        real(dp) :: Xnew1, Xnew2, Xnew3_1, Xnew3_2, Xfac ! mass fraction
        real(dp) :: xmg, dxmg, xo, dxo, xne1, xne2 ! number fractions
        real(dp) :: dx1_,dx2_
        integer :: i,j

        Xfac = X1 + X2 + X3_1 + X3_2

        xmg = (X1/24)/(X1/24 + X2/16 + X3_1/20 + X3_2/22)
        xo = (X2/16)/(X1/24 + X2/16 + X3_1/20 + X3_2/22)
        
        xne1 = (X3_1/20)/(X1/24 + X2/16 + X3_1/20 + X3_2/22)
        xne2 = (X3_2/22)/(X1/24 + X2/16 + X3_1/20 + X3_2/22)

        call tab_interp_medin_cumming_dx1(xmg,xo,'NeOMg',dx1_)
        call tab_interp_medin_cumming_dx2(xmg,xo,'NeOMg',dx2_)
        dxmg=dx1_
        dxo=dx2_     
        
        xmg = xmg + dxmg
        xo = xo + dxo

        ! convert deltas in number fraction to mass fraction

        Xnew1 = 24*xmg/(24*xmg + 16*xo + 20*(1-xmg-xo)*(xne1)/(xne1+xne2)+22*(1-xmg-xo)*(xne2)/(xne1+xne2))
        Xnew2 = 16*xo/(24*xmg + 16*xo + 20*(1-xmg-xo)*(xne1)/(xne1+xne2)+22*(1-xmg-xo)*(xne2)/(xne1+xne2))

        Xnew3_1 = (20*(1-xmg-xo)*(xne1)/(xne1+xne2))/(24*xmg + 16*xo + 20*(1-xmg-xo)*(xne1)/(xne1+xne2)+22*(1-xmg-xo)*(xne2)/(xne1+xne2))
        Xnew3_2 = (22*(1-xmg-xo)*(xne2)/(xne1+xne2))/(24*xmg + 16*xo + 20*(1-xmg-xo)*(xne1)/(xne1+xne2)+22*(1-xmg-xo)*(xne2)/(xne1+xne2))
         
         Dd=[0,0,0,0]
         Dd(1)= Xnew1 - X1
         Dd(2)= Xnew2 - X2
         Dd(3)= Xnew3_1 - X3_1
         Dd(4)= Xnew3_2 - X3_2

         !write(*,*) 'delta_XMg: ',Dd(1),' delta_XO: ', Dd(2), 'delta_XNe:', Dd(3)+Dd(4) 

      end subroutine medin_cumming_3p_d_neomg

      subroutine medin_cumming_3p_d_onena(X1,X2,X3_1,X3_2,Dd)
        real(dp), intent(in) :: X1, X2, X3_1, X3_2 ! mass fraction
        real(dp), dimension(4),intent(out) :: Dd
        real(dp) :: Xnew1, Xnew2, Xnew3_1, Xnew3_2, Xfac ! mass fraction
        real(dp) :: xna, dxna, xo, dxo, xne1, xne2 ! number fractions
        real(dp) :: dx1_,dx2_
        integer :: i,j

        Xfac = X1 + X2 + X3_1 + X3_2
      
        xna = (X1/23)/(X1/23 + X2/16 + X3_1/20 + X3_2/22)
        xo = (X2/16)/(X1/23 + X2/16 + X3_1/20 + X3_2/22)
        
        xne1 = (X3_1/20)/(X1/23 + X2/16 + X3_1/20 + X3_2/22)
        xne2 = (X3_2/22)/(X1/23 + X2/16 + X3_1/20 + X3_2/22)

        call tab_interp_medin_cumming_dx1(xna,xo,'ONeNa',dx1_)
        call tab_interp_medin_cumming_dx2(xna,xo,'ONeNa',dx2_)
        dxna=dx1_
        dxo=dx2_

        !write(*,*) xna,xo
        !write(*,*) 'delta_xna: ',dxna,' delta_xo: ', dxo      

        xna = xna + dxna
        xo = xo + dxo
        
        ! convert deltas in number fraction to mass fraction

        Xnew1 = 23*xna/(23*xna + 16*xo + 20*(1-xna-xo)*(xne1)/(xne1+xne2)+22*(1-xna-xo)*(xne2)/(xne1+xne2))
        Xnew2 = 16*xo/(23*xna + 16*xo + 20*(1-xna-xo)*(xne1)/(xne1+xne2)+22*(1-xna-xo)*(xne2)/(xne1+xne2))

        Xnew3_1 = (20*(1-xna-xo)*(xne1)/(xne1+xne2))/(23*xna + 16*xo + 20*(1-xna-xo)*(xne1)/(xne1+xne2)+22*(1-xna-xo)*(xne2)/(xne1+xne2))
        Xnew3_2 = (22*(1-xna-xo)*(xne2)/(xne1+xne2))/(23*xna + 16*xo + 20*(1-xna-xo)*(xne1)/(xne1+xne2)+22*(1-xna-xo)*(xne2)/(xne1+xne2))
         
         Dd=[0,0,0,0]
         Dd(1)= Xnew1 - X1
         Dd(2)= Xnew2 - X2
         Dd(3)= Xnew3_1 - X3_1
         Dd(4)= Xnew3_2 - X3_2

         !write(*,*) 'delta_XNa: ',Dd(1),' delta_XO: ', Dd(2), 'delta_XNe:', Dd(3)+Dd(4) 

      end subroutine medin_cumming_3p_d_onena


      subroutine medin_cumming_3p_d_comg(X1,X2,X3,Dd)
        real(dp), intent(in) :: X1, X2, X3 ! mass fraction
        real(dp), dimension(4),intent(out) :: Dd
        real(dp) :: Xnew1, Xnew2, Xfac ! mass fraction
        real(dp) :: xc, dxc, xmg, dxmg, xo ! number fractions
        real(dp) :: dx1_,dx2_
        integer :: i,j

        Xfac = X1 + X2 + X3
 
        xc = (X1/12)/(X1/12 + X2/24 + X3/16)
        xmg = (X2/24)/(X1/12 + X2/24 + X3/16)
        
        xo = (X3/16)/(X1/12 + X2/24 + X3/16)

        call tab_interp_medin_cumming_dx1(xc,xmg,'COMg',dx1_)
        call tab_interp_medin_cumming_dx2(xc,xmg,'COMg',dx2_)
        dxc=dx1_
        dxmg=dx2_    

        xc = xc + dxc
        xmg = xmg + dxmg
        
        ! convert deltas in number fraction to mass fraction

        Xnew1 = 12*xc/(12*xc + 24*xmg + 16*(1-xc-xmg))
        Xnew2 = 24*xmg/(12*xc + 24*xmg + 16*(1-xc-xmg))
         
         Dd=[0,0,0,0]
         Dd(1)= Xnew1 - X1
         Dd(2)= Xnew2 - X2

      end subroutine medin_cumming_3p_d_comg

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
             skip_eos=.TRUE., skip_net=.TRUE., skip_neu=.TRUE., skip_kap=.TRUE., ierr=ierr)
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
