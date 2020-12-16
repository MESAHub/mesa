! ***********************************************************************
!
!   Copyright (C) 2013-2019  Bill Paxton, Pablo Marchant & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module binary_tides

      use star_lib
      use star_def
      use const_def
      use utils_lib
      use math_lib
      use binary_def
      
      implicit none


      contains
      
      
      subroutine sync_spin_orbit_torque(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: osep ! orbital separation (cm)
         real(dp) :: qratio ! mass_other_star/mass_this_star
         real(dp) :: rlr ! roche lobe radius (cm)
         real(dp) :: dt_next ! next timestep
         real(dp) :: Ftid  ! efficiency of tidal synchronization. (time scale × FSYNC). 
         character (len=strlen) :: sync_type
         character (len=strlen) :: sync_mode 
         type (binary_info), pointer :: b
         ierr = 0

         call star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_ptr'
            return
         end if

         call binary_ptr(s% binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         if (.not. b% do_tidal_sync) return

         call star_ptr(id, s, ierr)

         osep = b% separation
         qratio = b% m(b% a_i) / b% m(b% d_i)
         if (is_donor(b, s)) then
            rlr = b% rl(b% d_i)
         else
            qratio = 1d0/qratio
            rlr = b% rl(b% a_i)
         end if
         dt_next = s% dt

         if (b% point_mass_i /= 1 .and. id == b% s1% id) then
            Ftid = b% Ftid_1
            sync_type = b% sync_type_1
            sync_mode = b% sync_mode_1
         else
            Ftid = b% Ftid_2
            sync_type = b% sync_type_2
            sync_mode = b% sync_mode_2
         end if

         if (b% use_other_sync_spin_to_orbit) then
            call b% other_sync_spin_to_orbit(s% id, s% nz, osep, qratio, rlr, dt_next, Ftid, sync_type, sync_mode, ierr)
         else
            call sync_spin_to_orbit(s% id, s% nz, osep, qratio, rlr, dt_next, Ftid, sync_type, sync_mode, ierr)
         end if
         
      end subroutine sync_spin_orbit_torque
      
      subroutine sync_spin_to_orbit(id, nz, osep, qratio, rl, dt_next, Ftid, sync_type, sync_mode, ierr)
         ! initially based on spiba.f kindly provided by Norbert Langer and group.
         integer, intent(in) :: id
         integer, intent(in) :: nz
         real(dp), intent(in) :: osep ! orbital separation (cm)
         real(dp), intent(in) :: qratio ! mass_other_star/mass_this_star
         real(dp), intent(in) :: rl ! roche lobe radius (cm)
         real(dp), intent(in) :: dt_next ! next timestep
         real(dp), intent(in) :: Ftid ! efficiency of tidal synchronization. (time scale / Ftid ). 
         
         character (len=strlen), intent(in) :: sync_type ! synchronization timescale
         character (len=strlen), intent(in) :: sync_mode ! where to put/take angular momentum
         integer, intent(out) :: ierr
      
         type (star_info), pointer :: s
         real(dp) :: G, m, m_lim, t_sync, r_phot, delta_total_J, &
            sum_J_sync, sum_J_non_sync, tdyn, tkh, rho_face, cv_face, &
            T_face, csound_face, ff, omega_orb
            
         real(dp), dimension(nz) :: j_sync, delta_j, tdyn_div_tkh         
         integer, dimension(nz) :: layers_in_sync
         integer :: k, k_rl, k_xm, num_sync_layers
         type (binary_info), pointer :: b

         real(dp) :: a1,a2
      
         include 'formats'
      
         ierr = 0

         call star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_ptr'
            return
         end if

         if (osep <= 0) return
         if (qratio <= 0) return
         if (rl <= 0) return
         if (dt_next <= 0) return
         if (Ftid <= 0) return
         if (sync_type == "None") return

         call binary_ptr(s% binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         
         t_sync = 0
      
         G = standard_cgrav

         if (is_donor(b, s)) then
            m = b% m(b% d_i)
            r_phot = b% r(b% d_i)
         else
            m = b% m(b% a_i)
            r_phot = b% r(b% a_i)
         end if
         
         omega_orb = 2d0*pi/b% period
         do k=1,nz
            j_sync(k) = omega_orb*s% i_rot(k)
         end do
      
         if (sync_type == "Instantaneous") then ! instantaneous synchronisation
            do k=1,nz
               delta_j(k) = s% j_rot(k) - j_sync(k)
            end do
         else
            ! get synchronization timescale, i.e.
            ! 1/t_sync = \dot(omega_star)/(omega_orb-omega_star)
            ! to order e (e=eccentricity).
            if (.not. b% use_other_tsync) then
               call get_tsync(s% id, sync_type, Ftid, qratio, m, r_phot, osep, t_sync, ierr)
               if (ierr/=0) return
            else
               call b% other_tsync(s% id, sync_type, Ftid, qratio, m, r_phot, osep, t_sync, ierr)
               if (ierr/=0) return
            end if
            if (sync_mode == "Uniform") then
               a1 = f2(b% eccentricity)
               a2 = pow(1-pow2(b% eccentricity), 1.5d0)*f5(b% eccentricity)
               do k=1,nz
                  delta_j(k) = (1d0 - exp(-a2*dt_next/t_sync))*(s% j_rot(k) - a1/a2*j_sync(k))
               end do
            else if (sync_mode == "tdyn_div_tkh") then
               ! set local timescale following Gautschy & Glatzel 1990 MNRAS 245, 597-613.
               do k=1,nz
                  !tdyn_div_tkh := local dynamical time-scale / local thermal time-scale
                  rho_face = star_interp_val_to_pt(s% rho,k,nz,s% dq,'binary_tides')
                  cv_face = star_interp_val_to_pt(s% cv,k,nz,s% dq,"binary_tides")
                  T_face = star_interp_val_to_pt(s% T,k,nz,s% dq,"binary_tides")
                  csound_face = star_interp_val_to_pt(s% csound,k,nz,s% dq,"binary_tides")
                  tkh = 4*pi*s% r(k)*s% r(k)*rho_face*cv_face*T_face/s% L(k) ! (4.4)
                  tdyn = 1/csound_face
                  tdyn_div_tkh(k) = tdyn/tkh
               end do

               ! j_sync = synchronized specific angular momentum
               delta_total_J = 0
               do k=1,nz
                  delta_total_J = delta_total_J + abs(s% j_rot(k) - j_sync(k))*s% dm_bar(k)
                  delta_j(k) = 0
                  layers_in_sync(k) = 0
               end do
               delta_total_J = delta_total_J*(1d0 - exp(-dt_next/t_sync))
      
               ! Iteratively solve the scaling factor ff to add (or remove) delta_total_J.
               ! At each iteration, ff is solved such that each zone k has a change on
               ! its angular momentum J_k of the form:
               !
               ! Delta J_k = (J_k-J_{k,sync})*tdyn_div_tkh(k)**2*ff,
               !
               ! and the total change adds up to delta_total_J.
               ! Since tides can at most drive a layer into sync, ff cannot be
               ! solved right away. Then, each iteration solves ff ignoring layers
               ! going over sync, spreads the angular momentum to see which ones
               ! become synced, and then recalculates taking this into account
               ! until it converges.
               num_sync_layers = 0
               !write(*,*) "Entering tide iteration", delta_total_J
               do while (.true.)
                  sum_J_sync = 0
                  sum_J_non_sync = 0
                  do k=1,nz
                     delta_j(k) = s% j_rot(k) - j_sync(k)
                     if (layers_in_sync(k) /= 1) then
                        sum_J_non_sync = sum_J_non_sync + &
                            s% dm_bar(k) * abs(delta_j(k)) * pow2(tdyn_div_tkh(k))
                     else
                        sum_J_sync = sum_J_sync + s% dm_bar(k) * abs(delta_j(k))
                     end if
                  end do

                  ff = (delta_total_J - sum_J_sync) / sum_J_non_sync

                  do k=1,nz
                     if (layers_in_sync(k) == 1) cycle
                     if (delta_j(k) >= 0) then
                        delta_j(k) = min(delta_j(k),delta_j(k)*pow2(tdyn_div_tkh(k))*ff)
                     else
                        delta_j(k) = max(delta_j(k),delta_j(k)*pow2(tdyn_div_tkh(k))*ff)
                     end if
                     if (delta_j(k) == s% j_rot(k) - j_sync(k)) layers_in_sync(k) = 1
                  end do
                  if (num_sync_layers == sum(layers_in_sync)) exit
                  num_sync_layers = sum(layers_in_sync)

               end do
               !write(*,*) "tide iteration", num_sync_layers, &
               !    delta_total_J, sum_J_sync, sum_J_non_sync, ff
            else
               write(*,*) "sync_mode = " , sync_mode , " not recognized"
               write(*,*) "not doing tides"
               do k=1,nz
                  delta_j(k) = 0d0
               end do
            end if
      
         end if

         if (b% point_mass_i /= 1 .and. b% s1% id == s% id) then
            b% t_sync_1 = t_sync
            if (b% model_twins_flag) then
               b% t_sync_2 = t_sync
            end if
         else
            b% t_sync_2 = t_sync
         end if

         if (.not. b% doing_first_model_of_run) then
            do k=1,nz
               s% extra_jdot(k) = s% extra_jdot(k) - delta_j(k)/dt_next
            end do
         end if
      
      end subroutine sync_spin_to_orbit


      real(dp) function f2(e)
         real(dp), intent(in) :: e
         
         f2 = 1d0

         ! Hut 1981, A&A, 99, 126, definition of f2 after eq. 11
         if (e > 0d0) then
             f2 = 1d0 + 15d0/2d0*pow2(e) + 45d0/8d0*pow4(e) + 5d0/16d0*pow6(e)
         end if
          
      end function f2
    
      real(dp) function f3(e)
         real(dp), intent(in) :: e
         
         f3 = 1d0

         ! Hut 1981, A&A, 99, 126, definition of f3 after eq. 11
         if (e > 0d0) then
             f3 = 1d0 + 15d0/4d0*pow2(e) + 15d0/8d0*pow4(e) + 5d0/64d0*pow6(e)
         end if
          
      end function f3
      
      
      real(dp) function f4(e)
         real(dp), intent(in) :: e
         
         f4 = 1d0

         ! Hut 1981, A&A, 99, 126, definition of f4 after eq. 11
         if (e > 0d0) then
             f4 = 1d0 + 3d0/2d0*pow2(e) + 1d0/8d0*pow4(e)
         end if
          
      end function f4
      

      real(dp) function f5(e)
         real(dp), intent(in) :: e
         
         f5 = 1d0

         ! Hut 1981, A&A, 99, 126, definition of f5 after eq. 11
         if (e > 0d0) then
             f5 = 1d0 + 3d0*pow2(e) + 3d0/8d0*pow4(e)
         end if
          
      end function f5
      
      subroutine get_tsync(id, sync_type, Ftid, qratio, m, r_phot, osep, t_sync, ierr)
         integer, intent(in) :: id
         character (len=strlen), intent(in) :: sync_type ! synchronization timescale
         real(dp), intent(in) :: Ftid ! efficiency of tidal synchronization. (time scale / Ftid ). 
         real(dp), intent(in) :: qratio ! mass_other_star/mass_this_star
         real(dp), intent(in) :: m
         real(dp), intent(in) :: r_phot
         real(dp), intent(in) :: osep ! orbital separation (cm)
         real(dp), intent(out) :: t_sync
         integer, intent(out) :: ierr
         real(dp) :: rGyr_squared, moment_of_inertia
         type (binary_info), pointer :: b
         type (star_info), pointer :: s
      
         include 'formats'
   
         ierr = 0

         call star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_ptr'
            return
         end if

         call binary_ptr(s% binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         ! calculate the gyration radius squared
         moment_of_inertia = dot_product(s% i_rot(:s% nz), s% dm_bar(:s% nz))
         rGyr_squared = (moment_of_inertia/(m*r_phot*r_phot))
         if (sync_type == "Hut_conv") then
            ! eq. (11) of Hut, P. 1981, A&A, 99, 126
            t_sync = 3d0*k_div_T(b, s, .true.)*(qratio*qratio/rGyr_squared)*pow6(r_phot/osep)
            ! invert it.
            t_sync = 1d0/t_sync
         else if (sync_type == "Hut_rad") then
            ! eq. (11) of Hut, P. 1981, A&A, 99, 126
            t_sync = 3d0*k_div_T(b, s,.false.)*(qratio*qratio/rGyr_squared)*pow6(r_phot/osep)
            ! invert it.
            t_sync = 1d0/t_sync
         else if (sync_type == "Orb_period") then ! sync on timescale of orbital period
            t_sync = b% period ! synchronize on timescale of orbital period
         else
            ierr = -1
            write(*,*) 'unrecognized sync_type', sync_type
            return
         end if
         t_sync = t_sync / Ftid
      end subroutine get_tsync

      real(dp) function k_div_T(b, s, has_convective_envelope)
         type(binary_info), pointer :: b
         type(star_info), pointer :: s
         logical, intent(in) :: has_convective_envelope

         integer :: k
         real(dp) osep, qratio, m, r_phot,porb, m_env, r_env, tau_conv, P_tid, f_conv

         ! k/T computed as in Hurley, J., Tout, C., Pols, O. 2002, MNRAS, 329, 897
         ! Kudos to Francesca Valsecchi for help implementing and testing this

         k_div_T = 0d0

         osep = b% separation
         qratio = b% m(b% a_i) / b% m(b% d_i)
         if (is_donor(b, s)) then
            m = b% m(b% d_i)
            r_phot = b% r(b% d_i)
         else
            qratio = 1d0/qratio
            m = b% m(b% a_i)
            r_phot = b% r(b% a_i)
         end if
         porb = b% period

         if (has_convective_envelope) then
            m_env = 0d0
            r_env = 0d0
            do k=1, s% nz
               if (s% mixing_type(k) /= convective_mixing .and. &
                   s% rho(k) > 1d5*s% rho(1)) then
                  r_env = (r_phot - s% r(k))/Rsun
                  m_env = (s% m(1) - s% m(k))/Msun
                  exit
               end if
            end do
            tau_conv = 0.431d0*pow(m_env*r_env* &
               (r_phot/Rsun-r_env/2d0)/3d0/s% L_phot,one_third) * secyer
            if (1d0/porb /= s% omega_avg_surf/(2d0*pi)) then
               P_tid = 1d0/abs(1d0/porb-s% omega_avg_surf/(2d0*pi))
               f_conv = min(1.0d0, pow(P_tid/(2d0*tau_conv),b% tidal_reduction))
            else
               f_conv = 1d0
            end if

            k_div_T = 2d0/21d0*f_conv/tau_conv*m_env/(m/Msun)
         else
            !NOTE:There is a typo in eq. (42) of Hurley+ 2002,
            !correct expression is given in footnote 3 of
            !Sepinsky+ 2007
            k_div_T = 1.9782d4*sqrt(m*r_phot*r_phot/pow5(osep)/(Msun/pow3(Rsun)))
            k_div_T = k_div_T*pow(1d0+qratio,5d0/6d0)
            k_div_T = k_div_T*1.592d-9*pow(m/Msun,2.84d0)/secyer
         end if
          
      end function k_div_T

      end module binary_tides

