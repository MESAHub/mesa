! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use chem_def
      
      implicit none

      ! s% xtra(x_old_wind) every step we save the used wind mass loss rate here,
      ! this is used to soften too large changes in wind mass loss rates 
      integer, parameter :: x_old_wind = 1

      ! s% xtra(x_time_thermal_eq) holds time the star has been in thermal equilibrium
      ! used to determine the start of ZAMS
      integer, parameter :: x_time_thermal_eq = 2

      ! s% lxtra(lx_pre_ZAMS) true if the star has not already properly settled in the ZAMS
      integer, parameter :: lx_pre_ZAMS = 1
	  
      ! s% xtra(x_old_wind) every step we save the used wind mass loss rate here,
      ! this is used to soften too large changes in wind mass loss rates 
      !integer, parameter :: x_old_wind = 1
      
      ! these routines are called by the standard run_star check_model
      contains
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)


         ! the extras functions in this file will not be called
         ! unless you set their function pointers as done below.
         ! otherwise we use a null_ version which does nothing (except warn).

         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

         s% other_wind => erupt_other_wind

      end subroutine extras_controls
      
      subroutine get_conv_regions_above_T(id, T_limit, ierr, num_conv_regions)
         ! use mlt_def, only: convective_mixing
         type (star_info), pointer :: s
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         real(dp) :: T_limit
         integer :: prev_type, cur_type, cur_top, n, k, num_conv_regions, max_num_conv_regions, n_limit
         include 'formats'

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ierr = 0
         cur_type = s% mixing_type(1)
         cur_top = 1
         n = 0
         n_limit = 0
         max_num_conv_regions = 4 ! Max number of convective regions allowed


         ! Find gridpoint corresponding to max temperature (select only outer layers)
         do k = 1, s% nz
            if (s% T(k) < T_limit) then
                  n_limit = k
            end if
         end do

         ! Find all convective regions in the outer layers down to T_limit
         do k = 2, n_limit
            prev_type = cur_type
            cur_type = s% mixing_type(k)
            if (cur_type == prev_type .and. k < n_limit) cycle
            ! change of type from k-1 to k
            if (prev_type == convective_mixing) then
               n = n + 1
               s% mixing_region_type(n) = prev_type
               s% mixing_region_top(n) = cur_top
               if (k == n_limit) then
                  s% mixing_region_bottom(n) = k
               else
                  s% mixing_region_bottom(n) = k-1
               end if
               if (n == max_num_conv_regions) exit
            end if
            cur_top = k
         end do

         num_conv_regions = n

      end subroutine get_conv_regions_above_T
	  

      subroutine get_convective_core(id, sc_convective_core,ierr)
           type (star_info), pointer :: s
           integer, intent(in) :: id
           integer, intent(out) :: ierr, sc_convective_core
           integer :: cur_type, k
           include 'formats'

           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return

           k = s% nz
           cur_type = s% mixing_type(k)
           !write(*,*) 'Convective type', cur_type
           do while (cur_type == 1 .and. k > 2)
             cur_type = s% mixing_type(k)
             k = k - 1
           end do
           sc_convective_core = k
      end subroutine get_convective_core

  	  subroutine classify_conv_region_above_T(id, ierr, sc_top, sc_bottom, sc_type)

           type (star_info), pointer :: s
           integer, intent(in) :: id
           integer, intent(out) :: ierr

           character (len=7) ::  sc_type
           real(dp), DIMENSION(2) :: T_HI, T_HeI, T_HeII, T_FeCZ
           integer :: n, k, sc_top, sc_bottom
           include 'formats'

           ! Pass upper and lower gridpoints of convective regions, check temperature and classify

           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return


           T_HI    = (/ 3000,11000 /)     ! Rough T range for H Conv Region
           T_HeI   = (/ 11000,35000 /)    ! Rough T range for HeI Conv Region
           T_HeII  = (/ 35000,100000 /)   ! Rough T range for HeII Conv Region
           T_FeCZ  = (/ 100000,500000 /)  ! Rough T range for FeCZ Conv Region

           !write(*,*)   T_HI(1), T_HI(2), MAXVAL(s% T(sc_top:sc_bottom))

           ! Find gridpoint corresponding to max temperature (select only outer layers)

           sc_type = 'UNKNOWN'
           if ( sc_top > 0 ) then
             if (s% T(sc_top) < T_HI(2)) then
               sc_type = 'HI'
             else
               do k = sc_top, sc_bottom
                  if  (s% T(k) > T_HeI(1) .AND. s% T(k) < T_HeI(2)) then
                    sc_type = 'HeI'
              	else if (s% T(k) > T_HeII(1) .AND. s% T(k) < T_HeII(2)) then
                    sc_type = 'HeII'
              	else if (s% T(k) > T_FeCZ(1) .AND. s% T(k) < T_FeCZ(2)) then
                    sc_type = 'FeCZ'
              	else
                    sc_type = 'UNKNOWN'
              	end if
           	 end do
             end if
           end if
           write(*,*) 'Type: ', s% T(sc_top), s% T(sc_bottom), sc_type
        end subroutine classify_conv_region_above_T
		
		
    	  subroutine get_conv_radii(id, ierr, sc_top, sc_bottom, r_top, r_bottom)
             type (star_info), pointer :: s
             integer, intent(in) :: id
             integer, intent(out) :: ierr
             real(dp) :: r_top, r_bottom
             integer :: n, k, sc_top, sc_bottom
             include 'formats'

             ierr = 0
             call star_ptr(id, s, ierr)
             if (ierr /= 0) return
             r_top = 0d0
             r_bottom = 0d0

             if ( sc_top > 0 ) then
             	r_top = s% r(sc_top)
             	r_bottom = s% r(sc_bottom)
             end if
          end subroutine get_conv_radii
		  
		  

		  
		  
		  

          subroutine get_hp_radii(id, ierr, hp, k_hp)
             type (star_info), pointer :: s
             integer, intent(in) :: id
             integer, intent(out) :: ierr
             integer :: k_hp
             real(dp) :: hp
             include 'formats'

             ierr = 0
             call star_ptr(id, s, ierr)
             if (ierr /= 0) return
             k_hp = 1
    !        write(*,*) s% lnP
    !         do while (((s% lnP(k_hp) - s% lnP(1)) < LOG10(hp)) .and. (k_hp < s% nz))
             do while ((LOG(EXP(s% lnPgas(k_hp))/EXP(s% lnPgas(1))) < hp).and.(k_hp < s% nz))
             !           write(*,*)  s% lnP(k_hp),  s% lnP(1), LOG10(hp)
             	k_hp = k_hp + 1
             end do
             !write(*,*) 'Hp: ',hp,'P/P0: ', EXP(s% lnP(1))/EXP(s% lnP(k_hp))
          end subroutine get_hp_radii
		  
		  
          subroutine get_max_fc(id, ierr, fcmax, sc_top, sc_bottom)
             type (star_info), pointer :: s
             integer, intent(in) :: id
             integer, intent(out) :: ierr
             real(dp) :: fcmax
             integer :: n, k, sc_top, sc_bottom
             include 'formats'

             ierr = 0
             call star_ptr(id, s, ierr)
             if (ierr /= 0) return
             fcmax = 0
             ! fcmax = MAXVAL((s% conv_vel(sc_top:sc_bottom)**3.0) * s% rho(sc_top:sc_bottom))
             fcmax = MAXVAL(s% L_conv(sc_top:sc_bottom)/s% L(sc_top:sc_bottom)) ! Max of Lconv / Lstar
          end subroutine get_max_fc
		  
          subroutine get_max_p_turb_over_p(id, ierr, p_turb_over_p, sc_top, sc_bottom)
             type (star_info), pointer :: s
             integer, intent(in) :: id
             integer, intent(out) :: ierr
             real(dp) :: p_turb_over_p, p_turb
             integer :: n, k, sc_top, sc_bottom
             include 'formats'

             ierr = 0
             call star_ptr(id, s, ierr)
             if (ierr /= 0) return
             p_turb_over_p = 0d0
             p_turb = 0d0
             ! Calculate Pturb = (1d0/3d0) * rho vc^2 ! Follow Eq (1) of Grassitelli et al. 2015
             ! Calculate max Pturb/P
             ! p_turb =  (1d0/3d0) * s%rho(sc_top:sc_bottom) * pow(s% conv_vel(sc_top:sc_bottom),2)
             p_turb_over_p =  MAXVAL((1d0/3d0) * s%rho(sc_top:sc_bottom) * pow(s% conv_vel(sc_top:sc_bottom),2)/EXP(s% lnPgas(sc_top:sc_bottom)))
             !write(*,*) 'p_turb_over_p: ',p_turb_over_p
          end subroutine get_max_p_turb_over_p
		  
		  
          subroutine get_max_p_turb_over_ptopHp(id, ierr, p_turb_over_ptopHp, sc_top,sc_bottom)
             type (star_info), pointer :: s
             integer, intent(in) :: id
             integer, intent(out) :: ierr
			 integer :: Hpidxarr(1)
			 integer :: Hpidx
             real(dp) :: p_turb_over_ptopHp, p_turb
             integer :: n, k, sc_top, sc_bottom
             include 'formats'

             ierr = 0
             call star_ptr(id, s, ierr)
             if (ierr /= 0) return
             p_turb_over_ptopHp = 0d0
             p_turb = 0d0
			 Hpidxarr = minloc(abs(s% Pgas(sc_top:s% nz)-(s% Pgas(sc_top)*2.71828d0)))
			 Hpidx = Hpidxarr(1)+sc_top
			 !write(*,*) 'Pgas at top', s% Pgas(sc_top)
			 !write(*,*) 'Pgas at top div e', s% Pgas(sc_top)*2.71828d0
			 
			 !write(*,*) 'Hpidxarr', Hpidxarr
			 !write(*,*) 'Hpidx', Hpidx
			 !write(*,*) 'sc_top', sc_top
			 !write(*,*) 'Pgasfull', s% Pgas(sc_top:Hpidx)
			 
             ! Calculate Pturb = (1d0/3d0) * rho vc^2 ! Follow Eq (1) of Grassitelli et al. 2015
             ! Calculate max Pturb/P
             ! p_turb =  (1d0/3d0) * s%rho(sc_top:sc_bottom) * pow(s% conv_vel(sc_top:sc_bottom),2)
             p_turb_over_ptopHp =  MAXVAL((1d0/3d0) * s%rho(sc_top:Hpidx) * pow(s% conv_vel(sc_top:Hpidx),2)/EXP(s% lnPgas(sc_top:Hpidx)))
             !write(*,*) 'p_turb_over_p: ',p_turb_over_p
          end subroutine get_max_p_turb_over_ptopHp
	  
	  
          subroutine get_conv_mass(id, ierr, sc_top, sc_bottom, total_mass)
          type (star_info), pointer :: s
             integer, intent(in) :: id
             integer, intent(out) :: ierr
             real(dp) :: total_mass
             integer :: n, k, sc_top, sc_bottom
             include 'formats'

             ierr = 0
             call star_ptr(id, s, ierr)
             if (ierr /= 0) return

             total_mass = 0d0
             if ( sc_top > 0 ) then
             	total_mass = SUM(s% dm(sc_top:sc_bottom))
             end if

          end subroutine  get_conv_mass
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
          subroutine get_conv_velocities(id, ierr, v_max, v_aver, sc_top, sc_bottom,b_eq,b_max,rho_aver)
             type (star_info), pointer :: s
             integer, intent(in) :: id
             integer, intent(out) :: ierr
             real(dp) :: v_max, v_aver, b_eq, b_max, rho_aver, rho_v_max, pi
             integer :: n, k, sc_top, sc_bottom, i_v_max
             include 'formats'

             ierr = 0
             call star_ptr(id, s, ierr)
             if (ierr /= 0) return

             pi = 3.1415
             v_max = 0d0
             v_aver = 0d0
             rho_aver = 0d0
             b_eq = 0d0
             b_max = 0d0
             rho_aver = 0d0
             i_v_max = 0
             rho_v_max = 0d0
             ! Calculate Max and average V_conv in Conv region (average in dr, see Eq. 6 in Cantiello et al. 2009)


             if ( sc_top > 0 ) then
             	do k = sc_top, sc_bottom
                	v_aver = v_aver + (s% r(k) - s% r(k+1)) * s% conv_vel(k)
                	rho_aver = rho_aver + (s% r(k) - s% r(k+1)) * s% rho(k)
                !  write(*,*) 'DRHO: ', (s% r(k) - s% r(k+1)) * s% rho(k), s% rho(k)
              !    write(*,*) 'DV: ',(s% r(k) - s% r(k+1))*s% conv_vel(k),s% conv_vel(k)
             	end do
             	v_max = MAXVAL(s% conv_vel(sc_top:sc_bottom))
              i_v_max = MAXLOC (s% conv_vel(sc_top:sc_bottom), DIM=1)
              rho_v_max = s% rho(i_v_max)
             	v_aver = v_aver /( s% r(sc_top) - s% r(sc_bottom) )
             	rho_aver = rho_aver /( s% r(sc_top) - s% r(sc_bottom) )


             end if
             ! Calculate B_equipartition and B_max

             b_eq = (v_aver)*(4.0*pi*rho_aver)**(0.5)
             b_max = (v_max)*(4.0*pi*rho_v_max)**(0.5)
             !b_max = (v_max)*(4.0*pi*rho_aver)**(0.5) ! For convective core this would work better

             !write(*,*) v_aver, v_max, rho_aver, b_eq, b_max

          end subroutine get_conv_velocities
		  
		  
		  
          subroutine get_conv_velocitiestopHp(id, ierr, v_max, v_aver, sc_top, sc_bottom,b_eq,b_max,rho_aver)
             type (star_info), pointer :: s
             integer, intent(in) :: id
             integer, intent(out) :: ierr
			 integer :: Hpidxarr(1)
			 integer :: Hpidx
             real(dp) :: v_max, v_aver, b_eq, b_max, rho_aver, rho_v_max, pi
             integer :: n, k, sc_top, sc_bottom, i_v_max
             include 'formats'

             ierr = 0
             call star_ptr(id, s, ierr)
             if (ierr /= 0) return

             pi = 3.1415
             v_max = 0d0
             v_aver = 0d0
             rho_aver = 0d0
             b_eq = 0d0
             b_max = 0d0
             rho_aver = 0d0
             i_v_max = 0
             rho_v_max = 0d0
             ! Calculate Max and average V_conv in Conv region (average in dr, see Eq. 6 in Cantiello et al. 2009)
			 Hpidxarr = minloc(abs(s% Pgas(sc_top:s% nz)-(s% Pgas(sc_top)*2.71828d0)))
			 Hpidx = Hpidxarr(1)+sc_top

             if ( sc_top > 0 ) then
             	do k = sc_top, Hpidx
                	v_aver = v_aver + (s% r(k) - s% r(k+1)) * s% conv_vel(k)
                	rho_aver = rho_aver + (s% r(k) - s% r(k+1)) * s% rho(k)
                !  write(*,*) 'DRHO: ', (s% r(k) - s% r(k+1)) * s% rho(k), s% rho(k)
              !    write(*,*) 'DV: ',(s% r(k) - s% r(k+1))*s% conv_vel(k),s% conv_vel(k)
             	end do
             	v_max = MAXVAL(s% conv_vel(sc_top:Hpidx))
              i_v_max = MAXLOC (s% conv_vel(sc_top:Hpidx), DIM=1)
              rho_v_max = s% rho(i_v_max)
             	v_aver = v_aver /( s% r(sc_top) - s% r(Hpidx) )
             	rho_aver = rho_aver /( s% r(sc_top) - s% r(Hpidx) )


             end if
             ! Calculate B_equipartition and B_max

             b_eq = (v_aver)*(4.0*pi*rho_aver)**(0.5)
             b_max = (v_max)*(4.0*pi*rho_v_max)**(0.5)
             !b_max = (v_max)*(4.0*pi*rho_aver)**(0.5) ! For convective core this would work better

             !write(*,*) v_aver, v_max, rho_aver, b_eq, b_max

          end subroutine get_conv_velocitiestopHp
		  
		  
		  
          subroutine get_pressure_eq_field(id, ierr, sc_top, sc_bottom,b_p_eq,b_p_max)
             type (star_info), pointer :: s
             integer, intent(in) :: id
             integer, intent(out) :: ierr
             real(dp) :: p_max, p_aver, b_p_eq, b_p_max
             integer :: n, k, sc_top, sc_bottom
             include 'formats'

             ierr = 0
             call star_ptr(id, s, ierr)
             if (ierr /= 0) return


             b_p_eq = 0d0
             b_p_max = 0d0
             p_aver = 0d0

             if ( sc_top > 0 ) then
             	do k = sc_top, sc_bottom
                	p_aver = p_aver + (s% r(k) - s% r(k+1)) * EXP(s% lnPgas(k))
             	end do
             	p_max = EXP(MAXVAL(s% lnPgas(sc_top:sc_bottom)))
             	p_aver = p_aver /( s% r(sc_top) - s% r(sc_bottom) )
             end if
             ! Calculate B_Pressure_equipartition and B_max

             b_p_eq = (p_aver*8.0*pi)**(0.5)
             b_p_max = (p_max*8*pi)**(0.5)

    !         write(*,*) b_p_eq, b_p_max

          end subroutine get_pressure_eq_field
		  
		  

          subroutine get_average_hp(id, ierr, sc_top, sc_bottom, hp_aver)
             type (star_info), pointer :: s
             integer, intent(in) :: id
             integer, intent(out) :: ierr
             real(dp) :: hp_aver
             integer :: n, k, sc_top, sc_bottom
             include 'formats'

             ierr = 0
             call star_ptr(id, s, ierr)
             if (ierr /= 0) return

             hp_aver = 0d0
             ! Calculate average hp in convective region (dr weighted)

             if ( sc_top > 0 ) then
             	do k = sc_top, sc_bottom
                	hp_aver = hp_aver + (s% r(k) - s% r(k+1)) *s% scale_height(k)
             	end do
                hp_aver = hp_aver /( s% r(sc_top) - s% r(sc_bottom) )
             end if

          end subroutine get_average_hp
		  
		  
		  
          subroutine get_microturb(mach1_aver_ahp, rho1_aver, rho_surf,v1_aver_ahp, v1_surf_aver)
             real(dp) :: v1_surf_aver, v1_aver_ahp, mach1_aver_ahp, rho1_aver, rho_surf

             v1_surf_aver = 0d0
             if (rho_surf /= 0) then
    	         v1_surf_aver = v1_aver_ahp * (mach1_aver_ahp * rho1_aver/rho_surf)**(1./2.)
             end if

          end subroutine get_microturb
		  
		  
		  
      	  subroutine get_turnover(mixing_length_alpha,v_HI_aver, HI_hp_aver, turnover_HI)
               real(dp) :: v_HI_aver, HI_hp_aver, turnover_HI, mixing_length_alpha


                turnover_HI = 0d0
                if (v_HI_aver /= 0) then
                  turnover_HI = mixing_length_alpha*HI_hp_aver/v_HI_aver
               endif


            end subroutine get_turnover
			
			
            subroutine get_F0(rho, radius, bruntN2, F0)
                 real(dp) ::  F0, rho, radius, bruntN2, pi
                 pi = 3.1415
                 ! F0 = 0.5d0 * rho * (bruntN2**0.5)/(2d0*pi) * radius**3d0 * (1/turnover)**2d0
                 F0 = 0.5d0 * rho * (bruntN2**0.5)/(2d0*pi) * radius**3d0
            end subroutine get_F0
			
			
            subroutine get_bsurf(rho_surf, rho_HI_aver, b_HI_aver, b_HI_max, b_HI_surf, b_HI_surf_max)
               real(dp) :: rho_surf, rho_HI_aver, b_HI_aver, b_HI_max, b_HI_surf, b_HI_surf_max
               b_HI_surf = 0d0
               b_HI_surf_max = 0d0
                if (rho_HI_aver /= 0) then
                  b_HI_surf = b_HI_aver * (rho_surf/rho_HI_aver)**(2./3.)
                  b_HI_surf_max = b_HI_max * (rho_surf/rho_HI_aver)**(2./3.)
               endif


            end subroutine get_bsurf
			
			
		    pure function integrate(x, y) result(r)
		      !! Calculates the integral of an array y with respect to x using the trapezoid
		      !! approximation. Note that the mesh spacing of x does not have to be uniform.
		      real(dp), intent(in)  :: x(:)         !! Variable x
		      real(dp), intent(in)  :: y(size(x))   !! Function y(x)
		      real(dp)              :: r            !! Integral ∫y(x)·dx
			  
		      ! Integrate using the trapezoidal rule
		      associate(n => size(x))
		        r = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
		      end associate
		    end function

			
            subroutine get_mloss_forwind(id,ierr, testsize,mdot_erupt)
				 type (star_info), pointer :: s
				 integer, intent(in) :: testsize, id
				 integer, intent(out) :: ierr
				 integer :: taucritidxarr(1)
				 integer :: taucritidx
				 integer :: i,  mlost_loc
				 real(dp) :: mlost_final, tdyn_final, masslossres, tdyn_final_alt, masslossres_alt, mdot_erupt
                 real(dp), dimension(testsize) ::  t_dyn_test, L_excess_test, E_excess_test, E_excess_pos, masked, masked2, E_diff, E_excess_final, E_binding_final, E_diff_masked, M_above, mdot_loss, Mloss
				 include 'formats'
				 
                 ierr = 0
                 call star_ptr(id, s, ierr)
                 if (ierr /= 0) return
				 
				 !Find the index at which critical optical depth occurs
				 taucritidxarr = minloc(abs(s% tau(1:s% nz)-(clight/(s% csound(1:s% nz)))))
				 taucritidx = taucritidxarr(1)
				 
				 !Initialize arrays
				 do i  =1,s% nz
				    t_dyn_test(i) = 0d0
					L_excess_test(i) = 0d0
					E_excess_test(i) = 0d0
				 end do
                 
				 !Calculate dynamical time, excess Luminosity, and excess Energy. 
				 t_dyn_test = ((s% r(1:s% nz))**3/(6.674e-8*s% m(1:s% nz)))**(1.0/2.0)
				 L_excess_test = (s% L(1:s% nz)) - s% L_conv(1:s% nz) - (pi4*clight*s% cgrav(1:s% nz)*s% m_grav(1:s% nz)/(s% opacity(1:s% nz)))
				 E_excess_test = L_excess_test*t_dyn_test
				 
				 !Mask out negative excess Energy
				 where (E_excess_test < 0)
					 masked2 = 0
				elsewhere
					masked2 = E_excess_test
				end where
				
				!Mask out where tau < critical tau
				where ((s% tau(1:s% nz)) > (clight/(s% csound(1:s% nz))))
					masked = 0 
				elsewhere 
					masked = masked2
				end where
				
				!Calculate total excess energy from surface to critical tau 
				do i  =1,s% nz-1
					E_excess_final(s% nz+1-i) = -1*integrate(s% m(s% nz-i:s% nz), masked(s% nz-i:s% nz))/ABS(s% m(s% nz+1-i)-s% m(taucritidx)) !mass method
				end do
				
				!Find binding energy, mass above
				do i = 1, s% nz-1
					E_binding_final(i) = SUM(6.674e-8*s% m(1:i)*s% dm(1:i)/s% r(1:i)) !still right way, surface first
					M_above(i) = (s% m(1) - s% m(i))/1.9884098706980504d33
				end do
				
				!Patch up first element of arrays
				E_excess_final(1) = E_excess_final(2)
				E_binding_final(s% nz) = E_binding_final(s% nz-1)
				M_above(s% nz) = (s% m(1) - s% m(s% nz))/1.9884098706980504d33
				
				!Define energy difference between excess (due to local supereddington L) and binding
				E_diff = E_excess_final-E_binding_final
				
				!Define total lost mass. Hack-like filtering for where E_diff < 0.
			 	where (E_diff > 0)
				 	Mloss = M_above
				elsewhere
					Mloss = 0d0
				end where
				
				Mloss(taucritidx) = 0d0
				
				!Find indices
				mlost_loc = MAXLOC(Mloss, DIM=1) !E_diff < 0 was set to 0d0 earlier, so can use max() now. Location at first E_diff > 0 location more shallow than critical tau.
				mlost_final = MAXVAL(Mloss) !Final mass lost
				
				tdyn_final = t_dyn_test(mlost_loc)*3.17098d-8 !years

				masslossres = mlost_final/tdyn_final !Msun/yr
				
			 	mdot_erupt = masslossres
			
			 s% xtra(7) = M_above(taucritidx)
			 s% xtra(8) = mlost_final
			 s% xtra(9) = tdyn_final
			 s% xtra(10) = mlost_loc
			 
			
            end subroutine get_mloss_forwind
			
			
			
	        subroutine erupt_other_wind(id, Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z, wind, ierr)
	           use star_def
	           integer, intent(in) :: id
	           real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z ! surface values (cgs)
	           ! NOTE: surface is outermost cell. not necessarily at photosphere.
	           ! NOTE: don't assume that vars are set at this point.
	           ! so if you want values other than those given as args,
	           ! you should use values from s% xh(:,:) and s% xa(:,:) only.
	           ! rather than things like s% Teff or s% lnT(:) which have not been set yet.
	           real(dp), intent(out) :: wind ! wind in units of Msun/year (value is >= 0)
	           integer, intent(out) :: ierr
			   
			   
	           integer :: h1, he4
	           real(dp) :: w, Xs, Ys, Z_div_Z_solar, Teff_jump, alfa, L1, M1, R1, T1, &
	              vink_wind, nieu_wind, hamann_wind, highT_w, lowT_w, Twindow, T_high, T_low, surface_h1, w1, w2, T_decin
			   real(dp), parameter :: Zsolar = 0.019d0 ! for Vink et al formula
			   
			   type (star_info), pointer :: s
	           ierr = 0
	           call star_ptr(id, s, ierr)
	           if (ierr /= 0) return
			   
	           L1 = Lsurf
	           M1 = Msurf
	           R1 = Rsurf
	           T1 = Tsurf
			   
	           h1 = s% net_iso(ih1)
	           he4 = s% net_iso(ihe4)
	           Xs = s% xa(h1,1)
	           Ys = s% xa(he4,1)
			   surface_h1 = s% xa(h1,1)
			   
               T_high = 11000
               T_low = 10000
               if (s% Dutch_scaling_factor == 0) then
                  wind = 0
               else if (T1 <= T_low) then
                  call eval_lowT_Dutch(wind)
               else if (T1 >= T_high) then
                  call eval_highT_Dutch(wind)
               else ! transition
                  call eval_lowT_Dutch(w1)
                  call eval_highT_Dutch(w2)
                  alfa = (T1 - T_low)/(T_high - T_low)
                  wind = (1-alfa)*w1 + alfa*w2
               end if

			   s% xtra(1) = s% Dutch_scaling_factor*wind

			   w = 0
			   
			   call get_mloss_forwind(id, ierr, s% nz, w)
			   
			   s% xtra(5) = w * s% x_ctrl(18)
			   
			   w = w * s% x_ctrl(18)
			   
			   wind = w + s% Dutch_scaling_factor*wind
			   
			   s% xtra(3) = wind
		   
		       contains
			   
	           subroutine eval_highT_Dutch(w)
	              real(dp), intent(out) :: w
	              include 'formats'
	              if (surface_h1 < 0.4d0) then ! helium rich Wolf-Rayet star: Nugis & Lamers
	                 w = 1d-11 * pow(L1/Lsun,1.29d0) * pow(Y,1.7d0) * sqrt(Z)
	                 !if (dbg) write(*,1) 'Dutch_wind = Nugis & Lamers', log10(wind)
	              else
	                 call eval_Vink_wind(w)
	              end if
	           end subroutine eval_highT_Dutch


	           subroutine eval_lowT_Dutch(w)
	              real(dp), intent(out) :: w
	              include 'formats'
	              if (s% Dutch_wind_lowT_scheme == 'de Jager') then
	                 call eval_de_Jager_wind(w)
	                 !if (dbg) write(*,1) 'Dutch_wind = de Jager', safe_log10(wind), T1, T_low, T_high
	              !else if (s% Dutch_wind_lowT_scheme == 'van Loon') then
	                 !call eval_van_Loon_wind(w)
	                 !if (dbg) write(*,1) 'Dutch_wind = van Loon', safe_log10(wind), T1, T_low, T_high
	              !else if (s% Dutch_wind_lowT_scheme == 'Nieuwenhuijzen') then
	                 !call eval_Nieuwenhuijzen_wind(w)
	                 !if (dbg) write(*,1) 'Dutch_wind = Nieuwenhuijzen', safe_log10(wind), T1, T_low, T_high
				  !else if (s% Dutch_wind_lowT_scheme == 'Decin') then
						 !call eval_decin_wind(w)
	              else
	                 write(*,*) 'unknown value for Dutch_wind_lowT_scheme ' // &
	                    trim(s% Dutch_wind_lowT_scheme)
	                 w = 0
	              end if
	           end subroutine eval_lowT_Dutch
			   
	           subroutine eval_Vink_wind(w)
	              real(dp), intent(inout) :: w
	              real(dp) :: alfa, w1, w2, Teff_jump, logMdot, dT, vinf_div_vesc

	              ! alfa = 1 for hot side, = 0 for cool side
	              if (T1 > 27500d0) then
	                 alfa = 1
	              else if (T1 < 22500d0) then
	                 alfa = 0
	              else ! use Vink et al 2001, eqns 14 and 15 to set "jump" temperature
	                 Teff_jump = 1d3*(61.2d0 + 2.59d0*(-13.636d0 + 0.889d0*log10(Z/Zsolar)))
	                 dT = 100d0
	                 if (T1 > Teff_jump + dT) then
	                    alfa = 1
	                 else if (T1 < Teff_jump - dT) then
	                    alfa = 0
	                 else
	                    alfa = 0.5d0*(T1 - (Teff_jump - dT)) / dT
	                 end if
	              end if

	              if (alfa > 0) then ! eval hot side wind (eqn 24)
	                 vinf_div_vesc = 2.6d0 ! this is the hot side galactic value
	                 vinf_div_vesc = vinf_div_vesc*pow(Z/Zsolar,0.13d0) ! corrected for Z
	                 logMdot = &
	                    - 6.697d0 &
	                    + 2.194d0*log10(L1/Lsun/1d5) &
	                    - 1.313d0*log10(M1/Msun/30d0) &
	                    - 1.226d0*log10(vinf_div_vesc/2d0) &
	                    + 0.933d0*log10(T1/4d4) &
	                    - 10.92d0*pow2(log10(T1/4d4)) &
	                    + 0.85d0*log10(Z/Zsolar)
	                 w1 = exp10(logMdot)
	              else
	                 w1 = 0
	              end if

	              if (alfa < 1) then ! eval cool side wind (eqn 25)
	                 vinf_div_vesc = 1.3d0 ! this is the cool side galactic value
	                 vinf_div_vesc = vinf_div_vesc*pow(Z/Zsolar,0.13d0) ! corrected for Z
	                 logMdot = &
	                    - 6.688d0 &
	                    + 2.210d0*log10(L1/Lsun/1d5) &
	                    - 1.339d0*log10(M1/Msun/30d0) &
	                    - 1.601d0*log10(vinf_div_vesc/2d0) &
	                    + 1.07d0*log10(T1/2d4) &
	                    + 0.85d0*log10(Z/Zsolar)
	                 w2 = exp10(logMdot)
	              else
	                 w2 = 0
	              end if

	              w = alfa*w1 + (1 - alfa)*w2

	              !if (dbg) write(*,*) 'vink wind', w

	           end subroutine eval_Vink_wind
			   

	           subroutine eval_de_Jager_wind(w)
	              ! de Jager, C., Nieuwenhuijzen, H., & van der Hucht, K. A. 1988, A&AS, 72, 259.
	              real(dp), intent(out) :: w
	              real(dp) :: log10w
	              include 'formats'
	              log10w = 1.769d0*log10(L1/Lsun) - 1.676d0*log10(T1) - 8.158d0
	              w = exp10(log10w)
	              !if (dbg) then
	                 !write(*,1) 'de_Jager log10 wind', log10w
	              !end if
	           end subroutine eval_de_Jager_wind
			   
	           subroutine eval_decin_wind(w)
	              ! Decin, L. et al. 2024, A&A
	              real(dp), intent(out) :: w
	              real(dp) :: log10w, term1, term2, term3
	              include 'formats'
				  term1 = 1.77d0+0.49d0
				  term2 = -1.68d0 * (s% initial_mass/10d0)
				  term3 = 3.5d0*log10(L1/(1d5*Lsun))
	              
				  log10w = term1+term2+term3
				  
				  w = exp10(log10w)*(1d-5)
				  ! 1.769d0*log10(L1/Lsun) - 1.676d0*log10(T1) - 8.158d0
	           end subroutine eval_decin_wind
			   
			   
	        end subroutine erupt_other_wind
			
			

        	  subroutine get_conv_ahp(id, ierr, sc_top, sc_bottom, v_aver_ahp, mach_top_cz, mach_aver_ahp, rho_aver_ahp)
                 type (star_info), pointer :: s
                 integer, intent(in) :: id
                 integer, intent(out) :: ierr
                 real(dp) ::  v_aver_ahp, mach_aver_ahp, rho_aver_ahp, cs_aver_ahp, delta_r, mach_top_cz
                 integer :: n, k, sc_top, sc_bottom, kk
                 include 'formats'

                 ierr = 0
                 call star_ptr(id, s, ierr)
                 if (ierr /= 0) return

                 v_aver_ahp = 0d0
                 rho_aver_ahp = 0d0
                 cs_aver_ahp = 0d0
                 mach_aver_ahp = 0d0
                 mach_top_cz = 0d0
                 kk = 0

                 ! Calculate rho_aver and v_aver in top alpha*Hp of convective zone ! Follows Eq.6 in Cantiello et al. 2009

                 if ( sc_top > 0 ) then
                 	do k = sc_top, sc_bottom
                    	if (s% r(k) > s% r(sc_top) - s% mixing_length_alpha * s% scale_height(sc_top)) then
                       		rho_aver_ahp = rho_aver_ahp + (s% r(k) - s% r(k+1)) * s% rho(k)
                       		v_aver_ahp =  v_aver_ahp + (s% r(k) - s% r(k+1)) * s%  conv_vel(k)
                       		cs_aver_ahp = cs_aver_ahp + (s% r(k) - s% r(k+1)) * s% csound(k)
                       		kk = k
                    	end if
                 	end do
                 end if
                 rho_aver_ahp = rho_aver_ahp / ( s% r(sc_top) - s% r(kk) )
                 v_aver_ahp = v_aver_ahp / ( s% r(sc_top) - s% r(kk) )
                 cs_aver_ahp = cs_aver_ahp / ( s% r(sc_top) - s% r(kk) )

                 if (cs_aver_ahp /=0) then
                 	mach_aver_ahp = v_aver_ahp/cs_aver_ahp
                 end if

                 if (s% csound(sc_top) /=0) then
                 	mach_top_cz = s%  conv_vel(sc_top) / s% csound(sc_top)
                 end if


              end subroutine get_conv_ahp
	  
	  
              subroutine get_conv_aver(id, ierr, sc_top, sc_bottom, mach_max, mach_aver)
                   type (star_info), pointer :: s
                   integer, intent(in) :: id
                   integer, intent(out) :: ierr
                   real(dp) ::  v_aver, mach_aver, rho_aver, cs_aver, delta_r, mach_max
                   integer :: n, k, sc_top, sc_bottom, kk
                   include 'formats'

                   ierr = 0
                   call star_ptr(id, s, ierr)
                   if (ierr /= 0) return

                   v_aver = 0d0
                   rho_aver = 0d0
                   cs_aver = 0d0
                   mach_aver = 0d0
                   mach_max = 0d0
                   kk = 0

                   ! Calculate rho_aver and v_aver in the full convective zone ! Radial average

                   if ( sc_top > 0 ) then
                   	do k = sc_top, sc_bottom
                        rho_aver = rho_aver + (s% r(k) - s% r(k+1)) * s% rho(k)
                        v_aver =  v_aver + (s% r(k) - s% r(k+1)) * s%  conv_vel(k)
                        cs_aver = cs_aver + (s% r(k) - s% r(k+1)) * s% csound(k)
                        kk = k
                   	end do
                   end if
                   rho_aver = rho_aver / ( s% r(sc_top) - s% r(kk) )
                   v_aver  = v_aver / ( s% r(sc_top) - s% r(kk) )
                   cs_aver = cs_aver / ( s% r(sc_top) - s% r(kk) )

                   if (cs_aver /=0) then
                    mach_max = MAXVAL(s% conv_vel(sc_top:sc_bottom))
                    mach_max = mach_max/cs_aver
                   	mach_aver = v_aver/cs_aver
                   end if

                end subroutine get_conv_aver
	  
	  
      subroutine brott_wind(id, Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z, w, ierr)
         use star_def
         integer, intent(in) :: id
         real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z ! surface values (cgs)
         ! NOTE: surface is outermost cell. not necessarily at photosphere.
         ! NOTE: don't assume that vars are set at this point.
         ! so if you want values other than those given as args,
         ! you should use values from s% xh(:,:) and s% xa(:,:) only.
         ! rather than things like s% Teff or s% lnT(:) which have not been set yet.
         real(dp), intent(out) :: w ! wind in units of Msun/year (value is >= 0)
         integer, intent(out) :: ierr

         integer :: h1, he4
         real(dp) :: Xs, Ys, Z_div_Z_solar, Teff_jump, alfa, L1, M1, R1, T1, &
            vink_wind, nieu_wind, hamann_wind, highT_w, lowT_w, Twindow
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         L1 = Lsurf
         M1 = Msurf
         R1 = Rsurf
         T1 = Tsurf

         h1 = s% net_iso(ih1)
         he4 = s% net_iso(ihe4)
         Xs = s% xa(h1,1)
         Ys = s% xa(he4,1)
         ! For wind scaling we use the ratio of iron abundance to the A09 value
         Z_div_Z_solar = s% kap_rq% Zbase/0.0142d0
         ! use Vink et al 2001, eqns 14 and 15 to set "jump" temperature
         Teff_jump = 1d3*(61.2d0 + 2.59d0*(-13.636d0 + 0.889d0*log10(Z_div_Z_solar)))

         vink_wind = 0d0
         nieu_wind = 0d0
         hamann_wind = 0d0
         w = 0

         call eval_Vink_wind(vink_wind)
         call eval_Nieuwenhuijzen_wind(nieu_wind)
         call eval_Hamann_wind(hamann_wind)

         ! use 1/10 hamann
         hamann_wind = hamann_wind/10d0

         lowT_w = max(vink_wind, nieu_wind)

         alfa = 0d0
         if (Xs > 0.7d0) then
            alfa = 1d0
         else if (Xs > 0.4d0 .and. Xs < 0.7d0) then
            alfa = (Xs - 0.4d0)/0.3d0
         end if
         highT_w = alfa * vink_wind + (1d0-alfa) * hamann_wind

         ! have a 10% Teff_jump window to switch from the lowT to the highT wind
         Twindow = Teff_jump*0.10d0
         alfa = 0d0
         if (T1 < Teff_jump - Twindow/2d0) then
            alfa = 1d0
         else if (T1 > Teff_jump - Twindow/2d0 .and. T1 < Teff_jump + Twindow/2d0) then
            alfa = ((Teff_jump + Twindow/2d0)-T1)/Twindow
         end if
         w = alfa * lowT_w + (1d0-alfa) * highT_w

         ! soften change in wind to avoid things going bad
         if (s% xtra(x_old_wind) /= 0) then
            if(abs(w) > abs(s% xtra(x_old_wind))*1.05) then
               w = s% xtra(x_old_wind)*1.05
            else if(abs(w) < abs(s% xtra(x_old_wind))*0.95) then
               w = s% xtra(x_old_wind)*0.95
            end if
         end if
         s% xtra(x_old_wind) = w

         ierr = 0

         contains

         subroutine eval_Vink_wind(w)
            real(dp), intent(inout) :: w
            real(dp) :: alfa, w1, w2, logMdot, dT, vinf_div_vesc

            ! alfa = 1 for hot side, = 0 for cool side
            if (T1 > 27500d0) then
               alfa = 1
            else if (T1 < 22500d0) then
               alfa = 0
            else
               dT = 100d0
               if (T1 > Teff_jump + dT) then
                  alfa = 1
               else if (T1 < Teff_jump - dT) then
                  alfa = 0
               else
                  alfa = (T1 - (Teff_jump - dT)) / (2*dT)
               end if
            end if
            
            if (alfa > 0) then ! eval hot side wind (eqn 24)
               vinf_div_vesc = 2.6d0 ! this is the hot side galactic value
               vinf_div_vesc = vinf_div_vesc*pow(Z_div_Z_solar,0.13d0) ! corrected for Z
               logMdot = &
                  - 6.697d0 &
                  + 2.194d0*log10(L1/Lsun/1d5) &
                  - 1.313d0*log10(M1/Msun/30) &
                  - 1.226d0*log10(vinf_div_vesc/2d0) &
                  + 0.933d0*log10(T1/4d4) &
                  - 10.92d0*pow2(log10(T1/4d4)) &
                  + 0.85d0*log10(Z_div_Z_solar)
               w1 = 10**(logMdot)
            else
               w1 = 0
            end if
            
            if (alfa < 1) then ! eval cool side wind (eqn 25)
               vinf_div_vesc = 1.3d0 ! this is the cool side galactic value
               vinf_div_vesc = vinf_div_vesc*pow(Z_div_Z_solar,0.13d0) ! corrected for Z
               logMdot = &
                  - 6.688d0 &
                  + 2.210d0*log10(L1/Lsun/1d5) &
                  - 1.339d0*log10(M1/Msun/30) &
                  - 1.601d0*log10(vinf_div_vesc/2d0) &
                  + 1.07d0*log10(T1/2d4) &
                  + 0.85d0*log10(Z_div_Z_solar)
               w2 = 10**(logMdot)
            else
               w2 = 0
            end if
            
            w = alfa*w1 + (1 - alfa)*w2
            
         end subroutine eval_Vink_wind

         subroutine eval_Nieuwenhuijzen_wind(w)
            ! Nieuwenhuijzen, H.; de Jager, C. 1990, A&A, 231, 134 (eqn 2)
            real(dp), intent(out) :: w
            real(dp) :: log10w
            include 'formats'
            log10w = -14.02d0 &
                     +1.24d0*log10(L1/Lsun) &
                     +0.16d0*log10(M1/Msun) &
                     +0.81d0*log10(R1/Rsun) !&
             !+0.85d0*log10(Z_div_Z_solar) ! we do not apply the Vink Z scaling here
            w = 10**(log10w)
         end subroutine eval_Nieuwenhuijzen_wind

         subroutine eval_Hamann_wind(w)
            ! Hamann, W.-R.; Koesterke, L.; Wessolowski, U. 1995, A&A, 299, 151
            real(dp), intent(out) :: w
            real(dp) :: log10w
            include 'formats'
            log10w = -11.95d0 &
                     +1.5d0*log10(L1/Lsun) &
                     -2.85d0*Xs &
                     + 0.85d0*log10(Z_div_Z_solar)
            w = 10**(log10w)
         end subroutine eval_Hamann_wind

      end subroutine brott_wind
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         if (.not. restart) then
            s% xtra(x_old_wind) = 0d0
            s% xtra(x_time_thermal_eq) = 0d0
            s% lxtra(lx_pre_ZAMS) = .true.
         end if
      end subroutine extras_startup
      

      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         real(dp) :: omega_crit
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0

         ! check if the star is not yet at ZAMS. If that's the case, keep rotation fixed
         ! s% xtra(x_time_thermal_eq) is set in extras_check_model
         if (s% lxtra(lx_pre_ZAMS)) then
            write(*,*) "Check thermal timescale", standard_cgrav*(s% m(1))**2/(s% r(1)*s% L(1))/secyer
            if (s% xtra(x_time_thermal_eq) > standard_cgrav*(s% m(1))**2/(s% r(1)*s% L(1))) then
               write(*,*) "thermal equilibrium for GM^2/RL, normal evolution from here"
               s% lxtra(lx_pre_ZAMS) = .false.
            else if (s% center_h1 < 0.69d0) then
               write(*,*) "Have burned significant hydrogen, normal evolution from here"
               s% lxtra(lx_pre_ZAMS) = .false.
            else if (s% star_age*secyer > 5*standard_cgrav*(s% m(1))**2/(s% r(1)*s% L(1))) then
               write(*,*) "WARNING: no equilibrium found after evolving for 5GM^2/RL"
               write(*,*) "Switching to normal evolution"
               s% lxtra(lx_pre_ZAMS) = .false.
            else
               ! keep rotation fixed
               write(*,*) "Not at ZAMS yet, keeping omega_div_omega_crit fixed"
               omega_crit = star_surface_omega_crit(id,ierr)
               if (ierr /= 0) then
                  write(*,*) "Error in star_surface_omega_crit"
                  return
               end if
               call star_set_uniform_omega(s% id, omega_crit*s% job% new_omega_div_omega_crit, ierr)
               if (ierr /= 0) then
                  write(*,*) "Error in star_set_uniform_omega"
                  return
               end if
            end if
         end if

         !!slowly turn on superad reduction to make it easier to produce pre-MS models
         !if (s% star_age >= 0.01d0) then
         !   s% superad_reduction_diff_grads_limit = 1d-2
         !else
         !   s% superad_reduction_diff_grads_limit = 10**(2-4d0*(s% star_age)/0.01d0)
         !end if

         !write(*,*) "check superad_reduction_diff_grads_limit", s% superad_reduction_diff_grads_limit

         if (s% center_he4 < 1d-3 .and. s% center_c12 < 1d-3) then
             ! use stricter timestep control on evolution of center density
             s% delta_lgRho_cntr_limit = 0.005d0 
             s% delta_lgRho_cntr_hard_limit = 0.01d0 
             s% convergence_ignore_equL_residuals = .true.
         else
             s% delta_lgRho_cntr_limit = 0.02d0 
             s% delta_lgRho_cntr_hard_limit = 0.04d0
             s% convergence_ignore_equL_residuals = .false.
         end if

      end function extras_start_step


      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         

         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model
	  


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 138
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
		 integer :: nz_test
		 integer :: tauoneidxarr(1)
		 integer :: tauoneidx
		 integer :: gamidxarr(1)
		 integer :: gamidx
		 integer :: Hpidxtotarr(1)
		 integer :: Hpidxtot
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         
		 real(dp) :: t_dyn_test, L_excess_test, E_excess_test
		 
		 
         real(dp) :: v_HI_max, v_HeI_max, v_HeII_max, v_FeCZ_max, v_HI_maxtopHp
         real(dp) :: v_HI_aver, v_HeI_aver, v_HeII_aver, v_FeCZ_aver, v_HI_avertopHp
         real(dp) :: b_HI_max, b_HeI_max, b_HeII_max, b_FeCZ_max, b_HI_maxtopHp
         real(dp) :: b_HI_aver, b_HeI_aver, b_HeII_aver, b_FeCZ_aver, b_HI_avertopHp
         real(dp) :: b_HI_surf, b_HeI_surf, b_HeII_surf, b_FeCZ_surf
         real(dp) :: b_HI_surf_max, b_HeI_surf_max, b_HeII_surf_max, b_FeCZ_surf_max
         real(dp) :: HI_hp_aver, HeI_hp_aver, HeII_hp_aver, FeCZ_hp_aver
         real(dp) :: mach_HI_top, mach_HeI_top, mach_HeII_top, mach_FeCZ_top
         real(dp) :: mach_HI_max, mach_HI_aver
         real(dp) :: mach_HeI_max, mach_HeI_aver
         real(dp) :: mach_HeII_max, mach_HeII_aver
         real(dp) :: mach_FeCZ_max, mach_FeCZ_aver
         real(dp) :: rho_HI_aver, rho_HeI_aver, rho_HeII_aver, rho_FeCZ_aver, rho_HI_avertopHp
         real(dp) :: turnover_HI, turnover_HeI, turnover_HeII, turnover_FeCZ
         real(dp) :: mach_HI_aver_ahp, mach_HeI_aver_ahp, mach_HeII_aver_ahp, mach_FeCZ_aver_ahp
         real(dp) :: v_HI_aver_ahp, v_HeI_aver_ahp, v_HeII_aver_ahp, v_FeCZ_aver_ahp
         real(dp) :: v_HI_surf, v_HeI_surf, v_HeII_surf, v_FeCZ_surf
         real(dp) :: HI_r_top, HI_r_bottom, HeI_r_top, HeI_r_bottom
         real(dp) :: HeII_r_top, HeII_r_bottom, FeCZ_r_top, FeCZ_r_bottom
         real(dp) :: HI_mass, HeI_mass, HeII_mass, FeCZ_mass
         real(dp) :: r_hp_1, r_hp_2, r_hp_3, r_hp_4, r_hp_5, r_hp_6, r_hp_7, r_hp_8
         real(dp) :: r_hp_10, r_hp_15, r_hp_20, r_hp_30, r_hp_50, r_hp_100
         real(dp) :: HI_fcmax, HeI_fcmax, HeII_fcmax, FeCZ_fcmax
         real(dp) :: HI_b_p_eq, HI_b_p_max, HeI_b_p_eq, HeI_b_p_max, HeII_b_p_eq, HeII_b_p_max
         real(dp) :: HI_B_shutoff_conv, HeI_B_shutoff_conv, HeII_B_shutoff_conv, FeCZ_B_shutoff_conv
         real(dp) :: FeCZ_b_p_eq, FeCZ_b_p_max
         real(dp) :: mixing_length_alpha, rho_surf
         real(dp) :: v_max_core, v_aver_core,b_eq_core,b_max_core,rho_aver_core, hp_aver_core, turnover_core
         real(dp) :: hp_core_top, r_core, m_core, mach_top_cz_core, mach_aver_ahp_core, rho_aver_ahp_core, v_aver_ahp_core
         real(dp) :: mach_max_core, mach_aver_core
         real(dp) :: F0
         real(dp) :: HI_p_turb_over_p, HeI_p_turb_over_p, HeII_p_turb_over_p, FeCZ_p_turb_over_p, HI_p_turb_over_ptopHp
		 real(dp) :: mdot_erupt, mdot_dutch
		 
		 integer, intent(out) :: ierr
         integer ::  i, j, k, m, num_conv_regions, sc1_top, sc2_top, sc3_top, sc4_top
         integer ::  sc1_bottom, sc2_bottom, sc3_bottom, sc4_bottom, col_count, sc_convective_core
         integer ::  hp_1, hp_2, hp_3, hp_4, hp_5, hp_6, hp_7
         integer ::  hp_8, hp_10, hp_15, hp_20, hp_30, hp_50, hp_100
         character (len=100) :: col_name
         character (len=10) :: str
         character (len=7) ::  sc1_type
         character (len=7), dimension(4) :: sc_type ! Fixed max number of cz. Can be improved
         integer, dimension(4) :: sc_top, sc_bottom
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
		 
         
         ! note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in history column options.
         ! it must not include the new column names you are adding here.
		 
		 
         ! 1) Need to initialize to zero all the columns!!!
		 
		nz_test = s% nz
		 
		t_dyn_test = 0d0
		L_excess_test = 0d0
		E_excess_test = 0d0

        b_HI_aver = 0d0
        b_HI_max= 0d0
		b_HI_maxtopHp = 0d0
		b_HI_avertopHp = 0d0
        b_HI_surf= 0d0
        b_HI_surf_max= 0d0
        HI_hp_aver= 0d0
        mach_HI_aver_ahp= 0d0
        turnover_HI= 0d0
        v_HI_surf= 0d0
        HI_r_top = 0d0
        HI_r_bottom = 0d0
        HI_mass = 0d0
        HI_fcmax = 0d0
        HI_b_p_eq = 0d0
        HI_b_p_max = 0d0
        v_HI_aver = 0d0
		v_HI_avertopHp = 0d0
        mach_HI_max = 0d0
        mach_HI_aver = 0d0

        b_HeI_aver= 0d0
        b_HeI_max= 0d0
        b_HeI_surf= 0d0
        b_HeI_surf_max= 0d0
        HeI_hp_aver= 0d0
        mach_HeI_aver_ahp= 0d0
        turnover_HeI= 0d0
        v_HeI_surf= 0d0
        HeI_r_top = 0d0
        HeI_r_bottom = 0d0
        HeI_mass = 0d0
        HeI_fcmax = 0d0
        HeI_b_p_eq = 0d0
        HeI_b_p_max = 0d0
        v_HeI_aver = 0d0
        mach_HeI_max = 0d0
        mach_HeI_aver = 0d0

        b_HeII_aver= 0d0
        b_HeII_max= 0d0
        b_HeII_surf= 0d0
        b_HeII_surf_max= 0d0
        HeII_hp_aver= 0d0
        mach_HeII_aver_ahp= 0d0
        turnover_HeII= 0d0
        v_HeII_surf= 0d0
        HeII_r_top = 0d0
        HeII_r_bottom = 0d0
        HeII_mass = 0d0
        HeII_fcmax = 0d0
        HeII_b_p_eq = 0d0
        HeII_b_p_max = 0d0
        v_HeII_aver = 0d0
        mach_HeII_max = 0d0
        mach_HeII_aver = 0d0

        b_FeCZ_aver= 0d0
        b_FeCZ_max= 0d0
        b_FeCZ_surf= 0d0
        b_FeCZ_surf_max= 0d0
        FeCZ_hp_aver= 0d0
        mach_FeCZ_aver_ahp= 0d0
        turnover_FeCZ= 0d0
        v_FeCZ_surf= 0d0
        FeCZ_r_top = 0d0
        FeCZ_r_bottom = 0d0
        FeCZ_mass = 0d0
        FeCZ_fcmax = 0d0
        FeCZ_b_p_eq = 0d0
        FeCZ_b_p_max = 0d0
        v_FeCZ_aver = 0d0
        mach_FeCZ_max = 0d0
        mach_FeCZ_aver = 0d0

        HI_B_shutoff_conv = 0d0
        HeI_B_shutoff_conv = 0d0
        HeII_B_shutoff_conv = 0d0
        FeCZ_B_shutoff_conv = 0d0

        v_max_core = 0d0
        v_aver_core = 0d0
        b_eq_core = 0d0
        b_max_core = 0d0
        rho_aver_core = 0d0
        hp_aver_core = 0d0
        hp_core_top = 0d0
        turnover_core = 0d0
        m_core = 0d0
        r_core = 0d0
        v_aver_ahp_core = 0d0
        mach_top_cz_core = 0d0
        mach_aver_ahp_core = 0d0
        rho_aver_ahp_core = 0d0
        mach_max_core = 0d0
        mach_aver_core = 0d0

        HI_p_turb_over_p = 0d0
        HeI_p_turb_over_p = 0d0
        HeII_p_turb_over_p = 0d0
        FeCZ_p_turb_over_p = 0d0
		HI_p_turb_over_ptopHp = 0d0
		
		mdot_erupt = 0d0
		mdot_dutch = 0d0
	

        F0 = 0d0 ! F0/omega_c^2 at the surface (For IGWs)
        call get_F0(s% rho(1), s% r(1), s% brunt_N2(1), F0)



         mixing_length_alpha = s% mixing_length_alpha
		 
         ! Identify top of convective core (center has singular values of e.g. density. Use s% nz -1 )
         call get_convective_core(id, sc_convective_core, ierr)
         if (sc_convective_core < s% nz) then
            !write(*,*) 'Mass Convective Core: ', s% m(sc_convective_core)/Msun
            !write(*,*) 'sc_convective_core, s nz', sc_convective_core, s% nz
            call get_conv_velocities(id, ierr, v_max_core, v_aver_core, sc_convective_core, &
                 s% nz - 1, b_eq_core,b_max_core,rho_aver_core)
            !write(*,*) 'CORE:', v_max_core/1e5, v_aver_core/1e5, b_eq_core, b_max_core,rho_aver_core
            call get_average_hp(id, ierr, sc_convective_core, s% nz - 1, hp_aver_core)
            !write(*,*) 'HP average, boundary:', hp_aver_core/Rsun, s% scale_height(sc_convective_core)/Rsun
            call get_turnover(mixing_length_alpha, v_aver_core, hp_aver_core, turnover_core)
            !write(*,*) 'Turnover V_aver+Hp_aver:', turnover_core/(3600*24)
            m_core = s% m(sc_convective_core)
            r_core = s% r(sc_convective_core)
            hp_core_top = s% scale_height(sc_convective_core)
            call get_conv_ahp(id, ierr, sc_convective_core, s% nz - 1, v_aver_ahp_core, &
                              mach_top_cz_core, mach_aver_ahp_core, rho_aver_ahp_core)
            call get_conv_aver(id, ierr, sc_convective_core, s% nz - 1,  mach_max_core, mach_aver_core)
            !write(*,*) 'v_aver_ahp, mach_top_cz, mach_aver_ahp, rho_aver_ahp,rho_aver:', &
            !                     v_aver_ahp_core, mach_top_cz_core, mach_aver_ahp_core, &
            !                     rho_aver_ahp_core, rho_aver_core
            !call get_turnover(mixing_length_alpha, v_aver_core, s% scale_height(sc_convective_core), turnover_core)
            !write(*,*) 'Turnover core V_aver+Hp_top:', turnover_core/(3600*24)
            !call get_turnover(mixing_length_alpha, v_max_core, hp_aver_core, turnover_core)
            !write(*,*) 'Turnover core Vmax+Hp_aver:', turnover_core/(3600*24)
            !call get_turnover(mixing_length_alpha, v_max_core, s% scale_height(sc_convective_core), turnover_core)
            !write(*,*) 'Turnover core Vmax+Hp_top:', turnover_core/(3600*24)
         end if


         ! Identify number of convective regions above a certain temperature  (Max 4, HI, HeI, HeII, FeCZ)

         call get_conv_regions_above_T(id,1d6,ierr,num_conv_regions)
         names(1) = 'subsurface_convective_regions'
         vals(1)  = num_conv_regions

         rho_surf = s% rho(1)
         names(2) = 'rho_surf'
         vals(2)  = rho_surf



         ! Calculate relevant column values
         do k = 1, num_conv_regions ! num_conv_regions should always be >= 1
           sc_top(k) = s% mixing_region_top(k)
           sc_bottom(k) = s% mixing_region_bottom(k)
           if (sc_top(k) .NE. 0) then
             call classify_conv_region_above_T(id, ierr, sc_top(k), sc_bottom(k), sc_type(k))
             if ( sc_type(k) == 'HI' ) then
				call get_conv_velocitiestopHp(id, ierr, v_HI_maxtopHp, v_HI_avertopHp, sc_top(k), sc_bottom(k), b_HI_avertopHp, b_HI_maxtopHp, rho_HI_avertopHp)
				call get_conv_velocities(id, ierr, v_HI_max, v_HI_aver, sc_top(k), sc_bottom(k), b_HI_aver, b_HI_max, rho_HI_aver)
                call get_average_hp(id, ierr, sc_top(k),  sc_bottom(k), HI_hp_aver)
                call get_conv_ahp(id, ierr, sc_top(k),  sc_bottom(k), v_HI_aver_ahp, mach_HI_top, mach_HI_aver_ahp, rho_HI_aver)
                call get_conv_aver(id, ierr, sc_top(k), sc_bottom(k),  mach_HI_max, mach_HI_aver)
                call get_microturb(mach_HI_aver_ahp, rho_HI_aver, rho_surf,v_HI_aver_ahp, v_HI_surf)
                call get_turnover(mixing_length_alpha, v_HI_aver, HI_hp_aver, turnover_HI)
                call get_bsurf(rho_surf, rho_HI_aver, b_HI_aver, b_HI_max, b_HI_surf, b_HI_surf_max)
                call get_conv_radii(id, ierr, sc_top(k), sc_bottom(k), HI_r_top, HI_r_bottom)
                call get_conv_mass(id, ierr, sc_top(k), sc_bottom(k), HI_mass)
                call get_max_fc(id, ierr, HI_fcmax, sc_top(k), sc_bottom(k))
                call get_pressure_eq_field(id, ierr, sc_top(k), sc_bottom (k),HI_b_p_eq,HI_b_p_max)
                !call get_B_shutoff_conv_region_above_T(id, ierr, sc_top(k), sc_bottom (k), HI_B_shutoff_conv)
                call get_max_p_turb_over_p(id, ierr, HI_p_turb_over_p, sc_top(k), sc_bottom(k))
				call get_max_p_turb_over_ptopHp(id, ierr, HI_p_turb_over_ptopHp, sc_top(k), sc_bottom(k))
				
                !write(*,*) sc_top(k), sc_bottom(k), sc_type(k), v_HI_aver, rho_HI_aver, b_HI_max, b_HI_surf, v_HI_surf
             else if ( sc_type(k) == 'HeI' ) then
                call get_conv_velocities(id, ierr, v_HeI_max, v_HeI_aver, sc_top(k), sc_bottom(k), &
                b_HeI_aver, b_HeI_max, rho_HeI_aver)
                call get_average_hp(id, ierr, sc_top(k),  sc_bottom(k), HeI_hp_aver)
                call get_conv_ahp(id, ierr, sc_top(k),  sc_bottom(k), v_HeI_aver_ahp, mach_HeI_top, &
                mach_HeI_aver_ahp, rho_HeI_aver)
                call get_conv_aver(id, ierr, sc_top(k), sc_bottom(k),  mach_HeI_max, mach_HeI_aver)
                call get_microturb(mach_HeI_aver_ahp, rho_HeI_aver, rho_surf,v_HeI_aver_ahp, v_HeI_surf)
                call get_turnover(mixing_length_alpha, v_HeI_aver, HeI_hp_aver, turnover_HeI)
                call get_bsurf(rho_surf, rho_HeI_aver, b_HeI_aver, b_HeI_max, b_HeI_surf, b_HeI_surf_max)
                call get_conv_radii(id, ierr, sc_top(k), sc_bottom(k), HeI_r_top, HeI_r_bottom)
                call get_conv_mass(id, ierr, sc_top(k), sc_bottom(k), HeI_mass)
                call get_max_fc(id, ierr, HeI_fcmax, sc_top(k), sc_bottom(k))
                call get_pressure_eq_field(id, ierr, sc_top(k), sc_bottom (k),HeI_b_p_eq,HeI_b_p_max)
                !call get_B_shutoff_conv_region_above_T(id, ierr, sc_top(k), sc_bottom (k), HeI_B_shutoff_conv)
                call get_max_p_turb_over_p(id, ierr, HeI_p_turb_over_p, sc_top(k), sc_bottom(k))
                !write(*,*) sc_top(k), sc_bottom(k), sc_type(k), v_HeI_aver, rho_HeI_aver, b_HeI_max, b_HeI_surf, v_HeI_surf
             else if ( sc_type(k) == 'HeII' ) then
                call get_conv_velocities(id, ierr, v_HeII_max, v_HeII_aver, sc_top(k), sc_bottom(k), &
                b_HeII_aver, b_HeII_max, rho_HeII_aver)
                call get_average_hp(id, ierr, sc_top(k),  sc_bottom(k), HeII_hp_aver)
                call get_conv_ahp(id, ierr, sc_top(k),  sc_bottom(k), v_HeII_aver_ahp, mach_HeII_top, &
                mach_HeII_aver_ahp, rho_HeII_aver)
                call get_conv_aver(id, ierr, sc_top(k), sc_bottom(k),  mach_HeII_max, mach_HeII_aver)
                call get_microturb(mach_HeII_aver_ahp, rho_HeII_aver, rho_surf,v_HeII_aver_ahp, v_HeII_surf)
                call get_turnover(mixing_length_alpha, v_HeII_aver, HeII_hp_aver, turnover_HeII)
                call get_bsurf(rho_surf, rho_HeII_aver, b_HeII_aver, b_HeII_max, b_HeII_surf, b_HeII_surf_max)
                call get_conv_radii(id, ierr, sc_top(k), sc_bottom(k), HeII_r_top, HeII_r_bottom)
                call get_conv_mass(id, ierr, sc_top(k), sc_bottom(k), HeII_mass)
                call get_max_fc(id, ierr, HeII_fcmax, sc_top(k), sc_bottom(k))
                call get_pressure_eq_field(id, ierr, sc_top(k), sc_bottom (k),HeII_b_p_eq,HeII_b_p_max)
                !call get_B_shutoff_conv_region_above_T(id, ierr, sc_top(k), sc_bottom (k), HeII_B_shutoff_conv)
                call get_max_p_turb_over_p(id, ierr, HeII_p_turb_over_p, sc_top(k), sc_bottom(k))
                !write(*,*) sc_top(k), sc_bottom(k), sc_type(k), v_HeII_aver, rho_HeII_aver, b_HeII_max, b_HeII_surf, v_HeII_surf
             else if ( sc_type(k) == 'FeCZ' ) then
             	  call get_conv_velocities(id, ierr, v_FeCZ_max, v_FeCZ_aver, sc_top(k), sc_bottom(k), &
                b_FeCZ_aver, b_FeCZ_max, rho_FeCZ_aver)
                call get_average_hp(id, ierr, sc_top(k),  sc_bottom(k), FeCZ_hp_aver)
                call get_conv_ahp(id, ierr, sc_top(k),  sc_bottom(k), v_FeCZ_aver_ahp, mach_FeCZ_top, &
                mach_FeCZ_aver_ahp, rho_FeCZ_aver)
                call get_conv_aver(id, ierr, sc_top(k), sc_bottom(k), mach_FeCZ_max, mach_FeCZ_aver)
                call get_microturb(mach_FeCZ_aver_ahp, rho_FeCZ_aver, rho_surf,v_FeCZ_aver_ahp, v_FeCZ_surf)
                call get_turnover(mixing_length_alpha, v_FeCZ_aver, FeCZ_hp_aver, turnover_FeCZ)
                call get_bsurf(rho_surf, rho_FeCZ_aver, b_FeCZ_aver, b_FeCZ_max, b_FeCZ_surf, b_FeCZ_surf_max)
                call get_conv_radii(id, ierr, sc_top(k), sc_bottom(k), FeCZ_r_top, FeCZ_r_bottom)
                call get_conv_mass(id, ierr, sc_top(k), sc_bottom(k), FeCZ_mass)
                call get_max_fc(id, ierr, FeCZ_fcmax, sc_top(k), sc_bottom(k))
                call get_pressure_eq_field(id, ierr, sc_top(k), sc_bottom (k),FeCZ_b_p_eq,FeCZ_b_p_max)
                !call get_B_shutoff_conv_region_above_T(id, ierr, sc_top(k), sc_bottom (k), FeCZ_B_shutoff_conv)
                call get_max_p_turb_over_p(id, ierr, FeCZ_p_turb_over_p, sc_top(k), sc_bottom(k))
                ! write(*,*) 'F0/1d15, brunt/muhz , omega/muhz', F0/1d15,(((s% brunt_N2(1))**0.5)*1d6)/6.28, (1/turnover_FeCZ)*1d6
                !write(*,*) sc_top(k), sc_bottom(k), sc_type(k), v_FeCZ_aver, rho_FeCZ_aver, b_FeCZ_max, b_FeCZ_surf, v_FeCZ_surf
             end if
          end if
         end do
         ! Store relevant column values (8x4) = 32 columns

         names(3) = 'v_HI_surf'
         vals(3)  = v_HI_surf
         names(4) = 'b_HI_surf'
         vals(4)  = b_HI_surf
         names(5) = 'b_HI_surf_max'
         vals(5)  = b_HI_surf_max
         names(6) = 'b_HI_aver'
         vals(6)  = b_HI_aver
         names(7) = 'b_HI_max'
         vals(7)  = b_HI_max
         names(8) = 'HI_hp_aver'
         vals(8)  = HI_hp_aver
         names(9) = 'mach_HI_aver_ahp'
         vals(9)  = mach_HI_aver_ahp
         names(10) = 'turnover_HI'
         vals(10)  = turnover_HI


         names(11) = 'v_HeI_surf'
         vals(11)  = v_HeI_surf
         names(12) = 'b_HeI_surf'
         vals(12)  = b_HeI_surf
         names(13) = 'b_HeI_surf_max'
         vals(13)  = b_HeI_surf_max
         names(14) = 'b_HeI_aver'
         vals(14)  = b_HeI_aver
         names(15) = 'b_HeI_max'
         vals(15)  = b_HeI_max
         names(16) = 'HeI_hp_aver'
         vals(16)  = HeI_hp_aver
         names(17) = 'mach_HeI_aver_ahp'
         vals(17)  = mach_HeI_aver_ahp
         names(18) = 'turnover_HeI'
         vals(18)  = turnover_HeI

         names(19) = 'v_HeII_surf'
         vals(19)  = v_HeII_surf
         names(20) = 'b_HeII_surf'
         vals(20)  = b_HeII_surf
         names(21) = 'b_HeII_surf_max'
         vals(21)  = b_HeII_surf_max
         names(22) = 'b_HeII_aver'
         vals(22)  = b_HeII_aver
         names(23) = 'b_HeII_max'
         vals(23)  = b_HeII_max
         names(24) = 'HeII_hp_aver'
         vals(24)  = HeII_hp_aver
         names(25) = 'mach_HeII_aver_ahp'
         vals(25)  = mach_HeII_aver_ahp
         names(26) = 'turnover_HeII'
         vals(26)  = turnover_HeII

         names(27) = 'v_FeCZ_surf'
         vals(27)  = v_FeCZ_surf
         names(28) = 'b_FeCZ_surf'
         vals(28)  = b_FeCZ_surf
         names(29) = 'b_FeCZ_surf_max'
         vals(29)  = b_FeCZ_surf_max
         names(30) = 'b_FeCZ_aver'
         vals(30)  = b_FeCZ_aver
         names(31) = 'b_FeCZ_max'
         vals(31)  = b_FeCZ_max
         names(32) = 'FeCZ_hp_aver'
         vals(32)  = FeCZ_hp_aver
         names(33) = 'mach_FeCZ_aver_ahp'
         vals(33)  = mach_FeCZ_aver_ahp
         names(34) = 'turnover_FeCZ'
         vals(34)  = turnover_FeCZ

         names(35) = 'HI_r_top'
         vals(35)  = HI_r_top
         names(36) = 'HI_r_bottom'
         vals(36)  = HI_r_bottom

         names(37) = 'HeI_r_top'
         vals(37)  = HeI_r_top
         names(38) = 'HeI_r_bottom'
         vals(38)  = HeI_r_bottom

         names(39) = 'HeII_r_top'
         vals(39)  = HeII_r_top
         names(40) = 'HeII_r_bottom'
         vals(40)  = HeII_r_bottom

         names(41) = 'FeCZ_r_top'
         vals(41)  = FeCZ_r_top
         names(42) = 'FeCZ_r_bottom'
         vals(42)  = FeCZ_r_bottom

         names(43) = 'HI_mass'
         vals(43)  = HI_mass
         names(44) = 'HeI_mass'
         vals(44)  = HeI_mass
         names(45) = 'HeII_mass'
         vals(45)  = HeII_mass
         names(46) = 'FeCZ_mass'
         vals(46)  = FeCZ_mass

         names(47) = 'HI_Fc_max'
         vals(47)  = HI_fcmax
         names(48) = 'HeI_Fc_max'
         vals(48)  = HeI_fcmax
         names(49) = 'HeII_Fc_max'
         vals(49)  = HeII_fcmax
         names(50) = 'FeCZ_Fc_max'
         vals(50)  = FeCZ_fcmax



!        Pressure scale Heigths (0,1,2,3,4,5,6,7,8)
         call get_hp_radii(id, ierr, 1d0, hp_1)
         call get_hp_radii(id, ierr, 2d0, hp_2)
         call get_hp_radii(id, ierr, 3d0, hp_3)
         call get_hp_radii(id, ierr, 4d0, hp_4)
         call get_hp_radii(id, ierr, 5d0, hp_5)
         call get_hp_radii(id, ierr, 6d0, hp_6)
         call get_hp_radii(id, ierr, 7d0, hp_7)
         call get_hp_radii(id, ierr, 8d0, hp_8)
         call get_hp_radii(id, ierr, 10d0, hp_10)
         call get_hp_radii(id, ierr, 15d0, hp_15)
         call get_hp_radii(id, ierr, 20d0, hp_20)
         call get_hp_radii(id, ierr, 30d0, hp_30)
         call get_hp_radii(id, ierr, 50d0, hp_50)
         call get_hp_radii(id, ierr, 100d0, hp_100)

         r_hp_1 = s% r(hp_1)
         r_hp_2 = s% r(hp_2)
         r_hp_3 = s% r(hp_3)
         r_hp_4 = s% r(hp_4)
         r_hp_5 = s% r(hp_5)
         r_hp_6 = s% r(hp_6)
         r_hp_7 = s% r(hp_7)
         r_hp_8 = s% r(hp_8)
         r_hp_10 = s% r(hp_10)
         r_hp_15 = s% r(hp_15)
         r_hp_20 = s% r(hp_20)
         r_hp_30 = s% r(hp_30)
         r_hp_50 = s% r(hp_50)
         r_hp_100 = s% r(hp_100)

         names(51) = 'r_hp_1'
         vals(51)  = r_hp_1
         names(52) = 'r_hp_2'
         vals(52)  = r_hp_2
         names(53) = 'r_hp_3'
         vals(53)  = r_hp_3
         names(54) = 'r_hp_4'
         vals(54)  = r_hp_4
         names(55) = 'r_hp_5'
         vals(55)  = r_hp_5
         names(56) = 'r_hp_6'
         vals(56)  = r_hp_6
         names(57) = 'r_hp_7'
         vals(57)  = r_hp_7
         names(58) = 'r_hp_8'
         vals(58)  = r_hp_8
         names(59) = 'r_hp_10'
         vals(59)  = r_hp_10
         names(60) = 'r_hp_15'
         vals(60)  = r_hp_15
         names(61) = 'r_hp_20'
         vals(61)  = r_hp_20
         names(62) = 'r_hp_30'
         vals(62)  = r_hp_30
         names(63) = 'r_hp_50'
         vals(63)  = r_hp_50
         names(64) = 'r_hp_100'
         vals(64)  = r_hp_100

         names(65) = 'HI_b_p_eq'
         vals(65) = HI_b_p_eq
         names(66) = 'HI_b_p_max'
         vals(66) = HI_b_p_max

         names(67) = 'HeI_b_p_eq'
         vals(67) = HeI_b_p_eq
         names(68) = 'HeI_b_p_max'
         vals(68) = HeI_b_p_max

         names(69) = 'HeII_b_p_eq'
         vals(69) = HeII_b_p_eq
         names(70) = 'HeII_b_p_max'
         vals(70) = HeII_b_p_max

         names(71) = 'FeCZ_b_p_eq'
         vals(71) = FeCZ_b_p_eq
         names(72) = 'FeCZ_b_p_max'
         vals(72) = FeCZ_b_p_max


         names(73) = 'v_max_core'
         names(74) = 'v_aver_core'
         names(75) = 'b_eq_core'
         names(76) = 'b_max_core'
         names(77) = 'rho_aver_core'
         names(78) = 'hp_aver_core'
         names(79) = 'hp_core_top'
         names(80) = 'turnover_core'
         names(81) = 'm_core'
         names(82) = 'r_core'

         vals(73) = v_max_core
         vals(74) = v_aver_core
         vals(75) = b_eq_core
         vals(76) = b_max_core
         vals(77) = rho_aver_core
         vals(78) = hp_aver_core
         vals(79) = hp_core_top
         vals(80) = turnover_core
         vals(81) = m_core
         vals(82) = r_core



         names(83) = 'v_aver_ahp_core'
         names(84) = 'mach_top_cz_core'
         names(85) = 'mach_aver_ahp_core'
         names(86) = 'rho_aver_ahp_core'

         vals(83) = v_aver_ahp_core
         vals(84) = mach_top_cz_core
         vals(85) = mach_aver_ahp_core
         vals(86) = rho_aver_ahp_core

         names(87) = 'HI_B_shutoff_conv'
         vals(87) = HI_B_shutoff_conv
         names(88) = 'HeI_B_shutoff_conv'
         vals(88) = HeI_B_shutoff_conv
         names(89) = 'HeII_B_shutoff_conv'
         vals(89) = HeII_B_shutoff_conv
         names(90) = 'FeCZ_B_shutoff_conv'
         vals(90) = FeCZ_B_shutoff_conv

         names(91) = 'F0_div_omega_c' ! 1/2 rho(R) N(R) R^3. Multiply by omega_c^2 to obtain F0
         vals(91) = F0

         names(92) = 'v_HI_aver'
         names(93) = 'mach_HI_max'
         names(94) = 'mach_HI_aver'

         vals(92) = v_HI_aver
         vals(93) = mach_HI_max
         vals(94) = mach_HI_aver

         names(94) = 'v_HeI_aver'
         names(95) = 'mach_HeI_max'
         names(96) = 'mach_HeI_aver'

         vals(94) = v_HeI_aver
         vals(95) = mach_HeI_max
         vals(96) = mach_HeI_aver

         names(97) = 'v_HeII_aver'
         names(98) = 'mach_HeII_max'
         names(99) = 'mach_HeII_aver'

         vals(97) = v_HeII_aver
         vals(98) = mach_HeII_max
         vals(99) = mach_HeII_aver

         names(100) = 'v_FeCZ_aver'
         names(101) = 'mach_FeCZ_max'
         names(102) = 'mach_FeCZ_aver'

         vals(100) = v_FeCZ_aver
         vals(101) = mach_FeCZ_max
         vals(102) = mach_FeCZ_aver

         names(103) = 'mach_max_core'
         vals(103) = mach_max_core
         names(104) = 'mach_aver_core'
         vals(104) = mach_aver_core

         names(105) = 'HI_p_turb_over_p'
         vals(105) = HI_p_turb_over_p
         names(106) = 'HeI_p_turb_over_p'
         vals(106) = HeI_p_turb_over_p

         names(107) = 'HeII_p_turb_over_p'
         vals(107) = HeII_p_turb_over_p
         names(108) = 'FeCZ_p_turb_over_p'
         vals(108) = FeCZ_p_turb_over_p
         
         names(109) = 'pre_ZAMS'
         names(110) = 'max_superad_reduction'
		 
		 names(111) = "Gamma_tot_max"
		 names(112) = "Gamma_tot_max_idx"
		 names(113) = "Gamma_rad_max"
		 names(114) = "Gamma_rad_max_idx"
		 names(115) = "Gamma_tot_surf"
		 names(116) = "Gamma_rad_surf"
		 names(117) = 'testtau'
		 names(118) = "Gamma_tot_tau1"
		 names(119) = "Gamma_rad_tau1"
		 names(120) = "testLedd"
		 names(121) = "testT1"
		 names(122) = "testT2"
		 names(123) = "testT1_1"
		 names(124) = "testT1_2"
		 
         if (s% lxtra(lx_pre_ZAMS)) then
            vals(109) = 1d0
         else
            vals(109) = 0d0
         end if
         vals(110) = maxval(s% superad_reduction_factor(1:s% nz))
		 
		 write(*,*) "Check max superad reduction", vals(110)
		
		tauoneidxarr = minloc(abs(s% tau(1:s% nz)-1.0d0))
		
		tauoneidx = tauoneidxarr(1)
		
		write(*,*) "tauoneidx", tauoneidx
		
		
		gamidxarr = minloc(abs(s% T(1:s% nz)-1.0d6))
		
		gamidx = gamidxarr(1)
		
		Hpidxtotarr = minloc(abs(s% Pgas(1:s% nz)-(s% Pgas(1)*2.71828d0)))
		Hpidxtot = Hpidxtotarr(1)
		
		
		call get_mloss_forwind(id, ierr, s% nz, mdot_erupt)
		
		t_dyn_test = ((s% r(1))**3/(6.674e-8*s% m(1)))**(1.0/2.0)
		L_excess_test = (s% L(1)) - (pi4*clight*s% cgrav(1)*s% m_grav(1)/(s% opacity(1)))
		E_excess_test = L_excess_test*t_dyn_test
		
		vals(111) = MAXVAL(s% L(1:gamidx) / (pi4*clight*s% cgrav(1:gamidx)*s% m_grav(1:gamidx)/(s% opacity(1:gamidx))))
		vals(112) = MAXLOC(s% L(1:gamidx) / (pi4*clight*s% cgrav(1:gamidx)*s% m_grav(1:gamidx)/(s% opacity(1:gamidx))), DIM=1)
		vals(113) = MAXVAL(4d0*crad/3d0*pow4(s% T(1:gamidx))/s% Peos(1:gamidx)*s%gradT(1:gamidx))
		vals(114) = MAXLOC(4d0*crad/3d0*pow4(s% T(1:gamidx))/s% Peos(1:gamidx)*s%gradT(1:gamidx), DIM=1)
		
		vals(115) = s% L(1) / (pi4*clight*s% cgrav(1)*s% m_grav(1)/(s% opacity(1)))
		vals(116) = 4d0*crad/3d0*pow4(s% T(1))/s% Peos(1)*s%gradT(1)
		vals(117) = s%tau(tauoneidx)
		vals(118) = s% L(tauoneidx) / (pi4*clight*s% cgrav(tauoneidx)*s% m_grav(tauoneidx)/(s% opacity(tauoneidx)))
		vals(119) = 4d0*crad/3d0*pow4(s% T(tauoneidx))/s% Peos(tauoneidx)*s%gradT(tauoneidx)
		vals(120) = pi4*clight*s% cgrav(1)*s% m_grav(1)/(s% opacity(1))
		vals(121) = s% T(1)
		vals(122) = s% T(gamidx)
		vals(123) = s% T(1)
		vals(124) = s% T(1)

        names(125) = 'HI_p_turb_over_ptopHp'
        vals(125) = HI_p_turb_over_ptopHp
		
		names(126) = 'v_HI_maxtopHp'
		vals(126) = v_HI_maxtopHp
		
		names(127) = 'v_HI_max'
		vals(127) = v_HI_max
		
		names(128) = 'v_max_topHp'
		vals(128) = MAXVAL(s% conv_vel(1:Hpidxtot))
		
		names(129) = 't_dyn_test'
		vals(129) = t_dyn_test
		
		names(130) = 'L_excess_test'
		vals(130) = L_excess_test
		
		names(131) = 'E_excess_test'
		vals(131) = E_excess_test
		
		names(132) = 'mdot_erupt'
		vals(132) = s% xtra(5)
		
		names(133) = 'mdot_dutch'
		vals(133) = s% xtra(1)
		
		names(134) = 'mdot_test_total'
		vals(134) = s% xtra(3)
		
		names(135) = 'Mabove_crit'
		vals(135) = s% xtra(7)
		
		names(136) = 'mlost_final'
		vals(136) = s% xtra(8)
		
		names(137) = 'tdyn_final'
		vals(137) = s% xtra(9)
		
		names(138) = 'mlost_loc'
		vals(138) = s% xtra(10)

		write(*,*) "mdot_erupt", s% xtra(5)
		write(*,*) "mdot_dutch", s% xtra(1)
		write(*,*) "mdot_total", s% xtra(3)
		write(*,*) "mdot_total_real", ABS(s% star_mdot)
		
		write(*,*) 'Mabove_crit', s% xtra(7)
		write(*,*) 'mlost_final', s% xtra(8)
		write(*,*) 'tdyn_final', s% xtra(9)
		write(*,*) 'mlost_loc', s% xtra(10)

	
      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 7
		 !how_many_extra_profile_columns = 2
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, j
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.
         names(1) = "Lrad_div_Leddd"
         names(2) = "Lrad_div_Leddd2"
		 names(3) = "Ledd_opacity"
		 names(4) = "Ledd_electron"
		 names(5) = "Gamma_opacity"
		 names(6) = "Gamma_electron"
		 !names(7) = 'opacity_test'
         do k=1,s% nz
         vals(k,1) = 4d0*crad/3d0*pow4(s% T(k))/s% Peos(k)*s%gradT(k)
         vals(k,2) = 4d0*crad/3d0*pow4(s% T(k))/s% Peos(k)*s%mlt_gradT(k)
		 vals(k,3) = pi4*clight*s% cgrav(k)*s% m_grav(k)/(s% opacity(k)*Lsun)
		 vals(k,4) = pi4*clight*s% cgrav(k)*s% m_grav(k)/(0.2d0*(1+s% x(k))*Lsun)
		 vals(k,5) = s% L(k) / vals(k,3)
		 vals(k,6) = s% L(k) / vals(k,4)
		 ! vals(k,7) = pi4*clight*s% cgrav(k)*s% m_grav(k)/(*Lsun)
		 
         end do
         
      end subroutine data_for_extra_profile_columns


      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
		 integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
		 integer :: k
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra history header item
         ! also set how_many_extra_history_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra profile header item
         ! also set how_many_extra_profile_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going

         ! to save a profile, 
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! check for Lnuc ~ L, this is used in star step to determine if we
         ! keep rotation fixed
         if (abs(log10(abs(s% L_nuc_burn_total * Lsun / s% L(1)))) < 0.005) then
            s% xtra(x_time_thermal_eq) = s% xtra(x_time_thermal_eq) + s% dt_old
         else
            s% xtra(x_time_thermal_eq) = 0d0
         end if

         !if (s% center_he4 < 1d-3 .and. s% center_c12 < 1d-3) then
         !    extras_finish_step = terminate
         !    write(*,*) "Terminate due to carbon depletion"
         !end if
         if (s% center_he4 < 1d-2 .and. s% center_c12 < 1d-2 .and. s% center_o16 < 1d-2) then
             extras_finish_step = terminate
             write(*,*) "Terminate due to oxygen depletion"
         end if

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve

      end module run_star_extras
      
