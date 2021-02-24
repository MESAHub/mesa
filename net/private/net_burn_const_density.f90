! ***********************************************************************
!
!   Copyright (C) 2016-2019  Bill Paxton & The MESA Team
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

      module net_burn_const_density
      use const_def
      use math_lib
      use chem_def
      use net_def
         
      use utils_lib, only: is_bad
      
      use net_burn_support, only: netint
         
      implicit none
      
      
      !logical, parameter :: use_ludcmp = .true.
      logical, parameter :: use_ludcmp = .false.
      
      !logical, parameter :: show_mesa_rates = .true.
      logical, parameter :: show_mesa_rates = .false.


      

      contains
      
      
      integer function burn_const_density_1_zone_work_size(g) result(sz)
         use net_initialize, only: work_size
         use net_approx21, only: approx21_nrat
         type (Net_General_Info), pointer :: g
         integer :: net_lwork, num_reactions, species, nvar
         include 'formats'
         num_reactions = g% num_reactions
         if (g% doing_approx21) num_reactions = approx21_nrat
         species = g% num_isos
         nvar = species + 1
         net_lwork = work_size(g)
         sz = net_lwork + 2*num_reactions + 2*nvar*nvar         
         ! dratdumdy1, dratdumdy2, dens_dfdy1, dmat1
      end function burn_const_density_1_zone_work_size
         

      subroutine get_pointers( &
            g, burn_lwork, burn_work_array, species, nvar, num_reactions, &
            dratdumdy1, dratdumdy2, dens_dfdy, dmat, i, ierr)
         use net_approx21, only: approx21_nrat
         type (Net_General_Info), pointer :: g
         real(dp), pointer :: burn_work_array(:)
         integer, intent(in) :: burn_lwork, species, nvar, num_reactions
         real(dp), pointer, dimension(:) :: dratdumdy1, dratdumdy2
         real(dp), pointer, dimension(:,:) :: dens_dfdy, dmat
         integer, intent(inout) :: i
         integer, intent(out) :: ierr
         integer :: sz
         include 'formats'
         ierr = 0
         sz = num_reactions
         if (g% doing_approx21) sz = approx21_nrat
         dratdumdy1(1:sz) => burn_work_array(i+1:i+sz); i = i+sz
         dratdumdy2(1:sz) => burn_work_array(i+1:i+sz); i = i+sz
         
         sz = nvar*nvar
         dens_dfdy(1:nvar,1:nvar) => burn_work_array(i+1:i+sz); i = i+sz
         dmat(1:nvar,1:nvar) => burn_work_array(i+1:i+sz); i = i+sz
         
      end subroutine get_pointers

      
      subroutine burn_const_density_1_zone( &
            net_handle, eos_handle, species, nvar, num_reactions, t_start, t_end, &
            starting_x, starting_log10T, log10Rho, &
            get_eos_info_for_burn_at_const_density, &
            rate_factors, weak_rate_factor, reaction_Qs, reaction_neuQs, &
            screening_mode, &
            stptry_in, max_steps, eps, odescal, &
            use_pivoting, trace, burn_dbg, burner_finish_substep, &
            burn_lwork, burn_work_array, net_lwork, net_work_array, &
            ending_x, eps_nuc_categories, ending_log10T, avg_eps_nuc, eps_neu_total, &
            nfcn, njac, nstep, naccpt, nrejct, ierr)
         use num_def
         use num_lib 
         use mtx_lib
         use mtx_def
         use rates_def, only: rates_reaction_id_max, reaction_Name
         use rates_lib, only: rates_reaction_id
         use net_initialize, only: set_rate_ptrs, setup_net_info, work_size
         use chem_lib, only: basic_composition_info, get_Q
         use net_initialize, only: work_size
         use net_approx21, only: approx21_nrat
         
         integer, intent(in) :: net_handle, eos_handle, species, nvar, num_reactions
         real(dp), intent(in) :: t_start, t_end, starting_x(:) ! (species)
         real(dp), intent(in) :: starting_log10T, log10Rho
         interface
            include 'burner_const_density_get_eos_info.inc'
         end interface
         real(dp), intent(in), pointer :: rate_factors(:) ! (num_reactions)
         real(dp), intent(in) :: weak_rate_factor
         real(dp), pointer, intent(in) :: reaction_Qs(:) ! (rates_reaction_id_max)
         real(dp), pointer, intent(in) :: reaction_neuQs(:) ! (rates_reaction_id_max)
         integer, intent(in) :: screening_mode ! see screen_def
         real(dp), intent(in) :: stptry_in ! try this for 1st step.  0 means try in 1 step.
         integer, intent(in) :: max_steps ! maximal number of allowed steps.
         real(dp), intent(in) :: eps, odescal ! tolerances.  e.g., set both to 1d-6
         logical, intent(in) :: use_pivoting ! for matrix solves
         logical, intent(in) :: trace, burn_dbg
         interface
            include 'burner_finish_substep.inc'
         end interface
         integer, intent(in) :: net_lwork, burn_lwork
         real(dp), intent(inout), pointer :: burn_work_array(:) ! (burn_lwork)
         real(dp), intent(inout), pointer :: net_work_array(:) ! (net_lwork)
         real(dp), intent(inout) :: ending_x(:) ! (species)
         real(dp), intent(inout) :: eps_nuc_categories(:) ! (num_categories)
         real(dp), intent(out) :: ending_log10T, avg_eps_nuc, eps_neu_total
         integer, intent(out) :: nfcn    ! number of function evaluations
         integer, intent(out) :: njac    ! number of jacobian evaluations
         integer, intent(out) :: nstep   ! number of computed steps
         integer, intent(out) :: naccpt  ! number of accepted steps
         integer, intent(out) :: nrejct  ! number of rejected steps
         integer, intent(out) :: ierr
         
         type (Net_General_Info), pointer :: g
         integer :: ijac, lrd, lid, lout, i, j, ir, idid, sz
         logical :: okay
         real(dp) :: temp, rho, lgT, lgRho, r, prev_lgRho, prev_lgT
         
         integer :: stpmax, imax_dydx, nstp
         real(dp) :: &
            h, start, stptry, stpmin, stopp, max_dydx, abs_max_dydx, lnT, &
            eta, d_eta_dlnT, d_eta_dlnRho, Cv, d_Cv_dlnT, burn_ergs, dx
                  
         real(dp), dimension(nvar), target :: starting_y_a, ending_y_a, save_x_a
         real(dp), dimension(:), pointer :: starting_y, ending_y, save_x
         real(dp), dimension(:), pointer :: dratdumdy1, dratdumdy2

         real(dp), dimension(:,:), pointer :: dens_dfdy, dmat

         real(dp) :: xh, xhe, z, abar, zbar, z2bar, z53bar, ye, mass_correction, sumx
         real(dp) :: aion(species)
      
         logical :: dbg
         
         type (Net_Info), target :: net_info_target
         type (Net_Info), pointer :: n
         
         real(dp), dimension(:), pointer :: &
            rate_raw, rate_raw_dT, rate_raw_dRho, &
            rate_screened, rate_screened_dT, rate_screened_dRho
         integer :: iwork, cid
         
         include 'formats'
         
         !dbg = .true.
         dbg = burn_dbg
         
         if (dbg) then
            do i=1,species
               write(*,2) 'starting_x', i, starting_x(i)
            end do
            write(*,1) 'starting_log10T', starting_log10T
            write(*,1) 'log10Rho', log10Rho
            write(*,1) 't_start', t_start
            write(*,1) 't_end', t_end
            write(*,1) 'stptry_in', stptry_in
            write(*,1) 'eps', eps
            write(*,1) 'odescal', odescal
            write(*,2) 'max_steps', max_steps
            !stop 'burn_const_density_1_zone'
         end if

         starting_y => starting_y_a
         ending_y => ending_y_a
         save_x => save_x_a
         n => net_info_target

         d_eta_dlnRho = 0d0
         
         lgT = starting_log10T
         lnT = lgT*ln10
         temp = exp10(lgT)
         lgRho = log10Rho
         rho = exp10(lgRho)
         prev_lgT = -1d99
         prev_lgRho = -1d99
         
         ierr = 0
         call get_net_ptr(net_handle, g, ierr)
         if (ierr /= 0) then
            write(*,*) 'invalid handle for burn_const_density_1_zone'
            return
         end if

         if (g% num_isos /= species) then
            write(*,*) 'invalid species', species
            return
         end if
         
         if (g% num_reactions /= num_reactions) then
            write(*,*) 'invalid num_reactions', num_reactions
            return
         end if
         
         nfcn = 0
         njac = 0
         nstep = 0
         naccpt = 0
         nrejct = 0
         
         do i=1,species
            cid = g% chem_id(i)
            if (cid < 0) cid = g% approx21_ye_iso
            aion(i) = chem_isos% Z_plus_N(cid)
            save_x(i) = starting_x(i)
            ending_x(i) = starting_x(i)
            starting_y(i) = starting_x(i)/aion(i)
            ending_y(i) = starting_y(i)
         end do
         starting_y(nvar) = lnT
         ending_y(nvar) = lnT
         
         start = t_start
         stptry = stptry_in
         if (stptry == 0d0) stptry = t_end
         
         !write(*,1) 'stptry', stptry
         
         stpmin = min(t_end*1d-20,stptry*1d-6)
         stopp = t_end
         stpmax = max_steps
         
         if (dbg) write(*,2) 'call set_rate_ptrs', burn_lwork
         call set_rate_ptrs(g, &
            rate_screened, rate_screened_dT, rate_screened_dRho, &
            rate_raw, rate_raw_dT, rate_raw_dRho, &
            burn_lwork, burn_work_array, iwork, ierr)
         if (ierr /= 0) then
            if (dbg) write(*,*) 'failed in set_ptrs_in_work'
            return
         end if
         
         if (dbg) write(*,2) 'call setup_net_info', iwork
         call setup_net_info( &
            g, n, eps_nuc_categories,  &
            screening_mode, &
            rate_screened, rate_screened_dT, rate_screened_dRho, &
            rate_raw, rate_raw_dT, rate_raw_dRho, burn_lwork, burn_work_array, &
            .false., .false., iwork, ierr)
         if (dbg) write(*,*) 'done setup_net_info'
         if (ierr /= 0) then
            if (dbg) write(*,*) 'failed in setup_net_info'
            return
         end if
      
         if (dbg) write(*,*) 'call get_pointers'
         call get_pointers( &
            g, burn_lwork, burn_work_array, species, nvar, num_reactions, &
            dratdumdy1, dratdumdy2, dens_dfdy, dmat, iwork, ierr)
         if (ierr /= 0) return

         call basic_composition_info( &
            species, g% chem_id, starting_x, xh, xhe, z, &
            abar, zbar, z2bar, z53bar, ye, mass_correction, sumx)

         !write(*,*) 'stptry for netint', stptry
         !if (stptry == 0) then
            !stptry = 1d-16 ! TESTING
            !stptry = max((stopp - start) * 1.0d-10,1.0d-12)
            !stpmin = stptry * 1.0d-12
         !else
         !   stpmin = stptry * (stopp - start)
         !end if
         !write(*,*) 'stpmin dt for netint', stpmin, stopp - start
         !stop
         stptry = max(start * 1.0d-10,1.0d-16)
         stpmin = stptry * 1.0d-12
         
         if (dbg) write(*,*) 'call netint'
         call netint( &
            start,stptry,stpmin,max_steps,stopp,ending_y, &
            eps,species,nvar,naccpt,nrejct,nstep,odescal,dens_dfdy,dmat, &
            burner_derivs,burner_jakob,burner_finish_substep,ierr)
         if (dbg) write(*,*) 'done netint'
         if (ierr /= 0) then
            return
            write(*,*) 'netint ierr'
            stop
         end if

         ! set burn_ergs according to change in abundances
         burn_ergs = 0
         do i=1,species
            ending_x(i) = ending_y(i)*aion(i)
            dx = ending_x(i) - save_x(i)
            !write(*,2) 'dx aion end_x', i, dx, aion(i), ending_x(i)
            cid = g% chem_id(i)             
            burn_ergs = burn_ergs + &
               (get_Q(chem_isos,cid))*dx/chem_isos% Z_plus_N(cid)
         end do
         burn_ergs = burn_ergs*Qconv
         avg_eps_nuc = burn_ergs/(t_end - t_start) - eps_neu_total
         ending_log10T = ending_y(nvar)/ln10
         
         !write(*,*) 'starting, ending logT', starting_y(nvar)/ln10, ending_log10T
      
      contains
         
         subroutine burner_derivs(x,h,y,f,nvar,ierr)
            integer, intent(in) :: nvar
            real(dp) :: x, h, y(:), f(:)
            integer, intent(out) :: ierr
            integer, parameter :: ld_dfdx = 0
            real(dp), target :: dfdx_arry(ld_dfdx,nvar)
            real(dp), pointer :: dfdx(:,:)
            real(dp) :: dxdt_sum, dxdt_sum_aprox21, &
               Z_plus_N, xsum, r, r1, r2
            integer :: i, ir, ci, j, k, ibad
            logical :: okay, reuse_rates

            real(dp), target :: f21_a(nvar)
            real(dp), pointer :: f21(:)
            
            include 'formats'
            
            ierr = 0
            nfcn = nfcn + 1
            dfdx => dfdx_arry
            reuse_rates = .false. ! use_aprox21_rates .or. set_aprox21_rates
            call jakob_or_derivs(x,h,y,f,dfdx,reuse_rates,ierr)
            if (ierr /= 0) return            
         
         end subroutine burner_derivs

         subroutine burner_jakob(x,h,y,dfdy,nvar,ierr)
            integer, intent(in) :: nvar
            real(dp) :: x, h, y(:)
            real(dp), pointer :: dfdy(:,:)
            integer, intent(out) :: ierr
            real(dp), target :: f_arry(0)
            real(dp), pointer :: f(:)
            
            real(dp), target :: dfdy21_a(nvar,nvar)
            real(dp), pointer :: dfdy21(:,:)
            real(dp) :: Z_plus_N, df_t, df_m
            integer :: i, ci, j, cj
            logical :: okay, reuse_rates
            include 'formats'

            ierr = 0
            njac = njac + 1
            f => f_arry

            reuse_rates = .false.
            call jakob_or_derivs(x,h,y,f,dfdy,reuse_rates,ierr)
            if (ierr /= 0) return
                     
         end subroutine burner_jakob

         subroutine jakob_or_derivs(time,dt,y,f,dfdy,reuse_rates,ierr)
            use chem_lib, only: basic_composition_info
            use net_eval, only: eval_net
            use rates_def, only: rates_reaction_id_max, i_rate, i_rate_dT, i_rate_dRho
            use interp_1d_lib, only: interp_value
         
            real(dp) :: time, dt, y(:), f(:)
            real(dp), pointer :: dfdy(:,:)
            logical, intent(in) :: reuse_rates
            integer, intent(out) :: ierr
         
            real(dp) :: T, lgT, rate_limit, rat, dratdt, dratdd
            real(dp) :: eps_nuc
            real(dp) :: d_eps_nuc_dT
            real(dp) :: d_eps_nuc_dRho
            real(dp) :: d_eps_nuc_dx(nvar) 
            real(dp) :: dxdt(nvar)
            real(dp) :: d_dxdt_dRho(nvar)
            real(dp) :: d_dxdt_dT(nvar)
            real(dp) :: d_dxdt_dx(nvar, nvar)
            
            logical ::  rates_only, dxdt_only, okay
            integer :: i, j, k, ir

            real(dp), pointer, dimension(:) :: actual_Qs, actual_neuQs
            real(dp) :: xsum
            logical, pointer :: from_weaklib(:)
            
            real(dp), target :: x_a(nvar), dfdx_a(nvar,nvar)
            real(dp), pointer :: x(:), dfdx(:,:)
            real(dp) :: d_eta_dlnRho
         
            include 'formats'
         
            ierr = 0

            x => x_a
            dfdx => dfdx_a
            
            actual_Qs => null()
            actual_neuQs => null()
            from_weaklib => null()
            
            xsum = 0 
            do i=1,species
               if (is_bad(y(i))) then
                  ierr = -1
                  write(*,2) 'net_burn_const_density failed in jakob_or_derivs: bad y(i)', i, y(i)
                  return
                  stop
               end if               
               y(i) = min(1.0d0, max(y(i),1.0d-30))
               x(i) = y(i)*aion(i)
               xsum = xsum + x(i)
            end do
            if (trace .and. xsum > 2) write(*,*) 'sum_x, time', xsum, time
         
            lgT = y(nvar)/ln10
            T = exp10(lgT)
                
            call basic_composition_info( &
               species, g% chem_id, x, xh, xhe, z, &
               abar, zbar, z2bar, z53bar, ye, mass_correction, sumx)
               
            call get_eos_info_for_burn_at_const_density( &
               eos_handle, species, g% chem_id, g% net_iso, x, &
               Rho, lgRho, T, lgT, &
               Cv, d_Cv_dlnT, eta, d_eta_dlnT, ierr)
            if (ierr /= 0) then
               !write(*,*) 'failed in get_eos_info_for_burn_at_const_density'
               return
               do j=1,species
                  write(*,2) 'x', j, x(j)
               end do
               write(*,1) 'sum(x)', sum(x(1:species))
               write(*,1) 'T', T
               write(*,1) 'lgT', lgT
               write(*,1) 'rho', rho
               write(*,1) 'lgRho', lgRho
               write(*,1) 'z', z
               write(*,1) 'xh', xh
               write(*,1) 'abar', abar
               write(*,1) 'zbar', zbar
               stop 'net burn const density'
            end if
            
            rates_only = .false.
            dxdt_only = (size(dfdy,dim=1) == 0)
            
            !write(*,1) 'eval_net lgT', lgT
            
            d_eta_dlnRho = 0d0
            call eval_net( &
               n, g, rates_only, dxdt_only, &
               species, num_reactions, g% num_wk_reactions, &
               x, T, lgT, rho, lgRho, &
               abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
               rate_factors, weak_rate_factor, &
               reaction_Qs, reaction_neuQs, &
               .false., .false., &
               eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
               dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx,  &
               screening_mode,  &
               eps_nuc_categories, eps_neu_total, &
               net_lwork, net_work_array, &
               actual_Qs, actual_neuQs, from_weaklib, .false., ierr)
            if (ierr /= 0) then
               !write(*,*) 'failed in eval_net'
               return
               do j=1,species
                  write(*,2) 'x', j, x(j)
               end do
               write(*,1) 'sum(x)', sum(x(1:species))
               write(*,1) 'T', T
               write(*,1) 'lgT', lgT
               write(*,1) 'rho', rho
               write(*,1) 'lgRho', lgRho
               write(*,1) 'abar', abar
               write(*,1) 'zbar', zbar
               write(*,1) 'z2bar', z2bar
               write(*,1) 'eta', eta
               write(*,1) 'd_eta_dlnT', d_eta_dlnT
               stop 'net burn const density'
            end if
            
            if (size(f,dim=1) > 0) then
               do j = 1, species
                  f(j) = dxdt(j)/aion(j)
               end do
               f(nvar) = eps_nuc/(Cv*T) ! dlnT_dt
               !f(nvar) = 0d0 ! TESTING
            end if

            if (.not. dxdt_only) then
               do j = 1, species
                  do i = 1, species
                     dfdy(i,j) = d_dxdt_dx(i,j)*aion(j)/aion(i)
                     if (.false. .and. is_bad(dfdy(i,j))) then
                        write(*,1) 'x', x
                        write(*,3) 'dfdy(i,j)', i, j, dfdy(i,j)
                        stop 'jakob_or_derivs'
                     end if
                  end do
                  !dfdy(j,nvar) = T*d_dxdt_dT(j) ! d_dxdt_dlnT
                  dfdy(j,nvar) = 0d0 ! TESTING
                  !if (j /= 22) dfdy(j,nvar) = 0d0 ! TESTING
                  !if (j==5 .or. j==13 .or. j==14 .or. j==19 .or. j==20 .or. j==22) dfdy(j,nvar) = 0d0 ! TESTING
                  !if (j== 11) dfdy(j,nvar) = 0d0 ! TESTING
               end do
               do j = 1,species
                  dfdy(nvar,j) = d_eps_nuc_dx(j)/(Cv*T) ! d_lnT_dx(j)
                  !dfdy(nvar,j) = 0d0 ! TESTING
               end do               
               dfdy(nvar,nvar) = &
                  d_eps_nuc_dT/Cv - (1d0 + d_Cv_dlnT/Cv)*eps_nuc/(Cv*T)
               !dfdy(nvar,nvar) = -1d0 ! TESTING
               if (is_bad(dfdy(nvar,nvar))) then
                  ierr = -1
                  return
                  write(*,1) 'dfdy(nvar,nvar)', dfdy(nvar,nvar)
                  write(*,1) 'd_eps_nuc_dT', d_eps_nuc_dT
                  write(*,1) 'Cv', Cv
                  stop 'jakob_or_derivs'
               end if
            end if
         
         end subroutine jakob_or_derivs


      end subroutine burn_const_density_1_zone
      





      end module net_burn_const_density

