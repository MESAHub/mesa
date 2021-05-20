! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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

      module net_burn
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
      
      !logical, parameter :: report_ierr = .true.
      logical, parameter :: report_ierr = .false.


      

      contains
      
      
      integer function burn_1_zone_work_size(g) result(sz)
         use net_initialize, only: work_size
         use net_approx21, only: approx21_nrat
         type (Net_General_Info), pointer :: g
         integer :: net_lwork, num_reactions, species
         num_reactions = g% num_reactions
         if (g% doing_approx21) num_reactions = approx21_nrat
         species = g% num_isos
         net_lwork = work_size(g)
         sz = net_lwork + 2*num_reactions + 2*species*species
         ! dratdumdy1, dratdumdy2, dens_dfdy1, dmat1
      end function burn_1_zone_work_size
         

      subroutine get_pointers( &
            g, burn_lwork, burn_work_array, species, num_reactions, &
            dratdumdy1, dratdumdy2, dens_dfdy, dmat, i, ierr)
         use net_approx21, only: approx21_nrat
         type (Net_General_Info), pointer :: g
         real(dp), pointer :: burn_work_array(:)
         integer, intent(in) :: burn_lwork, species, num_reactions
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
         
         sz = species*species
         dens_dfdy(1:species,1:species) => burn_work_array(i+1:i+sz); i = i+sz
         dmat(1:species,1:species) => burn_work_array(i+1:i+sz); i = i+sz
         
      end subroutine get_pointers

      
      subroutine burn_1_zone( &
            net_handle, eos_handle, species, num_reactions, t_start, t_end, starting_x, &
            ntimes, times, log10Ts_f1, log10Rhos_f1, etas_f1, dxdt_source_term, &
            rate_factors, weak_rate_factor, &
            reaction_Qs, reaction_neuQs, &
            screening_mode,  &
            stptry_in, max_steps, eps, odescal, &
            use_pivoting, trace, burn_dbg, burner_finish_substep, &
            burn_lwork, burn_work_array, &
            net_lwork, net_work_array, &
            ending_x, eps_nuc_categories, avg_eps_nuc, eps_neu_total, &
            nfcn, njac, nstep, naccpt, nrejct, ierr)
         use num_def
         use num_lib 
         use mtx_lib
         use mtx_def
         use rates_def, only: rates_reaction_id_max, reaction_Name, reaction_categories
         use rates_lib, only: rates_reaction_id
         use net_initialize, only: set_rate_ptrs, setup_net_info, work_size
         use chem_lib, only: basic_composition_info, get_Q
         use net_initialize, only: work_size
         use net_approx21, only: approx21_nrat
         
         integer, intent(in) :: net_handle, eos_handle
         integer, intent(in) :: species
         integer, intent(in) :: num_reactions
         real(dp), intent(in) :: t_start, t_end, starting_x(:) ! (species)
         integer, intent(in) :: ntimes ! ending time is times(num_times); starting time is 0
         real(dp), pointer, intent(in) :: times(:) ! (num_times) 
         real(dp), pointer, intent(in) :: log10Ts_f1(:) 
            ! =(4,numtimes) interpolant for log10T(time)
         real(dp), pointer, intent(in) :: log10Rhos_f1(:)
            ! =(4,numtimes) interpolant for log10Rho(time)
         real(dp), pointer, intent(in) :: etas_f1(:)
            ! =(4,numtimes) interpolant for eta(time)
         real(dp), pointer, intent(in) :: dxdt_source_term(:)
            ! (species)  or null if no source term.
         real(dp), intent(in), pointer :: rate_factors(:) ! (num_reactions)
         real(dp), intent(in) :: weak_rate_factor
         real(dp), pointer, intent(in) :: reaction_Qs(:) ! (rates_reaction_id_max)
         real(dp), pointer, intent(in) :: reaction_neuQs(:) ! (rates_reaction_id_max)
         integer, intent(in) :: screening_mode ! see screen_def
         real(dp), intent(in) :: stptry_in
         integer, intent(in) :: max_steps ! maximal number of allowed steps.
         real(dp), intent(in) :: eps, odescal ! tolerances.  e.g., set both to 1d-6
         logical, intent(in) :: use_pivoting
         logical, intent(in) :: trace, burn_dbg
         interface
            include 'burner_finish_substep.inc'
         end interface
         integer, intent(in) :: burn_lwork, net_lwork
         real(dp), intent(inout), pointer :: burn_work_array(:)
         real(dp), intent(inout), pointer :: net_work_array(:)
         real(dp), intent(inout) :: ending_x(:) ! (species)
         real(dp), intent(inout) :: eps_nuc_categories(:) ! (num_categories)
         real(dp), intent(out) :: avg_eps_nuc, eps_neu_total
         integer, intent(out) :: nfcn    ! number of function evaluations
         integer, intent(out) :: njac    ! number of jacobian evaluations
         integer, intent(out) :: nstep   ! number of computed steps
         integer, intent(out) :: naccpt  ! number of accepted steps
         integer, intent(out) :: nrejct  ! number of rejected steps
         integer, intent(out) :: ierr
         
         type (Net_General_Info), pointer :: g
         integer :: ijac, lrd, lid, lout, i, j, ir, idid, sz
         logical :: okay, have_set_rate_screened
         real(dp) :: temp, rho, eta, lgT, lgRho, r, prev_lgRho, prev_lgT
         
         integer :: stpmax, imax_dydx, nstp
         real(dp) :: &
            h, start, stptry, stpmin, stopp, max_dydx, abs_max_dydx, &
            burn_ergs, dx
                  
         real(dp), dimension(species), target :: starting_y_a, ending_y_a, save_x_a
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
         end if

         starting_y => starting_y_a
         ending_y => ending_y_a
         save_x => save_x_a
         n => net_info_target
         
         have_set_rate_screened = .false.
         
         lgT = log10Ts_f1(1)
         temp = exp10(lgT)
         lgRho = log10Rhos_f1(1)
         rho = exp10(lgRho)
         eta = etas_f1(1)
         prev_lgT = -1d99
         prev_lgRho = -1d99
         
         ierr = 0
         call get_net_ptr(net_handle, g, ierr)
         if (ierr /= 0) then
            if (report_ierr) write(*,*) 'invalid handle for burn_1_zone'
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
            iwork, ierr)
         if (dbg) write(*,*) 'done setup_net_info'
         if (ierr /= 0) then
            if (dbg) write(*,*) 'failed in setup_net_info'
            return
         end if
      
         if (dbg) write(*,*) 'call get_pointers'
         call get_pointers( &
            g, burn_lwork, burn_work_array, species, num_reactions, &
            dratdumdy1, dratdumdy2, dens_dfdy, dmat, iwork, ierr)
         if (ierr /= 0) return

         call basic_composition_info( &
            species, g% chem_id, starting_x, xh, xhe, z, &
            abar, zbar, z2bar, z53bar, ye, mass_correction, sumx)

         stptry = max(start * 1.0d-10,1.0d-16)
         stpmin = stptry * 1.0d-12

         if (dbg) write(*,*) 'call netint'
         call netint( &
            start,stptry,stpmin,max_steps,stopp,ending_y, &
            eps,species,species,naccpt,nrejct,nstep,odescal,dens_dfdy,dmat, &
            burner_derivs,burner_jakob,burner_finish_substep,ierr)
         if (dbg) write(*,*) 'done netint'
         if (ierr /= 0) then
            return
            !write(*,*) 'netint ierr'
            !stop
         end if

         do i=1,species
            ending_x(i) = ending_y(i)*aion(i)
         end do

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
         !write(*,1) 'burn_ergs', burn_ergs
         avg_eps_nuc = burn_ergs/(t_end - t_start) - eps_neu_total
      
      contains
         
         subroutine burner_derivs(x,y,f,species,ierr)
            integer, intent(in) :: species
            real(dp) :: x, y(:), f(:)
            integer, intent(out) :: ierr
            integer, parameter :: ld_dfdx = 0
            real(dp), target :: dfdx_arry(ld_dfdx,species)
            real(dp), pointer :: dfdx(:,:)
            real(dp) :: dxdt_sum, dxdt_sum_aprox21, &
               Z_plus_N, xsum, r, r1, r2
            integer :: i, ir, ci, j, k, ibad
            logical :: okay

            real(dp), target :: f21_a(species)
            real(dp), pointer :: f21(:)
            
            include 'formats'
            
            ierr = 0
            nfcn = nfcn + 1
            dfdx => dfdx_arry
            call jakob_or_derivs(x,y,f,dfdx,ierr)
            if (ierr /= 0) return            
         
         end subroutine burner_derivs

         subroutine burner_jakob(x,y,dfdy,species,ierr)
            integer, intent(in) :: species
            real(dp) :: x, y(:)
            real(dp), pointer :: dfdy(:,:)
            integer, intent(out) :: ierr
            real(dp), target :: f_arry(0)
            real(dp), pointer :: f(:)
            
            real(dp), target :: dfdy21_a(species,species)
            real(dp), pointer :: dfdy21(:,:)
            real(dp) :: Z_plus_N, df_t, df_m
            integer :: i, ci, j, cj
            logical :: okay
            include 'formats'

            ierr = 0
            njac = njac + 1
            f => f_arry

            call jakob_or_derivs(x,y,f,dfdy,ierr)
            if (ierr /= 0) return
                     
         end subroutine burner_jakob

         subroutine jakob_or_derivs(time,y,f,dfdy,ierr)
            use chem_lib, only: basic_composition_info
            use chem_def, only: chem_isos, num_categories, category_name, ih1
            use net_eval, only: eval_net
            use rates_def, only: rates_reaction_id_max, i_rate, i_rate_dT, i_rate_dRho
            use interp_1d_lib, only: interp_value
            use eos_def, only: num_eos_basic_results, num_eos_d_dxa_results, i_eta
            use eos_lib, only: eosDT_get
         
            real(dp) :: time, y(:), f(:)
            real(dp), pointer :: dfdy(:,:)
            integer, intent(out) :: ierr
         
            real(dp) :: rho, lgRho, T, lgT, rate_limit, rat, dratdt, dratdd
            real(dp) :: eta, d_eta_dlnT, d_eta_dlnRho
            real(dp) :: eps_nuc
            real(dp) :: d_eps_nuc_dT
            real(dp) :: d_eps_nuc_dRho
            real(dp) :: d_eps_nuc_dx(species) 
            real(dp) :: dxdt(species)
            real(dp) :: d_dxdt_dRho(species)
            real(dp) :: d_dxdt_dT(species)
            real(dp) :: d_dxdt_dx(species, species)
            
            logical :: rates_only, dxdt_only, okay
            integer :: i, j, k, ir

            real(dp), pointer, dimension(:) :: actual_Qs, actual_neuQs
            real(dp) :: xsum
            logical, pointer :: from_weaklib(:)
            
            real(dp), target :: x_a(species), dfdx_a(species,species)
            real(dp), pointer :: x(:), dfdx(:,:)

            real(dp), dimension(num_eos_basic_results) :: res, d_dlnd, d_dlnT
            real(dp) :: d_dxa(num_eos_d_dxa_results,species)
         
            include 'formats'
         
            ierr = 0

            x => x_a
            dfdx => dfdx_a
            
            actual_Qs => null()
            actual_neuQs => null()
            from_weaklib => null()
         
            if (ntimes == 1) then
         
               lgT = log10Ts_f1(1)
               lgRho = log10Rhos_f1(1)
            
            else
         
               call interp_value(times, ntimes, log10Ts_f1, time, lgT, ierr)
               if (ierr /= 0) then
                  if (report_ierr) &
                     write(*,1) 'interp_value for lgT failed in jakob for 1 zone burn', time
                  return
               end if

               call interp_value(times, ntimes, log10Rhos_f1, time, lgRho, ierr)
               if (ierr /= 0) then
                  if (report_ierr) &
                     write(*,1) 'interp_value for lgRho failed in jakob for 1 zone burn', time
                  return
               end if

            end if
            
            xsum = 0 
            do i=1,species
               if (is_bad(y(i))) then
                  ierr = -1
                  if (report_ierr) &
                     write(*,2) 'net_burn failed in jakob_or_derivs: bad y(i) lgT lgRho', i, y(i), lgT, lgRho
                  return
                  stop
               end if               
               y(i) = min(1.0d0, max(y(i),1.0d-30))
               x(i) = y(i)*aion(i)
               xsum = xsum + x(i)
            end do
            if (trace .and. xsum > 2) write(*,*) 'sum_x, time', xsum, time

            rho = exp10(lgRho)
            T = exp10(lgT)
                
            call basic_composition_info( &
               species, g% chem_id, x, xh, xhe, z, &
               abar, zbar, z2bar, z53bar, ye, mass_correction, sumx)
            

            call eosDT_get( &
               eos_handle, species, g% chem_id, g% net_iso, x, &
               Rho, lgRho, T, lgT, &
               res, d_dlnd, d_dlnT, d_dxa, ierr)
            if (ierr /= 0) then
               if (report_ierr) write(*,*) 'failed in eos_get_helm_results'
               return
            end if
            eta = res(i_eta)
            d_eta_dlnT = d_dlnT(i_eta)
            d_eta_dlnRho = d_dlnd(i_eta)

            
            rates_only = .false.
            dxdt_only = (size(dfdy,dim=1) == 0)
            
            call eval_net( &
               n, g, rates_only, dxdt_only, &
               species, num_reactions, g% num_wk_reactions, &
               x, T, lgT, rho, lgRho, &
               abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
               rate_factors, weak_rate_factor, &
               reaction_Qs, reaction_neuQs, &
               eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
               dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx,  &
               screening_mode, &
               eps_nuc_categories, eps_neu_total, &
               net_lwork, net_work_array, &
               actual_Qs, actual_neuQs, from_weaklib, .false., ierr)
            
            if (size(f,dim=1) > 0) then
               do j = 1, species
                  f(j) = dxdt(j)/aion(j)
                  if (.false. .and. is_bad(f(j))) then
                     write(*,1) 'x', x
                     write(*,2) 'f(j)', j, f(j)
                     stop 'jakob_or_derivs'
                  end if
               end do
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
               end do
            end if
            
         
         end subroutine jakob_or_derivs
         
         
         

      end subroutine burn_1_zone
      





      end module net_burn

