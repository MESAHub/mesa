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

      module net_burn_const_P
      use const_def
      use chem_def
      use math_lib
      use net_def
      use rates_def, only: num_rvs
      use mtx_def
      
      implicit none

      contains

      subroutine burn_1_zone_const_P( &
            net_handle, eos_handle, num_isos, num_reactions, &
            which_solver, starting_temp, starting_x, clip, &
            ntimes, times, log10Ps_f1, &
            rate_factors, weak_rate_factor, &
            reaction_Qs, reaction_neuQs, screening_mode, &
            h, max_step_size, max_steps, rtol, atol, itol, x_min, x_max, which_decsol, &
            caller_id, solout, iout, &
            ending_x, ending_temp, ending_rho, ending_lnS, initial_rho, initial_lnS, &
            nfcn, njac, nstep, naccpt, nrejct, time_doing_net, time_doing_eos, ierr)
         use num_def
         use num_lib 
         use mtx_lib
         use rates_def, only: rates_reaction_id_max
         use net_initialize, only: work_size
         
         integer, intent(in) :: net_handle, eos_handle
         integer, intent(in) :: num_isos
         integer, intent(in) :: num_reactions
         real(dp), pointer, intent(in) :: starting_x(:) ! (num_isos)
         real(dp), intent(in) :: starting_temp
         logical, intent(in) :: clip ! if true, set negative x's to zero during burn.
         
         integer, intent(in) :: which_solver ! as defined in num_def.f
         integer, intent(in) :: ntimes ! ending time is times(num_times); starting time is 0
         real(dp), pointer, intent(in) :: times(:) ! (num_times) 
         real(dp), pointer, intent(in) :: log10Ps_f1(:) ! =(4,numtimes) interpolant for log10P(time)

         real(dp), intent(in) :: rate_factors(:) ! (num_reactions)
         real(dp), intent(in) :: weak_rate_factor
         real(dp), pointer, intent(in) :: reaction_Qs(:) ! (rates_reaction_id_max)
         real(dp), pointer, intent(in) :: reaction_neuQs(:) ! (rates_reaction_id_max)
         integer, intent(in) :: screening_mode
         
         ! args to control the solver -- see num/public/num_isolve.dek
         real(dp), intent(inout) :: h 
         real(dp), intent(in) :: max_step_size ! maximal step size.
         integer, intent(in) :: max_steps ! maximal number of allowed steps.
         ! absolute and relative error tolerances
         real(dp), intent(inout) :: rtol(*) ! relative error tolerance(s)
         real(dp), intent(inout) :: atol(*) ! absolute error tolerance(s)
         integer, intent(in) :: itol ! switch for rtol and atol
         real(dp), intent(in) :: x_min, x_max ! bounds on allowed values
         integer, intent(in) :: which_decsol ! from mtx_def
         integer, intent(in) :: caller_id
         interface ! subroutine called after each successful step
            include "num_solout.dek"
         end interface
         integer, intent(in)  :: iout
         real(dp), intent(inout) :: ending_x(:)
         real(dp), intent(out) :: ending_temp, ending_rho, ending_lnS, initial_rho, initial_lnS
         integer, intent(out) :: nfcn    ! number of function evaluations
         integer, intent(out) :: njac    ! number of jacobian evaluations
         integer, intent(out) :: nstep   ! number of computed steps
         integer, intent(out) :: naccpt  ! number of accepted steps
         integer, intent(out) :: nrejct  ! number of rejected steps
         real(dp), intent(inout) :: time_doing_net
            ! if < 0, then ignore
            ! else on return has input value plus time spent doing eval_net
         real(dp), intent(inout) :: time_doing_eos
            ! if < 0, then ignore
            ! else on return has input value plus time spent doing eos
         integer, intent(out) :: ierr
         
         type (Net_Info) :: n
         type (Net_General_Info), pointer :: g
         integer :: ijac, nzmax, isparse, mljac, mujac, imas, mlmas, mumas, lrd, lid, &
               lout, liwork, lwork, i, j, lrpar, lipar, idid, nvar
         integer, pointer :: ipar(:), iwork(:), ipar_burn_const_P_decsol(:)
         real(dp), pointer :: rpar(:), work(:), v(:), rpar_burn_const_P_decsol(:)
         real(dp) :: t, lgT, lgRho, tend
         
         include 'formats'
                  
         ending_x = 0
         ending_temp = 0
         ending_rho = 0
         ending_lnS = 0
         nfcn = 0
         njac = 0
         nstep = 0
         naccpt = 0
         nrejct = 0
         ierr = 0
         
         nvar = num_isos + 1

         call get_net_ptr(net_handle, g, ierr)
         if (ierr /= 0) then
            return
         end if
         
         if (g% num_isos /= num_isos) then
            write(*,*) 'invalid num_isos', num_isos
            return
         end if
         
         if (g% num_reactions /= num_reactions) then
            write(*,*) 'invalid num_reactions', num_reactions
            return
         end if
         
         if (which_decsol == lapack) then
            nzmax = 0
            isparse = 0
            call lapack_work_sizes(nvar, lrd, lid)
         else
            write(*,1) 'net 1 zone burn const P: unknown value for which_decsol', which_decsol
            call mesa_error(__FILE__,__LINE__)
         end if
         
         ijac = 1
         mljac = nvar ! square matrix
         mujac = nvar

         imas = 0
         mlmas = 0
         mumas = 0        
         
         lout = 0
                  
         call isolve_work_sizes(nvar, nzmax, imas, mljac, mujac, mlmas, mumas, liwork, lwork)

         lipar = burn_lipar
         lrpar = burn_const_P_lrpar 
         
         allocate(v(nvar), iwork(liwork), work(lwork), rpar(lrpar), ipar(lipar), &
               ipar_burn_const_P_decsol(lid), rpar_burn_const_P_decsol(lrd), stat=ierr)
         if (ierr /= 0) then
            write(*, *) 'allocate ierr', ierr
            return
         end if
         
         i = burn_const_P_lrpar
            
         ipar(i_burn_caller_id) = caller_id
         ipar(i_net_handle) = net_handle
         ipar(i_eos_handle) = eos_handle
         ipar(i_screening_mode) = screening_mode
         ipar(i_sparse_format) = isparse
         if (clip) then
            ipar(i_clip) = 1
         else
            ipar(i_clip) = 0
         end if

         iwork = 0
         work = 0
         
         t = 0
         tend = times(ntimes)
         
         rpar(r_burn_const_P_pressure) = exp10(log10Ps_f1(1)) ! no interpolation yet
         rpar(r_burn_const_P_temperature) = starting_temp
         rpar(r_burn_const_P_init_rho) = -1d99
         rpar(r_burn_const_P_init_lnS) = -1d99
         rpar(r_burn_const_P_time_net) = time_doing_net
         rpar(r_burn_const_P_time_eos) = time_doing_eos

         v(1:num_isos) = starting_x(1:num_isos)
         v(nvar) = log(starting_temp)
                     
         if (which_decsol == lapack) then
            call do_isolve(lapack_decsol, null_decsols, ierr)
         else
            write(*,*) 'unknown value for which_decsol', which_decsol
            call mesa_error(__FILE__,__LINE__)
         end if
         if (ierr /= 0) then
            call dealloc
            return
         end if

         nfcn = iwork(14)
         njac = iwork(15)
         nstep = iwork(16)
         naccpt = iwork(17)
         nrejct = iwork(18)            
         time_doing_net = rpar(r_burn_const_P_time_net)
         time_doing_eos = rpar(r_burn_const_P_time_eos)

         ending_x(1:num_isos) = v(1:num_isos)
         ending_temp = exp(v(nvar))
         ending_rho = rpar(r_burn_const_P_rho)
         ending_lnS = rpar(r_burn_const_P_lnS)
         initial_rho = rpar(r_burn_const_P_init_rho)
         initial_lnS = rpar(r_burn_const_P_init_lnS)
         
         if (ierr /= 0) then
            write(*, '(a30,i10)') 'nfcn', nfcn
            write(*, '(a30,i10)') 'njac', njac
            write(*, '(a30,i10)') 'nstep', nstep
            write(*, '(a30,i10)') 'naccpt', naccpt
            write(*, '(a30,i10)') 'nrejct', nrejct
            call mesa_error(__FILE__,__LINE__)
         end if
         
         call dealloc
      
         
         contains
         
         
         subroutine dealloc
            deallocate(iwork, work, rpar, ipar, ipar_burn_const_P_decsol, rpar_burn_const_P_decsol)
         end subroutine dealloc
         
         
         subroutine do_isolve(decsol, decsols, ierr)
            interface
               include "mtx_decsol.dek"
               include "mtx_decsols.dek"
            end interface
            integer, intent(out) :: ierr
            integer :: caller_id, nvar_blk, nz_blk
            real(dp), dimension(:), pointer :: &
               lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk
         
            nullify(lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk)
            caller_id = 0
            nvar_blk = 0
            nz_blk = 0
            include 'formats'
            ierr = 0
            ! NOTE: don't use x_max for const_P since x includes other variables
            call isolve( &
               which_solver, nvar, burn_derivs, t, v, tend, &
               h, max_step_size, max_steps, &
               rtol, atol, itol, x_min, 1d99, &
               burn_jacob, ijac, null_sjac, nzmax, isparse, mljac, mujac, &
               null_mas, imas, mlmas, mumas, &
               solout, iout, &
               decsol, decsols, null_decsolblk, &
               lrd, rpar_burn_const_P_decsol, lid, ipar_burn_const_P_decsol, &
               caller_id, nvar_blk, nz_blk, &
               lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk, & 
               null_fcn_blk_dble, null_jac_blk_dble, &
               work, lwork, iwork, liwork, &
               lrpar, rpar, lipar, ipar, &
               lout, idid)
            if (idid < 0) ierr = -1
         end subroutine do_isolve


         subroutine burn_derivs(nvar, t, h, v, f, lrpar, rpar, lipar, ipar, ierr)
            integer, intent(in) :: nvar, lrpar, lipar
            real(dp), intent(in) :: t, h
            real(dp), intent(inout) :: v(:) ! (nvar)
            real(dp), intent(inout) :: f(:) ! (nvar) ! dvdt
            integer, intent(inout), pointer :: ipar(:) ! (lipar)
            real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
            integer, intent(out) :: ierr
            integer, parameter :: ld_dfdv = 0
            real(dp) :: dfdv(ld_dfdv,nvar)
            ierr = 0
            call burn_jacob(nvar, t, h, v, f, dfdv, ld_dfdv, lrpar, rpar, lipar, ipar, ierr)
         end subroutine burn_derivs


         subroutine burn_jacob( &
               nvar, time, h, v, f, dfdv, ld_dfdv, lrpar, rpar, lipar, ipar, ierr)
            use chem_lib, only: basic_composition_info
            use net_eval, only: eval_net
            use eos_def
            use eos_lib, only: Radiation_Pressure, eosPT_get
            use rates_def, only: rates_reaction_id_max

            integer, intent(in) :: nvar, ld_dfdv, lrpar, lipar
            real(dp), intent(in) :: time, h
            real(dp), intent(inout) :: v(:)
            real(dp), intent(inout) :: f(:), dfdv(:,:)
            integer, intent(inout), pointer :: ipar(:) ! (lipar)
            real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
            integer, intent(out) :: ierr
         
            integer :: net_handle, num_reactions, eos_handle
            real(dp) :: &
                  abar, zbar, z2bar, z53bar, ye, mass_correction, sumx, T, logT, &
                  rho, logRho, pressure, Pgas, Prad, lgPgas, &
                  dlnT_dt, x(nvar-1)
            real(dp) :: eta, d_eta_dlnT, d_eta_dlnRho
            real(dp) :: eps_neu_total, eps_nuc
            real(dp) :: d_eps_nuc_dT
            real(dp) :: d_eps_nuc_dRho
            real(dp) :: d_eps_nuc_dx(nvar-1) 
            real(dp) :: dxdt(nvar-1)
            real(dp) :: d_dxdt_dRho(nvar-1)
            real(dp) :: d_dxdt_dT(nvar-1)
            real(dp) :: d_dxdt_dx(nvar-1, nvar-1)
            real(dp), target :: eps_nuc_categories(num_categories)
            logical :: rates_only, skip_jacobian
            integer :: screening_mode, i, num_isos
            integer(8) :: time0, time1, clock_rate

            real(dp) :: xh, Y, z, Cp, rate_limit
            real(dp) :: dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas
            real(dp) :: dlnRho_dlnT_const_P, d_epsnuc_dlnT_const_P, d_Cp_dlnT
            real(dp) :: res(num_eos_basic_results)
            real(dp) :: d_dlnRho_const_T(num_eos_basic_results) 
            real(dp) :: d_dlnT_const_Rho(num_eos_basic_results) 
            real(dp) :: d_dxa_const_TRho(num_eos_d_dxa_results, nvar-1) 
            integer, pointer :: net_iso(:), chem_id(:)

            type (Net_General_Info), pointer :: g
            real(dp), pointer, dimension(:) :: actual_Qs, actual_neuQs
            logical, pointer :: from_weaklib(:)
            actual_Qs => null()
            actual_neuQs => null()
            from_weaklib => null()
         
            include 'formats'
         
            num_isos = nvar-1
         
            ierr = 0
            f = 0
            dfdv = 0
         
            eos_handle = ipar(i_eos_handle)
         
            net_handle = ipar(i_net_handle)
            call get_net_ptr(net_handle, g, ierr)
            if (ierr /= 0) then
               write(*,*) 'invalid handle for eval_net -- did you call alloc_net_handle?'
               return
            end if
         
            v(1:num_isos) = max(1d-30, min(1d0, v(1:num_isos)))
               ! positive definite mass fractions
            v(1:num_isos) = v(1:num_isos)/sum(v(1:num_isos))
            x(1:num_isos) = v(1:num_isos)
         
            num_reactions = g% num_reactions

            i = burn_const_P_lrpar
         
            if (ipar(i_clip) /= 0) then
               do i=1,num_isos
                  x(i) = max(0d0, min(1d0, x(i)))
               end do
            end if

            call basic_composition_info( &
               num_isos, g% chem_id, x, xh, Y, z, abar, zbar, z2bar, z53bar, ye, mass_correction, sumx)
     
            logT = v(nvar)/ln10
            T = exp10(logT)

            pressure = rpar(r_burn_const_P_pressure)         
            Prad = Radiation_Pressure(T)
            Pgas = pressure - Prad
            if (Pgas <= 0) then
               write(*,1) 'Pgas <= 0 in burn at const P', Pgas
               ierr = -1
               return
            end if
            lgPgas = log10(Pgas)

            chem_id => g% chem_id
            net_iso => g% net_iso
                  
            if (rpar(r_burn_const_P_time_eos) >= 0) then
               call system_clock(time0,clock_rate)
            else
               time0 = 0
            endif
         
            call eosPT_get( &
               eos_handle, &
               num_isos, chem_id, net_iso, x, &
               Pgas, lgPgas, T, logT, &
               Rho, logRho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
               res, d_dlnRho_const_T, d_dlnT_const_Rho,  &
               d_dxa_const_TRho, ierr)

            Cp = res(i_Cp)
            eta = res(i_eta)
            d_eta_dlnRho = d_dlnRho_const_T(i_eta)
            d_eta_dlnT = d_dlnT_const_Rho(i_eta)
            screening_mode = ipar(i_screening_mode)
            rpar(r_burn_const_P_rho) = Rho
            rpar(r_burn_const_P_temperature) = T
            rpar(r_burn_const_P_lnS) = res(i_lnS)
            if (rpar(r_burn_const_P_init_rho) < -1d90) &
               rpar(r_burn_const_P_init_rho) = Rho
            if (rpar(r_burn_const_P_init_lnS) < -1d90) &
               rpar(r_burn_const_P_init_lnS) = res(i_lnS)
            if (ierr /= 0 .or. Cp <= 0) then
               ierr = -1
               return

               write(*,*) 'eosPT_get failed'
               write(*,1) 'xh', xh
               write(*,1) 'Y', Y
               write(*,1) 'Z', 1 - (xh + Y)
               write(*,1) 'abar', abar
               write(*,1) 'zbar', zbar
               write(*,1) 'pressure', pressure
               write(*,1) 'Prad', Prad
               write(*,1) 'Pgas', Pgas
               write(*,1) 'lgPgas', lgPgas
               write(*,1) 'T', T
               write(*,1) 'logT', logT
               write(*,1) 'Rho', Rho
               write(*,1) 'logRho', logRho
               write(*,1) 'Cp', Cp
               return
            end if

            if (rpar(r_burn_const_P_time_eos) >= 0) then
               call system_clock(time1,clock_rate)
               rpar(r_burn_const_P_time_eos) = &
                  rpar(r_burn_const_P_time_eos) + dble(time1 - time0) / clock_rate
               if (rpar(r_burn_const_P_time_net) >= 0) time0 = time1
            else if (rpar(r_burn_const_P_time_net) >= 0) then
               call system_clock(time0,clock_rate)
            end if

            rates_only = .false.
            skip_jacobian = .false.

            call eval_net( &
                  n, g, rates_only, skip_jacobian, &
                  num_isos, num_reactions, g% num_wk_reactions, &
                  x, T, logT, rho, logRho, &
                  abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
                  rate_factors, weak_rate_factor, &
                  reaction_Qs, reaction_neuQs, &
                  eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx, &
                  dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx, &
                  screening_mode, &
                  eps_nuc_categories, eps_neu_total, &
                  actual_Qs, actual_neuQs, from_weaklib, .false., ierr)

            if (rpar(r_burn_const_P_time_net) >= 0) then
               call system_clock(time1,clock_rate)
               rpar(r_burn_const_P_time_net) = &
                  rpar(r_burn_const_P_time_net) + dble(time1 - time0) / clock_rate
            end if

            if (ierr /= 0) then
               return
            
            
            
               write(*,*) 'eval_net failed'
               write(*,1) 'xh', xh
               write(*,1) 'Y', Y
               write(*,1) 'Z', 1 - (xh + Y)
               write(*,1) 'abar', abar
               write(*,1) 'zbar', zbar
               write(*,1) 'pressure', pressure
               write(*,1) 'Prad', Prad
               write(*,1) 'Pgas', Pgas
               write(*,1) 'lgPgas', lgPgas
               write(*,1) 'T', T
               write(*,1) 'logT', logT
               write(*,1) 'Rho', Rho
               write(*,1) 'logRho', logRho
               write(*,1) 'Cp', Cp
               ierr = -1
            
            
               call mesa_error(__FILE__,__LINE__,'net_burn_const_P')
            
               return
            end if
         
            f(1:num_isos) = dxdt
            dlnT_dt = eps_nuc/(Cp*T)
            f(nvar) = dlnT_dt
         
            if (ld_dfdv > 0) then

               dlnRho_dlnT_const_P = -res(i_chiT)/res(i_chiRho)
               d_epsnuc_dlnT_const_P = d_eps_nuc_dT*T + d_eps_nuc_dRho*Rho*dlnRho_dlnT_const_P
               d_Cp_dlnT = d_dlnT_const_Rho(i_Cp) + d_dlnRho_const_T(i_Cp)*dlnRho_dlnT_const_P
            
               dfdv(1:num_isos,1:num_isos) = d_dxdt_dx

               dfdv(nvar,nvar) = d_epsnuc_dlnT_const_P/(Cp*T) - dlnT_dt*(1 + d_Cp_dlnT/Cp)
            
               ! d_dxdt_dlnT
               dfdv(1:num_isos,nvar) = &
                  d_dxdt_dT(1:num_isos)*T + d_dxdt_dRho(1:num_isos)*Rho*dlnRho_dlnT_const_P
            
               ! d_dlnTdt_dx
               dfdv(nvar,1:num_isos) = d_eps_nuc_dx(1:num_isos)/(Cp*T)

            end if
         
         end subroutine burn_jacob


         subroutine burn_sjac(n,time,h,y,f,nzmax,ia,ja,values,lrpar,rpar,lipar,ipar,ierr)  
            use mtx_lib, only: dense_to_sparse_with_diag
            integer, intent(in) :: n, nzmax, lrpar, lipar
            real(dp), intent(in) :: time, h
            real(dp), intent(inout) :: y(n)
            integer, intent(out) :: ia(:), ja(:)
            real(dp), intent(inout) :: f(:), values(:)
            integer, intent(inout), pointer :: ipar(:) ! (lipar)
            real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
            integer, intent(out) :: ierr ! nonzero means terminate integration
            real(dp), pointer :: dfdv(:,:) ! (n,n)
            integer :: ld_dfdv, nz, i, j, cnt, nnz
            include 'formats'
            !write(*,1) 'burn_sjac', x
            ierr = 0
            ld_dfdv = n
            allocate(dfdv(n,n),stat=ierr)
            if (ierr /= 0) return
            call burn_jacob(n,time,h,y,f,dfdv,ld_dfdv,lrpar,rpar,lipar,ipar,ierr)
            if (ierr /= 0) then
               deallocate(dfdv)
               return
            end if
            ! remove entries with abs(value) < 1d-16
            cnt = 0; nnz = 0
            do i=1,n
               do j=1,n
                  if (dfdv(i,j) /= 0) then
                     nnz = nnz + 1
                     if (abs(dfdv(i,j)) < 1d-16) then
                        cnt = cnt+1; dfdv(i,j) = 0
                     end if
                  end if
               end do
            end do
            call dense_to_sparse_with_diag( &
               ipar(i_sparse_format),n,n,dfdv,nzmax,nz,ia,ja,values,ierr)
            deallocate(dfdv)
         end subroutine burn_sjac
      
         
      end subroutine burn_1_zone_const_P
      

      end module net_burn_const_P

