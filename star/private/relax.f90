! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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

      module relax

      use star_private_def
      use const_def
      use utils_lib

      implicit none


      real(dp), parameter :: min_dlnz = -12
      real(dp) :: min_z = 1d-12

      ! some relax routines depend on things such as other_energy and other_torque
      ! to which interpolation parameters cannot be passed directly. So for simplicity
      ! use these two global variables instead.
      integer :: relax_num_pts
      real(dp), pointer :: relax_work_array(:)

      contains


      subroutine do_relax_to_limit(id, restore_at_end, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restore_at_end
         integer, intent(out) :: ierr

         integer, parameter ::  lipar=1, lrpar=1

         type (star_info), pointer :: s

         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:)
         real(dp), target :: rpar_ary(lrpar)
         real(dp), pointer :: rpar(:)
         rpar => rpar_ary
         ipar => ipar_ary

         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call do_internal_evolve( &
            id, before_relax_to_limit, relax_to_limit_adjust_model, &
            relax_to_limit_check_model, null_finish_model, &
            restore_at_end, lipar, ipar, lrpar, rpar, ierr)
      end subroutine do_relax_to_limit


      subroutine before_relax_to_limit(s, id, lipar, ipar, lrpar, rpar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine before_relax_to_limit

      integer function relax_to_limit_adjust_model(s, id, lipar, ipar, lrpar, rpar)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         relax_to_limit_adjust_model = keep_going
      end function relax_to_limit_adjust_model

      integer function relax_to_limit_check_model(s, id, lipar, ipar, lrpar, rpar)
         use do_one_utils, only: do_bare_bones_check_model, do_check_limits
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         relax_to_limit_check_model = do_bare_bones_check_model(id)
         if (relax_to_limit_check_model == keep_going) then
            relax_to_limit_check_model = do_check_limits(id)
            if (relax_to_limit_check_model == terminate) &
               s% termination_code = t_relax_finished_okay
         end if
      end function relax_to_limit_check_model


      subroutine do_relax_mass(id, new_mass, lg_max_abs_mdot, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: new_mass, lg_max_abs_mdot
         integer, intent(out) :: ierr
         integer, parameter ::  lipar=1, lrpar=3
         integer :: max_model_number
         character (len=32) :: cool_wind_AGB_scheme, cool_wind_RGB_scheme
         real(dp) :: starting_dt_next, varcontrol_target, eps_mdot_factor
         logical :: adding_mass
         type (star_info), pointer :: s
         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:)
         real(dp), target :: rpar_ary(lrpar)
         real(dp), pointer :: rpar(:)
         rpar => rpar_ary
         ipar => ipar_ary
         include 'formats'
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (abs(new_mass - s% star_mass) < 1d-12*new_mass) then
            s% star_mass = new_mass
            s% mstar = new_mass*Msun
            s% xmstar = s% mstar - s% M_center
            return
         end if
         write(*,*)
         write(*,1) 'current mass', s% mstar/Msun
         write(*,1) 'relax to new_mass', new_mass
         write(*,1) 'lg_max_abs_mdot', lg_max_abs_mdot
         write(*,*)
         if (new_mass <= 0) then
            ierr = -1
            write(*,*) 'invalid new mass'
            return
         end if

         rpar(1) = new_mass*Msun
         rpar(2) = lg_max_abs_mdot
         rpar(3) = s% mstar

         adding_mass = (new_mass > s% star_mass)

         cool_wind_AGB_scheme = s% cool_wind_AGB_scheme
         s% cool_wind_AGB_scheme = ''

         cool_wind_RGB_scheme = s% cool_wind_RGB_scheme
         s% cool_wind_RGB_scheme = ''

         varcontrol_target = s% varcontrol_target
         s% varcontrol_target = 1d-3

         eps_mdot_factor = s% eps_mdot_factor
         s% eps_mdot_factor = 0

         max_model_number = s% max_model_number
         s% max_model_number = -1111
         starting_dt_next = s% dt_next
         call do_internal_evolve( &
            id, before_evolve_relax_mass, relax_mass_adjust_model, relax_mass_check_model, &
            null_finish_model, .true., lipar, ipar, lrpar, rpar, ierr)
         s% max_model_number = max_model_number
         s% dt_next = min(s% dt_next, starting_dt_next) * 1d-1
         s% initial_mass = new_mass

         s% cool_wind_AGB_scheme = cool_wind_AGB_scheme
         s% cool_wind_RGB_scheme = cool_wind_RGB_scheme
         s% varcontrol_target = varcontrol_target
         s% eps_mdot_factor = eps_mdot_factor
         s% star_mass = new_mass
         s% mstar = new_mass*Msun
         s% xmstar = s% mstar - s% M_center

         call error_check('relax mass',ierr)

         write(*,2) 's% max_model_number', s% max_model_number


      end subroutine do_relax_mass


      subroutine before_evolve_relax_mass(s, id, lipar, ipar, lrpar, rpar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         ierr = 0
         call setup_before_relax(s)
         s% mass_change = 0
         s% dt_next = min(s% dt_next, 1d4*secyer)
      end subroutine before_evolve_relax_mass

      integer function relax_mass_adjust_model(s, id, lipar, ipar, lrpar, rpar)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         relax_mass_adjust_model = keep_going
      end function relax_mass_adjust_model

      integer function relax_mass_check_model(s, id, lipar, ipar, lrpar, rpar)
         use do_one_utils, only:do_bare_bones_check_model
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer :: ramp
         real(dp) :: &
            new_mass, old_mass, mdot, max_abs_mdot, abs_diff, lg_max_abs_mdot

         logical, parameter :: dbg = .false.

         include 'formats'

         if (dbg) write(*,*) 'relax_mass_check_model'

         relax_mass_check_model = do_bare_bones_check_model(id)
         if (relax_mass_check_model /= keep_going) then
            write(*,*) 'forced termination'
            return
         end if

         new_mass = rpar(1)
         lg_max_abs_mdot = rpar(2)
         if (lg_max_abs_mdot <= -100) then ! use default
            if (s% star_mass < 0.003d0) then
               lg_max_abs_mdot = -7d0
            else if (s% star_mass < 0.006d0) then
               lg_max_abs_mdot = -6.3d0
            else if (s% star_mass < 0.01d0) then
               lg_max_abs_mdot = -6d0
            else if (s% star_mass < 0.1d0) then
               lg_max_abs_mdot = -4d0
            else
               lg_max_abs_mdot = -3d0
            end if
         end if
         max_abs_mdot = exp10(lg_max_abs_mdot)*(Msun/secyer)
         old_mass = rpar(3)

         if (s% model_number >= s% max_model_number .and. s% max_model_number > 0) then
            write(*,3) 's% model_number >= s% max_model_number', &
               s% model_number, s% max_model_number
            relax_mass_check_model = terminate
            return
         end if

         abs_diff = abs(s% mstar - new_mass)
         mdot = (new_mass - s% mstar)/s% dt_next
         if (abs(mdot) > max_abs_mdot) then
            mdot = sign(max_abs_mdot, mdot)
         end if

         if (abs_diff < 1d-4*new_mass .or. abs(mdot) < 1d-50) then
            s% mass_change = 0
            s% star_mass = new_mass
            s% mstar = new_mass*Msun
            s% xmstar = s% mstar - s% M_center
            s% mass_change = 0
            write(*,1) 's% tau_base =', s% tau_base
            write(*,1) 's% tau_factor =', s% tau_factor
            write(*,1) 'atm_option = ' // trim(s% atm_option)
            relax_mass_check_model = terminate
            s% termination_code = t_relax_finished_okay
            return
         end if

         s% max_timestep = abs(new_mass - s% mstar)/mdot

         if (dbg) then
            write(*,1) 'new_mass/Msun', new_mass/Msun
            write(*,1) 'old_mass/Msun', old_mass/Msun
            write(*,1) 'abs_diff/Msun', abs_diff/Msun
         end if

         ramp = 12
         if (s% model_number < ramp) then
            mdot = mdot * pow(1.1d0,dble(s% model_number-ramp))
         end if

         if (abs(mdot)*s% dt > 0.05d0*s% mstar) mdot = sign(0.05d0*s% mstar/s% dt,mdot)

         s% mass_change = mdot/(Msun/secyer)

         if (dbg) write(*,1) 's% mass_change', s% mass_change

      end function relax_mass_check_model


      subroutine do_relax_composition(  &
            id, num_steps_to_use, num_pts, species, xa, xq, ierr)
         integer, intent(in) :: id
         integer, intent(in) :: num_steps_to_use ! use this many steps to do conversion
         integer, intent(in) :: num_pts
            ! length of composition vector; need not equal nz for current model
         integer, intent(in) :: species
            ! per point; must = number of species for current model
         real(dp), intent(in) :: xa(species,num_pts) ! desired composition profile
         real(dp), intent(in) :: xq(num_pts)
         integer, intent(out) :: ierr

         integer, parameter ::  lipar=3
         integer :: lrpar, max_model_number
         real(dp), pointer :: rpar(:)
         real(dp) :: starting_dt_next, mix_factor, dxdt_nuc_factor
         logical :: do_element_diffusion, include_composition_in_eps_grav
         type (star_info), pointer :: s
         real(dp), pointer :: x(:), f1(:), f(:,:,:)
         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:)
         ipar => ipar_ary

         include 'formats'

         ierr = 0

         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         if (species /= s% species) then
            ierr = -1
            return
         end if

         ipar(1) = num_pts
         ipar(2) = num_steps_to_use
         ipar(3) = s% model_number
         lrpar = (1 + 4*species)*num_pts
         allocate(rpar(lrpar), stat=ierr)
         if (ierr /= 0) return

         x(1:num_pts) => rpar(1:num_pts)
         f1(1:4*num_pts*species) => rpar(num_pts+1:lrpar)
         f(1:4,1:num_pts,1:species) => f1(1:4*num_pts*species)

         call store_rpar(species, num_pts, ierr)
         if (ierr /= 0) return

         max_model_number = s% max_model_number
         s% max_model_number = num_steps_to_use + s% model_number + 1
         write(*,*) 'relax_composition: num_steps_to_use', num_steps_to_use

         dxdt_nuc_factor = s% dxdt_nuc_factor
         s% dxdt_nuc_factor = 0 ! turn off composition change by nuclear burning
         mix_factor = s% mix_factor
         s% mix_factor = 0 ! turn off mixing
         do_element_diffusion = s% do_element_diffusion
         s% do_element_diffusion = .false. ! turn off diffusion
         include_composition_in_eps_grav = s% include_composition_in_eps_grav
         s% include_composition_in_eps_grav = .false. ! don't need energetic effects of artificial changes
         starting_dt_next = s% dt_next
         call do_internal_evolve( &
               id, before_evolve_relax_composition, &
               relax_composition_adjust_model, relax_composition_check_model, &
               null_finish_model, .true., lipar, ipar, lrpar, rpar, ierr)

         s% max_model_number = max_model_number
         s% dt_next = starting_dt_next
         s% dxdt_nuc_factor = dxdt_nuc_factor
         s% mix_factor = mix_factor
         s% do_element_diffusion = do_element_diffusion
         s% include_composition_in_eps_grav = include_composition_in_eps_grav
         
         call error_check('relax composition',ierr)

         deallocate(rpar)


         contains


         subroutine store_rpar(species, num_pts, ierr) ! get interpolation info
            use interp_1d_def, only: pm_work_size
            use interp_1d_lib, only: interp_pm
            integer, intent(in) :: species, num_pts
            integer, intent(out) :: ierr
            integer :: j, op_err
            real(dp), pointer :: interp_work(:), work(:), p(:)
            allocate(interp_work(num_pts*pm_work_size*species), stat=ierr)
            if (ierr /= 0) return
            x(:) = xq(:)
            do j=1, species ! make piecewise monotonic cubic interpolants
               op_err = 0
               f(1,1:num_pts,j) = xa(j,1:num_pts)
               work(1:num_pts*pm_work_size) => &
                  interp_work(1+num_pts*pm_work_size*(j-1):num_pts*pm_work_size*j)
               p(1:4*num_pts) => f1(1+4*num_pts*(j-1):4*num_pts*j)
               call interp_pm(x, num_pts, p, pm_work_size, work, &
                  'do_relax_composition', op_err)
               if (op_err /= 0) ierr = op_err
            end do
            deallocate(interp_work)
         end subroutine store_rpar

      end subroutine do_relax_composition


      subroutine before_evolve_relax_composition(s, id, lipar, ipar, lrpar, rpar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         ierr = 0
         call setup_before_relax(s)
      end subroutine before_evolve_relax_composition

      integer function relax_composition_adjust_model(s, id, lipar, ipar, lrpar, rpar)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer :: num_pts, num_steps_to_use, species, starting_model_number, ierr
         real(dp) :: lambda, avg_err

         real(dp), pointer :: x(:) ! =(num_pts)
         real(dp), pointer :: f1(:) ! =(4, num_pts, species)

         include 'formats'

         num_pts = ipar(1)
         num_steps_to_use = ipar(2)
         starting_model_number = ipar(3)
         species = s% species

         if (lrpar /= (1 + 4*species)*num_pts) then
            write(*,*) 'bad lrpar for relax_composition_check_model'
            relax_composition_adjust_model = terminate
            return
         end if

         ierr = 0
         if (s% job% timescale_for_relax_composition < 0d0) then
            lambda = min(1d0, max(0d0, (s% model_number - starting_model_number - 1) &
                 / max(1d0, dble(num_steps_to_use) - 1)))
         else
            lambda = min(1d0, s% dt/s% job% timescale_for_relax_composition)
         end if

         x(1:num_pts) => rpar(1:num_pts)
         f1(1:4*num_pts*species) => rpar(num_pts+1:lrpar)
         call adjust_xa(species, num_pts, avg_err, ierr)
         if (ierr /= 0) relax_composition_adjust_model = terminate

         write(*,1) 'avg remaining difference, lambda', avg_err, lambda

         relax_composition_adjust_model = keep_going

         contains


         subroutine adjust_xa(species, num_pts, avg_err, ierr)
            use interp_1d_lib, only: interp_values
            integer, intent(in) :: species, num_pts
            real(dp), intent(out) :: avg_err
            integer, intent(out) :: ierr
            integer :: j, k, nz, op_err
            real(dp), pointer :: vals(:,:), xq(:), f(:)
            ierr = 0
            nz = s% nz
            allocate(vals(species, nz), xq(nz), stat=ierr)
            if (ierr /= 0) return
            xq(1) = 0.5d0*s% dq(1) ! xq for cell center
            do k = 2, nz
               xq(k) = xq(k-1) + 0.5d0*(s% dq(k) + s% dq(k-1))
            end do
!$OMP PARALLEL DO PRIVATE(j,op_err,f) SCHEDULE(dynamic,2)
            do j=1, species ! interpolate target composition
               f(1:4*num_pts) => f1(1+(j-1)*4*num_pts:j*4*num_pts)
               call interp_values(x, num_pts, f, nz, xq, vals(j,:), op_err)
               if (op_err /= 0) ierr = op_err
               s% xa(j,1:nz) = (1d0-lambda)*s% xa(j,1:nz) + lambda*vals(j,1:nz)
            end do
!$OMP END PARALLEL DO
            avg_err = sum(abs(s% xa(1:species,1:nz)-vals(1:species,1:nz)))/(nz*species)
            deallocate(vals, xq)
         end subroutine adjust_xa

      end function relax_composition_adjust_model

      integer function relax_composition_check_model(s, id, lipar, ipar, lrpar, rpar)
         use do_one_utils, only:do_bare_bones_check_model
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)

         integer :: num_steps_to_use, starting_model_number, ierr

         include 'formats'

         relax_composition_check_model = do_bare_bones_check_model(id)
         if (relax_composition_check_model /= keep_going) return

         num_steps_to_use = ipar(2)
         starting_model_number = ipar(3)

         if (s% job% timescale_for_relax_composition < 0d0) then
            if (s% model_number - starting_model_number >= num_steps_to_use) then
               relax_composition_check_model = terminate
               s% termination_code = t_relax_finished_okay
               return
            end if
         else
            if (s% generations > 1 .and. &
               s% dt > s% job% timescale_for_relax_composition .and. &
               s% dt_old > s% job% timescale_for_relax_composition) then
               relax_composition_check_model = terminate
               s% termination_code = t_relax_finished_okay
               return
            end if
         end if


      end function relax_composition_check_model


      subroutine do_relax_to_xaccrete(id, num_steps_to_use, ierr)
         use adjust_xyz, only: get_xa_for_accretion

         integer, intent(in) :: id
         integer, intent(in) :: num_steps_to_use ! use this many steps to do conversion
         integer, intent(out) :: ierr

         integer, parameter :: lipar=2
         integer :: lrpar, max_model_number, species
         real(dp), pointer :: rpar(:)
         real(dp) :: starting_dt_next, mix_factor, dxdt_nuc_factor
         logical :: do_element_diffusion
         type (star_info), pointer :: s
         real(dp), pointer :: xa(:), f1(:), f(:,:,:)
         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:)
         ipar => ipar_ary

         include 'formats'

         ierr = 0

         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         species = s% species
         ipar(1) = s% model_number
         ipar(2) = num_steps_to_use
         if (num_steps_to_use <= 0) then
            ierr = -1
            write(*,2) 'invalid num_steps_to_use to relax_to_xaccrete', num_steps_to_use
            return
         end if

         lrpar = species
         allocate(rpar(lrpar), stat=ierr)
         if (ierr /= 0) return

         xa(1:species) => rpar(1:species)

         call get_xa_for_accretion(s, xa, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) &
               write(*, *) 'get_xa_for_accretion failed in relax_to_xaccrete'
            deallocate(rpar)
            return
         end if

         max_model_number = s% max_model_number
         s% max_model_number = num_steps_to_use + 1
         write(*,*) 'num_steps_to_use', num_steps_to_use

         dxdt_nuc_factor = s% dxdt_nuc_factor
         s% dxdt_nuc_factor = 0 ! turn off composition change by nuclear burning
         mix_factor = s% mix_factor
         s% mix_factor = 0 ! turn off mixing
         do_element_diffusion = s% do_element_diffusion
         do_element_diffusion = .false. ! turn off diffusion
         starting_dt_next = s% dt_next

         call do_internal_evolve( &
               id, before_evolve_relax_to_xaccrete, &
               relax_to_xaccrete_adjust_model, relax_to_xaccrete_check_model, &
               null_finish_model, .true., lipar, ipar, lrpar, rpar, ierr)

         s% max_model_number = max_model_number
         s% dt_next = starting_dt_next
         s% dxdt_nuc_factor = dxdt_nuc_factor
         s% mix_factor = mix_factor
         s% do_element_diffusion = do_element_diffusion

         call error_check('relax to xacrrete',ierr)

         deallocate(rpar)

      end subroutine do_relax_to_xaccrete


      subroutine before_evolve_relax_to_xaccrete(s, id, lipar, ipar, lrpar, rpar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         ierr = 0
         call setup_before_relax(s)
      end subroutine before_evolve_relax_to_xaccrete

      integer function relax_to_xaccrete_adjust_model(s, id, lipar, ipar, lrpar, rpar)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         relax_to_xaccrete_adjust_model = keep_going
      end function relax_to_xaccrete_adjust_model

      integer function relax_to_xaccrete_check_model(s, id, lipar, ipar, lrpar, rpar)
         use do_one_utils, only:do_bare_bones_check_model
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)

         integer :: num_steps_to_use, starting_model_number, species, k, j
         real(dp), pointer :: xa(:)
         real(dp) :: frac

         include 'formats'

         relax_to_xaccrete_check_model = do_bare_bones_check_model(id)
         if (relax_to_xaccrete_check_model /= keep_going) return

         starting_model_number = ipar(1)
         num_steps_to_use = ipar(2)
         species = s% species

         frac = dble(s% model_number - starting_model_number)/dble(num_steps_to_use)
         frac = frac*frac

         if (lrpar /= species) then
            write(*,*) 'bad lrpar for relax_to_xaccrete_check_model'
            relax_to_xaccrete_check_model = terminate
         end if

         xa(1:species) => rpar(1:species)

         do k=1,s% nz
            do j=1,species
               s% xa(j,k) = (1d0-frac)*s% xa(j,k) + frac*xa(j)
            end do
         end do

         if (mod(s% model_number, s% terminal_interval) == 0) &
            write(*,2) 'relax to xaccrete: fraction', s% model_number, frac

         if (s% model_number - starting_model_number >= num_steps_to_use) then
            relax_to_xaccrete_check_model = terminate
            s% termination_code = t_relax_finished_okay
            return
         end if


      end function relax_to_xaccrete_check_model

      subroutine do_relax_entropy(  &
            id, max_steps_to_use, num_pts, entropy, xq, ierr)
         use alloc, only: alloc_extras, dealloc_extras
         integer, intent(in) :: id
         integer, intent(in) :: max_steps_to_use ! use this many steps to do conversion
         integer, intent(in) :: num_pts
            ! length of entropy vector; need not equal nz for current model
         real(dp), intent(in) :: entropy(num_pts) ! desired entropy profile
         real(dp), intent(in) :: xq(num_pts)
         integer, intent(out) :: ierr

         integer, parameter ::  lipar=2
         integer :: lrpar, max_model_number
         real(dp), pointer :: rpar(:)
         real(dp) :: starting_dt_next, mix_factor, dxdt_nuc_factor, max_years_for_timestep
         logical :: do_element_diffusion, use_other_energy
         type (star_info), pointer :: s
         real(dp), pointer :: x(:), f1(:), f(:,:)
         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:)
         procedure (other_energy_interface), pointer :: &
            other_energy => null()

         ipar => ipar_ary

         include 'formats'

         ierr = 0

         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         max_model_number = s% max_model_number
         s% max_model_number = max_steps_to_use + s% model_number + 1
         write(*,*) 'relax_entropy: max_steps_to_use', max_steps_to_use

         ipar(1) = num_pts
         ipar(2) = s% max_model_number
         lrpar = 5*num_pts
         allocate(rpar(lrpar), stat=ierr)
         if (ierr /= 0) return

         x(1:num_pts) => rpar(1:num_pts)
         f1(1:4*num_pts) => rpar(num_pts+1:lrpar)
         f(1:4,1:num_pts) => f1(1:4*num_pts)

         call store_rpar(num_pts, ierr)
         if (ierr /= 0) return
         
         ! need to use global variables, as relax_entropy uses
         ! the other_energy routine to which it can't pass rpar
         relax_num_pts = num_pts
         relax_work_array => rpar

         dxdt_nuc_factor = s% dxdt_nuc_factor
         s% dxdt_nuc_factor = 0 ! turn off composition change by nuclear burning
         mix_factor = s% mix_factor
         s% mix_factor = 0 ! turn off mixing
         do_element_diffusion = s% do_element_diffusion
         s% do_element_diffusion = .false. ! turn off diffusion
         starting_dt_next = s% dt_next
         max_years_for_timestep = s% max_years_for_timestep
         s% max_years_for_timestep = s% job% max_dt_for_relax_entropy
         s% dt_next = min(s% dt_next, s% job% max_dt_for_relax_entropy * secyer)
         use_other_energy = s% use_other_energy
         s% use_other_energy = .true.
         other_energy => s% other_energy
         s% other_energy => entropy_relax_other_energy

         call do_internal_evolve( &
               id, before_evolve_relax_entropy, &
               relax_entropy_adjust_model, relax_entropy_check_model, &
               null_finish_model, .true., lipar, ipar, lrpar, rpar, ierr)

         s% max_model_number = max_model_number
         s% dt_next = starting_dt_next
         s% dxdt_nuc_factor = dxdt_nuc_factor
         s% mix_factor = mix_factor
         s% do_element_diffusion = do_element_diffusion
         s% max_years_for_timestep = max_years_for_timestep
         s% use_other_energy = use_other_energy
         s% other_energy => other_energy

         call error_check('relax entropy',ierr)

         deallocate(rpar)


         contains


         subroutine store_rpar(num_pts, ierr) ! get interpolation info
            use interp_1d_def, only: pm_work_size
            use interp_1d_lib, only: interp_pm
            integer, intent(in) :: num_pts
            integer, intent(out) :: ierr
            integer :: op_err
            real(dp), pointer :: interp_work(:), work(:), p(:)
            allocate(interp_work(num_pts*pm_work_size), stat=ierr)
            if (ierr /= 0) return
            x(:) = xq(:)
            op_err = 0
            f(1,1:num_pts) = entropy(1:num_pts)
            work(1:num_pts*pm_work_size) => &
               interp_work(1:num_pts*pm_work_size)
            p(1:4*num_pts) => f1(1:4*num_pts)
            call interp_pm(x, num_pts, p, pm_work_size, work, &
               'do_relax_entropy', op_err)
            if (op_err /= 0) ierr = op_err
            deallocate(interp_work)
         end subroutine store_rpar

      end subroutine do_relax_entropy

      subroutine before_evolve_relax_entropy(s, id, lipar, ipar, lrpar, rpar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         ierr = 0
         call setup_before_relax(s)
      end subroutine before_evolve_relax_entropy

      integer function relax_entropy_adjust_model(s, id, lipar, ipar, lrpar, rpar)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         relax_entropy_adjust_model = keep_going
      end function relax_entropy_adjust_model

      integer function relax_entropy_check_model(s, id, lipar, ipar, lrpar, rpar)
         use do_one_utils, only: do_bare_bones_check_model
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)

         integer :: num_pts, ierr, max_model_number
         real(dp) :: lambda, avg_err
         real(dp), pointer :: x(:) ! =(num_pts)
         real(dp), pointer :: f1(:) ! =(4, num_pts)

         include 'formats'

         relax_entropy_check_model = do_bare_bones_check_model(id)
         if (relax_entropy_check_model /= keep_going) return

         if (lipar /= 2) then
            write(*,*) 'bad lipar for relax_entropy_check_model'
            relax_entropy_check_model = terminate
            return
         end if

         num_pts = ipar(1)
         max_model_number = ipar(2)

         if (lrpar /= 5*num_pts) then
            write(*,*) 'bad lrpar for relax_entropy_check_model'
            relax_entropy_check_model = terminate
         end if

         ierr = 0
         x(1:num_pts) => rpar(1:num_pts)
         f1(1:4*num_pts) => rpar(num_pts+1:lrpar)
         call adjust_entropy(num_pts, avg_err, ierr)
         if (ierr /= 0) relax_entropy_check_model = terminate

         if (mod(s% model_number, s% terminal_interval) == 0) &
            write(*,*) 'relax_entropy avg rel err, dt, model', avg_err, s% dt/secyer, s% model_number

         if (s% star_age >= s% job% timescale_for_relax_entropy*s% job% num_timescales_for_relax_entropy) then
            relax_entropy_check_model = terminate
            s% termination_code = t_relax_finished_okay
            return
         end if

         if (s% model_number >= max_model_number) then
            write(*,*) "Terminated relax because of max_model_number instead of relax_time"
            relax_entropy_check_model = terminate
            s% termination_code = t_relax_finished_okay
            return
         end if


         contains


         subroutine adjust_entropy(num_pts, avg_err, ierr)
            use interp_1d_lib, only: interp_values
            integer, intent(in) :: num_pts
            real(dp), intent(out) :: avg_err
            integer, intent(out) :: ierr
            integer :: k, nz, op_err
            real(dp) :: dentropy_sum
            real(dp), pointer :: vals(:), xq(:), f(:), entropy(:)
            ierr = 0
            nz = s% nz
            allocate(vals(nz), xq(nz), stat=ierr)
            if (ierr /= 0) return
            xq(1) = s% dq(1)/2 ! xq for cell center
            do k = 2, nz
               xq(k) = xq(k-1) + (s% dq(k) + s% dq(k-1))/2
            end do
            dentropy_sum = 0
            f(1:4*num_pts) => f1(1:4*num_pts)
            call interp_values(x, num_pts, f, nz, xq, vals(:), op_err)
            if (op_err /= 0) ierr = op_err
            allocate(entropy(1:nz))
            do k=1,nz
               entropy(k) = exp(s% lnS(k))
            end do
            
            dentropy_sum = sum(abs((entropy(1:nz)-vals(1:nz))/vals(1:nz)))
            avg_err = dentropy_sum/nz
            deallocate(vals, xq, entropy)
         end subroutine adjust_entropy

      end function relax_entropy_check_model

      subroutine entropy_relax_other_energy(id, ierr)
         use interp_1d_lib, only: interp_values
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, nz, num_pts, op_err
         real(dp), pointer :: vals(:), xq(:), x(:), f(:)
         ierr = 0
         call star_ptr(id, s, ierr)

         nz = s% nz
         num_pts = relax_num_pts
         allocate(vals(nz), xq(nz), stat=ierr)
         f(1:4*num_pts) => relax_work_array(num_pts+1:5*num_pts)
         x(1:num_pts) => relax_work_array(1:num_pts)
         xq(1) = s% dq(1)/2 ! xq for cell center
         do k = 2, nz
            xq(k) = xq(k-1) + (s% dq(k) + s% dq(k-1))/2
         end do
         call interp_values(x, num_pts, f, nz, xq, vals(:), op_err)
         if (op_err /= 0) ierr = op_err

         if (ierr /= 0) return
         s% extra_heat(:) = 0d0
         do k = 1, s% nz
            s% extra_heat(k) =  ( 1d0 - exp(s%lnS(k))/vals(k) ) * exp(s%lnE(k))
            s% extra_heat(k) = s% extra_heat(k) / (s% job% timescale_for_relax_entropy * secyer)
         end do
      end subroutine entropy_relax_other_energy

      subroutine do_relax_angular_momentum(  &
            id, max_steps_to_use, num_pts, angular_momentum, xq, ierr)
         use alloc, only: alloc_extras, dealloc_extras
         integer, intent(in) :: id
         integer, intent(in) :: max_steps_to_use ! use this many steps to do conversion
         integer, intent(in) :: num_pts
            ! length of angular momentum vector; need not equal nz for current model
         real(dp), intent(in) :: angular_momentum(num_pts) ! desired angular momentum profile
         real(dp), intent(in) :: xq(num_pts)
         integer, intent(out) :: ierr

         integer, parameter ::  lipar=2
         integer :: lrpar, max_model_number
         real(dp), pointer :: rpar(:)
         real(dp) :: starting_dt_next, mix_factor, dxdt_nuc_factor, max_timestep, &
            am_D_mix_factor
         logical :: do_element_diffusion, use_other_torque
         type (star_info), pointer :: s
         real(dp), pointer :: x(:), f1(:), f(:,:)
         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:)
         procedure (other_torque_interface), pointer :: &
            other_torque => null()

         ipar => ipar_ary

         include 'formats'

         ierr = 0

         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         max_model_number = s% max_model_number
         s% max_model_number = max_steps_to_use + s% model_number + 1
         write(*,*) 'relax_angular_momentum: max_steps_to_use', max_steps_to_use

         ipar(1) = num_pts
         ipar(2) = s% max_model_number
         lrpar = 5*num_pts
         allocate(rpar(lrpar), stat=ierr)
         if (ierr /= 0) return

         x(1:num_pts) => rpar(1:num_pts)
         f1(1:4*num_pts) => rpar(num_pts+1:lrpar)
         f(1:4,1:num_pts) => f1(1:4*num_pts)

         call store_rpar(num_pts, ierr)
         if (ierr /= 0) return
         
         ! need to use global variables, as relax_angular_momentum uses
         ! the other_torque routine to which it can't pass rpar
         relax_num_pts = num_pts
         relax_work_array => rpar

         dxdt_nuc_factor = s% dxdt_nuc_factor
         s% dxdt_nuc_factor = 0 ! turn off composition change by nuclear burning
         mix_factor = s% mix_factor
         s% mix_factor = 0 ! turn off mixing
         am_D_mix_factor = s% am_D_mix_factor
         s% am_D_mix_factor = 0d0
         do_element_diffusion = s% do_element_diffusion
         s% do_element_diffusion = .false. ! turn off diffusion
         starting_dt_next = s% dt_next
         max_timestep = s% max_timestep
         s% max_timestep = s% job% max_dt_for_relax_angular_momentum * secyer
         s% dt_next = min(s% dt_next, s% job% max_dt_for_relax_angular_momentum * secyer)
         use_other_torque = s% use_other_torque
         s% use_other_torque = .true.
         other_torque => s% other_torque
         s% other_torque => angular_momentum_relax_other_torque

         call do_internal_evolve( &
               id, before_evolve_relax_angular_momentum, &
               relax_angular_momentum_adjust_model, relax_angular_momentum_check_model, &
               null_finish_model, .true., lipar, ipar, lrpar, rpar, ierr)

         s% max_model_number = max_model_number
         s% dt_next = starting_dt_next
         s% dxdt_nuc_factor = dxdt_nuc_factor
         s% mix_factor = mix_factor
         s% am_D_mix_factor = am_D_mix_factor
         s% do_element_diffusion = do_element_diffusion
         s% max_timestep = max_timestep
         s% use_other_torque = use_other_torque
         s% other_torque => other_torque

         call error_check('relax angular momentum',ierr)

         deallocate(rpar)


         contains


         subroutine store_rpar(num_pts, ierr) ! get interpolation info
            use interp_1d_def, only: pm_work_size
            use interp_1d_lib, only: interp_pm
            integer, intent(in) :: num_pts
            integer, intent(out) :: ierr
            integer :: op_err
            real(dp), pointer :: interp_work(:), work(:), p(:)
            allocate(interp_work(num_pts*pm_work_size), stat=ierr)
            if (ierr /= 0) return
            x(:) = xq(:)
            op_err = 0
            f(1,1:num_pts) = angular_momentum(1:num_pts)
            work(1:num_pts*pm_work_size) => &
               interp_work(1:num_pts*pm_work_size)
            p(1:4*num_pts) => f1(1:4*num_pts)
            call interp_pm(x, num_pts, p, pm_work_size, work, &
               'do_relax_entropy', op_err)
            if (op_err /= 0) ierr = op_err
            deallocate(interp_work)
         end subroutine store_rpar

      end subroutine do_relax_angular_momentum

      subroutine before_evolve_relax_angular_momentum(s, id, lipar, ipar, lrpar, rpar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         ierr = 0
         call setup_before_relax(s)
      end subroutine before_evolve_relax_angular_momentum

      integer function relax_angular_momentum_adjust_model(s, id, lipar, ipar, lrpar, rpar)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         relax_angular_momentum_adjust_model = keep_going
      end function relax_angular_momentum_adjust_model

      integer function relax_angular_momentum_check_model(s, id, lipar, ipar, lrpar, rpar)
         use do_one_utils, only: do_bare_bones_check_model
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)

         integer :: num_pts, ierr, max_model_number
         real(dp) :: lambda, avg_err
         real(dp), pointer :: x(:) ! =(num_pts)
         real(dp), pointer :: f1(:) ! =(4, num_pts)

         include 'formats'

         relax_angular_momentum_check_model = do_bare_bones_check_model(id)
         if (relax_angular_momentum_check_model /= keep_going) return

         if (lipar /= 2) then
            write(*,*) 'bad lipar for relax_angular_momentum_check_model'
            relax_angular_momentum_check_model = terminate
            return
         end if

         num_pts = ipar(1)
         max_model_number = ipar(2)

         if (lrpar /= 5*num_pts) then
            write(*,*) 'bad lrpar for relax_angular_momentum_check_model'
            relax_angular_momentum_check_model = terminate
         end if

         ierr = 0
         x(1:num_pts) => rpar(1:num_pts)
         f1(1:4*num_pts) => rpar(num_pts+1:lrpar)
         call adjust_angular_momentum(num_pts, avg_err, ierr)
         if (ierr /= 0) relax_angular_momentum_check_model = terminate

         if (mod(s% model_number, s% terminal_interval) == 0) &
            write(*,*) 'relax_angular_momentum avg rel err, dt, model', avg_err, s% dt/secyer, s% model_number

         if (s% star_age >= &
            s% job% timescale_for_relax_angular_momentum*&
               s% job% num_timescales_for_relax_angular_momentum) then
            relax_angular_momentum_check_model = terminate
            s% termination_code = t_relax_finished_okay
            return
         end if

         if (s% model_number >= max_model_number) then
            write(*,*) "Terminated relax because of max_model_number instead of relax_time"
            relax_angular_momentum_check_model = terminate
            s% termination_code = t_relax_finished_okay
            return
         end if


         contains


         subroutine adjust_angular_momentum(num_pts, avg_err, ierr)
            use interp_1d_lib, only: interp_values
            integer, intent(in) :: num_pts
            real(dp), intent(out) :: avg_err
            integer, intent(out) :: ierr
            integer :: k, nz, op_err
            real(dp) :: dangular_momentum_sum
            real(dp), pointer :: vals(:), xq(:), f(:)
            ierr = 0
            nz = s% nz
            allocate(vals(nz), xq(nz), stat=ierr)
            if (ierr /= 0) return
            xq(1) = s% dq(1)/2 ! xq for cell center
            do k = 2, nz
               xq(k) = xq(k-1) + (s% dq(k) + s% dq(k-1))/2
            end do
            dangular_momentum_sum = 0
            f(1:4*num_pts) => f1(1:4*num_pts)
            call interp_values(x, num_pts, f, nz, xq, vals(:), op_err)
            if (op_err /= 0) ierr = op_err
            dangular_momentum_sum = sum(abs((s% j_rot(1:nz)-vals(1:nz))/vals(1:nz)))
            avg_err = dangular_momentum_sum/nz
            deallocate(vals, xq)
         end subroutine adjust_angular_momentum

      end function relax_angular_momentum_check_model

      subroutine angular_momentum_relax_other_torque(id, ierr)
         use interp_1d_lib, only: interp_values
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, k_inner, nz, num_pts, op_err
         real(dp), pointer :: vals(:), xq(:), x(:), f(:)
         real(dp) :: omega_target, j_rot_target
         ierr = 0
         call star_ptr(id, s, ierr)

         nz = s% nz
         num_pts = relax_num_pts
         allocate(vals(nz), xq(nz), stat=ierr)
         f(1:4*num_pts) => relax_work_array(num_pts+1:5*num_pts)
         x(1:num_pts) => relax_work_array(1:num_pts)
         xq(1) = s% dq(1)/2 ! xq for cell center
         do k = 2, nz
            xq(k) = xq(k-1) + (s% dq(k) + s% dq(k-1))/2
         end do
         call interp_values(x, num_pts, f, nz, xq, vals(:), op_err)
         if (op_err /= 0) ierr = op_err

         if (ierr /= 0) return
         s% extra_jdot(:) = 0d0
         do k = 1, s% nz
            if (s% j_rot_flag) then
               s% extra_jdot(k) =  &
                  (1d0 - exp(-s% dt/(s% job% timescale_for_relax_angular_momentum*secyer))) * &
                  (vals(k) - s% j_rot_start(k))/s% dt
            else
               s% extra_jdot(k) =  &
                  (1d0 - exp(-s% dt/(s% job% timescale_for_relax_angular_momentum*secyer))) * &
                  (vals(k) - s% j_rot(k))/s% dt
            end if
         end do
         ! for cells near center use a constant omega, prevents extrapolating j_rot and causing artificially large core rotation
         k_inner = -1
         do k = 2, s% nz
            if (xq(k) > x(num_pts)) then
               k_inner = k-1
            end if
         end do
         if (s% job% relax_angular_momentum_constant_omega_center) then
            if (k_inner > 0) then
               omega_target = vals(k_inner)/s% i_rot(k_inner)
               do k=k_inner+1, s% nz
                  j_rot_target = omega_target*s% i_rot(k)
                  if (s% j_rot_flag) then
                     s% extra_jdot(k) =  &
                        (1d0 - exp(-s% dt/(s% job% timescale_for_relax_angular_momentum*secyer))) * &
                        (j_rot_target - s% j_rot_start(k))/s% dt
                  else
                     s% extra_jdot(k) =  &
                        (1d0 - exp(-s% dt/(s% job% timescale_for_relax_angular_momentum*secyer))) * &
                        (j_rot_target - s% j_rot(k))/s% dt
                  end if
               end do
            end if
         end if
         deallocate(vals, xq)
      end subroutine angular_momentum_relax_other_torque

      subroutine do_relax_uniform_omega( &
            id, kind_of_relax, target_value, num_steps_to_relax_rotation,&
            relax_omega_max_yrs_dt, ierr)
         integer, intent(in) :: id, kind_of_relax
         real(dp), intent(in) :: target_value,relax_omega_max_yrs_dt
         integer, intent(in) :: num_steps_to_relax_rotation
         integer, intent(out) :: ierr

         integer, parameter :: lipar=3, lrpar=2
         integer :: max_model_number, k
         real(dp) ::  max_years_for_timestep, mix_factor, dxdt_nuc_factor, am_D_mix_factor
         type (star_info), pointer :: s
         real(dp), target :: rpar_ary(lrpar)
         real(dp), pointer :: rpar(:)
         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:)
         logical :: okay

         include 'formats'

         ierr = 0

         rpar => rpar_ary
         ipar => ipar_ary

         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         if (.not. s% rotation_flag) return

         rpar(1) = target_value
         okay = .true.
         do k=1,s% nz
            if (s% omega(k) /= target_value) then
               okay = .false.
               exit
            end if
         end do
         if (okay) return

         rpar(2) = s% omega(1)
         ipar(1) = s% model_number
         ipar(2) = num_steps_to_relax_rotation
         ipar(3) = kind_of_relax
         if (num_steps_to_relax_rotation <= 0) then
            ierr = -1
            write(*,2) 'invalid num_steps_to_relax_rotation', num_steps_to_relax_rotation
            return
         end if

         dxdt_nuc_factor = s% dxdt_nuc_factor
         s% dxdt_nuc_factor = 0d0 ! turn off composition change by nuclear burning
         mix_factor = s% mix_factor
         am_D_mix_factor = s% am_D_mix_factor
         s% mix_factor = 0d0
         s% am_D_mix_factor = 0d0
         max_model_number = s% max_model_number
         s% max_model_number = num_steps_to_relax_rotation + 1
         max_years_for_timestep = s% max_years_for_timestep
         s% max_years_for_timestep = relax_omega_max_yrs_dt
         write(*,*) 'num_steps_to_relax_rotation', num_steps_to_relax_rotation

         call do_internal_evolve( &
               id, before_evolve_relax_omega, &
               relax_omega_adjust_model, relax_omega_check_model, &
               null_finish_model, .true., lipar, ipar, lrpar, rpar, ierr)

         s% dxdt_nuc_factor = dxdt_nuc_factor
         s% mix_factor = mix_factor
         s% am_D_mix_factor = am_D_mix_factor
         s% max_model_number = max_model_number
         s% max_years_for_timestep = max_years_for_timestep
         call error_check('relax uniform omega',ierr)
         s% D_omega(1:s% nz) = 0
         s% am_nu_rot(1:s% nz) = 0

      end subroutine do_relax_uniform_omega


      subroutine before_evolve_relax_omega(s, id, lipar, ipar, lrpar, rpar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         ierr = 0
         call setup_before_relax(s)
      end subroutine before_evolve_relax_omega

      integer function relax_omega_adjust_model(s, id, lipar, ipar, lrpar, rpar)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         relax_omega_adjust_model = keep_going
      end function relax_omega_adjust_model

      integer function relax_omega_check_model(s, id, lipar, ipar, lrpar, rpar)
         use do_one_utils, only: do_bare_bones_check_model
         use hydro_rotation, only: set_uniform_omega, set_i_rot
         use star_utils, only: set_surf_avg_rotation_info
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)

         integer :: num_steps_to_use, starting_model_number, kind_of_relax, ierr
         real(dp) :: frac, target_value, starting_omega, new_omega, this_step_omega

         include 'formats'

         relax_omega_check_model = do_bare_bones_check_model(id)
         if (relax_omega_check_model /= keep_going) return

         starting_model_number = ipar(1)
         num_steps_to_use = ipar(2)
         kind_of_relax = ipar(3)
         target_value = rpar(1)
         starting_omega = rpar(2)

         frac = max(0d0, min(1d0, &
            dble(s% model_number - starting_model_number)/dble(num_steps_to_use)))

         if (kind_of_relax == relax_to_new_omega) then
            new_omega = target_value
         else if (kind_of_relax == relax_to_new_omega_div_omega_crit) then
            call set_i_rot(s, .false.)
            call set_surf_avg_rotation_info(s)
            new_omega = target_value*s% omega_crit_avg_surf
         else if (kind_of_relax == relax_to_new_surface_rotation_v) then
            new_omega = target_value*1d5/(s% photosphere_r*Rsun)
         else
            write(*,2) 'bad value for kind_of_relax', kind_of_relax
            stop 'relax_omega_check_model'
         end if

         this_step_omega = frac*new_omega + (1 - frac)*starting_omega

         if (s% model_number > starting_model_number + num_steps_to_use + 50 .or. &
               abs(s% omega(1) - new_omega) < 1d-4*max(1d-10,abs(new_omega))) then
            write(*,2) 'final step: wanted-current, current, wanted', &
               s% model_number, new_omega-s% omega(1), s% omega(1), new_omega
            relax_omega_check_model = terminate
            s% termination_code = t_relax_finished_okay
         else if (mod(s% model_number, s% terminal_interval) == 0) then
            write(*,2) 'relax to omega: wanted-current, current, wanted', &
               s% model_number, new_omega-s% omega(1), s% omega(1), new_omega
         end if

         ierr = 0
         call set_uniform_omega(id, this_step_omega, ierr)
         if (ierr /= 0) then
            write(*,*) 'set_uniform_omega failed'
            relax_omega_check_model = terminate
            return
         end if

      end function relax_omega_check_model


      subroutine do_relax_tau_factor(id, new_tau_factor, dlogtau_factor, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: new_tau_factor, dlogtau_factor
         integer, intent(out) :: ierr
         integer, parameter ::  lipar=1, lrpar=2
         integer :: max_model_number
         real(dp) :: tau_factor
         type (star_info), pointer :: s
         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:)
         real(dp), target :: rpar_ary(lrpar)
         real(dp), pointer :: rpar(:)
         rpar => rpar_ary
         ipar => ipar_ary
         include 'formats'
         ierr = 0
         if (new_tau_factor <= 0) then
            ierr = -1
            write(*,*) 'invalid new_tau_factor', new_tau_factor
            return
         end if
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         tau_factor = s% tau_factor
         if (abs(new_tau_factor - tau_factor) <= 1d-6) then
            s% tau_factor = new_tau_factor
            return
         end if
         write(*,*)
         write(*,1) 'current tau_factor', tau_factor
         write(*,1) 'relax to new tau_factor', new_tau_factor
         write(*,*)
         write(*,1) 'dlogtau_factor', dlogtau_factor
         write(*,*)
         rpar(1) = new_tau_factor
         rpar(2) = dlogtau_factor
         max_model_number = s% max_model_number
         s% max_model_number = -1111
         call do_internal_evolve( &
               id, before_evolve_relax_tau_factor, &
               relax_tau_factor_adjust_model, relax_tau_factor_check_model, &
               null_finish_model, .true., lipar, ipar, lrpar, rpar, ierr)
         s% max_model_number = max_model_number
         if (ierr == 0) s% force_tau_factor = s% tau_factor
         call error_check('relax tau factor',ierr)
      end subroutine do_relax_tau_factor


      subroutine before_evolve_relax_tau_factor(s, id, lipar, ipar, lrpar, rpar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         ierr = 0
         call turn_off_winds(s)
         s% max_model_number = -111
      end subroutine before_evolve_relax_tau_factor

      integer function relax_tau_factor_adjust_model(s, id, lipar, ipar, lrpar, rpar)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         relax_tau_factor_adjust_model = keep_going
      end function relax_tau_factor_adjust_model

      integer function relax_tau_factor_check_model(s, id, lipar, ipar, lrpar, rpar)
         use do_one_utils, only:do_bare_bones_check_model
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         real(dp) :: new_tau_factor, dlogtau_factor, current_tau_factor, next
         logical, parameter :: dbg = .false.

         include 'formats'

         relax_tau_factor_check_model = do_bare_bones_check_model(id)
         if (relax_tau_factor_check_model /= keep_going) return

         new_tau_factor = rpar(1)
         dlogtau_factor = rpar(2)
         current_tau_factor = s% tau_factor

         if (mod(s% model_number, s% terminal_interval) == 0) &
            write(*,1) 'tau_factor target current', new_tau_factor, current_tau_factor

         if (abs(current_tau_factor-new_tau_factor) < 1d-15) then
            s% tau_factor = new_tau_factor
            s% termination_code = t_relax_finished_okay
            relax_tau_factor_check_model = terminate
            return
         end if

         if (new_tau_factor < current_tau_factor) then
            next = exp10(safe_log10(current_tau_factor) - dlogtau_factor)
            if (next < new_tau_factor) next = new_tau_factor
         else
            next = exp10(safe_log10(current_tau_factor) + dlogtau_factor)
            if (next > new_tau_factor) next = new_tau_factor
         end if

         if (dbg) write(*,1) 'next', next, log10(next)

         s% tau_factor = next
         s% max_timestep = secyer*s% time_step

      end function relax_tau_factor_check_model


      subroutine do_relax_opacity_factor(id, new_opacity_factor, dopacity_factor, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: new_opacity_factor, dopacity_factor
         integer, intent(out) :: ierr
         integer, parameter ::  lipar=1, lrpar=2
         integer :: max_model_number
         real(dp) :: opacity_factor
         type (star_info), pointer :: s
         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:)
         real(dp), target :: rpar_ary(lrpar)
         real(dp), pointer :: rpar(:)
         rpar => rpar_ary
         ipar => ipar_ary
         include 'formats'
         ierr = 0
         if (new_opacity_factor <= 0) then
            ierr = -1
            write(*,*) 'invalid new_opacity_factor', new_opacity_factor
            return
         end if
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         opacity_factor = s% opacity_factor
         if (abs(new_opacity_factor - opacity_factor) <= 1d-6) then
            s% opacity_factor = new_opacity_factor
            return
         end if
         write(*,*)
         write(*,1) 'current opacity_factor', opacity_factor
         write(*,1) 'relax to new opacity_factor', new_opacity_factor
         write(*,*)
         write(*,1) 'dopacity_factor', dopacity_factor
         write(*,*)
         rpar(1) = new_opacity_factor
         rpar(2) = dopacity_factor
         max_model_number = s% max_model_number
         s% max_model_number = -1111
         call do_internal_evolve( &
               id, before_evolve_relax_opacity_factor, &
               relax_opacity_factor_adjust_model, relax_opacity_factor_check_model, &
               null_finish_model, .true., lipar, ipar, lrpar, rpar, ierr)
         s% max_model_number = max_model_number
         if (ierr == 0) s% force_opacity_factor = s% opacity_factor
         call error_check('relax opacity factor',ierr)
      end subroutine do_relax_opacity_factor


      subroutine before_evolve_relax_opacity_factor(s, id, lipar, ipar, lrpar, rpar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         ierr = 0
         call turn_off_winds(s)
         s% max_model_number = -111
      end subroutine before_evolve_relax_opacity_factor

      integer function relax_opacity_factor_adjust_model(s, id, lipar, ipar, lrpar, rpar)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         relax_opacity_factor_adjust_model = keep_going
      end function relax_opacity_factor_adjust_model

      integer function relax_opacity_factor_check_model(s, id, lipar, ipar, lrpar, rpar)
         use do_one_utils, only:do_bare_bones_check_model
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         real(dp) :: new_opacity_factor, dopacity_factor, current_opacity_factor, next
         logical, parameter :: dbg = .false.

         include 'formats'

         relax_opacity_factor_check_model = do_bare_bones_check_model(id)
         if (relax_opacity_factor_check_model /= keep_going) return

         new_opacity_factor = rpar(1)
         dopacity_factor = rpar(2)
         current_opacity_factor = s% opacity_factor

         if (dbg) then
            write(*,1) 'new_opacity_factor', new_opacity_factor
            write(*,1) 'current_opacity_factor', current_opacity_factor
         end if

         if (abs(current_opacity_factor-new_opacity_factor) < 1d-15) then
            s% opacity_factor = new_opacity_factor
            s% termination_code = t_relax_finished_okay
            relax_opacity_factor_check_model = terminate
            return
         end if

         if (new_opacity_factor < current_opacity_factor) then
            next = current_opacity_factor - dopacity_factor
            if (next < new_opacity_factor) next = new_opacity_factor
         else
            next = current_opacity_factor + dopacity_factor
            if (next > new_opacity_factor) next = new_opacity_factor
         end if

         if (mod(s% model_number, s% terminal_interval) == 0) &
            write(*,1) 'next opacity_factor', next ! provide terminal feedback to show working.

         s% opacity_factor = next
         s% max_timestep = secyer*s% time_step

      end function relax_opacity_factor_check_model


      subroutine do_relax_Tsurf_factor(id, new_Tsurf_factor, dlogTsurf_factor, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: new_Tsurf_factor, dlogTsurf_factor
         integer, intent(out) :: ierr
         integer, parameter ::  lipar=1, lrpar=2
         integer :: max_model_number
         real(dp) :: Tsurf_factor
         type (star_info), pointer :: s
         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:)
         real(dp), target :: rpar_ary(lrpar)
         real(dp), pointer :: rpar(:)
         rpar => rpar_ary
         ipar => ipar_ary
         include 'formats'
         ierr = 0
         if (new_Tsurf_factor <= 0) then
            ierr = -1
            write(*,*) 'invalid new_Tsurf_factor', new_Tsurf_factor
            return
         end if
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         Tsurf_factor = s% Tsurf_factor
         if (abs(new_Tsurf_factor - Tsurf_factor) <= 1d-6) then
            s% Tsurf_factor = new_Tsurf_factor
            return
         end if
         write(*,*)
         write(*,1) 'current Tsurf_factor', Tsurf_factor
         write(*,1) 'relax to new Tsurf_factor', new_Tsurf_factor
         write(*,*)
         write(*,1) 'dlogTsurf_factor', dlogTsurf_factor
         write(*,*)
         rpar(1) = new_Tsurf_factor
         rpar(2) = dlogTsurf_factor
         max_model_number = s% max_model_number
         s% max_model_number = -1111
         call do_internal_evolve( &
               id, before_evolve_relax_Tsurf_factor, &
               relax_Tsurf_factor_adjust_model, relax_Tsurf_factor_check_model, &
               null_finish_model, .true., lipar, ipar, lrpar, rpar, ierr)
         s% max_model_number = max_model_number
         if (ierr /= 0) then
            write(*,*)
            write(*,1) 'ERROR: failed doing relax Tsurf_factor', new_Tsurf_factor
            write(*,*)
            stop 'do_relax_Tsurf_factor'
         end if
         
         if (new_Tsurf_factor == 1d0) then
            s% force_Tsurf_factor = 0d0
         else
            s% force_Tsurf_factor = s% Tsurf_factor
         end if
         
         call error_check('relax tsurf factor',ierr)

      end subroutine do_relax_Tsurf_factor


      subroutine before_evolve_relax_Tsurf_factor(s, id, lipar, ipar, lrpar, rpar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         ierr = 0
         call turn_off_winds(s)
         s% max_model_number = -111
      end subroutine before_evolve_relax_Tsurf_factor

      integer function relax_Tsurf_factor_adjust_model(s, id, lipar, ipar, lrpar, rpar)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         relax_Tsurf_factor_adjust_model = keep_going
      end function relax_Tsurf_factor_adjust_model

      integer function relax_Tsurf_factor_check_model(s, id, lipar, ipar, lrpar, rpar)
         use do_one_utils, only:do_bare_bones_check_model
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         real(dp) :: new_Tsurf_factor, dlogTsurf_factor, current_Tsurf_factor, next
         logical, parameter :: dbg = .false.

         include 'formats'

         relax_Tsurf_factor_check_model = do_bare_bones_check_model(id)
         if (relax_Tsurf_factor_check_model /= keep_going) return

         new_Tsurf_factor = rpar(1)
         dlogTsurf_factor = rpar(2)
         current_Tsurf_factor = s% Tsurf_factor

         if (mod(s% model_number, s% terminal_interval) == 0) &
            write(*,1) 'Tsurf_factor target current', new_Tsurf_factor, current_Tsurf_factor

         if (abs(current_Tsurf_factor-new_Tsurf_factor) < 1d-15) then
            s% Tsurf_factor = new_Tsurf_factor
            s% termination_code = t_relax_finished_okay
            relax_Tsurf_factor_check_model = terminate
            return
         end if

         if (new_Tsurf_factor < current_Tsurf_factor) then
            next = exp10(safe_log10(current_Tsurf_factor) - dlogTsurf_factor)
            if (next < new_Tsurf_factor) next = new_Tsurf_factor
         else
            next = exp10(safe_log10(current_Tsurf_factor) + dlogTsurf_factor)
            if (next > new_Tsurf_factor) next = new_Tsurf_factor
         end if

         if (dbg) write(*,1) 'next Tsurf_factor', next, log10(next)

         s% Tsurf_factor = next
         s% max_timestep = secyer*s% time_step

      end function relax_Tsurf_factor_check_model


      subroutine do_relax_irradiation(id, &
            min_steps, new_irrad_flux, new_irrad_col_depth, relax_irradiation_max_yrs_dt, ierr)
         integer, intent(in) :: id, min_steps
         real(dp), intent(in) :: &
            new_irrad_flux, new_irrad_col_depth, relax_irradiation_max_yrs_dt
         integer, intent(out) :: ierr

         integer, parameter ::  lipar=2, lrpar=4
         integer :: max_model_number, i
         real(dp) :: max_years_for_timestep
         type (star_info), pointer :: s
         logical :: all_same
         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:)
         real(dp), target :: rpar_ary(lrpar)
         real(dp), pointer :: rpar(:)
         rpar => rpar_ary
         ipar => ipar_ary
         include 'formats'
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ipar(1) = s% model_number
         ipar(2) = min_steps

         rpar(1) = new_irrad_flux
         rpar(2) = new_irrad_col_depth
         rpar(3) = s% irradiation_flux
         rpar(4) = s% column_depth_for_irradiation

         all_same = .true.
         do i = 1, 2
            if (abs(rpar(i)-rpar(i+2)) > 1d-10) then
               all_same = .false.; exit
            end if
         end do
         if (all_same) return

         write(*,*)
         write(*,2) 'relax to new irradiation -- min steps', min_steps
         write(*,*)

         max_model_number = s% max_model_number
         s% max_model_number = -1111
         max_years_for_timestep = s% max_years_for_timestep
         s% max_years_for_timestep = relax_irradiation_max_yrs_dt
         call do_internal_evolve( &
               id, before_evolve_relax_irradiation, &
               relax_irradiation_adjust_model, relax_irradiation_check_model, &
               null_finish_model, .true., lipar, ipar, lrpar, rpar, ierr)
         s% max_model_number = max_model_number
         s% max_years_for_timestep = max_years_for_timestep
         call error_check('relax irradiation',ierr)
      end subroutine do_relax_irradiation


      subroutine before_evolve_relax_irradiation(s, id, lipar, ipar, lrpar, rpar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         ierr = 0
         call turn_off_winds(s)
         s% max_model_number = -111
      end subroutine before_evolve_relax_irradiation

      integer function relax_irradiation_adjust_model(s, id, lipar, ipar, lrpar, rpar)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         relax_irradiation_adjust_model = keep_going
      end function relax_irradiation_adjust_model

      integer function relax_irradiation_check_model(s, id, lipar, ipar, lrpar, rpar)
         use do_one_utils, only:do_bare_bones_check_model
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)

         integer :: adjust_model, max_num_steps, num_steps
         real(dp) :: old_irrad_flux, old_irrad_col_depth
         real(dp) :: new_irrad_flux, new_irrad_col_depth, frac
         logical, parameter :: dbg = .false.

         include 'formats'

         relax_irradiation_check_model = do_bare_bones_check_model(id)
         if (relax_irradiation_check_model /= keep_going) return

         adjust_model = ipar(1)
         max_num_steps = ipar(2)

         new_irrad_flux = rpar(1)
         new_irrad_col_depth = rpar(2)
         old_irrad_flux = rpar(3)
         old_irrad_col_depth = rpar(4)
         num_steps = s% model_number - adjust_model
         frac = dble(num_steps)/dble(max_num_steps)

         if (s% dt < s% max_years_for_timestep*secyer) then
            ipar(1) = adjust_model + 1
            write(*,'(a60,2i6,3x,99e12.3)') 'relax irradiation, model, step, frac, flux, wait for dt', &
               s% model_number, num_steps, frac, s% irradiation_flux
            return
         end if

         if (frac >= 1d0) then
            s% irradiation_flux = new_irrad_flux
            s% column_depth_for_irradiation = new_irrad_col_depth
            relax_irradiation_check_model = terminate
            s% termination_code = t_relax_finished_okay
            write(*,'(a60,2i6,3x,99e12.3)') &
               'DONE: relax irradiation, model, step, fraction done, flux', &
               s% model_number, num_steps, frac, s% irradiation_flux
            return
         end if

         s% irradiation_flux = &
            frac*new_irrad_flux + (1-frac)*old_irrad_flux
         s% column_depth_for_irradiation = &
            frac*new_irrad_col_depth + (1-frac)*old_irrad_col_depth

         if (mod(s% model_number, s% terminal_interval) == 0) &
            write(*,'(a60,2i6,3x,99e12.3)') 'relax irradiation, model, step, fraction done, flux', &
               s% model_number, num_steps, frac, s% irradiation_flux

      end function relax_irradiation_check_model


      subroutine do_relax_mass_change( &
            id, min_steps, initial_mass_change, final_mass_change, relax_mass_change_max_yrs_dt, ierr)
         integer, intent(in) :: id, min_steps
         real(dp), intent(in) :: initial_mass_change, final_mass_change, relax_mass_change_max_yrs_dt
         integer, intent(out) :: ierr

         integer, parameter ::  lipar=2, lrpar=2
         integer :: max_model_number, i
         real(dp) :: max_years_for_timestep
         type (star_info), pointer :: s
         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:)
         real(dp), target :: rpar_ary(lrpar)
         real(dp), pointer :: rpar(:)
         logical :: adding_mass

         rpar => rpar_ary
         ipar => ipar_ary
         include 'formats'

         ierr = 0
         if (abs(initial_mass_change - final_mass_change) < 1d-10) return

         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return


         ipar(1) = s% model_number
         ipar(2) = min_steps
         rpar(1) = initial_mass_change
         rpar(2) = final_mass_change

         write(*,*)
         write(*,2) 'relax_mass_change -- min steps, init, final, max dt', &
            min_steps, initial_mass_change, final_mass_change, relax_mass_change_max_yrs_dt
         write(*,*)

         max_model_number = s% max_model_number
         s% max_model_number = -1111
         max_years_for_timestep = s% max_years_for_timestep
         s% max_years_for_timestep = relax_mass_change_max_yrs_dt

         call do_internal_evolve( &
               id, before_evolve_relax_mass_change, &
               relax_mass_change_adjust_model, relax_mass_change_check_model, &
               null_finish_model, .true., lipar, ipar, lrpar, rpar, ierr)

         s% max_model_number = max_model_number
         s% max_years_for_timestep = max_years_for_timestep

         if (ierr == 0) s% mass_change = final_mass_change

         call error_check('relax mass change',ierr)

      end subroutine do_relax_mass_change


      subroutine before_evolve_relax_mass_change(s, id, lipar, ipar, lrpar, rpar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         ierr = 0
         call turn_off_winds(s)
         s% max_model_number = -111
      end subroutine before_evolve_relax_mass_change

      integer function relax_mass_change_adjust_model(s, id, lipar, ipar, lrpar, rpar)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         relax_mass_change_adjust_model = keep_going
      end function relax_mass_change_adjust_model

      integer function relax_mass_change_check_model(s, id, lipar, ipar, lrpar, rpar)
         use do_one_utils, only:do_bare_bones_check_model
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)

         integer :: adjust_model, max_num_steps, num_steps
         real(dp) :: init_mass_change, final_mass_change, mass_change, frac
         logical, parameter :: dbg = .false.

         include 'formats'

         relax_mass_change_check_model = do_bare_bones_check_model(id)
         if (relax_mass_change_check_model /= keep_going) return

         adjust_model = ipar(1)
         max_num_steps = ipar(2)
         num_steps = s% model_number - adjust_model

         init_mass_change = rpar(1)
         final_mass_change = rpar(2)
         frac = dble(num_steps)/dble(max_num_steps)

         if (s% dt < s% max_years_for_timestep*secyer) then
            ipar(1) = adjust_model + 1 ! don't count this one
            write(*,'(a60,2i6,3x,99e12.3)') 'relax_mass_change wait for dt: model, step, frac', &
               s% model_number, num_steps, frac
            return
         end if

         if (frac >= 1d0) then
            s% mass_change = final_mass_change
            relax_mass_change_check_model = terminate
            s% termination_code = t_relax_finished_okay
            write(*,'(a60,2i6,3x,99e12.3)') 'DONE: relax_mass_change'
            return
         end if

         s% mass_change = frac*final_mass_change + (1-frac)*init_mass_change

         if (mod(s% model_number, s% terminal_interval) == 0) &
            write(*,'(a60,2i6,3x,99e12.3)') 'relax_mass_change, model, step, fraction done, mass_change', &
               s% model_number, num_steps, frac, s% mass_change

      end function relax_mass_change_check_model


      subroutine do_relax_core( &
            id, new_core_mass, dlg_core_mass_per_step, &
            relax_core_years_for_dt, core_avg_rho, core_avg_eps, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: new_core_mass ! in Msun units
         real(dp), intent(in) :: dlg_core_mass_per_step, relax_core_years_for_dt
         real(dp), intent(in) :: core_avg_rho, core_avg_eps
            ! adjust R_center according to core_avg_rho (g cm^-3)
            ! adjust L_center according to core_avg_eps (erg g^-1 s^-1)
         integer, intent(out) :: ierr

         integer, parameter ::  lipar=1, lrpar=5
         integer :: max_model_number
         real(dp) :: max_years_for_timestep
         real(dp) :: current_core_mass
         type (star_info), pointer :: s
         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:)
         real(dp), target :: rpar_ary(lrpar)
         real(dp), pointer :: rpar(:)
         rpar => rpar_ary
         ipar => ipar_ary
         include 'formats'
         ierr = 0
         if (new_core_mass <= 0) then
            ierr = -1
            write(*,*) 'invalid new_core_mass', new_core_mass
            return
         end if
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         current_core_mass = s% M_center/Msun
         if (abs(new_core_mass - current_core_mass) <= 1d-12*new_core_mass) then
            call do1_relax_core(s, new_core_mass, core_avg_rho, core_avg_eps, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in do1_relax_core'
            end if
            return
         end if
         write(*,*)
         write(*,1) 'current core mass', current_core_mass
         write(*,1) 'relax to new_core_mass', new_core_mass
         write(*,*)
         write(*,1) 'dlg_core_mass_per_step', dlg_core_mass_per_step
         write(*,*)
         rpar(1) = new_core_mass
         rpar(2) = dlg_core_mass_per_step
         rpar(3) = relax_core_years_for_dt
         rpar(4) = core_avg_rho
         rpar(5) = core_avg_eps
         max_model_number = s% max_model_number
         s% max_model_number = -1111
         max_years_for_timestep = s% max_years_for_timestep
         s% max_years_for_timestep = relax_core_years_for_dt
         call do_internal_evolve( &
               id, before_evolve_relax_core, relax_core_adjust_model, &
               relax_core_check_model, null_finish_model,  &
               .true., lipar, ipar, lrpar, rpar, ierr)
         s% max_model_number = max_model_number
         s% max_years_for_timestep = max_years_for_timestep
         call error_check('relax core',ierr)
      end subroutine do_relax_core


      subroutine before_evolve_relax_core(s, id, lipar, ipar, lrpar, rpar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         real(dp) :: relax_core_years_for_dt
         ierr = 0
         call setup_before_relax(s)
         s% max_model_number = -111
         relax_core_years_for_dt = rpar(3)
         s% dt_next = secyer*relax_core_years_for_dt
      end subroutine before_evolve_relax_core

      integer function relax_core_adjust_model(s, id, lipar, ipar, lrpar, rpar)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         relax_core_adjust_model = keep_going
      end function relax_core_adjust_model

      integer function relax_core_check_model(s, id, lipar, ipar, lrpar, rpar)
         use do_one_utils, only:do_bare_bones_check_model
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         real(dp) :: new_core_mass, dlg_core_mass_per_step, next
         real(dp) :: relax_core_dt, relax_core_years_for_dt
         real(dp) :: core_avg_rho, core_avg_eps, current_core_mass
         integer :: ierr
         logical :: end_now

         logical, parameter :: dbg = .false.

         include 'formats'

         relax_core_check_model = do_bare_bones_check_model(id)
         if (relax_core_check_model /= keep_going) return

         new_core_mass = rpar(1)
         dlg_core_mass_per_step = rpar(2)
         relax_core_years_for_dt = rpar(3)
         core_avg_rho = rpar(4)
         core_avg_eps = rpar(5)
         ierr = 0

         current_core_mass = s% M_center/Msun
         if (abs(new_core_mass - current_core_mass) <= 1d-12*new_core_mass) then
            call do1_relax_core(s, new_core_mass, core_avg_rho, core_avg_eps, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in do1_relax_core'
            end if
            relax_core_check_model = terminate
            s% termination_code = t_relax_finished_okay
            return
         end if

         if (mod(s% model_number, s% terminal_interval) == 0) then
            write(*,1) 'current & target core masses', &
               current_core_mass, new_core_mass, &
               (new_core_mass - current_core_mass)/new_core_mass
         end if

         relax_core_dt = secyer*relax_core_years_for_dt

         if (s% dt < relax_core_dt*0.9d0) then
            write(*,1) 's% dt < relax_core_dt*0.9d0', s% dt, relax_core_dt*0.9d0
            write(*,1) 's% max_timestep', s% max_timestep
            write(*,1) 's% max_years_for_timestep*secyer', s% max_years_for_timestep*secyer
            write(*,*)
            return ! give a chance to stabilize
         end if

         end_now=.false.
         if (new_core_mass < current_core_mass) then
            next = exp10(safe_log10(current_core_mass) - dlg_core_mass_per_step)
            if (next < new_core_mass) then
               next = new_core_mass
               end_now = .true.
            end if
         else
            next = exp10(log10(max(1d-8,current_core_mass)) + dlg_core_mass_per_step)
            if (next > new_core_mass) then
               next = new_core_mass
               end_now = .true.
            end if
         end if

         if (dbg) write(*,1) 'next', next, log10(next)

         call do1_relax_core(s, next, core_avg_rho, core_avg_eps, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in do1_relax_core'
            relax_core_check_model = terminate
         end if

         if (.true.) then
            write(*,1) 's% M_center', s% M_center
            write(*,1) 's% L_center', s% L_center
            write(*,1) 's% R_center', s% R_center
            write(*,1) 's% xmstar', s% xmstar
            write(*,*)
         end if
         
         if(end_now) then
            relax_core_check_model = terminate
            s% termination_code = t_relax_finished_okay
            return         
         end if

      end function relax_core_check_model


      subroutine do1_relax_core(s, next, core_avg_rho, core_avg_eps, ierr)
         type (star_info), pointer :: s
         real(dp), intent(in) :: next, core_avg_rho, core_avg_eps
         integer, intent(out) :: ierr
         real(dp) :: next_M_center, next_R_center, next_L_center
         ierr = 0
         next_M_center = next*Msun
         s% M_center = next_M_center
         s% xmstar = s% mstar - s% M_center
         next_R_center = pow(s% M_center/(core_avg_rho*four_thirds_pi),one_third)
         call do1_relax_R_center(s, next_R_center, ierr)
         if (ierr /= 0) return
         next_L_center = s% M_center*core_avg_eps
         call do1_relax_L_center(s, next_L_center, ierr)
      end subroutine do1_relax_core


      subroutine do_relax_mass_scale( &
            id, new_mass, dlgm_per_step, change_mass_years_for_dt, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: new_mass, dlgm_per_step, change_mass_years_for_dt
         integer, intent(out) :: ierr
         integer, parameter ::  lipar=1, lrpar=3
         integer :: max_model_number
         real(dp) :: max_years_for_timestep, relax_mass_scale_dt, eps_mdot_factor
         logical :: adding_mass
         type (star_info), pointer :: s
         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:)
         real(dp), target :: rpar_ary(lrpar)
         real(dp), pointer :: rpar(:)
         rpar => rpar_ary
         ipar => ipar_ary
         include 'formats'
         ierr = 0
         if (new_mass <= 0) then
            ierr = -1
            write(*,*) 'invalid new_mass', new_mass
            return
         end if
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (abs(new_mass - s% star_mass) <= 1d-12*new_mass) then
            s% star_mass = new_mass
            s% mstar = new_mass*Msun
            s% xmstar = s% mstar - s% M_center
            return
         end if
         write(*,*)
         write(*,1) 'relax_mass_scale'
         write(*,1) 'current star_mass', s% star_mass
         write(*,1) 'relax to new_mass', new_mass
         write(*,1) 'dlgm_per_step', dlgm_per_step
         write(*,*)
         rpar(1) = new_mass
         rpar(2) = dlgm_per_step
         rpar(3) = change_mass_years_for_dt
         max_model_number = s% max_model_number
         s% max_model_number = -1111

         adding_mass = (new_mass > s% star_mass)

         eps_mdot_factor = s% eps_mdot_factor
         s% eps_mdot_factor = 0

         max_years_for_timestep = s% max_years_for_timestep
         relax_mass_scale_dt = secyer*change_mass_years_for_dt
         s% max_years_for_timestep = relax_mass_scale_dt/secyer
         call do_internal_evolve( &
               id, before_evolve_relax_mass_scale, &
               relax_mass_scale_adjust_model, relax_mass_scale_check_model, null_finish_model, &
               .true., lipar, ipar, lrpar, rpar, ierr)
         s% max_model_number = max_model_number
         s% star_mass = new_mass
         s% mstar = new_mass*Msun
         s% xmstar = s% mstar - s% M_center
         s% eps_mdot_factor = eps_mdot_factor
         s% max_years_for_timestep = max_years_for_timestep

         call error_check('relax mass scale',ierr)
      end subroutine do_relax_mass_scale


      subroutine before_evolve_relax_mass_scale(s, id, lipar, ipar, lrpar, rpar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         real(dp) :: change_mass_years_for_dt
         ierr = 0
         call setup_before_relax(s)
         s% max_model_number = -111
         change_mass_years_for_dt = rpar(3)
         s% dt_next = secyer*change_mass_years_for_dt
      end subroutine before_evolve_relax_mass_scale

      integer function relax_mass_scale_adjust_model(s, id, lipar, ipar, lrpar, rpar)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         relax_mass_scale_adjust_model = keep_going
      end function relax_mass_scale_adjust_model

      integer function relax_mass_scale_check_model(s, id, lipar, ipar, lrpar, rpar)
         use do_one_utils, only:do_bare_bones_check_model
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         real(dp) :: new_mass, dlgm_per_step, next
         real(dp) :: relax_mass_scale_dt, change_mass_years_for_dt
         logical, parameter :: dbg = .false.

         include 'formats'

         relax_mass_scale_check_model = do_bare_bones_check_model(id)
         if (relax_mass_scale_check_model /= keep_going) return

         new_mass = rpar(1)
         dlgm_per_step = rpar(2)
         change_mass_years_for_dt = rpar(3)
         if (s% star_mass < 0.01d0) dlgm_per_step = dlgm_per_step*0.1d0
         !if (s% star_mass < 0.001d0) dlgm_per_step = dlgm_per_step*0.1d0

         if (dbg) then
            write(*,1) 'new_mass', new_mass
            write(*,1) 'current mass', s% star_mass
         end if

         if (abs(s% star_mass-new_mass) < 1d-12*new_mass) then
            s% star_mass = new_mass
            s% mstar = new_mass*Msun
            s% xmstar = s% mstar - s% M_center
            s% termination_code = t_relax_finished_okay
            relax_mass_scale_check_model = terminate
            return
         end if

         relax_mass_scale_dt = secyer*change_mass_years_for_dt

         if (s% dt < relax_mass_scale_dt*0.9d0) return ! give a chance to stabilize

         if (new_mass < s% star_mass) then
            next = exp10(safe_log10(s% star_mass) - dlgm_per_step)
            if (next < new_mass) next = new_mass
         else
            next = exp10(safe_log10(s% star_mass) + dlgm_per_step)
            if (next > new_mass) next = new_mass
         end if

         if (dbg) write(*,1) 'next', next, log10(next)

         s% star_mass = next
         s% mstar = next*Msun
         s% xmstar = s% mstar - s% M_center

      end function relax_mass_scale_check_model


      subroutine do_relax_M_center(id, new_mass, dlgm_per_step, relax_M_center_dt, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: new_mass, dlgm_per_step, relax_M_center_dt
         integer, intent(out) :: ierr
         integer, parameter ::  lipar=1, lrpar=3
         integer :: max_model_number
         real(dp) :: max_years_for_timestep
         type (star_info), pointer :: s
         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:)
         real(dp), target :: rpar_ary(lrpar)
         real(dp), pointer :: rpar(:)
         rpar => rpar_ary
         ipar => ipar_ary
         include 'formats'
         ierr = 0
         if (new_mass <= 0) then
            ierr = -1
            write(*,*) 'invalid new_mass', new_mass
            return
         end if
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (abs(new_mass - s% star_mass) <= 1d-6) then
            s% star_mass = new_mass
            s% mstar = new_mass*Msun
            s% xmstar = s% mstar - s% M_center
            return
         end if
         write(*,*)
         write(*,1) 'current star_mass', s% star_mass
         write(*,1) 'relax to new_mass', new_mass
         write(*,*)
         write(*,1) 'dlgm_per_step', dlgm_per_step
         write(*,*)
         rpar(1) = new_mass
         rpar(2) = dlgm_per_step
         rpar(3) = relax_M_center_dt
         max_model_number = s% max_model_number
         max_years_for_timestep = s% max_years_for_timestep
         s% max_model_number = -1111
         s% max_years_for_timestep = relax_M_center_dt/secyer
         call do_internal_evolve( &
               id, before_evolve_relax_M_center, &
               relax_M_center_adjust_model, relax_M_center_check_model, &
               null_finish_model, .true., lipar, ipar, lrpar, rpar, ierr)

         s% max_model_number = max_model_number
         s% max_years_for_timestep = max_years_for_timestep
         call error_check('relax M center',ierr)
      end subroutine do_relax_M_center


      subroutine before_evolve_relax_M_center(s, id, lipar, ipar, lrpar, rpar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         ierr = 0
         call setup_before_relax(s)
         s% max_model_number = -111
         s% dt_next =  rpar(3) ! relax_M_center_dt
      end subroutine before_evolve_relax_M_center

      integer function relax_M_center_adjust_model(s, id, lipar, ipar, lrpar, rpar)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         relax_M_center_adjust_model = keep_going
      end function relax_M_center_adjust_model

      integer function relax_M_center_check_model(s, id, lipar, ipar, lrpar, rpar)
         use do_one_utils, only:do_bare_bones_check_model
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer :: ierr
         real(dp) :: new_mass, dlgm_per_step, relax_M_center_dt, next
         logical, parameter :: dbg = .false.
         logical :: end_now

         include 'formats'

         relax_M_center_check_model = do_bare_bones_check_model(id)
         if (relax_M_center_check_model /= keep_going) return

         new_mass = rpar(1)
         dlgm_per_step = rpar(2)
         relax_M_center_dt = rpar(3)

         if (mod(s% model_number, s% terminal_interval) == 0) &
            write(*,1) 'relax_M_center target/current', new_mass/(s% M_center/Msun)
         
         end_now=.false.
         if (new_mass < s% star_mass) then
            next = exp10(safe_log10(s% star_mass) - dlgm_per_step)
            if (next < new_mass) then
               next = new_mass
               end_now=.true.
            end if
         else
            next = exp10(safe_log10(s% star_mass) + dlgm_per_step)
            if (next > new_mass) then
               next = new_mass
               end_now=.true.
            end if
         end if

         if (abs(s% star_mass - new_mass) < 1d-15 .or. end_now) then
            call set_new_mass_for_relax_M_center(s, new_mass, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed to set mass for relax_M_center'
               relax_M_center_check_model = terminate
               return
            end if
            write(*,1) 'final mass', s% star_mass, s% mstar, s% M_center, s% xmstar
            relax_M_center_check_model = terminate
            s% termination_code = t_relax_finished_okay
            return
         end if

         if (s% dt < relax_M_center_dt*0.9d0) return ! give a chance to stabilize

         if (dbg) write(*,1) 'next', next, log10(next)

         call set_new_mass_for_relax_M_center(s, next, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to set mass for relax_M_center'
            relax_M_center_check_model = terminate
            return
         end if

      end function relax_M_center_check_model


      subroutine set_new_mass_for_relax_M_center(s, new_mass, ierr)
         use star_utils, only: set_qs
         type (star_info), pointer :: s
         real(dp), intent(in) :: new_mass ! Msun
         integer, intent(out) :: ierr
         include 'formats'
         ierr = 0
         s% star_mass = new_mass
         s% mstar = new_mass*Msun
         s% M_center = s% mstar - s% xmstar
      end subroutine set_new_mass_for_relax_M_center


      subroutine do_relax_R_center(id, new_R_center, dlgR_per_step, relax_R_center_dt, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: new_R_center, dlgR_per_step, relax_R_center_dt
         integer, intent(out) :: ierr
         integer, parameter ::  lipar=1, lrpar=3
         integer :: max_model_number
         real(dp) :: max_years_for_timestep
         type (star_info), pointer :: s
         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:)
         real(dp), target :: rpar_ary(lrpar)
         real(dp), pointer :: rpar(:)
         rpar => rpar_ary
         ipar => ipar_ary
         include 'formats'
         ierr = 0
         if (new_R_center < 0) then
            ierr = -1
            write(*,*) 'invalid new_R_center', new_R_center
            return
         end if
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (abs(new_R_center - s% R_center) <= 1d-6) then
            s% R_center = new_R_center
            return
         end if
         write(*,*)
         write(*,1) 'current R_center', s% R_center
         write(*,1) 'relax to new_R_center', new_R_center
         write(*,*)
         write(*,1) 'dlgR_per_step', dlgR_per_step
         write(*,*)
         rpar(1) = new_R_center
         rpar(2) = dlgR_per_step
         rpar(3) = relax_R_center_dt
         max_model_number = s% max_model_number
         max_years_for_timestep = s% max_years_for_timestep
         s% max_model_number = -1111
         s% max_years_for_timestep = relax_R_center_dt/secyer
         call do_internal_evolve( &
               id, before_evolve_relax_R_center, &
               relax_R_center_adjust_model, relax_R_center_check_model, &
               null_finish_model, .true., lipar, ipar, lrpar, rpar, ierr)

         s% max_model_number = max_model_number
         s% max_years_for_timestep = max_years_for_timestep
         call error_check('relax R center',ierr)
      end subroutine do_relax_R_center


      subroutine before_evolve_relax_R_center(s, id, lipar, ipar, lrpar, rpar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         ierr = 0
         call setup_before_relax(s)
         s% max_model_number = -111
         s% dt_next =  rpar(3) ! relax_R_center_dt
      end subroutine before_evolve_relax_R_center

      integer function relax_R_center_adjust_model(s, id, lipar, ipar, lrpar, rpar)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         relax_R_center_adjust_model = keep_going
      end function relax_R_center_adjust_model

      integer function relax_R_center_check_model(s, id, lipar, ipar, lrpar, rpar)
         use do_one_utils, only:do_bare_bones_check_model
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer :: ierr
         real(dp) :: new_R_center, dlgR_per_step, relax_R_center_dt, next
         logical, parameter :: dbg = .false.

         include 'formats'

         relax_R_center_check_model = do_bare_bones_check_model(id)
         if (relax_R_center_check_model /= keep_going) return

         new_R_center = rpar(1)
         dlgR_per_step = rpar(2)
         relax_R_center_dt = rpar(3)

         if (mod(s% model_number, s% terminal_interval) == 0) &
            write(*,1) 'relax_R_center target/current', new_R_center/s% R_center

         if (abs(s% R_center - new_R_center) < 1d-15) then
            call do1_relax_R_center(s, new_R_center, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in relax_R_center'
            end if
            relax_R_center_check_model = terminate
            s% termination_code = t_relax_finished_okay
            return
         end if

         if (s% dt < relax_R_center_dt*0.9d0) return ! give a chance to stabilize

         if (new_R_center < s% R_center) then
            next = exp10(safe_log10(s% R_center) - dlgR_per_step)
            if (next < new_R_center) next = new_R_center
         else if (s% R_center < 1) then
            next = 1
         else
            next = exp10(safe_log10(s% R_center) + dlgR_per_step)
            if (next > new_R_center) next = new_R_center
         end if

         if (dbg) write(*,1) 'next', next

         call do1_relax_R_center(s, next, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in relax_R_center'
            relax_R_center_check_model = terminate
         end if

      end function relax_R_center_check_model


      subroutine do1_relax_R_center(s, new_Rcenter, ierr)
         ! adjust all lnR's to keep same density for each cell as 1st guess for next model
         use star_utils, only: set_qs
         type (star_info), pointer :: s
         real(dp), intent(in) :: new_Rcenter ! cm
         integer, intent(out) :: ierr
         real(dp) :: dm, rho, dr3, rp13
         integer :: k
         include 'formats'
         ierr = 0
         s% R_center = new_Rcenter
         ! adjust lnR's
         rp13 = s% R_center*s% R_center*s% R_center
         do k = s% nz, 1, -1
            dm = s% dm(k)
            rho = s% rho(k)
            dr3 = dm/(rho*four_thirds_pi) ! dm/rho is cell volume
            s% xh(s% i_lnR,k) = log(rp13 + dr3)*one_third
            rp13 = rp13 + dr3
         end do
      end subroutine do1_relax_R_center


      subroutine do_relax_v_center(id, new_v_center, dv_per_step, relax_v_center_dt, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: new_v_center, dv_per_step, relax_v_center_dt
         integer, intent(out) :: ierr
         integer, parameter ::  lipar=1, lrpar=3
         integer :: max_model_number
         real(dp) :: max_years_for_timestep
         type (star_info), pointer :: s
         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:)
         real(dp), target :: rpar_ary(lrpar)
         real(dp), pointer :: rpar(:)
         rpar => rpar_ary
         ipar => ipar_ary
         include 'formats'
         ierr = 0
         if (new_v_center < 0) then
            ierr = -1
            write(*,*) 'invalid new_v_center', new_v_center
            return
         end if
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (abs(s% v_center - new_v_center) < &
               1d-6*max(1d-6,abs(s% v_center),abs(new_v_center))) then
            s% v_center = new_v_center
            return
         end if
         write(*,*)
         write(*,1) 'current v_center', s% v_center
         write(*,1) 'relax to new_v_center', new_v_center
         write(*,*)
         write(*,1) 'dv_per_step', dv_per_step
         write(*,*)
         rpar(1) = new_v_center
         rpar(2) = dv_per_step
         rpar(3) = relax_v_center_dt
         max_model_number = s% max_model_number
         max_years_for_timestep = s% max_years_for_timestep
         s% max_model_number = -1111
         s% max_years_for_timestep = relax_v_center_dt/secyer
         call do_internal_evolve( &
               id, before_evolve_relax_v_center, &
               relax_v_center_adjust_model, relax_v_center_check_model, &
               null_finish_model, .true., lipar, ipar, lrpar, rpar, ierr)

         s% max_model_number = max_model_number
         s% max_years_for_timestep = max_years_for_timestep
         call error_check('relax V center',ierr)
      end subroutine do_relax_v_center


      subroutine before_evolve_relax_v_center(s, id, lipar, ipar, lrpar, rpar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         ierr = 0
         call setup_before_relax(s)
         s% max_model_number = -111
         s% dt_next =  rpar(3) ! relax_v_center_dt
      end subroutine before_evolve_relax_v_center

      integer function relax_v_center_adjust_model(s, id, lipar, ipar, lrpar, rpar)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         relax_v_center_adjust_model = keep_going
      end function relax_v_center_adjust_model

      integer function relax_v_center_check_model(s, id, lipar, ipar, lrpar, rpar)
         use do_one_utils, only:do_bare_bones_check_model
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer :: ierr
         real(dp) :: new_v_center, dv_per_step, relax_v_center_dt, next
         logical, parameter :: dbg = .false.

         include 'formats'

         relax_v_center_check_model = do_bare_bones_check_model(id)
         if (relax_v_center_check_model /= keep_going) return

         new_v_center = rpar(1)
         dv_per_step = rpar(2)
         relax_v_center_dt = rpar(3)

         if (mod(s% model_number, s% terminal_interval) == 0) &
            write(*,1) 'target v_center current', new_v_center, s% v_center

         if (abs(s% v_center - new_v_center) < &
               1d-6*max(1d-6,abs(s% v_center),abs(new_v_center))) then
            s% v_center = new_v_center
            relax_v_center_check_model = terminate
            s% termination_code = t_relax_finished_okay
            return
         end if

         if (s% dt < relax_v_center_dt*0.9d0) return ! give a chance to stabilize

         if (new_v_center < s% v_center) then
            next = s% v_center - dv_per_step
            if (next < new_v_center) next = new_v_center
         else
            next = s% v_center + dv_per_step
            if (next > new_v_center) next = new_v_center
         end if

         if (dbg) write(*,1) 'next', next

         s% v_center = next

      end function relax_v_center_check_model


      subroutine do_relax_L_center(id, new_L_center, dlgL_per_step, relax_L_center_dt, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: new_L_center, dlgL_per_step, relax_L_center_dt
         integer, intent(out) :: ierr
         integer, parameter ::  lipar=1, lrpar=3
         integer :: max_model_number
         real(dp) :: max_years_for_timestep
         type (star_info), pointer :: s
         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:)
         real(dp), target :: rpar_ary(lrpar)
         real(dp), pointer :: rpar(:)
         rpar => rpar_ary
         ipar => ipar_ary
         include 'formats'
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (abs(new_L_center - s% L_center) <= &
               1d-10*max(abs(new_L_center),abs(s% L_center),1d0)) then
            s% L_center = new_L_center
            return
         end if
         write(*,*)
         write(*,1) 'current L_center', s% L_center
         write(*,1) 'relax to new_L_center', new_L_center
         write(*,*)
         write(*,1) 'dlgL_per_step', dlgL_per_step
         write(*,*)
         rpar(1) = new_L_center
         rpar(2) = dlgL_per_step*(new_L_center - s% L_center)
         rpar(3) = relax_L_center_dt
         max_model_number = s% max_model_number
         max_years_for_timestep = s% max_years_for_timestep
         s% max_model_number = -1111
         s% max_years_for_timestep = relax_L_center_dt/secyer
         call do_internal_evolve( &
               id, before_evolve_relax_L_center, &
               relax_L_center_adjust_model, relax_L_center_check_model, &
               null_finish_model, .true., lipar, ipar, lrpar, rpar, ierr)

         s% max_model_number = max_model_number
         s% max_years_for_timestep = max_years_for_timestep
         call error_check('relax L center',ierr)
      end subroutine do_relax_L_center


      subroutine before_evolve_relax_L_center(s, id, lipar, ipar, lrpar, rpar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         ierr = 0
         call setup_before_relax(s)
         s% max_model_number = -111
         s% dt_next =  rpar(3) ! relax_L_center_dt
      end subroutine before_evolve_relax_L_center

      integer function relax_L_center_adjust_model(s, id, lipar, ipar, lrpar, rpar)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         relax_L_center_adjust_model = keep_going
      end function relax_L_center_adjust_model

      integer function relax_L_center_check_model(s, id, lipar, ipar, lrpar, rpar)
         use do_one_utils, only:do_bare_bones_check_model
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer :: ierr
         real(dp) :: new_L_center, dL, relax_L_center_dt, next
         logical, parameter :: dbg = .false.

         include 'formats'

         relax_L_center_check_model = do_bare_bones_check_model(id)
         if (relax_L_center_check_model /= keep_going) return

         new_L_center = rpar(1)
         dL = rpar(2)
         relax_L_center_dt = rpar(3)


         if (mod(s% model_number, s% terminal_interval) == 0) &
            write(*,1) 'relax_L_center target/current', new_L_center/s% L_center

         if (abs(new_L_center - s% L_center) < abs(dL)) then
            call do1_relax_L_center(s, new_L_center, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in relax_L_center'
            end if
            relax_L_center_check_model = terminate
            s% termination_code = t_relax_finished_okay
            return
         end if

         if (s% dt < relax_L_center_dt*0.9d0) return ! give a chance to stabilize

         next = s% L_center + dL
         if (dbg) write(*,1) 'next', next

         call do1_relax_L_center(s, next, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in relax_L_center'
            relax_L_center_check_model = terminate
         end if

      end function relax_L_center_check_model


      subroutine do1_relax_L_center(s, new_Lcenter, ierr)
         type (star_info), pointer :: s
         real(dp), intent(in) :: new_Lcenter
         integer, intent(out) :: ierr
         real(dp) :: L_center_prev, dL
         integer :: i_lum, nz, k
         include 'formats'
         ierr = 0
         nz = s% nz
         L_center_prev = s% L_center
         s% L_center = new_Lcenter
         i_lum = s% i_lum
         if (i_lum == 0) return
         dL = new_Lcenter - L_center_prev
         do k=1,nz
            s% xh(i_lum,k) = s% xh(i_lum,k) + dL
         end do
      end subroutine do1_relax_L_center


      subroutine do_relax_dxdt_nuc_factor(id, new_value, per_step_multiplier, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: new_value, per_step_multiplier
         integer, intent(out) :: ierr
         integer, parameter ::  lipar=1, lrpar=2
         integer :: max_model_number
         type (star_info), pointer :: s
         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:)
         real(dp), target :: rpar_ary(lrpar)
         real(dp), pointer :: rpar(:)
         rpar => rpar_ary
         ipar => ipar_ary
         include 'formats'
         ierr = 0
         if (new_value <= 0) then
            ierr = -1
            write(*,*) 'invalid new_value', new_value
            return
         end if
         if (per_step_multiplier <= 0 .or. per_step_multiplier == 1) then
            ierr = -1
            write(*,*) 'invalid per_step_multiplier', per_step_multiplier
            return
         end if
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (abs(new_value - s% dxdt_nuc_factor) <= 1d-6) then
            s% dxdt_nuc_factor = new_value
            return
         end if
         write(*,*)
         write(*,1) 'current dxdt_nuc_factor', s% dxdt_nuc_factor
         write(*,1) 'relax to new_value', new_value
         write(*,*)
         write(*,1) 'per_step_multiplier', per_step_multiplier
         write(*,*)
         rpar(1) = new_value
         rpar(2) = per_step_multiplier
         max_model_number = s% max_model_number
         s% max_model_number = -1111
         call do_internal_evolve( &
               id, before_evolve_relax_dxdt_nuc_factor, &
               relax_dxdt_nuc_factor_adjust_model, relax_dxdt_nuc_factor_check_model, &
               null_finish_model, .true., lipar, ipar, lrpar, rpar, ierr)
         s% max_model_number = max_model_number
         s% dxdt_nuc_factor = new_value
         call error_check('relax dxdt nuc factor',ierr)
      end subroutine do_relax_dxdt_nuc_factor


      subroutine before_evolve_relax_dxdt_nuc_factor(s, id, lipar, ipar, lrpar, rpar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         ierr = 0
         call turn_off_winds(s)
         s% max_model_number = -111
         s% dt_next = secyer*1d-3
      end subroutine before_evolve_relax_dxdt_nuc_factor

      integer function relax_dxdt_nuc_factor_adjust_model(s, id, lipar, ipar, lrpar, rpar)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         relax_dxdt_nuc_factor_adjust_model = keep_going
      end function relax_dxdt_nuc_factor_adjust_model

      integer function relax_dxdt_nuc_factor_check_model(s, id, lipar, ipar, lrpar, rpar)
         use do_one_utils, only:do_bare_bones_check_model
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer :: ierr
         real(dp) :: new_value, per_step_multiplier
         logical, parameter :: dbg = .false.

         include 'formats'

         relax_dxdt_nuc_factor_check_model = do_bare_bones_check_model(id)
         if (relax_dxdt_nuc_factor_check_model /= keep_going) return

         new_value = rpar(1)
         per_step_multiplier = rpar(2)

         if (mod(s% model_number, s% terminal_interval) == 0) &
            write(*,1) 'new_value current', new_value, s% dxdt_nuc_factor

         s% dxdt_nuc_factor = s% dxdt_nuc_factor * per_step_multiplier

         if ((per_step_multiplier < 1 .and. s% dxdt_nuc_factor < new_value) .or. &
             (per_step_multiplier > 1 .and. s% dxdt_nuc_factor > new_value)) then
            s% dxdt_nuc_factor = new_value
            relax_dxdt_nuc_factor_check_model = terminate
            s% termination_code = t_relax_finished_okay
            return
         end if

      end function relax_dxdt_nuc_factor_check_model


      subroutine do_relax_eps_nuc_factor(id, new_value, per_step_multiplier, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: new_value, per_step_multiplier
         integer, intent(out) :: ierr
         integer, parameter ::  lipar=1, lrpar=2
         integer :: max_model_number
         type (star_info), pointer :: s
         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:)
         real(dp), target :: rpar_ary(lrpar)
         real(dp), pointer :: rpar(:)
         rpar => rpar_ary
         ipar => ipar_ary
         include 'formats'
         ierr = 0
         if (new_value <= 0) then
            ierr = -1
            write(*,*) 'invalid new_value', new_value
            return
         end if
         if (per_step_multiplier <= 0 .or. per_step_multiplier == 1) then
            ierr = -1
            write(*,*) 'invalid per_step_multiplier', per_step_multiplier
            return
         end if
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (abs(new_value - s% eps_nuc_factor) <= 1d-6) then
            s% eps_nuc_factor = new_value
            return
         end if
         write(*,*)
         write(*,1) 'current eps_nuc_factor', s% eps_nuc_factor
         write(*,1) 'relax to new_value', new_value
         write(*,*)
         write(*,1) 'per_step_multiplier', per_step_multiplier
         write(*,*)
         rpar(1) = new_value
         rpar(2) = per_step_multiplier
         max_model_number = s% max_model_number
         s% max_model_number = -1111
         call do_internal_evolve( &
               id, before_evolve_relax_eps_nuc_factor, &
               relax_eps_nuc_factor_adjust_model, relax_eps_nuc_factor_check_model, &
               null_finish_model, .true., lipar, ipar, lrpar, rpar, ierr)
         s% max_model_number = max_model_number
         s% eps_nuc_factor = new_value
         call error_check('relax eps nuc factor',ierr)
      end subroutine do_relax_eps_nuc_factor


      subroutine before_evolve_relax_eps_nuc_factor(s, id, lipar, ipar, lrpar, rpar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         ierr = 0
         call turn_off_winds(s)
         s% max_model_number = -111
         s% dt_next = secyer*1d-3
      end subroutine before_evolve_relax_eps_nuc_factor

      integer function relax_eps_nuc_factor_adjust_model(s, id, lipar, ipar, lrpar, rpar)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         relax_eps_nuc_factor_adjust_model = keep_going
      end function relax_eps_nuc_factor_adjust_model

      integer function relax_eps_nuc_factor_check_model(s, id, lipar, ipar, lrpar, rpar)
         use do_one_utils, only:do_bare_bones_check_model
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer :: ierr
         real(dp) :: new_value, per_step_multiplier
         logical, parameter :: dbg = .false.

         include 'formats'

         relax_eps_nuc_factor_check_model = do_bare_bones_check_model(id)
         if (relax_eps_nuc_factor_check_model /= keep_going) return

         new_value = rpar(1)
         per_step_multiplier = rpar(2)

         if (mod(s% model_number, s% terminal_interval) == 0) &
            write(*,1) 'new_value, current', new_value, s% eps_nuc_factor

         s% eps_nuc_factor = s% eps_nuc_factor * per_step_multiplier

         if ((per_step_multiplier < 1 .and. s% eps_nuc_factor < new_value) .or. &
             (per_step_multiplier > 1 .and. s% eps_nuc_factor > new_value)) then
            s% eps_nuc_factor = new_value
            relax_eps_nuc_factor_check_model = terminate
            s% termination_code = t_relax_finished_okay
            return
         end if

      end function relax_eps_nuc_factor_check_model


      subroutine do_relax_opacity_max(id, new_value, per_step_multiplier, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: new_value, per_step_multiplier
         integer, intent(out) :: ierr
         integer, parameter ::  lipar=1, lrpar=2
         integer :: max_model_number
         type (star_info), pointer :: s
         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:)
         real(dp), target :: rpar_ary(lrpar)
         real(dp), pointer :: rpar(:)
         rpar => rpar_ary
         ipar => ipar_ary
         include 'formats'
         ierr = 0
         if (new_value <= 0) then
            ierr = -1
            write(*,*) 'invalid new_value', new_value
            return
         end if
         if (per_step_multiplier <= 0 .or. per_step_multiplier == 1) then
            ierr = -1
            write(*,*) 'invalid per_step_multiplier', per_step_multiplier
            return
         end if
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (s% opacity_max <= 0) then
            ierr = -1
            write(*,*) 'invalid opacity_max', s% opacity_max
            return
         end if
         if (abs(new_value - s% opacity_max) <= 1d-6) then
            s% opacity_max = new_value
            return
         end if
         write(*,*)
         write(*,1) 'current opacity_max', s% opacity_max
         write(*,1) 'relax to new_value', new_value
         write(*,*)
         write(*,1) 'per_step_multiplier', per_step_multiplier
         write(*,*)
         rpar(1) = new_value
         rpar(2) = per_step_multiplier
         max_model_number = s% max_model_number
         s% max_model_number = -1111
         call do_internal_evolve( &
               id, before_evolve_relax_opacity_max, &
               relax_opacity_max_adjust_model, relax_opacity_max_check_model, &
               null_finish_model, .true., lipar, ipar, lrpar, rpar, ierr)
         
         s% max_model_number = max_model_number
         s% opacity_max = new_value
         s% dt_next = rpar(1) ! keep dt from relax
         call error_check('relax opacity max',ierr)
      end subroutine do_relax_opacity_max


      subroutine before_evolve_relax_opacity_max(s, id, lipar, ipar, lrpar, rpar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         ierr = 0
         s% max_model_number = -111
         s% dt_next = secyer*1d-3
         call turn_off_winds(s)
      end subroutine before_evolve_relax_opacity_max

      integer function relax_opacity_max_adjust_model(s, id, lipar, ipar, lrpar, rpar)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         relax_opacity_max_adjust_model = keep_going
      end function relax_opacity_max_adjust_model

      integer function relax_opacity_max_check_model(s, id, lipar, ipar, lrpar, rpar)
         use do_one_utils, only:do_bare_bones_check_model
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         real(dp) :: new_value, per_step_multiplier
         logical, parameter :: dbg = .false.

         include 'formats'

         relax_opacity_max_check_model = do_bare_bones_check_model(id)
         if (relax_opacity_max_check_model /= keep_going) return

         new_value = rpar(1)
         per_step_multiplier = rpar(2)

         s% opacity_max = s% opacity_max * per_step_multiplier

         if (mod(s% model_number, s% terminal_interval) == 0) &
            write(*,1) 'relax opacity', s% opacity_max, new_value

         if ((per_step_multiplier < 1 .and. s% opacity_max < new_value) .or. &
             (per_step_multiplier > 1 .and. s% opacity_max > new_value)) then
            s% opacity_max = new_value
            relax_opacity_max_check_model = terminate
            s% termination_code = t_relax_finished_okay
            rpar(1) = s% dt
            return
         end if

      end function relax_opacity_max_check_model


      subroutine do_relax_fixed_L(id, steps, ierr)
         integer, intent(in) :: id
         integer, intent(in) :: steps
         integer, intent(out) :: ierr
         integer, parameter ::  lipar=1, lrpar=1
         integer :: max_model_number
         type (star_info), pointer :: s
         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:)
         real(dp), target :: rpar_ary(lrpar)
         real(dp), pointer :: rpar(:)
         rpar => rpar_ary
         ipar => ipar_ary
         include 'formats'
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (.not. s% use_fixed_L_for_BB_outer_BC) return
         if (steps <= 0) return
         ipar(1) = steps
         max_model_number = s% max_model_number
         s% max_model_number = -1111
         call do_internal_evolve( &
               id, before_evolve_relax_fixed_L, &
               relax_fixed_L_adjust_model, relax_fixed_L_check_model, &
               null_finish_model, .true., lipar, ipar, lrpar, rpar, ierr)

         s% max_model_number = max_model_number
         s% use_fixed_L_for_BB_outer_BC = .false.
         s% dt_next = rpar(1) ! keep dt from relax
         call error_check('relax fixed L',ierr)
      end subroutine do_relax_fixed_L


      subroutine before_evolve_relax_fixed_L(s, id, lipar, ipar, lrpar, rpar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         ierr = 0
         s% max_model_number = -111
         s% dt_next = secyer*1d-3
         call turn_off_winds(s)
      end subroutine before_evolve_relax_fixed_L

      integer function relax_fixed_L_adjust_model(s, id, lipar, ipar, lrpar, rpar)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         relax_fixed_L_adjust_model = keep_going
      end function relax_fixed_L_adjust_model

      integer function relax_fixed_L_check_model(s, id, lipar, ipar, lrpar, rpar)
         use do_one_utils, only: do_bare_bones_check_model
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer :: steps_remaining
         real(dp) :: f
         logical, parameter :: dbg = .false.

         include 'formats'

         relax_fixed_L_check_model = do_bare_bones_check_model(id)
         if (relax_fixed_L_check_model /= keep_going) return

         steps_remaining = ipar(1)
         
         if (steps_remaining <= 0) then
            relax_fixed_L_check_model = terminate
            s% termination_code = t_relax_finished_okay
            rpar(1) = s% dt
            return
         end if
         
         f = 1d0/dble(steps_remaining)
         steps_remaining = steps_remaining - 1
         ipar(1) = steps_remaining         
         s% fixed_L_for_BB_outer_BC = f*s% L(1) + (1d0 - f)*s% fixed_L_for_BB_outer_BC
         if (mod(s% model_number, s% terminal_interval) == 0) &
            write(*,1) 's% fixed_L_for_BB_outer_BC/Lsun', s% fixed_L_for_BB_outer_BC/Lsun

      end function relax_fixed_L_check_model


      subroutine do_relax_max_surf_dq(id, new_value, per_step_multiplier, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: new_value, per_step_multiplier
         integer, intent(out) :: ierr
         integer, parameter ::  lipar=1, lrpar=2
         integer :: max_model_number
         type (star_info), pointer :: s
         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:)
         real(dp), target :: rpar_ary(lrpar)
         real(dp), pointer :: rpar(:)
         rpar => rpar_ary
         ipar => ipar_ary
         include 'formats'
         ierr = 0
         if (new_value <= 0) then
            ierr = -1
            write(*,*) 'invalid new_value', new_value
            return
         end if
         if (per_step_multiplier <= 0 .or. per_step_multiplier == 1) then
            ierr = -1
            write(*,*) 'invalid per_step_multiplier', per_step_multiplier
            return
         end if
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (s% max_surface_cell_dq <= 0) then
            ierr = -1
            write(*,*) 'invalid max_surf_dq', s% max_surface_cell_dq
            return
         end if
         if (abs(new_value - s% max_surface_cell_dq) <= &
               1d-6*min(new_value,s% max_surface_cell_dq)) then
            s% max_surface_cell_dq = new_value
            return
         end if
         write(*,*)
         write(*,1) 'current max_surf_dq', s% max_surface_cell_dq
         write(*,1) 'relax to new_value', new_value
         write(*,*)
         write(*,1) 'per_step_multiplier', per_step_multiplier
         write(*,*)
         rpar(1) = new_value
         rpar(2) = per_step_multiplier
         max_model_number = s% max_model_number
         s% max_model_number = -1111
         call do_internal_evolve( &
               id, before_evolve_relax_max_surf_dq, &
               relax_max_surf_dq_adjust_model, relax_max_surf_dq_check_model, &
               null_finish_model, .true., lipar, ipar, lrpar, rpar, ierr)
         s% max_model_number = max_model_number
         s% max_surface_cell_dq = new_value
         s% dt_next = rpar(1) ! keep dt from relax
         call error_check('relax max surf dq',ierr)
      end subroutine do_relax_max_surf_dq


      subroutine before_evolve_relax_max_surf_dq(s, id, lipar, ipar, lrpar, rpar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         ierr = 0
         s% max_model_number = -111
         s% dt_next = secyer*1d-3
         call turn_off_winds(s)
      end subroutine before_evolve_relax_max_surf_dq

      integer function relax_max_surf_dq_adjust_model(s, id, lipar, ipar, lrpar, rpar)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         relax_max_surf_dq_adjust_model = keep_going
      end function relax_max_surf_dq_adjust_model

      integer function relax_max_surf_dq_check_model(s, id, lipar, ipar, lrpar, rpar)
         use do_one_utils, only:do_bare_bones_check_model
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         real(dp) :: new_value, per_step_multiplier
         logical, parameter :: dbg = .false.

         include 'formats'

         relax_max_surf_dq_check_model = do_bare_bones_check_model(id)
         if (relax_max_surf_dq_check_model /= keep_going) return

         new_value = rpar(1)
         per_step_multiplier = rpar(2)

         s% max_surface_cell_dq = s% max_surface_cell_dq * per_step_multiplier

         if (mod(s% model_number, s% terminal_interval) == 0) &
           write(*,1) 'relax max_surface_cell_dq', s% max_surface_cell_dq, new_value

         if ((per_step_multiplier < 1 .and. s% max_surface_cell_dq < new_value) .or. &
             (per_step_multiplier > 1 .and. s% max_surface_cell_dq > new_value)) then
            s% max_surface_cell_dq = new_value
            relax_max_surf_dq_check_model = terminate
            s% termination_code = t_relax_finished_okay
            rpar(1) = s% dt
            return
         end if

      end function relax_max_surf_dq_check_model


      subroutine do_relax_num_steps(id, num_steps, max_timestep, ierr)
         integer, intent(in) :: id, num_steps
         real(dp), intent(in) :: max_timestep
         integer, intent(out) :: ierr
         integer, parameter ::  lipar=1, lrpar=1
         integer :: max_model_number, model_number
         real(dp) :: save_max_timestep, save_max_years_for_timestep
         type (star_info), pointer :: s
         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:)
         real(dp), target :: rpar_ary(lrpar)
         real(dp), pointer :: rpar(:)
         rpar => rpar_ary
         ipar => ipar_ary

         include 'formats'
         ierr = 0
         if (num_steps <= 0) return
         

         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         write(*,*)
         write(*,2) 'relax_num_steps', num_steps
         write(*,*)
         ipar(1) = num_steps
         if (max_timestep <= 0) then
            rpar(1) = secyer
         else
            rpar(1) = max_timestep
         end if
         max_model_number = s% max_model_number
         model_number = s% model_number
         save_max_timestep = s% max_timestep
         save_max_years_for_timestep = s% max_years_for_timestep
         s% model_number = 0
         call do_internal_evolve( &
               id, before_evolve_relax_num_steps, &
               relax_num_steps_adjust_model, relax_num_steps_check_model, &
               null_finish_model, .true., lipar, ipar, lrpar, rpar, ierr)
         s% max_model_number = max_model_number
         s% model_number = model_number
         s% max_timestep = save_max_timestep
         s% max_years_for_timestep = save_max_years_for_timestep
         call error_check('relax num steps',ierr)

      end subroutine do_relax_num_steps


      subroutine before_evolve_relax_num_steps(s, id, lipar, ipar, lrpar, rpar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         ierr = 0
         call setup_before_relax(s)
         s% max_timestep = rpar(1)
         s% max_years_for_timestep = s% max_timestep/secyer
         s% dt_next = s% max_timestep
         s% max_model_number = ipar(1)
      end subroutine before_evolve_relax_num_steps

      integer function relax_num_steps_adjust_model(s, id, lipar, ipar, lrpar, rpar)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         relax_num_steps_adjust_model = keep_going
      end function relax_num_steps_adjust_model

      integer function relax_num_steps_check_model(s, id, lipar, ipar, lrpar, rpar)
         use do_one_utils, only:do_bare_bones_check_model
         
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer :: ierr, klo, khi
         real(dp) :: lnbc_target, new_pre_ms, new_lnbc, dlnbc, lnbc, &
            current_pre_ms, next_pre_ms

         logical, parameter :: dbg = .false.

         include 'formats'

         relax_num_steps_check_model = do_bare_bones_check_model(id)
         if (relax_num_steps_check_model /= keep_going) return
         if (s% model_number >= ipar(1)) then
            relax_num_steps_check_model = terminate
            s% termination_code = t_relax_finished_okay
         end if

      end function relax_num_steps_check_model


      subroutine do_relax_to_radiative_core(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         integer, parameter ::  lipar=1, lrpar=2
         integer :: max_model_number, model_number
         real(dp) :: save_max_timestep, save_max_years_for_timestep
         type (star_info), pointer :: s
         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:)
         real(dp), target :: rpar_ary(lrpar)
         real(dp), pointer :: rpar(:)
         real(dp) :: max_timestep
         rpar => rpar_ary
         ipar => ipar_ary

         include 'formats'
         ierr = 0

         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         if (s% star_mass < s% job% pre_ms_check_radiative_core_min_mass) then
            write(*,*) 'stop relax to begin radiative core because star_mass < pre_ms_check_radiative_core_min_mass'
            return
         end if
         
         max_timestep = 1d3*secyer  ! can provide a parameter for this if necessary
         
         write(*,*)
         write(*,1) 'relax_to_radiative_core'
         write(*,*)
         if (max_timestep <= 0) then
            rpar(1) = secyer
         else
            rpar(1) = max_timestep
         end if
         rpar(2) = 1d99 ! min_conv_mx1_bot
         max_model_number = s% max_model_number
         model_number = s% model_number
         save_max_timestep = s% max_timestep
         save_max_years_for_timestep = s% max_years_for_timestep
         s% model_number = 0
         call do_internal_evolve( &
               id, before_evolve_relax_to_radiative_core, &
               relax_to_radiative_core_adjust_model, relax_to_radiative_core_check_model, &
               null_finish_model, .true., lipar, ipar, lrpar, rpar, ierr)
         s% max_model_number = max_model_number
         s% model_number = model_number
         s% max_timestep = save_max_timestep
         s% max_years_for_timestep = save_max_years_for_timestep
         call error_check('relax_to_radiative_core',ierr)

      end subroutine do_relax_to_radiative_core

      subroutine before_evolve_relax_to_radiative_core(s, id, lipar, ipar, lrpar, rpar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         ierr = 0
         call setup_before_relax(s)
         !s% max_timestep = rpar(1)
         !s% max_years_for_timestep = s% max_timestep/secyer
         !s% dt_next = s% max_timestep
      end subroutine before_evolve_relax_to_radiative_core

      integer function relax_to_radiative_core_adjust_model(s, id, lipar, ipar, lrpar, rpar)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         relax_to_radiative_core_adjust_model = keep_going
      end function relax_to_radiative_core_adjust_model

      integer function relax_to_radiative_core_check_model(s, id, lipar, ipar, lrpar, rpar)
         use do_one_utils, only:do_bare_bones_check_model
         use report, only: set_power_info
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         logical, parameter :: dbg = .false.
         integer :: ierr
         real(dp) :: min_conv_mx1_bot
         include 'formats'
         relax_to_radiative_core_check_model = do_bare_bones_check_model(id)
         if (relax_to_radiative_core_check_model /= keep_going) return
         ierr = 0
         call set_power_info(s)
         if (s% L_nuc_burn_total/s% L_phot > s% job% pre_ms_check_radiative_core_Lnuc_div_L_limit) then
            write(*,*) 'reached pre_ms_check_radiative_core_Lnuc_div_L_limit in relax to begin radiative core'
            relax_to_radiative_core_check_model = terminate
            s% termination_code = t_relax_finished_okay
            return
         end if
         min_conv_mx1_bot = rpar(2)
         if (s% conv_mx1_bot < min_conv_mx1_bot) then
            min_conv_mx1_bot = s% conv_mx1_bot
            rpar(2) = min_conv_mx1_bot
         end if
         if (min_conv_mx1_bot < s% job% pre_ms_check_radiative_core_start) then
            if (s% conv_mx1_bot > s% job% pre_ms_check_radiative_core_stop) then
               write(*,2) 'finished relax to begin radiative core', s% model_number, s% conv_mx1_bot
               relax_to_radiative_core_check_model = terminate
               s% termination_code = t_relax_finished_okay
            else if (mod(s% model_number, s% terminal_interval) == 0) then
               write(*,2) 'relative mass of radiative core still tiny', s% model_number, s% conv_mx1_bot
            end if
         else if (mod(s% model_number, s% terminal_interval) == 0) then
            write(*,2) 'waiting for fully convective core to develop', s% model_number, s% conv_mx1_bot
         end if
      end function relax_to_radiative_core_check_model


      subroutine do_relax_Z(id, new_z, dlnz, minq, maxq, ierr)
         use star_utils, only: eval_current_z
         use adjust_xyz, only:set_z
         integer, intent(in) :: id
         real(dp), intent(in) :: new_z, dlnz, minq, maxq
         integer, intent(out) :: ierr
         integer, parameter ::  lipar=1, lrpar=5
         integer :: max_model_number
         real(dp) :: z
         type (star_info), pointer :: s
         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:)
         real(dp), target :: rpar_ary(lrpar)
         real(dp), pointer :: rpar(:)
         rpar => rpar_ary
         ipar => ipar_ary
         include 'formats'
         ierr = 0
         if (new_Z < 0 .or. new_Z > 1) then
            ierr = -1
            write(*,*) 'invalid new_Z', new_Z
            return
         end if
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         z = eval_current_z(s, 1, s% nz, ierr)
         if (ierr /= 0) return
         if (abs(new_z - z) <= 1d-6*z) return
         if (max(new_z, z) > 1d-6) then
            if (abs(new_z - z) <= 1d-3*new_z) then
               call set_z(s, new_z, 1, s% nz, ierr)
               return
            end if
         end if
         write(*,*)
         write(*,1) 'current Z', z
         write(*,1) 'relax to new Z', new_z
         write(*,1) '(new - current) / current', (new_z - z) / z
         write(*,*)
         write(*,1) 'dlnz per step', dlnz
         write(*,*)
         rpar(1) = log(max(min_z,new_z))
         rpar(2) = new_z
         rpar(3) = abs(dlnz)
         rpar(4) = minq
         rpar(5) = maxq
         max_model_number = s% max_model_number
         s% max_model_number = -1111
         s% initial_z = z
         call do_internal_evolve( &
               id, before_evolve_relax_Z, &
               relax_Z_adjust_model, relax_Z_check_model, &
               null_finish_model, .true., lipar, ipar, lrpar, rpar, ierr)
         s% max_model_number = max_model_number
         call error_check('relax Z',ierr)
      end subroutine do_relax_Z


      subroutine before_evolve_relax_Z(s, id, lipar, ipar, lrpar, rpar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         ierr = 0
         call setup_before_relax(s)
         s% max_model_number = -111
         s% max_timestep = secyer
         s% dt_next = s% max_timestep
      end subroutine before_evolve_relax_Z

      integer function relax_Z_adjust_model(s, id, lipar, ipar, lrpar, rpar)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         relax_Z_adjust_model = keep_going
      end function relax_Z_adjust_model

      integer function relax_Z_check_model(s, id, lipar, ipar, lrpar, rpar)
         use adjust_xyz, only: set_z
         use star_utils, only: k_for_q, eval_current_z
         use do_one_utils, only:do_bare_bones_check_model
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer :: ierr, klo, khi
         real(dp) :: lnz_target, new_z, new_lnz, dlnz, lnz, current_z, next_z, &
            min_q_for_relax_Z, max_q_for_relax_Z

         logical, parameter :: zdbg = .true.

         include 'formats'

         relax_Z_check_model = do_bare_bones_check_model(id)
         if (relax_Z_check_model /= keep_going) return

         lnz_target = rpar(1)
         new_z = rpar(2)
         dlnz = rpar(3)
         min_q_for_relax_Z = rpar(4)
         max_q_for_relax_Z = rpar(5)

         khi = k_for_q(s, min_q_for_relax_Z)
         klo = k_for_q(s, max_q_for_relax_Z)
         if (zdbg) write(*,2) 'klo', klo, max_q_for_relax_Z
         if (zdbg) write(*,2) 'khi', khi, min_q_for_relax_Z

         current_z = eval_current_z(s, klo, khi, ierr)
         if (ierr /= 0) return

         if (mod(s% model_number, s% terminal_interval) == 0) &
            write(*,1) 'new_z, current', new_z, current_z

         if (abs(current_z-new_z) <= 1d-6*new_z) then
            relax_Z_check_model = terminate
            s% termination_code = t_relax_finished_okay
            return
         end if

         lnz = log(max(min_z,current_z))

         if (zdbg) then
            write(*,1) 'lnz_target', lnz_target
            write(*,1) 'lnz', lnz
            write(*,1) 'lnz - lnz_target', lnz - lnz_target
            write(*,1) 'dlnz', dlnz
         end if

         if (abs(lnz - lnz_target) < dlnz) then
            dlnz = abs(lnz - lnz_target)
            if (zdbg) write(*,1) 'reduced dlnz', dlnz
         end if

         if (lnz_target < lnz) then
            new_lnz = lnz - dlnz
         else
            new_lnz = lnz + dlnz
         end if

         if (zdbg) write(*,1) 'new_lnz', new_lnz

         if (new_lnz >= min_dlnz) then
            next_z = exp(new_lnz)
         else
            next_z = new_z
         end if

         if (zdbg) write(*,1) 'next_z', next_z

         ierr = 0
         call set_z(s, next_z, klo, khi, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*, *) 'relax_Z_check_model ierr', ierr
            relax_Z_check_model = terminate
            s% result_reason = nonzero_ierr
            return
         end if

         write(*,1) 'relax Z, z diff, new, current', new_z - current_z, new_z, current_z

         if (klo == 1 .and. khi == s% nz) s% initial_z = next_z
         s% max_timestep = secyer*s% time_step

      end function relax_Z_check_model




      subroutine do_relax_Y(id, new_Y, dY, minq, maxq, ierr)
         use star_utils, only: k_for_q, eval_current_y
         integer, intent(in) :: id
         real(dp), intent(in) :: new_Y, dY, minq, maxq
         integer, intent(out) :: ierr
         integer, parameter ::  lipar=1, lrpar=4
         integer :: max_model_number, khi, klo
         real(dp) :: y
         type (star_info), pointer :: s
         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:)
         real(dp), target :: rpar_ary(lrpar)
         real(dp), pointer :: rpar(:)
         rpar => rpar_ary
         ipar => ipar_ary
         include 'formats'
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         khi = k_for_q(s, minq)
         klo = k_for_q(s, maxq)
         
         y = eval_current_y(s, klo, khi, ierr)
         if (ierr /= 0) return
         if (is_bad(y)) then
            write(*,1) 'y', y
            stop 'do_relax_Y'
         end if
         
         if (abs(new_Y - y) <= 1d-6*new_Y) return
         if (new_Y < 0 .or. new_Y > 1) then
            ierr = -1
            write(*,*) 'invalid new_Y', new_Y
            return
         end if
         write(*,*)
         write(*,1) 'current Y', Y
         write(*,1) 'relax to new_Y', new_Y
         write(*,1) 'dY per step', dY
         write(*,*)
         rpar(1) = new_Y
         rpar(2) = abs(dY)
         rpar(3) = minq
         rpar(4) = maxq
         max_model_number = s% max_model_number
         s% max_model_number = -1111
         s% initial_y = y
         call do_internal_evolve( &
               id, before_evolve_relax_Y, &
               relax_Y_adjust_model, relax_Y_check_model, &
               null_finish_model, .true., lipar, ipar, lrpar, rpar, ierr)
         s% max_model_number = max_model_number
         call error_check('relax Y',ierr)
      end subroutine do_relax_Y


      subroutine before_evolve_relax_Y(s, id, lipar, ipar, lrpar, rpar, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         ierr = 0
         call setup_before_relax(s)
         s% max_model_number = -111
         s% max_timestep = secyer
         s% dt_next = s% max_timestep
      end subroutine before_evolve_relax_Y

      integer function relax_Y_adjust_model(s, id, lipar, ipar, lrpar, rpar)
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         relax_Y_adjust_model = keep_going
      end function relax_Y_adjust_model

      integer function relax_Y_check_model(s, id, lipar, ipar, lrpar, rpar)
         use star_utils, only: k_for_q, eval_current_y
         use adjust_xyz, only: set_y
         use do_one_utils, only: do_bare_bones_check_model
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer :: ierr, klo, khi
         real(dp) :: new_y, dy, current_y, next_y, minq, maxq, actual_next_y
         logical, parameter :: ydbg = .true.
         logical :: end_now

         include 'formats'

         end_now=.false.

         relax_Y_check_model = do_bare_bones_check_model(id)
         if (relax_Y_check_model /= keep_going) return

         new_y = rpar(1)
         dy = rpar(2)
         minq = rpar(3)
         maxq = rpar(4)
         
         khi = k_for_q(s, minq)
         klo = k_for_q(s, maxq)

         if (ydbg) then
            write(*,4) 'klo, khi nz', klo, khi, s% nz
         end if
         
         current_y = eval_current_y(s, klo, khi, ierr)
         if (is_bad(current_y)) then
            write(*,1) 'current_y', current_y
            stop 'relax_y_check_model'
         end if
         if (ierr /= 0) return

         if (mod(s% model_number, s% terminal_interval) == 0) then
            write(*,1) 'new_y', new_y
            write(*,1) 'dy', dy
            write(*,1) 'current_y', current_y
            write(*,1) 'current_y - new_y', current_y - new_y
         end if

         if (abs(current_y - new_y) < 1d-10) then
            relax_Y_check_model = terminate
            s% termination_code = t_relax_finished_okay
            return
         end if

         if (abs(current_y - new_y) < dY) then
            dY = abs(current_y - new_y)
            end_now = .true.
            if (ydbg) write(*,1) 'reduced dY', dY
         end if

         if (new_y >= current_y) then
            next_y = current_y + dy
         else
            next_y = current_y - dy
         end if

         if (ydbg) write(*,1) 'next_y', next_y

         ierr = 0

         call set_y(s, next_y, klo, khi, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*, *) 'relax_Y_check_model ierr', ierr
            relax_Y_check_model = terminate
            s% result_reason = nonzero_ierr
            return
         end if

         actual_next_y = eval_current_y(s, klo, khi, ierr)

         write(*,1) 'relax Y, y diff, new, current', new_y - current_y, new_y, actual_next_y

         if (ydbg) write(*,1) 'actual_next_y', actual_next_y
         if (ydbg) write(*,1) 'actual_next_y - next_y', actual_next_y - next_y
         if (ydbg) write(*,1) 'y diff', actual_next_y - new_y

         if (abs(actual_next_y - new_y) < 1d-10 .or. end_now) then
            relax_Y_check_model = terminate
            s% termination_code = t_relax_finished_okay
            return
         end if

         if (klo == 1 .and. khi == s% nz) s% initial_y = next_y
         s% max_timestep = secyer*s% time_step

      end function relax_Y_check_model


      subroutine setup_before_relax(s)
         type (star_info), pointer :: s
         s% dxdt_nuc_factor = 0
         s% max_age = 1d50
         s% max_age_in_seconds = 1d50
         s% max_timestep_factor = 2
         s% max_model_number = -1111
         call turn_off_winds(s)
      end subroutine setup_before_relax


      subroutine turn_off_winds(s)
         type (star_info), pointer :: s
         s% mass_change = 0
         s% Reimers_scaling_factor = 0d0
         s% Blocker_scaling_factor = 0d0
         s% de_Jager_scaling_factor = 0d0
         s% van_Loon_scaling_factor = 0d0
         s% Nieuwenhuijzen_scaling_factor = 0d0
         s% Vink_scaling_factor = 0d0
         s% Dutch_scaling_factor = 0d0
         s% use_other_wind = .false.
      end subroutine turn_off_winds


      subroutine do_internal_evolve( &
            id, before_evolve, adjust_model, check_model,&
            finish_model, restore_at_end, &
            lipar, ipar, lrpar, rpar, ierr)
         use evolve
         use star_utils, only: yrs_for_init_timestep
         use pgstar
         use history_specs, only: set_history_columns
         use profile, only: set_profile_columns
         integer, intent(in) :: id, lipar, lrpar
         logical, intent(in) :: restore_at_end
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         interface
            subroutine before_evolve(s, id, lipar, ipar, lrpar, rpar, ierr)
               use const_def, only: dp
               use star_private_def, only:star_info
               type (star_info), pointer :: s
               integer, intent(in) :: id, lipar, lrpar
               integer, intent(inout), pointer :: ipar(:) ! (lipar)
               real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
               integer, intent(out) :: ierr
            end subroutine before_evolve
            integer function adjust_model(s, id, lipar, ipar, lrpar, rpar)
               use const_def, only: dp
               use star_private_def, only:star_info
               type (star_info), pointer :: s
               integer, intent(in) :: id, lipar, lrpar
               integer, intent(inout), pointer :: ipar(:) ! (lipar)
               real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
            end function adjust_model
            integer function check_model(s, id, lipar, ipar, lrpar, rpar)
               use const_def, only: dp
               use star_private_def, only:star_info
               type (star_info), pointer :: s
               integer, intent(in) :: id, lipar, lrpar
               integer, intent(inout), pointer :: ipar(:) ! (lipar)
               real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
            end function check_model
            integer function finish_model(s)
               use star_def, only:star_info
               type (star_info), pointer :: s
            end function finish_model
         end interface
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: result, model_number, model_number_for_last_retry, &
            recent_log_header, num_retries, &
            photo_interval, profile_interval, priority_profile_interval, &
            model_number_old, max_number_retries, &
            solver_iters_timestep_limit, iter_for_resid_tol2, iter_for_resid_tol3, &
            steps_before_use_gold_tolerances, steps_before_use_gold2_tolerances
         real(dp) :: star_age, time, max_age, max_age_in_seconds, max_timestep, &
            Reimers_scaling_factor, Blocker_scaling_factor, de_Jager_scaling_factor, Dutch_scaling_factor, &
            van_Loon_scaling_factor, Nieuwenhuijzen_scaling_factor, Vink_scaling_factor, &
            dxdt_nuc_factor, tol_correction_norm, tol_max_correction, warning_limit_for_max_residual, &
            tol_residual_norm1, tol_max_residual1, &
            tol_residual_norm2, tol_max_residual2, &
            tol_residual_norm3, tol_max_residual3, &
            max_timestep_factor, mass_change, varcontrol_target, dt_next, &
            time_old, maxT_for_gold_tolerances
         logical :: do_history_file, write_profiles_flag, first_try, use_other_wind, &
            use_gold_tolerances, use_gold2_tolerances

         procedure(integer), pointer :: tmp_ptr1 => null(), tmp_ptr3 => null()
         procedure(), pointer :: tmp_ptr2 => null(), tmp_ptr4 => null()
         
         include 'formats'

         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         call save_stuff

         s% do_history_file = .false.
         s% write_profiles_flag = .false.
         s% warning_limit_for_max_residual = 1d99
         s% recent_log_header = -1
         s% max_number_retries = s% relax_max_number_retries
         s% use_gold_tolerances = s% relax_use_gold_tolerances
         s% steps_before_use_gold_tolerances = -1
         s% use_gold2_tolerances = .false.
         s% steps_before_use_gold2_tolerances = -1
                  
         if (s% relax_solver_iters_timestep_limit /= 0) &
            s% solver_iters_timestep_limit = s% relax_solver_iters_timestep_limit

         if (s% relax_tol_correction_norm /= 0) &
            s% tol_correction_norm = s% relax_tol_correction_norm
         if (s% relax_tol_max_correction /= 0) &
            s% tol_max_correction = s% relax_tol_max_correction
            
         if (s% relax_iter_for_resid_tol2 /= 0) &
            s% iter_for_resid_tol2 = s% relax_iter_for_resid_tol2
         if (s% relax_tol_residual_norm1 /= 0) &
            s% tol_residual_norm1 = s% relax_tol_residual_norm1
         if (s% relax_tol_max_residual1 /= 0) &
            s% tol_max_residual1 = s% relax_tol_max_residual1
         
         if (s% relax_iter_for_resid_tol3 /= 0) &
            s% iter_for_resid_tol3 = s% relax_iter_for_resid_tol3
         if (s% relax_tol_residual_norm2 /= 0) &
            s% tol_residual_norm2 = s% relax_tol_residual_norm2
         if (s% relax_tol_max_residual2 /= 0) &
            s% tol_max_residual2 = s% relax_tol_max_residual2
            
         if (s% relax_tol_residual_norm3 /= 0) &
            s% tol_residual_norm3 = s% relax_tol_residual_norm3
         if (s% relax_tol_max_residual3 /= 0) &
            s% tol_max_residual3 = s% relax_tol_max_residual3
            
         if (s% relax_maxT_for_gold_tolerances /= 0) &
            s% maxT_for_gold_tolerances = s% relax_maxT_for_gold_tolerances

         if (s% doing_first_model_of_run) then
            s% num_retries = 0
            s% time = 0
            s% star_age = 0
            s% model_number_for_last_retry = 0
            s% photo_interval = 0
            s% profile_interval = 0
            s% priority_profile_interval = 0
         end if

         if( s% job% pgstar_flag) then
            ! Can't use the star_lib versions otherwise we have a circular dependency in the makefile
            call set_history_columns(id,s% job% history_columns_file, .true., ierr)
            if (ierr /= 0) return
            call set_profile_columns(id, s% job% profile_columns_file, .true., ierr)
            if (ierr /= 0) return
         end if

         if (s% doing_first_model_of_run .and. s% job% pgstar_flag) then
            call do_start_new_run_for_pgstar(s, ierr)
            if (ierr /= 0) return
         end if

         call before_evolve(s, id, lipar, ipar, lrpar, rpar, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in before_evolve'
            return
         end if

         s% termination_code = -1
         s% doing_relax = .true.
         s% need_to_setvars = .true.

         evolve_loop: do ! evolve one step per loop

            first_try = .true.

            step_loop: do ! may need to repeat this loop for retry
            
               result = do_evolve_step_part1(id, first_try)
               if (result == keep_going) &
                  result = adjust_model(s, id, lipar, ipar, lrpar, rpar)
               if (result == keep_going) &
                  result = do_evolve_step_part2(id, first_try)

               if (result == keep_going) result = check_model(s, id, lipar, ipar, lrpar, rpar)
               if (result == keep_going) result = pick_next_timestep(id)
               if (result == keep_going) exit step_loop

               if (result == retry) then
               else if (s% result_reason /= result_reason_normal) then
                  write(*, *) model_number, 'terminate reason: ' // trim(result_reason_str(s% result_reason))
               end if

               if (result == redo) result = prepare_to_redo(id)
               if (result == retry) result = prepare_to_retry(id)
               if (result == terminate) exit evolve_loop
               first_try = .false.

            end do step_loop

            if (.false. .and. s% job% pgstar_flag) then
               ! Can't use the star_lib versions otherwise we have a circular dependency in the makefile
               write(*,2) 'after step_loop: call update_pgstar_data', s% model_number
               call update_pgstar_data(s, ierr)
               if (failed()) return
               call do_read_pgstar_controls(s, s% inlist_fname, ierr) 
               if (failed()) return
               call do_pgstar_plots( s, .false., ierr)
               if (failed()) return
            end if
            
            result = finish_model(s)
            if (result /= keep_going) exit evolve_loop

            tmp_ptr1 => s% how_many_extra_history_columns
            tmp_ptr2 => s% data_for_extra_history_columns
            tmp_ptr3 => s% how_many_extra_profile_columns
            tmp_ptr4 => s% data_for_extra_profile_columns

            s% how_many_extra_profile_columns => no_extra_profile_columns
            s% data_for_extra_profile_columns => none_for_extra_profile_columns
            s% how_many_extra_history_columns => no_extra_history_columns
            s% data_for_extra_history_columns => none_for_extra_history_columns

            result = finish_step(id, .false., ierr)
            s% how_many_extra_history_columns => tmp_ptr1
            s% data_for_extra_history_columns => tmp_ptr2
            s% how_many_extra_profile_columns => tmp_ptr3
            s% data_for_extra_profile_columns => tmp_ptr4
            nullify(tmp_ptr1,tmp_ptr2,tmp_ptr3,tmp_ptr4)

            if (result /= keep_going) exit evolve_loop

            if (associated(s% finish_relax_step)) then
               result = s% finish_relax_step(id)
               if (result /= keep_going) exit evolve_loop
            end if

         end do evolve_loop

         if (s% job% pgstar_flag) then
         ! Can't use the star_lib versions otherwise we have a circular dependency in the makefile
         write(*,2) 'after evolve_loop: call update_pgstar_data', s% model_number
            call update_pgstar_data(s, ierr)
            if (ierr /= 0) return
            call do_read_pgstar_controls(s, s% inlist_fname, ierr)
            if (ierr /= 0) return
            call do_pgstar_plots(s, s% job% save_pgstar_files_when_terminate, ierr)
            if (ierr /= 0) return
         end if

         s% doing_relax = .false.
         s% need_to_setvars = .true. ! just to be safe

         if (.not. (s% termination_code == t_relax_finished_okay .or. &
                    s% termination_code == t_extras_check_model)) ierr = -1

         s% termination_code = -1

         if (associated(s% finished_relax)) call s% finished_relax(id)

         if (restore_at_end) call restore_stuff
         
         if (s% job% set_cumulative_energy_error_each_relax) &
            s% cumulative_energy_error = s% job% new_cumulative_energy_error

         s% dt = 0
         s% dt_old = 0

         s% timestep_hold = -100
         s% model_number_for_last_retry = -100

         contains
         
         logical function failed()
            failed = .false.
            if (ierr == 0) return
            failed = .true.
         end function failed

         subroutine save_stuff

            warning_limit_for_max_residual = s% warning_limit_for_max_residual
            do_history_file = s% do_history_file
            write_profiles_flag = s% write_profiles_flag
            recent_log_header = s% recent_log_header
            mass_change = s% mass_change
            Reimers_scaling_factor = s% Reimers_scaling_factor
            Blocker_scaling_factor = s% Blocker_scaling_factor
            de_Jager_scaling_factor = s% de_Jager_scaling_factor
            van_Loon_scaling_factor = s% van_Loon_scaling_factor
            Nieuwenhuijzen_scaling_factor = s% Nieuwenhuijzen_scaling_factor
            Vink_scaling_factor = s% Vink_scaling_factor
            Dutch_scaling_factor = s% Dutch_scaling_factor
            use_other_wind = s% use_other_wind

            num_retries = s% num_retries
            star_age = s% star_age
            time = s% time
            model_number = s% model_number
            dxdt_nuc_factor = s% dxdt_nuc_factor
            max_age = s% max_age
            max_age_in_seconds = s% max_age_in_seconds
            max_timestep_factor = s% max_timestep_factor
            varcontrol_target = s% varcontrol_target
            max_timestep = s% max_timestep
            model_number_for_last_retry = s% model_number_for_last_retry
            photo_interval = s% photo_interval
            profile_interval = s% profile_interval
            priority_profile_interval = s% priority_profile_interval
            dt_next = s% dt_next
            max_number_retries = s% max_number_retries
            
            use_gold2_tolerances = s% use_gold2_tolerances
            steps_before_use_gold2_tolerances = s% steps_before_use_gold2_tolerances
            use_gold_tolerances = s% use_gold_tolerances
            steps_before_use_gold_tolerances = s% steps_before_use_gold_tolerances
            solver_iters_timestep_limit = s% solver_iters_timestep_limit
            tol_correction_norm = s% tol_correction_norm
            tol_max_correction = s% tol_max_correction            
            iter_for_resid_tol2 = s% iter_for_resid_tol2
            tol_residual_norm1 = s% tol_residual_norm1
            tol_max_residual1 = s% tol_max_residual1         
            iter_for_resid_tol3 = s% iter_for_resid_tol3
            tol_residual_norm2 = s% tol_residual_norm2
            tol_max_residual2 = s% tol_max_residual2         
            tol_residual_norm3 = s% tol_residual_norm3
            tol_max_residual3 = s% tol_max_residual3
            maxT_for_gold_tolerances = s% maxT_for_gold_tolerances

            ! selected history
            time_old = s% time_old
            model_number_old = s% model_number_old

         end subroutine save_stuff

         subroutine restore_stuff
            s% warning_limit_for_max_residual = warning_limit_for_max_residual
            s% do_history_file = do_history_file
            s% write_profiles_flag = write_profiles_flag
            s% recent_log_header = recent_log_header
            s% mass_change = mass_change
            s% Reimers_scaling_factor = Reimers_scaling_factor
            s% Blocker_scaling_factor = Blocker_scaling_factor
            s% de_Jager_scaling_factor = de_Jager_scaling_factor
            s% van_Loon_scaling_factor = van_Loon_scaling_factor
            s% Nieuwenhuijzen_scaling_factor = Nieuwenhuijzen_scaling_factor
            s% Vink_scaling_factor = Vink_scaling_factor
            s% Dutch_scaling_factor = Dutch_scaling_factor
            s% use_other_wind = use_other_wind
            s% num_retries = num_retries
            s% star_age = star_age
            s% time = time
            s% model_number = model_number
            s% dxdt_nuc_factor = dxdt_nuc_factor
            s% max_age = max_age
            s% max_age_in_seconds = max_age_in_seconds
            s% max_timestep_factor = max_timestep_factor
            s% varcontrol_target = varcontrol_target
            s% max_timestep = max_timestep
            s% model_number_for_last_retry = model_number_for_last_retry
            s% photo_interval = photo_interval
            s% profile_interval = profile_interval
            s% priority_profile_interval = priority_profile_interval
            s% dt_next = dt_next
            s% max_number_retries = max_number_retries
            
            s% use_gold2_tolerances = use_gold2_tolerances
            s% steps_before_use_gold2_tolerances = steps_before_use_gold2_tolerances
            s% use_gold_tolerances = use_gold_tolerances
            s% steps_before_use_gold_tolerances = steps_before_use_gold_tolerances
            s% solver_iters_timestep_limit = solver_iters_timestep_limit
            s% tol_correction_norm = tol_correction_norm
            s% tol_max_correction = tol_max_correction
            s% iter_for_resid_tol2 = iter_for_resid_tol2
            s% tol_residual_norm1 = tol_residual_norm1
            s% tol_max_residual1 = tol_max_residual1
            s% iter_for_resid_tol3 = iter_for_resid_tol3
            s% tol_residual_norm2 = tol_residual_norm2
            s% tol_max_residual2 = tol_max_residual2
            s% tol_residual_norm3 = tol_residual_norm3
            s% tol_max_residual3 = tol_max_residual3
            s% maxT_for_gold_tolerances = maxT_for_gold_tolerances

            ! selected history
            s% time_old = time_old
            s% model_number_old = model_number_old

         end subroutine restore_stuff

      end subroutine do_internal_evolve


      integer function null_finish_model(s)
         use star_def, only:star_info
         type (star_info), pointer :: s
         null_finish_model = keep_going
      end function null_finish_model


      integer function no_extra_history_columns(id)
         integer, intent(in) :: id
         no_extra_history_columns = 0
      end function no_extra_history_columns


      subroutine none_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine none_for_extra_history_columns


      integer function no_extra_profile_columns(id)
         integer, intent(in) :: id
         no_extra_profile_columns = 0
      end function no_extra_profile_columns


      subroutine none_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine none_for_extra_profile_columns

      subroutine error_check(name,ierr)
         character(len=*), intent(in) :: name
         integer, intent(in) :: ierr
         include 'formats'

         if (ierr /= 0) then
            write(*,*) 'failed in ', name
            return
         end if

         write(*,*)
         write(*,*) 'finished doing ', name
         write(*,*)
      end subroutine error_check


      end module relax


