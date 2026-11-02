program test_tdc_mass_remap
   use const_def, only: dp, qp
   use utils_lib, only: mesa_error
   use adjust_mass, only: remap_tdc_turbulent_energy
   use mass_utils, only: make_compressed_intersect
   implicit none

   integer :: i, ierr
   real(qp), parameter :: mass_unit = 1.0e30_qp
   real(dp), parameter :: velocity_unit = 1.0e6_dp
   real(dp), parameter :: tolerance = 5.0e-14_dp
   real(qp) :: old_mass(4), new_mass(4), removed, mass_fraction
   real(dp) :: old_vc(4), new_vc(4), repeated_vc(4), checkpoint_vc(4)
   real(dp) :: reference_vc

   write(*,'(a)') 'case,relative_energy_error,max_relative_velocity_change'

   old_mass = mass_unit*[1.0_qp, 2.0_qp, 4.0_qp, 1.0_qp]
   old_vc = velocity_unit*[0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
   new_mass = old_mass
   call evaluate('identity')
   call require(all(new_vc == old_vc), 'identity changed a starting velocity')
   call require(removed == 0.0_qp, 'identity removed turbulent energy')

   old_vc = 0.0_dp
   new_mass(1:3) = 0.8_qp*old_mass(1:3)
   call evaluate('zero_turbulence_loss')
   call require(all(new_vc == 0.0_dp), 'zero turbulence generated a velocity')
   call require(removed == 0.0_qp, 'zero turbulence removed energy')
   new_mass(1:3) = 1.2_qp*old_mass(1:3)
   call evaluate('zero_turbulence_gain')

   old_vc = velocity_unit
   new_mass(1:3) = 0.8_qp*old_mass(1:3)
   call evaluate('uniform_velocity_loss')
   call require(maxval(abs(new_vc/velocity_unit - 1.0_dp)) < tolerance, &
      'uniform wind changed the retained turbulent velocity')
   call require(abs(real(removed/ &
      (1.5_qp*velocity_unit*velocity_unit*sum(old_mass - new_mass)) - 1.0_qp,dp)) < tolerance, &
      'uniform wind has the wrong escaping turbulent energy')

   old_vc = velocity_unit*[1.0_dp, 2.0_dp, 3.0_dp, 0.0_dp]
   checkpoint_vc = old_vc
   new_mass(1:3) = 0.6_qp*old_mass(1:3)
   call evaluate('nonuniform_velocity_loss')
   call require(all(old_vc == checkpoint_vc), 'remap changed its input checkpoint')
   repeated_vc = new_vc
   call remap_tdc_turbulent_energy(4, old_mass, new_mass, old_vc, new_vc, removed, ierr)
   call require(ierr == 0, 'repeat call returned an error')
   call require(all(new_vc == repeated_vc), 'repeating the same remap changed the result')

   new_mass(1:3) = 1.4_qp*old_mass(1:3)
   call evaluate('nonuniform_velocity_gain')
   call require(removed == 0.0_qp, 'cold accretion removed turbulent energy')

   old_mass = mass_unit
   old_vc = velocity_unit*[1.0_dp, 2.0_dp, 3.0_dp, 0.0_dp]
   mass_fraction = 1.0_qp
   do i=1,3
      mass_fraction = mass_fraction*1.0e-6_qp
      new_mass = old_mass
      new_mass(1:3) = old_mass(1:3)*(1.0_qp + mass_fraction)
      call evaluate('vanishing_mass_gain')
      ! Energy is linear in the moved mass; a zero velocity can grow as sqrt(dM).
      call require(maxval(abs(new_vc*new_vc - old_vc*old_vc))/maxval(old_vc*old_vc) < &
         20.0_dp*real(mass_fraction,dp) + tolerance, 'small mass gain changed the energy profile discontinuously')
      new_mass(1:3) = old_mass(1:3)*(1.0_qp - mass_fraction)
      call evaluate('vanishing_mass_loss')
      call require(maxval(abs(new_vc*new_vc - old_vc*old_vc))/maxval(old_vc*old_vc) < &
         20.0_dp*real(mass_fraction,dp) + tolerance, 'small mass loss changed the energy profile discontinuously')
      call require(abs(real(removed/ &
         (1.5_qp*velocity_unit*velocity_unit*sum(old_mass - new_mass)) - 1.0_qp,dp)) < tolerance, &
         'small mass loss has the wrong escaping turbulent energy')
   end do

   ! Face 3 shares changed cell 2 even though cells 3 and 4 keep their masses.
   ! Its new dual volume contains all of the old face-3 volume plus zero energy.
   old_mass = mass_unit
   new_mass = mass_unit*[1.0_qp, 2.0_qp, 1.0_qp, 1.0_qp]
   old_vc = velocity_unit*[0.0_dp, 0.0_dp, 2.0_dp, 0.0_dp]
   call evaluate('changed_dual_boundary')
   reference_vc = real(sqrt(8.0_qp/3.0_qp),dp)*velocity_unit
   call require(abs(new_vc(3)/reference_vc - 1.0_dp) < tolerance, &
      'the shared boundary face was not transported with its changed dual mass')

   call check_single_cell_accretion
   call check_intersections('single_cell_identity', [1.0_qp], [1.0_qp])
   call check_intersections('single_cell_gain', [1.0_qp], [2.0_qp])
   call check_intersections('single_cell_loss', [2.0_qp], [1.0_qp])
   call check_intersections('multiple_cell_identity', &
      [0.5_qp, 1.0_qp, 1.0_qp, 1.0_qp], [0.5_qp, 1.0_qp, 1.0_qp, 1.0_qp])
   call check_intersections('surface_cell_gain', &
      [0.5_qp, 1.0_qp, 1.0_qp, 1.0_qp], [0.500001_qp, 1.000001_qp, 1.000001_qp, 1.0_qp])
   call check_intersections('multiple_cell_gain', &
      [0.5_qp, 1.0_qp, 1.0_qp, 1.0_qp], [0.7_qp, 1.4_qp, 1.4_qp, 1.0_qp])
   call check_intersections('surface_cell_loss', &
      [0.5_qp, 1.0_qp, 1.0_qp, 1.0_qp], [0.499999_qp, 0.999999_qp, 0.999999_qp, 1.0_qp])
   call check_intersections('multiple_cell_loss', &
      [0.5_qp, 1.0_qp, 1.0_qp, 1.0_qp], [0.3_qp, 0.6_qp, 0.6_qp, 1.0_qp])
   call check_intersections('fixed_total_mass', &
      [0.5_qp, 1.0_qp, 1.0_qp, 1.0_qp], [0.6_qp, 0.9_qp, 1.0_qp, 1.0_qp])

   write(*,'(a)') 'All native turbulent remap checks passed.'

contains

   subroutine evaluate(label)
      character(len=*), intent(in) :: label
      real(qp) :: energy_before, energy_after, relative_error
      real(dp) :: change

      call remap_tdc_turbulent_energy(4, old_mass, new_mass, old_vc, new_vc, removed, ierr)
      call require(ierr == 0, trim(label)//': remap returned an error')
      call require(all(new_vc >= 0.0_dp), trim(label)//': negative turbulent velocity')
      call require(removed >= 0.0_qp, trim(label)//': negative escaping energy')
      energy_before = total_turbulent_energy(old_mass, old_vc)
      energy_after = total_turbulent_energy(new_mass, new_vc)
      relative_error = abs(energy_before - energy_after - removed)/max(1.0_qp, energy_before)
      change = maxval(abs(new_vc - old_vc))/max(1.0_dp, maxval(old_vc))
      write(*,'(a,2(a,es16.8))') trim(label), ',', real(relative_error,dp), ',', change
      call require(relative_error < real(tolerance,qp), trim(label)//': turbulent energy budget did not close')
   end subroutine evaluate


   subroutine check_single_cell_accretion
      real(qp) :: mass_old(1), mass_new(1), escaped
      real(dp) :: vc_old(1), vc_new(1)
      integer :: remap_ierr

      mass_old = 2.0_qp*mass_unit
      mass_new = 4.0_qp*mass_unit
      vc_old = velocity_unit
      call remap_tdc_turbulent_energy(1, mass_old, mass_new, vc_old, vc_new, escaped, remap_ierr)
      call require(remap_ierr == 0, 'single-cell accretion returned an error')
      call require(abs(vc_new(1)/(sqrt(0.5_dp)*velocity_unit) - 1.0_dp) < tolerance, &
         'single-cell cold accretion did not dilute turbulent energy')
      call require(escaped == 0.0_qp, 'single-cell accretion removed turbulent energy')
      call require(abs(real(total_turbulent_energy(mass_new,vc_new)/ &
         total_turbulent_energy(mass_old,vc_old) - 1.0_qp,dp)) < tolerance, &
         'single-cell cold accretion did not conserve turbulent energy')
   end subroutine check_single_cell_accretion


   subroutine check_intersections(label, mass_old, mass_new)
      character(len=*), intent(in) :: label
      real(qp), intent(in) :: mass_old(:), mass_new(:)
      real(qp) :: mesh_intersects(2*size(mass_old)), covered_old(size(mass_old)), &
         covered_new(size(mass_old)), incoming, outgoing, expected_incoming, expected_outgoing, scale
      real(qp), parameter :: intersection_tolerance = 64.0_qp*epsilon(1.0_qp)
      integer :: ranges(2*size(mass_old),2), interval, j, k, nz

      nz = size(mass_old)
      call require(size(mass_new) == nz, trim(label)//': inconsistent test dimensions')
      call make_compressed_intersect(mass_new, mass_old, nz, mesh_intersects, ranges)
      covered_old = 0.0_qp
      covered_new = 0.0_qp
      incoming = 0.0_qp
      outgoing = 0.0_qp
      scale = max(sum(mass_old),sum(mass_new))
      call require(all(mesh_intersects >= -intersection_tolerance*scale), &
         trim(label)//': negative intersection width')
      do interval=1,2*nz
         if (mesh_intersects(interval) <= 0.0_qp) cycle
         j = ranges(interval,1)
         k = ranges(interval,2)
         call require(j >= 0 .and. j <= nz .and. k >= 0 .and. k <= nz, &
            trim(label)//': invalid intersection indices')
         if (j > 0) covered_new(j) = covered_new(j) + mesh_intersects(interval)
         if (k > 0) covered_old(k) = covered_old(k) + mesh_intersects(interval)
         if (j == 0) outgoing = outgoing + mesh_intersects(interval)
         if (k == 0) incoming = incoming + mesh_intersects(interval)
      end do
      call require(all(abs(covered_old - mass_old) <= intersection_tolerance*mass_old), &
         trim(label)//': old-cell intersection widths do not sum to their masses')
      call require(all(abs(covered_new - mass_new) <= intersection_tolerance*mass_new), &
         trim(label)//': new-cell intersection widths do not sum to their masses')
      expected_incoming = max(0.0_qp,sum(mass_new) - sum(mass_old))
      expected_outgoing = max(0.0_qp,sum(mass_old) - sum(mass_new))
      call require(abs(incoming - expected_incoming) <= intersection_tolerance*scale, &
         trim(label)//': incorrect incoming mass')
      call require(abs(outgoing - expected_outgoing) <= intersection_tolerance*scale, &
         trim(label)//': incorrect outgoing mass')
   end subroutine check_intersections


   real(qp) function total_turbulent_energy(mass, velocity) result(energy)
      real(qp), intent(in) :: mass(:)
      real(dp), intent(in) :: velocity(:)
      integer :: k
      real(qp) :: left, right

      energy = 0.0_qp
      do k=1,size(mass)
         left = real(velocity(k),qp)
         right = 0.0_qp
         if (k < size(mass)) right = real(velocity(k+1),qp)
         energy = energy + 0.75_qp*mass(k)*(left*left + right*right)
      end do
   end function total_turbulent_energy


   subroutine require(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (condition) return
      write(*,'(a)') trim(message)
      call mesa_error(__FILE__,__LINE__,'native turbulent remap regression failed')
   end subroutine require

end program test_tdc_mass_remap
