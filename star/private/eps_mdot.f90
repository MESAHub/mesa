



      module eps_mdot

      use star_private_def
      use const_def
      use star_utils
      use accurate_sum  ! Provides the accurate_real type, which enables us to do
                        !sums and differences without much loss of precision.
      use mass_utils    ! Helper methods for accurately computiong mass differences
                        ! and intersections between different meshes.

      implicit none

      private
      public :: calculate_eps_mdot



      contains


      !> We choose as a convention to fix F on the faces between cells and to be
      !! positive when mass is flowing downward. With this, the mass flux may be
      !! related to the change in the mass of the cell as
      !!
      !! Delta(dm_j) = F_{j} - F_{j+1}.
      !!
      !! The left-hand side is known: we have
      !!
      !! Delta(dm_j) = (dm_j_after_adust_mass - dm_j_before_adjust_mass)
      !!
      !! This is what we call change_in_dm.
      !!
      !! Furthermore we fix the flux through the centre (F_{nz+1}) to zero.
      !! Hence
      !!
      !! F_{nz} = Delta(dm_nz)
      !! F_{nz-1} = F_{nz} + Delta(dm_{nz-1})
      !!
      !! and so on.
      !!
      !! With a little rearranging we therefore find_mass_flux
      !!
      !! F_{j} = F_{j+1} + (dm_j_after_adust_mass - dm_j_before_adjust_mass)
      !!
      !! @param change_in_dm Change in cell masses.
      !! @param nz Number of cells.
      !! @param mass_flux Mass flux between cells.
      subroutine find_mass_flux(nz, change_in_dm, mass_flux)
         ! Inputs
         real(qp), dimension(:), intent(in) :: change_in_dm
         integer, intent(in) :: nz

         ! Intermediates
         integer j

         real(qp), dimension(:), intent(out) :: mass_flux

         mass_flux(nz+1) = 0 ! Fixes the inner boundary.

         do j=nz,1,-1
            mass_flux(j) =  mass_flux(j+1) + change_in_dm(j)
         end do

      end subroutine find_mass_flux


      real(dp) function interpolate_onto_faces(vec, dm, nz, j)
         ! Inputs
         real(dp), dimension(:) :: vec
         real(qp), dimension(:) :: dm
         integer nz, j

         ! Intermediates
         real(dp) alpha, beta, dmbar1

         ! Return
         real(dp) interp

         !!! High-level explanation

         ! Returns the value of vec interpolated onto face j.
         ! This works for j running from 1 through length(vec) + 1 == nz + 1,
         ! so that we include the 'bounding' faces.

         !!! Mechanics explanation

         ! The interpolation is done so that the finite difference
         !
         ! (interp_vec(j+1)-interp_vec(j) / dm(j))
         !
         ! Equals the average of the finite differences
         !
         ! (vec(j+1)-vec(j))/(dm-bar(j+1))
         !
         ! and
         !
         ! (vec(j)-vec(j-1))/(dm-bar(j)),
         !
         ! where
         ! 
         ! dm-bar(j) = (1/2)(dm(j-1) + dm(j)).
         !
         ! This is done because these finite differences are what
         ! MESA is using elsewhere, so in order to ensure consistency
         ! we want our interpolated vector to have derivatives consistent
         ! with these.
         ! 
         ! When j == 1 we can't do this because we don't know vec(j-1), so
         ! we take vec(j-1) == vec(j) and return vec(1).
         ! When j == length(vec) + 1 we likewise need to assume vec(j-1) == vec(j)
         ! so we return vec(nz).

         interp = 0
         if (j == 1) interp = vec(1)
         if (j == nz + 1) interp = vec(nz)
         if (j > 1 .and. j < nz + 1) then
            dmbar1 = (dm(j-1) + dm(j)) / 2.0d0
            beta = dm(j-1) / (2.0d0 * dmbar1)
            alpha = 1 - beta
            interp = alpha * vec(j-1) + beta * vec(j)
         end if
         interpolate_onto_faces = interp

      end function interpolate_onto_faces

      subroutine set_leak_frac(nz, L, dm, dt, specific_entropy, grad_r_sub_grad_a, mass_flux, leak_frac, eps_mdot_leak_frac_factor)
         ! Inputs
         real(qp), dimension(:) :: mass_flux, dm
         real(dp), dimension(:) :: L, grad_r_sub_grad_a, specific_entropy, leak_frac
         real(dp) eps_mdot_leak_frac_factor
         real(dp) dt
         integer nz

         ! Intermediates
         integer k, km1
         real(dp) mass_flux_bar

         !!! High-level explanation

         ! This subroutine calculates the maximum fractional change in entropy that thermal
         ! leakage can produce. When this is zero all changes in state must be adiabatic.
         ! When this is one there is enough time for arbitrarily large changes in entropy.

         ! The derivation is as follows.
         ! The perturbation to the luminosity L of a fluid element due to a perturbation
         ! in its specific entropy S follows dL/L ~ d(grad S)/grad S. In radiative zones this
         ! is ~dS/S. In convective zones it's ~ dS / (S * grad_r_sub_grad_ad).
         ! Because grad_r_sub_grad_ad is of order unity in radiative zones we can use this correction everywhere.
         ! The takes to leak out is therefore dm * TdS/dL ~ dm * grad_r_sub_grad_ad * TS/L, where dm is the mass of the element.
         ! This is the thermal time-scale of the fluid element. By comparison, in time dm*dt/mass_flux the
         ! element will no longer be in the same cell, so while in this cell a fraction L*dt/(TS mass_flux)
         ! of this perturbation can leak out. The fraction which leaks out cannot exceed one, and so
         ! we obtain the expression below.

         do k=1,nz
            if (k == 1) then
               mass_flux_bar = mass_flux(1)
            else
               km1 = k-1
               mass_flux_bar = (mass_flux(km1) + mass_flux(k)) / 2.0d0
            end if
            if (abs(grad_r_sub_grad_a(k)) == 0d0 .or. mass_flux_bar == 0d0) then
               leak_frac(k) = 1.0d0
            else
               leak_frac(k) = eps_mdot_leak_frac_factor *&
                   abs(L(k) * dt / (grad_r_sub_grad_a(k) * specific_entropy(k) * mass_flux_bar))
               leak_frac(k) = min(1.0d0, leak_frac(k))
            end if
         end do

      end subroutine set_leak_frac

      subroutine leak_control(nz, mass_flux, dm, mesh_intersects, ranges,&
                              total_mass_through_cell, eps_mdot_per_total_mass,&
                              accumulated, mdot_adiabatic_surface, leak_frac)
         ! Inputs
         integer nz
         integer, dimension(:,:) :: ranges
         real(dp) mdot_adiabatic_surface
         real(qp), dimension(:) :: mass_flux, dm, mesh_intersects, total_mass_through_cell
         real(dp), dimension(:) :: eps_mdot_per_total_mass, accumulated, leak_frac

         ! Intermediates
         integer i, j, k
         integer i0, i1
         integer i_start, i_end
         integer, dimension(:), allocatable :: i_min, i_max, j_min, j_max
         logical do_now
         real(qp) delta_m
         real(dp) sgn
         real(dp), dimension(:), allocatable :: excess
         type(non_rect_array), dimension(:), allocatable :: pf

         !!! High-level explanation

         ! When the thermal time-scale is long, material changes adiabatically as it descends.
         ! When the thermal time-scale is short, material bleeds entropy along the way.
         ! The adjust_mass routine assumes that we are in the latter limit and so produces a state
         ! which assumes entropy has been lost along the way.
         ! If this is correct, we should count up the entropy lost in each cell and use that to set eps_mdot.
         ! If this is not correct, we should instead follow parcels of material as they descend, track
         ! the entropy they were (incorrectly) assumed to have lost, and restore it via eps_mdot.
         ! We connect these limits using the leak fraction, as defined in set_leak_frac.

         ! The calculation of what energy to restore and what to leak is done in the leak subroutine.
         ! This (leak_control) subroutine precomputes some useful quantities for that calculation and
         ! then divides the star into contiguous blocks within which the mass flows monotonically
         ! up (surface-ward) or down (center-ward). It then calls leak on each such block.

         ! Initialize
         accumulated = 0
         mdot_adiabatic_surface = 0
         delta_m = mass_flux(1)
         if (mass_flux(1) == 0) then
            ! Should never happen, because we bail out earlier if mdot == 0.
            ! Nevertheless the logic that follows should work without issue.
            ! We just need to set delta_m to be something non-zero for the purposes
            ! of this subroutine because -delta_m is used to normalize the pass fraction.
            ! That is ALL it is used to do in this routine however, and the normalization
            ! is arbitrary (we just choose to use delta_m for convenience) so we can set it
            ! to anything we like.
            delta_m = -1.
         end if

         ! Calculate pass fraction:
         ! pf holds the fraction of mass in cell j which passed through cell i for
         ! all (i,j) for which this is neither 0 nor 1. The method compute_pass_fraction
         ! takes this as an input and does the relevant tests to determine which of
         ! (0, 1, pf(j)%arr(i)) is the appropriate fraction.
         ! Note that the pass fraction for cell j == 0 is either zero (if delta_m > 0) or
         ! else is normalized against the total mass which leaves the model (-delta_m == -mass_flux(1)).
         allocate(i_min(0:nz), i_max(0:nz), pf(0:nz))
         call prepare_pass_fraction(nz, delta_m, dm, mesh_intersects, ranges, i_min, i_max, pf)


         ! We determine for each cell i the set of cells {j} whose ending material
         ! passed through i during adjust_mass. This is a contiguous set because if
         ! material in cells k1 and k2 (k1 < k2) passed through cell i then because
         ! the material between cells k1 and k2 is between them it must also have passed
         ! through cell i. Hence we just need to find the first and last cell which
         ! contains material which passed through cell i. j_min and j_max. These are
         ! computed in the mass_utils subroutine find_j_ranges.
         allocate(j_min(nz), j_max(nz))
         call find_j_ranges(nz, ranges, mass_flux, j_min, j_max)

         ! Starting state
         i_start = 1 ! Where we'll begin leaking
         i_end = 1 ! Where we'll end leaking
         if (mass_flux(1) == 0) then
            sgn = 0
         else
            sgn = mass_flux(1) / abs(mass_flux(1))
         end if

         ! The following loop cuts the star up into segments in which
         ! mass is either all flowing or all flowing down. It then calls
         ! the leak subroutine on each of these in turn.
         do while (i_start <= nz .and. i_end <= nz)
            if (sgn > 0) then
               if (mass_flux(i_end + 1) > 0) then
                  i_end = i_end + 1
               else if (mass_flux(i_end + 1) < 0) then
                  call leak(nz, i_start, i_end, i_min, i_max, j_min, j_max, pf,&
                      dm, delta_m, accumulated, eps_mdot_per_total_mass, leak_frac,&
                      total_mass_through_cell, mdot_adiabatic_surface)
                  i_start = i_end
                  sgn = -1
               else
                  call leak(nz, i_start, i_end, i_min, i_max, j_min, j_max, pf,&
                      dm, delta_m, accumulated, eps_mdot_per_total_mass, leak_frac,&
                      total_mass_through_cell, mdot_adiabatic_surface)
                  i_start = i_end
                  sgn = 0
               end if
            else if (sgn < 0) then
               if (mass_flux(i_start + 1) < 0) then
                  i_start = i_start + 1
               else if (mass_flux(i_start + 1) > 0) then
                  call leak(nz, i_start, i_end, i_min, i_max, j_min, j_max, pf,&
                      dm, delta_m, accumulated, eps_mdot_per_total_mass, leak_frac,&
                      total_mass_through_cell, mdot_adiabatic_surface)
                  i_end = i_start
                  sgn = 1
               else
                  call leak(nz, i_start, i_end, i_min, i_max, j_min, j_max, pf,&
                      dm, delta_m, accumulated, eps_mdot_per_total_mass, leak_frac,&
                      total_mass_through_cell, mdot_adiabatic_surface)
                  i_end = i_start
                  sgn = 0
               end if
            else
               if (mass_flux(i_start + 1) > 0) then
                  i_start = i_start + 1
                  i_end = i_end + 1
                  sgn = 1
               else if (mass_flux(i_start + 1) < 0) then
                  i_start = i_start + 1
                  i_end = i_end + 1
                  sgn = -1
               else
                  call leak(nz, i_start, i_end, i_min, i_max, j_min, j_max, pf,&
                      dm, delta_m, accumulated, eps_mdot_per_total_mass, leak_frac,&
                      total_mass_through_cell, mdot_adiabatic_surface)
                  i_start = i_start + 1
                  i_end = i_end + 1
               end if
            end if
         end do

         do j=1,nz
            accumulated(j) = accumulated(j) / dm(j)
         end do

         ! This handles the contribution of material which starts and ends in the same cell,
         ! or which starts outside the star and ends in the top cell.
         do j=1,nz
            if (leak_frac(j) == 1) then
               ! Means all energy leaks into this cell, so we don't do the normal leak calculation.
               accumulated(j) = accumulated(j) + eps_mdot_per_total_mass(j) * total_mass_through_cell(j)
            else
               ! Leakage happens so we just need this piece.
               accumulated(j) = accumulated(j) + eps_mdot_per_total_mass(j) * dm(j) * compute_pass_fraction(j, j, i_min, i_max, pf)
            end if
         end do

      end subroutine leak_control


      subroutine leak(nz, i_start, i_end, i_min, i_max, j_min, j_max, pf,&
                      dm, delta_m, accumulated, eps_mdot_per_total_mass, leak_frac,&
                      total_mass_through_cell, mdot_adiabatic_surface)

         ! Inputs
         integer nz, i_start, i_end
         integer, dimension(:) :: i_min, i_max, j_min, j_max
         type(non_rect_array), dimension(:) :: pf
         real(qp) delta_m
         real(dp), dimension(:) :: leak_frac, accumulated, eps_mdot_per_total_mass
         real(qp), dimension(:) :: dm, total_mass_through_cell
         real(dp) mdot_adiabatic_surface


         ! Intermediates
         integer i, j, k, direction, ii
         real(qp) pass_frac, next, pass_mass
         real(dp), dimension(:), allocatable :: excess

         !!! High-level explanation

         ! The leak fraction determines the fraction of any attempted entropy
         ! change which leaks out (versus being held adiabatically).
         ! excess is the amount of heat the material has which it would deposit
         ! if given the time.
         ! Hence if we follow the material that ultimately ends in cell j as it
         ! passes through cell i we find
         !
         ! excess(j) -> excess(j) + (heat picked up while passing through cell i)
         ! (deposited heat) = leak_frac(i) * excess(j)
         ! excess(j) -> (1 - leak_frac) * excess(j)
         !
         ! That is, we first calculate the entropy change associated with moving the material
         ! through the next cell (heat picked up). We then add that to the current excess.
         ! A fraction leak_frac of the cumulative entropy excess is then deposited as heat and
         ! decremented from the excess.
         ! When the material reaches whatever cell it ends in (i == j) the excess is deposited
         ! in that cell. If material exits the star the excess it leaves with is accounted for
         ! in mdot_adiabatic_surface. 



         ! There are only two ways for i_start to equal i_end:
         ! 1. All mass which begins in this cell ends in this cell so there's nothing to do.
         ! 2. This is the top cell and all mass which ends in this cell begins in either this cell or outside the star.
         !    Because no eps_mdot is accumulated outside the star this is equivalent to it all
         !    beginning in the top cell, so there's again nothing to do.
         ! In either case the loop at the end of leak_control will ensure that
         ! accumulated(i_start) = eps_mdot_per_total_mass(i) * dm(i),
         ! which is just the previously calculated eps_mdot(i).
         if (i_start == i_end) return

         ! Determine flow direction
         if (i_start < i_end) then
            direction = 1
         else
            direction = -1
         end if

         allocate(excess(0:nz))
         excess = 0
         i = i_start
         do while (i /= i_end + direction)
            if (leak_frac(i) == 1) then
               ! If this is the first cell with full leakage we dump all excess here.
               ! If not it means we've done this before (so there is no excess)
               ! or we're at the very start (so there is no excess).
               if (i /= i_start) then
                  if (leak_frac(i - direction) < 1) then
                     do j=j_min(i), j_max(i)
                        accumulated(i) = accumulated(i) + excess(j)
                        excess(j) = 0
                     end do
                  end if
               end if
            else
               do j=j_min(i), j_max(i)
                  if ((i == i_end .and. j > 0) .or. j == i) then
                     ! This material ends here, so it drops its heat here.
                     ! The first condition means that this is the last cell (i)
                     ! to be considered in this flow direction and that cell (j)
                     ! does not represent material that exits the star.
                     ! The second condition means the cell under consideration (i)
                     ! is the same as that of the material under consideration (j)
                     ! and so that material must end up here.

                     ! Note that we do not include the eps_mdot contribution from
                     ! this cell here because that could cause double counting in
                     ! cells across which mass_flux changes sign. To ensure single
                     ! counting that contribution is accounted for in the loop
                     ! at the end of leak_frac.
                     accumulated(i) = accumulated(i) + excess(j)
                     excess(j) = 0                    
                  else if (i == i_end .and. i == 1 .and. j == 0) then
                     ! Material with j == 0 exits the star. Note that this implies direction == -1.
                     ! For i > 1 this material can be handled by the 'just passing through' else 
                     ! clause, so we only need to think about the i == 1 case.

                     ! Because this material isn't in the star at the end, we have to account
                     ! for the heat it picks up on its way out:
                     pass_frac = compute_pass_fraction(1, 0, i_min, i_max, pf)
                     excess(j) = excess(j) - delta_m * pass_frac * eps_mdot_per_total_mass(1)

                     ! Because our accounting assumes that material flows through the
                     ! surface of the star with energy equal to the surface specific energy
                     ! we have to account for any excess this material carries on its way out.
                     mdot_adiabatic_surface = mdot_adiabatic_surface + (1 - leak_frac(i_end)) * excess(j)
                     accumulated(i) = accumulated(i) + leak_frac(i_end) * excess(j)

                     excess(j) = 0
                  else
                     ! The material in cell j is just passing through cell i.
                     pass_frac = compute_pass_fraction(i, j, i_min, i_max, pf)
                     if (j > 0) then
                        pass_mass = pass_frac * dm(j)
                     else
                        pass_mass = -pass_frac * delta_m
                     end if
                     next = excess(j) + dm(i) * eps_mdot_per_total_mass(i) * pass_mass

                     accumulated(i) = accumulated(i) + leak_frac(i) * next
                     excess(j) = (1 - leak_frac(i)) * next
                  end if
               end do
            end if
            i = i + direction
         end do

      end subroutine leak

      subroutine calculate_eps_mdot(s, dt, ierr)
         use adjust_mass, only: compute_prev_mesh_dm

         ! Inputs
         type (star_info), pointer :: s
         real(dp) dt
         integer ierr
         
         ! Intermediates
         logical :: dbg = .false.
         integer nz, j, k, l, n
         real(dp) delta_m, sgn, change_sum, leak_sum, err, abs_err, mdot_adiabatic_surface, gradT_mid
         real(dp), dimension(:), allocatable :: &
            p_bar, rho_bar, te_bar, te, curr_m, &
            leak_frac, thermal_energy, density_weighted_flux, eps_mdot_per_total_mass,&
            accumulated, grad_r_sub_grad_a
         real(qp), dimension(:), allocatable :: change_in_dm, mass_flux, dm, prev_mesh_dm,&
             total_mass_through_cell
         type(accurate_real) sum
         integer, dimension(:,:), allocatable :: ranges
         real(qp), dimension(:), allocatable :: remainders, mesh_intersects
         real(qp) m

         if (s% mstar_dot == 0d0 .or. dt <= 0d0) then
            s% eps_mdot(1:s%nz) = 0d0
            s% mdot_adiabatic_surface = 0d0
            s% mdot_acoustic_surface = 0d0
            s% total_internal_energy_after_adjust_mass = 0d0
            s% total_gravitational_energy_after_adjust_mass = 0d0
            s% total_radial_kinetic_energy_after_adjust_mass = 0d0
            s% total_turbulent_energy_after_adjust_mass = 0d0
            s% total_rotational_kinetic_energy_after_adjust_mass = 0d0
            s% total_energy_after_adjust_mass = 0d0
            return
         end if
         
         s% need_to_setvars = .true.

         ! Stellar properties
         nz = s%nz
         delta_m = s%mstar_dot * dt

         allocate(prev_mesh_dm(nz), dm(nz), change_in_dm(nz))
         call compute_prev_mesh_dm(s, prev_mesh_dm, dm, change_in_dm)

         allocate(mass_flux(nz+1))
         call find_mass_flux(nz, change_in_dm, mass_flux)

         ! Tabulate cell intersection widths between the new mesh and the old
         allocate(mesh_intersects(2*nz))               
         allocate(ranges(2*nz,2))
         call make_compressed_intersect(dm, prev_mesh_dm, nz, mesh_intersects, ranges)


         ! Weighted dm/dt by rho, used for calculating PdV work.
         allocate(density_weighted_flux(nz+1))
         density_weighted_flux(nz+1) = 0
         do j=nz,1,-1
            density_weighted_flux(j) = density_weighted_flux(j+1) + change_in_dm(j) / s%rho(j) 
         end do

         ! We attribute eps_mdot evenly to all of the mass which is at any point within a cell.
         ! This is a zeroth-order accurate scheme, but it is good enough for our purposes.
         ! The next best approximation would be to track the amount of "time" material spends
         ! in each cell and assign contributions proportional to that, but unless there is evidence
         ! suggesting this is necessary it seems overly complicated.
         allocate(total_mass_through_cell(nz))
         !$OMP PARALLEL DO
         do j=1,nz
            ! This equals Sum_j pass_fraction(i,j) * dm(j)
            total_mass_through_cell(j) = prev_mesh_dm(j) &
               + max(mass_flux(j),0.0_qp) - min(mass_flux(j+1),0.0_qp)
         end do
         !$OMP END PARALLEL DO

         ! Puts total_energy in specific terms and interpolates that as well as p and rho onto faces.
         ! For convenience we include the bottom face (face nz+1) and the top face (face 1), for which the
         ! interpolation just returns the value in the adjacent cell.
         allocate(te(nz))
         allocate(p_bar(nz+1))
         allocate(rho_bar(nz+1))
         allocate(te_bar(nz+1))

         !$OMP PARALLEL DO
         do j=1,nz
            te(j) = s% total_energy_profile_before_adjust_mass(j) / prev_mesh_dm(j)
         end do
         !$OMP END PARALLEL DO
         call eval_total_energy_profile(s, s% total_energy_profile_after_adjust_mass)


         !$OMP PARALLEL DO
         do j=1,nz+1
            ! We use the previous mesh for interpolation because that's the one for which our derivatives were calculated.
            p_bar(j) = interpolate_onto_faces(s%p, prev_mesh_dm, nz, j) 
            rho_bar(j) = interpolate_onto_faces(s%rho, prev_mesh_dm, nz, j) 
            te_bar(j) = interpolate_onto_faces(te, prev_mesh_dm, nz, j) 
         end do
         !$OMP END PARALLEL DO

         if (s% use_other_accreting_state) then
            call s%other_accreting_state(s%id, te_bar(1), p_bar(1), rho_bar(1), ierr)
            s% surface_cell_specific_total_energy_old = te_bar(1)
         end if

         ! Calculate grad_r_sub_grad_a
         allocate(grad_r_sub_grad_a(nz))
         !$OMP PARALLEL DO PRIVATE(gradT_mid)
         do j=1,nz-1
            gradT_mid = 0.5d0*(s% gradT(j) + s% gradT(j+1))
            grad_r_sub_grad_a(j) = s% grada(j) - gradT_mid
         end do
         !$OMP END PARALLEL DO

         grad_r_sub_grad_a(nz) = grad_r_sub_grad_a(nz-1)

         ! Get thermal energy, which is not quite the same as the specific internal energy
         ! because there may be an offset on that (e.g. due to a high Fermi level).
         ! What we want, more specifically, is a quantity which, when perturbed, perturbs
         ! the luminosity to a comparable degree. S*T is a good approximation of this.
         allocate(thermal_energy(nz))
         !$OMP PARALLEL DO
         do j=1,nz
            thermal_energy(j) = s%T(j) * s%cp(j)
         end do
         !$OMP END PARALLEL DO

         ! Our PV and TE methods can only handle a situation in which all cells increase in mass.
         ! We handle the case in which all cells decrease in mass by swapping dm and prev_mesh_dm.
         ! In effect this let's us treat mass loss as mass accretion in reverse. None of the thermodynamic
         ! routines care about this ordering so we can just negate the output. The leakage routine
         ! does care, because it undergoes an irreversible process, so to make sure that is handled correctly
         ! we explicitly pass the direction of the change to that routine.

         allocate(eps_mdot_per_total_mass(nz))
         do j=1,nz
            eps_mdot_per_total_mass(j) = 0
            s%eps_mdot(j) = 0
         end do

         change_sum = 0

         ! Calculate change in heat.
         s % mdot_acoustic_surface = delta_m * p_bar(1) / rho_bar(1) / dt
         do j=1,nz
            ! We calculate eps_mdot using accurate reals to reduce roundoff errors.
            sum % sum = 0.0d0
            sum % compensator = 0.0d0

            ! More numerically stable version of
            ! te_bar(j) * mass_flux(j) - te_bar(j+1) * mass_flux(j+1) - change_in_dm(j) * te(j)
            sum = sum + real(change_in_dm(j) * te_bar(j),qp)
            sum = sum + real(mass_flux(j+1) * (te_bar(j) - te_bar(j+1)),qp)
            sum = sum - real(change_in_dm(j) * te(j),qp)

            ! (PV F_m - PV_0)_{k} - (PV F_m - PV_0)_{k+1}
            sum = sum + real(p_bar(j) * (mass_flux(j) / rho_bar(j) - density_weighted_flux(j)),qp)
            sum = sum - real(p_bar(j+1) * (mass_flux(j+1) / rho_bar(j+1) - density_weighted_flux(j+1)),qp)

            ! Change in specific total energy in the cell with the final cell mass. The
            ! In the absence of composition changes this is just dm_new * (change in gravitational potential),
            ! but when composition changes this term also captures that.
            sum = sum - real(s% total_energy_profile_after_adjust_mass(j) - s%dm(j) * te(j), qp)

            change_sum = change_sum + sum%value() / (dt)

            ! Multiplicative factors
            eps_mdot_per_total_mass(j) = sum % value() / (s%dm(j) * dt) / total_mass_through_cell(j)

         end do

         err = 0d0
         do j=1,nz
            err = err + eps_mdot_per_total_mass(j) * s%dm(j) * dt * total_mass_through_cell(j)
         end do
         err = err - s%mdot_acoustic_surface
         err = err - te_bar(1) * delta_m

         ! Calculate leak fraction
         allocate(leak_frac(nz))
         call set_leak_frac(nz, s%L, dm, dt, thermal_energy, grad_r_sub_grad_a,&
                            mass_flux, leak_frac, s%eps_mdot_leak_frac_factor)

         ! Now we calculate what ought to have happened accounting for energy retention/leakage.
         allocate(accumulated(nz))
         call leak_control(nz, mass_flux, dm, mesh_intersects, ranges,&
                           total_mass_through_cell, eps_mdot_per_total_mass,&
                           accumulated, mdot_adiabatic_surface, leak_frac)
         do j=1,nz
            s%eps_mdot(j) = accumulated(j)
         end do
         s% mdot_adiabatic_surface = -mdot_adiabatic_surface

         ! Apply user-specified scaling factor.
         ! This non-physical, but can be useful for stellar engineering.
         s% eps_mdot(1:nz) = s% eps_mdot_factor * s% eps_mdot(1:nz)

         ! Evaluate total energies for later use outside of eps_mdot.
            call eval_deltaM_total_energy_integrals( &
               s, 1, nz, s% mstar, .true., &
               s% total_energy_profile_after_adjust_mass, &
               s% total_internal_energy_after_adjust_mass, &
               s% total_gravitational_energy_after_adjust_mass, &
               s% total_radial_kinetic_energy_after_adjust_mass, &
               s% total_rotational_kinetic_energy_after_adjust_mass, &
               s% total_turbulent_energy_after_adjust_mass, &
               s% total_energy_after_adjust_mass)

         if (dbg) then
            write(*,*) 'Mdot',s%mstar_dot
            leak_sum = mdot_adiabatic_surface
            do j=1,nz
               leak_sum = leak_sum + s%eps_mdot(j) * s%dm(j)
            end do
            write(*,*) 'Leak Difference:', (leak_sum - change_sum) / change_sum, dt * change_sum, dt * leak_sum, &
                     dt * (leak_sum - change_sum)

            err = 0
            abs_err = 0
            do j=1,nz
               err = err + s% total_energy_profile_after_adjust_mass(j) - te(j) * prev_mesh_dm(j)
               err = err + s%eps_mdot(j) * dm(j) * dt
!               abs_err = abs_err + abs(s% total_energy_profile_after_adjust_mass(j) - te(j) * prev_mesh_dm(j))
               abs_err = abs_err + abs(s% total_energy_profile_after_adjust_mass(j))
!               write(*,"(5ES16.8)") s% total_energy_profile_after_adjust_mass(j), te(j), prev_mesh_dm(j), &
!               s% total_energy_profile_after_adjust_mass(j) - te(j) * prev_mesh_dm(j), s%eps_mdot(j) * dm(j) * dt
            end do

            err = err - te(1) * delta_m - (s%mdot_acoustic_surface + s%mdot_adiabatic_surface) * dt
            write(*,*) 'Err:', err, err / abs_err, err / delta_m, delta_m / s%m(1), abs_err, s%mdot_acoustic_surface*dt
         end if


      end subroutine calculate_eps_mdot


      end module eps_mdot
