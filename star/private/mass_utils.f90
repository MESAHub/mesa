



      module mass_utils

      use const_def
      use accurate_sum  ! Provides the accurate_real type, which enables us to do
                        !sums and differences without much loss of precision.
 

      implicit none

      private
      public :: accurate_mass_difference, reconstruct_m, reconstruct_xm, &
                make_compressed_intersect, non_rect_array, &
                prepare_pass_fraction, compute_pass_fraction, find_j_ranges, integrate_conserved


      ! A derived array type for storing non-rectangular 2D arrays
      type non_rect_array
         real(qp), dimension(:), allocatable :: arr
      end type non_rect_array


      contains


      ! Reconstructs m(j) given dm.
      ! Not currently used, but helpful for debugging.
      real(qp) function reconstruct_m(dm, nz, j)
         ! Inputs
         real(qp), dimension(:) :: dm
         integer nz, j

         ! Intermediates
         real(qp) sum, compensator
         integer l

         sum = 0.0
         compensator = 0.0
         do l=nz,j,-1
            call neumaier_sum(sum, compensator, dm(l))
         end do

         reconstruct_m = sum + compensator
      end function reconstruct_m


      ! Reconstructs xm(j) given dm.
      ! Not currently used, but helpful for debugging.
      real(qp) function reconstruct_xm(dm, nz, j)
         ! Inputs
         real(qp), dimension(:) :: dm
         integer nz, j

         ! Intermediates
         real(qp) sum, compensator
         integer l

         sum = 0.0
         compensator = 0.0
         do l=1,j
            call neumaier_sum(sum, compensator, dm(l))
         end do

         reconstruct_xm = sum + compensator
      end function reconstruct_xm


      ! Returns
      !
      ! m1(j) - m2(k)
      !
      ! where
      !
      ! m_i(j) = Sum_{l >= j} dm_i(l)
      !
      ! The reason we specify dm's rather than m's is to avoid
      ! having to subtract large numbers and incur the resultant
      ! penalty in precision. 
      !
      ! We require length(dm1) == length(dm2) == nz.
      ! We permit j and k to range from 1 ... nz+1 inclusive.
      ! By the above definition m_i(nz+1) == 0.
      !
      ! With the above we begin by writing
      !
      ! m1(j) - m2(k) = Sum_{l >= j} dm1(l) - Sum_{l >= k} dm2(l)
      !
      ! We then apply Neumaier's algorithm (see accurate_real docs)
      ! to perform this sum accurately. As a preprocessing step we present
      ! terms to this algorithm with alternating signs and
      ! descending magnitudes for as long as possible.
      ! Hence we actually begin with the 
      !
      ! Sum_{l = nz ... max(j,k)} dm1(l) - dm2(l) ()
      !
      ! If j < k we then add Sum_{l=k-1 ... j} dm1(l).
      ! If j > k we then subtract Sum_{l=j-1 ... k} dm2(l). 
      !
      ! Not currently used, but helpful for debugging.
      real(qp) function accurate_mass_difference(dm1, dm2, j, k, nz)
         ! Inputs
         real(qp), dimension(:) :: dm1, dm2
         integer j, k, nz

         ! Intermediates
         type(accurate_real) sum
         real(qp) summand
         integer l

         sum % sum = 0.0
         sum % compensator = 0.0

         if (max(j,k) <= nz) then
            do l=nz,max(j,k),-1
               sum = sum + dm1(l)
               sum = sum - dm2(l)
            end do 
         end if           

         if (j /= k) then
            do l=max(j,k) - 1,min(j,k),-1
               if (j < k) then
                  summand = dm1(l)
               else
                  summand = -dm2(l)
               end if

               sum = sum + summand

            end do
         end if

         accurate_mass_difference = sum % value()
      end function accurate_mass_difference


      ! Returns the size of the overlap between the intervals
      !
      ! [m1(j+1), m1(j)]
      !
      ! and
      !
      ! [m2(k+1), m2(k)]
      !
      ! where
      !
      ! m_i(i) = Sum_{l >= i} dm_i(l)
      !
      ! The reason we specify dm's rather than m's is to avoid
      ! having to subtract large numbers and incur the resultant
      ! penalty in precision.
      !
      ! We assume dm1 and dm2 have the same length nz.
      ! We allow the indices j and k to range between 0 and nz-1
      ! inclusive. Here m(nz) is just the mass of the centre.
      !
      ! The width of the overlap is given by
      !
      ! min(m1(j), m2(k)) - max(m1(j+1), m2(k+1))
      !
      ! if positive. Otherwise the intervals do not intersect.
      ! We first compute the differences between the arguments of each
      ! min and max to determine the final difference.
      !
      ! There are nz*nz combinations of (j,k), but only <= 2*nz of them
      ! are non-zero. To see this assume  withouut loss of generality that
      ! m1(1) < m2(1). Augment m1 with m1(0) = m2(1). It follows that the
      ! overlap (0,1) is non-zero, as is the overlap (nz,nz). Because
      ! mass is monotonic, the (j,k) with non-zero overlaps form a sequence
      ! which is monotonic in both j and k. There are at most 2*nz points in
      ! such a sequence. To avoid allocating a quadratic amount of memory in nz 
      ! we therefore store just 2*nz overlaps. The locations of these overlaps are
      ! given by j == ranges(:,1), k == ranges(:,2). Some of these j's or k's
      ! will generally be zero, corresponding to the padding described above.
      subroutine make_compressed_intersect(dm1, dm2, nz, mesh_intersects, ranges)
         ! Inputs
         real(qp), dimension(:) :: dm1, dm2
         real(qp), dimension(:) :: mesh_intersects
         integer, dimension(:,:) :: ranges
         integer nz

         ! Intermediates
         integer side
         real(qp), dimension(:), allocatable :: remainders1, remainders2
         type(accurate_real) diff, dBottom, dTop, dBottomTop, dTopBottom
         integer i, j, k, counter

         allocate(remainders1(nz), remainders2(nz))

         j = nz
         k = nz
         counter = 2 * nz

         mesh_intersects = 0

         ! j tracks our position in dm1, k tracks our position in dm2.
         ! We begin with both equal to nz.
         ! We first compute the intersection between the nz cells.
         ! After that our algorithm is:
         ! 1. Begin with whichever side has its current cell has its top face deepest down.
         !    Decrement that index. If it cannot be decremented
         !    terminate, because all remaining overlaps will be zero.
         ! 2. Calculate the overlap between the new cell and the old one on the other side.
         ! 3. Return to 1.
         ! 
         ! Along the way we store overlaps in the order in which we encounter them.
         ! As discussed in the comments above, we know there will be precisely 2*nz of these.
         ! The indices of these are stored as ranges(counter,1) = j, ranges(counter,2) = k.
         ranges = -1

         ! remainders1(j) is the amount of mass in dm1(j) above all cells of dm2.
         ! remainders2(k) is the amount of mass in dm2(k) above all cells of dm1.
         ! Only one of these may have non-zero entries, and those entries form a contiguous
         ! block which, if non-empty, includes the surface cell.
         remainders1 = dm1
         remainders2 = dm2

         ! For convenience we track
         ! dTop = m1(j) - m2(k)
         ! dBottom = m1(j+1) - m2(k+1)
         ! dTopBottom = m1(j) - m2(k+1)
         ! dBottomTop = m1(j+1) - m2(k)
         dBottom = 0.0_qp
         dTop = dm1(nz) - dm2(nz)
         dBottomTop = -dm2(nz)
         dTopBottom = dm1(nz)

         if (dm1(j) > dm2(j)) then
            side = 1

            ranges(counter, 1) = j
            ranges(counter, 2) = k
            mesh_intersects(counter) = dm2(k)
            counter = counter - 1

            remainders1(j) = remainders1(j) - dm2(k)
            remainders2(j) = remainders2(j) - dm2(k)
            diff = dm1(j) - dm2(k)
         else
            side = 2

            ranges(counter, 1) = j
            ranges(counter, 2) = k
            mesh_intersects(counter) = dm1(j)
            counter = counter - 1

            remainders1(j) = remainders1(j) - dm1(j)
            remainders2(j) = remainders2(j) - dm1(j)
            diff = dm2(k) - dm1(j)
         end if      

         do while (k > 0 .and. j > 0 .and. ((side == 1 .and. k > 1) .or. (side == 2 .and. j > 1)))
            if (side == 1) then
               k = k - 1
               if (dm2(k) < diff) then
                  ranges(counter, 1) = j
                  ranges(counter, 2) = k
                  mesh_intersects(counter) = dm2(k)
                  counter = counter - 1

                  remainders1(j) = remainders1(j) - dm2(k)
                  remainders2(k) = remainders2(k) - dm2(k)
                  diff = diff - dm2(k)
               else
                  ranges(counter, 1) = j
                  ranges(counter, 2) = k
                  mesh_intersects(counter) = diff
                  counter = counter - 1

                  remainders1(j) = remainders1(j) - diff
                  remainders2(k) = remainders2(k) - diff
                  diff = dm2(k) - diff
                  side = 2
               end if

               dBottom = dBottom - dm2(k+1)
               dTopBottom = dTopBottom - dm2(k+1)
               dTop = dTop - dm2(k)
               dBottomTop = dBottomTop - dm2(k)

            else
               j = j -1
               if (dm1(j) < diff) then
                  ranges(counter, 1) = j
                  ranges(counter, 2) = k
                  mesh_intersects(counter) = dm1(j)
                  counter = counter - 1

                  remainders1(j) = remainders1(j) - dm1(j)
                  remainders2(k) = remainders2(k) - dm1(j)
                  diff = diff - dm1(j)
               else
                  ranges(counter, 1) = j
                  ranges(counter, 2) = k
                  mesh_intersects(counter) = diff
                  counter = counter - 1

                  remainders1(j) = remainders1(j) - diff
                  remainders2(k) = remainders2(k) - diff
                  diff = dm1(j) - diff
                  side = 1
               end if

               dBottom = dBottom + dm1(j+1)
               dBottomTop = dBottomTop + dm1(j+1)
               dTopBottom = dTopBottom + dm1(j)
               dTop = dTop + dm1(j)

            end if

         end do

         if (j > 1) then   
            do i=j,1,-1
               ranges(counter,1) = i
               ranges(counter,2) = 0
               mesh_intersects(counter) = remainders1(i)
               counter = counter - 1
            end do
         else
            do i=k,1,-1
               ranges(counter,1) = 0
               ranges(counter,2) = i
               mesh_intersects(counter) = remainders2(i)
               counter = counter - 1
            end do
         end if

      end subroutine make_compressed_intersect

      ! Computes for each final cell j the first (i_min(j)) and
      ! last (i_max(j)) initial cells whose overlap with j is non-zero.
      subroutine find_i_ranges(nz, ranges, i_min, i_max)
         ! Inputs
         integer nz
         integer, dimension(:,:) :: ranges
         integer, dimension(0:nz) :: i_min, i_max

         ! Intermediates
         integer i, j, counter

         ! Ensures that anything is less than initial i_min
         ! and anything is morre than initial i_max.
         do j=0,nz
            i_min(j) = nz+1
            i_max(j) = -1
         end do

         ! Find correct answers
         do counter=1,2*nz
            j = ranges(counter,1)
            i = ranges(counter,2)
            if (i < i_min(j)) i_min(j) = i
            if (i > i_max(j)) i_max(j) = i
         end do

      end subroutine find_i_ranges

      ! Computes the fraction of the mass which ends in cell j that passed through cell i
      ! for all (i,j) such that mesh_intersects(i,j) /= 0. These are the (i,j) for which
      ! the pass_fraction is not required to be 0 or 1, so we only need to precompute them.
      ! To do this we note that these (i,j) form a monotonic path from (0,0) to (nz,nz),
      ! and that within a row or column the path must be contiguous because cells represent
      ! contiguous mass intervals. This allows us to just sum along j in the direction the
      ! material flows to get the cumulative fraction of each cell which begins to one side
      ! of each face. Note that for cells where the flow converges (i.e. mass_flux(j) > 0,
      ! mass_flux(j+1) < 0) we perform a cumulative sum in both directions, meeting at
      ! cell j, at which point the two sums are added together.
      subroutine prepare_pass_fraction(nz, delta_m, dm, mesh_intersects, ranges, i_min, i_max, pf)
         ! Inputs
         integer nz
         real(qp) delta_m
         real(qp), dimension(:) :: dm, mesh_intersects
         integer, dimension(:,:) :: ranges
         integer, dimension(0:nz) :: i_min, i_max
         type(non_rect_array), dimension(0:nz) :: pf

         ! Intermediates
         integer i, j, l, counter, start, end, direction
         logical work_backwards
         real(qp), dimension(:), allocatable :: test

         ! Outputs
         call find_i_ranges(nz, ranges, i_min, i_max)

         do j=0,nz
            allocate(pf(j)%arr(i_min(j):i_max(j)))
            pf(j)%arr = 0
         end do

         ! Store mesh_intersects in pf
         counter = 1
         do j=0,nz
            do i=i_min(j),i_max(j)
               pf(j) % arr(i) = mesh_intersects(counter)
               counter = counter + 1
            end do
         end do

         ! Cumulatively sum mesh_intersects
         ! from i_min -> min(j, i_max) and from i_max -> max(j, i-min).
         do j=0,nz
            do i=i_min(j)+1,min(j,i_max(j))
               pf(j)%arr(i) = pf(j)%arr(i) + pf(j)%arr(i - 1)
            end do
            do i=i_max(j)-1,max(j,i_min(j)),-1
               pf(j)%arr(i) = pf(j)%arr(i) + pf(j)%arr(i + 1)               
            end do
         end do

         ! Normalize by cell mass.
         ! Note that j == 0 only has a non-zero pass fraction when delta_m < 0,
         ! in which case it sums over l to -delta_m.
         ! So in all cases we may normalize the j == 0 case by -delta_m.
         do j=0,nz
            do l=i_min(j), i_max(j)
               if (j > 0) then
                  pf(j)%arr(l) = pf(j)%arr(l) / dm(max(1,j)) ! bp: get rid of bogus compiler warning
               else
                  pf(j)%arr(l) = -pf(j)%arr(l) / delta_m
               end if
            end do
         end do

      end subroutine prepare_pass_fraction


      ! Computes the fraction of mass in final cell j which was at some point during adjust_mass
      ! inside cell i.
      real(qp) function compute_pass_fraction(i, j, i_min, i_max, pf)
         ! Inputs
         integer i, j
         integer, dimension(0:) :: i_min, i_max
         type(non_rect_array), dimension(0:) :: pf

         if (i >= i_min(j) .and. i <= i_max(j)) then
            compute_pass_fraction = pf(j) % arr(i) ! In the pre-computed portion
         else if (i > max(j,i_max(j)) .or. i < min(j,i_min(j))) then
            compute_pass_fraction = 0 ! Outside of the flow
         else
            compute_pass_fraction = 1 ! In the flow but not pre-computed
         end if

      end function compute_pass_fraction

      ! Computes for each i the minimum and maximum j such that pass_fraction(i,j) > 0.
      ! The set of non-zero pass_fraction(i,:) is guaranteed to be contiguous so this
      ! is equivalent to enumerating all (i,j) with pass_fraction(i,j) > 0.
      subroutine find_j_ranges(nz, ranges, mass_flux, j_min, j_max)
         ! Inputs
         integer nz
         integer, dimension(:,:) :: ranges
         integer, dimension(:) :: j_min, j_max
         real(qp), dimension(:) :: mass_flux

         ! Intermediates
         integer i, j, counter

         ! Ensures that anything is less than initial i_min
         ! and anything is morre than initial i_max.
         do i=1,nz
            j_min(i) = nz+1
            j_max(i) = -1
         end do

         ! If the flow is downward then j_min = i (all lesser j have i > j).
         ! If the flow is downward then j_max is the greatest j such that (i,j) appears
         ! in ranges. All greater j have vanishing pass_fraction.

         ! If the flow is upward then j_max = i (all greater j have i < j).
         ! If the flow is upward then j_min is the least j such that (i,j) appears in ranges.
         ! All lesser j have vanishing pass_fraction.

         ! If material enters a cell from both faces then all material which touches
         ! the cell ends in the cell, so pass_fraction(i,:) is non-zero only for j == i.
         ! Hence j_min = j_max = i.

         ! If material leaves a cell from both faces then j_min is,
         ! the least j such that (i,j) appears in ranges, and j_max is
         ! the greatest j such that (i,j) appears in ranges.

         ! Note that for all of the above we may use inclusive inequalities
         ! because when one or more of the faces of a cell have zero mass flux
         ! these different cases agree.

         ! Set the easy ones
         do i=1,nz
            if (mass_flux(i) >= 0 .and. mass_flux(i+1) >= 0) then
               j_min(i) = i
            end if

            if (mass_flux(i) <= 0 .and. mass_flux(i+1) <= 0) then
               j_max(i) = i
            end if

            if (mass_flux(i) >= 0 .and. mass_flux(i+1) <= 0) then
               j_min(i) = i
               j_max(i) = i
            end if

         end do

         ! Now go through ranges
         do counter=1,2*nz
            j = ranges(counter,1)
            i = ranges(counter,2)

            if (i > 0) then
               if (mass_flux(i) >= 0 .and. mass_flux(i+1) >= 0) then
                  if (j > j_max(i)) j_max(i) = j
               end if

               if (mass_flux(i) <= 0 .and. mass_flux(i+1) <= 0) then
                  if (j < j_min(i)) j_min(i) = j
               end if

               if (mass_flux(i) <= 0 .and. mass_flux(i+1) >= 0) then
                  if (j < j_min(i)) j_min(i) = j
                  if (j > j_max(i)) j_max(i) = j
               end if
            end if
         end do

      end subroutine find_j_ranges

      ! Integrates conserved quantities given in vals_old from the mesh specified by dm_old onto the mesh given by dm_new.
      ! vals_new is set equal to \int vals_old dm / dm_new, where the integral runs over the mass range of the new cell.
      ! vals_outside is the value to use in the integral when the new cell runs outside the old mesh.
      subroutine integrate_conserved(vals_new, vals_old, vals_outside, dm_new, dm_old, nz, nvars)

         ! Inputs
         integer, intent(in) :: nz, nvars
         real(dp), intent(in), dimension(:,:) :: vals_old ! (cell index, variable index)
         real(dp), intent(in), dimension(:) :: vals_outside ! (variable index)
         real(qp), intent(in), dimension(:) :: dm_new, dm_old

         ! Intermediates
         integer i, j, k, l
         real(qp) :: mesh_intersects(2*nz)
         integer :: ranges(2*nz,2)

         ! Output
         real(dp), intent(out), dimension(:,:) :: vals_new

         ranges = 0
         mesh_intersects = 0d0

         ! Intersect meshes
         call make_compressed_intersect(dm_new, dm_old, nz, mesh_intersects, ranges)

         vals_new = 0d0
         do i=1,2*nz
            j = ranges(i,1)
            k = ranges(i,2)

            if (j == 0) cycle

            do l=1,nvars
               if (k > 0) then
                  vals_new(j,l) = vals_new(j,l) + mesh_intersects(i) * vals_old(k,l)
               else
                  vals_new(j,l) = vals_new(j,l) + mesh_intersects(i) * vals_outside(l)
               end if
            end do

         end do

         do j=1,nz
            do l=1,nvars
               vals_new(j,l) = vals_new(j,l) / dm_new(j)
            end do
         end do

      end subroutine integrate_conserved


  end module mass_utils
