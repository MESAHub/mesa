module eos_blend
      use const_def, only: dp
      use math_lib
      use auto_diff
      use eos_def

      implicit none
      private
      public :: is_contained, min_distance_to_polygon

      contains

      !! Determines which quadrant the given point is in.
      !!
      !! @param p The coordinates of the point (x,y).
      !! @param q The quadrant index.
      integer function quadrant(p) result(q)
         real(dp), intent(in) :: p(2)

         if (p(1) >= 0d0 .and. p(2) >= 0d0) then
            q = 1
         else if (p(1) < 0d0 .and. p(2) >= 0d0) then
            q = 2
         else if (p(1) < 0d0 .and. p(2) < 0d0) then
            q = 3
         else
            q = 4
         end if

      end function quadrant

      !! Determines the winding number of a polygon around the origin.
      !! 
      !! Implements the winding number algorithm of
      !! Moscato, Titolo, Feliu, and Munoz (https://shemesh.larc.nasa.gov/people/cam/publications/FM2019-draft.pdf)
      !!
      !! @param num_points The number of points specifying the polygon.
      !! @param coords The coordinates of the polygon. Given as an array of shape (num_coords, 2) storing (x,y) pairs.
      !! @param w The winding number of the polygon about the origin.
      integer function winding_number(num_points, coords) result(w)
         integer, intent(in) :: num_points
         real(dp), intent(in) :: coords(num_points, 2)
         real(dp) :: determinant, x_start(2), x_end(2), dist(2)
         integer :: i, q_start, q_end, increment

         w = 0
         do i=1,num_points
            x_start = coords(i,1:2)

            if (i == num_points) then
               x_end = coords(1,1:2)
            else
               x_end = coords(i+1,1:2)
            end if

            dist(1) = x_end(1) - x_start(1)
            dist(2) = x_end(2) - x_start(2)

            determinant = x_start(2) * dist(1) - x_start(1) * dist(2)

            q_start = quadrant(x_start)
            q_end = quadrant(x_end)

            if (q_start == q_end) then
               increment = 0
            else if (q_end == mod(q_start,4) + 1) then
               increment = 1
            else if (q_start == mod(q_end, 4) + 1) then
               increment = -1
            else if (determinant < 0d0) then
               increment = 2
            else
               increment = -2
            end if

            w = w + increment
         end do
      end function winding_number

      !! Determines if a point is contained within a polygon specified by connecting a sequence of points.
      !! This is done by checking the winding number of the polygon around the point.
      !!
      !! @param num_points The number of points specifying the polygon.
      !! @param coords The coordinates of the polygon. An array of shape (num_points, 2) storing (x,y) pairs.
      !! @param p The point whose distance to compute.
      !! @param contained Whether or not the point p is contained in the polygon.
      logical function is_contained(num_points, coords, p) result(contained)
         integer, intent(in) :: num_points
         real(dp), intent(in) :: coords(num_points,2)
         type(auto_diff_real_2var_order1), intent(in) :: p(2)

         integer :: i, w
         real(dp) :: diff(num_points,2)

         do i=1,num_points
            diff(i,1) = coords(i,1) - p(1)%val
            diff(i,2) = coords(i,2) - p(2)%val
         end do

         w = winding_number(num_points, diff)
         if (w /= 0) then
            contained = .true.
         else
            contained = .false.
         end if

      end function is_contained

      !! Computes the minimum distance from a given point to a given line segment.
      !! 
      !! @param line_start The coordinates of the start of the line segment (x,y).
      !! @param line_end The coordinates of the end of the line segment (x,y).
      !! @param p The point whose distance to compute.
      type(auto_diff_real_2var_order1) function min_distance_from_point_to_line_segment(line_start, line_end, p) result(d)
         real(dp), intent(in) :: line_start(2), line_end(2)
         type(auto_diff_real_2var_order1), intent(in) :: p(2)
         real(dp), parameter :: eps = 1e-10 ! To avoid singularity in the derivatives near edges and corners.

         type(auto_diff_real_2var_order1) :: diff_line(2), diff_start(2), diff_end(2)
         type(auto_diff_real_2var_order1) :: length_squared, lambda, nearest_point_on_line(2)


         diff_start(1) = p(1) - line_start(1)
         diff_start(2) = p(2) - line_start(2)

         diff_end(1) = p(1) - line_end(1)
         diff_end(2) = p(2) - line_end(2)

         diff_line(1) = line_end(1) - line_start(1)
         diff_line(2) = line_end(2) - line_start(2)

         length_squared = pow2(diff_line(1)) + pow2(diff_line(2))

         ! First, find the nearest_point_on_line. We parameterize this by
         ! (nearest) = (start) + lambda (end - start)
         !
         ! We can then pretend the line is infinite, solve for lambda, and then restrict it to lie in [0,1].
         lambda = (diff_start(1) * diff_line(1) + diff_start(2) * diff_line(2)) / length_squared

         if (lambda < 0d0) then ! Nearest point is line_start
            d = sqrt(pow2(diff_start(1)) + pow2(diff_start(2)))
         else if (lambda > 1d0) then ! Nearest point is line_end
            d = sqrt(pow2(diff_end(1)) + pow2(diff_end(2)))
         else ! Nearest point is interior to the line segment
            nearest_point_on_line(1) = line_start(1) + lambda * diff_line(1)
            nearest_point_on_line(2) = line_start(2) + lambda * diff_line(2)
            nearest_point_on_line(1) = nearest_point_on_line(1) - p(1)
            nearest_point_on_line(2) = nearest_point_on_line(2) - p(2)
            d = sqrt(pow2(nearest_point_on_line(1)) + pow2(nearest_point_on_line(2)) + eps**2)
         end if

      end function min_distance_from_point_to_line_segment

      !! Computes the distance to the nearest line segment.
      !! This is done by looping over segments, computing the minimum distance to each, and 
      !! returning the smallest of those differences.
      !! 
      !! @param num_points The number of points specifying the polygon.
      !! @param coords The coordinates of the polygon. An array of shape (num_points, 2) storing (x,y) pairs.
      !! @param p The point whose distance to compute.
      type(auto_diff_real_2var_order1) function min_distance_to_polygon(num_points, coords, p) result(d)
         integer, intent(in) :: num_points
         real(dp), intent(in) :: coords(num_points, 2)
         type(auto_diff_real_2var_order1), intent(in) :: p(2)

         real(dp) :: line_start(2), line_end(2)
         type(auto_diff_real_2var_order1) :: line_dist
         integer :: i

         d = 1d99
         do i=1,num_points
            line_start(1) = coords(i,1)
            line_start(2) = coords(i,2)

            if (i == num_points) then
               line_end(1) = coords(1,1)
               line_end(2) = coords(1,2)
            else
               line_end(1) = coords(i+1,1)
               line_end(2) = coords(i+1,2)
            end if

            line_dist = min_distance_from_point_to_line_segment(line_start, line_end, p)
            d = min(d, line_dist)
         end do

      end function min_distance_to_polygon

end module eos_blend
