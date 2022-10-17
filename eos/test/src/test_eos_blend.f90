module test_eos_blend
   
   use const_def, only : dp
   use auto_diff
   use eos_blend
   
   implicit none
   private
   public :: do_test_eos_blend

contains
   
   subroutine do_test_eos_blend()
      integer :: w
      real(dp) :: coords(4, 2)
      type(auto_diff_real_2var_order1) :: p(2), d
      
      coords(1, 1) = -1d0
      coords(1, 2) = -1d0
      coords(2, 1) = -1d0
      coords(2, 2) = 1d0
      coords(3, 1) = 1d0
      coords(3, 2) = 1d0
      coords(4, 1) = 1d0
      coords(4, 2) = -1d0
      
      p(1) = 0d0
      p(2) = 1d0
      write(*, *) is_contained(4, coords, p)
      d = min_distance_to_polygon(4, coords, p)
      write(*, '(99(1pd26.16))') d
      p(1) = -1d0
      p(2) = 1d0
      write(*, *) is_contained(4, coords, p)
      d = min_distance_to_polygon(4, coords, p)
      write(*, '(99(1pd26.16))') d
      p(1) = 0d0
      p(2) = 0d0
      write(*, *) is_contained(4, coords, p)
      d = min_distance_to_polygon(4, coords, p)
      write(*, '(99(1pd26.16))') d
      p(1) = -2d0
      p(2) = 0d0
      write(*, *) is_contained(4, coords, p)
      d = min_distance_to_polygon(4, coords, p)
      write(*, '(99(1pd26.16))') d
      p(1) = -1d0
      p(2) = 3d0
      write(*, *) is_contained(4, coords, p)
      d = min_distance_to_polygon(4, coords, p)
      write(*, '(99(1pd26.16))') d
   
   end subroutine do_test_eos_blend

end module test_eos_blend
