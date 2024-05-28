!  
! lib_util
! francois hebert, 10/08/11   
!
! this subroutine contains "utility" routines
!

      module lib_util

      use def_type

      implicit none

      contains

      subroutine reallocate_dp(array, newsize)
         integer, intent(in) :: newsize
         real (dp), intent(inout), dimension(:), allocatable :: array

         integer :: oldsize
         real (dp), dimension(:), allocatable :: brray

         ! if array had already been allocated, then we:
         !  -save its data
         !  -deallocate and reallocate to new length
         !  -copy data within limitations of size
         if (allocated(array)) then
            oldsize = size(array)
            if (oldsize == newsize) then
               return
            else
               allocate(brray(oldsize))
               brray(:) = array(:)

               deallocate(array)
               allocate(array(newsize))

               if (oldsize >= newsize) then
                  array(:) = brray(1:newsize)
               else
                  array(1:oldsize) = brray(:)
               end if
            end if
         ! if array was not previously allocated, allocate it
         else
            allocate(array(newsize))
         end if
      end subroutine reallocate_dp


      end module lib_util
