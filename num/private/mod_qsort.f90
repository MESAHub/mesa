

      module mod_qsort
      use const_def, only: dp
      
      
      implicit none


      contains


  
           ! FILE: sort.f
         ! PURPOSE: demonstrate the use of "qsort_inline.inc" and
         ! "qsort_inline_index.inc". These can be used as specific
         ! sort procedures under a common SORT generic name.
         !---------------------------------------------------------------
         ! Sort a string array, with any string length.
         subroutine sortp_string(array_size,index,string)
           integer, intent(in) :: array_size
           integer, intent(out) :: index(:) ! (array_size)
           character(len=*), intent(in) :: string(:) ! (array_size)
#include "qsort_inline.inc"
         contains
         ! set up initial index:
           subroutine init()
             integer :: i
             do i=1,array_size
               index(i)=i
             end do
           end subroutine init

         ! swap indices a,b
           subroutine swap(a,b)
             integer, intent(in) :: a,b
             integer :: hold
             hold=index(a)
             index(a)=index(b)
             index(b)=hold
           end subroutine swap

         ! circular shift-right by one:
           subroutine rshift(left,right)
             implicit none
             integer, intent(in) :: left, right
             integer :: hold, i
             hold=index(right)
             ! This sytnax is valid, but has poor optimization in GFortran:
             ! index(left+1:right)=index(left:right-1)
             do i=right,left+1,-1
               index(i)=index(i-1)
             end do
             index(left)=hold
           end subroutine rshift
  
           logical &
           function less_than(a,b)
             integer, intent(in) :: a,b
             if ( string(index(a)) == string(index(b))  ) then
               less_than = ( index(a) < index(b) )
             else
               less_than = ( string(index(a)) < string(index(b)) )
             end if
           end function less_than
           
         end subroutine sortp_string
         !---------------------------------------------------------------
         ! Sort an array of indices into a string array, with any string length.
         subroutine sortp_string_index(array_size,index,str_index,string)
           integer, intent(in) :: array_size
           integer, intent(out) :: index(:) ! (array_size)
           integer, intent(in) :: str_index(:) ! (array_size)
           character(len=*), intent(in) :: string(:) ! 1..maxval(str_index)
#include "qsort_inline.inc"
         contains
         ! set up initial index:
           subroutine init()
             integer :: i
             do i=1,array_size
               index(i)=i
             end do
           end subroutine init

         ! swap indices a,b
           subroutine swap(a,b)
             integer, intent(in) :: a,b
             integer :: hold
             hold=index(a)
             index(a)=index(b)
             index(b)=hold
           end subroutine swap

         ! circular shift-right by one:
           subroutine rshift(left,right)
             implicit none
             integer, intent(in) :: left, right
             integer :: hold, i
             hold=index(right)
             ! This sytnax is valid, but has poor optimization in GFortran:
             ! index(left+1:right)=index(left:right-1)
             do i=right,left+1,-1
               index(i)=index(i-1)
             end do
             index(left)=hold
           end subroutine rshift
           
           logical &
           function less_than(a,b)
             integer, intent(in) :: a,b
             if ( string(str_index(index(a))) == string(str_index(index(b)))  ) then
               less_than = ( str_index(index(a)) < str_index(index(b)) )
             else
               less_than = ( string(str_index(index(a))) < string(str_index(index(b))) )
             end if
           end function less_than
           
         end subroutine sortp_string_index
         !---------------------------------------------------------------
         ! Sort a double-precision array by index
         subroutine sortp_dp(array_size,index,value)
           integer, intent(in) :: array_size
           integer, intent(inout) :: index(:) ! (array_size)
           real(dp), intent(in) :: value(:) ! (array_size)
#include "qsort_inline.inc"
         contains
         ! set up initial index:
           subroutine init()
             integer :: i
             do i=1,array_size
               index(i)=i
             end do
           end subroutine init

         ! swap indices a,b
           subroutine swap(a,b)
             integer, intent(in) :: a,b
             integer :: hold
             hold=index(a)
             index(a)=index(b)
             index(b)=hold
           end subroutine swap

         ! circular shift-right by one:
           subroutine rshift(left,right)
             implicit none
             integer, intent(in) :: left, right
             integer :: hold, i
             hold=index(right)
             ! This sytnax is valid, but has poor optimization in GFortran:
             ! index(left+1:right)=index(left:right-1)
             do i=right,left+1,-1
               index(i)=index(i-1)
             end do
             index(left)=hold
           end subroutine rshift
  
           logical &
           function less_than(a,b)
             integer, intent(in) :: a,b
             less_than = value(index(a)) < value(index(b))
           end function less_than
           
         end subroutine sortp_dp


      
      
      end module mod_qsort
      
