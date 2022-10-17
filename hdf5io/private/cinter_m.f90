! Module   : cinter_m
! Purpose  : C interoperability support
!
! Copyright 2021 Rich Townsend
!
! This file is part of the ForUM (Fortran Utility Modules)
! package. ForUM is free software: you can redistribute it and/or
! modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, version 3.
!
! ForUM is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
! License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
























module cinter_m
   
   ! Uses
   
   use kinds_m
   
   use ISO_FORTRAN_ENV
   use ISO_C_BINDING
   
   ! No implicit typing
   
   implicit none
   
   ! Interfaces
   
   interface c_f_len
      module procedure c_f_len_arr_
      module procedure c_f_len_ptr_
   end interface c_f_len
   
   interface c_f_string
      module procedure c_f_string_arr_
      module procedure c_f_string_ptr_
   end interface c_f_string
   
   ! Access specifiers
   
   private
   
   public :: c_f_len
   public :: c_f_string
   
   ! Procedures

contains
   
   function c_f_len_arr_(c_str) result (f_len)
      
      character(C_CHAR) :: c_str(*)
      integer :: f_len
      
      integer :: i
      
      ! Determine the Fortran length of the C string (expressed as an
      ! assumed-size array of C_CHARs, terminated by a NULL)
      
      i = 1
      
      len_loop : do
         if (c_str(i) == C_NULL_CHAR) exit len_loop
         i = i + 1
      end do len_loop
      
      f_len = i - 1
      
      ! Finish
      
      return
   
   end function c_f_len_arr_
   
   !****
   
   function c_f_len_ptr_(c_str) result (f_len)
      
      type(C_PTR), value :: c_str
      integer :: f_len
      
      character(C_CHAR), pointer :: p(:)
      integer :: i
      
      ! Determine the Fortran length of the C string (expressed as a
      ! pointer to a NULL-termninated sequence)
      
      if (C_ASSOCIATED(c_str)) then
         
         call c_f_pointer(c_str, p, [HUGE(0)])
         
         i = 1
         
         len_loop : do
            if (p(i) == C_NULL_CHAR) exit len_loop
            i = i + 1
         end do len_loop
         
         f_len = i - 1
      
      else
         
         f_len = 0
      
      endif
      
      ! Finish
      
      return
   
   end function c_f_len_ptr_
   
   !****
   
   function c_f_string_arr_(c_str) result (f_str)
      
      character(C_CHAR) :: c_str(*)
      character(:), allocatable :: f_str
      
      integer :: n
      integer :: i
      
      ! Convert the C string (expressed as an assumed-size array of
      ! C_CHARs, terminated by a NULL) into a Fortran string
      
      n = c_f_len(c_str)
      
      allocate(character(LEN = c_f_len(c_str)) :: f_str)
      
      copy_loop : do i = 1, n
         f_str(i:i) = c_str(i)
      end do copy_loop
      
      ! Finish
      
      return
   
   end function c_f_string_arr_
   
   !****
   
   function c_f_string_ptr_(c_str) result (f_str)
      
      type(C_PTR), value :: c_str
      character(:), allocatable :: f_str
      
      character(C_CHAR), pointer :: p(:)
      integer :: n
      integer :: i
      
      ! Convert the C string (expressed as a pointer to a
      ! NULL-terminated sequence) into a Fortran string
      
      if (C_ASSOCIATED(c_str)) then
         
         call c_f_pointer(c_str, p, [HUGE(0)])
         
         n = c_f_len(p)
         
         allocate(character(LEN = n) :: f_str)
         
         copy_loop : do i = 1, n
            f_str(i:i) = p(i)
         end do copy_loop
      
      else
         
         allocate(character(LEN = 0) :: f_str)
      
      endif
      
      ! Finish
      
      return
   
   end function c_f_string_ptr_

end module cinter_m
