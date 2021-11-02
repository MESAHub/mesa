! ***********************************************************************
!
!   Copyright (C) 2018 Robert Farmer
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

module utils_system
   implicit none

   interface 
      function f_mkdir_p(folder) bind(C,name='c_mkdir_p')
         use, intrinsic :: ISO_C_BINDING, only: C_CHAR, C_INT
         integer(C_INT) :: f_mkdir_p
         character(kind=C_CHAR) :: folder(*)
      end function f_mkdir_p
   
      function f_mv(src, dest) bind(C,name='c_mv')
         use, intrinsic :: ISO_C_BINDING, only: C_CHAR, C_INT
         integer(C_INT) :: f_mv
         character(kind=C_CHAR) :: src(*), dest(*)
      end function f_mv
      
      function f_cp(src, dest) bind(C,name='c_cp')
         use, intrinsic :: ISO_C_BINDING, only: C_CHAR, C_INT
         integer(C_INT) :: f_cp
         character(kind=C_CHAR) :: src(*), dest(*)
      end function f_cp

      function f_is_dir(folder) bind(C,name='is_dir')
         use, intrinsic :: ISO_C_BINDING, only: C_CHAR, C_INT
         integer(C_INT) :: f_is_dir
         character(kind=C_CHAR) :: folder(*)
      end function f_is_dir

   end interface

   private 
   public :: mkdir_p, mv, cp, is_dir


   contains
   
   
   ! Converts a fortran string to a NULL terminated string 
   pure function f_c_string (f_str) result (c_str)
      use, intrinsic :: ISO_C_BINDING, only: C_CHAR, C_NULL_CHAR
      character(len=*), intent(in) :: f_str
      character(len=1,kind=C_CHAR) :: c_str(len_trim(f_str)+1)
      integer                      :: n, i
      
      n = len_trim(f_str)
      do i = 1, n
         c_str(i) = f_str(i:i)
      end do
      c_str(n + 1) = C_NULL_CHAR
   
   end function f_c_string 
   
   ! Makes a directory, potentially making any needed parent directories
   integer function mkdir_p(folder)
      character(len=*), intent(in) :: folder

      mkdir_p = f_mkdir_p(f_c_string(folder))
   
   end function mkdir_p
   
   ! Moves src to dest, if dest is on a different filesystem, do a cp
   ! to the same filesystem then mv to dest
   integer function mv(src,dest)
      character(len=*), intent(in) :: src, dest
      
      mv = f_mv(f_c_string(src),f_c_string(dest))
   
   end function mv
   
   ! Copies src to dest
   integer function cp(src,dest)
      character(len=*), intent(in) :: src, dest
      
      cp = f_cp(f_c_string(src),f_c_string(dest))
   
   end function cp

   ! Checks if folder exists or not
   logical function is_dir(folder)
      character(len=*), intent(in) :: folder

      is_dir = f_is_dir(folder) == 1

   end function is_dir


end module utils_system


! Left for testing

!program sys
!   use system_utils
!   implicit none
!   integer :: num, res
!   character(len=256) :: f1, f2
   
!   num = command_argument_count()
!   call get_command_argument(1,f1)
   
!   if(num==1) then
!      write(*,*) "Test mkdir_p ",trim(f1)
!      res = mkdir_p(f1)
!   else
!      call get_command_argument(2,f2)
!!      write(*,*) "Test mv ",trim(f1)," * ",trim(f2)
!!      res = mv(f1,f2)
!      write(*,*) "Test cp ",trim(f1)," * ",trim(f2)
!      res = cp(f1,f2)
!   end if
   
!   write(*,*) "Result: ", res

!end program sys
