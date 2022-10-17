! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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

module pulse_gr1d
   
   ! Uses
   
   use star_private_def
   use const_def
   use utils_lib
   
   use pulse_utils
   
   ! No implicit typing
   
   implicit none
   
   ! Access specifiers
   
   private
   
   public :: get_gr1d_data
   public :: write_gr1d_data

contains
   
   subroutine get_gr1d_data (id, global_data, point_data, ierr)
      
      integer, intent(in) :: id
      real(dp), allocatable, intent(out) :: global_data(:)
      real(dp), allocatable, intent(out) :: point_data(:, :)
      integer, intent(out) :: ierr
      
      type(star_info), pointer :: s
      integer :: nn
      integer :: j
      integer :: k
      integer :: sg
      
      ! Get model data for GR1D output
      
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) then
         write(*, *) 'bad star id for get_gyre_data'
         return
      end if
      
      ! Determine data dimensiones
      
      nn = s%nz
      
      ! Store global data
      
      allocate(global_data(0))
      
      ! Store point data
      
      allocate(point_data(7, nn))
      
      j = 1
      
      ! Envelope
      
      sg = 1
      
      env_loop : do k = s%nz, 1, -1
         
         call store_point_data_env(j, k)
      
      end do env_loop
      
      ! Finish
      
      return
   
   contains
      
      subroutine store_point_data_env (j, k)
         
         integer, intent(in) :: j
         integer, intent(in) :: k
         
         ! Store data associated with envelope face k into the point_data
         ! array at position j
         
         associate (&
            m => point_data(1, j), &
            r => point_data(2, j), &
            T => point_data(3, j), &
            rho => point_data(4, j), &
            v => point_data(5, j), &
            ye => point_data(6, j), &
            omega => point_data(7, j))
            
            m = s% m(k) - 0.5d0 * s% dm(k)
            
            T = s%T(k)
            rho = s%rho(k)
            ye = s%ye(k)
            
            if (s% rotation_flag) then
               if (k == s%nz) then
                  omega = s%omega(k)
               else
                  omega = 0.5d0 * (s%omega(k) + s%omega(k + 1))
               end if
            else
               omega = 0d0
            endif
            
            if (k == s% nz) then
               r = 0.5d0 * s%r(k)
               v = 0.5d0 * s%v(k)
            else
               r = 0.5d0 * (s%r(k) + s%r(k + 1))
               v = 0.5d0 * (s%v(k) + s%v(k + 1))
            end if
         
         end associate
         ! Finish
         
         return
      
      end subroutine store_point_data_env
   
   end subroutine get_gr1d_data
   
   !****
   
   subroutine write_gr1d_data (id, filename, global_data, point_data, ierr)
      
      integer, intent(in) :: id
      character(*), intent(in) :: filename
      real(dp), intent(in) :: global_data(:)
      real(dp), intent(in) :: point_data(:, :)
      integer, intent(out) :: ierr
      
      type(star_info), pointer :: s
      integer :: iounit
      integer :: nn
      integer :: j
      
      ! Write GR1D data to file
      
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) then
         write(*, *) 'bad star id for write_gr1d_data'
         return
      end if
      
      ! Open the file
      
      open(newunit = iounit, file = TRIM(filename), status = 'REPLACE', iostat = ierr)
      if (ierr /= 0) then
         write(*, *) 'failed to open ' // TRIM(filename)
         return
      end if
      
      ! Write the data
      
      nn = SIZE(point_data, 2)
      
      write(iounit, 100) nn
      100 format(I6)
      
      do j = 1, nn
         write(iounit, 100) j, point_data(:, j)
         110    format(I6, 1P, 99E20.10)
      end do
      
      ! Close the file
      
      close(iounit)
      
      ! Finish
      
      return
   
   end subroutine write_gr1d_data

end module pulse_gr1d
