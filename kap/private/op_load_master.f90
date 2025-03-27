! ***********************************************************************
!
!   Copyright (C) 2013-2019  The MESA Team
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   https://http://mesastar.org/
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

module op_load_master

   use const_def, only: dp

   implicit none

   private
   public :: load_op_master

   logical :: loaded_op_master = .false.

contains

   subroutine load_op_master(emesh_data_for_op_mono_path, iz, ite, jne, epatom, amamu, sig, eumesh, ierr)

      character(len=*), intent(in) :: emesh_data_for_op_mono_path

      integer, intent(inout) :: ierr
      integer, pointer, intent(out) :: iz(:), ite(:), jne(:)
      real(dp), pointer, intent(out) :: sig(:, :, :)
      real(dp), pointer, intent(out):: epatom(:, :), amamu(:), eumesh(:, :, :)

      integer :: n, m, ke
      CHARACTER(LEN=72) :: FMT
      integer :: nel, nptot, np
      parameter(nel=17, nptot=10000, np=1648)  ! number of elements and number of u-mesh points
      real(dp), allocatable :: amamu_f(:, :)
      integer, allocatable  :: iz_f(:, :)

      if (loaded_op_master) return

      allocate (iz_f(nel, np), iz(nel), ite(np), jne(np), stat=ierr)
      allocate (sig(nel, np, nptot), stat=ierr)
      allocate (epatom(nel, np), amamu_f(nel, np), amamu(nel), eumesh(nel, np, nptot), stat=ierr)

      FMT = '(i2,1x,i3,1x,i3,1x,F14.10,1x,F14.10,10000(1x,E12.6E3),10000(1x,E13.6E3))'

      write (*, *) 'Opening file...'
      open (1, file=emesh_data_for_op_mono_path, form='formatted', action='read')
      write (*, *) 'Loading OP mono data...'

      do ke = 1, nel
         do n = 1, np
            read (1, FMT)iz_f(ke,n),ite(n),jne(n),epatom(ke,n),amamu_f(ke,n),(sig(ke,n,m), m=1,nptot),(eumesh(ke,n,m), m=1,nptot)
         end do
      end do

      close (1)

      do ke = 1, nel
         amamu(ke) = amamu_f(ke, 1)
         iz(ke) = iz_f(ke, 1)
      end do

      write (*, *) 'OP mono data loaded.'
      ierr = 0
      loaded_op_master = .true.

   end subroutine load_op_master

end module op_load_master
