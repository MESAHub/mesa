! ***********************************************************************
!
!   Copyright (C) 2011-2019  The MESA Team
!
!   This file is part of MESA.
!
!   MESA is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

program eos_correction
   use eos_def
   use eos_lib
   use chem_def
   use chem_lib
   use const_lib
   use math_lib
   
   implicit none
   
   integer :: handle
   real(dp) :: X, Z, Y, abar, zbar, z2bar, z53bar, ye, WoA
   integer, parameter :: species = 2
   integer, parameter :: h1 = 1, c12 = 2
   integer, pointer, dimension(:) :: net_iso, chem_id
   real(dp) :: xa(species)
   
   integer, parameter :: num_lgRhos = 601
   integer, parameter :: num_lgTs = 321
   real(dp), parameter :: lg_Tmin = 6.0_dp, lg_Tmax = 9.2_dp, &
      lg_Rhomin = 4.0_dp, lg_Rhomax = 10.0_dp
   real(dp), dimension(num_lgTs) :: lg_Ts
   real(dp), dimension(num_lgRhos) :: lg_Rhos
   real(dp), dimension(num_lgTs, num_lgRhos) :: tab, Ytab, Etab
   integer :: i
   character (len = 256) :: my_mesa_dir
   
   call setup
   
   do i = 1, num_lgts
      lg_Ts(i) = lg_Tmax + real(i - 1) * (lg_Tmin - lg_Tmax) / real(num_lgTs - 1)
   end do
   do i = 1, num_lgrhos
      lg_Rhos(i) = lg_Rhomin + real(i - 1) * (lg_Rhomax - lg_Rhomin) / real(num_lgRhos - 1)
   end do
   call write_axes_to_file
   
   X = 0.0_dp; Y = 0.0_dp; Z = 1.0_dp
   call make_correction_table(X, Z, tab, Ytab, Etab)
   call write_correction_table('correction_C12.data', tab)
   call write_correction_table('Yfree_C12.data', Ytab)
   call write_correction_table('EoC2_C12.data', Etab)
   
   call shutdown

contains
   
   subroutine setup()
      use math_lib
      integer :: ierr
      
      ierr = 0
      my_mesa_dir = '..' ! if empty string, uses environment variable MESA_DIR
      call const_init(my_mesa_dir, ierr)
      if (ierr /= 0) then
         write(*, *) 'const_init failed'
         call mesa_error(__FILE__, __LINE__)
      end if
      
      call math_init()
      
      call chem_init('isotopes.data', ierr)
      if (ierr /= 0) then
         write(*, *) 'failed in chem_init'
         call mesa_error(__FILE__, __LINE__)
      end if
      
      ! allocate and initialize the eos tables
      call Setup_eos(handle)
      
      allocate(net_iso(num_chem_isos), chem_id(species), stat = ierr)
      if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'allocate failed')
   end subroutine setup
   
   subroutine make_correction_table(X, Z, Ecorr, Ytab, Etab)
      implicit none
      real(dp), intent(in) :: X, Z
      real(dp), dimension(num_lgTs, num_lgRhos) :: Ecorr, Ytab, Etab
      real(dp) :: Rho, T, log10Rho, log10T
      real(dp), dimension(num_eos_basic_results) :: res, d_dlnd, d_dlnT
      real(dp), dimension(num_eos_d_dxa_results, species) :: d_dxa
      real(dp) :: Yplus, Yfree, Eoc2
      integer :: i, j, ierr
      
      call Init_Composition(X, Z)
      
      do i = 1, num_lgRhos
         log10Rho = lg_Rhos(i)
         Rho = exp10(log10Rho)
         do j = 1, num_lgTs
            log10T = lg_Ts(j)
            T = exp10(log10T)
            
            call eosDT_get(&
               handle, &
               species, chem_id, net_iso, xa, &
               Rho, log10Rho, T, log10T, res, d_dlnd, d_dlnT, &
               d_dxa, ierr)
            
            Yfree = exp(res(i_lnfree_e))
            Yplus = max(Yfree - ye, 0.0_dp)
            Eoc2 = exp(res(i_lnE)) / (clight * clight)
            
            Ecorr(j, i) = WoA + Yplus * me / amu + Eoc2
            Ytab(j, i) = Yplus
            Etab(j, i) = Eoc2
         end do
      end do
      write (*, '(a16,"=",2f13.6)') 'Ecorr min, max', minval(Ecorr(1:num_lgTs, 1:num_lgRhos)), &
         maxval(Ecorr(1:num_lgTs, 1:num_lgRhos))
      write (*, '(a16,"=",2f13.6)') 'Yfree min, max', minval(Ytab(1:num_lgTs, 1:num_lgRhos)), &
         maxval(Ytab(1:num_lgTs, 1:num_lgRhos))
      write (*, '(a16,"=",2f13.6)') 'E/c**2 min, max', minval(Etab(1:num_lgTs, 1:num_lgRhos)), &
         maxval(Etab(1:num_lgTs, 1:num_lgRhos))
   end subroutine make_correction_table
   
   
   subroutine shutdown()
      ! deallocate the eos tables
      call Shutdown_eos(handle)
      
      deallocate(net_iso, chem_id)
   
   end subroutine shutdown
   
   
   subroutine Setup_eos(handle)
      ! allocate and load the eos tables
      integer, intent(out) :: handle
      
      integer :: ierr
      logical, parameter :: use_cache = .true.
      
      call eos_init('', use_cache, ierr)
      if (ierr /= 0) then
         write(*, *) 'eos_init failed in Setup_eos'
         call mesa_error(__FILE__, __LINE__)
      end if
      
      write(*, *) 'loading eos tables'
      
      handle = alloc_eos_handle(ierr)
      if (ierr /= 0) then
         write(*, *) 'failed trying to allocate eos handle'
         call mesa_error(__FILE__, __LINE__)
      end if
   
   end subroutine Setup_eos
   
   
   subroutine Shutdown_eos(handle)
      use eos_def
      use eos_lib
      integer, intent(in) :: handle
      call free_eos_handle(handle)
      call eos_shutdown
   end subroutine Shutdown_eos
   
   
   subroutine Init_Composition(X, Z)
      use chem_lib
      real(dp), intent(in) :: X, Z
      real(dp) :: dabar_dx(species), dzbar_dx(species), dmc_dx(species), sumx, xh, xhe, xz
      
      xa(h1) = X
      xa(c12) = Z
      
      net_iso(:) = 0
      chem_id(h1) = ih1; net_iso(ih1) = h1
      chem_id(c12) = ic12; net_iso(ic12) = c12
      
      call composition_info(&
         & species, chem_id, xa, xh, xhe, xz, abar, zbar, z2bar, z53bar, ye, &
         & WoA, sumx, dabar_dx, dzbar_dx, dmc_dx)
   
   end subroutine Init_Composition
   
   subroutine write_axes_to_file()
      use utils_lib
      character(len = *), parameter :: form = '(f6.2)'
      integer :: iounit, ierr
      
      open(newunit = iounit, file = 'data/lgTs', iostat = ierr, action = "write")
      if (ierr /= 0) then
         write(*, *) "Error opening file data/lgTs"
         stop
      end if
      
      write(iounit, form) lg_Ts
      close(iounit)
      
      open(newunit = iounit, file = 'data/lgRhos', iostat = ierr, action = "write")
      if (ierr /= 0) then
         write(*, *) "Error opening file data/lgRhos"
         stop
      end if
      
      write(iounit, form) lg_Rhos
      close(iounit)
   
   end subroutine write_axes_to_file
   
   subroutine write_correction_table(filename, tab)
      use utils_lib
      character(len = *), intent(in) :: filename
      real(dp), dimension(num_lgTs, num_lgRhos), intent(in) :: tab
      integer :: iounit, ierr
      character(len = 32) :: form
      
      open(newunit = iounit, file = 'data/' // trim(filename), iostat = ierr, action = "write")
      if (ierr /= 0) then
         write (*, *) "Error opening file data/Ecorrection.data"
         stop
      end if
      
      write(form, '("(",i0,"f11.6)")') num_lgRhos
      
      do i = 1, num_lgTs
         write(iounit, form) tab(i, :)
      end do
      
      close(iounit)
   end subroutine write_correction_table

end program eos_correction
   
