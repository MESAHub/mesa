! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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

      module mod_ionization

      use const_def, only: one_third, two_thirds, avo, dp, mesa_data_dir
      use math_lib
      use ionization_def
      use utils_lib, only: mesa_error

      implicit none
      

      logical, parameter :: dbg = .false.
     

      contains


      subroutine do_init_ionization(ionization_cache_dir_in, use_cache, ierr)      
         character (len=*), intent(in) :: ionization_cache_dir_in
         logical, intent(in) :: use_cache
         integer, intent(out) :: ierr
         ierr = 0
         table_is_initialized = .false.
      end subroutine do_init_ionization


      subroutine do_load(ierr)
         use ionization_def
         integer, intent(out) :: ierr
      
         integer :: io_log_ne, io_logT, io_z
         integer, pointer :: ibound(:,:), tmp_version(:)
         integer, parameter :: num_log_ne_fe56_he4 = 105, num_logT_fe56_he4 = 30
      
         ierr = 0

         fe_he_ptr => fe_he_info  
      
         call load_table_summary( &
            'log_ne_fe56_he4.data', 'logT_fe56_he4.data', 'z_fe56_he4.data', &
            num_log_ne_fe56_he4, num_logT_fe56_he4, fe_he_ptr, ierr)
         if (ierr /= 0) return         

         call create_interpolants(fe_he_ptr,num_log_ne_fe56_he4,num_logT_fe56_he4,ierr)
         if (ierr /= 0) return         
         
         table_is_initialized = .true.
         
         contains
         
         subroutine openfile(filename, iounit, ierr)
            character(len=*) :: filename
            integer, intent(inout) :: iounit
            integer, intent(out) :: ierr
            if (dbg) write(*,*) 'read ' // trim(filename)            
            ierr = 0
            open(newunit=iounit,file=trim(filename),action='read',status='old',iostat=ierr)
            if (ierr/= 0) then
               write(*,*) 'table_ionization_init: missing ionization data'
               write(*,*) filename
               write(*,*)
               write(*,*)
               write(*,*)
               write(*,*)
               write(*,*)
               write(*,*) 'FATAL ERROR: missing or bad ionization data.'
               write(*,*) 'Please update by removing the directory mesa/data/ionization_data,'
               write(*,*) 'and rerunning the mesa ./install script.'
               write(*,*)
               call mesa_error(__FILE__,__LINE__)
            endif
         end subroutine openfile
      
      
         subroutine load_table_summary( &
               log_ne_fname, logT_fname, z_fname, num_log_ne, num_logT, p, ierr)
            character(len=*), intent(in) :: log_ne_fname, logT_fname, z_fname
            integer, intent(in) :: num_log_ne, num_logT
            type (Ionization_Info), pointer :: p
            integer, intent(out) :: ierr
      
            character(len=256) :: filename
            real(dp), pointer :: f(:,:,:)
            integer :: i, j
            
            ierr = 0
            p% have_interpolation_info = .false.
            p% num_log_ne = num_log_ne
            p% num_logT = num_logT
            allocate( &
               p% log_ne(num_log_ne), p% logT(num_logT), &
               p% f1(4*num_log_ne*num_logT), stat=ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in allocate for ionization tables'
               call mesa_error(__FILE__,__LINE__)
            end if
            f(1:4,1:num_log_ne,1:num_logT) => p% f1(1:4*num_log_ne*num_logT)
            
            filename = trim(mesa_data_dir) // '/ionization_data/' // trim(z_fname)            
            call openfile(filename, io_z, ierr)
            if (ierr /= 0) return
            do i=1,num_logT
               read(io_z,fmt=*,iostat=ierr) p% log_ne(1:num_log_ne)
               if (ierr /= 0) then
                  write(*,*) 'failed in reading ionization z ' // trim(filename)
                  call mesa_error(__FILE__,__LINE__)
               end if
               !p% f(1,1:num_log_ne,i) = p% log_ne(1:num_log_ne)  << segfault on UBUNTU
               do j=1,num_log_ne
                  f(1,j,i) = p% log_ne(j) ! sets p% f1
               end do
            end do
            close(io_z)
            
            filename = trim(mesa_data_dir) // '/ionization_data/' // trim(log_ne_fname)            
            call openfile(filename, io_log_ne, ierr)
            if (ierr /= 0) return
            do i=1,num_log_ne
               read(io_log_ne,fmt=*,iostat=ierr) p% log_ne(i)
               if (ierr /= 0) then
                  write(*,*) 'failed in reading ionization log_ne ' // trim(filename)
                  call mesa_error(__FILE__,__LINE__)
               end if
            end do
            close(io_log_ne)
            
            filename = trim(mesa_data_dir) // '/ionization_data/' // trim(logT_fname)            
            call openfile(filename, io_logT, ierr)
            if (ierr /= 0) return
            do i=1,num_logT
               read(io_logT,fmt=*,iostat=ierr) p% logT(i)
               if (ierr /= 0) then
                  write(*,*) 'failed in reading ionization logT ' // trim(filename)
                  call mesa_error(__FILE__,__LINE__)
               end if
            end do
            close(io_logT)
         
         end subroutine load_table_summary


      end subroutine do_load

      
      subroutine create_interpolants(p,nx,ny,ierr)
         use interp_2d_lib_db
         type (Ionization_Info), pointer :: p
         integer, intent(in) :: nx, ny
         integer, intent(out) :: ierr
         integer :: ibcxmin, ibcxmax, ibcymin, ibcymax
         real(dp) :: bcxmin(ny), bcxmax(ny), bcymin(nx), bcymax(nx)
         ! use "not a knot" bc's
         ibcxmin = 0; bcxmin(:) = 0d0
         ibcxmax = 0; bcxmax(:) = 0d0
         ibcymin = 0; bcymin(:) = 0d0
         ibcymax = 0; bcymax(:) = 0d0
         call interp_mkbicub_db( &
            p% log_ne, p% num_log_ne, p% logT, p% num_logT, &
            p% f1, p% num_log_ne, &
            ibcxmin, bcxmin, ibcxmax, bcxmax, &
            ibcymin, bcymin, ibcymax, bcymax, &
            p% ilinx, p% iliny, ierr )
         if (ierr /= 0) then
            if (dbg) write(*,*) 'interp_mkbicub_db failed for ionization interpolant'
            return
         end if
         p% have_interpolation_info = .true.
      end subroutine create_interpolants

      real(dp) function charge_of_Fe56_in_He4(log_ne, logT, ierr)
         use interp_2d_lib_db
         real(dp), intent(in) :: log_ne ! ne=avo*rho*free_e
         real(dp), intent(in) :: logT
         integer, intent(out) :: ierr

         integer :: ict(6) ! code specifying output desired
         real(dp) :: fval(6) ! output data
         type (Ionization_Info), pointer :: p
         
         ierr = 0
         charge_of_Fe56_in_He4 = 0
         
         if (.not. table_is_initialized) then
!$omp critical (ionization_table)
            if (.not. table_is_initialized) call do_load(ierr)
!$omp end critical (ionization_table)
            if (ierr /= 0) return
         endif
         
         ict = 0; ict(1) = 1 ! just the result; no partials
         p => fe_he_ptr
         call interp_evbicub_db( &
            log_ne, logT, p% log_ne, p% num_log_ne, p% logT, p% num_logT, &
            p% ilinx, p% iliny, p% f1, p% num_log_ne, ict, fval, ierr)
         
         charge_of_Fe56_in_He4 = fval(1)
         
      end function charge_of_Fe56_in_He4

      subroutine chi_info(a1, z1, T, log_T, rho, log_rho, chi, c0, c1, c2)
         real(dp), intent(in) :: a1, z1, T, log_T, rho, log_rho
         real(dp), intent(out) :: chi, c0, c1, c2
         chi = 1.987d-4*T*(-8.392d0 - log_rho + 1.5d0*log_T - log10(z1/a1)) ! eqn 20
         ! coef's used in eqn 21
         c0 = 1.085d-4*rho*T/a1
         c1 = 1.617d4*sqrt(rho*(z1*z1 + z1)/(T*a1))
         c2 = 29.38d0*z1*pow(rho/a1,one_third)
         ! c2 had a typo in eqn 21, now corrected to match Dupuis et al. (1992) eqn 3
      end subroutine chi_info
      
      real(dp) function chi_effective(chi, c0, c1, c2, z1, z2)
         real(dp), intent(in) :: chi, c0, c1, c2, z1, z2
         chi_effective = chi + c0/(z2*z2*z2) + &
            min(c1*z2, c2*(pow(z2/z1,two_thirds) + 0.6d0))
      end function chi_effective
      
      real(dp) function get_typical_charge(cid, a1, z1, T, log_T, rho, log_rho)
         use mod_tables
         use chem_def
         integer, intent(in) :: cid
         real(dp), intent(in) :: a1, z1, T, log_T, rho, log_rho      
         real(dp) :: chi, c0, c1, c2, z2, chi_eff
         integer :: i, izmax
         include 'formats'
         if (.not. ionization_tables_okay) then
            call set_ionization_potentials
            ionization_tables_okay = .true.
         end if
         izmax = int(chem_isos% Z(cid))
         get_typical_charge = dble(izmax)
         if (izmax > 30) return
         call chi_info(a1, z1, T, log_T, rho, log_rho, chi, c0, c1, c2)
         do i=1, izmax-1
            z2 = dble(i)
            chi_eff = chi_effective(chi, c0, c1, c2, z1, z2+1)
            if (chi_eff < ip(izmax,i+1)) then
               get_typical_charge = z2
               return
            end if
         end do
      end function get_typical_charge
      
      



      end module mod_ionization

