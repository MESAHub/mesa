! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
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
!
! ***********************************************************************

      module eos_initialize

      implicit none

      contains
            
      
      subroutine Init_eos( &
           eosDT_cache_dir_in, eosPT_cache_dir_in, &
           use_cache, ierr)
         use eos_def
         use helm_alloc
         use utils_lib, only : mkdir
         use const_def, only: mesa_data_dir, mesa_caches_dir, mesa_temp_caches_dir
         character(*), intent(IN) :: &
            eosDT_cache_dir_in, eosPT_cache_dir_in
         logical, intent(in) :: use_cache
         integer, intent(OUT) :: ierr ! 0 means AOK.
         !integer, parameter :: imax = 261, jmax = 101  
            ! dimensions of small version of helm table
         !integer, parameter :: imax = 1081, jmax = 401
            ! dimensions of medium version of helm table; 40 points per decade
         integer, parameter :: imax = 2701, jmax = 1001
            ! dimensions of large version of helm table; 100 points per decade
         ! helm table lives in eosDT_data
         character (len=256) :: eosDT_data_dir, eosPT_data_dir
         ierr = 0
         if (eos_root_is_initialized) return
         use_cache_for_eos = use_cache
         eosDT_data_dir = trim(mesa_data_dir) // '/eosDT_data'
         eosPT_data_dir = trim(mesa_data_dir) // '/eosPT_data'
         if (use_cache_for_eos) then
            if (len_trim(eosDT_cache_dir_in) > 0) then
               eosDT_cache_dir = eosDT_cache_dir_in
            else if (len_trim(mesa_caches_dir) > 0) then
               eosDT_cache_dir = trim(mesa_caches_dir) // '/eosDT_cache'
            else
               eosDT_cache_dir = trim(eosDT_data_dir) // '/cache'
            end if
            call mkdir(eosDT_cache_dir)
            eosDT_temp_cache_dir = trim(mesa_temp_caches_dir) // '/eosDT_cache'
            if(use_mesa_temp_cache) call mkdir(eosDT_temp_cache_dir)
            if (len_trim(eosPT_cache_dir_in) > 0) then
               eosPT_cache_dir = eosPT_cache_dir_in
            else if (len_trim(mesa_caches_dir) > 0) then
               eosPT_cache_dir = trim(mesa_caches_dir) // '/eosPT_cache'
            else
               eosPT_cache_dir = trim(eosPT_data_dir) // '/cache'
            end if
            call mkdir(eosPT_cache_dir)
            eosPT_temp_cache_dir = trim(mesa_temp_caches_dir) // '/eosPT_cache'
            if(use_mesa_temp_cache) call mkdir(eosPT_temp_cache_dir)
         end if
         
         call alloc_helm_table(eos_ht, imax, jmax, ierr)
         if (ierr /= 0) return
         
         call read_helm_table(eos_ht, &
            eosDT_data_dir, eosDT_cache_dir, eosDT_temp_cache_dir, use_cache_for_eos, ierr)
         if (ierr /= 0) return

         call eos_def_init
         ! replace defaults from eos_def_init by argument
         
         call theta_e_interp_init(ierr)
         if (ierr /= 0) return
         
         eos_root_is_initialized = .true.
      
      end subroutine Init_eos
      
      
      subroutine theta_e_interp_init(ierr)
         use const_def, only: dp
         use interp_1d_def, only : pm_work_size
         use interp_1d_lib, only : interp_pm
         use eos_def
         integer, intent(out) :: ierr
         real(dp), parameter :: eta0 = -100
         real(dp), parameter :: deta0 = 10
         
         real(dp), parameter :: eta1 = -20
         real(dp), parameter :: deta1 = 1
         
         real(dp), parameter :: eta2 = -5
         real(dp), parameter :: deta2 = 0.1d0
         
         real(dp), parameter :: eta3 = 25
         real(dp), parameter :: deta3 = 1

         real(dp), parameter :: eta4 = 100
         real(dp), parameter :: deta4 = 10

         real(dp), parameter :: eta5 = 200
         real(dp), parameter :: deta5 = 100

         real(dp), parameter :: eta6 = 2000
         
         integer :: nx, pass, i, j
         real(dp) :: eta, theta_e
         real(dp), pointer :: work(:)
         
         ierr = 0
      
         do pass = 1, 2
            i = 0
            eta = eta0
            do while (eta < eta1)
               call do_eta
               eta = eta + deta0
            end do
            do while (eta < eta2)
               call do_eta
               eta = eta + deta1
            end do
            do while (eta < eta3)
               call do_eta
               eta = eta + deta2
            end do
            do while (eta < eta4)
               call do_eta
               eta = eta + deta3
            end do
            do while (eta < eta5)
               call do_eta
               eta = eta + deta4
            end do
            do while (eta < eta6)
               call do_eta
               eta = eta + deta5
            end do
            if (pass == 1) then
               nx = i
               theta_e_nx = nx
               allocate(f_theta_e1(4*nx), x_theta_e(nx), work(nx*pm_work_size), stat=ierr)
               if (ierr /= 0) return
               f_theta_e(1:4,1:nx) => f_theta_e1(1:4*nx)
            end if
         end do
         
         call interp_pm(x_theta_e, nx, f_theta_e1, pm_work_size, work, 'theta_e_interp_init', ierr)
         deallocate(work)
         
         if (ierr /= 0) then
            deallocate(f_theta_e1, x_theta_e)
            theta_e_nx = 0
         end if
         
         contains
         
         subroutine do_eta
            i = i+1
            if (pass == 1) return
            x_theta_e(i) = eta
            f_theta_e(1, i) = eval(eta)
         end subroutine do_eta
         
         real(dp) function eval(eta)
            use gauss_fermi, only: dfermi
           
            real(dp), intent(in) :: eta
            real(dp) :: theta, d_theta_e_deta, &
               imh, dimh_deta, dimh_dtheta,  &
               iph, diph_deta, diph_dtheta
            theta = 0
            call dfermi(-0.5d0, eta, theta, imh, dimh_deta, dimh_dtheta)
            call dfermi(0.5d0, eta, theta, iph, diph_deta, diph_dtheta)
            ! the factor of 0.5 comes from slightly different definitions of fd integrals
            eval = 0.5d0*imh/iph            
         end function eval      
      
      end subroutine theta_e_interp_init
      
      
      end module eos_initialize
