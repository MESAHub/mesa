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

module pulse_saio
   
   ! Uses
   
   use star_private_def
   use const_def
   use utils_lib
   use atm_def
   use atm_support
   
   use pulse_utils
   
   ! No implicit typing
   
   implicit none
   
   ! Access specifiers
   
   private
   
   public :: get_saio_data
   public :: write_saio_data

contains
   
   subroutine get_saio_data (id, &
      keep_surface_point, add_atmosphere, global_data, point_data, ierr)
      
      integer, intent(in) :: id
      logical, intent(in) :: keep_surface_point
      logical, intent(in) :: add_atmosphere
      real(dp), allocatable, intent(out) :: global_data(:)
      real(dp), allocatable, intent(out) :: point_data(:, :)
      integer, intent(out) :: ierr
      
      type(star_info), pointer :: s
      integer, allocatable :: k_a(:)
      integer, allocatable :: k_b(:)
      integer :: n_sg
      integer :: n_atm
      integer :: nn_atm
      integer :: n_env
      integer :: nn_env
      integer :: nn
      real(dp) :: r_outer
      real(dp) :: m_outer
      integer :: j
      integer :: k
      integer :: sg
      
      ! Get SAIO data
      
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) then
         write(*, *) 'bad star id for get_saio_data'
         return
      end if
      
      ! Set up segment indices
      
      call set_segment_indices(s, k_a, k_b, .FALSE.)
      
      n_sg = SIZE(k_a)
      
      ! Determine data dimensiones
      
      if (add_atmosphere) then
         call build_atm(s, s%L(1), s%r(1), s%Teff, s%m_grav(1), s%cgrav(1), ierr)
         if (ierr /= 0) then
            write(*, *) 'failed in build_atm'
            return
         end if
         n_atm = s%atm_structure_num_pts
         nn_atm = n_atm - 1
      else
         n_atm = 0
         nn_atm = 0
      end if
      
      n_env = s%nz
      
      if (keep_surface_point) then
         nn_env = n_env + n_sg - 1
      else
         nn_env = n_env - 1 + n_sg - 1
      endif
      
      nn = nn_env + nn_atm
      
      ! Store global data
      
      allocate(global_data(4))
      
      r_outer = Rsun * s%photosphere_r
      m_outer = s%m_grav(1)
      
      global_data(1) = m_outer
      global_data(2) = log10(s%L(1) / Lsun)
      global_data(3) = log10(r_outer / Rsun)
      global_data(4) = s%star_age
      
      ! Store point data
      
      allocate(point_data(20, nn))
      
      j = 1
      
      ! Atmosphere (we skip the point at the base of the atm to
      ! smooth the transition)
      
      atm_loop : do k = 1, n_atm - 1
         call store_saio_data_atm(j, k, k_a(1), k_b(1))
         j = j + 1
      end do atm_loop
      
      ! Envelope
      
      sg = 1
      
      env_loop : do k = 1, n_env
         
         if (k == 1 .AND. .NOT. keep_surface_point) cycle env_loop
         
         call store_saio_data_env(j, k, k_a(sg), k_b(sg))
         j = j + 1
         
         if (k == k_b(sg) + 1) then
            
            sg = sg + 1
            
            call store_saio_data_env(j, k, k_a(sg), k_b(sg))
            j = j + 1
         
         endif
      
      end do env_loop
      
      ! Reverse the data ordering (SAIO format is center-to-surface)
      
      point_data = point_data(:, nn:1:-1)
      
      ! Tidy up
      
      if (ASSOCIATED(s%atm_structure)) then
         deallocate(s%atm_structure)
      end if
      
      ! Finish
      
      return
   
   contains
      
      subroutine store_saio_data_atm (j, k, k_a, k_b)
         
         integer, intent(in) :: j
         integer, intent(in) :: k
         integer, intent(in) :: k_a
         integer, intent(in) :: k_b
         
         ! Store data associated with atmosphere point k into the
         ! point_data array at position j
         
         associate (&
            r => point_data(1, j), &
            m => point_data(2, j), &
            Lrad => point_data(3, j), &
            T => point_data(4, j), &
            rho => point_data(5, j), &
            P => point_data(6, j), &
            eps => point_data(7, j), &
            kap => point_data(8, j), &
            c_V => point_data(9, j), &
            chi_rho => point_data(10, j), &
            chi_T => point_data(11, j), &
            eps_rho => point_data(12, j), &
            eps_T => point_data(13, j), &
            kap_rho => point_data(14, j), &
            kap_T => point_data(15, j), &
            nabla => point_data(16, j), &
            nabla_ad => point_data(17, j), &
            X => point_data(18, j), &
            Y => point_data(19, j), &
            L => point_data(20, j))
            
            r = s% r(1) + s% atm_structure(atm_delta_r, k)
            m = s% m_grav(1)
            T = exp(s% atm_structure(atm_lnT, k))
            rho = exp(s% atm_structure(atm_lnd, k))
            P = exp(s% atm_structure(atm_lnP, k))
            eps = 0d0
            kap = s% atm_structure(atm_kap, k)
            c_V = s% atm_structure(atm_cv, k)
            chi_rho = s% atm_structure(atm_chiRho, k)
            chi_T = s% atm_structure(atm_chiT, k)
            eps_rho = 0d0
            eps_T = 0d0
            kap_rho = s% atm_structure(atm_dlnkap_dlnd, k)
            kap_T = s% atm_structure(atm_dlnkap_dlnT, k)
            nabla = s% atm_structure(atm_gradT, k)
            nabla_ad = s% atm_structure(atm_grada, k)
            X = eval_face(s%dq, s%X, 1, k_a, k_b, v_lo = 0d0, v_hi = 1d0)
            Y = eval_face(s%dq, s%Y, 1, k_a, k_b, v_lo = 0d0, v_hi = 1d0)
            L = s% L(1)
            
            if (s%atm_structure(atm_tau, k) >= 2d0 / 3d0) then
               Lrad = (16 * pi * clight * crad * s%cgrav(k) * m * T * T * T * T * nabla) / (3 * P * kap)
            else
               Lrad = L
            end if
         
         end associate
         
         ! Finish
         
         return
      
      end subroutine store_saio_data_atm
      
      !****
      
      subroutine store_saio_data_env (j, k, k_a, k_b)
         
         integer, intent(in) :: j
         integer, intent(in) :: k
         integer, intent(in) :: k_a
         integer, intent(in) :: k_b
         
         ! Store data for envelope face k into the point_data array at
         ! position j
         
         associate (&
            r => point_data(1, j), &
            m => point_data(2, j), &
            Lrad => point_data(3, j), &
            T => point_data(4, j), &
            rho => point_data(5, j), &
            P => point_data(6, j), &
            eps => point_data(7, j), &
            kap => point_data(8, j), &
            c_V => point_data(9, j), &
            chi_rho => point_data(10, j), &
            chi_T => point_data(11, j), &
            eps_rho => point_data(12, j), &
            eps_T => point_data(13, j), &
            kap_rho => point_data(14, j), &
            kap_T => point_data(15, j), &
            nabla => point_data(16, j), &
            nabla_ad => point_data(17, j), &
            X => point_data(18, j), &
            Y => point_data(19, j), &
            L => point_data(20, j))
            
            r = s%r(k)
            m = s%m_grav(k)
            T = eval_face(s%dq, s%T, k, 1, s%nz)
            if (s%interpolate_rho_for_pulse_data) then
               rho = eval_face(s%dq, s%rho, k, k_a, k_b)
            else
               rho = eval_face_rho(s, k, k_a, k_b)
            end if
            P = eval_face(s%dq, s%Peos, k, 1, s%nz)
            eps = eval_face(s%dq, s%eps_nuc, k, k_a, k_b) + eval_face(s%dq, s%eps_grav_ad%val, k, k_a, k_b)
            c_V = eval_face(s%dq, s%Cv, k, k_a, k_b)
            chi_rho = eval_face(s%dq, s%chiRho, k, k_a, k_b)
            chi_T = eval_face(s%dq, s%chiT, k, k_a, k_b)
            eps_rho = eval_face(s%dq, s%d_epsnuc_dlnd, k, k_a, k_b)
            eps_T = eval_face(s%dq, s%d_epsnuc_dlnT, k, k_a, k_b)
            kap_rho = eval_face(s%dq, s%d_opacity_dlnd, k, k_a, k_b) / kap
            kap_T = eval_face(s%dq, s%d_opacity_dlnT, k, k_a, k_b) / kap
            nabla = s%gradT(k) ! Not quite right; gradT can be discontinuous
            nabla_ad = eval_face(s%dq, s%grada, k, k_a, k_b)
            X = eval_face(s%dq, s%X, k, k_a, k_b, v_lo = 0d0, v_hi = 1d0)
            Y = eval_face(s%dq, s%Y, k, k_a, k_b, v_lo = 0d0, v_hi = 1d0)
            L = s%L(k)
            
            if (s%tau(k) >= 2d0 / 3d0) then
               Lrad = (16 * pi * clight * crad * s%cgrav(k) * m * T * T * T * T * nabla) / (3 * P * kap)
            else
               Lrad = L
            end if
         
         end associate
         
         ! Finish
         
         return
      
      end subroutine store_saio_data_env
   
   end subroutine get_saio_data
   
   !****
   
   subroutine write_saio_data (id, filename, global_data, point_data, ierr)
      
      integer, intent(in) :: id
      character(*), intent(in) :: filename
      real(dp), intent(in) :: global_data(:)
      real(dp), intent(in) :: point_data(:, :)
      integer, intent(out) :: ierr
      
      type(star_info), pointer :: s
      integer :: iounit
      integer :: nn
      integer :: j
      
      ! Write SAIO data to file
      
      call get_star_ptr(id, s, ierr)
      if (ierr /= 0) then
         write(*, *) 'bad star id for get_fgong_info'
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
      
      write(iounit, 100) nn, global_data
      100 format(I6, 99(1PE16.9))
      
      do j = 1, nn
         write(iounit, 110) point_data(:, j)
         110    format(5(1PE16.9))
      end do
      
      ! Close the file
      
      close(iounit)
      
      ! Finish
      
      return
   
   end subroutine write_saio_data

end module pulse_saio
