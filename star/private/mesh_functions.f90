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


      module mesh_functions

      use star_private_def
      use const_def
      use num_lib
      use utils_lib
      use chem_def

      implicit none

      private
      public :: num_mesh_functions, set_mesh_function_data, max_allowed_gvals

      integer, parameter :: max_allowed_gvals = 50

      contains


      integer function get_net_iso(s, species)
         use chem_lib, only: chem_get_iso_id
         type (star_info), pointer :: s
         character (len=*) :: species
         integer :: j
         get_net_iso = 0
         if (len_trim(species) == 0) return
         j = chem_get_iso_id(species)
         if (j <= 0) then
            write(*,*) 'unknown species name for mesh function: ' // trim(species)
            write(*,*) 'len_trim(species)', len_trim(species)
            return
         end if
         get_net_iso = s% net_iso(j) ! 0 if species not in current net
      end function get_net_iso


      logical function do_mass_function(s,species,weight,j)
         type (star_info), pointer :: s
         character (len=iso_name_length) :: species
         real(dp), intent(in) :: weight
         integer, intent(out) :: j
         j = 0
         if (weight > 0) j = get_net_iso(s, species)
         do_mass_function = (j > 0)
      end function do_mass_function


      integer function num_mesh_functions(s)
         type (star_info), pointer :: s
         integer :: i, j, k
         i = 0
         if (s% use_other_mesh_functions) &
            call s% how_many_other_mesh_fcns(s% id, i)
         if (s% E_function_weight > 0) i=i+1
         if (s% P_function_weight > 0) i=i+1
         if (s% T_function1_weight > 0) i=i+1
         if (s% T_function2_weight > 0) i=i+1
         if (s% R_function_weight > 0) i=i+1
         if (s% R_function2_weight > 0) i=i+1
         if (s% R_function3_weight > 0) i=i+1
         if (s% M_function_weight > 0) i=i+1
         if (s% gradT_function_weight > 0) i=i+1
         if (s% log_tau_function_weight > 0) i=i+1
         if (s% log_kap_function_weight > 0) i=i+1
         if (s% gam_function_weight > 0) i=i+1
         if (s% omega_function_weight > 0 .and. s% rotation_flag) i=i+1
         do k=1,num_xa_function
            if (do_mass_function( &
                  s, s% xa_function_species(k), s% xa_function_weight(k), j)) &
               i=i+1
         end do
         num_mesh_functions = i
      end function num_mesh_functions


      subroutine set_mesh_function_data( &
            s, nfcns, names, gval_is_xa_function, gval_is_logT_function, vals1, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nfcns
         character (len=32), intent(out) :: names(max_allowed_gvals)
         logical, intent(out), dimension(max_allowed_gvals) :: &
            gval_is_xa_function, gval_is_logT_function
         real(dp), pointer :: vals1(:) ! =(nz, nfcns)
         integer, intent(out) :: ierr

         integer :: i, nz, j, k, i_other
         logical, parameter :: dbg = .false.
         real(dp), dimension(:,:), pointer :: vals

         ierr = 0
         nz = s% nz

         vals(1:nz,1:nfcns) => vals1(1:nz*nfcns)
         gval_is_xa_function = .false.
         gval_is_logT_function = .false.

         i_other = 0
         if (s% use_other_mesh_functions) then
            call s% how_many_other_mesh_fcns(s% id, i_other)
            if (i_other > 0) then
               if (i_other > nfcns) then
                  ierr = -1
                  return
               end if
               call s% other_mesh_fcn_data( &
                  s% id, i_other, names, gval_is_xa_function, vals1, ierr)
               if (ierr /= 0) return
            end if
         end if

         i = i_other
         if (s% E_function_weight > 0) then
            i = i+1; names(i) = 'E_function'
         end if
         if (s% P_function_weight > 0) then
            i = i+1; names(i) = 'P_function'
         end if
         if (s% T_function1_weight > 0) then
            i = i+1; names(i) = 'T_function1'
            gval_is_logT_function(i) = .true.
         end if
         if (s% T_function2_weight > 0) then
            i = i+1; names(i) = 'T_function2'
            gval_is_logT_function(i) = .true.
         end if
         if (s% R_function_weight > 0) then
            i = i+1; names(i) = 'R_function'
         end if
         if (s% R_function2_weight > 0) then
            i = i+1; names(i) = 'R_function2'
         end if
         if (s% R_function3_weight > 0) then
            i = i+1; names(i) = 'R_function3'
         end if
         if (s% M_function_weight > 0) then
            i = i+1; names(i) = 'M_function'
         end if
         if (s% gradT_function_weight > 0) then
            i = i+1; names(i) = 'gradT_function'
         end if
         if (s% log_tau_function_weight > 0) then
            i = i+1; names(i) = 'log_tau_function'
         end if
         if (s% log_kap_function_weight > 0) then
            i = i+1; names(i) = 'log_kap_function'
         end if
         if (s% gam_function_weight > 0) then
            i = i+1; names(i) = 'gam_function'
         end if
         if (s% omega_function_weight > 0 .and. s% rotation_flag) then
            i = i+1; names(i) = 'omega_function'
         end if
         do k=1,num_xa_function
            if (do_mass_function(s, s% xa_function_species(k), s% xa_function_weight(k), j)) then
               i = i+1; names(i) = trim(s% xa_function_species(k))
               gval_is_xa_function(i) = .true.
               !write(names(i),'(a,i1)') 'xa_function_', k
            end if
         end do
         if (i /= nfcns) then
            write(*,*) 'error in set_mesh_function_names: incorrect nfcns'
            ierr = -1
         end if

         do i=i_other+1,nfcns

            if (ierr /= 0) cycle
            if (dbg) write(*,*) trim(names(i))

            if (names(i) == 'E_function') then
               do k=1,nz
                  vals(k,i) = s% E_function_weight * max(s% E_function_param, s% lnE(k)/ln10)
               end do

            else if (names(i) == 'P_function') then
               do k=1,nz
                  vals(k,i) = s% P_function_weight * s% lnPeos(k)/ln10
               end do

            else if (names(i) == 'T_function1') then
               do k=1,nz
                  vals(k,i) = s% T_function1_weight*s% lnT(k)/ln10
               end do

            else if (names(i) == 'T_function2') then
               do k=1,nz
                  vals(k,i) = &
                     s% T_function2_weight*log10(s% T(k) / (s% T(k) + s% T_function2_param))
               end do

            else if (names(i) == 'R_function') then
               do k=1,nz
                  vals(k,i) = &
                     s% R_function_weight*log10(1 + (s% r(k)/Rsun)/s% R_function_param)
               end do

            else if (names(i) == 'R_function2') then
               do k=1,nz
                  vals(k,i) = &
                     s% R_function2_weight * &
                     min(s% R_function2_param1, max(s% R_function2_param2,s% r(k)/s% r(1)))
               end do

            else if (names(i) == 'R_function3') then
               do k=1,nz
                  vals(k,i) = s% R_function3_weight*s% r(k)/s% r(1)
               end do

            else if (names(i) == 'M_function') then
               do k=1,nz
                  vals(k,i) = &
                  s% M_function_weight*log10(1 + (s% xmstar*s% q(k)/Msun)/s% M_function_param)
               end do

            else if (names(i) == 'gradT_function') then
               do k=1,nz
                  vals(k,i) = s% gradT_function_weight*s% gradT(k)
               end do

            else if (names(i) == 'log_tau_function') then
               do k=1,nz
                  vals(k,i) = s% log_tau_function_weight*log10(s% tau(k))
               end do

            else if (names(i) == 'log_kap_function') then
               do k=1,nz
                  vals(k,i) = s% log_kap_function_weight*log10(s% opacity(k))
               end do

            else if (names(i) == 'gam_function') then
               do k=1,nz
                  vals(k,i) = s% gam_function_weight* &
                     tanh((s% gam(k) - s% gam_function_param1)/s% gam_function_param2)
               end do

            else if (names(i) == 'omega_function') then
               do k=1,nz
                  vals(k,i) = s% omega_function_weight*log10(max(1d-99,abs(s% omega(k))))
               end do

            else
               do k=1,num_xa_function
                  call do1_xa_function(k,i)
               end do

            end if

         end do

         contains

         subroutine do1_xa_function(k,i)
            integer, intent(in) :: k,i
            real(dp) :: weight, param
            integer :: j, m
            if (len_trim(s% xa_function_species(k))==0) return
            if (trim(names(i)) /= trim(s% xa_function_species(k))) return
            j = get_net_iso(s, s% xa_function_species(k))
            if (j <= 0) then
               ierr = -1
            else
               weight = s% xa_function_weight(k)
               param = s% xa_function_param(k)
               do m=1,s% nz
                  vals(m,i) = weight*log10(s% xa(j,m) + param)
               end do
            end if
         end subroutine do1_xa_function


      end subroutine set_mesh_function_data


      end module mesh_functions
