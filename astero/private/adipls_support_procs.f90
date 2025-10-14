! ***********************************************************************
!
!   Copyright (C) 2013  The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! ***********************************************************************

      ! routines that are called by adipls
      ! uses adipls_support, so must compile after that

      subroutine spcout_adi(x,y,aa,data,nn,iy,iaa,ispcpr)
         ! must set ispcpr > 0 to get this called
         use astero_def, only: store_new_oscillation_results, &
            el, order, em, cyclic_freq, inertia, num_results
         use adipls_support, only: adipls_mode_info
         use const_def, only: dp, pi4
         use utils_lib, only: mesa_error

         implicit none

         integer :: nn, iy, iaa, ispcpr
         real(dp) :: x(1:nn), y(1:iy,1:nn), aa(1:iaa,1:nn), data(8)

         !  common for storage of model parameters
         !  degree, order, cyclic frequency (microHz), inertia
         common/cobs_param/ icobs_st, nobs_st, obs_st
         real(dp) :: csummm(50)
         common/csumma/ csummm

         integer :: icobs_st, nobs_st
         real(dp) :: obs_st(10,100000)  ! huge 2nd dimension to satisfy bounds checking

         integer :: ierr, new_el, new_order, new_em, n
         real(dp) :: new_inertia, new_cyclic_freq

         include 'formats'

         new_el = int(obs_st(1,nobs_st) + 0.5_dp)
         new_order = int(obs_st(2,nobs_st) + 0.5_dp)
         new_em = csummm(38)
         new_inertia = obs_st(4,nobs_st)*pi4
         new_cyclic_freq = obs_st(3,nobs_st)

         call store_new_oscillation_results( &
            new_el, new_order, new_em, new_inertia, new_cyclic_freq, 0._dp, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)

         n = num_results
         call adipls_mode_info( &
            el(n), order(n), em(n), cyclic_freq(n), inertia(n), &
            x, y, aa, data, nn, iy, iaa, ispcpr)

      end subroutine spcout_adi


      subroutine modmod(x,aa,data,nn,ivarmd,iaa,imdmod)
         use const_def, only: dp
         integer :: nn, ivarmd, iaa, imdmod
         real(dp) :: x(nn), aa(iaa,nn), data(8)
      end subroutine modmod


      subroutine resdif
         return
      end subroutine resdif
