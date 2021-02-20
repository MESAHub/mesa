! Module  : hdf5io
! Purpose : definition of hdf5io_t type, supporting HDF5 input/output
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


























module hdf5io_m

   ! Uses

   use kinds_m
   use cinter_m

   use hdf5

   use ISO_FORTRAN_ENV
   use ISO_C_BINDING

   ! No implicit typing

   implicit none

   ! Parameters

   integer, parameter :: CREATE_FILE = 1
   integer, parameter :: OPEN_FILE = 2

   integer, parameter :: TYPE_LEN = 31

   ! Derived-type definitions

   type hdf5io_t
      private
      integer(HID_T)   :: file_id = 0
      integer(HID_T)   :: group_id = 0
      integer          :: comp_level = 0
      integer, pointer :: ref_count => null()
   contains
      private
      procedure, public :: final => final_
         procedure, public :: group_exists
         procedure, public :: attr_exists
         procedure, public :: dset_exists
         procedure, public :: group_names
         procedure, public :: dset_names
         procedure, public :: attr_shape
         procedure, public :: dset_shape
                  procedure :: read_attr_is_0_
                  procedure :: read_attr_id_0_
                  procedure :: read_attr_rs_0_
                  procedure :: read_attr_rd_0_
                  procedure :: read_attr_cs_0_
                  procedure :: read_attr_cd_0_
                  procedure :: read_attr_l_0_
                  procedure :: read_attr_s_0_
                  procedure :: read_attr_is_1_
                  procedure :: read_attr_id_1_
                  procedure :: read_attr_rs_1_
                  procedure :: read_attr_rd_1_
                  procedure :: read_attr_cs_1_
                  procedure :: read_attr_cd_1_
                  procedure :: read_attr_l_1_
                  procedure :: read_attr_s_1_
                  procedure :: read_attr_is_2_
                  procedure :: read_attr_id_2_
                  procedure :: read_attr_rs_2_
                  procedure :: read_attr_rd_2_
                  procedure :: read_attr_cs_2_
                  procedure :: read_attr_cd_2_
                  procedure :: read_attr_l_2_
                  procedure :: read_attr_s_2_
                  procedure :: read_attr_is_3_
                  procedure :: read_attr_id_3_
                  procedure :: read_attr_rs_3_
                  procedure :: read_attr_rd_3_
                  procedure :: read_attr_cs_3_
                  procedure :: read_attr_cd_3_
                  procedure :: read_attr_l_3_
                  procedure :: read_attr_s_3_
                  procedure :: read_attr_is_4_
                  procedure :: read_attr_id_4_
                  procedure :: read_attr_rs_4_
                  procedure :: read_attr_rd_4_
                  procedure :: read_attr_cs_4_
                  procedure :: read_attr_cd_4_
                  procedure :: read_attr_l_4_
                  procedure :: read_attr_s_4_
                  procedure :: read_attr_is_5_
                  procedure :: read_attr_id_5_
                  procedure :: read_attr_rs_5_
                  procedure :: read_attr_rd_5_
                  procedure :: read_attr_cs_5_
                  procedure :: read_attr_cd_5_
                  procedure :: read_attr_l_5_
                  procedure :: read_attr_s_5_
                  procedure :: read_attr_is_6_
                  procedure :: read_attr_id_6_
                  procedure :: read_attr_rs_6_
                  procedure :: read_attr_rd_6_
                  procedure :: read_attr_cs_6_
                  procedure :: read_attr_cd_6_
                  procedure :: read_attr_l_6_
                  procedure :: read_attr_s_6_
                  procedure :: read_attr_is_7_
                  procedure :: read_attr_id_7_
                  procedure :: read_attr_rs_7_
                  procedure :: read_attr_rd_7_
                  procedure :: read_attr_cs_7_
                  procedure :: read_attr_cd_7_
                  procedure :: read_attr_l_7_
                  procedure :: read_attr_s_7_
            generic, public :: read_attr => read_attr_is_0_,read_attr_id_0_,read_attr_rs_0_,read_attr_rd_0_,read_attr_cs_0_,read_at&
                &tr_cd_0_,read_attr_l_0_,read_attr_s_0_,read_attr_is_1_,read_attr_id_1_,read_attr_rs_1_,read_attr_rd_1_,read_attr_c&
                &s_1_,read_attr_cd_1_,read_attr_l_1_,read_attr_s_1_,read_attr_is_2_,read_attr_id_2_,read_attr_rs_2_,read_attr_rd_2_&
                &,read_attr_cs_2_,read_attr_cd_2_,read_attr_l_2_,read_attr_s_2_,read_attr_is_3_,read_attr_id_3_,read_attr_rs_3_,rea&
                &d_attr_rd_3_,read_attr_cs_3_,read_attr_cd_3_,read_attr_l_3_,read_attr_s_3_,read_attr_is_4_,read_attr_id_4_,read_at&
                &tr_rs_4_,read_attr_rd_4_,read_attr_cs_4_,read_attr_cd_4_,read_attr_l_4_,read_attr_s_4_,read_attr_is_5_,read_attr_i&
                &d_5_,read_attr_rs_5_,read_attr_rd_5_,read_attr_cs_5_,read_attr_cd_5_,read_attr_l_5_,read_attr_s_5_,read_attr_is_6_&
                &,read_attr_id_6_,read_attr_rs_6_,read_attr_rd_6_,read_attr_cs_6_,read_attr_cd_6_,read_attr_l_6_,read_attr_s_6_,rea&
                &d_attr_is_7_,read_attr_id_7_,read_attr_rs_7_,read_attr_rd_7_,read_attr_cs_7_,read_attr_cd_7_,read_attr_l_7_,read_a&
                &ttr_s_7_
                  procedure :: read_dset_is_0_
                  procedure :: read_dset_id_0_
                  procedure :: read_dset_rs_0_
                  procedure :: read_dset_rd_0_
                  procedure :: read_dset_cs_0_
                  procedure :: read_dset_cd_0_
                  procedure :: read_dset_l_0_
                  procedure :: read_dset_s_0_
                  procedure :: read_dset_is_1_
                  procedure :: read_dset_id_1_
                  procedure :: read_dset_rs_1_
                  procedure :: read_dset_rd_1_
                  procedure :: read_dset_cs_1_
                  procedure :: read_dset_cd_1_
                  procedure :: read_dset_l_1_
                  procedure :: read_dset_s_1_
                  procedure :: read_dset_is_2_
                  procedure :: read_dset_id_2_
                  procedure :: read_dset_rs_2_
                  procedure :: read_dset_rd_2_
                  procedure :: read_dset_cs_2_
                  procedure :: read_dset_cd_2_
                  procedure :: read_dset_l_2_
                  procedure :: read_dset_s_2_
                  procedure :: read_dset_is_3_
                  procedure :: read_dset_id_3_
                  procedure :: read_dset_rs_3_
                  procedure :: read_dset_rd_3_
                  procedure :: read_dset_cs_3_
                  procedure :: read_dset_cd_3_
                  procedure :: read_dset_l_3_
                  procedure :: read_dset_s_3_
                  procedure :: read_dset_is_4_
                  procedure :: read_dset_id_4_
                  procedure :: read_dset_rs_4_
                  procedure :: read_dset_rd_4_
                  procedure :: read_dset_cs_4_
                  procedure :: read_dset_cd_4_
                  procedure :: read_dset_l_4_
                  procedure :: read_dset_s_4_
                  procedure :: read_dset_is_5_
                  procedure :: read_dset_id_5_
                  procedure :: read_dset_rs_5_
                  procedure :: read_dset_rd_5_
                  procedure :: read_dset_cs_5_
                  procedure :: read_dset_cd_5_
                  procedure :: read_dset_l_5_
                  procedure :: read_dset_s_5_
                  procedure :: read_dset_is_6_
                  procedure :: read_dset_id_6_
                  procedure :: read_dset_rs_6_
                  procedure :: read_dset_rd_6_
                  procedure :: read_dset_cs_6_
                  procedure :: read_dset_cd_6_
                  procedure :: read_dset_l_6_
                  procedure :: read_dset_s_6_
                  procedure :: read_dset_is_7_
                  procedure :: read_dset_id_7_
                  procedure :: read_dset_rs_7_
                  procedure :: read_dset_rd_7_
                  procedure :: read_dset_cs_7_
                  procedure :: read_dset_cd_7_
                  procedure :: read_dset_l_7_
                  procedure :: read_dset_s_7_
            generic, public :: read_dset => read_dset_is_0_,read_dset_id_0_,read_dset_rs_0_,read_dset_rd_0_,read_dset_cs_0_,read_ds&
                &et_cd_0_,read_dset_l_0_,read_dset_s_0_,read_dset_is_1_,read_dset_id_1_,read_dset_rs_1_,read_dset_rd_1_,read_dset_c&
                &s_1_,read_dset_cd_1_,read_dset_l_1_,read_dset_s_1_,read_dset_is_2_,read_dset_id_2_,read_dset_rs_2_,read_dset_rd_2_&
                &,read_dset_cs_2_,read_dset_cd_2_,read_dset_l_2_,read_dset_s_2_,read_dset_is_3_,read_dset_id_3_,read_dset_rs_3_,rea&
                &d_dset_rd_3_,read_dset_cs_3_,read_dset_cd_3_,read_dset_l_3_,read_dset_s_3_,read_dset_is_4_,read_dset_id_4_,read_ds&
                &et_rs_4_,read_dset_rd_4_,read_dset_cs_4_,read_dset_cd_4_,read_dset_l_4_,read_dset_s_4_,read_dset_is_5_,read_dset_i&
                &d_5_,read_dset_rs_5_,read_dset_rd_5_,read_dset_cs_5_,read_dset_cd_5_,read_dset_l_5_,read_dset_s_5_,read_dset_is_6_&
                &,read_dset_id_6_,read_dset_rs_6_,read_dset_rd_6_,read_dset_cs_6_,read_dset_cd_6_,read_dset_l_6_,read_dset_s_6_,rea&
                &d_dset_is_7_,read_dset_id_7_,read_dset_rs_7_,read_dset_rd_7_,read_dset_cs_7_,read_dset_cd_7_,read_dset_l_7_,read_d&
                &set_s_7_
                  procedure :: alloc_read_attr_is_0_
                  procedure :: alloc_read_attr_id_0_
                  procedure :: alloc_read_attr_rs_0_
                  procedure :: alloc_read_attr_rd_0_
                  procedure :: alloc_read_attr_cs_0_
                  procedure :: alloc_read_attr_cd_0_
                  procedure :: alloc_read_attr_l_0_
                  procedure :: alloc_read_attr_s_0_
                  procedure :: alloc_read_attr_is_1_
                  procedure :: alloc_read_attr_id_1_
                  procedure :: alloc_read_attr_rs_1_
                  procedure :: alloc_read_attr_rd_1_
                  procedure :: alloc_read_attr_cs_1_
                  procedure :: alloc_read_attr_cd_1_
                  procedure :: alloc_read_attr_l_1_
                  procedure :: alloc_read_attr_s_1_
                  procedure :: alloc_read_attr_is_2_
                  procedure :: alloc_read_attr_id_2_
                  procedure :: alloc_read_attr_rs_2_
                  procedure :: alloc_read_attr_rd_2_
                  procedure :: alloc_read_attr_cs_2_
                  procedure :: alloc_read_attr_cd_2_
                  procedure :: alloc_read_attr_l_2_
                  procedure :: alloc_read_attr_s_2_
                  procedure :: alloc_read_attr_is_3_
                  procedure :: alloc_read_attr_id_3_
                  procedure :: alloc_read_attr_rs_3_
                  procedure :: alloc_read_attr_rd_3_
                  procedure :: alloc_read_attr_cs_3_
                  procedure :: alloc_read_attr_cd_3_
                  procedure :: alloc_read_attr_l_3_
                  procedure :: alloc_read_attr_s_3_
                  procedure :: alloc_read_attr_is_4_
                  procedure :: alloc_read_attr_id_4_
                  procedure :: alloc_read_attr_rs_4_
                  procedure :: alloc_read_attr_rd_4_
                  procedure :: alloc_read_attr_cs_4_
                  procedure :: alloc_read_attr_cd_4_
                  procedure :: alloc_read_attr_l_4_
                  procedure :: alloc_read_attr_s_4_
                  procedure :: alloc_read_attr_is_5_
                  procedure :: alloc_read_attr_id_5_
                  procedure :: alloc_read_attr_rs_5_
                  procedure :: alloc_read_attr_rd_5_
                  procedure :: alloc_read_attr_cs_5_
                  procedure :: alloc_read_attr_cd_5_
                  procedure :: alloc_read_attr_l_5_
                  procedure :: alloc_read_attr_s_5_
                  procedure :: alloc_read_attr_is_6_
                  procedure :: alloc_read_attr_id_6_
                  procedure :: alloc_read_attr_rs_6_
                  procedure :: alloc_read_attr_rd_6_
                  procedure :: alloc_read_attr_cs_6_
                  procedure :: alloc_read_attr_cd_6_
                  procedure :: alloc_read_attr_l_6_
                  procedure :: alloc_read_attr_s_6_
                  procedure :: alloc_read_attr_is_7_
                  procedure :: alloc_read_attr_id_7_
                  procedure :: alloc_read_attr_rs_7_
                  procedure :: alloc_read_attr_rd_7_
                  procedure :: alloc_read_attr_cs_7_
                  procedure :: alloc_read_attr_cd_7_
                  procedure :: alloc_read_attr_l_7_
                  procedure :: alloc_read_attr_s_7_
            generic, public :: alloc_read_attr => alloc_read_attr_is_0_,alloc_read_attr_id_0_,alloc_read_attr_rs_0_,alloc_read_attr&
                &_rd_0_,alloc_read_attr_cs_0_,alloc_read_attr_cd_0_,alloc_read_attr_l_0_,alloc_read_attr_s_0_,alloc_read_attr_is_1_&
                &,alloc_read_attr_id_1_,alloc_read_attr_rs_1_,alloc_read_attr_rd_1_,alloc_read_attr_cs_1_,alloc_read_attr_cd_1_,all&
                &oc_read_attr_l_1_,alloc_read_attr_s_1_,alloc_read_attr_is_2_,alloc_read_attr_id_2_,alloc_read_attr_rs_2_,alloc_rea&
                &d_attr_rd_2_,alloc_read_attr_cs_2_,alloc_read_attr_cd_2_,alloc_read_attr_l_2_,alloc_read_attr_s_2_,alloc_read_attr&
                &_is_3_,alloc_read_attr_id_3_,alloc_read_attr_rs_3_,alloc_read_attr_rd_3_,alloc_read_attr_cs_3_,alloc_read_attr_cd_&
                &3_,alloc_read_attr_l_3_,alloc_read_attr_s_3_,alloc_read_attr_is_4_,alloc_read_attr_id_4_,alloc_read_attr_rs_4_,all&
                &oc_read_attr_rd_4_,alloc_read_attr_cs_4_,alloc_read_attr_cd_4_,alloc_read_attr_l_4_,alloc_read_attr_s_4_,alloc_rea&
                &d_attr_is_5_,alloc_read_attr_id_5_,alloc_read_attr_rs_5_,alloc_read_attr_rd_5_,alloc_read_attr_cs_5_,alloc_read_at&
                &tr_cd_5_,alloc_read_attr_l_5_,alloc_read_attr_s_5_,alloc_read_attr_is_6_,alloc_read_attr_id_6_,alloc_read_attr_rs_&
                &6_,alloc_read_attr_rd_6_,alloc_read_attr_cs_6_,alloc_read_attr_cd_6_,alloc_read_attr_l_6_,alloc_read_attr_s_6_,all&
                &oc_read_attr_is_7_,alloc_read_attr_id_7_,alloc_read_attr_rs_7_,alloc_read_attr_rd_7_,alloc_read_attr_cs_7_,alloc_r&
                &ead_attr_cd_7_,alloc_read_attr_l_7_,alloc_read_attr_s_7_
                  procedure :: alloc_read_dset_is_0_
                  procedure :: alloc_read_dset_id_0_
                  procedure :: alloc_read_dset_rs_0_
                  procedure :: alloc_read_dset_rd_0_
                  procedure :: alloc_read_dset_cs_0_
                  procedure :: alloc_read_dset_cd_0_
                  procedure :: alloc_read_dset_l_0_
                  procedure :: alloc_read_dset_s_0_
                  procedure :: alloc_read_dset_is_1_
                  procedure :: alloc_read_dset_id_1_
                  procedure :: alloc_read_dset_rs_1_
                  procedure :: alloc_read_dset_rd_1_
                  procedure :: alloc_read_dset_cs_1_
                  procedure :: alloc_read_dset_cd_1_
                  procedure :: alloc_read_dset_l_1_
                  procedure :: alloc_read_dset_s_1_
                  procedure :: alloc_read_dset_is_2_
                  procedure :: alloc_read_dset_id_2_
                  procedure :: alloc_read_dset_rs_2_
                  procedure :: alloc_read_dset_rd_2_
                  procedure :: alloc_read_dset_cs_2_
                  procedure :: alloc_read_dset_cd_2_
                  procedure :: alloc_read_dset_l_2_
                  procedure :: alloc_read_dset_s_2_
                  procedure :: alloc_read_dset_is_3_
                  procedure :: alloc_read_dset_id_3_
                  procedure :: alloc_read_dset_rs_3_
                  procedure :: alloc_read_dset_rd_3_
                  procedure :: alloc_read_dset_cs_3_
                  procedure :: alloc_read_dset_cd_3_
                  procedure :: alloc_read_dset_l_3_
                  procedure :: alloc_read_dset_s_3_
                  procedure :: alloc_read_dset_is_4_
                  procedure :: alloc_read_dset_id_4_
                  procedure :: alloc_read_dset_rs_4_
                  procedure :: alloc_read_dset_rd_4_
                  procedure :: alloc_read_dset_cs_4_
                  procedure :: alloc_read_dset_cd_4_
                  procedure :: alloc_read_dset_l_4_
                  procedure :: alloc_read_dset_s_4_
                  procedure :: alloc_read_dset_is_5_
                  procedure :: alloc_read_dset_id_5_
                  procedure :: alloc_read_dset_rs_5_
                  procedure :: alloc_read_dset_rd_5_
                  procedure :: alloc_read_dset_cs_5_
                  procedure :: alloc_read_dset_cd_5_
                  procedure :: alloc_read_dset_l_5_
                  procedure :: alloc_read_dset_s_5_
                  procedure :: alloc_read_dset_is_6_
                  procedure :: alloc_read_dset_id_6_
                  procedure :: alloc_read_dset_rs_6_
                  procedure :: alloc_read_dset_rd_6_
                  procedure :: alloc_read_dset_cs_6_
                  procedure :: alloc_read_dset_cd_6_
                  procedure :: alloc_read_dset_l_6_
                  procedure :: alloc_read_dset_s_6_
                  procedure :: alloc_read_dset_is_7_
                  procedure :: alloc_read_dset_id_7_
                  procedure :: alloc_read_dset_rs_7_
                  procedure :: alloc_read_dset_rd_7_
                  procedure :: alloc_read_dset_cs_7_
                  procedure :: alloc_read_dset_cd_7_
                  procedure :: alloc_read_dset_l_7_
                  procedure :: alloc_read_dset_s_7_
            generic, public :: alloc_read_dset => alloc_read_dset_is_0_,alloc_read_dset_id_0_,alloc_read_dset_rs_0_,alloc_read_dset&
                &_rd_0_,alloc_read_dset_cs_0_,alloc_read_dset_cd_0_,alloc_read_dset_l_0_,alloc_read_dset_s_0_,alloc_read_dset_is_1_&
                &,alloc_read_dset_id_1_,alloc_read_dset_rs_1_,alloc_read_dset_rd_1_,alloc_read_dset_cs_1_,alloc_read_dset_cd_1_,all&
                &oc_read_dset_l_1_,alloc_read_dset_s_1_,alloc_read_dset_is_2_,alloc_read_dset_id_2_,alloc_read_dset_rs_2_,alloc_rea&
                &d_dset_rd_2_,alloc_read_dset_cs_2_,alloc_read_dset_cd_2_,alloc_read_dset_l_2_,alloc_read_dset_s_2_,alloc_read_dset&
                &_is_3_,alloc_read_dset_id_3_,alloc_read_dset_rs_3_,alloc_read_dset_rd_3_,alloc_read_dset_cs_3_,alloc_read_dset_cd_&
                &3_,alloc_read_dset_l_3_,alloc_read_dset_s_3_,alloc_read_dset_is_4_,alloc_read_dset_id_4_,alloc_read_dset_rs_4_,all&
                &oc_read_dset_rd_4_,alloc_read_dset_cs_4_,alloc_read_dset_cd_4_,alloc_read_dset_l_4_,alloc_read_dset_s_4_,alloc_rea&
                &d_dset_is_5_,alloc_read_dset_id_5_,alloc_read_dset_rs_5_,alloc_read_dset_rd_5_,alloc_read_dset_cs_5_,alloc_read_ds&
                &et_cd_5_,alloc_read_dset_l_5_,alloc_read_dset_s_5_,alloc_read_dset_is_6_,alloc_read_dset_id_6_,alloc_read_dset_rs_&
                &6_,alloc_read_dset_rd_6_,alloc_read_dset_cs_6_,alloc_read_dset_cd_6_,alloc_read_dset_l_6_,alloc_read_dset_s_6_,all&
                &oc_read_dset_is_7_,alloc_read_dset_id_7_,alloc_read_dset_rs_7_,alloc_read_dset_rd_7_,alloc_read_dset_cs_7_,alloc_r&
                &ead_dset_cd_7_,alloc_read_dset_l_7_,alloc_read_dset_s_7_
                  procedure :: write_attr_is_0_
                  procedure :: write_attr_id_0_
                  procedure :: write_attr_rs_0_
                  procedure :: write_attr_rd_0_
                  procedure :: write_attr_cs_0_
                  procedure :: write_attr_cd_0_
                  procedure :: write_attr_l_0_
                  procedure :: write_attr_s_0_
                  procedure :: write_attr_is_1_
                  procedure :: write_attr_id_1_
                  procedure :: write_attr_rs_1_
                  procedure :: write_attr_rd_1_
                  procedure :: write_attr_cs_1_
                  procedure :: write_attr_cd_1_
                  procedure :: write_attr_l_1_
                  procedure :: write_attr_s_1_
                  procedure :: write_attr_is_2_
                  procedure :: write_attr_id_2_
                  procedure :: write_attr_rs_2_
                  procedure :: write_attr_rd_2_
                  procedure :: write_attr_cs_2_
                  procedure :: write_attr_cd_2_
                  procedure :: write_attr_l_2_
                  procedure :: write_attr_s_2_
                  procedure :: write_attr_is_3_
                  procedure :: write_attr_id_3_
                  procedure :: write_attr_rs_3_
                  procedure :: write_attr_rd_3_
                  procedure :: write_attr_cs_3_
                  procedure :: write_attr_cd_3_
                  procedure :: write_attr_l_3_
                  procedure :: write_attr_s_3_
                  procedure :: write_attr_is_4_
                  procedure :: write_attr_id_4_
                  procedure :: write_attr_rs_4_
                  procedure :: write_attr_rd_4_
                  procedure :: write_attr_cs_4_
                  procedure :: write_attr_cd_4_
                  procedure :: write_attr_l_4_
                  procedure :: write_attr_s_4_
                  procedure :: write_attr_is_5_
                  procedure :: write_attr_id_5_
                  procedure :: write_attr_rs_5_
                  procedure :: write_attr_rd_5_
                  procedure :: write_attr_cs_5_
                  procedure :: write_attr_cd_5_
                  procedure :: write_attr_l_5_
                  procedure :: write_attr_s_5_
                  procedure :: write_attr_is_6_
                  procedure :: write_attr_id_6_
                  procedure :: write_attr_rs_6_
                  procedure :: write_attr_rd_6_
                  procedure :: write_attr_cs_6_
                  procedure :: write_attr_cd_6_
                  procedure :: write_attr_l_6_
                  procedure :: write_attr_s_6_
                  procedure :: write_attr_is_7_
                  procedure :: write_attr_id_7_
                  procedure :: write_attr_rs_7_
                  procedure :: write_attr_rd_7_
                  procedure :: write_attr_cs_7_
                  procedure :: write_attr_cd_7_
                  procedure :: write_attr_l_7_
                  procedure :: write_attr_s_7_
            generic, public :: write_attr => write_attr_is_0_,write_attr_id_0_,write_attr_rs_0_,write_attr_rd_0_,write_attr_cs_0_,w&
                &rite_attr_cd_0_,write_attr_l_0_,write_attr_s_0_,write_attr_is_1_,write_attr_id_1_,write_attr_rs_1_,write_attr_rd_1&
                &_,write_attr_cs_1_,write_attr_cd_1_,write_attr_l_1_,write_attr_s_1_,write_attr_is_2_,write_attr_id_2_,write_attr_r&
                &s_2_,write_attr_rd_2_,write_attr_cs_2_,write_attr_cd_2_,write_attr_l_2_,write_attr_s_2_,write_attr_is_3_,write_att&
                &r_id_3_,write_attr_rs_3_,write_attr_rd_3_,write_attr_cs_3_,write_attr_cd_3_,write_attr_l_3_,write_attr_s_3_,write_&
                &attr_is_4_,write_attr_id_4_,write_attr_rs_4_,write_attr_rd_4_,write_attr_cs_4_,write_attr_cd_4_,write_attr_l_4_,wr&
                &ite_attr_s_4_,write_attr_is_5_,write_attr_id_5_,write_attr_rs_5_,write_attr_rd_5_,write_attr_cs_5_,write_attr_cd_5&
                &_,write_attr_l_5_,write_attr_s_5_,write_attr_is_6_,write_attr_id_6_,write_attr_rs_6_,write_attr_rd_6_,write_attr_c&
                &s_6_,write_attr_cd_6_,write_attr_l_6_,write_attr_s_6_,write_attr_is_7_,write_attr_id_7_,write_attr_rs_7_,write_att&
                &r_rd_7_,write_attr_cs_7_,write_attr_cd_7_,write_attr_l_7_,write_attr_s_7_
                  procedure :: write_dset_is_0_
                  procedure :: write_dset_id_0_
                  procedure :: write_dset_rs_0_
                  procedure :: write_dset_rd_0_
                  procedure :: write_dset_cs_0_
                  procedure :: write_dset_cd_0_
                  procedure :: write_dset_l_0_
                  procedure :: write_dset_s_0_
                  procedure :: write_dset_is_1_
                  procedure :: write_dset_id_1_
                  procedure :: write_dset_rs_1_
                  procedure :: write_dset_rd_1_
                  procedure :: write_dset_cs_1_
                  procedure :: write_dset_cd_1_
                  procedure :: write_dset_l_1_
                  procedure :: write_dset_s_1_
                  procedure :: write_dset_is_2_
                  procedure :: write_dset_id_2_
                  procedure :: write_dset_rs_2_
                  procedure :: write_dset_rd_2_
                  procedure :: write_dset_cs_2_
                  procedure :: write_dset_cd_2_
                  procedure :: write_dset_l_2_
                  procedure :: write_dset_s_2_
                  procedure :: write_dset_is_3_
                  procedure :: write_dset_id_3_
                  procedure :: write_dset_rs_3_
                  procedure :: write_dset_rd_3_
                  procedure :: write_dset_cs_3_
                  procedure :: write_dset_cd_3_
                  procedure :: write_dset_l_3_
                  procedure :: write_dset_s_3_
                  procedure :: write_dset_is_4_
                  procedure :: write_dset_id_4_
                  procedure :: write_dset_rs_4_
                  procedure :: write_dset_rd_4_
                  procedure :: write_dset_cs_4_
                  procedure :: write_dset_cd_4_
                  procedure :: write_dset_l_4_
                  procedure :: write_dset_s_4_
                  procedure :: write_dset_is_5_
                  procedure :: write_dset_id_5_
                  procedure :: write_dset_rs_5_
                  procedure :: write_dset_rd_5_
                  procedure :: write_dset_cs_5_
                  procedure :: write_dset_cd_5_
                  procedure :: write_dset_l_5_
                  procedure :: write_dset_s_5_
                  procedure :: write_dset_is_6_
                  procedure :: write_dset_id_6_
                  procedure :: write_dset_rs_6_
                  procedure :: write_dset_rd_6_
                  procedure :: write_dset_cs_6_
                  procedure :: write_dset_cd_6_
                  procedure :: write_dset_l_6_
                  procedure :: write_dset_s_6_
                  procedure :: write_dset_is_7_
                  procedure :: write_dset_id_7_
                  procedure :: write_dset_rs_7_
                  procedure :: write_dset_rd_7_
                  procedure :: write_dset_cs_7_
                  procedure :: write_dset_cd_7_
                  procedure :: write_dset_l_7_
                  procedure :: write_dset_s_7_
            generic, public :: write_dset => write_dset_is_0_,write_dset_id_0_,write_dset_rs_0_,write_dset_rd_0_,write_dset_cs_0_,w&
                &rite_dset_cd_0_,write_dset_l_0_,write_dset_s_0_,write_dset_is_1_,write_dset_id_1_,write_dset_rs_1_,write_dset_rd_1&
                &_,write_dset_cs_1_,write_dset_cd_1_,write_dset_l_1_,write_dset_s_1_,write_dset_is_2_,write_dset_id_2_,write_dset_r&
                &s_2_,write_dset_rd_2_,write_dset_cs_2_,write_dset_cd_2_,write_dset_l_2_,write_dset_s_2_,write_dset_is_3_,write_dse&
                &t_id_3_,write_dset_rs_3_,write_dset_rd_3_,write_dset_cs_3_,write_dset_cd_3_,write_dset_l_3_,write_dset_s_3_,write_&
                &dset_is_4_,write_dset_id_4_,write_dset_rs_4_,write_dset_rd_4_,write_dset_cs_4_,write_dset_cd_4_,write_dset_l_4_,wr&
                &ite_dset_s_4_,write_dset_is_5_,write_dset_id_5_,write_dset_rs_5_,write_dset_rd_5_,write_dset_cs_5_,write_dset_cd_5&
                &_,write_dset_l_5_,write_dset_s_5_,write_dset_is_6_,write_dset_id_6_,write_dset_rs_6_,write_dset_rd_6_,write_dset_c&
                &s_6_,write_dset_cd_6_,write_dset_l_6_,write_dset_s_6_,write_dset_is_7_,write_dset_id_7_,write_dset_rs_7_,write_dse&
                &t_rd_7_,write_dset_cs_7_,write_dset_cd_7_,write_dset_l_7_,write_dset_s_7_
   end type hdf5io_t

   ! Module variables
  
   integer, save :: ref_count = 0

      integer(HID_T), save :: mem_type_id_is
      integer(HID_T), save :: file_type_id_is
      integer(HID_T), save :: mem_type_id_id
      integer(HID_T), save :: file_type_id_id
      integer(HID_T), save :: mem_type_id_rs
      integer(HID_T), save :: file_type_id_rs
      integer(HID_T), save :: mem_type_id_rd
      integer(HID_T), save :: file_type_id_rd
      integer(HID_T), save :: mem_type_id_cs
      integer(HID_T), save :: file_type_id_cs
      integer(HID_T), save :: mem_type_id_cd
      integer(HID_T), save :: file_type_id_cd
      
   ! Interfaces

   interface hdf5io_t
      module procedure hdf5io_t_file_
      module procedure hdf5io_t_group_
   end interface hdf5io_t
   
   ! Access specifiers

   private

   public :: CREATE_FILE
   public :: OPEN_FILE
   public :: TYPE_LEN
   public :: hdf5io_t
   
   ! Procedures

contains

   function hdf5io_t_file_(file_name, access_type, comp_level) result (hi)

      character(*), intent(in)      :: file_name
      integer, intent(in)           :: access_type
      integer, intent(in), optional :: comp_level
      type(hdf5io_t)                :: hi

      integer        :: hdf_err
      integer(HID_T) :: file_id
      integer(HID_T) :: group_id

      ! If necessary, open the HDF5 library

      if (ref_count == 0) then
         call open_library_()
      endif

      ref_count = ref_count + 1

      ! Depending on the access_type, open or create the file

      select case(access_type)
      case(CREATE_FILE)
   call h5fcreate_f(file_name,H5F_ACC_TRUNC_F,file_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 132 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5fcreate_f'
   error stop
   endif
      case(OPEN_FILE)
   call h5fopen_f(file_name,H5F_ACC_RDWR_F,file_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 134 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5fopen_f'
   error stop
   endif
      case default
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 136 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'invalid access_type'
   error stop
      end select

      ! Open the root group

   call h5gopen_f(file_id,'/',group_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 141 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5gopen_f'
   error stop
   endif

      ! Construct the hdf5io_t

      hi%file_id = file_id
      hi%group_id = group_id
      if (PRESENT(comp_level)) then
         hi%comp_level = comp_level
      else
         hi%comp_level = 0
      endif

      allocate(hi%ref_count)
      hi%ref_count = 1

      ! Finish

      return

   end function hdf5io_t_file_

   !****

   function hdf5io_t_group_(hi_parent, group_name) result (hi)

      type(hdf5io_t), intent(inout) :: hi_parent
      character(*), intent(in)      :: group_name
      type(hdf5io_t)                :: hi

      integer        :: hdf_err
      integer(HID_T) :: group_id

      ! Depending on whether the group already exists, open or create it

      if (hi_parent%group_exists(group_name)) then
   call h5gopen_f(hi_parent%group_id,group_name,group_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 176 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5gopen_f'
   error stop
   endif
      else
   call h5gcreate_f(hi_parent%group_id,group_name,group_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 178 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5gcreate_f'
   error stop
   endif
      endif

      ! Construct the hdf5io_t

      hi%file_id = hi_parent%file_id
      hi%group_id = group_id
      hi%comp_level = hi_parent%comp_level

      hi%ref_count => hi_parent%ref_count
      hi%ref_count = hi%ref_count + 1

      ! Finish

      return

   end function hdf5io_t_group_

   !****

   subroutine final_(self)

      class(hdf5io_t), intent(inout) :: self

      integer :: hdf_err

      ! Close the group

   call h5gclose_f(self%group_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 206 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5gclose_f'
   error stop
   endif

      self%ref_count = self%ref_count - 1

      ! If necessary, close the file also

      if (self%ref_count == 0) then
   call h5fclose_f(self%file_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 213 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5fclose_f'
   error stop
   endif
         deallocate(self%ref_count)
         ref_count = ref_count - 1
      endif

      ! If necessary, close the HDF5 library

      if (ref_count == 0) then
         call close_library_()
      endif

      ! Finish

      return

   end subroutine final_

   !****

   subroutine open_library_()

      integer       :: hdf_err

      ! Open the HDF5 library
      
   call h5open_f(hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 238 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5open_f'
   error stop
   endif

   if (.NOT. (hdf_err == 0)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 240 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "hdf_err == 0" failed with message "Failed to open HDF5 library"'
   error stop
   end if

      ! Set up data types

      mem_type_id_is = h5kind_to_type(is, H5_INTEGER_KIND)
      mem_type_id_id = h5kind_to_type(id, H5_INTEGER_KIND)

      mem_type_id_rs = h5kind_to_type(rs, H5_REAL_KIND)
      mem_type_id_rd = h5kind_to_type(rd, H5_REAL_KIND)

      call create_complex_type_(mem_type_id_rs, mem_type_id_cs)
      call create_complex_type_(mem_type_id_rd, mem_type_id_cd)

      file_type_id_is = H5T_STD_I32LE
      file_type_id_id = H5T_STD_I64LE

      file_type_id_rs = H5T_IEEE_F32LE
      file_type_id_rd = H5T_IEEE_F64LE

      call create_complex_type_(file_type_id_rs, file_type_id_cs)
      call create_complex_type_(file_type_id_rd, file_type_id_cd)

      ! Finish

      return

   contains

      subroutine create_complex_type_(comp_type_id, type_id)

         integer(HID_T), intent(in)  :: comp_type_id
         integer(HID_T), intent(out) :: type_id

         integer         :: hdf_err
         integer(SIZE_T) :: comp_size

         ! Create a complex data type

   call h5tget_size_f(comp_type_id,comp_size,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 278 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tget_size_f'
   error stop
   endif

   call h5tcreate_f(H5T_COMPOUND_F,INT(2*comp_size, SIZE_T),type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 280 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcreate_f'
   error stop
   endif
   call h5tinsert_f(type_id,'re',INT(0, SIZE_T),comp_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 281 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tinsert_f'
   error stop
   endif
   call h5tinsert_f(type_id,'im',INT(comp_size, SIZE_T),comp_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 282 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tinsert_f'
   error stop
   endif

         ! Finish

         return

      end subroutine create_complex_type_

   end subroutine open_library_

   !****

   subroutine close_library_()

      integer :: hdf_err

      ! Close complex data types

   call h5tclose_f(mem_type_id_cs,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 301 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif
   call h5tclose_f(file_type_id_cs,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 302 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif
   call h5tclose_f(mem_type_id_cd,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 301 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif
   call h5tclose_f(file_type_id_cd,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 302 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

      ! Close the HDF5 library

   call h5close_f(hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 307 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5close_f'
   error stop
   endif

      ! Finish

      return

   end subroutine close_library_

   !****

   function group_exists(self, group_name)

      class(hdf5io_t), intent(inout) :: self
      character(*), intent(in)       :: group_name
      logical                        :: group_exists

      integer        :: hdf_err
      integer(HID_T) :: group_id

      ! Determine whether the named group already exists

      call h5eset_auto_f(0, hdf_err)

      call h5gopen_f(self%group_id, group_name, group_id, hdf_err)

      if (hdf_err >= 0) then
         group_exists = .TRUE.
         call h5gclose_f(group_id, hdf_err)
      else
         group_exists = .FALSE.
      endif

      call h5eset_auto_f(1, hdf_err)

      ! Finish

      return

   end function group_exists

   !****

   function attr_exists(self, attr_name)

      class(hdf5io_t), intent(inout) :: self
      character(*), intent(in)       :: attr_name
      logical                        :: attr_exists

      integer :: hdf_err

      ! Check if the attribute exists

   call h5aexists_f(self%group_id,attr_name,attr_exists,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 359 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aexists_f'
   error stop
   endif

      ! Finish

      return

   end function attr_exists

   !****

   function dset_exists(self, dset_name)

      class(hdf5io_t), intent(inout) :: self
      character(*), intent(in)       :: dset_name
      logical                        :: dset_exists

      integer          :: hdf_err
      logical          :: link_exists
      type(h5o_info_t) :: obj_info

      ! Check if the dataset exists

   call h5lexists_f(self%group_id,dset_name,link_exists,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 381 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5lexists_f'
   error stop
   endif

      if (link_exists) then

   call h5oget_info_by_name_f(self%group_id,dset_name,obj_info,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 385 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5oget_info_by_name_f'
   error stop
   endif

         dset_exists = obj_info%type == H5O_TYPE_DATASET_F

      else

         dset_exists = .FALSE.

      endif

      ! Finish

      return

   end function dset_exists

   !****


      function group_names(self) result(names)

         class(hdf5io_t), intent(inout) :: self
         character(:), allocatable      :: names(:)

         integer(HSIZE_T) :: i
         integer          :: n_names
         integer          :: max_len
         type(C_FUNPTR)   :: op
         integer          :: ret_value
         integer          :: hdf_err
         integer          :: i_name

         ! Iterate to get number of items and maximum name length

         i = 0

         n_names = 0
         max_len = 0

         op = C_FUNLOC(op_len_)

   call h5literate_f(self%group_id,H5_INDEX_NAME_F,H5_ITER_NATIVE_F,i,op,C_NULL_PTR,ret_value,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 427 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5literate_f'
   error stop
   endif

         ! Allocate space

         allocate(character(LEN=max_len)::names(n_names))

         ! Iterate again to store the names

         i = 0

         i_name = 0

         op = C_FUNLOC(op_name_)

   call h5literate_f(self%group_id,H5_INDEX_NAME_F,H5_ITER_NATIVE_F,i,op,C_NULL_PTR,ret_value,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 441 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5literate_f'
   error stop
   endif

         ! Finish

         return

      contains
      
         function op_len_(loc_id, c_name, c_info, c_data) result(ret_value) bind(C)

            integer(HID_T), value :: loc_id
            character(C_CHAR)     :: c_name(*)
            type(C_PTR), value    :: c_info
            type(C_PTR), value    :: c_data
            integer(C_INT)        :: ret_value

            type(H5O_INFO_T) :: info
            integer          :: hdf_err

            ! Check that the link points to the appropriate object type

   call h5oget_info_f(loc_id,info,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 462 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5oget_info_f'
   error stop
   endif

            if (info%type == H5O_TYPE_GROUP_F) then

               ! Update counters

               n_names = n_names + 1

               max_len = MAX(c_f_len(c_name), max_len)

            endif

            ! Finish

            ret_value = 0

         end function op_len_

         !****

         function op_name_(loc_id, c_name, c_info, c_data) result(ret_value) bind(C)

            integer(HID_T), value :: loc_id
            character(C_CHAR)     :: c_name(*)
            type(C_PTR), value    :: c_info
            type(C_PTR), value    :: c_data
            integer(C_INT)        :: ret_value

            type(H5O_INFO_T) :: info
            integer          :: hdf_err

            ! Check that the link points to the appropriate object type

   call h5oget_info_f(loc_id,info,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 495 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5oget_info_f'
   error stop
   endif

            if (info%type == H5O_TYPE_GROUP_F) then

               ! Store the name

               i_name = i_name + 1

               names(i_name) = c_f_string(c_name)

            end if

            ! Finish

            ret_value = 0

         end function op_name_

      end function group_names
         

      function dset_names(self) result(names)

         class(hdf5io_t), intent(inout) :: self
         character(:), allocatable      :: names(:)

         integer(HSIZE_T) :: i
         integer          :: n_names
         integer          :: max_len
         type(C_FUNPTR)   :: op
         integer          :: ret_value
         integer          :: hdf_err
         integer          :: i_name

         ! Iterate to get number of items and maximum name length

         i = 0

         n_names = 0
         max_len = 0

         op = C_FUNLOC(op_len_)

   call h5literate_f(self%group_id,H5_INDEX_NAME_F,H5_ITER_NATIVE_F,i,op,C_NULL_PTR,ret_value,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 427 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5literate_f'
   error stop
   endif

         ! Allocate space

         allocate(character(LEN=max_len)::names(n_names))

         ! Iterate again to store the names

         i = 0

         i_name = 0

         op = C_FUNLOC(op_name_)

   call h5literate_f(self%group_id,H5_INDEX_NAME_F,H5_ITER_NATIVE_F,i,op,C_NULL_PTR,ret_value,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 441 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5literate_f'
   error stop
   endif

         ! Finish

         return

      contains
      
         function op_len_(loc_id, c_name, c_info, c_data) result(ret_value) bind(C)

            integer(HID_T), value :: loc_id
            character(C_CHAR)     :: c_name(*)
            type(C_PTR), value    :: c_info
            type(C_PTR), value    :: c_data
            integer(C_INT)        :: ret_value

            type(H5O_INFO_T) :: info
            integer          :: hdf_err

            ! Check that the link points to the appropriate object type

   call h5oget_info_f(loc_id,info,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 462 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5oget_info_f'
   error stop
   endif

            if (info%type == H5O_TYPE_DATASET_F) then

               ! Update counters

               n_names = n_names + 1

               max_len = MAX(c_f_len(c_name), max_len)

            endif

            ! Finish

            ret_value = 0

         end function op_len_

         !****

         function op_name_(loc_id, c_name, c_info, c_data) result(ret_value) bind(C)

            integer(HID_T), value :: loc_id
            character(C_CHAR)     :: c_name(*)
            type(C_PTR), value    :: c_info
            type(C_PTR), value    :: c_data
            integer(C_INT)        :: ret_value

            type(H5O_INFO_T) :: info
            integer          :: hdf_err

            ! Check that the link points to the appropriate object type

   call h5oget_info_f(loc_id,info,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 495 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5oget_info_f'
   error stop
   endif

            if (info%type == H5O_TYPE_DATASET_F) then

               ! Store the name

               i_name = i_name + 1

               names(i_name) = c_f_string(c_name)

            end if

            ! Finish

            ret_value = 0

         end function op_name_

      end function dset_names
         

   !****


      function attr_shape(self, item_name) result(shape)

         class(hdf5io_t), intent(inout) :: self
         character(*), intent(in)       :: item_name
         integer(HSIZE_T), allocatable  :: shape(:)

         integer                       :: hdf_err
         integer(HID_T)                :: item_id
         integer(HID_T)                :: space_id
         integer                       :: rank
         integer(HSIZE_T), allocatable :: max_shape(:)

         ! Get the shape of the item

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 535 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif

   call h5aget_space_f(item_id,space_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 537 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aget_space_f'
   error stop
   endif

   call h5sget_simple_extent_ndims_f(space_id,rank,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 539 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sget_simple_extent_ndims_f'
   error stop
   endif

         allocate(shape(rank))
         allocate(max_shape(rank))

   call h5sget_simple_extent_dims_f(space_id,shape,max_shape,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 544 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sget_simple_extent_dims_f'
   error stop
   endif

   call h5sclose_f(space_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 546 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 547 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

         ! Finish

         return

      end function attr_shape


      function dset_shape(self, item_name) result(shape)

         class(hdf5io_t), intent(inout) :: self
         character(*), intent(in)       :: item_name
         integer(HSIZE_T), allocatable  :: shape(:)

         integer                       :: hdf_err
         integer(HID_T)                :: item_id
         integer(HID_T)                :: space_id
         integer                       :: rank
         integer(HSIZE_T), allocatable :: max_shape(:)

         ! Get the shape of the item

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 535 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif

   call h5dget_space_f(item_id,space_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 537 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dget_space_f'
   error stop
   endif

   call h5sget_simple_extent_ndims_f(space_id,rank,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 539 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sget_simple_extent_ndims_f'
   error stop
   endif

         allocate(shape(rank))
         allocate(max_shape(rank))

   call h5sget_simple_extent_dims_f(space_id,shape,max_shape,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 544 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sget_simple_extent_dims_f'
   error stop
   endif

   call h5sclose_f(space_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 546 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 547 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

         ! Finish

         return

      end function dset_shape


   !****



            subroutine read_attr_is_0_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(out)  :: data

               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


               ! Read the item

               mem_type_id = mem_type_id_is

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_is_0_


            subroutine read_attr_id_0_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(out)  :: data

               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


               ! Read the item

               mem_type_id = mem_type_id_id

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_id_0_


            subroutine read_attr_rs_0_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(out)  :: data

               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


               ! Read the item

               mem_type_id = mem_type_id_rs

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_rs_0_


            subroutine read_attr_rd_0_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(out)  :: data

               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


               ! Read the item

               mem_type_id = mem_type_id_rd

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_rd_0_


            subroutine read_attr_cs_0_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(out)  :: data

               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


               ! Read the item

               mem_type_id = mem_type_id_cs

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_cs_0_


            subroutine read_attr_cd_0_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(out)  :: data

               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


               ! Read the item

               mem_type_id = mem_type_id_cd

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_cd_0_


            
            subroutine read_attr_s_0_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(out)  :: data

               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

         
               ! Read the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 640 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 641 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 645 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 646 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 647 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 649 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_s_0_


            
            subroutine read_attr_l_0_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(out)  :: data

               integer, allocatable :: data_i

               ! Read the logical item

               call self%alloc_read_attr(item_name, data_i)

               data = data_i /= 0

               ! Finish

               return

            end subroutine read_attr_l_0_

            


            subroutine read_attr_is_1_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(out)  :: data(:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 1)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 1" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_is

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_is_1_


            subroutine read_attr_id_1_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(out)  :: data(:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 1)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 1" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_id

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_id_1_


            subroutine read_attr_rs_1_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(out)  :: data(:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 1)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 1" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_rs

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_rs_1_


            subroutine read_attr_rd_1_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(out)  :: data(:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 1)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 1" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_rd

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_rd_1_


            subroutine read_attr_cs_1_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(out)  :: data(:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 1)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 1" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_cs

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_cs_1_


            subroutine read_attr_cd_1_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(out)  :: data(:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 1)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 1" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_cd

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_cd_1_


            
            subroutine read_attr_s_1_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(out)  :: data(:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 1)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 633 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 1" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 634 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if

         
               ! Read the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 640 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 641 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 645 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 646 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 647 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 649 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_s_1_


            
            subroutine read_attr_l_1_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(out)  :: data(:)
                  contiguous :: data

               integer, allocatable :: data_i(:)

               ! Read the logical item

               call self%alloc_read_attr(item_name, data_i)

               data = data_i /= 0

               ! Finish

               return

            end subroutine read_attr_l_1_

            


            subroutine read_attr_is_2_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(out)  :: data(:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 2)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 2" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_is

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_is_2_


            subroutine read_attr_id_2_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(out)  :: data(:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 2)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 2" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_id

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_id_2_


            subroutine read_attr_rs_2_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(out)  :: data(:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 2)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 2" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_rs

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_rs_2_


            subroutine read_attr_rd_2_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(out)  :: data(:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 2)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 2" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_rd

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_rd_2_


            subroutine read_attr_cs_2_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(out)  :: data(:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 2)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 2" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_cs

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_cs_2_


            subroutine read_attr_cd_2_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(out)  :: data(:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 2)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 2" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_cd

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_cd_2_


            
            subroutine read_attr_s_2_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(out)  :: data(:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 2)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 633 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 2" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 634 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if

         
               ! Read the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 640 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 641 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 645 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 646 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 647 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 649 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_s_2_


            
            subroutine read_attr_l_2_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(out)  :: data(:,:)
                  contiguous :: data

               integer, allocatable :: data_i(:,:)

               ! Read the logical item

               call self%alloc_read_attr(item_name, data_i)

               data = data_i /= 0

               ! Finish

               return

            end subroutine read_attr_l_2_

            


            subroutine read_attr_is_3_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(out)  :: data(:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 3)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 3" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_is

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_is_3_


            subroutine read_attr_id_3_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(out)  :: data(:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 3)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 3" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_id

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_id_3_


            subroutine read_attr_rs_3_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(out)  :: data(:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 3)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 3" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_rs

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_rs_3_


            subroutine read_attr_rd_3_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(out)  :: data(:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 3)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 3" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_rd

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_rd_3_


            subroutine read_attr_cs_3_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(out)  :: data(:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 3)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 3" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_cs

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_cs_3_


            subroutine read_attr_cd_3_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(out)  :: data(:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 3)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 3" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_cd

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_cd_3_


            
            subroutine read_attr_s_3_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(out)  :: data(:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 3)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 633 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 3" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 634 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if

         
               ! Read the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 640 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 641 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 645 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 646 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 647 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 649 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_s_3_


            
            subroutine read_attr_l_3_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(out)  :: data(:,:,:)
                  contiguous :: data

               integer, allocatable :: data_i(:,:,:)

               ! Read the logical item

               call self%alloc_read_attr(item_name, data_i)

               data = data_i /= 0

               ! Finish

               return

            end subroutine read_attr_l_3_

            


            subroutine read_attr_is_4_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(out)  :: data(:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 4)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 4" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_is

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_is_4_


            subroutine read_attr_id_4_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(out)  :: data(:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 4)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 4" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_id

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_id_4_


            subroutine read_attr_rs_4_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(out)  :: data(:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 4)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 4" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_rs

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_rs_4_


            subroutine read_attr_rd_4_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(out)  :: data(:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 4)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 4" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_rd

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_rd_4_


            subroutine read_attr_cs_4_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(out)  :: data(:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 4)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 4" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_cs

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_cs_4_


            subroutine read_attr_cd_4_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(out)  :: data(:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 4)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 4" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_cd

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_cd_4_


            
            subroutine read_attr_s_4_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(out)  :: data(:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 4)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 633 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 4" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 634 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if

         
               ! Read the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 640 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 641 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 645 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 646 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 647 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 649 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_s_4_


            
            subroutine read_attr_l_4_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(out)  :: data(:,:,:,:)
                  contiguous :: data

               integer, allocatable :: data_i(:,:,:,:)

               ! Read the logical item

               call self%alloc_read_attr(item_name, data_i)

               data = data_i /= 0

               ! Finish

               return

            end subroutine read_attr_l_4_

            


            subroutine read_attr_is_5_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(out)  :: data(:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 5)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 5" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_is

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_is_5_


            subroutine read_attr_id_5_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(out)  :: data(:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 5)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 5" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_id

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_id_5_


            subroutine read_attr_rs_5_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(out)  :: data(:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 5)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 5" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_rs

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_rs_5_


            subroutine read_attr_rd_5_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(out)  :: data(:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 5)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 5" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_rd

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_rd_5_


            subroutine read_attr_cs_5_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(out)  :: data(:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 5)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 5" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_cs

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_cs_5_


            subroutine read_attr_cd_5_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(out)  :: data(:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 5)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 5" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_cd

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_cd_5_


            
            subroutine read_attr_s_5_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(out)  :: data(:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 5)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 633 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 5" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 634 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if

         
               ! Read the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 640 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 641 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 645 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 646 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 647 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 649 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_s_5_


            
            subroutine read_attr_l_5_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(out)  :: data(:,:,:,:,:)
                  contiguous :: data

               integer, allocatable :: data_i(:,:,:,:,:)

               ! Read the logical item

               call self%alloc_read_attr(item_name, data_i)

               data = data_i /= 0

               ! Finish

               return

            end subroutine read_attr_l_5_

            


            subroutine read_attr_is_6_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(out)  :: data(:,:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 6)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 6" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_is

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_is_6_


            subroutine read_attr_id_6_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(out)  :: data(:,:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 6)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 6" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_id

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_id_6_


            subroutine read_attr_rs_6_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(out)  :: data(:,:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 6)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 6" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_rs

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_rs_6_


            subroutine read_attr_rd_6_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(out)  :: data(:,:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 6)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 6" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_rd

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_rd_6_


            subroutine read_attr_cs_6_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(out)  :: data(:,:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 6)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 6" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_cs

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_cs_6_


            subroutine read_attr_cd_6_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(out)  :: data(:,:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 6)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 6" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_cd

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_cd_6_


            
            subroutine read_attr_s_6_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(out)  :: data(:,:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 6)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 633 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 6" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 634 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if

         
               ! Read the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 640 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 641 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 645 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 646 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 647 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 649 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_s_6_


            
            subroutine read_attr_l_6_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(out)  :: data(:,:,:,:,:,:)
                  contiguous :: data

               integer, allocatable :: data_i(:,:,:,:,:,:)

               ! Read the logical item

               call self%alloc_read_attr(item_name, data_i)

               data = data_i /= 0

               ! Finish

               return

            end subroutine read_attr_l_6_

            


            subroutine read_attr_is_7_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(out)  :: data(:,:,:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 7)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 7" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_is

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_is_7_


            subroutine read_attr_id_7_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(out)  :: data(:,:,:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 7)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 7" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_id

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_id_7_


            subroutine read_attr_rs_7_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(out)  :: data(:,:,:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 7)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 7" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_rs

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_rs_7_


            subroutine read_attr_rd_7_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(out)  :: data(:,:,:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 7)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 7" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_rd

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_rd_7_


            subroutine read_attr_cs_7_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(out)  :: data(:,:,:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 7)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 7" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_cs

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_cs_7_


            subroutine read_attr_cd_7_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(out)  :: data(:,:,:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 7)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 7" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_cd

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_cd_7_


            
            subroutine read_attr_s_7_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(out)  :: data(:,:,:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id


                  ! Check rank and shape

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 7)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 633 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 7" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 634 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if

         
               ! Read the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 640 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 641 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5aopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 645 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aopen_f'
   error stop
   endif
   call h5aread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 646 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aread_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 647 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 649 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_attr_s_7_


            
            subroutine read_attr_l_7_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(out)  :: data(:,:,:,:,:,:,:)
                  contiguous :: data

               integer, allocatable :: data_i(:,:,:,:,:,:,:)

               ! Read the logical item

               call self%alloc_read_attr(item_name, data_i)

               data = data_i /= 0

               ! Finish

               return

            end subroutine read_attr_l_7_

            


            subroutine read_dset_is_0_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(out)  :: data

               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


               ! Read the item

               mem_type_id = mem_type_id_is

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_is_0_


            subroutine read_dset_id_0_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(out)  :: data

               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


               ! Read the item

               mem_type_id = mem_type_id_id

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_id_0_


            subroutine read_dset_rs_0_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(out)  :: data

               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


               ! Read the item

               mem_type_id = mem_type_id_rs

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_rs_0_


            subroutine read_dset_rd_0_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(out)  :: data

               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


               ! Read the item

               mem_type_id = mem_type_id_rd

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_rd_0_


            subroutine read_dset_cs_0_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(out)  :: data

               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


               ! Read the item

               mem_type_id = mem_type_id_cs

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_cs_0_


            subroutine read_dset_cd_0_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(out)  :: data

               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


               ! Read the item

               mem_type_id = mem_type_id_cd

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_cd_0_


            
            subroutine read_dset_s_0_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(out)  :: data

               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

         
               ! Read the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 640 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 641 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 645 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 646 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 647 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 649 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_s_0_


            
            subroutine read_dset_l_0_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(out)  :: data

               integer, allocatable :: data_i

               ! Read the logical item

               call self%alloc_read_dset(item_name, data_i)

               data = data_i /= 0

               ! Finish

               return

            end subroutine read_dset_l_0_

            


            subroutine read_dset_is_1_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(out)  :: data(:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 1)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 1" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_is

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_is_1_


            subroutine read_dset_id_1_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(out)  :: data(:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 1)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 1" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_id

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_id_1_


            subroutine read_dset_rs_1_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(out)  :: data(:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 1)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 1" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_rs

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_rs_1_


            subroutine read_dset_rd_1_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(out)  :: data(:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 1)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 1" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_rd

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_rd_1_


            subroutine read_dset_cs_1_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(out)  :: data(:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 1)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 1" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_cs

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_cs_1_


            subroutine read_dset_cd_1_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(out)  :: data(:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 1)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 1" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_cd

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_cd_1_


            
            subroutine read_dset_s_1_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(out)  :: data(:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 1)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 633 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 1" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 634 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if

         
               ! Read the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 640 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 641 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 645 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 646 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 647 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 649 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_s_1_


            
            subroutine read_dset_l_1_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(out)  :: data(:)
                  contiguous :: data

               integer, allocatable :: data_i(:)

               ! Read the logical item

               call self%alloc_read_dset(item_name, data_i)

               data = data_i /= 0

               ! Finish

               return

            end subroutine read_dset_l_1_

            


            subroutine read_dset_is_2_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(out)  :: data(:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 2)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 2" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_is

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_is_2_


            subroutine read_dset_id_2_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(out)  :: data(:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 2)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 2" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_id

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_id_2_


            subroutine read_dset_rs_2_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(out)  :: data(:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 2)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 2" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_rs

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_rs_2_


            subroutine read_dset_rd_2_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(out)  :: data(:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 2)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 2" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_rd

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_rd_2_


            subroutine read_dset_cs_2_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(out)  :: data(:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 2)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 2" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_cs

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_cs_2_


            subroutine read_dset_cd_2_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(out)  :: data(:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 2)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 2" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_cd

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_cd_2_


            
            subroutine read_dset_s_2_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(out)  :: data(:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 2)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 633 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 2" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 634 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if

         
               ! Read the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 640 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 641 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 645 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 646 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 647 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 649 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_s_2_


            
            subroutine read_dset_l_2_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(out)  :: data(:,:)
                  contiguous :: data

               integer, allocatable :: data_i(:,:)

               ! Read the logical item

               call self%alloc_read_dset(item_name, data_i)

               data = data_i /= 0

               ! Finish

               return

            end subroutine read_dset_l_2_

            


            subroutine read_dset_is_3_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(out)  :: data(:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 3)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 3" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_is

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_is_3_


            subroutine read_dset_id_3_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(out)  :: data(:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 3)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 3" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_id

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_id_3_


            subroutine read_dset_rs_3_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(out)  :: data(:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 3)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 3" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_rs

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_rs_3_


            subroutine read_dset_rd_3_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(out)  :: data(:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 3)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 3" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_rd

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_rd_3_


            subroutine read_dset_cs_3_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(out)  :: data(:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 3)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 3" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_cs

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_cs_3_


            subroutine read_dset_cd_3_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(out)  :: data(:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 3)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 3" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_cd

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_cd_3_


            
            subroutine read_dset_s_3_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(out)  :: data(:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 3)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 633 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 3" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 634 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if

         
               ! Read the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 640 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 641 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 645 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 646 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 647 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 649 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_s_3_


            
            subroutine read_dset_l_3_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(out)  :: data(:,:,:)
                  contiguous :: data

               integer, allocatable :: data_i(:,:,:)

               ! Read the logical item

               call self%alloc_read_dset(item_name, data_i)

               data = data_i /= 0

               ! Finish

               return

            end subroutine read_dset_l_3_

            


            subroutine read_dset_is_4_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(out)  :: data(:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 4)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 4" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_is

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_is_4_


            subroutine read_dset_id_4_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(out)  :: data(:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 4)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 4" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_id

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_id_4_


            subroutine read_dset_rs_4_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(out)  :: data(:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 4)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 4" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_rs

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_rs_4_


            subroutine read_dset_rd_4_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(out)  :: data(:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 4)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 4" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_rd

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_rd_4_


            subroutine read_dset_cs_4_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(out)  :: data(:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 4)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 4" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_cs

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_cs_4_


            subroutine read_dset_cd_4_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(out)  :: data(:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 4)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 4" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_cd

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_cd_4_


            
            subroutine read_dset_s_4_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(out)  :: data(:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 4)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 633 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 4" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 634 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if

         
               ! Read the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 640 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 641 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 645 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 646 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 647 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 649 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_s_4_


            
            subroutine read_dset_l_4_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(out)  :: data(:,:,:,:)
                  contiguous :: data

               integer, allocatable :: data_i(:,:,:,:)

               ! Read the logical item

               call self%alloc_read_dset(item_name, data_i)

               data = data_i /= 0

               ! Finish

               return

            end subroutine read_dset_l_4_

            


            subroutine read_dset_is_5_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(out)  :: data(:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 5)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 5" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_is

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_is_5_


            subroutine read_dset_id_5_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(out)  :: data(:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 5)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 5" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_id

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_id_5_


            subroutine read_dset_rs_5_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(out)  :: data(:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 5)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 5" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_rs

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_rs_5_


            subroutine read_dset_rd_5_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(out)  :: data(:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 5)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 5" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_rd

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_rd_5_


            subroutine read_dset_cs_5_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(out)  :: data(:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 5)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 5" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_cs

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_cs_5_


            subroutine read_dset_cd_5_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(out)  :: data(:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 5)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 5" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_cd

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_cd_5_


            
            subroutine read_dset_s_5_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(out)  :: data(:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 5)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 633 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 5" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 634 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if

         
               ! Read the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 640 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 641 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 645 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 646 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 647 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 649 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_s_5_


            
            subroutine read_dset_l_5_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(out)  :: data(:,:,:,:,:)
                  contiguous :: data

               integer, allocatable :: data_i(:,:,:,:,:)

               ! Read the logical item

               call self%alloc_read_dset(item_name, data_i)

               data = data_i /= 0

               ! Finish

               return

            end subroutine read_dset_l_5_

            


            subroutine read_dset_is_6_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(out)  :: data(:,:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 6)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 6" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_is

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_is_6_


            subroutine read_dset_id_6_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(out)  :: data(:,:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 6)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 6" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_id

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_id_6_


            subroutine read_dset_rs_6_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(out)  :: data(:,:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 6)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 6" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_rs

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_rs_6_


            subroutine read_dset_rd_6_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(out)  :: data(:,:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 6)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 6" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_rd

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_rd_6_


            subroutine read_dset_cs_6_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(out)  :: data(:,:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 6)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 6" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_cs

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_cs_6_


            subroutine read_dset_cd_6_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(out)  :: data(:,:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 6)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 6" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_cd

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_cd_6_


            
            subroutine read_dset_s_6_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(out)  :: data(:,:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 6)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 633 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 6" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 634 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if

         
               ! Read the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 640 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 641 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 645 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 646 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 647 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 649 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_s_6_


            
            subroutine read_dset_l_6_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(out)  :: data(:,:,:,:,:,:)
                  contiguous :: data

               integer, allocatable :: data_i(:,:,:,:,:,:)

               ! Read the logical item

               call self%alloc_read_dset(item_name, data_i)

               data = data_i /= 0

               ! Finish

               return

            end subroutine read_dset_l_6_

            


            subroutine read_dset_is_7_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(out)  :: data(:,:,:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 7)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 7" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_is

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_is_7_


            subroutine read_dset_id_7_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(out)  :: data(:,:,:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 7)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 7" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_id

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_id_7_


            subroutine read_dset_rs_7_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(out)  :: data(:,:,:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 7)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 7" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_rs

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_rs_7_


            subroutine read_dset_rd_7_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(out)  :: data(:,:,:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 7)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 7" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_rd

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_rd_7_


            subroutine read_dset_cs_7_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(out)  :: data(:,:,:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 7)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 7" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_cs

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_cs_7_


            subroutine read_dset_cd_7_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(out)  :: data(:,:,:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer(HID_T) :: mem_type_id
               integer        :: hdf_err
               integer(HID_T) :: item_id
               type(C_PTR)    :: data_ptr


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 7)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 586 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 7" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 587 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if


               ! Read the item

               mem_type_id = mem_type_id_cd

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 597 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 598 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 599 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_cd_7_


            
            subroutine read_dset_s_7_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(out)  :: data(:,:,:,:,:,:,:)
                  contiguous :: data

                  integer(HSIZE_T), allocatable :: item_shape(:)
               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id


                  ! Check rank and shape

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 7)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 633 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 7" failed with message "rank mismatch"'
   error stop
   end if
   if (.NOT. (ALL(item_shape == SHAPE(data)))) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 634 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "ALL(item_shape == SHAPE(data))" failed with message "shape mismatch"'
   error stop
   end if

         
               ! Read the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 640 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 641 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5dopen_f(self%group_id,item_name,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 645 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dopen_f'
   error stop
   endif
   call h5dread_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 646 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dread_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 647 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 649 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine read_dset_s_7_


            
            subroutine read_dset_l_7_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(out)  :: data(:,:,:,:,:,:,:)
                  contiguous :: data

               integer, allocatable :: data_i(:,:,:,:,:,:,:)

               ! Read the logical item

               call self%alloc_read_dset(item_name, data_i)

               data = data_i /= 0

               ! Finish

               return

            end subroutine read_dset_l_7_

            

   !****

         
            subroutine alloc_read_attr_is_0_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(IS), allocatable, intent(out) :: data


               ! Allocate the item
                  

                  allocate(data)


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_is_0_

         
            subroutine alloc_read_attr_id_0_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(ID), allocatable, intent(out) :: data


               ! Allocate the item
                  

                  allocate(data)


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_id_0_

         
            subroutine alloc_read_attr_rs_0_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RS), allocatable, intent(out) :: data


               ! Allocate the item
                  

                  allocate(data)


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_rs_0_

         
            subroutine alloc_read_attr_rd_0_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RD), allocatable, intent(out) :: data


               ! Allocate the item
                  

                  allocate(data)


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_rd_0_

         
            subroutine alloc_read_attr_cs_0_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RS), allocatable, intent(out) :: data


               ! Allocate the item
                  

                  allocate(data)


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_cs_0_

         
            subroutine alloc_read_attr_cd_0_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RD), allocatable, intent(out) :: data


               ! Allocate the item
                  

                  allocate(data)


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_cd_0_

         
            subroutine alloc_read_attr_l_0_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               logical, allocatable, intent(out) :: data


               ! Allocate the item
                  

                  allocate(data)


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_l_0_

         
            subroutine alloc_read_attr_s_0_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               character(*), allocatable, intent(out) :: data


               ! Allocate the item
                  

                  allocate(data)


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_s_0_

         
            subroutine alloc_read_attr_is_1_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(IS), allocatable, intent(out) :: data(:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 1)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 1" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_is_1_

         
            subroutine alloc_read_attr_id_1_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(ID), allocatable, intent(out) :: data(:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 1)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 1" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_id_1_

         
            subroutine alloc_read_attr_rs_1_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RS), allocatable, intent(out) :: data(:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 1)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 1" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_rs_1_

         
            subroutine alloc_read_attr_rd_1_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RD), allocatable, intent(out) :: data(:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 1)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 1" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_rd_1_

         
            subroutine alloc_read_attr_cs_1_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RS), allocatable, intent(out) :: data(:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 1)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 1" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_cs_1_

         
            subroutine alloc_read_attr_cd_1_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RD), allocatable, intent(out) :: data(:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 1)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 1" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_cd_1_

         
            subroutine alloc_read_attr_l_1_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               logical, allocatable, intent(out) :: data(:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 1)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 1" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_l_1_

         
            subroutine alloc_read_attr_s_1_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               character(*), allocatable, intent(out) :: data(:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 1)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 1" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_s_1_

         
            subroutine alloc_read_attr_is_2_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(IS), allocatable, intent(out) :: data(:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 2)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 2" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_is_2_

         
            subroutine alloc_read_attr_id_2_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(ID), allocatable, intent(out) :: data(:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 2)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 2" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_id_2_

         
            subroutine alloc_read_attr_rs_2_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RS), allocatable, intent(out) :: data(:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 2)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 2" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_rs_2_

         
            subroutine alloc_read_attr_rd_2_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RD), allocatable, intent(out) :: data(:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 2)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 2" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_rd_2_

         
            subroutine alloc_read_attr_cs_2_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RS), allocatable, intent(out) :: data(:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 2)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 2" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_cs_2_

         
            subroutine alloc_read_attr_cd_2_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RD), allocatable, intent(out) :: data(:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 2)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 2" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_cd_2_

         
            subroutine alloc_read_attr_l_2_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               logical, allocatable, intent(out) :: data(:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 2)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 2" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_l_2_

         
            subroutine alloc_read_attr_s_2_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               character(*), allocatable, intent(out) :: data(:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 2)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 2" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_s_2_

         
            subroutine alloc_read_attr_is_3_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(IS), allocatable, intent(out) :: data(:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 3)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 3" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_is_3_

         
            subroutine alloc_read_attr_id_3_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(ID), allocatable, intent(out) :: data(:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 3)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 3" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_id_3_

         
            subroutine alloc_read_attr_rs_3_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RS), allocatable, intent(out) :: data(:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 3)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 3" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_rs_3_

         
            subroutine alloc_read_attr_rd_3_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RD), allocatable, intent(out) :: data(:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 3)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 3" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_rd_3_

         
            subroutine alloc_read_attr_cs_3_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RS), allocatable, intent(out) :: data(:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 3)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 3" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_cs_3_

         
            subroutine alloc_read_attr_cd_3_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RD), allocatable, intent(out) :: data(:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 3)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 3" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_cd_3_

         
            subroutine alloc_read_attr_l_3_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               logical, allocatable, intent(out) :: data(:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 3)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 3" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_l_3_

         
            subroutine alloc_read_attr_s_3_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               character(*), allocatable, intent(out) :: data(:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 3)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 3" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_s_3_

         
            subroutine alloc_read_attr_is_4_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(IS), allocatable, intent(out) :: data(:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 4)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 4" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_is_4_

         
            subroutine alloc_read_attr_id_4_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(ID), allocatable, intent(out) :: data(:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 4)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 4" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_id_4_

         
            subroutine alloc_read_attr_rs_4_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RS), allocatable, intent(out) :: data(:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 4)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 4" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_rs_4_

         
            subroutine alloc_read_attr_rd_4_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RD), allocatable, intent(out) :: data(:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 4)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 4" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_rd_4_

         
            subroutine alloc_read_attr_cs_4_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RS), allocatable, intent(out) :: data(:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 4)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 4" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_cs_4_

         
            subroutine alloc_read_attr_cd_4_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RD), allocatable, intent(out) :: data(:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 4)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 4" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_cd_4_

         
            subroutine alloc_read_attr_l_4_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               logical, allocatable, intent(out) :: data(:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 4)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 4" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_l_4_

         
            subroutine alloc_read_attr_s_4_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               character(*), allocatable, intent(out) :: data(:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 4)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 4" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_s_4_

         
            subroutine alloc_read_attr_is_5_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(IS), allocatable, intent(out) :: data(:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 5)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 5" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_is_5_

         
            subroutine alloc_read_attr_id_5_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(ID), allocatable, intent(out) :: data(:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 5)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 5" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_id_5_

         
            subroutine alloc_read_attr_rs_5_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RS), allocatable, intent(out) :: data(:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 5)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 5" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_rs_5_

         
            subroutine alloc_read_attr_rd_5_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RD), allocatable, intent(out) :: data(:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 5)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 5" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_rd_5_

         
            subroutine alloc_read_attr_cs_5_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RS), allocatable, intent(out) :: data(:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 5)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 5" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_cs_5_

         
            subroutine alloc_read_attr_cd_5_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RD), allocatable, intent(out) :: data(:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 5)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 5" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_cd_5_

         
            subroutine alloc_read_attr_l_5_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               logical, allocatable, intent(out) :: data(:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 5)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 5" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_l_5_

         
            subroutine alloc_read_attr_s_5_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               character(*), allocatable, intent(out) :: data(:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 5)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 5" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_s_5_

         
            subroutine alloc_read_attr_is_6_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(IS), allocatable, intent(out) :: data(:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 6)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 6" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_is_6_

         
            subroutine alloc_read_attr_id_6_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(ID), allocatable, intent(out) :: data(:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 6)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 6" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_id_6_

         
            subroutine alloc_read_attr_rs_6_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RS), allocatable, intent(out) :: data(:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 6)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 6" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_rs_6_

         
            subroutine alloc_read_attr_rd_6_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RD), allocatable, intent(out) :: data(:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 6)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 6" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_rd_6_

         
            subroutine alloc_read_attr_cs_6_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RS), allocatable, intent(out) :: data(:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 6)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 6" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_cs_6_

         
            subroutine alloc_read_attr_cd_6_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RD), allocatable, intent(out) :: data(:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 6)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 6" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_cd_6_

         
            subroutine alloc_read_attr_l_6_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               logical, allocatable, intent(out) :: data(:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 6)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 6" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_l_6_

         
            subroutine alloc_read_attr_s_6_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               character(*), allocatable, intent(out) :: data(:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 6)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 6" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_s_6_

         
            subroutine alloc_read_attr_is_7_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(IS), allocatable, intent(out) :: data(:,:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 7)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 7" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6),item_shape(7)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_is_7_

         
            subroutine alloc_read_attr_id_7_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(ID), allocatable, intent(out) :: data(:,:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 7)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 7" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6),item_shape(7)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_id_7_

         
            subroutine alloc_read_attr_rs_7_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RS), allocatable, intent(out) :: data(:,:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 7)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 7" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6),item_shape(7)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_rs_7_

         
            subroutine alloc_read_attr_rd_7_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RD), allocatable, intent(out) :: data(:,:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 7)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 7" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6),item_shape(7)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_rd_7_

         
            subroutine alloc_read_attr_cs_7_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RS), allocatable, intent(out) :: data(:,:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 7)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 7" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6),item_shape(7)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_cs_7_

         
            subroutine alloc_read_attr_cd_7_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RD), allocatable, intent(out) :: data(:,:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 7)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 7" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6),item_shape(7)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_cd_7_

         
            subroutine alloc_read_attr_l_7_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               logical, allocatable, intent(out) :: data(:,:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 7)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 7" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6),item_shape(7)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_l_7_

         
            subroutine alloc_read_attr_s_7_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               character(*), allocatable, intent(out) :: data(:,:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%attr_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 7)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 7" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6),item_shape(7)))


               ! Read the item

               call self%read_attr(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_attr_s_7_

         
            subroutine alloc_read_dset_is_0_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(IS), allocatable, intent(out) :: data


               ! Allocate the item
                  

                  allocate(data)


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_is_0_

         
            subroutine alloc_read_dset_id_0_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(ID), allocatable, intent(out) :: data


               ! Allocate the item
                  

                  allocate(data)


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_id_0_

         
            subroutine alloc_read_dset_rs_0_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RS), allocatable, intent(out) :: data


               ! Allocate the item
                  

                  allocate(data)


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_rs_0_

         
            subroutine alloc_read_dset_rd_0_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RD), allocatable, intent(out) :: data


               ! Allocate the item
                  

                  allocate(data)


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_rd_0_

         
            subroutine alloc_read_dset_cs_0_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RS), allocatable, intent(out) :: data


               ! Allocate the item
                  

                  allocate(data)


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_cs_0_

         
            subroutine alloc_read_dset_cd_0_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RD), allocatable, intent(out) :: data


               ! Allocate the item
                  

                  allocate(data)


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_cd_0_

         
            subroutine alloc_read_dset_l_0_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               logical, allocatable, intent(out) :: data


               ! Allocate the item
                  

                  allocate(data)


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_l_0_

         
            subroutine alloc_read_dset_s_0_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               character(*), allocatable, intent(out) :: data


               ! Allocate the item
                  

                  allocate(data)


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_s_0_

         
            subroutine alloc_read_dset_is_1_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(IS), allocatable, intent(out) :: data(:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 1)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 1" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_is_1_

         
            subroutine alloc_read_dset_id_1_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(ID), allocatable, intent(out) :: data(:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 1)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 1" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_id_1_

         
            subroutine alloc_read_dset_rs_1_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RS), allocatable, intent(out) :: data(:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 1)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 1" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_rs_1_

         
            subroutine alloc_read_dset_rd_1_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RD), allocatable, intent(out) :: data(:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 1)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 1" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_rd_1_

         
            subroutine alloc_read_dset_cs_1_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RS), allocatable, intent(out) :: data(:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 1)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 1" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_cs_1_

         
            subroutine alloc_read_dset_cd_1_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RD), allocatable, intent(out) :: data(:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 1)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 1" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_cd_1_

         
            subroutine alloc_read_dset_l_1_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               logical, allocatable, intent(out) :: data(:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 1)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 1" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_l_1_

         
            subroutine alloc_read_dset_s_1_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               character(*), allocatable, intent(out) :: data(:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 1)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 1" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_s_1_

         
            subroutine alloc_read_dset_is_2_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(IS), allocatable, intent(out) :: data(:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 2)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 2" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_is_2_

         
            subroutine alloc_read_dset_id_2_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(ID), allocatable, intent(out) :: data(:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 2)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 2" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_id_2_

         
            subroutine alloc_read_dset_rs_2_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RS), allocatable, intent(out) :: data(:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 2)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 2" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_rs_2_

         
            subroutine alloc_read_dset_rd_2_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RD), allocatable, intent(out) :: data(:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 2)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 2" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_rd_2_

         
            subroutine alloc_read_dset_cs_2_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RS), allocatable, intent(out) :: data(:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 2)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 2" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_cs_2_

         
            subroutine alloc_read_dset_cd_2_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RD), allocatable, intent(out) :: data(:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 2)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 2" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_cd_2_

         
            subroutine alloc_read_dset_l_2_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               logical, allocatable, intent(out) :: data(:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 2)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 2" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_l_2_

         
            subroutine alloc_read_dset_s_2_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               character(*), allocatable, intent(out) :: data(:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 2)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 2" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_s_2_

         
            subroutine alloc_read_dset_is_3_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(IS), allocatable, intent(out) :: data(:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 3)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 3" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_is_3_

         
            subroutine alloc_read_dset_id_3_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(ID), allocatable, intent(out) :: data(:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 3)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 3" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_id_3_

         
            subroutine alloc_read_dset_rs_3_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RS), allocatable, intent(out) :: data(:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 3)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 3" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_rs_3_

         
            subroutine alloc_read_dset_rd_3_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RD), allocatable, intent(out) :: data(:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 3)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 3" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_rd_3_

         
            subroutine alloc_read_dset_cs_3_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RS), allocatable, intent(out) :: data(:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 3)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 3" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_cs_3_

         
            subroutine alloc_read_dset_cd_3_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RD), allocatable, intent(out) :: data(:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 3)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 3" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_cd_3_

         
            subroutine alloc_read_dset_l_3_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               logical, allocatable, intent(out) :: data(:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 3)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 3" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_l_3_

         
            subroutine alloc_read_dset_s_3_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               character(*), allocatable, intent(out) :: data(:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 3)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 3" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_s_3_

         
            subroutine alloc_read_dset_is_4_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(IS), allocatable, intent(out) :: data(:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 4)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 4" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_is_4_

         
            subroutine alloc_read_dset_id_4_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(ID), allocatable, intent(out) :: data(:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 4)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 4" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_id_4_

         
            subroutine alloc_read_dset_rs_4_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RS), allocatable, intent(out) :: data(:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 4)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 4" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_rs_4_

         
            subroutine alloc_read_dset_rd_4_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RD), allocatable, intent(out) :: data(:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 4)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 4" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_rd_4_

         
            subroutine alloc_read_dset_cs_4_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RS), allocatable, intent(out) :: data(:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 4)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 4" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_cs_4_

         
            subroutine alloc_read_dset_cd_4_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RD), allocatable, intent(out) :: data(:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 4)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 4" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_cd_4_

         
            subroutine alloc_read_dset_l_4_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               logical, allocatable, intent(out) :: data(:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 4)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 4" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_l_4_

         
            subroutine alloc_read_dset_s_4_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               character(*), allocatable, intent(out) :: data(:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 4)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 4" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_s_4_

         
            subroutine alloc_read_dset_is_5_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(IS), allocatable, intent(out) :: data(:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 5)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 5" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_is_5_

         
            subroutine alloc_read_dset_id_5_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(ID), allocatable, intent(out) :: data(:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 5)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 5" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_id_5_

         
            subroutine alloc_read_dset_rs_5_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RS), allocatable, intent(out) :: data(:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 5)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 5" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_rs_5_

         
            subroutine alloc_read_dset_rd_5_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RD), allocatable, intent(out) :: data(:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 5)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 5" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_rd_5_

         
            subroutine alloc_read_dset_cs_5_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RS), allocatable, intent(out) :: data(:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 5)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 5" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_cs_5_

         
            subroutine alloc_read_dset_cd_5_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RD), allocatable, intent(out) :: data(:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 5)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 5" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_cd_5_

         
            subroutine alloc_read_dset_l_5_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               logical, allocatable, intent(out) :: data(:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 5)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 5" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_l_5_

         
            subroutine alloc_read_dset_s_5_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               character(*), allocatable, intent(out) :: data(:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 5)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 5" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_s_5_

         
            subroutine alloc_read_dset_is_6_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(IS), allocatable, intent(out) :: data(:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 6)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 6" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_is_6_

         
            subroutine alloc_read_dset_id_6_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(ID), allocatable, intent(out) :: data(:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 6)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 6" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_id_6_

         
            subroutine alloc_read_dset_rs_6_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RS), allocatable, intent(out) :: data(:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 6)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 6" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_rs_6_

         
            subroutine alloc_read_dset_rd_6_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RD), allocatable, intent(out) :: data(:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 6)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 6" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_rd_6_

         
            subroutine alloc_read_dset_cs_6_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RS), allocatable, intent(out) :: data(:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 6)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 6" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_cs_6_

         
            subroutine alloc_read_dset_cd_6_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RD), allocatable, intent(out) :: data(:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 6)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 6" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_cd_6_

         
            subroutine alloc_read_dset_l_6_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               logical, allocatable, intent(out) :: data(:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 6)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 6" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_l_6_

         
            subroutine alloc_read_dset_s_6_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               character(*), allocatable, intent(out) :: data(:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 6)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 6" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_s_6_

         
            subroutine alloc_read_dset_is_7_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(IS), allocatable, intent(out) :: data(:,:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 7)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 7" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6),item_shape(7)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_is_7_

         
            subroutine alloc_read_dset_id_7_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               integer(ID), allocatable, intent(out) :: data(:,:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 7)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 7" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6),item_shape(7)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_id_7_

         
            subroutine alloc_read_dset_rs_7_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RS), allocatable, intent(out) :: data(:,:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 7)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 7" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6),item_shape(7)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_rs_7_

         
            subroutine alloc_read_dset_rd_7_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               real(RD), allocatable, intent(out) :: data(:,:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 7)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 7" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6),item_shape(7)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_rd_7_

         
            subroutine alloc_read_dset_cs_7_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RS), allocatable, intent(out) :: data(:,:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 7)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 7" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6),item_shape(7)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_cs_7_

         
            subroutine alloc_read_dset_cd_7_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               complex(RD), allocatable, intent(out) :: data(:,:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 7)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 7" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6),item_shape(7)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_cd_7_

         
            subroutine alloc_read_dset_l_7_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               logical, allocatable, intent(out) :: data(:,:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 7)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 7" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6),item_shape(7)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_l_7_

         
            subroutine alloc_read_dset_s_7_(self, item_name, data)

               class(hdf5io_t), intent(inout)     :: self
               character(*), intent(in)           :: item_name
               character(*), allocatable, intent(out) :: data(:,:,:,:,:,:,:)

                  integer(HSIZE_T), allocatable :: item_shape(:)

               ! Allocate the item
                  

                  item_shape = self%dset_shape(item_name)
   if (.NOT. (SIZE(item_shape) == 7)) then
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 710 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'assertion "SIZE(item_shape) == 7" failed with message "rank mismatch"'
   error stop
   end if
                  
                  allocate(data(item_shape(1),item_shape(2),item_shape(3),item_shape(4),item_shape(5),item_shape(6),item_shape(7)))


               ! Read the item

               call self%read_dset(item_name, data)

               ! Finish

               return

            end subroutine alloc_read_dset_s_7_


   !****



            subroutine write_attr_is_0_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(in)   :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_f(H5S_SCALAR_F,dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 763 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_is
               file_type_id = file_type_id_is

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_is_0_


            subroutine write_attr_id_0_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(in)   :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_f(H5S_SCALAR_F,dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 763 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_id
               file_type_id = file_type_id_id

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_id_0_


            subroutine write_attr_rs_0_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(in)   :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_f(H5S_SCALAR_F,dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 763 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rs
               file_type_id = file_type_id_rs

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_rs_0_


            subroutine write_attr_rd_0_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(in)   :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_f(H5S_SCALAR_F,dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 763 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rd
               file_type_id = file_type_id_rd

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_rd_0_


            subroutine write_attr_cs_0_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(in)   :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_f(H5S_SCALAR_F,dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 763 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cs
               file_type_id = file_type_id_cs

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_cs_0_


            subroutine write_attr_cd_0_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(in)   :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_f(H5S_SCALAR_F,dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 763 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cd
               file_type_id = file_type_id_cd

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_cd_0_



            subroutine write_attr_s_0_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(in)   :: data

               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               integer(HID_T) :: dspace_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 820 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 821 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5tcopy_f(H5T_NATIVE_CHARACTER,file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 823 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(file_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 824 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5screate_f(H5S_SCALAR_F,dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 829 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 834 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif
   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 835 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 836 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 838 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 840 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif
   call h5tclose_f(file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 841 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_s_0_



            subroutine write_attr_l_0_ (self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(in)   :: data

               ! Write the logical item

               call self%write_attr(item_name, MERGE(1, 0, MASK=data))

               ! Finish

               return

            end subroutine write_attr_l_0_




            subroutine write_attr_is_1_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(in)   :: data(:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(1,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_is
               file_type_id = file_type_id_is

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_is_1_


            subroutine write_attr_id_1_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(in)   :: data(:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(1,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_id
               file_type_id = file_type_id_id

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_id_1_


            subroutine write_attr_rs_1_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(in)   :: data(:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(1,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rs
               file_type_id = file_type_id_rs

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_rs_1_


            subroutine write_attr_rd_1_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(in)   :: data(:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(1,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rd
               file_type_id = file_type_id_rd

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_rd_1_


            subroutine write_attr_cs_1_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(in)   :: data(:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(1,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cs
               file_type_id = file_type_id_cs

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_cs_1_


            subroutine write_attr_cd_1_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(in)   :: data(:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(1,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cd
               file_type_id = file_type_id_cd

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_cd_1_



            subroutine write_attr_s_1_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(in)   :: data(:)
                  contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               integer(HID_T) :: dspace_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 820 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 821 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5tcopy_f(H5T_NATIVE_CHARACTER,file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 823 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(file_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 824 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5screate_simple_f(1,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 827 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 834 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif
   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 835 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 836 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 838 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 840 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif
   call h5tclose_f(file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 841 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_s_1_



            subroutine write_attr_l_1_ (self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(in)   :: data(:)
                  contiguous :: data

               ! Write the logical item

               call self%write_attr(item_name, MERGE(1, 0, MASK=data))

               ! Finish

               return

            end subroutine write_attr_l_1_




            subroutine write_attr_is_2_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(in)   :: data(:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(2,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_is
               file_type_id = file_type_id_is

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_is_2_


            subroutine write_attr_id_2_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(in)   :: data(:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(2,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_id
               file_type_id = file_type_id_id

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_id_2_


            subroutine write_attr_rs_2_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(in)   :: data(:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(2,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rs
               file_type_id = file_type_id_rs

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_rs_2_


            subroutine write_attr_rd_2_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(in)   :: data(:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(2,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rd
               file_type_id = file_type_id_rd

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_rd_2_


            subroutine write_attr_cs_2_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(in)   :: data(:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(2,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cs
               file_type_id = file_type_id_cs

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_cs_2_


            subroutine write_attr_cd_2_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(in)   :: data(:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(2,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cd
               file_type_id = file_type_id_cd

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_cd_2_



            subroutine write_attr_s_2_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(in)   :: data(:,:)
                  contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               integer(HID_T) :: dspace_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 820 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 821 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5tcopy_f(H5T_NATIVE_CHARACTER,file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 823 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(file_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 824 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5screate_simple_f(2,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 827 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 834 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif
   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 835 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 836 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 838 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 840 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif
   call h5tclose_f(file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 841 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_s_2_



            subroutine write_attr_l_2_ (self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(in)   :: data(:,:)
                  contiguous :: data

               ! Write the logical item

               call self%write_attr(item_name, MERGE(1, 0, MASK=data))

               ! Finish

               return

            end subroutine write_attr_l_2_




            subroutine write_attr_is_3_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(in)   :: data(:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(3,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_is
               file_type_id = file_type_id_is

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_is_3_


            subroutine write_attr_id_3_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(in)   :: data(:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(3,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_id
               file_type_id = file_type_id_id

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_id_3_


            subroutine write_attr_rs_3_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(in)   :: data(:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(3,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rs
               file_type_id = file_type_id_rs

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_rs_3_


            subroutine write_attr_rd_3_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(in)   :: data(:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(3,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rd
               file_type_id = file_type_id_rd

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_rd_3_


            subroutine write_attr_cs_3_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(in)   :: data(:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(3,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cs
               file_type_id = file_type_id_cs

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_cs_3_


            subroutine write_attr_cd_3_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(in)   :: data(:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(3,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cd
               file_type_id = file_type_id_cd

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_cd_3_



            subroutine write_attr_s_3_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(in)   :: data(:,:,:)
                  contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               integer(HID_T) :: dspace_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 820 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 821 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5tcopy_f(H5T_NATIVE_CHARACTER,file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 823 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(file_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 824 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5screate_simple_f(3,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 827 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 834 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif
   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 835 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 836 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 838 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 840 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif
   call h5tclose_f(file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 841 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_s_3_



            subroutine write_attr_l_3_ (self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(in)   :: data(:,:,:)
                  contiguous :: data

               ! Write the logical item

               call self%write_attr(item_name, MERGE(1, 0, MASK=data))

               ! Finish

               return

            end subroutine write_attr_l_3_




            subroutine write_attr_is_4_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(in)   :: data(:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(4,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_is
               file_type_id = file_type_id_is

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_is_4_


            subroutine write_attr_id_4_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(in)   :: data(:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(4,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_id
               file_type_id = file_type_id_id

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_id_4_


            subroutine write_attr_rs_4_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(in)   :: data(:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(4,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rs
               file_type_id = file_type_id_rs

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_rs_4_


            subroutine write_attr_rd_4_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(in)   :: data(:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(4,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rd
               file_type_id = file_type_id_rd

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_rd_4_


            subroutine write_attr_cs_4_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(in)   :: data(:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(4,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cs
               file_type_id = file_type_id_cs

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_cs_4_


            subroutine write_attr_cd_4_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(in)   :: data(:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(4,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cd
               file_type_id = file_type_id_cd

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_cd_4_



            subroutine write_attr_s_4_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(in)   :: data(:,:,:,:)
                  contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               integer(HID_T) :: dspace_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 820 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 821 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5tcopy_f(H5T_NATIVE_CHARACTER,file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 823 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(file_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 824 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5screate_simple_f(4,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 827 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 834 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif
   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 835 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 836 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 838 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 840 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif
   call h5tclose_f(file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 841 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_s_4_



            subroutine write_attr_l_4_ (self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(in)   :: data(:,:,:,:)
                  contiguous :: data

               ! Write the logical item

               call self%write_attr(item_name, MERGE(1, 0, MASK=data))

               ! Finish

               return

            end subroutine write_attr_l_4_




            subroutine write_attr_is_5_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(in)   :: data(:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(5,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_is
               file_type_id = file_type_id_is

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_is_5_


            subroutine write_attr_id_5_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(in)   :: data(:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(5,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_id
               file_type_id = file_type_id_id

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_id_5_


            subroutine write_attr_rs_5_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(in)   :: data(:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(5,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rs
               file_type_id = file_type_id_rs

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_rs_5_


            subroutine write_attr_rd_5_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(in)   :: data(:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(5,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rd
               file_type_id = file_type_id_rd

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_rd_5_


            subroutine write_attr_cs_5_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(in)   :: data(:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(5,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cs
               file_type_id = file_type_id_cs

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_cs_5_


            subroutine write_attr_cd_5_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(in)   :: data(:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(5,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cd
               file_type_id = file_type_id_cd

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_cd_5_



            subroutine write_attr_s_5_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(in)   :: data(:,:,:,:,:)
                  contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               integer(HID_T) :: dspace_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 820 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 821 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5tcopy_f(H5T_NATIVE_CHARACTER,file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 823 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(file_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 824 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5screate_simple_f(5,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 827 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 834 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif
   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 835 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 836 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 838 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 840 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif
   call h5tclose_f(file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 841 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_s_5_



            subroutine write_attr_l_5_ (self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(in)   :: data(:,:,:,:,:)
                  contiguous :: data

               ! Write the logical item

               call self%write_attr(item_name, MERGE(1, 0, MASK=data))

               ! Finish

               return

            end subroutine write_attr_l_5_




            subroutine write_attr_is_6_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(in)   :: data(:,:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(6,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_is
               file_type_id = file_type_id_is

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_is_6_


            subroutine write_attr_id_6_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(in)   :: data(:,:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(6,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_id
               file_type_id = file_type_id_id

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_id_6_


            subroutine write_attr_rs_6_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(in)   :: data(:,:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(6,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rs
               file_type_id = file_type_id_rs

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_rs_6_


            subroutine write_attr_rd_6_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(in)   :: data(:,:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(6,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rd
               file_type_id = file_type_id_rd

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_rd_6_


            subroutine write_attr_cs_6_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(in)   :: data(:,:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(6,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cs
               file_type_id = file_type_id_cs

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_cs_6_


            subroutine write_attr_cd_6_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(in)   :: data(:,:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(6,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cd
               file_type_id = file_type_id_cd

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_cd_6_



            subroutine write_attr_s_6_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(in)   :: data(:,:,:,:,:,:)
                  contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               integer(HID_T) :: dspace_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 820 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 821 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5tcopy_f(H5T_NATIVE_CHARACTER,file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 823 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(file_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 824 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5screate_simple_f(6,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 827 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 834 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif
   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 835 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 836 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 838 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 840 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif
   call h5tclose_f(file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 841 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_s_6_



            subroutine write_attr_l_6_ (self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(in)   :: data(:,:,:,:,:,:)
                  contiguous :: data

               ! Write the logical item

               call self%write_attr(item_name, MERGE(1, 0, MASK=data))

               ! Finish

               return

            end subroutine write_attr_l_6_




            subroutine write_attr_is_7_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(in)   :: data(:,:,:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(7,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_is
               file_type_id = file_type_id_is

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_is_7_


            subroutine write_attr_id_7_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(in)   :: data(:,:,:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(7,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_id
               file_type_id = file_type_id_id

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_id_7_


            subroutine write_attr_rs_7_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(in)   :: data(:,:,:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(7,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rs
               file_type_id = file_type_id_rs

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_rs_7_


            subroutine write_attr_rd_7_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(in)   :: data(:,:,:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(7,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rd
               file_type_id = file_type_id_rd

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_rd_7_


            subroutine write_attr_cs_7_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(in)   :: data(:,:,:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(7,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cs
               file_type_id = file_type_id_cs

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_cs_7_


            subroutine write_attr_cd_7_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(in)   :: data(:,:,:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(7,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_ATTRIBUTE_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 767 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cd
               file_type_id = file_type_id_cd

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,acpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 782 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif

   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_cd_7_



            subroutine write_attr_s_7_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(in)   :: data(:,:,:,:,:,:,:)
                  contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               integer(HID_T) :: dspace_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 820 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 821 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5tcopy_f(H5T_NATIVE_CHARACTER,file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 823 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(file_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 824 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5screate_simple_f(7,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 827 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5acreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 834 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5acreate_f'
   error stop
   endif
   call h5awrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 835 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5awrite_f'
   error stop
   endif
   call h5aclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 836 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5aclose_f'
   error stop
   endif

   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 838 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 840 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif
   call h5tclose_f(file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 841 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_attr_s_7_



            subroutine write_attr_l_7_ (self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(in)   :: data(:,:,:,:,:,:,:)
                  contiguous :: data

               ! Write the logical item

               call self%write_attr(item_name, MERGE(1, 0, MASK=data))

               ! Finish

               return

            end subroutine write_attr_l_7_




            subroutine write_dset_is_0_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(in)   :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_f(H5S_SCALAR_F,dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 763 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_is
               file_type_id = file_type_id_is

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_is_0_


            subroutine write_dset_id_0_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(in)   :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_f(H5S_SCALAR_F,dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 763 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_id
               file_type_id = file_type_id_id

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_id_0_


            subroutine write_dset_rs_0_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(in)   :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_f(H5S_SCALAR_F,dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 763 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rs
               file_type_id = file_type_id_rs

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_rs_0_


            subroutine write_dset_rd_0_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(in)   :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_f(H5S_SCALAR_F,dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 763 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rd
               file_type_id = file_type_id_rd

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_rd_0_


            subroutine write_dset_cs_0_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(in)   :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_f(H5S_SCALAR_F,dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 763 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cs
               file_type_id = file_type_id_cs

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_cs_0_


            subroutine write_dset_cd_0_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(in)   :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_f(H5S_SCALAR_F,dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 763 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cd
               file_type_id = file_type_id_cd

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_cd_0_



            subroutine write_dset_s_0_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(in)   :: data

               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               integer(HID_T) :: dspace_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 820 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 821 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5tcopy_f(H5T_NATIVE_CHARACTER,file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 823 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(file_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 824 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5screate_f(H5S_SCALAR_F,dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 829 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 834 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif
   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 835 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 836 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 838 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 840 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif
   call h5tclose_f(file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 841 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_s_0_



            subroutine write_dset_l_0_ (self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(in)   :: data

               ! Write the logical item

               call self%write_dset(item_name, MERGE(1, 0, MASK=data))

               ! Finish

               return

            end subroutine write_dset_l_0_




            subroutine write_dset_is_1_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(in)   :: data(:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(1,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,1,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_is
               file_type_id = file_type_id_is

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_is_1_


            subroutine write_dset_id_1_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(in)   :: data(:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(1,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,1,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_id
               file_type_id = file_type_id_id

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_id_1_


            subroutine write_dset_rs_1_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(in)   :: data(:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(1,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,1,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rs
               file_type_id = file_type_id_rs

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_rs_1_


            subroutine write_dset_rd_1_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(in)   :: data(:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(1,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,1,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rd
               file_type_id = file_type_id_rd

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_rd_1_


            subroutine write_dset_cs_1_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(in)   :: data(:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(1,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,1,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cs
               file_type_id = file_type_id_cs

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_cs_1_


            subroutine write_dset_cd_1_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(in)   :: data(:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(1,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,1,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cd
               file_type_id = file_type_id_cd

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_cd_1_



            subroutine write_dset_s_1_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(in)   :: data(:)
                  contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               integer(HID_T) :: dspace_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 820 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 821 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5tcopy_f(H5T_NATIVE_CHARACTER,file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 823 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(file_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 824 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5screate_simple_f(1,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 827 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 834 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif
   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 835 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 836 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 838 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 840 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif
   call h5tclose_f(file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 841 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_s_1_



            subroutine write_dset_l_1_ (self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(in)   :: data(:)
                  contiguous :: data

               ! Write the logical item

               call self%write_dset(item_name, MERGE(1, 0, MASK=data))

               ! Finish

               return

            end subroutine write_dset_l_1_




            subroutine write_dset_is_2_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(in)   :: data(:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(2,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,2,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_is
               file_type_id = file_type_id_is

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_is_2_


            subroutine write_dset_id_2_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(in)   :: data(:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(2,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,2,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_id
               file_type_id = file_type_id_id

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_id_2_


            subroutine write_dset_rs_2_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(in)   :: data(:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(2,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,2,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rs
               file_type_id = file_type_id_rs

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_rs_2_


            subroutine write_dset_rd_2_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(in)   :: data(:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(2,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,2,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rd
               file_type_id = file_type_id_rd

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_rd_2_


            subroutine write_dset_cs_2_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(in)   :: data(:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(2,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,2,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cs
               file_type_id = file_type_id_cs

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_cs_2_


            subroutine write_dset_cd_2_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(in)   :: data(:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(2,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,2,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cd
               file_type_id = file_type_id_cd

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_cd_2_



            subroutine write_dset_s_2_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(in)   :: data(:,:)
                  contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               integer(HID_T) :: dspace_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 820 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 821 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5tcopy_f(H5T_NATIVE_CHARACTER,file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 823 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(file_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 824 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5screate_simple_f(2,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 827 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 834 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif
   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 835 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 836 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 838 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 840 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif
   call h5tclose_f(file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 841 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_s_2_



            subroutine write_dset_l_2_ (self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(in)   :: data(:,:)
                  contiguous :: data

               ! Write the logical item

               call self%write_dset(item_name, MERGE(1, 0, MASK=data))

               ! Finish

               return

            end subroutine write_dset_l_2_




            subroutine write_dset_is_3_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(in)   :: data(:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(3,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,3,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_is
               file_type_id = file_type_id_is

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_is_3_


            subroutine write_dset_id_3_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(in)   :: data(:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(3,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,3,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_id
               file_type_id = file_type_id_id

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_id_3_


            subroutine write_dset_rs_3_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(in)   :: data(:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(3,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,3,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rs
               file_type_id = file_type_id_rs

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_rs_3_


            subroutine write_dset_rd_3_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(in)   :: data(:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(3,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,3,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rd
               file_type_id = file_type_id_rd

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_rd_3_


            subroutine write_dset_cs_3_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(in)   :: data(:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(3,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,3,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cs
               file_type_id = file_type_id_cs

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_cs_3_


            subroutine write_dset_cd_3_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(in)   :: data(:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(3,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,3,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cd
               file_type_id = file_type_id_cd

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_cd_3_



            subroutine write_dset_s_3_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(in)   :: data(:,:,:)
                  contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               integer(HID_T) :: dspace_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 820 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 821 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5tcopy_f(H5T_NATIVE_CHARACTER,file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 823 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(file_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 824 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5screate_simple_f(3,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 827 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 834 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif
   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 835 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 836 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 838 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 840 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif
   call h5tclose_f(file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 841 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_s_3_



            subroutine write_dset_l_3_ (self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(in)   :: data(:,:,:)
                  contiguous :: data

               ! Write the logical item

               call self%write_dset(item_name, MERGE(1, 0, MASK=data))

               ! Finish

               return

            end subroutine write_dset_l_3_




            subroutine write_dset_is_4_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(in)   :: data(:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(4,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,4,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_is
               file_type_id = file_type_id_is

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_is_4_


            subroutine write_dset_id_4_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(in)   :: data(:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(4,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,4,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_id
               file_type_id = file_type_id_id

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_id_4_


            subroutine write_dset_rs_4_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(in)   :: data(:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(4,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,4,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rs
               file_type_id = file_type_id_rs

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_rs_4_


            subroutine write_dset_rd_4_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(in)   :: data(:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(4,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,4,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rd
               file_type_id = file_type_id_rd

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_rd_4_


            subroutine write_dset_cs_4_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(in)   :: data(:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(4,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,4,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cs
               file_type_id = file_type_id_cs

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_cs_4_


            subroutine write_dset_cd_4_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(in)   :: data(:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(4,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,4,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cd
               file_type_id = file_type_id_cd

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_cd_4_



            subroutine write_dset_s_4_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(in)   :: data(:,:,:,:)
                  contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               integer(HID_T) :: dspace_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 820 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 821 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5tcopy_f(H5T_NATIVE_CHARACTER,file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 823 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(file_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 824 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5screate_simple_f(4,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 827 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 834 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif
   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 835 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 836 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 838 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 840 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif
   call h5tclose_f(file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 841 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_s_4_



            subroutine write_dset_l_4_ (self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(in)   :: data(:,:,:,:)
                  contiguous :: data

               ! Write the logical item

               call self%write_dset(item_name, MERGE(1, 0, MASK=data))

               ! Finish

               return

            end subroutine write_dset_l_4_




            subroutine write_dset_is_5_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(in)   :: data(:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(5,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,5,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_is
               file_type_id = file_type_id_is

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_is_5_


            subroutine write_dset_id_5_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(in)   :: data(:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(5,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,5,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_id
               file_type_id = file_type_id_id

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_id_5_


            subroutine write_dset_rs_5_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(in)   :: data(:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(5,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,5,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rs
               file_type_id = file_type_id_rs

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_rs_5_


            subroutine write_dset_rd_5_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(in)   :: data(:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(5,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,5,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rd
               file_type_id = file_type_id_rd

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_rd_5_


            subroutine write_dset_cs_5_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(in)   :: data(:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(5,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,5,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cs
               file_type_id = file_type_id_cs

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_cs_5_


            subroutine write_dset_cd_5_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(in)   :: data(:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(5,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,5,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cd
               file_type_id = file_type_id_cd

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_cd_5_



            subroutine write_dset_s_5_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(in)   :: data(:,:,:,:,:)
                  contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               integer(HID_T) :: dspace_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 820 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 821 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5tcopy_f(H5T_NATIVE_CHARACTER,file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 823 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(file_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 824 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5screate_simple_f(5,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 827 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 834 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif
   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 835 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 836 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 838 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 840 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif
   call h5tclose_f(file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 841 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_s_5_



            subroutine write_dset_l_5_ (self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(in)   :: data(:,:,:,:,:)
                  contiguous :: data

               ! Write the logical item

               call self%write_dset(item_name, MERGE(1, 0, MASK=data))

               ! Finish

               return

            end subroutine write_dset_l_5_




            subroutine write_dset_is_6_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(in)   :: data(:,:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(6,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,6,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_is
               file_type_id = file_type_id_is

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_is_6_


            subroutine write_dset_id_6_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(in)   :: data(:,:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(6,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,6,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_id
               file_type_id = file_type_id_id

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_id_6_


            subroutine write_dset_rs_6_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(in)   :: data(:,:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(6,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,6,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rs
               file_type_id = file_type_id_rs

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_rs_6_


            subroutine write_dset_rd_6_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(in)   :: data(:,:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(6,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,6,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rd
               file_type_id = file_type_id_rd

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_rd_6_


            subroutine write_dset_cs_6_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(in)   :: data(:,:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(6,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,6,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cs
               file_type_id = file_type_id_cs

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_cs_6_


            subroutine write_dset_cd_6_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(in)   :: data(:,:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(6,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,6,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cd
               file_type_id = file_type_id_cd

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_cd_6_



            subroutine write_dset_s_6_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(in)   :: data(:,:,:,:,:,:)
                  contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               integer(HID_T) :: dspace_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 820 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 821 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5tcopy_f(H5T_NATIVE_CHARACTER,file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 823 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(file_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 824 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5screate_simple_f(6,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 827 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 834 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif
   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 835 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 836 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 838 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 840 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif
   call h5tclose_f(file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 841 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_s_6_



            subroutine write_dset_l_6_ (self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(in)   :: data(:,:,:,:,:,:)
                  contiguous :: data

               ! Write the logical item

               call self%write_dset(item_name, MERGE(1, 0, MASK=data))

               ! Finish

               return

            end subroutine write_dset_l_6_




            subroutine write_dset_is_7_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(IS), target, intent(in)   :: data(:,:,:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(7,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,7,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_is
               file_type_id = file_type_id_is

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_is_7_


            subroutine write_dset_id_7_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               integer(ID), target, intent(in)   :: data(:,:,:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(7,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,7,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_id
               file_type_id = file_type_id_id

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_id_7_


            subroutine write_dset_rs_7_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RS), target, intent(in)   :: data(:,:,:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(7,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,7,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rs
               file_type_id = file_type_id_rs

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_rs_7_


            subroutine write_dset_rd_7_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               real(RD), target, intent(in)   :: data(:,:,:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(7,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,7,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_rd
               file_type_id = file_type_id_rd

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_rd_7_


            subroutine write_dset_cs_7_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RS), target, intent(in)   :: data(:,:,:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(7,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,7,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cs
               file_type_id = file_type_id_cs

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_cs_7_


            subroutine write_dset_cd_7_(self, item_name, data)
            
               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               complex(RD), target, intent(in)   :: data(:,:,:,:,:,:,:)
               contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: dspace_id
               integer(HID_T) :: plist_id
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the item

   call h5screate_simple_f(7,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 761 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 769 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pcreate_f'
   error stop
   endif
   call h5pset_chunk_f(plist_id,7,INT(SHAPE(data), HSIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 771 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_chunk_f'
   error stop
   endif
   call h5pset_deflate_f(plist_id,self%comp_level,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 772 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pset_deflate_f'
   error stop
   endif

               mem_type_id = mem_type_id_cd
               file_type_id = file_type_id_cd

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err,dcpl_id=plist_id)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 784 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif

   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 787 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5pclose_f(plist_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 788 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5pclose_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 789 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif
   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 790 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_cd_7_



            subroutine write_dset_s_7_(self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               character(*), target, intent(in)   :: data(:,:,:,:,:,:,:)
                  contiguous :: data

               integer        :: hdf_err
               integer(HID_T) :: mem_type_id
               integer(HID_T) :: file_type_id
               integer(HID_T) :: dspace_id
               type(C_PTR)    :: data_ptr
               integer(HID_T) :: item_id

               ! Write the character item

   call h5tcopy_f(H5T_NATIVE_CHARACTER,mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 820 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(mem_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 821 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5tcopy_f(H5T_NATIVE_CHARACTER,file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 823 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tcopy_f'
   error stop
   endif
   call h5tset_size_f(file_type_id,LEN(data, SIZE_T),hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 824 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tset_size_f'
   error stop
   endif

   call h5screate_simple_f(7,INT(SHAPE(data), HSIZE_T),dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 827 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5screate_simple_f'
   error stop
   endif

               data_ptr = C_LOC(data)

   call h5dcreate_f(self%group_id,item_name,file_type_id,dspace_id,item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 834 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dcreate_f'
   error stop
   endif
   call h5dwrite_f(item_id,mem_type_id,data_ptr,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 835 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dwrite_f'
   error stop
   endif
   call h5dclose_f(item_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 836 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5dclose_f'
   error stop
   endif

   call h5sclose_f(dspace_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 838 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5sclose_f'
   error stop
   endif

   call h5tclose_f(mem_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 840 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif
   call h5tclose_f(file_type_id,hdf_err)
   if (hdf_err == -1) then
      call h5eprint_f(hdf_err)
   write(UNIT=ERROR_UNIT, FMT=*) 'ABORT at line 841 of hdf5io_m.fypp'
      write(UNIT=ERROR_UNIT, FMT=*) 'error in call to h5tclose_f'
   error stop
   endif

               ! Finish

               return

            end subroutine write_dset_s_7_



            subroutine write_dset_l_7_ (self, item_name, data)

               class(hdf5io_t), intent(inout) :: self
               character(*), intent(in)       :: item_name
               logical, target, intent(in)   :: data(:,:,:,:,:,:,:)
                  contiguous :: data

               ! Write the logical item

               call self%write_dset(item_name, MERGE(1, 0, MASK=data))

               ! Finish

               return

            end subroutine write_dset_l_7_


            
end module hdf5io_m
