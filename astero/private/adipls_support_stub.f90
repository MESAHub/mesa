! ***********************************************************************
!
!   Copyright (C) 2021  The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************


module adipls_support
   
   use astero_def
   use star_lib
   use star_def
   use const_def
   use utils_lib
   
   implicit none
   
   
   ! args for adipls
   integer :: i_paramset, ierr_param, i_inout, nn
   real(dp), pointer :: x(:) ! (nn)
   real(dp), pointer :: aa(:, :) ! (iaa_arg,nn)
   real(dp) :: data(8)
   
   integer, parameter :: ivarmd = 6, iaa_arg = 10
   
   integer :: iounit_dev_null = -1
   
   integer :: nn_redist ! set from redistrb.c input file
   
   real(dp), pointer :: x_arg(:), aa_arg(:, :)
   integer :: nn_arg
   real(dp) :: data_arg(8)
   
   logical, parameter :: ADIPLS_IS_ENABLED = .false.

contains
   
   
   ! this can be called from user run_star_extras check model routine
   subroutine do_adipls_get_one_el_info(&
      s, l, nu1, nu2, iscan, R, G, M, &
      add_center_point, keep_surface_point, add_atmosphere, &
      redist_mesh, store_for_adipls, &
      save_mode_info, order_to_save_in, save_mode_filename_in, &
      num, l_freq, l_inertia, l_order, l_em, ierr)
      use num_lib, only : qsort
      use star_def, only : star_info
      type (star_info), pointer :: s
      integer, intent(in) :: l, iscan
      real(dp), intent(in) :: nu1, nu2, R, G, M
      logical, intent(in) :: &
         add_center_point, keep_surface_point, add_atmosphere, &
         redist_mesh, store_for_adipls, save_mode_info
      integer, intent(in) :: order_to_save_in
      character (len = *), intent(in) :: save_mode_filename_in
      integer, intent(out) :: num
      real(dp), pointer, dimension(:) :: l_freq, l_inertia
      integer, pointer, dimension(:) :: l_order, l_em
      integer, intent(out) :: ierr
      
      real(dp) :: sig_fac
      integer :: nsel, itrsig, nsig, irotkr, nprtkr, igm1kr, npgmkr
      real(dp) :: els1, dels, sig1, sig2, dfsig
      integer :: k, i, j
      integer, pointer :: index(:)
      
      logical, parameter :: dbg = .false.
      
      include 'formats'
      
      ierr = -1
   
   end subroutine do_adipls_get_one_el_info
   
   
   subroutine adipls_mode_info(&
      l, order, em, freq, inertia, x, y, aa, data, nn, iy, iaa, ispcpr)
      integer, intent(in) :: l, order, em
      real(dp), intent(in) :: freq, inertia
      integer, intent(in) :: nn, iy, iaa, ispcpr
      real(dp), intent(in) :: x(1:nn), y(1:iy, 1:nn), aa(1:iaa, 1:nn), data(8)
      integer :: iounit, ierr, i, j, skip
      real(dp) :: y_r, y_h
   
   end subroutine adipls_mode_info
   
   
   subroutine store_model_for_adipls (s, add_atmosphere, do_redistribute_mesh, ierr)
      
      type (star_info), pointer :: s
      logical, intent(in) :: add_atmosphere
      logical, intent(in) :: do_redistribute_mesh
      integer, intent(out) :: ierr
      
      common/ccgrav/ cgrav
      real(dp) :: cgrav
      
      integer :: i, iriche, iturpr
      integer :: iconst, ivar, ivers
      real(dp), allocatable :: global_data(:) ! (iconst)
      real(dp), allocatable :: point_data(:, :) ! (ivar,nn)
      character (len = 2000) :: format_string, num_string, filename
      
      ierr = -1
   
   end subroutine store_model_for_adipls
   
   
   subroutine run_adipls(&
      s, first_time, store_model, &
      add_center_point, keep_surface_point, add_atmosphere, &
      do_redistribute_mesh, ierr)
      type (star_info), pointer :: s
      logical, intent(in) :: &
         first_time, store_model, &
         add_center_point, keep_surface_point, add_atmosphere, &
         do_redistribute_mesh
      integer, intent(out) :: ierr
      
      integer :: iounit, nn_arg_0
      integer(8) :: time0, time1, clock_rate
      real(dp) :: time, x_arg0(0), aa_arg0(0, 0)
      character (len = 256) :: filename
      common/cstdio/ istdin, istdou, istdpr, istder
      integer :: istdin, istdou, istdpr, istder
      
      logical, parameter :: dbg = .false.
      
      include 'formats'
      
      ierr = -1
   
   end subroutine run_adipls
   
   
   subroutine set_adipls_controls(&
      el, nsel, els1, dels, itrsig, iscan, sig1, sig2, dfsig, nsig, &
      irotkr, nprtkr, igm1kr, npgmkr)
      integer, intent(in) :: el, nsel, itrsig, iscan, nsig, &
         irotkr, nprtkr, igm1kr, npgmkr
      real(dp), intent(in) :: els1, dels, sig1, sig2, dfsig
   
   end subroutine set_adipls_controls
   
   
   ! this is called by modmod
   subroutine check_arg_data(nn, data, ldaa, aa, x, ierr)
      integer, intent(in) :: nn, ldaa
      real(dp), intent(in) :: data(8)
      real(dp) :: aa(ldaa, nn)
      real(dp) :: x(nn)
      integer, intent(out) :: ierr
      
      real(dp), parameter :: rtol = 1d-9, atol = 1d-9
      
      integer :: i, j
      
      ierr = -1
   
   end subroutine check_arg_data
   
   
   subroutine read_and_store(iriche, iturpr, cgrav)
      integer, intent(inout) :: iriche, iturpr
      real(dp), intent(in) :: cgrav
      character (len = 64) :: fname
      integer :: nn, iconst, ivar, ivers, ierr
      real(dp), pointer :: glob(:) ! (iconst)   will be allocated
      real(dp), pointer :: var(:, :) ! (ivar,nn)   will be allocated
      real(dp), pointer :: aa(:, :) ! (iaa_arg,nn)   will be allocated
      real(dp), pointer :: x(:) ! (nn)   will be allocated
      real(dp) :: data(8)
   
   end subroutine read_and_store
   
   
   subroutine store_amdl(nn_in, iriche, iturpr, data, aa, x, nn, ierr)
      ! derived from adipls readml.n.d.f
      integer, intent(in) :: nn_in, iriche
      integer, intent(inout) :: iturpr
      real(dp), intent(in) :: data(8)
      real(dp), pointer :: aa(:, :)
      real(dp), pointer :: x(:) ! (nn)     will be allocated
      ! nn can be less than nn_in
      integer, intent(out) :: nn, ierr
      
      ! local
      integer :: i, j, nsin, iggt, inp, in, nshift, nnr, n, n1, nstart, idata8
      logical :: sincen, sinsur
      real(dp), pointer :: aa1(:, :)
      real(dp) :: ggt
      
      ierr = -1
   
   end subroutine store_amdl
   
   
   subroutine fgong_amdl(&
      cgrav, nn_in, iconst, ivar, ivers, glob, var, data, aa, nn, ierr)
      ! derived from fgong-amdl.d.f
      real(dp), intent(in) :: cgrav
      integer, intent(in) :: nn_in, iconst, ivar, ivers
      real(dp), intent(inout) :: glob(:) ! (iconst)
      real(dp), intent(inout) :: var(:, :) ! (ivar,nn_in)
      real(dp), intent(inout) :: data(8)
      real(dp), pointer :: aa(:, :) ! (iaa_arg,nn)   will be allocated
      integer, intent(out) :: nn, ierr
      
      integer, parameter :: ireset(16) = &
         (/3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20/)
      integer :: nn1, i, n, ir
      real(dp) :: d2amax, var1(ivar, nn_in + 100), q(nn_in + 100), x(nn_in + 100)
      
      ierr = -1
   
   end subroutine fgong_amdl
   
   
   subroutine read_fgong_file(fin, nn, iconst, ivar, ivers, glob, var, ierr)
      character (len = *), intent(in) :: fin
      integer, intent(out) :: nn, iconst, ivar, ivers
      real(dp), pointer :: glob(:) ! (iconst)   will be allocated
      real(dp), pointer :: var(:, :) ! (ivar,nn)   will be allocated
      integer, intent(out) :: ierr
      
      real(dp), pointer :: var1(:, :) ! (ivar,nn)
      integer :: ios, iounit, i, n, ir, nn1
      character(80) :: head
      
      ierr = -1
   
   end subroutine read_fgong_file
   
   
   ! for testing
   subroutine dump(filename_for_dump, nn, glob, var, ierr)
      character (len = *), intent(in) :: filename_for_dump
      integer, intent(in) :: nn
      real(dp), pointer :: glob(:) ! (iconst)
      real(dp), pointer :: var(:, :) ! (ivar,nn)
      integer, intent(out) :: ierr
      
      ierr = -1
   end subroutine dump
   
   
   subroutine show_adipls_results
      integer :: k
   
   end subroutine show_adipls_results


end module adipls_support
