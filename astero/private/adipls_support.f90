! ***********************************************************************
!
!   Copyright (C) 2013  Bill Paxton
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
      real(dp), pointer :: aa(:,:) ! (iaa_arg,nn)
      real(dp) :: data(8)

      integer, parameter :: ivarmd = 6, iaa_arg = 10

      integer :: iounit_dev_null = -1

      integer :: nn_redist ! set from redistrb.c input file
      
      
      real(dp), pointer :: x_arg(:), aa_arg(:,:)
      integer :: nn_arg
      real(dp) :: data_arg(8)
      

      
      contains


      ! this can be called from user run_star_extras check model routine
      subroutine do_adipls_get_one_el_info( &
            s, l, nu1, nu2, iscan, R, G, M, &
            add_center_point, keep_surface_point, add_atmosphere, &
            redist_mesh, store_for_adipls, &
            save_mode_info, order_to_save_in, save_mode_filename_in, &
            num, l_freq, l_inertia, l_order, l_em, ierr)
         use num_lib, only: qsort
         use star_def, only: star_info
         type (star_info), pointer :: s
         integer, intent(in) :: l, iscan
         real(dp), intent(in) :: nu1, nu2, R, G, M
         logical, intent(in) :: &
            add_center_point, keep_surface_point, add_atmosphere, &
            redist_mesh, store_for_adipls, save_mode_info
         integer, intent(in) :: order_to_save_in
         character (len=*), intent(in) :: save_mode_filename_in
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
         
         ierr = 0
         sig_fac = (2*pi)**2*R**3/(G*M)
         nsel = 0
         dels = 1
         els1 = dble(l)
         itrsig = 1
         sig1 = sig_fac*(nu1*1d-6)*(nu1*1d-6)
         sig2 = sig_fac*(nu2*1d-6)*(nu2*1d-6)
         nsig = 2
         dfsig = sig_fac*delta_nu_model*delta_nu_model
         
         if (dbg) write(*,*) 'call set_adipls_controls'
         call set_adipls_controls( &
            l, nsel, els1, dels, itrsig, iscan, sig1, sig2, dfsig, nsig, &
            adipls_irotkr, adipls_nprtkr, adipls_igm1kr, adipls_npgmkr)
         if (dbg) write(*,*) 'done set_adipls_controls'

         el_to_save = l
         order_to_save = order_to_save_in
         save_mode_filename = save_mode_filename_in
      
         num_results = 0
         if (dbg) write(*,*) 'call run_adipls'
         call run_adipls(s, .false., store_for_adipls, &
            add_center_point, keep_surface_point, add_atmosphere, &
            redist_mesh, ierr)
         if (dbg) write(*,*) 'done run_adipls'
         if (ierr /= 0) then
            write(*,*) 'failed in run_adipls'
            return
         end if
         num = num_results
         
         if (num_results == 0) then
            write(*,*) 'failed to find any modes in specified frequency range'
            return
         end if
         
         ! sort results by increasing frequency
         allocate(index(num_results), stat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in allocate before calling qsort'
            return
         end if
         call qsort(index, num_results, cyclic_freq)

         if (.not. associated(l_order)) then
            allocate(l_order(num_results))
         else if (num_results >= size(l_order,dim=1)) then ! enlarge
            call realloc_integer(l_order,num_results,ierr)
            if (ierr /= 0) return
         end if

         if (.not. associated(l_em)) then
            allocate(l_em(num_results))
         else if (num_results >= size(l_em,dim=1)) then ! enlarge
            call realloc_integer(l_em,num_results,ierr)
            if (ierr /= 0) return
         end if

         if (.not. associated(l_freq)) then
            allocate(l_freq(num_results))
         else if (num_results >= size(l_freq,dim=1)) then ! enlarge
            call realloc_double(l_freq,num_results,ierr)
            if (ierr /= 0) return
         end if

         if (.not. associated(l_inertia)) then
            allocate(l_inertia(num_results))
         else if (num_results >= size(l_inertia,dim=1)) then ! enlarge
            call realloc_double(l_inertia,num_results,ierr)
            if (ierr /= 0) return
         end if
         
         do j = 1, num_results
            i = index(j)
            l_freq(j) = cyclic_freq(i)
            l_inertia(j) = inertia(i)
            l_order(j) = order(i)
            l_em(j) = em(i)
         end do
         
         deallocate(index)
         
      end subroutine do_adipls_get_one_el_info

      
      subroutine adipls_mode_info( &
            l, order, em, freq, inertia, x, y, aa, data, nn, iy, iaa, ispcpr)
         integer, intent(in) :: l, order, em
         real(dp), intent(in) :: freq, inertia
         integer, intent(in) :: nn, iy, iaa, ispcpr
         real(dp), intent(in) :: x(1:nn), y(1:iy,1:nn), aa(1:iaa,1:nn), data(8)
         integer :: iounit, ierr, i, j, skip
         real(dp) :: y_r, y_h
         include 'formats'
         if (use_other_adipls_mode_info) then
            call astero_other_procs% other_adipls_mode_info( &
               l, order, freq, inertia, x, y, aa, data, nn, iy, iaa, ispcpr, ierr)
         end if
         if (star_model_number /= save_mode_model_number) return
         !write(*,3) 'adipls_mode_info l order freq inertia', l, order, freq, inertia
         if (l /= el_to_save .or. order /= order_to_save) return
         if (len_trim(save_mode_filename) <= 0) save_mode_filename = 'save_mode.data'
         write(*,*) 'save eigenfunction info to file ' // trim(save_mode_filename)
         write(*,'(3a8,99a20)') 'el', 'order', 'em', 'freq (microHz)', 'inertia'
         write(*,'(3i8,f20.10,e20.10,i20)') l, order, em, freq, inertia
         ierr = 0
         iounit = alloc_iounit(ierr)
         if (ierr /= 0) return
         open(unit=iounit, file=trim(save_mode_filename), action='write', iostat=ierr)
         if (ierr /= 0) return
         if (abs(x(1)) < 1d-20) then
            skip = 1
         else
            skip = 0
         end if
         write(iounit,'(3a8,99a20)') 'el', 'order', 'em', 'freq (microHz)', 'inertia', 'nn'
         write(iounit,'(3i8,f20.10,e20.10,i20)') l, order, em, freq, inertia, nn-skip
         write(iounit,'(a)') &
            'x = r/R;  y_r = xi_r/R;  y_h = xi_h*l*(l+1)/R;' // &
            '  displacement xi normalized to xi_r = R at surface.'
         write(iounit,'(a6,4a26)') 'i', 'x', 'y_r', 'y_h'
         do i = 1+skip, nn
            y_r = y(1,i)
            if (l > 0) then
               y_h = y(2,i)
            else
               y_h = 0
            end if
            write(iounit,'(i6,4e26.16)') i-skip, x(i), y_r, y_h
         end do
         close(iounit)
         call free_iounit(iounit)         
      end subroutine adipls_mode_info
      
      
      subroutine store_model_for_adipls (s, add_atmosphere, do_redistribute_mesh, ierr)

        type (star_info), pointer :: s
        logical, intent(in)       :: add_atmosphere
        logical, intent(in)       ::do_redistribute_mesh
         integer, intent(out)     :: ierr

         common/ccgrav/ cgrav
         real(dp) :: cgrav

         integer :: i, iriche, iturpr
         integer :: iconst, ivar, ivers
         real(dp), allocatable :: global_data(:) ! (iconst)
         real(dp), allocatable :: point_data(:,:) ! (ivar,nn)
         character (len=2000) :: format_string, num_string, filename
         
         ierr = 0
         iriche = 0
         iturpr = 0
         
         if (associated(x)) deallocate(x)
         if (associated(aa)) deallocate(aa)

         ! Get the model data

         call star_get_pulse_data(s%id, 'FGONG', &
              add_center_point, keep_surface_point, add_atmosphere, global_data, point_data, ierr)

         if (ierr /= 0) then
            write(*,*) 'failed in star_get_pulse_data'
            call mesa_error(__FILE__,__LINE__)
         end if

         ! If necessary, write it
         
         if (write_fgong_for_each_model) then

            write(format_string,'( "(i",i2.2,".",i2.2,")" )') &
               model_num_digits, model_num_digits
            write(num_string,format_string) s% model_number

            filename = trim(fgong_prefix) // trim(num_string) // trim(fgong_postfix)

            call star_write_pulse_data(s%id, 'FGONG', filename, global_data, point_data, ierr)

            if (ierr /= 0) then
               write(*,*) 'failed in star_write_pulse_data'
               call mesa_error(__FILE__,__LINE__)
            end if

         end if

         ! Convert to AMDL format

         iconst = SIZE(global_data)
         ivar = SIZE(point_data, 1)
         nn = SIZE(point_data, 2)

         ivers = 0 ! It's not clear what this does in fgong_amdl

         call fgong_amdl( &
            cgrav, nn, iconst, ivar, ivers, global_data, point_data, data, aa, nn, ierr)
         deallocate(global_data, point_data)

         if (ierr /= 0) then
            write(*,*) 'failed in fgong_amdl'
            call mesa_error(__FILE__,__LINE__)
         end if
        
         call store_amdl(nn, iriche, iturpr, data, aa, x, nn, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in store_amdl'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         call redist_amdl(ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in redist_amdl'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         
         contains
         
         
         subroutine redist_amdl(ierr)
            integer, intent(out) :: ierr
            real(dp), pointer :: aa_new(:,:)
            real(dp), pointer :: x_new(:)
            integer :: nn_new
            include 'formats'
            ierr = 0
            if (.not. do_redistribute_mesh) return
            nn_new = nn_redist ! srdist uses nn from input file
            allocate(aa_new(iaa_arg,nn_new), x_new(nn_new))
            ierr_param = 0
            !write(*,2) 'call srdist: nn_redist', nn_new
            call srdist(i_paramset, ierr_param, i_inout, &
               x, aa, data, x_new, aa_new, &
               nn, nn_new, ivarmd, iaa_arg, iaa_arg)
            !write(*,1) 'done srdist'
            deallocate(aa, x)
            aa => aa_new
            x => x_new
            nn = nn_new
            if (ierr_param < 0) ierr = -1
         end subroutine redist_amdl
      
      
      end subroutine store_model_for_adipls
      
      
      subroutine run_adipls( &
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
         real(dp) :: time, x_arg0(0), aa_arg0(0,0)
         character (len=256) :: filename
         common/cstdio/ istdin, istdou, istdpr, istder
         integer :: istdin, istdou, istdpr, istder
         
         logical, parameter :: dbg = .false.
         
         include 'formats'
         
         ierr = 0
         
         i_inout = 0
         i_paramset = 1
         ierr_param = 0

         if (first_time) then
            nullify(aa, x)
            call setup_redist
            call setup_adipls
            return
         end if
         
         if (iounit_dev_null > 0) then
            close(iounit_dev_null)
         else
            iounit_dev_null = alloc_iounit(ierr)
            if (ierr /= 0) then
               write(*,*) 'adipls failed in alloc_iounit for iounit_dev_null'
               stop 'run_adipls'
            end if
         end if
         
         filename = 'adipls.stdout'
         open(unit=iounit_dev_null, file=trim(filename), iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'adipls failed to open ' // trim(filename)
            stop 'run_adipls'
         end if            
         istdou = iounit_dev_null
         istdpr = iounit_dev_null
         
         if (store_model) then
            if (dbg) write(*,*) 'call store_model_for_adipls'
            call store_model_for_adipls(s, add_atmosphere, do_redistribute_mesh, ierr)
            if (dbg) write(*,*) 'done store_model_for_adipls'
            if (ierr /= 0) return
         end if
         
         ! ivarmd and iaa_arg are defined in file with store_amdl         
         if (dbg) write(*,*) 'call adipls'
         
         if (trace_time_in_oscillation_code) then
            call system_clock(time0, clock_rate)
         end if

         call adipls(i_paramset, ierr_param, i_inout, &
               x, aa, data, nn, ivarmd, iaa_arg)
         if (dbg) write(*,*) 'done adipls'
         
         if (trace_time_in_oscillation_code) then
            call system_clock(time1, clock_rate)
            time = dble(time1-time0)/clock_rate
            total_time_in_oscillation_code = total_time_in_oscillation_code + time
            write(*,1) 'time_in_oscillation_code and total', time, total_time_in_oscillation_code
         end if

         if (ierr_param < 0) then
            ierr = ierr_param
            write(*,*) 'call to adipls failed'
            stop 'run_adipls'
         end if
         
         
         contains
         
         
         subroutine setup_adipls
            iounit = alloc_iounit(ierr)
            if (ierr /= 0) then
               write(*,*) 'setup_adipls failed in alloc_iounit'
               stop 'run_adipls'
            end if
            filename = 'adipls.c.pruned.in'
            open(unit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
            if (ierr /= 0) then
               write(*,*) 
               write(*,*) 
               write(*,*) 
               write(*,*) 
               write(*,*) 'ERROR: failed to open ' // trim(filename)
               write(*,*) 'please convert adipls.c.in to "pruned" form'
               write(*,*) 'e.g., you can run the get-input script from mesa/adipls/adipack.c/bin:'
               write(*,*) './../../adipls/adipack.c/bin/get-input adipls.c.in > adipls.c.pruned.in'
               write(*,*) 
               write(*,*) 
               write(*,*) 
               write(*,*) 
               stop 'run_adipls'
            end if            
            
            write(*,*)
            write(*,'(a)') 'call adipls to read ' // trim(filename)
            call setups_adi
            nn_arg_0 = 0
            istdin = iounit
            call adipls(i_paramset, ierr_param, i_inout, &
                  x_arg0, aa_arg0, data_arg, nn_arg_0, ivarmd, iaa_arg)
            close(iounit)
            
            call free_iounit(iounit)
            
            if (ierr_param < 0) then
               ierr = ierr_param
               write(*,*) '1st call to adipls failed in setup_adipls'
               stop 'run_adipls'
            end if
            
            write(*,*) 'back from 1st call on adipls'
            write(*,*)

         end subroutine setup_adipls

         
         subroutine setup_redist

            common/comgrp/ isprtp
            integer :: isprtp

            if (.not. do_redistribute_mesh) return
            
            iounit = alloc_iounit(ierr)
            if (ierr /= 0) then
               write(*,*) 'setup_redist failed in alloc_iounit'
               stop 'run_rdist'
            end if
            filename = 'redistrb.c.pruned.in'
            open(unit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
            if (ierr /= 0) then
               write(*,*) 'setup_redist failed to open ' // trim(filename)
               write(*,*) 
               write(*,*) 
               write(*,*) 
               write(*,*) 
               write(*,*) 'ERROR: failed to open ' // trim(filename)
               write(*,*) 'please convert redistrb.c.in to "pruned" form'
               write(*,*) 'e.g., you can run the get-input script from mesa/adipls/adipack.c/bin:'
               write(*,*) './../../adipls/adipack.c/bin/get-input redistrb.c.in > redistrb.c.pruned.in'
               write(*,*) 
               write(*,*) 
               write(*,*) 
               write(*,*) 
               stop 'run_adipls'
               stop 'run_rdist'
            end if 
            
            read(iounit,*,iostat=ierr) nn_redist
            if (ierr /= 0) then
               write(*,*) 'setup_redist failed to read nn_redist from ' // trim(filename)
               stop 'run_rdist'
            end if 
            write(*,*) 'nn_redist', nn_redist
            
            rewind(iounit)
                       
            write(*,*)
            write(*,'(a)') 'call srdist to read ' // trim(filename)
            
            istdin = iounit
            i_inout = 0
            i_paramset = 1
            ierr_param = 0
            isprtp = 0
            
            call srdist(i_paramset, ierr_param, i_inout, &
               x_arg0, aa_arg0, data_arg, x_arg0, aa_arg0, &
               nn_arg_0, nn_arg_0, ivarmd, iaa_arg, iaa_arg)
                  
            close(iounit)
            
            call free_iounit(iounit)
            
            if (ierr_param < 0) then
               ierr = ierr_param
               write(*,*) '1st call to srdist failed'
               stop 'run_rdist'
            end if
            
            write(*,*) 'back from 1st call on srdist'
            write(*,*)

         end subroutine setup_redist


      end subroutine run_adipls
      
      
      subroutine set_adipls_controls( &
            el, nsel, els1, dels, itrsig, iscan, sig1, sig2, dfsig, nsig, &
            irotkr, nprtkr, igm1kr, npgmkr)
         integer, intent(in) :: el, nsel, itrsig, iscan, nsig, &
            irotkr, nprtkr, igm1kr, npgmkr
         real(dp), intent(in) :: els1, dels, sig1, sig2, dfsig

         !  commons for parameter transmission
         common/cadr_param/ &
            para_el, para_els1, para_dels, para_dfsig1, para_dfsig2, &
            para_sig1, para_sig2, para_dfsig, para_eltrw1, para_eltrw2, &
            para_sgtrw1, para_sgtrw2
         real(dp) :: &
            para_el, para_els1, para_dels, para_dfsig1, para_dfsig2, &
            para_sig1, para_sig2, para_dfsig, para_eltrw1, para_eltrw2, &
            para_sgtrw1, para_sgtrw2
            
         common/cadi_param/ &
            ipara_nsel, ipara_nsig1, ipara_nsig2, ipara_itrsig, ipara_nsig, &
            ipara_istsig, ipara_inomd1, ipara_iscan
         integer :: &
            ipara_nsel, ipara_nsig1, ipara_nsig2, ipara_itrsig, ipara_nsig, &
            ipara_istsig, ipara_inomd1, ipara_iscan
            
         common/coutpt/ &
            ipara_nout, ipara_nprcen, ipara_iper, ipara_irotkr, ipara_nprtkr, &
            ipara_igm1kr, ipara_npgmkr, ipara_nfmode, ipara_nfmesh, ipara_ispcpr, &
            ipara_npout, ipara_nobs_stmx, ipara_nfmscn
         integer :: &
            ipara_nout, ipara_nprcen, ipara_iper, ipara_irotkr, ipara_nprtkr, &
            ipara_igm1kr, ipara_npgmkr, ipara_nfmode, ipara_nfmesh, ipara_ispcpr, &
            ipara_npout, ipara_nobs_stmx, ipara_nfmscn
         
         para_el = dble(el)
         ipara_nsel = nsel
         para_els1 = els1
         para_dels = dels
         ipara_itrsig = itrsig
         ipara_iscan = iscan
         para_sig1 = sig1
         para_sig2 = sig2
         para_dfsig = dfsig
         ipara_nsig = nsig
         ipara_irotkr = irotkr
         ipara_nprtkr = nprtkr
         ipara_igm1kr = igm1kr
         ipara_npgmkr = npgmkr
      
      end subroutine set_adipls_controls
      
      
      ! this is called by modmod
      subroutine check_arg_data(nn, data, ldaa, aa, x, ierr)
         integer, intent(in) :: nn, ldaa
         real(dp), intent(in) :: data(8)
         real(dp) :: aa(ldaa,nn)
         real(dp) :: x(nn)
         integer, intent(out) :: ierr
         
         real(dp), parameter :: rtol = 1d-9, atol = 1d-9
         
         integer :: i, j
         
         ierr = 0
         
         if (ldaa /= iaa_arg) then
            write(*,*) 'ldaa /= iaa_arg', ldaa, iaa_arg
            ierr = -1
            call mesa_error(__FILE__,__LINE__)
         end if
         
         if (nn /= nn_arg) then
            write(*,*) 'nn /= nn_arg', nn, nn_arg
            ierr = -1
            call mesa_error(__FILE__,__LINE__)
         end if
         
         do i=1,8
            if (is_bad(data(i),data_arg(i))) then
               write(*,'(a40,i6,99e26.16)') 'data(i) /= data_arg(i)', i, data(i), data_arg(i)
               ierr = -1
               call mesa_error(__FILE__,__LINE__)
            end if
         end do
         
         do j=1,nn
            if (is_bad(x(j),x_arg(j))) then
               write(*,'(a40,i6,99e26.16)') 'x(j) /= x_arg(j)', j, x(j), x_arg(j)
               ierr = -1
               call mesa_error(__FILE__,__LINE__)
            end if
            do i=1,iaa_arg
               if (is_bad(aa(i,j),aa_arg(i,j))) then
                  write(*,'(a40,2i6,99e26.16)') 'aa(i,j) /= aa_arg(i,j)', i, j, aa(i,j), aa_arg(i,j)
                  call mesa_error(__FILE__,__LINE__)
                  ierr = -1
               end if
            end do
         end do
         
         if (ierr /= 0) stop 'check_arg_data'
         
         
         contains
         
         logical function is_bad(v1,v2)
            real(dp), intent(in) :: v1, v2
            real(dp) :: err
            err = abs(v1-v2)/(atol + rtol*max(abs(v1),abs(v2)))
            is_bad = (err > 1d0)
         end function is_bad
         
      
      end subroutine check_arg_data
      
      
      subroutine read_and_store(iriche, iturpr, cgrav)
         integer, intent(inout) :: iriche, iturpr
         real(dp), intent(in) :: cgrav
         character (len=64) :: fname
         integer :: nn, iconst, ivar, ivers, ierr
         real(dp), pointer :: glob(:) ! (iconst)   will be allocated
         real(dp), pointer :: var(:,:) ! (ivar,nn)   will be allocated
         real(dp), pointer :: aa(:,:) ! (iaa_arg,nn)   will be allocated
         real(dp), pointer :: x(:) ! (nn)   will be allocated
         real(dp) :: data(8)
         
         ierr = 0
         fname = 'test.fgong'
         call read_fgong_file(fname, nn, iconst, ivar, ivers, glob, var, ierr)
         if (ierr /= 0) then
            write(*,*) 'read_and_store failed in read_fgong_file'
            call mesa_error(__FILE__,__LINE__)
         end if
         call fgong_amdl( &
            cgrav, nn, iconst, ivar, ivers, glob, var, data, aa, nn, ierr)
         if (ierr /= 0) then
            write(*,*) 'read_and_store failed in fgong_amdl'
            call mesa_error(__FILE__,__LINE__)
         end if
         deallocate(glob, var)
         call store_amdl(nn, iriche, iturpr, data, aa, x, nn, ierr)
         if (ierr /= 0) then
            write(*,*) 'read_and_store failed in store_amdl'
            call mesa_error(__FILE__,__LINE__)
         end if
         
      end subroutine read_and_store
      
      
      subroutine store_amdl(nn_in, iriche, iturpr, data, aa, x, nn, ierr)
         ! derived from adipls readml.n.d.f
         integer, intent(in) :: nn_in, iriche
         integer, intent(inout) :: iturpr
         real(dp), intent(in) :: data(8)
         real(dp), pointer :: aa(:,:)
         real(dp), pointer :: x(:) ! (nn)     will be allocated
         ! nn can be less than nn_in
         integer, intent(out) :: nn, ierr
         
         ! local
         integer :: i, j, nsin, iggt, inp, in, nshift, nnr, n, n1, nstart, idata8
         logical :: sincen, sinsur
         real(dp), pointer :: aa1(:,:)
         real(dp) :: ggt
         
         ierr = 0
         nn = nn_in
         
         allocate(aa1(iaa_arg,nn))
         do i=1,nn
            do j=1,ivarmd
               aa1(j,i) = aa(j,i)
            end do
         end do
      
         ! test for singular centre and/or surface

         sincen=aa1(1,1).eq.0
         sinsur=data(7).ge.0
         nsin=0
         if (sincen) nsin=nsin+1
         if (sinsur) nsin=nsin+1

         ! test for inclusion of g/(g tilde)

         idata8 = int(data(8)+0.1)
         if (mod(idata8/10,10).eq.2) then
            iggt = 1
            iturpr=8
         else
            iggt=0
         end if

         ! we always take every point in model
         
         ! test for number of nonsingular points

         if (iriche.ne.1.or.mod(nn-nsin,2).eq.1) then
            nshift=0
         else
            nshift=1
         end if
         nnr=nn
         if (nshift.ne.0) then
            nn=nn-nshift
         end if
         
         allocate(x(nn))
         
         if (sincen) then
            x(1)=aa1(1,1)
            do i=1,ivarmd
               aa(i,1)=aa1(i+1,1)
            end do
            do n=2,nnr
               n1=n+nshift
               x(n)=aa1(1,n1)
               do i=1,ivarmd
                  aa(i,n)=aa1(i+1,n1)
               end do
            end do
         else
            do n=1,nnr
               if (n.eq.1) then
                  n1=1
               else
                  n1=n+nshift
               end if
               x(n)=aa1(1,n1)
               do i=1,ivarmd
                  aa(i,n)=aa1(i+1,n1)
               end do
            end do
         end if
         
         deallocate(aa1)

         ! set g/gtilde (=1 in models without turbulent pressure)

         if (iturpr.eq.1) then
            do n=1,nn
               if (x(n).lt.0.999) then
                  ggt=1
               else
                  ggt=1./(x(n)*x(n)*x(n)*aa(1,n))
               end if
               aa(10,n)=ggt
            end do
         else if (iggt.eq.1) then
            do n=1,nn
               aa(10,n)=aa(6,n)
            end do
         else
            do n=1,nn
               aa(10,n)=1
            end do
         end if
         
         x_arg => x
         aa_arg => aa
         nn_arg = nn
         data_arg(:) = data(:)

      end subroutine store_amdl
            
      
      subroutine fgong_amdl( &
            cgrav, nn_in, iconst, ivar, ivers, glob, var, data, aa, nn, ierr)
         ! derived from fgong-amdl.d.f
         real(dp), intent(in) :: cgrav
         integer, intent(in) :: nn_in, iconst, ivar, ivers
         real(dp), intent(inout) :: glob(:) ! (iconst)
         real(dp), intent(inout) :: var(:,:) ! (ivar,nn_in)
         real(dp), intent(inout) :: data(8)
         real(dp), pointer :: aa(:,:) ! (iaa_arg,nn)   will be allocated
         integer, intent(out) :: nn, ierr
         
         integer, parameter :: ireset(16) = &
            (/3,4,5,6,8,9,10,11,12,13,14,16,17,18,19,20/)
         integer :: nn1, i, n, ir
         real(dp) :: d2amax, var1(ivar,nn_in+100), q(nn_in+100), x(nn_in+100)
         real(dp), parameter :: pi4 = 4d0*3.14159265358979323846d0
      
         ierr = 0
         nn = nn_in
      
         if (var(1,1).gt.var(1,nn)) then 
            nn1=nn+1
            do i=1,ivar
               do n=1,nn
                  var1(i,n)=var(i,nn1-n)
               end do
               do n=1,nn
                  var(i,n)=var1(i,n)
               end do
            end do
         end if
         
         if (var(1,1).gt.1.d6) then 
            do i=1,ivar
               do n=1,nn
                  var1(i,n+1)=var(i,n)
               end do
            end do
         
            do i=1,ivar
               var1(i,1)=0
            end do
         
            do ir=1,16
               i=ireset(ir)
               var1(i,1)=var1(i,2)
            end do
         
            nn=nn+1 
            do i=1,ivar
               do n=1,nn
                  var(i,n)=var1(i,n)
               end do
            end do
         end if
         
         do n=1,nn
            q(n)=exp(var(2,n))
            x(n)=var(1,n)/glob(2)
         end do
         
         x(1)=0
         q(1)=0
         
         allocate(aa(iaa_arg,nn))
         
         do n=2,nn
            aa(1,n)=x(n)
            aa(2,n)=q(n)/x(n)**3
            aa(3,n)=cgrav*glob(1)*q(n)*var(5,n)/(var(10,n)*var(4,n)*var(1,n))
            aa(4,n)=var(10,n)
            aa(5,n)=var(15,n)
            aa(6,n)=pi4*var(5,n)*var(1,n)**3/(glob(1)*q(n))
         end do
         
         aa(1,1)=0
         aa(2,1)=pi4/3.d0*var(5,1)*glob(2)**3/glob(1)
         aa(3,1)=0
         aa(4,1)=var(10,1)
         aa(5,1)=0
         aa(6,1)=3.d0
         if (aa(5,nn).le.10) then 
            nn=nn-1 
            !write(6,*) 'Chop off outermost point' 
         end if
         data(1)=glob(1)
         data(2)=glob(2)
         data(3)=var(4,1)
         data(4)=var(5,1)
         if (glob(11).lt.0.and.glob(11).gt.-10000) then 
            data(5)=-glob(11)/var(10,1)
            data(6)=-glob(12) 
         else 
            data(5)=pi4/3.d0*cgrav*(var(5,1)*glob(2))**2/(var(4,1)*var(10,1))
            d2amax=0.d0
            do n=2,nn
               d2amax=max(d2amax,aa(5,n)/x(n)**2)
               if (x(n).ge.0.05d0) exit
            end do
            data(6)=d2amax+data(5)
            !write(6,140) data(5), data(6)
         end if
         data(7)=-1.d0
         data(8)=0.d0
      
      end subroutine fgong_amdl
      

      subroutine read_fgong_file(fin, nn, iconst, ivar, ivers, glob, var, ierr)
         character (len=*), intent(in) :: fin
         integer, intent(out) :: nn, iconst, ivar, ivers
         real(dp), pointer :: glob(:) ! (iconst)   will be allocated
         real(dp), pointer :: var(:,:) ! (ivar,nn)   will be allocated
         integer, intent(out) :: ierr
      
         real(dp), pointer :: var1(:,:) ! (ivar,nn)
         integer :: ios, iounit, i, n, ir, nn1
         character(80) :: head
      
  120 format(4i10)
  130 format(5e16.9)

         ierr = 0
         iounit = 55
         if (ierr /= 0) then
            write(*,*) 'failed in read_fgong_file'
            return
         end if
      
         ios = 0
         open(iounit,file=trim(fin),status='old', iostat=ios)
         if (ios /= 0) then
            write(*,*) 'failed to open ' // trim(fin)
            return
         end if

         do i=1,4
            read(iounit,'(a)', iostat=ios) head
            if (ios /= 0) then
               write(*,*) 'failed to read header line ', i
               return
            end if
         end do

         read(iounit,120, iostat=ios) nn, iconst, ivar, ivers
         if (ios /= 0) then
            write(*,*) 'failed to read dimensions'
            return
         end if
      
         allocate(glob(iconst), var(ivar,nn+10))
      
         read(iounit,130, iostat=ios) (glob(i),i=1,iconst)
         if (ios /= 0) then
            write(*,*) 'failed to read globals'
            return
         end if      

         do n=1,nn
            read(iounit,130, iostat=ios) (var(i,n),i=1,ivar)
            if (ios /= 0) exit
         end do
         close(iounit)
      
         if (ios /= 0) then
            write(*,*) 'failed to read vars'
            return
         end if

      end subroutine read_fgong_file
      
      
      ! for testing
      subroutine dump(filename_for_dump,nn,glob,var,ierr)
         character (len=*), intent(in) :: filename_for_dump
         integer, intent(in) :: nn
         real(dp), pointer :: glob(:) ! (iconst)
         real(dp), pointer :: var(:,:) ! (ivar,nn)
         integer, intent(out) :: ierr
      
         real(dp), parameter :: Msun = 1.9892d33, Rsun = 6.9598d10, Lsun = 3.8418d33
         integer :: iounit, k, offset
      
         ierr = 0
         if (len_trim(filename_for_dump) == 0) return

         iounit = 55
         if (ierr /= 0) then
            write(*,*) 'failed in alloc_iounit for dump fgong'
            return
         end if
      
         open(iounit, file=trim(filename_for_dump),  iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'dump fgong failed to open ' // trim(filename_for_dump)
            return
         end if

         write(*,*) 'dump fgong data to ' // trim(filename_for_dump)
      
         if (VAR(1,1) <= 1) then ! skip tny r
            offset = 1
         else
            offset = 0
         end if
      
         write(iounit,'(99a24)') &
            'num_zones', 'star_mass', 'star_radius', 'star_L', 'initial_z',  &
            'mlt_alpha', 'star_age', 'star_Teff'
         write(iounit,fmt='(i24,99e24.12)') &
            nn-offset, GLOB(1)/Msun, GLOB(2)/Rsun, GLOB(3)/Lsun, GLOB(4),  &
            GLOB(6), GLOB(13), GLOB(14)      
      
         write(iounit,'(a5,99a24)')  &
            'i', 'r', 'm', 'temperature', 'pressure', 'density', &
            'xh1', 'luminosity', 'opacity', 'eps', 'gamma1', &
            'grada', 'chiT_div_chiRho', 'cp', 'free_e', 'brunt_A', &
            'dxdt_nuc_h1', 'z', 'dr_to_surf', 'eps_grav', 'xhe3', &
            'xc12', 'xc13', 'xn14', 'xo16', 'xh2', 'xhe4', 'xli7', &
            'xbe7', 'xn15', 'xo17', 'xo18', 'xne20'                               
      
         do k=1+offset,nn
            write(iounit,fmt='(i5,99e24.12)') k-offset, &
               VAR(1,k), &
               exp(VAR(2,k))*GLOB(1), &
               VAR(3,k), &
               VAR(4,k), &
               VAR(5,k), &
               VAR(6,k), &
               VAR(7,k), &
               VAR(8,k), &
               VAR(9,k), &
               VAR(10,k), &
               VAR(11,k), &
               VAR(12,k), &
               VAR(13,k), &
               VAR(14,k), &
               VAR(15,k), &
               VAR(16,k), &
               VAR(17,k), &
               VAR(18,k), &
               VAR(19,k), &
               VAR(21,k), &
               VAR(22,k), &
               VAR(23,k), &
               VAR(24,k), &
               VAR(25,k), &
               VAR(29,k), &
               VAR(30,k), &
               VAR(31,k), &
               VAR(32,k), &
               VAR(33,k), &
               VAR(34,k), &
               VAR(35,k), &
               VAR(36,k)                               
         end do
         close(iounit)
      
      end subroutine dump


      subroutine show_adipls_results
         integer :: k
         include 'formats'
         do k = 1, num_results
            write(*,4) 'ADIPLS', k, el(k), order(k), cyclic_freq(k), inertia(k)
         end do
         write(*,*)
      end subroutine show_adipls_results


      end module adipls_support
