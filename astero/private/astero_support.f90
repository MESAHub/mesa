! ***********************************************************************
!
!   Copyright (C) 2013  The MESA Team
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
 

      module astero_support
      
      use astero_def
      use star_lib
      use star_def
      use const_def
      use math_lib
      use utils_lib
      use auto_diff
      
      implicit none


      
      contains
      
      subroutine check_search_controls(ierr)
         integer, intent(out) :: ierr
         integer :: i, l
         include 'formats'
         ierr = 0

         do l = 0, 3
            do i = 2, nl(l)
               if (freq_target(l,i) <= freq_target(l,i-1)) then
                  write(*,4) 'freq_target values out of order', l, i-1, i, freq_target(l,i-1), freq_target(l,i)
                  ierr = -1
               end if
            end do
         end do

         if (ierr /= 0) &
            write(*,1) 'please put frequency values in ascending order'      
      end subroutine check_search_controls
      

      subroutine get_one_el_info( &
            s, l, nu1, nu2, iscan, i1, i2, store_model, code, ierr)
         use num_lib, only: qsort
         use adipls_support
         use gyre_support
         type (star_info), pointer :: s
         integer, intent(in) :: l, iscan, i1, i2
         real(dp), intent(in) :: nu1, nu2
         logical, intent(in) :: store_model
         character (len=*), intent(in) :: code
         integer, intent(out) :: ierr
                  
         real(dp) :: nu_obs, dist_j, nu, dist, min_dist, min_freq, &
            R, G, M, sig_fac, b, sum_1, sum_2, sum_3, empty(0)
         integer :: min_dist_j, min_order, n, cnt, int_empty(0)
         integer :: nsel, itrsig, nsig
         real(dp) :: els1, dels, sig1, sig2, dfsig
         integer :: num_l0_terms, k, i, j
         integer, pointer :: index(:) 

         include 'formats'
         
         ierr = 0
         
         if (code == 'gyre') then
         
            if (.not. gyre_is_enabled) then
               ierr = -1
               write(*,'(A)')
               write(*,'(a)') 'gyre is not currently enabled in your configuration of mesa.'
               write(*,'(a)') 'check that your utils/makefile_header has USE_GYRE = YES'
               write(*,'(A)')
               return
            end if

            b = correction_b

            num_results = 0
            call do_gyre_get_modes(s, l, store_model, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in do_gyre_get_modes'
               call mesa_error(__FILE__,__LINE__,'get_one_el_info')
            end if
         
         else if (code == 'adipls') then 

            if (.not. adipls_is_enabled) then
               ierr = -1
               write(*,'(A)')
               write(*,'(a)') 'adipls is not currently enabled in your configuration of mesa.'
               write(*,'(a)') 'check that your utils/makefile_header has USE_ADIPLS = YES'
               write(*,'(A)')
               return
            end if

            R = Rsun*s% photosphere_r
            G = standard_cgrav
            M = s% m_grav(1)
            sig_fac = (2*pi)*(2*pi)*R*R*R/(G*M)
            b = correction_b
         
            ! set controls for adipls
            nsel = 0
            dels = 1
            els1 = dble(l)
            itrsig = 1
            sig1 = sig_fac*(nu1*1d-6)*(nu1*1d-6)
            sig2 = sig_fac*(nu2*1d-6)*(nu2*1d-6)
            dfsig = sig_fac*delta_nu_model*delta_nu_model
            nsig = 2
         
            call set_adipls_controls( &
               l, nsel, els1, dels, itrsig, iscan, sig1, sig2, dfsig, nsig, &
               adipls_irotkr, adipls_nprtkr, adipls_igm1kr, &
               adipls_npgmkr)
      
            num_results = 0
            call run_adipls( &
                 s, .false., store_model, &
                 add_center_point, keep_surface_point, add_atmosphere, &
                 do_redistribute_mesh, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in run_adipls'
               call mesa_error(__FILE__,__LINE__,'get_one_el_info')
            end if
            
         else
         
            write(*,'(a)') 'invalid oscillation_code: ' // trim(oscillation_code)
            ierr = -1
            return
         
         end if

         ! sort results by increasing frequency
         allocate(index(num_results), stat=ierr)
         if (ierr /= 0) then
            call mesa_error(__FILE__,__LINE__,'failed in allocate before calling qsort')
         end if
         call qsort(index, num_results, cyclic_freq)

         if (l == 0) then
            call set_to_closest(freq_target(0,:), &
               model_freq(0,:), empty, empty, &
               model_inertia(0,:), empty, empty, &
               model_order(0,:), int_empty, int_empty, ierr)

            model_freq_alt_up(0,:) = model_freq(0,:)
            model_inertia_alt_up(0,:) = model_inertia(0,:)
            model_order_alt_up(0,:) = model_order(0,:)

            model_freq_alt_down(0,:) = model_freq(0,:)
            model_inertia_alt_down(0,:) = model_inertia(0,:)
            model_order_alt_down(0,:) = model_order(0,:)
         else if (0 < l .and. l <= 3) then
            call set_to_closest(freq_target(l,:), &
               model_freq(l,:), model_freq_alt_up(l,:), model_freq_alt_down(l,:), &
               model_inertia(l,:), model_inertia_alt_up(l,:), model_inertia_alt_down(l,:), &
               model_order(l,:), model_order_alt_up(l,:), model_order_alt_down(l,:), ierr)
         else
            call mesa_error(__FILE__,__LINE__,'bad value for l in get_one_el_info')
         end if
         if (ierr /= 0) then
            !write(*,2) 'failed to find frequency for matching for l =', l
            write(*,2) 'failed to match frequencies for l =', l
            return
         end if
         
         if (l == 0 .and. correction_factor > 0 .and. nl(0) > 0 .and. &
               delta_nu > 0 .and. nu_max > 0 .and. avg_nu_obs > 0) then 
            ! calculate surface correction info
            
            cnt = 0
            sum_1 = 0
            do i=1,nl(0)
               if (freq_target(0,i) < 0) cycle
               cnt = cnt + 1
               sum_1 = sum_1 + model_freq(0,i)
            end do
            if (cnt == 0) return
            avg_nu_model = sum_1/cnt
            
            sum_1 = 0
            sum_2 = 0
            sum_3 = 0
            do i=1,nl(0)
               if (freq_target(0,i) < 0) cycle
               sum_1 = sum_1 + &
                  (model_freq(0,i) - avg_nu_model)*(l0_n_obs(i) - avg_radial_n)
               sum_2 = sum_2 + pow2(l0_n_obs(i) - avg_radial_n)
               sum_3 = sum_3 + pow(freq_target(0,i)/nu_max,b)
            end do
            if (sum_2 == 0 .or. sum_3 == 0) return
            delta_nu_model = sum_1/sum_2
            correction_r = & ! K08 eqn 6
               (b-1)/(b*avg_nu_model/avg_nu_obs - delta_nu_model/delta_nu)
            if (correction_r <= 0) return
            correction_a = & ! K08 eqn 10
               min(0d0, avg_nu_obs - correction_r*avg_nu_model)*nl(0)/sum_3
            a_div_r = correction_a/correction_r
            
         end if
            
         deallocate(index)
         
         
         contains
         

         subroutine set_to_closest( &
               l_obs, &
               l_freq, l_freq_up, l_freq_down, &
               l_inertia, l_inertia_up, l_inertia_down, &
               l_order, l_order_up, l_order_down, &
               ierr)
            real(dp), intent(in) :: l_obs(:)
            real(dp), intent(inout), dimension(:) :: &
               l_freq, l_freq_up, l_freq_down, &
               l_inertia, l_inertia_up, l_inertia_down
            integer, intent(out), dimension(:) :: &
               l_order, l_order_up, l_order_down
            integer, intent(out) :: ierr
            integer :: i, j, jprev, j_alt
            include 'formats'
            jprev = 0
            ierr = 0
            do i = i1, i2
               l_freq(i) = 0
               l_inertia(i) = 0
               l_order(i) = 0
               if (l > 0) then
                  l_freq_up(i) = 0
                  l_inertia_up(i) = 0
                  l_order_up(i) = 0
                  l_freq_down(i) = 0
                  l_inertia_down(i) = 0
                  l_order_down(i) = 0
               end if
               j = find_closest(l_obs(i),jprev)
               if (j <= 0) then
                  ierr = -1
               else
                  l_freq(i) = cyclic_freq(j)
                  l_inertia(i) = inertia(j)
                  l_order(i) = order(j)
                  if (l > 0) then
                     j_alt = find_next_down(j)
                     if (j_alt > 0) then
                        l_freq_down(i) = cyclic_freq(j_alt)
                        l_inertia_down(i) = inertia(j_alt)
                        l_order_down(i) = order(j_alt)
                     end if
                     j_alt = find_next_up(j)
                     if (j_alt > 0) then
                        l_freq_up(i) = cyclic_freq(j_alt)
                        l_inertia_up(i) = inertia(j_alt)
                        l_order_up(i) = order(j_alt)
                     end if
                  end if
                  jprev = j
               end if
            end do            
         end subroutine set_to_closest
            
            
         integer function find_closest(nu,jprev) ! find closest model frequency 
            real(dp), intent(in) :: nu
            integer, intent(in) :: jprev
            min_dist = 1d99; min_dist_j = -1
            do j = jprev+1, num_results
               if (el(j) /= l) cycle
               dist = abs(cyclic_freq(j) - nu)
               if (min_dist_j <= 0 .or. dist < min_dist) then
                  min_dist = dist; min_dist_j = j
               end if
               if (cyclic_freq(j) > nu) exit
            end do
            find_closest = min_dist_j
         end function find_closest
            
            
         integer function find_next_down(j) result(j_down) ! same l, next lower freq 
            integer, intent(in) :: j
            do j_down = j-1, 1, -1
               if (el(j_down) /= l) cycle
               return
            end do
            j_down = 0
         end function find_next_down
            
            
         integer function find_next_up(j) result(j_up) ! same l, next higher freq 
            integer, intent(in) :: j
            do j_up = j+1, num_results
               if (el(j_up) /= l) cycle
               return
            end do
            j_up = 0
         end function find_next_up
         

      end subroutine get_one_el_info
      
      
      subroutine get_frequency_ratios( &
            init, nl0, l0, nl1, l1, n, l0_first, l1_first, r01, r10)
         logical, intent(in) :: init
         integer, intent(in) :: nl0, nl1
         real(dp), intent(in) :: l0(:), l1(:)
         integer, intent(out) :: n, l0_first, l1_first
         real(dp), intent(out) :: r01(:), r10(:)
         
         integer :: l0_seq_n, l0_last, l1_seq_n, l1_last, i, i0, i1
         real(dp) :: d01, d10, sd01, sd10, dnu, sdnu
         
         logical :: dbg
            
         include 'formats'
         
         dbg = .false.
         call fill_with_NaNs(r01)
         call fill_with_NaNs(r10)
         
         n = 0
         
         if (nl1 <= 0) return
      
         call get_max_sequence(nl0, l0, l0_first, l0_seq_n)
         l0_last = l0_first + l0_seq_n - 1
         if (dbg) write(*,4) 'l0_first l0_last l0_seq_n', l0_first, l0_last, l0_seq_n

         call get_max_sequence(nl1, l1, l1_first, l1_seq_n)
         l1_last = l1_first + l1_seq_n - 1
         if (dbg) write(*,4) 'l1_first l1_last l1_seq_n', l1_first, l1_last, l1_seq_n
         
         do ! trim high end of l0 until < last l1
            if (l0(l0_last) < l1(l1_last)) exit
            l0_last = l0_last - 1
            if (l0_last < l0_first) then ! no overlap
               if (dbg) then
                  write(*,*) 'l0_last < l0_first', l0_last, l0_first
               end if
               return
            end if
         end do
         if (dbg) write(*,2) 'l0_last after trim', l0_last
         
         do ! trim low end of l1 until > first l0
            if (l1(l1_first) > l0(l0_first)) exit
            l1_first = l1_first + 1
            if (l1_first > l1_last) then ! no overlap
               return
            end if
         end do
         if (dbg) write(*,2) 'l1_first after trim', l1_first
         
         do ! trim low end of l0 until only 1 < 1st l1
            if (l0_first == l0_last) exit
            if (l0(l0_first+1) >= l1(l1_first)) exit
            l0_first = l0_first + 1
         end do
         if (dbg) write(*,2) 'l0_first after trim', l0_first
         
         do ! trim high end of l1 until only 1 > last l0
            if (l1_last == l1_first) exit
            if (l1(l1_last-1) <= l0(l0_last)) exit
            l1_last = l1_last - 1
         end do
         if (dbg) write(*,2) 'l1_last after trim', l1_last
         
         l0_seq_n = l0_last - l0_first + 1
         l1_seq_n = l1_last - l1_first + 1
         n = l0_seq_n - 2
         if (dbg) write(*,2) 'n', n
         
         if (l0_seq_n /= l1_seq_n .or. n < 1) then
            return
         end if
         
         do i = 1, n         
            i0 = i + l0_first
            i1 = i + l1_first
            d01 = (l0(i0-1) - 4*l1(i1-1) + 6*l0(i0) - 4*l1(i1) + l0(i0+1))/8d0
            r01(i) = d01/(l1(i1) - l1(i1-1))
            d10 = -(l1(i1-1) - 4*l0(i0) + 6*l1(i1) - 4*l0(i0+1) + l1(i1+1))/8d0
            r10(i) = d10/(l0(i0+1) - l0(i0))
            if (.not. init) cycle
            ! set ratio sigmas
            sd01 = sqrt(pow2(freq_sigma(0,i0-1)) + pow2(4*freq_sigma(1,i1-1)) + &
               pow2(6*freq_sigma(0,i0)) + pow2(4*freq_sigma(1,i1)) + pow2(freq_sigma(0,i0+1)))/8d0
            dnu = l1(i1) - l1(i1-1)
            sdnu = sqrt(pow2(freq_sigma(1,i1)) + pow2(freq_sigma(1,i1-1)))
            sigmas_r01(i) = sqrt(pow2(sd01/dnu) + pow2(sdnu*d01/(dnu*dnu)))
            sd10 = sqrt(pow2(freq_sigma(1,i1-1)) + pow2(4*freq_sigma(0,i0)) + &
               pow2(6*freq_sigma(1,i1)) + pow2(4*freq_sigma(0,i0+1)) + pow2(freq_sigma(1,i1+1)))/8d0
            dnu = l0(i0+1) - l0(i0)
            sdnu = sqrt(pow2(freq_sigma(0,i0+1)) + pow2(freq_sigma(0,i0)))
            sigmas_r10(i) = sqrt(pow2(sd10/dnu) + pow2(sdnu*d10/(dnu*dnu)) )
            if (trace_chi2_seismo_ratios_info) then
               write(*,'(a30,i4,99f16.6)') 'r01 r10 sigmas_r01 sigmas_r10', &
                  i, r01(i), r10(i), sigmas_r01(i), sigmas_r10(i)
            end if
         end do
         
      end subroutine get_frequency_ratios
      
      
      subroutine get_r02_frequency_ratios(init, nl0, l0, nl1, l1, nl2, l2, r02)
         logical, intent(in) :: init
         integer, intent(in) :: nl0, nl1, nl2
         real(dp), intent(in) :: l0(:), l1(:), l2(:)
         real(dp), intent(out) :: r02(:)
         
         integer :: i, i0, i1, i2, jmin, j
         real(dp) :: d02, sd02, dnu, sdnu, df, f0, f2, fmin, fmax, dfmin
         
         logical :: dbg
            
         include 'formats'
         
         dbg = .false.
         call fill_with_NaNs(r02)
         
         if (init) then ! set i2_for_r02
            do i = 1, ratios_n         
               i0 = i + ratios_l0_first
               i1 = i + ratios_l1_first            
               dnu = l1(i1) - l1(i1-1)
               df = 0.25*dnu
               f0 = l0(i0)
               fmin = f0 - df
               fmax = f0 + df
               dfmin = 1d99
               jmin = 0
               do j = 1, nl2
                  f2 = l2(j)
                  !if (i==1) write(*,2) 'f2', j, fmin, f0, f2, fmax, abs(f2 - f0), dfmin
                  if (f2 <= fmax .and. f2 >= fmin .and. &
                        abs(f2 - f0) < dfmin) then
                     dfmin = abs(f2 - f0)
                     jmin = j
                  end if
                  if (f2 > f0) exit
               end do
               !if (.true.) write(*,3) 'i2_for_r02', i, jmin, fmin, f0, fmax, dnu
               i2_for_r02(i) = jmin
               sigmas_r02(i) = 0d0
            end do
         end if
         !write(*,2) 'ratios_n', ratios_n
         !stop
         
         do i = 1, ratios_n
            if ((.not. init) .and. sigmas_r02(i) == 0d0) cycle
            i2 = i2_for_r02(i)
            if (i2 == 0) cycle
            i0 = i + ratios_l0_first
            i1 = i + ratios_l1_first            
            d02 = l0(i0) - l2(i2)
            dnu = l1(i1) - l1(i1-1)
            r02(i) = d02/dnu
            if (.not. init) cycle
            ! set ratio sigmas
            sd02 = sqrt(pow2(freq_sigma(0,i0)) + pow2(freq_sigma(2,i2)))
            sdnu = sqrt(pow2(freq_sigma(1,i1)) + pow2(freq_sigma(1,i1-1)))
            sigmas_r02(i) = sqrt(pow2(sd02/dnu) + pow2(sdnu*d02/(dnu*dnu)))
            if (trace_chi2_seismo_ratios_info) then
               write(*,'(a30,i4,99f16.6)') 'r02 sigmas_r02', &
                  i, r02(i), sigmas_r02(i)
            end if
         end do
         
      end subroutine get_r02_frequency_ratios
      
      
      real(dp) function interpolate_ratio_r010( &
            freq, first, model_freqs, model_ratios) result(ratio)
         real(dp), intent(in) :: freq
         integer, intent(in) :: first
         real(dp), intent(in), dimension(:) :: model_freqs, model_ratios
         integer :: i, j
         real(dp) :: alfa, beta
         ratio = 0
         if (ratios_n == 0) return
         j = 1 + first
         if (freq <= model_freqs(j)) then
            ratio = model_ratios(j)
            return
         end if
         j = ratios_n + first
         if (freq >= model_freqs(j)) then
            ratio = model_ratios(j)
            return
         end if
         do i=2,ratios_n
            j = i+first
            if (freq < model_freqs(j)) then
               alfa = (freq - model_freqs(j-1))/(model_freqs(j) - model_freqs(j-1))
               beta = 1d0 - alfa
               ratio = alfa*model_ratios(j) + beta*model_ratios(j-1)
               return
            end if
         end do
      end function interpolate_ratio_r010

      
      real(dp) function interpolate_ratio_r02( &
            freq, model_freqs, model_ratios) result(ratio)
         real(dp), intent(in) :: freq
         real(dp), intent(in), dimension(:) :: model_freqs, model_ratios
         integer :: i, i_lo, i_hi
         real(dp) :: alfa, beta
         ratio = 0
         i_lo = 0
         do i=1,nl(0)
            if (sigmas_r02(i) == 0) cycle
            i_lo = i
            exit
         end do
         if (i_lo == 0) return
         if (freq <= model_freqs(i_lo)) then
            ratio = model_ratios(i_lo)
            return
         end if
         i_hi = i_lo
         do i=i_lo+1,nl(0)
            if (sigmas_r02(i) == 0) cycle
            i_hi = i
            if (freq > model_freqs(i)) cycle
            alfa = (freq - model_freqs(i_lo))/ &
               (model_freqs(i_hi) - model_freqs(i_lo))
            beta = 1d0 - alfa
            ratio = alfa*model_ratios(i_hi) + beta*model_ratios(i_lo)
            return
         end do
         ratio = model_ratios(i_hi)
      end function interpolate_ratio_r02


      subroutine get_max_sequence(nl, l_obs, max_seq_i, max_seq_n)
         integer, intent(in) :: nl
         real(dp), intent(in) :: l_obs(:)
         integer, intent(out) :: max_seq_i, max_seq_n
         
         integer :: i, j, seq_i, seq_n
         
         max_seq_i = 0
         max_seq_n = 0
         seq_i = 0
         seq_n = 0
      
         do 
            i = seq_i + seq_n + 1 ! start of next sequence
            if (i >= nl) exit
            seq_i = i
            seq_n = 1
            do j = seq_i, nl-1 ! j is in series; try to add j+1
               if (l_obs(j+1) - l_obs(j) > 1.5*delta_nu) then ! end of series
                  if (seq_n > max_seq_n) then
                     max_seq_i = seq_i
                     max_seq_n = seq_n
                  end if
                  exit
               end if
               seq_n = seq_n + 1
            end do
         end do
      
         if (seq_n > max_seq_n) then
            max_seq_i = seq_i
            max_seq_n = seq_n
         end if
         
      end subroutine get_max_sequence
               

      subroutine init_obs_data(ierr)
         integer, intent(out) :: ierr
         
         integer :: i, cnt, norders
         integer, dimension(max_nl) :: orders
         real(dp) :: sum_1, sum_2, sum_3, range, nmax
         real(dp) :: x, y, isig2, sum_xy, sum_x, sum_y, sum_x2, sum_isig2, d
         
         logical, parameter :: dbg = .false.
         
         include 'formats'
         
         ierr = 0
         
         !call test_get_frequency_ratios
         
         if (nl(0) <= 0) return
         
         sigmas_r02 = 0d0
         ratios_r02 = 0d0
         
         if (chi2_seismo_r_010_fraction > 0 .or. &
             chi2_seismo_r_02_fraction > 0) then
            call get_frequency_ratios( &
               .true., nl(0), freq_target(0,:), nl(1), freq_target(1,:), &
               ratios_n, ratios_l0_first, ratios_l1_first, &
               ratios_r01, ratios_r10)
         end if
         if (chi2_seismo_r_02_fraction > 0) then
            call get_r02_frequency_ratios(.true., &
               nl(0), freq_target(0,:), &
               nl(1), freq_target(1,:), &
               nl(2), freq_target(2,:), ratios_r02)
         end if
         
         if (delta_nu <= 0 .and. nl(0) > 1 .and. l0_n_obs(1) > 0) then
            sum_xy = 0
            sum_x = 0
            sum_y = 0
            sum_x2 = 0
            sum_isig2 = 0
            do i=1,nl(0)
               isig2 = 1d0/pow2(freq_sigma(0,i))
               x = dble(l0_n_obs(i))
               y = freq_target(0,i)
               sum_xy = sum_xy + x*y*isig2
               sum_x = sum_x + x*isig2
               sum_y = sum_y + y*isig2
               sum_x2 = sum_x2 + x*x*isig2
               sum_isig2 = sum_isig2 + isig2
            end do
            d = sum_isig2*sum_x2 - sum_x*sum_x
            delta_nu = (sum_isig2*sum_xy - sum_x*sum_y)/d
            if (delta_nu_sigma <= 0) delta_nu_sigma = sqrt(sum_isig2/d)            
         end if
         
         ! if (correction_factor <= 0) return
         if (correction_scheme /= 'kjeldsen') return
         
         if (l0_n_obs(1) <= 0) then
            if (delta_nu <= 0) then
               write(*,*) 'must supply value for delta_nu'
               ierr = -1
               return
            end if
            ! set l0_n_obs(i) to order of freq_target(0,i)
            range = freq_target(0,nl(0)) - freq_target(0,1)
            norders = int(range/delta_nu + 0.5d0) + 1
            nmax = (nu_max/delta_nu)*(delta_nu_sun/nu_max_sun)*22.6 - 1.6         
            l0_n_obs(1) = int(nmax - (norders-1)/2)
            if (dbg) write(*,3) 'l0_n_obs(i)', 1, l0_n_obs(1), freq_target(0,1)
            do i=2,norders
               l0_n_obs(i) = l0_n_obs(1) + &
                  int((freq_target(0,i) - freq_target(0,1))/delta_nu + 0.5)
               if (dbg) write(*,3) 'l0_n_obs(i)', i, l0_n_obs(i), freq_target(0,i)
            end do
            if (dbg) then
               write(*,1) 'range', range
               write(*,2) 'norders', norders
               write(*,1) 'nmax', nmax
               write(*,2) '(norders+1)/2', (norders+1)/2
               write(*,2) 'l0_n_obs(1)', l0_n_obs(1)
               write(*,'(A)')
               !stop
            end if
         end if 
         
         
         cnt = 0
         sum_1 = 0
         sum_2 = 0
         do i=1,nl(0)
            if (freq_target(0,i) < 0) cycle
            cnt = cnt + 1
            sum_1 = sum_1 + freq_target(0,i)
            sum_2 = sum_2 + l0_n_obs(i)
         end do
         avg_nu_obs = sum_1/cnt
         avg_radial_n = sum_2/cnt
         
         if (dbg) then
            write(*,1) 'avg_nu_obs', avg_nu_obs
            write(*,1) 'avg_radial_n', avg_radial_n
            write(*,2) 'cnt', cnt
            write(*,1) 'sum_1', sum_1
            write(*,1) 'sum_2', sum_2
            write(*,'(A)')
            call mesa_error(__FILE__,__LINE__,'init_obs_data')
         end if
            
      end subroutine init_obs_data
      
      
      real(dp) function interpolate_l0_inertia(freq) result(inertia)
         real(dp), intent(in) :: freq
         integer :: i
         real(dp) :: alfa, beta
         inertia = 0
         if (nl(0) == 0) return
         if (freq <= model_freq(0,1)) then
            inertia = model_inertia(0,1)
            return
         end if
         if (freq >= model_freq(0,nl(0))) then
            inertia = model_inertia(0,nl(0))
            return
         end if
         do i=2,nl(0)
            if (freq < model_freq(0,i)) then
               alfa = (freq - model_freq(0,i-1))/(model_freq(0,i) - model_freq(0,i-1))
               beta = 1d0 - alfa
               inertia = alfa*model_inertia(0,i) + beta*model_inertia(0,i-1)
               return
            end if
         end do
      end function interpolate_l0_inertia
      
      
      subroutine get_kjeldsen_radial_freq_corr( &
            a_div_r, b, nu_max, correction_factor, check_obs, &
            nl0, l0_obs, l0_freq, l0_freq_corr, l0_inertia)
         real(dp), intent(in) :: a_div_r, b, nu_max, correction_factor
         logical, intent(in) :: check_obs ! if false, then l0_obs is not used
         integer, intent(in) :: nl0
         real(dp), intent(in), dimension(:) :: &
            l0_obs, l0_freq, l0_inertia
         real(dp), intent(inout) :: l0_freq_corr(:)
         integer :: i
         real(dp) :: Qnl
         do i = 1, nl0
            if (check_obs) then
               if (l0_obs(i) < 0) cycle
            end if
            Qnl = 1
            l0_freq_corr(i) = l0_freq(i)
            if (b > 0 .and. correction_factor > 0 .and. l0_freq_corr(i) > 0) &
               l0_freq_corr(i) = l0_freq_corr(i) + &
                  correction_factor*(a_div_r/Qnl)*pow(l0_freq(i)/nu_max,b)
         end do
      end subroutine get_kjeldsen_radial_freq_corr
      
      
      subroutine get_kjeldsen_nonradial_freq_corr( &
            a_div_r, b, nu_max, correction_factor, check_obs, &
            nl1, l1_obs, l1_freq, l1_freq_corr, l1_inertia, l0_inertia)
         real(dp), intent(in) :: a_div_r, b, nu_max, correction_factor
         logical, intent(in) :: check_obs ! if false, then l1_obs is not used
         integer, intent(in) :: nl1
         real(dp), intent(in), dimension(:) :: &
            l1_obs, l1_freq, l1_inertia, l0_inertia
         real(dp), intent(inout) :: l1_freq_corr(:)
         integer :: i
         real(dp) :: Qnl, interp_l0_inertia
         include 'formats'
         do i = 1, nl1
            if (check_obs) then
               if (l1_obs(i) < 0) cycle
            end if
            l1_freq_corr(i) = l1_freq(i)
            if (b > 0 .and. correction_factor > 0 .and. l1_freq_corr(i) > 0 .and. nl(0) > 0) then
               interp_l0_inertia = interpolate_l0_inertia(l1_freq(i))
               !write(*,2) 'l1_freq_corr: l0_inertia interp prev', i, interp_l0_inertia, &
               !   (l0_inertia(min(nl(0),i)) + l0_inertia(min(nl(0),i+1)))/2
               Qnl = l1_inertia(i)/interp_l0_inertia
               l1_freq_corr(i) = l1_freq_corr(i) + &
                  correction_factor*(a_div_r/Qnl)*pow(l1_freq(i)/nu_max,b)
            end if
         end do
      end subroutine get_kjeldsen_nonradial_freq_corr

      
      subroutine get_kjeldsen_nonradial_freq_corr_alt_up
         !call mesa_error(__FILE__,__LINE__,'get_kjeldsen_nonradial_freq_corr_alt_up')
      end subroutine get_kjeldsen_nonradial_freq_corr_alt_up

      
      subroutine get_kjeldsen_nonradial_freq_corr_alt_down 
         !call mesa_error(__FILE__,__LINE__,'get_kjeldsen_nonradial_freq_corr_alt_down')
      end subroutine get_kjeldsen_nonradial_freq_corr_alt_down

      
      subroutine get_kjeldsen_freq_corr
         integer :: l

         call get_kjeldsen_radial_freq_corr( &
            a_div_r, correction_b, nu_max, correction_factor, .true., &
            nl(0), freq_target(0,:), model_freq(0,:), model_freq_corr(0,:), model_inertia(0,:))

         do l = 1, 3
            call get_kjeldsen_nonradial_freq_corr( &
               a_div_r, correction_b, nu_max, correction_factor, .true., &
               nl(l), freq_target(l,:), model_freq(l,:), model_freq_corr(l,:), model_inertia(l,:), model_inertia(0,:))
         end do

      end subroutine get_kjeldsen_freq_corr

      
      subroutine get_kjeldsen_freq_corr_alt_up
         integer :: l

         do l = 1, 3
            call get_kjeldsen_nonradial_freq_corr( &
               a_div_r, correction_b, nu_max, correction_factor, .true., &
               nl(l), freq_target(l,:), model_freq_alt_up(l,:), model_freq_corr_alt_up(l,:), &
               model_inertia_alt_up(l,:), model_inertia(0,:))
         end do

      end subroutine get_kjeldsen_freq_corr_alt_up

      
      subroutine get_kjeldsen_freq_corr_alt_down 
         integer :: l

         do l = 1, 3
            call get_kjeldsen_nonradial_freq_corr( &
               a_div_r, correction_b, nu_max, correction_factor, .true., &
               nl(l), freq_target(l,:), model_freq_alt_down(l,:), model_freq_corr_alt_down(l,:), &
               model_inertia_alt_down(l,:), model_inertia(0,:))
         end do

      end subroutine get_kjeldsen_freq_corr_alt_down
      
      
      subroutine get_no_freq_corr
         integer :: i, l

         do l = 0, 3
            do i = 1, nl(l)
               model_freq_corr(l,i) = model_freq(l,i)
            end do
         end do
      end subroutine get_no_freq_corr
      
      
      subroutine get_no_freq_corr_alt_up
         integer :: i, l

         do l = 0, 3
            do i = 1, nl(l)
               model_freq_corr_alt_up(l,i) = model_freq_alt_up(l,i)
            end do
         end do
      end subroutine get_no_freq_corr_alt_up
      
      
      subroutine get_no_freq_corr_alt_down
         integer :: i, l

         do l = 0, 3
            do i = 1, nl(l)
               model_freq_corr_alt_down(l,i) = model_freq_alt_down(l,i)
            end do
         end do
      end subroutine get_no_freq_corr_alt_down


      subroutine get_cubic_all_freq_corr(a3, radial_only, &
            nl, obs, sigma, freq, freq_corr, inertia)
           
        integer, intent(in) :: nl(0:)
        real(dp), intent(in), dimension(0:,:) :: &
             obs, sigma, freq, inertia
        real(dp), intent(inout), dimension(0:,:) :: freq_corr
        real(dp), intent(out) :: a3

        logical :: radial_only
        real(dp) :: X, y, XtX, Xty
        integer :: i, l

        XtX = 0d0
        Xty = 0d0

        do i = 1, nl(0)
           X = pow3(freq(0,i))/inertia(0,i)/sigma(0,i)
           y = (obs(0,i) - freq(0,i))/sigma(0,i)

           XtX = XtX + X*X
           Xty = Xty + X*y
        end do

        if (.not. radial_only) then
           do l = 1, 3
              do i = 1, nl(l)
                 X = pow3(freq(l,i))/inertia(l,i)/sigma(l,i)
                 y = (obs(l,i) - freq(l,i))/sigma(l,i)

                 XtX = XtX + X*X
                 Xty = Xty + X*y
              end do
           end do
        end if

        a3 = Xty/XtX

        do l = 0, 3
           do i = 1, nl(l)
              freq_corr(l,i) = freq(l,i) + correction_factor*a3*pow3(freq(l,i))/inertia(l,i)
           end do
        end do

      end subroutine get_cubic_all_freq_corr


      subroutine get_cubic_freq_corr(radial_only)
         logical, intent(in) :: radial_only
         call get_cubic_all_freq_corr(a3, radial_only, &
            nl, freq_target, freq_sigma, model_freq, model_freq_corr, model_inertia)
      end subroutine get_cubic_freq_corr
      
      
      subroutine get_cubic_freq_corr_alt_up(radial_only)
         logical, intent(in) :: radial_only
         call get_cubic_all_freq_corr(a3, radial_only, &
            nl, freq_target, freq_sigma, model_freq_alt_up, model_freq_corr, model_inertia_alt_up)
      end subroutine get_cubic_freq_corr_alt_up
      
      
      subroutine get_cubic_freq_corr_alt_down(radial_only)
         logical, intent(in) :: radial_only
         call get_cubic_all_freq_corr(a3, radial_only, &
            nl, freq_target, freq_sigma, model_freq_alt_down, model_freq_corr_alt_down, model_inertia_alt_down)
      end subroutine get_cubic_freq_corr_alt_down


      subroutine get_combined_all_freq_corr(a3, a1, radial_only, &
            nl, obs, sigma, freq, freq_corr, inertia)

        integer, intent(in) :: nl(0:)
        real(dp), intent(in), dimension(0:,:) :: &
             obs, sigma, freq, inertia
        real(dp), intent(inout), dimension(0:,:) :: freq_corr
        real(dp), intent(out) :: a3, a1
        logical :: radial_only
        
        integer :: i, l
        real(dp) :: X(2), XtX(2,2), XtXi(2,2), Xty(2), y
        real(dp) :: detXtX

        XtX = 0d0
        Xty = 0d0

        do i = 1, nl(0)
           X(1) = powm1(freq(0,i))/inertia(0,i)/sigma(0,i)
           X(2) = pow3(freq(0,i))/inertia(0,i)/sigma(0,i)
           y = (obs(0,i) - freq(0,i))/sigma(0,i)

           XtX(1,1) = XtX(1,1) + X(1)*X(1)
           XtX(1,2) = XtX(1,2) + X(1)*X(2)
           XtX(2,2) = XtX(2,2) + X(2)*X(2)
           Xty(1) = Xty(1) + X(1)*y
           Xty(2) = Xty(2) + X(2)*y
        end do

        if (.not. radial_only) then
           do l = 1, 3
              do i = 1, nl(l)
                 X(1) = powm1(freq(l,i))/inertia(l,i)/sigma(l,i)
                 X(2) = pow3(freq(l,i))/inertia(l,i)/sigma(l,i)
                 y = (obs(l,i) - freq(l,i))/sigma(l,i)

                 XtX(1,1) = XtX(1,1) + X(1)*X(1)
                 XtX(1,2) = XtX(1,2) + X(1)*X(2)
                 XtX(2,2) = XtX(2,2) + X(2)*X(2)
                 Xty(1) = Xty(1) + X(1)*y
                 Xty(2) = Xty(2) + X(2)*y
              end do
           end do
        end if

        XtX(2,1) = XtX(1,2)

        XtXi(1,1) = XtX(2,2)
        XtXi(2,2) = XtX(1,1)
        XtXi(1,2) = -XtX(1,2)
        XtXi(2,1) = -XtX(2,1)

        detXtX = XtX(1,1)*XtX(2,2) - XtX(1,2)*XtX(2,1)
        XtXi = XtXi/detXtX

        a1 = XtXi(1,1)*Xty(1) + XtXi(1,2)*Xty(2)
        a3 = XtXi(2,1)*Xty(1) + XtXi(2,2)*Xty(2)

        do l = 0, 3
           do i = 1, nl(l)
              freq_corr(l,i) = freq(l,i) + &
                    correction_factor*(a1*powm1(freq(l,i))+a3*pow3(freq(l,i)))/inertia(l,i)
           end do
        end do

      end subroutine get_combined_all_freq_corr


      subroutine get_combined_freq_corr(radial_only)
         logical, intent(in) :: radial_only
         call get_combined_all_freq_corr(a3, a1, radial_only, &
            nl, freq_target, freq_sigma, model_freq, model_freq_corr, model_inertia)
      end subroutine get_combined_freq_corr
      
      
      subroutine get_combined_freq_corr_alt_up(radial_only)
         logical, intent(in) :: radial_only
         call get_combined_all_freq_corr(a3, a1, radial_only, &
            nl, freq_target, freq_sigma, model_freq_alt_up, model_freq_corr_alt_up, model_inertia_alt_up)
      end subroutine get_combined_freq_corr_alt_up
      
      
      subroutine get_combined_freq_corr_alt_down(radial_only)
         logical, intent(in) :: radial_only
         call get_combined_all_freq_corr(a3, a1, radial_only, &
            nl, freq_target, freq_sigma, model_freq_alt_down, model_freq_corr_alt_down, model_inertia_alt_down)
      end subroutine get_combined_freq_corr_alt_down
      
      
      type(auto_diff_real_2var_order1) function power_law(freq, freq_ref, a, b)
        real(dp), intent(in) :: freq, freq_ref, a, b
        type(auto_diff_real_2var_order1) :: a_ad, b_ad

        a_ad = a
        a_ad%d1val1 = 1.0_dp

        b_ad = b
        b_ad%d1val2 = 1.0_dp
        
        power_law = a_ad*pow(freq/freq_ref, b_ad)
      end function power_law
      
      
      subroutine get_power_law_all_freq_corr(a, b, radial_only, freq_ref, &
            nl, obs, sigma, freq, freq_corr, inertia)

        integer, intent(in) :: nl(0:)
        real(dp), intent(in) :: freq_ref
        real(dp), intent(in), dimension(0:,:) :: &
             obs, sigma, freq, inertia
        real(dp), intent(out), dimension(0:,:) :: freq_corr
        real(dp), intent(out) :: a, b

        logical :: radial_only
        integer :: i, l, iter
        real(dp) :: X(2), XtX(2,2), XtXi(2,2), Xty(2), y
        real(dp) :: detXtX, da, db
        real(dp) :: Q(0:3,max_nl)
        type(auto_diff_real_2var_order1) :: power_law_ad

        ! Power_Law's solar values happen to be the same as MESA's but the
        ! commented expression is there in case it changes
        ! nu_max = s% nu_max ! *(3100d0/3100d0)/sqrt(5777d0/5777d0)

        ! initial guesses are solar-calibrated values from Ball et al. (2016)
        a = -5.25d0
        b = 5.37d0

        do iter=1,1000
           XtX = 0d0
           Xty = 0d0

           do i = 1, nl(0)
              Q(0,i) = 1

              power_law_ad = power_law(freq(0,i), freq_ref, a, b)

              X(1) = -power_law_ad%d1val1/sigma(0,i) ! dpower_law/da
              X(2) = -power_law_ad%d1val2/sigma(0,i) ! dpower_law/db
              y = (obs(0,i) - freq(0,i) - power_law_ad%val)/sigma(0,i)
              
              XtX(1,1) = XtX(1,1) + X(1)*X(1)
              XtX(1,2) = XtX(1,2) + X(1)*X(2)
              XtX(2,2) = XtX(2,2) + X(2)*X(2)
              Xty(1) = Xty(1) + X(1)*y
              Xty(2) = Xty(2) + X(2)*y
           end do

           if (.not. radial_only) then
              do l = 1, 3
                 do i = 1, nl(l)
                    Q(l,i) = inertia(l,i)/interpolate_l0_inertia(freq(l,i))

                    power_law_ad = power_law(freq(l,i), freq_ref, a, b)

                    X(1) = -power_law_ad%d1val1/sigma(l,i)
                    X(2) = -power_law_ad%d1val2/sigma(l,i)
                    y = ((obs(l,i) - freq(l,i))*Q(l,i) - power_law_ad%val)/sigma(l,i)
                 
                    XtX(1,1) = XtX(1,1) + X(1)*X(1)
                    XtX(1,2) = XtX(1,2) + X(1)*X(2)
                    XtX(2,2) = XtX(2,2) + X(2)*X(2)
                    Xty(1) = Xty(1) + X(1)*y
                    Xty(2) = Xty(2) + X(2)*y
                 end do
              end do
           end if

           XtX(2,1) = XtX(1,2)
           
           XtXi(1,1) = XtX(2,2)
           XtXi(2,2) = XtX(1,1)
           XtXi(1,2) = -XtX(1,2)
           XtXi(2,1) = -XtX(2,1)
              
           detXtX = XtX(1,1)*XtX(2,2) - XtX(1,2)*XtX(2,1)
           XtXi = XtXi/detXtX
           
           da = XtXi(1,1)*Xty(1) + XtXi(1,2)*Xty(2)
           db = XtXi(2,1)*Xty(1) + XtXi(2,2)*Xty(2)
           
           if ((da /= da) .or. (db /= db)) exit
           
           a = a - da
           b = b - db

           if ((abs(da)<1d-8) .and. (abs(db)<1d-6)) exit
        end do

        ! if ((da /= da) .or. (db /= db)) then
        !    write(*,*) 'NaN in Power_Law surface correction', iter
        !    write(*,*) freq(0,1:nl(0))
        !    if (nl(1) > 0) write(*,*) freq(1,1:nl(1))
        !    if (nl(2) > 0) write(*,*) freq(2,1:nl(2))
        !    if (nl(3) > 0) write(*,*) freq(3,1:nl(3))
        ! end if

        do l = 0, 3
           do i = 1, nl(l)
              power_law_ad = power_law(freq(l,i), freq_ref, a, b)
              freq_corr(l,i) = freq(l,i) + correction_factor*power_law_ad%val/Q(l,i)
           end do
        end do

      end subroutine get_power_law_all_freq_corr


      subroutine get_power_law_freq_corr(radial_only, freq_ref)
         logical, intent(in) :: radial_only
         real(dp), intent(in) :: freq_ref
         call get_power_law_all_freq_corr(power_law_a, power_law_b, radial_only, freq_ref, &
            nl, freq_target, freq_sigma, model_freq, model_freq_corr, model_inertia)
      end subroutine get_power_law_freq_corr
      
      
      subroutine get_power_law_freq_corr_alt_up(radial_only, freq_ref)
         logical, intent(in) :: radial_only
         real(dp), intent(in) :: freq_ref
         call get_power_law_all_freq_corr(power_law_a, power_law_b, radial_only, freq_ref, &
            nl, freq_target, freq_sigma, model_freq_alt_up, model_freq_corr_alt_up, model_inertia_alt_up)
      end subroutine get_power_law_freq_corr_alt_up
      
      
      subroutine get_power_law_freq_corr_alt_down(radial_only, freq_ref)
         logical, intent(in) :: radial_only
         real(dp), intent(in) :: freq_ref
         call get_power_law_all_freq_corr(power_law_a, power_law_b, radial_only, freq_ref, &
            nl, freq_target, freq_sigma, model_freq_alt_down, model_freq_corr_alt_down, model_inertia_alt_down)
      end subroutine get_power_law_freq_corr_alt_down
      
      
      type(auto_diff_real_2var_order1) function sonoi(freq, freq_ref, a, b)
        real(dp), intent(in) :: freq, freq_ref, a, b
        type(auto_diff_real_2var_order1) :: a_ad, b_ad

        a_ad = a
        a_ad%d1val1 = 1.0_dp

        b_ad = b
        b_ad%d1val2 = 1.0_dp
        
        sonoi = a_ad*freq_ref*(1d0 - 1d0/(1d0+pow(freq/freq_ref, b_ad)))
      end function sonoi


      subroutine get_sonoi_all_freq_corr(a, b, radial_only, freq_ref, &
            nl, obs, sigma, freq, freq_corr, inertia)

        integer, intent(in) :: nl(0:)
        real(dp), intent(in) :: freq_ref
        real(dp), intent(in), dimension(0:,:) :: &
             obs, sigma, freq, inertia
        real(dp), intent(out), dimension(0:,:) :: freq_corr
        real(dp), intent(out) :: a, b

        logical :: radial_only
        integer :: i, l, iter
        real(dp) :: X(2), XtX(2,2), XtXi(2,2), Xty(2), y
        real(dp) :: detXtX, da, db
        real(dp) :: Q(0:3,max_nl)
        type(auto_diff_real_2var_order1) :: sonoi_ad

        ! Sonoi's solar values happen to be the same as MESA's but the
        ! commented expression is there in case it changes
        ! nu_max = s% nu_max ! *(3100d0/3100d0)/sqrt(5777d0/5777d0)

        a = -3.59d-3
        b = 11.26d0

        do iter=1,1000
           XtX = 0d0
           Xty = 0d0

           do i = 1, nl(0)
              Q(0,i) = 1

              sonoi_ad = sonoi(freq(0,i), freq_ref, a, b)

              X(1) = -sonoi_ad%d1val1/sigma(0,i)
              X(2) = -sonoi_ad%d1val2/sigma(0,i)
              y = (obs(0,i) - freq(0,i) - sonoi_ad%val)/sigma(0,i)
              
              XtX(1,1) = XtX(1,1) + X(1)*X(1)
              XtX(1,2) = XtX(1,2) + X(1)*X(2)
              XtX(2,2) = XtX(2,2) + X(2)*X(2)
              Xty(1) = Xty(1) + X(1)*y
              Xty(2) = Xty(2) + X(2)*y
           end do

           if (.not. radial_only) then
              do l = 1, 3
                 do i = 1, nl(l)
                    Q(l,i) = inertia(l,i)/interpolate_l0_inertia(freq(l,i))

                    sonoi_ad = sonoi(freq(l,i), freq_ref, a, b)

                    X(1) = -sonoi_ad%d1val1/sigma(l,i)
                    X(2) = -sonoi_ad%d1val2/sigma(l,i)
                    y = ((obs(l,i) - freq(l,i))*Q(l,i) - sonoi_ad%val)/sigma(l,i)
                 
                    XtX(1,1) = XtX(1,1) + X(1)*X(1)
                    XtX(1,2) = XtX(1,2) + X(1)*X(2)
                    XtX(2,2) = XtX(2,2) + X(2)*X(2)
                    Xty(1) = Xty(1) + X(1)*y
                    Xty(2) = Xty(2) + X(2)*y
                 end do
              end do
           end if

           XtX(2,1) = XtX(1,2)
           
           XtXi(1,1) = XtX(2,2)
           XtXi(2,2) = XtX(1,1)
           XtXi(1,2) = -XtX(1,2)
           XtXi(2,1) = -XtX(2,1)
              
           detXtX = XtX(1,1)*XtX(2,2) - XtX(1,2)*XtX(2,1)
           XtXi = XtXi/detXtX
           
           da = XtXi(1,1)*Xty(1) + XtXi(1,2)*Xty(2)
           db = XtXi(2,1)*Xty(1) + XtXi(2,2)*Xty(2)
           
           if ((da /= da) .or. (db /= db)) exit
           
           a = a - da
           b = b - db

           if ((abs(da)<1d-8) .and. (abs(db)<1d-6)) exit
        end do

        ! if ((da /= da) .or. (db /= db)) then
        !    write(*,*) 'NaN in Sonoi surface correction', iter
        !    write(*,*) freq(0,1:nl(0))
        !    if (nl(1) > 0) write(*,*) freq(1:nl(1))
        !    if (nl(2) > 0) write(*,*) freq(1:nl(2))
        !    if (nl(3) > 0) write(*,*) freq(1:nl(3))
        ! end if

        if (b < 0d0) a = 0d0

        do l = 0, 3
           do i = 1, nl(l)
              sonoi_ad = sonoi(freq(l,i), freq_ref, a, b)
              freq_corr(l,i) = freq(l,i) + correction_factor*sonoi_ad%val/Q(l,i)
           end do
        end do

      end subroutine get_sonoi_all_freq_corr


      subroutine get_sonoi_freq_corr(radial_only, freq_ref)
         logical, intent(in) :: radial_only
         real(dp), intent(in) :: freq_ref
         call get_sonoi_all_freq_corr(sonoi_a, sonoi_b, radial_only, freq_ref, &
            nl, freq_target, freq_sigma, model_freq, model_freq_corr, model_inertia)
      end subroutine get_sonoi_freq_corr
      
      
      subroutine get_sonoi_freq_corr_alt_up(radial_only, freq_ref)
         logical, intent(in) :: radial_only
         real(dp), intent(in) :: freq_ref
         call get_sonoi_all_freq_corr(sonoi_a, sonoi_b, radial_only, freq_ref, &
            nl, freq_target, freq_sigma, model_freq_alt_up, model_freq_corr_alt_up, model_inertia_alt_up)
      end subroutine get_sonoi_freq_corr_alt_up
      
      
      subroutine get_sonoi_freq_corr_alt_down(radial_only, freq_ref)
         logical, intent(in) :: radial_only
         real(dp), intent(in) :: freq_ref
         call get_sonoi_all_freq_corr(sonoi_a, sonoi_b, radial_only, freq_ref, &
            nl, freq_target, freq_sigma, model_freq_alt_down, model_freq_corr_alt_down, model_inertia_alt_down)
      end subroutine get_sonoi_freq_corr_alt_down
      
      
      subroutine get_freq_corr(s, radial_only, ierr)      
         type (star_info), pointer :: s
         logical, intent(in) :: radial_only
         integer, intent(out) :: ierr
         ierr = 0
         if (s% use_other_astero_freq_corr) then
            call s% other_astero_freq_corr(s% id, ierr)
            return
         end if
         if (correction_scheme == 'kjeldsen') then
            call get_kjeldsen_freq_corr
            surf_coef1 = a_div_r
            surf_coef2 = correction_r
            
            if (save_next_best_at_higher_frequency) &
               call get_kjeldsen_freq_corr_alt_up
            if (save_next_best_at_lower_frequency) &
               call get_kjeldsen_freq_corr_alt_down
            call get_kjeldsen_freq_corr
         else if (correction_scheme == 'cubic') then
            call get_cubic_freq_corr(radial_only)
            surf_coef1 = a3*pow3(5000.*s%nu_max/s% nu_max_sun)
            surf_coef2 = 0
            
            if (save_next_best_at_higher_frequency) &
               call get_cubic_freq_corr_alt_up(radial_only)
            if (save_next_best_at_lower_frequency) &
               call get_cubic_freq_corr_alt_down(radial_only)
            call get_cubic_freq_corr(radial_only)
         else if (correction_scheme == 'combined') then
            call get_combined_freq_corr(radial_only)
            surf_coef1 = a3*pow3(5000.*s%nu_max/s% nu_max_sun)
            surf_coef2 = a1/(5000.*s%nu_max/s% nu_max_sun)
            
            if (save_next_best_at_higher_frequency) &
               call get_combined_freq_corr_alt_up(radial_only)
            if (save_next_best_at_lower_frequency) &
               call get_combined_freq_corr_alt_down(radial_only)
            call get_combined_freq_corr(radial_only)
         else if (correction_scheme == 'sonoi') then
            call get_sonoi_freq_corr(radial_only, s% nu_max)
            surf_coef1 = sonoi_a
            surf_coef2 = sonoi_b
            
            if (save_next_best_at_higher_frequency) &
               call get_sonoi_freq_corr_alt_up(radial_only, s% nu_max)
            if (save_next_best_at_lower_frequency) &
               call get_sonoi_freq_corr_alt_down(radial_only, s% nu_max)
            call get_sonoi_freq_corr(radial_only, s% nu_max)
         else if (correction_scheme == 'power_law') then
            call get_power_law_freq_corr(radial_only, s% nu_max)
            surf_coef1 = power_law_a
            surf_coef2 = power_law_b
            
            if (save_next_best_at_higher_frequency) &
               call get_power_law_freq_corr_alt_up(radial_only, s% nu_max)
            if (save_next_best_at_lower_frequency) &
               call get_power_law_freq_corr_alt_down(radial_only, s% nu_max)
            call get_power_law_freq_corr(radial_only, s% nu_max)
         else 
            call get_no_freq_corr
            surf_coef1 = 0
            surf_coef2 = 0
            
            if (save_next_best_at_higher_frequency) &
               call get_no_freq_corr_alt_up
            if (save_next_best_at_lower_frequency) &
               call get_no_freq_corr_alt_down
            call get_no_freq_corr
         end if
      end subroutine get_freq_corr
      

      ! chi2 = chi2_seismo*chi2_seismo_fraction &
      !      + chi2_spectro*(1 - chi2_seismo_fraction)
      real(dp) function get_chi2(s, max_el, trace_okay, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: max_el
         logical, intent(in) :: trace_okay
         integer, intent(out) :: ierr

         integer :: i, l, n, chi2N1, chi2N2
         real(dp) :: chi2term, Teff, logL, chi2sum1, chi2sum2, frac, &
            model_r01, model_r10, model_r02
         
         ! calculate chi^2 following Brandao et al, 2011, eqn 11
         include 'formats'
         
         ierr = 0
         chi2sum1 = 0
         chi2N1 = 0
         chi2_r_010_ratios = 0
         chi2_r_02_ratios = 0
         chi2_frequencies = 0
         
         if (chi2_seismo_freq_fraction > 0) then
         
            if (trace_okay .and. trace_chi2_seismo_frequencies_info) &
               write(*,'(4a6,99a20)') &
                  'model', 'i', 'l', 'n', 'chi2term', 'freq', 'corr', &
                  'obs', 'sigma', 'log E'

            do i = 1, nl(0)
               if (freq_target(0,i) < 0) cycle
               chi2term = &
                  pow2((model_freq_corr(0,i) - freq_target(0,i))/freq_sigma(0,i))
               if (trace_okay .and. trace_chi2_seismo_frequencies_info) &
                  write(*,'(4i6,99(1pe20.10))') &
                     s% model_number, i, 0, model_order(0,i), chi2term, model_freq(0,i), &
                     model_freq_corr(0,i), freq_target(0,i), freq_sigma(0,i), safe_log10(model_inertia(0,i))
               chi2sum1 = chi2sum1 + chi2term
               chi2N1 = chi2N1 + 1
            end do

            do l = 1, 3
               if (max_el >= l) then
                  do i = 1, nl(l)
                     if (freq_target(l,i) < 0) cycle
                     chi2term = &
                           pow2((model_freq_corr(l,i) - freq_target(l,i))/freq_sigma(l,i))
                     if (trace_okay .and. trace_chi2_seismo_frequencies_info) &
                           write(*,'(4i6,99(1pe20.10))') &
                           s% model_number, i, l, model_order(l,i), chi2term, model_freq(l,i), &
                           model_freq_corr(l,i), freq_target(l,i), freq_sigma(l,i), safe_log10(model_inertia(l,i))
                     chi2sum1 = chi2sum1 + chi2term
                     chi2N1 = chi2N1 + 1
                  end do
               end if
            end do

            num_chi2_seismo_terms = chi2N1
            if (normalize_chi2_seismo_frequencies) then
               chi2_frequencies = chi2sum1/max(1,chi2N1)
            else
               chi2_frequencies = chi2sum1
            end if

         end if
         
         if (chi2_seismo_r_010_fraction > 0 .and. max_el >= 1) then

            if (ratios_n == 0) then
               write(*,*) 'ERROR: chi2_seismo_r_010_fraction > 0 but cannot evaluate r_010'
               ierr = -1
               return
            end if
            
            chi2sum1 = 0
            do i=1,ratios_n
               model_r01 = interpolate_ratio_r010( &
                  freq_target(0,i + ratios_l0_first), ratios_l0_first, model_freq(0,:), model_ratios_r01)
               if (trace_okay .and. trace_chi2_seismo_ratios_info) &
                  write(*,2) 'r01 obs, model, interp, model - interp', &
                     i, ratios_r01(i), model_ratios_r01(i + ratios_l0_first), model_r01, &
                     model_ratios_r01(i + ratios_l0_first) - model_r01
               model_r10 = interpolate_ratio_r010( &
                  freq_target(1,i + ratios_l1_first), ratios_l1_first, model_freq(1,:), model_ratios_r10)
               if (trace_okay .and. trace_chi2_seismo_ratios_info) &
                  write(*,2) 'r10 obs, model, interp, model - interp', &
                     i, ratios_r10(i), model_ratios_r10(i + ratios_l1_first), model_r10, &
                     model_ratios_r10(i + ratios_l1_first) - model_r10
               chi2term = &
                  pow2((model_r01 - ratios_r01(i))/sigmas_r01(i)) + &
                  pow2((model_r10 - ratios_r10(i))/sigmas_r10(i))
               chi2sum1 = chi2sum1 + chi2term
               if (trace_okay .and. trace_chi2_seismo_ratios_info) &
                  write(*,2) 'chi2 ratios terms r01 r10', i, chi2term, &
                     pow2((model_r01 - ratios_r01(i))/sigmas_r01(i)), &
                     pow2((model_r10 - ratios_r10(i))/sigmas_r10(i)), &
                     model_r01, model_r10
            end do
            n = 2*ratios_n
            if (normalize_chi2_seismo_r_010) then
               chi2_r_010_ratios = chi2sum1/max(1,n)
            else
               chi2_r_010_ratios = chi2sum1
            end if
            
         end if
         
         if (chi2_seismo_r_02_fraction > 0 .and. max_el >= 2) then
            
            chi2sum1 = 0
            n = 0
            do i=1,nl(0)
               if (sigmas_r02(i) == 0d0) cycle
               model_r02 = interpolate_ratio_r02( &
                  freq_target(0,i + ratios_l0_first), model_freq(0,:), model_ratios_r02)
               if (trace_okay .and. trace_chi2_seismo_ratios_info) &
                  write(*,2) 'r02 obs, model, interp, model - interp', &
                     i, ratios_r02(i), model_ratios_r02(i), model_r02, &
                     model_ratios_r02(i) - model_r02
               chi2sum1 = chi2sum1 + &
                  pow2((model_r02 - ratios_r02(i))/sigmas_r02(i))
               n = n+1
            end do
            if (n == 0) then
               write(*,*) 'ERROR: chi2_seismo_r_02_fraction > 0 but cannot evaluate r_02'
               ierr = -1
               return
            end if
            if (normalize_chi2_seismo_r_02) then
               chi2_r_02_ratios = chi2sum1/max(1,n)
            else
               chi2_r_02_ratios = chi2sum1
            end if
                              
         end if

         chi2_seismo = &
            chi2_seismo_r_010_fraction*chi2_r_010_ratios + &
            chi2_seismo_r_02_fraction*chi2_r_02_ratios + &
            chi2_seismo_freq_fraction*chi2_frequencies + &
            chi2_seismo_delta_nu_fraction*chi2_delta_nu + &
            chi2_seismo_nu_max_fraction*chi2_nu_max
         
         chi2sum2 = 0
         chi2N2 = 0
         
         if (age_sigma > 0 .and. include_age_in_chi2_spectro) then
            chi2term = pow2((s% star_age - age_target)/age_sigma)
            if (trace_okay .and. trace_chi2_spectro_info) &
               write(*,2) 'chi2_spectro_term age', s% model_number, chi2term
            chi2sum2 = chi2sum2 + chi2term
            chi2N2 = chi2N2 + 1
         end if

         do i = 1, max_constraints
            if (constraint_sigma(i) > 0 .and. include_constraint_in_chi2_spectro(i)) then
               chi2term = pow2((constraint_value(i) - constraint_target(i))/constraint_sigma(i))
               if (trace_okay .and. trace_chi2_spectro_info) &
                  write(*,2) 'chi2_spectro_term ' // trim(constraint_name(i)), s% model_number, chi2term
               chi2sum2 = chi2sum2 + chi2term
               chi2N2 = chi2N2 + 1
            end if
         end do

         num_chi2_spectro_terms = chi2N2
         if (normalize_chi2_spectro) then
            chi2_spectro = chi2sum2/max(1,chi2N2)
         else
            chi2_spectro = chi2sum2
         end if
         
         frac = chi2_seismo_fraction
         chi2 = frac*chi2_seismo + (1-frac)*chi2_spectro         

         get_chi2 = chi2
         
         if (chi2_seismo_fraction < 0 .or. chi2_seismo_fraction > 1) then
            write(*,1) 'ERROR: bad chi2_seismo_fraction', chi2_seismo_fraction
            stop
         end if
         
         !if (is_bad(chi2)) call mesa_error(__FILE__,__LINE__,'get_chi2')
                  
      end function get_chi2


      end module astero_support
