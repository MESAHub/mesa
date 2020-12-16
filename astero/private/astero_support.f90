! ***********************************************************************
!
!   Copyright (C) 2013  Bill Paxton & The MESA Team
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
      
      implicit none


      
      contains
      
      subroutine check_search_controls(ierr)
         integer, intent(out) :: ierr
         integer :: i
         include 'formats'
         ierr = 0         
         do i=2,nl0
            if (l0_obs(i) <= l0_obs(i-1)) then
               write(*,3) 'l0_obs values out of order', i-1, i, l0_obs(i-1), l0_obs(i)
               ierr = -1
            end if
         end do         
         do i=2,nl1
            if (l1_obs(i) <= l1_obs(i-1)) then
               write(*,3) 'l1_obs values out of order', i-1, i, l1_obs(i-1), l1_obs(i)
               ierr = -1
            end if
         end do         
         do i=2,nl2
            if (l2_obs(i) <= l2_obs(i-1)) then
               write(*,3) 'l2_obs values out of order', i-1, i, l2_obs(i-1), l2_obs(i)
               ierr = -1
            end if
         end do         
         do i=2,nl3
            if (l3_obs(i) <= l3_obs(i-1)) then
               write(*,3) 'l3_obs values out of order', i-1, i, l3_obs(i-1), l3_obs(i)
               ierr = -1
            end if
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
               write(*,*)
               write(*,'(a)') 'gyre is not currently enabled in your configuration of mesa.'
               write(*,'(a)') 'check that your utils/makefile_header has USE_GYRE = YES'
               write(*,*)
               return
            end if

            b = correction_b

            num_results = 0
            call do_gyre_get_modes(s, l, store_model, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in do_gyre_get_modes'
               stop 'get_one_el_info'
            end if
         
         else if (code == 'adipls') then 

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
               stop 'get_one_el_info'
            end if
            
         else
         
            write(*,'(a)') 'invalid oscillation_code: ' // trim(oscillation_code)
            ierr = -1
            return
         
         end if

         ! sort results by increasing frequency
         allocate(index(num_results), stat=ierr)
         if (ierr /= 0) then
            stop 'failed in allocate before calling qsort'
         end if
         call qsort(index, num_results, cyclic_freq)

         if (l == 0) then
            call set_to_closest(l0_obs, &
               l0_freq, empty, empty, &
               l0_inertia, empty, empty, &
               l0_order, int_empty, int_empty, ierr)
         else if (l == 1) then
            call set_to_closest(l1_obs, &
               l1_freq, l1_freq_alt_up, l1_freq_alt_down, &
               l1_inertia, l1_inertia_alt_up, l1_inertia_alt_down, &
               l1_order, l1_order_alt_up, l1_order_alt_down, ierr)
         else if (l == 2) then
            call set_to_closest(l2_obs, &
               l2_freq, l2_freq_alt_up, l2_freq_alt_down, &
               l2_inertia, l2_inertia_alt_up, l2_inertia_alt_down, &
               l2_order, l2_order_alt_up, l2_order_alt_down, ierr)
         else if (l == 3) then
            call set_to_closest(l3_obs, &
               l3_freq, l3_freq_alt_up, l3_freq_alt_down, &
               l3_inertia, l3_inertia_alt_up, l3_inertia_alt_down, &
               l3_order, l3_order_alt_up, l3_order_alt_down, ierr)
         else
            stop 'bad value for l in get_one_el_info'
         end if
         if (ierr /= 0) then
            !write(*,2) 'failed to find frequency for matching for l =', l
            write(*,2) 'failed to match frequencies for l =', l
            return
         end if
         
         if (l == 0 .and. correction_factor > 0 .and. nl0 > 0 .and. &
               delta_nu > 0 .and. nu_max > 0 .and. avg_nu_obs > 0) then 
            ! calculate surface correction info
            
            cnt = 0
            sum_1 = 0
            do i=1,nl0
               if (l0_obs(i) < 0) cycle
               cnt = cnt + 1
               sum_1 = sum_1 + l0_freq(i)
            end do
            if (cnt == 0) return
            avg_nu_model = sum_1/cnt
            
            sum_1 = 0
            sum_2 = 0
            sum_3 = 0
            do i=1,nl0
               if (l0_obs(i) < 0) cycle
               sum_1 = sum_1 + &
                  (l0_freq(i) - avg_nu_model)*(l0_n_obs(i) - avg_radial_n)
               sum_2 = sum_2 + pow2(l0_n_obs(i) - avg_radial_n)
               sum_3 = sum_3 + pow(l0_obs(i)/nu_max,b)
            end do
            if (sum_2 == 0 .or. sum_3 == 0) return
            delta_nu_model = sum_1/sum_2
            correction_r = & ! K08 eqn 6
               (b-1)/(b*avg_nu_model/avg_nu_obs - delta_nu_model/delta_nu)
            if (correction_r <= 0) return
            correction_a = & ! K08 eqn 10
               min(0d0, avg_nu_obs - correction_r*avg_nu_model)*nl0/sum_3
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
         real(dp), intent(inout) :: r01(:), r10(:)
         
         integer :: l0_seq_n, l0_last, l1_seq_n, l1_last, i, i0, i1
         real(dp) :: d01, d10, sd01, sd10, dnu, sdnu
         
         logical :: dbg
            
         include 'formats'
         
         dbg = .false.
         
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
            sd01 = sqrt(pow2(l0_obs_sigma(i0-1)) + pow2(4*l1_obs_sigma(i1-1)) + &
               pow2(6*l0_obs_sigma(i0)) + pow2(4*l1_obs_sigma(i1)) + pow2(l0_obs_sigma(i0+1)))/8d0
            dnu = l1(i1) - l1(i1-1)
            sdnu = sqrt(pow2(l1_obs_sigma(i1)) + pow2(l1_obs_sigma(i1-1)))
            sigmas_r01(i) = sqrt(pow2(sd01/dnu) + pow2(sdnu*d01/(dnu*dnu)))
            sd10 = sqrt(pow2(l1_obs_sigma(i1-1)) + pow2(4*l0_obs_sigma(i0)) + &
               pow2(6*l1_obs_sigma(i1)) + pow2(4*l0_obs_sigma(i0+1)) + pow2(l1_obs_sigma(i1+1)))/8d0
            dnu = l0(i0+1) - l0(i0)
            sdnu = sqrt(pow2(l0_obs_sigma(i0+1)) + pow2(l0_obs_sigma(i0)))
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
         real(dp), intent(inout) :: r02(:)
         
         integer :: i, i0, i1, i2, jmin, j
         real(dp) :: d02, sd02, dnu, sdnu, df, f0, f2, fmin, fmax, dfmin
         
         logical :: dbg
            
         include 'formats'
         
         dbg = .false.
         
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
            sd02 = sqrt(pow2(l0_obs_sigma(i0)) + pow2(l2_obs_sigma(i2)))
            sdnu = sqrt(pow2(l1_obs_sigma(i1)) + pow2(l1_obs_sigma(i1-1)))
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
         do i=1,nl0
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
         do i=i_lo+1,nl0
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
         integer, dimension(max_nl0) :: orders
         real(dp) :: sum_1, sum_2, sum_3, range, nmax
         real(dp) :: x, y, isig2, sum_xy, sum_x, sum_y, sum_x2, sum_isig2, d
         
         logical, parameter :: dbg = .false.
         
         include 'formats'
         
         ierr = 0
         
         !call test_get_frequency_ratios
         
         if (nl0 <= 0) return
         
         sigmas_r02 = 0d0
         ratios_r02 = 0d0
         
         if (chi2_seismo_r_010_fraction > 0 .or. &
             chi2_seismo_r_02_fraction > 0) then
            call get_frequency_ratios( &
               .true., nl0, l0_obs, nl1, l1_obs, &
               ratios_n, ratios_l0_first, ratios_l1_first, &
               ratios_r01, ratios_r10)
         end if
         if (chi2_seismo_r_02_fraction > 0) then
            call get_r02_frequency_ratios( &
               .true., nl0, l0_obs, nl1, l1_obs, nl2, l2_obs, ratios_r02)
         end if
         
         if (delta_nu <= 0 .and. nl0 > 1 .and. l0_n_obs(1) > 0) then
            sum_xy = 0
            sum_x = 0
            sum_y = 0
            sum_x2 = 0
            sum_isig2 = 0
            do i=1,nl0
               isig2 = 1d0/pow2(l0_obs_sigma(i))
               x = dble(l0_n_obs(i))
               y = l0_obs(i)
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
            ! set l0_n_obs(i) to order of l0_obs(i)
            range = l0_obs(nl0) - l0_obs(1)
            norders = int(range/delta_nu + 0.5d0) + 1
            nmax = (nu_max/delta_nu)*(delta_nu_sun/nu_max_sun)*22.6 - 1.6         
            l0_n_obs(1) = int(nmax - (norders-1)/2)
            if (dbg) write(*,3) 'l0_n_obs(i)', 1, l0_n_obs(1), l0_obs(1)
            do i=2,norders
               l0_n_obs(i) = l0_n_obs(1) + &
                  int((l0_obs(i) - l0_obs(1))/delta_nu + 0.5)
               if (dbg) write(*,3) 'l0_n_obs(i)', i, l0_n_obs(i), l0_obs(i)
            end do
            if (dbg) then
               write(*,1) 'range', range
               write(*,2) 'norders', norders
               write(*,1) 'nmax', nmax
               write(*,2) '(norders+1)/2', (norders+1)/2
               write(*,2) 'l0_n_obs(1)', l0_n_obs(1)
               write(*,*)
               !stop
            end if
         end if 
         
         
         cnt = 0
         sum_1 = 0
         sum_2 = 0
         do i=1,nl0
            if (l0_obs(i) < 0) cycle
            cnt = cnt + 1
            sum_1 = sum_1 + l0_obs(i)
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
            write(*,*)
            stop 'init_obs_data'
         end if
            
      end subroutine init_obs_data
      
      
      real(dp) function interpolate_l0_inertia(freq) result(inertia)
         real(dp), intent(in) :: freq
         integer :: i
         real(dp) :: alfa, beta
         inertia = 0
         if (nl0 == 0) return
         if (freq <= l0_freq(1)) then
            inertia = l0_inertia(1)
            return
         end if
         if (freq >= l0_freq(nl0)) then
            inertia = l0_inertia(nl0)
            return
         end if
         do i=2,nl0
            if (freq < l0_freq(i)) then
               alfa = (freq - l0_freq(i-1))/(l0_freq(i) - l0_freq(i-1))
               beta = 1d0 - alfa
               inertia = alfa*l0_inertia(i) + beta*l0_inertia(i-1)
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
            if (b > 0 .and. correction_factor > 0 .and. l1_freq_corr(i) > 0 .and. nl0 > 0) then
               interp_l0_inertia = interpolate_l0_inertia(l1_freq(i))
               !write(*,2) 'l1_freq_corr: l0_inertia interp prev', i, interp_l0_inertia, &
               !   (l0_inertia(min(nl0,i)) + l0_inertia(min(nl0,i+1)))/2
               Qnl = l1_inertia(i)/interp_l0_inertia
               l1_freq_corr(i) = l1_freq_corr(i) + &
                  correction_factor*(a_div_r/Qnl)*pow(l1_freq(i)/nu_max,b)
            end if
         end do
      end subroutine get_kjeldsen_nonradial_freq_corr

      
      subroutine get_kjeldsen_nonradial_freq_corr_alt_up
         !stop 'get_kjeldsen_nonradial_freq_corr_alt_up'
      end subroutine get_kjeldsen_nonradial_freq_corr_alt_up

      
      subroutine get_kjeldsen_nonradial_freq_corr_alt_down 
         !stop 'get_kjeldsen_nonradial_freq_corr_alt_down'
      end subroutine get_kjeldsen_nonradial_freq_corr_alt_down

      
      subroutine get_kjeldsen_freq_corr
         call get_kjeldsen_radial_freq_corr( &
            a_div_r, correction_b, nu_max, correction_factor, .true., &
            nl0, l0_obs, l0_freq, l0_freq_corr, l0_inertia)
         call get_kjeldsen_nonradial_freq_corr( &
            a_div_r, correction_b, nu_max, correction_factor, .true., &
            nl1, l1_obs, l1_freq, l1_freq_corr, l1_inertia, l0_inertia)
         call get_kjeldsen_nonradial_freq_corr( &
            a_div_r, correction_b, nu_max, correction_factor, .true., &
            nl2, l2_obs, l2_freq, l2_freq_corr, l2_inertia, l0_inertia)
         call get_kjeldsen_nonradial_freq_corr( &
            a_div_r, correction_b, nu_max, correction_factor, .true., &
            nl3, l3_obs, l3_freq, l3_freq_corr, l3_inertia, l0_inertia)
      end subroutine get_kjeldsen_freq_corr

      
      subroutine get_kjeldsen_freq_corr_alt_up
         call get_kjeldsen_nonradial_freq_corr( &
            a_div_r, correction_b, nu_max, correction_factor, .true., &
            nl1, l1_obs, l1_freq_alt_up, l1_freq_corr_alt_up, l1_inertia_alt_up, l0_inertia)
         call get_kjeldsen_nonradial_freq_corr( &
            a_div_r, correction_b, nu_max, correction_factor, .true., &
            nl2, l2_obs, l2_freq_alt_up, l2_freq_corr_alt_up, l2_inertia_alt_up, l0_inertia)
         call get_kjeldsen_nonradial_freq_corr( &
            a_div_r, correction_b, nu_max, correction_factor, .true., &
            nl3, l3_obs, l3_freq_alt_up, l3_freq_corr_alt_up, l3_inertia_alt_up, l0_inertia)
      end subroutine get_kjeldsen_freq_corr_alt_up

      
      subroutine get_kjeldsen_freq_corr_alt_down 
         call get_kjeldsen_nonradial_freq_corr( &
            a_div_r, correction_b, nu_max, correction_factor, .true., &
            nl1, l1_obs, l1_freq_alt_down, l1_freq_corr_alt_down, l1_inertia_alt_down, l0_inertia)
         call get_kjeldsen_nonradial_freq_corr( &
            a_div_r, correction_b, nu_max, correction_factor, .true., &
            nl2, l2_obs, l2_freq_alt_down, l2_freq_corr_alt_down, l2_inertia_alt_down, l0_inertia)
         call get_kjeldsen_nonradial_freq_corr( &
            a_div_r, correction_b, nu_max, correction_factor, .true., &
            nl3, l3_obs, l3_freq_alt_down, l3_freq_corr_alt_down, l3_inertia_alt_down, l0_inertia)
      end subroutine get_kjeldsen_freq_corr_alt_down
      
      
      subroutine get_no_freq_corr   ! use correction_factor = 0d0 for this
         call get_kjeldsen_radial_freq_corr( &
            a_div_r, correction_b, nu_max, 0d0, .true., &
            nl0, l0_obs, l0_freq, l0_freq_corr, l0_inertia)
         call get_kjeldsen_nonradial_freq_corr( &
            a_div_r, correction_b, nu_max, 0d0, .true., &
            nl1, l1_obs, l1_freq, l1_freq_corr, l1_inertia, l0_inertia)
         call get_kjeldsen_nonradial_freq_corr( &
            a_div_r, correction_b, nu_max, 0d0, .true., &
            nl2, l2_obs, l2_freq, l2_freq_corr, l2_inertia, l0_inertia)
         call get_kjeldsen_nonradial_freq_corr( &
            a_div_r, correction_b, nu_max, 0d0, .true., &
            nl3, l3_obs, l3_freq, l3_freq_corr, l3_inertia, l0_inertia)
      end subroutine get_no_freq_corr
      
      
      subroutine get_no_freq_corr_alt_up
         integer :: i
         do i=1,nl1
            l1_freq_corr_alt_up(i) = l1_freq_alt_up(i)
         end do
         do i=1,nl2
            l2_freq_corr_alt_up(i) = l2_freq_alt_up(i)
         end do
         do i=1,nl3
            l3_freq_corr_alt_up(i) = l3_freq_alt_up(i)
         end do
      end subroutine get_no_freq_corr_alt_up
      
      
      subroutine get_no_freq_corr_alt_down
         integer :: i
         do i=1,nl1
            l1_freq_corr_alt_down(i) = l1_freq_alt_down(i)
         end do
         do i=1,nl2
            l2_freq_corr_alt_down(i) = l2_freq_alt_down(i)
         end do
         do i=1,nl3
            l3_freq_corr_alt_down(i) = l3_freq_alt_down(i)
         end do
      end subroutine get_no_freq_corr_alt_down


      subroutine get_cubic_radial_freq_corr(a3, &
            nl0, l0_obs, l0_obs_sigma, l0_freq, l0_freq_corr, l0_inertia)
        integer, intent(in) :: nl0
        real(dp), intent(in), dimension(:) :: &
             l0_obs, l0_obs_sigma, l0_freq, l0_inertia
        real(dp), intent(inout), dimension(:) :: l0_freq_corr
        real(dp), intent(out) :: a3

        real(dp) :: X, y, XtX, Xty
        integer :: i

        XtX = 0d0
        Xty = 0d0

        do i = 1, nl0
           X = l0_freq(i)**3/l0_inertia(i)/l0_obs_sigma(i)
           y = (l0_obs(i) - l0_freq(i))/l0_obs_sigma(i)

           XtX = XtX + X*X
           Xty = Xty + X*y
        end do

        a3 = Xty/XtX

        l0_freq_corr(1:nl0) = l0_freq(1:nl0) + correction_factor*a3*l0_freq(1:nl0)**3/l0_inertia(1:nl0)

      end subroutine get_cubic_radial_freq_corr


      subroutine get_cubic_all_freq_corr(a3, radial_only, &
            nl0, l0_obs, l0_obs_sigma, l0_freq, l0_freq_corr, l0_inertia, &
            nl1, l1_obs, l1_obs_sigma, l1_freq, l1_freq_corr, l1_inertia, &
            nl2, l2_obs, l2_obs_sigma, l2_freq, l2_freq_corr, l2_inertia, &
            nl3, l3_obs, l3_obs_sigma, l3_freq, l3_freq_corr, l3_inertia)
           
        integer, intent(in) :: nl0, nl1, nl2, nl3
        real(dp), intent(in), dimension(:) :: &
             l0_obs, l0_obs_sigma, l0_freq, l0_inertia, &
             l1_obs, l1_obs_sigma, l1_freq, l1_inertia, &
             l2_obs, l2_obs_sigma, l2_freq, l2_inertia, &
             l3_obs, l3_obs_sigma, l3_freq, l3_inertia
        real(dp), intent(inout), dimension(:) :: &
             l0_freq_corr, l1_freq_corr, l2_freq_corr, l3_freq_corr
        real(dp), intent(out) :: a3

        logical :: radial_only
        real(dp) :: X, y, XtX, Xty
        integer :: i

        XtX = 0d0
        Xty = 0d0

        do i = 1, nl0
           X = l0_freq(i)**3/l0_inertia(i)/l0_obs_sigma(i)
           y = (l0_obs(i) - l0_freq(i))/l0_obs_sigma(i)

           XtX = XtX + X*X
           Xty = Xty + X*y
        end do

        if (.not. radial_only) then
           do i = 1, nl1
              X = l1_freq(i)**3/l1_inertia(i)/l1_obs_sigma(i)
              y = (l1_obs(i) - l1_freq(i))/l1_obs_sigma(i)

              XtX = XtX + X*X
              Xty = Xty + X*y
           end do

           do i = 1, nl2
              X = l2_freq(i)**3/l2_inertia(i)/l2_obs_sigma(i)
              y = (l2_obs(i) - l2_freq(i))/l2_obs_sigma(i)

              XtX = XtX + X*X
              Xty = Xty + X*y
           end do

           do i = 1, nl3
              X = l3_freq(i)**3/l3_inertia(i)/l3_obs_sigma(i)
              y = (l3_obs(i) - l3_freq(i))/l3_obs_sigma(i)

              XtX = XtX + X*X
              Xty = Xty + X*y
           end do
        end if

        a3 = Xty/XtX

        l0_freq_corr(1:nl0) = l0_freq(1:nl0) + correction_factor*a3*l0_freq(1:nl0)**3/l0_inertia(1:nl0)
        l1_freq_corr(1:nl1) = l1_freq(1:nl1) + correction_factor*a3*l1_freq(1:nl1)**3/l1_inertia(1:nl1)
        l2_freq_corr(1:nl2) = l2_freq(1:nl2) + correction_factor*a3*l2_freq(1:nl2)**3/l2_inertia(1:nl2)
        l3_freq_corr(1:nl3) = l3_freq(1:nl3) + correction_factor*a3*l3_freq(1:nl3)**3/l3_inertia(1:nl3)

      end subroutine get_cubic_all_freq_corr


      subroutine get_cubic_freq_corr(radial_only)
         logical, intent(in) :: radial_only
         call get_cubic_all_freq_corr(a3, radial_only, &
            nl0, l0_obs, l0_obs_sigma, l0_freq, l0_freq_corr, l0_inertia, &
            nl1, l1_obs, l1_obs_sigma, l1_freq, l1_freq_corr, l1_inertia, &
            nl2, l2_obs, l2_obs_sigma, l2_freq, l2_freq_corr, l2_inertia, &
            nl3, l3_obs, l3_obs_sigma, l3_freq, l3_freq_corr, l3_inertia)
      end subroutine get_cubic_freq_corr
      
      
      subroutine get_cubic_freq_corr_alt_up(radial_only)
         logical, intent(in) :: radial_only
         call get_cubic_all_freq_corr(a3, radial_only, &
            nl0, l0_obs, l0_obs_sigma, l0_freq, l0_freq_corr, l0_inertia, &
            nl1, l1_obs, l1_obs_sigma, l1_freq_alt_up, l1_freq_corr_alt_up, l1_inertia_alt_up, &
            nl2, l2_obs, l2_obs_sigma, l2_freq_alt_up, l2_freq_corr_alt_up, l2_inertia_alt_up, &
            nl3, l3_obs, l3_obs_sigma, l3_freq_alt_up, l3_freq_corr_alt_up, l3_inertia_alt_up)
      end subroutine get_cubic_freq_corr_alt_up
      
      
      subroutine get_cubic_freq_corr_alt_down(radial_only)
         logical, intent(in) :: radial_only
         call get_cubic_all_freq_corr(a3, radial_only, &
            nl0, l0_obs, l0_obs_sigma, l0_freq, l0_freq_corr, l0_inertia, &
            nl1, l1_obs, l1_obs_sigma, l1_freq_alt_down, l1_freq_corr_alt_down, l1_inertia_alt_down, &
            nl2, l2_obs, l2_obs_sigma, l2_freq_alt_down, l2_freq_corr_alt_down, l2_inertia_alt_down, &
            nl3, l3_obs, l3_obs_sigma, l3_freq_alt_down, l3_freq_corr_alt_down, l3_inertia_alt_down)
      end subroutine get_cubic_freq_corr_alt_down

      
      subroutine get_combined_radial_freq_corr(a3, a1, &
           nl0, l0_obs, l0_obs_sigma, l0_freq, l0_freq_corr, l0_inertia)

        integer, intent(in) :: nl0
        real(dp), intent(in), dimension(:) :: &
             l0_obs, l0_obs_sigma, l0_freq, l0_inertia
        real(dp), intent(inout), dimension(:) :: l0_freq_corr
        real(dp), intent(out) :: a3, a1

        integer :: i
        real(dp) :: X(2), XtX(2,2), XtXi(2,2), Xty(2), y
        real(dp) :: detXtX

        XtX = 0d0
        Xty = 0d0

        do i = 1, nl0
           X(1) = l0_freq(i)**(-1)/l0_inertia(i)/l0_obs_sigma(i)
           X(2) = l0_freq(i)**3/l0_inertia(i)/l0_obs_sigma(i)
           y = (l0_obs(i) - l0_freq(i))/l0_obs_sigma(i)

           XtX(1,1) = XtX(1,1) + X(1)*X(1)
           XtX(1,2) = XtX(1,2) + X(1)*X(2)
           XtX(2,2) = XtX(2,2) + X(2)*X(2)
           Xty(1) = Xty(1) + X(1)*y
           Xty(2) = Xty(2) + X(2)*y
        end do

        XtX(2,1) = XtX(1,2)

        XtXi(1,1) = XtX(2,2)
        XtXi(2,2) = XtX(1,1)
        XtXi(1,2) = -XtX(1,2)
        XtXi(2,1) = -XtX(2,1)

        detXtX = XtX(1,1)*XtX(2,2) - XtX(1,2)*XtX(2,1)
        XtXi = XtXi/detXtX

        a1 = XtXi(1,1)*Xty(1) + XtXi(1,2)*Xty(2)
        a3 = XtXi(2,1)*Xty(1) + XtXi(2,2)*Xty(2)

        l0_freq_corr(1:nl0) = l0_freq(1:nl0) + &
             correction_factor*(a1*l0_freq(1:nl0)**(-1)+a3*l0_freq(1:nl0)**3)/l0_inertia(1:nl0)

      end subroutine get_combined_radial_freq_corr


      subroutine get_combined_all_freq_corr(a3, a1, radial_only, &
            nl0, l0_obs, l0_obs_sigma, l0_freq, l0_freq_corr, l0_inertia, &
            nl1, l1_obs, l1_obs_sigma, l1_freq, l1_freq_corr, l1_inertia, &
            nl2, l2_obs, l2_obs_sigma, l2_freq, l2_freq_corr, l2_inertia, &
            nl3, l3_obs, l3_obs_sigma, l3_freq, l3_freq_corr, l3_inertia)

        integer, intent(in) :: nl0, nl1, nl2, nl3
        real(dp), intent(in), dimension(:) :: &
             l0_obs, l0_obs_sigma, l0_freq, l0_inertia, &
             l1_obs, l1_obs_sigma, l1_freq, l1_inertia, &
             l2_obs, l2_obs_sigma, l2_freq, l2_inertia, &
             l3_obs, l3_obs_sigma, l3_freq, l3_inertia
        real(dp), intent(inout), dimension(:) :: &
             l0_freq_corr, l1_freq_corr, l2_freq_corr, l3_freq_corr
        real(dp), intent(out) :: a3, a1
        logical :: radial_only
        
        integer :: i
        real(dp) :: X(2), XtX(2,2), XtXi(2,2), Xty(2), y
        real(dp) :: detXtX

        XtX = 0d0
        Xty = 0d0

        do i = 1, nl0
           X(1) = l0_freq(i)**(-1)/l0_inertia(i)/l0_obs_sigma(i)
           X(2) = l0_freq(i)**3/l0_inertia(i)/l0_obs_sigma(i)
           y = (l0_obs(i) - l0_freq(i))/l0_obs_sigma(i)

           XtX(1,1) = XtX(1,1) + X(1)*X(1)
           XtX(1,2) = XtX(1,2) + X(1)*X(2)
           XtX(2,2) = XtX(2,2) + X(2)*X(2)
           Xty(1) = Xty(1) + X(1)*y
           Xty(2) = Xty(2) + X(2)*y
        end do

        if (.not. radial_only) then
           do i = 1, nl1
              X(1) = l1_freq(i)**(-1)/l1_inertia(i)/l1_obs_sigma(i)
              X(2) = l1_freq(i)**3/l1_inertia(i)/l1_obs_sigma(i)
              y = (l1_obs(i) - l1_freq(i))/l1_obs_sigma(i)

              XtX(1,1) = XtX(1,1) + X(1)*X(1)
              XtX(1,2) = XtX(1,2) + X(1)*X(2)
              XtX(2,2) = XtX(2,2) + X(2)*X(2)
              Xty(1) = Xty(1) + X(1)*y
              Xty(2) = Xty(2) + X(2)*y
           end do

           do i = 1, nl2
              X(1) = l2_freq(i)**(-1)/l2_inertia(i)/l2_obs_sigma(i)
              X(2) = l2_freq(i)**3/l2_inertia(i)/l2_obs_sigma(i)
              y = (l2_obs(i) - l2_freq(i))/l2_obs_sigma(i)

              XtX(1,1) = XtX(1,1) + X(1)*X(1)
              XtX(1,2) = XtX(1,2) + X(1)*X(2)
              XtX(2,2) = XtX(2,2) + X(2)*X(2)
              Xty(1) = Xty(1) + X(1)*y
              Xty(2) = Xty(2) + X(2)*y
           end do

           do i = 1, nl3
              X(1) = l3_freq(i)**(-1)/l3_inertia(i)/l3_obs_sigma(i)
              X(2) = l3_freq(i)**3/l3_inertia(i)/l3_obs_sigma(i)
              y = (l3_obs(i) - l3_freq(i))/l3_obs_sigma(i)

              XtX(1,1) = XtX(1,1) + X(1)*X(1)
              XtX(1,2) = XtX(1,2) + X(1)*X(2)
              XtX(2,2) = XtX(2,2) + X(2)*X(2)
              Xty(1) = Xty(1) + X(1)*y
              Xty(2) = Xty(2) + X(2)*y
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

        l0_freq_corr(1:nl0) = l0_freq(1:nl0) + &
             correction_factor*(a1*l0_freq(1:nl0)**(-1)+a3*l0_freq(1:nl0)**3)/l0_inertia(1:nl0)
        l1_freq_corr(1:nl1) = l1_freq(1:nl1) + &
             correction_factor*(a1*l1_freq(1:nl1)**(-1)+a3*l1_freq(1:nl1)**3)/l1_inertia(1:nl1)
        l2_freq_corr(1:nl2) = l2_freq(1:nl2) + &
             correction_factor*(a1*l2_freq(1:nl2)**(-1)+a3*l2_freq(1:nl2)**3)/l2_inertia(1:nl2)
        l3_freq_corr(1:nl3) = l3_freq(1:nl3) + &
             correction_factor*(a1*l3_freq(1:nl3)**(-1)+a3*l3_freq(1:nl3)**3)/l3_inertia(1:nl3)
      end subroutine get_combined_all_freq_corr


      subroutine get_combined_freq_corr(radial_only)
         logical, intent(in) :: radial_only
         call get_combined_all_freq_corr(a3, a1, radial_only, &
              nl0, l0_obs, l0_obs_sigma, l0_freq, l0_freq_corr, l0_inertia, &
              nl1, l1_obs, l1_obs_sigma, l1_freq, l1_freq_corr, l1_inertia, &
              nl2, l2_obs, l2_obs_sigma, l2_freq, l2_freq_corr, l2_inertia, &
              nl3, l3_obs, l3_obs_sigma, l3_freq, l3_freq_corr, l3_inertia)
      end subroutine get_combined_freq_corr
      
      
      subroutine get_combined_freq_corr_alt_up(radial_only)
         logical, intent(in) :: radial_only
         call get_combined_all_freq_corr(a3, a1, radial_only, &
            nl0, l0_obs, l0_obs_sigma, l0_freq, l0_freq_corr, l0_inertia, &
            nl1, l1_obs, l1_obs_sigma, l1_freq_alt_up, l1_freq_corr_alt_up, l1_inertia_alt_up, &
            nl2, l2_obs, l2_obs_sigma, l2_freq_alt_up, l2_freq_corr_alt_up, l2_inertia_alt_up, &
            nl3, l3_obs, l3_obs_sigma, l3_freq_alt_up, l3_freq_corr_alt_up, l3_inertia_alt_up)
      end subroutine get_combined_freq_corr_alt_up
      
      
      subroutine get_combined_freq_corr_alt_down(radial_only)
         logical, intent(in) :: radial_only
         call get_combined_all_freq_corr(a3, a1, radial_only, &
            nl0, l0_obs, l0_obs_sigma, l0_freq, l0_freq_corr, l0_inertia, &
            nl1, l1_obs, l1_obs_sigma, l1_freq_alt_down, l1_freq_corr_alt_down, l1_inertia_alt_down, &
            nl2, l2_obs, l2_obs_sigma, l2_freq_alt_down, l2_freq_corr_alt_down, l2_inertia_alt_down, &
            nl3, l3_obs, l3_obs_sigma, l3_freq_alt_down, l3_freq_corr_alt_down, l3_inertia_alt_down)
      end subroutine get_combined_freq_corr_alt_down
      
      
      real(dp) function power_law(freq, freq_ref, a, b)
        real(dp), intent(in) :: freq, freq_ref, a, b
        
        power_law = a*pow(freq/freq_ref, b)
      end function power_law

      
      real(dp) function dpower_law_da(freq, freq_ref, a, b)
        real(dp), intent(in) :: freq, freq_ref, a, b

        dpower_law_da = pow(freq/freq_ref, b)
      end function dpower_law_da


      real(dp) function dpower_law_db(freq, freq_ref, a, b)
        real(dp), intent(in) :: freq, freq_ref, a, b
        
        dpower_law_db = a*pow(freq/freq_ref, b)*log(freq/freq_ref)
      end function dpower_law_db
      
      
      subroutine get_power_law_all_freq_corr(a, b, radial_only, freq_ref, &
            nl0, l0_obs, l0_obs_sigma, l0_freq, l0_freq_corr, l0_inertia, &
            nl1, l1_obs, l1_obs_sigma, l1_freq, l1_freq_corr, l1_inertia, &
            nl2, l2_obs, l2_obs_sigma, l2_freq, l2_freq_corr, l2_inertia, &
            nl3, l3_obs, l3_obs_sigma, l3_freq, l3_freq_corr, l3_inertia)

        integer, intent(in) :: nl0, nl1, nl2, nl3
        real(dp), intent(in) :: freq_ref
        real(dp), intent(in), dimension(:) :: &
             l0_obs, l0_obs_sigma, l0_freq, l0_inertia, &
             l1_obs, l1_obs_sigma, l1_freq, l1_inertia, &
             l2_obs, l2_obs_sigma, l2_freq, l2_inertia, &
             l3_obs, l3_obs_sigma, l3_freq, l3_inertia
        real(dp), intent(out), dimension(:) :: &
             l0_freq_corr, l1_freq_corr, l2_freq_corr, l3_freq_corr
        real(dp), intent(out) :: a, b

        logical :: radial_only
        integer :: i, iter
        real(dp) :: X(2), XtX(2,2), XtXi(2,2), Xty(2), y
        real(dp) :: detXtX, da, db
        real(dp) :: Q1(nl1), Q2(nl2), Q3(nl3)

        ! Power_Law's solar values happen to be the same as MESA's but the
        ! commented expression is there in case it changes
        ! nu_max = s% nu_max ! *(3100d0/3100d0)/sqrt(5777d0/5777d0)

        ! initial guesses are solar-calibrated values from Ball et al. (2016)
        a = -5.25d0
        b = 5.37d0

        do iter=1,1000
           XtX = 0d0
           Xty = 0d0

           do i = 1, nl0
              X(1) = -dpower_law_da(l0_freq(i), freq_ref, a, b)/l0_obs_sigma(i)
              X(2) = -dpower_law_db(l0_freq(i), freq_ref, a, b)/l0_obs_sigma(i)
              y = (l0_obs(i) - l0_freq(i) - power_law(l0_freq(i), freq_ref, a, b))/l0_obs_sigma(i)
              
              XtX(1,1) = XtX(1,1) + X(1)*X(1)
              XtX(1,2) = XtX(1,2) + X(1)*X(2)
              XtX(2,2) = XtX(2,2) + X(2)*X(2)
              Xty(1) = Xty(1) + X(1)*y
              Xty(2) = Xty(2) + X(2)*y
           end do

           if (.not. radial_only) then
              do i = 1, nl1
                 Q1(i) = l1_inertia(i)/interpolate_l0_inertia(l1_freq(i))

                 X(1) = -dpower_law_da(l1_freq(i), freq_ref, a, b)/l1_obs_sigma(i)
                 X(2) = -dpower_law_db(l1_freq(i), freq_ref, a, b)/l1_obs_sigma(i)
                 y = ((l1_obs(i) - l1_freq(i))*Q1(i) - power_law(l1_freq(i), freq_ref, a, b))/l1_obs_sigma(i)
                 
                 XtX(1,1) = XtX(1,1) + X(1)*X(1)
                 XtX(1,2) = XtX(1,2) + X(1)*X(2)
                 XtX(2,2) = XtX(2,2) + X(2)*X(2)
                 Xty(1) = Xty(1) + X(1)*y
                 Xty(2) = Xty(2) + X(2)*y
              end do

              do i = 1, nl2
                 Q2(i) = l2_inertia(i)/interpolate_l0_inertia(l2_freq(i))

                 X(1) = -dpower_law_da(l2_freq(i), freq_ref, a, b)/l2_obs_sigma(i)
                 X(2) = -dpower_law_db(l2_freq(i), freq_ref, a, b)/l2_obs_sigma(i)
                 y = ((l2_obs(i) - l2_freq(i))*Q2(i) - power_law(l2_freq(i), freq_ref, a, b))/l2_obs_sigma(i)
                 
                 XtX(1,1) = XtX(1,1) + X(1)*X(1)
                 XtX(1,2) = XtX(1,2) + X(1)*X(2)
                 XtX(2,2) = XtX(2,2) + X(2)*X(2)
                 Xty(1) = Xty(1) + X(1)*y
                 Xty(2) = Xty(2) + X(2)*y
              end do

              do i = 1, nl3
                 Q3(i) = l3_inertia(i)/interpolate_l0_inertia(l3_freq(i))

                 X(1) = -dpower_law_da(l3_freq(i), freq_ref, a, b)/l3_obs_sigma(i)
                 X(2) = -dpower_law_db(l3_freq(i), freq_ref, a, b)/l3_obs_sigma(i)
                 y = ((l3_obs(i) - l3_freq(i))*Q3(i) - power_law(l3_freq(i), freq_ref, a, b))/l3_obs_sigma(i)
                 
                 XtX(1,1) = XtX(1,1) + X(1)*X(1)
                 XtX(1,2) = XtX(1,2) + X(1)*X(2)
                 XtX(2,2) = XtX(2,2) + X(2)*X(2)
                 Xty(1) = Xty(1) + X(1)*y
                 Xty(2) = Xty(2) + X(2)*y
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
        !    write(*,*) l0_freq(1:nl0)
        !    if (nl1 > 0) write(*,*) l1_freq(1:nl1)
        !    if (nl2 > 0) write(*,*) l2_freq(1:nl2)
        !    if (nl3 > 0) write(*,*) l3_freq(1:nl3)
        ! end if

        do i=1,nl0
           l0_freq_corr(i) = l0_freq(i) + correction_factor*power_law(l0_freq(i), freq_ref, a, b)
        end do
        do i=1,nl1
           l1_freq_corr(i) = l1_freq(i) + correction_factor*power_law(l1_freq(i), freq_ref, a, b)/Q1(i)
        end do
        do i=1,nl2
           l2_freq_corr(i) = l2_freq(i) + correction_factor*power_law(l2_freq(i), freq_ref, a, b)/Q2(i)
        end do
        do i=1,nl3
           l3_freq_corr(i) = l3_freq(i) + correction_factor*power_law(l3_freq(i), freq_ref, a, b)/Q3(i)
        end do

      end subroutine get_power_law_all_freq_corr


      subroutine get_power_law_freq_corr(radial_only, freq_ref)
         logical, intent(in) :: radial_only
         real(dp), intent(in) :: freq_ref
         call get_power_law_all_freq_corr(power_law_a, power_law_b, radial_only, freq_ref, &
              nl0, l0_obs, l0_obs_sigma, l0_freq, l0_freq_corr, l0_inertia, &
              nl1, l1_obs, l1_obs_sigma, l1_freq, l1_freq_corr, l1_inertia, &
              nl2, l2_obs, l2_obs_sigma, l2_freq, l2_freq_corr, l2_inertia, &
              nl3, l3_obs, l3_obs_sigma, l3_freq, l3_freq_corr, l3_inertia)
      end subroutine get_power_law_freq_corr
      
      
      subroutine get_power_law_freq_corr_alt_up(radial_only, freq_ref)
         logical, intent(in) :: radial_only
         real(dp), intent(in) :: freq_ref
         call get_power_law_all_freq_corr(power_law_a, power_law_b, radial_only, freq_ref, &
            nl0, l0_obs, l0_obs_sigma, l0_freq, l0_freq_corr, l0_inertia, &
            nl1, l1_obs, l1_obs_sigma, l1_freq_alt_up, l1_freq_corr_alt_up, l1_inertia_alt_up, &
            nl2, l2_obs, l2_obs_sigma, l2_freq_alt_up, l2_freq_corr_alt_up, l2_inertia_alt_up, &
            nl3, l3_obs, l3_obs_sigma, l3_freq_alt_up, l3_freq_corr_alt_up, l3_inertia_alt_up)
      end subroutine get_power_law_freq_corr_alt_up
      
      
      subroutine get_power_law_freq_corr_alt_down(radial_only, freq_ref)
         logical, intent(in) :: radial_only
         real(dp), intent(in) :: freq_ref
         call get_power_law_all_freq_corr(power_law_a, power_law_b, radial_only, freq_ref, &
            nl0, l0_obs, l0_obs_sigma, l0_freq, l0_freq_corr, l0_inertia, &
            nl1, l1_obs, l1_obs_sigma, l1_freq_alt_down, l1_freq_corr_alt_down, l1_inertia_alt_down, &
            nl2, l2_obs, l2_obs_sigma, l2_freq_alt_down, l2_freq_corr_alt_down, l2_inertia_alt_down, &
            nl3, l3_obs, l3_obs_sigma, l3_freq_alt_down, l3_freq_corr_alt_down, l3_inertia_alt_down)
      end subroutine get_power_law_freq_corr_alt_down
      
      
      real(dp) function sonoi(freq, freq_ref, a, b)
        real(dp), intent(in) :: freq, freq_ref, a, b
        
        sonoi = a*freq_ref*(1d0 - 1d0/(1d0+pow(freq/freq_ref,b)))
      end function sonoi


      real(dp) function dsonoi_da(freq, freq_ref, a, b)
        real(dp), intent(in) :: freq, freq_ref, a, b
        
        dsonoi_da = freq_ref*(1d0 - 1d0/(1d0+pow(freq/freq_ref,b)))
      end function dsonoi_da


      real(dp) function dsonoi_db(freq, freq_ref, a, b)
        real(dp), intent(in) :: freq, freq_ref, a, b
        
        dsonoi_db = a*freq_ref*pow(freq/freq_ref,b)*log(freq/freq_ref)/(1d0+pow(freq/freq_ref,b))**2d0
      end function dsonoi_db


      subroutine get_sonoi_all_freq_corr(a, b, radial_only, freq_ref, &
            nl0, l0_obs, l0_obs_sigma, l0_freq, l0_freq_corr, l0_inertia, &
            nl1, l1_obs, l1_obs_sigma, l1_freq, l1_freq_corr, l1_inertia, &
            nl2, l2_obs, l2_obs_sigma, l2_freq, l2_freq_corr, l2_inertia, &
            nl3, l3_obs, l3_obs_sigma, l3_freq, l3_freq_corr, l3_inertia)

        integer, intent(in) :: nl0, nl1, nl2, nl3
        real(dp), intent(in) :: freq_ref
        real(dp), intent(in), dimension(:) :: &
             l0_obs, l0_obs_sigma, l0_freq, l0_inertia, &
             l1_obs, l1_obs_sigma, l1_freq, l1_inertia, &
             l2_obs, l2_obs_sigma, l2_freq, l2_inertia, &
             l3_obs, l3_obs_sigma, l3_freq, l3_inertia
        real(dp), intent(out), dimension(:) :: &
             l0_freq_corr, l1_freq_corr, l2_freq_corr, l3_freq_corr
        real(dp), intent(out) :: a, b

        logical :: radial_only
        integer :: i, iter
        real(dp) :: X(2), XtX(2,2), XtXi(2,2), Xty(2), y
        real(dp) :: detXtX, da, db
        real(dp) :: Q1(nl1), Q2(nl2), Q3(nl3)

        ! Sonoi's solar values happen to be the same as MESA's but the
        ! commented expression is there in case it changes
        ! nu_max = s% nu_max ! *(3100d0/3100d0)/sqrt(5777d0/5777d0)

        a = -3.59d-3
        b = 11.26d0

        do iter=1,1000
           XtX = 0d0
           Xty = 0d0

           do i = 1, nl0
              X(1) = -dsonoi_da(l0_freq(i), freq_ref, a, b)/l0_obs_sigma(i)
              X(2) = -dsonoi_db(l0_freq(i), freq_ref, a, b)/l0_obs_sigma(i)
              y = (l0_obs(i) - l0_freq(i) - sonoi(l0_freq(i), freq_ref, a, b))/l0_obs_sigma(i)
              
              XtX(1,1) = XtX(1,1) + X(1)*X(1)
              XtX(1,2) = XtX(1,2) + X(1)*X(2)
              XtX(2,2) = XtX(2,2) + X(2)*X(2)
              Xty(1) = Xty(1) + X(1)*y
              Xty(2) = Xty(2) + X(2)*y
           end do

           if (.not. radial_only) then
              do i = 1, nl1
                 Q1(i) = l1_inertia(i)/interpolate_l0_inertia(l1_freq(i))

                 X(1) = -dsonoi_da(l1_freq(i), freq_ref, a, b)/l1_obs_sigma(i)
                 X(2) = -dsonoi_db(l1_freq(i), freq_ref, a, b)/l1_obs_sigma(i)
                 y = ((l1_obs(i) - l1_freq(i))*Q1(i) - sonoi(l1_freq(i), freq_ref, a, b))/l1_obs_sigma(i)
                 
                 XtX(1,1) = XtX(1,1) + X(1)*X(1)
                 XtX(1,2) = XtX(1,2) + X(1)*X(2)
                 XtX(2,2) = XtX(2,2) + X(2)*X(2)
                 Xty(1) = Xty(1) + X(1)*y
                 Xty(2) = Xty(2) + X(2)*y
              end do

              do i = 1, nl2
                 Q2(i) = l2_inertia(i)/interpolate_l0_inertia(l2_freq(i))

                 X(1) = -dsonoi_da(l2_freq(i), freq_ref, a, b)/l2_obs_sigma(i)
                 X(2) = -dsonoi_db(l2_freq(i), freq_ref, a, b)/l2_obs_sigma(i)
                 y = ((l2_obs(i) - l2_freq(i))*Q2(i) - sonoi(l2_freq(i), freq_ref, a, b))/l2_obs_sigma(i)
                 
                 XtX(1,1) = XtX(1,1) + X(1)*X(1)
                 XtX(1,2) = XtX(1,2) + X(1)*X(2)
                 XtX(2,2) = XtX(2,2) + X(2)*X(2)
                 Xty(1) = Xty(1) + X(1)*y
                 Xty(2) = Xty(2) + X(2)*y
              end do

              do i = 1, nl3
                 Q3(i) = l3_inertia(i)/interpolate_l0_inertia(l3_freq(i))

                 X(1) = -dsonoi_da(l3_freq(i), freq_ref, a, b)/l3_obs_sigma(i)
                 X(2) = -dsonoi_db(l3_freq(i), freq_ref, a, b)/l3_obs_sigma(i)
                 y = ((l3_obs(i) - l3_freq(i))*Q3(i) - sonoi(l3_freq(i), freq_ref, a, b))/l3_obs_sigma(i)
                 
                 XtX(1,1) = XtX(1,1) + X(1)*X(1)
                 XtX(1,2) = XtX(1,2) + X(1)*X(2)
                 XtX(2,2) = XtX(2,2) + X(2)*X(2)
                 Xty(1) = Xty(1) + X(1)*y
                 Xty(2) = Xty(2) + X(2)*y
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
        !    write(*,*) l0_freq(1:nl0)
        !    if (nl1 > 0) write(*,*) l1_freq(1:nl1)
        !    if (nl2 > 0) write(*,*) l2_freq(1:nl2)
        !    if (nl3 > 0) write(*,*) l3_freq(1:nl3)
        ! end if

        if (b < 0d0) a = 0d0

        do i=1,nl0
           l0_freq_corr(i) = l0_freq(i) + correction_factor*sonoi(l0_freq(i), freq_ref, a, b)
        end do
        do i=1,nl1
           l1_freq_corr(i) = l1_freq(i) + correction_factor*sonoi(l1_freq(i), freq_ref, a, b)/Q1(i)
        end do
        do i=1,nl2
           l2_freq_corr(i) = l2_freq(i) + correction_factor*sonoi(l2_freq(i), freq_ref, a, b)/Q2(i)
        end do
        do i=1,nl3
           l3_freq_corr(i) = l3_freq(i) + correction_factor*sonoi(l3_freq(i), freq_ref, a, b)/Q3(i)
        end do

      end subroutine get_sonoi_all_freq_corr


      subroutine get_sonoi_freq_corr(radial_only, freq_ref)
         logical, intent(in) :: radial_only
         real(dp), intent(in) :: freq_ref
         call get_sonoi_all_freq_corr(sonoi_a, sonoi_b, radial_only, freq_ref, &
              nl0, l0_obs, l0_obs_sigma, l0_freq, l0_freq_corr, l0_inertia, &
              nl1, l1_obs, l1_obs_sigma, l1_freq, l1_freq_corr, l1_inertia, &
              nl2, l2_obs, l2_obs_sigma, l2_freq, l2_freq_corr, l2_inertia, &
              nl3, l3_obs, l3_obs_sigma, l3_freq, l3_freq_corr, l3_inertia)
      end subroutine get_sonoi_freq_corr
      
      
      subroutine get_sonoi_freq_corr_alt_up(radial_only, freq_ref)
         logical, intent(in) :: radial_only
         real(dp), intent(in) :: freq_ref
         call get_sonoi_all_freq_corr(sonoi_a, sonoi_b, radial_only, freq_ref, &
            nl0, l0_obs, l0_obs_sigma, l0_freq, l0_freq_corr, l0_inertia, &
            nl1, l1_obs, l1_obs_sigma, l1_freq_alt_up, l1_freq_corr_alt_up, l1_inertia_alt_up, &
            nl2, l2_obs, l2_obs_sigma, l2_freq_alt_up, l2_freq_corr_alt_up, l2_inertia_alt_up, &
            nl3, l3_obs, l3_obs_sigma, l3_freq_alt_up, l3_freq_corr_alt_up, l3_inertia_alt_up)
      end subroutine get_sonoi_freq_corr_alt_up
      
      
      subroutine get_sonoi_freq_corr_alt_down(radial_only, freq_ref)
         logical, intent(in) :: radial_only
         real(dp), intent(in) :: freq_ref
         call get_sonoi_all_freq_corr(sonoi_a, sonoi_b, radial_only, freq_ref, &
            nl0, l0_obs, l0_obs_sigma, l0_freq, l0_freq_corr, l0_inertia, &
            nl1, l1_obs, l1_obs_sigma, l1_freq_alt_down, l1_freq_corr_alt_down, l1_inertia_alt_down, &
            nl2, l2_obs, l2_obs_sigma, l2_freq_alt_down, l2_freq_corr_alt_down, l2_inertia_alt_down, &
            nl3, l3_obs, l3_obs_sigma, l3_freq_alt_down, l3_freq_corr_alt_down, l3_inertia_alt_down)
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

         integer :: i, n, chi2N1, chi2N2
         real(dp) :: chi2term, Teff, logL, chi2sum1, chi2sum2, frac, &
            model_r01, model_r10, model_r02
         
         ! calculate chi^2 following Brandao et al, 2011, eqn 11
         include 'formats'
         
         ierr = 0
         chi2sum1 = 0
         chi2N1 = 0     
         chi2_r_010_ratios = -1    
         chi2_r_02_ratios = -1   
         chi2_frequencies = 0
         
         if (chi2_seismo_freq_fraction > 0) then
         
            if (trace_okay .and. trace_chi2_seismo_frequencies_info) &
               write(*,'(4a6,99a20)') &
                  'model', 'i', 'l', 'n', 'chi2term', 'freq', 'corr', &
                  'obs', 'sigma', 'log E'
            do i = 1, nl0
               if (l0_obs(i) < 0) cycle
               chi2term = &
                  pow2((l0_freq_corr(i) - l0_obs(i))/l0_obs_sigma(i))
               if (trace_okay .and. trace_chi2_seismo_frequencies_info) &
                  write(*,'(4i6,99(1pe20.10))') &
                     s% model_number, i, 0, l0_order(i), chi2term, l0_freq(i), &
                     l0_freq_corr(i), l0_obs(i), l0_obs_sigma(i), safe_log10(l0_inertia(i))
               chi2sum1 = chi2sum1 + chi2term
               chi2N1 = chi2N1 + 1
            end do
         
            if (max_el >= 1) then
               do i = 1, nl1
                  if (l1_obs(i) < 0) cycle
                  chi2term = &
                     pow2((l1_freq_corr(i) - l1_obs(i))/l1_obs_sigma(i))       
                  if (trace_okay .and. trace_chi2_seismo_frequencies_info) &
                     write(*,'(4i6,99(1pe20.10))') &
                        s% model_number, i, 1, l1_order(i), chi2term, l1_freq(i), &
                        l1_freq_corr(i), l1_obs(i), l1_obs_sigma(i), safe_log10(l1_inertia(i))
                  chi2sum1 = chi2sum1 + chi2term
                  chi2N1 = chi2N1 + 1
               end do
            end if
         
            if (max_el >= 2) then
               do i = 1, nl2
                  if (l2_obs(i) < 0) cycle
                  chi2term = &
                     pow2((l2_freq_corr(i) - l2_obs(i))/l2_obs_sigma(i))            
                  if (trace_okay .and. trace_chi2_seismo_frequencies_info) &
                     write(*,'(4i6,99(1pe20.10))') &
                        s% model_number, i, 2, l2_order(i), chi2term, l2_freq(i), &
                        l2_freq_corr(i), l2_obs(i), l2_obs_sigma(i), safe_log10(l2_inertia(i))
                  chi2sum1 = chi2sum1 + chi2term
                  chi2N1 = chi2N1 + 1
               end do
            end if
         
            if (max_el >= 3) then
               do i = 1, nl3
                  if (l3_obs(i) < 0) cycle
                  chi2term = &
                     pow2((l3_freq_corr(i) - l3_obs(i))/l3_obs_sigma(i))            
                  if (trace_okay .and. trace_chi2_seismo_frequencies_info) &
                     write(*,'(4i6,99(1pe20.10))') &
                        s% model_number, i, 3, l3_order(i), chi2term, l3_freq(i), &
                        l3_freq_corr(i), l3_obs(i), l3_obs_sigma(i), safe_log10(l3_inertia(i))
                  chi2sum1 = chi2sum1 + chi2term
                  chi2N1 = chi2N1 + 1
               end do
            end if
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
                  l0_obs(i + ratios_l0_first), ratios_l0_first, l0_freq, model_ratios_r01)
               if (trace_okay .and. trace_chi2_seismo_ratios_info) &
                  write(*,2) 'r01 obs, model, interp, model - interp', &
                     i, ratios_r01(i), model_ratios_r01(i + ratios_l0_first), model_r01, &
                     model_ratios_r01(i + ratios_l0_first) - model_r01
               model_r10 = interpolate_ratio_r010( &
                  l1_obs(i + ratios_l1_first), ratios_l1_first, l1_freq, model_ratios_r10)
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
            do i=2,nl0
               if (sigmas_r02(i) == 0d0) cycle
               model_r02 = interpolate_ratio_r02( &
                  l0_obs(i + ratios_l0_first), l0_freq, model_ratios_r02)
               if (trace_okay .and. trace_chi2_seismo_ratios_info) &
                  write(*,2) 'r02 obs, model, interp, model - interp', &
                     i, ratios_r02(i), model_ratios_r02(i), model_r02, &
                     model_ratios_r02(i) - model_r10
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
         
         if (Teff_sigma > 0 .and. include_Teff_in_chi2_spectro) then
            Teff = s% Teff
            chi2term = pow2((Teff - Teff_target)/Teff_sigma)
            if (trace_okay .and. trace_chi2_spectro_info) &
               write(*,2) 'chi2_spectro_term Teff', s% model_number, chi2term
            chi2sum2 = chi2sum2 + chi2term
            chi2N2 = chi2N2 + 1
         end if
         
         if (logL_sigma > 0 .and. include_logL_in_chi2_spectro) then
            logL = s% log_surface_luminosity
            chi2term = pow2((logL - logL_target)/logL_sigma)
            if (trace_okay .and. trace_chi2_spectro_info) &
               write(*,2) 'chi2_spectro_term logL', s% model_number, chi2term
            chi2sum2 = chi2sum2 + chi2term
            chi2N2 = chi2N2 + 1
         end if
         
         if (logg_sigma > 0 .and. include_logg_in_chi2_spectro) then
            chi2term = pow2((logg - logg_target)/logg_sigma)
            if (trace_okay .and. trace_chi2_spectro_info) &
               write(*,2) 'chi2_spectro_term logg', s% model_number, chi2term
            chi2sum2 = chi2sum2 + chi2term
            chi2N2 = chi2N2 + 1
         end if
         
         if (FeH_sigma > 0 .and. include_FeH_in_chi2_spectro) then
            chi2term = pow2((FeH - FeH_target)/FeH_sigma)
            if (trace_okay .and. trace_chi2_spectro_info) &
               write(*,2) 'chi2_spectro_term FeH', s% model_number, chi2term
            chi2sum2 = chi2sum2 + chi2term
            chi2N2 = chi2N2 + 1
         end if
         
         if (logR_sigma > 0 .and. include_logR_in_chi2_spectro) then
            chi2term = pow2((logR - logR_target)/logR_sigma)
            if (trace_okay .and. trace_chi2_spectro_info) &
               write(*,2) 'chi2_spectro_term logR', s% model_number, chi2term
            chi2sum2 = chi2sum2 + chi2term
            chi2N2 = chi2N2 + 1
         end if
         
         if (age_sigma > 0 .and. include_age_in_chi2_spectro) then
            chi2term = pow2((s% star_age - age_target)/age_sigma)
            if (trace_okay .and. trace_chi2_spectro_info) &
               write(*,2) 'chi2_spectro_term age', s% model_number, chi2term
            chi2sum2 = chi2sum2 + chi2term
            chi2N2 = chi2N2 + 1
         end if
         
         if (surface_Z_div_X_sigma > 0 .and. include_surface_Z_div_X_in_chi2_spectro) then
            chi2term = pow2((surface_Z_div_X - surface_Z_div_X_target)/surface_Z_div_X_sigma)
            if (trace_okay .and. trace_chi2_spectro_info) &
               write(*,2) 'chi2_spectro_term surface_Z_div_X', s% model_number, chi2term
            chi2sum2 = chi2sum2 + chi2term
            chi2N2 = chi2N2 + 1
         end if
         
         if (surface_He_sigma > 0 .and. include_surface_He_in_chi2_spectro) then
            chi2term = pow2((surface_He - surface_He_target)/surface_He_sigma)
            if (trace_okay .and. trace_chi2_spectro_info) &
               write(*,2) 'chi2_spectro_term surface_He', s% model_number, chi2term
            chi2sum2 = chi2sum2 + chi2term
            chi2N2 = chi2N2 + 1
         end if
         
         if (Rcz_sigma > 0 .and. include_Rcz_in_chi2_spectro) then
            chi2term = pow2((Rcz - Rcz_target)/Rcz_sigma)
            if (trace_okay .and. trace_chi2_spectro_info) &
               write(*,2) 'chi2_spectro_term Rcz', s% model_number, chi2term
            chi2sum2 = chi2sum2 + chi2term
            chi2N2 = chi2N2 + 1
         end if
         
         if (my_var1_sigma > 0 .and. include_my_var1_in_chi2_spectro) then
            chi2term = pow2((my_var1 - my_var1_target)/my_var1_sigma)
            if (trace_okay .and. trace_chi2_spectro_info) &
               write(*,2) 'chi2_spectro_term ' // trim(my_var1_name), s% model_number, chi2term
            chi2sum2 = chi2sum2 + chi2term
            chi2N2 = chi2N2 + 1
         end if
         
         if (my_var2_sigma > 0 .and. include_my_var2_in_chi2_spectro) then
            chi2term = pow2((my_var2 - my_var2_target)/my_var2_sigma)
            if (trace_okay .and. trace_chi2_spectro_info) &
               write(*,2) 'chi2_spectro_term ' // trim(my_var2_name), s% model_number, chi2term
            chi2sum2 = chi2sum2 + chi2term
            chi2N2 = chi2N2 + 1
         end if
         
         if (my_var3_sigma > 0 .and. include_my_var3_in_chi2_spectro) then
            chi2term = pow2((my_var3 - my_var3_target)/my_var3_sigma)
            if (trace_okay .and. trace_chi2_spectro_info) &
               write(*,2) 'chi2_spectro_term ' // trim(my_var3_name), s% model_number, chi2term
            chi2sum2 = chi2sum2 + chi2term
            chi2N2 = chi2N2 + 1
         end if

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
         
         !if (is_bad(chi2)) stop 'get_chi2'
                  
      end function get_chi2


      end module astero_support
