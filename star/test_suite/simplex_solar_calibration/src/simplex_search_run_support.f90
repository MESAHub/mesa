! ***********************************************************************
!
!   Copyright (C) 2020  Bill Paxton and MESA Team
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

      module simplex_search_run_support

      use star_lib
      use star_def
      use math_lib
      use utils_lib
      use const_def
      use simplex_search_data

      
      implicit none
      
      
      contains        

               
      subroutine do_run_star_simplex( &
            extras_controls, inlist_simplex_search_controls_fname)
         use run_star_support, only: do_read_star_job, id_from_read_star_job
         interface
            subroutine extras_controls(id, ierr)
               integer, intent(in) :: id
               integer, intent(out) :: ierr
            end subroutine extras_controls      
         end interface
         character (len=256) :: inlist_simplex_search_controls_fname
         optional inlist_simplex_search_controls_fname

         integer :: id, ierr
         character (len=256) :: inlist_fname
         
         include 'formats'

         ierr = 0
         call do_read_star_job('inlist', ierr) ! this does alloc_star
         ! and saves the id in id_from_read_star_job
         if (ierr /= 0) stop 1

         id = id_from_read_star_job
         id_from_read_star_job = 0
         star_id = id

         okay_to_restart = .true.
         
         call init_simplex_search_data(ierr)
         if (ierr /= 0) stop 1
         
         star_simplex_procs% extras_controls => extras_controls

         if (present(inlist_simplex_search_controls_fname)) then
            inlist_simplex_fname = inlist_simplex_search_controls_fname
         else
            inlist_simplex_fname = 'inlist_simplex_search_controls'
         end if
         write(*,*) 'read ' // trim(inlist_simplex_fname)
         call read_simplex_search_controls(inlist_simplex_fname, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in read_simplex_search_controls'
            stop 1
         end if
         
         if (Y_depends_on_Z .and. vary_Y) then
            vary_Y = .false.
            write(*,*) &
               'WARNING: vary_Y has been changed to false since Y_depends_on_Z is true'
         end if

         if (save_controls) then
            call write_simplex_search_controls(save_controls_filename, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in write_simplex_search_controls'
               stop 1
            end if
         end if
         
         next_Y_to_try = -1
         next_FeH_to_try = -1
         next_mass_to_try = -1      
         next_alpha_to_try = -1
         next_f_ov_to_try = -1 
         sample_number = 0
         max_num_samples = 0
         num_chi2_too_big = 0
         avg_age_top_samples = 1d99
         avg_age_sigma = 1d99
         avg_model_number_top_samples = 1d99
         avg_model_number_sigma = 1d99
         nvar = 0
         my_var1 = 0d0
         my_var2 = 0d0
         my_var3 = 0d0
         my_param1 = 0
         my_param2 = 0
         my_param3 = 0
         
         call init_sample_ptrs

         if (just_do_first_values) then
            vary_Y = .false.
            vary_FeH = .false.
            vary_mass = .false.
            vary_alpha = .false.
            vary_f_ov = .false.
            vary_my_param1 = .false.
            vary_my_param2 = .false.
            vary_my_param3 = .false.
            chi2 = eval1(id,ierr)
         else
            call do_simplex(ierr)
         end if

      end subroutine do_run_star_simplex
      
      
      subroutine do_simplex(ierr)
         use num_lib
         integer, intent(out) :: ierr
         
         real(dp) :: final_mass, final_alpha, final_Y, final_FeH
         real(dp), dimension(:), pointer :: x_first, x_lower, x_upper, x_final
         real(dp), pointer :: simplex(:,:), f(:)
         real(dp) :: f_final
         integer :: lrpar, lipar
         integer, pointer :: ipar(:) ! (lipar)
         real(dp), pointer :: rpar(:) ! (lrpar)
         integer :: num_iters, num_fcn_calls, &
            num_fcn_calls_for_ars, num_accepted_for_ars
         integer :: seed, i, j, k, num_samples
         logical :: start_from_given_simplex_and_f
         
         include 'formats'
         
         ierr = 0
         
         if (vary_mass) then
            nvar = nvar+1; i_mass = nvar
            if (min_mass >= max_mass) then
               write(*,1) 'min_mass >= max_mass', &
                  min_mass, max_mass
               ierr = -1
            end if
         end if

         if (vary_Y) then
            nvar = nvar+1; i_Y = nvar
            if (min_Y >= max_Y) then
               write(*,1) 'min_Y >= max_Y', min_Y, max_Y
               ierr = -1
            end if
         end if
         
         if (vary_FeH) then
            nvar = nvar+1; i_FeH = nvar
            if (min_FeH >= max_FeH) then
               write(*,1) 'min_FeH >= max_FeH', &
                  min_FeH, max_FeH
               ierr = -1
            end if
         end if
         
         if (vary_alpha) then
            nvar = nvar+1; i_alpha = nvar
            if (min_alpha >= max_alpha) then
               write(*,1) 'min_alpha >= max_alpha', &
                  min_alpha, max_alpha
               ierr = -1
            end if
         end if
         
         if (vary_f_ov) then
            nvar = nvar+1; i_f_ov = nvar
            if (min_f_ov >= max_f_ov) then
               write(*,1) 'min_f_ov >= max_f_ov', &
                  min_f_ov, max_f_ov
               ierr = -1
            end if
         end if
         
         if (vary_my_param1) then
            nvar = nvar+1; i_my_param1 = nvar
            if (min_my_param1 >= max_my_param1) then
               write(*,1) 'min_my_param1 >= max_my_param1', &
                  min_my_param1, max_my_param1
               ierr = -1
            end if
         end if
         
         if (vary_my_param2) then
            nvar = nvar+1; i_my_param2 = nvar
            if (min_my_param2 >= max_my_param2) then
               write(*,1) 'min_my_param2 >= max_my_param2', &
                  min_my_param2, max_my_param2
               ierr = -1
            end if
         end if
         
         if (vary_my_param3) then
            nvar = nvar+1; i_my_param3 = nvar
            if (min_my_param3 >= max_my_param3) then
               write(*,1) 'min_my_param3 >= max_my_param3', &
                  min_my_param3, max_my_param3
               ierr = -1
            end if
         end if

         if (ierr /= 0) return
         
         lrpar = 0; lipar = 0

         allocate( &
            rpar(lrpar), ipar(lipar), simplex(nvar,nvar+1), f(nvar+1), &
            x_lower(nvar), x_upper(nvar), x_first(nvar), x_final(nvar))
         
         if (.not. scale_simplex_params) then
            call set_xs
         else ! values are scaled to -1..1 with first at 0
            x_lower(1:nvar) = -1
            x_upper(1:nvar) = 1
            x_first(1:nvar) = 0
         end if
                  
         if (restart_simplex_from_file) then
            call read_samples_from_file(simplex_output_filename, ierr)
            if (ierr /= 0) return
            if (sample_number < nvar+1) then
               write(*,2) 'sorry: too few points. for simplex restart need at least', nvar+1
               ierr = -1
               return
            end if
            num_samples = sample_number
            call setup_simplex_and_f(ierr)
            if (ierr /= 0) return            
            start_from_given_simplex_and_f = .true.
            call set_sample_averages
         else
            start_from_given_simplex_and_f = .false.
         end if
         
         call NM_simplex( &
            nvar, x_lower, x_upper, x_first, x_final, f_final, &
            simplex, f, start_from_given_simplex_and_f, simplex_f, &
            simplex_x_atol, simplex_x_rtol, &
            simplex_itermax, simplex_fcn_calls_max, &
            simplex_centroid_weight_power, simplex_enforce_bounds, &
            simplex_adaptive_random_search, simplex_seed, &
            simplex_alpha, simplex_beta, &
            simplex_gamma, simplex_delta, &
            lrpar, rpar, lipar, ipar, &
            num_iters, num_fcn_calls, &
            num_fcn_calls_for_ars, num_accepted_for_ars, ierr)
         
         if (vary_Y) &
            final_Y = simplex_param( &
               x_final(i_Y), first_Y, min_Y, max_Y)

         if (vary_FeH) &
            final_FeH = simplex_param( &
               x_final(i_FeH), first_FeH, &
               min_FeH, max_FeH)
         
         if (vary_mass) &
            final_mass = simplex_param( &
               x_final(i_mass), first_mass, &
               min_mass, max_mass)
         
         if (vary_alpha) &
            final_alpha = simplex_param( &
               x_final(i_alpha), first_alpha, &
               min_alpha, max_alpha)
         
         if (vary_f_ov) &
            final_f_ov = simplex_param( &
               x_final(i_f_ov), first_f_ov, &
               min_f_ov, max_f_ov)
         
         if (vary_my_param1) &
            final_my_param1 = simplex_param( &
               x_final(i_my_param1), first_my_param1, &
               min_my_param1, max_my_param1)
         
         if (vary_my_param2) &
            final_my_param2 = simplex_param( &
               x_final(i_my_param2), first_my_param2, &
               min_my_param2, max_my_param2)
         
         if (vary_my_param3) &
            final_my_param3 = simplex_param( &
               x_final(i_my_param3), first_my_param3, &
               min_my_param3, max_my_param3)

         deallocate( &
            rpar, ipar, simplex, f, x_lower, x_upper, x_first, x_final)
            
            
         contains
         
         
         subroutine set_xs ! x_first, x_lower, x_upper
            if (vary_Y) then
               x_first(i_Y) = first_Y
               x_lower(i_Y) = min_Y
               x_upper(i_Y) = max_Y
            end if
            if (vary_FeH) then
               x_first(i_FeH) = first_FeH
               x_lower(i_FeH) = min_FeH
               x_upper(i_FeH) = max_FeH
            end if
            if (vary_mass) then
               x_first(i_mass) = first_mass
               x_lower(i_mass) = min_mass
               x_upper(i_mass) = max_mass
            end if
            if (vary_alpha) then
               x_first(i_alpha) = first_alpha
               x_lower(i_alpha) = min_alpha
               x_upper(i_alpha) = max_alpha
            end if
            if (vary_f_ov) then
               x_first(i_f_ov) = first_f_ov
               x_lower(i_f_ov) = min_f_ov
               x_upper(i_f_ov) = max_f_ov
            end if         
            if (vary_my_param1) then
               x_first(i_my_param1) = first_my_param1
               x_lower(i_my_param1) = min_my_param1
               x_upper(i_my_param1) = max_my_param1
            end if         
            if (vary_my_param2) then
               x_first(i_my_param2) = first_my_param2
               x_lower(i_my_param2) = min_my_param2
               x_upper(i_my_param2) = max_my_param2
            end if         
            if (vary_my_param3) then
               x_first(i_my_param3) = first_my_param3
               x_lower(i_my_param3) = min_my_param3
               x_upper(i_my_param3) = max_my_param3
            end if         
         end subroutine set_xs
         
         
         subroutine setup_simplex_and_f(ierr)
            use num_lib, only: qsort
            integer, intent(out) :: ierr
            
            integer :: j, i, k, max_i, jj
            integer, pointer :: index(:)
            ! sort results by increasing sample_chi2
            
            include 'formats'
            
            ierr = 0
            allocate(index(num_samples), stat=ierr)
            if (ierr /= 0) then
               stop 'failed in allocate before calling qsort from show_all_sample_results'
            end if
            call qsort(index, num_samples, sample_chi2)
            max_i = 0
            do j=1,nvar+1
               i = index(j)
               if (i > max_i) max_i = i ! max sample restored
               write(*,3) 'restore simplex', j, i
               f(j) = sample_chi2(i)
               write(*,3) 'chi2', j, i, f(j)
               if (vary_Y) then
                  simplex(i_Y,j) = &
                     simplex_inverse(sample_init_Y(i), first_Y, min_Y, max_Y)
                  write(*,3) 'Y', j, i, sample_init_Y(i)
               end if
               if (vary_FeH) then
                  simplex(i_FeH,j) = &
                     simplex_inverse(sample_init_FeH(i), first_FeH, min_FeH, max_FeH)
                  write(*,3) 'FeH', j, i, sample_init_FeH(i)
               end if         
               if (vary_mass) then
                  simplex(i_mass,j) = &
                     simplex_inverse(sample_mass(i), first_mass, min_mass, max_mass)
                  write(*,3) 'mass', j, i, sample_mass(i)
               end if        
               if (vary_alpha) then
                  simplex(i_alpha,j) = &
                     simplex_inverse(sample_alpha(i), first_alpha, min_alpha, max_alpha)
                  write(*,3) 'alpha', j, i, sample_alpha(i)
               end if         
               if (vary_f_ov) then
                  simplex(i_f_ov,j) = &
                     simplex_inverse(sample_f_ov(i), first_f_ov, min_f_ov, max_f_ov)
                  write(*,3) 'f_ov', j, i, sample_f_ov(i)
               end if
               if (vary_my_param1) then
                  simplex(i_my_param1,j) = &
                     simplex_inverse(sample_my_param1(i), first_my_param1, min_my_param1, max_my_param1)
                  write(*,3) 'my_param1', j, i, sample_my_param1(i)
               end if
               if (vary_my_param2) then
                  simplex(i_my_param2,j) = &
                     simplex_inverse(sample_my_param2(i), first_my_param2, min_my_param2, max_my_param2)
                  write(*,3) 'my_param2', j, i, sample_my_param2(i)
               end if
               if (vary_my_param3) then
                  simplex(i_my_param3,j) = &
                     simplex_inverse(sample_my_param3(i), first_my_param3, min_my_param3, max_my_param3)
                  write(*,3) 'my_param3', j, i, sample_my_param3(i)
               end if
                write(*,*)
            end do
            
            deallocate(index)

            write(*,2) 'num_samples', max_i
            j = 0
            do i=max_i+1,num_samples
               j = j+1
               do jj=1,j
                  write(*,'(a)',advance='no') ' damn, '
               end do
               write(*,'(a)',advance='no') 'will have to rerun sample '
               if (i < 10) then
                  write(*,'(i1)') i
               else if (i < 100) then
                  write(*,'(i2)') i
               else if (i < 1000) then
                  write(*,'(i3)') i
               else if (i < 10000) then
                  write(*,'(i4)') i
               else if (i < 100000) then
                  write(*,'(i5)') i
               else
                  write(*,'(i6)') i
               end if
            end do
            write(*,*)
            write(*,*)
            num_samples = max_i
            
         end subroutine setup_simplex_and_f
         
      end subroutine do_simplex


      ! for nelder-mead minimization
      real(dp) function simplex_f( &
            n, x, lrpar, rpar, lipar, ipar, op_code, ierr)
         use const_def, only: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: x(:) ! (n)
         integer, intent(in) :: lrpar, lipar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(in) :: op_code
         integer, intent(out) :: ierr
         
         integer :: prev_sample_number
         include 'formats'
         
         ierr = 0
         
         write(*,*)
         write(*,*)
         
         if (vary_Y) then
            next_Y_to_try = simplex_param( &
               x(i_Y), first_Y, min_Y, max_Y)
            write(*,1) 'next_Y_to_try', next_Y_to_try, x(i_Y)
            if (next_Y_to_try < 0) then
               write(*,1) 'ERROR: bad arg from simplex -- try again'
               simplex_f = 1d99
               return
            end if
         end if

         if (vary_FeH) then
            next_FeH_to_try = simplex_param( &
               x(i_FeH), first_FeH, min_FeH, max_FeH)
            write(*,1) 'next_FeH_to_try', &
               next_FeH_to_try, x(i_FeH)
         end if
         
         if (vary_mass) then
            next_mass_to_try = simplex_param(&
               x(i_mass), first_mass, min_mass, max_mass)
            write(*,1) 'next_mass_to_try', &
               next_mass_to_try, x(i_mass)
            if (next_mass_to_try < 0) then
               write(*,1) 'ERROR: bad arg from simplex -- try again'
               simplex_f = 1d99
               return
            end if
         end if
         
         if (vary_alpha) then
            next_alpha_to_try = simplex_param( &
               x(i_alpha), first_alpha, min_alpha, max_alpha)
            write(*,1) 'next_alpha_to_try', &
               next_alpha_to_try, x(i_alpha)
            if (next_alpha_to_try < 0) then
               write(*,1) 'ERROR: bad arg from simplex -- try again'
               simplex_f = 1d99
               return
            end if
         end if
         
         if (vary_f_ov) then
            next_f_ov_to_try = simplex_param( &
               x(i_f_ov), first_f_ov, min_f_ov, max_f_ov)
            write(*,1) 'next_f_ov_to_try', &
               next_f_ov_to_try, x(i_f_ov)
            if (next_f_ov_to_try < 0) then
               write(*,1) 'ERROR: bad arg from simplex -- try again'
               simplex_f = 1d99
               return
            end if
         end if
         
         if (vary_my_param1) then
            next_my_param1_to_try = simplex_param( &
               x(i_my_param1), first_my_param1, min_my_param1, max_my_param1)
            write(*,1) 'next_my_param1_to_try', &
               next_my_param1_to_try, x(i_my_param1)
            if (next_my_param1_to_try < 0) then
               write(*,1) 'ERROR: bad arg from simplex -- try again'
               simplex_f = 1d99
               return
            end if
         end if
         
         if (vary_my_param2) then
            next_my_param2_to_try = simplex_param( &
               x(i_my_param2), first_my_param2, min_my_param2, max_my_param2)
            write(*,1) 'next_my_param2_to_try', &
               next_my_param2_to_try, x(i_my_param2)
            if (next_my_param2_to_try < 0) then
               write(*,1) 'ERROR: bad arg from simplex -- try again'
               simplex_f = 1d99
               return
            end if
         end if
         
         if (vary_my_param3) then
            next_my_param3_to_try = simplex_param( &
               x(i_my_param3), first_my_param3, min_my_param3, max_my_param3)
            write(*,1) 'next_my_param3_to_try', &
               next_my_param3_to_try, x(i_my_param3)
            if (next_my_param3_to_try < 0) then
               write(*,1) 'ERROR: bad arg from simplex -- try again'
               simplex_f = 1d99
               return
            end if
         end if
         
         prev_sample_number = sample_number
         simplex_f = eval1(star_id, ierr)
         if (ierr /= 0) then
            write(*,*) 'got ierr from eval1 in simplex -- try again'
            ierr = 0
            simplex_f = 1d99
            return
         end if
         
         if (sample_number == prev_sample_number) then
            write(*,*) 'failed to get new chi^2 -- try again'
            simplex_f = 1d99
            return
         end if
         
         call save_best_for_sample(sample_number, op_code)

         call save_sample_results_to_file(-1, simplex_output_filename, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in save_sample_results_to_file'
            stop 1
         end if

         if (simplex_f < simplex_chi2_tol) then
            write(*,*) 'chi2 < simplex_chi2_tol; stopping further iteration'
            ierr = -1
            return
         endif
         
      end function simplex_f

      
      real(dp) function simplex_param(x, first, min, max) result(param)
         real(dp), intent(in) :: x, first, min, max
         if (.not. scale_simplex_params) then
            param = x
         else
            if (x > 0) then
               param = first + x*(max-first)
            else
               param = first + x*(first-min)
            end if
         end if
      end function simplex_param

      
      real(dp) function simplex_inverse(param, first, min, max) result(x)
         real(dp), intent(in) :: param, first, min, max
         if (.not. scale_simplex_params) then
            x = param
         else
            if (param > first) then
               if (max == first) then
                  x = 1d0
               else
                  x = (param-first)/(max-first)
               end if
            else
               if (first == min) then
                  x = -1d0
               else
                  x = (param-first)/(first-min)
               end if
            end if
         end if
      end function simplex_inverse

      
      subroutine save_best_for_sample(i, op_code)
         integer, intent(in) :: i, op_code
         integer :: ierr
                
         if (i <= 0) return
         if (i > max_num_samples) then
            call alloc_sample_ptrs(ierr)
            if (ierr /= 0) then
               write(*,*) 'ERROR -- failed to allocate for samples'
               stop 'save_best_for_sample'
               return
            end if
         end if
                  
         sample_op_code(i) = op_code
         sample_chi2(i) = best_chi2
         
         sample_age(i) = best_age
         sample_init_Y(i) = current_Y
         sample_init_FeH(i) = current_FeH
         sample_mass(i) = current_mass
         sample_alpha(i) = current_alpha
         sample_f_ov(i) = current_f_ov

         sample_init_h1(i) = current_h1
         sample_init_he3(i) = current_he3
         sample_init_he4(i) = current_he4
         sample_init_Z(i) = current_Z

         sample_radius(i) = best_radius
         sample_logL(i) = best_logL
         sample_Teff(i) = best_Teff
         sample_logg(i) = best_logg
         sample_FeH(i) = best_FeH
         
         sample_logR(i) = best_logR
         sample_surface_Z_div_X(i) = best_surface_Z_div_X
         sample_surface_He(i) = best_surface_He
         sample_Rcz(i) = best_Rcz
         sample_solar_cs_rms(i) = best_solar_cs_rms
         sample_my_var1(i) = best_my_var1
         sample_my_var2(i) = best_my_var2
         sample_my_var3(i) = best_my_var3
         sample_my_param1(i) = best_my_param1
         sample_my_param2(i) = best_my_param2
         sample_my_param3(i) = best_my_param3
         sample_model_number(i) = best_model_number

         call set_sample_averages

      end subroutine save_best_for_sample
      
      
      subroutine set_sample_averages
         integer :: ierr, jj, j, n
         real(dp) :: avg_age_top_samples2, avg_model_number_top_samples2
         
         include 'formats'
      
         call set_sample_index_by_chi2
         n = min(sample_number, max_num_samples_for_avg)
         if (n < max(2,min_num_samples_for_avg)) then
            avg_age_top_samples = 1d99
            avg_age_sigma = 1d99
            avg_model_number_top_samples = 1d99
            avg_model_number_sigma = 1d99
            return
         end if
         avg_age_top_samples = 0
         avg_model_number_top_samples = 0
         avg_age_top_samples2 = 0
         avg_model_number_top_samples2 = 0
         do jj=1,n
            j = sample_index_by_chi2(jj)
            !write(*,3) 'j, jj', j, jj
            avg_age_top_samples = &
               avg_age_top_samples + sample_age(j)
            avg_age_top_samples2 = &
               avg_age_top_samples2 + sample_age(j)*sample_age(j)
            avg_model_number_top_samples = &
               avg_model_number_top_samples + sample_model_number(j)
            avg_model_number_top_samples2 = &
               avg_model_number_top_samples2 + &
               sample_model_number(j)*sample_model_number(j)
            !write(*,2) 'avg_age_top_samples', j, avg_age_top_samples/j
            !write(*,2) 'avg_age_top_samples2', j, avg_age_top_samples2/j
            !write(*,2) 'avg_model_number_top_samples', j, avg_model_number_top_samples/j
            !write(*,2) 'avg_model_number_top_samples2', j, avg_model_number_top_samples2/j
         end do
         avg_age_sigma = &
            sqrt(max(0d0,(avg_age_top_samples2 - avg_age_top_samples*avg_age_top_samples/n)/(n-1)))
         avg_age_top_samples = avg_age_top_samples/n
         avg_model_number_sigma = &
            sqrt(max(0d0,(avg_model_number_top_samples2 - &
                  avg_model_number_top_samples*avg_model_number_top_samples/n)/(n-1)))
         avg_model_number_top_samples = avg_model_number_top_samples/n
         
         write(*,*)
         write(*,2) 'n for averages', n
         write(*,1) 'avg_age_top_samples', avg_age_top_samples
         write(*,1) 'avg_age_sigma', avg_age_sigma
         write(*,1) 'age limit', avg_age_top_samples + avg_age_sigma_limit*avg_age_sigma
         write(*,1) 'avg_model_number_top_samples', avg_model_number_top_samples
         write(*,1) 'avg_model_number_sigma', avg_model_number_sigma
         write(*,1) 'model number limit', &
            avg_model_number_top_samples + &
               avg_model_number_sigma_limit*avg_model_number_sigma
         write(*,*)
         !stop 'set_sample_averages'
         
      end subroutine set_sample_averages
      
      
      subroutine zero_best_info
         best_chi2 = 0
         best_init_h1 = 0
         best_init_he3 = 0
         best_init_he4 = 0
         best_init_Z = 0
         best_age = 0
         best_radius = 0
         best_logL = 0
         best_Teff = 0
         best_logg = 0
         best_FeH = 0
         best_logR = 0
         best_surface_Z_div_X = 0
         best_surface_He = 0
         best_Rcz = 0
         best_solar_cs_rms = 0
         best_my_var1 = 0
         best_my_var2 = 0
         best_my_var3 = 0
         best_my_param1 = 0
         best_my_param2 = 0
         best_my_param3 = 0
      end subroutine zero_best_info

      
      real(dp) function eval1(id_in,ierr)         
         use run_star_support, only: run1_star
         
         integer, intent(in) :: id_in
         integer, intent(out) :: ierr
         
         logical, parameter :: &
            do_alloc_star = .false., &
            do_free_star = .false.
            
         type (star_info), pointer :: s
         logical :: restart
         integer :: id, i      

         include 'formats'
         
         ierr = 0
         id = id_in
         
         call star_ptr(id, s, ierr) 
         if (ierr /= 0) return

         eval1 = -1
         
         ! init for start of run
         best_chi2 = -1    
         num_chi2_too_big = 0     
         simplex_max_dt_next = 1d99
         
         call run1_star( &
            do_alloc_star, do_free_star, &   ! note that these are both false
            okay_to_restart, &
            id, restart, &
            simplex_extras_controls, &
            ierr)
         if (ierr /= 0) return
         
         s% max_years_for_timestep = initial_max_years_for_timestep
         simplex_using_revised_max_yr_dt = .false.
         simplex_revised_max_yr_dt = s% max_years_for_timestep
                  
         okay_to_restart = .false. ! only allow restart on 1st call to run1_star
         
         eval1 = best_chi2
         
         if (simplex_just_call_my_extras_check_model) return
         
         if (best_chi2 < 0) then
            write(*,*) 'failed to find chi^2 for this run'
            call zero_best_info
            best_chi2 = 999999d0
            return
         end if
         
         sample_number = sample_number + 1
         write(*,*)         
         call show_best(6)
         
         if (write_best_model_data_for_each_sample) &
            call write_best(sample_number)
         
      end function eval1
      

      subroutine simplex_extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         real(dp) :: X, Y, Z, FeH, f_ov, a, b, c
         type (star_info), pointer :: s
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
          
         call star_simplex_procs% extras_controls(id, ierr)
         if (ierr /= 0) return

         initial_max_years_for_timestep = s% max_years_for_timestep
         
         if (simplex_just_call_my_extras_check_model) return

         ! overwrite various inlist controls

         if (vary_alpha) then
            s% mixing_length_alpha = next_alpha_to_try
         else
            s% mixing_length_alpha = first_alpha
         end if

         if (vary_my_param1) then
            call star_simplex_procs% will_set_my_param(id, 1, next_my_param1_to_try, ierr)
            if (ierr /= 0) return
            my_param1 = next_my_param1_to_try
         else
            my_param1 = first_my_param1
         end if

         if (vary_my_param2) then
            call star_simplex_procs% will_set_my_param(id, 2, next_my_param2_to_try, ierr)
            if (ierr /= 0) return
            my_param2 = next_my_param2_to_try
         else
            my_param2 = first_my_param2
         end if

         if (vary_my_param3) then
            call star_simplex_procs% will_set_my_param(id, 3, next_my_param3_to_try, ierr)
            if (ierr /= 0) return
            my_param3 = next_my_param3_to_try
         else
            my_param3 = first_my_param3
         end if

         if (vary_f_ov) then
            f_ov = next_f_ov_to_try
         else
            f_ov = first_f_ov
         end if
      
         if (vary_FeH) then
            FeH = next_FeH_to_try
         else
            FeH = first_FeH
         end if

         initial_FeH = FeH
         initial_Z_div_X = Z_div_X_solar*exp10(FeH)

         if (Y_depends_on_Z) then
            a = initial_Z_div_X
            b = dYdZ
            c = 1d0 + a*(1d0 + b)
            X = (1d0 - Y0)/c
            Y = (Y0 + a*(b + Y0))/c
            Z = 1d0 - (X + Y)
         else 
            if (vary_Y) then
               Y = next_Y_to_try
            else
               Y = first_Y
            end if
            X = (1d0 - Y)/(1d0 + initial_Z_div_X)
            Z = X*initial_Z_div_X
         end if

         if (vary_mass) then
            s% job% new_mass = next_mass_to_try
         else
            s% job% new_mass = first_mass
         end if
         
         s% job% relax_initial_mass = .true.
         s% initial_mass = s% job% new_mass
         
         initial_Y = Y
         !s% initial_Z = Z << don't do this. it interferes with use of zams file.
         
         s% job% initial_h1 = X
         s% job% initial_h2 = 0
         s% job% initial_he3 = Y_frac_he3*Y
         s% job% initial_he4 = Y - s% job% initial_he3
         s% job% set_uniform_initial_composition = .true. 
         
         current_Y = Y
         current_FeH = FeH
         current_mass = s% job% new_mass
         current_alpha = s% mixing_length_alpha
         current_f_ov = f_ov
         
         current_h1 = X
         current_he3 = s% job% initial_he3
         current_he4 = s% job% initial_he4
         current_Z = Z
         
         if (f_ov /= 0d0) then
            s% overshoot_scheme(1) = 'exponential'
            s% overshoot_zone_type(1) = 'any'
            s% overshoot_zone_loc(1) = 'any'
            s% overshoot_bdy_loc(1) = 'any'
            s% overshoot_f(1) = f_ov
            s% overshoot_f0(1) = f0_ov_div_f_ov*f_ov
         end if
         
         s% extras_check_model => simplex_extras_check_model
         s% extras_finish_step => simplex_extras_finish_step
         s% extras_after_evolve => simplex_extras_after_evolve
           
      end subroutine simplex_extras_controls


      integer function simplex_extras_check_model(id)            
         integer, intent(in) :: id
         integer :: other_check, ierr
         type (star_info), pointer :: s
         
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         if (simplex_just_call_my_extras_check_model) then
            simplex_extras_check_model = star_simplex_procs% extras_check_model(id)
            best_chi2 = 0
         else
            other_check = star_simplex_procs% extras_check_model(id)
            simplex_extras_check_model = &
                  do_simplex_extras_check_model(s, id)
            if (other_check > simplex_extras_check_model) &
               simplex_extras_check_model = other_check
         end if
         
         star_model_number = s% model_number
               
      end function simplex_extras_check_model

      
      integer function simplex_extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         simplex_extras_finish_step = star_simplex_procs% extras_finish_step(id)         
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         s% dt_next = min(s% dt_next, simplex_max_dt_next)
      end function simplex_extras_finish_step

      
      subroutine simplex_extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         integer :: iounit, ckm
         type (star_info), pointer :: s
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if ((eval_chi2_at_target_age_only .and. s% star_age >= age_target) .or. &
             (include_age_in_chi2 .and. s% star_age >= max_age_for_chi2) .or. &
             save_info_for_last_model) then
            !write(*,*) 'call do_simplex_extras_check_model before terminate'
            ckm = do_simplex_extras_check_model(s, id)
            !write(*,*) 'done do_simplex_extras_check_model before terminate'
         end if         
         call star_simplex_procs% extras_after_evolve(id, ierr)
         if (save_info_for_last_model) then
            write(*,1) 'chi2', chi2
            call store_best_info(s)
            iounit = alloc_iounit(ierr)
            if (ierr /= 0) return
            open(unit=iounit, file=trim(last_model_save_info_filename), &
               action='write', status='replace', iostat=ierr)
            if (ierr /= 0) then
               write(*,'(a)') 'failed to open last_model_save_info_filename ' // &
                  trim(last_model_save_info_filename)
               return
            end if
            write(*,*) 'write ' // trim(last_model_save_info_filename)
            write(*,*) 'call show_best'
            call show_best(iounit)
            write(*,*) 'done show_best'
            call free_iounit(iounit)
         end if
      end subroutine simplex_extras_after_evolve


      integer function do_simplex_extras_check_model(s, id)

         type (star_info), pointer :: s
         integer, intent(in) :: id
         
         integer :: ierr, i, j, n
         logical :: store_model, checking_age
         real(dp) :: age_limit, model_limit, err, X, Y, Z, &
            surface_X, surface_Z, remaining_years, prev_max_years, min_max
         
         include 'formats'
         
         do_simplex_extras_check_model = keep_going
         simplex_max_dt_next = 1d99
         chi2 = -1
         FeH = -1
         checking_age = &
            eval_chi2_at_target_age_only .or. include_age_in_chi2
         
         if (checking_age) then
            if (num_smaller_steps_before_age_target <= 0 .or. &
                dt_for_smaller_steps_before_age_target <= 0) then
               write(*,*) 'ERROR: must set num_smaller_steps_before_age_target'
               write(*,*) 'and dt_for_smaller_steps_before_age_target'
               stop 1
            end if
            if (age_target > s% star_age) then
               remaining_years = age_target - s% star_age
               if (simplex_using_revised_max_yr_dt) &
                  s% max_years_for_timestep = simplex_revised_max_yr_dt
               n = floor(remaining_years/s% max_years_for_timestep + 1d-6)
               j = num_smaller_steps_before_age_target
               if (remaining_years <= s% max_years_for_timestep) then
                  if (mod(s% model_number, s% terminal_interval) == 0) &
                     write(*,'(a40,i6,f20.10)') '(age_target - star_age)/age_sigma', &
                        s% model_number, (age_target - s% star_age)/age_sigma
                  s% max_years_for_timestep = remaining_years
                  simplex_using_revised_max_yr_dt = .true.
                  simplex_revised_max_yr_dt = s% max_years_for_timestep
                  simplex_max_dt_next = s% max_years_for_timestep*secyer
               else if (n <= j) then
                  if (mod(s% model_number, s% terminal_interval) == 0) &
                     write(*,'(a40,i6,f20.10)') '(age_target - star_age)/age_sigma', &
                        s% model_number, (age_target - s% star_age)/age_sigma
                  prev_max_years = s% max_years_for_timestep
                  i = floor(remaining_years/dt_for_smaller_steps_before_age_target + 1d-6)
                  if ((i+1d-9)*dt_for_smaller_steps_before_age_target < remaining_years) then
                     s% max_years_for_timestep = remaining_years/(i+1)
                  else
                     s% max_years_for_timestep = remaining_years/i
                  end if
                  min_max = prev_max_years*s% reduction_factor_for_max_timestep
                  if (s% max_years_for_timestep < min_max) &
                     s% max_years_for_timestep = min_max
                  if (.not. simplex_using_revised_max_yr_dt) then
                     simplex_using_revised_max_yr_dt = .true.
                     write(*,2) 'begin reducing max timestep prior to age target', &
                        s% model_number, remaining_years
                  else if (mod(s% model_number, s% terminal_interval) == 0) then
                     if (simplex_revised_max_yr_dt > s% max_years_for_timestep) then
                        write(*,2) 'reducing max timestep prior to age target', &
                           s% model_number, remaining_years
                     else if (s% max_years_for_timestep <= dt_for_smaller_steps_before_age_target) then
                        i = floor(remaining_years/s% max_years_for_timestep + 1d-6)
                        write(*,3) 'remaining steps and years until age target', &
                           s% model_number, i, remaining_years
                     else 
                        write(*,2) 'remaining_years until age target', &
                           s% model_number, remaining_years
                     end if
                  end if
                  simplex_revised_max_yr_dt = s% max_years_for_timestep
                  if (s% dt_next/secyer > s% max_years_for_timestep) &
                     simplex_max_dt_next = s% max_years_for_timestep*secyer
               end if
            else if (include_age_in_chi2) then
               if (abs(s% max_years_for_timestep - dt_for_smaller_steps_before_age_target) > &
                     dt_for_smaller_steps_before_age_target*1d-2) then
                  write(*,'(a40,i6,f20.10)') '(age_target - star_age)/age_sigma', &
                     s% model_number, (age_target - s% star_age)/age_sigma
                  write(*,1) 'dt_for_smaller_steps_before_age_target', &
                     dt_for_smaller_steps_before_age_target
                  write(*,1) 'max_years_for_timestep', &
                     s% max_years_for_timestep
                  stop 'bad max_years_for_timestep'
               else if (mod(s% model_number, s% terminal_interval) == 0) then
                  write(*,'(a40,i6,f20.10)') '(age_target - star_age)/age_sigma', &
                     s% model_number, (age_target - s% star_age)/age_sigma
               end if
            end if
         else if (s% star_age < min_age_for_chi2) then
            return
         end if

         if (include_age_in_chi2 .and. s% star_age < min_age_for_chi2) return         
         if (eval_chi2_at_target_age_only .and. s% star_age < age_target) return

         if (s% L_nuc_burn_total < s% L_phot*Lnuc_div_L_limit .or. &
               s% star_age < min_age_limit) then
            return
         end if
         
         if (.not. checking_age) then
         
            age_limit = avg_age_top_samples + avg_age_sigma_limit*avg_age_sigma
            if (s% star_age > age_limit) then
               write(*,1) 'star age > limit from top samples', s% star_age, age_limit
               write(*,1) 'avg_age_top_samples', avg_age_top_samples
               write(*,1) 'avg_age_sigma_limit', avg_age_sigma_limit
               write(*,1) 'avg_age_sigma', avg_age_sigma
               do_simplex_extras_check_model = terminate
               return
            end if
         
            model_limit = &
               avg_model_number_top_samples + &
                  avg_model_number_sigma_limit*avg_model_number_sigma
            if (dble(s% model_number) > model_limit) then
               write(*,2) 'model number > limit from top samples', &
                  s% model_number, model_limit
               write(*,1) 'avg_model_number_top_samples', avg_model_number_top_samples
               write(*,1) 'avg_model_number_sigma_limit', avg_model_number_sigma_limit
               write(*,1) 'avg_model_number_sigma', avg_model_number_sigma
               do_simplex_extras_check_model = terminate
               return
            end if

         end if
         
         surface_X = max(s% surface_h1, 1d-10)
         surface_He = s% surface_he3 + s% surface_he4
         surface_Z = max(1d-99, min(1d0, 1 - (surface_X + surface_He)))
         surface_Z_div_X = surface_Z/surface_X
         FeH = log10((surface_Z_div_X)/Z_div_X_solar)
         logg = log10(s% grav(1))
         logR = log10(s% photosphere_r)
         if (.not. include_Rcz_in_chi2) then
            Rcz = 0
         else
            do i = 1, s% nz-1 ! locate bottom of solar convective zone
               if (s% mixing_type(i+1) /= convective_mixing &
                     .and. s% mixing_type(i) == convective_mixing) then
                  if (s% r(i+1) > 0.25*Rsun .and. s% r(i) < 0.9*Rsun) then
                     Rcz = s% r(i)/Rsun
                     exit
                  end if
               end if
            end do
         end if
         if (include_solar_cs_rms_in_chi2 .or. report_solar_cs_rms) then
            solar_cs_rms = calc_current_rms(s, s% nz)
         else
            solar_cs_rms = 0
         end if
         
         call check_limits
         if (do_simplex_extras_check_model /= keep_going) return
         
         chi2 = get_chi2(s, .true., ierr)
         if (ierr /= 0) then
            write(*,'(a40,i6)') 'failed to calculate chi^2', s% model_number
            call check_too_many_bad
            return
         end if

         if (is_bad(chi2)) then
            write(*,1) 'bad chi2', chi2
            write(*,1) 'FeH', FeH
            write(*,1) 'surface_Z', surface_Z
            write(*,1) 'surface_X', surface_X
            write(*,1) 'Z_div_X_solar', Z_div_X_solar
            chi2 = 1d99
            do_simplex_extras_check_model = terminate
            return
            !stop
         end if

         if (chi2 > chi2_limit) then
            write(*,'(a50,i6,99f16.2)') 'chi2 > limit', &
                  s% model_number, chi2, chi2_limit
            call check_too_many_bad
            return
         end if

         write(*,'(a50,i6,99f16.2)') 'chi^2', s% model_number, chi2

         store_model = .true.
          
         if (checking_age) then
            ! leave max_years_for_timestep as is
         else if (chi2 <= chi2_limit_for_smallest_timesteps) then
            s% max_years_for_timestep = max_yrs_dt_chi2_smallest_limit 
            if (s% dt > max_yrs_dt_chi2_smallest_limit*secyer) then
               s% dt = max_yrs_dt_chi2_smallest_limit*secyer
               s% timestep_hold = s% model_number + 10
               write(*,'(a50,i6,2f16.2,1p,99e16.4)') 'redo timestep for "smallest" chi^2 limit', &
                  s% model_number, chi2, chi2_limit_for_smallest_timesteps, &
                  max_yrs_dt_chi2_smallest_limit
               do_simplex_extras_check_model = redo
               return
            end if         
         else if (chi2 <= chi2_limit_for_smaller_timesteps) then
            s% max_years_for_timestep = max_yrs_dt_chi2_smaller_limit 
            if (s% dt > max_yrs_dt_chi2_smaller_limit*secyer) then
               s% dt = max_yrs_dt_chi2_smaller_limit*secyer
               s% timestep_hold = s% model_number + 10
               write(*,'(a50,i6,2f16.2,1p,99e16.4)') 'redo timestep for "smaller" chi^2 limit', &
                  s% model_number, chi2, chi2_limit_for_smaller_timesteps, &
                  max_yrs_dt_chi2_smaller_limit
               do_simplex_extras_check_model = redo
               return
            end if         
         else if (chi2 <= chi2_limit_for_small_timesteps) then
            s% max_years_for_timestep = max_yrs_dt_chi2_small_limit 
            if (s% dt > max_yrs_dt_chi2_small_limit*secyer) then
               s% dt = max_yrs_dt_chi2_small_limit*secyer
               s% timestep_hold = s% model_number + 10
               write(*,'(a50,i6,2f16.2,1p,99e16.4)') 'redo timestep for "small" chi^2 limit', &
                  s% model_number, chi2, chi2_limit_for_small_timesteps, &
                  max_yrs_dt_chi2_small_limit
               do_simplex_extras_check_model = redo
               return
            end if         
         end if
         
         if (best_chi2 <= 0 .or. chi2 < best_chi2) then
            call save_best_info(s)
         end if
            
         call final_checks

         
         contains
         
         
         subroutine check_too_many_bad
            if (best_chi2 > 0) then
               num_chi2_too_big = num_chi2_too_big + 1
               if (num_chi2_too_big > limit_num_chi2_too_big) then
                  write(*,*) 'have reached too many bad chi2 limit'
                  do_simplex_extras_check_model = terminate
               end if
               return
            end if
            num_chi2_too_big = 0
         end subroutine check_too_many_bad
         
         
         subroutine final_checks
            if (include_age_in_chi2 .and. s% star_age >= max_age_for_chi2) then
               write(*,*) 'have reached max_age_for_chi2'
               do_simplex_extras_check_model = terminate
            end if
            if (eval_chi2_at_target_age_only .and. s% star_age >= age_target) then
               write(*,*) 'have reached age_target'
               do_simplex_extras_check_model = terminate
            end if
            if (best_chi2 > 0) then
               if (best_chi2 <= chi2_search_limit1 .and. &
                  chi2 >= chi2_search_limit2) then
                  write(*,*) 'have reached chi2_search_limit2'
                  do_simplex_extras_check_model = terminate
                  return
               end if
               if (chi2 >= chi2_relative_increase_limit*best_chi2) then
                  num_chi2_too_big = num_chi2_too_big + 1
                  if (num_chi2_too_big > limit_num_chi2_too_big) then
                     write(*,*) 'have reached too many bad chi2 limit'
                     do_simplex_extras_check_model = terminate
                  end if
                  return
               end if
               num_chi2_too_big = 0
            end if
         end subroutine final_checks
         
         
         subroutine check_limits
            real(dp) :: logg_limit, logL_limit, Teff_limit, &
               logR_limit, surface_Z_div_X_limit, surface_He_limit, solar_Rcz_limit, &
               solar_cs_rms_limit, my_var1_limit, my_var2_limit, my_var3_limit
            integer :: nz
            include 'formats'
            nz = s% nz
            
            if (s% star_age >= max_age_for_chi2) then
               write(*,*) 'have reached max_age_for_chi2'
               do_simplex_extras_check_model = terminate
               return
            end if

            if (sigmas_coeff_for_Teff_limit /= 0 .and. Teff_sigma > 0) then
               Teff_limit = Teff_target + Teff_sigma*sigmas_coeff_for_Teff_limit
               if ((sigmas_coeff_for_Teff_limit > 0 .and. s% Teff > Teff_limit) .or. &
                   (sigmas_coeff_for_Teff_limit < 0 .and. s% Teff < Teff_limit)) then
                  write(*,*) 'have reached Teff limit'
                  write(*,1) 'Teff', s% Teff
                  write(*,1) 'Teff_limit', Teff_limit
                  write(*,*)
                  do_simplex_extras_check_model = terminate
                  return
               end if    
               if (trace_limits) then
                  write(*,1) 'Teff', s% Teff
                  write(*,1) 'Teff_limit', Teff_limit
               end if
            end if
            
            if (sigmas_coeff_for_logg_limit /= 0 .and. logg_sigma > 0) then    
               logg_limit = logg_target + logg_sigma*sigmas_coeff_for_logg_limit
               if ((sigmas_coeff_for_logg_limit > 0 .and. logg > logg_limit) .or. &
                   (sigmas_coeff_for_logg_limit < 0 .and. logg < logg_limit)) then
                  write(*,*) 'have reached logg limit'
                  write(*,1) 'logg', logg
                  write(*,1) 'logg_limit', logg_limit
                  write(*,*)
                  do_simplex_extras_check_model = terminate
                  return
               end if
               if (trace_limits) then
                  write(*,1) 'logg', logg
                  write(*,1) 'logg_limit', logg_limit
               end if
            end if
            
            if (sigmas_coeff_for_logL_limit /= 0 .and. logL_sigma > 0) then
               logL_limit = logL_target + logL_sigma*sigmas_coeff_for_logL_limit
               if ((sigmas_coeff_for_logL_limit > 0 .and. s% log_surface_luminosity > logL_limit) .or. &
                   (sigmas_coeff_for_logL_limit < 0 .and. s% log_surface_luminosity < logL_limit)) then
                  write(*,*) 'have reached logL limit'
                  write(*,1) 'logL', s% log_surface_luminosity
                  write(*,1) 'logL_limit', logL_limit
                  write(*,*)
                  do_simplex_extras_check_model = terminate
                  return
               end if
               if (trace_limits) then
                  write(*,1) 'logL', s% log_surface_luminosity
                  write(*,1) 'logL_limit', logL_limit
               end if
            end if
            
            if (sigmas_coeff_for_logR_limit /= 0 .and. logR_sigma > 0) then
               logR_limit = logR_target + logR_sigma*sigmas_coeff_for_logR_limit
               if ((sigmas_coeff_for_logR_limit > 0 .and. logR > logR_limit) .or. &
                   (sigmas_coeff_for_logR_limit < 0 .and. logR < logR_limit)) then
                  write(*,*) 'have reached logR limit'
                  write(*,1) 'logR', logR
                  write(*,1) 'logR_limit', logR_limit
                  write(*,*)
                  do_simplex_extras_check_model = terminate
                  return
               end if
               if (trace_limits) then
                  write(*,1) 'logR', logR
                  write(*,1) 'logR_limit', logR_limit
               end if
            end if
            
            if (sigmas_coeff_for_surface_Z_div_X_limit /= 0 .and. surface_Z_div_X_sigma > 0) then
               surface_Z_div_X_limit = surface_Z_div_X_target + &
                  surface_Z_div_X_sigma*sigmas_coeff_for_surface_Z_div_X_limit
               if ((sigmas_coeff_for_surface_Z_div_X_limit > 0 .and. &
                     surface_Z_div_X > surface_Z_div_X_limit) .or. &
                     (sigmas_coeff_for_surface_Z_div_X_limit < 0 .and. &
                      surface_Z_div_X < surface_Z_div_X_limit)) then
                  write(*,*) 'have reached surface_Z_div_X limit'
                  write(*,1) 'surface_Z_div_X', surface_Z_div_X
                  write(*,1) 'surface_Z_div_X_limit', surface_Z_div_X_limit
                  write(*,*)
                  do_simplex_extras_check_model = terminate
                  return
               end if
               if (trace_limits) then
                  write(*,1) 'surface_Z_div_X', surface_Z_div_X
                  write(*,1) 'surface_Z_div_X_limit', surface_Z_div_X_limit
               end if
            end if
            
            if (sigmas_coeff_for_surface_He_limit /= 0 .and. surface_He_sigma > 0) then
               surface_He_limit = surface_He_target + &
                     surface_He_sigma*sigmas_coeff_for_surface_He_limit
               if ((sigmas_coeff_for_surface_He_limit > 0 .and. &
                     surface_He > surface_He_limit) .or. &
                   (sigmas_coeff_for_surface_He_limit < 0 .and. &
                     surface_He < surface_He_limit)) then
                  write(*,*) 'have reached surface_He limit'
                  write(*,1) 'surface_He', surface_He
                  write(*,1) 'surface_He_limit', surface_He_limit
                  write(*,*)
                  do_simplex_extras_check_model = terminate
                  return
               end if
               if (trace_limits) then
                  write(*,1) 'surface_He', surface_He
                  write(*,1) 'surface_He_limit', surface_He_limit
               end if
            end if
            
            if (sigmas_coeff_for_solar_Rcz_limit /= 0 .and. Rcz_sigma > 0) then
               solar_Rcz_limit = Rcz_target + Rcz_sigma*sigmas_coeff_for_solar_Rcz_limit
               if ((sigmas_coeff_for_solar_Rcz_limit > 0 .and. Rcz > solar_Rcz_limit) .or. &
                   (sigmas_coeff_for_solar_Rcz_limit < 0 .and. Rcz < solar_Rcz_limit)) then
                  write(*,*) 'have reached Rcz limit'
                  write(*,1) 'Rcz', Rcz
                  write(*,1) 'solar_Rcz_limit', solar_Rcz_limit
                  write(*,*)
                  do_simplex_extras_check_model = terminate
                  return
               end if
               if (trace_limits) then
                  write(*,1) 'Rcz', Rcz
                  write(*,1) 'solar_Rcz_limit', solar_Rcz_limit
               end if
            end if
            
            if (sigmas_coeff_for_solar_cs_rms_limit /= 0 .and. solar_cs_rms_sigma > 0) then
               solar_cs_rms_limit = solar_cs_rms_target + &
                     solar_cs_rms_sigma*sigmas_coeff_for_solar_cs_rms_limit
               if ((sigmas_coeff_for_solar_cs_rms_limit > 0 .and. &
                     solar_cs_rms > solar_cs_rms_limit) .or. &
                   (sigmas_coeff_for_solar_cs_rms_limit < 0 .and. &
                     solar_cs_rms < solar_cs_rms_limit)) then
                  write(*,*) 'have reached solar_cs_rms limit'
                  write(*,1) 'solar_cs_rms', solar_cs_rms
                  write(*,1) 'solar_cs_rms_limit', solar_cs_rms_limit
                  write(*,*)
                  do_simplex_extras_check_model = terminate
                  return
               end if
               if (trace_limits) then
                  write(*,1) 'solar_cs_rms', solar_cs_rms
                  write(*,1) 'solar_cs_rms_limit', solar_cs_rms_limit
               end if
            end if
            
            if (sigmas_coeff_for_my_var1_limit /= 0 .and. my_var1_sigma > 0) then
               my_var1_limit = &
                  my_var1_target + my_var1_sigma*sigmas_coeff_for_my_var1_limit
               if ((sigmas_coeff_for_my_var1_limit > 0 .and. &
                        my_var1 > my_var1_limit) .or. &
                   (sigmas_coeff_for_my_var1_limit < 0 .and. &
                        my_var1 < my_var1_limit)) then
                  write(*,*) 'have reached my_var1 limit'
                  write(*,1) 'my_var1', my_var1
                  write(*,1) 'my_var1_limit', my_var1_limit
                  write(*,*)
                  do_simplex_extras_check_model = terminate
                  return
               end if
               if (trace_limits) then
                  write(*,1) 'my_var1', my_var1
                  write(*,1) 'my_var1_limit', my_var1_limit
               end if
            end if
            
            if (sigmas_coeff_for_my_var2_limit /= 0 .and. my_var2_sigma > 0) then
               my_var2_limit = &
                  my_var2_target + my_var2_sigma*sigmas_coeff_for_my_var2_limit
               if ((sigmas_coeff_for_my_var2_limit > 0 .and. &
                        my_var2 > my_var2_limit) .or. &
                   (sigmas_coeff_for_my_var2_limit < 0 .and. &
                        my_var2 < my_var2_limit)) then
                  write(*,*) 'have reached my_var2 limit'
                  write(*,1) 'my_var2', my_var2
                  write(*,1) 'my_var2_limit', my_var2_limit
                  write(*,*)
                  do_simplex_extras_check_model = terminate
                  return
               end if
               if (trace_limits) then
                  write(*,1) 'my_var2', my_var2
                  write(*,1) 'my_var2_limit', my_var2_limit
               end if
            end if
            
            if (sigmas_coeff_for_my_var3_limit /= 0 .and. my_var3_sigma > 0) then
               my_var3_limit = &
                  my_var3_target + my_var3_sigma*sigmas_coeff_for_my_var3_limit
               if ((sigmas_coeff_for_my_var3_limit > 0 .and. &
                        my_var3 > my_var3_limit) .or. &
                   (sigmas_coeff_for_my_var3_limit < 0 .and. &
                        my_var3 < my_var3_limit)) then
                  write(*,*) 'have reached my_var3 limit'
                  write(*,1) 'my_var3', my_var3
                  write(*,1) 'my_var3_limit', my_var3_limit
                  write(*,*)
                  do_simplex_extras_check_model = terminate
                  return
               end if
               if (trace_limits) then
                  write(*,1) 'my_var3', my_var3
                  write(*,1) 'my_var3_limit', my_var3_limit
               end if
            end if
            
         end subroutine check_limits


         subroutine setup_solar_data_for_calc_rms(ierr)
            use const_def, only: mesa_data_dir
            integer, intent(out) :: ierr
         
            integer, parameter :: lines_to_skip = 11
            integer :: iounit, i, k
            real(dp) :: jnk
         
            character (len=256) :: fname
         
            have_sound_speed_data = .true.
            ierr = 0
            iounit = alloc_iounit(ierr)
            if (ierr /= 0) then
               return
            end if
            fname = trim(mesa_data_dir) // '/star_data/solar_csound.data'
            open(iounit, file=trim(fname), action='read', status='old', iostat=ierr)
            if (ierr /= 0) then
               write(*,*) 'failed to open ' // trim(fname)
               call free_iounit(iounit)
               return
            end if                  
         
            do i=1,lines_to_skip
               read(iounit,fmt=*,iostat=ierr)
               if (ierr /= 0) then
                  write(*,*) 'failed when skipping line', i
                  close(iounit)
                  call free_iounit(iounit)
                  return
               end if
            end do
         
            do i=1,npts
               read(iounit,fmt=*,iostat=ierr) &
                  data_r(i), jnk, data_csound(i), jnk, jnk, jnk, data_width(i)
               if (ierr /= 0) then
                  write(*,*) 'failed when reading data point', i
                  close(iounit)
                  call free_iounit(iounit)
                  return
               end if
            end do

            close(iounit)
            call free_iounit(iounit)
      
         end subroutine setup_solar_data_for_calc_rms


         real(dp) function calc_current_rms(s, nz) ! dR weighted
            use interp_1d_lib
            use interp_1d_def
            type (star_info), pointer :: s
            integer, intent(in) :: nz
         
            logical, parameter :: dbg = .false.
            real(dp), target :: calc_rms_f1_ary(4*nz)
            real(dp), pointer :: calc_rms_f1(:), calc_rms_f(:,:)
            real(dp) :: sumy2, sumdr, dr, y2, cs
            real(dp), parameter :: min_R = 0.094, max_R = 0.94
            real(dp), target :: pm_work_ary(nz*pm_work_size)
            real(dp), pointer :: pm_work(:)
            integer :: k, i, ierr
         
            include 'formats'
         
            calc_current_rms = -1  
            pm_work => pm_work_ary
            calc_rms_f1 => calc_rms_f1_ary
            calc_rms_f(1:4,1:nz) => calc_rms_f1(1:4*nz)

            ierr = 0
            if (.not. have_sound_speed_data) then
               call setup_solar_data_for_calc_rms(ierr)
               if (ierr /= 0) then
                  write(*,*) 'failed in setup_solar_data_for_calc_rms'
                  return
               end if
            end if
         
            do k=1,nz
               if (k == 1) then
                  calc_rms_f(1,k) = s% csound(k)
               else
                  calc_rms_f(1,k) = &
                     (s% csound(k)*s% dq(k-1) + s% csound(k-1)*s% dq(k))/(s% dq(k-1) + s% dq(k))
               end if
            end do

            call interp_pm( &
               s% r, nz, calc_rms_f1, pm_work_size, pm_work, 'calc_current_rms', ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in interp_pm'
               return
            end if

            sumy2 = 0
            sumdr = 0
            do i=1,npts
               if (data_r(i) < min_R .or. data_r(i) > max_R) cycle
               call interp_value(s% r, nz, calc_rms_f1, Rsun*data_r(i), cs, ierr)
               if (ierr /= 0) then
                  write(*,*) 'failed in interp_value', i
                  return
               end if
               if (i == 1) then
                  dr = data_r(2) - data_r(1)
               else
                  dr = data_r(i) - data_r(i-1)
               end if
               if (dr < 0) dr = -dr
               ! change to weigh by point rather than by dr
               dr = 1
            
               sumdr = sumdr + dr
               y2 = dr*pow2((cs - data_csound(i))/data_csound(i))
               sumy2 = sumy2 + y2
               if (dbg) write(*,2) 'rms cs, data_cs, reldiff, y2, dr', i, cs, data_csound(i), &
                  (cs - data_csound(i))/data_csound(i), y2, dr
            end do
         
            calc_current_rms = sqrt(sumy2/sumdr)
            if (dbg) write(*,1) 'calc_current_rms', calc_current_rms

         end function calc_current_rms
         
      end function do_simplex_extras_check_model


      real(dp) function get_chi2(s, trace_okay, ierr)
         type (star_info), pointer :: s
         logical, intent(in) :: trace_okay
         integer, intent(out) :: ierr

         integer :: i, n, chi2N
         real(dp) :: chi2term, Teff, logL, chi2sum
         
         ! calculate chi^2 following Brandao et al, 2011, eqn 11
         include 'formats'
         
         ierr = 0
         chi2sum = 0
         chi2N = 0
         
         call star_simplex_procs% set_my_vars(s% id, ierr)
         if (ierr /= 0) return
         
         if (Teff_sigma > 0 .and. include_Teff_in_chi2) then
            Teff = s% Teff
            chi2term = pow2((Teff - Teff_target)/Teff_sigma)
            if (trace_okay .and. trace_chi2_info) &
               write(*,2) 'chi2_term Teff', s% model_number, chi2term
            chi2sum = chi2sum + chi2term
            chi2N = chi2N + 1
         end if
         
         if (logL_sigma > 0 .and. include_logL_in_chi2) then
            logL = s% log_surface_luminosity
            chi2term = pow2((logL - logL_target)/logL_sigma)
            if (trace_okay .and. trace_chi2_info) &
               write(*,2) 'chi2_term logL', s% model_number, chi2term
            chi2sum = chi2sum + chi2term
            chi2N = chi2N + 1
         end if
         
         if (logg_sigma > 0 .and. include_logg_in_chi2) then
            chi2term = pow2((logg - logg_target)/logg_sigma)
            if (trace_okay .and. trace_chi2_info) &
               write(*,2) 'chi2_term logg', s% model_number, chi2term
            chi2sum = chi2sum + chi2term
            chi2N = chi2N + 1
         end if
         
         if (FeH_sigma > 0 .and. include_FeH_in_chi2) then
            chi2term = pow2((FeH - FeH_target)/FeH_sigma)
            if (trace_okay .and. trace_chi2_info) &
               write(*,2) 'chi2_term FeH', s% model_number, chi2term
            chi2sum = chi2sum + chi2term
            chi2N = chi2N + 1
         end if
         
         if (logR_sigma > 0 .and. include_logR_in_chi2) then
            chi2term = pow2((logR - logR_target)/logR_sigma)
            if (trace_okay .and. trace_chi2_info) &
               write(*,2) 'chi2_term logR', s% model_number, chi2term
            chi2sum = chi2sum + chi2term
            chi2N = chi2N + 1
         end if
         
         if (age_sigma > 0 .and. include_age_in_chi2) then
            chi2term = pow2((s% star_age - age_target)/age_sigma)
            if (trace_okay .and. trace_chi2_info) &
               write(*,2) 'chi2_term age', s% model_number, chi2term
            chi2sum = chi2sum + chi2term
            chi2N = chi2N + 1
         end if
         
         if (surface_Z_div_X_sigma > 0 .and. include_surface_Z_div_X_in_chi2) then
            chi2term = pow2((surface_Z_div_X - surface_Z_div_X_target)/surface_Z_div_X_sigma)
            if (trace_okay .and. trace_chi2_info) &
               write(*,2) 'chi2_term surface_Z_div_X', s% model_number, chi2term
            chi2sum = chi2sum + chi2term
            chi2N = chi2N + 1
         end if
         
         if (surface_He_sigma > 0 .and. include_surface_He_in_chi2) then
            chi2term = pow2((surface_He - surface_He_target)/surface_He_sigma)
            if (trace_okay .and. trace_chi2_info) &
               write(*,2) 'chi2_term surface_He', s% model_number, chi2term
            chi2sum = chi2sum + chi2term
            chi2N = chi2N + 1
         end if
         
         if (Rcz_sigma > 0 .and. include_Rcz_in_chi2) then
            chi2term = pow2((Rcz - Rcz_target)/Rcz_sigma)
            if (trace_okay .and. trace_chi2_info) &
               write(*,2) 'chi2_term Rcz', s% model_number, chi2term
            chi2sum = chi2sum + chi2term
            chi2N = chi2N + 1
         end if
         
         if (solar_cs_rms_sigma > 0 .and. include_solar_cs_rms_in_chi2) then
            chi2term = pow2((solar_cs_rms - solar_cs_rms_target)/solar_cs_rms_sigma)
            if (trace_okay .and. trace_chi2_info) &
               write(*,2) 'chi2_term solar_cs_rms', s% model_number, chi2term
            chi2sum = chi2sum + chi2term
            chi2N = chi2N + 1
         end if
         
         if (my_var1_sigma > 0 .and. include_my_var1_in_chi2) then
            chi2term = pow2((my_var1 - my_var1_target)/my_var1_sigma)
            if (trace_okay .and. trace_chi2_info) &
               write(*,2) 'chi2_term ' // trim(my_var1_name), s% model_number, chi2term
            chi2sum = chi2sum + chi2term
            chi2N = chi2N + 1
         end if
         
         if (my_var2_sigma > 0 .and. include_my_var2_in_chi2) then
            chi2term = pow2((my_var2 - my_var2_target)/my_var2_sigma)
            if (trace_okay .and. trace_chi2_info) &
               write(*,2) 'chi2_term ' // trim(my_var2_name), s% model_number, chi2term
            chi2sum = chi2sum + chi2term
            chi2N = chi2N + 1
         end if
         
         if (my_var3_sigma > 0 .and. include_my_var3_in_chi2) then
            chi2term = pow2((my_var3 - my_var3_target)/my_var3_sigma)
            if (trace_okay .and. trace_chi2_info) &
               write(*,2) 'chi2_term ' // trim(my_var3_name), s% model_number, chi2term
            chi2sum = chi2sum + chi2term
            chi2N = chi2N + 1
         end if

         num_chi2_terms = chi2N
         chi2 = chi2sum/max(1,chi2N)

         get_chi2 = chi2
                           
      end function get_chi2
      
      
      subroutine save_best_info(s)     
         type (star_info), pointer :: s
         integer :: ierr
         logical :: write_controls_info_with_profile
         
         include 'formats'
         
         if (save_model_for_best_model) then
            ierr = 0
            call star_write_model(s% id, best_model_save_model_filename, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in star_write_model'
               stop 1
            end if
            write(*, '(a,i7)') 'save ' // trim(best_model_save_model_filename), s% model_number
         end if

         if (write_profile_for_best_model) then
            ierr = 0
            write_controls_info_with_profile = s% write_controls_info_with_profile
            s% write_controls_info_with_profile = .false.
            call star_write_profile_info(s% id, best_model_profile_filename, ierr)
            s% write_controls_info_with_profile = write_controls_info_with_profile
            if (ierr /= 0) then
               write(*,*) 'failed in star_write_profile_info'
               stop 1
            end if
            call save_profile(s% id, 3, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in save_profile'
               stop 1
            end if
         end if         
         
         call store_best_info(s)
            
      end subroutine save_best_info
            
   
      subroutine init_sample_ptrs
         nullify( &
            sample_chi2, &
            sample_age, &
            sample_init_Y, &
            sample_init_FeH, &
            sample_init_h1, &
            sample_init_he3, &
            sample_init_he4, &
            sample_init_Z, &
            sample_mass, &
            sample_alpha, &
            sample_f_ov, &
            sample_radius, &
            sample_logL, &
            sample_Teff, &
            sample_logg, &
            sample_FeH, &
            sample_logR, &
            sample_surface_Z_div_X, &
            sample_surface_He, &
            sample_Rcz, &
            sample_solar_cs_rms, &
            sample_my_var1, &
            sample_my_var2, &
            sample_my_var3, &
            sample_model_number, &
            sample_index_by_chi2)
      end subroutine init_sample_ptrs
      
      
      subroutine alloc_sample_ptrs(ierr)
         use utils_lib
         integer, intent(out) :: ierr
         ierr = 0
         max_num_samples = 1.5*max_num_samples + 200
         
         call realloc_double(sample_chi2,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_age,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_init_Y,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_init_FeH,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_init_h1,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_init_he3,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_init_he4,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_init_Z,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_mass,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_alpha,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_f_ov,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_radius,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_logL,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_Teff,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_logg,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_FeH,max_num_samples,ierr); if (ierr /= 0) return
         
         call realloc_double(sample_logR,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_surface_Z_div_X,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_surface_He,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_Rcz,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_solar_cs_rms,max_num_samples,ierr); if (ierr /= 0) return
         
         call realloc_double(sample_my_var1,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_my_var2,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_my_var3,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_my_param1,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_my_param2,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_double(sample_my_param3,max_num_samples,ierr); if (ierr /= 0) return
         
         call realloc_integer(sample_index_by_chi2,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_integer(sample_op_code,max_num_samples,ierr); if (ierr /= 0) return
         call realloc_integer(sample_model_number,max_num_samples,ierr); if (ierr /= 0) return

      end subroutine alloc_sample_ptrs
   
   
      subroutine read_simplex_search_controls(filename, ierr)
         character (len=*), intent(in) :: filename
         integer, intent(out) :: ierr
         ! initialize controls to default values
         call set_simplex_search_defaults
         ierr = 0
         call read1_simplex_search_inlist(filename, 1, ierr)
      end subroutine read_simplex_search_controls
      
      
      subroutine set_simplex_search_defaults
            include 'simplex_solar.defaults'
      end subroutine set_simplex_search_defaults

         
      recursive subroutine read1_simplex_search_inlist(filename, level, ierr)
         character (len=*), intent(in) :: filename
         integer, intent(in) :: level  
         integer, intent(out) :: ierr
         
         logical :: read_extra1, read_extra2, read_extra3, read_extra4, read_extra5
         character (len=256) :: message, extra1, extra2, extra3, extra4, extra5
         integer :: unit
         
         if (level >= 10) then
            write(*,*) 'ERROR: too many levels of nested extra star_job inlist files'
            ierr = -1
            return
         end if
         
         ierr = 0
         unit=alloc_iounit(ierr)
         if (ierr /= 0) return
         
         open(unit=unit, file=trim(filename), action='read', delim='quote', iostat=ierr)
         if (ierr /= 0) then
            write(*, *) 'Failed to open simplex search inlist file ', trim(filename)
         else
            read(unit, nml=simplex_search_controls, iostat=ierr)  
            close(unit)
            if (ierr /= 0) then
               write(*, *) &
                  'Failed while trying to read simplex search inlist file ', trim(filename)
               write(*, '(a)') trim(message)
               write(*, '(a)') &
                  'The following runtime error message might help you find the problem'
               write(*, *) 
               open(unit=unit, file=trim(filename), &
                  action='read', delim='quote', status='old', iostat=ierr)
               read(unit, nml=simplex_search_controls)
               close(unit)
            end if  
         end if
         call free_iounit(unit)
         if (ierr /= 0) return
         
         ! recursive calls to read other inlists
         
         read_extra1 = read_extra_simplex_search_inlist1
         read_extra_simplex_search_inlist1 = .false.
         extra1 = extra_simplex_search_inlist1_name
         extra_simplex_search_inlist1_name = 'undefined'
         
         read_extra2 = read_extra_simplex_search_inlist2
         read_extra_simplex_search_inlist2 = .false.
         extra2 = extra_simplex_search_inlist2_name
         extra_simplex_search_inlist2_name = 'undefined'
         
         read_extra3 = read_extra_simplex_search_inlist3
         read_extra_simplex_search_inlist3 = .false.
         extra3 = extra_simplex_search_inlist3_name
         extra_simplex_search_inlist3_name = 'undefined'
         
         read_extra4 = read_extra_simplex_search_inlist4
         read_extra_simplex_search_inlist4 = .false.
         extra4 = extra_simplex_search_inlist4_name
         extra_simplex_search_inlist4_name = 'undefined'
         
         read_extra5 = read_extra_simplex_search_inlist5
         read_extra_simplex_search_inlist5 = .false.
         extra5 = extra_simplex_search_inlist5_name
         extra_simplex_search_inlist5_name = 'undefined'
         
         if (read_extra1) then
            !write(*,*) 'read extra simplex_search inlist1 from ' // trim(extra1)
            call read1_simplex_search_inlist(extra1, level+1, ierr)
            if (ierr /= 0) return
         end if
         
         if (read_extra2) then
            !write(*,*) 'read extra simplex_search inlist2 from ' // trim(extra2)
            call read1_simplex_search_inlist(extra2, level+1, ierr)
            if (ierr /= 0) return
         end if
         
         if (read_extra3) then
            !write(*,*) 'read extra simplex_search inlist3 from ' // trim(extra3)
            call read1_simplex_search_inlist(extra3, level+1, ierr)
            if (ierr /= 0) return
         end if
         
         if (read_extra4) then
            !write(*,*) 'read extra simplex_search inlist4 from ' // trim(extra4)
            call read1_simplex_search_inlist(extra4, level+1, ierr)
            if (ierr /= 0) return
         end if
         
         if (read_extra5) then
            write(*,*) 'read extra simplex_search inlist5 from ' // trim(extra5)
            call read1_simplex_search_inlist(extra5, level+1, ierr)
            if (ierr /= 0) return
         end if
         
      end subroutine read1_simplex_search_inlist


      subroutine write_simplex_search_controls(filename_in, ierr)
         use utils_lib
         character(*), intent(in) :: filename_in
         integer, intent(out) :: ierr
         character (len=256) :: filename
         integer :: unit
         ierr = 0
         filename = trim(filename_in)
         if (len_trim(filename) == 0) filename = 'simplex_search_controls.out'
         unit=alloc_iounit(ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to alloc iounit in write_simplex_search_controls'
            return
         end if
         ! NOTE: when open namelist file, must include delim='APOSTROPHE'
         open(unit=unit, file=trim(filename), action='write', delim='APOSTROPHE', iostat=ierr)
         if (ierr /= 0) then
            write(*, *) 'Failed to open ', trim(filename)
         else
            write(unit, nml=simplex_search_controls)
            close(unit)
         end if
         call free_iounit(unit)
         
         write(*,*)
         write(*,*) 'saved initial &simplex_search_controls to ' // trim(filename)
         write(*,*)
         write(*,*)

      end subroutine write_simplex_search_controls


      subroutine save_sample_results_to_file(i_total, results_fname, ierr)
         use utils_lib
         integer, intent(in) :: i_total
         character (len=*), intent(in) :: results_fname
         integer, intent(out) :: ierr
         integer :: iounit
         write(*,*) 'save_sample_results_to_file ' // trim(results_fname)
         iounit = alloc_iounit(ierr)
         if (ierr /= 0) stop 'alloc_iounit failed'
         open(unit=iounit, file=trim(results_fname), action='write', iostat=ierr)
         if (ierr /= 0) return
         call show_all_sample_results(iounit, i_total, ierr)
         close(iounit)
         call free_iounit(iounit)         
      end subroutine save_sample_results_to_file
      
      
      subroutine set_sample_index_by_chi2
         use num_lib, only: qsort
         if (sample_number <= 0) return
         if (sample_number == 1) then
            sample_index_by_chi2(1) = 1
            return
         end if
         call qsort(sample_index_by_chi2, sample_number, sample_chi2)
      end subroutine set_sample_index_by_chi2
      
      
      subroutine show_sample_header(iounit)
         integer, intent(in) ::iounit
         
         integer :: j
         character (len=10) :: str
      
         write(iounit,'(2x,a6,7a26,a16,99a26)') &
            'sample', &
            
            'chi2', &
            'mass', &
            'init_Y', &
            'init_FeH', &
            'alpha', &
            'f_ov', &
            'age', &
            
            'model_number', &
            
            'init_h1', &
            'init_he3', &
            'init_he4', &
            'init_Z', &
            'log_radius', &
            'logL', &
            'Teff', &
            'logg', &
            'Fe_H', &
            'logR', &
            'surface_Z_div_X', &
            'surface_He', &
            'Rcz', &
            'solar_cs_rms', &
            trim(my_var1_name), &
            trim(my_var2_name), &
            trim(my_var3_name), &
            trim(my_param1_name), &
            trim(my_param2_name), &
            trim(my_param3_name)
                                                   
      end subroutine show_sample_header
      
      
      subroutine show1_sample_results(i, iounit)
         use num_lib, only: simplex_info_str
         integer, intent(in) :: i, iounit
            
         integer :: j, k, op_code, ierr
         character (len=256) :: info_str
         
         ierr = 0

         op_code = sample_op_code(i) 
         if (op_code <= 0) then
            info_str = ''
         else
            call simplex_info_str(op_code, info_str, ierr)
            if (ierr /= 0) then
               info_str = ''
               ierr = 0
            end if
         end if
         
         write(iounit,'(3x,i5,7(1pd26.16),i16,99(1pd26.16))',advance='no') i, &
            sample_chi2(i), &
            sample_mass(i), &
            sample_init_Y(i), &
            sample_init_FeH(i), &
            sample_alpha(i), &
            sample_f_ov(i), &
            sample_age(i), &
            sample_model_number(i), &
            sample_init_h1(i), &
            sample_init_he3(i), &
            sample_init_he4(i), &
            sample_init_Z(i), &
            safe_log10(sample_radius(i)), &
            sample_logL(i), &
            sample_Teff(i), &
            sample_logg(i), &
            sample_FeH(i), &
            sample_logR(i), &
            sample_surface_Z_div_X(i), &
            sample_surface_He(i), &
            sample_Rcz(i), &
            sample_solar_cs_rms(i), &
            sample_my_var1(i), &
            sample_my_var2(i), &
            sample_my_var3(i), &
            sample_my_param1(i), &
            sample_my_param2(i), &
            sample_my_param3(i)
            
         if (iounit == 6) return

         write(iounit,'(a12)') trim(info_str)
      
      
      end subroutine show1_sample_results
      
      
      subroutine show_all_sample_results(iounit, i_total, ierr)
         integer, intent(in) :: iounit, i_total
         integer, intent(out) :: ierr
         integer :: i, j

         ierr = 0
         ! sort results by increasing sample_chi2
         call set_sample_index_by_chi2
         if (i_total > 0) then
            write(iounit,*) sample_number, ' of ', i_total
         else
            write(iounit,*) sample_number, ' samples'
         end if
         call show_sample_header(iounit)
         do j = 1, sample_number
            i = sample_index_by_chi2(j)
            call show1_sample_results(i, iounit)
         end do

         call show_sample_header(iounit)
         do i = 1, 3
            write(iounit,*)
         end do

      end subroutine show_all_sample_results

      
      subroutine show_best(io)
         integer, intent(in) :: io
         
         real(dp) :: chi2term
         include 'formats'

         if (Teff_sigma > 0 .and. include_Teff_in_chi2) then
            chi2term = pow2((best_Teff - Teff_target)/Teff_sigma)
            write(io,*)
            call write1('Teff chi2term', chi2term)
            call write1('Teff', best_Teff)
            call write1('Teff_obs', Teff_target)
            call write1('Teff_sigma', Teff_sigma)
         end if
         
         if (logL_sigma > 0 .and. include_logL_in_chi2) then
            chi2term = pow2((best_logL - logL_target)/logL_sigma)
            write(io,*)
            call write1('logL chi2term', chi2term)
            call write1('logL', best_logL)
            call write1('logL_obs', logL_target)
            call write1('logL_sigma', logL_sigma)
         end if
         
         if (logg_sigma > 0 .and. include_logg_in_chi2) then
            chi2term = pow2((best_logg - logg_target)/logg_sigma)
            write(io,*)
            call write1('logg chi2term', chi2term)
            call write1('logg', best_logg)
            call write1('logg_obs', logg_target)
            call write1('logg_sigma', logg_sigma)
         end if
         
         if (FeH_sigma > 0 .and. include_FeH_in_chi2) then
            chi2term = pow2((best_FeH - FeH_target)/FeH_sigma)
            write(io,*)
            call write1('FeH chi2term', chi2term)
            call write1('FeH', best_FeH)
            call write1('FeH_obs', FeH_target)
            call write1('FeH_sigma', FeH_sigma)
         end if
         
         if (logR_sigma > 0 .and. include_logR_in_chi2) then
            chi2term = pow2((best_logR - logR_target)/logR_sigma)
            write(io,*)
            call write1('logR chi2term', chi2term)
            call write1('logR', best_logR)
            call write1('logR_obs', logR_target)
            call write1('logR_sigma', logR_sigma)
         end if
         
         if (age_sigma > 0 .and. include_age_in_chi2) then
            chi2term = pow2((best_age - age_target)/age_sigma)
            write(io,*)
            write(io,'(a40,e20.10,99f20.10)') 'age chi2term', chi2term
            write(io,'(a40,1pd20.10)') 'age', best_age
            write(io,'(a40,1pd20.10)') 'age_target', age_target
            write(io,'(a40,1pd20.10)') 'age_sigma', age_sigma
         end if
         
         if (surface_Z_div_X_sigma > 0 .and. &
               include_surface_Z_div_X_in_chi2) then
            chi2term = &
               pow2((best_surface_Z_div_X - surface_Z_div_X_target)/surface_Z_div_X_sigma)
            write(io,*)
            write(io,'(a40,e20.10,99f20.10)') 'surface_Z_div_X chi2term', chi2term
            call write1('surface_Z_div_X', best_surface_Z_div_X)
            call write1('surface_Z_div_X_obs', surface_Z_div_X_target)
            call write1('surface_Z_div_X_sigma', surface_Z_div_X_sigma)
         end if
         
         if (surface_He_sigma > 0 .and. include_surface_He_in_chi2) then
            chi2term = pow2((best_surface_He - surface_He_target)/surface_He_sigma)
            write(io,*)
            call write1('surface_He chi2term', chi2term)
            call write1('surface_He', best_surface_He)
            call write1('surface_He_obs', surface_He_target)
            call write1('surface_He_sigma', surface_He_sigma)
         end if
         
         if (Rcz_sigma > 0 .and. include_Rcz_in_chi2) then
            chi2term = pow2((best_Rcz - Rcz_target)/Rcz_sigma)
            write(io,*)
            call write1('Rcz chi2term', chi2term)
            call write1('Rcz', best_Rcz)
            call write1('Rcz_obs', Rcz_target)
            call write1('Rcz_sigma', Rcz_sigma)
         end if
         
         if (solar_cs_rms_sigma > 0 .and. include_solar_cs_rms_in_chi2) then
            chi2term = pow2((best_solar_cs_rms - solar_cs_rms_target)/solar_cs_rms_sigma)
            write(io,*)
            call write1('solar_cs_rms chi2term', chi2term)
            call write1('solar_cs_rms', best_solar_cs_rms)
            call write1('solar_cs_rms_obs', solar_cs_rms_target)
            call write1('solar_cs_rms_sigma', solar_cs_rms_sigma)
         end if
         
         if (my_var1_sigma > 0 .and. include_my_var1_in_chi2) then
            chi2term = pow2( &
                  (best_my_var1 - my_var1_target)/my_var1_sigma)
            write(io,*)
            call write1(trim(my_var1_name) // ' chi2term', chi2term)
            call write1(trim(my_var1_name), best_my_var1)
            call write1(trim(my_var1_name) // '_obs', my_var1_target)
            call write1(trim(my_var1_name) // '_sigma', my_var1_sigma)
         end if
         
         if (my_var2_sigma > 0 .and. include_my_var2_in_chi2) then
            chi2term = pow2( &
                  (best_my_var2 - my_var2_target)/my_var2_sigma)
            write(io,*)
            call write1(trim(my_var2_name) // ' chi2term', chi2term)
            call write1(trim(my_var2_name), best_my_var2)
            call write1(trim(my_var2_name) // '_obs', my_var2_target)
            call write1(trim(my_var2_name) // '_sigma', my_var2_sigma)
         end if
         
         if (my_var3_sigma > 0 .and. include_my_var3_in_chi2) then
            chi2term = pow2( &
                  (best_my_var3 - my_var3_target)/my_var3_sigma)
            write(io,*)
            call write1(trim(my_var3_name) // ' chi2term', chi2term)
            call write1(trim(my_var3_name), best_my_var3)
            call write1(trim(my_var3_name) // '_obs', my_var3_target)
            call write1(trim(my_var3_name) // '_sigma', my_var3_sigma)
         end if
         
         write(io,*)
         call write1('R/Rsun', best_radius)
         call write1('logL/Lsun', best_logL)
         call write1('Teff', best_Teff)
         call write1('logg', best_logg)
         call write1('FeH', best_FeH)
         call write1('logR', best_logR)
         call write1('surface_Z_div_X', best_surface_Z_div_X)
         call write1('surface_He', best_surface_He)
         call write1('Rcz', best_Rcz)
         call write1('solar_cs_rms', best_solar_cs_rms)
         write(io,*)        
         call write1('initial h1', current_h1)
         call write1('initial he3', current_he3)
         call write1('initial he4', current_he4)
         call write1('initial Y', current_Y)
         call write1('initial Z', current_Z)
         call write1('initial FeH', current_FeH)
         write(io,*)        
         call write1('mass/Msun', current_mass)
         call write1('alpha', current_alpha)
         call write1('f_ov', current_f_ov)
         write(io,'(a40,1pd26.16)') 'age', best_age
         write(io,*)
         call write1('chi^2', best_chi2)
         write(io,*)
         write(io,'(a40,i16)') 'model number', best_model_number
         write(io,*)
         write(io,*)
         
         contains
         
         subroutine write1(str,x)
            character (len=*), intent(in) :: str
            real(dp), intent(in) :: x
            if (abs(x) < 1d6) then
               write(io,'(a40,99f20.10)') trim(str), x
            else
               write(io,'(a40,99e20.10)') trim(str), x
            end if
         end subroutine write1

      end subroutine show_best

      
      subroutine store_best_info(s)
         type (star_info), pointer :: s
         integer :: i
      
         best_chi2 = chi2
         
         best_age = s% star_age
         best_model_number = s% model_number
         best_radius = s% photosphere_r
         best_logL = s% log_surface_luminosity
         best_Teff = s% Teff
         best_logg = logg
         best_FeH = FeH
         
         best_logR = logR
         best_surface_Z_div_X = surface_Z_div_X
         best_surface_He = surface_He
         best_Rcz = Rcz
         best_solar_cs_rms = solar_cs_rms
         
         best_my_var1 = my_var1
         best_my_var2 = my_var2
         best_my_var3 = my_var3
         best_my_param1 = my_param1
         best_my_param2 = my_param2
         best_my_param3 = my_param3
      
      end subroutine store_best_info
      
      
      subroutine write_best(num)
         use utils_lib, only: mkdir
         integer, intent(in) :: num
         integer :: ierr, iounit
         character (len=256) :: format_string, num_string, filename
         integer, parameter :: max_len_out = 2000
         character (len=max_len_out) :: script         
         ierr = 0
         iounit = alloc_iounit(ierr)
         if (ierr /= 0) return
         write(format_string,'( "(i",i2.2,".",i2.2,")" )') num_digits, num_digits
         write(num_string,format_string) num
         call mkdir(sample_results_dir)
         filename = trim(sample_results_dir) // '/' // &
            trim(sample_results_prefix) // trim(num_string) // trim(sample_results_postfix)
         open(unit=iounit, file=trim(filename), action='write', status='replace', iostat=ierr)
         if (ierr == 0) then
            call show_best(iounit)
            close(iounit)
            write(*,*) 'save best model results to ' // trim(filename)
         else
            write(*,*) 'failed to open ' // trim(filename)
         end if
         call free_iounit(iounit)
      end subroutine write_best


      subroutine read_samples_from_file(results_fname, ierr)
         use utils_lib
         character (len=*), intent(in) :: results_fname
         integer, intent(out) :: ierr
         integer :: iounit, num, i, j, model_number
         character (len=100) :: line
         
         include 'formats'
         
         ierr = 0         
         write(*,*) 'read samples from file ' // trim(results_fname)
         
         iounit = alloc_iounit(ierr)
         if (ierr /= 0) stop 'alloc_iounit failed'
         open(unit=iounit, file=trim(results_fname), action='read', status='old', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(results_fname)
            call free_iounit(iounit) 
            return
         end if
         
         read(iounit, fmt=*, iostat=ierr) num
         if (ierr /= 0) then
            write(*,*) 'failed to read number of samples on 1st line of ' // trim(results_fname)
            call done
            return
         end if
         
         write(*,2) 'number of samples in file', num
         
         read(iounit, fmt='(a)', iostat=ierr) line
         if (ierr /= 0) then
            write(*,*) 'failed to read 2nd line of ' // trim(results_fname)
            write(*,'(a)') 'line <' // trim(line) // '>'
            call done
            return
         end if
         
         do while (max_num_samples < num)
            call alloc_sample_ptrs(ierr)
            if (ierr /= 0) then
               write(*,*) 'ERROR -- failed to allocate for samples'
               call done
               return
            end if
         end do
         
         do j = 1, num
            call read1_sample_from_file(j, iounit, ierr)
            if (ierr /= 0) then
               write(*,2) 'ERROR -- failed while reading sample on line', j+2
               call done
               return
            end if
         end do
                  
         sample_number = num
         write(*,2) 'number of samples read from file', num
         
         call done

         
         contains
         
         
         subroutine done
            close(iounit)
            call free_iounit(iounit)         
         end subroutine done
         

      end subroutine read_samples_from_file
      
      
      subroutine read1_sample_from_file(j, iounit, ierr)
         use num_lib, only: simplex_op_code
         integer, intent(in) :: j, iounit
         integer, intent(out) :: ierr
            
         integer :: i, k
         character (len=256) :: info_str
         real(dp) :: logR
         
         include 'formats'
         
         ierr = 0
         read(iounit,fmt='(i8)',advance='no',iostat=ierr) i
         if (ierr /= 0) return
         if (i <= 0 .or. i > size(sample_chi2,dim=1)) then
            write(*,2) 'invalid sample number', i
            ierr = -1
            return
         end if
         
         read(iounit,'(7(1pd26.16),i16,99(1pd26.16))',advance='no',iostat=ierr) &
            sample_chi2(i), &
            sample_mass(i), &
            sample_init_Y(i), &
            sample_init_FeH(i), &
            sample_alpha(i), &
            sample_f_ov(i), &
            sample_age(i), &
            sample_model_number(i), &
            sample_init_h1(i), &
            sample_init_he3(i), &
            sample_init_he4(i), &
            sample_init_Z(i), &
            logR, &
            sample_logL(i), &
            sample_Teff(i), &
            sample_logg(i), &
            sample_FeH(i), &
            sample_logR(i), &
            sample_surface_Z_div_X(i), &
            sample_surface_He(i), &
            sample_Rcz(i), &
            sample_solar_cs_rms(i), &
            sample_my_var1(i), &
            sample_my_var2(i), &
            sample_my_var3(i), &
            sample_my_param1(i), &
            sample_my_param2(i), &
            sample_my_param3(i)
         if (failed('results')) return
            
         sample_radius(i) = exp10(logR)

         read(iounit,'(a12)',iostat=ierr) info_str
         if (ierr /= 0) then
            ierr = 0
            sample_op_code(i) = 0
            return
         end if
      
         if (len_trim(info_str) == 0) then
            sample_op_code(i) = 0
         else
            sample_op_code(i) = simplex_op_code(info_str, ierr)
            if (ierr /= 0) then
               ierr = 0
               sample_op_code(i) = 0
               return
            end if
         end if
         
         
         contains
         
         
         logical function failed(str)
            character (len=*), intent(in) :: str
            include 'formats'
            failed = .false.
            if (ierr == 0) return
            write(*,2) 'failed reading ' // trim(str) // ' data for sample number', i
            failed = .true.
         end function failed
         
      
      end subroutine read1_sample_from_file
      


      end module simplex_search_run_support
