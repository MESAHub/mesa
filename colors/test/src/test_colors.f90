program test_colors
   use colors_lib, only: colors_init, colors_shutdown
   implicit none

   integer :: ierr
   logical :: use_cache

   ierr = 0
   use_cache = .false.

   ! TODO: add tests for colors module functionality here

   write(*,*) 'Testing colors module initialization...'

   ! Initialize colors module
   ! call colors_init(use_cache, '', ierr)
   if (ierr /= 0) then
      write(*,*) 'Error: colors_init failed with code', ierr
      stop 1
   end if

   ! enable photometry so how_many_colors_history_columns returns > 0
   cs%use_colors = .true.
   cs%mag_system = 'Vega'

   n_cols = how_many_colors_history_columns(handle)
   if (n_cols == 0) then
      write(*,*) 'how_many_colors_history_columns returned 0'
      stop 1
   end if

   allocate(col_names(n_cols), col_vals(n_cols))
   model_num = 0

   ! -----------------------------------------------------------------------
   ! group 1: representative stellar types
   ! -----------------------------------------------------------------------

   write(*,'(a)') '# group1  system=Vega  grid=Kurucz2003  filters=Johnson'
   do j = 1, n_cases
      model_num = model_num + 1
      call data_for_colors_history_columns( &
         test_teff(j), test_logg(j), test_R(j), test_meta(j), model_num, &
         handle, n_cols, col_names, col_vals, ierr)
      if (ierr /= 0) then
         write(*,*) 'data_for_colors_history_columns failed, group1 case', j, ', ierr =', ierr
         stop 1
      end if
      write(*,'(a, a)') '# case: ', trim(adjustl(labels(j)))
      do k = 1, n_cols
         write(*,'(a40, 1pe24.14)') trim(col_names(k)), col_vals(k)
      end do
   end do

   ! -----------------------------------------------------------------------
   ! group 2a: vary [M/H]
   ! -----------------------------------------------------------------------

   write(*,'(a)') '# group2a  vary_MH  Teff=5778  logg=4.44'
   do j = 1, n_meta
      model_num = model_num + 1
      call data_for_colors_history_columns( &
         5778d0, 4.44d0, rsun, sweep_meta(j), model_num, &
         handle, n_cols, col_names, col_vals, ierr)
      if (ierr /= 0) then
         write(*,*) 'data_for_colors_history_columns failed, group2a case', j, ', ierr =', ierr
         stop 1
      end if
      write(*,'(a, f6.2)') '# MH= ', sweep_meta(j)
      do k = 1, n_cols
         write(*,'(a40, 1pe24.14)') trim(col_names(k)), col_vals(k)
      end do
   end do

   ! -----------------------------------------------------------------------
   ! group 2b: vary log g
   ! -----------------------------------------------------------------------

   write(*,'(a)') '# group2b  vary_logg  Teff=5778  MH=0.0'
   do j = 1, n_logg
      model_num = model_num + 1
      call data_for_colors_history_columns( &
         5778d0, sweep_logg(j), rsun, 0.0d0, model_num, &
         handle, n_cols, col_names, col_vals, ierr)
      if (ierr /= 0) then
         write(*,*) 'data_for_colors_history_columns failed, group2b case', j, ', ierr =', ierr
         stop 1
      end if
      write(*,'(a, f6.2)') '# logg= ', sweep_logg(j)
      do k = 1, n_cols
         write(*,'(a40, 1pe24.14)') trim(col_names(k)), col_vals(k)
      end do
   end do

   ! -----------------------------------------------------------------------
   ! group 2c: vary Teff
   ! -----------------------------------------------------------------------

   write(*,'(a)') '# group2c  vary_Teff  logg=4.0  MH=0.0'
   do j = 1, n_teff
      model_num = model_num + 1
      call data_for_colors_history_columns( &
         sweep_teff(j), 4.0d0, rsun, 0.0d0, model_num, &
         handle, n_cols, col_names, col_vals, ierr)
      if (ierr /= 0) then
         write(*,*) 'data_for_colors_history_columns failed, group2c case', j, ', ierr =', ierr
         stop 1
      end if
      write(*,'(a, f10.1)') '# Teff= ', sweep_teff(j)
      do k = 1, n_cols
         write(*,'(a40, 1pe24.14)') trim(col_names(k)), col_vals(k)
      end do
   end do

   ! -----------------------------------------------------------------------
   ! SED comparison: solar case, n_sed_samples evenly spaced wavelength/flux
   ! -----------------------------------------------------------------------

   sed_filepath = trim(mesa_dir)//trim(cs%stellar_atm)
   call calculate_bolometric( &
      cs, test_teff(1), test_logg(1), test_meta(1), test_R(1), d_10pc, &
      bol_mag, bol_flux, wavelengths, fluxes, sed_filepath, interp_rad)

   n_wav = size(wavelengths)
   stride = max(1, n_wav / n_sed_samples)

   write(*,'(a)') '# SED sample  case=solar  columns=wavelength_AA  flux_erg_s_cm2_AA'
   do i = 1, n_wav, stride
      write(*,'(1pe24.14, 1x, 1pe24.14)') wavelengths(i), fluxes(i)
   end do

   ! -----------------------------------------------------------------------
   ! cleanup
   ! -----------------------------------------------------------------------

   deallocate(col_names, col_vals)
   if (allocated(wavelengths)) deallocate(wavelengths)
   if (allocated(fluxes))      deallocate(fluxes)

   call free_colors_handle(handle)
   call colors_shutdown()

   ! Clean up
   ! call colors_shutdown()

end program test_colors