program eos_composition_partials_check

   use const_def, only: dp, crad
   use const_lib, only: const_init
   use chem_def
   use chem_lib, only: chem_init
   use eos_def
   use eos_lib
   use math_lib

   implicit none

   integer, parameter :: species = 14
   integer, parameter :: h1=1, he3=2, he4=3, c12=4, n14=5, o16=6, ne20=7
   integer, parameter :: f20=8, o20=9, mg24=10, na24=11, ne24=12, si28=13, fe56=14
   integer, parameter :: sink = he4
   integer, parameter :: nsample_dirs = 3
   integer, parameter :: nsweeps = 4
   integer, parameter :: npoints = 31
   integer, parameter :: ncontour_logT = 150
   integer, parameter :: ncontour_logRho = 150
   integer, parameter :: ndfridr_logT = 60
   integer, parameter :: ndfridr_logRho = 60
   real(dp), parameter :: contour_logT_min = 2d0
   real(dp), parameter :: contour_logT_max = 10d0
   real(dp), parameter :: contour_logRho_min = -15d0
   real(dp), parameter :: contour_logRho_max = 10d0
   integer, parameter :: ncontour_rows = 13
   integer, parameter :: ncontour_dirs = 3
   integer, parameter :: contour_rows(ncontour_rows) = [ &
      i_mu, i_lnPgas, i_lnE, i_lnS, i_lnfree_e, i_eta, i_grad_ad, &
      i_chiRho, i_chiT, i_Cp, i_Cv, i_gamma1, i_gamma3]
   integer, parameter :: sample_dirs(nsample_dirs) = [h1, c12, o16]
   integer, parameter :: contour_dirs(ncontour_dirs) = [he4, c12, o16]
   integer, parameter :: ndfridr_rows = 3
   integer, parameter :: ndfridr_dirs = 3
   integer, parameter :: dfridr_rows(ndfridr_rows) = [0, i_mu, i_lnE]
   integer, parameter :: dfridr_dirs(ndfridr_dirs) = [he4, c12, o16]
   integer, parameter :: dfridr_sinks(ndfridr_dirs) = [h1, h1, h1]
   integer, parameter :: ridders_ntab = 20
   real(dp), parameter :: ridders_con2 = 2d0
   real(dp), parameter :: ridders_con = sqrt(ridders_con2)
   real(dp), parameter :: ridders_safe = 2d0
   real(dp), parameter :: ridders_big = 1d50
   real(dp), parameter :: ridders_initial_step = 1d-4
   integer, target :: chem_id_array(species)
   integer, target, allocatable :: net_iso_array(:)
   integer, pointer :: chem_id(:), net_iso(:)
   character(len=12) :: species_name(species)
   character(len=18) :: sweep_name(nsweeps)
   real(dp) :: xa(species)
   integer :: handle, ierr, unit

   call initialize(ierr)
   if (ierr /= 0) stop 1

   call set_composition
   call set_sweeps

   open(newunit=unit, file='data/eos_composition_partials.csv', status='replace', action='write')
   write(unit,'(a)') &
      'sweep,point,logT,logRho,direction,sink,row,row_name,' // &
      'analytic,finite_difference,finite_difference_error,abs_error,rel_error,' // &
      'frac_OPAL_SCVH,frac_HELM,frac_Skye,frac_PC,frac_FreeEOS,frac_CMS,frac_ideal'
   call write_samples(unit)
   close(unit)

   open(newunit=unit, file='data/eos_composition_partials_contours.csv', status='replace', action='write')
   write(unit,'(a)') &
      'logT,logRho,direction,sink,row,row_name,raw_partial,sink_projected_partial,' // &
      'frac_OPAL_SCVH,frac_HELM,frac_Skye,frac_PC,frac_FreeEOS,frac_CMS,frac_ideal'
   call write_contours(unit)
   close(unit)

   open(newunit=unit, file='data/eos_composition_partials_dfridr_contours.csv', &
      status='replace', action='write')
   write(unit,'(a)') &
      'logT,logRho,direction,sink,row,row_name,' // &
      'analytic,finite_difference,finite_difference_error,abs_error,rel_error,' // &
      'frac_OPAL_SCVH,frac_HELM,frac_Skye,frac_PC,frac_FreeEOS,frac_CMS,frac_ideal'
   call write_dfridr_contours(unit)
   close(unit)

   call free_eos_handle(handle)
   call eos_shutdown

   write(*,'(a)') 'wrote data/eos_composition_partials.csv'
   write(*,'(a)') 'wrote data/eos_composition_partials_contours.csv'
   write(*,'(a)') 'wrote data/eos_composition_partials_dfridr_contours.csv'

contains

   subroutine initialize(ierr)
      integer, intent(out) :: ierr
      logical :: use_cache

      ierr = 0
      species_name = ['h1        ', 'he3       ', 'he4       ', 'c12       ', 'n14       ', &
         'o16       ', 'ne20      ', 'f20       ', 'o20       ', 'mg24      ', 'na24      ', &
         'ne24      ', 'si28      ', 'fe56      ']

      call const_init('../../..', ierr)
      if (ierr /= 0) then
         write(*,*) 'const_init failed'
         return
      end if

      call math_init()

      call chem_init('isotopes.data', ierr)
      if (ierr /= 0) then
         write(*,*) 'chem_init failed'
         return
      end if

      allocate(net_iso_array(num_chem_isos))
      chem_id => chem_id_array
      net_iso => net_iso_array
      net_iso = 0

      chem_id(h1) = ih1; net_iso(ih1) = h1
      chem_id(he3) = ihe3; net_iso(ihe3) = he3
      chem_id(he4) = ihe4; net_iso(ihe4) = he4
      chem_id(c12) = ic12; net_iso(ic12) = c12
      chem_id(n14) = in14; net_iso(in14) = n14
      chem_id(o16) = io16; net_iso(io16) = o16
      chem_id(ne20) = ine20; net_iso(ine20) = ne20
      chem_id(f20) = if20; net_iso(if20) = f20
      chem_id(o20) = io20; net_iso(io20) = o20
      chem_id(mg24) = img24; net_iso(img24) = mg24
      chem_id(na24) = ina24; net_iso(ina24) = na24
      chem_id(ne24) = ine24; net_iso(ine24) = ne24
      chem_id(si28) = isi28; net_iso(isi28) = si28
      chem_id(fe56) = ife56; net_iso(ife56) = fe56

      use_cache = .true.
      call eos_init(' ', use_cache, ierr)
      if (ierr /= 0) then
         write(*,*) 'eos_init failed'
         return
      end if

      handle = alloc_eos_handle(ierr)
      if (ierr /= 0) write(*,*) 'alloc_eos_handle failed'
   end subroutine initialize

   subroutine set_composition
      xa = 0d0
      xa(h1) = 0.70d0
      xa(c12) = 0.01d0
      xa(o16) = 0.01d0
      xa(he4) = 1d0 - xa(h1) - xa(c12) - xa(o16)
   end subroutine set_composition

   subroutine set_sweeps
      sweep_name(1) = 'hot_T_sweep'
      sweep_name(2) = 'dense_T_sweep'
      sweep_name(3) = 'cool_T_sweep'
      sweep_name(4) = 'rho_sweep'
   end subroutine set_sweeps

   subroutine get_point(isweep, ipoint, logT, logRho)
      integer, intent(in) :: isweep, ipoint
      real(dp), intent(out) :: logT, logRho
      real(dp) :: f

      f = dble(ipoint - 1)/dble(npoints - 1)
      select case(isweep)
      case(1)
         logRho = 2.0d0
         logT = 6.5d0 + f*(9.0d0 - 6.5d0)
      case(2)
         logRho = 5.5d0
         logT = 6.5d0 + f*(8.5d0 - 6.5d0)
      case(3)
         logRho = -4.0d0
         logT = 4.0d0 + f*(6.0d0 - 4.0d0)
      case(4)
         logT = 7.0d0
         logRho = -6.0d0 + f*(6.0d0 + 6.0d0)
      end select
   end subroutine get_point

   subroutine write_samples(unit)
      integer, intent(in) :: unit

      integer :: isweep, ipoint, idir, idir_index, row, ierr
      real(dp) :: logT, logRho, T, Rho, analytic, abs_error, rel_error, denom
      real(dp) :: fd(num_eos_basic_results), fd_err(num_eos_basic_results)
      real(dp) :: res0(num_eos_basic_results), d_dlnd0(num_eos_basic_results)
      real(dp) :: d_dlnT0(num_eos_basic_results), d_dxa0(num_eos_basic_results,species)

      do isweep = 1, nsweeps
         do ipoint = 1, npoints
            call get_point(isweep, ipoint, logT, logRho)
            T = exp10(logT)
            Rho = exp10(logRho)
            call eval_eos(xa, Rho, logRho, T, logT, res0, d_dlnd0, d_dlnT0, d_dxa0, ierr)
            if (ierr /= 0) then
               write(*,'(a,1x,i0,1x,i0)') 'base eos failure', isweep, ipoint
               cycle
            end if

            do idir_index = 1, nsample_dirs
               idir = sample_dirs(idir_index)
               if (idir == sink) cycle
               call finite_difference_direction(idir, sink, Rho, logRho, T, logT, fd, fd_err, ierr)
               if (ierr /= 0) then
                  write(*,'(a,1x,i0,1x,i0,1x,i0)') 'finite-difference eos failure', &
                     isweep, ipoint, idir
                  cycle
               end if

               do row = 1, num_eos_basic_results
                  analytic = d_dxa0(row,idir) - d_dxa0(row,sink)
                  abs_error = abs(analytic - fd(row))
                  denom = max(1d-14, abs(analytic), abs(fd(row)))
                  rel_error = abs_error/denom
                  write(unit,'(a,",",i0,",",es24.16,",",es24.16,",",a,",",a,",",i0,",",a,",",12(es24.16,:,","))') &
                     trim(sweep_name(isweep)), ipoint, logT, logRho, trim(species_name(idir)), &
                     trim(species_name(sink)), row, trim(eosDT_result_names(row)), &
                     analytic, fd(row), fd_err(row), abs_error, rel_error, &
                     res0(i_frac_OPAL_SCVH), res0(i_frac_HELM), res0(i_frac_Skye), &
                     res0(i_frac_PC), res0(i_frac_FreeEOS), res0(i_frac_CMS), res0(i_frac_ideal)
               end do
            end do
         end do
      end do
   end subroutine write_samples

   subroutine write_contours(unit)
      integer, intent(in) :: unit

      integer :: iT, iRho, irow, idir, row, dir, ierr, num_failures
      real(dp) :: fT, fRho, logT, logRho, T, Rho, raw_partial, sink_projected_partial
      real(dp) :: Pgas, Prad, Peos
      real(dp) :: res0(num_eos_basic_results), d_dlnd0(num_eos_basic_results)
      real(dp) :: d_dlnT0(num_eos_basic_results), d_dxa0(num_eos_basic_results,species)

      num_failures = 0
      do iT = 1, ncontour_logT
         fT = dble(iT - 1)/dble(ncontour_logT - 1)
         logT = contour_logT_min + fT*(contour_logT_max - contour_logT_min)
         T = exp10(logT)
         do iRho = 1, ncontour_logRho
            fRho = dble(iRho - 1)/dble(ncontour_logRho - 1)
            logRho = contour_logRho_min + fRho*(contour_logRho_max - contour_logRho_min)
            Rho = exp10(logRho)
            call eval_eos(xa, Rho, logRho, T, logT, res0, d_dlnd0, d_dlnT0, d_dxa0, ierr)
            if (ierr /= 0) then
               num_failures = num_failures + 1
               cycle
            end if

            Pgas = exp(res0(i_lnPgas))
            Prad = crad*T*T*T*T/3d0
            Peos = Pgas + Prad

            do idir = 1, ncontour_dirs
               dir = contour_dirs(idir)
               do irow = 1, ncontour_rows
                  row = contour_rows(irow)
                  raw_partial = d_dxa0(row,dir)
                  sink_projected_partial = raw_partial - d_dxa0(row,sink)
                  write(unit,'(es24.16,",",es24.16,",",a,",",a,",",i0,",",a,",",9(es24.16,:,","))') &
                     logT, logRho, trim(species_name(dir)), trim(species_name(sink)), &
                     row, trim(eosDT_result_names(row)), raw_partial, sink_projected_partial, &
                     res0(i_frac_OPAL_SCVH), res0(i_frac_HELM), res0(i_frac_Skye), &
                     res0(i_frac_PC), res0(i_frac_FreeEOS), res0(i_frac_CMS), res0(i_frac_ideal)
               end do
               raw_partial = (Pgas/Peos)*d_dxa0(i_lnPgas,dir)
               sink_projected_partial = (Pgas/Peos)*(d_dxa0(i_lnPgas,dir) - d_dxa0(i_lnPgas,sink))
               write(unit,'(es24.16,",",es24.16,",",a,",",a,",",i0,",",a,",",9(es24.16,:,","))') &
                  logT, logRho, trim(species_name(dir)), trim(species_name(sink)), &
                  -1, 'chiX', raw_partial, sink_projected_partial, &
                  res0(i_frac_OPAL_SCVH), res0(i_frac_HELM), res0(i_frac_Skye), &
                  res0(i_frac_PC), res0(i_frac_FreeEOS), res0(i_frac_CMS), res0(i_frac_ideal)
            end do
         end do
      end do
      if (num_failures > 0) write(*,'(a,1x,i0)') 'contour eos failures', num_failures
   end subroutine write_contours

   subroutine write_dfridr_contours(unit)
      integer, intent(in) :: unit

      integer :: iT, iRho, idir, irow, dir, sink_in, row, ierr, num_base_failures, num_fd_failures
      real(dp) :: fT, fRho, logT, logRho, T, Rho, analytic, abs_error, rel_error, denom
      real(dp) :: Pgas
      real(dp) :: fd(ndfridr_rows), fd_err(ndfridr_rows)
      real(dp) :: res0(num_eos_basic_results), d_dlnd0(num_eos_basic_results)
      real(dp) :: d_dlnT0(num_eos_basic_results), d_dxa0(num_eos_basic_results,species)
      character(len=eos_name_length) :: row_name

      num_base_failures = 0
      num_fd_failures = 0
      do iT = 1, ndfridr_logT
         fT = dble(iT - 1)/dble(ndfridr_logT - 1)
         logT = contour_logT_min + fT*(contour_logT_max - contour_logT_min)
         T = exp10(logT)
         do iRho = 1, ndfridr_logRho
            fRho = dble(iRho - 1)/dble(ndfridr_logRho - 1)
            logRho = contour_logRho_min + fRho*(contour_logRho_max - contour_logRho_min)
            Rho = exp10(logRho)
            call eval_eos(xa, Rho, logRho, T, logT, res0, d_dlnd0, d_dlnT0, d_dxa0, ierr)
            if (ierr /= 0) then
               num_base_failures = num_base_failures + 1
               cycle
            end if

            Pgas = exp(res0(i_lnPgas))
            do idir = 1, ndfridr_dirs
               dir = dfridr_dirs(idir)
               sink_in = dfridr_sinks(idir)
               call finite_difference_dfridr_direction(dir, sink_in, Rho, logRho, T, logT, fd, fd_err, ierr)
               if (ierr /= 0) then
                  num_fd_failures = num_fd_failures + 1
                  cycle
               end if

               do irow = 1, ndfridr_rows
                  row = dfridr_rows(irow)
                  call get_dfridr_row_name(row, row_name)
                  if (row == 0) then
                     analytic = Pgas*(d_dxa0(i_lnPgas,dir) - d_dxa0(i_lnPgas,sink_in))
                  else
                     analytic = d_dxa0(row,dir) - d_dxa0(row,sink_in)
                  end if
                  abs_error = abs(analytic - fd(irow))
                  denom = max(1d-14, abs(analytic), abs(fd(irow)))
                  rel_error = abs_error/denom
                  write(unit,'(es24.16,",",es24.16,",",a,",",a,",",i0,",",a,",",12(es24.16,:,","))') &
                     logT, logRho, trim(species_name(dir)), trim(species_name(sink_in)), &
                     row, trim(row_name), analytic, fd(irow), fd_err(irow), abs_error, rel_error, &
                     res0(i_frac_OPAL_SCVH), res0(i_frac_HELM), res0(i_frac_Skye), &
                     res0(i_frac_PC), res0(i_frac_FreeEOS), res0(i_frac_CMS), res0(i_frac_ideal)
               end do
            end do
         end do
      end do
      if (num_base_failures > 0) write(*,'(a,1x,i0)') 'dfridr base eos failures', num_base_failures
      if (num_fd_failures > 0) write(*,'(a,1x,i0)') 'dfridr finite-difference eos failures', num_fd_failures
   end subroutine write_dfridr_contours

   subroutine get_dfridr_row_name(row, row_name)
      integer, intent(in) :: row
      character(len=eos_name_length), intent(out) :: row_name

      if (row == 0) then
         row_name = 'Pgas'
      else
         row_name = eosDT_result_names(row)
      end if
   end subroutine get_dfridr_row_name

   subroutine finite_difference_direction(idir, sink_in, Rho, logRho, T, logT, fd, fd_err, ierr)
      integer, intent(in) :: idir, sink_in
      real(dp), intent(in) :: Rho, logRho, T, logT
      real(dp), intent(out) :: fd(num_eos_basic_results), fd_err(num_eos_basic_results)
      integer, intent(out) :: ierr

      integer :: i, j, row
      real(dp) :: hh, fac, errt
      real(dp) :: a(num_eos_basic_results,ridders_ntab,ridders_ntab)
      real(dp) :: res_plus(num_eos_basic_results), res_minus(num_eos_basic_results)

      ierr = 0
      fd = 0d0
      fd_err = ridders_big
      a = 0d0

      if (idir == sink_in) then
         ierr = -1
         return
      end if

      hh = min(ridders_initial_step, 2d-1*min(xa(idir), xa(sink_in)))
      if (hh <= 0d0) then
         ierr = -1
         return
      end if

      call eval_shifted_eos(idir, sink_in, hh, Rho, logRho, T, logT, res_plus, ierr)
      if (ierr /= 0) return
      call eval_shifted_eos(idir, sink_in, -hh, Rho, logRho, T, logT, res_minus, ierr)
      if (ierr /= 0) return

      a(:,1,1) = (res_plus - res_minus)/(2d0*hh)
      fd = a(:,1,1)

      do i = 2, ridders_ntab
         hh = hh/ridders_con
         call eval_shifted_eos(idir, sink_in, hh, Rho, logRho, T, logT, res_plus, ierr)
         if (ierr /= 0) return
         call eval_shifted_eos(idir, sink_in, -hh, Rho, logRho, T, logT, res_minus, ierr)
         if (ierr /= 0) return

         a(:,1,i) = (res_plus - res_minus)/(2d0*hh)
         fac = ridders_con2
         do j = 2, i
            a(:,j,i) = (a(:,j-1,i)*fac - a(:,j-1,i-1))/(fac - 1d0)
            do row = 1, num_eos_basic_results
               errt = max(abs(a(row,j,i) - a(row,j-1,i)), abs(a(row,j,i) - a(row,j-1,i-1)))
               if (errt <= fd_err(row)) then
                  fd_err(row) = errt
                  fd(row) = a(row,j,i)
               end if
            end do
            fac = ridders_con2*fac
         end do

         if (all(abs(a(:,i,i) - a(:,i-1,i-1)) >= ridders_safe*fd_err)) return
      end do
   end subroutine finite_difference_direction

   subroutine finite_difference_dfridr_direction(idir, sink_in, Rho, logRho, T, logT, fd, fd_err, ierr)
      integer, intent(in) :: idir, sink_in
      real(dp), intent(in) :: Rho, logRho, T, logT
      real(dp), intent(out) :: fd(ndfridr_rows), fd_err(ndfridr_rows)
      integer, intent(out) :: ierr

      integer :: i, j, row
      real(dp) :: hh, fac, errt
      real(dp) :: a(ndfridr_rows,ridders_ntab,ridders_ntab)
      real(dp) :: vals_plus(ndfridr_rows), vals_minus(ndfridr_rows)

      ierr = 0
      fd = 0d0
      fd_err = ridders_big
      a = 0d0
      if (idir == sink_in) then
         ierr = -1
         return
      end if

      hh = min(ridders_initial_step, 2d-1*min(xa(idir), xa(sink_in)))
      if (hh <= 0d0) then
         ierr = -1
         return
      end if

      call eval_shifted_dfridr_values(idir, sink_in, hh, Rho, logRho, T, logT, vals_plus, ierr)
      if (ierr /= 0) return
      call eval_shifted_dfridr_values(idir, sink_in, -hh, Rho, logRho, T, logT, vals_minus, ierr)
      if (ierr /= 0) return

      a(:,1,1) = (vals_plus - vals_minus)/(2d0*hh)
      fd = a(:,1,1)

      do i = 2, ridders_ntab
         hh = hh/ridders_con
         call eval_shifted_dfridr_values(idir, sink_in, hh, Rho, logRho, T, logT, vals_plus, ierr)
         if (ierr /= 0) return
         call eval_shifted_dfridr_values(idir, sink_in, -hh, Rho, logRho, T, logT, vals_minus, ierr)
         if (ierr /= 0) return

         a(:,1,i) = (vals_plus - vals_minus)/(2d0*hh)
         fac = ridders_con2
         do j = 2, i
            a(:,j,i) = (a(:,j-1,i)*fac - a(:,j-1,i-1))/(fac - 1d0)
            do row = 1, ndfridr_rows
               errt = max(abs(a(row,j,i) - a(row,j-1,i)), abs(a(row,j,i) - a(row,j-1,i-1)))
               if (errt <= fd_err(row)) then
                  fd_err(row) = errt
                  fd(row) = a(row,j,i)
               end if
            end do
            fac = ridders_con2*fac
         end do

         if (all(abs(a(:,i,i) - a(:,i-1,i-1)) >= ridders_safe*fd_err)) return
      end do
   end subroutine finite_difference_dfridr_direction

   subroutine eval_shifted_eos(idir, sink_in, delta_x, Rho, logRho, T, logT, res, ierr)
      integer, intent(in) :: idir, sink_in
      real(dp), intent(in) :: delta_x, Rho, logRho, T, logT
      real(dp), intent(out) :: res(num_eos_basic_results)
      integer, intent(out) :: ierr

      real(dp) :: x(species)

      x = xa
      x(idir) = x(idir) + delta_x
      x(sink_in) = x(sink_in) - delta_x
      if (minval(x) < -1d-12) then
         ierr = -1
         return
      end if

      call eval_eos_values(x, Rho, logRho, T, logT, res, ierr)
   end subroutine eval_shifted_eos

   subroutine eval_shifted_dfridr_values(idir, sink_in, delta_x, Rho, logRho, T, logT, vals, ierr)
      integer, intent(in) :: idir, sink_in
      real(dp), intent(in) :: delta_x, Rho, logRho, T, logT
      real(dp), intent(out) :: vals(ndfridr_rows)
      integer, intent(out) :: ierr

      real(dp) :: x(species)

      x = xa
      x(idir) = x(idir) + delta_x
      x(sink_in) = x(sink_in) - delta_x
      if (minval(x) < -1d-12) then
         ierr = -1
         return
      end if

      call eval_dfridr_values(x, Rho, logRho, T, logT, vals, ierr)
   end subroutine eval_shifted_dfridr_values

   subroutine eval_dfridr_values(x, Rho, logRho, T, logT, vals, ierr)
      real(dp), intent(in) :: x(species), Rho, logRho, T, logT
      real(dp), intent(out) :: vals(ndfridr_rows)
      integer, intent(out) :: ierr

      real(dp) :: res(num_eos_basic_results)

      call eval_eos_values(x, Rho, logRho, T, logT, res, ierr)
      if (ierr /= 0) return
      vals(1) = exp(res(i_lnPgas))
      vals(2) = res(i_mu)
      vals(3) = res(i_lnE)
   end subroutine eval_dfridr_values

   subroutine eval_eos(x, Rho, logRho, T, logT, res, d_dlnd, d_dlnT, d_dxa, ierr)
      real(dp), intent(in) :: x(species), Rho, logRho, T, logT
      real(dp), intent(out) :: res(num_eos_basic_results), d_dlnd(num_eos_basic_results)
      real(dp), intent(out) :: d_dlnT(num_eos_basic_results)
      real(dp), intent(out) :: d_dxa(num_eos_basic_results,species)
      integer, intent(out) :: ierr

      res = 0d0
      d_dlnd = 0d0
      d_dlnT = 0d0
      d_dxa = 0d0
      call eosDT_get_full_dxa( &
         handle, species, chem_id, net_iso, x, &
         Rho, logRho, T, logT, res, d_dlnd, d_dlnT, d_dxa, ierr)
   end subroutine eval_eos

   subroutine eval_eos_values(x, Rho, logRho, T, logT, res, ierr)
      real(dp), intent(in) :: x(species), Rho, logRho, T, logT
      real(dp), intent(out) :: res(num_eos_basic_results)
      integer, intent(out) :: ierr

      real(dp) :: d_dlnd(num_eos_basic_results), d_dlnT(num_eos_basic_results)
      real(dp) :: d_dxa(num_eos_d_dxa_results,species)

      res = 0d0
      d_dlnd = 0d0
      d_dlnT = 0d0
      d_dxa = 0d0
      call eosDT_get( &
         handle, species, chem_id, net_iso, x, &
         Rho, logRho, T, logT, res, d_dlnd, d_dlnT, d_dxa, ierr)
   end subroutine eval_eos_values

end program eos_composition_partials_check
