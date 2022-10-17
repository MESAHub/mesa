module chem_support
   use const_def
   use const_lib
   use math_lib
   
   real(dp), parameter :: no_mass_table_entry = -999.d0
   character(len = 256) :: data_dir, output_dir, masstable_filename, winvn_filename
   integer :: masstable_header_length
   integer :: winvn_header_length
   integer :: number_nuclides
   ! integer, parameter :: number_masstable_header_lines = 15
   ! integer, parameter :: number_winvn_header_lines = 7857
   ! character(len=*), parameter :: masstable_filename = 'masslib_library_5.data'
   real(dp), dimension(:, :), allocatable :: mass_table
   
   
   ! al26 partition functions from Gupta & Meyer (2001)
   ! http://ftp.aip.org/epaps//phys_rev_c/E-PRVCAN-64-028108/
   ! Note that these values do not match their Table IV.
   
   real(dp), dimension(24) :: partiton_al26_1 = &
      (/ 1.00000E+00_dp, 1.00000E+00_dp, 1.00000E+00_dp, 1.00000E+00_dp, 1.00000E+00_dp, 1.00000E+00_dp, 1.00020E+00_dp, 1.00060E+00_dp, &
         1.00150E+00_dp, 1.00290E+00_dp, 1.00500E+00_dp, 1.02440E+00_dp, 1.03250E+00_dp, 1.01910E+00_dp, 1.01260E+00_dp, 1.00940E+00_dp, &
         1.00780E+00_dp, 1.00760E+00_dp, 1.00870E+00_dp, 1.01620E+00_dp, 1.03310E+00_dp, 1.06000E+00_dp, 1.09270E+00_dp, 1.12850E+00_dp /)
   
   real(dp), dimension(24) :: partiton_al26_2 = &
      (/ 1.00000E+00_dp, 1.00000E+00_dp, 1.00000E+00_dp, 1.00010E+00_dp, 1.00030E+00_dp, 1.00090E+00_dp, 1.00180E+00_dp, 1.00310E+00_dp, &
         1.00460E+00_dp, 1.00620E+00_dp, 1.00800E+00_dp, 1.06220E+00_dp, 2.02360E+00_dp, 3.38330E+00_dp, 4.19860E+00_dp, 4.81450E+00_dp, &
         5.36050E+00_dp, 5.88330E+00_dp, 6.40960E+00_dp, 7.52430E+00_dp, 8.75760E+00_dp, 1.01610E+01_dp, 1.18190E+01_dp, 1.37730E+00_dp /)
   
   namelist /chem/  &
      & data_dir, output_dir, masstable_filename, winvn_filename, &
      &  masstable_header_length, winvn_header_length, number_nuclides

contains
   
   subroutine read_input_parameters(inlist_fname)
      character(len = *), intent(in) :: inlist_fname
      integer :: iounit, ios
      
      ! set default values
      data_dir = 'chem_input_data'
      output_dir = 'data/chem_data'
      masstable_filename = 'masslib_library_5.data'
      winvn_filename = 'winvne_v2.data'
      masstable_header_length = 15
      winvn_header_length = 4
      number_nuclides = 7853
      
      open(newunit = iounit, file = trim(inlist_fname), iostat = ios, status = "old", action = "read", delim = 'quote')
      if (ios /= 0) then
         write(*, '(A)')
         write(*, '(A)')
         write(*, *) 'Failed to open namelist file ', trim(inlist_fname)
         write(*, '(A)')
         write(*, '(A)')
         stop
      end if
      read(iounit, nml = chem, iostat = ios)
      if (ios /= 0) then
         write (*, *)
         write (*, *)
         write (*, *) 'Failed to read namelist file ', trim(inlist_fname)
         write (*, *)
         write (*, *)
         stop
      end if
   end subroutine read_input_parameters
   
   subroutine init_preprocessor()
      integer :: ierr
      ierr = 0
      call const_init("../../", ierr)
      if(ierr/=0) stop "Error in const_init"
      
      call math_init()
   
   end subroutine init_preprocessor
   
   subroutine read_mass_table()
      integer :: mass_unit
      integer :: ios, i, nlines, zmin, zmax, nmin, nmax, Z, A, N, ierr
      real(dp), parameter :: keV_to_MeV = 1.0d-3
      real(dp) :: mass
      character(len = 256) :: filename, buf
      character(len = 24) :: eval, error
      
      write(filename, '(a)') trim(data_dir) // '/' // trim(masstable_filename)
      open(newunit = mass_unit, file = trim(filename), iostat = ios, status = "old", action = "read")
      
      if (ios /= 0) call mesa_error(__FILE__, __LINE__, 'unable to open mass table for reading')
      
      ! first pass
      call skip_header
      nlines = 0
      zmin = 999; zmax = -1; nmin = 999; nmax = -1
      do
         !read(mass_unit,*,iostat = ios) Z,A, eval, mass, error
         read(mass_unit, '(A)', iostat = ios) buf
         if (ios /= 0) exit
         call parse_line_mass_unit(buf, Z, A, eval, mass, error, ierr)
         nlines = nlines + 1
         N = A - Z
         if (Z < zmin) zmin = Z
         if (Z > zmax) zmax = Z
         if (N < nmin) nmin = N
         if (N > nmax) nmax = N
      end do
      
      ! now read the table
      allocate(mass_table(zmin:zmax, nmin:nmax))
      mass_table = no_mass_table_entry
      rewind(mass_unit)
      call skip_header
      do i = 1, nlines
         read(mass_unit, '(A)', iostat = ios) buf
         call parse_line_mass_unit(buf, Z, A, eval, mass, error, ierr)
         N = A - Z
         ! store mass and convert from keV to MeV
         mass_table(Z, N) = mass * keV_to_MeV
      end do
      close(mass_unit)
   
   contains
      subroutine skip_header()
         integer :: i, ios
         do i = 1, masstable_header_length
            read(mass_unit, *, iostat = ios)
            if (ios /= 0) call mesa_error(__FILE__, __LINE__, 'error while reading header of masstable')
         end do
      end subroutine skip_header
   end subroutine read_mass_table
   
   subroutine process_winvn_table()
      use iso_fortran_env, only : error_unit
      integer :: in_unit, out_unit
      integer :: ios, i, ierr
      character(len = 256) :: infile_name, outfile_name, buf
      real(dp) :: fac
      integer :: zmin, zmax, nmin, nmax
      integer :: Z, N
      real(dp) :: W, spin, mass_excess, pfcn(24)
      character(len = 8) :: name, ref
      logical :: ground_state
      
      fac = mev_to_ergs / amu / (clight * clight)
      
      ! the mass table must be allocated first
      if (.not.allocated(mass_table)) then
         write (error_unit, *) 'mass_table must be allocated first'
         return
      end if
      
      zmin = lbound(mass_table, dim = 1)
      zmax = ubound(mass_table, dim = 1)
      nmin = lbound(mass_table, dim = 2)
      nmax = ubound(mass_table, dim = 2)
      
      write(infile_name, '(a)') trim(data_dir) // '/' // trim(winvn_filename)
      write(outfile_name, '(a)') trim(output_dir) // '/isotopes.data'
      
      open(newunit = in_unit, file = trim(infile_name), iostat = ios, status = "old", action = "read")
      if (ios /= 0) stop "Error opening raw winvn file"
      
      open(newunit = out_unit, file = trim(outfile_name), iostat = ios, action = "write")
      if (ios /= 0) stop "Error opening processed winvn file"
      
      !Add blank line to file at the start
      write (out_unit, *)
      
      ! skim off the header
      do i = 1, winvn_header_length + number_nuclides
         read(in_unit, *)
      end do
      
      do i = 1, number_nuclides  ! 4 lines per nuclide
         read(in_unit, '(A)', iostat = ios) buf
         if (ios /= 0) exit
         call parse_line_nuclides_header(buf, name, w, z, n, spin, mass_excess, ref, ierr)
         read(in_unit, '(A)', iostat = ios) buf
         call parse_line_pfcn(buf, pfcn(1:8), ierr)
         if (ios /= 0) exit
         read(in_unit, '(A)', iostat = ios) buf
         call parse_line_pfcn(buf, pfcn(9:16), ierr)
         if (ios /= 0) exit
         read(in_unit, '(A)', iostat = ios) buf
         call parse_line_pfcn(buf, pfcn(17:24), ierr)
         if (ios /= 0) exit
         
         if (name /= 'al-6' .and. name /= 'al*6') then
            ground_state = .true.
         else
            ground_state = .false.
         end if
         
         ! lookup the mass information if it is available
         if (Z >= zmin .and. Z <= zmax .and. N >= nmin .and. N <= nmax .and. ground_state) then
            if (mass_table(Z, N) /= no_mass_table_entry)  mass_excess = mass_table(Z, N)
         end if
         
         ! set the atomic weight
         W = Z + N + mass_excess * fac
         
         ! convert the name
         if (name == 'n') name = 'neut'
         if (name == 'p') name = 'h1'
         if (name == 'd') name = 'h2'
         if (name == 't') name = 'h3'
         if (name == 'al-6') then
            name = 'al26-1'
            pfcn = partiton_al26_1
         end if
         if (name == 'al*6') then
            name = 'al26-2'
            pfcn = partiton_al26_2
         end if
         
         ! write to the processed datafile
         call write_entry
         
         ! duplicate entry for h1 as prot
         if (name == 'h1') then
            name = 'prot'
            call write_entry
         end if
      end do
      ! write the 'xtra' entries for the rates
      name = 'xtra1'; W = 100.d0; Z = 0; N = 100; spin = 0.d0; mass_excess = 0.d0; pfcn = 1.d0
      call write_entry
      
      name = 'xtra2'; W = 200.d0; Z = 0; N = 200; spin = 0.d0; mass_excess = 0.d0; pfcn = 1.d0
      call write_entry
      
      close(in_unit)
      close(out_unit)
   
   contains
      subroutine write_entry()
         write (out_unit, '(a8,f13.7,i5,i5,f6.1,f14.9)') name, W, Z, N, spin, mass_excess
         write (out_unit, '(8(es12.5,tr1))') pfcn(1:8)
         write (out_unit, '(8(es12.5,tr1))') pfcn(9:16)
         write (out_unit, '(8(es12.5,tr1))') pfcn(17:24)
      end subroutine write_entry
   end subroutine process_winvn_table
   
   subroutine cleanup()
      deallocate(mass_table)
   end subroutine cleanup
   
   
   subroutine parse_line_mass_unit(line, z, a, eval, mass, error, ierr)
      character(len = *), intent(in) :: line
      integer, intent(out) :: z, a
      character(len = *), intent(out) :: eval, error
      real(dp), intent(out) :: mass
      integer, intent(inout) :: ierr
      
      integer :: j, k
      integer, parameter :: num_cols = 5
      character(len = 256), dimension(num_cols) :: tmp
      
      call parse_line(line, num_cols, tmp)
      
      read(tmp(1), *) z
      read(tmp(2), *) a
      eval = trim(tmp(3))
      call str_to_double(tmp(4), mass, ierr)
      if (ierr /= 0) return
      error = trim(tmp(5))
   
   end subroutine parse_line_mass_unit
   
   
   subroutine parse_line_nuclides_header(line, name, w, z, n, spin, mass, ref, ierr)
      character(len = *), intent(in) :: line
      integer, intent(out) :: z, n
      character(len = *), intent(out) :: ref, name
      real(dp), intent(out) :: w, spin, mass
      integer, intent(inout) :: ierr
      integer, parameter :: num_cols = 7
      character(len = 256), dimension(num_cols) :: tmp
      
      call parse_line(line, num_cols, tmp)
      name = tmp(1)
      call str_to_double(tmp(2), W, ierr)
      if (ierr /= 0) return
      read(tmp(3), *) Z
      if (ierr /= 0) return
      read(tmp(4), *) N
      if (ierr /= 0) return
      call str_to_double(tmp(5), spin, ierr)
      if (ierr /= 0) return
      call str_to_double(tmp(6), mass, ierr)
      if (ierr /= 0) return
      ref = tmp(7)
   
   end subroutine parse_line_nuclides_header
   
   subroutine parse_line_pfcn(line, pfcn, ierr)
      character(len = *), intent(in) :: line
      integer :: i
      integer, intent(inout) :: ierr
      integer, parameter :: num_cols = 8
      character(len = 256), dimension(num_cols) :: tmp
      real(dp), dimension(:), intent(out) :: pfcn
      
      call parse_line(line, num_cols, tmp)
      
      do i = 1, num_cols
         call str_to_double(tmp(i), pfcn(i), ierr)
         if (ierr /= 0) return
      end do
   
   end subroutine parse_line_pfcn
   
   
   subroutine parse_line(line, num_cols, line_out)
      character(len = *), intent(in) :: line
      integer, intent(in) :: num_cols
      character(len = 256), dimension(num_cols), intent(out) :: line_out
      character(len = 256) :: tmp
      integer :: k, i, j
      
      k = 1
      tmp = ''
      line_out = ''
      do j = 1, len(line)
         if(line(j:j)==' ') cycle
         tmp = trim(tmp) // line(j:j)
         if(j<len(line))then
            if(line(j + 1:j + 1)==' ') then
               line_out(k) = tmp(1:len_trim(tmp))
               k = k + 1
               tmp = ''
            end if
         end if
         if(k==num_cols + 1) exit
      end do
   
   end subroutine parse_line


end module chem_support
