! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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
! ***********************************************************************

module load_kap

  use kap_def
  use load_CO_kap, only: Setup_Kap_CO_Tables, min_version
  use kapcn, only: kapcn_init
  use kap_aesopus, only: kap_aesopus_init
  use const_def, only: dp, use_mesa_temp_cache
  use math_lib
  use utils_lib, only: mesa_error, mv, switch_str

  implicit none


  logical, parameter :: dbg = .false.

  logical, parameter :: dbg_cache = .false.

contains


  subroutine Setup_Kap_Tables(rq, &
       use_cache, load_on_demand, ierr)
    use const_def, only: mesa_data_dir
    use condint, only: init_potekhin
    type (Kap_General_Info), pointer :: rq
    logical, intent(in) :: use_cache, load_on_demand
    integer, intent(out) :: ierr

    integer :: id, iz, ix, num_Xs
    real(dp) :: X, Z

    include 'formats'

    ierr = 0

    kap_dir = trim(mesa_data_dir) // '/kap_data'

    kap_use_cache = use_cache
    
    call Setup_Kap_CO_Tables(rq, kap_CO_z_tables(rq% kap_CO_option)% ar, use_cache, load_on_demand, ierr)
    if (ierr /= 0) return
    
    call setup_lowT(kap_lowT_z_tables(rq% kap_lowT_option)% ar)
    call setup(kap_z_tables(rq% kap_option)% ar)
    
    rq% logT_Compton_blend_hi = kap_z_tables(rq% kap_option)% ar(1)% x_tables(1)% logT_max - 0.01d0
      !rq% kap_z_tables(1)% x_tables(1)% logT_max - 0.01d0
    rq% logR_Compton_blend_lo = kap_z_tables(rq% kap_option)% ar(1)% x_tables(1)% logR_min + 0.01d0
      !rq% kap_z_tables(1)% x_tables(1)% logR_min + 0.01d0

    call init_potekhin(ierr)

  contains

    subroutine setup(z_tables)
      type (Kap_Z_Table), dimension(:), pointer :: z_tables
      integer :: iz
      include 'formats'
      if (.not. associated(z_tables)) then
         allocate(z_tables(num_kap_Zs(rq% kap_option)), STAT=ierr)
         if (ierr /= 0) return
         do iz = 1, num_kap_Zs(rq% kap_option)
            z_tables(iz)% lowT_flag = .false.
            z_tables(iz)% Z = kap_Zs(iz, rq% kap_option)
            allocate(z_tables(iz)% x_tables(num_kap_Xs(rq% kap_option)), STAT=ierr)
            if (ierr /= 0) return
         end do
         if (show_allocations) write(*,2) 'setup z_tables', &
             num_kap_Zs(rq% kap_option) + num_kap_Zs(rq% kap_option)*num_kap_Xs(rq% kap_option)
      end if
      do iz = 1, num_kap_Zs(rq% kap_option)
         Z = kap_Zs(iz,rq% kap_option)
         num_Xs = num_kap_Xs_for_this_Z(iz,rq% kap_option)
         z_tables(iz)% num_Xs = num_Xs
         do ix = 1, num_Xs
            X = kap_Xs(ix,rq% kap_option)
            if (X + Z > 1) X = 1-Z
            call load_one(rq, z_tables, iz, ix, X, Z, &
                 load_on_demand, ierr)
            if (ierr /= 0) return
         end do
      end do

    end subroutine setup

    subroutine setup_lowT(z_tables)
      type (Kap_Z_Table), dimension(:), pointer :: z_tables
      integer :: iz
      include 'formats'
      select case (rq% kap_lowT_option)
      case (kap_lowT_kapCN)
         if (rq% show_info) write(*,*) 'Loading lowT kapCN tables'
         call kapCN_init(rq% show_info, ierr)
         if (ierr /= 0) then
            write(*,*) 'Failed to load kapCN opacities'
            call mesa_error(__FILE__, __LINE__)
         end if
      case (kap_lowT_AESOPUS)
         if (rq% show_info) write(*,*) 'Loading lowT AESOPUS tables'
         call kap_aesopus_init(rq, ierr)
         if (ierr /= 0)  then
            write(*,*) 'Failed to load AESOPUS opacities'
            call mesa_error(__FILE__, __LINE__)
         end if
      case default
         if (rq% kap_lowT_option > kap_lowT_options_max .or. &
             rq% kap_lowT_option <= 0) then
           write(*,*) 'invalid lowT opacity option'
           call mesa_error(__FILE__, __LINE__)
         end if
         if (rq% show_info) then
            if (rq% kap_lowT_option == kap_lowT_test) then
               write(*,*) 'Loading lowT test tables'
            else if (rq% kap_lowT_option == kap_lowT_Freedman11) then
               write(*,*) 'Loading lowT Freedman tables'
            else
               write(*,*) 'Loading lowT Fergusson tables'
            end if
         end if
         if (.not. associated(z_tables)) then
            allocate(z_tables(num_kap_lowT_Zs(rq% kap_lowT_option)), STAT=ierr)
            if (ierr /= 0) return
            do iz = 1, num_kap_lowT_Zs(rq% kap_lowT_option)
               z_tables(iz)% lowT_flag = .true.
               z_tables(iz)% Z = kap_lowT_Zs(iz, rq% kap_lowT_option)
               allocate(z_tables(iz)% x_tables(num_kap_lowT_Xs(rq% kap_lowT_option)), STAT=ierr)
               if (ierr /= 0) return
            end do
            if (show_allocations) write(*,2) 'setup_lowT z_tables', &
                num_kap_lowT_Zs(rq% kap_lowT_option) + num_kap_lowT_Zs(rq% kap_lowT_option)*num_kap_lowT_Xs(rq% kap_lowT_option)
         end if
         do iz = 1, num_kap_lowT_Zs(rq% kap_lowT_option)
            Z = kap_lowT_Zs(iz, rq% kap_lowT_option)
            num_Xs = num_kap_lowT_Xs_for_this_Z(iz, rq% kap_lowT_option)
            z_tables(iz)% num_Xs = num_Xs
            do ix = 1, num_Xs
               X = kap_lowT_Xs(ix, rq% kap_lowT_option)
               if (X + Z > 1) X = 1-Z
               call load_one(rq, z_tables, iz, ix, X, Z, &
                    load_on_demand, ierr)
               if (ierr /= 0) then
                  if (rq% kap_lowT_option == kap_lowT_test) then
                     write(*,*) 'Failed to load lowT test tables'
                  else if (rq% kap_lowT_option == kap_lowT_Freedman11) then
                     write(*,*) 'Failed to load lowT Freedman tables'
                  else
                     write(*,*) 'Failed to load lowT Fergusson tables'
                  end if
                  call mesa_error(__FILE__, __LINE__)
               end if
            end do
         end do
      end select

    end subroutine setup_lowT


  end subroutine Setup_Kap_Tables


  subroutine load_one(rq, &
       z_tables, iz, ix, X, Z, read_later, ierr)
    type (Kap_General_Info), pointer :: rq
    type (Kap_Z_Table), dimension(:), pointer :: z_tables
    integer, intent(in) :: iz, ix
    real(dp), intent(in) :: X, Z
    logical, intent(in) :: read_later
    integer, intent(out) :: ierr

    character (len=256) :: fname, filename, cache_filename, temp_cache_filename

    ierr = 0

    call Get_Filenames(rq, z_tables, X, Z, &
         kap_dir, fname, filename, cache_filename, temp_cache_filename, ierr)
    if (ierr /= 0) return
    
    if (rq% show_info) write(*,*) 'read filename <' // trim(filename) // '>'

    call Prepare_Kap_X_Table(rq, &
         z_tables, iz, z_tables(iz)% x_tables, ix, X, Z, &
         read_later, &
         fname, filename, cache_filename, temp_cache_filename,&
         ierr)
    if (ierr /= 0) return

    if (.not. read_later) z_tables(iz)% x_tables(ix)% not_loaded_yet = .false.

  end subroutine load_one


  subroutine Prepare_Kap_X_Table(rq, &
       z_tables, iz, x_tables, ix, X_in, Z_in, &
       read_later, &
       fname, filename, cache_filename, temp_cache_filename, &
       ierr)
    type (Kap_General_Info), pointer :: rq
    type (Kap_Z_Table), dimension(:), pointer :: z_tables
    ! allocates the arrays and stores data
    type (Kap_X_Table), dimension(:), pointer :: x_tables
    integer, intent(in) :: iz, ix
    real(dp), intent(in) :: X_in, Z_in
    logical, intent(in) :: read_later
    character (*), intent(in) :: fname, filename, cache_filename, temp_cache_filename
    integer, intent(out) :: ierr

    integer :: io_unit, cache_io_unit

    real(dp) :: X, Z
    integer :: ios, form, version, n, num_logRs, num_logTs, status, nvec, i
    real(dp) :: xin, zz, logR_min, logR_max, logT_min, logT_max, xErr, zErr
    character (len=1000) :: message
    real(dp), target :: vec_ary(50)
    real(dp), pointer :: vec(:)
    real(dp), parameter :: tiny = 1d-6
    include 'formats'

    ierr = 0
    vec => vec_ary
    X = X_in
    Z = Z_in
    nvec=-1

    open(newunit=io_unit,file=trim(filename),action='read',status='old',iostat=ios)
    if (ios /= 0) then
       !write(*,*) 'load kap tables ' // trim(filename)
       write(*,'(A)')
       write(*,'(A)')
       write(*,'(A)')
       write(*,'(A)')
       write(*,'(A)')
       write(*,*) 'NOTICE: missing kap data ' // trim(filename)
       write(*,*) 
       write(*,*) 'Please check the validity of the kap_prefix string for this file.'
       write(*,*) 
       write(*,*) 'If it is okay, you may need to install new kap data.'
       write(*,*) 'To do that, remove the directory mesa/data/kap_data,'
       write(*,*) 'and rerun the mesa ./install script.'
       write(*,'(A)')
       call mesa_error(__FILE__,__LINE__)
    end if

    version = -1

    read(io_unit,*,iostat=ierr) ! skip
    if (ierr == 0) then
       read(io_unit,*,iostat=ierr) ! skip
       if (ierr == 0) then
          read(io_unit,'(a)',iostat=ierr) message
          if (ierr == 0) call str_to_vector(message, vec, nvec, ierr)
          if (nvec < 10) ierr = -1
          if (ierr == 0) then
             i = 0
             i = i+1; form = int(vec(i))
             i = i+1; version = int(vec(i))
             i = i+1; xin = vec(i)
             i = i+1; zz = vec(i)
             i = i+1; num_logRs = int(vec(i))
             i = i+1; logR_min = vec(i)
             i = i+1; logR_max = vec(i)
             i = i+1; num_logTs = int(vec(i))
             i = i+1; logT_min = vec(i)
             i = i+1; logT_max = vec(i)
          end if
       end if
    end if



    if (ierr /= 0 .or. version < min_version) then
       write(*,*) 'load kap tables ' // trim(filename)
       write(*,'(A)')
       write(*,'(A)')
       write(*,'(A)')
       write(*,'(A)')
       write(*,'(A)')
       write(*,*) 'NOTICE: you need to install a new verion of the kap data.'
       write(*,*) 'Please remove the directory mesa/data/kap_data,'
       write(*,*) 'and rerun the mesa ./install script.'
       write(*,'(A)')
       if (ierr /= 0) write(*,*) 'ierr', ierr
       if (version < min_version) &
            write(*,*) 'version < min_version', version, min_version
       write(*,*) 'form', form
       write(*,*) 'xin', xin
       write(*,*) 'zz', zz
       write(*,*) 'num_logRs', num_logRs
       write(*,*) 'logR_min', logR_min
       write(*,*) 'logR_max', logR_max
       write(*,*) 'num_logTs', num_logTs
       write(*,*) 'logT_min', logT_min
       write(*,*) 'logT_max', logT_max
       call mesa_error(__FILE__,__LINE__)
    end if

    if (form /= kap_table_fixed_metal_form) then
       call mesa_error(__FILE__,__LINE__,'form /= kap_table_fixed_metal_form')
    end if

    call Setup_Kap_X_Table(ierr)
    if (ierr /= 0) return

    if (read_later) then
       close(io_unit)
       return
    end if

    if (kap_use_cache) then
       ios = 0
       if (dbg_cache) then
          open(newunit=cache_io_unit,file=trim(cache_filename),action='read', &
               status='old',iostat=ios)
       else
          open(newunit=cache_io_unit,file=trim(cache_filename),action='read', &
               status='old',iostat=ios,form='unformatted')
       end if
       if (ios == 0) then ! try reading the cached data
          !write(*,'(a)') 'loading ' // trim(cache_filename)
          call Read_Kap_X_Table(cache_io_unit, .true., ierr)
          close(cache_io_unit)
          if (ierr == 0) then
             close(io_unit)
             return
          end if
          ierr = 0
       else
          !write(*,*) 'failed to open ' // trim(cache_filename)
       end if
    end if

    if (show_allocations) write(*,2) 'loading ' // trim(filename)
    call Read_Kap_X_Table(io_unit, .false., ierr)
    close(io_unit)
    if (ierr /= 0) then
       write(*,*) 'failed in Read_Kap_X_Table ' // trim(filename)
       call mesa_error(__FILE__,__LINE__)
       return
    end if

    if (.not. kap_use_cache) return

    if (dbg_cache) then
       open(newunit=cache_io_unit, file=trim(switch_str(temp_cache_filename, cache_filename, use_mesa_temp_cache)),&
            iostat=ios, action='write')
    else
       open(newunit=cache_io_unit, file=trim(switch_str(temp_cache_filename, cache_filename, use_mesa_temp_cache)),&
            iostat=ios, action='write', form='unformatted')
    end if

    if (ios == 0) then
       write(*,'(a)') 'write ' // trim(cache_filename)
       call Write_Kap_X_Table_Cache( &
            z_tables(iz)% x_tables, ix, cache_io_unit,  version)
       close(cache_io_unit)
       ! Atomic move temp cache file to cache folder
       if(use_mesa_temp_cache) call mv(temp_cache_filename,cache_filename,.true.)
       if (kap_read_after_write_cache) then
          open(newunit=cache_io_unit,file=trim(cache_filename),action='read', &
               status='old',iostat=ios,form='unformatted')
          if (ios == 0) then
             call Read_Kap_X_Table(cache_io_unit, .true., ierr)
             close(cache_io_unit)
          end if
       end if
    end if


  contains


    subroutine Setup_Kap_X_Table(ierr)
      integer, intent(out) :: ierr

      integer :: i

      xErr = abs(xin - X); zErr = abs(zz - Z)
      if (xErr > tiny .or. zErr > tiny) then
         ierr = -1
         write(*,*) 'bug in file ' // trim(filename), xErr, zErr
         write(*,*) 'xErr > tiny', xErr > tiny
         write(*,*) 'zErr > tiny', zErr > tiny
         write(*,*) 'xin', xin
         write(*,*) 'X', X
         write(*,*) 'zz', zz
         write(*,*) 'Z', Z
         return
      end if

      x_tables(ix)% not_loaded_yet = .true.            
      x_tables(ix)% X = X
      x_tables(ix)% Z = Z
      x_tables(ix)% logR_min = logR_min
      x_tables(ix)% logR_max = logR_max

      x_tables(ix)% num_logRs = num_logRs
      nullify(x_tables(ix)% logRs)

      x_tables(ix)% logT_min = logT_min
      x_tables(ix)% logT_max = logT_max
      x_tables(ix)% num_logTs = num_logTs
      nullify(x_tables(ix)% logTs)

      nullify(x_tables(ix)% kap1) ! allocate when read the data

    end subroutine Setup_Kap_X_Table


    subroutine Read_Kap_X_Table(io_unit, reading_cache, ierr)
      integer, intent(in) :: io_unit ! use this for file access
      logical, intent(in) :: reading_cache
      integer, intent(out) :: ierr

      character (len=1000) :: message
      character (len=1) :: char
      integer :: i, j, c_version, c_num_logRs, c_num_logTs
      real(dp) :: c_xin, c_zz, c_logR_min, c_logR_max, &
           c_logT_min, c_logT_max
      real(dp) :: kap_logKs(num_logRs), logT
      real(dp), pointer :: kap(:,:,:), kap1(:)
      real(dp), target :: vec_ary(100)
      real(dp), pointer :: vec(:)
      integer :: nvec

      include 'formats'

      vec => vec_ary
      nvec=-1
      
      if (reading_cache) then

         ios = 0
         if (dbg_cache) then
            write(*,*) 'io_unit', io_unit
            read(io_unit, *, iostat=ios) c_version, c_num_logRs, c_num_logTs
         else
            read(io_unit, iostat=ios) c_version, c_num_logRs, c_num_logTs
         end if
         if (ios /= 0 .or. c_version /= version) then
            ierr = 1
            if (ios /= 0) then
               write(*,*) 'cache failed in read 1'
            else if (c_version /= version) then
               write(*,*) 'cache failed for c_version /= version'
               write(*,*) 'c_version', c_version
               write(*,*) 'version', version
            else if (c_num_logRs /= num_logRs .or. c_num_logTs /= num_logTs) then
               if (c_num_logRs /= num_logRs) write(*,*) 'cache failed for c_num_logRs /= num_logRs'
               if (c_num_logTs /= num_logTs) write(*,*) 'cache failed for c_num_logTs /= num_logTs'
            end if
            return
         end if

         read(io_unit, iostat=ios) &
              c_xin, c_zz, c_logR_min, c_logR_max, &
              c_logT_min, c_logT_max
         if (ios /= 0) then
            ierr = 1
            if (ios /= 0) write(*,*) 'cache failed in read 2'
         end if

      end if

      xErr = abs(xin - X); zErr = abs(zz - Z)
      if (xErr > tiny .or. zErr > tiny) then
         if (reading_cache) then
            if (xErr > tiny) write(*,*) 'cache failed for xErr > tiny'
            if (zErr > tiny) write(*,*) 'cache failed for zErr > tiny'
            ierr = 1; return
         end if
         ierr = -1
         return
      end if
      
      if (show_allocations) write(*,2) 'x_tables ' // trim(filename), &
        num_logRs + num_logTs + sz_per_Kap_point*num_logRs*num_logTs
      allocate(x_tables(ix)% logRs(num_logRs), x_tables(ix)% logTs(num_logTs), &
           x_tables(ix)% kap1(sz_per_Kap_point*num_logRs*num_logTs), STAT=status)
      if (status .ne. 0) then
         ierr = -1
         return
      end if

      kap1 => x_tables(ix)% kap1
      kap(1:sz_per_Kap_point,1:num_logRs,1:num_logTs) => &
           kap1(1:sz_per_Kap_point*num_logRs*num_logTs)

      if (.not. reading_cache) then

         read(io_unit,*,iostat=ierr) ! skip
         if (ierr /= 0) return
         read(io_unit,*,iostat=ierr) ! skip
         if (ierr /= 0) return

         read(io_unit,'(a)',iostat=ierr) message
         if (ierr == 0) call str_to_vector(message, vec, nvec, ierr)
         if (nvec < num_logRs) ierr = -1
         if (ierr /= 0) return
         do j=1,num_logRs
            x_tables(ix)% logRs(j) = vec(j)
            !write(*,*) 'read logR', j, vec(j)
         end do
         !write(*,*) "input line: <" // trim(message) // ">"

         read(io_unit,*,iostat=ierr) ! skip
         if (ierr /= 0) return

         do i = 1, num_logTs
            read(io_unit,'(a)',iostat=ierr) message
            if (ierr /= 0) then
               write(*,*) 'failed in read for logT i', i
               return
               !stop
            end if
            call str_to_vector(message, vec, nvec, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in str_to_vector'
               write(*,*) 'nvec', nvec
               write(*,'(a)') 'message "' // trim(message) // '"'
               return
               !stop
            end if

            if (nvec < 1+num_logRs) ierr = -1
            if (ierr /= 0) return
            x_tables(ix)% logTs(i) = vec(1)
            do j=1,num_logRs
               kap_logKs(j) = vec(j+1)
               kap(1,j,i) = kap_logKs(j)
            end do

            if (.false.) then
               write(*,2) 'logT', i, x_tables(ix)% logTs(i)
               do j=1,num_logRs
                  write(*,3) 'kap', j, i, kap(1,j,i)
               end do
               write(*,'(a)') 'message "' // trim(message) // '"'
            end if

         end do

         call Make_Interpolation_Data( &
              kap1, x_tables(ix)% logRs, num_logRs, &
              x_tables(ix)% logTs, num_logTs, &
              x_tables(ix)% ili_logRs, x_tables(ix)% ili_logTs, ierr)

         if (ierr /= 0) write(*,*) 'Read_Kap_X_Table failed in Make_Interpolation_Data'

      else ! reading_cache

         read(io_unit, iostat=ierr) &
              x_tables(ix)% ili_logRs, x_tables(ix)% ili_logTs
         if (ierr /= 0) return

         read(io_unit, iostat=ierr) &
              x_tables(ix)% logRs(1:num_logRs), &
              x_tables(ix)% logTs(1:num_logTs)            
         if (ierr /= 0) return

         read(io_unit, iostat=ierr) kap1
         if (ierr /= 0) return

      end if

    end subroutine Read_Kap_X_Table


  end subroutine Prepare_Kap_X_Table


  subroutine Make_Interpolation_Data( &
       kap1, logRs, num_logRs, logTs, num_logTs, ili_logRs, ili_logTs, ierr)
    use interp_2d_lib_db
    real(dp), pointer :: kap1(:)
    integer, intent(in) :: num_logRs, num_logTs
    real(dp), intent(in), pointer :: logRs(:) ! (num_logRs)
    real(dp), intent(in), pointer :: logTs(:) ! (num_logTs)
    integer, intent(out) :: ili_logRs, ili_logTs, ierr

    character (len=256) :: message
    integer :: ibcxmin                   ! bc flag for x=xmin
    real(dp) :: bcxmin(num_logTs)               ! bc data vs. y at x=xmin
    integer :: ibcxmax                   ! bc flag for x=xmax
    real(dp) :: bcxmax(num_logTs)               ! bc data vs. y at x=xmax
    integer :: ibcymin                   ! bc flag for y=ymin
    real(dp) :: bcymin(num_logRs)               ! bc data vs. x at y=ymin
    integer :: ibcymax                   ! bc flag for y=ymax
    real(dp) :: bcymax(num_logRs)               ! bc data vs. x at y=ymax
    integer :: ier                       ! =0 on exit if there is no error.
    integer :: i, j
    real(dp) :: tol
    real(dp), pointer :: kap(:,:,:)
    real(dp), pointer :: x_out(:), y_out(:), f_out(:,:)

    kap(1:sz_per_kap_point,1:num_logRs,1:num_logTs) => &
         kap1(1:sz_per_kap_point*num_logRs*num_logTs)

    ! just use "not a knot" bc's at edges of tables
    ibcxmin = 0; bcxmin(1:num_logTs) = 0d0
    ibcxmax = 0; bcxmax(1:num_logTs) = 0d0
    ibcymin = 0; bcymin(1:num_logRs) = 0d0
    ibcymax = 0; bcymax(1:num_logRs) = 0d0

    call interp_mkbicub_db( &
         logRs, num_logRs, logTs, num_logTs, kap1, num_logRs, &
         ibcxmin,bcxmin,ibcxmax,bcxmax, &
         ibcymin,bcymin,ibcymax,bcymax, &
         ili_logRs,ili_logTs,ier)

    if (ier /= 0) then
       write(*,*) 'interp_mkbicub_db error happened for Make_Interpolation_Data for table'
       write(*,*) 'num_logRs', num_logRs
       do i=1,num_logRs
          write(*,*) 'logR', i, logRs(i)
       end do
       write(*,*) 'num_logTs', num_logTs
       do i=1,num_logTs
          write(*,*) 'logT', i, logTs(i)
       end do
       write(*,'(A)')
       call mesa_error(__FILE__,__LINE__,'kap interp error')
       ierr = -1
       return
    end if

    call Check_Interpolation_Data

    ierr = 0

  contains

    subroutine Check_Interpolation_Data
      use utils_lib,only:is_bad
      integer :: i, iR, jtemp
      real(dp) :: val

      do i = 1, sz_per_kap_point
         do iR = 1, num_logRs
            do jtemp = 1, num_logTs
               val = kap(i,iR,jtemp)
               if (is_bad(val)) then
                  if (.true.) then
                     write(*,*) 'bad value in xz', val, i, iR, jtemp
                     write(*,'(99(a15,3x,f15.8,3x))')  &
                          'logR', logRs(iR), 'logT', logTs(jtemp)
                  end if
                  kap(i,iR,jtemp) = 0
               end if
            end do
         end do
      end do

    end subroutine Check_Interpolation_Data


  end subroutine Make_Interpolation_Data


  subroutine Write_Kap_X_Table_Cache(x_tables, ix, io_unit,  version)
    type (Kap_X_Table), dimension(:), pointer :: x_tables
    integer, intent(in) :: ix, io_unit, version


    if (dbg_cache) then
       write(*,*) 'write cache plain text', io_unit
       write(io_unit,*) version, x_tables(ix)% num_logTs, x_tables(ix)% num_logRs
       write(io_unit,'(999(1pe26.16))') &
            x_tables(ix)% X, &
            x_tables(ix)% Z, &
            x_tables(ix)% logR_min, &
            x_tables(ix)% logR_max, &
            x_tables(ix)% logT_min, &
            x_tables(ix)% logT_max
       !write(io_unit) &
       !   x_tables(ix)% logRs(:), &
       !   x_tables(ix)% logTs(:)                            
       !write(io_unit) x_tables(ix)% kap(:,:,:)   
    else
       write(io_unit) version, x_tables(ix)% num_logTs, x_tables(ix)% num_logRs
       write(io_unit) &
            x_tables(ix)% X, &
            x_tables(ix)% Z, &
            x_tables(ix)% logR_min, &
            x_tables(ix)% logR_max, &
            x_tables(ix)% logT_min, &
            x_tables(ix)% logT_max                           
       write(io_unit) x_tables(ix)% ili_logRs, x_tables(ix)% ili_logTs 
       write(io_unit) &
            x_tables(ix)% logRs(:), &
            x_tables(ix)% logTs(:)                             
       write(io_unit) x_tables(ix)% kap1(:) 
    end if

  end subroutine Write_Kap_X_Table_Cache


  subroutine Get_Filenames(rq, &
       z_tables, X, Z, data_dir, &
       fname, filename, cache_filename, temp_cache_filename, ierr)
    type (Kap_General_Info), pointer :: rq
    type (Kap_Z_Table), dimension(:), pointer :: z_tables
    real(dp), intent(in) :: X, Z
    character (*), intent(in) :: data_dir
    character (*), intent(out) :: fname, filename, cache_filename, temp_cache_filename
    integer, intent(out) :: ierr
    character (len=256) :: cache_fname
    ierr=0
    call create_fname(rq, z_tables, X, Z, fname, cache_fname,ierr)
    filename = trim(data_dir) // '/' // fname
    cache_filename = trim(kap_cache_dir) // '/' // cache_fname
    temp_cache_filename=trim(kap_temp_cache_dir) // '/' // cache_fname
  end subroutine Get_Filenames


  ! this must match the preprocessor naming scheme
  subroutine create_fname(rq, z_tables, X, Z, fname, cache_fname, ierr)
    type (Kap_General_Info), pointer :: rq
    type (Kap_Z_Table), dimension(:), pointer :: z_tables
    real(dp), intent(in) :: X, Z
    character (len=*), intent(out) :: fname, cache_fname
    integer, intent(out) :: ierr
    character (len=256) :: zstr, xstr, prefix

    ierr=0
    
    if (z_tables(1)% lowT_flag .and. rq% kap_lowT_option == kap_lowT_Freedman11) then
       call get_output_string(Z,zstr,ierr)
       fname = trim(kap_lowT_option_str(rq% kap_lowT_option)) // '_z' // trim(zstr) // '.data'
       cache_fname = trim(kap_lowT_option_str(rq% kap_lowT_option))// '_z' // trim(zstr) // '.bin'
       return
    end if

    call get_output_string(Z,zstr,ierr)
    call get_output_string(X,xstr,ierr)
    if (z_tables(1)% lowT_flag) then
       prefix = kap_lowT_option_str(rq% kap_lowT_option)
    else
       prefix = kap_option_str(rq% kap_option)
    end if

    fname = &
         trim(prefix) // '_z' // trim(zstr) // '_x' // trim(xstr) // '.data'
    cache_fname = &
         trim(prefix) // '_z' // trim(zstr) // '_x' // trim(xstr) // '.bin'
  end subroutine create_fname


end module load_kap

