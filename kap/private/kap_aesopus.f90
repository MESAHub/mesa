! ***********************************************************************
!
!   Copyright (C) 2018-2020  Aaron Dotter, Josiah Schwab & The MESA Team
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

module kap_aesopus

  use hdf5
  use iso_c_binding

  use math_lib
  use kap_def

  implicit none

  logical, parameter :: debug = .false.

  ! for 2-D interpolation
  integer :: ibcx=0, ibcy=0, ilinx=1, iliny=1
  integer :: ict(6) = [ 1, 1, 1, 0, 0, 0 ]
  real(dp), parameter :: bc(100) = 0d0

contains

  subroutine kap_aesopus_init(rq, ierr)
    use const_def, only: mesa_data_dir
    use kap_def, only: kap_aesopus_is_initialized
    type (Kap_General_Info), pointer :: rq
    integer, intent(out) :: ierr

    ! real(dp) :: kap, dlnkap_dlnRho, dlnkap_dlnT

    call read_kap_aesopus_tables(rq, ierr)

    ! call AESOPUS_get(0.02d0, 0.7d0, 0d0, 0d0, 0d0, -10d0, 4d0, &
    !      kap, dlnkap_dlnRho, dlnkap_dlnT,ierr)

    ! write(*,*) kap, dlnkap_dlnRho, dlnkap_dlnT

    if (ierr == 0) kap_aesopus_is_initialized = .true.

  end subroutine kap_aesopus_init


  subroutine kap_aesopus_shutdown
    use kap_def, only: kap_aesopus_is_initialized
    kap_aesopus_is_initialized = .false.
  end subroutine kap_aesopus_shutdown


  subroutine read_kap_aesopus_tables(rq, ierr)

    use interp_2d_lib_db, only: interp_mkbicub_db

    type (Kap_General_Info), pointer :: rq
    integer, intent(out) :: ierr

    character(len=256) :: filename
    integer :: n, io
    integer :: iX, iCO, iC, iN

    integer(hid_t) :: file_id, subgroup_id, group_id, dspace_id, dset_id
    integer(hsize_t), dimension(6) :: data_dims, max_dims

    real(dp), dimension(:,:,:,:,:,:), allocatable :: table_data
    integer :: table_size

    integer(hsize_t) :: idx ! index
    integer(size_t) :: size ! size of group_name
    character(len=80) :: group_name ! output buffer

    character(len=30) :: efmt = '(A14, 99ES12.3)'
    character(len=30) :: ffmt = '(A14, 99F8.3)'
    character(len=30) :: ifmt = '(A14, I4)'

    logical :: file_exists 

    ! get the filename
    filename = trim(aesopus_filename)
    if (filename == '') then
       write(*,*) 'failed to specify AESOPUS_filename'
       ierr = -1
       return
    end if

    ! first try local directory
    inquire(file=trim(filename), exist=file_exists)
    if (.not. file_exists) then
       ! then try MESA directory
       filename = trim(kap_dir) // '/' // filename
       inquire(file=trim(filename), exist=file_exists)
       if (.not. file_exists) then
          write(*,*) 'failed to open AESOPUS file ' // trim(filename)
          ierr = -1
          return
       end if
    end if

    if (rq% show_info) write(*,*) 'read filename <' // trim(filename) // '>'

    ierr = 0

    ! open hdf5 interface
    call h5open_f(ierr)
    if (ierr /= 0) return

    ! open file (read-only)
    call h5fopen_f(filename, h5f_acc_rdonly_f, file_id, ierr)
    if (ierr /= 0) return

    ! open root group
    call h5gopen_f(file_id, "/", group_id, ierr)
    if (ierr /= 0) return

    if (rq% show_info) write(*,*) 'AESOPUS composition parameters'

    ! read composition parameters
    data_dims = 0
    call h5dopen_f(group_id, "Zsun", dset_id, ierr)
    call h5dread_f(dset_id, H5T_IEEE_F64LE, kA% Zsun, data_dims, ierr)
    call h5dclose_f(dset_id, ierr)
    if (rq% show_info) write(*,efmt) 'Zsun =', kA % Zsun

    call h5dopen_f(group_id, "fCO_ref", dset_id, ierr)
    call h5dread_f(dset_id, H5T_IEEE_F64LE, kA% fCO_ref, data_dims, ierr)
    call h5dclose_f(dset_id, ierr)
    if (rq% show_info) write(*,ffmt) 'fCO_ref =', kA % fCO_ref

    call h5dopen_f(group_id, "fC_ref", dset_id, ierr)
    call h5dread_f(dset_id, H5T_IEEE_F64LE, kA% fC_ref, data_dims, ierr)
    call h5dclose_f(dset_id, ierr)
    if (rq% show_info) write(*,ffmt) 'fC_ref =', kA % fC_ref

    call h5dopen_f(group_id, "fN_ref", dset_id, ierr)
    call h5dread_f(dset_id, H5T_IEEE_F64LE, kA% fN_ref, data_dims, ierr)
    call h5dclose_f(dset_id, ierr)
    if (rq% show_info) write(*,ffmt) 'fN_ref =', kA % fN_ref


    if (rq% show_info) then
       write(*,*)
       write(*,*) 'AESOPUS logT and logR range (logR = logRho - 3 * logT + 18)'
    end if

    ! read logT
    data_dims = 0
    call h5dopen_f(group_id, "logTs", dset_id, ierr)
    call h5dget_space_f(dset_id, dspace_id, ierr)
    call H5sget_simple_extent_dims_f(dspace_id, data_dims, max_dims, ierr)
    if (rq% show_info) write(*,ifmt) "num logTs =", data_dims(1)
    kA % num_logTs = data_dims(1)
    allocate(kA% logTs(kA % num_logTs))
    call h5dread_f(dset_id, H5T_IEEE_F64LE, kA% logTs, data_dims, ierr)
    call h5dclose_f(dset_id, ierr)

    kA% min_logT = minval(kA% logTs)
    kA% max_logT = maxval(kA% logTs)

    if (rq% show_info) then
       write(*,ffmt) 'logTs =', kA% logTs
    end if

    ! read logR
    data_dims = 0
    call h5dopen_f(group_id, "logRs", dset_id, ierr)
    call h5dget_space_f(dset_id, dspace_id, ierr)
    call H5sget_simple_extent_dims_f(dspace_id, data_dims, max_dims, ierr)
    if (rq% show_info) write(*,ifmt) "num logRs =", data_dims(1)
    kA % num_logRs = data_dims(1)
    allocate(kA% logRs(kA % num_logRs))
    call h5dread_f(dset_id, H5T_IEEE_F64LE, kA% logRs, data_dims, ierr)
    call h5dclose_f(dset_id, ierr)

    kA% min_logR = minval(kA% logRs)
    kA% max_logR = maxval(kA% logRs)

    if (rq% show_info) then
       write(*,ffmt) 'logRs =', kA% logRs
    end if

    if (rq% show_info) then
       write(*,*)
       write(*,*) 'AESOPUS metallicities'
       write(*,*) '(These Z values are the reference metallicities)'
    end if

    ! read Zs
    data_dims = 0
    call h5dopen_f(group_id, "Zs", dset_id, ierr)
    call h5dget_space_f(dset_id, dspace_id, ierr)
    call H5sget_simple_extent_dims_f(dspace_id, data_dims, max_dims, ierr)
    if (rq% show_info) write(*,ifmt) "num Zs =", data_dims(1)
    kA % num_Zs = data_dims(1)
    allocate(kA% Zs(kA % num_Zs))
    call h5dread_f(dset_id, H5T_IEEE_F64LE, kA% Zs, data_dims, ierr)
    call h5dclose_f(dset_id, ierr)
    if (debug) write(*,*) 'Zs', kA% Zs

    if (rq% show_info) then
       write(*,efmt) "Zs =", kA% Zs
    end if


    ! pre-compute logZs
    allocate(kA% logZs(kA % num_Zs))
    do n = 1, kA % num_Zs
       kA% logZs(n) = safe_log10(kA% Zs(n))
    end do

    ! now, there is one group for each Z
    ! walk through these and construct a AESOPUS_TableSet for each one

    allocate(kA% ts(kA% num_Zs))

    do n = 1, kA% num_Zs

       associate(ts => kA% ts(n))

       ! get group name and open group
       idx = n - 1
       call h5lget_name_by_idx_f(file_id, ".", H5_INDEX_NAME_F, H5_ITER_INC_F, idx, group_name, ierr, size)
       if (rq% show_info) then
          write(*,*)
          write(*,'(3A, ES9.3)') "Table ",  trim(group_name), ": Z = ", kA% Zs(n)
       end if
       call h5gopen_f(group_id, group_name, subgroup_id, ierr)

       ! read Xs
       data_dims = 0
       call h5dopen_f(subgroup_id, "Xs", dset_id, ierr)
       call h5dget_space_f(dset_id, dspace_id, ierr)
       call H5sget_simple_extent_dims_f(dspace_id, data_dims, max_dims, ierr)
       if (rq% show_info) write(*,ifmt) "num Xs =", data_dims(1)
       ts % num_Xs = data_dims(1)
       allocate(ts% Xs(ts % num_Xs))
       call h5dread_f(dset_id, H5T_IEEE_F64LE, ts% Xs, data_dims, ierr)
       call h5dclose_f(dset_id, ierr)
       if (debug) write(*,*) "Xs", ts% Xs

       if (rq% show_info) then
          write(*,ffmt) "Xs =", ts% Xs
       end if


       ! read fCOs
       data_dims = 0
       call h5dopen_f(subgroup_id, "fCOs", dset_id, ierr)
       call h5dget_space_f(dset_id, dspace_id, ierr)
       call H5sget_simple_extent_dims_f(dspace_id, data_dims, max_dims, ierr)
       if (rq% show_info) write(*,ifmt) "num fCOs =", data_dims(1)
       ts % num_fCOs = data_dims(1)
       allocate(ts% fCOs(ts % num_fCOs))
       call h5dread_f(dset_id, H5T_IEEE_F64LE, ts% fCOs, data_dims, ierr)
       call h5dclose_f(dset_id, ierr)
       if (debug) write(*,*) "fCOs", ts% fCOs

       if (rq% show_info) then
          write(*,ffmt) "fCOs =", ts% fCOs
       end if


       ! read fCs
       data_dims = 0
       call h5dopen_f(subgroup_id, "fCs", dset_id, ierr)
       call h5dget_space_f(dset_id, dspace_id, ierr)
       call H5sget_simple_extent_dims_f(dspace_id, data_dims, max_dims, ierr)
       if (rq% show_info) write(*,ifmt) "num fCs =", data_dims(1)
       ts % num_fCs = data_dims(1)
       allocate(ts% fCs(ts % num_fCs))
       call h5dread_f(dset_id, H5T_IEEE_F64LE, ts% fCs, data_dims, ierr)
       call h5dclose_f(dset_id, ierr)
       if (debug) write(*,*) "fCs", ts% fCs

       if (rq% show_info) then
          write(*,ffmt) "fCs =", ts% fCs
       end if

       ! read fNs
       data_dims = 0
       call h5dopen_f(subgroup_id, "fNs", dset_id, ierr)
       call h5dget_space_f(dset_id, dspace_id, ierr)
       call H5sget_simple_extent_dims_f(dspace_id, data_dims, max_dims, ierr)
       if (rq% show_info) write(*,ifmt) "num fNs =", data_dims(1)
       ts % num_fNs = data_dims(1)
       allocate(ts% fNs(ts % num_fNs))
       call h5dread_f(dset_id, H5T_IEEE_F64LE, ts% fNs, data_dims, ierr)
       call h5dclose_f(dset_id, ierr)
       if (debug) write(*,*) "fNs", ts% fNs

       if (rq% show_info) then
          write(*,ffmt) "fNs =", ts% fNs
       end if


       ! read opacities
       data_dims = 0
       call h5dopen_f(subgroup_id, "kap", dset_id, ierr)
       call h5dget_space_f(dset_id, dspace_id, ierr)
       call H5sget_simple_extent_dims_f(dspace_id, data_dims, max_dims, ierr)
       if (debug) write(*,*) "data_dims", data_dims
       allocate(table_data(data_dims(1), data_dims(2), data_dims(3), &
                           data_dims(4), data_dims(5), data_dims(6)))
       call h5dread_f(dset_id, H5T_IEEE_F64LE, table_data, data_dims, ierr)
       call h5dclose_f(dset_id, ierr)


       allocate(ts% t(ts % num_Xs, ts % num_fCOs, ts % num_fCs, ts % num_fNs))

       do iX = 1, ts % num_Xs
          do iCO = 1, ts % num_fCOs
             do iC = 1, ts % num_fCs
                do iN = 1, ts % num_fNs

                   associate(t=> ts% t(iX, iCO, iC, iN))

                     !allocate(t% kap(4, kA% num_logRs, kA% num_logTs))

                     table_size = kA% num_logRs*kA% num_logTs
                     allocate(t% kap(4*table_size))

                     ! insert data
                     t% kap(1:4*table_size:4) = reshape(table_data(iN, iC, iCO, iX, :, :), [table_size])
                     t% X = ts% Xs(iX)
                     t% Z = kA% Zs(n)
                     t% fCO = ts% fCOs(iCO)
                     t% fC = ts% fCs(iC)
                     t% fN = ts% fNs(iN)

                     ! create interpolant
                     call interp_mkbicub_db(kA% logRs, kA% num_logRs, &
                          ka% logTs, kA% num_logTs, &
                          t% kap, ka% num_logRs, &
                          ibcx, bc, ibcx, bc, &
                          ibcy, bc, ibcy, bc, &
                          ilinx, iliny, ierr)

                     ! if (debug) call interp_evbicub_db(0d0, 3.85d0, &
                     !      kA% logRs, kA% num_logRs, &
                     !      ka% logTs, kA% num_logTs, &
                     !      ilinx, iliny, &
                     !      t% kap, ka% num_logRs, &
                     !      ict, res, ierr)
                     ! if (debug) write(*,'(10F10.4)') kA% logRs(17), kA% logTs(54), &
                     !      t% X, t% fCO, t% fC, t% fN, &
                     !      table_data(iN, iC, iCO, iX, 17, 54), res(1)

                   end associate

                end do
             end do
          end do
       end do

       deallocate(table_data)

       call h5gclose_f(subgroup_id, ierr)

       end associate

    end do

    ! close file
    call h5fclose_f(file_id, ierr)
    if (ierr /= 0) return

    ! close interface
    call h5close_f(ierr)
    if (ierr /= 0) return


    if (rq% show_info) then
       write(*,*)
       write(*,*) 'Finished reading AESOPUS tables'
       write(*,*)
    end if


  end subroutine read_kap_aesopus_tables


  subroutine AESOPUS_interp(Zref, X, XC, XN, XO, logR, logT, res, ierr)
    !result=(logKappa, dlogKappa/dlogR, dlogKappa/dlogT)
    use interp_1d_def, only: pm_work_size
    use interp_1d_lib, only: interpolate_vector, interp_pm
    use num_lib, only: binary_search
    real(dp), intent(in) :: Zref, X, XC, XN, XO, logR, logT
    real(dp), intent(out) :: res(3)
    integer, intent(out) :: ierr
    real(dp), pointer :: Z_ary(:), work1(:)
    integer, parameter :: nZ = 3
    integer :: i,iZ
    real(dp) :: my_Z, raw_res(3, nZ), x_new(1), v_new(1)
    character(len=32) :: junk
    logical :: clipped

    ! result=0d0; iZ=0; ierr=0

    if (outside_R_and_T_bounds(logR,logT)) then
       write(*,*) 'AESOPUS_interp: logR, logT outside of table bounds'
       ierr = -1
       return
    endif

    ! restrict to range
    clipped = .false.
    if (Zref .le. kA% Zs(1)) then
       my_Z = kA% Zs(1)
       iZ = 1
       clipped = .true.
    else if (Zref .ge. kA% Zs(kA% num_Zs)) then
       my_Z = kA% Zs(kA% num_Zs)
       iZ = kA% num_Zs
       clipped = .true.
    endif

    ! if clipped in Z, then just need one call
    if (clipped) then
       call AESOPUS_interp_fixedZref(iZ, X, XC, XN, XO, logR, logT, res, ierr)
       return
    endif

    my_Z = Zref

    ! it might be easier just to do linear interpolation
    ! but for now, we use an adapted version of what kapCN does

    ! require at least 3 Zs for interpolation
    if (nZ .gt. kA% num_Zs) then
       write(*,*) 'AESOPUS_interp: insufficient number of Z values for interpolation'
       write(*,'(I2, A, I2, A)') nZ, ' values are required; ', kA% num_Zs, ' were provided'
       ierr = -1
       return
    endif

    ! binary_search returns iZ between 1 and num_Zs-1
    ! such that Zs(iZ) <= Z < Zs(iZ+1)
    iZ = binary_search(kA% num_Zs, kA% Zs, 0, my_Z)

    ! make sure this is in an acceptable range
    ! unless that would go off the ends
    ! this check is hard-coded assuming nZ = 3
    iZ = min(kA% num_Zs-1, max(iZ, 2))

    ! want to call at [iZ-1, iZ, iZ+1]
    do i = 1, nZ
       call AESOPUS_interp_fixedZref(iZ+i-2, X, XC, XN, XO, logR, logT, raw_res(:,i), ierr)
    enddo

    nullify(work1)
    allocate(work1(nZ*pm_work_size))
    x_new(1)=safe_log10(my_Z)

    ! do the Z interpolation
    ! loop does over kap and its derivatives
    do i = 1, 3
       call interpolate_vector(nZ, kA% logZs(iZ-1:iZ+1), 1, &
            x_new, raw_res(i,:), v_new, &
            interp_pm, pm_work_size, &
            work1, junk, ierr)
       res(i) = v_new(1)
    enddo

    deallocate(work1)

  end subroutine AESOPUS_interp


  logical function outside_R_and_T_bounds(logR,logT)
    real(dp), intent(in) :: logR, logT
    outside_R_and_T_bounds = &
         logR < kA% min_logR .or. logR > kA% max_logR .or. &
         logT < kA% min_logT .or. logT > kA% max_logT
  end function outside_R_and_T_bounds


  subroutine AESOPUS_interp_fixedZref(iZ, X, XC, XN, XO, logR, logT, res, ierr)
    ! simple interpolation in each of X, fCO, fC, fN
    ! at present, interpolation is linear (order = 2)
    integer, intent(in) :: iZ
    real(dp), intent(in) :: X, XC, XN, XO, logR, logT
    real(dp), intent(out) :: res(3)
    integer, intent(out) :: ierr

    integer, parameter :: npts = 2
    real(dp), dimension(npts) :: wX, wfCO, wfC, wfN
    real(dp) :: raw_res(3)
    integer :: iX, ifCO, ifC, ifN
    integer :: i, j, k, l

    real(dp) :: Zref, fCO, fC, fN
    logical :: clipped_X, clipped_fCO, clipped_fC, clipped_fN

    ierr = 0
    res = 0d0
    iX=0; ifCO=0; ifC=0; ifN=0

    ! AESOPUS defines these quantities as follows

    ! fCO=log10(XC/XO)-log10(XC/XO)ref
    ! fC=log10(XC)-log10(XC)ref
    ! fN=log10(XN)-log10(XN)ref

    ! Note that the reference values are different for different solar abundance patterns

    Zref = kA% Zs(iZ)
    fCO = safe_log10(XC/XO) - kA% fCO_ref
    fC = safe_log10(XC/Zref) - kA% fC_ref
    fN = safe_log10(XN/Zref) - kA% fN_ref

    if (debug) then
       write(*,*) 'call to AESOPUS_interp_RT'
       write(*,*) 'logR = ', logR
       write(*,*) 'logT = ', logT
       write(*,*) 'Zref = ', Zref
       write(*,*) 'X = ', X
       write(*,*) 'fCO = ', fCO
       write(*,*) 'fC = ', fC
       write(*,*) 'fN = ', fN
    end if

    associate(ts => kA% ts(iZ))

      ! get weights for a clipped linear interpolation in each parameter
      call clipped_linear_weights(X, ts% num_Xs, ts% Xs, iX, wX, clipped_X)
      call clipped_linear_weights(fCO, ts% num_fCOs, ts% fCOs, ifCO, wfCO, clipped_fCO)
      call clipped_linear_weights(fC, ts% num_fCs, ts% fCs, ifC, wfC, clipped_fC)
      call clipped_linear_weights(fN, ts% num_fNs, ts% fNs, ifN, wfN, clipped_fN)

      ! cycles prevent wastefully calling interp_RT with zero weights
      do i = 1, npts
         if (wX(i) .eq. 0) cycle
         do j = 1, npts
            if (wfCO(j) .eq. 0) cycle
            do k = 1, npts
               if (wfC(k) .eq. 0) cycle
               do l = 1, npts
                  if (wfN(l) .eq. 0) cycle

                  if (debug) then
                     write(*,*) 'call to AESOPUS_interp_RT'
                     write(*,*) iX+i-1, ifCO+j-1, ifC+k-1,ifN+l-1
                     write(*,*) 'X = ', X, ts% Xs(iX+i-1), wX(i)
                     write(*,*) 'fCO = ', fCO, ts% fCOs(ifCO+j-1), wfCO(j)
                     write(*,*) 'fC = ', fC, ts% fCs(ifC+k-1), wfC(k)
                     write(*,*) 'fN = ', fN, ts% fNs(ifN+l-1), wfN(l)
                  end if

                  ! now do the call and collect the results

                  call AESOPUS_interp_RT(ts% t(iX+i-1,ifCO+j-1,ifC+k-1,ifN+l-1), logR, logT, raw_res, ierr)
                  if (ierr .ne. 0) return

                  res = res + wX(i)*wfCO(j)*wfC(k)*wfN(l) * raw_res

               end do
            end do
         end do
      end do

    end associate

  contains

    subroutine clipped_linear_weights(val, len, vec, loc, weights, clipped)

      ! calculate the weights for a linear interpolation
      ! clip to table edges

      use num_lib, only: binary_search

      real(dp), intent(in) :: val ! value
      integer, intent(in) :: len ! number of tabulated values
      real(dp), dimension(:), intent(in) :: vec ! tabulated values
      integer,  intent(out) :: loc ! vec(loc) <= val <= vec(loc+1)
      real(dp), dimension(2), intent(out) :: weights ! for linear interpolation
      logical, intent(out) :: clipped ! did we clip? if so, only locs(1)/weights(1) matter

      integer :: i

      weights = 0

      ! clip to range, if needed
      clipped = .false.
      if (val .le. vec(1)) then
         loc = 1
         weights(1) = 1d0
         weights(2) = 0d0
         clipped = .true.
      else if (val .ge. vec(len)) then
         loc = len
         weights(1) = 1d0
         weights(2) = 0d0
         clipped = .true.
      endif

      ! find location and calculate linear weights
      if (.not. clipped) then
         ! binary_search returns k between 1 and n-1 such that vec(k) <= val < vec(k+1)
         loc = binary_search(len, vec, len/2, val)
         weights(2) = (val - vec(loc)) / (vec(loc+1) - vec(loc))
         weights(1) = 1d0 - weights(2)
      endif

    end subroutine clipped_linear_weights

  end subroutine AESOPUS_interp_fixedZref


  subroutine AESOPUS_interp_RT(t, logR, logT, res, ierr)
    use interp_2d_lib_db, only: interp_evbicub_db
    type(AESOPUS_Table) :: t
    real(dp), intent(in) :: logR, logT
    real(dp), intent(out) :: res(3)
    integer, intent(out) :: ierr
    real(dp) :: raw_res(6)

    if (debug) then
       write(*,*) 'inside call to AESOPUS_interp_RT'
       write(*,*) 'X = ', t% X
       write(*,*) 'Z = ', t% Z
       write(*,*) 'fCO = ', t% fCO
       write(*,*) 'fC = ', t% fC
       write(*,*) 'fN = ', t% fN
    end if

    call interp_evbicub_db(logR, logT, &
         kA% logRs, kA% num_logRs, &
         ka% logTs, kA% num_logTs, &
         ilinx, iliny, &
         t% kap, ka% num_logRs, &
         ict, raw_res, ierr)

    res = raw_res(1:3)

  end subroutine AESOPUS_interp_RT


  subroutine AESOPUS_get(Zref, X, XC, XN, XO, logRho, logT, kap, &
       dlnkap_dlnRho, dlnkap_dlnT,ierr)
    real(dp), intent(in) :: Zref, X, XC, XN, XO
    real(dp), intent(in) :: logRho, logT
    real(dp), intent(out) :: kap, dlnkap_dlnRho, dlnkap_dlnT
    integer, intent(out) :: ierr
    real(dp) :: logR, res(3)

    ierr = 0
    logR = logRho - 3d0*logT + 18d0

    if (debug) write(*,*) Zref, X, XC, XN, XO, logR, logT
    call AESOPUS_interp(Zref, X, XC, XN, XO, logR, logT, res, ierr)
    if (ierr == 0) then
       kap = exp10(res(1))
       dlnkap_dlnRho = res(2)
       dlnkap_dlnT = res(3) - 3d0*res(2)
    else
       kap = -1d0
       dlnkap_dlnRho = 0d0
       dlnkap_dlnT = 0d0
    endif

  end subroutine AESOPUS_get


end module kap_aesopus
