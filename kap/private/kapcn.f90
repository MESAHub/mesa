! kapCN by Aaron Dotter 
! This is a MESA module that reads in and interpolates the low-T
! opacities by M.T. Lederer & B. Aringer, 2009, A&A, 494, 403
! see doc/ReadMe for details of the data files; it uses MESA
! modules for interpolation and a few other things
! kapCN_ZC=0.1644, kapCN_ZN=0.0532 

module kapcn

  use kap_def
  use num_lib, only: binary_search
  use const_def, only: dp
  use math_lib

  implicit none

  !local stuff
  logical, parameter :: debug = .false.
  character(len=32), parameter :: kapCN_param_file = 'kR_Z_fCN.data'

  !for kap_def

  !for 2-D interpolation
  integer :: ibcx=0, ibcy=0, ilinx=1, iliny=1
  integer :: ict(6) = [ 1, 1, 1, 0, 0, 0 ]
  real(dp), parameter :: bc(kapCN_num_logT) = 0d0

  
contains

  !for kapCN
  subroutine kapCN_init(show_info, ierr)
    use const_def, only: mesa_data_dir
    use kap_def, only: kapcn_is_initialized
    logical, intent(in) :: show_info
    integer, intent(out) :: ierr
    integer :: iZ
    type(kapCN_set), pointer :: k

    do iZ=1,num_kapCN_Zs
       k => kCN(iZ)
       k% not_loaded_yet = .true.
       k% Zbase = kapCN_Z(iZ)
       k% fC = kapCN_fC(:,iZ)
       k% fN = kapCN_fN(:,iZ)
    enddo

    call read_kapCN_tables(ierr)

    if(ierr==0) kapCN_is_initialized = .true.

  end subroutine kapCN_init

  subroutine kapCN_shutdown
    use kap_def, only: kapcn_is_initialized
    kapCN_is_initialized = .false.
  end subroutine kapCN_shutdown

  subroutine read_kapCN_tables(ierr)
    integer, intent(out) :: ierr
    character(len=256) :: filename
    character(len=4) :: prefix, string_Z(num_kapCN_Zs)
    character(len=6) :: tmp_Z, suffix
    integer :: i, io

    prefix = 'kR_Z'
    suffix = '.data'

    string_Z = ''

    filename = trim(kap_dir) // '/' // trim(kapCN_param_file)

    if(debug) write(*,*) '   parameter file = ', trim(filename)

    open(newunit=io,file=trim(filename))
    read(io,*) !skip first line
    do i=num_kapCN_Zs,1,-1 !read backwards into array so Z is increasing
       read(io,'(f7.5,11f9.1)') kapCN_Z(i), kapCN_fC(:,i), kapCN_fN(:,i)
    enddo
    close(io)

    do i=1,num_kapCN_Zs
       write(tmp_Z,'(1p,1e6.1e1)') kapCN_Z(i)
       string_Z(i) = tmp_Z(1:1) // tmp_Z(4:6)
    enddo

    kCN(:)% filename = prefix // string_Z // trim(suffix)

    do i=1,num_kapCN_Zs
       call read_one_kapCN_table(i,ierr)
    enddo

    kapCN_min_logR = minval(kapCN_logR)
    kapCN_max_logR = maxval(kapCN_logR)
    kapCN_min_logT = minval(kapCN_logT)
    kapCN_max_logT = maxval(kapCN_logT)

  end subroutine read_kapCN_tables


  subroutine read_one_kapCN_table(iZ,ierr)
    use interp_2d_lib_db, only: interp_mkbicub_db
    integer, intent(in) :: iZ
    integer, intent(out) :: ierr
    character(len=256) :: ascii_file, cache_file
    character(len=10) :: my_Z
    integer :: i, io, j
    real(dp) :: c_div_o, y
    real(dp) :: table(kapCN_num_logR,kapCN_num_logT)
    real(dp), pointer :: logR(:), logT(:)
    logical :: have_cache
    type(kapCN_set), pointer :: k

    ierr=0

    k => kCN(iZ)
    logR => kapCN_logR
    logT => kapCN_logT

    ! data files
    ascii_file = trim(kap_dir) // '/' // trim(k% filename)
    cache_file = trim(kap_dir) // '/cache/' // trim(k% filename(1:9)) // 'bin'

    if(debug) write(*,*) ' ascii file = ', trim(ascii_file)
    if(debug) write(*,*) ' cache file = ', trim(cache_file)

    !check for cached bin file; if it exists, read it and exit
    if(kapCN_use_cache)then
       inquire(file=trim(cache_file),exist=have_cache)
       if(have_cache)then
          !read cache
          open(newunit=io,file=trim(cache_file),form='unformatted',iostat=ierr)
          if(ierr/=0) write(*,*) &
               ' read_one_kapCN_table: problem opening cache file for read'
          read(io) logR, logT
          read(io) k% t(:)% X
          read(io) k% t(:)% Z
          read(io) k% t(:)% fC
          read(io) k% t(:)% fN
          read(io) k% t(:)% falpha
          do i=1,kapCN_num_tbl
             allocate(k% t(i)% kap(4*kapCN_tbl_size))
             read(io) k% t(i)% kap
          enddo
          close(io)
          !table is now fully loaded
          k% not_loaded_yet = .false.
          return
       endif
    endif

    !no cache file, so read ascii file
    open(newunit=io,file=trim(ascii_file),iostat=ierr)
    if(ierr/=0) write(*,*) &
         ' read_one_kapCN_table: problem opening ascii file for read'
    do i=1,4
       read(io,*) !skip header
    enddo
    read(io,'(6x,a10)') my_Z
    if(debug) write(*,*) ' has Z = ', my_Z, ' Zbase = ', k% Zbase

    do i=1,14
       read(io,*)
    enddo

    do i=1,kapCN_num_tbl
       read(io,*) j, k% t(i)% x, y, k% t(i)% z, c_div_o, &
            k% t(i)% fc, k% t(i)% fn, k% t(i)% falpha
    enddo

    read(io,*) !last line of #'s

    do i=1,kapCN_num_tbl
       do j=1,14
          read(io,*) !skip individual headers
       enddo
       read(io,'(5x,17f7.3)') logR
       do j=1,kapCN_num_logT
          read(io,'(f5.3,17f7.3)') logT(j), table(:,j)
       enddo
       allocate(k% t(i)% kap(4*kapCN_tbl_size))
       k% t(i)% kap(1:4*kapCN_tbl_size:4) = reshape(table,[kapCN_tbl_size])
    enddo

    close(io)

    !create interpolants for tables
    do i=1,kapCN_num_tbl
       call interp_mkbicub_db( logR, kapCN_num_logR, &
            logT, kapCN_num_logT, &
            k% t(i)% kap, kapCN_num_logR, &
            ibcx, bc, ibcx, bc, &
            ibcy, bc, ibcy, bc, &
            ilinx, iliny, ierr)
    end do

    !table is now fully loaded
    k% not_loaded_yet = .false.


    !write cache file
    open(newunit=io,file=trim(cache_file),form='unformatted',iostat=ierr)
    if(ierr/=0) write(*,*) &
         ' read_one_kapCN_table: problem opening cache file for write'
    write(io) logR, logT
    write(io) k% t(:)% x
    write(io) k% t(:)% z
    write(io) k% t(:)% fC
    write(io) k% t(:)% fN
    write(io) k% t(:)% falpha
    do i=1,kapCN_num_tbl
       write(io) k% t(i)% kap
    enddo
    close(io)

  end subroutine read_one_kapCN_table


  !this one is specifically for use with MESA


  subroutine kapCN_interp(Z,X,fC,fN,logR,logT,result,ierr)
    !result=(logKappa, dlogKappa/dlogR, dlogKappa/dlogT)
    use interp_1d_def, only: pm_work_size
    use interp_1d_lib, only: interpolate_vector, interp_pm
    real(dp), intent(in) :: Z, X, fC, fN, logR, logT
    real(dp), intent(out) :: result(3)
    integer, intent(out) :: ierr
    real(dp), pointer :: Z_ary(:), work1(:), logZ_ary(:)
    integer, parameter :: nZ = 4
    integer :: i,iZ
    real(dp) :: my_Z, res(3,nZ), x_new(1), v_new(1)
    character(len=32) :: junk

    result=0d0; iZ=0; ierr=0

    if(.not.kapCN_is_initialized)then
       write(*,*) ' kapCN is not initialized; call kapCN_init()'
       ierr=-99
       return
    endif

    if(outside_R_and_T_bounds(logR,logT))then
       !write(*,*) 'kapCN_interp: logR, logT outside of table bounds'
       ierr=-1
       return
    endif

    Z_ary => kapCN_Z
    my_Z = max(min(Z,Z_ary(num_kapCN_Zs)),Z_ary(1))
    iZ=binary_search(num_kapCN_Zs,Z_ary, iZ, Z)
    iZ = max(nz/2,min(iZ,num_kapCN_Zs-nz/2))

    !check to see if exact match in Z, then just need one call
    if(Z==Z_ary(iZ))then
       call kapCN_interp_fixedZ(iZ,X,fC,fN,logR,logT,result)
       return
    endif

    !else do the full Z interpolation
    do i=1,nZ
       call kapCN_interp_fixedZ(iZ+i-2,X,fC,fN,logR,logT,res(:,i))
    enddo

    nullify(work1)
    allocate(work1(nZ*pm_work_size))
    x_new(1)=log10(my_Z)

    allocate(logZ_ary(iZ-1:iZ+2))
    
    do i=iZ-1,iZ+2
        logZ_ary(i) = log10(z_ary(i))
    end do

    do i=1,3
       call interpolate_vector(nZ,logZ_ary, 1, &
            x_new, res(i,:), v_new, &
            interp_pm, pm_work_size, &
            work1, junk, ierr )
       result(i)=v_new(1)
    enddo

    deallocate(work1, logz_ary)

  end subroutine kapCN_interp


  logical function outside_R_and_T_bounds(logR,logT)
    real(dp), intent(in) :: logR, logT     
    outside_R_and_T_bounds = &
         logR < kapCN_min_logR .or. logR > kapCN_max_logR .or. &
         logT < kapCN_min_logT .or. logT > kapCN_max_logT
  end function outside_R_and_T_bounds


  subroutine kapCN_interp_fixedZ(iZ,X,fC,fN,logR,logT,result)
    !using simple quadratic interpolation in each of X, fC, fN
    integer, intent(in) :: iZ
    real(dp), intent(in) :: X, fC, fN, logR, logT
    real(dp), intent(out) :: result(3)
    integer, parameter :: my_num_kapCN_fCs = 3, f1=num_kapCN_fCs*num_kapCN_Xs, f2=num_kapCN_fCs
    real(dp) :: my_X, my_fC, my_fN, wX(num_kapCN_Xs), wfN(num_kapCN_fNs), wfC(my_num_kapCN_fCs)
    real(dp) :: res(3)
    integer :: iX, ifC, ifN, i,j,tbl,n,nlo,nhi,ierr
    real(dp), pointer :: X_ary(:), fC_ary(:), fN_ary(:)

    result = 0d0; iX=0; ifC=0; ifN=0
    X_ary => kapCN_X; fC_ary => kapCN_fC(:,iZ); fN_ary => kapCN_fN(:,iZ)

    !limit input quantities to lie within tabulated range
    my_X = max(min(X, X_ary(num_kapCN_Xs)),X_ary(1))
    my_fC = max(min(fC,fC_ary(num_kapCN_fCs)),fC_ary(1))
    my_fN = max(min(fN,fN_ary(num_kapCN_fNs)),fN_ary(1))

    !locate inputs in arrays
    iX = binary_search(num_kapCN_Xs, X_ary, iX, my_X)
    ifC= binary_search(num_kapCN_fCs, fC_ary, ifC, my_fC)
    ifN= binary_search(num_kapCN_fNs, fN_ary, ifN, my_fN)

    !make sure 1 < ifC < num_kapCN_fCs
    ifC = max(2,min(ifC,num_kapCN_fCs-1))

    !interpolation coefficients in X
    call interp(X_ary,wX,my_X,num_kapCN_Xs)

    !interpolation coefficients in fN
    call interp(fN_ary(:),wfN,my_fN,num_kapCN_fNs)

    !interpolation coefficients in fC
    call interp(fC_ary(ifC-1:ifC+1),wfC,my_fC,my_num_kapCN_fCs)

    do i=1,num_kapCN_fNs
       do j=1,num_kapCN_Xs
          do n=1,my_num_kapCN_fCs
             res=0d0
             tbl = f1*(i-1) + f2*(j-1) + (ifC+n-2)
             nlo = 3*(n-1)+1
             nhi = nlo+2
             call kapCN_interp_RT(iZ,tbl,logR,logT,res,ierr)
             result = result + wfN(i)*wX(j)*wfC(n)*res
          enddo
       enddo
    enddo

    if(debug) write(*,'(1p9e12.4)') my_X, my_fC, my_fN, result

  contains

    subroutine interp(a,b,x,n)
      ! {a} are the tabulated values for use in interpolation
      ! {b} are coefficients of the interpolating polynomial
      !  x  is the abscissa to be interpolated
      !  n  is the number of points to be used, interpolating polynomial
      !     has order n-1 
      integer, intent(in) :: n
      real(dp), intent(in) :: a(n), x
      real(dp), intent(out) :: b(n)
      integer :: i,j
      do i=1,n
         b(i)=1d0
         do j=1,n
            if(j/=i) b(i)=b(i)*(x-a(j))/(a(i)-a(j))
         enddo
      enddo
    end subroutine interp

  end subroutine kapCN_interp_fixedZ


  subroutine kapCN_interp_RT(iZ,tbl,logR,logT,result,ierr)
    use interp_2d_lib_db, only: interp_evbicub_db
    integer, intent(in) :: iZ, tbl
    real(dp), intent(in) :: logR, logT
    real(dp), intent(out) :: result(3)
    integer, intent(out) :: ierr
    real(dp) :: res(6)
    real(dp), pointer :: logR_ary(:), logT_ary(:)
    type(kapCN_set), pointer :: k

    k => kCN(iZ)
    logR_ary => kapCN_logR
    logT_ary => kapCN_logT

    call interp_evbicub_db( logR, logT, &
         logR_ary, kapCN_num_logR, &
         logT_ary, kapCN_num_logT, &
         ilinx, iliny, &
         k% t(tbl)% kap, kapCN_num_logR, &
         ict, res, ierr)

    result = res(1:3)

  end subroutine kapCN_interp_RT

  subroutine kapCN_get(Z,X,fC,fN,logRho,logT,kappa, &
       dlnkap_dlnRho,dlnkap_dlnT,ierr)
    real(dp), intent(in) :: Z, X, fC, fN, logRho, logT
    real(dp), intent(out) :: kappa
    real(dp), intent(out) :: dlnkap_dlnRho, dlnkap_dlnT
    integer, intent(out) :: ierr
    real(dp) :: logR, result(3)

    ierr = 0
    logR = logRho - 3d0*logT + 18d0

    call kapCN_interp(Z,X,fC,fN,logR,logT,result,ierr)
    if(ierr==0)then
       kappa = exp10(result(1))
       dlnkap_dlnRho = result(2)
       dlnkap_dlnT = result(3) - 3*result(2)
    else
       kappa = -1d0
       dlnkap_dlnRho = 0d0
       dlnkap_dlnT = 0d0
    endif

  end subroutine kapCN_get


end module kapcn
