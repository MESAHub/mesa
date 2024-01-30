program turb_plotter

  use utils_lib
  use const_lib
  use math_lib
  use turb
  
  implicit none

  integer :: ierr, op_err
  character (len=32) :: my_mesa_dir

  integer               :: nR0, nks, j, jmax, spectral_resolution
  real(dp), allocatable :: ks(:), res(:,:)
  real(dp)              :: tau, Pr, Pm, HB1, HB2, DB, R0
  real(dp)              :: HB(3)
  logical               :: FRG_withTC
  integer               :: iounit

  real(dp), parameter :: UNSET = -999
  
  namelist /plotter/ &
       tau, Pr, Pm, HB1, HB2, nR0, nks, spectral_resolution, FRG_withTC

  include 'formats'

  ierr = 0

  my_mesa_dir = '../..'
  call const_init(my_mesa_dir,ierr)
  if (ierr /= 0) then
     write(*,*) 'const_init failed'
     call mesa_error(__FILE__,__LINE__)
  end if

  call math_init()

  tau = UNSET
  Pr = UNSET
  Pm = UNSET
  HB1 = UNSET
  HB2 = UNSET
  nR0 = UNSET
  nks = UNSET
  spectral_resolution = UNSET
  
  ! get info from namelist
  open(newunit=iounit, file='inlist_plotter')
  read(iounit, nml=plotter)
  close(iounit)

  HB = [0d0, HB1, HB2]
  DB = Pr/Pm
  
  ! spectral_resolution must be odd integer, so promote to next odd number if even
  spectral_resolution = (spectral_resolution/2)*2 + 1

  allocate(ks(nks*2))
  allocate(res(nR0,7))
  res(:,:) = 0d0
  
  do j = 1,nks
     ! first nks entries are log space from 1e-6 to 0.1 (don't include endpoint)
     ks(j) = pow(10d0,-6d0 + (j-1)*(-1d0 + 6d0)/nks)

     ! last nks entries are linear space from 0.1 to 2
     ks(j+nks) = 0.1d0 + (j-1)*(2d0 - 0.1d0)/(nks - 1)
  end do
  
  ! loop stays interior to interval 1 < R0 < 1/tau,
  ! so need two extra points to define the endpoints.
  jmax = nR0 + 2
  
!$OMP PARALLEL DO PRIVATE(j, R0, op_err) SCHEDULE(dynamic,2)
  do j = 2,jmax-1
     R0 = 1d0 + (1d0/tau - 1d0)*(j - 1)/(jmax - 1)
     call set_res_for(j-1, R0, op_err)
     if(op_err /= 0) ierr = op_err
  end do
!$OMP END PARALLEL DO

  ! file for output
  open(newunit=iounit, file='turb_plotter.dat')
  ! header for use by numpy genfromtxt in plotter.py
  write(iounit,*) "# index R0 Dth_HG19_HB0 Dth_HG19_HB1 Dth_HG19_HB2 Dth_FRG24_HB0 Dth_FRG24_HB1 Dth_FRG24_HB2"

  ! write out results
  do j = 1,nR0
     write(iounit,*) j, res(j,:)
  end do
  
  deallocate(ks)
  deallocate(res)
  
contains

  subroutine set_res_for(j, R0, ierr)
     integer, intent(in)  :: j
     real(dp), intent(in) :: R0
     integer, intent(out) :: ierr

     real(dp) :: l2hat, lhat, lamhat, w
     integer  :: i
  
     ! calculate lamhat and l2hat
     call thermohaline_mode_properties(Pr, tau, R0, lamhat, l2hat, ierr)
     if (ierr /= 0) then
        write(*,*) 'thermohaline_mode_properties failed'
        call mesa_error(__FILE__,__LINE__)
     end if

     lhat = sqrt(l2hat)

     res(j,1) = R0

     ! Calculate results for 3 different magnetic field strengths (0,HB1,HB2)
     do i = 1,3
        ! use mode properties to calculate velocity w (Harrington model for now)
        call calc_hg19_w(HB(i), l2hat, lhat, lamhat, w, ierr)
        if (ierr /= 0) then
           write(*,*) 'calc_hg19_w failed'
           write(*,*) 'R0', R0
           write(*,*) '1/tau', 1/tau
           write(*,*) 'HB', HB(i)
           write(*,*) 'l2hat', l2hat
           write(*,*) 'lhat', lhat
           write(*,*) 'lamhat', lamhat
           write(*,*) 'w', w
           call mesa_error(__FILE__,__LINE__)
        end if
        
        ! caclulate resulting D/kappa_T as function of tau, Pr, R0
        ! KB = 1.24 for Harrington model
        res(j,i+1) = thermohaline_nusseltC(tau, w, lamhat, l2hat, 1.24d0) - 1d0

     end do

     ! Now calculate again for full FRG24 model (adds in Pm dependence)
     do i = 1,3
        call calc_frg24_w(Pr, tau, R0, HB(i), DB, ks, spectral_resolution, w, FRG_withTC, ierr, lamhat, l2hat)
        if (ierr /= 0) then
           write(*,*) 'calc_frg24_w failed'
           write(*,*) 'R0', R0
           write(*,*) '1/tau', 1/tau
           write(*,*) 'HB', HB(i)
           write(*,*) 'DB', DB
           write(*,*) 'l2hat', l2hat
           write(*,*) 'lhat', lhat
           write(*,*) 'lamhat', lamhat
           write(*,*) 'w', w
           write(*,*) 'ierr', ierr
           call mesa_error(__FILE__,__LINE__)
        end if
        write(*,*) "calc_frg24_w, R0, HB, w", R0, HB(i), w
        ! KB = 0.62 for Fraser model
        res(j,i+4) = thermohaline_nusseltC(tau, w, lamhat, l2hat, 0.62d0) - 1d0
     end do
     
   end subroutine set_res_for

end program turb_plotter
