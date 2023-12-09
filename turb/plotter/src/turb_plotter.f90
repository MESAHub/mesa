program turb_plotter

  use utils_lib
  use const_lib
  use math_lib
  use turb
  
  implicit none

  integer :: ierr
  character (len=32) :: my_mesa_dir

  integer  :: nR0, i, j, jmax, spectral_resolution
  real(dp) :: ks(100)
  real(dp) :: tau, Pr, Pm, HB1, HB2, DB, R0
  real(dp) :: l2hat, lhat, lamhat, w
  real(dp) :: HB(3), res(5)
  integer  :: iounit

  real(dp), parameter :: UNSET = -999
  
  namelist /plotter/ &
       tau, Pr, Pm, HB1, HB2, nR0

  include 'formats'

  ierr = 0

  spectral_resolution = 17 ! must be odd. Should promote to inlist option eventually

  do i = 1,50
     ! first 50 entries are log space from 1e-6 to 0.1 (don't include endpoint)
     ks(i) = pow(10d0,-6d0 + (i-1)*(-1d0 + 6d0)/50d0)

     ! last 50 entries are linear space from 0.1 to 2
     ks(i+50) = 0.1d0 + (i-1)*(2d0 - 0.1d0)/49d0
  end do
  
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
  
  ! get info from namelist
  open(newunit=iounit, file='inlist_plotter')
  read(iounit, nml=plotter)
  close(iounit)

  HB = [0d0, HB1, HB2]
  res(:) = 0d0
  
  ! file for output
  open(newunit=iounit, file='turb_plotter.dat')

  ! loop stays interior to interval 1 < R0 < 1/tau,
  ! so need two extra points to define the endpoints.
  jmax = nR0 + 2
  
  ! loop for making calls to calculate things
  do j = 2,jmax-1
     R0 = 1d0 + (1d0/tau - 1d0)*(j - 1)/(jmax - 1)

     ! calculate lamhat and l2hat
     call thermohaline_mode_properties(Pr, tau, R0, lamhat, l2hat, ierr)
     if (ierr /= 0) then
        write(*,*) 'thermohaline_mode_properties failed'
        call mesa_error(__FILE__,__LINE__)
     end if

     lhat = sqrt(l2hat)

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
        res(i) = thermohaline_nusseltC(tau, w, lamhat, l2hat) - 1d0

     end do

     ! Now calculate again for full FRG24 model (adds in Pm dependence)
     DB = Pr/Pm
     do i = 2,3
        call calc_frg24_w(Pr, tau, R0, HB(i), DB, ks, spectral_resolution, w, ierr, lamhat, l2hat)
        if (ierr /= 0) then
           write(*,*) 'calc_frg24_w failed'
           write(*,*) 'R0', R0
           write(*,*) '1/tau', 1/tau
           write(*,*) 'HB', HB(i)
           write(*,*) 'l2hat', l2hat
           write(*,*) 'lhat', lhat
           write(*,*) 'lamhat', lamhat
           write(*,*) 'w', w
           call mesa_error(__FILE__,__LINE__)
        end if
        print *, "calc_frg24_w, R0, HB, w", R0, HB(i), w
        res(i+2) = thermohaline_nusseltC(tau, w, lamhat, l2hat) - 1d0
     end do
     
     write(iounit,*) j-1, R0, res
  end do
  
end program turb_plotter
