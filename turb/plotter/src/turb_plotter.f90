program turb_plotter

  use utils_lib
  use const_lib
  use math_lib
  
  implicit none

  integer :: ierr
  character (len=32) :: my_mesa_dir

  integer :: nR0, j
  real(dp) :: tau, Pr, Pm, R0, res
  integer :: iounit

  real(dp), parameter :: UNSET = -999
  
  namelist /plotter/ &
       tau, Pr, Pm, nR0

  include 'formats'

  ierr = 0
  res = 0d0
  
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
  nR0 = UNSET
  
  ! get info from namelist
  open(newunit=iounit, file='inlist_plotter')
  read(iounit, nml=plotter)
  close(iounit)
  
  ! file for output
  open(newunit=iounit, file='turb_plotter.dat')

  ! loop for making calls to calculate things
  do j = 1,nR0
     R0 = 1 + (1/tau)*(j-1)/(nR0 -1)

     ! caclulate res as function of tau, Pr, Pm, R0 here
     res = 37d0
     
     write(iounit,*) j, R0, res
  end do
  
end program turb_plotter
