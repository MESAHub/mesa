program turb_plotter

   use utils_lib
   use const_lib
   use math_lib
   use turb_lib, only: set_info_HG19, set_info_FRG24
   use turb_def

   implicit none

   integer :: ierr, op_err
   character (len=32) :: my_mesa_dir

   integer                :: n_R_0, nks, j, jmax, N, safety
   real(dp), allocatable  :: k_z(:), res(:,:)
   real(dp)               :: tau, Pr, Pm, H_B_1, H_B_2, D_B, R_0
   real(dp)               :: H_B(3)
   integer                :: iounit

   real(dp), parameter :: UNSET = -999

   namelist /plotter/ &
      tau, Pr, Pm, H_B_1, H_B_2, n_R_0, nks, N, safety

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
   H_B_1 = UNSET
   H_B_2 = UNSET
   n_R_0 = UNSET
   nks = UNSET
   N = UNSET

   ! get info from namelist
   open(newunit=iounit, file='inlist_plotter')
   read(iounit, nml=plotter)
   close(iounit)

   H_B = [0d0, H_B_1, H_B_2]
   D_B = Pr/Pm

   allocate(res(n_R_0,7))

   ! loop stays interior to interval 1 < R_0 < 1/tau,
   ! so need two extra points to define the endpoints.
   jmax = n_R_0 + 2

   !$OMP PARALLEL DO PRIVATE(j, R_0, op_err) SCHEDULE(dynamic,2)
   do j = 2,jmax-1

      R_0 = 1d0 + (1d0/tau - 1d0)*(j - 1)/(jmax - 1)
      call set_res_for(j-1, R_0, op_err)
      if(op_err /= 0) ierr = op_err

   end do
   !$OMP END PARALLEL DO

   ! file for output
   open(newunit=iounit, file='turb_plotter.dat')
   ! header for use by numpy genfromtxt in plotter.py
   write(iounit,*) "# index R_0 Dth_HG19_HB0 Dth_HG19_HB1 Dth_HG19_HB2 Dth_FRG24_HB0 Dth_FRG24_HB1 Dth_FRG24_HB2"

   ! write out results
   do j = 1, n_R_0
      write(iounit,*) j, res(j,:)
   end do

contains

   subroutine set_res_for(j, R_0, ierr)

      integer, intent(in)  :: j
      real(dp), intent(in) :: R_0
      integer, intent(out) :: ierr

      type(th_info_t) :: th_info
      integer  :: i

      ! Set parameters in th_info

      th_info%Pr = Pr
      th_info%tau = tau
      th_info%R_0 = R_0
      th_info%r = (th_info%R_0 - 1._dp)/(1._dp/th_info%tau - 1._dp)
      th_info%D_B = D_B
      th_info%K_C = 1._dp ! Required so that D_thrm = D_thrm/K_C

      res(j, 1) = R_0

      ! Calculate HG19 results for 3 different magnetic field strengths (0,HB1,HB2)

      do i = 1,3

         th_info%H_B = H_B(i)

         call set_info_HG19(th_info, ierr)
         if (ierr /= 0) then
            call mesa_error(__FILE__,__LINE__)
         end if

         res(j,i+1) = th_info%D_thrm

      end do

      ! Now calculate again for full FRG24 model (adds in Pm dependence)

      do i = 1,3

         th_info%H_B = H_B(i)

         call set_info_FRG24(safety, nks, N, th_info, ierr)
         if (ierr /= 0) then
            call mesa_error(__FILE__,__LINE__)
         end if

         res(j,i+4) = th_info%D_thrm

      end do

   end subroutine set_res_for

end program turb_plotter
