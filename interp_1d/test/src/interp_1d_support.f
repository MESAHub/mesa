      module interp_1d_support
      use const_lib
      use math_lib
      use interp_1d_def
      use interp_1d_lib
      use utils_lib, only: mesa_error

      implicit none


      integer, parameter :: num_xpts = 31
      integer, parameter :: num_ypts = 35
      integer, parameter :: sz_per_pt = 4

      contains
      
      
      subroutine do_test

         write(*, *)
         call test_min(.true.)      

         write(*, *)
         call test_min(.false.)  

         write(*, *)
         call test1(.true.)

         write(*, *)
         call test1(.false.)
      
         write(*, *)
         call test2(.true.)
      
         write(*, *)
         call test2(.false.)
      
         write(*, *)
         call test3(.true.)
      
         write(*, *)
         call test3(.false.)
      
         write(*, *)
         call test_4pt

         write(*, *)
         call test_interp_3_to_2
         
      end subroutine do_test
      
      
      subroutine test_data(n,U)
         integer, intent(in) :: n
         real(dp), intent(out) :: U(n)
         
        U(  1)=          0.3691382918814334D0
        U(  2)=          0.3695899829574669D0
        U(  3)=          0.3699633032242006D0
        U(  4)=          0.3709860746618878D0
        U(  5)=          0.3724175838188959D0
        U(  6)=          0.3742482756438662D0
        U(  7)=          0.3764768199547272D0
        U(  8)=          0.3790910321656262D0
        U(  9)=          0.3820871425096383D0
        U( 10)=          0.3854500131142883D0
        U( 11)=          0.3891727128876379D0
        U( 12)=          0.3932377816410259D0
        U( 13)=          0.3976355585696444D0
        U( 14)=          0.4023461172079230D0
        U( 15)=          0.4073570798427897D0
        U( 16)=          0.4126461990878889D0
        U( 17)=          0.4181985492965749D0
        U( 18)=          0.4239897218739848D0
        U( 19)=          0.4300024071530404D0
        U( 20)=          0.4362102270543894D0
        U( 21)=          0.4425936978527515D0
        U( 22)=          0.4491247203356548D0
        U( 23)=          0.4557819260544076D0
        U( 24)=          0.4625358769752868D0
        U( 25)=          0.4693638888826394D0
        U( 26)=          0.4762362700345437D0
        U( 27)=          0.4831315999154387D0
        U( 28)=          0.4900124769496845D0
        U( 29)=          0.4968509249135122D0
        U( 30)=          0.5036231656397859D0
        U( 31)=          0.5102978108992456D0
        U( 32)=          0.5168503329341820D0
        U( 33)=          0.5232493058916751D0
        U( 34)=          0.5294708374154181D0
        U( 35)=          0.5354843215442489D0
        U( 36)=          0.5412671682413143D0
        U( 37)=          0.5467901937111218D0
        U( 38)=          0.5520326792853251D0
        U( 39)=          0.5569674077513700D0
        U( 40)=          0.5615760727117703D0
        U( 41)=          0.5658339305000317D0
        U( 42)=          0.5697255981615711D0
        U( 43)=          0.5732292848999179D0
        U( 44)=          0.5763330504836860D0
        U( 45)=          0.5790185013283228D0
        U( 46)=          0.5812774415045446D0
        U( 47)=          0.5830949206727917D0
        U( 48)=          0.5844671751637984D0
        U( 49)=          0.5853828024764719D0
        U( 50)=          0.5858432670435806D0
        U( 51)=          0.5858432670435451D0
        U( 52)=          0.5853828024764413D0
        U( 53)=          0.5844671751637710D0
        U( 54)=          0.5830949206727681D0
        U( 55)=          0.5812774415045240D0
        U( 56)=          0.5790185013283052D0
        U( 57)=          0.5763330504836706D0
        U( 58)=          0.5732292848999045D0
        U( 59)=          0.5697255981615594D0
        U( 60)=          0.5658339305000216D0
        U( 61)=          0.5615760727117616D0
        U( 62)=          0.5569674077513626D0
        U( 63)=          0.5520326792853185D0
        U( 64)=          0.5467901937111161D0
        U( 65)=          0.5412671682413093D0
        U( 66)=          0.5354843215442447D0
        U( 67)=          0.5294708374154145D0
        U( 68)=          0.5232493058916718D0
        U( 69)=          0.5168503329341793D0
        U( 70)=          0.5102978108992431D0
        U( 71)=          0.5036231656397836D0
        U( 72)=          0.4968509249135103D0
        U( 73)=          0.4900124769496828D0
        U( 74)=          0.4831315999154373D0
        U( 75)=          0.4762362700345427D0
        U( 76)=          0.4693638888826384D0
        U( 77)=          0.4625358769752859D0
        U( 78)=          0.4557819260544067D0
        U( 79)=          0.4491247203356541D0
        U( 80)=          0.4425936978527509D0
        U( 81)=          0.4362102270543888D0
        U( 82)=          0.4300024071530398D0
        U( 83)=          0.4239897218739841D0
        U( 84)=          0.4181985492965744D0
        U( 85)=          0.4126461990878886D0
        U( 86)=          0.4073570798427895D0
        U( 87)=          0.4023461172079228D0
        U( 88)=          0.3976355585696443D0
        U( 89)=          0.3932377816410258D0
        U( 90)=          0.3891727128876379D0
        U( 91)=          0.3854500131142883D0
        U( 92)=          0.3820871425096385D0
        U( 93)=          0.3790910321656264D0
        U( 94)=          0.3764768199547278D0
        U( 95)=          0.3742482756438669D0
        U( 96)=          0.3724175838188968D0
        U( 97)=          0.3709860746618888D0
        U( 98)=          0.3699633032242018D0
        U( 99)=          0.3695899829574566D0
        U(100)=          0.3691382918814347D0
      
      end subroutine test_data
      
      
      
      subroutine test1(increasing)
         logical, intent(in) :: increasing
         integer, parameter :: n = 100
         real(dp) :: U(n)
         
         integer, parameter :: nvals = 4
         integer, parameter :: nwork = max(pm_work_size, mp_work_size)
         real(dp), target :: work_ary(n*nwork)
         real(dp), pointer :: work(:)
         real(dp) :: init_xs(n), xs(nvals), vals(nvals), dx
         integer :: i, ierr
         real(dp), target :: f_ary(4*n)
         real(dp), pointer :: f1(:), f(:,:)
         f1 => f_ary
         f(1:4,1:n) => f1(1:4*n)
         
         include 'formats'
         
         write(*, *) 'test1', increasing
         
         work => work_ary
         dx = 0.594059405940594d0
         
         do i = 1, n
            init_xs(i) = dx*(i-1)
         end do
         
         if (.not. increasing) init_xs(:) = -init_xs(:)
      
         call test_data(n, U)
         
         f(1, 1:n) = U(1:n)
         call interp_m3q(init_xs, n, f1, nwork, work, 'test1', ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         write(*, 1) 'f(2, 1)', f(2, 1)
         write(*, 1) 'f(2, 2)', f(2, 2)
         write(*, 1) 'f(2, 3)', f(2, 3)
         write(*, *)
         
         xs = (/ 10d0, 10.4d0, 10.8d0, 11d0 /)
         if (.not. increasing) xs(:) = -xs(:)
         call interp_values(init_xs, n, f1, nvals, xs, vals, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         write(*, 1) 'z(10.0)', vals(1)
         write(*, 1) 'z(10.4)', vals(2)
         write(*, 1) 'z(10.8)', vals(3)
         write(*, 1) 'z(11.0)', vals(4)
         write(*, *)
      
      end subroutine test1
      
      
      subroutine test2(increasing)
         logical, intent(in) :: increasing
         integer, parameter :: n = 43
         real(dp) :: U(n)
         
         logical, parameter :: show_errors = .false.
         integer, parameter :: nvals = 7
         real(dp) :: init_xs(n), xs(nvals), vals(nvals), dx
         integer, parameter :: nwork = max(pm_work_size, mp_work_size)
         real(dp), target :: work_ary(n*nwork)
         real(dp), pointer :: work(:)
         integer :: i, ierr
         real(dp), target :: f_ary(4*n)
         real(dp), pointer :: f1(:), f(:,:)
         f1 => f_ary
         f(1:4,1:n) => f1(1:4*n)
         
         include 'formats'

         write(*, *) 'test2', increasing
         
         work => work_ary
         dx = 1d0 / (n-1)
         do i = 1, n
            init_xs(i) = 4*dx*(i-1)
            if (.not. increasing) init_xs(i) = -init_xs(i)
            U(i) = sin(init_xs(i))
         end do
         
         f(1, 1:n) = U(1:n)
         call interp_m3q(init_xs, n, f1, nwork, work, 'test2', ierr)
         if (ierr /= 0) then
            call mesa_error(__FILE__,__LINE__)
         end if
         write(*, *)
         
         xs = (/ 0d0, pi/4, pi/3, pi/2, 2*pi/3, 3*pi/4, 19*pi/20 /)
         if (.not. increasing) xs(:) = -xs(:)
         call interp_values(init_xs, n, f1, nvals, xs, vals, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         
         if (show_errors) then
            do i = 1, nvals
               write(*, 2) 'x val exact err', i, xs(i), vals(i), sin(xs(i)), vals(i) - sin(xs(i))
            end do
         else
            do i = 1, nvals
               write(*, 2) 'x val exact', i, xs(i), vals(i), sin(xs(i))
            end do
         end if
         write(*, *)
      
         call integrate_values(init_xs, n, f1, nvals, xs, vals, ierr)
         if (ierr /= 0) then
            call mesa_error(__FILE__,__LINE__)
         end if

         if (show_errors) then
            do i = 2, nvals
               write(*, 2) 'x val exact err', i, xs(i), vals(i), cos(xs(i-1)) - cos(xs(i)), 
     >                     vals(i) - (cos(xs(i-1)) - cos(xs(i)))
            end do
         else
            do i = 2, nvals
               write(*, 2) 'x val exact', i, xs(i), vals(i), cos(xs(i-1)) - cos(xs(i))
            end do
         end if
         write(*, *)
      
      end subroutine test2
      
      
      subroutine test3(increasing)
         logical, intent(in) :: increasing
         integer, parameter :: n = 100
         real(dp) :: U(n)
         
         integer, parameter :: nvals = 3
         integer, parameter :: nwork = max(pm_work_size, mp_work_size)
         real(dp), target :: work_ary(n*nwork)
         real(dp), pointer :: work(:)
         real(dp) :: init_xs(n), xs(nvals), vals(nvals), dx
         integer :: i, ierr

         real(dp), target :: f_ary(4*n)
         real(dp), pointer :: f1(:), f(:,:)
         f1 => f_ary
         f(1:4,1:n) => f1(1:4*n)
         
         include 'formats'
         
         write(*, *) 'test3', increasing

         work => work_ary
         dx = 0.594059405940594d0
         
         do i = 1, n
            init_xs(i) = dx*(i-1)
         end do
         
         if (.not. increasing) init_xs(:) = -init_xs(:)
      
         call test_data(n, U)
         
         f(1, 1:n) = U(1:n)
         call interp_pm(init_xs, n, f1, nwork, work, 'test3', ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         write(*, 1) 'f(2, 1)', f(2, 1)
         write(*, 1) 'f(2, 2)', f(2, 2)
         write(*, 1) 'f(2, 3)', f(2, 3)
         write(*, *)
         
         xs = (/ 10d0, 10.4d0, 10.8d0 /)
         if (.not. increasing) xs(:) = -xs(:)
         call interp_values(init_xs, n, f1, nvals, xs, vals, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         write(*, 1) 'z(10.0)', vals(1)
         write(*, 1) 'z(10.4)', vals(2)
         write(*, 1) 'z(10.8)', vals(3)
         write(*, *)
      
      end subroutine test3
      
      
      subroutine test_min(increasing)
         logical, intent(in) :: increasing
         integer, parameter :: n = 100
         real(dp) :: U(n)
         
         integer, parameter :: nvals = 2
         integer, parameter :: nwork = max(pm_work_size, mp_work_size)
         real(dp), target :: work_ary(n*nwork)
         real(dp), pointer :: work(:)
         real(dp) :: init_xs(n), xs(nvals), vals(nvals), dx
         integer :: i, ierr

         real(dp), target :: f_ary(4*n)
         real(dp), pointer :: f1(:), f(:,:)
         f1 => f_ary
         f(1:4,1:n) => f1(1:4*n)
         
         include 'formats'
         
         
         write(*, *) 'test3', increasing

         dx = 0.594059405940594d0
         work => work_ary
         
         do i = 1, n
            init_xs(i) = dx*(i-1)
         end do
         
         if (.not. increasing) init_xs(:) = -init_xs(:)
      
         call test_data(n, U)
         
         f(1, 1:n) = U(1:n)
         call interp_pm(init_xs, n, f1, nwork, work, 'test_min', ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         write(*, 1) 'f(2, 1)', f(2, 1)
         write(*, 1) 'f(2, 2)', f(2, 2)
         write(*, *)
         
         xs = (/ 10d0, 10.8d0 /)
         if (.not. increasing) xs(:) = -xs(:)
         call interp_values(init_xs, n, f1, nvals, xs, vals, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         write(*, 1) 'z(10.0)', vals(1)
         write(*, 1) 'z(10.8)', vals(2)
         write(*, *)
      
      end subroutine test_min
      
      
      subroutine test_interp_3_to_2
         real(dp) :: pdqm1, pdq00, ndqm1, ndq00, pfm1, pf00, pfp1, nf00, nfp1
         integer :: ierr
         
         include 'formats'
         write(*, *) 'test_interp_3_to_2'
         pdqm1  =  2.0114947208182310D-03
         pdq00  =  2.0097307373083619D-03
         ndqm1  =  1.5784933404892154D-03
         ndq00  =  1.3270582904786616D-03
         pfm1  =  1.8984458478374221D+30
         pf00  =  1.4579738233866260D+30
         pfp1  =  1.1136297079131214D+30
         call interp_3_to_2(
     >      pdqm1, pdq00, ndqm1, ndq00, pfm1, pf00, pfp1, nf00, nfp1, 'test_interp_3_to_2', ierr)
         write(*,*) 'ierr', ierr
         write(*,1) 'nf00', nf00
         write(*,1) 'nfp1', nfp1
      end subroutine test_interp_3_to_2
      
      
      subroutine test_4pt
         integer, parameter :: n = 6
         real(dp) :: x(n), y(n), a(3), dx, result, exact
         integer :: ierr, i
         integer, parameter :: nwork = max(pm_work_size, mp_work_size)
         real(dp), target :: work_ary(n*nwork)
         real(dp), pointer :: work(:)
         real(dp), target :: f_ary(4*n)
         real(dp), pointer :: f1(:)
         
         include 'formats'
         f1 => f_ary
         work => work_ary
         x = (/ pi/4, pi/3, pi/2, 2*pi/3, 3*pi/4, 19*pi/20 /)
         do i=1, n
            y(i) = sin(x(i))
         end do
         call interp_4pt_pm(x(2:5), y(2:5), a)
         dx = (x(4)-x(3))/2
         result = y(3) + dx*(a(1) + dx*(a(2) + dx*a(3)))
         exact = sin(x(3)+dx)
         write(*, *) 'test_4pt'
         write(*, 1) 'res exact err', result, exact, abs((result-exact)/exact)
      
      end subroutine test_4pt


      end module interp_1d_support




