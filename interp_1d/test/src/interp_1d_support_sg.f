      module interp_1d_support_sg
      use const_lib
      use math_lib
      use interp_1d_def
      use interp_1d_lib_sg
      use utils_lib, only: mesa_error

      implicit none


      integer, parameter :: num_xpts = 31
      integer, parameter :: num_ypts = 35
      integer, parameter :: sz_per_pt = 4

      contains
      
      
      subroutine do_test_sg

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
         
      end subroutine do_test_sg
      
      
      subroutine test_data(n,U)
         integer, intent(in) :: n
         real, intent(out) :: U(n)
         
        U(  1)=          0.3691382918814334
        U(  2)=          0.3695899829574669
        U(  3)=          0.3699633032242006
        U(  4)=          0.3709860746618878
        U(  5)=          0.3724175838188959
        U(  6)=          0.3742482756438662
        U(  7)=          0.3764768199547272
        U(  8)=          0.3790910321656262
        U(  9)=          0.3820871425096383
        U( 10)=          0.3854500131142883
        U( 11)=          0.3891727128876379
        U( 12)=          0.3932377816410259
        U( 13)=          0.3976355585696444
        U( 14)=          0.4023461172079230
        U( 15)=          0.4073570798427897
        U( 16)=          0.4126461990878889
        U( 17)=          0.4181985492965749
        U( 18)=          0.4239897218739848
        U( 19)=          0.4300024071530404
        U( 20)=          0.4362102270543894
        U( 21)=          0.4425936978527515
        U( 22)=          0.4491247203356548
        U( 23)=          0.4557819260544076
        U( 24)=          0.4625358769752868
        U( 25)=          0.4693638888826394
        U( 26)=          0.4762362700345437
        U( 27)=          0.4831315999154387
        U( 28)=          0.4900124769496845
        U( 29)=          0.4968509249135122
        U( 30)=          0.5036231656397859
        U( 31)=          0.5102978108992456
        U( 32)=          0.5168503329341820
        U( 33)=          0.5232493058916751
        U( 34)=          0.5294708374154181
        U( 35)=          0.5354843215442489
        U( 36)=          0.5412671682413143
        U( 37)=          0.5467901937111218
        U( 38)=          0.5520326792853251
        U( 39)=          0.5569674077513700
        U( 40)=          0.5615760727117703
        U( 41)=          0.5658339305000317
        U( 42)=          0.5697255981615711
        U( 43)=          0.5732292848999179
        U( 44)=          0.5763330504836860
        U( 45)=          0.5790185013283228
        U( 46)=          0.5812774415045446
        U( 47)=          0.5830949206727917
        U( 48)=          0.5844671751637984
        U( 49)=          0.5853828024764719
        U( 50)=          0.5858432670435806
        U( 51)=          0.5858432670435451
        U( 52)=          0.5853828024764413
        U( 53)=          0.5844671751637710
        U( 54)=          0.5830949206727681
        U( 55)=          0.5812774415045240
        U( 56)=          0.5790185013283052
        U( 57)=          0.5763330504836706
        U( 58)=          0.5732292848999045
        U( 59)=          0.5697255981615594
        U( 60)=          0.5658339305000216
        U( 61)=          0.5615760727117616
        U( 62)=          0.5569674077513626
        U( 63)=          0.5520326792853185
        U( 64)=          0.5467901937111161
        U( 65)=          0.5412671682413093
        U( 66)=          0.5354843215442447
        U( 67)=          0.5294708374154145
        U( 68)=          0.5232493058916718
        U( 69)=          0.5168503329341793
        U( 70)=          0.5102978108992431
        U( 71)=          0.5036231656397836
        U( 72)=          0.4968509249135103
        U( 73)=          0.4900124769496828
        U( 74)=          0.4831315999154373
        U( 75)=          0.4762362700345427
        U( 76)=          0.4693638888826384
        U( 77)=          0.4625358769752859
        U( 78)=          0.4557819260544067
        U( 79)=          0.4491247203356541
        U( 80)=          0.4425936978527509
        U( 81)=          0.4362102270543888
        U( 82)=          0.4300024071530398
        U( 83)=          0.4239897218739841
        U( 84)=          0.4181985492965744
        U( 85)=          0.4126461990878886
        U( 86)=          0.4073570798427895
        U( 87)=          0.4023461172079228
        U( 88)=          0.3976355585696443
        U( 89)=          0.3932377816410258
        U( 90)=          0.3891727128876379
        U( 91)=          0.3854500131142883
        U( 92)=          0.3820871425096385
        U( 93)=          0.3790910321656264
        U( 94)=          0.3764768199547278
        U( 95)=          0.3742482756438669
        U( 96)=          0.3724175838188968
        U( 97)=          0.3709860746618888
        U( 98)=          0.3699633032242018
        U( 99)=          0.3695899829574566
        U(100)=          0.3691382918814347
      
      end subroutine test_data
      
      
      
      subroutine test1(increasing)
         logical, intent(in) :: increasing
         integer, parameter :: n = 100
         real :: U(n)
         
         integer, parameter :: nvals = 4
         integer, parameter :: nwork = max(pm_work_size, mp_work_size)
         real, target :: work_ary(n*nwork)
         real, pointer :: work(:)
         real :: init_xs(n), xs(nvals), vals(nvals), dx
         integer :: i, ierr

         real, target :: f_ary(4*n)
         real, pointer :: f1(:), f(:,:)
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
         call interp_m3q_sg(init_xs, n, f1, nwork, work, 'test1', ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         write(*, 1) 'f(2, 1)', f(2, 1)
         write(*, 1) 'f(2, 2)', f(2, 2)
         write(*, 1) 'f(2, 3)', f(2, 3)
         write(*, *)
         
         xs = (/ 10d0, 10.4d0, 10.8d0, 11d0 /)
         if (.not. increasing) xs(:) = -xs(:)
         call interp_values_sg(init_xs, n, f1, nvals, xs, vals, ierr)
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
         real :: U(n)
         
         integer, parameter :: nvals = 7
         real :: init_xs(n), xs(nvals), vals(nvals), dx
         integer, parameter :: nwork = max(pm_work_size, mp_work_size)
         real, target :: work_ary(n*nwork)
         real, pointer :: work(:)
         integer :: i, ierr
         
         logical, parameter :: show_errors = .false.

         real, target :: f_ary(4*n)
         real, pointer :: f1(:), f(:,:)
         f1 => f_ary
         f(1:4,1:n) => f1(1:4*n)
         
         include 'formats'
        

         write(*, *) 'test2', increasing
         
         work => work_ary
         dx = 1d0 / (n-1)
         do i = 1, n
            init_xs(i) = 4*dx*(i-1)
            if (.not. increasing) init_xs(i) = -init_xs(i)
            U(i) = real(sin(dble(init_xs(i))))
         end do
         
         f(1, 1:n) = U(1:n)
         call interp_m3q_sg(init_xs, n, f1, nwork, work, 'test2', ierr)
         !call interp_m3q_on_uniform_grid(dx, n, f, ierr)
         if (ierr /= 0) then
            call mesa_error(__FILE__,__LINE__)
         end if
         write(*, *)
         
         xs = (/ 0d0, pi/4, pi/3, pi/2, 2*pi/3, 3*pi/4, 19*pi/20 /)
         if (.not. increasing) xs(:) = -xs(:)
         call interp_values_sg(init_xs, n, f1, nvals, xs, vals, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         
         if (show_errors) then
            do i = 1, nvals
               write(*, 2) 'x val exact err', 
     >               i, xs(i), vals(i), real(sin(dble(xs(i)))), vals(i) - real(sin(dble(xs(i))))
            end do
         else
            do i = 1, nvals
               write(*, 2) 'x val exact', i, xs(i), vals(i), real(sin(dble(xs(i))))
            end do
         end if
         write(*, *)
      
         call integrate_values_sg(init_xs, n, f1, nvals, xs, vals, ierr)
         if (ierr /= 0) then
            call mesa_error(__FILE__,__LINE__)
         end if

         if (show_errors) then
            do i = 2, nvals
               write(*, 2) 'x val exact err', i, xs(i), vals(i), 
     >               real(cos(dble(xs(i-1)))) - real(cos(dble(xs(i)))), 
     >                     vals(i) - (real(cos(dble(xs(i-1)))) - real(cos(dble(xs(i)))))
            end do
         else
            do i = 2, nvals
               write(*, 2) 'x val exact', i, xs(i), vals(i), 
     >               real(cos(dble(xs(i-1)))) - real(cos(dble(xs(i))))
            end do
         end if
         write(*, *)
      
      end subroutine test2
      
      
      subroutine test3(increasing)
         logical, intent(in) :: increasing
         integer, parameter :: n = 100
         real :: U(n)
         
         integer, parameter :: nvals = 3
         integer, parameter :: nwork = max(pm_work_size, mp_work_size)
         real, target :: work_ary(n*nwork)
         real, pointer :: work(:)
         real :: init_xs(n), xs(nvals), vals(nvals), dx
         integer :: i, ierr

         real, target :: f_ary(4*n)
         real, pointer :: f1(:), f(:,:)
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
         call interp_pm_sg(init_xs, n, f1, nwork, work, 'test3', ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         write(*, 1) 'f(2, 1)', f(2, 1)
         write(*, 1) 'f(2, 2)', f(2, 2)
         write(*, 1) 'f(2, 3)', f(2, 3)
         write(*, *)
         
         xs = (/ 10d0, 10.4d0, 10.8d0 /)
         if (.not. increasing) xs(:) = -xs(:)
         call interp_values_sg(init_xs, n, f1, nvals, xs, vals, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         write(*, 1) 'z(10.0)', vals(1)
         write(*, 1) 'z(10.4)', vals(2)
         write(*, 1) 'z(10.8)', vals(3)
         write(*, *)
      
      end subroutine test3
      
      
      subroutine test_min(increasing)
         logical, intent(in) :: increasing
         integer, parameter :: n = 100
         real :: U(n)
         
         integer, parameter :: nvals = 2
         integer, parameter :: nwork = max(pm_work_size, mp_work_size)
         real, target :: work_ary(n*nwork)
         real, pointer :: work(:)
         real :: init_xs(n), xs(nvals), vals(nvals), dx
         integer :: i, ierr

         real, target :: f_ary(4*n)
         real, pointer :: f1(:), f(:,:)
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
         call interp_pm_sg(init_xs, n, f1, nwork, work, 'test_min', ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         write(*, 1) 'f(2, 1)', f(2, 1)
         write(*, 1) 'f(2, 2)', f(2, 2)
         write(*, *)
         
         xs = (/ 10d0, 10.8d0 /)
         if (.not. increasing) xs(:) = -xs(:)
         call interp_values_sg(init_xs, n, f1, nvals, xs, vals, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         write(*, 1) 'z(10.0)', vals(1)
         write(*, 1) 'z(10.8)', vals(2)
         write(*, *)
      
      end subroutine test_min
      
      
      subroutine test_interp_3_to_2
         real :: pdqm1, pdq00, ndqm1, ndq00, pfm1, pf00, pfp1, nf00, nfp1
         integer :: ierr
         
         include 'formats'
         write(*, *) 'test_interp_3_to_2'
         pdqm1  =  2.0114947208182310E-03
         pdq00  =  2.0097307373083619E-03
         ndqm1  =  1.5784933404892154E-03
         ndq00  =  1.3270582904786616E-03
         pfm1  =  1.8984458478374221E+30
         pf00  =  1.4579738233866260E+30
         pfp1  =  1.1136297079131214E+30
         call interp_3_to_2_sg(
     >            pdqm1, pdq00, ndqm1, ndq00, pfm1, pf00, pfp1, nf00, nfp1, 'test_interp_3_to_2', ierr)
         write(*,*) 'ierr', ierr
         write(*,1) 'nf00', nf00
         write(*,1) 'nfp1', nfp1
      end subroutine test_interp_3_to_2
      
      
      subroutine test_4pt
         integer, parameter :: n = 6
         real :: x(n), y(n), a(3), dx, result, exact, f(4, n)
         integer :: ierr, i
         integer, parameter :: nwork = max(pm_work_size, mp_work_size)
         real, target :: work_ary(n*nwork)
         real, pointer :: work(:)
         
         include 'formats'
         work => work_ary
         x = (/ pi/4, pi/3, pi/2, 2*pi/3, 3*pi/4, 19*pi/20 /)
         do i=1, n
            y(i) = real(sin(dble(x(i))))
         end do
         call interp_4pt_pm_sg(x(2:5), y(2:5), a)
         dx = (x(4)-x(3))/2
         result = y(3) + dx*(a(1) + dx*(a(2) + dx*a(3)))
         exact = real(sin(dble(x(3)+dx)))
         write(*, *) 'test_4pt'
         write(*, 1) 'res exact err', result, exact, abs((result-exact)/exact)      
      end subroutine test_4pt
      

      end module interp_1d_support_sg




