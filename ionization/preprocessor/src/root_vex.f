!
! root_vex
! francois hebert, jul 17 2010
!
! subroutine to converge on the correct central potential value
!

      module root_vex

      use def_args
      use def_type, only: dp
      use mod_fast_tfdh

      implicit none

      contains


      subroutine rootfind_vex(args)

         type(datastruct), intent(inout) :: args

         integer, parameter :: max_iter = 100
         real (dp), parameter :: eps = 1.0e-12_dp

         integer :: i
         logical :: flag
         real (dp) :: v0, v1, vr, dv

         dv = 5.0_dp
         v0 = 0.0_dp
         v1 = 0.0_dp

         ! search for v0,v1 pair that bracket the root
         do, i=1, max_iter
            args % dv0 = v0

            flag = fasttfdh(args)

            if (.not. flag) exit
            v1 = v0
            v0 = v0 - dv
         end do

         ! now close in on root
         do, i=1, max_iter
            vr = (v0 + v1)/2.0_dp
            if (vr==v0 .or. vr==v1) exit !underflow
            args % dv0 = vr

            flag = fasttfdh(args)

            if (flag) then
               v1 = vr
            else
               v0 = vr
            end if
            if (abs((v1-v0)/v1) < eps) exit
         end do

         args % dv0 = v1

      end subroutine rootfind_vex


      end module root_vex

   