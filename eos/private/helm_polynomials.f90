module helm_polynomials
   use const_def, only : dp
   
   implicit none

contains
   
   !..quintic hermite polynomial statement functions
   !..psi0 and its derivatives
   pure real(dp) function psi0(z)
      real(dp), intent(in) :: z
      psi0 = z * z * z * (z * (-6.0d0 * z + 15.0d0) - 10.0d0) + 1.0d0
   end function psi0
   
   pure real(dp) function dpsi0(z)
      real(dp), intent(in) :: z
      dpsi0 = z * z * (z * (-30.0d0 * z + 60.0d0) - 30.0d0)
   end function dpsi0
   
   pure real(dp) function ddpsi0(z)
      real(dp), intent(in) :: z
      ddpsi0 = z * (z * (-120.0d0 * z + 180.0d0) - 60.0d0)
   end function ddpsi0
   
   pure real(dp) function dddpsi0(z)
      real(dp), intent(in) :: z
      dddpsi0 = z * (-360.0d0 * z + 360.0d0) - 60.0d0
   end function dddpsi0
   
   
   !..psi1 and its derivatives
   
   pure real(dp) function psi1(z)
      real(dp), intent(in) :: z
      psi1 = z * (z * z * (z * (-3.0d0 * z + 8.0d0) - 6.0d0) + 1.0d0)
   end function psi1
   
   pure real(dp) function dpsi1(z)
      real(dp), intent(in) :: z
      dpsi1 = z * z * (z * (-15.0d0 * z + 32.0d0) - 18.0d0) + 1.0d0
   end function dpsi1
   
   pure real(dp) function ddpsi1(z)
      real(dp), intent(in) :: z
      ddpsi1 = z * (z * (-60.0d0 * z + 96.0d0) - 36.0d0)
   end function ddpsi1
   
   pure real(dp) function dddpsi1(z)
      real(dp), intent(in) :: z
      dddpsi1 = z * (-180.0d0 * z + 192.0d0) - 36.0d0
   end function dddpsi1
   
   !..psi2  and its derivatives
   
   pure real(dp) function psi2(z)
      real(dp), intent(in) :: z
      psi2 = 0.5d0 * z * z * (z * (z * (-z + 3.0d0) - 3.0d0) + 1.0d0)
   end function psi2
   
   pure real(dp) function dpsi2(z)
      real(dp), intent(in) :: z
      dpsi2 = 0.5d0 * z * (z * (z * (-5.0d0 * z + 12.0d0) - 9.0d0) + 2.0d0)
   end function dpsi2
   
   pure real(dp) function ddpsi2(z)
      real(dp), intent(in) :: z
      ddpsi2 = 0.5d0 * (z * (z * (-20.0d0 * z + 36.0d0) - 18.0d0) + 2.0d0)
   end function ddpsi2
   
   pure real(dp) function dddpsi2(z)
      real(dp), intent(in) :: z
      dddpsi2 = 0.5d0 * (z * (-60.0d0 * z + 72.0d0) - 18.0d0)
   end function dddpsi2
   
   !..biquintic hermite polynomial statement function
   
   pure real(dp) function h5(i, j, fi, &
      w0t, w1t, w2t, w0mt, w1mt, w2mt, &
      w0d, w1d, w2d, w0md, w1md, w2md)
      integer, intent(in) :: i, j
      real(dp), intent(in) :: fi(36)
      real(dp), intent(in) :: w0t, w1t, w2t, w0mt, w1mt, w2mt
      real(dp), intent(in) :: w0d, w1d, w2d, w0md, w1md, w2md
      
      h5 = fi(1) * w0d * w0t + fi(2) * w0md * w0t &
         + fi(3) * w0d * w0mt + fi(4) * w0md * w0mt &
         + fi(5) * w0d * w1t + fi(6) * w0md * w1t &
         + fi(7) * w0d * w1mt + fi(8) * w0md * w1mt &
         + fi(9) * w0d * w2t + fi(10) * w0md * w2t &
         + fi(11) * w0d * w2mt + fi(12) * w0md * w2mt &
         + fi(13) * w1d * w0t + fi(14) * w1md * w0t &
         + fi(15) * w1d * w0mt + fi(16) * w1md * w0mt &
         + fi(17) * w2d * w0t + fi(18) * w2md * w0t &
         + fi(19) * w2d * w0mt + fi(20) * w2md * w0mt &
         + fi(21) * w1d * w1t + fi(22) * w1md * w1t &
         + fi(23) * w1d * w1mt + fi(24) * w1md * w1mt &
         + fi(25) * w2d * w1t + fi(26) * w2md * w1t &
         + fi(27) * w2d * w1mt + fi(28) * w2md * w1mt &
         + fi(29) * w1d * w2t + fi(30) * w1md * w2t &
         + fi(31) * w1d * w2mt + fi(32) * w1md * w2mt &
         + fi(33) * w2d * w2t + fi(34) * w2md * w2t &
         + fi(35) * w2d * w2mt + fi(36) * w2md * w2mt
   end function h5
   
   !..cubic hermite polynomial statement functions
   !..psi0 & derivatives
   
   pure real(dp) function xpsi0(z)
      real(dp), intent(in) :: z
      xpsi0 = z * z * (2.0d0 * z - 3.0d0) + 1.0d0
   end function xpsi0
   
   pure real(dp) function xdpsi0(z)
      real(dp), intent(in) :: z
      xdpsi0 = z * (6.0d0 * z - 6.0d0)
   end function xdpsi0
   
   pure real(dp) function xddpsi0(z)
      real(dp), intent(in) :: z
      xddpsi0 = 12.0d0 * z - 6.0d0
   end function xddpsi0
   
   pure real(dp) function xdddpsi0(z)
      real(dp), intent(in) :: z
      xdddpsi0 = 12.0d0
   end function xdddpsi0
   
   !..psi1 & derivatives
   
   pure real(dp) function xpsi1(z)
      real(dp), intent(in) :: z
      xpsi1 = z * (z * (z - 2.0d0) + 1.0d0)
   end function xpsi1
   
   pure real(dp) function xdpsi1(z)
      real(dp), intent(in) :: z
      xdpsi1 = z * (3.0d0 * z - 4.0d0) + 1.0d0
   end function xdpsi1
   
   pure real(dp) function xddpsi1(z)
      real(dp), intent(in) :: z
      xddpsi1 = 6.0d0 * z - 4.0d0
   end function xddpsi1
   
   pure real(dp) function xdddpsi1(z)
      real(dp), intent(in) :: z
      xdddpsi1 = 6.0d0
   end function xdddpsi1
   
   !..bicubic hermite polynomial statement function
   pure real(dp) function h3(i, j, fi, &
      w0t, w1t, w0mt, w1mt, w0d, w1d, w0md, w1md)
      integer, intent(in) :: i, j
      real(dp), intent(in) :: fi(36)
      real(dp), intent(in) :: w0t, w1t, w0mt, w1mt, w0d, w1d, w0md, w1md
      h3 = fi(1) * w0d * w0t + fi(2) * w0md * w0t  &
         + fi(3) * w0d * w0mt + fi(4) * w0md * w0mt &
         + fi(5) * w0d * w1t + fi(6) * w0md * w1t  &
         + fi(7) * w0d * w1mt + fi(8) * w0md * w1mt &
         + fi(9) * w1d * w0t + fi(10) * w1md * w0t  &
         + fi(11) * w1d * w0mt + fi(12) * w1md * w0mt &
         + fi(13) * w1d * w1t + fi(14) * w1md * w1t  &
         + fi(15) * w1d * w1mt + fi(16) * w1md * w1mt
   end function h3

end module helm_polynomials
