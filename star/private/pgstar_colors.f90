! ***********************************************************************
!
!   Copyright (C) 2010-2025  The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! ***********************************************************************

module pgstar_colors
   implicit none
   public

   integer, parameter :: clr_Background = 0
   integer, parameter :: clr_Foreground = 1
   ! Values are set for backwards compat
   integer, private, parameter :: clr_RedAlt = 2
   integer, private, parameter :: clr_GreenAlt = 3
   integer, private, parameter :: clr_BlueAlt = 4
   integer, parameter :: clr_Cyan = 5
   integer, parameter :: clr_Magenta = 6
   integer, parameter :: clr_Yellow = 7
   integer, parameter :: clr_Orange = 8
   integer, parameter :: clr_LimeGreen = 9
   integer, parameter :: clr_GreenYellow = 10
   integer, parameter :: clr_DodgerBlue = 11
   integer, parameter :: clr_MagentaDark = 12
   integer, parameter :: clr_Plum = 13
   integer, parameter :: clr_SandyBrown = 14
   integer, parameter :: clr_Salmon = 15
   integer, parameter :: clr_Grey59 = 16
   integer, parameter :: clr_Grey30 = 17
   integer, parameter :: clr_Black = 18
   integer, parameter :: clr_Blue = 19
   integer, parameter :: clr_BrightBlue = 20
   integer, parameter :: clr_Goldenrod = 21
   integer, parameter :: clr_Lilac = 22
   integer, parameter :: clr_Coral = 23
   integer, parameter :: clr_FireBrick = 24
   integer, parameter :: clr_RoyalPurple = 25
   integer, parameter :: clr_Gold = 26
   integer, parameter :: clr_Crimson = 27
   integer, parameter :: clr_SlateGray = 28
   integer, parameter :: clr_Teal = 29
   integer, parameter :: clr_LightSteelBlue = 30
   integer, parameter :: clr_MediumSlateBlue = 31
   integer, parameter :: clr_MediumSpringGreen = 32
   integer, parameter :: clr_MediumBlue = 33
   integer, parameter :: clr_RoyalBlue = 34
   integer, parameter :: clr_LightGray = 35
   integer, parameter :: clr_Silver = 36
   integer, parameter :: clr_DarkGray = 37
   integer, parameter :: clr_Gray = 38
   integer, parameter :: clr_LightSkyBlue = 39
   integer, parameter :: clr_LightSkyGreen = 40
   integer, parameter :: clr_SeaGreen = 41
   integer, parameter :: clr_Tan = 42
   integer, parameter :: clr_IndianRed = 43
   integer, parameter :: clr_LightOliveGreen = 44
   integer, parameter :: clr_CadetBlue = 45
   integer, parameter :: clr_Beige = 46

   integer, parameter :: colormap_offset = 46
   integer, parameter :: colormap_length = 101

   integer, parameter :: clr_no_mixing = clr_SeaGreen
   integer, parameter :: clr_convection = clr_LightSkyBlue
   integer, parameter :: clr_leftover_convection = clr_BrightBlue
   integer, parameter :: clr_semiconvection = clr_SlateGray
   integer, parameter :: clr_thermohaline = clr_Lilac
   integer, parameter :: clr_overshoot = clr_Beige
   integer, parameter :: clr_rotation = clr_LightSkyGreen
   integer, parameter :: clr_rayleigh_taylor = clr_IndianRed
   integer, parameter :: clr_minimum = clr_Coral
   integer, parameter :: clr_anonymous = clr_Tan
contains
   subroutine set_device_colors(white_on_black_flag)
      logical, intent(in) :: white_on_black_flag
      integer :: i
      real, dimension(3, colormap_length) :: colormap

      if (white_on_black_flag) then
         call pgscr(clr_Background, 0., 0., 0.)
         call pgscr(clr_Foreground, 1., 1., 1.)
      else
         call pgscr(clr_Foreground, 0., 0., 0.)
         call pgscr(clr_Background, 1., 1., 1.)
      end if
      ! These are duplicated later, but we keep them for backwards compatibility
      call pgscr(clr_RedAlt, 1., 0., 0.)
      call pgscr(clr_GreenAlt, 0., 1., 0.)
      call pgscr(clr_BlueAlt, 0., 0., 1.)

      call pgscr(clr_Cyan, 0., 1., 1.)
      call pgscr(clr_Magenta, 1., 0., 1.)
      call pgscr(clr_Yellow, 1., 1., 0.)
      call pgscr(clr_Orange, 1., 0.65, 0.)
      call pgscr(clr_LimeGreen, 0.2, 0.8, 0.2)
      call pgscr(clr_GreenYellow, 0.68, 1., 0.18)
      call pgscr(clr_DodgerBlue, 0.12, 0.56, 1.0)
      call pgscr(clr_MagentaDark, 0.55, 0., 0.55)
      call pgscr(clr_Plum, 0.87, 0.63, 0.87)
      call pgscr(clr_SandyBrown, 0.96, 0.64, 0.38)
      call pgscr(clr_Salmon, 0.91, 0.59, 0.48)
      call pgscr(clr_Grey59, 0.59, 0.59, 0.59)
      call pgscr(clr_Grey30, 0.3, 0.3, 0.3)

      ! Tioga colors
      call pgscr(clr_Black, 0.0, 0.0, 0.0)
      call pgscr(clr_Blue, 0.0, 0.0, 1.0)
      call pgscr(clr_BrightBlue, 0.0, 0.4, 1.0)
      call pgscr(clr_LightSkyBlue, 0.53, 0.808, 0.98)
      call pgscr(clr_LightSkyGreen, 0.125, 0.698, 0.668)
      call pgscr(clr_MediumSpringGreen, 0.0, 0.98, 0.604)
      call pgscr(clr_Goldenrod, 0.855, 0.648, 0.125)
      call pgscr(clr_Lilac, 0.8, 0.6, 1.0)
      call pgscr(clr_Coral, 1.0, 0.498, 0.312)
      call pgscr(clr_FireBrick, 0.698, 0.132, 0.132)
      call pgscr(clr_RoyalPurple, 0.4, 0.0, 0.6)
      call pgscr(clr_Gold, 1.0, 0.844, 0.0)
      call pgscr(clr_Crimson, 0.8, 0.0, 0.2)
      call pgscr(clr_SlateGray, 0.44, 0.5, 0.565)
      call pgscr(clr_SeaGreen, 0.18, 0.545, 0.34)
      call pgscr(clr_Teal, 0.0, 0.5, 0.5)
      call pgscr(clr_LightSteelBlue, 0.69, 0.77, 0.87)
      call pgscr(clr_MediumSlateBlue, 0.484, 0.408, 0.932)
      call pgscr(clr_MediumBlue, 0.0, 0.0, 0.804)
      call pgscr(clr_RoyalBlue, 0.255, 0.41, 0.884)
      call pgscr(clr_LightGray, 0.828, 0.828, 0.828)
      call pgscr(clr_Silver, 0.752, 0.752, 0.752)
      call pgscr(clr_DarkGray, 0.664, 0.664, 0.664)
      call pgscr(clr_Gray, 0.5, 0.5, 0.5)
      call pgscr(clr_IndianRed, 0.804, 0.36, 0.36)
      call pgscr(clr_Tan, 0.824, 0.705, 0.55)
      call pgscr(clr_LightOliveGreen, 0.6, 0.8, 0.6)
      call pgscr(clr_CadetBlue, 0.372, 0.62, 0.628)
      call pgscr(clr_Beige, 0.96, 0.96, 0.864)

      colormap(1:3, 1) = [ 0.0, 0.0, 1.0 ]
      colormap(1:3, 2) = [ 0.0196078431372549, 0.0196078431372549, 1.0 ]
      colormap(1:3, 3) = [ 0.0352941176470588, 0.0352941176470588, 1.0 ]
      colormap(1:3, 4) = [ 0.0588235294117647, 0.0588235294117647, 1.0 ]
      colormap(1:3, 5) = [ 0.0705882352941176, 0.0705882352941176, 1.0 ]
      colormap(1:3, 6) = [ 0.0941176470588235, 0.0941176470588235, 1.0 ]
      colormap(1:3, 7) = [ 0.105882352941176, 0.105882352941176, 1.0 ]
      colormap(1:3, 8) = [ 0.129411764705882, 0.129411764705882, 1.0 ]
      colormap(1:3, 9) = [ 0.141176470588235, 0.141176470588235, 1.0 ]
      colormap(1:3, 10) = [ 0.164705882352941, 0.164705882352941, 1.0 ]
      colormap(1:3, 11) = [ 0.184313725490196, 0.184313725490196, 1.0 ]
      colormap(1:3, 12) = [ 0.2, 0.2, 1.0 ]
      colormap(1:3, 13) = [ 0.219607843137255, 0.219607843137255, 1.0 ]
      colormap(1:3, 14) = [ 0.235294117647059, 0.235294117647059, 1.0 ]
      colormap(1:3, 15) = [ 0.254901960784314, 0.254901960784314, 1.0 ]
      colormap(1:3, 16) = [ 0.270588235294118, 0.270588235294118, 1.0 ]
      colormap(1:3, 17) = [ 0.294117647058824, 0.294117647058824, 1.0 ]
      colormap(1:3, 18) = [ 0.305882352941176, 0.305882352941176, 1.0 ]
      colormap(1:3, 19) = [ 0.329411764705882, 0.329411764705882, 1.0 ]
      colormap(1:3, 20) = [ 0.341176470588235, 0.341176470588235, 1.0 ]
      colormap(1:3, 21) = [ 0.364705882352941, 0.364705882352941, 1.0 ]
      colormap(1:3, 22) = [ 0.384313725490196, 0.384313725490196, 1.0 ]
      colormap(1:3, 23) = [ 0.4, 0.4, 1.0 ]
      colormap(1:3, 24) = [ 0.419607843137255, 0.419607843137255, 1.0 ]
      colormap(1:3, 25) = [ 0.435294117647059, 0.435294117647059, 1.0 ]
      colormap(1:3, 26) = [ 0.454901960784314, 0.454901960784314, 1.0 ]
      colormap(1:3, 27) = [ 0.470588235294118, 0.470588235294118, 1.0 ]
      colormap(1:3, 28) = [ 0.490196078431373, 0.490196078431373, 1.0 ]
      colormap(1:3, 29) = [ 0.505882352941176, 0.505882352941176, 1.0 ]
      colormap(1:3, 30) = [ 0.529411764705882, 0.529411764705882, 1.0 ]
      colormap(1:3, 31) = [ 0.549019607843137, 0.549019607843137, 1.0 ]
      colormap(1:3, 32) = [ 0.564705882352941, 0.564705882352941, 1.0 ]
      colormap(1:3, 33) = [ 0.584313725490196, 0.584313725490196, 1.0 ]
      colormap(1:3, 34) = [ 0.6, 0.6, 1.0 ]
      colormap(1:3, 35) = [ 0.619607843137255, 0.619607843137255, 1.0 ]
      colormap(1:3, 36) = [ 0.635294117647059, 0.635294117647059, 1.0 ]
      colormap(1:3, 37) = [ 0.654901960784314, 0.654901960784314, 1.0 ]
      colormap(1:3, 38) = [ 0.670588235294118, 0.670588235294118, 1.0 ]
      colormap(1:3, 39) = [ 0.690196078431373, 0.690196078431373, 1.0 ]
      colormap(1:3, 40) = [ 0.705882352941177, 0.705882352941177, 1.0 ]
      colormap(1:3, 41) = [ 0.725490196078431, 0.725490196078431, 1.0 ]
      colormap(1:3, 42) = [ 0.749019607843137, 0.749019607843137, 1.0 ]
      colormap(1:3, 43) = [ 0.764705882352941, 0.764705882352941, 1.0 ]
      colormap(1:3, 44) = [ 0.784313725490196, 0.784313725490196, 1.0 ]
      colormap(1:3, 45) = [ 0.8, 0.8, 1.0 ]
      colormap(1:3, 46) = [ 0.831372549019608, 0.831372549019608, 1.0 ]
      colormap(1:3, 47) = [ 0.854901960784314, 0.854901960784314, 1.0 ]
      colormap(1:3, 48) = [ 0.890196078431372, 0.890196078431372, 1.0 ]
      colormap(1:3, 49) = [ 0.913725490196078, 0.913725490196078, 1.0 ]
      colormap(1:3, 50) = [ 0.949019607843137, 0.949019607843137, 1.0 ]
      colormap(1:3, 51) = [ 1.0, 0.972549019607843, 0.972549019607843 ]
      colormap(1:3, 52) = [ 1.0, 0.949019607843137, 0.949019607843137 ]
      colormap(1:3, 53) = [ 1.0, 0.913725490196078, 0.913725490196078 ]
      colormap(1:3, 54) = [ 1.0, 0.890196078431372, 0.890196078431372 ]
      colormap(1:3, 55) = [ 1.0, 0.854901960784314, 0.854901960784314 ]
      colormap(1:3, 56) = [ 1.0, 0.831372549019608, 0.831372549019608 ]
      colormap(1:3, 57) = [ 1.0, 0.8, 0.8 ]
      colormap(1:3, 58) = [ 1.0, 0.784313725490196, 0.784313725490196 ]
      colormap(1:3, 59) = [ 1.0, 0.764705882352941, 0.764705882352941 ]
      colormap(1:3, 60) = [ 1.0, 0.749019607843137, 0.749019607843137 ]
      colormap(1:3, 61) = [ 1.0, 0.725490196078431, 0.725490196078431 ]
      colormap(1:3, 62) = [ 1.0, 0.705882352941177, 0.705882352941177 ]
      colormap(1:3, 63) = [ 1.0, 0.690196078431373, 0.690196078431373 ]
      colormap(1:3, 64) = [ 1.0, 0.670588235294118, 0.670588235294118 ]
      colormap(1:3, 65) = [ 1.0, 0.654901960784314, 0.654901960784314 ]
      colormap(1:3, 66) = [ 1.0, 0.635294117647059, 0.635294117647059 ]
      colormap(1:3, 67) = [ 1.0, 0.619607843137255, 0.619607843137255 ]
      colormap(1:3, 68) = [ 1.0, 0.6, 0.6 ]
      colormap(1:3, 69) = [ 1.0, 0.584313725490196, 0.584313725490196 ]
      colormap(1:3, 70) = [ 1.0, 0.564705882352941, 0.564705882352941 ]
      colormap(1:3, 71) = [ 1.0, 0.541176470588235, 0.541176470588235 ]
      colormap(1:3, 72) = [ 1.0, 0.529411764705882, 0.529411764705882 ]
      colormap(1:3, 73) = [ 1.0, 0.505882352941176, 0.505882352941176 ]
      colormap(1:3, 74) = [ 1.0, 0.490196078431373, 0.490196078431373 ]
      colormap(1:3, 75) = [ 1.0, 0.470588235294118, 0.470588235294118 ]
      colormap(1:3, 76) = [ 1.0, 0.454901960784314, 0.454901960784314 ]
      colormap(1:3, 77) = [ 1.0, 0.435294117647059, 0.435294117647059 ]
      colormap(1:3, 78) = [ 1.0, 0.419607843137255, 0.419607843137255 ]
      colormap(1:3, 79) = [ 1.0, 0.4, 0.4 ]
      colormap(1:3, 80) = [ 1.0, 0.384313725490196, 0.384313725490196 ]
      colormap(1:3, 81) = [ 1.0, 0.364705882352941, 0.364705882352941 ]
      colormap(1:3, 82) = [ 1.0, 0.341176470588235, 0.341176470588235 ]
      colormap(1:3, 83) = [ 1.0, 0.329411764705882, 0.329411764705882 ]
      colormap(1:3, 84) = [ 1.0, 0.305882352941176, 0.305882352941176 ]
      colormap(1:3, 85) = [ 1.0, 0.294117647058824, 0.294117647058824 ]
      colormap(1:3, 86) = [ 1.0, 0.270588235294118, 0.270588235294118 ]
      colormap(1:3, 87) = [ 1.0, 0.254901960784314, 0.254901960784314 ]
      colormap(1:3, 88) = [ 1.0, 0.235294117647059, 0.235294117647059 ]
      colormap(1:3, 89) = [ 1.0, 0.219607843137255, 0.219607843137255 ]
      colormap(1:3, 90) = [ 1.0, 0.2, 0.2 ]
      colormap(1:3, 91) = [ 1.0, 0.176470588235294, 0.176470588235294 ]
      colormap(1:3, 92) = [ 1.0, 0.164705882352941, 0.164705882352941 ]
      colormap(1:3, 93) = [ 1.0, 0.141176470588235, 0.141176470588235 ]
      colormap(1:3, 94) = [ 1.0, 0.129411764705882, 0.129411764705882 ]
      colormap(1:3, 95) = [ 1.0, 0.105882352941176, 0.105882352941176 ]
      colormap(1:3, 96) = [ 1.0, 0.0941176470588235, 0.0941176470588235 ]
      colormap(1:3, 97) = [ 1.0, 0.0705882352941176, 0.0705882352941176 ]
      colormap(1:3, 98) = [ 1.0, 0.0588235294117647, 0.0588235294117647 ]
      colormap(1:3, 99) = [ 1.0, 0.0352941176470588, 0.0352941176470588 ]
      colormap(1:3, 100) = [ 1.0, 0.0196078431372549, 0.0196078431372549 ]
      colormap(1:3, 101) = [ 1.0, 0.0, 0.0 ]

      do i = 1, colormap_length
         call pgscr(colormap_offset + i, colormap(1, i), colormap(2, i), colormap(3, i))
      end do
   end subroutine set_device_colors
end module pgstar_colors

