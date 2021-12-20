! ***********************************************************************
!
!   Copyright (C) 2020  The MESA Team
!
!   This file is part of MESA.
!
!   MESA is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

      module create_FSCRliq8_table
      
      use const_def
      use chem_def
      use utils_lib, only: is_bad
      use math_lib
      
      implicit none

      public :: do_create_FSCRliq8_table
      private

      
   contains
   
   
   subroutine do_create_FSCRliq8_table(fname,iZion)
      character (len=*), intent(in) :: fname
      integer, intent(in) :: iZion
      real(dp) :: logRS, logRS_min, logRS_max, dlogRS
      real(dp) :: logGAME, logGAME_min, logGAME_max, dlogGAME
      integer :: nlogRS, nlogGAME, io_unit, i, j, ierr
      real(dp) :: Zion, RS, GAME, FSCR, USCR, PSCR, CVSCR, PDTSCR, PDRSCR
      include 'formats'
      
      Zion = iZion
      
      logRS_min = -3.5d0
      logRS_max = 0.0d0
      dlogRS = 1d-2
      nlogRS = (logRS_max - logRS_min)/dlogRS + 1
      
      logGAME_min = -2d0
      logGAME_max = 4d0
      dlogGAME = 1d-2
      nlogGAME = (logGAME_max - logGAME_min)/dlogGAME + 1

      io_unit = 40
         
      !write(*,'(a)') 'create ' // trim(fname)

      open(unit=io_unit,file=trim(fname))
      
      write(io_unit,'(99(a14))') 'num logRS', 'logRS min', 'logRS max', 'del logRS', &
      		'num logGAME', 'logGAME min', 'logGAME max', 'del logGAME', 'Zion'
      
      write(io_unit,'(2(i10,4x,3(f14.4)),i14)') &
         nlogRS, logRS_min, logRS_max, dlogRS, &
         nlogGAME, logGAME_min, logGAME_max, dlogGAME, iZion

      do i = 1, nlogRS
         logRS = logRS_min + (i-1) * dlogRS
         RS = exp10(logRS)
         write(io_unit,'(/,7x,a)') 'logRS' 
         write(io_unit,'(2x,f14.6/)') logRS
         write(io_unit,'(99(a26))') &
            'logGAME', 'FSCR', 'USCR', 'PSCR', 'CVSCR', 'PDTSCR', 'PDRSCR', 'RS', 'GAME'
         do j = 1, nlogGAME
            logGAME = logGAME_min + (j-1) * dlogGAME
            GAME = exp10(logGAME)
            if (is_bad(GAME)) then
               write(*,3) 'GAME', i, j, GAME, logGAME, logGAME_min, dlogGAME
            end if
            ierr = 0
            call FSCRliq8(RS, GAME, Zion, FSCR, USCR, PSCR, CVSCR, PDTSCR, PDRSCR, ierr)
            if (ierr /= 0) then
               write(*,*) 'FSCRliq8 error'
               write(*,1) 'RS', RS
               write(*,1) 'GAME', GAME
               write(*,1) 'Zion', Zion
               call mesa_error(__FILE__,__LINE__,'do_create_FSCRliq8_table')
            end if
            write(io_unit,'(99(1pd26.16))') &
               logGAME, FSCR, USCR, PSCR, CVSCR, PDTSCR, PDRSCR, RS, GAME
         end do
         write(io_unit,*)
      end do

      write(io_unit,*)
      write(io_unit,*)
      close(io_unit)
         
   end subroutine do_create_FSCRliq8_table


      subroutine FSCRliq8(RS,GAME,Zion, &
           FSCR,USCR,PSCR,CVSCR,PDTSCR,PDRSCR,ierr) ! fit to the el.-ion scr.
!                                                       Version 11.09.08
!                                                       cleaned 16.06.09
! Stems from FSCRliq7 v. 09.06.07. Included a check for RS=0.
!   INPUT: RS - density parameter, GAME - electron Coulomb parameter,
!          Zion - ion charge number,
!   OUTPUT: FSCR - screening (e-i) free energy per kT per 1 ion,
!           USCR - internal energy per kT per 1 ion (screen.contrib.)
!           PSCR - pressure divided by (n_i kT) (screen.contrib.)
!           CVSCR - heat capacity per 1 ion (screen.contrib.)
!           PDTSCR,PDRSCR = PSCR + d PSCR / d ln(T,\rho)
      real(dp), intent(in) :: RS,GAME
      real(dp), intent(in) :: Zion
      real(dp), intent(out) :: FSCR,USCR,PSCR,CVSCR,PDTSCR,PDRSCR
      integer, intent(out) :: ierr

      real(dp) :: SQG, SQR, SQZ1, SQZ, CDH0, CDH, ZLN, Z13
      real(dp) :: X, CTF, P01, P03, PTX
      real(dp) :: TX, TXDG, TXDGG, TY1, TY1DG, TY1DGG
      real(dp) :: TY2, TY2DX, TY2DXX, TY, TYX, TYDX, TYDG, P1
      real(dp) :: COR1, COR1DX, COR1DG, COR1DXX, COR1DGG, COR1DXG
      real(dp) :: U0, U0DX, U0DG, U0DXX, U0DGG, U0DXG
      real(dp) :: D0DG, D0, D0DX, D0DXX
      real(dp) :: COR0, COR0DX, COR0DG, COR0DXX, COR0DGG, COR0DXG
      real(dp) :: RELE, Q1, Q2, H1U, H1D, H1, H1X, H1DX, H1DXX
      real(dp) :: UP, UPDX, UPDG, UPDXX, UPDGG, UPDXG
      real(dp) :: DN1, DN1DX, DN1DG, DN1DXX, DN1DGG, DN1DXG
      real(dp) :: DN, DNDX, DNDG, DNDXX, DNDGG, DNDXG
      real(dp) :: FX, FXDG, FDX, FG, FDG, FDGDH, FDXX, FDGG, FDXG
      
      real(dp), parameter :: XRS=.0140047d0
      real(dp), parameter :: TINY=1.d-19

      ierr = 0
      if (RS.lt.0d0) then
         ierr = -1
         return
         !call mesa_error(__FILE__,__LINE__,'FSCRliq8: RS < 0')
      end if
      if (RS.lt.TINY) then
         FSCR=0.d0
         USCR=0.d0
         PSCR=0.d0
         CVSCR=0.d0
         PDTSCR=0.d0
         PDRSCR=0.d0
         return
      endif
      SQG=sqrt(GAME)
      SQR=sqrt(RS)
      SQZ1=sqrt(1d0+Zion)
      SQZ=sqrt(Zion)
      CDH0=Zion/1.73205d0 ! 1.73205=sqrt(3.)
      CDH=CDH0*(SQZ1*SQZ1*SQZ1-SQZ*SQZ*SQZ-1d0)
      SQG=sqrt(GAME)
      ZLN=log(Zion)
      Z13=exp(ZLN/3.d0) ! Zion**(1./3.)
      X=XRS/RS ! relativity parameter
      CTF=Zion*Zion*.2513d0*(Z13-1d0+.2d0/sqrt(Z13))
! Thomas-Fermi constant; .2513=(18/175)(12/\pi)^{2/3}
      P01=1.11d0*exp(0.475d0*ZLN)
      P03=0.2d0+0.078d0*ZLN*ZLN
      PTX=1.16d0+0.08d0*ZLN
      TX=pow(GAME,PTX)
      TXDG=PTX*TX/GAME
      TXDGG=(PTX-1.d0)*TXDG/GAME
      TY1=1d0/(1.d-3*Zion*Zion+2d0*GAME)
      TY1DG=-2d0*TY1*TY1
      TY1DGG=-4d0*TY1*TY1DG
      TY2=1d0+6d0*RS*RS
      TY2DX=-12d0*RS*RS/X
      TY2DXX=-3d0*TY2DX/X
      TY=RS*RS*RS/TY2*(1d0+TY1)
      TYX=3d0/X+TY2DX/TY2
      TYDX=-TY*TYX
      TYDG=RS*RS*RS*TY1DG/TY2
      P1=(Zion-1d0)/9d0
      COR1=1d0+P1*TY
      COR1DX=P1*TYDX
      COR1DG=P1*TYDG
      COR1DXX=P1*(TY*(3d0/(X*X)+(TY2DX/TY2)*(TY2DX/TY2)-TY2DXX/TY2)-TYDX*TYX)
      COR1DGG=P1*RS*RS*RS*TY1DGG/TY2
      COR1DXG=-P1*TYDG*TYX
      U0=0.78d0*sqrt(GAME/Zion)*RS*RS*RS
      U0DX=-3d0*U0/X
      U0DG=.5d0*U0/GAME
      U0DXX=-4.d0*U0DX/X
      U0DGG=-.5d0*U0DG/GAME
      U0DXG=-3.d0*U0DG/X
      D0DG=Zion*Zion*Zion
      D0=GAME*D0DG+21d0*RS*RS*RS
      D0DX=-63d0*RS*RS*RS/X
      D0DXX=252d0*RS*RS*RS/(X*X)
      COR0=1d0+U0/D0
      COR0DX=(U0DX-U0*D0DX/D0)/D0
      COR0DG=(U0DG-U0*D0DG/D0)/D0
      COR0DXX=(U0DXX-(2d0*U0DX*D0DX+U0*D0DXX)/D0+2d0*(D0DX/D0)*(D0DX/D0))/D0
      COR0DGG=(U0DGG-2d0*U0DG*D0DG/D0+2d0*U0*(D0DG/D0)*(D0DG/D0))/D0
      COR0DXG=(U0DXG-(U0DX*D0DG+U0DG*D0DX)/D0+2d0*U0*D0DX*D0DG/(D0*D0))/D0
! Relativism:
      RELE=sqrt(1.d0+X*X)
      Q1=0.18d0/sqrt(sqrt(Zion))
      Q2=0.2d0+0.37d0/sqrt(Zion)
      H1U=1d0+X*X/5.d0
      H1D=1d0+Q1*X+Q2*X*X
      H1=H1U/H1D
      H1X=0.4d0*X/H1U-(Q1+2d0*Q2*X)/H1D
      H1DX=H1*H1X
      H1DXX=H1DX*H1X &
         + H1*(0.4d0/H1U-(0.4d0*X/H1U)*(0.4d0*X/H1U)-2d0*Q2/H1D+pow2((Q1+2d0*Q2*X)/H1D))
      UP=CDH*SQG+P01*CTF*TX*COR0*H1
      UPDX=P01*CTF*TX*(COR0DX*H1+COR0*H1DX)
      UPDG=.5d0*CDH/SQG+P01*CTF*(TXDG*COR0+TX*COR0DG)*H1
      UPDXX=P01*CTF*TX*(COR0DXX*H1+2d0*COR0DX*H1DX+COR0*H1DXX)
      UPDGG=-.25d0*CDH/(SQG*GAME) &
         + P01*CTF*(TXDGG*COR0+2d0*TXDG*COR0DG+TX*COR0DGG)*H1
      UPDXG=P01*CTF*(TXDG*(COR0DX*H1+COR0*H1DX) &
        + TX*(COR0DXG*H1+COR0DG*H1DX))
      DN1=P03*SQG+P01/RS*TX*COR1
      DN1DX=P01*TX*(COR1/XRS+COR1DX/RS)
      DN1DG=.5d0*P03/SQG+P01/RS*(TXDG*COR1+TX*COR1DG)
      DN1DXX=P01*TX/XRS*(2d0*COR1DX+X*COR1DXX)
      DN1DGG=-.25d0*P03/(GAME*SQG) &
         + P01/RS*(TXDGG*COR1+2d0*TXDG*COR1DG+TX*COR1DGG)
      DN1DXG=P01*(TXDG*(COR1/XRS+COR1DX/RS)+TX*(COR1DG/XRS+COR1DXG/RS))
      DN=1d0+DN1/RELE
      DNDX=DN1DX/RELE-X*DN1/(RELE*RELE*RELE)
      DNDXX=(DN1DXX-((2d0*X*DN1DX+DN1)-3.d0*X*X*DN1/(RELE*RELE))/(RELE*RELE))/RELE
      DNDG=DN1DG/RELE
      DNDGG=DN1DGG/RELE
      DNDXG=DN1DXG/RELE-X*DN1DG/(RELE*RELE*RELE)
      FSCR=-UP/DN*GAME
      FX=(UP*DNDX/DN-UPDX)/DN
      FXDG=((UPDG*DNDX+UPDX*DNDG+UP*DNDXG-2d0*UP*DNDX*DNDG/DN)/DN-UPDXG)/DN
      FDX=FX*GAME
      FG=(UP*DNDG/DN-UPDG)/DN
      FDG=FG*GAME-UP/DN
      FDGDH=SQG*DNDG/(DN*DN) ! d FDG / d CDH
      FDXX=((UP*DNDXX+2d0*(UPDX*DNDX-UP*DNDX*DNDX/DN))/DN-UPDXX)/DN*GAME
      FDGG=2d0*FG+GAME*((2d0*DNDG*(UPDG-UP*DNDG/DN)+UP*DNDGG)/DN-UPDGG)/DN
      FDXG=FX+GAME*FXDG
      USCR=GAME*FDG
      CVSCR=-GAME*GAME*FDGG
      PSCR=(X*FDX+GAME*FDG)/3.d0
      PDTSCR=-GAME*GAME*(X*FXDG+FDGG)/3.d0
      PDRSCR=(12d0*PSCR+X*X*FDXX+2d0*X*GAME*FDXG+GAME*GAME*FDGG)/9d0
      return
      end subroutine FSCRliq8


      end module create_FSCRliq8_table
