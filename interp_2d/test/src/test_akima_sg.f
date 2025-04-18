      subroutine test_akima_sg
      use interp_2d_lib_sg

! Test for the RGBI3P_sg/RGSF3P_sg subroutine package

! Hiroshi Akima
! U.S. Department of Commerce, NTIA/ITS
! Version of 1995/08

! This program calls the RGBI3P_sg and RGSF3P_sg subroutines.

! This program requires no input data files.

! This program creates the TPRG3P data file.  All elements of
! the DZI array in the data file are expected to be zero.

! Specification statements
!     .. Parameters ..
      integer :: NXD,NYD,NXI
      real :: XIMN,XIMX
      integer :: NYI
      real :: YIMN,YIMX
      parameter (NXD=9,NYD=11,NXI=19,XIMN=-0.5,XIMX=8.5,NYI=23,YIMN=-0.5,YIMX=10.5)

!     .. Local Scalars ..
      real :: ANXIM1,ANYIM1,DXI,DYI
      integer :: IER,ISEC,IXD,IXI,IXIMN,IXIMX,IYD,IYDR,IYI,IYIR,MD,NYDO2
      character(len=9) :: NMPR
      character(len=6) :: NMWF

!     .. Local Arrays ..
      real :: DZI(NXI,NYI),WK(3,NXD,NYD),XD(NXD),XI(NXI),YD(NYD),YI(NYI),ZD(NXD,NYD),ZI(NXI,NYI),ZIE(NXI,NYI)
      character(len=9) :: NMSR(2)
      character(len=20) :: LBL(2)

!     .. Intrinsic Functions ..
      INTRINSIC MOD,real

! Data statements
      DATA   NMPR/'TPRG3P_sg'/,NMWF/'WFRG3P'/,NMSR/'RGBI3P_sg','RGSF3P_sg'/,LBL/'Calculated ZI Values','Differences         '/
      DATA   XD/0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0/
      DATA   YD/0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0/
      DATA   ZD/9*0.0,9*0.0,9*0.0,3.2,0.7,7*0.0,7.4,4.8,1.4,
     &       0.1,5*0.0,12.0,8.0,5.3,2.9,0.6,4*0.0,16.8,14.4,
     &       8.1,6.9,6.2,0.6,0.1,2*0.0,21.8,20.5,12.8,17.6,
     &       5.8,7.6,0.8,0.6,0.6,22.4,22.5,14.6,22.5,4.7,7.2,
     &       1.8,2.1,2.1,37.2,40.0,27.0,41.3,14.1,24.5,17.3,
     &       20.2,20.8,58.2,61.5,47.9,62.3,34.6,45.5,38.2,
     &       41.2,41.7/
      DATA   ((ZIE(IXI,IYI),IXI=1,NXI),IYI=1,5)/-.847,-.533,
     &       -.274,-.117,-.031,.000,.000,.000,.000,.000,.000,
     &       .000,.000,.000,.000,.000,.000,.000,.000,.000,
     &       .000,.000,.000,.000,.000,.000,.000,.000,.000,
     &       .000,.000,.000,.000,.000,.000,.000,.000,.000,
     &       .401,.250,.119,.043,.011,.000,.000,.000,.000,
     &       .000,.000,.000,.000,.000,.000,.000,.000,.000,
     &       .000,.000,.000,.000,.000,.000,.000,.000,.000,
     &       .000,.000,.000,.000,.000,.000,.000,.000,.000,
     &       .000,.000,-.665,-.376,-.143,-.033,-.007,.000,
     &       .000,.000,.000,.000,.000,.000,.000,.000,.000,
     &       .000,.000,.000,.000/
      DATA   ((ZIE(IXI,IYI),IXI=1,NXI),IYI=6,10)/.000,.000,
     &       .000,.000,.000,.000,.000,.000,.000,.000,.000,
     &       .000,.000,.000,.000,.000,.000,.000,.000,2.449,
     &       1.368,.537,.149,.025,.000,.000,.000,.000,.000,
     &       .000,.000,.000,.000,.000,.000,.000,.000,.000,
     &       5.083,3.200,1.642,.700,.187,.000,.000,.000,.000,
     &       .000,.000,.000,.000,.000,.000,.000,.000,.000,
     &       .000,6.588,5.234,3.878,2.542,1.188,.253,.026,
     &       .026,.007,.000,.000,.000,.000,.000,.000,.000,
     &       .000,.000,.000,8.017,7.400,6.400,4.800,2.963,
     &       1.400,.457,.100,.027,.000,.000,.000,.000,.000,
     &       .000,.000,.000,.000,.000/
      DATA   ((ZIE(IXI,IYI),IXI=1,NXI),IYI=11,15)/11.055,
     &       9.670,8.083,6.305,4.786,3.421,2.043,1.112,.565,
     &       .131,-.019,.000,.000,.000,.000,.000,.000,.000,
     &       .000,14.492,12.000,9.746,8.000,6.594,5.300,4.081,
     &       2.900,1.697,.600,.059,.000,.000,.000,.000,.000,
     &       .000,.000,.000,15.999,14.376,12.657,10.774,8.620,
     &       6.659,5.291,4.392,3.926,3.005,1.223,.139,.051,
     &       .025,.009,.000,.000,.000,-.005,15.525,16.800,
     &       16.749,14.400,10.956,8.100,6.735,6.900,7.298,
     &       6.200,3.010,.600,.248,.100,.024,.000,.006,.000,
     &       -.025,15.876,19.280,20.563,17.856,13.242,10.219,
     &       10.577,11.999,10.170,7.053,5.198,3.543,1.831,
     &       .350,-.130,.168,.408,.168,-.224/
      DATA   ((ZIE(IXI,IYI),IXI=1,NXI),IYI=16,20)/17.700,
     &       21.800,23.531,20.500,15.087,12.800,15.817,17.600,
     &       11.477,5.800,6.988,7.600,4.410,.800,-.392,.600,
     &       1.261,.600,-.417,17.913,22.788,24.944,21.881,
     &       16.302,14.382,18.557,20.807,11.916,4.561,7.327,
     &       8.518,5.133,1.284,-.013,1.201,1.998,1.200,-.065,
     &       16.383,22.400,25.330,22.500,16.796,14.600,19.172,
     &       22.500,13.159,4.700,6.689,7.200,4.392,1.800,
     &       1.150,2.100,2.734,2.100,1.025,18.109,26.756,
     &       31.311,28.143,21.004,18.237,24.236,28.979,17.970,
     &       7.469,10.467,11.985,9.022,6.833,6.901,8.292,
     &       9.186,8.524,7.101,24.667,37.200,44.007,40.000,
     &       30.508,27.000,34.974,41.300,27.136,14.100,20.473,
     &       24.500,20.557,17.300,17.639,20.200,21.826,20.800,
     &       18.458/
      DATA   ((ZIE(IXI,IYI),IXI=1,NXI),IYI=21,23)/33.414,
     &       48.009,56.017,51.561,40.817,36.922,45.856,52.860,
     &       37.376,23.200,30.839,36.192,31.969,28.037,28.437,
     &       31.604,33.579,32.332,29.561,44.842,58.200,65.537,
     &       61.500,51.657,47.900,55.899,62.300,47.891,34.600,
     &       41.239,45.500,41.479,38.200,38.591,41.200,42.823,
     &       41.700,39.192,58.284,68.917,74.644,71.333,63.413,
     &       60.125,66.293,71.400,59.129,47.725,52.451,54.592,
     &       50.842,48.483,48.639,50.142,51.089,50.200,48.268/

! Calculation
! Opens the output file and writes the input data.
!      OPEN (6,FILE=NMWF)
      NYDO2 = NYD/2
      write (6,FMT=9000) NMPR
      write (6,FMT=9010) XD
      do IYDR = 1,NYD
          if (MOD(IYDR-1,NYDO2) <= 1) write (6,FMT='(1X)')
          IYD = NYD + 1 - IYDR
          write (6,FMT=9020) YD(IYD), (ZD(IXD,IYD),IXD=1,NXD)
      end do
! Program check for the RGBI3P_sg subroutine
! - Performs interpolation and calculates the differences.
      DXI = XIMX - XIMN
      ANXIM1 = NXI - 1
      do IXI = 1,NXI
          XI(IXI) = XIMN + DXI*real(IXI-1)/ANXIM1
      end do
      DYI = YIMX - YIMN
      ANYIM1 = NYI - 1
      do IYI = 1,NYI
          YI(IYI) = YIMN + DYI*real(IYI-1)/ANYIM1
      end do
      do IYI = 1,NYI
          do IXI = 1,NXI
              if (IXI == 1 .and. IYI == 1) then
                  MD = 1
              else
                  MD = 2
              end if
              CALL interp_RGBI3P_sg(MD,NXD,NYD,XD,YD,ZD,1,XI(IXI),YI(IYI),ZI(IXI,IYI),IER,WK)
              if (IER > 0) STOP 1
              DZI(IXI,IYI) = ZI(IXI,IYI) - ZIE(IXI,IYI)
          end do
      end do
! - Writes the calculated results.
      write (6,FMT=9030) NMPR,NMSR(1),LBL(1)
      do ISEC = 1,2
          if (ISEC == 1) then
              IXIMN = 1
              IXIMX = 11
          else
              IXIMN = 9
              IXIMX = NXI
          end if
          write (6,FMT=9040) (XI(IXI),IXI=IXIMN,IXIMX)
          do IYIR = 1,NYI
              IYI = NYI + 1 - IYIR
              write (6,FMT=9050) YI(IYI), (ZI(IXI,IYI),IXI=IXIMN,IXIMX)
           end do
      end do
! - Writes the differences.
      write (6,FMT=9030) NMPR,NMSR(1),LBL(2)
      do ISEC = 1,2
          if (ISEC == 1) then
              IXIMN = 1
              IXIMX = 11
          else
              IXIMN = 9
              IXIMX = NXI
          end if
          write (6,FMT=9060) (XI(IXI),IXI=IXIMN,IXIMX)
          do IYIR = 1,NYI
              IYI = NYI + 1 - IYIR
              write (6,FMT=9050) YI(IYI),(DZI(IXI,IYI),IXI=IXIMN,IXIMX)
          end do
      end do
! Program check for the RGSF3P_sg subroutine
! - Performs surface fitting and calculates the differences.
      MD = 1
      call interp_RGSF3P_sg(MD,NXD,NYD,XD,YD,ZD,NXI,XI,NYI,YI,ZI,IER,WK)
      if (IER > 0) STOP 1
      do IYI = 1,NYI
          do IXI = 1,NXI
              DZI(IXI,IYI) = ZI(IXI,IYI) - ZIE(IXI,IYI)
          end do
      end do
! - Writes the calculated results.
      write (6,FMT=9030) NMPR,NMSR(2),LBL(1)
      do ISEC = 1,2
          if (ISEC == 1) then
              IXIMN = 1
              IXIMX = 11
          else
              IXIMN = 9
              IXIMX = NXI
          end if
          write (6,FMT=9040) (XI(IXI),IXI=IXIMN,IXIMX)
          do IYIR = 1,NYI
              IYI = NYI + 1 - IYIR
              write (6,FMT=9050) YI(IYI), (ZI(IXI,IYI),IXI=IXIMN,IXIMX)
          end do
      end do
! - Writes the differences.
      write (6,FMT=9030) NMPR,NMSR(2),LBL(2)
      do ISEC = 1,2
          if (ISEC == 1) then
              IXIMN = 1
              IXIMX = 11
          else
              IXIMN = 9
              IXIMX = NXI
          end if
          write (6,FMT=9060) (XI(IXI),IXI=IXIMN,IXIMX)
          do IYIR = 1,NYI
              IYI = NYI + 1 - IYIR
              write (6,FMT=9050) YI(IYI),(DZI(IXI,IYI),IXI=IXIMN,IXIMX)
          end do
      end do
      return
! Format statements
 9000 FORMAT (A9,7X,'Original Data',/,/,/,/,35X,'ZD(XD,YD)')
 9010 FORMAT (4X,'YD    XD=',/,7X,F8.1,2 (1X,3F6.1,F7.1),/)
 9020 FORMAT (1X,F6.1,F8.1,2 (1X,3F6.1,F7.1))
 9030 FORMAT (/,A9,3X,'Program Check for ',A9,3X,A20)
 9040 FORMAT (1X,/,38X,'ZI(XI,YI)',/,2X,'YI',3X,'XI=',/,5X,3F15.10,2F15.10,2F15.10,2F15.10,2F15.10,/)
 9050 FORMAT (F5.2,3F15.10,2F15.10,2F15.10,2F15.10,2F15.10)
 9060 FORMAT (1X,/,38X,'DZI(XI,YI)',/,2X,'YI',3X,'XI=',/,5X,3F15.10,2F15.10,2F15.10,2F15.10,2F15.10,/)
      end subroutine test_akima_sg
