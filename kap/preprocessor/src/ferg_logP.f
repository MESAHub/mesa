! ***********************************************************************
!
!   Copyright (C) 2009-2019  Bill Paxton, Frank Timmes & The MESA Team
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

      module ferg_logP
      use interp_2D_lib_sg
      use math_lib
      use utils_lib, only: mesa_error

      implicit none
      
      ! range of (logR,logT) in logP data file
      real, parameter :: logP_logT_min = 2.70
      real, parameter :: logP_logT_max = 4.00
      real, parameter :: logP_dlogT = 0.05
      integer, parameter :: logP_num_logTs = 27
      
      real, parameter :: logP_logR_min =-6.00
      real, parameter :: logP_logR_max = 7.00
      real, parameter :: logP_dlogR = 1.00
      integer, parameter :: logP_num_logRs = 14
      
      real :: logP_logTs(logP_num_logTs)
      real :: logP_logRs(logP_num_logRs)
      real :: logPs(4, logP_num_logRs, logP_num_logTs)
            
      ! range of (logP,logT) in logK data file
      real, parameter :: logK_logT_min = 2.70
      real, parameter :: logK_logT_max = 4.00
      real, parameter :: logK_dlogT = 0.05
      integer, parameter :: logK_num_logTs = 27
      
      real, parameter :: logK_logP_min =-9.00
      real, parameter :: logK_logP_max = 10.00
      real, parameter :: logK_dlogP = 1.00
      integer, parameter :: logK_num_logPs = 20
      
      real :: logK_logTs(logK_num_logTs)
      real :: logK_logPs(logK_num_logPs)
      real :: logKs(4, logK_num_logPs, logK_num_logTs)
            
      contains

      
      
      subroutine Init_ferg_logP(Z, X, dir, logP_dir)
         double precision, intent(IN) :: Z, X
         character (len=256), intent(IN) :: dir, logP_dir
         character (len=256) :: fname
         ! read the table of logP values for the given Z and X
         call Read_logPs(Z, X, logP_dir)
         call Get_logK_filename(Z, X, dir, fname)
         call Read_logKs(fname)
      end subroutine Init_ferg_logP


      
      subroutine Get_logK_filename(Z, X, dir, fname)
         save
         double precision, intent(IN) :: Z, X
         character (len=*), intent(IN) :: dir
         character (len=*), intent(OUT) :: fname
         character (len=256) :: zstr, xstr
         integer :: iz, ix

         iz = int(z * 1d3 + 0.5d0)
         ix = int(x * 1d2 + 0.5d0)
         zstr = ''
         if (iz == 0) zstr = '0'
         if (iz == 1) zstr = '001'
         if (iz == 4) zstr = '004'
         if (iz == 10) zstr = '01'
         if (iz == 20) zstr = '02'
         if (iz == 30) zstr = '03'
         if (iz == 50) zstr = '05'
         if (iz == 100) zstr = '1'
         xstr = ''
         if (ix == 0) xstr = '0'
         if (ix == 3) xstr = '03'
         if (ix == 10) xstr = '1'
         if (ix == 35) xstr = '35'
         if (ix == 70) xstr = '7'
         
         write(fname,'(6A)') trim(dir), '/gs98.hd.', trim(xstr), '.', trim(zstr), '.tron'
         
      end subroutine Get_logK_filename
      
   
      subroutine Read_logKs(fname)
         character (len=*), intent(IN) :: fname
      
         integer :: io_logK, i, ios       
         real :: logT, lKs(logK_num_logPs)
      
         integer :: ibcxmin                   ! bc flag for x=xmin
         real :: bcxmin(logK_num_logTs)               ! bc data vs. y at x=xmin
         integer :: ibcxmax                   ! bc flag for x=xmax
         real :: bcxmax(logK_num_logTs)               ! bc data vs. y at x=xmax
         integer :: ibcymin                   ! bc flag for y=ymin
         real :: bcymin(logK_num_logPs)               ! bc data vs. x at y=ymin
         integer :: ibcymax                   ! bc flag for y=ymax
         real :: bcymax(logK_num_logPs)               ! bc data vs. x at y=ymax
         integer :: ili_logPs                 
         integer :: ili_logTs                 
         integer :: ier                       ! =0 on exit if there is no error.
      
         io_logK = 40
         open(unit=io_logK, file=trim(fname), action='read', status='old', iostat=ios)
         if (ios /= 0) then
            write(*,*) 'failed to open', trim(fname)
            call mesa_error(__FILE__,__LINE__)
         end if
   
         do i = 1, 3
            read(io_logK, *)
         end do

         read(io_logK, '(6x,99(f7.3))') logK_logPs(1:logK_num_logPs)
            
         i = logK_num_logTs
         logT = logK_logT_max
         do 
            read(io_logK,'(f5.3,1x,99f7.3)') logK_logTs(i), lKs(1:logK_num_logPs)
            logKs(1, 1:logK_num_logPs, i) = lKs(1:logK_num_logPs)
            if (abs(logT - logK_logTs(i)) < 0.01) then
               logT = logK_logT_max - (logK_num_logTs-i+1)*logK_dlogT
               i = i - 1
               if (i < 1) exit
            end if
         end do
            
         close(io_logK)
            
         ! just use "not a knot" bc's at edges of tables
         ibcxmin = 0; bcxmin(1:logK_num_logTs) = 0
         ibcxmax = 0; bcxmax(1:logK_num_logTs) = 0
         ibcymin = 0; bcymin(1:logK_num_logPs) = 0
         ibcymax = 0; bcymax(1:logK_num_logPs) = 0

         call interp_mkbicub_sg(logK_logPs, logK_num_logPs, logK_logTs, logK_num_logTs, logKs, logK_num_logPs, &
                  ibcxmin,bcxmin,ibcxmax,bcxmax, &
                  ibcymin,bcymin,ibcymax,bcymax, &
                  ili_logPs,ili_logTs,ier)
         if (ier /= 0) then
            write(*,*) 'interp_mkbicub_sg error happened for logKs'
            call mesa_error(__FILE__,__LINE__)
         end if

      end subroutine Read_logKs


      real function Find_logP(logR, logT)
         real, intent(IN) :: logR, logT
         real :: fval(6)
         integer :: ict(6), ier
         
         ict = 0; ict(1) = 1
         call evbicub(logR,logT,logP_logRs,logP_num_logRs, logP_logTs,logP_num_logTs, &
                             0, 0, logPs, logP_num_logRs, ict, fval, ier)
         Find_logP = fval(1)
         
      end function Find_logP
      

      real function Eval_lowTemp_logK(logP, logT)
         real, intent(IN) :: logP, logT
         real :: fval(6)
         integer :: ict(6), ier
         
         ict = 0; ict(1) = 1
         call evbicub(logP, logT, logK_logPs, logK_num_logPs, logK_logTs, logK_num_logTs, &
                             0, 0, logKs, logK_num_logPs, ict, fval, ier)
         Eval_lowTemp_logK = fval(1)
         
      end function Eval_lowTemp_logK
      

      subroutine Read_logPs(Z, X, dir)
         double precision, intent(IN) :: Z, X
         character (len=256), intent(IN) :: dir
      
         integer :: ix, iz, io_logP, i, j
         character (len=256) :: fname, xstr, zstr
         
         integer :: first_logR(logP_num_logTs), last_logR(logP_num_logTs)
         real :: logP, dlogP, lPs(logP_num_logRs)
      
         integer :: ibcxmin                   ! bc flag for x=xmin
         real :: bcxmin(logP_num_logTs)               ! bc data vs. y at x=xmin
         integer :: ibcxmax                   ! bc flag for x=xmax
         real :: bcxmax(logP_num_logTs)               ! bc data vs. y at x=xmax
         integer :: ibcymin                   ! bc flag for y=ymin
         real :: bcymin(logP_num_logRs)               ! bc data vs. x at y=ymin
         integer :: ibcymax                   ! bc flag for y=ymax
         real :: bcymax(logP_num_logRs)               ! bc data vs. x at y=ymax
         integer :: ili_logRs                     ! =1: logRho grid is "nearly" equally spaced
         integer :: ili_logTs                     ! =1: logT grid is "nearly" equally spaced
         integer :: ier                       ! =0 on exit if there is no error.

         iz = int(Z * 1d3 + 0.1d0)
         ix = int(X * 1d2 + 0.1d0)
         zstr = ''
         if (iz == 0) zstr = '000'
         if (iz == 1) zstr = '001'
         if (iz == 4) zstr = '004'
         if (iz == 10) zstr = '010'
         if (iz == 20) zstr = '020'
         if (iz == 30) zstr = '030'
         if (iz == 50) zstr = '050'
         if (iz == 100) zstr = '100'
         xstr = ''
         if (ix == 0) xstr = '00'
         if (ix == 3) xstr = '03'
         if (ix == 10) xstr = '10'
         if (ix == 35) xstr = '35'
         if (ix == 70) xstr = '70'
         
         write(fname,'(6A)') trim(dir), '/logP_z', trim(zstr), '_x', trim(xstr), '.data'
      
         io_logP = 40
         open(unit=io_logP,file=trim(fname))
   
         read(io_logP, *)
         read(io_logP, *)
         read(io_logP, *) logP_logRs(1:logP_num_logRs)
         read(io_logP, *)
            
         do i = 1, logP_num_logTs
            read(io_logP,*) logP_logTs(i), lPs(1:logP_num_logRs)
            logPs(1, 1:logP_num_logRs, i) = lPs(1:logP_num_logRs)
         end do
            
         close(io_logP)
         
         ! replace fillers by extrapolation
         do j = 1, logP_num_logTs
            first_logR(j) = 0
            do i = 1, logP_num_logRs
               if (logPs(1, i, j) > -9.0) then
                  if (first_logR(j) == 0) first_logR(j) = i
                  last_logR(j) = i
               end if
            end do
            
            if (first_logR(j) > 1) then ! extrapolate down
               i = first_logR(j)
               logP = logPs(1, i, j)
               dlogP = logPs(1, i+1, j) - logP
               do i = first_logR(j)-1, 1, -1
                  logP = logP - dlogP
                  logPs(1, i, j) = logP
               end do
            end if
               
            if (last_logR(j) < logP_num_logRs) then ! extrapolate up
               i = last_logR(j)
               logP = logPs(1, i, j)
               dlogP = logP - logPs(1, i-1, j)
               do i = last_logR(j)+1, logP_num_logRs
                  logP = logP + dlogP
                  logPs(1, i, j) = logP
               end do
            end if
         end do
            
         if (.false.) then
            write(*,'(7x,99(F7.3))') logP_logRs(1:logP_num_logRs)
            do i = 1, logP_num_logTs
               write(*,'(99(F7.3))') logP_logTs(i), logPs(1, 1:logP_num_logRs, i)
            end do
         end if
            
         ! just use "not a knot" bc's at edges of tables
         ibcxmin = 0; bcxmin(1:logP_num_logTs) = 0
         ibcxmax = 0; bcxmax(1:logP_num_logTs) = 0
         ibcymin = 0; bcymin(1:logP_num_logRs) = 0
         ibcymax = 0; bcymax(1:logP_num_logRs) = 0

         call interp_mkbicub_sg(logP_logRs, logP_num_logRs, logP_logTs, logP_num_logTs, logPs, logP_num_logRs, &
                  ibcxmin,bcxmin,ibcxmax,bcxmax, &
                  ibcymin,bcymin,ibcymax,bcymax, &
                  ili_logRs,ili_logTs,ier)
         if (ier /= 0) then
            write(*,*) 'interp_mkbicub_sg error happened for logPs'
            call mesa_error(__FILE__,__LINE__)
         end if
            
      end subroutine Read_logPs
      
      
      subroutine eval_ferg_logP (z_in,xh_in,xxc_in,xxo_in,t6_in,r_in,logKap)
         double precision z_in,xh_in,xxc_in,xxo_in,t6_in,r_in,logKap
         real z,xh,xxc,xxo,t6,r
         real logT, logR, logP
      
         z = real(z_in); xh = real(xh_in); xxc = real(xxc_in) 
         xxo = real(xxo_in); t6 = real(t6_in); r = real(r_in)

         xxc = xxc; xxo = xxo; xh = xh; z=z ! for now, we are ignoring these

         logT = log10(t6*1e6);
         if (logT > logP_logT_max) then
            logT = logK_logT_max
         end if
         if (logT < logP_logT_min) then
            logT = logK_logT_min
         end if
      
         logR = log10(R)
         if (logR > logP_logR_max) then
            logR = logP_logR_max
         end if
         if (logR < logP_logR_min) then
            logR = logP_logR_min
         end if
      
         logP = Find_logP(logR,logT)
         if (logP > logK_logP_max) logP = logK_logP_max
         
         logKap = Eval_lowTemp_logK(logP, logT)

      end subroutine eval_ferg_logP
      
      
      end module ferg_logP
