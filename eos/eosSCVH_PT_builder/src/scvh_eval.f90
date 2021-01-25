! ***********************************************************************
!
!   Copyright (C) 2009  Bill Paxton, Frank Timmes
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

      module scvh_eval
      use utils_lib, only: mesa_error
      implicit none
      
      logical, parameter :: dbg = .false.
      
      double precision, parameter :: max_scvh_logP = 19d0, min_scvh_logP = -0.6d0

      integer :: num_h_pts, num_he_pts
      double precision, pointer, dimension(:) :: X_h, X_he, Y_h, Y_he
      double precision, pointer, dimension(:) :: XD, YD
      
      type akima_info
         logical :: needs_initialization
         double precision, pointer, dimension(:,:) :: ZD ! (NXD,NYD)
         double precision, pointer, dimension(:,:,:) :: WK ! (3,NXD,NYD)
      end type akima_info
      
      type (akima_info), target :: &
            compv1_h_ai, compv2_h_ai, denlog_h_ai, slog_h_ai, ulog_h_ai,  &
            dddt_tab_h_ai, dddp_tab_h_ai, dsdt_tab_h_ai, dsdp_tab_h_ai, dtdp_tab_h_ai, &
            compv1_he_ai, compv2_he_ai, denlog_he_ai, slog_he_ai, ulog_he_ai,  &
            dddt_tab_he_ai, dddp_tab_he_ai, dsdt_tab_he_ai, dsdp_tab_he_ai, dtdp_tab_he_ai
            
            
      ! scvh data
      
      integer, parameter :: nt=63, np=99 
         ! nt is number of different temperatures; np is number of different pressures
      integer :: npoint(nt)
      double precision :: tlog(nt),plog(np)
      double precision, target :: compv1_h1a(4*nt*np),compv2_h1a(4*nt*np), &
                       denlog_h1a(4*nt*np),slog_h1a(4*nt*np), &
                       ulog_h1a(4*nt*np),dddt_tab_h1a(4*nt*np), &
                       dddp_tab_h1a(4*nt*np),dsdt_tab_h1a(4*nt*np), &
                       dsdp_tab_h1a(4*nt*np),dtdp_tab_h1a(4*nt*np)
      double precision, target :: compv1_he1a(4*nt*np),compv2_he1a(4*nt*np), &
                       denlog_he1a(4*nt*np),slog_he1a(4*nt*np), &
                       ulog_he1a(4*nt*np),dddt_tab_he1a(4*nt*np), &
                       dddp_tab_he1a(4*nt*np),dsdt_tab_he1a(4*nt*np), &
                       dsdp_tab_he1a(4*nt*np),dtdp_tab_he1a(4*nt*np)
      double precision, pointer, dimension(:) :: compv1_h1,compv2_h1, &
                       denlog_h1,slog_h1, &
                       ulog_h1,dddt_tab_h1, &
                       dddp_tab_h1,dsdt_tab_h1, &
                       dsdp_tab_h1,dtdp_tab_h1
      double precision, pointer, dimension(:) :: compv1_he1,compv2_he1, &
                       denlog_he1,slog_he1, &
                       ulog_he1,dddt_tab_he1, &
                       dddp_tab_he1,dsdt_tab_he1, &
                       dsdp_tab_he1,dtdp_tab_he1
      double precision, pointer, dimension(:,:,:) :: compv1_h,compv2_h, &
                       denlog_h,slog_h, &
                       ulog_h,dddt_tab_h, &
                       dddp_tab_h,dsdt_tab_h, &
                       dsdp_tab_h,dtdp_tab_h
      double precision, pointer, dimension(:,:,:) :: compv1_he,compv2_he, &
                       denlog_he,slog_he, &
                       ulog_he,dddt_tab_he, &
                       dddp_tab_he,dsdt_tab_he, &
                       dsdp_tab_he,dtdp_tab_he
      
      
      
      
      
      
      contains
      
      
      
      subroutine setup_scvh(data_dir)
         character (len=*), intent(IN) :: data_dir
         call setup_AKIMA(data_dir)
         
         compv1_h1 => compv1_h1a
         compv2_h1 => compv2_h1a
         denlog_h1 => denlog_h1a
         slog_h1 => slog_h1a
         ulog_h1 => ulog_h1a
         dddt_tab_h1 => dddt_tab_h1a
         dddp_tab_h1 => dddp_tab_h1a
         dsdt_tab_h1 => dsdt_tab_h1a
         dsdp_tab_h1 => dsdp_tab_h1a
         dtdp_tab_h1 => dtdp_tab_h1a

         compv1_he1 => compv1_he1a
         compv2_he1 => compv2_he1a
         denlog_he1 => denlog_he1a
         slog_he1 => slog_he1a
         ulog_he1 => ulog_he1a
         dddt_tab_he1 => dddt_tab_he1a
         dddp_tab_he1 => dddp_tab_he1a
         dsdt_tab_he1 => dsdt_tab_he1a
         dsdp_tab_he1 => dsdp_tab_he1a
         dtdp_tab_he1 => dtdp_tab_he1a

         compv1_h(1:4,1:nt,1:np) => compv1_h1(1:4*nt*np)
         compv2_h(1:4,1:nt,1:np) => compv2_h1(1:4*nt*np)
         denlog_h(1:4,1:nt,1:np) => denlog_h1(1:4*nt*np)
         slog_h(1:4,1:nt,1:np) => slog_h1(1:4*nt*np)
         ulog_h(1:4,1:nt,1:np) => ulog_h1(1:4*nt*np)
         dddt_tab_h(1:4,1:nt,1:np) => dddt_tab_h1(1:4*nt*np)
         dddp_tab_h(1:4,1:nt,1:np) => dddp_tab_h1(1:4*nt*np)
         dsdt_tab_h(1:4,1:nt,1:np) => dsdt_tab_h1(1:4*nt*np)
         dsdp_tab_h(1:4,1:nt,1:np) => dsdp_tab_h1(1:4*nt*np)
         dtdp_tab_h(1:4,1:nt,1:np) => dtdp_tab_h1(1:4*nt*np)

         compv1_he(1:4,1:nt,1:np) => compv1_he1(1:4*nt*np)
         compv2_he(1:4,1:nt,1:np) => compv2_he1(1:4*nt*np)
         denlog_he(1:4,1:nt,1:np) => denlog_he1(1:4*nt*np)
         slog_he(1:4,1:nt,1:np) => slog_he1(1:4*nt*np)
         ulog_he(1:4,1:nt,1:np) => ulog_he1(1:4*nt*np)
         dddt_tab_he(1:4,1:nt,1:np) => dddt_tab_he1(1:4*nt*np)
         dddp_tab_he(1:4,1:nt,1:np) => dddp_tab_he1(1:4*nt*np)
         dsdt_tab_he(1:4,1:nt,1:np) => dsdt_tab_he1(1:4*nt*np)
         dsdp_tab_he(1:4,1:nt,1:np) => dsdp_tab_he1(1:4*nt*np)
         dtdp_tab_he(1:4,1:nt,1:np) => dtdp_tab_he1(1:4*nt*np)
         
      end subroutine setup_scvh
      
      
      
      subroutine interp_vals( &
            den_h,ener_h,entr_h,dddt_cp_h,dddp_ct_h,dsdt_cp_h,dsdp_ct_h,dtdp_cs_h, &
            den_he,ener_he,entr_he,dddt_cp_he,dddp_ct_he,dsdt_cp_he,dsdp_ct_he,dtdp_cs_he, &
            xnh,xnh2,xnhe,xnhep,logT,logP,P,info)
         use interp_2D_lib_db
         double precision, intent(out) ::  &
            den_h,ener_h,entr_h,dddt_cp_h,dddp_ct_h,dsdt_cp_h,dsdp_ct_h,dtdp_cs_h, &
            den_he,ener_he,entr_he,dddt_cp_he,dddp_ct_he,dsdt_cp_he,dsdp_ct_he,dtdp_cs_he, &
            xnh,xnh2,xnhe,xnhep
         double precision, intent(in) :: logT,logP,P
         integer, intent(out) :: info
         info = 0
         call AKIMA_interp_vals( &
                     den_h,ener_h,entr_h,dddt_cp_h,dddp_ct_h,dsdt_cp_h,dsdp_ct_h,dtdp_cs_h, &
                     den_he,ener_he,entr_he,dddt_cp_he,dddp_ct_he,dsdt_cp_he,dsdp_ct_he,dtdp_cs_he, &
                     xnh,xnh2,xnhe,xnhep,logT,logP,P,info)
      end subroutine interp_vals
      
      
      subroutine setup_AKIMA(data_dir)
         character (len=*), intent(IN) :: data_dir
         call alloc_for_AKIMA
         call read_data_for_AKIMA(data_dir)
         call setup_for_interp_AKIMA
      end subroutine setup_AKIMA
      
      
      subroutine alloc_for_AKIMA
         
!..bring in the common block
         integer :: info
         integer, parameter :: NXD=np, NYD=nt
         info = 0
         
         allocate(XD(NXD),YD(NYD),stat=info)
         if (info /= 0) call do_stop('allocate failed in alloc_for_AKIMA')
         
         call alloc(compv1_h_ai,compv1_he_ai)
         call alloc(compv2_h_ai,compv2_he_ai)
         call alloc(denlog_h_ai,denlog_he_ai)
         call alloc(slog_h_ai,slog_he_ai)
         call alloc(ulog_h_ai,ulog_he_ai)
         call alloc(dddt_tab_h_ai,dddt_tab_he_ai)
         call alloc(dddp_tab_h_ai,dddp_tab_he_ai)
         call alloc(dsdt_tab_h_ai,dsdt_tab_he_ai)
         call alloc(dsdp_tab_h_ai,dsdp_tab_he_ai)
         call alloc(dtdp_tab_h_ai,dtdp_tab_he_ai)
     
         contains
         
         subroutine alloc(info_h,info_he)
            ! set up the interpolation information
            type (akima_info) :: info_h, info_he
            integer :: ierr
            ierr = 0
            allocate(info_h% ZD(NXD,NYD), info_h% WK(3,NXD,NYD), stat=ierr)
            if (ierr /= 0) call do_stop('allocate failed in alloc_for_AKIMA')
            info_h% ZD = 0
            info_h% needs_initialization = .true.
            allocate(info_he% ZD(NXD,NYD), info_he% WK(3,NXD,NYD), stat=ierr)
            if (ierr /= 0) call do_stop('allocate failed in alloc_for_AKIMA')
            info_he% ZD = 0
            info_he% needs_initialization = .true.
         end subroutine alloc

      end subroutine alloc_for_AKIMA
      
      
      subroutine read_data_for_AKIMA(data_dir)
         character (len=*), intent(IN) :: data_dir
         
!..bring in the common block
         integer :: j

         if (np /= 99) call do_stop('read_data_for_AKIMA expects np == 99')
         if (nt /= 63) call do_stop('read_data_for_AKIMA expects n6 == 63')
         
         ! NOTE: XD and YD must be in monotonic increasing order

         do j=1,np
            XD(j) = -0.60d0 + (j-1)*0.20d0 ! store logP values in XD
         end do
         if (abs(XD(np)-19d0) > 1d-4) call do_stop('read_data_for_AKIMA expects max logP of 19.0')
         do j=1,nt
            YD(j) = 2.10d0 + (j-1)*0.08d0 ! store logT values in YD
         end do
         if (abs(YD(nt)-7.06d0) > 1d-4) call do_stop('read_data_for_AKIMA expects max logT of 7.06')
         
         call read_file_for_AKIMA(data_dir,'scvh/h_tab.dat', num_h_pts, np, nt, X_h, Y_h, &
            compv1_h_ai, compv2_h_ai, denlog_h_ai, slog_h_ai, ulog_h_ai,  &
            dddt_tab_h_ai, dddp_tab_h_ai, dsdt_tab_h_ai, dsdp_tab_h_ai, dtdp_tab_h_ai, .false.)
     
         call read_file_for_AKIMA(data_dir,'scvh/he_tab.dat', num_he_pts, np, nt, X_he, Y_he, &
            compv1_he_ai, compv2_he_ai, denlog_he_ai, slog_he_ai, ulog_he_ai,  &
            dddt_tab_he_ai, dddp_tab_he_ai, dsdt_tab_he_ai, dsdp_tab_he_ai, dtdp_tab_he_ai, .false.)
     
      end subroutine read_data_for_AKIMA
      
      
      subroutine read_file_for_AKIMA(data_dir,filename,n,np,nt,X,Y, &
            compv1, compv2, denlog, slog, ulog,  &
            dddt_tab, dddp_tab, dsdt_tab, dsdp_tab, dtdp_tab, show_bounds)
         character (len=*) :: data_dir,filename
         integer, intent(in) :: n,np,nt
         double precision, pointer, dimension(:) :: X, Y ! X is logP, Y is logT
         type (akima_info) ::compv1, compv2, denlog, slog, ulog,  &
            dddt_tab, dddp_tab, dsdt_tab, dsdp_tab, dtdp_tab
         integer :: i, j, jj, nps, k
         double precision :: tlg,plg,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10
         character (len=256) :: full_name
         logical :: show_bounds
         include 'formats'
         write(full_name,'(3a)') trim(data_dir), '/', trim(filename)
         write(*,*) 'read ' // trim(full_name)
         open(unit=1,file=trim(full_name),action='read')
         do i=1,nt
            read(1,*) tlg, nps
            do jj=1,nps
               read(1,*) plg,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10
               k = 0
               do j=1,np
                  if (abs(XD(j)-plg) < 1d-6) then
                     k = j; exit
                  end if
               end do
               if (k == 0) then
                  write(*,*) 'failed to find pressure in XD', plg
                  write(*,*) 'tlg', i, tlg
                  call do_stop('bad pressure in table')
               end if
               compv1% ZD(k,i) = c1
               compv2% ZD(k,i) = c2
               denlog% ZD(k,i) = c3
               slog% ZD(k,i) = c4
               ulog% ZD(k,i) = c5
               dddt_tab% ZD(k,i) = c6
               dddp_tab% ZD(k,i) = c7
               dsdt_tab% ZD(k,i) = c8
               dsdp_tab% ZD(k,i) = c9
               dtdp_tab% ZD(k,i) = c10
               if (show_bounds) then
                  if (jj == 1) write(*,'(99f15.8)', advance = 'no') tlg, plg
                  if (jj == nps) write(*,'(99f15.8)') plg
               end if
            enddo
         enddo
         if (show_bounds) then
            write(*,*) trim(filename)
            stop 'read_file_for_AKIMA'
         end if
         close(unit=1)
         write(*,*) 'close ' // trim(filename)
      
      end subroutine read_file_for_AKIMA
      
      
      subroutine setup_for_interp_AKIMA
         ! no setup required for AKIMA
      end subroutine setup_for_interp_AKIMA
      
      
      subroutine AKIMA_interp_vals( &
            den_h,ener_h,entr_h,dddt_cp_h,dddp_ct_h,dsdt_cp_h,dsdp_ct_h,dtdp_cs_h, &
            den_he,ener_he,entr_he,dddt_cp_he,dddp_ct_he,dsdt_cp_he,dsdp_ct_he,dtdp_cs_he, &
            xnh,xnh2,xnhe,xnhep,logT,logP,pres,info)
         use interp_2D_lib_db
         double precision, intent(out) ::  &
            den_h,ener_h,entr_h,dddt_cp_h,dddp_ct_h,dsdt_cp_h,dsdp_ct_h,dtdp_cs_h, &
            den_he,ener_he,entr_he,dddt_cp_he,dddp_ct_he,dsdt_cp_he,dsdp_ct_he,dtdp_cs_he, &
            xnh,xnh2,xnhe,xnhep
         double precision, intent(in) :: logT,logP,pres
         
!..bring in the common block

         integer, parameter :: NIP = 1
         double precision :: XI(NIP), YI(NIP), ZI(NIP)
         integer :: info, NXD, NYD
         
         include 'formats'
         
         info = 0
         NXD = np
         NYD = nt
         
         XI(1) = logP
         YI(1) = logT

         den_h = 10.0**interp_value(num_h_pts, X_h, Y_h, denlog_h_ai); if (info /= 0) return
         den_he = 10.0**interp_value(num_he_pts, X_he, Y_he, denlog_he_ai); if (info /= 0) return
         
         !..hydrogen section

         !..get d (log den) / d (log p) for hydrogen
         dddp_ct_h  = den_h/pres * interp_value(num_h_pts, X_h, Y_h, dddp_tab_h_ai); if (info /= 0) return

         !..number concentration of h2 molecules
         xnh2 = interp_value(num_h_pts, X_h, Y_h, compv1_h_ai); if (info /= 0) return

         !..number concentration of neutral h atoms 
         xnh = interp_value(num_h_pts, X_h, Y_h, compv2_h_ai); if (info /= 0) return

         !.. entropy in erg/g/k
         entr_h = 10.0**interp_value(num_h_pts, X_h, Y_h, slog_h_ai); if (info /= 0) return

         !..internal energy in erg/g
         ener_h = 10.0**interp_value(num_h_pts, X_h, Y_h, ulog_h_ai); if (info /= 0) return

         !..d log rho/d log t|p
         !..logarithmic derivative of the density with respect 
         !..to the temperature at constant p.
         dddt_cp_h = interp_value(num_h_pts, X_h, Y_h, dddt_tab_h_ai); if (info /= 0) return

         !..d log rho/d log p|t
         !..logarithmic derivative of the density with respect 
         !..to the pressure at constant t.
         dddp_ct_h = interp_value(num_h_pts, X_h, Y_h, dddp_tab_h_ai); if (info /= 0) return

         !..d log s/d log t|p
         !..logarithmic derivative of the entropy with respect 
         !..to the temperature at constant p.
         dsdt_cp_h = interp_value(num_h_pts, X_h, Y_h, dsdt_tab_h_ai); if (info /= 0) return

         !..d log s/d log p|t
         !..logarithmic derivative of the entropy with respect 
         !..to the pressure at constant t.
         dsdp_ct_h = interp_value(num_h_pts, X_h, Y_h, dsdp_tab_h_ai); if (info /= 0) return
         
         !..d log t/d log p|s  adiabatic gradiant
         dtdp_cs_h = interp_value(num_h_pts, X_h, Y_h, dtdp_tab_h_ai); if (info /= 0) return
         

         !..helium section

         !..get d (log den) / d (log p) for helium
         dddp_ct_he  = den_he/pres * interp_value(num_he_pts, X_he, Y_he, dddp_tab_he_ai); if (info /= 0) return

         !..number concentration of neutral helium atoms
         xnhe = interp_value(num_he_pts, X_he, Y_he, compv1_he_ai); if (info /= 0) return

         !..number concentration of he+ ions
         xnhep = interp_value(num_he_pts, X_he, Y_he, compv2_he_ai); if (info /= 0) return

         !.. entropy in erg/g/k
         entr_he = 10.0**interp_value(num_he_pts, X_he, Y_he, slog_he_ai); if (info /= 0) return

         !..internal energy in erg/g
         ener_he = 10.0**interp_value(num_he_pts, X_he, Y_he, ulog_he_ai); if (info /= 0) return

         !..d log rho/d log t|p
         !..logarithmic derivative of the density with respect 
         !..to the temperature at constant p.
         dddt_cp_he = interp_value(num_he_pts, X_he, Y_he, dddt_tab_he_ai); if (info /= 0) return

         !..d log rho/d log p|t
         !..logarithmic derivative of the density with respect 
         !..to the pressure at constant t.
         dddp_ct_he = interp_value(num_he_pts, X_he, Y_he, dddp_tab_he_ai); if (info /= 0) return

         !..d log s/d log t|p
         !..logarithmic derivative of the entropy with respect 
         !..to the temperature at constant p.
         dsdt_cp_he = interp_value(num_he_pts, X_he, Y_he, dsdt_tab_he_ai); if (info /= 0) return

         !..d log s/d log p|t
         !..logarithmic derivative of the entropy with respect 
         !..to the pressure at constant t.
         dsdp_ct_he = interp_value(num_he_pts, X_he, Y_he, dsdp_tab_he_ai); if (info /= 0) return

         !..d log t/d log p|s  adiabatic gradiant
         dtdp_cs_he = interp_value(num_he_pts, X_he, Y_he, dtdp_tab_he_ai); if (info /= 0) return
      
      
         contains

         double precision function interp_value(n,X,Y,ai)
            integer, intent(in) :: n
            double precision, pointer, dimension(:) :: X, Y
            type (akima_info) :: ai
            integer :: MD
            if (ai% needs_initialization) then
               MD = 1
            else
               MD = 2
            end if
            call interp_RGBI3P_db(MD,NXD,NYD,XD,YD,ai% ZD,NIP,XI,YI,ZI,info,ai% WK)
            interp_value = ZI(1)
            ai% needs_initialization = .false.
         end function interp_value
      
      end subroutine AKIMA_interp_vals



      integer function num_scvh_P_T_pairs(fname)
         implicit none
         character (len=*), intent(IN) :: fname
         double precision :: t
         integer :: iounit, cnt, n, i, j
!..bring in the common block
         iounit = 40
         cnt = 0
         open(unit=iounit,file=trim(fname),action='read')
         do i=1,nt
            read(iounit,*) t, n
            cnt = cnt + n
            do j=1,n
               read(iounit,*)
            end do
         end do
         close(iounit)
         num_scvh_P_T_pairs = cnt
      end function num_scvh_P_T_pairs









      subroutine read_scvh(data_dir) ! read and create interpolation info
      use interp_2D_lib_db
      implicit none
      character (len=*), intent(IN) :: data_dir
      save

!..this routine read the scvh data files

!..set the arrary sizes 
!..for the hydrogen table
      integer          i,j,jj,k
      double precision plo,pdelta,pdum,dp_fac,plg,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10


!  input bdy condition data:
         integer :: ibcxmin                   ! bc flag for x=xmin
         double precision :: bcxmin(nt)               ! bc data vs. y at x=xmin
         integer :: ibcxmax                   ! bc flag for x=xmax
         double precision :: bcxmax(nt)               ! bc data vs. y at x=xmax

         integer :: ibcymin                   ! bc flag for y=ymin
         double precision :: bcymin(np)               ! bc data vs. x at y=ymin
         integer :: ibcymax                   ! bc flag for y=ymax
         double precision :: bcymax(np)               ! bc data vs. x at y=ymax
      
      integer ilinx, iliny, info, ier
      
      integer ict(6)                    ! code specifying output desired
      double precision fval(6)                      ! output data
      character (len=256) :: filename
      

      double precision :: output_values(np,nt)
      integer :: i_lo, i_hi, j_lo, j_hi, j_mx



!..open the hydrogen file and read it
      
      if (nt==63) then
         write(filename,'(2a)') trim(data_dir), '/scvh/h_tab.dat'
         do j=1,np
            plog(j) = -0.60d0 + (j-1)*0.20d0
         end do
      else if (nt==621) then
         write(*,*) 'need to set plog'
         write(filename,'(2a)') trim(data_dir), '/scvh/H_TAB_I_HD.dat'
         !write(filename,'(2a)') trim(data_dir), '/SCVH_Akima/H_TAB_I_Akima.dat'
         !write(filename,'(2a)') trim(data_dir), '/SCVH_Renka/H_TAB_I_Renka.dat'
         call do_stop('read hires hydrogen file')
      else
         call do_stop('bad value for nt')
      end if
      write(*,*) 'read ' // trim(filename)
      open(unit=1,file=trim(filename),action='read')
      compv1_h = 0
      compv2_h = 0
      denlog_h = 0
      slog_h = 0
      ulog_h = 0
      dddt_tab_h = 0
      dddp_tab_h = 0
      dsdt_tab_h = 0
      dsdp_tab_h = 0
      dtdp_tab_h = 0
      do i=1,nt
         read(1,*) tlog(i), npoint(i)
         do jj=1,npoint(i)
            read(1,*) plg,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10
            k = 0
            do j=1,np
               if (abs(plog(j)-plg) < 1d-6) then
                  k = j; exit
               end if
            end do
            if (k == 0) then
               write(*,*) 'failed to find pressure in plog table', plg
               write(*,*) 'tlog(i)', i, tlog(i)
               call do_stop('bad pressure for h table')
            end if
            compv1_h(1,i,j) = c1
            compv2_h(1,i,j) = c2
            denlog_h(1,i,j) = c3
            slog_h(1,i,j) = c4
            ulog_h(1,i,j) = c5
            dddt_tab_h(1,i,j) = c6
            dddp_tab_h(1,i,j) = c7
            dsdt_tab_h(1,i,j) = c8
            dsdp_tab_h(1,i,j) = c9
            dtdp_tab_h(1,i,j) = c10
         enddo
      enddo
      close(unit=1)
      write(*,*) 'close ' // trim(filename)

!..open the helium file and read it

      if (nt==63) then
         write(filename,'(2a)') trim(data_dir), '/scvh/he_tab.dat'
      else if (nt==621) then
         write(filename,'(2a)') trim(data_dir), '/scvh/He_TAB_I_HD.dat'
         !write(filename,'(2a)') trim(data_dir), '/SCVH_Akima/He_TAB_I_Akima.dat'
         !write(filename,'(2a)') trim(data_dir), '/SCVH_Renka/He_TAB_I_Renka.dat'
      else
         call do_stop('bad value for nt')
      end if
      write(*,*) 'read ' // trim(filename)
      open(unit=1,file=trim(filename),action='read')
      compv1_he = 0
      compv2_he = 0
      denlog_he = 0
      slog_he = 0
      ulog_he = 0
      dddt_tab_he = 0
      dddp_tab_he = 0
      dsdt_tab_he = 0
      dsdp_tab_he = 0
      dtdp_tab_he = 0
      do i=1,nt
         read(1,*) tlog(i),npoint(i)
         do jj=1,npoint(i)
            read(1,*) plg,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10
            k = 0
            do j=1,np
               if (abs(plog(j)-plg) < 1d-6) then
                  k = j; exit
               end if
            end do
            if (k == 0) then
               write(*,*) 'failed to find pressure in plog table', plg
               write(*,*) 'tlog(i)', i, tlog(i)
               call do_stop('bad pressure for he table')
            end if
            compv1_he(1,i,j) = c1
            compv2_he(1,i,j) = c2
            denlog_he(1,i,j) = c3
            slog_he(1,i,j) = c4
            ulog_he(1,i,j) = c5
            dddt_tab_he(1,i,j) = c6
            dddp_tab_he(1,i,j) = c7
            dsdt_tab_he(1,i,j) = c8
            dsdp_tab_he(1,i,j) = c9
            dtdp_tab_he(1,i,j) = c10
         enddo
      enddo
      close(unit=1)
      write(*,*) 'close ' // trim(filename)


!..use "not a knot" bc's
      ibcxmin = 0; bcxmin(1:nt) = 0
      ibcxmax = 0; bcxmax(1:nt) = 0
      ibcymin = 0; bcymin(1:np) = 0
      ibcymax = 0; bcymax(1:np) = 0
      
      call mk_interp(compv1_h1); if (info /= 0) return         
      call mk_interp(compv2_h1); if (info /= 0) return      
      call mk_interp(denlog_h1); if (info /= 0) return
      call mk_interp(slog_h1); if (info /= 0) return        
      call mk_interp(ulog_h1); if (info /= 0) return        
      call mk_interp(dddt_tab_h1); if (info /= 0) return         
      call mk_interp(dddp_tab_h1); if (info /= 0) return         
      call mk_interp(dsdt_tab_h1); if (info /= 0) return         
      call mk_interp(dsdp_tab_h1); if (info /= 0) return         
      call mk_interp(dtdp_tab_h1); if (info /= 0) return

      call mk_interp(compv1_he1); if (info /= 0) return         
      call mk_interp(compv2_he1); if (info /= 0) return         
      call mk_interp(denlog_he1); if (info /= 0) return
      call mk_interp(slog_he1); if (info /= 0) return         
      call mk_interp(ulog_he1); if (info /= 0) return         
      call mk_interp(dddt_tab_he1); if (info /= 0) return         
      call mk_interp(dddp_tab_he1); if (info /= 0) return         
      call mk_interp(dsdt_tab_he1); if (info /= 0) return         
      call mk_interp(dsdp_tab_he1); if (info /= 0) return         
      call mk_interp(dtdp_tab_he1); if (info /= 0) return
      
      
      contains
      
      subroutine mk_interp(table1)
         double precision, pointer :: table1(:)
         call interp_mkbicub_db( &
            tlog,nt,plog,np,table1,nt, &
            ibcxmin,bcxmin,ibcxmax,bcxmax, &
            ibcymin,bcymin,ibcymax,bcymax, &
            ilinx,iliny,info)
      end subroutine mk_interp
      
      end subroutine read_scvh


      subroutine do_stop(str)
         character (len=*) :: str
         write(*,*) trim(str)
         call mesa_error(__FILE__,__LINE__)
      end subroutine do_stop
      







      subroutine interpolate_scvh( &
         include_radiation,logT,T,logPgas,Pgas,xmassh1, &
         logRho, logE, logS, chiRho, chiT, &
         Cp, Cv, dE_dRho, dS_dT, dS_dRho, &
         mu, log_free_e, gamma1, gamma3, grad_ad, eta, &
         info)
      
      use interp_2D_lib_db
      use eos_lib, only: eos_fermi_dirac_integral  
      implicit none
      save

!..this routine interpolates the scvh tables

!..declare the pass  
      integer, intent(in) :: include_radiation
      double precision, intent(in) :: logT, T, logPgas, Pgas, xmassh1
      double precision, intent(out) :: &
         logRho, logE, logS, chiRho, chiT, &
         Cp, Cv, dE_dRho, dS_dT, dS_dRho, &
         mu, log_free_e, gamma1, gamma3, grad_ad, eta
      integer info ! returned = 0 if AOK
      
      double precision xnh,xnh2,xnhe,xnhep,dpdd,dpdt,ener,dedd, &
         Rho,din,xnhp,xnhepp,logtin,P,tin,pres,gam3,entr,dsdd,dsdt,xtra, &
         dxdd,dxdt,Ne,logNe,log_free_e0,log_free_e1,kt,theta,free_e
     
      
      double precision xmasshe4, dedt

!..local variables
      
      double precision prad,erad,srad

      double precision logpin,inv_din,inv_tin,inv_pres,small_value,eta_ele, &
         f, new_eta, alfa, fd, fdeta, fdtheta,eta_min,eta_max,coef
      parameter        (small_value = 1.0d-16)
      double precision, parameter :: tiny_eta_ele = -20d0

!..for hydrogen
      double precision den_h,ener_h,entr_h,dddt_cp_h,dddp_ct_h,dsdt_cp_h, &
         dsdp_ct_h,dtdp_cs_h,dpdd_h,dpdt_h,dedd_h,dedt_h,dsdd_h,dsdt_h


!..for helium
      double precision den_he,ener_he,entr_he,dddt_cp_he,dddp_ct_he,dsdt_cp_he, &
                       dsdp_ct_he,dtdp_cs_he,dpdd_he,dpdt_he,dedd_he,dedt_he,dsdd_he,dsdt_he


!..for the mixture
      double precision beta,gama,delt,smix,dddt_cp_hhe,dddp_ct_hhe,dsdt_cp_hhe, &
                       dsdp_ct_hhe,dtdp_cs_hhe,dpdd_hhe,dpdt_hhe,dedd_hhe,dedt_hhe,dsdd_hhe,dsdt_hhe
      double precision tiny
      parameter        (tiny = 1.0d-30)


!..for the search
      integer          ii,itmax
      parameter        (itmax = 200)
      double precision pmax,pmin,den_calc,func,denom,eostol
      parameter        (eostol = 1.0d-6)


!..for the photons
      integer          radmult
      parameter        (radmult = 1)
      double precision dpraddd,dpraddt,deraddd,deraddt,dsraddd,dsraddt

!..constants
      double precision clight,ssol,asol,asoli3,avo,kerg,xka,mh1,mhe4,third
      parameter        (clight  = 2.99792458d10, ssol    = 5.67051d-5, &
                        asol    = 4.0d0 * ssol / clight, asoli3  = asol/3.0d0, &
                        avo     = 6.0221367d23,kerg    = 1.380658d-16, &
                        xka     = kerg*avo,mh1     = 1.67357d-24, &
                        mhe4    = 6.646442d-24,third   = 1.0d0/3.0d0)
      double precision, parameter ::  pi = 3.1415926535897932384d0
      double precision, parameter ::  logAvo = 23.7797506095929d0
      double precision, parameter ::  h = 6.6260755d-27
      double precision, parameter ::  me = 9.1093897d-28
      double precision, parameter ::  mecc = me * clight * clight


      logical, parameter :: pure_splines = .true.

!..for the interpolations
      integer          iat, jat, i, j
      double precision dt,dt2,dti,dt2i,dp2,dpi,dp2i,xt,xp,mxt,mxp, &
                       si0t,si1t,si2t,si0mt,si1mt,si2mt,si0p,si1p,si2p,si0mp,si1mp,si2mp, &
                       z,fi_h(36),fi_he(36),w0t,w1t,w2t,w0mt,w1mt,w2mt,w0p,w1p,w2p,w0mp,w1mp,w2mp

      include 'formats'
      
      info = 0
      
      z=0; w0t=0; w1t=0; w2t=0; w0mt=0; w1mt=0; w2mt=0
      w0p=0; w1p=0; w2p=0; w0mp=0; w1mp=0; w2mp=0; i=0; j=0

!..initialize
      
      xnh2    = small_value
      xnh     = small_value
      xnhe    = small_value
      xnhep   = small_value
      dpdd    = small_value
      dpdt    = small_value
      dedd    = small_value
      dedt    = small_value
      dsdd    = small_value
      dsdt    = small_value
      pres    = small_value
      ener    = small_value
      entr    = small_value

      xmasshe4 = 1 - xmassh1

      logtin = logT
      tin = T
      inv_tin = 1.0d0/T


!..interpolate

      logpin   = logPgas
      pres     = Pgas

      call interp_vals( &
         den_h,ener_h,entr_h,dddt_cp_h,dddp_ct_h,dsdt_cp_h,dsdp_ct_h,dtdp_cs_h, &
         den_he,ener_he,entr_he,dddt_cp_he,dddp_ct_he,dsdt_cp_he,dsdp_ct_he,dtdp_cs_he, &
         xnh,xnh2,xnhe,xnhep,logtin,logpin,pres,info)
      if (info /= 0) then
         write(*,*) 'scvh interp_vals error'
         write(*,1) 'logT', logT
         write(*,1) 'logPgas', logPgas
         write(*,1) 'X', xmassh1
         write(*,*)
         return
      end if
      
      if (dsdt_cp_h <= 0 .or. dsdt_cp_he <= 0) then ! indicates off table
         info = 1 ! indicates off table
         !write(*,1) 'off table for scvh: logPgas, logW, logT, X', logPgas, logPgas - 4*logT, logT, xmassh1         
         return
      end if
      
      ! density
      din = 1.0d0/(xmassh1/den_h + xmasshe4/den_he)
      inv_din = 1.0d0/din

!..radiation section:
      prad    = 0.0
      dpraddd = 0.0
      dpraddt = 0.0

      erad    = 0.0
      deraddd = 0.0
      deraddt = 0.0

      srad    = 0.0
      dsraddd = 0.0
      dsraddt = 0.0

      if (radmult .ne. 0) then
       prad    = asoli3 * tin * tin * tin * tin
       dpraddd = 0.0
       dpraddt = 4.0d0 * prad * inv_tin

       erad    = 3.0d0 * prad * inv_din
       deraddd = -erad * inv_din
       deraddt = 3.0d0 * dpraddt * inv_din

       srad    = (prad*inv_din + erad) * inv_tin
       dsraddd = ((dpraddd - prad*inv_din)*inv_din + deraddd)*inv_tin
       dsraddt = (dpraddt*inv_din + deraddt - srad)*inv_tin
      end if

!..form the mixture totals
!..first get an approximation to the mixture entropy, eq 53 of scvh95
!..avoid beta being zero or infinite
!..avoid delt being <= zero

 
      beta = (mh1 * (xmasshe4+tiny)) / (mhe4 * (xmassh1 + tiny)) 

      gama = 1.5d0 * (1.0d0 + xnh + 3.0d0*xnh2) / (1.0d0 + 2.0d0*xnhe + xnhep)

      delt = 1.5d0 * beta * gama * (2.0d0 - 2.0d0*xnhe - xnhep) / (1.0d0 - xnh2 - xnh)

      if (delt <= 1d-10) delt = 1d-10 ! ADDED BY BP

      smix = xmassh1/mh1 * (2.0d0/(1.0d0 + xnh + 3.0d0*xnh2)) *  &
              (log(1.0d0 + beta*gama) - 0.5d0*(1.0d0 - xnh2 - xnh)*log(1d0 + delt) +  &
               beta*gama*(log(1.0d0 + 1.0d0/(beta*gama)) -  &
               third*(2.0d0 - 2.0d0*xnhe - xnhep) * log(1.0d0 + 1.0d0/delt)))
      smix = kerg * smix

!..equations 40 and 41 of scvh95
!..[BP] we're missing smix in entr
      ener = xmassh1*ener_h + xmasshe4*ener_he  ! eqn 40
      entr = xmassh1*entr_h + xmasshe4*entr_he  ! eqn 41
      
!..form the mixture derivatives
!..equations 43 to 47 of scvh95
!..[BP] we're missing the partials of smix in the these expressions
! FXT says we should just skip smix since the partials would involve differencing the tables.
 
      dddt_cp_hhe = din*(xmassh1/den_h*dddt_cp_h + xmasshe4/den_he*dddt_cp_he)  ! eqn 43

      dddp_ct_hhe = din*(xmassh1/den_h*dddp_ct_h + xmasshe4/den_he*dddp_ct_he)  ! eqn 44

      dsdt_cp_hhe = din*(xmassh1/den_h*dsdt_cp_h + xmasshe4/den_he*dsdt_cp_he)  ! eqn 45
         ! need smix/entr*dlogSmix/dlogT|P

      dsdp_ct_hhe = din*(xmassh1/den_h*dsdp_ct_h + xmasshe4/den_he*dsdp_ct_he)  ! eqn 46
         ! need smix/entr*dlogSmix/dlogP|T

      dtdp_cs_hhe = -dsdp_ct_hhe/dsdt_cp_hhe  ! eqn 47

!..delog the derivatives
      inv_pres = 1.0d0/pres
      dddt_cp_h  = din * inv_tin     * dddt_cp_h
      dddp_ct_h  = din * inv_pres    * dddp_ct_h
      dsdt_cp_h  = entr_h * inv_tin  * dsdt_cp_h
      dsdp_ct_h  = entr_h * inv_pres * dsdp_ct_h
      dtdp_cs_h  = tin * inv_pres    * dtdp_cs_h

      dddt_cp_he = din * inv_tin      * dddt_cp_he
      dddp_ct_he = din * inv_pres     * dddp_ct_he
      dsdt_cp_he = entr_he * inv_tin  * dsdt_cp_he
      dsdp_ct_he = entr_he * inv_pres * dsdp_ct_he
      dtdp_cs_he = tin * inv_pres     * dtdp_cs_he

      dddt_cp_hhe = din * inv_tin   * dddt_cp_hhe
      dddp_ct_hhe = din * inv_pres  * dddp_ct_hhe
      dsdt_cp_hhe = entr * inv_tin  * dsdt_cp_hhe
      dsdp_ct_hhe = entr * inv_pres * dsdp_ct_hhe
      dtdp_cs_hhe = tin * inv_pres  * dtdp_cs_hhe

!..form the usual thermodynamic derivatives 
!..for hydrogen

!..d(pres)/ d(den)|t
      dpdd_h  = 1.0d0/dddp_ct_h

!..d(entr)/d(den)|t
      dsdd_h = dddt_cp_h/(din*din*dddp_ct_h)

!..d(pres)/d(temp)|d
      dpdt_h = -din*din*dsdd_h

!..d(ener)/d(rho)|t
      dedd_h = (pres - tin * dpdt_h)*inv_din*inv_din

!..d(entr)/d(temp)|d
      dsdt_h = dsdp_ct_h * dpdt_h + dsdt_cp_h

!..d(ener)/d(temp)|d
      dedt_h = tin * dsdt_h



!..for helium
!..d(pres)/ d(den)|t
      dpdd_he  = 1.0d0/dddp_ct_he

!..d(entr)/d(den)|t
      dsdd_he = dddt_cp_he/(din*din*dddp_ct_he)

!..d(pres)/d(temp)|d
      dpdt_he = -din*din*dsdd_he

!..d(ener)/d(rho)|t
      dedd_he = (pres - tin * dpdt_he)*inv_din*inv_din

!..d(entr)/d(temp)|d
      dsdt_he = dsdp_ct_he * dpdt_he + dsdt_cp_he

!..d(ener)/d(temp)|d
      dedt_he = tin * dsdt_he



!..for the mixture
!..d(pres)/ d(den)|t
      dpdd_hhe  = 1.0d0/dddp_ct_hhe

!..d(entr)/d(den)|t
      dsdd_hhe = dddt_cp_hhe/(din*din*dddp_ct_hhe)

!..d(pres)/d(temp)|d
      dpdt_hhe = -din*din*dsdd_hhe

!..d(ener)/d(rho)|t
      dedd_hhe = (pres - tin * dpdt_hhe)*inv_din*inv_din

!..d(entr)/d(temp)|d
      dsdt_hhe = dsdp_ct_hhe * dpdt_hhe + dsdt_cp_hhe

!..d(ener)/d(temp)|d
      dedt_hhe = tin * dsdt_hhe

!..sum the components

      if (include_radiation /= 0) then
         pres = pres     + prad
         dpdd = dpdd_hhe + dpraddd
         dpdt = dpdt_hhe + dpraddt
   
         ener = ener     + erad
         dedd = dedd_hhe + deraddd
         dedt = dedt_hhe + deraddt
         
         entr = entr     + srad
         dsdd = dsdd_hhe + dsraddd
         dsdt = dsdt_hhe + dsraddt
      else
         dpdd = dpdd_hhe
         dpdt = dpdt_hhe
   
         dedd = dedd_hhe
         dedt = dedt_hhe
   
         dsdd = dsdd_hhe
         dsdt = dsdt_hhe
      end if
      
      gam3 = 1 + dpdt / (din * dedt) ! C&G 9.93      

      ! store the output

      Rho = din
      logE = log10(ener)
      logS = log10(entr)
      
      P = pres
      
      logRho = log10(Rho)
      
      gamma3 = gam3
      grad_ad = dtdp_cs_hhe/(tin * inv_pres)
      gamma1 = (gamma3 - 1) / grad_ad ! C&G 9.88 & 9.89
      
      chiRho = dpdd * Rho / P
      chiT = dpdt * T / P
      
      Cv = chiT * P / (rho * T * (gamma3 - 1)) ! C&G 9.93
      Cp = Cv + P * chiT*chiT / (Rho * T * chiRho) ! C&G 9.86
      
      dE_dRho = dedd
      
      dS_dT = dsdt
      dS_dRho = dsdd


      xnhp = max(0d0, min(1d0, 1 - (xnh2 + xnh)))
      xnhepp = max(0d0, min(1d0, 1 - (xnhep + xnhe)))
      mu = xmassh1 * (xnhp / 2 + xnh + 2 * xnh2) + 4 * (1 - xmassh1) * (xnhepp / 3 + xnhep / 2 + xnhe)
      
      Ne = Rho * avo * (xmassh1 * xnhp + (1 - xmassh1) * (xnhep + 2 * xnhepp) / 4 )
      if (Ne < 1d-99) Ne = 1d-99
      logNe = log10(Ne)
      if (logNe - logRho < 17d0) logNe = logRho + 17d0
         ! put a lower limit on the electron abundance to avoid spurious 
         ! bicubic interpolations at very low T
      free_e = 10**logNe / (avo * Rho) ! convert to mean number of free electrons per nucleon      
      
      ! calculate eta_ele as function of Ne and T
      ! get guess to start the process
      
      log_free_e = logNe - log10(avo) - logRho
      log_free_e0 = -2d0
      log_free_e1 = -4d0
      ! when number of free electrons per nucleon gets too low, eta becomes meaningless
      ! so cut it off

      kt      = kerg * T
      theta = kt/mecc

      if (log_free_e < log_free_e1) then
      
         eta_ele = tiny_eta_ele
         call eos_fermi_dirac_integral(0.5d0, eta_ele, theta, fd, fdeta, fdtheta)
         
      else
      
         eta_min = -1200d0
         eta_max = 1200d0
         eta_ele = 0d0
   
         coef    = 4 * pi * (2 * me * kt) ** 1.5d0 / h**3
         ne      = 10**logNe
      
         do i = 1, 100
   
            ! get fermi dirac integral
            call eos_fermi_dirac_integral(0.5d0, eta_ele, theta, fd, fdeta, fdtheta)
      
            if (fdeta < 0) then
               write(*,*) fd, fdeta
               write(*,*) 'expected fdeta > 0'
               call mesa_error(__FILE__,__LINE__)
            end if
         
            f = coef * fd - ne
            if (f > 0) then ! fd too large, make eta smaller to reduce it
               eta_max = eta_ele
            else  ! fd too small, make eta larger to reduce it
               eta_min = eta_ele
            end if
            new_eta = 0.5d0 * (eta_min + eta_max)
            if (abs(eta_ele - new_eta) < eostol) exit
            if (i == 100) then
               write(*,*) 'logT', logT
               write(*,*) 'logRho', logRho
               write(*,*) 'failed to find eta_ele by root solve'
               call mesa_error(__FILE__,__LINE__)
            end if
      
            eta_ele = new_eta
            if (eta_ele > 1000d0) then
               eta_ele = 1000d0; exit
            end if
            if (eta_ele < -1000d0) then
               eta_ele = -1000d0; exit
            end if
   
         end do
      
         if (log_free_e < log_free_e0) then
            alfa = (log_free_e - log_free_e1) / (log_free_e0 - log_free_e1)
            alfa = 0.5d0 * (1 - cos(pi*alfa))
            beta = 1 - alfa
            eta_ele = alfa * eta_ele + beta * tiny_eta_ele
         end if

      end if
      
      eta = eta_ele
      
      if (.false.) then
      write(*,*) 'include_radiation', include_radiation
      write(*,1) 'xmassh1', xmassh1
      write(*,1) 'den_h', den_h
      write(*,1) 'dddt_cp_h', dddt_cp_h
      write(*,1) 'xmasshe4', xmasshe4
      write(*,1) 'den_he', den_he
      write(*,1) 'dddt_cp_he', dddt_cp_he
      write(*,1) 'gam3', gam3
      write(*,1) 'dtdp_cs_hhe', dtdp_cs_hhe
      write(*,1) 'tin', tin
      write(*,1) 'inv_pres', inv_pres
      write(*,1) 'grad_ad', grad_ad
      write(*,1) 'logRho', logRho
      write(*,1) 'logT', logT
      write(*,*)
      stop 'scvh_eval'
      end if
      
      end subroutine interpolate_scvh



      end module scvh_eval

