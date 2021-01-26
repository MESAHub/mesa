! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton,Frank Timmes
!
!   This file is part of MESA.
!
!   MESA is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,or
!   (at your option) any later version.
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not,write to the Free Software
!   Foundation,Inc.,59 Temple Place,Suite 330,Boston,MA 02111-1307 USA
!
! ***********************************************************************

      module scvh_core
      implicit none
      
      logical,parameter :: dbg = .false.

      integer :: num_h_pts,num_he_pts
      double precision,pointer,dimension(:) :: logPs,logTs
      integer,pointer,dimension(:) :: k_for_MIN_logP,k_for_MAX_logP
      
      integer,parameter :: NlogPs = 99,NlogTs = 63
      double precision,parameter :: dlogP = 0.20d0
      double precision,parameter :: logP_min = -0.60d0
      double precision,parameter :: logP_max = 19d0
      double precision,parameter :: dlogT = 0.08d0
      double precision,parameter :: logT_min = 2.10d0
      double precision,parameter :: logT_max = 7.06d0

      type scvh_info
         logical :: needs_initialization
         double precision,pointer,dimension(:) :: Z1, f1
         double precision,pointer,dimension(:,:) :: Z ! (NlogPs,NlogTs)
         ! for alternative evaluation using bicubic splines
         double precision,pointer,dimension(:,:,:) :: f ! (4,NlogPs,NlogTs)
         integer :: ilinx ! =1: x grid is "nearly" equally spaced
         integer :: iliny ! =1: y grid is "nearly" equally spaced
         character (len=64) :: name
      end type scvh_info
      
      type (scvh_info),target ::
     >      compv1_h_si,compv2_h_si,denlog_h_si,slog_h_si,ulog_h_si, 
     >      compv1_he_si,compv2_he_si,denlog_he_si,slog_he_si,ulog_he_si
     
      
      
      contains


      subroutine setup_scvh(data_dir)
         character (len=*),intent(IN) :: data_dir
         call alloc_for_scvh
         call read_data_for_scvh(data_dir)
      end subroutine setup_scvh


      subroutine alloc_for_scvh

         integer :: info
         info = 0

         allocate(logPs(NlogPs),logTs(NlogTs),k_for_MAX_logP(NlogTs),k_for_MIN_logP(NlogTs),stat=info)
         if (info /= 0) call do_stop('allocate failed in alloc_for_scvh')
         
         call alloc(compv1_h_si,compv1_he_si,'compv1_h_si','compv1_he_si')
         call alloc(compv2_h_si,compv2_he_si,'compv2_h_si','compv2_he_si')
         call alloc(denlog_h_si,denlog_he_si,'denlog_h_si','denlog_he_si')
         call alloc(slog_h_si,slog_he_si,'slog_h_si','slog_he_si')
         call alloc(ulog_h_si,ulog_he_si,'ulog_h_si','ulog_he_si')
     
         contains
         
         subroutine alloc(info_h,info_he,name_h,name_he)
            ! set up the interpolation information
            type (scvh_info) :: info_h,info_he
            character (len=*),intent(in) :: name_h,name_he
            integer :: ierr
            ierr = 0
            allocate(info_h% Z1(NlogPs*NlogTs),info_h% f1(4*NlogPs*NlogTs),stat=ierr)
            if (ierr /= 0) call do_stop('allocate failed in alloc_for_scvh')
            info_h% Z(1:NlogPs,1:NlogTs) => info_h% Z1(1:NlogPs*NlogTs)
            info_h% f(1:4,1:NlogPs,1:NlogTs) => info_h% f1(1:4*NlogPs*NlogTs)
            info_h% Z = 0
            info_h% needs_initialization = .true.
            info_h% name = name_h
            allocate(info_he% Z1(NlogPs*NlogTs),info_he% f1(4*NlogPs*NlogTs),stat=ierr)
            if (ierr /= 0) call do_stop('allocate failed in alloc_for_scvh')
            info_he% Z(1:NlogPs,1:NlogTs) => info_he% Z1(1:NlogPs*NlogTs)
            info_he% f(1:4,1:NlogPs,1:NlogTs) => info_he% f1(1:4*NlogPs*NlogTs)
            info_he% Z = 0
            info_he% needs_initialization = .true.
            info_he% name = name_he
         end subroutine alloc

      end subroutine alloc_for_scvh


      subroutine read_data_for_scvh(data_dir)
         character (len=*),intent(IN) :: data_dir

         integer :: j
         
         include 'formats'
         
         ! NOTE: logPs and logTs must be in monotonic increasing order

         do j=1,NlogPs
            logPs(j) = logP_min + (j-1)*dlogP ! store logP values in logPs
         end do
         if (abs(logPs(NlogPs)-logP_max) > 1d-4) call do_stop('read_data_for_scvh expects max logP of 19.0')
         do j=1,NlogTs
            logTs(j) = logT_min + (j-1)*dlogT ! store logT values in logTs
         end do
         if (abs(logTs(NlogTs)-logT_max) > 1d-4) call do_stop('read_data_for_scvh expects max logT of 7.06')
         
         call read_file_for_scvh(data_dir,'scvh/h_tab.asc.data',num_h_pts,NlogPs,NlogTs,
     >         compv1_h_si,compv2_h_si,denlog_h_si,slog_h_si,ulog_h_si)
         call read_file_for_scvh(data_dir,'scvh/he_tab.asc.data',num_he_pts,NlogPs,NlogTs,
     >         compv1_he_si,compv2_he_si,denlog_he_si,slog_he_si,ulog_he_si)
         
         
         if (.false.) then
            call write_file_for_scvh(data_dir,'scvh/h_tab.new.data',num_h_pts,NlogPs,NlogTs,
     >         compv1_h_si,compv2_h_si,denlog_h_si,slog_h_si,ulog_h_si)
            call write_file_for_scvh(data_dir,'scvh/he_tab.new.data',num_he_pts,NlogPs,NlogTs,
     >         compv1_he_si,compv2_he_si,denlog_he_si,slog_he_si,ulog_he_si)
         
         end if
         
         
         call fill_and_smooth(compv1_h_si)
         call fill_and_smooth(compv2_h_si)
         call fill_and_smooth(denlog_h_si,6d0,-15d0)
         call fill_and_smooth(slog_h_si)
         call fill_and_smooth(ulog_h_si)
         
         call fill_and_smooth(compv1_he_si)
         call fill_and_smooth(compv2_he_si)
         call fill_and_smooth(denlog_he_si,6d0,-15d0)
         call fill_and_smooth(slog_he_si)
         call fill_and_smooth(ulog_he_si)
         
         ! do extra smoothing of compv1_he for logT < 4.5 and logP > 11.0 in 
         ! this is where we get problems from the He+ to He++ pressure ionization
         !call extra_smoothing(compv1_he_si,1,logTs(1),4.5d0,10d0,logPs(NlogPs))
         
         !call extra_smoothing(compv1_h_si,1,logTs(1),3.5d0,logPs(1),8d0)
         !call extra_smoothing(compv2_h_si,100,logTs(1),3.5d0,2d0,8d0)
         
                  
         contains
         
         
         subroutine fill_and_smooth(si,val_lr,val_ul)
            type (scvh_info) :: si
            double precision,intent(in),optional :: val_lr,val_ul
            integer :: j,kmax,kmin,k,i,kk
            double precision :: logP0,logP1,logT0,logT1,d0,d1,vlr,vul
            if (present(val_lr)) then
               vlr = val_lr
            else
               vlr = si% Z(k_for_MAX_logP(20)-4,20)
            end if
            if (present(val_ul)) then
               vul = val_ul
            else
               vul = si% Z(k_for_MIN_logP(NlogTs/2),NlogTs/2)
            end if
            do j = 1,NlogTs
               kmax = k_for_MAX_logP(j)
               do k = kmax+1,NlogPs
                  logP0 = logPs(k); logP1 = logP0
                  logT0 = logTs(j); logT1 = logT0
                  d0 = sqrt((logP_max+1-logP0)**2 + (logT_min-1-logT0)**2)
                  do i = j+1,NlogTs
                     logT1 = logT1 + dlogT
                     logP1 = logP1 - 4*dlogT
                     kk = k_for_MAX_logP(i)
                     if (logP1 <= logPs(kk) .or. i == NlogTs) then
                        d1 = sqrt((logP1-logP0)**2 + (logT1-logT0)**2)
                        si% Z(k,j) = si% Z(min(k,kk),i)*d0/(d0+d1) + vlr*d1/(d0+d1)
                        exit
                     end if
                  end do
               end do


               kmin = k_for_MIN_logP(j)
               do k = 1,kmin-1
                  logP0 = logPs(k); logP1 = logP0
                  logT0 = logTs(j); logT1 = logT0
                  d0 = sqrt((logP_min-1-logP0)**2 + (logT_max+1-logT0)**2)
                  do i = j-1,1,-1
                     logT1 = logT1 - dlogT
                     logP1 = logP1 + 4*dlogT
                     kk = k_for_MIN_logP(i)
                     if (logP1 >= logPs(kk) .or. i == 1) then
                        d1 = sqrt((logP1-logP0)**2 + (logT1-logT0)**2)
                        si% Z(k,j) = si% Z(max(k,kk),i)*d0/(d0+d1) + vul*d1/(d0+d1)
                        exit
                     end if
                  end do
               end do
            end do
            
            
            ! smooth the added entries
            do i = 1,20
               do j = 2,NlogTs-1
               
                  kmax = k_for_MAX_logP(j)
                  do k = kmax,NlogPs-1
                     si% Z(k,j) = 
     >                  (si% Z(k,j) 
     >                  + si% Z(k,j-1) 
     >                  + si% Z(k,j+1) 
     >                  + si% Z(k-1,j) 
     >                  + si% Z(k+1,j)) / 5
                  end do
                  si% Z(NlogPs,j) = si% Z(NlogPs-1,j)

                  kmin = k_for_MIN_logP(j)
                  do k = 2,kmin-1
                     si% Z(k,j) = 
     >                  (si% Z(k,j) 
     >                  + si% Z(k,j-1) 
     >                  + si% Z(k,j+1) 
     >                  + si% Z(k-1,j) 
     >                  + si% Z(k+1,j)) / 5
                  end do
                  si% Z(1,j) = si% Z(2,j)
                  
               end do
            end do
            
            si% f(1,:,:) = si% Z(:,:)

         end subroutine fill_and_smooth
         
         
         subroutine extra_smoothing(si,n,logT_lo,logT_hi,logP_lo,logP_hi)
            type (scvh_info) :: si
            integer,intent(in) :: n
            double precision,intent(in) :: logT_lo,logT_hi,logP_lo,logP_hi
            double precision :: Z(NlogPs,NlogTs)
            integer :: i,j,k
            Z(:,:) = si% f(1,:,:)
            do i = 1,n
               do j = 2,NlogTs-1
                  if (logTs(j) < logT_lo .or. logTs(j) > logT_hi) cycle
                  do k = 2,NlogPs-1
                     if (logPs(k) < logP_lo .or. logPs(k) > logP_hi) cycle
                     Z(k,j) = 
     >                  (Z(k,j) 
     >                  + Z(k,j-1) 
     >                  + Z(k,j+1) 
     >                  + Z(k-1,j) 
     >                  + Z(k+1,j)) / 5
                  end do
               end do
            end do                        
            si% f(1,:,:) = Z(:,:)            
         end subroutine extra_smoothing
                  
     
      end subroutine read_data_for_scvh


      subroutine read_file_for_scvh(data_dir,filename,n,np,nt,
     >      compv1,compv2,denlog,slog,ulog)
         character (len=*) :: data_dir,filename
         integer,intent(in) :: n,np,nt
         type (scvh_info) ::compv1,compv2,denlog,slog,ulog
         integer :: i,j,jj,nps,k,kk
         double precision :: tlg,plg,c1,c2,c3,c4,c5
         character (len=256) :: full_name
         write(full_name,'(3a)') trim(data_dir),'/',trim(filename)
         write(*,*) 'read ' // trim(full_name)
         open(unit=1,file=trim(full_name),action='read')
         do i=1,nt
            kk = 0
            read(1,*) tlg,nps
            if (abs(tlg - logTs(i)) > 1d-6) then
               write(*,*) 'tlg',i,tlg
               call do_stop('bad logT in table')
            end if
            do jj=1,nps
               read(1,*) plg,c1,c2,c3,c4,c5
               k = 0
               do j=1,np
                  if (abs(logPs(j)-plg) < 1d-6) then
                     k = j; exit
                  end if
               end do
               if (jj == 1) k_for_MIN_logP(i) = k
               if (k == 0) then
                  write(*,*) 'SCVH: failed to find pressure in logPs',plg
                  write(*,*) 'tlg',i,tlg
                  call do_stop('bad pressure in table')
               end if
               if (kk > 0 .and. k > kk+1) then ! check for internal gaps in pressures
                  write(*,'(a30,2i6,4x,e20.10)') 'SCVH: mistake in sequence of pressures',k,kk,plg
                  write(*,*) 'tlg',i,tlg
                  call do_stop('mistake in pressure in table')
               end if
               if (c1 < 1d-8) then
                  compv1% Z(k,i) = 0
               else
                  compv1% Z(k,i) = c1
               end if
               if (c2 < 1d-8) then
                  compv2% Z(k,i) = 0
               else
                  compv2% Z(k,i) = c2
               end if
               denlog% Z(k,i) = c3
               slog% Z(k,i) = c4
               ulog% Z(k,i) = c5
               kk = k
            enddo
            k_for_MAX_logP(i) = kk
         enddo
         close(unit=1)
         write(*,*) 'close ' // trim(filename)

      end subroutine read_file_for_scvh




      subroutine write_file_for_scvh(data_dir,filename,n,np,nt,
     >      compv1,compv2,denlog,slog,ulog)
         character (len=*) :: data_dir,filename
         integer,intent(in) :: n,np,nt
         type (scvh_info) ::compv1,compv2,denlog,slog,ulog
         integer :: i,j,jj,nps,k,kk
         double precision :: tlg,plg,c1,c2,c3,c4,c5
         logical :: lowT
         character (len=256) :: full_name
         write(full_name,'(3a)') trim(data_dir),'/',trim(filename)
         write(*,*) 'write ' // trim(full_name)
         open(unit=1,file=trim(full_name),action='write')
         do i=1,nt
            lowT = (logTs(i) < 4.05d0)
            if (lowT) then
               nps = 99
            else
               nps = 76
            end if
            write(1,'(f5.2,i4)') logTs(i),nps
            do k=NlogPs-nps+1,NlogPs
               if (lowT .and. logPs(k) < 3.99d0) then
                  c2 = compv1% Z(k,i)
                  c1 = compv2% Z(k,i)
               else
                  c1 = compv1% Z(k,i)
                  c2 = compv2% Z(k,i)
               end if
               c3 = denlog% Z(k,i)
               c4 = slog% Z(k,i)
               c5 = ulog% Z(k,i)
               write(1,'(f6.2,2(1pe14.5),3(0pf10.4))') logPs(k),c1,c2,c3,c4,c5
            enddo
         enddo
         close(unit=1)
         write(*,*) 'close ' // trim(filename)

      end subroutine write_file_for_scvh


      integer function locate_logT(logT)
         double precision,intent(in) :: logT
         integer :: jT,j
         jT = NlogTs-1
         do j=2,NlogTs-1
            if (logT < logTs(j)) then
               jT = j-1; exit
            end if
         end do
         locate_logT = jT
      end function locate_logT


      subroutine interp_densities(logT,logP,pres,den_h,den_he,info)
         double precision,intent(in) :: logT,logP,pres
         double precision,intent(out) :: den_h,den_he
         integer,intent(out) :: info
         double precision :: 
     1      ener_h,entr_h,dddt_cp_h,dddpress_ct_h,dsdt_cp_h,dsdpress_ct_h,
     1      ener_he,entr_he,dddt_cp_he,dddpress_ct_he,dsdt_cp_he,dsdpress_ct_he,
     1      xnh,dxnh_dlogT,dxnh_dlogP,
     >      xnh2,dxnh2_dlogT,dxnh2_dlogP,
     >      xnhe,dxnhe_dlogT,dxnhe_dlogP,
     >      xnhep,dxnhep_dlogT,dxnhep_dlogP
         logical,parameter :: only_densities = .true.,search_for_SCVH = .false.
         call interp_vals_bicub(
     1      only_densities,search_for_SCVH,
     1      den_h,ener_h,entr_h,dddt_cp_h,dddpress_ct_h,dsdt_cp_h,dsdpress_ct_h,
     1      den_he,ener_he,entr_he,dddt_cp_he,dddpress_ct_he,dsdt_cp_he,dsdpress_ct_he,
     1      xnh,dxnh_dlogT,dxnh_dlogP,
     >      xnh2,dxnh2_dlogT,dxnh2_dlogP,
     >      xnhe,dxnhe_dlogT,dxnhe_dlogP,
     >      xnhep,dxnhep_dlogT,dxnhep_dlogP,
     >      logT,logP,pres,info)
      end subroutine interp_densities
      

      subroutine interp_vals_bicub(
     1      only_densities,search_for_SCVH,
     1      den_h,ener_h,entr_h,dddt_cp_h,dddpress_ct_h,dsdt_cp_h,dsdpress_ct_h,
     1      den_he,ener_he,entr_he,dddt_cp_he,dddpress_ct_he,dsdt_cp_he,dsdpress_ct_he,
     1      xnh,dxnh_dlogT,dxnh_dlogP,
     >      xnh2,dxnh2_dlogT,dxnh2_dlogP,
     >      xnhe,dxnhe_dlogT,dxnhe_dlogP,
     >      xnhep,dxnhep_dlogT,dxnhep_dlogP,
     >      logT,logP,pres,info)
         logical,intent(in) :: only_densities,search_for_SCVH
         double precision,intent(out) :: 
     1      den_h,ener_h,entr_h,dddt_cp_h,dddpress_ct_h,dsdt_cp_h,dsdpress_ct_h,
     1      den_he,ener_he,entr_he,dddt_cp_he,dddpress_ct_he,dsdt_cp_he,dsdpress_ct_he,
     1      xnh,dxnh_dlogT,dxnh_dlogP,
     >      xnh2,dxnh2_dlogT,dxnh2_dlogP,
     >      xnhe,dxnhe_dlogT,dxnhe_dlogP,
     >      xnhep,dxnhep_dlogT,dxnhep_dlogP
         double precision,intent(in) :: logT,logP,pres
         integer,intent(out) :: info

         integer,parameter :: NIP = 1
         double precision :: XI(NIP),YI(NIP),ZI(NIP),max_logP,d_dlogT,d_dlogP
         integer :: j,jT,kP
         
         
         include 'formats'

         info = 0
         
         XI(1) = logP
         YI(1) = logT
         
         jT = locate_logT(logT)
         kP = k_for_max_logP(jT)
         max_logP = logPs(kP)
         if (logPs(k_for_max_logP(jT+1)) < max_logP) then
            kP = k_for_max_logP(jT+1)
            max_logP = logPs(kP)
         end if

         den_h = 10.0**interp_value(num_h_pts,denlog_h_si); if (info /= 0) return
         dddpress_ct_h = d_dlogP
         dddt_cp_h = d_dlogT
                
         den_he = 10.0**interp_value(num_he_pts,denlog_he_si); if (info /= 0) return
         dddpress_ct_he = d_dlogP
         dddt_cp_he = d_dlogT
         
         if (only_densities) return
         
         !..hydrogen section

         !..number concentration of h2 molecules
         xnh2 = interp_value(num_h_pts,compv1_h_si); if (info /= 0) return
         if (xnh2 > 1) then
            xnh2 = 1
            dxnh2_dlogP = 0
            dxnh2_dlogT = 0
         else if (xnh2 < 0) then
            xnh2 = 0
            dxnh2_dlogP = 0
            dxnh2_dlogT = 0
         else
            dxnh2_dlogP = d_dlogP
            dxnh2_dlogT = d_dlogT
         end if

         !..number concentration of neutral h atoms 
         xnh = interp_value(num_h_pts,compv2_h_si); if (info /= 0) return
         if (xnh > 1) then
            xnh = 1
            dxnh_dlogP = 0
            dxnh_dlogT = 0
         else if (xnh < 0) then
            xnh = 0
            dxnh_dlogP = 0
            dxnh_dlogT = 0
         else
            dxnh_dlogP = d_dlogP
            dxnh_dlogT = d_dlogT
         end if
         
         if (xnh + xnh2 > 1) then
            if (xnh > xnh2) then
               xnh = 1d0 - xnh2
               dxnh_dlogP = -dxnh2_dlogP
               dxnh_dlogT = -dxnh2_dlogT
            else
               xnh2 = 1d0 - xnh
               dxnh2_dlogP = -dxnh_dlogP
               dxnh2_dlogT = -dxnh_dlogT
            end if
         end if

         !.. entropy in erg/g/k
         entr_h = 10.0**interp_value(num_h_pts,slog_h_si); if (info /= 0) return
         dsdpress_ct_h = d_dlogP
         dsdt_cp_h = d_dlogT
         
         !..internal energy in erg/g
         ener_h = 10.0**interp_value(num_h_pts,ulog_h_si); if (info /= 0) return


         !..helium section

         !..number concentration of neutral helium atoms
         xnhe = interp_value(num_he_pts,compv1_he_si); if (info /= 0) return
         if (xnhe > 1) then
            xnhe = 1
            dxnhe_dlogP = 0
            dxnhe_dlogT = 0
         else if (xnhe < 0) then
            xnhe = 0
            dxnhe_dlogP = 0
            dxnhe_dlogT = 0
         else
            dxnhe_dlogP = d_dlogP
            dxnhe_dlogT = d_dlogT
         end if

         !..number concentration of he+ ions
         xnhep = interp_value(num_he_pts,compv2_he_si); if (info /= 0) return
         if (xnhep > 1) then
            xnhep = 1
            dxnhep_dlogP = 0
            dxnhep_dlogT = 0
         else if (xnhep < 0) then
            xnhep = 0
            dxnhep_dlogP = 0
            dxnhep_dlogT = 0
         else
            dxnhep_dlogP = d_dlogP
            dxnhep_dlogT = d_dlogT
         end if
         
         if (xnhe + xnhep > 1) then
            if (xnhe > xnhep) then
               xnhe = 1d0 - xnhep
               dxnhe_dlogP = -dxnhep_dlogP
               dxnhe_dlogT = -dxnhep_dlogT
            else
               xnhep = 1d0 - xnhe
               dxnhep_dlogP = -dxnhe_dlogP
               dxnhep_dlogT = -dxnhe_dlogT
            end if
         end if

         !.. entropy in erg/g/k
         entr_he = 10.0**interp_value(num_he_pts,slog_he_si); if (info /= 0) return
         dsdpress_ct_he = d_dlogP
         dsdt_cp_he = d_dlogT

         !..internal energy in erg/g
         ener_he = 10.0**interp_value(num_he_pts,ulog_he_si); if (info /= 0) return

               
         
         contains


         double precision function interp_value(n,si)
            use interp_2d_lib_db,only: interp_mkbicub_db,interp_evbicub_db
            use num_lib
            integer,intent(in) :: n
            type (scvh_info) :: si
            integer :: MD,ierr,iP
            double precision :: frac,v0,v1
            integer :: ibcxmin                   ! bc flag for x=xmin
            double precision :: bcxmin(NlogTs)       ! bc data vs. y at x=xmin
            integer :: ibcxmax                   ! bc flag for x=xmax
            double precision :: bcxmax(NlogTs)       ! bc data vs. y at x=xmax
            integer :: ibcymin                   ! bc flag for y=ymin
            double precision :: bcymin(NlogPs)       ! bc data vs. x at y=ymin
            integer :: ibcymax                   ! bc flag for y=ymax
            double precision :: bcymax(NlogPs)       ! bc data vs. x at y=ymax
            integer :: ict(6)                    ! code specifying output desired
            double precision :: fval(6)          ! output data
            include 'formats'
            
            ierr = 0
            if (si% needs_initialization) then
               ibcxmin = 0; bcxmin(:) = 0
               ibcxmax = 0; bcxmax(:) = 0
               ibcymin = 0; bcymin(:) = 0
               ibcymax = 0; bcymax(:) = 0
               call interp_mkbicub_db(
     >               logPs,NlogPs,logTs,NlogTs,si% f1,NlogPs,
     >               ibcxmin,bcxmin,ibcxmax,bcxmax,
     >               ibcymin,bcymin,ibcymax,bcymax,
     >               si% ilinx,si% iliny,ierr)
               if (ierr /= 0) then
                  return
                  write(*,*) 'scvh: failed in interp_mkbicub_db ' // si% name
                  stop 1
               end if
            end if
            ict(1:3) = 1; ict(4:6) = 0
            call interp_evbicub_db(
     >            XI(1),YI(1),logPs,NlogPs,logTs,NlogTs,si% ilinx,si% iliny,si% f1,NlogPs,ict,fval,info)
            interp_value = fval(1)
            d_dlogP = fval(2)
            d_dlogT = fval(3)
            if (info /= 0) then
               return
               write(*,*)
               write(*,2) 'logT',jT,logT
               write(*,1) 'logP',logP
               write(*,2) 'min_logP',1,logPs(1)
               write(*,2) 'max_logP',kP,max_logP
               write(*,*) 'scvh: failed in interp_evbicub_db ' // si% name
               stop 1
               interp_value = 0
            end if
            si% needs_initialization = .false.
            
         end function interp_value
      
      
      end subroutine interp_vals_bicub


      subroutine do_stop(str)
         character (len=*) :: str
         write(*,*) trim(str)
         stop 1
      end subroutine do_stop


      subroutine interpolate_scvh(
     >               include_radiation,search_for_SCVH,
     >               logT,logRho,T,Rho,xmassh1,
     >               logPgas,logE,logS,chiRho,chiT,
     >               Cp,Cv,dE_dRho,dS_dT,dS_dRho,
     >               mu,gamma1,gamma3,grad_ad,logNe,
     >               info)

      use num_lib,only: safe_root_without_brackets
      implicit none
      save

!..this routine interpolates the scvh tables

      logical,intent(in) :: include_radiation
      logical,intent(in) :: search_for_SCVH
      double precision,intent(inout) :: logT,logRho,T,Rho,xmassh1
      double precision,intent(out) ::  
     >               logPgas,logE,logS,chiRho,chiT,
     >               Cp,Cv,dE_dRho,dS_dT,dS_dRho,
     >               mu,gamma1,gamma3,grad_ad,logNe
      integer,intent(out) :: info ! returned = 0 if AOK
     
      double precision dpressdd,dpressdt,ener,dedd,
     >      xnh,dxnh_dlogT,dxnh_dlogP,
     >      xnh2,dxnh2_dlogT,dxnh2_dlogP,
     >      xnhe,dxnhe_dlogT,dxnhe_dlogP,
     >      xnhep,dxnhep_dlogT,dxnhep_dlogP,
     >      xnhp,xnhepp,P,entr,dsdd,dsdt,xtra,
     >      dxdd,dxdt,Ne,log_free_e,log_free_e0,log_free_e1,kt,theta
      
      double precision xmasshe4,dedt

!..local variables

      double precision prad,erad,srad

      double precision logP,inv_Rho,inv_T,inv_P,small_value,Rho_min
      parameter        (small_value = 1.0d-16)
      parameter        (Rho_min = 1.0d-10)

!..for hydrogen
      double precision den_h,ener_h,entr_h,dddt_cp_h,dddpress_ct_h,dsdt_cp_h,
     >      dsdpress_ct_h,dpressdd_h,dpressdt_h,dedd_h,dedt_h,dsdd_h,dsdt_h


!..for helium
      double precision den_he,ener_he,entr_he,dddt_cp_he,dddpress_ct_he,dsdt_cp_he,
     >      dsdpress_ct_he,dpressdd_he,dpressdt_he,dedd_he,dedt_he,dsdd_he,dsdt_he


!..for the mixture
      double precision beta,gama,delt,smix,d_smix_dT,d_smix_dP,dddt_cp_hhe,dddpress_ct_hhe,dsdt_cp_hhe,
     1                 dsdpress_ct_hhe,dtdpress_cs_hhe,dpressdd_hhe,dpressdt_hhe,dedd_hhe,dedt_hhe,dsdd_hhe,dsdt_hhe
      double precision tiny
      parameter        (tiny = 1.0d-30)


!..for the search
      integer          ii,itmax
      parameter        (itmax = 200)
      double precision logP_new,pmax,pmin,den_calc,func,denom,eostol
      parameter        (eostol = 1.0d-9)


!..for the photons
      integer          radmult
      parameter        (radmult = 1)
      double precision dpressraddd,dpressraddt,deraddd,deraddt,dsraddd,dsraddt

!..constants
      double precision clight,ssol,asol,asoli3,avo,kerg,xka,mh1,mhe4,third
      parameter        (clight  = 2.99792458d10,ssol    = 5.67051d-5,
     1                  asol    = 4.0d0 * ssol / clight,asoli3  = asol/3.0d0,
     1                  avo     = 6.0221367d23,kerg    = 1.380658d-16,
     1                  xka     = kerg*avo,mh1     = 1.67357d-24,
     1                  mhe4    = 6.646442d-24,third   = 1.0d0/3.0d0)

      logical,parameter :: pure_splines = .true.,DT_flag = .true.
      
      integer,parameter :: lrpar=3,lipar=0
      integer, target :: ipar_array(lipar)
      double precision, target :: rpar_array(lrpar)
      integer, pointer :: ipar(:)
      double precision, pointer :: rpar(:)
      integer          iat,jat,i,j,imax
      double precision dt,dt2,dti,dt2i,dpress,dpress2,dpressi,dpress2i,xt,xp,mxt,mxp,logRho_new,dlogRho_dlogP,del_logP,under,
     1                 si0t,si1t,si2t,si0mt,si1mt,si2mt,si0p,si1p,si2p,si0mp,si1mp,si2mp,logRho_test,
     1                 logP_guess,epslogP,epslogRho,x1,x3,y1,y3,dfdlogP,
     1                 z,fi_h(36),fi_he(36),w0t,w1t,w2t,w0mt,w1mt,w2mt,w0p,w1p,w2p,w0mp,w1mp,w2mp

      include 'formats'
      info = 0
      ipar => ipar_array
      rpar => rpar_array

      z=0; w0t=0; w1t=0; w2t=0; w0mt=0; w1mt=0; w2mt=0
      w0p=0; w1p=0; w2p=0; w0mp=0; w1mp=0; w2mp=0; i=0; j=0

!..initialize

      if (Rho < Rho_min) then
         Rho = Rho_min; logRho = log10(Rho)
      end if
      
      xnh2    = small_value
      xnh     = small_value
      xnhe    = small_value
      xnhep   = small_value
      dpressdd    = small_value
      dpressdt    = small_value
      dedd    = small_value
      dedt    = small_value
      dsdd    = small_value
      dsdt    = small_value
      P    = small_value
      ener    = small_value
      entr    = small_value
      inv_Rho = 1.0d0/Rho
      inv_T = 1.0d0/T

      xmasshe4 = 1 - xmassh1

!..radiation section:
      prad    = 0.0
      dpressraddd = 0.0
      dpressraddt = 0.0

      erad    = 0.0
      deraddd = 0.0
      deraddt = 0.0

      srad    = 0.0
      dsraddd = 0.0
      dsraddt = 0.0

      if (radmult .ne. 0) then
       prad    = asoli3 * T**4
       dpressraddd = 0.0
       dpressraddt = 4.0d0 * prad * inv_T

       erad    = 3.0d0 * prad * inv_Rho
       deraddd = -erad * inv_Rho
       deraddt = 3.0d0 * dpressraddt * inv_Rho

       srad    = (prad*inv_Rho + erad) * inv_T
       dsraddd = ((dpressraddd - prad*inv_Rho)*inv_Rho + deraddd)*inv_T
       dsraddt = (dpressraddt*inv_Rho + deraddt - srad)*inv_T
      end if

!..find pressure that gives the desired density

      logP_guess = min(logP_max - 0.1d0,
     >         log10(kerg/mh1) - log10(4/(3+5*xmassh1)) + logRho + logT) ! classical perfect gas value
      rpar(1) = logT
      rpar(2) = logRho
      rpar(3) = xmassh1
      imax = 100
      epslogP = 1d-14
      epslogRho = 1d-10
      logP_new = safe_root_without_brackets(
     >      fscvh,logP_guess,dlogP,imax/2,imax,epslogP,epslogRho,lrpar,rpar,lipar,ipar,info)
      if (dbg) write(*,1) 'logP_new',logP_new
      if (info /= 0) then
         info = -2
         if (dbg) write(*,*) 'search for P to give rho failed.  reject and retry smaller rho.'
         return
      end if

      if (abs(logP_new-logP_max) < 1d-5) then ! have hit edge of table.  reject it.
         info = -2
         if (dbg) write(*,*) 'search for P to give rho hit max logP.  reject and retry smaller rho.'
         return
      end if

      logP   = logP_new
      P     = 10.0**logP

!..now get the rest of the hydrogen and helium quantities
      call get_values(.false.)
      if (info /= 0) then
         return
         write(*,1) 'scvh: failed in get_values,logRho,logT',logRho,logT
         stop 1
      end if
      
      if (dsdt_cp_h <= 0 .or. dsdt_cp_he <= 0) then ! indicates off table
         info = -1
         return

         write(*,*) 'off table'
         write(*,*)
         write(*,1) 'den_h',den_h
         write(*,1) 'ener_h',ener_h
         write(*,1) 'dddt_cp_h',dddt_cp_h
         write(*,1) 'dddpress_ct_h',dddpress_ct_h
         write(*,1) 'dsdt_cp_h',dsdt_cp_h
         write(*,1) 'dsdpress_ct_h',dsdpress_ct_h
         write(*,*)
         write(*,1) 'den_he',den_he
         write(*,1) 'ener_he',ener_he
         write(*,1) 'dddt_cp_he',dddt_cp_he
         write(*,1) 'dddpress_ct_he',dddpress_ct_he
         write(*,1) 'dsdt_cp_he',dsdt_cp_he
         write(*,1) 'dsdpress_ct_he',dsdpress_ct_he
         write(*,*) 'scvh'
         stop 1
         info = -1
         return
      end if

      logRho_test = log10(1.0d0/(xmassh1/den_h + xmasshe4/den_he))
      if (abs(logRho - logRho_test) > 1d-6) then
         info = -1
         return
         write(*,1) 'abs(logRho - logRho_test)',abs(logRho - logRho_test),logRho,logT
         write(*,*) 'bad match for logRho'
         write(*,*)
         write(*,1) 'den_h',den_h
         write(*,1) 'den_he',den_he
         write(*,1) 'logRho - logRho_test',logRho - logRho_test
         write(*,1) 'logRho_test',logRho_test
         write(*,1) 'logRho',logRho
         write(*,1) 'logT',logT
         write(*,*) 'scvh'
         stop 1
      end if

      call entropy_of_mixing(
     >     xmassh1,xmasshe4,T,P,
     >     xnh,dxnh_dlogT,dxnh_dlogP,
     >     xnh2,dxnh2_dlogT,dxnh2_dlogP,
     >     xnhe,dxnhe_dlogT,dxnhe_dlogP,
     >     xnhep,dxnhep_dlogT,dxnhep_dlogP,
     >     smix,d_smix_dT,d_smix_dP)
      
      !smix = 0; d_smix_dT = 0; d_smix_dP = 0
      
      ener = xmassh1*ener_h + xmasshe4*ener_he  ! eqn 40
      entr = xmassh1*entr_h + xmasshe4*entr_he + smix ! eqn 41
      if (.false.) then
         write(*,1) 'entr_h',entr_h
         write(*,1) 'entr_he',entr_he
         write(*,1) 'smix',smix
         write(*,1) 'd_smix_dP',d_smix_dP
         write(*,1) 'd_smix_dT',d_smix_dT
         write(*,*)
      end if

! form logarithmic derivatives for the mixture

      dddt_cp_hhe = Rho*(xmassh1/den_h*dddt_cp_h + xmasshe4/den_he*dddt_cp_he)  ! eqn 43
      dddpress_ct_hhe = Rho*(xmassh1/den_h*dddpress_ct_h + xmasshe4/den_he*dddpress_ct_he)  ! eqn 44

      dsdt_cp_hhe=(xmassh1*entr_h*dsdt_cp_h + xmasshe4*entr_he*dsdt_cp_he + T*d_smix_dT)/entr ! eqn 45
      dsdpress_ct_hhe=(xmassh1*entr_h*dsdpress_ct_h + xmasshe4*entr_he*dsdpress_ct_he + P*d_smix_dP)/entr ! eqn 46

      dtdpress_cs_hhe = -dsdpress_ct_hhe/dsdt_cp_hhe  ! eqn 47


!..delog the derivatives
      inv_P = 1.0d0/P
      dddt_cp_h  = Rho * inv_T     * dddt_cp_h
      dddpress_ct_h  = Rho * inv_P    * dddpress_ct_h
      dsdt_cp_h  = entr_h * inv_T  * dsdt_cp_h
      dsdpress_ct_h  = entr_h * inv_P * dsdpress_ct_h

      dddt_cp_he = Rho * inv_T      * dddt_cp_he
      dddpress_ct_he = Rho * inv_P     * dddpress_ct_he
      dsdt_cp_he = entr_he * inv_T  * dsdt_cp_he
      dsdpress_ct_he = entr_he * inv_P * dsdpress_ct_he

      dddt_cp_hhe = Rho * inv_T   * dddt_cp_hhe
      dddpress_ct_hhe = Rho * inv_P  * dddpress_ct_hhe
      dsdt_cp_hhe = entr * inv_T  * dsdt_cp_hhe
      dsdpress_ct_hhe = entr * inv_P * dsdpress_ct_hhe
      dtdpress_cs_hhe = T * inv_P  * dtdpress_cs_hhe


!..form the usual thermodynamic derivatives 
!..for hydrogen

!..d(P)/ d(den)|t
      dpressdd_h  = 1.0d0/dddpress_ct_h

!..d(entr)/d(den)|t
      dsdd_h = dddt_cp_h/(Rho*Rho*dddpress_ct_h)

!..d(P)/d(temp)|d
      dpressdt_h = -Rho*Rho*dsdd_h

!..d(ener)/d(rho)|t
      dedd_h = (P - T * dpressdt_h)*inv_Rho*inv_Rho

!..d(entr)/d(temp)|d
      dsdt_h = dsdpress_ct_h * dpressdt_h + dsdt_cp_h

!..d(ener)/d(temp)|d
      dedt_h = T * dsdt_h



!..for helium
!..d(P)/ d(den)|t
      dpressdd_he  = 1.0d0/dddpress_ct_he

!..d(entr)/d(den)|t
      dsdd_he = dddt_cp_he/(Rho*Rho*dddpress_ct_he)

!..d(P)/d(temp)|d
      dpressdt_he = -Rho*Rho*dsdd_he

!..d(ener)/d(rho)|t
      dedd_he = (P - T * dpressdt_he)*inv_Rho*inv_Rho

!..d(entr)/d(temp)|d
      dsdt_he = dsdpress_ct_he * dpressdt_he + dsdt_cp_he

!..d(ener)/d(temp)|d
      dedt_he = T * dsdt_he


!..for the mixture
!..d(P)/ d(den)|t
      dpressdd_hhe  = 1.0d0/dddpress_ct_hhe

!..d(entr)/d(den)|t
      dsdd_hhe = dddt_cp_hhe/(Rho*Rho*dddpress_ct_hhe)

!..d(P)/d(temp)|d
      dpressdt_hhe = -Rho*Rho*dsdd_hhe

!..d(ener)/d(rho)|t
      dedd_hhe = (P - T * dpressdt_hhe)*inv_Rho*inv_Rho

!..d(entr)/d(temp)|d
      dsdt_hhe = dsdpress_ct_hhe * dpressdt_hhe + dsdt_cp_hhe

!..d(ener)/d(temp)|d
      dedt_hhe = T * dsdt_hhe

!..sum the components

      logPgas = log10(P) ! store this before add prad

      if (include_radiation) then
         P = P     + prad
         dpressdd = dpressdd_hhe + dpressraddd
         dpressdt = dpressdt_hhe + dpressraddt
   
         ener = ener     + erad
         dedd = dedd_hhe + deraddd
         dedt = dedt_hhe + deraddt
         
         entr = entr     + srad
         dsdd = dsdd_hhe + dsraddd
         dsdt = dsdt_hhe + dsraddt
      else
         dpressdd = dpressdd_hhe
         dpressdt = dpressdt_hhe
   
         dedd = dedd_hhe
         dedt = dedt_hhe
   
         dsdd = dsdd_hhe
         dsdt = dsdt_hhe
      end if

      ! store the output

      logE = log10(ener)
      logS = log10(entr)
            
      Cv = dedt
      dE_dRho = dedd
      dS_dT = dsdt
      dS_dRho = dsdd

      chiRho = dpressdd * Rho / P
      chiT = dpressdt * T / P
      
      gamma3 = 1 + dpressdt / (Rho * dedt)   
      grad_ad = dtdpress_cs_hhe/(T * inv_P)
      gamma1 = (gamma3 - 1) / grad_ad ! C&G 9.88 & 9.89
      
      Cp = Cv + P * chiT**2 / (Rho * T * chiRho) ! C&G 9.86
      


      xnhp = max(0d0,min(1d0,1 - (xnh2 + xnh)))
      xnhepp = max(0d0,min(1d0,1 - (xnhep + xnhe)))
      mu = xmassh1 * (xnhp / 2 + xnh + 2 * xnh2) + 4 * (1 - xmassh1) * (xnhepp / 3 + xnhep / 2 + xnhe)
      
      Ne = Rho * avo * (xmassh1 * xnhp + (1 - xmassh1) * (xnhep + 2 * xnhepp) / 4 )
      if (Ne < 1d-99) Ne = 1d-99
      logNe = log10(Ne)
      if (logNe - logRho < 17d0) logNe = logRho + 17d0
         ! put a lower limit on the electron abundance to avoid spurious
         ! bicubic interpolations at very low T

      if (.false. .and. gamma1 <= 0) then
         write(*,1) 'scvh_core gamma1',gamma1
         write(*,1) 'logT',logT
         write(*,1) 'logRho',logRho
         write(*,1) 'logPgas',logPgas
         write(*,1) 'xmassh1',xmassh1
         write(*,*)
      end if
      

      if (.false.) then
         write(*,1) 'logT',logT
         write(*,1) 'logRho',logRho
         write(*,*)


         !dsdd_hhe = dddt_cp_hhe/(Rho*Rho*dddpress_ct_hhe)
         write(*,1) 'dddpress_ct_hhe',dddpress_ct_hhe
         write(*,1) 'dddt_cp_hhe',dddt_cp_hhe
         write(*,1) 'dsdd_hhe',dsdd_hhe
         write(*,*)
         write(*,1) 'gamma3',gamma3
         write(*,1) 'dpressdt',dpressdt
         write(*,1) 'dedt',dedt
         write(*,1) 'Rho',Rho
         write(*,*)
         write(*,1) 'T',T
         write(*,1) 'dsdt_hhe',dsdt_hhe
         write(*,*)
         write(*,1) 'dsdpress_ct_hhe',dsdpress_ct_hhe
         write(*,1) 'dpressdt_hhe',dpressdt_hhe
         write(*,1) 'dsdt_cp_hhe',dsdt_cp_hhe
         write(*,*)
         write(*,1) 'dedd',dedd
         write(*,1) 'dpressdt_hhe',dpressdt_hhe
         write(*,1) 'P',P
         write(*,1) 'T',T
         write(*,1) 'inv_Rho',inv_Rho
         write(*,*)
      end if
            
      contains
      
      subroutine get_values(only_densities)
         logical,intent(in) :: only_densities
         info = 0
         call interp_vals_bicub(
     1      only_densities,search_for_SCVH,
     1      den_h,ener_h,entr_h,dddt_cp_h,dddpress_ct_h,dsdt_cp_h,dsdpress_ct_h,
     1      den_he,ener_he,entr_he,dddt_cp_he,dddpress_ct_he,dsdt_cp_he,dsdpress_ct_he,
     1      xnh,dxnh_dlogT,dxnh_dlogP,
     >      xnh2,dxnh2_dlogT,dxnh2_dlogP,
     >      xnhe,dxnhe_dlogT,dxnhe_dlogP,
     >      xnhep,dxnhep_dlogT,dxnhep_dlogP,
     >      logT,logP,P,info)
      end subroutine get_values
      
      end subroutine interpolate_scvh
      


      double precision function fscvh(logP,dfdlogP,lrpar,rpar,lipar,ipar,ierr)
         ! returns with ierr = 0 if was able to evaluate f and df/dx at x
         ! if df/dx not available,it is okay to set it to 0
         integer,intent(in) :: lrpar,lipar
         double precision,intent(in) :: logP
         double precision,intent(out) :: dfdlogP
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         double precision, intent(inout), pointer :: rpar(:) ! (lrpar)
         integer,intent(out) :: ierr
         
         double precision :: logT,P,Rho,logRho_new,logRho_target,dddpress,
     1      den_h,ener_h,entr_h,dddt_cp_h,dddpress_ct_h,dsdt_cp_h,dsdpress_ct_h,
     1      den_he,ener_he,entr_he,dddt_cp_he,dddpress_ct_he,dsdt_cp_he,dsdpress_ct_he,
     1      xnh,dxnh_dlogT,dxnh_dlogP,
     >      xnh2,dxnh2_dlogT,dxnh2_dlogP,
     >      xnhe,dxnhe_dlogT,dxnhe_dlogP,
     >      xnhep,dxnhep_dlogT,dxnhep_dlogP,xmassh1,xmasshe4
         logical :: only_densities,search_for_SCVH
         include 'formats'
         if (logP < logP_min .or. logP > logP_max) then
            ierr = -1
            if (dbg) write(*,*)
            if (dbg) write(*,1) 'logP guess out of bounds',logP
            return
         end if
         logT = rpar(1)
         logRho_target = rpar(2)
         xmassh1 = rpar(3)
         xmasshe4 = 1 - xmassh1
         P = 10**logP
         ierr = 0
         only_densities = .false.
         search_for_SCVH = .false.
         if (dbg) write(*,*)
         if (dbg) write(*,1) 'logP guess',logP
         call interp_vals_bicub(
     1      only_densities,search_for_SCVH,
     1      den_h,ener_h,entr_h,dddt_cp_h,dddpress_ct_h,dsdt_cp_h,dsdpress_ct_h,
     1      den_he,ener_he,entr_he,dddt_cp_he,dddpress_ct_he,dsdt_cp_he,dsdpress_ct_he,
     1      xnh,dxnh_dlogT,dxnh_dlogP,
     >      xnh2,dxnh2_dlogT,dxnh2_dlogP,
     >      xnhe,dxnhe_dlogT,dxnhe_dlogP,
     >      xnhep,dxnhep_dlogT,dxnhep_dlogP,
     >      logT,logP,P,ierr)
         if (ierr /= 0) then
            return
            write(*,*) 'fscvh failed in interp_vals'
            stop 1
         end if
         Rho = 1.0d0/(xmassh1/den_h + xmasshe4/den_he)
         logRho_new = log10(Rho)
         if (dbg) write(*,1) 'new logRho',logRho_new
         fscvh = logRho_new - logRho_target
         dfdlogP = Rho*(xmassh1/den_h*dddpress_ct_h + xmasshe4/den_he*dddpress_ct_he)
         if (dbg) write(*,1) 'DDDPT',dfdlogP
         if (dbg) write(*,1) 'DP',fscvh/dfdlogP
         if (dbg) write(*,1) 'DP/PR',fscvh/dfdlogP/logP
         if (dbg) write(*,1) 'fscvh',fscvh
      end function fscvh
      


      subroutine entropy_of_mixing(
     >     xmassh1,xmasshe4,T,P,
     >     xnh_in,dxnh_dlogT,dxnh_dlogP,
     >     xnh2_in,dxnh2_dlogT,dxnh2_dlogP,
     >     xnhe_in,dxnhe_dlogT,dxnhe_dlogP,
     >     xnhep_in,dxnhep_dlogT,dxnhep_dlogP,
     >     smix,d_smix_dT,d_smix_dP)
         implicit none
         double precision,intent(in) :: xmassh1,xmasshe4,T,P,
     >     xnh_in,dxnh_dlogT,dxnh_dlogP,
     >     xnh2_in,dxnh2_dlogT,dxnh2_dlogP,
     >     xnhe_in,dxnhe_dlogT,dxnhe_dlogP,
     >     xnhep_in,dxnhep_dlogT,dxnhep_dlogP
         double precision,intent(out) :: smix,d_smix_dT,d_smix_dP

         double precision,parameter :: tiny = 1d-14,ln10 = 2.30258509299405d0,
     >      mh1 = 1.67357d-24,mhe4 = 6.646442d-24,kerg = 1.3806504D-16,small = 1d-8
         double precision ::
     >      xnh,dxnh_dP,dxnh_dT,xnh2,dxnh2_dP,dxnh2_dT,
     >      xnhe,dxnhe_dP,dxnhe_dT,xnhep,dxnhep_dP,dxnhep_dT,
     >      beta,
     >      num,dnum_dP,dnum_dT,
     >      denom,ddenom_dP,ddenom_dT,
     >      gama,dgama_dP,dgama_dT,
     >      delt,ddelt_dP,ddelt_dT,        
     >      a,da_dP,da_dT,
     >      b1,db1_dT,db1_dP,
     >      b21,db21_dT,db21_dP,
     >      b22,db22_dT,db22_dP,
     >      b2,db2_dT,db2_dP,
     >      b31,db31_dT,db31_dP,
     >      b32,db32_dT,db32_dP,
     >      b331,db331_dT,db331_dP,
     >      b332,db332_dT,db332_dP,
     >      b33,db33_dT,db33_dP,
     >      b3,db3_dT,db3_dP,
     >      b,db_dT,db_dP

         include 'formats'
         
         if (xnh_in > 1d0) then
            xnh = 1d0
            dxnh_dP = 0
            dxnh_dT = 0
         else if (xnh_in > small) then
            xnh = xnh_in
            dxnh_dP = dxnh_dlogP/(P*ln10)
            dxnh_dT = dxnh_dlogT/(T*ln10)
         else
            xnh = 0
            dxnh_dP = 0
            dxnh_dT = 0
         end if
         
         if (xnh2_in > 1d0) then
            xnh2 = 1d0
            dxnh2_dP = 0
            dxnh2_dT = 0
         else if (xnh2_in > small) then
            xnh2 = xnh2_in
            dxnh2_dP = dxnh2_dlogP/(P*ln10)
            dxnh2_dT = dxnh2_dlogT/(T*ln10)
         else
            xnh2 = 0
            dxnh2_dP = 0
            dxnh2_dT = 0
         end if

         if (xnhe_in > 1d0) then
            xnhe = 1d0
            dxnhe_dP = 0
            dxnhe_dT = 0
         else if (xnhe_in > small) then
            xnhe = xnhe_in
            dxnhe_dP = dxnhe_dlogP/(P*ln10)
            dxnhe_dT = dxnhe_dlogT/(T*ln10)
         else
            xnhe = 0
            dxnhe_dP = 0
            dxnhe_dT = 0
         end if
         
         if (xnhep_in > 1d0) then
            xnhep = 1d0
            dxnhep_dP = 0
            dxnhep_dT = 0
         else if (xnhep_in > small) then
            xnhep = xnhep_in
            dxnhep_dP = dxnhep_dlogP/(P*ln10)
            dxnhep_dT = dxnhep_dlogT/(T*ln10)
         else
            xnhep = 0
            dxnhep_dP = 0
            dxnhep_dT = 0
         end if
      
         beta = (mh1 * (xmasshe4+tiny)) / (mhe4 * (xmassh1 + tiny)) 
      
         num = 1.5d0*(1 + xnh + 3*xnh2)
         dnum_dP = 1.5d0*(dxnh_dP + 3*dxnh2_dP)
         dnum_dT = 1.5d0*(dxnh_dT + 3*dxnh2_dT)

         denom = 1 + 2*xnhe + xnhep
         ddenom_dP = 2*dxnhe_dP + dxnhep_dP
         ddenom_dT = 2*dxnhe_dT + dxnhep_dT

         gama = num/denom
         dgama_dP = dnum_dP/denom - ddenom_dP*gama/denom
         dgama_dT = dnum_dT/denom - ddenom_dT*gama/denom

         num = 1.5d0*beta*gama*(2-2*xnhe-xnhep)
         dnum_dP = 1.5d0*beta*(dgama_dP*(2-2*xnhe-xnhep) + gama*(-2*dxnhe_dP-dxnhep_dP))
         dnum_dT = 1.5d0*beta*(dgama_dT*(2-2*xnhe-xnhep) + gama*(-2*dxnhe_dT-dxnhep_dT))

         denom = 1 - xnh2 - xnh
         if (denom <= 1d-5) then
            denom = 1d-5
            ddenom_dP = 0
            ddenom_dT = 0
         else
            ddenom_dP = - dxnh2_dP - dxnh_dP
            ddenom_dT = - dxnh2_dT - dxnh_dT
         end if

         delt = num/denom
         if (delt <= tiny) then
            delt = tiny
            ddelt_dP = 0
            ddelt_dT = 0
         else
            ddelt_dP = dnum_dP/denom - ddenom_dP*delt/denom
            ddelt_dT = dnum_dT/denom - ddenom_dT*delt/denom
         end if
      
         a = xmassh1/mh1*(2/(1 + xnh + 3*xnh2))
         da_dP = -a*(dxnh_dP + 3*dxnh2_dP)/(1 + xnh + 3*xnh2)
         da_dT = -a*(dxnh_dT + 3*dxnh2_dT)/(1 + xnh + 3*xnh2)
      
         b1 = log(1 + beta*gama)
         db1_dT = beta*dgama_dT/(1 + beta*gama)
         db1_dP = beta*dgama_dP/(1 + beta*gama)
      
         b21 = 0.5d0*(1-xnh2-xnh)
         db21_dT = 0.5d0*(-dxnh2_dT-dxnh_dT)
         db21_dP = 0.5d0*(-dxnh2_dP-dxnh_dP)
      
         b22 = log(1 + delt)
         db22_dT = ddelt_dT/(1 + delt)
         db22_dP = ddelt_dP/(1 + delt)
      
         b2 = b21*b22
         db2_dT = db21_dT*b22 + b21*db22_dT
         db2_dP = db21_dP*b22 + b21*db22_dP
      
         b31 = beta*gama
         db31_dT = beta*dgama_dT
         db31_dP = beta*dgama_dP
      
         b32 = log(1 + 1/(beta*gama))
         db32_dT = -dgama_dT/(gama*(1 + beta*gama))
         db32_dP = -dgama_dP/(gama*(1 + beta*gama))
      
         b331 = log(1 + 1/delt)
         db331_dT = -ddelt_dT/(delt*(1 + delt))
         db331_dP = -ddelt_dP/(delt*(1 + delt))
      
         b332 = (2-2*xnhe-xnhep)/3
         db332_dT = (-2*dxnhe_dT-dxnhep_dT)/3
         db332_dP = (-2*dxnhe_dP-dxnhep_dP)/3
      
         b33 = b331*b332
         db33_dT = db331_dT*b332 + b331*db332_dT
         db33_dP = db331_dP*b332 + b331*db332_dP
      
         b3 = b31*(b32 - b33)
         db3_dT = db31_dT*(b32 - b33) + b31*(db32_dT - db33_dT)
         db3_dP = db31_dP*(b32 - b33) + b31*(db32_dP - db33_dP)
      
         b = b1 - b2 + b3
         db_dT = db1_dT - db2_dT + db3_dT
         db_dP = db1_dP - db2_dP + db3_dP

         smix = kerg*a*b
         d_smix_dT = kerg*(a*db_dT + da_dT*b)
         d_smix_dP = kerg*(a*db_dP + da_dP*b)
         
         return
         
         if (P < 1d12) return
         write(*,1) 'T',T
         write(*,1) 'P',P
         write(*,1) 'smix',smix
         write(*,1) 'd_smix_dT',d_smix_dT
         write(*,1) 'a',a
         write(*,1) 'b',b
         write(*,1) 'db_dT',db_dT
         write(*,1) 'da_dT',da_dT
         write(*,1) 'a*db_dT',a*db_dT
         write(*,1) 'da_dT*b',da_dT*b
         write(*,1) 'db1_dT',db1_dT
         write(*,1) 'db2_dT',db2_dT
         write(*,1) 'db3_dT',db3_dT
         write(*,1) 'db21_dT*b22',db21_dT*b22
         write(*,1) 'b21*db22_dT',b21*db22_dT
         write(*,1) 'b21',b21
         write(*,1) 'b22',b22
         write(*,1) 'db21_dT',db21_dT
         write(*,1) 'db22_dT',db22_dT
         write(*,1) 'delt',delt
         write(*,1) 'xnh',xnh
         write(*,1) 'xnh2',xnh2
         write(*,1) 'xnhe',xnhe
         write(*,1) 'xnhep',xnhep
         !write(*,1) '',
         !write(*,1) '',
         !write(*,1) '',
         
!   >     ,dxnh_dlogT,dxnh_dlogP,
!   >     ,dxnh2_dlogT,dxnh2_dlogP,
!   >     ,dxnhe_dlogT,dxnhe_dlogP,
!   >     ,dxnhep_dlogT,dxnhep_dlogP
         write(*,*)

      end subroutine entropy_of_mixing




      end module scvh_core

