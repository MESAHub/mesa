! ***********************************************************************
!
!   Copyright (C) 2010  Aaron Dotter
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
      program table_merge
      use const_def, only: dp
      implicit none
!Combines three types of surface BC tables [ Pgas(Teff,logg) ] as listed below.
!Resulting tables have 14 columns (1=Teff; 2-14=Pgas@logg = -0.5 - 5.5 dex) and
!82 rows of Teff values + 2 header lines.  logg is provided from -0.5 (0.5) 5.5
!and Teff is provided as 2,000(250)13,000 + 14,000(1,000)50,000 following the C&K
!tables' Teff distribution. Only the 4 highest logg values are defined at all Teff
!values.  Other Teff values have lower maximum temperatures; empty table points
!are filled with zeros.  The program weights points in each table such that if a
!point is near the edge of the table, it is more "transparent" (has less weight).
!Eddington T-tau atm integrations are used as the first "layer" in this scheme as
!they provide the best coverage of the Teff-logg plane. C&K tables are layered on
!top of the Eddington table and then PHOENIX (only at cool Teff) are layered on
!top of the C&K tables.
!Three types of tables to work with are: 
!   1. Eddington T-tau to fill in holes   [ T:2,500-50,000K, logg: -0.5 - 5.5 (nT_Ttau=80) ]
!   2. Castelli & Kurucz model atmosphere [ T:3,500-50,000K, logg:  0.0 - 5.0 (nT_ck  =76) ]
!   3. PHOENIX model atmosphere           [ T:2,000-10,000K, logg: -0.5 - 5.5 (nT_phx =33) ]
      integer :: i, j, ii, j_max, i_max, jj
      integer, parameter :: ng=13, nT_phx=33, nT_ck=76, nT_Ttau=80, nT_all=82 !table dimensions
      integer, parameter :: table_version = 4
      integer :: zero_phx(ng), zero_ck(ng), ibound(ng) !locate the edges of the tables at constant logg
      character(len=64) :: f_phx, f_ck, f_Ttau, f_output, head_phx, head_Ttau !filenames, table headers
      real(dp) :: T_phx(nT_phx), T_ck(nT_ck), T_Ttau(nT_Ttau), T_all(nT_all),tran_tbl, tran_atm
      real(dp) :: P_phx(ng,nT_phx), P_ck(ng,nT_ck), P_Ttau(ng,nT_Ttau), P_all(ng,nT_all)
      real(dp) :: Pinterp, Pinterp_max, d0, d1
      real(dp), parameter :: Pmin = 1d-5 !used to set very small P to zero
      real(dp), parameter :: Pfill = 10 !used for interpolation in filling missing values

      if(COMMAND_ARGUMENT_COUNT()/=4) stop 'usage: ./table_merge [phoenix] [C&K] [T-tau] [output]'
      !filenames of input tables and output file'
      call GET_COMMAND_ARGUMENT(1,f_phx)
      call GET_COMMAND_ARGUMENT(2,f_ck)
      call GET_COMMAND_ARGUMENT(3,f_Ttau)
      call GET_COMMAND_ARGUMENT(4,f_output)

      !read PHX table
      open(1,file=trim(f_phx))
      read(1,'(a40)') head_phx
      read(1,*)
      do i=1,nT_phx
         read(1,*) T_phx(i), (P_phx(j,i),j=1,ng)
      enddo
      close(1)

      !read C&K table
      open(1,file=trim(f_ck))
      read(1,*)
      do i=1,nT_ck
         read(1,*) T_ck(i), (P_ck(j,i),j=1,ng)
      enddo
      close(1)
      P_ck(6,40)=0d0 !C&K tables have an inconvenient data point here

      !read T-tau table
      open(1,file=trim(f_Ttau))
      read(1,'(a40)') head_Ttau
      read(1,*)
      do i=1,nT_Ttau
         read(1,*) T_Ttau(i), (P_Ttau(j,i),j=1,ng)
      enddo
      close(1)

      !locate index of first zero for each logg, used in transparency function
      do j=1,ng
         i=1
         do while(P_ck(j,i) > Pmin .and. i < nT_ck)
            i=i+1
         enddo
         zero_ck(j)=i
         i=1
         do while(P_phx(j,i) > Pmin .and. i < nT_phx)
            i=i+1
         enddo
         zero_phx(j)=i
      enddo
      !FUDGE WARNING--reset zero for j=4,7 to the point at which T=7000K
      !atm integration gives very different result from model atmosphere
      !for T > 7000K.  Elsewhere the difference is not noticeable.
      do j=4,7
         i=1
         do while(T_ck(i)<7d3)
            i=i+1
         enddo
         zero_ck(j)=i
         i=1
         do while(T_phx(i)<7d3)
            i=i+1
         enddo
         zero_phx(j)=i
      enddo
      !end of FUDGE

      !initialize final table T, P with Eddington T-tau values, offset is 2
      P_all(:,:)=0d0
      T_all(1)=T_phx(1)
      T_all(2)=T_phx(2)
      do i=1,nT_Ttau
         ii=i+2
         T_all(ii)=T_Ttau(i)
         do j=1,ng
            P_all(j,ii) = P_Ttau(j,i)
         enddo
      enddo

      !add C&K pressures to the T-tau table, T offset is 4
      do i=1,nT_ck
         ii=i+6
         if(T_ck(i) /= T_all(ii)) stop 'T_ck != Tall'
         do j=1,ng
            tran_tbl = transparency(i,zero_ck(j))
            tran_atm = 1d0 - tran_tbl
            if(P_ck(j,i) > Pmin) P_all(j,ii) = tran_tbl*P_ck(j,i) + tran_atm*P_all(j,ii)
         enddo
      enddo

      !add PHX pressures to the combined C&K/T-tau table, T offset is -2
      do i=1,nT_phx
         ii = i
         if(T_phx(i) /= T_all(ii)) stop 'T_phx != T_all'
         do j=1,ng
            tran_tbl = transparency(i,zero_phx(j))
            tran_atm = 1d0-tran_tbl
            if(P_phx(j,i) > Pmin) P_all(j,ii) = tran_tbl*P_phx(j,i) + tran_atm*P_all(j,ii)
         enddo
      enddo

      !set validity range of combined table
      do j=1,ng
         i=1
         do while(P_all(j,i) > Pfill .and. i < nT_all)
            i=i+1
         enddo
         ibound(j)=i
      enddo
      
      ! fill and smooth so bicubic splines will be happy
      P_all(1,nT_all) = Pfill
      Pinterp_max = 0
      i_max = 0
      j_max = 0
      do i = 1, ng
         do j = max(1,ibound(i)-1), nT_all
            d0 = sqrt(dble((j-nT_all)**2 + (i-1)**2))
            jj = j
            do ii = i, ng
               if (ibound(ii)-1 > jj) then
                  d1 = sqrt(dble((j-jj)**2 + (i-ii)**2))
                  Pinterp = Pfill*d1/(d0+d1) + P_all(ii,jj)*d0/(d0+d1)
                  if (Pinterp > Pinterp_max .and. i == 1) then
                     Pinterp_max = Pinterp
                     i_max = i
                     j_max = j
                  end if
                  P_all(i,j) = Pinterp
                  exit
               end if
               jj = max(1,jj-1)
            end do
         end do
      end do
   
      ! smooth
      do i=1,20
         do ii=1,ng
            if (ibound(ii) == nT_all) cycle
            do jj=ibound(ii)-1,nT_all
               P_all(ii,jj) = ( &
                  P_all(ii,max(1,jj-1)) + &
                  P_all(ii,min(nT_all,jj+1)) + &
                  P_all(ii,jj) + &
                  P_all(max(1,ii-1),jj) + &
                  P_all(min(ng,ii+1),jj)) / 5
            end do
         end do
      end do
         
      
      

      !write out combined table
      open(4,file=f_output)
      write(4,'("#Table Version",i4)') table_version
      write(4,'(a40,a15,20i4)') head_phx,'| VALID RANGE: ', ibound
      write(4,'("#Teff(K)| Pgas@",13("  log g =",f5.2," "))') (-0.5+0.5*(i-1),i=1,ng)
      do i=1,nT_all
         write(4,'(1p,20e15.7)') T_all(i), P_all(:,i)
      enddo
      close(4)

      contains
      
      real(dp) function transparency(i,zero)
      integer :: i,zero         !transparency tells what fraction of the data point
      if(i>=zero) then          !in question to apply to the final table.  if the 
         transparency=0d0       !point is near the edge of the table, it is more
      else if(i == zero-1) then !"transparent" i.e. receives less weight
         transparency=0.25d0
      else if(i == zero-2) then
         transparency=0.5d0
      else if(i == zero-3) then
         transparency=0.75d0
      else
         transparency=1d0
      endif
      end function transparency
      
      end program table_merge
