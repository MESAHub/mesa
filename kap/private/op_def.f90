! ***********************************************************************
!
!   Copyright (C) 2013  Haili Hu and Bill Paxton
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************
 

MODULE op_def
   IMPLICIT NONE

   integer, parameter :: nptot = 10000
   integer, parameter :: nrad = 17
   integer, parameter :: ipe = 17      

   integer,dimension(140:320),parameter :: JS=(/14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 15, 18, 19, 22, 23,&
   26, 27, 30, 31, 34, 34, 34, 34, 36, 36, 36, 36, 36, 36, 38, 38,&
   38, 38, 38, 39, 40, 40, 40, 40, 40, 41, 42, 42, 42, 42, 42, 43,&
   44, 44, 44, 44, 44, 45, 46, 46, 46, 46, 46, 47, 48, 48, 48, 48,&
   48, 49, 50, 50, 50, 50, 52, 52, 52, 52, 52, 52, 54, 54, 54, 54,&
   54, 54, 56, 56, 56, 56, 56, 56, 58, 58, 58, 58, 58, 58, 60, 60,&
   60, 60, 60, 60, 62, 62, 62, 62, 62, 63, 64, 64, 64, 64, 64, 65,&
   66, 66, 66, 66, 66, 67, 68, 68, 68, 68, 68, 69, 70, 70, 70, 70,&
   70, 71, 72, 72, 72, 72, 72, 73, 74, 74, 74, 74, 76, 76, 76, 76,&
   76, 76, 78, 78, 78, 78, 78, 78, 80, 80, 80, 80, 80, 80, 82, 82,&
   82, 82, 82, 82, 84, 84, 84, 84, 84, 84, 86, 86, 86, 86, 86, 87,&
   88, 88, 88, 88, 88/)

   integer,dimension(140:320),parameter :: JE=(/&
   52, 56, 56, 58, 58, 60, 60, 60, 60, 62, 62, 64, 64, 66, 66, 68,&
   68, 70, 70, 72, 72, 74, 74, 74, 74, 76, 76, 78, 78, 80, 80, 80,&
   80, 81, 82, 82, 82, 82, 82, 83, 84, 84, 84, 84, 84, 86, 86, 86,&
   86, 86, 86, 88, 88, 88, 88, 89, 90, 90, 90, 88, 88, 88, 88, 92,&
   92, 92, 92, 94, 94, 94, 94, 94, 94, 89, 90, 90, 90, 90, 90, 91,&
   92, 92, 92, 92, 92, 94, 94, 90, 90, 92, 92, 92, 92, 92, 92, 93,&
   94, 94, 94, 94, 94, 94, 94, 96, 96, 96, 96, 96, 96, 98, 98, 98,&
   98, 98, 98, 99,100,100,100,100,100,100,100,102,102,102,102,102,&
   102,103,104,104,104,104,104,104,104,106,106,106,106,106,106,107,&
   108,108,108,108,108,108,108,110,110,110,110,110,110,111,112,112,&
   112,112,112,112,112,114,114,114,114,114,114,115,116,116,116,116,&
   116,116,116,118,118/)

!        
   INTEGER,DIMENSION(17),parameter :: kz=(/1,2,6,7,8,10,11,12,13,14,16,18,20,24,25,26,28/)
   character(len=2), dimension(17),parameter :: name=(/'H ','He','C ','N ','O ','Ne','Na',&
                                       'Mg','Al','Si','S ','Ar','Ca','Cr','Mn',&
                                       'Fe','Ni'/)
   REAL,DIMENSION(17),PARAMETER :: AMASS=(/1.0080,4.0026,12.0111,14.0067,15.9994,20.179,&
                                 22.9898,24.305,26.9815,28.086,32.06,39.948,  &
                                 40.08,51.996,54.9380,55.847,58.71/)
   
   integer,save :: ite1, ite2, ite3, jne3 , ntotp, nc, nf
   integer,dimension(91),save :: jn1, jn2 
   integer,dimension(17),save :: int    
   real,save :: umin, umax
   real,dimension(17,91,25),save :: epatom, oplnck
   integer,dimension(17,91,25),save :: ne1p, ne2p,np,kp1,kp2,kp3,npp
   real,dimension(-1:28,28,91,25),save :: fionp
   real,allocatable,DIMENSION(:),save :: yy2,yx
   INTEGER,allocatable,DIMENSION(:),save :: nx

   integer, parameter :: op_cache_version = 1
   
END module op_def
