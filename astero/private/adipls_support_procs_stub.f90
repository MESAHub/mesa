! ***********************************************************************
!
!   Copyright (C) 2021  The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      
      

      ! routines that are called by adipls   
      ! uses adipls_support, so must compile after that   
      
      subroutine spcout_adi(x,y,aa,data,nn,iy,iaa,ispcpr)
         ! must set ispcpr > 0 to get this called
         use astero_def, only: store_new_oscillation_results, &
            el, order, em, cyclic_freq, inertia, num_results
         use adipls_support, only: adipls_mode_info
         use const_def, only: dp, pi4
         use utils_lib, only: mesa_error

         implicit none
         
         integer :: nn, iy, iaa, ispcpr
         real(dp) :: x(1:nn), y(1:iy,1:nn), aa(1:iaa,1:nn), data(8)

         

      end subroutine spcout_adi      
      
      
      subroutine modmod(x,aa,data,nn,ivarmd,iaa,imdmod)
         use const_def, only: dp
         integer :: nn, ivarmd, iaa, imdmod
         real(dp) :: x(nn), aa(iaa,nn), data(8)
      end subroutine modmod
      
      
      subroutine resdif
         return
      end subroutine resdif
