! ***********************************************************************
!
!   Copyright (C) 2009-2019  Aaron Dotter, Bill Paxton & The MESA Team
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
      program sample_kap
      use kap_lib
      use kap_def
      use chem_lib
      use chem_def
      use const_def
      use const_lib
      use math_lib
      use utils_lib, only: mesa_error
      
      implicit none
      !this program demonstrates how to use mesa/kap in a stellar structure code
      !it reads in a mesa/star model of an AGB star so that it makes full use
      !of mesa/kap's capabilities. modules kap_lib and kap_def contain access to
      !all of the pieces required to set up and use mesa/kap.

      character(len=256) :: model_file, output_file
      logical :: use_cache, show_info
      integer :: handle1, handle2, i, ii, iounit, Npts, Nspec
      type (Kap_General_Info), pointer :: rq1, rq2
      integer, parameter :: maxpts = 2000, maxspec = 31
      integer :: ierr
      
      integer, parameter :: h1 = 1
      integer, parameter :: h2 = 2
      integer, parameter :: he3 = 3
      integer, parameter :: he4 = 4
      integer, parameter :: li7 = 5
      integer, parameter :: be7 = 6
      integer, parameter :: b8 = 7
      integer, parameter :: c12 = 8
      integer, parameter :: c13 = 9
      integer, parameter :: n13 = 10
      integer, parameter :: n14 = 11
      integer, parameter :: n15 = 12
      integer, parameter :: o16 = 13
      integer, parameter :: o17 = 14
      integer, parameter :: o18 = 15
      integer, parameter :: f19 = 16
      integer, parameter :: ne20 = 17
      integer, parameter :: ne21 = 18
      integer, parameter :: ne22 = 19
      integer, parameter :: na22 = 20
      integer, parameter :: na23 = 21
      integer, parameter :: mg24 = 22
      integer, parameter :: mg25 = 23
      integer, parameter :: mg26 = 24
      integer, parameter :: al26 = 25
      integer, parameter :: al27 = 26
      integer, parameter :: si28 = 27
      integer, parameter :: si29 = 28
      integer, parameter :: si30 = 29
      integer, parameter :: p31 = 30
      integer, parameter :: s32 = 31
      
      real(dp) :: Mstar, Xc, Xn, Xo, Xne, xc_base, xn_base, xo_base, xne_base, &
         zbar, frac_Type2, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
      real(dp) :: Z_init, Z
      real(dp) :: lnRho(maxpts), lnT(maxpts), logRho(maxpts), logT(maxpts), X(maxspec,maxpts)
      real(dp) :: lnR(maxpts), L, dq(maxpts)
      real(dp) :: kappa(maxpts), dlnkap_dlnRho(maxpts), dlnkap_dlnT(maxpts)
      real(dp) :: kappaCO(maxpts), dlnkapCO_dlnRho(maxpts), dlnkapCO_dlnT(maxpts)
      character (len=32) :: my_mesa_dir
      
      use_cache = .false.
      show_info = .false.
      model_file = 'sample_kap_agb.model'
      output_file = 'kap_test.data'
      ierr = 0

      ! initialization and setup

      my_mesa_dir = '../..'         
      call const_init(my_mesa_dir,ierr)     
      if (ierr /= 0) then
         write(*,*) 'const_init failed'
         call mesa_error(__FILE__,__LINE__)
      end if        
      
      call math_init()

      call chem_init('isotopes.data', ierr)
      call kap_init(use_cache, '', ierr) 
      if(ierr/=0) stop 'problem in kap_init'

      !next it is necessary to create a 'handle' for the general kap structure
      !using handles, it is possible to simultaneously access more than one
      !working copy of the opacity subroutines with different settings
      handle1 = alloc_kap_handle_using_inlist('inlist_sample', ierr)
      call kap_ptr(handle1, rq1, ierr)
      rq1% use_Type2_opacities = .false.
      
      handle2 = alloc_kap_handle_using_inlist('inlist_sample', ierr)
      call kap_ptr(handle2, rq2, ierr)
      rq2% use_Type2_opacities = .true.

      if(ierr/=0) stop 'problem in alloc_kap_handle'

      !read in AGB model
      iounit=99
      open(unit=iounit,file=trim(model_file),status='old',iostat=ierr)
      if(ierr/=0) stop 'problem opening agb.mod file'
      read(iounit,*)
      read(iounit,*)            !skip 3 header lines
      read(iounit,*)
      read(iounit,1) Mstar      !read stellar mass
      read(iounit,1) Z_init     !read initial Z
      read(iounit,2) Npts       !read number of points in model
      read(iounit,*)            !skip
      read(iounit,2) Nspec      !read number of chemical species in model
      read(iounit,*)            !skip 2 lines
      read(iounit,*)
      
      write(*,*) ' Npts', Npts
      write(*,*) 'Nspec', Nspec
      
      do i=1,Npts               !read model
         read(iounit,*) ii, lnRho(i), lnT(i), lnR(i), L, dq(i), X(1:Nspec,i)
         if (ii /= i) then
            write(*,*) 'bad data for zone', i
            stop
         end if
      enddo
      close(iounit)

      write(*,*)
      write(*,*) 'Z_init', Z_init
      write(*,*)

      rq1% Zbase = Z_init
      rq2% Zbase = Z_init
      
      XC_base = GN93_element_zfrac(e_C)*Z_init
      XN_base = GN93_element_zfrac(e_N)*Z_init
      XO_base = GN93_element_zfrac(e_O)*Z_init
      XNe_base = GN93_element_zfrac(e_Ne)*Z_init

      logRho(:) = lnRho(:)/ln10 !convert ln's to log10's
      logT(:)   = lnT(:)  /ln10   
      zbar = 1.5d0              ! needed for electron conduction at high rho
      lnfree_e = 0d0            ! needed for Compton at high T
      d_lnfree_e_dlnRho = 0d0
      d_lnfree_e_dlnT = 0d0

!$omp parallel do private(i,ierr,XC,XN,XO,XNe,Z) schedule(dynamic)
      do i=1,Npts            
         XC = X(c12,i) + X(c13,i)
         XN = X(n14,i) + X(n15,i)
         XO = X(o16,i) + X(o17,i) + X(o18,i)
         XNe = X(ne20,i) + X(ne21,i) + X(ne22,i)
         Z = 1d0 - (X(h1,i)+X(he3,i)+X(he4,i))

         call kap_get( &
            handle1, zbar, X(h1,i), Z, XC, XN, XO, XNe, logRho(i), logT(i), &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            frac_Type2, kappa(i), dlnkap_dlnRho(i), dlnkap_dlnT(i), ierr)
         if(ierr/=0) write(*,*) 'kap_get (Type 1) failed at i=', i

         call kap_get( &
            handle2, zbar, X(h1,i), Z, XC, XN, XO, XNe, logRho(i), logT(i), &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            frac_Type2, kappaCO(i), dlnkapCO_dlnRho(i), dlnkapCO_dlnT(i), ierr)
         if(ierr/=0) write(*,*) 'kap_get (Type 2) failed at i=', i

      enddo
!$omp end parallel do

      open(unit=iounit,file=trim(output_file),iostat=ierr)
      if(ierr/=0) stop 'problem opening kap_test.data file'
      
      write(*,*) 'write ' // trim(output_file)

      write(iounit,3) 'grid', 'log_T', 'log_Rho', 'kappa', 'kappa_CO', &
      'dlnK_dlnRho', 'dlnK_dlnT'

      do i=1,Npts
         write(iounit,4) i, logT(i), logRho(i), kappa(i), kappaCO(i), &
         dlnkap_dlnRho(i), dlnkap_dlnT(i)
       enddo
       
       close(iounit)

      !all finished? then deallocate the handle and unload the opacity tables
      call free_kap_handle(handle1)
      call free_kap_handle(handle2)
      call kap_shutdown

 1    format(37x,e23.16)
 2    format(37x,i6)
 3    format(a28,99(a26,1x))
 4    format(i28,99(1pes26.16e3,1x))
      
      end program
