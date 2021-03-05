! ***********************************************************************
!
!   Copyright (C) 2010-2019  Aaron Dotter, Bill Paxton & The MESA Team
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

module mod_tau100

  use atm_def
  use atm_lib
  use chem_def
  use chem_lib, only: chem_init, basic_composition_info
  use const_def
  use const_lib
  use eos_def
  use eos_lib
  use kap_lib
  use utils_lib
  use math_lib
  use interp_2d_lib_db

  implicit none

  logical, parameter :: write_plot_files = .true.

  real(dp) :: tau_base ! tau at base of atm (1, 10, or 100)
  integer :: which_cond_layer, which_ck_layer
  character(len=256) :: my_mesa_dir, output_file1, output_file2  

  integer, parameter :: &
       cond_layer_tau1m2 = 26, cond_layer_tau1m1 = 32, cond_layer_tau1 = 38, &
       cond_layer_tau10 = 44, cond_layer_tau100 = 50
  integer, parameter :: &
       ck_layer_tau1m2 = 40, ck_layer_tau1m1 = 48, ck_layer_tau1 = 56, &
       ck_layer_tau10 = 64, ck_layer_tau100 = 72


  integer :: eos_handle, kap_handle

contains

  subroutine build_tau_tables

    integer, parameter :: num_isos = 7, nmet = 5
    integer, pointer :: chem_id(:), net_iso(:)
    real(dp), pointer :: xa(:)

    integer :: i, j, k, ierr, i_Teff, i_logg, logg_i_lo, logg_i_hi, &
         io_out1, io_out2, io, ii, jj
    include 'tau100_Ts.dek'
    character(len=256) :: filename     
    integer, parameter :: ng = 13, nT = num_tau100_Ts, num_layers = 50
    real(dp) :: loggs(ng), Teffs(nT), vals(num_layers), jnk1, jnk2
    real(dp) :: T(ng,nT), Pgas(ng,nT), logT(ng,nT), logPgas(ng,nT)
    real(dp) :: logT_plot(nT,ng), logPgas_plot(nT,ng)
    real(dp) :: X, Y, Z, XC, XN, XO, XNe, XMg, abar, zbar, z2bar, z53bar, logg, Teff
    integer, parameter :: nt_for_CK = 76, max_ng_for_CK = 11
    integer :: ng_for_CK_Teff(nt_for_CK)

    include 'formats'

    write(*,1) 'tau_base', tau_base     

    ng_for_CK_Teff(1:11) = 11
    ng_for_CK_Teff(12:17)= 10
    ng_for_CK_Teff(18:20)=  9
    ng_for_CK_Teff(21:23)=  8
    ng_for_CK_Teff(24:34)=  7
    ng_for_CK_Teff(35:39)=  6
    ng_for_CK_Teff(40)   =  7
    ng_for_CK_Teff(41:45)=  6
    ng_for_CK_Teff(46:52)=  5
    ng_for_CK_Teff(53:57)=  4
    ng_for_CK_Teff(58:65)=  3
    ng_for_CK_Teff(66:75)=  2
    ng_for_CK_Teff(76)   =  1

    ierr = 0
    io = 33 ! for reading

    ! setup for T(tau) calculations
    call set_table_composition

    loggs(1:ng) = (/ &
         0d0, 0.5d0, 1d0, 1.5d0, 2d0, 2.5d0, 3d0, 3.5d0, 4d0, 4.5d0, 5d0, 5.5d0, 6d0 /)
    Teffs(1:nT) = tau100_Ts(1:nT)

    Pgas = -1
    T = -1

    call read_COND
    call read_CK

    ! set logg = 0.0 to 2.0 for Teff = 100 to 3400 using T(tau)
    write(*,*) 'do T(tau)'
    do i = 1, ng
       logg = loggs(i)
       if (logg < 0d0) cycle
       if (logg > 2.01d0) exit
       do j = 1, nT
          Teff = Teffs(j)
          if (Teff < 99) cycle
          if (Teff > 3401) exit
          !write(*,*) 'logg Teff', logg, Teff
          call eval1_Ttau(logg, Teff, Pgas(i,j), T(i,j), ierr)
          if (ierr /= 0) then
             write(*,*) 'failed in eval1_Ttau', logg, Teff, i, j
             call mesa_error(__FILE__,__LINE__)
          end if
       end do
       !do j = 1, 4
       !   write(*,*)
       !end do
    end do

    write(*,*) 'smooth transitions'

    ! extrapolate logg = 5.5 and 6.0 for Teff = 3500 to 50000 using logg = 4.50 and 5.0
    do i = 1, ng
       logg = loggs(i)
       if (logg < 5.49d0) cycle
       if (logg > 5.51d0) exit
       do j = 1, nT
          Teff = Teffs(j)
          if (Teff < 3499) cycle
          if (Teff > 50001) exit
          Pgas(i,j) = (3*Pgas(i-1,j) - Pgas(i-2,j))/2
          T(i,j) = (3*T(i-1,j) - T(i-2,j))/2
          Pgas(i+1,j) = Pgas(i,j)
          T(i+1,j) = T(i,j)
       end do
    end do

    ! smooth the T(tau) region
    call smooth(0d0, 2.01d0, 99d0, 3401d0, 10)

    call fill_and_smooth_CK_corner()

    ! set Teff = 3400 and 3500 for logg 0.0 to 6.0 to linear interpolation of 3300 and 3600
    do i = 1, ng
       logg = loggs(i)
       do j = 1, nT
          Teff = Teffs(j)
          if (Teff < 3399) cycle
          Pgas(i,j) = (2*Pgas(i,j-1) + Pgas(i,j+2))/3
          T(i,j) = (2*T(i,j-1) + T(i,j+2))/3
          Pgas(i,j+1) = (Pgas(i,j-1) + 2*Pgas(i,j+2))/3
          T(i,j+1) = (T(i,j-1) + 2*T(i,j+2))/3
          exit
       end do
    end do

    call write_output()

    if (write_plot_files) call write_plots()

  contains

    subroutine fill_and_smooth_CK_corner()
      integer :: ibound(ng), j, i, jj, ii, j1
      real(dp) :: d0, d1, Pinterp, Tinterp
      real(dp), parameter :: Pfill = 10
      real(dp), parameter :: Tfill = 1d5
      include 'formats'

      do i=1,ng
         j=1
         do while (Teffs(j) < 4000)
            j=j+1
         end do
         do while(Pgas(i,j) > 0 .and. j < nT)
            j=j+1
         enddo
         ibound(i)=j
         !write(*,3) 'ibound(i) = j', i, j, loggs(i), Teffs(j)
      enddo

      Pgas(1,nT) = Pfill
      T(1,nT) = Tfill
      do i = 1, ng
         j1 = max(1,ibound(i)-1)
         do j = j1, nT
            d0 = sqrt(dble((j-nT)**2 + (i-1)**2))
            jj = j
            do ii = i, ng
               if (ibound(ii)-1 > jj) then
                  d1 = sqrt(dble((j-jj)**2 + (i-ii)**2))
                  Pinterp = Pfill*d1/(d0+d1) + Pgas(ii,jj)*d0/(d0+d1)
                  Pgas(i,j) = Pinterp
                  Tinterp = Tfill*d1/(d0+d1) + T(ii,jj)*d0/(d0+d1)
                  T(i,j) = Tinterp
                  !if (j == j1) then
                  !   write(*,3) 'T(i,j)', i, j, T(i,j), loggs(i), Teffs(j)
                  !end if
                  exit
               end if
               jj = max(1,jj-1)
            end do
         end do
      end do
      do i=1,20
         do ii=1,ng
            if (ibound(ii) == nT) cycle
            do jj=ibound(ii)-1,nT
               Pgas(ii,jj) = ( &
                    Pgas(ii,max(1,jj-1)) + &
                    Pgas(ii,min(nT,jj+1)) + &
                    Pgas(ii,jj) + &
                    Pgas(max(1,ii-1),jj) + &
                    Pgas(min(ng,ii+1),jj)) / 5
               T(ii,jj) = ( &
                    T(ii,max(1,jj-1)) + &
                    T(ii,min(nT,jj+1)) + &
                    T(ii,jj) + &
                    T(max(1,ii-1),jj) + &
                    T(min(ng,ii+1),jj)) / 5
            end do
         end do
      end do
    end subroutine fill_and_smooth_CK_corner

    !****

    subroutine smooth(logg_lo, logg_hi, Teff_lo, Teff_hi, n)
      real(dp), intent(in) :: logg_lo, logg_hi, Teff_lo, Teff_hi
      integer, intent(in) :: n
      real(dp) :: logg, Teff
      integer :: k, ii, jj
      include 'formats'
      do k = 1, n
         do ii = 1, ng
            logg = loggs(ii)
            if (logg < logg_lo) cycle
            if (logg > logg_hi) exit
            do jj = 1, nT
               Teff = Teffs(jj)
               if (Teff < Teff_lo) cycle
               if (Teff > Teff_hi) exit
               Pgas(ii,jj) = ( &
                    Pgas(ii,max(1,jj-1)) + &
                    Pgas(ii,min(nT,jj+1)) + &
                    Pgas(ii,jj) + &
                    Pgas(max(1,ii-1),jj) + &
                    Pgas(min(ng,ii+1),jj)) / 5
               T(ii,jj) = ( &
                    T(ii,max(1,jj-1)) + &
                    T(ii,min(nT,jj+1)) + &
                    T(ii,jj) + &
                    T(max(1,ii-1),jj) + &
                    T(min(ng,ii+1),jj)) / 5
               !write(*,3) 'smooth', ii, jj, logg, Teff, Pgas(ii,jj), T(ii,jj)
            end do
         end do
      end do
    end subroutine smooth

    !****
    
    subroutine write_plots()
      
      integer :: i, j, ii, jj

      ! rearrange for plotting
      do i=1,ng
         ii = ng-i+1
         do j=1,nT
            jj = nT-j+1
            logT_plot(j,i) = safe_log10(T(i,j))
            logPgas_plot(j,i) = safe_log10(Pgas(i,j))
         end do
      end do

      open(io,file='plot_tau100/logT.data',action='write')
      write(io,'(e20.10)') logT_plot(:,:)
      close(io)

      open(io,file='plot_tau100/logPgas.data',action='write')
      write(io,'(e20.10)') logPgas_plot(:,:)
      close(io)

      open(io,file='plot_tau100/Teff.data',action='write')
      do j=1,nT
         write(io,'(e20.10)') Teffs(j)
      end do
      close(io)

      open(io,file='plot_tau100/logg.data',action='write')
      do i=1,ng
         write(io,'(e20.10)') loggs(i)
      end do
      close(io)

    end subroutine write_plots

    !****

    subroutine write_output
      write(*,*) trim(output_file1)
      io_out1 = 33
      open(io_out1, file=output_file1)         
      write(*,*) trim(output_file2)
      io_out2 = 34
      open(io_out2, file=output_file2)         
      write(io_out1,'(a15)',advance='no') "#Teff(K)| Pgas@" 
      write(io_out2,'(a15)',advance='no') "#Teff(K)|    T@" 
      do i = 1, ng
         write(io_out1,fmt='("  log g =",f5.2," ")',advance='no') loggs(i)
         write(io_out2,fmt='("  log g =",f5.2," ")',advance='no') loggs(i)
      end do
      write(io_out1,*)
      write(io_out2,*)
      do j = 1, nT
         write(io_out1,fmt='(e15.7)',advance='no') Teffs(j)
         write(io_out2,fmt='(e15.7)',advance='no') Teffs(j)
         do i = 1, ng
            write(io_out1,fmt='(e15.7)',advance='no') Pgas(i,j)
            write(io_out2,fmt='(e15.7)',advance='no') T(i,j)
         end do
         write(io_out1,*)
         write(io_out2,*)
      end do
      close(io_out1)
      close(io_out2)
    end subroutine write_output

    !****
    
    subroutine read_COND()
      ! read COND; store logg = 2.5 to 6.0, Teff = 100 to 3600
      write(*,*) 'read COND'
      call read_phx(2.49d0, 6.01d0, 99d0, 3601d0, &
           'atm_input_data/cond/lte', '-0.0.AMES-Cond-2000.20')
    end subroutine read_COND

    !****

    subroutine read_phx(logg_lo, logg_hi, Teff_lo, Teff_hi, prefix, suffix)
      real(dp), intent(in) :: logg_lo, logg_hi, Teff_lo, Teff_hi
      character (len=*), intent(in) :: prefix, suffix
      integer :: T2
      do i = 1, ng
         logg = loggs(i)
         if (logg < logg_lo) cycle
         if (logg > logg_hi) exit
         do j = 1, nT
            Teff = Teffs(j)
            if (Teff < Teff_lo) cycle
            if (Teff > Teff_hi) exit
            T2 = floor(tau100_Ts(j))/100
            if (T2 > 99) then
               write(filename,'(a,i3,a,f3.1,a)') &
                    trim(prefix), T2, '-', loggs(i), trim(suffix)
            else if (T2 > 9) then
               write(filename,'(a,i2,a,f3.1,a)') &
                    trim(prefix), T2, '-', loggs(i), trim(suffix)
            else
               write(filename,'(a,i1,a,f3.1,a)') &
                    trim(prefix) // '0', T2, '-', loggs(i), trim(suffix)
            end if
            !write(*,*) trim(filename)
            open(io,file=trim(filename),action='read',status='old',iostat=ierr)
            if (ierr /= 0) then
               write(*,*) 'failed to open ' // trim(filename)
               call mesa_error(__FILE__,__LINE__)
            end if
            call read_file(i,j)
            close(io)
         end do
         !write(*,*)
      end do
    end subroutine read_phx

    !****

    subroutine read_file(i,j)
      integer, intent(in) :: i, j
      integer :: layers
      real(dp) :: Pg, Pe
      read(io,*)
      read(io,'(i5)') layers
      if (layers /= num_layers) then
         write(*,*) 'layers /= n', i, j, layers
         call mesa_error(__FILE__,__LINE__)
      end if
      call read_vals(.false.) ! tau
      call read_vals(.false.) ! temperature
      T(i,j) = vals(num_layers)
      call read_vals(.true.) ! flxrad
      call read_vals(.true.) ! terad
      call read_vals(.true.) ! bkmean
      call read_vals(.true.) ! jkmean
      call read_vals(.true.) ! fkmean
      call read_vals(.true.) ! rkmean
      call read_vals(.true.) ! pgas
      Pg = vals(which_cond_layer)
      call read_vals(.true.) ! pe
      Pe = vals(which_cond_layer)
      Pgas(i,j) = Pg + Pe
    end subroutine read_file

    !****
    
    subroutine read_vals(skip)
      logical, intent(in) :: skip
      integer :: k
      if (skip) read(io,*) ! skip line
      do k=1,16
         read(io,*) vals(1+(k-1)*3:k*3)
      end do
      read(io,*) vals(49:50)
    end subroutine read_vals

    !****

    subroutine read_CK
      ! read C&K; Teff = 3500 to 50000
      integer :: i, j1, j, k, npi
      write(*,*) 'read C&K'
      filename = 'atm_input_data/ck03/ck03_structures/ap00k2odfnew.dat'
      open(io,file=trim(filename),action='read',status='old',iostat=ierr)
      if (ierr /= 0) then
         write(*,*) 'failed to open ' // trim(filename)
         call mesa_error(__FILE__,__LINE__)
      end if
      do i=1,nt_for_CK
         j1 = max_ng_for_CK - ng_for_CK_Teff(i) + 1
         do j = j1, max_ng_for_CK
            read(io,'(4x,f8.0,9x,f8.5)') Teff, logg
            call locate(Teff, logg, ii, jj)
            do k=1,21
               read(io,*)
            enddo
            read(io,'(10x,i3)') npi
            ! 72 plane parallel layers from logtau = -6.875 to logtau = 2.0 by 0.125
            ! tau = 100 in line 72
            ! tau = 10 in line 64
            ! tau = 1 in line 56
            ! tau = 0.1 in line 48
            ! tau = 0.01 in line 40
            if (npi /= 72) then
               write(*,*) 'unexpected value for npi', npi, i, j
               call mesa_error(__FILE__,__LINE__)
            end if
            do k=1,which_ck_layer-1
               read(io,*)
            enddo
            read(io,'(15x,f9.1,e10.3)') T(ii,jj), Pgas(ii,jj)
            do k=which_ck_layer+1,npi
               read(io,*)
            enddo
            read(io,*)
            read(io,*)
         enddo
      enddo
      close(io)
    end subroutine read_CK

    !****

    subroutine locate(Teff, logg, i, j)
      real(dp), intent(in) :: Teff, logg
      integer, intent(out) :: i, j
      real(dp), parameter :: eps = 1d-6
      do i = 1, ng
         if (loggs(i) < logg-eps) cycle
         do j = 1, nT
            if (Teffs(j) < Teff-eps) cycle
            return
         end do
      end do
      write(*,*) 'failed in locate'
      call mesa_error(__FILE__,__LINE__)
    end subroutine locate

    !****

    subroutine eval1_Ttau(logg, Teff, Pgas, T, ierr)
      real(dp), intent(in) :: logg, Teff
      real(dp), intent(out) :: Pgas, T
      integer, intent(out) :: ierr

      integer :: iters
      real(dp) :: M, R, L, Teff_out, kap, err, Pextra_factor
      real(dp) :: lnT, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap
      real(dp) :: lnP, dlnP_dL, dlnP_dlnR, dlnP_dlnm, dlnP_dlnkap

      logical, parameter :: skip_partials = .true.
      integer, parameter :: max_iters = 100
      real(dp), parameter :: errtol = 1d-6
      M = Msun
      R = sqrt ( standard_cgrav*M / 10d0**logg )               
      L = pi*crad*clight * R**2 * Teff**4
      ierr = 0
      Pextra_factor = 1
      call atm_eval_T_tau_uniform( &
           tau_base, L, R, M, standard_cgrav, 0.2d0*(1 + X), Pextra_factor, &
           ATM_T_TAU_EDDINGTON, eos_proc, kap_proc, errtol, max_iters, skip_partials, &
           Teff_out, kap, &
           lnT, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
           lnP, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
           ierr)    
      if (ierr /= 0 .or. abs(Teff-Teff_out) > 1d-2) then
         Pgas = -1
         T = -1
         write(*,*) 'failed in atm_get_grey_and_kap', logg, Teff
      else
         T = exp(lnT)
         Pgas = max(1d-99, exp(lnP) - Radiation_Pressure(T))
      end if
    end subroutine eval1_Ttau

    !****

    subroutine set_table_composition()
      integer :: i
      real(dp) :: dabar_dx(num_isos), dzbar_dx(num_isos), ye, sumx, &
           xh, xhe, xz, approx_abar, approx_zbar, mass_correction
      real(dp), parameter :: &
           Xbbn   = 0.75d0, &
           Ybbn   = 0.25d0, &           
           Xsun   = 7.348d-1, &
           Ysun   = 2.348d-1, &
           Zsun   = 1.685d-2

      Z = Zsun
      Y = 1.5d0*Z + Ybbn
      X = 1d0 - (Y+Z)

      XC  = 0.1721d0*Z
      XN  = 0.0504d0*Z
      XO  = 0.4680d0*Z
      XNe = 0.1050d0*Z
      XMg = 0.2450d0*Z

      allocate(xa(num_isos), chem_id(num_isos), net_iso(num_chem_isos))

      chem_id(:) = (/ ih1, ihe4, ic12, in14, io16, ine20, img24 /)
      net_iso(:) = 0
      do i=1,num_isos
         net_iso(chem_id(i)) = i
      end do

      xa(:) = (/ X, Y, xc, xn, xo, xne, xmg /)
      xa(num_isos) = 1 - sum(xa(:))
      call basic_composition_info( &
           num_isos, chem_id, xa, xh, xhe, xz, abar, zbar, z2bar, z53bar, ye, &
           mass_correction, sumx)

    end subroutine set_table_composition

    !****
    
    subroutine eos_proc( &
         lnP, lnT, &
         lnRho, res, dres_dlnRho, dres_dlnT, &
         ierr)

      use eos_def, only: num_eos_basic_results
      use eos_lib, only: eosPT_get, radiation_pressure

      real(dp), intent(in)  :: lnP
      real(dp), intent(in)  :: lnT
      real(dp), intent(out) :: lnRho
      real(dp), intent(out) :: res(:)
      real(dp), intent(out) :: dres_dlnRho(:)
      real(dp), intent(out) :: dres_dlnT(:)
      integer, intent(out)  :: ierr

      real(dp) :: T, P, Prad, Pgas, logPgas, rho
      real(dp) :: logRho, dlnRho_dlnPgas, dlnRho_dlnT
      real(dp), dimension(num_eos_basic_results) :: &
           dres_dabar, dres_dzbar

      T = exp(lnT)
      P = exp(lnP)

      Prad = radiation_pressure(T)
      Pgas = max(1E-99_dp, P - Prad)
      logPgas = log10(Pgas)

      call eosPT_get( &
           eos_handle, Z, X, abar, zbar, &
           num_isos, chem_id, net_iso, xa, &
           Pgas, logPgas, T, lnT/ln10, &
           Rho, logRho, dlnRho_dlnPgas, dlnRho_dlnT, &
           res, dres_dlnRho, dres_dlnT, dres_dabar, dres_dzbar, ierr)

      lnRho = logRho*ln10

    end subroutine eos_proc

    !****

    subroutine kap_proc( &
         lnRho, lnT, res, dres_dlnRho, dres_dlnT, &
         kap, dlnkap_dlnRho, dlnkap_dlnT, &
         ierr)

      use kap_def, only: num_kap_fracs
      use kap_lib, only: kap_get
      use eos_def, only: i_lnfree_e, i_eta

      real(dp), intent(in)  :: lnRho
      real(dp), intent(in)  :: lnT
      real(dp), intent(in)  :: res(:)
      real(dp), intent(in)  :: dres_dlnRho(:)
      real(dp), intent(in)  :: dres_dlnT(:)
      real(dp), intent(out) :: kap
      real(dp), intent(out) :: dlnkap_dlnRho
      real(dp), intent(out) :: dlnkap_dlnT
      integer, intent(out)  :: ierr

      real(dp) :: kap_fracs(num_kap_fracs), dlnkap_dxa(num_isos)

      call kap_get( &
           kap_handle, num_isos, chem_id, net_iso, xa, &
           lnRho/ln10, lnT/ln10, res(i_lnfree_e), dres_dlnRho(i_lnfree_e), dres_dlnT(i_lnfree_e), &
           res(i_eta), dres_dlnRho(i_eta), dres_dlnT(i_eta), &
           kap_fracs, kap, dlnkap_dlnRho, dlnkap_dlnT, dlnkap_dxa, ierr)

    end subroutine kap_proc

  end subroutine build_tau_tables

  !****

  subroutine mesa_init(ierr)
    integer, intent(out) :: ierr
    logical, parameter :: use_cache = .true.
    call const_init(my_mesa_dir, ierr)      
    if (ierr /= 0) call mesa_error(__FILE__,__LINE__)         
    call math_init()
    call chem_init('isotopes.data', ierr)
    if (ierr /= 0) call mesa_error(__FILE__,__LINE__)         
    call eos_init(' ', use_cache, ierr)
    if (ierr /= 0) call mesa_error(__FILE__,__LINE__)            
    eos_handle = alloc_eos_handle(ierr)
    if (ierr /= 0) call mesa_error(__FILE__,__LINE__)            
    call kap_init(use_cache, ' ', ierr) 
    if (ierr /= 0) call mesa_error(__FILE__,__LINE__)            
    kap_handle = alloc_kap_handle(ierr)
    if (ierr /= 0) call mesa_error(__FILE__,__LINE__)            
    call atm_init(.false., ierr)
    if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
  end subroutine mesa_init


  subroutine mesa_shutdown
    call eos_shutdown
    call free_eos_handle(eos_handle)
    call kap_shutdown
    call free_kap_handle(kap_handle)
    call atm_shutdown
  end subroutine mesa_shutdown

end module mod_tau100


program create_tau100

  use mod_tau100

  integer :: ierr
  ierr = 0

  my_mesa_dir = '../..'
  call mesa_init(ierr)
  if (ierr /= 0) call mesa_error(__FILE__,__LINE__)

  tau_base = 0.1d0
  which_cond_layer = cond_layer_tau1m1
  which_ck_layer = ck_layer_tau1m1
  output_file1 = 'atm_data/tau1m1_Pgas.data'
  output_file2 = 'atm_data/tau1m1_T.data'
  call build_tau_tables

  tau_base = 1d0
  which_cond_layer = cond_layer_tau1
  which_ck_layer = ck_layer_tau1
  output_file1 = 'atm_data/tau1_Pgas.data'
  output_file2 = 'atm_data/tau1_T.data'
  call build_tau_tables

  tau_base = 10d0
  which_cond_layer = cond_layer_tau10
  which_ck_layer = ck_layer_tau10
  output_file1 = 'atm_data/tau10_Pgas.data'
  output_file2 = 'atm_data/tau10_T.data'
  call build_tau_tables

  tau_base = 100d0
  which_cond_layer = cond_layer_tau100
  which_ck_layer = ck_layer_tau100
  output_file1 = 'atm_data/tau100_Pgas.data'
  output_file2 = 'atm_data/tau100_T.data'
  call build_tau_tables

  call mesa_shutdown

end program create_tau100
