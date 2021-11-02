! ***********************************************************************
!
!   Copyright (C) 2020  Aaron Dotter
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

!takes a CMS P,T table and converts to a Rho,T table
!also calculates extra quantities wanted by MESA EOS

program cms_resample
   
   use const_def
   use const_lib
   use interp_1d_def
   use interp_1d_lib
   use math_lib
   
   implicit none
   
   integer, parameter :: version = 1
   integer, parameter :: NT = 121
   integer, parameter :: NP = 441
   integer, parameter :: NRho = 281 ! -8 <= logRho <= +6 by 0.05
   integer :: ierr, i, j, io
   
   character(len=256) :: input, output

   real(dp) :: H_mass_fraction, He_mass_fraction
   
   !old EOS table
   real(dp), dimension(NP,NT) :: logT(NP,NT)
   real(dp), dimension(NP,NT) :: logRho, logP, logU, logS, dlnRho_dlnT_P, &
      dlnRho_dlnP_T, dlnS_dlnT_P, dlnS_dlnP_T, grad_ad

   !new EOS quantities
   real(dp), parameter :: logRho_min_for_interp = -8.1_dp
   real(dp), parameter :: logRho_min = -8.0_dp
   real(dp), parameter :: logRho_max =  6.0_dp
   real(dp), parameter :: delta_logRho = 0.05_dp
   real(dp), parameter :: delta_logT   = 0.05_dp
   real(dp) :: new_logT, new_logRho, new_logP, new_logU, new_logS, new_dlnS_dlnT, new_dlnS_dlnP
   real(dp) :: new_dlnRho_dlnT, new_dlnRho_dlnP, new_grad_ad, dU_dP_T, dU_dT_P, dS_dRho, mu, dP_dT
   real(dp) :: dS_dT, dU_dRho, Cv, Cp, gamma1, gamma3, eta, P, U, S, T, Rho, chiRho, chiT, lnfree_e
   real(dp) :: dS_dP_T, dS_dT_P, dse, dsp, dpe !for consistency check

   ierr=0
   
   if(command_argument_count()<2)then
      write(*,*) './cms_resample [input] [output]'
      stop
   endif
   
   call get_command_argument(1,input)
   call get_command_argument(2,output)

   call mesa_init

   open(newunit=io,file=trim(input),action='read',status='old',iostat=ierr)
   read(io,'(2x,f5.3)') H_mass_fraction
   read(io,*) !header
   do i=1,NT
      do j=1,NP
         read(io,'(1p99e15.6)') logT(j,i), logP(j,i), logRho(j,i), logU(j,i), logS(j,i), &
            dlnRho_dlnT_P(j,i), dlnRho_dlnP_T(j,i), dlnS_dlnT_P(j,i), dlnS_dlnP_T(j,i), &
            grad_ad(j,i)
      enddo
   enddo
   close(io)

   !set composition variables
   He_mass_fraction = 1.0_dp - H_mass_fraction

   write(*,*) 'io=', io
   open(newunit=io,file=trim(output),action='write',status='unknown',iostat=ierr)
   if(ierr/=0) stop
   write(io,'(99a14)') 'version', 'X', 'Z', 'num logTs', 'logT min', 'logT max', 'del logT', &
      'num logRhos', 'logRho min', 'logRho max', 'del logRho'
   write(io,'(i14,2f14.4,2(i14,3f14.4))') version, H_mass_fraction, 0.0_dp, NT, logT(1,1), &
      logT(1,NT), delta_logT, NRho, logRho_min, logRho_max, delta_logRho
   write(io,'(A)')

   write(io,'(99a15)') 'logT', 'logRho', 'logPgas', 'logU', 'logS', 'chiRho', 'chiT', 'Cp', 'Cv', &
      'dE_dRho', 'dS_dT', 'dS_dRho', 'mu', 'lnfree_e', 'gamma1', 'gamma3', 'grad_ad', 'eta', &
      'dsp', 'dse'
   !loop through logT and logRho points and interpolate
   do i=1,NT
      new_logT = logT(1,i)
      do j=1,NRho
         new_logRho = logRho_min + real(j-1,kind=dp)*delta_logRho
         
         call do_stuff(i, new_logRho, new_logP, new_logU, new_logS, &
            new_dlnRho_dlnT, new_dlnRho_dlnP, new_dlnS_dlnT, new_dlnS_dlnP, new_grad_ad)  

         P = exp10(new_logP)
         U = exp10(new_logU)
         S = exp10(new_logS)
         T = exp10(new_logT)
         Rho = exp10(new_logRho)

         !CMS eqn 5
         chiRho = 1.0_dp / new_dlnRho_dlnP ! dlnP/dlnRho at const T
         chiT   = -new_dlnRho_dlnT / new_dlnRho_dlnP !dlnP/dlnT at const Rho
         Cp = S * new_dlnS_dlnT ! at const P
         Cv = Cp - (P*chiT*chiT)/(Rho*T*chiRho) ! at const Rho (V)

         dU_dP_T = new_dlnRho_dlnP/Rho + T*S*new_dlnS_dlnP/P
         dU_dT_P = (P/(rho*T))*new_dlnRho_dlnT + S*new_dlnS_dlnT

         dS_dP_T = (S/P) * new_dlnS_dlnP
         dS_dT_P = (S/T) * new_dlnS_dlnT
         
         dU_dRho = dU_dP_T /( (rho/P) * new_dlnRho_dlnP)
         
         mu = 4.0_dp / (6.0_dp*H_mass_fraction + He_mass_fraction + 2.0_dp)
         lnfree_e = log(0.5_dp*(1.0_dp + H_mass_fraction))

         dS_dT = (S/T)*new_dlnS_dlnP * chiT !dS/dT_|Rho
         dS_dRho = new_dlnS_dlnT / new_dlnRho_dlnT * (S/Rho)

         !gamma3 = (P*chiT)/(rho*T*Cv) + 1.0_dp
         gamma3 = 1.0_dp - new_dlnS_dlnP / (new_dlnS_dlnT * new_dlnRho_dlnP)
         !gamma3 = new_dlnRho_dlnT
         gamma1 = chiT*(gamma3 - 1.0_dp) + chiRho

         eta = 0.0_dp

         dP_dT = P*chiT/T

         !dpe
         dse = T*(dS_dT/Cv) - 1.0_dp
         dsp = -rho*rho*(dS_dRho/dP_dT) - 1.0_dp
         dpe = 0.0_dp !not avaiable as yet...

         write(io,'(1p99e15.6)') new_logT, new_logRho, new_logP, new_logU, new_logS, &
            chiRho, chiT, Cp, Cv, dU_dRho, dS_dT, dS_dRho, mu, lnfree_e, gamma1, gamma3, &
            new_grad_ad, eta, dsp, dse
      enddo
   enddo
   close(io)
   
contains

   subroutine mesa_init
      call const_init(' ',ierr)     
      if (ierr /= 0) then
         write(0,*) 'const_init failed'
         stop 1
      end if
      
      call math_init()
      
   end subroutine mesa_init

   subroutine do_stuff(iT, new_logRho, new_logP, new_logU, new_logS, new_dlnRho_dlnT, &
            new_dlnRho_dlnP, new_dlnS_dlnT, new_dlnS_dlnP, new_grad_ad)
      integer, intent(in) :: iT
      real(dp), intent(in) :: new_logRho
      real(dp), intent(out) :: new_logP, new_logU, new_logS, new_dlnRho_dlnT, new_dlnRho_dlnP
      real(dp), intent(out) :: new_dlnS_dlnT, new_dlnS_dlnP, new_grad_ad
      real(dp), allocatable :: x_old(:)
      real(dp), allocatable :: y_old(:)
      real(dp) :: tmp_logRho(NP), tmp_logP(NP), tmp_logU(NP), tmp_logS(NP)
      real(dp) :: tmp_dlnRho_dlnT(NP), tmp_dlnRho_dlnP(NP), tmp_grad_ad(NP)
      real(dp) :: tmp_dlnS_dlnT(NP), tmp_dlnS_dlnP(NP)
      real(dp) :: x_new(1), y_new(1)
      integer :: j, count, num_pts, ierr
      integer, parameter :: nwork = NP*mp_work_size
      real(dp), target :: work_ary(nwork)
      real(dp), pointer :: work(:)

      work => work_ary
      
      count = 0
      do j=1,NP
         if(logRho(j,iT) >= logRho_min_for_interp)then
            if(logRho(j,iT) /= -8.0_dp)then
               count = count + 1
               tmp_logRho(count) = logRho(j,iT)
               tmp_logP(count) = logP(j,iT)
               tmp_logU(count) = logU(j,iT)
               tmp_logS(count) = logS(j,iT)
               tmp_dlnRho_dlnT(count) = dlnRho_dlnT_P(j,iT)
               tmp_dlnRho_dlnP(count) = dlnRho_dlnP_T(j,iT)
               tmp_dlnS_dlnT(count) = dlnS_dlnT_P(j,iT)
               tmp_dlnS_dlnP(count) = dlnS_dlnP_T(j,iT)
               tmp_grad_ad(count) = grad_ad(j,iT)
            endif
         endif
      enddo
      
      num_pts = count
      allocate(x_old(num_pts), y_old(num_pts))
      
      !get logP for logRho input
      x_old(1:num_pts) = tmp_logRho(1:num_pts)
      y_old(1:num_pts) = tmp_logP(1:num_pts)
      x_new(1) = new_logRho
      call interpolate_vector(num_pts, x_old, 1, x_new, y_old, y_new, &
         interp_m3a, nwork, work, 'sigh', ierr)
      new_logP = y_new(1)

      !now do other interpolations to get other EOS quantities
      !logU
      x_old(1:num_pts) = tmp_logP(1:num_pts)
      y_old(1:num_pts) = tmp_logU(1:num_pts)
      x_new(1) = new_logP
      call interpolate_vector(num_pts, x_old, 1, x_new, y_old, y_new, &
         interp_m3a, nwork, work, 'sigh', ierr)
      new_logU = y_new(1)

      !logS
      x_old(1:num_pts) = tmp_logP(1:num_pts)
      y_old(1:num_pts) = tmp_logS(1:num_pts)
      x_new(1) = new_logP
      call interpolate_vector(num_pts, x_old, 1, x_new, y_old, y_new, &
         interp_m3a, nwork, work, 'sigh', ierr)
      new_logS = y_new(1)      

      !dlnRho_dlnT_constP
      x_old(1:num_pts) = tmp_logP(1:num_pts)
      y_old(1:num_pts) = tmp_dlnRho_dlnT(1:num_pts)
      x_new(1) = new_logP
      call interpolate_vector(num_pts, x_old, 1, x_new, y_old, y_new, &
         interp_m3a, nwork, work, 'sigh', ierr)
      new_dlnRho_dlnT = y_new(1)      

      !dlnRho_dlnP_constT
      x_old(1:num_pts) = tmp_logP(1:num_pts)
      y_old(1:num_pts) = tmp_dlnRho_dlnP(1:num_pts)
      x_new(1) = new_logP
      call interpolate_vector(num_pts, x_old, 1, x_new, y_old, y_new, &
         interp_m3a, nwork, work, 'sigh', ierr)
      new_dlnRho_dlnP = y_new(1)      
      
      !dlnS_dlnT_constP
      x_old(1:num_pts) = tmp_logP(1:num_pts)
      y_old(1:num_pts) = tmp_dlnS_dlnT(1:num_pts)
      x_new(1) = new_logP
      call interpolate_vector(num_pts, x_old, 1, x_new, y_old, y_new, &
         interp_m3a, nwork, work, 'sigh', ierr)
      new_dlnS_dlnT = y_new(1)      

      !dlnS_dlnP_constT
      x_old(1:num_pts) = tmp_logP(1:num_pts)
      y_old(1:num_pts) = tmp_dlnS_dlnP(1:num_pts)
      x_new(1) = new_logP
      call interpolate_vector(num_pts, x_old, 1, x_new, y_old, y_new, &
         interp_m3a, nwork, work, 'sigh', ierr)
      new_dlnS_dlnP = y_new(1)      

      !grad_ad
      x_old(1:num_pts) = tmp_logP(1:num_pts)
      y_old(1:num_pts) = tmp_grad_ad(1:num_pts)
      x_new(1) = new_logP
      call interpolate_vector(num_pts, x_old, 1, x_new, y_old, y_new, &
         interp_m3a, nwork, work, 'sigh', ierr)
      new_grad_ad = min(0.5_dp, max(0.1_dp, y_new(1)))
      
   end subroutine do_stuff
   
end program cms_resample
