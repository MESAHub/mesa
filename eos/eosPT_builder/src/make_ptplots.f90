! ***********************************************************************
!
!   Copyright (C) 2009-2019  Bill Paxton & The MESA Team
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

      module make_PTplots
      
      use eos_def
      use eos_lib
      use const_def
      use chem_def
      use math_lib
      
      implicit none
      
      integer, parameter :: eosPT = 1
      integer, parameter :: scvh_eos = 2
      integer, parameter :: mesa_eos = 3
      
      contains
      
      
      subroutine Do_Test
         use eval_eosPT
         use utils_lib, only: is_bad
         
         integer :: which_eos, call_number, ierr
         double precision :: &
            X, Z, abar, zbar, logT, T, logPgas, Pgas, logW, logRho_guess, &
            logRho, logE, logS, chiRho, chiT, &
            Cp, Cv, dE_dRho, dS_dT, dS_dRho, &
            mu, log_free_e, gamma1, gamma3, grad_ad, eta, P, Prad, E, S
         character (len=256) :: dir
         integer, parameter :: io_unit0 = 40
         
         include 'formats'
         
         write(*,*) 'do1_test for eosPT'

         which_eos = mesa_eos
         ierr = 0
         
         X = 1.00d0  ! 0.00d0, 0.20d0, 0.40d0, 0.60d0, 0.80d0
         Z = 0.00d0  ! 0.00, 0.02, 0.04

         call Init_Composition(X, Z, abar, zbar)                  

         logT = 3.5d0

         !logPgas = 15
         !logW = logPgas - 4*logT
         
         logW = -14.5d0
         logPgas = logW + 4*logT
         
         logW =    -1.4800000000000001d+01
         logT =     3.5397959183673473d+00
               
         T = 10 ** logT
         Pgas = 10 ** logPgas

         call_number = 1
         logRho_guess = 0
         select case(which_eos)
            case (scvh_eos)
               call scvh_4_PT( &
                  logW, logT, T, logPgas, Pgas, &
                  abar, zbar, X, Z, logRho_guess, &
                  logRho, logE, logS, chiRho, chiT, &
                  Cp, Cv, dE_dRho, dS_dT, dS_dRho, &
                  mu, log_free_e, gamma1, gamma3, grad_ad, eta, call_number, ierr)
            case (mesa_eos)
               call mesa_4_PT( &
                  logW, logT, T, logPgas, Pgas, &
                  abar, zbar, X, Z, logRho_guess, &
                  logRho, logE, logS, chiRho, chiT, &
                  Cp, Cv, dE_dRho, dS_dT, dS_dRho, &
                  mu, log_free_e, gamma1, gamma3, grad_ad, eta, call_number, ierr)
            case (eosPT)
               call eos_4_Pgas_T( &
                  logW, logT, T, logPgas, Pgas, &
                  abar, zbar, X, Z, logRho_guess, &
                  logRho, logE, logS, chiRho, chiT, &
                  Cp, Cv, dE_dRho, dS_dT, dS_dRho, &
                  mu, log_free_e, gamma1, gamma3, grad_ad, eta, call_number, ierr)
            case default
               stop 'bad value for which_eos'
         end select
         
         write(*,1) 'logW', logW
         write(*,1) 'logT', logT
         write(*,1) 'logPgas', log10(Pgas)
         write(*,1) 'abar', abar
         write(*,1) 'zbar', zbar
         write(*,1) 'X', X
         write(*,1) 'Z', Z
         write(*,*)
         write(*,1) 'logRho', logRho
         write(*,1) 'logS', logS
         write(*,1) 'logE', logE
         write(*,1) 'chiRho', chiRho
         write(*,1) 'chiT', chiT
         write(*,1) 'Cp', Cp
         write(*,1) 'Cv', Cv
         write(*,1) 'dE_dRho', dE_dRho
         write(*,1) 'dS_dT', dS_dT
         write(*,1) 'dS_dRho', dS_dRho
         write(*,1) 'mu', mu
         write(*,1) 'log_free_e', log_free_e
         write(*,1) 'gamma1', gamma1
         write(*,1) 'gamma3', gamma3
         write(*,1) 'grad_ad', grad_ad
         write(*,1) 'eta', eta
         write(*,*)
         Prad = crad*T**4/3
         P = Pgas + Prad
         write(*,1) 'logP', log10(P)
         write(*,1) 'Pgas/P', Pgas/P
         write(*,*)
         stop 'Do_Test'
               
      
      end subroutine Do_Test

      
      subroutine make_plot_files
         use eval_eosPT
         use utils_lib, only: is_bad
         
         integer :: which_eos, logT_points, logW_points, io_params, io_logW, io_logT, &
            io_first, io_last, io, num_vals, j, i, k, call_number, ierr
         double precision, pointer :: output_values(:,:,:)
         double precision :: &
            X, Z, abar, zbar, logT_max, logT_min, logW_min, logW_max, dlogT, dlogW, &
            logT, T, logPgas, Pgas, logW, logRho_guess, &
            logRho, logE, logS, chiRho, chiT, &
            Cp, Cv, dE_dRho, dS_dT, dS_dRho, &
            mu, log_free_e, gamma1, gamma3, grad_ad, eta, energy, S, Rho
         character (len=256) :: dir
         integer, parameter :: io_unit0 = 40
         
         include 'formats'
         
         which_eos = scvh_eos
         
         X = 1d0  ! 0.00d0, 0.20d0, 0.40d0, 0.60d0, 0.80d0
         Z = 0d0  ! 0.00, 0.02, 0.04

         call Init_Composition(X, Z, abar, zbar)
                  
         logT_points = 100
         logW_points = 100
         
         logT_min = 2.9d0
         logT_max = 5.5d0

         logW_min = -3.0d0
         logW_max = -1.5d0


         io_params = io_unit0
         io_logW = io_unit0+1
         io_logT = io_unit0+2
         io_first = io_unit0+3

         dir = 'plot_data'
         call Open_Plot_Outfiles(io_first, io_last, io_params, io_logW, io_logT, dir)
         write(io_params, '(2(f10.6),2(i7))') Z, X, logW_points, logT_points
         close(io_params)
         num_vals  = io_last - io_first + 1
         allocate(output_values(logW_points,logT_points,num_vals))
         
         dlogT = (logT_max - logT_min)/(logT_points-1)
         dlogW = (logW_max - logW_min)/(logW_points-1)
         
         call_number = 0
      
         logRho_guess = 0
         do j=1, logT_points
            logT = logT_min + dlogT*(j-1)
            T = 10 ** logT

            do i=1,logW_points
               logW = logW_min + dlogW*(i-1)
               
               logPgas = logW + 4*logT
               Pgas = 10 ** logPgas

               call_number = call_number + 1
               select case(which_eos)
                  case (scvh_eos)
                     call scvh_4_PT( &
                        logW, logT, T, logPgas, Pgas, &
                        abar, zbar, X, Z, logRho_guess, &
                        logRho, logE, logS, chiRho, chiT, &
                        Cp, Cv, dE_dRho, dS_dT, dS_dRho, &
                        mu, log_free_e, gamma1, gamma3, grad_ad, eta, call_number, ierr)
                  case (mesa_eos)
                     call mesa_4_PT( &
                        logW, logT, T, logPgas, Pgas, &
                        abar, zbar, X, Z, logRho_guess, &
                        logRho, logE, logS, chiRho, chiT, &
                        Cp, Cv, dE_dRho, dS_dT, dS_dRho, &
                        mu, log_free_e, gamma1, gamma3, grad_ad, eta, call_number, ierr)
                  case (eosPT)
                     call eos_4_Pgas_T( &
                        logW, logT, T, logPgas, Pgas, &
                        abar, zbar, X, Z, logRho_guess, &
                        logRho, logE, logS, chiRho, chiT, &
                        Cp, Cv, dE_dRho, dS_dT, dS_dRho, &
                        mu, log_free_e, gamma1, gamma3, grad_ad, eta, call_number, ierr)
                  case default
                     stop 'bad value for which_eos'
               end select
                              
               k = 0
               k = k+1; output_values(i,j,k) = logRho
               k = k+1; output_values(i,j,k) = logPgas
               k = k+1; output_values(i,j,k) = logE
               k = k+1; output_values(i,j,k) = logS
               k = k+1; output_values(i,j,k) = chiRho
               k = k+1; output_values(i,j,k) = chiT
               k = k+1; output_values(i,j,k) = safe_log10(Cp)
               k = k+1; output_values(i,j,k) = safe_log10(Cv)
               energy = 10**logE
               S = 10**logS
               Rho = 10**logRho
               k = k+1; output_values(i,j,k) = dE_dRho*Rho/energy
               k = k+1; output_values(i,j,k) = dS_dT*T/S
               k = k+1; output_values(i,j,k) = dS_dRho*Rho/S
               k = k+1; output_values(i,j,k) = mu
               k = k+1; output_values(i,j,k) = log_free_e
               k = k+1; output_values(i,j,k) = gamma1
               k = k+1; output_values(i,j,k) = gamma3
               k = k+1; output_values(i,j,k) = grad_ad
               k = k+1; output_values(i,j,k) = eta
               
               if (is_bad(Cv)) then
                  write(*,1) 'logW', logW
                  write(*,1) 'logT', logT
                  write(*,1) 'logPgas', logPgas
                  write(*,1) 'logRho', logRho
                  write(*,1) 'Cv', Cv
                  stop
               end if
               
               if (is_bad(logRho) .or. ierr /= 0) then
                  output_values(i,j,1:k) = -1d-99
               else
                  logRho_guess = logRho
               end if
            
            enddo
         
         enddo
   
!$OMP PARALLEL DO PRIVATE(k)
         do k = 1, num_vals
            write(*,*) k
            write(io_first+k-1,'(e14.6)') output_values(1:logW_points,1:logT_points,k)
         end do
!$OMP END PARALLEL DO

         do i = 1, logT_points
            logT = logT_min + dlogT*(i-1)
            write(io_logT,*) logT
         end do
         close(io_logT)
      
         do j=1,logW_points
            logW = logW_min + dlogW*(j-1)
            write(io_logW,*) logW
         end do
         close(io_logW)
   
         do io=io_first,io_last
            close(io)
         end do
   
         deallocate(output_values)
      
      end subroutine make_plot_files


      subroutine Open_Plot_Outfiles(io_first, io_last, io_params, io_logW, io_logT, dir)
         integer, intent(in) :: io_first, io_params, io_logW, io_logT
         integer, intent(out) :: io_last
         character (len=*), intent(in) :: dir
         character (len=256) :: fname
         integer :: io
         
         fname = trim(dir) // '/params.data'
         open(unit=io_params,action='write',file=trim(fname))
         
         fname = trim(dir) // '/logW.data'
         open(unit=io_logW,action='write',file=trim(fname))
         
         fname = trim(dir) // '/logT.data'
         open(unit=io_logT,action='write',file=trim(fname))
         
         io = io_first-1

         fname = trim(dir) // '/logRho.data'
         io = io+1; open(unit=io,action='write',file=trim(fname))

         fname = trim(dir) // '/logPgas.data'
         io = io+1; open(unit=io,action='write',file=trim(fname))

         fname = trim(dir) // '/logE.data'
         io = io+1; open(unit=io,action='write',file=trim(fname))

         fname = trim(dir) // '/logS.data'
         io = io+1; open(unit=io,action='write',file=trim(fname))

         fname = trim(dir) // '/chiRho.data'
         io = io+1; open(unit=io,action='write',file=trim(fname))

         fname = trim(dir) // '/chiT.data'
         io = io+1; open(unit=io,action='write',file=trim(fname))

         fname = trim(dir) // '/logCp.data'
         io = io+1; open(unit=io,action='write',file=trim(fname))

         fname = trim(dir) // '/logCv.data'
         io = io+1; open(unit=io,action='write',file=trim(fname))

         fname = trim(dir) // '/dlnE_dlnRho.data'
         io = io+1; open(unit=io,action='write',file=trim(fname))

         fname = trim(dir) // '/dlnS_dlnT.data'
         io = io+1; open(unit=io,action='write',file=trim(fname))

         fname = trim(dir) // '/dlnS_dlnRho.data'
         io = io+1; open(unit=io,action='write',file=trim(fname))

         fname = trim(dir) // '/mu.data'
         io = io+1; open(unit=io,action='write',file=trim(fname))

         fname = trim(dir) // '/log_free_e.data'
         io = io+1; open(unit=io,action='write',file=trim(fname))

         fname = trim(dir) // '/gamma1.data'
         io = io+1; open(unit=io,action='write',file=trim(fname))

         fname = trim(dir) // '/gamma3.data'
         io = io+1; open(unit=io,action='write',file=trim(fname))

         fname = trim(dir) // '/grad_ad.data'
         io = io+1; open(unit=io,action='write',file=trim(fname))

         fname = trim(dir) // '/eta.data'
         io = io+1; open(unit=io,action='write',file=trim(fname))

         io_last = io
      
      end subroutine Open_Plot_Outfiles


      end module make_PTplots
