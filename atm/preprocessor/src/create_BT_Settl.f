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
      module mod_BT_Settl
         use chem_def
         use chem_lib
         use const_def
         use const_lib
         use num_lib, only: binary_search

         implicit none
         
         logical, parameter :: write_plot_files = .false., verbose = .false.

         integer, parameter :: fixed_tau=1, photosphere=2
         
         real(dp) :: tau_base, BT_Settl_M ! tau at base of atm (1, 10, or 100); [M/H]

         character(len=256) :: data_dir, output_file1, output_file2, prefix, suffix, BT_Settl_Mstr
                           
         integer :: eos_handle, kap_handle, table_type

         contains
         
         
         subroutine build_BT_Settl_tables

            integer, parameter :: num_isos = 7, nmet = 5
            integer, pointer :: chem_id(:), net_iso(:)
            real(dp), pointer :: xa(:)            
            integer :: i, j, k, ierr, i_Teff, i_logg, logg_i_lo, logg_i_hi, &
               io_out1, io_out2, io, ii, jj
            character(len=256) :: filename
            integer, parameter :: ng = 13, nT = 142, num_layers = 128 !nT=158 for jmin=4
            real(dp) :: loggs(ng), Teffs(nT), vals(num_layers), jnk1, jnk2
            real(dp) :: T(ng,nT), Pgas(ng,nT), logT(ng,nT), logPgas(ng,nT)
            real(dp) :: logT_plot(nT,ng), logPgas_plot(nT,ng), BT_Settl_Ts(nT)
            real(dp) :: X, Y, Z, XC, XN, XO, XNe, XMg, abar, zbar, z2bar
            integer, parameter :: nt_for_CK = 76, max_ng_for_CK = 11
            integer :: ng_for_CK_Teff(nt_for_CK)
            
            include 'formats.dek'
            
            ! set the Teff values
            i = 0
            do j=20,69 !4,69
               i = i+1
               BT_Settl_Ts(i) = dble(j*100)
               !write(*,3) 'T2', i, j
            end do
            do j=70,118,2
               i = i+1
               BT_Settl_Ts(i) = dble(j*100)
               !write(*,3) 'T2', i, j
            end do
            do j=120,195,5
               i = i+1
               BT_Settl_Ts(i) = dble(j*100)
               !write(*,3) 'T2', i, j
            end do
            do j=200,700,10
               i = i+1
               BT_Settl_Ts(i) = dble(j*100)
               !write(*,3) 'T2', i, j
            end do
            if (i /= nT) then
               write(*,*) 'error in setting Teff for BT_Settl'
               call mesa_error(__FILE__,__LINE__)
            end if
            
            if(table_type == fixed_tau)then
               write(*,*)
               write(*,*) 'table_type = fixed_tau'
               write(*,1) 'tau_base', tau_base     
            else
               write(*,*)
               write(*,*) 'table_type = photosphere'
            endif
      
            ierr = 0
            io = 33 ! for reading
            
            loggs(1:ng) = (/ &
               -0.5d0, 0.0d0, 0.5d0, 1.0d0, 1.5d0, 2.0d0, 2.5d0, 3.0d0, 3.5d0, 4.0d0, 4.5d0, 5.0d0, 5.5d0 /)
               
            Teffs(1:nT) = BT_Settl_Ts(1:nT)

            Pgas = -1d2
            T = -1d2
            
            call read_BT_Settl
                        
            call write_output
         
            if (write_plot_files) call write_plots
         


            contains
            
            subroutine write_plots
            use math_lib
               
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
         
               open(io,file='plot_BT_Settl/logT.data',action='write')
               write(io,'(e20.10)') logT_plot(:,:)
               close(io)
         
               open(io,file='plot_BT_Settl/logPgas.data',action='write')
               write(io,'(e20.10)') logPgas_plot(:,:)
               close(io)
         
               open(io,file='plot_BT_Settl/Teff.data',action='write')
               do j=1,nT
                  write(io,'(e20.10)') Teffs(j)
               end do
               close(io)
         
               open(io,file='plot_BT_Settl/logg.data',action='write')
               do i=1,ng
                  write(io,'(e20.10)') loggs(i)
               end do
               close(io)

            end subroutine write_plots
            
            
            subroutine write_output

               write(*,*) trim(output_file1)
               io_out1 = 33
               open(io_out1, file=output_file1)         
               write(io_out1,'(a15)',advance='no') "#Teff(K)| Pgas@" 

               if(table_type == fixed_tau)then
                  write(*,*) trim(output_file2)
                  io_out2 = 34
                  open(io_out2, file=output_file2)
                  write(io_out2,'(a15)',advance='no') "#Teff(K)|    T@" 
               endif

               do i = 1, ng
                  write(io_out1,fmt='("  log g =",f5.2," ")',advance='no') loggs(i)
                  if(table_type == fixed_tau) then
                     write(io_out2,fmt='("  log g =",f5.2," ")',advance='no') loggs(i)
                  endif
               end do !'
               write(io_out1,*)
               if(table_type == fixed_tau) write(io_out2,*)

               do j = 1, nT
                  write(io_out1,fmt='(1pe15.7)',advance='no') Teffs(j)
                  if(table_type == fixed_tau) then
                     write(io_out2,fmt='(1pe15.7)',advance='no') Teffs(j)
                  endif
                  do i = 1, ng
                     write(io_out1,fmt='(1pe15.7)',advance='no') Pgas(i,j)
                     if(table_type == fixed_tau)then
                        write(io_out2,fmt='(1pe15.7)',advance='no') T(i,j)
                     endif
                  end do
                  write(io_out1,*)
                  if(table_type == fixed_tau) write(io_out2,*)
               end do
               close(io_out1)
               if(table_type == fixed_tau) close(io_out2)
            end subroutine write_output
            
            
            subroutine read_BT_Settl
               integer :: T2
               character(len=1) :: logg_sign
               if(verbose) write(*,*) 'read BT_Settl'
               do i = 1, ng
                  !satisfy backwards naming convention
                  if(loggs(i) < 0d0)then
                     logg_sign = '+'
                  else
                     logg_sign = '-'
                  endif

                  do j = 1, nT
                     T2 = floor(BT_Settl_Ts(j))/100
                     if (T2 > 99) then
                        write(filename,'(a,i3,a,f3.1,a)') &
                           trim(prefix), T2, logg_sign, abs(loggs(i)), trim(suffix)
                        if (verbose) write(*,'(i3,a,f3.1)') trim(prefix), T2, '-', loggs(i)
                     else if (T2 > 9) then
                        write(filename,'(a,i2,a1,f3.1,a)') &
                           trim(prefix) // '0', T2, logg_sign, abs(loggs(i)), trim(suffix)
                        if (verbose) write(*,'(i2,a,f3.1)') T2, '-', loggs(i)
                     else
                        write(filename,'(a,i1,a,f3.1,a)') &
                           trim(prefix) // '00', T2, logg_sign, abs(loggs(i)), trim(suffix)
                        if (verbose) write(*,'(a,i1,a,f3.1)') '0', T2, '-', loggs(i)
                     end if
                     
                     open(io,file=trim(filename),action='read',status='old',iostat=ierr)
                     if (ierr /= 0) then
                        if(verbose) write(*,*) 'failed to open ' // trim(filename)
                        close(io)
                        cycle
                     end if
                     call read_file(i,j,ierr)
                     if (ierr /= 0) then
                        write(*,*) 'failed while reading ' // trim(filename)
                        close(io)
                        call mesa_error(__FILE__,__LINE__)
                     end if
                     close(io)
                  end do
               end do
            end subroutine read_BT_Settl
         
         
            subroutine read_file(i,j,ierr)
               integer, intent(in) :: i, j
               integer, intent(out) :: ierr
               integer :: layers, base_layer=1
               real(dp) :: Pg, Pe, interp_factor(2)
               ierr = 0
               read(io,*) !skip header info
               read(io,'(i5)') layers
               if (layers /= num_layers) then
                  write(*,*) 'layers /= n', i, j, layers
                  call mesa_error(__FILE__,__LINE__)
               end if
               read(io,*) !skip header info
               
               !read tau array
               if(verbose) print *, 'tau'
               call read_vals(.false.,ierr) ! tau
               if (ierr /= 0) return
               if(table_type == fixed_tau)then
               !locate tau_base layer and set linear interpolation factors
                  base_layer = binary_search(num_layers, base_layer, vals, tau_base)
                  if(base_layer >= num_layers .or. base_layer <= 1) then
                     !try one quick fix:
                     base_layer = binary_search(num_layers-10, num_layers/2, vals(6:num_layers-5), tau_base)
                     print *, 'got here!'
                     print *, 'base_layer = ', base_layer
                  endif
                  if(base_layer >= num_layers .or. base_layer <= 1) then
                     write(*,*) 'problem with ', trim(filename)
                     return
                  endif
                  interp_factor(2) = (tau_base - vals(base_layer))/(vals(base_layer+1) - vals(base_layer))
                  interp_factor(1) = 1d0 - interp_factor(2)
                  if(verbose)then
                     print *, ''
                     print *, 'filename = ', trim(filename)
                     print *, 'table_type = fixed_tau'
                     print *, 'base_layer = ', base_layer
                     print *, 'tau_base = ', tau_base
                     print *, 'tau(base_layer) = ', vals(base_layer)
                     print *, 'tau(base_layer+1) = ', vals(base_layer+1)
                     print *, ''
                  endif
               endif

               !read T array
               if(verbose) print *, 'T'
               call read_vals(.false.,ierr) ! temperature
               if (ierr /= 0) return
               if(table_type == fixed_tau)then
                  T(i,j) = interp_factor(1)*vals(base_layer) + interp_factor(2)*vals(base_layer+1)
               else !doing photosphere, find T=Teff
                  base_layer = binary_search(num_layers, base_layer, vals, Teffs(j))
                  !in a very few cases there is something funny in the temperature profile that
                  !makes it hard for the search to find the right T range for Teff
                  ! experience suggests that the problem areas are near the beginning or end of the array
                  ! so we try to look again after removing the inner and outer 5 points -- same test is 
                  ! employed above in the tau search, though it does not seem to experience the same problem
                  if(base_layer >= num_layers .or. base_layer <= 1) then
                     !try one quick fix:
                     base_layer = binary_search(num_layers-10, num_layers/2, vals(6:num_layers-5), Teffs(j))
                     print *, 'got here!'
                     print *, 'base_layer = ', base_layer
                  endif
                  if(base_layer >= num_layers .or. base_layer <= 1) then
                     write(*,*) 'problem with ', trim(filename)
                     return
                  endif
                  interp_factor(2) = (Teffs(j) - vals(base_layer))/(vals(base_layer+1)-vals(base_layer))
                  interp_factor(1) = 1d0 - interp_factor(2)
                  if(verbose)then
                     print *, ''
                     print *, 'filename = ', trim(filename)
                     print *, 'table_type = photosphere'
                     print *, 'base_layer = ', base_layer
                     print *, 'Teff = ', Teffs(j)
                     print *, 'T(base_layer) = ', vals(base_layer)
                     print *, 'T(base_layer+1) = ', vals(base_layer+1)
                     print *, ''
                  endif
               endif

               if(verbose) print *, 'flxrad'
               call read_vals(.true.,ierr) ! flxrad
               if (ierr /= 0) return
               if(verbose) print *, 'terad'
               call read_vals(.true.,ierr) ! terad
               if (ierr /= 0) return
               if(verbose) print *, 'bkmean'
               call read_vals(.true.,ierr) ! bkmean
               if (ierr /= 0) return
               if(verbose) print *, 'jkmean'
               call read_vals(.true.,ierr) ! jkmean
               if (ierr /= 0) return
               if(verbose) print *, 'fkmean'
               call read_vals(.true.,ierr) ! fkmean
               if (ierr /= 0) return
               if(verbose) print *, 'rkmean'
               call read_vals(.true.,ierr) ! rkmean
               if (ierr /= 0) return
               if(verbose) print *, 'ul_db'
               call read_vals(.true.,ierr) ! ul_db
               if (ierr /= 0) return
               if(verbose) print *, 'ul_dbc'
               call read_vals(.true.,ierr) ! ul_dbc
               if (ierr /= 0) return
               if(verbose) print *, 'ul_gam'
               call read_vals(.true.,ierr) ! ul_gam
               if (ierr /= 0) return
               if(verbose) print *, 'ul_bet'
               call read_vals(.true.,ierr) ! ul_bet
               if (ierr /= 0) return

               if(verbose) print *, 'pgas'
               call read_vals(.true.,ierr) ! pgas
               if (ierr /= 0) return
               Pg = interp_factor(1)*vals(base_layer) + interp_factor(2)*vals(base_layer+1)

               if(verbose) print *, 'pe'
               call read_vals(.true.,ierr) ! pe
               if (ierr /= 0) return
               Pe = interp_factor(1)*vals(base_layer) + interp_factor(2)*vals(base_layer+1)
               Pgas(i,j) = Pg + Pe


               if(verbose) print *, 'done with read_file'


            end subroutine read_file
         
         
            subroutine read_vals(skip,ierr) !'
               logical, intent(in) :: skip
               integer, intent(out) :: ierr
               integer :: k
               ierr = 0
               if (skip) read(io,*,iostat=ierr) ! skip line
               if (ierr /= 0) return
               read(io,'(3d24.17)',iostat=ierr) vals(1:num_layers)
               if(verbose) then
                  do k=1,num_layers,4
                     write(*,'(i4,1p4e20.10)') k, vals(k), vals(k+1), vals(k+2), vals(k+3)
                  enddo
               endif
            end subroutine read_vals       
            
         end subroutine build_BT_Settl_tables

      end module mod_BT_Settl
      

      program create_BT_Settl
      use mod_BT_Settl
      
      integer :: ierr, i
      integer, parameter :: num_Zs = 11
      real(dp) :: MH(num_Zs)
      character(len=8) :: MH_Str(num_Zs) 
      
      ierr = 0

      MH = (/ -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.3, 0.5 /)
      MH_Str = (/ '-4.0', '-3.5', '-3.0', '-2.5', '-2.0', '-1.5', '-1.0', '-0.5', '-0.0', '+0.3', '+0.5' /)

      data_dir = '../../data'

      call const_init
      call chem_init(data_dir, 'isotopes.data_approx', ierr)
      if (ierr /= 0) stop 'chem_init'
     
      do i=1,num_Zs

         write(*,*)
         write(*,*)
         write(*,*) 'doing [M/H] = ', trim(MH_Str(i))
         BT_Settl_M = MH(i)  ! log10(Z/Zsun)
         BT_Settl_Mstr = MH_Str(i)

         prefix = 'atm_input_data/BT-Settl-M' // trim(BT_Settl_Mstr) // '/lte'
         !prefix = '/astro/dotter/PHOENIX/BT-Settl-M' // trim(BT_Settl_Mstr) // '/lte'
         suffix = trim(BT_Settl_Mstr) // '.BT-Settl.20'

         ! do fixed tau tables
         table_type = fixed_tau
 
         tau_base = 1d2
         output_file1 = 'atm_data/BT_Settl_M' // trim(BT_Settl_Mstr) // '_tau100_Pgas.data'
         output_file2 = 'atm_data/BT_Settl_M' // trim(BT_Settl_Mstr) // '_tau100_T.data'
         call build_BT_Settl_tables
      
         tau_base = 1d1
         output_file1 = 'atm_data/BT_Settl_M' // trim(BT_Settl_Mstr) // '_tau10_Pgas.data'
         output_file2 = 'atm_data/BT_Settl_M' // trim(BT_Settl_Mstr) // '_tau10_T.data'
         call build_BT_Settl_tables
      
         tau_base = 1d0
         output_file1 = 'atm_data/BT_Settl_M' // trim(BT_Settl_Mstr) // '_tau1_Pgas.data'
         output_file2 = 'atm_data/BT_Settl_M' // trim(BT_Settl_Mstr) // '_tau1_T.data'
         call build_BT_Settl_tables
     
         ! do photosphere table
         table_type = photosphere
         output_file1 = 'atm_data/BT_Settl_M' // trim(BT_Settl_Mstr) // '_phot_Pgas.data'
         call build_BT_Settl_tables

      enddo

      end program create_BT_Settl
