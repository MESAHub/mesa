! ***********************************************************************
!
!   Copyright (C) 2009-2019  Bill Paxton, Frank Timmes & The MESA Team
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

      module test_net_support
      use chem_def
      use chem_lib, only: chem_init
      use math_lib
      use net_def
      use net_lib
      use const_def
      use rates_def
      use test_net_do_one
      use utils_lib, only: mesa_error
      
      implicit none
      
      integer, parameter :: max_files = 20, max_cnt = 100000      
      
      
      character (len=256) :: cache_suffix
      integer :: num_reactions
      
      integer, dimension(:), pointer :: net_iso, chem_id, isos_to_show
      
      integer, pointer :: reaction_table(:)
      integer, pointer :: rates_to_show(:)

      real(dp), dimension(:), pointer :: rho_vector, T_vector

      integer :: nrates_to_show, nisos_to_show, net_handle
      


      contains
      

      subroutine do_test_net(do_plots, symbolic)
         logical, intent(in) :: do_plots, symbolic
         call set_composition(species, xin)
         eta = 0
         if (do_plots) then
            call Create_Plot_Files(species, xin)
         else
            call Do_One_Net(symbolic)
         end if
      end subroutine do_test_net 
      
      
      subroutine load_libs
         use const_lib
         use const_def, only: mesa_dir
         use chem_lib
         use rates_lib, only: rates_init, rates_warning_init
         integer :: ierr
         character (len=64) :: my_mesa_dir
         
         my_mesa_dir = '../..'         
         call const_init(my_mesa_dir,ierr)     
         if (ierr /= 0) then
            write(*,*) 'const_init failed'
            call mesa_error(__FILE__,__LINE__)
         end if        

         call math_init()

         ierr = 0
         call chem_init('isotopes.data', ierr)
         if (ierr /= 0) then
            write(*,*) 'chem_init failed'
            call mesa_error(__FILE__,__LINE__)
         end if
      
         call rates_init('reactions.list', '', 'rate_tables', .false., .false., &
                        '', '', '',ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         
         call rates_warning_init(.true., 10d0)
      
      end subroutine load_libs
      
      
      subroutine test_net_setup(net_file_in)
         character (len=*), intent(in) :: net_file_in
         integer, pointer :: r_id(:)
         type(Net_General_Info), pointer :: g

         integer :: info, i, ierr
         
         include 'formats'
         
         net_file = net_file_in

         allocate(net_iso(num_chem_isos), isos_to_show(num_chem_isos), chem_id(num_chem_isos))
         
         call net_init(ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         
         net_handle = alloc_net_handle(ierr)
         if (ierr /= 0) then
            write(*,*) 'alloc_net_handle failed'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         call net_start_def(net_handle, ierr)
         if (ierr /= 0) then
            write(*,*) 'net_start_def failed'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         call read_net_file(net_file, net_handle, ierr)
         if (ierr /= 0) then
            write(*,*) 'read_net_file failed ', trim(net_file)
            call mesa_error(__FILE__,__LINE__)
         end if
         
         call net_finish_def(net_handle, ierr)
         if (ierr /= 0) then
            write(*,*) 'net_finish_def failed'
            call mesa_error(__FILE__,__LINE__)
         end if
      
         allocate(reaction_id(rates_reaction_id_max), reaction_table(rates_reaction_id_max))
         allocate(rates_to_show(rates_reaction_id_max))
               
         cache_suffix = ''
         call net_setup_tables(net_handle, cache_suffix, info)
         if (ierr /= 0) then
            write(*,*) 'net_setup_tables failed'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         call net_ptr(net_handle, g, ierr)
         if (ierr /= 0) then
            write(*,*) 'net_ptr failed'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         species = g% num_isos
         num_reactions = g% num_reactions
         
         call get_chem_id_table(net_handle, species, chem_id, ierr)
         if (ierr /= 0) then
            write(*,*) 'get_chem_id_table failed'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         call get_net_iso_table(net_handle, net_iso, ierr)
         if (ierr /= 0) then
            write(*,*) 'get_net_iso_table failed'
            call mesa_error(__FILE__,__LINE__)
         end if
                  
         call get_reaction_id_table(net_handle, num_reactions, reaction_id, ierr)
         if (ierr /= 0) then
            write(*,*) 'get_reaction_id_table failed'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         call get_net_reaction_table(net_handle, reaction_table, ierr)
         if (ierr /= 0) then
            write(*,*) 'get_net_reaction_table failed'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         call do_test_net_alloc(species)

      end subroutine test_net_setup        
              
      subroutine do_test_net_alloc(species)
         integer, intent(in) :: species
         allocate( &
            xin(species), xin_copy(species), d_eps_nuc_dx(species), dxdt(species), &
            d_dxdt_dRho(species), d_dxdt_dT(species), d_dxdt_dx(species, species))
      end subroutine do_test_net_alloc


      subroutine Do_One_Net(symbolic)
         logical, intent(in) :: symbolic
         integer :: i, id
         call do1_net(net_handle, symbolic)
      end subroutine Do_One_Net

      subroutine Setup_eos(handle)
 ! allocate and load the eos tables
         use eos_def
         use eos_lib
         integer, intent(out) :: handle
         integer :: ierr
         logical, parameter :: use_cache = .true.
         call eos_init('', use_cache, ierr)
         if (ierr /= 0) then
            write(*,*) 'eos_init failed in Setup_eos'
            call mesa_error(__FILE__,__LINE__)
         end if         
         handle = alloc_eos_handle(ierr)
         if (ierr /= 0) then
            write(*,*) 'failed trying to allocate eos handle'
            call mesa_error(__FILE__,__LINE__)
         end if      
      end subroutine Setup_eos
      
      
      subroutine test_net_cleanup
         call do_test_net_cleanup
         call free_net_handle(net_handle)
      end subroutine test_net_cleanup
           
      subroutine do_test_net_cleanup
         deallocate(xin)
         deallocate(d_eps_nuc_dx)
         deallocate(dxdt)
         deallocate(d_dxdt_dRho)
         deallocate(d_dxdt_dT)
         deallocate(d_dxdt_dx)
      end subroutine do_test_net_cleanup
      
      
      subroutine change_net(new_net_file)
         character (len=*), intent(in) :: new_net_file
         call test_net_cleanup
         call test_net_setup(new_net_file)
      end subroutine change_net


      subroutine Create_Plot_Files(species, xin)
         use utils_lib, only: mkdir
         integer, intent(in) :: species
         real(dp) :: xin(species)
         type (net_info)  :: n

         integer:: i, j, k, ierr, num_reactions
         real(dp) :: T, logT, logT_min, logT_max
         real(dp) :: logRho_min, logRho_max, dlogT, dlogRho
         integer :: logT_points, logRho_points
         integer :: io, io_first, io_last, io_rho, io_tmp, io_params, num_out
         character (len=256) :: fname, dir

         real(dp), allocatable :: output_values(:, :, :)
         type(Net_General_Info), pointer :: g



 ! full range
         logT_max = 9.1d0
         logT_min = 6d0
         logRho_min = -3d0
         logRho_max = 10d0

 ! oxygen burning range
         logT_max = 9.5d0
         logT_min = 8d0
         logRho_min = -3d0
         logRho_max = 12d0

 ! test FL
         logT_max = 9d0
         logT_min = 7d0
         logRho_min = 5d0
         logRho_max = 10.2d0

 ! test FL
         logT_max = 8.4d0
         logT_min = 7.8d0
         logRho_min = 2d0
         logRho_max = 6d0

 ! test C+C
         logT_max = 10d0
         logT_min = 7.5d0
         logRho_min = 6d0
         logRho_max = 12d0

         logT_points = 251
         logRho_points = 251
         
         dir = 'plot_data'
         call mkdir(dir)
         write(*, *) trim(dir)
         
 01   format(E30.22)

         dlogT = (logT_max-logT_min)/(logT_points-1)
         dlogRho = (logRho_max-logRho_min)/(logRho_points-1)         

         io_params = 40
         io_rho = 41
         io_tmp = 42
         io_first = 43

         fname = trim(dir) // '/' // 'params.data'
         open(unit=io_params, file=trim(fname))
         write(io_params, '(6f16.6)')  &
               xin(net_iso(ih1)), xin(net_iso(ihe4)), xin(net_iso(ic12)),   &
               xin(net_iso(in14)), xin(net_iso(io16))
         close(io_params)

         fname = trim(dir) // '/' // 'rho.data'
         open(unit=io_rho, file=trim(fname))

         fname = trim(dir) // '/' // 'tmp.data'
         open(unit=io_tmp, file=trim(fname))

         io = io_first-1
         io_last = Open_Files(io, dir)
         num_out = io_last - io_first + 1

         call get_net_ptr(net_handle, g, ierr)
         if(ierr/=0) return
         
         num_reactions = g% num_reactions

         allocate(output_values(logRho_points, logT_points, num_out))
         
!xx$OMP PARALLEL DO PRIVATE(logT, T, j)
         do j=1, logT_points
            logT = logT_min + dlogT*(j-1)
            T = exp10(logT)

            call do_inner_loop(species, num_reactions, logT, T, j, output_values, n, xin,  &
                     logRho_points, logRho_min, dlogRho)
            
         end do
!xx$OMP END PARALLEL DO

         write(*, *) 'write the files'


 ! write out the results
         do j=1, logRho_points
            write(io_rho, 01) logRho_min + dlogRho*(j-1)
         end do
         close(io_rho)

         do i=1, logT_points
            write(io_tmp, 01) logT_min + dlogT*(i-1)
         enddo
         close(io_tmp)
         
!$OMP PARALLEL DO PRIVATE(k)
         do k = 1, num_out
            write(*, *) k
            write(io_first+k-1, '(e14.6)') output_values(:, :, k)
         end do
!$OMP END PARALLEL DO
         
         do io=io_first, io_last
            close(io)
         end do
         
      end subroutine Create_Plot_Files


      subroutine do_inner_loop(species, num_reactions, logT, T, j, output_values, n, xin,  &
               logRho_points, logRho_min, dlogRho)
         integer, intent(in) :: species, num_reactions
         real(dp), intent(in) :: logT, T
         integer, intent(in) :: j, logRho_points
         type (Net_info) :: n
         real(dp), intent(OUT) :: output_values(:, :, :)
         real(dp), intent(in) :: xin(species), logRho_min, dlogRho

         integer :: i
         real(dp) :: logRho, Rho

         do i=1, logRho_points
            logRho = logRho_min + dlogRho*(i-1)
            Rho = exp10(logRho)
            call do_one_net_eval(species, num_reactions, logT, T, logRho, Rho,  &
                        i, j, output_values, n, xin)
         enddo
         
      end subroutine do_inner_loop
      
      
      subroutine do_one_net_eval(species, num_reactions, logT, T, logRho, Rho,  &
               i, j, output_values, n, xin)
         use chem_lib, only:composition_info
         integer, intent(in) :: species, num_reactions
         real(dp), intent(in) :: logT, T, logRho, Rho
         integer, intent(in) :: i, j
         type (Net_General_Info), pointer :: g => null()
         type(net_info) :: n
         real(dp), intent(OUT) :: output_values(:, :, :)
         real(dp), intent(in) :: xin(species)
      
         real(dp) :: z, abar, zbar, z2bar, z53bar, ye, sum, mx, weak_rate_factor
         real(dp), target :: rate_factors_a(num_reactions)
         real(dp), pointer :: rate_factors(:)

         real(dp) :: eps_nuc
         real(dp) :: d_eps_nuc_dT
         real(dp) :: d_eps_nuc_dRho
         real(dp) :: d_eps_nuc_dx(species) 
 ! partial derivatives wrt mass fractions
         
         real(dp) :: dxdt(species)  
 ! rate of change of mass fractions caused by nuclear reactions
         real(dp) :: d_dxdt_dRho(species)
         real(dp) :: d_dxdt_dT(species)
         real(dp) :: d_dxdt_dx(species, species)  
         real(dp) :: eps_nuc_categories(num_categories)  

         integer :: info, k, h1, he4, chem_id(species)
         real(dp) :: xh, xhe, mass_correction
         real(dp), dimension(species) :: dabar_dx, dzbar_dx, dmc_dx
         logical :: skip_jacobian
         
         rate_factors => rate_factors_a
         
         h1 = net_iso(ih1)
         he4 = net_iso(ihe4)

         g => n% g
      
         call get_chem_id_table(net_handle, species, chem_id, info)
         if (info /= 0) call mesa_error(__FILE__,__LINE__)

         call composition_info( &
            species, chem_id, xin, xh, xhe, z, abar, zbar, z2bar, z53bar, ye,  &
            mass_correction, sum, dabar_dx, dzbar_dx, dmc_dx)
     
         rate_factors(:) = 1
         weak_rate_factor = 1
         skip_jacobian = .false.
         
         call net_get(net_handle, skip_jacobian, n, species, num_reactions,  &
                  xin, T, logT, Rho, logRho,  &
                  abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
                  rate_factors, weak_rate_factor, &
                  std_reaction_Qs, std_reaction_neuQs, &
                  eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
                  dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx,  &
                  screening_mode,  &
                  eps_nuc_categories, eps_neu_total, &
                  info)
         if (info /= 0) then
            write(*, *) 'bad result from net_get'
            call mesa_error(__FILE__,__LINE__)
         end if

         sum = 0d0
         mx = 0d0
         do k = species, 1, -1
            if (abs(dxdt(k)) > mx) mx = abs(dxdt(k))
            sum = sum + dxdt(k)
         end do
         if (mx <= 0) mx = 1d0
         if (sum-sum /= 0) then
            write(*, *) logRho, logT
            call mesa_error(__FILE__,__LINE__)
         end if

         k =   1; output_values(i, j, k) = safe_log10(eps_nuc)
         
         k = k+1; output_values(i, j, k) =  sum / max(1d-20, mx)

         k = k+1; output_values(i, j, k) =  d_eps_nuc_dT * T / max(1d-20, eps_nuc)
         k = k+1; output_values(i, j, k) =  d_eps_nuc_dRho * Rho / max(1d-20, eps_nuc)
         k = k+1; output_values(i, j, k) =  safe_log10(abs(d_eps_nuc_dx(h1)))
         k = k+1; output_values(i, j, k) =  safe_log10(abs(d_eps_nuc_dx(he4)))
      
         if (k > max_files) then
            write(*, *) 'need to enlarge max_files'
            call mesa_error(__FILE__,__LINE__)
         end if

      end subroutine do_one_net_eval
      
      integer function Open_Files(io_start, dir)
         integer, intent(in) :: io_start
         character (len=256), intent(in) :: dir
         character (len=256) fname
         integer :: io
         io = io_start
         
         fname = trim(dir) // '/' // 'log_net_eps.data'
         io = io + 1; open(unit=io, file=trim(fname))

         fname = trim(dir) // '/' // 'sum_dxdt_div_max_dxdt.data'
         io = io+1; open(unit=io, file=trim(fname))

         fname = trim(dir) // '/' // 'd_lneps_dlnT.data'
         io = io + 1; open(unit=io, file=trim(fname))

         fname = trim(dir) // '/' // 'd_lneps_dlnRho.data'
         io = io + 1; open(unit=io, file=trim(fname))

         fname = trim(dir) // '/' // 'd_eps_nuc_dxh1.data'
         io = io + 1; open(unit=io, file=trim(fname))

         fname = trim(dir) // '/' // 'd_eps_nuc_dxhe4.data'
         io = io + 1; open(unit=io, file=trim(fname))
         
         Open_Files = io
         
      end function Open_Files


      subroutine set_composition(species, xin)
         integer, intent(in) :: species
         real(dp), intent(OUT) :: xin(species)
   
         real(dp) :: sum
         integer :: i, adjustment_iso
         
         eta = 0
         
         adjustment_iso = net_iso(img24)

         if (net_file == 'basic.net') then   
         
            adjustment_iso = net_iso(img24)
            
            xin = 0
            xin(net_iso(ih1))  =   0.655186E+00 ! h1   
            xin(net_iso(ihe4)) =   0.31002164D+00 ! he4  
            xin(net_iso(ic12)) =   0.002725D-01 ! c12
            xin(net_iso(in14)) =   0.203101D-01 + 0.612124D-06  + 0.109305D-02  + 0.356004D-04   
            xin(net_iso(io16)) =   0.094000D-01 ! o16 
            xin(net_iso(ine20)) =  0.162163D-02 ! ne20 
            xin(net_iso(img24)) =  0.658226D-25 ! mg24
            xin(net_iso(ihe3)) =   0.201852D-02 ! he3  
            
         else if (net_file == 'o18_and_ne22.net') then
         
            adjustment_iso = net_iso(img24)
            
            xin = 0
            xin(net_iso(ih1))  =   0.655186E+00 
            xin(net_iso(ihe4)) =   0.31002164D+00
            xin(net_iso(ic12)) =   0.002725D-01 
            xin(net_iso(in14)) =   0.203101D-01 + 0.612124D-06  + 0.109305D-02
            xin(net_iso(io16)) =   0.094000D-01
            xin(net_iso(io18)) =   1d-20
            xin(net_iso(ine20)) =  0.162163D-02  
            xin(net_iso(img24)) =  0.658226D-25 
            xin(net_iso(ine22)) =   0.201852D-02
            
         else if (net_file == 'pp_extras.net') then
         
            adjustment_iso = net_iso(img24)
            
            xin = 0
            xin(net_iso(ih1))  =   0.655186E+00 ! h1   
            xin(net_iso(ihe4)) =   0.31002164D+00 ! he4  
            xin(net_iso(ic12)) =   0.002725D-01 ! c12
            xin(net_iso(in14)) =   0.203101D-01 + 0.612124D-06  + 0.109305D-02  + 0.356004D-04   
            xin(net_iso(io16)) =   0.094000D-01 ! o16 
            xin(net_iso(ine20)) =  0.162163D-02 ! ne20 
            xin(net_iso(img24)) =  0.658226D-25 ! mg24
    
            xin(net_iso(ih2))  =   0.632956D-17 ! h2   
            xin(net_iso(ihe3)) =   0.201852D-02 ! he3  
            xin(net_iso(ili7)) =   0.664160D-15 ! li7  
            xin(net_iso(ibe7)) =   0.103866D-15 ! be7  
      
         else if (net_file == 'cno_extras.net') then

            adjustment_iso = net_iso(img24)
            xin = 0
            xin(net_iso(ih1))   =   0.173891680788D-01 ! h1   
            xin(net_iso(ihe4))  =   0.963245225401D+00 ! he4  
            xin(net_iso(ic12))  =   0.238935745993D-03 ! c12  
            xin(net_iso(in14))  =   0.134050688300D-01 ! n14  
            xin(net_iso(io16))  =   0.268791618452D-03 ! o16  
            xin(net_iso(ine20)) =   0.180001692845D-02 ! ne20 
            xin(net_iso(img24)) =   0.353667702698D-02 ! mg24 

            xin(net_iso(ic13))  =   0.717642727071D-04 ! c13  
            xin(net_iso(in13))  =   0.370732258156D-09 ! n13  
            xin(net_iso(in15))  =   0.450484708137D-06 ! n15  
            xin(net_iso(io14))  =   0.100000000000D-49 ! o14  
            xin(net_iso(io15))  =   0.874815374966D-10 ! o15  
            xin(net_iso(if17))  =   0.100000000000D-49 ! f17  
            xin(net_iso(if18))  =   0.100000000000D-49 ! f18  
            xin(net_iso(ine18)) =   0.100000000000D-49 ! ne18 
            xin(net_iso(ine19)) =   0.100000000000D-49 ! ne19 
            xin(net_iso(img22)) =   0.439011547696D-04 ! mg22             
            
         else if (net_file == 'pp_cno_extras_o18_ne22.net') then
         
            adjustment_iso = net_iso(img24)
            xin = 0
            xin(net_iso(ih1))   =   0.173891680788D-01 ! h1   
            xin(net_iso(ihe4))  =   0.963245225401D+00 ! he4  
            xin(net_iso(ic12))  =   0.238935745993D-03 ! c12  
            xin(net_iso(in14))  =   0.134050688300D-01 ! n14  
            xin(net_iso(io16))  =   0.268791618452D-03 ! o16  
            xin(net_iso(ine20)) =   0.180001692845D-02 ! ne20 
            xin(net_iso(img24)) =   0.353667702698D-02 ! mg24 

            xin(net_iso(ic13))  =   0.717642727071D-04 ! c13  
            xin(net_iso(in13))  =   0.370732258156D-09 ! n13  
            xin(net_iso(in15))  =   0.450484708137D-06 ! n15  
            xin(net_iso(io14))  =   0.100000000000D-49 ! o14  
            xin(net_iso(io15))  =   0.874815374966D-10 ! o15  
            xin(net_iso(if17))  =   0.100000000000D-49 ! f17  
            xin(net_iso(if18))  =   0.100000000000D-49 ! f18  
            xin(net_iso(ine18)) =   0.100000000000D-49 ! ne18 
            xin(net_iso(ine19)) =   0.100000000000D-49 ! ne19 
            xin(net_iso(img22)) =   0.439011547696D-04 ! mg22             
            
         else if (net_file == 'approx21_cr60_plus_co56.net') then
         
            adjustment_iso = net_iso(img24)
            xin = 0
            xin(net_iso(ihe4))=     3.4555392534813939D-01
            xin(net_iso(io16))=     1.9367778420430937D-01
            xin(net_iso(ih1))=     1.9052931367172501D-01
            xin(net_iso(isi28))=     9.1023936032064601D-02
            xin(net_iso(ini56))=     5.8492459653295005D-02
            xin(net_iso(is32))=     4.1416642476999908D-02
            xin(net_iso(ine20))=     2.0927782332477558D-02
            xin(net_iso(ic12))=     2.0589617104448312D-02
            xin(net_iso(img24))=     1.8975914108406104D-02
            xin(net_iso(iar36))=     7.0600338785939956D-03
            xin(net_iso(ica40))=     4.9182171981353726D-03
            xin(net_iso(ife52))=     3.1618820806005280D-03
            xin(net_iso(in14))=     1.9878460082760952D-03
            xin(net_iso(ife56))=     7.1668776621130190D-04
            xin(net_iso(ico56))=     4.4974971099921181D-04
            xin(net_iso(ife54))=     3.3432423890988422D-04
            xin(net_iso(icr48))=     1.7602626907769762D-04
            xin(net_iso(ihe3))=     5.1650097191631545D-06
            xin(net_iso(iti44))=     2.6452422203389906D-06
            xin(net_iso(icr60))=     4.7653353095937982D-08
            xin(net_iso(iprot))=     1.2739022407246250D-11
            xin(net_iso(ineut))=     0.0000000000000000D+00
            
         else if (net_file == 'approx21_cr60_plus_fe53_fe55_co56.net') then
         
            adjustment_iso = net_iso(img24)
            xin = 0
            xin(net_iso(ihe4))=     3.4555392534813939D-01
            xin(net_iso(io16))=     1.9367778420430937D-01
            xin(net_iso(ih1))=     1.9052931367172501D-01
            xin(net_iso(isi28))=     9.1023936032064601D-02
            xin(net_iso(ini56))=     5.8492459653295005D-02
            xin(net_iso(is32))=     4.1416642476999908D-02
            xin(net_iso(ine20))=     2.0927782332477558D-02
            xin(net_iso(ic12))=     2.0589617104448312D-02
            xin(net_iso(img24))=     1.8975914108406104D-02
            xin(net_iso(iar36))=     7.0600338785939956D-03
            xin(net_iso(ica40))=     4.9182171981353726D-03
            xin(net_iso(ife52))=     3.1618820806005280D-03
            xin(net_iso(in14))=     1.9878460082760952D-03
            xin(net_iso(ife56))=     7.1668776621130190D-04
            xin(net_iso(ico56))=     4.4974971099921181D-04
            xin(net_iso(ife54))=     3.3432423890988422D-04
            xin(net_iso(icr48))=     1.7602626907769762D-04
            xin(net_iso(ihe3))=     5.1650097191631545D-06
            xin(net_iso(iti44))=     2.6452422203389906D-06
            xin(net_iso(icr60))=     4.7653353095937982D-08
            xin(net_iso(iprot))=     1.2739022407246250D-11
            xin(net_iso(ife53))=     0.0000000000000000D+00
            xin(net_iso(ife55))=     0.0000000000000000D+00
            xin(net_iso(ineut))=     0.0000000000000000D+00
            
         else if (net_file == 'approx21_new.net' .or. &
                     net_file == 'approx21_old.net' .or. &
                     net_file == 'approx21_plus_co56.net' .or. &
                     net_file == 'approx21.net') then
         
            adjustment_iso = net_iso(img24)
            xin = 0
             xin(net_iso(ife56))=     8.0387021484318166D-01
             xin(net_iso(ife54))=     1.6096648736760832D-01
             xin(net_iso(icr56))=     2.9480945535920525D-02
              xin(net_iso(ihe4))=     4.8624637161320565D-03
             xin(net_iso(ini56))=     2.8376270731890360D-04
             xin(net_iso(isi28))=     1.5018628906135739D-04
              xin(net_iso(is32))=     1.1613271635573457D-04
             xin(net_iso(iprot))=     1.1139431633673653D-04
             xin(net_iso(ica40))=     5.3688377473494185D-05
             xin(net_iso(iar36))=     5.2702831822567062D-05
             xin(net_iso(ife52))=     3.7866504131935185D-05
             xin(net_iso(icr48))=     5.8401123037667974D-06
             xin(net_iso(ineut))=     4.9141227703118397D-06
             xin(net_iso(iti44))=     1.5085038561154746D-06
              xin(net_iso(io16))=     9.5384019609049255D-07
             xin(net_iso(img24))=     6.3808207717725580D-07
              xin(net_iso(ic12))=     2.9048656673991868D-07
             xin(net_iso(ine20))=     9.6468865023609427D-09
              xin(net_iso(ihe3))=     4.6203603862263096D-80
              xin(net_iso(in14))=     7.5867472225841235D-99

         else
            
            write(*,*) 'net_file ' // trim(net_file)
            call mesa_error(__FILE__,__LINE__,'set_composition: do not recognize net_file')
      
         end if

         sum = 0d0
         do i=1, species
            if (xin(i) < 1d-99) xin(i) = 1d-99
            sum = sum + xin(i)
         end do
         if (abs(1d0-sum) > 1d-4) write(*, *) 'change abundance sum by', 1d0-sum
         xin(adjustment_iso) = xin(adjustment_iso) + (1d0 - sum)
         if (xin(adjustment_iso) < 0d0) call mesa_error(__FILE__,__LINE__,'error in sum of abundances')

      end subroutine set_composition
      
      
      subroutine read_test_data(filename, n, rho_vec, T_vec, ierr)
 ! the data files have columns of mass, radius, density, temp
         use utils_lib
         character (len=*), intent(in) :: filename
         integer, intent(out) :: n
         real(dp), dimension(:), pointer :: rho_vec, T_vec ! to be allocated and filled
         integer, intent(out) :: ierr
         
         integer :: iounit, i
         real(dp) :: junk
         
         ierr = 0
         open(newunit=iounit, file=trim(filename), action='read', iostat=ierr)
         if (ierr /= 0) then
            write(*, *) 'failed to open ', trim(filename)
            return
         end if
         
         i = 0
         do
            read(unit=iounit, fmt=*, iostat=ierr) junk, junk, junk, junk
            if (ierr == 0) then
               i = i+1; cycle
            end if
            ierr = 0
            n = i
            exit
         end do
         rewind(iounit)
         
         allocate(rho_vec(n), T_vec(n), stat=ierr); if (ierr /= 0) return
         
         do i=1, n
            read(iounit, *) junk, junk, rho_vec(i), T_vec(i)
         end do
      
         close(iounit)
         
      end subroutine read_test_data

      
      subroutine Do_One_Test(net_file, do_timing)
         use chem_lib, only:composition_info
         use rates_lib
         character (len=*), intent(in) :: net_file
         logical, intent(in) :: do_timing
         call Do_One_Testcase(net_file, do_timing, .false.)
      end subroutine Do_One_Test

      
      subroutine Do_One_Test_and_show_Qs(net_file, do_timing)
         use chem_lib, only:composition_info
         use rates_lib
         character (len=*), intent(in) :: net_file
         logical, intent(in) :: do_timing
         call Do_One_Testcase(net_file, do_timing, .true.)
      end subroutine Do_One_Test_and_show_Qs
      
      
      subroutine Do_One_Testcase(net_file, do_timing, show_Qs)
         use chem_lib, only:composition_info
         use rates_lib
         character (len=*), intent(in) :: net_file
         logical, intent(in) :: do_timing, show_Qs
         
         real(dp) :: logRho, logT, Rho, T, xsum, Q1, Q2, &
           eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, weak_rate_factor, &
            dvardx, dvardx_0, dx_0, err, var_0, xdum, &
           eps_nuc_categories(num_categories), xh, xhe, mass_correction !approx_abar, approx_zbar
         integer :: i, j, k, info, ierr
         integer :: j_dx, j_dx_sink
         integer :: adjustment_iso, ir_c12_c12_to_he4_ne20, ir_he4_ne20_to_c12_c12
         real(dp), dimension(:), pointer :: d_eps_nuc_dx, dabar_dx, dzbar_dx, dmc_dx
         real(dp), pointer :: rate_factors(:), &
            actual_Qs(:), actual_neuQs(:)       
         logical, pointer :: from_weaklib(:)
         logical :: skip_jacobian, doing_d_dlnd, doing_dx
         real(dp), dimension(:), pointer :: &
            rate_raw, rate_raw_dT, rate_raw_dRho, &
            rate_screened, rate_screened_dT, rate_screened_dRho
         type(net_info) :: n
         
         include 'formats'
         
         write(*,*) 'Do_One_Test ' // trim(net_file)
         
         if (do_timing) call mesa_error(__FILE__,__LINE__,'no support for do_timing')
                  
         call test_net_setup(net_file)
            
         ierr = 0
         
         write(*,*) 'species', species
         
         allocate(d_eps_nuc_dx(species), dabar_dx(species), dzbar_dx(species), dmc_dx(species))
         
         info = 0

         
         allocate(&
               rate_factors(num_reactions), &
               actual_Qs(num_reactions), &
               actual_neuQs(num_reactions), &
               from_weaklib(num_reactions), &
               stat=info)
         if (info /= 0) call mesa_error(__FILE__,__LINE__)
         
         rate_factors(:) = 1
         weak_rate_factor = 1
         
         if (.false.) then ! get neutrino Q
         
            Q1 = eval_neutrino_Q(img22, is30)
            write(*,1) 'Qneu mg22->s30', Q1
         
            Q1 = eval_neutrino_Q(is30, ini56)
            write(*,1) 'Qneu s30->ni56', Q1
         
            Q1 = eval_neutrino_Q(ica38, ini56)
            write(*,1) 'Qneu ca38->ni56', Q1
         
            Q1 = eval_neutrino_Q(ini56, ige64)
            write(*,1) 'Qneu ni56->ge64', Q1

            Q1 = eval_neutrino_Q(ige64, ise68)
            write(*,1) 'Qneu ge64->se68', Q1

            Q1 = eval_neutrino_Q(ise68, ikr72)
            write(*,1) 'Qneu se68->kr72', Q1

            Q1 = eval_neutrino_Q(ikr72, isr76)
            write(*,1) 'Qneu kr72->sr76', Q1

            Q1 = eval_neutrino_Q(isr76, isn104)
            write(*,1) 'Qneu sr76->sn104', Q1

            stop
         end if
         
         if (.false.) then ! get reaction Q
 ! co55 -> fe55
            Q1 = isoB(ife55) - isoB(ico55)
            write(*,1) 'Q co55->fe55', Q1
            stop
         end if

         xin = 0
         eta = 0

         if (net_file == 'mesa_201.net') then
            
            nrates_to_show = 2
         
            rates_to_show(1:nrates_to_show) = (/  &
               ir_ar36_ag_ca40, &
               ir_ca40_ga_ar36  /)

                                  xin(net_iso(ife56))=     6.3551174304779179D-01
                                  xin(net_iso(icr52))=     1.0518849309100423D-01
                                  xin(net_iso(ini60))=     6.8579809304547934D-02
                                  xin(net_iso(ife55))=     3.0800958229563902D-02
                                  xin(net_iso(imn55))=     2.1373990699858084D-02
                                  xin(net_iso(ico57))=     2.0785551124804885D-02
                                  xin(net_iso(ife57))=     1.6083543077536507D-02
                                  xin(net_iso(imn53))=     1.5959192021165327D-02
                                  xin(net_iso(ico59))=     1.3346191916961030D-02
                                  xin(net_iso(ife58))=     1.2810477227656202D-02
                                  xin(net_iso(ife54))=     1.2428871023625330D-02
                                  xin(net_iso(ini62))=     1.1831182533246063D-02
                                  xin(net_iso(icr53))=     8.1041573656346795D-03
                                  xin(net_iso(ini61))=     5.1462544919734727D-03
                                  xin(net_iso(imn54))=     4.8562102089927551D-03
                                  xin(net_iso(icr54))=     3.4271335428296490D-03
                                  xin(net_iso(ini59))=     2.9811021846265560D-03
                                  xin(net_iso(ini58))=     2.8713349650366189D-03
                                   xin(net_iso(iv51))=     2.0988150935152424D-03
                                  xin(net_iso(ico58))=     2.0282210582857861D-03
                                  xin(net_iso(icr51))=     1.2727926750247761D-03
                                  xin(net_iso(icr50))=     3.5421790633561727D-04
                                  xin(net_iso(icu63))=     3.0335040211002022D-04
                                  xin(net_iso(iti48))=     2.8512639104984498D-04
                                  xin(net_iso(iti50))=     2.8394752519368671D-04
                                  xin(net_iso(ico56))=     2.4310701594699283D-04
                                  xin(net_iso(ico60))=     1.5319574812323961D-04
                                   xin(net_iso(ihe4))=     1.3780017854631601D-04
                                  xin(net_iso(imn56))=     1.1066944504563417D-04
                                   xin(net_iso(iv50))=     8.7723703257792468D-05
                                   xin(net_iso(iv49))=     8.0539507322740820D-05
                                  xin(net_iso(icu61))=     6.0898214714637941D-05
                                  xin(net_iso(iti49))=     5.8026224632063015D-05
                                  xin(net_iso(imn52))=     4.2516037186521133D-05
                                  xin(net_iso(ico61))=     4.1332712278653746D-05
                                  xin(net_iso(ini63))=     3.9285182071536010D-05
                                  xin(net_iso(ico55))=     3.3984664667118780D-05
                                  xin(net_iso(ife59))=     3.1874001320244153D-05
                                  xin(net_iso(izn64))=     2.7904240321969555D-05
                                  xin(net_iso(izn66))=     2.3416686728782276D-05
                                  xin(net_iso(icu62))=     2.1144678609352833D-05
                                  xin(net_iso(ini57))=     1.1679099363693715D-05
                                   xin(net_iso(iv52))=     1.0016381776869645D-05
                                  xin(net_iso(ini64))=     8.9564547480773243D-06
                                  xin(net_iso(iti47))=     8.7562205406651916D-06
                                  xin(net_iso(icu65))=     8.3830901771893844D-06
                                  xin(net_iso(iti46))=     6.5019528109988749D-06
                                  xin(net_iso(icu64))=     6.1556773376380492D-06
                                  xin(net_iso(iar38))=     5.0886468271422554D-06
                                  xin(net_iso(ife53))=     4.4893983930970075D-06
                                   xin(net_iso(is34))=     4.2369862519232918D-06
                                  xin(net_iso(izn65))=     4.1809380230467465D-06
                                  xin(net_iso(icr55))=     3.5647078269917385D-06
                                  xin(net_iso(imn51))=     1.3849206094166064D-06
                                  xin(net_iso(ife60))=     1.0919469947619407D-06
                                    xin(net_iso(ih1))=     9.8211241034612710D-07
                                  xin(net_iso(ica44))=     6.4721481111312556D-07
                                  xin(net_iso(isi30))=     6.0482126630410841D-07
                                  xin(net_iso(ica42))=     5.8926684360031471D-07
                                  xin(net_iso(isi28))=     5.5342250413192408D-07
                                   xin(net_iso(iv48))=     5.4911502105605942D-07
                                   xin(net_iso(iv53))=     5.4571869339657211D-07
                                  xin(net_iso(icu60))=     4.5761289432308086D-07
                                  xin(net_iso(izn63))=     4.4652045082610858D-07
                                  xin(net_iso(isc47))=     4.4099241256913694D-07
                                  xin(net_iso(ini56))=     3.9547421607331126D-07
                                  xin(net_iso(iti51))=     3.6358450091693021D-07
                                  xin(net_iso(izn62))=     3.0268870369662295D-07
                                   xin(net_iso(is32))=     3.0184105591459131D-07
                                  xin(net_iso(icr49))=     2.7706803242086906D-07
                                  xin(net_iso(isc45))=     2.5466905171154329D-07
                                   xin(net_iso(ik39))=     1.9906885458096793D-07
                                  xin(net_iso(icl35))=     1.5826438708200070D-07
                                   xin(net_iso(is33))=     1.2794750349676913D-07
                                   xin(net_iso(ip31))=     1.1805533268957108D-07
                                  xin(net_iso(icl37))=     1.0108890802543261D-07
                                  xin(net_iso(ica43))=     9.5766190955127930D-08
                                  xin(net_iso(iar36))=     8.3287350202069049D-08
                                  xin(net_iso(isi29))=     7.6157023642649312D-08
                                  xin(net_iso(isc49))=     7.2645544618960897D-08
                                  xin(net_iso(icu59))=     6.9036259006119691D-08
                                  xin(net_iso(ica40))=     6.4387988887989985D-08
                                  xin(net_iso(isc46))=     6.1294200120692868D-08
                                  xin(net_iso(iar37))=     5.0049708949178339D-08
                                  xin(net_iso(ineut))=     4.7236564548991118D-08
                                  xin(net_iso(isc48))=     3.7392395931746378D-08
                                  xin(net_iso(ife52))=     3.0982331086000098D-08
                                  xin(net_iso(icr56))=     2.9477366114308358D-08
                                  xin(net_iso(ica41))=     2.7100443580454983D-08
                                   xin(net_iso(is35))=     2.6527306613430685D-08
                                  xin(net_iso(ica45))=     2.5532508173462626D-08
                                  xin(net_iso(ico62))=     2.3608043613616947D-08
                                  xin(net_iso(iar39))=     2.3503192609361190D-08
                                  xin(net_iso(ica46))=     2.2455595751043146D-08
                                   xin(net_iso(iv47))=     1.9789662501752072D-08
                                  xin(net_iso(icu66))=     1.9285467280104737D-08
                                  xin(net_iso(icl36))=     1.8710358753573163D-08
                                   xin(net_iso(is36))=     1.6433794658525574D-08
                                   xin(net_iso(ik41))=     1.1021812817088955D-08
                                  xin(net_iso(ini65))=     1.0432196346423500D-08
                                   xin(net_iso(ip33))=     8.9625871046594568D-09
                                   xin(net_iso(ik40))=     7.2106087354006743D-09
                                  xin(net_iso(iar40))=     7.1174073521523365D-09
                                  xin(net_iso(iti45))=     5.6870435962784455D-09
                                   xin(net_iso(ip32))=     4.9537776441010461D-09
                                  xin(net_iso(icr48))=     3.0709436939798925D-09
                                  xin(net_iso(isc44))=     2.7862949914482315D-09
                                  xin(net_iso(isc43))=     1.4283233379580965D-09
                                  xin(net_iso(isi31))=     1.3747008133379580D-09
                                  xin(net_iso(iti52))=     1.2866447687101849D-09
                                  xin(net_iso(ico63))=     1.0456195723572430D-09
                                  xin(net_iso(iti44))=     7.1579364938842059D-10
                                   xin(net_iso(io16))=     5.6818053126812287D-10
                                  xin(net_iso(ial27))=     5.6117840295305725D-10
                                  xin(net_iso(ica47))=     5.4798226854137171D-10
                                  xin(net_iso(img24))=     3.8679049700264050D-10
                                  xin(net_iso(ini66))=     2.8640604416758339D-10
                                  xin(net_iso(izn61))=     2.4995195933116682D-10
                                  xin(net_iso(ica48))=     1.9626419275146166D-10
                                  xin(net_iso(ife61))=     1.8711288748163173D-10
                                   xin(net_iso(ik42))=     1.4923884785773920D-10
                                  xin(net_iso(isi32))=     1.4443919428884894D-10
                                   xin(net_iso(ip30))=     1.3908837795296053D-10
                                   xin(net_iso(ic12))=     1.3185610207511085D-10
                                   xin(net_iso(ik43))=     1.1518951971920992D-10
                                   xin(net_iso(iv54))=     8.8440974843643898D-11
                                  xin(net_iso(img26))=     8.5260453209638504D-11
                                  xin(net_iso(icl38))=     2.5853881572793256D-11
                                  xin(net_iso(isc50))=     1.7292118708137885D-11
                                  xin(net_iso(iar41))=     1.0665700916421081D-11
                                  xin(net_iso(img25))=     9.2051027558630546D-12
                                  xin(net_iso(ial28))=     8.9602126799277990D-12
                                  xin(net_iso(izn60))=     7.0658622923470296D-12
                                   xin(net_iso(ip34))=     4.2003981867966896D-12
                                  xin(net_iso(icr57))=     2.8195149736768269D-12
                                  xin(net_iso(ine20))=     1.1829006560529443D-12
                                  xin(net_iso(ife62))=     1.0149415562457171D-12
                                   xin(net_iso(ik44))=     5.5654687896999145D-13
                                   xin(net_iso(is31))=     4.0845373817884839D-13
                                   xin(net_iso(iv55))=     2.8610531204397973D-13
                                  xin(net_iso(ina23))=     2.5122330767236196D-13
                                   xin(net_iso(is37))=     2.1529809224524208D-13
                                  xin(net_iso(iti53))=     1.4129020276964862D-13
                                  xin(net_iso(iar35))=     1.2934962278023034D-13
                                  xin(net_iso(ial26))=     1.2407881461251650D-13
                                   xin(net_iso(in15))=     1.1945765698183132D-13
                                  xin(net_iso(ico64))=     8.1590671587313667D-14
                                  xin(net_iso(img27))=     6.7969591881380018D-14
                                    xin(net_iso(ih2))=     5.6613884110319004D-14
                                  xin(net_iso(ini67))=     5.1506647175988674D-14
                                  xin(net_iso(ica39))=     3.8945192023978035D-14
                                  xin(net_iso(ini55))=     3.0883773851523300D-14
                                  xin(net_iso(ine22))=     1.3808913342478939D-14
                                  xin(net_iso(ica49))=     1.0542673783683999D-14
                                  xin(net_iso(isc51))=     9.2083686560605166D-15
                                  xin(net_iso(isi27))=     8.3464191225256989D-15
                                  xin(net_iso(ine21))=     6.6697557229646203D-15
                                  xin(net_iso(ife51))=     6.5787879505834196D-15
                                   xin(net_iso(io17))=     3.8661065919908615D-15
                                   xin(net_iso(ic13))=     2.2604411883637647D-15
                                  xin(net_iso(icr58))=     2.1539938862813166D-15
                                  xin(net_iso(isi33))=     1.4888165334138879D-15
                                  xin(net_iso(ina24))=     7.2312125147882022D-16
                                  xin(net_iso(ial25))=     5.3576177397850227D-16
                                  xin(net_iso(ico65))=     4.9566812556376405D-16
                                  xin(net_iso(icr47))=     3.0194388465619186D-16
                                   xin(net_iso(in14))=     2.6891785216435022D-16
                                  xin(net_iso(ina22))=     2.5931749891034273D-16
                                   xin(net_iso(io15))=     2.5261003710407275D-16
                                  xin(net_iso(ini68))=     2.2419710841076782D-16
                                   xin(net_iso(io18))=     1.8206042351246347D-16
                                  xin(net_iso(iti43))=     9.9895985562055471D-17
                                   xin(net_iso(if19))=     7.0445166521693986D-17
                                   xin(net_iso(ihe3))=     6.8591249543794866D-17
                                  xin(net_iso(iti54))=     3.3955810406819116D-17
                                  xin(net_iso(ife63))=     3.0222399040043984D-17
                                   xin(net_iso(in13))=     2.5097133179904447D-17
                                  xin(net_iso(izn59))=     2.3782224732332168D-17
                                  xin(net_iso(img23))=     2.2784900370197838D-17
                                   xin(net_iso(if17))=     1.0052207505944401D-17
                                   xin(net_iso(if18))=     6.3013547340360942D-18
                                   xin(net_iso(iv56))=     2.1677062516769839D-18
                                  xin(net_iso(ina21))=     1.9109919103480182D-18
                                  xin(net_iso(ine23))=     1.2143041459039416D-18
                                   xin(net_iso(ili6))=     1.4203908924439747D-19
                                   xin(net_iso(ib11))=     1.1320699694157424D-19
                                   xin(net_iso(if20))=     4.4936577179455616D-20
                                  xin(net_iso(ine19))=     2.6387620496069204D-20
                                  xin(net_iso(ife64))=     1.1841283303587735D-20
                                   xin(net_iso(ib10))=     1.0494357182634838D-20
                                   xin(net_iso(ibe9))=     7.8292004126082962D-21
                                  xin(net_iso(ico66))=     6.4514234785580282D-21
                                   xin(net_iso(ili7))=     6.4231341717191173D-21
                                   xin(net_iso(in16))=     4.7443318701592232D-21
                                   xin(net_iso(ibe7))=     3.3532048357084064D-22
                                   xin(net_iso(io19))=     3.2062683754970905D-22
                                  xin(net_iso(ibe10))=     6.5331631260779713D-23
                                  xin(net_iso(ico67))=     4.9629777598135805D-24
                                  xin(net_iso(ife65))=     3.7335326045384740D-26
                                  xin(net_iso(ife66))=     8.6772104841406824D-30
                                    xin(net_iso(ib8))=     4.1439050548710376D-31
                              
                              write(*,*) 'sum xin', sum(xin(:))

                                                 logT =    9.6532818288064650D+00
                                               logRho =    7.9479966082179185D+00
                                                  eta =    2.7403163311838425D+00
                                   
                                   
               screening_mode = extended_screening
               
               call net_set_logTcut(net_handle, 0d0, 0d0, info)
               if (info /= 0) then
                  write(*,*) 'failed in net_set_logTcut'
                  call mesa_error(__FILE__,__LINE__)
               end if

         else if (net_file == 'approx21_cr60_plus_co56.net') then
            
            nrates_to_show = 2
         
            rates_to_show(1:nrates_to_show) = (/  &
                 irco56ec_to_fe56, &
                 irni56ec_to_co56 /)

                                    xin(net_iso(ineut))=    1.9615092621881698D-06
                                      xin(net_iso(ih1))=    0.0000000000000000D+00
                                    xin(net_iso(iprot))=    3.3544919533130298D-04
                                     xin(net_iso(ihe3))=    0.0000000000000000D+00
                                     xin(net_iso(ihe4))=    9.6491906054115388D-03
                                     xin(net_iso(ic12))=    4.4455299753403973D-07
                                     xin(net_iso(in14))=    0.0000000000000000D+00
                                     xin(net_iso(io16))=    7.0439210492179401D-07
                                    xin(net_iso(ine20))=    4.5369892002182717D-09
                                    xin(net_iso(img24))=    8.0018732551576689D-07
                                    xin(net_iso(isi28))=    3.3697411043955434D-04
                                     xin(net_iso(is32))=    2.7239844165292807D-04
                                    xin(net_iso(iar36))=    1.4240894544492388D-04
                                    xin(net_iso(ica40))=    1.5156321095874170D-04
                                    xin(net_iso(iti44))=    4.3979509377388036D-06
                                    xin(net_iso(icr48))=    2.8113145151721771D-05
                                    xin(net_iso(icr60))=    4.1810609055173498D-03
                                    xin(net_iso(ife52))=    1.8928784597602586D-04
                                    xin(net_iso(ife54))=    2.6996971418369603D-01
                                    xin(net_iso(ife56))=    7.0702963220409232D-01
                                    xin(net_iso(ico56))=    6.8333792315189530D-03
                                    xin(net_iso(ini56))=    8.7251484519146338D-04

                              
                              write(*,*) 'sum xin', sum(xin(:))

                                                   logT =    9.8200000000000003D+00
                                                 logRho =    8.2586740078478176D+00
                                                 eta = 0d0
                                   
                                   
               screening_mode = extended_screening
               
               call net_set_logTcut(net_handle, 0d0, 0d0, info)
               if (info /= 0) then
                  write(*,*) 'failed in net_set_logTcut'
                  call mesa_error(__FILE__,__LINE__)
               end if

         else if (net_file == 'approx21_cr60_plus_fe53_fe55_co56.net') then
            
            nrates_to_show = 2
         
            rates_to_show(1:nrates_to_show) = (/  &
                 irco56ec_to_fe56, &
                 irni56ec_to_co56 /)

                                   xin(net_iso(ihe4))=     3.4555392534813939D-01
                                   xin(net_iso(io16))=     1.9367778420430937D-01
                                    xin(net_iso(ih1))=     1.9052931367172501D-01
                                  xin(net_iso(isi28))=     9.1023936032064601D-02
                                  xin(net_iso(ini56))=     5.8492459653295005D-02
                                   xin(net_iso(is32))=     4.1416642476999908D-02
                                  xin(net_iso(ine20))=     2.0927782332477558D-02
                                   xin(net_iso(ic12))=     2.0589617104448312D-02
                                  xin(net_iso(img24))=     1.8975914108406104D-02
                                  xin(net_iso(iar36))=     7.0600338785939956D-03
                                  xin(net_iso(ica40))=     4.9182171981353726D-03
                                  xin(net_iso(ife52))=     3.1618820806005280D-03
                                   xin(net_iso(in14))=     1.9878460082760952D-03
                                  xin(net_iso(ife56))=     7.1668776621130190D-04
                                  xin(net_iso(ico56))=     4.4974971099921181D-04
                                  xin(net_iso(ife54))=     3.3432423890988422D-04
                                  xin(net_iso(icr48))=     1.7602626907769762D-04
                                   xin(net_iso(ihe3))=     5.1650097191631545D-06
                                  xin(net_iso(iti44))=     2.6452422203389906D-06
                                  xin(net_iso(icr60))=     4.7653353095937982D-08
                                  xin(net_iso(iprot))=     1.2739022407246250D-11
                                  xin(net_iso(ife53))=     0.0000000000000000D+00
                                  xin(net_iso(ife55))=     0.0000000000000000D+00
                                  xin(net_iso(ineut))=     0.0000000000000000D+00

                              
                              write(*,*) 'sum xin', sum(xin(:))

                                                 logT =    4.6233007922659333D+00
                                               logRho =   -1.0746410107891649D+01
                                                  eta =   -2.2590260158215202D+01
                                   
                                   
               screening_mode = extended_screening
               
               call net_set_logTcut(net_handle, 0d0, 0d0, info)
               if (info /= 0) then
                  write(*,*) 'failed in net_set_logTcut'
                  call mesa_error(__FILE__,__LINE__)
               end if

         else if (net_file == 'basic.net') then
            
            nrates_to_show = 1
            
            if (rates_reaction_id('rc12_to_n14') <= 0) call mesa_error(__FILE__,__LINE__,'bad reaction')
            write(*,*) 'rc12_to_n14', rates_reaction_id('rc12_to_n14')
         
            rates_to_show(1:nrates_to_show) = (/  &
            rates_reaction_id('rc12_to_n14') /)
     
                                   xin(net_iso(ihe4))=     9.8119124177708650D-01
                                   xin(net_iso(in14))=     9.8369547495994036D-03
                                   xin(net_iso(io16))=     2.9223115895360822D-03
                                  xin(net_iso(ine20))=     2.0337034688681288D-03
                                   xin(net_iso(ihe3))=     0.0000000000000000D+00
                                    xin(net_iso(ih1))=     0.0000000000000000D+00
                                    
                                    write(*,*) 'when 1st'
                                   xin(net_iso(ic12))=     2.3551202768735954D-04 ! 2.3551179217556737D-04
                                  xin(net_iso(img24))=     3.7802763872225075D-03 ! 3.7802766227342998D-03
                                    
                                    !write(*,*) 'when 2nd'
                                   !xin(net_iso(ic12))=     2.3551179217556737D-04
                                  !xin(net_iso(img24))=     3.7802766227342998D-03

                                                 logT =    7.8820961200821831D+00
                                               logRho =    5.6124502722887746D+00
                                                  eta =    1.2629003275927920D+01

                              write(*,*) 'sum xin', sum(xin(:))
                                   
                                   
               screening_mode = extended_screening
               
               call net_set_logTcut(net_handle, 0d0, 0d0, info)
               if (info /= 0) then
                  write(*,*) 'failed in net_set_logTcut'
                  call mesa_error(__FILE__,__LINE__)
               end if

         else if (net_file == 'agb.net') then
            
            nrates_to_show = 4
         
            rates_to_show(1:nrates_to_show) = (/  &
            ir_h1_h1_wk_h2, &
            ir_c13_an_o16, &
            ir_f19_ap_ne22, &
            ir_he3_ag_be7 /)
     
                     xin(net_iso(ih1))= 1
                              
                              write(*,*) 'sum xin', sum(xin(:))

                                  logT =    8.6864273893515023D+00
                                logRho =    2.0591020210828619D+00
                                   eta =   -1.4317150417353590D+01
                                   
                                   
               screening_mode = extended_screening
               
               call net_set_logTcut(net_handle, 0d0, 0d0, info)
               if (info /= 0) then
                  write(*,*) 'failed in net_set_logTcut'
                  call mesa_error(__FILE__,__LINE__)
               end if
               
               if (.false.) then
                  if (.false.) then
                     rate_factors(:) = 0
                     i = reaction_table(ir_h2_pg_he3)
                     if (i == 0) call mesa_error(__FILE__,__LINE__)
                     rate_factors(i) = 1
                  else
                     i = reaction_table(ir_ne20_ag_mg24)
                     if (i == 0) call mesa_error(__FILE__,__LINE__)
                     rate_factors(i) = 0
                     i = reaction_table(ir_o16_ag_ne20)
                     if (i == 0) call mesa_error(__FILE__,__LINE__)
                     rate_factors(i) = 0
                  end if
               end if

         

         else if (net_file == 'pp_and_cno_extras.net') then
            
            nrates_to_show = 8
         
            rates_to_show(1:nrates_to_show) = (/  &
            rates_reaction_id('r_n13_wk_c13'),                &
            rates_reaction_id('r_o15_wk_n15'),                &
            rates_reaction_id('r_f17_wk_o17'),                &
            rates_reaction_id('r_f18_wk_o18'),                &
            rates_reaction_id('r_o14_wk_n14'),                &
            rates_reaction_id('r_ne18_wk_f18'),                &
            rates_reaction_id('r_ne19_wk_f19'),                &
            ir_he4_he4_he4_to_c12 /)
     
         xin = 0
                     xin(net_iso(ih1))=     7.2265805432969643D-01
                    xin(net_iso(ihe3))=     6.7801726921522655D-04
                    xin(net_iso(ihe4))=     2.6667042876319019D-01
                    xin(net_iso(ic12))=     1.9056849943017622D-03
                    xin(net_iso(in13))=     1.4791107757148081D-04
                    xin(net_iso(in14))=     6.0253770632534619D-04
                    xin(net_iso(in15))=     1.9263919190065488D-07
                    xin(net_iso(io14))=     9.5977897394247582D-09
                    xin(net_iso(io15))=     7.6666060610240002D-08
                    xin(net_iso(io16))=     5.5886173952684358D-03
                    xin(net_iso(io17))=     2.0449560316023342D-06
                    xin(net_iso(io18))=     1.1313355448541673D-05
                    xin(net_iso(if19))=     2.1343308067891499D-07
                   xin(net_iso(ine20))=     1.2563102679570999D-03
                   xin(net_iso(img24))=     4.7858754879924638D-04
 
                                  logT =    9.6d0
                                logRho =    6.0d0
                                   rho =    7.8571498592117219D+00
                                     T =    8.5648111120065376D+06
                                  abar =    1.2655060647252907D+00
                                  zbar =    1.0901664301076275D+00
                                 z2bar =    1.3036906023574921D+00
                                    ye =    8.6144702146826535D-01
                                   eta =   -3.4387570967781595D+00                                   
                                   
               screening_mode = extended_screening

         else if (net_file == 'approx21_plus_co56.net') then
            

            nrates_to_show = 5
            
            rates_to_show(1:nrates_to_show) = (/  &
            ir_v47_pa_ti44,            &
            ir_mn51_pa_cr48, &
            ir_ar36_ap_k39, &
            ir_co55_pa_fe52,            &
            ir_s32_ap_cl35             &
            /)
                     xin = 0
            
                                  xin(net_iso(ife54))=     7.8234742556602999D-01
                                  xin(net_iso(isi28))=     7.8210084821085060D-02
                                  xin(net_iso(ini56))=     5.2306555890846963D-02
                                   xin(net_iso(is32))=     5.1718599300602117D-02
                                  xin(net_iso(iar36))=     1.3861579457866108D-02
                                  xin(net_iso(ica40))=     1.2347894507489101D-02
                                  xin(net_iso(ife56))=     4.8313904464514874D-03
                                  xin(net_iso(ife52))=     3.9677414859767722D-03
                                  xin(net_iso(icr48))=     2.9870241199304463D-04
                                   xin(net_iso(ihe4))=     4.0616997047809087D-05
                                  xin(net_iso(iti44))=     3.3724322385639036D-05
                                  xin(net_iso(iprot))=     1.7517084938626562D-05
                                  xin(net_iso(img24))=     1.1248868356795997D-05
                                   xin(net_iso(io16))=     6.1628768742062784D-06
                                   xin(net_iso(ic12))=     7.4857838132808694D-07
                                  xin(net_iso(ine20))=     5.6045205217740520D-09
                                  xin(net_iso(icr56))=     1.7593106525983871D-09
                                  xin(net_iso(ineut))=     1.9978632764577629D-11
                                   xin(net_iso(in14))=     0.0000000000000000D+00
                                   xin(net_iso(ihe3))=     0.0000000000000000D+00
                                    xin(net_iso(ih1))=     0.0000000000000000D+00
                              write(*,*) 'test case sum xin', sum(xin(1:species))
                              

                                                 logT =    9.5806070583042597D+00
                                               logRho =    7.1251356937727763D+00
                                                  eta =    6.3504500633747440D-01

                                                    T =    3.8072119794259882D+09
                                                  rho =    1.3339381513098311D+07
                                                 abar =    4.8254908015575154D+01
                                                 zbar =    2.3420437256420026D+01
                                                z2bar =    5.7180065624152712D+02
                                                   ye =    4.8534829345982072D-01
                        screening_mode = extended_screening

         else
            
            write(*, *) 'need to define setup for net_file ', trim(net_file)
            call mesa_error(__FILE__,__LINE__,'Do_One_Test')
         
         end if
         
         Rho = exp10(logRho)
         T = exp10(logT)
         
         write(*, *)
         write(*, *)
         
         info = 0
         
         ierr = 0
         call composition_info( &
            species, chem_id, xin, xh, xhe, z, abar, zbar, z2bar, z53bar, ye,  &
            mass_correction, xsum, dabar_dx, dzbar_dx, dmc_dx)
         
         write(*,'(a40,d26.16)') 'xh', xh
         write(*,'(a40,d26.16)') 'xhe', xhe
         write(*,'(a40,d26.16)') 'abar', abar
         write(*,'(a40,d26.16)') 'zbar', zbar
         write(*,'(a40,d26.16)') 'z2bar', z2bar
         write(*,'(a40,d26.16)') 'ye', ye
         do i = 1, species
            write(*,'(a40,i6,d26.16)')  'init x ' // trim(chem_isos% name(chem_id(i))), i, xin(i)
         end do
         write(*,'(A)')
         write(*,'(a40,d26.16)') 'logT', logT
         write(*,'(a40,d26.16)') 'logRho', logRho
         write(*,'(a40,d26.16)') 'eta', eta
         
         skip_jacobian = .false.

         n% screening_mode = screening_mode

         if (.false.) then
            write(*,*) 'call net_get_rates_only'
            call net_get_rates_only( &
               net_handle, n, species, num_reactions,  &
               xin, T, logT, Rho, logRho,  &
               abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
               rate_factors, weak_rate_factor, &
               std_reaction_Qs, std_reaction_neuQs, &
               screening_mode, &
               ierr)
            call mesa_error(__FILE__,__LINE__,'net_get_rates_only')
         end if

         call net_get_with_Qs( &
               net_handle, skip_jacobian, n, species, num_reactions,  &
               xin, T, logT, Rho, logRho,  &
               abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
               rate_factors, weak_rate_factor, &
               std_reaction_Qs, std_reaction_neuQs, &
               eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
               dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx,  &
               screening_mode,  &
               eps_nuc_categories, eps_neu_total,  &
               actual_Qs, actual_neuQs, from_weaklib, info)
         if (info /= 0) then
            write(*,1) 'logT', logT
            write(*,1) 'logRho', logRho
            write(*, *) 'bad return from net_get'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         
         
         if (.true.) then ! dfridr tests for partials
         
         
            !doing_d_dlnd = .true.
            doing_d_dlnd = .false.
            doing_dx = .false.

            j_dx = 22
            var_0 = dxdt(j_dx)
         
            if (doing_dx) then
               j_dx = 20 ! fe56 
               j_dx_sink = 17 ! cr56
               dx_0 = xin(j_dx)*1d-6
               dvardx_0 = d_eps_nuc_dx(j_dx)
               write(*,1) 'xin(fe56)', xin(j_dx)
               write(*,1) 'xin(cr56)', xin(j_dx_sink)
               write(*,1) 'eps_nuc', eps_nuc
               write(*,'(A)')
               !call mesa_error(__FILE__,__LINE__,'testing')
            else if (doing_d_dlnd) then
               dx_0 = max(1d-14, abs(logRho*ln10*1d-6))
               dvardx_0 = d_eps_nuc_dRho*Rho ! d_dxdt_dRho(1)*Rho ! analytic value of partial
            else
               dx_0 = max(1d-14, abs(logT*ln10*1d-6))
               !dvardx_0 = d_eps_nuc_dT*T
               dvardx_0 = d_dxdt_dT(j_dx)*T
            end if
            err = 0d0
            dvardx = dfridr(dx_0,err)
            xdum = (dvardx - dvardx_0)/max(abs(dvardx_0),1d-50)
            write(*,1) 'analytic, numeric, est err in numeric, rel diff', &
                  dvardx_0, dvardx, err, xdum
            if (doing_dx) then
               write(*,*) 'doing d_dx ' // trim(chem_isos% name(chem_id(j_dx)))
            else if (doing_d_dlnd) then
               write(*,*) 'doing dlnd'
            else ! doing d_dlnT
               write(*,*) 'doing dlnT'
            end if
            write(*,*) 'test net'
            write(*,'(A)')
         
            call mesa_error(__FILE__,__LINE__,'test net')
         
         
         
         
         
         end if
         
         
         

         if (show_Qs) then
            write(*,'(A)')
            write(*,1) 'logT', logT
            write(*,1) 'logRho', logRho
            write(*,'(A)')
            write(*,'(30x,4a20)') 'Q total', 'Q neutrino', 'Q total-neutrino'
            do i = 1, num_reactions
               if (from_weaklib(i)) then
                  write(*,'(a30,99f20.10)') 'weaklib ' // trim(reaction_Name(reaction_id(i))),  &
                     actual_Qs(i), actual_neuQs(i), actual_Qs(i) - actual_neuQs(i)
               else
                  write(*,'(a30,99f20.10)') trim(reaction_Name(reaction_id(i))),  &
                     actual_Qs(i), actual_neuQs(i), actual_Qs(i) - actual_neuQs(i)
               end if
            end do
            write(*,'(A)')
            stop
         end if
                  
         write(*,2) 'screening_mode', screening_mode
            
         if (.true.) then
            write(*,1) 'logT', logT
            write(*,1) 'logRho', logRho
            write(*,'(A)')
            
            write(*,1) 'eps_nuc', eps_nuc
            write(*,1) 'd_epsnuc_dlnd', d_eps_nuc_dRho*Rho
            write(*,1) 'd_epsnuc_dlnT', d_eps_nuc_dT*T
            write(*,'(A)')
            
            if (eps_nuc > 0) then
               write(*,1) 'log eps_nuc', log10(eps_nuc)
               write(*,1) 'd_lnepsnuc_dlnd', d_eps_nuc_dRho*Rho/eps_nuc
               write(*,1) 'd_lnepsnuc_dlnT', d_eps_nuc_dT*T/eps_nuc
               write(*,'(A)')
            end if

            
            !stop


            if (.false.) then

               do i = 1, species
                  write(*,1)  'd_eps_nuc_dx ' // trim(chem_isos% name(chem_id(i))), d_eps_nuc_dx(i)
               end do
               write(*,'(A)')


               do i = 1, species
                  write(*,1)  'd_dxdt_dlnRho ' // trim(chem_isos% name(chem_id(i))), d_dxdt_dRho(i)*rho
               end do
               write(*,'(A)')
            
               do i = 1, species
                  write(*,1)  'd_dxdt_dlnT ' // trim(chem_isos% name(chem_id(i))), d_dxdt_dT(i)*T
               end do
               write(*,'(A)')
            
               if (.false.) then
                  do i = 1, species
                     write(*,1)  'd_dxdt_dx(:,neut) ' // &
                        trim(chem_isos% name(chem_id(i))), d_dxdt_dx(i, net_iso(ineut))
                  end do
                  write(*,'(A)')
               end if


               do i = 1, species
                  write(*,1)  'dxdt ' // trim(chem_isos% name(chem_id(i))), dxdt(i)
               end do
               write(*,1) 'sum(dxdt)', sum(dxdt(1:species))
            
               do i=1,nrates_to_show
                  j = rates_to_show(i)
                  if (j == 0) cycle
                  write(*,1) 'd_rate_raw_dT ' // trim(reaction_Name(j)),  &
                                    rate_raw_dT(reaction_table(j))
               end do
               write(*,'(A)')
            
            end if
         
            do i=1,nrates_to_show
               j = rates_to_show(i)
               if (j == 0) cycle
               write(*,1) 'rate_raw ' // trim(reaction_Name(j)),  &
                                 rate_raw(reaction_table(j))
            end do
            write(*,'(A)')
            
            do i=1,nrates_to_show
               j = rates_to_show(i)
               if (j == 0) cycle
               write(*,1) 'rate_screened ' // trim(reaction_Name(j)), &
                                 rate_screened(reaction_table(j))
            end do
            write(*,'(A)')

            do i=1,nrates_to_show
               j = rates_to_show(i)
               if (j == 0) cycle
               write(*,1) 'Q total ' // trim(reaction_Name(j)), &
                                 actual_Qs(reaction_table(j))
            end do
            write(*,'(A)')
            
            do i=1,nrates_to_show
               j = rates_to_show(i)
               if (j == 0) cycle
               write(*,1) 'Q neutrino ' // trim(reaction_Name(j)), &
                                 actual_neuQs(reaction_table(j))
            end do
            write(*,'(A)')
            
            do i=1,nrates_to_show
               j = rates_to_show(i)
               if (j == 0) cycle
               write(*,1) 'Q total-neutrino ' // trim(reaction_Name(j)), &
                                 actual_Qs(reaction_table(j)) - actual_neuQs(reaction_table(j))
            end do
            write(*,'(A)')
            
            if (.false.) then
            
               do i = 1, species
                  write(*,1)  'x ' // trim(chem_isos% name(chem_id(i))), xin(i)
               end do
               write(*,'(A)')
            
               do i = 1, species
                  if (-dxdt(i) > 1d-90)  &
                     write(*,1)  'x/dxdt ' // trim(chem_isos% name(chem_id(i))), xin(i)/dxdt(i)
               end do
               write(*,'(A)')
         
               do i = 1, num_categories
                  if (abs(eps_nuc_categories(i)) < 1d-20) cycle
                  write(*,1)  'eps_nuc_cat ' // trim(category_name(i)), eps_nuc_categories(i)
               end do
               write(*,'(A)')
               
            end if
            
            stop
         end if

         write(*,'(A)')
         write(*,'(A)')
         write(*,*) 'net_name ', trim(net_file)
         write(*,*) 'species', species
         write(*,1) 'abar =', abar
         write(*,1) 'zbar =', zbar
         write(*,1) 'z2bar =', z2bar
         write(*,1) 'ye =', ye
         write(*, *)
         do i=1,nrates_to_show
            j = rates_to_show(i)
            if (j == 0) cycle
            if (reaction_table(j) == 0) then
               write(*,*) 'missing reaction_table(j) for ' // trim(reaction_Name(j))
               stop
            end if
         end do
         write(*,'(A)')
         write(*,1) 'eps_nuc', eps_nuc
         write(*,1) 'd_eps_nuc_dRho', d_eps_nuc_dRho
         write(*,1) 'd_eps_nuc_dT', d_eps_nuc_dT
         do j=1,species
            write(*,1) 'd_eps_nuc_dx ' // trim(chem_isos% name(chem_id(j))), d_eps_nuc_dx(j)
         end do
         write(*,'(A)')
         do j=1,species
            write(*,1) 'dxdt ' // trim(chem_isos% name(chem_id(j))), dxdt(j)
         end do
         write(*,'(A)')
         do j=1,species
            write(*,1) 'd_dxdt_dRho ' // trim(chem_isos% name(chem_id(j))), d_dxdt_dRho(j)
         end do
         write(*,'(A)')
         do j=1,species
            write(*,1) 'd_dxdt_dT ' // trim(chem_isos% name(chem_id(j))), d_dxdt_dT(j)
         end do
         write(*,'(A)')
         do j=1,species
            write(*,1) 'd_dxdt_dx(1,:) ' // trim(chem_isos% name(chem_id(j))), d_dxdt_dx(1,j)
         end do
         write(*, *)
         do j=1,num_categories
            write(*,1) trim(category_name(j)), eps_nuc_categories( j)
         end do
         write(*,'(A)')
         write(*,1) 'eta =', eta
         write(*,1) 'logT =', logT
         write(*,1) 'logRho =', logRho
         write(*,*) 'screening_mode =', screening_mode
         write(*,'(A)')
         do j=1,species
            write(*,1) 'xin(net_iso(i' // trim(chem_isos% name(chem_id(j))) // '))= ', xin(j)
         end do
         write(*,'(A)')
         write(*,1) 'sum(xin(1:species))', sum(xin(1:species))
         write(*,1) '1 - sum(xin(1:species))', 1 - sum(xin(1:species))
         write(*,'(A)')
         write(*,1) 'eps_nuc', eps_nuc
         write(*,1) 'eps_nuc_neu_total', eps_neu_total


         deallocate(rate_factors, actual_Qs, actual_neuQs, from_weaklib)

            
         contains
         
         
         real(dp) function dfridr_func(delta_x) result(val)
            real(dp), intent(in) :: delta_x
            integer :: ierr
            real(dp) :: pgas, prad, energy, entropy, var, log_var
            include 'formats'
            ierr = 0
            
            if (doing_dx) then
            
               xin_copy = xin
               xin_copy(j_dx) = xin_copy(j_dx) + delta_x
               xin_copy(j_dx_sink) = xin_copy(j_dx_sink) - delta_x
               
               call net_get_with_Qs( &
                  net_handle, skip_jacobian, n, species, num_reactions,  &
                  xin_copy, T, logT, Rho, logRho,  &
                  abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
                  rate_factors, weak_rate_factor, &
                  std_reaction_Qs, std_reaction_neuQs, &
                  eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
                  dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx,  &
                  screening_mode,    &
                  eps_nuc_categories, eps_neu_total,  &
                  actual_Qs, actual_neuQs, from_weaklib, ierr)

               write(*,1) 'xin(j)', xin_copy(j_dx)
               write(*,1) 'xin(j_sink)', xin_copy(j_dx_sink)
               write(*,1) 'eps_nuc', eps_nuc
               write(*,'(A)')
            
            else if (doing_d_dlnd) then
            
               log_var = logRho + delta_x/ln10
               var = exp10(log_var)
               
               call net_get_with_Qs( &
                  net_handle, skip_jacobian, n, species, num_reactions,  &
                  xin, T, logT, var, log_var,  &
                  abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
                  rate_factors, weak_rate_factor, &
                  std_reaction_Qs, std_reaction_neuQs, &
                  eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
                  dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx,  &
                  screening_mode,   &
                  eps_nuc_categories, eps_neu_total,  &
                  actual_Qs, actual_neuQs, from_weaklib, ierr)

            else
            
               log_var = logT + delta_x/ln10
               var = exp10(log_var)
               
               call net_get_with_Qs( &
                  net_handle, skip_jacobian, n, species, num_reactions,  &
                  xin, var, log_var, Rho, logRho,  &
                  abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
                  rate_factors, weak_rate_factor, &
                  std_reaction_Qs, std_reaction_neuQs, &
                  eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
                  dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx,  &
                  screening_mode,  &
                  eps_nuc_categories, eps_neu_total,  &
                  actual_Qs, actual_neuQs, from_weaklib, ierr)
                  
            end if
            
            if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'failed in call on net_get_with_Qs')
            !val = eps_nuc ! dxdt(1)
            val = dxdt(j_dx)
            
         end function dfridr_func


         real(dp) function dfridr(hx,err) ! from Frank
            real(dp), intent(in) :: hx
            real(dp), intent(out) :: err
            !  this routine returns the first derivative of a function func(x)
            !  at the point x, by ridders method of polynomial extrapolation.
            !  value hx is the initial step size;
            !  it should be an increment for which func changes substantially.
            !  an estimate of the error in the first derivative is returned in err.
            integer, parameter :: ntab = 20
            integer :: i,j
            real(dp) :: x,errt,fac,hh,a(ntab,ntab),xdum,ydum
            real(dp), parameter :: con2=2d0, con=sqrt(con2), big=1d50, safe=2d0
            include 'formats'
            dfridr = 0d0
            hh = hx
            ! 2nd order central difference
            a(1,1) = (dfridr_func(hh) - dfridr_func(-hh))/(2d0*hh)
            write(*,2) 'dfdx hh', 1, a(1,1), hh
            err = big
            ! succesive columns in the neville tableu will go to smaller stepsizes
            ! and higher orders of extrapolation
            do i=2,ntab
               hh = hh/con
               a(1,i) = (dfridr_func(hh) - dfridr_func(-hh))/(2d0*hh)
               write(*,2) 'dfdx hh', i, a(1,i), hh
               ! compute extrapolations of various orders; the error stratagy is to compare
               ! each new extrapolation to one order lower but both at the same stepsize
               ! and at the previous stepsize
               fac = con2
               do j=2,i
                  a(j,i) = (a(j-1,i)*fac - a(j-1,i-1))/(fac-1d0)
                  fac = con2*fac
                  errt = max(abs(a(j,i)-a(j-1,i)),abs(a(j,i)-a(j-1,i-1)))
                 ! write(*,3) 'a(j,i)', j, i, a(j,i), errt
                  if (errt <= err) then
                     err = errt
                     dfridr = a(j,i)
                     !write(*,3) 'dfridr err', i, j, dfridr, err
                  end if
               end do
               ! if higher order is worse by a significant factor safe, then bail
               if (abs(a(i,i) - a(i-1,i-1)) >= safe*err) then
                  write(*,1) 'higher order is worse', err, a(i,i), a(i-1,i-1)
                  return
               end if
            end do
         end function dfridr


      end subroutine Do_One_Testcase

      
      subroutine test_neutrino_Q
         real(dp), parameter :: Qnu_n13 = 0.714440d0 !..13n(e+nu)13c
         real(dp), parameter :: Qnu_o15 = 1.005513d0 !..15o(e+nu)15n
         real(dp), parameter :: Qnu_f17 = 1.009145d0 !..17f(e+nu)17o
         real(dp), parameter :: Qnu_f18 = 0.393075d0 !..18f(e+nu)18o   
         real(dp), parameter :: Qnu_o14 = 2.22d0 !..14o(e+nu)14n
         real(dp), parameter :: Qnu_ne18 = 1.87d0 !..18ne(e+nu)18f
         real(dp), parameter :: Qnu_ne19 = 1.25d0 !..19ne(e+nu)19f
 !real(dp), parameter :: Qnu_mg21 = 6.2d0 !..mg21(e+nu)na21
         real(dp), parameter :: Qnu_mg22 = 2.1d0 !..mg22(e+nu)na22
         
 1       format(a40, 1pd26.16)
         
         write(*, 1) 'expected Q for 13n(e+nu)13c', Qnu_n13
         write(*, 1) 'calculated Q for 13n(e+nu)13c', eval_neutrino_Q(in13, ic13)
         write(*, *)
         write(*, 1) 'expected Q for 15o(e+nu)15n', Qnu_o15
         write(*, 1) 'calculated Q for 15o(e+nu)15n', eval_neutrino_Q(io15, in15)
         write(*, *)
         write(*, 1) 'expected Q for 17f(e+nu)17o', Qnu_f17
         write(*, 1) 'calculated Q for 17f(e+nu)17o', eval_neutrino_Q(if17, io17)
         write(*, *)
         write(*, 1) 'expected Q for 18f(e+nu)18o', Qnu_f18
         write(*, 1) 'calculated Q for 18f(e+nu)18o', eval_neutrino_Q(if18, io18)
         write(*, *)
         write(*, 1) 'expected Q for 14o(e+nu)14n', Qnu_o14
         write(*, 1) 'calculated Q for 14o(e+nu)14n', eval_neutrino_Q(io14, in14)
         write(*, *)
         write(*, 1) 'expected Q for 18ne(e+nu)18f', Qnu_ne18
         write(*, 1) 'calculated Q for 18ne(e+nu)18f', eval_neutrino_Q(ine18, if18)
         write(*, *)
         write(*, 1) 'expected Q for 19ne(e+nu)19f', Qnu_ne19
         write(*, 1) 'calculated Q for 19ne(e+nu)19f', eval_neutrino_Q(ine19, if19)
         write(*, *)
         write(*, 1) 'expected Q for mg22(e+nu)na22', Qnu_mg22
         write(*, 1) 'calculated Q for mg22(e+nu)na22', eval_neutrino_Q(img22, ina22)
         write(*, *)
         stop
      
      end subroutine test_neutrino_Q
      

      
      end module test_net_support




