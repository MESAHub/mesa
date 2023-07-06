      module test_ionization_support
      
      
      use ionization_lib
      use chem_def
      use chem_lib, only: chem_init
      use const_lib
      use math_lib
      use utils_lib, only: mesa_error




      implicit none
      
      
      
      
      integer :: ierr
      logical, parameter :: use_cache = .true.
      character (len=256) :: my_mesa_dir, file_prefix, Z1_suffix, &
         in_dir, out_dir_ion, out_dir_eosDT, out_dir_eosPT

      contains
      
      
      subroutine do_test(quietly)
         logical, intent(in) :: quietly

         file_prefix = 'ion'
         Z1_suffix = '_CO_1'
         ierr = 0
         
         my_mesa_dir = '../..'         
         call const_init(my_mesa_dir,ierr)     
         if (ierr /= 0) then
            write(*,*) 'const_init failed'
            call mesa_error(__FILE__,__LINE__)
         end if        

         call math_init()
         
         call chem_init('isotopes.data', ierr)
         if (ierr /= 0) then
            write(*,*) 'FATAL ERROR: failed in chem_init'
            call mesa_error(__FILE__,__LINE__)
         end if
      
         call ionization_init(file_prefix, Z1_suffix, '', use_cache, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
      
         if (.false.) then
            in_dir = 'eos_macdonald'
            out_dir_ion = '' ! 'ionization_data'
            out_dir_eosDT = 'eosDT_data'
            out_dir_eosPT = 'eosPT_data'
            call create_ion_table_files( &
               in_dir, out_dir_ion, out_dir_eosDT, out_dir_eosPT)
            stop
         end if
      
         !call create_table_plot_files; stop
         !call Build_Plots; stop
      
         call do_test_Paquette_ionization(quietly)
         call test_fe56_in_he4(quietly)
         call do_test_eval_ionization(quietly)
         
      end subroutine do_test

      subroutine test_fe56_in_he4(quietly)
         logical, intent(in) :: quietly
         integer :: ierr
         1 format(a40,f16.7)
         ierr = 0
         if (ierr /= 0) then
            write(*,*) 'eval_charge_of_Fe56_in_He4 failed'
            call mesa_error(__FILE__,__LINE__)
         end if
         if (.not. quietly) then
            write(*,*)
            write(*,*) 'test_fe56_in_he4'
         end if
         call do1_fe56_in_he4(30.5d0, 7.9d0, quietly)
         call do1_fe56_in_he4(26d0, 7d0, quietly)
         call do1_fe56_in_he4(24.5d0, 6.1d0, quietly)         
      end subroutine test_fe56_in_he4

      subroutine do1_fe56_in_he4(log10_ne, log10_T, quietly)
         real(dp), intent(in) :: log10_ne, log10_T
         logical, intent(in) :: quietly
         real(dp) :: z
         integer :: ierr
         1 format(a40,3f16.7)
         ierr = 0
         z = eval_charge_of_Fe56_in_He4(log10_ne, log10_T, ierr)
         if (ierr /= 0) then
            write(*,1) 'eval_charge_of_Fe56_in_He4 failed', log10_ne, log10_T
            call mesa_error(__FILE__,__LINE__)
         end if
         if (quietly) return
         write(*,1) 'z for Fe56 in He4', z, log10_ne, log10_T
      end subroutine do1_fe56_in_he4

      subroutine do_test_eval_ionization(quietly)
         use ionization_def, only: num_ion_vals, ion_result_names
         logical, intent(in) :: quietly
         real(dp) :: Z, X, Rho, log10Rho, T, log10T
         real(dp) :: res(num_ion_vals)
         integer :: ierr, i

         include 'formats'

         ierr = 0
         Z = 0.0188d0
         X = 0.725d0
         T = 3.2345529591587989D+04
         log10T = log10(T)
         rho = 9.0768775206067858D-06
         log10Rho = log10(rho)

         if (.not. quietly) then
            write(*,*)
            write(*,*) 'do_test_eval_ionization'
            write(*,1) 'log10T', log10T
            write(*,1) 'log10rho', log10rho
            write(*,1) 'Z', Z
            write(*,1) 'X', X
         end if

         call eval_ionization(Z, X, Rho, log10Rho, T, log10T, res, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         
         if (quietly) return

         write(*,*)
         do i=1,num_ion_vals
            write(*,1) trim(ion_result_names(i)), res(i)
         end do
         write(*,*)
         write(*,*) 'done'
         write(*,*)
         
      end subroutine do_test_eval_ionization

      subroutine do_test_Paquette_ionization(quietly)
         logical, intent(in) :: quietly
         real(dp) :: abar, free_e, T, log10_T, rho, log10_rho, &
            typical_charge, actual
         integer :: cid, ierr
         
         include 'formats'

         abar = 1.4641872501488922D+00
         free_e = 7.7809739936525557D-01
         T = 3.2345529591587989D+05
         log10_T = log10(T)
         rho = 9.0768775206067858D-05
         log10_rho = log10(rho)
         
         if (.not. quietly) then
            write(*,*)
            write(*,*) 'do_test_ionization'
            write(*,1) 'log10_T =', log10_T
            write(*,1) 'log10_rho =', log10_rho
         end if



         if (.false.) then
            cid = ini58
            typical_charge = 1.0000000000000000D+01
            actual = eval_typical_charge(cid, abar, free_e, T, log10_T, rho, log10_rho)
            if (.not. quietly) write(*,1) trim(chem_isos% name(cid)) // ' typical charge = ', actual
            if (abs(actual - typical_charge) > 1d-6) call mesa_error(__FILE__,__LINE__)
            stop 'test ionization'
         end if



         cid = ih1
         typical_charge = 1.0000000000000000D+00
         actual = eval_typical_charge(cid, abar, free_e, T, log10_T, rho, log10_rho)
         if (.not. quietly) write(*,1) trim(chem_isos% name(cid)) // ' typical charge = ', actual
         if (abs(actual - typical_charge) > 1d-6) call mesa_error(__FILE__,__LINE__)

         cid = ihe4
         typical_charge = 2.0000000000000000D+00
         actual = eval_typical_charge(cid, abar, free_e, T, log10_T, rho, log10_rho)
         if (.not. quietly) write(*,1) trim(chem_isos% name(cid)) // ' typical charge = ', actual
         if (abs(actual - typical_charge) > 1d-6) call mesa_error(__FILE__,__LINE__)

         cid = io16
         typical_charge = 6.0000000000000000D+00
         actual = eval_typical_charge(cid, abar, free_e, T, log10_T, rho, log10_rho)
         if (.not. quietly) write(*,1) trim(chem_isos% name(cid)) // ' typical charge = ', actual
         if (abs(actual - typical_charge) > 1d-6) call mesa_error(__FILE__,__LINE__)

         cid = ife52
         typical_charge = 1.0000000000000000D+01
         actual = eval_typical_charge(cid, abar, free_e, T, log10_T, rho, log10_rho)
         if (.not. quietly) write(*,1) trim(chem_isos% name(cid)) // ' typical charge = ', actual
         if (abs(actual - typical_charge) > 1d-6) call mesa_error(__FILE__,__LINE__)
         
         if (quietly) return
         
         write(*,*)
         write(*,*) 'done'
         write(*,*)
         
      end subroutine do_test_Paquette_ionization 


      subroutine do_other_test
         real(dp) :: abar, free_e, T, log10_T, rho, log10_rho, &
            typical_charge, actual
         integer :: cid, ierr
         
         1 format(a40,1pe26.16)
         
         call chem_init('isotopes.data_approx', ierr)
         if (ierr /= 0) then
            write(*,*) 'FATAL ERROR: failed in chem_init'
            call mesa_error(__FILE__,__LINE__)
         end if

         abar = 1.4641872501488922D+00
         free_e = 7.7809739936525557D-01
         T = 3.2345529591587989D+05
         log10_T = 5.5098142662145326D+00
         rho = 9.0768775206067858D-05
         log10_rho =   -4.0420635246937922D+00
         
         write(*,*)
         write(*,1) 'log10_T =', log10_T
         write(*,1) 'log10_rho =', log10_rho


         cid = ih1
         typical_charge = 1.0000000000000000D+00
         actual = eval_typical_charge(cid, abar, free_e, T, log10_T, rho, log10_rho)
         write(*,1) trim(chem_isos% name(cid)) // ' typical charge = ', actual
         if (abs(actual - typical_charge) > 1d-6) call mesa_error(__FILE__,__LINE__)

         cid = ihe4
         typical_charge = 2.0000000000000000D+00
         actual = eval_typical_charge(cid, abar, free_e, T, log10_T, rho, log10_rho)
         write(*,1) trim(chem_isos% name(cid)) // ' typical charge = ', actual
         if (abs(actual - typical_charge) > 1d-6) call mesa_error(__FILE__,__LINE__)

         cid = io16
         typical_charge = 6.0000000000000000D+00
         actual = eval_typical_charge(cid, abar, free_e, T, log10_T, rho, log10_rho)
         write(*,1) trim(chem_isos% name(cid)) // ' typical charge = ', actual
         if (abs(actual - typical_charge) > 1d-6) call mesa_error(__FILE__,__LINE__)

         cid = ife52
         typical_charge = 1.0000000000000000D+01
         actual = eval_typical_charge(cid, abar, free_e, T, log10_T, rho, log10_rho)
         write(*,1) trim(chem_isos% name(cid)) // ' typical charge = ', actual
         if (abs(actual - typical_charge) > 1d-6) call mesa_error(__FILE__,__LINE__)
         
         write(*,*)
         write(*,*) 'done'
         write(*,*)

         
      end subroutine do_other_test
      
      
      subroutine Build_Plots
         use const_lib
         use utils_lib, only: mkdir
         integer :: ierr
         character(len=256) :: data_dir, dir
      
         real(dp) :: log_ne_min,log_ne_max,logT_min,logT_max,dlog_ne,dlogT,lgRho,log_ne,logT
      
         integer log_ne_points, logT_points
         integer i,j,k,info,io,io_first,io_last,io_log_ne,io_logT,num_vals

         integer, parameter :: io_unit0 = 40
      
         !real(dp) :: M, X, Z, kap

         real(dp), allocatable :: output_values(:,:,:)
         
         data_dir = '../../data'    

   !..set the sample size
         log_ne_points = 300
         logT_points = 300
            
            
   !..set the ranges

         log_ne_max = 31
         log_ne_min = 24
         logT_min = 6
         logT_max = 8
            
   !..open the output files
         io_log_ne = io_unit0
         io_logT = io_unit0+2
         io_first = io_unit0+3

         dir = 'plot_data'
         call mkdir(dir)
         call Open_Plot_Outfiles(io_first, io_last, io_log_ne, io_logT, dir)
         num_vals  = io_last - io_first + 1
         allocate(output_values(logT_points,log_ne_points,num_vals))
         
   !..get the results

         dlog_ne = (log_ne_max - log_ne_min)/(log_ne_points-1)
         dlogT = (logT_max - logT_min)/(logT_points-1)

         do j=1,logT_points
            logT = logT_min + dlogT*(j-1)
            do i=1, log_ne_points
               log_ne = log_ne_min + dlog_ne*(i-1)
               call Plot_one( &
                  i, j, log_ne, logT, &
                  output_values, num_vals, logT_points, log_ne_points, info)
            end do
         end do
   
         write(*,*) 'write files'
         do k = 1, num_vals
            write(*,*) k
            write(io_first+k-1,'(e24.16)') output_values(1:logT_points,1:log_ne_points,k)
         end do

         do i = 1, log_ne_points
            log_ne = log_ne_min + dlog_ne*(i-1)
            write(io_log_ne,*) log_ne
         end do
         close(io_log_ne)
      
         do j=1,logT_points
            logT = logT_min + dlogT*(j-1)
            write(io_logT,*) logT
         end do
         close(io_logT)
   
         do io=io_first,io_last
            close(io)
         end do
   
         deallocate(output_values)

      end subroutine Build_Plots


      subroutine Plot_one( &
            i, j, log_ne, logT, &
            output_values, num_vals, logT_points, log_ne_points, ierr)
         integer, intent(in) :: i, j, num_vals, logT_points, log_ne_points
         real(dp), intent(in) :: log_ne, logT
         real(dp), intent(out) :: output_values(logT_points,log_ne_points,num_vals)
         integer, intent(out) :: ierr

         real(dp) :: z
         integer :: k
         
         include 'formats'
                     
         ierr = 0
          
         z = eval_charge_of_Fe56_in_He4(log_ne, logT, ierr)
         if (ierr /= 0) then
            write(*,1) 'eval_charge_of_Fe56_in_He4 failed log_ne, logT', log_ne, logT
            return
         end if
         
         k = 0
         k = k+1; output_values(i,j,k) = z

      end subroutine Plot_one
      
      
      subroutine Open_Plot_Outfiles(io_first, io_last, io_log_ne, io_logT, dir)
         integer, intent(in) :: io_first, io_log_ne, io_logT
         integer, intent(out) :: io_last
         character (len=256), intent(in) :: dir
         character (len=256) :: fname
         integer :: io
         
         fname = trim(dir) // '/log_ne.data'
         open(unit=io_log_ne,file=trim(fname))
         
         fname = trim(dir) // '/logT.data'
         open(unit=io_logT,file=trim(fname))
         
         io = io_first-1
         
         fname = trim(dir) // '/z_fe56_he4.data'
         io = io+1; open(unit=io,file=trim(fname))  
         
         io_last = io
      
      end subroutine Open_Plot_Outfiles


      end module test_ionization_support




