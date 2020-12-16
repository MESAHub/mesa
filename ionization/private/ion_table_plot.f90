      module ion_table_plot
      
      use const_def
      use ion_tables_eval
      use math_lib
      use utils_lib, only: mesa_error, mkdir

      implicit none

      
      contains
      
            
      subroutine do_create_table_plot_files

         character (len=256) :: dir
   
         real(dp) :: lgT_min, lgT_max, lgRho_min, lgRho_max, dlgT, &
            dlgRho, lgRho, lgT, Rho, T, Z, X, lgQ_min, lgQ_max
   
         integer :: lgT_points, lgRho_points
         integer :: i, j, k, ierr, io, io_first, io_last, io_params, io_rho, io_tmp, num_vals

         integer, parameter :: io_unit0 = 40

         real(dp), allocatable :: output_values(:,:,:)
         
         Z = 0.018
         X = 0.72

         !..set the sample size
         lgT_points = 300
         lgRho_points = 300
      
         !lgT_points = 100
         !lgRho_points = 100
      
         !lgT_points = 2
         !lgRho_points = 2
            
         !..set the ranges
               
         ! check opal/scvh
         lgT_max = 7.7d0
         lgT_min = 2.0d0
         lgRho_min = -5d0
         lgRho_max = 5.5d0
                  
         ! table full range
         lgT_max = 8.2
         lgT_min = 2.1
         lgQ_min = -10
         lgQ_max = 5.69
         lgRho_min = -9 ! lgQ_min + 2*lgT_min - 12
         lgRho_max = 8 ! lgQ_max + 2*lgT_max - 12
         
         ! test
         lgT_max = 7.5d0
         lgT_min = 3d0
         lgRho_max = 3d0
         lgRho_min = -7d0
            
         io_params = io_unit0
         io_rho = io_unit0+1
         io_tmp = io_unit0+2
         io_first = io_unit0+3

         dir = 'plot_data'
         call mkdir(dir)
         call open_plot_files(io_first, io_last, io_params, io_rho, io_tmp, dir)
         write(io_params, '(2(f10.6),2(i7))') Z, X, lgRho_points, lgT_points
         close(io_params)
         num_vals  = io_last - io_first + 1
         allocate(output_values(lgRho_points,lgT_points,num_vals))
         
         dlgT = (lgT_max - lgT_min)/(lgT_points-1)
         dlgRho = (lgRho_max - lgRho_min)/(lgRho_points-1)

         do j=1, lgT_points
            lgT = lgT_min + dlgT*(j-1)
            T = exp10(lgT)
            do i=1,lgRho_points
               lgRho = lgRho_min + dlgRho*(i-1)
               Rho = exp10(lgRho)
               call plot_one( &
                  i, j, Z, X, lgT, T, lgRho, Rho, output_values, &
                  num_vals, lgRho_points, lgT_points, ierr)
               if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
            end do
         end do
   
         write(*,*) 'write plot files'
         do k = 1, num_vals
            write(*,*) k
            write(io_first+k-1,'(e24.16)') output_values(1:lgRho_points,1:lgT_points,k)
         end do

         do i = 1, lgT_points
            lgT = lgT_min + dlgT*(i-1)
            write(io_tmp,*) lgT
         end do
         close(io_tmp)
      
         do j=1,lgRho_points
            lgRho = lgRho_min + dlgRho*(j-1)
            write(io_rho,*) lgRho
         end do
         close(io_rho)
   
         do io=io_first,io_last
            close(io)
         end do
   
         deallocate(output_values)

      end subroutine do_create_table_plot_files


      subroutine plot_one( &
            i, j, Z, X, lgT, T, lgRho, Rho, output_values, num_vals, lgRho_points, lgT_points, ierr)
         integer, intent(in) :: i, j, num_vals, lgRho_points, lgT_points
         real(dp), intent(in) :: Z, X, lgT, T, lgRho, Rho
         real(dp), intent(inout) :: output_values(lgRho_points,lgT_points,num_vals)
         integer, intent(out) :: ierr

         real(dp), dimension(num_ion_vals) :: res
         integer :: k
         
         include 'formats.dek'
                     
         ierr = 0
         call Get_ion_Results(Z, X, Rho, lgRho, T, lgT, res, ierr)
         if (ierr /= 0) then
            write(*,1) 'Get_ion_Results failed for lgRho, lgT', lgRho, lgT
            stop
         end if

         if (ierr == 0) then

            k = 0
            call save1('fneut_H', res(ion_ifneut_H))
            call save1('fneut_He', res(ion_ifneut_He))
            call save1('fneut_C', res(ion_ifneut_C))
            call save1('fneut_N', res(ion_ifneut_N))
            call save1('fneut_O', res(ion_ifneut_O))
            call save1('fneut_Ne', res(ion_ifneut_Ne))
            call save1('fneut_Mg', res(ion_ifneut_Mg))
            call save1('fneut_Si', res(ion_ifneut_Si))
            call save1('fneut_Fe', res(ion_ifneut_Fe))
            call save1('z_H', res(ion_iZ_H))
            call save1('z_He', res(ion_iZ_He))
            call save1('z_C', res(ion_iZ_C))
            call save1('z_N', res(ion_iZ_N))
            call save1('z_O', res(ion_iZ_O))
            call save1('z_Ne', res(ion_iZ_Ne))
            call save1('z_Mg', res(ion_iZ_Mg))
            call save1('z_Si', res(ion_iZ_Si))
            call save1('z_Fe', res(ion_iZ_Fe))
            call save1('logpp_H', res(ion_ilogpp_H))
            call save1('logpp_He', res(ion_ilogpp_He))
            call save1('logpp_C', res(ion_ilogpp_C))
            call save1('logpp_N', res(ion_ilogpp_N))
            call save1('logpp_O', res(ion_ilogpp_O))
            call save1('logpp_Ne', res(ion_ilogpp_Ne))
            call save1('logpp_Mg', res(ion_ilogpp_Mg))
            call save1('logpp_Si', res(ion_ilogpp_Si))
            call save1('logpp_Fe', res(ion_ilogpp_Fe))
         
         end if
         
         contains
         
         subroutine save1(str,v)
            character (len=*), intent(in) :: str
            real(dp), intent(in) :: v
            k = k+1; output_values(i,j,k) = v
         end subroutine save1

      end subroutine plot_one
      
      
      subroutine open_plot_files(io_first, io_last, io_params, io_rho, io_tmp, dir)
         integer, intent(IN) :: io_first, io_params, io_rho, io_tmp
         integer, intent(OUT) :: io_last
         character (len=256), intent(IN) :: dir
         character (len=256) :: fname
         integer :: io
         
         fname = trim(dir) // '/params.data'
         open(unit=io_params,file=trim(fname))
         
         fname = trim(dir) // '/logRho.data'
         open(unit=io_rho,file=trim(fname))
         
         fname = trim(dir) // '/logT.data'
         open(unit=io_tmp,file=trim(fname))
         
         io = io_first-1
         !call open1('logP')
         !call open1('logPgas')
         !call open1('logS')
         !call open1('chiRho')
         !call open1('chiT')
         !call open1('logCp')
         !call open1('logCv')
         !call open1('dlnE_dlnRho')
         !call open1('dlnS_dlnT')
         !call open1('dlnS_dlnRho')
         !call open1('mu')
         !call open1('log_free_e')
         !call open1('gamma1')
         !call open1('gamma3')
         !call open1('grad_ad')
         !call open1('eta')
         call open1('fneut_H')
         call open1('fneut_He')
         call open1('fneut_C')
         call open1('fneut_N')
         call open1('fneut_O')
         call open1('fneut_Ne')
         call open1('fneut_Mg')
         call open1('fneut_Si')
         call open1('fneut_Fe')
         call open1('z_H')
         call open1('z_He')
         call open1('z_C')
         call open1('z_N')
         call open1('z_O')
         call open1('z_Ne')
         call open1('z_Mg')
         call open1('z_Si')
         call open1('z_Fe')
         call open1('logpp_H')
         call open1('logpp_He')
         call open1('logpp_C')
         call open1('logpp_N')
         call open1('logpp_O')
         call open1('logpp_Ne')
         call open1('logpp_Mg')
         call open1('logpp_Si')
         call open1('logpp_Fe')
         !call open1('logE')
         !call open1('logW')
         io_last = io
         
         
         contains 
         
         
         subroutine open1(name)
            character (len=*), intent(in) :: name
            fname = trim(dir) // '/' // trim(name) // '.data'
            io = io+1; open(unit=io,file=trim(fname))
         end subroutine open1
      
      end subroutine open_plot_files


         
      
      

      end module ion_table_plot


