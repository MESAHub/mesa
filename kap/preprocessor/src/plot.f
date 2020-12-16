      program plot_kap_builder
      
      use kap_support

      implicit none
      
      call Build_Plots

      
      contains


      subroutine Build_Plots
         
         !logical, parameter :: co_enhanced = .true.
         logical, parameter :: co_enhanced = .false.



         integer, parameter :: io_unit0 = 40
         integer :: info_out, i, j, k
         character (len=256) :: dir
         
         double precision :: X, Zbase, dXC, dXO, Rho, logRho, T, logT, logK, logR
         double precision logT_min,logT_max,logR_min,logR_max,dlogT,dlogR
         integer logT_points, logR_points
         integer io,io_first,io_last,io_params,ioR,io_tmp,num_out

         double precision, allocatable :: output_values(:,:,:)
         
         !call read_namelist
         
         Zbase = 0.02d0
         X = 0.00d0
         dXC = 0.00d0
         dXO = 0.00d0
   
         call setup(co_enhanced,data_dir,type1_table,Zbase,X,dXC,dXO)
         
         dir = 'plot_data'
         write(*,*) 'write data for opacity plots to ' // trim(dir)
         
         if (co_enhanced) then
            write(*,*) 'co_enhanced'
         else
            write(*,*) 'not co_enhanced'
         end if
         
         if (.true.) then
            logT_points = 351
            logR_points = 351

            logT_min = 2.60
            logT_max = 3.8d0

            logR_min = -5d0
            logR_max = 1d0
            
         else
            logT_points = 2
            logR_points = 2

            logT_min = 7d0
            logT_max = 7.3d0

            logR_min = -1d0
            logR_max = 1d0
         end if

   !..open the output files
         io_params = io_unit0
         ioR = io_unit0+1
         io_tmp = io_unit0+2
         io_first = io_unit0+3

         call Open_Plot_Outfiles(io_first, io_last, io_params, ioR, io_tmp, dir)
         write(io_params, '(4f16.6,6x,2i10)') Zbase, X, dXC, dXO, logR_points, logT_points
         close(io_params)
         num_out = io_last - io_first + 1
         allocate(output_values(logR_points,logT_points,num_out))
         
   !..get the results

         dlogT = (logT_max - logT_min)/(logT_points-1)
         dlogR = (logR_max - logR_min)/(logR_points-1)
       
         do j=1, logT_points
            logT = logT_min + dlogT*(j-1)
            T = 10 ** logT

            do i=1,logR_points
               logR = logR_min + dlogR*(i-1)
               logRho = logR + 3*logT - 18
               Rho = 10**logRho
               
               !write(*,*) 'call Get_Results', i, j, logT, logR
               info_out = 0
               call Get_Results( &
                     Zbase, X, dXC, dXO, Rho, logRho, T, logT,  &
                     logK, co_enhanced, data_dir, type1_table, .false., info_out)
                           
               if (info_out == 0) then

                  k =   1; output_values(i,j,k) = logK
                              
               else
            
                  k =   1; output_values(i,j,k) = -1d99
               
               end if
            
            enddo
         
         enddo

 01   format(E30.22)
 
 
 
         ! write out the results
         do j=1,logT_points
            write(io_tmp,01) logT_min + dlogT*(j-1)
         end do
         close(io_tmp)

         do i=1,logR_points
            write(ioR,01) logR_min + dlogR*(i-1)
         enddo
         close(ioR)
         
         do k = 1, num_out
            write(*,*) k
            write(io_first+k-1,'(e14.6)') output_values(1:logR_points,1:logT_points,k)
         end do
      
         do io=io_first,io_last
            close(io)
         end do
   
      end subroutine Build_Plots


      subroutine Open_Plot_Outfiles(io_first, io_last, io_params, ioR, io_tmp, dir)
         integer, intent(IN) :: io_first, io_params, ioR, io_tmp
         integer, intent(OUT) :: io_last
         character (len=256), intent(IN) :: dir
         character (len=256) :: fname
         integer :: io
         
         fname = trim(dir) // '/params.data'
         open(unit=io_params,file=trim(fname))
         
         fname = trim(dir) // '/logR.data'
         open(unit=ioR,file=trim(fname))
         
         fname = trim(dir) // '/logT.data'
         open(unit=io_tmp,file=trim(fname))
         
         io = io_first
         fname = trim(dir) // '/logK.data'
         open(unit=io,file=trim(fname))
            
         io_last = io
      
      end subroutine Open_Plot_Outfiles


      end   

