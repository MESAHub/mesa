      program test_colors

      use colors_def
      use colors_lib
      use math_lib
      use utils_lib, only: mesa_error, mkdir
      
      implicit none
      
      call do_test_colors

      contains


      subroutine do_test_colors
         use const_lib

         character (len=256) :: my_mesa_dir
         integer :: info,ierr
         
         logical, parameter :: do_one = .true.
         integer, parameter :: n_colors=11
   !     logical, parameter :: do_one = .false.
         
         my_mesa_dir = '../..'         
         call const_init(my_mesa_dir,info)     
         if (info /= 0) then
            write(*,*) 'const_init failed'
            call mesa_error(__FILE__,__LINE__)
         end if        
         call math_init()

         call colors_init(1,(/'../data/lcb98cor.dat'/),(/n_colors/),info)
         
         if (info /= 0) then
            write(*,*) 'colors_init failed during initialization'
            return
         end if

         if (do_one) then
            call do_one_colors
         else
             call create_plot_files
         end if
         call colors_shutdown()
               
      end subroutine do_test_colors 
      
      
      subroutine do_one_colors

         integer, parameter :: num_results=16
         real(dp) :: log_teff ! log10 of surface temp
         real(dp) :: log_l ! log10 of luminosity in solar units
         real(dp) :: mass ! mass in solar units
         real(dp) :: M_div_h ! [M/h], or as an approximation, log10[z/zsun]
         real(dp) :: boloMag
         real(dp), dimension(num_results) :: results
         real(dp) ::log_g
         integer :: info, i,total_num_colors
         character (len=8) :: vname
         character(len=strlen),dimension(num_results):: colors_name
         logical, parameter :: doing_solar = .true.
         real(dp), dimension(num_results) :: solar_expected_results = (/ &
         4.75d0, -0.11510d0, -0.14211d0, -0.61768d0, -0.36199d0, -0.68894d0, -1.46926d0, -0.32695d0, -0.78032d0, -0.39024d0, &
         0.05223d0, -0.10512d0,  -0.33801d0, -0.44312d0, -0.44123d0, -0.43080d0 /) 
         character(len=256),dimension(20) :: names

         
         ! solar values
         log_teff = log10(5780d0)
         log_l    = 0d0
         mass     = 1d0
         M_div_h     = 0d0
         log_g    = log10(mass) + 4.0D0*log_Teff - log_L - 10.6071D0 

         boloMag=get_abs_bolometric_mag(10**log_l) 
         
         !Store answers in results array
         
         ! These are bololmetric correction differences NOT magnitude differences
         ! Thus they are -1*mag diff
         results(1)=boloMag
         results(2)=get_bc_by_name('V',log_Teff,log_g, M_div_h, info)
         results(3)=get_bc_by_name('U',log_Teff,log_g, M_div_h, info)-get_bc_by_name('B',log_Teff,log_g, M_div_h, info)
         results(4)=get_bc_by_name('B',log_Teff,log_g, M_div_h, info)-get_bc_by_name('V',log_Teff,log_g, M_div_h, info)
         results(5)=get_bc_by_name('V',log_Teff,log_g, M_div_h, info)-get_bc_by_name('R',log_Teff,log_g, M_div_h, info)
         results(6)=get_bc_by_name('V',log_Teff,log_g, M_div_h, info)-get_bc_by_name('I',log_Teff,log_g, M_div_h, info)
         results(7)=get_bc_by_name('V',log_Teff,log_g, M_div_h, info)-get_bc_by_name('K',log_Teff,log_g, M_div_h, info)
         results(8)=get_bc_by_name('R',log_Teff,log_g, M_div_h, info)-get_bc_by_name('I',log_Teff,log_g, M_div_h, info)
         results(9)=get_bc_by_name('I',log_Teff,log_g, M_div_h, info)-get_bc_by_name('K',log_Teff,log_g, M_div_h, info)
         results(10)=get_bc_by_name('J',log_Teff,log_g, M_div_h, info)-get_bc_by_name('H',log_Teff,log_g, M_div_h, info)
         results(11)=get_bc_by_name('H',log_Teff,log_g, M_div_h, info)-get_bc_by_name('K',log_Teff,log_g, M_div_h, info)
         results(12)=get_bc_by_name('K',log_Teff,log_g, M_div_h, info)-get_bc_by_name('L',log_Teff,log_g, M_div_h, info)
         results(13)=get_bc_by_name('J',log_Teff,log_g, M_div_h, info)-get_bc_by_name('K',log_Teff,log_g, M_div_h, info)
         results(14)=get_bc_by_name('J',log_Teff,log_g, M_div_h, info)-get_bc_by_name('L',log_Teff,log_g, M_div_h, info)
         results(15)=get_bc_by_name('J',log_Teff,log_g, M_div_h, info)-get_bc_by_name('Lprime',log_Teff,log_g, M_div_h, info)
         results(16)=get_bc_by_name('K',log_Teff,log_g, M_div_h, info)-get_bc_by_name('M',log_Teff,log_g, M_div_h, info)
         
         if (info /= 0) then
            stop 'bad return from colors_get'
         end if
         
         write(*,*)
         write(*,*) 'color magnitude results'
         write(*,*)
         write(*,'(6a12)') 'teff', 'log_teff', 'log_l', 'mass', '[M_div_h]', 'log_g'
         write(*,'(i12,5f12.2)') floor(exp10(log_teff) + 0.5d0), log_teff, log_l, mass, M_div_h, log_g
         write(*,*)
         
         !call get_all_bc_names(colors_name,total_num_colors,info)
         colors_name=(/'bol ','bcv ','U-B ','B-V ','V-R ','V-I ','V-K ','R-I ',&
                        'I-K ','J-H ','H-K ','K-L ','J-K ','J-L ','J-Lp','K-M '/)
         
         do i=1,num_results
            write(*,'(9x,a8,f10.5)') colors_name(i), results(i)
         end do
         write(*,*)
         vname = 'vcolors'
         write(*,'(9x,a8,f10.5)') vname, results(1)-results(2)
         write(*,*)
         
         if (doing_solar) then
         
            do i=1,num_results
               if (abs(results(i) - solar_expected_results(i)) > 0.02d0) then
                  write(*,*)
                  write(*,*) "warning colors don't match solar expected values"
                  write(*,*) "Color: ",trim(colors_name(i))
                  write(*,*) "Got ", results(i), ' but expected ', solar_expected_results(i)
                  write(*,*)
                  stop
               end if
            end do
            
            write(*,*) 'matches expected solar results'
            write(*,*)
            
         end if         
         
      end subroutine do_one_colors


      subroutine create_plot_files

         character (len=256) fname, dir
         integer, parameter :: max_num_masses = 100
         integer, parameter :: num_results=16
         real(dp), dimension(max_num_masses) :: mass, logl, logteff
         real(dp) :: read_junk, log_g, M_div_h, boloMag
         real(dp),dimension(num_results) :: results 
         integer :: info, i, num_masses, io_unit, ios, iread_junk
         
         M_div_h = 0d0
         fname = 'zams_data/z02.log'
         io_unit = 40
         open(unit=io_unit, file=trim(fname), action='read', status='old', iostat=ios)
         if (ios /= 0) then
            write(*,*) 'failed to open the zams data'
            return
         end if
         
         num_masses = 0
         do i = 1, max_num_masses
            read(io_unit,fmt=*,iostat=ios) iread_junk, mass(i), logl(i), read_junk, logteff(i)
            if (ios /= 0) then
               num_masses = i-1
               exit
            end if
         end do
         read_junk = read_junk; iread_junk = iread_junk ! to keep g95 quiet
         
         close(io_unit)

         dir = 'plot_data'
         call mkdir(dir)
         fname = trim(dir) // '/' // 'colors.data'
         open(unit=io_unit,file=trim(fname))
            
         do i = 1, num_masses

            !call colors_get(logteff(i), logl(i), mass(i), M_div_h, results, log_g, info)
            
            boloMag=get_abs_bolometric_mag(10**logl(i)) 
            
            log_g=log10(mass(i)) + 4.0D0*logteff(i) - logl(i) - 10.6071D0  
            !Store answers in results array
            results(1)=boloMag
            results(2)=get_bc_by_name('V',logteff(i),log_g, M_div_h, info)
            results(3)=get_bc_by_name('U',logteff(i),log_g, M_div_h, info)-get_bc_by_name('B',logteff(i),log_g, M_div_h, info)
            results(4)=get_bc_by_name('B',logteff(i),log_g, M_div_h, info)-get_bc_by_name('V',logteff(i),log_g, M_div_h, info)
            results(5)=get_bc_by_name('V',logteff(i),log_g, M_div_h, info)-get_bc_by_name('R',logteff(i),log_g, M_div_h, info)
            results(6)=get_bc_by_name('V',logteff(i),log_g, M_div_h, info)-get_bc_by_name('I',logteff(i),log_g, M_div_h, info)
            results(7)=get_bc_by_name('V',logteff(i),log_g, M_div_h, info)-get_bc_by_name('K',logteff(i),log_g, M_div_h, info)
            results(8)=get_bc_by_name('R',logteff(i),log_g, M_div_h, info)-get_bc_by_name('I',logteff(i),log_g, M_div_h, info)
            results(9)=get_bc_by_name('I',logteff(i),log_g, M_div_h, info)-get_bc_by_name('K',logteff(i),log_g, M_div_h, info)
            results(10)=get_bc_by_name('J',logteff(i),log_g, M_div_h, info)-get_bc_by_name('H',logteff(i),log_g, M_div_h, info)
            results(11)=get_bc_by_name('H',logteff(i),log_g, M_div_h, info)-get_bc_by_name('K',logteff(i),log_g, M_div_h, info)
            results(12)=get_bc_by_name('K',logteff(i),log_g, M_div_h, info)-get_bc_by_name('L',logteff(i),log_g, M_div_h, info)
            results(13)=get_bc_by_name('J',logteff(i),log_g, M_div_h, info)-get_bc_by_name('K',logteff(i),log_g, M_div_h, info)
            results(14)=get_bc_by_name('J',logteff(i),log_g, M_div_h, info)-get_bc_by_name('L',logteff(i),log_g, M_div_h, info)
            results(15)=get_bc_by_name('J',logteff(i),log_g, M_div_h, info)-get_bc_by_name('Lprime',logteff(i),log_g, M_div_h, info)
            results(16)=get_bc_by_name('K',logteff(i),log_g, M_div_h, info)-get_bc_by_name('M',logteff(i),log_g, M_div_h, info)
            
            if (info == 0) write(io_unit,'(99f15.8)') mass(i), logl(i), logteff(i), log_g, results

         end do
         
         close(io_unit)
         
         write(*,*) 'finished creating plot files'

      end subroutine create_plot_files


      end program




