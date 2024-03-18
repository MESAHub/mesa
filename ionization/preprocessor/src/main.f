!
! main
! francois hebert, jul 20 2010
!
! initialize values of parameters and call tfdh code for calculations
!



      program main
   
      use def_args
      use def_type
      use def_const
      use lib_util
      use mod_run_tfdh

      implicit none

      !
      ! the mode parameter allows you to chose what the program will do on running:
      !
      ! 1: run a single model and print radial profiles to file
      ! 2: scan across binary composition mix for a specified ne,T pair
      ! 3: scan across binding energy ('cut' parameter) for a specified ne,T, composition
      ! 4: scan across ne at specified values of T :: good for making plots
      ! 5: do a full ne,T table calculation :: good for lookup interpolation
      !
      integer, parameter :: mode = 1
   
      !======================================================================
      !======================================================================

      real (dp) :: ztr, temp, ne, rho
      real (dp) :: tmin, tmax, nmin, nmax
      real (dp) :: cut, cut_max, cut_min
      real (dp) :: frac, tanh_input, tanh_range
      real (dp) :: z_real, z_elem1, z_elem2
      real (dp) :: zbar_m, zbar_p, z_interp_m, z_interp_p

      real (dp), dimension(:), allocatable :: as, xs, zs, t_list
      integer :: i, j, ns, nn, nt, nc, nx, jn, jt, jc, jx
      type (datastruct) :: args
      logical :: return_profile

      integer, parameter :: f_unit = 1, g_unit = 2
      character(len=20) :: tempstring
      character(len=64) :: dir, file_z, file_r, file_n, file_t, file_rho
      character(len=64) :: file_profile, file_cut, file_mix

      !======================================================================
      !======================================================================

      !
      ! read in directory and assign file names for output
      !
      if (iargc() /= 1) call alert(1, '(main) tfdh takes 1 argument - an output directory')
      call getarg(1, dir)

      file_profile = trim(dir) // '/' // 'profile'
      file_z       = trim(dir) // '/' // 'z_table'
      file_r       = trim(dir) // '/' // 'r_table'
      file_t       = trim(dir) // '/' // 't_list'
      file_n       = trim(dir) // '/' // 'n_list'
      file_rho     = trim(dir) // '/' // 'r_list'
      file_cut     = trim(dir) // '/' // 'cut_scan'
      file_mix     = trim(dir) // '/' // 'mix_scan'

      !======================================================================
      !======================================================================

      !
      ! run program in chosen mode
      !
      ! in each case, the physical conditions for the calculation are set up
      ! these quantities are:
      !
      !> ne : electron density of plasma
      !> temp : temperature of plasma
      !> xs, as, zs : composition of plasma, given as a list of mass fractions,
      !               mass numbers, and atomic numbers, for all species present
      !
      !> ztr : the atomic number of the trace species for which we are trying
      !        to calculate the average ionization
      !
      ! then, the function run_tfdh is called. this will use input values to 
      ! initialize and perform the TFDH calculation
      !
      
      select case (mode)

         case (1)

            !
            ! set up conditions for calculation
            !
            ztr = 26_dp
            ne  = 1e28_dp
            temp = 1e7_dp

            ns = 2
            call reallocate_dp(xs, ns)
            call reallocate_dp(zs, ns)
            call reallocate_dp(as, ns)
            !xs = (/1.0_dp/)
            !zs = (/2.0_dp/)
            !as = (/4.0_dp/)
            xs = (/0.5_dp, 0.5_dp/)
            zs = (/1.0_dp, 2.0_dp/)
            as = (/1.0_dp, 4.0_dp/)

            return_profile = .true.
            
            write (*,*) 'initialization complete'

            !
            ! run code
            !
            write (*,*) 'running..'
            call run_tfdh(ztr, temp, ne, xs, zs, as, return_profile, args)

            !
            ! print summary output to stdout
            !
            write (*,*) 'ne   = ', args % neinf
            write (*,*) 'mu   = ', args % mu
            write (*,*) 'rho  = ', args % mu * args % neinf * mp
            write (*,*) 'dv0  = ', args % dv0
            write (*,*) 'znet = ', args % znet
            do i=1, ns
               write (*,*) 'rex_i = ', args % rex(i)
            end do
            write (*,*) 'rws  = ', args % rws
            write (*,*) 'rfinal  = ', args % rfinal

            !
            ! write profile to file
            ! NOTE: tempstring formatting is unable to support ns > 9
            !
            write (*,*) 'writing profiles to file..'

            open(unit = f_unit, file = file_profile)
            call write_header(ztr, temp, ne, xs, zs, as)

            write (f_unit, '(8i20)', advance='no') 0,1,2,3,4,5,6,7
            do i=1, ns
               write (f_unit, '(i20)', advance='no') 7+i
            end do
            write (f_unit, *) '' ! new line
            write (f_unit, '(8a20)', advance='no') 'r/a0', 'v(r)', 'z(r)', 'zb(r)', 'ne,tot', 'ne,b', 'ne,free', 'qi_tot'
            do i=1, ns
               write (tempstring, '(16x,a3,i1)') 'qi_', i
               write (f_unit, '(a20)', advance='no') tempstring
            end do
            write (f_unit, *) '' ! new line

            do j=1, args%output_length
               write (f_unit, '(8es20.9)', advance='no') args%r(j)/rbohr, args%v(j), args%z(j), args%zb(j), &
               args%ne(j), args%nb(j), args%ne(j)-args%nb(j), args%qitot(j)
               do i=1, ns
                  write (f_unit, '(es20.9)', advance='no') args%qi(i + ns*(j-1))
               end do
               write (f_unit, *) '' ! new line
            end do

            close(f_unit)
            write (*,*) 'done'

      !======================================================================
      !======================================================================

         case (2)

            !
            ! setup
            !
            ztr = 26_dp
            ne  = 1e28_dp
            temp = 1e7_dp

            ns = 2
            call reallocate_dp(xs, ns)
            call reallocate_dp(zs, ns)
            call reallocate_dp(as, ns)
            xs = (/0.5_dp, 0.5_dp/)
            zs = (/2.0_dp, 6.0_dp/)
            as = (/4.0_dp, 12.0_dp/)

            nx = 25
            tanh_range = 3_dp

            return_profile = .false.

            write (*,*) 'initialization complete'

            !
            ! run code, write to file as you go
            !
            if (ns /= 2) call alert(1, '(main) mode 2 only implement binary mixtures [for now]')

            open(unit=f_unit, file=file_mix, status='replace')
            write (f_unit, '(3a20)') 'Z_trace', 'log ne', 'log T'
            write (f_unit, '(3es20.9)') ztr, log10(ne), log10(temp)
            write (f_unit, *) '' 
            write (f_unit, '(3a20)') 'element #', 'Z', 'A'
            write (f_unit, '(i20,3es20.9)') 1, zs(1), as(1)
            write (f_unit, '(i20,3es20.9)') 2, zs(2), as(2)
            write (f_unit, *) ''
            write (f_unit, '(7i20)') 0,1,2,3,4,5,6
            write (f_unit, '(7a20)') 'xs(1)', 'xs(2)', 'z_net', 'zbar_m', 'z_interp_m', 'zbar_p', 'z_interp_p'

            !
            ! begin by calculating charge for pure OCP of elements 1, 2
            !
            write (*,*) 'running z_elem1..'
            xs(1) = 1.0_dp
            xs(2) = 0.0_dp
            call run_tfdh(ztr, temp, ne, xs, zs, as, return_profile, args)
            z_elem1 = args%znet

            write (*,*) 'running z_elem2..'
            xs(1) = 0.0_dp
            xs(2) = 1.0_dp
            call run_tfdh(ztr, temp, ne, xs, zs, as, return_profile, args)
            z_elem2 = args%znet


            write (f_unit, '(7es20.9)') 1.0_dp, 0.0_dp, z_elem1, zs(1), z_elem1, zs(1), z_elem1
            flush(f_unit)
            
            !
            ! then loop to vary mixture of elements 1,2 and compare real value
            ! to interpolation on the pure OCP cases calculated above
            !
            do, jx=1, nx
               write (*,'(a20,4x,i3,a1,i3)') 'running iteration..', jx, '/', nx

               ! use a tanh dependence for x1,x2
               tanh_input = (tanh_range/(nx-1.0_dp))*(2*jx-nx-1.0_dp)
               frac = (tanh(tanh_input) + 1.0_dp)/2.0_dp

               ! assign mass fractions to element 1+2
               ! start out with elem1, then progressively have more of elem2
               xs(1) = 1.0_dp - frac
               xs(2) = frac

               call run_tfdh(ztr, temp, ne, xs, zs, as, return_profile, args)
               z_real = args%znet

               ! zbar as defined over mass or particle number
               zbar_m = xs(1)*zs(1) + xs(2)*zs(2)
               zbar_p = as(1)*as(2) / (args%mu * (as(2)*xs(1) + as(1)*xs(2)))

               z_interp_m = z_elem1 + (z_elem2-z_elem1)/(zs(2)-zs(1))*(zbar_m-zs(1))
               z_interp_p = z_elem1 + (z_elem2-z_elem1)/(zs(2)-zs(1))*(zbar_p-zs(1))

               write (f_unit, '(7es20.9)') 1.0_dp-frac, frac, z_real, zbar_m, z_interp_m, zbar_p, z_interp_p
               flush(f_unit)
            end do

            write (f_unit, '(7es20.9)') 0.0_dp, 1.0_dp, z_elem2, zs(2), z_elem2, zs(2), z_elem2
            close(f_unit)
            write (*,*) 'done'

      !======================================================================
      !======================================================================

         case (3)

            !
            ! setup
            !
            ztr = 26_dp
            ne  = 1e28_dp
            temp = 1e7_dp

            ns = 1
            call reallocate_dp(xs, ns)
            call reallocate_dp(zs, ns)
            call reallocate_dp(as, ns)
            xs = (/1.0_dp/)
            zs = (/2.0_dp/)
            as = (/4.0_dp/)

            nc = 33
            cut_min = 1e-2_dp
            cut_max = 1e2_dp

            return_profile = .false.

            write (*,*) 'initialization complete'

            !
            ! run code, write to file as you go
            !
            open(unit=f_unit, file=file_cut, status='replace')
            call write_header(ztr, temp, ne, xs, zs, as)
            write (f_unit, '(3a20)') 'E_cut/kT', 'Z_cut', 'N_bound'
            flush(f_unit)
            
            do, jc=1, nc
               write (*,'(a20,4x,i3,a1,i3)') 'running iteration..', jc, '/', nc

               args % cut = - cut_min * (cut_max/cut_min)**((jc-1.0_dp)/(nc-1.0_dp))

               call run_tfdh(ztr, temp, ne, xs, zs, as, return_profile, args)

               write (f_unit, '(3es20.9)') args%cut, args%znet, args%ztr-args%znet
               flush(f_unit)
            end do

            close(f_unit)
            write (*,*) 'done'

      !======================================================================
      !======================================================================

         case (4)

            !
            ! setup
            !
            ztr = 26_dp

            nt = 1 !5
            call reallocate_dp(t_list, nt)
            t_list = (/1e7_dp/) !(/1e6_dp, 3e6_dp, 1e7_dp, 3e7_dp, 1e8_dp/)

            nn = 57
            nmin = 1e24_dp
            nmax = 1e31_dp

            ns = 1
            call reallocate_dp(xs, ns)
            call reallocate_dp(zs, ns)
            call reallocate_dp(as, ns)
            xs = (/1.0_dp/)
            zs = (/6.0_dp/)
            as = (/12.0_dp/)

            return_profile = .false.

            write (*,*) 'initialization complete'

            !
            ! run code, print output
            !
            open(unit=f_unit, file=file_t)
            do, jt=1, nt
               write (f_unit, '(es20.9)') t_list(jt)
            end do
            close(f_unit)

            open(unit=f_unit, file=file_n)
            do, jn=1, nn
               write (f_unit, '(es20.9)') nmin * (nmax/nmin)**((jn-1.0_dp)/(nn-1.0_dp))
            end do
            close(f_unit)

            open(unit=f_unit, file=file_rho)
            do, jn=1, nn
               write (f_unit, '(es20.9)') mp/sum(zs(:)*xs(:)/as(:)) * nmin * (nmax/nmin)**((jn-1.0_dp)/(nn-1.0_dp))
            end do
            close(f_unit)


            open(unit=f_unit, file=file_z, status='replace')
            open(unit=g_unit, file=file_r, status='replace')
            do, jt=1, nt
               do, jn=1, nn
                  write (*,'(a20,4x,i3,a1,i3,4x,a3,4x,i3,a1,i3)') 'running iteration..', jt, '/', nt, 'and', jn, '/', nn

                  temp = t_list(jt)
                  ne = nmin * (nmax/nmin)**((jn-1.0_dp)/(nn-1.0_dp))

                  call run_tfdh(ztr, temp, ne, xs, zs, as, return_profile, args)

                  write (f_unit, '(es20.9)', advance='no') args % znet 
                  write (g_unit, '(es20.9)', advance='no') args % rex 
               end do
               write (f_unit, *) ''; flush (f_unit)
               write (g_unit, *) ''; flush (g_unit)
            end do
            close(f_unit)
            close(g_unit)
            
            write (*,*) 'done'

      !======================================================================
      !======================================================================

         case (5)

            !
            ! setup
            !
            ztr = 26_dp

            nt = 3
            tmin = 1e6_dp
            tmax = 1e8_dp

            nn = 10
            nmin = 1e24_dp
            nmax = 1e30_dp

            ns = 1
            call reallocate_dp(xs, ns)
            call reallocate_dp(zs, ns)
            call reallocate_dp(as, ns)
            xs = (/1.0_dp/)
            zs = (/2.0_dp/)
            as = (/4.0_dp/)

            return_profile = .false.

            write (*,*) 'initialization complete'

            !
            ! run code, print output
            !
            open(unit=f_unit, file=file_t)
            do, jt=1, nt
               write (f_unit, '(es20.9)') tmin * (tmax/tmin)**((jt-1.0_dp)/(nt-1.0_dp))
            end do
            close(f_unit)

            open(unit=f_unit, file=file_n)
            do, jn=1, nn
               write (f_unit, '(es20.9)') nmin * (nmax/nmin)**((jn-1.0_dp)/(nn-1.0_dp))
            end do
            close(f_unit)

            open(unit=f_unit, file=file_rho)
            do, jn=1, nn
               write (f_unit, '(es20.9)') mp/sum(zs(:)*xs(:)/as(:)) * nmin * (nmax/nmin)**((jn-1.0_dp)/(nn-1.0_dp))
            end do
            close(f_unit)


            open(unit=f_unit, file=file_z, status='replace')
            do, jt=1, nt
               do, jn=1, nn
                  write (*,'(a20,4x,i3,a1,i3,4x,a3,4x,i3,a1,i3)') 'running iteration..', jt, '/', nt, 'and', jn, '/', nn

                  temp = tmin * (tmax/tmin)**((jt-1.0_dp)/(nt-1.0_dp))
                  ne  = nmin * (nmax/nmin)**((jn-1.0_dp)/(nn-1.0_dp))

                  call run_tfdh(ztr, temp, ne, xs, zs, as, return_profile, args)

                  write (f_unit, '(es20.9)', advance='no') args % znet 
               end do
               write (f_unit, *) ''; flush(f_unit)
            end do
            close(f_unit)

            !
            ! print t, n, rho list files
            !
            write (*,*) 'done'

      !======================================================================
      !======================================================================

         case default
            write (*,*) 'invalid mode parameter in main'

      end select

      ! fortran should automatically deallocate the arrays (?)

      !======================================================================
      !======================================================================

      contains

      subroutine write_header(ztr, temp, ne, xs, zs, as)

         real (dp), intent(in) :: ztr, temp, ne
         real (dp), intent(in), dimension(:) :: xs, zs, as

         integer :: ns, i

         ns = size(xs)
         if (ns /= size(as))  call alert(1, '(main: write_header) composition vectors xs,as of different length')
         if (ns /= size(zs))  call alert(1, '(main: write_header) composition vectors xs,zs of different length')

         write (f_unit, '(3a20)') 'Z_trace', 'log ne', 'log T'
         write (f_unit, '(3es20.9)') ztr, log10(ne), log10(temp)
         write (f_unit, *) '' 

         write (f_unit, '(4a20)') 'element #', 'X', 'Z', 'A'
         do i=1, ns
            write (f_unit, '(i20,3es20.9)') i, xs(i), zs(i), as(i)
         end do
         write (f_unit, *) ''
   
      end subroutine write_header
   
      end program main



   
   