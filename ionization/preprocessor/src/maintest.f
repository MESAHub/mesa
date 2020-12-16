!
! maintest
! francois hebert, dec 26 2010
!
! generate simple file with summary info from a single tfdh calculation.. 
! with the purpose of comparing to a "reference" file to determine compilation success
!
! NOTE: use caution in modifying! you may no longer be able to meaningfully 
! compare with provided reference file
!


      program maintest
   
      use def_args
      use def_type
      use def_const
      use lib_util
      use mod_run_tfdh

      implicit none

      !======================================================================
      !======================================================================

      real (dp) :: ztr, tmp, ne, rho
      real (dp), dimension(:), allocatable :: as, xs, zs
      integer :: i, j, ns
      type (datastruct) :: args
      logical :: return_profile

      integer, parameter :: f_unit  = 1
      character(len=20) :: tmpstring
      character(len=64) :: dir, testfile

      !======================================================================
      !======================================================================

      !
      ! read in directory and assign file name for output
      !
      if (iargc() /= 1) call alert(1, '(main) tfdh takes 1 argument - an output directory')
      call getarg(1, dir)

      testfile = trim(dir) // '/' // 'test_file'
      
      !======================================================================
      !======================================================================

      !
      ! setup
      !
      ztr = 26_dp
      ne  = 1e27_dp
      tmp = 1e7_dp

      ns = 2
      call reallocate_dp(xs, ns)
      call reallocate_dp(zs, ns)
      call reallocate_dp(as, ns)
      xs = (/0.5_dp, 0.5_dp/)
      zs = (/1.0_dp, 2.0_dp/)
      as = (/1.0_dp, 4.0_dp/)

      return_profile = .false.

      write (*,*) 'initialization complete'

      !
      ! run code
      !
      write (*,*) 'running..'
      call run_tfdh(ztr, tmp, ne, xs, zs, as, return_profile, args)

      !
      ! write info to file
      ! NOTE: tmpstring formatting is unable to support ns > 9
      !
      write (*,*) 'writing summary to file..'

      open(unit = f_unit, file = testfile)
      call write_header(ztr, tmp, ne, xs, zs, as)
      write (f_unit,*) 'kt     = ', args % kt
      write (f_unit,*) 'ne     = ', args % neinf
      write (f_unit,*) 'mu     = ', args % mu
      write (f_unit,*) 'rho    = ', args % mu * args % neinf * mp
      write (f_unit,*) 'dv0    = ', args % dv0
      write (f_unit,*) 'znet   = ', args % znet
      write (f_unit,*) 'rex1   = ', args % rex(1)
      write (f_unit,*) 'rex2   = ', args % rex(2)
      write (f_unit,*) 'rws    = ', args % rws
      write (f_unit,*) 'nerws  = ', args % nerws
      close(f_unit)
      
      !======================================================================
      !======================================================================

      contains

      subroutine write_header(ztr, tmp, ne, xs, zs, as)

         real (dp), intent(in) :: ztr, tmp, ne
         real (dp), intent(in), dimension(:) :: xs, zs, as

         integer :: ns, i

         ns = size(xs)
         if (ns /= size(as))  call alert(1, '(main: write_header) composition vectors xs,as of different length')
         if (ns /= size(zs))  call alert(1, '(main: write_header) composition vectors xs,zs of different length')

         write (f_unit, '(3a20)') 'Z_trace', 'log ne', 'log T'
         write (f_unit, '(3es20.9)') ztr, log10(ne), log10(tmp)
         write (f_unit, *) '' 

         write (f_unit, '(4a20)') 'element #', 'X', 'Z', 'A'
         do i=1, ns
            write (f_unit, '(i20,3es20.9)') i, xs(i), zs(i), as(i)
         end do
         write (f_unit, *) ''
   
      end subroutine write_header
   
      end program maintest



   
   