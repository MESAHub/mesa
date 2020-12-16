
      program test_num
      
      use test_support
      use test_brent
      use test_newuoa
      use test_bobyqa
      use test_newton
      !use test_radau5_pollu, only: do_test_radau5_pollu
      !use test_radau5_hires, only: do_test_radau5_hires
      use test_int_support
      
      use test_beam
      use test_chemakzo
      use test_medakzo
      use test_vdpol
      
      use test_diffusion
      use test_simplex
            
      use const_def
      use const_lib
      use num_def
      use mtx_lib
      use mtx_def
      use utils_lib, only: mesa_error
      
      implicit none


      logical, parameter :: show_all = .false.  ! false for releases


      integer :: i, j, k, solver, decsol, omp_get_thread_num, ierr
      logical :: do_numerical_jacobian, m_band, j_band, quiet
      character (len=32) :: my_mesa_dir

      my_mesa_dir = '../..'         
      call const_init(my_mesa_dir,ierr)     
      if (ierr /= 0) then
         write(*,*) 'const_init failed'
         call mesa_error(__FILE__,__LINE__)
      end if        
      call math_init()

      quiet = .false.   
      m_band = .false.   
      j_band = .false.   
      do_numerical_jacobian = .false.   
      decsol = lapack

      ! newton solver
      do_numerical_jacobian = .false.
 
      write(*,*) 'call do_test_newton lapack'
      call do_test_newton(do_numerical_jacobian, lapack)
 
      write(*,*) 'call do_test_newton block_thomas_dble'
      call do_test_newton(do_numerical_jacobian, block_thomas_dble)

      write(*,*) 'call test_find0_quadratic'
      call test_find0_quadratic

      write(*,*) 'call test_find_max_quadratic'
      call test_find_max_quadratic

      write(*,*) 'call test_qsort'
      call test_qsort

      write(*,*) 'call do_test_newuoa'
      call do_test_newuoa

      write(*,*) 'call do_test_bobyqa'
      call do_test_bobyqa

      write(*,*) 'call do_test_simplex'
      call do_test_simplex

      write(*,*) 'call do_test_brent'
      call do_test_brent

      write(*,*) 'call test_binary_search'
      call test_binary_search

      write(*,*) 'call test_root routines'
      call test_root_with_brackets
      call test_root2
      call test_root3
      
      ! explicit solvers
      call test_cash_karp(show_all)
      call test_dopri(.false.,show_all)
      call test_dopri(.true.,show_all)


      ! test each implicit solver with dense matrix
      !   ijob      M     J           test
      !     1       I     F           vdpol
      write(*,*) 'ijob 1'
      decsol = lapack
      write(*,*) 'numerical jacobians'
      do_numerical_jacobian = .true.
      do i=1,num_solvers
         call do_test_vdpol(i,decsol,do_numerical_jacobian,show_all,quiet)
      end do    

      write(*,*) 'analytical jacobians'
      do_numerical_jacobian = .false.
      do i=1,num_solvers
         call do_test_vdpol(i,decsol,do_numerical_jacobian,show_all,quiet)
      end do
      
      ! test each implicit solver with banded matrix
      !   ijob      M     J           test
      !     2       I     B           medakzo
      write(*,*) 'ijob 2'
      decsol = lapack
      write(*,*) 'numerical jacobians'
      do_numerical_jacobian = .true.
      do i=1,num_solvers
         if (i <= ros3p_solver) cycle
         call do_test_medakzo(i,decsol,do_numerical_jacobian,show_all,quiet)
      end do
      
      
      write(*,*) 'analytical jacobians'
      do_numerical_jacobian = .false.
      do i=1,num_solvers
         if (i <= ros3p_solver) cycle
         call do_test_medakzo(i,decsol,do_numerical_jacobian,show_all,quiet)
      end do
      

! as of dec, 2013, non-identity mass matrix causes diff results with ifort vs gfortran
!      ! test each implicit solver with banded implicit ODE system and dense matrix
!      !   ijob      M     J           test
!      !     3       B     F           chemakzo
!      write(*,*) 'ijob 3'
!      decsol = lapack
!      do i=1,num_solvers
!         m_band = .true.
!         if (i <= ros3p_solver) cycle
!         call do_test_chemakzo(i,decsol,m_band,do_numerical_jacobian,show_all,quiet)
!      end do
!
!
!      ! each implicit solver with full implicit ODE system
!      !   ijob      M     J           test
!      !     5       F     F           chemakzo
!      write(*,*) 'ijob 5'
!      decsol = lapack
!      do i=1,num_solvers
!         m_band = .false.
!         if (i <= ros3p_solver) cycle
!         call do_test_chemakzo(i,decsol,m_band,do_numerical_jacobian,show_all,quiet)
!      end do
            
            
      ! test with m1 /= 0
      !   ijob      M     J           test
      !    11       I     F     x     beam
      write(*,*) 'ijob 11'
      decsol = lapack
      if (show_all) then
         do i=1,num_solvers
               if (i <= ros3pl_solver) cycle ! beam is too hard for these
            call do_test_beam(i,decsol,.true.,show_all,quiet)
         end do
      end if


      ! test solvers with tridiagonal jacobian
      decsol = lapack
      do_numerical_jacobian = .false.
      quiet = .false.
      do i=1,num_solvers
         call do_test_diffusion(i,decsol,do_numerical_jacobian,show_all,quiet)
      end do


      end program
