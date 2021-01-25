      module test_int_support
      use num_def
      use num_lib
      use utils_lib, only: mesa_error

      implicit none
            
      integer, parameter :: ipar_sparse_format = 1
         ! =0 means compressed row format; else, compressed column format.
      integer, parameter :: i_nfcn=2
      integer, parameter :: i_njac=3

      contains
      
      
      subroutine do_test_stiff_int( &
            which_solver, which_decsol, numerical_jacobian, &
            fcn, jac, sjac, solout, iout_input, &
            fcn_blk_dble,jac_blk_dble, &
            caller_id,nvar_blk_dble,nz_blk_dble,lblk,dblk,ublk,uf_lblk,uf_dblk,uf_ublk, &
            n,ndisc,mljac,mujac,matrix_type_spec, &
            mas,imas,mlmas,mumas,m1,m2,t, &
            rtol,atol,itol,h0,y,nstep,lrpar,rpar,lipar,ipar,quiet,ierr)
         use test_support,only:show_results,show_statistics
         use mtx_lib
         use mtx_def
         integer, intent(in) :: which_solver, which_decsol
         interface
            include 'num_fcn.dek'
            include 'num_fcn_blk_dble.dek'
            include 'num_jac.dek'
            include 'num_jac_blk_dble.dek'
            include 'num_sjac.dek'
            include 'num_solout.dek'
            include 'num_mas.dek'
         end interface
         integer :: caller_id, nvar_blk_dble, nz_blk_dble
         real(dp), dimension(:), pointer :: lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk ! =(nvar,nvar,nz)
         integer, intent(in) :: imas, mlmas, mumas, m1, m2, iout_input
         integer, intent(in) :: n, ndisc, mljac, mujac, matrix_type_spec, lrpar, lipar, itol
         logical, intent(in) :: numerical_jacobian, quiet
         real(dp), intent(inout) :: t(0:ndisc+1), rtol(*), atol(*), h0
         real(dp), pointer :: y(:) ! (n)
         integer, intent(inout) :: nstep
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr

         integer ::  i, k, nsteps, lout, iout, idid, ijac, max_cols_exptrap, &
            max_steps, xnstp, nfcn, njac, naccept, nreject, ndec, nsol, nrdens, &
            nzmax, lrd, lid, nrow, ncol, ndns, ndim, lfil, maxits, isparse, liwork, lwork
         real(dp) :: h, droptol, eps, max_step_size, y_min, y_max

         integer, pointer :: iwork(:) !(liwork)
         real(dp), pointer :: work(:) !(lwork)
         integer, pointer :: ipar_decsol(:) !(lid)
         real(dp), pointer :: rpar_decsol(:) !(lrd)
         
         iout = iout_input
         if (quiet) iout = 0
         max_steps = 500000
         max_step_size = 0 
         isparse = 0
         lout = 6
         
         y_min = -1d199
         y_max = 1d199
         
         if (numerical_jacobian) then
            ijac = 0
         else
            ijac = 1
         end if

         ipar = 0
         rpar = 0         

         nrdens = n
         max_cols_exptrap = 0 ! use default

         lid = 0; lrd = 0
         if (which_decsol == lapack) then
            nzmax = 0
            call lapack_work_sizes(n,lrd,lid)
         else
            if (mljac == n) then
               nzmax = n*n
            else
               nzmax = n*(mljac + mujac + 1)
            end if
            write(*,*) 'test_int_support: bad which_decsol', which_decsol
            call mesa_error(__FILE__,__LINE__) ! test_int_support
         end if

         call isolve_work_sizes(n,nzmax,imas,mljac,mujac,mlmas,mumas,liwork,lwork)
         
         allocate(iwork(liwork),work(lwork),ipar_decsol(lid),rpar_decsol(lrd),stat=ierr)
         if (ierr /= 0) then
            write(*,*) 'allocate ierr', ierr
            call mesa_error(__FILE__,__LINE__) ! test_int_support
         end if
      
         iwork = 0
         work = 0
      
         iwork(9) = m1
         iwork(10) = m2

         nstep = 0
         eps = rtol(1)
         do i=0,ndisc
            ierr = 0             
            h = h0
            select case(which_solver)
            case (ros2_solver)
               if (i==0 .and. .not. quiet) write(*,*) 'ros2'
            case (rose2_solver)
               if (i==0 .and. .not. quiet) write(*,*) 'rose2'
            case (ros3p_solver)
               if (i==0 .and. .not. quiet) write(*,*) 'ros3p'
            case (ros3pl_solver)
               if (i==0 .and. .not. quiet) write(*,*) 'ros3pl'
            case (rodas3_solver)
               if (i==0 .and. .not. quiet) write(*,*) 'rodas3'
            case (rodas4_solver)
               if (i==0 .and. .not. quiet) write(*,*) 'rodas4'
            case (rodasp_solver)
               if (i==0 .and. .not. quiet) write(*,*) 'rodasp'
            case default
               write(*,*) 'unknown value for which_solver'
               call mesa_error(__FILE__,__LINE__)
            end select

            if (which_decsol == lapack) then
               if (i==0 .and. .not. quiet) write(*,*) 'lapack_decsol'
               call do_isolve(lapack_decsol, null_decsols, null_decsolblk)
            else
               write(*,*) 'unknown value for which_decsol', which_decsol
               call mesa_error(__FILE__,__LINE__)
            end if

            if (idid /= 1) ierr = -1
            if (ierr /= 0) then
               write(*,*) 'solver returned ierr /= 0', idid
               call mesa_error(__FILE__,__LINE__)
            end if
            nstep = nstep + iwork(16) ! nsteps
         end do

         deallocate(iwork,work,ipar_decsol,rpar_decsol)
            
         contains
         
         
         subroutine do_isolve(decsol, decsols, decsolblk)
            interface
               include "mtx_decsol.dek"
               include "mtx_decsols.dek"
               include "mtx_decsolblk.dek"
            end interface
            integer :: j
            include 'formats'
            call isolve( &
               which_solver, n, fcn, t(i), y, t(i+1), & 
               h, max_step_size, max_steps, & 
               rtol, atol, itol, y_min, y_max, & 
               jac, ijac, sjac, nzmax, isparse, mljac, mujac, & 
               mas, imas, mlmas, mumas, & 
               solout, iout, & 
               decsol, decsols, decsolblk, &
               lrd, rpar_decsol, lid, ipar_decsol, &  
               caller_id, nvar_blk_dble, nz_blk_dble, &
               lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk, & 
               fcn_blk_dble, jac_blk_dble, &
               work, lwork, iwork, liwork, & 
               lrpar, rpar, lipar, ipar, & 
               lout, idid)
            return
            do j=1,n
               write(*,2) 'y(j)', j, y(j)
            end do
         end subroutine do_isolve
         

      end subroutine do_test_stiff_int


      end module test_int_support
