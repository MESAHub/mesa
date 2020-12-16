

      module test_bcyclic
      use cyclic_red, only: check_solver, bcyclic_init
      use mtx_lib, only: bcyclic_mt_lapack_decsolblk, bcyclic_mt_lapack_work_sizes
      use mtx_def, only: lapack
      use utils_lib, only: mesa_error
      implicit none

#ifdef mpi_opt
      include 'mpif.h'                                       !mpi stuff
#endif

      integer :: timeon, timeoff, countrate, ns0, nsn
      integer, parameter :: rprec = 8
      
      contains
      
      
      subroutine do_test_bcyclic(quiet, tol, fname_in, fname_out)
      logical, intent(in) :: quiet
      real*8, intent(in) :: tol
      character (len=*), intent(in) :: fname_in, fname_out
      
      integer                   :: ns, mblock, istat, k, nmin, irhs
      integer, pointer          :: ipivot(:,:)
      real*8, pointer, dimension(:,:,:) ::                     &
                   lblk, dblk, ublk, lblk1, dblk1, ublk1
      real*8, pointer  :: brhs(:,:), gc(:,:)
      real*8  :: t1, rms_error
      
      
      integer, parameter :: lid = 0, lrd = 0
      real*8, target :: rpar_decsol(lrd)
      integer, target :: ipar_decsol(lid)
      integer :: ierr, lid_check, lrd_check, iop
      
      include 'formats.dek'
      
!******************************************
      ierr = 0
      call read_data(ierr) ! sets mblock and ns; allocates and reads lblk, dblk, ublk
      if (ierr /= 0) return
!
!  mpi setup calls: 
!
      call bcyclic_init(ns,lapack,ierr)
      if (ierr /= 0) stop 'failed in bcyclic_init'

      write(*,3) 'test_bcyclic ' // trim(fname_in), mblock, ns

      call bcyclic_mt_lapack_work_sizes(mblock,ns,lrd_check,lid_check)
      if (lrd_check /= lrd .or. lid_check /= lid) then
         write(*,*) 'unexpected result from bcyclic_mt_lapack_work_sizes', lid_check, lrd_check
         call mesa_error(__FILE__,__LINE__)
      end if
!******************************************
!
!     do cyclic reduction
!     important: the blocks and data are assumed to be distributed consecutively
!                on nodes as described in the paper ....

!*****************************************

!     note: in production run, skip this since ipivot, lblk, dblk, ublk, brhs are
!           known inputs
     allocate (ipivot(mblock, ns0:nsn), stat=istat)
     if (istat .ne. 0) stop 'allocation error!'

!     store lblk1, etc only for back-solver check
      allocate (lblk1(mblock,mblock,ns0:nsn), dblk1(mblock,mblock,ns0:nsn),   &
                ublk1(mblock,mblock,ns0:nsn), gc(mblock,ns0:nsn),             &
                stat=istat)
      if (istat .ne. 0) stop 'allocation error!'

      lblk1 = lblk; dblk1 = dblk; ublk1 = ublk
      gc = brhs
!*****************************************

!     simulate multiple rhs
      multiple_rhs: do irhs = 1, 2

!     solve using cyclic reduction
!     solution overwrite brhs on each processor
      if (.not. quiet) then
#ifdef mpi_opt
      if (rank .eq. 0) print '(a,i1)',' solution rhs #',irhs
#else
      print '(a,i1)',' solution rhs #',irhs
#endif
      end if
      !call bcyclic_solver (lblk, dblk, ublk, ipivot, brhs, mblock, ns)
      iop = 0 ! factor
      call bcyclic_mt_lapack_decsolblk(iop,lblk,dblk,ublk,brhs,ipivot,lrd,rpar_decsol,lid,ipar_decsol,ierr)
      if (ierr /= 0) then
         write(*,*) 'failed in bcyclic_mt_lapack_decsolblk solve', ierr
         call mesa_error(__FILE__,__LINE__)
      end if
      iop = 1 ! solve
      call bcyclic_mt_lapack_decsolblk(iop,lblk,dblk,ublk,brhs,ipivot,lrd,rpar_decsol,lid,ipar_decsol,ierr)
      if (ierr /= 0) then
         write(*,*) 'failed in bcyclic_mt_lapack_decsolblk solve', ierr
         call mesa_error(__FILE__,__LINE__)
      end if

!     check solution
      call check_solver (lblk1, dblk1, ublk1, brhs, gc, mblock, ns, irhs, quiet, tol)

!     important: do not reset factored blocks
!     set brhs (solution)
      gc(:,ns0:nsn)=gc(:,ns0:nsn)*2 + 1.1d0 ! change from previous
      brhs(:,ns0:nsn)=gc(:,ns0:nsn)

      end do multiple_rhs

!     clean up allocated arrays
      !call clearstorage
      iop = 2 ! deallocate
      call bcyclic_mt_lapack_decsolblk(iop,lblk,dblk,ublk,brhs,ipivot,lrd,rpar_decsol,lid,ipar_decsol,ierr)
      if (ierr /= 0) then
         write(*,*) 'failed in bcyclic_mt_lapack_decsolblk deallocate', ierr
         call mesa_error(__FILE__,__LINE__)
      end if

      deallocate (ipivot, lblk, dblk, ublk, brhs, stat=istat)
      deallocate (lblk1, ublk1, dblk1, gc, stat=istat)

!     shut down mpi 
#ifdef mpi_opt
      call mpi_barrier(mpi_comm_world, ierr)
      call mpi_finalize(ierr)
#endif
      
      write(*,*) 'done test_bcyclic'
      write(*,*)
      
     
      contains
     
     
      subroutine read_data(ierr)
         integer, intent(out) :: ierr
         integer :: k, i, j, iounit
         ierr = 0
         iounit = 33
         write(*,*) 'read test matrix info ' // trim(fname_in)
         open(iounit,file=trim(fname_in), status='old', action='read', iostat=ierr)
         if (ierr /= 0) return
         read(iounit,*,iostat=ierr) mblock, ns
         if (ierr /= 0) return
         allocate (lblk(mblock,mblock,ns), dblk(mblock,mblock,ns), ublk(mblock,mblock,ns),brhs(mblock,ns))
         do k=1,ns
            do i=1,mblock
               read(iounit,*,iostat=ierr) brhs(i,k)
               if (ierr /= 0) return
               do j=1,mblock
                  read(iounit,*,iostat=ierr) lblk(i,j,k), dblk(i,j,k), ublk(i,j,k)
                  if (ierr /= 0) return
               end do
            end do
         end do
         close(iounit)
      end subroutine read_data
      

      end subroutine do_test_bcyclic



      end module test_bcyclic
