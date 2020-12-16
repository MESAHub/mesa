
      module test_newton
      use const_def, only: dp     
      use num_def
      use num_lib
      use mtx_def
      use mtx_lib
      use math_lib
      use utils_lib, only: mesa_error
      
      implicit none

      real(dp), parameter :: one=1
      
      integer, parameter :: nz = 1001, nvar = 2 !use odd number of zones for problem symmetry
      integer, parameter :: nsec = 0 ! number of secondaries per zone
      integer, parameter :: ldy = nz 

      integer, parameter :: i_conc=1, i_flux=2, equ_conc=1, equ_flux=2
      
      integer :: matrix_type
      real(dp), pointer, dimension(:) :: equ1, x1, xold1, dx1, xscale1, y1
      real(dp), pointer, dimension(:,:) :: equ, x, xold, dx, xscale, y
      real(dp), pointer, dimension(:,:,:) :: ublk, dblk, lblk

      logical, parameter :: dbg = .false.
      

      contains

      
      subroutine do_test_newton( &
         do_numerical_jacobian, which_decsol_in)
         logical, intent(in) :: do_numerical_jacobian
         integer, intent(in) :: which_decsol_in

         real(dp) :: dt, kappa, alphat, tmax, time, actual
         real(dp), pointer, dimension (:) :: concentration, fluxes
         real(dp) :: xmin, xmax, delx
         integer :: ierr, which_decsol, numsteps, midpt, maxsteps, neq
         character (len=64) :: decsol_option_name
         
         real(dp), parameter :: expected = 2.9347118120566711D-02 ! using lapack
         
         include 'formats.dek'

         which_decsol = which_decsol_in
         call decsol_option_str(which_decsol, decsol_option_name, ierr)
         if (ierr /= 0) return

         write(*,*) 'do_test_newton using ' // trim(decsol_option_name)

         ! diffusion problem setup
         kappa  = 1.0
         alphat = 1.d1 ! multiply explicit stability time step by this factor
         time   = 0.0
         xmin   = 0.0
         xmax   = 1.0
         delx   = (xmax - xmin)/float(nz) !use uniform spatial mesh
         tmax   = pow2(10.0*delx)/kappa !maximum evolution time in units of stability time step
      
         allocate(concentration(nz), fluxes(nz), stat=ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         
         neq = nvar*nz
         allocate( &
            equ1(neq), x1(neq), xold1(neq), dx1(neq), &
            xscale1(neq), y1(ldy*nsec), stat=ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         
         x(1:nvar,1:nz) => x1(1:neq)
         xold(1:nvar,1:nz) => xold1(1:neq)
         dx(1:nvar,1:nz) => dx1(1:neq)
         equ(1:nvar,1:nz) => equ1(1:neq)
         xscale(1:nvar,1:nz) => xscale1(1:neq)
         y(1:ldy,1:nsec) => y1(1:ldy*nsec)
         
         concentration = 0.0
         fluxes        = 0.0
         midpt = ceiling(float(nz)/2)
         concentration(midpt) = 1.0 !delta function spike 
         numsteps = 0
         maxsteps = 500
         dt = alphat*(delx*delx)/kappa ! explicit stability time step multiplied by alphat
         do while(time < tmax .and. numsteps < maxsteps)
            numsteps = numsteps + 1
            call do_1step_diffuse( &
               do_numerical_jacobian, which_decsol, &
               dt, kappa, nz, nvar, concentration, fluxes, delx, ierr)
            if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
            time = time + dt
            !write(*,2) 'diffusion step', numsteps, concentration(midpt), time/tmax
         end do
         
         actual = concentration(midpt)
         write(*,1) 'expected', expected
         write(*,1) 'actual', actual
         !write(*,1) '(actual - expected)/expected', (actual - expected)/expected
         write(*,*)
      
         deallocate(concentration, fluxes)
         deallocate(equ1, x1, xold1, dx1, xscale1, y1)
      
      end subroutine do_test_newton
      
      
      subroutine do_1step_diffuse( &
            do_numerical_jacobian, which_decsol, &
            dt, kappa, nz, nvar, concentration, fluxes, delx, ierr)
         integer i
         logical, intent(in) :: do_numerical_jacobian
         integer, intent(in) :: which_decsol
         integer, intent(in) :: nz, nvar
         real(dp), intent(inout) :: dt, kappa  
         real(dp), pointer, dimension(:), intent(inout) :: concentration, fluxes
         real(dp), intent(in) :: delx
         integer, intent(out) :: ierr

         integer :: liwork, lwork
         integer, dimension(:), pointer :: iwork
         real(dp), dimension(:), pointer :: work
         
         integer, parameter :: lrpar = 3, lipar = 1
         integer, target :: ipar_target(lipar)
         real(dp), target :: rpar_target(lrpar)
         integer, pointer :: ipar(:)
         real(dp), pointer :: rpar(:)
         
         integer :: lrd, lid
         integer, pointer :: ipar_decsol(:) ! (lid)
         real(dp), pointer :: rpar_decsol(:) ! (lrd)
         
         integer :: mljac, mujac
         real(dp) :: tol_correction_norm, tol_max_correction, tol_residual_norm
         logical :: nonconv
         
         include 'formats.dek'
         
         ierr = 0

         ipar => ipar_target
         rpar => rpar_target         

         rpar(1) = dt
         rpar(2) = kappa
         rpar(3) = delx
         
         if (do_numerical_jacobian) then
            ipar(1) = 1
         else
            ipar(1) = 0
         end if
         
         call set_fluxes(concentration, fluxes, kappa, delx)
         xold(i_conc,:) = concentration ! starting model
         xold(i_flux,:) = fluxes
         dx = 0d0
         x = xold
         
         tol_correction_norm = 1d-9 ! upper limit on magnitude of average scaled correction
         tol_max_correction = 1d99
         tol_residual_norm = 1d99

         mljac = 2*nvar-1 ! number of subdiagonals
         mujac = mljac ! number of superdiagonals
         if (which_decsol == lapack) then
            call lapack_work_sizes(nz*nvar, lrd, lid)
            if (do_numerical_jacobian) then
               matrix_type = square_matrix_type
               mljac = nz*nvar-1
               mujac = mljac
               if (nz*nvar > 51) then
                  write(*,*) 'numerical jac is very slow for large nz*nvar'
                  call mesa_error(__FILE__,__LINE__)
               end if
            else
               matrix_type = banded_matrix_type
            end if
         else
            write(*,*) 'bad value for which_decsol'
            call mesa_error(__FILE__,__LINE__)
         end if

         allocate(rpar_decsol(lrd), ipar_decsol(lid), stat=ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         
         call newton_work_sizes( &
            mljac, mujac, nvar, nz, nsec, matrix_type, lwork, liwork, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         
         allocate(work(lwork), iwork(liwork), stat=ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         
         work = 0
         iwork = 0
         
         iwork(i_try_really_hard) = 1 ! try really hard for first model
         iwork(i_model_number) = 1
                  
         !iwork(i_debug) = 1
         
         if (which_decsol == lapack) then
            call do_newt( &
               lapack_decsol, null_decsolblk, &
               nonconv, ierr)
         else
            stop 'bad which_decsol'
         end if
         
         if (nonconv) then
            write(*,*) 'failed to converge'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         concentration = x(i_conc,:)
         fluxes = x(i_flux,:)

         deallocate(iwork, work, rpar_decsol, ipar_decsol)
         
         
         contains         
         
         
         subroutine do_newt(decsol, decsolblk, nonconv, ierr)
            interface
               include "mtx_decsol.dek"
               include "mtx_decsolblk_dble.dek"
            end interface
            logical, intent(out) :: nonconv
            integer, intent(out) :: ierr
            
            real(dp), pointer :: AF1(:)
            AF1 => null()
            call newton( &
               nz, nvar, x1, xold1, matrix_type, mljac, mujac, &
               decsol, decsolblk, &
               lrd, rpar_decsol, lid, ipar_decsol, which_decsol, tol_correction_norm, &
               default_set_primaries, default_set_secondaries, diffusion_set_xscale, &
               default_Bdomain, default_xdomain, eval_equations, &
               default_size_equ, default_sizeB, default_inspectB, &
               enter_setmatrix, exit_setmatrix, failed_in_setmatrix, &
               default_force_another_iter, xscale1, equ1, ldy, nsec, y1, &
               work, lwork, iwork, liwork, &
               AF1, lrpar, rpar, lipar, ipar,  nonconv, ierr)
            if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
            if (associated(AF1)) deallocate(AF1)
         end subroutine do_newt
         
         
      end subroutine do_1step_diffuse
         
      
      subroutine set_fluxes(concentration, fluxes, kappa, delx)
         real(dp), intent(in) :: delx, kappa
         real(dp), pointer, dimension(:), intent(inout) :: concentration, fluxes
         integer :: i
         do i=2,nz
            fluxes(i) = -kappa*(concentration(i)-concentration(i-1))/delx
         end do      
         fluxes(1)  = 0
      end subroutine set_fluxes


      subroutine eval_equations(iter, nvar, nz, x, xscale, equ, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: iter, nvar, nz
         real(dp), pointer, dimension(:,:) :: x, xscale, equ ! (nvar, nz)
         integer, intent(in) :: lrpar, lipar
         real(dp), intent(inout) :: rpar(:) ! (lrpar)
         integer, intent(inout) :: ipar(:) ! (lipar)
         integer, intent(out) :: ierr

         integer :: k
         real(dp) :: dt, kappa, delx, fnzp1
         ierr = 0

         dt    = rpar(1)
         kappa = rpar(2)
         delx  = rpar(3)

         equ = 0.0

         ! Equation 1 : dc/dt = -div(fluxes)
         do k=1,nz-1
            equ(equ_conc,k) = &
               x(i_conc,k) - xold(i_conc,k) + dt*(x(i_flux,k+1) - x(i_flux,k))/delx
         end do
         fnzp1 = kappa*x(i_conc,nz)/delx
         equ(equ_conc,nz) = &
            x(i_conc,nz) - xold(i_conc,nz) + dt*(fnzp1 - x(i_flux,k))/delx

         ! Equation 2 : flux = -kappa*gradient(c)
         do k=2,nz
            equ(equ_flux,k) = x(i_flux,k) + kappa*(x(i_conc,k) - x(i_conc,k-1))/delx
         end do
         equ(equ_flux,1) = x(i_flux,1)

      end subroutine eval_equations


      subroutine eval_jacobian(ldA, A1, idiag, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: ldA ! leading dimension of A
         real(dp), pointer :: A1(:) ! (ldA, nvar*nz) ! the jacobian matrix
         ! A(idiag+q-v, v) = partial of equation(q) wrt variable(v)
         integer, intent(inout) :: idiag 
         integer, intent(in) :: lrpar, lipar
         real(dp), intent(inout) :: rpar(:) ! (lrpar)
         integer, intent(inout) :: ipar(:) ! (lipar)
         integer, intent(out) :: ierr         

         integer :: k
         real(dp) :: dt, kappa, delx
         real(dp), pointer :: blk3(:, :, :, :)
         ierr = 0
         
         blk3(1:nvar,1:nvar,1:nz,1:3) => A1(1:nvar*nvar*nz*3)
         ublk => blk3(:,:,:,1)
         dblk => blk3(:,:,:,2)
         lblk => blk3(:,:,:,3)

         dt = rpar(1)
         kappa = rpar(2)
         delx = rpar(3)

         A1 = 0.0
         do k=1,nz-1 ! concentration equ
            call e00(A1, equ_conc, i_conc, k, idiag, ldA, one)
            call e00(A1, equ_conc, i_flux, k, idiag, ldA, -dt/delx)
            call ep1(A1, equ_conc, i_flux, k, idiag, ldA, dt/delx)
         end do 
         call e00(A1, equ_conc, i_conc, nz, idiag, ldA, one + kappa*dt/delx**2)
            
         do k=2,nz ! flux equ
            call e00(A1, equ_flux, i_flux, k, idiag, ldA, one)
            call e00(A1, equ_flux, i_conc, k, idiag, ldA, kappa/delx)
            call em1(A1, equ_flux, i_conc, k, idiag, ldA, -kappa/delx)          
         end do 
         call e00(A1, equ_flux, i_flux, 1, idiag, ldA, one)

      end subroutine eval_jacobian
      
      
      subroutine e00(A1,i,j,k,idiag,ldA,v) ! partial of equ(i,k) wrt var(j,k)
         real(dp), pointer :: A1(:)
         real(dp) :: v
         real(dp), pointer :: A(:,:)
         integer, intent(in) :: i, j, k, idiag, ldA
         integer :: b, q, v00, sz, neq
         sz = size(A1,dim=1)
         neq = sz/ldA
         A(1:ldA,1:neq) => A1(1:sz)
         if (matrix_type == square_matrix_type) then
            b = nvar*(k-1)
            A(b+i,b+j) = A(b+i,b+j) + v
         else if (matrix_type == banded_matrix_type) then
            b = nvar*(k-1)
            q = idiag + b + i
            v00 = b + j
            A(q-v00,v00) = A(q-v00,v00) + v
         else if (matrix_type == block_tridiag_dble_matrix_type) then
            dblk(i,j,k) = dblk(i,j,k) + v
         else
            stop 'bad matrix_type'
         end if
      end subroutine e00
      
      
      subroutine em1(A1,i,j,k,idiag,ldA,v) ! partial of equ(i,k) wrt var(j,k-1)
         real(dp), pointer :: A1(:)
         real(dp) :: v
         real(dp), pointer :: A(:,:)
         integer, intent(in) :: i, j, k, idiag, ldA
         integer :: b, q, vm1, sz, neq
         sz = size(A1,dim=1)
         neq = sz/ldA
         A(1:ldA,1:neq) => A1(1:sz)
         if (k <= 0) stop 'bad k for em1'
         if (matrix_type == square_matrix_type) then
            b = nvar*(k-1)
            A(b+i,b+j-nvar) = A(b+i,b+j-nvar) + v
         else if (matrix_type == banded_matrix_type) then
            b = nvar*(k-1)
            q = idiag + b + i
            vm1 = b + j - nvar
            A(q-vm1,vm1) = A(q-vm1,vm1) + v
         else if (matrix_type == block_tridiag_dble_matrix_type) then
            lblk(i,j,k) = lblk(i,j,k) + v
         else
            stop 'bad matrix_type'
         end if
      end subroutine em1
      
      
      subroutine ep1(A1,i,j,k,idiag,ldA,v) ! partial of equ(i,k) wrt var(j,k+1)
         real(dp), pointer :: A1(:)
         real(dp) :: v
         real(dp), pointer :: A(:,:)
         integer, intent(in) :: i, j, k, idiag, ldA
         integer :: b, q, vp1, sz, neq
         sz = size(A1,dim=1)
         neq = sz/ldA
         A(1:ldA,1:neq) => A1(1:sz)
         if (k >= nz) stop 'bad k for ep1'
         if (matrix_type == square_matrix_type) then
            b = nvar*(k-1)
            A(b+i,b+j+nvar) = A(b+i,b+j+nvar) + v
         else if (matrix_type == banded_matrix_type) then
            b = nvar*(k-1)
            q = idiag + b + i
            vp1 = b + j + nvar
            A(q-vp1,vp1) = A(q-vp1,vp1) + v
         else if (matrix_type == block_tridiag_dble_matrix_type) then
            ublk(i,j,k) = ublk(i,j,k) + v
         else
            stop 'bad matrix_type'
         end if
      end subroutine ep1


      subroutine enter_setmatrix( &
            iter, nvar, nz, neqs, x, xold, xscale, xder, need_solver_to_eval_jacobian, &
            ldA, A1, idiag, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: iter, nvar, nz, neqs
         real(dp), pointer, dimension(:,:) :: x, xold, xscale, xder ! (nvar, nz)
         logical, intent(out) :: need_solver_to_eval_jacobian
         integer, intent(in) :: ldA ! leading dimension of A
         real(dp), pointer, dimension(:) :: A1 ! =(ldA, neqs)
         integer, intent(inout) :: idiag 
         integer, intent(in) :: lrpar, lipar
         real(dp), intent(inout) :: rpar(:) ! (lrpar)
         integer, intent(inout) :: ipar(:) ! (lipar)
         integer, intent(out) :: ierr
         real(dp) :: epsder
         include 'formats.dek'
         if (dbg) write(*, '(/, a)') 'enter_setmatrix'
         if (ipar(1) /= 0) then ! do numerical jacobian
            epsder = 1d-6 ! relative variation to compute numerical derivatives
            xder = epsder*(xscale+abs(xold))
            need_solver_to_eval_jacobian = .true.
         else
            call eval_jacobian(ldA, A1, idiag, lrpar, rpar, lipar, ipar, ierr)
            need_solver_to_eval_jacobian = .false.
         end if
      end subroutine enter_setmatrix


      subroutine exit_setmatrix( &
            iter, nvar, nz, neqs, dx, ldA, A1, idiag, xscale, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: ldA ! leading dimension of A
         integer, intent(in) :: iter, nvar, nz, neqs ! number of equations, 2nd dimension of A
         integer, intent(inout) :: idiag ! row of A with the matrix diagonal entries
         real(dp), pointer, dimension(:,:) :: dx
         real(dp), pointer, dimension(:) :: A1
         real(dp), pointer, dimension(:,:) :: xscale ! (nvar, nz)
         integer, intent(in) :: lrpar, lipar
         real(dp), intent(inout) :: rpar(:) ! (lrpar)
         integer, intent(inout) :: ipar(:) ! (lipar)
         integer, intent(out) :: ierr
         integer :: i, j
         if (dbg) write(*, '(a, /)') 'exit_setmatrix'
         ierr = 0
      end subroutine exit_setmatrix


      subroutine failed_in_setmatrix(j, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: j
         integer, intent(in) :: lrpar, lipar
         real(dp), intent(inout) :: rpar(:) ! (lrpar)
         integer, intent(inout) :: ipar(:) ! (lipar)
         integer, intent(out) :: ierr
         if (dbg) write(*, '(a, /)') 'failed_in_setmatrix'
         ierr = 0
      end subroutine failed_in_setmatrix


      ! you might want to use a different value of xscale_min for this
      subroutine diffusion_set_xscale(nvar, nz, xold, xscale, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: nvar, nz
         real(dp), pointer :: xold(:,:) ! (nvar, nz)
         real(dp), pointer :: xscale(:,:) ! (nvar, nz)
         integer, intent(in) :: lrpar, lipar
         real(dp), intent(inout) :: rpar(:) ! (lrpar)
         integer, intent(inout) :: ipar(:) ! (lipar)
         integer, intent(out) :: ierr
         real(dp), parameter :: xscale_min = 1d0
         xscale = 1.d0 ! max(xscale_min, abs(xold))
         ierr = 0
      end subroutine diffusion_set_xscale
            
      
      end module test_newton
