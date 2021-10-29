! ***********************************************************************
!
!   Copyright (C) 2016-2019  The MESA Team
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

      module net_burn_support
      use const_def, only: dp
      use math_lib
      use utils_lib, only: is_bad, mesa_error
      use mtx_def, only: lapack
      
      implicit none
      

      integer, parameter :: qcol_imax=13, kmaxx = 7, stifbs_imax  = kmaxx+1


      logical, parameter :: dbg = .false.



      contains
      
      
            
      subroutine netint( &
            start,stptry,stpmin,max_steps,stopp,y, &
            eps,species,nvar,nok,nbad,nstp,odescal,dens_dfdy,dmat, &
            derivs,jakob,burner_finish_substep,ierr)
      
! input:
! start    = beginning integration point
! stptry   = suggested first step size
! stpmin   = minimum allowable step size
! stopp    = ending integration point
! bc       = initial conditions, array of physical dimension yphys
! eps      = desired fraction error during the integration
! odescal  = error scaling factor
! derivs   = name of the routine that contains the odes
! jakob    = name of the routine that contains the jacobian of the odes
! steper   = name of the routine that will take a single step

! nvar >= species.  abundances come first in args,
! optionally followed by anything else such as lnT.

! output:
! nok      = number of succesful steps taken
! nbad     = number of bad steps taken, bad but retried and then succesful
! nstp     = total number of steps taken

      real(dp), pointer :: dens_dfdy(:,:),dmat(:,:)
      integer, intent(out) :: ierr
      
      interface
         include 'burner_derivs.inc'
      end interface         
      interface
         include 'burner_jakob.inc'
      end interface         
      interface
         include 'burner_finish_substep.inc'
      end interface

      integer    ::  species,nvar,nok,nbad,nstp,max_steps
      real(dp)   ::  start,stptry,stpmin,stopp,y(:),eps,odescal

      real(dp) :: yscal(nvar),dydx(nvar),cons,x,h,hdid,hnext,xx
      real(dp), parameter  :: tiny=1.0d-15

      
      real(dp) y0(nvar),a(stifbs_imax),alf(kmaxx,kmaxx),epsold,xnew,scale,red
      integer i,kmax,kopt,nseq(stifbs_imax),nvold
      logical first

      include 'formats'

! initialize
      first = .true.
      epsold = -1d0
      nvold = -1
      nseq = (/ 2, 6, 10, 14, 22, 34, 50, 70 /)
      x      = start
      h      = sign(stptry,stopp-start)
      nok    = 0
      nbad   = 0
      ierr = 0
      
      do i=1,nvar
         y0(i) = y(i)
      end do

! take at most max_steps steps
      do nstp=1,max_steps

! positive definite abundance fractions
       do i=1,species
         if (dbg) then
            if (is_bad(y(i))) then
               write(*,*) 'bad y for nstp', nstp, i, y(i)
               call mesa_error(__FILE__,__LINE__,'netint')
            end if
         end if
         y(i) = min(1.0d0, max(y(i),1.0d-30))
       end do
               
       call burner_finish_substep(nstp, x, y, ierr)
       if (ierr /= 0) return

! get the right hand sides and form scaling vector
       call derivs(x,y,dydx,nvar,ierr)
       if (ierr /= 0) then
         return
         write(*,*) 'derivs failed in netint'
       end if
       
      do i=1,nvar
         yscal(i) = max(odescal,abs(y(i)))
      end do

! if the step can overshoot the stop point cut it
       if ((x+h-stopp)*(x+h-start) .gt. 0.0d0) h = stopp - x

! do an integration step
        call stifbs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext, &
                    a,alf,epsold,first,kmax,kopt,nseq,nvold,xnew,scale,red, &
                    dens_dfdy,dmat,derivs,jakob,nstp,ierr)
        if (ierr /= 0) then
          write(*,*) 'stifbs ierr'
          return
        end if

! tally if we did or did not do the step
       if (hdid.eq.h) then
        nok = nok+1
       else
        nbad = nbad+1
       end if

! this is the normal exit point
       if ( (nstp .eq. max_steps) .or. &
            (x-stopp)*(stopp-start).ge. 0.0d0) then
        call burner_finish_substep(nstp, x, y, ierr)
        return
       end if

! normal timestep choice
       h = hnext

! die 
       if (abs(h).lt.stpmin) then
        write(*,*) 'netint failed: abs(h).lt.stpmin', abs(h), stpmin
        ierr = -1
        return

        write(6,210) 'nstp=',nstp
        write(6,210) 'nok nbad',nok,nbad
        write(6,220) 'attempted time step =',stptry
        write(6,220) 'current time step =',h
        write(6,220) 'input composition:'
        write(6,230) (y0(i), i=1,nvar)
        write(6,220) 'current composition:'
        write(6,230) (y(i), i=1,nvar)
        
        !call mesa_error(__FILE__,__LINE__,'h < stpmin in netint')

210     format(1x,a,4i6)
220     format(1x,a,1p3e24.16)
230     format(1x,1p3e24.16)

       end if


! back for another iteration or death
      enddo
      ierr = -1
      return
      call mesa_error(__FILE__,__LINE__,'more than max_steps steps required in netint')
      end subroutine netint


      

      subroutine stifbs(y,dydx,nvar,x,htry,eps,yscal,hdid,hnext, &
            a,alf,epsold,first,kmax,kopt,nseq,nvold,xnew,scale,red, &
            dens_dfdy,dmat,derivs,jakob,nstp,ierr)
!
! for dense analytic jacobians, lu decomposition linear algebra
!
! semi-implicit extrapolation step for integrating stiff ode's with monitoring
! of local truncation error to adjust stepsize. inputs are the dependent
! variable vector y(1:nvar) and its derivative dydx(1:nvar) at the starting of the
! independent variable x. also input are the stepsize to be attempted htry,
! the required accuracy eps, and the vector yscal against which the error is
! scaled. on output, y and x are replaced by their new values, hdid is the
! stepsize actually accomplished, and hnext is the estimated next stepsize.
! dervs is a user supplied function that computes the right hand side of
! the equations.
!
! declare
      real(dp) :: a(stifbs_imax),alf(kmaxx,kmaxx),epsold,xnew,scale,red
      integer :: kmax,kopt,nseq(stifbs_imax),nvold
      logical :: first
      integer :: ierr


      real(dp), pointer :: dens_dfdy(:,:),dmat(:,:)
      
      interface
         include 'burner_derivs.inc'
      end interface         
      interface
         include 'burner_jakob.inc'
      end interface         

      logical          reduct
      integer          nvar
      integer          i,iq,j,k,kk,km,i_errmax
      real(dp) y(:),dydx(:),x,htry,eps,yscal(:),hdid,hnext, &
                       eps1,errmax,fact,h,work,wrkmin,xest,err(kmaxx), &
                       yerr(nvar),ysav(nvar),yseq(nvar),safe1,safe2,dum, &
                       redmax,redmin,tiny,scalmx,qcol(nvar,qcol_imax),x_pzextr(qcol_imax)
      parameter        (safe1 = 0.25d0, safe2 = 0.7d0, redmax=1.0d-5, &
                        redmin = 0.7d0, tiny = 1.0d-30, &
                        scalmx = 0.1d0)

      integer          nstp, ierr2

      ierr = 0
      
! a new tolerance or a new number, so reinitialize
      if (eps .ne. epsold  .or.  nvar .ne. nvold) then
       hnext = -1.0d29
       xnew  = -1.0d29
       eps1  = safe1 * eps

! compute the work coefficients a_k
       a(1)  = nseq(1) + 1
       do k=1,kmaxx
        a(k+1) = a(k) + nseq(k+1)
       enddo

! compute alf(k,q)
       do iq=2,kmaxx
        do k=1,iq-1
         alf(k,iq) = pow(eps1,((a(k+1) - a(iq+1)) / &
                        ((a(iq+1) - a(1) + 1.0d0) * (2*k + 1))))
        enddo
       enddo
       epsold = eps
       nvold  = nvar

! add cost of jacobians to work coefficients
       a(1)   = nvar + a(1)
       do k=1,kmaxx
        a(k+1) = a(k) + nseq(k+1)
       enddo

! determine optimal row number for convergence
       do kopt=2,kmaxx-1
        if (a(kopt+1) .gt. a(kopt)*alf(kopt-1,kopt)) go to 01
       enddo
01     kmax = kopt
      end if

! save the starting values
      h    = htry
      do i=1,nvar
       ysav(i)  = y(i)
      enddo

! get the dense jacobian in dens_dfdy
      call jakob(x,y,dens_dfdy,nvar,ierr)
      if (ierr /= 0) then
         if (dbg) write(*,*) 'jakob failed in stifbs'
         return
       end if

! a new stepsize or a new integration, re-establish the order window
      if (h .ne. hnext  .or.  x .ne. xnew) then
       first = .true.
       kopt = kmax
      end if
      reduct = .false.

! evaluate the sequence of semi implicit midpoint rules
02    do 18 k=1,kmax

!       write(6,119) 'xnew x and h',xnew,x,h
! 119   format(1x,a,' ',1p3e12.4)

       xnew = x + h

!       write(6,119) 'xnew x and h',xnew,x,h
!       read(5,*)

       if (xnew .eq. x) then
         ierr = -1
         if (dbg) write(*,*) 'ierr: stepsize too small in routine stiffbs'
         return
         call mesa_error(__FILE__,__LINE__,'stepsize too small in routine stiffbs')
       end if

       call simpr( &
          ysav,dydx,nvar,x,h,nseq(k),yseq,dens_dfdy,dmat, &
          derivs,ierr)
       if (ierr /= 0) then


         h      = h * 0.1d0
         i_errmax   = 0
         reduct = .true.
         ierr = 0
         if (dbg) write(*,*) 'ierr: simpr failed in stifbs'
         go to 2


         write(*,*) 'simpr failed in stifbs'
         return
         
         
         
       end if
       xest = (h/nseq(k))*(h/nseq(k))
       call net_pzextr(k,xest,yseq,y,yerr,nvar,qcol,x_pzextr)

! compute normalized error estimate
       if (k .ne. 1) then
        errmax = tiny
        i_errmax   = 0
        do i=1,nvar
!        errmax = max(errmax,abs(yerr(i)/yscal(i)))
         dum = abs(yerr(i)/yscal(i))
         if (dum .ge. errmax) then
          errmax = dum
          i_errmax = i
         end if
        enddo

        errmax   = errmax/eps
        km = k - 1
        err(km) = pow(errmax/safe1,1.0d0/(2*km+1))
       end if

! in order window
       if (k .ne. 1  .and. (k .ge. kopt-1  .or. first)) then

! converged
        if (errmax .lt. 1.0d0) go to 04

! possible step size reductions
        if (k .eq. kmax  .or.  k .eq. kopt + 1) then
         red = safe2/err(km)
         go to 03
        else if (k .eq. kopt) then
         if (alf(kopt-1,kopt) .lt. err(km)) then
          red = 1.0d0/err(km)
          go to 03
         end if
        else if (kopt .eq. kmax) then
         if (alf(km,kmax-1) .lt. err(km)) then
          red = alf(km,kmax-1) * safe2/err(km)
          go to 03
         end if
        else if (alf(km,kopt) .lt. err(km)) then
         red = alf(km,kopt-1)/err(km)
         go to 03
        end if
       end if
18    continue

! reduce stepsize by at least redmin and at most redmax
03    red    = min(red,redmin)
      red    = max(red,redmax)
      h      = h * red
      i = i_errmax
      if (dbg) write(*,*) 'reduce stepsize', i, errmax, yerr(i), yscal(i), red, h
      i_errmax   = 0
      reduct = .true.
      go to 2


! successful step; get optimal row for convergence and corresponding stepsize
04    x = xnew
      hdid = h
      first = .false.
      wrkmin = 1.0d35
      do kk=1,km
       fact = max(err(kk),scalmx)
       work = fact * a(kk+1)
       if (work .lt. wrkmin) then
        scale  = fact
        wrkmin = work
        kopt   = kk + 1
       end if
      enddo
!
! check for possible order increase, but not if stepsize was just reduced
      hnext = h/scale
      if (kopt .ge. k  .and.  kopt .ne. kmax  .and.  .not.reduct) then
       fact = max(scale/alf(kopt-1,kopt),scalmx)
       if (a(kopt+1)*fact .le. wrkmin) then
        hnext = h/fact
        kopt = kopt + 1
       end if
      end if
      return
      end subroutine stifbs


      subroutine simpr( &
         y,dydx,nvar,xs,htot,nstep,yout,dens_dfdy,dmat, &
         derivs,ierr)
!
      real(dp), pointer :: dens_dfdy(:,:),dmat(:,:)
      
      interface
         include 'burner_derivs.inc'
      end interface
      
      integer, intent(out) :: ierr
      
      integer          nvar,nstep
      integer          i,j,nn,ii
      real(dp) y(:),dydx(:),xs,htot, &
                       yout(:),h,x,del(nvar),ytemp(nvar)

!..for the linear algebra
      integer, target :: indx_a(nvar)
      integer, pointer :: indx(:)
      
      include 'formats'

      indx => indx_a

! stepsize this trip, and make the a matrix
      h = htot/nstep
      do j=1,nvar
       do i=1,nvar
        dmat(i,j) = -h * dens_dfdy(i,j) 
       enddo
      enddo
      do i=1,nvar
       dmat(i,i) = 1.0d0 + dmat(i,i)
      end do

!..factor the matrix 
      call my_getf2(nvar, dmat, nvar, indx, ierr)  
      if (ierr /= 0) then
         if (dbg) write(*,*) 'my_getf2 failed in simpr'
         return
      end if    

! use yout as temporary storage; the first step
      do i=1,nvar 
         yout(i) = h * dydx(i)
         if (dbg) then
            if (is_bad(yout(i))) then
               write(*,*) 'bad yout in simpr nstep i yout', nstep, i, yout(i), dydx(i)
               call mesa_error(__FILE__,__LINE__,'simpr')
            end if
         end if
      enddo
      
      call my_getrs1(nvar, dmat, nvar, indx, yout, nvar, ierr)  
      if (ierr /= 0) then
         if (dbg) write(*,*) 'my_getrs1 failed in simpr'
         return
      end if     

      do i=1,nvar
         del(i)   = yout(i)
         ytemp(i) = y(i) + del(i)
         if (dbg) then
            if (is_bad(ytemp(i))) then
         
               do j=1,nvar
                do ii=1,nvar
                 if (dens_dfdy(ii,j) /= 0) write(*,3) 'dens_dfdy(ii,j)', ii, j, dens_dfdy(ii,j) 
                enddo
               enddo
            
               do ii=1,nvar
                 if (dydx(ii) /= 0) write(*,2) 'dydx(ii)', ii, dydx(ii)
               enddo
            
               do j=1,nvar
                do ii=1,nvar
                 if (dmat(ii,j) /= 0) write(*,3) 'dmat(ii,j)', ii, j, dmat(ii,j) 
                enddo
               enddo
            
               do ii=1,nvar
                 if (yout(ii) /= 0) write(*,2) 'yout(ii)', ii, yout(ii)
               enddo
            
               write(*,*) 'first step: bad ytemp in simpr nstep i ytemp', nstep, i, ytemp(i), del(i), y(i)
               call mesa_error(__FILE__,__LINE__,'simpr')
            
            end if
            
         end if
      enddo

      x = xs + h
      call derivs(x,ytemp,yout,nvar,ierr)
      if (ierr /= 0) then
         if (dbg) write(*,*) 'init derivs failed in simpr'
         return
       end if

! use yout as temporary storage; general step

      do nn=2,nstep
       do 15 i=1,nvar 
        yout(i) = h*yout(i) - del(i) 
         if (dbg) then
            if (is_bad(yout(i))) then
               write(*,*) 'bad yout in simpr nn i yout', nn, i, yout(i)
               call mesa_error(__FILE__,__LINE__,'simpr')
            end if
         end if
15     continue
      call my_getrs1(nvar, dmat, nvar, indx, yout, nvar, ierr)  
      if (ierr /= 0) then
         if (dbg) write(*,*) 'my_getrs1 failed in simpr'
         return
      end if     
       do i=1,nvar
        del(i)   = del(i) + 2.0d0 * yout(i)
        ytemp(i) = ytemp(i) + del(i)
         if (dbg) then
            if (is_bad(ytemp(i))) then
               write(*,*) 'general step: bad ytemp in simpr nn i yout', nn, i, ytemp(i)
               call mesa_error(__FILE__,__LINE__,'simpr')
            end if
         end if
       enddo

       x = x + h
       call derivs(x,ytemp,yout,nvar,ierr)
       if (ierr /= 0) then
         if (dbg) write(*,*) 'derivs failed in simpr general step'
         return
       end if
      enddo

! take the last step
      do 18 i=1,nvar 
       yout(i) = h * yout(i) - del(i)  
         if (dbg) then 
            if (is_bad(yout(i))) then
               write(*,*) 'bad yout in simpr last step: nstep i yout', nstep, i, yout(i)
               call mesa_error(__FILE__,__LINE__,'simpr')
            end if
         end if
18    continue
      call my_getrs1(nvar, dmat, nvar, indx, yout, nvar, ierr)  
      if (ierr /= 0) then
         write(*,*) 'my_getrs1 failed in simpr'
         return
      end if     
      
      do i=1,nvar
         yout(i) = ytemp(i) + yout(i)
         if (dbg) then
            if (is_bad(yout(i))) then
               write(*,*) 'bad yout in simpr result: nstep i yout', nstep, i, yout(i)
               call mesa_error(__FILE__,__LINE__,'simpr')
            end if
         end if
      end do

      return
      end subroutine simpr
      
      
      subroutine net_pzextr(iest,xest,yest,yz,dy,nvar,qcol,x)
! use polynomial extrapolation to evaluate nvar functions at x=0 by fitting
! a polynomial to a sequence of estimates with progressively smaller values
! x=xest, and corresponding function vectors yest(1:nvar). the call is number
! iest in the sequence of calls. extrapolated function values are output as
! yz(1:nvar), and their estimated error is output as dy(1:nvar)

! declare the pass
     integer :: iest,nvar
     real(dp)  :: xest,dy(:),yest(:),yz(:),qcol(:,:),x(:)

! locals; qcol and x must be "saved" between successive calls

      integer            :: j,k1
      real(dp)           :: delta,f1,f2,q,d(nvar)


! sanity checks

      if (iest .gt. qcol_imax) call mesa_error(__FILE__,__LINE__,'iest > qcol_imax in net_pzextr')

! save current independent variables
      x(iest) = xest
      do j=1,nvar
       dy(j) = yest(j)
       yz(j) = yest(j)
      enddo

! store first estimate in first column
      if (iest .eq. 1) then
       do j=1,nvar
        qcol(j,1) = yest(j)
       enddo
      else
       do j=1,nvar
        d(j) = yest(j)
       enddo
       do k1=1,iest-1
        delta = 1.0d0/(x(iest-k1) - xest)
        f1    = xest * delta
        f2    = x(iest-k1) * delta

! propagate tableu 1 diagonal more
        do j=1,nvar
         q          = qcol(j,k1)
         qcol(j,k1) = dy(j)
         delta      = d(j) - q
         dy(j)      = f1*delta
         d(j)       = f2*delta
         yz(j)      = yz(j) + dy(j)
        enddo
       enddo
       do j=1,nvar
        qcol(j,iest) = dy(j)
       enddo
      end if
      return
      end subroutine net_pzextr

      
      include 'mtx_solve_routines.inc'


      

      end module net_burn_support

