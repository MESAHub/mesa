! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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

! Ernst Hairer's copyright for rodas can be found at the end of this file.


      module mod_rosenbrock
      use mod_dc_decsol
      use utils_lib
      use const_def, only: dp
      use math_lib
      
      
      logical, parameter :: dbg = .false.
      
      integer, parameter :: ns_max = 8 ! current max allowed value for number of stages
         ! okay to increase this if necessary.
            
      
      contains


      subroutine null_mas(n,am,lmas,lrpar,rpar,lipar,ipar)
         integer, intent(in) :: n, lmas, lrpar, lipar
         real(dp), intent(inout) :: am(lmas,n)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         am = 0
      end subroutine null_mas

      
      subroutine do_ros2(
     >      n,fcn,ifcn,x,y,xend,
     >      h,max_step_size,max_steps,
     >      rtol,atol,itol,y_min,y_max,
     >      jac,ijac,sjac,nzmax,isparse,
     >      mljac_in,mujac_in,dfx,idfx,
     >      mas,imas,mlmas,mumas,
     >      solout,iout,
     >      decsol, decsols, decsolblk, 
     >      lrd, rpar_decsol, lid, ipar_decsol,  
     >      caller_id, nvar, nz, lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk,
     >      fcn_blk_dble, jac_blk_dble,
     >      work,lwork,iwork,liwork,
     >      lrpar,rpar,lipar,ipar,
     >      lout,idid)
         implicit real(dp) (a-h,o-z)
#include "rodas_args.dek"
         integer, parameter :: ns = 2 ! number of stages
         call do_rodas(
     >      ns,contro3,ros2_coeffs,n,fcn,ifcn,x,y,xend,
     >      h,max_step_size,max_steps,
     >      rtol,atol,itol,y_min,y_max,
     >      jac,ijac,sjac,nzmax,isparse,
     >      mljac_in,mujac_in,dfx,idfx,
     >      mas,imas,mlmas,mumas,
     >      solout,iout,
     >      decsol, decsols, decsolblk, 
     >      lrd, rpar_decsol, lid, ipar_decsol,  
     >      caller_id, nvar, nz, lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk,  
     >      fcn_blk_dble, jac_blk_dble,
     >      work,lwork,iwork,liwork,
     >      lrpar,rpar,lipar,ipar,
     >      lout,idid)  
      end subroutine do_ros2

      
      subroutine do_rose2(
     >      n,fcn,ifcn,x,y,xend,
     >      h,max_step_size,max_steps,
     >      rtol,atol,itol,y_min,y_max,
     >      jac,ijac,sjac,nzmax,isparse,
     >      mljac_in,mujac_in,dfx,idfx,
     >      mas,imas,mlmas,mumas,
     >      solout,iout,
     >      decsol, decsols, decsolblk, 
     >      lrd, rpar_decsol, lid, ipar_decsol,  
     >      caller_id, nvar, nz, lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk,  
     >      fcn_blk_dble, jac_blk_dble,
     >      work,lwork,iwork,liwork,
     >      lrpar,rpar,lipar,ipar,
     >      lout,idid)
         implicit real(dp) (a-h,o-z)
#include "rodas_args.dek"
         integer, parameter :: ns = 3 ! number of stages
         call do_rodas(
     >      ns,contro3,rose2_coeffs,n,fcn,ifcn,x,y,xend,
     >      h,max_step_size,max_steps,
     >      rtol,atol,itol,y_min,y_max,
     >      jac,ijac,sjac,nzmax,isparse,
     >      mljac_in,mujac_in,dfx,idfx,
     >      mas,imas,mlmas,mumas,
     >      solout,iout,
     >      decsol, decsols, decsolblk, 
     >      lrd, rpar_decsol, lid, ipar_decsol,  
     >      caller_id, nvar, nz, lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk,  
     >      fcn_blk_dble, jac_blk_dble,
     >      work,lwork,iwork,liwork,
     >      lrpar,rpar,lipar,ipar,
     >      lout,idid)  
      end subroutine do_rose2

      
      subroutine do_ros3p(
     >      n,fcn,ifcn,x,y,xend,
     >      h,max_step_size,max_steps,
     >      rtol,atol,itol,y_min,y_max,
     >      jac,ijac,sjac,nzmax,isparse,
     >      mljac_in,mujac_in,dfx,idfx,
     >      mas,imas,mlmas,mumas,
     >      solout,iout,
     >      decsol, decsols, decsolblk, 
     >      lrd, rpar_decsol, lid, ipar_decsol,  
     >      caller_id, nvar, nz, lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk,  
     >      fcn_blk_dble, jac_blk_dble,
     >      work,lwork,iwork,liwork,
     >      lrpar,rpar,lipar,ipar,
     >      lout,idid)
         implicit real(dp) (a-h,o-z)
#include "rodas_args.dek"
         integer, parameter :: ns = 3 ! number of stages
         call do_rodas(
     >      ns,contro3,ros3p_coeffs,n,fcn,ifcn,x,y,xend,
     >      h,max_step_size,max_steps,
     >      rtol,atol,itol,y_min,y_max,
     >      jac,ijac,sjac,nzmax,isparse,
     >      mljac_in,mujac_in,dfx,idfx,
     >      mas,imas,mlmas,mumas,
     >      solout,iout,
     >      decsol, decsols, decsolblk, 
     >      lrd, rpar_decsol, lid, ipar_decsol,  
     >      caller_id, nvar, nz, lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk,  
     >      fcn_blk_dble, jac_blk_dble,
     >      work,lwork,iwork,liwork,
     >      lrpar,rpar,lipar,ipar,
     >      lout,idid)    
      end subroutine do_ros3p
      

      subroutine do_ros3pl(
     >      n,fcn,ifcn,x,y,xend,
     >      h,max_step_size,max_steps,
     >      rtol,atol,itol,y_min,y_max,
     >      jac,ijac,sjac,nzmax,isparse,
     >      mljac_in,mujac_in,dfx,idfx,
     >      mas,imas,mlmas,mumas,
     >      solout,iout,
     >      decsol, decsols, decsolblk, 
     >      lrd, rpar_decsol, lid, ipar_decsol,  
     >      caller_id, nvar, nz, lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk,  
     >      fcn_blk_dble, jac_blk_dble,
     >      work,lwork,iwork,liwork,
     >      lrpar,rpar,lipar,ipar,
     >      lout,idid)
         implicit real(dp) (a-h,o-z)
#include "rodas_args.dek"
         integer, parameter :: ns = 4 ! number of stages
         call do_rodas(
     >      ns,contro3,ros3pl_coeffs,n,fcn,ifcn,x,y,xend,
     >      h,max_step_size,max_steps,
     >      rtol,atol,itol,y_min,y_max,
     >      jac,ijac,sjac,nzmax,isparse,
     >      mljac_in,mujac_in,dfx,idfx,
     >      mas,imas,mlmas,mumas,
     >      solout,iout,
     >      decsol, decsols, decsolblk, 
     >      lrd, rpar_decsol, lid, ipar_decsol,  
     >      caller_id, nvar, nz, lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk,  
     >      fcn_blk_dble, jac_blk_dble,
     >      work,lwork,iwork,liwork,
     >      lrpar,rpar,lipar,ipar,
     >      lout,idid)     
      end subroutine do_ros3pl
      



      
      subroutine do_rodas3(
     >      n,fcn,ifcn,x,y,xend,
     >      h,max_step_size,max_steps,
     >      rtol,atol,itol,y_min,y_max,
     >      jac,ijac,sjac,nzmax,isparse,
     >      mljac_in,mujac_in,dfx,idfx,
     >      mas,imas,mlmas,mumas,
     >      solout,iout,
     >      decsol, decsols, decsolblk, 
     >      lrd, rpar_decsol, lid, ipar_decsol,  
     >      caller_id, nvar, nz, lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk,  
     >      fcn_blk_dble, jac_blk_dble,
     >      work,lwork,iwork,liwork,
     >      lrpar,rpar,lipar,ipar,
     >      lout,idid)
         implicit real(dp) (a-h,o-z)
#include "rodas_args.dek"
         integer, parameter :: ns = 4 ! number of stages
         call do_rodas(
     >      ns,contro3,rodas3_coeffs,n,fcn,ifcn,x,y,xend,
     >      h,max_step_size,max_steps,
     >      rtol,atol,itol,y_min,y_max,
     >      jac,ijac,sjac,nzmax,isparse,
     >      mljac_in,mujac_in,dfx,idfx,
     >      mas,imas,mlmas,mumas,
     >      solout,iout,
     >      decsol, decsols, decsolblk, 
     >      lrd, rpar_decsol, lid, ipar_decsol,  
     >      caller_id, nvar, nz, lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk,  
     >      fcn_blk_dble, jac_blk_dble,
     >      work,lwork,iwork,liwork,
     >      lrpar,rpar,lipar,ipar,
     >      lout,idid)     
      end subroutine do_rodas3


      subroutine do_rodas4(
     >      n,fcn,ifcn,x,y,xend,
     >      h,max_step_size,max_steps,
     >      rtol,atol,itol,y_min,y_max,
     >      jac,ijac,sjac,nzmax,isparse,
     >      mljac_in,mujac_in,dfx,idfx,
     >      mas,imas,mlmas,mumas,
     >      solout,iout,
     >      decsol, decsols, decsolblk, 
     >      lrd, rpar_decsol, lid, ipar_decsol,  
     >      caller_id, nvar, nz, lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk,  
     >      fcn_blk_dble, jac_blk_dble,
     >      work,lwork,iwork,liwork,
     >      lrpar,rpar,lipar,ipar,
     >      lout,idid)
         implicit real(dp) (a-h,o-z)
#include "rodas_args.dek"
         integer, parameter :: ns = 6 ! number of stages
         call do_rodas(
     >      ns,contro4,rodas4_coeffs,n,fcn,ifcn,x,y,xend,
     >      h,max_step_size,max_steps,
     >      rtol,atol,itol,y_min,y_max,
     >      jac,ijac,sjac,nzmax,isparse,
     >      mljac_in,mujac_in,dfx,idfx,
     >      mas,imas,mlmas,mumas,
     >      solout,iout,
     >      decsol, decsols, decsolblk, 
     >      lrd, rpar_decsol, lid, ipar_decsol,  
     >      caller_id, nvar, nz, lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk,  
     >      fcn_blk_dble, jac_blk_dble,
     >      work,lwork,iwork,liwork,
     >      lrpar,rpar,lipar,ipar,
     >      lout,idid)
      end subroutine do_rodas4


      subroutine do_rodasp(
     >      n,fcn,ifcn,x,y,xend,
     >      h,max_step_size,max_steps,
     >      rtol,atol,itol,y_min,y_max,
     >      jac,ijac,sjac,nzmax,isparse,
     >      mljac_in,mujac_in,dfx,idfx,
     >      mas,imas,mlmas,mumas,
     >      solout,iout,
     >      decsol, decsols, decsolblk, 
     >      lrd, rpar_decsol, lid, ipar_decsol,  
     >      caller_id, nvar, nz, lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk,  
     >      fcn_blk_dble, jac_blk_dble,
     >      work,lwork,iwork,liwork,
     >      lrpar,rpar,lipar,ipar,
     >      lout,idid)
         implicit real(dp) (a-h,o-z)
#include "rodas_args.dek"
         integer, parameter :: ns = 6 ! number of stages
         call do_rodas(
     >      ns,contro4,rodasp_coeffs,n,fcn,ifcn,x,y,xend,
     >      h,max_step_size,max_steps,
     >      rtol,atol,itol,y_min,y_max,
     >      jac,ijac,sjac,nzmax,isparse,
     >      mljac_in,mujac_in,dfx,idfx,
     >      mas,imas,mlmas,mumas,
     >      solout,iout,
     >      decsol, decsols, decsolblk, 
     >      lrd, rpar_decsol, lid, ipar_decsol,  
     >      caller_id, nvar, nz, lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk,  
     >      fcn_blk_dble, jac_blk_dble,
     >      work,lwork,iwork,liwork,
     >      lrpar,rpar,lipar,ipar,
     >      lout,idid)
      end subroutine do_rodasp

      
      subroutine do_rodas(
     >      ns,contro,coeffs,n,fcn,ifcn,x,y,xend,
     >      h,max_step_size,max_steps,
     >      rtol,atol,itol,y_min,y_max,
     >      jac,ijac,sjac,nzmax,isparse,
     >      mljac_in,mujac_in,dfx,idfx,
     >      mas,imas,mlmas,mumas,
     >      solout,iout,
     >      decsol, decsols, decsolblk, 
     >      lrd, rpar_decsol, lid, ipar_decsol,  
     >      caller_id, nvar, nz, lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk,  
     >      fcn_blk_dble, jac_blk_dble,
     >      work,lwork,iwork,liwork,
     >      lrpar,rpar,lipar,ipar,
     >      lout,idid)
         implicit real(dp) (a-h,o-z)
         integer, intent(in) :: ns ! number of stages
         interface
            real(dp) function contro(i,x,rwork,iwork,ierr)
               use const_def, only: dp
               integer, intent(in) :: i
               real(dp), intent(in) :: x
               real(dp), intent(inout), target :: rwork(*)
               integer, intent(inout), target :: iwork(*)
               integer, intent(out) :: ierr
            end function contro
            subroutine coeffs (ns,ra,rc,rd,ros_d,ros_m,ros_e,ros_alpha,
     &                ros_gamma,ros_newf,ros_elo,no_aux_in_error,ros_name)
               use const_def, only: dp
               integer, intent(in) :: ns
               real(dp), intent(inout) :: 
     >               ros_m(ns), ros_e(ns), ros_d(ns,ns), ra(ns,ns), rc(ns,ns), rd(ns,ns)
               real(dp), intent(inout) :: ros_alpha(ns), ros_gamma(ns)
               integer, intent(out) :: ros_elo
               logical, intent(out) :: no_aux_in_error, ros_newf(ns)
               character*12, intent(out) :: ros_name
            end subroutine coeffs
         end interface
#include "rodas_args.dek"

      logical autnms,implct,jband,arret,pred
      integer mljac, mujac, needed_lwork, needed_liwork
      real(dp), intent(in) :: y_min, y_max
      real(dp), pointer, dimension(:) :: p1, p2, p3, p4, p5
      integer, pointer, dimension(:) :: ip1
      mljac = mljac_in; mujac = mujac_in
! *** *** *** *** *** *** ***
!        setting the parameters 
! *** *** *** *** *** *** ***
      nfcn=0
      naccpt=0
      nrejct=0
      nstep=0
      njac=0
      ndec=0
      nsol=0
      arret=.false.
! -------- nmax , the maximal number of steps -----
      if(max_steps.eq.0)then
         nmax=100000
      else
         nmax=max_steps
         if(nmax.le.0)then
            if (lout > 0) write(lout,*)' wrong input max_steps=',max_steps
            arret=.true.
         end if
      end if
! -------- meth   coefficients of the method
      if(iwork(2).eq.0)then
         meth=1
      else
         meth=iwork(2)
         if(meth.le.0.or.meth.ge.4)then
            if (lout > 0) write(lout,*)' curious input iwork(2)=',iwork(2)
            arret=.true.
         end if
      end if
! -------- pred   step size control
      if(iwork(3).le.1)then
         pred=.true.
      else
         pred=.false.
      end if
! -------- parameter for second order equations
      m1=iwork(9)
      m2=iwork(10)
      nm1=n-m1
      if (m1.eq.0) m2=n
      if (m2.eq.0) m2=m1
      if (m1.lt.0.or.m2.lt.0.or.m1+m2.gt.n) then
       if (lout > 0) write(lout,*)' curious input for iwork(9,10)=',m1,m2
       arret=.true.
      end if
      nerror=iwork(11) ! number of variables to use for tolerances
      if (nerror.eq.0) nerror=n
! -------- uround   smallest number satisfying 1.d0+uround>1.d0  
      if(work(1).eq.0.d0)then
         uround=1.d-16
      else
         uround=work(1)
         if(uround.lt.1.d-16.or.uround.ge.1.d0)then
            if (lout > 0) write(lout,*)' coefficients have 16 digits, uround=',work(1)
            arret=.true.
         end if
      end if
! -------- maximal step size
      if(max_step_size.eq.0.d0)then
         hmax=xend-x
      else
         hmax=max_step_size
      end if
! -------  fac1,fac2     parameters for step size selection
      if(work(3).eq.0.d0)then
         fac1=5.d0
      else
         fac1=1.d0/work(3)
      end if
      if(work(4).eq.0.d0)then
         fac2=1.d0/3d0 ! 6.0d0
         ! originally default was 1/6, but that gave poor results on the diffusion test
      else
         fac2=1.d0/work(4)
      end if
      if (fac1.lt.1.0d0.or.fac2.gt.1.0d0) then
            if (lout > 0) write(lout,*)' curious input work(3,4)=',work(3),work(4)
            arret=.true.
         end if
! --------- safe     safety factor in step size prediction
      if (work(5).eq.0.0d0) then
         safe=0.9d0
      else
         safe=work(5)
         if (safe.le.0.001d0.or.safe.ge.1.0d0) then
            if (lout > 0) write(lout,*)' curious input for work(5)=',work(5)
            arret=.true.
         end if
      end if
! --------- check if tolerances are o.k.
      if (itol.eq.0) then
          if (atol(1).le.0.d0.or.rtol(1).le.10.d0*uround) then
              if (lout > 0) write(lout,*) ' tolerances are too small'
              arret=.true.
          end if
      else
          do i=1,n
             if (atol(i).le.0.d0.or.rtol(i).le.10.d0*uround) then
                if (lout > 0) write(lout,*) ' tolerances(',i,') are too small'
                arret=.true.
             end if
          end do
      end if
! *** *** *** *** *** *** *** *** *** *** *** *** ***
!         computation of array entries
! *** *** *** *** *** *** *** *** *** *** *** *** ***
! ---- autonomous, implicit, banded or not ?
      autnms=ifcn.eq.0
      implct=imas.ne.0
      jband=mljac.lt.nm1 .and. nzmax == 0
      if ((nzmax > 0) .and. (jband .or. ijac==0 .or. m1 /= 0)) then
         if (lout > 0) write(lout,*) 'sparse matrix -- nzmax > 0 -- requires ijac=1, m1=0, mljac=n'
          arret=.true.
      end if
! -------- computation of the row-dimensions of the 2-arrays ---
! -- jacobian and matrix e
      if(jband)then
         ldjac=mljac+mujac+1
         lde=mljac+ldjac
      else
         mljac=nm1
         mujac=nm1
         ldjac=nm1
         lde=nm1
      end if
! -- mass matrix
      if (implct) then
          if (mlmas.ne.nm1) then
              ldmas=mlmas+mumas+1
              if (nzmax > 0) then ! sparse jacobian
                 ijob=9
              else if (jband) then
                 ijob=4
              else
                 ijob=3
              end if
          else
              ldmas=nm1
              ijob=5
          end if
! ------ bandwith of "mas" not larger than bandwith of "jac"
          if (mlmas.gt.mljac.or.mumas.gt.mujac) then
              if (lout > 0) then
                  write(lout,*) 'bandwith of "mas" must not be larger than bandwith of "jac"'
                  write(lout,*) 'mlmas', mlmas
                  write(lout,*) 'mljac', mljac
                  write(lout,*) 'mumas', mumas 
                  write(lout,*) 'mujac', mujac
              end if
              arret=.true.
          end if
      else
          ldmas=0
          if (nzmax > 0) then ! sparse jacobian
             ijob=8
          else if (jband) then
             ijob=2
          else
             ijob=1
          end if
      end if
      ldmas2=max(1,ldmas)

      call calculate_work_sizes(
     >      n, ns_max, ldjac, nm1, ldmas, lde, nzmax,
     >      needed_lwork, needed_liwork, ieynew, iedy1, iedy, ieak, iefx, iecon, 
     >      iejac, iemas, iee, iesj, iesa, ieip, ieia, ieja)
     
      if(needed_lwork.gt.lwork)then
         ierr = 0
         call realloc_double(work,needed_lwork,ierr)
         if (ierr /= 0) then
            write(lout,*)
     >         ' insufficient storage for work, min. lwork=',needed_lwork
            arret=.true.
         end if
      end if
      
      if(needed_liwork.gt.liwork)then
         ierr = 0
         call realloc_integer(iwork,needed_liwork,ierr)
         if (ierr /= 0) then
            write(lout,*)
     >         ' insufficient storage for iwork, min. liwork=',needed_liwork
            arret=.true.
         end if
      end if
      
! ------ when a fail has occured, we return with idid=-1
      if (arret) then
         idid=-1
         return
      end if
      
! -------- call to core integrator ------------
      p1(1:n*ns) => work(ieak:ieak+n*ns-1)   ! arg ak1
      p2(1:lde*nm1) => work(iee:iee+lde*nm1-1) ! arg e1 
      p3(1:n) => work(ieynew:ieynew+n-1) ! ynew
      p4(1:n) => work(iedy1:iedy1+n-1) ! dy1
      p5(1:n) => work(iedy:iedy+n-1) ! dy
      
      ip1(1:nm1) => iwork(ieip:ieip+nm1-1)
      call roscor(
     &   ns,contro,coeffs,n,fcn,x,y,xend,hmax,h,rtol,atol,itol,y_min,y_max,
     &   jac,ijac,sjac,nzmax,isparse,
     &   mljac,mujac,dfx,idfx,mas,mlmas,mumas,solout,iout,idid,nmax,
     &   uround,meth,ijob,fac1,fac2,safe,autnms,implct,jband,pred,ldjac,
     &   lde,ldmas2,p3,p4,p5,p1,
     &   work(iefx:lwork),work(iejac:lwork),p2,work(iemas:lwork),
     &   ip1,work(iecon:lwork),
     &   decsol,decsols,decsolblk,
     &   caller_id, nvar, nz, lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk,
     &   fcn_blk_dble, jac_blk_dble,
     &   iwork(ieia:liwork),iwork(ieja:liwork),
     &   work(iesj:lwork),work(iesa:lwork),lrd,rpar_decsol,lid,ipar_decsol,m1,m2,nm1,nerror,
     &   nfcn,njac,nstep,naccpt,nrejct,ndec,nsol,lout,lrpar,rpar,lipar,ipar)
      iwork(14)=nfcn
      iwork(15)=njac
      iwork(16)=nstep
      iwork(17)=naccpt
      iwork(18)=nrejct
      iwork(19)=ndec
      iwork(20)=nsol
! ----------- return -----------
      return
      end subroutine do_rodas
!




      subroutine calculate_work_sizes(
     >      n, ns_max, ldjac, nm1, ldmas, lde, nzmax,
     >      lwork, liwork, ieynew, iedy1, iedy, ieak, iefx, iecon, 
     >      iejac, iemas, iee, iesj, iesa, ieip, ieia, ieja)
         integer, intent(in) :: n, ns_max, ldjac, nm1, ldmas, lde, nzmax
         integer, intent(out) :: 
     >      lwork, liwork, ieynew, iedy1, iedy, ieak, iefx, iecon, 
     >      iejac, iemas, iee, iesj, iesa, ieip, ieia, ieja
         ieynew=21
         iedy1=ieynew+n
         iedy=iedy1+n
         ieak=iedy+n      
         iefx=ieak+n*ns_max
         iecon=iefx+n
         iejac=iecon+2+5*n
         iemas=iejac+n*ldjac
         iee=iemas+nm1*ldmas
         iesj=iee+nm1*lde
         iesa=iesj+nzmax
         lwork=iesa+nzmax-1  
         ieip=21
         ieia=ieip+nm1
         ieja=ieia+n+1
         liwork=ieja+nzmax-1
      end subroutine calculate_work_sizes






!
!
!  ----- ... and here is the core integrator  ----------
!
      subroutine roscor(
     &  ns,contro,coeffs,n,fcn,x,y,xend,hmax,h,rtol,atol,itol,y_min,y_max,
     &  jac,ijac,sjac,nzmax,isparse,mljac,mujac,dfx,idfx,mas,mlmas,mumas,
     &  solout,iout,idid,nmax,uround,meth,ijob,fac1,fac2,safe,autnms,implct,banded,
     &  pred,ldjac,lde,ldmas,ynew,dy1,dy,ak1,fx,fjac,e1,fmas,ip,rwork,
     &  decsol,decsols,decsolblk,
     &  caller_id, nvar, nz, lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk,
     &  fcn_blk_dble, jac_blk_dble,
     &  ia,ja,sparse_jac,sa,
     &  lrd,rpar_decsol,lid,ipar_decsol,m1,m2,nm1,nerror,
     &  nfcn,njac,nstep,naccpt,nrejct,ndec,nsol,lout,lrpar,rpar,lipar,ipar)
! ----------------------------------------------------------
!     core integrator for rodas4
!     parameters same as in rodas4 with workspace added 
! ---------------------------------------------------------- 
!         declarations 
! ---------------------------------------------------------- 
      implicit real(dp) (a-h,o-z)
       interface
         real(dp) function contro(i,x,rwork,iwork,ierr)
            use const_def, only: dp
            integer, intent(in) :: i
            real(dp), intent(in) :: x
            real(dp), intent(inout), target :: rwork(*)
            integer, intent(inout), target :: iwork(*)
            integer, intent(out) :: ierr
         end function contro
         subroutine coeffs (ns,ra,rc,rd,ros_d,ros_m,ros_e,ros_alpha,
     &                ros_gamma,ros_newf,ros_elo,no_aux_in_error,ros_name)
            use const_def, only: dp
            implicit none
            integer, intent(in) :: ns
            real(dp), intent(inout) :: ros_m(ns), ros_e(ns), ros_d(ns,ns), ra(ns,ns), rc(ns,ns), rd(ns,ns)
            real(dp), intent(inout) :: ros_alpha(ns), ros_gamma(ns)
            integer, intent(out) :: ros_elo
            logical, intent(out) :: no_aux_in_error, ros_newf(ns)
            character*12, intent(out) :: ros_name
         end subroutine coeffs
#include "num_solout.dek"
#include "num_mas.dek"
#include "num_dfx.dek"
#include "num_fcn.dek"
#include "num_fcn_blk_dble.dek"
#include "num_jac.dek"
! for double block tridiagonal matrix
#include "num_jac_blk_dble.dek" 
#include "num_sjac.dek"
#include "mtx_decsol.dek"
#include "mtx_decsols.dek"
#include "mtx_decsolblk.dek"
       end interface
      integer, intent(in) :: nzmax, lrpar, lipar, lrd, lid
      integer, intent(out) :: ia(n+1), ja(nzmax)
      real(dp), intent(inout) :: sparse_jac(nzmax), sa(nzmax)
      integer, intent(in) :: caller_id, nvar, nz
      real(dp), dimension(:), pointer, intent(inout) :: lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk ! =(nvar,nvar,nz)
      ! the rosenbrock method parameters      
      integer, intent(in) :: ns ! number of stages
      real(dp) :: ros_m(ns), ros_e(ns) 
      real(dp) :: ros_alpha(ns), ros_gamma(ns)
      integer :: ros_elo
      logical :: ros_newf(ns)
      character*12 :: ros_name
      
      ! args
      integer, intent(inout), pointer :: ipar(:) ! (lipar)
      real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
      real(dp), intent(inout), pointer :: y(:) ! (n)
      real(dp), intent(inout), pointer :: dy(:) ! (n)
      real(dp), intent(inout), pointer :: ynew(:) ! (n)
      real(dp), intent(inout), pointer :: dy1(:) ! (n)
      
      dimension fx(n),fjac(ldjac,n),fmas(ldmas,nm1),atol(*),rtol(*)
     
      integer iwork(2)
      real(dp), target :: rwork(2+5*n)
      real(dp), intent(inout), pointer :: rpar_decsol(:) ! (lrd)
      integer, intent(inout), pointer :: ipar_decsol(:) ! (lid)
      real(dp), intent(in) :: y_min, y_max
      real(dp), pointer :: ak1(:) ! =(n,ns)
      real(dp), pointer :: e1(:) ! =(lde,nm1)
      integer, pointer :: ip(:) ! (nm1)
      
      ! locals
      dimension hd(ns), ra(ns,ns), rhc(ns,ns), rc(ns,ns), rd(ns,ns), ros_d(ns,ns)
      logical reject,autnms,implct,banded,last,pred,not_stage1,no_aux_in_error,need_free
      real(dp), pointer :: cont(:), dy2(:)
      real(dp) :: hprev, rd32
      integer :: i, j
      
      real(dp), pointer :: ak(:,:), e(:,:), p1(:)
      ak(1:n,1:ns) => ak1(1:n*ns)
      e(1:lde,1:nm1) => e1(1:lde*nm1)

! *** *** *** *** *** *** ***
!  initialisations
! *** *** *** *** *** *** ***    
      rc = 0
      md = 0 !n/2-1 ! for dbg
      need_free = .false.
      cont => rwork(1:5*n)
      dy2 => cont(n+1:2*n)
      nn=n 
      nn2=2*n
      nn3=3*n
      nn4=4*n
      lrc=5*n
! ------- compute mass matrix for implicit case ----------
      if (implct) call mas (nm1,fmas,ldmas,lrpar,rpar,lipar,ipar)
! ------ set the parameters of the method -----
      call coeffs(ns,ra,rc,rd,ros_d,ros_m,ros_e,ros_alpha,
     &                     ros_gamma,ros_newf,ros_elo,no_aux_in_error,ros_name)
      gamma = ros_gamma(1)
      if (ns >=3) then
         rd32 = rd(3,2)
      else
         rd32 = 0
      end if
      if (no_aux_in_error .and. nerror > m2) nerror = m2
! --- initial preparations
      if (m1.gt.0) ijob=ijob+10
       
      if (dbg) write(*,*) trim(ros_name) // ' ijob', ijob
       
      posneg=sign(1.d0,xend-x)
      hmaxn=min(abs(hmax),abs(xend-x))
      if (abs(h).le.10.d0*uround) h=1.0d-6
      h=min(abs(h),hmaxn) 
      h=sign(h,posneg) 
      reject=.false.
      last=.false.
      nsing=0
      ierr=0
      irtrn=1
      hd=0
! -------- prepare band-widths --------
      mbdiag=mumas+1
      if (banded) then
          mle=mljac
          mue=mujac
          mbjac=mljac+mujac+1
          mbb=mlmas+mumas+1
          mdiag=mle+mue+1
          mdiff=mle+mue-mumas
      end if
      if (iout.ne.0) then 
          xold=x
          irtrn=1
          hout=h
          iwork(1) = lrc
          iwork(2) = n
          j=1+lrc
          rwork(j) = xold; j=j+1
          rwork(j) = h
          call solout(naccpt,xold,x,n,y,rwork,iwork,contro,lrpar,rpar,lipar,ipar,irtrn)
          if (irtrn.lt.0) goto 179
      end if
      
! --- basic integration step  
 1    if (nstep.gt.nmax) goto 178
      hprev = h
      if (0.1d0*abs(h).le.abs(x)*uround) goto 177
      if (last) then
          h=hopt
          idid=1
          goto 191
      end if
      hopt=h
      if ((x+h*1.0001d0-xend)*posneg.ge.0.d0) then
         h=xend-x
         last=.true.
      end if
! *** *** *** *** *** *** ***
!  computation of the jacobian
! *** *** *** *** *** *** ***      
      njac=njac+1
      if (ijac.eq.0) then
         if (nvar > 0) then
            call fcn_blk_dble(n,caller_id,nvar,nz,x,h,y,dy1,lrpar,rpar,lipar,ipar,ierr)
         else
            call fcn(n,x,h,y,dy1,lrpar,rpar,lipar,ipar,ierr)
         endif
         if (ierr /= 0) goto 180
         nfcn=nfcn+1
! --- compute jacobian matrix numerically
         if (banded) then
! --- jacobian is banded
            mujacp=mujac+1
            md=min(mbjac,n)
            do mm=1,m1/m2+1
               do k=1,md
                  j=k+(mm-1)*m2
 12               ak(j,2)=y(j)
                  ak(j,3)=dsqrt(uround*max(1.d-5,abs(y(j))))
                  y(j)=y(j)+ak(j,3)
                  j=j+md
                  if (j.le.mm*m2) goto 12 
                  
                  p1(1:n) => ak(1:n,1)
                  if (nvar > 0) then
                     call fcn_blk_dble(n,caller_id,nvar,nz,x,h,y,p1,lrpar,rpar,lipar,ipar,ierr)
                  else
                     call fcn(n,x,h,y,p1,lrpar,rpar,lipar,ipar,ierr)
                  endif
                  
                  if (ierr /= 0) goto 180
                  j=k+(mm-1)*m2
                  j1=k
                  lbeg=max(1,j1-mujac)+m1
 14               lend=min(m2,j1+mljac)+m1
                  y(j)=ak(j,2)
                  mujacj=mujacp-j1-m1
                  do l=lbeg,lend
                     fjac(l+mujacj,j)=(ak(l,1)-dy1(l))/ak(j,3) 
                  end do
                  j=j+md
                  j1=j1+md
                  lbeg=lend+1
                  if (j.le.mm*m2) goto 14
               end do
            end do
         else
! --- jacobian is full
            do i=1,n
               ysafe=y(i)
               delt=dsqrt(uround*max(1.d-5,abs(ysafe)))
               y(i)=ysafe+delt
               p1(1:n) => ak(1:n,1)
               if (nvar > 0) then
                  call fcn_blk_dble(n,caller_id,nvar,nz,x,h,y,p1,lrpar,rpar,lipar,ipar,ierr)
               else
                  call fcn(n,x,h,y,p1,lrpar,rpar,lipar,ipar,ierr)
               endif
               
               if (ierr /= 0) goto 180
               do j=m1+1,n
                 fjac(j-m1,i)=(ak(j,1)-dy1(j))/delt
               end do
               y(i)=ysafe
            end do
         end if
      else
! --- compute jacobian matrix analytically
         if (nzmax == 0) then
            if (dbg) write(*,11) 'jac y(:)', y(1:min(4,n))
            
            if (nvar > 0) then
               call jac_blk_dble(n,caller_id,nvar,nz,x,h,y,dy1,uf_lblk,uf_dblk,uf_ublk,lrpar,rpar,lipar,ipar,ierr)
            else
               call jac(n,x,h,y,dy1,fjac,ldjac,lrpar,rpar,lipar,ipar,ierr)
            endif
            
            if (dbg) write(*,11) 'jac dy1(:)', dy1(1:min(4,n))
         else
            call sjac(n,x,h,y,dy1,nzmax,ia,ja,sparse_jac,lrpar,rpar,lipar,ipar,ierr)
         end if
         if (ierr /= 0) goto 180
      end if
      if (.not.autnms) then
         if (idfx.eq.0) then
! --- compute numerically the derivative with respect to x
            delt=sqrt(uround*max(1.d-5,abs(x)))
            xdelt=x+delt
            p1(1:n) => ak(1:n,1)
            if (nvar > 0) then
               call fcn_blk_dble(n,caller_id,nvar,nz,xdelt,h,y,p1,lrpar,rpar,lipar,ipar,ierr)
            else
               call fcn(n,xdelt,h,y,p1,lrpar,rpar,lipar,ipar,ierr)
            endif
            
            if (ierr /= 0) goto 180
            do j=1,n
               fx(j)=(ak(j,1)-dy1(j))/delt
            end do
         else
! --- compute analytically the derivative with respect to x
            call dfx(n,x,y,fx,lrpar,rpar,lipar,ipar,ierr)
            if (ierr /= 0) goto 180
         end if
      end if
   2  continue
      if (h <= hprev*1d-30) goto 190
      ierr=0
! *** *** *** *** *** *** ***
!  compute the stages
! *** *** *** *** *** *** ***
      fac=1.d0/(h*gamma)
      if (need_free) then
         call decsol_done(n,fjac,ldjac,fmas,ldmas,mlmas,mumas,
     &            m1,m2,nm1,fac,e1,lde,ip,ak1,ier,ijob,implct,ip,
     &            mle,mue,mbjac,mbb,mdiag,mdiff,mbdiag,
     &            decsol,decsols,decsolblk,
     &            caller_id, nvar, nz, lblk, dblk, ublk, 
     &            sparse_jac,nzmax,isparse,ia,ja,sa,lrd,rpar_decsol,lid,ipar_decsol)
      end if
      call decomr(n,fjac,ldjac,fmas,ldmas,mlmas,mumas,
     &            m1,m2,nm1,fac,e1,lde,ip,ak1,ier,ijob,implct,ip,
     &            mle,mue,mbjac,mbb,mdiag,mdiff,mbdiag,
     &            decsol,decsols,decsolblk,
     &            caller_id, nvar, nz, lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk,
     &            sparse_jac,nzmax,isparse,ia,ja,sa,lrd,rpar_decsol,lid,ipar_decsol)
      need_free = .true.
      if (ier.ne.0) goto 80
      if (dbg .and. .false.) then
         write(*,11) 'fac', fac
         write(*,*) 'n', n
         write(*,*) 'ldjac', ldjac
         write(*,*) 'ldmas', ldmas
         write(*,*) 'mlmas', mlmas
         write(*,*) 'mumas', mumas
         write(*,*) 'm1', m1
         write(*,*) 'm2', m2
         write(*,*) 'nm1', nm1
         write(*,*) 'lde', lde
         do i=1,lde
            write(*,11) 'e', e(i,1), e(i,2), e(i,min(3,n)), e(i,min(4,n))
         end do
         write(*,*)  
      end if
      ndec=ndec+1
! --- prepare for the computation of the stages
      if (.not.autnms) hd = h*ros_gamma
      if (h <= 0d0 .or. is_bad(h)) goto 177
      rhc = rc/h
     
! ------------ stages ---------------
      if (dbg) then
         write(*,*)
         write(*,*) 'nstep', nstep
         write(*,11) 'x', x
         write(*,11) 'h', h
         if (nstep > 30) stop 'rosenbrock nstep > 1'
      end if

      if (is_bad(h)) goto 177
      
      do is=1,ns
         if (dbg) write(*,*) 'stage', is
         if (is == 1) then
            dy = dy1
            not_stage1 = .false.
         else
            if (ros_newf(is)) then
               do j=1,n
                  ynew(j) = y(j) + sum(ra(is,1:is-1)*ak(j,1:is-1))
                  if (ynew(j) < y_min .or. ynew(j) > y_max) then
                     if (dbg) write(*,*) 'stage ynew(j) < y_min .or. ynew(j) > y_max', 
     >                     is, j, ynew(j), y(j), sum(ra(is,1:is-1)*ak(j,1:is-1))
                     goto 82
                  end if
               end do
               if (dbg) write(*,11) 'fcn y(:)', ynew(1:min(4,n))
               
               if (nvar > 0) then
                  call fcn_blk_dble(n,caller_id,nvar,nz,x+ros_alpha(is)*h,h,ynew,dy,lrpar,rpar,lipar,ipar,ierr)
               else
                  call fcn(n,x+ros_alpha(is)*h,h,ynew,dy,lrpar,rpar,lipar,ipar,ierr)
               endif
               
               if (dbg) write(*,11) 'fcn dy(:)', dy(1:min(4,n))
               if (ierr /= 0) goto 81
               if (rd32 /= 0) dy2(1:n) = dy(1:n)
            end if
            do j=1,n
               cont(j) = sum(rhc(is,1:is-1)*ak(j,1:is-1))
            end do
            if (is == 3 .and. rd32 /= 0) then
               do j=1,n
                  cont(j) = cont(j) + rd32*dy2(j) + rd(3,1)*dy1(j)
               end do
            end if
            not_stage1 = .true.
         end if
         p1(1:n) => ak1(1+(is-1)*n:n*is)
         call slvrod(n,fjac,ldjac,mljac,mujac,fmas,ldmas,mlmas,mumas,
     &      m1,m2,nm1,fac,e1,lde,ip,dy,p1,fx,cont,hd(is),ijob,not_stage1,
     &      mle,mue,mbjac,mbb,mdiag,mdiff,mbdiag,
     &      decsol,decsols,decsolblk,
     &      caller_id, nvar, nz, lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk,
     &      nzmax,isparse,ia,ja,sa,lrd,rpar_decsol,lid,ipar_decsol,ierr)
         if (ierr /= 0) goto 81
         if (dbg) write(*,11) 'ak(:,is)', ak(1,is), ak(2,is), ak(min(3,n),is), ak(min(4,n),is)
         if (dbg) write(*,*)
      end do

! ------------ new solution ---------------
      do j=1,n
         ynew(j) = y(j) + sum(ros_m(1:ns)*ak(j,1:ns))
         if (ynew(j) < y_min .or. ynew(j) > y_max) then
            if (dbg) write(*,*) 'new solution ynew(j) < y_min .or. ynew(j) > y_max', j, ynew(j)
            goto 82
         end if
      end do
      nsol=nsol+ns
      nfcn=nfcn+ns-1 
      
! ------------ error estimation ----------
      nstep=nstep+1
      ! put the error vector in dy
      do i=1,n
         dy(i) = sum(ros_e(1:ns)*ak(i,1:ns))
      end do
      
      err=0.d0
      do i=1,nerror
         if (itol.eq.0) then
            sk=atol(1)+rtol(1)*max(abs(y(i)),abs(ynew(i)))
         else
            sk=atol(i)+rtol(i)*max(abs(y(i)),abs(ynew(i)))
         !if (dbg) write(*,*) 'sk', i, sk
         end if
         err=err+(dy(i)/sk)**2
         !if (dbg) write(*,*) 'err', i, err
      end do
      
      err=sqrt(err/nerror)
      if (is_bad(err)) goto 81
      
      
 11    format(a20,2x,99(e26.16,1x))
 
      if (dbg) then
      
         write(*,11) 'y(:)', y(1:min(4,n))
         write(*,11) 'ynew(:)', ynew(1:min(4,n))
         write(*,11) 'err(:)', dy(1:min(4,n))
         write(*,11) 'err=', err
         write(*,*)
         write(*,*)
      
         if (nstep > 20) stop
      
         if (is_bad(y(1)) .or. is_bad(y(2))) stop
      
      end if
      

! --- computation of hnew
! --- we require .2<=hnew/h<=6.
      eloi = 1d0/ros_elo ! inverse of estimated local order
      fac=max(fac2,min(fac1,pow(err,eloi)/safe))


      
      if (minval(ynew(1:n)) < y_min) then
         if (dbg) write(*,*) 'reject < y_min', minval(ynew(1:n)), y_min
         fac = 100
         err = 100
      else if (maxval(ynew(1:n)) > y_max) then
         if (dbg) write(*,*) 'reject > y_max', maxval(ynew(1:n)), y_max
         fac = 100
         err = 100
      end if


      hnew=h/fac 
      if (is_bad(hnew)) goto 81

      hopt=hnew 
      
! *** *** *** *** *** *** ***
!  is the error small enough ?
! *** *** *** *** *** *** ***
      if (err.le.1.d0) then
! --- step is accepted  
         naccpt=naccpt+1
         if (pred) then
c       --- predictive controller of gustafsson
            if (naccpt.gt.1) then
               facgus=(hacc/h)*pow(err**2/erracc,eloi)/safe
               facgus=max(fac2,min(fac1,facgus))
               fac=max(fac,facgus)
               hnew=h/fac
            end if
            hacc=h
            erracc=max(1.0d-2,err)
         end if
         if (iout.ne.0) then
            do i=1,n
               cont(i)=y(i)
            end do
         end if
         do i=1,n 
            y(i)=ynew(i)
         end do
         xold=x 
         x=x+h
         if (iout.ne.0) then 
            do i=1,n
               cont(i+nn)=ynew(i)
            end do
            if (ros_elo >= 4) then
               do i=1,n
                  cont(i+nn2)=sum(ros_d(2,1:ns-1)*ak(i,1:ns-1))
                  cont(i+nn3)=sum(ros_d(3,1:ns-1)*ak(i,1:ns-1))
               end do
               if (ros_elo == 5) then
                  do i=1,n
                     cont(i+nn4)=sum(ros_d(4,1:ns-1)*ak(i,1:ns-1))
                  end do
               end if
            else
               do i=1,n
                  cont(i+nn2) = dy1(i)
               end do
            end if
            irtrn=1
            hout=h
            iwork(1) = lrc
            iwork(2) = n
            j=1+lrc
            rwork(j) = xold; j=j+1
            rwork(j) = h
            call solout(naccpt,xold,x,n,y,rwork,iwork,contro,lrpar,rpar,lipar,ipar,irtrn)
            if (irtrn.lt.0) goto 179
         end if
         if (abs(hnew).gt.hmaxn) hnew=posneg*hmaxn
         if (reject) hnew=posneg*min(abs(hnew),abs(h)) 
         reject=.false.
         h=hnew
         goto 1
      else
! --- step is rejected  
         reject=.true.
         last=.false.
         h=hnew         
         if (naccpt.ge.1) nrejct=nrejct+1
         goto 2
      end if
! --- singular matrix
  80  nsing=nsing+1
      if (nsing.ge.5) goto 176
      h=h*0.5d0
      reject=.true.
      last=.false.
      goto 2
! --- step rejected
  81  h=h*0.5d0
      reject=.true.
      last=.false.
      goto 2
! --- step rejected for y_min or y_max
  82  h=h*0.1d0
      reject=.true.
      last=.false.
      goto 2
! --- fail exit
 176  continue
      if (lout > 0) write(lout,979)x   
      if (lout > 0) write(lout,*) ' matrix is repeatedly singular, ier=',ier
      idid=-4
      goto 191
 177  continue
      if (lout > 0) write(lout,979)x   
      if (lout > 0) write(lout,*) ' step size too small, h=',h
      idid=-3
      goto 191
 178  continue
      if (lout > 0) write(lout,979)x   
      if (lout > 0) write(lout,*) ' more than nmax =',nmax,'steps are needed' 
      idid=-2
      goto 191
! --- solout exit
 179  continue
      if (lout > 0) write(lout,979)x
 979  format(' exit at x=',e18.4) 
      idid=2
      goto 191
! --- forced exit because of ierr /= 0
 180  continue
      idid=-5
      goto 191
! --- too many reductions in stepsize and still not okay error
 190  continue
      idid=-8


 191  continue
 
      if (need_free) then
         p1(1:n) => ak(1:n,1)
         call decsol_done(n,fjac,ldjac,fmas,ldmas,mlmas,mumas,
     &            m1,m2,nm1,fac,e1,lde,ip,p1,ier,ijob,implct,ip,
     &            mle,mue,mbjac,mbb,mdiag,mdiff,mbdiag,
     &            decsol,decsols,decsolblk,
     &            caller_id, nvar, nz, lblk, dblk, ublk, 
     &            sparse_jac,nzmax,isparse,ia,ja,sa,lrd,rpar_decsol,lid,ipar_decsol)
      end if
      
      return
      end subroutine roscor


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
      subroutine ros2_coeffs (ns,ra,rc,rd,ros_d,ros_m,ros_e,ros_alpha,
     &                ros_gamma,ros_newf,ros_elo,no_aux_in_error,ros_name)

!       Rosenbrock Method "ROS2"
!       CWI, MAS-R9717
!       J.G.Verwer, E.J.Spee, J.G.Blom, W.H.Hundsdorfer

! --- an L-stable method, 2 stages, order 2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
      implicit none
      integer, intent(in) :: ns
      real(dp), intent(inout) :: ros_m(ns), ros_e(ns)
      real(dp), intent(inout) :: ros_d(ns,ns), ra(ns,ns), rc(ns,ns), rd(ns,ns)
      real(dp), intent(inout) :: ros_alpha(ns), ros_gamma(ns)
      integer, intent(out) :: ros_elo
      logical, intent(out) :: no_aux_in_error, ros_newf(ns)
      character*12, intent(out) :: ros_name
      
      real(dp) :: g
       no_aux_in_error = .true.
       g = 1.0d0 + 1.0d0/sqrt(2.0d0)
      
!~~~> name of the method
       ros_name = 'ros2'      
!~~~> number of stages
      if (ns /= 2) stop 'bad ns arg for ros2_coeffs'
       
       rd = 0
       ros_d = 0
      
       ra(2,1) = (1.d0)/g
       rc(2,1) = (-2.d0)/g
       
!~~~> does the stage i require a new function evaluation (ros_newf(i)=true)
!   or does it re-use the function evaluation from stage i-1 (ros_newf(i)=false)
       ros_newf(1) = .true.
       ros_newf(2) = .true.
!~~~> m_i = coefficients for new step solution
       ros_m(1)= (3.d0)/(2.d0*g)
       ros_m(2)= (1.d0)/(2.d0*g)
! e_i = coefficients for error estimator       
       ros_e(1) = 1.d0/(2.d0*g)
       ros_e(2) = 1.d0/(2.d0*g)
!~~~> ros_elo = estimator of local order
       ros_elo = 2   
!~~~> y_stage_i ~ y( t + h*alpha_i )
       ros_alpha(1) = 0.0d0
       ros_alpha(2) = 1.0d0 
!~~~> gamma_i = \sum_j  gamma_{i,j}       
       ros_gamma(1) = g
       ros_gamma(2) =-g
      
      return
      end subroutine ros2_coeffs


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
      subroutine rose2_coeffs (ns,ra,rc,rd,ros_d,ros_m,ros_e,ros_alpha,
     &                ros_gamma,ros_newf,ros_elo,no_aux_in_error,ros_name)

!       Rosenbrock Method "ROSE2"
!       Shampine & Reichelt, SIAM J Sci. Comput., 18, (1997) 1-22.
!       ode23s from matlab ode suite.
! --- an l-stable method, 3 stages, order 2
! --- final function evaluation uses the new state.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
      implicit none
      integer, intent(in) :: ns
      real(dp), intent(inout) :: ros_m(ns), ros_e(ns), ros_d(ns,ns), ra(ns,ns), rc(ns,ns), rd(ns,ns)
      real(dp), intent(inout) :: ros_alpha(ns), ros_gamma(ns)
      integer, intent(out) :: ros_elo
      logical, intent(out) :: no_aux_in_error, ros_newf(ns)
      character*12, intent(out) :: ros_name
      real(dp) :: g, e32
      real(dp), parameter :: sqrt2 = 1.4142135623731d0 ! sqrt(2d0)
       no_aux_in_error = .true.
       g = 1.0d0/(2d0 + sqrt2)       
       e32 = 6d0 + sqrt2
!~~~> name of the method
      ros_name = 'rose2'      
!~~~> number of stages
      if (ns /= 3) stop 'bad ns arg for rose2_coeffs'
      ros_d = 0
      
!~~~> the coefficient matrices a and c are strictly lower triangular.
     
      ra(2,1)= 0.5d0/g
      ra(3,1)= 1.0d0/g
      ra(3,2)= 1.0d0/g

      rc(2,1) = -1.0d0/g
      rc(3,1) = -(e32+2d0)/g
      rc(3,2) = -e32/g

      rd(2,1) = 0d0
      rd(3,1) = 2d0
      rd(3,2) = e32
      
!~~~> does the stage i require a new function evaluation (ros_newf(i)=true)
!   or does it re-use the function evaluation from stage i-1 (ros_newf(i)=false)
      ros_newf(1) = .true.
      ros_newf(2) = .true.
      ros_newf(3) = .true.
      
!~~~> m_i = coefficients for new step solution
      ros_m(1) = 1/g
      ros_m(2) = 1/g
      ros_m(3) = 0
      
! e_i = coefficients for error estimator
      ros_e(1) = -1/(6*g)
      ros_e(2) = -1/(3*g)
      ros_e(3) = 1/(6*g)
      
!~~~> ros_elo = estimator of local order
      ros_elo = 2
      
!~~~> y_stage_i ~ y( t + h*alpha_i )
      ros_alpha(1)= 0.0d0
      ros_alpha(2)= 0.5d0
      ros_alpha(3)= 1.0d0
      
!~~~> gamma_i = \sum_j  gamma_{i,j}
      ros_gamma(1)= g
      ros_gamma(2)= 0
      ros_gamma(3)= g
      
      return
      end subroutine rose2_coeffs

      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
      subroutine ros3p_coeffs (ns,ra,rc,rd,ros_d,ros_m,ros_e,ros_alpha,
     &               ros_gamma,ros_newf,ros_elo,no_aux_in_error,ros_name)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
       
!       Rosenbrock Method "ROS3P"
!       Visit at CWI, January 2000
!       J.Lang, J.G.Verwer

! --- A-stable, 3 stages, order 3, 2 function evaluations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
      implicit none
      integer, intent(in) :: ns
      real(dp), intent(inout) :: ros_m(ns), ros_e(ns), ros_d(ns,ns), ra(ns,ns), rc(ns,ns), rd(ns,ns)
      real(dp), intent(inout) :: ros_alpha(ns), ros_gamma(ns)
      integer, intent(out) :: ros_elo
      logical, intent(out) :: no_aux_in_error, ros_newf(ns)
      character*12, intent(out) :: ros_name
      no_aux_in_error = .true.
!~~~> name of the method
      ros_name = 'ros3p'      
!~~~> number of stages
      if (ns /= 3) stop 'bad ns arg for ros3p_coeffs'
      ra = 0
      rc = 0
      rd = 0
      ros_d = 0
      
!~~~> the coefficient matrices a and c are strictly lower triangular.
     
      ra(2,1)= 1.26794919243112273d0
      ra(3,1)= 1.26794919243112273d0
      ra(3,2)= 0.d0

      rc(2,1) = -1.60769515458673630d0
      rc(3,1) = -3.46410161513775469d0
      rc(3,2) = -1.73205080756887734d0
      
!~~~> does the stage i require a new function evaluation (ros_newf(i)=true)
!   or does it re-use the function evaluation from stage i-1 (ros_newf(i)=false)
      ros_newf(1) = .true.
      ros_newf(2) = .true.
      ros_newf(3) = .false.
!~~~> m_i = coefficients for new step solution
      ros_m(1) =  2.0d0
      ros_m(2) =  0.57735026918962578d0
      ros_m(3) = 0.42264973081037424d0
! e_i = coefficients for error estimator
      ros_e(1) = ros_m(1) - 2.11324865405187124d0
      ros_e(2) = ros_m(2) - 1.0d0
      ros_e(3) = ros_m(3) - 0.42264973081037424d0
!~~~> ros_elo = estimator of local order
      ros_elo = 3
!~~~> y_stage_i ~ y( t + h*alpha_i )
      ros_alpha(1)= 0.0d+00
      ros_alpha(2)= 1d+00
      ros_alpha(3)= 1d+00
!~~~> gamma_i = \sum_j  gamma_{i,j}
      ros_gamma(1)= 0.78867513459481286d0
      ros_gamma(2)= -0.21132486540518713d0
      ros_gamma(3)= -1.07735026918962573d0
      return
      end subroutine ros3p_coeffs

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      


      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
      subroutine ros3pl_coeffs (ns,ra,rc,rd,ros_d,ros_m,ros_e,ros_alpha,
     &               ros_gamma,ros_newf,ros_elo,no_aux_in_error,ros_name)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
! --- a stiffly-stable method for parabolic equations; 4 stages, order 3, 3 function evaluations.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
      implicit none
      integer, intent(in) :: ns
      real(dp), intent(inout) :: ros_m(ns), ros_e(ns), ros_d(ns,ns), ra(ns,ns), rc(ns,ns), rd(ns,ns)
      real(dp), intent(inout) :: ros_alpha(ns), ros_gamma(ns)
      integer, intent(out) :: ros_elo
      logical, intent(out) :: no_aux_in_error, ros_newf(ns)
      character*12, intent(out) :: ros_name
      no_aux_in_error = .true.
!~~~> name of the method
      ros_name = 'ros3pl'      
!~~~> number of stages
      if (ns /= 4) stop 'bad ns arg for ros3pl_coeffs'
      ra = 0
      rc = 0
      rd = 0
      ros_d = 0
      
!~~~> the coefficient matrices a and c are strictly lower triangular.
     
      ! matA (alpha21, alpha31, ...)
      ra(2,1) = 1.147140180139521d0
      ra(3,1) = 2.463070773030053d0
      ra(3,2) = ra(2,1)
      ra(4,1) = ra(3,1)
      ra(4,2) = ra(2,1)
      ra(4,3) = 0.d0
      
      ! -matC  (c21, c31, ...)
      rc(2,1) = -2.631861185781065d0
      rc(3,1) = -1.302364158113095d0
      rc(3,2) =  2.769432022251304d0
      rc(4,1) = -1.552568958732400d0
      rc(4,2) =  2.587743501215153d0
      rc(4,3) = -1.416993298352020d0
      
!~~~> does the stage i require a new function evaluation (ros_newf(i)=true)
!   or does it re-use the function evaluation from stage i-1 (ros_newf(i)=false)
      ros_newf(1)  = .true.
      ros_newf(2)  = .true.
      ros_newf(3)  = .true.
      ros_newf(4)  = .false.
      
!~~~> m_i = coefficients for new step solution
      ! vecB1  (m1, m2, ...)
      ros_m(1) = 2.463070773030053d0
      ros_m(2) = 1.147140180139521d0
      ros_m(3) = 0
      ros_m(4) = 1
      
! e_i = coefficients for error estimator
      ! vecB1-vecB2  (m1-mhat1, ...)
      ros_e(1) = ros_m(1) - 2.346947683513665d0
      ros_e(2) = ros_m(2) - 0.456530569451895d0
      ros_e(3) = ros_m(3) - 0.056949243945495d0
      ros_e(4) = ros_m(4) - 0.738684936166224d0
      
!~~~> ros_elo = estimator of local order
      ros_elo = 3
      
!~~~> y_stage_i ~ y( t + h*alpha_i )
      ! vecA (alpha1, alpha2, ...)
      ros_alpha(1)= 0
      ros_alpha(2)= 0.5d0
      ros_alpha(3)= 1
      ros_alpha(4)= 1
      
!~~~> gamma_i = \sum_j  gamma_{i,j}
      ! vecG (gamma1, gamma2,...)
      ros_gamma(1)= 0.435866521508459d0
      ros_gamma(2)= -0.064133478491541d0
      ros_gamma(3)= 0.111028172512505d0
      ros_gamma(4)= 0
      
      return
      end subroutine ros3pl_coeffs

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
      subroutine rodas3_coeffs (ns,ra,rc,rd,ros_d,ros_m,ros_e,ros_alpha,
     &                ros_gamma,ros_newf,ros_elo,no_aux_in_error,ros_name)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
! --- a stiffly-stable method; 4 stages, order 3, 3 function evaluations.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
      implicit none
      integer, intent(in) :: ns
      real(dp), intent(inout) :: ros_m(ns), ros_e(ns), ros_d(ns,ns), ra(ns,ns), rc(ns,ns), rd(ns,ns)
      real(dp), intent(inout) :: ros_alpha(ns), ros_gamma(ns)
      integer, intent(out) :: ros_elo
      logical, intent(out) :: no_aux_in_error, ros_newf(ns)
      character*12, intent(out) :: ros_name
      ra = 0
      rc = 0
      rd = 0
      ros_d = 0
      no_aux_in_error = .false.
!~~~> name of the method
      ros_name = 'rodas3'      
!~~~> number of stages
      if (ns /= 4) stop 'bad ns arg for rodas3_coeffs'
      
!~~~> the coefficient matrices a and c are strictly lower triangular.
      ra(2,1) = 0.0d+00
      ra(3,1) = 2.0d+00
      ra(3,2) = 0.0d+00
      ra(4,1) = 2.0d+00
      ra(4,2) = 0.0d+00
      ra(4,3) = 1.0d+00

      rc(2,1) = 4.0d+00
      rc(3,1) = 1.0d+00
      rc(3,2) =-1.0d+00
      rc(4,1) = 1.0d+00
      rc(4,2) =-1.0d+00 
      rc(4,3) =-(8.0d+00/3.0d+00) 
               
!~~~> does the stage i require a new function evaluation (ros_newf(i)=true)
!   or does it re-use the function evaluation from stage i-1 (ros_newf(i)=false)
      ros_newf(1)  = .true.
      ros_newf(2)  = .false.
      ros_newf(3)  = .true.
      ros_newf(4)  = .true.
!~~~> m_i = coefficients for new step solution
      ros_m(1) = 2.0d+00
      ros_m(2) = 0.0d+00
      ros_m(3) = 1.0d+00
      ros_m(4) = 1.0d+00
!~~~> e_i  = coefficients for error estimator       
      ros_e(1) = 0.0d+00
      ros_e(2) = 0.0d+00
      ros_e(3) = 0.0d+00
      ros_e(4) = 1.0d+00
!~~~> ros_elo  = estimator of local order
      ros_elo  = 3     
!~~~> y_stage_i ~ y( t + h*alpha_i )
      ros_alpha(1) = 0.0d+00
      ros_alpha(2) = 0.0d+00
      ros_alpha(3) = 1.0d+00
      ros_alpha(4) = 1.0d+00
!~~~> gamma_i = \sum_j  gamma_{i,j}       
      ros_gamma(1) = 0.5d+00
      ros_gamma(2) = 1.5d+00
      ros_gamma(3) = 0.0d+00
      ros_gamma(4) = 0.0d+00
      return
      end subroutine rodas3_coeffs
      
      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
      subroutine rodas4_coeffs (ns,ra,rc,rd,ros_d,ros_m,ros_e,ros_alpha,
     &                 ros_gamma,ros_newf,ros_elo,no_aux_in_error,ros_name)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!     stiffly-stable rosenbrock method of order 4, with 6 stages
!
!         e. hairer and g. wanner, solving ordinary differential
!         equations ii. stiff and differential-algebraic problems.
!         springer series in computational mathematics,
!         springer-verlag (1996)     
!          
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

      implicit none
      integer, intent(in) :: ns
      real(dp), intent(inout) :: 
     >      ros_m(ns), ros_e(ns), ros_d(ns,ns), ra(ns,ns), rc(ns,ns), rd(ns,ns)
      real(dp), intent(inout) :: ros_alpha(ns), ros_gamma(ns)
      integer, intent(out) :: ros_elo
      logical, intent(out) :: no_aux_in_error, ros_newf(ns)
      character*12, intent(out) :: ros_name
      
      no_aux_in_error = .false.
      rd = 0
!~~~> name of the method
       ros_name = 'rodas4'      
!~~~> number of stages
       if (ns /= 6) stop 'bad ns arg for rodas4_coeffs'

!~~~> y_stage_i ~ y( t + h*alpha_i )
       ros_alpha(1) = 0.000d0
       ros_alpha(2) = 0.386d0
       ros_alpha(3) = 0.210d0 
       ros_alpha(4) = 0.630d0
       ros_alpha(5) = 1.000d0
       ros_alpha(6) = 1.000d0
   
!~~~> gamma_i = \sum_j  gamma_{i,j}       
       ros_gamma(1) = 0.2500000000000000d+00
       ros_gamma(2) =-0.1043000000000000d+00
       ros_gamma(3) = 0.1035000000000000d+00
       ros_gamma(4) =-0.3620000000000023d-01
       ros_gamma(5) = 0.0d0
       ros_gamma(6) = 0.0d0

!~~~> the coefficient matrices a and c are strictly lower triangular.
     
       ra(2,1) = 0.1544000000000000d+01
       ra(3,1) = 0.9466785280815826d+00
       ra(3,2) = 0.2557011698983284d+00
       ra(4,1) = 0.3314825187068521d+01
       ra(4,2) = 0.2896124015972201d+01
       ra(4,3) = 0.9986419139977817d+00
       ra(5,1) = 0.1221224509226641d+01
       ra(5,2) = 0.6019134481288629d+01
       ra(5,3) = 0.1253708332932087d+02
       ra(5,4) =-0.6878860361058950d+00
       
       ra(6,1:4) = ra(5,1:4)
       ra(6,5) = 1

       rc(2,1) =-0.5668800000000000d+01
       rc(3,1) =-0.2430093356833875d+01
       rc(3,2) =-0.2063599157091915d+00
       rc(4,1) =-0.1073529058151375d+00
       rc(4,2) =-0.9594562251023355d+01
       rc(4,3) =-0.2047028614809616d+02
       rc(5,1) = 0.7496443313967647d+01
       rc(5,2) =-0.1024680431464352d+02
       rc(5,3) =-0.3399990352819905d+02
       rc(5,4) = 0.1170890893206160d+02
       rc(6,1) = 0.8083246795921522d+01
       rc(6,2) =-0.7981132988064893d+01
       rc(6,3) =-0.3152159432874371d+02
       rc(6,4) = 0.1631930543123136d+02
       rc(6,5) =-0.6058818238834054d+01

!~~~> m_i = coefficients for new step solution
       ros_m(1:4) = ra(5,1:4)
       ros_m(5) = 1
       ros_m(6) = 1

!~~~> e_i  = coefficients for error estimator       
       ros_e(1) = 0.0d+00
       ros_e(2) = 0.0d+00
       ros_e(3) = 0.0d+00
       ros_e(4) = 0.0d+00
       ros_e(5) = 0.0d+00
       ros_e(6) = 1.0d+00

!~~~> does the stage i require a new function evaluation (ros_newf(i)=true)
!   or does it re-use the function evaluation from stage i-1 (ros_newf(i)=false)
       ros_newf(1) = .true.
       ros_newf(2) = .true.
       ros_newf(3) = .true.
       ros_newf(4) = .true.
       ros_newf(5) = .true.
       ros_newf(6) = .true.
       
      ! dense output
      ros_d(2,1)= 0.1012623508344586d+02
      ros_d(2,2)=-0.7487995877610167d+01
      ros_d(2,3)=-0.3480091861555747d+02
      ros_d(2,4)=-0.7992771707568823d+01
      ros_d(2,5)= 0.1025137723295662d+01
      ros_d(3,1)=-0.6762803392801253d+00
      ros_d(3,2)= 0.6087714651680015d+01
      ros_d(3,3)= 0.1643084320892478d+02
      ros_d(3,4)= 0.2476722511418386d+02
      ros_d(3,5)=-0.6594389125716872d+01
     
!~~~> ros_elo  = estimator of local order
       ros_elo = 4
     
      return
      end subroutine rodas4_coeffs
      
      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
      subroutine rodasp_coeffs (ns,ra,rc,rd,ros_d,ros_m,ros_e,ros_alpha,
     &                 ros_gamma,ros_newf,ros_elo,no_aux_in_error,ros_name)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!     stiffly-stable rosenbrock method of order 4, with 6 stages
!     for parabolic equations (G.Steinbach,1993).
!
!         e. hairer and g. wanner, solving ordinary differential
!         equations ii. stiff and differential-algebraic problems.
!         springer series in computational mathematics,
!         springer-verlag (1996)     
!          
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

      implicit none
      integer, intent(in) :: ns
      real(dp), intent(inout) :: ros_m(ns), ros_e(ns), ros_d(ns,ns), ra(ns,ns), rc(ns,ns), rd(ns,ns)
      real(dp), intent(inout) :: ros_alpha(ns), ros_gamma(ns)
      integer, intent(out) :: ros_elo
      logical, intent(out) :: no_aux_in_error, ros_newf(ns)
      character*12, intent(out) :: ros_name
      no_aux_in_error = .false.
      rd = 0
!~~~> name of the method
       ros_name = 'rodasp'      
!~~~> number of stages
       if (ns /= 6) stop 'bad ns arg for rodasp_coeffs'

!~~~> y_stage_i ~ y( t + h*alpha_i )
       ros_alpha(1) = 0.000d0
       ros_alpha(2) = 0.750d0
       ros_alpha(3) = 0.210d0 
       ros_alpha(4) = 0.630d0
       ros_alpha(5) = 1.000d0
       ros_alpha(6) = 1.000d0
   
!~~~> gamma_i = \sum_j  gamma_{i,j}       
       ros_gamma(1) = 0.2500000000000000d+00
       ros_gamma(2) =-0.5000000000000000d+00
       ros_gamma(3) =-0.2350400000000000d-01
       ros_gamma(4) =-0.3620000000000000d-01
       ros_gamma(5) = 0.0d0
       ros_gamma(6) = 0.0d0

!~~~> the coefficient matrices a and c are strictly lower triangular.
       ra(2,1) = 0.3000000000000000d+01
       ra(3,1) = 0.1831036793486759d+01
       ra(3,2) = 0.4955183967433795d+00
       ra(4,1) = 0.2304376582692669d+01
       ra(4,2) =-0.5249275245743001d-01
       ra(4,3) =-0.1176798761832782d+01
       ra(5,1) =-0.7170454962423024d+01
       ra(5,2) =-0.4741636671481785d+01
       ra(5,3) =-0.1631002631330971d+02
       ra(5,4) =-0.1062004044111401d+01
       
       ra(6,1:4) = ra(5,1:4)
       ra(6,5) = 1

       rc(2,1)=-0.1200000000000000d+02
       rc(3,1)=-0.8791795173947035d+01
       rc(3,2)=-0.2207865586973518d+01
       rc(4,1)= 0.1081793056857153d+02
       rc(4,2)= 0.6780270611428266d+01
       rc(4,3)= 0.1953485944642410d+02
       rc(5,1)= 0.3419095006749676d+02
       rc(5,2)= 0.1549671153725963d+02
       rc(5,3)= 0.5474760875964130d+02
       rc(5,4)= 0.1416005392148534d+02
       rc(6,1)= 0.3462605830930532d+02
       rc(6,2)= 0.1530084976114473d+02
       rc(6,3)= 0.5699955578662667d+02
       rc(6,4)= 0.1840807009793095d+02
       rc(6,5)=-0.5714285714285717d+01

!~~~> m_i = coefficients for new step solution
       ros_m(1:4) = ra(5,1:4)
       ros_m(5) = 1
       ros_m(6) = 1

!~~~> e_i  = coefficients for error estimator       
       ros_e(1) = 0.0d+00
       ros_e(2) = 0.0d+00
       ros_e(3) = 0.0d+00
       ros_e(4) = 0.0d+00
       ros_e(5) = 0.0d+00
       ros_e(6) = 1.0d+00

!~~~> does the stage i require a new function evaluation (ros_newf(i)=true)
!   or does it re-use the function evaluation from stage i-1 (ros_newf(i)=false)
       ros_newf(1) = .true.
       ros_newf(2) = .true.
       ros_newf(3) = .true.
       ros_newf(4) = .true.
       ros_newf(5) = .true.
       ros_newf(6) = .true.
       
      ! dense output
      ros_d(2,1)= 0.2509876703708589d+02
      ros_d(2,2)= 0.1162013104361867d+02
      ros_d(2,3)= 0.2849148307714626d+02
      ros_d(2,4)=-0.5664021568594133d+01
      ros_d(2,5)= 0.0000000000000000d+00
      ros_d(3,1)= 0.1638054557396973d+01
      ros_d(3,2)=-0.7373619806678748d+00
      ros_d(3,3)= 0.8477918219238990d+01
      ros_d(3,4)= 0.1599253148779520d+02
      ros_d(3,5)=-0.1882352941176471d+01
     
!~~~> ros_elo  = estimator of local order
       ros_elo = 4
     
      return
      end subroutine rodasp_coeffs


      
      ! continuous output routines


      real(dp) function contro3(i,x,rwork,iwork,ierr)
         integer, intent(in) :: i
         real(dp), intent(in) :: x
         real(dp), intent(inout), target :: rwork(*)
         integer, intent(inout), target :: iwork(*)
         integer, intent(out) :: ierr
      
         real(dp) :: xold, h, dx, s, s2
         real(dp), pointer :: cont(:)
         integer :: lrc, n, j
         ierr=0
         lrc = iwork(1)
         n = iwork(2)
         cont => rwork(1:lrc); j=lrc+1
         xold = rwork(j); j=j+1
         h = rwork(j)
         dx=x-xold
         s = dx/h
         s2 = s**2
         contro3 = (1-s2)*cont(i) + dx*(1-s)*cont(i+2*n) + s2*cont(i+n)
         return
      end function contro3


      real(dp) function contro4(i,x,rwork,iwork,ierr)
         integer, intent(in) :: i
         real(dp), intent(in) :: x
         real(dp), intent(inout), target :: rwork(*)
         integer, intent(inout), target :: iwork(*)
         integer, intent(out) :: ierr
      
         real(dp) :: xold, h, s
         real(dp), pointer :: cont(:)
         integer :: lrc, n, j
         ierr=0
         lrc = iwork(1)
         n = iwork(2)

         cont => rwork(1:lrc); j=1+lrc
         xold = rwork(j); j=j+1
         h = rwork(j)
         s=(x-xold)/h 
         contro4=cont(i)*(1-s)+s*(cont(i+n)+(1-s)*(cont(i+n*2)+s*cont(i+n*3)))
         return
      end function contro4
      
      
      end module mod_rosenbrock
      
      
      

! copyright (c) 2004, ernst hairer

! redistribution and use in source and binary forms, with or without 
! modification, are permitted provided that the following conditions are 
! met:

! - redistributions of source code must retain the above copyright 
! notice, this list of conditions and the following disclaimer.

! - redistributions in binary form must reproduce the above copyright 
! notice, this list of conditions and the following disclaimer in the 
! documentation and/or other materials provided with the distribution.

! this software is provided by the copyright holders and contributors “as 
! is” and any express or implied warranties, including, but not limited 
! to, the implied warranties of merchantability and fitness for a 
! particular purpose are disclaimed. in no event shall the regents or 
! contributors be liable for any direct, indirect, incidental, special, 
! exemplary, or consequential damages (including, but not limited to, 
! procurement of substitute goods or services; loss of use, data, or 
! profits; or business interruption) however caused and on any theory of 
! liability, whether in contract, strict liability, or tort (including 
! negligence or otherwise) arising in any way out of the use of this 
! software, even if advised of the possibility of such damage.

