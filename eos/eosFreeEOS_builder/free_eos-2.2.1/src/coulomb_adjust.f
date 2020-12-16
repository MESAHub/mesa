C*******************************************************************************
C       Start of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C       Copyright (C) 1996, 1998, 2000, 2001, 2004, 2005, 2006 Alan W. Irwin
C
C       $Id: coulomb_adjust.f 352 2006-04-13 02:12:47Z airwin $
C
C       For the latest version of this source code, please contact
C       Alan W. Irwin
C       Department of Physics and Astronomy
C       University of Victoria,
C       Box 3055
C       Victoria, B.C., Canada
C       V8W 3P6
C       e-mail: irwin@beluga.phys.uvic.ca.
C
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; either version 2 of the License, or
C    (at your option) any later version.
C
C    This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program; if not, write to the Free Software
C    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
C
C       End of the header for a fortran source file for a subroutine
C       of the Free_EOS stellar interior equation of state code
C*******************************************************************************
      subroutine coulomb_adjust(xdh, xmocp,
     &  acon, bcon, scon, occon, odcon, ccon, dcon)
C      subroutine to smoothly join Y(X) = ln(DH solution) for all
C      X = ln Gamma < xdh with a Y(X) = ln(modified OCP solution)
C      for all X = ln Gamma > xmocp using a cubic
C      polynomial that satisfies the total of 6 continuity constraints
C      for Y, Y', and Y'' at xdh and xmocp.  The method used is to determine
C      the cubic using the 4 continuity conditions for Y and Y'', and use
C      the two continuity conditions for Y' to determine dccon = ccon - occon
C      and ddcon = dcon - odcon.  The solution of the two non-linear
C      equations in two unknowns is obtained by NR iteration.
C      input quantities:
C      xdh, xmocp define the boundaries of the dh and modified ocp region
C        as explained above.
C      acon, bcon, scon, occon, and odcon define the original ocp solution:
C       ocp = acon*gamma + (bcon*gamma**scon + occon*X + odcon)
C      output quantities:
C      ccon, dcon are the changed coefficients for the modified OCP result:
C       mocp = acon*gamma + (bcon*gamma**scon + ccon*X + dcon)
      implicit none
      double precision
     &  gamma,
     &  acon, bcon, scon, ccon, dcon, !econ, fcon, gcon,
     &  occon, dccon, odcon, ddcon,
     &  dh, ddh, ddh2,
     &  ocp, docp, docp2,
     &  xdh, xmocp, xvalue, dx
C      Variables used for delta c and delta d solution of two simultaneous
C      equations.
      integer iter
      double precision maxrel, f(2), df(2), df2(2), fsolve(2), 
     &  jacobian(2,2),
     &  docpdc, ddocpdc, ddocp2dc, docpdd, ddocpdd, ddocp2dd,
     &  ddocpdocp, ddocp2docp
C       lapack variables
      double precision lu_lapack(2,2),
     &  row_lapack(2), col_lapack(2), sol_lapack(2),
     &  rcond_lapack, ferr_lapack, berr_lapack,
     &  work_lapack(4*2)
      integer ipiv_lapack(2), iwork_lapack(2), info_lapack
      character*1 equed_lapack

C      find dccon and ddcon (modifications to OCP) such that there is a
C      smooth cubic polynomial transition in gap region.
C
C      Ensure smooth result by satisfying two first derivative requirements
C      for gap region given by eq. 3.3.5 of Press et al.  Must solve
C      these two equations for two unknowns, dccon and ddcon.

C      initial estimates:
      dccon = 0.d0
      ddcon = 0.d0
      ccon = occon + dccon
      dcon = odcon + ddcon
      maxrel = 1.d0
      iter = 1
      xvalue = xdh
      dh = 1.5d0*xvalue - 0.5d0*log(3.d0)
      ddh = 1.5d0
      ddh2 = 0.d0
      xvalue = xmocp
      do while(maxrel.gt.1.d-10.and.iter.lt.10)
        gamma = exp(xvalue)
        ocp = acon*gamma +
     &    (bcon*gamma**scon + ccon*xvalue + dcon)
        docp = acon + 
     &    (scon*bcon*gamma**(-1.d0+scon) + ccon/gamma)
        docp2 = 
     &    ((-1.d0+scon)*scon*bcon*gamma**(-2.d0+scon) -
     &    ccon/(gamma*gamma))
        docpdc = xvalue
        docpdd = 1.d0
        ddocpdc = 1.d0/gamma
        ddocp2dc = -1.d0/(gamma*gamma)
C         transform to ln(ocp) as a function of x
        ddocp2docp = -(docp2*gamma*gamma +
     &    docp*gamma)/
     &    (ocp*ocp)
     &    + 2.d0*docp*gamma*
     &    docp*gamma/(ocp*ocp*ocp)
        docp2 = (docp2*gamma*gamma +
     &    docp*gamma*(1.d0 - docp*gamma/ocp))
     &    /ocp
        ddocp2dd = ddocp2docp*docpdd
        ddocp2dc = (ddocp2dc*gamma*gamma +
     &    ddocpdc*gamma*(1.d0 - 2.d0*docp*gamma/ocp))
     &    /ocp + ddocp2docp*docpdc
        docp = docp*gamma/ocp
        ddocpdocp = -docp/ocp
        ddocpdd = ddocpdocp*docpdd
        ddocpdc = ddocpdc*gamma/ocp + ddocpdocp*docpdc
        docpdc = docpdc/ocp
        docpdd = docpdd/ocp
        ocp = log(ocp)
        f(1) = dh
        df(1) = ddh
        df2(1) = ddh2
        f(2) = ocp
        df(2) = docp
        df2(2) = docp2
        dx = xmocp - xdh
C        From Press et al eq. 3.3.5.
C        A = 1, B = 0.
        fsolve(1) = (f(2)-f(1))/dx
     &    - dx*df2(1)/3.d0
     &    - dx*df2(2)/6.d0
     &    - df(1)
C        A = 0, B = 1.
        fsolve(2) = (f(2)-f(1))/dx
     &    + dx*df2(1)/6.d0
     &    + dx*df2(2)/3.d0
     &    - df(2)
C        Actually store negative of Jacobian in jacobian.  Note only
C        f(2)= ocp, df(2)= docp, and df2(2)=docp2 are functions of
C        ccon and dcon.
        jacobian(1,1) = -(
     &    docpdc/dx
     &    - dx*ddocp2dc/6.d0
     &    )
        jacobian(1,2) = -(
     &    docpdd/dx
     &    - dx*ddocp2dd/6.d0
     &    )
        jacobian(2,1) = -(
     &    docpdc/dx
     &    + dx*ddocp2dc/3.d0
     &    - ddocpdc
     &    )
        jacobian(2,2) = -(
     &    docpdd/dx
     &    + dx*ddocp2dd/3.d0
     &    - ddocpdd
     &    )
!        test section for comparing numerical and analytical derivatives
!        write(90,*) iter+1
!        write(90,*) f(2), docpdc, docpdd
!        write(90,*) df(2), ddocpdc, ddocpdd
!        write(90,*) df2(2), ddocp2dc, ddocp2dd
!        jacobian(1,1) = 10000.d0
!        jacobian(2,2) = 10000.d0
!        jacobian(2,1) = 0.d0
!        jacobian(1,2) = 0.d0
!        fsolve(1) = 0.d0
!        fsolve(2) = 1.d0
C        compute change in dccon and ddcon using lapack routine.
        call dgesvx('E', 'N',
     &    2, 1, jacobian, 2,
     &    lu_lapack, 2, ipiv_lapack, equed_lapack,
     &    row_lapack, col_lapack, fsolve, 2,
     &    sol_lapack, 2,
     &    rcond_lapack, ferr_lapack, berr_lapack,
     &    work_lapack, iwork_lapack, info_lapack)
        if (info_lapack.ne.0) stop 'something wrong with dgesvx setup'
        dccon = dccon + sol_lapack(1)
        ddcon = ddcon + sol_lapack(2)
        iter = iter + 1
        maxrel = max(
     &    abs((occon + dccon - ccon)/ccon),
     &    abs((odcon + ddcon - dcon)/dcon))
        ccon = occon + dccon
        dcon = odcon + ddcon
!        write(90,*) iter, maxrel, dccon, ddcon
      enddo
C      at this stage have found definitive dccon, ddcon, ccon, and dcon
C      using solution of two non-linear equations in two unknowns.
      end
