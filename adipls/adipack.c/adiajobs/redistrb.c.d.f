      program main  
c   
c  redistributes mesh in physical model 
c   
c  modified 9/1/1985 to allow use of central expansion to reset 
c  points in new model interior to first non-zero point in old  
c  model. old expansion may still be used by setting ioldex = 1.
c   
c  modified 17/5/1985. now interpolates in x**(-2)*aa(i,n)  
c  for i = 2 and 4 (these tend to a constant at x = 0), 
c  and sets aa(2,n) interior to first mesh point in original
c  model from expansion of p and other aa-s.
c   
c  modified 21/5/1985 in treatment of convection zone boundary  
c  (see s/r stcvzb) 
c   
c  modified 21/5/1985 to include contribution in stretching 
c  from superadiabatic gradient, to ensure adequate resolution very 
c  near surface.
c   
c  modified 29/3 - 1/4 1989 to include various fixes for stretching
c  polytropic model
c
c  modified 30/5/92, to include term in change in buoyancy frequency,
c  and smoothing with running mean.
c
c  modified 25/7/92, to smooth A4 near centre, particularly 
c  near hydrogen exhaustion
c
c  modified 17/2/98, to allow reading general-format amdl files
c  (including aa(6,.), say)
c
c  Modified 5/7/01, fixing criterion for discontinuity in rho, now
c  based on aa(4,.). This is, for the time being, flagged by a
c  negative ndisc. The old criterion is obtained by selecting
c  a posive ndisc.
c
c  Modified 28/1/02, adding contribution from cx at x = 0.
c
c  Modified 9/6/03, adding term to A(4,.)/x**2 in g-mode term to correct 
c  for under-resolution in very compact convective cores.
c
c  Modified 10/6/03, changing treatment of region near density
c  discontinuity.
c  Also test and possibly correct behaviour of A(1,.) at central
c  point. (Note that expansion in evolution code requires checking
c  with 4He burning!)
c
c  Modified 16/2/05, fixing the new definition of the boundaries
c  of a discontinuity region (applied when ndisc .lt. 0)
c
c  Modified 9/3/05, introducing option to ensure a minimum fraction of
c  points in p-mode region, flagged by cg lt 0.
c
c  Major modification 6/7/05, separating bulk of code into a
c  subroutine for later integration with evolution package.
c
c  Modified 7/3/08, using iresa4 = -1 to reset data(6)
c
c  Modified 22/9/08, introducing test of monotonicity of xsi after vinta
c  integration in s/r srdist, with s/r test_mono. If xsi is not monotonic, 
c  xsi is reset on relevant interval with trapezoidal integration.
c
c  Double precision version
c  ++++++++++++++++++++++++
c
c  Dated 10/3/90
c
      implicit double precision (a-h, o-z)
      include 'adipr.incl'
      parameter (iaa=11)
      dimension x(nnmax),aa(iaa,nnmax),xn(nnmax),
     *  an(iaa,nnmax), data(8)
      common/comgrp/ isprtp, irotcp, omgrtp(nnmax)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  Input parameters (read from s/r srdist):
c
c  nn: number of points in new mesh
c  Default: nn=601
c  Default: icnmsh=0  
c  icase: for icase = 1, 2, 3, 11, 12, 21, 22, 23, 24 set standard parameters.
c     icase = 1, 2, 3 give old parameters, as used up until June 1988.
c     icase = 11 gives new, optimized p mode parameters
c     icase = 1  corresponds to p modes, old parameters
c     icase = 2  corresponds to g modes, old parameters
c     icase = 3  corresponds to intermediate modes, old parameters.
c     icase = 11 corresponds to p modes, new parameters.
c     icase = 12 corresponds to g modes, new parameters.
c  Revision 3/5/10, including cg .lt. 0 (see below).
c  For icase .ge. 20, set new standard nsmth, ndisc, dlxdsc
c     icase = 21  corresponds to p modes, revised parameters
c     icase = 22  corresponds to g modes, revised parameters
c     icase = 23  corresponds to intermediate modes, revised parameters.
c     icase = 24  corresponds to modes in red giants
c
c  Note that cg .lt. 0 is used to flag for ensuring that a reasonable
c  fraction of points correspond to p modes (fraction hardcoded in
c  call of s/r testcg).
c  Default: icase=0   
c  icvzbn: for icvzbn .gt. 0 ensure that there is a mesh point at   
c     the lower boundary of the convective envelope.
c     for icvzbn = 1 interpolate normally in aa(i,.) for i .ne. 4,   
c     and reset aa(4,.) by linear interpolation near boundary.  
c     for icvzbn = 2 interpolate separately below and above  
c     boundary. 
c  Default: icvzbn=0  
c  nsmth: for nsmth .gt. 1, smooth integrand with Gaussian-weighted
c     running mean over nsmth points
c  Default: nsmth=0
c  ndisc: for ndisc ne 0 add abs(ndisc) points where there is a
c     discontinuity in rho
c     NOTE: ndisc lt 0 selects for new (post-05/07/01) criterion
c  Default: ndisc=0
c  dlxdsc: range in ln x over which discontinuity is distributed
c  Default: dlxdsc=1.e-3
c  dlgrmx: minimum (d log rho/d log p) for discontinuity
c  Default: dlgrmx=10
c  cacvzb: amplitude of point added at base of convection zone
c     should be used only in conjunction with smoothing
c  Default: cacvzb=0
c
c  cg: abs(cg) weight for g-mode part of mesh
c     if cg .lt. 0 ensure a minimum numbers of meshpoint in p-mode region
c     The fraction is hard-coded to 0.1 in ppfrac below.
c  Default: cg=1  
c  cx: constant part of mesh function
c  Default: cx=40.
c  ca: weight to near-surface term
c  Default: ca=0.01   
c  cdgr, cddgr: additional weight to term depending on A4 (of somewhat
c     uncertain significance for now)
c  Default: cdgr=0
c  Default: cddgr=0   
c  cdg: weight for gradient in A in interior of model
c     if cdg .gt.0, limit contribution (before scaling by cdg)
c     to partial sum of other terms. If cdg .lt. 0, do not limit
c     (and use abs(cdg) as weight).
c  Default: cdg=0
c  alphsf: cuts off the singularity at the surface of a polytropic
c     model. In setting mesh, A2 is replaced by A2/(1+alphsf*A2).
c  Default: alphsf=0
c
c  adda4: scaling of term to add to A(4,.)/x**2 in g-mode term
c     Actual term added is adda4*A(4,nr)/x(nr)**2, where nr-2 is
c     the innermost convectively stable point.
c  Default: adda4=0.d0
c  accrmn: Minimum value of A(4,.) used in cg term (set .gt. 0 to
c     increase number of points in convective core
c  Default: accrmn=0
c  nout: number of meshpoints in printed output
c  Default: nout=100  
c  cn: control of innermost non-zero meshpoint
c  Default: cn=2
c  irsu, unew: if irsu .gt. 0, reset U at innermost non-zero meshpoint to
c  unew
c  Default: irsu=0
c  Default: unew=3
c  iresa4: if iresa4 = 1, reset A4 before redistribution,
c  to correct possible errors in A4 resulting from near exhaustion
c  of hydrogen
c  if iresa4 = -1, reset data(6) on input, to correct for apparent
c  errors in setting amdl file (added 7/3/08)
c  Default: iresa4=0
c  nmodel: for nmodel .gt. 1 read and redistribute nmodel models.   
c          otherwise only 1.
c  Default: nmodel=0  
c  kmodel: for kmodel .gt. 0, start at model no. kmodel 
c  Default: kmodel=0  
c  itsaml: if itsaml = 1 test input model for conversion errors 
c  Default: itsaml=0  
c  ioldex: if ioldex = 1, old expansion (not using second derivatives   
c     in data) is used. 
c  Default: ioldex=0  
c idsin, idsout: input and output units for file I/O with i_inout = 1
c  Default: idsin=2
c  Default: idsout=3
c
c  setup for simple run, assuming no rotation (isprtp = 0)
c
      i_paramset=0
      i_inout=1
      isprtp=0
c
      call srdist(i_paramset, ierr_param, i_inout, 
     *  x, aa, data, xn, an, nh, nn, ivar, iaa, iaa)
      stop
      end
