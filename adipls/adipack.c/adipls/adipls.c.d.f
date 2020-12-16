      subroutine adipls(i_paramset, ierr_param, i_inout, 
     *  x_arg, aa_arg, data_arg, nn_arg, ivarmd, iaa_arg)
c
c  main subroutine for adiabatic oscillation integration
c  *****************************************************
c
c  Version with parameter passing. Dated 7/7/05
c
c  Arguments with parameter passing:
c
c  i_paramset: Controls parameter setting [in a way to be described]
c  ierr_param: Is returned as 0 for successful completion, as .lt. 0
c     in case of error.
c  i_inout: If i_inout = 1 read model from file, as usual.
c     If i_inout = 0, model quantities must be provided in
c     x_arg(1:nn_arg), aa(arg(1:ivarmd,1:nn_arg), data_arg(1:8)
c     If i_inout = -1, assume that model is already stored,
c     regardless of possible model parameters read in.
c
      implicit double precision (a-h,o-z)
      include 'adipls.c.d.incl'
      parameter(iaa = 10, iaa1 = 10, iy = 8)
      logical fulmod,frsrun,radial,nonrad
      character tail*5, cntrd*50, ctrdsf*1, in_stat*2, out_stat*2
      character*280 file, filess, obsrot_file, strcompr
c
      dimension x_arg(*), aa_arg(iaa_arg,*), data_arg(*)
      dimension csummm(50),icsumm(8),ctrdsf(3)
c
c  commons for parameter transmission
c
      common/cadr_param/
     *  para_el, para_els1, para_dels, para_dfsig1, para_dfsig2,
     *  para_sig1, para_sig2, para_dfsig, para_eltrw1, para_eltrw2,
     *  para_sgtrw1, para_sgtrw2 
      common/cadi_param/
     *  ipara_nsel, ipara_nsig1, ipara_nsig2, ipara_itrsig, ipara_nsig, 
     *  ipara_istsig, ipara_inomd1, ipara_iscan, ipara_irotkr
c
      common/csumma/ xmod1,datsum(8),datmd1(2),xtrnc1,dum13,xfit1,
     *  fsbcsm,fcbcsm,albsum,
     *  elsum,ordsum,sigsum,sigc,dmx,xmx,ekin,per,perv,frqv,
     *  ddsig,ddsol,ysum(4),dosdum(5),
     *  in,nnw,mdints,ivarfs,icase,iorign,idum7,idum8,mlname(4)
      common/rhsdat/ el,ell,alb,els,el1,sig,anres,perfac,data(8),
     *  aa(iaa,1)
      common/rotdat/ em, irotsl
      common/bcsdat/ fctsbc, fcttbc, istsbc, ibotbc
      common/nrmchk/ nnww,ncfac,irsdif
      common/cincnt/ xmnevn,xfit, fcnorm, eps, epssol, dsigmx, fsig, 
     *  dsigre,icow,iturpr,iplneq,iriche,nnwwin, moddet, itmax, irsevn, 
     *  nftmax,irsord,iekinr,itsord,imissl,imjssl,imstsl,inomd1
      common/cincnt_new/ fsig0
      common/cnrkbc/ cbc(4,4),cavc(4,4),cavs(4,4),isngbt,isngsf,ibcnr1,
     *  nnq1,nnq2,idgnrk
      common/sysord/ ii
      common/ccgrav/ cgrav
      common/xarra/ x(1)
      common/wrklir/ wwwww(300)
      common/worksp/ aa1(iaa1,1)
      common/yyyyyy/ y(iy,1)
      common/varcon/ ivarf,npvarf,kvarfc
      common/cmdtst/ iordtr,icaswn, sigwn1, sigwn2, frqwn1, frqwn2, 
     *  iorwn1, iorwn2, frlwn1, frlwn2
c
      common/cprcns/ epsprc,epsufl,epsofl,eprufl
      common/cofile/ nfiles, idsfil(99), file(99), iopen(99), filess(99)
      common/cwindw/ eltrw1, eltrw2, sgtrw1, sgtrw2
c
c  common for storage of modal parameters (degree, order, cyclic frequency,
c  inertia)
c
      common/cobs_param/ icobs_st, nobs_st, obs_st(10,1)
c
c  common controlling output
c
      common/coutpt/ nout,nprcen,iper,irotkr,nprtkr,igm1kr,npgmkr,
     *  nfmode,nfmesh,ispcpr,npout,nobs_stmx,nfmscn
c
c  common containing diagnostic flags possibly set in various parts of the 
c  code
c
      common/cdgsbr/ kdgbcs
c
c  common controlling diagnostics
c
      common/cdiagn/ idgrhs, idgrh1, iprdet, itssol, idgtss
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      common/cstdio_def/ istdin_def, istdou_def, istdpr_def, istder_def
c
c  common defining file input and output
c
      common/cdadsg/ imlds,idslog,idsgsm,idsssm,idsefn,
     *   idsrkr,idsgkr,itrds
c
      external adirhs,rhsnrm,bcnrk,blstio
c
      equivalence (xmod1,csummm(1))
      equivalence (csummm(39),icsumm(1))
c
      save
c
c
c  set initial flags for unit numbers
c
      data imldsp, itrdsp, idsgsp, idsssp,  idsefp, idsrkp, 
     *     idsgkp, istdpp
     *     /  -1,     -1,     -1,     -1,     -1,     -1,
     *        -1,     -1  /
c
      data tail /'    @'/
      data ctrdsf/ 'u', 'u', 'f'/
c
c  initialize flag for parameter setting
c
      data init_paramset /0/
c
c  initialize flag for input of obsrot_file
c
      data init_obsrot /0/
c
c  precision parameters. values used here are suitable for a machine
c  with a nn1-bit matissa and nn2-bit exponent (as the univac 1100,
c  in single precision).
c
c  set in common/cprcns/ epsprc,epsufl,epsofl,eprufl,
c  where
c
c  epsprc is accuracy
c  epsufl is smallest real variable
c  epsofl is largest real variable
c  eprufl is sqrt(epsufl)
c
      epsprc=1.d-13
      epsufl=1.d-37
      epsofl=1.d0/epsufl
      eprufl=sqrt(epsufl)
c
      pi = 4.d0*atan(1.0d0)
c
c  ********************************************************************
c
c  Notes on equations and integration procedures:
c  ---------------------------------------------
c
c  for non-radial oscillations a factor x**(l-1) is taken out of
c  the solution out to the point x(nev1). the location of this
c  point is close to, but slightly below, the first turning point,
c  determined from the asymptotic relations using the trial frequency.
c  this is corrected for in the right hand side subroutine rhs
c  by setting the parameter el1 (passed in common/rhsdat/) to el - 1.
c  the unmodified equations are obtained by setting el1 = 0.
c
c  the equations are solved either by a shooting technique or by
c  a proper boundary value technique, as determined by the parameter
c  mdintg.
c
c  for mdintg = 1, 2 or 4 a shooting technique is used. the solutions
c  are integrated from the centre and from the surface, and matched
c  at the point x(nfit) closest to xfit, where xfit is input
c  in namelist /exec/. the condition that the matching be continuous
c  determines the eigenfrequency. the position of the matching point
c  should be chosen such that, as far as possible, the integration 
c  is not along a decreasing exponential in an evanescent region.
c
c  for mdintg = 1 the equations are approximated by centred differences.
c
c  for mdintg = 2 or 4 the equations are approximated on each mesh interval
c  by a system of constant coefficient equations, and these are
c  solved exactly in terms of harmonic and exponential functions.
c  thus modes of arbitrary high radial order can be computed
c  (see gabriel and noels, a & a, vol. 53, p. 149, 1976).
c  mdintg = 2 assumes a second-order system, i.e.
c  radial oscillations or non-radial oscillations in the cowling
c  approximation.
c  mdintg = 4 is for the full fourth-order system.
c
c  for mdintg = 3 the equations, together with some of, or all, the
c  boundary conditions are solved using a newton-raphson technique,
c  implemented in the subroutine nrkm. the eigenfrequency is
c  found by requiring that one of the boundary conditions, or a matching
c  condition, be satisfied.
c
c  the action taken depends on the value of xfit and hence nfit, 
c  relative to the end-points nibc and nn of the range of integration.
c
c  for nibc .lt. nfit .lt. nn nrkm is used with two regions, matched at
c  n = nfit. in this case all inner and outer boundary conditions are
c  satisfied, as well as the matching conditions at nfit for
c  y(1), y(3) and y(4). (in the cowling approximation, or in the
c  radial case, only continuity of y(1) is required). det is the
c  discontinuity in y(2) at the fitting point nfit.
c
c  for nfit .le. nibc or nfit .ge. nn nrkm is called with one region.
c  here either one of the inner or one of the outer boundary conditions
c  are left unsatisfied, and are used to iterate for sig.
c  for nfit .ge. nn det is set to the surface pressure 
c  boundary condition.
c  for nfit .le. nibc  det is set to the inner displacement boundary
c  condition.
c
c  the shooting technique is the fastest and probably the most flexible.
c  however it gives problems in the full non-radial case, when the
c  gravitational potential perturbation has almost no effect (i.e.
c  for high radial order or high degree), as
c  the matching determinant may then become almost degenerate. here
c  mdintg = 3 should be used, preferably with interior fitting.
c
c  in the cowling approximation the correction to the eigenfrequency,
c  and the perturbation in the gravitational potential (y(3,.)) is
c  calculated separately.
c
c  to obtain improved eigenfrequency the variational expression
c  may be used, as described by christensen-dalsgaard (monthly
c  notices, vol. 199, p. 735). as described there two formulations
c  are possible, one suitable for p modes and one suitable for
c  g modes. they are obtained by setting ivarf = 1 and 2, respectively.
c  for high-order modes it is essential that the correct formulation
c  is chosen.
c
c  notice that for radial modes only, ivarf may be set to 3.
c  this has the effect of using essentially the nonradial
c  formulation of the variational integral (in its p-mode form)
c  for radial modes, and may be useful to get strict consistency
c  between the calculation for radial and non-radial modes.
c  however this formulation apparently leads to fairly severe
c  cancellation and consequent loss of accuracy. thus it
c  probably cannot be recommended, although the modes for
c  which it may be of use must be investigated.
c
c  ********************************************************************
c
c  Notes on output to file:
c  ------------------------
c
c  In the current version, unit numbers are set up below, in the
c  source. This may be made more flexible later.
c
c  On d/s idslog a log is provided of modes where the iteration
c  failed or was questionable. (Note that inn most cases the 
c  "questionable" modes are probably fine.)
c  (Default value of idslog: 20; if idslog is not assigned to
c  a file in the input, the default name adipls-status.log is
c  used.)
c
c  on d/s idsgsm a complete summary of results for each mode is written,
c  defined by csummm(1-50). this is equivalenced to common/csumma/
c  above. the write is unformatted, in one record per mode.
c  see also separate notes (default value of idsgsm: 11)
c
c  on d/s idsssm a minimal summary of results for each mode is written,
c  defined by ssummm(1-7). the write is unformatted, in one record
c  per mode. see also separate notes (default value of idsssm: 15)
c
c  on d/s idsefn the eigenfunction is written when nfmode .ge. 1. 
c  The write is unformatted, the form depending on nfmode 
c  (default value of idsefn: 4)
c
c  nfmode = 1:
c     write(idsefn) csummm,nnw,(x(n),(y(i,n),i=1,4),(y(i,n),i=7,8),n=nw1,nn)
c
c  nfmode = 2:
c     First record: write(idsefn) nnw,(x(n),n=nw1,nn)
c     Subsequent records: write(idsefn) csummm,((y(i,n),i=1,2),n=nw1,nn)
c
c  nfmode = 3:
c     First record: write(idsefn) nnw,(x(n),n=nw1,nn)
c     Subsequent records: write(idsefn) csummm,((y(i,n),i=7,8),n=nw1,nn)
c
c  where nnw = nn - nw1 + 1 is the total number of points.
c  thus the extensive summary is written together with 
c  the eigenfunction.
c  see separate notes for the meaning of the quantities output.
c
c  on d/s idsrkr the rotational kernel is written for irotkr .ge. 1. 
c  The write is unformatted, of the form
c
c     write(idsrkr) csummm,nnw,(x(n),rotker(n),n=nw1,nn)
c  (default value of idsrkr: 12)
c
c  on d/s idsgkr the gamma1 kernel is written for igm1kr = 1. 
c  The write is unformatted, of the form
c
c     write(idsgkr) csummm,nnw,(x(n),gm1ker(n),n=nw1,nn)
c  (default value of idsgkr: 13)
c  ********************************************************************
c
c  notes on corrections and modifications:
c
c  23/1/1985: an error was found in s/r leq. this apparently
c     only affected the determination of maxima in s/r aramax,
c     and only on certain computers. as a result all values
c     of ximax and xmax (and of location and value of maximum
c     of energy normalized eigenfunction) computed at recku
c     since 1983 are wrong.
c
c  28/1/1985: allow possibility of no region where x**(l-1) is
c    taken out. this is required when sigma**2 is very small, and
c    the inner evanescent region consequently only occupies a few
c    mesh points.
c    to correct for it, a flag is passed to s/r setbcs to determine
c    whether to set y(i,.) or x**(-l+1)*y(i,.).
c
c  2/4/1985: the arrays ds and ds1 used in determining
c     matching determinant have been made double precision,
c     and the call of leqdet has been changed to dleq.
c     this is an attempt to get round cancellation errors when
c     solving for high-order modes in the full case.
c     **** but see below, note of 11/7/1985.
c
c     in addition, the calculation of the matching coefficients
c     has been modified, to use the original equations.
c     the parameters imstsl, imissl and imjssl have been added
c     to /exec/ to provide flexibility as to how this is done.
c
c  3/4/1985: the definition of the printed 'calcellation factor'
c     has been changed from being the minimum to being the mean
c     of the quantity abs(yy1+yy2)/(abs(yy1)+abs(yy2)), where
c     yy1 and yy2 are the two, hopefully, independent solutions.
c
c  30/5/1985: parameter irsord added to namelist /exec/.
c     if irsord = 1 and icow = 0, increment order for l = 1 modes, 
c     original order = 0 or 1, by 1. this is suitable for normal
c     models of  p r e s e n t  sun.
c
c  11/7/1985: the matching determinant and fitting have been restored
c     to single precision.
c
c  12/7/1985: sigc, in grand summary, is now set to sig ( = sigma**2)
c     when cowling approximation is not used.
c
c     frfit replaced by xfit in namelist /exec/.
c
c  17/8/1985: modifications developed in cambridge in the spring of 1984
c     incorporated in this version. specifically these involve
c    - the option of having ivarf = 3 for radial modes in s/r varfrq
c      (see notes above).
c    - the option of using a short summary to set trial solution.
c    - the option of producing the gamma1-kernel (for igm1kr = 1).
c     in addition the programme now produces a short summary 
c     (on d/s 15).
c
c  18/8/1985: change meaning of ddsol when integrating with relaxation
c     technique and interior matching.
c     allow extra iteration if eigenfunction is not converged when
c     solution is.
c
c  19/3/1986: s/r varfrq modified to use tabulated differentiation and
c     integration coefficients in s/r derivk and s/r vintk.
c
c  25/3/1986:  meaning of parameter fctsbc in s/r setbcs changed.
c     previously used for an approximation to beta+, beta- mixing.
c     now enables transition from delta p = 0 to p prime = 0.
c
c  12/8/87: version from RECKU and version from Boulder (1985) 
c     merged.
c
c  13/8/87: begin modifications to make programme run under UNIX.
c     Unit number for formatted summary changed from 9 to 16.
c     Introduction of output units istdou and istdpr.
c     Argument list of s/r dmpsum changed.
c
c  19/2/88: include cowling approximation, and variable lambda, 
c     for the radial case.
c
c  23/3/88: include consistent setting of dimensional values,
c     from gravitational constant cgrav. cgrav is read in in 
c     namelist, if cntrd contains cst
c
c  30/11/88: include option ibotbc = 2, corresponding
c     to bottom condition d(xir)/dx = 0, matching Stein & Nordlund.
c
c  13/2/89: include option for Richardson extrapolation.
c     Note that the extrapolated cyclic frequency is set in
c     cs(37). All other quantities output are not extrapolated.
c     When Cowling approximation is set with icow = 2, Richardson
c     extrapolation is based on corrected eigenfrequencies.
c     When icow = 3, uncorrected frequencies are used instead.
c     This is flagged by setting mod(icase,10) = 2, instead of 1.
c
c  10/5/89: modified to test for same mesh etc. before using tabulated
c     coefficients in s/r derivk and vintk in s/r varfrq.
c
c  11/5/89: modified to allow radial oscillations of truncated model
c     without Cowling approximation. 
c
c  11/5/89: modified to calculate correction to frequency
c     in Cowling approximation for truncated model.
c
c  18/2/90: correct setting of turbulent pressure ratio in aa(6,.)
c     Previously aa(6,n) was set to 1/(x**3*aa(1,n)) for x .ge. 0.999.
c     Now aa(6,n) is set to 1, unless iturpr = 1
c
c  12/7/90: correct setting gravitational potential correction in 
c     s/r gravpo on submesh during Richardson extrapolation, for
c     model with singular surface.
c
c  16/7/90: change value of epsprc from 1.e-6 to 1.e-10, more suitable for
c     double precision version of programme.
c     Add check for same order in Richardson extrapolation.
c
c  18/8/90: introduce options itrsig = 6 or 7 to read trial frequency
c     from file of observed data.
c
c  8/8/91: add parameter dsigre, as decrement for trial sig when doing
c     Richardson extrapolation (previously dfsig was used also in
c     this capacity).
c
c  16/12/93: major revision, to improve treatment of singular surface:
c     - Include additional terms in expansion in setbcs
c     - Allow integration of inhomogeneous equation for 
c       0 lt npol lt 1, in rhs, shtint and nrkint (note:
c       this should be looked into; strange statement
c       anres = 0 after call of setbcs).
c     - Fix treatment of near-surface contribution in varfrq
c
c  23/12/93: reset constants in period-to-frequency conversions
c     to be consistent double precision.
c
c  21/7/94:  Add diagnostic output to file 'adipls-status.log'
c     (unit idslog; set to 20) in cases of no or questionable
c     convergence.
c  21/7/94:  Add test on eigenfunction difference if Richardson
c     extrapolation encounters different orders.
c  22/7/94:  Add automatic change of fitting point if
c     iteration fails to converge. Controlled by new
c     parameter nftmax in int group
c  27/10/94: Fix resetting of nwmod during frequency scan.
c     Error affected writing of eigenfunctions for nfmode .gt. 1
c  19/1/95: Introduce possibility of treating turbulent pressure
c     as in CR formulation, flagged by iturpr = 2
c     Move turbulent pressure ratio from aa(6,.) to aa(10,.)
c     Add 1 000 000*iturpr to icase.
c  10/8/96: Introduce use of istsig as flag for choice of frequency
c     limits when scanning in frequency for iscan .gt. 1
c
c  17/2/98: Set iturpr to 8 in readml, if model contains g/(g tilde)
c     in aa(6,.) (as is the case with models including the
c     spherically symmetric effect of rotation).
c
c  17/7/02: Add option imdmod = -1, to add point near centre if 
c     required for Richardson extrapolation. REMOVED 14/2/06, as
c     being apparently inactive.
c
c  29/3/05: Add option mdintg = 5, to use fourth-order integrating method
c     in shooting technique.
c
c  19/12/07: Increase size of arrays in common/cofile/
c
c  28/2/08: Increase size of internal storage for frequency results in
c     common/cobs_param/, and check to avoid overwriting if more modes
c     are computed.
c     Introduce option inomd1 .ge. 2 to stop scan if frequency increases
c     beyond acoustical cut-off frequency (or for other diagnostics from 
c     setbcs in sigscn.
c
c  24/7/09: Allow output of unconverged eigenfunctions during simple scan,
c     for setting core phase, flagged by setting nfmode.
c     Also allow output of unconverged eigenfunctions for itmax = 0
c     (typically when reading already computed, or observed, frequencies),
c     flagged by setting nfmode/10 = 1.
c
c  6/5/10: Include option irsord = 20 to use the Takata (2006; ESA SP-624)
c     scheme for computing the order of dipolar modes.
c
c  ********************************************************************
c
      write(istdou,'(//
     *  '' --------------------------------------------------''/
     *  '' Entering adipls with parameter passing''/
     *  '' i_paramset ='',i2/)') i_paramset
      if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,'(//
     *  '' --------------------------------------------------''/
     *  '' Entering adipls with parameter passing''/
     *  '' i_paramset ='',i2/)') i_paramset
c..      write(0,*) 'Enter adipls. istdou, istdpr =',istdou, istdpr
c
c  initialize diagnostic flags
c
      kdgbcs=0
c
      ierr_param=-1
c
c  test for skipping to continuing run with parameter setting
c
      if(i_paramset.eq.1.and.init_paramset.ne.0) go to 10000
c
c  as initilization zero all elements in csummm
c
      call zero(csummm,50)
c
      init_paramset=0
c
c  ********************************************************************
c
c  defaults in exec
c  ****************
c
c  cntrd: string determining which control fields are read.
c      contains some or all of the following 
c      dsn: change output data set designators
c      mod: read model controls
c      osc: read controls for mode selection
c      cst: read new values of fundamental constants (currently
c           just gravitational constant)
c      int: read controls for type of equations or integration
c           procedure
c      out: read controls for output
c      dgn: read controls for diagnostics
c      fields may be separated by, e.g., "." (but  n o t comma
c      or blank).
      cntrd='mod.osc.int.out.dgn'
c
c  output file unit numbers
c  ************************
c
c  idsgsm: grand summary output
      idsgsm=11
c  idsssm: short summary output
      idsssm=15
c  idsefn: eigenfunction output
      idsefn=4
c  idsrkr: output of kernels for spherical rotation
c  if idsrkr = 0, no file output is produced
      idsrkr=12
c  idsgkr: gamma-kernel output
      idsgkr=13
c  idslog: logfile
      idslog=20
c
c  model controls
c  **************
c
c  ifind, xmod: determines which model is used
c  ifind .lt.0: do not read new model, except if no model has been
c               read so far.
c  ifind   = 0: rewind data set before reading model
c  ifind   = 1: read next model on data set
c  ifind   = 2: read model no. xmod on dataset. xmod is currently
c               assumed to be integer-valued. non-integer xmod may
c               be implemented later, with interpolation between
c               models.
c  For ifind = -2, possibly reset model with s/r modmod or set
c  new truncation, even if no new model was read.
      ifind=-1
      xmod=0
c  imlds: models are read from d/s imlds
      imlds=2
c  in: after reading model, take every in-th point
      in=1
c  nprmod: if nprmod .gt. 0 print the model at nprmod points.
      nprmod=0
c  xtrnct: if xtrnct .gt. 0 truncate the model at fractional radius
c          xtrnct. This is only implemented for radial
c          oscillations, or in the Cowling approximation. Thus
c          if l .gt. 0 and the model is truncated, icow = 2 is
c          forced.
      xtrnct=0
c  ntrnsf: if ntrnsf .gt. 1 truncate model ntrnsf points from the
c          surface.
      ntrnsf=0
c  imdmod: when imdmod .ne. 0 call s/r modmod, which may be 
c          user specified to modify the model
c  imdmod = -1: With Richardson extrapolation, add a mesh point near
c  lower boundary if required to ensure even number of mesh intervals.
c  imdmod = -1 option REMOVED 14/2/06.
      imdmod=0
c
c  controls for mode selection:
c  ***************************
c
c  el, nsel, els1, dels: controls for determining the value 
c      of the degree
c  el: when nsel .le. 0 (and itrsig .ne. 2, see below) use el as input
c      in exec.
      el=1
c  nsel: when nsel .ge. 1 step through el = els1 + i*dels,
c  i = 0, ..., nsel-1. when iscan .le. 1 (see below) sig is given
c  the initial value in sig1, and is incremented by dfsig for each
c  new value of el.
      nsel=0
      els1=10
      dels=1
c  dfsig1, dfsig2, nsig1, nsig2: allows resetting of the sig1 and sig2
c  during step in el and scan in sig. For each new value of el,
c  sig1 and sig2 are incremented as determined by (nsig1, dfsig1),
c  (nsig2, dfsig2), following the definition of nsig and dfsig below.
c  Note: currently dfsig1, dfsig2 refer to the scan limits in
c  squared dimensionless frequency, regardless of the value of istsig.
c
      dfsig1=0
      dfsig2=0
      nsig1=1
      nsig2=1
c
c  sig1, sig2, itrsig, istsig, itrds, dfsig, nsig: 
c       determine trial frequency.
c  itrsig = 0: trial frequency taken from sig1
c  itrsig = 1: sig is found from previous value sigp and dfsig.
c      meaning of dfsig depends on nsig:
c       nsig = 1: sig = sigp + dfsig
c       nsig = 2: dfsig is increment in frequency (i.e. in sqrt(sig))
c       nsig = 3: dfsig is increment in 1/(frequency) 
c       (i.e. in 1/sqrt(sig))
c       nsig = 4: uniform in 1/(frequency) at low frequency and
c       in frequency at high frequency.
c       (note that itrsig and nsig has a special meaning 
c       when iscan .gt. 1, see below).
c  itrsig = 2 or 3: find trial frequency from values on
c       double precision grand summary residing on d/s itrds.
c       if itrsig = 2 mode no inomd1 + jstsig - 1 is taken,
c       and el is reset to the value for this mode.
c       if itrsig = 3 the mode with the value of el input in exec
c       and of order inomd1 + jstsig - 1 is used.
c       here jstsig is stepped through 1,...,istsig.
c  **** The following appears not to be implemented!
c  **** when itrsig .gt. 1 and sig2 .gt. sig1 only values of sig between
c  **** sig1 and sig2 are used. when in addition itrsig = 3 the
c  **** modes in the given range and with the given el are used,
c  **** starting from order inomd1 and up to order inomd1 + istsig - 1.
c  itrsig = 4 or 5: find trial frequency from values on
c       double precision short summary residing on d/s itrds.
c       otherwise itrsig = 4 and 5 corresponds to itrsig = 2 and 3 
c       above.
c  itrsig = 6 or 7: find trial frequency from values of 
c       cyclic frequencies (in microHz) from file residing on d/s itrds.
c       otherwise itrsig = 6 and 7 corresponds to itrsig = 2 and 3 
c       above.
c  itrsig = -2 to -5: find trial frequency from single precision
c       grand or short summary, as above.
c  other values of itrsig are currently not implemented.
c  Note that istsig and inomd1 have special meaning when iscan .gt. 1;
c  see below.)
c
      itrsig=0
      sig1=10
      sig2=0
      dfsig=0.d0
      nsig=1
      istsig=1
      inomd1=1
      itrds=10
c  iscan: flag for scan in sig.
c  when iscan .gt. 1 step in sig with iscan steps between 
c  frequency limits sig1 and sig2.
c  When istsig .le. 1, sig1 and sig2 are limits in squared dimensionless
c  frequency. Otherwise, sig1 and sig2 are taken as limits in
c  cyclic frequency, in mHz.
c  step is uniform in squared frequency, frequency or period, for
c  nsig = 1, 2, and 3, respectively.
c  when itrsig = 1 check for change of sign in the matching determinant
c  and iterate for eigenfrequency at each change of sign.
c  if inomd1 .ge. 2 stop scan after diagnostics from integration.
c  (Implemented 28/2/08; currently only for diagnostics from setbcs).
c

      iscan=1
c
c  notes on preference in determining trial el and sig: 
c     nsel .ge. 1 may be combined with iscan .gt. 1 to produce 
c     a scan between sig1 and sig2 in sig for each 
c     el = els1, els1+dels, ..., els1+(nsel-1)*dels.
c
c     test on iscan takes place before test on itrsig.
c
c     nsel .gt. 1 cannot currently be combined with itrsig .gt. 1.
c     when itrsig .gt. 1 (or itrsig = 1 for iscan .le. 1) the
c     value of the degree in el is always used.
c  eltrw1, eltrw2, sgtrw1, sgtrw2: windows applied to trial degree or sig.
c  Only takes effect if eltrw1 .le. eltrw2, 
c  respectively sgtrw1 .le. sgtrw2
      eltrw1=0
      eltrw2=-1
      sgtrw1=0
      sgtrw2=-1
c
c  controls for rotational solution
c  ********************************
c
c  irotsl: if irotsl = 1, calculate solution including first-order
c  rotational effects as in Soufi et al., for azimuthal order em.
      irotsl=0
      em=0
c  nsem, ems1, dems: Controls step in azimuthal order, for
c  irotsl = 1.
c  nsem = -1: step with step dems, from -el to el
c  nsem .gt. 0: take nsem steps, starting from ems1, with step dems
      nsem=0
      ems1=0
      dems=0
c
c  fundamental constants:
c  *********************
c
c  cgrav: value of gravitational constant (default is value hardcoded
c     before 23/3/88)
      cgrav = 6.6732d-8
c
c  controls for equations and integration:
c  **************************************
c
c  iplneq: when iplneq = 1 use equations for plane-parallel layer
c          (note that model coefficients should then also correspond
c           to plane-parallel case).
      iplneq=0
c  iturpr: flag for turbulent pressure in equilibrium model
c     when iturpr = 1, set aa(10,.) to turbulent pressure ratio in
c     subroutine readml. otherwise aa(10,.) = 1 always.
c     When iturpr = 2, assume that model is given with 6 variables,
c     the sixth being correction to A2.
c     If model contains g/(g tilde) in aa(6,.), iturpr is set to 8
c     in s/r readml.
      iturpr=0
c  icow: flag for cowling approximation
c  icow = 0: solve full equations
c       = 1: solve equations in cowling approximation and go back
c            and solve full equations with cowling result as trial
c       = 2: solve equations in cowling approximation only. Also sets
c            frequencies corrected by perturbation technique.
c       = 3: As icow = 2, except that Richardson extrapolation of
c            cyclic frequency is based on uncorrected eigenfrequency.
      icow=0
c  alb: fudge factor in poisson's equation. should be 1, except
c       when studying gradual transition between cowling approximation
c       and full case.
      alb=1.d0
c
c  istsbc, fctsbc: determines surface pressure boundary condition.
c  istsbc = 1: find condition by matching to exponentially
c     decaying solution in an isothermal atmosphere matched to
c     the outermost mesh point. this assumes that the frequency is
c     below the acoustical cut-off frequency at that point.
c     otherwise a message is printed and the condition for
c     istsbc = 0 is used.
c  istsbc = 0: use (lagrangian pressure perturbation) = 0 on surface.
c     fctsbc: determines condition when istsbc = 0.
c     for fctsbc = 0 use delta p = 0.
c     for fctsbc = 1 use p prime = 0. note that all intermediate cases
c     are allowed.
c  istsbc = 9: use delta r = 0 on surface.
      istsbc=1
      fctsbc=0
c  ibotbc, fcttbc: determines bottom boundary condition in truncated
c     model.
c  ibotbc = 0: set relation between y1 and y2 at bottom to isolate
c     solution that decreases exponentially towards the interior
c     this assumes that the bottom is in an evanescent region for
c     the l-value and frequency used. 
c     otherwise the condition corresponding to ibotbc = 1 is used.
c  ibotbc = 1: set y1 = alpha*(1-fcttbc), y2 = alpha*fcttbc at bottom,
c     where alpha is an arbitrary proportionality constant.
c  ibotbc = 2: use bottom condition d(xir)/dx = 0, to match 
c     Stein & Nordlund calculation.
c  ibotbc = 3: use bottom condition d/dx(xir/r) = 0, as suggested
c     by JOP.
      ibotbc=0
      fcttbc=0
c
c  mdintg: determines type of integration used
c  mdintg = 1: use centred differences, integrate from centre and
c              surface and get frequency from matching
c  mdintg = 2: use constant coefficient integration on each mesh
c              interval, for second-order systems,
c              i.e. radial or cowling approximation. if mdintg = 2
c              and l .gt. 0 mdintg = 4 is forced.
c  mdintg = 3: use nrkm. find frequency by iterating on interior
c              boundary condition for xfit = 0, outer boundary
c              condition for xfit = 1, or at an internal point for
c              0 .lt. xfit .lt. 1.
c  mdintg = 4: use constant coefficient integration on each mesh
c              interval, for the full fourth-order system.
c              in the case of a second-order system,
c              i.e. radial or cowling approximation, 
c              mdintg is reset to 2.
c  mdintg = 5: use fourth-order scheme (Cash & Moore 1980; BIT 20, 44),
c              integrate from centre and surface and get frequency 
c              from matching
      mdintg=1
c  iriche: if iriche = 1, Richardson extrapolation is used to
c     improve the computed frequency (but not the other quantities),
c     by solving initially the equations on every second meshpoint.
c     Note: trial frequency is decreased by dsigre, below.
      iriche=0
c  xfit: match solution at point xfit. 
c     Note that if irsevn = 2 the setting of matching point by xfit
c     may be overridden by the location of the evanescent 
c     region.
c     For xfit = -1 match solution at edge of inner evanescent region.
      xfit=0.5d0
c  fcnorm: for mdintg = 3 and xfit = 0 or 1, normalize boundary 
c          condition by solution at point fcnorm*nn.
      fcnorm=0.5d0
c  eps: convergence of sig assumed when relative change in sig between
c       two iterations is less than eps.
      eps=2*epsprc
c  epssol: when mdintg = 1, 2 or 4 convergence criterion for 
c          eigenfunction is assumed to be that relative discontinuity 
c          at matching point is less than epssol. 
c          for mdintg = 3 convergence criterion is
c          that the mean relative change in the eigenfunction between
c          two iterations is less than epssol.
      epssol=1.d-5
c  itmax: maximum number of iterations
c         when itmax = 0 just integrates once and, when 
c         mdintg = 1, 2 or 4, tests the continuity of the eigenfunction.
c         this is useful for re-computing eigenfunctions when 
c         the eigenfrequency is known.
      itmax=8
c  dsigre: when using Richardson extrapolation, make relative change in
c  trial sig by dsigre before iteration on full mesh
c  In addition, dsigre .le. -1 is used to flag for frequency resetting
c  if orders do not agree with Richardson extrapolation.
      dsigre=0
c  fsig: in the secant iteration for sig, the second trial value is
c        (1+fsig)*(the first trial value)
      fsig=0.001d0
c  dsigmx: the relative change in sig during the iteration is limited
c          to be less than dsigmx in absolute value.
      dsigmx=0.1d0
c  irsevn: for irsevn .ge. 1, when iterating for eigenfrequency during
c     scan (i.e. for iscan .gt. 1, itrsig = 1), reset transition
c     point between modified and standard equations (at boundary
c     of evanescent inner region for the trial sig) before 
c     each iteration.
c     For irsevn = 2, reset fitting point if it is deeper
c     than the evanescent transition (default is to shift
c     evanescent transition to fitting point).
c     For irsenv = 3, in addition do not rescale inner boundary values
c     with x**(l-1), even if that factor is not taken out in the
c     solution. This is useful to ensure a smooth scanning, when the
c     onset of the rescaling takes place during the scan.
c     for irsevn = -1, do not use transformation in evanescent region.
      irsevn=2
c  xmnevn: search for evanescent transition only outside x = xmnevn
      xmnevn=0
c  nftmax: controls resetting of fitting point, if .gt. 1.
c     nftmax = nftmx1 + 10*nftrcs.
c     attempt changing fitting point up to nftmx1 times, 
c     if iteration does not converge
c     or the correct order is not obtained for itsord = 1
c     mode of resetting is controlled by nftrcs. For nftrcs = 0,
c     base resetting on mesh-point number, and for nftrcs = 1,
c     base resetting on mesh in x.
      nftmax=3
c  itsord: if itsord = 1, test that computed order corresponds to order
c     of trial mode. Otherwise try resetting nfit
c     if itsord = -1, drop test of order in Richardson extrapolation
c     (as a rather desperate measure)
      itsord = 0
c
c  controls for output:
c  *******************
c
c  istdpr: unit number for printed output. default set in block data
c      subprogram blstio
c
c  nout: when nout .gt. 0 print solution at nout points
      nout=50
c  nprcen: if nprcen .gt. 1 print solution at all the nprcen points
c          closest to the centre.
      nprcen=0
c  irsord: if irsord gt 0 and irsord le 10, and icow = 0, increment 
c  order for l = 1 modes with Scuflaire order between 0 and irsord 
c     irsord = 1 is appropriate in normal models of present sun.
c  for iabs(irsord) = 11, use instead Lee classification in all cases
c  for iabs(irsord) = 12, use instead Lee classification only for l = 1
c  for irsord = 20 use the Takata (2006; ESA SP-624) scheme for computing 
c     the order of modes with l = 1.
      irsord=0
c  iekinr: For iekinr = 0, normalize mode energy with surface
c  vertical displacement. For iekinr = 1, normalize with total
c  photospheric displacement
      iekinr=0
c  iper, ivarf, kvarf, npvarf: controls for calculation of the
c  variational frequency.
c  iper: when iper = 1 calculate variational frequency.
      iper=0
c  ivarf = 1: use p mode formulation.
c  ivarf = 2: use g mode formulation.
c  ivarf = 3: (only applicable for radial modes) use non-radial
c              formulation (see note above).
      ivarf=1
c  kvarf: numerical differentiation and integration in calculation
c         of variational frequency uses polynomials of degree 2*kvarf
      kvarf=2
c  npvarf: for npvarf .gt. 0 print integrands and integrals in
c          variational calculation at npvarf points.
      npvarf=0
c  nfmode: Control of eigenfunction output. 
c  On input nfmode = nfmod0 + 10*nfmod1.
c  Internally, nfmode is set to nfmod0.
c  for nfmode .ge. 1 write eigenfunctions to file on d/s 4.
c  The variables output and the format depend on whether 
c  nfmode = 1,2 or 3, as described above.
c  nfmod1 = 1 flags for output of possibly discontinuous eigenfunctions
c  when itmax = 0.
      nfmode=0
c  irotkr, nprtkr: for irotkr .ge. 1 calculate rotational kernels.
c    For irotkr = 11 in addition compute rotational splitting from
c    angular velocity in common/comgrp/. 
c update: 18/05/07 KDB
c    For irotkr = 21 in addition compute second-order rotational splitting
c    from angular velocity in common/comgrp/.
c end update
c    When the model is passed
c    as an argument to this routine (i.e., with i_inout .eq. 0)
c    the angular velocity must be set in common/comgrp/ before the
c    call of adipls; in this case it will often be transferred
c    from evolution code but could also be set by an external call
c    of s/r set_rotation. When the model is read in during the calculation,
c    set_rotation is called after the model read. In this case the
c    parameter input file should contain the file name of the
c    file describing the angular velocity (if required in set_rotation)
c    and the file name for `obs' output including rotational results,
c    after the parameter input.
c    If in addition nprtkr .gt. 1 rotational kernel is printed
c    at nprtkr points.
      irotkr=0
      nprtkr=0
c  igm1kr, npgmkr: for igm1kr = 1 calculate gamma1 kernels.
c    if in addition npgmkr .gt. 1 gamma1 kernel is printed
c    at npgmkr points.
      igm1kr=0
      npgmkr=50
c  ispcpr: if ispcpr .ne. 0 special output may be produced by
c          call of user-supplied routine spcout.
c  Also used to define type of frequency stored in common/cobs_param/
c  (see s/r setobs_st for details):
c  ispcpr = 1: variational frequency.
c  ispcpr = 4: from eigenfrequency in cs(20).
c              Note that this allows setting Cowling
c              approximation frequency.
c  ispcpr = 5: from Richardson extrapolation frequency
c              in cs(37), if this is set.
c              Otherwise variational frequency is used.
c  ispcpr = 6: from (possibly corrected) eigenfrequency in cs(21).
      ispcpr=0
c  icaswn, sigwn1, sigwn2, frqwn1, frqwn2, iorwn1, iorwn2, frlwn1, frlwn2:
c  controls which modes are output to file.
c  For windows defined by pairs of parameters (e.g. sigwn1, sigwn2)
c  windowing is applied only if first is .le. second 
c  (e.g. sigwn1 .le. sigwn2). Defaults are now windowing.
c    icaswn: if .ge. 0, select only modes with icase = icaswn
c    sigwn1, sigwn2: Window in sigma**2
c    frqwn1, frqwn2: Window in cyclic frequency nu (in microHz) 
c    iorwn1, iorwn2: Window in mode order
c    frlwn1, frlwn2: Window in nu/(l + 1/2) (nu in microHz)
c  Defaults:
      icaswn=-1
      sigwn1=0
      sigwn2=-1
      frqwn1=0
      frqwn2=-1
      iorwn1=0
      iorwn2=-1
      frlwn1=0
      frlwn2=-1
c
c  controls for diagnostics:
c  ************************
c
c  itssol: when itssol = 1 test solution with nrkm right hand side
c          routine. for idgtss = 1 additional details about solution
c          etc. are printed.
      itssol=0
      idgtss=0
c  moddet, iprdet: flags for modifying and printing matching 
c     determinant.
      moddet=0
      iprdet=0
c  npout: when finding eigenfrequency by matching and npout .gt. 0
c         print the separate solutions at npout points.
      npout=0
c  imstsl, imissl, imjssl (added 2/4/1985):
c     parameters controlling the determination of the matching
c     coefficients in the full case.
c     the equation no. imissl is ignored, and the coefficient no.
c     imjssl is given a fixed value.
c     if imstsl .ne. 1, imissl and imjssl are taken as input.
c     otherwise imissl and imjssl are determined to minimize
c     the error.
      imstsl=1
      imissl=2
      imjssl=2
c  idgnrk: determines level of diagnostics in nrkint solution,
c     for mdintg = 3
      idgnrk=0
c
c
c  ********************************************************************
c
c  notes on indices for model and solution:
c  ---------------------------------------
c
c  the model is set up (in call of s/r readml) in x(n), aa(i,n),
c  i = 1, ..., 6, n = 1, ..., nn.
c  if the model is truncated at the surface (for ntrnsf .gt. 0) nn is
c  reset in the present subroutine.
c  the solution is set into y(i,n), i = 1, ..., ii, n = nw1, ..., nn.
c  here nw1 = 1 if the model is not truncated in the interior,
c  otherwise nw1 = ntrnct. note that nnw = nn - nw1 + 1 is used for
c  the total number of points in the solution.
c
c  the inner boundary conditions are applied at the point nibc.
c  when nw1 = 1, and x(1) is at the centre of the model, nibc = 2.
c  otherwise nibc = nw1.
c
c  when integrating the equations with matching 
c  (i.e. for mdintg = 1, 2 or 4 ) the solution outside 
c  the matching point is temporarily shifted by one,
c  so that y(i,n+1) is the solution at point x(n). the outer boundary
c  condition is applied on y(i,ny) where ny = nn+1 if the surface is
c  non-singular, and ny = nn if the surface is singular. ny (and nd,
c  used for storage information for the integration routine genint)
c  are set in s/r setbcs.
c
c  ********************************************************************
c
c  note: the setting and resetting of the case number icase is unclear
c        for some cases, and has to be checked and cleared up.
c
c  ********************************************************************
c
c  ....................................................................
c
      sig=0
c
      ntrncp=-1
      nnp=0
c
      idoutd=11
c
c  store input value of nnww in nnwwin
c
      nnwwin=nnww
c
      frsrun=.true.
c
c  initialize previous kvarf
c
      kvarfp=-1
c
c              ----------------------------------------
c
c  set file names, and fixed output files
c
      write(istdou,*) ' Files needed:'
      write(istdou,*) ' Model input on unit imlds (default 2)'
      write(istdou,*) 
     *           ' Trial frequency input on unit itrds (default 10).',
     *           ' Optional.'
      write(istdou,*) ' Grand summary on unit idsgsm (default 11).'
      write(istdou,*) ' Short summary on unit idsssm (default 15).'
      write(istdou,*) 
     *  ' Eigenfunction on unit idsefn (default 4) (optional).'
      write(istdou,*) 
     *  ' Rotational kernel on unit idsrkr (default 12) (optional).'
      write(istdou,*) 
     *  ' Gamma1 kernel on unit idsgkr (default 13) (optional).'
      write(istdou,*) ' Printed output on unit istdpr, if not default'
c
      call ofiles
c
      nfmesh=1
c                     ----------------------------
c
c  read exec
c
    6 xtrncp=xtrnc1
c
      write(istdou,
     *  '(//'' ----------------------------------------------''//)')
c
      write(istdou,*) ' cntrd?'
      write(istdou,*) cntrd
      read(istdin,'(a)',end=85,err=90) cntrd
c
c  decode cntrd
c
      ircdsn=index(cntrd,'dsn')
      ircmod=index(cntrd,'mod')
      ircosc=index(cntrd,'osc')
      ircrot=index(cntrd,'rot')
      irccst=index(cntrd,'cst')
      ircint=index(cntrd,'int')
      ircout=index(cntrd,'out')
      ircdgn=index(cntrd,'dgn')
c
c  test whether any input has been selected
c
      if(ircmod+ircosc+irccst+ircint+ircout+ircdgn.eq.0) then
	write(istdou,103) cntrd
	go to 90
      end if
c
      if(ircdsn.gt.0) then
        write(istdou,*) 
     *    'idsgsm, idsssm, idsefn, idsrkr, idsgkr, idslog?'
        write(istdou,*) idsgsm, idsssm, idsefn, idsrkr, idsgkr, idslog
        read(istdin,*,end=85,err=90) 
     *    idsgsm, idsssm, idsefn, idsrkr, idsgkr, idslog
      end if
c
      if(ircmod.gt.0) then
c
        write(istdou,*) ' ifind,xmod,imlds,in,nprmod?'
        write(istdou,*) ifind,xmod,imlds,in,nprmod
        read(istdin,*,end=85,err=90) ifind,xmod,imlds,in,nprmod
        write(istdou,*) ' xtrnct,ntrnsf,imdmod?'
        write(istdou,*) xtrnct,ntrnsf,imdmod
        read(istdin,*,end=85,err=90) xtrnct,ntrnsf,imdmod
c
      end if
c
      if(ircosc.gt.0) then
c
        write(istdou,*) ' el,nsel,els1,dels,dfsig1,dfsig2,nsig1,nsig2?'
        write(istdou,*) el,nsel,els1,dels,dfsig1,dfsig2,nsig1,nsig2
        read(istdin,*,end=85,err=90) 
     *    el,nsel,els1,dels,dfsig1,dfsig2,nsig1,nsig2
        write(istdou,*) ' itrsig,sig1,istsig,inomd1,itrds?'
        write(istdou,*) itrsig,sig1,istsig,inomd1,itrds
        read(istdin,*,end=85,err=90) itrsig,sig1,istsig,inomd1,itrds
        write(istdou,*) ' dfsig,nsig,iscan,sig2?'
        write(istdou,*) dfsig,nsig,iscan,sig2
        read(istdin,*,end=85,err=90) dfsig,nsig,iscan,sig2
        write(istdou,*) ' eltrw1, eltrw2, sgtrw1, sgtrw2?'
        write(istdou,*) eltrw1, eltrw2, sgtrw1, sgtrw2
        read(istdin,*,end=85,err=90) eltrw1, eltrw2, sgtrw1, sgtrw2
c
      end if
c
      if(ircrot.gt.0) then
c
        write(istdou,*) ' irotsl, em, nsem, ems1, dems?'
        write(istdou,*)   irotsl, em, nsem, ems1, dems
        read(istdin,*,end=85,err=90) irotsl, em, nsem, ems1, dems
c
      end if
c
      if(irccst.gt.0) then
c
        write(istdou,*) ' cgrav?'
        write(istdou,*) cgrav
        read(istdin,*,end=85,err=90) cgrav
c
      end if
c
      if(ircint.gt.0) then
c
        write(istdou,*) ' iplneq,iturpr,icow,alb?'
        write(istdou,*) iplneq,iturpr,icow,alb
        read(istdin,*,end=85,err=90) iplneq,iturpr,icow,alb
        write(istdou,*) ' istsbc,fctsbc,ibotbc,fcttbc?'
        write(istdou,*) istsbc,fctsbc,ibotbc,fcttbc
        read(istdin,*,end=85,err=90) istsbc,fctsbc,ibotbc,fcttbc
        write(istdou,*) 
     *    ' mdintg,iriche,xfit,fcnorm,eps,epssol,itmax,dsigre?'
        write(istdou,*) 
     *    mdintg,iriche,xfit,fcnorm,eps,epssol,itmax,dsigre
        read(istdin,*,end=85,err=90) 
     *    mdintg,iriche,xfit,fcnorm,eps,epssol,itmax,dsigre
        write(istdou,*) ' fsig,dsigmx,irsevn,xmnevn,nftmax,itsord?'
        write(istdou,*) fsig,dsigmx,irsevn,xmnevn,nftmax,itsord
        read(istdin,*,end=85,err=90) 
     *    fsig,dsigmx,irsevn,xmnevn,nftmax,itsord
c
      end if
c
      if(ircout.gt.0) then
c
        write(istdou,*) ' istdpr,nout,nprcen,irsord,iekinr?'
        write(istdou,*) istdpr,nout,nprcen,irsord,iekinr
        read(istdin,*,end=85,err=90) istdpr,nout,nprcen,irsord,iekinr
        write(istdou,*) ' iper,ivarf,kvarf,npvarf,nfmode?'
        write(istdou,*) iper,ivarf,kvarf,npvarf,nfmode
        read(istdin,*,end=85,err=90) iper,ivarf,kvarf,npvarf,nfmode
        write(istdou,*) ' irotkr,nprtkr,igm1kr,npgmkr,ispcpr?'
        write(istdou,*) irotkr,nprtkr,igm1kr,npgmkr,ispcpr
        read(istdin,*,end=85,err=90) irotkr,nprtkr,igm1kr,npgmkr,ispcpr
        write(istdou,*) 'icaswn, sigwn1, sigwn2, frqwn1, frqwn2, ',
     *                'iorwn1, iorwn2, frlwn1, frlwn2?'
        write(istdou,*) icaswn, sigwn1, sigwn2, frqwn1, frqwn2, 
     *                iorwn1, iorwn2, frlwn1, frlwn2
        read(istdin,*,end=85,err=90)
     *  	      icaswn, sigwn1, sigwn2, frqwn1, frqwn2, 
     *                iorwn1, iorwn2, frlwn1, frlwn2
c
      end if
c
      if(ircdgn.gt.0) then
c
        write(istdou,*) ' itssol,idgtss,moddet,iprdet,npout?'
        write(istdou,*) itssol,idgtss,moddet,iprdet,npout
        read(istdin,*,end=85,err=90) 
     *    itssol,idgtss,moddet,iprdet,npout
        write(istdou,*) ' imstsl,imissl,imjssl,idgnrk?'
        write(istdou,*) imstsl,imissl,imjssl,idgnrk
        read(istdin,*,end=85,err=90) imstsl,imissl,imjssl,idgnrk
c
      end if
c
      istdpr_adi=istdpr
c
c
c     *********************  End of parameter input *********************
c
c  In initial call with parameter setting, return after input of
c  parameters. Actual calculation is carried out with the subsequent
c  calls.
c
      if(i_paramset.eq.1.and.init_paramset.eq.0) then
        write(istdou,'(//'' End parameter input'')')
        if(istdpr.gt.0.and.istdpr.ne.istdou) 
     *    write(istdpr,'(//'' End parameter input'')')
C
c  store parameters in common /cadr_param/ and /cadi_param/
c
        i_param=1
        call res_adipar(i_param,
     *  el, els1, dels, dfsig1, dfsig2, sig1, sig2, dfsig,
     *  eltrw1, eltrw2, sgtrw1, sgtrw2, 
     *  nsel, nsig1, nsig2, itrsig, nsig, istsig,
     *  inomd1, iscan, irotkr)
        init_paramset=1
c
        ierr_param=1
        return
c
      end if
c
c  ***************************************************************
c
c  Entry point for repeated call with parameter variations.
c
c  ***************************************************************
10000 continue
c
c  begin by opening file for printed output, if different from
c  standard output hardcoded in istdou
c
c  Test for including trailer in output file names
c
      if(i_paramset.eq.0) then
        in_stat ='o'
        out_stat='u'
      else
        in_stat ='ot'
        out_stat='ut'
      end if
c
      istdpr=istdpr_adi
      if(istdpr.ne.istdpr_def) then
        write(istdpr,'(//'' Change istdpr to'',i3)') istdpr
	call stfile(istdpr,nstdpr)
        call openfc(istdpr,istdpp,out_stat,'f')
	write(istdou,'(/'' Opening standard print on unit'',i3,
     *    '' file '',a40)') istdpr, filess(nstdpr)
      end if
c
c  test for resetting parameters from common
c
      if(i_paramset.ne.0.and.init_paramset.ne.0) then
c
	i_param=2
        call res_adipar(i_param,
     *   el, els1, dels, dfsig1, dfsig2, sig1, sig2, dfsig,
     *   eltrw1, eltrw2, sgtrw1, sgtrw2, 
     *   nsel, nsig1, nsig2, itrsig, nsig, istsig,
     *   inomd1, iscan, irotkr)
      end if
c
c  test for resetting parameters due to inconsistencies
c
      if(iscan.gt.1) then
	if(itrsig.gt.1) then
	  if(istdpr.gt.0) write(istdpr,105) iscan, itrsig
	  if(istdpr.ne.istdou) write(istdou,105) iscan, itrsig
	  itrsig = 1
        end if
	if(istsig.gt.1) then
	  if(istdpr.gt.0) write(istdpr,107) iscan, istsig
	  if(istdpr.ne.istdou) write(istdou,107) iscan, istsig
        end if
      end if
c
      if(istsig.gt.1.and.mod(itrsig,2).eq.0.and.nsel.gt.1
     *    .and.iscan.le.0) then
	if(istdpr.gt.0) write(istdpr,110) istsig, itrsig, nsel
	if(istdpr.ne.istdou) write(istdou,110) istsig, itrsig, nsel
	nsel = 0
      end if
c
      if(istsig.gt.1.and.itrsig.eq.0.and.iscan.le.0) then
	if(istdpr.gt.0) write(istdpr,112) istsig
	if(istdpr.ne.istdou) write(istdou,112) istsig
	if(i_paramset.le.0) then
	  go to 6
	else
	  go to 90
        end if
      end if
c
      if(itsord.eq.1.and.iabs(itrsig).lt.2) then
	if(istdpr.gt.0) write(istdpr,114) itrsig
	if(istdpr.ne.istdou) write(istdou,114) itrsig
	itsord = 0
      end if
c
c  test for resetting eps
c
      if(eps.le.epsprc) then
        eps=1.1d0*epsprc
        if(istdpr.gt.0) write(istdou,120) eps
        if(istdou.ne.istdpr) write(istdpr,120) eps
      end if
c
c  test for prohibiting Richardson extrapolation if mdintg = 2, 4 or 5
c
      if((mdintg.eq.2.or.mdintg.eq.4.or.mdintg.eq.5).and.iriche.eq.1)
     *  then
        iriche=0
        write(istdou,934)
        if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,934)
      end if
c                   -----------------------------------
c  store parameters in common /cadr_param/ and /cadi_param/
c
      i_param=1
      call res_adipar(i_param,
     * el, els1, dels, dfsig1, dfsig2, sig1, sig2, dfsig,
     * eltrw1, eltrw2, sgtrw1, sgtrw2, 
     * nsel, nsig1, nsig2, itrsig, nsig, istsig,
     * inomd1, iscan, irotkr)
      init_paramset=1
c
c  unpack nfmode, reset to nfmode mod 10 for consistency
c
      nfmod1=nfmode/10
      nfmode=mod(nfmode,10)
c
c  set flag for output of eigenfunctions during scan
c
      if((iscan.gt.1.and.itrsig.ne.1.and.nfmode.gt.0).or.
     *  (itmax.eq.0.and.nfmod1.eq.1)) then
        nfmscn=1
        write(istdou,'(/'' Output eigenfunctions during scan''/)')
        if(istdpr.gt.0.and.istdpr.ne.istdou) 
     *    write(istdpr,'(/'' Output eigenfunctions during scan''/)')
      else
        nfmscn=0
      end if
c
c  write exec
c  **********
c
      write(istdou,
     *  '(//'' ----------------------------------------------''//)')
      write(istdou,*) ' cntrd'
      write(istdou,*) cntrd,tail
c
      if(ircmod.gt.0) then
c
        write(istdou,*) ' ifind,xmod,imlds,in,nprmod'
        write(istdou,*) ifind,xmod,imlds,in,nprmod,tail
        write(istdou,*) ' xtrnct,ntrnsf,imdmod'
        write(istdou,*) xtrnct,ntrnsf,imdmod,tail
c
      end if
c
      if(ircosc.gt.0) then
c
        write(istdou,*) ' el,nsel,els1,dels,dfsig1,dfsig2,nsig1,nsig2'
        write(istdou,*) 
     *    el,nsel,els1,dels,dfsig1,dfsig2,nsig1,nsig2,tail
        write(istdou,*) ' itrsig,sig1,istsig,inomd1,itrds'
        write(istdou,*) itrsig,sig1,istsig,inomd1,itrds,tail
        write(istdou,*) ' dfsig,nsig,iscan,sig2'
        write(istdou,*) dfsig,nsig,iscan,sig2,tail
        write(istdou,*) ' eltrw1, eltrw2, sgtrw1, sgtrw2'
        write(istdou,*) eltrw1, eltrw2, sgtrw1, sgtrw2,tail
c
      end if
c
      if(ircrot.gt.0) then
        write(istdou,*) ' irotsl, em, nsem, ems1, dems'
        write(istdou,*)   irotsl, em, nsem, ems1, dems, tail
c
      end if
c
      if(irccst.gt.0) then
c
        write(istdou,*) ' cgrav'
        write(istdou,*) cgrav,tail
c
      end if
c
      if(ircint.gt.0) then
c
        write(istdou,*) ' iplneq,iturpr,icow,alb'
        write(istdou,*) iplneq,iturpr,icow,alb,tail
        write(istdou,*) ' istsbc,fctsbc,ibotbc,fcttbc'
        write(istdou,*) istsbc,fctsbc,ibotbc,fcttbc,tail
        write(istdou,*) 
     *    ' mdintg,iriche,xfit,fcnorm,eps,epssol,itmax,dsigre'
        write(istdou,*) 
     *    mdintg,iriche,xfit,fcnorm,eps,epssol,itmax,dsigre,tail
        write(istdou,*) ' fsig,dsigmx,irsevn,xmnevn,nftmax,itsord'
        write(istdou,*) fsig,dsigmx,irsevn,xmnevn,nftmax,itsord,tail
c
      end if
c
      if(ircout.gt.0) then
c
        write(istdou,*) ' istdpr,nout,nprcen,irsord,iekinr'
        write(istdou,*) istdpr,nout,nprcen,irsord,iekinr,tail
        write(istdou,*) ' iper,ivarf,kvarf,npvarf,nfmode'
        write(istdou,*) iper,ivarf,kvarf,npvarf,nfmode,tail
        write(istdou,*) ' irotkr,nprtkr,igm1kr,npgmkr,ispcpr'
        write(istdou,*) irotkr,nprtkr,igm1kr,npgmkr,ispcpr,tail
        write(istdou,*) 'icaswn, sigwn1, sigwn2, frqwn1, frqwn2, ',
     *                'iorwn1, iorwn2, frlwn1, frlwn2'
        write(istdou,*) icaswn, sigwn1, sigwn2, frqwn1, frqwn2, 
     *                iorwn1, iorwn2, frlwn1, frlwn2, tail
c
      end if
c
      if(ircdgn.gt.0) then
c
        write(istdou,*) ' itssol,idgtss,moddet,iprdet,npout'
        write(istdou,*) itssol,idgtss,moddet,iprdet,npout,tail
        write(istdou,*) ' imstsl,imissl,imjssl,idgnrk'
        write(istdou,*) imstsl,imissl,imjssl,idgnrk,tail
c
      end if
c
c  if istdpr.ne.istdou, also print exec
c
      if(istdpr.ne.istdou.and.istdpr.gt.0) then
c
        write(istdpr,*) 'cntrd'
        write(istdpr,*) cntrd,tail
        write(istdpr,*) 'ifind,xmod,imlds,in,nprmod'
        write(istdpr,*) ifind,xmod,imlds,in,nprmod,tail
        write(istdpr,*) 'xtrnct,ntrnsf,imdmod'
        write(istdpr,*) xtrnct,ntrnsf,imdmod,tail
        write(istdpr,*) 'el,nsel,els1,dels,dfsig1,dfsig2,nsig1,nsig2'
        write(istdpr,*) 
     *    el,nsel,els1,dels,dfsig1,dfsig2,nsig1,nsig2,tail
        write(istdpr,*) 'itrsig,sig1,istsig,inomd1,itrds'
        write(istdpr,*) itrsig,sig1,istsig,inomd1,itrds,tail
        write(istdpr,*) 'dfsig,nsig,iscan,sig2'
        write(istdpr,*) dfsig,nsig,iscan,sig2,tail
        write(istdpr,*) 'eltrw1, eltrw2, sgtrw1, sgtrw2'
        write(istdpr,*) eltrw1, eltrw2, sgtrw1, sgtrw2,tail
        write(istdpr,*) ' cgrav'
        write(istdpr,*) cgrav,tail
        write(istdpr,*) 'iplneq,iturpr,icow,alb'
        write(istdpr,*) iplneq,iturpr,icow,alb,tail
        write(istdpr,*) 'istsbc,fctsbc,ibotbc,fcttbc'
        write(istdpr,*) istsbc,fctsbc,ibotbc,fcttbc,tail
        write(istdpr,*) 
     *    'mdintg,iriche,xfit,fcnorm,eps,epssol,itmax,dsigre'
        write(istdpr,*) 
     *    mdintg,iriche,xfit,fcnorm,eps,epssol,itmax,dsigre,tail
        write(istdpr,*) 'fsig,dsigmx,irsevn,xmnevn,nftmax,itsord'
        write(istdpr,*) fsig,dsigmx,irsevn,xmnevn,nftmax,itsord,tail
        write(istdpr,*) 'nout,nprcen,irsord,iekinr'
        write(istdpr,*) nout,nprcen,irsord,iekinr,tail
        write(istdpr,*) 'iper,ivarf,kvarf,npvarf,nfmode'
        write(istdpr,*) iper,ivarf,kvarf,npvarf,nfmode,tail
        write(istdpr,*) 'irotkr,nprtkr,igm1kr,npgmkr,ispcpr'
        write(istdpr,*) irotkr,nprtkr,igm1kr,npgmkr,ispcpr,tail
        write(istdpr,*) 'icaswn, sigwn1, sigwn2, frqwn1, frqwn2, ',
     *                'iorwn1, iorwn2, frlwn1, frlwn2'
        write(istdpr,*) icaswn, sigwn1, sigwn2, frqwn1, frqwn2, 
     *                iorwn1, iorwn2, frlwn1, frlwn2, tail
        write(istdpr,*) 'itssol,idgtss,moddet,iprdet,npout'
        write(istdpr,*) itssol,idgtss,moddet,iprdet,npout,tail
        write(istdpr,*) 'imstsl,imissl,imjssl,idgnrk'
        write(istdpr,*) imstsl,imissl,imjssl,idgnrk,tail
c
        write(istdpr,125)
      end if
c
c  set itrsga to test for type of trial from file, regardless
c  of sign (and hence precision)
c
      itrsga=iabs(itrsig)
c
c  store original fsig
c
      fsig0=fsig
c
c  open variable or conditional files (note: model file is opened 
c  in sr/ readml)
c
      if(itrsga.ge.2) then
	itr=itrsga/2
	call openfc(itrds,itrdsp,'o',ctrdsf(itr))
      end if
c
c  open status log. Use standard name if not set in ofiles
c
      call stfile(idslog,ndslog)
      if(ndslog.eq.-1) then
	write(istdou,102)
        open(idslog,file='adipls-status.log',status='unknown')
      else
	call openf(idslog,'u','f')
      end if
      write(idslog,*) "  "
c
      call openfc(idsgsm,idsgsp,out_stat,'u')
      call openfc(idsssm,idsssp,out_stat,'u')
      if(nfmode.ge.1) call openfc(idsefn,idsefp,out_stat,'u')
      if(irotkr.ge.1.and.idsrkr.gt.0) 
     *  call openfc(idsrkr,idsrkp,out_stat,'u')
      if(igm1kr.eq.1) call openfc(idsgkr,idsgkp,out_stat,'u')
c
c  reset kvarfc if kvarf is changed
c
      if(kvarf.ne.kvarfp) then
        kvarfc=kvarf
        kvarfp=kvarf
      end if
c
      nwmod=0
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'set nwmod no. 1 to',nwmod
c 
c  for i_inout = 0, assume that model quantities are provided as 
c  arguments
c
      if(i_inout.eq.0) then
c
c  store variables in x, aa, possibly resetting number of mesh points 
c  or iturpr, depending on iriche or data_arg
c
	call res_adimod(x_arg, aa_arg, data_arg, nn_arg, ivarmd, 
     *    in, data, x, aa, iaa_arg, iaa, nn0, icry)
c
	nwmod=1
c
c  Otherwise, in first run read model
c
      else if(i_inout.gt.0.and.(frsrun.or.ifind.ge.0)) then
c
        xmod1=xmod
c
        call readml(ifind,xmod1,mdintg,imlds,imldsp,in,data,x,aa,iaa,
     *    nn0,nprmod,ivarmd,mlname,nwmod,icryml)
c
c  test for success
c
        if(icryml.eq.-1) then
	  go to 6
        else if(icryml.eq.2) then
          write(istdou,910) xmod
          if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,910) xmod
        end if
c
c  possibly set rotation rate and read in output file for obs file
c
	if(irotkr.eq.11.or.irotkr.eq.21) then
	  call set_rotation(x, nn0, icontr)
	  if(init_obsrot.eq.0) read(istdin,'(a)') obsrot_file
	  init_obsrot=1
        end if
      end if
c
c  test whether model was reset
c
      if(nwmod.ne.0) then
c
c  reset frsrun
c
        frsrun=.false.
c
        nwmod=1
        nfmesh=1
        an=data(7)
c
c  set total number of points (note: this may be reset by s/r trnmod
c  below)
c
	nn = nn0
c
c  factor for period calculation (period in minutes)
c
        perfac=(pi/30)*sqrt(data(2)**3/(data(1)*cgrav))
c
c  reset kvarfc for new model:
c
        kvarfc=kvarf
c
c  reset counter for storage of mode quantities in common/cobs_param/
c
	nobs_st=0  
c
      end if
c
c  end of model input section
c  **************************
c
c                         -----------------------------------
c
c  test for modifications to the model
c
      if(nwmod.eq.1.or.ifind.eq.-2) then
c
        if(imdmod.ne.0) then
	  call modmod(x,aa,data,nn,ivarmd,iaa,imdmod)
	  nn0 = nn
        end if
c
c  test for truncation of model
c
        if(xtrnct.gt.0.or.ntrnsf.ge.1) then
c
          call trnmod(x,nn0,xtrnct,ntrnsf,ntrnct,xtrnc1,ntrns1,
     *      nw1,nibc,nn,fulmod)
c
c  reset flags for new model
c
          if(ntrncp.ge.0.and.(ntrncp.ne.ntrnct.or.nnp.ne.nn)) then
            nwmod=1
	    nfmesh=1
	  end if
c
	else
c
c  set parameters for normal model, depending on whether or not
c  this is an envelope model
c
	  nw1 = 1
	  if(x(1).eq.0) then
	    nibc = 2
	    fulmod = .true.
	    xtrnc1 = 0
          else
	    nibc = 1
	    fulmod = .false.
	    xtrnc1 = x(1)
	  end if
	  nn = nn0
	  ntrnct = 0
        end if
c
	nnp = nn
	ntrncp = ntrnct
      end if
c
c  set surface values into datmd1
c
      datmd1(1)=aa(2,nn)
      datmd1(2)=aa(5,nn)
c
c  output steps
c
      nnw=nn-nw1+1
c
c  fitting point
c  *************
c
      if(xfit.eq.-1) then
c
c  for xfit = -1 temporarily set nfit to nn
c
        nfit=nn
c
      else
        nfit=-1
        do 15 n=1,nn
        if(x(n).lt.xfit) nfit=n
   15   continue
      end if
c
c  reset nfit because of truncation when mdintg .ne. 3
c
      if(nfit.le.nibc.and.mdintg.ne.3) then
c
c  note: in the following two lines, frfit is used purely
c  as a temporary variable
c
        frfit=nfit/float(nn)
        nfit=nw1-1+frfit*(nn+1-nw1)
        write(istdou,920) ntrnct,nfit
        if(istdou.ne.istdpr.and.istdpr.gt.0) 
     *    write(istdpr,920) ntrnct,nfit
      end if
c
      xfit1=x(nfit)
c
c  set limits in frequency scan as squared dimensionless frequency
c
      if(iscan.gt.1) then
	if(istsig.le.1) then
	  sigs1=sig1
	  sigs2=sig2
        else
	  sigfct=4.d-6*pi**2*data(2)**3/(data(1)*cgrav)
	  sigs1=sigfct*sig1**2
	  sigs2=sigfct*sig2**2
        end if
      end if
c
c  test for step in el
c  *******************
c
      nel=nsel
c
      if(nsel.ge.1) then
        el=els1
        nel=1
      end if
c
c  entry point for continuing step in el
c
   20 continue
c
c  carry out various tests depending on the degree
c
c  test for radial oscillation
c
      radial=el.le.1.e-6
c
c  for radial oscillations, exclude xfit = -1
c
      idg930=0
      if(radial.and.xfit.eq.-1) then
        idg930=1
        nfit=0.5*nn
      end if
c
c  test for resetting nfit with Richardson extrapolation
c
      if(mod(nfit-nibc,2).ne.0.and.iriche.eq.1) then
        nfit=nfit+1
      end if
c
      xfit1=x(nfit)
c
c  test whether cowling approximation must be used, and if so reset
c  icow if necessary
c
c  cowling approximation is forced for:
c  1) iplneq = 1
c  2) truncated model for nonradial oscillations
c
      idg924=0
      nonrad=.not.radial.or.nsel.gt.1
      if((iplneq.eq.1
     *  .or.(nonrad.and..not.fulmod))
     *  .and.icow.ne.2.and.icow.ne.3) then
        icow=2
        idg924=1
      end if
c
c  test for proper use of mdintg
c
      idg926=0
      idg928=0
      if(mdintg.eq.2.and..not.radial.and.icow.ne.2) then
        mdintg=4
        idg926=1
      else if(mdintg.eq.4.and.(radial.or.icow.eq.2)) then
        mdintg=2
        idg928=1
      end if
c
c  skip radial oscillations in plane-parallel case
c
      if(radial.and.iplneq.eq.1) then
        write(istdou,914)
        if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,914)
        go to 80
      end if
c
c  test for step in em
c  *******************
c
      if(irotsl.eq.1.and.nsem.ne.0) then
        if(nsem.eq.-1) then
	  em=-el
        else
	  em=ems1
        end if
	nem=1
      end if
c
c  entry point for continuing step in em
c
   25 continue
c
c  test for consistent el and em
c
      if(irotsl.eq.1.and.em.gt.el) then
        write(istdou,127) em, el
        if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,127) em, el
	go to 80
      end if
c
c  new sigma**2
c  ************
c
      sigp=sig
      if(iscan.gt.1) then
        sig=sigs1
      else
        sig=sig1
      end if
      jstsig=istsig
c
c  test for step in sig
c
      if(istsig.gt.1.and.iscan.le.1) then
	jstsig=1
      else
	jstsig=istsig
      end if
c
c  set sig
c  Entry point for continuing step in sig
c  Test for using step in sig (note: for first step with itrsig = 1
c  use frequency given in sig1)
c
   30 if(iscan.le.1) then
	if(istsig.ge.1.and.
     *    (itrsig.ge.2.or.itrsig.eq.1.and.jstsig.ge.2)) then
c
          inomde=inomd1+jstsig-1
	  sigp=sig
          sig=signew(itrsig,nsig,1,inomde,sigp,dfsig,itrds,el,icry)
          if(icry.eq.-1) then
	    go to 75
          else if(icry.eq.-2) then
	    go to 6
          end if
        end if
      end if
c
c  Ready to start for new mode. Write el
c
      if(irotsl.ne.1) then
        write(istdou,180) jstsig, el, sig
      else
        write(istdou,182) jstsig, el, em, sig
      end if
      if(istdpr.gt.0.and.istdpr.ne.istdou) then
        write(istdpr,125)
        if(irotsl.ne.1) then
          write(istdpr,130) el
        else
          write(istdpr,131) el, em
        end if
      end if
c
c  zero part of /csumma/ containing results from previous mode
c
      call zero(csummm(18),21)
c
c  test for diagnostic output
c
      if(idg924.eq.1) write(istdou,924)
      if(idg926.eq.1) write(istdou,926)
      if(idg928.eq.1) write(istdou,928)
      if(idg930.eq.1) write(istdou,930)
c
      if(istdpr.gt.0.and.istdpr.ne.istdou) then
        if(idg924.eq.1) write(istdpr,924)
        if(idg926.eq.1) write(istdpr,926)
        if(idg928.eq.1) write(istdpr,928)
        if(idg930.eq.1) write(istdpr,930)
      end if
c
      ell=el*(el+1)
      if(iplneq.eq.1) then
        ell=el*el
        write(istdou,132)
        if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,132)
      end if
c
c  This needs a little attention later
c
      istsb1=istsbc
c
      if(istdpr.gt.0) then
c
        if(alb.ne.1) write(istdpr,135) alb
c
c  inclusion of beta+ solution
c
        if(fctsbc.ne.0) write(istdpr,138) fctsbc
c
        if(istsbc.eq.1) write(istdpr,136)
c
        write(istdpr,137)
      end if
c
      if(.not.fulmod) then
        write(istdou,139) xtrnc1
        if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,139) xtrnc1
      end if
c
      if(ntrnsf.eq.nn) then
        if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,140) datmd1
        write(istdou,140) datmd1
      end if
c
c  in case of Richardson extrapolation, decrease trial frequency
c
      if(iriche.eq.1) then
        sig=sig-dsigre
      end if
c
c  test for scan 
c
      if(iscan.le.1) then
c
c  normal iteration
c
        call sigsol(sig,x,y,iy,nw1,nibc,nn,nnw,mdintg,nev1,nfit,istsb1,
     *    iscan,iord,icry,isolcv,isigcv,sigtst)
c
        if(icry.lt.0) then
	  if(i_paramset.le.0) then
	    go to 6
          else
	    go to 90
          end if
        end if
c
c  call output routine
c
        call sigout(sig,x,y,iy,nw1,nibc,nn,nnw,mdintg,nev1,nfit,
     *    istsb1,iord,isolcv,isigcv)
c
      else
c
        call sigscn(sigs1,sigs2,iscan,nsig,itrsig,x,y,iy,
     *    nw1,nibc,nn,nnw,mdintg,nev1,nfit,istsb1,icry)
c
        if(icry.lt.0) then
	  if(i_paramset.le.0) then
	    go to 6
          else
	    go to 90
          end if
        end if
c
      end if
c
c  test for continued step in sig, for itrsig .ge. 1.
c
   75 if(jstsig.lt.istsig.and.itrsga.ge.1) then
	jstsig=jstsig+1
	go to 30
      end if
c
c  continue step in em?
c
      if(irotsl.eq.1.and.nsem.ne.0) then
        if((nsem.eq.-1.or.nem.lt.nsem).and.em.le.el-dems+1.d-10) then
	  nem=nem+1
	  em=em+dems
c
c  test for resetting of sigs1, sigs2 
c
          if(iscan.gt.1) then
            if(dfsig1.ne.0) then
              sigs1p=sigs1
              sigs1=signew(1,nsig1,1,0,sigs1p,dfsig1,10,el,icry)
            end if
            if(dfsig2.ne.0) then
              sigs2p=sigs2
              sigs2=signew(1,nsig2,1,0,sigs2p,dfsig2,10,el,icry)
            end if
          end if
	  go to 25
        end if
      end if
c
c  continue step in el?
c
   80 if(nel.lt.nsel) then
        nel=nel+1
        sig=sig+dfsig
        el=el+dels
c
c  test for resetting of sigs1, sigs2 
c
	if(iscan.gt.1) then
          if(dfsig1.ne.0) then
	    sigs1p=sigs1
	    sigs1=signew(1,nsig1,1,0,sigs1p,dfsig1,10,el,icry)
          end if
          if(dfsig2.ne.0) then
	    sigs2p=sigs2
	    sigs2=signew(1,nsig2,1,0,sigs2p,dfsig2,10,el,icry)
          end if
        end if
c
	go to 20
      end if
c
c  read new exec
c  *************
c
      if(i_paramset.le.0) go to 6
c
c                                   ---------------------------------
c  Normal end
c
   85 ierr_param=0
c
c  close output files
c
   90 close(idsgsm)
      close(idsssm)
c
      if(idsefp.gt.0) close(idsefn)
      if(idsrkp.gt.0) close(idsrkr)
      if(idsgkp.gt.0) close(idsgkr)
c
c  test for output of rotation results
c
      if(init_obsrot.eq.1) call dump_obs(strcompr(obsrot_file))
      if(istdpr.gt.0.and.istdou.ne.istdpr) write(istdpr,'(//
     *  '' Exiting adipls''/
     *  '' --------------------------------------------------''/)')
      write(istdou,'(//
     *  '' Exiting adipls''/
     *  '' --------------------------------------------------'')')
c
      if(istdpr.ne.istdpr_def) then
        close(istdpr)
      end if
      return
  102 format(//' Log in standard name: adipls-status.log')
  103 format(//' ***** Warning. Empty cntrd: ',a35/
     *         '       Execution terminated.')
  105 format(/' **** Warning. For iscan =',i5,' itrsig was reset from',
     *  i5,' to 1')
  107 format(/' For iscan =',i5,
     *  ' setting frequency limits for scan with istsig =', i5)
  110 format(/' **** Warning. For istsig =',i5,' itrsig = ',i2,
     *  ' nsel reset from ',i5,' to 0')
  112 format(/
     *   ' **** Error. istsig =',i5,' not allowed with itrsig = 0'/
     *   '      Case skipped.')
  114 format(/
     *  ' **** Warning. For itrsig = ',i2,' itsord = 1 not allowed'/
     *  '      itsord reset to 0')
  120 format(//' **** Warning. eps has been reset to',1pe13.5)
  125 format(///10x,60(1h*))
  127 format(//' ***** Error in s/r adipls. em =',1pe13.5,' gt el =',
     *   e13.5)
  130 format(//'     l   =',f12.7)
  131 format(//'     l   =',f12.7,'  m =',f12.7)
  132 format(//' plane-parallel equations')
  135 format(//'  lambda=',f10.7)
  136 format(/' surface mechanical boundary condition assumes',
     *  ' isothermal atmosphere')
  137 format(/)
  138 format(///' at outer boundary beta+/beta- =',f10.5)
  139 format(//' model truncated at x =',f10.7)
  140 format(//' at surface model stops at v/gamma1 =',1pe13.5,
     *  ' and u =',e13.5)
  180 format(/i5,'. Start iteration for l =',f10.3,' trial sigma**2 =',
     *  1pe12.4)
  182 format(/i5,'. Start iteration for l =',f10.3,' m =',f10.3,
     *  ' trial sigma**2 =', 1pe12.4)
  910 format(/1x,5(1h*),' xmod =',f10.5,' is non-integral.',
     *  ' integer part used.')
  914 format(//' ***** error. programme cannot solve radial',
     *  ' equations in plane-parallel case')
  920 format(//'  model truncated at n =',i5,'  nfit reset to',i5)
  924 format(/' ***** warning. cowling approximation is forced.',
     *  ' icow has been reset to 2')
  926 format(/' ***** warning. mdintg = 2 for full nonradial case'/
     *  ' mdintg has been reset to 4')
  928 format(/
     *  ' ***** warning. mdintg = 4 for radial case or Cowling appr.'/
     *  ' mdintg has been reset to 2')
  930 format(/' ***** warning. xfit = -1 not allowed for radial',
     *  ' oscillations. xfit reset to x(0.5*nn)')
  934 format(//' ***** error. Richardson extrapolation not allowed',
     *  ' for mdintg = 2 or 4'/
     *         '              iriche reset to 0')
      end
      subroutine setlog(idslog,istatus,el,iord,iordtr,sig,
     *  eps,epssol,ddsig,ddsol,iord1,iord2,sigrc1,sigrc2,dsolrc,
     *  dsig1,dsig2,xfit1,xfitst)
c
c  outputs diagnistics to log file in unit 99, if iteration
c  has failed to converge. 
c  istatus is determined by the convergence problem:
c    istatus = 1: frequency iteration failed to converge
c    istatus = 2: eigenfunction iteration failed to converge
c    istatus = 3: change of xfit to obtain convergence
c    istatus = 4: different orders and eigenfunctions
c                 in Richardson extrapolation
c    istatus = 5: different orders but similar eigenfunctions
c                 in Richardson extrapolation. Solution accepted.
c    istatus = 6: orders in Richardson extrapolation 
c                 brought to agree through change of xfit.
c    istatus = 7: order of solution differs from order of trial mode
c    istatus = 8: order of solution brought to agree with
c                 order of trial mode through change of xfit.
c    istatus = 9: different gravitational-potential corrections
c                 to frequency with Richardson extrapolation
c
c  Original version: 21/7/94
c
      implicit double precision (a-h, o-z)
      character*280 file, filess
      parameter(iaa=10)
      common/cofile/ nfiles, ids(99), file(99), iopen(99), filess(99)
      common/rhsdat/ elcom,ell,alb,els,el1,sigcom,anres,perfac,data(8),
     *  aacom(iaa,1)
      common/rotdat/ em, irotsl
      common/cprcns/ epsprc,epsufl,epsofl,eprufl
c
c  common defining file input and output
c
      common/cdadsg/ imlds,idslg1,idsgsm,idsssm,idsefn,
     *  idsrkr,idsgkr,itrds
c
      data init/1/
c..      write(6,*) 'Enter setlog with istatus,el,iord,per,sig,',
c..     *  ' eps,epssol,ddsig,ddsol,iord1,iord2,sigrc1,sigrc2,dsolrc,',
c..     *  ' xfit1,xfitst ='
c..      write(6,*) istatus,el,iord,per,sig,
c..     *  eps,epssol,ddsig,ddsol,iord1,iord2,sigrc1,sigrc2,dsolrc,
c..     *  xfit1,xfitst
c
c  In initial call, set model and grand summary file names
c
      if(init.eq.1) then
	init=0
	call stfile(imlds,nmlds)
	call stfile(idsgsm,ndsgsm)
	write(idslog,100) file(nmlds),file(ndsgsm)
      end if
c
      per=perfac/sqrt(max(abs(sig),epsufl))
      anu=16666.667d0/per
      write(idslog,105) el, iord, anu, sig
      if(istatus.eq.1) then
	write(idslog,110) ddsig, eps
      else if(istatus.eq.2) then
	write(idslog,120) ddsol, epssol
      else if(istatus.eq.3) then
	write(idslog,130) xfit1, xfitst
      else if(istatus.eq.4) then
	write(idslog,140) iord1, iord2, dsolrc
      else if(istatus.eq.5) then
	write(idslog,150) iord1, iord2, dsolrc
      else if(istatus.eq.6) then
	write(idslog,160) xfit1, xfitst
      else if(istatus.eq.7) then
	write(idslog,170) iord, iordtr
      else if(istatus.eq.8) then
	write(idslog,180) xfit1, xfitst
      else if(istatus.eq.9) then
	write(idslog,190) dsig1, dsig2
      end if
      return
  100 format(' Model file:   ', a60/
     *       ' Grand summary:', a60)
  105 format(/' ----------------------------------------------------'/
     * ' Mode l =',f6.1,' order =',i5,' nu =',
     *  1pe12.3,' microHz, sigma**2 =',e13.5)
  110 format(' Error. Frequency iteration unconverged.'/
     *  ' Last correction =',1pe13.5,' exceeds limit =', e13.5)
  120 format(' Error. Eigenfunction iteration unconverged.'/
     *  ' Last normalized misfit =',1pe13.5,' exceeds limit =', e13.5)
  130 format(' Warning: xfit changed to ',f10.5,
     *  ' from standard value ',f10.5/
     *       '          to ensure convergence')
  140 format(' Error in Richardson extrapolation. Orders =',
     *  2i5,' differ'/
     *  ' Relative eigenfunction difference =',1pe13.5,
     *  ' is excessive')
  150 format(' Warning in Richardson extrapolation. Orders =',
     *  2i5,' differ'/
     *  ' Relative eigenfunction difference =',1pe13.5,
     *  ' is acceptable')
  160 format(' Warning: xfit changed to ',f10.5,
     *  ' from standard value ',f10.5/
     *       '          to obtain correct order in Richardson ',
     *  'extrapolation')
  170 format(' Error: computed order',i5,' disagrees with trial order',
     *  i5)
  180 format(' Warning: xfit changed to ',f10.5,
     *  ' from standard value ',f10.5/
     *       '          to obtain agreement in order')
  190 format(' Error: With Richardson extrapolation,'/
     *       ' gravitational potential corrections =',1p2e11.3,
     *       ' differ excessively')
      end
      subroutine testri(x,y,yri,iy,iyri,nn,iasn,aa,ia,
     *      el,iord1,iord2,sigrc1,sigrc2,dsgrcp,dsolrc,isigcv,istatus,
     *      idgtri)
c
c  Carry out various tests to see whether Richardson extrapolation must
c  rejected when the two mode orders differ.
c
c  Returns istatus = 4, isigcv = -2 for rejection.
c  If mode is accepted, returns isigcv = 1 but istatus = 4 for warning.
c
c  If idgtri .gt. 1, extensive diagnostics is printed
c
c  Original version: 21/7/94
c
c  Modified 26/4/99, changing dsol to dsolrc 
c
      implicit double precision (a-h, o-z)
      parameter(iw=10)
      dimension x(1),y(iy,1),yri(iyri,1),aa(ia,1)
      common/xarr1/ x1(1)
      common/worksp/ wrk(iw,1)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  set normalized norm of difference between eigenfunctions
c
      if(x(1).eq.0) then
        nc=2
      else
        nc=1
      end if
c
      n1=0
      if(el.gt.1.e-6) then
        elli=1./(el*(el+1))
        ellis=sqrt(elli)
        do 20 n=nc,nn,iasn
        n1=n1+1
        wrk(3,n1)=x(n)*sqrt(aa(1,n)*aa(5,n))
	wrk(4,n1)=wrk(3,n1)*(y(1,n)-yri(1,n))
	wrk(5,n1)=ellis*wrk(3,n1)*(y(2,n)-yri(2,n))
	wrk(6,n1)=wrk(3,n1)*y(1,n)
	wrk(7,n1)=ellis*wrk(3,n1)*y(2,n)
	x1(n1)=x(n)
	wrk(1,n1)=wrk(6,n1)*wrk(6,n1)+wrk(7,n1)*wrk(7,n1)
	wrk(2,n1)=wrk(4,n1)*wrk(4,n1)+wrk(5,n1)*wrk(5,n1)
        if(idgtri.gt.1.and.istdpr.gt.0) 
     *    write(istdpr,'(i5,f10.5,1p2e15.7)') 
     *    n, x(n), wrk(6,n1), wrk(3,n1)*yri(1,n)
   20   continue
      else
        do 25 n=nc,nn,iasn
        n1=n1+1
        wrk(3,n1)=x(n)*sqrt(aa(1,n)*aa(5,n))
	wrk(4,n1)=wrk(3,n1)*(y(1,n)-yri(1,n))
	wrk(6,n1)=wrk(3,n1)*y(1,n)
	x1(n1)=x(n)
	wrk(1,n1)=wrk(6,n1)*wrk(6,n1)
	wrk(2,n1)=wrk(4,n1)*wrk(4,n1)
        if(idgtri.gt.1.and.istdpr.gt.0) 
     *    write(istdpr,'(i5,f10.5,1p2e15.7)') 
     *    n, x(n), wrk(6,n1), wrk(3,n1)*yri(1,n)
   25   continue
      end if
c
      nn1=n1
c
c  integrate norm and difference
c
      call vinta(x1,wrk(1,1),wrk(3,1),nn1,iw,iw)
      call vinta(x1,wrk(2,1),wrk(4,1),nn1,iw,iw)
      if(idgtri.gt.1.and.istdpr.gt.0) 
     *  write(istdpr,'(i5,0pf10.5,1p4e13.5)') 
     *  (n1, x1(n1), (wrk(i,n1), i=1,4),n1=1,nn1)
c
      dsolrc=wrk(4,nn1)/wrk(3,nn1)
c
c  carry out test. Note: these may have to be fine-tuned
c
      if(idgtri.gt.1.and.istdpr.gt.0) 
     *  write(istdpr,*) 'sigrc1, sigrc2', sigrc1, sigrc2
      dsig=abs(sigrc1-sigrc2)
      if(idgtri.gt.1.and.istdpr.gt.0) 
     *  write(istdpr,*) 'dsig, dsolrc, dsgrcp',dsig, dsolrc, 
     *  dsgrcp
      dsigr=dsig/sigrc1
      if(dsgrcp.gt.0) then
	dsgrc1=dsgrcp
      else
	dsgrc1=sigrc1
      end if
c
      if(dsig.gt.10*dsgrcp.or.dsolrc.gt.5.*dsigr.or.dsolrc.gt.0.1) 
     *  then
	istatus=4
	isigcv=-2
        write(istdou,110) iord1, iord2, dsolrc, dsigr
        if(istdpr.gt.0.and.istdpr.ne.istdou)
     *    write(istdpr,110) iord1, iord2, dsolrc, dsigr
      else
	istatus=5
	isigcv=1
        write(istdou,120) iord1, iord2, dsolrc, dsigr
        if(istdou.ne.istdpr.and.istdpr.gt.0) 
     *    write(istdpr,120) iord1, iord2, dsolrc, dsigr
      end if
      return
  110 format(//' *******  Error in Richardson extrapolation.',
     *                  ' Orders =',2i5,' differ'/
     *         '          rel. change in solution =',1pe13.5/
     *         '          rel. change in sigma =',e13.5/
     *         '          Solution rejected')
  120 format(//' *******  Warning in Richardson extrapolation',
     *                  ' orders =',2i5,' differ'/
     *         '          rel. change in solution =',1pe13.5/
     *         '          rel. change in sigma =',e13.5/
     *         '          Solution accepted')
      end
      subroutine dmpfls(i)
c
c  dump contents of commons /cofile/ and /cdadsg/
c
      character*280 file,filess
      common/cofile/ nfiles, ids(99), file(99), iopen(99), filess(99)
c
c  common defining file input and output
c
      common/cdadsg/ imlds,idslog,idsgsm,idsssm,idsefn,
     *  idsrkr,idsgkr,itrds
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      if(istdpr.gt.0) then
        write(istdpr,*) ' Call of dmpfls with i =',i
        write(istdpr,*) ' nfiles =',nfiles
        if(i.gt.0) then
          write(istdpr,*) ' ids =',ids
          write(istdpr,'(a60)') 'file:',file
        end if
      end if
      return
      end
      subroutine dump_csumma(indicator)
      implicit double precision (a-h, o-z)
      dimension csummm(50)
      character*(*) indicator
c
      common/csumma/ xmod1,datsum(8),datmd1(2),xtrnc1,dum13,xfit1,
     *  fsbcsm,fcbcsm,albsum,
     *  elsum,ordsum,sigsum,sigc,dmx,xmx,ekin,per,perv,frqv,
     *  ddsig,ddsol,ysum(4),dosdum(5),
     *  in,nnw,mdints,ivarfs,icase,iorign,idum7,idum8,mlname(4)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      equivalence(xmod1,csummm(1))
      write(istder,*) indicator,' in',in
      return
      end
