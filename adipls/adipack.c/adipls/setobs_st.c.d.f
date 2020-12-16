      subroutine setobs_st(cs,ics,iriche)
c
c  sets mode quantities in common/cobs_param/ from grand summary
c  and possibly also rotational results. Output is controlled by
c  parameters in common/coutpt/.
c
c  The frequency is set, depending on ispcpr:
c
c  ispcpr = 1: variational frequency.
c  ispcpr = 4: from eigenfrequency in cs(20).
c              Note that this allows setting Cowling
c              approximation frequency.
c  ispcpr = 5: from Richardson extrapolation frequency
c              in cs(37), if this is set.
c              Otherwise variational frequency is used.
c  ispcpr = 6: from (possibly corrected) eigenfrequency
c              in cs(21).
c
c  If ispcpr = 0, the frequency is set in order of priority, depending
c  on the availability of the relevant quantity: 5, 1, 6
c  The same order of priority is used if the quantitity selected by
c  another value of ispcpr is not available.
c
c  If irotkr gt 0 in addition sets rotational-splitting quantities.
c  Only spherical rotation is considered.
c  So far only implemented for irotkr = 11, for first-order splitting.
c  Here obs_st(5,.) = beta and obs_st(6,.) = split(1) (the uniform splitting
c  between adjacent m values).
c
c  For irotkr = 21 second-order splitting is computed in addition to the 
c  first-order terms above. Here obs_st(7,.) = split(2) (the m^0 coefficient)
c  and obs_st(8,.) = split(3) (the m^2 coefficient).
c
c  The actual quantities returned are flagged by 
c  icobs_st = icobs_st0 + 10*icobs_st1. This is set when 
c  nobs_st = 0 (and assumed valid for the full set for the given model).
c
c  icobs_st0 flags for the choice of cyclic frequency and icobs_st1 = 1
c  flags for the presence of rotational splitting quantities.
c  (Later to be expanded to allow for higher-order rotation also.)
c  icobs_st1 = 11 flags for the presence of second-order rotational 
c  splitting quantities.
c
c  original version 8/7/05
c
c  Modified 3/8/05 to include rotational results.
c
c  Modified 18/5/07 to include second-order rotational results - KDB
c
c
      implicit double precision (a-h,o-z)
      dimension cs(*),ics(*)
      common/rhsdat/ el,ell,alb,els,el1,sig,anres,perfac
      common/crot_split/ beta, split(4)
c
c  common controlling output
c
      common/coutpt/ nout,nprcen,iper,irotkr,nprtkr,igm1kr,npgmkr,
     *  nfmode,nfmesh,ispcpr,npout,nobs_stmx,nfmscn
c
c  common for storage of modal parameters 
c  degree, order, cyclic frequency (microHz), inertia
c
      common/cobs_param/ icobs_st, nobs_st, obs_st(10,1)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      data idiag_er /0/
c
      save 
c
c  test for selecting the relevant quantity
c
      if(nobs_st.eq.0) then
        if(ispcpr.eq.0.or.ispcpr.eq.5) then
	  if(iriche.eq.1) then
	    icobs_st0=5
	  else if(iper.eq.1) then
	    icobs_st0=1
	  else
	    icobs_st0=6
	  end if
        else if(ispcpr.eq.1) then
	  if(iper.eq.1) then
	    icobs_st0=1
	  else
	    icobs_st0=6
	  end if
        else 
	  icobs_st0=ispcpr
        end if
c
c  test for rotational results
c
	if(irotkr.eq.11) then
	  icobs_st1=1
        else if(irotkr.eq.21) then
          icobs_st1=11
        else
	  icobs_st1=0
        end if
	icobs_st=icobs_st0+10*icobs_st1
      end if
c
c  test for excessive number of modes
c
      if(nobs_st.ge.nobs_stmx) then
        if(idiag_er.eq.0) write(istder,'(//
     *   '' ***** Error in setobs_st.'',
     *   '' Number of modes exceeds maximum = '',i7)') nobs_stmx
        idiag_er=1
        write(istdou,'(//
     *   '' ***** Error in setobs_st.'',
     *   '' Number of modes exceeds maximum = '',i7)') nobs_stmx
        if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,'(//
     *   '' ***** Error in setobs_st.'',
     *   '' Number of modes exceeds maximum = '',i7)') nobs_stmx
        return
      end if
c
      nobs_st=nobs_st+1
      obs_st(1,nobs_st)=cs(18)
      obs_st(2,nobs_st)=cs(19)
      obs_st(4,nobs_st)=cs(24)
      if(icobs_st0.eq.5) then
        obs_st(3,nobs_st)=1000d0*cs(37)
      else if(icobs_st0.eq.1) then
        obs_st(3,nobs_st)=1000d0*cs(27)
      else if(icobs_st0.eq.6) then
        obs_st(3,nobs_st)=sqrt(cs(21))/(60.d-6*perfac)
      else 
        obs_st(3,nobs_st)=sqrt(cs(20))/(60.d-6*perfac)
      end if
c
c  test for rotational results
c
      if(irotkr.eq.11) then
	obs_st(5,nobs_st)=beta
	obs_st(6,nobs_st)=split(1)
      end if
      if(irotkr.eq.21) then
	obs_st(5,nobs_st)=beta
	obs_st(6,nobs_st)=split(1)
        obs_st(7,nobs_st)=split(2)
        obs_st(8,nobs_st)=split(3)
      end if      
      return
      end
