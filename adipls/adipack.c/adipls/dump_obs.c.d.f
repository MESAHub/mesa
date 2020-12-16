      subroutine dump_obs(obs_file)
c
c  dumps mode quantities in common/cobs_param/  to file outfile
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
c  The actual quantity returned is flagged by icobs_st, set when 
c  nobs_st = 0 (and assumed valid for the full set for the given
c  model).
c
c  original version 8/7/05
c
c  Modified 3/8/05, allowing for rotational results
c
      implicit double precision (a-h,o-z)
      character obs_file*(*)
c
c  common for storage of modal parameters 
c  degree, order, cyclic frequency (microHz), inertia
c
      common/cobs_param/ icobs_st, nobs_st, obs_st(10,1)
c
c  common defining standard input and output, standard error
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  test for no data
c
      if(nobs_st.le.0) then
        write(istdou,105)
        if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,105)
        return
      end if
c
      open(99,status='unknown',file=obs_file,form='formatted')
      write(99,110)
      icobs_st0=mod(icobs_st,10)
      icobs_st1=icobs_st/10
      if(icobs_st0 .eq. 1) then
        write(99,120) 'Variational frequency'
      else if(icobs_st0 .eq. 4) then
        write(99,120) 'Raw eigenfrequency'
      else if(icobs_st0 .eq. 5) then
        write(99,120) 'Richardson extrapolation frequency'
      else 
        write(99,120) 'Eigenfrequency (possibly corrected)'
      end if
c
      if(icobs_st1.eq.0) then
	write(99,125)
	imax=4
      else if(icobs_st1.gt.1) then
         write(99,129)
         imax=8
      else
	write(99,127)
	imax=6
      end if
c
      do n=1,nobs_st
        ll=intgpt(obs_st(1,n)+0.5)
        nord=intgpt(obs_st(2,n)+0.5)
        write(99,130) ll, nord, (obs_st(i,n),i=3,imax)
      end do
      close(99)
c
      return
  105 format(//' ***** Error in s/r dump_obs. No data stored'/)
  110 format('#  Dump contents of common/cobs_param/'/'#')
  120 format('#  Frequency from: ',a/'#')
  125 format('# degree, order, frequency (microHz), inertia'/'#')
  127 format('# degree, order, frequency (microHz), inertia, ',
     *   ' beta, splitting (microHz)'/'#')
  129 format('# degree, order, frequency (microHz), inertia, ',
     *   ' beta, splitting1 (microHz), splitting2 (microHz), ',
     *   ' splitting3 (microHz) '/'#')
  130 format(2i5,f10.3,1pe13.5,0pf10.6,1pe13.5,1pe13.5,1pe13.5)
      end
