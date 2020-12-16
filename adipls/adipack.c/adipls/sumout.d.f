      subroutine sumout(iw,ih,nwmod,an,xtrcnt,icase,data,
     *  id,nd,nc1)
c  writes summary of calculation on d/s iw
c
c  summaries are numbered consecutively, starting at nc1
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
      implicit double precision (a-h,o-z)
      dimension data(id,nd)
      common/cicode/ icow,iper,intmod,ial,istsbc,icc(4)
c
      save
c
      data icasep/-1/
c  unpack case number
      call decicd(icase,iel,iekin)
c
      ih1=ih
c  test for change in icase since last call
      if(icase.ne.icasep) ih1=1
      icasep=icase
c  new model?
      if(nwmod.ne.1) go to 6
      ih1=1
c  model type
      if(an.ge.0) go to 3
      write(iw,100)
      go to 4
    3 write(iw,105) an
c  truncated?
    4 if(xtrcnt.gt.0) write(iw,107) xtrcnt
      go to 7
c  first part of heading
    6 if(ih1.ne.1) go to 10
c  isothermal atmosphere?
    7 if(istsbc.ne.0) write(iw,108)
      if(icow.gt.0) then
        write(iw,110)
      else
        write(iw,115)
      end if
c  test for lambda.ne.1
   10 if(icow.eq.0) then
        inw=7
      else
        inw=8
      end if
      iord=2
      if(ial.eq.0) go to 15
      if(ih1.eq.1) write(iw,120)
      iord=3
      inw=inw+1
      go to 20
c  test for interpolation in model
   15 if(intmod.ne.1) go to 20
      if(ih1.eq.1) write(iw,125)
      inw=inw+1
      iord=3
c  test for periods
   20 if(iper.ne.1) go to 30
      if(ih1.ne.1) go to 25
      if(icow.gt.0) then
        write(iw,130)
      else
        write(iw,135)
      end if
   25 inw=inw+2
c  now heading is written
   30 if(nd.le.0) return
c  write data
      if(icow.le.0) then
        ifc=2*(iord-2)+1
      else
        ifc=2*(iord-2)+2
      end if
      if(ih1.eq.1) write(iw,140)
c
c  step in modes
c
      ncc=nc1
c
      if(iper.ne.1) go to 50
c
c  summaries with variational periods and frequencies
c
      do 40 n=1,nd
      iorder=data(iord,n)
      go to (32,34,36,38), ifc
   32 write(iw,142) data(1,n),iorder,(data(i,n),i=3,inw),ncc
      go to 40
   34 write(iw,144) data(1,n),iorder,(data(i,n),i=3,inw),ncc
      go to 40
   36 write(iw,146) data(1,n),data(2,n),iorder,(data(i,n),i=4,inw),ncc
      go to 40
   38 write(iw,148) data(1,n),data(2,n),iorder,(data(i,n),i=4,inw),ncc
c
   40 ncc=ncc+1
      return
c
c  summaries without variational periods and frequencies
c
   50 do 60 n=1,nd
      iorder=data(iord,n)
      go to (52,54,56,58), ifc
   52 write(iw,152) data(1,n),iorder,(data(i,n),i=3,inw),ncc
      go to 60
   54 write(iw,154) data(1,n),iorder,(data(i,n),i=3,inw),ncc
      go to 60
   56 write(iw,156) data(1,n),data(2,n),iorder,(data(i,n),i=4,inw),ncc
      go to 60
   58 write(iw,158) data(1,n),data(2,n),iorder,(data(i,n),i=4,inw),ncc
c
   60 ncc=ncc+1
c
      return
  100 format(//' physical model')
  105 format(//' polytrope with n=',f10.7)
  107 format(1h+,28x,', truncated at x =',f10.7)
  108 format(//' isothermal atmosphere')
  110 format(//8x,'l,order,(sig**2)sol.,(sig**2)corr.,ximax,',
     *  'xmax,ekin,period(sol.)')
  115 format(//8x,'l,order,sig**2,ximax,xmax,ekin,period(sol.)')
  120 format('+lambda,')
  125 format('+  zeta,')
  130 format(1h+,70x,',period(v.i.),freq.(v.i.)')
  135 format(1h+,50x,',period(v.i.),freq.(v.i.)')
  140 format(/)
  142 format(f14.8,i5,1pe18.10,e13.5,0pf8.5,1pe13.5,0p2f11.3,
     *  f11.7,22x,i5)
  144 format(f14.8,i5,1p2e18.10,e13.5,0pf8.5,1pe13.5,0p2f11.3,
     *  f11.7,4x,i5)
  146 format(f12.8,f12.6,i5,1pe18.10,e13.5,0pf10.7,1pe13.5,0p2f11.3,
     *  f11.7,10x,i5)
  148 format(f11.7,f11.5,i5,1p2e15.7,e13.5,0pf10.7,1pe13.5,0p2f11.3,
     *  f11.7,i5)
  152 format(f14.8,i5,1pe18.10,e13.5,0pf8.5,1pe13.5,0pf11.3,44x,i5)
  154 format(f14.8,i5,1p2e18.10,e13.5,0pf8.5,1pe13.5,0pf11.3,26x,i5)
  156 format(f12.8,f12.6,i5,1pe18.10,e13.5,0pf10.7,1pe13.5,0pf11.3,
     *  32x,i5)
  158 format(f11.7,f11.5,i5,1p2e15.7,e13.5,0pf10.7,1pe13.5,0pf11.3,
     *  22x,i5)
      end
c
c
