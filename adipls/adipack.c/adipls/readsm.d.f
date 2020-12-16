      subroutine readsm(iw,ih,nwmod,an,xtrnct,icas1,data,id,nd,ndmax,
     *  idiag)
c  reads summary from d/s iw and sets case number in icas1(n), data in
c  data(i,n), n=1,...,nd, where nd is the number of data records or
c  ndmax (set in call), whichever is the smallest.
c  if heading for new model is encountered, nwmod is set to 1,
c  and if any heading is encountered, ih is set to 1.
c  for idiag = 1  each line is printed.
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
      implicit double precision (a-h,o-z)
      integer dot,blank,a,peh,el
      dimension a(132),data(id,ndmax),icas1(ndmax),iaper(4),inw1(4)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      data dot,blank,peh,el /1h.,1h ,1hp,1hl/
      data iaper,inw1 /87,105,99,111,6,7,7,8/
c
      nr=0
      n=0
      ih=0
      nwmod=0
c  for the moment do not set an,xtrnct from data
      an=-1
      xtrnct=0
c  new record. empty a
   10 do 15 j=1,132
   15 a(j)=blank
c  read record
      read(iw,100,end=70,err=70) a
      nr=nr+1
c  find case
   20 if(a(6).ne.dot) go to 25
      if(a(53).ne.dot) go to 22
      ifc=1
      icow=0
      ial=0
      go to 30
   22 if(a(59).ne.dot) go to 90
      ifc=2
      icow=1
      ial=0
      go to 30
c
   25 if(a(4).ne.dot) go to 50
      if(a(18).ne.dot) go to 27
      ifc=3
      icow=0
      ial=1
      go to 30
   27 if(a(17).ne.dot) go to 90
      ifc=4
      icow=1
      ial=1
c  now case is determined. periods?
   30 iper=0
      j=iaper(ifc)
      if(a(j).eq.dot) iper=1
      n=n+1
      icas1(n)=1000*ial+10*iper+icow
c  number of variables
      inw=inw1(ifc)+3*iper
c  read data
      backspace iw
      go to (31,32,33,34), ifc
   31 read(iw,111) data(1,n),iorder,(data(i,n),i=3,inw)
      data(2,n)=iorder
      go to 40
   32 read(iw,112) data(1,n),iorder,(data(i,n),i=3,inw)
      data(2,n)=iorder
      go to 40
   33 read(iw,113) data(1,n),data(2,n),iorder,(data(i,n),i=4,inw)
      data(3,n)=iorder
      go to 40
   34 read(iw,114) data(1,n),data(2,n),iorder,(data(i,n),i=4,inw)
      data(3,n)=iorder
c  now data is set. testoutput
   40 if(idiag.eq.1.and.istdpr.gt.0) 
     *  write(istdpr,120) nr,n,icas1(n),(data(i,n),i=1,inw)
      if(n.lt.ndmax) go to 10
c  n = ndmax. finish
      nd=n
      return
c  not a data card. print card
   50 if(idiag.eq.1.and.istdpr.gt.0) write(istdpr,130) (a(i),i=1,120),nr
c  set nwmod or ih?
      if(a(2).eq.peh) nwmod=1
      if(a(9).eq.el) ih=1
      go to 10
c  end of data reached. test that data has been read
   70 nd=n
      if(n.gt.0) return
      if(istdpr.gt.0) write(istdpr,140) iw,nr
      return
c  internal inconsistency detected. write diagnostic
   90 if(istdpr.gt.0) write(istdpr,150) nr,iw,(a(i),i=1,131)
c  try next card
      go to 10
  100 format(132a1)
  111 format(f14.8,i5,1p e18.10,e13.5,0pf8.5,1pe13.5,0p3f11.6)
  112 format(f14.8,i5,1p2e18.10,e13.5,0pf8.5,1pe13.5,0p3f11.6)
  113 format(f12.8,f12.6,i5,1pe18.10,e13.5,0pf10.7,1pe13.5,0p3f11.6)
  114 format(f11.7,f11.5,i5,1p2e16.8,e13.5,0pf10.7,1pe13.5,0p3f11.6)
  120 format(3i6,1p11e10.2)
  130 format(1x,120a1,i6)
  140 format(//1x,10(1h*),' no data on d/s',i3,'.  ',i6,
     *  ' records read')
  150 format(//1x,10(1h*),i6,'-th record on d/s',i3,
     *  ' inconsistent. the record:'//1x,131a1/1x,
     *  13(10h1234567890))
      end
c
c
