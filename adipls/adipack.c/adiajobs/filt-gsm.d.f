      program main
c
c  filter grand summary. Various possibilities:
c
c  a) eliminate modes with excessive
c     difference between variational and Richardson frequencies.
c
c  b) select only modes with a specific case number (e.g. to
c     exclude modes that do not use isothermal atmosphere
c     solution.
c  c) exclude modes with a specific case number.
c
c  Double precision version
c  ++++++++++++++++++++++++
c
c  Dated 11/4/90
c
      implicit double precision (a-h, o-z)
      character*180 fin, fout
      dimension cs(50),ics(7),iss(2)
      equivalence (cs(39),ics(1)),(cs(6),iss(2))
c
c  set up data sets
c
      write(6,*) 'Enter input file'
      read(5,'(a)') fin
      open(2,file=fin,status='old',form='unformatted')
c
      write(6,*) 'Enter output file'
      read(5,'(a)') fout
      open(12,file=fout,status='unknown',form='unformatted')
c
c  set up case
c
      write(6,110)
      read(5,*) icase
c
      if(icase.eq.1) then
        write(6,*) 'Enter maximum allowed difference (microHz)'
        read(5,*) dfrmax
      else
        write(6,*) 'Enter case number for selection'
        read(5,*) icssel
      end if
c
c  start stepping through dataset
c
      nmd=0
c
   30 read(2,end=90,err=90) cs
      nmd=nmd+1
c
      if(icase.eq.1) then
        if(abs(cs(27)-cs(37)).gt.0.001d0*dfrmax) then
c
c  error message
c
          l=cs(18)+0.5d0
          nord=cs(19)+0.5d0
          write(6,120) l, nord, 1000*cs(27), 1000*cs(37)
c
        else
c
c  output mode
c
          write(12) cs
        end if
c
      else if((icase.eq.2.and.icssel.ne.ics(5))
     *    .or.(icase.eq.3.and.icssel.eq.ics(5))) then
c
c  error message
c
        l=cs(18)+0.5d0
        nord=cs(19)+0.5d0
        write(6,130) l, nord, 1000*cs(27), ics(5)
c
      else
c
c  output mode
c
        write(12) cs
c
      end if
c
      go to 30
c
   90 continue
c
      stop
  110 format(/
     *  ' Options:'/
     *  ' 1: Select on difference between ',
     *  'variational and Richardson frequencies'/
     *  ' 2: Select on case number to be included'/
     *  ' 3: Select on case number to be excluded'/
     *  '    Enter option:')
  120 format(' mode skipped. l, order =',2i4,' nu.var , nu.rich =',
     *  2f11.4,' microHz')
  130 format(' mode skipped. l, order, nu =',2i4,f8.2,
     *  ' microHz. case =',i10)
      end
