      program main
c
c  make fit as in Elsworth et al. to frequency separations.
c
      character fin*280
      parameter(nnmax=200)
      dimension data(4,nnmax), xfit(nnmax), yfit(nnmax)
c
      idiag=2
c
c  set data
c
      write(6,*) 'Enter input data file'
      read(5,'(a)') fin
      open(2,file=fin,status='old')
      call skpcom(2)
c
c  open output files
c
      open(10,file='ttt.fit-dnl',status='unknown')
      open(11,file='ttt.fit-dnl1',status='unknown')
c
      n=1
   10 read(2,*,end=20,err=20) (data(i,n),i=1,4)
      n=n+1
      go to 10
c
   20 nn=n-1
c
c  set up data for fit, do fit
c
      write(6,*) 'Enter limits in order, central order for fit'
      read(5,*) nord1, nord2, nordc
      lfin=length(fin)
      write(10,105) nord1, nord2, nordc, fin(1:lfin)
      write(11,107) fin(1:lfin), nord1, nord2, nordc 
c
      lp=-1
      nl=0
c
      do 30 n=1,nn
      l=data(1,n)+0.5
      nord=data(2,n)+0.5
c
      if(nord.ge.nord1.and.nord.le.nord2) then
	if(l.eq.lp) then
	  nl=nl+1
        else
c
c  new l-value. Test for making fit and printing output
c
	  if(nl.gt.0) then
	    nnl=nl
            call lstsqr(nnl,xfit,yfit,a,b,rms,sa,sb,idiag)
	    write(6,110) lp, a, sa, b, sb
	    write(10,110) lp, a, sa, b, sb
	    write(11,115) lp, a, sa, b, sb
          end if
          nl=1
	  lp=l
        end if
c
	xfit(nl)=nord-nordc
	yfit(nl)=data(4,n)
      end if
c
   30 continue
c
c  test for doing final analysis
c
      if(nl.gt.0) then
        nnl=nl
        call lstsqr(nnl,xfit,yfit,a,b,rms,sa,sb,idiag)
        write(6,110) lp, a, sa, b, sb
        write(10,110) lp, a, sa, b, sb
	write(11,115) lp, a, sa, b, sb
      end if
      stop
  105 format(' Least square fit to frequency separation.'//
     *       ' Range in order =',i3,'  -  ',i3,
     *       '. Central order in fit =',i3//
     *       ' Input file: ',a60/)
  107 format('# Input file, nord1, nord2, nordc'/
     *  '#  l, intercept, error, slope, error:'/ a, 3i5)
  110 format(/' l =',i2,': intercept =',f8.3,'(',f5.3,') microHz'/
     *        '    ',2x,'  slope     =',f8.3,'(',f5.3,') microHz')
  115 format(i3,1p4e13.5)
      end
