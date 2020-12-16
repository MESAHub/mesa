      subroutine rdamdl(ids,x,aa,data,nn,nmod,ivar,icry,ia)
c
c  reads amdl file from unit ids checking (if this is the first read
c  from the given unit) for number of variables from value of
c  data(8)
c
c  ivar returns the number of variables (which is also used in 
c  subsequent reads). icry is returned as 0 for successful read,
c  as 1 on end of file and as 2 on error in read.
c
c  Original version: 2/8/95
c
      implicit double precision(a-h, o-z)
      dimension x(*), aa(ia,*), data(*)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      data idsp /-1/
c
      save
c
c  open file for model input
c
      if(ids.ne.idsp) then 
	if(idsp.gt.0) close(idsp)
	call openf(ids,'o','u')
c
c  read start of file, to check for number of variables
c
        read(ids) nmod,nn,(data(i),i=1,8)
        close(ids)
	call openf(ids,'o','u')
	idata8 = idint(data(8)+0.1)
	if(idata8.ge.100) then
	  ivar = 8
        else if(idata8.ge.10) then
	  ivar = 6
        else
	  ivar = 5
        end if
	ivarc = ivar
	idsp = ids
c
      else
c
	ivar = ivarc
c
      end if
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'In rdamdl ivar =',ivar
c
      read(ids,end=80, err=90) nmod,nn,(data(i),i=1,8),
     *  (x(n),(aa(i,n),i=1,ivar), n=1,nn)
c
      icry = 0
      return
c
   80 write(istdpr,110) ids
      icry = 1
      return
c
   90 write(istdpr,120) ids
      icry = 2
      return
  110 format(/' ***** EOF reached on unit',i4,' in s/r rdamdl')
  120 format(/' ***** Error on unit',i4,' in s/r rdamdl')
      end
