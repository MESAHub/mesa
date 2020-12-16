      subroutine skpcom(ids)
c
c  skips initial comment lines on unit ids.
c  a comment is flagged by a "#" in column no 1.
c
c  on exit file is positioned at first non-comment line
c
      character*1 i,ic
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      data ic /'#'/
c
   10 read(ids,100,end=90,err=90) i
      if(i.eq.ic) go to 10
c
      backspace ids
      return
c
c  diagnostics for errors
c
   90 if(istdpr.gt.0) 
     *  write(istdpr,*) 'error or EOF detected in s/r skpcom on unit',
     *  ids
      return
  100 format(a)
      end
