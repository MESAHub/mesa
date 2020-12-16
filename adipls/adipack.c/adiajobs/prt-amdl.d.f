      program main
c   
c  print adiabatic model variables
c   
      implicit double precision (a-h,o-z)
      parameter(nnmax=5000)
      character fin*280
      dimension x(nnmax),aa(6,nnmax),data(8)
c
      write(6,*) 'Enter file name'
      read(5,'(a)') fin
c   
      open(2,file=fin,status='old',
     *  form='unformatted')
c   
   10 read(2,end=90,err=90) nmod,nn,data,(x(n),(aa(i,n),i=1,5),n=1,nn)
c
      write(6,100) fin,(data(i),i=1,7)
      if(data(8).ne.1) then
	write(6,102)
      else
	write(6,103)
      end if
      write(6,105)
c   
      write(6,110) (n,x(n),(aa(i,n),i=1,5),n=1,nn) 
c
      go to 10
c
   90 continue
      stop
c   
  100 format('# file: ',a/'#'/'# data:'/'#',1p7e11.3/'#')
  102 format('# old, error-prone setting of r/R')
  103 format('# new setting of r/R')
  105 format('#'/'#  n, x, aa(1-5):'/'#')
  110 format(i5,0pf12.8,1p5e12.5)
      end   
