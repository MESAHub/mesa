      program main
c
c  reset file of adiabatic modes. the variables included depend on icase
c  (input in namelist /exec/). when icase = 1 (the default) the variables
c  included are y(1) and y(2). when icase = 2 the variables are
c  zhat(1) and zhat(2), defined such that the energy is proportional to
c
c      integral((zhat(1)**2 + zhat(2)**2) dr/r)
c
c  the first record of the output dataset contains the mesh in x
c  which is assumed to be the same for all modes.
c  the subsequent records contain cs (the grand summary) and
c  the chosen variables.
c
c  Double precision version
c  ++++++++++++++++++++++++
c
c  Dated 14/5/90
c
      implicit double precision (a-h, o-z)
      character*280 ccase, fin, fout
      dimension x(2410),cs(50),yy(6,2410),ccase(2)
      data ccase /'set y(1) and y(2)','set zhat(1) and zhat(2)'/
c
c$nl      namelist /exec/ icase
c
      icase=1
c
c$nl      read(5,exec,end=8)
c
      write(6,*) 'Enter input file'
      read(5,'(a)') fin
      write(6,*) 'Enter output file'
      read(5,'(a)') fout
c
      open(4,file=fin,status='old',form='unformatted')
      open(20,file=fout,status='unknown',form='unformatted')
c
      write(6,*) 'Enter case number (1 for y1, y2; 2 for zh1, zh2)'
      read(5,*) icase
c
      if(icase.ne.2) icase=1
c
    8 nmde=0
c
c  set range in variables
c
      iv1=1+4*(icase-1)
      iv2=iv1+1
c
      write(6,100) ccase(icase)
c
   10 read(4,end=30,err=30) cs,nn,(x(n),(yy(i,n),i=1,6),
     *  n=1,nn)
c
      nmde=nmde+1
c
c  test for output of nn and x
c
      if(nmde.eq.1) then
        write(20) nn,(x(n),n=1,nn)
      end if
c
c  output of restricted variables
c
      write(20) cs,((yy(i,n),i=iv1,iv2),n=1,nn)
      write(6,110) nmde,(cs(i),i=18,20),cs(27)
      go to 10
c
   30 continue
      stop
  100 format('1 ',a//' nmode, l, order, sigma**2, nu (millihz):'/)
  110 format(i5,2f10.1,f12.3,f12.5)
      end
