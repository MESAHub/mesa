c this is lineexpop.f but for linronvvf.dat (ron + Verner-Ferland)
c ------------------------------------------------------------------------ 
c                s u b r o u t i n e   l i n e e x p o p 
c ------------------------------------------------------------------------ 
c The routine computes the opacity contribution due to transitions between 
c levels in ions whose population are computed using a Saha-Boltzmann 
c equation of state. The opacity computed is an `expansion opacity',  
C computed using the formulae given in Eastman and Pinto (1993).  
 
c As is, the opacities generated are purely scattering--contributions are 
c  added to SCATOP, and OPACITY is returned unchanged. If the routine 
c  is preprocessed with ABSORPTIVE_LINE_OPAC defined, then it will use the 
c  approximation described by L.S. Anderson, ApJ. 339, p558, 1989, which 
c  uses the Van Regemorter formula to compute the probability that an 
c  upper level will be collisionally depopulated before spontaneously 
c  decaying. In this case, the fraction of the line opacity which is 
c  absorptive and scattering is estimated, and contributions are added 
c  to both opacities. The emissivity is also incremented appropriately. 
 
c The routine can also be used to generate a `short' list of  
c  lines from a much longer list, by giving it a appropriate 
c  atomic model and a minimum Sobolev optical depth which a line 
c  must possess in order to be included in the shortened list of lines. 
 
      subroutine lineexpop(linelist, readdump, nwave, mxwave, 
     .     hfreq3c2, dlnfreq, wlgrd, scatop, opacity, emis, 
     .     plnkfnc, init, mkshtlst, taumin, longlist, 
     .     nradii, nrdim, eden, temp, dvdravg, iondenz, ionstg1, 
     .     nr1, nr2) 
 
      implicit real*8 (a-h,o-z) 
 
 

      PARAMETER(MAXIONM1=(6  - 1))
 

      PARAMETER(NSBINTVL=3*10 )
 


      data pi/3.141592653589793d+00/, fourpi/12.5637061d+00/
      data bc/1.380626d-16/, h/6.626205d-27/, hev/4.1357d-15/
      data c/2.997925d+10/, elecxsec/6.6524d-25/
      data evtoerg/1.6022d-12/, a/7.56464d-15/
      data stefbltz/5.66956d-05/, hmass/1.67352d-24/
      data esu/4.80298d-10/, emass/9.1091d-28/
      data srpi/1.772453851d+00/, bcev/8.617064d-05/


 
c Calling arguments: 
 
c (character*80 ) LINELIST is the name of an ascii file containing the line  
c list. 
 
c (logical) READDUMP specifies whether the line list should be read 
c  from a previously written binary image by the name of linelist. 
c  T=> read line data from linelist. F => write a binary image to 
c  `linedata.dump'. 
 
c (integer) NWAVE is the number of wavelengths in the wavelength/frequency 
c  grid. 
 
c (integer) MXWAVE is the first (frequency grid) dimension of SCATOP, 
c OPACITY and EMIS. 
 
c (real*8) HFREQ3C2(NWAVE) = h * freq**3 / c**2 on the frequency grid. 
 
c (real*8) DLNFREQ(NWAVE) = 2 (freq(k+1)-freq(k))/(freq(k)+freq(k+1)). 
 
c (real*8) WLGRD(NWAVE) contains the wavelength grid, in angstroms, with 
c  WLGRD(K+1)<WLGRD(K). 
 
c (real*8) SCATOP(MXWAVE,NRADII), OPACITY(MXWAVE,NRADII) and  
c  EMIS(MXWAVE,NRADII) contain the scattering opacity, absorptive opacity, 
c  and emissivity, respectively, tabulated as a function of frequency and 
c  radius. 
 
c (real*8) PLNKFNC(MXWAVE,NRADII) must contain the Planck function, B-NU. 
 
c (integer) INIT should be set to zero the first time LINEEXPOP is called. 
c  This causes the line list to be read and sorted, and various other 
c  quantities to be initialized. If the frequency grid is changed, then 
c  upon the very next call to LINEEXPOP, init should be set to zero. It is 
c  set to 1 upon return. 
 
c (logical) MKSHTLST is a flag which, if true, causes LINEEXPOP to read in 
c  the line list contained in the ascii file, LONGLIST, and write out to  
c  the file LINELIST all lines which have TAUSOB > TAUMIN, based on the  
c  current model. 
 
c (real*8) TAUMIN is the minimum Sobolev optical depth a line should have, 
c  somewhere in the star, in order to be included in the file LINELIST (used 
c  when MKSHTLST is true. 
 
c (character*80 ) LONGLIST is the name of a (large) file containing data 
c  on more lines than we want to know about. 
 
c (integer) NRADII is the number of radial grid points. 
 
c (integer) NRDIM is the first dimension of IONDENZ and IONSTG1. 
 
c (real*8) EDEN(NRADII) should contain the electron density on the 
c  radial grid. 
 
c (real*8) TEMP(NRADII) should contain the temperature on the radial 
c  grid. 
 
c (real*8) DVDRAVG(NRADII) should contain the angle averaged velocity  
c  gradient on the radial grid. 
 
c (real*8) IONDENZ(NRDIM,6,99) should contain the ion density divided by 
c  the partition function of the six most populated successive ionization 
c  stages of every element up to Z=99. 
 
c (integer) IONSTG1(NRDIM,99) contains the index of the first ionization 
c  stage, as a function of zone and Z, which is stored in IONDENZ. Thus, 
c  in zone NR for element IZ, IONDENZ(NR,1,IZ) contains the IONSTG1(NR,IZ)th 
c  ionization stage, IONDENZ(NR,2,IZ) contains the (IONSTG1(NR,IZ)+1)th 
c  ionization stage, etc. 
 
c (integer) NR1 and NR2 are the indecies of the first and second zone over which 
c  to compute the opacity. 
 
      character*80  linelist, longlist 
 
      logical readdump, mkshtlst 
 
      integer ionstg1(nrdim,99) 
 
      logical wrtdump, fexist 
 
      integer nwave, mxwave, nradii, nrdim, init, nr1, nr2 
 
      real*8 wlgrd(mxwave), taumin 
      real*8 dlnfreq(mxwave), hfreq3c2(mxwave) 
      real*8 scatop(mxwave,nradii), opacity(mxwave,nradii) 
      real*8 emis(mxwave,nradii) 
      real*8 plnkfnc(mxwave,nradii) 
      real*8 eden(nradii), temp(nradii), dvdravg(nradii) 
      real*8 iondenz(nrdim,6,99) 
 
c Local variables 
      pointer (pdvdrlog, dvdrlog) 
      pointer (ptrwavln, wavelen) 
      pointer (ptgf, gf) 
      pointer (ptelow, elower) 
      pointer (pteup, eupper) 
      pointer (ptenrwt, enrwt) 
      pointer (ptznucln, znucline) 
      pointer (ptionstg, ionstage) 
      pointer (ptindexv, indexv) 
      pointer (ptlteflg, lteflgx) 
      pointer (ptenrndx, enrindex) 
      pointer (ptwlindx, wlindex) 
      pointer (ptenrvec, enrvec) 
      pointer (ptrlvndx, levindex) 
      pointer (ptrhckt, hckt) 
      pointer (ptrlabund, labund) 
      pointer (ptrexpvc, expvecx) 
      real*8 expvecx(75 *200 ) 

 
      pointer (ptrlscat, lscat) 

      pointer (ptrlopac, lopac) 

      pointer (ptrlemis, lemis) 

      real*8 lscat(100000000), lopac(100000000) 
      real*8 lemis(100000000) 
      pointer (ptridenz, idenz) 
      real*8 idenz(31*99 *75 ) 

 

      pointer (ptrlcountr, linecountr) 
      integer linecountr(4950) 

      pointer (ptrlcountw, linecountw) 
      integer linecountw(4950) 
 
      real*8 enrvec(200 ) 
 

 
      real*4 wavelen(112219 ), gf(112219 ) 
      real*4 elower(112219 ), eupper(112219 ) 
      real*4 enrwt(112219 ) 
      real*4 dumvec(112219 ), element(112219 ) 
      real*8 dvdrlog(75 ), hckt(75 ) 
      logical labund(99 ) 
 
      integer idumvec(112219 ), levindex(112219 ) 
      integer znucline(112219 ), ionstage(112219 ) 
      integer indexv(112219 ), lteflgx(112219 ) 
      integer enrindex(112219 ), wlindex(112219 ) 
 
      pointer (pdumvec, dumvec) 
      pointer (pidmvec, idumvec) 
      pointer (pelement, element) 
 
      real*4 jupper, jlower, tenlog 
 
      character*4 ref 
      character*8 lnamelw, lnameup 
 
      integer ptrinit 
 
      common /ltelnz/ptrwavln, pelement, pidmvec, pdumvec , ptznucln, ptenrwt, ptenrvec, 
     .     ptgf, ptelow, ptenrndx, ptionstg, ptwlindx, 
     .     ptrlvndx, nlines
 
      real*8 quench_fac(2) 
      integer discard1, discard2, type1, type2, type3, type4 
 
      save ptrinit 
      data c18/2.997925d+18/, ptrinit/0/, smallnum/1.d-70/ 
      data ediffmax/1.d-01/ 
 
c QENCHFAC is the electron quenching correction factor computed by 
c L. Anderson, ApJ. 399, p558, 1989. I have rescaled Anderson's values  
c so as to be mulitiplied by 2 * H * FREQ(NF)**2 / C**2, with everything  
c in cgs units. 
      data quench_fac/2.9d+17, 3.4d+16/ 
 
c  ------------------------------------------------------------------ 
      iiii=1  
      saha = 0.5d+00 * ((h / emass) * (h / bc) / (2.d+00 * pi))**1.5 
      cross0 = pi * (esu / emass) * (esu / c) * 1.d-08 
 
        
 
      if (mkshtlst) call mzalloc('pdvdrlog ',pdvdrlog, nradii * 8 ) 
      if (ptrinit .eq. 0) then 
         call mzalloc('ptrhckt ', ptrhckt, nradii * 8 ) 
      end if 
      call mzalloc('ptridenz ', ptridenz, 31 * nradii * 99 * 8 ) 
 
C Creating a full array of 31 possible ionization stages for every element 
c obviates the necessity for evaluating 50 million conditionals inside of 
c the opacity calculation do loop. 
 
      do iz = 1, 99 
c Zero idenz. 
         do nr = 1, nradii 
            idenz(nr+(1-1)*nradii+(iz-1)*31*nradii)  = 0.d+00 
         end do 
         do i = 2, 31 
          call dcopy(nradii, 
     *idenz(1+(i-1-1)*nradii+(iz-1)*31*nradii) ,1, 
     *idenz(1+(i-1)*nradii+(iz-1)*31*nradii) ,1)  
         end do 
 
         do nr = 1, nradii 
            do ii = 1, 6 
               i =  ionstg1(nr,iz) + ii - 1 
               idenz(nr+(i-1)*nradii+(iz-1)*31*nradii)  
     * = iondenz(nr,ii,iz) 
            end do 
         end do 
 
      end do 
c          <- iz. 
 
        
 
c Compute hc / kt 
 
      do nr = 1, nradii 
         hckt(nr) = h * c / (bc * temp(nr)) 
      end do 
 
c If a dump file isn't being read in, write one out. 
      wrtdump = (.not. readdump) 
 
      if (mkshtlst) then 
c Make a `short' list of the strongest lines. 
 
         inquire(file=longlist,exist=fexist) 
         print *,'LongList=',longlist 
         print *,'linelist=',linelist 
         if (.not. fexist) then 
            write(7,'(3a)') ' lineexpop: unable to', 
     .           ' find file ', linelist(1:lnblnk(linelist)) 
            stop 
         end if 
 
         open(file=longlist,unit=20,status='old') 
         open(file=linelist,unit=21,status='unknown') 
 
         tsminln = dlog(cross0 * 10.d+00 / taumin) 
 
 
 
 
         do 3000 iz = 1, 99 
 
            do 2620 i = 1, min(iz,6) 
 
               do 2610 nr = 1, nradii 
                  idenz(nr+(i-1)*nradii+(iz-1)*31*nradii)  = dlog(max( 
     *idenz(NR+(I-1)*nradii+(IZ-1)*31*nradii) , 1.d-99)) 
 2610          continue 
 
 2620       continue 
 
 3000    continue 
 
         do 3010 nr = 1, nradii 
            dvdrlog(nr) = dlog(dabs(dvdravg(nr))) 
 3010    continue 
 
 
c Allocate space for linecountr and linecountw, used to keep statistics 
c on the number of lines in each ionization stage read and written. 
         call mzalloc('ptrlcountr ', ptrlcountr, 4950 * 8 ) 
         call mzalloc('ptrlcountw ', ptrlcountw, 4950 * 8 ) 
 
         do l = 1, 4950 
            linecountr(l) = 0 
         end do 
 
         do l = 1, 4950 
            linecountw(l) = 0 
         end do 
 
         nlcount = 0 
 
         do 3130 il = 1, 9999999 
 
c The first entry in the ``long'' line list is read in first, and the log of 
c it's wavelength is computed by a call to WLLOG. 
 
            read(20,30,end=3140) wl, ref, gflog, 
     .           jlower, el, jupper, eu, code, lnamelw, lnameup 
            iz = int(code + 0.001) 
            i = int(100.d+00 * (code - dfloat(iz)) + 0.1d+00) + 1 
            linecountr(i+((iz-1)*iz)/2)  
     * = linecountr(i+((iz-1)*iz)/2)  + 1 
 
            el = dabs(el) 
            eu = dabs(eu) 
 
            if (el .gt. eu) then 
               oldeu= eu 
               eu = el 
               el = oldeu 
               oldju= jupper 
               jupper = jlower 
               jlower = oldju 
            end if 
 
c Currently, the ionization is computed (elsewhere) by routine SAHAEQN, 
c which can only handle the first six ionization stages. Thus, here the 
c decision is made to limit the line list to transitions of the first 
c six ionization stages of any element. 
 
            wllog = dlog(wl) 
            denmaxln = idenz(1+(i-1)*nradii+(iz-1)*31*nradii)  - el * hc
     ~kt(1) - dvdrlog(1) 
 
            do nr = 2, nradii 
               denmaxln = max(denmaxln, idenz(nr+(i-1)*nradii+(iz-1)*31*
     ~nradii)  - 
     .              el * hckt(nr) - dvdrlog(nr))  
            end do 
 
            if ((tsminln + denmaxln + gflog + wllog) .gt. 0.d+00) then 
               nlcount = nlcount + 1 
               write(21,30) wl, ref, gflog, jlower, el, jupper, 
     .              eu, code, lnamelw, lnameup 
               linecountw(i+((iz-1)*iz)/2)  = linecountw(i+((iz-1)*iz)/2
     ~)  + 1 
            end if 
 
            wlo = wl 
 
c The rest of the entries are read in one at a time, in the same fashion 
c as above, but the log of the wavelength is approximated without the call 
c to dlog. 
 
            do nl= 1, 300 
               read(20,30,end=3140) wl, ref, gflog, 
     .              jlower, el, jupper, eu, code, lnamelw, lnameup 
               dlnwl = 2.d+00 * (wl - wlo) / (wl + wlo) 
               wllog = wllog + dlnwl 
               wlo = wl 
               iz = int(code + 0.001) 
               i = int(100.d+00 * (code - dfloat(iz)) + 0.1d+00) + 1 
               denmaxln = idenz(1+(i-1)*nradii+(iz-1)*31*nradii)  - el *
     ~ hckt(1) - dvdrlog(1) 
 
               do nr = 2, nradii 
                  denmaxln = max(denmaxln, idenz(nr+(i-1)*nradii+(iz-1)*
     ~31*nradii)  - 
     .                 el * hckt(nr) - dvdrlog(nr))  
               end do 
 
              if ((tsminln + denmaxln + gflog + wllog) .gt. 0.d+00) then
     ~ 
                  nlcount = nlcount + 1 
                  write(21,30) wl, ref, gflog, jlower, el, jupper, 
     .                 eu, code, lnamelw, lnameup 
               end if 
 
            end do 
c                <- nl. 
 3130    continue 
 
c --------------- 
 3140    continue 
c --------------- 
 
         write(7,'(a,a,a,a,a,i8)')  
     .        ' Number of lines written from ',  
     .        longlist(1:lnblnk(longlist)), 
     .        ' to ', linelist(1:lnblnk(linelist)), 
     .        ' equals ', nlcount 
         write(7,*) 
         write(7,'(a)') ' Lines per ionization stage read and written:' 
     ~
         write(7,'(a)') '    z    i      read   written' 
         do iz = 1, 99 
            do i = 1, iz 
               if (linecountr(i+((iz-1)*iz)/2)  .gt. 0) then 
                  write(7,'(1x,i4,1x,i4,1x,i9,1x,i9)') iz, i, 
     .                 linecountr(i+((iz-1)*iz)/2) , linecountw(i+((iz-1
     ~)*iz)/2)  
               end if 
            end do 
         end do 
 
         write(7,*) 
c        ---- 
         stop 
c        ---- 
 
      end if 
c          <- if (mkshtlst) 
 
c Allocate seperate space for temporary storage of scattering opacity.  
c Scattering opacity is computed by spacial zone, but will be transfered 
c to arrays where it is stored by frequency and then spacial zone. The 
c same holds for absorptive opacity and emissivity, below. 
          call mzalloc('ptrlscat ', ptrlscat, 8  * nwave * nradii) 
          call mzalloc('ptrlopac ', ptrlopac, 8  * nwave * nradii) 
 
      do nr = 1, nradii 
         lscat(nr+(1-1)*nradii)  = 0.d+00 
         lopac(nr+(1-1)*nradii)  = 0.d+00
      end do
 
      do nf = 2, nwave 
         call dcopy(nradii, lscat(1+(nf-1-1)*nradii) , 1, lscat(1+(nf-1)
     ~*nradii) , 1) 
         call dcopy(nradii, lopac(1+(nf-1-1)*nradii) , 1, lopac(1+(nf-1)
     ~*nradii) , 1) 
      end do 

 
        
       
c     -------------------------> 
      if (init .eq. 0) then 
 
         if (readdump) then 
 
c Read lines from binary dump file of Kurucz line list. 
 
            inquire(file=linelist, exist=fexist) 
            print * , ' Linelist ' , linelist 
            if (.not. fexist) then 
               write(7,'(3a)') ' lineexpop: readdump=true: unable to', 
     .              ' find file: ', linelist(1:lnblnk(linelist)) 
             print * ,' program stoped' 
               stop 
            end if 
             print * ,' program continues , file exists ',linelist 
 
c                <- if (.not. fexist) 
 
            open(file=linelist,form='unformatted',unit=35, 
     .           status='old') 
            read(35) nlines 
 
            if (ptrinit .eq. 0) then 
               call mzalloc('ptrwavln ',ptrwavln, nlines * 8 ) 
               call mzalloc('ptgf ',ptgf, nlines * 8 ) 
               call mzalloc('ptelow ',ptelow, nlines * 8 ) 
               call mzalloc('ptenrwt ',ptenrwt, nlines * 8 ) 
               call mzalloc('ptznucln ',ptznucln, nlines * 8 ) 
               call mzalloc('ptrlvndx ', ptrlvndx, nlines * 8 ) 
               call mzalloc('ptionstg ',ptionstg, nlines * 8 ) 
               call mzalloc('ptenrndx ',ptenrndx, nlines * 8 ) 
               call mzalloc('ptwlindx ',ptwlindx, nlines * 8 ) 
               call mzalloc('ptenrvec ',ptenrvec, 200  * 8 ) 
               call mzalloc('ptrexpvc ',ptrexpvc, 200  * nradii * 8 ) 
            end if 
c                <- if (ptrinit .eq. 0) 
 
            call mzalloc('ptindexv ',ptindexv, nlines * 8 ) 
            call mzalloc('ptlteflg ',ptlteflg, nlines * 8 ) 
            call mzalloc('pteup ',pteup, nlines * 8 ) 
            call mzalloc('pdumvec ',pdumvec, nlines * 8 ) 
            call mzalloc('pidmvec ',pidmvec, nlines * 8 ) 
            call mzalloc('pelement ',pelement, nlines * 8 ) 
 
            call readbinr(35, elower, nlines) 
            call readbinr(35, eupper, nlines) 
            call readbinr(35, gf, nlines) 
            call readbini(35, ionstage, nlines) 
            call readbinr(35, wavelen, nlines) 
            call readbini(35, znucline, nlines) 
            call readbinr(35, element, nlines) 
            close(unit=35) 
 
c This shouldn't be necessary, but Kurucz tagged some of his energy 
c levels by making them negative, and these have made there way into 
c binary dumps. 
 
            do nl = 1, nlines 
               elower(nl) = abs(elower(nl)) 
               eupper(nl) = abs(eupper(nl)) 
            end do 
 
              
 
         else 
 
c Read in line data from Kurucz line list (ascii file). 
 
            inquire(file=linelist,exist=fexist) 
 
            if (.not. fexist) then 
               write(7,'(3a)') ' lineexpop: unable to', 
     .              ' find file ', linelist(1:lnblnk(linelist)) 
               stop 
            end if 
c                <- if (.not. fexist) 
 
            open(file=linelist,unit=20,status='old') 
 19         format(i10) 
            read(20,19) nlines 
            write(*,*) ' nlines=',nlines
 
            if (ptrinit .eq. 0) then 
               call mzalloc('ptrwavln ', ptrwavln, nlines * 8 ) 
               call mzalloc('ptgf ', ptgf, nlines * 8 ) 
               call mzalloc('ptelow ', ptelow, nlines * 8 ) 
               call mzalloc('ptenrwt ', ptenrwt, nlines * 8 ) 
               call mzalloc('ptznucln ', ptznucln, nlines * 8 ) 
               call mzalloc('ptrlvndx ', ptrlvndx, nlines * 8 ) 
               call mzalloc('ptionstg ', ptionstg, nlines * 8 ) 
               call mzalloc('ptenrndx ', ptenrndx, nlines * 8 ) 
               call mzalloc('ptwlindx ', ptwlindx, nlines * 8 ) 
               call mzalloc('ptenrvec ', ptenrvec, 200  * 8 ) 
               call mzalloc('ptrexpvc ', ptrexpvc, 200  * nradii * 8 ) 
            end if 
c                <- if (ptrinit .eq. 0) 
 
            call mzalloc('ptindexv ', ptindexv, nlines * 8 ) 
            call mzalloc('ptlteflg ', ptlteflg, nlines * 8 ) 
            call mzalloc('pteup ',pteup, nlines * 8 ) 
            call mzalloc('pdumvec ',pdumvec, nlines * 8 ) 
            call mzalloc('pidmvec ',pidmvec, nlines * 8 ) 
            call mzalloc('pelement ',pelement, nlines * 8 ) 
 
            do 20 nl= 1, nlines 
               read(20,30) wavelen(nl), ref, gf(nl), 
     .              jlower, elower(nl), jupper, eupper(nl), 
     .              element(nl), lnamelw, lnameup 
 30            format(f10.4,1x,a4,1x,f7.3,f4.1,f11.3,f4.1, 
     .              1x,f11.3,f7.2,a8,2x,a8) 
c ZNUCLINE and IONSTAGE are the nuclear charge and ionization stage, 
c respectively, to which the transition belongs. 
               znucline(nl) = int(element(nl) + 0.001) 
               ionstage(nl) = int(100.d+00 * (element(nl) -  
     .              dfloat(znucline(nl))) + 0.1d+00) + 1 
               eupper(nl) = abs(eupper(nl)) 
               elower(nl) = abs(elower(nl)) 
 20         continue 
c                  <- end sum over input records, nl. 
 
            close(20) 
 
            if (wrtdump) then 
               open(file='linedata.dump',form='unformatted',unit=35, 
     .              status='unknown') 
               write(35) nlines 
               call writbinr(35, elower, nlines) 
               call writbinr(35, eupper, nlines) 
               call writbinr(35, gf, nlines) 
               call writbini(35, ionstage, nlines) 
               call writbinr(35, wavelen, nlines) 
               call writbini(35, znucline, nlines) 
               call writbinr(35, element, nlines) 
               close(unit=35) 
            end if 
c                <- if (wrtdump) 
         end if 
c              <- if (readdump) 
 
         write(7,*) ' lineexpop: number of lines read in = ',nlines 
 
         enrmax = 0.d+00 
c Determine maximum energy of lower level, of all lines. 
 
         do nl = 1, nlines 
            if (elower(nl) .ge. enrmax) enrmax = elower(nl) + 1.d+00 
         end do 
 
         tenlog = log(10.e+00) 
 
         do nl = 1, nlines 
c Convert GF from log(gf) to actual gf-value. 
            gf(nl) = exp(tenlog * gf(nl)) 
c Convert wavelength from nanometers to angstroms. 
            wavelen(nl) = 10.e+00 * wavelen(nl) 
         end do 
 
c Throw out those lines which don't fit inside the included wavelength 
c range. 
 
         wlmin = wlgrd(nwave) 
         wlmax = wlgrd(1) 
 
         if (dble(wavelen(1)) .gt. wlmin) then 
            nl1 = 1 
         else 
            nl1 = ndexs(sngl(wlmin), wavelen, nlines) + 1 
         end if 
c             <- if ((wavelen(1)) .gt. wlmin) 
 
         if (dble(wavelen(nlines)) .lt. wlmax) then 
            nl2 = nlines 
         else 
            nl2 = ndexs(sngl(wlmax), wavelen, nlines) 
         end if 
c             <- if ((wavelen(nlines)) .lt. wlmax) 
 
         if (nl1 .gt. 1) then 
c Some lines are outside included wavelength range, so shift those which 
c are included over so that they start in the first position. 
            do 45 nnl = nl1, nl2 
               nl = nnl - nl1 + 1 
               elower(nl) = elower(nnl) 
               eupper(nl) = eupper(nnl) 
               gf(nl) = gf(nnl) 
               ionstage(nl) = ionstage(nnl) 
               wavelen(nl) = wavelen(nnl) 
               znucline(nl) = znucline(nnl) 
               element(nl) = element(nnl) 
 45         continue 
 
            nlines = nl2 - nl1 + 1 
         else 
 
            if (nl2 .lt. nlines) then 
               nlines = nl2 - nl1 + 1 
            end if 
c                <- if (nl2 .lt. nlines) 
         end if 
c             <- if (nl1 .gt. 1) 
 
         write(7,'(a,a,i10)') ' lineexpop: lines inside ', 
     .        'included wavelength interval: ', nlines 
 
c labund is a logical flag which is set to false for elements which have 
c zero total column density. 
         call mzalloc('ptrlabund ', ptrlabund, 99  * 8 ) 
 
         do iz = 1, 99 
            labund(iz) = .false. 
 
            do ii = 1, 6 
               labund(iz) = (labund(iz) .or. 
     .              (dasum(nradii, iondenz(1,ii,iz), 1) .gt. 1.d-10)) 
            end do 
         end do 
 
c Throw out lines of elements which have zero column density. 
 
         do nl = 1, nlines 
            iz = znucline(nl) 
            if (labund(iz)) then 
               lteflgx(nl) = 0 
            else 
               lteflgx(nl) = 10 
            end if 
         end do 
 
c Allocate space for linecountr, used to keep statistics on 
c the number of lines in each ionization stage read in. 
         call mzalloc('ptrlcountr ', ptrlcountr, 4950 * 8 ) 
 
         do l = 1, 4950 
            linecountr(l) = 0 
         end do 
 
c Count number of lines per ion which are being used by this routine. 
 
         do nl = 1, nlines 
            if (lteflgx(nl) .lt. 10) then 
               iz = znucline(nl) 
               i = ionstage(nl) 
               linecountr(i+((iz-1)*iz)/2)  = linecountr(i+((iz-1)*iz)/2
     ~)  + 1 
            end if 
         end do 
 
         write(7,*) 
         write(7,'(a)') ' Lines per ionization stage read in and used:' 
     ~
         write(7,'(a)') '    z    i     count' 
         do iz = 1, 99 
            do i = 1, iz 
               if (linecountr(i+((iz-1)*iz)/2)  .gt. 0) then 
                  write(7,'(1x,i4,1x,i4,1x,i9)') iz, i, 
     .                 linecountr(i+((iz-1)*iz)/2)  
               end if 
            end do 
         end do 
 
         write(7,*) 
         call mzdalloc('ptrlcountr ', ptrlcountr, 0) 
 
c     Sort lines according to the value stored in lteflgx. 
 
           
 
         call itablsrt(nlines, lteflgx, indexv, 1) 
 
         call sgthr(nlines, wavelen, dumvec, indexv) 
         call scopy(nlines, dumvec, 1, wavelen, 1) 
 
         call sgthr(nlines, elower, dumvec, indexv) 
         call scopy(nlines, dumvec, 1, elower, 1) 
 
         call igather(nlines, idumvec, znucline, indexv) 
         call icopy(nlines, idumvec, 1, znucline, 1) 
 
         call igather(nlines, idumvec, ionstage, indexv) 
         call icopy(nlines, idumvec, 1, ionstage, 1) 
 
         call igather(nlines, idumvec, levindex, indexv) 
         call icopy(nlines, idumvec, 1, levindex, 1) 
 
         call igather(nlines, idumvec, lteflgx, indexv) 
         call icopy(nlines, idumvec, 1, lteflgx, 1) 
 
         call sgthr(nlines, gf, dumvec, indexv) 
         call scopy(nlines, dumvec, 1, gf, 1) 
 
         call sgthr(nlines, element, dumvec, indexv) 
         call scopy(nlines, dumvec, 1, element, 1) 
 
         do nl = 1, nlines 
c                                   --------- 
            if (lteflgx(nl) .gt. 0) go to 600 
c                                   --------- 
         end do 
 
c   ------------- 
 600     continue 
c   ------------- 
         nlines = nl - 1 
 
c Sort LTE lines by decreasing value of lower energy level. 
 
           
 
         call stblsort(nlines, elower, indexv, 1) 
 
         call sgthr(nlines, wavelen, dumvec, indexv) 
         call scopy(nlines, dumvec, 1, wavelen, 1) 
 
         call sgthr(nlines, elower, dumvec, indexv) 
         call scopy(nlines, dumvec, 1, elower, 1) 
 
         call igather(nlines, idumvec, znucline, indexv) 
         call icopy(nlines, idumvec, 1, znucline, 1) 
 
         call igather(nlines, idumvec, ionstage, indexv) 
         call icopy(nlines, idumvec, 1, ionstage, 1) 
 
         call sgthr(nlines, gf, dumvec, indexv) 
         call scopy(nlines, dumvec, 1, gf, 1) 
 
         call sgthr(nlines, element, dumvec, indexv) 
         call scopy(nlines, dumvec, 1, element, 1) 
 
           
 
c Compute energy levels used in evaluation of the Boltzmann equation. 
c The Boltzmann factors, exp(-E/kt), are computed on a grid. Values 
c for specific energies are determined by interpolation on this grid. 
 
         denergy = enrmax / dfloat(200  - 1) 
         enrvec(1) = 0.d+00 
 
         do i = 1, 200  
            enrvec(i) = dfloat(i-1) * denergy 
         end do 
 
c Determine energy interval index. 
         enrindex(1) = max(ndexd(dble(elower(1)), enrvec, 200 ), 1) 
 
         do 650 nl = 2, nlines 
 
            if (elower(nl) .eq. elower(nl-1)) then 
               enrindex(nl) = enrindex(nl-1) 
            else 
               enrindex(nl) = max(ndexd(dble(elower(nl)), enrvec, 
     .              200 ), 1) 
            end if 
 
            enrwt(nl) = (elower(nl) - enrvec(enrindex(nl))) / denergy 
 650     continue 
 
      end if 
c     -------------------------> if (init .eq. 0) 
 
        
 
c Compute exponential of energy Boltzmann factors. 
 
      do 420 i = 1, 200  
 
         do nr = nr1, nr2 
            expvecx(nr+(i-1)*nradii)  = dexp(-enrvec(i) * hckt(nr)) 
         end do 
 
 420  continue 
 
c     Compute frequency index of line. 
 
      if (init .eq. 0) then 
 
         do 660 nl = 1, nlines 
            wlindex(nl) = ndexd(dble(wavelen(nl)), wlgrd, nwave) 
 660     continue 
 
      end if 
 
      ndtfac = 128 
      dtfac = 1.d+00 / float(ndtfac) 
 
        
 
c Add on contributions from lte lines to scattering opacity. 
 
      do 920 nl = 1, nlines 
         ne = enrindex(nl) 
         nf = wlindex(nl) 
         scatx = cross0 * (gf(nl) * wavelen(nl)) 
         iz = znucline(nl) 
         io = ionstage(nl) 

         do 910 nr = nr1, nr2 
            expterm = enrwt(nl) * expvecx(nr+(ne+1-1)*nradii)  + 
     .           (1.d+00 - enrwt(nl)) * expvecx(nr+(ne-1)*nradii)  
 
               qfactor = eden(nr) / (eden(nr) + quench_fac(min(io,2)) *
     .                   hfreq3c2(nf))
            tausob = scatx * expterm * idenz(nr+(io-1)*nradii+(iz-1)*31*
     ~nradii)  / dvdravg(nr) 
 
c The following expression approximates exp(-tausob) as 
c 1. / (1. + tausob * dtfac)**ndtfac, with dtfac=1/float(ndtfac). 
cccc DEC-alpha failed 
c            opx = dvdravg(nr) * (1.d+00 - 1.d+00 / (1.d+00 + 
c     .           tausob * dtfac)**ndtfac) / (dlnfreq(nf) * c) 
            opx = dvdravg(nr) * (1.d+00 - 
     .            (1.d+00 / (1.d+00 + tausob * dtfac))**ndtfac) 
     .            / (dlnfreq(nf) * c) 
            lscat(nr+(nf-1)*nradii)  = lscat(nr+(nf-1)*nradii)  + 
     .                                 (1.d+00 - qfactor) * opx 
            lopac(nr+(nf-1)*nradii)  = lopac(nr+(nf-1)*nradii)  +
     .                                 qfactor * opx 
 910     continue 
 
 920  continue 
 
        
 
      if (init .eq. 0) then 
         call mzdalloc('pdumvec ',pdumvec, 0) 
         call mzdalloc('pidmvec ',pidmvec, 0) 
         call mzdalloc('pelement ',pelement, 0) 
         call mzdalloc('pteup ', pteup, 0) 
         call mzdalloc('ptindexv ', ptindexv, 0) 
         call mzdalloc('ptlteflg ', ptlteflg, 0) 
         call mzdalloc('ptrlabund ', ptrlabund, 0) 
      end if 
 
      call mzdalloc('ptridenz ', ptridenz, 0) 
 
      do nr = nr1, nr2 
         do nw = 1, nwave 
c      	    print *,'s>',nw,scatop(nw,nr), lscat(nr+(nw-1)*nradii)
            scatop(nw,nr) = scatop(nw,nr) + lscat(nr+(nw-1)*nradii)  
c            print *,'a>',nw,opacity(nw,nr), lopac(nr+(nw-1)*nradii)
            opacity(nw,nr) = opacity(nw,nr) + lopac(nr+(nw-1)*nradii)
         end do 
      end do 
 
      call mzdalloc('ptrlscat ', ptrlscat, 0) 
      call mzdalloc('ptrlopac ', ptrlopac, 0) 
 

 
      ptrinit = 1 
      init = 1 
 
      return 
 
      end 
      subroutine readbini(iunit, x, lenx) 
      integer x(lenx) 
      read(iunit) x 
      return 
      end 
      subroutine readbinr(iunit, x, lenx) 
      real*4 x(lenx) 
      read(iunit) x 
      return 
      end 
      subroutine writbini(iunit, x, lenx) 
      integer x(lenx) 
      write(iunit) x 
      return 
      end 
      subroutine writbinr(iunit, x, lenx) 
      real*4 x(lenx) 
      write(iunit) x 
      return 
      end 
       
