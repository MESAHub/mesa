c This is a front end routine for computing absorptive opacity
c and scattering opacity.

      subroutine opacity(rho, temp, abund, nfreq, freq,
     .     dvdr, linedatafile, rdlndump, mkshortlist, longlist,
     .     taumin, xsecdatadir, opac, scatop, eden)
      implicit real*8 (a-h, o-z)

c Calling parameters:
c <input>
c     rho = gas density in g/cc.
c     temp = temperature.
c     abund(99) = mass fractions through z=99.
c     nfreq = number of frequency grid points.
c     freq(nfreq) = frequency grid (freq(nf+1)>freq(nf)).
c     dvdr = velocity gradient (in inverse seconds).
c     linedatafile = name of file holding linelist.
c     rdlndump = logical flag which should be set to true if linedatafile
c      is a binary file.
c     mkshortlist = logical flag, set to true if you want routine
c      lineexpop to compute a short list of lines.
c     longlist = name of file to containing long list of lines.
c     taumin = minimum Sobolev optical depth a lines must have to
c      be included in the shortlist.
c     xsecdatadir = path to directory which contains the Hartree-Fock-Dirac
c      photoionization data.
c <output>
c     opac(nfreq) = array used to store computed opacity.
c     scatop(nfreq) = array used to store computed scattering opacity.
c     eden = computed electron density.

      integer nfreq
      character*(*) linedatafile, longlist, xsecdatadir
      logical rdlndump, mkshortlist
      real*8 rho, temp, abund(99), freq(nfreq), taumin
      real*8 opac(nfreq), scatop(nfreq)

      real*8 iondenz(6,99), iondenz_part(6,99)
      integer ionstg1(99)
      real*8 frac(31), part(31), nucden


C      pointer (ptrexpnu, expnu)
      parameter(iptrexpnu=9999)
      real*8 expnu(iptrexpnu)
C      pointer (ptrplnkfnc, plnkfnc)
      parameter(iptrplnkfnc=9999)
      real*8 plnkfnc(iptrplnkfnc)
C      pointer (ptrwave, wave)
      parameter(iptrwave=9999)
      real*8 wave(iptrwave)
C      pointer (ptrhfreq3c2, hfreq3c2)
      parameter(iptrhfreq3c2=9999)
      real*8 hfreq3c2(iptrhfreq3c2)
C      pointer (ptrdlnfreq, dlnfreq)
      parameter(iptrdlnfreq=9999)
      real*8 dlnfreq(iptrdlnfreq)
C      pointer (ptrfreqm3, freqm3)
      parameter(iptrfreqm3=9999)
      real*8 freqm3(iptrfreqm3)
C      pointer (ptrsigma, sigma)
      parameter(iptrsigma=9999)
      real*8 sigma(iptrsigma)

C      pointer (ptrgff, gff)
      parameter(iptrgff=9999*30)
      real*8 gff(iptrgff)

      common/Uex/SUex

      save ptrplnkfnc, ptrexpnu, ptrwave, ptrhfreq3c2
      save ptrdlnfreq, ptrgff, ptrfreqm3, ptrsigma
      save init, lopinit


      data pi/3.141592653589793d+00/, fourpi/12.5637061d+00/
      data bc/1.380626d-16/, h/6.626205d-27/, hev/4.1357d-15/
      data c/2.997925d+10/, elecxsec/6.6524d-25/
      data evtoerg/1.6022d-12/, a/7.56464d-15/
      data stefbltz/5.66956d-05/, hmass/1.67352d-24/
      data esu/4.80298d-10/, emass/9.1091d-28/
      data srpi/1.772453851d+00/, bcev/8.617064d-05/

      data init/1/, ffop0/3.6914403278312d+08/,  mode/2/
      data rydnu/3.2880513d+15/, n_hyd_max/9/, lopinit/0/
c - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      if (init .eq. 1) then
c Get space for holding auxiliary frequency grid quantities and
c compute them.
         call mzalloc('ptrexpnu ', ptrexpnu, nfreq * 8 ,
     ~ iptrexpnu*8,'iptrexpnu     ','opacity  ')
         call mzalloc('ptrplnkfnc ', ptrplnkfnc, nfreq * 8 ,
     ~ iptrplnkfnc*8,'iptrplnkfnc     ','opacity  ')
         call mzalloc('ptrwave ', ptrwave, nfreq * 8 ,
     ~ iptrwave*8,'iptrwave     ','opacity  ')
         call mzalloc('ptrhfreq3c2 ', ptrhfreq3c2, nfreq * 8 ,
     ~ iptrhfreq3c2*8,'iptrhfreq3c2     ','opacity  ')
         call mzalloc('ptrdlnfreq ', ptrdlnfreq, nfreq * 8 ,
     ~ iptrdlnfreq*8,'iptrdlnfreq     ','opacity  ')
         call mzalloc('ptrfreqm3 ', ptrfreqm3, nfreq * 8 ,
     ~ iptrfreqm3*8,'iptrfreqm3     ','opacity  ')
         call mzalloc('ptrsigma ', ptrsigma, nfreq * 8 ,
     ~ iptrsigma*8,'iptrsigma     ','opacity  ')
         call mzalloc('ptrgff ', ptrgff, nfreq * 8  * 30,
     ~ iptrgff*8,'iptrgff     ','opacity  ')

c Compute rest of stuff needed for frequency grid.

         do nf = 1, nfreq
            wave(nf) = c * 1.d+08 / freq(nf)
            hfreq3c2(nf) = h * freq(nf)**3 / c**2
            freqm3(nf) = freq(nf)**(-3)
         end do

         do nf = 1, nfreq - 1
            dlnfreq(nf) = dlog(freq(nf+1) / freq(nf))
         end do

         dlnfreq(nfreq) = 1.d+00
         init = 0
      end if

      bct = bc * temp
      twopi = 2.d+00 * pi

c Compute Planck function.
      do nf = 1, nfreq
         expnu(nf) = dexp(-h * freq(nf) / bct)
         plnkfnc(nf) = hfreq3c2(nf) * expnu(nf) * (1.d+00 - expnu(nf))
      end do

c total nuclear density.
      nucden = abund(1) * rho / hmass

      do iz = 2, 99
         nucden = nucden + abund(iz) * rho / (dfloat(2 * iz) * hmass)
      end do

c Compute the electron density.
      maxiter = 30
      acc = 1.d-03
      numelems = 99
      call edensol(rho, temp, abund, eden, maxiter, acc,
     .     numelems)

      if (eden .ne. eden) print *, ' Bardak:',abund
c ****************************************************
      write(0,3) rho, temp, eden / nucden
 3    format(' rho: ', 1pe10.3, ' temp: ', 1pe10.3,
     .     ' eden / nucden: ', 1pe10.3)
c ****************************************************

c the pressure
      press = bct * (eden + nucden)

c compute ionization fractions.

      SUex=0.d0
      do iz = 1, 99
         if (iz .eq. 1) then
            totden = abund(iz) * rho / hmass
         else
            totden = abund(iz) * rho / (dfloat(2 * iz) * hmass)
         end if

c Solve Saha-Boltzmann equation to compute the ionization.
         call sahaeqn(1, .true., press, temp, eden, iz,
     .        6, ionstg1(iz), frac, part)

c    if(iz.lt.4)
c     +  write(0,'(a,1p6g11.4)')' frac=',(frac(iqq),iqq=1,6)
         do i = 1, min(6, iz+1)
            iondenz(i,iz) = totden * frac(i)
        iondenz_part(i,iz) = iondenz(i,iz) / part(i)
           SUex=SUex+(abund(iz)/dfloat(2*iz))*log(part(i))
         end do
c    if(iz.lt.4)
c     +  write(0,'(a,1p6g11.4)')' iondenz=',(iondenz(iqq,iz),iqq=1,6)
      end do

c zero things.
      do nf = 1, nfreq
         scatop(nf) = eden * elecxsec
         opac(nf) = 0.d+00
      end do

c Compute the line scattering opacity.


      call lineexpop(linedatafile, rdlndump, nfreq, nfreq,
     .     hfreq3c2, dlnfreq, wave, scatop, opac, dummy, plnkfnc,
     .     lopinit, mkshortlist, taumin, longlist, 1, 1,
     .     eden, temp, dvdr, iondenz_part, ionstg1, 1, 1)


c Compute free-free opacity.

      call gffcalc(gff, temp, freq, nfreq, nfreq, 30)

      ffgcoeff = ffop0 * eden / dsqrt(temp)

      do iz = 1, 99
         if (abund(iz) .gt. 0.d+00) then
            i1 = ionstg1(iz)
c            write(7,'(a,2i5)')' ionstg1 iz ',i1,iz

            do k = 1, min(6, iz+1)
c seb changed:
               i = k + ionstg1(iz) - 1
               ieff = i-1
               z2 = dfloat(ieff)**2
               if (ieff .gt. 0) then
                 do nf = 1, nfreq
                    opterm = ffgcoeff * freqm3(nf)
                    ffop = iondenz(k,iz) * opterm * z2 *
     *                     gff(nf+(ieff-1)*nfreq)
                    opac(nf) = opac(nf) + ffop *(1.d+00 -  expnu(nf))
                 end do
c                      <- end sum over freq, nf
               endif
            end do
c                <- sum over abundant ions.
         end if
c             <- if (abund(iz) .gt. 0)
      end do
c          <- sum over elements.

c Compute bound-free opacity.

      do iz = 1, 99
         if (abund(iz) .gt. 0.d+00) then
            i1 = ionstg1(iz)

            do k = 1, min(6, iz)
               i = ionstg1(iz) + k - 1
               ne = iz - i + 1

               if (ne .gt. 1) then
                  call valence_nl(iz, i, n_princ, l_ang)
                  if ((n_princ .lt. 4) .or. ((n_princ .eq. 4) .and.
     .                 (l_ang .eq. 0))) then
                     call gshfdxsec(iz, ne, n_princ, l_ang, nfreq, xsecd
     ~atadir,
     .                    mode, freq, e_thresh, sigma)

                     do nf = 1, nfreq
                        opac(nf) = opac(nf) + iondenz_part(k,iz)
     .                       * sigma(nf) * dfloat(2 * l_ang + 1) *
     .                       (1.d+00 - expnu(nf))
                     end do
                  end if
c                      <- if ((n_princ .lt. 4) .or. ((n_princ .eq. 4)...
               else if (ne .eq. 1) then
                  freq0 = rydnu * dfloat(iz**2)

                  do n_princ = 1, n_hyd_max
                     freqn = rydnu * dfloat(iz**2) / dfloat(n_princ)**2

                     if (freqn .lt. freq(nfreq)) then
                        expfac = 2.d+00 * dfloat(n_princ**2) *
     .                       dexp(-h * (freq0 - freqn) / bct)

                        do nf = 1, nfreq
                           phot = freq(nf) / freqn
                           if (phot .ge. 1.d+00) then
                              eps = dsqrt(phot - 1.d+00) + 1.d-06
                              sig = 6.3d-18 * (dexp(4. * (1. - atan(eps)
     ~ / eps))
     .                             / (1. - dexp(-twopi / eps))) /
     .                             (phot**4 * dfloat(iz**2))
c                              opac(nf) = opac(nf) + iondenz(i,iz) *
c seb changed: iondenz(i,iz) to iondenz(k,iz)
c !!! check also iondenz_part(i,iz) on this subj
c           if ( i .ne. k ) then
c               print *, ' in opacity.f  i,k',i,k
c               pause
c           endif
                              opac(nf) = opac(nf) + iondenz_part(k,iz) *
     .                             expfac * sig * (1.d+00 - expnu(nf))
                           end if
                        end do
c                            <- nf
                     end if
c                         <- if (freqn .lt. freq(nfreq))
                  end do
c                      <- n_princ
               end if
c                   <- if (ne .gt. 1)
            end do
c                <- ion stages, k.
         end if
c             <- if (abund(iz) .gt. 0)
      end do
c          <- iz.

      return

      end
