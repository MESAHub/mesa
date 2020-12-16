c ------------------------------------------------------------------------
c                 s u b r o u t i n e    g s h f d x s e c
c ------------------------------------------------------------------------
c This routines returns the Hartree-Fock-Dirac ground state photoionization 
c cross sections computed by D. A. Verner, D. G. Yakovlev, I. M. Band and 
c M. B. Trzhaskovskaya, 1993 (Atomic Data and Nuclear Data Tables preprint).
c Ron Eastman, UCSC, 10 Nov 1993.

      subroutine gshfdxsec(z, ne, n_princ, l_ang, npts, datadir, mode,
     .     x, e_thresh, sigma)
      implicit real*8 (a-h, o-z)

c Calling parameters:
c [input]
c     (int) z is charge on the nucleus (z<30).
c     (int) ne is the total number of electrons.
c     (int) n_princ is the orbital principal quantum number.
c     (int) l_ang is the orbital angular momentum quantum number.
c     (int) npts is the number of points in y and sigma.
c     (chr) datadir is a string containing the path to the directory
c           where the photoionzation fitting parameters can be found.
c     (int) mode specifies how the energy grid is given:
c               mode = 0 => x contains e/e_thresh.
c               mode = 1 => x contains e [eV].
c               mode = 2 => x contains frequency [Hertz].
c               mode = 3 => x contains wavelength [Angstroms],
c           where e_thresh is the threshold crossection.
c     (r*8) x(npts) contains the energy grid at which the cross sections
c           should be computed, as specified by the mode parameter.
c  [output]
c     (r*8) e_thresh is the threshold energy in eV.
c     (r*8) sigma(npts) is the photoionization crossection in cm^2.

      integer z, ne, n_princ, l_ang, npts, mode
      real*8 x(npts), sigma(npts)
      character*(*) datadir

      integer init
      logical fo_search

c Arrays to hold first order fitting parameters.
      parameter (n1s=337)
      parameter(n2s=308)
      parameter (n2p=260)
      parameter (n3s=165)
      parameter (n3p=128)
      parameter (n3d=57)
      parameter (n4s=18)

      integer z_1s(n1s), ne_1s(n1s)
      real*8 eth_1s(n1s), e0_1s(n1s), sig0_1s(n1s)
      real*8 ya_1s(n1s), p_1s(n1s), yw_1s(n1s)

      integer z_2s(n2s), ne_2s(n2s)
      real*8 eth_2s(n2s), e0_2s(n2s), sig0_2s(n2s)
      real*8 ya_2s(n2s), p_2s(n2s), yw_2s(n2s)

      integer z_2p(n2p), ne_2p(n2p)
      real*8 eth_2p(n2p), e0_2p(n2p), sig0_2p(n2p)
      real*8 ya_2p(n2p), p_2p(n2p), yw_2p(n2p)

      integer z_3s(n3s), ne_3s(n3s)
      real*8 eth_3s(n3s), e0_3s(n3s), sig0_3s(n3s)
      real*8 ya_3s(n3s), p_3s(n3s), yw_3s(n3s)

      integer z_3p(n3p), ne_3p(n3p)
      real*8 eth_3p(n3p), e0_3p(n3p), sig0_3p(n3p)
      real*8 ya_3p(n3p), p_3p(n3p), yw_3p(n3p)

      integer z_3d(n3d), ne_3d(n3d)
      real*8 eth_3d(n3d), e0_3d(n3d), sig0_3d(n3d)
      real*8 ya_3d(n3d), p_3d(n3d), yw_3d(n3d)

      integer z_4s(n4s), ne_4s(n4s)
      real*8 eth_4s(n4s), e0_4s(n4s), sig0_4s(n4s)
      real*8 ya_4s(n4s), p_4s(n4s), yw_4s(n4s)

c Save all these quantities!

      save z_1s, ne_1s, eth_1s, e0_1s, sig0_1s
      save ya_1s, p_1s, yw_1s

      save z_2s, ne_2s, eth_2s, e0_2s, sig0_2s
      save ya_2s, p_2s, yw_2s

      save z_2p, ne_2p, eth_2p, e0_2p, sig0_2p
      save ya_2p, p_2p, yw_2p

      save z_3s, ne_3s, eth_3s, e0_3s, sig0_3s
      save ya_3s, p_3s, yw_3s

      save z_3p, ne_3p, eth_3p, e0_3p, sig0_3p
      save ya_3p, p_3p, yw_3p

      save z_3d, ne_3d, eth_3d, e0_3d, sig0_3d
      save ya_3d, p_3d, yw_3d

      save z_4s, ne_4s, eth_4s, e0_4s, sig0_4s
      save ya_4s, p_4s, yw_4s


c Arrays to hold second order fitting parameters.
      parameter (ngroups=20)
      real*8 p_g(ngroups), yw_g(ngroups)
      integer ya_g(ngroups), zc_g(ngroups), alpha_g(ngroups)
      integer z0_g(ngroups), i0_g(ngroups)
      real*8 a_g(6,ngroups), b_g(6,ngroups), c_g(6,ngroups)
      character*4 group_name(ngroups)

      common /vybtso/p_g, yw_g, a_g, b_g, c_g, ya_g, zc_g,
     .     alpha_g, z0_g, i0_g, group_name

      data init/1/
c - - - - - - - - - - - - - - - - - - - - - - - - - - -

      if (init .eq. 1) then
c Read in first order fitting parameters.
         call read_vybt_fo_fits(datadir, '1sfits ', n1s, z_1s,
     .        ne_1s, eth_1s, e0_1s, sig0_1s, ya_1s, p_1s, yw_1s)

         call read_vybt_fo_fits(datadir, '2sfits ', n2s, z_2s,
     .        ne_2s, eth_2s, e0_2s, sig0_2s, ya_2s, p_2s, yw_2s)

         call read_vybt_fo_fits(datadir, '2pfits ', n2p, z_2p,
     .        ne_2p, eth_2p, e0_2p, sig0_2p, ya_2p, p_2p, yw_2p)

         call read_vybt_fo_fits(datadir, '3sfits ', n3s, z_3s,
     .        ne_3s, eth_3s, e0_3s, sig0_3s, ya_3s, p_3s, yw_3s)

         call read_vybt_fo_fits(datadir, '3pfits ', n3p, z_3p,
     .        ne_3p, eth_3p, e0_3p, sig0_3p, ya_3p, p_3p, yw_3p)

         call read_vybt_fo_fits(datadir, '3dfits ', n3d, z_3d,
     .        ne_3d, eth_3d, e0_3d, sig0_3d, ya_3d, p_3d, yw_3d)

         call read_vybt_fo_fits(datadir, '4sfits ', n4s, z_4s,
     .        ne_4s, eth_4s, e0_4s, sig0_4s, ya_4s, p_4s, yw_4s)

c Read in second order fitting parameters.
         call read_vybt_so_fits(datadir)

         init = 0

      end if

c Zero xsecs.
      do k = 1, npts
         sigma(k) = 0.d+00
      end do

      e_thresh = 0.d+00

c Check for invalid or inapplicable calling parameters.
      if (z .gt. 29) return
      if (n_princ .gt. 4) return
      if ((n_princ .eq. 4) .and. (l_ang .gt. 0)) return

c Here I am assuming that shells are filled in order of increasing
c principal quantum number and angular momentum quantum number.
c ne_lim is the minimum number of orbital electrons which must be
c present for there to be any population of the orbital.
      ne_lim = ((n_princ - 1) * n_princ * (2 * n_princ - 1)) / 3
     .     + 2 * l_ang**2
c ne_shell is the number of electrons in this shell.
      ne_shell = min(ne - ne_lim, 2 * (2 * l_ang + 1))

      if (ne_shell .lt. 1) return

      if (l_ang .gt. (n_princ - 1)) then
         write(0,10) n_princ, l_ang
 10      format('gshfdxsec: bad quantum numbers: n = ', i5, ' l = ',
     .        i5, '.')
         stop
      end if

c Look through the appropriate list of first order fitting parameters
c to find one which goes with this ion. If none is found, use second
c order fit. The method of searching is not very efficient, but for
c this routine it doesn't need to be. Presumably, photoionization crossections
c are not going to be computed over and over again.

      if (n_princ .eq. 1) then
c 1s.
         if (fo_search(z, ne, l_ang, n1s, z_1s, ne_1s, eth_1s,
     .        p_1s, e0_1s, sig0_1s, yw_1s, ya_1s, mode, 
     .        npts, x, sigma, e_thresh)) return

c No first order fit found, so determine correct group for second
c order fit.

         if (ne .eq. 2) then
            ng = 1
         else if (ne .gt. 10) then
            ng = 3
         else
            ng = 2
         end if

      else if ((n_princ .eq. 2) .and. (l_ang .eq. 0)) then
c 2s.
         if (fo_search(z, ne, l_ang, n2s, z_2s, ne_2s, eth_2s,
     .        p_2s, e0_2s, sig0_2s, yw_2s, ya_2s, mode, 
     .        npts, x, sigma, e_thresh)) return

c Determine 2nd order group.

         if (ne .lt. 11) then
            ng = 4
         else
            ng = 5
         end if

      else if ((n_princ .eq. 2) .and. (l_ang .eq. 1)) then
c 2p.
         if (fo_search(z, ne, l_ang, n2p, z_2p, ne_2p, eth_2p,
     .        p_2p, e0_2p, sig0_2p, yw_2p, ya_2p, mode, 
     .        npts, x, sigma, e_thresh)) return

c Determine 2nd order group.

         if (ne .lt. 11) then
            if (ne .eq. z) then
               ng = 6
            else
               ng = 7
            end if
         else
            ng = 8
         end if

      else if ((n_princ .eq. 3) .and. (l_ang .eq. 0)) then
c 3s.
         if (fo_search(z, ne, l_ang, n3s, z_3s, ne_3s, eth_3s,
     .        p_3s, e0_3s, sig0_3s, yw_3s, ya_3s, mode, 
     .        npts, x, sigma, e_thresh)) return

c Determine 2nd order group.

         if ((ne .lt. 19) .and. ((z - ne) .lt. 3)) then
            ng = 9
         else if (((z .eq. 23) .and. (ne .eq. 22)) .or.
     .           ((z .eq. 24) .and. ((ne .eq. 23) .or. (ne .eq. 24)))
     .           .or. ((z .eq. 27) .and. (ne .eq. 26)) .or.
     .           ((z .eq. 28) .and. (ne .eq. 27)) .or.
     .           ((z .eq. 29) .and. ((ne .eq. 28) .or. (ne .eq. 29))))
     .           then
            ng = 10
         else if ((ne .gt. 18) .and. ((z - ne) .lt. 3)) then
            ng = 11
         else
            ng = 12
         end if

      else if ((n_princ .eq. 3) .and. (l_ang .eq. 1)) then
c 3p.
         if (fo_search(z, ne, l_ang, n3p, z_3p, ne_3p, eth_3p,
     .        p_3p, e0_3p, sig0_3p, yw_3p, ya_3p, mode, 
     .        npts, x, sigma, e_thresh)) return

c Determine 2nd order group.

         if ((ne .lt. 19) .and. ((z - ne) .lt. 3)) then
            ng = 13
         else if (((z .eq. 23) .and. (ne .eq. 22)) .or.
     .           ((z .eq. 24) .and. ((ne .eq. 23) .or. (ne .eq. 24)))
     .           .or. ((z .eq. 27) .and. (ne .eq. 26)) .or.
     .           ((z .eq. 28) .and. (ne .eq. 27)) .or.
     .           ((z .eq. 29) .and. ((ne .eq. 28) .or. (ne .eq. 29))))
     .           then
            ng = 14
         else if ((ne .gt. 18) .and. ((z - ne) .lt. 3)) then
            ng = 15
         else
            ng = 16
         end if

      else if ((n_princ .eq. 3) .and. (l_ang .eq. 2)) then
c 3d.
         if (fo_search(z, ne, l_ang, n3d, z_3d, ne_3d, eth_3d,
     .        p_3d, e0_3d, sig0_3d, yw_3d, ya_3d, mode, 
     .        npts, x, sigma, e_thresh)) return

c Determine 2nd order group.

         if (((z .eq. 23) .and. (ne .eq. 22)) .or.
     .        ((z .eq. 24) .and. ((ne .eq. 23) .or. (ne .eq. 24)))
     .        .or. ((z .eq. 27) .and. (ne .eq. 26)) .or.
     .        ((z .eq. 28) .and. (ne .eq. 27)) .or.
     .        ((z .eq. 29) .and. ((ne .eq. 28) .or. (ne .eq. 29))))
     .        then
            ng = 17
         else if ((z - ne) .lt. 3) then
            ng = 18
         else
            ng = 19
         end if

      else if ((n_princ .eq. 4) .and. (l_ang .eq. 0)) then
c 4s.
         if (fo_search(z, ne, l_ang, n4s, z_4s, ne_4s, eth_4s,
     .        p_4s, e0_4s, sig0_4s, yw_4s, ya_4s, mode, 
     .        npts, x, sigma, e_thresh)) return

         ng = 20

      end if

c If here, no first order fit was found -> use second order fit.

c `i' is the net charge on the ion.
         i = z - ne
         ir = i - i0_g(ng)
         zr = z - z0_g(ng)
         e_thresh = c_g(1,ng) + c_g(2,ng) * zr + c_g(3,ng) * ir +
     .        c_g(4,ng) * zr**2 + c_g(5,ng) * ir * zr + 
     .        c_g(6,ng) * ir**2

         if ((ng .eq. 1) .or. (ng .eq. 6)) then
            e0 = a_g(1,ng) + a_g(2,ng) * z + a_g(3,ng) * z**2
            sig0 = ne_shell * (1.d+00 - dble(zc_g(ng)) /
     .           dble(z))**(alpha_g(ng)) / (b_g(1,ng) + b_g(2,ng) * z
     .           + b_g(3,ng) * z**2)
         else
            e0 = a_g(1,ng) + a_g(2,ng) * z + a_g(3,ng) * ne + a_g(4,ng) 
     ~*
     .           z**2 + a_g(5,ng) * z * ne + a_g(6,ng) * ne**2
            sig0 = ne_shell * (1.d+00 - dble(zc_g(ng)) /
     .           dble(z))**(alpha_g(ng)) /(b_g(1,ng) + 
     .           b_g(2,ng) * z + b_g(3,ng) * ne + b_g(4,ng) *
     .           z**2 + b_g(5,ng) * z * ne + b_g(6,ng) * ne**2)
         end if

         call gshfd_sigma_calc(l_ang, e_thresh, p_g(ng), e0,
     .        sig0, yw_g(ng), dble(ya_g(ng)), npts, x, sigma, mode)

         return

         end

c ---------------------------------------------------------------
c This subroutine reads in the first order fitting parameters.
      subroutine read_vybt_fo_fits(datadir, fitfile, n, z, ne, eth,
     .     e0, sig0, ya, p, yw)
      implicit real*8 (a-h, o-z)

      character*(*) datadir, fitfile
      integer z(n), ne(n)
      real*8 eth(n), e0(n), sig0(n)
      real*8 ya(n), p(n), yw(n)


      character*160  path
      logical fexist
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      len_datadir = lnblnk(datadir)
      len_fitfile = lnblnk(fitfile)
      path(1:160 ) = ' '

      if (len_datadir .eq. 0) then
         len_path = len_fitfile + 2
         path(1:2) = './'
         path(3:len_fitfile+2) = fitfile(1:len_fitfile)
      else if (datadir(len_datadir:len_datadir) .ne. '/') then
         len_path = len_fitfile + len_datadir + 1
         path(1:len_datadir) = datadir(1:len_datadir)
         path(len_datadir+1:len_datadir+1) = '/'
         path(len_datadir+2:len_datadir+len_fitfile+1) = fitfile(1:len_f
     ~itfile)
      else
         len_path = len_fitfile + len_datadir
         path(1:len_datadir) = datadir(1:len_datadir)
         path(len_datadir+1:len_datadir+len_fitfile) = fitfile(1:len_fit
     ~file)
      end if

      inquire(file=path, exist=fexist)

      if (.not. fexist) then
         write(0,10) path(1:len_path)
 10      format(' gshfdxsec: bad path or file name. file ', a,
     .        ' not found.')
         stop
      end if

      open(file=path, unit=55, status='old')

      read(55,'(///)')
      do k = 1, n
         read(55,20) z(k), ne(k), eth(k), e0(k), sig0(k), 
     .        ya(k), p(k), yw(k)
      end do

 20   format(9x,i3,i3,6e10.4)

      close(unit=55)

      return

      end

c ----------------------------------------------------------
c This subroutine reads in the second order fitting parameters.
      subroutine read_vybt_so_fits(datadir)
      implicit real*8 (a-h, o-z)

      character*(*) datadir
c Arrays to hold second order fitting parameters.
      parameter (ngroups=20)
      real*8 p_g(ngroups), yw_g(ngroups)
      integer ya_g(ngroups), zc_g(ngroups), alpha_g(ngroups)
      integer z0_g(ngroups), i0_g(ngroups)
      real*8 a_g(6,ngroups), b_g(6,ngroups), c_g(6,ngroups)
      character*4 group_name(ngroups)

      common /vybtso/p_g, yw_g, a_g, b_g, c_g, ya_g, zc_g,
     .     alpha_g, z0_g, i0_g, group_name

      character*160  path
      logical fexist

      character*29  sofitfile
      data sofitfile/'second.fit '/
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      len_datadir = lnblnk(datadir)
      len_sofitfile = 29 
      path(1:160 ) = ' '

      if (len_datadir .eq. 0) then
         len_path = len_sofitfile + 2
         path(1:2) = './'
         path(3:len_sofitfile+2) = sofitfile(1:len_sofitfile)
      else if (datadir(len_datadir:len_datadir) .ne. '/') then
         len_path = len_sofitfile + len_datadir + 1
         path(1:len_datadir) = datadir(1:len_datadir)
         path(len_datadir+1:len_datadir+1) = '/'
         path(len_datadir+2:len_datadir+len_sofitfile+1) = sofitfile(1:l
     ~en_sofitfile)
      else
         len_path = len_sofitfile + len_datadir
         path(1:len_datadir) = datadir(1:len_datadir)
         path(len_datadir+1:len_datadir+len_sofitfile) = sofitfile(1:len
     ~_sofitfile)
      end if

      inquire(file=path, exist=fexist)

      if (.not. fexist) then
         write(0,10) path(1:len_path)
 10      format(' gshfdxsec: bad path or file name. file ', a,
     .        ' not found.')
         stop
      end if

      open(file=path, unit=55, status='old')

      read(55,'(/////)')

      do ng = 1, ngroups
         read(55,20) group_name(ng), ya_g(ng), p_g(ng), yw_g(ng),
     .        zc_g(ng), alpha_g(ng), z0_g(ng), i0_g(ng)
 20      format(a4,5x,i3,9x,f5.2,8x,f5.3,7x,i2,11x,i1,10x,i2,5x,i1)
         read(55,30) a_g(1,ng), a_g(2,ng), a_g(3,ng), a_g(4,ng), a_g(5,n
     ~g),
     .        a_g(6,ng)
         read(55,30) b_g(1,ng), b_g(2,ng), b_g(3,ng), b_g(4,ng), b_g(5,n
     ~g),
     .        b_g(6,ng)
         read(55,30) c_g(1,ng), c_g(2,ng), c_g(3,ng), c_g(4,ng), c_g(5,n
     ~g),
     .        c_g(6,ng)
 30      format(6x,e11.4,5(2x,e10.3))

      end do

      close(unit=55)
      return
      end

c ---------------------------------------------------------------------
c This subroutine does a front to back search of all elements and ionization
c stages which cross section fits are given for the relevant n and l.
c If a match is found, it is used to compute the cross section, and the
c threshold energy is returned. Otherwise, zero is returned.

      logical function fo_search(z, ne, l_ang, nfit, zf, nef, eth,
     .     p, e0, sig0, yw, ya, mode, npts, x, sigma, e_thresh)
      implicit real*8 (a-h, o-z)

      integer z, ne, l_ang, nfit, npts, mode
      integer zf(nfit), nef(nfit)
      real*8 e_thresh, eth(nfit), e0(nfit), sig0(nfit)
      real*8 yw(nfit), ya(nfit), p(nfit)

      real*8 x(npts), sigma(npts)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      fo_search = .false.

      do k = 1, nfit
         if ((z .eq. zf(k)) .and. (ne .eq. nef(k))) then
c First order fit found.               
            e_thresh = eth(k)
            fo_search = .true.
            call gshfd_sigma_calc(l_ang, eth(k), p(k), e0(k),
     .           sig0(k), yw(k), ya(k), npts, x, sigma, mode)
            return
         end if
      end do

      return

      end

c -------------------------------------------------------------------
c This routine does the acutal calculation of the photoionization
c crossection.

      subroutine gshfd_sigma_calc(l_ang, eth, p, e0, sig0, yw, ya,
     .        npts, x, sigma, mode)
      implicit real*8 (a-h, o-z)

      integer l_ang, npts, mode
      real*8 eth, e0, sig0, yw, ya, p

      real*8 x(npts), sigma(npts)

      data hev/4.1357d-15/, c18/2.997925d+18/
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      hevc18 = hev * c18
      q = 5.5 + dble(l_ang) - 0.5 * p

      if (mode .eq. 0) then
c x contains e/eth:
         do n = 1, npts
            if (x(n) .ge. 1.d+00) then
               y = x(n) * eth / e0
               sigma(n) = 1.d-18 * sig0 * ((y - 1.d+00)**2 + yw**2)
     .              / (y**q * (1.d+00 + dsqrt(y / ya))**p)
            end if
         end do
      else if (mode .eq. 1) then
c x contains e [eV]:
         do n = 1, npts
            if (x(n) .ge. eth) then
               y = x(n) / e0
               sigma(n) = 1.d-18 * sig0 * ((y - 1.d+00)**2 + yw**2)
     .              / (y**q * (1.d+00 + dsqrt(y / ya))**p)
            end if
         end do
      else if (mode .eq. 2) then
c x contains frequency:
         do n = 1, npts
            xx = x(n) * hev
            if (xx .ge. eth) then
               y = xx / e0
               sigma(n) = 1.d-18 * sig0 * ((y - 1.d+00)**2 + yw**2)
     .              / (y**q * (1.d+00 + dsqrt(y / ya))**p)
            end if
         end do
      else if (mode .eq. 3) then
c x contains wavelength in angstroms:
         do n = 1, npts
            xx = hevc18 / x(n)
            if (xx .ge. eth) then
               y = xx / e0
               sigma(n) = 1.d-18 * sig0 * ((y - 1.d+00)**2 + yw**2)
     .              / (y**q * (1.d+00 + dsqrt(y / ya))**p)
            end if
         end do
      else
         write(0,10) mode
 10      format(' gshfdxsec: ', i6, ' is and invalid value for',
     .        ' mode parameter.')
         stop
      end if

      return

      end
