subroutine read_inlist_other_winds(ierr)
   integer, intent(out) :: ierr
   character (len = 256) :: filename, message
   integer :: unit
   filename = 'inlist_other_winds'
   
   write(*, *) 'read_inlist_other_winds'
   
   ! set defaults
   wind_scheme = ''
   Kudritzki_scaling_factor = 0d0
   Grafener_scaling_factor = 0d0
   Stern51_scaling_factor = 0d0
   
   open(newunit = unit, file = trim(filename), action = 'read', delim = 'quote', iostat = ierr)
   if (ierr /= 0) then
      write(*, *) 'Failed to open control namelist file ', trim(filename)
   else
      read(unit, nml = other_winds, iostat = ierr)
      close(unit)
      if (ierr /= 0) then
         write(*, *) 'Failed while trying to read control namelist file ', trim(filename)
         write(*, '(a)') &
               'The following runtime error message might help you find the problem'
         write(*, *)
         open(newunit = unit, file = trim(filename), action = 'read', delim = 'quote', status = 'old', iostat = ierr)
         read(unit, nml = other_winds)
         close(unit)
      end if
   end if

end subroutine read_inlist_other_winds


subroutine other_wind(id, L, M, R, Tsurf, X, Y, Z, w, ierr)
   use star_def
   integer, intent(in) :: id
   real(dp), intent(in) :: L, M, R, Tsurf, X, Y, Z ! surface values (cgs)
   ! NOTE: surface is outermost cell. not necessarily at photosphere.
   ! NOTE: don't assume that vars are set at this point.
   ! so if you want values other than those given as args,
   ! you should use values from s% xh(:,:) and s% xa(:,:) only.
   ! rather than things like s% Teff or s% lnT(:) which have not been set yet.
   real(dp), intent(out) :: w ! wind in units of Msun/year (value is >= 0)
   integer, intent(out) :: ierr
   type(star_info), pointer :: s
   w = 0
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return
   if (wind_scheme == 'Stern51') then
      if (s% model_number == 1) write(*, *) 'use Stern51 wind'
      w = eval_stern51_wind(s, L, M, R, Tsurf, X, Y, Z)
   else if (wind_scheme == 'Grafener') then
      Grafener_scaling_factor = 1d0
      if (s% model_number == 1) write(*, *) 'use Grafener wind'
      call eval_Grafener_wind(L, M, R, Tsurf, X, Y, Z, w)
   else if (wind_scheme == 'Kudritzki') then
      Kudritzki_scaling_factor = 1d0
      if (s% model_number == 1) write(*, *) 'use Kudritzki wind'
      call eval_Kudritzki_wind(L, M, R, Tsurf, X, Y, Z, w)
   else
      write(*, *) 'other_wind: invalid wind_scheme: <' // trim(wind_scheme) // '>'
      ierr = -1
   end if
end subroutine other_wind


subroutine eval_Grafener_wind(L, M, R, Tsurf, X, Y, Z, w)
   ! Grafener, G. & Hamann, W.-R. 2008, A&A 482, 945
   ! routine contributed by Nilou Afsari
   real(dp), intent(in) :: L, M, R, Tsurf, X, Y, Z
   real(dp), intent(inout) :: w
   real(dp) :: w1, logMdot, gamma_edd, xsurf, beta, gammazero, lgZ
   real(dp), parameter :: Zsolar = 0.019d0
   include 'formats'
   xsurf = X
   gamma_edd = exp10(-4.813d0) * (1 + xsurf) * (L / Lsun) * (Msun / M)
   lgZ = log10(Z / Zsolar)
   beta = 1.727d0 + 0.250d0 * lgZ
   gammazero = 0.326d0 - 0.301d0 * lgZ - 0.045d0 * lgZ * lgZ
   if (gammazero > gamma_edd) then
      w = 0d0
      return
   end if
   logMdot = &
         + 10.046d0 &
               - 3.5d0 * log10(Tsurf) &
               + beta * log10(gamma_edd - gammazero) &
               + 0.42d0 * log10(L / Lsun) &
               - 0.45d0 * xsurf
   w = exp10(logMdot)
   w = w * Grafener_scaling_factor
end subroutine eval_Grafener_wind


subroutine eval_Kudritzki_wind(L, M, R, Tsurf, X, Y, Z, w)
   real(dp), intent(in) :: L, M, R, Tsurf, X, Y, Z
   real(dp), intent(inout) :: w
   real(dp) :: xlogl, xlteff, stmass, xsurf, xmdfic, vinfkm
   xlogl = log10(L / Lsun)
   xlteff = log10(Tsurf)
   stmass = M / Msun
   xsurf = X
   ! output xmdfic = mass loss rate [Msun/yr]
   ! output vinfkm = v_infinity
   call eval_kuma_wind(xlogl, xlteff, stmass, xsurf, xmdfic, vinfkm)
   w = xmdfic * Kudritzki_scaling_factor
   if (.false.) write(*, *) 'lg eval_Kudritzki_wind', log10(w)
end subroutine eval_Kudritzki_wind


real(dp) function eval_stern51_wind(s, L, M, R, Teff, X, Y, Z) result(dms)
   use chem_def
   type (star_info), pointer :: s
   real(dp), intent(in) :: L, M, R, Teff, X, Y, Z ! cgs
   
   real(dp) :: xhc, xhs, xhs2, znow, fHeavy, ZHeavy, Zlight, Zcno, &
         dms1, dms2, fff, logl, logm, logmd, logr, logz, &
         lteff, vinfkm, xxx, zfac, XC, XN, XO, XNe, XMg, &
         Zbase, ZFe, fC, fN, fO, fNe, fFe
   integer :: i, j, cid, h1, c12, n14, o16, ne20, mg24
   
   include 'formats'
   
   h1 = s% net_iso(ih1)
   Zbase = s% kap_rq% Zbase
   if (Zbase < 0) then
      write(*, *) 'must set Zbase to use stern51 wind'
      stop 'error in eval_stern51_wind'
   end if
   fFe = 0.07d0 ! approx metal fraction by mass of Fe
   ZFe = Zbase * fFe
   Zcno = Z - ZFe ! this is the Stern definition
   
   if (h1 > 0) then
      xhs = s% xa(h1, 1)
      xhc = s% xa(h1, s% nz)
   else
      xhs = X
      xhc = 0
   end if
   xhs2 = X ! surface avg
   
   if (xhs2 > 0.3d0) then ! still have H in envelope
      
      znow = max(0d0, min(1d0, &
            0.02d0 * (0.5d0 * (ZFe / 1.32d-3) + 0.5d0 * (Zcno / 0.01868d0))))
      
      if (xhc < 1d-6) then ! center h gone
         logr = log10(R / Rsun)
         logl = log10(L / Lsun)
         logm = log10(M / Msun)
         logz = log10(znow / 0.02d0)
         logmd = 1.24d0 * logl + 0.16d0 * logm + 0.81d0 * logr - 14.016d0
         dms = exp10(logmd)
      else ! if (xhc >= 1d-6) then ! still have H at center
         logl = log10(L / Lsun)
         lteff = log10(Teff)
         call eval_kuma_wind(logl, lteff, M / Msun, xhs, dms, vinfkm)
      end if
      
      dms = dms * pow(znow / 0.02d0, 0.69d0)
      
      if (xhs2 <= 0.45d0) then ! envelope H is partially depleted
         znow = max(Zbase, 1d-6)
         dms1 = dms
         call wr_hamann(L, xhs, znow, dms2)
         dms2 = dms2 / 10d0
         xxx = (xhs - 0.375d0) * 12.d0
         fff = 0.25d0 * (2d0 - 3d0 * xxx + xxx * xxx * xxx)
         dms = (1.d0 - fff) * dms1 + fff * dms2
         dms = max(dms1, dms)
      end if
   
   else ! if (xhs2 <= 0.3d0) then ! envelope H is gone or going
      
      znow = max(Zbase, 1d-6)
      call wr_hamann(L, xhs, znow, dms)
      dms = dms / 10.d0
      if (Zcno >= znow) then
         zfac = 1.d0 + 19.d0 * (Zcno - znow) / (1.d0 - znow)
         dms = dms * zfac
      end if
   
   end if

end function eval_stern51_wind


subroutine wr_hamann(xsl, xhs, znow, dms)
   real(dp) :: xsl, xhs, znow, dms
   !-------------------------------------------------------------------------
   !     l and x dep. mass loss for WNL and WNE/WC/WO stars
   !     see Hamann et al. 1995, A&A 299, 151
   !     for log(L)<4.5, mass loss is strongly reduced according
   !     to Hamann et al. 1982, A&A 116, 273 (very low rates for
   !     Helium stars with low L)
   !--------------------------------------------------------------------------
   real(dp) :: xlogl, xlogmd, logz
   include 'formats'
   xlogl = log10(xsl / Lsun)
   logz = log10(max(1d-20, znow) / 0.02d0)
   if (xlogl > 4.5d0) then
      xlogmd = 1.5d0 * xlogl - 2.85d0 * xhs - 11.95d0 + 0.86d0 * logz
   else
      xlogmd = 6.8d0 * xlogl - 2.85d0 * xhs - 35.8d0 + 0.86d0 * logz
   end if
   dms = exp10(xlogmd)
   return
   write(*, 1) 'wr_hamann dms', dms, xlogmd
   write(*, 1) 'wr_hamann xlogl', xlogl
   write(*, 1) 'wr_hamann xhs', xhs
   write(*, 1) 'wr_hamann znow', znow
   write(*, 1) 'wr_hamann logz', logz
end subroutine wr_hamann


subroutine eval_kuma_wind(xlogl, xlteff, stmass, xsurf, xmdfic, vinfkm)
   real(dp), intent(in) :: xlogl, xlteff, stmass, xsurf
   real(dp), intent(out) :: xmdfic, vinfkm
   ! input: xlogl= log L/Lsun; xlteff= log Teff, stmass=M/Msun,
   !        xsurf=surface hydrogen mass fraction
   ! output: xmdfic=mass loss rate [Msun/yr], vinfkm=v_infinity
   ! ANALYTICAL APPROXIMATION FOR RADIATION DRIVEN WINDS
   ! PROGRAM WRITTEN BY R.P. KUDRITZKI
   real(dp) :: teff, stlum, ttt, rst, yps, xihe, sig, ga, ganeu, &
         vesckm, vesc, vsound, vsquar, vsoukm, xk, alpha, beta, delta, &
         ucrit, xx, vc, vckms, xmdca, xdbaco, xdjage, xdlame
   
   include 'formats'
   
   teff = exp10(xlteff)
   stlum = exp10(xlogl)
   ttt = teff / 5780.d0
   ttt = ttt * ttt * ttt * ttt
   rst = sqrt(stlum / ttt)
   call mdote(xlogl, rst, stmass, teff, xdjage, xdlame, xdbaco)
   yps = (0.98d0 - xsurf) / (4.d0 * xsurf)
   xihe = 2.d0
   sig = sigmae(yps, xihe)
   ga = gammma(stlum, stmass, sig)
   vesckm = stmass * (1.d0 - ga) / rst
   vesckm = 617.d0 * sqrt(vesckm)
   vesc = vesckm * 1.d05
   vsound = (2.d0 + (1.d0 + xihe) * yps) / (1.d0 + 4.d0 * yps)
   vsound = vsound * teff
   vsound = 9.085d3 * sqrt(vsound)
   vsquar = vsound / vesc
   vsquar = vsquar * vsquar
   vsoukm = vsound * 1.d-5
   
   ! Pauldrach et al. 1993 (MPA 740) for Zeta Pup
   xk = 0.085d0
   alpha = 0.657d0
   beta = 1.0d0
   delta = 0.095d0
   ucrit = uc(alpha, beta, delta, vsquar)
   
   ! Approx. wind integral
   xx = xintap(alpha, beta, delta, ucrit)
   vc = vaucr(vsound, vesc, alpha)
   xx = xx + (1.d0 - alpha) * vc * vc / (alpha * vesc * vesc)
   vckms = vc * 1.d-5
   vinfkm = vesckm * sqrt(xx)
   vinfkm = vinfkm * sqrt(alpha / (1.d0 - alpha))
   
   ! === Approx. Mdot formulae ====
   xmdca = xmdcak(xk, alpha, stlum, stmass, yps, teff)
   xmdfic = xmdfc(xk, alpha, beta, delta, xx, ucrit, rst, vesc, yps, xmdca)

end subroutine eval_kuma_wind


real(dp) function cf(u, alpha, beta)
   real(dp), intent(in) :: u, alpha, beta
   ! cf IS THE FINITE CONE ANGLE CORRECTION FACTOR
   real(dp) :: xl, absxl, eps
   include 'formats'
   xl = (beta + 1.) * u - 1.d0
   xl = xl * u / beta
   absxl = abs(xl)
   eps = 1.d-3
   if (absxl <= eps) then
      cf = 1.d0
      return
   end if
   cf = 1.d0 - pow(1.d0 - xl, 1.d0 + alpha)
   cf = cf / xl
   cf = cf / (1.d0 + alpha)
   !write(*,1) 'cf', cf
end function cf


real(dp) function gdel(u, beta, delta)
   real(dp), intent(in) :: u, beta, delta
   ! COMPUTES FUNCTION G(DELTA)
   real(dp) :: q
   include 'formats'
   q = qq(beta, delta)
   gdel = q * u * u + 1.d0
   !write(*,1) 'gdel', gdel
end function gdel


real(dp) function uc(alpha, beta, delta, vsquar)
   real(dp), intent(in) :: alpha, beta, delta, vsquar
   ! ESTIMATE OF RECIPROCAL CRITICAL POINT
   real(dp) :: sigma, xphi, bbb, q, xc, hilf
   include 'formats'
   sigma = alpha * alpha * (1.d0 - alpha) / vsquar
   xphi = phicr(vsquar)
   xphi = xphi * xphi
   sigma = sigma / xphi
   bbb = betcr(alpha, delta)
   q = qq(bbb, delta)
   xc = -2d0 * q / (sigma * 3d0)
   hilf = 1.46d0 / sigma
   hilf = hilf * (alpha + 1.) * pow(alpha, 0.6d0)
   hilf = pow(hilf, 1d0 / 3d0)
   xc = hilf + xc
   if (xc < 1.03d0) xc = 1.03d0
   uc = 1.d0 / xc
   !write(*,1) 'uc', uc
end function uc


real(dp) function phicr(vsquar)
   real(dp), intent(in) :: vsquar
   ! CALCULATES PAULDRACH PHI FUNCTION AT UCRIT
   real(dp) :: vsq, xl
   include 'formats'
   vsq = sqrt(vsquar)
   xl = log10(vsq)
   xl = 0.36d0 + xl
   xl = 0.3d0 * xl
   xl = pow(vsq, xl)
   phicr = 3.0d0 * xl
   !write(*,1) 'phicr', phicr
end function phicr


real(dp) function vaucr(vsound, vesc, alpha)
   real(dp), intent(in) :: vsound, vesc, alpha
   ! COMPUTES VELOCITY AT CRITICAL POINT ACC. TO PAULDRACH
   real(dp) :: vsq, phi, hilf, h
   include 'formats'
   vsq = vsound * vsound / (vesc * vesc)
   phi = phicr(vsq)
   hilf = sqrt(phi / (1.d0 - alpha))
   h = 1.d0 - pow(1.d0 - alpha, 2.d0 / alpha)
   hilf = hilf / h
   vaucr = vsound * hilf
end function vaucr


real(dp) function betcr(alpha, delta)
   real(dp), intent(in) :: alpha, delta
   ! CALCULATES BETA VALUE FOR UCRIT ACC. TO PAULDRACH
   ! LVERS=2 MODIFIED BETACRIT ACC. TO SECOND VERSION
   include 'formats'
   if (alpha <= 0.7d0) then
      if (delta < 0.03d0) then
         betcr = 2d0
      else if (delta < 0.055d0) then
         betcr = 1d0
      else if (delta < 0.085d0) then
         betcr = 0.7d0
      else if (delta < 0.095d0) then
         betcr = 0.5d0
      else
         betcr = 0.25d0
      end if
   else if (delta >= 0.03d0) then
      if (delta >= 0.095d0) then
         betcr = 0.7d0
      else if (delta <= 0.055d0) then
         betcr = 2d0
      else
         betcr = 1d0
      end if
   else
      betcr = 2d0
   end if
   !write(*,1) 'betcr', betcr
end function betcr


real(dp) function qq(beta, delta)
   real(dp), intent(in) :: beta, delta
   ! CALCULATES Q AS FUNCTION OF BETA AND DELTA
   ! Q=A(BETA)**DELTA-1.
   ! A(BETA) LINEAR BETWEEN A(2)=22.1, A(1)=7.5
   !   A(0.7)=4.0, A(0.5)=2.5, A(0.25)=1.18
   real(dp) :: a, b
   include 'formats'
   b = min(2.5d0, max(0.2d0, beta))
   if (b >= 1.0d0) then
      a = 15d0 * (b - 1.) + 7.5d0
   else if (b >= 0.7d0) then
      a = 11.66667d0 * (b - 0.7d0) + 4.0d0
   else if (b >= 0.5d0) then
      a = 7.5d0 * (b - 0.5d0) + 2.5d0
   else
      a = 5.28d0 * (b - 0.25d0) + 1.18d0
   end if
   qq = pow(a, delta - 1)
   !write(*,1) 'qq', qq
end function qq


real(dp) function xintap(alpha, beta, delta, ucrit)
   real(dp), intent(in) :: alpha, beta, delta, ucrit
   ! APPROX> WIND INTEGRAL = EQ. (60) +
   real(dp) :: q, bet1, bet2, bet3, bet4, bet5, &
         uc1, uc2, uc3, uc4, hilf, x, g, z, a0, a1, b1, b2
   include 'formats'
   q = qq(beta, delta)
   call coeff(alpha, beta, a0, a1, b1, b2)
   bet1 = 1.d0 / (1.d0 + beta)
   bet2 = bet1 * bet1
   bet3 = bet2 * bet1
   bet4 = bet3 * bet1
   bet5 = bet4 * bet1
   uc1 = ucrit
   uc2 = uc1 * uc1
   uc3 = uc2 * uc1
   uc4 = uc3 * uc1
   hilf = 1.d0 / (1.d0 - alpha)
   x = a0 * (uc1 - bet1 + hilf * q * (uc3 - bet3) / 3.d0)
   x = x - a1 * ((uc2 - bet2) / 2.d0 + hilf * q * (uc4 - bet4) / 4.d0)
   x = x + bet1
   x = x + b1 * (bet2 / 2.d0 + hilf * q * bet4 / 4.d0)
   x = x - b2 * (bet3 / 3.d0 + hilf * q * bet5 / 5.d0)
   x = x + hilf * q * bet3 / 3.d0
   g = 1.d0 / (a0 - a1 * uc1)
   g = g * pow(1.d0 / (q * uc2 + 1.d0), hilf)
   z = 1.d0 / g
   z = 2.d0 / alpha * (1.d0 - z)
   z = g * (1.d0 + sqrt(z))
   xintap = z * x
   !write(*,1) 'xintap', xintap
end function xintap


subroutine coeff(alpha, beta, a0, a1, b1, b2)
   real(dp), intent(in) :: alpha, beta
   real(dp), intent(out) :: a0, a1, b1, b2
   ! COEFFICIENTS A0 A1 B0 B1 FOR APPROX. TERMINAL VELOCITY +
   real(dp) :: hilfex, alaaf
   include 'formats'
   hilfex = 1.d0 - alpha
   hilfex = 1.d0 / hilfex
   alaaf = pow(1.d0 / (1.d0 + alpha), hilfex)
   a0 = (1.d0 + beta - alaaf) / beta
   a1 = a0 - alaaf
   b1 = alpha * hilfex / (2.d0 * beta)
   b2 = b1 * (beta + 1.d0)
end subroutine coeff


real(dp) function xmdcak(xk, alpha, stlum, stmass, yps, teff)
   real(dp), intent(in) :: xk, alpha, stlum, stmass, yps, teff
   ! MDOT AFTER CAK
   real(dp) :: hilfex, hilf, hiover, xmue, xihe, sig, ga, hi
   include 'formats'
   hilfex = 1.d0 / alpha
   hilf = 1.2762d23 * xk * stlum
   hilf = hilf * 1.d-10
   hilf = pow(hilf, hilfex)
   hiover = pow(1.d10, hilfex)
   xmue = 1.d0
   hilf = hilf * alpha / vtherm(teff, xmue)
   xihe = 2.d0
   sig = sigmae(yps, xihe)
   ga = gammma(stlum, stmass, sig)
   hilfex = hilfex * (1.d0 - alpha)
   hi = sig * (1.d0 - alpha) / (1.6671d27 * stmass * (1.d0 - ga))
   hi = pow(hi, hilfex)
   xmdcak = hilf * hi * hiover
   !write(*,1) 'xmdcak', xmdcak
end function xmdcak


real(dp) function sigmae(yps, xihe)
   real(dp), intent(in) :: yps, xihe
   ! SIGMA = THOMSON ABS. COEFF. DIVIDED BY DENSITY
   ! XIHE  = NUMBER OF e- PROVIDIED BY He NUCLEUS
   ! YPS   = N(He)/N(H)
   sigmae = (1.d0 + xihe * yps) / (1.d0 + 4.d0 * yps)
   sigmae = sigmae * 0.3978d0
end function sigmae


real(dp) function gammma(stlum, stmass, sig)
   real(dp), intent(in) :: stlum, stmass, sig
   ! EDDINGTON LIMIT GAMMA, SIG COMPUTED BY SIGMAE
   gammma = 7.655d-5 * sig * stlum / stmass
end function gammma


real(dp) function vtherm(t, xmue)
   real(dp), intent(in) :: t, xmue
   ! THERMAL VELOCITY OF ATOM WITH ATOMIC WHEIGHT XMUE
   vtherm = 1.2848d4 * sqrt(t / xmue)
end function vtherm


real(dp) function xmdfc(xk, alpha, beta, delta, xint, ucrit, rst, vesc, yps, xmdcak)
   real(dp), intent(in) :: xk, alpha, beta, delta, xint, ucrit, rst, vesc, yps, xmdcak
   ! MDOT FOR FINITE CONE ANGLE [Msun/yr],  EQ. (61)
   ! XMDCAK = CAK MDOT
   ! XINT   = F.C. VELOCITY INTEGRAL
   real(dp) :: xihe, d, alaaf, hiexa, hiexb, hiexc, a, b, c
   include 'formats'
   xihe = 2d0
   d = 3.92859d1 * (1d0 + xihe * yps) / (1d0 + 4d0 * yps)
   d = d / (rst * rst)
   d = d / vesc
   alaaf = (1.d0 - alpha) / alpha
   d = d * 1d-11 * sqrt(alaaf)
   hiexa = 1.d0 / (alpha - delta)
   hiexb = hiexa * delta
   hiexc = hiexa * alpha
   a = d / sqrt(xint)
   a = pow(a, hiexb)
   b = cf(ucrit, alpha, beta)
   b = b * (1.d0 + qq(beta, delta) * ucrit * ucrit)
   b = pow(b, hiexa)
   c = pow(xmdcak, hiexc)
   ! NORMALISATION TO [Msun/yr] +
   c = c * 1.58661d-26
   xmdfc = a * b * c
end function xmdfc


subroutine mdote(xlogl, rst, stmass, teff, xdjage, xdlame, xdbaco)
   real(dp), intent(in) :: xlogl, rst, stmass, teff
   real(dp), intent(out) :: xdjage, xdlame, xdbaco
   ! COMPUTE EMPIRICAL MASS LOSS RATES ACCORDING TO:              +
   ! 1) DE JAGER ET AL.,                          ---> xdjage     +
   ! 2) Lamers, 1981, Ap.J. 245, 593              ---> xdlame     +
   ! 3) Barlow & Cohen, 1977, Ap.J. 213, 737      ---> xdbaco     +
   real(dp) :: x, y
   real(dp), parameter :: &
         a1 = 6.3168d0, a2 = 0.1104d0, a3 = -0.4311d0, a4 = 3.579d0, a5 = -1.571d0, &
         a6 = -0.0109d0, a7 = -0.2175d0, a8 = -0.8381d0, a9 = -1.2487d0, a10 = 1.5822d0
   include 'formats'
   x = log10(teff) - 4.d0
   y = xlogl - 5.d0
   xdjage = a1 + a2 * x + a3 * x * x + a4 * x * x * x + a5 * y + a6 * y * y + a7 * y * y * y + a8 * x * y + a9 * x * x * y + a10 * x * y * y
   xdjage = -xdjage
   xdlame = 1.42d0 * xlogl + 0.61d0 * log10(rst) - 0.99d0 * log10(stmass) - 12.71d0
   xdbaco = 1.10d0 * xlogl + log10(1.35d-12)
end subroutine mdote
      
