c -----------------------------------------------------------------------
c                   f u n c t i o n   h y d x s e c n
c -----------------------------------------------------------------------
c Compute hydrogenic crossection for hydrogenic state of charge Z and
c radial quantum number N.
      function hydxsecn(phot, z, n)
      implicit real*8 (a-h, o-z)

      real*8 phot, z
      integer n

c PHOT is the photon energy in units of the ionization threshold energy.
c Z is the charge on the nucleus after electron removal.
c N is the principal quantum number.

       
      hydxsecn = 0.d+00

      if (phot .lt. 1.d+00) return

      do 10 l = 0, n - 1
        swt = dble(2 * l + 1)
        hydxsecn = hydxsecn + hydxsecl(phot, z, n, l) * swt
  10  continue

      hydxsecn = hydxsecn / dble(n**2)

      return

      end
