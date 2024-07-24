c -----------------------------------------------------------------------
c                 f u n c t i o n   h y d x s e c l
c -----------------------------------------------------------------------
c This version of HYDXSECL calls 
      function hydxsecl(phot, z, n, l)
      implicit real*8 (a-h,o-z)
      real*8 phot, z
      integer n, l
      logical hyphinit

      data hyphinit/.true./

       

      if (hyphinit) then
         call hypho(-1.d+00, 0, 0, 0, 0, 0.0d+00, xsec)
         hyphinit = .false.
      end if

      call hypho(z, n, l, l, 1, phot / dble(n**2), xsec)
      hydxsecl = xsec
      return
      end
