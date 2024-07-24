c -------------------------------------------------------------------------
c               s u b r o u t i n e    g f f c a l c
c -------------------------------------------------------------------------
c This subroutine computes free-free gaunt factors using a fit to the
c calculations of Karzas and Latter (1961), by E. Gronenschild and 
c R. Mewe 1978, A&A Suppl, 32, 283-305.

      subroutine gffcalc(gauntff, t, freq, nfreq, freqdim, iondim)
      implicit real*8 (a-h, o-z)


      PARAMETER(MAXIONM1=(6  - 1))
 

      PARAMETER(NSBINTVL=3*10 )
 

      integer freqdim
      dimension gauntff(freqdim, iondim), freq(freqdim)

      data bc/1.38054d-16/, pi/3.141592653589793/
      data h/6.6262d-27/
c ----------------------------------------------------------------------
       
      bctinv = 1.d+00 / (bc * t)
      sr3pi = dsqrt(3.d+00 / pi)

      do i = 1, iondim
         z2t = dble(i**2) / t
         gamma2 = 1.578d+05 * z2t
         xcoef = 0.5d+00 * (1.d+00 + dsqrt(10.d+00 * gamma2))
         g2log = dlog10(gamma2)
         a = 1.2 * dexp(-((g2log - 1.d+00) / 3.7d+00)**2)
         b = 0.37 * dexp(-((g2log + 1.d+00) * 0.5d+00)**2)

         do nf = 1, nfreq
            u = h * freq(nf) * bctinv
            ulog = dlog10(u)
            gauntff(nf,i) = dsqrt((sr3pi * bessk0ex(xcoef * u))**2
     .           + (a - b * ulog)**2)
         end do

      end do

      return

      end
