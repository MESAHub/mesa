

      subroutine azbar(xmass,aion,zion,ionmax,
     1                 ymass,abar,zbar)
      implicit none
      save

c..this routine calculates composition variables for an eos routine

c..input:
c..mass fractions     = xmass(1:ionmax)
c..number of nucleons = aion(1:ionmax)
c..charge of nucleus  = zion(1:ionmax)
c..number of isotopes = ionmax
c..
c..output:
c..molar abundances        = ymass(1:ionmax), 
c..mean number of nucleons = abar
c..mean nucleon charge     = zbar

c..declare
      integer          i,ionmax
      double precision xmass(ionmax),aion(ionmax),zion(ionmax),
     1                 ymass(ionmax),abar,zbar,zbarxx,ytot1

      zbarxx  = 0.0d0
      ytot1   = 0.0d0
      do i=1,ionmax
       ymass(i) = xmass(i)/aion(i)
       ytot1    = ytot1 + ymass(i)
       zbarxx   = zbarxx + zion(i) * ymass(i)
      enddo
      abar   = 1.0d0/ytot1
      zbar   = zbarxx * abar
      return
      end
