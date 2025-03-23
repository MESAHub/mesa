
subroutine azbar(xmass, aion, zion, ionmax, ymass, abar, zbar)
   implicit none
   save

!..this routine calculates composition variables for an eos routine

!..input:
!..mass fractions     = xmass(1:ionmax)
!..number of nucleons = aion(1:ionmax)
!..charge of nucleus  = zion(1:ionmax)
!..number of isotopes = ionmax
!..
!..output:
!..molar abundances        = ymass(1:ionmax),
!..mean number of nucleons = abar
!..mean nucleon charge     = zbar

!..declare
   integer :: i, ionmax
   double precision :: xmass(ionmax), aion(ionmax), zion(ionmax), ymass(ionmax), abar, zbar, zbarxx, ytot1

   zbarxx = 0.0d0
   ytot1 = 0.0d0
   do i = 1, ionmax
      ymass(i) = xmass(i)/aion(i)
      ytot1 = ytot1 + ymass(i)
      zbarxx = zbarxx + zion(i)*ymass(i)
   end do
   abar = 1.0d0/ytot1
   zbar = zbarxx*abar
end subroutine azbar
