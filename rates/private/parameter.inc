C *** TYP declaration, Routinen mit IMPLICIT NONE statement
      INTEGER NSP, NRE, MSL, ndim_rl, nnpartdim
C *** dimension global arrays
C *** for NNN=0
C      PARAMETER( NSP=36, NRE=65, MSL=1 )
C *** vitalwrapper for EVOL vital
      PARAMETER( NSP=300, NRE=109, MSL=1 )
C *** for NNN > 0
C      PARAMETER( NSP=1999, NRE=25000, MSL=1 )
! *** reaclib parameters --> Basel:
      PARAMETER(ndim_rl=61721, nnpartdim=5427) 
! *** reaclib parameters --> Jina (without species beyond At):
!      PARAMETER(ndim_rl=74313, nnpartdim=5456) 
      
