      FUNCTION BESSI0(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 Y,P1,P2,P3,P4,P5,P6,P7,
     *         Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9
      DATA 
     *   P1,P2,P3,P4,P5,P6,P7/1.0D0,3.5156229D0,3.0899424D0,1.2067492D0,
     ~
     *         0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,
     *         0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,
     *         0.2635537D-1,-0.1647633D-1,0.392377D-2/
      IF (ABS(X).LT.3.75) THEN
         Y=(X/3.75)**2
         BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
         AX=ABS(X)
         Y=3.75/AX
         BESSI0=(EXP(AX)/SQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4
     *            +Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
      RETURN
      END
