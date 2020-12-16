c This function is the routine bessk0 from `Numerical Recipes', modified
c so that it returns the expression exp(x) * K0(x). - R. Eastman 11/15/93

      function bessk0ex(x)
      implicit real*8 (a-h,o-z)
      real*8 y,p1,p2,p3,p4,p5,p6,p7,
     .     q1,q2,q3,q4,q5,q6,q7
      data p1,p2,p3,p4,p5,p6,p7/-0.57721566d0,0.42278420d0,0.23069756d0,
     ~
     .     0.3488590d-1,0.262698d-2,0.10750d-3,0.74d-5/
      data q1,q2,q3,q4,q5,q6,q7/1.25331414d0,-0.7832358d-1,0.2189568d-1,
     ~
     .     -0.1062446d-1,0.587872d-2,-0.251540d-2,0.53208d-3/

      if (x.le.2.0) then
         y=x*x/4.0
         bessk0ex=((-log(x/2.0)*bessi0(x))+(p1+y*(p2+y*(p3+
     .        y*(p4+y*(p5+y*(p6+y*p7))))))) * exp(x)
      else
         y=(2.0/x)
         bessk0ex=(q1+y*(q2+y*(q3+
     .        y*(q4+y*(q5+y*(q6+y*q7)))))) / sqrt(x)
      endif

      return

      end
