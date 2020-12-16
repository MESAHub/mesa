      real*8 cr, compton, avogadro, pi, c_e, cd, clight,
     &  electron_mass, h_mass, ct, cpe, c2, alpha20, alpha2,
     &  boltzmann, echarge, ergsperev, stefan_boltzmann, prad_const,
     &  planck, rydberg, bohr, ergspercmm1, proton_mass
      parameter (pi = 3.141592653589793d0)
      parameter (clight = 2.99792458d10)    !c
      parameter (planck = 6.6260693d-27)    !h
      parameter (cr = 8.314472d7)      !R (gas constant)
      parameter (avogadro = 6.0221415d23)    !Avogadro number N_0
      parameter (rydberg = 109737.31568525d0)  !R(infinity)
      parameter (proton_mass = 1.00727646688d0) ! units are AMU
      parameter (echarge = 1.60217653d-19*1.d-1*clight)
      parameter (ergsperev = (echarge/(1.d-1*clight))*1.d7)
      parameter (electron_mass = avogadro*rydberg*
     &    (planck/(echarge*echarge))*
     &    (planck/echarge)*(planck/echarge)*clight/(2.d0*pi*pi))
      parameter (h_mass = proton_mass + electron_mass)
      parameter (boltzmann = cr/avogadro)
      parameter (compton = planck*avogadro/(electron_mass*clight))
      parameter (c_e = 8.d0*pi/(compton*compton*compton))  !n_e = c_e*rhostar
      parameter (cd = c_e/avogadro)  !8 pi H/lambda_c^3
      parameter (ct = cr/(electron_mass*clight*clight))  !k/m c^2
      parameter (cpe = cr*cd/ct)  !8 pi m c^2/lambda_c^3
      parameter (c2 = planck*clight/boltzmann)
      parameter (alpha20 = c2*c2*cr/(2.d0*pi*clight))
      parameter (alpha2 = alpha20*alpha20*alpha20*
     &  (avogadro/(clight*clight))*(avogadro/clight))
      parameter (stefan_boltzmann = 2.d0*pi*pi*pi*pi*pi/
     &  (15.d0*c2*c2*c2)*clight*boltzmann)
      parameter (prad_const = 4.d0*stefan_boltzmann/(3.d0*clight))
      parameter (bohr = 0.5d0*echarge*echarge/(rydberg*clight*planck))
      parameter (ergspercmm1 = c2*boltzmann)
