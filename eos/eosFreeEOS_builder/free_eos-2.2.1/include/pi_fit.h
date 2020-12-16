      real*8 xdh10, xmocp10
      data xdh10, xmocp10/-0.4d0, 0.d0/
      integer*4 nelements_pi_fit
      parameter (nelements_pi_fit = 20)
      real*8 pi_fit_neutral_ln_original(nelements_pi_fit+2),
     &  pi_fit_neutral_ln(nelements_pi_fit+2),
     &  pi_fit_neutral_ln_saumon(nelements_pi_fit+2)
      data pi_fit_neutral_ln_original/0.d0, 0.d0, 18*0.d0, 2*0.d0/
      data pi_fit_neutral_ln/-0.4d0, -0.2d0, 18*0.d0, -4.0d0, -40.d0/
      data pi_fit_neutral_ln_saumon/-0.4d0, -0.2d0, 18*0.d0,
     &  -4.0d0, -40.d0/
      integer*4 nions_pi_fit
      parameter (nions_pi_fit = 295)
      real*8 pi_fit_ion_ln_original(nions_pi_fit+2),
     &  pi_fit_ion_ln(nions_pi_fit+2),
     &  pi_fit_ion_ln_saumon(nions_pi_fit+2)
      data pi_fit_ion_ln_original/0.d0, 0.d0, 0.d0, 292*0.d0, 0.d0,
     &  0.d0/
      data pi_fit_ion_ln/
     &  -1.8d0, !H
     &  -1.8d0, -1.8d0, !He
     &  0.d0, 0.d0,  4*6.d0, !C
     &  0.d0, 0.d0,  5*6.d0, !N
     &  0.d0, 0.d0,  6*6.d0, !O
     &  0.d0, 0.d0,  8*6.d0, !Ne
     &  0.d0, 0.d0,  9*6.d0, !Na
     &  0.d0, 0.d0, 10*6.d0, !Mg
     &  0.d0, 0.d0, 11*6.d0, !Al
     &  0.d0, 0.d0, 12*6.d0, !Si
     &  0.d0, 0.d0, 13*6.d0, !P
     &  0.d0, 0.d0, 14*6.d0, !S
     &  0.d0, 0.d0, 15*6.d0, !Cl
     &  0.d0, 0.d0, 16*6.d0, !A
     &  0.d0, 0.d0, 18*6.d0, !Ca
     &  0.d0, 0.d0, 20*6.d0, !Ti
     &  0.d0, 0.d0, 22*6.d0, !Cr
     &  0.d0, 0.d0, 23*6.d0, !Mn
     &  0.d0, 0.d0, 24*6.d0, !Fe
     &  0.d0, 0.d0, 26*6.d0, !Ni
     &  -1.8d0, -1.8d0/ !H2, H2+
      data pi_fit_ion_ln_saumon/-40.d0, -40.d0, -40.d0, 294*-40.d0/

      real*8 quad_original, quad_modified
      parameter(quad_original = 10.d0)
      parameter(quad_modified = 10.d0)

      real*8 pi_fitx_neutral_ln_original, pi_fitx_ion_ln_original
      real*8 pi_fitx_neutral_ln, pi_fitx_ion_ln
      real*8 pi_fitx_neutral_ln_saumon, pi_fitx_ion_ln_saumon
      data pi_fitx_neutral_ln_original, pi_fitx_ion_ln_original
     &  /2*0.d0/
      data pi_fitx_neutral_ln, pi_fitx_ion_ln/0.d0,-1.8d0/
      data pi_fitx_neutral_ln_saumon, pi_fitx_ion_ln_saumon
     &  /0.d0,-40.d0/
