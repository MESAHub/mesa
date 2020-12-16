      logical*4 ifhe1_special
      parameter(ifhe1_special = .true.)
      integer*4 nx, max_nmin_max, max_izhi, nions_excitation
      parameter (nx = 5)
      parameter(max_nmin_max = 10)
      parameter(max_izhi = 29)
      parameter(nions_excitation = 295)
      real*8 c2t, x(nx), x_old_excitation(nx),
     &  qh2, qh2t, qh2t2, qh2plus, qh2plust, qh2plust2,
     &  qstar(max_nmin_max, max_izhi),
     &  qstart(max_nmin_max, max_izhi),
     &  qstarx(nx, max_nmin_max, max_izhi),
     &  qstart2(max_nmin_max, max_izhi),
     &  qstartx(nx, max_nmin_max, max_izhi),
     &  qstarx2(nx, nx, max_nmin_max, max_izhi),
     &  qmhd_he1, qmhd_he1t, qmhd_he1x(nx),
     &  qmhd_he1t2, qmhd_he1tx(nx), qmhd_he1x2(nx,nx),
     &  psum, psumf, psumt, psum_dv(nions_excitation+2),
     &  ssum, ssumf, ssumt, usum,
     &  free_sum, free_sumf, free_sum_dv(nions_excitation+2)
      real*8 tl_old_excitation
      integer ifpi_fit_old_excitation,
     &  ifh2_old_excitation, ifh2plus_old_excitation
      logical ifpl_logical_old_excitation,
     &  ifmhd_logical_old_excitation,
     &  ifapprox_old_excitation, ifdiff_x_excitation
      common /excitation_block/ c2t, x, x_old_excitation,
     &  qh2, qh2t, qh2t2, qh2plus, qh2plust, qh2plust2,
     &  qstar, qstart, qstarx, qstart2, qstartx, qstarx2,
     &  qmhd_he1, qmhd_he1t, qmhd_he1x,
     &  qmhd_he1t2, qmhd_he1tx, qmhd_he1x2,
     &  psum, psumf, psumt, psum_dv,
     &  ssum, ssumf, ssumt, usum,
     &  free_sum, free_sumf, free_sum_dv,
     &  tl_old_excitation, ifpi_fit_old_excitation,
     &  ifh2_old_excitation, ifh2plus_old_excitation,
     &  ifpl_logical_old_excitation,
     &  ifmhd_logical_old_excitation,
     &  ifapprox_old_excitation, ifdiff_x_excitation
