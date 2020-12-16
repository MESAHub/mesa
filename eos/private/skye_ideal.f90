module skye_ideal
   use math_lib
   use auto_diff
   use const_def

   implicit none

   private

   public :: compute_F_rad, compute_F_ideal_ion, compute_xne, compute_ideal_ele

   real(dp), parameter :: sifac  = planck_h * planck_h * planck_h / (2d0 * pi * amu * sqrt(2d0 * pi * amu))

   logical, parameter :: dbg = .false.

   contains

   type(auto_diff_real_2var_order3) function compute_F_rad(temp, den) result(Frad)
      type(auto_diff_real_2var_order3), intent(in) :: temp, den

      ! F = - P / rho
      Frad = -crad * pow4(temp) / (3d0 * den)

   end function compute_F_rad

   type(auto_diff_real_2var_order3) function compute_F_ideal_ion(temp, den, abar, species, weights, ya) result(F_ideal_ion)
      type(auto_diff_real_2var_order3), intent(in) :: temp, den
      integer, intent(in) :: species
      real(dp), intent(in) :: weights(species), ya(species), abar

      integer :: j
      type(auto_diff_real_2var_order3) :: n, nj, nQ, nQj
      type(auto_diff_real_2var_order3) :: y, yy, z
      type(auto_diff_real_2var_order3) :: kt, xni, etaion

      ! Helper quantity
      kt = kerg * temp
      n = den / (amu * abar) ! Number density of ions
      nQ = pow(kT, 1.5d0) / sifac ! Quantum density for a 1-amu species

      F_ideal_ion = 0d0
      do j=1,species
         nj = ya(j) * n
         nQj = nQ * pow(weights(j), 1.5d0)
         F_ideal_ion = F_ideal_ion + ya(j) * (log(nj / nQj) - 1d0)
      end do
      F_ideal_ion = F_ideal_ion * kt / (amu * abar)

   end function compute_F_ideal_ion

   type(auto_diff_real_2var_order3) function compute_xne(den, ytot1, zbar) result(xne)
      type(auto_diff_real_2var_order3), intent(in) :: den
      real(dp), intent(in) :: ytot1, zbar
      type(auto_diff_real_2var_order3) :: xni

      ! xne is the electron density due to matter, and so is proportional to density and independent of temperature
      ! for a fully ionized system.
      xni     = avo * ytot1 * den
      xne    = xni * zbar

   end function compute_xne

   subroutine compute_ideal_ele(temp_in, den_in, din_in, logtemp_in, logden_in, zbar, ytot1, ye, ht, s, e, p, adr_etaele, adr_xnefer, ierr)
      use helm_polynomials
      use eos_def
      real(dp), intent(in) :: temp_in, den_in, din_in, logtemp_in, logden_in, zbar, ytot1, ye
      type (Helm_Table), pointer, intent(in) :: ht

      ! Intermediates
      real(dp) :: x, y, z, ww, zz, xni
      real(dp) :: temp, den, din, logtemp, logden
      integer :: k

      !..for the interpolations
      integer          iat, jat
      real(dp) dth, dt2, dti, dt2i, dt3i, dd, dd2, ddi, dd2i, dd3i, &
                       xt, xd, mxt, mxd, fi(36), &
                       dindd, dinda, dindz, dindda, dinddz, dindaa, &
                       dindaz, dindzz, dinddaa, dinddaz, &
                       w0t, w1t, w2t, w0mt, w1mt, w2mt, &
                       w0d, w1d, w2d, w0md, w1md, w2md, &
                       dpepdd_in, dpepddd_in, dpepddt_in

      real(dp) si0t, si1t, si2t, si0mt, si1mt, si2mt, &
                       si0d, si1d, si2d, si0md, si1md, si2md, &
                       dsi0t, dsi1t, dsi2t, dsi0mt, dsi1mt, dsi2mt, &
                       dsi0d, dsi1d, dsi2d, dsi0md, dsi1md, dsi2md, &
                       ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, &
                       ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md, &
                       dddsi0t, dddsi1t, dddsi2t, &
                       dddsi0mt, dddsi1mt, dddsi2mt, &
                       dddsi0d, dddsi1d, dddsi2d, &
                       dddsi0md, dddsi1md, dddsi2md

      ! For some results we don't need, but which helm_electron_positron.dek computes
      integer :: elemult

      real(dp) :: detada,detadz
      real(dp) :: detadda,detaddz
      real(dp) :: detadta,detadtz,detadaa,detadaz,detadzz
      real(dp) :: detadddd, detadddt, detaddtt, detadttt

      real(dp) :: dpepda,dpepdz
      real(dp) :: dpepdda,dpepddz
      real(dp) :: dpepdta,dpepdtz,dpepdaa
      real(dp) :: dpepdaz,dpepdzz
      real(dp) :: deepda,deepdz
      real(dp) :: deepdda,deepddz
      real(dp) :: deepdta,deepdtz,deepdaa
      real(dp) :: deepdaz,deepdzz
      real(dp) :: dsepda,dsepdz
      real(dp) :: dsepdda,dsepddz
      real(dp) :: dsepdta,dsepdtz,dsepdaa
      real(dp) :: dsepdaz,dsepdzz

      real(dp) :: etapos,zeff

      real(dp) :: dpeledd,dpeledt,dpeleda,dpeledz
      real(dp) :: dpeleddd,dpeleddt,dpeledda,dpeleddz
      real(dp) :: dpeledtt,dpeledta,dpeledtz,dpeledaa
      real(dp) :: dpeledaz,dpeledzz
      real(dp) :: deeledd,deeledt,deeleda,deeledz
      real(dp) :: deeleddd,deeleddt,deeledda,deeleddz
      real(dp) :: deeledtt,deeledta,deeledtz,deeledaa
      real(dp) :: deeledaz,deeledzz
      real(dp) :: dseledd,dseledt,dseleda,dseledz
      real(dp) :: dseleddd,dseleddt,dseledda,dseleddz
      real(dp) :: dseledtt,dseledta,dseledtz,dseledaa
      real(dp) :: dseledaz,dseledzz

      real(dp) :: ppos,dpposdd,dpposdt,dpposda,dpposdz
      real(dp) :: dpposddd,dpposddt,dpposdda,dpposddz
      real(dp) :: dpposdtt,dpposdta,dpposdtz,dpposdaa
      real(dp) :: dpposdaz,dpposdzz
      real(dp) :: epos,deposdd,deposdt,deposda,deposdz
      real(dp) :: deposddd,deposddt,deposdda,deposddz
      real(dp) :: deposdtt,deposdta,deposdtz,deposdaa
      real(dp) :: deposdaz,deposdzz
      real(dp) :: spos,dsposdd,dsposdt,dsposda,dsposdz
      real(dp) :: dsposddd,dsposddt,dsposdda,dsposddz
      real(dp) :: dsposdtt,dsposdta,dsposdtz,dsposdaa
      real(dp) :: dsposdaz,dsposdzz

      real(dp) :: xne
      real(dp) :: dxnedd,dxnedt,dxneda,dxnedz
      real(dp) :: dxneddd,dxneddt,dxnedda,dxneddz
      real(dp) :: dxnedtt,dxnedta,dxnedtz,dxnedaa
      real(dp) :: dxnedaz,dxnedzz

      real(dp) :: xnefer,dxneferdd,dxneferdt,dxneferda,dxneferdz
      real(dp) :: dxneferddd,dxneferddt,dxneferdda,dxneferddz
      real(dp) :: dxneferdtt,dxneferdta,dxneferdtz,dxneferdaa
      real(dp) :: dxneferdaz,dxneferdzz
      real(dp) :: xnpfer,dxnpferdd,dxnpferdt,dxnpferda,dxnpferdz
      real(dp) :: dxnpferddd,dxnpferddt,dxnpferdda,dxnpferddz
      real(dp) :: dxnpferdtt,dxnpferdta,dxnpferdtz,dxnpferdaa
      real(dp) :: dxnpferdaz,dxnpferdzz
      real(dp) :: dxneferdddd, dxneferdddt, dxneferddtt, dxneferdttt
      real(dp) :: dxneferddda, dxneferdtta, dxneferddta
      real(dp) :: dxneferdddz, dxneferdttz, dxneferddtz
      real(dp) :: detaddda, detadtta, detaddta
      real(dp) :: detadddz, detadttz, detaddtz

      ! Results we do need
      real(dp) :: free, df_d, df_t, df_dd, df_tt, df_dt, &
                       df_ttt, df_dtt, df_ddt, df_ddd
      real(dp) :: etaele,detadd,detadt
      real(dp) :: detaddd,detaddt,detadtt
      type(auto_diff_real_2var_order3), intent(out) :: adr_etaele, adr_xnefer, s, e, p
      real(dp) :: sele, dsepdt, dsepdd, dsepddt, dsepdtt, dsepddd, &
                 eele, deepdt, deepdd, deepddt, deepdtt, deepddd, &
                 pele, dpepdt, dpepdd, dpepddt, dpepdtt, dpepddd

      integer, intent(out) :: ierr


      temp = temp_in
      logtemp = logtemp_in
      den = den_in
      logden = logden_in
      din = din_in

      include 'helm_electron_positron.dek'

      ! Pack quantities into auto_diff types

      ! Electron chemical potential
      adr_etaele = etaele

      adr_etaele%d1val1 = detadt
      adr_etaele%d1val2 = detadd
!      adr_etaele%d1val3 = detada ! d/dabar
!      adr_etaele%d1val4 = detadz ! d/dzbar

      adr_etaele%d2val1 = detadtt
      adr_etaele%d2val2 = detaddd
      adr_etaele%d1val1_d1val2 = detaddt
!      adr_etaele%d1val1_d1val3 = detadta ! d/dabar
!      adr_etaele%d1val1_d1val4 = detadtz ! d/dzbar
!      adr_etaele%d1val2_d1val3 = detadda ! d/dabar
!      adr_etaele%d1val2_d1val4 = detaddz ! d/dzbar

      adr_etaele%d3val1 = detadttt
      adr_etaele%d2val1_d1val2 = detaddtt
      adr_etaele%d1val1_d2val2 = detadddt
      adr_etaele%d3val2 = detadddd

!      adr_etaele%d2val1_d1val3 = detadtta
!      adr_etaele%d2val2_d1val3 = detaddda
!      adr_etaele%d1val1_d1val2_d1val3 = detaddta
!      adr_etaele%d2val1_d1val4 = detadttz
!      adr_etaele%d2val2_d1val4 = detadddz
!      adr_etaele%d1val1_d1val2_d1val4 = detaddtz

      ! Electron density
      adr_xnefer = xnefer

      adr_xnefer%d1val1 = dxnedt
      adr_xnefer%d1val2 = dxnedd
!      adr_xnefer%d1val3 = dxneda ! d/dabar
!      adr_xnefer%d1val4 = dxnedz ! d/dzbar

      adr_xnefer%d2val1 = dxnedtt
      adr_xnefer%d2val2 = dxneddd
      adr_xnefer%d1val1_d1val2 = dxneddt
!      adr_xnefer%d1val1_d1val3 = dxnedta ! d/dabar
!      adr_xnefer%d1val1_d1val4 = dxnedtz ! d/dzbar
!      adr_xnefer%d1val2_d1val3 = dxnedda ! d/dabar
!      adr_xnefer%d1val2_d1val4 = dxneddz ! d/dzbar

      adr_xnefer%d3val1 = dxneferdttt
      adr_xnefer%d2val1_d1val2 = dxneferddtt
      adr_xnefer%d1val1_d2val2 = dxneferdddt
      adr_xnefer%d3val2 = dxneferdddd

!      adr_xnefer%d2val1_d1val3 = dxneferdtta
!      adr_xnefer%d2val2_d1val3 = dxneferddda
!      adr_xnefer%d1val1_d1val2_d1val3 = dxneferddta
!      adr_xnefer%d2val1_d1val4 = dxneferdttz
!      adr_xnefer%d2val2_d1val4 = dxneferdddz
!      adr_xnefer%d1val1_d1val2_d1val4 = dxneferddtz


      ! Entropy
      s%val = sele

      s%d1val1 = dsepdt
      s%d1val2 = dsepdd
!      s%d1val3 = dsepda ! d/dabar
!      s%d1val4 = dsepdz ! d/dzbar

      s%d2val1 = dsepdtt
      s%d2val2 = dsepddd
      s%d1val1_d1val2 = dsepddt

!      s%d1val1_d1val3 = dsepdta ! d/dabar
!      s%d1val1_d1val4 = dsepdtz ! d/dzbar

!      s%d1val2_d1val3 = dsepdda ! d/dabar
!      s%d1val2_d1val4 = dsepddz ! d/dzbar

      ! Energy
      e%val = eele

      e%d1val1 = deepdt
      e%d1val2 = deepdd
!      e%d1val3 = deepda ! d/dabar
!      e%d1val4 = deepdz ! d/dzbar

      e%d2val1 = deepdtt
      e%d2val2 = deepddd
      e%d1val1_d1val2 = deepddt

!      e%d1val1_d1val3 = deepdta ! d/dabar
!      e%d1val1_d1val4 = deepdtz ! d/dzbar

!      e%d1val2_d1val3 = deepdda ! d/dabar
!      e%d1val2_d1val4 = deepddz ! d/dzbar

      ! Pressure
      p%val = pele

      p%d1val1 = dpepdt
      p%d1val2 = dpepdd
!      p%d1val3 = dpepda ! d/dabar
!      p%d1val4 = dpepdz ! d/dzbar

      p%d2val1 = dpepdtt
      p%d2val2 = dpepddd
      p%d1val1_d1val2 = dpepddt

!      p%d1val1_d1val3 = dpepdta ! d/dabar
!      p%d1val1_d1val4 = dpepdtz ! d/dzbar

!      p%d1val2_d1val3 = dpepdda ! d/dabar
!      p%d1val2_d1val4 = dpepddz ! d/dzbar

   end subroutine compute_ideal_ele

end module skye_ideal
