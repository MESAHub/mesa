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

   type(auto_diff_real_2var_order3_1var_order2) function compute_F_rad(temp, den) result(Frad)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: temp, den

      ! F = - P / rho
      Frad = -crad * pow4(temp) / (3d0 * den)

   end function compute_F_rad

   type(auto_diff_real_2var_order3_1var_order2) function compute_F_ideal_ion(temp, den, abar, species, weights, ya) result(F_ideal_ion)
      integer, intent(in) :: species
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: temp, den
      real(dp), intent(in) :: weights(species), ya(species), abar

      integer :: j
      type(auto_diff_real_2var_order3_1var_order2) :: n, nj, nQ, nQj
      type(auto_diff_real_2var_order3_1var_order2) :: y, yy, z
      type(auto_diff_real_2var_order3_1var_order2) :: kt, xni, etaion

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

   type(auto_diff_real_2var_order3_1var_order2) function compute_xne(den, ytot1, zbar) result(xne)
      type(auto_diff_real_2var_order3_1var_order2), intent(in) :: den
      real(dp), intent(in) :: zbar, ytot1
      type(auto_diff_real_2var_order3_1var_order2) :: xni

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
      type(auto_diff_real_2var_order3_1var_order2), intent(out) :: adr_etaele, adr_xnefer, s, e, p
      real(dp) :: sele, dsepdt, dsepdd, dsepddt, dsepdtt, dsepddd, &
                 eele, deepdt, deepdd, deepddt, deepdtt, deepddd, &
                 pele, dpepdt, dpepdd, dpepddt, dpepdtt, dpepddd

      integer, intent(out) :: ierr


      temp = temp_in
      logtemp = logtemp_in
      den = den_in
      logden = logden_in
      din = din_in


!..density derivs
        dindd   = ye
        dinda   = -din*ytot1
        dindz   = den*ytot1
        dindda  = -ye*ytot1
        dinddz  = ytot1
        ww      = ytot1*ytot1
        dindaa  = 2.0d0*din*ww
        dindaz  = -den*ww
        dinddaa = 2.0d0*ye*ww
        dinddaz = -ww
        dindzz = 0.0d0 


!..hash locate this temperature and density
        jat = int((logtemp - ht% logtlo)*ht% logtstpi) + 1
        jat = max(1,min(jat,ht% jmax-1))
        iat = int((log10(din) - ht% logdlo)*ht% logdstpi) + 1
        iat = max(1,min(iat,ht% imax-1))


!..access the table locations only once
        fi(1)  = ht% f(iat,jat)
        fi(2)  = ht% f(iat+1,jat)
        fi(3)  = ht% f(iat,jat+1)
        fi(4)  = ht% f(iat+1,jat+1)
        fi(5)  = ht% ft(iat,jat)
        fi(6)  = ht% ft(iat+1,jat)
        fi(7)  = ht% ft(iat,jat+1)
        fi(8)  = ht% ft(iat+1,jat+1)
        fi(9)  = ht% ftt(iat,jat)
        fi(10) = ht% ftt(iat+1,jat)
        fi(11) = ht% ftt(iat,jat+1)
        fi(12) = ht% ftt(iat+1,jat+1)
        fi(13) = ht% fd(iat,jat)
        fi(14) = ht% fd(iat+1,jat)
        fi(15) = ht% fd(iat,jat+1)
        fi(16) = ht% fd(iat+1,jat+1)
        fi(17) = ht% fdd(iat,jat)
        fi(18) = ht% fdd(iat+1,jat)
        fi(19) = ht% fdd(iat,jat+1)
        fi(20) = ht% fdd(iat+1,jat+1)
        fi(21) = ht% fdt(iat,jat)
        fi(22) = ht% fdt(iat+1,jat)
        fi(23) = ht% fdt(iat,jat+1)
        fi(24) = ht% fdt(iat+1,jat+1)
        fi(25) = ht% fddt(iat,jat)
        fi(26) = ht% fddt(iat+1,jat)
        fi(27) = ht% fddt(iat,jat+1)
        fi(28) = ht% fddt(iat+1,jat+1)
        fi(29) = ht% fdtt(iat,jat)
        fi(30) = ht% fdtt(iat+1,jat)
        fi(31) = ht% fdtt(iat,jat+1)
        fi(32) = ht% fdtt(iat+1,jat+1)
        fi(33) = ht% fddtt(iat,jat)
        fi(34) = ht% fddtt(iat+1,jat)
        fi(35) = ht% fddtt(iat,jat+1)
        fi(36) = ht% fddtt(iat+1,jat+1)

!..various differences
        xt  = max( (temp - ht% t(jat)) * ht% dti_sav(jat), 0.0d0)
        xd  = max( (din - ht% d(iat)) * ht% ddi_sav(iat), 0.0d0)
        mxt = 1.0d0 - xt
        mxd = 1.0d0 - xd

!..the six density and six temperature basis functions
        si0t =   psi0(xt)
        si1t =   psi1(xt) * ht% dt_sav(jat)
        si2t =   psi2(xt) * ht% dt2_sav(jat)

        si0mt =  psi0(mxt)
        si1mt = -psi1(mxt) * ht% dt_sav(jat)
        si2mt =  psi2(mxt) * ht% dt2_sav(jat)

        si0d =   psi0(xd)
        si1d =   psi1(xd) * ht% dd_sav(iat)
        si2d =   psi2(xd) * ht% dd2_sav(iat)

        si0md =  psi0(mxd)
        si1md = -psi1(mxd) * ht% dd_sav(iat)
        si2md =  psi2(mxd) * ht% dd2_sav(iat)

!..first derivatives of the weight functions
        dsi0t =   dpsi0(xt) * ht% dti_sav(jat)
        dsi1t =   dpsi1(xt)
        dsi2t =   dpsi2(xt) * ht% dt_sav(jat)

        dsi0mt = -dpsi0(mxt) * ht% dti_sav(jat)
        dsi1mt =  dpsi1(mxt)
        dsi2mt = -dpsi2(mxt) * ht% dt_sav(jat)

        dsi0d =   dpsi0(xd) * ht% ddi_sav(iat)
        dsi1d =   dpsi1(xd)
        dsi2d =   dpsi2(xd) * ht% dd_sav(iat)

        dsi0md = -dpsi0(mxd) * ht% ddi_sav(iat)
        dsi1md =  dpsi1(mxd)
        dsi2md = -dpsi2(mxd) * ht% dd_sav(iat)

!..second derivatives of the weight functions
        ddsi0t =   ddpsi0(xt) * ht% dt2i_sav(jat)
        ddsi1t =   ddpsi1(xt) * ht% dti_sav(jat)
        ddsi2t =   ddpsi2(xt)
 
        ddsi0mt =  ddpsi0(mxt) * ht% dt2i_sav(jat)
        ddsi1mt = -ddpsi1(mxt) * ht% dti_sav(jat)
        ddsi2mt =  ddpsi2(mxt)

        ddsi0d =   ddpsi0(xd) * ht% dd2i_sav(iat)
        ddsi1d =   ddpsi1(xd) * ht% ddi_sav(iat)
        ddsi2d =   ddpsi2(xd)

        ddsi0md =  ddpsi0(mxd) * ht% dd2i_sav(iat)
        ddsi1md = -ddpsi1(mxd) * ht% ddi_sav(iat)
        ddsi2md =  ddpsi2(mxd)

!..third derivatives of the weight functions
        dddsi0t =   dddpsi0(xt) * ht% dt3i_sav(jat)
        dddsi1t =   dddpsi1(xt) * ht% dt2i_sav(jat)
        dddsi2t =   dddpsi2(xt) * ht% dti_sav(jat)
 
        dddsi0mt = -dddpsi0(mxt) * ht% dt3i_sav(jat)
        dddsi1mt =  dddpsi1(mxt) * ht% dt2i_sav(jat)
        dddsi2mt = -dddpsi2(mxt) * ht% dti_sav(jat)

        dddsi0d =   dddpsi0(xd) * ht% dd3i_sav(iat)
        dddsi1d =   dddpsi1(xd) * ht% dd2i_sav(iat)
        dddsi2d =   dddpsi2(xd) * ht% ddi_sav(iat)

        dddsi0md = -dddpsi0(mxd) * ht% dd3i_sav(iat)
        dddsi1md =  dddpsi1(mxd) * ht% dd2i_sav(iat)
        dddsi2md = -dddpsi2(mxd) * ht% ddi_sav(iat)


!..the free energy
        free  = h5(iat,jat,fi, &
                si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
                si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

!..first derivative with respect to density
        df_d  = h5(iat,jat,fi, &
                si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
                dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)

!..first derivative with respect to temperature
        df_t = h5(iat,jat,fi, &
                dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
                si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

!..second derivative with respect to density**2
         df_dd = h5(iat,jat,fi, &
                si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
                ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md)

!..second derivative with respect to temperature**2
        df_tt = h5(iat,jat,fi, &
              ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, &
                si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

!..second derivative with respect to temperature and density
        df_dt = h5(iat,jat,fi, &
                dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
                dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)

!..third derivative with respect to temperature**3
        df_ttt = h5(iat,jat,fi, &
                dddsi0t, dddsi1t, dddsi2t, dddsi0mt, dddsi1mt, dddsi2mt, &
                si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

!..third derivative with respect to density temperature**2
        df_dtt = h5(iat,jat,fi, &
                ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, &
                dsi0d,   dsi1d,   dsi2d,   dsi0md,   dsi1md,   dsi2md)

!..third derivative with respect to density**2 temperature
        df_ddt = h5(iat,jat,fi, &
                dsi0t, dsi1t, dsi2t, dsi0mt, dsi1mt, dsi2mt, &
                ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md)

!..third derivative with respect to density**3
        df_ddd = h5(iat,jat,fi, &
                si0t, si1t, si2t, si0mt, si1mt, si2mt, &
                dddsi0d, dddsi1d, dddsi2d, dddsi0md, dddsi1md, dddsi2md)

!..now get the pressure derivative with density, chemical potential, and 
!..electron positron number densities
!..get the cubic interpolation weight functions

        si0t   =  xpsi0(xt)
        si1t   =  xpsi1(xt) * ht% dt_sav(jat)
        si0mt  =  xpsi0(mxt)
        si1mt  =  -xpsi1(mxt) * ht% dt_sav(jat)

        si0d   =  xpsi0(xd)
        si1d   =  xpsi1(xd) * ht% dd_sav(iat)
        si0md  =  xpsi0(mxd)
        si1md  =  -xpsi1(mxd) * ht% dd_sav(iat)

!..first derivatives of weight functions
        dsi0t  = xdpsi0(xt) * ht% dti_sav(jat)
        dsi1t  = xdpsi1(xt)
        dsi0mt = -xdpsi0(mxt) * ht% dti_sav(jat)
        dsi1mt = xdpsi1(mxt)

        dsi0d  = xdpsi0(xd) * ht% ddi_sav(iat)
        dsi1d  = xdpsi1(xd)
        dsi0md = -xdpsi0(mxd) * ht% ddi_sav(iat)
        dsi1md = xdpsi1(mxd)

!..second derivatives of weight functions
        ddsi0t  = xddpsi0(xt) * ht% dt2i_sav(jat)
        ddsi1t  = xddpsi1(xt) * ht% dti_sav(jat)
        ddsi0mt = xddpsi0(mxt) * ht% dt2i_sav(jat)
        ddsi1mt = -xddpsi1(mxt) * ht% dti_sav(jat)

        ddsi0d  = xddpsi0(xd) * ht% dd2i_sav(iat)
        ddsi1d  = xddpsi1(xd) * ht% ddi_sav(iat)
        ddsi0md = xddpsi0(mxd) * ht% dd2i_sav(iat)
        ddsi1md = -xddpsi1(mxd) * ht% ddi_sav(iat)

!..third derivatives of weight functions
        dddsi0t  = xdddpsi0(xt) * ht% dt3i_sav(jat)
        dddsi1t  = xdddpsi1(xt) * ht% dt2i_sav(jat)
        dddsi0mt = -xdddpsi0(mxt) * ht% dt3i_sav(jat)
        dddsi1mt = xdddpsi1(mxt) * ht% dt2i_sav(jat)

        dddsi0d  = xdddpsi0(xd) * ht% dd3i_sav(iat)
        dddsi1d  = xdddpsi1(xd) * ht% dd2i_sav(iat)
        dddsi0md = -xdddpsi0(mxd) * ht% dd3i_sav(iat)
        dddsi1md = xdddpsi1(mxd) * ht% dd2i_sav(iat)      

!..look in the pressure derivative only once
        fi(1)  = ht% dpdf(iat,jat)
        fi(2)  = ht% dpdf(iat+1,jat)
        fi(3)  = ht% dpdf(iat,jat+1)
        fi(4)  = ht% dpdf(iat+1,jat+1)
        fi(5)  = ht% dpdft(iat,jat)
        fi(6)  = ht% dpdft(iat+1,jat)
        fi(7)  = ht% dpdft(iat,jat+1)
        fi(8)  = ht% dpdft(iat+1,jat+1)
        fi(9)  = ht% dpdfd(iat,jat)
        fi(10) = ht% dpdfd(iat+1,jat)
        fi(11) = ht% dpdfd(iat,jat+1)
        fi(12) = ht% dpdfd(iat+1,jat+1)
        fi(13) = ht% dpdfdt(iat,jat)
        fi(14) = ht% dpdfdt(iat+1,jat)
        fi(15) = ht% dpdfdt(iat,jat+1)
        fi(16) = ht% dpdfdt(iat+1,jat+1)         

!..pressure derivative with density
        dpepdd_in = h3(iat,jat,fi, &
                       si0t,   si1t,   si0mt,   si1mt, &
                       si0d,   si1d,   si0md,   si1md)
        dpepdd_in = max(dpepdd_in,1.0d-30)
        dpepdd    = dpepdd_in * dindd

!..second pressure derivative with density**2
        dpepddd_in = h3(iat,jat,fi, &
                       si0t,   si1t,   si0mt,   si1mt, &
                       dsi0d,   dsi1d, dsi0md,  dsi1md)
        dpepddd  = dpepddd_in * dindd * dindd

!..second pressure derivative with density temperature
        dpepddt_in = h3(iat,jat,fi, &
                       dsi0t,   dsi1t,   dsi0mt,   dsi1mt, &
                       si0d,   si1d, si0md,  si1md)
        dpepddt  = dpepddt_in * dindd

!..look in the electron chemical potential table only once
        fi(1)  = ht% ef(iat,jat)
        fi(2)  = ht% ef(iat+1,jat)
        fi(3)  = ht% ef(iat,jat+1)
        fi(4)  = ht% ef(iat+1,jat+1)
        fi(5)  = ht% eft(iat,jat)
        fi(6)  = ht% eft(iat+1,jat)
        fi(7)  = ht% eft(iat,jat+1)
        fi(8)  = ht% eft(iat+1,jat+1)
        fi(9)  = ht% efd(iat,jat)
        fi(10) = ht% efd(iat+1,jat)
        fi(11) = ht% efd(iat,jat+1)
        fi(12) = ht% efd(iat+1,jat+1)
        fi(13) = ht% efdt(iat,jat)
        fi(14) = ht% efdt(iat+1,jat)
        fi(15) = ht% efdt(iat,jat+1)
        fi(16) = ht% efdt(iat+1,jat+1)



!..electron chemical potential etaele
        etaele  = h3(iat,jat,fi, &
                     si0t,   si1t,   si0mt,   si1mt, &
                     si0d,   si1d,   si0md,   si1md)

!..first derivatives
        x       = h3(iat,jat,fi, &
                     si0t,   si1t,   si0mt,   si1mt, &
                    dsi0d,  dsi1d,  dsi0md,  dsi1md)
        detadd  = ye * x  
        detadt  = h3(iat,jat,fi, &
                    dsi0t,  dsi1t,  dsi0mt,  dsi1mt, &
                     si0d,   si1d,   si0md,   si1md)
       detada = -x * din * ytot1
       detadz =  x * den * ytot1

!..second derivatives
        y       = h3(iat,jat,fi, &
                    si0t,   si1t,  si0mt,  si1mt, &
                    ddsi0d,  ddsi1d,  ddsi0md,  ddsi1md)
        detaddd = ye * ye * y
        z       = h3(iat,jat,fi, &
                    dsi0t,   dsi1t,   dsi0mt,   dsi1mt, &
                    dsi0d,  dsi1d,  dsi0md,  dsi1md)
        detaddt = ye * z
        detadda = -ye*ytot1*x + ye*y*dinda
        detaddz = ytot1*x + ye*y*dindz
        detadtt   = h3(iat,jat,fi, &
                    ddsi0t,   ddsi1t,   ddsi0mt,   ddsi1mt, &
                    si0d,  si1d,  si0md,  si1md)
        detadta = z * dinda
        detadtz = z * dindz
        detadaa = (y*din + 2.0d0*x)*din*ww
        detadaz = -(y*dindz*din + x*den*ytot1)*ytot1 
        detadzz = y*dindz*den*ytot1

!..third derivatives
        y       = h3(iat,jat,fi, &
                    si0t,   si1t,  si0mt,  si1mt, &
                    dddsi0d,  dddsi1d,  dddsi0md,  dddsi1md)        
        detadddd = ye * ye * ye * y ! Actual interpolation variable is ye * rho, so we multiply by ye to get d/d(density)  
                                    !  ! d/drho^3

        detadttt = h3(iat,jat,fi, &
                    dddsi0t,   dddsi1t,  dddsi0mt,  dddsi1mt, &
                    si0d,  si1d,  si0md,  si1md)       ! d/dT^3

        y       = h3(iat,jat,fi, &
                    dsi0t,   dsi1t,  dsi0mt,  dsi1mt, &
                    ddsi0d,  ddsi1d,  ddsi0md,  ddsi1md)    
        detadddt = ye * ye * y ! d/drho^2 d/dT

        y       = h3(iat,jat,fi, &
                    ddsi0t,   ddsi1t,  ddsi0mt,  ddsi1mt, &
                    dsi0d,  dsi1d,  dsi0md,  dsi1md)    

        detaddtt = ye * y ! d/drho d/dT^2

        ! dg/da = dg/d(din) d(din)/da = dg/d(din) d(ye*rho)/da
        ! = dg/d(din) d(z rho / a)/da = -(z rho / a^2) dg/d(din)
        ! = -(z rho / a^2) dg/d(Rho) d(Rho)/d(din)
        ! = -(z rho / a^2) dg/d(Rho) (1 / ye)
        ! = -(rho / a) dg/d(Rho) = - ytot1 * rho * dg/d(Rho)
        ! = -ytot1 * den * dg/d(Rho)
        detaddda = -den * ytot1 * detadddd
        detadtta = -den * ytot1 * detaddtt
        detaddta = -den * ytot1 * detadddt

        ! dg/dz = dg/d(din) d(din)/dz = dg/d(din) d(ye*rho)/dz
        ! = dg/d(din) d(z rho / a)/dz = (rho / a) dg/d(din)
        ! = (rho / a) dg/d(Rho) d(Rho)/d(din)
        ! = (rho / a) dg/d(Rho) (1 / ye)
        ! = (rho / z) dg/d(Rho)
        ! = ytot1 * din * dg/d(Rho)
        detadddz = din * ytot1 * detadddd
        detadttz = din * ytot1 * detaddtt
        detaddtz = din * ytot1 * detadddt

!..look in the number density table only once
        fi(1)  = ht% xf(iat,jat)
        fi(2)  = ht% xf(iat+1,jat)
        fi(3)  = ht% xf(iat,jat+1)
        fi(4)  = ht% xf(iat+1,jat+1)
        fi(5)  = ht% xft(iat,jat)
        fi(6)  = ht% xft(iat+1,jat)
        fi(7)  = ht% xft(iat,jat+1)
        fi(8)  = ht% xft(iat+1,jat+1)
        fi(9)  = ht% xfd(iat,jat)
        fi(10) = ht% xfd(iat+1,jat)
        fi(11) = ht% xfd(iat,jat+1)
        fi(12) = ht% xfd(iat+1,jat+1)
        fi(13) = ht% xfdt(iat,jat)
        fi(14) = ht% xfdt(iat+1,jat)
        fi(15) = ht% xfdt(iat,jat+1)
        fi(16) = ht% xfdt(iat+1,jat+1)

!..electron + positron number densities
       xnefer   = h3(iat,jat,fi, &
                     si0t,   si1t,   si0mt,   si1mt, &
                     si0d,   si1d,   si0md,   si1md)

!..first derivatives
       x        = h3(iat,jat,fi, &
                     si0t,   si1t,   si0mt,   si1mt, &
                    dsi0d,  dsi1d,  dsi0md,  dsi1md)
       x = max(x,1.0d-30)
       dxnedd   = ye * x
       dxnedt   = h3(iat,jat,fi, &
                    dsi0t,  dsi1t,  dsi0mt,  dsi1mt, &
                     si0d,   si1d,   si0md,   si1md)
       dxneda = -x * din * ytot1
       dxnedz =  x * den * ytot1

!..second derivatives
        y       = h3(iat,jat,fi, &
                    si0t,   si1t,  si0mt,  si1mt, &
                    ddsi0d,  ddsi1d,  ddsi0md,  ddsi1md)
        dxneddd = ye * ye * y  
        z       = h3(iat,jat,fi, &
                    dsi0t,   dsi1t,   dsi0mt,   dsi1mt, &
                    dsi0d,  dsi1d,  dsi0md,  dsi1md)
        dxneddt = ye * z
        dxnedda = -ye*ytot1*x + ye*y*dinda
        dxneddz = ytot1*x + ye*y*dindz
        dxnedtt = h3(iat,jat,fi, &
                    ddsi0t,   ddsi1t,   ddsi0mt,   ddsi1mt, &
                    si0d,  si1d,  si0md,  si1md)
        dxnedta = z * dinda
        dxnedtz = z * dindz
        dxnedaa = (y*din + 2.0d0*x)*din*ytot1*ytot1
        dxnedaz = -(y*dindz*din + x*den*ytot1)*ytot1 
        dxnedzz = y*dindz*den*ytot1

!..third derivatives
        y       = h3(iat,jat,fi, &
                    si0t,   si1t,  si0mt,  si1mt, &
                    dddsi0d,  dddsi1d,  dddsi0md,  dddsi1md)        
        dxneferdddd = ye * ye * ye * y ! Actual interpolation variable is ye * rho, so we multiply by ye to get d/d(density)  
                                    !  ! d/drho^3

        dxneferdttt = h3(iat,jat,fi, &
                    dddsi0t,   dddsi1t,  dddsi0mt,  dddsi1mt, &
                    si0d,  si1d,  si0md,  si1md)       ! d/dT^3


        y       = h3(iat,jat,fi, &
                    dsi0t,   dsi1t,  dsi0mt,  dsi1mt, &
                    ddsi0d,  ddsi1d,  ddsi0md,  ddsi1md)    
        dxneferdddt = ye * ye * y ! d/drho^2 d/dT

        y       = h3(iat,jat,fi, &
                    ddsi0t,   ddsi1t,  ddsi0mt,  ddsi1mt, &
                    dsi0d,  dsi1d,  dsi0md,  dsi1md)    

        dxneferddtt = ye * y ! d/drho d/dT^2


        ! dg/da = dg/d(din) d(din)/da = dg/d(din) d(ye*rho)/da
        ! = dg/d(din) d(z rho / a)/da = -(z rho / a^2) dg/d(din)
        ! = -(z rho / a^2) dg/d(Rho) d(Rho)/d(din)
        ! = -(z rho / a^2) dg/d(Rho) (1 / ye)
        ! = -(rho / a) dg/d(Rho) = - ytot1 * rho * dg/d(Rho)
        ! = -ytot1 * den * dg/d(Rho)
        dxneferddda = -den * ytot1 * dxneferdddd
        dxneferdtta = -den * ytot1 * dxneferddtt
        dxneferddta = -den * ytot1 * dxneferdddt

        ! dg/dz = dg/d(din) d(din)/dz = dg/d(din) d(ye*rho)/dz
        ! = dg/d(din) d(z rho / a)/dz = (rho / a) dg/d(din)
        ! = (rho / a) dg/d(Rho) d(Rho)/d(din)
        ! = (rho / a) dg/d(Rho) (1 / ye)
        ! = (rho / z) dg/d(Rho)
        ! = ytot1 * din * dg/d(Rho)
        dxneferdddz = din * ytot1 * dxneferdddd
        dxneferdttz = din * ytot1 * dxneferddtt
        dxneferddtz = din * ytot1 * dxneferdddt


!..the desired electron-positron thermodynamic quantities

!..evaluating the pressure, entropy, and energy in that order
!..guarantees thermodymic consistency

!..dpepdd at high temperatures and low densities is below the
!..floating point limit of the subtraction of two large terms.
!..since dpresdd doesn't enter the maxwell relations at all, use the
!..bicubic interpolation done above instead of the formally correct expression,
!..which are commented out below.

!..pressure in erg/vol
        x       = din * din
        pele    = x * df_d

        ww = x*df_dd + 2.0d0*din*df_d
        zz = x*df_ddd + 4.0d0*din*df_dd + 2.0d0*df_d
        y = x * df_ddt + 2.0d0*din*df_dt

!..first derivatives
        dpepdt  = x * df_dt
        dpepdd  = ww*dindd
        dpepda  = ww*dinda
        dpepdz  = ww*dindz

!..second derivatives
        dpepddd = zz * dindd * dindd
        dpepddt = y*dindd
        dpepdda = (zz*dindd + ww)*dindda
        dpepddz = (zz*dindd + ww)*dinddz
        dpepdtt = x * df_dtt
        dpepdta = y * dinda
        dpepdtz = y * dindz
        dpepdaa = (zz*dinda + ww)*dindaa
        dpepdaz = (zz*dinda + ww)*dindaz
        dpepdzz = zz*dindz*dindz 



!..entropy in erg/g/k
        x       = ye * ye
        y       = ye * ytot1
        sele    = -df_t * ye

!..first derivatives
        dsepdt  = -df_tt * ye
        if (dsepdt < 0d0) then
         ierr = -1
         if (dbg) then
            write(*,*) 'dsepdt < 0d0 in helm_electron_positron'
            write(*,*) 'dsepdt', dsepdt
            write(*,*)  'df_tt', df_tt
            write(*,*) '   ye', ye
            stop
         end if
         return
        end if
        dsepdd  = -df_dt * x
        dsepda  = -df_dt*dinda*ye + df_t*y
        dsepdz  = -df_dt*dindz*ye - df_t*ytot1

!..second derivatives
        dsepddd = -df_ddt * x * ye
        dsepddt = -df_dtt * x
        dsepdda = -df_ddt * dinda * x + 2.0d0*df_dt*x*ytot1
        dsepddz = -df_ddt * dindz * x - 2.0d0*df_dt*y
        dsepdtt = -df_ttt * ye
        dsepdta = -df_dtt*dinda*ye + df_tt*y
        dsepdtz = -df_dtt * dindz * ye - df_tt*ytot1
        dsepdaa = -df_ddt*dinda*dinda*ye  &
                  - df_dt*(dindaa*ye - 2.0d0*dinda*ye*ytot1) &
                  - 2.0d0*df_t*y*ytot1
        dsepdaz = -df_ddt*dindz*dinda*ye  &
                  - df_dt*(dindaz*ye + dinda*ytot1  + dindz*dindda) &
                  - df_t*dinddaz
        dsepdzz = -df_ddt*dindz*dindz*ye  &
                  - df_dt*(dindzz*ye + dindz*ytot1 + dindz*dinddz)



!..energy in erg/g
        eele    = ye*free + temp*sele

!..first derivatives
        deepdt  = temp * dsepdt
        deepdd  = ye*df_d*dindd + temp*dsepdd
        deepda  = -y*free +  ye*df_d*dinda + temp*dsepda
        deepdz  = ytot1*free + ye*df_d*dindz + temp*dsepdz

!..second derivatives
        deepddd = ye*df_dd*dindd*dindd + temp*dsepddd
        deepddt = ye*df_dt*dindd + dsepdd + temp*dsepddt
        deepdda = -y*df_d*dindd + ye*df_dd*dinda*dindd  &
                  + ye*df_d*dindda + temp*dsepdda
        deepddz = ytot1*df_d*dindd + ye*df_dd*dindz*dindd  &
                  + ye*df_d*dinddz + temp*dsepddz
        deepdtt = dsepdt + temp*dsepdtt
        deepdta = temp * dsepdta
        deepdtz = temp * dsepdtz
        deepdaa  = 2.0d0*y*(ytot1*free - df_d*dinda)  &
                    + ye*(df_d*dindaa + df_dd*dinda*dinda)  &
                    + temp*dsepdaa
        deepdaz  = -ytot1*ytot1*free - y*df_d*dindz +  ytot1*df_d*dinda  &
                  + ye*df_dd*dindz*dinda + ye*df_d*dindaz + temp*dsepdaz
        deepdzz  = (2.0d0*ytot1*df_d + ye*df_dd*dindz)*dindz  &
                   + temp*dsepdzz


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
