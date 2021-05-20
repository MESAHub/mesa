! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

      module mod_neu
      use const_def
      use neu_def
      use math_lib
      use utils_lib, only: mesa_error

      implicit none


      !..various numerical constants
      
      real(dp), parameter :: con1   = 1.0d0/5.9302d0

      !..cv and ca are the vector and axial currents.

      real(dp), parameter :: cv     = 0.5d0 + 2.0d0 * weinberg_theta
      real(dp), parameter :: cvp    = 1.0d0 - cv
      real(dp), parameter :: ca     = 0.5d0
      real(dp), parameter :: cap    = 1.0d0 - ca
      real(dp), parameter :: tfac1  = cv*cv + ca*ca + (num_neu_fam-1.0d0) * (cvp*cvp+cap*cap)
      real(dp), parameter :: tfac2  = cv*cv - ca*ca + (num_neu_fam-1.0d0) * (cvp*cvp - cap*cap)
      real(dp), parameter :: tfac3  = tfac2/tfac1
      real(dp), parameter :: tfac4  = 0.5d0 * tfac1
      real(dp), parameter :: tfac5  = 0.5d0 * tfac2
      real(dp), parameter :: tfac6  = cv*cv + 1.5d0*ca*ca + (num_neu_fam - 1.0d0)*(cvp*cvp + 1.5d0*cap*cap)



      real(dp), parameter :: fac1   = 5.0d0 * pi / 3.0d0
      real(dp), parameter :: fac2   = 10.0d0 * pi
      real(dp), parameter :: fac3   = pi / 5.0d0


      type t8s
         real(dp) :: t8,t812,t832,t82,t83,t85,t86,t8m1,t8m2,t8m3,t8m5,t8m6
      end type t8s


      type inputs
         real(dp) :: temp,logtemp,den,logden,den6
         real(dp) :: ye,deni,tempi,abari,zbari,abar,zbar
         real(dp) :: t9,xl,xldt,xlp5,xl2,xl3,xl4,xl5,xl6,xl7,xl8,xl9,xlmp5,xlm1,xlm2,xlm3,xlm4
         real(dp) :: rm,rmdd,rmda,rmdz,rmi
         real(dp) :: zeta,zetadt,zetadd,zetada,zetadz,zeta2,zeta3
      end type inputs


      contains

      real(dp) function ifermi12(f)

!..this routine applies a rational function expansion to get the inverse
!..fermi-dirac integral of order 1/2 when it is equal to f.
!..maximum error is 4.19d-9.   reference: antia apjs 84,101 1993

!..declare
      integer          i,m1,k1,m2,k2
      real(dp) f,an,a1(12),b1(12),a2(12),b2(12),rn,den,ff


!..load the coefficients of the expansion
      data  an,m1,k1,m2,k2 /0.5d0, 4, 3, 6, 5/
      data  (a1(i),i=1,5)/ 1.999266880833d4,   5.702479099336d3, &
           6.610132843877d2,   3.818838129486d1, &
           1.0d0/
      data  (b1(i),i=1,4)/ 1.771804140488d4,  -2.014785161019d3, &
           9.130355392717d1,  -1.670718177489d0/
      data  (a2(i),i=1,7)/-1.277060388085d-2,  7.187946804945d-2,  &
                          -4.262314235106d-1,  4.997559426872d-1, &
                          -1.285579118012d0,  -3.930805454272d-1, &
           1.0d0/
      data  (b2(i),i=1,6)/-9.745794806288d-3,  5.485432756838d-2, &
                          -3.299466243260d-1,  4.077841975923d-1, &
                          -1.145531476975d0,  -6.067091689181d-2/


      if (f .lt. 4.0d0) then
         rn  = f + a1(m1)
         do i=m1-1,1,-1
            rn  = rn*f + a1(i)
         enddo
         den = b1(k1+1)
         do i=k1,1,-1
            den = den*f + b1(i)
         enddo
         ifermi12 = log(f * rn/den)

      else
         ff = 1.0d0/pow(f,1.0d0/(1.0d0 + an))
         rn = ff + a2(m2)
         do i=m2-1,1,-1
            rn = rn*ff + a2(i)
         enddo
         den = b2(k2+1)
         do i=k2,1,-1
            den = den*ff + b2(i)
         enddo
         ifermi12 = rn/(den*ff)
      end if

      end function ifermi12






      real(dp) function zfermim12(x)

!..this routine applies a rational function expansion to get the fermi-dirac
!..integral of order -1/2 evaluated at x. maximum error is 1.23d-12.
!..reference: antia apjs 84,101 1993

!..declare
      integer          i,m1,k1,m2,k2
      real(dp) x,an,a1(12),b1(12),a2(12),b2(12),rn,den,xx

!..load the coefficients of the expansion
      data  an,m1,k1,m2,k2 /-0.5d0, 7, 7, 11, 11/
      data  (a1(i),i=1,8)/ 1.71446374704454d7,    3.88148302324068d7, &
                           3.16743385304962d7,    1.14587609192151d7, &
                           1.83696370756153d6,    1.14980998186874d5, &
                           1.98276889924768d3,    1.0d0/
      data  (b1(i),i=1,8)/ 9.67282587452899d6,    2.87386436731785d7, &
                           3.26070130734158d7,    1.77657027846367d7, &
                           4.81648022267831d6,    6.13709569333207d5, &
                           3.13595854332114d4,    4.35061725080755d2/
      data (a2(i),i=1,12)/-4.46620341924942d-15, -1.58654991146236d-12, &
                          -4.44467627042232d-10, -6.84738791621745d-8, &
                          -6.64932238528105d-6,  -3.69976170193942d-4, &
                          -1.12295393687006d-2,  -1.60926102124442d-1, &
                          -8.52408612877447d-1,  -7.45519953763928d-1, &
                           2.98435207466372d0,    1.0d0/
      data (b2(i),i=1,12)/-2.23310170962369d-15, -7.94193282071464d-13, &
                          -2.22564376956228d-10, -3.43299431079845d-8, &
                          -3.33919612678907d-6,  -1.86432212187088d-4, &
                          -5.69764436880529d-3,  -8.34904593067194d-2, &
                          -4.78770844009440d-1,  -4.99759250374148d-1, &
                           1.86795964993052d0,    4.16485970495288d-1/


      if (x .lt. 2.0d0) then
         xx = exp(x)
         rn = xx + a1(m1)
         do i=m1-1,1,-1
            rn = rn*xx + a1(i)
         enddo
         den = b1(k1+1)
         do i=k1,1,-1
            den = den*xx + b1(i)
         enddo
         zfermim12 = xx * rn/den
!..
      else
         xx = 1.0d0/(x*x)
         rn = xx + a2(m2)
         do i=m2-1,1,-1
            rn = rn*xx + a2(i)
         enddo
         den = b2(k2+1)
         do i=k2,1,-1
            den = den*xx + b2(i)
         enddo
         zfermim12 = sqrt(x)*rn/den
      end if

      end function zfermim12

      

      subroutine neutrinos(T, logT, Rho, logRho, abar, zbar, log10_Tlim,  &
                  flags, loss, sources, info)
      use utils_lib, only: is_bad
      
      !..this routine computes neutrino losses from the analytic fits of
      !..itoh et al. apjs 102, 411, 1996, and also returns their derivatives. 
      
      ! provide T or logT or both (the code needs both, so pass 'em if you've got 'em!)
      ! same for Rho and logRho
      
      real(dp), intent(in) :: T ! temperature
      real(dp), intent(in) :: logT ! log10 of temperature
      real(dp), intent(in) :: Rho ! density
      real(dp), intent(in) :: logRho ! log10 of density
      real(dp), intent(in) :: abar ! mean atomic weight
      real(dp), intent(in) :: zbar ! mean charge
      real(dp), intent(in) :: log10_Tlim ! start to cutoff at this temperature
      logical, intent(in) :: flags(num_neu_types) ! true if should include the type
      ! in most cases for stellar evolution, you may want to include brem, plas, pair, and phot
      ! but skip reco.  it is fairly expensive to compute and typically makes only a small contribution.

      real(dp), intent(inout) :: loss(num_neu_rvs) ! total from all sources
      real(dp), intent(inout) :: sources(num_neu_types, num_neu_rvs)
      integer, intent(out) :: info ! 0 means AOK.

!..local variables
      type(inputs) ::  input

      real(dp) :: temp,logtemp,den,logden,tcutoff_factor

      real(dp) :: snu, dsnudt, dsnudd, dsnuda, dsnudz, dtcutoff_factordt, dtlim
      real(dp) spair,spairdt,spairdd,spairda,spairdz, &
                       splas,splasdt,splasdd,splasda,splasdz, &
                       sphot,sphotdt,sphotdd,sphotda,sphotdz, &
                       sbrem,sbremdt,sbremdd,sbremda,sbremdz, &
                       sreco,srecodt,srecodd,srecoda,srecodz


      info = 0
   
      if ((T /= arg_not_provided .and. T .le. Tmin_neu) .or. &
            (logT /= arg_not_provided .and. logT .le. log10Tmin_neu)) then
         loss = 0d0
         sources(1:num_neu_types, 1:num_neu_rvs) = 0d0
         return
      end if

      if (T == arg_not_provided) then
         if (logT == arg_not_provided) then
            info = -2
            return
         end if
         temp = exp10(logT)
      else
         temp = T
      end if
         
      if (T <= 0) then
         info = -1
         return
      end if
      
      if (logT == arg_not_provided) then
         logtemp = log10(T)
      else
         logtemp = logT
      end if
      
      if (logtemp > 20) then
         info = -1
         return
      end if

      if (Rho == arg_not_provided) then
         if (logRho == arg_not_provided) then
            info = -3
            return
         end if
         den = exp10(logRho)
      else
         den = Rho
      end if
         
      if (Rho <= 0) then
         info = -1
         return
      end if

      if (logRho == arg_not_provided) then
         logden = log10(Rho)
      else
         logden = logRho
      end if
      
      if (logden > 20) then
         info = -1
         return
      end if

!..initialize 
      spair   = 0.0d0
      spairdt = 0.0d0
      spairdd = 0.0d0
      spairda = 0.0d0
      spairdz = 0.0d0

      splas   = 0.0d0
      splasdt = 0.0d0
      splasdd = 0.0d0
      splasda = 0.0d0
      splasdz = 0.0d0

      sphot   = 0.0d0
      sphotdt = 0.0d0
      sphotdd = 0.0d0
      sphotda = 0.0d0
      sphotdz = 0.0d0

      sbrem   = 0.0d0
      sbremdt = 0.0d0
      sbremdd = 0.0d0
      sbremda = 0.0d0
      sbremdz = 0.0d0

      sreco   = 0.0d0
      srecodt = 0.0d0
      srecodd = 0.0d0
      srecoda = 0.0d0
      srecodz = 0.0d0

      snu     = 0.0d0
      dsnudt  = 0.0d0
      dsnudd  = 0.0d0
      dsnuda  = 0.0d0
      dsnudz  = 0.0d0


      call set_inputs(input,temp,logtemp,den,logden,abar,zbar)
      
!..do the requested types

      if (flags(pair_neu_type)) call pair_neu(spair,spairdt,spairdd,spairda,spairdz, input)
      if (flags(plas_neu_type)) call plas_neu(splas,splasdt,splasdd,splasda,splasdz, input)
      if (flags(phot_neu_type)) call phot_neu(sphot,sphotdt,sphotdd,sphotda,sphotdz, input)
      if (flags(brem_neu_type)) call brem_neu(sbrem,sbremdt,sbremdd,sbremda,sbremdz, input)
      if (flags(reco_neu_type)) call reco_neu(sreco,srecodt,srecodd,srecoda,srecodz, input)

!..convert from erg/cm^3/s to erg/g/s 
!..comment these out to duplicate the itoh et al plots

      spair   = spair   * input% deni
      spairdt = spairdt * input% deni
      spairdd = spairdd * input% deni - spair * input% deni
      spairda = spairda * input% deni
      spairdz = spairdz * input% deni  

      splas   = splas   * input% deni
      splasdt = splasdt * input% deni
      splasdd = splasdd * input% deni - splas * input% deni
      splasda = splasda * input% deni
      splasdz = splasdz * input% deni  

      sphot   = sphot   * input% deni
      sphotdt = sphotdt * input% deni
      sphotdd = sphotdd * input% deni - sphot * input% deni
      sphotda = sphotda * input% deni
      sphotdz = sphotdz * input% deni  

      sbrem   = sbrem   * input% deni
      sbremdt = sbremdt * input% deni
      sbremdd = sbremdd * input% deni - sbrem * input% deni
      sbremda = sbremda * input% deni
      sbremdz = sbremdz * input% deni  

      sreco   = sreco   * input% deni
      srecodt = srecodt * input% deni
      srecodd = srecodd * input% deni - sreco * input% deni
      srecoda = srecoda * input% deni
      srecodz = srecodz * input% deni  

!..calculate temperature cutoff factor

      if (input% logtemp <= log10_Tlim .and. log10_Tlim > log10Tmin_neu) then
          dtlim = log10_Tlim - log10Tmin_neu
         tcutoff_factor = 0.5d0* &
            (1 - cospi((input% logtemp - log10Tmin_neu)/(log10_Tlim - log10Tmin_neu)))
     
     
         dtcutoff_factordt = 0.5d0 * pi * sinpi((input% logtemp - log10Tmin_neu)/dtlim) * &
                             1.d0/(dtlim * temp * ln10)

         splasdt = tcutoff_factor * splasdt + dtcutoff_factordt * splas
         splasdd = tcutoff_factor * splasdd
         splasda = tcutoff_factor * splasda
         splasdz = tcutoff_factor * splasdz
         splas   = tcutoff_factor * splas
     
         spairdt = tcutoff_factor * spairdt + dtcutoff_factordt * spair
         spairdd = tcutoff_factor * spairdd
         spairda = tcutoff_factor * spairda
         spairdz = tcutoff_factor * spairdz
         spair   = tcutoff_factor * spair     
     
         sphotdt = tcutoff_factor * sphotdt + dtcutoff_factordt * sphot
         sphotdd = tcutoff_factor * sphotdd
         sphotda = tcutoff_factor * sphotda
         sphotdz = tcutoff_factor * sphotdz
         sphot   = tcutoff_factor * sphot
     
         sbremdt = tcutoff_factor * sbremdt + dtcutoff_factordt * sbrem
         sbremdd = tcutoff_factor * sbremdd
         sbremda = tcutoff_factor * sbremda
         sbremdz = tcutoff_factor * sbremdz
         sbrem   = tcutoff_factor * sbrem
     
         srecodt = tcutoff_factor * srecodt + dtcutoff_factordt * sreco
         srecodd = tcutoff_factor * srecodd
         srecoda = tcutoff_factor * srecoda
         srecodz = tcutoff_factor * srecodz
         sreco   = tcutoff_factor * sreco

         
      end if


!..the total neutrino loss rate
      snu    =  splas   + spair   + sphot   + sbrem   + sreco
      dsnudt =  splasdt + spairdt + sphotdt + sbremdt + srecodt  
      dsnudd =  splasdd + spairdd + sphotdd + sbremdd + srecodd 
      dsnuda =  splasda + spairda + sphotda + sbremda + srecoda 
      dsnudz =  splasdz + spairdz + sphotdz + sbremdz + srecodz 
      
      if (is_bad(snu)) then
         info = -1
         return
      end if

!..packup the results

      loss(1) = snu
      loss(2) = dsnudt
      loss(3) = dsnudd
      loss(4) = dsnuda
      loss(5) = dsnudz
      
      call store(pair_neu_type, spair, spairdt, spairdd, spairda, spairdz)
      call store(plas_neu_type, splas, splasdt, splasdd, splasda, splasdz)
      call store(phot_neu_type, sphot, sphotdt, sphotdd, sphotda, sphotdz)
      call store(brem_neu_type, sbrem, sbremdt, sbremdd, sbremda, sbremdz)
      call store(reco_neu_type, sreco, srecodt, srecodd, srecoda, srecodz)
      

      contains
      
      
         subroutine store(neu_type, s, sdt, sdd, sda, sdz)
            integer, intent(in) :: neu_type
            real(dp), intent(in) :: s, sdt, sdd, sda, sdz
            sources(neu_type,1) = s
            sources(neu_type,2) = sdt
            sources(neu_type,3) = sdd
            sources(neu_type,4) = sda
            sources(neu_type,5) = sdz
         end subroutine store
         
      end subroutine neutrinos


      subroutine set_inputs(input,temp,logtemp,den,logden,abar,zbar)
         type(inputs),intent(out) :: input

         real(dp), intent(in) :: temp,logtemp,den,logden,abar, zbar
         real(dp) ::a0,a1,a2

         input% temp = temp
         input% logtemp = logtemp
         input% den = den
         input% logden = logden
         input% abar = abar
         input% zbar = zbar
         input% den6 = input% den * 1.0d-6

!..to avoid lots of divisions
         input% deni  = 1.0d0 / input% den
         input% tempi = 1.0d0 / input% temp
         input% abari = 1.0d0 / input% abar
         input% zbari  = 1.0d0 / input% zbar
   
   
   !..some composition variables
         input% ye    = input% zbar * input% abari
         !xmue  = abar * zbari
   
   
   !..some frequent factors
         input% t9     = input% temp * 1.0d-9
         input% xl     = input% t9 * con1
         input% xldt   = 1.0d-9 * con1
         input% xlp5   = sqrt(input% xl)
         input% xl2    = input% xl  * input% xl
         input% xl3    = input% xl2 * input% xl
         input% xl4    = input% xl3 * input% xl
         input% xl5    = input% xl4 * input% xl
         input% xl6    = input% xl5 * input% xl
         input% xl7    = input% xl6 * input% xl
         input% xl8    = input% xl7 * input% xl
         input% xl9    = input% xl8 * input% xl
         input% xlmp5  = 1.0d0 / input% xlp5
         input% xlm1   = 1.0d0 / input% xl
         input% xlm2   = input% xlm1*input% xlm1
         input% xlm3   = input% xlm1*input% xlm2
         input% xlm4   = input% xlm1*input% xlm3      
         
         
         input% rm     = input% den*input% ye
         input% rmdd   = input% ye
         input% rmda   = -input% rm*input% abari
         input% rmdz   = input% den *input% abari
         input% rmi    = 1.0d0/input% rm


         a0     = input% rm * 1.0d-9
         a1     = pow(a0,one_third)
         input% zeta   = a1 * input% xlm1
         input% zetadt = -a1 * input% xlm2 * input% xldt
         a2     = one_third * a1*input%rmi * input% xlm1
         input% zetadd = a2 * input%rmdd 
         input% zetada = a2 * input%rmda
         input% zetadz = a2 * input%rmdz
         
         input% zeta2 = input%zeta * input%zeta
         input% zeta3 = input%zeta2 * input%zeta


      end subroutine set_inputs


      subroutine phot_neu(sphot,sphotdt,sphotdd,sphotda,sphotdz, input)
         real(dp), intent(out) :: sphot,sphotdt,sphotdd,sphotda,sphotdz
         type(inputs), intent(in) :: input
   
         real(dp) tau,taudt,cos1,cos2,cos3,cos4,cos5,sin1,sin2, &
            sin3,sin4,sin5,last,xast, &
            fphot,fphotdt,fphotdd,fphotda,fphotdz, &
            qphot,qphotdt,qphotdd,qphotda,qphotdz

         real(dp) :: a0,a1,a2,a3,c00,c01,c02,c03,c04,c05,c06,cc, &
            c10,c11,c12,c13,c14,c15,c16,c20,c21,c22,c23,c24, &
            c25,c26,dd01,dd02,dd03,dd04,dd05,dd11,dd12, &
            dd13,dd14,dd15,dd21,dd22,dd23,dd24,dd25,f0,f1,f2,z, &
            dum,dumdt,dumdd,dumda,dumdz

         real(dp) :: xnum,xnumdt,xnumdd,xnumda,xnumdz, &
            xden,xdendt,xdendd,xdenda,xdendz

         real(dp) :: dccdt

   !..photoneutrino process section  
   !..for reactions like e- + gamma => e- + nu_e + nubar_e
   !..                   e+ + gamma => e+ + nu_e + nubar_e
   !..equation 3.8 for tau, equation 3.6 for cc,
   !..and table 2 written out for speed
         dccdt = 0d0

         sphot   = 0.0d0
         sphotdt = 0.0d0
         sphotdd = 0.0d0
         sphotda = 0.0d0
         sphotdz = 0.0d0

         if(input% temp .le.1.0d7) then
            return
         else if (input% temp .ge. 1.0d7  .and. input% temp .lt. 1.0d8) then
            tau  =  input% logtemp - 7d0
            cc   =  0.5654d0 + tau
            dccdt = 1d0/(ln10*input% temp)
            c00  =  1.008d11
            c01  =  0.0d0
            c02  =  0.0d0
            c03  =  0.0d0
            c04  =  0.0d0
            c05  =  0.0d0
            c06  =  0.0d0
            c10  =  8.156d10
            c11  =  9.728d8
            c12  = -3.806d9
            c13  = -4.384d9
            c14  = -5.774d9
            c15  = -5.249d9
            c16  = -5.153d9
            c20  =  1.067d11
            c21  = -9.782d9 
            c22  = -7.193d9
            c23  = -6.936d9
            c24  = -6.893d9
            c25  = -7.041d9
            c26  = -7.193d9
            dd01 =  0.0d0
            dd02 =  0.0d0
            dd03 =  0.0d0
            dd04 =  0.0d0
            dd05 =  0.0d0
            dd11 = -1.879d10
            dd12 = -9.667d9
            dd13 = -5.602d9
            dd14 = -3.370d9
            dd15 = -1.825d9
            dd21 = -2.919d10
            dd22 = -1.185d10
            dd23 = -7.270d9
            dd24 = -4.222d9
            dd25 = -1.560d9

         else if (input% temp .ge. 1.0d8  .and. input% temp .lt. 1.0d9) then
            tau   =  input% logtemp - 8d0
            cc   =  1.5654d0
            dccdt = 0d0
            c00  =  9.889d10 
            c01  = -4.524d8
            c02  = -6.088d6 
            c03  =  4.269d7 
            c04  =  5.172d7 
            c05  =  4.910d7 
            c06  =  4.388d7
            c10  =  1.813d11
            c11  = -7.556d9 
            c12  = -3.304d9  
            c13  = -1.031d9
            c14  = -1.764d9  
            c15  = -1.851d9
            c16  = -1.928d9
            c20  =  9.750d10
            c21  =  3.484d10
            c22  =  5.199d9  
            c23  = -1.695d9  
            c24  = -2.865d9  
            c25  = -3.395d9  
            c26  = -3.418d9
            dd01 = -1.135d8   
            dd02 =  1.256d8   
            dd03 =  5.149d7   
            dd04 =  3.436d7   
            dd05 =  1.005d7
            dd11 =  1.652d9  
            dd12 = -3.119d9  
            dd13 = -1.839d9  
            dd14 = -1.458d9  
            dd15 = -8.956d8
            dd21 = -1.548d10
            dd22 = -9.338d9  
            dd23 = -5.899d9  
            dd24 = -3.035d9  
            dd25 = -1.598d9

         else if (input% temp .ge. 1.0d9) then
            tau  =  input% logtemp - 9d0
            cc   =  1.5654d0
            dccdt = 0d0
            c00  =  9.581d10
            c01  =  4.107d8
            c02  =  2.305d8   
            c03  =  2.236d8   
            c04  =  1.580d8   
            c05  =  2.165d8   
            c06  =  1.721d8
            c10  =  1.459d12
            c11  =  1.314d11
            c12  = -1.169d11  
            c13  = -1.765d11  
            c14  = -1.867d11  
            c15  = -1.983d11  
            c16  = -1.896d11
            c20  =  2.424d11
            c21  = -3.669d9
            c22  = -8.691d9  
            c23  = -7.967d9  
            c24  = -7.932d9  
            c25  = -7.987d9  
            c26  = -8.333d9
            dd01 =  4.724d8
            dd02 =  2.976d8   
            dd03 =  2.242d8   
            dd04 =  7.937d7   
            dd05 =  4.859d7
            dd11 = -7.094d11
            dd12 = -3.697d11
            dd13 = -2.189d11  
            dd14 = -1.273d11  
            dd15 = -5.705d10
            dd21 = -2.254d10
            dd22 = -1.551d10
            dd23 = -7.793d9
            dd24 = -4.489d9
         dd25 = -2.185d9
         end if

         taudt = iln10*input% tempi


   !..equation 3.7, compute the expensive trig functions only one time
         cos1 = cos(fac1*tau)
         cos2 = cos(fac1*2.0d0*tau)
         cos3 = cos(fac1*3.0d0*tau)
         cos4 = cos(fac1*4.0d0*tau)
         cos5 = cos(fac1*5.0d0*tau)
         last = cos(fac2*tau)

         sin1 = sin(fac1*tau)
         sin2 = sin(fac1*2.0d0*tau)
         sin3 = sin(fac1*3.0d0*tau)
         sin4 = sin(fac1*4.0d0*tau)
         sin5 = sin(fac1*5.0d0*tau)
         xast = sin(fac2*tau)

         a0 = 0.5d0*c00  &
            + c01*cos1 + dd01*sin1 + c02*cos2 + dd02*sin2 &
            + c03*cos3 + dd03*sin3 + c04*cos4 + dd04*sin4 &
            + c05*cos5 + dd05*sin5 + 0.5d0*c06*last

         f0 =  taudt*fac1*(-c01*sin1 + dd01*cos1 - c02*sin2*2.0d0  &
            + dd02*cos2*2.0d0 - c03*sin3*3.0d0 + dd03*cos3*3.0d0  &
            - c04*sin4*4.0d0 + dd04*cos4*4.0d0 &
            - c05*sin5*5.0d0 + dd05*cos5*5.0d0)  &
            - 0.5d0*c06*xast*fac2*taudt

         a1 = 0.5d0*c10  &
            + c11*cos1 + dd11*sin1 + c12*cos2 + dd12*sin2 &
            + c13*cos3 + dd13*sin3 + c14*cos4 + dd14*sin4 &
            + c15*cos5 + dd15*sin5 + 0.5d0*c16*last

         f1 = taudt*fac1*(-c11*sin1 + dd11*cos1 - c12*sin2*2.0d0  &
            + dd12*cos2*2.0d0 - c13*sin3*3.0d0 + dd13*cos3*3.0d0  &
            - c14*sin4*4.0d0 + dd14*cos4*4.0d0 - c15*sin5*5.0d0  &
            + dd15*cos5*5.0d0) - 0.5d0*c16*xast*fac2*taudt

         a2 = 0.5d0*c20  &
            + c21*cos1 + dd21*sin1 + c22*cos2 + dd22*sin2 &
            + c23*cos3 + dd23*sin3 + c24*cos4 + dd24*sin4 &
            + c25*cos5 + dd25*sin5 + 0.5d0*c26*last

         f2 = taudt*fac1*(-c21*sin1 + dd21*cos1 - c22*sin2*2.0d0  &
            + dd22*cos2*2.0d0 - c23*sin3*3.0d0 + dd23*cos3*3.0d0  &
            - c24*sin4*4.0d0 + dd24*cos4*4.0d0 - c25*sin5*5.0d0  &
            + dd25*cos5*5.0d0) - 0.5d0*c26*xast*fac2*taudt

   !..equation 3.4
         dum   = a0 + a1*input% zeta + a2*input% zeta2
         dumdt = f0 + f1*input% zeta + a1*input% zetadt + f2*input% zeta2 + 2.0d0*a2*input% zeta*input% zetadt
         dumdd = a1*input% zetadd + 2.0d0*a2*input% zeta*input% zetadd
         dumda = a1*input% zetada + 2.0d0*a2*input% zeta*input% zetada
         dumdz = a1*input% zetadz + 2.0d0*a2*input% zeta*input% zetadz

         z      = exp(-cc*input% zeta)

         xnum   = dum*z
         xnumdt = dumdt*z - dum*z*(dccdt*input% zeta + input% zetadt*cc)
         xnumdd = dumdd*z - dum*z*cc*input% zetadd
         xnumda = dumda*z - dum*z*cc*input% zetada
         xnumdz = dumdz*z - dum*z*cc*input% zetadz

         xden   = input% zeta3 + 6.290d-3*input% xlm1 + 7.483d-3*input% xlm2 + 3.061d-4*input% xlm3

         dum    = 3.0d0*input% zeta2
         xdendt = dum*input% zetadt - input% xldt*(6.290d-3*input% xlm2  &
                  + 2.0d0*7.483d-3*input% xlm3 + 3.0d0*3.061d-4*input% xlm4)
         xdendd = dum*input% zetadd
         xdenda = dum*input% zetada
         xdendz = dum*input% zetadz

         dum      = 1.0d0/xden
         fphot   = xnum*dum
         fphotdt = (xnumdt - fphot*xdendt)*dum
         fphotdd = (xnumdd - fphot*xdendd)*dum
         fphotda = (xnumda - fphot*xdenda)*dum
         fphotdz = (xnumdz - fphot*xdendz)*dum
   

   !..equation 3.3
         a0     = 1.0d0 + 2.045d0 * input% xl
         xnum   = 0.666d0*pow(a0,-2.066d0)
         xnumdt = -2.066d0*xnum/a0 * 2.045d0*input% xldt

         dum    = 1.875d8*input% xl + 1.653d8*input% xl2 + 8.499d8*input% xl3 - 1.604d8*input% xl4
         dumdt  = input% xldt*(1.875d8 + 2.0d0*1.653d8*input% xl + 3.0d0*8.499d8*input% xl2  &
                  - 4.0d0*1.604d8*input% xl3)

         z      = 1.0d0/dum
         xden   = 1.0d0 + input% rm*z
         xdendt =  -input% rm*z*z*dumdt
         xdendd =  input% rmdd*z
         xdenda =  input% rmda*z
         xdendz =  input% rmdz*z

         z      = 1.0d0/xden
         qphot = xnum*z
         qphotdt = (xnumdt - qphot*xdendt)*z
         dum      = -qphot*z
         qphotdd = dum*xdendd
         qphotda = dum*xdenda
         qphotdz = dum*xdendz

   !..equation 3.2
         sphot   = input% xl5 * fphot
         sphotdt = 5.0d0*input% xl4*input% xldt*fphot + input% xl5*fphotdt
         sphotdd = input% xl5*fphotdd
         sphotda = input% xl5*fphotda
         sphotdz = input% xl5*fphotdz

         a1      = sphot
         sphot   = input% rm*a1
         sphotdt = input% rm*sphotdt  
         sphotdd = input% rm*sphotdd + input% rmdd*a1  
         sphotda = input% rm*sphotda + input% rmda*a1  
         sphotdz = input% rm*sphotdz + input% rmdz*a1  

         a1      = tfac4*(1.0d0 - tfac3 * qphot)
         a2      = -tfac4*tfac3

         a3      = sphot
         sphot   = a1*a3
         sphotdt = a1*sphotdt + a2*qphotdt*a3
         sphotdd = a1*sphotdd + a2*qphotdd*a3
         sphotda = a1*sphotda + a2*qphotda*a3
         sphotdz = a1*sphotdz + a2*qphotdz*a3

            
         if (.false.) then
            write(*,*) 'logT', input% logtemp
            write(*,*) 'logRho', input% logden
            write(*,*) 'sphot', sphot
            write(*,*) 
         end if

         if (sphot .le. 0.0d0) then
            sphot   = 0.0d0
            sphotdt = 0.0d0
            sphotdd = 0.0d0
            sphotda = 0.0d0
            sphotdz = 0.0d0
         end if

      end subroutine phot_neu
      
      
      subroutine brem_neu_weak_degen(sbrem,sbremdt,sbremdd,sbremda,sbremdz,t8, input)
         type(t8s), intent(in) :: t8
         real(dp), intent(out) :: sbrem,sbremdt,sbremdd,sbremda,sbremdz
         type(inputs), intent(in) :: input

         real(dp) :: a0,a1,a2,c00,c01,c02,c03,c04, &
         dd00,dd01,dd02,&
         f0,z, &
         dum,dumdt,dumdd,dumda,dumdz

         real(dp) :: xnum,xnumdt,xnumdd,xnumda,xnumdz, &
         xden,xdendt,xdendd,xdenda,xdendz

   ! brem

         real(dp) :: eta,etadt,etadd,etada,etadz,etam1,etam2,etam3, &
         fbrem,fbremdt,fbremdd,fbremda,fbremdz, &
         gbrem,gbremdt,gbremdd,gbremda,gbremdz
         
         real(dp) :: p

         sbrem=0d0
         sbremdt=0d0
         sbremdd=0d0
         sbremda=0d0
         sbremdz=0d0
         
   !..equation 5.3
         dum   = 7.05d6 * t8% t832 + 5.12d4 * t8% t83
         dumdt = (1.5d0*7.05d6*t8% t812 + 3.0d0*5.12d4*t8% t82)*1.0d-8

         z     = 1.0d0/dum
         eta   = input% rm*z
         etadt = -input% rm*z*z*dumdt
         etadd = input% rmdd*z
         etada = input% rmda*z
         etadz = input% rmdz*z

         etam1 = 1.0d0/eta
         etam2 = etam1 * etam1
         etam3 = etam2 * etam1


   !..equation 5.2
         a0    = 23.5d0 + 6.83d4*t8% t8m2 + 7.81d8*t8% t8m5
         f0    = (-2.0d0*6.83d4*t8% t8m3 - 5.0d0*7.81d8*t8% t8m6)*1.0d-8
         xnum  = 1.0d0/a0

         dum   = 1.0d0 + 1.47d0*etam1 + 3.29d-2*etam2
         z     = -1.47d0*etam2 - 2.0d0*3.29d-2*etam3
         dumdt = z*etadt
         dumdd = z*etadd
         dumda = z*etada
         dumdz = z*etadz

         c00   = 1.26d0 * (1.0d0+etam1)
         z     = -1.26d0*etam2
         c01   = z*etadt
         c02   = z*etadd
         c03   = z*etada
         c04   = z*etadz
         
         z      = 1.0d0/dum
         xden   = c00*z
         xdendt = (c01 - xden*dumdt)*z
         xdendd = (c02 - xden*dumdd)*z
         xdenda = (c03 - xden*dumda)*z
         xdendz = (c04 - xden*dumdz)*z

         fbrem   = xnum + xden
         fbremdt = -xnum*xnum*f0 + xdendt
         fbremdd = xdendd
         fbremda = xdenda
         fbremdz = xdendz



   !..equation 5.9
         a0    = 230.0d0 + 6.7d5*t8% t8m2 + 7.66d9*t8% t8m5
         f0    = (-2.0d0*6.7d5*t8% t8m3 - 5.0d0*7.66d9*t8% t8m6)*1.0d-8

         z     = 1.0d0 + input% rm*1.0d-9 
         dum   = a0*z
         dumdt = f0*z
         z     = a0*1.0d-9
         dumdd = z*input% rmdd
         dumda = z*input% rmda
         dumdz = z*input% rmdz

         xnum   = 1.0d0/dum
         z      = -xnum*xnum
         xnumdt = z*dumdt
         xnumdd = z*dumdd
         xnumda = z*dumda
         xnumdz = z*dumdz
         
         p = pow(t8% t8,3.85d0)
         c00   = 7.75d5*t8% t832 + 247.0d0*p
         dd00  = (1.5d0*7.75d5*t8% t812 + 3.85d0*247.0d0*p/t8% T8)*1.0d-8
         
         p = pow(t8% t8,1.4d0)
         c01   = 4.07d0 + 0.0240d0 * p
         dd01  = 1.4d0*0.0240d0*(p/t8% T8)*1.0d-8
         
         p = pow(t8% t8,-0.110d0)
         c02   = 4.59d-5 * p
         dd02  = -0.11d0*4.59d-5*(p/t8% T8)*1.0d-8

         z     = pow(input% den,0.656d0)
         dum   = c00*input% rmi  + c01  + c02*z 
         dumdt = dd00*input% rmi + dd01 + dd02*z
         z     = -c00*input% rmi*input% rmi
         dumdd = z*input% rmdd + 0.656d0*c02*pow(input% den,-0.344d0)
         dumda = z*input% rmda 
         dumdz = z*input% rmdz 

         xden  = 1.0d0/dum
         z      = -xden*xden
         xdendt = z*dumdt
         xdendd = z*dumdd
         xdenda = z*dumda
         xdendz = z*dumdz


         gbrem   = xnum + xden
         gbremdt = xnumdt + xdendt
         gbremdd = xnumdd + xdendd
         gbremda = xnumda + xdenda
         gbremdz = xnumdz + xdendz



   !..equation 5.1
         dum    = 0.5738d0*input% zbar*input% ye*t8% t86*input% den
         dumdt  = 0.5738d0*input% zbar*input% ye*6.0d0*t8% t85*input% den*1.0d-8
         dumdd  = 0.5738d0*input% zbar*input% ye*t8% t86
         dumda  = -dum*input% abari
         dumdz  = 0.5738d0*2.0d0*input% ye*t8% t86*input% den

         z       = tfac4*fbrem - tfac5*gbrem
         sbrem   = dum * z
         sbremdt = dumdt*z + dum*(tfac4*fbremdt - tfac5*gbremdt) 
         sbremdd = dumdd*z + dum*(tfac4*fbremdd - tfac5*gbremdd) 
         sbremda = dumda*z + dum*(tfac4*fbremda - tfac5*gbremda) 
         sbremdz = dumdz*z + dum*(tfac4*fbremdz - tfac5*gbremdz) 

         
      end subroutine brem_neu_weak_degen
      
      
      subroutine brem_neu_liquid_metal(sbrem,sbremdt,sbremdd,sbremda,sbremdz,t8, input)
         type(t8s), intent(in) :: t8
         real(dp), intent(out) :: sbrem,sbremdt,sbremdd,sbremda,sbremdz
         type(inputs), intent(in) :: input

         real(dp) :: a0,a1,a2,c00,c01,c02,c03, &
         z, &
         dum,dumdt,dumdd,dumda,dumdz

   ! brem

         real(dp) :: u,gm1,gm2,gm13,gm23,gm43,gm53,v,w,fb,gt,gb, &
         fliq,fliqdt,fliqdd,fliqda,fliqdz,  &
         gliq,gliqdt,gliqdd,gliqda,gliqdz 

         real(dp) :: cos1,cos2,cos3,cos4,cos5,sin1,sin2, &
         sin3,sin4,sin5

         real(dp) :: ft

         sbrem=0d0
         sbremdt=0d0
         sbremdd=0d0
         sbremda=0d0
         sbremdz=0d0

         
   !..liquid metal with c12 parameters (not too different for other elements)
   !..equation 5.18 and 5.16
         u     = fac3 * (input% logden - 3.0d0)
         a0    = iln10*fac3 * input% deni

   !..compute the expensive trig functions of equation 5.21 only once
         cos1 = cos(u)
         cos2 = cos(2.0d0*u)
         cos3 = cos(3.0d0*u)
         cos4 = cos(4.0d0*u)
         cos5 = cos(5.0d0*u)

         sin1 = sin(u)
         sin2 = sin(2.0d0*u)
         sin3 = sin(3.0d0*u)
         sin4 = sin(4.0d0*u)
         sin5 = sin(5.0d0*u)

   !..equation 5.21
         fb =  0.5d0 * 0.17946d0  + 0.00945d0*u + 0.34529d0    &
               - 0.05821d0*cos1 - 0.04969d0*sin1 &
               - 0.01089d0*cos2 - 0.01584d0*sin2 &
               - 0.01147d0*cos3 - 0.00504d0*sin3 &
               - 0.00656d0*cos4 - 0.00281d0*sin4 &
               - 0.00519d0*cos5 

         c00 =  a0*(0.00945d0  &
               + 0.05821d0*sin1       - 0.04969d0*cos1 &
               + 0.01089d0*sin2*2.0d0 - 0.01584d0*cos2*2.0d0 &
               + 0.01147d0*sin3*3.0d0 - 0.00504d0*cos3*3.0d0 &
               + 0.00656d0*sin4*4.0d0 - 0.00281d0*cos4*4.0d0 &
               + 0.00519d0*sin5*5.0d0) 

      
   !..equation 5.22
         ft =  0.5d0 * 0.06781d0 - 0.02342d0*u + 0.24819d0 &
               - 0.00944d0*cos1 - 0.02213d0*sin1 &
               - 0.01289d0*cos2 - 0.01136d0*sin2 &
               - 0.00589d0*cos3 - 0.00467d0*sin3 &
               - 0.00404d0*cos4 - 0.00131d0*sin4 &
               - 0.00330d0*cos5 

         c01 = a0*(-0.02342d0   &
               + 0.00944d0*sin1       - 0.02213d0*cos1 &
               + 0.01289d0*sin2*2.0d0 - 0.01136d0*cos2*2.0d0 &
               + 0.00589d0*sin3*3.0d0 - 0.00467d0*cos3*3.0d0 &
               + 0.00404d0*sin4*4.0d0 - 0.00131d0*cos4*4.0d0 &
               + 0.00330d0*sin5*5.0d0) 


   !..equation 5.23
         gb =  0.5d0 * 0.00766d0 - 0.01259d0*u + 0.07917d0 &
               - 0.00710d0*cos1 + 0.02300d0*sin1 &
               - 0.00028d0*cos2 - 0.01078d0*sin2 &
               + 0.00232d0*cos3 + 0.00118d0*sin3 &
               + 0.00044d0*cos4 - 0.00089d0*sin4 &
               + 0.00158d0*cos5

         c02 = a0*(-0.01259d0 &
               + 0.00710d0*sin1       + 0.02300d0*cos1 &
               + 0.00028d0*sin2*2.0d0 - 0.01078d0*cos2*2.0d0 &
               - 0.00232d0*sin3*3.0d0 + 0.00118d0*cos3*3.0d0 &
               - 0.00044d0*sin4*4.0d0 - 0.00089d0*cos4*4.0d0 &
               - 0.00158d0*sin5*5.0d0)


   !..equation 5.24
         gt =  -0.5d0 * 0.00769d0  - 0.00829d0*u + 0.05211d0 &
               + 0.00356d0*cos1 + 0.01052d0*sin1 &
               - 0.00184d0*cos2 - 0.00354d0*sin2 &
               + 0.00146d0*cos3 - 0.00014d0*sin3 &
               + 0.00031d0*cos4 - 0.00018d0*sin4 &
               + 0.00069d0*cos5 

         c03 = a0*(-0.00829d0 &
               - 0.00356d0*sin1       + 0.01052d0*cos1 &
               + 0.00184d0*sin2*2.0d0 - 0.00354d0*cos2*2.0d0 &
               - 0.00146d0*sin3*3.0d0 - 0.00014d0*cos3*3.0d0 &
               - 0.00031d0*sin4*4.0d0 - 0.00018d0*cos4*4.0d0 &
               - 0.00069d0*sin5*5.0d0) 


         dum   = 2.275d-1 * input% zbar * input% zbar*t8% t8m1 * pow(input% den6*input% abari, one_third)
         dumdt = -dum*input% tempi
         dumdd = one_third*dum * input% deni
         dumda = -one_third*dum*input% abari
         dumdz = 2.0d0*dum*input% zbari
      
         gm1   = 1.0d0/dum
         gm2   = gm1*gm1
         gm13  = pow(gm1,one_third)
         gm23  = gm13 * gm13
         gm43  = gm13*gm1
         gm53  = gm23*gm1


   !..equation 5.25 and 5.26
         v  = -0.05483d0 - 0.01946d0*gm13 + 1.86310d0*gm23 - 0.78873d0*gm1
         a0 = one_third*0.01946d0*gm43 - two_thirds*1.86310d0*gm53 + 0.78873d0*gm2

         w  = -0.06711d0 + 0.06859d0*gm13 + 1.74360d0*gm23 - 0.74498d0*gm1
         a1 = -one_third*0.06859d0*gm43 - two_thirds*1.74360d0*gm53 + 0.74498d0*gm2


   !..equation 5.19 and 5.20
         fliq   = v*fb + (1.0d0 - v)*ft
         fliqdt = a0*dumdt*(fb - ft)
         fliqdd = a0*dumdd*(fb - ft) + v*c00 + (1.0d0 - v)*c01 
         fliqda = a0*dumda*(fb - ft)
         fliqdz = a0*dumdz*(fb - ft)

         gliq   = w*gb + (1.0d0 - w)*gt
         gliqdt = a1*dumdt*(gb - gt)
         gliqdd = a1*dumdd*(gb - gt) + w*c02 + (1.0d0 - w)*c03
         gliqda = a1*dumda*(gb - gt)
         gliqdz = a1*dumdz*(gb - gt)


   !..equation 5.17
         dum    = 0.5738d0*input% zbar*input% ye*t8% t86*input% den
         dumdt  = 0.5738d0*input% zbar*input% ye*6.0d0*t8% t85*input% den*1.0d-8
         dumdd  = 0.5738d0*input% zbar*input% ye*t8% t86
         dumda  = -dum*input% abari
         dumdz  = 0.5738d0*2.0d0*input% ye*t8% t86*input% den

         z       = tfac4*fliq - tfac5*gliq
         sbrem   = dum * z
         sbremdt = dumdt*z + dum*(tfac4*fliqdt - tfac5*gliqdt) 
         sbremdd = dumdd*z + dum*(tfac4*fliqdd - tfac5*gliqdd) 
         sbremda = dumda*z + dum*(tfac4*fliqda - tfac5*gliqda) 
         sbremdz = dumdz*z + dum*(tfac4*fliqdz - tfac5*gliqdz) 
      
      end subroutine brem_neu_liquid_metal
      
      
      subroutine brem_neu(sbrem,sbremdt,sbremdd,sbremda,sbremdz, input)
         real(dp), intent(out) :: sbrem,sbremdt,sbremdd,sbremda,sbremdz
         type(inputs), intent(in) :: input

         real(dp) :: a0,a1,a2,tfermi

         ! brem

         type(t8s) :: t8

         
   !..bremsstrahlung neutrino section 
   !..for reactions like e- + (z,a) => e- + (z,a) + nu + nubar
   !..                   n  + n     => n + n + nu + nubar
   !..                   n  + p     => n + p + nu + nubar

         real(dp) :: alfa, beta, sb, sbdt, sbdd, sbda, sbdz
         real(dp) :: sb2, sbdt2, sbdd2, sbda2, sbdz2
         real(dp), parameter :: tflo = 0.05d0, tfhi = 0.3d0
         real(dp) :: dtfracdt, dtfracdd,  dtfracda,  dtfracdz
         real(dp) :: dalfadt, dalfadd,  dalfada,  dalfadz
         real(dp) :: dbetadt, dbetadd,  dbetada,  dbetadz
         real(dp) :: dtfermidd,  dtfermida,  dtfermidz
         real(dp) :: A, B, C, D, U, dtfermidu, dtfracdtfermi, dtf, tfrac
         real(dp) :: dudv, dudd, duda, dudz

         sbrem=0d0
         sbremdt=0d0
         sbremdd=0d0
         sbremda=0d0
         sbremdz=0d0

   !..equation 4.3

         t8% t8     = input% temp * 1.0d-8
         t8% t812   = sqrt(t8% t8)
         t8% t832   = t8% t8 * t8% t812
         t8% t82    = t8% t8*t8% t8
         t8% t83    = t8% t82*t8% t8 
         t8% t85    = t8% t82*t8% t83
         t8% t86    = t8% t85*t8% t8
         t8% t8m1   = 1.0d0/t8% t8
         t8% t8m2   = t8% t8m1*t8% t8m1
         t8% t8m3   = t8% t8m2*t8% t8m1
         t8% t8m5   = t8% t8m3*t8% t8m2
         t8% t8m6   = t8% t8m5*t8% t8m1


         A = 5.9302d9
         B = 1.d0
         C = 1.018d0
         D = 1.0d0
         
         U = pow(input% den6*input% ye,two_thirds)
         tfermi = A * (sqrt(U) - D)
         
         if (input% temp .ge. tfhi * tfermi) then
         
            call brem_neu_weak_degen(sbrem,sbremdt,sbremdd,sbremda,sbremdz,t8, input)

         else if (input% temp .le. tflo * tfermi) then

            call brem_neu_liquid_metal(sbrem,sbremdt,sbremdd,sbremda,sbremdz,t8, input)
         
         else ! blend
         
            call brem_neu_weak_degen(sbrem,sbremdt,sbremdd,sbremda,sbremdz,t8, input)
            sb   = sbrem
            sbdt = sbremdt
            sbdd = sbremdd
            sbda = sbremda
            sbdz = sbremdz
         
            call brem_neu_liquid_metal(sbrem,sbremdt,sbremdd,sbremda,sbremdz,t8, input)
            sb2   = sbrem
            sbdt2 = sbremdt
            sbdd2 = sbremdd
            sbda2 = sbremda
            sbdz2 = sbremdz
            
            dtf = tfhi - tflo
            tfrac = (input% temp / tfermi - tflo) / dtf
            alfa = 0.5d0 * (1d0 - cospi(tfrac))
            beta = 1d0 - alfa
            

            dtfermidu = (1d0/2d0) * A * pow(U,-1d0/2d0)
            dtfracdtfermi = -input% temp/(tfermi * tfermi * dtf )
            
            ! v = den6* ye  = den *10**-6 * ye
            dudv = two_thirds * pow(input% den6 * input% ye,-1d0/3d0)
            
            dudd = dudv * 1d-6 * input% ye
            duda = dudv * input% den6 * input% ye * input% abari * (-1d0)
            dudz = dudv * input% den6 * input% abari
            
   
            dtfermidd = dtfermidu * dudd
            dtfermida = dtfermidu * duda
            dtfermidz = dtfermidu * dudz
            
            dtfracdt = 1.0d0/(tfermi * dtf) 
            
            dtfracdd = dtfermidd * dtfracdtfermi
            dtfracda = dtfermida * dtfracdtfermi
            dtfracdz = dtfermidz * dtfracdtfermi
            
            dalfadt = dtfracdt * 0.5d0 * pi * sinpi(tfrac)
            dalfadd = dtfracdd * 0.5d0 * pi * sinpi(tfrac) 
            dalfada = dtfracda * 0.5d0 * pi * sinpi(tfrac) 
            dalfadz = dtfracdz * 0.5d0 * pi * sinpi(tfrac) 
            
            dbetadt = -dalfadt
            dbetadd = -dalfadd
            dbetada = -dalfada
            dbetadz = -dalfadz
            
            sbrem   = alfa * sb   + beta * sb2
            sbremdt = alfa * sbdt + beta * sbdt2 + dalfadt * sb + dbetadt * sb2
            sbremdd = alfa * sbdd + beta * sbdd2 + dalfadd * sb + dbetadd * sb2
            sbremda = alfa * sbda + beta * sbda2 + dalfada * sb + dbetada * sb2
            sbremdz = alfa * sbdz + beta * sbdz2 + dalfadz * sb + dbetadz * sb2

   
         end if

      
      end subroutine brem_neu
      
      
      subroutine reco_neu(sreco,srecodt,srecodd,srecoda,srecodz, input)
         real(dp), intent(out) :: sreco,srecodt,srecodd,srecoda,srecodz
         type(inputs), intent(in) :: input

         real(dp) :: a0,a1,a2,a3,c00,c01, &
         dd00,dd01, &
         b,c,d,f1,f2,f3,z, &
         dum,dumdt,dumdd,dumda,dumdz, &
         gum,gumdt,gumdd,gumda,gumdz

         real(dp) :: zeta,zetadt,zetadd,zetada,zetadz
         real(dp) :: xnum,xnumdt,xnumdd,xnumda,xnumdz

         real(dp) :: nu,nudt,nudd,nuda,nudz, &
         nu2,nu3,bigj,bigjdt,bigjdd,bigjda,bigjdz

         sreco=0d0
         srecodt=0d0
         srecodd=0d0
         srecoda=0d0
         srecodz=0d0
         
   !..recombination neutrino section
   !..for reactions like e- (continuum) => e- (bound) + nu_e + nubar_e

   !..equation 6.11 solved for nu
         xnum   = 1.10520d8 * input% den * input% ye /(input% temp*sqrt(input% temp))
         xnumdt = -1.50d0*xnum*input% tempi
         xnumdd = xnum * input% deni
         xnumda = -xnum*input% abari
         xnumdz = xnum*input% zbari

   !..the chemical potential
         nu   = ifermi12(xnum)

   !..a0 is d(nu)/d(xnum)
         a0 = 1.0d0/(0.5d0*zfermim12(nu))
         nudt = a0*xnumdt
         nudd = a0*xnumdd
         nuda = a0*xnumda
         nudz = a0*xnumdz

         nu2  = nu * nu
         nu3  = nu2 * nu

   !..table 12
         if (nu .ge. -20.0d0  .and. nu .lt. 0.0d0) then
            a1 = 1.51d-2
            a2 = 2.42d-1
            a3 = 1.21d0
            b  = 3.71d-2
            c  = 9.06d-1
            d  = 9.28d-1
            f1 = 0.0d0
            f2 = 0.0d0
            f3 = 0.0d0
         else if (nu .ge. 0.0d0  .and. nu .le. 10.0d0) then
            a1 = 1.23d-2
            a2 = 2.66d-1
            a3 = 1.30d0
            b  = 1.17d-1
            c  = 8.97d-1
            d  = 1.77d-1
            f1 = -1.20d-2
            f2 = 2.29d-2
            f3 = -1.04d-3
         end if


   !..equation 6.7, 6.13 and 6.14
         if (nu .ge. -20.0d0  .and.  nu .le. 10.0d0) then

            zeta   = 1.579d5*input% zbar*input% zbar*input% tempi
            zetadt = -zeta*input% tempi
            zetadd = 0.0d0
            zetada = 0.0d0
            zetadz = 2.0d0*zeta*input% zbari

            c00    = 1.0d0/(1.0d0 + f1*nu + f2*nu2 + f3*nu3)  
            c01    = f1 + f2*2.0d0*nu + f3*3.0d0*nu2
            dum    = zeta*c00
            dumdt  = zetadt*c00 -1d0 * c00 * c00 * zeta*c01*nudt
            dumdd  = -1d0 * c00 * c00 * zeta*c01*nudd
            dumda  = -1d0 * c00 * c00 * zeta*c01*nuda
            dumdz  = zetadz*c00 -1d0 * c00 *c00 * zeta*c01*nudz
         
            z      = 1.0d0/dum
            dd00   = pow(dum,-2.25d0) 
            dd01   = pow(dum,-4.55d0)
            c00    = a1*z + a2*dd00 + a3*dd01
            c01    = -(a1*z + 2.25d0*a2*dd00 + 4.55d0*a3*dd01)*z
         
            z      = exp(c*nu)  
            dd00   = b*z*(1.0d0 + d*dum)        
            gum    = 1.0d0 + dd00
            gumdt  = dd00*c*nudt + b*z*d*dumdt  
            gumdd  = dd00*c*nudd + b*z*d*dumdd  
            gumda  = dd00*c*nuda + b*z*d*dumda  
            gumdz  = dd00*c*nudz + b*z*d*dumdz  

            z   = exp(nu)  
            a1  = 1.0d0/gum

            bigj   = c00 * z * a1
            bigjdt = c01*dumdt*z*a1 + c00*z*nudt*a1 - c00*z*a1*a1 * gumdt
            bigjdd = c01*dumdd*z*a1 + c00*z*nudd*a1 - c00*z*a1*a1 * gumdd
            bigjda = c01*dumda*z*a1 + c00*z*nuda*a1 - c00*z*a1*a1 * gumda
            bigjdz = c01*dumdz*z*a1 + c00*z*nudz*a1 - c00*z*a1*a1 * gumdz

            !..equation 6.5
            z     = exp(zeta + nu)
            dum   = 1.0d0 + z
            a1    = 1.0d0/dum
            a2    = 1.0d0/bigj

            sreco   = tfac6 * 2.649d-18 * input% ye * pow(input% zbar,13) * input% den * bigj*a1
            srecodt = sreco*(bigjdt*a2 - z*(zetadt + nudt)*a1)
            srecodd = sreco*(1.0d0 * input% deni + bigjdd*a2 - z*(zetadd + nudd)*a1)
            srecoda = sreco*(-1.0d0*input% abari + bigjda*a2 - z*(zetada+nuda)*a1)
            srecodz = sreco*(14.0d0*input% zbari + bigjdz*a2 - z*(zetadz+nudz)*a1)

         end if 

      
      end subroutine reco_neu
      
      subroutine plas_neu(splas,splasdt,splasdd,splasda,splasdz, input)
         real(dp), intent(out) :: splas,splasdt,splasdd,splasda,splasdz
         type(inputs), intent(in) :: input

         real(dp) :: a0,a1,a2,a3,b1,b2,c00,c01,c02,c03,c04,xlnt,cc, &
                     c,d,f1,&
                     dumdt,dumdd,dumda,dumdz, gum

         real(dp) :: xnum,xnumdt,xnumdd,xnumda,xnumdz, &
         xden,xdendt,xdendd,xdenda,xdendz

         real(dp) :: gl,gl2,gl2dt,gl2dd,gl2da,gl2dz,gl12,gl32,gl72,gl6, &
         ft,ftdt,ftdd,ftda,ftdz,fl,fldt,fldd,flda,fldz, &
         fxy,fxydt,fxydd,fxyda,fxydz

         splas=0d0
         splasdt=0d0
         splasdd=0d0
         splasda=0d0
         splasdz=0d0
         
   !..plasma neutrino section 
   !..for collective reactions like gamma_plasmon => nu_e + nubar_e
   !..equation 4.6

         a1   = 1.019d-6*input% rm
         a2   = pow(a1,two_thirds)
         a3   = two_thirds*a2/a1

         b1   =  sqrt(1.0d0 + a2)
         b2   = 1.0d0/b1
   
         c00  = 1.0d0/(input% temp*input% temp*b1)

         gl2   = 1.1095d11 * input% rm * c00

         gl2dt = -2.0d0*gl2*input% tempi
         d     = input% rm*c00*b2*0.5d0*b2*a3*1.019d-6
         gl2dd = 1.1095d11 * (input% rmdd*c00  - d*input% rmdd)
         gl2da = 1.1095d11 * (input% rmda*c00  - d*input% rmda)
         gl2dz = 1.1095d11 * (input% rmdz*c00  - d*input% rmdz)
         

         gl    = sqrt(gl2)
         gl12  = sqrt(gl)
         gl32  = gl * gl12
         gl72  = gl2 * gl32
         gl6   = gl2 * gl2 * gl2


   !..equation 4.7
         ft   = 2.4d0 + 0.6d0*gl12 + 0.51d0*gl + 1.25d0*gl32
         gum  = 1.0d0/gl2
         a1   =(0.25d0*0.6d0*gl12 +0.5d0*0.51d0*gl +0.75d0*1.25d0*gl32)*gum
         ftdt = a1*gl2dt
         ftdd = a1*gl2dd
         ftda = a1*gl2da
         ftdz = a1*gl2dz


   !..equation 4.8
         a1   = 8.6d0*gl2 + 1.35d0*gl72
         a2   = 8.6d0 + 1.75d0*1.35d0*gl72*gum

         b1   = 225.0d0 - 17.0d0*gl + gl2
         b2   = -0.5d0*17.0d0*gl*gum + 1.0d0

         c    = 1.0d0/b1
         fl   = a1*c

         d    = (a2 - fl*b2)*c       
         fldt = d*gl2dt
         fldd = d*gl2dd
         flda = d*gl2da
         fldz = d*gl2dz
      

   !..equation 4.9 and 4.10
         cc   = log10(2.0d0*input% rm)
         xlnt = input% logtemp

         xnum   = one_sixth * (17.5d0 + cc - 3.0d0*xlnt)
         xnumdt = -iln10*0.5d0*input% tempi
         a2     = iln10*one_sixth*input% rmi
         xnumdd = a2*input% rmdd
         xnumda = a2*input% rmda 
         xnumdz = a2*input% rmdz 

         xden   = one_sixth * (-24.5d0 + cc + 3.0d0*xlnt)
         xdendt = iln10*0.5d0*input% tempi
         xdendd = a2*input% rmdd
         xdenda = a2*input% rmda 
         xdendz = a2*input% rmdz 


   !..equation 4.11
         if (abs(xnum) .gt. 0.7d0  .or.  xden .lt. 0.0d0) then
            fxy   = 1.0d0
            fxydt = 0.0d0
            fxydd = 0.0d0
            fxydz = 0.0d0
            fxyda = 0.0d0

         else 

            a1  = 0.39d0 - 1.25d0*xnum - 0.35d0*sin(4.5d0*xnum)
            a2  = -1.25d0 - 4.5d0*0.35d0*cos(4.5d0*xnum)

            b1  = 0.3d0 * exp(-1.0d0*pow((4.5d0*xnum + 0.9d0),2))
            b2  = -b1*2.0d0*(4.5d0*xnum + 0.9d0)*4.5d0

            c   = min(0.0d0, xden - 1.6d0 + 1.25d0*xnum)
         if (c .eq. 0.0d0) then
            dumdt = 0.0d0
            dumdd = 0.0d0
            dumda = 0.0d0
            dumdz = 0.0d0
         else
            dumdt = xdendt + 1.25d0*xnumdt
            dumdd = xdendd + 1.25d0*xnumdd
            dumda = xdenda + 1.25d0*xnumda
            dumdz = xdendz + 1.25d0*xnumdz
         end if

            d   = 0.57d0 - 0.25d0*xnum
            a3  = c/d
            c00 = exp(-1.0d0*a3*a3)

            f1  = -c00*2.0d0*a3/d
            c01 = f1*(dumdt + a3*0.25d0*xnumdt)
            c02 = f1*(dumdd + a3*0.25d0*xnumdd)
            c03 = f1*(dumda + a3*0.25d0*xnumda)
            c04 = f1*(dumdz + a3*0.25d0*xnumdz)

            fxy   = 1.05d0 + (a1 - b1)*c00
            fxydt = (a2*xnumdt -  b2*xnumdt)*c00 + (a1-b1)*c01
            fxydd = (a2*xnumdd -  b2*xnumdd)*c00 + (a1-b1)*c02
            fxyda = (a2*xnumda -  b2*xnumda)*c00 + (a1-b1)*c03
            fxydz = (a2*xnumdz -  b2*xnumdz)*c00 + (a1-b1)*c04

         end if



   !..equation 4.1 and 4.5
         splas   = (ft + fl) * fxy
         splasdt = (ftdt + fldt)*fxy + (ft+fl)*fxydt
         splasdd = (ftdd + fldd)*fxy + (ft+fl)*fxydd
         splasda = (ftda + flda)*fxy + (ft+fl)*fxyda
         splasdz = (ftdz + fldz)*fxy + (ft+fl)*fxydz

         a2      = exp(-gl)
         a3      = -0.5d0*a2*gl*gum

         a1      = splas
         splas   = a2*a1
         splasdt = a2*splasdt + a3*gl2dt*a1
         splasdd = a2*splasdd + a3*gl2dd*a1
         splasda = a2*splasda + a3*gl2da*a1
         splasdz = a2*splasdz + a3*gl2dz*a1

         a2      = gl6
         a3      = 3.0d0*gl6*gum

         a1      = splas
         splas   = a2*a1
         splasdt = a2*splasdt + a3*gl2dt*a1
         splasdd = a2*splasdd + a3*gl2dd*a1
         splasda = a2*splasda + a3*gl2da*a1
         splasdz = a2*splasdz + a3*gl2dz*a1


         a2      = 0.93153d0 * 3.0d21 * input% xl9
         a3      = 0.93153d0 * 3.0d21 * 9.0d0*input% xl8*input% xldt

         a1      = splas
         splas   = a2*a1
         splasdt = a2*splasdt + a3*a1
         splasdd = a2*splasdd 
         splasda = a2*splasda 
         splasdz = a2*splasdz 

      
      end subroutine plas_neu



      subroutine pair_neu(spair,spairdt,spairdd,spairda,spairdz, input)
         real(dp), intent(out) :: spair,spairdt,spairdd,spairda,spairdz
         type(inputs), intent(in) :: input

         real(dp) :: a0,a1,a2,a3,b1,b2,c,d, gl,gldt
         
         real(dp) :: xnum,xnumdt,xnumdd,xnumda,xnumdz, &
         xden,xdendt,xdendd,xdenda,xdendz

         real(dp) :: fpair,fpairdt,fpairdd,fpairda,fpairdz, &
         qpair,qpairdt,qpairdd,qpairda,qpairdz

   !..pair neutrino section
   !..for reactions like e+ + e- => nu_e + nubar_e 

   !..equation 2.8 
         gl   = 1.0d0 - 13.04d0*input% xl2 +133.5d0*input% xl4 +1534.0d0*input% xl6 +918.6d0*input% xl8
         gldt = input% xldt*(-26.08d0*input% xl +534.0d0*input% xl3 +9204.0d0*input% xl5 +7348.8d0*input% xl7)

   !..equation 2.7

         a1     = 6.002d19 + 2.084d20*input% zeta + 1.872d21*input% zeta2
         a2     = 2.084d20 + 2.0d0*1.872d21*input% zeta

         if (input% t9 .lt. 10.0d0) then
            b1     = exp(-5.5924d0*input% zeta)
            b2     = -b1*5.5924d0
         else
            b1     = exp(-4.9924d0*input% zeta)
            b2     = -b1*4.9924d0
         end if
         
         xnum   = a1 * b1
         c      = a2*b1 + a1*b2
         xnumdt = c*input% zetadt
         xnumdd = c*input% zetadd
         xnumda = c*input% zetada
         xnumdz = c*input% zetadz

         if (input% t9 .lt. 10.0d0) then
            a1   = 9.383d-1*input% xlm1 - 4.141d-1*input% xlm2 + 5.829d-2*input% xlm3
            a2   = -9.383d-1*input% xlm2 + 2.0d0*4.141d-1*input% xlm3 - 3.0d0*5.829d-2*input% xlm4
         else
            a1   = 1.2383d0*input% xlm1 - 8.141d-1*input% xlm2 
            a2   = -1.2383d0*input% xlm2 + 2.0d0*8.141d-1*input% xlm3 
         end if

         b1   = 3.0d0*input% zeta2

         xden   = input% zeta3 + a1
         xdendt = b1*input% zetadt + a2*input% xldt
         xdendd = b1*input% zetadd
         xdenda = b1*input% zetada
         xdendz = b1*input% zetadz

         a1      = 1.0d0/xden
         fpair   = xnum*a1
         fpairdt = (xnumdt - fpair*xdendt)*a1
         fpairdd = (xnumdd - fpair*xdendd)*a1
         fpairda = (xnumda - fpair*xdenda)*a1
         fpairdz = (xnumdz - fpair*xdendz)*a1


   !..equation 2.6
         a1     = 10.7480d0*input% xl2 + 0.3967d0*input% xlp5 + 1.005d0
         a2     = input% xldt*(2.0d0*10.7480d0*input% xl + 0.5d0*0.3967d0*input% xlmp5) 
         xnum   = 1.0d0/a1
         xnumdt = -xnum*xnum*a2

         a1     = 7.692d7*input% xl3 + 9.715d6*input% xlp5
         a2     = input% xldt*(3.0d0*7.692d7*input% xl2 + 0.5d0*9.715d6*input% xlmp5)

         c      = 1.0d0/a1
         b1     = 1.0d0 + input% rm*c

         xden   = pow(b1,-0.3d0)

         d      = -0.3d0*xden/b1
         xdendt = -d*input% rm*c*c*a2
         xdendd = d*input% rmdd*c 
         xdenda = d*input% rmda*c 
         xdendz = d*input% rmdz*c 

         qpair   = xnum*xden
         qpairdt = xnumdt*xden + xnum*xdendt
         qpairdd = xnum*xdendd
         qpairda = xnum*xdenda
         qpairdz = xnum*xdendz

   !..equation 2.5
         a1    = exp(-2.0d0*input% xlm1)
         a2    = a1*2.0d0*input% xlm2*input% xldt

         spair   = a1*fpair
         spairdt = a2*fpair + a1*fpairdt
         spairdd = a1*fpairdd
         spairda = a1*fpairda
         spairdz = a1*fpairdz

         a1      = spair
         spair   = gl*a1
         spairdt = gl*spairdt + gldt*a1
         spairdd = gl*spairdd
         spairda = gl*spairda
         spairdz = gl*spairdz

         a1      = tfac4*(1.0d0 + tfac3 * qpair)
         a2      = tfac4*tfac3

         a3      = spair
         spair   = a1*a3
         spairdt = a1*spairdt + a2*qpairdt*a3
         spairdd = a1*spairdd + a2*qpairdd*a3
         spairda = a1*spairda + a2*qpairda*a3
         spairdz = a1*spairdz + a2*qpairdz*a3
      
      end subroutine pair_neu


      end module mod_neu
   


