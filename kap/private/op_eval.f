! ***********************************************************************
!
!   Copyright (C) 2013-2019  Haili Hu & The MESA Team
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
 
c FORTRAN 90 module for calculation of radiative accelerations,
c based on the Opacity Project (OP) code "OPserver".
c See CHANGES_HU for changes made to the original code.
c
c Haili Hu 2010
c
      module op_eval
      use op_load
      use const_def, only: dp
      use math_lib
      use kap_def, only: kap_test_partials, kap_test_partials_val, kap_test_partials_dval_dx
      
      logical, parameter :: dbg = .false.

      contains




c HH: Based on "op_ax.f"
c Input:   kk = number of elements to calculate g_rad for
c          iz1(kk) = charge of element to calculate g_rad for
c          nel = number of elements in mixture
c          izzp(nel) = charge of elements
c          fap(nel) = number fractions of elements
c          fac(nel) = scale factors for element opacity   
c          flux = local radiative flux (Lrad/4*pi*r^2)
c          fltp = log T
c          flrhop = log rho
c          screening   if true, use screening corrections
c Output: g1 = log kappa
c         gx1 = d(log kappa)/d(log T)
c         gy1 = d(log kappa)/d(log rho)
c         gp1(kk) = d(log kappa)/d(log xi) 
c         grl1(kk) = log grad
c         fx1(kk) = d(log grad)/d(log T) 
c         fy1(kk) = d(log grad)/d(log rho)
c         grlp1(kk) = d(log grad)/d(log xi)
c         meanZ(nel) = average ionic charge of elements
c         zetx1(nel) = d(meanZ)/d(log T) 
c         zety1(nel) = d(meanZ)/d(log rho)
c         ierr = 0 for correct use 
      subroutine eval_op_radacc(
     > kk, izk, nel, izzp, fap, fac, flux, fltp, flrhop, screening,
     : g1, grl1,
     > umesh, semesh, ff, ta, rs, ierr)
      use op_radacc
      use op_load, only: msh
      use op_common
      implicit none
      integer, intent(in) :: kk, nel
      integer, intent(in) :: izk(kk), izzp(nel)
      real(dp), intent(in) :: fap(nel), fac(nel)
      real(dp), intent(in) :: flux, fltp, flrhop
      logical, intent(in) :: screening
      real(dp), intent(out) :: g1
      real(dp), intent(inout) :: grl1(kk)
      real, pointer :: umesh(:), semesh(:), ff(:,:,:,:), ta(:,:,:,:), rs(:,:,:)
         ! umesh(nptot)
         ! semesh(nptot)
         ! ff(nptot, ipe, 4, 4)
         ! ta(nptot, nrad, 4, 4), 
         ! rs(nptot, 4, 4)
      integer,intent(out) :: ierr
c local variables      
      integer :: n, i, k2, i3, ntot, jhmin, jhmax
      integer :: ih(4), jh(4), ilab(4), kzz(nrad), nkz(ipe), izz(ipe), iz1(nrad)
      real :: const, gx, gy, flt, flrho, flmu, dscat, dv, xi, flne,
     : epa, eta, ux, uy, g
      real :: uf(0:100), rion(28,4,4), rossl(4,4), flr(4,4),
     : rr(28, ipe, 4, 4), fa(ipe), 
     : gaml(4, 4, nrad), f(nrad), am1(nrad), 
     : fmu1(nrad)
c
c  Initialisations
      include 'formats'
      
      ierr = 0      
      if(nel.le.0.or.nel.gt.ipe) then
         write(6,*)'OP - NUMBER OF ELEMENTS OUT OF RANGE:', nel
         ierr = 1
         return
      endif

      if (dbg) write(*,*) 'start eval_op_radacc'
      


c  Get i3 for mesh type q='m'
      i3=2      
c      
c HH: k2 loops over elements for which to calculate grad.
      do k2 = 1, kk
         do n = 1, ipe
            if(izk(k2).eq.kz(n)) then
               iz1(k2) = izk(k2)
               exit
            endif   
            if(n.eq.ipe) then
               write(6,*) 'OP - SELECTED ELEMENT CANNOT BE TREATED: Z = ', izk(k2)
               write(6,*) 'OP supports C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Cr, Mn, Fe, and Ni'
               ierr = 5
               return
            endif
         enddo
      enddo
c      
      outer: do i = 1, nel
         inner: do n = 1, ipe
            if(izzp(i).eq.kz(n)) then
               izz(i) = izzp(i)
               fa(i) = fap(i)
               if(fa(i).lt.0.0) then
                  write(6,*)'OP - NEGATIVE FRACTIONAL ABUNDANCE:',fa(i)
                  ierr = 7
                  return
               endif
               cycle outer
            endif
         enddo inner
         write(6,*)'OP - CHEM. ELEMENT CANNOT BE INCLUDED: Z = ', izzp(i)
         ierr = 8
         return
      enddo outer
c
c Calculate mean atomic weight (flmu) and 
c array kzz indicating elements for which to calculate g_rad
      if (dbg) write(*,*) 'call abund'
      call abund(nel, izz, kk, iz1, fa,     !input variables 
     :   kzz, flmu, am1, fmu1, nkz)         !output variables
c           
c  Other initialisations
c       dv = interval in frequency variable v
c       ntot=number of frequency points
c       umesh, values of u=(h*nu/k*T) on mesh points
c       semesh, values of 1-exp(-u) on mesh points
c       uf, dscat used in scattering correction
      if (dbg) write(*,*) 'call msh'
      call msh(dv, ntot, umesh, semesh, uf, dscat)  !output variables
c
c  Start loop on temperature-density points
c  flt=log10(T, K)
c  flrho=log10(rho, cgs)
c
      flt = fltp
      flrho = flrhop
c     Get temperature indices
c       Let ite(i) be temperature index used in mono files
c       Put ite(i)=2*ih(i)
c       Use ih(i), i=1 to 4
c       xi=interpolation variable
c       log10(T)=flt=0.025*(ite(1)+xi+3)
c       ilab(i) is temperature label
      if (dbg) write(*,*) 'call xindex'
      call xindex(flt, ilab, xi, ih, i3, ierr) 
      if (ierr /= 0) then
         write(*,*) "xindex errored in radacc"
         return
      endif
c      
c     Get density indices
c       Let jne(j) be density index used  in mono files
c       Put jne(j)=2*jh(j)
c       Use jh(j), j=1 to 4
c       Get extreme range for jh
      if (dbg) write(*,*) 'call jrange'
      call jrange(ih, jhmin, jhmax, i3)
c
c     Get electron density flne=log10(Ne) for specified mass density flrho
c       Also:  UY=0.25*[d log10(rho)]/[d log10(Ne)] 
c              epa=electrons per atom
      if (dbg) write(*,*) 'call findne'
      call findne(ilab, fa, nel, nkz, jhmin, jhmax, ih, flrho, flt, xi,
     : flne, flmu, flr, epa, uy, i3, ierr)
      if (ierr .ne. 0 ) then
         write(*,*) "findne encountered error in radacc"
         return
      endif
c
c     Get density indices jh(j), j=1 to 4,
c       Interpolation variable eta
c       log10(Ne)=flne=0.25*(jne(1)+eta+3)
      if (dbg) write(*,*) 'call yindex'
      call yindex(jhmin, jhmax, flne, jh, i3, eta)
c        
c     Get ux=0.025*[d log10(rho)]/[d log10(T)]
      if (dbg) write(*,*) 'call findux'
      call findux(flr, xi, eta, ux)
            
c    rossl(i,j)=log10(Rosseland mean) on mesh points (i,j)
c     Get new mono opacities, ff(n,k,i,j)
      if (dbg) write(*,*) 'call rd'
      call rd(i3, kk, kzz, nel, nkz, izz, ilab, jh, ntot, umesh,
     : semesh, ff, rr, ta, fac)

c     Get rs = weighted sum of monochromatic opacity cross sections
      if (dbg) write(*,*) 'call mix'
      call mix(kk, kzz, ntot, nel, fa, ff, rr, rs, rion)
c
c     Screening corrections      
      if (screening) then      
c        Get Boercker scattering correction
         if (dbg) write(*,*) 'call scatt'
         call scatt(ih, jh, rion, uf, rs, umesh, semesh, dscat, ntot, epa, ierr)
         if (ierr .ne. 0 ) return
c        Get correction for Debye screening
         if (dbg) write(*,*) 'call screen1'
         call screen1(ih, jh, rion, umesh, ntot, epa, rs)
      endif               
c      
c     Get rossl, array of log10(Rosseland mean in cgs)    
      if (dbg) write(*,*) 'call ross'
      call ross(kk, flmu, fmu1, dv, ntot, rs, rossl, gaml, ta)
c
c     Interpolate to required flt, flrho
c     g=log10(ross, cgs)
      if (dbg) write(*,*) 'call interp'
      call interp(nel, kk, rossl, gaml, xi, eta, g, i3, f)
      if (dbg) write(*,*) 'done interp'
c      
c Write grad in terms of local radiative flux instead of (Teff, r/R*):
      const = 13.30295d0 + log10(flux) ! = -log10(c) - log10(amu) + log10(flux)
      do k2 = 1, kk 
         grl1(k2) = const + flmu - log10(dble(am1(k2))) + f(k2) + g    ! log g_rad 
      enddo   !k2
c      
      g1 = g                ! log kappa

c
      if (dbg) write(*,*) 'done eval_op_radacc'
      return
c
      end subroutine eval_op_radacc
c***********************************************************************
c HH: Based on "op_mx.f", opacity calculations to be used for stellar evolution calculations 
c Input:   nel = number of elements in mixture
c          izzp(nel) = charge of elements
c          fap(nel) = number fractions of elements
c          fac(nel) = scale factors for element opacity   
c          fltp = log (temperature)
c          flrhop = log (mass density) 
c          screening   if true, use screening corrections
c Output: g1 = log kappa
c         gx1 = d(log kappa)/d(log T)
c         gy1 = d(log kappa)/d(log rho)
c         ierr = 0 for correct use 
      subroutine eval_op_ev(
     >         nel, izzp, fap, fac, fltp, flrhop, screening, g1, gx1, gy1, 
     >         umesh, semesh, ff, rs, ierr)
      use op_ev
      use op_load, only: msh
      use op_common
      implicit none
      integer, intent(in) :: nel
      integer, intent(in) :: izzp(nel)
      real(dp), intent(in) :: fap(nel), fac(nel)
      real(dp), intent(in) :: fltp, flrhop
      logical, intent(in) :: screening
      real(dp), intent(inout) :: g1, gx1, gy1
      real, pointer :: umesh(:), semesh(:), ff(:,:,:,:), rs(:,:,:)
         ! umesh(nptot)
         ! semesh(nptot)
         ! ff(nptot, ipe, 4, 4)
         ! rs(nptot, 4, 4)
         ! s(nptot, nrad, 4, 4)
      integer,intent(out) :: ierr
c local variables      
      integer :: n, i, i3, jhmin, jhmax, ntot
      integer :: ih(4), jh(4), ilab(4), izz(ipe), nkz(ipe) 
      real :: flt, flrho, flmu, flne, dv, dscat, const, gx, gy, g,
     : eta, epa, xi, ux, uy
      real :: uf(0:100), rion(28, 1:4, 1:4), rossl(4, 4), flr(4, 4), 
     : fa(ipe), rr(28, ipe, 4, 4),
     : fmu1(nrad)
c
c  Initialisations
      ierr=0      
      if(nel.le.0.or.nel.gt.ipe) then
         write(6,*)'OP - NUMBER OF ELEMENTS OUT OF RANGE:',nel
         ierr=1
         return
      endif
c      
c  Get i3 for mesh type q='m'
      i3=2
c      
        outer: do i=1,nel
          inner: do n=1,ipe
            if(izzp(i).eq.kz(n)) then
              izz(i)=izzp(i)
              fa(i)=fap(i)
              if(fa(i).lt.0.d0) then
                write(6,*)'OP - NEGATIVE FRACTIONAL ABUNDANCE:',fa(i)
                ierr=7
                return
              endif
              cycle outer
            endif
          enddo inner
          write(6,*)'OP - CHEM. ELEMENT CANNOT BE INCLUDED: Z = ',
     +      izzp(i)
          ierr=8
          return
       enddo outer

c Calculate mean atomic weight (flmu) 
      call abund(nel, izz, fa, flmu, fmu1, nkz)
c          
c  Other initialisations
c       dv = interval in frequency variable v
c       ntot=number of frequency points
c       umesh, values of u=(h*nu/k*T) on mesh points
c       semesh, values of 1-exp(-u) on mesh points
c       uf, dscat used in scattering correction
      call msh(dv, ntot, umesh, semesh, uf, dscat)
c
c  Start loop on temperature-density points
c  flt=log10(T, K)
c  flrho=log10(rho, cgs)
      flt = fltp
      flrho = flrhop
c     Get temperature indices
c       Let ite(i) be temperature index used in mono files
c       Put ite(i)=2*ih(i)
c       Use ih(i), i=1 to 4
c       xi=interpolation variable
c       log10(T)=flt=0.025*(ite(1)+xi+3)
c       ilab(i) is temperature label
      call xindex(flt, ilab, xi, ih, i3, ierr) 
      if (ierr /= 0) then
         write(*,*) "eval_op_ev failed in xindex"
         return
      endif
c      
c     Get density indices
c       Let jne(j) be density index used  in mono files
c       Put jne(j)=2*jh(j)
c       Use jh(j), j=1 to 4
c       Get extreme range for jh
      call jrange(ih, jhmin, jhmax, i3)   
c
c     Get electron density flne=log10(Ne) for specified mass density flrho
c       Also:  UY=0.25*[d log10(rho)]/[d log10(Ne)] 
c              epa=electrons per atom
      call findne(ilab, fa, nel, nkz, jhmin, jhmax, ih, flrho, flt,
     + xi, flne, flmu, flr, epa, uy, i3, ierr)
      if (ierr /= 0) then
         write(*,*) "eval_op_ev failed in findne"
         return
      endif
c
c     Get density indices jh(j), j=1 to 4,
c       Interpolation variable eta
c       log10(Ne)=flne=0.25*(jne(1)+eta+3)
      call yindex(jhmin, jhmax, flne, jh, i3, eta)
c        
c     Get ux=0.025*[d log10(rho)]/[d log10(T)]
      call findux(flr, xi, eta, ux)
            
c    rossl(i,j)=log10(Rosseland mean) on mesh points (i,j)
c     Get new mono opacities, ff(n,k,i,j)
      call rd(nel, nkz, izz, ilab, jh, ntot, ff, rr, i3, umesh, fac)

c     Up-date mixture
      call mix(ntot, nel, fa, ff, rs, rr, rion)  
c
      if(screening) then 
c        Get Boercker scattering correction
         call scatt(ih, jh, rion, uf, rs, umesh, semesh, dscat, ntot, epa, ierr)
         if (ierr /= 0) then
            write(*,*) "scattering correction failed in op_eval_ev"
            return
         endif
c        Get correction for Debye screening
         call screen1(ih, jh, rion, umesh, ntot, epa, rs)
      endif     
c      
c     Get rossl, array of log10(Rosseland mean in cgs)
      call ross(flmu, fmu1, dv, ntot, rs, rossl)
c      
c     Interpolate to required flt, flrho
c     g=log10(ross, cgs)
      call interp(nel, rossl, xi, eta, 
     : g, i3, ux, uy, gx, gy)
c            
      g1 = g                ! log kappa
      gx1 = gx              ! dlogkappa/dt
      gy1 = gy              ! dlogkappa/drho   
c
      return
c
      end subroutine eval_op_ev
***********************************************************************

c     HH: Based on "op_mx.f", opacity calculations to be used for non-adiabatic pulsation calculations
c Special care is taken to ensure smoothness of opacity derivatives
c Input:   nel = number of elements in mixture
c          izzp(nel) = charge of elements
c          fap(nel) = number fractions of elements
c          fac(nel) = scale factors for element opacity   
c          fltp = log (temperature)
c          flrhop = log (mass density) 
c          screening   if true, use screening corrections
c Output: g1 = log kappa
c         gx1 = d(log kappa)/d(log T)
c         gy1 = d(log kappa)/d(log rho)
c         ierr = 0 for correct use 
      subroutine eval_alt_op(
     >         nel, izzp, fap, fac, fltp, flrhop, screening, g1, gx1, gy1,
     >         umesh, semesh, ff, rs, ierr)
      use op_osc
      use op_load, only: msh
      implicit none
      integer, intent(in) :: nel
      integer, intent(in) :: izzp(nel)
      real(dp), intent(in) :: fap(nel), fac(nel)
      real(dp), intent(in) :: fltp, flrhop
      logical, intent(in) :: screening
      real(dp), intent(out) :: g1, gx1, gy1
!      real(dp), intent(out) :: meanZ(nel)
      real, pointer :: umesh(:), semesh(:), ff(:,:,:,:), rs(:,:,:)
         ! umesh(nptot)
         ! semesh(nptot)
         ! ff(nptot, ipe, 0:5, 0:5)
         ! rs(nptot, 0:5, 0:5)
      integer,intent(out) :: ierr
c local variables      
      integer :: n, i, i3, jhmin, jhmax, ntot
      integer :: ih(0:5), jh(0:5), ilab(0:5), izz(ipe), nkz(ipe) 
      real :: flt, flrho, flmu, flne, dv, dscat, const, gx, gy, g,
     : eta, epa, xi, ux, uy
      real :: uf(0:100), rion(28, 0:5, 0:5), rossl(0:5, 0:5), flr(4, 4), 
     : fa(ipe), rr(28, ipe, 0:5, 0:5)
c
c  Initialisations
      ierr=0      
      if(nel.le.0.or.nel.gt.ipe) then
         write(6,*)'OP - NUMBER OF ELEMENTS OUT OF RANGE:',nel
         ierr=1
         return
      endif
c      
c  Get i3 for mesh type q='m'
      i3=2
c      
        outer: do i=1,nel
          inner: do n=1,ipe
            if(izzp(i).eq.kz(n)) then
              izz(i)=izzp(i)
              fa(i)=fap(i)
              if(fa(i).lt.0.0) then
                write(6,*)'OP - NEGATIVE FRACTIONAL ABUNDANCE:',fa(i)
                ierr=7
                return
              endif
              cycle outer
            endif
          enddo inner
          write(6,*)'OP - CHEM. ELEMENT CANNOT BE INCLUDED: Z = ',
     +      izzp(i)
          ierr=8
          return
       enddo outer

c Calculate mean atomic weight (flmu) 
      call abund(nel, izz, fa, flmu, nkz)
c          
c  Other initialisations
c       dv = interval in frequency variable v
c       ntot=number of frequency points
c       umesh, values of u=(h*nu/k*T) on mesh points
c       semesh, values of 1-exp(-u) on mesh points
c       uf, dscat used in scattering correction
      call msh(dv, ntot, umesh, semesh, uf, dscat)
c
c  Start loop on temperature-density points
c  flt=log10(T, K)
c  flrho=log10(rho, cgs)
      flt = fltp
      flrho = flrhop
c     Get temperature indices
c       Let ite(i) be temperature index used in mono files
c       Put ite(i)=2*ih(i)
c       Use ih(i), i=1 to 4
c       xi=interpolation variable
c       log10(T)=flt=0.025*(ite(1)+xi+3)
c       ilab(i) is temperature label
      call xindex(flt, ilab, xi, ih, i3, ierr) 
      if (ierr /= 0) return
c      
c     Get density indices
c       Let jne(j) be density index used  in mono files
c       Put jne(j)=2*jh(j)
c       Use jh(j), j=1 to 4
c       Get extreme range for jh
      call jrange(ih, jhmin, jhmax, i3)   
c
c     Get electron density flne=log10(Ne) for specified mass density flrho
c       Also:  UY=0.25*[d log10(rho)]/[d log10(Ne)] 
c              epa=electrons per atom
      call findne(ilab, fa, nel, nkz, jhmin, jhmax, ih, flrho, flt,
     + xi, flne, flmu, flr, epa, uy, i3, ierr)
      if (ierr /= 0) return
c
c     Get density indices jh(j), j=1 to 4,
c       Interpolation variable eta
c       log10(Ne)=flne=0.25*(jne(1)+eta+3)
      call yindex(jhmin, jhmax, flne, jh, i3, eta)
c        
c     Get ux=0.025*[d log10(rho)]/[d log10(T)]
      call findux(flr, xi, eta, ux)
            
c    rossl(i,j)=log10(Rosseland mean) on mesh points (i,j)
c     Get new mono opacities, ff(n,k,i,j)
      call rd(nel, nkz, izz, ilab, jh, ntot, ff, rr, i3, umesh, fac)

c     Up-date mixture
      call mix(ntot, nel, fa, ff, rs, rr, rion)  
c
      if(screening) then 
c        Get Boercker scattering correction
         call scatt(ih, jh, rion, uf, rs, umesh, semesh, dscat, ntot, epa, ierr)
         if (ierr /= 0) return
c        Get correction for Debye screening
         call screen1(ih, jh, rion, umesh, ntot, epa, rs)
      endif     
c      
c     Get rossl, array of log10(Rosseland mean in cgs)
      call ross(flmu, dv, ntot, rs, rossl)
c      
c     Interpolate to required flt, flrho
c     g=log10(ross, cgs)
      call interp(nel, rossl, xi, eta, g, i3, ux, uy, gx, gy)
c      
      g1 = g                ! log kappa
      gx1 = gx              ! dlogkappa/dt
      gy1 = gy              ! dlogkappa/drho   
c        
      return
c
      end subroutine eval_alt_op
c***********************************************************************


      end module op_eval
