! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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

      module screen5
      use rates_def
      use const_def, only: ln10, pi, two_13
      use math_lib
      
      implicit none

      contains
      
      subroutine screen5_init_AZ_info( &
               zs13, zhat, zhat2, lzav, aznut, zs13inv, a1, z1, a2, z2, ierr)
         !..compute and store things that only depend on reaction A's and Z's
         real(dp), intent(out) :: zs13, zhat, zhat2, lzav, aznut, zs13inv
         ! zs13 = (z1+z2)**(1./3.)
         ! zhat = combination of z1 and z2 raised to the 5/3 power
         ! zhat2 = combination of z1 and z2 raised to the 5/12 power
         ! lzav = log of effective charge
         ! aznut = combination of a1, z1, a2, z2 raised to 1/3 power
         ! zs13inv = 1 / zs13
         real(dp), intent(in) :: a1, z1, a2, z2
         integer, intent(out) :: ierr
         
         real(dp), parameter :: x13   = 1.0d0/3.0d0 
         real(dp), parameter :: x14   = 1.0d0/4.0d0
         real(dp), parameter :: x53   = 5.0d0/3.0d0
         real(dp), parameter :: x532  = 5.0d0/32.0d0
         real(dp), parameter :: x512  = 5.0d0/12.0d0
         
         ierr = 0
         
         if (z1 <= 0 .or. z2 <= 0) then
            zs13    = 0
            zs13inv = 0
            zhat    = 0
            zhat2   = 0
            lzav    = 0
            aznut   = 0
            return
         end if

         zs13    = pow(z1 + z2, x13)
         zs13inv = 1.0d0/zs13
         zhat    = pow(z1 + z2, x53)  - pow(z1,x53) - pow(z2,x53)
         zhat2   = pow(z1 + z2, x512) - pow(z1,x512) - pow(z2,x512)
         lzav    = x53 * log(max(1d-99,z1*z2/(z1 + z2)))
         aznut   = pow(z1*z1*z2*z2*a1*a2/(a1 + a2), x13)
               
      end subroutine screen5_init_AZ_info
      
      subroutine fxt_screen5(sc, zs13, zhat, zhat2, lzav, aznut, zs13inv,  &
                           a1, z1, a2, z2, scor, scordt, scordd, ierr)
!..this subroutine calculates screening factors and their derivatives
!..for nuclear reaction rates in the weak, intermediate and strong regimes. 

!..based on graboske, dewit, grossman and cooper apj 181 457 1973 for weak screening. 

!..based on alastuey and jancovici apj 226 1034 1978, 
!..with plasma parameters from itoh et al apj 234 1079 1979, for strong screening. 

!..input:
!..temp    = temperature
!..den     = density
!..zbar    = mean charge per nucleus
!..abar    = mean number of nucleons per nucleus 
!..z2bar   = mean square charge per nucleus
!..z1 a1   = charge and number in the entrance channel
!..z2 a2   = charge and number in the exit channel
!..jscreen = counter of which reaction is being calculated 
!..init    = flag to compute the more expensive functions just once

!..output:
!..scor    = screening correction
!..scordt  = derivative of screening correction with temperature
!..scordd  = derivative of screening correction with density


!..declare the pass
      type (Screen_Info), pointer :: sc
      real(dp), intent(in) :: zs13, zhat, zhat2, lzav, aznut, zs13inv
      ! zs13 = (z1+z2)**(1./3.)
      ! zhat = combination of z1 and z2 raised to the 5/3 power
      ! zhat2 = combination of z1 and z2 raised to the 5/12 power
      ! lzav = log of effective charge
      ! aznut = combination of a1, z1, a2, z2 raised to 1/3 power
      ! zs13inv = 1 / zs13
      real(dp), intent(in) :: a1, z1, a2, z2
      real(dp), intent(out) :: scor, scordt, scordd
      integer, intent(out) :: ierr

!..local variables
      real(dp) aa, daadt, daadd, bb, cc, dccdt, dccdd,  &
                       qq, dqqdt, dqqdd, rr, drrdt, drrdd,  &
                       ss, dssdt, dssdd, tt, dttdt, dttdd, uu, duudt, duudd,  &
                       vv, dvvdt, dvvdd, a3, da3, &
                       qlam0z, qlam0zdt, qlam0zdd,  &
                       h12w, dh12wdt, dh12wdd, h12, dh12dt, dh12dd,  &
                       taufac, taufacdt, gamp, gampdt, gampdd,  &
                       gamef, gamefdt, gamefdd,  &
                       tau12, tau12dt, alph12, alph12dt, alph12dd,  &
                       xlgfac, dxlgfacdt, dxlgfacdd,  &
                       gamp14, gamp14dt, gamp14dd, h12x, dh12xdt, dh12xdd,  &
                       gamefx, gamefs, alfa, beta, temp, den, zbar, abar, z2bar,&
                       dgamma

!..screening variables
!..zs13    = (z1+z2)**(1./3.)
!..zhat    = combination of z1 and z2 raised to the 5/3 power
!..zhat2   = combination of z1 and z2 raised to the 5/12 power
!..lzav    = log of effective charge
!..aznut   = combination of a1, z1, a2, z2 raised to 1/3 power


      real(dp), parameter :: alph12_lim = 1.6d0 ! ln(10)
      real(dp), parameter :: h12_max = 300d0
      real(dp), parameter :: x13   = 1.0d0/3.0d0 
      real(dp), parameter :: x14   = 1.0d0/4.0d0
      real(dp), parameter :: x53   = 5.0d0/3.0d0
      real(dp), parameter :: x532  = 5.0d0/32.0d0
      real(dp), parameter :: x512  = 5.0d0/12.0d0
      real(dp), parameter :: fact  = two_13 ! the cube root of 2
      real(dp), parameter :: co2   = x13 * 4.248719d3
      
      logical, parameter :: debug = .false.
      !logical :: debug
      
      include 'formats'


      !debug = (abs(a1 - 4d0) < 1d-14 .and. abs(a2 - 12d0) < 1d-14)
      
      ierr = 0
      
      if (z1 <= 0d0 .or. z2 <= 0d0) then
         scor = 1d0
         scordt = 0d0
         scordd = 0d0
         return
      end if

      zbar  = sc% zbar
      abar  = sc% abar
      z2bar = sc% z2bar
      temp  = sc% temp
      den   = sc% den

!..unload the screen info
      !ytot     = sc% ytot
      rr       = sc% rr
      !tempi    = sc% tempi
      !dtempi   = sc% dtempi
      !deni     = sc% deni
      
      !pp       = sc% pp 
      !dppdt    = sc% dppdt
      !dppdd    = sc% dppdd
      
      qlam0z   = sc% qlam0z
      qlam0zdt = sc% qlam0zdt
      qlam0zdd = sc% qlam0zdd
      
      taufac   = sc% taufac
      taufacdt = sc% taufacdt
      
      !xni      = sc% xni
      !dxnidd   = sc% dxnidd
      
      aa       = sc% aa
      daadt    = sc% daadt
      daadd    = sc% daadd


!..calculate individual screening factors 
      bb       = z1 * z2
      gamp     = aa
      gampdt   = daadt
      gampdd   = daadd

      qq       = fact * bb * zs13inv       
      gamef    = qq * gamp 
      gamefdt  = qq * gampdt
      gamefdd  = qq * gampdd

      tau12    = taufac * aznut 
      tau12dt  = taufacdt * aznut
      
      qq       = 1.0d0/tau12 
      alph12   = gamef * qq
      alph12dt = (gamefdt - alph12*tau12dt) * qq
      alph12dd = gamefdd * qq

!..limit alph12 to alph12_lim to prevent unphysical behavior.  
!..this should really be replaced by a pycnonuclear reaction rate formula 
      if (alph12 .gt. alph12_lim) then 

         alph12   = alph12_lim 
         alph12dt = 0.0d0
         alph12dd = 0.0d0

         gamef    = alph12 * tau12 
         gamefdt  = alph12 * tau12dt
         gamefdd  = 0.0d0 
         
         qq       = zs13/(fact * bb) 
         gamp     = gamef * qq
         gampdt   = gamefdt * qq
         gampdd   = 0.0d0
         
      end if 


!..weak screening regime 
      h12w    = bb * qlam0z 
      dh12wdt = bb * qlam0zdt
      dh12wdd = bb * qlam0zdd

      h12     = h12w 
      dh12dt  = dh12wdt
      dh12dd  = dh12wdd

!..intermediate and strong sceening regime

      if (debug) write(*, 1) 'gamef', gamef
      if (debug) write(*, 1) 'fact', fact
      if (debug) write(*, 1) 'z1', z1
      if (debug) write(*, 1) 'z2', z2
      if (debug) write(*, 1) 'gamp', gamp
      if (debug) write(*, 1) 'sc% aa', sc% aa
      if (debug) write(*, 1) 'sc% tempi', sc% tempi
      if (debug) write(*, 1) 'sc% xni', sc% xni
      
      gamefx = 0.3d0
      if (gamef .gt. gamefx) then
         if (debug) write(*,1) 'intermediate and strong sceening regime'

         gamp14   = pow(gamp,x14)
         rr       = 1.0d0/gamp
         qq       = 0.25d0*gamp14*rr 
         gamp14dt = qq * gampdt
         gamp14dd = qq * gampdd
   
         cc       = 0.896434d0 * gamp * zhat  &
                  - 3.44740d0  * gamp14 * zhat2   &
                  - 0.5551d0   * (log(gamp) + lzav)  &
                  - 2.996d0 

         dccdt    = 0.896434d0 * gampdt * zhat  &
                  - 3.44740d0  * gamp14dt * zhat2   &
                  - 0.5551d0*rr*gampdt

         dccdd    = 0.896434d0 * gampdd * zhat  &
                  - 3.44740d0  * gamp14dd * zhat2   &
                  - 0.5551d0*rr*gampdd

         a3     = alph12 * alph12 * alph12 
         da3    = 3.0d0 * alph12 * alph12
   
         qq     = 0.014d0 + 0.0128d0*alph12
         dqqdt  = 0.0128d0*alph12dt
         dqqdd  = 0.0128d0*alph12dd
   
         rr     = x532 - alph12*qq
         drrdt  = -(alph12dt*qq + alph12*dqqdt) 
         drrdd  = -(alph12dd*qq + alph12*dqqdd) 
   
         ss     = tau12*rr
         dssdt  = tau12dt*rr + tau12*drrdt
         dssdd  = tau12*drrdd
   
         tt     = -0.0098d0 + 0.0048d0*alph12
         dttdt  = 0.0048d0*alph12dt
         dttdd  = 0.0048d0*alph12dd
         
         uu     =  0.0055d0 + alph12*tt
         duudt  = alph12dt*tt + alph12*dttdt
         duudd  = alph12dd*tt + alph12*dttdd
   
         vv   = gamef * alph12 * uu  
         dvvdt= gamefdt*alph12*uu + gamef*alph12dt*uu + gamef*alph12*duudt  
         dvvdd= gamefdd*alph12*uu + gamef*alph12dd*uu + gamef*alph12*duudd  
         
         h12     = cc - a3 * (ss + vv)
         rr      = da3 * (ss + vv)
         dh12dt  = dccdt - rr*alph12dt - a3*(dssdt + dvvdt)
         dh12dd  = dccdd - rr*alph12dd - a3*(dssdd + dvvdd)
   
         rr     =  1.0d0 - 0.0562d0*a3
         ss     =  -0.0562d0*da3
         drrdt  = ss*alph12dt
         drrdd  = ss*alph12dd
   

         if (debug) write(*, 1) 'rr', rr

         if (rr .ge. 0.77d0) then
            xlgfac    = rr
            dxlgfacdt = drrdt 
            dxlgfacdd = drrdd 
         else
            xlgfac    = 0.77d0
            dxlgfacdt = 0.0d0
            dxlgfacdd = 0.0d0
         end if 

         h12    = log(xlgfac) + h12 
         rr     = 1.0d0/xlgfac
         dh12dt = rr*dxlgfacdt + dh12dt 
         dh12dd = rr*dxlgfacdd + dh12dd 
   
   

         if (debug) write(*, 1) 'gamef', gamef

         gamefs = 0.8d0
         if (gamef .le. gamefs) then 
             dgamma  = 1.0d0/(gamefs - gamefx)
      
             rr     =  dgamma*(gamefs - gamef)
             drrdt  = -dgamma*gamefdt
             drrdd  = -dgamma*gamefdd
      
             ss     = dgamma*(gamef - gamefx)
             dssdt  = dgamma*gamefdt
             dssdd  = dgamma*gamefdd
      
            vv     = h12 

            h12x    = h12
            dh12xdt = dh12dt
            dh12xdd = dh12dd
            
            h12    = h12w*rr + vv*ss
            dh12dt = dh12wdt*rr + h12w*drrdt + dh12dt*ss + vv*dssdt
            dh12dd = dh12wdd*rr + h12w*drrdd + dh12dd*ss + vv*dssdd

            if (debug) write(*, 1) 'gamef', gamef
            if (debug) write(*, 1) 'gamefx', gamefx

            if (debug) write(*,*) 'intermediate screening'
         
         else
            if (debug) write(*,*) 'strong screening'

         end if 
         
         !..end of intermediate and strong screening if
      
      else if (debug) then
         write(*,*) 'weak screening'

      end if 

      if (debug) write(*, 1) 'h12', h12
      if (debug) write(*, 1) 'h12/ln10', h12/ln10
      if (debug) write(*, *)

      !..machine limit the output
      h12    = max(min(h12, h12_max), 0.0d0) 
      scor   = exp(h12) 
      if (h12 .eq. h12_max) then
         scordt = 0.0d0
         scordd = 0.0d0
      else 
         scordt = scor * dh12dt
         scordd = scor * dh12dd
      end if

      if (debug) write(*, 1) 'scor', scor
      if (debug) write(*, 1) 'scordt', scordt
      if (debug) write(*, 1) 'scordd', scordd
      if (debug) write(*, *)
      if (debug) write(*, *)
      if (debug) write(*, *)

      end subroutine fxt_screen5


      
      



      end module screen5

