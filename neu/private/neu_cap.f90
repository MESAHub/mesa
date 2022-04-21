! ***********************************************************************
!
!   Copyright (C) 2022  The MESA Team
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

   module neu_cap
      use const_def
      use neu_def
      use num_lib
      use num_def
      use math_lib
      use utils_lib, only: mesa_error, is_bad

      implicit none

      ! Based of Lippuner & Roberts 2017 ApJS https://ui.adsabs.harvard.edu/abs/2017ApJS..233...18L/abstract


      real(dp), parameter :: Qec = (mn-mp) * clight2 ! ergs
      real(dp), parameter :: Qpc = -Qec ! ergs
      real(dp), parameter :: me_erg = me * clight2 ! ergs

      real(dp), parameter :: B = 5.76d0
      real(dp), parameter :: K = 6144.0d0
      real(dp), parameter :: C = (B*ln2)/(K*me_erg*me_erg*me_erg*me_erg*me_erg)

      ! Indexs to args
      integer, parameter :: IND_MUE = 1 ! electron chemical potential
      integer, parameter :: IND_TEMP = 2 ! temperature (kelvin)
      integer, parameter :: IND_Q = 3 ! Energy (ev)
      integer, parameter :: IND_NU = 4 ! neutrino distribution functions normally 0

      integer, parameter :: NUM_IND = 5

      ! Integration tolerances
      real(dp), parameter :: atol=1d-5, rtol=1d-5
      integer,  parameter :: MAX_STEPS = 20, MIN_STEPS=15
      real(dp), parameter :: MAX_EXP = 200d0



      contains
   
      real(dp) function fermidirac(temp, E, mu)
         real(dp),intent(in) :: temp ! Kelvin
         real(dp),intent(in) :: E ! erg
         real(dp),intent(in) :: mu ! ergs?
         real(dp) :: ex

         ex = (E/(boltzm*temp)) - mu

         if(ex>MAX_EXP) then
            fermidirac = 0d0
            return
         end if

         fermidirac = 1.d0/(exp(ex)+1d0)

         !write(90,*) log10(temp),mu,E,fermidirac,(E- mu),(E- mu)/(boltzm*temp),ex,(boltzm*temp)

      end function fermidirac


      real(dp) function w(Q) ! Lower bound on E integral
         real(dp), intent(in) :: Q ! erg

         w = max(me_erg,Q)

      end function w

      real(dp) function hwm() ! weak magnetism and recoil corrections Horowitz 2002
         !TODO: compute it

         hwm=0d0

      end function hwm

      real(dp) function lambda_ec(temp ,mu, ierr)
         real(dp),intent(in) :: temp, mu ! Temp in K mu electron chemical potential
         integer, intent(inout) :: ierr
         real(dp) :: args(NUM_IND)

         ierr = 0

         args(IND_TEMP) = temp
         args(IND_MUE) = mu + me_erg/(boltzm*temp)
         args(IND_NU) = 0
         args(IND_Q) = Qec

         lambda_ec = C * integrate_infinity(lambda1, w(args(IND_Q)), args, atol, rtol,MIN_STEPS,MAX_STEPS,ierr)

         ! write(70,*) log10(temp),mu, args(IND_MUE),w(args(IND_Q)),lambda_ec,i
         ! flush(70)


      end function lambda_ec

      real(dp) function lambda_pc(temp ,mu, ierr)
         real(dp),intent(in) :: temp, mu ! Temp in K mu electron chemical potential
         integer, intent(inout) :: ierr
         real(dp) :: args(NUM_IND)
         integer :: i

         ierr = 0

         args(IND_TEMP) = temp
         args(IND_MUE) = -mu - 2*me_erg/(boltzm*temp)
         args(IND_NU) = 0
         args(IND_Q) = Qpc

         lambda_pc = C * integrate_infinity(lambda1,w(args(IND_Q)),args,atol,rtol,MIN_STEPS,MAX_STEPS,ierr)

         !write(70,*) log10(temp),mu, args(IND_MUE),w(args(IND_Q)),lambda_pc,i
         !flush(70)

      end function lambda_pc

      real(dp) function lambda_nue(temp ,mu, ierr)
         real(dp),intent(in) :: temp, mu ! Temp in K mu electron chemical potential
         integer, intent(inout) :: ierr
         real(dp) :: args(NUM_IND)

         ierr = 0

         args(IND_TEMP) = temp
         args(IND_MUE) = mu + me_erg/(boltzm*temp)
         args(IND_NU) = 0
         args(IND_Q) = Qec

         lambda_nue = C * integrate_infinity(lambda2,w(args(IND_Q)),args,atol,rtol,MIN_STEPS,MAX_STEPS,ierr)

      end function lambda_nue

      real(dp) function lambda_anue(temp ,mu, ierr)
         real(dp),intent(in) :: temp, mu ! Temp in K mu electron chemical potential
         integer, intent(inout) :: ierr
         real(dp) :: args(NUM_IND)

         ierr = 0

         args(IND_TEMP) = temp
         args(IND_MUE) = -mu - 2*me_erg/(boltzm*temp)
         args(IND_NU) = 0
         args(IND_Q) = Qpc

         lambda_anue = C * integrate_infinity(lambda2,w(args(IND_Q)),args,atol,rtol,MIN_STEPS,MAX_STEPS,ierr)

      end function lambda_anue

      real(dp) function lambda1(E, args, ierr)
         real(dp),intent(in) :: E ! erg
         real(dp), intent(in) :: args(:)
         integer, intent(inout) :: ierr

         real(dp) :: temp,mu,nu,q
         real(dp) :: EmQ

         ierr = 0

         temp = args(IND_TEMP)
         mu= args(IND_MUE)
         nu = args(IND_NU)
         Q = args(IND_Q)

         EmQ = E-Q

         if(abs(E)<abs(Q)) then
            lambda1=0d0
            return
         end if 

         lambda1 = E * sqrt(E*E-Q*Q) * EmQ*EmQ * (1+E*hwm()) * &
                  fermidirac(temp, E, mu) * (1d0-fermidirac(temp, EmQ, nu))


        ! write(80,*) log10(temp),E,mu,sqrt(E*E-Q*Q),EmQ*EmQ,fermidirac(temp, E, mu),fermidirac(temp, EmQ, nu),lambda1

      end function lambda1

      real(dp) function lambda2(E, args, ierr)
         real(dp),intent(in) :: E ! erg
         real(dp), intent(in) :: args(:)
         integer, intent(inout) :: ierr

         real(dp) :: temp,mu,nu,q
         real(dp) :: EmQ

         ierr = 0

         temp = args(IND_TEMP)
         mu= args(IND_MUE)
         nu = args(IND_NU)
         Q = args(IND_Q)

         EmQ = E-Q

         if(abs(E)<abs(Q)) then
            lambda2=0d0
            return
         end if 

         lambda2 = E * sqrt(E*E-Q*Q) * EmQ*EmQ * (1+E*hwm()) * &
                     (1.d0-fermidirac(temp, E, mu)) * fermidirac(temp, EmQ, nu)

         !write(80,*) log10(temp),E,mu,Q,sqrt(E*E-Q*Q),EmQ*EmQ,fermidirac(temp, E, mu),fermidirac(temp, EmQ, nu),lambda2

      end function lambda2



   end module neu_cap
