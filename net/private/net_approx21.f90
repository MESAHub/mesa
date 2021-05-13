! ***********************************************************************
!
!   Copyright (C) 2014  Bill Paxton, Frank Timmes
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

      module net_approx21
      use const_def, only: dp, qp, avo, clight
      use utils_lib, only: is_bad, mesa_error


      implicit none

      

      !logical :: plus_co56 ! Must now be passed as an argument
      
      logical, parameter :: reduced_net_for_testing = .true.
      !logical, parameter :: reduced_net_for_testing = .false.

      integer, parameter :: species_21 = 21, species_co56 = 22


      integer :: iso_cid(species_co56) ! these are corresponding chem ids for the isos
         ! e.g., iso_cid(ife52) is = the iso number for fe52 as defined in mesa/chem
         ! Define as largest possible array
      integer :: & ! these are indices in y vector
         ih1, &
         ihe3, &
         ihe4, &
         ic12, &
         in14, &
         io16, &
         ine20, &
         img24, &
         isi28, &
         is32, &
         iar36, &
         ica40, &
         iti44, &
         icr48, &
         icrx, &
         ife52, &
         ife53, &
         ife54, &
         ife55, &
         ife56, &
         ico56, &
         ini56, &
         ineut, &
         iprot
      
      integer, parameter :: approx21_num_mesa_reactions_21 = 93, approx21_nrat = 116
      integer, parameter :: approx21_num_mesa_reactions_co56 = approx21_num_mesa_reactions_21+1, &
                              approx21_plus_co56_nrat = approx21_nrat+1
      
      ! integer :: num_mesa_reactions 
      ! integer :: num_reactions
      
      integer :: rate_id(approx21_num_mesa_reactions_co56) ! rate ids for the mesa reactions
         ! e.g., rate_id(ir3a) is reaction id for triple alpha as defined in mesa/rates
         ! Define as largest possible array
      integer :: & ! these are indices in rates arrays
         ir3a, &
         irg3a, &
         ircag, &
         ir1212, &
         ir1216, &
         ir1616, &
         iroga, &
         iroag, &
         irnega, &
         irneag, &
         irmgga, &
         irmgag, &
         irsiga, &
         irmgap, &
         iralpa, &
         iralpg, &
         irsigp, &
         irsiag, &
         irsga, &
         irsiap, &
         irppa, &
         irppg, &
         irsgp, &
         irsag, &
         irarga, &
         irsap, &
         irclpa, &
         irclpg, &
         irargp, &
         irarag, &
         ircaga, &
         irarap, &
         irkpa, &
         irkpg, &
         ircagp, &
         ircaag, &
         irtiga, &
         ircaap, &
         irscpa, &
         irscpg, &
         irtigp, &
         irtiag, &
         ircrga, &
         irtiap, &
         irvpa , &
         irvpg, &
         ircrgp, &
         ircrag, &
         irfega, &
         ircrap, &
         irmnpa, &
         irmnpg, &
         irfegp, &
         irfeag, &
         irniga, &
         irfeap, &
         ircopa, &
         ircopg, &
         irnigp, &

         ! for fe54 photodisintegration
         ir52ng, &
         ir53gn, &
         ir53ng, &
         ir54gn, &
         irfepg, &
         ircogp, &

         ! for he4 photodisintegration
         irheng, &
         irhegn, &
         irhng, &
         irdgn, &
         irdpg, &
         irhegp, &

         ! weak reactions
         irpen, &
         irnep, &
         irn56ec, &
         irco56ec, &

         ! ppchain
         irpp, &
         ir33, &
         irhe3ag, &
         ir_be7_wk_li7, &
         ir_be7_pg_b8, &

         ! cno cycles
         ircpg, &
         irnpg, &
         iropg, &
         irnag, &

         ! for reactions to fe56 
         ir54ng, &
         ir55gn, &
         ir55ng, &
         ir56gn, &
         irfe54ap, &
         irco57pa, &
         irfe56pg, &
         irco57gp, &

         ! for n15 branching
         irn15pa, &
         irn15pg, &

         ! the equilibrium links
         ifa, &
         ifg, &

         irr1, &
         irs1, &
         irt1, &
         iru1, &
         irv1, &
         irw1, &
         irx1, &

         ir1f54, &
         ir2f54, &
         ir3f54, &
         ir4f54, &
         ir5f54, &
         ir6f54, &
         ir7f54, &
         ir8f54, &

         iralf1, &
         iralf2, &

         irfe56_aux1, &
         irfe56_aux2, &
         irfe56_aux3, &
         irfe56_aux4

      ! names
      character (len=40) :: ratnam(approx21_plus_co56_nrat)
         ! Define as largest possible array

      contains



      ! call this after get raw rates
         subroutine approx21_pa_pg_fractions( &
            ratraw,dratrawdt,dratrawdd,ierr)
         real(dp), dimension(:) :: ratraw,dratrawdt,dratrawdd
         integer, intent(out) :: ierr
         
         include 'formats'
         
         ierr = 0

         call set1(ifa,irn15pg,irn15pa)
         ratraw(ifg)    = 1.0d0 - ratraw(ifa)
         dratrawdt(ifg) = -dratrawdt(ifa)
         dratrawdd(ifg) = -dratrawdd(ifa)

         call set1(irr1,iralpg,iralpa) ! al27
         call set1(irs1,irppg,irppa)   ! p31
         call set1(irt1,irclpg,irclpa) ! cl35
         call set1(iru1,irkpg,irkpa)   ! k39
         call set1(irv1,irscpg,irscpa) ! sc43
         call set1(irw1,irvpg,irvpa)   ! v47
         call set1(irx1,irmnpg,irmnpa) ! mn51


         contains

         subroutine set1(ifa,irn15pg,irn15pa)
            integer, intent(in) :: ifa,irn15pg,irn15pa
            real(dp) :: ff1,dff1dt,dff1dd,ff2,dff2dt,dff2dd, &
               tot,dtotdt,dtotdd,invtot

            ff1 = ratraw(irn15pg)
            dff1dt = dratrawdt(irn15pg)
            dff1dd = dratrawdd(irn15pg)

            ff2 = ratraw(irn15pa)
            dff2dt = dratrawdt(irn15pa)
            dff2dd = dratrawdd(irn15pa)

            tot            = ff1 + ff2
            dtotdt         = dff1dt + dff2dt
            dtotdd         = dff1dd + dff2dd

            if (tot > 1d-30) then
               invtot         = 1.0d0/tot
               ratraw(ifa)    = ff2 * invtot
               dratrawdt(ifa) = dff2dt * invtot - ff2 * invtot*invtot * dtotdt
               dratrawdd(ifa) = dff2dd * invtot - ff2 * invtot*invtot * dtotdd
            else
               ratraw(ifa)    = 0.0d0
               dratrawdt(ifa) = 0.0d0
               dratrawdd(ifa) = 0.0d0
            end if

         end subroutine set1

      end subroutine approx21_pa_pg_fractions

         
         ! call this before screening
         subroutine approx21_weak_rates( &
               y, ratraw, dratrawdt, dratrawdd, &
               temp, den, ye, eta, zbar, &
               weak_rate_factor,  plus_co56, ierr)
            use rates_lib, only: eval_ecapnuc_rate
            use net_derivs, only: eval_ni56_ec_rate, eval_co56_ec_rate
            
            real(dp), dimension(:) :: y, ratraw, dratrawdt, dratrawdd
            real(dp), intent(in) :: temp, den, ye, eta, zbar, weak_rate_factor
            logical, intent(in) :: plus_co56
            integer, intent(out) :: ierr
            
            real(dp) :: rpen, rnep, spen, snep, &
               rate, d_rate_dlnT, d_rate_dlnRho, Q, Qneu
            include 'formats'
            
            ierr = 0
            
            call eval_ecapnuc_rate(eta, temp, rpen, rnep, spen, snep)
            
            ratraw(irpen) = rpen
            dratrawdt(irpen) = 0
            dratrawdd(irpen) = 0
            if (rpen > 0) then
               Qneu = spen/rpen
            else
               Qneu = 0
            end if
            
            ratraw(irnep) = rnep
            dratrawdt(irnep) = 0
            dratrawdd(irnep) = 0
            if (rnep > 0) then
               Qneu = snep/rnep
            else
               Qneu = 0
            end if
            
            call eval_ni56_ec_rate( &
               temp, den, ye, eta, zbar, weak_rate_factor, &
               rate, d_rate_dlnT, d_rate_dlnRho, Q, Qneu, &
               ierr)
            if (ierr /= 0) then
               !write(*,*) 'failed in eval_ni56_ec_rate'
               return
            end if
            ratraw(irn56ec) = rate
            dratrawdt(irn56ec) = 0 
            dratrawdd(irn56ec) = 0
            
            if (plus_co56) then         
               call eval_co56_ec_rate( &
                  temp, den, ye, eta, zbar, weak_rate_factor, &
                  rate, d_rate_dlnT, d_rate_dlnRho, Q, Qneu, &
                  ierr)
               if (ierr /= 0) then
                  !write(*,*) 'failed in eval_co56_ec_rate'
                  return
               end if
               ratraw(irco56ec) = rate
               dratrawdt(irco56ec) = 0 
               dratrawdd(irco56ec) = 0           
            end if

         end subroutine approx21_weak_rates


         ! call this after screening -- depends on y, so don't reuse results.
         subroutine approx21_special_reactions( &
               btemp, bden, abar, zbar, y, &
               use_3a_FL, conv_eps_3a, &
               ratdum, dratdumdt, dratdumdd, dratdumdy1, dratdumdy2, &
               plus_co56, ierr)
            use math_lib
            use rates_lib, only: eval_FL_epsnuc_3alf
            use utils_lib, only: mesa_error
            real(dp), intent(in) :: &
               btemp, bden, abar, zbar, y(:), conv_eps_3a
            logical, intent(in) :: use_3a_FL
            real(dp), dimension(:) :: &
               ratdum,dratdumdt,dratdumdd,dratdumdy1,dratdumdy2
            logical, intent(in) ::  plus_co56
            integer, intent(out) :: ierr
            
            real(dp) :: denom, denomdt, denomdd, zz, xx, eps, deps_dT, deps_dRho
            real(dp), parameter :: tiny_denom = 1d-50, tiny_y = 1d-30
            integer :: i
            logical :: okay
            include 'formats'
            
            ierr = 0
            
            
         if (use_3a_FL) then
            ! Fushiki and Lamb, Apj, 317, 368-388, 1987
            if (y(ihe4) < tiny_y) then
               ratdum(ir3a)     = 0.0d0
               dratdumdt(ir3a)  = 0.0d0
               dratdumdd(ir3a)  = 0.0d0
            else
               call eval_FL_epsnuc_3alf( &
                  btemp, bden, 4*y(ihe4), abar/zbar, eps, deps_dT, deps_dRho)
               ! convert from eps back to rate
               xx = conv_eps_3a*y(ihe4)*y(ihe4)*y(ihe4)/6d0
               ratdum(ir3a) = eps/xx
               dratdumdt(ir3a) = deps_dT/xx
               dratdumdd(ir3a) = deps_dRho/xx    
            end if
         end if
            
            
            okay = .true.
            do i=1,num_mesa_reactions(plus_co56) 
               if (ratdum(i) < 0d0) then
                  write(*,2) 'approx21 missing rate for ' // ratnam(i), i, ratdum(i), &
                     btemp, log10(btemp), bden, log10(bden)
                  okay = .false.
               end if
            end do
            if (.not. okay) call mesa_error(__FILE__,__LINE__)

            ! for debugging: sum(cat)/eps_nuc
            
            if (reduced_net_for_testing) then 

               !if (.true.) then 
               
               !end if

               if (.false.) then 

                  ! turn off PP
                  call turn_off_reaction(irpp)
                  call turn_off_reaction(ir33)
                  call turn_off_reaction(irhe3ag)

                  !turn off CNO
                  call turn_off_reaction(ircpg)
                  call turn_off_reaction(irnpg)
                  call turn_off_reaction(iropg)

                  ! turn off n14ag + 3alpha
                  call turn_off_reaction(irnag)
                  call turn_off_reaction(ir3a)
                  call turn_off_reaction(irg3a)

                  !c12+c12, c12ag, o16ag
                  call turn_off_reaction(ir1212)
                  call turn_off_reaction(iroag)
                  call turn_off_reaction(irnega)
                  call turn_off_reaction(ircag)
                  call turn_off_reaction(iroga)

                  !Ne/O burn
                  call turn_off_reaction(ir1216)
                  call turn_off_reaction(ir1616) 
                  call turn_off_reaction(irneag)
                  call turn_off_reaction(irmgga)

                  !alpha links
                  call turn_off_reaction(irmgag)
                  call turn_off_reaction(irsiga)
                  call turn_off_reaction(irmgap)
                  call turn_off_reaction(iralpa)
                  call turn_off_reaction(iralpg)
                  call turn_off_reaction(irsigp)
                  call turn_off_reaction(irsiag)
                  call turn_off_reaction(irsga)

                  call turn_off_reaction(irppa)
                  call turn_off_reaction(irppg)
                  call turn_off_reaction(irsiap)
                  call turn_off_reaction(irsgp)

                  call turn_off_reaction(irsag)
                  call turn_off_reaction(irarga)
                  call turn_off_reaction(irsap)
                  call turn_off_reaction(irclpa)
                  call turn_off_reaction(irclpg)
                  call turn_off_reaction(irargp)

                  call turn_off_reaction(irarag)
                  call turn_off_reaction(ircaga)
                  call turn_off_reaction(irarap)
                  call turn_off_reaction(irkpa)
                  call turn_off_reaction(irkpg)
                  call turn_off_reaction(ircagp)

                  call turn_off_reaction(ircaag)
                  call turn_off_reaction(irtiga)
                  call turn_off_reaction(ircaap)
                  call turn_off_reaction(irscpa)
                  call turn_off_reaction(irscpg)
                  call turn_off_reaction(irtigp)
                  call turn_off_reaction(irtiag)

                  call turn_off_reaction(ircrga)
                  call turn_off_reaction(irtiap)
                  call turn_off_reaction(irvpa )
                  call turn_off_reaction(irvpg)
                  call turn_off_reaction(ircrgp)
                  call turn_off_reaction(ircrag)
                  call turn_off_reaction(irfega)
                  call turn_off_reaction(ircrap)
                  call turn_off_reaction(irmnpa)
                  call turn_off_reaction(irmnpg)
                  call turn_off_reaction(irfegp)
                  call turn_off_reaction(irfeag)

                  !iron group
                  call turn_off_reaction(irniga)                                  
                  call turn_off_reaction(irfeap)
                  call turn_off_reaction(ircopa)
               
                  call turn_off_reaction(irnigp)
                  call turn_off_reaction(irfepg)
                  call turn_off_reaction(ircogp)
               
                  call turn_off_reaction(irheng)
                  call turn_off_reaction(irhegn)
                  
                  call turn_off_reaction(irhng)
                  call turn_off_reaction(irdgn)
                  
                  call turn_off_reaction(irdpg)
                  call turn_off_reaction(irhegp)
                  
                  call turn_off_reaction(irpen)
                  call turn_off_reaction(irnep)
            
                  call turn_off_reaction(ircopg)
               
                  call turn_off_reaction(ir54ng)
                  call turn_off_reaction(ir55gn)
                  call turn_off_reaction(ir55ng)
                  call turn_off_reaction(ir56gn)

                  call turn_off_reaction(ir52ng)
                  call turn_off_reaction(ir53gn)
                  call turn_off_reaction(ir53ng)
                  call turn_off_reaction(ir54gn)

                  call turn_off_reaction(irfe56pg)
                  call turn_off_reaction(irco57gp)

                  call turn_off_reaction(irfe54ap)
                  call turn_off_reaction(irco57pa)
                  
                  call turn_off_reaction(irco56ec)
                  call turn_off_reaction(irn56ec)
               
               end if

               !if (.true.) then 
                  

               !end if 

            end if
            
   ! fe52(n,g)fe53(n,g)fe54 equilibrium links
         ratdum(ir1f54)     = 0.0d0
         dratdumdy1(ir1f54) = 0.0d0
         dratdumdt(ir1f54)  = 0.0d0
         dratdumdd(ir1f54)  = 0.0d0

         ratdum(ir2f54)     = 0.0d0
         dratdumdy1(ir2f54) = 0.0d0
         dratdumdt(ir2f54)  = 0.0d0
         dratdumdd(ir2f54)  = 0.0d0

         denom   = ratdum(ir53gn) + y(ineut)*ratdum(ir53ng)
         denomdt = dratdumdt(ir53gn) + y(ineut)*dratdumdt(ir53ng)
         denomdd = dratdumdd(ir53gn) + y(ineut)*dratdumdd(ir53ng)

         if (denom > tiny_denom .and. btemp .gt. 1.5d9) then
         zz      = 1.0d0/denom

         ratdum(ir1f54)     = ratdum(ir54gn)*ratdum(ir53gn)*zz
         dratdumdy1(ir1f54) = -ratdum(ir1f54)*zz * ratdum(ir53ng)
         dratdumdt(ir1f54)  = dratdumdt(ir54gn)*ratdum(ir53gn)*zz &
                              + ratdum(ir54gn)*dratdumdt(ir53gn)*zz &
                              - ratdum(ir1f54)*zz*denomdt
         dratdumdd(ir1f54) = dratdumdd(ir54gn)*ratdum(ir53gn)*zz &
                           + ratdum(ir54gn)*dratdumdd(ir53gn)*zz &
                           - ratdum(ir1f54)*zz*denomdd

         ratdum(ir2f54)     = ratdum(ir52ng)*ratdum(ir53ng)*zz
         dratdumdy1(ir2f54) = -ratdum(ir2f54)*zz * ratdum(ir53ng)
         dratdumdt(ir2f54)  = dratdumdt(ir52ng)*ratdum(ir53ng)*zz &
                              + ratdum(ir52ng)*dratdumdt(ir53ng)*zz &
                              - ratdum(ir2f54)*zz*denomdt
         dratdumdd(ir2f54) = dratdumdd(ir52ng)*ratdum(ir53ng)*zz &
                           + ratdum(ir52ng)*dratdumdd(ir53ng)*zz &
                           - ratdum(ir2f54)*zz*denomdd
         end if

   ! fe54(n,g)fe55(n,g)fe56 equilibrium links
         ratdum(irfe56_aux1)     = 0.0d0
         dratdumdy1(irfe56_aux1) = 0.0d0
         dratdumdt(irfe56_aux1)  = 0.0d0
         dratdumdd(irfe56_aux1)  = 0.0d0

         ratdum(irfe56_aux2)     = 0.0d0
         dratdumdy1(irfe56_aux2) = 0.0d0
         dratdumdt(irfe56_aux2)  = 0.0d0
         dratdumdd(irfe56_aux2)  = 0.0d0

         denom   = ratdum(ir55gn)    + y(ineut)*ratdum(ir55ng)
         denomdt = dratdumdt(ir55gn) + y(ineut)*dratdumdt(ir55ng)
         denomdd = dratdumdd(ir55gn) + y(ineut)*dratdumdd(ir55ng)

         if (denom > tiny_denom .and. btemp .gt. 1.5d9) then
         zz      = 1.0d0/denom

         ratdum(irfe56_aux1)     = ratdum(ir56gn)*ratdum(ir55gn)*zz
         dratdumdy1(irfe56_aux1) = -ratdum(irfe56_aux1)*zz * ratdum(ir55ng)
         dratdumdt(irfe56_aux1)  = dratdumdt(ir56gn)*ratdum(ir55gn)*zz &
                                 + ratdum(ir56gn)*dratdumdt(ir55gn)*zz &
                                 - ratdum(irfe56_aux1)*zz*denomdt
         dratdumdd(irfe56_aux1)  = dratdumdd(ir56gn)*ratdum(ir55gn)*zz &
                                 + ratdum(ir56gn)*dratdumdd(ir55gn)*zz &
                                 - ratdum(irfe56_aux1)*zz*denomdd

         ratdum(irfe56_aux2)     = ratdum(ir54ng)*ratdum(ir55ng)*zz
         dratdumdy1(irfe56_aux2) = -ratdum(irfe56_aux2)*zz * ratdum(ir55ng)
         dratdumdt(irfe56_aux2)  = dratdumdt(ir54ng)*ratdum(ir55ng)*zz &
                                 + ratdum(ir54ng)*dratdumdt(ir55ng)*zz &
                                 - ratdum(irfe56_aux2)*zz*denomdt
         dratdumdd(irfe56_aux2) = dratdumdd(ir54ng)*ratdum(ir55ng)*zz &
                                 + ratdum(ir54ng)*dratdumdd(ir55ng)*zz &
                                 - ratdum(irfe56_aux2)*zz*denomdd
         
         end if

   ! fe54(a,p)co57(g,p)fe56 equilibrium links 

         ratdum(irfe56_aux3)     = 0.0d0
         dratdumdy1(irfe56_aux3) = 0.0d0
         dratdumdt(irfe56_aux3)  = 0.0d0
         dratdumdd(irfe56_aux3)  = 0.0d0

         ratdum(irfe56_aux4)     = 0.0d0
         dratdumdy1(irfe56_aux4) = 0.0d0
         dratdumdt(irfe56_aux4)  = 0.0d0
         dratdumdd(irfe56_aux4)  = 0.0d0

         denom   = ratdum(irco57gp)    + y(iprot)*ratdum(irco57pa)
         denomdt = dratdumdt(irco57gp) + y(iprot)*dratdumdt(irco57pa)
         denomdd = dratdumdd(irco57gp) + y(iprot)*dratdumdd(irco57pa)

         if (denom > tiny_denom .and. btemp .gt. 1.5d9) then
         zz      = 1.0d0/denom

         ratdum(irfe56_aux3)     = ratdum(irfe56pg) * ratdum(irco57pa) * zz
         dratdumdy1(irfe56_aux3) = -ratdum(irfe56_aux3) * zz * ratdum(irco57pa)
         dratdumdt(irfe56_aux3)  = dratdumdt(irfe56pg) * ratdum(irco57pa) * zz &
                                 + ratdum(irfe56pg) * dratdumdt(irco57pa) * zz &
                                 - ratdum(irfe56_aux3) * zz * denomdt
         dratdumdd(irfe56_aux3)  = dratdumdd(irfe56pg) * ratdum(irco57pa) * zz &
                                 + ratdum(irfe56pg) * dratdumdd(irco57pa) * zz &
                                 - ratdum(irfe56_aux3) * zz * denomdd

         ratdum(irfe56_aux4)     = ratdum(irfe54ap) * ratdum(irco57gp) * zz
         dratdumdy1(irfe56_aux4) = -ratdum(irfe56_aux4) * zz * ratdum(irco57pa)
         dratdumdt(irfe56_aux4)  = dratdumdt(irfe54ap) * ratdum(irco57gp) * zz &
                                 + ratdum(irfe54ap) * dratdumdt(irco57gp) * zz &
                                 - ratdum(irfe56_aux4) * zz * denomdt
         dratdumdd(irfe56_aux4)  = dratdumdd(irfe54ap) * ratdum(irco57gp) * zz &
                                 + ratdum(irfe54ap) * dratdumdd(irco57gp) * zz &
                                 - ratdum(irfe56_aux4) * zz * denomdd
         end if


   ! fe54(p,g)co55(p,g)ni56 equilibrium links r3f54 r4f54
   ! fe52(a,p)co55(g,p)fe54 equilibrium links r5f54 r6f54
   ! fe52(a,p)co55(p,g)ni56 equilibrium links r7f54 r8f54

         ratdum(ir3f54)     = 0.0d0
         dratdumdy1(ir3f54) = 0.0d0
         dratdumdt(ir3f54)  = 0.0d0
         dratdumdd(ir3f54)  = 0.0d0

         ratdum(ir4f54)     = 0.0d0
         dratdumdy1(ir4f54) = 0.0d0
         dratdumdt(ir4f54)  = 0.0d0
         dratdumdd(ir4f54)  = 0.0d0

         ratdum(ir5f54)     = 0.0d0
         dratdumdy1(ir5f54) = 0.0d0
         dratdumdt(ir5f54)  = 0.0d0
         dratdumdd(ir5f54)  = 0.0d0

         ratdum(ir6f54)     = 0.0d0
         dratdumdy1(ir6f54) = 0.0d0
         dratdumdt(ir6f54)  = 0.0d0
         dratdumdd(ir6f54)  = 0.0d0

         ratdum(ir7f54)     = 0.0d0
         dratdumdy1(ir7f54) = 0.0d0
         dratdumdt(ir7f54)  = 0.0d0
         dratdumdd(ir7f54)  = 0.0d0

         ratdum(ir8f54)     = 0.0d0
         dratdumdy1(ir8f54) = 0.0d0
         dratdumdt(ir8f54)  = 0.0d0
         dratdumdd(ir8f54)  = 0.0d0

         denom   = ratdum(ircogp)+y(iprot)*(ratdum(ircopg)+ratdum(ircopa))

         if (denom > tiny_denom .and. btemp .gt. 1.5d9) then

         denomdt = dratdumdt(ircogp) &
                  + y(iprot)*(dratdumdt(ircopg) + dratdumdt(ircopa))
         denomdd = dratdumdd(ircogp) &
                  + y(iprot)*(dratdumdd(ircopg) + dratdumdd(ircopa))

         zz      = 1.0d0/denom

         ratdum(ir3f54)     = ratdum(irfepg) * ratdum(ircopg) * zz
         dratdumdy1(ir3f54) = -ratdum(ir3f54) * zz * &
                              (ratdum(ircopg) + ratdum(ircopa))
         dratdumdt(ir3f54)  = dratdumdt(irfepg) * ratdum(ircopg) * zz &
                  + ratdum(irfepg) * dratdumdt(ircopg) * zz &
                  - ratdum(ir3f54)*zz*denomdt
         dratdumdd(ir3f54)  = dratdumdd(irfepg) * ratdum(ircopg) * zz &
                  + ratdum(irfepg) * dratdumdd(ircopg) * zz &
                  - ratdum(ir3f54)*zz*denomdd

         ratdum(ir4f54)     = ratdum(irnigp) * ratdum(ircogp) * zz
         dratdumdy1(ir4f54) = -ratdum(ir4f54) * zz * &
                              (ratdum(ircopg)+ratdum(ircopa))
         dratdumdt(ir4f54)  =  dratdumdt(irnigp) * ratdum(ircogp) * zz &
                  + ratdum(irnigp) * dratdumdt(ircogp) * zz &
                  - ratdum(ir4f54)*zz*denomdt
         dratdumdd(ir4f54)  = dratdumdd(irnigp) * ratdum(ircogp) * zz &
                  + ratdum(irnigp) * dratdumdd(ircogp) * zz &
                  - ratdum(ir4f54)*zz*denomdd

         ratdum(ir5f54)     = ratdum(irfepg) * ratdum(ircopa) * zz
         dratdumdy1(ir5f54) = -ratdum(ir5f54) * zz * &
                              (ratdum(ircopg)+ratdum(ircopa))
         dratdumdt(ir5f54)  = dratdumdt(irfepg) * ratdum(ircopa) * zz &
                  + ratdum(irfepg) * dratdumdt(ircopa) * zz &
                  - ratdum(ir5f54) * zz * denomdt
         dratdumdd(ir5f54)  = dratdumdd(irfepg) * ratdum(ircopa) * zz &
                  + ratdum(irfepg) * dratdumdd(ircopa) * zz &
                  - ratdum(ir5f54) * zz * denomdd

         ratdum(ir6f54)     = ratdum(irfeap) * ratdum(ircogp) * zz
         dratdumdy1(ir6f54) = -ratdum(ir6f54) * zz * &
                              (ratdum(ircopg)+ratdum(ircopa))
         dratdumdt(ir6f54)  = dratdumdt(irfeap) * ratdum(ircogp) * zz &
                  + ratdum(irfeap) * dratdumdt(ircogp) * zz &
                  - ratdum(ir6f54) * zz * denomdt
         dratdumdd(ir6f54)  = dratdumdd(irfeap) * ratdum(ircogp) * zz &
                  + ratdum(irfeap) * dratdumdd(ircogp) * zz &
                  - ratdum(ir6f54) * zz * denomdd

         ratdum(ir7f54)     = ratdum(irfeap) * ratdum(ircopg) * zz

         dratdumdy1(ir7f54) = -ratdum(ir7f54) * zz * &
                              (ratdum(ircopg)+ratdum(ircopa))
         dratdumdt(ir7f54)  = dratdumdt(irfeap) * ratdum(ircopg) * zz &
                  + ratdum(irfeap) * dratdumdt(ircopg) * zz &
                  - ratdum(ir7f54) * zz * denomdt
         dratdumdd(ir7f54)  = dratdumdd(irfeap) * ratdum(ircopg) * zz &
                  + ratdum(irfeap) * dratdumdd(ircopg) * zz &
                  - ratdum(ir7f54) * zz * denomdd

         ratdum(ir8f54)     = ratdum(irnigp) * ratdum(ircopa) * zz

         dratdumdy1(ir8f54) = -ratdum(ir8f54) * zz * &
                              (ratdum(ircopg)+ratdum(ircopa))
         dratdumdt(ir8f54)  = dratdumdt(irnigp) * ratdum(ircopa) * zz &
                  + ratdum(irnigp) * dratdumdt(ircopa) * zz &
                  - ratdum(ir8f54) * zz * denomdt
         dratdumdd(ir8f54)  = dratdumdd(irnigp) * ratdum(ircopa) * zz &
                  + ratdum(irnigp) * dratdumdd(ircopa) * zz &
                  - ratdum(ir8f54) * zz * denomdd
                  

         end if


   ! p(n,g)h2(n,g)3h(p,g)he4   photodisintegrated n and p back to he4 equilibrium links
   ! p(n,g)h2(p,g)he3(n,g)he4

         ratdum(iralf1)     = 0.0d0
         dratdumdy1(iralf1) = 0.0d0
         dratdumdy2(iralf1) = 0.0d0
         dratdumdt(iralf1)  = 0.0d0
         dratdumdd(iralf1)  = 0.0d0

         ratdum(iralf2)     = 0.0d0
         dratdumdy1(iralf2) = 0.0d0
         dratdumdy2(iralf2) = 0.0d0
         dratdumdt(iralf2)  = 0.0d0
         dratdumdd(iralf2)  = 0.0d0

         denom  = ratdum(irhegp)*ratdum(irdgn) + &
                  y(ineut)*ratdum(irheng)*ratdum(irdgn) + &
                  y(ineut)*y(iprot)*ratdum(irheng)*ratdum(irdpg)

         if (denom > tiny_denom .and. btemp .gt. 1.5d9) then

         denomdt  = dratdumdt(irhegp)*ratdum(irdgn) &
                  + ratdum(irhegp)*dratdumdt(irdgn) &
                  +  y(ineut) * (dratdumdt(irheng)*ratdum(irdgn) &
                              + ratdum(irheng)*dratdumdt(irdgn)) &
                  +  y(ineut)*y(iprot) * (dratdumdt(irheng)*ratdum(irdpg) &
                                       + ratdum(irheng)*dratdumdt(irdpg))

         denomdd  = dratdumdd(irhegp)*ratdum(irdgn) &
                  + ratdum(irhegp)*dratdumdd(irdgn) &
                  +  y(ineut) * (dratdumdd(irheng)*ratdum(irdgn) &
                              + ratdum(irheng)*dratdumdd(irdgn)) &
                  +  y(ineut)*y(iprot) * (dratdumdd(irheng)*ratdum(irdpg) &
                                       + ratdum(irheng)*dratdumdd(irdpg))

         zz = 1.0d0/denom

         ratdum(iralf1)     = ratdum(irhegn) * ratdum(irhegp)* &
                              ratdum(irdgn) * zz
         dratdumdy1(iralf1) = -ratdum(iralf1) * zz * &
                              (ratdum(irheng)*ratdum(irdgn) + &
                              y(iprot)*ratdum(irheng)*ratdum(irdpg))
         dratdumdy2(iralf1) = -ratdum(iralf1) * zz * y(ineut) * &
                              ratdum(irheng) * ratdum(irdpg)
         dratdumdt(iralf1)  = dratdumdt(irhegn)*ratdum(irhegp)* &
                              ratdum(irdgn) * zz &
                  + ratdum(irhegn)*dratdumdt(irhegp)*ratdum(irdgn)*zz &
                  + ratdum(irhegn)*ratdum(irhegp)*dratdumdt(irdgn)*zz &
                  - ratdum(iralf1)*zz*denomdt
         dratdumdd(iralf1)  = dratdumdd(irhegn) * ratdum(irhegp)* &
                              ratdum(irdgn) * zz &
                  + ratdum(irhegn)*dratdumdd(irhegp)*ratdum(irdgn)*zz &
                  + ratdum(irhegn)*ratdum(irhegp)*dratdumdd(irdgn)*zz &
                  - ratdum(iralf1)*zz*denomdt


         ratdum(iralf2)     = ratdum(irheng)*ratdum(irdpg)* &
                              ratdum(irhng)*zz
         dratdumdy1(iralf2) = -ratdum(iralf2) * zz * &
                              (ratdum(irheng)*ratdum(irdgn) + &
                                 y(iprot)*ratdum(irheng)*ratdum(irdpg))


         denom  = ratdum(irhegp)*ratdum(irdgn) + &
                  y(ineut)*ratdum(irheng)*ratdum(irdgn) + &
                  y(ineut)*y(iprot)*ratdum(irheng)*ratdum(irdpg)
                                 
         if (is_bad(dratdumdy1(iralf2))) then
            write(*,1) 'denom', denom
            write(*,1) 'zz', zz
            write(*,1) 'dratdumdy1(iralf2)', dratdumdy1(iralf2)
            write(*,1) 'ratdum(irhegp)*ratdum(irdgn)', ratdum(irhegp)*ratdum(irdgn)
            write(*,1) 'y(ineut)*ratdum(irheng)*ratdum(irdgn)', y(ineut)*ratdum(irheng)*ratdum(irdgn)
            write(*,1) 'y(ineut)*y(iprot)*ratdum(irheng)*ratdum(irdpg)', &
               y(ineut)*y(iprot)*ratdum(irheng)*ratdum(irdpg)
            write(*,1) 'ratdum(irhegp)', ratdum(irhegp)
            write(*,1) 'ratdum(irdgn)', ratdum(irdgn)
            write(*,1) 'ratdum(irheng)', ratdum(irheng)
            write(*,1) 'ratdum(irdgn)', ratdum(irdgn)
            write(*,1) 'y(ineut)', y(ineut)
            write(*,1) 'y(iprot)', y(iprot)
            stop
         end if
         
         
         dratdumdy2(iralf2) = -ratdum(iralf2) * zz * y(ineut)* &
                              ratdum(irheng) * ratdum(irdpg)
         dratdumdt(iralf2)  = dratdumdt(irheng)*ratdum(irdpg) * &
                              ratdum(irhng) * zz &
                  + ratdum(irheng)*dratdumdt(irdpg)*ratdum(irhng)*zz &
                  + ratdum(irheng)*ratdum(irdpg)*dratdumdt(irhng)*zz &
                  - ratdum(iralf2)*zz*denomdt
         dratdumdd(iralf2)  = dratdumdd(irheng)*ratdum(irdpg)* &
                              ratdum(irhng)*zz &
                  + ratdum(irheng)*dratdumdd(irdpg)*ratdum(irhng)*zz &
                  + ratdum(irheng)*ratdum(irdpg)*dratdumdd(irhng)*zz &
                  - ratdum(iralf2)*zz*denomdd

         end if



   ! he3(a,g)be7(p,g)8b(e+nu)8be(2a)
   ! beta limit he3+he4 by the 8B decay half life
         if (y(ihe4) > tiny_y) then
         xx            = 0.896d0/y(ihe4)
         ratdum(irhe3ag)  = min(ratdum(irhe3ag),xx)
         if (ratdum(irhe3ag) .eq. xx) then
         dratdumdy1(irhe3ag) = -xx/y(ihe4)
         dratdumdt(irhe3ag)  = 0.0d0
         dratdumdd(irhe3ag)  = 0.0d0
         else
         dratdumdy1(irhe3ag) = 0.0d0
         endif
         endif


   ! beta limit n14(p,g)o15(enu)o16  and o16(p,g)f17(e+nu)17o(p,a)n14
         if (y(ih1) > tiny_y) then

            xx = 5.68d-03/(y(ih1)*1.57d0)
            ratdum(irnpg) = min(ratdum(irnpg),xx)
            if (ratdum(irnpg) .eq. xx) then
            dratdumdy1(irnpg) = -xx/y(ih1)
            dratdumdt(irnpg)  = 0.0d0
            dratdumdd(irnpg)  = 0.0d0
            else
            dratdumdy1(irnpg) = 0.0d0
            end if

            xx = 0.0105d0/y(ih1)
            ratdum(iropg) = min(ratdum(iropg),xx)
            if (ratdum(iropg) .eq. xx) then
            dratdumdy1(iropg) = -xx/y(ih1)
            dratdumdt(iropg)  = 0.0d0
            dratdumdd(iropg)  = 0.0d0
            else
            dratdumdy1(iropg) = 0.0d0
            end if

         end if
         
            
         contains
                  
         subroutine turn_off_reaction(i)
            integer, intent(in) :: i
            if (i == 0) return
            ratdum(i) = 0
            dratdumdt(i) = 0
            dratdumdd(i) = 0
            dratdumdy1(i) = 0
            dratdumdy2(i) = 0
         end subroutine turn_off_reaction         

         end subroutine approx21_special_reactions
         

         subroutine approx21_dydt( &
            y, rate, ratdum, dydt, deriva, &
            fe56ec_fake_factor_in, min_T, fe56ec_n_neut, temp, plus_co56, ierr)
         logical, intent(in) :: deriva ! false for dydt, true for partials wrt T, Rho
         real(dp), dimension(:), intent(in) :: y, rate, ratdum
         integer, intent(in) :: fe56ec_n_neut
         real(dp), dimension(:), intent(out) :: dydt
         real(dp), intent(in) :: fe56ec_fake_factor_in, temp
         real(dp) :: fe56ec_fake_factor, min_T
         logical, intent(in) ::  plus_co56
         integer, intent(out) :: ierr

         integer :: i

   ! quad precision dydt sums
         real(qp) :: a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,&
            a11,a12,a13,a14,a15,a16,a17,a18,a19,a20
         real(qp) :: qray(species(plus_co56))

         logical :: okay
         
         include 'formats'
         
         ierr = 0

         ! Turn on special fe56ec rate above some temperature
         fe56ec_fake_factor=eval_fe56ec_fake_factor(fe56ec_fake_factor_in, min_T, temp)
         
         dydt(1:species(plus_co56)) = 0.0d0
         qray(1:species(plus_co56)) = 0.0_qp

   ! hydrogen reactions
         a1 = -1.5d0 * y(ih1) * y(ih1) * rate(irpp)
         a2 =  y(ihe3) * y(ihe3) * rate(ir33) 
         a3 = -y(ihe3) * y(ihe4) * rate(irhe3ag) 
         a4 = -2.0d0 * y(ic12) * y(ih1) * rate(ircpg) 
         a5 = -2.0d0 * y(in14) * y(ih1) * rate(irnpg) 
         a6 = -2.0d0 * y(io16) * y(ih1) * rate(iropg) 
         a7 = -3.0d0 * y(ih1) * rate(irpen)

         qray(ih1) = qray(ih1) + a1 + a2 + a3 + a4 + a5 + a6 + a7      

   ! he3 reactions

         a1  =  0.5d0 * y(ih1) * y(ih1) * rate(irpp) 
         a2  = -y(ihe3) * y(ihe3) * rate(ir33) 
         a3  = -y(ihe3) * y(ihe4) * rate(irhe3ag) 
         a4  =  y(ih1) * rate(irpen)

         qray(ihe3) = qray(ihe3) + a1 + a2 + a3 + a4


   ! he4 reactions
   ! heavy ion reactions
         a1  = 0.5d0 * y(ic12) * y(ic12) * rate(ir1212) 
         a2  = 0.5d0 * y(ic12) * y(io16) * rate(ir1216) 
         a3  = 0.56d0 * 0.5d0 * y(io16) * y(io16) * rate(ir1616)
         a4 = -y(ihe4) * y(in14) * rate(irnag) * 1.5d0 ! n14 + 1.5 alpha => ne20
         qray(ihe4) =  qray(ihe4) + a1 + a2 + a3 + a4


   ! (a,g) and (g,a) reactions

         a1  = -0.5d0 * y(ihe4) * y(ihe4) * y(ihe4) * rate(ir3a) 
         a2  =  3.0d0 * y(ic12) * rate(irg3a) 
         a3  = -y(ihe4) * y(ic12) * rate(ircag) 
         a4  =  y(io16) * rate(iroga) 
         a5  = -y(ihe4) * y(io16) * rate(iroag) 
         a6  =  y(ine20) * rate(irnega) 
         a7  = -y(ihe4) * y(ine20) * rate(irneag) 
         a8  =  y(img24) * rate(irmgga) 
         a9  = -y(ihe4) * y(img24)* rate(irmgag) 
         a10 =  y(isi28) * rate(irsiga) 
         a11 = -y(ihe4) * y(isi28)*rate(irsiag) 
         a12 =  y(is32) * rate(irsga)

         qray(ihe4) =  qray(ihe4) + &
            a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12
               
         a1  = -y(ihe4) * y(is32) * rate(irsag) 
         a2  =  y(iar36) * rate(irarga) 
         a3  = -y(ihe4) * y(iar36)*rate(irarag) 
         a4  =  y(ica40) * rate(ircaga) 
         a5  = -y(ihe4) * y(ica40)*rate(ircaag) 
         a6  =  y(iti44) * rate(irtiga) 
         a7  = -y(ihe4) * y(iti44)*rate(irtiag) 
         a8  =  y(icr48) * rate(ircrga) 
         a9  = -y(ihe4) * y(icr48)*rate(ircrag) 
         a10 =  y(ife52) * rate(irfega) 
         a11 = -y(ihe4) * y(ife52) * rate(irfeag) 
         a12 =  y(ini56) * rate(irniga)

         qray(ihe4) =  qray(ihe4) + &
            a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12


   ! (a,p)(p,g) and (g,p)(p,a) reactions

         if (.not.deriva) then
         a1  =  0.34d0*0.5d0*y(io16)*y(io16)*rate(irs1)*rate(ir1616) 
         a2  = -y(ihe4) * y(img24) * rate(irmgap)*(1.0d0-rate(irr1))
         a3  =  y(isi28) * rate(irsigp) * rate(irr1) 
         a4  = -y(ihe4) * y(isi28) * rate(irsiap)*(1.0d0-rate(irs1)) 
         a5  =  y(is32) * rate(irsgp) * rate(irs1) 
         a6  = -y(ihe4) * y(is32) * rate(irsap)*(1.0d0-rate(irt1)) 
         a7  =  y(iar36) * rate(irargp) * rate(irt1) 
         a8  = -y(ihe4) * y(iar36) * rate(irarap)*(1.0d0-rate(iru1)) 
         a9  =  y(ica40) * rate(ircagp) * rate(iru1) 
         a10 = -y(ihe4) * y(ica40) * rate(ircaap)*(1.0d0-rate(irv1)) 
         a11 =  y(iti44) * rate(irtigp) * rate(irv1)
         a12 = -y(ihe4) * y(iti44) * rate(irtiap)*(1.0d0-rate(irw1)) 
         a13 =  y(icr48) * rate(ircrgp) * rate(irw1) 
         a14 = -y(ihe4) * y(icr48) * rate(ircrap)*(1.0d0-rate(irx1)) 
         a15 =  y(ife52) * rate(irfegp) * rate(irx1) 

         qray(ihe4) = qray(ihe4) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + &
            a8 + a9 + a10 + a11 + a12 + a13 + a14 + a15

         else
         a1  =  0.34d0*0.5d0*y(io16)*y(io16) * ratdum(irs1)*rate(ir1616) 
         a2  =  0.34d0*0.5d0*y(io16)*y(io16) * rate(irs1) * ratdum(ir1616)
         a3  = -y(ihe4)*y(img24) * rate(irmgap)*(1.0d0 - ratdum(irr1)) 
         a4  =  y(ihe4)*y(img24) * ratdum(irmgap)*rate(irr1)
         a5  =  y(isi28) * ratdum(irsigp) * rate(irr1) 
         a6  =  y(isi28) * rate(irsigp) * ratdum(irr1)
         a7  = -y(ihe4)*y(isi28) * rate(irsiap)*(1.0d0 - ratdum(irs1)) 
         a8  =  y(ihe4)*y(isi28) * ratdum(irsiap) * rate(irs1)
         a9  =  y(is32) * ratdum(irsgp) * rate(irs1)
         a10 =  y(is32) * rate(irsgp) * ratdum(irs1)

         qray(ihe4) =  qray(ihe4) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10

         a1  = -y(ihe4)*y(is32) * rate(irsap)*(1.0d0 - ratdum(irt1)) 
         a2  =  y(ihe4)*y(is32) * ratdum(irsap)*rate(irt1)
         a3  =  y(iar36) * ratdum(irargp) * rate(irt1) 
         a4  =  y(iar36) * rate(irargp) * ratdum(irt1)
         a5  = -y(ihe4)*y(iar36) * rate(irarap)*(1.0d0 - ratdum(iru1))
         a6  =  y(ihe4)*y(iar36) * ratdum(irarap)*rate(iru1)
         a7  =  y(ica40) * ratdum(ircagp) * rate(iru1)
         a8  =  y(ica40) * rate(ircagp) * ratdum(iru1)
         a9  = -y(ihe4)*y(ica40) * rate(ircaap)*(1.0d0-ratdum (irv1)) 
         a10 =  y(ihe4)*y(ica40) * ratdum(ircaap)*rate(irv1)
         a11 =  y(iti44) * ratdum(irtigp) * rate(irv1)
         a12 =  y(iti44) * rate(irtigp) * ratdum(irv1)

         qray(ihe4) =  qray(ihe4) + &
            a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11 + a12

         a1  = -y(ihe4)*y(iti44) * rate(irtiap)*(1.0d0 - ratdum(irw1))
         a2  =  y(ihe4)*y(iti44) * ratdum(irtiap)*rate(irw1)
         a3  =  y(icr48) * ratdum(ircrgp) * rate(irw1) 
         a4  =  y(icr48) * rate(ircrgp) * ratdum(irw1)
         a5  = -y(ihe4)*y(icr48) * rate(ircrap)*(1.0d0 - ratdum(irx1))
         a6  =  y(ihe4)*y(icr48) * ratdum(ircrap)*rate(irx1)
         a7  =  y(ife52) * ratdum(irfegp) * rate(irx1) 
         a8  =  y(ife52) * rate(irfegp) * ratdum(irx1)

         qray(ihe4) = qray(ihe4) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8
         end if


   ! photodisintegration reactions
         a1 =  y(ife54) * y(iprot) * y(iprot) * rate(ir5f54) 
         a2 = -y(ife52) * y(ihe4) * rate(ir6f54) 
         a3 = -y(ife52) * y(ihe4) * y(iprot) * rate(ir7f54) 
         a4 =  y(ini56) * y(iprot) * rate(ir8f54) 
         a5 = -y(ihe4) * rate(iralf1) 
         a6 =  y(ineut)*y(ineut) * y(iprot)*y(iprot) * rate(iralf2) 
         a7 =  y(ife56) * y(iprot) * y(iprot) * rate(irfe56_aux3) 
         a8 = -y(ife54) * y(ihe4) * rate(irfe56_aux4) 

         qray(ihe4) =  qray(ihe4) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8


   ! ppchain
         a1 = 0.5d0 * y(ihe3) * y(ihe3) * rate(ir33)
         a2 = y(ihe3) * y(ihe4) * rate(irhe3ag)

         qray(ihe4) =  qray(ihe4) + a1 + a2 


   ! cno cycles
         a1 = y(io16) * y(ih1) * rate(iropg) 

         qray(ihe4) =  qray(ihe4) + a1 + a2
         
         if (.not. deriva) then
            a1 = y(in14) * y(ih1) * rate(ifa) * rate(irnpg) 
            qray(ihe4) =  qray(ihe4) + a1
         else
            a1 = y(in14) * y(ih1) * rate(ifa) * ratdum(irnpg) 
            a2 = y(in14) * y(ih1) * ratdum(ifa) * rate(irnpg) 
            qray(ihe4) =  qray(ihe4) + a1 + a2
         end if


   ! c12 reactions
         a1 = -y(ic12) * y(ic12) * rate(ir1212) 
         a2 = -y(ic12) * y(io16) * rate(ir1216) 
         a3 =  (1d0/6d0) * y(ihe4) * y(ihe4) * y(ihe4) * rate(ir3a) 
         a4 = -y(ic12) * rate(irg3a) 
         a5 = -y(ic12) * y(ihe4) * rate(ircag) 
         a6 =  y(io16) * rate(iroga) 
         a7 = -y(ic12) * y(ih1) * rate(ircpg) 

         qray(ic12) =  qray(ic12) + a1 + a2 + a3 + a4 + a5 + a6 + a7
               
         if (.not. deriva) then
            a1 =  y(in14) * y(ih1) * rate(ifa) * rate(irnpg)
            qray(ic12) =  qray(ic12) + a1
         else
            a1 =  y(in14) * y(ih1) * rate(ifa) * ratdum(irnpg)
            a2 =  y(in14) * y(ih1) * ratdum(ifa) * rate(irnpg)
            qray(ic12) =  qray(ic12) + a1 + a2         
         end if


   ! n14 reactions
         a1 =  y(ic12) * y(ih1) * rate(ircpg) 
         a2 = -y(in14) * y(ih1) * rate(irnpg) 
         a3 =  y(io16) * y(ih1) * rate(iropg) 
         a4 = -y(ihe4) * y(in14) * rate(irnag) ! n14 + 1.5 alpha => ne20

         qray(in14) =  qray(in14) + a1 + a2 + a3 + a4 


   ! o16 reactions
         a1 = -y(ic12) * y(io16) * rate(ir1216) 
         a2 = -y(io16) * y(io16) * rate(ir1616) 
         a3 =  y(ic12) * y(ihe4) * rate(ircag) 
         a4 = -y(io16) * y(ihe4) * rate(iroag) 
         a5 = -y(io16) * rate(iroga) 
         a6 =  y(ine20) * rate(irnega) 
         a7 = -y(io16) * y(ih1) * rate(iropg)

         qray(io16) =  qray(io16) + a1 + a2 + a3 + a4 + a5 + a6 + a7
         
         if (.not. deriva) then
            a1 =  y(in14) * y(ih1) * rate(ifg) * rate(irnpg) 
            qray(io16) =  qray(io16) + a1
         else
            a1 =  y(in14) * y(ih1) * rate(ifg) * ratdum(irnpg) 
            a2 =  y(in14) * y(ih1) * ratdum(ifg) * rate(irnpg) 
            qray(io16) =  qray(io16) + a1 + a2
         end if


   ! ne20 reactions
         a1 =  0.5d0 * y(ic12) * y(ic12) * rate(ir1212) 
         a2 =  y(io16) * y(ihe4) * rate(iroag) 
         a3 = -y(ine20) * y(ihe4) * rate(irneag) 
         a4 = -y(ine20) * rate(irnega) 
         a5 =  y(img24) * rate(irmgga) 
         a6 =  y(in14) * y(ihe4) * rate(irnag) ! n14 + 1.5 alpha => ne20

         qray(ine20) =  qray(ine20) + a1 + a2 + a3 + a4 + a5 + a6


   ! mg24 reactions
         a1 =  0.5d0 * y(ic12) * y(io16) * rate(ir1216) 
         a2 =  y(ine20) * y(ihe4) * rate(irneag) 
         a3 = -y(img24) * y(ihe4) * rate(irmgag) 
         a4 = -y(img24) * rate(irmgga) 
         a5 =  y(isi28) * rate(irsiga)
         
         qray(img24) =  qray(img24) + a1 + a2 + a3 + a4 + a5 

         if (.not.deriva) then
         a1 = -y(img24) * y(ihe4) * rate(irmgap)*(1.0d0-rate(irr1)) 
         a2 =  y(isi28) * rate(irr1) * rate(irsigp)

         qray(img24) =  qray(img24) + a1 + a2

         else
         a1 = -y(img24)*y(ihe4) * rate(irmgap)*(1.0d0 - ratdum(irr1))
         a2 =  y(img24)*y(ihe4) * ratdum(irmgap)*rate(irr1)
         a3 =  y(isi28) * ratdum(irr1) * rate(irsigp) 
         a4 =  y(isi28) * rate(irr1) * ratdum(irsigp)

         qray(img24) =  qray(img24) + a1 + a2 + a3 + a4
         end if


   ! si28 reactions
         a1 =  0.5d0 * y(ic12) * y(io16) * rate(ir1216) 
         a2 =  0.56d0 * 0.5d0*y(io16) * y(io16) * rate(ir1616) 
         a3 =  y(img24) * y(ihe4) * rate(irmgag) 
         a4 = -y(isi28) * y(ihe4) * rate(irsiag) 
         a5 = -y(isi28) * rate(irsiga) 
         a6 =  y(is32) * rate(irsga)

         qray(isi28) =  qray(isi28) + a1 + a2 + a3 + a4 + a5 + a6

         if (.not.deriva) then
         
         a1 =  0.34d0*0.5d0*y(io16)*y(io16)*rate(irs1)*rate(ir1616) 
         a2 =  y(img24) * y(ihe4) * rate(irmgap)*(1.0d0-rate(irr1)) 
         a3 = -y(isi28) * rate(irr1) * rate(irsigp) 
         a4 = -y(isi28) * y(ihe4) * rate(irsiap)*(1.0d0-rate(irs1)) 
         a5 =  y(is32) * rate(irs1) * rate(irsgp)

         qray(isi28) =  qray(isi28) + a1 + a2 + a3 + a4 + a5

         else
         a1  =  0.34d0*0.5d0*y(io16)*y(io16) * ratdum(irs1)*rate(ir1616) 
         a2  =  0.34d0*0.5d0*y(io16)*y(io16) * rate(irs1)*ratdum(ir1616)
         a3  =  y(img24)*y(ihe4) * rate(irmgap)*(1.0d0 - ratdum(irr1))
         a4  = -y(img24)*y(ihe4) * ratdum(irmgap)*rate(irr1)
         a5  = -y(isi28) * ratdum(irr1) * rate(irsigp) 
         a6  = -y(isi28) * rate(irr1) * ratdum(irsigp)
         a7  = -y(isi28)*y(ihe4) * rate(irsiap)*(1.0d0 - ratdum(irs1))
         a8  =  y(isi28)*y(ihe4) * ratdum(irsiap)*rate(irs1)
         a9  = y(is32) * ratdum(irs1) * rate(irsgp)
         a10 = y(is32) * rate(irs1) * ratdum(irsgp)

         qray(isi28) =  qray(isi28) + &
            a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10
         end if


   ! s32 reactions
         a1 =  0.1d0 * 0.5d0*y(io16) * y(io16) * rate(ir1616) 
         a2 =  y(isi28) * y(ihe4) * rate(irsiag) 
         a3 = -y(is32) * y(ihe4) * rate(irsag) 
         a4 = -y(is32) * rate(irsga) 
         a5 =  y(iar36) * rate(irarga)

         qray(is32) =  qray(is32) + a1 + a2 + a3 + a4 + a5

         if (.not.deriva) then

         a1 =  0.34d0*0.5d0*y(io16)*y(io16)* rate(ir1616)*(1.0d0-rate(irs1)) 
         a2 =  y(isi28) * y(ihe4) * rate(irsiap)*(1.0d0-rate(irs1)) 
         a3 = -y(is32) * rate(irs1) * rate(irsgp) 
         a4 = -y(is32) * y(ihe4) * rate(irsap)*(1.0d0-rate(irt1)) 
         a5 =  y(iar36) * rate(irt1) * rate(irargp)

         qray(is32) =  qray(is32) + a1 + a2 + a3 + a4 + a5

         else
         a1  =  0.34d0*0.5d0*y(io16)*y(io16) * rate(ir1616)*(1.0d0-ratdum(irs1))
         a2  = -0.34d0*0.5d0*y(io16)*y(io16) * ratdum(ir1616)*rate(irs1)
         a3  =  y(isi28)*y(ihe4) * rate(irsiap)*(1.0d0-ratdum(irs1))
         a4  = -y(isi28)*y(ihe4) * ratdum(irsiap)*rate(irs1) 
         a5  = -y(is32) * ratdum(irs1) * rate(irsgp) 
         a6  = -y(is32) * rate(irs1) * ratdum(irsgp)
         a7  = -y(is32)*y(ihe4) * rate(irsap)*(1.0d0-ratdum(irt1))
         a8  =  y(is32)*y(ihe4) * ratdum(irsap)*rate(irt1)
         a9  =  y(iar36) * ratdum(irt1) * rate(irargp)
         a10 =  y(iar36) * rate(irt1) * ratdum(irargp)

         qray(is32) =  qray(is32) + &
            a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10
         end if


   ! ar36 reactions
         a1 =  y(is32) * y(ihe4) * rate(irsag) 
         a2 = -y(iar36) * y(ihe4) * rate(irarag)
         a3 = -y(iar36) * rate(irarga) 
         a4 =  y(ica40) * rate(ircaga)

         qray(iar36) =  qray(iar36) + a1 + a2 + a3 + a4

         if (.not.deriva) then
         a1 = y(is32) * y(ihe4) * rate(irsap)*(1.0d0-rate(irt1)) 
         a2 = -y(iar36) * rate(irt1) * rate(irargp) 
         a3 = -y(iar36) * y(ihe4) * rate(irarap)*(1.0d0-rate(iru1)) 
         a4 =  y(ica40) * rate(ircagp) * rate(iru1)

         qray(iar36) =  qray(iar36) + a1 + a2 + a3 + a4

         else
         a1 =  y(is32)*y(ihe4) * rate(irsap)*(1.0d0 - ratdum(irt1)) 
         a2 = -y(is32)*y(ihe4) * ratdum(irsap)*rate(irt1)
         a3 = -y(iar36) * ratdum(irt1) * rate(irargp) 
         a4 = -y(iar36) * rate(irt1) * ratdum(irargp)
         a5 = -y(iar36)*y(ihe4) * rate(irarap)*(1.0d0-ratdum(iru1))
         a6 =  y(iar36)*y(ihe4) * ratdum(irarap)*rate(iru1)
         a7 =  y(ica40) * ratdum(ircagp) * rate(iru1) 
         a8 =  y(ica40) * rate(ircagp) * ratdum(iru1)

         qray(iar36) =  qray(iar36) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 
         end if


   ! ca40 reactions
         a1 =  y(iar36) * y(ihe4) * rate(irarag)
         a2 = -y(ica40) * y(ihe4) * rate(ircaag)
         a3 = -y(ica40) * rate(ircaga)
         a4 =  y(iti44) * rate(irtiga)

         qray(ica40) =  qray(ica40) + a1 + a2 + a3 + a4

         if (.not.deriva) then

         a1 =  y(iar36) * y(ihe4) * rate(irarap)*(1.0d0-rate(iru1)) 
         a2 = -y(ica40) * rate(ircagp) * rate(iru1) 
         a3 = -y(ica40) * y(ihe4) * rate(ircaap)*(1.0d0-rate(irv1)) 
         a4 =  y(iti44) * rate(irtigp) * rate(irv1)

         qray(ica40) =  qray(ica40) + a1 + a2 + a3 + a4

         else
         a1 =  y(iar36)*y(ihe4) * rate(irarap)*(1.0d0-ratdum(iru1))
         a2 = -y(iar36)*y(ihe4) * ratdum(irarap)*rate(iru1)
         a3 = -y(ica40) * ratdum(ircagp) * rate(iru1) 
         a4 = -y(ica40) * rate(ircagp) * ratdum(iru1)
         a5 = -y(ica40)*y(ihe4) * rate(ircaap)*(1.0d0-ratdum(irv1)) 
         a6 =  y(ica40)*y(ihe4) * ratdum(ircaap)*rate(irv1)
         a7 =  y(iti44) * ratdum(irtigp) * rate(irv1) 
         a8 =  y(iti44) * rate(irtigp) * ratdum(irv1)

         qray(ica40) =  qray(ica40) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 
         end if


   ! ti44 reactions
         a1 =  y(ica40) * y(ihe4) * rate(ircaag) 
         a2 = -y(iti44) * y(ihe4) * rate(irtiag) 
         a3 = -y(iti44) * rate(irtiga) 
         a4 =  y(icr48) * rate(ircrga)

         qray(iti44) =  qray(iti44) + a1 + a2 + a3 + a4

         if (.not.deriva) then
         a1 =  y(ica40) * y(ihe4) * rate(ircaap)*(1.0d0-rate(irv1)) 
         a2 = -y(iti44) * rate(irv1) * rate(irtigp) 
         a3 = -y(iti44) * y(ihe4) * rate(irtiap)*(1.0d0-rate(irw1)) 
         a4 =  y(icr48) * rate(irw1) * rate(ircrgp)

         qray(iti44) =  qray(iti44) + a1 + a2 + a3 + a4

         else
         a1 =  y(ica40)*y(ihe4) * rate(ircaap)*(1.0d0-ratdum(irv1)) 
         a2 = -y(ica40)*y(ihe4) * ratdum(ircaap)*rate(irv1)
         a3 = -y(iti44) * ratdum(irv1) * rate(irtigp) 
         a4 = -y(iti44) * rate(irv1) * ratdum(irtigp)
         a5 = -y(iti44)*y(ihe4) * rate(irtiap)*(1.0d0-ratdum(irw1))
         a6 =  y(iti44)*y(ihe4) * ratdum(irtiap)*rate(irw1)
         a7 =  y(icr48) * ratdum(irw1) * rate(ircrgp)
         a8 =  y(icr48) * rate(irw1) * ratdum(ircrgp)

         qray(iti44) =  qray(iti44) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 
         end if


   ! cr48 reactions
         a1 =  y(iti44) * y(ihe4) * rate(irtiag) 
         a2 = -y(icr48) * y(ihe4) * rate(ircrag) 
         a3 = -y(icr48) * rate(ircrga) 
         a4 =  y(ife52) * rate(irfega)

         qray(icr48) =  qray(icr48) + a1 + a2 + a3 + a4

         if (.not.deriva) then
         a1 =  y(iti44) * y(ihe4) * rate(irtiap)*(1.0d0-rate(irw1)) 
         a2 = -y(icr48) * rate(irw1) * rate(ircrgp) 
         a3 = -y(icr48) * y(ihe4) * rate(ircrap)*(1.0d0-rate(irx1)) 
         a4 =  y(ife52) * rate(irx1) * rate(irfegp)

         qray(icr48) =  qray(icr48) + a1 + a2 + a3 + a4

         else
         a1 =  y(iti44)*y(ihe4) * rate(irtiap)*(1.0d0-ratdum(irw1)) 
         a2 = -y(iti44)*y(ihe4) * ratdum(irtiap)*rate(irw1)
         a3 = -y(icr48) * ratdum(irw1) * rate(ircrgp) 
         a4 = -y(icr48) * rate(irw1) * ratdum(ircrgp)
         a5 = -y(icr48)*y(ihe4) * rate(ircrap)*(1.0d0-ratdum(irx1))
         a6 =  y(icr48)*y(ihe4) * ratdum(ircrap)*rate(irx1)
         a7 =  y(ife52) * ratdum(irx1) * rate(irfegp) 
         a8 =  y(ife52) * rate(irx1) * ratdum(irfegp)

         qray(icr48) =  qray(icr48) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 
         end if


   ! crx reactions
         a1 = y(ife56) * fe56ec_fake_factor * rate(irn56ec)
         
         qray(icrx) = qray(icrx) + a1

   ! fe52 reactions
         a1 =  y(icr48) * y(ihe4) * rate(ircrag) 
         a2 = -y(ife52) * y(ihe4) * rate(irfeag) 
         a3 = -y(ife52) * rate(irfega) 
         a4 =  y(ini56) * rate(irniga)

         qray(ife52) =  qray(ife52) + a1 + a2 + a3 + a4

         if (.not.deriva) then
         a1 =  y(icr48) * y(ihe4) * rate(ircrap)*(1.0d0-rate(irx1)) 
         a2 = -y(ife52) * rate(irx1) * rate(irfegp) 

         qray(ife52) =  qray(ife52) + a1 + a2

         else
         a1 =  y(icr48)*y(ihe4) * rate(ircrap)*(1.0d0-ratdum(irx1))
         a2 = -y(icr48)*y(ihe4) * ratdum(ircrap)*rate(irx1)
         a3 = -y(ife52) * ratdum(irx1) * rate(irfegp) 
         a4 = -y(ife52) * rate(irx1) * ratdum(irfegp)

         qray(ife52) =  qray(ife52) + a1 + a2 + a3 + a4
         end if

         a1 =  y(ife54) * rate(ir1f54) 
         a2 = -y(ife52) * y(ineut) * y(ineut) * rate(ir2f54) 
         a3 =  y(ife54) * y(iprot) * y(iprot) * rate(ir5f54) 
         a4 = -y(ife52) * y(ihe4) * rate(ir6f54) 
         a5 = -y(ife52) * y(ihe4) * y(iprot) * rate(ir7f54) 
         a6 =  y(ini56) * y(iprot) * rate(ir8f54) 

         qray(ife52) =  qray(ife52) + a1 + a2 + a3 + a4 + a5 + a6 


   ! fe54 reactions
         a1  = -y(ife54) * rate(ir1f54)
         a2  =  y(ife52) * y(ineut) * y(ineut) * rate(ir2f54) 
         a3  = -y(ife54) * y(iprot) * y(iprot) * rate(ir3f54) 
         a4  =  y(ini56) * rate(ir4f54) 
         a5  = -y(ife54) * y(iprot) * y(iprot) * rate(ir5f54) 
         a6  =  y(ife52) * y(ihe4) * rate(ir6f54) 
         a7  =  y(ife56) * rate(irfe56_aux1) 
         a8  = -y(ife54) * y(ineut) * y(ineut) * rate(irfe56_aux2) 
         a9  =  y(ife56) * y(iprot) * y(iprot) * rate(irfe56_aux3) 
         a10 = -y(ife54) * y(ihe4) * rate(irfe56_aux4) 

         qray(ife54) =  qray(ife54) + &
            a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10


   ! fe56 reactions
         if (plus_co56) then
            a1 =  y(ico56) * rate(irco56ec) 
         else
            a1 =  y(ini56) * rate(irn56ec)  
         end if
         a2 = -y(ife56) * fe56ec_fake_factor * rate(irn56ec) 
         a3 = -y(ife56) * rate(irfe56_aux1) 
         a4 =  y(ife54) * y(ineut) * y(ineut) * rate(irfe56_aux2)  
         a5 = -y(ife56) * y(iprot) * y(iprot) * rate(irfe56_aux3) 
         a6 =  y(ife54) * y(ihe4) * rate(irfe56_aux4) 

         qray(ife56) =  qray(ife56) + a1 + a2 + a3 + a4 + a5 + a6 

         if (plus_co56) then
      ! co56 reactions
            a1 =  y(ini56) * rate(irn56ec)  
            a2 = -y(ico56) * rate(irco56ec) 
            
            qray(ico56) =  qray(ico56) + a1 + a2
         end if

   ! ni56 reactions
         a1 =  y(ife52) * y(ihe4) * rate(irfeag) 
         a2 = -y(ini56) * rate(irniga) 
         a3 = -y(ini56) * rate(irn56ec) 
         
         qray(ini56) =  qray(ini56) + a1 + a2 + a3

         a1 =  y(ife54) * y(iprot) * y(iprot) * rate(ir3f54) 
         a2 = -y(ini56) * rate(ir4f54) 
         a3 =  y(ife52) * y(ihe4)* y(iprot) * rate(ir7f54) 
         a4 = -y(ini56) * y(iprot) * rate(ir8f54)

         qray(ini56) =  qray(ini56) + a1 + a2 + a3 + a4

   ! neutrons
         a1 =  2.0d0 * y(ife54) * rate(ir1f54) 
         a2 = -2.0d0 * y(ife52) * y(ineut) * y(ineut) * rate(ir2f54) 
         a3 =  2.0d0 * y(ihe4) * rate(iralf1) 
         a4 = -2.0d0 * y(ineut)*y(ineut) * y(iprot)*y(iprot) * rate(iralf2) 
         a5 =  y(iprot) * rate(irpen) 
         a6 = -y(ineut) * rate(irnep) 
         a7 =  2.0d0 * y(ife56) * rate(irfe56_aux1)
         a8 = -2.0d0 * y(ife54) * y(ineut) * y(ineut) * rate(irfe56_aux2)
         a9 = -fe56ec_n_neut * y(ife56) * fe56ec_fake_factor * rate(irn56ec)

         qray(ineut) =  qray(ineut) + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9

   ! photodisintegration protons
         a1  = -2.0d0 * y(ife54) * y(iprot) * y(iprot) * rate(ir3f54) 
         a2  =  2.0d0 * y(ini56) * rate(ir4f54) 
         a3  = -2.0d0 * y(ife54) * y(iprot) * y(iprot) * rate(ir5f54) 
         a4  =  2.0d0 * y(ife52) * y(ihe4) * rate(ir6f54) 
         a5  =  2.0d0 * y(ihe4) * rate(iralf1) 
         a6  = -2.0d0 * y(ineut)*y(ineut) * y(iprot)*y(iprot) * rate(iralf2) 
         a7  = -y(iprot) * rate(irpen) 
         a8  =  y(ineut) * rate(irnep)  
         a9  = -2.0d0 * y(ife56) * y(iprot) * y(iprot) * rate(irfe56_aux3) 
         a10 =  2.0d0 * y(ife54) * y(ihe4) * rate(irfe56_aux4)

         qray(iprot) =  qray(iprot) + &
            a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10

   ! now set the real(dp) return argument dydt
         okay = .true.
         do i=1,species(plus_co56)
            dydt(i) = qray(i)
            if (is_bad(dydt(i))) then
               write(*,*) 'dydt(i)', i, dydt(i), y(i)
               okay = .false.
            end if
         end do
         if (.not. okay) then
            do i=1,num_reactions(plus_co56)
               write(*,*) trim(ratnam(i)), i, rate(i)
            end do
            stop 'approx21_dydt'
         end if

         end subroutine approx21_dydt
         
         
         real(dp) function approx21_eval_PPII_fraction(y, rate) result(fII)
            real(dp), dimension(:), intent(in) :: y, rate
            real(dp) :: rateII, rateIII, rsum
            include 'formats'
            rateII  = rate(ir_be7_wk_li7)
            rateIII = y(ih1) * rate(ir_be7_pg_b8)
            rsum = rateII + rateIII
            if (rsum < 1d-50) then
               fII = 0.5d0
            else
               fII = rateII / rsum
            end if               
         end function approx21_eval_PPII_fraction
               

         subroutine approx21_eps_info( &
               n, y, mion, dydt, rate, fII, &
               Qtotal_rpp, Qneu_rpp, &
               Qr33, &
               Qtotal_rpp2, Qneu_rpp2, &
               Qtotal_rpp3, Qneu_rpp3, &
               Qtotal_rcpg, Qneu_rcpg, &
               Qtotal_rnpg, Qneu_rnpg, &
               Qtotal_ropg, Qneu_ropg, &
               Qrn14_to_o16, &
               Qtotal_rpen, Qneu_rpen, &
               Qtotal_rnep, Qneu_rnep, &
               Qtotal_rni56ec, Qneu_rni56ec, &
               Qtotal_rco56ec, Qneu_rco56ec, &
               Qtotal_rfe56ec, Qneu_rfe56ec, &
               Qrn14ag, &
               Qr3alf, &
               Qrc12ag, Qro16ag,&
               Qr1212, &
               Qr1216_to_mg24, Qr1216_to_si28, &
               Qr1616a, Qr1616g, &
               Qrne20ag, &
               Qrmg24ag, &
               Qrsi28ag, &
               Qrs32ag, &
               Qrar36ag, &
               Qrca40ag, &
               Qrti44ag, &
               Qrcr48ag, &
               Qrfe52ag, &
               Qrfe52ng, & 
               Qrfe53ng, & 
               Qrfe54ng, & 
               Qrfe55ng, & 
               Qrfe52neut_to_fe54, &
               Qrfe52aprot_to_fe54, &
               Qrfe54ng_to_fe56, &
               Qrfe54aprot_to_fe56, &
               Qrfe52aprot_to_ni56, &
               Qrfe54prot_to_ni56, &     
               Qrhe4_breakup, &
               Qrhe4_rebuild, &      
               eps_total, eps_neu, &
               do_eps_nuc_categories, eps_nuc_categories, &
               dbg, &
               plus_co56, &
               ierr)
            use const_def, only: Qconv
            use chem_def, only: num_categories, category_name, chem_isos, &
               ipp, icno, i3alf, i_burn_c, i_burn_n, i_burn_o, i_burn_ne, i_burn_na, &
               i_burn_mg, i_burn_si, i_burn_s, i_burn_ar, i_burn_ca, i_burn_ti, i_burn_cr, &
               i_burn_fe, icc, ico, ioo, ipnhe4, iphoto, i_ni56_co56, i_co56_fe56, iother
            use net_def, only: Net_Info
            type (Net_Info), pointer :: n
            real(dp), dimension(:), intent(in) :: y, mion, dydt, rate
            real(dp), intent(in) :: fII, &
               Qtotal_rpp, Qneu_rpp, Qr33, &
               Qtotal_rpp2, Qneu_rpp2, &
               Qtotal_rpp3, Qneu_rpp3, &
               Qtotal_rcpg, Qneu_rcpg, &
               Qtotal_rnpg, Qneu_rnpg, &
               Qtotal_ropg, Qneu_ropg, &
               Qrn14_to_o16, Qtotal_rpen, Qneu_rpen, &
               Qtotal_rnep, Qneu_rnep, &
               Qtotal_rni56ec, Qneu_rni56ec, &
               Qtotal_rco56ec, Qneu_rco56ec, &
               Qtotal_rfe56ec, Qneu_rfe56ec, &
               Qrn14ag, &
               Qr3alf, &
               Qrc12ag, Qro16ag, &
               Qr1212, &
               Qr1216_to_mg24, Qr1216_to_si28, &
               Qr1616a, Qr1616g, &
               Qrne20ag, &
               Qrmg24ag, &
               Qrsi28ag, &
               Qrs32ag, &
               Qrar36ag, &
               Qrca40ag, &
               Qrti44ag, &
               Qrcr48ag, &
               Qrfe52ag, & 
               Qrfe52ng, & 
               Qrfe53ng, & 
               Qrfe54ng, & 
               Qrfe55ng, & 
               Qrfe52neut_to_fe54, &
               Qrfe52aprot_to_fe54, &
               Qrfe54ng_to_fe56, &
               Qrfe54aprot_to_fe56, &
               Qrfe52aprot_to_ni56, &
               Qrfe54prot_to_ni56, &           
               Qrhe4_breakup, &
               Qrhe4_rebuild
            logical, intent(in) :: do_eps_nuc_categories
            real(dp), intent(out) :: eps_total, eps_neu, eps_nuc_categories(:)
            logical, intent(in) :: dbg
            integer, intent(out) :: ierr

            integer :: i, fe56ec_n_neut
            real(qp) :: a1, a2, xx, eps_neu_q, eps_nuc_cat(num_categories), eps_total_q, &
               eps_nuc_q, sum_categories_q
            real(dp) :: enuc_conv2, sum_categories, eps_nuc, fe56ec_fake_factor
            logical, intent(in) ::  plus_co56
            
            include 'formats'
               
            !write(*,1) 'reaction_Qs(irn14_to_o16) Qrn14_to_o16*Qconv', Qrn14_to_o16*Qconv
            
            ierr = 0

            xx = 0.0_qp
            do i=1,species(plus_co56)
               a1 = dydt(i) 
               a2 = mion(i)
               xx = xx + a1*a2
            end do
            eps_total_q = -m3(avo,clight,clight) * xx
            eps_total = eps_total_q
            
            fe56ec_fake_factor = eval_fe56ec_fake_factor( &
               n% g% fe56ec_fake_factor, n% g% min_T_for_fe56ec_fake_factor, n% temp)
            fe56ec_n_neut = n% g% fe56ec_n_neut
            
            eps_neu_q = &
               m5(Qneu_rpp, 0.5d0, y(ih1), y(ih1), rate(irpp)) + &
               m5(Qneu_rpp2, y(ihe3), y(ihe4), rate(irhe3ag), fII) + &
               m5(Qneu_rpp3, y(ihe3), y(ihe4), rate(irhe3ag), (1d0-fII)) + &
               m4(Qneu_rcpg, y(ic12), y(ih1), rate(ircpg)) + &
               m5(Qneu_rnpg, y(in14), y(ih1), rate(irnpg), rate(ifa)) + &
               m4(Qneu_ropg, y(io16), y(ih1), rate(iropg)) + &
               m3(Qneu_rpen, y(ih1), rate(irpen)) + &
               m3(Qneu_rpen, y(iprot), rate(irpen)) + &
               m3(Qneu_rnep, y(ineut), rate(irnep)) + &
               m4(Qneu_rfe56ec, y(ife56), fe56ec_fake_factor, rate(irn56ec))

            if (plus_co56) then
               eps_neu_q = eps_neu_q + &
                  m3(Qneu_rni56ec, y(ini56), rate(irn56ec)) + &
                  m3(Qneu_rco56ec, y(ico56), rate(irco56ec))
            else
               eps_neu_q = eps_neu_q + &
                  m3(Qneu_rni56ec + Qneu_rco56ec, y(ini56), rate(irn56ec))
            end if
            eps_neu_q = eps_neu_q * Qconv
            eps_neu = eps_neu_q

            eps_nuc_q = eps_total_q - eps_neu_q
            eps_nuc = eps_nuc_q
            
            if (.not. do_eps_nuc_categories) return

            do i=1,num_categories
               eps_nuc_cat(i) = 0d0
            end do
            
            eps_nuc_cat(ipp) = &
               m5(Qtotal_rpp - Qneu_rpp, 0.5d0, y(ih1), y(ih1), rate(irpp)) + &
               m5(Qr33, 0.5d0, y(ihe3), y(ihe3), rate(ir33)) + &
               m4(y(ihe3), y(ihe4), rate(irhe3ag), ( &
                  (Qtotal_rpp2 - Qneu_rpp2)*fII + &
                  (Qtotal_rpp3 - Qneu_rpp3)*(1d0 - fII)) )
            eps_nuc_cat(ipp) = eps_nuc_cat(ipp) * Qconv
               
            eps_nuc_cat(icno) = &
               m4(Qtotal_rcpg - Qneu_rcpg, y(ic12), y(ih1), rate(ircpg)) + &
               m5(Qtotal_rnpg - Qneu_rnpg, y(in14), y(ih1), rate(ifa), rate(irnpg)) + &
               m4(Qtotal_ropg - Qneu_ropg, y(io16), y(ih1), rate(iropg)) 
            eps_nuc_cat(icno) = eps_nuc_cat(icno) * Qconv
               
            eps_nuc_cat(i3alf) = m5((1d0/6d0)*Qr3alf, y(ihe4), y(ihe4), y(ihe4), rate(ir3a))
            eps_nuc_cat(i3alf) = eps_nuc_cat(i3alf) * Qconv
            
            eps_nuc_cat(i_burn_c) = m4(Qrc12ag, y(ic12), y(ihe4), rate(ircag))
            eps_nuc_cat(i_burn_c) = eps_nuc_cat(i_burn_c) * Qconv
               
            eps_nuc_cat(i_burn_n) = &
               m4(Qrn14ag, y(ihe4), y(in14), rate(irnag)) + &
               m5(Qrn14_to_o16, y(in14), y(ih1), rate(ifg), rate(irnpg)) 
            eps_nuc_cat(i_burn_n) = eps_nuc_cat(i_burn_n) * Qconv

            eps_nuc_cat(i_burn_o) = m4(Qro16ag, y(io16), y(ihe4), rate(iroag) )
            eps_nuc_cat(i_burn_o) = eps_nuc_cat(i_burn_o) * Qconv

            eps_nuc_cat(icc) = m5(Qr1212, 0.5d0, y(ic12), y(ic12), rate(ir1212))
            eps_nuc_cat(icc) = eps_nuc_cat(icc) * Qconv
            
            eps_nuc_cat(ico) = m4(0.5d0*(Qr1216_to_mg24 + Qr1216_to_si28), &
               y(ic12), y(io16), rate(ir1216))
            eps_nuc_cat(ico) = eps_nuc_cat(ico) * Qconv

            eps_nuc_cat(ioo) = m5(0.5d0, y(io16), y(io16), rate(ir1616), &
               ! these make he4 + si28
               Qr1616a * (0.56d0 + 0.34d0*rate(irs1)) + &
               ! these make s32
               Qr1616g * (0.1d0 + 0.34d0*(1d0 - rate(irs1))))
            eps_nuc_cat(ioo) = eps_nuc_cat(ioo) * Qconv
            
            eps_nuc_cat(i_burn_ne) = m4(Qrne20ag, y(ihe4), y(ine20), rate(irneag))
            eps_nuc_cat(i_burn_ne) = eps_nuc_cat(i_burn_ne) * Qconv
               
            eps_nuc_cat(i_burn_mg) = m4(Qrmg24ag, y(ihe4), y(img24), ( &
               rate(irmgag) + rate(irmgap)*(1.0d0-rate(irr1))) )
            eps_nuc_cat(i_burn_mg) = eps_nuc_cat(i_burn_mg) * Qconv
               
            eps_nuc_cat(i_burn_si) = m4(Qrsi28ag, y(ihe4), y(isi28), ( &
               rate(irsiag) + rate(irsiap)*(1.0d0-rate(irs1))) )
            eps_nuc_cat(i_burn_si) = eps_nuc_cat(i_burn_si) * Qconv
               
            eps_nuc_cat(i_burn_s) = m4(Qrs32ag, y(ihe4), y(is32), ( &
               rate(irsag) + rate(irsap)*(1.0d0-rate(irt1))) )
            eps_nuc_cat(i_burn_s) = eps_nuc_cat(i_burn_s) * Qconv
               
            eps_nuc_cat(i_burn_ar) = m4(Qrar36ag, y(ihe4), y(iar36), ( &
               rate(irarag) + rate(irarap)*(1.0d0-rate(iru1))) )
            eps_nuc_cat(i_burn_ar) = eps_nuc_cat(i_burn_ar) * Qconv
               
            eps_nuc_cat(i_burn_ca) = m4(Qrca40ag, y(ihe4), y(ica40), ( &
               rate(ircaag) + rate(ircaap)*(1.0d0-rate(irv1))) )
            eps_nuc_cat(i_burn_ca) = eps_nuc_cat(i_burn_ca) * Qconv
               
            eps_nuc_cat(i_burn_ti) = m4(Qrti44ag, y(ihe4), y(iti44), ( &
               rate(irtiag) + rate(irtiap)*(1.0d0-rate(irw1)) ))
            eps_nuc_cat(i_burn_ti) = eps_nuc_cat(i_burn_ti) * Qconv
               
            eps_nuc_cat(i_burn_cr) = m4(Qrcr48ag, y(ihe4), y(icr48), ( &
               rate(ircrag) + rate(ircrap)*(1.0d0-rate(irx1)) ))
            eps_nuc_cat(i_burn_cr) = eps_nuc_cat(i_burn_cr) * Qconv
               
            eps_nuc_cat(i_burn_fe) = &
               m4(Qrfe52ag, y(ihe4), y(ife52), rate(irfeag)) + &
               m5(Qrfe52aprot_to_ni56, y(ife52), y(ihe4), y(iprot), rate(ir7f54)) + &
               m5(Qrfe52neut_to_fe54, y(ife52), y(ineut), y(ineut), rate(ir2f54)) + &
               m5(Qrfe54ng_to_fe56, y(ife54), y(ineut), y(ineut), rate(irfe56_aux2)) + &
               m4(Qrfe52aprot_to_fe54, y(ife52), y(ihe4), rate(ir6f54)) + &
               m4(Qrfe54aprot_to_fe56, y(ife54), y(ihe4), rate(irfe56_aux4)) + &
               m5(Qrfe54prot_to_ni56, y(ife54), y(iprot), y(iprot), rate(ir3f54)) + &
               m4(Qtotal_rfe56ec - Qneu_rfe56ec, y(ife56), fe56ec_fake_factor, rate(irn56ec))
            eps_nuc_cat(i_burn_fe) = eps_nuc_cat(i_burn_fe) * Qconv
            
            eps_nuc_cat(ipnhe4) = &
               m4(Qrhe4_rebuild, y(ineut)*y(ineut), y(iprot)*y(iprot), rate(iralf2))
            eps_nuc_cat(ipnhe4) = eps_nuc_cat(ipnhe4) * Qconv

            eps_nuc_cat(iother) = &
               m3(Qtotal_rpen, y(ih1), rate(irpen)) + &
               m3(Qtotal_rpen, y(iprot), rate(irpen)) + &
               m3(Qtotal_rnep, y(ineut), rate(irnep))
            eps_nuc_cat(iother) = eps_nuc_cat(iother) * Qconv
               
            eps_nuc_cat(iphoto) = &
               m3(Qrhe4_breakup, y(ihe4), rate(iralf1)) - ( & ! note: Qrhe4_breakup < 0
               m3(Qrc12ag, y(io16), rate(iroga)) + & ! all the rest are > 0 Q's for forward reactions
               m3(Qr3alf, y(ic12), rate(irg3a)) + &
               m3(Qro16ag, y(ine20), rate(irnega)) + & 
               m3(Qrne20ag, y(img24), rate(irmgga)) + &
               m3(Qrmg24ag, y(isi28), rate(irsiga)) + &
               m4(Qrmg24ag, y(isi28), rate(irsigp), rate(irr1)) + &
               m3(Qrsi28ag, y(is32), rate(irsga)) + &
               m4(Qrsi28ag, y(is32), rate(irsgp), rate(irs1)) + &
               m3(Qrs32ag, y(iar36), rate(irarga)) + &
               m4(Qrs32ag, y(iar36), rate(irargp), rate(irt1)) + &
               m3(Qrar36ag, y(ica40), rate(ircaga)) + &
               m4(Qrar36ag, y(ica40), rate(ircagp), rate(iru1)) + &
               m3(Qrca40ag, y(iti44), rate(irtiga)) + &
               m4(Qrca40ag, y(iti44), rate(irtigp), rate(irv1)) + &
               m3(Qrti44ag, y(icr48), rate(ircrga)) + &
               m4(Qrti44ag, y(icr48), rate(ircrgp), rate(irw1)) + &
               m3(Qrcr48ag, y(ife52), rate(irfega)) + &
               m4(Qrcr48ag, y(ife52), rate(irfegp), rate(irx1)) + &
               m4(Qrfe52aprot_to_ni56, y(ini56), y(iprot), rate(ir8f54)) + &
               m5(Qrfe52aprot_to_fe54, y(ife54), y(iprot), y(iprot), rate(ir5f54)) + &
               m3(Qrfe52ag, y(ini56), rate(irniga)) + &
               m3(Qrfe52neut_to_fe54, y(ife54), rate(ir1f54)) + &
               m3(Qrfe54ng_to_fe56, y(ife56), rate(irfe56_aux1)) + &
               m5(Qrfe54aprot_to_fe56, y(ife56), y(iprot), y(iprot), rate(irfe56_aux3)) + &
               m3(Qrfe54prot_to_ni56, y(ini56), rate(ir4f54)))
            eps_nuc_cat(iphoto) = eps_nuc_cat(iphoto) * Qconv
            
            eps_nuc_cat(i_ni56_co56) = &
               m4(Qtotal_rni56ec - Qneu_rni56ec, y(ini56), rate(irn56ec), Qconv)

            if (plus_co56) then
               eps_nuc_cat(i_co56_fe56) = &
                  m4(Qtotal_rco56ec - Qneu_rco56ec, y(ico56), rate(irco56ec), Qconv)
            else
               eps_nuc_cat(i_co56_fe56) = &
                  m4(Qtotal_rco56ec - Qneu_rco56ec, y(ini56), rate(irn56ec), Qconv)
            end if
            
            do i=1,num_categories
               eps_nuc_categories(i) = eps_nuc_cat(i)
            end do
            
            ! check eps_nuc vs sum(eps_nuc_cat)
            
            sum_categories_q = sum(eps_nuc_cat)
            sum_categories = sum_categories_q
            
            if (.false. .and. &
               abs(eps_nuc) > 1d-10*abs(eps_nuc_cat(iphoto)) .and. abs(eps_nuc) > 1d0 .and. &
               abs(sum_categories - eps_nuc) > 1d-2*min(abs(sum_categories),abs(eps_nuc))) then
            !$OMP critical (net21_crit1)
               !write(*,*) '>>>>> problem in net_approx21_procs.inc, approx21_eps_info <<<<<<'
               !write(*,*) ' please report it.  can edit the file in eos/private to remove this test. '
               write(*,1) 'logT', n% logT
               write(*,1) 'approx21_eps_info'
               write(*,*)
               do i=1,num_categories
                  if (abs(eps_nuc_categories(i)) > 1d-6) then
                     write(*,1) trim(category_name(i)), eps_nuc_categories(i)
                  end if
               end do
               write(*,*)         
               write(*,1) 'eps_total', eps_total
               write(*,1) 'eps_neu', eps_neu
               write(*,1) 'eps_nuc', eps_nuc
               write(*,1) 'sum(cat)', sum_categories_q
               write(*,1) 'sum(cat) - eps_nuc', sum_categories_q - eps_nuc_q
               write(*,1) 'sum(cat)/eps_nuc - 1', (sum_categories_q - eps_nuc_q)/eps_nuc_q
               write(*,*)
               stop 1
            !$OMP end critical (net21_crit1)
            end if
            
            ! for debugging use reduced_net_for_testing

            if (.false. .and. n% logT >= 9.220336900d0 .and. n% logT <= 9.2203369009d0 .and. &
               abs(eps_nuc) > 1d-10*abs(eps_nuc_cat(iphoto)) .and. abs(eps_nuc) > 1d0 .and. &
               abs(sum_categories - eps_nuc) > 1d-2*min(abs(sum_categories),abs(eps_nuc))) then
               write(*,1) 'logT', n% logT
               write(*,2) 'icrx ' // trim(chem_isos% name(n% g% approx21_ye_iso)), icrx
               write(*,2) 'fe56ec_n_neut', fe56ec_n_neut
               write(*,1) 'fe56ec_fake_factor', fe56ec_fake_factor
               write(*,*)
               do i=1,num_categories
                  write(*,1) trim(category_name(i)), eps_nuc_categories(i)
               end do
               write(*,*)         
               write(*,1) 'eps_total', eps_total
               write(*,1) 'eps_neu', eps_neu
               write(*,1) 'eps_nuc', eps_nuc
               write(*,1) 'sum(cat)', sum_categories
               write(*,1) 'sum(cat) - eps_nuc', sum_categories_q - eps_nuc_q
               write(*,1) 'sum(cat)/eps_nuc - 1', (sum_categories_q - eps_nuc_q)/eps_nuc_q
               stop 'approx21_eps_info'
            end if
            
            
            contains
            
            real(qp) function m2(a1,a2)
               real(dp), intent(in) :: a1, a2
               real(qp) :: q1, q2
               q1 = a1; q2 = a2; m2 = q1*q2
            end function m2
            
            real(qp) function m3(a1,a2,a3)
               real(dp), intent(in) :: a1, a2, a3
               real(qp) :: q1, q2, q3
               q1 = a1; q2 = a2; q3 = a3; m3 = q1*q2*q3
            end function m3
            
            real(qp) function m4(a1,a2,a3,a4)
               real(dp), intent(in) :: a1, a2, a3, a4
               real(qp) :: q1, q2, q3, q4
               q1 = a1; q2 = a2; q3 = a3; q4 = a4; m4 = q1*q2*q3*q4
            end function m4
            
            real(qp) function m5(a1,a2,a3,a4,a5)
               real(dp), intent(in) :: a1, a2, a3, a4, a5
               real(qp) :: q1, q2, q3, q4, q5
               q1 = a1; q2 = a2; q3 = a3; q4 = a4; q5 = a5; m5 = q1*q2*q3*q4*q5
            end function m5

         end subroutine approx21_eps_info


         subroutine approx21_d_epsneu_dy( &
               y, rate, &
               Qneu_rpp, Qneu_rpp2, Qneu_rpp3, &
               Qneu_rcpg, Qneu_rnpg, Qneu_ropg, &
               d_epsneu_dy, plus_co56, ierr)
            use const_def, only: Qconv
            real(dp), dimension(:), intent(in) :: y, rate
            real(dp), intent(in) :: &
               Qneu_rpp, Qneu_rpp2, Qneu_rpp3, &
               Qneu_rcpg, Qneu_rnpg, Qneu_ropg
            real(dp), intent(inout) :: d_epsneu_dy(:)
            logical, intent(in) ::  plus_co56
            integer, intent(out) :: ierr
            
            real(dp) :: fII
            
            ierr = 0
            
            fII = 0.5d0 ! fix this
            
            d_epsneu_dy(1:species(plus_co56)) = 0d0
            
            d_epsneu_dy(ih1) = Qconv*( &
               Qneu_rpp * y(ih1) * rate(irpp) + & ! rpp_to_he3
               Qneu_rcpg * y(ic12) * rate(ircpg) + & ! C of CNO
               Qneu_rnpg * y(in14) * rate(irnpg) + & ! N of CNO
               Qneu_ropg * y(io16) * rate(iropg)) ! O of CNO
               
            d_epsneu_dy(ihe3) = Qconv*( &
               Qneu_rpp2 * y(ihe4) * rate(irhe3ag) * fII + & ! r34_pp2
               Qneu_rpp3 * y(ihe4) * rate(irhe3ag) * (1d0-fII)) ! r34_pp3
               
            d_epsneu_dy(ihe4) = Qconv*( &
               Qneu_rpp2 * y(ihe3) * rate(irhe3ag) * fII + & ! r34_pp2
               Qneu_rpp3 * y(ihe3) * rate(irhe3ag) * (1d0-fII)) ! r34_pp3
               
            d_epsneu_dy(ic12) = Qconv* &
               Qneu_rcpg * y(ih1) * rate(ircpg) ! C of CNO
               
            d_epsneu_dy(in14) = Qconv* &
               Qneu_rnpg * y(ih1) * rate(irnpg)  ! N of CNO
               
            d_epsneu_dy(io16) = Qconv* &
               Qneu_ropg * y(ih1) * rate(iropg) ! O of CNO

         end subroutine approx21_d_epsneu_dy


         subroutine approx21_dfdy( &
            y, dfdy, fe56ec_fake_factor_in, min_T, fe56ec_n_neut, &
            ratdum, dratdumdt, dratdumdd, dratdumdy1, dratdumdy2, btemp, plus_co56, ierr)
         real(dp), intent(in) :: fe56ec_fake_factor_in, min_T, btemp
         integer, intent(in) :: fe56ec_n_neut
         real(dp), intent(in), dimension(:) :: &
            y, ratdum, dratdumdt, dratdumdd, dratdumdy1, dratdumdy2
         real(dp), intent(inout) :: dfdy(:,:)
         logical, intent(in) ::  plus_co56
         integer, intent(out) :: ierr

         integer :: i,j
         real(dp) abar,zbar,ye,taud,taut, b1, &
               snuda,snudz,enuc,velx,posx,zz
         real(dp) :: fe56ec_fake_factor
               
         ierr = 0
         
         ! Turn on special fe56ec rate above some temperature
         fe56ec_fake_factor=eval_fe56ec_fake_factor(fe56ec_fake_factor_in,min_T,btemp)
            
         ! NOTE: use of quad precision for dfdy doesn't make a difference.

         dfdy(1:species(plus_co56),1:species(plus_co56)) = 0.0d0

   ! h1 jacobian elements
         dfdy(ih1,ih1)  = -3.0d0 * y(ih1) * ratdum(irpp) &
                        - 2.0d0 * y(ic12) * ratdum(ircpg) &
                        - 2.0d0 * y(in14) * ratdum(irnpg) &
                        - 2.0d0 * y(in14) * y(ih1) * dratdumdy1(irnpg) &
                        - 2.0d0 * y(io16) * ratdum(iropg) &
                        - 2.0d0 * y(io16) * y(ih1) * dratdumdy1(iropg) &
                        - 3.0d0 * ratdum(irpen)

         dfdy(ih1,ihe3) = 2.0d0 * y(ihe3) * ratdum(ir33) &
                        - y(ihe4) * ratdum(irhe3ag)

         dfdy(ih1,ihe4) = -y(ihe3) * ratdum(irhe3ag) &
                        - y(ihe3) * y(ihe4) * dratdumdy1(irhe3ag)

         dfdy(ih1,ic12) = -2.0d0 * y(ih1) * ratdum(ircpg)

         dfdy(ih1,in14) = -2.0d0 * y(ih1) * ratdum(irnpg)

         dfdy(ih1,io16) = -2.0d0 * y(ih1) * ratdum(iropg)


   ! he3 jacobian elements
         dfdy(ihe3,ih1)  =  y(ih1) * ratdum(irpp) &
                        + ratdum(irpen)

         dfdy(ihe3,ihe3) = -2.0d0 * y(ihe3) * ratdum(ir33) &
                        - y(ihe4) * ratdum(irhe3ag)

         dfdy(ihe3,ihe4) = -y(ihe3) * ratdum(irhe3ag) &
                        - y(ihe3) * y(ihe4) * dratdumdy1(irhe3ag)


   ! he4 jacobian elements
         dfdy(ihe4,ih1)  = y(in14) * ratdum(ifa) * ratdum(irnpg) &
                        + y(in14) * y(ih1) * ratdum(ifa) * dratdumdy1(irnpg) &
                        + y(io16) * ratdum(iropg) &
                        + y(io16) * y(ih1) * dratdumdy1(iropg)

         dfdy(ihe4,ihe3)  = y(ihe3) * ratdum(ir33) &
                        + y(ihe4) * ratdum(irhe3ag)


         dfdy(ihe4,ihe4)  = -1.5d0 * y(ihe4) * y(ihe4) * ratdum(ir3a) &
                           - y(ic12) * ratdum(ircag) &
                           - y(io16) * ratdum(iroag) &
                           - y(ine20) * ratdum(irneag) &
                           - y(img24) * ratdum(irmgag) &
                           - y(isi28) * ratdum(irsiag) &
                           - y(is32) * ratdum(irsag) &
                           - y(iar36) * ratdum(irarag) &
                           - y(ica40) * ratdum(ircaag) &
                           - y(iti44) * ratdum(irtiag) &
                           - y(icr48) * ratdum(ircrag) &
                           - y(ife52) * ratdum(irfeag)

         dfdy(ihe4,ihe4)  = dfdy(ihe4,ihe4) &
                           - y(img24) * ratdum(irmgap) * (1.0d0-ratdum(irr1)) &
                           - y(isi28) * ratdum(irsiap) * (1.0d0-ratdum(irs1)) &
                           - y(is32) * ratdum(irsap) * (1.0d0-ratdum(irt1)) &
                           - y(iar36) * ratdum(irarap) * (1.0d0-ratdum(iru1)) &
                           - y(ica40) * ratdum(ircaap) * (1.0d0-ratdum(irv1)) &
                           - y(iti44) * ratdum(irtiap) * (1.0d0-ratdum(irw1)) &
                           - y(icr48) * ratdum(ircrap) * (1.0d0-ratdum(irx1))

         dfdy(ihe4,ihe4)  = dfdy(ihe4,ihe4) &
                           - y(ife52) * ratdum(ir6f54) &
                           - y(ife52) * y(iprot) * ratdum(ir7f54) &
                           - ratdum(iralf1) & 
                           - y(ife54) * ratdum(irfe56_aux4) 


         dfdy(ihe4,ihe4)  = dfdy(ihe4,ihe4) &
                        + y(ihe3) * ratdum(irhe3ag) &
                        + y(ihe3) * y(ihe4) * dratdumdy1(irhe3ag) &
                        - y(in14) * ratdum(irnag) * 1.5d0


         dfdy(ihe4,ic12)  = y(ic12) * ratdum(ir1212) &
                           + 0.5d0 * y(io16) * ratdum(ir1216) &
                           + 3.0d0 * ratdum(irg3a) &
                           - y(ihe4) * ratdum(ircag)


         dfdy(ihe4,in14)  = y(ih1) * ratdum(ifa) * ratdum(irnpg) &
                        - y(ihe4) * ratdum(irnag) * 1.5d0


         dfdy(ihe4,io16)  = 0.5d0 * y(ic12) * ratdum(ir1216) &
                           + 1.12d0 * 0.5d0*y(io16) * ratdum(ir1616) &
                           + 0.68d0 * ratdum(irs1) * 0.5d0*y(io16) * ratdum(ir1616) &
                           + ratdum(iroga) &
                           - y(ihe4) * ratdum(iroag) &
                           + y(ih1) * ratdum(iropg)

         dfdy(ihe4,ine20) =  ratdum(irnega) &
                        - y(ihe4) * ratdum(irneag)

         dfdy(ihe4,img24) =   ratdum(irmgga) &
                           - y(ihe4) * ratdum(irmgag) &
                           - y(ihe4) * ratdum(irmgap) * (1.0d0-ratdum(irr1))

         dfdy(ihe4,isi28) =   ratdum(irsiga) &
                           - y(ihe4) * ratdum(irsiag) &
                           - y(ihe4) * ratdum(irsiap) * (1.0d0-ratdum(irs1)) &
                           + ratdum(irr1) * ratdum(irsigp)

         dfdy(ihe4,is32)  =   ratdum(irsga) &
                           - y(ihe4) * ratdum(irsag) &
                           - y(ihe4) * ratdum(irsap) * (1.0d0-ratdum(irt1)) &
                           + ratdum(irs1) * ratdum(irsgp)

         dfdy(ihe4,iar36) =   ratdum(irarga) &
                           - y(ihe4) * ratdum(irarag) &
                           - y(ihe4) * ratdum(irarap) * (1.0d0-ratdum(iru1)) &
                           + ratdum(irt1) * ratdum(irargp)

         dfdy(ihe4,ica40) =   ratdum(ircaga) &
                           - y(ihe4) * ratdum(ircaag) &
                           - y(ihe4) * ratdum(ircaap) * (1.0d0-ratdum(irv1)) &
                           + ratdum(iru1) * ratdum(ircagp)

         dfdy(ihe4,iti44) =   ratdum(irtiga) &
                           - y(ihe4) * ratdum(irtiag) &
                           - y(ihe4) * ratdum(irtiap) * (1.0d0-ratdum(irw1)) &
                           + ratdum(irv1) * ratdum(irtigp)

         dfdy(ihe4,icr48) =   ratdum(ircrga) &
                           - y(ihe4) * ratdum(ircrag) &
                           - y(ihe4) * ratdum(ircrap) * (1.0d0-ratdum(irx1)) &
                           + ratdum(irw1) * ratdum(ircrgp)

         dfdy(ihe4,ife52) =   ratdum(irfega) &
                           - y(ihe4) * ratdum(irfeag) &
                           + ratdum(irx1) * ratdum(irfegp) &
                           - y(ihe4) * ratdum(ir6f54) &
                           - y(ihe4) * y(iprot) * ratdum(ir7f54)

         dfdy(ihe4,ife54) =   y(iprot) * y(iprot) * ratdum(ir5f54) &
                           - y(ihe4) * ratdum(irfe56_aux4) 

         dfdy(ihe4,ife56) =   y(iprot) * y(iprot) * ratdum(irfe56_aux3)

         dfdy(ihe4,ini56) =   ratdum(irniga) &
                           + y(iprot) * ratdum(ir8f54)


         dfdy(ihe4,ineut) = -y(ihe4) * dratdumdy1(iralf1) &
                        + 2.0d0 * y(ineut) * y(iprot)*y(iprot) * ratdum(iralf2) &
                        + y(ineut)*y(ineut) * y(iprot)*y(iprot) * dratdumdy1(iralf2)
                        
         include 'formats'

         dfdy(ihe4,iprot) =   2.0d0 * y(ife54) * y(iprot) * ratdum(ir5f54) &
                           + y(ife54) * y(iprot) * y(iprot) * dratdumdy1(ir5f54) &
                           - y(ihe4) * y(ife52) * dratdumdy1(ir6f54) &
                           - y(ife52) * y(ihe4) * ratdum(ir7f54) &
                           - y(ife52) * y(ihe4) * y(iprot) * dratdumdy1(ir7f54) &
                           + y(ini56) * ratdum(ir8f54) &
                           + y(ini56) * y(iprot) * dratdumdy1(ir8f54) &
                           - y(ihe4) * dratdumdy2(iralf1) &
                           + 2.0d0 * y(ineut)*y(ineut) * y(iprot) * ratdum(iralf2) &
                           + y(ineut)*y(ineut) * y(iprot)*y(iprot) * dratdumdy2(iralf2) &
                           + 2.0d0 * y(ife56) * y(iprot) * ratdum(irfe56_aux3) &
                           + y(ife56) * y(iprot) * y(iprot) * dratdumdy1(irfe56_aux3) &
                           - y(ihe4) * y(ife54) * dratdumdy1(irfe56_aux4) 



   ! c12 jacobian elements
         dfdy(ic12,ih1)  = -y(ic12) * ratdum(ircpg) &
                        + y(in14) * ratdum(ifa) * ratdum(irnpg) &
                        + y(in14) * y(ih1) * ratdum(ifa) * dratdumdy1(irnpg)

         dfdy(ic12,ihe4) = 0.5d0 * y(ihe4) * y(ihe4) * ratdum(ir3a) &
                        - y(ic12) * ratdum(ircag)

         dfdy(ic12,ic12) = -2.0d0 * y(ic12) * ratdum(ir1212) &
                        - y(io16) * ratdum(ir1216) &
                        - ratdum(irg3a) &
                        - y(ihe4) * ratdum(ircag) &
                        - y(ih1) * ratdum(ircpg)

         dfdy(ic12,in14) = y(ih1) * ratdum(ifa) * ratdum(irnpg)

         dfdy(ic12,io16) = -y(ic12) * ratdum(ir1216) &
                        + ratdum(iroga)


   ! n14 jacobian elements
         dfdy(in14,ih1)  = y(ic12) * ratdum(ircpg) &
                        - y(in14) * ratdum(irnpg) &
                        - y(in14) * y(ih1) * dratdumdy1(irnpg) &
                        + y(io16) * ratdum(iropg) &
                        + y(io16) * y(ih1) * dratdumdy1(iropg)

         dfdy(in14,ihe4) = -y(in14) * ratdum(irnag)

         dfdy(in14,ic12) = y(ih1) * ratdum(ircpg)

         dfdy(in14,in14) = -y(ih1) * ratdum(irnpg) &
                        - y(ihe4) * ratdum(irnag)

         dfdy(in14,io16) = y(ih1) * ratdum(iropg)



   ! o16 jacobian elements
         dfdy(io16,ih1) = y(in14) * ratdum(ifg) * ratdum(irnpg) &
                     + y(in14) * y(ih1) * ratdum(ifg) * dratdumdy1(irnpg) &
                     - y(io16) * ratdum(iropg) &
                     - y(io16) * y(ih1) * dratdumdy1(iropg)

         dfdy(io16,ihe4) = y(ic12)*ratdum(ircag) &
                        - y(io16)*ratdum(iroag)

         dfdy(io16,ic12) = -y(io16)*ratdum(ir1216) &
                        + y(ihe4)*ratdum(ircag)

         dfdy(io16,in14) = y(ih1) * ratdum(ifg) * ratdum(irnpg)

         dfdy(io16,io16) = - y(ic12) * ratdum(ir1216) &
                        - 2.0d0 * y(io16) * ratdum(ir1616) &
                        - y(ihe4) * ratdum(iroag) &
                        - ratdum(iroga) &
                        - y(ih1) * ratdum(iropg)

         dfdy(io16,ine20) = ratdum(irnega)


   ! ne20 jacobian elements
         dfdy(ine20,ihe4)  = y(io16) * ratdum(iroag) &
                        - y(ine20) * ratdum(irneag) &
                        + y(in14) * ratdum(irnag)

         dfdy(ine20,ic12)  = y(ic12) * ratdum(ir1212)

         dfdy(ine20,in14)  = y(ihe4) * ratdum(irnag)

         dfdy(ine20,io16)  = y(ihe4) * ratdum(iroag)

         dfdy(ine20,ine20) = -y(ihe4) * ratdum(irneag) &
                           - ratdum(irnega)

         dfdy(ine20,img24) = ratdum(irmgga)



   ! mg24 jacobian elements
         dfdy(img24,ihe4)  = y(ine20) * ratdum(irneag) &
                           -y(img24) * ratdum(irmgag) &
                           -y(img24) * ratdum(irmgap) * (1.0d0-ratdum(irr1))

         dfdy(img24,ic12)  = 0.5d0 * y(io16) * ratdum(ir1216)

         dfdy(img24,io16)  = 0.5d0 * y(ic12) * ratdum(ir1216)

         dfdy(img24,ine20) = y(ihe4) * ratdum(irneag)

         dfdy(img24,img24) = -y(ihe4) * ratdum(irmgag) &
                           - ratdum(irmgga) &
                           - y(ihe4)*ratdum(irmgap)*(1.0d0-ratdum(irr1))

         dfdy(img24,isi28) = ratdum(irsiga) &
                           + ratdum(irr1) * ratdum(irsigp)



   ! si28 jacobian elements
         dfdy(isi28,ihe4)  = y(img24) * ratdum(irmgag) &
                        - y(isi28) * ratdum(irsiag) &
                        + y(img24) * ratdum(irmgap) * (1.0d0-ratdum(irr1)) &
                        - y(isi28) * ratdum(irsiap) * (1.0d0-ratdum(irs1))

         dfdy(isi28,ic12)  = 0.5d0 * y(io16) * ratdum(ir1216)

         dfdy(isi28,io16)  =   0.5d0 * y(ic12) * ratdum(ir1216) &
                           + 1.12d0 * 0.5d0*y(io16) * ratdum(ir1616) &
                           + 0.68d0 * 0.5d0*y(io16) * ratdum(irs1) * ratdum(ir1616)

         dfdy(isi28,img24) = y(ihe4) * ratdum(irmgag) &
                        + y(ihe4) * ratdum(irmgap) * (1.0d0-ratdum(irr1))

         dfdy(isi28,isi28) = -y(ihe4) * ratdum(irsiag) &
                           - ratdum(irsiga) &
                           - ratdum(irr1) * ratdum(irsigp) &
                           - y(ihe4) * ratdum(irsiap) * (1.0d0-ratdum(irs1))

         dfdy(isi28,is32)  = ratdum(irsga) &
                        + ratdum(irs1) * ratdum(irsgp)



   ! s32 jacobian elements
         dfdy(is32,ihe4)  = y(isi28) * ratdum(irsiag) &
                        - y(is32) * ratdum(irsag) &
                        + y(isi28) * ratdum(irsiap) * (1.0d0-ratdum(irs1)) &
                        - y(is32) * ratdum(irsap) * (1.0d0-ratdum(irt1))

         dfdy(is32,io16)  = &
                        + 0.68d0*0.5d0*y(io16)*ratdum(ir1616)*(1.0d0-ratdum(irs1)) &
                           + 0.2d0 * 0.5d0*y(io16) * ratdum(ir1616)

         dfdy(is32,isi28) = y(ihe4) * ratdum(irsiag) &
                        + y(ihe4) * ratdum(irsiap) * (1.0d0-ratdum(irs1))

         dfdy(is32,is32)  = -y(ihe4) * ratdum(irsag) &
                        - ratdum(irsga) &
                        - ratdum(irs1) * ratdum(irsgp) &
                        - y(ihe4) * ratdum(irsap) * (1.0d0-ratdum(irt1))

         dfdy(is32,iar36) = ratdum(irarga) &
                        + ratdum(irt1) * ratdum(irargp)



   ! ar36 jacobian elements
         dfdy(iar36,ihe4)  = y(is32) * ratdum(irsag) &
                        - y(iar36) * ratdum(irarag) &
                        + y(is32) * ratdum(irsap) * (1.0d0-ratdum(irt1)) &
                        - y(iar36) * ratdum(irarap) * (1.0d0-ratdum(iru1))

         dfdy(iar36,is32)  = y(ihe4) * ratdum(irsag) &
                           + y(ihe4) * ratdum(irsap) * (1.0d0-ratdum(irt1))

         dfdy(iar36,iar36) = -y(ihe4) * ratdum(irarag) &
                           - ratdum(irarga) &
                           - ratdum(irt1) * ratdum(irargp) &
                           - y(ihe4) * ratdum(irarap) * (1.0d0-ratdum(iru1))

         dfdy(iar36,ica40) = ratdum(ircaga) &
                        + ratdum(ircagp) * ratdum(iru1)



   ! ca40 jacobian elements
         dfdy(ica40,ihe4)   = y(iar36) * ratdum(irarag) &
                           - y(ica40) * ratdum(ircaag) &
                           + y(iar36) * ratdum(irarap)*(1.0d0-ratdum(iru1)) &
                           - y(ica40) * ratdum(ircaap)*(1.0d0-ratdum(irv1))

         dfdy(ica40,iar36)  = y(ihe4) * ratdum(irarag) &
                           + y(ihe4) * ratdum(irarap)*(1.0d0-ratdum(iru1))

         dfdy(ica40,ica40)  = -y(ihe4) * ratdum(ircaag) &
                           - ratdum(ircaga) &
                           - ratdum(ircagp) * ratdum(iru1) &
                           - y(ihe4) * ratdum(ircaap)*(1.0d0-ratdum(irv1))

         dfdy(ica40,iti44)  = ratdum(irtiga) &
                           + ratdum(irtigp) * ratdum(irv1)
            


   ! ti44 jacobian elements
         dfdy(iti44,ihe4)   = y(ica40) * ratdum(ircaag) &
                           - y(iti44) * ratdum(irtiag) &
                           + y(ica40) * ratdum(ircaap)*(1.0d0-ratdum(irv1)) &
                           - y(iti44) * ratdum(irtiap)*(1.0d0-ratdum(irw1))

         dfdy(iti44,ica40)  = y(ihe4) * ratdum(ircaag) &
                           + y(ihe4) * ratdum(ircaap)*(1.0d0-ratdum(irv1))

         dfdy(iti44,iti44)  = -y(ihe4) * ratdum(irtiag) &
                           - ratdum(irtiga) &
                           - ratdum(irv1) * ratdum(irtigp) &
                           - y(ihe4) * ratdum(irtiap)*(1.0d0-ratdum(irw1))

         dfdy(iti44,icr48)  = ratdum(ircrga) &
                           + ratdum(irw1) * ratdum(ircrgp)



   ! cr48 jacobian elements
         dfdy(icr48,ihe4)  = y(iti44) * ratdum(irtiag) &
                        - y(icr48) * ratdum(ircrag) &
                        + y(iti44) * ratdum(irtiap)*(1.0d0-ratdum(irw1)) &
                        - y(icr48) * ratdum(ircrap)*(1.0d0-ratdum(irx1))

         dfdy(icr48,iti44) = y(ihe4) * ratdum(irtiag) &
                        + y(ihe4) * ratdum(irtiap)*(1.0d0-ratdum(irw1))

         dfdy(icr48,icr48) = -y(ihe4) * ratdum(ircrag) &
                           - ratdum(ircrga) &
                           - ratdum(irw1) * ratdum(ircrgp) &
                           - y(ihe4) * ratdum(ircrap)*(1.0d0-ratdum(irx1))

         dfdy(icr48,ife52) = ratdum(irfega) &
                        + ratdum(irx1) * ratdum(irfegp)


   ! crx jacobian elements
         dfdy(icrx,ife56)  = fe56ec_fake_factor * ratdum(irn56ec)


   ! fe52 jacobian elements
         dfdy(ife52,ihe4)  = y(icr48) * ratdum(ircrag) &
                        - y(ife52) * ratdum(irfeag) &
                        + y(icr48) * ratdum(ircrap) * (1.0d0-ratdum(irx1)) &
                        - y(ife52) * ratdum(ir6f54) &
                        - y(ife52) * y(iprot) * ratdum(ir7f54)

         dfdy(ife52,icr48) = y(ihe4) * ratdum(ircrag) &
                        + y(ihe4) * ratdum(ircrap) * (1.0d0-ratdum(irx1))

         dfdy(ife52,ife52) = - y(ihe4) * ratdum(irfeag) &
                           - ratdum(irfega) &
                           - ratdum(irx1) * ratdum(irfegp) &
                           - y(ineut) * y(ineut) * ratdum(ir2f54) &
                           - y(ihe4) * ratdum(ir6f54) &
                           - y(ihe4) * y(iprot) * ratdum(ir7f54)

         dfdy(ife52,ife54) = ratdum(ir1f54) + &
                           y(iprot) * y(iprot) * ratdum(ir5f54)

         dfdy(ife52,ini56) = ratdum(irniga) &
                        + y(iprot) * ratdum(ir8f54)

         dfdy(ife52,ineut) = &
                           y(ife54) * dratdumdy1(ir1f54) &
                           - 2.0d0 * y(ife52) * y(ineut) * ratdum(ir2f54) &
                           - y(ife52) * y(ineut) * y(ineut) * dratdumdy1(ir2f54)

         dfdy(ife52,iprot) = 2.0d0 * y(ife54) * y(iprot) * ratdum(ir5f54) &
                        + y(ife54) * y(iprot) * y(iprot) * dratdumdy1(ir5f54) &
                        - y(ihe4) * y(ife52) * dratdumdy1(ir6f54) &
                        - y(ife52) * y(ihe4) * ratdum(ir7f54) &
                        - y(ife52) * y(ihe4) * y(iprot) * dratdumdy1(ir7f54) &
                        + y(ini56) * ratdum(ir8f54) &
                        + y(ini56) * y(iprot) * dratdumdy1(ir8f54)
                        

   ! fe54 jacobian elements
         dfdy(ife54,ihe4)  = y(ife52) * ratdum(ir6f54) & 
                           - y(ife54) * ratdum(irfe56_aux4)
                           
         dfdy(ife54,ife52) = &
                              y(ineut) * y(ineut) * ratdum(ir2f54) + &
                              y(ihe4) * ratdum(ir6f54)

         dfdy(ife54,ife54) = &
                           - ratdum(ir1f54) &
                           - y(ineut) * y(ineut) * ratdum(irfe56_aux2) & 
                           - y(iprot) * y(iprot) * ratdum(ir3f54) &
                           - y(iprot) * y(iprot) * ratdum(ir5f54) &
                           - y(ihe4) * ratdum(irfe56_aux4)

         dfdy(ife54,ife56) = &
                           ratdum(irfe56_aux1) + &
                           y(iprot) * y(iprot) * ratdum(irfe56_aux3)
   
         dfdy(ife54,ini56) = ratdum(ir4f54) 

         dfdy(ife54,ineut) = &
                           - y(ife54) * dratdumdy1(ir1f54) &
                           + 2.0d0 * y(ife52) * y(ineut) * ratdum(ir2f54) &
                           + y(ife52) * y(ineut) * y(ineut) * dratdumdy1(ir2f54) &
                           + y(ife56) * dratdumdy1(irfe56_aux1) & 
                           - 2.0d0 * y(ife54) * y(ineut) * ratdum(irfe56_aux2) &
                           - y(ife54) * y(ineut) * y(ineut) * dratdumdy1(irfe56_aux2) 

         dfdy(ife54,iprot) = -2.0d0 * y(ife54) * y(iprot) * ratdum(ir3f54) &
                           - y(ife54) * y(iprot) * y(iprot) * dratdumdy1(ir3f54) &
                           + y(ini56) * dratdumdy1(ir4f54) &
                           - 2.0d0 * y(ife54) * y(iprot) * ratdum(ir5f54) &
                           - y(ife54) * y(iprot) * y(iprot) * dratdumdy1(ir5f54) &
                           + y(ihe4) * y(ife52) * dratdumdy1(ir6f54) & 
                           + 2.0d0 * y(ife56) * y(iprot) * ratdum(irfe56_aux3) &
                           + y(ife56) * y(iprot) * y(iprot) * dratdumdy1(irfe56_aux3) &
                           - y(ihe4) * y(ife54) * dratdumdy1(irfe56_aux4)


   ! fe56 jacobian elements

         dfdy(ife56,ihe4)  = y(ife54) * ratdum(irfe56_aux4)


         dfdy(ife56,ife54) = &
                           y(ineut) * y(ineut) * ratdum(irfe56_aux2) + & 
                           y(ihe4) * ratdum(irfe56_aux4)

         dfdy(ife56,ife56)  = - fe56ec_fake_factor * ratdum(irn56ec) &
                              - ratdum(irfe56_aux1) & 
                              - y(iprot) * y(iprot) * ratdum(irfe56_aux3)

         if (plus_co56) then
            dfdy(ife56,ico56)  = ratdum(irco56ec)
         else
            dfdy(ife56,ini56)  = ratdum(irn56ec)
         end if


         dfdy(ife56,ineut) =  &
                           -y(ife56) * dratdumdy1(irfe56_aux1) &
                           + 2.0d0 * y(ife54) * y(ineut) * ratdum(irfe56_aux2) &
                           + y(ife54) * y(ineut) * y(ineut) * dratdumdy1(irfe56_aux2)
                           

         dfdy(ife56,iprot) = -2.0d0 * y(ife56) * y(iprot) * ratdum(irfe56_aux3) &
                           - y(ife56) * y(iprot) * y(iprot) * dratdumdy1(irfe56_aux3) &
                           + y(ihe4) * y(ife54) * dratdumdy1(irfe56_aux4)

         if (plus_co56) then
   ! co56 jacobian elements      
            dfdy(ico56,ini56) =  ratdum(irn56ec)
            dfdy(ico56,ico56) = -ratdum(irco56ec)
         end if


   ! ni56 jacobian elements
         dfdy(ini56,ihe4)  = y(ife52) * ratdum(irfeag) &
                        + y(ife52) * y(iprot) * ratdum(ir7f54)

         dfdy(ini56,ife52) = y(ihe4) * ratdum(irfeag) &
                        + y(ihe4)* y(iprot) * ratdum(ir7f54)

         dfdy(ini56,ife54) = y(iprot) * y(iprot) * ratdum(ir3f54)

         dfdy(ini56,ini56) = -ratdum(irniga) &
                           - ratdum(ir4f54) &
                           - y(iprot) * ratdum(ir8f54) &
                           - ratdum(irn56ec)

         dfdy(ini56,iprot) = 2.0d0 * y(ife54) * y(iprot) * ratdum(ir3f54) &
                        + y(ife54) * y(iprot) * y(iprot) * dratdumdy1(ir3f54) &
                        - y(ini56) * dratdumdy1(ir4f54) &
                        + y(ife52) * y(ihe4)* ratdum(ir7f54) &
                        + y(ife52) * y(ihe4)* y(iprot) * dratdumdy1(ir7f54) &
                        - y(ini56) * ratdum(ir8f54) &
                        - y(ini56) * y(iprot) * dratdumdy1(ir8f54)


   ! photodisintegration neutrons jacobian elements
         dfdy(ineut,ihe4)  = 2.0d0 * ratdum(iralf1)

         dfdy(ineut,ife52) = -2.0d0 * y(ineut) * y(ineut) * ratdum(ir2f54)                    
                              
         dfdy(ineut,ife54) =  2.0d0 * ratdum(ir1f54) &
                           - 2.0d0 * y(ineut) * y(ineut) * ratdum(irfe56_aux2)
                           
         dfdy(ineut,ife56) = 2.0d0 * ratdum(irfe56_aux1) &
                           - fe56ec_n_neut * fe56ec_fake_factor * ratdum(irn56ec)

         dfdy(ineut,ineut) = &
                           2.0d0 * y(ife54) * dratdumdy1(ir1f54) &
                           - 4.0d0 * y(ife52) * y(ineut) * ratdum(ir2f54) &
                           - 2.0d0 * y(ife52) * y(ineut) * y(ineut) * dratdumdy1(ir2f54) &
                           + 2.0d0 * y(ihe4) * dratdumdy1(iralf1) &
                           - 4.0d0 * y(ineut) * y(iprot)*y(iprot) * ratdum(iralf2) &
                           - 2.0d0 * y(ineut)*y(ineut) * y(iprot)*y(iprot) * dratdumdy1(iralf2) &
                           - ratdum(irnep) &
                           + 2.0d0 * y(ife56) * dratdumdy1(irfe56_aux1) & 
                           - 4.0d0 * y(ife54) * y(ineut) * ratdum(irfe56_aux2) &
                           - 2.0d0 * y(ife54) * y(ineut) * y(ineut) * dratdumdy1(irfe56_aux2)

         dfdy(ineut,iprot) = 2.0d0 * y(ihe4) * dratdumdy2(iralf1) &
                        - 4.0d0 * y(ineut)*y(ineut) * y(iprot) * ratdum(iralf2) &
                        - 2.0d0 * y(ineut)*y(ineut) * y(iprot)*y(iprot) * dratdumdy2(iralf2) &
                        + ratdum(irpen)

   ! photodisintegration protons jacobian elements
         dfdy(iprot,ihe4)  = 2.0d0 * y(ife52) * ratdum(ir6f54) &
                           + 2.0d0 * ratdum(iralf1) & 
                           + 2.0d0 * y(ife54) * ratdum(irfe56_aux4)

         dfdy(iprot,ife52) = 2.0d0 * y(ihe4) * ratdum(ir6f54)

         dfdy(iprot,ife54) = -2.0d0 * y(iprot) * y(iprot) * ratdum(ir3f54) &
                           - 2.0d0 * y(iprot) * y(iprot) * ratdum(ir5f54) & 
                           + 2.0d0 * y(ihe4) * ratdum(irfe56_aux4)

         dfdy(iprot,ife56) = -2.0d0 * y(iprot) * y(iprot) * ratdum(irfe56_aux3)

         dfdy(iprot,ini56) = 2.0d0 * ratdum(ir4f54)

         dfdy(iprot,ineut) = 2.0d0 * y(ihe4) * dratdumdy1(iralf1) &
                        - 4.0d0 * y(ineut) * y(iprot)*y(iprot) * ratdum(iralf2) &
                        - 2.0d0 * y(ineut)*y(ineut) * y(iprot)*y(iprot) * dratdumdy1(iralf2) &
                        + ratdum(irnep)

         dfdy(iprot,iprot) = -4.0d0 * y(ife54) * y(iprot) * ratdum(ir3f54) &
                           - 2.0d0 * y(ife54) * y(iprot)*y(iprot)*dratdumdy1(ir3f54) &
                           + 2.0d0 * y(ini56) * dratdumdy1(ir4f54) &
                           - 4.0d0 * y(ife54) * y(iprot) * ratdum(ir5f54) &
                           - 2.0d0 * y(ife54) * y(iprot)*y(iprot)*dratdumdy1(ir5f54) &
                           + 2.0d0 * y(ihe4) * y(ife52) * dratdumdy1(ir6f54) &
                           + 2.0d0 * y(ihe4) * dratdumdy2(iralf1) &
                           - 4.0d0 * y(ineut)*y(ineut) * y(iprot) * ratdum(iralf2) &
                           - 2.0d0 * y(ineut)*y(ineut) * y(iprot)*y(iprot) * dratdumdy2(iralf2) &
                           - ratdum(irpen) & 
                           - 4.0d0 * y(ife56) * y(iprot) * ratdum(irfe56_aux3) &
                           - 2.0d0 * y(ife56) * y(iprot) * y(iprot) * dratdumdy1(irfe56_aux3) &
                           + 2.0d0 * y(ihe4) * y(ife54) * dratdumdy1(irfe56_aux4)

         end subroutine approx21_dfdy
         
         
         subroutine approx21_dfdT_dfdRho( & ! epstotal includes neutrinos
               y, mion, dfdy, ratdum, dratdumdt, dratdumdd, &
               fe56ec_fake_factor, min_T, fe56ec_n_neut, temp, &
               dfdT, dfdRho, d_epstotal_dy,  plus_co56, ierr)
            real(dp), intent(in), dimension(:) :: &
               y, mion, ratdum, dratdumdt, dratdumdd
            real(dp), intent(in) :: fe56ec_fake_factor, min_T, temp, dfdy(:,:)
            integer, intent(in) :: fe56ec_n_neut
            real(dp), intent(inout), dimension(:) :: d_epstotal_dy, dfdT, dfdRho
            logical, intent(in) ::  plus_co56
            integer, intent(out) :: ierr
            
            integer :: i, j
            real(dp) :: enuc_conv2
            logical, parameter :: deriva = .true.
            
            ! temperature dependence of the rate equations            
            dfdT(1:species(plus_co56)) = 0d0
            call approx21_dydt( &
               y,dratdumdt,ratdum,dfdT,deriva,&
               fe56ec_fake_factor,min_T,fe56ec_n_neut,temp,plus_co56,ierr)
            if (ierr /= 0) return

            ! density dependence of the rate equations
            dfdRho(1:species(plus_co56)) = 0d0
            call approx21_dydt( &
               y,dratdumdd,ratdum,dfdRho,deriva,&
               fe56ec_fake_factor,min_T,fe56ec_n_neut,0d0,plus_co56,ierr)
            if (ierr /= 0) return

            ! energy generation rate partials (total energy; do neutrinos elsewhere)
            enuc_conv2 = -avo*clight*clight
            d_epstotal_dy(1:species(plus_co56)) = 0d0
            do j=1,species(plus_co56)
               do i=1,species(plus_co56)
                  d_epstotal_dy(j) = d_epstotal_dy(j) + dfdy(i,j)*mion(i)
               enddo
               d_epstotal_dy(j) = d_epstotal_dy(j) * enuc_conv2
            enddo
         
         end subroutine approx21_dfdT_dfdRho
      
      
         subroutine mark_approx21(handle, ierr)
            use net_def, only: Net_General_Info, get_net_ptr
            use chem_def, only: chem_isos
            integer, intent(in) :: handle
            integer, intent(out) :: ierr
            type (Net_General_Info), pointer :: g
            include 'formats'
            call get_net_ptr(handle, g, ierr)
            if (ierr /= 0) then
               write(*,*) 'invalid handle for do_mark_approx21_on_coprocessor'
               return
            end if
            call mark_approx21_isos( &
               g% net_iso, chem_isos% name(g% approx21_ye_iso),g% add_co56_to_approx21, ierr)
            if (ierr /= 0) return
            call mark_approx21_reactions(g% net_reaction,g% add_co56_to_approx21, ierr)
            if (ierr /= 0) return
         end subroutine mark_approx21


         subroutine set_approx21(handle, ierr)
            use net_def, only: Net_General_Info, get_net_ptr
            use chem_def, only: chem_isos
            integer, intent(in) :: handle
            integer, intent(out) :: ierr
            type (Net_General_Info), pointer :: g
            include 'formats'
            call get_net_ptr(handle, g, ierr)
            if (ierr /= 0) then
               write(*,*) 'invalid handle for do_mark_approx21_on_coprocessor'
               return
            end if
            call set_approx21_isos( &
               g% net_iso, chem_isos% name(g% approx21_ye_iso),g% add_co56_to_approx21,ierr)
            if (ierr /= 0) return
            call set_approx21_reactions(g% net_reaction,g% add_co56_to_approx21,ierr)
            if (ierr /= 0) return
         end subroutine set_approx21


         subroutine mark_approx21_isos(itab, ye_iso_name,plus_co56, ierr)
            use chem_lib, only: chem_get_iso_id
            integer :: itab(:)
            character (len=*), intent(in) :: ye_iso_name
            logical, intent(in) :: plus_co56
            integer, intent(out) :: ierr
            integer :: i, cid
            ierr = 0
            
            call do1('h1')
            call do1('he3')
            call do1('he4')
            call do1('c12')
            call do1('n14')
            call do1('o16')
            call do1('ne20')
            call do1('mg24')
            call do1('si28')
            call do1('s32')
            call do1('ar36')
            call do1('ca40')
            call do1('ti44')
            call do1('cr48')
            call do1('fe52')
            call do1('fe54')
            call do1('fe56')
            if (plus_co56) call do1('co56')
            call do1('ni56')
            call do1('neut')
            call do1('prot')
            call do1(ye_iso_name)
         
            contains
         
            subroutine do1(str)
               use utils_lib, only: mesa_error
               character (len=*), intent(in) :: str
               integer :: cid
               cid = chem_get_iso_id(str)
               if (cid <= 0) then
                  ierr = -1
                  write(*,*) 'mark_approx21_isos failed for ' // trim(str)
                  call mesa_error(__FILE__,__LINE__)
               end if
               itab(cid) = 1
            end subroutine do1

         end subroutine mark_approx21_isos      
         
         
         subroutine set_approx21_isos(itab, ye_iso_name, plus_co56, ierr)
            use chem_lib, only: chem_get_iso_id
            use const_def, only: ev2erg, clight
            integer :: itab(:)
            character (len=*), intent(in) :: ye_iso_name
            logical, intent(in) :: plus_co56
            integer, intent(out) :: ierr
            integer :: i, cid
            ierr = 0
            
            ih1   = do1('h1')
            ihe3  = do1('he3')
            ihe4  = do1('he4')
            ic12  = do1('c12')
            in14  = do1('n14')
            io16  = do1('o16')
            ine20 = do1('ne20')
            img24 = do1('mg24')
            isi28 = do1('si28')
            is32  = do1('s32')
            iar36 = do1('ar36')
            ica40 = do1('ca40')
            iti44 = do1('ti44')
            icr48 = do1('cr48')
            ife52 = do1('fe52')
            ife54 = do1('fe54')
            ife56 = do1('fe56')
            if (plus_co56) ico56 = do1('co56')
            ini56 = do1('ni56')
            ineut = do1('neut')
            iprot = do1('prot')
            icrx = do1(ye_iso_name)
            iso_cid(icrx) = -1 ! different for different approx21 nets
         
            contains
         
            integer function do1(str)
               use chem_def, only: chem_isos
               use utils_lib, only: mesa_error
               character (len=*), intent(in) :: str
               integer :: cid
               cid = chem_get_iso_id(str)
               if (cid <= 0) then
                  write(*,*) 'set_approx21_isos failed for ' // trim(str)
                  call mesa_error(__FILE__,__LINE__)
               end if
               do1 = itab(cid)
               iso_cid(do1) = cid
            end function do1
            
         end subroutine set_approx21_isos
         
         
         subroutine mark_approx21_reactions(rtab, plus_co56, ierr)
            use rates_lib, only: rates_reaction_id
            integer :: rtab(:)
            logical, intent(in) :: plus_co56
            integer, intent(out) :: ierr
            integer :: i, ir
            include 'formats'
            ierr = 0
            
            call do1('r_he4_he4_he4_to_c12')
            call do1('r_c12_to_he4_he4_he4')
            call do1('r_c12_ag_o16')
            call do1('r1212')
            call do1('r1216')
            call do1('r1616')
            call do1('r_o16_ga_c12')
            call do1('r_o16_ag_ne20')
            call do1('r_ne20_ga_o16')
            call do1('r_ne20_ag_mg24')
            call do1('r_mg24_ga_ne20')
         
            call do1('r_mg24_ag_si28')
            call do1('r_si28_ga_mg24')
            call do1('r_mg24_ap_al27')
            call do1('r_al27_pa_mg24')
            call do1('r_al27_pg_si28')
            call do1('r_si28_gp_al27')
            call do1('r_si28_ag_s32')
            call do1('r_s32_ga_si28')
            call do1('r_si28_ap_p31')
            call do1('r_p31_pa_si28')
            call do1('r_p31_pg_s32')
            call do1('r_s32_gp_p31')
            call do1('r_s32_ag_ar36')
            call do1('r_ar36_ga_s32')
            call do1('r_s32_ap_cl35')
            call do1('r_cl35_pa_s32')
            call do1('r_cl35_pg_ar36')
            call do1('r_ar36_gp_cl35')
            call do1('r_ar36_ag_ca40')
            call do1('r_ca40_ga_ar36')
            call do1('r_ar36_ap_k39')
            call do1('r_k39_pa_ar36')
            call do1('r_k39_pg_ca40')
            call do1('r_ca40_gp_k39')
            call do1('r_ca40_ag_ti44')
            call do1('r_ti44_ga_ca40')
            call do1('r_ca40_ap_sc43')
            call do1('r_sc43_pa_ca40')
            call do1('r_sc43_pg_ti44')
            call do1('r_ti44_gp_sc43')
            call do1('r_ti44_ag_cr48')
            call do1('r_cr48_ga_ti44')
            call do1('r_ti44_ap_v47')
            call do1('r_v47_pa_ti44')
            call do1('r_v47_pg_cr48')
            call do1('r_cr48_gp_v47')
            call do1('r_cr48_ag_fe52')
            call do1('r_fe52_ga_cr48')
            call do1('r_cr48_ap_mn51')
            call do1('r_mn51_pa_cr48')
            call do1('r_mn51_pg_fe52')
            call do1('r_fe52_gp_mn51')
            call do1('r_fe52_ag_ni56')
            call do1('r_ni56_ga_fe52')
            call do1('r_fe52_ap_co55')
            call do1('r_co55_pa_fe52')
            call do1('r_co55_pg_ni56')
            call do1('r_ni56_gp_co55')

            ! for fe54 photodisintegration
            call do1('r_fe52_ng_fe53')
            call do1('r_fe53_gn_fe52')
            call do1('r_fe53_ng_fe54')
            call do1('r_fe54_gn_fe53')
            call do1('r_fe54_pg_co55')
            call do1('r_co55_gp_fe54')

            ! for he4 photodisintegration
            call do1('r_he3_ng_he4')
            call do1('r_he4_gn_he3')
            call do1('r_h1_ng_h2')
            call do1('r_h2_gn_h1')
            call do1('r_h2_pg_he3')
            call do1('r_he3_gp_h2')

            ! for weak reactions
            call do1('rprot_to_neut')
            call do1('rneut_to_prot')

            if (plus_co56) then
               call do1('rni56ec_to_co56')
               call do1('rco56ec_to_fe56')
            else
               call do1('rni56ec_to_fe56')
            end if

            ! ppchain
            call do1('rpp_to_he3')
            call do1('r_he3_he3_to_h1_h1_he4')
            call do1('r_he3_ag_be7')
            call do1('r_be7_wk_li7')
            call do1('r_be7_pg_b8')

            ! cno cycles
            call do1('r_c12_pg_n13')
            call do1('r_n14_pg_o15')
            call do1('r_o16_pg_f17')      
            call do1('r_n15_pg_o16')
            call do1('r_n15_pa_c12')      
            call do1('r_n14_ag_f18')

            ! for reactions to fe56 
            call do1('r_fe54_ng_fe55')
            call do1('r_fe55_gn_fe54')
            call do1('r_fe55_ng_fe56')
            call do1('r_fe56_gn_fe55')
            call do1('r_fe54_ap_co57')
            call do1('r_co57_pa_fe54')
            call do1('r_fe56_pg_co57')
            call do1('r_co57_gp_fe56')
            
            contains
         
            subroutine do1(str)
               use utils_lib, only: mesa_error
               character (len=*), intent(in) :: str
               integer :: ir
               ir = rates_reaction_id(str)
               if (ir <= 0) then
                  ierr = -1
                  write(*,*) 'mark_approx21_reactions failed for ' // trim(str)
                  call mesa_error(__FILE__,__LINE__)
               end if
               rtab(ir) = 1
            end subroutine do1
            
         end subroutine mark_approx21_reactions


         subroutine set_approx21_reactions(rtab, plus_co56, ierr)
            use rates_lib, only: rates_reaction_id
            use utils_lib, only: mesa_error
            integer :: rtab(:)
            logical, intent(in) :: plus_co56
            integer, intent(out) :: ierr
            ierr = 0
            
            ir3a = do1('r_he4_he4_he4_to_c12')
            irg3a = do1('r_c12_to_he4_he4_he4')
            ircag = do1('r_c12_ag_o16')
            ir1212 = do1('r1212')
            ir1216 = do1('r1216')
            ir1616 = do1('r1616')
            iroga = do1('r_o16_ga_c12')
            iroag = do1('r_o16_ag_ne20')
            irnega = do1('r_ne20_ga_o16')
            irneag = do1('r_ne20_ag_mg24')
            irmgga = do1('r_mg24_ga_ne20')
         
            irmgag = do1('r_mg24_ag_si28')
            irsiga = do1('r_si28_ga_mg24')
            irmgap = do1('r_mg24_ap_al27')
            iralpa = do1('r_al27_pa_mg24')
            iralpg = do1('r_al27_pg_si28')
            irsigp = do1('r_si28_gp_al27')
            irsiag = do1('r_si28_ag_s32')
            irsga = do1('r_s32_ga_si28')
            irsiap = do1('r_si28_ap_p31')
            irppa = do1('r_p31_pa_si28')
            irppg = do1('r_p31_pg_s32')
            irsgp = do1('r_s32_gp_p31')
            irsag = do1('r_s32_ag_ar36')
            irarga = do1('r_ar36_ga_s32')
            irsap = do1('r_s32_ap_cl35')
            irclpa = do1('r_cl35_pa_s32')
            irclpg = do1('r_cl35_pg_ar36')
            irargp = do1('r_ar36_gp_cl35')
            irarag = do1('r_ar36_ag_ca40')
            ircaga = do1('r_ca40_ga_ar36')
            irarap = do1('r_ar36_ap_k39')
            irkpa = do1('r_k39_pa_ar36')
            irkpg = do1('r_k39_pg_ca40')
            ircagp = do1('r_ca40_gp_k39')
            ircaag = do1('r_ca40_ag_ti44')
            irtiga = do1('r_ti44_ga_ca40')
            ircaap = do1('r_ca40_ap_sc43')
            irscpa = do1('r_sc43_pa_ca40')
            irscpg = do1('r_sc43_pg_ti44')
            irtigp = do1('r_ti44_gp_sc43')
            irtiag = do1('r_ti44_ag_cr48')
            ircrga = do1('r_cr48_ga_ti44')
            irtiap = do1('r_ti44_ap_v47')
            irvpa = do1('r_v47_pa_ti44')
            irvpg = do1('r_v47_pg_cr48')
            ircrgp = do1('r_cr48_gp_v47')
            ircrag = do1('r_cr48_ag_fe52')
            irfega = do1('r_fe52_ga_cr48')
            ircrap = do1('r_cr48_ap_mn51')
            irmnpa = do1('r_mn51_pa_cr48')
            irmnpg = do1('r_mn51_pg_fe52')
            irfegp = do1('r_fe52_gp_mn51')
            irfeag = do1('r_fe52_ag_ni56')
            irniga = do1('r_ni56_ga_fe52')
            irfeap = do1('r_fe52_ap_co55')
            ircopa = do1('r_co55_pa_fe52')
            ircopg = do1('r_co55_pg_ni56')
            irnigp = do1('r_ni56_gp_co55')

            ! for fe54 photodisintegration
            ir52ng = do1('r_fe52_ng_fe53')
            ir53gn = do1('r_fe53_gn_fe52')
            ir53ng = do1('r_fe53_ng_fe54')
            ir54gn = do1('r_fe54_gn_fe53')
            irfepg = do1('r_fe54_pg_co55')
            ircogp = do1('r_co55_gp_fe54')

            ! for he4 photodisintegration
            irheng = do1('r_he3_ng_he4')
            irhegn = do1('r_he4_gn_he3')
            irhng = do1('r_h1_ng_h2')
            irdgn = do1('r_h2_gn_h1')
            irdpg = do1('r_h2_pg_he3')
            irhegp = do1('r_he3_gp_h2')

            ! for weak reactions
            irpen = do1('rprot_to_neut')
            irnep = do1('rneut_to_prot')

            if (plus_co56) then
               irn56ec = do1('rni56ec_to_co56')
               irco56ec = do1('rco56ec_to_fe56')
            else
               irn56ec = do1('rni56ec_to_fe56')
            end if

            ! ppchain
            irpp = do1('rpp_to_he3')
            ir33 = do1('r_he3_he3_to_h1_h1_he4')
            irhe3ag = do1('r_he3_ag_be7')
            ir_be7_wk_li7 = do1('r_be7_wk_li7')
            ir_be7_pg_b8 = do1('r_be7_pg_b8')

            ! cno cycles
            ircpg = do1('r_c12_pg_n13')
            irnpg = do1('r_n14_pg_o15')
            iropg = do1('r_o16_pg_f17')      
            irn15pg = do1('r_n15_pg_o16')
            irn15pa = do1('r_n15_pa_c12')      
            irnag = do1('r_n14_ag_f18')

            ! for reactions to fe56 
            ir54ng = do1('r_fe54_ng_fe55')
            ir55gn = do1('r_fe55_gn_fe54')
            ir55ng = do1('r_fe55_ng_fe56')
            ir56gn = do1('r_fe56_gn_fe55')
            irfe54ap = do1('r_fe54_ap_co57')
            irco57pa = do1('r_co57_pa_fe54')
            irfe56pg = do1('r_fe56_pg_co57')
            irco57gp = do1('r_co57_gp_fe56')

            ! the equilibrium links come after the mesa reactions
            ifa = num_mesa_reactions(plus_co56)+1
            ifg = ifa+1

            irr1 = ifg+1
            irs1 = irr1+1
            irt1 = irs1+1
            iru1 = irt1+1
            irv1 = iru1+1
            irw1 = irv1+1
            irx1 = irw1+1

            ir1f54 = irx1+1
            ir2f54 = ir1f54+1
            ir3f54 = ir2f54+1
            ir4f54 = ir3f54+1
            ir5f54 = ir4f54+1
            ir6f54 = ir5f54+1
            ir7f54 = ir6f54+1
            ir8f54 = ir7f54+1

            iralf1 = ir8f54+1
            iralf2 = iralf1+1

            irfe56_aux1 = iralf2+1
            irfe56_aux2 = irfe56_aux1+1
            irfe56_aux3 = irfe56_aux2+1
            irfe56_aux4 = irfe56_aux3+1
            
            if( (plus_co56 .and. irfe56_aux4 /= num_reactions(plus_co56)) .or. &
               (.not.plus_co56 .and. irfe56_aux4 /= num_reactions(plus_co56))) then
               write(*,*) 'set_approx21_reactions found bad num_reactions'
               write(*,*) plus_co56,irfe56_aux4,num_reactions(plus_co56)
               call mesa_error(__FILE__,__LINE__)
            end if
            
            call init_approx21(plus_co56)
            
            contains
         
            integer function do1(str)
               use utils_lib, only: mesa_error
               character (len=*), intent(in) :: str
               integer :: ir
               ir = rates_reaction_id(str)
               if (ir <= 0) then
                  write(*,*) 'set_approx21_reactions failed for ' // trim(str)
                  call mesa_error(__FILE__,__LINE__)
               end if
               do1 = rtab(ir)
               if (do1 <= 0) then
                  write(*,*) 'set_approx21_reactions failed to find rate for ' // trim(str)
                  call mesa_error(__FILE__,__LINE__)
               end if
               rate_id(do1) = ir
            end function do1
                  
         end subroutine set_approx21_reactions
         
         
         ! call this after have set rate numbers
         subroutine init_approx21(plus_co56)
            integer :: i
            logical, intent(in) :: plus_co56
            include 'formats'

            ! set the names of the reaction rates (use mesa standard names)

            ratnam(ir3a)   = 'r_he4_he4_he4_to_c12'
            ratnam(irg3a)  = 'r_c12_to_he4_he4_he4'
            ratnam(ircag)  = 'r_c12_ag_o16'
            ratnam(ir1212) = 'r1212'
            ratnam(ir1216) = 'r1216'
            ratnam(ir1616) = 'r1616'
            ratnam(iroga)  = 'r_o16_ga_c12'
            ratnam(iroag)  = 'r_o16_ag_ne20'
            ratnam(irnega) = 'r_ne20_ga_o16'
            ratnam(irneag) = 'r_ne20_ag_mg24'
            ratnam(irmgga) = 'r_mg24_ga_ne20'

            ratnam(irmgag) = 'r_mg24_ag_si28'
            ratnam(irsiga) = 'r_si28_ga_mg24'
            ratnam(irmgap) = 'r_mg24_ap_al27'
            ratnam(iralpa) = 'r_al27_pa_mg24'
            ratnam(iralpg) = 'r_al27_pg_si28'
            ratnam(irsigp) = 'r_si28_gp_al27'
            ratnam(irsiag) = 'r_si28_ag_s32'
            ratnam(irsga)  = 'r_s32_ga_si28'
            ratnam(irsiap) = 'r_si28_ap_p31'
            ratnam(irppa)  = 'r_p31_pa_si28'
            ratnam(irppg)  = 'r_p31_pg_s32'
            ratnam(irsgp)  = 'r_s32_gp_p31'
            ratnam(irsag)  = 'r_s32_ag_ar36'
            ratnam(irarga) = 'r_ar36_ga_s32'
            ratnam(irsap)  = 'r_s32_ap_cl35'
            ratnam(irclpa) = 'r_cl35_pa_s32'
            ratnam(irclpg) = 'r_cl35_pg_ar36'
            ratnam(irargp) = 'r_ar36_gp_cl35'
            ratnam(irarag) = 'r_ar36_ag_ca40'
            ratnam(ircaga) = 'r_ca40_ga_ar36'
            ratnam(irarap) = 'r_ar36_ap_k39'
            ratnam(irkpa)  = 'r_k39_pa_ar36'
            ratnam(irkpg)  = 'r_k39_pg_ca40'
            ratnam(ircagp) = 'r_ca40_gp_k39'
            ratnam(ircaag) = 'r_ca40_ag_ti44'
            ratnam(irtiga) = 'r_ti44_ga_ca40'
            ratnam(ircaap) = 'r_ca40_ap_sc43'
            ratnam(irscpa) = 'r_sc43_pa_ca40'
            ratnam(irscpg) = 'r_sc43_pg_ti44'
            ratnam(irtigp) = 'r_ti44_gp_sc43'
            ratnam(irtiag) = 'r_ti44_ag_cr48'
            ratnam(ircrga) = 'r_cr48_ga_ti44'
            ratnam(irtiap) = 'r_ti44_ap_v47'
            ratnam(irvpa)  = 'r_v47_pa_ti44'
            ratnam(irvpg)  = 'r_v47_pg_cr48'
            ratnam(ircrgp) = 'r_cr48_gp_v47'
            ratnam(ircrag) = 'r_cr48_ag_fe52'
            ratnam(irfega) = 'r_fe52_ga_cr48'
            ratnam(ircrap) = 'r_cr48_ap_mn51'
            ratnam(irmnpa) = 'r_mn51_pa_cr48'
            ratnam(irmnpg) = 'r_mn51_pg_fe52'
            ratnam(irfegp) = 'r_fe52_gp_mn51'
            ratnam(irfeag) = 'r_fe52_ag_ni56'
            ratnam(irniga) = 'r_ni56_ga_fe52'
            ratnam(irfeap) = 'r_fe52_ap_co55'
            ratnam(ircopa) = 'r_co55_pa_fe52'
            ratnam(ircopg) = 'r_co55_pg_ni56'
            ratnam(irnigp) = 'r_ni56_gp_co55'

            ! for fe54 photodisintegration
            ratnam(ir52ng) = 'r_fe52_ng_fe53'
            ratnam(ir53gn) = 'r_fe53_gn_fe52'
            ratnam(ir53ng) = 'r_fe53_ng_fe54'
            ratnam(ir54gn) = 'r_fe54_gn_fe53'
            ratnam(irfepg) = 'r_fe54_pg_co55'
            ratnam(ircogp) = 'r_co55_gp_fe54'

            ! for he4 photodisintegration
            ratnam(irheng)  = 'r_he3_ng_he4'
            ratnam(irhegn)  = 'r_he4_gn_he3'
            ratnam(irhng)   = 'r_h1_ng_h2'
            ratnam(irdgn)   = 'r_h2_gn_h1'
            ratnam(irdpg)   = 'r_h2_pg_he3'
            ratnam(irhegp)  = 'r_he3_gp_h2'

            ! for weak reactions
            ratnam(irpen)   = 'rprot_to_neut'
            ratnam(irnep)   = 'rneut_to_prot'

            if (plus_co56) then
               ratnam(irn56ec) = 'r_ni56_wk_co56'
               ratnam(irco56ec) = 'r_co56_wk_fe56'
            else
               ratnam(irn56ec) = 'rni56ec_to_fe56'
            end if

            ! ppchain
            ratnam(irpp)    = 'rpp_to_he3'
            ratnam(ir33)    = 'r_he3_he3_to_h1_h1_he4'
            ratnam(irhe3ag) = 'r_he3_ag_be7'
            ratnam(ir_be7_wk_li7) = 'r_be7_wk_li7'
            ratnam(ir_be7_pg_b8) = 'r_be7_pg_b8'

            ! cno cycles
            ratnam(ircpg)   = 'r_c12_pg_n13'
            ratnam(irnpg)   = 'r_n14_pg_o15'
            ratnam(iropg)   = 'r_o16_pg_f17'

            ratnam(irn15pg) = 'r_n15_pg_o16'
            ratnam(irn15pa) = 'r_n15_pa_c12'

            ratnam(irnag)   = 'r_n14_ag_f18'
            
            ratnam(ir54ng)   = 'r_fe54_ng_fe55'
            ratnam(ir55gn)   = 'r_fe55_gn_fe54'
            ratnam(ir55ng)   = 'r_fe55_ng_fe56'
            ratnam(ir56gn)   = 'r_fe56_gn_fe55'
            ratnam(irfe54ap) = 'r_fe54_ap_co57'
            ratnam(irco57pa) = 'r_co57_pa_fe54'
            ratnam(irfe56pg) = 'r_fe56_pg_co57'
            ratnam(irco57gp) = 'r_co57_gp_fe56'


            ! the combo links

            ratnam(ifa)     ='fa'  ! this is fraction of n15 that goes to c12 by pa
            ratnam(ifg)     ='fg'  ! this is fraction of n15 that goes to o16 by pg

            ratnam(irr1)   = 'r1'
            ratnam(irs1)   = 's1'
            ratnam(irt1)   = 't1'
            ratnam(iru1)   = 'u1'
            ratnam(irv1)   = 'v1'
            ratnam(irw1)   = 'w1'
            ratnam(irx1)   = 'x1'

            ratnam(ir1f54) = 'r1f54'
            ratnam(ir2f54) = 'r2f54'
            ratnam(ir3f54) = 'r3f54' ! rfe54prot_to_ni56
            ratnam(ir4f54) = 'r4f54' ! rni56gprot_to_fe54
            ratnam(ir5f54) = 'r5f54'
            ratnam(ir6f54) = 'r6f54'
            ratnam(ir7f54) = 'r7f54' ! rfe52aprot_to_ni56
            ratnam(ir8f54) = 'r8f54' ! rni56gprot_to_fe52

            ratnam(iralf1) = 'ralf1'
            ratnam(iralf2) = 'ralf2'

            ratnam(irfe56_aux1) = 'rfe56aux1'
            ratnam(irfe56_aux2) = 'rfe56aux2'
            ratnam(irfe56_aux3) = 'rfe56aux3'
            ratnam(irfe56_aux4) = 'rfe56aux4'
            
            return
            
            do i=1,num_mesa_reactions(plus_co56)
               write(*,2) trim(ratnam(i)), i
            end do
            write(*,*) ''
            do i=num_mesa_reactions(plus_co56)+1,num_reactions(plus_co56)
               write(*,2) 'extra ' // trim(ratnam(i)), i
            end do
            stop 'init_approx21'

         end subroutine init_approx21
         
         real(dp) function eval_fe56ec_fake_factor(fe56ec_fake_factor,min_T,temp)
            real(dp), intent(in) :: fe56ec_fake_factor,min_T,temp
         
            eval_fe56ec_fake_factor = 0.d0
            if(temp >= min_T)then
               eval_fe56ec_fake_factor = fe56ec_fake_factor
            end if
      
         end function eval_fe56ec_fake_factor
      

         pure integer function num_reactions(plus_co56)
            logical, intent(in) :: plus_co56

            if(plus_co56) then
               num_reactions = approx21_plus_co56_nrat
            else
               num_reactions = approx21_nrat
            end if

         end function num_reactions


         pure integer function num_mesa_reactions(plus_co56)
            logical, intent(in) :: plus_co56

            if(plus_co56) then
               num_mesa_reactions = approx21_num_mesa_reactions_co56
            else
               num_mesa_reactions = approx21_num_mesa_reactions_21
            end if

         end function num_mesa_reactions

         pure integer function species(plus_co56)
            logical, intent(in) :: plus_co56

            if(plus_co56) then
               species = species_co56
            else
               species = species_21
            end if

         end function species


      end module net_approx21
