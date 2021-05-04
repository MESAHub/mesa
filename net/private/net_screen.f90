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

      module net_screen
      
      use const_def,only: dp, pi, ln10
      use math_lib
      use chem_def, only: chem_isos, ih1, num_chem_isos
      use net_def, only: Net_General_Info, Net_Info
      use rates_def
      
      implicit none
      

      contains


      subroutine make_screening_tables(n, ierr)
         type (Net_Info), pointer :: n
         integer, intent(out) :: ierr
         real(dp) :: y(num_chem_isos)
         real(dp), dimension(3,0), target :: screen_h1, screen_he4
         y = 0
         screen_h1 = 0
         screen_he4 = 0
         call screen_net( &
            n% g, num_chem_isos, y, 1d0, 1d0, 0d0, 0d0, .true., &
            n% rate_raw, n% rate_raw_dT, n% rate_raw_dRho, &
            n% rate_screened, n% rate_screened_dT, n% rate_screened_dRho, &
            n% screening_mode,  &
            screen_h1, screen_he4, 0d0, 0d0, 0d0, 1d0, ierr)
      end subroutine make_screening_tables
      

      subroutine screen_net( &
            g, num_isos, y, temp, den, logT, logRho, init,  &
            rate_raw, rate_raw_dT, rate_raw_dRho, &
            rate_screened, rate_screened_dT, rate_screened_dRho, &
            screening_mode, &
            screen_h1, screen_he4, zbar, abar, z2bar, ye, ierr)

         use rates_def, only: Screen_Info, reaction_name
         use rates_lib, only: screen_set_context
         
         type (Net_General_Info), pointer  :: g
         integer, intent(in) :: num_isos, screening_mode
         real(dp), intent(in) :: y(:), temp, den, logT, logRho, &
            zbar, abar, z2bar, ye
         real(dp), intent(inout), dimension(:) :: &
            rate_raw, rate_raw_dT, rate_raw_dRho, &
            rate_screened, rate_screened_dT, rate_screened_dRho
         real(dp), intent(inout), dimension(:,:) :: screen_h1, screen_he4
         logical, intent(in) :: init
         integer, intent(out) :: ierr

         type (Screen_Info) :: sc
         integer :: num_reactions, i, ir, j, op_err
         real(dp) :: Tfactor, dTfactordt
         logical :: all_okay
         
         include 'formats'
         
         ierr = 0

         if (.not. init) then
            call screen_set_context( &
               sc, temp, den, logT, logRho, zbar, abar, z2bar,  &
               screening_mode, num_isos, y, g% z158)
         end if

         num_reactions = g% num_reactions
         
         do i = 1, num_reactions
            ir = g% reaction_id(i)
            if (ir == 0) then
               write(*,*) 'g% reaction_id(i) == 0', i, num_reactions
               stop 'screen_net'
            end if
            if (reaction_screening_info(3,ir) > 0) then
               call eval_screen_triple(  &
                  init, i, &
                  reaction_screening_info(1,ir),  &
                  reaction_screening_info(2,ir),   &
                  reaction_screening_info(3,ir),   &
                  i, sc, ir, ierr)
               if (ierr /= 0) then
                  write(*,*) 'screen_net failed in eval_screen_triple ' // &
                     trim(reaction_name(ir))                  
                  return
               end if
            else if (reaction_screening_info(2,ir) > 0) then
               call eval_screen_pair(  &
                  init, i, &
                  reaction_screening_info(1,ir),  &
                  reaction_screening_info(2,ir),   &
                  i, sc, ir, ierr)
               if (ierr /= 0) then
                  write(*,*) 'screen_net failed in eval_screen_pair ' // &
                     trim(reaction_name(ir))                  
                  return
               end if
            else
               rate_screened(i) = rate_raw(i)
               rate_screened_dT(i) = rate_raw_dT(i)
               rate_screened_dRho(i) = rate_raw_dRho(i)
            end if
         end do
         if (ierr /= 0) return
         
         call set_combo_screen_rates(num_isos, y, sc, ierr)
         if (ierr /= 0) then
            write(*,*) 'screen_net failed in set_combo_screen_rates'
            return
         end if
               
         if (nrattab > 1 .and. (logT < g% logTcut_lim .or. logT <= g% logTcut_lo)) then
            ! strong rates cutoff smoothly for logT < logTcut_lim
            if (logT <= g% logTcut_lo) then
               do i = 1, num_reactions
                  if (g% weak_reaction_index(i) > 0) cycle
                  rate_screened(i) = 0
                  rate_screened_dT(i) = 0
                  rate_screened_dRho(i) = 0
               end do
            else
               Tfactor = (logT - g% logTcut_lo)/(g% logTcut_lim - g% logTcut_lo)
               Tfactor = 0.5d0*(1 - cospi(Tfactor*Tfactor))
               dTfactordt = 0.5d0 * pi * sinpi(Tfactor*Tfactor) * &
                             2.d0/((g% logTcut_lim - g% logTcut_lo) * temp * ln10)
               do i = 1, num_reactions
                  if (g% weak_reaction_index(i) > 0) cycle
                  rate_screened_dT(i) = Tfactor * rate_screened_dT(i) + dTfactordt * rate_screened(i)
                  rate_screened_dRho(i) = Tfactor * rate_screened_dRho(i)
                  rate_screened(i) = Tfactor * rate_screened(i)
               end do
            end if
         end if
         
         
         contains
      
         subroutine screening_pair( &
               init, ir, jscr, sc, cid1, a1, z1, cid2, a2, z2, scor, scordt, scordd, ierr)
            use rates_lib, only: screen_init_AZ_info, screen_pair
            use rates_def, only: Screen_Info
            use chem_def, only: ih1, ih2, ihe4, ico55, ico57
            logical, intent(in) :: init
            integer, intent(in) :: ir, jscr
            type (Screen_Info) :: sc
            integer, intent(in) :: cid1, cid2
            real(dp), intent(in) :: a1, z1, a2, z2
            real(dp), intent(out) :: scor, scordt, scordd
            integer, intent(out) :: ierr

            integer :: i1, i2
            include 'formats'
            ierr = 0
            if (init) then
               call screen_init_AZ_info( &
                  a1, z1, a2, z2, &
                  g% zs13(jscr), g% zhat(jscr), g% zhat2(jscr), g% lzav(jscr), &
                  g% aznut(jscr), g% zs13inv(jscr), &
                  ierr)
               if (ierr /= 0) write(*,*) 'screen_init_AZ_info failed in screening_pair ' // &
                     trim(reaction_name(ir))    
            else
               if (cid1 > 0 .and. cid2 > 0) then
                  i1 = g% net_iso(cid1) 
                  i2 = g% net_iso(cid2)
                  if (i1 == 0 .or. i2 == 0) then ! not in current net
                     if (g% doing_approx21 .and. &
                              .not. (cid1 == ico55 .or. cid1 == ico57 .or. &
                                     cid2 == ico55 .or. cid2 == ico57 .or. &
                                     cid2 == ih2 .or. cid2 == ih2)) then
                        ! this skips screening things like al27 + p for approx21
                        scor = 1d0
                        scordt = 0d0
                        scordd = 0d0
                        return
                     end if
                  else
                     if (cid1 == ih1 .and. screen_h1(1,i2) > 0) then
                        scor = screen_h1(1,i2)
                        scordt = screen_h1(2,i2)
                        scordd = screen_h1(3,i2)
                        return
                     else if (cid1 == ihe4 .and. screen_he4(1,i2) > 0) then
                        scor = screen_he4(1,i2)
                        scordt = screen_he4(2,i2)
                        scordd = screen_he4(3,i2)
                        return
                     else if (cid2 == ih1 .and. screen_h1(1,i1) > 0) then
                        scor = screen_h1(1,i1)
                        scordt = screen_h1(2,i1)
                        scordd = screen_h1(3,i1)
                        return
                     else if (cid2 == ihe4 .and. screen_he4(1,i1) > 0) then
                        scor = screen_he4(1,i1)
                        scordt = screen_he4(2,i1)
                        scordd = screen_he4(3,i1)
                        return
                     end if
                  end if
               end if 
               call screen_pair( &
                  sc, a1, z1, a2, z2, screening_mode, &
                  g% zs13(jscr), g% zhat(jscr), g% zhat2(jscr), g% lzav(jscr), &
                  g% aznut(jscr), g% zs13inv(jscr), g% logTcut_lo, &
                  scor, scordt, scordd, ierr) 
               if (ierr /= 0) write(*,*) 'screen_pair failed in screening_pair ' // &
                     trim(reaction_name(ir)) 
               
               if (cid1 > 0 .and. cid2 > 0) then
                  i1 = g% net_iso(cid1) 
                  i2 = g% net_iso(cid2) 
                  if (i1 /= 0 .and. i2 /= 0) then
                     if (cid1 == ih1) then
                        screen_h1(1,i2) = scor
                        screen_h1(2,i2) = scordt
                        screen_h1(3,i2) = scordd
                        return
                     else if (cid1 == ihe4) then
                        screen_he4(1,i2) = scor
                        screen_he4(2,i2) = scordt
                        screen_he4(3,i2) = scordd
                        return
                     else if (cid2 == ih1) then
                        screen_h1(1,i1) = scor
                        screen_h1(2,i1) = scordt
                        screen_h1(3,i1) = scordd
                        return
                     else if (cid2 == ihe4) then
                        screen_he4(1,i1) = scor
                        screen_he4(2,i1) = scordt
                        screen_he4(3,i1) = scordd
                        return
                     end if
                  end if
               end if 
            end if         
         end subroutine screening_pair
     
         subroutine set_rate_screening(i, sc1a, sc1adt, sc1add)
            integer, intent(in) :: i
            real(dp), intent(in) :: sc1a, sc1adt, sc1add
            include 'formats'
            if (i == 0) return         
            rate_screened(i) = rate_raw(i)*sc1a
            rate_screened_dT(i) = rate_raw_dT(i)*sc1a + rate_raw(i)*sc1adt
            rate_screened_dRho(i) = rate_raw_dRho(i)*sc1a + rate_raw(i)*sc1add
         end subroutine set_rate_screening      
      
         subroutine eval_screen_pair(init, jscr, i1, i2, i, sc, ir, ierr)
            use rates_def, only: Screen_Info
            logical, intent(in) :: init
            integer, intent(in) :: jscr
            type (Screen_Info) :: sc
            integer, intent(in) :: i1, i2 ! chem id's for the isotopes
            integer, intent(in) :: i ! rate number
            integer, intent(in) :: ir
            integer, intent(out) :: ierr
            real(dp) :: sc1a, sc1adt, sc1add, a1, z1, a2, z2
            include 'formats'
            ierr = 0
            a1 = chem_isos% Z_plus_N(i1)
            z1 = dble(chem_isos% Z(i1))
            a2 = chem_isos% Z_plus_N(i2)
            z2 = dble(chem_isos% Z(i2))
            !if (z1 == 0d0 .or. z2 == 0d0) return
               ! with this, get bad burn result, but okay for restart
               ! without it, reversed. get good burn, bad restart.
               ! bad restart max diff ~ 5e-15 in abundances of n17, n18, o14
               ! perhaps tiny difference in some screening factor?
            call screening_pair( &
               init, ir, jscr, sc, i1, a1, z1, i2, a2, z2, sc1a, sc1adt, sc1add, ierr)
            if (ierr /= 0) return
            if (init) return
            call set_rate_screening(i, sc1a, sc1adt, sc1add)         
         end subroutine eval_screen_pair
      
         subroutine eval_screen_triple(init, jscr, i1_in, i2_in, i3_in, i, sc, ir, ierr)
            use rates_def, only: Screen_Info
            logical, intent(in) :: init
            integer, intent(in) :: jscr
            type (Screen_Info) :: sc
            integer, intent(in) :: i1_in, i2_in, i3_in ! chem id's for the isotopes
            integer, intent(in) :: i ! rate number
            integer, intent(in) :: ir
            integer, intent(out) :: ierr
            integer :: i1, i2, i3, ii
            real(dp) :: sc1, sc1dt, sc1dd
            real(dp) :: sc2, sc2dt, sc2dd
            real(dp) :: scor, scordt, scordd
            real(dp) :: a1, z1, a2, z2, a3, z3
            include 'formats'
            ierr = 0
            i1 = i1_in; i2 = i2_in; i3 = i3_in
            a1 = chem_isos% Z_plus_N(i1)
            z1 = dble(chem_isos% Z(i1))
            a2 = chem_isos% Z_plus_N(i2)
            z2 = dble(chem_isos% Z(i2))
            a3 = chem_isos% Z_plus_N(i3)
            z3 = dble(chem_isos% Z(i3))
            if (z2 == 0) then
               if (z1 == 0) return ! n + n + A
               ! have A + n + B
               ! swap 1 and 2 so have n + A + B
               ii = i2; i2 = i1; i1 = ii
               a1 = chem_isos% Z_plus_N(i1)
               z1 = dble(chem_isos% Z(i1))
               a2 = chem_isos% Z_plus_N(i2)
               z2 = dble(chem_isos% Z(i2))
            end if
            if (z3 == 0) then ! have A + B + n
               ! swap 1 and 3 so have n + A + B
               ii = i1; i1 = i3; i3 = ii
               a1 = chem_isos% Z_plus_N(i1)
               z1 = dble(chem_isos% Z(i1))
               a3 = chem_isos% Z_plus_N(i3)
               z3 = dble(chem_isos% Z(i3))
            end if
            call screening_pair( &
               init, ir, jscr, sc, i2, a2, z2, i3, a3, z3, sc2, sc2dt, sc2dd, ierr)
            if (ierr /= 0) return
            if (z1 == 0) then
               if (init) return
               call set_rate_screening(i, sc2, sc2dt, sc2dd)
               return ! n + (A + B)
            end if
            i2 = 0 ! 0 for cid to disable caching
            a2 = a2 + a3
            z2 = z2 + z3
            call screening_pair( &
               init, ir, jscr, sc, i1, a1, z1, i2, a2, z2, sc1, sc1dt, sc1dd, ierr)
            if (init) return
            scor = sc1*sc2
            scordt = sc1*sc2dt + sc1dt*sc2
            scordd = sc1*sc2dd + sc1dd*sc2
            call set_rate_screening(i, scor, scordt, scordd)
         
            if (.false.) write(*,2) 'scr 3 ' // trim(reaction_Name(ir)) &
                     // ' ' // trim(chem_isos% name(i1)) &
                     // ' ' // trim(chem_isos% name(i2)) &
                     // ' ' // trim(chem_isos% name(i3)),  &
                  ir, scor
         
         end subroutine eval_screen_triple

         subroutine set_combo_screen_rates(num_isos, y, sc, ierr)
            use rates_def, only: Screen_Info
            integer, intent(in) :: num_isos
            real(dp), intent(in) :: y(:)
            type (Screen_Info) :: sc
            integer, intent(out) :: ierr

            integer, pointer :: rtab(:)
            real(dp) :: rateII, rateIII, rsum, fII, fIII

            include 'formats'
          
            rtab => g% net_reaction
            ierr = 0
         
            if (rtab(ir34_pp2) /= 0 .and. rtab(ir34_pp3) /= 0) then
               if (rate_screened(rtab(ir34_pp2)) /= &
                     rate_screened(rtab(ir34_pp3))) then
                  ierr = -1
                  return
               end if
               if (rtab(ir_be7_wk_li7) /= 0) then
                  rateII  = rate_screened(rtab(ir_be7_wk_li7))
               else if (rtab(irbe7ec_li7_aux) /= 0) then
                  rateII  = rate_screened(rtab(irbe7ec_li7_aux))
               else
                  write(*,*) 'need either r_be7_wk_li7 or rbe7ec_li7_aux'
                  stop 'set_combo_screen_rates'
               end if
               if (rtab(ir_be7_pg_b8) /= 0) then
                  rateIII = y(g% net_iso(ih1)) * rate_screened(rtab(ir_be7_pg_b8))
               else if (rtab(irbe7pg_b8_aux) /= 0) then
                  rateIII = y(g% net_iso(ih1)) * rate_screened(rtab(irbe7pg_b8_aux))
               else
                  write(*,*) 'need either r_be7_pg_b8 or rbe7pg_b8_aux'
                  stop 'set_combo_screen_rates'
               end if
               rsum = rateII + rateIII
               if (rsum < 1d-50) then
                  fII = 0.5d0
               else
                  fII = rateII / rsum
               end if
               fIII = 1d0 - fII
               
               rate_screened(rtab(ir34_pp2)) = fII*rate_screened(rtab(ir34_pp2))
               rate_screened_dT(rtab(ir34_pp2)) = fII*rate_screened_dT(rtab(ir34_pp2))
               rate_screened_dRho(rtab(ir34_pp2)) = fII*rate_screened_dRho(rtab(ir34_pp2))

               rate_screened(rtab(ir34_pp3)) = fIII*rate_screened(rtab(ir34_pp3))
               rate_screened_dT(rtab(ir34_pp3)) = fIII*rate_screened_dT(rtab(ir34_pp3))
               rate_screened_dRho(rtab(ir34_pp3)) = fIII*rate_screened_dRho(rtab(ir34_pp3))

            end if

            if (rtab(irn14_to_c12) /= 0)  &
               call rate_for_pg_pa_branches( &
                        rtab(irn14pg_aux), rtab(irn15pg_aux), rtab(irn15pa_aux),  &
                        0, rtab(irn14_to_c12))         

            if (rtab(irn14_to_o16) /= 0) &
               call rate_for_pg_pa_branches( &
                        rtab(irn14pg_aux), rtab(irn15pg_aux), rtab(irn15pa_aux),  &
                        rtab(irn14_to_o16), 0)               
      
            if (rtab(ir1616ppa) /= 0)  &
               call rate_for_pg_pa_branches( &
                        rtab(ir1616p_aux), rtab(irp31pg_aux), rtab(irp31pa_aux),  &
                        0, rtab(ir1616ppa))         
      
            if (rtab(ir1616ppg) /= 0)  &
               call rate_for_pg_pa_branches( &
                        rtab(ir1616p_aux), rtab(irp31pg_aux), rtab(irp31pa_aux),  &
                        rtab(ir1616ppg), 0)         

            call rate_for_alpha_ap( &
                        irc12ap_aux, irn15pg_aux, irn15pa_aux,  &
                        irc12ap_to_o16)    

            call rate_for_alpha_gp( &
                        iro16gp_aux, irn15pg_aux, irn15pa_aux,  &
                        iro16gp_to_c12)         

            call rate_for_alpha_ap( &
                        iro16ap_aux, irf19pg_aux, irf19pa_aux,  &
                        iro16ap_to_ne20)    

            call rate_for_alpha_gp( &
                        irne20gp_aux, irf19pg_aux, irf19pa_aux,  &
                        irne20gp_to_o16)         
              
            call rate_for_alpha_ap( &
                        irne20ap_aux, irna23pg_aux, irna23pa_aux,  &
                        irne20ap_to_mg24)         
                                       
            call rate_for_alpha_gp( &
                        irmg24gp_aux, irna23pg_aux, irna23pa_aux,  &
                        irmg24gp_to_ne20)               
      
            call rate_for_alpha_ap( &
                        irmg24ap_aux, iral27pg_aux, iral27pa_aux,  &
                        irmg24ap_to_si28)         
                                       
            call rate_for_alpha_gp( &
                        irsi28gp_aux, iral27pg_aux, iral27pa_aux,  &
                        irsi28gp_to_mg24)                      
                           
            call rate_for_alpha_ap( &
                        irsi28ap_aux, irp31pg_aux, irp31pa_aux,  &
                        irsi28ap_to_s32)         

            call rate_for_alpha_gp( &
                        irs32gp_aux, irp31pg_aux, irp31pa_aux,  &
                        irs32gp_to_si28)         
            
            call rate_for_alpha_ap( &
                        irs32ap_aux, ircl35pg_aux, ircl35pa_aux,  &
                        irs32ap_to_ar36)         
               
            call rate_for_alpha_gp( &
                        irar36gp_aux, ircl35pg_aux, ircl35pa_aux,  &
                        irar36gp_to_s32)         
                    
            call rate_for_alpha_ap( &
                        irar36ap_aux, irk39pg_aux, irk39pa_aux,  &
                        irar36ap_to_ca40)         

            call rate_for_alpha_gp( &
                        irca40gp_aux, irk39pg_aux, irk39pa_aux,  &
                        irca40gp_to_ar36)         

            call rate_for_alpha_ap( &
                        irca40ap_aux, irsc43pg_aux, irsc43pa_aux,  &
                        irca40ap_to_ti44)         

            call rate_for_alpha_gp( &
                        irti44gp_aux, irsc43pg_aux, irsc43pa_aux,  &
                        irti44gp_to_ca40)         

            call rate_for_alpha_ap( &
                        irti44ap_aux, irv47pg_aux, irv47pa_aux,  &
                        irti44ap_to_cr48)         
            
            call rate_for_alpha_gp( &
                        ircr48gp_aux, irv47pg_aux, irv47pa_aux,  &
                        ircr48gp_to_ti44)         

            call rate_for_alpha_ap( &
                        ircr48ap_aux, irmn51pg_aux, irmn51pa_aux,  &
                        ircr48ap_to_fe52)         
            
            call rate_for_alpha_gp( &
                        irfe52gp_aux, irmn51pg_aux, irmn51pa_aux,  &
                        irfe52gp_to_cr48)                          
         

         end subroutine set_combo_screen_rates

         subroutine rate_for_alpha_ap(ir_start, irpg, irpa, ir_with_pg)
            integer, intent(in) :: ir_start, irpg, irpa, ir_with_pg
            integer, pointer :: rtab(:)
            include 'formats'
            if (ir_start == 0) return
            rtab => g% net_reaction
            if (rtab(ir_with_pg) == 0) return
            call rate_for_pg_pa_branches( &
               rtab(ir_start), rtab(irpg), rtab(irpa), rtab(ir_with_pg), 0)
         end subroutine rate_for_alpha_ap

         subroutine rate_for_alpha_gp(ir_start, irpg, irpa, ir_with_pa)
            integer, intent(in) :: ir_start, irpg, irpa, ir_with_pa         
            integer, pointer :: rtab(:)
            if (ir_start == 0) return
            rtab => g% net_reaction
            if (rtab(ir_with_pa) == 0) return
            call rate_for_pg_pa_branches( &
                  rtab(ir_start), rtab(irpg), rtab(irpa), 0, rtab(ir_with_pa))
         end subroutine rate_for_alpha_gp
         
         subroutine rate_for_pg_pa_branches(ir_start, irpg, irpa, ir_with_pg, ir_with_pa)
            integer, intent(in) :: ir_start, irpg, irpa, ir_with_pg, ir_with_pa
            
            real(dp) :: pg_raw_rate, pa_raw_rate, pg_frac, pa_frac
            real(dp) :: d_pg_frac_dT, d_pg_frac_dRho, d_pa_frac_dT, d_pa_frac_dRho
            real(dp) :: r, drdT, drdd, x
         
            if (ir_start == 0) then
               write(*,*) 'ir_start', ir_start
               if (irpg /= 0) write(*,*) trim(reaction_Name(g% reaction_id(irpg))) // ' irpg'
               if (irpa /= 0) write(*,*) trim(reaction_Name(g% reaction_id(irpa))) // ' irpa'
               if (ir_with_pg /= 0) write(*,*) trim(reaction_Name(g% reaction_id(ir_with_pg))) // ' ir_with_pg'
               if (ir_with_pa /= 0) write(*,*) trim(reaction_Name(g% reaction_id(ir_with_pa))) // ' ir_with_pa'
               stop 'rate_for_pg_pa_branches'
            end if
         
            if (irpg == 0) then
               write(*,*) 'irpg', irpg
               if (ir_with_pg /= 0) write(*,*) trim(reaction_Name(g% reaction_id(ir_with_pg))) // ' ir_with_pg'
               if (ir_with_pa /= 0) write(*,*) trim(reaction_Name(g% reaction_id(ir_with_pa))) // ' ir_with_pa'
               stop 'rate_for_pg_pa_branches'
            end if
         
            if (irpa == 0) then
               write(*,*) 'irpg', irpg
               if (ir_with_pg /= 0) write(*,*) trim(reaction_Name(g% reaction_id(ir_with_pg))) // ' ir_with_pg'
               if (ir_with_pa /= 0) write(*,*) trim(reaction_Name(g% reaction_id(ir_with_pa))) // ' ir_with_pa'
               stop 'rate_for_pg_pa_branches'
            end if
         
            pg_raw_rate = rate_raw(irpg)
            pa_raw_rate = rate_raw(irpa)
         
            if (pg_raw_rate + pa_raw_rate < 1d-99) then ! avoid divide by 0
               pg_raw_rate = 1; pa_raw_rate = 1
            end if
         
            pg_frac = pg_raw_rate / (pg_raw_rate + pa_raw_rate)
            pa_frac = 1 - pg_frac
         
            x = pg_raw_rate + pa_raw_rate
            d_pg_frac_dT =  &
               (pa_raw_rate*rate_raw_dT(irpg) - pg_raw_rate*rate_raw_dT(irpa)) / (x*x)
            d_pa_frac_dT = -d_pg_frac_dT
         
            d_pg_frac_dRho =  &
               (pa_raw_rate*rate_raw_dRho(irpg) - pg_raw_rate*rate_raw_dRho(irpa)) / (x*x)
            d_pa_frac_dRho = -d_pg_frac_dRho
         
            r    = rate_screened(ir_start)
            drdT = rate_screened_dT(ir_start)
            drdd = rate_screened_dRho(ir_start)
         
            if (ir_with_pg /= 0) then
               rate_screened(ir_with_pg) = r*pg_frac
               rate_screened_dT(ir_with_pg) = r*d_pg_frac_dT + drdT*pg_frac
               rate_screened_dRho(ir_with_pg) = r*d_pg_frac_dRho + drdd*pg_frac
            end if
         
            if (ir_with_pa /= 0) then
               rate_screened(ir_with_pa)  = r*pa_frac
               rate_screened_dT(ir_with_pa) = r*d_pa_frac_dT + drdT*pa_frac
               rate_screened_dRho(ir_with_pa) = r*d_pa_frac_dRho + drdd*pa_frac
            end if
               
         end subroutine rate_for_pg_pa_branches
         

      end subroutine screen_net



      end module net_screen

