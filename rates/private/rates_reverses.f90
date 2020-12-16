! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
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

      module rates_reverses
      
      implicit none
      
      
      contains
      
      subroutine set_reaction_reverses
         use rates_def
      
         reverse_reaction_id(ir_al27_pa_mg24) = ir_mg24_ap_al27
         reverse_reaction_id(ir_mg24_ap_al27) = ir_al27_pa_mg24

         reverse_reaction_id(ir_ar36_ag_ca40) = ir_ca40_ga_ar36
         reverse_reaction_id(ir_ca40_ga_ar36) = ir_ar36_ag_ca40

         reverse_reaction_id(ir_ar36_ga_s32) = ir_s32_ag_ar36
         reverse_reaction_id(ir_s32_ag_ar36) = ir_ar36_ga_s32

         reverse_reaction_id(ir_b8_gp_be7) = ir_be7_pg_b8
         reverse_reaction_id(ir_be7_pg_b8) = ir_b8_gp_be7
         
         reverse_reaction_id(ir_c12_ag_o16) = ir_o16_ga_c12
         reverse_reaction_id(ir_o16_ga_c12) = ir_c12_ag_o16

         reverse_reaction_id(ir_c12_ap_n15) = ir_n15_pa_c12
         reverse_reaction_id(ir_n15_pa_c12) = ir_c12_ap_n15

         reverse_reaction_id(ir_c12_pg_n13) = ir_n13_gp_c12
         reverse_reaction_id(ir_n13_gp_c12) = ir_c12_pg_n13

         reverse_reaction_id(ir_c12_to_he4_he4_he4) = ir_he4_he4_he4_to_c12
         reverse_reaction_id(ir_he4_he4_he4_to_c12) = ir_c12_to_he4_he4_he4

         reverse_reaction_id(ir_c13_an_o16) = ir_n14_gp_c13
         reverse_reaction_id(ir_n14_gp_c13) = ir_c13_an_o16
         
         reverse_reaction_id(ir_ca40_ag_ti44) = ir_ti44_ga_ca40
         reverse_reaction_id(ir_ti44_ga_ca40) = ir_ca40_ag_ti44

         reverse_reaction_id(ir_cr48_ag_fe52) = ir_fe52_ga_cr48
         reverse_reaction_id(ir_fe52_ga_cr48) = ir_cr48_ag_fe52

         reverse_reaction_id(ir_cr48_ga_ti44) = ir_ti44_ag_cr48
         reverse_reaction_id(ir_ti44_ag_cr48) = ir_cr48_ga_ti44

         reverse_reaction_id(ir_f17_ap_ne20) = ir_o16_pg_f17
         reverse_reaction_id(ir_o16_pg_f17) = ir_f17_ap_ne20
         
         reverse_reaction_id(ir_f17_gp_o16) = ir_o14_ap_f17
         reverse_reaction_id(ir_o14_ap_f17) = ir_f17_gp_o16

         reverse_reaction_id(ir_f18_gp_o17) = ir_o17_pg_f18
         reverse_reaction_id(ir_o17_pg_f18) = ir_f18_gp_o17

         reverse_reaction_id(ir_f17_pg_ne18) = ir_ne18_gp_f17
         reverse_reaction_id(ir_ne18_gp_f17) = ir_f17_pg_ne18

         reverse_reaction_id(ir_f18_pg_ne19) = ir_ne19_gp_f18
         reverse_reaction_id(ir_ne19_gp_f18) = ir_f18_pg_ne19

         reverse_reaction_id(ir_f18_pa_o15) = ir_o15_ap_f18
         reverse_reaction_id(ir_o15_ap_f18) = ir_f18_pa_o15         
         
         reverse_reaction_id(ir_f19_gp_o18) = ir_o18_pg_f19
         reverse_reaction_id(ir_o18_pg_f19) = ir_f19_gp_o18

         reverse_reaction_id(ir_f19_pa_o16) = ir_o16_ap_f19
         reverse_reaction_id(ir_o16_ap_f19) = ir_f19_pa_o16

         reverse_reaction_id(ir_f19_pg_ne20) = ir_ne20_gp_f19
         reverse_reaction_id(ir_ne20_gp_f19) = ir_f19_pg_ne20

         reverse_reaction_id(ir_fe52_ag_ni56) = ir_ni56_ga_fe52
         reverse_reaction_id(ir_ni56_ga_fe52) = ir_fe52_ag_ni56

         reverse_reaction_id(ir_mg24_ag_si28) = ir_si28_ga_mg24
         reverse_reaction_id(ir_si28_ga_mg24) = ir_mg24_ag_si28

         reverse_reaction_id(ir_mg24_ga_ne20) = ir_ne20_ag_mg24
         reverse_reaction_id(ir_ne20_ag_mg24) = ir_mg24_ga_ne20
         
         reverse_reaction_id(ir_n13_pg_o14) = ir_o14_gp_n13
         reverse_reaction_id(ir_o14_gp_n13) = ir_n13_pg_o14
         
         reverse_reaction_id(ir_n14_ap_o17) = ir_o17_pa_n14
         reverse_reaction_id(ir_o17_pa_n14) = ir_n14_ap_o17

         reverse_reaction_id(ir_n14_pg_o15) = ir_o15_gp_n14
         reverse_reaction_id(ir_o15_gp_n14) = ir_n14_pg_o15
         
         reverse_reaction_id(ir_n15_ap_o18) = ir_o18_pa_n15
         reverse_reaction_id(ir_o18_pa_n15) = ir_n15_ap_o18

         reverse_reaction_id(ir_n15_pg_o16) = ir_o16_gp_n15
         reverse_reaction_id(ir_o16_gp_n15) = ir_n15_pg_o16

         reverse_reaction_id(ir_na23_pa_ne20) = ir_ne20_ap_na23
         reverse_reaction_id(ir_ne20_ap_na23) = ir_na23_pa_ne20
         
         reverse_reaction_id(ir_ne19_ga_o15) = ir_o15_ag_ne19
         reverse_reaction_id(ir_o15_ag_ne19) = ir_ne19_ga_o15

         reverse_reaction_id(ir_ne20_ga_o16) = ir_o16_ag_ne20
         reverse_reaction_id(ir_o16_ag_ne20) = ir_ne20_ga_o16
         
         reverse_reaction_id(ir_s32_ga_si28) = ir_si28_ag_s32
         reverse_reaction_id(ir_si28_ag_s32) = ir_s32_ga_si28
         
         reverse_reaction_id(irc12ap_to_o16) = iro16gp_to_c12
         reverse_reaction_id(iro16gp_to_c12) = irc12ap_to_o16
         
         reverse_reaction_id(iro16ap_to_ne20) = irne20gp_to_o16
         reverse_reaction_id(irne20gp_to_o16) = iro16ap_to_ne20
         
         reverse_reaction_id(irne20ap_to_mg24) = irmg24gp_to_ne20
         reverse_reaction_id(irmg24gp_to_ne20) = irne20ap_to_mg24
         
         reverse_reaction_id(irmg24ap_to_si28) = irsi28gp_to_mg24
         reverse_reaction_id(irsi28gp_to_mg24) = irmg24ap_to_si28
         
         reverse_reaction_id(irsi28ap_to_s32) = irs32gp_to_si28
         reverse_reaction_id(irs32gp_to_si28) = irsi28ap_to_s32
         
         reverse_reaction_id(irs32ap_to_ar36) = irar36gp_to_s32
         reverse_reaction_id(irar36gp_to_s32) = irs32ap_to_ar36
         
         reverse_reaction_id(irar36ap_to_ca40) = irca40gp_to_ar36
         reverse_reaction_id(irca40gp_to_ar36) = irar36ap_to_ca40
         
         reverse_reaction_id(irca40ap_to_ti44) = irti44gp_to_ca40
         reverse_reaction_id(irti44gp_to_ca40) = irca40ap_to_ti44
         
         reverse_reaction_id(irti44ap_to_cr48) = ircr48gp_to_ti44
         reverse_reaction_id(ircr48gp_to_ti44) = irti44ap_to_cr48
         
         reverse_reaction_id(ircr48ap_to_fe52) = irfe52gp_to_cr48
         reverse_reaction_id(irfe52gp_to_cr48) = ircr48ap_to_fe52
         
         reverse_reaction_id(irfe54ng_to_fe56) = irfe56gn_to_fe54
         reverse_reaction_id(irfe56gn_to_fe54) = irfe54ng_to_fe56
         
         reverse_reaction_id(irfe52neut_to_fe54) = irfe54g_to_fe52
         reverse_reaction_id(irfe54g_to_fe52) = irfe52neut_to_fe54

      end subroutine set_reaction_reverses

         
      end module rates_reverses


