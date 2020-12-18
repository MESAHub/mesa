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

      module raw_rates
      
      use rates_def
      use const_def !, only: missing_value, dp
      
      implicit none

      abstract interface
         subroutine rate_fcn(tf, temp, fr, rr)
           use const_def, only: dp
           use ratelib, only: T_factors
           type (T_Factors), pointer :: tf
           real(dp), intent(in) :: temp
           real(dp), intent(out) :: fr, rr
         end subroutine rate_fcn
      end interface
      
      logical, parameter :: show_rates = .false.

      
      contains
      
      
      subroutine set_raw_rates(n, irs, which_rates, temp, tf, rates, ierr)
         use rates_def, only : T_Factors
         integer, intent(in) :: n
         integer, intent(in) :: irs(:) ! (n) maps 1..n to reaction id
         integer, intent(in) :: which_rates(:) ! (rates_reaction_id_max)
         real(dp), intent(in) :: temp
         type (T_Factors), pointer :: tf
         real(dp), intent(inout) :: rates(:)
         integer, intent(out) :: ierr
         integer :: i, ir, op_err
         include 'formats.dek'
         ierr = 0

!x$OMP PARALLEL DO PRIVATE(i,ir,op_err)
         do i=1,n
            ir = irs(i)
            if (ir <= 0) cycle
            op_err = 0
            call set_raw_rate(ir, which_rates(ir), temp, tf, rates(i), op_err)
            if (op_err /= 0) then
               ierr = op_err
               cycle
            end if
         end do
!x$OMP END PARALLEL DO
      end subroutine set_raw_rates
      

      subroutine set_raw_rate(ir, which_rate, temp, tf, raw_rate, ierr)
         use ratelib
         use reaclib_eval, only: do_reaclib_lookup
         integer, intent(in) :: ir, which_rate
         real(dp), intent(in) :: temp
         type (T_Factors), pointer :: tf
         real(dp), intent(out) :: raw_rate
         integer, intent(out) :: ierr
         integer :: rir, reaclib_id_ir

         real(dp) :: rr

         include 'formats.dek'
         
         ierr = 0

         ! See if the rate or its reverse is being loaded from a rate_table
         rir = reverse_reaction_id(ir)
         if(rir/=0) then
            if (raw_rates_records(ir)% use_rate_table) then
               ! We want a rate from its table
               call eval_table(ir, tf, temp, raw_rate, rr, ierr)
               return
            end if

            ! We want to compute a rate from the "other" table
            if (raw_rates_records(rir)% use_rate_table) then
               reaclib_id_ir = do_reaclib_lookup(reaction_name(ir), reaclib_rates% reaction_dict)
               ! if ir == a reverse rate (e.g r_o16_ga_c12) will have reaclib_id_ir == 0
               ! if ir == a forward rate (e.g r_c12_ag_o16) will have reaclib_id_ir /= 0 

               if(reaclib_id_ir == 0 ) then
                  ! We want a reverse rate from a forward table
                  call eval_table_reverse(ir, rir, tf, temp, raw_rate, rr, ierr)
                  return
               else
                  ! We want a forward rate from a reverse table
                  write(*,*)
                  write(*,*) "ERROR: Can not evalute ",trim(reaction_name(ir)), &
                             " from detailed balance with ",trim(reaction_name(rir))
                  write(*,*) "Provide either both rates or only provide ",trim(reaction_name(ir))
                  write(*,*)
                  call mesa_error(__FILE__,__LINE__)
                  return
               end if
            end if            
         else if (raw_rates_records(ir)% use_rate_table) then
            ! Only ir is set as a table rate and rate does not have a reverse
            call eval_table(ir, tf, temp, raw_rate, rr, ierr)
            return
         end if


         select case(ir)

            case(ir_he4_he4_he4_to_c12) ! triple alpha to c12 
               call do1_of_3( &
                  rate_tripalf_nacre, rate_tripalf_jina, rate_tripalf_fxt)
      
            case(ir_c12_to_he4_he4_he4) ! c12 to 3 alpha
               call do1_of_3_reverse( &
                  rate_tripalf_nacre, rate_tripalf_jina, rate_tripalf_fxt)

            case(ir_c12_ag_o16)
               call do1_of_4( &
                  rate_c12ag_nacre, rate_c12ag_jina, rate_c12ag_kunz, rate_c12ag_fxt)

            case(ir_o16_ga_c12) ! o16(g, a)c12
               call do1_of_4_reverse( &
                  rate_c12ag_nacre, rate_c12ag_jina, rate_c12ag_kunz, rate_c12ag_fxt)

            case(ir1212) ! c12(c12,n)mg23, c12(c12,p)na23, c12(c12,a)ne20
               call do1_of_2(rate_c12c12_fxt_multi, rate_c12c12_fxt_basic)
               ! NOTE: Gasques option for c12+c12 is implemented in net, not in rates.
      
            case(ir1216)
               call do1(rate_c12o16_fxt)

            case(ir1616) ! o16 + o16 -> si28 + he4
               call do1_of_2(rate_o16o16_fxt, rate_o16o16_jina)

            case(ir1216_to_mg24) ! ! c12 + o16 -> mg24 + he4
               call do1_of_2(rate_c12o16_to_mg24_fxt, rate_c12o16a_jina)

            case(ir1216_to_si28) ! ! c12 + o16 -> si28
               call do1_of_2(rate_c12o16_to_si28_fxt, rate_c12o16p_jina)

            case(ir1616a) ! o16(o16, a)si28 
               call do1_of_2(rate_o16o16a_fxt, rate_o16o16a_jina)

            case(ir1616g) ! o16(o16, g)s32 
               ! no jina rate
               call do1(rate_o16o16g_fxt)

            case(ir1616p_aux) ! o16(o16, p)p31
               call do1_of_2(rate_o16o16p_fxt, rate_o16o16p_jina)

            case(ir1616ppa) ! o16(o16, p)p31(p, a)si28 
               call do1(rate_o16o16p_jina)

            case(ir1616ppg) ! o16(o16, p)p31(p, g)s32 
               call do1(rate_o16o16p_jina)

            case(ir_he3_ag_be7) ! he3(he4, g)be7
               call do1_of_2(rate_he3he4_nacre, rate_he3he4_jina)

            case(ir34_pp2) ! he4(he3, g)be7(e-, nu)li7(p, a)he4
               call do1_of_2(rate_he3he4_nacre, rate_he3he4_jina)

            case(ir34_pp3) ! he4(he3, g)be7(p, g)b8(e+, nu)be8(, a)he4
               call do1_of_2(rate_he3he4_nacre, rate_he3he4_jina)

            case(ir_b8_gp_be7) ! be7(p, g)b8
                ! no jina rate
               call do1_of_2_reverse(rate_be7pg_nacre, rate_be7pg_nacre)

            case(ir_be7_wk_li7)      ! be7(e-, nu)li7
               call do1_of_2(rate_be7em_fxt, rate_be7em_jina)

            case(ir_c12_ap_n15) ! c12(a, p)n15
               call do1_reverse(rate_n15pa_jina)

            case(ir_c12_pg_n13)
               call do1_of_2(rate_c12pg_nacre, rate_c12pg_jina)

            case(ir_c13_pg_n14)
               call do1_of_2(rate_c13pg_nacre, rate_c13pg_jina)

            case(ir_f17_ap_ne20) ! f17(a, p)ne20
               call do1(rate_f17ap_jina)

            case(ir_f17_gp_o16) ! f17(g, p)o16
               call do1_of_2_reverse(rate_o16pg_nacre, rate_o16pg_jina)

            case(ir_f18_gp_o17) ! f18(g, p)o17
               call do1_reverse(rate_o17pg_jina)

            case(ir_f19_gp_o18) ! f19(g, p)o18
               call do1_reverse(rate_o18pg_jina)

            case(ir_f19_pa_o16)  ! f19(p, a)o16
               call do1_of_2(rate_f19pa_nacre, rate_f19pa_jina)

            case(ir_f19_pg_ne20)
               call do1(rate_f19pg_jina)

            case(ir_h2_be7_to_h1_he4_he4)      ! be7(d, p)2he4
               call do1(rate_be7dp_jina)

            case(ir_h2_h2_to_he4)        ! h2(h2, g)he4
               call do1(rate_ddg_jina)

            case(ir_h2_he3_to_h1_he4)  ! he3(d, p)he4
               call do1(rate_he3d_jina)

            case(ir_h2_pg_he3)
               call do1_of_2(rate_dpg_nacre, rate_dpg_jina)

            case(ir_he3_be7_to_h1_h1_he4_he4)      ! be7(he3, 2p)2he4
               call do1(rate_be7he3_jina)

            case(ir_he3_he3_to_h1_h1_he4) ! he3(he3, 2p)he4
               call do1_of_2(rate_he3he3_nacre, rate_he3he3_jina)

            case(ir_h1_h1_he4_to_he3_he3) ! he4(2p, he3)he3
               call do1_of_2_reverse(rate_he3he3_nacre, rate_he3he3_jina)

            case(ir_li7_pa_he4) ! li7(p, a)he4
               call do1_of_2(rate_li7pa_nacre, rate_li7pa_jina)

            case(ir_mg24_ga_ne20) ! mg24(g, a)ne20
               call do1_reverse(rate_ne20ag_jina)

            case(ir_n13_gp_c12) ! n13(g, p)c12
               call do1_of_2_reverse(rate_c12pg_nacre, rate_c12pg_jina)

            case(ir_n13_pg_o14)
               call do1_of_2(rate_n13pg_nacre, rate_n13pg_jina)

            case(ir_n14_ag_f18) ! n14(a, g)f18
               call do1(rate_n14ag_jina)

            case(ir_n14_ap_o17) ! n14(a, p)o17
               call do1_reverse(rate_o17pa_jina)

            case(ir_n14_gp_c13) ! n14(g, p)c13
               call do1_of_2_reverse(rate_c13pg_nacre, rate_c13pg_jina)
               
            case(ir_n14_pg_o15)
               call do1_of_2(rate_n14pg_nacre, rate_n14pg_jina)

            case(ir_n15_ap_o18) ! n15(a, p)o18
               call do1_of_2_reverse(rate_o18pa_nacre, rate_o18pa_jina)

            case(ir_n15_pa_c12) ! n15(p, a)c12
               call do1(rate_n15pa_jina)

            case(ir_n15_pg_o16)
               call do1(rate_n15pg_jina)

            case(ir_ne20_ag_mg24) ! ne20(a, g)mg24
               call do1(rate_ne20ag_jina)

            case(ir_ne20_ap_na23) ! ne20(a, p)na23
               call do1(rate_ne20ap_jina)

            case(ir_ne20_ga_o16) ! ne20(g, a)o16
               call do1_of_2_reverse(rate_o16ag_nacre, rate_o16ag_jina)

            case(ir_ne20_gp_f19) ! ne20(g, p)f19
               call do1_reverse(rate_f19pg_jina)

            case(ir_ne22_pg_na23)
               call do1(rate_ne22pg_jina)

            case(ir_o14_gp_n13) ! o14(g, p)n13
               call do1_of_2_reverse(rate_n13pg_nacre, rate_n13pg_jina)

            case(ir_o15_gp_n14)  ! o15(g, p)n14
               call do1_of_2_reverse(rate_n14pg_nacre, rate_n14pg_jina)

            case(ir_o16_ag_ne20) ! o16(a, g)ne20
               call do1_of_2(rate_o16ag_nacre, rate_o16ag_jina)

            case(ir_o16_ap_f19)  ! o16(a, p)f19
               call do1_of_2_reverse(rate_f19pa_nacre, rate_f19pa_jina)

            case(ir_o16_gp_n15)  ! o16(g, p)n15
               call do1_reverse(rate_n15pg_jina)

            case(ir_o16_pg_f17)
               call do1_of_2(rate_o16pg_nacre, rate_o16pg_jina)

            case(ir_o17_pa_n14) ! o17(p, a)n14
               call do1( rate_o17pa_jina)

            case(ir_o17_pg_f18)
               call do1(rate_o17pg_jina)

            case(ir_o18_ag_ne22) ! o18(a, g)ne22 
               call do1(rate_o18ag_jina)

            case(ir_o18_pa_n15) ! o18(p, a)n15 and n15(a, p)o18
               call do1_of_2(rate_o18pa_nacre, rate_o18pa_jina)

            case(ir_o18_pg_f19)
               call do1(rate_o18pg_jina)

            case(iral27pa_aux) ! al27(p, a)mg24 
               call do1_reverse(rate_mg24ap_jina)

            case(iral27pg_aux) ! al27(p, g)si28
               call do1(rate_al27pg_jina)

            case(irar36ap_aux) ! ar36(a, p)k39
               call do1(rate_ar36ap_jina)

            case(irar36ap_to_ca40)  
               call do1(rate_ar36ap_jina)

            case(irar36gp_aux) ! ar36(g, p)cl35
               call do1_reverse(rate_cl35pg_jina)

            case(irar36gp_to_s32) 
               call do1_reverse(rate_cl35pg_jina)

            case(irbe7ec_li7_aux)      ! be7(e-, nu)li7(p, a)he4 
               call do1(rate_be7em_fxt)

            case(irbe7pg_b8_aux) ! be7(p, g)b8(e+, nu)be8(, a)he4
               call do1_of_2(rate_be7pg_nacre, rate_be7pg_jina)

            case(irc12_to_c13) ! c12(p, g)n13(e+nu)c13
               call do1_of_2(rate_c12pg_nacre, rate_c12pg_jina)

            case(irc12_to_n14) ! c12(p, g)n13(e+nu)c13(p, g)n14
               call do1_of_2(rate_c12pg_nacre, rate_c12pg_jina)

            case(irc12ap_aux)  ! c12(a, p)n15
               call do1_reverse(rate_n15pa_jina)

            case(irc12ap_to_o16) ! c12(a, p)n15(p, g)o16
               call do1_reverse(rate_n15pa_jina)

            case(irca40ap_aux) ! ca40(a, p)sc43 
               call do1(rate_ca40ap_jina)

            case(irca40ap_to_ti44)
               call do1(rate_ca40ap_jina)

            case(irca40gp_aux) ! ca40(g, p)k39 
               call do1_reverse(rate_k39pg_jina)

            case(irca40gp_to_ar36)  
               call do1_reverse(rate_k39pg_jina)

            case(ircl35pa_aux) ! cl35(p, a)s32 
               call do1_reverse(rate_s32ap_jina)

            case(ircl35pg_aux) ! cl35(p, g)ar36
               call do1(rate_cl35pg_jina)

            case(irco55gprot_aux) ! co55(g, prot)fe54 
               call do1_reverse(rate_fe54pg_jina)

            case(irco55pg_aux) ! co55(p, g)ni56 
               call do1(rate_co55pg_jina)

            case(irco55protg_aux) ! co55(prot, g)ni56 
               call do1(rate_co55pg_jina)

            case(ircr48ap_aux) ! cr48(a, p)mn51 
               call do1(rate_cr48ap_jina)

            case(ircr48ap_to_fe52)
               call do1(rate_cr48ap_jina)

            case(ircr48gp_aux) ! cr48(g, p)v47 
               call do1_reverse(rate_v47pg_jina)

            case(ircr48gp_to_ti44)
               call do1_reverse(rate_v47pg_jina)

            case(irf19pg_aux) ! f19(p, g)ne20
               call do1(rate_f19pg_jina)

            case(irfe52ap_aux) ! fe52(a, p)co55 
               call do1(rate_fe52ap_jina)

            case(irfe52ap_to_ni56) ! fe52(a, p)co55(p, g)ni56
               call do1(rate_fe52ap_jina)

            case(irfe52aprot_aux) ! fe52(a, prot)co55 
               call do1(rate_fe52ap_jina)

            case(irfe52aprot_to_fe54) ! fe52(a, prot)co55(g, prot)fe54
               call do1(rate_fe52ap_jina)

            case(irfe52aprot_to_ni56) ! fe52(a, prot)co55(prot, g)ni56
               call do1(rate_fe52ap_jina)

            case(irfe52gp_aux) ! fe52(g, p)mn51 
               call do1_reverse(rate_mn51pg_jina)

            case(irfe52gp_to_cr48)
               call do1_reverse(rate_mn51pg_jina)

            case(irfe52neut_to_fe54) ! fe52(neut, g)fe53(neut, g)fe54
               raw_rate = -1 ! rate calculated by special routine.

            case(irfe52ng_aux) ! fe52(n, g)fe53
               call do1(rate_fe52ng_jina)

            case(irfe53gn_aux) ! fe53(g, n)fe52
               call do1_reverse(rate_fe52ng_jina)

            case(irfe53ng_aux) ! fe53(n, g)fe54
               call do1(rate_fe53ng_jina)

            case(irfe54a_to_ni56) ! fe54 + alpha -> ni56 + 2 neut
               call do1(rate_fe54a_jina)

            case(irfe54an_aux) ! fe54(a,n)ni57                        
               call do1(rate_fe54an_jina)

            case(irfe54an_to_ni56) ! fe54(a,n)ni57(g,n)ni56
               call do1(rate_fe54an_jina)            

            case(irfe54aprot_to_fe56) ! fe54(a, prot)co57(g, prot)fe56
               call do1(rate_fe54ap_jina)            

            case(irfe54g_to_fe52) ! fe54(g, neut)fe53(g, neut)fe52
               call do1_reverse(rate_fe53ng_jina)

            case(irfe54ng_aux) ! fe54(neut, g)fe55                        
               call do1(rate_fe54ng_jina)

            case(irfe54ng_to_fe56) ! fe54(neut, g)fe55(neut, g)fe56
               raw_rate = -1 ! rate calculated by special routine.

            case(irfe54prot_to_fe52) ! fe54(prot, g)co55(prot, a)fe52
               raw_rate = -1 ! rate calculated by special routine.

            case(irfe54prot_to_ni56) ! fe54(prot, g)co55(prot, g)ni56
               raw_rate = -1 ! rate calculated by special routine.                       

            case(irfe54protg_aux) ! fe54(prot, g)co55
               call do1(rate_fe54pg_jina)

            case(irfe55gn_aux) ! fe55(g, neut)fe54                            
               call do1_reverse(rate_fe54ng_jina)

            case(irfe55ng_aux) ! fe55(neut, g)fe56                            
               call do1(rate_fe55ng_jina)

            case(irfe56ec_fake_to_mn56)
               raw_rate = -1 ! rate calculated by special routine.

            case(irfe56ec_fake_to_mn57)
               raw_rate = -1 ! rate calculated by special routine.

            case(irfe56ec_fake_to_cr56)
               raw_rate = -1 ! rate calculated by special routine.

            case(irfe56ec_fake_to_cr57)
               raw_rate = -1 ! rate calculated by special routine.

            case(irfe56ec_fake_to_cr58)
               raw_rate = -1 ! rate calculated by special routine.

            case(irfe56ec_fake_to_cr59)
               raw_rate = -1 ! rate calculated by special routine.

            case(irfe56ec_fake_to_cr60)
               raw_rate = -1 ! rate calculated by special routine.

            case(irfe56ec_fake_to_cr61)
               raw_rate = -1 ! rate calculated by special routine.

            case(irfe56ec_fake_to_cr62)
               raw_rate = -1 ! rate calculated by special routine.

            case(irfe56ec_fake_to_cr63)
               raw_rate = -1 ! rate calculated by special routine.

            case(irfe56ec_fake_to_cr64)
               raw_rate = -1 ! rate calculated by special routine.

            case(irfe56ec_fake_to_cr65)
               raw_rate = -1 ! rate calculated by special routine.

            case(irfe56ec_fake_to_cr66)
               raw_rate = -1 ! rate calculated by special routine.

            case(irfe56ee_to_ni56)
               raw_rate = -1 ! rate calculated by special routine.

            case(irfe56gn_aux) ! fe56(g, neut)fe55                            
               call do1_reverse(rate_fe55ng_jina)

            case(irfe56gn_to_fe54) ! fe56(g, neut)fe55(g, neut)fe54
               raw_rate = -1 ! rate calculated by special routine.

            case(irfe56prot_to_fe54) ! fe56(prot, g)co57(prot, a)fe54
               raw_rate = -1 ! rate calculated by special routine.

            case(irh2_protg_aux) ! h2(prot, g)he3
               call do1(rate_dpg_fxt)

            case(irh2g_neut_aux) ! h2(g, neut)prot
               call do1_reverse(rate_png_fxt)

            case(irhe3_neutg_aux) ! he3(neut, g)he4
               call do1(rate_he3ng_fxt)

            case(irhe3gprot_aux) ! he3(g, prot)h2
               call do1_of_2_reverse(rate_dpg_fxt, rate_dpg_jina)

            case(irhe4_breakup) ! he4(g, neut)he3(g, prot)h2(g, neut)prot
               raw_rate = -1 ! rate calculated by special routine.

            case(irhe4_rebuild) ! prot(neut, g)h2(prot, g)he3(neut, g)he4
               raw_rate = -1 ! rate calculated by special routine.

            case(irhe4g_neut_aux) ! he4(g, neut)he3
               call do1_reverse(rate_he3ng_fxt) ! no jina rate

            case(irk39pa_aux) ! k39(p, a)ar36
               call do1_reverse(rate_ar36ap_jina)

            case(irk39pg_aux) ! k39(p, g)ca40 
               call do1(rate_k39pg_jina)

            case(irmg24ap_aux) ! mg24(a, p)al27 
               call do1(rate_mg24ap_jina)

            case(irmg24ap_to_si28)
               call do1(rate_mg24ap_jina)

            case(irmg24gp_aux) ! mg24(g, p)na23 
               call do1_reverse(rate_na23pg_jina)

            case(irmg24gp_to_ne20) ! mg24(g, p)na23(p, a)ne20 
               call do1_reverse(rate_na23pg_jina)

            case(irmn51pg_aux) ! mn51(p, g)fe52 
               call do1(rate_mn51pg_jina)

            case(irn14_to_c12) ! n14(p, g)o15(e+nu)n15(p, a)c12
               call do1_of_3(rate_n14pg_nacre, rate_n14pg_jina, rate_n14pg_fxt)

            case(irn14_to_n15) ! n14(p, g)o15(e+nu)n15
               call do1_of_3(rate_n14pg_nacre, rate_n14pg_jina, rate_n14pg_fxt)

            case(irn14_to_o16) ! n14(p, g)o15(e+nu)n15(p, g)o16
               call do1_of_3(rate_n14pg_nacre, rate_n14pg_jina, rate_n14pg_fxt)

            case(irn14ag_lite) ! n14 + 1.5 alpha => ne20
               call do1(rate_n14ag_jina)

            case(irn14ag_to_o18) ! n14(a, g)f18(e+nu)o18
               call do1(rate_n14ag_jina)

            case(irn14gc12) ! n14 => c12 + neut + prot
               call do1_of_2_reverse(rate_c13pg_nacre, rate_c13pg_jina)

            case(irn14pg_aux) ! n14(p, g)o15
               call do1_of_3(rate_n14pg_nacre, rate_n14pg_jina, rate_n14pg_fxt)

            case(irn15pa_aux) ! n15(p, a)c12
               call do1( rate_n15pa_jina)

            case(irn15pg_aux) ! n15(p, g)o16
               call do1(rate_n15pg_jina)

            case(irna23pa_aux) ! na23(p, a)ne20 
               call do1(rate_na23pa_jina)

            case(irna23pg_aux) ! na23(p, g)mg24 
               call do1(rate_na23pg_jina)

            case(irne18ag_to_mg24) ! ne18(a, g)mg22 -> mg24
               call do1(rate_ne18ag_jina)

            case(irne18ap_to_mg22) ! ne18(a, p)na21(p, g)mg22
               call do1(rate_ne18ap_jina)

            case(irne18ap_to_mg24) ! ne18(a, p)na21(p, g)mg22 -> mg24
               call do1(rate_ne18ap_jina)

            case(irne19pg_to_mg22) ! ne19(p, g)na20(p, g)mg21(e+nu)na21(p, g)mg22 
               call do1(rate_ne19pg_jina)

            case(irne19pg_to_mg24) ! ne19(p, g)na20(p, g)mg21(e+nu)na21(p, g)mg22 -> mg24
               call do1(rate_ne19pg_jina)

            case(irne20ap_aux) ! ne20(a, p)na23
               call do1(rate_ne20ap_jina)

            case(irne20ap_to_mg24) ! ne20(a, p)na23(p, g)mg24
               call do1(rate_ne20ap_jina)

            case(irne20gp_aux) ! ne20(g, p)f19
               call do1_reverse(rate_f19pg_jina)

            case(irne20gp_to_o16) ! ne20(g, p)f19(p, a)o16
               call do1_reverse(rate_f19pg_jina)

            case(irne20pg_to_mg22) ! ne20(p, g)na21(p, g)mg22 
               call do1(rate_ne20pg_nacre)

            case(irne20pg_to_mg24) ! ne20(p, g)na21(p, g)mg22 -> mg24
               call do1_of_2(rate_ne20pg_nacre, rate_ne20pg_jina)

            case(irneut_to_prot) ! neut(e+nu)prot
               raw_rate = -1 ! rate calculated by special routine.

            case(irni56ec_to_fe54) ! ni56 + 2 e- => 56/54*fe54
               raw_rate = -1 ! rate calculated by special routine.

            case(irni56ec_to_fe56) ! ni56 + 2 e- => fe56
               raw_rate = -1 ! rate calculated by special routine.

            case(irni56ec_to_co56)
               raw_rate = -1 ! rate calculated by special routine.

            case(irco56ec_to_fe56)
               raw_rate = -1 ! rate calculated by special routine.

            case(irni56gp_aux) ! ni56(g, p)co55 
               call do1_reverse(rate_co55pg_jina)

            case(irni56gp_to_fe52) ! ni56(g, p)co55(p, a)fe52 
               raw_rate = -1 ! rate calculated by special routine.

            case(irni56gprot_aux) ! ni56(g, prot)co55 
               call do1_reverse(rate_co55pg_jina)

            case(irni56gprot_to_fe52) ! ni56(g, prot)co55(prot, a)fe52 
               raw_rate = -1 ! rate calculated by special routine.

            case(irni56gprot_to_fe54) ! ni56(g, prot)co55(g, prot)fe54 
               raw_rate = -1 ! rate calculated by special routine.

            case(irni56ng_to_fe54) ! ni56(n,g)ni57(n,a)fe54
               raw_rate = -1 ! rate calculated by special routine.

            case(irni57na_aux) ! ni57(n,a)fe54                            
               call do1_reverse(rate_fe54an_jina)

            case(iro16_to_n14) ! o16(p, g)f17(e+nu)o17(p, a)n14
               call do1_of_2(rate_o16pg_nacre, rate_o16pg_jina)

            case(iro16_to_o17) ! o16(p, g)f17(e+nu)o17
               call do1_of_2(rate_o16pg_nacre, rate_o16pg_jina)

            case(iro16ap_aux)  ! o16(a, p)f19
               call do1_of_2_reverse(rate_f19pa_nacre, rate_f19pa_jina)

            case(iro16ap_to_ne20)  ! o16(a, p)f19(p, a)ne20
               call do1_of_2_reverse(rate_f19pa_nacre, rate_f19pa_jina)

            case(iro16gp_aux)  ! o16(g, p)n15 
               call do1_reverse(rate_n15pg_jina)

            case(iro16gp_to_c12)  ! o16(g, p)n15(p, a)c12
               call do1_reverse(rate_n15pg_jina)

            case(iro17_to_o18) ! o17(p, g)f18(e+nu)o18
               call do1(rate_o17pg_jina)

            case(irp31pa_aux) ! p31(p, a)si28 
               call do1_reverse(rate_si28ap_jina)

            case(irp31pg_aux) ! p31(p, g)s32 
               call do1(rate_p31pg_jina)

            case(irpep_to_he3)        ! p(e-p, nu)h2(p, g)he3
               call do1_of_2(rate_pep_fxt, rate_pep_jina)

            case(irpp_to_he3) ! p(p, e+nu)h2(p, g)he3
               call do1_of_2(rate_pp_nacre, rate_pp_jina)

            case(irprot_neutg_aux) ! prot(neut, g)h2
               call do1(rate_png_fxt)

            case(irprot_to_neut) ! prot(e-nu)neut
               raw_rate = -1 ! rate calculated by special routine.

            case(irs32ap_aux) ! s32(a, p)cl35 
               call do1(rate_s32ap_jina)

            case(irs32ap_to_ar36) 
               call do1(rate_s32ap_jina)

            case(irs32gp_aux) ! s32(g, p)p31 
               call do1_reverse(rate_p31pg_jina)

            case(irs32gp_to_si28)    
               call do1_reverse(rate_p31pg_jina)

            case(irsc43pa_aux) ! sc43(p, a)ca40 
               call do1_reverse(rate_ca40ap_jina)

            case(irsc43pg_aux) ! sc43(p, g)ti44 
               call do1(rate_sc43pg_jina)

            case(irsi28ap_aux) ! si28(a, p)p31 
               call do1(rate_si28ap_jina)

            case(irsi28ap_to_s32) 
               call do1(rate_si28ap_jina)

            case(irsi28gp_aux) ! si28(g, p)al27
               call do1_reverse(rate_al27pg_jina)

            case(irsi28gp_to_mg24) 
               call do1_reverse(rate_al27pg_jina)

            case(irti44ap_aux) ! ti44(a, p)v47 
               call do1(rate_ti44ap_jina)

            case(irti44ap_to_cr48) 
               call do1(rate_ti44ap_jina)

            case(irti44gp_aux) ! ti44(g, p)sc43 
               call do1_reverse(rate_sc43pg_jina)

            case(irti44gp_to_ca40)
               call do1_reverse(rate_sc43pg_jina)

            case(irv47pa_aux) ! v47(p, a)ti44 
               call do1_reverse(rate_ti44ap_jina)

            case(irv47pg_aux) ! v47(p, g)cr48 
               call do1(rate_v47pg_jina)

            case(ir_h1_h1_wk_h2) ! p(p, e+nu)h2
               call do1_of_2(rate_pp_nacre, rate_pp_jina)

            case(ir_h1_h1_ec_h2)        ! p(e-p, nu)h2
               call do1_of_2(rate_pep_fxt, rate_pep_jina)

            case(irn14ag_to_ne22) ! n14(a, g)f18(e+nu)o18(a, g)ne22
               call do1(rate_n14ag_jina)

            case(irf19pa_aux)  ! f19(p, a)o16
               call do1_of_2(rate_f19pa_nacre, rate_f19pa_jina)

            case(ir_be7_pg_b8)
               call do1_of_2(rate_be7pg_nacre, rate_be7pg_jina)

            case(ir_b8_wk_he4_he4) ! b8(p=>n)be8=>2 he4
               call do1_of_2(rate_b8ep, rate_b8_wk_he4_he4_jina)               

            case(irmn51pa_aux) ! mn51(p, a)cr48 
               call do1_reverse(rate_cr48ap_jina)

            case(irfe54gn_aux) ! fe54(g, n)fe53
               call do1_reverse(rate_fe53ng_jina)

            case(irco55pa_aux) ! co55(p, a)fe52 
               call do1_reverse(rate_fe52ap_jina)

            case(irco55prota_aux) ! co55(prot, a)fe52 
               call do1_reverse(rate_fe52ap_jina)

            case(ir_h1_he3_wk_he4)  ! he3(p, e+nu)he4
               call do1_of_2(rate_hep_fxt, rate_hep_jina)
                        
            case(ir_he3_ng_he4)
               call do1(rate_he3ng_fxt)
            
            case(ir_he4_gn_he3)
               call do1_reverse(rate_he3ng_fxt)
            
            case(ir_h1_ng_h2)
               call do1(rate_png_fxt)
            
            case(ir_h2_gn_h1)
               call do1_reverse(rate_png_fxt)
            
            case(ir_he3_gp_h2)
               call do1_reverse(rate_dpg_nacre)
            
            case(ir_c12_c12_to_h1_na23)
               call do1(rate_c12_c12_to_h1_na23_jina)
            
            case(ir_he4_ne20_to_c12_c12)
               call do1(rate_he4_ne20_to_c12_c12_jina)
            
            case(ir_c12_c12_to_he4_ne20)
               call do1_reverse(rate_he4_ne20_to_c12_c12_jina)
            
            case(ir_he4_mg24_to_c12_o16)
               call do1(rate_he4_mg24_to_c12_o16_jina)

            case default
               call do_default(ierr)
               
         end select
         
         
         contains


         subroutine do_default(ierr) 
            integer, intent(out) :: ierr ! set ierr to -1 if cannot find rate                      
            real(dp) :: lambda, dlambda_dlnT, rlambda, drlambda_dlnT             
            include 'formats.dek'           
            ierr = 0                     
            ! look for rate in reaclib
            call get_reaclib_rate_and_dlnT( &
               ir, temp, lambda, dlambda_dlnT, rlambda, drlambda_dlnT, ierr)
            raw_rate = lambda
         end subroutine do_default


         subroutine do1(rate_fcn1)
            procedure(rate_fcn) :: rate_fcn1
            call eval_raw_rate(ir, rate_fcn1, tf, temp, raw_rate, rr, ierr)
         end subroutine do1
         
         
         subroutine do1_reverse(rate_fcn1)
            procedure(rate_fcn) :: rate_fcn1
            call eval_raw_rate(ir, rate_fcn1, tf, temp, rr, raw_rate, ierr)
         end subroutine do1_reverse
         
         
         subroutine do1_of_2(rate_fcn1, rate_fcn2)
            procedure(rate_fcn) :: rate_fcn1, rate_fcn2
            call eval_which_raw_rate(  &
               ir, min(2,which_rate),  &
               rate_fcn1, rate_fcn2, rate_fcn_null, rate_fcn_null,  &
               tf, temp, raw_rate, rr, ierr)
         end subroutine do1_of_2
         
         
         subroutine do1_of_2_reverse(rate_fcn1, rate_fcn2)
            procedure(rate_fcn) :: rate_fcn1, rate_fcn2
            real(dp) :: r
            call do1_of_2(rate_fcn1, rate_fcn2)
            r = raw_rate; raw_rate = rr; rr = r
         end subroutine do1_of_2_reverse
         
         
         subroutine do1_of_3(rate_fcn1, rate_fcn2, rate_fcn3)
            procedure(rate_fcn) :: rate_fcn1, rate_fcn2, rate_fcn3
            call eval_which_raw_rate(  &
               ir, min(3,which_rate),  &
               rate_fcn1, rate_fcn2, rate_fcn3, rate_fcn_null,  &
               tf, temp, raw_rate, rr, ierr)
         end subroutine do1_of_3
         
         
         subroutine do1_of_3_reverse(rate_fcn1, rate_fcn2, rate_fcn3)
            procedure(rate_fcn) :: rate_fcn1, rate_fcn2, rate_fcn3
            real(dp) :: r
            call do1_of_3(rate_fcn1, rate_fcn2, rate_fcn3)
            r = raw_rate; raw_rate = rr; rr = r
         end subroutine do1_of_3_reverse
         
         
         subroutine do1_of_4(rate_fcn1, rate_fcn2, rate_fcn3, rate_fcn4)
            procedure(rate_fcn) :: rate_fcn1, rate_fcn2, rate_fcn3, rate_fcn4
            call eval_which_raw_rate(  &
               ir, min(4,which_rate),  &
               rate_fcn1, rate_fcn2, rate_fcn3, rate_fcn4,  &
               tf, temp, raw_rate, rr, ierr)
         end subroutine do1_of_4
         
         
         subroutine do1_of_4_reverse(rate_fcn1, rate_fcn2, rate_fcn3, rate_fcn4)
            procedure(rate_fcn) :: rate_fcn1, rate_fcn2, rate_fcn3, rate_fcn4
            real(dp) :: r
            call do1_of_4(rate_fcn1, rate_fcn2, rate_fcn3, rate_fcn4)
            r = raw_rate; raw_rate = rr; rr = r
         end subroutine do1_of_4_reverse


      end subroutine set_raw_rate


      subroutine rate_fcn_null(tf, temp, fr, rr)
         use ratelib, only: T_factors
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         fr = -1; rr = -1
      end subroutine rate_fcn_null


      subroutine eval_which_raw_rate(  &
            ir, which_rate,  &
            rate_fcn1, rate_fcn2, rate_fcn3, rate_fcn4,  &
            tf, temp, fr, rr, ierr)
         use interp_1d_lib, only: interp_values
         use ratelib
         use const_def, only: pi
         use math_lib
         integer, intent(in) :: ir, which_rate
         procedure(rate_fcn) :: rate_fcn1, rate_fcn2, rate_fcn3, rate_fcn4
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         integer, intent(out) :: ierr
         real(dp) :: fr1, rr1, alfa, beta
         include 'formats.dek'
         ierr = 0
         if (which_rate == 2 .and. temp < JR_T_full_on) then
            call eval_raw_rate(ir, rate_fcn1, tf, temp, fr1, rr1, ierr)
            if (ierr /= 0) return
            if (temp <= JR_T_full_off) then
               fr = fr1; rr = rr1; return
            end if
         end if
         select case(which_rate)
            case (1)
               call eval_raw_rate(ir, rate_fcn1, tf, temp, fr, rr, ierr)
            case (2)
               call eval_raw_rate(ir, rate_fcn2, tf, temp, fr, rr, ierr)
            case (3)
               call eval_raw_rate(ir, rate_fcn3, tf, temp, fr, rr, ierr)
            case (4)
               call eval_raw_rate(ir, rate_fcn4, tf, temp, fr, rr, ierr)
            case default         
               write(*,*) 'bad which rate ' // trim(reaction_Name(ir)), which_rate
               ierr = -1
               return
         end select
         if (which_rate == 2 .and. temp < JR_T_full_on) then ! blend 1 and 2
            ! alfa = fraction from rate_fcn2; beta = 1-alfa = fraction from rate_fcn1
            alfa = (temp - JR_T_full_off) / (JR_T_full_on - JR_T_full_off)
            alfa = 0.5d0*(1d0 - cospi(alfa))
            beta = 1d0 - alfa
            fr = beta*fr1 + alfa*fr
            rr = beta*rr1 + alfa*rr
         end if
      end subroutine eval_which_raw_rate


      subroutine eval_raw_rate(ir, rate_fcn1, tf, temp, fr, rr, ierr)
         use ratelib
         integer, intent(in) :: ir ! reaction id
         procedure(rate_fcn) :: rate_fcn1
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         integer, intent(out) :: ierr
         include 'formats.dek'
         ierr = 0
         call rate_fcn1(tf, temp, fr, rr)
         if (fr < 0 .or. rr < 0) then
            write(*,1) 'invalid which rate for ' // trim(reaction_Name(ir)), fr, rr
            ierr = -1
            return
         end if
         if (fr == missing_value .or. rr == missing_value) then
            write(*,1) 'missing value for ' // trim(reaction_Name(ir)), fr, rr
            ierr = -1
            return
         end if
      end subroutine eval_raw_rate
      
      subroutine eval_table(ir, tf, temp, fr, rr, ierr)
         use interp_1d_lib, only: interp_values
         use ratelib
         integer, intent(in) :: ir ! reaction id
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         integer, intent(out) :: ierr
         integer, parameter :: nv = 1
         real(dp) :: x(nv), vals(nv)
         type (rate_table_info), pointer :: ri
         include 'formats.dek'
         ierr = 0

         ri => raw_rates_records(ir)
         if (.not. ri% use_rate_table) then
            ierr = -1
            return
         end if    
         if (ri% need_to_read) then
!$omp critical (load_rate_table)
            if (ri% need_to_read) then
               call get_interp_table(ri% rate_fname, ri% nT8s, ri% T8s, ri% f1, ierr)
               ri% need_to_read = .false.
            end if
!$omp end critical (load_rate_table)
            if (ierr /= 0) then
               write(*,*) 'failed to load table ' // trim(ri% rate_fname)
               return
            end if
         end if        
         x(1) = temp*1d-8
         call interp_values(ri% T8s, ri% nT8s, ri% f1, nv, x, vals, ierr)
         fr = vals(1)
         rr = 0 ! no reverse rates for tables
      end subroutine eval_table


      subroutine eval_table_reverse(ir, rir, tf, temp, fr, rr, ierr)
         use interp_1d_lib, only: interp_values
         use ratelib
         use reaclib_eval, only: compute_some_inverse_lambdas,&
                                 do_reaclib_indices_for_reaction
         integer, intent(in) :: ir, rir ! reaction id, reverse reaction id
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         integer, intent(out) :: ierr
         real(dp), dimension(1):: ln_lambda, lambda, dlambda_dlnT, &
                                    inv_lambda, dinv_lambda_dlnT
         integer :: lo, hi, num_lambdas
         real(dp) :: fr_table
         include 'formats.dek'
         ierr = 0

         call eval_table(rir, tf, temp, fr_table, rr, ierr)

         lambda = 1d0
         ln_lambda = 0d0
         dlambda_dlnT = 0d0

         call do_reaclib_indices_for_reaction(reaction_name(rir), reaclib_rates, lo, hi, ierr)
         if(ierr/=0) return 

         inv_lambda = 0d0
         dinv_lambda_dlnT = 0d0

         call compute_some_inverse_lambdas(1, lo, lo, &
                                          tf%T9, reaclib_rates, &
                                          ln_lambda, lambda, dlambda_dlnT, &
                                          inv_lambda, dinv_lambda_dlnT)

         fr = inv_lambda(1) * fr_table 

         rr = 0 
      end subroutine eval_table_reverse


      subroutine get_interp_table(f_name, nT8s, T8s_out, f1_out, ierr)
         use interp_1d_def, only: pm_work_size
         use interp_1d_lib, only: interp_pm
         use utils_lib
         use math_lib, only: str_to_vector
         character (len=*), intent(in) :: f_name
         integer, intent(out) :: nT8s
         real(dp), pointer :: T8s_out(:) ! will be allocated.  (nT8s)
         real(dp), pointer :: f1_out(:) ! will be allocated.  (4,nT8s)
         integer, intent(out) :: ierr
         
         integer :: iounit, j, nvec
         real(dp) :: tmp
         real(dp), pointer :: work(:)
         real(dp), pointer :: T8s(:)
         real(dp), pointer :: f1(:), f(:,:)
         character (len=256) :: line, rate_file
         real(dp), target :: vec_ary(20)
         real(dp), pointer :: vec(:)
         
         ierr = 0
         vec => vec_ary

         rate_file = trim(rates_table_dir) // '/' // trim(f_name)
         !write(*,*) 'load table ' // trim(rate_file)
         
         open(newunit=iounit,file=trim(rate_file),action='read',status='old',iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'ERROR: cannot open rate info file ' // trim(rate_file)
            !return
            call mesa_error(__FILE__,__LINE__)
         end if

         do ! read until reach line starting with an integer (nT8s)
            ierr = 0
            read(iounit, fmt=*, iostat=ierr) nT8s
            if (ierr == 0) exit
         end do
         if (failed('skip to find line starting with an integer for nT8s')) return

         allocate(T8s(nT8s), f1(4*nT8s), stat=ierr)
         if (failed('allocate')) return
         f(1:4,1:nT8s) => f1(1:4*nT8s)
         
         do j=1,nT8s
            read(iounit,'(a)',iostat=ierr) line
            if (ierr == 0) call str_to_vector(line, vec, nvec, ierr)
            if (failed('read table')) return
            if (nvec < 2) then
               ierr = -1
               if (failed('read table')) return
            end if
            T8s(j) = vec(1)
            f(1,j) = vec(2)
         end do
         
         allocate(work(nT8s*pm_work_size), stat=ierr)
         if (failed('allocate')) return

         call interp_pm(T8s, nT8s, f1, pm_work_size, work,  &
                  'rates get_interp_table', ierr)
         deallocate(work)
         
         if (failed('interp_pm')) return
         
         close(iounit)
         
         ! don't set the pointers until have finished setting up the data
         
         if (associated(T8s_out)) deallocate(T8s_out)
         if (associated(f1_out)) deallocate(f1_out)

         T8s_out => T8s
         f1_out => f1
         
         contains

         logical function failed(str)
            character (len=*), intent(in) :: str
            failed = .false.
            if (ierr == 0) return
            close(iounit)
            write(*,*) trim(str) // ' failed in reading ' // trim(rate_file)
            failed = .true.
            return
         end function failed
      
      end subroutine get_interp_table
      

      subroutine get_reaclib_rate_and_dlnT( &
            ir, temp, lambda, dlambda_dlnT, rlambda, drlambda_dlnT, ierr)
         use ratelib
         integer, intent(in) :: ir
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: lambda, dlambda_dlnT, rlambda, drlambda_dlnT
         integer, intent(out) :: ierr
         integer :: reverse
         
         include 'formats'
         ierr = 0
         
         reverse = reaction_is_reverse(ir)
         if (reverse == 0) then ! that means don't know
            reverse = reaclib_reverse(reaction_Name(ir))
            if (reverse > 0) then
               reaction_is_reverse(ir) = reverse
            else
               reaction_is_reverse(ir) = -1
            end if
         else if (reverse == -1) then ! that means reaclib_reverse returned 0
            reverse = 0
         end if
         if (reverse > 0) then
            !write(*,1) 'call do_jina_reaclib_reverse ' // trim(reaction_Name(ir))
            call do_jina_reaclib_reverse(reaclib_rates% reaction_handle(reverse))
         else
            !write(*,1) 'call do_jina_reaclib ' // trim(reaction_Name(ir))
            call do_jina_reaclib
         end if 
         
         return
         
         write(*,1) 'temp', temp
         write(*,1) 'lambda', lambda
         write(*,1) 'dlambda_dlnT', dlambda_dlnT
         write(*,1) 'rlambda', rlambda
         write(*,1) 'drlambda_dlnT', drlambda_dlnT
         
         stop 'get_reaclib_rate_and_dlnT'
         
         contains
         
         
         subroutine get_reaclib_lo_hi(ir, handle, lo, hi, ierr)
         use reaclib_eval, only: do_reaclib_indices_for_reaction
            integer, intent(in) :: ir
            character (len=*) :: handle
            integer, intent(out) :: lo, hi, ierr
            include 'formats'
            ierr = 0
            lo = reaction_reaclib_lo(ir)
            hi = reaction_reaclib_hi(ir)
            if (lo > 0 .and. hi > 0) return
            call do_reaclib_indices_for_reaction( &
               handle, reaclib_rates, lo, hi, ierr)
            if (ierr /= 0) return
            if (lo <= 0 .or. hi <= 0) ierr = -1
            reaction_reaclib_lo(ir) = lo
            reaction_reaclib_hi(ir) = hi
         end subroutine get_reaclib_lo_hi
         
         
         subroutine do_jina_reaclib
            integer :: ierr, lo, hi
            include 'formats'
            ierr = 0
            call get_reaclib_lo_hi(ir, reaction_Name(ir), lo, hi, ierr)
            if (ierr /= 0) then
               write(*,'(a,3x,i5)')  &
                  trim(reaction_Name(ir)) // ' failed in do_jina_reaclib', ir
               !stop 'raw_rates'
               return
            end if
            !write(*,3) trim(reaction_Name(ir)) // ' lo hi', lo, hi
            call reaclib_rate_and_dlnT( &
               lo, hi, reaction_Name(ir), temp*1d-9,  &
               lambda, dlambda_dlnT, rlambda, drlambda_dlnT, ierr)
            !write(*,2) trim(reaction_Name(ir)) // ' lambda', ir, lambda
            if (ierr /= 0) then
               write(*,'(a,3x,i5)')  &
                  trim(reaction_Name(ir)) // ' failed in get_reaclib_rate_and_dlnT', ir
               return
            end if
         end subroutine do_jina_reaclib
         
         
         subroutine do_jina_reaclib_reverse(reverse_handle)
            character (len=*) :: reverse_handle
            integer :: ierr, lo, hi, r_id
            include 'formats.dek'
            ierr = 0
            r_id = reverse_reaction_id(ir)
            if (r_id == 0) then ! don't know
               r_id = get_rates_reaction_id(reverse_handle)
               if (r_id == 0) then 
                  write(*,'(a,3x,i5)')  &
                     trim(reverse_handle) // ' failed in reaclib_index', r_id
                  !stop 'raw_rates'
               end if
               reverse_reaction_id(ir) = r_id
            end if
            call get_reaclib_lo_hi(r_id, reverse_handle, lo, hi, ierr)
            if (ierr /= 0) then
               write(*,'(a,3x,i5)')  &
                  trim(reverse_handle) // ' failed in do_jina_reaclib_reverse', r_id
               stop 'raw_rates'
               return
            end if
            call reaclib_rate_and_dlnT( &
               lo, hi, reverse_handle, temp*1d-9,  &
               rlambda, drlambda_dlnT, lambda, dlambda_dlnT, ierr)
            if (ierr /= 0) then
               write(*,'(a)') trim(reverse_handle) // ' failed in get_reaclib_rate_and_dlnT'
               return
            end if
         end subroutine do_jina_reaclib_reverse

                
      end subroutine get_reaclib_rate_and_dlnT




      end module raw_rates

