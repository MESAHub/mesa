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




! last checked for consistency with Frank's ratlib.f on Sept 18, 2008

      module ratelib
      use rates_def
      use utils_lib
      use chem_def !, only: nuclide_data, chem_isos
      use chem_lib, only: chem_get_iso_id
      use const_def, only: dp, ln2
      use math_lib
      
      implicit none

      real(dp), parameter :: lowT9_cutoff = 1d-3 ! all non-pp rates except decays go to 0 below this


      
      real(dp), parameter :: lowT9pp_cutoff = 1d-5 ! all pp rates except decays go to 0 below this

      real(dp)  oneth, twoth, fourth, fiveth, elvnth, fivfour, onesix, &
                        fivsix, sevsix, onefif, sixfif, onesev, twosev, foursev
      parameter        (oneth   = 1.0d0/3.0d0,  &
                        twoth   = 2.0d0/3.0d0,  &
                        fourth  = 4.0d0/3.0d0,  &
                        fiveth  = 5.0d0/3.0d0,  &
                        elvnth  = 11.0d0/3.0d0,  &
                        fivfour = 1.25d0,  &
                        onesix  = 1.0d0/6.0d0,  &
                        fivsix  = 5.0d0/6.0d0,  &
                        sevsix  = 7.0d0/6.0d0,  &
                        onefif  = 0.2d0,  &
                        sixfif  = 1.2d0,  &
                        onesev  = 1.0d0/7.0d0,  &
                        twosev  = 2.0d0/7.0d0,  &
                        foursev = 4.0d0/7.0d0)

      integer, parameter :: nTs_rf18ap = 13
      logical :: have_f_rf18ap = .false.
      real(dp), target :: f_rf18ap_ary(4*nTs_rf18ap)
      real(dp), pointer :: f_rf18ap1(:), f_rf18ap(:,:)

      contains
   ! sources:
   ! cf88      Caughlin, G. R. & Fowler, W. A. 1988, Atom. Data and Nuc. Data Tables, 40, 283
   ! nacre     C. Angulo et al., Nucl. Phys. A656 (1999)3-187
   ! wk82      wiescher and kettner, ap. j., 263, 891 (1982)
   ! c96       champagne 1996

      

! Hydrogen

! rpp, p(p,e+nu)h2
      
      subroutine rate_pp_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) :: term,aa,bb

         if (tf% t9 < lowT9pp_cutoff) then
            fr = 0; rr = 0; return
         end if 
                    
         if (tf% t9 .le. 3d0) then
            aa   = 4.01d-15 * tf% t9i23 * exp(-3.380d0*tf% t9i13) 
            bb   = 1.0d0 + 0.123d0*tf% t913 + 1.09d0*tf% t923 + 0.938d0*tf% t9 
            term    = aa * bb
         else
            term    = 1.1581136d-15
         end if
         fr    = term 
         rr    = 0.0d0
      end subroutine rate_pp_fxt
      
      
      subroutine rate_pp_nacre(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) :: term
         ! a0 T9i23 exp(-a1 T9i13 - (T9*a2)^2) 
         !     * (1 + b0 T9 + b1 T92 + b2 T93 + b3 T94 + b4 T95) 
         ! + c0 T9i32 exp(-c1/T9)
         ! + d0 T9i32 exp(-d1/T9)
         ! + e0 T9^e1 exp(-e2/T9)

         if (tf% t9 < lowT9pp_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         call rnacre(tf,  &
            4.08d-15, 3.381d0, 0d0, & ! a0, a1, a2
            3.82d0, 1.51d0, 0.144d0, -1.14d-02, 0d0, & ! b0, b1, b2, b3, b4
            0d0, 0d0, &  ! c0, c1
            0d0, 0d0, & ! d0, d1
            0d0, 0d0, 0d0, & ! e0, e1, e2
            term)         
         fr    = term 
         rr    = 0.0d0
      end subroutine rate_pp_nacre
      

      subroutine rate_pp_jina(tf, temp, fr, rr) ! cf88
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         integer :: ierr
         include 'formats.dek'
         ierr = 0
!         p    p    d                       bet+w     1.44206d+00          
         call reaclib_rate_for_handle('r_h1_h1_wk_h2', tf% T9, fr, rr, ierr)
         if (ierr /= 0) then
            write(*,'(a)') 'failed to get reaclib rate r_h1_h1_wk_h2'
            call rate_pp_fxt(tf, temp, fr, rr)
         end if
      end subroutine rate_pp_jina
              
              
! rpep, p(e-p, nu)h2
   
      subroutine rate_pep_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, aa, bb

         if (tf% t9 < lowT9pp_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         if ((tf% T9)  <=  3d0) then
            aa   = 1.36d-20 * (tf% T9i76) * exp(-3.380d0*(tf% T9i13))
            bb   = (1.0d0 - 0.729d0*(tf% T913) + 9.82d0*(tf% T923))
            term    = aa * bb
         else
            term    = 7.3824387d-21
         end if
         fr    = term 
         rr    = 0.0d0
      end subroutine rate_pep_fxt
      

      subroutine rate_pep_jina(tf, temp, fr, rr) ! cf88
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         integer :: ierr
         ierr = 0
!         p    p    d                         ecw     1.44206d+00          
         call reaclib_rate_for_handle('r_h1_h1_ec_h2', tf% T9, fr, rr, ierr)
         if (ierr /= 0) then
            write(*,'(a)') 'failed to get reaclib rate r_h1_h1_ec_h2'
            call rate_pep_fxt(tf, temp, fr, rr)
         end if
      end subroutine rate_pep_jina
              
              
! rdpg, h2(p,g)he3


      subroutine rate_dpg_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) :: term, rev, aa, bb
         include 'formats.dek'

         if (tf% t9 < lowT9pp_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         aa      = 2.24d+03 * tf% t9i23 * exp(-3.720d0*tf% t9i13) 
         bb      = 1.0d0 + 0.112d0*tf% t913 + 3.38d0*tf% t923 + 2.65d0*tf% t9
         term    = aa * bb 
         fr    = term 
         rev      = 1.63d+10 * tf% t932 * exp(-63.750d0*tf% t9i)
         rr    = rev * term
         !if (temp > 3.1d6 .and. temp < 3.2d6) write(*,1) 'rates dpg', fr, temp
      end subroutine rate_dpg_fxt
      

      subroutine rate_dpg_nacre(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, dtermdt, rev, drevdt, aa, daa, bb, dbb

         if (tf% t9 < lowT9pp_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         if (tf% T9 <= 0.11d0) then
            ! a0 T9i23 exp(-a1 T9i13 - (T9*a2)^2) 
            !     * (1 + b0 T9 + b1 T92 + b2 T93 + b3 T94 + b4 T95) 
            ! + c0 T9i32 exp(-c1/T9)
            ! + d0 T9i32 exp(-d1/T9)
            ! + e0 T9^e1 exp(-e2/T9)
            call rnacre(tf,  &
               1.81d3, 3.721d0, 0d0, & ! a0, a1, a2
               14.3d0, -90.5d0, 395d0, 0d0, 0d0, & ! b0, b1, b2, b3, b4
               0d0, 0d0, &  ! c0, c1
               0d0, 0d0, & ! d0, d1
               0d0, 0d0, 0d0, & ! e0, e1, e2
               term)         
         else
            ! a0 T9i23 exp(-a1 T9i13 - (T9*a2)^2) 
            !     * (1 + b0 T9 + b1 T92 + b2 T93 + b3 T94 + b4 T95) 
            ! + c0 T9i32 exp(-c1/T9)
            ! + d0 T9i32 exp(-d1/T9)
            ! + e0 T9^e1 exp(-e2/T9)
            call rnacre(tf,  &
               2.58d3, 3.721d0, 0d0, & ! a0, a1, a2
               3.96d0, 0.116d0, 0d0, 0d0, 0d0, & ! b0, b1, b2, b3, b4
               0d0, 0d0, &  ! c0, c1
               0d0, 0d0, & ! d0, d1
               0d0, 0d0, 0d0, & ! e0, e1, e2
               term)               
         end if
         call rnacre_rev(tf, &  ! a0 T932 exp(-a1/T9)
            1.63d10, 63.749d0, &  ! a0, a1
            rev)     
         fr    = term 
         rr    = rev * term
      end subroutine rate_dpg_nacre
      

      subroutine rate_dpg_jina(tf, temp, fr, rr) ! cf88
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!         p    d  he3                       de04      5.49300d+00          
         call jina_reaclib_2_1(ih1, ih2, ihe3, tf, fr, rr, 'rate_dpg_jina')
      end subroutine rate_dpg_jina


      subroutine rate_png_fxt(tf, temp, fr, rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
      real(dp) term, rev, aa

         if (tf% t9 < lowT9pp_cutoff) then
            fr = 0; rr = 0; return
         end if 
          


! p(n, g)d
! smith, kawano, malany 1992

       aa      = 1.0d0 - 0.8504d0*(tf% T912) + 0.4895d0*(tf% T9) &
                 - 0.09623d0*(tf% T932) + 8.471d-3*(tf% T92)  &
                 - 2.80d-4*(tf% T952)
  
       term    = 4.742d4 * aa

! wagoner, schramm 1977
!      aa      = 1.0d0 - 0.86d0*(tf% T912) + 0.429d0*(tf% T9)
!      daa     =  -0.5d0*0.86d0*(tf% T9i12) + 0.429

!      term    = 4.4d4 * aa
!      dtermdt = 4.4d4 * daa

      fr    = term 
      rev      = 4.71d+09 * (tf% T932) * exp(-25.82d0*(tf% T9i))
      rr    = rev * term

      end subroutine rate_png_fxt
      

      

      subroutine rate_ddg_jina(tf, temp, fr, rr) ! cf88
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!         d    d  he4                       cf88n     2.38470d+01          
         call jina_reaclib_2_1(ih2, ih2, ihe4, tf, fr, rr, 'rate_ddg_jina')
      end subroutine rate_ddg_jina


! Helium

! rhe3p, he3(p,e+nu)he4
      

      subroutine rate_hep_jina(tf, temp, fr, rr) ! cf88
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         integer :: ierr
         ierr = 0
!         p  he3  he4                       bet+w     1.97960d+01          
         call reaclib_rate_for_handle('r_h1_he3_wk_he4', tf% T9, fr, rr, ierr)
         if (ierr /= 0) then
            write(*,'(a)') 'failed to get reaclib rate r_h1_he3_wk_he4'
            call rate_hep_fxt(tf, temp, fr, rr)
         end if
      end subroutine rate_hep_jina


      subroutine rate_hep_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, aa

         if (tf% t9 < lowT9pp_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         if ((tf% T9)  <=  3d0) then
            aa   = 8.78d-13 * (tf% T9i23) * exp(-6.141d0*(tf% T9i13))
            term    = aa
         else
            term    = 5.9733434d-15
         end if
         fr    = term 
         rr    = 0.0d0
      end subroutine rate_hep_fxt


! rhe3d     he3(d,p)he4    de04

      subroutine rate_he3d_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!         d  he3    p  he4                  de04      1.83530d+01          
         call jina_reaclib_2_2(ih2, ihe3, ih1, ihe4, tf, fr, rr, 'rate_he3d_jina')
      end subroutine rate_he3d_jina



! r33, he3(he3, 2p)he4       

      subroutine rate_he3he3_nacre(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, rev
         ! a0 T9i23 exp(-a1 T9i13 - (T9*a2)^2) 
         !     * (1 + b0 T9 + b1 T92 + b2 T93 + b3 T94 + b4 T95) 
         ! + c0 T9i32 exp(-c1/T9)
         ! + d0 T9i32 exp(-d1/T9)
         ! + e0 T9^e1 exp(-e2/T9)

         if (tf% t9 < lowT9pp_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         call rnacre(tf,  &
            5.59d10, 12.277d0, 0d0, & ! a0, a1, a2
            -0.135d0, 2.54d-2, -1.29d-03, 0d0, 0d0, & ! b0, b1, b2, b3, b4
            0d0, 0d0, &  ! c0, c1
            0d0, 0d0, & ! d0, d1
            0d0, 0d0, 0d0, & ! e0, e1, e2
            term)         
         call rnacre_rev(tf, &  ! a0 T932 exp(-a1/T9)
            3.392d-10, 149.23d0, &  ! a0, a1
            rev)     
         fr    = term 
         rr    = rev * term 
      end subroutine rate_he3he3_nacre

      subroutine rate_he3he3_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, rev
!       he3  he3  h1 h1 he4         
         call jina_reaclib_2_3(ihe3, ihe3, ih1, ih1, ihe4, tf, fr, rr, 'rate_he3he3_jina')
      end subroutine rate_he3he3_jina


! r34, he4(he3,g)be7 


      subroutine rate_he3he4_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!       he4  he3  be7                       de04      1.58700d+00          
         call jina_reaclib_2_1(ihe4, ihe3, ibe7, tf, fr, rr, 'rate_he3he4_jina')
      end subroutine rate_he3he4_jina
      

      subroutine rate_he3he4_nacre(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, rev
         ! a0 T9i23 exp(-a1 T9i13 - (T9*a2)^2) 
         !     * (1 + b0 T9 + b1 T92 + b2 T93 + b3 T94 + b4 T95) 
         ! + c0 T9i32 exp(-c1/T9)
         ! + d0 T9i32 exp(-d1/T9)
         ! + e0 T9^e1 exp(-e2/T9)

         if (tf% t9 < lowT9pp_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         call rnacre(tf,  &
            5.46d6, 12.827d0, 0d0, & ! a0, a1, a2
            -0.307d0, 8.81d-2, -1.06d-2, 4.46d-4, 0d0, & ! b0, b1, b2, b3, b4
            0d0, 0d0, &  ! c0, c1
            0d0, 0d0, & ! d0, d1
            0d0, 0d0, 0d0, & ! e0, e1, e2
            term)        
         call rnacre_rev(tf, &  ! a0 T932 exp(-a1/T9)
            1.113d10, 18.412d0, &  ! a0, a1
            rev)     
         fr    = term
         rr    = rev * term 
      end subroutine rate_he3he4_nacre



! r3a, triple alpha


      subroutine rate_tripalf_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
      real(dp) :: fr1, rr1
      include 'formats.dek'
      
!       he4  he4  he4  c12                  fy05r     7.27500d+00          
         call jina_reaclib_3_1(ihe4, ihe4, ihe4, ic12, tf, fr, rr, 'rate_tripalf_jina')
         
         return 
         call rate_tripalf_reaclib(tf, temp, fr1, rr1)
         
         write(*,1) 'fr', fr
         write(*,1) 'fr1', fr1
         write(*,1) 'rr', rr
         write(*,1) 'rr1', rr1
         write(*,*)
         stop 'rate_tripalf_jina' 
      end subroutine rate_tripalf_jina



      subroutine rate_tripalf_reaclib(tf, temp, fr, rr)
      !      HE4(2A,G)C12    reaclib JINA - Fynbo et al. 2005 Nature 433, 136-139
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
      real(dp) :: fr1, fr2, fr3, rev
      rr = 0
      if (temp < 1d7) then
         fr = 0; return
      end if
      call do_reaclib(tf,   &
                  -9.710520d-01,  0.000000d+00, -3.706000d+01,  2.934930d+01,                        &
                  -1.155070d+02, -1.000000d+01, -1.333330d+00,                                     &
                  fr1)
      call do_reaclib(tf,   &
                  -2.435050d+01, -4.126560d+00, -1.349000d+01,  2.142590d+01,                        &
                  -1.347690d+00,  8.798160d-02, -1.316530d+01,                                     &
                  fr2)
      call do_reaclib(tf,   &
                  -1.178840d+01, -1.024460d+00, -2.357000d+01,  2.048860d+01,                        &
                  -1.298820d+01, -2.000000d+01, -2.166670d+00,                                     &
                  fr3)
         fr = fr1 + fr2 + fr3         
         ! use the fxt reverse rate term
         rev    = 2.00d+20*(tf% t93)*exp(-84.424d0*(tf% t9i))
         rr = fr * rev         
      end subroutine rate_tripalf_reaclib

      

      subroutine rate_tripalf_nacre(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) :: r2abe, rbeac, bb, term, rev
         ! he4(a, g)be8
         ! a0 T9i23 exp(-a1 T9i13 - (T9*a2)^2) 
         !     * (1 + b0 T9 + b1 T92 + b2 T93 + b3 T94 + b4 T95) 
         ! + c0 T9i32 exp(-c1/T9)
         ! + d0 T9i32 exp(-d1/T9)
         ! + e0 T9^e1 exp(-e2/T9)

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         call rnacre(tf,  &
               2.43d9, 13.490d0, 1d0/0.15d0, & ! a0, a1, a2
               74.5d0, 0d0, 0d0, 0d0, 0d0, & ! b0, b1, b2, b3, b4
               6.09d5, 1.054d0, &  ! c0, c1
               0d0, 0d0, & ! d0, d1
               0d0, 0d0, 0d0, & ! e0, e1, e2
               r2abe)         
        ! be8(a, g)c12
         ! a0 T9i23 exp(-a1 T9i13 - (T9*a2)^2) 
         !     * (1 + b0 T9 + b1 T92 + b2 T93 + b3 T94 + b4 T95) 
         ! + c0 T9i32 exp(-c1/T9)
         ! + d0 T9i32 exp(-d1/T9)
         ! + e0 T9^e1 exp(-e2/T9)
         call rnacre(tf,  &
               2.76d7, 23.570d0, 1d0/0.4d0, & ! a0, a1, a2
               5.47d0, 326d0, 0d0, 0d0, 0d0, & ! b0, b1, b2, b3, b4
               130.7d0, 3.338d0, &  ! c0, c1
               2.51d4, 20.307d0, & ! d0, d1
               0d0, 0d0, 0d0, & ! e0, e1, e2
               rbeac)               
         if (tf% T9 <= 0.03d0) then
            bb    = 3.07d-16*(1d0 - 29.1d0*(tf% T9) + 1308d0*(tf% T92))
            if (bb < 0) then
               bb = 0
            end if      
         else
            bb    = 3.44d-16*(1 +     0.0158d0*pow(tf% T9,-0.65d0))
         end if      
         term    = r2abe * rbeac * bb
         call rnacre_rev(tf, &  ! a0 T932 exp(-a1/T9)
            2.003d20, 84.415d0, &  ! a0, a1
            rev)     
         fr    = term
         rr    = rev * term
      end subroutine rate_tripalf_nacre


      subroutine rate_tripalf_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr

      real(dp) term, dtermdt, rev, drevdt, r2abe, dr2abedt, rbeac,  &
                       drbeacdt, aa, daa, bb, dbb, cc, dcc, dd, ddd, ee, dee,  &
                       ff, dff, xx, dxx, yy, dyy, zz, dzz, uu, vv, f1, df1, rc28,  &
                       q1, q2
      parameter        (rc28   = 0.1d0,  &
                        q1     = 1.0d0/0.009604d0,  &
                        q2     = 1.0d0/0.055225d0) 


         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
! this is a(a, g)be8
         aa    = 7.40d+05 * (tf% t9i32) * exp(-1.0663d0*(tf% t9i)) 

         bb    = 4.164d+09 * (tf% t9i23) * exp(-13.49d0*(tf% t9i13) - (tf% t92)*q1)

         cc    = 1.0d0 + 0.031d0*(tf% t913) + 8.009d0*(tf% t923) + 1.732d0*(tf% t9)   &
              + 49.883d0*(tf% t943) + 27.426d0*(tf% t953)

         r2abe    = aa + bb * cc

! this is be8(a, g)c12
         dd    = 130.0d0 * (tf% t9i32) * exp(-3.3364d0*(tf% t9i))  

         ee    = 2.510d+07 * (tf% t9i23) * exp(-23.57d0*(tf% t9i13) - (tf% t92)*q2)

         ff    = 1.0d0 + 0.018d0*(tf% t913) + 5.249d0*(tf% t923) + 0.650d0*(tf% t9) +  &
              19.176d0*(tf% t943) + 6.034d0*(tf% t953)

         rbeac    = dd + ee * ff

! a factor
         xx    = rc28 * 1.35d-07 * (tf% t9i32) * exp(-24.811d0*(tf% t9i))


! high temperature rate
         if ((tf% t9).gt.0.08d0) then
          term    = 2.90d-16 * r2abe * rbeac + xx


! low temperature rate
         else
          uu   = 0.8d0*exp(-pow(0.025d0*(tf% t9i),3.263d0)) 
          yy   = 0.2d0 + uu  ! fixes a typo in Frank's original
          vv   = 4.0d0*exp(-pow((tf% t9)/0.025d0,9.227d0)) 
          zz   = 1.0d0 + vv
          aa   = 1.0d0/zz
          f1   = 0.01d0 + yy * aa  ! fixes a typo in Frank's original
          term = 2.90d-16 * r2abe * rbeac * f1 +  xx 
         end if

         rev    = 2.00d+20*(tf% t93)*exp(-84.424d0*(tf% t9i))

!      term    = 1.2d0 * term
!      dtermdt = 1.2d0 * term

      fr    = term

      rr    = rev * term

      end subroutine rate_tripalf_fxt


      subroutine rate_he3ng_fxt(tf, temp, fr, rr)

      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr

      real(dp) term, rev


         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
! he3(n, g)he4
      term    = 6.62d0 * (1.0d0 + 905.0d0*(tf% T9))
      fr    = term 
      rev      = 2.61d+10 * (tf% T932) * exp(-238.81d0*(tf% T9i))
      rr    = rev * term
      end subroutine rate_he3ng_fxt
      
      


! Lithium  


! rli7pa, li7(p,a)he4

      subroutine rate_li7pa_nacre(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, rev
         ! a0 T9i23 exp(-a1 T9i13 - (T9*a2)^2) 
         !     * (1 + b0 T9 + b1 T92 + b2 T93 + b3 T94 + b4 T95) 
         ! + c0 T9i32 exp(-c1/T9)
         ! + d0 T9i32 exp(-d1/T9)
         ! + e0 T9^e1 exp(-e2/T9)

         if (tf% t9 < lowT9pp_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         call rnacre(tf,  &
            7.20d8, 8.473d0, 1d0/6.5d0, & ! a0, a1, a2
            1.05d0, -0.653d0, 0.185d0, -2.12d-2, 9.30d-4, & ! b0, b1, b2, b3, b4
            0d0, 0d0, & ! c0, c1
            0d0, 0d0, & ! d0, d1
            9.85d6, 0.576d0, 10.415d0, & ! e0, e1, e2
            term)              
         call rnacre_rev(tf, &  ! a0 T932 exp(-a1/T9)
            4.676d0, 201.30d0, &  ! a0, a1
            rev)     
         fr    = term
         rr    = rev * term
      end subroutine rate_li7pa_nacre


      
      
      subroutine rate_li7pa_jina(tf, temp, fr, rr) ! jina reaclib
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         call jina_reaclib_2_2(ih1, ili7, ihe4, ihe4, tf, fr, rr, 'rate_li7pa_jina')
      end subroutine rate_li7pa_jina


! rli7pg, li7(p,g)be8 => 2 he4  


! Beryllium 

! rbe7ec, be7(e-, nu)li7

      subroutine rate_be7em_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, dtermdt, aa, daa, bb, dbb

         if (tf% t9 < lowT9pp_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         if (tf% T9 <= 3d0 .and. tf% T9 >= 1d-3) then
            aa  = 0.0027d0*(tf% T9i) * exp(2.515d-3*(tf% T9i)) 
            bb  = 1.0d0 - 0.537d0*(tf% T913) + 3.86d0*(tf% T923) + aa
            term    = 1.34d-10 * (tf% T9i12) * bb
         else
            term    = 0.0d0
         endif
         fr    = term
         rr    = 0.0d0
      end subroutine rate_be7em_fxt


      subroutine rate_be7em_jina(tf, temp, fr, rr)        
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         integer :: ierr
         ierr = 0
         call reaclib_rate_for_handle('r_be7_wk_li7', tf% T9, fr, rr, ierr)
         if (ierr /= 0) then
            write(*,'(a)') 'failed to get reaclib rate r_be7_wk_li7'
            call rate_b8ep(tf, temp, fr, rr)
         end if
      end subroutine rate_be7em_jina

       
! rbe7pg, be7(p,g)b8

      subroutine rate_be7pg_nacre(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, rev, aa, bb, cc
         ! a0 T9i23 exp(-a1 T9i13 - (T9*a2)^2) 
         !     * (1 + b0 T9 + b1 T92 + b2 T93 + b3 T94 + b4 T95) 
         ! + c0 T9i32 exp(-c1/T9)
         ! + d0 T9i32 exp(-d1/T9)
         ! + e0 T9^e1 exp(-e2/T9)

         if (tf% t9 < lowT9pp_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         call rnacre(tf,  &
            2.61d5, 10.264d0, 0d0, & ! a0, a1, a2
            -5.11d-2, 4.68d-2, -6.60d-3, 3.12d-4, 0d0, & ! b0, b1, b2, b3, b4
            2.05d3, 7.345d0, &  ! c0, c1
            0d0, 0d0, & ! d0, d1
            0d0, 0d0, 0d0, & ! e0, e1, e2
            term)              
         call rnacre_rev(tf, &  ! a0 T932 exp(-a1/T9)
            1.306d10, 1.594d0, &  ! a0, a1
            rev)     
         fr    = term 
         rr    = rev * term
      end subroutine rate_be7pg_nacre
      



      subroutine rate_be7pg_jina(tf, temp, fr, rr) ! jina reaclib   cf88
!         p  be7   b8                       cf88n     1.37000d-01          
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         call jina_reaclib_2_1(ih1, ibe7, ib8, tf, fr, rr, 'rate_be7pg_jina')
      end subroutine rate_be7pg_jina


! rbe7dp    be7(d,p)2he4

      subroutine rate_be7dp_jina(tf, temp, fr, rr)
!         d  be7    p  he4  he4             cf88n     1.67660d+01          
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!
         call jina_reaclib_2_3( &
                  ih2, ibe7, ih1, ihe4, ihe4, tf, fr, rr, 'rate_be7pg_jina')
      end subroutine rate_be7dp_jina


! rbe7dp    be7(he3,2p)2he4    


      subroutine rate_be7he3_jina(tf, temp, fr, rr)
!       he3  be7    p    p  he4  he4        mafon     1.12721d+01          
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!
         call jina_reaclib_2_4( &
               ihe3, ibe7, ih1, ih1, ihe4, ihe4, tf, fr, rr, 'rate_be7he3_jina')
      end subroutine rate_be7he3_jina



! be9(p,d)be8 => 2a



! Boron

! rb8ep, b8(e+, nu)be8 => 2a    

      subroutine rate_b8ep(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp), parameter :: halflife = 0.77d0 ! 770 ms
         rr    = 0.0d0
         fr    = ln2/halflife
         rr = 0
      end subroutine rate_b8ep

      subroutine rate_b8_wk_he4_he4_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         integer :: ierr
         ierr = 0
         call reaclib_rate_for_handle('r_b8_wk_he4_he4', tf% T9, fr, rr, ierr)
         if (ierr /= 0) then
            write(*,'(a)') 'failed to get reaclib rate r_b8_wk_he4_he4'
            call rate_b8ep(tf, temp, fr, rr)
         end if
      end subroutine rate_b8_wk_he4_he4_jina
      

! rb8gp, b8(g,p)be7
      ! see rbe7pg
      

! Carbon

! rc12pg, c12(p,g)n13


      subroutine rate_c12pg_nacre(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, rev
         ! a0 T9i23 exp(-a1 T9i13 - (T9*a2)^2) 
         !     * (1 + b0 T9 + b1 T92 + b2 T93 + b3 T94 + b4 T95) 
         ! + c0 T9i32 exp(-c1/T9)
         ! + d0 T9i32 exp(-d1/T9)
         ! + e0 T9^e1 exp(-e2/T9)

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         call rnacre(tf,  &
            2.00d7, 13.692d0, 1d0/0.46d0, & ! a0, a1, a2
            9.89d0, -59.8d0, 266d0, 0d0, 0d0, & ! b0, b1, b2, b3, b4
            1.00d5, 4.913d0, &  ! c0, c1
            4.24d5, 21.62d0, & ! d0, d1
            0d0, 0d0, 0d0, & ! e0, e1, e2
            term)      
         call rnacre_rev(tf, &  ! a0 T932 exp(-a1/T9)
            8.847d9, 22.553d0, &  ! a0, a1
            rev)     
         fr    = term 
         rr    = rev * term 
      end subroutine rate_c12pg_nacre


      subroutine rate_c12pg_jina(tf, temp, fr, rr) ! jina reaclib   nacre
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!         p  c12  n13                       nacrn     1.94300d+00
         call jina_reaclib_2_1(ih1, ic12, in13, tf, fr, rr, 'rate_c12pg_jina')
      end subroutine rate_c12pg_jina
    
! rc12ap, c12(a,p)n15


      subroutine rate_n15pa_nacre(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, rev

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         if (tf% T9 <= 2.5d0) then
            ! a0 T9i23 exp(-a1 T9i13 - (T9*a2)^2) 
            !     * (1 + b0 T9 + b1 T92 + b2 T93 + b3 T94 + b4 T95) 
            ! + c0 T9i32 exp(-c1/T9)
            ! + d0 T9i32 exp(-d1/T9)
            ! + e0 T9^e1 exp(-e2/T9)
            call rnacre(tf,  &
               1.12d12, 15.253d0, 1d0/0.28d0, & ! a0, a1, a2
               4.95d0, 143d0, 0d0, 0d0, 0d0, & ! b0, b1, b2, b3, b4
               1.01d8, 3.643d0, &  ! c0, c1
               1.19d9, 7.406d0, & ! d0, d1
               0d0, 0d0, 0d0, & ! e0, e1, e2
               term)               
         else
            term     = 4.17d7 * pow((tf% T9),0.917d0) * exp(-3.292d0*(tf% T9i)) 
         end if
         call rnacre_rev(tf, &  ! a0 T932 exp(-a1/T9)
            7.06d-1, 57.622d0, &  ! a0, a1
            rev)     
         fr    = term 
         rr    = rev * term 
      end subroutine rate_n15pa_nacre

      subroutine rate_n15pa_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!         p  n15  he4  c12                  nacrr     4.96600d+00          
         call jina_reaclib_2_2(ih1, in15, ihe4, ic12, tf, fr, rr, 'rate_n15pa_jina')
      end subroutine rate_n15pa_jina


! rc12ag, c12(a,g)o16          
      subroutine rate_c12ag_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr

      real(dp) term, dtermdt, rev, drevdt, aa, daa, bb, dbb, cc, dcc,  &
                       dd, ddd, ee, dee, ff, dff, gg, dgg, hh, dhh, f1, df1, f2, df2,  &
                       zz, q1, termE1, dtermE1, termE2, dtermE2, termRes, dtermRes
      parameter        (q1 = 1.0d0/12.222016d0)

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
      aa   = 1.0d0 + 0.0489d0*(tf% t9i23)
      bb   = (tf% t92)*aa*aa
      cc   = exp(-32.120d0*(tf% t9i13) - (tf% t92)*q1)
      dd   = 1.0d0 + 0.2654d0*(tf% t9i23)
      ee   = (tf% t92)*dd*dd
      ff   = exp(-32.120d0*(tf% t9i13))
      gg   = 1.25d3 * (tf% t9i32) * exp(-27.499d0*(tf% t9i))
      hh   = 1.43d-2 * (tf% t95) * exp(-15.541d0*(tf% t9i))
      zz   = 1.0d0/bb
      f1   = cc*zz
      zz   = 1.0d0/ee
      f2   = ff*zz
      term    = 1.04d8*f1  + 1.76d8*f2 + gg + hh
! 1.7 times cf88 value
      term     = 1.7d0 * term
      rev    = 5.13d10 * (tf% t932) * exp(-83.111d0*(tf% t9i))
      fr    = term
      rr     = rev * term
      end subroutine rate_c12ag_fxt




      subroutine rate_c12ag_nacre(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, rev, termE1, termE2, termRes, aa, bb, cc
         include 'formats.dek'

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         ! note: uses T9i2 instead of T9i23, so special case it.
         aa   = 6.66d7 * (tf% T9i2) * exp(-32.123d0*(tf% T9i13) - (tf% T92)/(4.6d0*4.6d0)) 
         bb   = 1 + 2.54d0*(tf% T9) + 1.04d0*(tf% T92) -  0.226d0*(tf% T93) 
         if (bb < 0) bb = 0
         cc    = 1.39d3 * (tf% T9i32) * exp(-28.930d0*(tf% T9i)) 
         termE1    = aa * bb + cc
         aa   = 6.56d7 * (tf% T9i2) * exp(-32.123d0*(tf% T9i13) - (tf% T92)/(1.3d0*1.3d0)) 
         bb   = 1 + 9.23d0*(tf% T9) - 13.7d0*(tf% T92) +  7.4d0*(tf% T93) 
         termE2    = aa * bb
         termRes   = 19.2d0 * (tf% T92) * exp(-26.9d0*(tf% T9i)) 
         term    = termE1 + termE2 + termRes
         rev      = 5.132d10 * (tf% T932) * exp(-83.109d0*(tf% T9i))
         fr    = term
         rr     = rev * term
      end subroutine rate_c12ag_nacre
      

      subroutine rate_c12ag_kunz(tf, temp, fr, rr)
         ! kunz et al (2002)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) :: term, rev, aa, bb, cc, dd, ee
         real(dp), parameter :: &
            a0 = 1.21d8, a1 = 6.06d-2, a2 = 32.12d0, a3 = 1.7d0, a4 = 7.4d8, &
            a5 = 0.47d0, a6 = 32.12d0, a9tilda = 3.06d10, a11 = 38.534d0   
         include 'formats.dek'

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         aa   = a0 * (tf% T9i2) * exp(-a2*(tf% T9i13) - (tf% T92)/(a3*a3))
         bb   = 1 / pow(1 + a1*(tf% T9i23),2)
         cc   = a4 * (tf% T9i2) * exp(-a6*(tf% T9i13)) 
         dd   = 1 / pow(1 + a5*(tf% T9i23),2)
         ee   = a9tilda * (tf% T9i13) * exp(-a11*(tf% T9i13)) 
         term    = aa*bb + cc*dd + ee         
         rev      = 5.132d10 * (tf% T932) * exp(-83.109d0*(tf% T9i))         
         fr    = term         
         rr     = rev * term
      end subroutine rate_c12ag_kunz
      

      subroutine rate_c12ag_jina(tf, temp, fr, rr) 
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!       he4  c12  o16                       bu96n     7.16192d+00          
         call jina_reaclib_2_1(ihe4, ic12, io16, tf, fr, rr, 'rate_c12ag_jina')
      end subroutine rate_c12ag_jina


! r1212p, c12(c12,p)na23

      subroutine rate_c12c12p_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!       c12  c12    p na23                  cf88r     2.24200d+00          
         call jina_reaclib_2_2(ic12, ic12, ih1, ina23, tf, fr, rr, 'rate_c12c12p_jina')
      end subroutine rate_c12c12p_jina

! r1212a, c12(c12,a)ne20

      subroutine rate_c12c12_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) :: term, T9a, dt9a, T9a13, T9a56, aa, zz
         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
         aa      = 1.0d0 + 0.0396d0*tf% t9
         zz      = 1.0d0/aa
         t9a     = tf% t9*zz
         dt9a    = (1.0d0 -  t9a*0.0396d0)*zz
         zz      = dt9a/t9a
         t9a13   = pow(t9a,oneth)
         t9a56   = pow(t9a,fivsix)
         term    = 4.27d+26 * t9a56 * tf% t9i32 * exp(-84.165d0/t9a13 - 2.12d-03*tf% t93)
         fr    = term
         rr    = 0.0d0
         return
         end subroutine rate_c12c12_fxt


         subroutine rate_c12c12_fxt_basic(tf, temp, fr, rr)
            type (T_Factors), pointer :: tf
            real(dp), intent(in) :: temp
            real(dp), intent(out) :: fr, rr
            real(dp) :: term,t9a,t9a13,t9a56,aa,zz
            
            include 'formats.dek'

            if (tf% t9 < lowT9_cutoff) then
               fr = 0; rr = 0; return
            end if 
               
            aa      = 1.0d0 + 0.0396d0*tf% t9
            zz      = 1.0d0/aa
            t9a     = tf% t9*zz
            t9a13   = pow(t9a,oneth)
            t9a56   = pow(t9a,fivsix)
            term    = 4.27d+26 * t9a56 * tf% t9i32 * exp(-84.165d0/t9a13 - 2.12d-03*tf% t93)
            fr    = term
            rr    = 0.0d0

         end subroutine rate_c12c12_fxt_basic


         subroutine rate_c12c12_fxt_multi(tf, temp, fr, rr)
            type (T_Factors), pointer :: tf
            real(dp), intent(in) :: temp
            real(dp), intent(out) :: fr, rr
            real(dp) :: fr1, rr1, fr2, rr2
            call rate_c12c12npa(tf, temp, fr1, rr1, fr2, rr2, fr, rr)
            fr = fr + fr1 + fr2
            rr = rr + rr1 + rr2
         end subroutine rate_c12c12_fxt_multi

         subroutine rate_c12c12a_fxt(tf, temp, fr, rr)
            type (T_Factors), pointer :: tf
            real(dp), intent(in) :: temp
            real(dp), intent(out) :: fr, rr
            real(dp) :: fr1, rr1, fr2, rr2
            call rate_c12c12npa(tf, temp, fr1, rr1, fr2, rr2, fr, rr)
         end subroutine rate_c12c12a_fxt

      subroutine rate_c12c12a_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!       c12  c12  he4 ne20                  cf88r     4.62100d+00          
         call jina_reaclib_2_2(ic12, ic12, ihe4, ine20, tf, fr, rr, 'rate_c12c12a_jina')
      end subroutine rate_c12c12a_jina

      subroutine rate_c12c12npa(tf, temp,  &
         fr1, rr1, & ! c12(c12,n)mg23
         fr2, rr2, & ! c12(c12,p)na23
         fr3, rr3) ! c12(c12,a)ne20

         type (T_Factors), pointer :: tf
         real(dp) temp, fr1, rr1, fr2, rr2, fr3, rr3

         real(dp) term, rev, T9a, T9a13,  &
                  T9a56, aa, bb, cc, dd, &
                  b24n, b24p, b24a


         if (tf% t9 < lowT9_cutoff) then
         fr1 = 0; rr1 = 0
         fr2 = 0; rr2 = 0
         fr3 = 0; rr3 = 0
         return
         end if 

         aa      = 1.0d0 + 0.0396d0*(tf% T9)

         T9a     = (tf% T9)/aa

         T9a13   = pow(T9a,oneth)

         T9a56   = pow(T9a,fivsix)

         aa = 4.27d+26 * T9a56 * (tf% T9i32) *  &
         exp(-84.165d0/T9a13 - 2.12d-03*(tf% T93))

         ! neutron branching from dayras switkowski and woosley 1976
         if ((tf% T9) .ge. 1.5d0) then

         bb    =  0.055d0 * exp(0.976d0 - 0.789d0*(tf% T9))

         b24n  = 0.055d0  - bb

         else 

         bb    = 1.0d0 + 0.0789d0*(tf% T9) + 7.74d0*(tf% T92)

         cc    = 0.766d0*(tf% T9i3)

         dd    = bb * cc

         b24n  = 0.859d0*exp(-dd)

         end if

         ! proton branching ratio
         if ((tf% T9).gt.3d0) then

         b24p  = oneth*(1.0d0 - b24n)

         b24a  = 2.0d0 * b24p

         else

         b24p  = 0.5d0*(1.0d0 - b24n)

         b24a  = b24p

         end if

         ! c12(c12, n)mg23
         term    = aa * b24n
         fr1     = term

         rev    = 0.0d0
         if ((tf% T9) .gt. 0.1d0) then
         rev    = 3.93d0 * exp(30.16100515d0*(tf% T9i))
         end if
         rr1    = rev * term


         ! c12(c12, p)na23
         term    = aa * b24p
         fr2     = term

         rev    = 0.0d0
         if ((tf% T9) .gt. 0.1d0) then
         rev    = 3.93d0 * exp(-25.98325915d0*(tf% T9i))
         end if
         rr2    = rev * term


         ! c12(c12, a)ne20
         term    = aa * b24a
         fr3     = term

         rev    = 0.0d0
         if ((tf% T9) .gt. 0.1d0) then
         rev    = 2.42d0 * exp(-53.576110995d0*(tf% T9i))
         end if
         rr3    = rev * term

         end subroutine rate_c12c12npa


! r1216g, c12(o16,g)si28

      subroutine rate_c12o16_to_mg24_fxt(tf, temp, fr, rr)
      ! this is a combined rate for reactions to mg24 and si28
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
      call rate_c12o16_fxt(tf, temp, fr, rr)
      fr = 0.5d0*fr
      end subroutine rate_c12o16_to_mg24_fxt


      subroutine rate_c12o16_to_si28_fxt(tf, temp, fr, rr)
      ! this is a combined rate for reactions to mg24 and si28
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
      call rate_c12o16_fxt(tf, temp, fr, rr)
      fr = 0.5d0*fr
      end subroutine rate_c12o16_to_si28_fxt


      subroutine rate_c12o16_fxt(tf, temp, fr, rr)
      ! this is a combined rate for reactions to mg24 and si28
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
      real(dp) term, T9a, T9a13, T9a23, T9a56, aa, bb, cc

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
!  c12 + o16 reaction; see cf88 references 47-4
      if ((tf% T9).ge.0.5d0) then
       aa     = 1.0d0 + 0.055d0*(tf% T9)
       T9a    = (tf% T9)/aa
       T9a13  = pow(T9a,oneth)
       T9a23  = T9a13*T9a13
       T9a56  = pow(T9a,fivsix)
       aa      = exp(-0.18d0*T9a*T9a) 
       bb      = 1.06d-03*exp(2.562d0*T9a23)
       cc      = aa + bb
       term    = 1.72d+31 * T9a56 * (tf% T9i32) * exp(-106.594d0/T9a13)/cc
      else
!       term    = 2.6288035d-29
       term    = 0.0d0
      endif
      fr    = term
      rr    = 0.0d0      
      end subroutine rate_c12o16_fxt

      subroutine rate_c12o16_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) :: fr1, rr1
         call rate_c12o16a_jina(tf, temp, fr, rr)
         call rate_c12o16p_jina(tf, temp, fr1, rr1)
         fr = fr + fr1
         rr = rr + rr1
      end subroutine rate_c12o16_jina


      subroutine rate_c12o16p_fxt(tf, temp, fr, rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
      call rate_c12o16_fxt(tf, temp, fr, rr)
      fr = 0.5d0*fr
      end subroutine rate_c12o16p_fxt


! r1216n, c12(o16,n)si27
      subroutine rate_c12o16n(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) :: fr2, rr2, fr3, rr3
         call rate_c12o16npa(tf, temp, fr, rr, fr2, rr2, fr3, rr3)
      end subroutine rate_c12o16n

! r1216p, c12(o16,p)al27
      subroutine rate_c12o16p(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) :: fr1, rr1, fr3, rr3
         call rate_c12o16npa(tf, temp, fr1, rr1, fr, rr, fr3, rr3)
      end subroutine rate_c12o16p

      subroutine rate_c12o16p_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!       c12  o16    p al27                  cf88r     5.17100d+00          
         call jina_reaclib_2_2(ic12, io16, ih1, ial27, tf, fr, rr, 'rate_c12o16p_jina')
      end subroutine rate_c12o16p_jina

      subroutine rate_c12o16a_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!       c12  o16  he4 mg24                  cf88r     6.77100d+00          
         call jina_reaclib_2_2(ic12, io16, ihe4, img24, tf, fr, rr, 'rate_c12o16a_jina')
      end subroutine rate_c12o16a_jina


! r1216n, c12(o16,n)si27
! r1216p, c12(o16,p)al27
! r1216a, c12(o16,a)mg24

      subroutine rate_c12o16npa(tf, temp, &
         fr1, rr1, & ! c12(o16,n)si27
         fr2, rr2, & ! c12(o16,p)al27 
         fr3, rr3) ! c12(o16,a)mg24

         type (T_Factors), pointer :: tf
         real(dp) temp, fr1, rr1, fr2, rr2, fr3, rr3

         real(dp) term, rev, T9a, T9a13,  &
               T9a23, T9a56, aa, bb, cc,  &
               dd, b27n, b27p, b24a


         if (tf% t9 < lowT9_cutoff) then
         fr1 = 0; rr1 = 0
         fr2 = 0; rr2 = 0
         fr3 = 0; rr3 = 0
         return
         end if 

         if ((tf% T9).ge.0.5d0) then
         aa     = 1.0d0 + 0.055d0*(tf% T9)
         T9a    = (tf% T9)/aa
         t9a13  = pow(t9a,oneth)
         T9a23  = T9a13*T9a13
         t9a56  = pow(t9a,fivsix)
         aa     = exp(-0.18d0*T9a*T9a) 
         bb     = 1.06d-03*exp(2.562d0*T9a23)
         cc     = aa + bb
         dd     = 1.72d+31 * T9a56 * (tf% T9i32) * exp(-106.594d0/T9a13)/cc
         else
         !       dd     = 2.6288035d-29
         dd     = 0.0d0
         endif

         ! branching ratios from pwnsz data
         b27n = 0.1d0
         b27p = 0.5d0
         b24a = 0.4d0

         ! c12(o16,n)si27
         term    = dd * b27n
         fr1     = term

         rev    = 0.0d0
         if ((tf% T9) .gt. 0.1d0) then
         rev    = 1.58d0 * exp(4.8972467d0*(tf% T9i))
         end if
         rr1    = rev * term

         ! c12(o16,p)al27
         term    = dd * b27p
         fr2     = term

         rev    = 0.0d0
         if ((tf% T9) .gt. 0.1d0) then
         rev    = 1.58d0 * exp(-59.9970745d0*(tf% T9i))
         end if
         rr2    = rev * term

         ! c12(o16,a)mg24
         term    = dd * b24a
         fr3     = term
         rev    = 0.0d0
         if ((tf% T9) .gt. 0.1d0) then
         rev    = 2.83d0 * exp(-78.5648345d0*(tf% T9i))
         end if
         rr3    = rev * term

         end subroutine rate_c12o16npa

! rc13pg, c13(p,g)n14 

      subroutine rate_c13pg_nacre(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, rev, bb, gs
         ! a0 T9i23 exp(-a1 T9i13 - (T9*a2)^2) 
         !     * (1 + b0 T9 + b1 T92 + b2 T93 + b3 T94 + b4 T95) 
         ! + c0 T9i32 exp(-c1/T9)
         ! + d0 T9i32 exp(-d1/T9)
         ! + e0 T9^e1 exp(-e2/T9)

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         call rnacre(tf,  &
            9.57d7, 13.72d0, 1d0, & ! a0, a1, a2
            3.56d0, 0d0, 0d0, 0d0, 0d0, & ! b0, b1, b2, b3, b4
            1.50d6, 5.930d0, &  ! c0, c1
            0d0, 0d0, & ! d0, d1
            6.83d5, -0.864d0, 12.057d0, & ! e0, e1, e2
            gs)               
         bb   = 2.070d0 * exp(-37.938d0*(tf% T9i))
         if (bb > 1) then ! guard against rate going negative
            bb  = 1
         end if      
         term    = gs * (1 - bb)
         call rnacre_rev(tf, &  ! a0 T932 exp(-a1/T9)
            1.190d10, 87.619d0, &  ! a0, a1
            rev)     
         fr    = term 
         rr    = rev * term 
      end subroutine rate_c13pg_nacre


      subroutine rate_c13pg_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!         p  c13  n14                       nacrr     7.55100d+00          
         call jina_reaclib_2_1(ih1, ic13, in14, tf, fr, rr, 'rate_c13pg_jina')
      end subroutine rate_c13pg_jina

! rc13an, c13(a,n)o16
      subroutine rate_c13an_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!       he4  c13    n  o16                  nacrn     2.21600d+00          
         call jina_reaclib_2_2(ihe4, ic13, ineut, io16, tf, fr, rr, 'rate_c13an_jina')
      end subroutine rate_c13an_jina

      subroutine rate_c13an_fxt(tf, temp, fr, rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
      real(dp) term, rev, aa, bb, cc, dd, ee, ff, gg, q1
      parameter        (q1 = 1.0d0/1.648656d0)

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
      aa  = 6.77d+15 * (tf% T9i23) * exp(-32.329d0*(tf% T9i13) - (tf% T92)*q1)
      bb  = 1.0d0 + 0.013d0*(tf% T913) + 2.04d0*(tf% T923) + 0.184d0*(tf% T9)
      cc   = aa * bb
      dd   = 3.82d+05 * (tf% T9i32) * exp(-9.373d0*(tf% T9i))
      ee   = 1.41d+06 * (tf% T9i32) * exp(-11.873d0*(tf% T9i))
      ff   = 2.0d+09 * (tf% T9i32) * exp(-20.409d0*(tf% T9i))
      gg   = 2.92d+09 * (tf% T9i32) * exp(-29.283d0*(tf% T9i))
      term    = cc + dd + ee + ff + gg
      fr    = term 
      rev    = 5.79d+00 * exp(-25.711d0*(tf% T9i))
      rr    = rev * term
      end subroutine rate_c13an_fxt


! Nitrogen
      

! rn13pg, n13(p,g)o14       

      subroutine rate_n13pg_nacre(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, rev
         ! a0 T9i23 exp(-a1 T9i13 - (T9*a2)^2) 
         !     * (1 + b0 T9 + b1 T92 + b2 T93 + b3 T94 + b4 T95) 
         ! + c0 T9i32 exp(-c1/T9)
         ! + d0 T9i32 exp(-d1/T9)
         ! + e0 T9^e1 exp(-e2/T9)

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         call rnacre(tf,  &
            4.02d7, 15.205d0, 1d0/0.54d0, & ! a0, a1, a2
            3.81d0, 18.6d0, 32.3d0, 0d0, 0d0, & ! b0, b1, b2, b3, b4
            0d0, 0d0, &  ! c0, c1
            0d0, 0d0, & ! d0, d1
            3.25d5, -1.35d0, 5.926d0, & ! e0, e1, e2
            term)        
         call rnacre_rev(tf, &  ! a0 T932 exp(-a1/T9)
            3.571d10, 53.705d0, &  ! a0, a1
            rev)     
         fr    = term 
         rr    = rev * term 
      end subroutine rate_n13pg_nacre

      subroutine rate_n13pg_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!         p  n13  o14                       lg06n     4.62797d+00          
         call jina_reaclib_2_1(ih1, in13, io14, tf, fr, rr, 'rate_n13pg_jina')
      end subroutine rate_n13pg_jina


! rn13ap n13(a,p)o16     cf88
      subroutine rate_n13ap_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!       he4  n13    p  o16                  cf88n     5.21800d+00          
         call jina_reaclib_2_2(ihe4, in13, ih1, io16, tf, fr, rr, 'rate_n13ap_jina')
      end subroutine rate_n13ap_jina

! rn13gp, n13(g,p)c12
   ! see c12pg

! n14(p,g)o15          

      subroutine rate_n14pg_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!         p  n14  o15                       im05n     7.29680d+00          
         call jina_reaclib_2_1(ih1, in14, io15, tf, fr, rr, 'rate_n14pg_jina')
      end subroutine rate_n14pg_jina


      subroutine rate_n14pg_fxt(tf, temp, fr, rr)
! rn14pg ro15gp
! n14(p, g)o15
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
      real(dp) term, rev, aa, bb, cc, dd, ee, q1
      parameter        (q1 = 1.0d0/10.850436d0)

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
      aa  = 4.90d+07 * (tf% T9i23) * exp(-15.228d0*(tf% T9i13) - (tf% T92)*q1)
      bb   = 1.0d0 + 0.027d0*(tf% T913) - 0.778d0*(tf% T923) - 0.149d0*(tf% T9)  &
             + 0.261d0*(tf% T943) + 0.127d0*(tf% T953)
      cc   = aa * bb
      dd   = 2.37d+03 * (tf% T9i32) * exp(-3.011d0*(tf% T9i))
      ee   = 2.19d+04 * exp(-12.530d0*(tf% T9i))
      term    = cc + dd + ee
      rev    = 2.70d+10 * (tf% T932) * exp(-84.678d0*(tf% T9i))
      fr    = term 
      rr    = rev * term 
      end subroutine rate_n14pg_fxt
      

      subroutine rate_n14pg_nacre(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, rev
         ! a0 T9i23 exp(-a1 T9i13 - (T9*a2)^2) 
         !     * (1 + b0 T9 + b1 T92 + b2 T93 + b3 T94 + b4 T95) 
         ! + c0 T9i32 exp(-c1/T9)
         ! + d0 T9i32 exp(-d1/T9)
         ! + e0 T9^e1 exp(-e2/T9)

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         call rnacre(tf,  &
            4.83d7, 15.231d0, 1d0/0.8d0, & ! a0, a1, a2
            -2.00d0, 3.41d0, -2.43d0, 0d0, 0d0, & ! b0, b1, b2, b3, b4
            2.36d3, 3.010d0, &  ! c0, c1
            0d0, 0d0, & ! d0, d1
            6.72d3, 0.380d0, 9.530d0, & ! e0, e1, e2
            term)        
         call rnacre_rev(tf, &  ! a0 T932 exp(-a1/T9)
            2.699d10, 84.677d0, &  ! a0, a1
            rev)     
         fr    = term 
         rr    = rev * term 
      end subroutine rate_n14pg_nacre
      
! rn14ap, n14(a,p)o17
   ! see ro17pa
        
! rn14gp, n14(g,p)c13     
   ! see rc13pg
        
! rn14ag, n14(a,g)f18          
         

      subroutine rate_n14ag_nacre(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, rev, bb, gs

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         if (tf% T9 <= 2d0) then  
            ! a0 T9i23 exp(-a1 T9i13 - (T9*a2)^2) 
            !     * (1 + b0 T9 + b1 T92 + b2 T93 + b3 T94 + b4 T95) 
            ! + c0 T9i32 exp(-c1/T9)
            ! + d0 T9i32 exp(-d1/T9)
            ! + e0 T9^e1 exp(-e2/T9)
            call rnacre(tf,  &
               7.93d11, 36.035d0, 1d0/0.07d0, & ! a0, a1, a2
               0d0, 0d0, 0d0, 0d0, 0d0, & ! b0, b1, b2, b3, b4
               1.85d-10, 2.750d0, &  ! c0, c1
               2.62d0, 5.045d0, & ! d0, d1
               2.93d3, 0.344d0, 10.561d0, & ! e0, e1, e2
               gs)                  
         else   
            gs   = 1.52d2 * pow((tf% T9),1.567d0) * exp(-6.315d0*(tf% T9i)) 
         end if   
         bb   = 0.340d0 * exp(-26.885d0*(tf% T9i) - 0.012d0*(tf% T9))
         if (bb > 1) then ! guard against rate going negative
            bb  = 1
         end if      
         term    = gs * (1 - bb)
         call rnacre_rev(tf, &  ! a0 T932 exp(-a1/T9)
            5.420d10, 51.231d0, &  ! a0, a1
            rev)     
         fr    = term 
         rr    = rev * term 
      end subroutine rate_n14ag_nacre
 
! n14(a,g)f18      
      subroutine rate_n14ag_jina(tf, temp, fr, rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
!       he4  n14  f18                       ga00r     4.41500d+00          
         call jina_reaclib_2_1(ihe4, in14, if18, tf, fr, rr, 'rate_n14ag_jina')
      end subroutine rate_n14ag_jina


! rn15pg, n15(p,g)o16           


      subroutine rate_n15pg_nacre(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, rev
         if (tf% T9 <= 3.5d0) then      
            ! a0 T9i23 exp(-a1 T9i13 - (T9*a2)^2) 
            !     * (1 + b0 T9 + b1 T92 + b2 T93 + b3 T94 + b4 T95) 
            ! + c0 T9i32 exp(-c1/T9)
            ! + d0 T9i32 exp(-d1/T9)
            ! + e0 T9^e1 exp(-e2/T9)
            call rnacre(tf,  &
               1.08d9, 15.254d0, 1d0/0.34d0, & ! a0, a1, a2
               6.15d0, 16.4d0, 0d0, 0d0, 0d0, & ! b0, b1, b2, b3, b4
               9.23d3, 3.597d0, &  ! c0, c1
               3.27d6, 11.024d0, & ! d0, d1
               0d0, 0d0, 0d0, & ! e0, e1, e2
            term)        
         else
            term     = 3.54d4 * pow((tf% T9),0.095d0) * exp(-2.306d0*(tf% T9i)) 
         end if   
         call rnacre_rev(tf, &  ! a0 T932 exp(-a1/T9)
            3.622d10, 140.73d0, &  ! a0, a1
            rev)     
         fr    = term 
         rr    = rev * term 
      end subroutine rate_n15pg_nacre


      subroutine rate_n15pg_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         call jina_reaclib_2_1(ih1, in15, io16, tf, fr, rr, 'rate_n15pg_jina')
      end subroutine rate_n15pg_jina


! rn15pa, n15(p,a)c12
   ! see rc12ap
   
! rn15ap, n15(a,p)o18 
   ! see ro18pa
   

! Oxygen


! ro14ap, o14(a,p)f17          

      subroutine rate_o14ap_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, rev, aa, bb, cc, dd, ee, ff, q1
         parameter        (q1 = 1.0d0/0.514089d0)

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         aa  = 1.68d+13 * (tf% T9i23) * exp(-39.388d0*(tf% T9i13)- (tf% T92)*q1)
         bb  = 1.0d0 + 0.011d0*(tf% T913) + 13.117d0*(tf% T923) + 0.971d0*(tf% T9)  &
            + 85.295d0*(tf% T943) + 16.061d0*(tf% T953)
         cc  = aa * bb
         dd  = 3.31d+04 * (tf% T9i32) * exp(-11.733d0*(tf% T9i))
         ee  = 1.79d+07 * (tf% T9i32) * exp(-22.609d0*(tf% T9i)) 
         ff  = 9.00d+03 * (tf% T9113) * exp(-12.517d0*(tf% T9i))
         term    = cc + dd + ee + ff
         fr    = term 
         rev      = 4.93d-01*exp(-13.820d0*(tf% T9i))
         rr    = rev * term 
      end subroutine rate_o14ap_fxt


      subroutine rate_o14ap_jina(tf, temp, fr, rr) !  Hahn 1996    PhRvC 54, 4, p1999-2013
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!       he4  o14    p  f17                  Ha96r     1.19200d+00          
         call jina_reaclib_2_2(ihe4, io14, ih1, if17, tf, fr, rr, 'rate_o14ap_jina')
      end subroutine rate_o14ap_jina
      
      
! ro14ag, o14(a,g)ne18  
      subroutine rate_o14ag_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!       he4  o14 ne18                       wh87n     5.11400d+00          
         call jina_reaclib_2_1(ihe4, io14, ine18, tf, fr, rr, 'rate_o14ag_jina')
      end subroutine rate_o14ag_jina




! ro14gp, o14(g,p)n13
   ! see rn13pg


! ro15ap, o15(a,p)f18  
   ! see rf18pa
              
! ro15ag, o15(a,g)ne19          

      subroutine rate_o15ag_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, rev, aa, bb, cc, dd, ee, ff, gg, hh, q1, q2, q3
         parameter        (q1 = 1.0d0/9.0d0,  &
                           q2 = 1.0d0/3.751969d0,  &
                           q3 = 1.0d0/64.0d0)

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         aa  = 3.57d+11 * (tf% T9i23) * exp(-39.584d+0*(tf% T9i13) - (tf% T92)*q1)
         bb  = 1.0d0 + 0.011d0*(tf% T913) - 0.273d0*(tf% T923) - 0.020d0*(tf% T9)
         cc  = aa*bb
         dd  = 5.10d+10 * (tf% T9i23) * exp(-39.584d+0*(tf% T9i13) - (tf% T92)*q2)
         ee  = 1.0d0 + 0.011d0*(tf% T913) + 1.59d0*(tf% T923) + 0.117d0*(tf% T9) &
            + 1.81d0*(tf% T943) + 0.338d0*(tf% T953)
         ff  = dd*ee
         gg  = 3.95d-1 * (tf% T9i32) * exp(-5.849d0*(tf% T9i))
         hh  = 1.90d+1 * pow((tf% T9),2.85d0) * exp(-7.356d0*(tf% T9i) - (tf% T92)*q3)
         term    = cc + ff + gg + hh
         fr    = term
         rev      = 5.54d+10 * (tf% T932) * exp(-40.957d0*(tf% T9i))
         rr    = rev * term
      end subroutine rate_o15ag_fxt


      subroutine rate_o15ag_jina(tf, temp, fr, rr) 
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!       he4  o15 ne19                       Ha96n     3.52900d+00          
         call jina_reaclib_2_1(ihe4, io15, ine19, tf, fr, rr, 'rate_o15ag_jina')
      end subroutine rate_o15ag_jina


! ro15gp, o15(g,p)n14
   ! see rn14pg                     

! ro16pg, o16(p,g)f17          
      subroutine rate_o16pg_nacre(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, rev, aa, bbm1, bb

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         aa   = 7.37d7 * pow((tf% T9),-0.82d0) * exp(-16.696d0*(tf% T9i13))          
         bbm1 = 202d0 * exp(-70.348d0*(tf% T9i) - 0.161d0*(tf% T9))
         bb   = 1 + bbm1         
         term    = aa * bb         
         call rnacre_rev(tf, &  ! a0 T932 exp(-a1/T9)
            3.037d9, 6.966d0, &  ! a0, a1
            rev)     
         fr    = term 
         rr    = rev * term 
      end subroutine rate_o16pg_nacre
         

      subroutine rate_o16pg_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!         p  o16  f17                       nacrn     6.00000d-01          
         call jina_reaclib_2_1(ih1, io16, if17, tf, fr, rr, 'rate_o16pg_jina')
      end subroutine rate_o16pg_jina
     
! ro16ap, o16(a,p)f19
   ! see rf19pa
                        
! ro16ag, o16(a,g)ne20                  

      subroutine rate_o16ag_nacre(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, rev

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         ! a0 T9i23 exp(-a1 T9i13 - (T9*a2)^2) 
         !     * (1 + b0 T9 + b1 T92 + b2 T93 + b3 T94 + b4 T95) 
         ! + c0 T9i32 exp(-c1/T9)
         ! + d0 T9i32 exp(-d1/T9)
         ! + e0 T9^e1 exp(-e2/T9)
         call rnacre(tf,  &
            2.68d10, 39.760d0, 1d0/1.6d0, & ! a0, a1, a2
            0d0, 0d0, 0d0, 0d0, 0d0, & ! b0, b1, b2, b3, b4
            51.1d0, 10.32d0, &  ! c0, c1
            616.1d0, 12.200d0, & ! d0, d1
            0.41d0, 2.966d0, 11.900d0, & ! e0, e1, e2
            term)        
         call rnacre_rev(tf, &  ! a0 T932 exp(-a1/T9)
            5.653d10, 54.886d0, &  ! a0, a1
            rev)     
         fr    = term
         rr    = rev * term
      end subroutine rate_o16ag_nacre

      subroutine rate_o16ag_jina(tf, temp, fr, rr) ! jina reaclib -- nacre
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!       he4  o16 ne20                       nacrr     4.73000d+00          
         call jina_reaclib_2_1(ihe4, io16, ine20, tf, fr, rr, 'rate_o16ag_jina')
      end subroutine rate_o16ag_jina



! ro16gp, o16(g,p)n15
   ! see rn15pg
               
! ro16ga, o16(g,a)c12   
   ! see rc12ag
                        
! r1616 cf88 fxt

      subroutine rate_o16o16_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         term  = 7.10d36 * (tf% T9i23) * &
              exp(-135.93d0*(tf% T9i13) - 0.629d0*(tf% T923)  &
                   - 0.445d0*(tf% T943) + 0.0103d0*(tf% T9)*(tf% T9))
         fr    = term
         rr    = 0.0d0      
      end subroutine rate_o16o16_fxt

      subroutine rate_o16o16g_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!
         call rate_o16o16_fxt(tf, temp, fr, rr)
         fr    = fr*0.10d0
      end subroutine rate_o16o16g_fxt

! r1616n, o16(o16,n)s31


! r1616p, o16(o16, p)p31
      subroutine rate_o16o16p_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!       o16  o16    p  p31                  cf88r     7.67800d+00          
         call jina_reaclib_2_2(io16, io16, ih1, ip31, tf, fr, rr, 'rate_o16o16p_jina')
      end subroutine rate_o16o16p_jina


! r1616a, o16(o16, a)si28
      subroutine rate_o16o16a_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!       o16  o16  he4 si28                  cf88r     9.59300d+00          
         call jina_reaclib_2_2(io16, io16, ihe4, isi28, tf, fr, rr, 'rate_o16o16a_jina')
      end subroutine rate_o16o16a_jina

      subroutine rate_o16o16a(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) :: fr1, rr1, fr2, rr2, fr4, rr4
         call rate_o16o16npad(tf, temp, fr1, rr1, fr2, rr2, fr, rr, fr4, rr4)
      end subroutine rate_o16o16a

      subroutine rate_o16o16a_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!
         call rate_o16o16_fxt(tf, temp, fr, rr)
         fr    = fr*0.56d0
      end subroutine rate_o16o16a_fxt

      subroutine rate_o16o16p_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!
         call rate_o16o16_fxt(tf, temp, fr, rr)
         fr    = fr*0.34d0
      end subroutine rate_o16o16p_fxt

      subroutine rate_o16o16_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) :: fr1, rr1, fr2, rr2
         call rate_o16o16a_jina(tf, temp, fr1, rr1)
         call rate_o16o16p_jina(tf, temp, fr2, rr2)
         fr = fr1 + fr2
         rr = rr1 + rr2
      end subroutine rate_o16o16_jina

      subroutine rate_o16o16npad(tf, temp, &
         fr1, rr1, &       ! o16(o16, n)s31 &
         fr2, rr2, &       ! o16(o16, p)p31 &
         fr3, rr3, &       ! o16(o16, a)si28 &
         fr4, rr4)       ! o16(o16, d)p30

         type (T_Factors), pointer :: tf
         real(dp) temp, fr1, rr1, fr2, rr2, fr3, rr3, fr4, rr4

         real(dp) term, rev, aa, daa,  &
                  b32n, b32p, b32a, b32d, ezro, dlt, xxt, thrs

         if (tf% t9 < lowT9_cutoff) then
         fr1 = 0; rr1 = 0
         fr2 = 0; rr2 = 0
         fr3 = 0; rr3 = 0
         fr4 = 0; rr4 = 0
         return
         end if 


         aa  = 7.10d36 * (tf% T9i23) * &
         exp(-135.93d0*(tf% T9i13) - 0.629d0*(tf% T923)  &
               - 0.445d0*(tf% T943) + 0.0103d0*(tf% T9)*(tf% T9))

         daa = -twoth*aa*(tf% T9i) &
         + aa * (oneth*135.93d0*(tf% T9i43) - twoth*0.629d0*(tf% T9i13) &
                     - fourth*0.445d0*(tf% T913) + 0.0206d0*(tf% T9))


         ! branching ratios highly uncertain;  guessed using fcz 1975
         ! deuteron channel is endoergic. apply error function cut-off.
         ezro = 3.9d0*(tf% T923)
         dlt  = 1.34d0*pow(tf% T9,fivsix)
         xxt  = 2.0d0*(2.406d0 - ezro)/dlt
         call fowthrsh(xxt, thrs)
         b32d  = 0.05d0*thrs
         b32n  = 0.1d0
         b32a  = 0.25d0
         b32p  = 1.0d0 - b32d - b32a - b32n

         ! o16(o16, n)s31
         term    = aa * b32n
         fr1     = term
         rev    = 0.0d0
         if ((tf% T9) .gt. 0.1d0) then
         rev    = 5.92d0 * exp(-16.8038228d0*(tf% T9i))
         end if
         rr1    = rev * term

         ! o16(o16, p)p31
         term    = aa * b32p
         fr2     = term
         rev    = 0.0d0
         if ((tf% T9) .gt. 0.1d0) then
         rev    = 5.92d0*exp(-89.0788286d0*(tf% T9i))
         end if
         rr2    = rev * term

         ! o16(o16, a)si28
         term    = aa * b32a
         fr3     = term
         rev    = 0.0d0
         if ((tf% T9) .gt. 0.1d0) then
         rev    = 3.46d0*exp(-111.3137212d0*(tf% T9i))
         end if
         rr3    = rev * term

         ! o16(o16, d)p30
         term    = aa * b32d
         fr4     = term
         rev    = 0.0d0
         if ((tf% T9) .gt. 0.1d0) then
         rev    = 0.984d0*exp(27.9908982d0*(tf% T9i))
         end if
         rr4    = rev * term
         end subroutine rate_o16o16npad


      subroutine fowthrsh(x, thrs)

! fowler threshold fudge function. 
! err func rational (abramowitz p.299)7.1.25 and its derivative

! declare
      real(dp) x, thrs, ag, z, z2, t, t2, t3, tt, er, aa
      ag   = sign(1.0d0, x)
      z    = abs(x)
      z2   = z*z
      aa   = 1.0d0 + 0.47047d0*z
      t    = 1.0d0/aa
      t2   = t*t
      t3   = t2*t
      tt   = 0.3480242d0*t - 0.0958798d0*t2 + 0.7478556d0*t3
      thrs  = 0.5d0
      if (z .ne. 0) then
       aa   = exp(-z2)
       er   = 1.0d0 - tt * aa
       thrs  = 0.5d0 * (1.0d0 - ag*er)
      end if
      end subroutine fowthrsh

      
      
! o17(a,g)ne21
      
      
! ro17pa, o17(p,a)n14                    

      subroutine rate_o17pa_nacre(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, rev, bb, cc, gs 
         include 'formats.dek'

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         if (tf% T9 <= 6d0) then      
            ! a0 T9i23 exp(-a1 T9i13 - (T9*a2)^2) 
            !     * (1 + b0 T9 + b1 T92 + b2 T93 + b3 T94 + b4 T95) 
            ! + c0 T9i32 exp(-c1/T9)
            ! + d0 T9i32 exp(-d1/T9)
            ! + e0 T9^e1 exp(-e2/T9)
            call rnacre(tf,  &
            9.20d8, 16.715d0, 1d0/0.06d0, & ! a0, a1, a2
            -80.31d0, 2211d0, 0d0, 0d0, 0d0, & ! b0, b1, b2, b3, b4
            9.13d-4, 0.7667d0, &  ! c0, c1
            9.68d0, 2.083d0, & ! d0, d1
            1.85d6, 1.591d0, 4.848d0, & ! e0, e1, e2
            gs)              
            cc   = 8.13d6 * (tf% T9i32) * exp(-5.685d0*(tf% T9i)) 
            gs = gs + cc
         else
            gs   = 8.73d6 * pow(tf% T9,0.950d0) * exp(-7.508d0*(tf% T9i)) 
         end if   
         bb   = 1.033d0 * exp(-10.034d0*(tf% T9i) - 0.165d0*(tf% T9))
         term    = gs * (1 + bb)
         call rnacre_rev(tf, &  ! a0 T932 exp(-a1/T9)
            6.759d-1, 13.829d0, &  ! a0, a1
            rev)     
         fr    = term 
         rr    = rev * term 
         
         
         return
         
         
         call show_nacre_terms(tf,  &
            9.20d8, 16.715d0, 1d0/0.06d0, & ! a0, a1, a2
            -80.31d0, 2211d0, 0d0, 0d0, 0d0, & ! b0, b1, b2, b3, b4
            9.13d-4, 0.7667d0, &  ! c0, c1
            9.68d0, 2.083d0, & ! d0, d1
            1.85d6, 1.591d0, 4.848d0)
         write(*,1) 'gs', gs
         write(*,1) 'bb', bb
         write(*,1) 'term', term
         write(*,*) 'rate_o17pa_nacre' 
         
      end subroutine rate_o17pa_nacre
      
      
      subroutine rate_o17pa_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!         p  o17  he4  n14                  ct07r     1.19164d+00          
         call jina_reaclib_2_2(ih1, io17, ihe4, in14, tf, fr, rr, 'rate_o17pa_jina')
      end subroutine rate_o17pa_jina


! ro17pg, o17(p,g)f18                           

      subroutine rate_o17pg_nacre(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, rev, bb, gs

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         if (tf% T9 <= 3d0) then  
            ! a0 T9i23 exp(-a1 T9i13 - (T9*a2)^2) 
            !     * (1 + b0 T9 + b1 T92 + b2 T93 + b3 T94 + b4 T95) 
            ! + c0 T9i32 exp(-c1/T9)
            ! + d0 T9i32 exp(-d1/T9)
            ! + e0 T9^e1 exp(-e2/T9)
            call rnacre(tf,  &
            1.50d8, 16.710d0, 1d0/0.2d0, & ! a0, a1, a2
            0d0, 0d0, 0d0, 0d0, 0d0, & ! b0, b1, b2, b3, b4
            9.79d-6, 0.7659d0, &  ! c0, c1
            4.15d0, 2.083d0, & ! d0, d1
            7.74d4, 1.16d0, 6.342d0, & ! e0, e1, e2
            gs)                     
         else   
            gs   = 1.74d3 * pow(tf% T9,0.700d0) * exp(-1.072d0*(tf% T9i)) 
          end if   
         bb   = 0.287d0 * exp(-10.011d0*(tf% T9i) - 0.062d0*(tf% T9))
         if (bb > 1) then ! guard against rate going negative
            bb  = 1
         end if      
         term    = gs * (1 - bb)
         call rnacre_rev(tf, &  ! a0 T932 exp(-a1/T9)
            3.663d10, 65.060d0, &  ! a0, a1
            rev)     
         fr    = term 
         rr    = rev * term 
      end subroutine rate_o17pg_nacre
      
      
      subroutine rate_o17pg_jina(tf, temp, fr, rr) ! jina reaclib   Chafa et al. (2007)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!         p  o17  f18                       ct07n     5.60650d+00          
         call jina_reaclib_2_1(ih1, io17, if18, tf, fr, rr, 'rate_o17pg_jina')
      end subroutine rate_o17pg_jina

! ro18pa, o18(p,a)n15          

      subroutine rate_o18pa_nacre(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, rev, bb, gs

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         ! a0 T9i23 exp(-a1 T9i13 - (T9*a2)^2) 
         !     * (1 + b0 T9 + b1 T92 + b2 T93 + b3 T94 + b4 T95) 
         ! + c0 T9i32 exp(-c1/T9)
         ! + d0 T9i32 exp(-d1/T9)
         ! + e0 T9^e1 exp(-e2/T9)
         call rnacre(tf,  &
            5.58d11, 16.732d0, 1d0/0.51d0, & ! a0, a1, a2
            3.2d0, 21.8d0, 0d0, 0d0, 0d0, & ! b0, b1, b2, b3, b4
            9.91d-14, 0.232d0, &  ! c0, c1
            2.58d4, 1.665d0, & ! d0, d1
            3.24d8, -0.378d0, 6.395d0, & ! e0, e1, e2
            gs)               
         bb   = 1.968d0 * exp(-25.673d0*(tf% T9i) - 0.083d0*(tf% T9))
         if (bb > 1) then ! guard against rate going negative
            bb  = 1
         end if      
         term    = gs * (1 - bb)
         call rnacre_rev(tf, &  ! a0 T932 exp(-a1/T9)
            1.660d-1, 46.192d0, &  ! a0, a1
            rev)     
         fr    = term 
         rr    = rev * term 
      end subroutine rate_o18pa_nacre


      subroutine rate_o18pa_jina(tf, temp, fr, rr) ! jina reaclib    nacre
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!         p  o18  he4  n15                  nacrn     3.98100d+00          
         call jina_reaclib_2_2(ih1, io18, ihe4, in15, tf, fr, rr, 'rate_o18pa_jina')
      end subroutine rate_o18pa_jina
      

! ro18pg, o18(p,g)f19                     

      subroutine rate_o18pg_nacre(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, rev, bb, dd, gs

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         if (tf% T9 <= 2d0) then   
            ! a0 T9i23 exp(-a1 T9i13 - (T9*a2)^2) 
            !     * (1 + b0 T9 + b1 T92 + b2 T93 + b3 T94 + b4 T95) 
            ! + c0 T9i32 exp(-c1/T9)
            ! + d0 T9i32 exp(-d1/T9)
            ! + e0 T9^e1 exp(-e2/T9)
            call rnacre(tf,  &
               4.59d8, 16.732d0, 1d0/0.15d0, & ! a0, a1, a2
               -9.02d0, 506d0, -2400d0, 0d0, 0d0, & ! b0, b1, b2, b3, b4
               9.91d-17, 0.232d0, &  ! c0, c1
               3.30d-3, 1.033d0, & ! d0, d1
               1.25d4, 0.458d0, 5.297d0, & ! e0, e1, e2
               gs)        
            dd   = 1.61d2 * (tf% T9i32) * exp(-1.665d0*(tf% T9i)) 
            gs   = gs + dd
         else   
            gs   = 1.38d4 * pow(tf% T9,0.829d0) * exp(-5.919d0*(tf% T9i)) 
          end if   
         bb   = 0.475d0 * exp(-15.513d0*(tf% T9i) - 0.102d0*(tf% T9))
         if (bb > 1) then ! guard against rate going negative
            bb  = 1
         end if      
         term    = gs * (1 - bb)
         call rnacre_rev(tf, &  ! a0 T932 exp(-a1/T9)
            9.201d9, 92.769d0, &  ! a0, a1
            rev)     
         fr    = term 
         rr    = rev * term 
      end subroutine rate_o18pg_nacre


      subroutine rate_o18pg_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!         p  o18  f19                       nacrr     7.99400d+00          
         call jina_reaclib_2_1(ih1, io18, if19, tf, fr, rr, 'rate_o18pg_jina')
      end subroutine rate_o18pg_jina


! ro18ag, o18(a,g)ne22                       

      subroutine rate_o18ag_nacre(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) :: gs, term, term1, bb, cc, dd, rev    

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         if (tf% T9 <= 6d0) then   
            ! a0 T9i23 exp(-a1 T9i13 - (T9*a2)^2) 
            !     * (1 + b0 T9 + b1 T92 + b2 T93 + b3 T94 + b4 T95) 
            ! + c0 T9i32 exp(-c1/T9)
            ! + d0 T9i32 exp(-d1/T9)
            ! + e0 T9^e1 exp(-e2/T9)
            call rnacre(tf,  &
               0d0, 0d0, 0d0, & ! a0, a1, a2
               0d0, 0d0, 0d0, 0d0, 0d0, & ! b0, b1, b2, b3, b4
               1.95d-13, 2.0690d0, &  ! c0, c1
               1.56d-2, 4.462d0, & ! d0, d1
               3.44d5, -0.5d0, 22.103d0, & ! e0, e1, e2
            gs)             
            cc   = 1.01d0 * (tf% T9i32) * exp(-6.391d0*(tf% T9i)) 
            dd   = 44.1d0 * (tf% T9i32) * exp(-7.389d0*(tf% T9i)) 
            gs    = gs + cc + dd
         else   
            gs   = 3.31d5 * pow(tf% T9,-0.221d0) * exp(-24.990d0*(tf% T9i)) 
         end if   
         bb   = 1.411d0 * exp(-20.533d0*(tf% T9i) - 0.0382d0*(tf% T9))
         if (bb > 1) then ! guard against rate going negative
            bb  = 1
         end if      
         term    = gs * (1 - bb)
         call rnacre_rev(tf, &  ! a0 T932 exp(-a1/T9)
            5.847d10, 112.18d0, &  ! a0, a1
            rev)     
         fr    = term 
         rr    = rev * term 
      end subroutine rate_o18ag_nacre


      subroutine rate_o18ag_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!       he4  o18 ne22                       dh03r     9.66900d+00          
         call jina_reaclib_2_1(ihe4, io18, ine22, tf, fr, rr, 'rate_o18ag_jina')
      end subroutine rate_o18ag_jina

       
! Fluorine 


! rf17pa, f17(p,a)o14    
   ! see ro14ap
         


! rf17gp, f17(g,p)o16    
   ! see ro16pg 
   
   
! rf17ap    f17(a,p)ne20     
      subroutine rate_f17ap_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!       he4  f17    p ne20                  nacr      4.13000d+00          
         call jina_reaclib_2_2(ihe4, if17, ih1, ine20, tf, fr, rr, 'rate_f17ap_jina')
      end subroutine rate_f17ap_jina

! rf18pa, f18(p,a)o15              

      subroutine rate_f18pa_wk82(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, rev, aa, bb, cc, dd, ee, ff

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         aa  = 1.66d-10 * (tf% T9i32) * exp(-0.302d0*(tf% T9i))
         bb  = 1.56d+05 * (tf% T9i32) * exp(-3.84d0*(tf% T9i))
         cc  = 1.36d+06 * (tf% T9i32) * exp(-5.22d0*(tf% T9i))
         dd  = 8.1d-05 * (tf% T9i32) * exp(-1.05d0*(tf% T9i))
         ee  = 8.9d-04 * (tf% T9i32) * exp(-1.51d0*(tf% T9i))
         ff  = 3.0d+05 * (tf% T9i32) * exp(-4.29d0*(tf% T9i))
         term    = aa + bb + cc + dd + ee + ff
         fr    = term 
         rev      = 4.93d-01 * exp(-33.433d0*(tf% T9i))
         rr    = rev * term 
      end subroutine rate_f18pa_wk82


      subroutine rate_f18pa_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!         p  f18  he4  o15                  sh03r     2.88215d+00          
         call jina_reaclib_2_2(ih1, if18, ihe4, io15, tf, fr, rr, 'rate_f18pa_jina')
      end subroutine rate_f18pa_jina


! rf18gp, f18(g,p)o17
   ! see ro17pg                  

! rf19pg, f19(p,g)ne20                      

      subroutine rate_f19pg_nacre(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, rev, bb, gs

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         if (tf% T9 <= 1.5d0) then
            ! a0 T9i23 exp(-a1 T9i13 - (T9*a2)^2) 
            !     * (1 + b0 T9 + b1 T92 + b2 T93 + b3 T94 + b4 T95) 
            ! + c0 T9i32 exp(-c1/T9)
            ! + d0 T9i32 exp(-d1/T9)
            ! + e0 T9^e1 exp(-e2/T9)
            call rnacre(tf,  &
               6.37d7, 18.116d0, 0d0, & ! a0, a1, a2
               0.775d0, 36.1d0, 0d0, 0d0, 0d0, & ! b0, b1, b2, b3, b4
               8.27d2, 3.752d0, &  ! c0, c1
               0d0, 0d0, & ! d0, d1
               1.28d6, -3.667d0, 9.120d0, & ! e0, e1, e2
            gs)        
         else   
            gs   = 3.66d3 * pow(tf% T9,0.947d0) * exp(-2.245d0*(tf% T9i)) 
         end if   
         bb   = 0.990d0 * exp(-1.207d0*(tf% T9i) - 0.0886d0*(tf% T9))
         if (bb > 1) then ! guard against rate going negative
            bb  = 1
         end if      
         term    = gs * (1 - bb)
         call rnacre_rev(tf, &  ! a0 T932 exp(-a1/T9)
            3.696d10, 149.04d0, &  ! a0, a1
            rev)      
         fr    = term 
         rr    = rev * term
      end subroutine rate_f19pg_nacre


      subroutine rate_f19pg_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!         p  f19 ne20                       cf88r     1.28480d+01          
         call jina_reaclib_2_1(ih1, if19, ine20, tf, fr, rr, 'rate_f19pg_jina')
      end subroutine rate_f19pg_jina


! rf19pa, f19(p,a)o16                         

      subroutine rate_f19pa_nacre(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, rev, bb, dd, gs

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         ! a0 T9i23 exp(-a1 T9i13 - (T9*a2)^2) 
         !     * (1 + b0 T9 + b1 T92 + b2 T93 + b3 T94 + b4 T95) 
         ! + c0 T9i32 exp(-c1/T9)
         ! + d0 T9i32 exp(-d1/T9)
         ! + e0 T9^e1 exp(-e2/T9)
         call rnacre(tf,  &
            2.62d11, 18.116d0, 1d0/0.185d0, & ! a0, a1, a2
            6.26d-2, 0.285d0, 4.94d-3, 11.5d0, 7.40d4, & ! b0, b1, b2, b3, b4
            3.80d6, 3.752d0, &  ! c0, c1
            0d0, 0d0, & ! d0, d1
            3.27d7, -0.193d0, 6.587d0, & ! e0, e1, e2
            gs)              
         dd   = 7.30d8 * pow(tf% T9,-0.201d0) * exp(-16.249d0*(tf% T9i)) 
         gs    = gs + dd
         bb   = 0.755d0 * exp(-1.755d0*(tf% T9i) - 0.174d0*(tf% T9))
         term    = gs * (1 + bb)
         call rnacre_rev(tf, &  ! a0 T932 exp(-a1/T9)
            6.538d-1, 94.154d0, &  ! a0, a1
            rev)     
         fr    = term 
         rr    = rev * term 
      end subroutine rate_f19pa_nacre


      subroutine rate_f19pa_jina(tf, temp, fr, rr) ! jina reaclib    
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         call jina_reaclib_2_2(ih1, if19, ihe4, io16, tf, fr, rr, 'rate_f19pa_jina')
      end subroutine rate_f19pa_jina


! rf19gp, f19(g,p)o18 
   ! see ro18pg
   
   
! rf19ap, f19(a,p)ne22 
      subroutine rate_f19ap_cf88(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) :: term
         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
         term = 4.50d18*tf% T9i23*exp(-43.467d0*tf% T9i13-pow(tf% T9/0.637d0,2))+ &
                7.98d04*tf% T932*exp(-12.760d0*tf% T9i)
         fr = term*6.36d00*exp(-19.439d0*tf% T9i)
         rr    = 0.0d0
      end subroutine rate_f19ap_cf88

      subroutine rate_f19ap_jina(tf, temp, fr, rr) ! jina reaclib    
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!       he4  f19    p ne22                  cf88r     1.67500d+00          
         call jina_reaclib_2_2(ihe4, if19, ih1, ine22, tf, fr, rr, 'rate_f19ap_jina')
      end subroutine rate_f19ap_jina
            

! Neon

      
! rne18ap, ne18(a,p)na21

      subroutine rate_ne18ap_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, rev, aa, bb, cc, dd, ee, zz
         real(dp) z1, a1, ztot, ared, r, c1, c2, c3, c4
         parameter     (z1   = 10.0d0,  &
                        a1   = 18.0d0,  &
                        ztot = 2.0d0 * z1,  &
                        ared = 4.0d0*a1/(4.0d0 + a1),  &
                        r    = 5.1566081196876965d0,  &
                        c1   = 4.9080044545315392d10,  &
                        c2   = 4.9592784569936502d-2,  &
                        c3   = 1.9288564401521285d1,  &
                        c4   = 4.6477847042196437d1)

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         ! note:
         !      r    = 1.09 * a1**oneth + 2.3
         !      c1   = 7.833e9 * 0.31 * ztot**fourth/(ared**fivsix)
         !      c2   = 0.08617 * 0.1215 * sqrt(ared*r**3/ztot)
         !      c3   = 2.0d0 * 0.52495 * sqrt(ared*r*ztot)
         !      c4   = 4.2487 * (ztot**2*ared)**oneth
         ! ne18ap(a, p)na21
         ! was a call to aprate
         aa  = 1.0d0 + c2*(tf% T9)
         zz  = c2/aa
         bb  = pow(aa,fivsix)
         cc  = (tf% T923) * bb
         dd = pow(aa,oneth)
         ee  = (tf% T9i13) * dd
         term    = c1*exp(c3 - c4*ee)/cc 
         fr    = term 
         rev    = 0.0d0
         rr    = 0.0d0
      end subroutine rate_ne18ap_fxt

      subroutine rate_ne18ap_jina(tf, temp, fr, rr) ! jina reaclib    
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!       he4 ne18    p na21                  GW95r     2.62700d+00          
         call jina_reaclib_2_2(ihe4, ine18, ih1, ina21, tf, fr, rr, 'rate_ne18ap_jina')
      end subroutine rate_ne18ap_jina


! rne18ag   ne18(a,g)mg22  rath
      subroutine rate_ne18ag_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!       he4 ne18 mg22                       rath      8.14100d+00          
         call jina_reaclib_2_1(ihe4, ine18, img22, tf, fr, rr, 'rate_ne18ag_jina')
      end subroutine rate_ne18ag_jina
      


! rne18gp, ne18(g,p)f17
   ! see rf17pg   


! rne19pg, ne19(p,g)na20

      subroutine rate_ne19pg_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, rev, aa, bb, cc, dd, ee, ff, gg, q1
         parameter        (q1 = 1.0d0/1.304164d0)

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         aa  = 1.71d+6 * (tf% T9i23) * exp(-19.431d0*(tf% T9i13))
         bb  = 1.0d0 + 0.021d0*(tf% T913) + 0.130d0*(tf% T923) + 1.95d-2*(tf% T9) &
            + 3.86d-2*(tf% T943) + 1.47d-02*(tf% T953) 
         cc  = aa*bb
         dd  = 1.89d+5 * (tf% T9i23) * exp(-19.431d0*(tf% T9i13) - (tf% T92)*q1)
         ee  = 1.0d0 + 0.021d0*(tf% T913) + 2.13d0*(tf% T923) + 0.320d0*(tf% T9)  &
            + 2.80d0*(tf% T943) + 1.07d0*(tf% T953)
         ff  = dd*ee
         gg  = 8.45d+3 * (tf% T9i54) * exp(-7.64d0*(tf% T9i))
         term    = cc + ff + gg
         fr    = term 
         rev      = 7.39d+09 * (tf% T932) * exp(-25.519d0*(tf% T9i))
         rr    = rev * term 
      end subroutine rate_ne19pg_fxt
      
      
      subroutine rate_ne19pg_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!         p ne19 na20                       cf88r     2.19900d+00          
         call jina_reaclib_2_1(ih1, ine19, ina20, tf, fr, rr, 'rate_ne19pg_jina')
      end subroutine rate_ne19pg_jina


! rne19ga, ne19(g,a)o15
   ! see r016ag
            
! rne19gp, ne19(g,p)f18
   ! see rf18pg          

! rne20pg, ne20(p,g)na21 

      subroutine rate_ne20pg_nacre(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, rev, aa, bb, cc, dd, ff, gg, zz, gs

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         aa   = 2.35d7 * pow(tf% T9,-1.84d0) * exp(-19.451d0*(tf% T9i13)) * (1 + 10.80d0*(tf% T9))
         gs   = aa
         aa   = 18.0d0 * (tf% T9i32) * exp(-4.247d0*(tf% T9i)) 
         gs   = gs + aa
         aa   = 9.83d0 * (tf% T9i32) * exp(-4.619d0*(tf% T9i)) 
         gs   = gs + aa
         aa   = 6.76d4 * pow(tf% T9,-0.641d0) * exp(-11.922d0*(tf% T9i)) 
         gs   = gs + aa
         bb   = 7.929d0 * exp(-20.108d0*(tf% T9i) - 0.327d0*(tf% T9))
         if (bb > 1) then ! guard against rate going negative
            bb  = 1
         end if         
         term    = gs * (1 - bb)
         rev      = 4.637d9 * (tf% T932) * exp(-28.214d0*(tf% T9i))
         fr    = term 
         rr    = rev * term
      end subroutine rate_ne20pg_nacre
      
      
      subroutine rate_ne20pg_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!         p ne20 na21                       nacrr     2.43100d+00          
         call jina_reaclib_2_1(ih1, ine20, ina21, tf, fr, rr, 'rate_ne20pg_jina')
      end subroutine rate_ne20pg_jina
      
      
! ne20(a,p)na23      
      subroutine rate_ne20ap_jina(tf, temp, fr, rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
      real(dp) :: fr1, rr1
      include 'formats.dek'
!       he4 ne20    p na23                  ha04rv   -2.37900d+00          
      call jina_reaclib_2_2(ih1, ina23, ihe4, ine20, tf, rr, fr, 'rate_ne20ap_jina')
      end subroutine rate_ne20ap_jina
    
! rne20ag, ne20(a,g)mg24                    

      subroutine rate_ne20ag_nacre(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, rev, bb, gs
         
         include 'formats'

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         if (tf% T9 <= 1d0) then         
            gs   = 8.72d0 * pow(tf% T9,-0.532d0) * exp(-8.995d0*(tf% T9i)) 
         else   
            gs   = 3.74d2 * pow(tf% T9,2.229d0) * exp(-12.681d0*(tf% T9i)) 
         end if      
         bb   = 7.787d0 * exp(-19.821d0*(tf% T9i) - 0.114d0*(tf% T9))
         if (bb > 1) then ! guard against rate going negative
            bb  = 1
         end if      
         term    = gs * (1 - bb)
         call rnacre_rev(tf, &  ! a0 T932 exp(-a1/T9)
            6.010d10, 108.11d0, &  ! a0, a1
            rev)     
         fr    = term
         rr    = rev * term
         
      end subroutine rate_ne20ag_nacre

      subroutine rate_ne20ag_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!       he4 ne20 mg24                       nacrr     9.31600d+00          
         call jina_reaclib_2_1(ihe4, ine20, img24, tf, fr, rr, 'rate_ne20ag_jina')
      end subroutine rate_ne20ag_jina



! rne20ga, ne20(g,a)o16  
   ! see ro16ag
            
! rne20gp, ne20(g,p)f19       
   ! see rf19pg
   
! rne22pg, ne22(p,g)na23

      subroutine rate_ne22pg_nacre(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, rev, bb, cc, dd, gs

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         if (tf% T9 <= 12d0) then   
            ! a0 T9i23 exp(-a1 T9i13 - (T9*a2)^2) 
            !     * (1 + b0 T9 + b1 T92 + b2 T93 + b3 T94 + b4 T95) 
            ! + c0 T9i32 exp(-c1/T9)
            ! + d0 T9i32 exp(-d1/T9)
            ! + e0 T9^e1 exp(-e2/T9)
            call rnacre(tf,  &
               0d0, 0d0, 0d0, & ! a0, a1, a2
               0d0, 0d0, 0d0, 0d0, 0d0, & ! b0, b1, b2, b3, b4
               1.11d-9, 0.422d0, &  ! c0, c1
               6.83d-5, 0.810d0, & ! d0, d1
               8.51d4, 0.725d0, 4.315d0, & ! e0, e1, e2
            gs)        
            cc   = 9.76d-3 * (tf% T9i32) * exp(-1.187d0*(tf% T9i)) 
            dd   = 1.06d-1 * (tf% T9i32) * exp(-1.775d0*(tf% T9i)) 
            gs   = gs + cc + dd
         else   
            gs   = 6.30d4 * pow(tf% T9, 0.816d0) * exp(-3.910d0*(tf% T9i)) 
         end if  
         bb   = 1.410d0 * exp(-14.651d0*(tf% T9i) - 0.020d0*(tf% T9))
         if (bb > 1) then ! guard against rate going negative
            bb  = 1
         end if    
         term    = gs * (1 - bb)
         call rnacre_rev(tf, &  ! a0 T932 exp(-a1/T9)
            4.668d9, 102.05d0, &  ! a0, a1
            rev)      
         fr    = term 
         rr    = rev * term
      end subroutine rate_ne22pg_nacre
      
      
      subroutine rate_ne22pg_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!         p ne22 na23                       ha01r     8.79400d+00          
         call jina_reaclib_2_1(ih1, ine22, ina23, tf, fr, rr, 'rate_ne22pg_jina')
      end subroutine rate_ne22pg_jina
      
! ne22(n,g)ne23      

      subroutine rate_ne22ag_fxt(tf, temp, fr, rr)

      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr

      real(dp) term, rev, aa, bb, cc, dd, res1,  &
                       fT9a, fpT9a, gT9x,  &
                       T9a, T9a13, T9a56,  &
                       rdmass, res2, zz
      parameter        (rdmass = 22.0d0*4.0d0/26.0d0,  &
                        res2   = -11.604d0 * 22.0d0/26.0d0)

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          

! ne22(a, g)mg26
! kappeler 1994 apj 437, 396 

      if ((tf% T9) .lt. 1.25d0) then
       res1 = 1.54d-01*pow((tf% T9)*rdmass,-1.5d0)
       aa    = 1.7d-36 * res1 * exp(res2*(tf% T9i)*0.097d0)
       bb    = 1.5d-7 * res1 * exp(res2*(tf% T9i)*0.400d0)
       cc    = 0.5d0 * res1 * 3.7d-2 * exp(res2*(tf% T9i)*0.633d0)
       dd    = res1 * 3.6d+1 * exp(res2*(tf% T9i)*0.828d0)
       term    = aa + bb + cc + dd
! cf88
      else
       aa    = 1.0d0 + 0.0548d0*(tf% T9)
       zz    = 1.0d0/aa
       T9a   = (tf% T9) *zz
       t9a13  = pow(t9a,oneth)
       t9a56  = pow(t9a,fivsix)
       aa     = 0.197d0/T9a
       bb     = pow(aa,4.82d0)
       fT9a   = exp(-bb)
       aa     = T9a/0.249d0
       bb     = pow(aa,2.31d0)
       fpT9a  = exp(-bb)
       aa     = 5.0d0*exp(-14.791d0*(tf% T9i))
       gT9x   = 1.0d0 + aa
       zz     = 1.0d0/gT9x
       aa     = 4.16d19 * fpT9a*zz
       bb     = 2.08d16 * fT9a*zz
       term    = (aa+bb) * T9a56 * (tf% T9i32) * exp(-47.004d0/T9a13)
      end if
      fr    = term 
      rev    = 6.15d+10 * (tf% T932) * exp(-123.151d0*(tf% T9i))
      rr    = rev * term
      end subroutine rate_ne22ag_fxt


      subroutine rate_ne22ag_jina(tf, temp, fr, rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
!       he4 ne22 mg26                       nacr      1.06150d+01          
         call jina_reaclib_2_1(ihe4, ine22, img26, tf, fr, rr, 'rate_ne22ag_jina')
      end subroutine rate_ne22ag_jina
      
      subroutine rate_na23pa_fxt(tf, temp, fr, rr)

      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr

      real(dp) term, dtermdt, rev, drevdt, aa, daa, bb, dbb, cc, dcc,  &
                       dd, ddd, ee, dee, ff, dff, gg, dgg, hh, dhh, theta, q1, q2
      parameter        (theta = 0.1d0,  &
                        q1    = 1.0d0/0.0169d0,  &
                        q2    = 1.0d0/0.017161d0)
         
         include 'formats.dek'

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
! na23(p, a)ne20
! el eid & champagne 1995 
      if ((tf% T9) <= 2.0d0) then
       aa  = 1.26d+10 * (tf% T9i23) * exp(-20.758d0*(tf% T9i13) - (tf% T92)*q1)
       bb  = 1.0d0  + 0.02d0*(tf% T913) - 13.8d0*(tf% T923) - 1.93d0*(tf% T9)  &
             + 234.0d0*(tf% T943) + 83.6d0*(tf% T953)
       cc   = aa * bb
       dd   = 4.38d0*(tf% T9i32) * exp(-1.979d0*(tf% T9i))
       ee   = 6.50d+06 * pow(tf% T9,-1.366d0) * exp(-6.490d0*(tf% T9i))
       ff   = 1.19d+08 * pow(tf% T9,1.055d0) * exp(-11.411d0*(tf% T9i))
       gg   = theta * 9.91d-14 * (tf% T9i32) * exp(-0.418d0*(tf% T9i))
       term    = cc + dd + ee + ff + gg 

! cf88 + one term from gorres, wiesher & rolfs 1989, apj 343, 365
      else 
       aa  = 8.56d+09 * (tf% T9i23) * exp(-20.766d0*(tf% T9i13) - (tf% T92)*q2)
       bb  = 1.0d0  + 0.02d0*(tf% T913) + 8.21d0*(tf% T923) + 1.15d0*(tf% T9)  &
             + 44.36d0*(tf% T943) + 15.84d0*(tf% T953)
       cc   = aa * bb
       dd   = 4.02d0*(tf% T9i32) * exp(-1.99d0*(tf% T9i))
       ee   = 1.18d+04 * (tf% T9i54) * exp(-3.148d0*(tf% T9i))
       ff   = 8.59d+05 * (tf% T943) * exp(-4.375d0*(tf% T9i))
       gg   = theta * 3.06d-12 * (tf% T9i32) * exp(-0.447d0*(tf% T9i))
       hh   = theta * 0.820d0*(tf% T9i32) * exp(-1.601d0*(tf% T9i))
       term    = cc + dd + ee + ff + gg + hh
      end if

      fr    = term 
      rev      = 1.25d0 * exp(-27.606d0*(tf% T9i))
      rr    = rev * term
      
      if (.false. .and. tf% t9 > 1.3d-2 .and. tf% t9 < 1.8d-2) then
         write(*,1) 'rate_na23pa_fxt', fr, rr, tf% t9
      end if
      
      end subroutine rate_na23pa_fxt


      subroutine rate_na23pa_jina(tf, temp, fr, rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
!         p na23  he4 ne20                  ha04n     2.37900d+00          
         call jina_reaclib_2_2(ih1, ina23, ihe4, ine20, tf, fr, rr, 'rate_na23pa_jina')
      end subroutine rate_na23pa_jina


      subroutine rate_na23pg_fxt(tf, temp, fr, rr)

      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr

      real(dp) term, dtermdt, rev, drevdt, aa, daa, bb, dbb, cc, dcc,  &
                       dd, ddd, ee, dee, ff, dff, gg, dgg, hh, hhi, xx, dxx,  &
                       theta, q1
      parameter        (theta = 0.1d0,  &
                        q1    = 1.0d0/0.088209d0)


         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          

! na23(p, g)mg24

! el eid & champagne 1995 
      if ((tf% T9) <= 2.0d0) then
       aa  = 2.47d+09 * (tf% T9i23) * exp(-20.758d0*(tf% T9i13))
       bb  = 9.19d+01 * (tf% T9i32) * exp(-2.789d0*(tf% T9i))
       cc  = 1.72d+04 * (tf% T9i32) * exp(-3.433d0*(tf% T9i))
       dd  = 3.44d+04 * pow(tf% T9, 0.323d0) * exp(-5.219d0*(tf% T9i))
       ee   = theta * 2.34d-04 * (tf% T9i32) * exp(-1.590d0*(tf% T9i))
       term    = aa + bb + cc + dd + ee 

! cf88 + gorres, wiesher & rolfs 1989, apj 343, 365
      else 
 
       aa  = 2.93d+08 * (tf% T9i23) * exp(-20.766d0*(tf% T9i13) - (tf% T92)*q1)
       bb  = 1.0d0 + 0.02d0*(tf% T913) + 1.61d0*(tf% T923) + 0.226d0*(tf% T9)  &
             + 4.94d0*(tf% T943) + 1.76d0*(tf% T953)
       xx  = aa * bb
       cc   = 9.34d+01 * (tf% T9i32) * exp(-2.789d0*(tf% T9i))
       dd   = 1.89d+04 * (tf% T9i32) * exp(-3.434d0*(tf% T9i))
       ee   = 5.1d+04 * (tf% T915) * exp(-5.51d0*(tf% T9i))
       ff   = theta * 0.820d0*(tf% T9i32) * exp(-1.601d0*(tf% T9i))
       gg   = 1.5d0 * exp(-5.105d0*(tf% T9i))
       hh   = 1.0d0 + gg
       hhi  = 1.0d0/hh
       term    = (xx + cc + dd + ee + ff) * hhi
      end if

      fr    = term 
      rev      = 7.49d+10 * (tf% T932) * exp(-135.665d0*(tf% T9i))
      rr    = rev * term
      end subroutine rate_na23pg_fxt
      
      
      subroutine rate_na23pg_jina(tf, temp, fr, rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
      real(dp) :: fr1
      rr = 0
!         p na23 mg24                       ha04r     1.16910d+01          
         call jina_reaclib_2_1(ih1, ina23, img24, tf, fr, rr, 'rate_na23pg_jina')
      end subroutine rate_na23pg_jina

! Magnesium

! rmg24ag, mg24(a,g)si28

      subroutine rate_mg24ag_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, dtermdt, aa, daa, bb, dbb, cc, dcc, dd, ddd, ee, dee,  &
                       ff, dff, gg, dgg, hh, hhi, rev, drevdt, rc121
      parameter        (rc121 = 0.1d0)

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         aa    = 4.78d+01 * (tf% T9i32) * exp(-13.506d0*(tf% T9i)) 
         bb    =  2.38d+03 * (tf% T9i32) * exp(-15.218d0*(tf% T9i))
         cc    = 2.47d+02 * (tf% T932) * exp(-15.147d0*(tf% T9i)) 
         dd    = rc121 * 1.72d-09 * (tf% T9i32) * exp(-5.028d0*(tf% T9i))
         ee    = rc121* 1.25d-03 * (tf% T9i32) * exp(-7.929d0*(tf% T9i))
         ff    = rc121 * 2.43d+01 * (tf% T9i) * exp(-11.523d0*(tf% T9i))
         gg    = 5.0d0*exp(-15.882d0*(tf% T9i))
         hh    = 1.0d0 + gg
         hhi   = 1.0d0/hh
         term    = (aa + bb + cc + dd + ee + ff) * hhi
         fr    = term
         rev      = 6.27d+10 * (tf% T932) * exp(-115.862d0*(tf% T9i))
         rr    = rev * term
      end subroutine rate_mg24ag_fxt


      subroutine rate_mg24ag_jina(tf, temp, fr, rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
!       he4 mg24 si28                       cf88r     9.98400d+00          
         call jina_reaclib_2_1(ihe4, img24, isi28, tf, fr, rr, 'rate_mg24ag_jina')
      end subroutine rate_mg24ag_jina


! rmg24ga, mg24(g,a)ne20
   ! see rne20ag
             
! rmg24ap, mg24(a,p)al27

      subroutine rate_mg24ap_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, dtermdt, aa, daa, bb, dbb, cc, dcc, dd, ddd, ee, dee,  &
                       ff, dff, gg, dgg, term1, dterm1, term2, dterm2,  &
                       rev, drevdt, rc148, q1
         parameter        (rc148 = 0.1d0,  &
                        q1    = 1.0d0/0.024649d0)                    

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         aa     = 1.10d+08 * (tf% T9i23) * exp(-23.261d0*(tf% T9i13) - (tf% T92)*q1)
         bb     =  1.0d0 + 0.018d0*(tf% T913) + 12.85d0*(tf% T923) + 1.61d0*(tf% T9)   &
               + 89.87d0*(tf% T943) + 28.66d0*(tf% T953)
         term1  = aa * bb
         aa     = 129.0d0 * (tf% T9i32) * exp(-2.517d0*(tf% T9i)) 
         bb     = 5660.0d0 * (tf% T972) * exp(-3.421d0*(tf% T9i)) 
         cc     = rc148 * 3.89d-08 * (tf% T9i32) * exp(-0.853d0*(tf% T9i))  
         dd     = rc148 * 8.18d-09 * (tf% T9i32) * exp(-1.001d0*(tf% T9i))
         term2  = aa + bb + cc + dd
         ee     = oneth*exp(-9.792d0*(tf% T9i))
         ff     =  twoth * exp(-11.773d0*(tf% T9i))
         gg     = 1.0d0 + ee + ff
         term    = (term1 + term2)/gg
         rev      = 1.81d0 * exp(-18.572d0*(tf% T9i))
         fr    = rev * term
         rr    = term
      end subroutine rate_mg24ap_fxt


      subroutine rate_mg24ap_jina(tf, temp, fr, rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
      real(dp) :: fr1
      rr = 0
!       he4 mg24    p al27                  il01rv   -1.60060d+00          
         call jina_reaclib_2_2(ih1, ial27, ihe4, img24, tf, rr, fr, 'rate_mg24ap_jina')
      end subroutine rate_mg24ap_jina

! Aluminum 

! ral27pg, al27(p,g)si28   

      subroutine rate_al27pg_c96(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, rev, aa, bb, cc, dd, ee, ff, gg

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         aa  = 1.32d+09 * (tf% T9i23) * exp(-23.26d0*(tf% T9i13))
         bb  = 3.22d-10 * (tf% T9i32) * exp(-0.836d0*(tf% T9i))*0.17d0
         cc  = 1.74d+00 * (tf% T9i32) * exp(-2.269d0*(tf% T9i))
         dd  = 9.92d+00 * (tf% T9i32) * exp(-2.492d0*(tf% T9i))
         ee  = 4.29d+01 * (tf% T9i32) * exp(-3.273d0*(tf% T9i))
         ff  = 1.34d+02 * (tf% T9i32) * exp(-3.654d0*(tf% T9i))
         gg  = 1.77d+04 * pow(tf% T9, 0.53d0) * exp(-4.588d0*(tf% T9i))
         term    = aa + bb + cc + dd + ee + ff + gg
         fr    = term 
         rev   = 1.13d+11 * (tf% T932) * exp(-134.434d0*(tf% T9i))
         rr    = rev * term
      end subroutine rate_al27pg_c96
      
      subroutine rate_al27pg_jina(tf, temp, fr, rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
      real(dp) :: fr1
!         p al27 si28                       il01r     1.15860d+01          
         call jina_reaclib_2_1(ih1, ial27, isi28, tf, fr, rr, 'rate_al27pg_jina')
      end subroutine rate_al27pg_jina

! Silicon 
 
! rsi28ag, si28(a,g)s32      
      subroutine rate_si28ag_jina(tf, temp, fr, rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
!       he4 si28  s32                       rath      6.94800d+00          
         call jina_reaclib_2_1(ihe4, isi28, is32, tf, fr, rr, 'rate_si28ag_jina')
         !if (abs(temp - 3.0097470376051402D+09) < 1d2) then
         !   include 'formats'
         !   write(*,1) 'rate_si28ag_jina', fr, rr, temp
         !end if
      end subroutine rate_si28ag_jina

      subroutine rate_si28ag_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, dtermdt, aa, daa, rev, drevdt, z, z2, z3
         
         include 'formats'

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         z     = min((tf% T9), 10.0d0)
         z2    = z*z
         z3    = z2*z
         aa    = 1.0d0 + 6.340d-2*z + 2.541d-3*z2 - 2.900d-4*z3
         term    = 4.82d+22 * (tf% T9i23) * exp(-61.015d0*(tf% T9i13) * aa)
         fr    = term
         rev      = 6.461d+10 * (tf% T932) * exp(-80.643d0*(tf% T9i))
         rr    = rev * term

         !if (abs(temp - 3.0097470376051402D+09) < 1d2) then
         !   write(*,1) 'rate_si28ag_fxt', fr, rr, temp
         !end if
         
      end subroutine rate_si28ag_fxt


! rsi28ga, si28(g,a)mg24
   ! see rmg24ag
          
! rsi28ap, si28(a,p)p31   

      subroutine rate_si28ap_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, dtermdt, aa, daa, rev, drevdt, z, z2, z3

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         z     = min((tf% T9), 10.0d0)
         z2    = z*z
         z3    = z2*z
         aa    = 1.0d0 + 2.798d-3*z + 2.763d-3*z2 - 2.341d-4*z3
         term    = 4.16d+13 * (tf% T9i23) * exp(-25.631d0*(tf% T9i13) * aa)
         rev      = 0.5825d0 * exp(-22.224d0*(tf% T9i))
         fr    = rev * term
         rr    = term
      end subroutine rate_si28ap_fxt


      subroutine rate_si28ap_jina(tf, temp, fr, rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
      real(dp) :: fr1
!       he4 si28    p  p31                  il01rv   -1.91710d+00          
         call jina_reaclib_2_2(ih1, ip31, ihe4, isi28, tf, rr, fr, 'rate_si28ap_jina')
      end subroutine rate_si28ap_jina


! rsi28gp, si28(g,p)al27
   ! see ral27pg

! Phosphorus 

! rp31pg, p31(p,g)s32  

      subroutine rate_p31pg_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, dtermdt, aa, daa, rev, drevdt, z, z2, z3
         include 'formats.dek'

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         z     = min((tf% T9), 10.0d0)
         z2    = z*z
         z3    = z2*z
         aa    = 1.0d0 + 1.928d-1*z - 1.540d-2*z2 + 6.444d-4*z3
         term    = 1.08d+16 * (tf% T9i23) * exp(-27.042d0*(tf% T9i13) * aa)
         fr    = term
         rev      = 3.764d+10 * (tf% T932) * exp(-102.865d0*(tf% T9i))
         rr    = rev * term
      end subroutine rate_p31pg_fxt


      subroutine rate_p31pg_jina(tf, temp, fr, rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
!         p  p31  s32                       il01n     8.86400d+00          
         call jina_reaclib_2_1(ih1, ip31, is32, tf, fr, rr, 'rate_p31pg_jina')
      end subroutine rate_p31pg_jina


! rp31pa, p31(p,a)si28  
   ! see rsi28ap
   

! Sulfur 
      
      
! rs32ag, s32(a,g)ar36      
      subroutine rate_s32ag_jina(tf, temp, fr, rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
!       he4  s32 ar36                       rath      6.63900d+00          
         call jina_reaclib_2_1(ihe4, is32, iar36, tf, fr, rr, 'rate_s32ag_jina')
      end subroutine rate_s32ag_jina

! rs32ag, s32(a,g)ar36  

      subroutine rate_s32ag_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, dtermdt, aa, daa, rev, drevdt, z, z2, z3

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         z     = min((tf% T9), 10.0d0)
         z2    = z*z
         z3    = z2*z
         aa    = 1.0d0 + 4.913d-2*z + 4.637d-3*z2 - 4.067d-4*z3
         term    = 1.16d+24 * (tf% T9i23) * exp(-66.690d0*(tf% T9i13) * aa)
         fr    = term
         rev      = 6.616d+10 * (tf% T932) * exp(-77.080d0*(tf% T9i))
         rr    = rev * term
      end subroutine rate_s32ag_fxt


! rs32ga, s32(g,a)si28          
   ! see rsi28ag
   
! rs32ap, s32(a,p)cl35    

      subroutine rate_s32ap_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, dtermdt, aa, daa, rev, drevdt, z, z2, z3

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         z     = min((tf% T9), 10.0d0)
         z2    = z*z
         z3    = z2*z
         aa    = 1.0d0 + 1.041d-1*z - 1.368d-2*z2 + 6.969d-4*z3
         term    = 1.27d+16 * (tf% T9i23) * exp(-31.044d0*(tf% T9i13) * aa)
         rev      = 1.144d0 * exp(-21.643d0*(tf% T9i))
         fr    = rev * term
         rr    = term
      end subroutine rate_s32ap_fxt


      subroutine rate_s32ap_jina(tf, temp, fr, rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
!       he4  s32    p cl35                  il01rv   -1.86700d+00          
         call jina_reaclib_2_2(ih1, icl35, ihe4, is32, tf, rr, fr, 'rate_s32ap_jina')
      end subroutine rate_s32ap_jina

! rs32gp, s32(g,p)p31     
   ! see rp31pg


! Chlorine

! rcl35pg, cl35(p,g)ar36

      subroutine rate_cl35pg_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, dtermdt, aa, daa, rev, drevdt

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         aa    = 1.0d0 + 1.761d-1*(tf% T9) - 1.322d-2*(tf% T92) + 5.245d-4*(tf% T93)
         term    =  4.48d+16 * (tf% T9i23) * exp(-29.483d0*(tf% T9i13) * aa)
         fr    = term
         rev      = 7.568d+10*(tf% T932)*exp(-98.722d0*(tf% T9i))
         rr    = rev * term
      end subroutine rate_cl35pg_fxt


      subroutine rate_cl35pg_jina(tf, temp, fr, rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
!         p cl35 ar36                       il01r     8.50600d+00          
         call jina_reaclib_2_1(ih1, icl35, iar36, tf, fr, rr, 'rate_cl35pg_jina')
      end subroutine rate_cl35pg_jina

! rcl35pa, cl35(p,a)s32
   ! see rs32ap     

! Argon 
      
      
! rar36ag, ar36(a,g)ca40      
      subroutine rate_ar36ag_jina(tf, temp, fr, rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr      
!       he4 ar36 ca40                       rath      7.04000d+00          
         call jina_reaclib_2_1(ihe4, iar36, ica40, tf, fr, rr, 'rate_ar36ag_jina')
      end subroutine rate_ar36ag_jina

! rar36ag, ar36(a,g)ca40   

      subroutine rate_ar36ag_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, dtermdt, aa, daa, rev, drevdt, z, z2, z3

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         z     = min((tf% T9), 10.0d0)
         z2    = z*z
         z3    = z2*z
         aa    = 1.0d0 + 1.458d-1*z - 1.069d-2*z2 + 3.790d-4*z3
         term    = 2.81d+30 * (tf% T9i23) * exp(-78.271d0*(tf% T9i13) * aa)
         fr    = term
         rev      = 6.740d+10 * (tf% T932) * exp(-81.711d0*(tf% T9i))
         rr    = rev * term
      end subroutine rate_ar36ag_fxt

! rar36ga, ar36(g,a)s32             
   ! see rs32ag
   
! rar36ap, ar36(a,p)k39  

      subroutine rate_ar36ap_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, dtermdt, aa, daa, rev, drevdt, z, z2, z3

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         z     = min((tf% T9), 10.0d0)
         z2    = z*z
         z3    = z2*z
         aa    = 1.0d0 + 4.826d-3*z - 5.534d-3*z2 + 4.021d-4*z3
         term    = 2.76d+13 * (tf% T9i23) * exp(-34.922d0*(tf% T9i13) * aa)
         rev      = 1.128d0*exp(-14.959d0*(tf% T9i))
         fr    = rev * term
         rr    = term
      end subroutine rate_ar36ap_fxt


      subroutine rate_ar36ap_jina(tf, temp, fr, rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
!       he4 ar36    p  k39                  rath v   -1.28800d+00          
         call jina_reaclib_2_2(ih1, ik39, ihe4, iar36, tf, rr, fr, 'rate_ar36ap_jina')
      end subroutine rate_ar36ap_jina

! rar36gp, ar36(g,p)cl35  
   ! see rcl35pg
   
! Potassium

! rk39pg, k39(p,g)ca40

      subroutine rate_k39pg_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, dtermdt, aa, daa, rev, drevdt, z, z2, z3

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         z     = min((tf% T9), 10.0d0)
         z2    = z*z
         z3    = z2*z
         aa    = 1.0d0 + 1.622d-1*z - 1.119d-2*z2 + 3.910d-4*z3
         term    = 4.09d+16 * (tf% T9i23) * exp(-31.727d0*(tf% T9i13) * aa)
         fr    = term
         rev      = 7.600d+10 * (tf% T932) * exp(-96.657d0*(tf% T9i))
         rr    = rev * term
      end subroutine rate_k39pg_fxt


      subroutine rate_k39pg_jina(tf, temp, fr, rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
!         p  k39 ca40                       rath      8.32800d+00          
         call jina_reaclib_2_1(ih1, ik39, ica40, tf, fr, rr, 'rate_k39pg_jina')
      end subroutine rate_k39pg_jina

! rk39pa, k39(p,a)ar36
   ! see rar36ap

! Calcium 
      
      
! rca40ag, ca40(a,g)ti44      
      subroutine rate_ca40ag_jina(tf, temp, fr, rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
      include 'formats.dek'
!       he4 ca40 ti44                       rath      5.12700d+00          
         call jina_reaclib_2_1(ihe4, ica40, iti44, tf, fr, rr, 'rate_ca40ag_jina')
      end subroutine rate_ca40ag_jina

! rca40ag, ca40(a,g)ti44    

      subroutine rate_ca40ag_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, dtermdt, aa, daa, rev, drevdt, z, z2, z3

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         z     = min((tf% T9), 10.0d0)
         z2    = z*z
         z3    = z2*z
         aa    = 1.0d0 + 1.650d-2*z + 5.973d-3*z2 - 3.889d-04*z3
         term    = 4.66d+24 * (tf% T9i23) * exp(-76.435d0*(tf% T9i13) * aa)
         fr    = term
         rev      = 6.843d+10 * (tf% T932) * exp(-59.510d0*(tf% T9i))
         rr    = rev * term
      end subroutine rate_ca40ag_fxt


! rca40ga, ca40(g,a)ar36      
   ! see rar36ag

! rca40ap, ca40(a,p)sc43(p,g)ti44        

      subroutine rate_ca40ap_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, dtermdt, aa, daa, rev, drevdt, z, z2, z3

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         z     = min((tf% T9), 10.0d0)
         z2    = z*z
         z3    = z2*z
         aa    = 1.0d0 - 1.206d-2*z + 7.753d-3*z2 - 5.071d-4*z3
         term    = 4.54d+14 * (tf% T9i23) * exp(-32.177d0*(tf% T9i13) * aa)
         rev      = 2.229d0 * exp(-40.966d0*(tf% T9i))
         fr    = rev * term
         rr    = term
      end subroutine rate_ca40ap_fxt


      subroutine rate_ca40ap_jina(tf, temp, fr, rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
!       he4 ca40    p sc43                  rath v   -3.52300d+00          
         call jina_reaclib_2_2(ih1, isc43, ihe4, ica40, tf, rr, fr, 'rate_ca40ap_jina')
      end subroutine rate_ca40ap_jina


! rca40ap, ca40(a,p)sc43 
   ! see rsc43pa
   
! rca40gp, ca40(g,p)k39  
   ! see rk39pg
   
! Scandium 

! rsc43pg, sc43(p,g)ti44     

      subroutine rate_sc43pg_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, dtermdt, aa, daa, rev, drevdt, z, z2, z3

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         z     = min((tf% T9), 10.0d0)
         z2    = z*z
         z3    = z2*z
         aa    = 1.0d0 + 1.023d-1*z - 2.242d-3*z2 - 5.463d-5*z3
         term    = 3.85d+16 * (tf% T9i23) * exp(-33.234d0*(tf% T9i13) * aa)
         fr    = term
         rev      = 1.525d+11 * (tf% T932) * exp(-100.475d0*(tf% T9i))
         rr    = rev * term
      end subroutine rate_sc43pg_fxt


      subroutine rate_sc43pg_jina(tf, temp, fr, rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
!         p sc43 ti44                       rath      8.65000d+00          
         call jina_reaclib_2_1(ih1, isc43, iti44, tf, fr, rr, 'rate_sc43pg_jina')
      end subroutine rate_sc43pg_jina


! rsc43pa, sc43(p,a)ca40
   ! see rca40ap
        

! Titanium 
      
      
! rti44ag, ti44(a,g)cr48      
      subroutine rate_ti44ag_jina(tf, temp, fr, rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
!       he4 ti44 cr48                       rath      7.69200d+00          
         call jina_reaclib_2_1(ihe4, iti44, icr48, tf, fr, rr, 'rate_ti44ag_jina')
      end subroutine rate_ti44ag_jina

! rti44ag, ti44(a,g)cr48   

      subroutine rate_ti44ag_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, dtermdt, aa, daa, rev, drevdt, z, z2, z3

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         z     = min((tf% T9), 10.0d0)
         z2    = z*z
         z3    = z2*z
         aa    = 1.0d0 + 1.066d-1*z - 1.102d-2*z2 + 5.324d-4*z3
         term    = 1.37d+26 * (tf% T9i23) * exp(-81.227d0*(tf% T9i13) * aa)
         fr    = term
         rev      = 6.928d+10*(tf% T932)*exp(-89.289d0*(tf% T9i))
         rr    = rev * term
      end subroutine rate_ti44ag_fxt


! rti44ga, ti44(g,a)ca40       
   ! see rca40ag
   
! rti44ap, ti44(a,p)v47

      subroutine rate_ti44ap_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, dtermdt, aa, daa, rev, drevdt, z, z2, z3

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         z     = min((tf% T9), 10.0d0)
         z2    = z*z
         z3    = z2*z
         aa    = 1.0d0 + 2.655d-2*z - 3.947d-3*z2 + 2.522d-4*z3
         term    = 6.54d+20 * (tf% T9i23) * exp(-66.678d0*(tf% T9i13) * aa)
         rev      = 1.104d0 * exp(-4.723d0*(tf% T9i))
         fr    = rev * term
         rr    = term
      end subroutine rate_ti44ap_fxt


      subroutine rate_ti44ap_jina(tf, temp, fr, rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
!       he4 ti44    p  v47                  chw0r    -4.10500d-01          
         call jina_reaclib_2_2(ihe4, iti44, ih1, iv47, tf, fr, rr, 'rate_ti44ap_jina')
      end subroutine rate_ti44ap_jina
      

! rti44gp, ti44(g,p)sc43 
   ! see rsc43pg
   

! Vanadium 

! rv47pg, v47(p,g)cr48     

      subroutine rate_v47pg_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, dtermdt, aa, daa, rev, drevdt, z, z2, z3

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         z     = min((tf% T9), 10.0d0)
         z2    = z*z
         z3    = z2*z
         aa    = 1.0d0 + 9.979d-2*z - 2.269d-3*z2 - 6.662d-5*z3
         term    = 2.05d+17 * (tf% T9i23) * exp(-35.568d0*(tf% T9i13) * aa)
         fr    = term
         rev      = 7.649d+10*(tf% T932)*exp(-93.999d0*(tf% T9i))
         rr    = rev * term
      end subroutine rate_v47pg_fxt


      subroutine rate_v47pg_jina(tf, temp, fr, rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
!         p  v47 cr48                       nfisn     8.10607d+00          
         call jina_reaclib_2_1(ih1, iv47, icr48, tf, fr, rr, 'rate_v47pg_jina')
      end subroutine rate_v47pg_jina


! rv47pa, v47(p,a)ti44 
   ! see rti44ap
   

! Chromium 
      
      
! rcr48ag, cr48(a,g)fe52      
      subroutine rate_cr48ag_jina(tf, temp, fr, rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
!       he4 cr48 fe52                       rath      7.93900d+00          
         call jina_reaclib_2_1(ihe4, icr48, ife52, tf, fr, rr, 'rate_cr48ag_jina')
      end subroutine rate_cr48ag_jina

! rcr48ag, cr48(a,g)fe52    

      subroutine rate_cr48ag_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, dtermdt, aa, daa, rev, drevdt, z, z2, z3

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         z     = min((tf% T9), 10.0d0)
         z2    = z*z
         z3    = z2*z
         aa    = 1.0d0 + 6.325d-2*z - 5.671d-3*z2 + 2.848d-4*z3
         term    = 1.04d+23 * (tf% T9i23) * exp(-81.420d0*(tf% T9i13) * aa)
         fr    = term
         rev      = 7.001d+10 * (tf% T932) * exp(-92.177d0*(tf% T9i))
         rr    = rev * term
      end subroutine rate_cr48ag_fxt


! rcr48ga, cr48(g,a)ti44     
   ! see rti44ag

! rcr48ap, cr48(a,p)mn51 

      subroutine rate_cr48ap_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, dtermdt, aa, daa, rev, drevdt, z, z2, z3

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         z     = min((tf% T9), 10.0d0)
         z2    = z*z
         z3    = z2*z
         aa    = 1.0d0 + 1.384d-2*z + 1.081d-3*z2 - 5.933d-5*z3
         term    = 1.83d+26 * (tf% T9i23) * exp(-86.741d0*(tf% T9i13) * aa)
         fr    = term
         rev      = 0.6087d0*exp(-6.510d0*(tf% T9i))
         rr    = rev * term
      end subroutine rate_cr48ap_fxt


      subroutine rate_cr48ap_jina(tf, temp, fr, rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
!       he4 cr48    p mn51                  rath      5.58000d-01          
         call jina_reaclib_2_2(ihe4, icr48, ih1, imn51, tf, fr, rr, 'rate_cr48ap_jina')
      end subroutine rate_cr48ap_jina


! rcr48gp, cr48(g,p)v47  
   ! see rv47pg


! Manganese 

! rmn51pg, mn51(p,g)fe52     

      subroutine rate_mn51pg_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, dtermdt, aa, daa, rev, drevdt, z, z2, z3

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         z     = min((tf% T9), 10.0d0)
         z2    = z*z
         z3    = z2*z
         aa    = 1.0d0 + 8.922d-2*z - 1.256d-3*z2 - 9.453d-5*z3
         term    = 3.77d+17 * (tf% T9i23) * exp(-37.516d0*(tf% T9i13) * aa)
         fr    = term
         rev      = 1.150d+11*(tf% T932)*exp(-85.667d0*(tf% T9i))
         rr    = rev * term
      end subroutine rate_mn51pg_fxt


      subroutine rate_mn51pg_jina(tf, temp, fr, rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
!         p mn51 fe52                       rath      7.38100d+00          
         call jina_reaclib_2_1(ih1, imn51, ife52, tf, fr, rr, 'rate_mn51pg_jina')
      end subroutine rate_mn51pg_jina


! rmn51pa, mn51(p,a)cr48 
    ! see rcr48ap


! rfe52ag, fe52(a,g)ni56   

      subroutine rate_fe52ag_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, dtermdt, aa, daa, rev, drevdt, z, z2, z3

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         z     = min((tf% T9), 10.0d0)
         z2    = z*z
         z3    = z2*z
         aa    = 1.0d0 + 7.846d-2*z - 7.430d-3*z2 + 3.723d-4*z3
         term    = 1.05d+27 * (tf% T9i23) * exp(-91.674d0*(tf% T9i13) * aa)
         fr    = term
         rev      = 7.064d+10*(tf% T932)*exp(-92.850d0*(tf% T9i))
         rr    = rev * term
      end subroutine rate_fe52ag_fxt


! rfe52ga, fe52(g,a)cr48       
   ! see rcr48ag

! rfe52ap, fe52(a,p)co55 

      subroutine rate_fe52ap_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, dtermdt, aa, daa, rev, drevdt, z, z2, z3

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         z     = min((tf% T9), 10.0d0)
         z2    = z*z
         z3    = z2*z
         aa    = 1.0d0 + 1.367d-2*z + 7.428d-4*z2 - 3.050d-5*z3
         term    = 1.30d+27 * (tf% T9i23) * exp(-91.674d0*(tf% T9i13) * aa)
         fr    = term
         rev      = 0.4597d0*exp(-9.470d0*(tf% T9i))
         rr    = rev * term
      end subroutine rate_fe52ap_fxt


! rfe52gp, fe52(g,p)mn51 
   ! see mg51pg
   


      subroutine rate_fe52ng_jina(tf, temp, fr,  rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
!         n fe52 fe53                       rath      1.06840d+01          
         call jina_reaclib_2_1(ineut, ife52, ife53, tf, fr, rr, 'rate_fe52ng_jina')
      end subroutine rate_fe52ng_jina   


      subroutine rate_fe52ng_fxt(tf, temp, fr,  rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
      real(dp) term, dtermdt, rev, drevdt, tq2

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
! fe52(n, g)fe53
      tq2     = (tf% T9) - 0.348d0
      term    = 9.604d+05 * exp(-0.0626d0*tq2)
      fr    = term
      rev      = 2.43d+09 * (tf% T932) * exp(-123.951d0*(tf% T9i)) 
      rr    = rev * term
      end subroutine rate_fe52ng_fxt


      subroutine rate_fe53ng_fxt(tf, temp, fr, rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
      real(dp) term, dtermdt, rev, drevdt, tq1, tq10, dtq10, tq2

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
! fe53(n, g)fe54
      tq1   = (tf% T9)/0.348d0
      tq10  = pow(tq1, 0.10d0)
      tq2   = (tf% T9) - 0.348d0
      term    = 1.817d+06 * tq10 * exp(-0.06319d0*tq2)
      fr    = term
      rev      = 1.56d+11 * (tf% T932) * exp(-155.284d0*(tf% T9i))
      rr    = rev * term
      end subroutine rate_fe53ng_fxt


      subroutine rate_fe53ng_jina(tf, temp, fr,  rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
!         n fe53 fe54                       rath      1.33780d+01          
         call jina_reaclib_2_1(ineut, ife53, ife54, tf, fr, rr, 'rate_fe53ng_jina')
      end subroutine rate_fe53ng_jina   


      subroutine rate_fe54pg_fxt(tf, temp, fr, rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
      real(dp) term, dtermdt, rev, drevdt, aa, daa, z, z2, z3

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
! fe54(p, g)co55
      z     = min((tf% T9), 10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 9.593d-2*z - 3.445d-3*z2 + 8.594d-5*z3
      term    = 4.51d+17 * (tf% T9i23) * exp(-38.483d0*(tf% T9i13) * aa)
      fr    = term
      rev      = 2.400d+09 * (tf% T932) * exp(-58.605d0*(tf% T9i))
      rr    = rev * term
      end subroutine rate_fe54pg_fxt


      subroutine rate_fe54a_jina(tf, temp, fr,  rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
      real(dp) :: fr1, fr2, fr3
      include 'formats.dek'
      call rate_fe54ag_jina(tf, temp, fr1,  rr)
      call rate_fe54an_jina(tf, temp, fr2,  rr)
      call rate_fe54ap_jina(tf, temp, fr3,  rr)
      fr = fr1 + fr2 + fr3
      end subroutine rate_fe54a_jina

      subroutine rate_fe54an_jina(tf, temp, fr,  rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
            call jina_reaclib_2_2(ineut, ini57, ihe4, ife54, tf, rr, fr, 'rate_fe54an_jina')
      end subroutine rate_fe54an_jina


      subroutine rate_fe54ng_jina(tf, temp, fr,  rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
         call jina_reaclib_2_1(ineut, ife54, ife55, tf, fr, rr, 'rate_fe54ng_jina')
      end subroutine rate_fe54ng_jina



      subroutine rate_fe55ng_jina(tf, temp, fr,  rr)
      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: temp
      real(dp), intent(out) :: fr, rr
         call jina_reaclib_2_1(ineut, ife55, ife56, tf, fr, rr, 'rate_fe55ng_jina')
      end subroutine rate_fe55ng_jina


! Cobalt 

! rco55pg, co55(p,g)ni56     

      subroutine rate_co55pg_fxt(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) term, dtermdt, aa, daa, rev, drevdt, z, z2, z3

         if (tf% t9 < lowT9_cutoff) then
            fr = 0; rr = 0; return
         end if 
          
         z     = min((tf% T9), 10.0d0)
         z2    = z*z
         z3    = z2*z
         aa    = 1.0d0 + 9.894d-2*z - 3.131d-3*z2 - 2.160d-5*z3
         term    = 1.21d+18 * (tf% T9i23) * exp(-39.604d0*(tf% T9i13) * aa)
         fr    = term
         rev      = 1.537d+11*(tf% T932)*exp(-83.382d0*(tf% T9i))
         rr    = rev * term
      end subroutine rate_co55pg_fxt


! rco55pa, co55(p,a)fe52  
   ! see rfe52ap

! Nickel 

! rni56ga, ni56(g,a)fe52   
   ! see rfe52ag
   
! rni56gp, ni56(g,p)co55 
   ! see rco55pg
   
      subroutine rate_v44pg_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         real(dp) :: fr1
!         p  v44 cr45                       rath      3.10000d+00          
         call jina_reaclib_2_1(ih1, iv44, icr45, tf, fr, rr, 'rate_v44pg_jina')
      end subroutine rate_v44pg_jina


      subroutine rate_v45pg_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!         p  v45 cr46                       rath      4.88600d+00          
         call jina_reaclib_2_1(ih1, iv45, icr46, tf, fr, rr, 'rate_v45pg_jina')
      end subroutine rate_v45pg_jina

      subroutine rate_co53pg_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!         p co53 ni54                       rath      3.85600d+00          
         call jina_reaclib_2_1(ih1, ico53, ini54, tf, fr, rr, 'rate_co53pg_jina')
      end subroutine rate_co53pg_jina


      subroutine rate_co54pg_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!         p co54 ni55                       rath      4.61400d+00          
         call jina_reaclib_2_1(ih1, ico54, ini55, tf, fr, rr, 'rate_co54pg_jina')
      end subroutine rate_co54pg_jina

      subroutine rate_ga62pg_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!         p ga62 ge63                       nfisn     2.23867d+00          
         call jina_reaclib_2_1(ih1, iga62, ige63, tf, fr, rr, 'rate_ga62pg_jina')
      end subroutine rate_ga62pg_jina


      subroutine rate_ga63pg_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
!         p ga63 ge64                       rath      5.02500d+00          
         call jina_reaclib_2_1(ih1, iga63, ige64, tf, fr, rr, 'rate_ga63pg_jina')
      end subroutine rate_ga63pg_jina
      
      ! ni56
      
      subroutine rate_fe52ag_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
            call jina_reaclib_2_1(ihe4, ife52, ini56, tf, fr, rr, 'rate_fe52ag_jina')
      end subroutine rate_fe52ag_jina


      subroutine rate_fe52ap_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         call jina_reaclib_2_2(ihe4, ife52, ih1, ico55, tf, fr, rr, 'rate_fe52ap_jina')
      end subroutine rate_fe52ap_jina


      subroutine rate_co55pg_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         call jina_reaclib_2_1(ih1, ico55, ini56, tf, fr, rr, 'rate_co55pg_jina')
      end subroutine rate_co55pg_jina


      subroutine rate_fe54pg_jina(tf, temp, fr,  rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
            call jina_reaclib_2_1(ih1, ife54, ico55, tf, fr, rr, 'rate_fe54pg_jina')
      end subroutine rate_fe54pg_jina   

      ! ni58
      
      subroutine rate_fe54ag_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
            call jina_reaclib_2_1(ihe4, ife54, ini58, tf, fr, rr, 'rate_fe54ag_jina')
      end subroutine rate_fe54ag_jina


      subroutine rate_fe54ap_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         call jina_reaclib_2_2(ihe4, ife54, ih1, ico57, tf, fr, rr, 'rate_fe54ap_jina')
      end subroutine rate_fe54ap_jina

      subroutine rate_fe56pg_jina(tf, temp, fr,  rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
            call jina_reaclib_2_1(ih1, ife56, ico57, tf, fr, rr, 'rate_fe56pg_jina')
      end subroutine rate_fe56pg_jina   

 
      subroutine rate_c12_c12_to_h1_na23_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         call jina_reaclib_2_2( &
            ic12, ic12, ih1, ina23, tf, fr, rr, 'rate_c12_c12_to_h1_na23_jina')
      end subroutine rate_c12_c12_to_h1_na23_jina
      

      subroutine rate_he4_ne20_to_c12_c12_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         call jina_reaclib_2_2( &
            ihe4, ine20, ic12, ic12, tf, fr, rr, 'rate_he4_ne20_to_c12_c12_jina')
      end subroutine rate_he4_ne20_to_c12_c12_jina
      

      subroutine rate_he4_mg24_to_c12_o16_jina(tf, temp, fr, rr)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: temp
         real(dp), intent(out) :: fr, rr
         call jina_reaclib_2_2( &
            ihe4, img24, ic12, io16, tf, fr, rr, 'rate_he4_mg24_to_c12_o16_jina')
      end subroutine rate_he4_mg24_to_c12_o16_jina
      
      

      subroutine tfactors(tf, logT_in, temp_in)
      use const_def, only: ln10
! sets various popular temperature factors
! this routine must be called before any of the rates are called

      type (T_Factors), pointer :: tf
      real(dp), intent(in) :: logT_in, temp_in

      real(dp) :: logT, temp
      
      logT = max(logT_in, 0d0)
      temp = max(temp_in, 1d0)
      
      tf% lnT9 = (logT - 9)*ln10
      tf% T9    = temp * 1.0d-9
      tf% T92   = tf% T9 * tf% T9
      tf% T93   = tf% T9 * tf% T92
      tf% T94   = tf% T9 * tf% T93
      tf% T95   = tf% T9 * tf% T94
      tf% T96   = tf% T9 * tf% T95

      tf% T912  = sqrt(tf% T9)
      tf% T932  = tf% T9 * tf% T912
      tf% T952  = tf% T9 * tf% T932
      tf% T972  = tf% T9 * tf% T952

      tf% T913  = pow(tf% T9,oneth)
      tf% T923  = tf% T913 * tf% T913
      tf% T943  = tf% T9 * tf% T913
      tf% T953  = tf% T9 * tf% T923
      tf% T973  = tf% T953 * tf% T923
      tf% T9113 = tf% T973 * tf% T943

      tf% T914  = pow(tf% T9, 0.25d0)
      tf% T934  = tf% T914 * tf% T914 * tf% T914
      tf% T954  = tf% T9 * tf% T914
      tf% T974  = tf% T9 * tf% T934

      tf% T915  = pow(tf% T9,onefif)
      tf% T935  = tf% T915 * tf% T915 * tf% T915
      tf% T945  = tf% T915 * tf% T935
      tf% T965  = tf% T9 * tf% T915

      tf% T916  = pow(tf% T9, onesix)
      tf% T976  = tf% T9 * tf% T916
      tf% T9i76 = 1.0d0 / tf% T976

      tf% T917  = pow(tf% T9, onesev)
      tf% T927  = tf% T917 * tf% T917
      tf% T947  = tf% T927 * tf% T927

      tf% T918  = sqrt(tf% T914)
      tf% T938  = tf% T918 * tf% T918 * tf% T918
      tf% T958  = tf% T938 * tf% T918 * tf% T918

      tf% T9i   = 1.0d0 / tf% T9
      tf% T9i2  = tf% T9i * tf% T9i
      tf% T9i3  = tf% T9i2 * tf% T9i

      tf% T9i12 = 1.0d0 / tf% T912
      tf% T9i32 = tf% T9i * tf% T9i12
      tf% T9i52 = tf% T9i * tf% T9i32
      tf% T9i72 = tf% T9i * tf% T9i52

      tf% T9i13 = 1.0d0 / tf% T913
      tf% T9i23 = tf% T9i13 * tf% T9i13
      tf% T9i43 = tf% T9i * tf% T9i13
      tf% T9i53 = tf% T9i * tf% T9i23

      tf% T9i14 = 1.0d0 / tf% T914
      tf% T9i34 = tf% T9i14 * tf% T9i14 * tf% T9i14
      tf% T9i54 = tf% T9i * tf% T9i14

      tf% T9i15 = 1.0d0 / tf% T915
      tf% T9i35 = tf% T9i15 * tf% T9i15 * tf% T9i15
      tf% T9i45 = tf% T9i15 * tf% T9i35
      tf% T9i65 = tf% T9i * tf% T9i15

      tf% T9i17 = 1.0d0 / tf% T917
      tf% T9i27 = tf% T9i17 * tf% T9i17 
      tf% T9i47 = tf% T9i27 * tf% T9i27

      tf% T9i18 = 1.0d0 / tf% T918
      tf% T9i38 = tf% T9i18 * tf% T9i18 * tf% T9i18
      tf% T9i58 = tf% T9i38 * tf% T9i18 * tf% T9i18

      end subroutine tfactors

      
      subroutine show_nacre_terms( &
               tf, a0, a1, a2, b0, b1, b2, b3, b4, c0, c1, d0, d1, e0, e1, e2)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: a0, a1, a2, b0, b1, b2, b3, b4,  &
               c0, c1, d0, d1, e0, e1, e2
         
         include 'formats.dek'
         real(dp) :: aa, bb, cc, dd, ee
         aa   = a0 * (tf% T9i23) * exp(-a1*(tf% T9i13) - (tf% T92)*(a2*a2)) 
         bb   = 1 + b0*(tf% T9) + b1*(tf% T92) +  b2*(tf% T93) +   b3*(tf% T94) +   b4*(tf% T95)
         cc   = c0 * (tf% T9i32) * exp(-c1*(tf% T9i)) 
         dd   = d0 * (tf% T9i32) * exp(-d1*(tf% T9i)) 
         ee   = e0 * pow(tf% T9,e1) * exp(-e2*(tf% T9i)) 

         write(*,1) 'aa', aa
         write(*,1) 'bb', bb
         write(*,1) 'aa*bb', aa*bb
         write(*,1) 'cc', cc
         write(*,1) 'dd', dd
         write(*,1) 'ee', ee
         
      end subroutine show_nacre_terms

     
      subroutine rnacre( &
         ! a0 T9i23 exp(-a1 T9i13 - (T9*a2)^2) 
         !     * (1 + b0 T9 + b1 T92 + b2 T93 + b3 T94 + b4 T95) 
         ! + c0 T9i32 exp(-c1/T9)
         ! + d0 T9i32 exp(-d1/T9)
         ! + e0 T9^e1 exp(-e2/T9) &
               tf, a0, a1, a2, b0, b1, b2, b3, b4, c0, c1, d0, d1, e0, e1, e2, term)
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: a0, a1, a2, b0, b1, b2, b3, b4,  &
               c0, c1, d0, d1, e0, e1, e2
         real(dp), intent(out) :: term         
         real(dp) :: aa, daa, bb, dbb, cc, dcc, dd, ddd, ee, dee         
         aa   = a0 * (tf% T9i23) * exp(-a1*(tf% T9i13) - (tf% T92)*(a2*a2)) 
         bb   = 1 + b0*(tf% T9) + b1*(tf% T92) +  b2*(tf% T93) +   b3*(tf% T94) +   b4*(tf% T95)
         if (bb < 0) then
            bb = 0
         end if
         cc   = c0 * (tf% T9i32) * exp(-c1*(tf% T9i)) 
         dd   = d0 * (tf% T9i32) * exp(-d1*(tf% T9i)) 
         ee   = e0 * pow(tf% T9,e1) * exp(-e2*(tf% T9i)) 
         term    = aa * bb + cc + dd + ee
      end subroutine rnacre

      
      subroutine rnacre_rev(tf, a0, a1, rev) ! a0 T932 exp(-a1/T9)
         real(dp), intent(in) :: a0, a1
         real(dp), intent(out) :: rev
         type (T_Factors), pointer :: tf
         rev    = a0 * (tf% T932) * exp(-a1*(tf% T9i))
      end subroutine rnacre_rev      
      
      
      subroutine jina_reaclib_1_1(i1, o1, tf, fr, rr, str)
         integer, intent(in) :: i1, o1
         type (T_Factors), pointer :: tf
         real(dp), intent(out) :: fr, rr
         character (len=*), intent(in) :: str
         integer :: ierr
         call try1_reaclib_1_1(i1, o1, tf, fr, rr, str, ierr)
         if (ierr == 0) return
         ierr = 0
         call try1_reaclib_1_1(o1, i1, tf, rr, fr, str, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in jina_reaclib_1_1 ' // trim(str), tf% T9
            call mesa_error(__FILE__,__LINE__)
         end if
      end subroutine jina_reaclib_1_1
      
      
      subroutine jina_reaclib_1_2(i1, o1, o2, tf, fr, rr, str)
         integer, intent(in) :: i1, o1, o2
         type (T_Factors), pointer :: tf
         real(dp), intent(out) :: fr, rr
         character (len=*), intent(in) :: str
         integer :: ierr
         call try1_reaclib_1_2(i1, o1, o2, tf, fr, rr, str, ierr)
         if (ierr == 0) return
         ierr = 0
         call try1_reaclib_2_1(o1, o2, i1, tf, rr, fr, str, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in jina_reaclib_1_2 ' // trim(str), tf% T9
            call mesa_error(__FILE__,__LINE__)
         end if
      end subroutine jina_reaclib_1_2
      
      
      subroutine jina_reaclib_1_3(i1, o1, o2, o3, tf, fr, rr, str)
         integer, intent(in) :: i1, o1, o2, o3
         type (T_Factors), pointer :: tf
         real(dp), intent(out) :: fr, rr
         character (len=*), intent(in) :: str
         integer :: ierr
         call try1_reaclib_1_3(i1, o1, o2, o3, tf, fr, rr, str, ierr)
         if (ierr == 0) return
         ierr = 0
         call try1_reaclib_3_1(o1, o2, o3, i1, tf, rr, fr, str, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in jina_reaclib_1_3 ' // trim(str), tf% T9
            call mesa_error(__FILE__,__LINE__)
         end if
      end subroutine jina_reaclib_1_3
      
      
      subroutine jina_reaclib_1_4(i1, o1, o2, o3, o4, tf, fr, rr, str)
         integer, intent(in) :: i1, o1, o2, o3, o4
         type (T_Factors), pointer :: tf
         real(dp), intent(out) :: fr, rr
         character (len=*), intent(in) :: str
         integer :: ierr
         call try1_reaclib_1_4(i1, o1, o2, o3, o4, tf, fr, rr, str, ierr)
         !if (ierr == 0) return
         !ierr = 0
         !call try1_reaclib_4_1(o1, o2, o3, o4, i1, tf, rr, fr, str, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in jina_reaclib_1_4 ' // trim(str), tf% T9
            call mesa_error(__FILE__,__LINE__)
         end if
      end subroutine jina_reaclib_1_4
      
      
      subroutine jina_reaclib_2_1(i1, i2, o1, tf, fr, rr, str)
         integer, intent(in) :: i1, i2, o1
         type (T_Factors), pointer :: tf
         real(dp), intent(out) :: fr, rr
         character (len=*), intent(in) :: str
         integer :: ierr
         call try1_reaclib_2_1(i1, i2, o1, tf, fr, rr, str, ierr)
         if (ierr == 0) return
         ierr = 0
         call try1_reaclib_1_2(o1, i1, i2, tf, rr, fr, str, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in jina_reaclib_2_1 ' // trim(str), tf% T9
            call mesa_error(__FILE__,__LINE__)
         end if
      end subroutine jina_reaclib_2_1
      
      
      subroutine jina_reaclib_2_2(i1, i2, o1, o2, tf, fr, rr, str)
         integer, intent(in) :: i1, i2, o1, o2
         type (T_Factors), pointer :: tf
         real(dp), intent(out) :: fr, rr
         character (len=*), intent(in) :: str
         integer :: ierr
         call try1_reaclib_2_2(i1, i2, o1, o2, tf, fr, rr, str, ierr)
         if (ierr == 0) return
         ierr = 0
         call try1_reaclib_2_2(o1, o2, i1, i2, tf, rr, fr, str, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in jina_reaclib_2_2 ' // trim(str), tf% T9
            call mesa_error(__FILE__,__LINE__)
         end if
      end subroutine jina_reaclib_2_2
      
      
      subroutine jina_reaclib_2_3(i1, i2, o1, o2, o3, tf, fr, rr, str)
         integer, intent(in) :: i1, i2, o1, o2, o3
         type (T_Factors), pointer :: tf
         real(dp), intent(out) :: fr, rr
         character (len=*), intent(in) :: str
         integer :: ierr
         call try1_reaclib_2_3(i1, i2, o1, o2, o3, tf, fr, rr, str, ierr)
         if (ierr == 0) return
         ierr = 0
         call try1_reaclib_3_2(o1, o2, o3, i1, i2, tf, rr, fr, str, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in jina_reaclib_3_2 ' // trim(str), tf% T9
            call mesa_error(__FILE__,__LINE__)
         end if
      end subroutine jina_reaclib_2_3
      
      
      subroutine jina_reaclib_2_4(i1, i2, o1, o2, o3, o4, tf, fr, rr, str)
         integer, intent(in) :: i1, i2, o1, o2, o3, o4
         type (T_Factors), pointer :: tf
         real(dp), intent(out) :: fr, rr
         character (len=*), intent(in) :: str
         integer :: ierr
         call try1_reaclib_2_4(i1, i2, o1, o2, o3, o4, tf, fr, rr, str, ierr)
         if (ierr == 0) return
         ierr = 0
         call try1_reaclib_4_2(o1, o2, o3, o4, i1, i2, tf, rr, fr, str, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in jina_reaclib_2_4 ' // trim(str), tf% T9
            call mesa_error(__FILE__,__LINE__)
         end if
      end subroutine jina_reaclib_2_4
      
      
      subroutine jina_reaclib_3_1(i1, i2, i3, o1, tf, fr, rr, str)
         integer, intent(in) :: i1, i2, i3, o1
         type (T_Factors), pointer :: tf
         real(dp), intent(out) :: fr, rr
         character (len=*), intent(in) :: str
         integer :: ierr
         call try1_reaclib_3_1(i1, i2, i3, o1, tf, fr, rr, str, ierr)
         if (ierr == 0) return
         ierr = 0
         call try1_reaclib_3_1(o1, i1, i2, i3, tf, rr, fr, str, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in jina_reaclib_3_1 ' // trim(str), tf% T9
            call mesa_error(__FILE__,__LINE__)
         end if
      end subroutine jina_reaclib_3_1
      
      
      subroutine jina_reaclib_3_2(i1, i2, i3, o1, o2, tf, fr, rr, str)
         integer, intent(in) :: i1, i2, i3, o1, o2
         type (T_Factors), pointer :: tf
         real(dp), intent(out) :: fr, rr
         character (len=*), intent(in) :: str
         integer :: ierr
         call try1_reaclib_3_2(i1, i2, i3, o1, o2, tf, fr, rr, str, ierr)
         if (ierr == 0) return
         ierr = 0
         call try1_reaclib_2_3(o1, o2, i1, i2, i3, tf, rr, fr, str, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in jina_reaclib_3_2 ' // trim(str), tf% T9
            call mesa_error(__FILE__,__LINE__)
         end if
      end subroutine jina_reaclib_3_2
      
      
      subroutine jina_reaclib_4_2(i1, i2, i3, i4, o1, o2, tf, fr, rr, str)
         integer, intent(in) :: i1, i2, i3, i4, o1, o2
         type (T_Factors), pointer :: tf
         real(dp), intent(out) :: fr, rr
         character (len=*), intent(in) :: str
         integer :: ierr
         call try1_reaclib_4_2(i1, i2, i3, i4, o1, o2, tf, fr, rr, str, ierr)
         if (ierr == 0) return
         ierr = 0
         call try1_reaclib_2_4(o1, o2, i1, i2, i3, i4, tf, rr, fr, str, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in jina_reaclib_4_2 ' // trim(str), tf% T9
            call mesa_error(__FILE__,__LINE__)
         end if
      end subroutine jina_reaclib_4_2
      
      
      subroutine try1_reaclib_1_1(i1, o1, tf, fr, rr, str, ierr)
         integer, intent(in) :: i1, o1
         type (T_Factors), pointer :: tf
         real(dp), intent(out) :: fr, rr
         character (len=*), intent(in) :: str
         integer, intent(out) :: ierr
         integer :: num_in, num_out
         integer, dimension(4) :: nuclides_in, nuclides_out
         logical, parameter :: dbg = .false.
         include 'formats.dek'
         ierr = 0
         num_in = 1
         nuclides_in(1) = i1
         num_out = 1
         nuclides_out(1) = o1
         call reaclib_rate(  &
            str, num_in, nuclides_in, num_out, nuclides_out, tf% T9, &
            fr, rr, ierr)
      end subroutine try1_reaclib_1_1
      
      
      subroutine try1_reaclib_1_2(i1, o1, o2, tf, fr, rr, str, ierr)
         integer, intent(in) :: i1, o1, o2
         type (T_Factors), pointer :: tf
         real(dp), intent(out) :: fr, rr
         character (len=*), intent(in) :: str
         integer, intent(out) :: ierr
         integer :: num_in, num_out
         integer, dimension(4) :: nuclides_in, nuclides_out
         logical, parameter :: dbg = .false.
         include 'formats.dek'
         ierr = 0
         num_in = 1
         nuclides_in(1) = i1
         num_out = 2
         nuclides_out(1) = o1
         nuclides_out(2) = o2
         call reaclib_rate(  &
            str, num_in, nuclides_in, num_out, nuclides_out, tf% T9, &
            fr, rr, ierr)
      end subroutine try1_reaclib_1_2
      
      
      subroutine try1_reaclib_1_3(i1, o1, o2, o3, tf, fr, rr, str, ierr)
         integer, intent(in) :: i1, o1, o2, o3
         type (T_Factors), pointer :: tf
         real(dp), intent(out) :: fr, rr
         character (len=*), intent(in) :: str
         integer, intent(out) :: ierr
         integer :: num_in, num_out
         integer, dimension(4) :: nuclides_in, nuclides_out
         logical, parameter :: dbg = .false.
         include 'formats.dek'
         ierr = 0
         num_in = 1
         nuclides_in(1) = i1
         num_out = 3
         nuclides_out(1) = o1
         nuclides_out(2) = o2
         nuclides_out(3) = o3
         call reaclib_rate(  &
            str, num_in, nuclides_in, num_out, nuclides_out, tf% T9, &
            fr, rr, ierr)
      end subroutine try1_reaclib_1_3
      
      
      subroutine try1_reaclib_1_4(i1, o1, o2, o3, o4, tf, fr, rr, str, ierr)
         integer, intent(in) :: i1, o1, o2, o3, o4
         type (T_Factors), pointer :: tf
         real(dp), intent(out) :: fr, rr
         character (len=*), intent(in) :: str
         integer, intent(out) :: ierr
         integer :: num_in, num_out
         integer, dimension(4) :: nuclides_in, nuclides_out
         logical, parameter :: dbg = .false.
         include 'formats.dek'
         ierr = 0
         num_in = 1
         nuclides_in(1) = i1
         num_out = 4
         nuclides_out(1) = o1
         nuclides_out(2) = o2
         nuclides_out(3) = o3
         nuclides_out(4) = o4
         call reaclib_rate(  &
            str, num_in, nuclides_in, num_out, nuclides_out, tf% T9, &
            fr, rr, ierr)
      end subroutine try1_reaclib_1_4
      
      
      subroutine try1_reaclib_2_1(i1, i2, o1, tf, fr, rr, str, ierr)
         integer, intent(in) :: i1, i2, o1
         type (T_Factors), pointer :: tf
         real(dp), intent(out) :: fr, rr
         character (len=*), intent(in) :: str
         integer, intent(out) :: ierr
         integer :: num_in, num_out
         integer, dimension(4) :: nuclides_in, nuclides_out
         logical, parameter :: dbg = .false.
         include 'formats.dek'
         ierr = 0
         num_in = 2
         nuclides_in(1) = i1
         nuclides_in(2) = i2
         num_out = 1
         nuclides_out(1) = o1
         call reaclib_rate(  &
            str, num_in, nuclides_in, num_out, nuclides_out, tf% T9, &
            fr, rr, ierr)
      end subroutine try1_reaclib_2_1
      
      
      subroutine try1_reaclib_2_2(i1, i2, o1, o2, tf, fr, rr, str, ierr)
         integer, intent(in) :: i1, i2, o1, o2
         type (T_Factors), pointer :: tf
         real(dp), intent(out) :: fr, rr
         character (len=*), intent(in) :: str
         integer, intent(out) :: ierr
         integer :: num_in, num_out
         integer, dimension(4) :: nuclides_in, nuclides_out
         logical, parameter :: dbg = .false.
         include 'formats.dek'
         ierr = 0
         num_in = 2
         nuclides_in(1) = i1
         nuclides_in(2) = i2
         num_out = 2
         nuclides_out(1) = o1
         nuclides_out(2) = o2
         call reaclib_rate(  &
            str, num_in, nuclides_in, num_out, nuclides_out, tf% T9, &
            fr, rr, ierr)
      end subroutine try1_reaclib_2_2
      
      
      subroutine try1_reaclib_2_3(i1, i2, o1, o2, o3, tf, fr, rr, str, ierr)
         integer, intent(in) :: i1, i2, o1, o2, o3
         type (T_Factors), pointer :: tf
         real(dp), intent(out) :: fr, rr
         character (len=*), intent(in) :: str
         integer, intent(out) :: ierr
         integer :: num_in, num_out
         integer, dimension(4) :: nuclides_in, nuclides_out
         logical, parameter :: dbg = .false.
         include 'formats.dek'
         ierr = 0
         num_in = 2
         nuclides_in(1) = i1
         nuclides_in(2) = i2
         num_out = 3
         nuclides_out(1) = o1
         nuclides_out(2) = o2
         nuclides_out(3) = o3
         call reaclib_rate(  &
            str, num_in, nuclides_in, num_out, nuclides_out, tf% T9, &
            fr, rr, ierr)
      end subroutine try1_reaclib_2_3
      
      
      subroutine try1_reaclib_2_4(i1, i2, o1, o2, o3, o4, tf, fr, rr, str, ierr)
         integer, intent(in) :: i1, i2, o1, o2, o3, o4
         type (T_Factors), pointer :: tf
         real(dp), intent(out) :: fr, rr
         character (len=*), intent(in) :: str
         integer, intent(out) :: ierr
         integer :: num_in, num_out
         integer, dimension(4) :: nuclides_in, nuclides_out
         logical, parameter :: dbg = .false.
         include 'formats.dek'
         ierr = 0
         num_in = 2
         nuclides_in(1) = i1
         nuclides_in(2) = i2
         num_out = 4
         nuclides_out(1) = o1
         nuclides_out(2) = o2
         nuclides_out(3) = o3
         nuclides_out(4) = o4
         call reaclib_rate(  &
            str, num_in, nuclides_in, num_out, nuclides_out, tf% T9, &
            fr, rr, ierr)
      end subroutine try1_reaclib_2_4
      
      
      subroutine try1_reaclib_3_1(i1, i2, i3, o1, tf, fr, rr, str, ierr)
         integer, intent(in) :: i1, i2, i3, o1
         type (T_Factors), pointer :: tf
         real(dp), intent(out) :: fr, rr
         character (len=*), intent(in) :: str
         integer, intent(out) :: ierr
         integer :: num_in, num_out
         integer, dimension(4) :: nuclides_in, nuclides_out
         logical, parameter :: dbg = .false.
         include 'formats.dek'
         ierr = 0
         num_in = 3
         nuclides_in(1) = i1
         nuclides_in(2) = i2
         nuclides_in(3) = i3
         num_out = 1
         nuclides_out(1) = o1
         call reaclib_rate(  &
            str, num_in, nuclides_in, num_out, nuclides_out, tf% T9, &
            fr, rr, ierr)
      end subroutine try1_reaclib_3_1
      
      
      subroutine try1_reaclib_3_2(i1, i2, i3, o1, o2, tf, fr, rr, str, ierr)
         integer, intent(in) :: i1, i2, i3, o1, o2
         type (T_Factors), pointer :: tf
         real(dp), intent(out) :: fr, rr
         character (len=*), intent(in) :: str
         integer, intent(out) :: ierr
         integer :: num_in, num_out
         integer, dimension(4) :: nuclides_in, nuclides_out
         logical, parameter :: dbg = .false.
         include 'formats.dek'
         ierr = 0
         num_in = 3
         nuclides_in(1) = i1
         nuclides_in(2) = i2
         nuclides_in(3) = i3
         num_out = 2
         nuclides_out(1) = o1
         nuclides_out(2) = o2
         call reaclib_rate(  &
            str, num_in, nuclides_in, num_out, nuclides_out, tf% T9, &
            fr, rr, ierr)
      end subroutine try1_reaclib_3_2
      
      
      subroutine try1_reaclib_4_2(i1, i2, i3, i4, o1, o2, tf, fr, rr, str, ierr)
         integer, intent(in) :: i1, i2, i3, i4, o1, o2
         type (T_Factors), pointer :: tf
         real(dp), intent(out) :: fr, rr
         character (len=*), intent(in) :: str
         integer, intent(out) :: ierr
         integer :: num_in, num_out
         integer, dimension(4) :: nuclides_in, nuclides_out
         logical, parameter :: dbg = .false.
         include 'formats.dek'
         ierr = 0
         num_in = 4
         nuclides_in(1) = i1
         nuclides_in(2) = i2
         nuclides_in(3) = i3
         nuclides_in(4) = i4
         num_out = 2
         nuclides_out(1) = o1
         nuclides_out(2) = o2
         call reaclib_rate(  &
            str, num_in, nuclides_in, num_out, nuclides_out, tf% T9, &
            fr, rr, ierr)
      end subroutine try1_reaclib_4_2
      
      
      subroutine reaclib_rate( &
            str, num_in, nuclides_in, num_out, nuclides_out, T9, &
            lambda, rlambda, ierr)
         use reaclib_support, only: reaction_handle
         character (len=*), intent(in) :: str
         integer, intent(in) :: num_in, nuclides_in(:) 
         integer, intent(in) :: num_out, nuclides_out(:)
         real(dp), intent(in) :: T9
         real(dp), intent(out) :: lambda, rlambda
         integer, intent(out) :: ierr        
         character (len=max_id_length) :: handle
         integer :: lo, hi, i
         integer :: iso_ids(num_in+num_out)
         include 'formats.dek'
         ierr = 0
         iso_ids(1:num_in) = nuclides_in(1:num_in)
         iso_ids(1+num_in:num_out+num_in) = nuclides_out(1:num_out)
         call reaction_handle(num_in, num_out, iso_ids, '-', handle)
         call reaclib_rate_for_handle(handle, T9, lambda, rlambda, ierr)
      end subroutine reaclib_rate
      
      
      subroutine reaclib_rate_for_handle(handle, T9, lambda, rlambda, ierr)
         use reaclib_eval, only: do_reaclib_indices_for_reaction, do_reaclib_reaction_rates
         character (len=*), intent(in) :: handle
         real(dp), intent(in) :: T9
         real(dp), intent(out) :: lambda, rlambda
         integer, intent(out) :: ierr                 
         real(dp) :: dlambda_dlnT, drlambda_dlnT         
         integer :: lo, hi, i
         logical, parameter :: forward_only = .false.
         include 'formats.dek'         
         ierr = 0
         call do_reaclib_indices_for_reaction(handle, reaclib_rates, lo, hi, ierr)
         if (ierr /= 0) then
            return
            write(*,'(a)')  &
               'reaclib_rate_and_dlnT_for_handle: failed in reaclib_indices_for_reaction '  &
               // trim(handle)
            stop 'reaclib_rate_and_dlnT_for_handle'
         end if
         if (T9 < reaclib_min_T9 .and. reaclib_rates% reaction_flag(lo) /= 'w' .and. &
             reaclib_rates% reaction_flag(lo) /= 'e') then ! w or ec
            lambda = 0; rlambda = 0
            return
         end if
         call do_reaclib_reaction_rates(  &
            lo, hi, T9, reaclib_rates, chem_isos, forward_only, &
            lambda, dlambda_dlnT,  &
            rlambda, drlambda_dlnT,  &
            ierr)      
         if (ierr /= 0) then
            write(*,'(a)')  &
               'reaclib_rate_for_handle: failed in reaclib_reaction_rates '  &
               // trim(handle)
            stop 'reaclib_rate_for_handle'
            return
         end if
      end subroutine reaclib_rate_for_handle
      
      
      subroutine reaclib_rate_and_dlnT_for_handle( &
            handle, T9, lambda, dlambda_dlnT, rlambda, drlambda_dlnT, ierr)
         use reaclib_eval, only: do_reaclib_indices_for_reaction, do_reaclib_reaction_rates
         character (len=*), intent(in) :: handle
         real(dp), intent(in) :: T9
         real(dp), intent(out) :: lambda, dlambda_dlnT, rlambda, drlambda_dlnT
         integer, intent(out) :: ierr                 
         integer :: lo, hi
         include 'formats.dek'         
         ierr = 0
         call do_reaclib_indices_for_reaction(handle, reaclib_rates, lo, hi, ierr)
         if (ierr /= 0) then
            return
            write(*,'(a)')  &
               'reaclib_rate_and_dlnT_for_handle: failed in reaclib_indices_for_reaction '  &
               // trim(handle)
            stop 'reaclib_rate_and_dlnT_for_handle'
         end if
         call reaclib_rate_and_dlnT( &
            lo, hi, handle, T9, lambda, dlambda_dlnT, rlambda, drlambda_dlnT, ierr)
      end subroutine reaclib_rate_and_dlnT_for_handle

         
      subroutine reaclib_rate_and_dlnT( &
            lo, hi, handle, T9, lambda, dlambda_dlnT, rlambda, drlambda_dlnT, ierr)
         use reaclib_eval, only: do_reaclib_reaction_rates
         integer, intent(in) :: lo,hi
         character (len=*), intent(in) :: handle
         real(dp), intent(in) :: T9
         real(dp), intent(out) :: lambda, dlambda_dlnT, rlambda, drlambda_dlnT
         integer, intent(out) :: ierr                 
         logical, parameter :: forward_only = .false.
         include 'formats.dek'         
         ierr = 0
         if (T9 < reaclib_min_T9 .and. reaclib_rates% reaction_flag(lo) /= 'w') then
            lambda = 0; dlambda_dlnT = 0
            rlambda = 0; drlambda_dlnT = 0
            return
         end if
         call do_reaclib_reaction_rates(  &
            lo, hi, T9, reaclib_rates, chem_isos, forward_only, &
            lambda, dlambda_dlnT, rlambda, drlambda_dlnT,  &
            ierr)      
         if (ierr /= 0) then
            write(*,*) 'failed in reaclib_reaction_rates ' // trim(handle)
            return
         end if
         return
         write(*,3) 'reaclib_rate_and_dlnT_for_handle ' // trim(handle), lo, hi
         write(*,1) 'T9', T9
         write(*,1) 'lambda', lambda
         write(*,1) 'dlambda_dlnT', dlambda_dlnT
         write(*,1) 'rlambda', rlambda
         write(*,1) 'drlambda_dlnT', drlambda_dlnT
         write(*,*)
      end subroutine reaclib_rate_and_dlnT




      subroutine do_reaclib(tf, a1, a2, a3, a4, a5, a6, a7, term)
         use const_def
         type (T_Factors), pointer :: tf
         real(dp), intent(in) :: a1, a2, a3, a4, a5, a6, a7
         real(dp), intent(out) :: term
         real(dp) :: exponent
         include 'formats.dek'
         !rate = exp(a1 + a2/T9 + a3/T913 + a4*T913 + a5*T9 + a6*T953 + a7*ln(T9))
         if (tf% T9 < reaclib_min_T9) then
            term = 0; return
         end if
         exponent = a1 + a2*(tf% T9i) + a3*(tf% T9i13) + a4*(tf% T913) +  &
               a5*(tf% T9) + a6*(tf% T953) + a7*(tf% lnT9)
         term = exp(exponent)
         if (is_bad(term)) then
            write(*,1) 'exponent', exponent
            write(*,1) 'term', term
            write(*,1) 'a1', a1
            write(*,1) 'a2*(tf% T9i)', a2*(tf% T9i)
            write(*,1) 'a3*(tf% T9i13)', a3*(tf% T9i13)
            write(*,1) 'a4*(tf% T913)', a4*(tf% T913)
            write(*,1) 'a5*(tf% T9)', a5*(tf% T9)
            write(*,1) 'a6*(tf% T953)', a6*(tf% T953)
            write(*,1) 'a7*(tf% lnT9)', a7*(tf% lnT9)
            write(*,1) 'tf% lnT9/ln10', tf% lnT9/ln10
            stop 'reaclib'
         end if
      end subroutine do_reaclib


      ! returns the indx corresponding to Tpart just less than T9
      ! T9 is the temperature in units of GK
      ! returns a value of 0 or npart if value is less than the minimum or maximum of the partition function
      ! temperature array Tpart
      function get_partition_fcn_indx(T9) result(indx)
         real(dp), intent(in) :: T9
         integer :: indx
         integer, parameter :: max_iterations = 8
         integer :: low, high, mid, i
         low = 1
         high = npart
         if (T9 < Tpart(low)) then
            indx = low-1
            return
         end if
         if (T9 > Tpart(high)) then
            indx = high + 1
         end if
         do i = 1, max_iterations
            if (high-low <= 1) then
               indx = low
               return
            end if
            mid = (high+low)/2
            if (T9 < Tpart(mid)) then
               high = mid
            else
               low = mid
            end if
         end do
         ! should never get here
         indx = low-1
      end function get_partition_fcn_indx


      subroutine mazurek_init(ierr)
         use rates_def, only: tv,rv,rfdm,rfd0,rfd1,rfd2,tfdm,tfd0,tfd1,tfd2
         integer, intent(out) :: ierr
         integer :: k,j
         rv(:) = (/ 6D0, 7D0, 8D0, 9D0, 10D0, 11D0 /)
         tv(:) = (/ 2D0, 4D0, 6D0, 8D0, 10D0, 12D0, 14D0 /)
         ierr = 0
         do k=2,4 
            rfdm(k)=1.d0/((rv(k-1)-rv(k))*(rv(k-1)-rv(k+1))*(rv(k-1)-rv(k+2))) 
            rfd0(k)=1.d0/((rv(k)-rv(k-1))*(rv(k)-rv(k+1))*(rv(k)-rv(k+2))) 
            rfd1(k)=1.d0/((rv(k+1)-rv(k-1))*(rv(k+1)-rv(k))*(rv(k+1)-rv(k+2))) 
            rfd2(k)=1.d0/((rv(k+2)-rv(k-1))*(rv(k+2)-rv(k))*(rv(k+2)-rv(k+1))) 
         enddo
         do j=2,5 
            tfdm(j)=1.d0/((tv(j-1)-tv(j))*(tv(j-1)-tv(j+1))*(tv(j-1)-tv(j+2))) 
            tfd0(j)=1.d0/((tv(j)-tv(j-1))*(tv(j)-tv(j+1))*(tv(j)-tv(j+2))) 
            tfd1(j)=1.d0/((tv(j+1)-tv(j-1))*(tv(j+1)-tv(j))*(tv(j+1)-tv(j+2))) 
            tfd2(j)=1.d0/((tv(j+2)-tv(j-1))*(tv(j+2)-tv(j))*(tv(j+2)-tv(j+1))) 
         enddo
      end subroutine mazurek_init

      subroutine mazurek(btemp,bden,y56,ye,rn56ec,sn56ec)       
      use rates_def, only: tv,rv,rfdm,rfd0,rfd1,rfd2,tfdm,tfd0,tfd1,tfd2
      real(dp), intent(in) :: btemp,bden,y56,ye
      real(dp), intent(out) :: rn56ec,sn56ec

!  this routine evaluates mazurek's 1973 fits for the ni56 electron 
!  capture rate rn56ec and neutrino loss rate sn56ec 

!  input: 
!  y56 = nickel56 molar abundance
!  ye  = electron to baryon number, zbar/abar

!  output:
!  rn56ec = ni56 electron capture rate
!  sn56ec = ni56 neutrino loss rate

!  declare 
      integer          ifirst,jp,kp,jr,jd,ii,ik,ij,j,k 
      real(dp) rnt(2),rne(2,7),datn(2,6,7),  &
                       t9,r,rfm,rf0,rf1,rf2,dfacm,dfac0,dfac1,dfac2,  &
                       tfm,tf0,tf1,tf2,tfacm,tfac0,tfac1,tfac2

!  initialize 
      data (((datn(ii,ik,ij),ik=1,6),ij=1,7),ii=1,1) /  &
          -3.98d0, -2.84d0, -1.41d0,  0.20d0,  1.89d0,  3.63d0,  &
          -3.45d0, -2.62d0, -1.32d0,  0.22d0,  1.89d0,  3.63d0,  &
          -2.68d0, -2.30d0, -1.19d0,  0.27d0,  1.91d0,  3.62d0,  &
          -2.04d0, -1.87d0, -1.01d0,  0.34d0,  1.94d0,  3.62d0,  &
          -1.50d0, -1.41d0, -0.80d0,  0.45d0,  1.99d0,  3.60d0,  &
          -1.00d0, -0.95d0, -0.54d0,  0.60d0,  2.06d0,  3.58d0,  &
          -0.52d0, -0.49d0, -0.21d0,  0.79d0,  2.15d0,  3.55d0 / 
      data (((datn(ii,ik,ij),ik=1,6),ij=1,7),ii=2,2) /  &
          -3.68d0, -2.45d0, -0.80d0,  1.12d0,  3.13d0,  5.19d0,  &
          -2.91d0, -2.05d0, -0.64d0,  1.16d0,  3.14d0,  5.18d0,  &
          -1.95d0, -1.57d0, -0.40d0,  1.24d0,  3.16d0,  5.18d0,  &
          -1.16d0, -0.99d0, -0.11d0,  1.37d0,  3.20d0,  5.18d0,  &
          -0.48d0, -0.40d0,  0.22d0,  1.54d0,  3.28d0,  5.16d0,  &
           0.14d0,  0.19d0,  0.61d0,  1.78d0,  3.38d0,  5.14d0,  &
           0.75d0,  0.78d0,  1.06d0,  2.07d0,  3.51d0,  5.11d0 / 

!  calculate ni56 electron capture and neutrino loss rates 
      rn56ec = 0.0d0
      sn56ec = 0.0d0

      if (btemp*1d-9 < lowT9_cutoff) return
          
      if ( (btemp .lt. 2.0d9) .or. (bden*ye .lt. 1.0d6)) return
      t9    = min(btemp, 1.4d10) * 1.0d-9
      r     = max(6.0d0,min(11.0d0,log10(bden*ye))) 
      jp    = min(max(2,int(0.5d0*t9)), 5) 
      kp    = min(max(2,int(r)-5), 4) 
      rfm   = r - rv(kp-1) 
      rf0   = r - rv(kp) 
      rf1   = r - rv(kp+1) 
      rf2   = r - rv(kp+2) 
      dfacm = rf0*rf1*rf2*rfdm(kp) 
      dfac0 = rfm*rf1*rf2*rfd0(kp) 
      dfac1 = rfm*rf0*rf2*rfd1(kp) 
      dfac2 = rfm*rf0*rf1*rfd2(kp) 
      tfm   = t9 - tv(jp-1) 
      tf0   = t9 - tv(jp) 
      tf1   = t9 - tv(jp+1) 
      tf2   = t9 - tv(jp+2) 
      tfacm = tf0*tf1*tf2*tfdm(jp) 
      tfac0 = tfm*tf1*tf2*tfd0(jp) 
      tfac1 = tfm*tf0*tf2*tfd1(jp) 
      tfac2 = tfm*tf0*tf1*tfd2(jp) 

!  evaluate the spline fits
      do jr = 1,2 
       do jd = jp-1,jp+2 
        rne(jr,jd) =   dfacm*datn(jr,kp-1,jd) + dfac0*datn(jr,kp,jd)  &
                     + dfac1*datn(jr,kp+1,jd) + dfac2*datn(jr,kp+2,jd) 
       enddo
       rnt(jr) =  tfacm*rne(jr,jp-1) + tfac0*rne(jr,jp)  &
                + tfac1*rne(jr,jp+1) + tfac2*rne(jr,jp+2) 
      enddo

!  set the output
      rn56ec = exp10(rnt(1))
      sn56ec = 6.022548d+23 * 8.18683d-7 * y56 * exp10(rnt(2))
      return 
      end subroutine mazurek

      
      subroutine n14_electron_capture_rate(T,Rho,UE,rate)
         real(dp), intent(in) :: T ! temperature
         real(dp), intent(in) :: Rho ! density
         real(dp), intent(in) :: UE ! electron molecular weight
         real(dp), intent(out) :: rate ! (s^-1)
         
         real(dp) :: Q, AMC2, AMULTIP, AL92, T8, X, XFER, EF, Y, AA, GUESS, ELCAP
         
         ! from Lars
         
      
!      Inputs are T in K, rho in gr/cm^3, and UE=electron mean mol. weight
!
!     Gives a reasonable estimate (i.e. within factor of 50% or so) of the 
!     electron capture rate for electrons on 14N in a plasma assumed to be quite 
!     degenerate. 
!
!         x=KT/Q, y=E_FERMI/Q 
!   
!      ELCAP is the rate in 1/seconds 
!
!
!     Let's start by putting in the Q value, electron rest mass and 
!     temperature in units of keV.
!
!
         Q = 667.479d0
         AMC2 = 510.999d0
         AMULTIP = 0.693d0/1.104d9
         AL92 = LOG(9d0/2d0)
         T8 = T/1d8
         X = 8.617d0*T8/Q
!
!     For this value of the density, find the electron fermi momentum 
!     assuming that the KT corrections to the electron EOS are not
!     important. 
!
      XFER = pow(RHO/(0.9739D6*UE),1d0/3d0) 
!
!      The parameter we need that is used in the fitting formula is
!      the electron Fermi energy 
!
      EF = AMC2*SQRT(1.0D0 + XFER*XFER)
      Y = EF/Q
      IF(Y .LT. (1.0D0 + AL92*X)) THEN 
          AA = (Y-1.0D0)/X
          GUESS = 2.0D0*X*X*X*exp(AA)
      ELSE
          GUESS = pow3(Y-1.0D0+(3.0D0-AL92)*X)/3.0D0 
      ENDIF
!
!     Now multiply by the prefactors .. .
!
      ELCAP = GUESS*AMULTIP


         rate = ELCAP
         
      end subroutine n14_electron_capture_rate


      subroutine ecapnuc(etakep,temp,rpen,rnep,spen,snep)
         use const_def
      real(dp), intent(in) :: etakep,temp
      real(dp), intent(out) :: rpen,rnep,spen,snep

!  given the electron degeneracy parameter etakep (chemical potential
!  without the electron's rest mass divided by kt) and the temperature temp,
!  this routine calculates rates for 
!  electron capture on protons rpen (captures/sec/proton),
!  positron capture on neutrons rnep (captures/sec/neutron), 
!  and their associated neutrino energy loss rates 
!  spen (ergs/sec/proton) and snep (ergs/sec/neutron)

!  declare
      integer          iflag
      real(dp), parameter :: c2me = me * clight * clight
      real(dp), parameter :: cmk5 = (kerg / c2me) * (kerg / c2me) * (kerg / c2me) * (kerg / c2me) * (kerg / c2me)
      real(dp), parameter :: cmk6 = cmk5 * (kerg / c2me)
      real(dp) t9,t5,qn,etaef,etael,zetan,eta,etael2, &
                       etael3,etael4,f1l,f2l,f3l,f4l,f5l,f1g, &
                       f2g,f3g,f4g,f5g,exmeta,eta2,eta3,eta4, &
                       fac0,fac1,fac2,fac3,rie1,rie2,facv0,facv1, &
                       facv2,facv3,facv4,rjv1,rjv2,spenc,snepc, &
                       exeta,zetan2,f0,etael5, &
                       qn1,ft,qn2, &
                       qndeca,tmean
      parameter        (qn1    = -2.0716446d-06, &
                        ft     = 1083.9269d0, &
                        qn2    = 2.0716446d-06, &
                        qndeca = 1.2533036d-06, &
                        tmean  = 886.7d0)
      


!  tmean and qndeca are the mean lifetime and decay energy of the neutron
!  c2me is the constant used to convert the neutrino energy
!  loss rate from mec2/s (as in the paper) to ergs/particle/sec.

!  initialize
      rpen  = 0.0d0
      rnep  = 0.0d0
      spen  = 0.0d0
      snep  = 0.0d0
      t9    = temp * 1.0d-9

      if (t9 < lowT9_cutoff) return

      iflag = 0
      qn    = qn1
          

!  chemical potential including the electron rest mass
      etaef = etakep + c2me/kerg/temp


!  iflag=1 is for electrons,  iflag=2 is for positrons
502   iflag = iflag + 1
      if (iflag.eq.1) etael = qn2/kerg/temp
      if (iflag.eq.2) etael = c2me/kerg/temp
      if (iflag.eq.2) etaef = -etaef

      t5    = temp*temp*temp*temp*temp
      zetan = qn/kerg/temp
      eta   = etaef - etael

!  protect from overflowing with large eta values
      if (eta .le. 6.8d+02) then
       exeta = exp(eta)
      else 
       exeta = 0.0d0
      end if
      etael2 = etael*etael
      etael3 = etael2*etael
      etael4 = etael3*etael
      etael5 = etael4*etael
      zetan2 = zetan*zetan
      if (eta .le. 6.8d+02) then
       f0 = log1p(exeta)
      else
       f0 = eta
      end if

!  if eta le. 0., the following fermi integrals apply
      f1l = exeta
      f2l = 2.0d0   * f1l
      f3l = 6.0d0   * f1l
      f4l = 24.0d0  * f1l
      f5l = 120.0d0 * f1l

!  if eta gt. 0., the following fermi integrals apply:
      f1g = 0.0d0
      f2g = 0.0d0
      f3g = 0.0d0
      f4g = 0.0d0
      f5g = 0.0d0
      if (eta .gt. 0.0d0) then
       exmeta = exp(-eta)
       eta2   = eta*eta
       eta3   = eta2*eta
       eta4   = eta3*eta
       f1g = 0.5d0*eta2 + 2.0d0 - exmeta
       f2g = eta3/3.0d0 + 4.0d0*eta + 2.0d0*exmeta
       f3g = 0.25d0*eta4 + 0.5d0*pi2*eta2 + 12.0d0 - 6.0d0*exmeta
       f4g = 0.2d0*eta4*eta + 2.0d0*pi2/3.0d0*eta3 + 48.0d0*eta &
             + 24.0d0*exmeta
       f5g = eta4*eta2/6.0d0 + 5.0d0/6.0d0*pi2*eta4  &
             + 7.0d0/6.0d0*pi2*eta2  + 240.0d0 -120.d0*exmeta
       end if

!  factors which are multiplied by the fermi integrals
      fac3 = 2.0d0*zetan + 4.0d0*etael
      fac2 = 6.0d0*etael2 + 6.0d0*etael*zetan + zetan2
      fac1 = 4.0d0*etael3 + 6.0d0*etael2*zetan + 2.0d0*etael*zetan2
      fac0 = etael4 + 2.0d0*zetan*etael3 + etael2*zetan2

!  electron capture rates onto protons with no blocking
      rie1 = f4l + fac3*f3l + fac2*f2l + fac1*f1l + fac0*f0
      rie2 = f4g + fac3*f3g + fac2*f2g + fac1*f1g + fac0*f0

!  neutrino emission rate for electron capture:
      facv4 = 5.0d0*etael + 3.0d0*zetan
      facv3 = 10.0d0*etael2 + 12.0d0*etael*zetan + 3.0d0*zetan2
      facv2 = 10.0d0*etael3 + 18.0d0*etael2*zetan &
              + 9.0d0*etael*zetan2 + zetan2*zetan
      facv1 = 5.0d0*etael4 + 12.0d0*etael3*zetan  &
              + 9.0d0*etael2*zetan2 + 2.0d0*etael*zetan2*zetan
      facv0 = etael5 + 3.0d0*etael4*zetan &
              + 3.0d0*etael3*zetan2 + etael2*zetan2*zetan
      rjv1  = f5l + facv4*f4l + facv3*f3l &
              + facv2*f2l + facv1*f1l + facv0*f0
      rjv2  = f5g + facv4*f4g + facv3*f3g &
              + facv2*f2g + facv1*f1g + facv0*f0

!  for electrons capture onto protons
      if (iflag .eq. 2) go to 503
      if (eta .gt. 0d0) go to 505
      rpen  = ln2*cmk5*t5*rie1/ft
      spen  = ln2*cmk6*t5*temp*rjv1/ft
      spenc = ln2*cmk6*t5*temp*rjv1/ft*c2me
      go to 504
505   rpen = ln2*cmk5*t5*rie2/ft
      spen = ln2*cmk6*t5*temp*rjv2/ft
      spenc = ln2*cmk6*t5*temp*rjv2/ft*c2me
504   continue
      qn = qn2
      go to 502

!  for positrons capture onto neutrons
503   if (eta.gt.0d0) go to 507
      rnep  = ln2*cmk5*t5*rie1/ft
      snep  = ln2*cmk6*t5*temp*rjv1/ft
      snepc = ln2*cmk6*t5*temp*rjv1/ft*c2me
!      if (rho.lt.1.0d+06) snep=snep+qndeca*xn(9)/mn/tmean
      go to 506
507   rnep  = ln2*cmk5*t5*rie2/ft
      snep  = ln2*cmk6*t5*temp*rjv2/ft
      snepc = ln2*cmk6*t5*temp*rjv2/ft*c2me
!      if (rho.lt.1.0d+06) snep=snep+qndeca*xn(9)/mn/tmean
506   continue
      return
      end subroutine ecapnuc
    
      end module ratelib
      
      


