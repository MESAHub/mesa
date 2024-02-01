! ***********************************************************************
!
!   Copyright (C) 2011-2019  Bill Paxton & The MESA Team
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
!
! ***********************************************************************

      module ion_tables_eval
      
      use const_def, only: dp, one_sixth
      use ionization_def
      use math_lib
      use utils_lib, only: mesa_error

      implicit none
      
      

      logical, parameter :: dbg = .false.
      
      
      contains


      subroutine Get_ion_Results(&
               Z_in, X_in, arho, alogrho, atemp, alogtemp, &
               res, ierr)
         use const_def
         use utils_lib, only: is_bad
         use ion_tables_load, only: Load_ion_Table
         
         ! INPUT

         real(dp), intent(in) :: Z_in ! the desired Z
         real(dp), intent(in) :: X_in ! the desired X
         
         real(dp), intent(in) :: arho, alogrho ! the density
            ! provide both if you have them.  else pass one and set the other to = arg_not_provided
            
         real(dp), intent(in) :: atemp, alogtemp ! the temperature
            ! provide both if you have them.  else pass one and set the other to = arg_not_provided

         
         ! OUTPUT    
              
         real(dp), intent(inout) :: res(num_ion_vals)
         integer, intent(out) :: ierr
         
         real(dp), dimension(num_ion_vals) :: res1, d_dlnRho_c_T1, d_dlnT_c_Rho1
         real(dp), dimension(num_ion_vals) :: res2, d_dlnRho_c_T2, d_dlnT_c_Rho2
         real(dp) :: X, Z, Rho, logRho, T, logT
         integer :: iregion
         character (len=256) :: message
                  
         real(dp) :: alfa, beta, c_dx, c_dy, logRho_lo, logRho_hi, &
            logT1, logT2, logT7, logT8, logRho3, logRho4
         real(dp) :: logQ, A, B, sin_pi_alfa, dA_dlnT, dA_dlnRho, dB_dlnT, dB_dlnRho
         real(dp) :: d_dx_dlogT, d_dx_dlogRho, d_dy_dlogT, d_dy_dlogRho
         real(dp) ::&
            lnPgas, Pgas, dlnPgas_dlnRho, dlnPgas_dlnT, dPgas_dlnRho, dPgas_dlnT, &
            Prad, dPrad_dlnT, P, dP_dlnRho, dP_dlnT, dlnP_dlnRho, dlnP_dlnT, &
            lnE, energy, dlnE_dlnRho, dlnE_dlnT, &
            lnS, entropy, dlnS_dlnRho, dlnS_dlnT, &
            d_alfa_dlogT, d_alfa_dlogRho, d_beta_dlogT, d_beta_dlogRho
         real(dp), parameter :: dZ_transition = 0.01d0 ! must be > 0
         real(dp), parameter :: logT_margin = 0.1d0
         real(dp), parameter :: tiny = 1d-20
         
         logical :: debug
         
         include 'formats'
         
         ierr = 0
         debug = dbg
         
         if (.not. ion_is_initialized) then
            call Load_ion_Table(ierr)
            if (ierr /= 0) then
               if (dbg) write(*,*) 'Load_ion_Table ierr', ierr
               return
            end if
         end if
         
         if (is_bad(X_in) .or. is_bad(Z_in)) then
            ierr = -1
            return
         end if
         
         X = X_in; Z = Z_in
         if (X < tiny) X = 0
         if (Z < tiny) Z = 0
         
         !..get temp and rho args
         T = atemp; logT = alogtemp
         if (atemp == arg_not_provided .and. alogtemp == arg_not_provided) then
            ierr = -2; return
         end if
         if (alogtemp == arg_not_provided) logT = log10(T)
         if (atemp == arg_not_provided) T = exp10(logT)
         
         if (T <= 0) then
            ierr = -1
            return
         end if
         
         Rho = arho; logrho = alogrho
         if (arho == arg_not_provided .and. alogrho == arg_not_provided) then
            ierr = -3; return
         end if
         if (alogrho == arg_not_provided) logRho = log10(Rho)
         if (arho == arg_not_provided) Rho = exp10(logRho)
         
         if (Rho <= 0) then
            ierr = -1
            return
         end if
         
         if (is_bad(Rho) .or. is_bad(T)) then
            ierr = -1
            return
         end if
         call Get_ion_ZResults(Z, X, Rho, logRho, T, logT, res, ierr)

      end subroutine Get_ion_Results

      
      subroutine Get_ion_ZResults(&
               Z, X, Rho, logRho, T, logT, &
               res, ierr)
         use chem_def
         real(dp), intent(in) :: Z, X, Rho, logRho, T, logT
         real(dp), intent(inout) :: res(num_ion_vals)
         integer, intent(out) :: ierr

         real(dp), dimension(num_ion_vals, num_ion_Zs) :: &
            res_zx, d_dlnRho_c_T_zx, d_dlnT_c_Rho_zx, d_res_dX
         real(dp), dimension(num_ion_vals) :: d_dX, d_dZ
         real(dp) :: d_res_dZ(num_ion_vals), dZ, denom, c(3),&
            dP_dZ, dlnP_dZ
         real(dp), parameter :: tiny = 1d-6

         character (len=256) :: message
         integer :: iz, j, ci
         logical, parameter :: ion_dbg = dbg

         ierr = 0
         
         if (num_ion_Zs < 3) then
            write(*, *) 'error: Get_ion_ZResults assumes num_ion_Zs >= 3'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         if (ion_Zs(1) /= 0) then
            write(*, *) 'error: Get_ion_ZResults assumes ion_Zs(1) == 0'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         if (abs(ion_Zs(1) - 2*ion_Zs(2) + ion_Zs(3)) > tiny) then
            write(*, *) 'error: Get_ion_ZResults assumes equal spaced Zs(1:3)'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         if (Z < tiny) then
            call Get_ion_for_X(1, X, Rho, logRho, T, logT, res, ierr)
            return
         end if
         
         if (Z > ion_Zs(3)) then
         
            if (Z <= ion_Zs(4)) then
               call do_interp2(3,4,ierr)
               if (ierr /= 0) return
            else if (Z <= ion_Zs(5) - tiny) then
               call do_interp2(4,5,ierr)
               if (ierr /= 0) return
            else
               call Get_ion_for_X(5, X, Rho, logRho, T, logT, res, ierr)
               return
            end if
            
         else
         
            do iz = 1, 3
               call Get_ion_for_X(iz, X, Rho, logRho, T, logT, res_zx(:, iz), ierr)
               if (ierr /= 0) return
            end do

            dZ = ion_Zs(2) - ion_Zs(1)
            denom = 2*dZ*dZ
            c(1) = (2*dZ*dZ - 3*dZ*Z + Z*Z)/denom
            c(2) = 2*(2*dZ-Z)*Z/denom
            c(3) = Z*(Z-dZ)/denom         

         end if
         
         res(:) = c(1)*res_zx(:, 1) + c(2)*res_zx(:, 2) + c(3)*res_zx(:, 3)
     
     
         contains
         
         subroutine do_interp2(iz1, iz2, ierr)
            use const_def, only: pi
            integer, intent(in) :: iz1, iz2
            integer, intent(out) :: ierr
            real(dp) :: Z1, Z2, alfa, beta
            ierr = 0
            Z1 = ion_Zs(iz1)
            Z2 = ion_Zs(iz2)
            alfa = (Z - Z1) / (Z2 - Z1)
            alfa = (1 - cospi(alfa))/2 ! smooth the transitions
            beta = 1 - alfa
            call Get_ion_for_X(iz1, X, Rho, logRho, T, logT, res_zx(:, 1), ierr)
            if (ierr /= 0) return
            call Get_ion_for_X(iz2, X, Rho, logRho, T, logT, res_zx(:, 2), ierr)
            if (ierr /= 0) return
            c(1) = beta
            c(2) = alfa
            c(3) = 0
            res_zx(:,3) = 0
         end subroutine do_interp2
     
      end subroutine Get_ion_ZResults

      
      subroutine Get_ion_for_X(iz, X, Rho, logRho, T, logT,res, ierr)
         integer, intent(in) :: iz ! the index in ion_Zs
         real(dp), intent(in) :: X, Rho, logRho, T, logT
         real(dp), intent(inout) :: res(num_ion_vals)
         integer, intent(out) :: ierr

         real(dp), dimension(num_ion_vals, 4) :: res_zx
         real(dp) :: dX, c(4), denom, delX, coef
         character (len=256) :: message
         integer :: ix, ix_lo, ix_hi, j, num_Xs
         real(dp), parameter :: tiny = 1d-6
         logical, parameter :: dbg_for_X = dbg ! .or. .true.

         ierr = 0
         
         if (num_ion_Xs /= 6) then
            write(*, *) 'error: Get_ion_for_X assumes num_ion_Xs == 6'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         if (ion_Xs(1) /= 0) then
            write(*, *) 'error: Get_ion_for_X assumes ion_Xs(1) == 0'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         num_Xs = num_ion_Xs_for_Z(iz)
         
         if (X < tiny .or. num_Xs == 1) then
            call Get_ion_XTable_Results(1, iz, Rho, logRho, T, logT, res, ierr)
            return
         end if
         
         dX = ion_Xs(2)-ion_Xs(1)
         
         do ix = 3, num_Xs
            if (abs(dX - (ion_Xs(ix) - ion_Xs(ix-1))) > tiny) then
               write(*, *) 'error: Get_ion_for_X assumes equal spaced Xs'
               call mesa_error(__FILE__,__LINE__)
            end if
         end do
         
         ix_hi = -1
         if (X <= ion_Xs(2)) then
            ix_lo = 1; ix_hi = 3
         else if (X >= ion_Xs(num_Xs-1)) then
            ix_lo = num_Xs-2; ix_hi = num_Xs
         else
            do ix = 3, num_Xs-1
               if (X <= ion_Xs(ix)) then
                  ix_lo = ix-2; ix_hi = ix+1; exit
               end if
            end do
         end if
         
         if (ix_hi < 0) then
            write(*, *) 'X', X
            write(*, *) 'ix_lo', ix_lo
            write(*, *) 'ix_hi', ix_hi
            write(*, *) 'error: Get_ion_for_X logic bug'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         if (dbg_for_X) then
            write(*, *) 'X', X
            write(*, *) 'ix_lo', ix_lo
            write(*, *) 'ix_hi', ix_hi
         end if
         
         do ix=ix_lo, ix_hi
            j = ix-ix_lo+1
            call Get_ion_XTable_Results(ix, iz, Rho, logRho, T, logT, &
               res_zx(:, j), ierr)
            if (ierr /= 0) return
         end do
         
         delX = X - ion_Xs(ix_lo)
         if (ix_hi-ix_lo==2) then
         
            denom = 2*dX*dX
            c(1) = (2*dX*dX - 3*dX*delX + delX*delX)/denom
            c(2) = 2*(2*dX-delX)*delX/denom
            c(3) = delX*(delX-dX)/denom
            res(:) = c(1)*res_zx(:, 1) + c(2)*res_zx(:, 2) + c(3)*res_zx(:, 3)

 1          format(a30, e25.15)         
            if (dbg_for_X) then
            end if
            
         else
         
            coef = (X-ion_Xs(ix_lo+1))/dX 
            ! coef = fractional location of X between 2nd and 3rd X's for fit.
            ! coef is weight for the quadratic based on points 2, 3, 4 of fit.
            ! (1-coef) is weight for quadratic based on points 1, 2, 3 of fit.
            if (coef < 0 .or. coef > 1) then
               write(*,*) 'logic bug in Get_ion_for_X'
               call mesa_error(__FILE__,__LINE__)
            end if
            c(1) = -coef*(coef-1)*(coef-1)/2
            c(2) = (2 + coef*coef*(-5 + 3*coef))/2
            c(3) = (coef + coef*coef*(4 - 3*coef))/2
            c(4) = coef*coef*(coef-1)/2
            res(:) = c(1)*res_zx(:, 1) + c(2)*res_zx(:, 2) &
                  + c(3)*res_zx(:, 3) + c(4)*res_zx(:, 4)
          
         end if
         
         
      end subroutine Get_ion_for_X
         
      
      subroutine Get_ion_XTable_Results(ix, iz, Rho, logRho_in, T, logT_in, &
               res, ierr)
         integer, intent(in) :: ix, iz
         real(dp), intent(in) :: Rho, logRho_in, T, logT_in
         real(dp), intent(inout) :: res(num_ion_vals)
         integer, intent(out) :: ierr
         
         real(dp) :: logQ0, logQ1, logT0, logT1
         integer :: iQ, jtemp
         
         real(dp) :: logRho, logT, logQ
         
         include 'formats'

         logRho = logRho_in
         logT = logT_in
         logQ = logRho - 2*logT + 12

         ierr = 0
      
         call Locate_logQ(logQ, iQ, logQ0, logQ1, ierr)
         if (ierr /= 0) return
      
         call Locate_logT(logT, jtemp, logT0, logT1, ierr)
         if (ierr /= 0) return
         
         call Do_ion_Interpolations(&
                  ion_num_logQs, ion_logQs, ion_num_logTs, ion_logTs, &
                  ion_tbl(:, :, :, :, ix, iz), &
                  iQ, jtemp, logQ0, logQ, logQ1, logT0, logT, logT1, &
                  res, ierr)
                  

      end subroutine Get_ion_XTable_Results
      
      
      subroutine Locate_logQ(logQ, iQ, logQ0, logQ1, ierr)
         real(dp), intent(inout) :: logQ
         integer, intent(out) :: iQ
         real(dp), intent(out) :: logQ0, logQ1
         integer, intent(out) :: ierr
      
         ierr = 0
         iQ = int((logQ - ion_logQ_min) / ion_del_logQ + 1d-4) + 1
         
         if (iQ < 1 .or. iQ >= ion_num_logQs) then
            
            if (iQ < 1) then
               iQ = 1
               logQ0 = ion_logQ_min
               logQ1 = logQ0 + ion_del_logQ
               logQ = logQ0
            else
               iQ = ion_num_logQs-1
               logQ0 = ion_logQ_min + (iQ-1) * ion_del_logQ
               logQ1 = logQ0 + ion_del_logQ
               logQ = logQ1
            end if
            
         else
         
            logQ0 = ion_logQ_min + (iQ-1) * ion_del_logQ
            logQ1 = logQ0 + ion_del_logQ

         end if

      end subroutine Locate_logQ
      
      
      subroutine Locate_logT(logT, iT, logT0, logT1, ierr)
         real(dp), intent(inout) :: logT
         integer, intent(out) :: iT
         real(dp), intent(out) :: logT0, logT1
         integer, intent(out) :: ierr
      
         ierr = 0
         iT = int((logT - ion_logT_min) / ion_del_logT + 1d-4) + 1
         
         if (iT < 1 .or. iT >= ion_num_logTs) then
            
            if (iT < 1) then
               iT = 1
               logT0 = ion_logT_min
               logT1 = logT0 + ion_del_logT
               logT = logT0
            else
               iT = ion_num_logTs-1
               logT0 = ion_logT_min + (iT-1) * ion_del_logT
               logT1 = logT0 + ion_del_logT
               logT = logT1
            end if
            
         else
         
            logT0 = ion_logT_min + (iT-1) * ion_del_logT
            logT1 = logT0 + ion_del_logT

         end if

      end subroutine Locate_logT
      
      
      subroutine Do_ion_Interpolations(&
                   nx, x, ny, y, fin, i, j, &
                   x0, xget, x1, y0, yget, y1, &
                   res, ierr)
!     >             fval, df_dx, df_dy, ierr)
         use interp_2d_lib_sg, only: interp_evbipm_sg
         integer, intent(in) :: nx, ny
         real(dp), intent(in) :: x(nx), y(ny), fin(4,num_ion_vals,nx,ny)
         integer, intent(in) :: i, j           ! target cell in f
         real(dp), intent(in) :: x0, xget, x1      ! x0 <= xget <= x1;  x0 = xs(i), x1 = xs(i+1)
         real(dp), intent(in) :: y0, yget, y1      ! y0 <= yget <= y1;  y0 = ys(j), y1 = ys(j+1)
         real(dp), intent(inout) :: res(num_ion_vals)
         integer, intent(out) :: ierr
   
         real(dp), parameter :: z36th = 1d0/36d0
         real(dp) :: xp, xpi, xp2, xpi2, ax, axbar, bx, bxbar, cx, cxi, hx2, cxd, cxdi, hx, hxi
         real(dp) :: yp, ypi, yp2, ypi2, ay, aybar, by, bybar, cy, cyi, hy2, cyd, cydi, hy, hyi
         real(dp) :: sixth_hx2, sixth_hy2, z36th_hx2_hy2
         real(dp) :: sixth_hx, sixth_hxi_hy2, z36th_hx_hy2
         real(dp) :: sixth_hx2_hyi, sixth_hy, z36th_hx2_hy
         integer :: k, ip1, jp1
         real(dp) :: f(4,nx,ny)
         
         include 'formats'
         
         ierr = 0
         hx=x1-x0
         hxi=1.0/hx
         hx2=hx*hx
   
         xp=(xget-x0)*hxi
         xpi=1.0-xp
         xp2=xp*xp
         xpi2=xpi*xpi

         ax=xp2*(3.0-2.0*xp)
         axbar=1.0-ax
         
         bx=-xp2*xpi
         bxbar=xpi2*xp
   
         cx=xp*(xp2-1.0)
         cxi=xpi*(xpi2-1.0)
         cxd=3.0*xp2-1.0
         cxdi=-3.0*xpi2+1.0
   
         hy=y1-y0
         hyi=1.0/hy
         hy2=hy*hy
   
         yp=(yget-y0)*hyi
         ypi=1.0-yp
         yp2=yp*yp
         ypi2=ypi*ypi

         ay=yp2*(3.0-2.0*yp)
         aybar=1.0-ay
         
         by=-yp2*ypi
         bybar=ypi2*yp
   
         cy=yp*(yp2-1.0)
         cyi=ypi*(ypi2-1.0)
         cyd=3.0*yp2-1.0
         cydi=-3.0*ypi2+1.0
                  
         sixth_hx2 = one_sixth*hx2
         sixth_hy2 = one_sixth*hy2
         z36th_hx2_hy2 = z36th*hx2*hy2
         
         sixth_hx = one_sixth*hx
         sixth_hxi_hy2 = one_sixth*hxi*hy2
         z36th_hx_hy2 = z36th*hx*hy2
         
         sixth_hx2_hyi = one_sixth*hx2*hyi
         sixth_hy = one_sixth*hy
         z36th_hx2_hy = z36th*hx2*hy
         
         ip1 = i+1
         jp1 = j+1
         
         do k=1,num_ion_vals
            ! bicubic spline interpolation
            res(k) = dble(&
                  xpi*(&
                     ypi*fin(1,k,i,j)  +yp*fin(1,k,i,jp1))&
                     +xp*(ypi*fin(1,k,ip1,j)+yp*fin(1,k,ip1,jp1))&
                  +sixth_hx2*(&
                     cxi*(ypi*fin(2,k,i,j) +yp*fin(2,k,i,jp1))+&
                     cx*(ypi*fin(2,k,ip1,j)+yp*fin(2,k,ip1,jp1)))&
                  +sixth_hy2*(&
                     xpi*(cyi*fin(3,k,i,j) +cy*fin(3,k,i,jp1))+&
                     xp*(cyi*fin(3,k,ip1,j)+cy*fin(3,k,ip1,jp1)))&
                  +z36th_hx2_hy2*(&
                     cxi*(cyi*fin(4,k,i,j) +cy*fin(4,k,i,jp1))+&
                     cx*(cyi*fin(4,k,ip1,j)+cy*fin(4,k,ip1,jp1))))

            ! derivatives of bicubic splines
!            df_dx(k) =
!     >            hxi*(
!     >               -(ypi*fin(1,k,i,j)  +yp*fin(1,k,i,jp1))
!     >               +(ypi*fin(1,k,ip1,j)+yp*fin(1,k,ip1,jp1)))
!     >            +sixth_hx*(
!     >               cxdi*(ypi*fin(2,k,i,j) +yp*fin(2,k,i,jp1))+
!     >               cxd*(ypi*fin(2,k,ip1,j)+yp*fin(2,k,ip1,jp1)))
!     >            +sixth_hxi_hy2*(
!     >               -(cyi*fin(3,k,i,j)  +cy*fin(3,k,i,jp1))
!     >               +(cyi*fin(3,k,ip1,j)+cy*fin(3,k,ip1,jp1)))
!     >            +z36th_hx_hy2*(
!     >               cxdi*(cyi*fin(4,k,i,j) +cy*fin(4,k,i,jp1))+
!     >               cxd*(cyi*fin(4,k,ip1,j)+cy*fin(4,k,ip1,jp1)))
!
!            df_dy(k) =
!     >            hyi*(
!     >               xpi*(-fin(1,k,i,j) +fin(1,k,i,jp1))+
!     >               xp*(-fin(1,k,ip1,j)+fin(1,k,ip1,jp1)))
!     >            +sixth_hx2_hyi*(
!     >               cxi*(-fin(2,k,i,j) +fin(2,k,i,jp1))+
!     >               cx*(-fin(2,k,ip1,j)+fin(2,k,ip1,jp1)))
!     >            +sixth_hy*(
!     >               xpi*(cydi*fin(3,k,i,j) +cyd*fin(3,k,i,jp1))+
!     >               xp*(cydi*fin(3,k,ip1,j)+cyd*fin(3,k,ip1,jp1)))
!     >            +z36th_hx2_hy*(
!     >               cxi*(cydi*fin(4,k,i,j) +cyd*fin(4,k,i,jp1))+
!     >               cx*(cydi*fin(4,k,ip1,j)+cyd*fin(4,k,ip1,jp1)))
     
         end do
         
      end subroutine Do_ion_Interpolations


      end module ion_tables_eval
      
