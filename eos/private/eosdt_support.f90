! ***********************************************************************
!
!   Copyright (C) 2017-2019  Bill Paxton & The MESA Team
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

      module eosdt_support
      use eos_def
      use const_def, only: avo, crad, ln10, arg_not_provided, mp, kerg, dp, qp, one_sixth
      use utils_lib, only: is_bad, mesa_error
      use math_lib
      
      implicit none
         
      integer, parameter :: sz = sz_per_eos_point
      
      contains

      
      subroutine Do_EoS_Interpolations( &
             nvlo, nvhi, n, nx, x, ny, y, fin1, i, j, &
             x0, xget, x1, y0, yget, y1, &
             fval, df_dx, df_dy, ierr)
         integer, intent(in) :: nvlo, nvhi, n, nx, ny
         real(dp), intent(in) :: x(:) ! (nx)
         real(dp), intent(in) :: y(:) ! (ny)
         real(dp), intent(in), pointer :: fin1(:) ! =(4,n,nx,ny)
         integer, intent(in) :: i, j           ! target cell in f
         real(dp), intent(in) :: x0, xget, x1      ! x0 <= xget <= x1;  x0 = xs(i), x1 = xs(i+1)
         real(dp), intent(in) :: y0, yget, y1      ! y0 <= yget <= y1;  y0 = ys(j), y1 = ys(j+1)
         real(dp), intent(inout), dimension(nv) :: fval, df_dx, df_dy
         integer, intent(out) :: ierr
   
         real(dp) :: xp, xpi, xp2, xpi2, ax, axbar, bx, bxbar, cx, cxi, hx2, cxd, cxdi, hx, hxi
         real(dp) :: yp, ypi, yp2, ypi2, ay, aybar, by, bybar, cy, cyi, hy2, cyd, cydi, hy, hyi
         real(dp) :: sixth_hx2, sixth_hy2, z36th_hx2_hy2
         real(dp) :: sixth_hx, sixth_hxi_hy2, z36th_hx_hy2
         real(dp) :: sixth_hx2_hyi, sixth_hy, z36th_hx2_hy
         integer :: k, ip1, jp1
         real(dp), pointer :: fin(:,:,:,:)
         
         include 'formats'
         
         ierr = 0
         
         fin(1:sz_per_eos_point,1:n,1:nx,1:ny) => &
            fin1(1:sz_per_eos_point*n*nx*ny)
         
         hx=x1-x0
         hxi=1d0/hx
         hx2=hx*hx
   
         xp=(xget-x0)*hxi

         xpi=1d0-xp
         xp2=xp*xp
         xpi2=xpi*xpi

         ax=xp2*(3d0-2d0*xp)
         axbar=1d0-ax
         
         bx=-xp2*xpi
         bxbar=xpi2*xp
   
         cx=xp*(xp2-1d0)
         cxi=xpi*(xpi2-1d0)
         cxd=3d0*xp2-1d0
         cxdi=-3d0*xpi2+1d0
   
         hy=y1-y0
         hyi=1d0/hy
         hy2=hy*hy
   
         yp=(yget-y0)*hyi
         
         ypi=1d0-yp
         yp2=yp*yp
         ypi2=ypi*ypi

         ay=yp2*(3d0-2d0*yp)
         aybar=1d0-ay
         
         by=-yp2*ypi
         bybar=ypi2*yp
   
         cy=yp*(yp2-1d0)
         cyi=ypi*(ypi2-1d0)
         cyd=3d0*yp2-1d0
         cydi=-3d0*ypi2+1d0
                  
         sixth_hx2 = one_sixth*hx2
         sixth_hy2 = one_sixth*hy2
         z36th_hx2_hy2 = sixth_hx2*sixth_hy2
         
         sixth_hx = one_sixth*hx
         sixth_hxi_hy2 = sixth_hy2*hxi
         z36th_hx_hy2 = sixth_hx*sixth_hy2
         
         sixth_hx2_hyi = sixth_hx2*hyi
         sixth_hy = one_sixth*hy
         z36th_hx2_hy = sixth_hx2*sixth_hy
         
         ip1 = i+1
         jp1 = j+1
         
         !$omp simd
         do k = nvlo, nvhi
            ! bicubic spline interpolation
            
            ! f(1,i,j) = f(x(i),y(j))
            ! f(2,i,j) = d2f/dx2(x(i),y(j))
            ! f(3,i,j) = d2f/dy2(x(i),y(j))
            ! f(4,i,j) = d4f/dx2dy2(x(i),y(j))

            fval(k) = &
                  xpi*( &
                     ypi*fin(1,k,i,j)  +yp*fin(1,k,i,jp1)) &
                     +xp*(ypi*fin(1,k,ip1,j)+yp*fin(1,k,ip1,jp1)) &
                  +sixth_hx2*( &
                     cxi*(ypi*fin(2,k,i,j) +yp*fin(2,k,i,jp1))+ &
                     cx*(ypi*fin(2,k,ip1,j)+yp*fin(2,k,ip1,jp1))) &
                  +sixth_hy2*( &
                     xpi*(cyi*fin(3,k,i,j) +cy*fin(3,k,i,jp1))+ &
                     xp*(cyi*fin(3,k,ip1,j)+cy*fin(3,k,ip1,jp1))) &
                  +z36th_hx2_hy2*( &
                     cxi*(cyi*fin(4,k,i,j) +cy*fin(4,k,i,jp1))+ &
                     cx*(cyi*fin(4,k,ip1,j)+cy*fin(4,k,ip1,jp1)))
            
            ! derivatives of bicubic splines
            df_dx(k) = &
                  hxi*( &
                     -(ypi*fin(1,k,i,j)  +yp*fin(1,k,i,jp1)) &
                     +(ypi*fin(1,k,ip1,j)+yp*fin(1,k,ip1,jp1))) &
                  +sixth_hx*( &
                     cxdi*(ypi*fin(2,k,i,j) +yp*fin(2,k,i,jp1))+ &
                     cxd*(ypi*fin(2,k,ip1,j)+yp*fin(2,k,ip1,jp1))) &
                  +sixth_hxi_hy2*( &
                     -(cyi*fin(3,k,i,j)  +cy*fin(3,k,i,jp1)) &
                     +(cyi*fin(3,k,ip1,j)+cy*fin(3,k,ip1,jp1))) &
                  +z36th_hx_hy2*( &
                     cxdi*(cyi*fin(4,k,i,j) +cy*fin(4,k,i,jp1))+ &
                     cxd*(cyi*fin(4,k,ip1,j)+cy*fin(4,k,ip1,jp1)))

            df_dy(k) = &
                  hyi*( &
                     xpi*(-fin(1,k,i,j) +fin(1,k,i,jp1))+ &
                     xp*(-fin(1,k,ip1,j)+fin(1,k,ip1,jp1))) &
                  +sixth_hx2_hyi*( &
                     cxi*(-fin(2,k,i,j) +fin(2,k,i,jp1))+ &
                     cx*(-fin(2,k,ip1,j)+fin(2,k,ip1,jp1))) &
                  +sixth_hy*( &
                     xpi*(cydi*fin(3,k,i,j) +cyd*fin(3,k,i,jp1))+ &
                     xp*(cydi*fin(3,k,ip1,j)+cyd*fin(3,k,ip1,jp1))) &
                  +z36th_hx2_hy*( &
                     cxi*(cydi*fin(4,k,i,j) +cyd*fin(4,k,i,jp1))+ &
                     cx*(cydi*fin(4,k,ip1,j)+cyd*fin(4,k,ip1,jp1)))
         
         end do
         
      end subroutine Do_EoS_Interpolations


      subroutine Do_Blend( &
            rq, species, Rho, logRho, T, logT, &
            alfa_in, d_alfa_dlogT_in, d_alfa_dlogRho_in, linear_blend, &
            res_1, d_dlnd_1, d_dlnT_1, d_dxa_1, &
            res_2, d_dlnd_2, d_dlnT_2, d_dxa_2, &
            res, dlnd, dlnT, d_dxa)
         type (EoS_General_Info), pointer :: rq
         integer, intent(in) :: species
         real(dp), intent(in) :: Rho, logRho, T, logT, &
            alfa_in, d_alfa_dlogT_in, d_alfa_dlogRho_in
         logical, intent(in) :: linear_blend
         real(dp), intent(in), dimension(nv) :: res_1, d_dlnd_1, d_dlnT_1, res_2, d_dlnd_2, d_dlnT_2
         real(dp), intent(in), dimension(nv, species) :: d_dxa_1, d_dxa_2
         real(dp), intent(out), dimension(nv) :: res, dlnd, dlnT
         real(dp), intent(out), dimension(nv, species) :: d_dxa

         real(dp) :: alfa0, d_alfa0_dlnT, d_alfa0_dlnd
         real(dp) :: alfa, beta, A, &
            d_alfa_dlnd, d_alfa_dlnT, &
            d_beta_dlnd, d_beta_dlnT
         integer :: j, k
         
         if (.not. linear_blend) then 
         
            ! smooth the transitions near alfa = 0.0 and 1.0
            ! quintic smoothing function with 1st and 2nd derivs = 0 at ends

            alfa0 = alfa_in
            d_alfa0_dlnT = d_alfa_dlogT_in/ln10
            d_alfa0_dlnd = d_alfa_dlogRho_in/ln10
            alfa = -alfa0*alfa0*alfa0*(-10d0 + alfa0*(15d0 - 6d0*alfa0))
            A = 30d0*(alfa0 - 1d0)*(alfa0 - 1d0)*alfa0*alfa0            
            d_alfa_dlnd = A*d_alfa0_dlnd
            d_alfa_dlnT = A*d_alfa0_dlnT
            
         else
            
            alfa = alfa_in
            d_alfa_dlnT = d_alfa_dlogT_in/ln10
            d_alfa_dlnd = d_alfa_dlogRho_in/ln10
         
         end if            
   
         beta = 1d0 - alfa
         d_beta_dlnT = -d_alfa_dlnT
         d_beta_dlnd = -d_alfa_dlnd

         do j=1,nv
            res(j) = alfa*res_1(j) + beta*res_2(j) 
            dlnd(j) = &
               alfa*d_dlnd_1(j) + beta*d_dlnd_2(j) + &
               d_alfa_dlnd*res_1(j) + d_beta_dlnd*res_2(j)
            dlnT(j) = &
               alfa*d_dlnT_1(j) + beta*d_dlnT_2(j) + &
               d_alfa_dlnT*res_1(j) + d_beta_dlnT*res_2(j)
         end do

         do k = 1, species
            do j = 1, nv
               d_dxa(j,k) = &
                  alfa*d_dxa_1(j,k) + beta*d_dxa_2(j,k) ! ignores blend derivatives
            end do
         end do

      end subroutine Do_Blend


      end module eosdt_support
      
