! ***********************************************************************
!
!   Copyright (C) 2010-2024  The MESA Team
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

module parasite_model_matrices

   use const_def, only: dp

   implicit none

   private
   public :: build_parasite_matrix
   public :: build_parasite_matrix_LPN
   public :: build_parasite_matrix_LPN_QS

contains

   subroutine build_parasite_matrix(w, k_z, Pr, tau, R_0, H_B, D_B, lam_hat, l_hat, N, L, parity)

      real(dp), intent(in)               :: w
      real(dp), intent(in)               :: k_z
      real(dp), intent(in)               :: Pr
      real(dp), intent(in)               :: tau
      real(dp), intent(in)               :: R_0
      real(dp), intent(in)               :: H_B
      real(dp), intent(in)               :: D_B
      real(dp), intent(in)               :: lam_hat
      real(dp), intent(in)               :: l_hat
      integer, intent(in)                :: N
      real(dp), allocatable, intent(out) :: L(:,:)
      character(*), intent(in), optional :: parity

      real(dp) :: E_psi
      real(dp) :: E_T
      real(dp) :: E_C
      logical  :: parity_switch
      integer  :: i
      integer  :: m
      real(dp) :: B(4,4)

      ! Build the parasite model matrix

      E_psi = w/(2*l_hat)
      E_T   = -l_hat*E_psi/(lam_hat + l_hat**2)
      E_C   = -l_hat*E_psi/(R_0*(lam_hat + tau*l_hat**2))

      if (PRESENT(parity)) then

         ! Parity-split case

         select case (parity)
         case ('EVEN')
            parity_switch = .TRUE.
         case ('ODD')
            parity_switch = .FALSE.
         case default
            stop '** invalid parity in build_parasite_matrix'
         end select

         allocate(L(4*N+2,4*N+2))
         L = 0._dp

         i = 1

         split_block_loop : do m = 0, N

            call eval_4x4_block_m_(m, B)

            if (m == 0) then
               if (parity_switch) then
                  L(i:i+1,i:i+1) = B([1,4],[1,4])
               else
                  L(i:i+1,i:i+1) = B([2,3],[2,3])
               end if
            else
               L(i:i+3,i:i+3) = B
            end if

            if (m > 0) then

               call eval_4x4_block_mm1_(m, B)

               if (m == 1) then
                  if (parity_switch) then
                     L(i:i+3,i-2:i-1) = B(:,[1,4])
                  else
                     L(i:i+3,i-2:i-1) = B(:,[2,3])
                  end if
               else
                  L(i:i+3,i-4:i-1) = B
               end if

            end if

            if (m < N) then

               call eval_4x4_block_mp1_(m, B)

               if (m == 0) then
                  if (parity_switch) then
                     L(i:i+1,i+2:i+5) = 2*B([1,4],:)
                  else
                     L(i:i+1,i+2:i+5) = 2*B([2,3],:)
                  end if
               else
                  L(i:i+3,i+4:i+7) = B
               end if
            end if

            if (m == 0) then
               i = i + 2
            else
               i = i + 4
            end if

         end do split_block_loop

      else

         allocate(L(8*N+4,8*N+4))
         L = 0._dp

         i = 1

         full_block_loop : do m = -N, N

            call eval_4x4_block_m_(m, L(i:i+3,i:i+3))
 
            if (m > -N) then
               call eval_4x4_block_mm1_(m, L(i:i+3,i-4:i-1))
            endif
            
            if (m < N) then
               call eval_4x4_block_mp1_(m, L(i:i+3,i+4:i+7))
            end if
            
            i = i + 4

         end do full_block_loop

      end if

   contains

      subroutine eval_4x4_block_m_(m, B)

         integer, intent(in)   :: m
         real(dp), intent(out) :: B(:,:)

         real(dp) :: k2_m

         ! Evaluate the (m,m) 4x4 block

         k2_m = m**2*l_hat**2 + k_z**2

         if (m /= 0) then

            B(1,1) = -Pr*k2_m
            B(1,2) =  Pr*m*l_hat/k2_m
            B(1,3) = -Pr*m*l_hat/k2_m
            B(1,4) = -H_B*k_z

            B(2,1) = -m*l_hat
            B(2,2) = -k2_m
            B(2,3) = 0
            B(2,4) = 0

            B(3,1) = -m*l_hat/R_0
            B(3,2) = 0
            B(3,3) = -tau*k2_m
            B(3,4) = 0

            B(4,1) = k_z
            B(4,2) = 0
            B(4,3) = 0
            B(4,4) = -D_B*k2_m

         else

            B(1,1) = -Pr*k2_m
            B(1,2) = 0
            B(1,3) = 0
            B(1,4) = -H_B*k_z

            B(2,1) = 0
            B(2,2) = -k2_m
            B(2,3) = 0
            B(2,4) = 0

            B(3,1) = 0
            B(3,2) = 0
            B(3,3) = -tau*k2_m
            B(3,4) = 0
            
            B(4,1) = k_z
            B(4,2) = 0
            B(4,3) = 0
            B(4,4) = -D_B*k2_m

         end if

      end subroutine eval_4x4_block_m_

      !****
      
      subroutine eval_4x4_block_mm1_(m, B)

         integer, intent(in)   :: m
         real(dp), intent(out) :: B(:,:)

         real(dp) :: k2_m
         real(dp) :: k2_m_m

         ! Evaluate the (m,m-1) 4x4 block

         k2_m = m**2*l_hat**2 + k_z**2
         k2_m_m = (m-1)**2*l_hat**2 + k_z**2

         if (E_psi /= 0._dp) then

            B(1,1) = l_hat*k_z*E_psi/k2_m*(l_hat**2 - k2_m_m)
            B(1,2) = 0
            B(1,3) = 0
            B(1,4) = 0

            B(2,1) = l_hat*k_z*E_T
            B(2,2) = -l_hat*k_z*E_psi
            B(2,3) = 0
            B(2,4) = 0

            B(3,1) = l_hat*k_z*E_C
            B(3,2) = 0
            B(3,3) = -l_hat*k_z*E_psi
            B(3,4) = 0

            B(4,1) = 0
            B(4,2) = 0
            B(4,3) = 0
            B(4,4) = -l_hat*k_z*E_psi

         else

            B = 0

         end if

      end subroutine eval_4x4_block_mm1_

      !****

      subroutine eval_4x4_block_mp1_(m, B)

         integer, intent(in)   :: m
         real(dp), intent(out) :: B(:,:)

         real(dp) :: k2_m
         real(dp) :: k2_m_p
      
         ! Evaluate the (m,m+1) 4x4 matrix block

         k2_m = m**2*l_hat**2 + k_z**2
         k2_m_p = (m+1)**2*l_hat**2 + k_z**2

         if (E_psi /= 0._dp) then

            B(1,1) = -l_hat*k_z*E_psi/k2_m*(l_hat**2 - k2_m_p)
            B(1,2) = 0
            B(1,3) = 0
            B(1,4) = 0

            B(2,1) = l_hat*k_z*E_T
            B(2,2) = l_hat*k_z*E_psi
            B(2,3) = 0
            B(2,4) = 0

            B(3,1) = l_hat*k_z*E_C
            B(3,2) = 0
            B(3,3) = l_hat*k_z*E_psi
            B(3,4) = 0

            B(4,1) = 0
            B(4,2) = 0
            B(4,3) = 0
            B(4,4) = l_hat*k_z*E_psi

         else

            B = 0

         end if

      end subroutine eval_4x4_block_mp1_

   end subroutine build_parasite_matrix

   !****

   subroutine build_parasite_matrix_LPN(w, k_z, Pr, tau, R_0, H_B, D_B, lam_hat, l_hat, N, L, parity)

      real(dp), intent(in)               :: w
      real(dp), intent(in)               :: k_z
      real(dp), intent(in)               :: Pr
      real(dp), intent(in)               :: tau
      real(dp), intent(in)               :: R_0
      real(dp), intent(in)               :: H_B
      real(dp), intent(in)               :: D_B
      real(dp), intent(in)               :: lam_hat
      real(dp), intent(in)               :: l_hat
      integer, intent(in)                :: N
      real(dp), allocatable, intent(out) :: L(:,:)
      character(*), intent(in), optional :: parity

      real(dp) :: E_psi
      real(dp) :: E_C
      logical  :: parity_switch
      integer  :: i
      integer  :: m
      real(dp) :: B(3,3)

      if (.NOT. PRESENT(parity)) then
         stop 'unsplit option not implemented for build_parasite_matrix_LPN'
      end if

      ! Build the parasite model matrix in the low Peclet number (LPN)
      ! limit tau << 1

      E_psi = w/(2*l_hat)
      E_C   = -l_hat*E_psi/(R_0*(lam_hat + tau*l_hat**2))

      select case (parity)
      case ('EVEN')
         parity_switch = .TRUE.
      case ('ODD')
         parity_switch = .FALSE.
      case default
         stop '** invalid parity in build_parasite_matrix_LPN'
      end select

      if (parity_switch) then
         allocate(L(3*N+2,3*N+2))
      else
         allocate(L(3*N+1,3*N+1))
      end if

      L = 0._dp

      i = 1

      block_loop : do m = 0, N

         call eval_3x3_block_m_(m, B)

         if (m == 0) then
            if (parity_switch) then
               L(i:i+1,i:i+1) = B([1,3],[1,3])
            else
               L(i:i,i:i) = B([2],[2])
            end if
         else
            L(i:i+2,i:i+2) = B
         end if

         if (m > 0) then

            call eval_3x3_block_mm1_(m, B)

            if (m == 1) then
               if (parity_switch) then
                  L(i:i+2,i-2:i-1) = B(:,[1,3])
               else
                  L(i:i+2,i-1:i-1) = B(:,[2])
               end if
            else
               L(i:i+2,i-3:i-1) = B
            end if
            
         end if

         if (m < N) then

            call eval_3x3_block_mp1_(m, B)

            if (m == 0) then
               if (parity_switch) then
                  L(i:i+1,i+2:i+4) = 2*B([1,3],:)
               else
                  L(i:i,i+2:i+4) = 2*B([2],:)
               end if
            else
               L(i:i+2,i+3:i+5) = B
            end if
         end if

         if (m == 0) then
            if (parity_switch) then
               i = i + 2
            else
               i = i + 1
            end if
         else
            i = i + 3
         end if
         
      end do block_loop

   contains

      subroutine eval_3x3_block_m_(m, B)

         integer, intent(in)   :: m
         real(dp), intent(out) :: B(:,:)

         real(dp) :: k2_m

         ! Evaluate the (m,m) 3x3 block

         k2_m = m**2*l_hat**2 + k_z**2

         if (m /= 0) then

            B(1,1) = -Pr*(k2_m + m**2*l_hat**2/k2_m**2)
            B(1,2) = -Pr*m*l_hat/k2_m
            B(1,3) = -H_B*k_z

            B(2,1) = -m*l_hat/R_0
            B(2,2) = -tau*k2_m
            B(2,3) = 0

            B(3,1) = k_z
            B(3,2) = 0
            B(3,4) = -D_B*k2_m

         else

            B(1,1) = -Pr*k2_m
            B(1,2) = 0
            B(1,3) = -H_B*k_z

            B(2,1) = 0
            B(2,2) = -tau*k2_m
            B(2,3) = 0
            
            B(3,1) = k_z
            B(3,2) = 0
            B(3,3) = -D_B*k2_m

         end if

      end subroutine eval_3x3_block_m_

      !****
      
      subroutine eval_3x3_block_mm1_(m, B)

         integer, intent(in)   :: m
         real(dp), intent(out) :: B(:,:)

         real(dp) :: k2_m
         real(dp) :: k2_m_m

         ! Evaluate the (m,m-1) 3x3 block

         k2_m = m**2*l_hat**2 + k_z**2
         k2_m_m = (m-1)**2*l_hat**2 + k_z**2

         if (E_psi /= 0._dp) then

            B(1,1) = l_hat*k_z*E_psi/k2_m*(l_hat**2 - k2_m_m)
            B(1,2) = 0
            B(1,3) = 0

            B(2,1) = l_hat*k_z*E_C
            B(2,2) = -l_hat*k_z*E_psi
            B(2,3) = 0

            B(3,1) = 0
            B(3,2) = 0
            B(3,3) = -l_hat*k_z*E_psi

         else

            B = 0

         end if

      end subroutine eval_3x3_block_mm1_

      !****

      subroutine eval_3x3_block_mp1_(m, B)

         integer, intent(in)   :: m
         real(dp), intent(out) :: B(:,:)

         real(dp) :: k2_m
         real(dp) :: k2_m_p
      
         ! Evaluate the (m,m+1) 3x3 matrix block

         k2_m = m**2*l_hat**2 + k_z**2
         k2_m_p = (m+1)**2*l_hat**2 + k_z**2

         if (E_psi /= 0._dp) then

            B(1,1) = -l_hat*k_z*E_psi/k2_m*(l_hat**2 - k2_m_p)
            B(1,2) = 0
            B(1,3) = 0

            B(2,1) = l_hat*k_z*E_C
            B(2,2) = l_hat*k_z*E_psi
            B(2,3) = 0

            B(3,1) = 0
            B(3,2) = 0
            B(3,3) = l_hat*k_z*E_psi

         else

            B = 0

         end if

      end subroutine eval_3x3_block_mp1_

   end subroutine build_parasite_matrix_LPN

   !****

   subroutine build_parasite_matrix_LPN_QS(w, k_z, Pr, tau, R_0, H_B, D_B, lam_hat, l_hat, N, L, parity)

      real(dp), intent(in)               :: w
      real(dp), intent(in)               :: k_z
      real(dp), intent(in)               :: Pr
      real(dp), intent(in)               :: tau
      real(dp), intent(in)               :: R_0
      real(dp), intent(in)               :: H_B
      real(dp), intent(in)               :: D_B
      real(dp), intent(in)               :: lam_hat
      real(dp), intent(in)               :: l_hat
      integer, intent(in)                :: N
      real(dp), allocatable, intent(out) :: L(:,:)
      character(*), intent(in), optional :: parity

      real(dp) :: E_psi
      real(dp) :: E_C
      logical  :: parity_switch
      integer  :: i
      integer  :: m
      real(dp) :: B(2,2)

      if (.NOT. PRESENT(parity)) then
         stop 'unsplit option not implemented for build_parasite_matrix_LPN_QS'
      end if

      ! Build the parasite model matrix in the low Peclet number (LPN)
      ! limit tau << 1 and quasi-static limit Rm << 1

      E_psi = w/(2*l_hat)
      E_C   = -l_hat*E_psi/(R_0*(lam_hat + tau*l_hat**2))

      select case (parity)
      case ('EVEN')
         parity_switch = .TRUE.
      case ('ODD')
         parity_switch = .FALSE.
      case default
         stop '** invalid parity in build_parasite_matrix_LPN_QS'
      end select

      allocate(L(2*N+2,2*N+2))

      L = 0._dp

      i = 1

      block_loop : do m = 0, N

         call eval_2x2_block_m_(m, B)

         if (m == 0) then
            if (parity_switch) then
               L(i:i,i:i) = B([1],[1])
            else
               L(i:i,i:i) = B([2],[2])
            end if
         else
            L(i:i+1,i:i+1) = B
         end if

         if (m > 0) then

            call eval_2x2_block_mm1_(m, B)

            if (m == 1) then
               if (parity_switch) then
                  L(i:i+1,i-1:i-1) = B(:,[1])
               else
                  L(i:i+1,i-1:i-1) = B(:,[2])
               end if
            else
               L(i:i+1,i-2:i-1) = B
            end if
            
         end if

         if (m < N) then

            call eval_2x2_block_mp1_(m, B)

            if (m == 0) then
               if (parity_switch) then
                  L(i:i,i+2:i+3) = 2*B([1],:)
               else
                  L(i:i,i+2:i+3) = 2*B([2],:)
               end if
            else
               L(i:i+1,i+2:i+3) = B
            end if
         end if

         if (m == 0) then
            i = i + 1
         else
            i = i + 2
         end if
         
      end do block_loop

   contains

      subroutine eval_2x2_block_m_(m, B)

         integer, intent(in)   :: m
         real(dp), intent(out) :: B(:,:)

         real(dp) :: k2_m

         ! Evaluate the (m,m) 2x2 block

         k2_m = m**2*l_hat**2 + k_z**2

         if (m /= 0) then

            B(1,1) = -Pr*(k2_m + m**2*l_hat**2/k2_m**2) - H_B*k_z**2/(D_B*k2_m)
            B(1,2) = -Pr*m*l_hat/k2_m

            B(2,1) = -m*l_hat/R_0
            B(2,2) = -tau*k2_m
      
         else

            B(1,1) = -Pr*k2_m - H_B*k_z**2/(D_B*k2_m)
            B(1,2) = 0

            B(2,1) = 0
            B(2,2) = -tau*k2_m

         end if

      end subroutine eval_2x2_block_m_

      !****
      
      subroutine eval_2x2_block_mm1_(m, B)

         integer, intent(in)   :: m
         real(dp), intent(out) :: B(:,:)

         real(dp) :: k2_m
         real(dp) :: k2_m_m

         ! Evaluate the (m,m-1) 2x2 block

         k2_m = m**2*l_hat**2 + k_z**2
         k2_m_m = (m-1)**2*l_hat**2 + k_z**2

         if (E_psi /= 0._dp) then

            B(1,1) = l_hat*k_z*E_psi/k2_m*(l_hat**2 - k2_m_m)
            B(1,2) = 0

            B(2,1) = l_hat*k_z*E_C
            B(2,2) = -l_hat*k_z*E_psi

         else

            B = 0

         end if

      end subroutine eval_2x2_block_mm1_

      !****

      subroutine eval_2x2_block_mp1_(m, B)

         integer, intent(in)   :: m
         real(dp), intent(out) :: B(:,:)

         real(dp) :: k2_m
         real(dp) :: k2_m_p
      
         ! Evaluate the (m,m+1) 2x2 matrix block

         k2_m = m**2*l_hat**2 + k_z**2
         k2_m_p = (m+1)**2*l_hat**2 + k_z**2

         if (E_psi /= 0._dp) then

            B(1,1) = -l_hat*k_z*E_psi/k2_m*(l_hat**2 - k2_m_p)
            B(1,2) = 0

            B(2,1) = l_hat*k_z*E_C
            B(2,2) = l_hat*k_z*E_psi

         else

            B = 0

         end if

      end subroutine eval_2x2_block_mp1_

   end subroutine build_parasite_matrix_LPN_QS

end module parasite_model_matrices
