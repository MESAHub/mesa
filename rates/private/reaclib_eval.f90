! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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
 
      module reaclib_eval
      use rates_def
      use math_lib
      
      implicit none

      real(dp), parameter :: &
         lam_max = 1d99, &
         ln1_max = 227.955924206411d0 ! log(lam_max)


      contains
      
      
      subroutine reaction_rates( &
            num_lambdas, lo, hi, T9, rates, nuclides, forward_only, &
            lambda, dlambda_dlnT, &
            rlambda, drlambda_dlnT, &
            ierr)
         integer, intent(in) :: num_lambdas, lo, hi
         real(dp), intent(in) :: T9
         type(reaction_data), intent(in) :: rates
         type(nuclide_data), intent(in) :: nuclides
         logical, intent(in) :: forward_only
         real(dp), intent(out) :: lambda, dlambda_dlnT
         real(dp), intent(out) :: rlambda, drlambda_dlnT
         integer, intent(out) :: ierr
         
         real(dp), dimension(num_lambdas), target :: &
            ln_lambdas_ar, lambdas_ar, dlambdas_dlnT_ar, &
            ln_rlambdas_ar, rlambdas_ar, drlambdas_dlnT_ar
         real(dp), dimension(:), pointer :: &
            ln_lambdas, lambdas, dlambdas_dlnT, &
            ln_rlambdas, rlambdas, drlambdas_dlnT
         
         include 'formats'
         
         ierr = 0
         
         ln_lambdas => ln_lambdas_ar
         lambdas => lambdas_ar
         dlambdas_dlnT => dlambdas_dlnT_ar
         ln_rlambdas => ln_rlambdas_ar
         rlambdas => rlambdas_ar
         drlambdas_dlnT => drlambdas_dlnT_ar
         
         call compute_some_lambdas( &
            num_lambdas, lo, hi, T9, rates, ln_lambdas, lambdas, dlambdas_dlnT)
         lambda = sum(lambdas(1:num_lambdas))
         dlambda_dlnT = sum(dlambdas_dlnT(1:num_lambdas))
         
         if (forward_only) then
            rlambda = 0
            drlambda_dlnT = 0
         else
            call compute_some_inverse_lambdas( &
               num_lambdas, lo, hi, T9, rates, &
               ln_lambdas, lambdas, dlambdas_dlnT, &
               rlambdas, drlambdas_dlnT)
            rlambda = sum(rlambdas(1:num_lambdas))
            drlambda_dlnT = sum(drlambdas_dlnT(1:num_lambdas))
         end if

      end subroutine reaction_rates
      

      subroutine compute_some_lambdas( &
            num_lambdas, lo, hi, T9, rates, ln_lambda, lambda, dlambda_dlnT)
         use utils_lib, only: is_bad
         use const_def
         integer, intent(in) :: num_lambdas, lo, hi ! range of rates to do
         real(dp), intent(in) :: T9
         type(reaction_data), intent(in) :: rates
         real(dp), dimension(:), intent(out) :: ln_lambda, lambda, dlambda_dlnT
         
         real(dp) :: T9inv, logT, ln1
         real(dp), dimension(7) :: T9fac, dT9fac_dT9, dT9fac_dlnT
         integer :: i, j
         
         include 'formats'

         T9inv = 1d0/T9
         
         T9fac(1) = 1d0
         dT9fac_dT9(1) = 0d0
         
         T9fac(2) = T9inv
         dT9fac_dT9(2) = -T9inv*T9inv
         
         T9fac(3) = pow(T9inv,one_third)
         dT9fac_dT9(3) = -one_third*pow(T9inv,four_thirds)
         
         T9fac(4) = pow(T9,one_third)
         dT9fac_dT9(4) = one_third*pow(T9inv,two_thirds)
         
         T9fac(5) = T9
         dT9fac_dT9(5) = 1d0
         
         T9fac(6) = pow(T9,five_thirds)
         dT9fac_dT9(6) = five_thirds*pow(T9,two_thirds)
         
         T9fac(7) = log(T9)
         dT9fac_dT9(7) = T9inv
         
         dT9fac_dlnT = T9*dT9fac_dT9
         
         do i = lo, hi
            j = i+1-lo
            ln1 = dot_product(T9fac(:), rates% coefficients(:,i))
            if (ln1 > ln1_max) then
               ln_lambda(j) = ln1_max
               lambda(j) = lam_max ! == exp(ln1_max)
               dlambda_dlnT(j) = 0
            else
               ln_lambda(j) = ln1
               lambda(j) = exp(ln_lambda(j))
               dlambda_dlnT(j) = &
                  dot_product(dT9fac_dlnT(:), rates% coefficients(:,i))*lambda(j)
            end if
         end do
         
      end subroutine compute_some_lambdas
      

      subroutine compute_some_inverse_lambdas( &
            num_lambdas, lo, hi, T9, rates, &
            ln_lambda, lambda, dlambda_dlnT, &
            inv_lambda, dinv_lambda_dlnT)
         use utils_lib, only: is_bad
         use chem_def, only: Tpart, npart
         use chem_lib, only: get_partition_fcn_indx
         integer, intent(in) :: num_lambdas, lo, hi ! range of rates to do
         real(dp), intent(in) :: T9
         type(reaction_data), intent(in) :: rates
         real(dp), dimension(:), intent(in) :: ln_lambda, lambda, dlambda_dlnT
         real(dp), dimension(:), intent(out) :: inv_lambda, dinv_lambda_dlnT
         
         integer :: indx,indxp
         integer :: rstart, rend, i, j
         real(dp), dimension(num_lambdas) :: A, Qratio, dQratio_dlnT
         real(dp) :: tfac, dtfac_dlnT, lnT9, T9i, dT9i_dlnT, ln1, fac1, dfac1_dlnT, dln1_dlnT,blurp
         
         include 'formats'
   
         ! find index of partition function and logarithmically interpolate
         indx = get_partition_fcn_indx(T9)
         if (indx >= npart) indx = npart-1
         if (indx < 1) then
            Qratio(1:num_lambdas) = 1d0
            dQratio_dlnT(1:num_lambdas) = 0d0
         else
            indxp = indx+1
            tfac = (T9-Tpart(indx))/(Tpart(indxp)-Tpart(indx))
            dtfac_dlnT = T9/(Tpart(indxp)-Tpart(indx))
            do i = lo, hi
               j = i+1-lo
               A(j) = rates% inverse_part(indxp,i)/rates% inverse_part(indx,i)
               blurp = pow(A(j),tfac)
               Qratio(j) = rates% inverse_part(indx,i) * blurp
               dQratio_dlnT(j) = Qratio(j)*log(A(j))*dtfac_dlnT
            end do
         end if
      
         lnT9 = log(T9)
         T9i = 1d0/T9
         dT9i_dlnT = -T9i
         
         do i = lo, hi
         
            j = i+1-lo
         
            ln1 = ln_lambda(j) + &
                  rates% inverse_coefficients(1,i) + &
                  rates% inverse_coefficients(2,i)*T9i + &
                  1.5d0*rates% inverse_exp(i)*lnT9
                  
            if (ln1 < ln1_max) then
               fac1 = exp(ln1)
               
               dln1_dlnT = dlambda_dlnT(j)/max(1d-99,lambda(j)) + &
                     rates% inverse_coefficients(2,i)*dT9i_dlnT + &
                     1.5d0*rates% inverse_exp(i)
                           
               dfac1_dlnT = dln1_dlnT*fac1
               
            else
               ln1 = ln1_max
               fac1 = lam_max ! == exp(ln1_max)
               dln1_dlnT = 0
               dfac1_dlnT = 0
            end if
            
            inv_lambda(j) = fac1*Qratio(j)            
            if (lambda(j) < 1d-99) then
               dinv_lambda_dlnT(j) = 0
               cycle
            end if
            
            dinv_lambda_dlnT(j) = dfac1_dlnT*Qratio(j) + fac1*dQratio_dlnT(j)

         end do
         
      end subroutine compute_some_inverse_lambdas
      
      
      integer function do_reaclib_lookup(handle, rates_dict) result(indx)
         ! returns first reaction index that matches handle. 
         ! there may be several following that one having the same handle.
         ! returns 0 if handle doesn't match any of the reactions
         use utils_lib, only: integer_dict_lookup
         character(len=*), intent(in) :: handle ! as in rates% reaction_handle
         type (integer_dict), pointer :: rates_dict ! from create_reaclib_rates_dict
         integer :: ierr
         ierr = 0
         call integer_dict_lookup(rates_dict, handle, indx, ierr)
         if (ierr /= 0) indx = 0
      end function do_reaclib_lookup
      
      
      subroutine do_reaclib_indices_for_reaction(handle, rates, lo, hi, ierr)
         character(len=*), intent(in) :: handle ! as in rates% reaction_handle
         type(reaction_data), intent(in) :: rates
         integer, intent(out) :: lo, hi
         integer, intent(out) :: ierr
         integer :: ir
         include 'formats'
         ierr = 0
         lo = 0
         hi = 0
         lo = do_reaclib_lookup(handle, rates% reaction_dict)
         if (lo == 0) then
            ierr = -1
            return
         end if
         hi = rates% nreactions
         do ir = lo+1, rates% nreactions
            if (rates% reaction_handle(ir) /= handle) then
               hi = ir-1
               exit
            end if
         end do
         do ir = lo-1, 1, -1
            if (rates% reaction_handle(ir) /= handle) then
               lo = ir+1
               exit
            end if
         end do
      end subroutine do_reaclib_indices_for_reaction
      
      
      subroutine do_reaclib_reaction_rates( &
            lo, hi, T9, rates, nuclides, forward_only, &
            lambda, dlambda_dlnT, &
            rlambda, drlambda_dlnT, &
            ierr)
         integer, intent(in) :: lo, hi ! from reaclib_indices_for_reaction
         real(dp), intent(in) :: T9
         type(reaction_data), intent(in) :: rates
         type(nuclide_data), intent(in) :: nuclides
         logical, intent(in) :: forward_only
         real(dp), intent(out) :: lambda, dlambda_dlnT
         real(dp), intent(out) :: rlambda, drlambda_dlnT
         integer, intent(out) :: ierr
         call reaction_rates( &
            hi-lo+1, lo, hi, T9, rates, nuclides, forward_only, &
            lambda, dlambda_dlnT, &
            rlambda, drlambda_dlnT, &
            ierr)
      end subroutine do_reaclib_reaction_rates
      


      end module reaclib_eval
