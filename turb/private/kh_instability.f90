module kh_instability

    use const_def, only: dp
    use math_lib
    use utils_lib
    ! use f95_lapack ! Might need to link properly from SDK
  
    
    implicit none
  
    private
    public :: khparams_from_fingering
    public :: deln
    public :: lmat
    public :: gamfromL
    public :: gamma_over_k
    public :: gammax_kscan
    public :: gammax_kscan_withTC
    public :: gammax_minus_lambda
    public :: gammax_minus_lambda_withTC
  
    real(dp), parameter :: CH = 1.66_dp
    real(dp), parameter :: C2 = 0.33d0
  
  contains
  
    subroutine khparams_from_fingering(w, lhat, hb, pr, db, hb_star, re, rm) ! db = \eta / k_T with \eta the magnetic diffusivity hb_star is the alfven mach number
      real(dp), intent(in) :: w, lhat, hb, pr, db 
      real(dp), intent(out) :: hb_star, re, rm
  
      hb_star = hb / pow2(w)
      re = w / (pr * lhat)
      rm = w / (db * lhat)
  
    end subroutine khparams_from_fingering
  
    function deln(k, n, delta)
      real(dp), intent(in) :: k  
      integer, intent(in) :: n  
      real(dp), intent(in) :: delta
      real(dp) :: deln
  
      deln = pow2(k) + pow2(n + delta)
  
    end function deln
  
    subroutine lmat(delta, m2, re, rm, k, n, ideal, l)
      real(dp), intent(in) :: delta, m2, re, rm, k  
      integer, intent(in) :: n
      logical, intent(in) :: ideal
      complex(dp), intent(out) :: l(2*n,2*n)
  
      real(dp) :: diss
      integer :: i, j, m, ns(n), ms(2*n)
      real(dp) :: delns(n), delns_m(2*n)
      real(dp) :: deltan, deltanm1, deltanp1
      complex(dp), parameter :: j_imag = (0d0,1d0) ! shorthand for imaginary unit

      if(ideal) then
         diss = 0d0
      else
         diss = 1d0
      end if
  
      do i = -n/2, n/2
         ns(i+n/2+1) = i
      end do
  
      do i = -n+1, n
         ms(i+n) = i
      end do
  
      do i = 1, n
         delns(i) = deln(k, ns(i), delta) 
      end do
  
      do i = 1, 2*n
         if (mod(ms(i),2) == 0) then
            delns_m(i) = deln(k, ms(i)/2, delta)
         else
            delns_m(i) = deln(k, (ms(i)-1)/2, delta)
         end if
      end do
  
      do i = 1, 2*n
         do j = 1, 2*n
            l(i,j) = (0.0_dp, 0.0_dp)
         end do
      end do
  
      do i = 1, 2*n-1
         if (i > 1 .and. i < 2*n-1) then
            deltan = delns_m(i)
            if (mod(ms(i),2) == 0) then
               l(i,i) = j_imag * (diss/re) * deltan
            else
               l(i,i) = j_imag * diss/rm * deltan
            end if
         end if
      end do

      ! pm 1 entries
      do i = 2, 2*n-2
         if (mod(ms(i),2) == 0) then
            l(i,i+1) = m2 * k
         else
            l(i,i-1) = k
         end if
      end do
      
      ! pm 2 entries
      do i = 3, 2*n-2
         deltan = delns_m(i)
         deltanp1 = delns_m(i+2)
         deltanm1 = delns_m(i-2)
         if (mod(ms(i),2) == 0) then
            l(i,i+2) = -k * (1 - deltanp1) / (2.0_dp*j_imag*deltan)
            l(i,i-2) = k * (1 - deltanm1) / (2.0_dp*j_imag*deltan)
         else
            l(i,i+2) = k/(2.0_dp*j_imag)
            l(i,i-2) = -k/(2.0_dp*j_imag)
         end if
      end do

  
      ! Boundary conditions
      i = 1
      l(i,i) = j_imag * delns_m(i) * diss / re
      l(i,i+1) = m2 * k
      l(i,i+2) = -k * (1 - delns_m(i+2)) / (2.0_dp*j_imag*delns_m(i))
  
      i = 2
      l(i,i) = j_imag * delns_m(i) * diss / rm
      l(i,i-1) = k
      l(i,i+2) = k/(2.0_dp*j_imag)
  
      i = 2*n-1
      l(i,i) = j_imag * delns_m(i) * diss / re
      l(i,i+1) = m2*k
      l(i,i-2) = k * (1 - delns_m(i-2)) / (2.0_dp*j_imag*delns_m(i))
  
      i = 2*n
      l(i,i) = j_imag * delns_m(i) * diss / rm
      l(i,i-1) = k
      l(i,i-2) = -k/(2.0_dp*j_imag)
  
    end subroutine lmat

    subroutine sams_lmat(n, f, k, m, A_Psi, A_T, A_C, Pr, tau, R0, Pm, HB, l)
      integer, intent(in)      :: n
      real(dp), intent(in)     :: f, k, m, A_psi, A_T, A_C, Pr, tau, R0, Pm, HB
      complex(dp), intent(out) :: l((2*n+1)*4,(2*n+1)*4)
  
      integer  :: i, j, dim, index, stride
      real(dp) :: D_b, k_mode, k_x, k_z, P_k_mode, N_k_mode
      complex(dp) :: PsiPsi_P, PsiPsi_N, PsiPsi, PsiT, PsiC, PsiA, &
                     TT_P, TT_N, TPsi_P, TPsi_N, TPsi, TT, &
                     CC_P, CC_N, CPsi_P, CPsi_N, CPsi, CC, &
                     AA_P, AA_N, AA, APsi
      complex(dp), parameter :: j_imag = (0d0,1d0) ! shorthand for imaginary unit

      dim = (2*n+1)*4
      D_b = Pr/Pm

      l(:,:) = (0.0_dp, 0.0_dp)

      do i = -n, n
         k_mode = sqrt(pow2(i*k + f) + pow2(m)) ! $k_m$ in latex
         k_x = i*k + f ! $(f + (m+1)k_x)$ in latex
         k_z = m ! $k_z$ in latex
        
         P_k_mode = sqrt(pow2((i + 1)*k + f) + pow2(m))  ! $k_{m+1}$ in latex
         N_k_mode = sqrt(pow2((i - 1)*k + f) + pow2(m))  ! $k_{m-1}$ in latex

         ! \psi field         
         ! -\lambda k_m^2 \psi_m -
         ! 	  i k_x  k_z   E_{\psi} \left( k_{m+1}^2 \psi_{m+1} +  k_{m-1}^2 \psi_{m-1} \right)  +
         !  i k_x^3 k_z E_{\psi}\left( \psi_{m+1} + \psi_{m-1}\right)
         !  \end{equation}
         !  $$ =  \text{Pr} k_m^4 + i\text{Pr}(f + (m+1)k_x)(T_m -C_m) -
         ! 	i H_b  k_z k_m^2 A_m  $$
         
         ! (terms are in order that they appear in latex doc)
         
         PsiPsi_P = j_imag * (-pow2(P_k_mode) * A_Psi * k_z * k + A_Psi * pow3(k) * k_z) / pow2(k_mode)
         PsiPsi_N = j_imag * (-pow2(N_k_mode) * A_Psi * k_z * k + A_Psi * pow3(k) * k_z) / pow2(k_mode)

         PsiPsi = - Pr * pow2(k_mode)
         PsiT = -j_imag * Pr * k_x / pow2(k_mode)
         PsiC = j_imag * Pr * k_x / pow2(k_mode)

         PsiA = j_imag * k_z * HB

         ! T field
         ! \lambda T_m  +  i k_x  k_z E_{\psi}(T_{m+1} + T_{m-1})  + k_x k_z E_{T} (-\psi_{m+1} + \psi_{m-1}) + i(f + m k_x) \psi_m = -k_m^2 T_m

         TT_P = j_imag * (-k * A_Psi * k_z)
         TT_N = j_imag * (-k * A_Psi * k_z)

         TPsi_P = (k * A_T * k_z)
         TPsi_N = (-k * A_T * k_z)

         TPsi = -(j_imag * k_x)
         TT = -pow2(k_mode)

         ! C field
         ! \lambda C_m  +  i k_x  k_z E_{\psi} (C_{m+1} + C_{m-1})  + k_x k_z E_{C} (-\psi_{m+1} + \psi_{m-1})  + i(f + m k_x) \frac{1}{R_0}\psi_m = -\tau k_m^2 C_m
         
         CC_P = j_imag * (-k * A_Psi * k_z)
         CC_N = j_imag * (-k * A_Psi * k_z)
         
         CPsi_P = (k * A_C * k_z)
         CPsi_N = (-k * A_C * k_z)

         CPsi = -(j_imag * k_x / R0)
         CC = -tau * pow2(k_mode)

         ! A field
         ! \lambda A +    i k_x k_z E_{\psi}\left(  A_{m+1} +  A_{m-1} \right) = -D_B A_m + k_z i \psi_m

         AA_P = -j_imag * k * k_z * A_Psi
         AA_N = -j_imag * k * k_z * A_Psi

         AA = -D_b * pow2(k_mode)

         APsi = j_imag * k_z

         !  interactions

         stride = 4
         index = (n + i) * stride + 1

         !  same mode interactions

         l(index, index) = PsiPsi
         l(index, index + 1) = PsiT
         l(index, index + 2) = PsiC

         l(index + 1, index) = TPsi
         l(index + 1, index + 1) = TT

         l(index + 2, index) = CPsi
         l(index + 2, index + 2) = CC

         l(index + 3, index) = APsi
         l(index + 3, index + 3) = AA
         l(index, index + 3) = PsiA

         ! neighboring interactions

         if (index > stride) then
            l(index, index - stride) = PsiPsi_N
            l(index + 1, index - stride) = TPsi_N
            l(index + 2, index - stride) = CPsi_N
            l(index + 1, index - stride + 1) = TT_N
            l(index + 2, index - stride + 2) = CC_N

            l(index + 3, index - stride + 3) = AA_N
         end if

         if (index + stride + 3 <= dim) then
            l(index, index + stride) = PsiPsi_P
            l(index + 1, index + stride) = TPsi_P
            l(index + 2, index + stride) = CPsi_P
            l(index + 1, index + stride + 1) = TT_P
            l(index + 2, index + stride + 2) = CC_P

            l(index + 3, index + stride + 3) = AA_P
         end if
         
      end do
      
    end subroutine sams_lmat

    subroutine richs_lmat(n, k_z, R0, Pr, tau, l_f, E_Psi, E_T, E_C, HB, DB, l)
      integer, intent(in)      :: n
      real(dp), intent(in)     :: k_z, R0, Pr, tau, l_f, E_psi, E_T, E_C, HB, DB
      complex(dp), intent(out) :: l((2*n+1)*4,(2*n+1)*4)
  
      integer  :: m, j, dim, i_p, i_p_p, i_p_m, i_t, i_t_p, &
                  i_t_m, i_c, i_c_p, i_c_m, i_a, i_a_p, i_a_m
      real(dp) :: k2_m, k2_m_p, k2_m_m
      complex(dp), parameter :: j_imag = (0d0,1d0) ! shorthand for imaginary unit

      dim = (2*n+1)*4

      l(:,:) = (0.0_dp, 0.0_dp)

      do m = -n, n
         ! Set up block indices
         i_p = 4*(m + n) + 1
         i_p_p = i_p + 4
         i_p_m = i_p - 4

         i_t = i_p + 1
         i_t_p = i_t + 4
         i_t_m = i_t - 4

         i_c = i_t + 1
         i_c_p = i_c + 4
         i_c_m = i_c - 4

         i_a = i_c + 1
         i_a_p = i_a + 4
         i_a_m = i_a - 4

         ! Evaluate wavenumbers
         k2_m = pow2(m*l_f) + pow2(k_z)
         k2_m_p = pow2((m + 1)*l_f) + pow2(k_z)
         k2_m_m = pow2((m - 1)*l_f) + pow2(k_z)
         
         ! Set matrix elements
         
         ! Start with psi_m row
         ! psi_m column
         l(i_p, i_p) = -Pr * k2_m
         ! T_m column
         l(i_p, i_t) = -j_imag * Pr * m * l_f / k2_m
         ! C_m column
         l(i_p, i_c) = j_imag * Pr * m * l_f / k2_m
         ! A_m column
         l(i_p, i_a) = j_imag * HB * k_z

         if (m > -n) then
            ! psi_{m-1} column
            l(i_p, i_p_m) = j_imag * pow3(l_f) * k_z * E_psi / k2_m - j_imag * l_f * k_z * E_psi * k2_m_m / k2_m
         end if

         if (m < n) then
            ! psi_{m+1} column
            l(i_p, i_p_p) = j_imag * pow3(l_f) * k_z * E_psi / k2_m - j_imag * l_f * k_z * E_psi * k2_m_p / k2_m
         end if

         ! T_m row
         l(i_t, i_p) = -j_imag * m * l_f
         l(i_t, i_t) = -k2_m

         if (m > -n) then
            l(i_t, i_p_m) = -l_f * k_z * E_T
            l(i_t, i_t_m) = -j_imag * l_f * k_z * E_psi
         end if
         if (m < n) then
            l(i_t, i_p_p) = l_f * k_z * E_T
            l(i_t, i_t_p) = -j_imag * l_f * k_z * E_psi
         end if

         ! C_m row
         l(i_c, i_p) = -j_imag * m * l_f / R0
         l(i_c, i_c) = -tau * k2_m

         if (m > -n) then
            l(i_c, i_p_m) = -l_f * k_z * E_C
            l(i_c, i_c_m) = -j_imag * l_f * k_z * E_psi
         end if
         if (m < n) then
            l(i_c, i_p_p) = l_f * k_z * E_C
            l(i_c, i_c_p) = -j_imag * l_f * k_z * E_psi
         end if

         ! A_m row
         l(i_a, i_p) = j_imag * k_z
         l(i_a, i_a) = -DB * k2_m

         if (m > -n) then
            l(i_a, i_a_m) = -j_imag * l_f * k_z * E_psi
         end if
         if (m < n) then
            l(i_a, i_a_p) = -j_imag * l_f * k_z * E_psi
         end if
      
      end do
      
    end subroutine richs_lmat


    subroutine lmat_block_a(n, k_z, R0, Pr, tau, l_f, E_Psi, E_T, E_C, HB, DB, l)
      integer, intent(in)      :: n
      real(dp), intent(in)     :: k_z, R0, Pr, tau, l_f, E_psi, E_T, E_C, HB, DB
      complex(dp), intent(out) :: l(4*n + 2,4*n + 2)
  
      integer  :: m, j, dim, i_p, i_p_p, i_p_m, i_t, i_t_p, &
                  i_t_m, i_c, i_c_p, i_c_m, i_a, i_a_p, i_a_m
      real(dp) :: k2_m, k2_m_p, k2_m_m, k2_1, k2_2
      complex(dp), parameter :: j_imag = (0d0,1d0) ! shorthand for imaginary unit

      dim = 4*n + 2

      l(:,:) = (0.0_dp, 0.0_dp)

      ! The basis ordering is: psi_0, A_0, psi_1+, T_1-, C_1-, A_1+, psi_2+, T_2-, ...

      ! The m=0 and 1 cases are a little funny. While they can be implemented in 
      ! the loop over m, I find it easier to set them up here first, separately.
      
      ! first, set up the 0th row, corresponding to psi_0
      k2_1 = pow2(l_f) + pow2(k_z)
      l(1, 1) = -Pr * pow2(k_z)
      l(1, 2) = j_imag * HB * k_z
      l(1, 3) = j_imag * (l_f/k_z) * E_psi * (pow2(l_f) - k2_1)
      ! next, the 1st row, corresponding to A_0
      l(2, 1) = j_imag * k_z
      l(2, 2) = -DB * pow2(k_z)
      l(2, 6) = -j_imag * l_f * k_z * E_psi

      ! next, the psi_1+ row
      k2_2 = 4 * pow2(l_f) + pow2(k_z)
      ! psi_1+ equation:
      l(3, 1) = 2*j_imag * (l_f * k_z / k2_1) * E_psi * (pow2(l_f) - pow2(k_z))
      l(3, 3) = -Pr * k2_1
      if (n > 1) then
         ! coupling to psi_2+
         l(3, 7) = j_imag * (l_f * k_z / k2_1) * E_psi * (pow2(l_f) - k2_2)
      end if
      ! coupling to T_1-
      l(3, 4) = -j_imag * Pr * (l_f / k2_1)
      ! coupling to C_1-
      l(3, 5) = j_imag * Pr * (l_f / k2_1)
      ! coupling to A_1+
      l(3, 6) = j_imag * HB * k_z
      
      ! T_1- equation:
      ! coupling to psi_0
      l(4, 1) = -2 * l_f * k_z * E_T
      ! coupling to psi_1+
      l(4, 3) = -j_imag * l_f
      ! coupling to T_1-:
      l(4, 4) = -k2_1
      if (n > 1) then
         ! coupling to psi_2+:
         l(4, 7) = l_f * k_z * E_T
         ! coupling to T_2-:
         l(4, 8) = -j_imag * l_f * k_z * E_psi
      end if

      ! C_1- equation:
      ! coupling to psi_0
      l(5, 1) = -2 * l_f * k_z * E_C
      ! coupling to psi_1+
      l(5, 3) = -j_imag * l_f / R0
      ! coupling to C_1-:
      l(5, 5) = -tau * k2_1
      if (n > 1) then
         ! coupling to psi_2+:
         l(5, 7) = l_f * k_z * E_C
         ! coupling to C_2-:
         l(5, 9) = -j_imag * l_f * k_z * E_psi
      end if

      ! A_1+ equation:
      ! A_0
      l(6, 2) = -2 * j_imag * l_f * k_z * E_psi
      ! psi_1+
      l(6, 3) = j_imag * k_z
      ! A_1+:
      l(6, 6) = -DB * k2_1
      if (n > 1) then
         ! A_2+:
         l(6, 10) = -j_imag * l_f * k_z * E_psi
      end if
      
      if (n > 1) then
         do m = 2, n
            ! Set up block indices
            i_p = 4*m - 1
            i_p_p = i_p + 4
            i_p_m = i_p - 4

            i_t = i_p + 1
            i_t_p = i_t + 4
            i_t_m = i_t - 4

            i_c = i_t + 1
            i_c_p = i_c + 4
            i_c_m = i_c - 4

            i_a = i_c + 1
            i_a_p = i_a + 4
            i_a_m = i_a - 4

            ! Evaluate wavenumbers
            k2_m = pow2(m*l_f) + pow2(k_z)
            k2_m_p = pow2((m + 1)*l_f) + pow2(k_z)
            k2_m_m = pow2((m - 1)*l_f) + pow2(k_z)
            
            ! Set matrix elements
            
            ! Start with psi_m row
            ! psi_m column
            l(i_p, i_p) = -Pr * k2_m
            ! T_m column
            l(i_p, i_t) = -j_imag * Pr * m * l_f / k2_m
            ! C_m column
            l(i_p, i_c) = j_imag * Pr * m * l_f / k2_m
            ! A_m column
            l(i_p, i_a) = j_imag * HB * k_z
            ! psi_{m-1} column
            ! l(i_p, i_p_m) = j_imag * pow3(l_f) * k_z * E_psi / k2_m - j_imag * l_f * k_z * E_psi * k2_m_m / k2_m
            l(i_p, i_p_m) = j_imag * l_f * k_z * (E_psi / k2_m) * (pow2(l_f) - k2_m_m)
            if (m < n) then
               ! psi_{m+1} column
               ! l(i_p, i_p_p) = j_imag * pow3(l_f) * k_z * E_psi / k2_m - j_imag * l_f * k_z * E_psi * k2_m_p / k2_m
               l(i_p, i_p_p) = j_imag * l_f * k_z * (E_psi / k2_m) * (pow2(l_f) - k2_m_p)
            end if

            ! T_m row
            l(i_t, i_p) = -j_imag * m * l_f
            l(i_t, i_t) = -k2_m
            l(i_t, i_p_m) = -l_f * k_z * E_T
            l(i_t, i_t_m) = -j_imag * l_f * k_z * E_psi
            if (m < n) then
               l(i_t, i_p_p) = l_f * k_z * E_T
               l(i_t, i_t_p) = -j_imag * l_f * k_z * E_psi
            end if

            ! C_m row
            l(i_c, i_p) = -j_imag * m * l_f / R0
            l(i_c, i_c) = -tau * k2_m
            l(i_c, i_p_m) = -l_f * k_z * E_C
            l(i_c, i_c_m) = -j_imag * l_f * k_z * E_psi
            if (m < n) then
               l(i_c, i_p_p) = l_f * k_z * E_C
               l(i_c, i_c_p) = -j_imag * l_f * k_z * E_psi
            end if

            ! A_m row
            l(i_a, i_p) = j_imag * k_z
            l(i_a, i_a) = -DB * k2_m
            l(i_a, i_a_m) = -j_imag * l_f * k_z * E_psi
            if (m < n) then
               l(i_a, i_a_p) = -j_imag * l_f * k_z * E_psi
            end if
         
         end do
      end if
      
    end subroutine lmat_block_a


    subroutine lmat_block_b(n, k_z, R0, Pr, tau, l_f, E_Psi, E_T, E_C, HB, DB, l)
      integer, intent(in)      :: n
      real(dp), intent(in)     :: k_z, R0, Pr, tau, l_f, E_psi, E_T, E_C, HB, DB
      complex(dp), intent(out) :: l(4*n + 2,4*n + 2)
  
      integer  :: m, j, dim, i_p, i_p_p, i_p_m, i_t, i_t_p, &
                  i_t_m, i_c, i_c_p, i_c_m, i_a, i_a_p, i_a_m
      real(dp) :: k2_m, k2_m_p, k2_m_m, k2_1, k2_2
      complex(dp), parameter :: j_imag = (0d0,1d0) ! shorthand for imaginary unit

      dim = 4*n + 2

      l(:,:) = (0.0_dp, 0.0_dp)

      ! The basis ordering is: T_0, C_0, psi_1-, T_1+, C_1+, A_1-, psi_2-, T_2+, ...

      ! The m=0 and 1 cases are a little funny. While they can be implemented in 
      ! the loop over m, I find it easier to set them up here first, separately.
      
      ! first, set up the 0th row, corresponding to T_0
      l(1, 1) = -pow2(k_z)
      l(1, 3) = l_f * k_z * E_T
      l(1, 4) = -j_imag * l_f * k_z * E_psi
      ! next, the 1st row, corresponding to C_0
      l(2, 2) = -tau * pow2(k_z)
      l(2, 3) = l_f * k_z * E_C
      l(2, 5) = -j_imag * l_f * k_z * E_psi

      ! next, the psi_1- row
      k2_1 = pow2(l_f) + pow2(k_z)
      k2_2 = 4 * pow2(l_f) + pow2(k_z)
      ! psi_1- equation:
      l(3, 3) = -Pr * k2_1
      if (n > 1) then
         ! coupling to psi_2-
         l(3, 7) = j_imag * (l_f * k_z / k2_1) * E_psi * (pow2(l_f) - k2_2)
      end if
      ! coupling to T_1+
      l(3, 4) = -j_imag * Pr * (l_f / k2_1)
      ! coupling to C_1+
      l(3, 5) = j_imag * Pr * (l_f / k2_1)
      ! coupling to A_1-
      l(3, 6) = j_imag * HB * k_z
      
      ! T_1+ equation:
      ! coupling to T_0
      l(4, 1) = -2 * j_imag * l_f * k_z * E_psi
      ! coupling to psi_1-
      l(4, 3) = -j_imag * l_f
      ! coupling to T_1+:
      l(4, 4) = -k2_1
      if (n > 1) then
         ! coupling to psi_2-:
         l(4, 7) = l_f * k_z * E_T
         ! coupling to T_2+:
         l(4, 8) = -j_imag * l_f * k_z * E_psi
      end if

      ! C_1+ equation:
      ! coupling to C_0
      l(5, 2) = -2 * j_imag * l_f * k_z * E_psi
      ! coupling to psi_1-
      l(5, 3) = -j_imag * l_f / R0
      ! coupling to C_1+:
      l(5, 5) = -tau * k2_1
      if (n > 1) then
         ! coupling to psi_2-:
         l(5, 7) = l_f * k_z * E_C
         ! coupling to C_2+:
         l(5, 9) = -j_imag * l_f * k_z * E_psi
      end if

      ! A_1- equation:
      ! psi_1-
      l(6, 3) = j_imag * k_z
      ! A_1-:
      l(6, 6) = -DB * k2_1
      if (n > 1) then
         ! A_2-:
         l(6, 10) = -j_imag * l_f * k_z * E_psi
      end if
      
      ! the rest of the matrix is identical to the other sub-block
      if (n > 1) then
         do m = 2, n
            ! Set up block indices
            i_p = 4*m - 1
            i_p_p = i_p + 4
            i_p_m = i_p - 4

            i_t = i_p + 1
            i_t_p = i_t + 4
            i_t_m = i_t - 4

            i_c = i_t + 1
            i_c_p = i_c + 4
            i_c_m = i_c - 4

            i_a = i_c + 1
            i_a_p = i_a + 4
            i_a_m = i_a - 4

            ! Evaluate wavenumbers
            k2_m = pow2(m*l_f) + pow2(k_z)
            k2_m_p = pow2((m + 1)*l_f) + pow2(k_z)
            k2_m_m = pow2((m - 1)*l_f) + pow2(k_z)
            
            ! Set matrix elements
            
            ! Start with psi_m row
            ! psi_m column
            l(i_p, i_p) = -Pr * k2_m
            ! T_m column
            l(i_p, i_t) = -j_imag * Pr * m * l_f / k2_m
            ! C_m column
            l(i_p, i_c) = j_imag * Pr * m * l_f / k2_m
            ! A_m column
            l(i_p, i_a) = j_imag * HB * k_z
            ! psi_{m-1} column
            ! l(i_p, i_p_m) = j_imag * pow3(l_f) * k_z * E_psi / k2_m - j_imag * l_f * k_z * E_psi * k2_m_m / k2_m
            l(i_p, i_p_m) = j_imag * l_f * k_z * (E_psi / k2_m) * (pow2(l_f) - k2_m_m)
            if (m < n) then
               ! psi_{m+1} column
               ! l(i_p, i_p_p) = j_imag * pow3(l_f) * k_z * E_psi / k2_m - j_imag * l_f * k_z * E_psi * k2_m_p / k2_m
               l(i_p, i_p_p) = j_imag * l_f * k_z * (E_psi / k2_m) * (pow2(l_f) - k2_m_p)
            end if

            ! T_m row
            l(i_t, i_p) = -j_imag * m * l_f
            l(i_t, i_t) = -k2_m
            l(i_t, i_p_m) = -l_f * k_z * E_T
            l(i_t, i_t_m) = -j_imag * l_f * k_z * E_psi
            if (m < n) then
               l(i_t, i_p_p) = l_f * k_z * E_T
               l(i_t, i_t_p) = -j_imag * l_f * k_z * E_psi
            end if

            ! C_m row
            l(i_c, i_p) = -j_imag * m * l_f / R0
            l(i_c, i_c) = -tau * k2_m
            l(i_c, i_p_m) = -l_f * k_z * E_C
            l(i_c, i_c_m) = -j_imag * l_f * k_z * E_psi
            if (m < n) then
               l(i_c, i_p_p) = l_f * k_z * E_C
               l(i_c, i_c_p) = -j_imag * l_f * k_z * E_psi
            end if

            ! A_m row
            l(i_a, i_p) = j_imag * k_z
            l(i_a, i_a) = -DB * k2_m
            l(i_a, i_a_m) = -j_imag * l_f * k_z * E_psi
            if (m < n) then
               l(i_a, i_a_p) = -j_imag * l_f * k_z * E_psi
            end if
         
         end do
      end if
      
    end subroutine lmat_block_b


    subroutine lmat_block_a_LPN(n, k_z, R0, Pr, tau, l_f, E_Psi, E_C, HB, DB, l)
      ! same as lmat_block_a above, but implements the low peclet number limit
      ! to remove the T equation
      integer, intent(in)      :: n
      real(dp), intent(in)     :: k_z, R0, Pr, tau, l_f, E_psi, E_C, HB, DB
      complex(dp), intent(out) :: l(3*n + 2, 3*n + 2)
  
      integer  :: m, j, dim, i_p, i_p_p, i_p_m, &
                  i_c, i_c_p, i_c_m, i_a, i_a_p, i_a_m
      real(dp) :: k2_m, k2_m_p, k2_m_m, k2_1, k2_2
      complex(dp), parameter :: j_imag = (0d0,1d0) ! shorthand for imaginary unit

      dim = 3*n + 2

      l(:,:) = (0.0_dp, 0.0_dp)

      ! The basis ordering is: psi_0, A_0, psi_1+, C_1-, A_1+, psi_2+, C_2-, ...

      ! The m=0 and 1 cases are a little funny. While they can be implemented in 
      ! the loop over m, I find it easier to set them up here first, separately.
      
      ! first, set up the 0th row, corresponding to psi_0
      k2_1 = pow2(l_f) + pow2(k_z)
      l(1, 1) = -Pr * pow2(k_z)
      l(1, 2) = j_imag * HB * k_z
      l(1, 3) = j_imag * (l_f/k_z) * E_psi * (pow2(l_f) - k2_1)
      ! next, the 1st row, corresponding to A_0
      l(2, 1) = j_imag * k_z
      l(2, 2) = -DB * pow2(k_z)
      l(2, 5) = -j_imag * l_f * k_z * E_psi

      ! next, the psi_1+ row
      k2_2 = 4 * pow2(l_f) + pow2(k_z)
      ! psi_1+ equation:
      l(3, 1) = 2*j_imag * (l_f * k_z / k2_1) * E_psi * (pow2(l_f) - pow2(k_z))
      l(3, 3) = -Pr * (k2_1 + pow2(l_f / k2_1))
      if (n > 1) then
         ! coupling to psi_2+
         l(3, 6) = j_imag * (l_f * k_z / k2_1) * E_psi * (pow2(l_f) - k2_2)
      end if
      ! coupling to C_1-
      l(3, 4) = j_imag * Pr * (l_f / k2_1)
      ! coupling to A_1+
      l(3, 5) = j_imag * HB * k_z

      ! C_1- equation:
      ! coupling to psi_0
      l(4, 1) = -2 * l_f * k_z * E_C
      ! coupling to psi_1+
      l(4, 3) = -j_imag * l_f / R0
      ! coupling to C_1-:
      l(4, 4) = -tau * k2_1
      if (n > 1) then
         ! coupling to psi_2+:
         l(4, 6) = l_f * k_z * E_C
         ! coupling to C_2-:
         l(4, 7) = -j_imag * l_f * k_z * E_psi
      end if

      ! A_1+ equation:
      ! A_0
      l(5, 2) = -2 * j_imag * l_f * k_z * E_psi
      ! psi_1+
      l(5, 3) = j_imag * k_z
      ! A_1+:
      l(5, 5) = -DB * k2_1
      if (n > 1) then
         ! A_2+:
         l(5, 8) = -j_imag * l_f * k_z * E_psi
      end if
      
      if (n > 1) then
         do m = 2, n
            ! Set up block indices
            i_p = 3*m
            i_p_p = i_p + 3
            i_p_m = i_p - 3

            i_c = i_p + 1
            i_c_p = i_c + 3
            i_c_m = i_c - 3

            i_a = i_c + 1
            i_a_p = i_a + 3
            i_a_m = i_a - 3

            ! Evaluate wavenumbers
            k2_m = pow2(m*l_f) + pow2(k_z)
            k2_m_p = pow2((m + 1)*l_f) + pow2(k_z)
            k2_m_m = pow2((m - 1)*l_f) + pow2(k_z)
            
            ! Set matrix elements
            
            ! Start with psi_m row
            ! psi_m column
            l(i_p, i_p) = -Pr * (k2_m + pow2(m * l_f / k2_m))
            ! C_m column
            l(i_p, i_c) = j_imag * Pr * m * l_f / k2_m
            ! A_m column
            l(i_p, i_a) = j_imag * HB * k_z
            ! psi_{m-1} column
            ! l(i_p, i_p_m) = j_imag * pow3(l_f) * k_z * E_psi / k2_m - j_imag * l_f * k_z * E_psi * k2_m_m / k2_m
            l(i_p, i_p_m) = j_imag * l_f * k_z * (E_psi / k2_m) * (pow2(l_f) - k2_m_m)
            if (m < n) then
               ! psi_{m+1} column
               ! l(i_p, i_p_p) = j_imag * pow3(l_f) * k_z * E_psi / k2_m - j_imag * l_f * k_z * E_psi * k2_m_p / k2_m
               l(i_p, i_p_p) = j_imag * l_f * k_z * (E_psi / k2_m) * (pow2(l_f) - k2_m_p)
            end if

            ! C_m row
            l(i_c, i_p) = -j_imag * m * l_f / R0
            l(i_c, i_c) = -tau * k2_m
            l(i_c, i_p_m) = -l_f * k_z * E_C
            l(i_c, i_c_m) = -j_imag * l_f * k_z * E_psi
            if (m < n) then
               l(i_c, i_p_p) = l_f * k_z * E_C
               l(i_c, i_c_p) = -j_imag * l_f * k_z * E_psi
            end if

            ! A_m row
            l(i_a, i_p) = j_imag * k_z
            l(i_a, i_a) = -DB * k2_m
            l(i_a, i_a_m) = -j_imag * l_f * k_z * E_psi
            if (m < n) then
               l(i_a, i_a_p) = -j_imag * l_f * k_z * E_psi
            end if
         
         end do
      end if
      
    end subroutine lmat_block_a_LPN


    subroutine lmat_block_b_LPN(n, k_z, R0, Pr, tau, l_f, E_Psi, E_C, HB, DB, l)
      ! same as lmat_block_b above, but implements the low peclet number limit
      ! to remove the T equation
      integer, intent(in)      :: n
      real(dp), intent(in)     :: k_z, R0, Pr, tau, l_f, E_psi, E_C, HB, DB
      complex(dp), intent(out) :: l(3*n + 1, 3*n + 1)
  
      integer  :: m, j, dim, i_p, i_p_p, i_p_m, &
                  i_c, i_c_p, i_c_m, i_a, i_a_p, i_a_m
      real(dp) :: k2_m, k2_m_p, k2_m_m, k2_1, k2_2
      complex(dp), parameter :: j_imag = (0d0,1d0) ! shorthand for imaginary unit

      dim = 3*n + 1

      l(:,:) = (0.0_dp, 0.0_dp)

      ! The basis ordering is: C_0, psi_1-, C_1+, A_1-, psi_2-, C_2+, ...

      ! The m=0 and 1 cases are a little funny. While they can be implemented in 
      ! the loop over m, I find it easier to set them up here first, separately.
      
      ! first, set up the 0th row, corresponding to C_0
      l(1, 1) = -tau * pow2(k_z)
      l(1, 2) = l_f * k_z * E_C
      l(1, 3) = -j_imag * l_f * k_z * E_psi

      ! next, the psi_1- row
      k2_1 = pow2(l_f) + pow2(k_z)
      k2_2 = 4 * pow2(l_f) + pow2(k_z)
      ! psi_1- equation:
      l(2, 2) = -Pr * (k2_1 + pow2(l_f / k2_1))
      if (n > 1) then
         ! coupling to psi_2-
         l(2, 5) = j_imag * (l_f * k_z / k2_1) * E_psi * (pow2(l_f) - k2_2)
      end if
      ! coupling to C_1+
      l(2, 3) = j_imag * Pr * (l_f / k2_1)
      ! coupling to A_1-
      l(2, 4) = j_imag * HB * k_z

      ! C_1+ equation:
      ! coupling to C_0
      l(3, 1) = -2 * j_imag * l_f * k_z * E_psi
      ! coupling to psi_1-
      l(3, 2) = -j_imag * l_f / R0
      ! coupling to C_1+:
      l(3, 3) = -tau * k2_1
      if (n > 1) then
         ! coupling to psi_2-:
         l(3, 5) = l_f * k_z * E_C
         ! coupling to C_2+:
         l(3, 6) = -j_imag * l_f * k_z * E_psi
      end if

      ! A_1- equation:
      ! psi_1-
      l(4, 2) = j_imag * k_z
      ! A_1-:
      l(4, 4) = -DB * k2_1
      if (n > 1) then
         ! A_2-:
         l(4, 7) = -j_imag * l_f * k_z * E_psi
      end if
      
      ! the rest of the matrix is identical to the other sub-block
      if (n > 1) then
         do m = 2, n
            ! Set up block indices
            i_p = 3*m - 1
            i_p_p = i_p + 3
            i_p_m = i_p - 3

            i_c = i_p + 1
            i_c_p = i_c + 3
            i_c_m = i_c - 3

            i_a = i_c + 1
            i_a_p = i_a + 3
            i_a_m = i_a - 3

            ! Evaluate wavenumbers
            k2_m = pow2(m*l_f) + pow2(k_z)
            k2_m_p = pow2((m + 1)*l_f) + pow2(k_z)
            k2_m_m = pow2((m - 1)*l_f) + pow2(k_z)
            
            ! Set matrix elements
            
            ! Start with psi_m row
            ! psi_m column
            l(i_p, i_p) = -Pr * (k2_m + pow2(m * l_f / k2_m))
            ! C_m column
            l(i_p, i_c) = j_imag * Pr * m * l_f / k2_m
            ! A_m column
            l(i_p, i_a) = j_imag * HB * k_z
            ! psi_{m-1} column
            ! l(i_p, i_p_m) = j_imag * pow3(l_f) * k_z * E_psi / k2_m - j_imag * l_f * k_z * E_psi * k2_m_m / k2_m
            l(i_p, i_p_m) = j_imag * l_f * k_z * (E_psi / k2_m) * (pow2(l_f) - k2_m_m)
            if (m < n) then
               ! psi_{m+1} column
               ! l(i_p, i_p_p) = j_imag * pow3(l_f) * k_z * E_psi / k2_m - j_imag * l_f * k_z * E_psi * k2_m_p / k2_m
               l(i_p, i_p_p) = j_imag * l_f * k_z * (E_psi / k2_m) * (pow2(l_f) - k2_m_p)
            end if

            ! C_m row
            l(i_c, i_p) = -j_imag * m * l_f / R0
            l(i_c, i_c) = -tau * k2_m
            l(i_c, i_p_m) = -l_f * k_z * E_C
            l(i_c, i_c_m) = -j_imag * l_f * k_z * E_psi
            if (m < n) then
               l(i_c, i_p_p) = l_f * k_z * E_C
               l(i_c, i_c_p) = -j_imag * l_f * k_z * E_psi
            end if

            ! A_m row
            l(i_a, i_p) = j_imag * k_z
            l(i_a, i_a) = -DB * k2_m
            l(i_a, i_a_m) = -j_imag * l_f * k_z * E_psi
            if (m < n) then
               l(i_a, i_a_p) = -j_imag * l_f * k_z * E_psi
            end if
         
         end do
      end if
      
    end subroutine lmat_block_b_LPN


    subroutine lmat_block_b_LPN_QS(n, k_z, R0, Pr, tau, l_f, E_Psi, E_C, HB, DB, l)
      ! same as lmat_block_b_LPN above, but implements the low magnetic Reynolds number
      ! limit (AKA the quasistatic approximation/limit) to remove the A equation
      integer, intent(in)      :: n
      real(dp), intent(in)     :: k_z, R0, Pr, tau, l_f, E_psi, E_C, HB, DB
      complex(dp), intent(out) :: l(2*n + 1, 2*n + 1)
  
      integer  :: m, j, dim, i_p, i_p_p, i_p_m, &
                  i_c, i_c_p, i_c_m
      real(dp) :: k2_m, k2_m_p, k2_m_m, k2_1, k2_2
      complex(dp), parameter :: j_imag = (0d0,1d0) ! shorthand for imaginary unit

      dim = 2*n + 1

      l(:,:) = (0.0_dp, 0.0_dp)

      ! The basis ordering is: C_0, psi_1-, C_1+, psi_2-, C_2+, ...

      ! The m=0 and 1 cases are a little funny. While they can be implemented in 
      ! the loop over m, I find it easier to set them up here first, separately.
      
      ! first, set up the 0th row, corresponding to C_0
      l(1, 1) = -tau * pow2(k_z)
      l(1, 2) = l_f * k_z * E_C
      l(1, 3) = -j_imag * l_f * k_z * E_psi

      ! next, the psi_1- row
      k2_1 = pow2(l_f) + pow2(k_z)
      k2_2 = 4 * pow2(l_f) + pow2(k_z)
      ! psi_1- equation:
      l(2, 2) = -Pr * (k2_1 + pow2(l_f / k2_1)) - (HB / DB) * pow2(k_z) / k2_1
      if (n > 1) then
         ! coupling to psi_2-
         l(2, 4) = j_imag * (l_f * k_z / k2_1) * E_psi * (pow2(l_f) - k2_2)
      end if
      ! coupling to C_1+
      l(2, 3) = j_imag * Pr * (l_f / k2_1)

      ! C_1+ equation:
      ! coupling to C_0
      l(3, 1) = -2 * j_imag * l_f * k_z * E_psi
      ! coupling to psi_1-
      l(3, 2) = -j_imag * l_f / R0
      ! coupling to C_1+:
      l(3, 3) = -tau * k2_1
      if (n > 1) then
         ! coupling to psi_2-:
         l(3, 4) = l_f * k_z * E_C
         ! coupling to C_2+:
         l(3, 5) = -j_imag * l_f * k_z * E_psi
      end if
      
      ! the rest of the matrix
      if (n > 1) then
         do m = 2, n
            ! Set up block indices
            i_p = 2*m
            i_p_p = i_p + 2
            i_p_m = i_p - 2

            i_c = i_p + 1
            i_c_p = i_c + 2
            i_c_m = i_c - 2

            ! Evaluate wavenumbers
            k2_m = pow2(m*l_f) + pow2(k_z)
            k2_m_p = pow2((m + 1)*l_f) + pow2(k_z)
            k2_m_m = pow2((m - 1)*l_f) + pow2(k_z)
            
            ! Set matrix elements
            
            ! Start with psi_m row
            ! psi_m column
            l(i_p, i_p) = -Pr * (k2_m + pow2(m * l_f / k2_m)) - (HB / DB) * pow2(k_z) / k2_m
            ! C_m column
            l(i_p, i_c) = j_imag * Pr * m * l_f / k2_m
            ! psi_{m-1} column
            l(i_p, i_p_m) = j_imag * l_f * k_z * (E_psi / k2_m) * (pow2(l_f) - k2_m_m)
            if (m < n) then
               ! psi_{m+1} column
               l(i_p, i_p_p) = j_imag * l_f * k_z * (E_psi / k2_m) * (pow2(l_f) - k2_m_p)
            end if

            ! C_m row
            l(i_c, i_p) = -j_imag * m * l_f / R0
            l(i_c, i_c) = -tau * k2_m
            l(i_c, i_p_m) = -l_f * k_z * E_C
            l(i_c, i_c_m) = -j_imag * l_f * k_z * E_psi
            if (m < n) then
               l(i_c, i_p_p) = l_f * k_z * E_C
               l(i_c, i_c_p) = -j_imag * l_f * k_z * E_psi
            end if
         
         end do
      end if
      
    end subroutine lmat_block_b_LPN_QS
    
    
    function gamfromL(L, N, return_maxreal, withmode) result(gam)
      ! Note: withmode is not needed (is useful in the Python implementation for debugging/looking into things more closely)
      complex(dp), intent(in) :: L(:,:)
      integer, intent(in) :: N ! size of NxN matrix L
      logical, intent(in), optional :: return_maxreal, withmode
      real(dp) :: gam 
      !complex(dp), allocatable, optional :: mode(:)
  
      complex(dp) :: omega(N) ! eigenvalues
      complex(dp) :: VL(N,N), VR(N,N) ! left and right eigenvectors (not used, dummies for ZGEEV interface)
      complex(dp), allocatable :: work(:) ! workspace
      real(dp) :: rwork(2*N) ! workspace
      integer :: info, lwork

      ! Query ZGEEV for optimal workspace size
      lwork = -1 
      allocate(work(1))
      call ZGEEV('N','N',N,L,N,omega,VL,N,VR,N,work,lwork,rwork,info)
      
      ! Allocate optimal workspace
      lwork = INT(work(1)) 
      deallocate(work)
      allocate(work(lwork))
      
      ! lapack get eigenvalues (ZGEEV is for general complex double precision matrices)
      ! first two arguments specify not to compute ('N') either left or right eigenvectors
      call ZGEEV('N','N',N,L,N,omega,VL,N,VR,N,work,lwork,rwork,info)
      
      deallocate(work)
      
      if (info /= 0) then
         write(*,*) 'ZGEEV failed in gamfromL' 
      end if

!!$      if (present(withmode)) then
!!$         i = maxloc(aimag(omega),1)
!!$         gam = -aimag(omega(i))
!!$         allocate(mode(size(omega)))
!!$         mode = v(:,i)
!!$      else

      if(present(return_maxreal)) then
         gam = maxval(realpart(omega))
      else
         gam = maxval(-aimag(omega))
      end if

!!$      end if

    end function gamfromL
    
    function gamma_over_k(delta, m2, re, rm, ks, n, ideal) result(gamk)
      real(dp), intent(in) :: delta, m2, re, rm
      real(dp), intent(in) :: ks(:)
      integer, intent(in) :: n
      logical, intent(in) :: ideal
  
      complex(dp) :: l_result(2*n,2*n)
      real(dp) :: gamk(size(ks))
      integer :: i

      do i = 1, size(ks)
         call lmat(delta, m2, re, rm, ks(i), n, ideal, l_result)
         gamk(i) = gamfromL(l_result,2*n)
      end do
  
    end function gamma_over_k

    function gamma_over_k_withTC(w, hb, db, pr, tau, R0, lamhat, lhat, ks, n, safety) result(gamk)
      real(dp), intent(in) :: w, hb, db, pr, tau, R0, lamhat, lhat
      real(dp), intent(in) :: ks(:)
      integer, intent(in) :: n, safety
  
      ! complex(dp) :: l_result(4*n,4*n)
      ! complex(dp) :: l_result(2*n,2*n)  ! each block is half the size of the original matrix
      real(dp) :: gamk(size(ks))
      real(dp) :: kz, l2hat, A_psi, A_T, A_C, gam_a, gam_b
      integer :: i, n_Sam, ierr
      ! Now the two blocks have two different sizes; if this is ugly, could just be ok with one
      ! matrix having an extra row and column of zeroes
      ! complex(dp) :: l_result_a(3*n_Sam+2, 3*n_Sam+2), l_result_b(3*n_Sam+1, 3*n_Sam+1)
      ! complex(dp) :: l_result_a(3*((n-1)/2)+2, 3*((n-1)/2)+2), l_result_b(3*((n-1)/2)+1, 3*((n-1)/2)+1)
      ! complex(dp) :: l_result_a(3*((n-1)/2)+2, 3*((n-1)/2)+2), l_result_b(n, n)
      
      ! TODO: Evan, can you tell me what the reasonable thing to do here is? l_result's size depends on 
      ! the value of safety, so we could either make it a dynamic array as I've done here, or make it
      ! so each choice of "safety" calls a different equivalent of "gamma_over_k_withTC", where the
      ! appropriate size of l_result is declared
      complex(dp), allocatable :: l_result_a(:, :), l_result_b(:, :)

      n_Sam = (n-1)/2  ! Sam's definition of n is different than Adrian's

      if (safety == 0) then
         allocate(l_result_b(n, n))
         ! and no need to allocate l_result_a
      elseif (safety == 1) then
         allocate(l_result_a(3*n_Sam+2, 3*n_Sam+2))
         allocate(l_result_b(n, n))
      elseif (safety == 2) then
         allocate(l_result_a(3*n_Sam+2, 3*n_Sam+2))
         allocate(l_result_b(3*n_Sam+1, 3*n_Sam+1))
      elseif (safety == 3) then
         allocate(l_result_a(2*n, 2*n))
         ! and no need to allocate l_result_b: we can just reuse this array
      end if  ! TODO: there should be an "ELSE" statement here raising an error if safety isn't one of these values


      l2hat = pow2(lhat)
      
      A_psi = w / (2*lhat)
      A_T = -lhat * A_psi / (lamhat + l2hat)
      A_C = -lhat * A_psi / (R0 * (lamhat + tau * l2hat))

      do i = 1, size(ks)
         kz = ks(i)*lhat
         if (safety == 0) then
            call lmat_block_b_LPN_QS(n_Sam, kz, R0, pr, tau, lhat, A_Psi, A_C, hb, db, l_result_b)
            gamk(i) = gamfromL(l_result_b, n, .true.)
         elseif (safety == 1) then
            call lmat_block_a_LPN(n_Sam, kz, R0, pr, tau, lhat, A_Psi, A_C, hb, db, l_result_a)
            gam_a = gamfromL(l_result_a, 3*n_Sam+2, .true.)
            call lmat_block_b_LPN_QS(n_Sam, kz, R0, pr, tau, lhat, A_Psi, A_C, hb, db, l_result_b)
            gam_b = gamfromL(l_result_b, n, .true.)
            gamk(i) = MAX(gam_a, gam_b)
         elseif (safety == 2) then
            call lmat_block_a_LPN(n_Sam, kz, R0, pr, tau, lhat, A_Psi, A_C, hb, db, l_result_a)
            gam_a = gamfromL(l_result_a, 3*n_Sam+2, .true.)
            call lmat_block_b_LPN(n_Sam, kz, R0, pr, tau, lhat, A_Psi, A_C, hb, db, l_result_b)
            gam_b = gamfromL(l_result_b, 3*n_Sam+1, .true.)
            gamk(i) = MAX(gam_a, gam_b)
         elseif (safety == 3) then
            call lmat_block_a(n_Sam, kz, R0, pr, tau, lhat, A_Psi, A_T, A_C, hb, db, l_result_a)
            gam_a = gamfromL(l_result_a, 2*n, .true.)
            call lmat_block_b(n_Sam, kz, R0, pr, tau, lhat, A_Psi, A_T, A_C, hb, db, l_result_a)
            gam_b = gamfromL(l_result_a, 2*n, .true.)
            gamk(i) = MAX(gam_a, gam_b)
         end if
         ! ! sams_lmat has dimension (2*n_sam+1)*4 = 4*n
         ! ! call sams_lmat(n_Sam, 0d0, lhat, kz, A_Psi, A_T, A_C, pr, tau, R0, pr/db, hb, l_result)
         ! ! call richs_lmat(n_Sam, kz, R0, pr, tau, lhat, A_Psi, A_T, A_C, hb, db, l_result)
         ! ! gamk(i) = gamfromL(l_result,4*n,.true.)
         ! ! call lmat_block_a(n_Sam, kz, R0, pr, tau, lhat, A_Psi, A_T, A_C, hb, db, l_result)
         ! ! gam_a = gamfromL(l_result, 2*n, .true.)
         ! ! call lmat_block_b(n_Sam, kz, R0, pr, tau, lhat, A_Psi, A_T, A_C, hb, db, l_result)
         ! ! gam_b = gamfromL(l_result, 2*n, .true.)
         ! call lmat_block_a_LPN(n_Sam, kz, R0, pr, tau, lhat, A_Psi, A_C, hb, db, l_result_a)
         ! gam_a = gamfromL(l_result_a, 3*n_Sam+2, .true.)
         ! ! call lmat_block_b_LPN(n_Sam, kz, R0, pr, tau, lhat, A_Psi, A_C, hb, db, l_result_b)
         ! ! gam_b = gamfromL(l_result_b, 3*n_Sam+1, .true.)
         ! call lmat_block_b_LPN_QS(n_Sam, kz, R0, pr, tau, lhat, A_Psi, A_C, hb, db, l_result_b)
         ! gam_b = gamfromL(l_result_b, 2*n_Sam+1, .true.)
         ! gamk(i) = MAX(gam_a, gam_b)
         
      end do

      if (safety > 0) then
         deallocate(l_result_a)
      end if
      if (safety < 3) then
         deallocate(l_result_b)
      end if
  
    end function gamma_over_k_withTC

  
    function gammax_kscan(delta, m2, re, rm, ks, n, ideal, badks_except) result(gammax)
      real(dp), intent(in) :: delta, m2, re, rm
      real(dp), intent(in) :: ks(:)  
      integer, intent(in) :: n
      logical, intent(in) :: ideal, badks_except
  
      real(dp) :: gammax
      real(dp) :: gamk(size(ks))
      integer :: i
  
      gamk = gamma_over_k(delta, m2, re, rm, ks, n, ideal)
      i = maxloc(gamk,1)
      gammax = gamk(i)
  
      if (badks_except .and. gammax > 0.0_dp) then
         if (i == 1 .or. i == size(ks)) then
            write(*,*) 'WARNING: most unstable growth at edge of k range'
            stop ! may need to get rid of this stop statement
         end if
      end if
  
    end function gammax_kscan

   function gammax_kscan_withTC(w, hb, db, pr, tau, R0, lamhat, lhat, ks, n, badks_exception, safety) result(gammax)
      real(dp), intent(in) :: w, hb, db, pr, tau, R0, lamhat, lhat
      real(dp), intent(in) :: ks(:)
      integer, intent(in) :: n, safety
      logical, intent(in) :: badks_exception
  
      real(dp) :: gammax
      real(dp) :: gamk(size(ks))
      integer :: i
  
      gamk = gamma_over_k_withTC(w, hb, db, pr, tau, R0, lamhat, lhat, ks, n, safety)
      i = maxloc(gamk,1)
      gammax = gamk(i)
  
      if (badks_exception .and. gammax > 0.0_dp) then
         if (i == 1 .or. i == size(ks)) then
            write(*,*) 'WARNING: most unstable growth at edge of k range'
            write(*,*) 'have not fully implemented logic for badks_exception'
            call mesa_error(__FILE__,__LINE__)
         end if
      end if
  
    end function gammax_kscan_withTC
    
      
!!$        function sigma_from_fingering_params(delta, w, hb, db, pr, tau, r0, k_star, n, withmode) result(sigma)
!!$          real(dp), intent(in) :: delta, w, hb, db, pr, tau, r0, k_star 
!!$          integer, intent(in) :: n
!!$          logical, intent(in), optional :: withmode
!!$          real(dp) :: sigma
!!$          complex(dp), allocatable, optional :: mode(:)
!!$      
!!$          real(dp) :: lamhat, lhat, l2hat, m2, re, rm
!!$          complex(dp) :: l(2*n,2*n)
!!$      
!!$          call gaml2max(pr, tau, r0, lamhat, l2hat)
!!$          lhat = sqrt(l2hat)
!!$      
!!$          call khparams_from_fingering(w, lhat, hb, pr, db, m2, re, rm)
!!$      
!!$          call lmat(delta, m2, re, rm, k_star, n, l)
!!$          
!!$          if (present(withmode)) then
!!$             sigma = gamfromL(l, 2*n, .true., mode) 
!!$          else
!!$             sigma = gamfromL(l, 2*n)
!!$          end if
!!$      
!!$        end function sigma_from_fingering_params
!!$      
!!$        function omega_from_fingering_params(delta, w, hb, db, pr, tau, r0, k_star, n) result(omg)
!!$          real(dp), intent(in) :: delta, w, hb, db, pr, tau, r0, k_star
!!$          integer, intent(in) :: n
!!$          complex(dp) :: omg
!!$      
!!$          real(dp) :: lamhat, lhat, l2hat, m2, re, rm 
!!$          complex(dp) :: l(2*n,2*n)
!!$      
!!$          call gaml2max(pr, tau, r0, lamhat, l2hat)  
!!$          lhat = sqrt(l2hat)
!!$      
!!$          call khparams_from_fingering(w, lhat, hb, pr, db, m2, re, rm)
!!$      
!!$          call lmat(delta, m2, re, rm, k_star, n, l)
!!$          omg = omegafromL(l)
!!$      
!!$        end function omega_from_fingering_params
!!$      
!!$        function sigma_over_k_fingering_params(delta, w, hb, db, pr, tau, r0, k_stars, n) result(sigk)
!!$          real(dp), intent(in) :: delta, w, hb, db, pr, tau, r0 
!!$          real(dp), intent(in) :: k_stars(:)
!!$          integer, intent(in) :: n
!!$          real(dp) :: sigk(size(k_stars))
!!$      
!!$          integer :: i
!!$      
!!$          do i = 1, size(k_stars)
!!$             sigk(i) = sigma_from_fingering_params(delta, w, hb, db, pr, tau, r0, k_stars(i), n)
!!$          end do
!!$      
!!$        end function sigma_over_k_fingering_params
      
!!$        function gammax_kscan(delta, m2, re, rm, ks, n, ideal, badks_except, get_kmax) result(gammax)
!!$          ! Add get_kmax option
!!$      
!!$          real(dp), intent(in) :: delta, m2, re, rm 
!!$          real(dp), intent(in) :: ks(:)
!!$          integer, intent(in) :: n
!!$          logical, intent(in) :: ideal, badks_except, get_kmax
!!$          real(dp) :: gammax
!!$          real(dp), optional :: kmax
!!$      
!!$          ! Same as before but now optionally return kmax
!!$          
!!$          if (get_kmax) then
!!$             gammax = gammax_kscan(delta, m2, re, rm, ks, n, ideal, badks_except, kmax) 
!!$          else
!!$             gammax = gammax_kscan1(delta, m2, re, rm, ks, n, ideal, badks_except)
!!$          end if
!!$      
!!$        end function gammax_kscan
      
!!$        function sigma_max_kscan_fingering_params(delta, w, hb, db, pr, tau, r0, k_stars, n, &
!!$             badks_except, get_kmax) result(sigmax)
!!$      
!!$          ! Add get_kmax option
!!$      
!!$          real(dp), intent(in) :: delta, w, hb, db, pr, tau, r0
!!$          real(dp), intent(in) :: k_stars(:)
!!$          integer, intent(in) :: n
!!$          logical, intent(in) :: badks_except, get_kmax  
!!$          real(dp) :: sigmax
!!$          real(dp), optional :: kmax
!!$      
!!$          real(dp) :: lamhat, lhat, l2hat, m2, re, rm
!!$      
!!$          call gaml2max(pr, tau, r0, lamhat, l2hat)
!!$          lhat = sqrt(l2hat)
!!$      
!!$          call khparams_from_fingering(w, lhat, hb, pr, db, m2, re, rm)
!!$      
!!$          if (get_kmax) then
!!$             sigmax = gammax_kscan(delta, m2, re, rm, k_stars, n, .false., badks_except, kmax)
!!$          else 
!!$             sigmax = gammax_kscan(delta, m2, re, rm, k_stars, n, .false., badks_except)
!!$          end if
!!$      
!!$        end function sigma_max_kscan_fingering_params
      
        function gammax_minus_lambda(w, lamhat, lhat, hb, pr, db, delta, ks, n, ideal, badks_exception) result(f)
      
          ! Root finding helper function
      
          real(dp), intent(in) :: w, lamhat, lhat, hb, pr, db, delta    
          real(dp), intent(in) :: ks(:)
          integer, intent(in) :: n
          logical, intent(in) :: ideal, badks_exception
          real(dp) :: f
      
          real(dp) :: m2, re, rm
          real(dp) :: sigma
          
          m2 = hb / pow2(w)
          re = w / (pr * lhat) 
          rm = w / (db * lhat)

          sigma = gammax_kscan(delta, m2, re, rm, ks, n, ideal, badks_exception)

          f = sigma*w*lhat - CH*lamhat
      
        end function gammax_minus_lambda

        function gammax_minus_lambda_withTC(w, lamhat, lhat, hb, pr, tau, R0, db, ks, n, badks_exception, safety) result(f)
      
          ! Root finding helper function
      
          real(dp), intent(in) :: w, lamhat, lhat, hb, pr, tau, R0, db
          real(dp), intent(in) :: ks(:)
          integer, intent(in) :: n, safety
          logical, intent(in) :: badks_exception
          real(dp) :: gammax, f

          gammax = gammax_kscan_withTC(w, hb, db, pr, tau, R0, lamhat, lhat, ks, n, badks_exception, safety)
          f = gammax*C2 - lamhat

        end function gammax_minus_lambda_withTC
        
      end module kh_instability
