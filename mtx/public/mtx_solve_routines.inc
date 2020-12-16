
      
      subroutine my_getf2(m, a, lda, ipiv, info)
         integer :: info, lda, m
         integer :: ipiv(:)
         real(dp) :: a(:,:)
         real(dp), parameter :: one=1, zero=0
         integer :: i, j, jp, ii, jj, n, mm
         real(dp) :: tmp, da
         do j = 1, m
            info = 0
            jp = j - 1 + maxloc(abs(a(j:lda,j)),dim=1)
            ipiv( j ) = jp
            if( a( jp, j ).ne.zero ) then
               if( jp.ne.j ) then ! swap a(j,:) and a(jp,:)
                  !$omp simd private(tmp)
                  do i=1,m
                     tmp = a(j,i)
                     a(j,i) = a(jp,i)
                     a(jp,i) = tmp
                  end do
               end if
               if( j.lt.m ) then 
                  da = one / a( j, j )
                  !$omp simd
                  do i = 1, m-j 
                     a( j+i, j ) = da*a( j+i, j )
                  end do 
               end if 
            else if( info.eq.0 ) then
               info = j
            end if
            if( j.lt.m ) then
               !call dger( m-j, m-j, -one, a( j+1, j ), 1, a( j, j+1 ), lda, a( j+1, j+1 ), lda )
               do jj = j+1, m
                  !$omp simd
                  do ii = j+1, m
                     a(ii,jj) = a(ii,jj) - a(ii,j)*a(j,jj)
                  end do
               end do
            end if
         end do
      end subroutine my_getf2

      
      subroutine my_getf2_n4(a, ipiv, info)
         integer :: info
         integer :: ipiv(:)
         real(dp) :: a(:,:)
         real(dp), parameter :: one=1, zero=0
         integer, parameter :: m=4, lda=4
         integer :: i, j, jp, ii, jj, n, mm
         real(dp) :: tmp, da
         do j = 1, m
            info = 0
            jp = j - 1 + maxloc(abs(a(j:lda,j)),dim=1)
            ipiv( j ) = jp
            if( a( jp, j ).ne.zero ) then
               if( jp.ne.j ) then ! swap a(j,:) and a(jp,:)
                  !$omp simd private(tmp)
                  do i=1,m
                     tmp = a(j,i)
                     a(j,i) = a(jp,i)
                     a(jp,i) = tmp
                  end do
               end if
               if( j.lt.m ) then 
                  da = one / a( j, j )
                  !$omp simd
                  do i = 1, m-j 
                     a( j+i, j ) = da*a( j+i, j )
                  end do 
               end if 
            else if( info.eq.0 ) then
               info = j
            end if
            if( j.lt.m ) then
               !call dger( m-j, m-j, -one, a( j+1, j ), 1, a( j, j+1 ), lda, a( j+1, j+1 ), lda )
               do jj = j+1, m
                  !$omp simd
                  do ii = j+1, m
                     a(ii,jj) = a(ii,jj) - a(ii,j)*a(j,jj)
                  end do
               end do
            end if
         end do
      end subroutine my_getf2_n4

      
      subroutine my_getf2_n5(a, ipiv, info)
         integer :: info
         integer :: ipiv(:)
         real(dp) :: a(:,:)
         real(dp), parameter :: one=1, zero=0
         integer, parameter :: m=5, lda=5
         integer :: i, j, jp, ii, jj, n, mm
         real(dp) :: tmp, da
         do j = 1, m
            info = 0
            jp = j - 1 + maxloc(abs(a(j:lda,j)),dim=1)
            ipiv( j ) = jp
            if( a( jp, j ).ne.zero ) then
               if( jp.ne.j ) then ! swap a(j,:) and a(jp,:)
                  !$omp simd private(tmp)
                  do i=1,m
                     tmp = a(j,i)
                     a(j,i) = a(jp,i)
                     a(jp,i) = tmp
                  end do
               end if
               if( j.lt.m ) then 
                  da = one / a( j, j )
                  !$omp simd
                  do i = 1, m-j 
                     a( j+i, j ) = da*a( j+i, j )
                  end do 
               end if 
            else if( info.eq.0 ) then
               info = j
            end if
            if( j.lt.m ) then
               !call dger( m-j, m-j, -one, a( j+1, j ), 1, a( j, j+1 ), lda, a( j+1, j+1 ), lda )
               do jj = j+1, m
                  !$omp simd
                  do ii = j+1, m
                     a(ii,jj) = a(ii,jj) - a(ii,j)*a(j,jj)
                  end do
               end do
            end if
         end do
      end subroutine my_getf2_n5
            
      
      subroutine my_getrs( n, nrhs, a, lda, ipiv, b, ldb, info )
         integer :: info, lda, ldb, n, nrhs
         integer :: ipiv(:)
         real(dp) :: a(:,:), b(:,:) ! a( lda, * ), b( ldb, * )
         real(dp), parameter :: one=1, zero=0
         real(dp) :: temp
         integer :: i,j,k, n32, ix, ip
         info = 0
         call my_laswp(nrhs, b, ldb, 1, n, ipiv, 1 )
         !call dtrsm( 'left', 'lower', 'no transpose', 'unit', n, nrhs, one, a, lda, b, ldb )
         do j = 1,nrhs
            do k = 1,n
               if (b(k,j).ne.zero) then
                  !$omp simd
                  do i = k + 1,n
                     b(i,j) = b(i,j) - b(k,j)*a(i,k)
                  end do
               end if
            end do
         end do
         !call dtrsm( 'left', 'upper', 'no transpose', 'non-unit', n, nrhs, one, a, lda, b, ldb )
         do j = 1,nrhs
            do k = n,1,-1
               if (b(k,j).ne.zero) then
                  b(k,j) = b(k,j)/a(k,k)
                  !$omp simd
                  do i = 1,k - 1
                     b(i,j) = b(i,j) - b(k,j)*a(i,k)
                  end do
               end if
            end do
         end do
         
      end subroutine my_getrs
      
      
      subroutine my_laswp( n,   a, lda,  k1, k2, ipiv,  incx )
         integer :: incx, k1, k2, lda, n
         integer :: ipiv(:)
         real(dp) :: a(:,:) ! a( lda, * )
         integer :: i, i1, i2, inc, ip, ix, ix0, j, k, n32
         real(dp) :: temp
         ! interchange row i with row ipiv(i) for each of rows k1 through k2.
         if( incx.gt.0 ) then
            ix0 = k1
            i1 = k1
            i2 = k2
            inc = 1
         else if( incx.lt.0 ) then
            ix0 = 1 + ( 1-k2 )*incx
            i1 = k2
            i2 = k1
            inc = -1
         else
            return
         end if
         ix = ix0
         do i = i1, i2, inc
            ip = ipiv( ix )
            if( ip.ne.i ) then
               !$omp simd private(temp)
               do k = 1, n
                  temp = a( i, k )
                  a( i, k ) = a( ip, k )
                  a( ip, k ) = temp
               end do
            end if
            ix = ix + incx
         end do         
      end subroutine my_laswp
            
      
      subroutine my_getrs1( n, a, lda, ipiv, b, ldb, info )
         integer :: info, lda, ldb, n
         integer :: ipiv(:)
         real(dp) :: a(:,:), b(:)
         real(dp), parameter :: one=1, zero=0
         real(dp) :: temp
         integer :: i,k
         info = 0
         !$omp simd private(temp)
         do i = 1, n
            temp = b(i)
            b(i) = b(ipiv(i))
            b(ipiv(i)) = temp
         end do         
         do k = 1,n
            if (b(k).ne.zero) then
               !$omp simd
               do i = k + 1,n
                  b(i) = b(i) - b(k)*a(i,k)
               end do
            end if
         end do
         do k = n,1,-1
            if (b(k).ne.zero) then
               b(k) = b(k)/a(k,k)
               !$omp simd
               do i = 1,k - 1
                  b(i) = b(i) - b(k)*a(i,k)
               end do
            end if
         end do
      end subroutine my_getrs1
            
      
      subroutine my_getrs1_n4( a, ipiv, b, info )
         integer :: info
         integer :: ipiv(:)
         real(dp) :: a(:,:), b(:)
         real(dp), parameter :: one=1, zero=0
         integer, parameter :: n = 4
         real(dp) :: temp
         integer :: i,k
         info = 0
         do i = 1, n
            temp = b(i)
            b(i) = b(ipiv(i))
            b(ipiv(i)) = temp
         end do         
         do k = 1,n
            if (b(k).ne.zero) then
               !$omp simd
               do i = k + 1,n
                  b(i) = b(i) - b(k)*a(i,k)
               end do
            end if
         end do
         do k = n,1,-1
            if (b(k).ne.zero) then
               b(k) = b(k)/a(k,k)
               !$omp simd
               do i = 1,k - 1
                  b(i) = b(i) - b(k)*a(i,k)
               end do
            end if
         end do
      end subroutine my_getrs1_n4
            
      
      subroutine my_getrs1_n5( a, ipiv, b, info )
         integer :: info
         integer :: ipiv(:)
         real(dp) :: a(:,:), b(:)
         real(dp), parameter :: one=1, zero=0
         integer, parameter :: n = 5
         real(dp) :: temp
         integer :: i,k
         info = 0
         do i = 1, n
            temp = b(i)
            b(i) = b(ipiv(i))
            b(ipiv(i)) = temp
         end do         
         do k = 1,n
            if (b(k).ne.zero) then
               !$omp simd
               do i = k + 1,n
                  b(i) = b(i) - b(k)*a(i,k)
               end do
            end if
         end do
         do k = n,1,-1
            if (b(k).ne.zero) then
               b(k) = b(k)/a(k,k)
               !$omp simd
               do i = 1,k - 1
                  b(i) = b(i) - b(k)*a(i,k)
               end do
            end if
         end do
      end subroutine my_getrs1_n5
      
      
      subroutine my_gemm0_p1(m,n,k,a,lda,b,ldb,c,ldc) ! c := -a*b
         integer, intent(in) :: k,lda,ldb,ldc,m,n
         real(dp), dimension(:,:) :: a, b, c ! a(lda,*),b(ldb,*),c(ldc,*)
         integer :: j, i
         real(dp), parameter :: zero=0
         include 'formats.dek'
         ! transa = 'n'
         ! transb = 'n'
         ! alpha = -1
         ! beta = 0
         ! assumes other args are valid
         c(1:m,1:n) = zero
         call my_gemm_p1(m,n,k,a,lda,b,ldb,c,ldc)
      end subroutine my_gemm0_p1
      
      
      subroutine my_gemm_plus_mm(m,n,k,a,b,d,e,c) ! c := c + a*b + d*e
         integer, intent(in) :: k,m,n
         real(dp), dimension(:,:) :: a, b, c, d, e
         real(dp) :: tmp_b, tmp_e
         real(dp), parameter :: zero=0
         integer :: j, i, l
         do j = 1,n
            do l = 1,k
               tmp_b = b(l,j)
               tmp_e = e(l,j)
               if (tmp_b .ne. zero) then
                  if (tmp_e .ne. zero) then
                     !$omp simd
                     do i = 1,m
                        c(i,j) = c(i,j) + tmp_b*a(i,l) + tmp_e*d(i,l)
                     end do
                  else
                     !$omp simd
                     do i = 1,m
                        c(i,j) = c(i,j) + tmp_b*a(i,l)
                     end do
                  end if
               else if (tmp_e .ne. zero) then
                  !$omp simd
                  do i = 1,m
                     c(i,j) = c(i,j) + tmp_e*d(i,l)
                  end do
               end if
            end do
         end do      
      end subroutine my_gemm_plus_mm
      
      
      subroutine my_gemm(m,n,k,a,lda,b,ldb,c,ldc) ! c := c - a*b
         integer, intent(in) :: k,lda,ldb,ldc,m,n
         real(dp), dimension(:,:) :: a, b, c ! a(lda,*),b(ldb,*),c(ldc,*)
         real(dp) :: tmp
         real(dp), parameter :: zero=0
         integer :: j, i, l
         ! transa = 'n'
         ! transb = 'n'
         ! alpha = -1
         ! beta = 1
         ! assumes other args are valid
         do j = 1,n
            do l = 1,k
               tmp = b(l,j)
               if (tmp .ne. zero) then
                  !$omp simd
                  do i = 1,m
                     c(i,j) = c(i,j) - tmp*a(i,l)
                  end do
               end if
            end do
         end do      
      end subroutine my_gemm
      
      
      subroutine my_gemm_p1(m,n,k,a,lda,b,ldb,c,ldc) ! c := c + a*b
         integer, intent(in) :: k,lda,ldb,ldc,m,n
         real(dp), dimension(:,:) :: a, b, c ! a(lda,*),b(ldb,*),c(ldc,*)
         real(dp) :: tmp
         real(dp), parameter :: zero=0
         integer :: j, i, l
         ! transa = 'n'
         ! transb = 'n'
         ! alpha = 1
         ! beta = 1
         ! assumes other args are valid
         do j = 1,n
            do l = 1,k
               tmp = b(l,j)
               if (tmp .ne. zero) then
                  !$omp simd
                  do i = 1,m
                     c(i,j) = c(i,j) + tmp*a(i,l)
                  end do
               end if
            end do
         end do      
      end subroutine my_gemm_p1


      subroutine my_gemv_mv(m,n,a,x,b,z,y) ! y = y - a*x - b*z
         integer lda,m,n
         real(dp) :: a(:,:), b(:,:)
         real(dp) :: x(:), z(:), y(:)
         real(dp) :: tmp_x, tmp_z
         real(dp), parameter :: zero=0
         integer :: j, i
         do j = 1,n
            tmp_x = x(j)
            tmp_z = z(j)
            if (tmp_x.ne.zero) then
               if (tmp_z /= zero) then
                  !$omp simd
                  do i = 1,m
                     y(i) = y(i) - tmp_x*a(i,j) - tmp_z*b(i,j)
                  end do
               else
                  !$omp simd
                  do i = 1,m
                     y(i) = y(i) - tmp_x*a(i,j)
                  end do
               end if
            else if (tmp_z /= zero) then
               !$omp simd
               do i = 1,m
                  y(i) = y(i) - tmp_z*b(i,j)
               end do
            end if
         end do
      end subroutine my_gemv_mv


      subroutine my_gemv(m,n,a,lda,x,y) ! y = y - a*x
         integer lda,m,n
         real(dp) :: a(:,:) ! (lda,*)
         real(dp) :: x(:), y(:)
         real(dp) :: tmp
         real(dp), parameter :: zero=0
         ! trans = 'n'
         ! alpha = -1
         ! beta = 1
         ! incx = 1
         ! incy = 1
         integer :: j, i
         do j = 1,n
            tmp = x(j)
            if (tmp.ne.zero) then
               !$omp simd
               do i = 1,m
                  y(i) = y(i) - tmp*a(i,j)
               end do
            end if
         end do
      end subroutine my_gemv


      subroutine my_gemv_p_mv(m,n,a,x,b,z,y) ! y = y + a*x + b*z
         integer lda,m,n
         real(dp) :: a(:,:), b(:,:)
         real(dp) :: x(:), z(:), y(:)
         real(dp) :: tmp_x, tmp_z
         real(dp), parameter :: zero=0
         integer :: j, i
         do j = 1,n
            tmp_x = x(j)
            tmp_z = z(j)
            if (tmp_x.ne.zero) then
               if (tmp_z /= zero) then
                  !$omp simd
                  do i = 1,m
                     y(i) = y(i) + tmp_x*a(i,j) + tmp_z*b(i,j)
                  end do
               else
                  !$omp simd
                  do i = 1,m
                     y(i) = y(i) + tmp_x*a(i,j)
                  end do
               end if
            else if (tmp_z /= zero) then
               !$omp simd
               do i = 1,m
                  y(i) = y(i) + tmp_z*b(i,j)
               end do
            end if
         end do
      end subroutine my_gemv_p_mv


      subroutine my_gemv_p1(m,n,a,lda,x,y) ! y = y + a*x
         integer lda,m,n
         real(dp) :: a(:,:) ! (lda,*)
         real(dp) :: x(:), y(:)
         real(dp) :: tmp
         real(dp), parameter :: zero=0
         ! trans = 'n'
         ! alpha = -1
         ! beta = 1
         ! incx = 1
         ! incy = 1
         integer :: j, i
         do j = 1,n
            tmp = x(j)
            if (tmp.ne.zero) then
               !$omp simd
               do i = 1,m
                  y(i) = y(i) + tmp*a(i,j)
               end do
            end if
         end do
      end subroutine my_gemv_p1

      
      subroutine my_getf2_no_pivot(m, a, lda, info)
         integer :: info, lda, m
         real(dp) :: a(:,:)
         real(dp), parameter :: one=1, zero=0
         integer :: i, j, ii, jj, n, mm
         real(dp) :: tmp, da
         do j = 1, m
            info = 0
            if( a( j, j ).ne.zero ) then
               if( j.lt.m ) then 
                  da = one / a( j, j )
                  !$omp simd
                  do i = 1, m-j 
                     a( j+i, j ) = da*a( j+i, j )
                  end do 
               end if 
            else if( info.eq.0 ) then
               info = j
            end if
            if( j.lt.m ) then
               do jj = j+1, m
                  !$omp simd
                  do ii = j+1, m
                     a(ii,jj) = a(ii,jj) - a(ii,j)*a(j,jj)
                  end do
               end do
            end if
         end do
      end subroutine my_getf2_no_pivot
            
      
      subroutine my_getrs_no_pivot( n, nrhs, a, lda, b, ldb, info )
         integer :: info, lda, ldb, n, nrhs
         real(dp) :: a(:,:), b(:,:) ! a( lda, * ), b( ldb, * )
         real(dp), parameter :: one=1, zero=0
         real(dp) :: temp
         integer :: i,j,k, n32, ix, ip
         info = 0
         do j = 1,nrhs
            do k = 1,n
               if (b(k,j).ne.zero) then
                  !$omp simd
                  do i = k + 1,n
                     b(i,j) = b(i,j) - b(k,j)*a(i,k)
                  end do
               end if
            end do
         end do
         do j = 1,nrhs
            do k = n,1,-1
               if (b(k,j).ne.zero) then
                  b(k,j) = b(k,j)/a(k,k)
                  !$omp simd
                  do i = 1,k - 1
                     b(i,j) = b(i,j) - b(k,j)*a(i,k)
                  end do
               end if
            end do
         end do
      end subroutine my_getrs_no_pivot
            
      
      subroutine my_getrs1_no_pivot( n, a, lda, b, ldb, info )
         integer :: info, lda, ldb, n
         real(dp) :: a(:,:), b(:)
         real(dp), parameter :: one=1, zero=0
         real(dp) :: temp
         integer :: i,k
         info = 0
         do k = 1,n
            if (b(k).ne.zero) then
               !$omp simd
               do i = k + 1,n
                  b(i) = b(i) - b(k)*a(i,k)
               end do
            end if
         end do
         do k = n,1,-1
            if (b(k).ne.zero) then
               b(k) = b(k)/a(k,k)
               !$omp simd
               do i = 1,k - 1
                  b(i) = b(i) - b(k)*a(i,k)
               end do
            end if
         end do
      end subroutine my_getrs1_no_pivot



      subroutine qrdcmp(a,n,c,d,sing)
      ! constructs the qr decomposition of a(1:n,1:n).
      ! the upper triangular matrix r is returned in the upper triangle of a, except
      ! for the diagonal which is returned in d(1:n). the orthogonal matrix q is
      ! represented as a product of n-1 householder matrices q_1....q_{n-1}, where
      ! q_j = i - u_j * u_j/c_j. the ith component of u_j is zero for i=1..j-1,
      ! while the nonzero components are returned in a(i,j) for i=j....n. sing
      ! returns as true if singularity is encountered during the decomposition,
      ! but the decomposition is still completed.
      integer, intent(in) :: n
      real(dp), intent(inout) :: a(:,:)
      real(dp), intent(inout) :: c(:),d(:)
      logical, intent(out) :: sing
      integer :: i,j,k
      real(dp) :: scale,sigma,sum,tau,dum
      sing = .false.
      do k=1,n-1
       scale = 0.0d0
       do i=k,n
        scale = max(scale,abs(a(i,k)))
       enddo
       if (scale .eq. 0.0) then
        sing = .true.
        c(k) = 0.0d0
        d(k) = 0.0d0
       else
        dum = 1.0d0/scale
        do i=k,n
         a(i,k) = a(i,k) * dum
        enddo
        sum = 0.0d0
        do i=k,n
         sum = sum + a(i,k)*a(i,k)
        enddo
        sigma  = sign(sqrt(sum),a(k,k))
        a(k,k) = a(k,k) + sigma
        c(k)   = sigma * a(k,k)
        d(k)   = -scale * sigma
        do j=k+1,n
         sum = 0.0d0
         do i=k,n
          sum = sum + a(i,k)*a(i,j)
         enddo
         tau = sum/c(k)
         do i=k,n
          a(i,j) = a(i,j) - tau * a(i,k)
         enddo
        enddo
       end if
      enddo
      d(n) = a(n,n)
      if (d(n) .eq. 0.0) sing = .true.
      end subroutine qrdcmp


      subroutine qrsolv(a,n,c,d,b)
      ! solves the linear system a*x=b, given the qr decomposition from
      ! routine qrdcmp. b is input as the right hand side and returns with the
      ! solution x.
      integer, intent(in) :: n
      real(dp), intent(in) :: a(:,:),c(:),d(:)
      real(dp), intent(inout) :: b(:)
      integer :: i,j
      real(dp) :: sum,tau
      do j=1,n-1
       sum = 0.0d0
       do i=j,n
        sum = sum + a(i,j)*b(i)
       enddo
       tau = sum/c(j)
       do i=j,n
        b(i) = b(i) - tau*a(i,j)
       enddo
      enddo
      b(n) = b(n)/d(n)
      do i=n-1,1,-1
       sum = 0.0d0
       do j=i+1,n
        sum = sum + a(i,j)*b(j)
       enddo
       b(i) = (b(i) - sum)/d(i)
      enddo
      end subroutine qrsolv

