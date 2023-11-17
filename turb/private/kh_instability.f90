module kh_instability

    use fingering_modes, only: gaml2max
  
    
    implicit none
  
    private
    public :: khparams_from_fingering
    public :: deln
    public :: lmat
    public :: gamfromL
    public :: omegafromL
    public :: gamma_over_k
    public :: omega_over_k
    public :: gammax_kscan
    public :: gammax_minus_lambda
  
    integer, parameter :: dp = kind(1.d0)
    real(dp), parameter :: CH = 1.66_dp

  
  contains
  
    subroutine khparams_from_fingering(w, lhat, hb, pr, db, hb_star, re, rm)
      real(dp), intent(in) :: w, lhat, hb, pr, db 
      real(dp), intent(out) :: hb_star, re, rm
  
      hb_star = hb / w**2
      re = w / (pr * lhat)
      rm = w / (db * lhat)
  
    end subroutine khparams_from_fingering
  
    function deln(k, n, delta)
      real(dp), intent(in) :: k  
      integer, intent(in) :: n  
      real(dp), intent(in) :: delta
      real(dp) :: deln
  
      deln = k**2 + (n + delta)**2
  
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
               l(i,i) = (0.0_dp, 1.0_dp) * (diss/re) * deltan
            else
               l(i,i) = (0.0_dp, 1.0_dp) * diss/rm * deltan
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
      do i = 3, 2*n-3
         deltan = delns_m(i)
         deltanp1 = delns_m(i+2)
         deltanm1 = delns_m(i-2)
         if (mod(ms(i),2) == 0) then
            l(i,i+2) = -k * (1 - deltanp1) / (2.0_dp*deltan)
            l(i,i-2) = k * (1 - deltanm1) / (2.0_dp*deltan)
         else
            l(i,i+2) = k/2.0_dp 
            l(i,i-2) = -k/2.0_dp
         end if
      end do

  
      ! Boundary conditions
      i = 1
      l(i,i) = (0.0_dp,1.0_dp) * delns_m(i) * diss / re
      l(i,i+1) = m2 * k
      l(i,i+2) = -k * (1 - delns_m(i+2)) / (2.0_dp*delns_m(i))
  
      i = 2
      l(i,i) = (0.0_dp,1.0_dp) * delns_m(i) * diss / rm  
      l(i,i-1) = k
      l(i,i+2) = k/2.0_dp
  
      i = 2*n-1
      l(i,i) = (0.0_dp,1.0_dp) * delns_m(i) * diss / re
      l(i,i-1) = m2*k
      l(i,i-2) = k * (1 - delns_m(i-2)) / (2.0_dp*delns_m(i))
  
      i = 2*n
      l(i,i) = (0.0_dp,1.0_dp) * delns_m(i) * diss / rm
      l(i,i-2) = k
      l(i,i-3) = -k/2.0_dp
  
    end subroutine lmat
  
    function gamfromL(l, withmode) result(gam)
      complex(dp), intent(in) :: l(:,:)
      logical, intent(in), optional :: withmode
      real(dp) :: gam 
      !complex(dp), allocatable, optional :: mode(:)
  
      complex(dp) :: w(size(l,1))
      complex(dp), allocatable :: v(:,:)
      integer :: i
  
      call eigensystem(l, w, v)
  
!!$      if (present(withmode)) then
!!$         i = maxloc(aimag(w),1)
!!$         gam = -aimag(w(i))
!!$         allocate(mode(size(w)))
!!$         mode = v(:,i)
!!$      else
         gam = maxval(-aimag(w))
!!$      end if
  
    end function gamfromL
  
    function omegafromL(l) result(omg)
      complex(dp), intent(in) :: l(:,:)
      complex(dp) :: omg
  
      complex(dp) :: w(size(l,1)) 
      complex(dp), allocatable :: v(:,:)
      integer :: i
  
      call eigensystem(l, w, v)
  
      omg = w(maxloc(aimag(w)))
  
    end function omegafromL
  
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
         gamk(i) = gamfromL(l_result)
      end do
  
    end function gamma_over_k
  
    function omega_over_k(delta, m2, re, rm, ks, n, ideal) result(omgk)
      real(dp), intent(in) :: delta, m2, re, rm 
      real(dp), intent(in) :: ks(:)
      integer, intent(in) :: n
      logical, intent(in) :: ideal
  
      complex(dp) :: l_result(2*n,2*n)
      complex(dp) :: omgk(size(ks))
      integer :: i
  
      do i = 1, size(ks)
         call lmat(delta, m2, re, rm, ks(i), n, ideal, l_result)
         omgk(i) = omegafromL(l_result)
      end do
  
    end function omega_over_k
  
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
!!$             sigma = gamfromL(l, .true., mode) 
!!$          else
!!$             sigma = gamfromL(l)
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
          
          m2 = hb / w**2
          re = w / (pr * lhat) 
          rm = w / (db * lhat)
      
          sigma = gammax_kscan(delta, m2, re, rm, ks, n, ideal, badks_exception)
      
          f = sigma*w*lhat - CH*lamhat
      
        end function gammax_minus_lambda
      
      end module kh_instability
