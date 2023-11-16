module parasite_model

    use fingering_modes, only: gaml2max
    use kh_instability, only: khparams_from_fingering, gammax_kscan, gammax_minus_lambda
    
    implicit none
  
    private
    public :: wf
    public :: hg19_eq32
    public :: deq32dw
    public :: wf_hg19
    public :: nuc
    public :: nut
    public :: gamma_tot
    public :: parasite_results
    public :: results_vs_r0
  
    real(dp), parameter :: kb = 1.24_dp
    real(dp), parameter :: CH = 1.66_dp
  
  contains
  
    function wf(pr, tau, r0, hb, db, ks, n, delta, ideal, badks_exception, get_kmax, &
         lamhat, l2hat) result(w)
  
      ! Root finding for wf
  
      real(dp), intent(in) :: pr, tau, r0, hb, db 
      real(dp), intent(in) :: ks(:)
      integer, intent(in) :: n
      real(dp), intent(in) :: delta  
      logical, intent(in) :: ideal, badks_exception, get_kmax
      real(dp), intent(in), optional :: lamhat, l2hat
      real(dp) :: w
      real(dp), optional :: kmax
  
      real(dp) :: lamhat_, l2hat_
      real(dp) :: lhat, w1, w2, wb(2)
      real(dp) :: args(10)
      integer :: ierr
  
      if (present(lamhat)) then
         lamhat_ = lamhat
         l2hat_ = l2hat
      else
         call gaml2max(pr, tau, r0, lamhat_, l2hat_) 
      end if
  
      lhat = sqrt(l2hat_)
  
      ! Set bracketing interval
      w1 = 2.0_dp*pi*lamhat_/lhat
      w2 = sqrt(2.0_dp*hb)
      wb = [w1, w1 + w2]
  
      ! Arguments to pass
      args = [lamhat_, lhat, hb, pr, db, delta, ks, n, ideal, badks_exception, CH]
  
      ! Call bracketed root finder
      call brent(gammax_minus_lambda, w, wb, args, ierr)
  
      if (present(get_kmax)) then
         call khparams_from_fingering(w, lhat, hb, pr, db, m2, re, rm)
         call gammax_kscan(delta, m2, re, rm, ks, n, ideal, badks_exception, kmax) 
      end if
  
    end function wf
  
    function hg19_eq32(w, pr, tau, r0, hb, CH) result(f)
  
      ! Evaluates Eq. 32 in HG19
  
      real(dp), intent(in) :: w, pr, tau, r0, hb, CH  
      real(dp) :: f
  
      real(dp) :: lamhat, lhat, l2hat
  
      call gaml2max(pr, tau, r0, lamhat, l2hat)
      lhat = sqrt(l2hat)
  
      f = 0.5_dp*w**2 - hb - (CH*lamhat/(0.42_dp*lhat))**1.5 * sqrt(w)
  
    end function hg19_eq32
  
    function deq32dw(w, pr, tau, r0, hb, CH) result(dfdw)
  
      ! Derivative wrt w of hg19_eq32
  
      real(dp), intent(in) :: w, pr, tau, r0, hb, CH
      real(dp) :: dfdw
  
      real(dp) :: lamhat, lhat, l2hat
  
      call gaml2max(pr, tau, r0, lamhat, l2hat)
      lhat = sqrt(l2hat)
  
      dfdw = w - 0.5_dp*(CH*lamhat/(0.42_dp*lhat))**1.5/sqrt(w)
  
    end function deq32dw
  
    function wf_hg19(pr, tau, r0, hb, CH) result(w)
  
      ! Newton solver for Eq. 32 in HG19
  
      real(dp), intent(in) :: pr, tau, r0, hb, CH
      real(dp) :: w
  
      real(dp) :: w0
      integer :: ierr
  
      ! Initial guess
      call gaml2max(pr, tau, r0, lamhat, l2hat)
      lhat = sqrt(l2hat)
      w0 = max(sqrt(2.0_dp*hb), 2.0_dp*pi*lamhat/lhat)
  
      ! Call Newton solver
      call newt(hg19_eq32, deq32dw, w, w0, pr, tau, r0, hb, CH, ierr) 
  
    end function wf_hg19
  
    function nuc(tau, w, lamhat, l2hat, kb) result(nu)
  
      ! Nusselt number for composition
  
      real(dp), intent(in) :: tau, w, lamhat, l2hat, kb
      real(dp) :: nu
  
      nu = 1.0_dp + kb*w**2/(tau*(lamhat + tau*l2hat))
  
    end function nuc
  
    function nut(w, lamhat, l2hat, kb) result(nu)
  
      ! Nusselt number for temperature
  
      real(dp), intent(in) :: w, lamhat, l2hat, kb  
      real(dp) :: nu
  
      nu = 1.0_dp + kb*w**2/(lamhat + l2hat)
  
    end function nut
  
    function gamma_tot(tau, r0, w, lamhat, l2hat, kb) result(gamma)
  
      ! Total growth rate
  
      real(dp), intent(in) :: tau, r0, w, lamhat, l2hat, kb
      real(dp) :: gamma
  
      real(dp) :: nut_, nuc_
  
      nut_ = nut(w, lamhat, l2hat, kb)
      nuc_ = nuc(tau, w, lamhat, l2hat, kb)
  
      gamma = r0*nut_/(tau*nuc_)
  
    end function gamma_tot
  
    function parasite_results(r0, hb, pr, tau, db, ks, n, lamhat, l2hat, eq32, double_n, &
         delta, ideal, CH, badks_exception) result(results)
  
      ! Get full set of results for parasite model
  
      real(dp), intent(in) :: r0, hb, pr, tau, db
      real(dp), intent(in) :: ks(:)
      integer, intent(in) :: n
      real(dp), intent(in) :: lamhat, l2hat 
      logical, intent(in) :: eq32, double_n, ideal, badks_exception
      real(dp), intent(in) :: delta, CH  
      type(dict) :: results
  
      real(dp) :: wf, fc, ft, nu_c, nu_t, gamma_tot
      real(dp) :: m2, re, kmax
  
      if (eq32) then
         wf = wf_hg19(pr, tau, r0, hb, CH)
      else if (double_n) then
         call wf(pr, tau, r0, hb, db, ks, 2*n-1, delta, ideal, badks_exception, .true., &
              lamhat, l2hat, wf, kmax)
      else
         call wf(pr, tau, r0, hb, db, ks, n, delta, ideal, badks_exception, .true., & 
              lamhat, l2hat, wf, kmax)
      end if
  
      fc = kb*wf**2/(r0*(lamhat + tau*l2hat))
      ft = kb*wf**2/(lamhat + l2hat)
      nu_c = nuc(tau, wf, lamhat, l2hat)
      nu_t = nut(wf, lamhat, l2hat)
      gamma_tot = gamma_tot(tau, r0, wf, lamhat, l2hat, kb)
  
      call khparams_from_fingering(wf, sqrt(l2hat), hb, pr, db, m2, re, rm)
  
      results = dict(names=["FC", "FT", "NuC", "NuT", "gammatot", "wf", "Re-star", "HB-star"])
  
      if (.not. eq32) then
         results%names(size(results%names)) = "kmax-star"
         results%vals(size(results%vals)) = kmax  
      end if
  
      results%vals = [fc, ft, nu_c, nu_t, gamma_tot, wf, re, m2]
  
    end function parasite_results
  
    function results_vs_r0(r0s, hb, pr, tau, db, ks, n, lamhats, l2hats, eq32, double_n, &
         delta, ideal, CH, badks_exception) result(results)
  
      ! Scan over density ratio
  
      real(dp), intent(in) :: r0s(:), hb, pr, tau, db  
      real(dp), intent(in) :: ks(:), lamhats(:), l2hats(:)
      integer, intent(in) :: n
      logical, intent(in) :: eq32, double_n, ideal, badks_exception
      real(dp), intent(in) :: delta, CH
      type(dict) :: results
  
      integer :: i
      type(dict) :: result_ri
  
      results = dict(names=["FC", "FT", "NuC", "NuT", "gammatot", "wf", "Re-star", "HB-star"])
  
      if (.not. eq32) then
         results%names(size(results%names)) = "kmax-star"
      end if
  
      do i = 1, size(r0s)
         if (eq32) then
            result_ri = parasite_results(r0s(i), hb, pr, tau, db, ks, n, &
                 lamhats(i), l2hats(i), eq32=eq32, CH=CH)
         else
            result_ri = parasite_results(r0s(i), hb, pr, tau, db, ks, n, &
                 lamhats(i), l2hats(i), eq32=eq32, double_n=double_n, &
                 delta=delta, ideal=ideal, CH=CH, badks_exception=badks_exception)
         end if
  
         results%vals(1:size(result_ri%vals)) = result_ri%vals
  
      end do
  
    end function results_vs_r0
  
  end module parasite_model