module parasite_model

    use const_def, only: dp, pi
    use num_lib
    use fingering_modes, only: gaml2max
    use kh_instability, only: khparams_from_fingering, gammax_kscan, gammax_minus_lambda
    
    implicit none
  
    private
    public :: wf
    public :: nuc
    public :: nut
    public :: gamma_tot
  
    real(dp), parameter :: kb = 1.24_dp
    real(dp), parameter :: CH = 1.66_dp
  
  contains
  
   ! Adrian said we can get rid of ideal, badks_exception, delta. 
   ! delta = Floquet wavenumber. delta = 0 means we are doing the analyisis for just one finger (n=1). For n > 1, delta = 1/n (e.g. 1/2 for 2 fingers). 
   ! Adrian says for most (all) cases delta = 0 is enough. We're keeping it there in case this turns out to be false in some edge case. 
   ! For now ideal = 1 or 0. Ideal = 1 means ideal gas. Be sure it is an integer, not logical 
   ! badks_exception true should throw a warning telling the kz search domain does not contain sigma_max, so the search didn't complete. Changed that in kh_instability.
   ! Note we will need to decide what MESA should use as input for ks (the grid of kz values over which we do the search). Need to be robust

    function wf(pr, tau, r0, hb, db, ks, n, delta, ideal, badks_exception, get_kmax, &
         lamhat, l2hat) result(w)
  
      ! Root finding for wf
  
      real(dp), intent(in) :: pr, tau, r0, hb, db 
      real(dp), intent(in) :: ks(:)
      integer, intent(in) :: n, ideal ! n MUST be odd
      real(dp), intent(in) :: delta  
      logical, intent(in) ::  badks_exception, get_kmax
      real(dp), intent(in), optional :: lamhat, l2hat
      real(dp) :: w
      ! real(dp), optional :: kmax
  
      integer, parameter  :: IMAX = 50
      real(dp), parameter :: EPSX = 2e-12_dp
      real(dp), parameter :: EPSY = 0._dp

      real(dp) :: lamhat_, l2hat_
      real(dp) :: lhat, w1, w2, y1, y2, hb_star, re, rm, dfdx, gammax
      real(dp), pointer :: rpar(:)
      integer, pointer  :: ipar(:)
      integer :: lrpar, lipar
      integer :: ierr
      
      if(mod(n,2) == 0) then
         stop "thermohaline parasite wf: n must be odd"
      end if
      
      if (present(lamhat)) then
         lamhat_ = lamhat
         l2hat_ = l2hat
      else
         call gaml2max(pr, tau, r0, lamhat_, l2hat_, ierr)
      end if
  
      lhat = sqrt(l2hat_)
  
      ! Set bracketing interval
      w1 = 2.0_dp*pi*lamhat_/lhat
      w2 = w1 + sqrt(2.0_dp*hb)
      
      lrpar = 7 + size(ks)
      lipar = 1
      allocate(rpar(lrpar))
      allocate(ipar(lipar))
      
      ! Arguments to pass
      rpar = [lamhat_, lhat, hb, pr, db, delta, CH, ks]
      ipar = [n]
      dfdx = 0d0 ! unused for bracket search
      
      y1 = gx_m_lam(w1, dfdx, lrpar, rpar, lipar, ipar, ierr)
      y2 = gx_m_lam(w2, dfdx, lrpar, rpar, lipar, ipar, ierr)

      ! Call bracketed root finder
      w = safe_root_with_brackets( &
           gx_m_lam, w1, w2, y1, y2, imax, epsx, epsy, lrpar, rpar, lipar, ipar, ierr)

      ! needs more work on gammax_kscan
      if (get_kmax) then
         call khparams_from_fingering(w, lhat, hb, pr, db, hb_star, re, rm)
         gammax = gammax_kscan(delta, hb_star, re, rm, ks, n, .false., badks_exception) !, kmax) 
      end if

      deallocate(rpar)
      deallocate(ipar)
  
    end function wf

    ! wrapper for gammax_minus_lambda for use with safe_root_with_brackets
    real(dp) function gx_m_lam(x, dfdx, lrpar, rpar, lipar, ipar, ierr) result(f)
      
      real(dp), intent(in)             :: x
      real(dp), intent(out)            :: dfdx
      integer, intent(in)              :: lrpar
      real(dp), intent(inout), pointer :: rpar(:)
      integer, intent(in)              :: lipar
      integer, intent(inout), pointer  :: ipar(:)
      integer, intent(out)             :: ierr

      dfdx = 0d0 ! unused for bracket search
      
      if (lipar /= 1) then
         ierr = -1
         write(*,*) "lrpar, lipar", lrpar, lipar
         stop "wrong number of parameters for bracketed root solve"
      end if
      
      associate( &
           w => x, &
           lamhat => rpar(1), &
           lhat => rpar(2), &
           hb => rpar(3), &
           pr => rpar(4), &
           db => rpar(5), &
           delta => rpar(6), &
           CH => rpar(7), &
           ks => rpar(8:), & 
           n => ipar(1))

        f = gammax_minus_lambda(w, lamhat, lhat, hb, pr, db, delta, ks, n, .false., .false.)
        ! ideal = .false., and ignoring badks_exception right now
        
      end associate

      ierr = 0

      return
      
    end function gx_m_lam
  
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
    
  end module parasite_model
