module fingering_modes

    implicit none
  
    private
    public :: rfromR
    public :: lamguess
    public :: k2guess
    public :: eq1
    public :: eq2
    public :: fun
    public :: jac
    public :: gaml2max
  
  contains
  
    function rfromR(r0, tau) result(r)
      real(dp), intent(in) :: r0, tau
      real(dp) :: r
  
      r = (r0 - 1.0_dp) / (-1.0_dp + 1.0_dp/tau)
  
    end function rfromR
  
    function lamguess(pr, tau, r0) result(lam)
      real(dp), intent(in) :: pr, tau, r0
      real(dp) :: lam
      real(dp) :: r
  
      r = rfromR(r0, tau)
      if (r < tau) then
         lam = sqrt(pr) - pr*sqrt(1.0_dp + tau/pr) 
      else
         if (r > 0.5_dp) then
            lam = 2.0_dp*pr*(tau/pr)*((1.0_dp/3.0_dp)*(1.0_dp - r))**1.5_dp &
                 / (1.0_dp - (1.0_dp - r)*(1.0_dp + tau/pr)/3.0_dp)
         else
            lam = sqrt(pr*tau/r) - pr*sqrt(1.0_dp + tau/pr)
         end if
      end if
  
    end function lamguess
  
    function k2guess(pr, tau, r0) result(k2)
      real(dp), intent(in) :: pr, tau, r0
      real(dp) :: k2
      real(dp) :: r
  
      r = rfromR(r0, tau)
      if (r < tau) then
         k2 = (1.0_dp + tau/pr)**(-0.5_dp) - sqrt(pr) * (1.0_dp + (tau/pr)*(1.0_dp + tau/pr)**(-2.0_dp)) 
      else
         if (r > 0.5_dp) then
            k2 = sqrt((1.0_dp - r)/3.0_dp)
         else
            k2 = (1.0_dp + tau/pr)**(-0.5_dp) - 2.0_dp*sqrt(r*tau/pr) * (1.0_dp + tau/pr)**(-2.5_dp)
         end if
      end if
  
    end function k2guess
  
    function eq1(lam, k2, pr, tau, r0) result(eq)
      real(dp), intent(in) :: lam, k2, pr, tau, r0
      real(dp) :: eq
  
      real(dp) :: b2, b1, b0
  
      b2 = k2*(1.0_dp + pr + tau) 
      b1 = k2**2*(tau*pr + pr + tau) + pr*(1.0_dp - 1.0_dp/r0)
      b0 = k2**3*tau*pr + k2*pr*(tau - 1.0_dp/r0)
  
      eq = lam**3 + b2*lam**2 + b1*lam + b0
  
    end function eq1
  
    function eq2(lam, k2, pr, tau, r0) result(eq)
      real(dp), intent(in) :: lam, k2, pr, tau, r0  
      real(dp) :: eq
  
      real(dp) :: c2, c1, c0
  
      c2 = 1.0_dp + pr + tau
      c1 = 2.0_dp*k2*(tau*pr + tau + pr) 
      c0 = 3.0_dp*k2**2*tau*pr + pr*(tau - 1.0_dp/r0)
  
      eq = c2*lam**2 + c1*lam + c0
  
    end function eq2
  
    function fun(x, pr, tau, r0, passk1) result(f)
      real(dp), intent(in) :: x(2), pr, tau, r0
      logical, intent(in), optional :: passk1
      real(dp) :: f(2)
  
      if (present(passk1)) then
         f(1) = eq1(x(1), x(2)**2, pr, tau, r0)  
         f(2) = eq2(x(1), x(2)**2, pr, tau, r0)
      else
         f(1) = eq1(x(1), x(2), pr, tau, r0)
         f(2) = eq2(x(1), x(2), pr, tau, r0)
      end if
  
    end function fun
  
    function jac(x, pr, tau, r0, passk1) result(dfdx)
      real(dp), intent(in) :: x(2), pr, tau, r0
      logical, intent(in), optional :: passk1
      real(dp) :: dfdx(2,2)
  
      real(dp) :: lam, k2
      real(dp) :: b2, db2dk2, b1, db1dk2, b0, db0dk2
      real(dp) :: c2, c1, dc1dk2, c0, dc0dk2
  
      lam = x(1)
      if (present(passk1)) then
         k2 = x(2)**2
      else  
         k2 = x(2)
      end if
  
      b2 = k2*(1.0_dp + pr + tau)
      db2dk2 = 1.0_dp + pr + tau
  
      b1 = k2**2*(tau*pr + pr + tau) + pr*(1.0_dp - 1.0_dp/r0) 
      db1dk2 = 2.0_dp*k2*(tau*pr + pr + tau)
  
      b0 = k2**3*tau*pr + k2*pr*(tau - 1.0_dp/r0)
      db0dk2 = 3.0_dp*k2**2*tau*pr + pr*(tau - 1.0_dp/r0)
  
      dfdx(1,1) = 3.0_dp*lam**2 + 2.0_dp*b2*lam + b1
      dfdx(1,2) = lam**2*db2dk2 + lam*db1dk2 + db0dk2
      if (present(passk1)) dfdx(1,2) = dfdx(1,2)*2.0_dp*x(2)
  
      c2 = 1.0_dp + pr + tau
      c1 = 2.0_dp*k2*(tau*pr + tau + pr)
      dc1dk2 = c1/k2
      c0 = 3.0_dp*k2**2*tau*pr + pr*(tau - 1.0_dp/r0)
      dc0dk2 = 6.0_dp*k2*tau*pr
  
      dfdx(2,1) = 2.0_dp*c2*lam + c1
      dfdx(2,2) = lam*dc1dk2 + dc0dk2
      if (present(passk1)) dfdx(2,2) = dfdx(2,2)*2.0_dp*x(2)
  
    end function jac
  
    function gaml2max(pr, tau, r0) result(x)
      real(dp), intent(in) :: pr, tau, r0
      real(dp) :: x(2)
  
      real(dp) :: lam, k2
      real(dp) :: lamg, k2g  
      real(dp) :: f(2)
  
      lamg = lamguess(pr, tau, r0)
      k2g = k2guess(pr, tau, r0)
  
      ! First solve for lam, k2
      call hybrd1(fun, x, f, jac, 2, lamg, k2g, pr, tau, r0) 
  
      if (x(2) < 0.0_dp) then
         ! If k2 is negative, solve for lam, k instead
         lamg = lamguess(pr, tau, r0)
         k2g = sqrt(k2guess(pr, tau, r0))
         call hybrd1(fun, x, f, jac, 2, lamg, k2g, pr, tau, r0, .true.)
         ! Convert k to k2
         x(2) = x(2)**2
      end if
  
      if (any(abs(f) > 1e-12)) then
         print *, 'gaml2max did not converge!'
         print *, f
      end if
  
    contains
  
      ! Simple hybrd1 interface
      subroutine hybrd1(func, x, fvec, fjac, n, xtol, xguess, var1, var2, var3, passk1)
        external :: func, fjac
        integer, intent(in) :: n
        real(dp), intent(in) :: xtol, var1, var2, var3
        real(dp), intent(inout) :: x(n), xguess
        real(dp), intent(out) :: fvec(n)
        logical, intent(in), optional :: passk1
  
        call hybrd1_(func, n, x, fvec, fjac, xtol, maxfev, ml, mu, epsfcn, &
             diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf, &
             wa1, wa2, wa3, wa4, xguess, var1, var2, var3, passk1)
  
      end subroutine hybrd1
  
    end function gaml2max

    ! This is an attempt to implement Rich's improved gaml2max_new
      
        function gaml2max_new(pr, tau, r0) result(out)
          real(dp), intent(in) :: pr, tau, r0
          real(dp) :: out(2)
      
          real(dp) :: lam_til_min, lam_til_max
          real(dp) :: lam_til, l2, lam
      
          lam_til_min = 0.0_dp
          lam_til_max = (1.0_dp - r0*tau)/(r0 - tau)
      
          call minimize(fun, lam_til, lam_til_min, lam_til_max, pr, tau, r0)
      
          l2 = eval_l2(lam_til, pr, tau, r0)
          lam = l2*lam_til
      
          out = [lam, l2]
      
        end function gaml2max_new
      
        function fun(lam_til, pr, tau, r0) result(f)
          real(dp), intent(in) :: lam_til, pr, tau, r0
          real(dp) :: f
      
          real(dp) :: l2
      
          l2 = eval_l2(lam_til, pr, tau, r0)
          f = -l2*lam_til
      
        end function fun
      
        function eval_l2(lam_til, pr, tau, r0) result(l2)
          real(dp), intent(in) :: lam_til, pr, tau, r0
          real(dp) :: l2
      
          l2 = sqrt(pr*(1.0_dp + lam_til - r0*lam_til - r0*tau)/(r0*(1.0_dp + lam_til)*&
               (pr + lam_til)*(lam_til + tau)))
      
        end function eval_l2
      
      
        ! Minimal interface to scalar minimization
        
        subroutine minimize(func, x, xmin, xmax, var1, var2, var3)
          external :: func
          real(dp), intent(inout) :: x
          real(dp), intent(in) :: xmin, xmax, var1, var2, var3
          integer :: ierr
      
          call minimize_scalar(func, x, xmin, xmax, ierr, var1, var2, var3)
      
        end subroutine minimize
      

        subroutine minimize_scalar(func, x, xmin, xmax, ierr, var1, var2, var3)

            external func
            real(dp), intent(inout) :: x
            real(dp), intent(in) :: xmin, xmax
            integer, intent(out) :: ierr
            real(dp), intent(in), optional :: var1, var2, var3
          
            ! Parameters 
            integer, parameter :: ITMAX = 100
            real(dp), parameter :: EPS = 1.0e-5_dp
          
            ! Local variables
            integer :: iter
            real(dp) :: a,b,v,fv,w,fw
            logical :: bracketed
            external :: func
          
            a = xmin
            b = xmax  
            v = x  
            w = v
            fv = func(v, var1, var2, var3)  
            fw = fv
          
            do iter = 1, ITMAX
               call mnbrak(func, a, v, b, var1, var2, var3, bracketed)
               call brent(func, a, v, b, var1, var2, var3, w, fw)
               if (abs(x-w) < EPS*(abs(x)+abs(w))) exit
            end do
          
            x = w
            ierr = iter
          
            contains
          
              ! Local version of brent and mnbrak
              subroutine brent(func,ax,bx,cx,var1,var2,var3,xmin,fmin)
                real(dp), external :: func
                real(dp), intent(in) :: ax,bx,cx
                real(dp), intent(out) :: xmin,fmin 
                real(dp), intent(in), optional :: var1,var2,var3
                real(dp) :: a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol,t2,u,v,w,x,xm
                parameter (tol=3.0e-8_dp)
                a=min(ax,cx)
                b=max(ax,cx)
                v=bx
                w=v
                x=v
                e=0.0_dp
                fx=func(x,var1,var2,var3)
                fv=fx
                fw=fx
                do
                  xm=0.5_dp*(a+b)
                  tol1=tol*abs(x)+1.0e-5_dp
                  tol2=2.0_dp*tol1
                  if(abs(x-xm).le.(tol2-0.5_dp*(b-a))) then
                    xmin=x
                    fmin=fx
                    return
                  end if
                  if(abs(e).gt.tol1) then
                    r=(x-w)*(fx-fv)
                    q=(x-v)*(fx-fw)
                    p=(x-v)*q-(x-w)*r
                    q=2.0_dp*(q-r)
                    if(q.gt.0.0_dp) p=-p
                    q=abs(q)
                    etemp=e
                    e=d
                    if(abs(p).ge.abs(0.5_dp*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-x))then
                      e=merge(a-x,b-x, x-a.lt.x-b)
                      d=CG*e  
                    else
                      d=p/q
                      u=x+d
                      if((u-a).lt.tol2).or.(b-u).lt.tol2) then
                        d=sign(tol1,xm-x)
                      end if
                    end if
                  else
                    e=merge(a-x,b-x, x-a.lt.x-b)
                    d=CG*e
                  end if
                  if(abs(d).ge.tol1) then
                    u=x+d
                  else
                    u=x+sign(tol1,d)
                  end if
                  fu=func(u,var1,var2,var3)
                  if(fu.le.fx) then
                    if(u.ge.x) then
                      a=x
                    else
                      b=x
                    end if
                    v=w
                    fv=fw
                    w=x
                    fw=fx
                    x=u
                    fx=fu
                  else
                    if(u.lt.x) then
                      a=u
                    else
                      b=u
                    end if
                    if(fu.le.fw .or. w.eq.x) then
                      v=w 
                      fv=fw
                      w=u
                      fw=fu
                    else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
                      v=u
                      fv=fu
                    end if
                  end if
                end do
              end subroutine brent
                
              subroutine mnbrak(func,ax,bx,cx,var1,var2,var3,bracketed)
                real(dp), external :: func
                real(dp), intent(in) :: ax  
                real(dp), intent(inout) :: bx 
                real(dp), intent(out) :: cx
                logical, intent(out) :: bracketed
                real(dp), intent(in), optional :: var1, var2, var3
                real(dp) :: ulim,u,r,q,fu,dum    
                parameter (CGOLD=0.3819660_dp,GLIMIT=100.0_dp,TINY=1.0e-20_dp)
                fu=func(ax,var1,var2,var3)
                fb=func(bx,var1,var2,var3)
                if (fb.gt.fu) then
                  dum=ax
                  ax=bx
                  bx=dum
                  dum=fb
                  fb=fu
                  fu=dum
                end if
                cx=bx+CGOLD*(bx-ax)
                fc=func(cx,var1,var2,var3)
                bracketed = .true.
                if (fb.ge.fc) then
                  r=(bx-ax)*(fb-fc)
                  q=(bx-cx)*(fb-fa)
                  u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0_dp*sign(max(abs(q-r),TINY),q-r))
                  ulim=bx+GLIMIT*(cx-bx)
                  if ((bx-u)*(u-cx).gt.0.0_dp) then
                    fu=func(u,var1,var2,var3)
                    if (fu.lt.fc) then
                      ax=bx
                      bx=u
                      return
                    else if(fu.gt.fb) then
                      cx=u
                      return
                    end if
                    u=cx+GOLD*(cx-bx)
                    fu=func(u,var1,var2,var3)
                  else if((cx-u)*(u-ulim).gt.0.0_dp) then
                    fu=func(u,var1,var2,var3)
                    if(fu.lt.fc) then
                      bx=cx
                      cx=u
                      u=cx+GOLD*(cx-bx)
                      fb=fc
                      fc=fu
                      fu=func(u,var1,var2,var3)
                    end if
                  else if((u-ulim)*(ulim-cx).ge.0.0_dp) then
                    u=ulim
                    fu=func(u,var1,var2,var3)
                  else            
                    u=cx+GOLD*(cx-bx)
                    fu=func(u,var1,var2,var3)
                  end if
                  ax=bx
                  bx=cx
                  cx=u
                  fa=fb
                  fb=fc
                  fc=fu
                  return
                end if
              end subroutine mnbrak
          
          end subroutine minimize_scalar

  end module fingering_modes
