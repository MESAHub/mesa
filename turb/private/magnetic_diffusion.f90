! ***********************************************************************
!   magnetic_diffusion module
!   Part of MESA software suite
!
!   This module calculates various properties related to magnetic diffusion.
!
!   Contains:
!   - calc_eta: Computes magnetic diffusivity from electrical conductivity.
!   - calc_sige: Computes electrical conductivity for different regimes.
!   - sige1, sige2, sige3: Helper functions for electrical conductivity calculations.
!
!   References:
!   - S.-C. Yoon Oct. 10, 2003
!   - Spitzer 1962
!   - Wendell et al. 1987, ApJ 313:284
!   - Nandkumar & Pethick (1984)
!
!   This code is distributed under the terms of the MESA MANIFESTO and the GNU General Library Public License.
! ***********************************************************************

module magnetic_diffusion

    use const_def
    use num_lib
    use math_lib
    use utils_lib
    implicit none

    private
    public :: calc_eta, calc_sige, sige1, sige2, sige3

contains


!### Example of how to use this: 

! call kap_get_elect_cond_opacity( &
! s% kap_handle, s% zbar(k), log10(s% rho(k)),& 
! log10(s% T(k)), kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
! if (ierr /= 0) return
! kap_cond(k) = kap

! !## Then calculate the conductivity 

! sig = calc_sige(s% abar(k), s% zbar(k),s% rho(k), & 
! s% T(k),s% Cp(k), kap_cond(k),s% opacity(k))
! cond(k) = sig

! !## and the magnetic diffusivity

! eta(k) = clight * clight / (4d0 * pi * sig)




    !> Compute magnetic diffusivity from electrical conductivity.
    real(dp) function calc_eta(sig)
        real(dp), intent(in) :: sig
        calc_eta = (clight * clight / (4d0 * pi)) / sig
    end function calc_eta

    !> Compute electrical conductivity for various regimes.
    real(dp) function calc_sige(abar, zbar, rho, T, Cp, kap_cond, opacity)
        real(dp), intent(in) :: abar, zbar, rho, T, Cp, kap_cond, opacity
        real(dp) :: gamma, xlg, alpha, ffff, xxx, xsig1, xsig2, xsig3

        gamma = 0.2275d0 * pow2(zbar) * pow(rho * 1.d-6 / abar, one_third) * 1.d8 / T
        xlg = log10(gamma)
        alpha = 16d0 * boltz_sigma * pow3(T) / (3d0 * opacity * pow2(rho) * Cp)

        select case (true)
            case (xlg < -1.5d0)
                calc_sige = sige1(zbar, T, gamma)
            case (xlg >= -1.5d0 .and. xlg <= 0d0)
                xxx = (xlg + 0.75d0) * 4d0 / 3d0
                ffff = 0.25d0 * (2d0 - 3d0 * xxx + pow3(xxx))
                xsig1 = sige1(zbar, T, gamma)
                xsig2 = sige2(T, rho, kap_cond)
                calc_sige = (1d0 - ffff) * xsig2 + ffff * xsig1
            case (xlg > 0d0 .and. xlg < 0.5d0)
                xsig2 = sige2(T, rho, kap_cond)
                calc_sige = xsig2
            case (xlg >= 0.5d0 .and. xlg < 1d0)
                xxx = (xlg - 0.75d0) * 4d0
                ffff = 0.25d0 * (2d0 - 3d0 * xxx + pow3(xxx))
                xsig2 = sige2(T, rho, kap_cond)
                xsig3 = sige3(zbar, T, gamma)
                calc_sige = (1d0 - ffff) * xsig3 + ffff * xsig2
            case (xlg >= 1d0)
                calc_sige = sige3(zbar, T, gamma)
        end select
    end function calc_sige

    !> Computes one regime of the electrical conductivity.
    real(dp) function sige1(z, t, xgamma)
        real(dp), intent(in) :: z, t, xgamma
        real(dp) :: etan, xlambda, f
        if (t >= 4.2d5) then
            f = sqrt(4.2d5 / t)
        else
            f = 1.d0
        end if
        xlambda = sqrt(3d0 * z * z * z) * pow(xgamma, -1.5d0) * f + 1d0
        etan = 3.d11 * z * log(xlambda) * pow(t, -1.5d0)
        etan = etan / (1.d0 - 1.20487d0 * exp(-1.0576d0 * pow(z, 0.347044d0)))
        sige1 = clight * clight / (pi4 * etan)
    end function sige1

    !> Computes electrical conductivity using conductive opacity.
    real(dp) function sige2(t, rho, kap_cond)
        real(dp), intent(in) :: t, rho, kap_cond
        sige2 = 1.11d9 * t * t / (rho * kap_cond)
    end function sige2

    !> Computes electrical conductivity in degenerate matter.
    real(dp) function sige3(z, t, xgamma)
        real(dp), intent(in) :: z, t, xgamma
        real(dp) :: rme, rm23, ctmp, xi
        rme = 8.5646d-23 * t * t * t * xgamma * xgamma * xgamma / pow5(z)
        rm23 = pow(rme, 2d0 / 3d0)
        ctmp = 1d0 + 1.018d0 * rm23
        xi = sqrt(3.14159d0 / 3.) * log(z) / 3.d0 + 2.d0 * log(1.32d0 + 2.33d0 / sqrt(xgamma)) / 3.d0 - 0.484d0 * rm23 / ctmp
        sige3 = 8.630d21 * rme / (z * ctmp * xi)
    end function sige3

end module magnetic_diffusion