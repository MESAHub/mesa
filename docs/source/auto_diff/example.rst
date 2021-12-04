.. _autodiff example:

Example of auto_diff in run_star_extras
=======================================

.. toctree::
   :maxdepth: 1

The ``auto_diff`` module simplifies using ``run_star_extras`` hooks that
need to provide partial derivatives.
For example, the ``other_surface_PT`` hook can be used to specify an alternate
surface boundary condition by providing the surface temperature and pressure.
This hook has the following call signature::

    other_surface_PT(id, &
            skip_partials, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, ierr)

In addition to ``lnT_surf`` and ``lnP_surf``, we also need to provide four partial derivatives!
We could specify these manually.
For instance an Eddington atmosphere might have::

 subroutine other_surface_PT(id, &
            skip_partials, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, ierr)

 use const_def, only: dp
 integer, intent(in) :: id
 logical, intent(in) :: skip_partials
 real(dp), intent(out) :: &
    lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
    lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap
 integer, intent(out) :: ierr

 real(dp) :: g

 lnT_surf = log(pow(s%L(1) / (4d0 * pi * boltzsig * pow2(s%R(1))), 0.25d0))
 dlnT_dL = lnT_surf / s%L(1)
 dlnT_dlnR = -0.5d0
 dlnT_dlnM = 0d0
 dlnT_dlnkap = 0d0

 g = s%M(1) / pow2(s%R(1))
 lnP_surf = log(g / s%opacity(1))
 dlnP_dL = 0d0
 dlnP_dlnR = -2d0
 dlnP_dlnM = 1d0
 dlnP_dlnkap = -1d0

 end subroutine other_surface_PT

This relies on us calculating the derivatives correctly by hand and hard-coding them.
For this simple case that isn't too hard, but it is easy to imagine more complicated cases where doing it by hand is error-prone and tedious.
Instead, we could use ``auto_diff``::

 subroutine other_surface_PT(id, &
            skip_partials, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, ierr)

 use const_def, only: dp
 use auto_diff
 integer, intent(in) :: id
 logical, intent(in) :: skip_partials
 real(dp), intent(out) :: &
    lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
    lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap
 integer, intent(out) :: ierr

 type(auto_diff_real_4var_order1) :: L, lnR, lnM, lnkap, M, R, kap, g, lnT, lnP
 real(dp) :: g
 
 ! Set up auto_diff base variables.
 ! Note that the ordering of which variable goes in which d1val slot doesn't matter.
 ! We just need to use the same ordering when reading the results out.

 L = s%L(1)
 L%d1val1 = 1d0

 lnR = log(s%R(1))
 lnR%d1val2 = 1d0

 lnM = log(s%M(1))
 lnM%d1val3 = 1d0

 lnkap = log(s%opacity(1))
 lnkap%d1val4 = 1d0

 ! Calculate lnT and lnP

 M = exp(lnM)
 R = exp(lnR)
 kap = exp(lnkap)

 lnT = log(pow(L / (4d0 * pi * boltzsig * pow2(R))))

 g = M / pow2(R)
 lnP = log(g / kap)

 ! Read results.
 ! Note that the order matches the order we used to set up the base variables.
 lnT_surf = lnT%val
 dlnT_dL = lnT%d1val1
 dlnT_dlnR = lnT%d1val2
 dlnT_dlnM = lnT%d1val3
 dlnT_dlnkap = lnT%d1val4

 lnP_surf = lnP%val
 dlnP_dL = lnP%d1val1
 dlnP_dlnR = lnP%d1val2
 dlnP_dlnM = lnP%d1val3
 dlnP_dlnkap = lnP%d1val4


 end subroutine other_surface_PT

Note that ``auto_diff`` does not support the operation ``real variable = auto_diff variable``.
If you attempt to set a ``real`` equal to an ``auto_diff`` variable you will get an error that this operation is not implemented.
That's why we read the value of ``lnT`` and ``lnP`` above using ``lnT%val`` and ``lnP%val``.
The reason this operation is not supported is that ``real`` variables cannot store derivative information, so information would be lost in that conversion.
If this were done as an intermediate operation in a calculation (rather than during read out) the resulting derivatives that get read out at the end would be incorrect.
It is easy to make this conversion by mistake, and there is no valid use case for this operation other than at read out, so we choose to throw a compiler error rather than allow it.

Another common failure mode is not including ``use auto_diff`` at the start of the hook.
That will result in a compiler error saying that the relevant ``auto_diff`` types are not defined.

The above example with ``auto_diff`` looks more complicated than assigning derivatives by hand, but that's only because we used a simple example.
For instance the atmosphere model of Jermyn+2021_ could be implemented as::

   subroutine other_surface_PT(id, skip_partials, &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, ierr)

       use auto_diff
       integer, intent(in) :: id
       logical, intent(in) :: skip_partials
         real(dp), intent(out) :: &
            lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap
         integer, intent(out) :: ierr

      type (star_info), pointer :: s
      real(dp) :: minT, blend_minT, Ledd_rdp
      type(auto_diff_real_4var_order1) :: M, R, L, kappa
      type(auto_diff_real_4var_order1) :: R_ph, R_bondi, R_acc, vesc2, Mdot, loss, super_eddington_factor, L_incident
      type(auto_diff_real_4var_order1) :: H_disk, erf_arg, erf_f1, erf_f2, vortpar, R_AGN, R_Hill, tidepar
      type(auto_diff_real_4var_order1) :: L_stream, L_shock, Ledd, T_surf, P_surf, P_rad, P_acc, P_outflow, P_simple
      type(auto_diff_real_4var_order1) :: agn_Teff, simple_Teff, simple_lnT, simple_lnP, agn_lnT, agn_lnP
      type(auto_diff_real_4var_order1) :: lnT, lnP

      ! Basic setup
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      ! Setup autodiff variables
      M = s%m(1)
      M%d1val1 = 1d0

      R = s%r(1)
      R%d1val2 = 1d0

      L = s%L(1)
      L%d1val3 = 1d0

      kappa = const_kappa
      kappa%d1val4 = 0d0
      !kappa = s%opacity(1)
      !kappa%d1val4 = 1d0

      ! Basic quantities
      R_bondi = 2d0*standard_cgrav * M / pow2(const_csb) ! Bondi radius
      R_acc = R_bondi
      Mdot = pi * pow2(R_bondi) * const_rho * const_csb ! Mdot gain

      R_ph = pow2(kappa * const_rho) * pow3(R_bondi) ! Photosphere radius
      R_ph = R_ph + R
      R_ph = min(R_ph, R_bondi)

      ! Adjust gain down due to radiation pressure
      !Ledd_rdp = eval_Ledd(s, ierr)
      !Ledd = (Ledd_rdp / M%val) * M ! Gives Ledd the correct derivatives, because Ledd_rdp scales like M.
      Ledd = 4.0d0*pi*clight*s%cgrav(1)*M/kappa

      if (s% x_logical_ctrl(7)) then
         Ledd = Ledd * max(1e-2, (1d0 - pow2(s%job%extras_rpar(i_mu_J) / s%job%extras_rpar(i_J_crit))))
      end if

      if (s% x_logical_ctrl(1)) then
         if (L < Ledd) then
            Mdot = Mdot * (1d0 - L/Ledd)**2d0 ! Suggested by Alexander Dittmann
         else
            Mdot = 0d0
         end if
      else
         Mdot = Mdot * (1d0 - tanh(abs(L / Ledd)))
      end if

      ! Adjust gain down due to gap opening in the disk
      if (s% x_logical_ctrl(2)) then
         Mdot = Mdot * exp(-pow(M / Msun / cutoff_mass, 2d0/3d0)) ! Form suggested by Doug Lin
      end if


      if (s% x_logical_ctrl(3)) then
         H_disk = pow(2d0, 0.5d0) * const_csb / Omega_AGN
         erf_arg = -R_bondi / H_disk
         erf_f1 = (-2d0/pow(pi, 0.5d0) )*pow( 1d0 - exp(-erf_arg*erf_arg), 0.5d0)
         erf_f2 = pow(pi, 0.5d0)*0.5d0 + (3.1d1/2d2)*exp(-erf_arg*erf_arg) - (3.41d2/8d3)*exp(-2d0*erf_arg*erf_arg)
         Mdot = 0.5d0 * Mdot * pow( pi, 0.5d0 ) * erf_f1 * erf_f2 / erf_arg
      end if

      if (s% x_logical_ctrl(4)) then
            R_AGN = pow( standard_cgrav * Mass_AGN * Msun / (Omega_AGN * Omega_AGN) , 1d0/3d0 )
            vortpar = (R_bondi * R_bondi * Omega_AGN)/(2d0 * const_csb * R_AGN) !avg vorticity parameter at R_bondi
            vortpar = (2d0 / (pi * vortpar))*asinh( pow(2d0*vortpar, 1d0/3d0))  !Krumholz, McKee, Klein 2005 eqn. A7
            Mdot = Mdot * pow(1d0 + pow(vortpar, -10d0), -0.1d0) !pick minimum factor
      end if

      if (s% x_logical_ctrl(5)) then
            R_Hill = pow(standard_cgrav * s%m(1) / (3d0 * Omega_AGN*Omega_AGN), 1d0/3d0)
            tidepar = pow(R_Hill / R_bondi, 2d0)
            Mdot = Mdot * min(1d0, tidepar) !pick minimum factor
            R_acc = min(R_Hill, R_Bondi)
            ! Use min(Bondi, Hill) for accretion rate. Similar to Rosenthal et al. 2020
      end if

      if (s% x_logical_ctrl(6)) then
            H_disk = pow(2d0, 0.5d0) * const_csb / Omega_AGN
            R_Hill = pow(standard_cgrav * s%m(1) / (3d0 * Omega_AGN*Omega_AGN), 1d0/3d0)
            R_acc = min(R_Hill, R_Bondi)
            tidepar = min(1d0, R_Hill/R_bondi)
            tidepar = min(tidepar, H_disk/R_bondi)
            tidepar = tidepar * min(1d0, R_Hill/R_bondi)
            tidepar = tidepar * min(1d0, exp(2d0*(1d0-pow(R_Hill/H_disk, 2d0))))
            !Dobbs-Dixon, Li, Lin 2007, approx eqn. 28
            Mdot = Mdot * tidepar
            ! Adjust accretion rate to take into account the tidal barrier
      end if

      ! Calculate super-eddington mass loss
      vesc2 = 2d0*standard_cgrav*M / R
      super_eddington_factor = 0.1d0
      loss = super_eddington_mdot_autodiff(Ledd,M,L,R,super_eddington_factor)

      ! Set a minimum on Teff.
      minT = const_Tb

      ! Various luminosities
      L_shock = Mdot * M * standard_cgrav / R
      L_stream = L + L_shock
      L_incident = 4 * pi * pow2(R_ph) * boltz_sigma * pow4(minT)
      L_stream = L_stream + L_incident
      L_stream = max(L_stream, L_incident) ! Avoid NaN's due to bad near-surface grad(T)

      ! Evaluate simple photosphere
      simple_Teff = pow(L_stream / (4d0 * pi * R**2 * boltz_sigma), 0.25d0)
      simple_lnP = log(standard_cgrav * M / (R**2 * kappa))
      P_simple = exp(simple_lnP)
      simple_lnT = log(simple_Teff) ! First cell is tau=2/3 in this simple solution.

      ! Evaluate complicated photosphere
      agn_Teff = pow(L_stream / (4 * pi * boltz_sigma * pow2(R_ph)), 0.25d0)
      T_surf = agn_Teff * pow(R_ph / R, 5d0 / 8d0)
      agn_lnT = log(T_surf)

      ! P_surf
      P_acc = const_rho * pow2(const_csb) * pow(R_acc / R, 2.5d0)
      P_rad = crad * pow4(T_surf) / 3d0
      P_outflow = - loss * pow(vesc2, 0.5d0) / (4 * pi * pow2(R))
      P_surf = P_acc + P_rad + P_outflow + P_simple
      agn_lnP = log(P_surf)

      ! blend gets set in extras_start_step
      lnT = agn_lnT
      lnP = agn_lnP

      ! Format for output
      lnT_surf = lnT%val
      dlnT_dlnM = lnT%d1val1 * M%val
      dlnT_dlnR = lnT%d1val2 * R%val
      dlnT_dL = lnT%d1val3
      dlnT_dlnkap = lnT%d1val4 * kappa%val

      lnP_surf = lnP%val
      dlnP_dlnM = lnP%d1val1 * M%val
      dlnP_dlnR = lnP%d1val2 * R%val
      dlnP_dL = lnP%d1val3
      dlnP_dlnkap = lnP%d1val4 * kappa%val

   end subroutine other_surface_PT

Propagating the derivatives through this by hand would be tedious and error-prone.
With ``auto_diff`` the focus stays on the physics rather than differentiation.

.. _Jermyn+2021: https://arxiv.org/abs/2102.13114


