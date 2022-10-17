! ***********************************************************************
!
!   Copyright (C) 2011-2019  The MESA Team
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


module condint
   
   use const_def, only : dp
   use math_lib
   use utils_lib, only : mesa_error
   
   implicit none
   
   integer, parameter :: num_logTs = 29, num_logRhos = 71, num_logzs = 15
   !!! NB: These parameters must be consistent with the table "condtabl.d"!
   logical :: initialized = .false.
   
   real(dp) :: logTs(num_logTs), logRhos(num_logRhos), logzs(num_logzs)
   real(dp), target :: f_ary(4 * num_logRhos * num_logTs * num_logzs) ! for bicubic splines
   real(dp), pointer :: f(:, :, :, :)
   integer :: ilinx(num_logzs), iliny(num_logzs)


contains
   
   
   subroutine init_potekhin(ierr)
      use kap_def, only : kap_dir
      use interp_2d_lib_db, only : interp_mkbicub_db
      integer, intent(out) :: ierr
      
      character (len = 256) :: filename
      integer :: read_err, iz, it, ir, shift
      integer :: ibcxmin                   ! bc flag for x=xmin
      real(dp) :: bcxmin(num_logTs)               ! bc data vs. y at x=xmin
      integer :: ibcxmax                   ! bc flag for x=xmax
      real(dp) :: bcxmax(num_logTs)               ! bc data vs. y at x=xmax
      integer :: ibcymin                   ! bc flag for y=ymin
      real(dp) :: bcymin(num_logRhos)               ! bc data vs. x at y=ymin
      integer :: ibcymax                   ! bc flag for y=ymax
      real(dp) :: bcymax(num_logRhos)               ! bc data vs. x at y=ymax
      real(dp) :: Z
      real(dp), pointer :: f1(:)
      
      include 'formats'
      
      ierr = 0
      if (initialized) return
      
      shift = 4 * num_logRhos * num_logTs
      f(1:4, 1:num_logRhos, 1:num_logTs, 1:num_logzs) => f_ary(1:shift * num_logzs)
      
      filename = trim(kap_dir) // '/condtabl.data'
      open(1, file = trim(filename), status = 'OLD', iostat = ierr)
      if (ierr /= 0) then
         write(*, '(A)')
         write(*, '(A)')
         write(*, '(A)')
         write(*, '(A)')
         write(*, '(A)')
         write(*, *) 'NOTICE: missing ' // trim(filename)
         write(*, *) 'Please remove the directory mesa/data/kap_data,'
         write(*, *) 'and rerun the mesa ./install script.'
         write(*, '(A)')
         call mesa_error(__FILE__, __LINE__)
      end if
      !print*,'Reading thermal conductivity data...'
      read_err = 0
      read(1, '(A)', iostat = read_err) ! skip the first line
      if (read_err /= 0) ierr = read_err
      do iz = 1, num_logzs
         read (1, *, iostat = read_err) z, (logTs(it), it = 1, num_logTs)
         if (read_err /= 0) ierr = read_err
         if (z .eq. 1d0) then
            logzs(iz) = 0d0
         else
            logzs(iz) = log10(z)
         endif
         do ir = 1, num_logRhos
            read(1, *, iostat = read_err) logRhos(ir), (f(1, ir, it, iz), it = 1, num_logTs)
            if (read_err /= 0) ierr = read_err
         end do
      end do
      close(1)
      if (ierr /= 0) then
         write(*, '(A)')
         write(*, '(A)')
         write(*, '(A)')
         write(*, '(A)')
         write(*, '(A)')
         write(*, *) 'NOTICE: error trying to read ' // trim(filename)
         write(*, *) 'Please remove the directory mesa/data/kap_data,'
         write(*, *) 'and rerun the mesa ./install script.'
         write(*, '(A)')
         call mesa_error(__FILE__, __LINE__)
      end if
      ! boundary condition is slope=0 at edges of tables
      ! this ensures continuous derivatives when we clip
      ibcxmin = 3; bcxmin(1:num_logTs) = 0d0
      ibcxmax = 3; bcxmax(1:num_logTs) = 0d0
      ibcymin = 3; bcymin(1:num_logRhos) = 0d0
      ibcymax = 3; bcymax(1:num_logRhos) = 0d0
      do iz = 1, num_logzs
         f1(1:shift) => f_ary(1 + (iz - 1) * shift:iz * shift)
         call interp_mkbicub_db(&
            logRhos, num_logRhos, logTs, num_logTs, f1, num_logRhos, &
            ibcxmin, bcxmin, ibcxmax, bcxmax, &
            ibcymin, bcymin, ibcymax, bcymax, &
            ilinx(iz), iliny(iz), ierr)
         if (ierr /= 0) then
            write(*, '(A)')
            write(*, '(A)')
            write(*, '(A)')
            write(*, '(A)')
            write(*, '(A)')
            write(*, *) 'NOTICE: error in ' // trim(filename)
            write(*, *) 'Please report the problem.'
            write(*, '(A)')
            call mesa_error(__FILE__, __LINE__)
         end if
      end do
      initialized = .true.
   end subroutine init_potekhin
   
   
   subroutine do_electron_conduction_potekhin(&
      zbar, logRho_in, logT_in, kap, dlogkap_dlogRho, dlogkap_dlogT, ierr)
      
      use const_def, only : boltz_sigma
      real(dp), intent(in) :: zbar, logRho_in, logT_in
      real(dp), intent(out) :: kap, dlogkap_dlogRho, dlogkap_dlogT
      integer, intent(out) :: ierr
      
      integer :: iz, iz1, iz2, shift
      real(dp) :: zlog, logRho, logT
      real(dp) :: alfa, beta, &
         logK1, kap1, dlogK1_dlogRho, dlogK1_dlogT, &
         logK2, kap2, dlogK2_dlogRho, dlogK2_dlogT, &
         logK, logkap, dlogK_dlogRho, dlogK_dlogT
      real(dp), pointer :: f1(:)
      
      logical :: clipped_logRho, clipped_logT
      
      include 'formats'
      
      ierr = 0
      shift = 4 * num_logRhos * num_logTs
      
      if (logRho_in .lt. logRhos(1)) then
         logRho = logRhos(1)
         clipped_logRho = .true.
      else if (logRho_in .gt. logRhos(num_logRhos)) then
         logRho = logRhos(num_logRhos)
         clipped_logRho = .true.
      else
         logRho = logRho_in
         clipped_logRho = .false.
      end if
      
      if (logT_in .lt. logTs(1)) then
         logT = logTs(1)
         clipped_logT = .true.
      else if (logT_in .gt. logTs(num_logTs)) then
         logT = logTs(num_logTs)
         clipped_logT = .true.
      else
         logT = logT_in
         clipped_logT = .false.
      end if
      
      zlog = max(logzs(1), min(logzs(num_logzs), log10(max(1d-30, zbar))))
      
      if (zlog <= logzs(1)) then ! use 1st
         call get1(1, logK, dlogK_dlogRho, dlogK_dlogT, ierr)
      else if (zlog >= logzs(num_logzs)) then ! use last
         call get1(num_logzs, logK, dlogK_dlogRho, dlogK_dlogT, ierr)
      else ! interpolate
         iz1 = -1
         do iz = 2, num_logzs
            if (zlog >= logzs(iz - 1) .and. zlog <= logzs(iz)) then
               iz1 = iz - 1; iz2 = iz; exit
            end if
         end do
         if (iz1 < 0) then
            write(*, 2) 'num_logzs', num_logzs
            do iz = 1, num_logzs
               write(*, 2) 'logzs(iz)', iz, logzs(iz)
            end do
            write(*, 1) 'zlog', zlog
            write(*, *) 'confusion in do_electron_conduction'
            call mesa_error(__FILE__, __LINE__)
         end if
         
         call get1(iz1, logK1, dlogK1_dlogRho, dlogK1_dlogT, ierr)
         if (ierr /= 0) then
            write(*, *) 'interp failed for iz1 in do_electron_conduction', iz1, logRho, logT
            call mesa_error(__FILE__, __LINE__)
         end if
         
         call get1(iz2, logK2, dlogK2_dlogRho, dlogK2_dlogT, ierr)
         if (ierr /= 0) then
            write(*, *) 'interp failed for iz2 in do_electron_conduction', iz2, logRho, logT
            call mesa_error(__FILE__, __LINE__)
         end if
         
         ! linear interpolation in zlog
         alfa = (zlog - logzs(iz1)) / (logzs(iz2) - logzs(iz1))
         beta = 1d0 - alfa
         logK = alfa * logK2 + beta * logK1
         if (clipped_logRho) then
            dlogK_dlogRho = 0
         else
            dlogK_dlogRho = alfa * dlogK2_dlogRho + beta * dlogK1_dlogRho
         end if
         if (clipped_logT) then
            dlogK_dlogT = 0
         else
            dlogK_dlogT = alfa * dlogK2_dlogT + beta * dlogK1_dlogT
         end if
      end if
      
      ! chi = thermal conductivity, = 10**logK (cgs units)
      ! conduction opacity kappa = 16*boltz_sigma*T^3 / (3*rho*chi)
      ! logkap = 3*logT - logRho - logK + log10(16*boltz_sigma/3)
      
      logkap = 3d0 * logT_in - logRho_in - logK + log10(16d0 * boltz_sigma / 3d0)
      
      kap = exp10(logkap)
      dlogkap_dlogRho = -1d0 - dlogK_dlogRho
      dlogkap_dlogT = 3d0 - dlogK_dlogT
   
   contains
      
      
      subroutine get1(iz, logK, dlogK_dlogRho, dlogK_dlogT, ierr)
         use kap_eval_support, only : Do_Kap_Interpolations
         integer, intent(in) :: iz
         real(dp), intent(out) :: logK, dlogK_dlogRho, dlogK_dlogT
         integer, intent(out) :: ierr
         logical, parameter :: dbg = .false.
         real(dp) :: fval, df_dx, df_dy, &
            logRho0, logRho1, logT0, logT1
         integer :: i_logRho, j_logT, k
         include 'formats'
         ierr = 0
         f1(1:shift) => f_ary(1 + (iz - 1) * shift:iz * shift)
         
         if (logRho < logRhos(2)) then
            i_logRho = 1
         else
            i_logRho = num_logRhos - 1
            do k = 2, num_logRhos - 1
               if (logRho >= logRhos(k) .and. logRho < logRhos(k + 1)) then
                  i_logRho = k; exit
               end if
            end do
         end if
         logRho0 = logRhos(i_logRho)
         logRho1 = logRhos(i_logRho + 1)
         
         if (logT < logTs(2)) then
            j_logT = 1
         else
            j_logT = num_logTs - 1
            do k = 2, num_logTs - 1
               if (logT >= logTs(k) .and. logT < logTs(k + 1)) then
                  j_logT = k; exit
               end if
            end do
         end if
         logT0 = logTs(j_logT)
         logT1 = logTs(j_logT + 1)
         
         call Do_Kap_Interpolations(&
            f1, num_logRhos, num_logTs, i_logRho, j_logT, logRho0, &
            logRho, logRho1, logT0, logT, logT1, logK, dlogK_dlogRho, dlogK_dlogT)
         if (ierr /= 0) return
      
      end subroutine get1
   
   
   end subroutine do_electron_conduction_potekhin
   
   
   subroutine do_electron_conduction_blouin(&
      zbar, logRho, logT, &
      kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
      use const_def, only : dp
      use auto_diff
      real(dp), intent(in) :: zbar ! average ionic charge (for electron conduction)
      real(dp), intent(in) :: logRho ! the density
      real(dp), intent(in) :: logT ! the temperature
      real(dp), intent(out) :: kap ! electron conduction opacity
      real(dp), intent(out) :: dlnkap_dlnRho ! partial derivative at constant T
      real(dp), intent(out) :: dlnkap_dlnT   ! partial derivative at constant Rho
      integer, intent(out) :: ierr ! 0 means AOK.
      
      ! this implements the correction formulae from Blouin et al. (2020)
      ! https://ui.adsabs.harvard.edu/abs/2020ApJ...899...46B/abstract
      
      real(dp), parameter :: alpha_H = -0.52d0, alpha_He = -0.46d0
      real(dp), parameter :: a_H = 2.0d0, a_He = 1.25d0
      real(dp), parameter :: b_H = 10.0d0, b_He = 2.5d0
      real(dp), parameter :: logRho0_H = 5.45d0, logRho0_He = 6.50d0
      real(dp), parameter :: logT0_H = 8.40d0, logT0_He = 8.57d0
      real(dp), parameter :: sigRho_H = 5.14d0, sigRho_He = 6.20d0
      real(dp), parameter :: sigT_H = 0.45d0, sigT_He = 0.55d0
      
      type(auto_diff_real_2var_order1) :: logRho_auto, logT_auto
      type(auto_diff_real_2var_order1) :: Rhostar, Tstar
      type(auto_diff_real_2var_order1) :: g_H, g_He, H_H, H_He
      type(auto_diff_real_2var_order1) :: log_correction, log_correction_H, log_correction_He
      
      real(dp) :: alfa, frac_H, frac_He
      
      ! auto_diff
      ! var1: lnRho
      ! var2: lnT
      
      logRho_auto = logRho
      logRho_auto% d1val1 = iln10
      logT_auto = logT
      logT_auto% d1val2 = iln10
      
      ! call previous standard routines (Cassisi/Potekhin)
      call do_electron_conduction_potekhin(&
         zbar, logRho, logT, &
         kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
      
      ! combined correction
      !
      ! The thermal conductivity is tabulated at Zbar = {1,2,3,4,6,...}
      ! and linear interpolation in logK vs logZbar is applied.
      !
      ! Therefore, we apply the Blouin+ 2020 corrections in a manner
      ! equivalent to individually correcting the Zbar = {1,2} tables.
      
      if (Zbar .le. 1d0) then ! all H
         frac_H = 1d0
         frac_He = 0d0
      else if (Zbar .le. 2d0) then ! mix H and He
         alfa = (log10(Zbar) - log10(1d0)) / (log10(2d0) - log10(1d0))
         frac_H = 1d0 - alfa
         frac_He = alfa
      else if (Zbar .le. 3d0) then ! mix He and no correction
         alfa = (log10(Zbar) - log10(2d0)) / (log10(3d0) - log10(2d0))
         frac_H = 0d0
         frac_He = 1d0 - alfa
      else ! no correction
         frac_H = 0d0
         frac_He = 0d0
         return
      end if
      
      if (frac_H .gt. 0) then
         
         ! correction for H
         Rhostar = logRho_auto - logRho0_H
         Tstar = logT_auto - logT0_H
         
         g_H = a_H * exp(&
            -pow2(Tstar * cos(alpha_H) + Rhostar * sin(alpha_H)) / pow2(sigT_H)  &
               - pow2(Tstar * sin(alpha_H) - Rhostar * cos(alpha_H)) / pow2(sigRho_H))
         
         H_H = 0.5d0 * tanh(b_H * (g_H - 0.5d0)) + 0.5d0
         
         log_correction_H = -log(1d0 + g_H * H_H)
      
      else
         
         log_correction_H = 0d0
      
      end if
      
      if (frac_He .gt. 0) then
         
         ! correction for He
         Rhostar = logRho_auto - logRho0_He
         Tstar = logT_auto - logT0_He
         
         g_He = a_He * exp(&
            -pow2(Tstar * cos(alpha_He) + Rhostar * sin(alpha_He)) / pow2(sigT_He)  &
               - pow2(Tstar * sin(alpha_He) - Rhostar * cos(alpha_He)) / pow2(sigRho_He))
         
         H_He = 0.5d0 * tanh(b_He * (g_He - 0.5d0)) + 0.5d0
         
         log_correction_He = -log(1d0 + g_He * H_He)
      
      else
         
         log_correction_He = 0d0
      
      end if
      
      
      ! blend H correction, He correction, and other correction (none)
      log_correction = frac_H * log_correction_H + frac_He * log_correction_He
      
      ! apply correction factor
      kap = exp(log(kap) + log_correction% val)
      dlnkap_dlnRho = dlnkap_dlnRho + log_correction% d1val1
      dlnkap_dlnT = dlnkap_dlnT + log_correction% d1val2
   
   end subroutine do_electron_conduction_blouin
   
   subroutine do_electron_conduction(&
      rq, zbar, logRho, logT, &
      kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
      use kap_def, only : Kap_General_Info
      type (Kap_General_Info), pointer, intent(in) :: rq
      real(dp), intent(in) :: zbar ! average ionic charge (for electron conduction)
      real(dp), intent(in) :: logRho ! the density
      real(dp), intent(in) :: logT ! the temperature
      real(dp), intent(out) :: kap ! electron conduction opacity
      real(dp), intent(out) :: dlnkap_dlnRho ! partial derivative at constant T
      real(dp), intent(out) :: dlnkap_dlnT   ! partial derivative at constant Rho
      integer, intent(out) :: ierr ! 0 means AOK.
      
      if (rq% use_blouin_conductive_opacities) then
         call do_electron_conduction_blouin(&
            zbar, logRho, logT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
      else
         call do_electron_conduction_potekhin(&
            zbar, logRho, logT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
      end if
   
   end subroutine do_electron_conduction


end module condint
