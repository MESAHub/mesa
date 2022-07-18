      module op_eval_mombarg

      use const_def, only: dp
      !use crlibm_lib

      implicit none
      contains


      !!! Compute g_rad for a single cell. (Currently not used.)
      subroutine compute_grad(k, fk, logT_face, logRho_face,l, r,lkap_ross_cell, lgrad_cell, ierr,ite,jne,epatom,amamu,sig,eumesh)
        ! OP mono data for: H, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Cr, Mn, Fe, and Ni.
      use chem_def, only: chem_isos, ih1, ihe3, ihe4, ic12, in14, io16, ine20, ina23,img24, ial27, isi28, is32, iar40, ica40, icr52, imn55, ife56, ini58
      use interp_2d_lib_sg
      use cubic_interpolator, only: interpolator

      integer, intent(inout) :: ierr
      integer, intent(in) :: k
      real(dp), intent(in) :: fk(:), logT_face, logRho_face, l, r
      integer :: nel, nptot
      parameter(nel = 17, nptot = 10000) !number of elements and number of u-mesh points.
      real(dp), intent(out) :: lkap_ross_cell, lgrad_cell(nel)

      integer , pointer :: ite(:),jne(:)
      real(dp), pointer :: sig(:,:,:)
      real(dp), pointer  :: epatom(:,:),amamu(:),eumesh(:,:,:)


      ! ite: temperature index, logT = 0.025*ite, (1648)
      ! jne: electron density index, logNe = 0.25*jne, (1648)
      ! epatom: electrons per atom, (nel,nptot).
      ! amamu: mu = mean atomic weight (nel)
      ! sig: monochromatic cross sections in (a0^2), (nel,1648,nptot)
      ! eumesh: sigma_i*(1 - exp(-u(v))) - a_k,i, where u=h*nu/(kT). OP mono grid is equally spaced in variable v. a_k,i are the correction factors. (nel,1648,nptot)


      integer :: n, ke, nz, id, m, ik, i

      real(dp):: epa_mix_cell(1648), amu_mix_cell, logRho(1648),logT(1648) ! Number of electrons per atom, mean molecular weight, density and temperature as a function of ite (temp index) and jne (density index) from the OP mono data.
      real(dp) :: delta(1648)
      real(dp) :: log10_bohr_radius_sqr = -16.55280
      real(dp) :: lgamm_cell(nel) !interpolated log kappa Rossland and log gamma_k in each cell (k = element index).
      real(dp) :: logNa = 23.779750912481397 !log10_cr(6.0221409d23) Avogadro's number
      real(dp) :: dv = (1.0552976117319748 - 0.00010565516589892675)/nptot !v(u(1)) - v(u(nptot))/nptot
      real(dp) :: mH

      integer ::  delta_min_idx
      real(dp) :: lkap_Ross(4,4),gamma_k(4,4,nel),lgamm(4,4,nel),sig_Ross(4,4) ! Rossland cross section, log kappa, gamma_k, log gamma_k for interpolation.
      real(dp) :: sf, flux !'scale factor' for g_rad, local flux.
      real(dp), parameter :: pi = 3.141592653589793239
      real(dp), parameter :: log_c = 10.476820702927927 !log10_cr(dble(299792458e2)) c = speed of light
      real, allocatable :: sig_mix_cell(:,:,:),sig_int(:)
      integer :: ii, jj, ite_min, jne_min, ii_min, jj_min, ite_step, jne_step
      integer :: ite_i, jne_i, dite, djne, i_grid(4,4)
      real(dp) :: logT_min, logRho_min, logT_grid(4,4), logRho_grid(4,4)
      integer :: offset1, offset2, tries, missing_point(4,4)
      real(dp) :: log_amu_mix_cell, lkap_ross_face, gam
      integer :: imin, imax
      logical :: retry

      type(interpolator) :: rossl_interpolator
      type(interpolator) :: gaml_interpolators(nel)

      ierr = 0


      imin = 1 !-1
      imax = 1648 !-1
      ite_step = 2
      jne_step = 2

      !!! Compute the number of electrons per atom for the local mixture.
      epa_mix_cell = 0d0
        do i=imin,imax
          epa_mix_cell(i) = dot_product(fk,epatom(1:nel,i))
        enddo

      amu_mix_cell = dot_product(fk,amamu)

      mH = chem_isos% W(ih1) * 1.660538782d-24
      log_amu_mix_cell = log10(amu_mix_cell * mH)
      !logRho = 0.25*jne + log_amu_mix_cell - log10(epa_mix_cell)
      logRho = 0.25*jne + log10(amu_mix_cell) - log10(epa_mix_cell) -logNa
      logT   = 0.025*ite
      !write(*,*) 'rho computed'

      !!! Select nearest points in logT,logRho for interpolation.
      !!! First, find the nearest OP data point, and check which of the four possible positionings this minimum has wrt to (T,Rho)_cell.
      !!! Acquire the remaining 15 points of the 4x4 grid, where (T,Rho)_cell is located in the inner square.

      delta = sqrt((logRho - logRho_face)*(logRho - logRho_face)/0.25 + (logT-logT_face)*(logT-logT_face)/0.025)

      delta_min_idx = MINLOC(delta, DIM=1)
      ite_min   = ite(delta_min_idx)
      jne_min   = jne(delta_min_idx)
      logT_min   = logT(delta_min_idx)
      logRho_min = logRho(delta_min_idx)
      if (logT_min <= logT_face .and. logRho_min <= logRho_face) then
        ii_min = 2
        jj_min = 2
      else if (logT_min < logT_face .and. logRho_min > logRho_face) then
        ii_min = 2!3
        jj_min = 3!2
      else if (logT_min > logT_face .and. logRho_min < logRho_face) then
        ii_min = 3
        jj_min = 2
      else if (logT_min > logT_face .and. logRho_min > logRho_face) then
        ii_min = 3
        jj_min = 3
      endif
      !write(*,*) 'ii, jj min', ii_min, jj_min, logT_min, XC, logRho_min, logRho_cell

      offset1 = 0
      offset2 = 0
      missing_point = 1
      tries = 0
      retry = .true.
      !do while (point_not_found .ne. 0) !If (T,Rho)_cell lies to close to the edge of the OP mono grid (in Rho), shift the 4x4 grid up/down by 1 in jne.
      do while (retry) !If (T,Rho)_cell lies to close to the edge of the OP mono grid (in Rho), shift the 4x4 grid up/down by 1 in jne.
      missing_point = 1
      retry = .false.
      do jj=1,4
            do ii=1,4
              dite = (ii - ii_min + offset1)*ite_step !(ii - ii_min)*ite_step + offset2
              djne = (jj - jj_min + offset2)*jne_step !(jj - jj_min)*jne_step + offset1
             do i =imin,imax
               ite_i = ite(i)
               jne_i = jne(i)

                if (ite_i .eq. ite_min +  dite .and. jne_i .eq. jne_min + djne) THEN
                  logT_grid(ii,jj) = logT(i)
                  logRho_grid(ii,jj) = logRho(i)
                  i_grid(ii,jj) = i
                  !write(*,*) ite_i, jne_i, ii, jj, i_grid(ii,jj), ',',logT_grid(ii,jj),&
                  !',', logRho_grid(ii,jj)
                  missing_point(ii,jj) = 0
                endif
              enddo
            enddo
      enddo
      if (SUM(missing_point) > 0) then
        retry = .true.
        if (SUM(missing_point(2:4,1:3)) == 0) then
          offset1 = offset1 + 1
          offset2 = offset2 - 1
        else if (SUM(missing_point(1:3,1:3)) == 0) then
          offset1 = offset1 - 1
          offset2 = offset2 - 1
        else if (SUM(missing_point(2:4,2:4)) == 0) then
          offset1 = offset1 + 1
          offset2 = offset2 + 1
        else if (SUM(missing_point(1:3,2:4)) == 0) then
          offset1 = offset1 - 1
          offset2 = offset2 + 1
        else
          !write(*,*) 'warning interpolation quality in cell grad', k,ite_min,jne_min, &
          !logT_face, logRho_face, ii_min, jj_min
          if (ii_min == 3 .and. jj_min == 2) then
            offset1 = offset1 + 2
            offset2 = offset2 - 2
          else if (ii_min == 2 .and. jj_min == 3) then
            offset1 = offset1 - 2
            offset2 = offset2 + 2
          else if (ii_min == 2 .and. jj_min == 2) then
            offset1 = offset1 - 2
            offset2 = offset2 - 2
          else if (ii_min == 3 .and. jj_min == 3) then
            offset1 = offset1 + 2
            offset2 = offset2 - 1
          endif
        endif
      endif


      if (tries > 2) THEN ! To prevent loop from getting stuck.
        write(*,*) 'Cannot find points for interpolation compute_grad', k, ite_min, jne_min, logT_min, logRho_min, logT_face, logRho_face, missing_point, ii_min, jj_min, offset1, offset2, imin, imax, log_amu_mix_cell
        ierr = 1
        return
      endif
      tries = tries + 1
      enddo
      !write(*,*) 'Points for interpolation selected', k

      !!! Compute the monochromatic cross-section for the local mixture.
      allocate(sig_mix_cell(4,4,nptot),sig_int(nptot-1), stat=ierr)
      sig_mix_cell = 0d0
      do jj=1,4
          do ii=1,4
            ik = i_grid(ii,jj)!+ 1648*(ke-1)
            do m=1,nptot
              sig_mix_cell(ii,jj,m) = dot_product(fk,sig(:,ik,m))
            enddo
          enddo
      enddo

      !!! Compute the Rossland mean cross-section by integrating over variable v (mesh equally spaced in v).
      sig_Ross = 0d0
      do jj=1,4
        do ii=1,4
           sig_int(1:nptot-1) = (1/sig_mix_cell(ii,jj,1:nptot-1) + 1/sig_mix_cell(ii,jj,2:nptot))/2. * dv
           do m=1, nptot-1
                !sig_Ross(ii,jj) = sig_Ross(ii,jj) + (1/sig_mix_cell(ii,jj,m) + 1/sig_mix_cell(ii,jj,m+1))/2. * dv
                sig_Ross(ii,jj) = sig_Ross(ii,jj) + sig_int(m)
           enddo
        enddo
      enddo

      lkap_Ross =  log10_bohr_radius_sqr - log_amu_mix_cell - log10(sig_Ross)

      call rossl_interpolator% initialize()
      do jj = 1, 4
         do ii = 1, 4
            call rossl_interpolator% add_point(logT_grid(ii,jj), logRho_grid(ii,jj), lkap_Ross(ii,jj))
         enddo
      enddo

      lkap_ross_cell  = rossl_interpolator% evaluate(logT_face,logRho_face)


      !!! Compute gamma for each element for interpolation.
      gamma_k = 0d0

      do jj=1,4
        do ii=1,4
          ik = i_grid(ii,jj)! + 1648*(ke-1)
          do ke=3,nel
            gam = 0
            do m=1, nptot-1
            gam = gam + ((eumesh(ke,ik,m)/ sig_mix_cell(ii,jj,m)) +(eumesh(ke,ik,m+1)/ sig_mix_cell(ii,jj,m+1)))
            enddo
           gamma_k(ii,jj,ke) = gam/2. * dv
          enddo
        enddo
      enddo

      deallocate(sig_mix_cell,sig_int)

      where (gamma_k < 0.)
        gamma_k = 10d-30
      endwhere
      lgamm = log10(gamma_k)

      do ke = 3, nel
         call gaml_interpolators(ke)% initialize()
         do jj = 1, 4
            do ii = 1, 4
               call gaml_interpolators(ke)% add_point(logT_grid(ii,jj), logRho_grid(ii,jj), lgamm(ii,jj,ke))
            enddo
         enddo

        lgamm_cell(ke) = gaml_interpolators(ke)% evaluate(logT_face, logRho_face) !cell
      enddo
      lgamm_cell(1) = -30.
      lgamm_cell(2) = -30.
      !write(*,*) 'lgamm_cell computed'!, ke, lgamm_cell(ke)

      ! Compute log g_rad.
        flux  = l / (4d0*pi*r*r)
        sf = lkap_ross_cell  - log_c + log10(flux) + log10(amu_mix_cell)
        lgrad_cell = sf + lgamm_cell - log10(amamu)
      lgrad_cell(1) = -30.
      lgrad_cell(2) = -30.
      !write(*,*) 'log kappa + log g_rad Fe cell', k, lkap_ross_cell, lgrad_cell(16)
      !write(*,*) 'compute_grad done', k

      end subroutine compute_grad


      !!! Compute gamma factors and kappa_Ross for all OP mono data points for a given mixture.
      subroutine compute_gamma_grid(ngp, fk_all,lgamm_pcg, lkap_face_pcg, logT_pcg, logRho_pcg, ierr,ite,jne,epatom,amamu,sig,eumesh)
        ! OP mono data for: H, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Cr, Mn, Fe, and Ni.
      use chem_def, only: chem_isos, ih1, ihe3, ihe4, ic12, in14, io16, ine20, ina23,img24, ial27, isi28, is32, iar40, ica40, icr52, imn55, ife56, ini58
      use interp_2d_lib_sg
      use cubic_interpolator, only: interpolator

      integer, intent(inout) :: ierr
      integer, intent(in) :: ngp
      real(dp), intent(in) :: fk_all(:)
      integer :: nel, nptot
      parameter(nel = 17, nptot = 10000) !number of elements and number of u-mesh points.
      !real(dp), pointer, intent(out)::lgamm_pcg(:,:,:),lkap_face_pcg(:,:)
      real(dp), intent(out)::lgamm_pcg(nel,1648),lkap_face_pcg(1648)

      real(dp), intent(out) :: logT_pcg(1648), logRho_pcg(1648)

      integer , pointer :: ite(:),jne(:)
      real(dp), pointer :: sig(:,:,:)!, ak_f(:,:)
      real(dp), pointer :: epatom(:,:),amamu(:),eumesh(:,:,:)
      integer :: imin, imax

      integer :: n, ke, nz, id, m, ik, i, j

      !real(dp) :: fk_norm_fac !Local fractional abudance per element and normalization factor.
      real(dp):: epa_mix_cell(1648), amu_mix_cell, fk(nel) ! Number of electrons per atom, mean molecular weight, density and temperature as a function of ite (temp index) and jne (density index) from the OP mono data.
      !integer ::  eid(nel)
      real(dp) :: log10_bohr_radius_sqr = -16.55280
      real(dp) :: logNa = 23.779750912481397 !log10_cr(6.0221409d23) Avogadro's number
      real(dp) :: dv = (1.0552976117319748 - 0.00010565516589892675)/nptot !v(u(1)) - v(u(nptot))/nptot
      real(dp) :: mH

      !!!! For interpolator.


      real(dp) :: sig_Ross(1648), gamma_k(nel,1648)!, lgamm(nel,1648) ! Rossland cross section, log kappa, gamma_k, log gamma_k for interpolation.


      real, allocatable :: sig_mix_cell(:,:),sig_int(:)
      real(dp) :: log_amu_mix_cell, gam


      ierr = 0

      !!! Compute an estimated temperature range.
      imin = 1 !-1
      imax = 1648 !-1

      fk = fk_all!(j,:)
      !!! Compute the number of electrons per atom for the local mixture.
      epa_mix_cell = 0d0
        do i=imin,imax
          epa_mix_cell(i) = dot_product(fk,epatom(1:nel,i))
        enddo

      amu_mix_cell = dot_product(fk,amamu)

      mH = chem_isos% W(ih1) * 1.660538782d-24
      log_amu_mix_cell = log10(amu_mix_cell * mH)


      !!! Compute the monochromatic cross-section for the local mixture.
      allocate(sig_mix_cell(1648,nptot),sig_int(nptot-1), stat=ierr)
      sig_mix_cell = 0d0
!$OMP PARALLEL DO PRIVATE(i,m) SCHEDULE(guided)
          do i=1,1648
              do m=1,nptot
              sig_mix_cell(i,m) = dot_product(fk,sig(:,i,m))
              enddo
          enddo
!$OMP END PARALLEL DO

      !!! Compute the Rossland mean cross-section by integrating over variable v (mesh equally spaced in v).
      sig_Ross = 0d0
!$OMP PARALLEL DO PRIVATE(i,m,sig_int) SCHEDULE(guided)
        do i=1,1648
              sig_int(1:nptot-1) = (1/sig_mix_cell(i,1:nptot-1) + 1/sig_mix_cell(i,2:nptot))/2. * dv !inv_sig_mix_cell(ii,jj,1:nptot-1) + inv_sig_mix_cell(ii,jj,2:nptot)
              do m=1, nptot-1
                sig_Ross(i) = sig_Ross(i) + sig_int(m)
              enddo
        enddo
!$OMP END PARALLEL DO


      logT_pcg   = 0.025*ite
      logRho_pcg = 0.25*jne + log10(amu_mix_cell) - log10(epa_mix_cell) - logNa
      !if (j==1) allocate(lkap_face_pcg(ngp,1648),stat=ierr)
      lkap_face_pcg =  log10_bohr_radius_sqr - log_amu_mix_cell - log10(sig_Ross)
!$OMP PARALLEL DO PRIVATE(i,ke,m,gam) SCHEDULE(guided)
        do i=imin,imax
          do ke=3,nel
            gam = 0d0
            do m=1, nptot-1
            gam = gam + ((eumesh(ke,i,m)/ sig_mix_cell(i,m))+(eumesh(ke,i,m+1)/ sig_mix_cell(i,m+1)))
            enddo
           gamma_k(ke,i) = gam/2. * dv
          enddo
        enddo
!$OMP END PARALLEL DO

        where (gamma_k < 0.)
          gamma_k = 10d-30
        endwhere

        deallocate(sig_mix_cell,sig_int)

        lgamm_pcg = log10(gamma_k)
      end subroutine compute_gamma_grid


      !!! Compute g_rad from precomputeds for grid for gamma and kappa_Ross for the entire OP mono data.
      subroutine compute_grad_fast(k,fk, logT_face, logRho_face, l, r,lgrad_cell, ierr,ite,jne,epatom,amamu,logT_pcg,logRho_pcg,lgamm_pcg,lkap_face_pcg)

      use cubic_interpolator, only: interpolator

      integer, intent(inout) :: ierr
      integer, intent(in) :: k
      real(dp), intent(in) :: fk(:), logT_face, logRho_face, l, r
      real(dp), intent(in) :: logT_pcg(1648), logRho_pcg(1648)!, lgamm_pcg(:,:), lkap_face_pcg(:)
      integer :: nel, nptot
      parameter(nel = 17, nptot = 10000) !number of elements and number of u-mesh points.
      real(dp), intent(in) :: lgamm_pcg(nel,1648), lkap_face_pcg(1648)
      real(dp), intent(out) :: lgrad_cell(nel)

      integer , pointer :: ite(:),jne(:)
      real(dp), pointer :: epatom(:,:),amamu(:)

      integer :: n, ke, nz, id, m, ik, i

      real(dp):: epa_mix_cell(1648), amu_mix_cell ! Number of electrons per atom, mean molecular weight, density and temperature as a function of ite (temp index) and jne (density index) from the OP mono data.
      real(dp) :: delta(1648)
      integer ::  eid(nel)
      real(dp) :: log10_bohr_radius_sqr = -16.55280
      real(dp) :: lgamm_cell(nel) !interpolated log kappa Rossland and log gamma_k in each cell (k = element index).
      real(dp) :: logNa = 23.779750912481397 !log10_cr(6.0221409d23) Avogadro's number
      real(dp) :: dv = (1.0552976117319748 - 0.00010565516589892675)/nptot !v(u(1)) - v(u(nptot))/nptot
      real(dp) :: mH

      !!!! For interpolator.
      integer ::  delta_min_idx
      real(dp) :: lkap_Ross(4,4),sig_Ross(4,4) ! Rossland cross section, log kappa, gamma_k, log gamma_k for interpolation.
      real(dp) :: sf, flux !'scale factor' for g_rad, local flux.
      real(dp), parameter :: pi = 3.141592653589793239
      real(dp), parameter :: log_c = 10.476820702927927 !log10_cr(dble(299792458e2)) c = speed of light

      integer :: ii, jj, ite_min, jne_min, ii_min, jj_min, ite_step, jne_step
      integer :: ite_i, jne_i, dite, djne, i_grid(4,4)
      real(dp) :: logT_min, logRho_min, logT_grid(4,4),logRho_grid(4,4),lgamm(4,4,nel)
      integer ::  offset1, offset2, tries, missing_point(4,4)
      real(dp) :: log_amu_mix_cell, lkap_ross_cell
      integer :: imin, imax
      real(dp) :: logT(1648), logRho(1648)!, ite_grid(4,4), jne_grid(4,4)
      logical :: retry, do_difficult_point
      type(interpolator) :: rossl_interpolator
      type(interpolator) :: gaml_interpolators(nel)

      ierr = 0

      imin = 1
      imax = 1648



      amu_mix_cell = dot_product(fk,amamu)

      logT   = logT_pcg !0.025*ite
      logRho = logRho_pcg !0.25*jne + log10_cr(amu_mix_cell) - log10(epa_mix_cell) - logNa

      !!! Compute an estimated temperature range.
      imin = 1
      imax = 1648
      ite_step = 2
      jne_step = 2

      delta = sqrt((logRho - logRho_face)*(logRho - logRho_face)/0.25 +(logT-logT_face)*(logT-logT_face)/0.025)

      delta_min_idx = MINLOC(delta, DIM=1)
      ite_min   = ite(delta_min_idx)
      jne_min   = jne(delta_min_idx)
      logT_min   = logT(delta_min_idx)
      logRho_min = logRho(delta_min_idx)
      if (logT_min <= logT_face .and. logRho_min <= logRho_face) then
        ii_min = 2
        jj_min = 2
      else if (logT_min < logT_face .and. logRho_min > logRho_face) then
        ii_min = 2
        jj_min = 3
      else if (logT_min > logT_face .and. logRho_min < logRho_face) then
        ii_min = 3
        jj_min = 2
      else if (logT_min > logT_face .and. logRho_min > logRho_face) then
        ii_min = 3
        jj_min = 3
      endif


      offset1 = 0
      offset2 = 0
      missing_point = 1
      tries = 0
      retry = .true.
      do while (retry) !If (T,Rho)_cell lies to close to the edge of the OP mono grid (in Rho), shift the 4x4 grid up/down by 1 in jne.
      missing_point = 1
      retry = .false.
      do jj=1,4
            do ii=1,4
              dite = (ii - ii_min + offset1)*ite_step
              djne = (jj - jj_min + offset2)*jne_step
             do i =imin,imax
               ite_i = ite(i)
               jne_i = jne(i)

                if (ite_i .eq. ite_min +  dite .and. jne_i .eq. jne_min + djne) THEN
                  logT_grid(ii,jj) = logT(i)
                  logRho_grid(ii,jj) = logRho(i)
                  i_grid(ii,jj) = i
                  missing_point(ii,jj) = 0

                endif
              enddo
            enddo
      enddo

      if (SUM(missing_point) > 0) then
        retry = .true.
        if (SUM(missing_point(2:4,1:3)) == 0) then
          offset1 = offset1 + 1
          offset2 = offset2 - 1
        else if (SUM(missing_point(1:3,1:3)) == 0) then
          offset1 = offset1 - 1
          offset2 = offset2 - 1
        else if (SUM(missing_point(2:4,2:4)) == 0) then
          offset1 = offset1 + 1
          offset2 = offset2 + 1
        else if (SUM(missing_point(1:3,2:4)) == 0) then
          offset1 = offset1 - 1
          offset2 = offset2 + 1
        else
          if (ii_min == 3 .and. jj_min == 2) then
            offset1 = offset1 + 2
            offset2 = offset2 - 2
          else if (ii_min == 2 .and. jj_min == 3) then
            offset1 = offset1 - 2
            offset2 = offset2 + 2
          else if (ii_min == 2 .and. jj_min == 2) then
            offset1 = offset1 - 2
            offset2 = offset2 - 2
          else if (ii_min == 3 .and. jj_min == 3) then
            offset1 = offset1 + 2
            offset2 = offset2 - 1
          endif
        endif
      endif

      if (tries > 2) THEN ! To prevent loop from getting stuck.
       do_difficult_point = .false.
       if (ite_min == 226 .and. jne_min == 90 .and. ii_min == 2 .and. jj_min == 3) do_difficult_point = .true.
       if (ite_min == 226 .and. jne_min == 88 .and. ii_min == 2 .and. jj_min == 2) do_difficult_point = .true.
       if (ite_min == 230 .and. jne_min == 88 .and. ii_min == 3 .and. jj_min == 2) do_difficult_point = .true.
       if (ite_min == 230 .and. jne_min == 90 .and. ii_min == 3 .and. jj_min == 3) do_difficult_point = .true.

       if (do_difficult_point) then
          do jj=1,4
            do ii =1,4
              dite = (ii - ii_min)*2*ite_step
              djne = (jj - jj_min - 1)*jne_step
              do i=imin,imax
                ite_i = ite(i)
                jne_i = jne(i)
                if (ite_i == ite_min + dite .and. jne_i == jne_min + djne) then
                  logT_grid(ii,jj) = logT(i)
                  logRho_grid(ii,jj) = logRho(i)
                  i_grid(ii,jj) = i

                endif
              enddo
             enddo
            enddo
            retry = .false.
        else

        write(*,*) 'Cannot find points for interpolation compute_grad_fast', k, ite_min, jne_min, logT_min, logRho_min, logT_face, logRho_face, missing_point, ii_min, jj_min, offset1, offset2,imin, imax
        ierr = 1
        return
        endif
      endif
      tries = tries + 1
      enddo

      do jj=1,4
        do ii=1,4
          ik = i_grid(ii,jj)
          do ke=3,nel
            lgamm(ii,jj,ke) = lgamm_pcg(ke,ik)
          enddo
          lkap_Ross(ii,jj) = lkap_face_pcg(ik)
        enddo
      enddo


      call rossl_interpolator% initialize()
      do jj = 1, 4
         do ii = 1, 4
            call rossl_interpolator% add_point(logT_grid(ii,jj), logRho_grid(ii,jj), lkap_Ross(ii,jj))
         enddo
      enddo

      lkap_ross_cell  = rossl_interpolator% evaluate(logT_face,logRho_face)

      do ke = 3, nel
         call gaml_interpolators(ke)% initialize()
         do jj = 1, 4
            do ii = 1, 4
               call gaml_interpolators(ke)% add_point(logT_grid(ii,jj), logRho_grid(ii,jj), lgamm(ii,jj,ke))
            enddo
         enddo

        lgamm_cell(ke) = gaml_interpolators(ke)% evaluate(logT_face, logRho_face) !cell
      enddo
      lgamm_cell(1) = -30d0
      lgamm_cell(2) = -30d0

      flux  = l / (4d0*pi*r*r)
      sf = lkap_ross_cell  - log_c + log10(flux) + log10(amu_mix_cell)
      lgrad_cell = sf + lgamm_cell - log10(amamu)
      lgrad_cell(1) = -30.
      lgrad_cell(2) = -30.

      end subroutine compute_grad_fast


      !!! Compute Rossland opacity and derivatives for a single cell.
      subroutine compute_kappa(k,fk, logT_cntr, logRho_cntr, lkap_ross_cell, dlnkap_rad_dlnT, dlnkap_rad_dlnRho, ierr,ite,jne,epatom,amamu,sig,logT_grid,logRho_grid,lkap_Ross)
        ! OP mono data for: H, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Cr, Mn, Fe, and Ni.
      use chem_def, only: chem_isos, ih1, ihe3, ihe4, ic12, in14, io16, ine20, ina23,img24, ial27, isi28, is32, iar40, ica40, icr52, imn55, ife56, ini58
      use interp_2d_lib_sg
      use cubic_interpolator, only: interpolator

      integer, intent(inout) :: ierr
      integer, intent(in) :: k
      real(dp), intent(in) :: fk(:), logT_cntr, logRho_cntr
      integer :: nel, nptot
      parameter(nel = 17, nptot = 10000) !number of elements and number of u-mesh points.
      real(dp), intent(out) :: lkap_ross_cell, dlnkap_rad_dlnT, dlnkap_rad_dlnRho
      real(dp), intent(out), dimension(4,4) :: logT_grid, logRho_grid, lkap_Ross


      integer , pointer :: ite(:),jne(:)
      real(dp), pointer :: sig(:,:,:)
      real(dp), pointer :: epatom(:,:),amamu(:)

      integer :: n, ke, nz, id, m, ik, i

      real(dp):: epa_mix_cell(1648), amu_mix_cell, logRho(1648),logT(1648) ! Number of electrons per atom, mean molecular weight, density and temperature as a function of ite (temp index) and jne (density index) from the OP mono data.
      real(dp) :: delta(1648)
      real(dp) :: log10_bohr_radius_sqr = -16.55280
      real(dp) :: lgamm_cell(nel) !interpolated log kappa Rossland and log gamma_k in each cell (k = element index).
      real(dp) :: logNa = 23.779750912481397 !log10_cr(6.0221409d23) Avogadro's number
      real(dp) :: dv = (1.0552976117319748 - 0.00010565516589892675)/nptot !v(u(1)) - v(u(nptot))/nptot
      real(dp) :: mH


      integer ::  delta_min_idx
      real(dp) :: sig_Ross(4,4)!,lkap_Ross(4,4), ! Rossland cross section, log kappa, gamma_k, log gamma_k for interpolation.
      real(dp) :: sf, flux !'scale factor' for g_rad, local flux.
      real(dp), parameter :: pi = 3.141592653589793239
      real(dp), parameter :: log_c = 10.476820702927927 !log10_cr(dble(299792458e2)) c = speed of light

      real, allocatable :: sig_mix_cell(:,:,:),sig_int(:)
      integer :: ii, jj, ite_min, jne_min, ii_min, jj_min, ite_step, jne_step
      integer :: ite_i, jne_i, dite, djne, i_grid(4,4)
      real(dp) :: logT_min, logRho_min!, logT_grid(4,4), logRho_grid(4,4)
      integer :: offset1, offset2, tries, missing_point(4,4)
      real(dp) :: log_amu_mix_cell
      integer :: imin, imax
      real(dp) :: sig_grid(nel,4,4,nptot)
      logical :: retry

      type(interpolator) :: rossl_interpolator

      ierr = 0

      imin = 1 !-1
      imax = 1648 !-1
      ite_step = 2
      jne_step = 2

      !!! Compute the number of electrons per atom for the local mixture.
      epa_mix_cell = 0d0
        do i=imin,imax
          epa_mix_cell(i) = dot_product(fk,epatom(1:nel,i))
        enddo


      amu_mix_cell = dot_product(fk,amamu)

      mH = chem_isos% W(ih1) * 1.660538782d-24
      log_amu_mix_cell = log10(amu_mix_cell * mH)
      !logRho = 0.25*jne + log_amu_mix_cell - log10(epa_mix_cell)
      logRho = 0.25*jne + log10(amu_mix_cell) - log10(epa_mix_cell) -logNa

      logT   = 0.025*ite

      !!! Select nearest points in logT,logRho for interpolation.
      !!! First, find the nearest OP data point, and check which of the four possible positionings this minimum has wrt to (T,Rho)_cell.
      !!! Acquire the remaining 15 points of the 4x4 grid, where (T,Rho)_cell is located in the inner square.


      delta = sqrt((logRho - logRho_cntr)*(logRho - logRho_cntr)/0.25 +(logT-logT_cntr)*(logT-logT_cntr)/0.025)

      delta_min_idx = MINLOC(delta, DIM=1)
      ite_min   = ite(delta_min_idx)
      jne_min   = jne(delta_min_idx)
      logT_min   = logT(delta_min_idx)
      logRho_min = logRho(delta_min_idx)
      if (logT_min <= logT_cntr .and. logRho_min <= logRho_cntr) then
        ii_min = 2
        jj_min = 2
      else if (logT_min < logT_cntr .and. logRho_min > logRho_cntr) then
        ii_min = 2!3
        jj_min = 3!2
      else if (logT_min > logT_cntr .and. logRho_min < logRho_cntr) then
        ii_min = 3
        jj_min = 2
      else if (logT_min > logT_cntr .and. logRho_min > logRho_cntr) then
        ii_min = 3
        jj_min = 3
      endif

      offset1 = 0
      offset2 = 0
      missing_point = 1
      tries = 0
      retry = .true.
      do while (retry) !If (T,Rho)_cell lies to close to the edge of the OP mono grid (in Rho), shift the 4x4 grid up/down by 1 in jne.
      missing_point = 1
      retry = .false.
      do jj=1,4
            do ii=1,4
              dite = (ii - ii_min + offset1)*ite_step
              djne = (jj - jj_min + offset2)*jne_step
             do i =imin,imax
               ite_i = ite(i)
               jne_i = jne(i)

                if (ite_i .eq. ite_min +  dite .and. jne_i .eq. jne_min + djne) THEN
                  logT_grid(ii,jj) = logT(i)
                  logRho_grid(ii,jj) = logRho(i)
                  i_grid(ii,jj) = i

                  missing_point(ii,jj) = 0
                endif
              enddo
            enddo
      enddo

      if (SUM(missing_point) > 0) then
        retry = .true.
        if (SUM(missing_point(2:4,1:3)) == 0) then
          offset1 = offset1 + 1
          offset2 = offset2 - 1
        else if (SUM(missing_point(1:3,1:3)) == 0) then
          offset1 = offset1 - 1
          offset2 = offset2 - 1
        else if (SUM(missing_point(2:4,2:4)) == 0) then
          offset1 = offset1 + 1
          offset2 = offset2 + 1
        else if (SUM(missing_point(1:3,2:4)) == 0) then
          offset1 = offset1 - 1
          offset2 = offset2 + 1
        else
          if (ii_min == 3 .and. jj_min == 2) then
            offset1 = offset1 + 2
            offset2 = offset2 - 2
          else if (ii_min == 2 .and. jj_min == 3) then
            offset1 = offset1 - 2
            offset2 = offset2 + 2
          else if (ii_min == 2 .and. jj_min == 2) then
            offset1 = offset1 - 2
            offset2 = offset2 - 2
          else if (ii_min == 3 .and. jj_min == 3) then
            offset1 = offset1 + 2
            offset2 = offset2 - 1
          endif
        endif
      endif



      if (tries > 4) THEN ! To prevent loop from getting stuck.
        write(*,*) 'Cannot find points for interpolation compute_kappa', k, ite_min, jne_min, logT_min, logRho_min,logT_cntr, logRho_cntr, missing_point, ii_min, jj_min, offset1, offset2,imin, imax, log_amu_mix_cell
        ierr = 1
        return
      endif
      tries = tries + 1
      enddo


      !!! Compute the monochromatic cross-section for the local mixture.
      allocate(sig_mix_cell(4,4,nptot),sig_int(nptot-1), stat=ierr)


      sig_mix_cell = 0d0
      do jj=1,4
        do ii=1,4
              ik = i_grid(ii,jj)
              do m=1,nptot
                 sig_mix_cell(ii,jj,m) = dot_product(fk,sig(:,ik,m))
              enddo
        enddo
      enddo

      !!! Compute the Rossland mean cross-section by integrating over variable v (mesh equally spaced in v).
      sig_Ross = 0d0
      do jj=1,4
        do ii=1,4
              sig_int = (1/sig_mix_cell(ii,jj,1:nptot-1) + 1/sig_mix_cell(ii,jj,2:nptot))/2. * dv !(1:nptot-1) inv_sig_mix_cell(ii,jj,1:nptot-1) + inv_sig_mix_cell(ii,jj,2:nptot)
              do m=1, nptot-1
                sig_Ross(ii,jj) = sig_Ross(ii,jj) + sig_int(m)
              enddo
        enddo
      enddo

      lkap_Ross =  log10_bohr_radius_sqr - log_amu_mix_cell - log10(sig_Ross)


      call rossl_interpolator% initialize()
      do jj = 1, 4
         do ii = 1, 4
            call rossl_interpolator% add_point(logT_grid(ii,jj), logRho_grid(ii,jj), lkap_Ross(ii,jj))
         enddo
      enddo

      lkap_ross_cell  = rossl_interpolator% evaluate(logT_cntr,logRho_cntr)
      dlnkap_rad_dlnT = rossl_interpolator% evaluate_deriv(logT_cntr,logRho_cntr, .true., .false.)
      dlnkap_rad_dlnRho = rossl_interpolator% evaluate_deriv(logT_cntr, logRho_cntr, .false., .true.)


      deallocate(sig_mix_cell,sig_int)

      end subroutine compute_kappa



      !!! Compute the Rossland opacity for the entire OP mono data set given a mixture.
      subroutine compute_kappa_grid(fk,lkap_ross_pcg, ierr,ite,jne,epatom,amamu,sig)
        ! OP mono data for: H, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Cr, Mn, Fe, and Ni.
      use chem_def, only: chem_isos, ih1, ihe3, ihe4, ic12, in14, io16, ine20, ina23,img24, ial27, isi28, is32, iar40, ica40, icr52, imn55, ife56, ini58
      use interp_2d_lib_sg
      use cubic_interpolator, only: interpolator

      integer, intent(inout) :: ierr
      real(dp), intent(in) :: fk(:)
      integer :: nel, nptot
      parameter(nel = 17, nptot = 10000) !number of elements and number of u-mesh points.
      real(dp), pointer, intent(out) :: lkap_ross_pcg(:)!, fk_pcg(:)
      real(dp) :: logT_pcg(1648), logRho_pcg(1648)

      integer , pointer :: ite(:),jne(:)
      real(dp), pointer :: sig(:,:,:)!, ak_f(:,:)
      real(dp), pointer :: epatom(:,:),amamu(:)
      integer :: imin, imax

      integer :: n, ke, nz, id, m, ik, i

      real(dp):: epa_mix_cell(1648), amu_mix_cell ! Number of electrons per atom, mean molecular weight, density and temperature as a function of ite (temp index) and jne (density index) from the OP mono data.
      real(dp) :: log10_bohr_radius_sqr = -16.55280
      real(dp) :: lgamm_cell(nel) !interpolated log kappa Rossland and log gamma_k in each cell (k = element index).
      real(dp) :: logNa = 23.779750912481397 !log10_cr(6.0221409d23) Avogadro's number
      real(dp) :: dv = (1.0552976117319748 - 0.00010565516589892675)/nptot !v(u(1)) - v(u(nptot))/nptot
      real(dp) :: mH

      !!!! For interpolator.
      integer ::  delta_min_idx

      real(dp) :: sig_Ross(1648) ! Rossland cross section, log kappa, gamma_k, log gamma_k for interpolation.
      real(dp), parameter :: pi = 3.141592653589793239
      real(dp), parameter :: log_c = 10.476820702927927 !log10_cr(dble(299792458e2)) c = speed of light

      real, allocatable :: sig_mix_cell(:,:),sig_int(:)
      real(dp) :: log_amu_mix_cell

      type(interpolator) :: rossl_interpolator

      ierr = 0

      !!! Compute an estimated temperature range.
      imin = 1
      imax = 1648

      !!! Compute the number of electrons per atom for the local mixture.
      epa_mix_cell = 0d0
      do i=imin,imax
          epa_mix_cell(i) = dot_product(fk,epatom(1:nel,i))
      enddo

      amu_mix_cell = dot_product(fk,amamu)

      mH = chem_isos% W(ih1) * 1.660538782d-24
      log_amu_mix_cell = log10(amu_mix_cell * mH)

      !!! Compute the monochromatic cross-section for the local mixture.
      allocate(sig_mix_cell(1648,nptot),sig_int(nptot-1), stat=ierr)
      sig_mix_cell = 0d0
!$OMP PARALLEL DO PRIVATE(i,m) SCHEDULE(guided)

          do i=1,1648
              do m=1,nptot
              sig_mix_cell(i,m) = dot_product(fk,sig(:,i,m))
              enddo
          enddo
!$OMP END PARALLEL DO

      !!! Compute the Rossland mean cross-section by integrating over variable v (mesh equally spaced in v).
      sig_Ross = 0
!$OMP PARALLEL DO PRIVATE(i,m,sig_int) SCHEDULE(guided)
        do i=1,1648
              sig_int(1:nptot-1) = (1/sig_mix_cell(i,1:nptot-1) + 1/sig_mix_cell(i,2:nptot))/2. * dv !inv_sig_mix_cell(ii,jj,1:nptot-1) + inv_sig_mix_cell(ii,jj,2:nptot)
              do m=1, nptot-1
                sig_Ross(i) = sig_Ross(i) + sig_int(m)
              enddo
        enddo
!$OMP END PARALLEL DO

      deallocate(sig_mix_cell,sig_int)

      logT_pcg   = 0.025*ite
      logRho_pcg = 0.25*jne + log10(amu_mix_cell) - log10(epa_mix_cell) - logNa
      allocate(lkap_ross_pcg(1648),stat=ierr)
      lkap_ross_pcg =  log10_bohr_radius_sqr - log_amu_mix_cell - log10(sig_Ross)

      end subroutine compute_kappa_grid


      !!! Routine to compute Rossland opacity from a precomputed grid of the entire OP mono data.
      subroutine compute_kappa_fast(k,fk, logT_cntr, logRho_cntr,lkap_ross_cell, dlnkap_rad_dlnT, dlnkap_rad_dlnRho, ierr,ite,jne,epatom,amamu,sig,lkap_ross_pcg)
        ! OP mono data for: H, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Cr, Mn, Fe, and Ni.
      use chem_def, only: chem_isos, ih1, ihe3, ihe4, ic12, in14, io16, ine20, ina23,img24, ial27, isi28, is32, iar40, ica40, icr52, imn55, ife56, ini58
      use cubic_interpolator, only: interpolator

      integer, intent(inout) :: ierr
      integer, intent(in) :: k
      real(dp), intent(in) :: fk(:), logT_cntr, logRho_cntr
      real(dp), pointer, intent(in) :: lkap_ross_pcg(:)
      integer :: nel, nptot
      parameter(nel = 17, nptot = 10000) !number of elements and number of u-mesh points.
      real(dp), intent(out) :: lkap_ross_cell, dlnkap_rad_dlnT, dlnkap_rad_dlnRho
      integer , pointer :: ite(:),jne(:)
      real(dp), pointer :: sig(:,:,:)
      real(dp), pointer :: epatom(:,:),amamu(:)

      integer :: n, ke, nz, id, m, ik, i

      real(dp):: epa_mix_cell(1648), amu_mix_cell ! Number of electrons per atom, mean molecular weight, density and temperature as a function of ite (temp index) and jne (density index) from the OP mono data.
      real(dp) :: delta(1648)
      real(dp) :: log10_bohr_radius_sqr = -16.55280
      real(dp) :: lgamm_cell(nel) !interpolated log kappa Rossland and log gamma_k in each cell (k = element index).
      real(dp) :: logNa = 23.779750912481397 !log10(6.0221409d23) Avogadro's number
      real(dp) :: dv = (1.0552976117319748 - 0.00010565516589892675)/nptot !v(u(1)) - v(u(nptot))/nptot
      real(dp) :: mH

      !!!! For interpolator.
      integer ::  delta_min_idx
      real(dp) :: lkap_Ross(4,4),sig_Ross(4,4) ! Rossland cross section, log kappa, gamma_k, log gamma_k for interpolation.
      real(dp) :: sf, flux !'scale factor' for g_rad, local flux.
      real(dp), parameter :: pi = 3.141592653589793239
      real(dp), parameter :: log_c = 10.476820702927927 !log10_cr(dble(299792458e2)) c = speed of light

      integer :: ii, jj, ite_min, jne_min, ii_min, jj_min, ite_step, jne_step
      integer :: ite_i, jne_i, dite, djne, i_grid(4,4)
      real(dp) :: logT_min, logRho_min, logT_grid(4,4), logRho_grid(4,4)
      integer ::  offset1, offset2, tries, missing_point(4,4)
      real(dp) :: log_amu_mix_cell
      integer :: imin, imax
      real(dp) :: logT(1648), logRho(1648)!, ite_grid(4,4), jne_grid(4,4)
      logical :: retry, do_difficult_point

      type(interpolator) :: rossl_interpolator

      ierr = 0

      imin = 1
      imax = 1648

      !!! Compute the number of electrons per atom for the local mixture.
      epa_mix_cell = 0d0
        do i=imin,imax
          epa_mix_cell(i) = dot_product(fk,epatom(:,i))
        enddo

      amu_mix_cell = dot_product(fk,amamu)

      logT   = 0.025*ite
      logRho = 0.25*jne + log10(amu_mix_cell) - log10(epa_mix_cell) - logNa

      imin = 1
      imax = 1648
      ite_step = 2
      jne_step = 2

      !!! Select nearest points in logT,logRho for interpolation.
      !!! First, find the nearest OP data point, and check which of the four possible positionings this minimum has wrt to (T,Rho)_cell.
      !!! Acquire the remaining 15 points of the 4x4 grid, where (T,Rho)_cell is located in the inner square.

      delta = sqrt((logRho - logRho_cntr)*(logRho - logRho_cntr)/0.25 +(logT-logT_cntr)*(logT-logT_cntr)/0.025)

      delta_min_idx = MINLOC(delta, DIM=1)
      !delta_min(1)     = MINVAL(delta)(1)
      ite_min   = ite(delta_min_idx)
      jne_min   = jne(delta_min_idx)
      logT_min   = logT(delta_min_idx)
      logRho_min = logRho(delta_min_idx)
      if (logT_min <= logT_cntr .and. logRho_min <= logRho_cntr) then
        ii_min = 2
        jj_min = 2
      else if (logT_min < logT_cntr .and. logRho_min > logRho_cntr) then
        ii_min = 2!3
        jj_min = 3!2
      else if (logT_min > logT_cntr .and. logRho_min < logRho_cntr) then
        ii_min = 3
        jj_min = 2
      else if (logT_min > logT_cntr .and. logRho_min > logRho_cntr) then
        ii_min = 3
        jj_min = 3
      endif


      offset1 = 0
      offset2 = 0
      missing_point = 1
      tries = 0
      retry = .true.
      do while (retry) !If (T,Rho)_cell lies to close to the edge of the OP mono grid (in Rho), shift the 4x4 grid up/down by 1 in jne.
      missing_point = 1
      retry = .false.
      do jj=1,4
            do ii=1,4
              dite = (ii - ii_min + offset1)*ite_step !(ii - ii_min)*ite_step + offset2
              djne = (jj - jj_min + offset2)*jne_step !(jj - jj_min)*jne_step + offset1
             do i =imin,imax
               ite_i = ite(i)
               jne_i = jne(i)

                if (ite_i .eq. ite_min +  dite .and. jne_i .eq. jne_min + djne) THEN
                  logT_grid(ii,jj) = logT(i)
                  logRho_grid(ii,jj) = logRho(i)
                  i_grid(ii,jj) = i
                  missing_point(ii,jj) = 0

                endif
              enddo
            enddo
      enddo

      if (SUM(missing_point) > 0) then
        retry = .true.
        if (SUM(missing_point(2:4,1:3)) == 0) then
          offset1 = offset1 + 1
          offset2 = offset2 - 1
        else if (SUM(missing_point(1:3,1:3)) == 0) then
          offset1 = offset1 - 1
          offset2 = offset2 - 1
        else if (SUM(missing_point(2:4,2:4)) == 0) then
          offset1 = offset1 + 1
          offset2 = offset2 + 1
        else if (SUM(missing_point(1:3,2:4)) == 0) then
          offset1 = offset1 - 1
          offset2 = offset2 + 1
        else
          if (ii_min == 3 .and. jj_min == 2) then
            offset1 = offset1 + 2
            offset2 = offset2 - 2
          else if (ii_min == 2 .and. jj_min == 3) then
            offset1 = offset1 - 2
            offset2 = offset2 + 2
          else if (ii_min == 2 .and. jj_min == 2) then
            offset1 = offset1 - 2
            offset2 = offset2 - 2
          else if (ii_min == 3 .and. jj_min == 3) then
            offset1 = offset1 + 2
            offset2 = offset2 - 1
          endif

        endif
      endif

      if (tries > 2) THEN ! To prevent loop from getting stuck.
        do_difficult_point = .false.
        if (ite_min == 226 .and. jne_min == 90 .and. ii_min == 2 .and. jj_min == 3) do_difficult_point = .true.
        if (ite_min == 226 .and. jne_min == 88 .and. ii_min == 2 .and. jj_min == 2) do_difficult_point = .true.
        if (ite_min == 230 .and. jne_min == 88 .and. ii_min == 3 .and. jj_min == 2) do_difficult_point = .true.
        if (ite_min == 230 .and. jne_min == 90 .and. ii_min == 3 .and. jj_min == 3) do_difficult_point = .true.

        if (do_difficult_point) then
           do jj=1,4
             do ii =1,4
               dite = (ii - ii_min)*2*ite_step
               djne = (jj - jj_min - 1)*jne_step
               do i=imin,imax
                 ite_i = ite(i)
                 jne_i = jne(i)
                 if (ite_i == ite_min + dite .and. jne_i == jne_min + djne) then
                   logT_grid(ii,jj) = logT(i)
                   logRho_grid(ii,jj) = logRho(i)
                   i_grid(ii,jj) = i
                 endif
               enddo
              enddo
             enddo
             retry = .false.
         else
        write(*,*) 'Cannot find points for interpolation compute_kappa_fast', k, ite_min, jne_min, logT_min, logRho_min,logT_cntr, logRho_cntr, missing_point, ii_min, jj_min, offset1, offset2,imin, imax
        ierr = 1
        return
        !stop
        endif
      endif
      tries = tries + 1
      enddo

      do jj=1,4
        do ii=1,4
          ik = i_grid(ii,jj)
          lkap_Ross(ii,jj) = lkap_ross_pcg(ik)
        enddo
      enddo

      call rossl_interpolator% initialize()
      do jj = 1, 4
         do ii = 1, 4
            call rossl_interpolator% add_point(logT_grid(ii,jj), logRho_grid(ii,jj), lkap_Ross(ii,jj))
         enddo
      enddo

      lkap_ross_cell  = rossl_interpolator% evaluate(logT_cntr,logRho_cntr)
      dlnkap_rad_dlnT = rossl_interpolator% evaluate_deriv(logT_cntr, logRho_cntr, .true., .false.)
      dlnkap_rad_dlnRho = rossl_interpolator% evaluate_deriv(logT_cntr, logRho_cntr, .false., .true.)

      end subroutine compute_kappa_fast


      !!! Interpolate Rossland opacity and derivatives from a 4x4 grid computed last time step.
      subroutine interpolate_kappa(k,logT_cntr, logRho_cntr,lkap_ross_cell, dlnkap_rad_dlnT, dlnkap_rad_dlnRho, ierr,logT_grid, logRho_grid, lkap_grid)

      use cubic_interpolator, only: interpolator

      integer, intent(inout) :: ierr
      integer, intent(in) :: k
      real(dp), intent(in) :: logT_cntr, logRho_cntr
      real(dp), intent(in) :: logT_grid(4,4), logRho_grid(4,4), lkap_grid(4,4)
      real(dp), intent(out) :: lkap_ross_cell, dlnkap_rad_dlnT, dlnkap_rad_dlnRho

      integer :: ii,jj

      type(interpolator) :: rossl_interpolator

      ierr = 0

      call rossl_interpolator% initialize()
      do jj = 1, 4
         do ii = 1, 4
            call rossl_interpolator% add_point(logT_grid(ii,jj), logRho_grid(ii,jj), lkap_grid(ii,jj))
         enddo
      enddo

      lkap_ross_cell  = rossl_interpolator% evaluate(logT_cntr,logRho_cntr)
      dlnkap_rad_dlnT = rossl_interpolator% evaluate_deriv(logT_cntr,logRho_cntr, .true., .false.)
      dlnkap_rad_dlnRho = rossl_interpolator% evaluate_deriv(logT_cntr,logRho_cntr, .false., .true.)

      end subroutine interpolate_kappa
      end module op_eval_mombarg
