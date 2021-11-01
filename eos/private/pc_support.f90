   module pc_support
   use utils_lib, only: is_bad, mesa_error, mv, switch_str
   use eos_def
   use const_def, only: ln10, mesa_temp_caches_dir, one_sixth
   use math_lib

   implicit none
   
   public :: get_FITION9, get_EXCOR7, get_FSCRliq8
   private

   ! the file data for EXCOR7
   integer, parameter :: jFXC = 1
   integer, parameter :: jUXC = 2
   integer, parameter :: jPXC = 3
   integer, parameter :: jCVXC = 4
   integer, parameter :: jSXC = 5
   integer, parameter :: jPDTXC = 6
   integer, parameter :: jPDRXC = 7
   integer, parameter :: nvals_EXCOR7 = 7

   ! the file data for FSCRliq8   FSCR, USCR, PSCR, CVSCR, PDTSCR, PDRSCR
   integer, parameter :: iFSCR = 1
   integer, parameter :: iUSCR = 2
   integer, parameter :: iPSCR = 3
   integer, parameter :: iCVSCR = 4
   integer, parameter :: iPDTSCR = 5
   integer, parameter :: iPDRSCR = 6
   integer, parameter :: nvals_FSCRliq8 = 6

   
   contains
   
   
   subroutine get_FITION9(GAMI, &
         FION, dFION_dlnGAMI, &
         UION, dUION_dlnGAMI, &
         PION, dPION_dlnGAMI, &
         CVii, dCVii_dlnGAMI, &
         PDTii, dPDTii_dlnGAMI, &
         PDRii, dPDRii_dlnGAMI, &
         skip, ierr)
      real(dp), intent(in) :: GAMI
      real(dp), intent(out) :: &
         FION, dFION_dlnGAMI, &
         UION, dUION_dlnGAMI, &
         PION, dPION_dlnGAMI, &
         CVii, dCVii_dlnGAMI, &
         PDTii, dPDTii_dlnGAMI, &
         PDRii, dPDRii_dlnGAMI
      logical, intent(out) :: skip
      integer, intent(out) :: ierr

      integer, parameter :: iFION = 1
      integer, parameter :: iUION = 2
      integer, parameter :: iPION = 3
      integer, parameter :: iCVii = 4
      integer, parameter :: iPDTii = 5
      integer, parameter :: iPDRii = 6 
      real(dp) :: lnGAMI_min, lnGAMI_max, dlnGAMI, lnGAMI
      integer :: nlnGAMI, i, j, k_old
      real(dp) :: xk_old, xkp1_old, xk_new, delta
      real(dp), pointer, dimension(:,:) :: f_1, f_2, f_3, f_4, f_5, f_6
      
      include 'formats'

      ierr = 0
      skip = .false.
      
      lnGAMI_min = -2.0d0*ln10
      lnGAMI_max = 6.5d0*ln10
      dlnGAMI = 1d-2*ln10
      nlnGAMI = (lnGAMI_max - lnGAMI_min)/dlnGAMI + 1
      
      if (.not. FITION9_loaded) then
!$OMP CRITICAL(eosPC_support_FITION9_load)
         if (.not. FITION9_loaded) then
            call build_FITION9_data(ierr)
            if (ierr == 0) FITION9_loaded = .true.
         end if
!$OMP END CRITICAL(eosPC_support_FITION9_load)
         if (ierr /= 0) return
      end if
      
      lnGAMI = log(GAMI)
      if (lnGAMI < lnGAMI_min .or. lnGAMI > lnGAMI_max) then
         skip = .true.
         return
      end if
      
      f_1(1:4,1:nlnGAMI) => FITION_data(iFION)% f1
      f_2(1:4,1:nlnGAMI) => FITION_data(iUION)% f1
      f_3(1:4,1:nlnGAMI) => FITION_data(iPION)% f1
      f_4(1:4,1:nlnGAMI) => FITION_data(iCVii)% f1
      f_5(1:4,1:nlnGAMI) => FITION_data(iPDTii)% f1
      f_6(1:4,1:nlnGAMI) => FITION_data(iPDRii)% f1

      xk_new = lnGAMI
      k_old = int((lnGAMI - lnGAMI_min)/dlnGAMI + 1d-4) + 1  
      if (k_old < 1 .or. k_old >= nlnGAMI) then            
         if (k_old < 1) then
            k_old = 1
            xk_old = lnGAMI_min
            xkp1_old = xk_old + dlnGAMI
         else
            k_old = nlnGAMI-1
            xk_old = lnGAMI_min + (k_old-1)*dlnGAMI
            xkp1_old = xk_old + dlnGAMI
         end if            
      else         
         xk_old = lnGAMI_min + (k_old-1)*dlnGAMI
         xkp1_old = xk_old + dlnGAMI
      end if
      
      delta = xk_new - xk_old
      
      FION = &
            f_1(1, k_old) + delta*(f_1(2, k_old)  &
               + delta*(f_1(3, k_old) + delta*f_1(4, k_old)))
      dFION_dlnGAMI = &
            f_1(2, k_old) + 2*delta*(f_1(3, k_old) + 1.5d0*delta*f_1(4, k_old))
      
      UION = &
            f_2(1, k_old) + delta*(f_2(2, k_old)  &
               + delta*(f_2(3, k_old) + delta*f_2(4, k_old)))
      dUION_dlnGAMI = &
            f_2(2, k_old) + 2*delta*(f_2(3, k_old) + 1.5d0*delta*f_2(4, k_old))
      
      PION = &
            f_3(1, k_old) + delta*(f_3(2, k_old)  &
               + delta*(f_3(3, k_old) + delta*f_3(4, k_old)))
      dPION_dlnGAMI = &
            f_3(2, k_old) + 2*delta*(f_3(3, k_old) + 1.5d0*delta*f_3(4, k_old))
      
      CVii = &
            f_4(1, k_old) + delta*(f_4(2, k_old)  &
               + delta*(f_4(3, k_old) + delta*f_4(4, k_old)))
      dCVii_dlnGAMI = &
            f_4(2, k_old) + 2*delta*(f_4(3, k_old) + 1.5d0*delta*f_4(4, k_old))
      
      PDTii = &
            f_5(1, k_old) + delta*(f_5(2, k_old)  &
               + delta*(f_5(3, k_old) + delta*f_5(4, k_old)))
      dPDTii_dlnGAMI = &
            f_5(2, k_old) + 2*delta*(f_5(3, k_old) + 1.5d0*delta*f_5(4, k_old))
      
      PDRii = &
            f_6(1, k_old) + delta*(f_6(2, k_old)  &
               + delta*(f_6(3, k_old) + delta*f_6(4, k_old)))
      dPDRii_dlnGAMI = &
            f_6(2, k_old) + 2*delta*(f_6(3, k_old) + 1.5d0*delta*f_6(4, k_old))
      
      contains
   
         
      subroutine build_FITION9_data(ierr)
         use interp_1d_def, only : pm_work_size
         use interp_1d_lib, only : interp_pm
         integer, intent(out) :: ierr
         type (FITION_Info), pointer :: fi
         real(dp), pointer :: work(:)
         real(dp) :: lnGAMI, GAMI, FION, UION, PION, CVii, PDTii, PDRii
         integer :: j
         include 'formats'
         
         ierr = 0
         !write(*,'(a)') 'eosPC build_FITION9_data'
         
         allocate(FITION9_lnGAMIs(nlnGAMI), work(pm_work_size*nlnGAMI))
         do j=1,FITION_vals
            fi => FITION_data(j)
            allocate(fi% f1(4*nlnGAMI))
            fi% f(1:4, 1:nlnGAMI) => fi% f1(1:4*nlnGAMI)
         end do
         
         do j=1,nlnGAMI-1
            FITION9_lnGAMIs(j) = lnGAMI_min + (j-1)*dlnGAMI
         end do
         FITION9_lnGAMIs(nlnGAMI) = lnGAMI_max
         
         do j=1,nlnGAMI
            lnGAMI = FITION9_lnGAMIs(j)
            GAMI = exp(lnGAMI)
            call FITION9(GAMI,FION,UION,PION,CVii,PDTii,PDRii)
            fi => FITION_data(iFION); fi% f(1,j) = FION
            fi => FITION_data(iUION); fi% f(1,j) = UION
            fi => FITION_data(iPION); fi% f(1,j) = PION
            fi => FITION_data(iCVii); fi% f(1,j) = CVii
            fi => FITION_data(iPDTii); fi% f(1,j) = PDTii
            fi => FITION_data(iPDRii); fi% f(1,j) = PDRii
         end do
         
         do j=1,FITION_vals
            fi => FITION_data(j)
            call interp_pm(FITION9_lnGAMIs, nlnGAMI, &
               fi% f1, pm_work_size, work, 'build_FITION9_data', ierr)
            if (ierr /= 0) return
         end do
         
         deallocate(work)

         !write(*,'(a)') 'done eosPC build_FITION9_data'
         
      end subroutine build_FITION9_data
         
   end subroutine get_FITION9


   subroutine get_FSCRliq8(iZion, RS, GAME, &
         FSCR, dFSCR_dlnRS, dFSCR_dlnGAME, &
         USCR, dUSCR_dlnRS, dUSCR_dlnGAME, &
         PSCR, dPSCR_dlnRS, dPSCR_dlnGAME, &
         CVSCR, dCVSCR_dlnRS, dCVSCR_dlnGAME, &
         PDTSCR, dPDTSCR_dlnRS, dPDTSCR_dlnGAME, &
         PDRSCR, dPDRSCR_dlnRS, dPDRSCR_dlnGAME, &
         skip, ierr) 
      integer, intent(in) :: iZion
      real(dp), intent(in) :: RS, GAME
      real(dp), intent(out) :: &
         FSCR, dFSCR_dlnRS, dFSCR_dlnGAME, &
         USCR, dUSCR_dlnRS, dUSCR_dlnGAME, &
         PSCR, dPSCR_dlnRS, dPSCR_dlnGAME, &
         CVSCR, dCVSCR_dlnRS, dCVSCR_dlnGAME, &
         PDTSCR, dPDTSCR_dlnRS, dPDTSCR_dlnGAME, &
         PDRSCR, dPDRSCR_dlnRS, dPDRSCR_dlnGAME
      logical, intent(out) :: skip
      integer, intent(out) :: ierr
      
      integer :: iRS, jGAME, i
      real(dp) :: lnRS, lnRS0, lnRS1, lnGAME, lnGAME0, lnGAME1
      real(dp), dimension(nvals_FSCRliq8) :: fval, df_dlnRS, df_dlnGAME
      type (eosPC_Support_Info), pointer :: fq
      include 'formats'
      ierr = 0
      skip = .false.
      
      if (iZion < 1 .or. iZion > max_FSCRliq8_Zion) then
         write(*,2) 'invalid value for Z ion in get_FSCRliq8', iZion
         ierr = -1
         return
      end if
      
      fq => FSCRliq8_data(iZion)
      
      if (.not. FSCRliq8_Zion_loaded(iZion)) then
!$OMP CRITICAL(eosPC_support_FSCRliq8_load)
         if (.not. FSCRliq8_Zion_loaded(iZion)) then
            call load_eosPC_support_Info(fq,iZion,ierr)
            if (ierr == 0) FSCRliq8_Zion_loaded(iZion) = .true.
         end if
!$OMP END CRITICAL(eosPC_support_FSCRliq8_load)
         if (ierr /= 0) return
      end if

      ! get interpolation results
      lnRS = log(RS)
      lnGAME = log(GAME)
      if (lnRS < fq% lnRS_min .or. lnRS > fq% lnRS_max .or. &
          lnGAME < fq% lnGAME_min .or. lnGAME > fq% lnGAME_max) then
         skip = .true.
         return
      end if
      
      call Locate_lnRS(lnRS, &
         fq% nlnRS, fq% lnRS_min, fq% dlnRS, &
         iRS, lnRS0, lnRS1)
      call Locate_lnGAME(lnGAME, &
         fq% nlnGAME, fq% lnGAME_min, fq% dlnGAME, &
         jGAME, lnGAME0, lnGAME1)
      
      call Do_Interpolations( &
         1, nvals_FSCRliq8, nvals_FSCRliq8, nvals_FSCRliq8, &
         fq% nlnRS, fq% lnRSs, fq% nlnGAME, fq% lnGAMEs, fq% tbl1, &
         iRS, jGAME, lnRS0, lnRS, lnRS1, lnGAME0, lnGAME, lnGAME1, &
         fval, df_dlnRS, df_dlnGAME, ierr)               
      if (ierr /= 0) then
         write(*,1) 'Do_Interpolations failed in get_FSCRliq8'
         return
      end if
      
      FSCR = fval(iFSCR); dFSCR_dlnRS = df_dlnRS(iFSCR); dFSCR_dlnGAME = df_dlnGAME(iFSCR)
      USCR = fval(iUSCR); dUSCR_dlnRS = df_dlnRS(iUSCR); dUSCR_dlnGAME = df_dlnGAME(iUSCR)
      PSCR = fval(iPSCR); dPSCR_dlnRS = df_dlnRS(iPSCR); dPSCR_dlnGAME = df_dlnGAME(iPSCR)
      CVSCR = fval(iCVSCR); dCVSCR_dlnRS = df_dlnRS(iCVSCR); dCVSCR_dlnGAME = df_dlnGAME(iCVSCR)
      PDTSCR = fval(iPDTSCR); dPDTSCR_dlnRS = df_dlnRS(iPDTSCR); dPDTSCR_dlnGAME = df_dlnGAME(iPDTSCR)
      PDRSCR = fval(iPDRSCR); dPDRSCR_dlnRS = df_dlnRS(iPDRSCR); dPDRSCR_dlnGAME = df_dlnGAME(iPDRSCR)
   
   end subroutine get_FSCRliq8   
   
   
   subroutine get_EXCOR7(RS, GAME, &
         FXC, dFXC_dlnRS, dFXC_dlnGAME, &
         UXC, dUXC_dlnRS, dUXC_dlnGAME, &
         PXC, dPXC_dlnRS, dPXC_dlnGAME, &
         CVXC, dCVXC_dlnRS, dCVXC_dlnGAME, &
         SXC, dSXC_dlnRS, dSXC_dlnGAME, &
         PDTXC, dPDTXC_dlnRS, dPDTXC_dlnGAME, &
         PDRXC, dPDRXC_dlnRS, dPDRXC_dlnGAME, &
         skip, ierr) 
      real(dp), intent(in) :: RS, GAME
      real(dp), intent(out) :: &
         FXC, dFXC_dlnRS, dFXC_dlnGAME, &
         UXC, dUXC_dlnRS, dUXC_dlnGAME, &
         PXC, dPXC_dlnRS, dPXC_dlnGAME, &
         CVXC, dCVXC_dlnRS, dCVXC_dlnGAME, &
         SXC, dSXC_dlnRS, dSXC_dlnGAME, &
         PDTXC, dPDTXC_dlnRS, dPDTXC_dlnGAME, &
         PDRXC, dPDRXC_dlnRS, dPDRXC_dlnGAME
      logical, intent(out) :: skip
      integer, intent(out) :: ierr
      
      integer :: iRS, jGAME, i
      real(dp) :: lnRS, lnRS0, lnRS1, lnGAME, lnGAME0, lnGAME1
      real(dp), dimension(nvals_EXCOR7) :: fval, df_dlnRS, df_dlnGAME
      type (eosPC_Support_Info), pointer :: fq
      include 'formats'
      ierr = 0
      skip = .false.

      fq => EXCOR7_data
      
      if (.not. EXCOR7_table_loaded) then
!$OMP CRITICAL(eosPC_EXCOR7_support_load)
         if (.not. EXCOR7_table_loaded) then
            call load_eosPC_support_Info(fq,0,ierr)
            if (ierr == 0) EXCOR7_table_loaded = .true.
         end if
!$OMP END CRITICAL(eosPC_EXCOR7_support_load)
         if (ierr /= 0) return
      end if

      ! get interpolation results
      lnRS = log(RS)
      lnGAME = log(GAME)
      if (lnRS < fq% lnRS_min .or. lnRS > fq% lnRS_max .or. &
          lnGAME < fq% lnGAME_min .or. lnGAME > fq% lnGAME_max) then
         skip = .true.
         return
      end if
      
      call Locate_lnRS(lnRS, &
         fq% nlnRS, fq% lnRS_min, fq% dlnRS, iRS, lnRS0, lnRS1)
      call Locate_lnGAME(lnGAME, &
         fq% nlnGAME, fq% lnGAME_min, fq% dlnGAME, jGAME, lnGAME0, lnGAME1)
      
      call Do_Interpolations( &
         1, nvals_EXCOR7, nvals_EXCOR7, nvals_EXCOR7, &
         fq% nlnRS, fq% lnRSs, fq% nlnGAME, fq% lnGAMEs, &
         fq% tbl1, iRS, jGAME, lnRS0, lnRS, lnRS1, lnGAME0, lnGAME, lnGAME1, &
         fval, df_dlnRS, df_dlnGAME, ierr)               
      if (ierr /= 0) then
         write(*,1) 'Do_Interpolations failed in get_EXCOR7'
         return
      end if
      
      FXC = fval(jFXC); dFXC_dlnRS = df_dlnRS(jFXC); dFXC_dlnGAME = df_dlnGAME(jFXC)
      UXC = fval(jUXC); dUXC_dlnRS = df_dlnRS(jUXC); dUXC_dlnGAME = df_dlnGAME(jUXC)
      PXC = fval(jPXC); dPXC_dlnRS = df_dlnRS(jPXC); dPXC_dlnGAME = df_dlnGAME(jPXC)
      CVXC = fval(jCVXC); dCVXC_dlnRS = df_dlnRS(jCVXC); dCVXC_dlnGAME = df_dlnGAME(jCVXC)
      SXC = fval(jSXC); dSXC_dlnRS = df_dlnRS(jSXC); dSXC_dlnGAME = df_dlnGAME(jSXC)
      PDTXC = fval(jPDTXC); dPDTXC_dlnRS = df_dlnRS(jPDTXC); dPDTXC_dlnGAME = df_dlnGAME(jPDTXC)
      PDRXC = fval(jPDRXC); dPDRXC_dlnRS = df_dlnRS(jPDRXC); dPDRXC_dlnGAME = df_dlnGAME(jPDRXC)
   
   end subroutine get_EXCOR7

   
   subroutine Locate_lnRS( &
         lnRS, nlnRS, lnRS_min, dlnRS, iQ, lnRS0, lnRS1)
      real(dp), intent(inout) :: lnRS
      integer, intent(in) :: nlnRS
      real(dp), intent(in) :: lnRS_min, dlnRS
      integer, intent(out) :: iQ
      real(dp), intent(out) :: lnRS0, lnRS1
      iQ = int((lnRS - lnRS_min)/dlnRS + 1d-4) + 1         
      if (iQ < 1 .or. iQ >= nlnRS) then            
         if (iQ < 1) then
            iQ = 1
            lnRS0 = lnRS_min
            lnRS1 = lnRS0 + dlnRS
            lnRS = lnRS0
         else
            iQ = nlnRS-1
            lnRS0 = lnRS_min + (iQ-1)*dlnRS
            lnRS1 = lnRS0 + dlnRS
            lnRS = lnRS1
         end if            
      else         
         lnRS0 = lnRS_min + (iQ-1)*dlnRS
         lnRS1 = lnRS0 + dlnRS
      end if
   end subroutine Locate_lnRS

   
   subroutine Locate_lnGAME( &
         lnGAME, nlnGAME, lnGAME_min, dlnGAME, iQ, lnGAME0, lnGAME1)
      real(dp), intent(inout) :: lnGAME
      integer, intent(in) :: nlnGAME
      real(dp), intent(in) :: lnGAME_min, dlnGAME
      integer, intent(out) :: iQ
      real(dp), intent(out) :: lnGAME0, lnGAME1
      iQ = int((lnGAME - lnGAME_min)/dlnGAME + 1d-4) + 1         
      if (iQ < 1 .or. iQ >= nlnGAME) then            
         if (iQ < 1) then
            iQ = 1
            lnGAME0 = lnGAME_min
            lnGAME1 = lnGAME0 + dlnGAME
            lnGAME = lnGAME0
         else
            iQ = nlnGAME-1
            lnGAME0 = lnGAME_min + (iQ-1)*dlnGAME
            lnGAME1 = lnGAME0 + dlnGAME
            lnGAME = lnGAME1
         end if            
      else         
         lnGAME0 = lnGAME_min + (iQ-1)*dlnGAME
         lnGAME1 = lnGAME0 + dlnGAME
      end if
   end subroutine Locate_lnGAME
   
   
   subroutine Do_Interpolations( &
          nvlo, nvhi, nvals, n, nlnGAMI, x, ny, y, fin1, i, j, &
          x0, xget, x1, y0, yget, y1, &
          fval, df_dx, df_dy, ierr)
      integer, intent(in) :: nvlo, nvhi, nvals, n, nlnGAMI, ny
      real(dp), intent(in) :: x(:) ! (nlnGAMI)
      real(dp), intent(in) :: y(:) ! (ny)
      real(dp), intent(in), pointer :: fin1(:) ! = (4,n,nlnGAMI,ny)
      integer, intent(in) :: i, j           ! target cell in f
      real(dp), intent(in) :: x0, xget, x1      ! x0 <= xget <= x1;  x0 = xs(i), x1 = xs(i+1)
      real(dp), intent(in) :: y0, yget, y1      ! y0 <= yget <= y1;  y0 = ys(j), y1 = ys(j+1)
      real(dp), intent(inout), dimension(nvals) :: fval, df_dx, df_dy
      integer, intent(out) :: ierr

      real(dp), parameter :: z36th = 1d0/36d0
      real(dp) :: xp, xpi, xp2, xpi2, ax, axbar, bx, bxbar, cx, cxi, hx2, cxd, cxdi, hx, hxi
      real(dp) :: yp, ypi, yp2, ypi2, ay, aybar, by, bybar, cy, cyi, hy2, cyd, cydi, hy, hyi
      real(dp) :: sixth_hx2, sixth_hy2, z36th_hx2_hy2
      real(dp) :: sixth_hx, sixth_hxi_hy2, z36th_hx_hy2
      real(dp) :: sixth_hx2_hyi, sixth_hy, z36th_hx2_hy
      integer :: k, ip1, jp1
      real(dp), pointer :: fin(:,:,:,:)
      
      include 'formats'
      
      ierr = 0
      
      fin(1:4,1:n,1:nlnGAMI,1:ny) => &
         fin1(1:4*n*nlnGAMI*ny)
      
      hx=x1-x0
      hxi=1d0/hx
      hx2=hx*hx

      xp=(xget-x0)*hxi

      xpi=1d0-xp
      xp2=xp*xp
      xpi2=xpi*xpi

      ax=xp2*(3d0-2d0*xp)
      axbar=1d0-ax
      
      bx=-xp2*xpi
      bxbar=xpi2*xp

      cx=xp*(xp2-1d0)
      cxi=xpi*(xpi2-1d0)
      cxd=3d0*xp2-1d0
      cxdi=-3d0*xpi2+1d0

      hy=y1-y0
      hyi=1d0/hy
      hy2=hy*hy

      yp=(yget-y0)*hyi
      
      ypi=1d0-yp
      yp2=yp*yp
      ypi2=ypi*ypi

      ay=yp2*(3d0-2d0*yp)
      aybar=1d0-ay
      
      by=-yp2*ypi
      bybar=ypi2*yp

      cy=yp*(yp2-1d0)
      cyi=ypi*(ypi2-1d0)
      cyd=3d0*yp2-1d0
      cydi=-3d0*ypi2+1d0
               
      sixth_hx2 = one_sixth*hx2
      sixth_hy2 = one_sixth*hy2
      z36th_hx2_hy2 = z36th*hx2*hy2
      
      sixth_hx = one_sixth*hx
      sixth_hxi_hy2 = one_sixth*hxi*hy2
      z36th_hx_hy2 = z36th*hx*hy2
      
      sixth_hx2_hyi = one_sixth*hx2*hyi
      sixth_hy = one_sixth*hy
      z36th_hx2_hy = z36th*hx2*hy
      
      ip1 = i+1
      jp1 = j+1
      
      !$omp simd
      do k = nvlo, nvhi
         ! bicubic spline interpolation
         
         ! f(1,i,j) = f(x(i),y(j))
         ! f(2,i,j) = d2f/dx2(x(i),y(j))
         ! f(3,i,j) = d2f/dy2(x(i),y(j))
         ! f(4,i,j) = d4f/dx2dy2(x(i),y(j))

         fval(k) = &
               xpi*( &
                  ypi*fin(1,k,i,j)  +yp*fin(1,k,i,jp1)) &
                  +xp*(ypi*fin(1,k,ip1,j)+yp*fin(1,k,ip1,jp1)) &
               +sixth_hx2*( &
                  cxi*(ypi*fin(2,k,i,j) +yp*fin(2,k,i,jp1))+ &
                  cx*(ypi*fin(2,k,ip1,j)+yp*fin(2,k,ip1,jp1))) &
               +sixth_hy2*( &
                  xpi*(cyi*fin(3,k,i,j) +cy*fin(3,k,i,jp1))+ &
                  xp*(cyi*fin(3,k,ip1,j)+cy*fin(3,k,ip1,jp1))) &
               +z36th_hx2_hy2*( &
                  cxi*(cyi*fin(4,k,i,j) +cy*fin(4,k,i,jp1))+ &
                  cx*(cyi*fin(4,k,ip1,j)+cy*fin(4,k,ip1,jp1)))
         
         ! derivatives of bicubic splines
         df_dx(k) = &
               hxi*( &
                  -(ypi*fin(1,k,i,j)  +yp*fin(1,k,i,jp1)) &
                  +(ypi*fin(1,k,ip1,j)+yp*fin(1,k,ip1,jp1))) &
               +sixth_hx*( &
                  cxdi*(ypi*fin(2,k,i,j) +yp*fin(2,k,i,jp1))+ &
                  cxd*(ypi*fin(2,k,ip1,j)+yp*fin(2,k,ip1,jp1))) &
               +sixth_hxi_hy2*( &
                  -(cyi*fin(3,k,i,j)  +cy*fin(3,k,i,jp1)) &
                  +(cyi*fin(3,k,ip1,j)+cy*fin(3,k,ip1,jp1))) &
               +z36th_hx_hy2*( &
                  cxdi*(cyi*fin(4,k,i,j) +cy*fin(4,k,i,jp1))+ &
                  cxd*(cyi*fin(4,k,ip1,j)+cy*fin(4,k,ip1,jp1)))

         df_dy(k) = &
               hyi*( &
                  xpi*(-fin(1,k,i,j) +fin(1,k,i,jp1))+ &
                  xp*(-fin(1,k,ip1,j)+fin(1,k,ip1,jp1))) &
               +sixth_hx2_hyi*( &
                  cxi*(-fin(2,k,i,j) +fin(2,k,i,jp1))+ &
                  cx*(-fin(2,k,ip1,j)+fin(2,k,ip1,jp1))) &
               +sixth_hy*( &
                  xpi*(cydi*fin(3,k,i,j) +cyd*fin(3,k,i,jp1))+ &
                  xp*(cydi*fin(3,k,ip1,j)+cyd*fin(3,k,ip1,jp1))) &
               +z36th_hx2_hy*( &
                  cxi*(cydi*fin(4,k,i,j) +cyd*fin(4,k,i,jp1))+ &
                  cx*(cydi*fin(4,k,ip1,j)+cyd*fin(4,k,ip1,jp1)))
      
      end do
      
   end subroutine Do_Interpolations
   

   subroutine load_eosPC_support_Info(fq, iZion, ierr) 
      use const_def, only: mesa_data_dir
      use utils_lib, only: alloc_iounit, free_iounit
      use create_EXCOR7_table, only: do_create_EXCOR7_table
      use create_FSCRliq8_table, only: do_create_FSCRliq8_table
      type (eosPC_Support_Info), pointer :: fq
      integer, intent(in) :: iZion ! 0 means EXCOR7
      integer, intent(out) :: ierr
      
      integer :: io_unit, n, j, i, k, iQ, nparams, nvals
      character (len=256) :: filename, fname, cache_fname, temp_cache_fname
      character (len=1000) :: message
      real(dp), pointer :: tbl2_1(:), tbl2(:,:,:)
      real(dp), target :: vec_ary(50)
      real(dp), pointer :: vec(:)
      real(dp) :: lnGAME, lnRS
      
      include 'formats'
      ierr = 0
      vec => vec_ary
      
      if (iZion == 0) then
         filename = 'EXCOR7'
      else if (iZion < 10) then
         write(filename,'(a,i1)') 'FSCRliq8_Zion_0', iZion
      else
         write(filename,'(a,i2)') 'FSCRliq8_Zion_', iZion
      end if
      fname = trim(mesa_data_dir) // '/eosPC_support_data/' // trim(filename) // '.data'
      cache_fname = trim(mesa_data_dir) // '/eosPC_support_data/cache/' // trim(filename) // '.bin'
      temp_cache_fname = trim(mesa_temp_caches_dir) // '/' // trim(filename) // '.bin'
      
      io_unit = alloc_iounit(ierr)
      if (ierr /= 0) then
         write(*,'(a)') 'failed to alloc_iounit'
         call mesa_error(__FILE__,__LINE__)
      end if 
      open(unit=io_unit, FILE=trim(fname), ACTION='READ', STATUS='OLD', IOSTAT=ierr)
      if (ierr /= 0) then ! need to open the table file
         ierr = 0
         if (iZion == 0) then
            call do_create_EXCOR7_table(fname)
         else
            call do_create_FSCRliq8_table(fname,iZion)
         end if
         open(unit=io_unit, FILE=trim(fname), ACTION='READ', STATUS='OLD', IOSTAT=ierr)
         if (ierr /= 0) then 
            write(*,'(a)') 'failed to create ' // trim(fname)
            write(*,'(A)')
            call mesa_error(__FILE__,__LINE__)
         end if
      end if
         
      read(io_unit,*,iostat=ierr) ! skip line 1
      if (ierr /= 0) return
      
      if (iZion == 0) then
         nparams = 8
         nvals = nvals_EXCOR7
      else
         nparams = 9
         nvals = nvals_FSCRliq8
      end if
      read(io_unit,'(a)',iostat=ierr) message ! get parameters line
      if (ierr == 0) call str_to_vector(message, vec, n, ierr)
      if (ierr /= 0 .or. n < nparams) then
         write(*,'(a)') 'failed while reading ' // trim(fname)
         close(io_unit)
         ierr = -1
         return
      end if
      
      read(io_unit,*,iostat=ierr) ! line 3
      if (ierr /= 0) return
      
      fq% Zion = iZion
      fq% nvals = nvals
      fq% nlnRS = int(vec(1))
      fq% lnRS_min = vec(2)*ln10
      fq% lnRS_max = vec(3)*ln10
      fq% dlnRS = vec(4)*ln10
      
      fq% nlnGAME = int(vec(5))
      fq% lnGAME_min = vec(6)*ln10
      fq% lnGAME_max = vec(7)*ln10
      fq% dlnGAME = vec(8)*ln10
      
      if (iZion > 0 .and. int(vec(9)) /= iZion) then
         write(*,*) 'bad value for Zion in file', int(vec(9)), iZion
         close(io_unit)
         ierr = -1
         return
      end if
   
      if (show_allocations) write(*,2) 'allocate pc_support', &
          4*nvals*fq% nlnRS*fq% nlnGAME + fq% nlnRS + fq% nlnGAME
      allocate(fq% tbl1(4*nvals*fq% nlnRS*fq% nlnGAME), &
         fq% lnRSs(fq% nlnRS), fq% lnGAMEs(fq% nlnGAME), STAT=ierr)
      if (ierr /= 0) return      
      fq% tbl(1:4, 1:nvals, 1:fq% nlnRS, 1:fq% nlnGAME) =>  &
           fq% tbl1(1:4*nvals*fq% nlnRS*fq% nlnGAME)
            
      fq% lnRSs(1) = fq% lnRS_min
      do i = 2, fq% nlnRS-1
         fq% lnRSs(i) = fq% lnRSs(i-1) + fq% dlnRS
      end do
      fq% lnRSs(fq% nlnRS) = fq% lnRS_max
      
      fq% lnGAMEs(1) = fq% lnGAME_min
      do i = 2, fq% nlnGAME-1
         fq% lnGAMEs(i) = fq% lnGAMEs(i-1) + fq% dlnGAME
      end do
      fq% lnGAMEs(fq% nlnGAME) = fq% lnGAME_max

      if (use_cache_for_eos) then
         call Read_Cache(fq, cache_fname, ierr)
         if (ierr == 0) then ! got it from the cache
            close(io_unit)
            call free_iounit(io_unit)
            return
         end if
         ierr = 0
      end if
               
      allocate(tbl2_1(nvals*fq% nlnRS*fq% nlnGAME), STAT=ierr)
      if (ierr .ne. 0) return
      
      tbl2(1:nvals, 1:fq% nlnRS, 1:fq% nlnGAME) =>  &
            tbl2_1(1:nvals*fq% nlnRS*fq% nlnGAME)

      do iQ=1,fq% nlnRS
         
         read(io_unit,*,iostat=ierr) 
         if (failed('skip line1')) return
                     
         read(io_unit,'(a)',iostat=ierr) message
         if (ierr == 0) call str_to_double(message, vec(1), ierr)
         if (failed('read lnRS')) return
         lnRS = vec(1)

         read(io_unit,*,iostat=ierr)
         if (failed('skip line2')) return
         
         read(io_unit,*,iostat=ierr)
         if (failed('skip line3')) return
         
         do i=1,fq% nlnGAME
         
            read(io_unit,'(a)',iostat=ierr) message
            if (failed('read line')) then
               write(*,'(a)') trim(message)
               write(*,*) trim(fname)
               write(*,*) 'iQ, i', iQ, i
               write(*,*) 'lnRS', lnRS
               write(*,*) 'bad input line?'
               call mesa_error(__FILE__,__LINE__)
            end if
            
            call str_to_vector(message, vec, n, ierr)
            if (ierr /= 0 .or. n < 1+nvals) then
               write(*,'(a)') trim(message)
               write(*,*) trim(fname)
               write(*,*) 'iQ, i', iQ, i
               write(*,*) 'lnRS', lnRS
               write(*,*) 'bad input line?'
               call mesa_error(__FILE__,__LINE__)
            end if
            lnGAME = vec(1)
            do j=1,nvals
               tbl2(j,iQ,i) = vec(1+j)
            end do
            
         enddo
         
         if(iQ == fq% nlnRS) exit
         read(io_unit,*,iostat=ierr)
         if (failed('skip line4')) return
         read(io_unit,*,iostat=ierr)
         if (failed('skip line5')) return
         
      end do
         
      close(io_unit)
      
      call Make_Interpolation_Data(fq, nvals, tbl2_1, ierr)
      deallocate(tbl2_1)
      if (failed('Make_Interpolation_Data')) return
      
      call Check_Interpolation_Data
      
      if (.not. use_cache_for_eos) then
         call free_iounit(io_unit)
         return
      end if

      open(unit=io_unit, &
         file=trim(switch_str(temp_cache_fname, cache_fname, use_mesa_temp_cache)), &
         iostat=ierr, action='write', form='unformatted')
      if (ierr == 0) then
         write(*,'(a)') ' write ' // trim(cache_fname)
         write(io_unit)  &
            fq% Zion, fq% nlnGAME, fq% lnGAME_min, fq% lnGAME_max, fq% dlnGAME,  &
            fq% nlnRS, fq% lnRS_min, fq% lnRS_max, fq% dlnRS
         write(io_unit) fq% tbl1(1:4*nvals*fq% nlnRS*fq% nlnGAME)
         close(io_unit)
         if (use_mesa_temp_cache) call mv(temp_cache_fname, cache_fname, .true.)
      end if
      call free_iounit(io_unit)
      
      if (iZion == 0) then
         EXCOR7_table_loaded = .true.
      else
         FSCRliq8_Zion_loaded(iZion) = .true.
      end if
      
      contains
      
      subroutine Check_Interpolation_Data
         integer :: i, j, iQ, jtemp
         do i = 1, 4
            do j = 1, nvals
               do iQ = 1, fq% nlnRS
                  do jtemp = 1, fq% nlnGAME
                     if (is_bad(fq% tbl(i,j,iQ,jtemp))) then
                        fq% tbl(i,j,iQ,jtemp) = 0
                     end if
                  end do
               end do
            end do
         end do
      end subroutine Check_Interpolation_Data
      
      logical function failed(str)
         character (len=*), intent(in) :: str
         failed = (ierr /= 0)
         if (failed) write(*,*) 'load_eosPC_support_Info failed: ' // trim(str)
         if (failed) call mesa_error(__FILE__,__LINE__,'load_eosPC_support_Info')
      end function failed
      
   end subroutine load_eosPC_support_Info
      
      
   subroutine Make_Interpolation_Data(fq, nvals, tbl2_1, ierr)
      use interp_2d_lib_db
      type (eosPC_Support_Info), pointer :: fq
      integer, intent(in) :: nvals
      real(dp), pointer :: tbl2_1(:) ! =(nvals, nlnRS, nlnGAME)
      integer, intent(out) :: ierr

      real(dp), allocatable, target :: f1_ary(:) ! data & spline coefficients
      real(dp), pointer :: f1(:), f(:,:,:), tbl2(:,:,:)
      integer :: ibcxmin                   ! bc flag for x=xmin
      real(dp) :: bcxmin(fq% nlnGAME)    ! bc data vs. y at x=xmin
      integer :: ibcxmax                   ! bc flag for x=xmax
      real(dp) :: bcxmax(fq% nlnGAME)     ! bc data vs. y at x=xmax
      integer :: ibcymin                   ! bc flag for y=ymin
      real(dp) :: bcymin(fq% nlnRS)   ! bc data vs. x at y=ymin
      integer :: ibcymax                   ! bc flag for y=ymax
      real(dp) :: bcymax(fq% nlnRS)   ! bc data vs. x at y=ymax
      integer :: ili_lnRSs    ! =1: logRho grid is "nearly" equally spaced
      integer :: ili_lnGAMEs      ! =1: lnGAME grid is "nearly" equally spaced
      integer :: ier            ! =0 on exit if there is no error.
      integer :: iQ, jtemp, ilnGAME, ilnRS
      real(dp) :: fval(nvals), df_dx(nvals), df_dy(nvals)
            
      integer :: v, vlist(3), var, i, j, ii, jj
      character (len=256) :: message
      
      include 'formats'

      ierr = 0

      ! just use "not a knot" bc's at edges of tables
      ibcxmin = 0; bcxmin(:) = 0
      ibcxmax = 0; bcxmax(:) = 0
      ibcymin = 0; bcymin(:) = 0
      ibcymax = 0; bcymax(:) = 0
      
      tbl2(1:nvals, 1:fq% nlnRS, 1:fq% nlnGAME) =>  &
            tbl2_1(1:nvals*fq% nlnRS*fq% nlnGAME)
      ! copy file variables to internal interpolation tables
      do j=1,fq% nlnGAME
         do i=1,fq% nlnRS
            do v = 1, nvals
               fq% tbl(1,v,i,j) = tbl2(v,i,j)
            end do
         end do             
      end do
      
      allocate(f1_ary(4 * fq% nlnRS * fq% nlnGAME))
      
      f1 => f1_ary
      f(1:4, 1:fq% nlnRS, 1:fq% nlnGAME) => &
            f1_ary(1:4*fq% nlnRS*fq% nlnGAME)

      ! create tables for bicubic spline interpolation         
      do v = 1, nvals
         do i=1,fq% nlnRS
            do j=1,fq% nlnGAME
               f(1,i,j) = fq% tbl(1,v,i,j)
            end do
         end do
         call interp_mkbicub_db( &
               fq% lnRSs,fq% nlnRS,fq% lnGAMEs,fq% nlnGAME,f1,fq% nlnRS, &
               ibcxmin,bcxmin,ibcxmax,bcxmax, &
               ibcymin,bcymin,ibcymax,bcymax, &
               ili_lnRSs,ili_lnGAMEs,ier)
         if (ier /= 0) then
            write(*,*) 'interp_mkbicub_db error happened for EXCOR7 value', v
            ierr = 3
            return
         end if
         do i=1,fq% nlnRS
            do j=1,fq% nlnGAME
               fq% tbl(2,v,i,j) = f(2,i,j)
               fq% tbl(3,v,i,j) = f(3,i,j)
               fq% tbl(4,v,i,j) = f(4,i,j)
            end do
         end do
      end do
      
   end subroutine Make_Interpolation_Data
   
   
   subroutine Read_Cache(fq, cache_fname, ierr)
      use utils_lib, only: alloc_iounit, free_iounit
      type (eosPC_Support_Info), pointer :: fq
      character (*), intent(in) :: cache_fname
      integer, intent(out) :: ierr

      real(dp) :: lnGAME_min_in, lnGAME_max_in, dlnGAME_in,  &
            lnRS_min_in, lnRS_max_in, dlnRS_in
      integer :: Zion_in, nlnRS_in, nlnGAME_in, io_unit, i, j
      real(dp), parameter :: tiny = 1d-10
      
      include 'formats'
      
      ierr = 0
      
      io_unit = alloc_iounit(ierr)
      if (ierr /= 0) then
         write(*,'(a)') 'failed to alloc_iounit'
         call mesa_error(__FILE__,__LINE__)
      end if 
      open(unit=io_unit, file=trim(cache_fname), action='read', &
            status='old', iostat=ierr, form='unformatted')
      if (ierr /= 0) then
         call free_iounit(io_unit)
         return
      end if
      
      !write(*,'(a)') 'read ' // trim(cache_fname)
      
      read(io_unit, iostat=ierr)  &
            Zion_in, nlnGAME_in, lnGAME_min_in, lnGAME_max_in, dlnGAME_in,  &
            nlnRS_in, lnRS_min_in, lnRS_max_in, dlnRS_in
      if (ierr /= 0) then
         write(*,*) 'read cache failed'
      end if 

      if (fq% Zion /= Zion_in) then
         ierr = 1
         write(*,*) 'read cache failed for Zion'
      end if 
      if (fq% nlnRS /= nlnRS_in) then
         ierr = 1
         write(*,*) 'read cache failed for nlnRS'
      end if 
      if (fq% nlnGAME /= nlnGAME_in) then
         ierr = 1
         write(*,*) 'read cache failed for nlnGAME'
      end if
      if (abs(fq% lnGAME_min-lnGAME_min_in) > tiny) then
         ierr = 1
         write(*,*) 'read cache failed for eos_lnGAME_min'
      end if    
      if (abs(fq% lnGAME_max-lnGAME_max_in) > tiny) then
         ierr = 1
         write(*,*) 'read cache failed for eos_lnGAME_max'
      end if    
      if (abs(fq% dlnGAME-dlnGAME_in) > tiny) then
         ierr = 1
         write(*,*) 'read cache failed for eos_dlnGAME'
      end if    
      if (abs(fq% lnRS_min-lnRS_min_in) > tiny) then
         ierr = 1
         write(*,*) 'read cache failed for eos_lnRS_min'
      end if    
      if (abs(fq% lnRS_max-lnRS_max_in) > tiny) then
         ierr = 1
         write(*,*) 'read cache failed for eos_lnRS_max'
      end if
      if (abs(fq% dlnRS-dlnRS_in) > tiny) then
         ierr = 1
         write(*,*) 'read cache failed for eos_dlnRS'
      end if
      
      if (ierr /= 0) then
         write(*,*) 'read cache file failed 1'
         call mesa_error(__FILE__,__LINE__,'Read_Cache')
         close(io_unit); return
      end if

      read(io_unit, iostat=ierr) &
         fq% tbl1(1:4*fq% nvals*fq% nlnRS*fq% nlnGAME)
      if (ierr /= 0) then
         write(*,*) 'read cache file failed 2'
         call mesa_error(__FILE__,__LINE__,'Read_Cache')
      end if

      close(io_unit)
      call free_iounit(io_unit)

   end subroutine Read_Cache
   

! ==================  ELECTRON-ION COULOMB LIQUID  =================== *
      subroutine FITION9(GAMI, &
           FION,UION,PION,CVii,PDTii,PDRii)
!                                                       Version 11.09.08
! Non-ideal contributions to thermodynamic functions of classical OCP,
!       corrected at small density for a mixture.
!   Stems from FITION00 v.24.05.00.
! Input: GAMI - ion coupling parameter
! Output: FION - ii free energy / N_i kT
!         UION - ii internal energy / N_i kT
!         PION - ii pressure / n_i kT
!         CVii - ii heat capacity / N_i k
!         PDTii = PION + d(PION)/d ln T = (1/N_i kT)*(d P_{ii}/d ln T)
!         PDRii = PION + d(PION)/d ln\rho
!   Parameters adjusted to Caillol (1999).
      real(dp), intent(in) :: GAMI
      real(dp), intent(out) :: FION, UION, PION, CVii, PDTii, PDRii
      real(dp) :: F0, U0

      real(dp), parameter :: A1=-.907347d0
      real(dp), parameter :: A2=.62849d0
      real(dp) :: A3
      real(dp), parameter :: C1=.004500d0
      real(dp), parameter :: G1=170.0d0
      real(dp), parameter :: C2=-8.4d-5
      real(dp), parameter :: G2=.0037d0
      real(dp), parameter :: SQ32=.8660254038d0 ! SQ32=sqrt(3)/2

      A3=-SQ32-A1/sqrt(A2)
      F0=A1*(sqrt(GAMI*(A2+GAMI)) &
         - A2*log(sqrt(GAMI/A2)+sqrt(1d0+GAMI/A2))) &
         + 2d0*A3*(sqrt(GAMI)-atan(sqrt(GAMI)))
      U0=pow3(sqrt(GAMI))*(A1/sqrt(A2+GAMI)+A3/(1.d0+GAMI))
!   This is the zeroth approximation. Correction:
      UION=U0+C1*GAMI*GAMI/(G1+GAMI)+C2*GAMI*GAMI/(G2+GAMI*GAMI)
      FION=F0+C1*(GAMI-G1*log(1.d0+GAMI/G1)) &
         + C2/2d0*log(1.d0+GAMI*GAMI/G2)
      CVii=-0.5d0*pow3(sqrt(GAMI))*(A1*A2/pow3(sqrt(A2+GAMI)) &
         + A3*(1.d0-GAMI)/pow2(1.d0+GAMI)) &
         - GAMI*GAMI*(C1*G1/pow2(G1+GAMI)+C2*(G2-GAMI*GAMI)/pow2(G2+GAMI*GAMI))
      PION=UION/3.0d0
      PDRii=(4.0d0*UION-CVii)/9.0d0 ! p_{ii} + d p_{ii} / d ln\rho
      PDTii=CVii/3.0d0 ! p_{ii} + d p_{ii} / d ln T
      return
      end subroutine FITION9

   
   end module pc_support
