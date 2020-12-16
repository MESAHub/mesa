      module test_kap_support

      use kap_def
      use kap_lib
      use const_def, only: dp, ln10, arg_not_provided
      use math_lib
      use utils_lib, only: mesa_error
      
      implicit none

      logical, parameter :: use_shared_data_dir = .true. ! if false, then test using local version data
      !logical, parameter :: use_shared_data_dir = .false.
      logical, parameter :: use_cache = .true.
      logical, parameter :: show_info = .false.

      character (len=32) :: my_mesa_dir
      integer, parameter :: ionmax = 7
   
      integer :: nz
      real(dp), pointer, dimension(:) :: &
            x, y, z, c, n, o, ne, lgT, lgd, &
            old_lgKap, old_opacity, old_d_opacity_dlnd, old_d_opacity_dlnT
      real(dp), pointer, dimension(:,:) :: &
            new_lgKap, new_opacity, new_d_opacity_dlnd, new_d_opacity_dlnT
      integer, parameter :: num_new = 2
      
      ! for eos
      
      double precision :: abar, zbar, z2bar, z53bar, ye
      integer, parameter :: species = 7
      integer, parameter :: h1=1, he4=2, c12=3, n14=4, o16=5, ne20=6, mg24=7
      integer, pointer, dimension(:) :: net_iso, chem_id
      double precision :: xa(species)

      real(dp) :: Zbase

      contains
      
      
      subroutine Do_One(quietly)

         use kap_def, only: Kap_General_Info
        
         logical, intent(in) :: quietly
         integer :: ierr

         type (Kap_General_Info), pointer :: rq

         call setup(quietly)

         call test1(quietly, 1, 'fixed metals', ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)

         call test1(quietly, 2, 'C/O enhanced', ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)

         !call test1(quietly, 3, 'op_mono', ierr)
         !if (ierr /= 0) call mesa_error(__FILE__,__LINE__)

         call test1(quietly, 4, 'AESOPUS', ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)

         call kap_shutdown

      end subroutine Do_One


      subroutine test1_op_mono(quietly, test_str)
         logical, intent(in) :: quietly
         character (len=*), intent(in) :: test_str
         real(dp) :: &
            zbar, Z, xh, XC, XN, XO, XNe, frac, abar, kap1, &
            fC, fN, fO, fNe, dXC, dXO, xmass(ionmax), &
            frac_Type2, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            logT, logRho, logR, kap, log10kap, dlnkap_dlnRho, dlnkap_dlnT
            
         logical :: CO_enhanced
         logical, parameter :: dbg = .false.
         integer :: chem_id(ionmax), ierr, fe56
         real(dp) :: chem_factors(ionmax)
         
         include 'formats.dek'
         
         ierr = 0
         
         call setup_op_mono(ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in setup_op_mono'
            return
         end if
         
         lnfree_e=0; d_lnfree_e_dlnRho=0; d_lnfree_e_dlnT=0
         xc = 0d0
         xn = 0d0
         xo = 0d0
         xne = 0d0

            logT = 5.4d0
          logRho = -5.7d0
               z = 0.02d0
              xh = 0.65d0
                          
         
         call get_composition_info(Z, xh, abar, zbar, chem_id, xmass)
            
         chem_factors(:) = 1d0 ! scale factors for element opacity

         frac_Type2 = 0d0
         call test_op_mono(0,ierr)
         if (ierr /= 0) return

         log10kap = safe_log10(kap)
         
         if (.not. quietly) then
            write(*,*) trim(test_str)
            write(*,*)
            call show_args
            call show_results
         end if
         
         ! test element factors with pure Fe56                          
         write(*,1) 'pure fe56; factors all 1.0'
         call get_pure_fe56_composition_info(abar, zbar, chem_id, xmass, fe56)
         call test_op_mono(fe56,ierr)
         if (ierr /= 0) return
         write(*,*)
         kap1 = kap
            
         chem_factors(fe56) = 1.75d0
         write(*,1) 'pure fe56; fe56 factor increased', chem_factors(fe56)
         call test_op_mono(fe56,ierr)
         if (ierr /= 0) return
         write(*,*)
         write(*,1) 'new/old', kap/kap1, kap, kap1
         write(*,*)

      
         contains
         
         
         subroutine setup_op_mono(ierr)
            integer, intent(out) :: ierr
            character (len=256) :: op_mono_data_path, op_mono_data_cache_filename
            
            ierr = 0
            
            call GET_ENVIRONMENT_VARIABLE( &
               "MESA_OP_MONO_DATA_PATH", op_mono_data_path, status=ierr, trim_name=.true.)
            if (ierr /= 0) then
               write(*,*) 'failed to get environment variable MESA_OP_MONO_DATA_PATH'
               return
            end if
            call GET_ENVIRONMENT_VARIABLE( &
               "MESA_OP_MONO_DATA_CACHE_FILENAME", op_mono_data_cache_filename, &
               status=ierr, trim_name=.true.)
            if (ierr /= 0) then
               write(*,*) 'failed to get environment variable MESA_OP_MONO_DATA_CACHE_FILENAME'
               return
            end if

            write(*,*) 'call load_op_mono_data'
            call load_op_mono_data( &
               op_mono_data_path, op_mono_data_cache_filename, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in load_op_mono_data'
               write(*,*) 'my_mesa_dir ' // trim(my_mesa_dir)
               write(*,*) 'op_mono_data_path ' // trim(op_mono_data_path)
               write(*,*) 'op_mono_data_cache_filename ' // trim(op_mono_data_cache_filename)
               return
            end if
         
         end subroutine setup_op_mono
         
         
         subroutine test_op_mono(fe56,ierr)
            use const_def, only: Lsun, Rsun, pi
            integer, intent(in) :: fe56
            integer, intent(out) :: ierr
         
            real, pointer :: &
               umesh(:), semesh(:), ff(:,:,:,:), ta(:,:,:,:), rs(:,:,:)
            integer :: kk, nel, nptot, ipe, nrad, i, iz(ionmax), iZ_rad(ionmax)
            real(dp), dimension(ionmax) :: fap, fac, gp1, &
               lgrad
            real(dp) :: flux, L, r
            logical, parameter :: screening = .true.
            !logical, parameter :: screening = .false.
            
            include 'formats'
         
            ierr = 0

            !write(*,*) 'call get_op_mono_params'
            call get_op_mono_params(nptot, ipe, nrad)
            allocate( &
               umesh(nptot), semesh(nptot), ff(nptot,ipe,4,4), &
               ta(nptot,nrad,4,4), &
               rs(nptot,4,4), stat=ierr)
            if (ierr /= 0) return
            
            !write(*,*) 'call get_op_mono_args'
            call get_op_mono_args( &
               ionmax, xmass, 0d0, chem_id, chem_factors, &
               nel, iz, fap, fac, ierr)
            if (ierr /= 0) then
               write(*,*) 'error in get_op_mono_args, ierr = ',ierr
               return
            end if

            L = 11d0*Lsun
            r = 0.15d0*Rsun
            flux = L/(4*pi*r*r)
            if (fe56 > 0) then
               kk = 1
               iZ_rad(1) = iZ(fe56)
            else
               kk = ionmax
               iZ_rad(:) = iZ(:)
            end if
            
            call op_mono_get_radacc( &
               ! input
               kk, iZ_rad, ionmax, iZ, fap, fac, &
               flux, logT, logRho, screening, &
               ! output
               log10kap, &
               lgrad, &
               ! work arrays
               umesh, semesh, ff, ta, rs, &
               ierr)
               
            deallocate(umesh, semesh, ff, rs, ta)
            if (ierr /= 0) then
               write(*,*) 'error in op_mono_get_radacc, ierr = ',ierr
               return
            end if
            
            kap = exp10(log10kap)
            
            if (fe56 > 0) then
               write(*,*)
               write(*,1) 'grad', exp10(lgrad(kk))
               write(*,1) 'lgrad', lgrad(kk)
            end if
            
            write(*,*)
            write(*,1) 'kap', kap
            write(*,1) 'log10kap', log10kap
            write(*,1) 'dlnkap_dlnRho', dlnkap_dlnRho
            write(*,1) 'dlnkap_dlnT', dlnkap_dlnT
            write(*,1) 'sum fap', sum(fap(:))
            write(*,*)

         end subroutine test_op_mono
         
      
         subroutine show_args
            1 format(a40,1pe26.16)
            write(*,*) 'CO_enhanced', CO_enhanced
            write(*,1) 'logT', logT
            write(*,1) 'logRho', logRho
            write(*,1) 'Z', Z
            write(*,1) 'Zbase', Zbase
            write(*,1) 'zbar', zbar
            write(*,1) 'xh', xh
            write(*,1) 'xc', xc
            write(*,1) 'xn', xn
            write(*,1) 'xo', xo
            write(*,1) 'xne', xne
            write(*,1) 'lnfree_e', lnfree_e
            write(*,*)
         end subroutine show_args
         
      
         subroutine show_results
            use utils_lib
            1 format(a40,1pe26.16)
            write(*,1) 'log10kap', log10kap
            write(*,1) 'dlnkap_dlnRho', dlnkap_dlnRho
            write(*,1) 'dlnkap_dlnT', dlnkap_dlnT
            write(*,*)
            write(*,1) 'kap', kap
            write(*,1) 'dkap_dlnd', dlnkap_dlnRho*kap
            write(*,1) 'dkap_dlnT', dlnkap_dlnT*kap
            write(*,*)
            write(*,1) 'frac_Type2', frac_Type2
            write(*,*)
            if (is_bad(log10kap)) then
               write(*,*) 'bad log10kap'
            end if
         end subroutine show_results
         
         
      end subroutine test1_op_mono
      
   
      subroutine test1(quietly, which, test_str, ierr)
         use kap_def, only: Kap_General_Info
         logical, intent(in) :: quietly
         integer, intent(in) :: which
         character (len=*), intent(in) :: test_str
         integer, intent(out) :: ierr
         real(dp) :: &
            zbar, Z, xh, XC, XN, XO, XNe, frac, abar, xmass(ionmax), &
            fC, fN, fO, fNe, dXC, dXO, &
            frac_Type2, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            logT, logRho, logR, kap, log10kap, dlnkap_dlnRho, dlnkap_dlnT, &
            lnd, lnT, dx_0, dvardx_0, err, dvardx, xdum

         character(len=64) :: inlist
         integer :: handle
         type (Kap_General_Info), pointer :: rq
         
         logical :: CO_enhanced
         logical, parameter :: dbg = .false., test_partial = .true.
         integer :: chem_id(ionmax)
         logical :: doing_d_dlnd
         include 'formats.dek'
         
         ierr = 0

         lnfree_e=0; d_lnfree_e_dlnRho=0; d_lnfree_e_dlnT=0
         Zbase = 0d0
         xc = 0d0
         xn = 0d0
         xo = 0d0
         xne = 0d0
         
         select case(which)
         case (0) ! special test
            
         case (1) ! fixed

            inlist = 'inlist_test_fixed'
            
            CO_enhanced = .false.
            logT =    6d0
            logRho =   -6d0
            Z =    0.018d0
            xh =    0.7d0

         case (2) ! co

            inlist = 'inlist_test_co'
            
            CO_enhanced = .true.
            logT =    6d0
            logRho =   -6d0
            Zbase = 0.018d0
            dXC = 0.021d0
            dXO = 0.019d0
            xh = 0.0d0
            fC = 0.173312d0
            fN = 0.053152d0
            fO = 0.482398d0
            fNe = 0.098668d0
            Z = Zbase + dXC + dXO
            xc = dXC + fC*Zbase
            xn = fN*Zbase
            xo = dXO + fO*Zbase
            xne = fNe*Zbase

         case (3) ! OP

            call mesa_error(__FILE__,__LINE__)
            
         case (4) ! AESOPUS

            inlist = 'inlist_aesopus'
            CO_enhanced = .true.

            ! conditions from RG surface
            logT =    3.5571504546260235D+000
            logRho = -8.2496430699014667D+000
            Z =       2.0074120713487353D-002
            xh =      6.7888662523180188D-001
            xc =      2.9968458709806432D-003
            xn =      1.5270900145630591D-003
            xo =      9.3376263240514384D-003
            xc =      2.9968458709806432D-003
            xn =      1.5270900145630591D-003
            xo =      9.3376263240514384D-003
            
         end select

         handle = alloc_kap_handle_using_inlist(inlist, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)

         call kap_ptr(handle,rq,ierr)

         call get_composition_info(Z, xh, abar, zbar, chem_id, xmass)

         call kap_get( &
              handle, zbar, xh, Z, XC, XN, XO, XNe, logRho, logT, &
              lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
              frac_Type2, kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)

         log10kap = safe_log10(kap)
         
         if (.not. quietly) then
            write(*,*) 'test number', which
            write(*,*) trim(test_str)
            write(*,*)
            call show_args
            call show_results
         end if
         
         
         if (test_partial .and. which == 0) then
         
            !doing_d_dlnd = .true.
            doing_d_dlnd = .false.
         
            lnd = logRho*ln10
            lnT = logT*ln10

            if (doing_d_dlnd) then
               dx_0 = 1d-4*max(1d0, abs(lnd))
               dvardx_0 = dlnkap_dlnRho ! analytic value of partial
            else
               dx_0 = 1d-4*max(1d0, abs(lnT))
               dvardx_0 = dlnkap_dlnT ! analytic value of partial
            end if
            err = 0d0
            dvardx = dfridr(dx_0,err)
            xdum = (dvardx - dvardx_0)/max(abs(dvardx_0),1d-50)
            if (doing_d_dlnd) then
               write(*,1) 'dlnkap_dlnRho analytic, numeric, diff, rel diff', &
                     dvardx_0, dvardx, err, xdum
            else ! doing d_dlnT
               write(*,1) 'dlnkap_dlnT analytic, numeric, diff, rel diff', &
                     dvardx_0, dvardx, err, xdum
            end if
            write(*,*)
         
         end if
         
         
         
      
         contains
         
         
         real(dp) function dfridr_func(delta_x) result(val)
            real(dp), intent(in) :: delta_x
            integer :: ierr
            real(dp) :: log_var
            include 'formats'
            ierr = 0
            
            
            
            ! must call eos to get new lnfree_e info
            
            
            
            
            
            if (doing_d_dlnd) then
               log_var = (lnd + delta_x)/ln10
               call kap_get( &
                    handle, zbar, xh, Z, XC, XN, XO, XNe, log_var, logT, &
                    lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
                    frac_Type2, kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
            else
               log_var = (lnT + delta_x)/ln10
               call kap_get( &
                    handle, zbar, xh, Z, XC, XN, XO, XNe, logRho, log_var, &
                    lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
                    frac_Type2, kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
            end if
            val = log(kap)
         end function dfridr_func

         real(dp) function dfridr(hx,err) ! from Frank
            real(dp), intent(in) :: hx
            real(dp), intent(out) :: err
            !  this routine returns the first derivative of a function func(x)
            !  at the point x, by ridders method of polynomial extrapolation.
            !  value hx is the initial step size;
            !  it should be an increment for which func changes substantially.
            !  an estimate of the error in the first derivative is returned in err.
            integer, parameter :: ntab = 20
            integer :: i,j
            real(dp) :: x,errt,fac,hh,a(ntab,ntab),xdum,ydum
            real(dp), parameter :: con2=2d0, con=sqrt(con2), big=1d50, safe=2d0
            include 'formats'
            dfridr = 0d0
            hh = hx
            ! 2nd order central difference
            a(1,1) = (dfridr_func(hh) - dfridr_func(-hh))/(2d0*hh)
            write(*,2) 'dfdx hh', 1, a(1,1), hh
            err = big
            ! succesive columns in the neville tableu will go to smaller stepsizes
            ! and higher orders of extrapolation
            do i=2,ntab
               hh = hh/con
               a(1,i) = (dfridr_func(hh) - dfridr_func(-hh))/(2d0*hh)
               !write(*,2) 'dfdx hh', i, a(1,i), hh
               ! compute extrapolations of various orders; the error stratagy is to compare
               ! each new extrapolation to one order lower but both at the same stepsize
               ! and at the previous stepsize
               fac = con2
               do j=2,i
                  a(j,i) = (a(j-1,i)*fac - a(j-1,i-1))/(fac-1d0)
                  fac = con2*fac
                  errt = max(abs(a(j,i)-a(j-1,i)),abs(a(j,i)-a(j-1,i-1)))
                  if (errt <= err) then
                     err = errt
                     dfridr = a(j,i)
                     write(*,3) 'dfridr err', i, j, dfridr, err
                  end if
               end do
               ! if higher order is worse by a significant factor safe, then bail
               if (abs(a(i,i) - a(i-1,i-1)) >= safe*err) then
                  write(*,1) 'higher order is worse', err, a(i,i), a(i-1,i-1)
                  return
               end if
            end do
         end function dfridr

         subroutine show_args
            1 format(a40,1pe26.16)
            write(*,*) 'CO_enhanced', CO_enhanced
            write(*,1) 'logT', logT
            write(*,1) 'logRho', logRho
            write(*,1) 'Z', Z
            write(*,1) 'Zbase', rq% Zbase
            write(*,1) 'zbar', zbar
            write(*,1) 'xh', xh
            write(*,1) 'xc', xc
            write(*,1) 'xn', xn
            write(*,1) 'xo', xo
            write(*,1) 'xne', xne
            write(*,1) 'lnfree_e', lnfree_e
            write(*,*)
         end subroutine show_args

         subroutine show_results
            use utils_lib
            1 format(a40,1pe26.16)
            write(*,1) 'log10kap', log10kap
            write(*,1) 'dlnkap_dlnRho', dlnkap_dlnRho
            write(*,1) 'dlnkap_dlnT', dlnkap_dlnT
            write(*,*)
            write(*,1) 'kap', kap
            write(*,1) 'dkap_dlnd', dlnkap_dlnRho*kap
            write(*,1) 'dkap_dlnT', dlnkap_dlnT*kap
            write(*,*)
            write(*,1) 'frac_Type2', frac_Type2
            write(*,*)
            if (is_bad(log10kap)) then
               write(*,*) 'bad log10kap'
            end if
         end subroutine show_results

      end subroutine test1
      

      subroutine setup(quietly)
         use chem_lib
         use const_lib
         logical, intent(in) :: quietly
         
         character (len=256) :: kap_dir, opal_dir, cbeg_ferg
         integer :: ierr
         logical, parameter :: use_cache = .true.
         
         my_mesa_dir = '../..'
         
         call const_init(my_mesa_dir,ierr)
         if (ierr /= 0) then
            write(*,*) 'const_init failed'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         call math_init()

         call chem_init('isotopes.data', ierr)
         if (ierr /= 0) then
            write(*,*) 'chem_init failed'
            call mesa_error(__FILE__,__LINE__)
         end if

         call kap_init(use_cache, '', ierr) 
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         
      end subroutine setup
      
      
      subroutine get_composition_info(Z, X, abar, zbar, chem_id, xmass)
         use chem_def, only: ih1, ihe4, ic12, in14, io16, ine20, img24
         real(dp), intent(in) :: Z, X
         real(dp), intent(out) :: abar, zbar, xmass(ionmax)
         integer, intent(out) :: chem_id(ionmax)

         real(dp) :: Y
         real(dp) :: aion(ionmax),zion(ionmax),ymass(ionmax),zbarxx,ytot1
         integer :: iz(ionmax)

         real(dp), parameter :: Zfrac_C = 0.173312d0
         real(dp), parameter :: Zfrac_N = 0.053177d0
         real(dp), parameter :: Zfrac_O = 0.482398d0
         real(dp), parameter :: Zfrac_Ne = 0.098675d0
         real(dp), parameter :: Zfrac_Mg = 1d0 - (Zfrac_C + Zfrac_N + Zfrac_O + Zfrac_Ne)
         
         integer :: i
         integer :: h1,he4,c12,n14,o16,ne20,mg24

         
         Y = 1 - (X+Z)
         if (Y < 0) then ! adjust XC and XO
            write(*,*) 'bad args to get_composition_info'
            call mesa_error(__FILE__,__LINE__)
         end if
      
         h1        = 1
         iz(h1) = 1
         zion(h1)  = 1.0d0
         aion(h1)  = 1.0d0
         xmass(h1) = X
         chem_id(h1) = ih1

         he4        = 2
         iz(he4) = 2
         zion(he4)  = 2.0d0
         aion(he4)  = 4.0d0
         xmass(he4) = Y
         chem_id(he4) = ihe4

         c12        = 3
         iz(c12) = 6
         zion(c12)  = 6.0d0
         aion(c12)  = 12.0d0
         xmass(c12) = Z * Zfrac_C
         chem_id(c12) = ic12

         n14        = 4
         iz(n14) = 7
         zion(n14)  = 7.0d0
         aion(n14)  = 14.0d0
         xmass(n14) = Z * Zfrac_N
         chem_id(n14) = in14

         o16        = 5
         iz(o16) = 8
         zion(o16)  = 8.0d0
         aion(o16)  = 16.0d0
         xmass(o16) = Z * Zfrac_O
         chem_id(o16) = io16

         ne20       = 6
         iz(ne20) = 10
         zion(ne20)  = 10.0d0
         aion(ne20)  = 20.0d0
         xmass(ne20) = Z * Zfrac_Ne
         chem_id(ne20) = ine20

         mg24       = 7
         iz(mg24) = 12
         zion(mg24)  = 12.0d0
         aion(mg24)  = 24.0d0
         xmass(mg24) = Z * Zfrac_Mg
         chem_id(mg24) = img24

         zbarxx  = 0.0d0
         ytot1   = 0.0d0
         do i=1,ionmax
            ymass(i) = xmass(i)/aion(i)
            ytot1    = ytot1 + ymass(i)
            zbarxx   = zbarxx + zion(i) * ymass(i)
         enddo
         abar   = 1.0d0/ytot1
         zbar   = zbarxx * abar

      end subroutine get_composition_info
      
      
      subroutine get_pure_fe56_composition_info(abar, zbar, chem_id, xmass, fe56)
         use chem_def, only: ih1, ihe4, ic12, in14, io16, ine20, ife56
         real(dp), intent(out) :: abar, zbar, xmass(ionmax)
         integer, intent(out) :: chem_id(ionmax), fe56

         integer :: h1,he4,c12,n14,o16,ne20
      
         h1        = 1
         chem_id(h1) = ih1

         he4        = 2
         chem_id(he4) = ihe4

         c12        = 3
         chem_id(c12) = ic12

         n14        = 4
         chem_id(n14) = in14

         o16        = 5
         chem_id(o16) = io16

         ne20       = 6
         chem_id(ne20) = ine20

         fe56       = 7
         chem_id(fe56) = ife56

         xmass(:) = 0d0
         xmass(fe56) = 1d0

         abar   = 56.0d0
         zbar   = 26.0d0

      end subroutine get_pure_fe56_composition_info
      
      
      subroutine write_plot_data(handle)
         use chem_def
         use utils_lib, only: mkdir
         integer, intent(in) :: handle
         integer, parameter :: io_unit0 = 40
         character (len=256) :: dir
         real(dp), pointer, dimension(:,:,:) :: output_values, co_output_values
         real(dp) :: kap_elect, Z, xh, XC, XN, XO, XNe, abar, zbar, xmass(ionmax), &
            logT_max, logT_min, logRho_max, logRho_min, dlogT, dlogRho, logT, logRho
         integer :: logT_points, logRho_points, num_out, io_params, io_rho, io_tmp, &
            io_first, io_last, i, j, k, ierr, io, chem_id(ionmax)
         logical, parameter :: compare_to_CO = .false.
         
         dir = 'plot_data'
         call mkdir(dir)
         write(*,*) 'write data for opacity plots to ' // trim(dir)

         xh = 0.1d0
         Z = 0.7d0
         
         xh = 0.7d0
         Z = 0.02d0
         
         XC = 3.4416106108119959D-03
         XN = 0.0d0
         XO = 9.3604466091399743D-03
         XNe = 0.0d0

         logT_points = 201
         logRho_points = 201
         
         logT_max = 7d0 ! 5.74d0
         logT_min = 4d0 ! 5.63d0
         logRho_max = 0d0 ! -4.8d0
         logRho_min = -7d0 ! -5.1d0

         call get_composition_info(z, xh, abar, zbar, chem_id, xmass)
         
         kap_elect = 0.2d0*(1 + xh)
         
         io_params = io_unit0
         io_rho = io_unit0+1
         io_tmp = io_unit0+2
         io_first = io_unit0+3
         call Open_Plot_Outfiles(io_first, io_last, io_params, io_rho, io_tmp, dir, compare_to_CO)
         num_out = io_last - io_first + 1
         
         allocate(output_values(logRho_points,logT_points,num_out), &
            co_output_values(logRho_points,logT_points,num_out))

         write(io_params, '(99(f16.6,6x))') Z, xh, XC, XO
         write(io_params, '(99(i16,6x))') logRho_points, logT_points
         close(io_params)


         dlogT = (logT_max - logT_min)/(logT_points-1)
         dlogRho = (logRho_max - logRho_min)/(logRho_points-1)
         
         ierr = 0
         
         do j=1, logT_points
            logT = logT_min + dlogT*(j-1)
            do i=1,logRho_points
               logRho = logRho_min + dlogRho*(i-1)
               if (xc /= 0 .or. xo /= 0) then
                  call do1_CO_plot_data( &
                     handle, logRho, logT, zbar, xh, Z, XC, XN, XO, XNe, &
                     output_values, i, j, compare_to_CO, ierr)
               else
                  call do1_plot_data( &
                     handle, logRho, logT, zbar, xh, Z, &
                     output_values, i, j, compare_to_CO, ierr)
                  if (compare_to_CO) then
                     call do1_CO_plot_data( &
                        handle, logRho, logT, zbar, xh, Z, XC, XN, XO, XNe, &
                        co_output_values, i, j, compare_to_CO, ierr)
                  end if
               end if
               if (ierr /= 0) exit
            end do
         end do

 01   format(e30.22)
         ! write out the results
         do j=1,logT_points
            write(io_tmp,01) logT_min + dlogT*(j-1)
         end do
         close(io_tmp)

         do i=1,logRho_points
            write(io_rho,01) logRho_min + dlogRho*(i-1)
         enddo
         close(io_rho)
         
         if (compare_to_CO) then
            write(*,*) 1
            write(io_first,'(e14.6)') output_values(1:logRho_points,1:logT_points,1)
            write(*,*) 2
            write(io_first+1,'(e14.6)') co_output_values(1:logRho_points,1:logT_points,1)
            write(*,*) 3
            write(io_first+2,'(e14.6)') &
               output_values(1:logRho_points,1:logT_points,1) - &
               co_output_values(1:logRho_points,1:logT_points,1)
         else
            do k = 1, num_out
               write(*,*) k
               write(io_first+k-1,'(e14.6)') output_values(1:logRho_points,1:logT_points,k)
            end do
         end if
      
         do io=io_first,io_last
            close(io)
         end do
         
         deallocate(output_values, co_output_values)
         
      end subroutine write_plot_data
      

      subroutine Open_Plot_Outfiles( &
            io_first, io_last, io_params, io_rho, io_tmp, dir, compare_to_CO)
         integer, intent(IN) :: io_first, io_params, io_rho, io_tmp
         integer, intent(OUT) :: io_last
         character (len=256), intent(IN) :: dir
         logical, intent(in) :: compare_to_CO
         character (len=256) :: fname
         integer :: io
         
         fname = trim(dir) // '/params.data'
         open(unit=io_params,file=trim(fname))
         
         fname = trim(dir) // '/logRho.data'
         open(unit=io_rho,file=trim(fname))
         
         fname = trim(dir) // '/logT.data'
         open(unit=io_tmp,file=trim(fname))
         
         if (compare_to_CO) then
            io = io_first
            fname = trim(dir) // '/kap.data'
            open(unit=io,file=trim(fname))
            fname = trim(dir) // '/kapCO.data'
            io = io+1; open(unit=io,file=trim(fname))
            fname = trim(dir) // '/kap_sub_kapCO.data'
            io = io+1; open(unit=io,file=trim(fname))
            fname = trim(dir) // '/dfridr.data'
            io = io+1; open(unit=io,file=trim(fname))
         else
            io = io_first
            fname = trim(dir) // '/logK.data'
            open(unit=io,file=trim(fname))
            fname = trim(dir) // '/dlogK_dlogRho.data'
            io = io+1; open(unit=io,file=trim(fname))
            fname = trim(dir) // '/dlogK_dlogT.data'
            io = io+1; open(unit=io,file=trim(fname))
            fname = trim(dir) // '/logKec.data'
            io = io+1; open(unit=io,file=trim(fname))
            fname = trim(dir) // '/dlogKec_dlogRho.data'
            io = io+1; open(unit=io,file=trim(fname))
            fname = trim(dir) // '/dlogKec_dlogT.data'
            io = io+1; open(unit=io,file=trim(fname))
            fname = trim(dir) // '/dfridr.data'
            io = io+1; open(unit=io,file=trim(fname))
         end if
            
         io_last = io
      
      end subroutine Open_Plot_Outfiles
      
      
      subroutine do1_plot_data( &
            handle, lgd, lgT, zbar, xh, Z, &
            output_values, i, j, compare_to_CO, ierr)
         use utils_lib, only: is_bad
         real(dp), intent(in) :: lgd, lgT, zbar
         integer, intent(in) :: handle, i, j
         real(dp), intent(in) :: xh, Z
         real(dp), intent(out) :: output_values(:,:,:)
         logical, intent(in) :: compare_to_CO
         integer, intent(out) :: ierr

         real(dp) :: frac_Type2, XC, XN, XO, XNe
         
         real(dp) :: kap, dlnkap_dlnRho, dlnkap_dlnT, &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            lnd, lnT, dx_0, dvardx_0, err, xdum, dvardx
         logical :: doing_d_dlnd, kap_ec
         include 'formats'
         
         ierr = 0
         lnfree_e=0; d_lnfree_e_dlnRho=0; d_lnfree_e_dlnT=0
         
         call kap_get( &
              handle, zbar, xh, Z, XC, XN, XO, XNe, lgd, lgT, &
              lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
              frac_Type2, kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         
         if (ierr /= 0) then
            output_values(i,j,:) = -99
            ierr = 0
            return
         end if
         
         if (is_bad(kap)) then
            write(*,*) 'kap', kap
            stop 'do1_plot_data'
         end if

         if (compare_to_CO) then
            output_values(i,j,1) = kap
         else
            output_values(i,j,1) = safe_log10(kap)
         end if
         output_values(i,j,2) = dlnkap_dlnRho
         output_values(i,j,3) = dlnkap_dlnT
         
         if (.false. .and. i == 1) then
            if (j == 1) write(*,*) 'logRho, X, Z', lgd, xh, Z
            if (j == 1) write(*,*) 'logT logKap'
            write(*,*) lgT, output_values(i,j,1)
         end if
         
         !kap_ec = .true.
         kap_ec = .false.
         
         if (.not. kap_ec) then
         
            doing_d_dlnd = .true.
            !doing_d_dlnd = .false.
         
            lnd = lgd*ln10
            lnT = lgT*ln10

            if (doing_d_dlnd) then
               dx_0 = max(1d-4, abs(lnd*1d-4))
               dvardx_0 = dlnkap_dlnRho ! analytic value of partial
            else
               dx_0 = max(1d-4, abs(lnT*1d-4))
               dvardx_0 = dlnkap_dlnT ! analytic value of partial
            end if
            err = 0d0
            dvardx = dfridr(dx_0,err)
            xdum = (dvardx - dvardx_0)/max(abs(dvardx_0),1d-5)
            output_values(i,j,7) = safe_log10(abs(xdum))
         
         end if

         call kap_get_elect_cond_opacity( &
                  zbar, lgd, lgT, &
                  kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         if (ierr /= 0) then
            output_values(i,j,4) = -99d0
            output_values(i,j,5) = -99d0
            output_values(i,j,6) = -99d0
            output_values(i,j,7) = -99d0
            ierr = 0
            return
         end if
         
         if (is_bad(kap)) then
            write(*,*) 'kap_get_elect_cond_opacity kap', kap
            stop 'do1_plot_data'
         end if

         output_values(i,j,4) = safe_log10(kap)
         output_values(i,j,5) = dlnkap_dlnRho
         output_values(i,j,6) = dlnkap_dlnT
            
         if (kap_ec) then
         
            !doing_d_dlnd = .true.
            doing_d_dlnd = .false.
         
            lnd = lgd*ln10
            lnT = lgT*ln10

            if (doing_d_dlnd) then
               dx_0 = max(1d-4, abs(lnd*1d-4))
               dvardx_0 = dlnkap_dlnRho ! analytic value of partial
            else
               dx_0 = max(1d-4, abs(lnT*1d-4))
               dvardx_0 = dlnkap_dlnT ! analytic value of partial
            end if
            err = 0d0
            dvardx = dfridr(dx_0,err)
            xdum = (dvardx - dvardx_0)/max(abs(dvardx_0),1d-6)
            output_values(i,j,7) = safe_log10(abs(xdum))
         
         end if
         
            
         contains
         
         
         real(dp) function dfridr_func(delta_x) result(val)
            real(dp), intent(in) :: delta_x
            integer :: ierr
            real(dp) :: log_var
            include 'formats'
            ierr = 0
            if (doing_d_dlnd) then
               log_var = (lnd + delta_x)/ln10
               if (kap_ec) then
                  call kap_get_elect_cond_opacity( &
                     zbar, log_var, lgT, &
                     kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
               else
                  call kap_get( &
                       handle, zbar, xh, Z, XC, XN, XO, XNe, log_var, lgT, &
                       lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
                       frac_Type2, kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
               end if
            else
               log_var = (lnT + delta_x)/ln10
               if (kap_ec) then
                  call kap_get_elect_cond_opacity( &
                     zbar, lgd, log_var, &
                     kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
               else
                  call kap_get( &
                       handle, zbar, xh, Z, XC, XN, XO, XNe, lgd, log_var, &
                       lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
                       frac_Type2, kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
               end if
            end if
            val = log(kap)
         end function dfridr_func

         real(dp) function dfridr(hx,err) ! from Frank
            real(dp), intent(in) :: hx
            real(dp), intent(out) :: err
            !  this routine returns the first derivative of a function func(x)
            !  at the point x, by ridders method of polynomial extrapolation.
            !  value hx is the initial step size;
            !  it should be an increment for which func changes substantially.
            !  an estimate of the error in the first derivative is returned in err.
            integer, parameter :: ntab = 20
            integer :: i,j
            real(dp) :: x,errt,fac,hh,a(ntab,ntab),xdum,ydum
            real(dp), parameter :: con2=2d0, con=sqrt(con2), big=1d50, safe=2d0
            include 'formats'
            dfridr = 0d0
            hh = hx
            ! 2nd order central difference
            a(1,1) = (dfridr_func(hh) - dfridr_func(-hh))/(2d0*hh)
            !write(*,2) 'dfdx hh', 1, a(1,1), hh
            err = big
            ! succesive columns in the neville tableu will go to smaller stepsizes
            ! and higher orders of extrapolation
            do i=2,ntab
               hh = hh/con
               a(1,i) = (dfridr_func(hh) - dfridr_func(-hh))/(2d0*hh)
               !write(*,2) 'dfdx hh', i, a(1,i), hh
               ! compute extrapolations of various orders; the error stratagy is to compare
               ! each new extrapolation to one order lower but both at the same stepsize
               ! and at the previous stepsize
               fac = con2
               do j=2,i
                  a(j,i) = (a(j-1,i)*fac - a(j-1,i-1))/(fac-1d0)
                  fac = con2*fac
                  errt = max(abs(a(j,i)-a(j-1,i)),abs(a(j,i)-a(j-1,i-1)))
                 ! write(*,3) 'a(j,i)', j, i, a(j,i), errt
                  if (errt <= err) then
                     err = errt
                     dfridr = a(j,i)
                     !write(*,3) 'dfridr err', i, j, dfridr, err
                  end if
               end do
               ! if higher order is worse by a significant factor safe, then bail
               if (abs(a(i,i) - a(i-1,i-1)) >= safe*err) then
                  !write(*,1) 'higher order is worse', err, a(i,i), a(i-1,i-1)
                  return
               end if
            end do
         end function dfridr

      end subroutine do1_plot_data
         
      
      subroutine do1_CO_plot_data( &
            handle, lgd, lgT, zbar, xh, Z, XC, XN, XO, XNe, &
            output_values, i, j, compare_to_CO, ierr)
         use utils_lib, only: is_bad
         real(dp), intent(in) :: lgd, lgT, zbar
         integer, intent(in) :: handle, i, j
         real(dp), intent(in) :: xh, Z, XC, XN, XO, XNe
         real(dp), intent(out) :: output_values(:,:,:)
         logical, intent(in) :: compare_to_CO
         integer, intent(out) :: ierr
         real(dp) :: frac_Type2, kap, dlnkap_dlnRho, dlnkap_dlnT, &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            lnd, lnT, dx_0, dvardx_0, err, xdum, dvardx
         logical :: doing_d_dlnd
         ierr = 0
         lnfree_e=0; d_lnfree_e_dlnRho=0; d_lnfree_e_dlnT=0
         call kap_get( &
              handle, zbar, xh, Z, XC, XN, XO, XNe, lgd, lgT, &
              lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
              frac_Type2, kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         if (ierr /= 0) then
            output_values(i,j,1) = -99d0
            output_values(i,j,2) = -99d0
            output_values(i,j,3) = -99d0
            output_values(i,j,4) = -99d0
            output_values(i,j,5) = -99d0
            output_values(i,j,6) = -99d0
            output_values(i,j,7) = -99d0
            ierr = 0
            return
         end if
         
         if (is_bad(kap)) then
            write(*,*) 'kap', kap, lgd, lgT, zbar, xh, Z, xc, xo
            stop 'do1_CO_plot_data'
         end if
         
         if (is_bad(dlnkap_dlnRho)) then
            write(*,*) 'dlnkap_dlnRho', dlnkap_dlnRho, lgd, lgT, zbar, xh, Z, xc, xo
            stop 'do1_CO_plot_data'
         end if
         
         if (is_bad(dlnkap_dlnT)) then
            write(*,*) 'dlnkap_dlnT', dlnkap_dlnT, lgd, lgT, zbar, xh, Z, xc, xo
            stop 'do1_CO_plot_data'
         end if

         if (compare_to_CO) then
            output_values(i,j,1) = kap
         else
            output_values(i,j,1) = safe_log10(kap)
         end if
         output_values(i,j,2) = dlnkap_dlnRho
         output_values(i,j,3) = dlnkap_dlnT

         output_values(i,j,4) = 0d0
         output_values(i,j,5) = 0d0
         output_values(i,j,6) = 0d0

         doing_d_dlnd = .true.
         !doing_d_dlnd = .false.
         
         lnd = lgd*ln10
         lnT = lgT*ln10

         if (doing_d_dlnd) then
               dx_0 = max(1d-4, abs(lnd*1d-4))
            dvardx_0 = dlnkap_dlnRho ! analytic value of partial
         else
            dx_0 = max(1d-4, abs(lnT*1d-4))
            dvardx_0 = dlnkap_dlnT ! analytic value of partial
         end if
         err = 0d0
         dvardx = dfridr(dx_0,err)
         xdum = (dvardx - dvardx_0)/max(abs(dvardx_0),1d-50)
         output_values(i,j,7) = safe_log10(abs(xdum))
            
         contains
         
         
         real(dp) function dfridr_func(delta_x) result(val)
            real(dp), intent(in) :: delta_x
            integer :: ierr
            real(dp) :: log_var
            include 'formats'
            ierr = 0
            if (doing_d_dlnd) then
               log_var = (lnd + delta_x)/ln10
               call kap_get( &
                  handle, zbar, xh, Z, XC, XN, XO, XNe, log_var, lgT, &
                  lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
                  frac_Type2, kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
            else
               log_var = (lnT + delta_x)/ln10
               call kap_get( &
                  handle, zbar, xh, Z, XC, XN, XO, XNe, lgd, log_var, &
                  lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
                  frac_Type2, kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
            end if
            val = log(kap)
         end function dfridr_func

         real(dp) function dfridr(hx,err) ! from Frank
            real(dp), intent(in) :: hx
            real(dp), intent(out) :: err
            !  this routine returns the first derivative of a function func(x)
            !  at the point x, by ridders method of polynomial extrapolation.
            !  value hx is the initial step size;
            !  it should be an increment for which func changes substantially.
            !  an estimate of the error in the first derivative is returned in err.
            integer, parameter :: ntab = 20
            integer :: i,j
            real(dp) :: x,errt,fac,hh,a(ntab,ntab),xdum,ydum
            real(dp), parameter :: con2=2d0, con=sqrt(con2), big=1d50, safe=2d0
            include 'formats'
            dfridr = 0d0
            hh = hx
            ! 2nd order central difference
            a(1,1) = (dfridr_func(hh) - dfridr_func(-hh))/(2d0*hh)
            !write(*,2) 'dfdx hh', 1, a(1,1), hh
            err = big
            ! succesive columns in the neville tableu will go to smaller stepsizes
            ! and higher orders of extrapolation
            do i=2,ntab
               hh = hh/con
               a(1,i) = (dfridr_func(hh) - dfridr_func(-hh))/(2d0*hh)
               !write(*,2) 'dfdx hh', i, a(1,i), hh
               ! compute extrapolations of various orders; the error stratagy is to compare
               ! each new extrapolation to one order lower but both at the same stepsize
               ! and at the previous stepsize
               fac = con2
               do j=2,i
                  a(j,i) = (a(j-1,i)*fac - a(j-1,i-1))/(fac-1d0)
                  fac = con2*fac
                  errt = max(abs(a(j,i)-a(j-1,i)),abs(a(j,i)-a(j-1,i-1)))
                 ! write(*,3) 'a(j,i)', j, i, a(j,i), errt
                  if (errt <= err) then
                     err = errt
                     dfridr = a(j,i)
                     !write(*,3) 'dfridr err', i, j, dfridr, err
                  end if
               end do
               ! if higher order is worse by a significant factor safe, then bail
               if (abs(a(i,i) - a(i-1,i-1)) >= safe*err) then
                  !write(*,1) 'higher order is worse', err, a(i,i), a(i-1,i-1)
                  return
               end if
            end do
         end function dfridr
         
      end subroutine do1_CO_plot_data
      
      
      subroutine write_logP_plot_data(handle)
         use chem_def
         use utils_lib, only: mkdir
         use eos_lib
         use eos_def
         integer, intent(in) :: handle
         integer, parameter :: io_unit0 = 40
         character (len=256) :: dir
         real(dp), pointer, dimension(:,:,:) :: output_values
         real(dp) :: kap_elect, Z, xh, XC, XN, XO, XNe, abar, zbar, xmass(ionmax), &
            logT_max, logT_min, logP_max, logP_min, dlogT, dlogP, T, logT, logP, &
            logRho_tol, logP_tol, logRho_guess, logPgas, Pgas, Prad, P, rho, &
            logRho, logRho_result, logRho_bnd1, logRho_bnd2, logP_at_bnd1, logP_at_bnd2
         integer :: logT_points, logP_points, num_out, io_params, io_P, io_tmp, max_iter, &
            io_first, io_last, i, j, k, ierr, io, kap_chem_id(ionmax), eos_handle, eos_calls
         real(dp), dimension(:), allocatable :: res, d_dlnd, d_dlnT, &
            d_dlnRho_const_T, d_dlnT_const_Rho, &
            d_dabar_const_TRho, d_dzbar_const_TRho
         
         dir = 'plot_logP_data'
         call mkdir(dir)
         write(*,*) 'write data for logP opacity plots to ' // trim(dir)
         
         allocate( &
            res(num_eos_basic_results), d_dlnd(num_eos_basic_results), d_dlnT(num_eos_basic_results), &
            d_dlnRho_const_T(num_eos_basic_results), d_dlnT_const_Rho(num_eos_basic_results), &
            d_dabar_const_TRho(num_eos_basic_results), d_dzbar_const_TRho(num_eos_basic_results))

         xh = 0.68882d0
         Z = 0.02207d0
         
         call Setup_eos(eos_handle)
         write(*,*) 'done Setup_eos'
         
         call Init_Composition_for_eos(xh, Z)
         write(*,*) 'done Init_Composition_for_eos'
         
         XC = 0d0
         XN = 0.0d0
         XO = 0d0
         XNe = 0.0d0

         logT_points = 101
         logP_points = 101
         
         logT_max = 5.8d0
         logT_min = 5d0
         logP_max = 8d0
         logP_min = 6d0

         call get_composition_info(z, xh, abar, zbar, kap_chem_id, xmass)
         write(*,*) 'done get_composition_info'
         
         kap_elect = 0.2d0*(1 + xh)
         
         io_params = io_unit0
         io_P = io_unit0+1
         io_tmp = io_unit0+2
         io_first = io_unit0+3
         call Open_logP_Plot_Outfiles(io_first, io_last, &
            io_params, io_P, io_tmp, dir)
         num_out = io_last - io_first + 1
         write(*,*) 'done Open_logP_Plot_Outfiles'
         
         allocate(output_values(logP_points,logT_points,num_out))

         write(io_params, '(99(f16.6,6x))') Z, xh, XC, XO
         write(io_params, '(99(i16,6x))') logP_points, logT_points
         close(io_params)


         dlogT = (logT_max - logT_min)/(logT_points-1)
         dlogP = (logP_max - logP_min)/(logP_points-1)
         
         ierr = 0
         
         do j=1, logT_points
            logT = logT_min + dlogT*(j-1)
            write(*,*) 'j', j, logT_points, logT
            T = exp10(logT)
            Prad = Radiation_Pressure(T)
            do i=1,logP_points
               logP = logP_min + dlogP*(i-1)
               write(*,*) 'i', i, logP_points, logP
               P = exp10(logP)
               Pgas = P - Prad
               logPgas = log10(Pgas)
               
               logRho_guess = -6d0
               max_iter = 100
               logRho_bnd1 = arg_not_provided
               logP_at_bnd1 = arg_not_provided
               logRho_bnd2 = arg_not_provided
               logP_at_bnd2 = arg_not_provided

               logRho_tol = 5d-10
               logP_tol = 5d-13

               call eos_gamma_PT_get( &
                  eos_handle, abar, P, logP, T, logT, 5d0/3d0, &
                  rho, logRho_guess, res, d_dlnd, d_dlnT, &
                  ierr)
               if (ierr /= 0) then
                  write(*,*) 'ierr from eos_gamma_PT_get'
                  exit
               end if

               call eosDT_get_Rho_given_Ptotal( &
                  eos_handle, Z, xh, abar, zbar, &
                  species, chem_id, net_iso, xa, &
                  logT, logP, logRho_tol, logP_tol, max_iter, logRho_guess, &
                  logRho_bnd1, logRho_bnd2, logP_at_bnd1, logP_at_bnd2, &
                  logRho_result, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
                  d_dabar_const_TRho, d_dzbar_const_TRho, eos_calls, ierr)
               
               logRho = logRho_result
               write(*,*) 'logT logP logRho', logT, logP, logRho, logRho_guess

               if (ierr /= 0) then
                  write(*,*) 'ierr from eosDT_get_Rho_given_Ptotal'
                  exit
               end if
               
               call do1_plot_data( &
                  handle, logRho, logT, zbar, xh, Z, &
                  output_values, i, j, .false., ierr)
               if (ierr /= 0) then
                  write(*,*) 'ierr from do1_plot_data', logRho, logT, zbar, xh, Z
                  exit
               end if
            end do
         end do

 01   format(e30.22)
         ! write out the results
         do j=1,logT_points
            write(io_tmp,01) logT_min + dlogT*(j-1)
         end do
         close(io_tmp)

         do i=1,logP_points
            write(io_P,01) logP_min + dlogP*(i-1)
         enddo
         close(io_P)
         
         do k = 1, num_out
            write(*,*) k
            write(io_first+k-1,'(e14.6)') output_values(1:logP_points,1:logT_points,k)
         end do
      
         do io=io_first,io_last
            close(io)
         end do
         
         deallocate(output_values)
         
      end subroutine write_logP_plot_data
      

      subroutine Open_logP_Plot_Outfiles( &
            io_first, io_last, io_params, io_P, io_tmp, dir)
         integer, intent(IN) :: io_first, io_params, io_P, io_tmp
         integer, intent(OUT) :: io_last
         character (len=256), intent(IN) :: dir
         character (len=256) :: fname
         integer :: io
         
         fname = trim(dir) // '/params.data'
         open(unit=io_params,file=trim(fname))
         
         fname = trim(dir) // '/logP.data'
         open(unit=io_P,file=trim(fname))
         
         fname = trim(dir) // '/logT.data'
         open(unit=io_tmp,file=trim(fname))
         
         io = io_first
         fname = trim(dir) // '/logK.data'
         open(unit=io,file=trim(fname))
         fname = trim(dir) // '/dlogK_dlogRho.data'
         io = io+1; open(unit=io,file=trim(fname))
         fname = trim(dir) // '/dlogK_dlogT.data'
         io = io+1; open(unit=io,file=trim(fname))
         fname = trim(dir) // '/logKec.data'
         io = io+1; open(unit=io,file=trim(fname))
         fname = trim(dir) // '/dlogKec_dlogRho.data'
         io = io+1; open(unit=io,file=trim(fname))
         fname = trim(dir) // '/dlogKec_dlogT.data'
         io = io+1; open(unit=io,file=trim(fname))
         fname = trim(dir) // '/dfridr.data'
         io = io+1; open(unit=io,file=trim(fname))
            
         io_last = io
      
      end subroutine Open_logP_Plot_Outfiles
      

      subroutine Setup_eos(handle)
         ! allocate and load the eos tables
         use eos_def
         use eos_lib
         integer, intent(out) :: handle

         integer :: ierr
         logical, parameter :: use_cache = .true.

         call eos_init(' ', ' ', ' ', use_cache, ierr)
         if (ierr /= 0) then
            write(*,*) 'eos_init failed in Setup_eos'
            stop 1
         end if
         
         write(*,*) 'loading eos tables'
         
         handle = alloc_eos_handle(ierr)
         if (ierr /= 0) then
            write(*,*) 'failed trying to allocate eos handle'
            stop 1
         end if
      
      end subroutine Setup_eos
      
      
      subroutine Init_Composition_for_eos(X, Z)
         use chem_lib
         use chem_def
         real(dp), intent(inout) :: X, Z

         double precision, parameter :: Zfrac_C = 0.173312d0
         double precision, parameter :: Zfrac_N = 0.053177d0
         double precision, parameter :: Zfrac_O = 0.482398d0
         double precision, parameter :: Zfrac_Ne = 0.098675d0
         
         double precision :: xz, frac, dabar_dx(species), dzbar_dx(species), sumx, &
               mass_correction, dmc_dx(species), Y
         
         allocate(net_iso(num_chem_isos), chem_id(species))

         net_iso(:) = 0
         
         chem_id(h1) = ih1; net_iso(ih1) = h1
         chem_id(he4) = ihe4; net_iso(ihe4) = he4
         chem_id(c12) = ic12; net_iso(ic12) = c12
         chem_id(n14) = in14; net_iso(in14) = n14
         chem_id(o16) = io16; net_iso(io16) = o16
         chem_id(ne20) = ine20; net_iso(ine20) = ne20
         chem_id(mg24) = img24; net_iso(img24) = mg24
         
         Y = 1 - (X + Z)
               
         xa(h1) = X
         xa(he4) = Y
         xa(c12) = Z * Zfrac_C
         xa(n14) = Z * Zfrac_N
         xa(o16) = Z * Zfrac_O
         xa(ne20) = Z * Zfrac_Ne
         xa(species) = 1 - sum(xa(1:species-1))
         
         call composition_info( &
               species, chem_id, xa, X, Y, xz, abar, zbar, z2bar, z53bar, &
               ye, mass_correction, sumx, dabar_dx, dzbar_dx, dmc_dx)
         ! ! for now, we use the approx versions
         ! abar = approx_abar
         ! zbar = approx_zbar

      end subroutine Init_Composition_for_eos
      
      
      end module test_kap_support

