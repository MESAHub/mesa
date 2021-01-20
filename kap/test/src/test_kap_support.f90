      module test_kap_support

      use chem_def
      use chem_lib
      use eos_def
      use eos_lib
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
      integer, parameter :: ionmax = 8
   
      real(dp) :: abar, zbar, z2bar, z53bar, ye, mass_correction, sumx
      integer, parameter :: species = 8
      integer, parameter :: h1=1, he4=2, c12=3, n14=4, o16=5, ne20=6, mg24=7, fe56=8
      integer, pointer, dimension(:) :: net_iso, chem_id
      real(dp) :: xa(species)

      real(dp) :: X, Y, Z, Zbase

      contains
      
      
      subroutine Do_One(quietly)

         logical, intent(in) :: quietly
         integer :: ierr

         call setup(quietly)

         allocate(net_iso(num_chem_isos), chem_id(species))

         net_iso(:) = 0

         chem_id(h1) = ih1; net_iso(ih1) = h1
         chem_id(he4) = ihe4; net_iso(ihe4) = he4
         chem_id(c12) = ic12; net_iso(ic12) = c12
         chem_id(n14) = in14; net_iso(in14) = n14
         chem_id(o16) = io16; net_iso(io16) = o16
         chem_id(ne20) = ine20; net_iso(ine20) = ne20
         chem_id(mg24) = img24; net_iso(img24) = mg24
         chem_id(fe56) = ife56; net_iso(ife56) = fe56

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
         integer :: ierr
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
                          
         
         !call get_composition_info(Z, xh, abar, zbar, chem_id, xmass)
            
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
         !call get_pure_fe56_composition_info(abar, zbar, chem_id, xmass, fe56)
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
         use kap_def, only: Kap_General_Info, num_kap_fracs, i_frac_Type2
         logical, intent(in) :: quietly
         integer, intent(in) :: which
         character (len=*), intent(in) :: test_str
         integer, intent(out) :: ierr
         real(dp) :: &
            XC, XN, XO, XNe, &
            fC, fN, fO, fNe, dXC, dXO, &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            eta, d_eta_dlnRho, d_eta_dlnT, &
            logT, logRho, logR, kap, log10kap, dlnkap_dlnRho, dlnkap_dlnT
         real(dp) :: kap_fracs(num_kap_fracs), dlnkap_dxa(species)

         ! eos results
         real(dp), dimension(num_eos_basic_results) :: res, deos_dlnd, deos_dlnT
         real(dp), dimension(num_eos_d_dxa_results,species) :: deos_dxa

         character(len=64) :: inlist
         integer :: eos_handle, kap_handle
         type (Kap_General_Info), pointer :: rq
         
         logical :: CO_enhanced
         logical, parameter :: dbg = .false.
         include 'formats.dek'
         
         ierr = 0

         xa = 0
         X = 0; Z = 0; xc = 0; xn = 0; xo = 0; xne = 0
         
         select case(which)
         case (0) ! special test
            
         case (1) ! fixed

            inlist = 'inlist_test_fixed'
            
            CO_enhanced = .false.
            logT =    6d0
            logRho =   -6d0
            X = 0.7d0
            Z = 0.018d0
            
            xa(h1) = X
            xa(he4) = 1d0 - X - Z
            xa(fe56) = Z
            
         case (2) ! co

            inlist = 'inlist_test_co'
            
            CO_enhanced = .true.
            logT =    6d0
            logRho =   -6d0
            Zbase = 0.018d0
            dXC = 0.021d0
            dXO = 0.019d0
            X = 0.0d0
            fC = 0.173312d0
            fN = 0.053152d0
            fO = 0.482398d0
            fNe = 0.098668d0
            Z = Zbase + dXC + dXO
            xc = dXC + fC*Zbase
            xn = fN*Zbase
            xo = dXO + fO*Zbase
            xne = fNe*Zbase

            xa(h1) = X
            xa(he4) = 1d0 - X - Z
            xa(c12) = xc
            xa(n14) = xn
            xa(o16) = xo
            xa(ne20) = xne

         case (3) ! OP

            call mesa_error(__FILE__,__LINE__)
            
         case (4) ! AESOPUS

            inlist = 'inlist_aesopus'
            CO_enhanced = .true.

            ! conditions from RG surface
            logT =    3.5571504546260235D+000
            logRho = -8.2496430699014667D+000
            Z =       2.0074120713487353D-002
            X =       6.7888662523180188D-001
            xc =      2.9968458709806432D-003
            xn =      1.5270900145630591D-003
            xo =      9.3376263240514384D-003

            xa(h1) = X
            xa(he4) = 1d0 - X - Z
            xa(c12) = xc
            xa(n14) = xn
            xa(o16) = xo
            
         end select


         call basic_composition_info( &
            species, chem_id, xa, X, Y, Z, abar, zbar, z2bar, z53bar, &
            ye, mass_correction, sumx)

         eos_handle = alloc_eos_handle_using_inlist(inlist, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)

         call eosDT_get( &
               eos_handle, species, chem_id, net_iso, xa, &
               exp10(logRho), logRho, exp10(logT), logT, &
               res, deos_dlnd, deos_dlnT, deos_dxa, ierr)

         lnfree_e = res(i_lnfree_e)
         d_lnfree_e_dlnRho = deos_dlnd(i_lnfree_e)
         d_lnfree_e_dlnT = deos_dlnT(i_lnfree_e)

         eta = res(i_eta)
         d_eta_dlnRho = deos_dlnd(i_eta)
         d_eta_dlnT = deos_dlnT(i_eta)
         
         kap_handle = alloc_kap_handle_using_inlist(inlist, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)

         call kap_ptr(kap_handle,rq,ierr)

         call kap_get( &
              kap_handle, species, chem_id, net_iso, xa, logRho, logT, &
              lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
              eta, d_eta_dlnRho, d_eta_dlnT, &
              kap_fracs, kap, dlnkap_dlnRho, dlnkap_dlnT, dlnkap_dxa, ierr)

         log10kap = safe_log10(kap)
         
         if (.not. quietly) then
            write(*,*) 'test number', which
            write(*,*) trim(test_str)
            write(*,*)
            call show_args
            call show_results
         end if

         contains
         
         subroutine show_args
            1 format(a40,1pe26.16)
            write(*,*) 'CO_enhanced', CO_enhanced
            write(*,1) 'logT', logT
            write(*,1) 'logRho', logRho
            write(*,1) 'Z', Z
            write(*,1) 'Zbase', rq% Zbase
            write(*,1) 'zbar', zbar
            write(*,1) 'X', X
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
            write(*,1) 'frac_Type2', kap_fracs(i_frac_Type2)
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

         call eos_init('', '', '', use_cache, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         
         call kap_init(use_cache, '', ierr) 
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         
      end subroutine setup
      
      
      end module test_kap_support

