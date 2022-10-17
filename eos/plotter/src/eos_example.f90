program eos_plotter
   
   use eos_def
   use eos_lib, only : eos_ptr, eosDT_get, eosDT_get_T_given_Ptotal
   use chem_def
   use chem_lib
   use const_lib
   use math_lib
   use num_lib, only : dfridr
   use utils_lib, only : set_nan
   
   implicit none
   
   integer, parameter :: species = 14
   integer, parameter :: h1 = 1, he3 = 2, he4 = 3, c12 = 4, n14 = 5, o16 = 6, ne20 = 7, f20 = 8, o20 = 9, &
      mg24 = 10, na24 = 11, ne24 = 12, si28 = 13, fe56 = 14
   integer, pointer, dimension(:) :: net_iso, chem_id
   type (EoS_General_Info), pointer :: rq
   real(dp) :: xa(species)
   
   integer :: handle
   real(dp) :: Rho, T, log10Rho, log10T, X, Z
   real(dp), dimension(num_eos_basic_results) :: res, d_dlnd, d_dlnT
   real(dp), dimension(num_eos_d_dxa_results, species) :: d_dxa
   real(dp), dimension(num_eos_basic_results) :: res_other, d_dlnd_other, d_dlnT_other
   real(dp), dimension(num_eos_d_dxa_results, species) :: d_dxa_other
   integer :: ierr
   character (len = 32) :: my_mesa_dir
   
   real(dp) :: p, res1, res2
   
   integer :: nT, nRho, nX, nZ
   real(dp) :: logT_center, delta_logT, logRho_center, delta_logRho
   real(dp) :: logT_min, logT_max, logRho_min, logRho_max
   
   real(dp) :: X_center, delta_X, Z_center, delta_Z
   real(dp) :: X_min, X_max, Z_min, Z_max
   
   real(dp) :: logT_step, logRho_step, X_step, Z_step
   
   integer :: iounit
   
   integer :: i_var, i_max, i_eos, i_eos_other, i_cons, i
   integer :: j, k, njs, nks, eos_calls
   real(dp) :: jval, kval, logT_tol, logP_tol
   
   real(dp) :: var, dvardx_0, dvardx, err, dx_0, xdum, logT_guess, log10P
   logical :: doing_partial, doing_dfridr, doing_d_dlnd, doing_consistency, ignore_ierr
   logical :: only_blend_regions
   real(dp) :: abar, zbar, z2bar, z53bar, ye, mass_correction, sumx, xh, xhe
   integer, parameter :: MAX_ITER_FOR_SOLVE = 100
   
   character(len = 4) :: xname, yname
   
   real(dp), parameter :: UNSET = -999
   real(dp), parameter :: min_derivative_error = 1d-4
   
   namelist /plotter/ &
      nT, nRho, nX, nZ, &
      logT_center, delta_logT, logRho_center, delta_logRho, &
      logT_min, logT_max, logRho_min, logRho_max, &
      X_center, delta_X, Z_center, delta_Z, &
      X_min, X_max, Z_min, Z_max, &
      xname, yname, doing_partial, doing_dfridr, doing_d_dlnd, doing_consistency, &
      i_var, i_eos, i_eos_other, i_cons, ignore_ierr, only_blend_regions
   
   include 'formats'
   
   ierr = 0
   
   my_mesa_dir = '../..'
   call const_init(my_mesa_dir, ierr)
   if (ierr /= 0 .and. .not. ignore_ierr) then
      write(*, *) 'const_init failed'
      call mesa_error(__FILE__, __LINE__)
   end if
   
   call math_init()
   
   call chem_init('isotopes.data', ierr)
   if (ierr /= 0 .and. .not. ignore_ierr) then
      write(*, *) 'failed in chem_init'
      call mesa_error(__FILE__, __LINE__)
   end if
   
   allocate(net_iso(num_chem_isos), chem_id(species), stat = ierr)
   if (ierr /= 0 .and. .not. ignore_ierr) call mesa_error(__FILE__, __LINE__, 'allocate failed')
   
   ! allocate and initialize the eos tables
   call Setup_eos(handle)
   
   log10Rho = -8.0880854137850644D+00
   log10T = 3.9439201289513104D+00
   logT_guess = 3.9439201289513104D+00
   log10P = 1.0608529872259041D+02
   logT_tol = 9.9999999999999994D-12
   logP_tol = 9.9999999999999994D-12
   
   call Set_Composition
   call basic_composition_info(&
      species, chem_id, xa, xh, xhe, z, &
      abar, zbar, z2bar, z53bar, ye, mass_correction, sumx)
   X = xh
   
   ! get a set of results for given temperature and density
   call eosDT_get_T_given_Ptotal(handle, Z, X, abar, zbar, &
      species, chem_id, net_iso, xa, &
      log10Rho, log10P, logT_tol, logP_tol, MAX_ITER_FOR_SOLVE, logT_guess, &
      arg_not_provided, arg_not_provided, arg_not_provided, arg_not_provided, &
      log10T, res, d_dlnd, d_dlnT, &
      d_dlnd_other, d_dlnd_other, eos_calls, ierr)
   
   write(*, *) log10T
   write(*, *) eos_calls, ierr
   
   ! deallocate the eos tables
   call Shutdown_eos(handle)
   
   deallocate(net_iso, chem_id)
   
   if (ierr /= 0 .and. .not. ignore_ierr) then
      write(*, *) 'bad result from eos_get'
      call mesa_error(__FILE__, __LINE__)
   end if
   
   stop

contains
   
   subroutine Setup_eos(handle)
      ! allocate and load the eos tables
      use eos_def
      use eos_lib
      integer, intent(out) :: handle
      
      character (len = 256) :: eos_file_prefix
      integer :: ierr
      logical, parameter :: use_cache = .true.
      
      eos_file_prefix = 'mesa'
      
      call eos_init(' ', use_cache, ierr)
      if (ierr /= 0 .and. .not. ignore_ierr) then
         write(*, *) 'eos_init failed in Setup_eos'
         call mesa_error(__FILE__, __LINE__)
      end if
      
      write(*, *) 'loading eos tables'
      
      handle = alloc_eos_handle_using_inlist('inlist_plotter', ierr)
      if (ierr /= 0 .and. .not. ignore_ierr) then
         write(*, *) 'failed trying to allocate eos handle'
         call mesa_error(__FILE__, __LINE__)
      end if
   
   end subroutine Setup_eos
   
   
   subroutine Shutdown_eos(handle)
      use eos_def
      use eos_lib
      integer, intent(in) :: handle
      call free_eos_handle(handle)
      call eos_shutdown
   end subroutine Shutdown_eos
   
   
   subroutine Set_Composition
      
      net_iso(:) = 0
      
      chem_id(h1) = ih1; net_iso(ih1) = h1
      chem_id(he3) = ihe3; net_iso(ihe3) = he3
      chem_id(he4) = ihe4; net_iso(ihe4) = he4
      chem_id(c12) = ic12; net_iso(ic12) = c12
      chem_id(n14) = in14; net_iso(in14) = n14
      chem_id(o16) = io16; net_iso(io16) = o16
      chem_id(ne20) = ine20; net_iso(ine20) = ne20
      chem_id(f20) = if20; net_iso(if20) = f20
      chem_id(o20) = io20; net_iso(io20) = o20
      chem_id(mg24) = img24; net_iso(img24) = mg24
      chem_id(na24) = ina24; net_iso(ina24) = na24
      chem_id(ne24) = ine24; net_iso(ine24) = ne24
      chem_id(si28) = isi28; net_iso(isi28) = si28
      chem_id(fe56) = ife56; net_iso(ife56) = fe56
      
      xa = 0d0
      xa(h1) = 7.3571032624883048D-01
      xa(he3) = 4.7616221251102759D-05
      xa(he4) = 2.6121382354237399D-01
      xa(c12) = 3.6157725062562829D-04
      xa(n14) = 3.5477305122973770D-04
      xa(o16) = 1.2885109738934157D-03
      xa(ne20) = 2.8004524103287246D-04
      xa(mg24) = 7.4332747076278781D-04
      xa(si28) = 1.0256615760513371D-99
   
   end subroutine Set_Composition

end program eos_plotter
