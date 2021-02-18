program kap_plotter

   use eos_def
   use eos_lib, only: eosDT_get
   use kap_def
   use kap_lib
   use chem_def
   use chem_lib
   use const_lib
   use math_lib
   use num_lib, only : dfridr

   implicit none

   integer, parameter :: species = 14
   integer, parameter :: h1=1, he3=2, he4=3, c12=4, n14=5, o16=6, ne20=7, f20=8, o20=9, &
      mg24=10, na24=11, ne24=12, si28=13, fe56=14
   integer, pointer, dimension(:) :: net_iso, chem_id
   type (EoS_General_Info), pointer :: rq
   real(dp) :: xa(species)

   integer :: eos_handle, kap_handle
   real(dp) :: Rho, T, log10Rho, log10T, X, Z
   real(dp), dimension(num_eos_basic_results) :: res, d_dlnd, d_dlnT
   real(dp), dimension(num_eos_basic_results, species) :: d_dxa
   integer :: ierr
   character (len=32) :: my_mesa_dir

   real(dp) :: Y, z2bar, z53bar, ye, mass_correction, sumx
   real(dp) :: abar, zbar
   real(dp) :: kap, dlnkap_dlnRho, dlnkap_dlnT, frac_Type2
   real(dp) :: kap_fracs(num_kap_fracs), dlnkap_dxa(species)

   real(dp) :: res1

   integer :: nT, nRho, nX, nZ
   real(dp) :: logT_center, delta_logT, logRho_center, delta_logRho
   real(dp) :: logT_min, logT_max, logRho_min, logRho_max

   real(dp) :: X_center, delta_X, Z_center, delta_Z
   real(dp) :: X_min, X_max, Z_min, Z_max

   real(dp) :: logT_step, logRho_step, X_step, Z_step

   integer :: iounit

   integer :: i_var
   integer :: j,k, njs,nks
   real(dp) :: jval, kval

   real(dp) :: var, dvardx_0, dvardx, err, dx_0, xdum
   logical :: doing_partial, doing_dfridr, doing_d_dlnd

   character(len=4) :: xname, yname

   real(dp), parameter :: UNSET = -999
   real(dp), parameter :: min_derivative_error = 1d-4

   namelist /plotter/ &
      nT, nRho, nX, nZ, &
      logT_center, delta_logT, logRho_center, delta_logRho, &
      logT_min, logT_max, logRho_min, logRho_max, &
      X_center, delta_X, Z_center, delta_Z, &
      X_min, X_max, Z_min, Z_max, &
      xname, yname, doing_partial, doing_dfridr, doing_d_dlnd, &
      i_var


   include 'formats'

   ierr = 0

   my_mesa_dir = '../..'
   call const_init(my_mesa_dir,ierr)
   if (ierr /= 0) then
      write(*,*) 'const_init failed'
      stop 1
   end if

   call math_init()

   call chem_init('isotopes.data', ierr)
   if (ierr /= 0) then
      write(*,*) 'failed in chem_init'
      stop 1
   end if

   ! allocate and initialize the eos tables
   call Setup_eos(eos_handle)
   rq => eos_handles(eos_handle)

   ! setup kap
   call Setup_kap(kap_handle)

   allocate(net_iso(num_chem_isos), chem_id(species), stat=ierr)
   if (ierr /= 0) stop 'allocate failed'

   logRho_center = UNSET
   logT_center = UNSET
   X_center = UNSET
   Z_center = UNSET

   delta_logRho = UNSET
   delta_logT = UNSET
   delta_X = UNSET
   delta_Z = UNSET

   doing_dfridr = .false.
   doing_d_dlnd = .true.

   ! get info from namelist
   open(newunit=iounit, file='inlist_plotter')
   read(iounit, nml=plotter)
   close(iounit)

   if (trim(xname) == trim(yname)) then
      write(*,*) 'xname == yname'
      stop
   end if

   ! file for output
   open(newunit=iounit, file='kap_plotter.dat')

   if (doing_partial) then
      if (doing_d_dlnd) then
         write(*,*) 'plotting partial of log(kap) w.r.t. logRho'
         write(iounit,*) 'log(kap) (partial w.r.t. logRho)'
      else
         write(*,*) 'plotting partial of log(kap) w.r.t. logT'
         write(iounit,*) 'log(kap) (partial w.r.t. logT)'
      end if
   else
      write(*,*) 'plotting log(kap)'
      write(iounit,*) 'log(kap)'
   end if

   select case(xname)
   case('T')
      njs = nT
      write(iounit,*) 'log10(T)'
   case('Rho')
      njs = nRho
      write(iounit,*) 'log10(Rho)'
   case('X')
      njs = nX
      write(iounit,*) 'X'
   case('Z')
      njs = nZ
      write(iounit,*) 'Z'
   case default
      write(*,*) 'invalid xname'
      stop
   end select

   select case(yname)
   case('T')
      nks = nT
      write(iounit,*) 'log10(T)'
   case('Rho')
      nks = nRho
      write(iounit,*) 'log10(Rho)'
   case('X')
      nks = nX
      write(iounit,*) 'X'
   case('Z')
      nks = nZ
      write(iounit,*) 'Z'
   case default
      write(*,*) 'invalid yname'
      stop
   end select


   if ((logT_center == UNSET) .or. (delta_logT == UNSET)) then
      logT_center = 0.5d0 * (logT_max + logT_min)
      delta_logT = (logT_max - logT_min)
   else
      logT_min = logT_center - delta_logT * 0.5d0
      logT_max = logT_center - delta_logT * 0.5d0
   end if

   if ((logRho_center == UNSET) .or. (delta_logRho == UNSET)) then
      logRho_center = 0.5d0 * (logRho_max + logRho_min)
      delta_logRho = (logRho_max - logRho_min)
   else
      logRho_min = logRho_center - delta_logRho * 0.5d0
      logRho_max = logRho_center - delta_logRho * 0.5d0
   end if

   if ((X_center == UNSET) .or. (delta_X == UNSET)) then
      X_center = 0.5d0 * (X_max + X_min)
      delta_X = (X_max - X_min)
   else
      X_min = X_center - delta_X * 0.5d0
      X_max = X_center - delta_X * 0.5d0
   end if

   if ((Z_center == UNSET) .or. (delta_Z == UNSET)) then
      Z_center = 0.5d0*(Z_max + Z_min)
      delta_Z = (Z_max - Z_min)
   else
      Z_min = Z_center - delta_Z * 0.5d0
      Z_max = Z_center - delta_Z * 0.5d0
   end if


   if (nT .gt. 1) then
      logT_step = delta_logT / (nT-1d0)
   else
      logT_step = 0
   end if

   if (nRho .gt. 1) then
      logRho_step = delta_logRho / (nRho-1d0)
   else
      logRho_step = 0
   end if

   if (nX .gt. 1) then
      X_step = delta_X / (nX-1d0)
   else
      X_step = 0
   end if

   if (nZ .gt. 1) then
      Z_step = delta_Z / (nZ-1d0)
   else
      Z_step = 0
   end if


   write(iounit,*) nks, njs

   log10T = logT_center
   T = exp10(log10T)
   log10Rho = logRho_center
   Rho = exp10(log10Rho)
   X = X_center
   Z = Z_center

   do j=1,njs !x
      do k=1,nks !y

         select case(xname)
         case('T')
            log10T = logT_min + logT_step*(j - 1)
            T = exp10(log10T)
            jval = log10T
         case('Rho')
            log10Rho = logRho_min + logRho_step*(j - 1)
            rho = exp10(log10Rho)
            jval = log10Rho
         case('X')
            X = X_min + X_step*(j - 1)
            jval = X
         case('Z')
            Z = Z_min + Z_step*(j - 1)
            jval = Z
         end select

         select case(yname)
         case('T')
            log10T = logT_min + logT_step*(k - 1)
            T = exp10(log10T)
            kval = log10T
         case('Rho')
            log10Rho = logRho_min + logRho_step*(k - 1)
            rho = exp10(log10Rho)
            kval = log10Rho
         case('X')
            X = X_min + X_step*(k - 1)
            kval = X
         case('Z')
            Z = Z_min + Z_step*(k - 1)
            kval = Z
         end select

         call Set_Composition

         ! get a set of results for given temperature and density
         call eosDT_get(&
            eos_handle, species, chem_id, net_iso, xa, &
            Rho, log10Rho, T, log10T, &
            res, d_dlnd, d_dlnT, d_dxa, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in eosDT_get'
            write(*,1) 'log10Rho', log10Rho
            write(*,1) 'log10T', log10T
            stop 1
         end if
         
         call kap_get( &
              kap_handle, species, chem_id, net_iso, xa, log10Rho, log10T, &
              res(i_lnfree_e), d_dlnd(i_lnfree_e), d_dlnT(i_lnfree_e), &
              res(i_eta), d_dlnd(i_eta), d_dlnT(i_eta), &
              kap_fracs, kap, dlnkap_dlnRho, dlnkap_dlnT, dlnkap_dxa, ierr)

         frac_Type2 = kap_fracs(i_frac_Type2)

         if (doing_partial) then
            if (doing_d_dlnd) then
               res1 = dlnkap_dlnRho
            else
               res1 = dlnkap_dlnT
            end if
         else
            res1 = log(kap)
         end if

         if (doing_dfridr) then
            var = log(kap)
            if (doing_d_dlnd) then
               dvardx_0 = dlnkap_dlnRho
            else
               dvardx_0 = dlnkap_dlnT
            end if

            dx_0 = 1d-3
            err = 0d0
            dvardx = dfridr(dx_0,dfridr_func,err)
            xdum = (dvardx - dvardx_0)/max(abs(dvardx_0),min_derivative_error)
            res1 = safe_log10(abs(xdum))
         end if


         write(iounit,*) kval, jval, res1
      end do
   end do


   ! deallocate the eos tables
   call Shutdown_eos(eos_handle)

   deallocate(net_iso, chem_id)

   if (ierr /= 0) then
      write(*,*) 'bad result from kap_get'
      stop 1
   end if

contains

   subroutine Setup_eos(handle)
      ! allocate and load the eos tables
      use eos_def
      use eos_lib
      integer, intent(out) :: handle

      character (len=256) :: eos_file_prefix
      integer :: ierr
      logical, parameter :: use_cache = .true.

      eos_file_prefix = 'mesa'

      call eos_init(' ', ' ', ' ', use_cache, ierr)
      if (ierr /= 0) then
         write(*,*) 'eos_init failed in Setup_eos'
         stop 1
      end if

      write(*,*) 'loading eos tables'

      handle = alloc_eos_handle_using_inlist('inlist_plotter', ierr)
      if (ierr /= 0) then
         write(*,*) 'failed trying to allocate eos handle'
         stop 1
      end if

   end subroutine Setup_eos


   subroutine Shutdown_eos(handle)
      use eos_def
      use eos_lib
      integer, intent(in) :: handle
      call free_eos_handle(handle)
      call eos_shutdown
   end subroutine Shutdown_eos


   subroutine Setup_kap(handle)
      ! allocate and load the kap tables
      use kap_def
      use kap_lib
      integer, intent(out) :: handle

      integer :: ierr
      logical, parameter :: use_cache = .true.

      call kap_init(use_cache, ' ', ierr)
      if (ierr /= 0) then
         write(*,*) 'kap_init failed in Setup_eos'
         stop 1
      end if

      write(*,*) 'loading kap tables'

      handle = alloc_kap_handle_using_inlist('inlist_plotter', ierr)
      if (ierr /= 0) then
         write(*,*) 'failed trying to allocate kap handle'
         stop 1
      end if

   end subroutine Setup_kap


   subroutine Shutdown_kap(handle)
      use kap_def
      use kap_lib
      integer, intent(in) :: handle
      call free_kap_handle(handle)
      call kap_shutdown
   end subroutine Shutdown_kap


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

      ! xa = 0
      ! xa(o16) = 0.50d0
      ! xa(ne20) = 0.45d0
      ! xa(mg24) = 0.025d0
      ! xa(ne24) = 0.025d0

      ! xa(j,k) = 0.0000000000000000D+00
      ! xa(j,k) = 2.3070818686211134D-61
      ! xa(j,k) = 1.7933866546042857D-30
      ! xa(j,k) = 2.1886524317518550D-22
      ! xa(j,k) = 1.2133173318677428D-19
      ! xa(j,k) = 5.0000003959509065D-01
      ! xa(j,k) = 2.4380794843721077D-16
      ! xa(j,k) = 0.0000000000000000D+00
      ! xa(j,k) = 4.4999976093022009D-01
      ! xa(j,k) = 2.3466254632930512D-02
      ! xa(j,k) = 1.3104996850960785D-04
      ! xa(j,k) = 2.6402842204789112D-02
      ! xa(j,k) = 5.2668459755217583D-08

      xa = 0d0
      xa(h1) = X
      xa(c12) = 0.5*Z
      xa(o16) = 0.5*Z
      xa(fe56) = 0.0
      xa(he4) = 1d0 - xa(h1) - xa(c12) - xa(o16) - xa(fe56)

      call basic_composition_info( &
         species, chem_id, xa, X, Y, Z, &
         abar, zbar, z2bar, z53bar, ye, mass_correction, sumx)

   end subroutine Set_Composition

   real(dp) function dfridr_func(delta_x) result(val)
      real(dp), intent(in) :: delta_x
      integer :: ierr
      real(dp) :: var, log_var, lnT, lnd
      include 'formats'
      ierr = 0

      lnT = log10T*ln10
      lnd = log10Rho*ln10

      ! must call eos to get new lnfree_e info
      
      if (doing_d_dlnd) then
         log_var = (lnd + delta_x)/ln10
         var = exp10(log_var)
         call eosDT_get( &
            eos_handle, species, chem_id, net_iso, xa, &
            var, log_var, T, log10T, &
            res, d_dlnd, d_dlnT, d_dxa, ierr)
         call kap_get( &
            kap_handle, species, chem_id, net_iso, xa, log_var, log10T, &
            res(i_lnfree_e), d_dlnd(i_lnfree_e), d_dlnT(i_lnfree_e), &
            res(i_eta), d_dlnd(i_eta), d_dlnT(i_eta), &
            kap_fracs, kap, dlnkap_dlnRho, dlnkap_dlnT, dlnkap_dxa, ierr)
      else
         log_var = (lnT + delta_x)/ln10
         var = exp10(log_var)
         call eosDT_get( &
            eos_handle, species, chem_id, net_iso, xa, &
            Rho, log10Rho, var, log_var, &
            res, d_dlnd, d_dlnT, d_dxa, ierr)
         call kap_get( &
            kap_handle, species, chem_id, net_iso, xa, log10Rho, log_var, &
            res(i_lnfree_e), d_dlnd(i_lnfree_e), d_dlnT(i_lnfree_e), &
            res(i_eta), d_dlnd(i_eta), d_dlnT(i_eta), &
            kap_fracs, kap, dlnkap_dlnRho, dlnkap_dlnT, dlnkap_dxa, ierr)
      end if

      val = log(kap)
   end function dfridr_func

end program kap_plotter
