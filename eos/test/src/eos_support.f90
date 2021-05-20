      module eos_support
      use eos_def
      use eos_lib
      use const_def
      use chem_def
      use math_lib
      use utils_lib, only: is_bad_num
      implicit none




      logical, parameter :: use_shared_data_dir = .true.   ! MUST BE .true. FOR RELEASE
      !logical, parameter :: use_shared_data_dir = .false.


      
      real(dp) :: X, Z, Zinit, Y, dXC, dXO, XC, XO, abar, zbar, z2bar, z53bar, ye
      integer, parameter :: species = 7
      integer, parameter :: h1=1, he4=2, c12=3, n14=4, o16=5, ne20=6, mg24=7
      integer, target :: chem_id_array(species)
      integer, pointer, dimension(:) :: chem_id, net_iso
      real(dp) :: xa(species)
      
      
      real(dp), allocatable :: d_dxa(:,:) ! (num_d_dxa_basic_results,species)

      integer :: handle
      type (EoS_General_Info), pointer :: rq

      character (len=eos_name_length) :: eos_names(num_eos_basic_results)
      
         ! if false, then test using data from mesa/eos/data/eos_data
         ! if true, then test using data from mesa/data/eos_data
      
      contains

      

      subroutine Init_Composition(X_in, Zinit_in, XC_in, XO_in)
         use chem_lib
         real(dp), intent(IN) :: X_in, Zinit_in, XC_in, XO_in

         real(dp), parameter :: Zfrac_C = 0.173312d0
         real(dp), parameter :: Zfrac_N = 0.053177d0
         real(dp), parameter :: Zfrac_O = 0.482398d0
         real(dp), parameter :: Zfrac_Ne = 0.098675d0
         
         real(dp) :: Z, frac, dabar_dx(species), dzbar_dx(species), sumx, &
               mass_correction, dmc_dx(species)
         
         chem_id => chem_id_array
         
         allocate(net_iso(num_chem_isos))
         net_iso(:) = 0
         
         chem_id(h1) = ih1; net_iso(ih1) = h1
         chem_id(he4) = ihe4; net_iso(ihe4) = he4
         chem_id(c12) = ic12; net_iso(ic12) = c12
         chem_id(n14) = in14; net_iso(in14) = n14
         chem_id(o16) = io16; net_iso(io16) = o16
         chem_id(ne20) = ine20; net_iso(ine20) = ne20
         chem_id(mg24) = img24; net_iso(img24) = mg24
         
         X = X_in
         Zinit = Zinit_in
         XC = XC_in; XO = XO_in
         Y = 1 - (X + Zinit + XC + XO)
         if (Y < 0) then ! adjust XC and XO
            if (XC + XO <= 0) then
               write(*,*) 'bad args to Init_Composition'
               stop 1
            end if
            frac = (1 - X - Zinit) / (XC + XO)
            if (frac <= 0) stop 'bad args to Init_Composition'
            XC = frac*XC; XO = frac*XO
            Y = 1 - (X+Zinit+XC+XO)
            if (Y < -1d-10) then
               write(*,*) 'screw up in Init_Composition'
               stop 1
            end if
            if (Y < 0) Y = 0
         end if
      
         xa(h1) = X
         xa(he4) = Y
         xa(c12) = Zinit * Zfrac_C + XC
         xa(n14) = Zinit * Zfrac_N
         xa(o16) = Zinit * Zfrac_O + XO
         xa(ne20) = Zinit * Zfrac_Ne
         xa(species) = 1 - sum(xa(1:species-1))
         
         call composition_info( &
               species, chem_id, xa, X, Y, Z, abar, zbar, z2bar, z53bar, &
               ye, mass_correction, sumx, dabar_dx, dzbar_dx, dmc_dx)

      end subroutine Init_Composition


      subroutine Setup_eos
         use chem_lib
         use const_lib
         !..allocate and load the eos tables

         character (len=256) :: my_mesa_dir
         integer :: info
         real(dp) :: logT_all_HELM, logT_all_OPAL
         logical :: use_cache
         
         info = 0
         allocate(d_dxa(num_eos_d_dxa_results,species),stat=info)
         if (info /= 0) then
            write(*,*) 'allocate failed for Setup_eos'
            stop 1
         end if
         
         my_mesa_dir = '../..'
         call const_init(my_mesa_dir,info)
         if (info /= 0) then
            write(*,*) 'const_init failed'
            stop 1
         end if
         
         call math_init()

         call chem_init('isotopes.data', info)
         if (info /= 0) then
            write(*,*) 'chem_init failed'
            stop 1
         end if
         
         use_cache = .true.
         
         call eos_init(' ', use_cache, info)
         if (info /= 0) then
            write(*,*) 'failed in eos_init'
            stop 1
         end if
         eos_names = eosDT_result_names

         handle = alloc_eos_handle_using_inlist('inlist', info)
         if (info /= 0) then
            write(*,*) 'failed in alloc_eos_handle_using_inlist'
            stop 1
         end if

         call eos_ptr(handle,rq,info)
         if (info /= 0) then
            write(*,*) 'failed in eos_ptr'
            stop 1
         end if
      
      end subroutine Setup_eos
      
      
      subroutine Shutdown_eos
         call free_eos_handle(handle)
         call eos_shutdown
         deallocate(d_dxa)
      end subroutine Shutdown_eos

         
      
      end module eos_support

