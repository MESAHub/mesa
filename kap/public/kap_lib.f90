! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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


      module kap_lib
      ! library for calculating opacities

      ! the data interface for the library is defined in kap_def
      
      use const_def, only: dp
      use math_lib
      
      implicit none


      contains ! the procedure interface for the library
      ! client programs should only call these routines.
            
      
      ! call this routine to initialize the kap module. 
      ! only needs to be done once at start of run.
      ! Reads data from the 'kap' directory in the data_dir.
      ! If use_cache is true and there is a 'kap/cache' directory, it will try that first.
      ! If it doesn't find what it needs in the cache, 
      ! it reads the data and writes the cache for next time.    
      subroutine kap_init(use_cache, kap_cache_dir, ierr)
         use kap_def, only : kap_def_init, kap_use_cache, kap_is_initialized
         logical, intent(in) :: use_cache
         character (len=*), intent(in) :: kap_cache_dir ! blank means use default
         integer, intent(out) :: ierr ! 0 means AOK.         
         ierr = 0
         if (kap_is_initialized) return
         call kap_def_init(kap_cache_dir)
         kap_use_cache = use_cache
         kap_is_initialized = .true.
      end subroutine kap_init

      
      subroutine kap_shutdown
         use kap_def, only: do_Free_Kap_Tables, kap_is_initialized
         call do_Free_Kap_Tables()
         kap_is_initialized = .false.
      end subroutine kap_shutdown

      
      ! after kap_init has finished, you can allocate a "handle".
      
      integer function alloc_kap_handle(ierr) result(handle)
         integer, intent(out) :: ierr ! 0 means AOK.
         character (len=0) :: inlist
         handle = alloc_kap_handle_using_inlist(inlist, ierr)
      end function alloc_kap_handle

      integer function alloc_kap_handle_using_inlist(inlist,ierr) result(handle)
         use kap_def, only:do_alloc_kap,kap_is_initialized
         use kap_ctrls_io, only:read_namelist
         character (len=*), intent(in) :: inlist ! empty means just use defaults.
         integer, intent(out) :: ierr ! 0 means AOK.
         ierr = 0
         if (.not. kap_is_initialized) then
            ierr=-1
            return
         endif
         handle = do_alloc_kap(ierr)
         if (ierr /= 0) return
         call read_namelist(handle, inlist, ierr)
         if (ierr /= 0) return
         call kap_setup_tables(handle, ierr)
         call kap_setup_hooks(handle, ierr)
      end function alloc_kap_handle_using_inlist

      subroutine free_kap_handle(handle)
         ! frees the handle and all associated data
         use kap_def,only: Kap_General_Info,do_free_kap
         integer, intent(in) :: handle
         call do_free_kap(handle)
      end subroutine free_kap_handle
      
      
      subroutine kap_ptr(handle,rq,ierr)
         use kap_def,only:Kap_General_Info,get_kap_ptr,kap_is_initialized
         integer, intent(in) :: handle ! from alloc_kap_handle
         type (Kap_General_Info), pointer :: rq
         integer, intent(out):: ierr
         if (.not. kap_is_initialized) then
            ierr=-1
            return
         endif
         call get_kap_ptr(handle,rq,ierr)
      end subroutine kap_ptr


      subroutine kap_setup_tables(handle, ierr)
         use kap_def, only : Kap_General_Info, get_kap_ptr
         use load_kap, only : Setup_Kap_Tables
         integer, intent(in) :: handle
         integer, intent(out):: ierr
         
         type (Kap_General_Info), pointer :: rq
         logical, parameter :: use_cache = .true.
         logical, parameter :: load_on_demand = .true.

         ierr = 0
         call get_kap_ptr(handle,rq,ierr)
         call Setup_Kap_Tables(rq, use_cache, load_on_demand, ierr)

      end subroutine kap_setup_tables


      subroutine kap_setup_hooks(handle, ierr)
         use kap_def, only : Kap_General_Info, get_kap_ptr
         use other_elect_cond_opacity
         use other_compton_opacity
         integer, intent(in) :: handle
         integer, intent(out):: ierr

         type (Kap_General_Info), pointer :: rq

         ierr = 0
         call get_kap_ptr(handle,rq,ierr)

         rq% other_elect_cond_opacity => null_other_elect_cond_opacity
         rq% other_compton_opacity => null_other_compton_opacity

      end subroutine kap_setup_hooks


      ! kap evaluation
      ! you can call these routines after you've setup the tables for the handle.
      ! NOTE: the structures referenced via the handle are read-only
      ! for the evaulation routines, so you can do multiple evaluations in parallel
      ! using the same handle. 
      
      
      subroutine kap_get( &
         handle, species, chem_id, net_iso, xa, &
         logRho, logT, &
         lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
         eta, d_eta_dlnRho, d_eta_dlnT , &
         kap_fracs, kap, dlnkap_dlnRho, dlnkap_dlnT, dlnkap_dxa, ierr)
         use chem_def, only: chem_isos
         use chem_lib, only: basic_composition_info
         use kap_def, only : kap_is_initialized, Kap_General_Info, num_kap_fracs, i_frac_Type2
         use kap_eval, only : Get_kap_Results
         ! INPUT
         integer, intent(in) :: handle ! from alloc_kap_handle; in star, pass s% kap_handle
         integer, intent(in) :: species
         integer, pointer :: chem_id(:) ! maps species to chem id
         integer, pointer :: net_iso(:) ! maps chem id to species number
         real(dp), intent(in) :: xa(:) ! mass fractions
         real(dp), intent(in) :: logRho ! density
         real(dp), intent(in) :: logT ! temperature
         real(dp), intent(in) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
         ! free_e := total combined number per nucleon of free electrons and positrons
         real(dp), intent(in)  :: eta, d_eta_dlnRho, d_eta_dlnT
         ! eta := electron degeneracy parameter from eos

         ! OUTPUT
         real(dp), intent(out) :: kap_fracs(num_kap_fracs)
         real(dp), intent(out) :: kap ! opacity
         real(dp), intent(out) :: dlnkap_dlnRho ! partial derivative at constant T
         real(dp), intent(out) :: dlnkap_dlnT   ! partial derivative at constant Rho
         real(dp), intent(out) :: dlnkap_dxa(:) ! partial derivative w.r.t. species
         integer, intent(out) :: ierr ! 0 means AOK.

         type (Kap_General_Info), pointer :: rq

         real(dp) :: X, Y, Z, abar, zbar, z2bar, z53bar, ye, mass_correction, sumx
         real(dp) :: XC, XN, XO, XNe
         integer :: i, iz

         ierr = 0
         if (.not. kap_is_initialized) then
            ierr=-1
            return
         endif

         call kap_ptr(handle,rq,ierr)
         if (ierr /= 0) return

         call basic_composition_info( &
            species, chem_id, xa, X, Y, Z, &
            abar, zbar, z2bar, z53bar, ye, mass_correction, sumx)

         xc = 0; xn = 0; xo = 0; xne = 0
         do i=1, species
            iz = chem_isos% Z(chem_id(i))
            select case(iz)
            case (6)
               xc = xc + xa(i)
            case (7)
               xn = xn + xa(i)
            case (8)
               xo = xo + xa(i)
            case (10)
               xne = xne + xa(i)
            end select
         end do

         call Get_kap_Results( &
            rq, zbar, X, Z, XC, XN, XO, XNe, logRho, logT, &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            eta, d_eta_dlnRho, d_eta_dlnT, &
            kap_fracs, kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)

         ! composition derivatives not implemented
         dlnkap_dxa = 0

      end subroutine kap_get


      subroutine kap_get_elect_cond_opacity( &
            handle, zbar, logRho, logT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         use condint, only: do_electron_conduction
         use kap_def, only : kap_is_initialized, Kap_General_Info
         integer, intent(in) :: handle ! from alloc_kap_handle; in star, pass s% kap_handle
         real(dp), intent(in) :: zbar ! average ionic charge (for electron conduction)
         real(dp), intent(in) :: logRho ! the density
         real(dp), intent(in) :: logT ! the temperature
         real(dp), intent(out) :: kap ! electron conduction opacity
         real(dp), intent(out) :: dlnkap_dlnRho ! partial derivative at constant T
         real(dp), intent(out) :: dlnkap_dlnT   ! partial derivative at constant Rho
         integer, intent(out) :: ierr ! 0 means AOK.

         type (Kap_General_Info), pointer :: rq

         if (.not. kap_is_initialized) then
            ierr=-1
            return
         endif
         ierr = 0

         call kap_ptr(handle,rq,ierr)
         if (ierr /= 0) return

         call do_electron_conduction( &
            rq, zbar, logRho, logT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)

      end subroutine kap_get_elect_cond_opacity
      
      
      subroutine kap_get_compton_opacity( &
         handle, &
         Rho, T, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
         eta, d_eta_dlnRho, d_eta_dlnT, &
         kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         use kap_eval, only: Compton_Opacity
         use kap_def, only : kap_is_initialized, Kap_General_Info
         integer, intent(in) :: handle ! kap handle; from star, pass s% kap_handle
         real(dp), intent(in) :: Rho, T
         real(dp), intent(in) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
            ! free_e := total combined number per nucleon of free electrons and positrons
         real(dp), intent(in)  :: eta, d_eta_dlnT, d_eta_dlnRho
            ! eta := electron degeneracy parameter from eos
         real(dp), intent(out) :: kap ! electron conduction opacity
         real(dp), intent(out) :: dlnkap_dlnRho, dlnkap_dlnT
         integer, intent(out) :: ierr ! 0 means AOK.

         type (Kap_General_Info), pointer :: rq

         if (.not. kap_is_initialized) then
            ierr=-1
            return
         endif
         ierr = 0

         call kap_ptr(handle,rq,ierr)
         if (ierr /= 0) return

         call Compton_Opacity(rq, &
            Rho, T, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            eta, d_eta_dlnRho, d_eta_dlnT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)

      end subroutine kap_get_compton_opacity


      subroutine kap_get_radiative_opacity( &
         handle, &
         X, Z, XC, XN, XO, XNe, logRho, logT, &
         frac_lowT, frac_highT, frac_Type2, kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)

         use kap_eval, only: Get_kap_Results_blend_T
         use kap_def, only : kap_is_initialized, Kap_General_Info

         ! INPUT
         integer, intent(in) :: handle ! kap handle; from star, pass s% kap_handle
         real(dp), intent(in) :: X, Z, XC, XN, XO, XNe ! composition
         real(dp), intent(in) :: logRho ! density
         real(dp), intent(in) :: logT ! temperature

         ! OUTPUT
         real(dp), intent(out) :: frac_lowT, frac_highT, frac_Type2
         real(dp), intent(out) :: kap ! opacity
         real(dp), intent(out) :: dlnkap_dlnRho ! partial derivative at constant T
         real(dp), intent(out) :: dlnkap_dlnT   ! partial derivative at constant Rho
         integer, intent(out) :: ierr ! 0 means AOK.

         type (Kap_General_Info), pointer :: rq

         if (.not. kap_is_initialized) then
            ierr=-1
            return
         endif
         ierr = 0

         call kap_ptr(handle,rq,ierr)
         if (ierr /= 0) return

         call Get_kap_Results_blend_T( &
            rq, X, Z, XC, XN, XO, XNe, logRho, logT, &
            frac_lowT, frac_highT, frac_Type2, kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)

      end subroutine kap_get_radiative_opacity

      
      subroutine kap_get_op_mono( &
            handle, zbar, logRho, logT, &
            ! args for op_mono
            use_op_mono_alt_get_kap, &
            nel, izzp, fap, fac, screening, umesh, semesh, ff, rs, &
            ! output
            kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         use kap_def
         use kap_eval, only: combine_rad_with_conduction
         integer, intent(in) :: handle ! from alloc_kap_handle; in star, pass s% kap_handle
         real(dp), intent(in) :: zbar ! average ionic charge (for electron conduction)
         real(dp), intent(in) :: logRho ! the density
         real(dp), intent(in) :: logT ! the temperature
         ! args for op_mono_get_kap
         logical, intent(in) :: use_op_mono_alt_get_kap
         integer, intent(in) :: nel
         integer, intent(in) :: izzp(:) ! (nel)
         real(dp), intent(in) :: fap(:) ! (nel) number fractions of elements
         real(dp), intent(in) :: fac(:) ! (nel) scale factors for element opacity
         logical, intent(in) :: screening
         ! work arrays
         real, pointer :: umesh(:), semesh(:), ff(:,:,:,:), rs(:,:,:)
            ! umesh(nptot)
            ! semesh(nptot)
            ! ff(nptot, ipe, 4, 4)
            ! rs(nptot, 4, 4)
         ! output
         real(dp), intent(out) :: kap ! opacity
         real(dp), intent(out) :: dlnkap_dlnRho ! partial derivative at constant T
         real(dp), intent(out) :: dlnkap_dlnT   ! partial derivative at constant Rho
         integer, intent(out) :: ierr ! 0 means AOK.

         real(dp) :: g1, kap_rad, dlnkap_rad_dlnRho, dlnkap_rad_dlnT
         real(dp) :: Rho, T
         type (Kap_General_Info), pointer :: rq

         ierr = 0
         if (.not. kap_is_initialized) then
            ierr=-1
            return
         endif
         call kap_ptr(handle,rq,ierr)
         if (ierr /= 0) return

         Rho = exp10(logRho)
         T = exp10(logT)
         
         if (use_op_mono_alt_get_kap) then
            call op_mono_alt_get_kap( &
               nel, izzp, fap, fac, logT, logRho, screening, &
               g1, dlnkap_rad_dlnT, dlnkap_rad_dlnRho, &
               umesh, semesh, ff, rs, ierr)
         else
            call op_mono_get_kap( &
               nel, izzp, fap, fac, logT, logRho, screening, &
               g1, dlnkap_rad_dlnT, dlnkap_rad_dlnRho, &
               umesh, semesh, ff, rs, ierr)
         end if
         if (ierr /= 0) return
         
         kap_rad = exp10(g1)

         call combine_rad_with_conduction( &
            rq, rho, logRho, T, logT, zbar, &
            kap_rad, dlnkap_rad_dlnRho, dlnkap_rad_dlnT, &
            kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
         
      end subroutine kap_get_op_mono
      
      
      
      ! interface to OP routines as modified by Haili Hu for radiative levitation in diffusion
      
      ! ref: Hu et al MNRAS 418 (2011)
      
      subroutine load_op_mono_data(op_mono_data_path, op_mono_data_cache_filename, ierr)
         use kap_def
         use op_load, only: op_dload
         character (len=*), intent(in) :: op_mono_data_path, op_mono_data_cache_filename
         integer, intent(out) :: ierr
         if (.not. kap_is_initialized) then
            ierr=-1
            return
         endif
         call op_dload(op_mono_data_path, op_mono_data_cache_filename, ierr)
      end subroutine load_op_mono_data
      
      
      ! sizes for work arrays
      subroutine get_op_mono_params(op_nptot, op_ipe, op_nrad)
         use kap_def
         use op_def, only: nptot, ipe, nrad
         integer, intent(out) :: op_nptot, op_ipe, op_nrad
         op_nptot = nptot
         op_ipe = ipe
         op_nrad = nrad
      end subroutine get_op_mono_params
      
      
! HH: Based on "op_ax.f"
! Input:   kk = number of elements to calculate g_rad for   
!          iz1(kk) = charge of element to calculate g_rad for   
!          nel = number of elements in mixture   
!          izzp(nel) = charge of elements   
!          fap(nel) = number fractions of elements   
!          fac(nel) = scale factors for element opacity   
!          flux = local radiative flux (Lrad/4*pi*r^2)   
!          fltp = log10 T   
!          flrhop = log10 rho   
!          screening   if true, use screening corrections   
! Output: g1 = log10 kappa   
!         gx1 = d(log kappa)/d(log T)   
!         gy1 = d(log kappa)/d(log rho)   
!         gp1(kk) = d(log kappa)/d(log xi)    
!         grl1(kk) = log10 grad   
!         fx1(kk) = d(log grad)/d(log T)    
!         fy1(kk) = d(log grad)/d(log rho)   
!         grlp1(kk) = d(log grad)/d(log chi),
!              chi is the fraction with which the number fraction is varied, i.e.:
!                 chi = nf_new/nf_previous
!                 where nf is the number fraction
!         meanZ(nel) = average ionic charge of elements   
!         zetx1(nel) = d(meanZ)/d(log T)    
!         zety1(nel) = d(meanZ)/d(log rho)   
!         ierr = 0 for correct use
!         ierr = 101 for rho out of range for this T 
!         ierr = 102 for T out of range
      subroutine op_mono_get_radacc( &
            kk, izk, nel, izzp, fap, fac, flux, fltp, flrhop, screening, &
            g1, grl1, &
            umesh, semesh, ff, ta, rs, ierr)
         use kap_def
         use op_eval, only: eval_op_radacc
         integer, intent(in) :: kk, nel
         integer, intent(in) :: izk(kk), izzp(nel)
         real(dp), intent(in) :: fap(:) ! (nel) number fractions of elements
         real(dp), intent(in) :: fac(:) ! (nel) scale factors for element opacity
         real(dp), intent(in) :: flux, fltp
         real(dp), intent(in) :: flrhop
         logical, intent(in) :: screening
         real(dp), intent(out) :: g1
         real(dp), intent(inout) :: &
            grl1(kk)
         ! work arrays 
         real, pointer :: umesh(:), semesh(:), ff(:,:,:,:), ta(:,:,:,:), rs(:,:,:)
            ! umesh(nptot)
            ! semesh(nptot)
            ! ff(nptot, ipe, 4, 4)
            ! ta(nptot, nrad, 4, 4), 
            ! rs(nptot, 4, 4)
         integer,intent(out) :: ierr
         if (.not. kap_is_initialized) then
            ierr=-1
            return
         endif
         call eval_op_radacc( &
            kk, izk, nel, izzp, fap, fac, flux, fltp, flrhop, screening, &
            g1, grl1, &
            umesh, semesh, ff, ta, rs, ierr)
      end subroutine op_mono_get_radacc
      

! note: for op mono, elements must come from the set given in op_mono_element_Z in kap_def.
      
! HH: Based on "op_mx.f", opacity calculations to be used for stellar evolution calculations 
! Input:   nel = number of elements in mixture
!          izzp(nel) = charge of elements
!          fap(nel) = number fractions of elements
!          fac(nel) = scale factors for element opacity   
!          fltp = log10 (temperature)
!          flrhop = log10 (mass density) 
!          screening   if true, use screening corrections
! Output: g1 = log10 kappa
!         gx1 = d(log kappa)/d(log T)
!         gy1 = d(log kappa)/d(log rho)
!         ierr = 0 for correct use
!         ierr = 101 for rho out of range for this T 
!         ierr = 102 for T out of range
      subroutine op_mono_get_kap( &
            nel, izzp, fap, fac, fltp, flrhop, screening, &
            g1, gx1, gy1, &
            umesh, semesh, ff, rs, ierr)
         use kap_def
         use op_eval, only: eval_op_ev
         integer, intent(in) :: nel
         integer, intent(in) :: izzp(nel)
         real(dp), intent(in) :: fap(:) ! (nel) number fractions of elements
         real(dp), intent(in) :: fac(:) ! (nel) scale factors for element opacity
         real(dp), intent(in) :: fltp, flrhop
         logical, intent(in) :: screening
         real(dp), intent(inout) :: g1, gx1, gy1
         ! work arrays
         real, pointer :: umesh(:), semesh(:), ff(:,:,:,:), rs(:,:,:)
            ! umesh(nptot)
            ! semesh(nptot)
            ! ff(nptot, ipe, 4, 4)
            ! rs(nptot, 4, 4)
            ! s(nptot, nrad, 4, 4)
         integer,intent(out) :: ierr
         if (.not. kap_is_initialized) then
            ierr=-1
            return
         endif
         call eval_op_ev( &
            nel, izzp, fap, fac, fltp, flrhop, screening, &
            g1, gx1, gy1, &
            umesh, semesh, ff, rs, ierr)
      end subroutine op_mono_get_kap
      
      
! HH: Based on "op_mx.f", opacity calculations to be used for non-adiabatic pulsation calculations
! Special care is taken to ensure smoothness of opacity derivatives
! Input:   nel = number of elements in mixture
!          izzp(nel) = charge of elements
!          fap(nel) = number fractions of elements
!          fac(nel) = scale factors for element opacity   
!          fltp = log10 (temperature)
!          flrhop = log10 (mass density) 
!          screening   if true, use screening corrections
! Output: g1 = log10 kappa
!         gx1 = d(log kappa)/d(log T)
!         gy1 = d(log kappa)/d(log rho)
!         ierr = 0 for correct use
!         ierr = 101 for rho out of range for this T 
!         ierr = 102 for T out of range
      subroutine op_mono_alt_get_kap( &
            nel, izzp, fap, fac, fltp, flrhop, screening, &
            g1, gx1, gy1, &
            umesh, semesh, ff, rs, ierr)
         use kap_def
         use op_eval, only: eval_alt_op
         implicit none
         integer, intent(in) :: nel
         integer, intent(in) :: izzp(nel)
         real(dp), intent(in) :: fap(:) ! (nel) number fractions of elements
         real(dp), intent(in) :: fac(:) ! (nel) scale factors for element opacity
         real(dp), intent(in) :: fltp, flrhop
         logical, intent(in) :: screening
         real(dp), intent(out) :: g1, gx1, gy1
         ! work arrays
         real, pointer :: umesh(:), semesh(:), ff(:,:,:,:), rs(:,:,:)
            ! umesh(nptot)
            ! semesh(nptot)
            ! ff(nptot, ipe, 0:5, 0:5)
            ! rs(nptot, 0:5, 0:5)
         integer,intent(out) :: ierr
         if (.not. kap_is_initialized) then
            ierr=-1
            return
         endif
         call eval_alt_op( &
            nel, izzp, fap, fac, fltp, flrhop, screening, &
            g1, gx1, gy1, &
            umesh, semesh, ff, rs, ierr)
      end subroutine op_mono_alt_get_kap
      
      
      subroutine get_op_mono_args( &
            species, X, min_X_to_include, chem_id, chem_factors, &
            nel, izzp, fap, fac, ierr)
         use chem_def, only: chem_isos
         use kap_def
         integer, intent(in) :: species, chem_id(:)
         real(dp), intent(in) :: X(:) ! mass fractions (assumed baryonic)
         real(dp), intent(in) :: min_X_to_include ! skip iso if X < this
         real(dp), intent(in) :: chem_factors(:) ! (species)
         integer, intent(out) :: nel
         integer, intent(inout) :: izzp(:)
         real(dp), intent(inout) :: fap(:)
         real(dp), intent(inout) :: fac(:)
         integer,intent(out) :: ierr
         
         integer :: i, cid, j, Z, iel
         real(dp) :: tot
         
         ierr = 0
         if (.not. kap_is_initialized) then
            ierr=-1
            return
         endif
         
         nel = 0
         izzp(:) = 0
         fap(:) = 0d0
         
         do i=1,species
            if (X(i) < min_X_to_include) cycle
            cid = chem_id(i)
            Z = chem_isos% Z(cid)
            if (Z == 0) cycle
            ! change Z if necessary so that in op set
            do j=num_op_mono_elements,1,-1
               if (Z >= op_mono_element_Z(j)) then
                  Z = op_mono_element_Z(j)
                  exit
               end if
            end do
            iel = 0
            do j=1,nel
               if (izzp(j) == Z) then
                  iel = j
                  exit
               end if
            end do
            if (iel == 0) then
               nel = nel+1
               iel = nel
               izzp(nel) = Z
            end if
            fap(iel) = fap(iel) + X(i)/dble(chem_isos% Z_plus_N(cid))
            fac(iel) = chem_factors(i)
         end do
         
         tot = sum(fap(1:nel))
         if (tot <= 0d0) then
            ierr = -1
            return
         end if
         
         do j=1,nel
            fap(j) = fap(j)/tot ! number fractions
         end do
         
      end subroutine get_op_mono_args


      subroutine kap_get_control_namelist(handle, name, val, ierr)
         use kap_def
         use kap_ctrls_io, only: get_kap_controls
         integer, intent(in) :: handle ! kap handle; from star, pass s% kap_handle
         character(len=*),intent(in) :: name
         character(len=*),intent(out) :: val
         integer, intent(out) :: ierr
         type (kap_General_Info), pointer :: rq
         ierr = 0
         call kap_ptr(handle,rq,ierr)
         if(ierr/=0) return
         call get_kap_controls(rq, name, val, ierr)

      end subroutine kap_get_control_namelist

      subroutine kap_set_control_namelist(handle, name, val, ierr)
         use kap_def
         use kap_ctrls_io, only: set_kap_controls
         integer, intent(in) :: handle ! kap handle; from star, pass s% kap_handle
         character(len=*),intent(in) :: name
         character(len=*),intent(in) :: val
         integer, intent(out) :: ierr
         type (kap_General_Info), pointer :: rq
         ierr = 0
         call kap_ptr(handle,rq,ierr)
         if(ierr/=0) return
         call set_kap_controls(rq, name, val, ierr)

      end subroutine kap_set_control_namelist


      end module kap_lib


