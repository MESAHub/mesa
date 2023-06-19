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

      module alloc

      use star_private_def
      use const_def, only: ln10
      use utils_lib, only: &
         fill_with_NaNs, fill_with_NaNs_2D, fill_with_NaNs_3d, set_nan, &
         is_bad, mesa_error

      implicit none

      integer, parameter :: do_deallocate = 0
      integer, parameter :: do_allocate = 1
      integer, parameter :: do_check_size = 2
      integer, parameter :: do_remove_from_center = 3
      integer, parameter :: do_copy_pointers_and_resize = 4
      integer, parameter :: do_reallocate = 5
      integer, parameter :: do_fill_arrays_with_NaNs = 6

      logical, parameter :: work_array_debug = .false.
      logical, parameter :: work_array_trace = .false.

      logical, parameter :: quad_array_debug = .false.
      logical, parameter :: quad_array_trace = .false.


      ! working storage

      type work_array_pointer
         real(dp), dimension(:), pointer :: p
      end type work_array_pointer
      integer, parameter :: num_work_arrays = 250
      type (work_array_pointer), target :: &
         work_pointers(num_work_arrays)

      type quad_array_pointer
         real(qp), dimension(:), pointer :: p
      end type quad_array_pointer
      integer, parameter :: num_quad_arrays = 250
      type (quad_array_pointer), target :: &
         quad_pointers(num_quad_arrays)

      type int_work_array_pointer
         integer, dimension(:), pointer :: p
      end type int_work_array_pointer
      integer, parameter :: num_int_work_arrays = 250
      type (int_work_array_pointer), target :: &
         int_work_pointers(num_int_work_arrays)

      type logical_work_array_pointer
         logical, dimension(:), pointer :: p
      end type logical_work_array_pointer
      integer, parameter :: num_logical_work_arrays = 250
      type (logical_work_array_pointer), target :: &
         logical_work_pointers(num_logical_work_arrays)

      integer :: num_calls, num_returns
      integer :: num_allocs, num_deallocs



      contains

      
      subroutine init_alloc
         integer :: i
         num_calls=0; num_returns=0
         num_allocs=0; num_deallocs=0
         do i=1,num_work_arrays
            nullify(work_pointers(i)%p)
         end do
         do i=1,num_quad_arrays
            nullify(quad_pointers(i)%p)
         end do
         do i=1,num_int_work_arrays
            nullify(int_work_pointers(i)%p)
         end do
         do i=1,num_logical_work_arrays
            nullify(logical_work_pointers(i)%p)
         end do
      end subroutine init_alloc


      subroutine alloc_extras(id, liwork, lwork, ierr)
         integer, intent(in) :: liwork, lwork, id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return

         if (associated(s% extra_iwork)) deallocate(s% extra_iwork)
         if (associated(s% extra_iwork_old)) deallocate(s% extra_iwork_old)

         if (associated(s% extra_work)) deallocate(s% extra_work)
         if (associated(s% extra_work_old)) deallocate(s% extra_work_old)

         allocate( &
            s% extra_iwork(liwork), s% extra_work(lwork), &
            s% extra_iwork_old(liwork), s% extra_work_old(lwork), &
            stat=ierr)
         if (ierr == 0) then
            s% len_extra_iwork = liwork
            s% len_extra_work = lwork
            if (s% fill_arrays_with_NaNs) then
               call fill_with_NaNs(s% extra_work)
               call fill_with_NaNs(s% extra_work_old)
            else if (s% zero_when_allocate) then
               s% extra_work = 0
               s% extra_work_old = 0
            end if
         end if

      end subroutine alloc_extras


      subroutine dealloc_extras(s)
         type (star_info), pointer :: s
         if (associated(s% extra_iwork)) then
            deallocate(s% extra_iwork)
            nullify(s% extra_iwork)
         end if
         if (associated(s% extra_iwork_old)) then
            deallocate(s% extra_iwork_old)
            nullify(s% extra_iwork_old)
         end if
         s% len_extra_iwork = 0
         if (associated(s% extra_work)) then
            deallocate(s% extra_work)
            nullify(s% extra_work)
         end if
         if (associated(s% extra_work_old)) then
            deallocate(s% extra_work_old)
            nullify(s% extra_work_old)
         end if
         s% len_extra_work = 0
      end subroutine dealloc_extras


      subroutine update_nvar_allocs(s, nvar_hydro_old, nvar_chem_old, ierr)
         use utils_lib, only: realloc_double, realloc_double2, realloc_double3, realloc_quad2
         type (star_info), pointer :: s
         integer, intent(in) :: nvar_hydro_old, nvar_chem_old
         integer, intent(out) :: ierr

         integer :: nvar, nz, species, nvar_chem, nvar_hydro

         include 'formats'

         ierr = 0
         nvar_hydro = s% nvar_hydro
         nvar_chem = s% nvar_chem
         species = s% species

         if ((nvar_chem_old == nvar_chem) .and. (nvar_hydro_old == nvar_hydro)) return

         nvar = nvar_chem + nvar_hydro
         s% nvar_total = nvar
         nz = max(s% nz, s% prev_mesh_nz)

         if (nvar_chem_old == 0) return

         call realloc_double2(s% xh, nvar_hydro, nz + nz_alloc_extra, ierr)
         if (ierr /= 0) return

         call realloc_double2( &
            s% xh_old, nvar_hydro, s% nz_old + nz_alloc_extra, ierr)
         if (ierr /= 0) return

         call realloc_double(s% residual_weight1, nvar*(nz + nz_alloc_extra), ierr)
         if (ierr /= 0) return
         s% residual_weight(1:nvar,1:nz) => s% residual_weight1(1:nvar*nz)

         call realloc_double(s% correction_weight1, nvar*(nz + nz_alloc_extra), ierr)
         if (ierr /= 0) return
         s% correction_weight(1:nvar,1:nz) => s% correction_weight1(1:nvar*nz)

         call realloc_double(s% solver_dx1, nvar*(nz + nz_alloc_extra), ierr)
         if (ierr /= 0) return
         s% solver_dx(1:nvar,1:nz) => s% solver_dx1(1:nvar*nz)

         call realloc_double(s% x_scale1, nvar*(nz + nz_alloc_extra), ierr)
         if (ierr /= 0) return
         s% x_scale(1:nvar,1:nz) => s% x_scale1(1:nvar*nz)

         call realloc_double2(s% xh_start, nvar_hydro, (nz + nz_alloc_extra), ierr)
         if (ierr /= 0) return

         call realloc_double(s% xa_removed, species, ierr)
         if (ierr /= 0) return

         call realloc_double(s% op_mono_factors, species, ierr)
         if (ierr /= 0) return

         ! do prev_mesh arrays
         call realloc_double2(s% prev_mesh_xh, nvar_hydro, nz + nz_alloc_extra, ierr)
         if (ierr /= 0) return
         call realloc_double2(s% prev_mesh_xa, species, nz + nz_alloc_extra, ierr)
         if (ierr /= 0) return
         s% prev_mesh_species_or_nvar_hydro_changed = .true.

      end subroutine update_nvar_allocs


      subroutine free_star_data(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call free_arrays(s)
      end subroutine free_star_data


      subroutine free_arrays(s)
         use kap_lib, only: free_kap_handle
         use eos_lib, only: free_eos_handle
         use net_lib, only: free_net_handle
         use star_private_def, only: free_star
         use star_bcyclic, only: clear_storage
         type (star_info), pointer :: s
         integer :: ierr

         ! Free handles

         call free_net_handle(s% net_handle)
         s%net_handle = 0

         call free_eos_handle(s% eos_handle)
         s%eos_handle = 0

         call free_kap_handle(s% kap_handle)
         s%kap_handle = 0

         ! Free star_info arrays

         call star_info_arrays(s, null(), do_deallocate, ierr)
         call star_info_old_arrays(s, do_deallocate, ierr)

         if (ASSOCIATED(s%xa_removed)) then
            deallocate(s%xa_removed)
            nullify(s%xa_removed)
         end if
         if (ASSOCIATED(s%op_mono_factors)) then
            deallocate(s%op_mono_factors)
            nullify(s%op_mono_factors)
         end if
            
         call dealloc_extras(s)

         call free_other(s)

         if (ASSOCIATED(s%other_star_info)) then

            nullify(s%other_star_info%profile_extra)

            call star_info_arrays(s%other_star_info, null(), do_deallocate, ierr)

            deallocate(s%other_star_info)

         endif
         
         call dealloc_history(s)
         
         if (ASSOCIATED(s%bcyclic_odd_storage)) call clear_storage(s)

         ! Free the star handle itself

         call free_star(s)
         
      end subroutine free_arrays


      subroutine free_other (s)

        type (star_info), pointer :: s

         ! Pointer aliases (do not deallocate, as that's handled elsewhere)

         if (ASSOCIATED(s% chem_id)) nullify(s% chem_id)
         if (ASSOCIATED(s% net_iso)) nullify(s% net_iso)

         ! Misc arrays not handled in star_info_arrays and related
         ! routines. (These are "found" by judicious application of
         ! valgrind)

         if (ASSOCIATED(s% conv_bdy_q)) deallocate(s% conv_bdy_q)
         if (ASSOCIATED(s% conv_bdy_loc)) deallocate(s% conv_bdy_loc)
         if (ASSOCIATED(s% top_conv_bdy)) deallocate(s% top_conv_bdy)
         if (ASSOCIATED(s% burn_h_conv_region)) deallocate(s% burn_h_conv_region)
         if (ASSOCIATED(s% burn_he_conv_region)) deallocate(s% burn_he_conv_region)
         if (ASSOCIATED(s% burn_z_conv_region)) deallocate(s% burn_z_conv_region)

         if (ASSOCIATED(s% prev_mesh_xa)) deallocate(s% prev_mesh_xa)
         if (ASSOCIATED(s% prev_mesh_xh)) deallocate(s% prev_mesh_xh)
         if (ASSOCIATED(s% prev_mesh_j_rot)) deallocate(s% prev_mesh_j_rot)
         if (ASSOCIATED(s% prev_mesh_omega)) deallocate(s% prev_mesh_omega)
         if (ASSOCIATED(s% prev_mesh_dq)) deallocate(s% prev_mesh_dq)

         if (ASSOCIATED(s% mix_bdy_q)) deallocate(s% mix_bdy_q)
         if (ASSOCIATED(s% mix_bdy_loc)) deallocate(s% mix_bdy_loc)
         if (ASSOCIATED(s% top_mix_bdy)) deallocate(s% top_mix_bdy)
         if (ASSOCIATED(s% burn_h_mix_region)) deallocate(s% burn_h_mix_region)
         if (ASSOCIATED(s% burn_he_mix_region)) deallocate(s% burn_he_mix_region)
         if (ASSOCIATED(s% burn_z_mix_region)) deallocate(s% burn_z_mix_region)

         if (ASSOCIATED(s% rate_factors)) deallocate(s% rate_factors)

         if (ASSOCIATED(s% nameofvar)) deallocate(s% nameofvar)
         if (ASSOCIATED(s% nameofequ)) deallocate(s% nameofequ)

         if (ASSOCIATED(s% solver_work)) deallocate(s% solver_work)
         if (ASSOCIATED(s% solver_iwork)) deallocate(s% solver_iwork)

         if (ASSOCIATED(s% AF1)) deallocate(s% AF1)

      end subroutine free_other


      subroutine check_sizes(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         type (star_info), pointer :: c => null()
         call star_info_arrays(s, c, do_check_size, ierr)
         if (ierr /= 0) then
            write(*,*) 'check_sizes failed for s'
            return
         end if
         if (s% generations <= 1) return
         call star_info_old_arrays(s, do_check_size, ierr)
         if (ierr /= 0) then
            write(*,*) 'check_sizes failed for s old'
            return
         end if
      end subroutine check_sizes


      subroutine alloc_star_info_old_arrays(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         call star_info_old_arrays(s, do_allocate, ierr)
      end subroutine alloc_star_info_old_arrays


      subroutine star_info_old_arrays(s, action, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: action
         integer, intent(out) :: ierr

         integer :: nz, species, nvar_hydro
         
         include 'formats'

         nz = s% nz_old
         nvar_hydro = s% nvar_hydro
         species = s% species
         ierr = -1
         call do2D(s, s% xh_old, nvar_hydro, nz, action, ierr)
         if (failed('xh_old')) return
         call do2D(s, s% xa_old, species, nz, action, ierr)
         if (failed('xa_old')) return
         call do1D(s, s% q_old, nz, action, ierr)
         if (failed('q_old')) return
         call do1D(s, s% dq_old, nz, action, ierr)
         if (failed('dq_old')) return
         call do1D(s, s% mlt_vc_old, nz, action, ierr)
         if (failed('mlt_vc_old')) return
         call do1D(s, s% omega_old, nz, action, ierr)
         if (failed('omega_old')) return
         call do1D(s, s% j_rot_old, nz, action, ierr)
         if (failed('j_rot_old')) return
         ierr = 0

         contains

         logical function failed(str)
            character (len=*), intent(in) :: str
            failed = .false.
            if (ierr == 0) return
            write(*,*) 'star_info_old_arrays failed for ' // trim(str)
            failed = .true.
         end function failed

      end subroutine star_info_old_arrays


      subroutine free_star_info_arrays(s)
         type (star_info), pointer :: s
         integer :: ierr
         type (star_info), pointer :: c => null()
         call star_info_arrays(s, c, do_deallocate, ierr)
         if (ierr /= 0) write(*,*) 'free_star_info_arrays failed'
      end subroutine free_star_info_arrays


      subroutine allocate_star_info_arrays(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         type (star_info), pointer :: c => null()
         call star_info_arrays(s, c, do_allocate, ierr)
         if (ierr /= 0) write(*,*) 'allocate_star_info_arrays failed'
      end subroutine allocate_star_info_arrays


      subroutine prune_star_info_arrays(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         type (star_info), pointer :: c => null()
         call star_info_arrays(s, c, do_remove_from_center, ierr)
         if (ierr /= 0) write(*,*) 'prune_star_info_arrays failed'
      end subroutine prune_star_info_arrays


      subroutine resize_star_info_arrays(s, c, ierr)
         type (star_info), pointer :: s, c
         integer, intent(out) :: ierr
         call star_info_arrays(s, c, do_copy_pointers_and_resize, ierr)
         if (ierr /= 0) write(*,*) 'resize_star_info_arrays failed'
      end subroutine resize_star_info_arrays


      subroutine reallocate_star_info_arrays(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         call star_info_arrays(s, null(), do_reallocate, ierr)
         if (ierr /= 0) write(*,*) 'reallocate_star_info_arrays failed'
      end subroutine reallocate_star_info_arrays


      subroutine fill_star_info_arrays_with_NaNs(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         call star_info_arrays(s, null(), do_fill_arrays_with_NaNs, ierr)
         if (ierr /= 0) write(*,*) 'fill_star_info_arrays_with_NaNs failed'
         s% need_to_setvars = .true.
      end subroutine fill_star_info_arrays_with_NaNs


      subroutine star_info_arrays(s, c_in, action_in, ierr)
         use eos_def, only:num_eos_basic_results, num_eos_d_dxa_results
         use rates_def, only: num_rvs
         use chem_def, only: num_categories
         use const_def, only: standard_cgrav
         type (star_info), pointer :: s, c_in
         integer, intent(in) :: action_in
         integer, intent(out) :: ierr

         integer :: nz, species, num_reactions, &
            nvar, nvar_hydro, nvar_chem, sz_new, psz_new, action
         type (star_info), pointer :: c
         character (len=128) :: null_str

         include 'formats'

         ierr = 0
         null_str = '' ! avoid bogus compiler warnings 'array subscript 1 is above array bounds'
         
         
         species = s% species
         num_reactions = s% num_reactions
         nvar = s% nvar_total
         nvar_hydro = s% nvar_hydro
         nvar_chem = s% nvar_chem

         c => s
         action = action_in
         if (action == do_check_size) then
            nz = s% nz
            sz_new = nz
         else
            nz = max(s% prev_mesh_nz, s% nz)
            sz_new = nz + nz_alloc_extra
         end if

         if (action == do_copy_pointers_and_resize) then
            if (associated(c_in)) then
               c => c_in
            else ! nothing to copy, so switch to allocate
               action = do_allocate
            end if
         end if

         do ! just so can exit on failure

            if (action /= do_fill_arrays_with_NaNs) then
               ! these arrays must not be filled with NaNs  
               ! because they contain the inputs to the step          
               call do2(s% xh, c% xh, nvar_hydro, 'xh')
               if (failed('xh')) exit
               call do2(s% xa, c% xa, species, 'xa')
               if (failed('xa')) then
                  write(*,2) 'species', species
                  write(*,2) 'size(s% xa,dim=1)', size(s% xa,dim=1)
                  write(*,2) 's% nz', s% nz
                  write(*,2) 's% prev_mesh_nz', s% prev_mesh_nz
                  write(*,2) 'size(s% xa,dim=2)', size(s% xa,dim=2)
                  call mesa_error(__FILE__,__LINE__,'star_info_arrays')
                  exit
               end if
               call do1(s% dq, c% dq)
               if (failed('dq')) exit
               call do1(s% omega, c% omega)
               if (failed('omega')) exit
               call do1(s% j_rot, c% j_rot)
               if (failed('j_rot')) exit               
               call do1(s% mlt_vc, c% mlt_vc)
               if (failed('mlt_vc')) exit
               call do1(s% conv_vel, c% conv_vel)
               if (failed('conv_vel')) exit
            end if
            
            call do1(s% q, c% q)
            if (failed('q')) exit
            call do1(s% m, c% m)
            if (failed('m')) exit
            call do1(s% dm, c% dm)
            if (failed('dm')) exit
            call do1(s% dm_bar, c% dm_bar)
            if (failed('dm_bar')) exit   
            
            call do1(s% am_nu_rot, c% am_nu_rot)
            if (failed('am_nu_rot')) exit
            call do1(s% D_omega, c% D_omega)
            if (failed('D_omega')) exit
            call do1_ad(s% fp_rot, c% fp_rot)
            if (failed('fp_rot')) exit
            call do1_ad(s% ft_rot, c% ft_rot)
            if (failed('ft_rot')) exit
            call do1(s% am_nu_non_rot, c% am_nu_non_rot)
            if (failed('am_nu_non_rot')) exit
            call do1(s% am_nu_omega, c% am_nu_omega)
            if (failed('am_nu_omega')) exit
            call do1(s% am_nu_j, c% am_nu_j)
            if (failed('am_nu_j')) exit
            call do1(s% am_sig_omega, c% am_sig_omega)
            if (failed('am_sig_omega')) exit
            call do1(s% am_sig_j, c% am_sig_j)
            if (failed('am_sig_j')) exit
            call do1_ad(s% i_rot, c% i_rot)
            if (failed('i_rot')) exit
            call do1(s% w_div_w_crit_roche, c% w_div_w_crit_roche)
            if (failed('w_div_w_crit_roche')) exit
            call do1_ad(s% j_flux, c% j_flux)
            if (failed('j_flux')) exit

            call do2(s% xh_start, c% xh_start, nvar_hydro, 'xh_start')
            if (failed('xh_start')) exit
            
            call do1(s% r_polar, c% r_polar)
            if (failed('r_polar')) exit
            call do1(s% r_equatorial, c% r_equatorial)
            if (failed('r_equatorial')) exit               

            call do1(s% lnd, c% lnd)
            if (failed('lnd')) exit
            call do1(s% lnT, c% lnT)
            if (failed('lnT')) exit
            call do1(s% lnR, c% lnR)
            if (failed('lnR')) exit
            call do1(s% RSP_Et, c% RSP_Et)
            if (failed('RSP_Et')) exit
            call do1(s% L, c% L)
            if (failed('L')) exit
            call do1(s% v, c% v)
            if (failed('v')) exit
            call do1(s% u, c% u)
            if (failed('u')) exit
            call do1(s% alpha_RTI, c% alpha_RTI)
            if (failed('alpha_RTI')) exit
            call do1(s% w, c% w)
            if (failed('w')) exit
            call do1(s% w_start, c% w_start)
            if (failed('w_start')) exit
            call do1(s% Hp_face_start, c% Hp_face_start)
            if (failed('Hp_face_start')) exit
            
            call do1(s% dxh_lnR, c% dxh_lnR)
            if (failed('dxh_lnR')) exit
            call do1(s% dxh_lnd, c% dxh_lnd)
            if (failed('dxh_lnd')) exit
            call do1(s% dxh_lnT, c% dxh_lnT)
            if (failed('dxh_lnT')) exit
            call do1(s% dxh_v, c% dxh_v)
            if (failed('dxh_v')) exit
            call do1(s% dxh_u, c% dxh_u)
            if (failed('dxh_u')) exit
            call do1(s% dxh_alpha_RTI, c% dxh_alpha_RTI)
            if (failed('dxh_alpha_RTI')) exit

            call do1(s% dudt_RTI, c% dudt_RTI)
            if (failed('dudt_RTI')) exit
            call do1(s% dedt_RTI, c% dedt_RTI)
            if (failed('dedt_RTI')) exit

            call do1(s% T, c% T)
            if (failed('T')) exit
            call do1(s% rho, c% rho)
            if (failed('rho')) exit
            call do1(s% r, c% r)
            if (failed('r')) exit

            call do1(s% rmid, c% rmid)
            if (failed('rmid')) exit

            call do1(s% X, c% X)
            if (failed('X')) exit
            call do1(s% Y, c% Y)
            if (failed('Y')) exit
            call do1(s% Z, c% Z)
            if (failed('Z')) exit
            call do1(s% abar, c% abar)
            if (failed('abar')) exit
            call do1(s% zbar, c% zbar)
            if (failed('zbar')) exit
            call do1(s% z2bar, c% z2bar)
            if (failed('z2bar')) exit
            call do1(s% z53bar, c% z53bar)
            if (failed('z53bar')) exit
            call do1(s% ye, c% ye)
            if (failed('ye')) exit

            call do1(s% mass_correction, c% mass_correction)
            if (failed('mass_correction')) exit
            call do1(s% mass_correction_start, c% mass_correction_start)
            if (failed('mass_correction_start')) exit
            call do1(s% m_grav, c% m_grav)
            if (failed('m_grav')) exit
            call do1(s% m_grav_start, c% m_grav_start)
            if (failed('m_grav_start')) exit

            call do1(s% Peos, c% Peos)
            if (failed('Peos')) exit
            call do1(s% lnPeos, c% lnPeos)
            if (failed('lnPeos')) exit
            call do1(s% lnPgas, c% lnPgas)
            if (failed('lnPgas')) exit
            call do1(s% Pgas, c% Pgas)
            if (failed('Pgas')) exit
            call do1(s% Prad, c% Prad)
            if (failed('Prad')) exit
            call do1(s% energy, c% energy)
            if (failed('energy')) exit
            call do1(s% egas, c% egas)
            if (failed('egas')) exit
            call do1(s% erad, c% erad)
            if (failed('erad')) exit
            call do1(s% lnE, c% lnE)
            if (failed('lnE')) exit
            call do1(s% grada, c% grada)
            if (failed('grada')) exit
            call do1(s% dE_dRho, c% dE_dRho)
            if (failed('dE_dRho')) exit
            call do1(s% Cv, c% Cv)
            if (failed('Cv')) exit
            call do1(s% Cp, c% Cp)
            if (failed('Cp')) exit
            call do1(s% lnS, c% lnS)
            if (failed('lnS')) exit
            call do1(s% gamma1, c% gamma1)
            if (failed('gamma1')) exit
            call do1(s% gamma3, c% gamma3)
            if (failed('gamma3')) exit
            call do1(s% eta, c% eta)
            if (failed('eta')) exit
            call do1(s% gam, c% gam)
            if (failed('gam')) exit
            call do1(s% mu, c% mu)
            if (failed('mu')) exit
            call do1(s% lnfree_e, c% lnfree_e)
            if (failed('lnfree_e')) exit

            call do1(s% phase, c% phase)
            if (failed('phase')) exit
            call do1(s% latent_ddlnT, c% latent_ddlnT)
            if (failed('latent_ddlnT')) exit
            call do1(s% latent_ddlnRho, c% latent_ddlnRho)
            if (failed('latent_ddlnRho')) exit

            call do1(s% chiRho, c% chiRho)
            if (failed('chiRho')) exit
            call do1(s% chiT, c% chiT)
            if (failed('chiT')) exit

            call do1(s% eos_frac_OPAL_SCVH, c% eos_frac_OPAL_SCVH)
            if (failed('eos_frac_OPAL_SCVH')) exit
            call do1(s% eos_frac_HELM, c% eos_frac_HELM)
            if (failed('eos_frac_HELM')) exit
            call do1(s% eos_frac_Skye, c% eos_frac_Skye)
            if (failed('eos_frac_Skye')) exit
            call do1(s% eos_frac_PC, c% eos_frac_PC)
            if (failed('eos_frac_PC')) exit
            call do1(s% eos_frac_FreeEOS, c% eos_frac_FreeEOS)
            if (failed('eos_frac_FreeEOS')) exit
            call do1(s% eos_frac_CMS, c% eos_frac_CMS)
            if (failed('eos_frac_CMS')) exit
            call do1(s% eos_frac_ideal, c% eos_frac_ideal)
            if (failed('eos_frac_ideal')) exit

            call do1(s% QQ, c% QQ)
            if (failed('QQ')) exit
            call do2(s% d_eos_dlnd, c% d_eos_dlnd, num_eos_basic_results, 'd_eos_dlnd')
            if (failed('d_eos_dlnd')) exit
            call do2(s% d_eos_dlnT, c% d_eos_dlnT, num_eos_basic_results, 'd_eos_dlnT')
            if (failed('d_eos_dlnT')) exit
            call do3(s% d_eos_dxa, c% d_eos_dxa, num_eos_d_dxa_results, species)
            if (failed('d_eos_dxa')) exit

            call do1(s% chiRho_for_partials, c% chiRho_for_partials)
            if (failed('chiRho_for_partials')) exit
            call do1(s% chiT_for_partials, c% chiT_for_partials)
            if (failed('chiT_for_partials')) exit
            call do1(s% dE_dRho_for_partials, c% dE_dRho_for_partials)
            if (failed('dE_dRho_for_partials')) exit
            call do1(s% Cv_for_partials, c% Cv_for_partials)
            if (failed('Cv_for_partials')) exit
            call do1(s% dS_dRho_for_partials, c% dS_dRho_for_partials)
            if (failed('dS_dRho_for_partials')) exit
            call do1(s% dS_dT_for_partials, c% dS_dT_for_partials)
            if (failed('dS_dT_for_partials')) exit
            call do2(s% dlnE_dxa_for_partials, c% dlnE_dxa_for_partials, species, null_str)
            if (failed('dlnE_dxa_for_partials')) exit
            call do2(s% dlnPeos_dxa_for_partials, c% dlnPeos_dxa_for_partials, species, null_str)
            if (failed('dlnPeos_dxa_for_partials')) exit

            ! other model variables
            call do1(s% csound, c% csound)
            if (failed('csound')) exit
            call do1(s% csound_face, c% csound_face)
            if (failed('csound_face')) exit
            
            call do1(s% rho_face, c% rho_face)
            if (failed('rho_face')) exit
            
            call do1(s% scale_height, c% scale_height)
            if (failed('scale_height')) exit
            call do1(s% v_div_csound, c% v_div_csound)
            if (failed('v_div_csound')) exit
            call do1(s% entropy, c% entropy)
            if (failed('entropy')) exit
            call do1(s% grav, c% grav)
            if (failed('grav')) exit
            call do1(s% tau, c% tau)
            if (failed('tau')) exit
            call do1(s% dr_div_csound, c% dr_div_csound)
            if (failed('dr_div_csound')) exit

            call do1(s% ergs_error, c% ergs_error)
            if (failed('ergs_error')) exit            
            call do1(s% gradr_factor, c% gradr_factor)
            if (failed('gradr_factor')) exit
            call do1(s% adjust_mlt_gradT_fraction, c% adjust_mlt_gradT_fraction)
            if (failed('adjust_mlt_gradT_fraction')) exit
            call do1(s% gradT_excess_effect, c% gradT_excess_effect)
            if (failed('gradT_excess_effect')) exit
            call do1(s% superad_reduction_factor, c% superad_reduction_factor)
            if (failed('superad_reduction_factor')) exit

            call do1(s% domega_dlnR, c% domega_dlnR)
            if (failed('domega_dlnR')) exit
            call do1(s% richardson_number, c% richardson_number)
            if (failed('richardson_number')) exit

            call do1(s% D_mix_rotation, c% D_mix_rotation)
            if (failed('D_mix_rotation')) exit
            call do1(s% D_mix_non_rotation, c% D_mix_non_rotation)
            if (failed('D_mix_non_rotation')) exit
            call do1(s% D_visc, c% D_visc)
            if (failed('D_visc')) exit
            call do1(s% D_DSI, c% D_DSI)
            if (failed('D_DSI')) exit
            call do1(s% D_SH, c% D_SH)
            if (failed('D_SH')) exit
            call do1(s% D_SSI, c% D_SSI)
            if (failed('D_SSI')) exit
            call do1(s% D_ES, c% D_ES)
            if (failed('D_ES')) exit
            call do1(s% D_GSF, c% D_GSF)
            if (failed('D_GSF')) exit

            call do1(s% D_ST, c% D_ST)
            if (failed('D_ST')) exit
            call do1(s% nu_ST, c% nu_ST)
            if (failed('nu_ST')) exit
            call do1(s% omega_shear, c% omega_shear)
            if (failed('omega_shear')) exit

            call do1(s% dynamo_B_r, c% dynamo_B_r)
            if (failed('dynamo_B_r')) exit
            call do1(s% dynamo_B_phi, c% dynamo_B_phi)
            if (failed('dynamo_B_phi')) exit

            !for ST time smoothing
            call do1(s% D_ST_start, c% D_ST_start)
            if (failed('D_ST_start')) exit
            call do1(s% nu_ST_start, c% nu_ST_start)
            if (failed('nu_ST_start')) exit

            call do1(s% opacity, c% opacity)
            if (failed('opacity')) exit
            call do1(s% d_opacity_dlnd, c% d_opacity_dlnd)
            if (failed('d_opacity_dlnd')) exit
            call do1(s% d_opacity_dlnT, c% d_opacity_dlnT)
            if (failed('d_opacity_dlnT')) exit
            call do1(s% kap_frac_lowT, c% kap_frac_lowT)
            if (failed('kap_frac_lowT')) exit
            call do1(s% kap_frac_highT, c% kap_frac_highT)
            if (failed('kap_frac_highT')) exit
            call do1(s% kap_frac_Type2, c% kap_frac_Type2)
            if (failed('kap_frac_Type2')) exit
            call do1(s% kap_frac_Compton, c% kap_frac_Compton)
            if (failed('kap_frac_Compton')) exit
            call do1(s% kap_frac_op_mono, c% kap_frac_op_mono)
            if (failed('kap_frac_op_mono')) exit

            call do1(s% eps_nuc, c% eps_nuc)
            if (failed('eps_nuc')) exit
            call do1(s% d_epsnuc_dlnd, c% d_epsnuc_dlnd)
            if (failed('d_epsnuc_dlnd')) exit
            call do1(s% d_epsnuc_dlnT, c% d_epsnuc_dlnT)
            if (failed('d_epsnuc_dlnT')) exit
            call do2(s% d_epsnuc_dx, c% d_epsnuc_dx, species, null_str)
            if (failed('d_epsnuc_dx')) exit

            call do1(s% eps_nuc_neu_total, c% eps_nuc_neu_total)
            if (failed('eps_nuc_neu_total')) exit

            call do2(s% raw_rate, c% raw_rate, num_reactions, null_str)
            if (failed('raw_rates')) exit
            call do2(s% screened_rate, c% screened_rate, num_reactions, null_str)
            if (failed('screened_rates')) exit
            call do2(s% eps_nuc_rate, c% eps_nuc_rate, num_reactions, null_str)
            if (failed('eps_nuc_rates')) exit
            call do2(s% eps_neu_rate, c% eps_neu_rate, num_reactions, null_str)
            if (failed('eps_neu_rates')) exit

            call do2(s% dxdt_nuc, c% dxdt_nuc, species, null_str)
            if (failed('dxdt_nuc')) exit
            call do2(s% d_dxdt_nuc_dRho, c% d_dxdt_nuc_dRho, species, null_str)
            if (failed('d_dxdt_nuc_dRho')) exit
            call do2(s% d_dxdt_nuc_dT, c% d_dxdt_nuc_dT, species, null_str)
            if (failed('d_dxdt_nuc_dT')) exit
            call do3(s% d_dxdt_nuc_dx, c% d_dxdt_nuc_dx, species, species)
            if (failed('d_dxdt_nuc_dx')) exit

            call do2(s% dxdt_mix, c% dxdt_mix, species, null_str)
            if (failed('dxdt_mix')) exit
            call do1(s% d_dxdt_mix_dxp1, c% d_dxdt_mix_dxp1)
            if (failed('d_dxdt_mix_dRho')) exit
            call do1(s% d_dxdt_mix_dx00, c% d_dxdt_mix_dx00)
            if (failed('d_dxdt_mix_dx00')) exit
            call do1(s% d_dxdt_mix_dxm1, c% d_dxdt_mix_dxm1)
            if (failed('d_dxdt_mix_dxm1')) exit

            if (action /= do_check_size) then
               call do2(s% eps_nuc_categories, c% eps_nuc_categories, num_categories, null_str)
               if (failed('eps_nuc_categories')) exit
               call do2(s% luminosity_by_category, c% luminosity_by_category, num_categories, null_str)
               if (failed('luminosity_by_category')) exit
               call do2(s% luminosity_by_category_start, &
                  c% luminosity_by_category_start, num_categories, null_str)
               if (failed('luminosity_by_category_start')) exit
               call do2(s% L_nuc_by_category, c% L_nuc_by_category, num_categories, null_str)
               if (failed('L_nuc_by_category')) exit
            end if

            call do2(s% diffusion_D_self, c% diffusion_D_self, species, null_str)
            if (failed('diffusion_D_self')) exit
            call do2(s% extra_diffusion_factor, c% extra_diffusion_factor, species, null_str)
            if (failed('extra_diffusion_factor')) exit
            call do2(s% edv, c% edv, species, null_str)
            if (failed('edv')) exit
            call do2(s% v_rad, c% v_rad, species, null_str)
            if (failed('v_rad')) exit
            call do2(s% g_rad, c% g_rad, species, null_str)
            if (failed('g_rad')) exit
            call do2(s% typical_charge, c% typical_charge, species, null_str)
            if (failed('typical_charge')) exit
            call do2(s% diffusion_dX, c% diffusion_dX, species, null_str)
            if (failed('diffusion_dX')) exit
            call do1(s% E_field, c% E_field)
            if (failed('E_field')) exit
            call do1(s% eps_WD_sedimentation, c% eps_WD_sedimentation)
            if (failed('eps_WD_sedimentation')) exit
            call do1(s% eps_diffusion, c% eps_diffusion)
            if (failed('eps_diffusion')) exit
            call do1(s% g_field_element_diffusion, c% g_field_element_diffusion)
            if (failed('g_field_element_diffusion')) exit

            call do1(s% eps_phase_separation, c% eps_phase_separation)
            if (failed('eps_phase_separation')) exit

            call do1(s% non_nuc_neu, c% non_nuc_neu)
            if (failed('non_nuc_neu')) exit
            call do1(s% d_nonnucneu_dlnd, c% d_nonnucneu_dlnd)
            if (failed('d_nonnucneu_dlnd')) exit
            call do1(s% d_nonnucneu_dlnT, c% d_nonnucneu_dlnT)
            if (failed('d_nonnucneu_dlnT')) exit

            call do1(s% nonnucneu_plas, c% nonnucneu_plas)
            if (failed('nonnucneu_plas')) exit
            call do1(s% nonnucneu_brem, c% nonnucneu_brem)
            if (failed('nonnucneu_brem')) exit
            call do1(s% nonnucneu_phot, c% nonnucneu_phot)
            if (failed('nonnucneu_phot')) exit
            call do1(s% nonnucneu_pair, c% nonnucneu_pair)
            if (failed('nonnucneu_pair')) exit
            call do1(s% nonnucneu_reco, c% nonnucneu_reco)
            if (failed('nonnucneu_reco')) exit

            call do1(s% extra_opacity_factor, c% extra_opacity_factor)
            if (failed('extra_opacity_factor')) exit

            call do1_ad(s% extra_pressure, c% extra_pressure)
            if (failed('extra_pressure')) exit
            call do1(s% eps_heat, c% eps_heat)
            if (failed('eps_heat')) exit
            call do1(s% irradiation_heat, c% irradiation_heat)
            if (failed('irradiation_heat')) exit
            
            call do1_ad(s% extra_heat, c% extra_heat)
            if (failed('extra_heat')) exit
            call do1_ad(s% extra_grav, c% extra_grav)
            if (failed('extra_grav')) exit

            call do1(s% extra_jdot, c% extra_jdot)
            if (failed('extra_jdot')) exit
            call do1(s% extra_omegadot, c% extra_omegadot)
            if (failed('extra_omegadot')) exit

            call do1(s% d_extra_jdot_domega_m1, c% d_extra_jdot_domega_m1)
            if (failed('d_extra_jdot_domega_m1')) exit
            call do1(s% d_extra_omegadot_domega_m1, c% d_extra_omegadot_domega_m1)
            if (failed('d_extra_omegadot_domega_m1')) exit
            call do1(s% d_extra_jdot_domega_00, c% d_extra_jdot_domega_00)
            if (failed('d_extra_jdot_domega_00')) exit
            call do1(s% d_extra_omegadot_domega_00, c% d_extra_omegadot_domega_00)
            if (failed('d_extra_omegadot_domega_00')) exit
            call do1(s% d_extra_jdot_domega_p1, c% d_extra_jdot_domega_p1)
            if (failed('d_extra_jdot_domega_p1')) exit
            call do1(s% d_extra_omegadot_domega_p1, c% d_extra_omegadot_domega_p1)
            if (failed('d_extra_omegadot_domega_p1')) exit

            call do1(s% cgrav, c% cgrav)
            if (failed('cgrav')) exit
            if (action == do_allocate .or. &
                  action == do_copy_pointers_and_resize) &
               s% cgrav(1:nz) = standard_cgrav

            call do1(s% mesh_delta_coeff_factor, c% mesh_delta_coeff_factor)
            if (failed('mesh_delta_coeff_factor')) exit
            if (action == do_allocate .or. &
                  action == do_copy_pointers_and_resize) &
               s% mesh_delta_coeff_factor(1:nz) = 1d0

            call do1_logical(s% amr_split_merge_has_undergone_remesh, c% amr_split_merge_has_undergone_remesh)
            if (failed('amr_split_merge_has_undergone_remesh')) exit

            call do1(s% alpha_mlt, c% alpha_mlt)
            if (failed('alpha_mlt')) exit
            if (action == do_allocate .or. &
                  action == do_copy_pointers_and_resize) &
               s% alpha_mlt(1:nz) = s% mixing_length_alpha

            call do1(s% vc, c% vc)
            if (failed('vc')) exit

            call do1_ad(s% eps_grav_ad, c% eps_grav_ad)
            if (failed('eps_grav_ad')) exit
            call do2(s% d_eps_grav_dx, c% d_eps_grav_dx, species, null_str)
            if (failed('d_eps_grav_dx')) exit

            call do1(s% eps_grav_composition_term, c% eps_grav_composition_term)
            if (failed('eps_grav_composition_term')) exit

            call do1(s% eps_mdot, c% eps_mdot)
            if (failed('eps_mdot')) exit
            call do1(s% dm_before_adjust_mass, c% dm_before_adjust_mass)
            if (failed('dm_before_adjust_mass')) exit
            call do1(s% total_energy_profile_before_adjust_mass, c% total_energy_profile_before_adjust_mass)
            if (failed('total_energy_profile_before_adjust_mass')) exit
            call do1(s% total_energy_profile_after_adjust_mass, c% total_energy_profile_after_adjust_mass)
            if (failed('total_energy_profile_after_adjust_mass')) exit

            call do1(s% dL_dm, c% dL_dm)
            if (failed('dL_dm')) exit
            call do1(s% energy_sources, c% energy_sources)
            if (failed('energy_sources')) exit
            call do1(s% energy_others, c% energy_others)
            if (failed('energy_others')) exit
            call do1(s% dwork_dm, c% dwork_dm)
            if (failed('dwork_dm')) exit
            call do1(s% dkedt, c% dkedt)
            if (failed('dkedt')) exit
            call do1(s% dpedt, c% dpedt)
            if (failed('dpedt')) exit
            call do1(s% dedt, c% dedt)
            if (failed('dedt')) exit
            call do1(s% detrbdt, c% detrbdt)
            if (failed('detrbdt')) exit

            call do1_integer(s% mlt_mixing_type, c% mlt_mixing_type)
            if (failed('mlt_mixing_type')) exit
            call do1(s% mlt_mixing_length, c% mlt_mixing_length)
            if (failed('mlt_mixing_length')) exit
            call do1(s% mlt_D, c% mlt_D)
            if (failed('mlt_D')) exit

            call do1_logical(s% fixed_gradr_for_rest_of_solver_iters, c% fixed_gradr_for_rest_of_solver_iters)
            if (failed('fixed_gradr_for_rest_of_solver_iters')) exit
            
            call do1(s% mlt_Gamma, c% mlt_Gamma)
            if (failed('mlt_Gamma')) exit
            call do1(s% L_conv, c% L_conv)
            if (failed('L_conv')) exit

            call do1_integer(s% tdc_num_iters, c% tdc_num_iters)
            if (failed('tdc_num_iters')) exit

            call do1(s% gradT_sub_grada, c% gradT_sub_grada)
            if (failed('gradT_sub_grada')) exit
            call do1(s% grada_face, c% grada_face)
            if (failed('grada_face')) exit
            call do1(s% mlt_gradT, c% mlt_gradT)
            if (failed('mlt_gradT')) exit

            call do1(s% mlt_cdc, c% mlt_cdc)
            if (failed('mlt_cdc')) exit
            call do1(s% cdc, c% cdc)
            if (failed('cdc')) exit
            call do1(s% dvc_dt_TDC, c% dvc_dt_TDC)
            if (failed('dvc_dt_TDC')) exit

            call do1(s% D_mix, c% D_mix)
            if (failed('D_mix')) exit

            call do1_integer(s% mixing_type, c% mixing_type)
            if (failed('mixing_type')) exit
            call do1(s% cz_bdy_dq, c% cz_bdy_dq)
            if (failed('cz_bdy_dq')) exit
            
            call do1_ad(s% gradT_ad, c% gradT_ad)
            if (failed('gradT_ad')) exit
            call do1_ad(s% gradr_ad, c% gradr_ad)
            if (failed('gradr_ad')) exit
            call do1_ad(s% Y_face_ad, c% Y_face_ad)
            if (failed('Y_face_ad')) exit
            call do1_ad(s% grada_face_ad, c% grada_face_ad)
            if (failed('grada_face_ad')) exit
            call do1_ad(s% mlt_vc_ad, c% mlt_vc_ad)
            if (failed('mlt_vc_ad')) exit
            call do1_ad(s% gradL_ad, c% gradL_ad)
            if (failed('gradL_ad')) exit
            call do1_ad(s% scale_height_ad, c% scale_height_ad)
            if (failed('scale_height_ad')) exit
            call do1_ad(s% Lambda_ad, c% Lambda_ad)
            if (failed('Lambda_ad')) exit
            call do1_ad(s% mlt_D_ad, c% mlt_D_ad)
            if (failed('mlt_D_ad')) exit
            call do1_ad(s% mlt_Gamma_ad, c% mlt_Gamma_ad)
            if (failed('mlt_Gamma_ad')) exit
            
            call do1_ad(s% PII_ad, c% PII_ad)
            if (failed('PII_ad')) exit
            call do1_ad(s% Chi_ad, c% Chi_ad)
            if (failed('Chi_ad')) exit
            call do1_ad(s% Eq_ad, c% Eq_ad)
            if (failed('Eq_ad')) exit
            call do1_ad(s% COUPL_ad, c% COUPL_ad)
            if (failed('COUPL_ad')) exit
            call do1_ad(s% Lr_ad, c% Lr_ad)
            if (failed('Lr_ad')) exit
            call do1_ad(s% Lc_ad, c% Lc_ad)
            if (failed('Lc_ad')) exit
            call do1_ad(s% Lt_ad, c% Lt_ad)
            if (failed('Lt_ad')) exit

            call do1(s% gradT, c% gradT)
            if (failed('gradT')) exit
            call do1(s% gradr, c% gradr)
            if (failed('gradr')) exit

            call do1(s% grad_density, c% grad_density)
            if (failed('grad_density')) exit
            call do1(s% grad_temperature, c% grad_temperature)
            if (failed('grad_temperature')) exit
            call do1(s% gradL, c% gradL)
            if (failed('gradL')) exit
            call do1(s% gradL_composition_term, c% gradL_composition_term)
            if (failed('gradL_composition_term')) exit

            call do1_integer( &
               s% dominant_iso_for_thermohaline, c% dominant_iso_for_thermohaline)
            if (failed('dominant_iso_for_thermohaline')) exit

            call do1(s% sig, c% sig)
            if (failed('sig')) exit
            call do1(s% sig_raw, c% sig_raw)
            if (failed('sig_raw')) exit

            call do1(s% brunt_N2, c% brunt_N2)
            if (failed('brunt_N2')) exit
            call do1(s% brunt_N2_composition_term, c% brunt_N2_composition_term)
            if (failed('brunt_N2_composition_term')) exit
            call do1(s% brunt_B, c% brunt_B)
            if (failed('brunt_B')) exit
            call do1(s% unsmoothed_brunt_B, c% unsmoothed_brunt_B)
            if (failed('unsmoothed_brunt_B')) exit
            call do1(s% smoothed_brunt_B, c% smoothed_brunt_B)
            if (failed('smoothed_brunt_B')) exit
            
            call do1(s% RTI_du_diffusion_kick, c% RTI_du_diffusion_kick)
            if (failed('RTI_du_diffusion_kick')) exit

            call do1_ad(s% u_face_ad, c% u_face_ad)
            if (failed('u_face_ad')) exit
            call do1(s% u_face_start, c% u_face_start)
            if (failed('u_face_start')) exit
            call do1(s% u_face_val, c% u_face_val)
            if (failed('u_face_val')) exit
            call do1(s% d_uface_domega, c% d_uface_domega)
            if (failed('d_uface_domega')) exit

            call do1_ad(s% P_face_ad, c% P_face_ad)
            if (failed('P_face_ad')) exit
            call do1(s% P_face_start, c% P_face_start)
            if (failed('P_face_start')) exit
            call do1(s% abs_du_div_cs, c% abs_du_div_cs)
            if (failed('abs_du_div_cs')) exit
            call do1(s% abs_du_plus_cs, c% abs_du_plus_cs)
            if (failed('abs_du_plus_cs')) exit

            call do1(s% dPdr_dRhodr_info, c% dPdr_dRhodr_info)
            if (failed('dPdr_dRhodr_info')) exit
            call do1(s% dPdr_info, c% dPdr_info)
            if (failed('dPdr_info')) exit
            call do1(s% dRhodr_info, c% dRhodr_info)
            if (failed('dRhodr_info')) exit

            call do1(s% source_plus_alpha_RTI, c% source_plus_alpha_RTI)
            if (failed('source_plus_alpha_RTI')) exit
            call do1(s% source_minus_alpha_RTI, c% source_minus_alpha_RTI)
            if (failed('source_minus_alpha_RTI')) exit
            call do1(s% eta_RTI, c% eta_RTI)
            if (failed('eta_RTI')) exit
            call do1(s% etamid_RTI, c% etamid_RTI)
            if (failed('etamid_RTI')) exit
            call do1(s% boost_for_eta_RTI, c% boost_for_eta_RTI)
            if (failed('boost_for_eta_RTI')) exit

            call do1(s% sig_RTI, c% sig_RTI)
            if (failed('sig_RTI')) exit
            call do1(s% sigmid_RTI, c% sigmid_RTI)
            if (failed('sigmid_RTI')) exit

            call do1(s% L_nuc_burn, c% L_nuc_burn)
            if (failed('L_nuc_burn')) exit

            call do2(s% xa_start, c% xa_start, species, 'xa_start')
            if (failed('xa_start')) exit
            call do2(s% xa_sub_xa_start, c% xa_sub_xa_start, species, 'xa_sub_xa_start')
            if (failed('xa_sub_xa_start')) exit

            call do1(s% lnd_start, c% lnd_start)
            if (failed('lnd_start')) exit
            call do1(s% lnPgas_start, c% lnPgas_start)
            if (failed('lnPgas_start')) exit
            call do1(s% lnPeos_start, c% lnPeos_start)
            if (failed('lnPeos_start')) exit
            call do1(s% Peos_start, c% Peos_start)
            if (failed('Peos_start')) exit
            call do1(s% lnT_start, c% lnT_start)
            if (failed('lnT_start')) exit
            call do1(s% energy_start, c% energy_start)
            if (failed('energy_start')) exit
            call do1(s% egas_start, c% egas_start)
            if (failed('egas_start')) exit
            call do1(s% erad_start, c% erad_start)
            if (failed('erad_start')) exit
            call do1(s% Pgas_start, c% Pgas_start)
            if (failed('Pgas_start')) exit
            call do1(s% Prad_start, c% Prad_start)
            if (failed('Prad_start')) exit
            call do1(s% lnR_start, c% lnR_start)
            if (failed('lnR_start')) exit
            call do1(s% v_start, c% v_start)
            if (failed('v_start')) exit
            call do1(s% u_start, c% u_start)
            if (failed('u_start')) exit
            call do1(s% L_start, c% L_start)
            if (failed('L_start')) exit
            call do1(s% r_start, c% r_start)
            if (failed('r_start')) exit
            call do1(s% rmid_start, c% rmid_start)
            if (failed('rmid_start')) exit
            call do1(s% omega_start, c% omega_start)
            if (failed('omega_start')) exit
            call do1(s% ye_start, c% ye_start)
            if (failed('ye_start')) exit
            call do1(s% opacity_start, c% opacity_start)
            if (failed('opacity_start')) exit
            call do1(s% csound_start, c% csound_start)
            if (failed('csound_start')) exit
            call do1(s% alpha_RTI_start, c% alpha_RTI_start)
            if (failed('alpha_RTI_start')) exit

            call do1(s% j_rot_start, c% j_rot_start)
            if (failed('j_rot_start')) exit
            call do1(s% i_rot_start, c% i_rot_start)
            if (failed('i_rot_start')) exit
            call do1(s% eps_nuc_start, c% eps_nuc_start)
            if (failed('eps_nuc_start')) exit
            call do1(s% non_nuc_neu_start, c% non_nuc_neu_start)
            if (failed('non_nuc_neu_start')) exit
            call do2(s% dxdt_nuc_start, c% dxdt_nuc_start, species, null_str)
            if (failed('dxdt_nuc_start')) exit
            call do1(s% grada_start, c% grada_start)
            if (failed('grada_start')) exit
            call do1(s% chiT_start, c% chiT_start)
            if (failed('chiT_start')) exit
            call do1(s% chiRho_start, c% chiRho_start)
            if (failed('chiRho_start')) exit
            call do1(s% cp_start, c% cp_start)
            if (failed('cp_start')) exit
            call do1(s% Cv_start, c% Cv_start)
            if (failed('Cv_start')) exit
            call do1(s% dE_dRho_start, c% dE_dRho_start)
            if (failed('dE_dRho_start')) exit
            call do1(s% gam_start, c% gam_start)
            if (failed('gam_start')) exit
            call do1(s% rho_start, c% rho_start)
            if (failed('rho_start')) exit
            call do1(s% lnS_start, c% lnS_start)
            if (failed('lnS_start')) exit
            call do1(s% T_start, c% T_start)
            if (failed('T_start')) exit
            call do1(s% zbar_start, c% zbar_start)
            if (failed('zbar_start')) exit
            call do1(s% mu_start, c% mu_start)
            if (failed('mu_start')) exit

            call do1(s% phase_start, c% phase_start)
            if (failed('phase_start')) exit
            call do1(s% latent_ddlnT_start, c% latent_ddlnT_start)
            if (failed('latent_ddlnT_start')) exit
            call do1(s% latent_ddlnRho_start, c% latent_ddlnRho_start)
            if (failed('latent_ddlnRho_start')) exit

            call do1(s% max_burn_correction, c% max_burn_correction)
            if (failed('max_burn_correction')) exit
            call do1(s% avg_burn_correction, c% avg_burn_correction)
            if (failed('avg_burn_correction')) exit

            call do1(s% burn_avg_epsnuc, c% burn_avg_epsnuc)
            if (failed('burn_avg_epsnuc')) exit
            call do1_integer(s% burn_num_iters, c% burn_num_iters)
            if (failed('burn_num_iters')) exit

            call do1_neq(s% residual_weight1, c% residual_weight1)
            if (failed('residual_weight1')) exit
            if (action == do_remove_from_center .or. action == do_reallocate .or. &
                  (action /= do_check_size .and. action /= do_deallocate)) &
               s% residual_weight(1:nvar,1:nz) => s% residual_weight1(1:nvar*nz)

            call do1_neq(s% correction_weight1, c% correction_weight1)
            if (failed('correction_weight1')) exit
            if (action == do_remove_from_center .or. action == do_reallocate .or. &
                  (action /= do_check_size .and. action /= do_deallocate)) &
               s% correction_weight(1:nvar,1:nz) => s% correction_weight1(1:nvar*nz)

            call do1_neq(s% solver_dx1, c% solver_dx1)
            if (failed('solver_dx1')) exit
            if (action == do_remove_from_center .or. action == do_reallocate .or. &
                  (action /= do_check_size .and. action /= do_deallocate)) &
               s% solver_dx(1:nvar,1:nz) => s% solver_dx1(1:nvar*nz)

            call do1_neq(s% x_scale1, c% x_scale1)
            if (failed('x_scale1')) exit
            if (action == do_remove_from_center .or. action == do_reallocate .or. &
                  (action /= do_check_size .and. action /= do_deallocate)) &
               s% x_scale(1:nvar,1:nz) => s% x_scale1(1:nvar*nz)

            call do1(s% eps_pre_mix, c% eps_pre_mix)
            if (failed('eps_pre_mix')) exit

            call do1(s% max_abs_xa_corr, c% max_abs_xa_corr)
            if (failed('max_abs_xa_corr')) exit

            call do1(s% Hp_face, c% Hp_face); if (failed('Hp_face')) exit

            call do1(s% Y_face, c% Y_face); if (failed('Y_face')) exit
            call do1(s% Y_face_start, c% Y_face_start); if (failed('Y_face_start')) exit
            
            call do1(s% PII, c% PII); if (failed('PII')) exit

            call do1(s% Chi, c% Chi); if (failed('Chi')) exit
            call do1(s% Chi_start, c% Chi_start); if (failed('Chi_start')) exit

            call do1(s% Lr, c% Lr); if (failed('Lr')) exit
            call do1(s% Lc, c% Lc); if (failed('Lc')) exit
            call do1(s% Lc_start, c% Lc_start); if (failed('Lc_start')) exit
            call do1(s% Lt, c% Lt); if (failed('Lt')) exit
            call do1(s% Lt_start, c% Lt_start); if (failed('Lt_start')) exit
            
            call do1(s% Fr, c% Fr); if (failed('Fr')) exit
            call do1(s% Fr_start, c% Fr_start); if (failed('Fr_start')) exit
            call do1(s% Pvsc, c% Pvsc); if (failed('Pvsc')) exit
            call do1(s% Pvsc_start, c% Pvsc_start); if (failed('Pvsc_start')) exit
            call do1(s% Ptrb, c% Ptrb); if (failed('Ptrb')) exit
            call do1(s% Ptrb_start, c% Ptrb_start); if (failed('Ptrb_start')) exit
            call do1(s% Eq, c% Eq); if (failed('Eq')) exit
            call do1(s% SOURCE, c% SOURCE); if (failed('SOURCE')) exit
            call do1(s% DAMP, c% DAMP); if (failed('DAMP')) exit
            call do1(s% DAMPR, c% DAMPR); if (failed('DAMPR')) exit
            call do1(s% COUPL, c% COUPL); if (failed('COUPL')) exit
            call do1(s% COUPL_start, c% COUPL_start); if (failed('COUPL_start')) exit
            call do1(s% RSP_w, c% RSP_w); if (failed('w')) exit
            call do1(s% RSP_w_start, c% RSP_w_start); if (failed('w_start')) exit
            call do1(s% Vol, c% Vol); if (failed('Vol')) exit
            call do1(s% Vol_start, c% Vol_start); if (failed('Vol_start')) exit
            call do1(s% Uq, c% Uq); if (failed('Uq')) exit
            call do1(s% f_Edd, c% f_Edd); if (failed('f_Edd')) exit

            call do1(s% xtra1_array, c% xtra1_array)
            if (failed('xtra1_array')) exit
            call do1(s% xtra2_array, c% xtra2_array)
            if (failed('xtra2_array')) exit
            call do1(s% xtra3_array, c% xtra3_array)
            if (failed('xtra3_array')) exit
            call do1(s% xtra4_array, c% xtra4_array)
            if (failed('xtra4_array')) exit
            call do1(s% xtra5_array, c% xtra5_array)
            if (failed('xtra5_array')) exit
            call do1(s% xtra6_array, c% xtra6_array)
            if (failed('xtra6_array')) exit

            call do1_integer(s% ixtra1_array, c% ixtra1_array)
            if (failed('ixtra1_array')) exit
            call do1_integer(s% ixtra2_array, c% ixtra2_array)
            if (failed('ixtra2_array')) exit
            call do1_integer(s% ixtra3_array, c% ixtra3_array)
            if (failed('ixtra3_array')) exit
            call do1_integer(s% ixtra4_array, c% ixtra4_array)
            if (failed('ixtra4_array')) exit
            call do1_integer(s% ixtra5_array, c% ixtra5_array)
            if (failed('ixtra5_array')) exit
            call do1_integer(s% ixtra6_array, c% ixtra6_array)
            if (failed('ixtra6_array')) exit

            if (action_in /= do_check_size) then
               if (action_in /= do_copy_pointers_and_resize .and. &
                   action_in /= do_reallocate) then
                  call do2D_dim1(s, s% profile_extra, nz, max_num_profile_extras, action, ierr)
               else
                  deallocate(s% profile_extra)
                  allocate(s% profile_extra(nz, max_num_profile_extras), stat=ierr)
               end if
               if (failed('pstar extras')) exit
            end if

            call do2(s% prev_mesh_xh, c% prev_mesh_xh, nvar_hydro, 'prev_mesh_xh')
            if (failed('prev_mesh_xh')) exit
            call do2(s% prev_mesh_xa, c% prev_mesh_xa, species, 'prev_mesh_xa')
            if (failed('prev_mesh_xa')) exit
            call do1(s% prev_mesh_j_rot, c% prev_mesh_j_rot)
            if (failed('prev_mesh_j_rot')) exit
            call do1(s% prev_mesh_omega, c% prev_mesh_omega)
            if (failed('prev_mesh_omega')) exit
            call do1(s% prev_mesh_mlt_vc, c% prev_mesh_mlt_vc)
            if (failed('prev_mesh_mlt_vc')) exit
            call do1(s% prev_mesh_dq, c% prev_mesh_dq)
            if (failed('prev_mesh_dq')) exit
            ! These are needed for time-smoothing of ST mixing
            call do1(s% prev_mesh_D_ST_start, c% prev_mesh_D_ST_start)
            if (failed('prev_mesh_D_ST_start')) exit
            call do1(s% prev_mesh_nu_ST_start, c% prev_mesh_nu_ST_start)
            if (failed('prev_mesh_nu_ST_start')) exit

            if (s% fill_arrays_with_NaNs) s% need_to_setvars = .true.
            return
         end do
         ierr = -1


         contains


         subroutine do1_ad(ptr, other)
            type(auto_diff_real_star_order1), dimension(:), pointer :: ptr, other
            type(auto_diff_real_star_order1), dimension(:), pointer :: tmp
            if (action == do_fill_arrays_with_NaNs) then
               call fill_ad_with_NaNs(ptr,1,-1)
            else if (action == do_copy_pointers_and_resize) then
               ptr => other
               if (nz <= size(ptr,dim=1)) then
                  if (s% fill_arrays_with_NaNs) call fill_ad_with_NaNs(ptr,1,-1)
                  return
               end if
               deallocate(ptr)
               allocate(ptr(sz_new), stat=ierr)
               if (s% fill_arrays_with_NaNs) call fill_ad_with_NaNs(ptr,1,-1)
               if (s% zero_when_allocate) call fill_ad_with_zeros(ptr,1,-1)
            else
               if (action == do_reallocate) then
                  if (nz <= size(ptr,dim=1)) return
               end if
               call do1D_ad(s, ptr, sz_new, action, ierr)
               if (action == do_allocate) then
                  if (s% fill_arrays_with_NaNs) call fill_ad_with_NaNs(ptr,1,-1)
                  if (s% zero_when_allocate) call fill_ad_with_zeros(ptr,1,-1)
               end if
            end if
         end subroutine do1_ad


         subroutine do1(ptr, other)
            real(dp), dimension(:), pointer :: ptr, other
            real(dp), dimension(:), pointer :: tmp
            if (action == do_fill_arrays_with_NaNs) then
               call fill_with_NaNs(ptr)
            else if (action == do_copy_pointers_and_resize) then
               ptr => other
               if (.not. associated(ptr)) then
                  call mesa_error(__FILE__,__LINE__,'do1 ptr not associated')
               end if
               if (nz <= size(ptr,dim=1)) then
                  if (s% fill_arrays_with_NaNs) call fill_with_NaNs(ptr)
                  return
               end if
               deallocate(ptr)
               allocate(ptr(sz_new), stat=ierr)
               if (s% fill_arrays_with_NaNs) call fill_with_NaNs(ptr)
               if (s% zero_when_allocate) ptr(:) = 0
            else
               if (action == do_reallocate) then
                  if (nz <= size(ptr,dim=1)) return
               end if
               call do1D(s, ptr, sz_new, action, ierr)
               if (action == do_allocate) then
                  if (s% fill_arrays_with_NaNs) call fill_with_NaNs(ptr)
                  if (s% zero_when_allocate) ptr(:) = 0
               end if
            end if
         end subroutine do1


         subroutine do1_neq(ptr, other)
            real(dp), dimension(:), pointer :: ptr, other
            real(dp), dimension(:), pointer :: tmp
            if (action == do_fill_arrays_with_NaNs) then
               call fill_with_NaNs(ptr)
            else if (action == do_copy_pointers_and_resize) then
               ptr => other
               if (nvar*nz <= size(ptr,dim=1)) then
                  if (s% fill_arrays_with_NaNs) call fill_with_NaNs(ptr)
                  if (s% zero_when_allocate) ptr(:) = 0
                  return
               end if
               deallocate(ptr)
               allocate(ptr(nvar*sz_new), stat=ierr)
               if (s% fill_arrays_with_NaNs) call fill_with_NaNs(ptr)
               if (s% zero_when_allocate) ptr(:) = 0
            else
               if (action == do_reallocate) then
                  if (nvar*nz <= size(ptr,dim=1)) return
               end if
               call do1D(s, ptr, nvar*sz_new, action, ierr)
               if (action == do_allocate) then
                  if (s% fill_arrays_with_NaNs) call fill_with_NaNs(ptr)
                  if (s% zero_when_allocate) ptr(:) = 0
               end if
            end if
         end subroutine do1_neq


         subroutine do1_integer(ptr, other)
            integer, dimension(:), pointer :: ptr, other
            integer, dimension(:), pointer :: tmp
            if (action == do_copy_pointers_and_resize) then
               ptr => other
               if (nz <= size(ptr,dim=1)) return
               deallocate(ptr)
               allocate(ptr(sz_new), stat=ierr)
            else
               if (action == do_reallocate) then
                  if (nz <= size(ptr,dim=1)) return
               end if
               call do1D_integer(s, ptr, sz_new, action, ierr)
            end if
         end subroutine do1_integer


         subroutine do2_integer(ptr, other, sz1)
            integer, dimension(:,:), pointer :: ptr, other
            integer, intent(in) :: sz1
            real(dp), dimension(:,:), pointer :: tmp
            if (action == do_copy_pointers_and_resize) then
               ptr => other
               if (sz1 == size(ptr, dim=1) .and. nz <= size(ptr, dim=2)) return
               deallocate(ptr)
               allocate(ptr(sz1, sz_new), stat=ierr)
            else
               if (action == do_reallocate) then
                  if (sz1 == size(ptr, dim=1) .and. nz <= size(ptr, dim=2)) return
               end if
               call do2D_integer(s, ptr, sz1, sz_new, action, ierr)
            end if
         end subroutine do2_integer


         subroutine do1_logical(ptr, other)
            logical, dimension(:), pointer :: ptr, other
            logical, dimension(:), pointer :: tmp
            if (action == do_copy_pointers_and_resize) then
               ptr => other
               if (nz <= size(ptr,dim=1)) return
               deallocate(ptr)
               allocate(ptr(sz_new), stat=ierr)
            else
               if (action == do_reallocate) then
                  if (nz <= size(ptr,dim=1)) return
               end if
               call do1D_logical(s, ptr, sz_new, action, ierr)
            end if
         end subroutine do1_logical


         subroutine do2(ptr, other, sz1, str)
            real(dp), dimension(:,:), pointer :: ptr, other
            integer, intent(in) :: sz1
            character (len=*), intent(in) :: str
            real(dp), dimension(:,:), pointer :: tmp
            include 'formats'
            if (action == do_fill_arrays_with_NaNs) then
               call fill_with_NaNs_2d(ptr)
            else if (action == do_copy_pointers_and_resize) then
               ptr => other
               if (sz1 == size(ptr, dim=1) .and. nz <= size(ptr, dim=2)) then
                  if (s% fill_arrays_with_NaNs) call fill_with_NaNs_2d(ptr)
                  if (s% zero_when_allocate) ptr(:,:) = 0
                  return
               end if
               deallocate(ptr)
               allocate(ptr(sz1, sz_new), stat=ierr)
               if (s% fill_arrays_with_NaNs) call fill_with_NaNs_2d(ptr)
               if (s% zero_when_allocate) ptr(:,:) = 0
            else
               if (action == do_reallocate) then
                  if (sz1 == size(ptr, dim=1) .and. nz <= size(ptr, dim=2)) return
               end if
               call do2D(s, ptr, sz1, sz_new, action, ierr)
               if (action == do_allocate) then
                  if (s% fill_arrays_with_NaNs) call fill_with_NaNs_2d(ptr)
                  if (s% zero_when_allocate) ptr(:,:) = 0
               end if
            end if
         end subroutine do2


         subroutine do3(ptr, other, sz1, sz2)
            real(dp), dimension(:,:,:), pointer :: ptr, other
            integer, intent(in) :: sz1, sz2
            real(dp), dimension(:,:,:), pointer :: tmp
            if (action == do_fill_arrays_with_NaNs) then
               call fill_with_NaNs_3d(ptr)
            elseif (action == do_copy_pointers_and_resize) then
               ptr => other
               if (sz1 == size(ptr, dim=1) .and. sz2 == size(ptr, dim=2) &
                     .and. nz <= size(ptr, dim=3)) then
                  if (s% fill_arrays_with_NaNs) call fill_with_NaNs_3d(ptr)
                  if (s% zero_when_allocate) ptr(:,:,:) = 0
                  return
               end if
               deallocate(ptr)
               allocate(ptr(sz1, sz2, sz_new), stat=ierr)
               if (s% fill_arrays_with_NaNs) call fill_with_NaNs_3d(ptr)
               if (s% zero_when_allocate) ptr(:,:,:) = 0
            else
               if (action == do_reallocate) then
                   if (sz1 == size(ptr, dim=1) .and. &
                       sz2 == size(ptr, dim=2) .and. &
                       nz <= size(ptr, dim=3)) return
               end if
               call do3D(s, ptr, sz1, sz2, sz_new, action, ierr)
               if (action == do_allocate) then
                  if (s% fill_arrays_with_NaNs) call fill_with_NaNs_3d(ptr)
                  if (s% zero_when_allocate) ptr(:,:,:) = 0
               end if
            end if
         end subroutine do3


         subroutine do2_quad(ptr, other, sz1)
            real(qp), dimension(:,:), pointer :: ptr, other
            integer, intent(in) :: sz1
            real(qp), dimension(:,:), pointer :: tmp
            if (action == do_copy_pointers_and_resize) then
               ptr => other
               if (sz1 == size(ptr, dim=1) .and. nz <= size(ptr, dim=2)) return
               deallocate(ptr)
               allocate(ptr(sz1, sz_new), stat=ierr)
            else
               if (action == do_reallocate) then
                  if (sz1 == size(ptr, dim=1) .and. nz <= size(ptr, dim=2)) return
               end if
               call do2D_quad(s, ptr, sz1, sz_new, action, ierr)
            end if
         end subroutine do2_quad


         logical function failed(str)
            character (len=*), intent(in) :: str
            include 'formats'
            failed = .false.
            if (ierr == 0) return
            write(*,*) 'star_info_arrays failed for ' // trim(str)
            write(*,2) 'action', action
            write(*,2) 'species', species
            write(*,2) 'nz', nz
            failed = .true.
         end function failed


      end subroutine star_info_arrays
         
         
      subroutine fill_ad_with_NaNs(ptr, klo, khi_in)
         type(auto_diff_real_star_order1), dimension(:), pointer :: ptr
         integer, intent(in) :: klo, khi_in
         integer :: k, khi
         if (khi_in == -1) then
            khi = size(ptr,dim=1)
         else
            khi = khi_in
         end if
         do k=klo,khi
            call set_nan(ptr(k)% val)
            call fill_with_NaNs(ptr(k)% d1Array)
         end do
      end subroutine fill_ad_with_NaNs
      
      
      subroutine fill_ad_with_zeros(ptr, klo, khi_in)
         type(auto_diff_real_star_order1), dimension(:), pointer :: ptr
         integer, intent(in) :: klo, khi_in
         integer :: k, khi
         if (khi_in == -1) then
            khi = size(ptr,dim=1)
         else
            khi = khi_in
         end if
         do k=klo,khi
            ptr(k)% val = 0d0
            ptr(k)% d1Array(:) = 0d0
         end do
      end subroutine fill_ad_with_zeros


      subroutine do1D_ad(s, ptr, sz, action, ierr)
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1), dimension(:), pointer :: ptr
         integer, intent(in) :: sz, action
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1), dimension(:), pointer :: ptr2
         integer :: old_sz, j
         include 'formats'
         ierr = 0
         select case(action)
            case (do_deallocate)
               if (associated(ptr)) then
                  deallocate(ptr)
                  nullify(ptr)
               end if
            case (do_allocate)
               allocate(ptr(sz), stat=ierr)
               if (s% fill_arrays_with_NaNs) then
                  call fill_ad_with_NaNs(ptr,1,-1)
               else if (s% zero_when_allocate) then
                  call fill_ad_with_zeros(ptr,1,-1)
               end if
            case (do_check_size)
               if (size(ptr,dim=1) < sz) ierr = -1
            case (do_remove_from_center)
               allocate(ptr2(sz), stat=ierr)
               old_sz = size(ptr,dim=1)
               do j=1,min(old_sz,sz)
                  ptr2(j) = ptr(j)
               end do
               deallocate(ptr)
               if (ierr /= 0) return
               ptr => ptr2
            case (do_reallocate)
               if (associated(ptr)) then
                  if (size(ptr,dim=1) >= sz) return
               else
                  ierr = -1
                  return
               end if
               allocate(ptr2(sz), stat=ierr)
               old_sz = size(ptr,dim=1)
               do j=1,old_sz
                  ptr2(j) = ptr(j)
               end do
               if (s% fill_arrays_with_NaNs) then
                  call fill_ad_with_NaNs(ptr2,old_sz+1,sz)
               else if (s% zero_when_allocate) then
                  call fill_ad_with_zeros(ptr2,old_sz+1,sz)
               end if
               deallocate(ptr)
               if (ierr /= 0) return
               ptr => ptr2
            case (do_fill_arrays_with_NaNs)
               if (associated(ptr)) call fill_ad_with_NaNs(ptr,1,-1)
         end select
      end subroutine do1D_ad


      subroutine do1D(s, ptr, sz, action, ierr)
         type (star_info), pointer :: s
         real(dp), dimension(:), pointer :: ptr
         integer, intent(in) :: sz, action
         integer, intent(out) :: ierr
         real(dp), dimension(:), pointer :: ptr2
         integer :: old_sz, j
         include 'formats'
         ierr = 0
         select case(action)
            case (do_deallocate)
               if (associated(ptr)) then
                  deallocate(ptr)
                  nullify(ptr)
               end if
            case (do_allocate)
               allocate(ptr(sz), stat=ierr)
               if (s% fill_arrays_with_NaNs) then
                  call set_nan(ptr)
               else if (s% zero_when_allocate) then
                  ptr(1:sz) = 0
               end if
            case (do_check_size)
               if (size(ptr,dim=1) < sz) ierr = -1
            case (do_remove_from_center)
               allocate(ptr2(sz), stat=ierr)
               old_sz = size(ptr,dim=1)
               do j=1,min(old_sz,sz)
                  ptr2(j) = ptr(j)
               end do
               deallocate(ptr)
               if (ierr /= 0) return
               ptr => ptr2
            case (do_reallocate)
               if (associated(ptr)) then
                  if (size(ptr,dim=1) >= sz) return
               else
                  ierr = -1
                  return
               end if
               allocate(ptr2(sz), stat=ierr)
               old_sz = size(ptr,dim=1)
               do j=1,old_sz
                  ptr2(j) = ptr(j)
               end do
               if (s% fill_arrays_with_NaNs) then
                  do j=old_sz+1,sz
                     call set_nan(ptr2(j))
                  end do
               else if (s% zero_when_allocate) then
                  do j=old_sz+1,sz
                     ptr2(j) = 0
                  end do
               end if
               deallocate(ptr)
               if (ierr /= 0) return
               ptr => ptr2
            case (do_fill_arrays_with_NaNs)
               if (associated(ptr)) call set_nan(ptr)
         end select
      end subroutine do1D


      subroutine do2D(s, ptr, sz1, sz2, action, ierr)

         type (star_info), pointer :: s
         real(dp), dimension(:,:), pointer:: ptr
         integer, intent(in) :: sz1, sz2, action
         integer, intent(out) :: ierr
         real(dp), dimension(:,:), pointer :: ptr2
         integer :: old_sz2, j, i
         ierr = 0
         select case(action)
            case (do_deallocate)
               if (associated(ptr)) then
                  deallocate(ptr)
                  nullify(ptr)
               end if
            case (do_allocate)
               allocate(ptr(sz1, sz2), stat=ierr)
               if (s% fill_arrays_with_NaNs) then
                  call set_nan(ptr)
               else if (s% zero_when_allocate) then
                  ptr(1:sz1,1:sz2) = 0
               end if
            case (do_check_size)
               if (size(ptr,dim=1) /= sz1) ierr = -1
               if (size(ptr,dim=2) < sz2) ierr = -1
            case (do_remove_from_center)
               allocate(ptr2(sz1, sz2), stat=ierr)
               old_sz2 = size(ptr,dim=2)
               do j=1,min(old_sz2,sz2)
                  do i=1,sz1
                     ptr2(i,j) = ptr(i,j)
                  end do
               end do
               deallocate(ptr)
               if (ierr /= 0) return
               ptr => ptr2
            case (do_reallocate)
               if (associated(ptr)) then
                  if (size(ptr,dim=1) /= sz1) then
                     ierr = -1
                     return
                  end if
                  if (size(ptr,dim=2) >= sz2) return
               else
                  ierr = -1
                  return
               end if
               allocate(ptr2(sz1, sz2), stat=ierr)
               old_sz2 = size(ptr,dim=2)
               do j=1,old_sz2
                  do i=1,sz1
                     ptr2(i,j) = ptr(i,j)
                  end do
               end do
               if (s% fill_arrays_with_NaNs) then
                  do j=old_sz2+1,sz2
                     do i=1,sz1
                        call set_nan(ptr2(i,j))
                     end do
                  end do
               else if (s% zero_when_allocate) then
                  do j=old_sz2+1,sz2
                     do i=1,sz1
                        ptr2(i,j) = 0
                     end do
                  end do
               end if
               deallocate(ptr)
               if (ierr /= 0) return
               ptr => ptr2
            case (do_fill_arrays_with_NaNs)
               if (associated(ptr)) call set_nan(ptr)
         end select
      end subroutine do2D


      subroutine do2D_quad(s, ptr, sz1, sz2, action, ierr)
         type (star_info), pointer :: s
         real(qp), dimension(:,:), pointer:: ptr
         integer, intent(in) :: sz1, sz2, action
         integer, intent(out) :: ierr
         real(qp), dimension(:,:), pointer :: ptr2
         real(qp) :: nan
         integer :: old_sz2, j, i
         ierr = 0
         select case(action)
            case (do_deallocate)
               if (associated(ptr)) then
                  deallocate(ptr)
                  nullify(ptr)
               end if
            case (do_allocate)
               allocate(ptr(sz1, sz2), stat=ierr)
               if (s% zero_when_allocate) ptr = 0
            case (do_check_size)
               if (size(ptr,dim=1) /= sz1) ierr = -1
               if (size(ptr,dim=2) < sz2) ierr = -1
            case (do_remove_from_center)
               allocate(ptr2(sz1, sz2), stat=ierr)
               old_sz2 = size(ptr,dim=2)
               do i=1,sz1
                  do j=1,min(old_sz2,sz2)
                     ptr2(i,j) = ptr(i,j)
                  end do
               end do
               deallocate(ptr)
               if (ierr /= 0) return
               ptr => ptr2
            case (do_reallocate)
               if (associated(ptr)) then
                  if (size(ptr,dim=1) /= sz1) then
                     ierr = -1
                     return
                  end if
                  if (size(ptr,dim=2) >= sz2) return
               else
                  ierr = -1
                  return
               end if
               allocate(ptr2(sz1, sz2), stat=ierr)
               old_sz2 = size(ptr,dim=2)
               do j=1,old_sz2
                  do i=1,sz1
                     ptr2(i,j) = ptr(i,j)
                  end do
               end do
               deallocate(ptr)
               if (ierr /= 0) return
               ptr => ptr2
         end select
      end subroutine do2D_quad


      subroutine do2D_dim1(s, ptr, sz1, sz2, action, ierr)

         type (star_info), pointer :: s
         real(dp), dimension(:,:), pointer:: ptr
         integer, intent(in) :: sz1, sz2, action
         integer, intent(out) :: ierr
         real(dp), dimension(:,:), pointer :: ptr2
         integer :: old_sz1, j,i
         ierr = 0
         select case(action)
            case (do_deallocate)
               if (associated(ptr)) then
                  deallocate(ptr)
                  nullify(ptr)
               end if
            case (do_allocate)
               allocate(ptr(sz1, sz2), stat=ierr)
               if (s% fill_arrays_with_NaNs) then
                  call set_nan(ptr)
               else if (s% zero_when_allocate) then
                  ptr(1:sz1,1:sz2) = 0
               end if
            case (do_remove_from_center)
               allocate(ptr2(sz1, sz2), stat=ierr)
               old_sz1 = size(ptr,dim=1)
               do j=1,sz2
                  do i=1,min(old_sz1,sz1)
                     ptr2(i,j) = ptr(i,j)
                  end do
               end do
               deallocate(ptr)
               if (ierr /= 0) return
               ptr => ptr2
         end select
      end subroutine do2D_dim1


      subroutine do3D(s, ptr, sz1, sz2, sz3, action, ierr)

         type (star_info), pointer :: s
         real(dp), dimension(:,:,:), pointer:: ptr
         integer, intent(in) :: sz1, sz2, sz3, action
         integer, intent(out) :: ierr
         real(dp), dimension(:,:,:), pointer :: ptr2
         integer :: old_sz3, j, i, k
         ierr = 0
         select case(action)
            case (do_deallocate)
               if (associated(ptr)) then
                  deallocate(ptr)
                  nullify(ptr)
               end if
            case (do_allocate)
               allocate(ptr(sz1, sz2, sz3), stat=ierr)
               if (s% fill_arrays_with_NaNs) then
                  call set_nan(ptr)
               else if (s% zero_when_allocate) then
                  ptr(1:sz1,1:sz2,1:sz3) = 0
               end if
            case (do_check_size)
               if (size(ptr,dim=1) /= sz1) ierr = -1
               if (size(ptr,dim=2) /= sz2) ierr = -1
               if (size(ptr,dim=3) < sz3) ierr = -1
            case (do_remove_from_center)
               allocate(ptr2(sz1, sz2, sz3), stat=ierr)
               old_sz3 = size(ptr,dim=3)
               do k=1,min(old_sz3,sz3)
                  do j=1,sz2
                     do i=1,sz1
                        ptr2(i,j,k) = ptr(i,j,k)
                     end do
                  end do
               end do
               deallocate(ptr)
               if (ierr /= 0) return
               ptr => ptr2
            case (do_reallocate)
               if (associated(ptr)) then
                  if (size(ptr,dim=1) /= sz1 .or. size(ptr,dim=2) /= sz2) then
                     ierr = -1
                     return
                  end if
                  if (size(ptr,dim=3) >= sz3) return
               else
                  ierr = -1
                  return
               end if
               allocate(ptr2(sz1, sz2, sz3), stat=ierr)
               old_sz3 = size(ptr,dim=3)
               do k=1,old_sz3
                  do j=1,sz2
                     do i=1,sz1
                        ptr2(i,j,k) = ptr(i,j,k)
                     end do
                  end do
               end do
               if (s% fill_arrays_with_NaNs) then
                  do k=old_sz3+1,sz3
                     do j=1,sz2
                        do i=1,sz1
                           call set_nan(ptr2(i,j,k))
                        end do
                     end do
                  end do
               else if (s% zero_when_allocate) then
                  do k=old_sz3+1,sz3
                     do j=1,sz2
                        do i=1,sz1
                           ptr2(i,j,k) = 0
                        end do
                     end do
                  end do
               end if
               deallocate(ptr)
               if (ierr /= 0) return
               ptr => ptr2
            case (do_fill_arrays_with_NaNs)
               if (associated(ptr)) call set_nan(ptr)
         end select
      end subroutine do3D


      subroutine do4D(s, ptr, sz1, sz2, sz3, sz4, action, ierr)

         type (star_info), pointer :: s
         real(dp), dimension(:,:,:,:), pointer:: ptr
         integer, intent(in) :: sz1, sz2, sz3, sz4, action
         integer, intent(out) :: ierr
         real(dp), dimension(:,:,:,:), pointer :: ptr2
         integer :: old_sz4, i, j, k, m
         ierr = 0
         select case(action)
            case (do_deallocate)
               if (associated(ptr)) then
                  deallocate(ptr)
                  nullify(ptr)
               end if
            case (do_allocate)
               allocate(ptr(sz1, sz2, sz3, sz4), stat=ierr)
               if (s% fill_arrays_with_NaNs) then
                  call set_nan(ptr)
               else if (s% zero_when_allocate) then
                  ptr(1:sz1,1:sz2,1:sz3,1:sz4) = 0
               end if
            case (do_check_size)
               if (size(ptr,dim=1) /= sz1) ierr = -1
               if (size(ptr,dim=2) /= sz2) ierr = -1
               if (size(ptr,dim=3) /= sz3) ierr = -1
               if (size(ptr,dim=4) < sz4) ierr = -1
            case (do_remove_from_center)
               allocate(ptr2(sz1, sz2, sz3, sz4), stat=ierr)
               old_sz4 = size(ptr,dim=4)
               do m=1,min(old_sz4,sz4)
                  do k=1,sz3
                     do j=1,sz2
                        do i=1,sz1
                           ptr2(i,j,k,m) = ptr(i,j,k,m)
                        end do
                     end do
                  end do
               end do
               deallocate(ptr)
               if (ierr /= 0) return
               ptr => ptr2
            case (do_reallocate)
               if (associated(ptr)) then
                  if (size(ptr,dim=1) /= sz1 .or. &
                      size(ptr,dim=2) /= sz2 .or. &
                      size(ptr,dim=3) /= sz3) then
                     ierr = -1
                     return
                  end if
                  if (size(ptr,dim=4) >= sz4) return
               else
                  ierr = -1
                  return
               end if
               allocate(ptr2(sz1, sz2, sz3, sz4), stat=ierr)
               old_sz4 = size(ptr,dim=4)
               do m=1,old_sz4
                  do k=1,sz3
                     do j=1,sz2
                        do i=1,sz1
                           ptr2(i,j,k,m) = ptr(i,j,k,m)
                        end do
                     end do
                  end do
               end do
               if (s% fill_arrays_with_NaNs) then
                  do m=old_sz4+1,sz4
                     do k=1,sz3
                        do j=1,sz2
                           do i=1,sz1
                              call set_nan(ptr2(i,j,k,m))
                           end do
                        end do
                     end do
                  end do
               else if (s% zero_when_allocate) then
                  do m=old_sz4+1,sz4
                     do k=1,sz3
                        do j=1,sz2
                           do i=1,sz1
                              ptr2(i,j,k,m) = 0
                           end do
                        end do
                     end do
                  end do
               end if
               deallocate(ptr)
               if (ierr /= 0) return
               ptr => ptr2
            case (do_fill_arrays_with_NaNs)
               call set_nan(ptr)
         end select
      end subroutine do4D


      subroutine do1D_integer(s, ptr, sz, action, ierr)
         type (star_info), pointer :: s
         integer, dimension(:), pointer:: ptr
         integer, intent(in) :: sz, action
         integer, intent(out) :: ierr
         integer, dimension(:), pointer :: ptr2
         integer :: old_sz, j
         ierr = 0
         select case(action)
            case (do_deallocate)
               if (associated(ptr)) then
                  deallocate(ptr)
                  nullify(ptr)
               end if
            case (do_allocate)
               allocate(ptr(sz), stat=ierr)
               if (s% zero_when_allocate) ptr = 0
            case (do_check_size)
               if (size(ptr,dim=1) < sz) ierr = -1
            case (do_remove_from_center)
               allocate(ptr2(sz), stat=ierr)
               old_sz = size(ptr,dim=1)
               do j=1,min(old_sz,sz)
                  ptr2(j) = ptr(j)
               end do
               deallocate(ptr)
               if (ierr /= 0) return
               ptr => ptr2
            case (do_reallocate)
               if (associated(ptr)) then
                  if (size(ptr,dim=1) >= sz) return
               else
                  ierr = -1
                  return
               end if
               allocate(ptr2(sz), stat=ierr)
               old_sz = size(ptr,dim=1)
               do j=1,old_sz
                  ptr2(j) = ptr(j)
               end do
               deallocate(ptr)
               if (ierr /= 0) return
               ptr => ptr2
         end select
      end subroutine do1D_integer


      subroutine do2D_integer(s, ptr, sz1, sz2, action, ierr)
         type (star_info), pointer :: s
         integer, dimension(:, :), pointer:: ptr
         integer, intent(in) :: sz1, sz2, action
         integer, intent(out) :: ierr
         integer, dimension(:,:), pointer :: ptr2
         integer :: old_sz2, j, i
         ierr = 0
         select case(action)
            case (do_deallocate)
               if (associated(ptr)) then
                  deallocate(ptr)
                  nullify(ptr)
               end if
            case (do_allocate)
               allocate(ptr(sz1, sz2), stat=ierr)
               if (s% zero_when_allocate) ptr = 0
            case (do_check_size)
               if (size(ptr,dim=1) /= sz1) ierr = -1
               if (size(ptr,dim=2) < sz2) ierr = -1
            case (do_remove_from_center)
               allocate(ptr2(sz1, sz2), stat=ierr)
               old_sz2 = size(ptr,dim=2)
               do j=1,min(old_sz2,sz2)
                  do i=1,sz1
                     ptr2(i,j) = ptr(i,j)
                  end do
               end do
               deallocate(ptr)
               if (ierr /= 0) return
               ptr => ptr2
            case (do_reallocate)
               if (associated(ptr)) then
                  if (size(ptr,dim=1) /= sz1) then
                     ierr = -1
                     return
                  end if
                  if (size(ptr,dim=2) >= sz2) return
               else
                  ierr = -1
                  return
               end if
               allocate(ptr2(sz1, sz2), stat=ierr)
               old_sz2 = size(ptr,dim=2)
               do j=1,old_sz2
                  do i=1,sz1
                     ptr2(i,j) = ptr(i,j)
                  end do
               end do
               deallocate(ptr)
               if (ierr /= 0) return
               ptr => ptr2
         end select
      end subroutine do2D_integer


      subroutine do1D_logical(s, ptr, sz, action, ierr)
         type (star_info), pointer :: s
         logical, dimension(:), pointer:: ptr
         integer, intent(in) :: sz, action
         integer, intent(out) :: ierr
         logical, dimension(:), pointer :: ptr2
         integer :: old_sz, j
         ierr = 0
         select case(action)
            case (do_deallocate)
               if (associated(ptr)) then
                  deallocate(ptr)
                  nullify(ptr)
               end if
            case (do_allocate)
               allocate(ptr(sz), stat=ierr)
               if (s% zero_when_allocate) ptr = .false.
            case (do_check_size)
               if (size(ptr,dim=1) < sz) ierr = -1
            case (do_remove_from_center)
               allocate(ptr2(sz), stat=ierr)
               old_sz = size(ptr,dim=1)
               do j=1,min(old_sz,sz)
                  ptr2(j) = ptr(j)
               end do
               deallocate(ptr)
               if (ierr /= 0) return
               ptr => ptr2
            case (do_reallocate)
               if (associated(ptr)) then
                  if (size(ptr,dim=1) >= sz) return
               else
                  ierr = -1
                  return
               end if
               allocate(ptr2(sz), stat=ierr)
               old_sz = size(ptr,dim=1)
               do j=1,old_sz
                  ptr2(j) = ptr(j)
               end do
               deallocate(ptr)
               if (ierr /= 0) return
               ptr => ptr2
         end select
      end subroutine do1D_logical


      subroutine set_var_info(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         integer :: i
         
         include 'formats'

         ierr = 0
         i = 0
         
         ! first assign variable numbers
         i = i+1; s% i_lnd = i
         i = i+1; s% i_lnT = i
         i = i+1; s% i_lnR = i
      
         if (.not. s% RSP_flag) then
            i = i+1; s% i_lum = i
         else
            s% i_lum = 0
         end if
      
         if (s% v_flag) then
            i = i+1; s% i_v = i
         else
            s% i_v = 0
         end if

         if (s% u_flag) then
            i = i+1;s% i_u = i
         else
            s% i_u = 0
         end if

         if (s% RTI_flag) then
            i = i+1; s% i_alpha_RTI = i
         else
            s% i_alpha_RTI = 0
         end if

         if (s% RSP_flag) then
            i = i+1; s% i_Et_RSP = i
            i = i+1; s% i_erad_RSP = i
            i = i+1; s% i_Fr_RSP = i
         else
            s% i_Et_RSP = 0
            s% i_erad_RSP = 0
            s% i_Fr_RSP = 0
         end if
         
         if (s% RSP2_flag) then
            i = i+1; s% i_w = i
            i = i+1; s% i_Hp = i
         else 
            s% i_w = 0
            s% i_Hp = 0
         end if

         if (s% w_div_wc_flag) then
            i = i+1; s% i_w_div_wc = i
         else
            s% i_w_div_wc = 0
         end if

         if (s% j_rot_flag) then
            i = i+1; s% i_j_rot = i
         else
            s% i_j_rot = 0
         end if
         
         ! now assign equation numbers
         if (s% i_v /= 0 .or. s% i_u /= 0) then
            s% i_dlnd_dt = s% i_lnd
            s% i_dlnE_dt = s% i_lnT
            s% i_equL = s% i_lum
            s% i_dlnR_dt = s% i_lnR
            s% i_dv_dt = s% i_v
            s% i_du_dt = s% i_u
         else ! HSE is included in dv_dt, so drop dlnR_dt
            s% i_equL = s% i_lnd
            s% i_dv_dt = s% i_lnT
            s% i_dlnE_dt = s% i_lum
            s% i_dlnd_dt = s% i_lnR
            s% i_dlnR_dt = 0
            s% i_du_dt = 0
         end if
      
         s% i_detrb_dt = s% i_w
         s% i_equ_Hp = s% i_Hp
         s% i_dalpha_RTI_dt = s% i_alpha_RTI
         s% i_dEt_RSP_dt = s% i_Et_RSP
         s% i_derad_RSP_dt = s% i_erad_RSP
         s% i_dFr_RSP_dt = s% i_Fr_RSP
         s% i_equ_w_div_wc = s% i_w_div_wc
         s% i_dj_rot_dt = s% i_j_rot

         s% nvar_hydro = i
         s% nvar_total = s% nvar_hydro + s% nvar_chem

         ! Names of the variables
         if (s% i_lnd /= 0) s% nameofvar(s% i_lnd) = 'lnd'
         if (s% i_lnT /= 0) s% nameofvar(s% i_lnT) = 'lnT'
         if (s% i_lnR /= 0) s% nameofvar(s% i_lnR) = 'lnR'
         if (s% i_lum /= 0) s% nameofvar(s% i_lum) = 'L'
         if (s% i_v /= 0) s% nameofvar(s% i_v) = 'v'
         if (s% i_w /= 0) s% nameofvar(s% i_w) = 'w'
         if (s% i_Hp/= 0) s% nameofvar(s% i_Hp) = 'Hp'
         if (s% i_alpha_RTI /= 0) s% nameofvar(s% i_alpha_RTI) = 'alpha_RTI'
         if (s% i_Et_RSP /= 0) s% nameofvar(s% i_Et_RSP) = 'etrb_RSP'
         if (s% i_erad_RSP /= 0) s% nameofvar(s% i_erad_RSP) = 'erad_RSP'
         if (s% i_Fr_RSP /= 0) s% nameofvar(s% i_Fr_RSP) = 'Fr_RSP'
         if (s% i_w_div_wc /= 0) s% nameofvar(s% i_w_div_wc) = 'w_div_wc'
         if (s% i_j_rot /= 0) s% nameofvar(s% i_j_rot) = 'j_rot'
         if (s% i_u /= 0) s% nameofvar(s% i_u) = 'u' 

         ! Names of the equations
         if (s% i_dv_dt /= 0) s% nameofequ(s% i_dv_dt) = 'dv_dt'
         if (s% i_equL /= 0) s% nameofequ(s% i_equL) = 'equL'
         if (s% i_dlnd_dt /= 0) s% nameofequ(s% i_dlnd_dt) = 'dlnd_dt'
         if (s% i_dlnE_dt /= 0) s% nameofequ(s% i_dlnE_dt) = 'dlnE_dt'
         if (s% i_dlnR_dt /= 0) s% nameofequ(s% i_dlnR_dt) = 'dlnR_dt'
         if (s% i_detrb_dt /= 0) s% nameofequ(s% i_detrb_dt) = 'detrb_dt'
         if (s% i_equ_Hp /= 0) s% nameofequ(s% i_equ_Hp) = 'equ_Hp'
         if (s% i_dalpha_RTI_dt /= 0) s% nameofequ(s% i_dalpha_RTI_dt) = 'dalpha_RTI_dt'
         if (s% i_dEt_RSP_dt /= 0) s% nameofequ(s% i_dEt_RSP_dt) = 'dEt_RSP_dt'
         if (s% i_derad_RSP_dt /= 0) s% nameofequ(s% i_derad_RSP_dt) = 'derad_RSP_dt'
         if (s% i_dFr_RSP_dt /= 0) s% nameofequ(s% i_dFr_RSP_dt) = 'dFr_RSP_dt'
         if (s% i_equ_w_div_wc /= 0) s% nameofequ(s% i_equ_w_div_wc) = 'equ_w_div_wc'
         if (s% i_dj_rot_dt /= 0) s% nameofequ(s% i_dj_rot_dt) = 'dj_rot_dt'
         if (s% i_du_dt /= 0) s% nameofequ(s% i_du_dt) = 'du_dt'

         ! chem names are done later by set_chem_names when have set up the net
         

         s% need_to_setvars = .true.

      end subroutine set_var_info


      subroutine set_chem_names(s)
         use chem_def
         type (star_info), pointer :: s
         integer ::  old_size, i, j, cid

         include 'formats'

         if (s% nvar_hydro == 0) return ! not ready to set chem names yet

         old_size = size(s% nameofvar,dim=1)
         if (old_size < s% nvar_total) then
            call realloc(s% nameofvar)
            call realloc(s% nameofequ)
         end if
         do i=1, s% nvar_chem
            cid = s% chem_id(i)
            j = s% nvar_hydro+i
            s% nameofvar(j) = trim(chem_isos% name(cid))
            s% nameofequ(j) = 'equ_' // trim(chem_isos% name(cid))
         end do

         contains

         subroutine realloc(p)
            character (len=name_len), dimension(:), pointer :: p
            character (len=name_len), dimension(:), pointer :: old_p
            integer :: cpy_len, j
            old_p => p
            old_size = size(p,dim=1)
            allocate(p(s% nvar_total))
            cpy_len = min(old_size, s% nvar_total)
            do j=1,cpy_len
               p(j) = old_p(j)
            end do
            deallocate(old_p)
         end subroutine realloc

         subroutine realloc_logical(p)
            logical, dimension(:), pointer :: p
            logical, dimension(:), pointer :: old_p
            integer :: cpy_len, j
            old_p => p
            old_size = size(p,dim=1)
            allocate(p(s% nvar_total))
            cpy_len = min(old_size, s% nvar_total)
            do j=1,cpy_len
               p(j) = old_p(j)
            end do
            deallocate(old_p)
         end subroutine realloc_logical

      end subroutine set_chem_names


      subroutine realloc_work_array(s, crit, ptr, oldsz, newsz, extra, str, ierr)
         type (star_info), pointer :: s
         logical, intent(in) :: crit
         integer, intent(in) :: oldsz, newsz, extra
         real(dp), pointer :: ptr(:)
         character (len=*), intent(in) :: str
         integer, intent(out) :: ierr
         real(dp), pointer :: tmp(:)
         integer :: k
         tmp => ptr
         call work_array(s, .true., crit, ptr, newsz, extra, str, ierr)
         if (ierr /= 0) return
         do k=1,min(oldsz,newsz)
            ptr(k) = tmp(k)
         end do
         call work_array(s, .false., crit, tmp, newsz, extra, str, ierr)
      end subroutine realloc_work_array

      ! if alloc is false, then deallocate ptr
      ! if crit is false, then don't need reentrant allocation
      subroutine work_array(s, alloc, crit, ptr, sz, extra, str, ierr)
         type (star_info), pointer :: s
         logical, intent(in) :: alloc, crit
         integer, intent(in) :: sz, extra
         real(dp), pointer :: ptr(:)
         character (len=*), intent(in) :: str
         integer, intent(out) :: ierr
         ierr = 0
         if (alloc) then
            call do_get_work_array(s, crit, ptr, sz, extra, str, ierr)
         else
            call do_return_work_array(s, crit, ptr, str)
         end if
      end subroutine work_array


      subroutine do_get_work_array(s, crit, ptr, sz, extra, str, ierr)
        type (star_info), pointer :: s
         logical, intent(in) :: crit
         integer, intent(in) :: sz, extra
         real(dp), pointer :: ptr(:)
         character (len=*), intent(in) :: str
         integer, intent(out) :: ierr

         integer :: i
         logical :: okay

         ierr = 0

         if (work_array_debug) then
            allocate(ptr(sz + extra), stat=ierr)
            if (s% fill_arrays_with_NaNs) call set_nan(ptr)
            return
         end if

         okay = .false.

         if (crit) then
!$omp critical (alloc_work_array1)
            num_calls = num_calls + 1 ! not safe, but just for info
            do i = 1, num_work_arrays
               if (get1(i)) then
                  okay = .true.
                  exit
               end if
            end do
!$omp end critical (alloc_work_array1)
         else
            num_calls = num_calls + 1
            do i = 1, num_work_arrays
               if (get1(i)) then
                  okay = .true.
                  exit
               end if
            end do
         end if
         if (okay) return

         allocate(ptr(sz + extra), stat=ierr)
         num_allocs = num_allocs + 1
         if (s% fill_arrays_with_NaNs) then
            call set_nan(ptr)
         else if (s% zero_when_allocate) then
            ptr(:) = 0
         end if

         if (work_array_trace) write(*,*) 'allocate new work array'

         contains

         logical function get1(itry)
            integer, intent(in) :: itry
            real(dp), pointer :: p(:)
            include 'formats'
            if (.not. associated(work_pointers(itry)% p)) then
               get1 = .false.
               return
            end if
            p => work_pointers(itry)% p
            work_pointers(itry)% p => null()
            if (size(p,dim=1) < sz) then
               if (work_array_trace) &
                  write(*,4) 'enlarge work array ' // trim(str), &
                     itry, size(p,dim=1), sz + extra
               deallocate(p)
               allocate(p(sz + extra), stat=ierr)
            else
               if (work_array_trace) &
                  write(*,4) 'use work array ' // trim(str), itry, size(p,dim=1)
            end if
            ptr => p
            get1 = .true.
            if (s% fill_arrays_with_NaNs) then
               call set_nan(ptr)
            else if (s% zero_when_allocate) then
               ptr(:) = 0
            end if
         end function get1

      end subroutine do_get_work_array


      subroutine do_return_work_array(s, crit, ptr, str)
         type (star_info), pointer :: s
         logical, intent(in) :: crit
         real(dp), pointer :: ptr(:)
         character (len=*), intent(in) :: str

         integer :: i
         logical :: okay

         if (.not. associated(ptr)) then
            !write(*,*) 'bogus call on do_return_work_array with nil ptr ' // trim(str)
            !call mesa_error(__FILE__,__LINE__,'do_return_work_array')
            return
         end if

         if (work_array_debug) then
            deallocate(ptr)
            return
         end if

         okay = .false.
         if (crit) then
!$omp critical (alloc_work_array2)
            num_returns = num_returns + 1
            do i=1,num_work_arrays
               if (return1(i)) then
                  okay = .true.
                  exit
               end if
            end do
!$omp end critical (alloc_work_array2)
         else
            num_returns = num_returns + 1
            do i=1,num_work_arrays
               if (return1(i)) then
                  okay = .true.
                  exit
               end if
            end do
         end if
         if (okay) return

         deallocate(ptr)
         num_deallocs = num_deallocs + 1

         contains

         logical function return1(itry)
            integer, intent(in) :: itry
            include 'formats'
            if (associated(work_pointers(itry)% p)) then
               return1 = .false.
               return
            end if
            if (work_array_trace) &
               write(*,3) 'return work array ' // trim(str), itry, size(ptr,dim=1)
            work_pointers(itry)% p => ptr
            ptr => null()
            return1 = .true.
         end function return1

      end subroutine do_return_work_array


      subroutine get_quad_array(s, ptr, sz, extra, str, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: sz, extra
         real(qp), pointer :: ptr(:)
         character (len=*), intent(in) :: str
         integer, intent(out) :: ierr
         call do_get_quad_array(s, .true., ptr, sz, extra, str, ierr)
      end subroutine get_quad_array


      ! okay to use this if sure don't need reentrant allocation
      subroutine non_crit_get_quad_array(s, ptr, sz, extra, str, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: sz, extra
         real(qp), pointer :: ptr(:)
         character (len=*), intent(in) :: str
         integer, intent(out) :: ierr
         call do_get_quad_array(s, .false., ptr, sz, extra, str, ierr)
      end subroutine non_crit_get_quad_array


      subroutine do_get_quad_array(s, crit, ptr, sz, extra, str, ierr)
         type (star_info), pointer :: s
         logical, intent(in) :: crit
         integer, intent(in) :: sz, extra
         real(qp), pointer :: ptr(:)
         character (len=*), intent(in) :: str
         integer, intent(out) :: ierr

         integer :: i
         logical :: okay

         ierr = 0

         if (quad_array_debug) then
            allocate(ptr(sz + extra), stat=ierr)
            return
         end if

         okay = .false.
         if (crit) then
!$omp critical (alloc_quad_work_array1)
            num_calls = num_calls + 1
            do i = 1, num_quad_arrays
               if (get1(i)) then
                  okay = .true.
                  exit
               end if
            end do
!$omp end critical (alloc_quad_work_array1)
         else
            num_calls = num_calls + 1
            do i = 1, num_quad_arrays
               if (get1(i)) then
                  okay = .true.
                  exit
               end if
            end do
         end if
         if (okay) return

         allocate(ptr(sz + extra), stat=ierr)
         num_allocs = num_allocs + 1

         if (quad_array_trace) write(*,*) 'allocate new quad array'

         contains

         logical function get1(itry)
            integer, intent(in) :: itry
            real(qp), pointer :: p(:)
            include 'formats'
            if (.not. associated(quad_pointers(itry)% p)) then
               get1 = .false.
               return
            end if
            p => quad_pointers(itry)% p
            quad_pointers(itry)% p => null()
            if (size(p,dim=1) < sz) then
               if (quad_array_trace) &
                  write(*,4) 'enlarge quad array ' // trim(str), itry, size(p,dim=1), sz + extra
               deallocate(p)
               allocate(p(sz + extra), stat=ierr)
            else
               if (quad_array_trace) &
                  write(*,4) 'use quad array ' // trim(str), itry, size(p,dim=1)
            end if
            ptr => p
            get1 = .true.
         end function get1

      end subroutine do_get_quad_array


      subroutine return_quad_array(s, ptr, str)
         type (star_info), pointer :: s
         real(qp), pointer :: ptr(:)
         character (len=*), intent(in) :: str
         call do_return_quad_array(s, .true., ptr, str)
      end subroutine return_quad_array


      ! okay to use this if sure don't need reentrant allocation
      subroutine non_crit_return_quad_array(s, ptr, str)
         type (star_info), pointer :: s
         real(qp), pointer :: ptr(:)
         character (len=*), intent(in) :: str
         if (.not. associated(ptr)) return
         call do_return_quad_array(s, .false., ptr, str)
      end subroutine non_crit_return_quad_array


      subroutine do_return_quad_array(s, crit, ptr, str)
         type (star_info), pointer :: s
         logical, intent(in) :: crit
         real(qp), pointer :: ptr(:)
         character (len=*), intent(in) :: str

         integer :: i
         logical :: okay

         if (.not. associated(ptr)) return

         if (quad_array_debug) then
            deallocate(ptr)
            return
         end if

         okay = .false.
         if (crit) then
!$omp critical (alloc_quad_work_array2)
            num_returns = num_returns + 1
            do i=1,num_quad_arrays
               if (return1(i)) then
                  okay = .true.
                  exit
               end if
            end do
!$omp end critical (alloc_quad_work_array2)
         else
            num_returns = num_returns + 1
            do i=1,num_quad_arrays
               if (return1(i)) then
                  okay = .true.
                  exit
               end if
            end do
         end if
         if (okay) return

         deallocate(ptr)
         num_deallocs = num_deallocs + 1

         contains

         logical function return1(itry)
            integer, intent(in) :: itry
            include 'formats'
            if (associated(quad_pointers(itry)% p)) then
               return1 = .false.
               return
            end if
            if (quad_array_trace) &
               write(*,3) 'return quad array ' // trim(str), itry, size(ptr,dim=1)
            quad_pointers(itry)% p => ptr
            ptr => null()
            return1 = .true.
         end function return1

      end subroutine do_return_quad_array


      subroutine realloc_integer_work_array(s, ptr, oldsz, newsz, extra, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: oldsz, newsz, extra
         integer, pointer :: ptr(:)
         integer, intent(out) :: ierr
         integer, pointer :: itmp(:)
         integer :: k
         itmp => ptr
         call get_integer_work_array(s, ptr, newsz, extra, ierr)
         if (ierr /= 0) return
         do k=1,min(oldsz,newsz)
            ptr(k) = itmp(k)
         end do
         call return_integer_work_array(s, itmp)
      end subroutine realloc_integer_work_array


      subroutine get_integer_work_array(s, ptr, sz, extra, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: sz, extra
         integer, pointer :: ptr(:)
         integer, intent(out) :: ierr

         integer :: i
         logical :: okay

         ierr = 0

         if (work_array_debug) then
            allocate(ptr(sz + extra), stat=ierr)
            return
         end if

         okay = .false.
!$omp critical (alloc_integer_work_array1)
         num_calls = num_calls + 1
         do i=1,num_int_work_arrays
            if (get1(i)) then
               okay = .true.
               exit
            end if
         end do
!$omp end critical (alloc_integer_work_array1)
         if (okay) return

         allocate(ptr(sz + extra), stat=ierr)
         num_allocs = num_allocs + 1

         if (work_array_trace) &
            write(*,*) 'allocate new integer work array'

         contains

         logical function get1(itry)
            integer, intent(in) :: itry
            integer, pointer :: p(:)
            include 'formats'
            if (.not. associated(int_work_pointers(i)% p)) then
               get1 = .false.
               return
            end if
            p => int_work_pointers(i)% p
            int_work_pointers(i)% p => null()
            if (size(p,dim=1) < sz) then
               if (work_array_trace) &
                  write(*,3) 'enlarge integer work array', size(p,dim=1), sz + extra
               deallocate(p)
               allocate(p(sz + extra), stat=ierr)
            end if
            ptr => p
            get1 = .true.
         end function get1

      end subroutine get_integer_work_array


      subroutine return_integer_work_array(s, ptr)
         type (star_info), pointer :: s
         integer, pointer :: ptr(:)

         integer :: i
         logical :: okay

         if (.not. associated(ptr)) return

         if (work_array_debug) then
            deallocate(ptr)
            return
         end if

         okay = .false.
!$omp critical (alloc_integer_work_array2)
         num_returns = num_returns + 1
         do i=1,num_int_work_arrays
            if (return1(i)) then
               okay = .true.
               exit
            end if
         end do
!$omp end critical (alloc_integer_work_array2)
         if (okay) return

         deallocate(ptr)
         num_deallocs = num_deallocs + 1

         contains

         logical function return1(itry)
            integer, intent(in) :: itry
            integer, pointer :: p(:)
            if (associated(int_work_pointers(itry)% p)) then
               return1 = .false.
               return
            end if
            int_work_pointers(itry)% p => ptr
            ptr => null()
            return1 = .true.
         end function return1

      end subroutine return_integer_work_array


      subroutine get_logical_work_array(s, ptr, sz, extra, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: sz, extra
         logical, pointer :: ptr(:)
         integer, intent(out) :: ierr

         integer :: i
         logical :: okay

         ierr = 0

         if (work_array_debug) then
            allocate(ptr(sz + extra), stat=ierr)
            return
         end if

         okay = .false.
!$omp critical (alloc_logical_work_array1)
         num_calls = num_calls + 1
         do i=1,num_logical_work_arrays
            if (get1(i)) then
               okay = .true.
               exit
            end if
         end do
!$omp end critical (alloc_logical_work_array1)
         if (okay) return

         allocate(ptr(sz + extra), stat=ierr)
         num_allocs = num_allocs + 1

         if (work_array_trace) &
            write(*,*) 'allocate new logical work array'

         contains

         logical function get1(itry)
            integer, intent(in) :: itry
            logical, pointer :: p(:)
            include 'formats'
            if (.not. associated(logical_work_pointers(itry)% p)) then
               get1 = .false.
               return
            end if
            p => logical_work_pointers(itry)% p
            logical_work_pointers(itry)% p => null()
            if (size(p,dim=1) < sz) then
               if (work_array_trace) &
                  write(*,3) 'enlarge logical work array', size(p,dim=1), sz + extra
               deallocate(p)
               allocate(p(sz + extra), stat=ierr)
            end if
            ptr => p
            get1 = .true.
         end function get1

      end subroutine get_logical_work_array


      subroutine return_logical_work_array(s, ptr)
         type (star_info), pointer :: s
         logical, pointer :: ptr(:)

         integer :: i
         logical :: okay

         if (.not. associated(ptr)) return

         if (work_array_debug) then
            deallocate(ptr)
            return
         end if

         okay = .false.
!$omp critical (alloc_logical_work_array2)
         num_returns = num_returns + 1
         do i=1,num_logical_work_arrays
            if (return1(i)) then
               okay = .true.
               exit
            end if
         end do
!$omp end critical (alloc_logical_work_array2)
         if (okay) return

         deallocate(ptr)
         num_deallocs = num_deallocs + 1

         contains

         logical function return1(itry)
            integer, intent(in) :: itry
            logical, pointer :: p(:)
            if (associated(logical_work_pointers(itry)% p)) then
               return1 = .false.
               return
            end if
            logical_work_pointers(itry)% p => ptr
            ptr => null()
            return1 = .true.
         end function return1

      end subroutine return_logical_work_array

      
      subroutine shutdown_alloc ()

         integer :: i

         call free_work_arrays()

      end subroutine shutdown_alloc

      
      subroutine free_work_arrays ()

         integer :: i

         do i=1,num_work_arrays
            if (associated(work_pointers(i)%p)) then
               deallocate(work_pointers(i)%p)
               nullify(work_pointers(i)%p)
               num_deallocs = num_deallocs + 1
            endif
         enddo
         do i=1,num_int_work_arrays
            if (associated(int_work_pointers(i)%p)) then
               deallocate(int_work_pointers(i)%p)
               nullify(int_work_pointers(i)%p)
               num_deallocs = num_deallocs + 1
            endif
         enddo
         do i=1,num_logical_work_arrays
            if (associated(logical_work_pointers(i)%p)) then
               deallocate(logical_work_pointers(i)%p)
               nullify(logical_work_pointers(i)%p)
               num_deallocs = num_deallocs + 1
            endif
         enddo

      end subroutine free_work_arrays

      subroutine size_work_arrays
         integer :: sz, num, i

         include 'formats'
         sz = 0; num = 0
         do i=1,num_work_arrays
            sz = sz + get_size(i)
         end do
         do i=1,num_int_work_arrays
            sz = sz + get_size_i(i)
         end do
         do i=1,num_logical_work_arrays
            sz = sz + get_size_l(i)
         end do

         write(*,'(a,5x,99i8)') &
            'work_arrays: num sz calls returns diff', &
            num, sz, num_calls, num_returns, num_calls-num_returns, &
            num_allocs, num_deallocs, num_allocs-num_deallocs
         write(*,'(A)')

         contains

         integer function get_size(i)
            integer, intent(in) :: i
            if (associated(work_pointers(i)% p)) then
               get_size = size(work_pointers(i)% p,dim=1)
               num = num + 1
            else
               get_size = 0
            end if
         end function get_size

         integer function get_size_i(i)
            integer, intent(in) :: i
            if (associated(int_work_pointers(i)% p)) then
               get_size_i = size(int_work_pointers(i)% p,dim=1)
               num = num + 1
            else
               get_size_i = 0
            end if
         end function get_size_i

         integer function get_size_l(i)
            integer, intent(in) :: i
            if (associated(logical_work_pointers(i)% p)) then
               get_size_l = size(logical_work_pointers(i)% p,dim=1)
               num = num + 1
            else
               get_size_l = 0
            end if
         end function get_size_l

      end subroutine size_work_arrays
      
      ! Cleans array used by history.f90, cant think of better place?
      subroutine dealloc_history(s)
         use utils_lib, only: integer_dict_free
         type(star_info), pointer :: s
      
         if (associated(s% history_values)) then
            deallocate(s% history_values)
            nullify(s% history_values)
         end if
         if (associated(s% history_names)) then
            deallocate(s% history_names)
            nullify(s% history_names)
         end if
         if (associated(s% history_value_is_integer)) then
            deallocate(s% history_value_is_integer)
            nullify(s% history_value_is_integer)
         end if
         if (associated(s% history_names_dict)) then
            call integer_dict_free(s% history_names_dict)
            nullify(s% history_names_dict)
         end if

      end subroutine dealloc_history


      subroutine get_work_array(s, ptr, sz, extra, str, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: sz, extra
         real(dp), pointer :: ptr(:)
         character (len=*), intent(in) :: str
         integer, intent(out) :: ierr
         call do_get_work_array(s, .true., ptr, sz, extra, str, ierr)
      end subroutine get_work_array


      subroutine return_work_array(s, ptr, str)
         type (star_info), pointer :: s
         real(dp), pointer :: ptr(:)
         character (len=*), intent(in) :: str
         call do_return_work_array(s, .true., ptr, str)
      end subroutine return_work_array


      ! okay to use this if sure don't need reentrant allocation
      subroutine non_crit_return_work_array(s, ptr, str)
         type (star_info), pointer :: s
         real(dp), pointer :: ptr(:)
         character (len=*), intent(in) :: str
         call do_return_work_array(s, .false., ptr, str)
      end subroutine non_crit_return_work_array




      end module alloc


