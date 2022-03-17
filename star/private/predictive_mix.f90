! ***********************************************************************
!
!   Copyright (C) 2017-2019  Rich Townsend & The MESA Team
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

module predictive_mix

  ! Uses

  use const_def
  use star_private_def
  use chem_def
  use chem_lib
  use num_lib

  ! No implicit typing

  implicit none

  ! Access specifers

  private

  public :: add_predictive_mixing

  ! Procedures

contains

  subroutine add_predictive_mixing (s, ierr)

    type(star_info), pointer :: s
    integer, intent(out)     :: ierr

    logical, parameter :: DEBUG = .FALSE.

    integer  :: i
    integer  :: j
    logical  :: is_core_zone
    logical  :: is_surf_zone
    logical  :: match_zone_type
    logical  :: match_zone_loc
    logical  :: match_bdy_loc
    integer  :: k

    logical :: mix_mask(s%nz)

    include 'formats'

    ! Initialize

    ierr = 0

    if (DEBUG) then
       write(*, *) 'add_predictive_mixing; model, n_conv_bdy=', &
            s%model_number, s%num_conv_boundaries
    end if

    ! Loop over convective boundaries, from center to surface

    mix_mask = .FALSE.

    conv_bdy_loop : do i = 1, s%num_conv_boundaries

       ! Skip this boundary if it's at the surface, since we don't
       ! predictively mix there

       if (s%conv_bdy_loc(i) == 1) then
          if (DEBUG) then
             write(*,*) 'skip since s%conv_bdy_loc(i) == 1', i
          endif
          cycle conv_bdy_loop
       end if

       ! Loop over predictive mixing criteria

       criteria_loop : do j = 1, NUM_PREDICTIVE_PARAM_SETS

          if (.NOT. s% ctrl% predictive_mix(j)) cycle criteria_loop

          ! Check if the criteria match the current boundary

          select case (s% ctrl% predictive_zone_type(j))
          case ('burn_H')
             match_zone_type = s%burn_h_conv_region(i)
          case ('burn_He')
             match_zone_type = s%burn_he_conv_region(i)
          case ('burn_Z')
             match_zone_type = s%burn_z_conv_region(i)
          case ('nonburn')
             match_zone_type = .NOT. ( &
                  s%burn_h_conv_region(i) .OR. &
                  s%burn_he_conv_region(i) .OR. &
                  s%burn_z_conv_region(i) )              
          case ('any')
             match_zone_type = .TRUE.
          case default
             write(*,*) 'Invalid predictive_zone_type: j, s% ctrl% predictive_zone_type(j)=', j, s% ctrl% predictive_zone_type(j)
             ierr = -1
             return
          end select

          is_core_zone = (i == 1 .AND. s%R_center == 0d0 .AND. s%top_conv_bdy(i))

          if (s%top_conv_bdy(i)) then
             is_surf_zone = s%conv_bdy_loc(i) == 1
          else
             is_surf_zone = s%conv_bdy_loc(i+1) == 1
          endif
                
          select case (s% ctrl% predictive_zone_loc(j))
          case ('core')
             match_zone_loc = is_core_zone
          case ('shell')
             match_zone_loc = .NOT. (is_core_zone .OR. is_surf_zone)
          case ('surf')
             match_zone_loc = is_surf_zone
          case ('any')
             match_zone_loc = .TRUE.
          case default
             write(*,*) 'Invalid predictive_zone_loc: j, s% ctrl% predictive_zone_loc(j)=', j, s% ctrl% predictive_zone_loc(j)
             ierr = -1
             return
          end select

          select case (s% ctrl% predictive_bdy_loc(j))
          case ('bottom')
             match_bdy_loc = .NOT. s%top_conv_bdy(i)
          case ('top')
             match_bdy_loc = s%top_conv_bdy(i)
          case ('any')
             match_bdy_loc = .TRUE.
          case default
             write(*,*) 'Invalid predictive_bdy_loc: j, s% ctrl% predictive_bdy_loc(j)=', j, s% ctrl% predictive_bdy_loc(j)
             ierr = -1
             return
          end select

          if (.NOT. (match_zone_type .AND. match_zone_loc .AND. match_bdy_loc)) cycle criteria_loop

          if (s%conv_bdy_q(i) < s% ctrl% predictive_bdy_q_min(j) .OR. &
              s%conv_bdy_q(i) > s% ctrl% predictive_bdy_q_max(j)) cycle criteria_loop
          
          if (DEBUG) then
             write(*,*) 'Predictive mixing at convective boundary: i, j=', i, j
             write(*,*) '  s% ctrl% predictive_zone_type=', TRIM(s% ctrl% predictive_zone_type(j))
             write(*,*) '  s% ctrl% predictive_zone_loc=', TRIM(s% ctrl% predictive_zone_loc(j))
             write(*,*) '  s% ctrl% predictive_bdy_loc=', TRIM(s% ctrl% predictive_bdy_loc(j))
          endif

          ! Perform the predictive mixing for this boundary

          if (s% ctrl% do_conv_premix) then
             call mesa_error(__FILE__,__LINE__,'Predictive mixing and convective premixing cannot be enabled at the same time')
             stop
          end if

          call do_predictive_mixing(s, i, j, ierr, mix_mask)
          if (ierr /= 0) return

          ! Finish (we apply at most a single predictive mix to each boundary)

          exit criteria_loop

       end do criteria_loop

    end do conv_bdy_loop

    ! Perform a sanity check on D_mix

    check_loop : do k = 1, s%nz

       if (is_bad_num(s%D_mix(k))) then

          if (s% ctrl% stop_for_bad_nums) call mesa_error(__FILE__,__LINE__,'add_predictive_mixing')

       end if

    end do check_loop

    ! Finish

    return

  end subroutine add_predictive_mixing

  !*****

  subroutine do_predictive_mixing (s, i, j, ierr, mix_mask)

    type(star_info), pointer :: s
    integer, intent(in)      :: i
    integer, intent(in)      :: j
    integer, intent(out)     :: ierr
    logical, intent(inout)   :: mix_mask(:)

    logical, parameter :: DEBUG = .FALSE.
    logical, parameter :: DUMP_PREDICTIONS = .FALSE.

    real(dp)       :: superad_thresh 
    real(dp)       :: ingest_factor
    integer        :: iso_id
    integer        :: iso_r
    integer        :: iso_i
    integer        :: k_bot_cz
    integer        :: k_top_cz  
    integer        :: k_bot_ez
    integer        :: k_top_ez
    integer        :: k_bot_mz
    integer        :: k_top_mz
    real(dp)       :: xa_cz(s%species)
    real(dp)       :: xa_cz_burn(s%species)
    real(dp)       :: xa_ez(s%species)
    real(dp)       :: xa_ez_burn(s%species)
    real(dp)       :: xa_mz(s%species)
    real(dp)       :: xa_mz_burn(s%species)
    logical        :: outward
    logical        :: ledoux_extension
    integer        :: k_a
    integer        :: k_b
    real(dp)       :: D(s%nz)
    real(dp)       :: vc(s%nz)
    real(dp)       :: gradr(s%nz)
    real(dp)       :: grada(s%nz)
    real(dp)       :: m_ingest
    real(dp)       :: m_ingest_limit
    real(dp)       :: superad_min
    integer        :: dk
    integer        :: k
    real(dp)       :: rho
    real(dp)       :: cdc
    real(dp)       :: dg0
    real(dp)       :: dg1
    character(256) :: filename
    integer        :: unit

    ! Perform predictive mixing at the i'th convective boundary,
    ! using the j'th set of predictive mixing parameters.  This
    ! follows the scheme described in the MESA IV instrument paper

    ierr = 0

    ! Extract parameters

    superad_thresh = MAX(s% ctrl% predictive_superad_thresh(j), 0._dp)

    if (s% ctrl% predictive_avoid_reversal(j) /= '') then
       iso_id = chem_get_iso_id(s% ctrl% predictive_avoid_reversal(j))
       if (iso_id == nuclide_not_found) then
          write(*,*) 'Invalid isotope name in predictive_avoid_reversal'
          ierr = -1
          return
       end if
       iso_r = s%net_iso(iso_id)
    else
       iso_r = 0
    endif

    if (s% ctrl% predictive_limit_ingestion(j) /= '') then
       iso_id = chem_get_iso_id(s% ctrl% predictive_limit_ingestion(j))
       if (iso_id == nuclide_not_found) then
          write(*,*) 'Invalid isotope name in predictive_limit_ingestion'
          ierr = -1
          return
       end if
       iso_i = s%net_iso(iso_id)
       ingest_factor = s% ctrl% predictive_ingestion_factor(j)
    else
       iso_i = 0
    endif

    ! Determine cell indices spanning the initial convection zone,
    ! including the cells which contain the actual convective
    ! boundaries

    if (s%top_conv_bdy(i)) then

       if (i == 1) then
          k_bot_cz = s%nz
       else
          if (s%top_conv_bdy(i-1)) then
             write(*,*) 'Double top boundary in predictive_mix; i=', i
             ierr = -1
             return
          end if
          k_bot_cz = s%conv_bdy_loc(i-1) - 1
       endif
       
       k_top_cz = s%conv_bdy_loc(i)

    else

       k_bot_cz = s%conv_bdy_loc(i) - 1

       if (i == s%num_conv_boundaries) then
          k_top_cz = 1
       else
          if (.NOT. s%top_conv_bdy(i+1)) then
             write(*,*) 'Double bottom boundary in predictive_mix; i=', i
             ierr = -1
             return
          endif
          k_top_cz = s%conv_bdy_loc(i+1)
       endif

    end if

    if (DEBUG) then
       if (k_bot_cz < s%nz) then
          write(*,*) 'Predictive mixing: i, j, q_top, q_bot:', i, j, s%q(k_top_cz), s%q(k_bot_cz+1)
       else
          write(*,*) 'Predictive mixing: i, j, q_top, q_bot:', i, j, s%q(k_top_cz), 0._dp
       endif
    end if

    ! Determine average abundances of the initial convection zone

    call eval_abundances(s, k_bot_cz, k_top_cz, xa_cz, xa_cz_burn)
    
    ! Decide whether we are starting in the "Ledoux extension" phase,
    ! where the boundary moves to where it would be if the
    ! Schwarzschild (rather than Ledoux) criterion had been used in
    ! mix_info. During Ledoux extension, some of the predictive mixing
    ! checks are disabled

    ledoux_extension = s% ctrl% use_ledoux_criterion

    ! Initialize the extended-zone data

    k_bot_ez = k_bot_cz
    k_top_ez = k_top_cz

    call eval_abundances(s, k_bot_ez, k_top_ez, xa_ez, xa_ez_burn)
    
    ! Begin the predictive mixing search, expanding the extent of the
    ! mixed zone until one of a number of criteria are met
 
    outward = s%top_conv_bdy(i)

    k_bot_mz = k_bot_cz
    k_top_mz = k_top_cz

    search_loop : do

       ! Grow the mixed zone by one cell

       if (outward) then
          k_top_mz = k_top_mz - 1
       else
          k_bot_mz = k_bot_mz + 1
       endif

       ! Evaluate average abundance in the mixed zone

       call eval_abundances(s, k_bot_mz, k_top_mz, xa_mz, xa_mz_burn)

       ! Evaluate mixing coefficients D and vc, together with grad's,
       ! at faces inside the mixed zone

       call eval_mixing_coeffs(s, k_bot_mz, k_top_mz, xa_mz_burn, &
                               k_a, k_b, D, vc, grada, gradr, ierr)
       if (ierr /= 0) return

       ! Update the Ledoux extension flag

       if (ledoux_extension) then

          ! Check if we've reached the Schwarzschild boundary

          if (.NOT. ALL(s%gradr(k_a:k_b) > s%grada_face(k_a:k_b))) then

             ledoux_extension = .FALSE.

          else

             ! Otherwise, update the extended-zone data

             k_bot_ez = k_bot_mz
             k_top_ez = k_top_mz

             call eval_abundances(s, k_bot_ez, k_top_ez, xa_ez, xa_ez_burn)

          endif

       endif

       ! Perform checks that only occur after the Ledoux extension
       ! phase has completed

       if (.NOT. ledoux_extension) then

          ! Check whether the predictive mixing will lead to a
          ! reversal in the abundance evolution of isotope iso_r due
          ! to nuclear burning; if so, finish the search.
       
          if (iso_r /= 0) then

             if (SIGN(1._dp, xa_mz_burn(iso_r)-xa_ez(iso_r)) /= SIGN(1._dp, xa_ez_burn(iso_r)-xa_ez(iso_r))) then
                if (DEBUG) then
                   write(*,*) 'Exiting predictive search due to abundance reversal'
                end if
                exit search_loop
             endif

          end if

          ! Check whether the predictive mixing will cause the
          ! ingestion rate for isotope iso_i to exceed the limit
       
          if (iso_i /= 0) then

             ! Calculate the mass ingested

             m_ingest = xa_mz(iso_i)*SUM(s%dm(k_top_mz:k_bot_mz)) - &
                        xa_ez(iso_i)*SUM(s%dm(k_top_ez:k_bot_ez))

             ! Calculate the limiting ingestion mass

             if (outward) then
                call eval_ingest_limit(s, grada, gradr, ingest_factor, k_a, m_ingest_limit)
             else
                call eval_ingest_limit(s, grada, gradr, ingest_factor, k_b, m_ingest_limit)
             endif

             ! If the mass ingested exceeds the limit, finish the search

             if (m_ingest > m_ingest_limit) then
                if (DEBUG) then
                   write(*,*) 'Exiting predictive search due to ingestion limit exceeded'
                end if
                exit search_loop
             end if

          end if

       end if

       ! See if the growing boundary face of the mixed region is
       ! non-convective; this signals that we've gone too far, and so
       ! finish the search

       if ((      outward .AND. gradr(k_a) < grada(k_a)) .OR. &
           (.NOT. outward .AND. gradr(k_b) < grada(k_b))) then
          if (DEBUG) then
             write(*,*) 'Exiting predictive search due to non-convective growing boundary'
          endif
          exit search_loop
       endif

       ! See if any faces (apart from the growing boundary face) have
       ! a ratio gradr/grada below the superadiabaticity threshold;
       ! this signals we're too close to splitting, and so finish the
       ! search

       if (outward) then
          superad_min = MINVAL(gradr(k_a+1:k_b)/grada(k_a+1:k_b)) - 1._dp
       else
          superad_min = MINVAL(gradr(k_a:k_b-1)/grada(k_a:k_b-1)) - 1._dp
       endif

       if (superad_min <= superad_thresh) then
          if (DEBUG) then
             write(*,*) 'Exiting predictive search due to convection-zone split'
          endif
          exit search_loop
       endif

       ! See if the mixed region has reached the center or surface

       if ((      outward .AND. k_top_mz == 1) .OR. &
           (.NOT. outward .AND. k_bot_mz == s%nz-1)) then
          exit search_loop
       endif

    end do search_loop

    ! If necessary, dump out the mixing prediction for this boundary

    if (DUMP_PREDICTIONS) then
       write(filename, 100) i, j, s%model_number
100    format('pred-mix.i',I0,'.j',I0,'.n',I0)
       open(NEWUNIT=unit, FILE=TRIM(filename), STATUS='REPLACE')
       dump_loop : do k = k_a, k_b
          write(unit, *) k, s%q(k), D(k), gradr(k), grada(k)
       end do dump_loop
       close(unit)
       print *,'Writing prediction data to file:',TRIM(filename)
    end if
   
    ! Back off the mixing by one zone

    if (outward) then
       k_top_mz = k_top_mz + 1
    else
       k_bot_mz = k_bot_mz - 1
    endif

    ! Check the mix mask

    if (ANY(mix_mask(k_top_mz:k_bot_mz))) then
       do k = 1, s%nz
          print *,k,mix_mask(k),k <= k_bot_mz .AND. k >= k_top_mz
       end do
       call mesa_error(__FILE__,__LINE__,'Double predictive')
    else
       mix_mask(k_top_mz:k_bot_mz) = .FALSE.
    endif

    ! Return now if no additional mixing should occur

    if (outward .AND. k_top_mz == k_top_cz) then
       if (DEBUG) then
          write(*,*) 'No predictive mixing at top of zone; boundary i=', i
       endif
       return
    elseif (.NOT. outward .AND. k_bot_mz == k_bot_cz) then
       if (DEBUG) then
          write(*,*) 'No predictive mixing at bottom of zone; boundary i=', i
       endif
       return
    endif

    ! Re-calculate mixed abundances and mixing coefficients, since we
    ! backed off by one zone

    call eval_abundances(s, k_bot_mz, k_top_mz, xa_mz, xa_mz_burn)

    call eval_mixing_coeffs(s, k_bot_mz, k_top_mz, xa_mz_burn, &
                            k_a, k_b, D, vc, grada, gradr, ierr)
    if (ierr /= 0) then
       if (DEBUG) write(*,*) 'Non-zero return from eval_mixing_coeffs in do_predictive_mixing/predictive_mix'
       return
    endif

    if (outward) then
       superad_min = MINVAL(gradr(k_a+1:k_b)/grada(k_a+1:k_b)) - 1._dp
    else
       superad_min = MINVAL(gradr(k_a:k_b-1)/grada(k_a:k_b-1)) - 1._dp
    endif

    ! Update the model with the new D and vc

    if (outward) then
       k_b = k_a
       k_a = k_top_cz
       dk = -1
    else
       k_a = k_bot_cz + 1
       dk = 1
    endif

    face_loop : do k = k_a, k_b, dk

       ! See if we would overwrite an adjacent convection zone; if so,
       ! stop

       if (s%mixing_type(k) == convective_mixing) then
          k_b = k - dk
          exit face_loop
       endif

       if (k > 1) then
          rho = (s%dq(k-1)*s%rho(k) + s%dq(k)*s%rho(k-1))/ &
                (s%dq(k-1) + s%dq(k))
       else
          rho = s%rho(k)
       endif
       
       cdc = (pi4*s%r(k)*s%r(k)*rho)*(pi4*s%r(k)*s%r(k)*rho)*D(k) ! gm^2/sec

       s%cdc(k) = cdc
       s%D_mix(k) = D(k)
       s%conv_vel(k) = vc(k)
       s%mixing_type(k) = convective_mixing

    end do face_loop

    ! Store the mass-fractional location of the new convective
    ! boundary, and reset the locations for the old convective
    ! boundary (cf. set_cz_boundary_info)

    if (outward) then

       s%cz_bdy_dq(k_top_cz) = 0._dp

       dg0 = s%grada_face(k_b-1) - s%gradr(k_b-1)
       dg1 = grada(k_b) - gradr(k_b)

       if (dg0*dg1 < 0) then
          s%cz_bdy_dq(k_top_mz) = find0(0._dp, dg0, s%dq(k_top_mz), dg1)
          if (s%cz_bdy_dq(k_top_mz) < 0._dp .OR. s%cz_bdy_dq(k_top_mz) > s%dq(k_top_mz)) then
             write(2,*) 'bad cz_bdy_dq in do_predictive_mixing', k_top_mz, s%cz_bdy_dq(k_top_mz), s%dq(k_top_mz)
             ierr = -1
             return
          end if
       endif

    else

       s%cz_bdy_dq(k_bot_cz) = 0._dp

       dg0 = grada(k_b) - gradr(k_b)
       dg1 = s%grada_face(k_b+1) - s%gradr(k_b+1)
          
       if (dg0*dg1 < 0) then
          s%cz_bdy_dq(k_bot_mz) = find0(0._dp, dg0, s%dq(k_bot_mz), dg1)
          if (s%cz_bdy_dq(k_bot_mz) < 0._dp .or. s%cz_bdy_dq(k_bot_mz) > s%dq(k_bot_mz)) then
             write(2,*) 'bad cz_bdy_dq in do_predictive_mixing', k_bot_mz, s%cz_bdy_dq(k_bot_mz), s%dq(k_bot_mz)
             ierr = -1
             return
          end if
       endif

    end if

    if (DEBUG) then
       write(*,*) 'Predictive mixing: i, k_a, k_b, q_a, q_b, superad_min=', i, k_a, k_b, s%q(k_a), s%q(k_b), &
            superad_min
    endif

    ! Finish

    return

  end subroutine do_predictive_mixing

  !****

  subroutine eval_abundances (s, k_bot, k_top, xa, xa_burn)

    type(star_info), pointer :: s
    integer, intent(in)      :: k_bot
    integer, intent(in)      :: k_top
    real(dp), intent(out)    :: xa(:)
    real(dp), intent(out)    :: xa_burn(:)

    real(dp) :: mc
    integer  :: l

    ! Evaluate average abundances over cells k_top:k_bot, without and
    ! with a burning correction

    mc = SUM(s%dm(k_top:k_bot))

    do l = 1, s%species

       xa(l) = SUM(s%dm(k_top:k_bot)*s%xa(l,k_top:k_bot))/mc

       xa_burn(l) = xa(l) + SUM(s%dm(k_top:k_bot)*s%dxdt_nuc(l,k_top:k_bot)*s%dt)/mc

    end do

    ! Apply limiting

    xa = MAX(MIN(xa, 1._dp), 0._dp)
    xa = xa/SUM(xa)

    xa_burn = MAX(MIN(xa_burn, 1._dp), 0._dp)
    xa_burn = xa_burn/SUM(xa_burn)

    ! Finish

    return

  end subroutine eval_abundances
    
  !****
  
  subroutine eval_mixing_coeffs (s, k_bot_mz, k_top_mz, xa_mx, k_a, k_b, D, vc, grada, gradr, ierr)

    use eos_def
    use micro
    use turb_info, only: do1_mlt_2

    type(star_info), pointer :: s
    integer, intent(in)      :: k_bot_mz
    integer, intent(in)      :: k_top_mz
    real(dp), intent(in)     :: xa_mx(:)
    integer, intent(out)     :: k_a
    integer, intent(out)     :: k_b
    real(dp), intent(out)    :: D(:)
    real(dp), intent(out)    :: vc(:)
    real(dp), intent(out)    :: grada(:)
    real(dp), intent(out)    :: gradr(:)
    integer, intent(out)     :: ierr

    logical, parameter :: DEBUG = .FALSE.

    real(dp) :: xh
    real(dp) :: xhe
    real(dp) :: z
    real(dp) :: abar
    real(dp) :: zbar, z53bar
    real(dp) :: z2bar
    real(dp) :: ye
    real(dp) :: mass_correction
    real(dp) :: sumx
    integer  :: h1
    real(dp) :: lnd_save(s%nz)
    real(dp) :: Cp_save(s%nz)
    real(dp) :: Cv_save(s%nz)
    real(dp) :: gamma1_save(s%nz)
    real(dp) :: grada_save(s%nz)
    real(dp) :: chiRho_save(s%nz)
    real(dp) :: chiT_save(s%nz)
    real(dp) :: lnfree_e_save(s%nz)
    real(dp) :: d_eos_dlnd_save(num_eos_basic_results,s%nz)
    real(dp) :: d_eos_dlnT_save(num_eos_basic_results,s%nz)
    real(dp) :: xa_save(s%species,s%nz)
    real(dp) :: zbar_save(s%nz)
    integer  :: k
    real(dp) :: w
    real(dp) :: rho_face_save(s%nz)
    integer  :: op_err
    logical  :: make_gradr_sticky_in_solver_iters

    ! Evaluate abundance data

    call basic_composition_info(s%species, s%chem_id, xa_mx, &
         xh, xhe, z, abar, zbar, z2bar, z53bar, ye, mass_correction, sumx)

    h1 = s%net_iso(ih1)

    ! Update EOS data for cells spanning the mixed region

!$OMP PARALLEL DO PRIVATE(k,op_err) SCHEDULE(dynamic,2)
    update_cell_loop_eos : do k = k_top_mz, k_bot_mz

       op_err = 0

       lnd_save(k) = s%lnd(k)

       Cp_save(k) = s%Cp(k)
       Cv_save(k) = s%Cv(k)

       gamma1_save(k) = s%gamma1(k)
       grada_save(k) = s%grada(k)

       chiRho_save(k) = s%chiRho(k)
       chiT_save(k) = s%chiT(k)

       lnfree_e_save(k) = s%lnfree_e(k)

       d_eos_dlnd_save(:,k) = s%d_eos_dlnd(:,k)
       d_eos_dlnT_save(:,k) = s%d_eos_dlnT(:,k)

       call eval_eos(s, k, z, xh, abar, zbar, xa_mx, &
            s%lnd(k), s%Cp(k), s%Cv(k), s%gamma1(k), s%grada(k), &
            s%chiRho(k), s%chiT(k), s%lnfree_e(k), &
            s%d_eos_dlnd(:,k), s%d_eos_dlnT(:,k), op_err)
       if (op_err /= 0) ierr = op_err

       xa_save(:,k) = s%xa(:,k)
       zbar_save(k) = s%zbar(k)

       s%xa(:,k) = xa_mx
       s%zbar(k) = zbar

    enddo update_cell_loop_eos
!$OMP END PARALLEL DO
    if (ierr /= 0) then
       if (s% ctrl% report_ierr) write(*,*) 'Non-zero return from eval_eos in eval_mixing_coeffs/predictive_mix'
       return
    endif


!$OMP PARALLEL DO PRIVATE(k,op_err) SCHEDULE(dynamic,2)
    update_cell_loop_kap : do k = k_top_mz, k_bot_mz
       op_err = 0
       call do_kap_for_cell(s, k, op_err)
       if (op_err /= 0) ierr = op_err
    enddo update_cell_loop_kap
!$OMP END PARALLEL DO
    if (ierr /= 0) then
       if (s% ctrl% report_ierr) write(*,*) 'Non-zero return from do_kap_for_cells in eval_mixing_coeffs/predictive_mix'
       print *,'xa:',xa_mx
       return
    endif



    ! Evaluate mixing coefficients D and vc, together with grad's,
    ! at faces inside the mixed region

    k_a = k_top_mz + 1
    k_b = k_bot_mz

    eval_face_loop: do k = k_a, k_b

       ! Update the face density

       rho_face_save(k) = s%rho_face(k)

       w = s%dq(k-1)/(s%dq(k-1) + s%dq(k))

       s%rho_face(k) = w*exp(s%lnd(k)) + (1._dp-w)*exp(s%lnd(k-1))

       ! Evaluate mixing coefficients etc.
       ! Explicitly set gradL_composition_term to 0 in this call.
       call do1_mlt_2(s, k, make_gradr_sticky_in_solver_iters, op_err, &
            s% alpha_mlt(k), 0._dp)
       if (op_err /= 0) call mesa_error(__FILE__,__LINE__,'non-zero op_err')

       D(k) = s%mlt_D(k)
       vc(k) = s%mlt_vc(k)

       grada(k) = s%grada_face(k)
       gradr(k) = s%gradr(k)

    end do eval_face_loop

    ! Restore EOS and MLT data

!$OMP PARALLEL DO PRIVATE(k,op_err) SCHEDULE(dynamic,2)
    restore_cell_loop : do k = k_top_mz, k_bot_mz

       op_err = 0

       s%lnd(k) = lnd_save(k)

       s%Cp(k) = Cp_save(k)
       s%Cv(k) = Cv_save(k)

       s%gamma1(k) = gamma1_save(k)
       s%grada(k) = grada_save(k)

       s%chiRho(k) = chiRho_save(k)
       s%chiT(k) = chiT_save(k)

       s%lnfree_e(k) = lnfree_e_save(k)

       s%d_eos_dlnd(:,k) = d_eos_dlnd_save(:,k)
       s%d_eos_dlnT(:,k) = d_eos_dlnT_save(:,k)

       s%xa(:,k) = xa_save(:,k)
       s%zbar(k) = zbar_save(k)

       call do_kap_for_cell(s, k, op_err)
       if (op_err /= 0) ierr = op_err

    enddo restore_cell_loop
!$OMP END PARALLEL DO
    if (ierr /= 0) then
       if (s% ctrl% report_ierr) write(*,*) 'Non-zero return from do_kap_for_cells in eval_mixing_coeffs/predictive_mix'
       return
    endif

    restore_face_loop: do k = k_a, k_b

       s%rho_face(k) = rho_face_save(k)
       call do1_mlt_2(s, k, make_gradr_sticky_in_solver_iters, op_err)
       if (op_err /= 0) call mesa_error(__FILE__,__LINE__,'non-zero op_err')

    end do restore_face_loop

    ! Finish

    return

  end subroutine eval_mixing_coeffs

  !****

  subroutine eval_eos (s, k, z, x, abar, zbar, xa, &
       lnRho, Cp, Cv, gamma1, grada, chiRho, chiT, &
       lnfree_e, d_dlnd, d_dlnT, ierr)

    use eos_def
    use eos_support

    type(star_info), pointer :: s
    integer, intent(in)      :: k
    real(dp), intent(in)     :: z
    real(dp), intent(in)     :: x
    real(dp), intent(in)     :: abar
    real(dp), intent(in)     :: zbar
    real(dp), intent(in)     :: xa(:)
    real(dp), intent(out)    :: lnRho
    real(dp), intent(out)    :: Cp
    real(dp), intent(out)    :: Cv
    real(dp), intent(out)    :: gamma1
    real(dp), intent(out)    :: grada
    real(dp), intent(out)    :: chiRho
    real(dp), intent(out)    :: chiT
    real(dp), intent(out)    :: lnfree_e
    real(dp), intent(out)    :: d_dlnd(:)
    real(dp), intent(out)    :: d_dlnT(:)
    integer, intent(out)     :: ierr

    logical, parameter  :: DEBUG = .FALSE.
    real(dp), parameter :: LOGRHO_TOL = 1E-8_dp
    real(dp), parameter :: LOGPGAS_TOL = 1E-8_dp

    real(dp) :: logRho
    real(dp) :: res(num_eos_basic_results)
    real(dp) :: d_dxa(num_eos_d_dxa_results,s%species)

    ! Evaluate EOS data in cell k, assuming the cell's temperature and
    ! pressure are as specified in the model, but with abundances
    ! given by xa and other input abundance parameters

    ! (NEEDS FIXING TO HANDLE CASE WHEN LNPGAS_FLAG = .TRUE.)

    call solve_eos_given_PgasT( &
         s, k, xa, &
         s%lnT(k)/ln10, s%lnPgas(k)/ln10, s%lnd(k)/ln10, LOGRHO_TOL, LOGPGAS_TOL, &
         logRho, res, d_dlnd, d_dlnT, d_dxa, &
       ierr)
    if (ierr /= 0) then
       if (DEBUG) write(*,*) 'Non-zero return from solve_eos_given_PgasT in eval_eos/predictive_mix'
       return
    endif

    lnRho = logRho*ln10

    Cp = res(i_Cp)
    Cv = res(i_Cv)

    gamma1 = res(i_gamma1)
    grada = res(i_grad_ad)

    chiRho = res(i_chiRho)
    chiT = res(i_chiT)

    lnfree_e = res(i_lnfree_e)

    ! Finish
    
    return

  end subroutine eval_eos

  !****

  subroutine eval_ingest_limit (s, grada, gradr, ingest_factor, k, m_ingest_limit)

    type(star_info), pointer :: s
    real(dp), intent(in)     :: grada(:)
    real(dp), intent(in)     :: gradr(:)
    real(dp), intent(in)     :: ingest_factor
    integer, intent(in)      :: k
    real(dp), intent(out)    :: m_ingest_limit

    real(dp) :: alfa
    real(dp) :: beta
    real(dp) :: T_face

    ! Evaluate the mass ingestion limit at face k following Spruit
    ! (2015)

    ! Interpolate T to the face

    if (k == 1) then
       alfa = 1._dp
    else
       alfa = s%dq(k-1)/(s%dq(k-1) + s%dq(k))
    end if
    beta = 1._dp - alfa
    
    T_face = alfa*s%T(k) + beta*s%T(k-1)

    ! Evaluate the limit

    m_ingest_limit = ingest_factor*(amu*s%L(k)/(boltzm*T_face))*(1._dp - grada(k)/gradr(k))*s%dt

    ! Finish

    return

  end subroutine eval_ingest_limit

end module predictive_mix
