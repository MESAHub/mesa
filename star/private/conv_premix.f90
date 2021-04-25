! ***********************************************************************
!
!   Copyright (C) 2017-2018  Rich Townsend & The MESA Team
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

module conv_premix

  ! Uses

  use const_def
  use star_private_def
  use chem_def
  use chem_lib
  use num_lib

  ! No implicit typing

  implicit none

  ! Parameter definitions

  integer, parameter :: BURN_TYPE_NONE = 1
  integer, parameter :: BURN_TYPE_H = 2
  integer, parameter :: BURN_TYPE_HE = 3
  integer, parameter :: BURN_TYPE_Z = 4

  integer, parameter :: FIXED_PT_MODE = 5
  integer, parameter :: FIXED_DT_MODE = 6

  logical, parameter :: TRACE_ALL = .FALSE.

  ! Derived-type definitions

  ! The zone_info type stores information about the extent &
  ! properties of each convective zone

  type zone_info
     integer               :: kc_t = 0                   ! Cell index of top boundary
     integer               :: kc_b = 0                   ! Cell index of bot boundary
     real(dp)              :: vc_t = 0._dp               ! Convective velocity near top boundary
     real(dp)              :: vc_b = 0._dp               ! Convective velocity near bot boundary
     real(dp)              :: vp_t = 0._dp               ! Penetration velocity at top boundary
     real(dp)              :: vp_b = 0._dp               ! Penetration velocity at bot boundary
     real(dp)              :: dt_t = 0._dp               ! Clock for top boundary
     real(dp)              :: dt_b = 0._dp               ! Clock for bottom boundary
     integer               :: burn_type = BURN_TYPE_NONE ! Type of burning in zone
     logical               :: sel_t = .FALSE.            ! Flag to select top boundary for modification
     logical               :: sel_b = .FALSE.            ! Flag to select bot boundary for modification
     logical               :: split_retreat = .FALSE.    ! Flag to indicate the zone has split, or the trailing bdy retreated
     real(dp), allocatable :: avg_xa(:)                  ! Initial average abundances
     real(dp), allocatable :: davg_xa_dt(:)              ! Initial rate of change of average abundances (due to burning)
     integer               :: avoid_inc_iso = 0          ! Isotope id for which increases should be avoided
  end type zone_info

  ! The saved_data type saves cell and face data from star_info, to
  ! enable us to restore the latter to an earlier state. Only data in
  ! cells kc_t:kc_b are used, but the arrays are full-size

  type saved_data
     real(dp), allocatable :: xa(:,:)
     real(dp), allocatable :: lnd(:)
     real(dp), allocatable :: rho(:)
     real(dp), allocatable :: lnPgas(:)
     real(dp), allocatable :: Pgas(:)
     real(dp), allocatable :: gradL_composition_term(:)
     integer, allocatable  :: update_mode(:)
     integer               :: kc_t = 0
     integer               :: kc_b = 0
  end type saved_data

  ! Access specifers

  private

  public :: do_conv_premix

  ! Procedures

contains

  subroutine do_conv_premix (s, ierr)

    use star_utils, only: start_time, update_time

    type(star_info), pointer       :: s
    integer, intent(out)           :: ierr

    logical, parameter :: TRACE_CONV_PREMIX = TRACE_ALL

    integer                      :: update_mode(s%nz)
    type(saved_data)             :: sd
    type(zone_info), allocatable :: zi(:)
    integer                      :: j_it
    real(dp), allocatable        :: dr(:)
    real(dp), allocatable        :: dt_t(:)
    real(dp), allocatable        :: dt_b(:)
    integer                      :: i_bdy
    logical                      :: t_bdy

    integer(8)                   :: time0
    real(dp)                     :: total

    ierr = 0

    if (s% doing_timing) call start_time(s, time0, total)

    ! Initialize the update_mode values. These control how cells will
    ! be updated (i.e., how the thermodynamic state will be
    ! re-evaluated after abundances have changed). For untouched
    ! cells, update_mode is determined by the lnPgas_flag setting

    !if (s%lnPgas_flag) then
    !   update_mode = FIXED_PT_MODE
    !else
       update_mode = FIXED_DT_MODE
    !endif

    ! Allocate the saved data arrays

    call alloc_saved_data_(s, sd)

    ! Initialize the zone info table

    call init_zone_info_(s, zi)

    ! Loop until there are no boundaries left to advance

    j_it = 0

    iter_loop : do

       j_it = j_it + 1

       ! Validate the current zone info table

       call validate_zone_info_(s, zi)

       ! Check for completion

       if (.NOT. (ANY(zi%sel_t .OR. zi%sel_b))) exit iter_loop

       ! Calculate mixing timescales for each boundary (the explicit
       ! allocations of dt_t and dt_b are required to workaround a
       ! gfortran bug, where allocate-on-assign will raise bogus array
       ! bounds violations)

       dr = s%r(zi%kc_t+1) - s%r(zi%kc_b)

       allocate(dt_t(SIZE(zi)))
       allocate(dt_b(SIZE(zi)))

       dt_t = MERGE(dr/zi%vc_t, HUGE(0._dp), MASK=zi%sel_t)
       dt_b = MERGE(dr/zi%vc_b, HUGE(0._dp), MASK=zi%sel_b)

       ! Select the boundary with the shortest timescale

       i_bdy = MINLOC(MIN(dt_t, dt_b), DIM=1)
       t_bdy = dt_t(i_bdy) < dt_b(i_bdy)

       if (TRACE_CONV_PREMIX) then
          write(*,*) 'Selected i_bdy, t_bdy:', i_bdy, t_bdy
       end if

       ! Advance the selected boundary, modifying the model and the
       ! zone info table

       call advance_bdy_(s, update_mode, zi, i_bdy, t_bdy, sd, j_it)

       ! Loop around

       deallocate(dt_t)
       deallocate(dt_b)

    end do iter_loop

    ! Finish
    s% need_to_setvars = .true.

    if (s% doing_timing) call update_time(s, time0, total, s% time_conv_premix)
    
    return

  end subroutine do_conv_premix
  
  !****

  subroutine advance_bdy_ (s, update_mode, zi, i_bdy, t_bdy, sd, j_it)

    type(star_info), pointer                    :: s
    integer, intent(inout)                      :: update_mode(:)
    type(zone_info), allocatable, intent(inout) :: zi(:)
    integer, intent(inout)                      :: i_bdy
    logical, intent(in)                         :: t_bdy
    type(saved_data), intent(inout)             :: sd
    integer, intent(in)                         :: j_it

    logical, parameter :: TRACE_ADVANCE_BDY = TRACE_ALL

    integer        :: j_ad
    character(256) :: filename

    ! Advance the top (bottom) boundary of the zone zi(i_bdy). t_bdy
    ! is .true. for the top boundary, .false. for the bottom. Because
    ! the zone info table can potentially change during the process
    ! (as zones split/merge), so can i_bdy

    if (TRACE_ADVANCE_BDY) then
       write(*,*) 'Advancing zone boundary: i_bdy, t_bdy, j_it=', i_bdy, t_bdy, j_it
    end if

    if (s% conv_premix_dump_snapshots) then

       write(filename, 100) 'cpm-snap.', s%model_number, '.it_', j_it, '.ad_', 0, '.dat'
100    format(A,I5.5,A,I5.5,A,I5.5,A)

       call dump_snapshot_(s, filename)

    end if

    j_ad = 0

    advance_loop : do

       j_ad = j_ad + 1

       ! Check whether we've reached the edge of the grid

       if (t_bdy .AND. zi(i_bdy)%kc_t == 1) then
          zi(i_bdy)%sel_t = .FALSE.
       elseif (.NOT. t_bdy .AND. zi(i_bdy)%kc_b == s%nz) then
          zi(i_bdy)%sel_b = .FALSE.
       endif
       
       ! Check whether the advancing face is still selected

       if ((t_bdy .AND. .NOT. zi(i_bdy)%sel_t) .OR. &
            (.NOT. t_bdy .AND. .NOT. zi(i_bdy)%sel_b)) exit advance_loop

       ! Advance the boundary by one cell

       call advance_bdy_by_one_cell_(s, update_mode, zi, i_bdy, t_bdy, sd)

       ! If necessary, dump a snapshot

       if (s% conv_premix_dump_snapshots .AND. MOD(j_ad, 50) == 0) then

          write(filename, 100) 'cpm-snap.', s%model_number, '.it_', j_it, '.ad_', j_ad, '.dat'

          call dump_snapshot_(s, filename)

       end if

       ! Loop around

    end do advance_loop

    if (s% conv_premix_dump_snapshots) then

       write(filename, 100) 'cpm-snap.', s%model_number, '.it_', j_it, '.ad_', j_ad-1, '.dat'

       call dump_snapshot_(s, filename)

    end if

    if (TRACE_ADVANCE_BDY) then
       write(*,*) 'Done advancing zone boundary'
    end if

    ! Finish

    return

  end subroutine advance_bdy_

  !****

  subroutine advance_bdy_by_one_cell_ (s, update_mode, zi, i_bdy, t_bdy, sd)

    type(star_info), pointer                    :: s
    integer, intent(inout)                      :: update_mode(:)
    type(zone_info), allocatable, intent(inout) :: zi(:)
    integer, intent(inout)                      :: i_bdy
    logical, intent(in)                         :: t_bdy
    type(saved_data), intent(inout)             :: sd

    integer, parameter :: N_SUB_INI = 1
    integer, parameter :: N_SUB_MAX = 1024
    logical, parameter :: TRACE_MIX_CELL = TRACE_ALL
    logical, parameter :: TRACE_MIX_SUBCELL = TRACE_ALL

    real(dp)              :: dr
    real(dp)              :: vp
    real(dp)              :: delta_dt
    real(dp)              :: dt_limit
    integer               :: kc_new
    integer               :: n_sub
    integer               :: kc_t
    integer               :: kc_b
    integer               :: kf_t
    integer               :: kf_b
    real(dp), allocatable :: xa_sub(:,:)
    real(dp)              :: m_sub
    logical               :: has_split
    integer               :: kc_sub
    real(dp)              :: m_zone
    integer               :: l
    real(dp)              :: avg_xa(s%species)
    integer               :: n_nc
    integer               :: kf
    real(dp)              :: davg_xa_dt(s%species)

    ! Advance the top (bottom) boundary of the zone zi(i_bdy) by one
    ! cell

    ! Calculate the clock increment for adding the new cell

    if (t_bdy) then
       kc_t = zi(i_bdy)%kc_t
       dr = s%r(kc_t) - s%r(kc_t+1)
       vp = zi(i_bdy)%vp_t
    else
       kc_b = zi(i_bdy)%kc_b
       if (kc_b < s%nz) then
          dr = s%r(kc_b) - s%r(kc_b+1)
       else
          dr = s%r(kc_b)
       endif
       vp = zi(i_bdy)%vp_b
    endif

    delta_dt = dr/vp

    ! Check whether we have enough time; if not, return

    if (s%conv_premix_time_factor > 0._dp) then
       dt_limit = s%conv_premix_time_factor*s%dt
    else
       dt_limit = HUGE(0._dp)
    endif

    if (t_bdy) then

       if (zi(i_bdy)%dt_t + delta_dt > dt_limit) then
        
          zi(i_bdy)%sel_t = .FALSE.

          if (TRACE_MIX_CELL) then
             write(*,*) '  Outcome: out of time', zi(i_bdy)%dt_t, delta_dt, dt_limit
          end if

          return

       else

          zi(i_bdy)%dt_t = zi(i_bdy)%dt_t + delta_dt

       endif

    else

       if (zi(i_bdy)%dt_b + delta_dt > dt_limit) then

          zi(i_bdy)%sel_b = .FALSE.

          if (TRACE_MIX_CELL) then
             write(*,*) '  Outcome: out of time', zi(i_bdy)%dt_b, delta_dt, dt_limit
          end if

          return

       else

          zi(i_bdy)%dt_b = zi(i_bdy)%dt_b + delta_dt

       endif

    endif

    ! Determine the index of the new cell we're going to add

    if (t_bdy) then
       kc_new = zi(i_bdy)%kc_t - 1
    else
       kc_new = zi(i_bdy)%kc_b + 1
    endif

    ! Save the current model state over all cells we *might* modify

    if (t_bdy) then
       call save_model_(s, update_mode, kc_new, zi(i_bdy)%kc_b, sd)
    else
       call save_model_(s, update_mode, zi(i_bdy)%kc_t, kc_new, sd)
    endif

    ! Loop over subcell refinement levels

    n_sub = N_SUB_INI

    refine_loop : do

       ! Set initial cell and face indices

       kc_t = zi(i_bdy)%kc_t
       kc_b = zi(i_bdy)%kc_b

       kf_t = kc_t
       kf_b = kc_b + 1

       if (TRACE_MIX_CELL) then
          write(*,*) 'Attempting to add cell via n_sub subcells:', n_sub
          write(*,*) '  kc_t   :', kc_t
          write(*,*) '  kc_b   :', kc_b
          write(*,*) '  kc_new :', kc_new
       end if

       ! Set update_mode over the cells being mixed

       if (s%conv_premix_fix_pgas) then
          update_mode(kc_t:kc_b) = FIXED_PT_MODE
       else
          update_mode(kc_t:kc_b) = FIXED_DT_MODE
       endif

       ! Set gradL_composition_term to zero on interior faces plus the
       ! advancing face

       if (t_bdy) then
          s%gradL_composition_term(kf_t:kf_b-1) = 0._dp
       else
          s%gradL_composition_term(kf_t+1:kf_b) = 0._dp
       endif

       ! Divide the new cell into n_sub virtual subcells, and store the
       ! (initially uniform) abundances of these subcells

       if (ALLOCATED(xa_sub)) deallocate(xa_sub)
       allocate(xa_sub(s%species,n_sub))

       do l = 1, s%species
          xa_sub(l,:) = s%xa(l,kc_new)
       end do

       m_sub = s%dm(kc_new)/n_sub

       ! Mix subcells into the zone, one by one

       has_split = .FALSE.

       subcell_loop : do kc_sub = 1, n_sub

          if (TRACE_MIX_SUBCELL) then
             write(*,*) 'In subcell loop:', kc_sub, n_sub
          end if

          ! Determine average abundances across the zone and the 1:k_s
          ! subcells

          m_zone = SUM(s%dm(kc_t:kc_b))
   
          do l = 1, s%species
             avg_xa(l) = (SUM(s%dm(kc_t:kc_b)*s%xa(l,kc_t:kc_b)) + &
                          SUM(m_sub*xa_sub(l,1:kc_sub))) / (m_zone + m_sub*kc_sub)
          end do

          avg_xa = MAX(MIN(avg_xa, 1._dp), 0._dp)
          avg_xa = avg_xa/SUM(avg_xa)

          ! Update abundances using the average

          do l = 1, s%species
             s%xa(l,kc_t:kc_b) = avg_xa(l)
             xa_sub(l,1:kc_sub) = avg_xa(l)
          end do

          call update_model_(s, update_mode, kc_t, kc_b)

          ! Determine how many interior faces are now non-convective

          n_nc = COUNT(s%mlt_mixing_type(kf_t+1:kf_b-1) /= convective_mixing)

          if (n_nc > 1 .AND. n_sub < N_SUB_MAX) then

             ! More than a single one: double the number of subcells
             ! and restart (as long as we're below the subcell limit)

             n_sub = 2*n_sub

             call restore_model_(s, update_mode, sd)

             cycle refine_loop

          elseif (n_nc > 0) then

             ! A single one (or more, if we're over the subcell
             ! limit): move the trailing boundary to its new position
             ! (defined as the closest non-convective face to the
             ! advancing boundary)

             has_split = .TRUE.
             
             if (t_bdy) then

                search_down_loop : do kf = kf_t+1, kf_b-1
                   if (s%mlt_mixing_type(kf) /= convective_mixing) exit search_down_loop
                end do search_down_loop

                kc_b = kf - 1
                kf_b = kc_b + 1

                if (TRACE_MIX_SUBCELL) then
                   write(*,*) '  Moved lower boundary to', kc_b
                end if

             else

                search_up_loop : do kf = kf_b-1, kf_t+1, -1
                   if (s%mlt_mixing_type(kf) /= convective_mixing) exit search_up_loop
                end do search_up_loop

                kc_t = kf
                kf_t = kc_t
                
                if (TRACE_MIX_SUBCELL) then
                   write(*,*) '  Moved upper boundary to', kc_t
                end if

             endif

          end if

       end do subcell_loop

       ! If we reach this point, all went well; exit

       if (TRACE_MIX_SUBCELL) then
          write(*,*) '  Completed mixing subcell'
       end if

       exit refine_loop

    end do refine_loop

    ! Update the abundances in the new cell, and adjust the cell/face
    ! indices

    s%xa(:,kc_new) = avg_xa

    if (s%conv_premix_fix_pgas) then
       update_mode(kc_new) = FIXED_PT_MODE
    else
       update_mode(kc_new) = FIXED_DT_MODE
    endif

    if (t_bdy) then
       kc_t = kc_new
       kf_t = kc_t
    else
       kc_b = kc_new
       kf_b = kc_b + 1
    endif

    call update_model_(s, update_mode, kc_t, kc_b)

    ! Check whether an abundance increase has occurred; if so, revert
    ! back to the starting model and return

    if (zi(i_bdy)%avoid_inc_iso /= 0) then

       ! Evaluate the rate of change of the zone-average abundances

       davg_xa_dt = (avg_xa - zi(i_bdy)%avg_xa)/s%dt

       ! Check whether an increase will occur (davg_xa_dt is positive,
       ! and greater in magnitude than the rate of decrease due to
       ! burning)

       associate (iso => zi(i_bdy)%avoid_inc_iso)

         if (davg_xa_dt(iso) > MAX(-zi(i_bdy)%davg_xa_dt(iso), 0._dp)) then

            call restore_model_(s, update_mode, sd)

            if (t_bdy) then
               zi(i_bdy)%sel_t = .FALSE.
            else
               zi(i_bdy)%sel_b = .FALSE.
            endif

            if (TRACE_MIX_CELL) then
               write(*,*) '  Outcome: abundance reversal for isotope', iso
            end if

            return

         end if

       end associate

    end if

    ! Determine whether the face just inside the advancing boundary
    ! has become/remained convective; if not, revert back to the
    ! starting model and return
       
    if (t_bdy) then

       if (s%mlt_mixing_type(kf_t+1) /= convective_mixing) then

          call restore_model_(s, update_mode, sd)

          zi(i_bdy)%sel_t = .FALSE.

          if (TRACE_MIX_CELL) then
             write(*,*) '  Outcome: non-convective at top'
          end if

          return

       endif

    else
       
       if (s%mlt_mixing_type(kf_b-1) /= convective_mixing) then

          call restore_model_(s, update_mode, sd)

          zi(i_bdy)%sel_b = .FALSE.

          if (TRACE_MIX_CELL) then
             write(*,*) '  Outcome: non-convective at bottom'
          end if

          return

       endif

    end if

    ! Update the zone info table to reflect all the changes we made

    call update_zone_info_(s, i_bdy, t_bdy, has_split, zi)

    ! Finish

    return

  end subroutine advance_bdy_by_one_cell_

  !****

  subroutine init_zone_info_ (s, zi)

    type(star_info), pointer                  :: s
    type(zone_info), allocatable, intent(out) :: zi(:)

    integer :: i

    ! Initialize the zone info table for the whole star

    call create_zone_info_(s, 1, s%nz, zi)

    ! Perform additional initializations

    zone_loop : do i = 1, SIZE(zi)

       ! Set penetration velocities

       call set_penetration_velocities_(s, zi(i))

       ! Set burn types and initial average abundances       

       call set_burn_data_(s, zi(i))
       call set_abund_data_(s, zi(i))

       ! Initialize clocks

       zi(i)%dt_t = 0._dp
       zi(i)%dt_b = 0._dp

       ! Set flags

       zi(i)%sel_t = zi(i)%kc_t /= 1
       zi(i)%sel_b = zi(i)%kc_b /= s%nz

       zi(i)%split_retreat = .FALSE.

    end do zone_loop

    ! Finish

    return

  end subroutine init_zone_info_

  !****

  subroutine create_zone_info_ (s, kc_t, kc_b, zi)

    type(star_info), pointer                  :: s
    integer, intent(in)                       :: kc_t
    integer, intent(in)                       :: kc_b
    type(zone_info), allocatable, intent(out) :: zi(:)

    logical         :: in_conv
    type(zone_info) :: zi_new
    integer         :: kf

    ! Create a zone info table spanning the cells kc_t:kc_b. For the
    ! purposes of convective premixing, a zone is defined as a regions
    ! of two or more cells with (i) convective faces inside the
    ! region, (ii) non-convective faces at the top and bottom
    ! boundaries of the region.  The possible exceptions to these
    ! rules apply to zones which begin (end) at kc_t (kc_b); their
    ! upper (lower) boundary faces are not checked for being
    ! non-convective.

    if (kc_t >= kc_b) then
       write(*,*) 'conv_premix: Invalid cell range in create_zone_info_'
       stop
    endif

    ! Create the empty zone info table

    allocate(zi(0))

    ! Loop through faces, looking for zone boundaries and adding zones
    ! (cf. locate_convection_boundaries in mix_info.f90)

    kf = kc_b

    in_conv = s%mlt_mixing_type(kf) == convective_mixing

    if (in_conv) then

       zi_new%kc_b = kf
       zi_new%vc_b = s%mlt_vc(kf)

    endif

    face_loop : do kf = kc_b-1, kc_t+1, -1

       if (in_conv .AND. s%mlt_mixing_type(kf) /= convective_mixing) then

          ! Transitioning out of a convective zone; complete
          ! definition of new zone info and add to to the table

          zi_new%kc_t = kf
          zi_new%vc_t = s%mlt_vc(kf+1)

          zi = [zi,zi_new]

          in_conv = .FALSE.

       elseif (.NOT. in_conv .AND. s%mlt_mixing_type(kf) == convective_mixing) then

          ! Transitioning into a convective zone; begin definition of
          ! new zone info

          zi_new%kc_b = kf
          zi_new%vc_b = s%mlt_vc(kf)

          in_conv = .TRUE.

       endif

    end do face_loop

    kf = kc_t

    if (in_conv) then

       zi_new%kc_t = kf
       zi_new%vc_t = s%mlt_vc(kf+1)

       zi = [zi,zi_new]

    end if

    ! Finish

    return

  end subroutine create_zone_info_

  !****

  subroutine update_zone_info_ (s, i_bdy, t_bdy, has_split, zi)

    type(star_info), pointer                    :: s
    integer, intent(inout)                      :: i_bdy
    logical, intent(in)                         :: t_bdy
    logical, intent(in)                         :: has_split
    type(zone_info), allocatable, intent(inout) :: zi(:)

    logical, parameter :: TRACE_UPDATE_ZONE = TRACE_ALL

    type(zone_info), allocatable :: zi_new(:)
    integer                      :: i

    ! Update the zone info table to account for a single-cell
    ! advancement in zone i_bdy, direction t_bdy

    ! Advance the zone by one cell

    if (t_bdy) then
       zi(i_bdy)%kc_t = zi(i_bdy)%kc_t - 1
       zi(i_bdy)%vc_t = s%mlt_vc(zi(i_bdy)%kc_t+1)
    else
       zi(i_bdy)%kc_b = zi(i_bdy)%kc_b + 1
       zi(i_bdy)%vc_b = s%mlt_vc(zi(i_bdy)%kc_b)
    endif

    call set_penetration_velocities_(s, zi(i_bdy))

    ! If necessary, truncate (or even delete) the adjacent zone we've
    ! encroached into

    if (t_bdy .AND. i_bdy < SIZE(zi)) then

       if (zi(i_bdy)%kc_t == zi(i_bdy+1)%kc_b) then

          ! Encroached into zone above

          if (zi(i_bdy+1)%kc_b-zi(i_bdy+1)%kc_t > 1) then

             ! Truncate the zone

             zi(i_bdy+1)%kc_b = zi(i_bdy+1)%kc_b - 1
             zi(i_bdy+1)%vc_b = s%mlt_vc(zi(i_bdy+1)%kc_b)

             if (TRACE_UPDATE_ZONE) then
                write(*,*) 'Truncated zone above to', zi(i_bdy+1)%kc_b, zi(i_bdy+1)%vc_b, &
                     s%mlt_mixing_type(zi(i_bdy+1)%kc_b), zi(i_bdy+1)%kc_b-zi(i_bdy+1)%kc_t+1
             endif
             
          else

             ! Delete the zone

             zi = [zi(:i_bdy),zi(i_bdy+2:)]

             if (TRACE_UPDATE_ZONE) then
                write(*,*) 'Deleted zone above'
             endif

          end if

       end if

    elseif (.NOT. t_bdy .AND. i_bdy > 1) then

       if (zi(i_bdy)%kc_b == zi(i_bdy-1)%kc_t) then

          ! Encroached into zone below

          if (zi(i_bdy-1)%kc_b-zi(i_bdy-1)%kc_t > 1) then

             ! Truncate the zone
          
             zi(i_bdy-1)%kc_t = zi(i_bdy-1)%kc_t + 1
             zi(i_bdy-1)%vc_t = s%mlt_vc(zi(i_bdy-1)%kc_t+1)

             if (TRACE_UPDATE_ZONE) then
                write(*,*) 'Truncated zone below to', zi(i_bdy-1)%kc_t
             endif
             
          else

             ! Delete the zone
             
             zi = [zi(:i_bdy-2),zi(i_bdy:)]
             i_bdy = i_bdy - 1

             if (TRACE_UPDATE_ZONE) then
                write(*,*) 'Deleted zone below'
             endif

          end if

       end if

    endif

    ! If the zone has split (or its trailing boundary has retreated),
    ! then replace the zone with as many new zones as necessary

    if (has_split) then

       ! Create new zones to replace the zone

       call create_zone_info_(s, zi(i_bdy)%kc_t, zi(i_bdy)%kc_b, zi_new)

       ! Set up parameters for the new zones

       new_zone_loop : do i = 1, SIZE(zi_new)

          ! Clocks & selection flags

          if (zi_new(i)%kc_t == zi(i_bdy)%kc_t) then
             zi_new(i)%dt_t = zi(i_bdy)%dt_t
             zi_new(i)%sel_t = zi(i_bdy)%sel_t
          else
             if (t_bdy) then
                zi_new(i)%dt_t = zi(i_bdy)%dt_t
                zi_new(i)%sel_t = .FALSE.
             else
                zi_new(i)%dt_t = zi(i_bdy)%dt_b
                zi_new(i)%sel_t = .FALSE.
             endif
          endif
          
          if (zi_new(i)%kc_b == zi(i_bdy)%kc_b) then
             zi_new(i)%dt_b = zi(i_bdy)%dt_b
             zi_new(i)%sel_b = zi(i_bdy)%sel_b
          else
             if (t_bdy) then
                zi_new(i)%dt_b = zi(i_bdy)%dt_t
                zi_new(i)%sel_b = .FALSE.
             else
                zi_new(i)%dt_b = zi(i_bdy)%dt_b
                zi_new(i)%sel_b = .FALSE.
             endif
          endif

          ! Penetration velocities

          call set_penetration_velocities_(s, zi_new(i))

          ! Burn data

          call set_burn_data_(s, zi_new(i))

          ! Initial abundances 

          zi_new(i)%avg_xa = zi(i_bdy)%avg_xa
          zi_new(i)%davg_xa_dt = zi(i_bdy)%davg_xa_dt

       end do new_zone_loop

       ! Replace the zone with the new zones

       zi = [zi(:i_bdy-1),zi_new,zi(i_bdy+1:)]

       if (TRACE_UPDATE_ZONE) then
          write(*,*) 'Replaced zone with new zones', SIZE(zi), SIZE(zi_new)
       end if

    else

       if (TRACE_UPDATE_ZONE) then
          write(*,*) 'Zone did not split, updating indices',&
               COUNT(s%mlt_mixing_type(zi(i_bdy)%kc_t+1:zi(i_bdy)%kc_b) /= convective_mixing)
       end if

    end if

    ! Finish

    return

  end subroutine update_zone_info_

  !****

  subroutine set_penetration_velocities_ (s, zi)

    type(star_info), pointer       :: s
    type(zone_info), intent(inout) :: zi

    real(dp) :: mu_i
    real(dp) :: mu_e
    integer  :: kf
    real(dp) :: z_s

    ! Set the penetration velocities for the zone info. These are
    ! calculated by evaluating the molecular weight contrast between
    ! the boundary cell and the adjacent (radiative) cell, and
    ! calculating a characteristic penetration velocity using the
    ! approach by Castellani (1971)

    ! Bottom boundary

    if (zi%kc_b < s%nz) then

       mu_i = s%mu(zi%kc_b)
       mu_e = s%mu(zi%kc_b+1)

       kf = zi%kc_b

       if (mu_e - mu_i > EPSILON(0._dp)) then
          z_s = MIN(zi%vc_b*zi%vc_b/(2._dp*s%grav(kf)*(1._dp - mu_i/mu_e)), s%mlt_mixing_length(kf))
       else
          z_s = s%mlt_mixing_length(kf)
       endif

       zi%vp_b = zi%vc_b*z_s/s%mlt_mixing_length(kf)

    else

       zi%vp_b = zi%vc_b

    end if

    if (zi%vp_b == 0._dp) print *,'Bottom bdy has zero vp', zi%kc_b, zi%vp_b, zi%vc_b, z_s

    ! Top boundary

    if (zi%kc_t > 1) then

       mu_i = s%mu(zi%kc_t)
       mu_e = s%mu(zi%kc_t-1)

       kf = zi%kc_t + 1

       if (mu_i - mu_e > EPSILON(0._dp)) then
          z_s = MIN(zi%vc_t*zi%vc_t/(2._dp*s%grav(kf)*(1._dp - mu_e/mu_i)), s%mlt_mixing_length(kf))
       else
          z_s = s%mlt_mixing_length(kf)
       endif

       zi%vp_t = zi%vc_t*z_s/s%mlt_mixing_length(kf)

    else

       zi%vp_t = zi%vc_t

    end if

    if (zi%vp_t == 0._dp) print *,'Top bdy has zero vp', zi%kc_t, zi%vp_t, zi%vc_t, z_s

    ! Finish

    return

  end subroutine set_penetration_velocities_

  !****

  subroutine set_burn_data_ (s, zi)

    type(star_info), pointer       :: s
    type(zone_info), intent(inout) :: zi

    real(dp), parameter :: EPS_THRESHOLD = 1._dp ! Threshold eps for a zone to be classified as burning

    integer  :: iso_h1
    integer  :: iso_he4
    integer  :: iso_c12
    real(dp) :: eps_max
    real(dp) :: eps_h_max
    real(dp) :: eps_he_max
    integer  :: kc

    ! Set burning data for the zone info
    
    iso_h1 = s%net_iso(chem_get_iso_id('h1'))
    iso_he4 = s%net_iso(chem_get_iso_id('he4'))
    iso_c12 = s%net_iso(chem_get_iso_id('c12'))

    associate (kc_t => zi%kc_t, &
               kc_b => zi%kc_b)

      ! Find the maximum burning rate, and the H/He burning rate at
      ! the same location

      eps_max = -HUGE(0._dp)
      eps_h_max = -HUGE(0._dp)
      eps_he_max = -HUGE(0._dp)

      cell_loop : do kc = kc_t, kc_b

         if (s%eps_nuc(kc) > eps_max) then

            eps_max = s%eps_nuc(kc)

            eps_h_max = s%eps_nuc_categories(ipp, kc) + &
                        s%eps_nuc_categories(icno, kc)
            eps_he_max = s%eps_nuc_categories(i3alf, kc) + &
                         s%eps_nuc_categories(i_burn_c, kc)
            
         endif

      end do cell_loop

      ! Classify the principal burning type

      if (eps_max > EPS_THRESHOLD) then

         if (eps_h_max > 0.5_dp*eps_max) then

            ! Hydrogen burning zone

            zi%burn_type = BURN_TYPE_H

         elseif (eps_he_max > 0.5_dp*eps_max) then

            ! Helium burning zone

            zi%burn_type = BURN_TYPE_HE

         else

            ! Metal burning zone

            zi%burn_type = BURN_TYPE_Z

         end if

      else

         ! Non-burning zone

         zi%burn_type = BURN_TYPE_NONE

      endif

    end associate

    ! Based on the burning type, decide which isotope (if any)
    ! should be monitored for abundance reversal avoidance

    if (s%conv_premix_avoid_increase) then

       select case (zi%burn_type)
       case (BURN_TYPE_H)
          zi%avoid_inc_iso = iso_h1
       case (BURN_TYPE_HE)
          zi%avoid_inc_iso = iso_he4
       case default
          zi%avoid_inc_iso = 0
       end select

    endif

    ! Finish

    return

  end subroutine set_burn_data_

  !****

  subroutine set_abund_data_ (s, zi)

    type(star_info), pointer       :: s
    type(zone_info), intent(inout) :: zi

    real(dp) :: m_zone
    integer  :: l

    ! Set abundance data for the zone info
    
    allocate(zi%avg_xa(s%species))
    allocate(zi%davg_xa_dt(s%species))

    associate (kc_t => zi%kc_t, &
               kc_b => zi%kc_b, &
               avg_xa => zi%avg_xa, &
               davg_xa_dt => zi%davg_xa_dt)

      ! Calculate average abundances and their rates of change

      m_zone = SUM(s%dm(kc_t:kc_b))

      do l = 1, s%species
         avg_xa(l) = SUM(s%dm(kc_t:kc_b)*s%xa(l,kc_t:kc_b))/m_zone
         davg_xa_dt(l) = SUM(s%dm(kc_t:kc_b)*s%dxdt_nuc(l,kc_t:kc_b))/m_zone
      end do

      avg_xa = MAX(MIN(avg_xa, 1._dp), 0._dp)
      avg_xa = avg_xa/SUM(avg_xa)

    end associate

    ! Finish

    return

  end subroutine set_abund_data_

  !****

  subroutine validate_zone_info_ (s, zi)

    type(star_info), pointer    :: s
    type(zone_info), intent(in) :: zi(:)

    logical, parameter :: TRACE_VALIDATE_INFO = TRACE_ALL

    logical :: valid
    integer :: n_zi
    integer :: i

    integer :: kf

    ! Validate zone info data

    n_zi = SIZE(zi)

    valid = .TRUE.

    zone_info_loop : do i = 1, n_zi

       if (TRACE_VALIDATE_INFO) then
          write(*,*) 'Validating info for zone', i
          call print_zone_info_(s, zi(i))
       end if

       ! (i) Degenerate zone (just one cell)

       if (zi(i)%kc_t == zi(i)%kc_b) then
          write(*,*) 'conv_premix: Degenerate zone'
          valid = .FALSE.
       endif

       ! (iii) Out of bounds zone (cells outside grid)

       if (zi(i)%kc_t < 1) then
          write(*,*) 'conv_premix: Zone top outside grid'
          valid = .FALSE.
       end if

       if (zi(i)%kc_b > s%nz) then
          write(*,*) 'conv_premix: Zone bottom outside grid'
          valid = .FALSE.
       end if

       ! (iv) Overlapping zone (cells inside previous/next zone)

       if (i > 1) then
          if (zi(i)%kc_b >= zi(max(1,i-1))%kc_t) then 
             ! bp: max(1,i-1) to prevent bogus warning from gfortran
             write(*,*) 'conv_premix: Zone bottom inside previous zone'
             valid = .FALSE.
          endif
       endif

       if (i < n_zi) then
          if (zi(i)%kc_t <= zi(i+1)%kc_b) then
             write(*,*) 'conv_premix: Zone top inside next zone'
             valid = .FALSE.
          endif
       endif

       ! (vi) Convective velocities

       if (zi(i)%vc_t <= 0._dp) then
          write(*,*) 'conv_premix: Convective velocity is not greater than zero at zone top'
          valid = .FALSE.
       endif

       if (zi(i)%vc_b <= 0._dp) then
          write(*,*) 'conv_premix: Convective velocity is not greater than zero at zone bottom'
          valid = .FALSE.
       endif

       ! (vii) Mixing type

       if (.NOT. ALL(s%mlt_mixing_type(zi(i)%kc_t+1:zi(i)%kc_b) == convective_mixing)) then
          do kf = zi(i)%kc_t+1, zi(i)%kc_b
             if (s%mlt_mixing_type(kf) /= convective_mixing) then
                write(*,*) 'Mix type at kf=',kf, s%mlt_mixing_type(kf)
             end if
          end do
          write(*,*) 'conv_premix: Zone contains cells with mixing_type /= convective_mixing'
          valid = .FALSE.
       end if

    end do zone_info_loop

    if (.NOT. valid) stop

    ! Finish

    return

  end subroutine validate_zone_info_

  !****

  subroutine print_zone_info_ (s, zi)

    type(star_info), pointer    :: s
    type(zone_info), intent(in) :: zi

    ! Print out zone info

    write(*,*) '  nz            :', s%nz
    write(*,*) '  kc_t          :', zi%kc_t
    write(*,*) '  kc_b          :', zi%kc_b
    write(*,*) '  mass(kf_t)    :', (s%M_center + s%xmstar*s%q(zi%kc_t))/Msun
    if (zi%kc_b < s%nz) then
       write(*,*) '  mass(kf_b)    :', (s%M_center + s%xmstar*s%q(zi%kc_b+1))/Msun
    else
       write(*,*) '  mass(kf_b)    :', (s%M_center)/Msun
    endif
    write(*,*) '  vc_t          :', zi%vc_t
    write(*,*) '  vc_b          :', zi%vc_b
    write(*,*) '  dt_t          :', zi%dt_t
    write(*,*) '  dt_b          :', zi%dt_b
    write(*,*) '  burn_type     :', zi%burn_type
    write(*,*) '  sel_t         :', zi%sel_t
    write(*,*) '  sel_b         :', zi%sel_b
    write(*,*) '  avoid_inc_iso :', zi%avoid_inc_iso
    write(*,*)

    ! Finish

    return

  end subroutine print_zone_info_

  !****

  subroutine update_model_ (s, update_mode, kc_t, kc_b)

    use mlt_info_newer, only: set_mlt_vars_newer
    use micro

    type(star_info), pointer :: s
    integer, intent(in)      :: update_mode(:)
    integer, intent(in)      :: kc_t
    integer, intent(in)      :: kc_b

    logical, parameter :: TRACE_UPDATE_MODEL = .FALSE.

    logical  :: lnPgas_flag
    integer  :: ierr
    integer  :: kf_t
    integer  :: kf_b

    ! Update the model to reflect changes in the abundances across
    ! cells kc_t:kc_b

    if (TRACE_UPDATE_MODEL) then
       write(*,*) '  Updating model'
    end if

    ! Perform the EOS updates, in accordance with the update_mode
    ! values. We make two passes -- the first to handle cells with
    ! update_mode == FIXED_PT_MODE, the second to handle cells with
    ! update_mode == FIXED_DT_MODE. If/when the lnPgas_flag
    ! functionality is retired, this will need to be rewritten

    lnPgas_flag = .FALSE. ! s%lnPgas_flag

    !s%lnPgas_flag = .TRUE.

    !call set_eos_with_mask(s, kc_t, kc_b, update_mode==FIXED_PT_MODE, ierr)
    !if (ierr /= 0) then
    !   write(*,*) 'conv_premix: error from call to set_eos_with_mask'
    !   stop
    !end if

    !s%lnPgas_flag = .FALSE.

    call set_eos_with_mask(s, kc_t, kc_b, update_mode==FIXED_DT_MODE, ierr)
    if (ierr /= 0) then
       write(*,*) 'conv_premix: error from call to set_eos_with_mask'
       stop
    end if

    !s%lnPgas_flag = lnPgas_flag

    ! Update opacities across cells kc_t:kc_b (this also sets rho_face
    ! and related quantities on faces kc_t:kc_b)

    call set_micro_vars(s, kc_t, kc_b, &
         skip_eos=.TRUE., skip_net=.TRUE., skip_neu=.TRUE., skip_kap=.FALSE., ierr=ierr)
    if (ierr /= 0) then
       write(*,*) 'conv_premix: error from call to set_micro_vars'
       stop
    end if

    ! Finally update MLT for interior faces

    kf_t = kc_t
    kf_b = kc_b + 1

    call set_mlt_vars_newer(s, kf_t+1, kf_b-1, ierr)
    if (ierr /= 0) then
       write(*,*) 'conv_premix: failed in call to set_mlt_vars_newer during update_model_'
       stop
    endif

    ! Finish

    return

  end subroutine update_model_

  !****

  subroutine dump_snapshot_ (s, filename)

    use brunt
    use mlt_info_newer, only: set_grads_newer

    type(star_info), pointer :: s
    character(*), intent(in) :: filename

    real(dp), allocatable :: gradL_composition_term(:)
    integer               :: ierr
    integer               :: unit
    integer               :: iso_h1
    integer               :: iso_he4
    integer               :: iso_c12
    integer               :: iso_o16
    integer               :: k

    ! Dump a snapshot of the current model state to file

    ! First, re-calculate gradL_composition_term
    allocate(gradL_composition_term(1:s%nz))

    gradL_composition_term = s%gradL_composition_term(1:s%nz)

    call do_brunt_B(s, 1, s%nz, ierr)
    if (ierr /= 0) then
       write(*,*) 'conv_premix: error from call to do_brunt_B'
       stop
    end if

    call set_grads_newer(s, ierr)
    if (ierr /= 0) then
       write(*,*) 'conv_premix: error from call to set_grads_newer'
       stop
    end if
       
    ! Dump the snapshot

    open(NEWUNIT=unit, FILE=filename, STATUS='REPLACE')

    iso_h1 = s%net_iso(chem_get_iso_id('h1'))
    iso_he4 = s%net_iso(chem_get_iso_id('he4'))
    iso_c12 = s%net_iso(chem_get_iso_id('c12'))
    iso_o16 = s%net_iso(chem_get_iso_id('o16'))

    write(*,*) '  Dumping to file:', TRIM(filename)

    write(unit, *) '1 2'
    write(unit, *) 'star_age model_number'
    write(unit, *) s%star_age, s%model_number
    write(unit, *) ''
    write(unit, *) '1 2 3 4 5 6 7 8'
    write(unit, *) 'k mass gradr grada gradL_composition_term mixing_type x y'

    do k = 1, s%nz
       write(unit, *) k, s%m(k)/msun, s%gradr(k), s%grada_face(k), &
            s%gradL_composition_term(k), s%mlt_mixing_type(k), s%xa(iso_h1,k), s%xa(iso_he4,k)
    end do

    close(unit)

    ! Restore gradL_composition_term

    s%gradL_composition_term(1:s%nz) = gradL_composition_term

    ! Finish

    return

  end subroutine dump_snapshot_

  !****

  subroutine alloc_saved_data_ (s, sd)

    type(star_info), pointer :: s
    type(saved_data)         :: sd

    ! Allocate cell data arrays

    allocate(sd%xa(s%species,s%nz))

    allocate(sd%lnd(s%nz))
    allocate(sd%rho(s%nz))

    allocate(sd%lnPgas(s%nz))
    allocate(sd%Pgas(s%nz))

    allocate(sd%update_mode(s%nz))

    allocate(sd%gradL_composition_term(s%nz))

    ! Finish

  end subroutine alloc_saved_data_

  !****

  subroutine save_model_ (s, update_mode, kc_t, kc_b, sd)

    type(star_info), pointer        :: s
    integer, intent(in)             :: update_mode(:)
    integer, intent(in)             :: kc_t
    integer, intent(in)             :: kc_b
    type(saved_data), intent(inout) :: sd

    logical, parameter :: TRACE_SAVE_MODEL = .FALSE.

    integer :: kf_t
    integer :: kf_b

    ! Save data for cells kc_t:kc_b and interior faces

    if (TRACE_SAVE_MODEL) then
       write(*,*) '  Saving model'
    end if

    ! Save cell data

    if (TRACE_SAVE_MODEL) then
       write(*,*) '    kc_t:', kc_t
       write(*,*) '    kc_b:', kc_b
    end if

    sd%update_mode(kc_t:kc_b) = update_mode(kc_t:kc_b)
    
    sd%xa(:,kc_t:kc_b) = s%xa(:,kc_t:kc_b)

    sd%lnd(kc_t:kc_b) = s%lnd(kc_t:kc_b)
    sd%Rho(kc_t:kc_b) = s%Rho(kc_t:kc_b)

    sd%lnPgas(kc_t:kc_b) = s%lnPgas(kc_t:kc_b)
    sd%Pgas(kc_t:kc_b) = s%Pgas(kc_t:kc_b)

    ! Save face data

    kf_t = kc_t
    kf_b = kc_b + 1

    if (TRACE_SAVE_MODEL) then
       write(*,*) '    kf_t:', kf_t
       write(*,*) '    kf_b:', kf_b
    end if
    
    sd%gradL_composition_term(kf_t+1:kf_b-1) = s%gradL_composition_term(kf_t+1:kf_b-1)

    ! Store indices (used when we call restore_model_)

    sd%kc_t = kc_t
    sd%kc_b = kc_b

    ! Finish

    return

  end subroutine save_model_

  !****

  subroutine restore_model_ (s, update_mode, sd)

    type(star_info), pointer        :: s
    integer, intent(inout)          :: update_mode(:)
    type(saved_data), intent(inout) :: sd

    logical, parameter :: TRACE_RESTORE_MODEL = .FALSE.

    integer :: kc_t
    integer :: kc_b
    integer :: kf_t
    integer :: kf_b

    ! Restore data from previous call to save_model_

    if (TRACE_RESTORE_MODEL) then
       write(*,*) '  Restoring model'
    end if

    ! Restore cell data

    kc_t = sd%kc_t
    kc_b = sd%kc_b

    if (TRACE_RESTORE_MODEL) then
       write(*,*) '    kc_t:', kc_t
       write(*,*) '    kc_b:', kc_b
    end if

    update_mode(kc_t:kc_b) = sd%update_mode(kc_t:kc_b)
    
    s%xa(:,kc_t:kc_b) = sd%xa(:,kc_t:kc_b)
    
    s%lnd(kc_t:kc_b) = sd%lnd(kc_t:kc_b)
    s%rho(kc_t:kc_b) = sd%rho(kc_t:kc_b)
    
    s%lnPgas(kc_t:kc_b) = sd%lnPgas(kc_t:kc_b)
    s%Pgas(kc_t:kc_b) = sd%Pgas(kc_t:kc_b)

    ! Restore face data

    kf_t = kc_t
    kf_b = kc_b + 1

    if (TRACE_RESTORE_MODEL) then
       write(*,*) '    kf_t:', kf_t
       write(*,*) '    kf_b:', kf_b
    end if
    
    s%gradL_composition_term(kf_t+1:kf_b-1) = sd%gradL_composition_term(kf_t+1:kf_b-1)

    ! Update the model set those quantities that are not stored
    ! explicitly in sd

    call update_model_(s, update_mode, kc_t, kc_b)

    ! Finish

    return

  end subroutine restore_model_

end module conv_premix

