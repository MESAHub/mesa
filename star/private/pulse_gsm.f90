! ***********************************************************************
!
!   Copyright (C) 2021  Rich Townsend & The MESA Team
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

module pulse_gsm

  ! Uses

  use star_private_def
  use forum_m, only: hdf5io_t, CREATE_FILE

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: write_gsm_data

contains

  subroutine write_gsm_data (id, filename, global_data, point_data, ierr)

    integer, intent(in)      :: id
    character(*), intent(in) :: filename
    real(dp), intent(in)     :: global_data(:)
    real(dp), intent(in)     :: point_data(:,:)
    integer, intent(out)     :: ierr

    type(star_info), pointer :: s
    type(hdf5io_t)           :: hi

    ! Write GYRE data to a GSM (GYRE stellar model) file

    call get_star_ptr(id, s, ierr)
    if (ierr /= 0) then
       write(*,*) 'bad star id for write_gyre_data'
       return
    end if

    select case(s%gyre_data_schema)
    case(110,120)
    case default
       write(*,*) 'invalid gyre_data_schema'
       ierr = -1
       return
    end select

    ! Open the file

    hi = hdf5io_t(filename, CREATE_FILE)

    ! Write the data

    call hi%write_attr('n', SIZE(point_data, 2))

    call hi%write_attr('M_star', global_data(1))
    call hi%write_attr('R_star', global_data(2))
    call hi%write_attr('L_star', global_data(3))
    
    call hi%write_attr('version', s%gyre_data_schema)

    select case(s%gyre_data_schema)

    case(110)

       call hi%write_dset('r', point_data(1,:))
       call hi%write_dset('M_r', point_data(2,:))
       call hi%write_dset('L_r', point_data(3,:))
       call hi%write_dset('P', point_data(4,:))
       call hi%write_dset('T', point_data(5,:))
       call hi%write_dset('rho', point_data(6,:))
       call hi%write_dset('nabla', point_data(7,:))
       call hi%write_dset('N2', point_data(8,:))
       call hi%write_dset('Gamma_1', point_data(9,:))
       call hi%write_dset('nabla_ad', point_data(10,:))
       call hi%write_dset('delta', point_data(11,:))
       call hi%write_dset('kap', point_data(12,:))
       call hi%write_dset('kap_kap_T', point_data(13,:))
       call hi%write_dset('kap_kap_rho', point_data(14,:))
       call hi%write_dset('eps', point_data(15,:))
       call hi%write_dset('eps_eps_T', point_data(16,:))
       call hi%write_dset('eps_eps_rho', point_data(17,:))
       call hi%write_dset('Omega_rot', point_data(18,:))

    case(120)

       call hi%write_dset('r', point_data(1,:))
       call hi%write_dset('M_r', point_data(2,:))
       call hi%write_dset('L_r', point_data(3,:))
       call hi%write_dset('P', point_data(4,:))
       call hi%write_dset('T', point_data(5,:))
       call hi%write_dset('rho', point_data(6,:))
       call hi%write_dset('nabla', point_data(7,:))
       call hi%write_dset('N2', point_data(8,:))
       call hi%write_dset('Gamma_1', point_data(9,:))
       call hi%write_dset('nabla_ad', point_data(10,:))
       call hi%write_dset('delta', point_data(11,:))
       call hi%write_dset('kap', point_data(12,:))
       call hi%write_dset('kap_kap_T', point_data(13,:))
       call hi%write_dset('kap_kap_rho', point_data(14,:))
       call hi%write_dset('eps', point_data(15,:))
       call hi%write_dset('eps_eps_T', point_data(16,:))
       call hi%write_dset('eps_eps_rho', point_data(17,:))
       call hi%write_dset('eps_grav', point_data(18,:))
       call hi%write_dset('Omega_rot', point_data(19,:))

    end select

    ! Close the file

    call hi%final()
    
    ! Finish

    ierr = 0

    return

 end subroutine write_gsm_data

end module pulse_gsm
