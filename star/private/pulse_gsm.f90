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
  use hdf5io_lib

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

    type(hdf5io_t) :: hi

    ! Write GYRE data to a GSM (GYRE stellar model) file

    ! Open the file

    hi = hdf5io_t(filename, CREATE_FILE)

    ! Write the data

    call hi%write_attr('n', SIZE(point_data, 2))

    call hi%write_attr('M_star', global_data(1))
    call hi%write_attr('R_star', global_data(2))
    call hi%write_attr('L_star', global_data(3))
    
    call hi%write_attr('version', GSM_MODEL_VERSION)

    associate ( &
         r => point_data(1,:), &
         m => point_data(2,:), &
         L => point_data(3,:), &
         P => point_data(4,:), &
         T => point_data(5,:), &
         rho => point_data(6,:), &
         nabla => point_data(7,:), &
         N2 => point_data(8,:), &
         Gamma_1 => point_data(9,:), &
         nabla_ad => point_data(10,:), &
         delta => point_data(11,:), &
         kap => point_data(12,:), &
         kap_T => point_data(13,:), &
         kap_rho => point_data(14,:), &
         eps => point_data(15,:), &
         eps_T => point_data(16,:), &
         eps_rho => point_data(17,:), &
         omega => point_data(18,:))

      call hi%write_dset('r', r)
      call hi%write_dset('M_r', m)
      call hi%write_dset('L_r', L)

      call hi%write_dset('P', P)
      call hi%write_dset('rho', rho)
      call hi%write_dset('T', T)
      
      call hi%write_dset('N2', N2)
      call hi%write_dset('Gamma_1', Gamma_1)
      call hi%write_dset('nabla_ad', nabla_ad)
      call hi%write_dset('delta', delta)
      call hi%write_dset('nabla', nabla)

      call hi%write_dset('kap', kap)
      call hi%write_dset('kap_kap_T', kap_T)
      call hi%write_dset('kap_kap_rho', kap_rho)

      call hi%write_dset('eps', eps)
      call hi%write_dset('eps_eps_T', eps_T)
      call hi%write_dset('eps_eps_rho', eps_rho)

      call hi%write_dset('Omega_rot', omega)

    end associate

    ! Close the file

    call hi%final()
    
    ! Finish

    ierr = 0

    return

 end subroutine write_gsm_data

end module pulse_gsm
