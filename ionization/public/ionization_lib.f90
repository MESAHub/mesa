! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
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

      module ionization_lib
      
      use const_def, only: dp
      
      implicit none


      contains


      subroutine ionization_init( &
            file_prefix, Z1_suffix, ionization_cache_dir, use_cache, ierr)      
         use ionization_def, only: ion_def_init
         use ion_tables_load, only: Init_ion_tables
         use mod_ionization, only: do_init_ionization
         character(len=*), intent(in) :: file_prefix, Z1_suffix
         character(len=*), intent(in) :: ionization_cache_dir ! '' means use default
         logical, intent(in) :: use_cache
         integer, intent(out) :: ierr
         ierr = 0
         call ion_def_init(ionization_cache_dir)
         call Init_ion_tables(file_prefix, Z1_suffix, use_cache, ierr)
         if (ierr /= 0) return
         call do_init_ionization(ionization_cache_dir, use_cache, ierr) 
      end subroutine ionization_init

      ! EXPERIMENTAL
      ! This routine is currently undocumented, not recommended for publishable work.
      ! Element diffusion uses eval_typical_charge, not eval_ionization.
      subroutine eval_ionization(Z, X, Rho, log10Rho, T, log10T, res, ierr)
         use ionization_def, only: num_ion_vals
         use ion_tables_eval, only: Get_ion_Results
         real(dp), intent(in) :: Z ! metals mass fraction
         real(dp), intent(in) :: X ! hydrogen mass fraction
         real(dp), intent(in) :: Rho, log10Rho ! the density
            ! provide both if you have them. 
            ! else pass one and set the other to = arg_not_provided
         real(dp), intent(in) :: T, log10T ! the temperature
            ! provide both if you have them. 
            ! else pass one and set the other to = arg_not_provided              
         real(dp), intent(inout) :: res(num_ion_vals) ! see ionization_def
         integer, intent(out) :: ierr
         call Get_ion_Results(Z, X, Rho, log10Rho, T, log10T, res, ierr)
      end subroutine eval_ionization

      
      ! average ionic charges based on
      ! C. Paquette, C. Pelletier, G. Fontaine, G. Michaud,
      ! "Diffusion in White Dwarfs: New Results and Comparative Study",
      ! ApJ Supp. Series, 61:197-217, 1986.

      ! Originally used eqn 21 of above paper for depression of the continuum,
      ! but that expression was missing a rho^1/3, so now corrected to match
      ! eqn 3 of
      ! J. Dupuis, G. Fontaine, C. Pelletier, F. Wesemael,
      ! "A Study of Metal Abundance Patterns in Cool White Dwarfs I.
      ! Time-dependent Calculations of Gravitational Settling",
      ! Apj Supp. Series, 82:505-521, 1992

      ! uses ionization potentials from
      ! Allen, C.W., 1973, "Astrophysical Quantities", 3rd edition, pg 37-38.
      
      ! tables go from H to Zn.
      
      ! the routine is written so that it doesn't ever return 0 as a typical charge.
      ! these are "typical" charges rather than average.  the values are whole numbers.

      real(dp) function eval_typical_charge( &
            cid, abar, free_e, T, log10_T, rho, log10_rho)
         use mod_ionization, only: get_typical_charge
         integer, intent(in) :: cid ! chem id such as ihe4.  defined in chem_def.
         real(dp), intent(in) :: abar ! average atomic weight (from chem_lib)
         real(dp), intent(in) :: free_e 
            ! mean number of free electrons per nucleon (from eos_lib)
            ! abar*free_e = (nucleons/particle)*(charge/nucleon) = charge/particle = z1
         real(dp), intent(in) :: T, log10_T, rho, log10_rho
         eval_typical_charge = get_typical_charge( &
            cid, abar, abar*free_e, T, log10_T, rho, log10_rho)    
      end function eval_typical_charge

      ! EXPERIMENTAL
      real(dp) function eval_charge_of_Fe56_in_He4(log10_ne, log10_T, ierr)
         use mod_ionization, only: charge_of_Fe56_in_He4
         real(dp), intent(in) :: log10_ne ! ne=avo*rho*free_e
         real(dp), intent(in) :: log10_T
         integer, intent(out) :: ierr
         eval_charge_of_Fe56_in_He4 = charge_of_Fe56_in_He4(log10_ne, log10_T, ierr)
      end function eval_charge_of_Fe56_in_He4


      
      subroutine create_ion_table_files( &
            in_dir, out_dir_ion, out_dir_eosDT, out_dir_eosPT)
         !use mod_ion_create_tables, only: do_create_ion_table_files
         character (len=*), intent(in) :: &
            in_dir, out_dir_ion, out_dir_eosDT, out_dir_eosPT
         stop 'create_ion_table_files not currently supported -- need to fix interpolation calls'
         !call do_create_ion_table_files( &
         !   in_dir, out_dir_ion, out_dir_eosDT, out_dir_eosPT)
      end subroutine create_ion_table_files
      
      
      subroutine create_table_plot_files
         use ion_table_plot, only: do_create_table_plot_files
         call do_create_table_plot_files
      end subroutine create_table_plot_files
      
      end module ionization_lib

