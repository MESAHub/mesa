module ion_offset
      use const_def
      use math_lib

      implicit none

      logical, parameter :: dbg = .false.
      !logical, parameter :: dbg = .true.
      
      
      private
      public :: compute_ion_offset

      contains

      !> Computes the energy required to fully ionize
      !! the specified composition. This is an offset
      !! between Skye/HELM and FreeEOS.
      !!
      !! @param species Number of species
      !! @param xa Mass fractions (vector)
      !! @param offset Energy in erg/g
      real(dp) function compute_ion_offset(species, xa, chem_id) result(offset)
         use chem_def, only: chem_isos
         integer, intent(in) :: species
         real(dp), intent(in) :: xa(species)
         integer, pointer :: chem_id(:)

         ! First 28 ionization energies, same species as used in FreeEOS, obtained from
         ! NIST (https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html)
         ! using the query 'H-Ni' asking for 'Total binding energy'.
         real(dp), parameter :: ionization_table(28) = [13.598434599702,79.0051545,203.4861711,399.14864,670.9809,&
                                                   1030.1085,1486.058,2043.8428,2715.89,3511.696,4420.0,5451.06,6604.95,&
                                                   7888.53,9305.8,10859.7,12556.4,14400.8,16382.0,18510.0,20788.0,&
                                                   23221.0,25820.0,28582.0,31514.0,34619.0,37899.0,41356.0]

         integer :: k, Z(species)
         real(dp) :: A(species), W(species), ya(species), norm

         ! Get basic species info
         norm = 0d0
         do k=1,species
            A(k) = chem_isos% Z_plus_N(chem_id(k)) ! baryon number
            Z(k) = chem_isos% Z(chem_id(k)) ! charge
            ya(k) = xa(k) / A(k)            ! number fraction (not normalized)
            norm = norm + ya(k)             ! accumulate number fraction normalization
         end do

         ! Normalize number fraction
         do k=1,species
            ya(k) = ya(k) / norm
         end do

         ! Compute offset in eV/baryon
         offset = 0d0
         do k=1,species
            if (Z(k) <= 28 .and. Z(k) >= 1) then
               offset = offset + ionization_table(Z(k)) * ya(k)
            end if
         end do

         ! Translate that into erg/g
         offset = offset * ev2erg / amu

      end function compute_ion_offset

end module ion_offset
