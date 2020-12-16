! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton, Ed Brown & The MESA Team
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
!
! ***********************************************************************

      module chem_lib
      use chem_def, only: chem_has_been_initialized
      use const_def, only: dp
      use math_lib
      
      implicit none


      contains
      
      
      subroutine chem_init(isotopes_filename, ierr) 
         ! uses mesa_data_dir from const_def
         use chem_def
         use chem_isos_io, only: do_read_chem_isos
         use lodders_mod, only : read_lodders03_data
         character (len=*), intent(in) :: isotopes_filename
         integer, intent(out) :: ierr

         ierr = 0
         if (chem_has_been_initialized) return
         call do_read_chem_isos(isotopes_filename, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in do_read_chem_isos'
            return
         end if
         call init_chem_tables
         call read_lodders03_data('lodders03.data',ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in read_lodders03_data'
            return
         end if
         chem_has_been_initialized = .true.
         call set_some_isos

      end subroutine chem_init


      subroutine chem_shutdown ()

        use chem_def
        use utils_lib, only: integer_dict_free

        if (.NOT. chem_has_been_initialized) return

        call free_lodders03_table()
        call free_nuclide_data(chem_isos)

        ! Free dictionaries

        if (ASSOCIATED(Xsol_names_dict)) call integer_dict_free(Xsol_names_dict)
        if (ASSOCIATED(chem_element_names_dict)) call integer_dict_free(chem_element_names_dict)
        if (ASSOCIATED(chem_isos_dict)) call integer_dict_free(chem_isos_dict)
        if (ASSOCIATED(category_names_dict)) call integer_dict_free(category_names_dict)

        chem_has_been_initialized = .FALSE.

      end subroutine chem_shutdown


      function lodders03_element_atom_percent(nuclei) result(percent)
      ! Katharina Lodders,
      ! "Solar System Abundances and Condensation Temperatures of the Elements",
      ! ApJ, 591, 1220-1247 (2003).
      ! Table 6: Abundances of the Isotopes in the Solar System
      
      ! These are element atom percentages (i.e., by number, not by mass)
      
      ! NOTE: The data here stops at ge -- the table in the paper goes to 92U
      
      ! TO DO: add the rest of the info from the paper
            use chem_def
            use lodders_mod
            character(len=*), intent(in) :: nuclei
            real(dp) :: percent
            integer :: ierr
            
            if (.not. chem_has_been_initialized) then
               write(*,*) 'must call chem_init before calling any other routine in chem_lib'
              percent = -1.0d0
               return
            end if
            percent = get_lodders03_isotopic_abundance(nuclei, ierr)
            if (ierr /= 0) percent = 0.0d0
      end function lodders03_element_atom_percent
      

      ! returns the index of a particular nuclide in a particular set
      ! returns nuclide_not_found if name not found
      function get_nuclide_index_in_set(nuclei, set) result(indx)
         use chem_def
         use nuclide_set_mod
         character(len=*), intent(in) :: nuclei
         type(nuclide_set), dimension(:), intent(in) :: set
         integer :: indx
         character(len=iso_name_length) :: name
         name = nuclei
         indx = rank_in_set(name, set)
      end function get_nuclide_index_in_set


      subroutine generate_nuclide_set(names, set)
         use chem_def
         use nuclide_set_mod, only: sort_nuclide_set
         character(len=iso_name_length), dimension(:), intent(in) :: names
         type(nuclide_set), dimension(size(names)), intent(out) :: set
         integer :: i
         set = [(nuclide_set(names(i), i), i=1, size(names))]
         call sort_nuclide_set(set)
      end subroutine generate_nuclide_set
      
      
      subroutine basic_composition_info( &
            num_isos, chem_id, x, xh, xhe, z, &
            abar, zbar, z2bar, z53bar, ye, mass_correction, sumx)
         integer, intent(in) :: num_isos
         integer, intent(in) :: chem_id(:) ! (num_isos) ! the nuclide indices for the entries in x
         real(dp), intent(in)  :: x(:) ! (num_isos) ! baryon fractions.  should sum to 1.0
         real(dp), intent(out) ::  &
               xh, xhe, z, abar, zbar, z2bar, z53bar, ye, mass_correction, sumx
         real(dp), dimension(0) :: dabar_dx, dzbar_dx, dmc_dx
         call get_composition_info( &
            num_isos, chem_id, x, xh, xhe, z, &
            abar, zbar, z2bar, z53bar, ye, mass_correction, &
            sumx, .true., dabar_dx, dzbar_dx, dmc_dx)
      end subroutine basic_composition_info
      
         
      subroutine composition_info( &
            num_isos, chem_id, x, xh, xhe, z, &
            abar, zbar, z2bar, z53bar, ye, mass_correction, &
            sumx, dabar_dx, dzbar_dx, dmc_dx)
         integer, intent(in) :: num_isos
         integer, intent(in) :: chem_id(:) ! (num_isos) ! the nuclide indices for the entries in x
         real(dp), intent(in)  :: x(:) ! (num_isos) ! baryon fractions.  should sum to 1.0
         real(dp), intent(out) ::  &
               xh, xhe, z, abar, zbar, z2bar, z53bar, ye, mass_correction, sumx
         real(dp), dimension(:) :: dabar_dx, dzbar_dx, dmc_dx
         call get_composition_info( &
            num_isos, chem_id, x, xh, xhe, z, &
            abar, zbar, z2bar, z53bar, ye, mass_correction, &
            sumx, .false., dabar_dx, dzbar_dx, dmc_dx)
      end subroutine composition_info
      
         
      subroutine get_composition_info( &
            num_isos, chem_id, x, xh, xhe, xz, &
            abar, zbar, z2bar, z53bar, ye, mass_correction, &
            sumx, skip_partials, dabar_dx, dzbar_dx, dmc_dx)

         ! here's a reminder of definitions:
         ! X(i) ion baryon fraction (g/g)
             ! A(i) ion atomic mass number
         ! W(i) ion atomic weight (g/mole)
         ! Z(i) ion charge (number of protons)
         ! Y(i) = X(i)/A(i), ion abundance 
         ! n(i) = rho*avo*Y(i), ion number density (g/cm^3)*(#/mole)*(mole/g) -> (#/cm^3)
         
         ! abar = sum(n(i)*A(i))/sum(n(i)) -- average atomic mass number
         ! zbar = sum(n(i)*Z(i))/sum(n(i)) -- average charge number
         ! z2bar = sum(n(i)*Z(i)^2)/sum(n(i)) -- average charge^2
         ! n = rho*avo/abar -- total number density (#/cm^3)
         ! ye = zbar/abar -- average charge per baryon = proton fraction
         ! xh = mass fraction hydrogen
         ! xhe = mass fraction helium
         ! xz = mass fraction metals
         ! mass_correction = sum(n(i)*W(i))/sum(n(i)*A(i)) --
         ! (mass density) = (baryon density) * m_u * mass_correction

         use chem_def
         integer, intent(in) :: num_isos
         integer, intent(in) :: chem_id(:) ! (num_isos) ! the nuclide indices for the entries in x
         real(dp), intent(in)  :: x(:) ! (num_isos) ! baryon fractions.  should sum to 1.0
         real(dp), intent(out) ::  &
               xh, xhe, xz, abar, zbar, z2bar, z53bar, ye, mass_correction, sumx
         logical, intent(in) :: skip_partials
         real(dp), dimension(:) :: dabar_dx, dzbar_dx, dmc_dx
     
         real(dp), dimension(num_isos) :: y, z, w, a
         integer :: i, cid, iz
         if (.not. chem_has_been_initialized) then
            write(*,*) 'must call chem_init before calling any other routine in chem_lib'
            xh=0; xhe=0; abar=0; zbar=0; z2bar=0; z53bar=0; ye=0; mass_correction = 0; sumx=0
            return
         end if
         xh = 0d0; xhe = 0d0
         do i=1,num_isos
            cid = chem_id(i)
            iz = chem_isos% Z(cid)
            z(i) = iz
            y(i) = x(i)/dble(chem_isos% Z_plus_N(cid))
            w(i) = chem_isos% W(cid)
            a(i) = chem_isos% Z_plus_N(cid)
            select case(iz)
               case (1)
                  xh = xh + x(i)
               case (2)
                  xhe = xhe + x(i)
            end select
         end do
         
         xz = max(0d0, min(1d0, 1d0 - (xh + xhe)))
         sumx = sum(x(1:num_isos))     ! this should be one, always, since we define x as a baryon fraction
         abar = sumx / sum(y(1:num_isos))
         zbar = sum(y(1:num_isos)*z(1:num_isos)) * abar
         z2bar = sum(y(1:num_isos)*z(1:num_isos)*z(1:num_isos)) * abar
         z53bar = sum(y(1:num_isos)*z(1:num_isos)*z(1:num_isos)) * abar
         ye = zbar/abar ! assume complete ionization
         mass_correction = sum(y(1:num_isos)*w(1:num_isos))

         z53bar = 0d0
         do i=1,num_isos
            z53bar = z53bar + y(i) * chem_isos% Z53(chem_id(i))
         end do
         z53bar = z53bar * abar

         if (skip_partials) return
         
         do i=1,num_isos
            dabar_dx(i) = abar*(a(i)-abar)/a(i)/sumx
            dzbar_dx(i) = abar*(z(i)-zbar)/a(i)/sumx
            dmc_dx(i) = w(i)/a(i) - mass_correction
         end do
         
      end subroutine get_composition_info
      
      
      ! Q:  Is positron annihilation actually included in Qtotal?
      ! A: (from Frank Timmes)
      !
      !   the formula used in the code is the atomic mass excess - not the nuclear mass excess - 
      !   in terms of the binding energy and Deltas. as bill notes, this formulation
      !   implicitly takes care of the electron masses. i can see why some confusion may 
      !   arise at first blush because of the notation used in the code - its says “Mp” but 
      !   really means “Mh” where Mh = m_p + m_e.
      !
      !   p + p -> d + e^+ + nu
      !
      !   let’s do it by hand and then by code, neglecting binding energy terms of 
      !   the atom (13.6 eV and such) as they are a million times smaller than the nuclear terms.
      !
      !   by hand:
      !
      !   mass excess left-hand-side = twice hydrogen mass excess = 
      !   2 * (m_h - 1 amu) = 2 (m_p + m_e - 1 amu)
      !
      !   mass excess right-hand-side = 
      !   m_h + m_n - B(d) - 2 amu = m_p + m_e + m_n - B(d) - 2 amu
      !
      !   Q = left - right = m_p + m_e - m_n + B(d) 
      !
      !   which may be written (adding and subtracting m_p)
      !
      !   Q = 2 m_p + m_e - (m_p + m_n - B(d))
      !    = (twice nuclear mass of proton) + (electron mass) - (nuclear mass of deuterium)
      !
      !
      !   by code (with interpretation):
      !
      !   first loop
      !   Q = -del_Mp = -(m_h - 1 amu) = -(m_p + m_e - 1 amu)
      !   Qtotal = del_Mp = m_p + m_e - 1 amu
      !
      !   second loop
      !   Q = -del_Mp = -(m_h - 1 amu) = -(m_p + m_e - 1 amu)
      !   Qtotal = 2 del_Mp = 2 (m_p + m_e - 1 amu)     ! note this is the same as the left-hand side term above
      !
      !
      !   third loop
      !   Q = B(d) - del_Mp - del_Mn 
      !    = B(d) - (m_h - 1 amu) - (m_n - 1 amu) 
      !    = B(d) - (m_p + m_e - 1 amu) - (m_n - 1 amu)
      !    = B(d) - (mp + m_e + m_n - 2 amu)             ! note this is the negative of the right-hand-side term above
      !
      !   Qtotal = 2 del_Mp + [B(d) - del_Mp - del_Mn]
      !         = m_p + m_e - m_n + B(d)                 ! same as above
      !
      !
      !   i think the confusion lay between the code nomenclature
      !   and nuclear vs atomic mass excesses. 
      
      real(dp) function reaction_Qtotal(num_in,num_out,reactants,nuclides)
         use chem_def
         integer, intent(in) :: num_in,num_out,reactants(:)
         type(nuclide_data), intent(in) :: nuclides
         integer :: j, l
         real(dp) :: Q
         reaction_Qtotal = 0d0
         do j=1,num_in+num_out
            l = reactants(j)
            Q = get_Q(nuclides,l)
            if (j <= num_in) then
               reaction_Qtotal = reaction_Qtotal - Q
            else
               reaction_Qtotal = reaction_Qtotal + Q
            end if
  
         end do
      end function reaction_Qtotal
            
      
      integer function chem_get_element_id(cname) 
      ! NOTE: this is for elements like 'h', not for isotopes like 'h1'
      ! use chem_get_iso_id for looking up isotope names
         use chem_def
         use utils_lib
         character (len=*), intent(in)  :: cname 
         ! name of the element (e.g. 'h', 'he', 'ne')
         ! same names as in chem_element_Name
         ! returns id for the element if there is a matching name
         ! returns 0 otherwise.
         integer :: ierr, value
         if (.not. chem_has_been_initialized) then
            write(*,*) 'must call chem_init before calling any other routine in chem_lib'
            chem_get_element_id = -1
            return
         end if
         call integer_dict_lookup(chem_element_names_dict, cname, value, ierr)
         if (ierr /= 0) value = -1
         chem_get_element_id = value
      end function chem_get_element_id
      
      
      real(dp) function chem_Xsol(nam)
         character (len=*), intent(in)  :: nam 
            ! name of the isotope (e.g. 'h1', 'he4', 'ne20')
         real(dp) :: z, a, xelem
         integer :: ierr
         if (.not. chem_has_been_initialized) then
            write(*,*) 'must call chem_init before calling any other routine in chem_lib'
            chem_Xsol = -1
            return
         end if
         call chem_get_solar(nam, z, a, xelem, ierr)
         if (ierr /= 0) then
            chem_Xsol = 0d0
         else
            chem_Xsol = xelem
         end if
      end function chem_Xsol
      

      subroutine chem_get_solar(nam, z, a, xelem, ierr)
         use chem_def
         use utils_lib
         ! returns data from Anders and Grevesse, 1989
         character (len=*), intent(in)  :: nam 
            ! name of the isotope (e.g. 'h1', 'he4', 'ne20')
            ! note that these names match those in the nuclear net library iso_Names array
            ! but some net isotopes are not here (ex/ be7, n13, o14, o15, f17, f18, ... fe52, ni56 )
            ! and some isotopes are here that are not in the nets (ex/ hf176)
         real(dp), intent(out) :: z ! charge
         real(dp), intent(out) :: a ! number of nucleons (protons and neutrons)
         real(dp), intent(out) :: xelem ! elemental mass fraction associated with this isotope
         integer, intent(out) :: ierr ! == 0 means AOK; == -1 means did not find the requested name
         integer :: i
         ierr = 0
         if (.not. chem_has_been_initialized) then
            write(*,*) 'must call chem_init before calling any other routine in chem_lib'
            ierr = -1; z = 0; a = 0; xelem = 0; return
         end if
         call integer_dict_lookup(Xsol_names_dict, nam, i, ierr)
         if (ierr /= 0) then
            z = 0; a = 0; xelem = 0; return
         end if
         z = dble(izsol(i))
         a = dble(iasol(i))
         xelem = solx(i)
      end subroutine chem_get_solar

      
      
      ! given an array of Z, A, returns an array of names in chem_isos format
      subroutine generate_nuclide_names(Z, A, names)
         use chem_def
         integer, dimension(:), intent(in) :: Z, A
         character(len=iso_name_length), dimension(size(Z)), intent(out) :: names
         integer :: i, ierr, count_isomer
         character(len=80) :: message
         logical :: use_al26_isomers

         count_isomer = 0
         use_al26_isomers = .false.
         do i = 1, size(Z)
           if (A(i) == 26 .and. Z(i) == 13) count_isomer = count_isomer + 1
         enddo
         if (count_isomer > 1) use_al26_isomers = .true.

         count_isomer = 1
         do i = 1, size(Z)
            select case(Z(i))
            case (0)
               names(i) = el_name(Z(i))
            case (1:max_el_z)
               write(names(i), '(a, i0)') trim(el_name(Z(i))), A(i)
               if (A(i) == 26 .and. Z(i) == 13 .and. use_al26_isomers) then
                  count_isomer = count_isomer + 1
                  !if (count_isomer > 1) names(i) = al_isomers(count_isomer)
                  names(i) = al_isomers(count_isomer)
               end if
            case default
               write(*, '(a, i0, a, i0)') 'warning: ', Z(i), ' greater than Zmax = ', max_el_z
               names(i) = '*****'
               ierr = nuclide_not_found
            end select
            names(i) = adjustr(names(i))
         end do
      end subroutine generate_nuclide_names
      

      subroutine generate_long_nuclide_names(Z, A, long_names)
         use chem_def
         integer, dimension(:), intent(in) :: Z, A
         character(len=long_name_length), dimension(size(Z)), intent(out) :: long_names
         integer :: i, ierr, count_isomer
         character(len=80) :: message
         logical :: use_al26_isomers

         count_isomer = 0
         use_al26_isomers = .false.
         do i = 1, size(Z)
           if (A(i) == 26 .and. Z(i) == 13) count_isomer = count_isomer + 1
         enddo
         if (count_isomer > 1) use_al26_isomers = .true.

         count_isomer = 1
         do i = 1, size(Z)
            select case(Z(i))
            case (0) ! neutrons are special?
               long_names(i) = el_long_name(Z(i))
            case (1:max_el_z)
               write(long_names(i), '(a, "-", i0)') trim(el_long_name(Z(i))), A(i)
               if (A(i) == 26 .and. Z(i) == 13 .and. use_al26_isomers) then
                  count_isomer = count_isomer + 1
                  !if (count_isomer > 1) long_names(i) = long_al_isomers(count_isomer)
                                    long_names(i) = long_al_isomers(count_isomer)
               end if
            case default
               write(*, '(a, i0, a, i0)') 'warning: ', Z(i), ' greater than Zmax = ', max_el_z
               long_names(i) = '********'
               ierr = nuclide_not_found
            end select
         end do
      end subroutine generate_long_nuclide_names
      
      
      ! nuclide information comes from the chem_isos tables
      ! the storage container for the data is called 'chem_isos'
      ! it has name, A, Z, N, spin, and B for each nuclide
      ! use the function chem_get_iso_id to find the index given the name.      
      integer function chem_get_iso_id(cname)
         use chem_def, only: get_nuclide_index
         character(len=*), intent(in) :: cname
         chem_get_iso_id = get_nuclide_index(cname)
      end function chem_get_iso_id  
            
      
      integer function lookup_ZN(Z,N)
         integer, intent(in) :: Z, N
         lookup_ZN = lookup_ZN_isomeric_state(Z,N,0)
      end function lookup_ZN  
            
      
      integer function lookup_ZN_isomeric_state(Z,N,isomeric_state)
         use chem_def, only: chem_isos, num_chem_isos
         integer, intent(in) :: Z, N, isomeric_state
         integer :: cid, i       
         iso_loop: do cid = 1, num_chem_isos
            if (chem_isos% Z(cid) == Z .and. chem_isos% N(cid) == N) then
               if (chem_isos% isomeric_state(cid) == isomeric_state) then
                  lookup_ZN_isomeric_state = cid
                  return
               end if
               do i = cid+1, num_chem_isos
                  if (chem_isos% Z(i) /= Z .or. chem_isos% N(i) /= N) exit iso_loop
                  if (chem_isos% isomeric_state(i) == isomeric_state) then
                     lookup_ZN_isomeric_state = i
                     return
                  end if
               end do
            end if
         end do iso_loop
         lookup_ZN_isomeric_state = 0 ! indicating failure
      end function lookup_ZN_isomeric_state  
            
      
      integer function rates_category_id(cname)
         use chem_def, only: category_names_dict
         use utils_lib
         character (len=*), intent(in)  :: cname 
         ! returns id for the category if there is a matching name
         ! returns 0 otherwise.
         integer :: ierr, value
         call integer_dict_lookup(category_names_dict, cname, value, ierr)
         if (ierr /= 0) value = 0
         rates_category_id = value
      end function rates_category_id
      

      function binding_energy(nuclides, Y) result (B)
         use chem_def
         type(nuclide_data), intent(in) :: nuclides
         real(dp), dimension(size(nuclides% binding_energy)), intent(in) :: Y
         real(dp) :: B
         B = dot_product(nuclides% binding_energy, Y)
      end function binding_energy
      
      ! returns mass excess in MeV
      function get_mass_excess(nuclides,chem_id) result (mass_excess)
         use chem_def
         type(nuclide_data), intent(in) :: nuclides
         integer, intent(in) :: chem_id
         real(dp) :: mass_excess
         logical :: use_nuclides_mass_excess=.false.
         
         ! These should be identical but can have slight ~ulp difference
         ! due to floating point maths
         if(use_nuclides_mass_excess)then
            mass_excess = nuclides% mass_excess(chem_id)
         else
            mass_excess = nuclides% Z(chem_id)*del_Mp + nuclides% N(chem_id)*del_Mn -&
                        nuclides% binding_energy(chem_id)
         end if   
               
      end function
      
      function get_Q(nuclides,chem_id) result (q)
         use chem_def
         type(nuclide_data), intent(in) :: nuclides
         integer, intent(in) :: chem_id
         real(dp) :: q
         
         !Minus the mass excess
         q=-get_mass_excess(nuclides,chem_id)
               
      end function

      ! returns the indx corresponding to Tpart just less than T9
      ! T9 is the temperature in units of GK
      ! returns a value of 0 or npart if value is less than the minimum or maximum of the partition function
      ! temperature array Tpart
      function get_partition_fcn_indx(T9) result(indx)
         use chem_def, only: Tpart, npart
         real(dp), intent(in) :: T9
         integer :: indx
         integer, parameter :: max_iterations = 8
         integer :: low, high, mid, i

         low = 1
         high = npart

         if (T9 < Tpart(low)) then
            indx = low-1
            return
         end if
         if (T9 > Tpart(high)) then
            indx = high + 1
         end if

         do i = 1, max_iterations
            if (high-low <= 1) then
               indx = low
               return
            end if
            mid = (high+low)/2
            if (T9 < Tpart(mid)) then
               high = mid
            else
               low = mid
            end if
         end do
         ! should never get here
         indx = low-1
         
      end function get_partition_fcn_indx

      ! Given a the chem_id's and abundances for a set of isotopes
      ! return the abundances of the stable isotopes where the unstable ones
      ! have decayed to the stable versions.
      ! Code from Frank Timmes "decay.zip" 
      subroutine get_stable_mass_frac(chem_id,num_species,abun_in,abun_out) 
         use chem_def
         integer,intent(in),dimension(:) :: chem_id
         integer,intent(in) :: num_species
         real(dp),dimension(:),intent(in) :: abun_in
         real(dp),dimension(:),intent(out) :: abun_out
         integer :: i,j,a,z,n

         abun_out(1:solsiz)=0.d0

         outer: do i=1,num_species
            Z=chem_isos%Z(chem_id(i))
            a=z+chem_isos%n(chem_id(i))
            inner: do j=1,solsiz
               if (a .ne. iasol(j)) then
                  cycle inner
               end if
               if ((jcode(j) .eq. 0) .or. &
                   (z.ge.izsol(j) .and. jcode(j).eq.1) .or. &
                   (z.le.izsol(j) .and. jcode(j).eq.2) .or. &
                   (z.eq.izsol(j) .and. jcode(j).eq.3) ) then
                     abun_out(j) = abun_out(j) + abun_in(i)
               end if
            end do inner
         end do outer

         !Normalise results
         abun_out(1:solsiz)=abun_out(1:solsiz)/sum(abun_out(1:solsiz))
         
      end subroutine get_stable_mass_frac


      real(dp) function chem_M_div_h(x,z,zfrac_choice) ! Returns [M/H] 
         use chem_def
         use utils_lib, only: mesa_error
         real(dp), intent(in) :: x ! Hydrogen fraction
         real(dp), intent(in) :: z ! metal fraction
         integer, intent(in) :: zfrac_choice ! See chem_def, *_zfracs options
      
         real(dp) :: zsolar,ysolar
         
         zsolar = 0d0
         ysolar = 0d0
         select case(zfrac_choice)
            case(AG89_zfracs)
               zsolar = AG89_zsol
               ysolar = AG89_ysol
            case(GN93_zfracs)
               zsolar = GN93_zsol
               ysolar = GN93_ysol
            case(GS98_zfracs)
               zsolar = GS98_zsol
               ysolar = GS98_ysol
            case(L03_zfracs)
               zsolar = L03_zsol
               ysolar = L03_ysol
            case(AGS05_zfracs)
               zsolar = AGS05_zsol
               ysolar = AGS05_ysol
            case(AGSS09_zfracs)
               zsolar = AGSS09_zsol
               ysolar = AGSS09_ysol
            case(L09_zfracs)
               zsolar = L09_zsol
               ysolar = L09_ysol
            case(A09_Prz_zfracs)
               zsolar = A09_Prz_zsol
               ysolar = A09_Prz_ysol 
            case default
               call mesa_error(__FILE__,__LINE__,"Bad zfrac_choice")
         end select
         
         chem_M_div_h = log10(z/x)-log10(zsolar/(1.d0-zsolar-ysolar))
      
      
      end function chem_M_div_h


      end module chem_lib

