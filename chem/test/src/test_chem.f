      module mod_test_chem
      use chem_lib
      use chem_def
      use const_def, only: dp
      use utils_lib, only: mesa_error
      
      implicit none
      

      contains
      
      
      subroutine do_test
         use const_lib, only: const_init
         use math_lib, only: math_init
         integer :: ierr, i
         character (len=32) :: my_mesa_dir 
         my_mesa_dir = '../..'         
         call const_init(my_mesa_dir,ierr)     
         if (ierr /= 0) then
            write(*,*) 'const_init failed'
            call mesa_error(__FILE__,__LINE__)
         end if        
         call math_init()
         call chem_init('isotopes.data', ierr)
         if (ierr /= 0) then
            write(*,*) 'FATAL ERROR: failed in chem_init'
            call mesa_error(__FILE__,__LINE__)
         end if
         call do_test_lodders
         call do_tests
         do i=1,50
            call do_tests
            write(*,*) 'done i', i
         end do
      end subroutine do_test
      
      
      subroutine do_tests
      
         write(*,*) 'call do_test_chem'
         flush(6)
         call do_test_chem
         
         write(*,*) 'call do_test_composition_info'
         flush(6)
         call do_test_composition_info

         write(*,*) 'call do_test_Qtotal'
         flush(6)
         call do_test_Qtotal
         
         write(*,*) 'done do_tests'
         flush(6)
         
      end subroutine do_tests


      subroutine do_test_chem
         real(dp) :: c, n, o, cno, w
         integer :: ic12, in14, io16
         include 'formats'
         
         ic12 = get_nuclide_index('c12')
         in14 = get_nuclide_index('n14')
         io16 = get_nuclide_index('o16')

         write(*,*)
         write(*,*)
         write(*,1) 'chem_W(io16)', chem_isos% W(io16)
         flush(6)
         write(*,1) 'chem_Z(io16)', dble(chem_isos% Z(io16))
         flush(6)
         write(*,1) 'chem_binding_energy(io16)', chem_isos% binding_energy(io16)
         flush(6)
         write(*,*)
         write(*,*) 'chem_name(io16) ', chem_isos% name(io16)
         flush(6)
         write(*,*) 'get_nuclide_index("o16") == io16', get_nuclide_index("o16") == io16
         flush(6)
         write(*,*)
         write(*,*) 'chem_element_Name(e_he) ', chem_element_Name(e_he)
         flush(6)
         write(*,*) 'chem_get_element_id("he") == e_he', chem_get_element_id("he") == e_he
         flush(6)
         write(*,*)
         write(*,*)
         write(*,1) 'Anders & Grevesse 1989 zsol', zsol
         write(*,1) 'Anders & Grevesse 1989 yesol', yesol
         write(*,*)
         write(*,*)
         write(*,*) 'cno fraction by mass of Z'
         write(*,*)
         flush(6)
         
         c = chem_Xsol('c12')
         n = chem_Xsol('n14')
         o = chem_Xsol('o16')
         !write(*,1) 'c', c
         !write(*,1) 'n', n
         !write(*,1) 'o', o
         cno = c + n + o
         write(*,1) 'Anders & Grevesse 1989', cno/zsol
         write(*,*)         
         cno = GN93_element_zfrac(e_c) + GN93_element_zfrac(e_n) + GN93_element_zfrac(e_o)
         write(*,1) 'Grevesse and Noels 1993', cno
         write(*,*)         
         cno = GS98_element_zfrac(e_c) + GS98_element_zfrac(e_n) + GS98_element_zfrac(e_o)
         write(*,1) 'Grevesse and Sauval 1998', cno
         write(*,*)         
         cno = L03_element_zfrac(e_c) + L03_element_zfrac(e_n) + L03_element_zfrac(e_o)
         write(*,1) 'Lodders 2003', cno
         write(*,*)         
         cno = AGS05_element_zfrac(e_c) + AGS05_element_zfrac(e_n) + AGS05_element_zfrac(e_o)
         write(*,1) 'Asplund, Grevesse and Sauval 2005', cno
         write(*,*)
         return
         
         write(*,1) 'Grevesse and Sauval 1998 C', GS98_element_zfrac(e_c)
         write(*,*)         
         write(*,1) 'Grevesse and Sauval 1998 N', GS98_element_zfrac(e_n)
         write(*,*)         
         write(*,1) 'Grevesse and Sauval 1998 O', GS98_element_zfrac(e_o)
         write(*,*)         
         write(*,1) 'Grevesse and Sauval 1998 Ne', GS98_element_zfrac(e_ne)
         write(*,*)         
         write(*,1) 'Grevesse and Sauval 1998 Mg', GS98_element_zfrac(e_mg)
         write(*,*)         
         write(*,1) 'Grevesse and Sauval 1998 Si', GS98_element_zfrac(e_si)
         write(*,*)         
         write(*,1) 'Grevesse and Sauval 1998 Fe', GS98_element_zfrac(e_fe)
         write(*,*)  
         return
         
                
         stop
         
      end subroutine do_test_chem 
      
      
         subroutine do_test_composition_info
            integer, parameter :: num_species = 2
            integer, dimension(num_species) :: ids
            real(dp), dimension(num_species) :: X
            real(dp) :: xh, xhe, xz, abar, zbar, z2bar, z53bar, ye, mass_correction, sumx
            real(dp) :: dabar_dx(num_species), dzbar_dx(num_species), dmc_dx(num_species)
            character (len=*), parameter :: form1 = '(a,t12,"=",f11.6)'
            ids(1) = get_nuclide_index('c12')
            ids(2) = get_nuclide_index('fe56')
            X = 0.5d0
            call composition_info(num_species, ids, X, xh, xhe, xz, abar, zbar, z2bar, z53bar, ye, 
     >         mass_correction, sumx, dabar_dx, dzbar_dx, dmc_dx)
            write (*,*) 'test moments of composition and the mass correction factor'
            write (*,*) 'for a C12-Fe56 mixture (both 50% by mass)'
            write (*,*)
            write (*,form1) 'Abar',abar
            write (*,form1) 'Zbar',zbar
            write (*,form1) 'Z2bar',z2bar
            write (*,form1) 'Ye',ye
            write (*,form1) 'sum(X*W/A)',mass_correction
         end subroutine do_test_composition_info
      
      subroutine do_test_Qtotal
      
         integer, parameter :: num_in = 2, num_out = 3
         integer :: reactants(num_in + num_out)
         real(dp) :: Qtotal
         
         include 'formats'
         
         write(*,*) 'test Qtotal for 2 he3 => 2 h1 + he4'
         write(*,*)
         
         reactants(1) = ihe3
         reactants(2) = ihe3
         reactants(3) = ihe4
         reactants(4) = ih1
         reactants(5) = ih1
         Qtotal = reaction_Qtotal(num_in,num_out,reactants,chem_isos)
         
         write(*,1) 'Qtotal', Qtotal   ! expect 12.86
         write(*,*)
      
      end subroutine do_test_Qtotal
      
      subroutine do_test_lodders
            integer :: ierr, i
            real(dp) :: percent

            write (*,*)
            write (*,'(a,/,72("="))') 'output of solar abundances: compare with Lodders (2003) table'
            write (*,'(a7,tr3,a11)') 'isotope','% abundance'
            do i = 1, size(namsol)
               percent = lodders03_element_atom_percent(namsol(i))
               write (*,'(a7,tr3,f11.6)') namsol(i), percent
            end do

         end subroutine do_test_lodders
         
      subroutine write_chem_ids_file
         integer, parameter :: iounit = 33
         integer :: ierr
         ierr = 0
         open(unit=iounit, file=trim('chem_ids.list'), action='write', status='replace', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open file for write_chem_ids'
            call mesa_error(__FILE__,__LINE__)
         end if
         call write_chem_ids(iounit)
         close(iounit)
      end subroutine write_chem_ids_file
      

      subroutine write_chem_ids(iounit)
         integer, intent(in) :: iounit
         integer :: i
         do i = 1, num_chem_isos
            write(iounit,'(5x,i5,3x,a5)') i, chem_isos% name(i)
         end do
         write(iounit,*)
      end subroutine write_chem_ids
      

      subroutine write_element_ids
         integer :: i
         do i = 0, max_el_z
            write(*,*) el_name(i), el_long_name(i)
         end do
      end subroutine write_element_ids
      
      
      end module mod_test_chem
      
      
      
      program test_chem
      use mod_test_chem
      !call write_element_ids; stop
      call do_test
      end program




