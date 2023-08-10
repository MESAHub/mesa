! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
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

module kap_def
  use const_def, only: dp, strlen

  implicit none

  ! interfaces for procedure pointers
  abstract interface

     subroutine other_elect_cond_opacity_interface( &
        handle, &
        zbar, logRho, logT, &
        kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
        use const_def, only: dp
        integer, intent(in) :: handle ! kap handle; from star, pass s% kap_handle
        real(dp), intent(in) :: zbar ! average ionic charge (for electron conduction)
        real(dp), intent(in) :: logRho ! the density
        real(dp), intent(in) :: logT ! the temperature
        real(dp), intent(out) :: kap ! electron conduction opacity
        real(dp), intent(out) :: dlnkap_dlnRho ! partial derivative at constant T
        real(dp), intent(out) :: dlnkap_dlnT   ! partial derivative at constant Rho
        integer, intent(out) :: ierr ! 0 means AOK.
     end subroutine other_elect_cond_opacity_interface

     subroutine other_compton_opacity_interface( &
        handle, &
        Rho, T, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
        eta, d_eta_dlnRho, d_eta_dlnT, &
        kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
        use const_def, only: dp
        integer, intent(in) :: handle ! kap handle; from star, pass s% kap_handle
        real(dp), intent(in) :: Rho, T
        real(dp), intent(in) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
        ! free_e := total combined number per nucleon of free electrons and positrons
        real(dp), intent(in) :: eta, d_eta_dlnRho, d_eta_dlnT
        ! eta := electron degeneracy parameter from eos
        real(dp), intent(out) :: kap ! electron conduction opacity
        real(dp), intent(out) :: dlnkap_dlnRho, dlnkap_dlnT
        integer, intent(out) :: ierr ! 0 means AOK.
     end subroutine other_compton_opacity_interface

     subroutine other_radiative_opacity_interface( &
        handle, &
        X, Z, XC, XN, XO, XNe, logRho, logT, &
        frac_lowT, frac_highT, frac_Type2, kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
        use const_def, only: dp
        integer, intent(in) :: handle ! kap handle; from star, pass s% kap_handle
        real(dp), intent(in) :: X, Z, XC, XN, XO, XNe ! composition
        real(dp), intent(in) :: logRho ! density
        real(dp), intent(in) :: logT ! temperature
        real(dp), intent(out) :: frac_lowT, frac_highT, frac_Type2
        real(dp), intent(out) :: kap ! opacity
        real(dp), intent(out) :: dlnkap_dlnRho ! partial derivative at constant T
        real(dp), intent(out) :: dlnkap_dlnT   ! partial derivative at constant Rho
        integer, intent(out) :: ierr ! 0 means AOK.
     end subroutine other_radiative_opacity_interface

  end interface

  
  logical, parameter :: show_allocations = .false.  ! for debugging memory usage

  ! for kap output
  integer, parameter :: num_kap_fracs = 4
  integer, parameter :: i_frac_lowT = 1
  integer, parameter :: i_frac_highT = i_frac_lowT + 1
  integer, parameter :: i_frac_Type2 = i_frac_highT + 1
  integer, parameter :: i_frac_Compton = i_frac_Type2 + 1

  ! info about op_mono elements
  integer, parameter :: num_op_mono_elements = 17
  integer :: op_mono_element_Z(num_op_mono_elements) 
  character(len=2) :: op_mono_element_name(num_op_mono_elements) 
  real(dp) :: op_mono_element_mass(num_op_mono_elements)


  integer, parameter :: kap_table_fixed_metal_form = 1
  integer, parameter :: kap_table_co_enhanced_form = 2


  ! for fixed metal tables (no enhancements)
  type Kap_Z_Table ! holds pointers to all the X tables for a particular Z
     logical :: lowT_flag
     real(dp) :: Z
     integer :: num_Xs ! number of X's for this Z
     type (Kap_X_Table), dimension(:), pointer :: x_tables  => null()! in order of increasing X
  end type Kap_Z_Table


  ! note:  logR = logRho - 3*logT + 18

  type Kap_X_Table
     logical :: not_loaded_yet
     real(dp) :: X
     real(dp) :: Z
     real(dp) :: logR_min
     real(dp) :: logR_max
     integer :: num_logRs
     integer :: ili_logRs ! =1 if logRs are evenly spaced
     real(dp), pointer :: logRs(:) => null() ! indexed from 1 to num_logRs
     real(dp) :: logT_min
     real(dp) :: logT_max
     integer :: num_logTs
     integer :: ili_logTs ! =1 if logTs are evenly spaced
     real(dp), pointer :: logTs(:) => null() ! indexed from 1 to num_logTs
     real(dp), pointer :: kap1(:) => null()
  end type Kap_X_Table


  ! for C/O enhanced tables
  type Kap_CO_Z_Table
     real(dp) :: Zbase, Zfrac_C, Zfrac_N, Zfrac_O, Zfrac_Ne
     real(dp) :: log10_Zbase ! log10(Zbase)
     type (Kap_CO_X_Table), dimension(:), pointer :: x_tables  => null()! stored in order of increasing X
     ! the X tables need not be equally spaced
  end type Kap_CO_Z_Table

  integer, parameter :: num_kap_CO_dXs = 8
  real(dp), parameter :: kap_CO_dXs(num_kap_CO_dXs) =  &
       [ 0.00d0, 0.01d0, 0.03d0, 0.10d0, 0.20d0, 0.40d0, 0.60d0, 1.0d0 ]

  type Kap_CO_Table
     integer :: table_num ! the table number from the data file
     real(dp) :: X
     real(dp) :: Z
     real(dp) :: dXC
     real(dp) :: dXO
     real(dp) :: dXC_lookup
     real(dp) :: dXO_lookup
     real(dp), dimension(:), pointer :: kap1  => null()
  end type Kap_CO_Table


  integer, parameter :: max_num_CO_tables = 70

  ! standard number of CO tables for each X+Z combo
  !           X  
  ! Z         0.0   0.1   0.35  0.7
  ! 0.0       58    58    51    32
  ! 0.001     58    58    51    30
  ! 0.004     58    58    51    30
  ! 0.010     58    58    51    30
  ! 0.020     58    58    51    30
  ! 0.030     58    58    49    30
  ! 0.100     58    53    43    26

  ! 1362 tables total



  type Kap_CO_X_Table

     logical :: not_loaded_yet
     real(dp) :: X
     real(dp) :: Z
     real(dp) :: logR_min
     real(dp) :: logR_max
     integer :: num_logRs
     integer :: ili_logRs ! =1 if logRs are evenly spaced
     real(dp), dimension(:), pointer :: logRs => null() ! indexed from 1 to num_logRs
     real(dp) :: logT_min
     real(dp) :: logT_max
     integer :: num_logTs
     integer :: ili_logTs ! =1 if logTs are evenly spaced
     real(dp), dimension(:), pointer :: logTs => null() ! indexed from 1 to num_logTs

     integer :: num_CO_tables
     ! the tables are in 3 groups
     ! 1) tables with dXC > dXO, ordered by increasing dXO, and by increasing dXC within same dXO.
     ! 2) tables with dXC = dXO, ordered by increasing value.
     ! 3) tables with dXC < dXO, ordered by increasing dXC, and by increasing dXO within same dXC.
     ! the spacing of dXC's is the same as dXO's, so there are as many tables in 3) as in 1).
     integer :: num_dXC_gt_dXO ! the number of tables with dXC > dXO
     integer :: CO_table_numbers(num_kap_CO_dXs,num_kap_CO_dXs) 
     ! entry (i,j) is the co_index for table with dXC=Xs(i) and dXO=Xs(j), or -1 if no such table.
     integer :: next_dXO_table(max_num_CO_tables) 
     ! entry (i) is the co_index for the table with same dXC and next larger dXO, or -1 if none such.
     integer :: next_dXC_table(max_num_CO_tables) 
     ! entry (i) is the co_index for the table with same dXO and next larger dXC, or -1 if none such.
     type (Kap_CO_Table), dimension(:), pointer :: co_tables => null()

  end type Kap_CO_X_Table

  type Kap_General_Info

      real(dp) :: Zbase
     
      integer :: kap_option, kap_CO_option, kap_lowT_option
      
      ! blending in T is done between the following limits
      real(dp) :: kap_blend_logT_upper_bdy ! = 3.88d0 ! old value was 4.1d0
      real(dp) :: kap_blend_logT_lower_bdy ! = 3.80d0 ! old value was 4.0d0
      ! last time I looked, the table bottom for the higher T tables was logT = 3.75
      ! while max logT for the lower T Freeman tables was 4.5
      ! so for those, you need to keep kap_blend_logT_upper_bdy < 4.5
      ! and kap_blend_logT_lower_bdy > 3.75 
      ! it is also probably a good idea to keep the blend away from H ionization
      ! logT upper of about 3.9 or a bit less will do that.

      logical :: cubic_interpolation_in_X
      logical :: cubic_interpolation_in_Z

      ! conductive opacities
      logical :: include_electron_conduction
      logical :: use_blouin_conductive_opacities

      logical :: use_Zbase_for_Type1
      logical :: use_Type2_opacities

      ! switch to Type1 if X too large
      real(dp) :: kap_Type2_full_off_X ! Type2 full off for X >= this
      real(dp) :: kap_Type2_full_on_X ! Type2 full on for X <= this

      ! switch to Type1 if dZ too small (dZ = Z - Zbase)
      real(dp) :: kap_Type2_full_off_dZ ! Type2 is full off for dZ <= this
      real(dp) :: kap_Type2_full_on_dZ ! Type2 can be full on for dZ >= this

      real(dp) :: logT_Compton_blend_hi, logR_Compton_blend_lo

      logical :: show_info

      ! debugging info
      logical :: dbg
      real(dp) :: logT_lo, logT_hi
      real(dp) :: logRho_lo, logRho_hi
      real(dp) :: X_lo, X_hi
      real(dp) :: Z_lo, Z_hi

      ! bookkeeping
      integer :: handle
      logical :: in_use

      ! User supplied inputs
      real(dp) :: kap_ctrl(10)
      integer :: kap_integer_ctrl(10)
      logical :: kap_logical_ctrl(10)
      character(len=strlen) :: kap_character_ctrl(10)

      ! other hooks

      logical :: use_other_elect_cond_opacity
      procedure (other_elect_cond_opacity_interface), pointer, nopass :: &
         other_elect_cond_opacity => null()

      logical :: use_other_compton_opacity
      procedure (other_compton_opacity_interface), pointer, nopass :: &
         other_compton_opacity => null()

      logical :: use_other_radiative_opacity
      procedure (other_radiative_opacity_interface), pointer, nopass :: &
         other_radiative_opacity => null()
      
  end type Kap_General_Info

  ! kap_options
  integer, parameter :: &
     kap_gn93 = 1, & 
     kap_gs98 = 2, & 
     kap_a09 = 3, & 
     kap_OP_gs98 = 4, & 
     kap_OP_a09 = 5, &
     kap_user = 6, &
     kap_test = 7, &
     kap_options_max = 7


  integer, parameter :: kap_max_dim = 20 !change this to make even larger grids in X and/or Z

  integer, dimension(kap_options_max) :: num_kap_Xs = 0
  real(dp), dimension(kap_max_dim, kap_options_max) :: kap_Xs = -1d0

  integer, dimension(kap_options_max) :: num_kap_Zs = 0
  real(dp), dimension(kap_max_dim, kap_options_max) :: kap_Zs= -1d0

  integer, dimension(kap_max_dim, kap_options_max) :: num_kap_Xs_for_this_Z = 0


  ! kap_CO_options
  integer, parameter :: &
     kap_CO_gn93 = 1, & 
     kap_CO_gs98 = 2, & 
     kap_CO_a09 = 3, &
     kap_CO_user = 4, &
     kap_CO_test = 5, &
     kap_CO_options_max = 5


  integer, dimension(kap_CO_options_max) :: num_kap_CO_Xs = 0
  real(dp), dimension(kap_max_dim, kap_CO_options_max) :: kap_CO_Xs = -1d0

  integer, dimension(kap_CO_options_max) :: num_kap_CO_Zs = 0
  real(dp), dimension(kap_max_dim, kap_CO_options_max) :: kap_CO_Zs= -1d0

  integer, dimension(kap_max_dim, kap_CO_options_max) :: num_kap_CO_Xs_for_this_Z = 0

  
  ! kap_lowT_options
  integer, parameter :: &
     kap_lowT_Freedman11 = 1, & 
     kap_lowT_fa05_gs98 = 2, & 
     kap_lowT_fa05_gn93 = 3, & 
     kap_lowT_fa05_a09p = 4, & 
     kap_lowT_af94_gn93 = 5, &
     kap_lowT_rt14_ag89 = 6, &
     kap_lowT_kapCN = 7, & 
     kap_lowT_AESOPUS = 8, &
     kap_lowT_user = 9, &
     kap_lowT_test = 10, &
     kap_lowT_options_max = 10

  
  integer, dimension(kap_lowT_options_max) :: num_kap_lowT_Xs = 0
  real(dp), dimension(kap_max_dim, kap_lowT_options_max) :: kap_lowT_Xs = -1d0

  integer, dimension(kap_lowT_options_max) :: num_kap_lowT_Zs = 0
  real(dp), dimension(kap_max_dim, kap_lowT_options_max) :: kap_lowT_Zs= -1d0

  integer, dimension(kap_max_dim, kap_lowT_options_max) :: num_kap_lowT_Xs_for_this_Z = 0

  
  character (len=256) :: &
     kap_option_str(kap_options_max), &
     kap_CO_option_str(kap_CO_options_max), &
     kap_lowT_option_str(kap_lowT_options_max)

  type Kap_Z_Table_Array ! in order of increasing Z
      type (Kap_Z_Table), dimension(:), pointer :: ar
  end type Kap_Z_Table_Array

  type Kap_CO_Z_Table_Array ! in order of increasing Z
      type (Kap_CO_Z_Table), dimension(:), pointer :: ar
  end type Kap_CO_Z_Table_Array

  type (Kap_Z_Table_Array), dimension(kap_options_max) :: kap_z_tables
  type (Kap_Z_Table_Array), dimension(kap_lowT_options_max) :: kap_lowT_z_tables
  type (Kap_CO_Z_Table_Array), dimension(kap_CO_options_max) :: kap_co_z_tables


  ! NOTE: in the following, "log" means base 10, "ln" means natural log, and units are cgs.

  integer, parameter :: sz_per_kap_point = 4
  !
  ! function f(x,y) with samples f(i,j) has bicubic spline fit s(x,y).
  ! compact representation of spline fit uses 4 entries as follows:
  !
  ! d(1,i,j) = s(i,j)
  ! d(2,i,j) = d2s_dx2(i,j)
  ! d(3,i,j) = d2s_dy2(i,j)
  ! d(4,i,j) = d4s_dx2_dy2(i,j)
  ! 
  ! given f(i,j), the spline fitting code can compute the other entries
  !
  ! given d(1:4,i,j), spline interpolation code can compute s(x,y)
  ! and also the partials ds_dx(x,y) and ds_dy(x,y)
  !


  logical :: kap_is_initialized = .false.

  character (len=1000) :: kap_dir, kap_cache_dir, kap_temp_cache_dir
  logical :: kap_use_cache = .true.
  logical :: kap_read_after_write_cache = .true.
  logical :: clip_to_kap_table_boundaries = .true. ! typically, this should be set true.
   ! if this is set true, then temperature and density args are
   ! clipped to the boundaries of the table.
  real(dp) :: kap_min_logRho = -40d0
   ! below this, clip logRho and set partials wrt logRho to zero


  integer, parameter :: max_kap_handles = 10
  type (Kap_General_Info), target :: kap_handles(max_kap_handles)

  !for kapCN
  logical, parameter :: kapCN_use_cache = .true.
  integer, parameter :: num_kapCN_Xs =    3
  integer, parameter :: num_kapCN_Zs =   14
  integer, parameter :: num_kapCN_fCs =   7
  integer, parameter :: num_kapCN_fNs =   3
  integer, parameter :: kapCN_num_logT = 18
  integer, parameter :: kapCN_num_logR = 17
  integer, parameter :: kapCN_tbl_size = kapCN_num_logR*kapCN_num_logT           ! 306
  integer, parameter :: kapCN_num_tbl = num_kapCN_Xs*num_kapCN_fCs*num_kapCN_fNs !  63
  
  real(dp), target :: kapCN_Z(num_kapCN_Zs)
  real(dp), target :: kapCN_fN(num_kapCN_fNs,num_kapCN_Zs)
  real(dp), target :: kapCN_fC(num_kapCN_fCs,num_kapCN_Zs)
  real(dp), target :: kapCN_X(num_kapCN_Xs)  = [ 0.5d0, 0.7d0, 0.8d0 ]
  real(dp), target :: kapCN_logT(kapCN_num_logT), kapCN_logR(kapCN_num_logR)
  real(dp) :: kapCN_min_logR, kapCN_max_logR, kapCN_min_logT, kapCN_max_logT
  logical :: kapCN_is_initialized = .false.

  real(dp), parameter :: kapCN_ZC=0.1644d0, kapCN_ZN=0.0532d0

  !stores a table for one {Z,X,fC,fN}
  type KapCN_Table
     real(dp) :: X
     real(dp) :: Z, Zbase
     real(dp) :: fC
     real(dp) :: fN
     real(dp) :: falpha
     real(dp), pointer :: kap(:) => null()
  end type KapCN_Table

  !stores a set of tables for one Z
  type KapCN_Set
     character(len=13) :: filename
     logical :: not_loaded_yet
     real(dp) :: Zbase
     real(dp) :: fC(num_kapCN_fCs)
     real(dp) :: fN(num_kapCN_fNs)
     type(KapCN_Table) :: t(kapCN_num_tbl)
  end type KapCN_Set
  type(KapCN_Set), target :: kCN(num_kapCN_Zs)


  logical :: kap_aesopus_is_initialized = .false.
  character(len=256) :: aesopus_filename

  ! stores a table for one {Z,X,fCO,fC,fN}
  type AESOPUS_Table
     real(dp) :: X
     real(dp) :: Z, Zbase
     real(dp) :: fCO
     real(dp) :: fC
     real(dp) :: fN

     ! this has to be a pointer because of the interp2d routines
     real(dp), dimension(:), pointer :: kap => null()

  end type AESOPUS_Table

  ! stores a set of tables for one Z
  type AESOPUS_TableSet

     integer :: num_Xs
     integer :: num_fCOs
     integer :: num_fCs
     integer :: num_fNs

     real(dp), dimension(:), allocatable :: Xs, fCOs, fCs, fNs

     type(AESOPUS_Table), dimension(:,:,:,:), allocatable :: t
     
  end type AESOPUS_TableSet

  
  type kapAESOPUS

     integer :: num_logTs
     integer :: num_logRs

     real(dp), dimension(:), allocatable :: logTs(:), logRs(:)

     real(dp) :: min_logR, max_logR
     real(dp) :: min_logT, max_logT
     
     integer :: num_Zs
     real(dp), dimension(:), allocatable :: Zs, logZs

     real(dp) :: Zsun
     real(dp) :: fCO_ref, fC_ref, fN_ref

     type(AESOPUS_TableSet), dimension(:), allocatable :: ts

  end type kapAESOPUS

  type(kapAESOPUS) :: kA

  logical :: kap_test_partials
  real(dp) :: kap_test_partials_val, kap_test_partials_dval_dx

contains


  subroutine kap_def_init(kap_cache_dir_in)
    use chem_def
    use utils_lib, only : mkdir
    use const_def, only: mesa_data_dir, mesa_caches_dir, mesa_temp_caches_dir, use_mesa_temp_cache
    character (*), intent(in) :: kap_cache_dir_in
    integer :: ierr, i
    
    kap_test_partials = .false.
    
    if (len_trim(kap_cache_dir_in) > 0) then
       kap_cache_dir = kap_cache_dir_in
    else if (len_trim(mesa_caches_dir) > 0) then
       kap_cache_dir = trim(mesa_caches_dir) // '/kap_cache'
    else
       kap_cache_dir = trim(mesa_data_dir) // '/kap_data/cache'
    end if
    call mkdir(kap_cache_dir)

    clip_to_kap_table_boundaries = .true.

    do i=1,max_kap_handles
       kap_handles(i)% handle = i
       kap_handles(i)% in_use = .false.
    end do

    op_mono_element_Z(1:num_op_mono_elements) = [ &
         1, 2, 6, 7, 8, 10, 11, 12, 13, 14, 16, 18, 20, 24, 25, 26, 28 ]
    op_mono_element_name(1:num_op_mono_elements) = [ &
         'h ', 'he', 'c ', 'n ', 'o ', 'ne', 'na', 'mg', 'al', &
         'si', 's ', 'ar', 'ca', 'cr', 'mn', 'fe', 'ni' ]
    op_mono_element_mass(1:num_op_mono_elements) = [ &
         1.0080d0, 4.0026d0, 12.0111d0, 14.0067d0, 15.9994d0, 20.179d0, &
         22.9898d0, 24.305d0, 26.9815d0, 28.086d0, 32.06d0, 39.948d0, &
         40.08d0, 51.996d0, 54.9380d0, 55.847d0, 58.71d0 ]

    kap_temp_cache_dir=trim(mesa_temp_caches_dir)//'/kap_cache'
    if(use_mesa_temp_cache) call mkdir(kap_temp_cache_dir)
    
    kap_option_str(kap_gn93) = 'gn93'
    kap_option_str(kap_gs98) = 'gs98'
    kap_option_str(kap_a09) = 'a09'
    kap_option_str(kap_OP_gs98) = 'OP_gs98'
    kap_option_str(kap_OP_a09) = 'OP_a09_nans_removed_by_hand'
    kap_option_str(kap_test) = 'test'
    
    kap_CO_option_str(kap_CO_gn93) = 'gn93_co'
    kap_CO_option_str(kap_CO_gs98) = 'gs98_co'
    kap_CO_option_str(kap_CO_a09) = 'a09_co'
    kap_CO_option_str(kap_CO_test) = 'test_co'
    
    kap_lowT_option_str(kap_lowT_Freedman11) = 'lowT_Freedman11'
    kap_lowT_option_str(kap_lowT_fa05_gs98) = 'lowT_fa05_gs98'
    kap_lowT_option_str(kap_lowT_fa05_gn93) = 'lowT_fa05_gn93'
    kap_lowT_option_str(kap_lowT_fa05_a09p) = 'lowT_fa05_a09p'
    kap_lowT_option_str(kap_lowT_af94_gn93) = 'lowT_af94_gn93'
    kap_lowT_option_str(kap_lowT_rt14_ag89) = 'lowT_rt14_ag89'
    kap_lowT_option_str(kap_lowT_kapCN) = 'kapCN'
    kap_lowT_option_str(kap_lowT_AESOPUS) = 'AESOPUS'
    kap_lowT_option_str(kap_lowT_test) = 'lowT_test'
    
    do i=1,kap_options_max
      nullify(kap_z_tables(i)% ar)
    end do
    do i=1,kap_lowT_options_max
      nullify(kap_lowT_z_tables(i)% ar)
    end do
    do i=1,kap_CO_options_max
      nullify(kap_co_z_tables(i)% ar)
    end do

    do i=1, kap_options_max
       select case (i)
       case DEFAULT
          num_kap_Xs(i) = 10
          kap_Xs(1:num_kap_Xs(i), i) = [0.00d0, 0.10d0, 0.20d0, &
             0.35d0, 0.50d0, 0.70d0, 0.80d0, 0.90d0, 0.95d0, 1d0]
          num_kap_Zs(i) = 13
          kap_Zs(1:num_kap_Zs(i), i) = [0.000d0, 0.0001d0, 0.0003d0, &
             0.001d0, 0.002d0, 0.004d0, 0.01d0, 0.02d0, 0.03d0, &
             0.04d0, 0.06d0, 0.08d0, 0.100d0 ]
          num_kap_Xs_for_this_Z(1:num_kap_Zs(i), i) = [ num_kap_Xs(i), num_kap_Xs(i), num_kap_Xs(i), &
             num_kap_Xs(i), num_kap_Xs(i), num_kap_Xs(i), num_kap_Xs(i), num_kap_Xs(i), num_kap_Xs(i), &
             num_kap_Xs(i), num_kap_Xs(i)-1, num_kap_Xs(i)-1, num_kap_Xs(i)-2 ]
       end select
    end do

    do i=1, kap_lowT_options_max
       select case (i)
       case(kap_lowT_Freedman11)
          num_kap_lowT_Xs(i) = 1
          kap_lowT_Xs(1:num_kap_lowT_Xs(i), i) = [ 0.00d0 ]
          num_kap_lowT_Zs(i) = 7
          kap_lowT_Zs(1:num_kap_lowT_Zs(i), i) = &
             [ 0.01d0, 0.02d0, 0.04d0, 0.100d0, 0.200d0, 0.63d0, 1.00d0 ]
          num_kap_lowT_Xs_for_this_Z(1:num_kap_lowT_Zs(i), i) = num_kap_lowT_Xs(i)
       case DEFAULT
          num_kap_lowT_Xs(i) = 10
          kap_lowT_Xs(1:num_kap_lowT_Xs(i), i) = [0.00d0, 0.10d0, 0.20d0, &
             0.35d0, 0.50d0, 0.70d0, 0.80d0, 0.90d0, 0.95d0, 1d0]
          num_kap_lowT_Zs(i) = 13
          kap_lowT_Zs(1:num_kap_lowT_Zs(i), i) = [0.000d0, 0.0001d0, 0.0003d0, &
             0.001d0, 0.002d0, 0.004d0, 0.01d0, 0.02d0, 0.03d0, &
             0.04d0, 0.06d0, 0.08d0, 0.100d0 ]
          num_kap_lowT_Xs_for_this_Z(1:num_kap_lowT_Zs(i), i) = [ num_kap_lowT_Xs(i), num_kap_lowT_Xs(i), num_kap_lowT_Xs(i), &
             num_kap_lowT_Xs(i), num_kap_lowT_Xs(i), num_kap_lowT_Xs(i), num_kap_lowT_Xs(i), num_kap_lowT_Xs(i), num_kap_lowT_Xs(i), &
             num_kap_lowT_Xs(i), num_kap_lowT_Xs(i)-1, num_kap_lowT_Xs(i)-1, num_kap_lowT_Xs(i)-2 ]
       end select
    end do

    do i=1, kap_CO_options_max
       select case (i)
       case DEFAULT
          num_kap_CO_Xs(i) = 5
          kap_CO_Xs(1:num_kap_CO_Xs(i), i) =  &
             [ 0.00d0, 0.03d0, 0.10d0, 0.35d0, 0.70d0 ]
          num_kap_CO_Zs(i) = 8
          kap_CO_Zs(1:num_kap_CO_Zs(i),i) =  &
             [ 0.000d0, 0.001d0, 0.004d0, 0.010d0, 0.020d0, 0.030d0, 0.050d0, 0.100d0 ]
          num_kap_CO_Xs_for_this_Z(1:num_kap_CO_Zs(i), i) = num_kap_CO_Xs(i)
       end select
    end do

    
  end subroutine kap_def_init



  integer function do_alloc_kap(ierr)
    integer, intent(out) :: ierr
    integer :: i
    ierr = 0
    do_alloc_kap = -1
    !$omp critical (kap_handle)
    do i = 1, max_kap_handles
       if (.not. kap_handles(i)% in_use) then
          kap_handles(i)% in_use = .true.
          do_alloc_kap = i
          exit
       end if
    end do
    !$omp end critical (kap_handle)
    if (do_alloc_kap == -1) then
       ierr = -1
       return
    end if
    if (kap_handles(do_alloc_kap)% handle /= do_alloc_kap) then
       ierr = -1
       return
    end if
  end function do_alloc_kap


  subroutine do_free_kap(handle)
    integer, intent(in) :: handle
    if (handle >= 1 .and. handle <= max_kap_handles) &
       kap_handles(handle)% in_use = .false.
  end subroutine do_free_kap


  subroutine get_kap_ptr(handle,rq,ierr)
    integer, intent(in) :: handle
    type (Kap_General_Info), pointer :: rq
    integer, intent(out):: ierr         
    if (handle < 1 .or. handle > max_kap_handles) then
       ierr = -1
       return
    end if
    rq => kap_handles(handle)
    ierr = 0
  end subroutine get_kap_ptr


  subroutine get_output_string(x,xstr,ierr) !works with X and Z
    real(dp), intent(in) :: x
    character(len=*), intent(out) :: xstr
    integer, intent(out) :: ierr
    character(len=1), parameter :: c(9)=['1','2','3','4','5','6','7','8','9']
    character(len=8) :: str
    integer :: i, j, k

    if(X < 0d0.or.X>1d0) then
       xstr='bad'
       ierr=-1
       return     
    endif
    ierr=0
    write(str,'(f8.6)') X
    k=0
    do i=1,9
       j=index(str,c(i),back=.true.)
       k=max(k,j)
    enddo
    xstr=str(1:max(k,3))
  end subroutine get_output_string
  

  subroutine do_Free_Kap_Tables
    integer :: i
    
    do i=1,kap_options_max
      if (associated(kap_z_tables(i)% ar)) &
        call free_z_tables(kap_z_tables(i)% ar)
    end do
    do i=1,kap_lowT_options_max
      if (associated(kap_lowT_z_tables(i)% ar)) &
        call free_z_tables(kap_lowT_z_tables(i)% ar)
    end do
    do i=1,kap_CO_options_max
      if (associated(kap_co_z_tables(i)% ar)) &
        call free_co_z_tables(kap_co_z_tables(i)% ar)
    end do

    contains

    subroutine free_z_tables(z_tables)
      type (Kap_Z_Table), dimension(:), pointer :: z_tables
      integer :: num_Zs
      integer :: iz
      num_Zs = size(z_tables,dim=1)
      do iz = 1, num_Zs
         if (associated(z_tables(iz)% x_tables)) &
            call free_x_tables(z_tables(iz)% x_tables, z_tables(iz)% num_Xs)
      end do
      deallocate(z_tables)
      nullify(z_tables) 
    end subroutine free_z_tables

    subroutine free_x_tables(x_tables, num_Xs)
      type (Kap_X_Table), dimension(:), pointer :: x_tables       
      integer, intent(in) :: num_Xs      
      integer :: ix
      do ix = 1, num_Xs  
         if (associated(x_tables(ix)% logRs)) then
          deallocate(x_tables(ix)% logRs)
          nullify(x_tables(ix)% logRs)
         end if
         if (associated(x_tables(ix)% logTs)) then
          deallocate(x_tables(ix)% logTs)       
          nullify(x_tables(ix)% logTs)    
         end if
         if (associated(x_tables(ix)% kap1)) then
          deallocate(x_tables(ix)% kap1)
          nullify(x_tables(ix)% kap1)
         end if
      end do
      if (associated(x_tables)) then
        deallocate(x_tables)
        nullify(x_tables)
      end if
    end subroutine free_x_tables

    subroutine free_co_z_tables(co_z_tables)
      type (Kap_CO_Z_Table), dimension(:), pointer :: co_z_tables
      integer :: num_Zs
      integer :: iz
      num_Zs = size(co_z_tables,dim=1)
      do iz = 1, num_Zs
         if (associated(co_z_tables(iz)% x_tables)) then
            call free_co_x_tables(co_z_tables(iz)% x_tables)
         end if
      end do
      deallocate(co_z_tables)
      nullify(co_z_tables) 
    end subroutine free_co_z_tables

    subroutine free_co_x_tables(x_tables)
      type (Kap_CO_X_Table), dimension(:), pointer :: x_tables 
      ! stored in order of increasing X
      integer :: num_Xs
      integer :: ix
      num_Xs = size(x_tables,dim=1)
      do ix = 1, num_Xs
         call free_co_table(x_tables(ix)% co_tables, x_tables(ix)% num_CO_tables)
         if (associated(x_tables(ix)% logRs)) then
          deallocate(x_tables(ix)% logRs)
          nullify(x_tables(ix)% logRs)
         end if
         if (associated(x_tables(ix)% logTs)) then
          deallocate(x_tables(ix)% logTs)    
          nullify(x_tables(ix)% logTs)    
         end if
      end do
      if (associated(x_tables))then
        deallocate(x_tables)    
        nullify(x_tables)
      end if
    end subroutine free_co_x_tables

    subroutine free_co_table(co_tables, num_COs)
      type (Kap_CO_Table), dimension(:), pointer :: co_tables 
      integer, intent(in) :: num_COs            
      integer :: ico            
      do ico = 1, num_COs               
         if (associated(co_tables(ico)% kap1)) then
          deallocate(co_tables(ico)% kap1)    
          nullify(co_tables(ico)% kap1) 
         end if
      end do
      if (associated(co_tables)) then
        deallocate(co_tables)      
        nullify(co_tables)
      end if
    end subroutine free_co_table

  end subroutine do_Free_Kap_Tables


end module kap_def

