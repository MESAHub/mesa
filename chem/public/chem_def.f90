! ***********************************************************************
!
!   Copyright (C) 2011-2019  Bill Paxton, Aaron Dotter, Ed Brown & The MESA Team
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

      module chem_def
      
      use utils_def, only: integer_dict
      use const_def, only: dp
      use math_lib, only: exp10
      use utils_lib, only: mesa_error
      
      implicit none

      ! Some notes on solar abundance scales:
      !
      ! published tabulations typically include both photospheric (spectroscopic)
      ! and meteoritic abundance measurements. many elements can be measured
      ! in both. some elements can only be measured one way (or the other). The
      ! following table, compiled by Frank Timmes, gives the breakdown of the
      ! elemental abundances for each of the solar abundance patterns included in
      ! MESA.
      !
      ! m = meteoritic, generally chondrites
      ! s = spectroscopic photosphere
      ! a = average of m and s
      ! r = recommended
      !
      !name    what it is    what mesa adopts
      !AG89     m,s          m, with some s
      !GS98     m,s          m, with some s
      !L03      m,s,a        as tabulated
      !AGS05    m,s          s, with some m
      !AGSS09   m,s          s, with some m
      !L09      m,s,r        r
      !A09_Prz               modified AGSS09 -  he, c, n, o, ne, mg, al, si, s, ar, fe.

      ! storage for Lodders (2003) isotopic percentages
      type isotopic_abundance_table_type
         type (integer_dict), pointer :: name_dict
         real(dp), dimension(:), allocatable :: isotopic_percent
      end type isotopic_abundance_table_type

      type(isotopic_abundance_table_type) :: lodders03_tab6


      integer, parameter :: iso_name_length = 8 ! no. characters in an nuclide name
         ! 2 for element, 3 for A, and 3 for isomeric state as '_' + 1 or 2 digit integer 0 to 99)
      integer, parameter :: long_name_length = 16 ! no. characters in a long name
      integer, parameter :: npart = 24 ! no. entries in partition fcn table
      integer, parameter :: max_el_z = 112   ! max. Z in the winvn database
      integer, parameter :: nhname = 3 ! no. isotopes for h in the winvn database

      logical, parameter :: convert_mass_excess_to_binding_energy = .true.
      integer, parameter :: nuclide_not_found = -1 ! warning flag

      integer, dimension(0:max_el_z) :: element_min_N, element_max_N 
         ! for isos included in chem_isos

      ! element names
      character(len=iso_name_length), dimension(0:max_el_z) :: &
         el_name = [character(len=iso_name_length) :: &
            'neut','h','he','li','be','b','c','n','o','f','ne',  &
            'na','mg','al','si','p','s','cl','ar','k','ca',  &
            'sc','ti','v','cr','mn','fe','co','ni','cu','zn',  &
            'ga','ge','as','se','br','kr','rb','sr','y','zr',  &
            'nb','mo','tc','ru','rh','pd','ag','cd','in','sn',  &
            'sb','te','i','xe','cs','ba','la','ce','pr','nd',  &
            'pm','sm','eu','gd','tb','dy','ho','er','tm','yb',  &
            'lu','hf','ta','w','re','os','ir','pt','au','hg',  &
            'tl','pb','bi','po','at','rn','fr','ra','ac','th', &
            'pa','u','np','pu','am','cm','bk','cf','es','fm','md', &
            'no','lr','rf','db','sg','bh','hs','mt','ds','rg','cn']
            
      character(len=long_name_length), dimension(0:max_el_z) ::  &
            el_long_name = [character(len=long_name_length) ::  &
            'neutron','hydrogen','helium','lithium','beryllium', &
            'boron','carbon','nitrogen','oxygen','fluorine','neon',  &
            'sodium','magnesium','aluminum','silicon','phosphorus', &
            'sulfur','chlorine','argon','potassium','calcium',  &
            'scandium','titanium','vanadium','chromium','manganese', &
            'iron','cobalt','nickel','copper','zinc','gallium',  &
            'germanium','arsenic','selenium','bromine','krypton', &
            'rubidium','strontium','yttrium','zirconium','niobium',  &
            'molybdenum','technetium','ruthenium','rhodium', &
            'palladium','silver','cadmium','indium','tin','antimony',  &
            'tellurium','iodine','xenon','cesium','barium', &
            'lanthanum','cerium','praseodymium','neodymium','promethium',  &
            'samarium','europium','gadolinium','terbium','dysprosium', &
            'holmium','erbium','thulium','ytterbium','lutetium',  &
            'hafnium','tantalum','tungsten','rhenium','osmium', &
            'iridium','platinum','gold','mercury','thallium','lead',  &
            'bisumth','polonium','astatine','radon','francium', &
            'radium','actinium','thorium','protactinium','uranium', &
            'neptunium','plutonium','americium','curium','berkelium', &
            'californium','einsteinium','fermium','mendelevium', &
            'nobelium','lawrencium','rutherfordium','dubnium', &
            'seaborgium','bohrium','hassium','meitnerium','darmstadtium', &
            'roentgenium','copernicum' ]

      ! aluminum isomers
      character(len=iso_name_length), dimension(2:3) :: &
         al_isomers = [character(len=iso_name_length) ::'al-6','al*6']
      character(len=long_name_length), dimension(2:3) :: &
         long_al_isomers = [character(len=long_name_length) ::  &
            'Aluminum-gs','Aluminum-ex']
      
      
      ! chem element id numbers (up to Cn)
      ! note: for isotope i, the element id number = chem_Z(i)

      !periodic table, row 1
      integer, parameter :: e_h = 1
      integer, parameter :: e_he = 2

      !periodic table, row 2
      integer, parameter :: e_li = 3
      integer, parameter :: e_be = 4
      integer, parameter :: e_b = 5
      integer, parameter :: e_c = 6
      integer, parameter :: e_n = 7
      integer, parameter :: e_o = 8
      integer, parameter :: e_f = 9
      integer, parameter :: e_ne = 10

      !periodic table, row 3
      integer, parameter :: e_na = 11
      integer, parameter :: e_mg = 12
      integer, parameter :: e_al = 13
      integer, parameter :: e_si = 14
      integer, parameter :: e_p = 15
      integer, parameter :: e_s = 16
      integer, parameter :: e_cl = 17
      integer, parameter :: e_ar = 18

      !periodic table, row 4
      integer, parameter :: e_k = 19
      integer, parameter :: e_ca = 20
      integer, parameter :: e_sc = 21
      integer, parameter :: e_ti = 22
      integer, parameter :: e_v = 23
      integer, parameter :: e_cr = 24
      integer, parameter :: e_mn = 25
      integer, parameter :: e_fe = 26
      integer, parameter :: e_co = 27
      integer, parameter :: e_ni = 28
      integer, parameter :: e_cu = 29
      integer, parameter :: e_zn = 30
      integer, parameter :: e_ga = 31
      integer, parameter :: e_ge = 32
      integer, parameter :: e_as = 33
      integer, parameter :: e_se = 34
      integer, parameter :: e_br = 35
      integer, parameter :: e_kr = 36
      
      !periodic table, row 5
      integer, parameter :: e_rb = 37
      integer, parameter :: e_sr = 38
      integer, parameter :: e_y = 39
      integer, parameter :: e_zr = 40
      integer, parameter :: e_nb = 41
      integer, parameter :: e_mo = 42
      integer, parameter :: e_tc = 43
      integer, parameter :: e_ru = 44
      integer, parameter :: e_rh = 45
      integer, parameter :: e_pd = 46
      integer, parameter :: e_ag = 47
      integer, parameter :: e_cd = 48
      integer, parameter :: e_in = 49
      integer, parameter :: e_sn = 50
      integer, parameter :: e_sb = 51
      integer, parameter :: e_te = 52
      integer, parameter :: e_i = 53
      integer, parameter :: e_xe = 54

      !periodic table, row 6
      integer, parameter :: e_cs = 55
      integer, parameter :: e_ba = 56      
      integer, parameter :: e_la = 57
      integer, parameter :: e_ce = 58
      integer, parameter :: e_pr = 59
      integer, parameter :: e_nd = 60
      integer, parameter :: e_pm = 61
      integer, parameter :: e_sm = 62
      integer, parameter :: e_eu = 63
      integer, parameter :: e_gd = 64
      integer, parameter :: e_tb = 65
      integer, parameter :: e_dy = 66
      integer, parameter :: e_ho = 67
      integer, parameter :: e_er = 68
      integer, parameter :: e_tm = 69
      integer, parameter :: e_yb = 70
      integer, parameter :: e_lu = 71
      integer, parameter :: e_hf = 72
      integer, parameter :: e_ta = 73
      integer, parameter :: e_w  = 74
      integer, parameter :: e_re = 75
      integer, parameter :: e_os = 76
      integer, parameter :: e_ir = 77
      integer, parameter :: e_pt = 78
      integer, parameter :: e_au = 79
      integer, parameter :: e_hg = 80
      integer, parameter :: e_tl = 81
      integer, parameter :: e_pb = 82
      integer, parameter :: e_bi = 83
      integer, parameter :: e_po = 84
      integer, parameter :: e_at = 85
      integer, parameter :: e_rn = 86
      
      !periodic table, row 7
      integer, parameter :: e_fr = 87
      integer, parameter :: e_ra = 88
      integer, parameter :: e_ac = 89
      integer, parameter :: e_th = 90
      integer, parameter :: e_pa = 91
      integer, parameter :: e_u  = 92
      integer, parameter :: e_np = 93
      integer, parameter :: e_pu = 94
      integer, parameter :: e_am = 95
      integer, parameter :: e_cm = 96
      integer, parameter :: e_bk = 97
      integer, parameter :: e_cf = 98
      integer, parameter :: e_es = 99
      integer, parameter :: e_fm = 100
      integer, parameter :: e_md = 101
      integer, parameter :: e_no = 102
      integer, parameter :: e_lr = 103
      integer, parameter :: e_rf = 104
      integer, parameter :: e_db = 105
      integer, parameter :: e_sg = 106
      integer, parameter :: e_bh = 107
      integer, parameter :: e_hs = 108
      integer, parameter :: e_mt = 109
      integer, parameter :: e_ds = 110
      integer, parameter :: e_rg = 111
      integer, parameter :: e_cn = 112
         
      integer, parameter :: num_chem_elements = max_el_z

      
      ! anders & grevesse 1989
      integer, parameter :: solsiz = 286
      integer, parameter :: solnamelen = 5
      character (len=solnamelen) :: namsol(solsiz)
      integer :: izsol(solsiz),iasol(solsiz),jcode(solsiz)
      real(dp) :: solx(solsiz), zsol, yesol ! according to AG89
      type (integer_dict), pointer :: Xsol_names_dict
      
      
      ! various values for current solar Z and Y (at photosphere)
      ! note that these have been reduced by diffusion from pre-MS values.
      ! values updated from Asplund et al. ARAA, 2009, 47, 481

      real(dp), parameter :: AG89_zsol   = 0.0201d0
      real(dp), parameter :: GN93_zsol   = 0.0179d0
      real(dp), parameter :: GS98_zsol   = 0.0169d0
      real(dp), parameter :: L03_zsol    = 0.0133d0
      real(dp), parameter :: AGS05_zsol  = 0.0122d0
      real(dp), parameter :: AGSS09_zsol = 0.0134d0
      real(dp), parameter :: L09_zsol    = 0.0154d0
      real(dp), parameter :: A09_Prz_zsol = 0.014d0
      
      real(dp), parameter :: AG89_ysol   = 0.2485d0
      real(dp), parameter :: GN93_ysol   = 0.2485d0
      real(dp), parameter :: GS98_ysol   = 0.2485d0
      real(dp), parameter :: L03_ysol    = 0.2377d0
      real(dp), parameter :: AGS05_ysol  = 0.2485d0
      real(dp), parameter :: AGSS09_ysol = 0.2485d0
      real(dp), parameter :: L09_ysol    = 0.2751d0
      real(dp), parameter :: A09_Prz_ysol = 0.276d0
      
      character(len=iso_name_length) :: chem_element_main_iso_name(num_chem_elements)
      integer, parameter :: chem_element_name_len = iso_name_length
      character (len=chem_element_name_len)  :: chem_element_Name(num_chem_elements) 
         ! names for elements
      
      
      ! identifiers for different Z fractions.
      integer, parameter :: AG89_zfracs = 1
      integer, parameter :: GN93_zfracs = 2
      integer, parameter :: GS98_zfracs = 3
      integer, parameter :: L03_zfracs = 4
      integer, parameter :: AGS05_zfracs = 5
      integer, parameter :: AGSS09_zfracs = 6
      integer, parameter :: L09_zfracs = 7
      integer, parameter :: A09_Prz_zfracs = 8
      
         
      real(dp) :: AG89_element_zfrac(num_chem_elements) ! fraction by mass of total Z
         ! Anders & Grevesse 1989
         
      real(dp) :: GN93_element_zfrac(num_chem_elements) ! fraction by mass of total Z
         ! Grevesse and Noels 1993 abundances
         
      real(dp) :: GS98_element_zfrac(num_chem_elements) ! fraction by mass of total Z
         ! Grevesse and Sauval 1998 abundances
         
      real(dp) :: L03_element_zfrac(num_chem_elements) ! fraction by mass of total Z
         ! Lodders 2003 abundances
         
      real(dp) :: AGS05_element_zfrac(num_chem_elements) ! fraction by mass of total Z
         ! Asplund, Grevesse, and Sauval 2005 abundances
         
      real(dp) :: AGSS09_element_zfrac(num_chem_elements) ! fraction by mass of total Z
         ! Asplund, Grevesse, Sauval, and Scott 2009 abundances
         ! Annu. Rev. Astron. Astrophys. 2009. 47:481–522
         
      real(dp) :: L09_element_zfrac(num_chem_elements) ! fraction by mass of total Z
         ! Lodders and Palme, 2009.  (http://adsabs.harvard.edu/abs/2009M%26PSA..72.5154L)
         
      real(dp) :: A09_Prz_zfrac(num_chem_elements)  ! fraction by mass of the total Z
         ! Abundances are abased on Asplund, Grevesse, Sauval, and Scott 2009, ARA&A, 47:481–522
         ! but that of the some key elements are updated based on:
         ! "Present-day cosmic abundances ..." 
            ! by Nieva, M.-F. & Przybilla, N. 2012, A&A, 539, 143 
         ! and the proceeding paper
            ! "Hot Stars and Cosmic Abundances" by
            ! Przybilla N. Nieva M. F., Irrgang A. and Butler K. 2013, EAS Publ. Ser.
            ! in preparation
         ! The modified abundances w.r.t. A09 are (eps_El = log(El/H)+12.0) 
            ! eps_He = 10.99
            ! eps_C  = 8.33
            ! eps_N  = 7.79
            ! eps_O  = 8.76
            ! eps_Ne = 8.09    ! Cunha et al. (2006) give 8.11
            ! eps_Mg = 7.56
            ! eps_Al = 6.30
            ! eps_Si = 7.50
            ! eps_S  = 7.14
            ! eps_Ar = 6.50
            ! eps_Fe = 7.52   
         

      type (integer_dict), pointer :: chem_element_names_dict

         
      real(dp) :: element_atomic_weight(num_chem_elements)
         ! de Laeter et al, Pure and Applied Chemistry 75(6), 683–799, 2003.
         ! (IUPAC Technical Report)
      


      ! temperature values at which partition function is defined
      real(dp), dimension(npart) :: Tpart
         

      ! mass excess of proton, neutron in MeV (for calculating binding energies)
      ! should be consistent with the mass excess of the prot and neut from the isotopes.data file
      ! Set in chem_isos_io.f90
      real(dp) :: del_Mp, del_Mn

      type nuclide_data
         integer :: nnuclides
         character(len=iso_name_length), dimension(:), pointer :: name  ! name of nuclide
         integer, dimension(:), pointer :: chem_id ! (nnuclides)
            ! gives chem_id for member of nuclide_data
         integer, dimension(:), pointer :: nuclide ! (num_chem_isos)
            ! gives index in nuclide_data 1 .. nnuclides or 0 if not included
         real(dp), dimension(:), pointer :: W ! atomic weight (mass in amu units)
         integer, dimension(:), pointer :: Z ! number of protons
         integer, dimension(:), pointer :: N ! number of neutrons
         integer, dimension(:), pointer :: Z_plus_N ! number of baryons
         integer, dimension(:), pointer :: isomeric_state ! 0 is default
         real(dp), dimension(:), pointer :: spin   ! ground-state spin
         real(dp), dimension(:), pointer :: binding_energy
            ! the binding energy is B = Z*del_Mp + N*del_Mn - mass_excess
         real(dp), dimension(:,:), pointer :: pfcn ! table of partition function
         real(dp), dimension(:), pointer :: mass_excess
         real(dp), dimension(:), pointer :: Z53 ! cache expensive Z^5/3 result
      end type nuclide_data      

      type (nuclide_data) :: chem_isos ! from winvn
      type (integer_dict), pointer :: chem_isos_dict
      integer :: num_chem_isos   ! no. entries in isotopes database

      ! storage container for a set of nuclides,
      ! used for extracting a subset of the full winvn database
      ! a set of nuclides is actually an array of these items
      type nuclide_set
         character(len=iso_name_length) :: nuclide
         integer :: rank
      end type nuclide_set
               
      
      ! reaction categories

      integer, parameter :: ipp = 1 ! pp chains
      integer, parameter :: icno = 2 ! cno cycles
      integer, parameter :: i3alf = 3 ! triple alpha  
      
      ! "burn" in the following means decays or captures of protons, alphas, or neutrons
      integer, parameter :: i_burn_c = 4
      integer, parameter :: i_burn_n = 5
      integer, parameter :: i_burn_o = 6
      integer, parameter :: i_burn_ne = 7
      integer, parameter :: i_burn_na = 8
      integer, parameter :: i_burn_mg = 9
      integer, parameter :: i_burn_si = 10
      integer, parameter :: i_burn_s = 11
      integer, parameter :: i_burn_ar = 12
      integer, parameter :: i_burn_ca = 13
      integer, parameter :: i_burn_ti = 14
      integer, parameter :: i_burn_cr = 15
      integer, parameter :: i_burn_fe = 16
      
      integer, parameter :: icc = 17 ! c12 + c12
      integer, parameter :: ico = 18 ! c12 + o16
      integer, parameter :: ioo = 19 ! o16 + o16
      
      integer, parameter :: ipnhe4 = 20 ! 2prot + 2neut -> he4
      
      integer, parameter :: iphoto = 21 ! photodisintegration
         ! note: for photodisintegrations, eps_nuc will be negative.
         
      integer, parameter :: i_ni56_co56 = 22 ! ni56 -> co56
      integer, parameter :: i_co56_fe56 = 23 ! co56 -> fe56
      integer, parameter :: iother = 24 ! misc.
   
      integer, parameter :: num_categories = iother

      integer, parameter :: maxlen_category_name = 16
      character (len=maxlen_category_name) :: category_name(num_categories)      
      type (integer_dict), pointer :: category_names_dict

      
      ! some commonly used values of get_nuclide_index      
      integer :: &
         ih1, ih2, ih3,  &
         ihe3, ihe4,  &
         ili6, ili7, ili8,  &
         ibe7, ibe8, ibe9, ibe10, ibe11,  &
         ib8, ib10, ib11, ib12, ib13, ib14,  &
         ic9, ic10, ic11, ic12, ic13, ic14, ic15, ic16,  &
         in12, in13, in14, in15, in16, in17, in18, in19, in20,  &
         io13, io14, io15, io16, io17, io18, io19, io20,  &
         if15, if16, if17, if18, if19, if20, if21, if22, if23, if24,  &
         ine17, ine18, ine19, ine20, ine21, ine22, ine23, ine24, ine25,  &
            ine26, ine27, ine28,  &
         ina20, ina21, ina22, ina23, ina24, ina25, ina26, ina27, ina28,  &
            ina29, ina30, ina31,  &
         img20, img21, img22, img23, img24, img25, img26, img27, img28,  &
            img29, img30, img31, img32, img33,  &
         ial22, ial23, ial24, ial25, ial26, ial27, ial28, ial29, ial30,  &
            ial31, ial32, ial33, ial34, ial35,  &
         isi22, isi23, isi24, isi25, isi26, isi27, isi28, isi29, isi30,  &
            isi31, isi32, isi33, isi34, isi35, isi36, isi37, isi38,  &
         ip26, ip27, ip28, ip29, ip30, ip31, ip32, ip33, ip34, ip35,  &
            ip36, ip37, ip38, ip39, ip40,  &
         is27, is28, is29, is30, is31, is32, is33, is34, is35, is36,  &
            is37, is38, is39, is40, is41, is42,  &
         icl31, icl32, icl33, icl34, icl35, icl36, icl37, icl38, icl39,  &
            icl40, icl41, icl42, icl43, icl44,  &
         iar31, iar32, iar33, iar34, iar35, iar36, iar37, iar38, iar39,  &
            iar40, iar41, iar42, iar43, iar44, iar45, iar46, iar47,  &
         ik35, ik36, ik37, ik38, ik39, ik40, ik41, ik42, ik43, ik44, ik45, ik46, ik47, &
         ica35, ica36, ica37, ica38, ica39, ica40, ica41, ica42, ica43,  &
            ica44, ica45, ica46, ica47, ica48, ica49, ica50, ica51, ica52, ica53,  &
         isc40, isc41, isc42, isc43, isc44, isc45, isc46, isc47, isc48,  &
            isc49, isc50, isc51, isc52, isc53,  &
         iti39, iti40, iti41, iti42, iti43, iti44, iti45, iti46, iti47,  &
            iti48, iti49, iti50, iti51, iti52, iti53, iti54, iti55,  &
         iv43, iv44, iv45, iv46, iv47, iv48, iv49, iv50, iv51, iv52, iv53,  &
            iv54, iv55, iv56, iv57,  &
         icr43, icr44, icr45, icr46, icr47, icr48, icr49, icr50, icr51,  &
            icr52, icr53, icr54, icr55, icr56, icr57, icr58, icr59, &
            icr60, icr61, icr62, icr63, icr64, icr65, icr66,  &
         imn46, imn47, imn48, imn49, imn50, imn51, imn52, imn53, imn54,  &
            imn55, imn56, imn57, imn58, imn59, imn60, imn61, imn62, imn63,  &
         ife46, ife47, ife48, ife49, ife50, ife51, ife52, ife53, ife54, ife55, ife56,  &
            ife57, ife58, ife59, ife60, ife61, ife62, ife63, ife64, ife65, ife66, ife68, &
         ico50, ico51, ico52, ico53, ico54, ico55, ico56, ico57, ico58,  &
            ico59, ico60, ico61, ico62, ico63, ico64, ico65, ico66, ico67,  &
         ini50, ini51, ini52, ini53, ini54, ini55, ini56, ini57, ini58,  &
            ini59, ini60, ini61, ini62, ini63, ini64, ini65, ini66, ini67,  &
            ini68, ini69, ini70, ini71, ini72, ini73,  &
         icu56, icu57, icu58, icu59, icu60, icu61, icu62, icu63, icu64,  &
            icu65, icu66, icu67, icu68, icu69, icu70, icu71, icu72,  &
         izn55, izn56, izn57, izn58, izn59, izn60, izn61, izn62, izn63, izn64,  &
            izn65, izn66, izn67, izn68, izn69, izn70, izn71, izn72, izn73, izn74, &
         iga60, iga61, iga62, iga63, iga64, iga65, iga66, iga67, iga68,  &
            iga69, iga70, iga71, iga72, iga73, iga74, iga75,  &
         ige59, ige60, ige61, ige62, ige63, ige64, ige65, ige66, ige67,  &
            ige68, ige69, ige70, ige71, ige72, ige73, ige74, ige75, ige76,  &
         ias71,ias72, ias73, ias74, ias75, ias76, ias77, ias78, ias79, &
         ise68, ise69, ise70, ise71, ise72, ise73, ise74, ise75, ise76,  &
            ikr70, ikr72, isr74, isr75, isr76, izr77, izr80, imo82, imo84,  &
            iru86, iru87, iru88, ipd89, ipd91, ipd92, icd93, icd96,  &
            isn98, isn100, isn102, isn104, &
         ineut, iprot


      logical :: chem_has_been_initialized = .false.


      
      contains
      
      
      subroutine init_chem_tables
         use utils_lib, only: integer_dict_define, integer_dict_create_hash
      
         integer :: i, ierr
         
         Tpart = (/  &
            0.10d0, 0.15d0, 0.20d0, 0.30d0, 0.40d0, 0.50d0, &
            0.60d0, 0.70d0, 0.80d0, 0.90d0, 1.00d0, 1.50d0, &
            2.00d0, 2.50d0, 3.00d0, 3.50d0, 4.00d0, 4.50d0, &
            5.00d0, 6.00d0, 7.00d0, 8.00d0, 9.00d0, 10.0d0 /)
            
         call init_ag_data
      
         call init_chem_element_names
      
         call init_chem_element_main_iso_names
      
         call init_element_atomic_weights
         
         call init_AG89_data
         
         call init_GN93_data
         
         call init_GS98_data
         
         call init_L03_data
         
         call init_AGS05_data
         
         call init_AGSS09_data
         
         call init_A09_Przybilla_data
         
         call init_L09_data
         
         nullify(chem_element_names_dict)
         do i=1,num_chem_elements
            call integer_dict_define(chem_element_names_dict, chem_element_Name(i), i, ierr)
            if (ierr /= 0) then
               write(*,*) 'FATAL ERROR: init_chem_tables failed in integer_dict_define'
               flush(6)
               call mesa_error(__FILE__,__LINE__)
            end if
         end do
         call integer_dict_create_hash(chem_element_names_dict, ierr)
         if (ierr /= 0) then
            write(*,*) 'FATAL ERROR: init_chem_tables failed in integer_dict_create_hash'
            flush(6)
            call mesa_error(__FILE__,__LINE__)
         end if
         
         call set_category_names
         nullify(category_names_dict)
         do i=1,num_categories
            call integer_dict_define(category_names_dict, category_name(i), i, ierr)
            if (ierr /= 0) then
               write(*,*) 'FATAL ERROR: rates_def_init failed in integer_dict_define'
               flush(6)
               call mesa_error(__FILE__,__LINE__)
            end if
         end do
         call integer_dict_create_hash(category_names_dict, ierr)
         if (ierr /= 0) then
            write(*,*) 'FATAL ERROR: rates_def_init failed in integer_dict_create_hash'
            flush(6)
            call mesa_error(__FILE__,__LINE__)
         end if
         
      end subroutine init_chem_tables
      
      
      subroutine init_ag_data
         use utils_lib
      
         real(dp) :: sum
         integer :: i, j, ierr

!..names of the stable isotopes
         namsol(1:120) = (/ &
       'h1   ','h2   ','he3  ','he4  ','li6  ','li7  ','be9  ','b10  ', &
       'b11  ','c12  ','c13  ','n14  ','n15  ','o16  ','o17  ','o18  ', &
       'f19  ','ne20 ','ne21 ','ne22 ','na23 ','mg24 ','mg25 ','mg26 ', &
       'al27 ','si28 ','si29 ','si30 ','p31  ','s32  ','s33  ','s34  ', &
       's36  ','cl35 ','cl37 ','ar36 ','ar38 ','ar40 ','k39  ','k40  ', &
       'k41  ','ca40 ','ca42 ','ca43 ','ca44 ','ca46 ','ca48 ','sc45 ', &
       'ti46 ','ti47 ','ti48 ','ti49 ','ti50 ','v50  ','v51  ','cr50 ', &
       'cr52 ','cr53 ','cr54 ','mn55 ','fe54 ','fe56 ','fe57 ','fe58 ', &
       'co59 ','ni58 ','ni60 ','ni61 ','ni62 ','ni64 ','cu63 ','cu65 ', &
       'zn64 ','zn66 ','zn67 ','zn68 ','zn70 ','ga69 ','ga71 ','ge70 ', &
       'ge72 ','ge73 ','ge74 ','ge76 ','as75 ','se74 ','se76 ','se77 ', &
       'se78 ','se80 ','se82 ','br79 ','br81 ','kr78 ','kr80 ','kr82 ', &
       'kr83 ','kr84 ','kr86 ','rb85 ','rb87 ','sr84 ','sr86 ','sr87 ', &
       'sr88 ','y89  ','zr90 ','zr91 ','zr92 ','zr94 ','zr96 ','nb93 ', &
       'mo92 ','mo94 ','mo95 ','mo96 ','mo97 ','mo98 ','mo100','ru96 '  /)

         namsol(121:240) = (/ &
       'ru98 ','ru99 ','ru100','ru101','ru102','ru104','rh103','pd102', &
       'pd104','pd105','pd106','pd108','pd110','ag107','ag109','cd106', &
       'cd108','cd110','cd111','cd112','cd113','cd114','cd116','in113', &
       'in115','sn112','sn114','sn115','sn116','sn117','sn118','sn119', &
       'sn120','sn122','sn124','sb121','sb123','te120','te122','te123', &
       'te124','te125','te126','te128','te130','i127 ','xe124','xe126', &
       'xe128','xe129','xe130','xe131','xe132','xe134','xe136','cs133', &
       'ba130','ba132','ba134','ba135','ba136','ba137','ba138','la138', &
       'la139','ce136','ce138','ce140','ce142','pr141','nd142','nd143', &
       'nd144','nd145','nd146','nd148','nd150','sm144','sm147','sm148', &
       'sm149','sm150','sm152','sm154','eu151','eu153','gd152','gd154', &
       'gd155','gd156','gd157','gd158','gd160','tb159','dy156','dy158', &
       'dy160','dy161','dy162','dy163','dy164','ho165','er162','er164', &
       'er166','er167','er168','er170','tm169','yb168','yb170','yb171', &
       'yb172','yb173','yb174','yb176','lu175','lu176','hf174','hf176' /)

         namsol(241:286) = (/ &
       'hf177','hf178','hf179','hf180','ta180','ta181','w180 ','w182 ', &
       'w183 ','w184 ','w186 ','re185','re187','os184','os186','os187', &
       'os188','os189','os190','os192','ir191','ir193','pt190','pt192', &
       'pt194','pt195','pt196','pt198','au197','hg196','hg198','hg199', &
       'hg200','hg201','hg202','hg204','tl203','tl205','pb204','pb206', &
       'pb207','pb208','bi209','th232','u235 ','u238 ' /)


!..anders & grevesse 1989 solar mass fractions
         solx(1:45) = (/ &
           7.0573D-01, 4.8010D-05, 2.9291D-05, 2.7521D-01, 6.4957D-10, &
           9.3490D-09, 1.6619D-10, 1.0674D-09, 4.7301D-09, 3.0324D-03, &
           3.6501D-05, 1.1049D-03, 4.3634D-06, 9.5918D-03, 3.8873D-06, &
           2.1673D-05, 4.0515D-07, 1.6189D-03, 4.1274D-06, 1.3022D-04, &
           3.3394D-05, 5.1480D-04, 6.7664D-05, 7.7605D-05, 5.8052D-05, &
           6.5301D-04, 3.4257D-05, 2.3524D-05, 8.1551D-06, 3.9581D-04, &
           3.2221D-06, 1.8663D-05, 9.3793D-08, 2.5320D-06, 8.5449D-07, &
           7.7402D-05, 1.5379D-05, 2.6307D-08, 3.4725D-06, 4.4519D-10, &
           2.6342D-07, 5.9898D-05, 4.1964D-07, 8.9734D-07, 1.4135D-06 /)

         solx(46:90)  = (/ &
             2.7926D-09, 1.3841D-07, 3.8929D-08, 2.2340D-07, 2.0805D-07, &
             2.1491D-06, 1.6361D-07, 1.6442D-07, 9.2579D-10, 3.7669D-07, &
             7.4240D-07, 1.4863D-05, 1.7160D-06, 4.3573D-07, 1.3286D-05, &
             7.1301D-05, 1.1686D-03, 2.8548D-05, 3.6971D-06, 3.3579D-06, &
             4.9441D-05, 1.9578D-05, 8.5944D-07, 2.7759D-06, 7.2687D-07, &
             5.7528D-07, 2.6471D-07, 9.9237D-07, 5.8765D-07, 8.7619D-08, &
             4.0593D-07, 1.3811D-08, 3.9619D-08, 2.7119D-08, 4.3204D-08, &
             5.9372D-08, 1.7136D-08, 8.1237D-08, 1.7840D-08, 1.2445D-08, &
             1.0295D-09, 1.0766D-08, 9.1542D-09, 2.9003D-08, 6.2529D-08 /)

         solx(91:135)  = (/ &
             1.1823D-08, 1.1950D-08, 1.2006D-08, 3.0187D-10, 2.0216D-09, &
             1.0682D-08, 1.0833D-08, 5.4607D-08, 1.7055D-08, 1.1008D-08, &
             4.3353D-09, 2.8047D-10, 5.0468D-09, 3.6091D-09, 4.3183D-08, &
             1.0446D-08, 1.3363D-08, 2.9463D-09, 4.5612D-09, 4.7079D-09, &
             7.7706D-10, 1.6420D-09, 8.7966D-10, 5.6114D-10, 9.7562D-10, &
             1.0320D-09, 5.9868D-10, 1.5245D-09, 6.2225D-10, 2.5012D-10, &
             8.6761D-11, 5.9099D-10, 5.9190D-10, 8.0731D-10, 1.5171D-09, &
             9.1547D-10, 8.9625D-10, 3.6637D-11, 4.0775D-10, 8.2335D-10, &
             1.0189D-09, 1.0053D-09, 4.5354D-10, 6.8205D-10, 6.4517D-10 /)

         solx(136:180)  = (/ &
             5.3893D-11, 3.9065D-11, 5.5927D-10, 5.7839D-10, 1.0992D-09, &
             5.6309D-10, 1.3351D-09, 3.5504D-10, 2.2581D-11, 5.1197D-10, &
             1.0539D-10, 7.1802D-11, 3.9852D-11, 1.6285D-09, 8.6713D-10, &
             2.7609D-09, 9.8731D-10, 3.7639D-09, 5.4622D-10, 6.9318D-10, &
             5.4174D-10, 4.1069D-10, 1.3052D-11, 3.8266D-10, 1.3316D-10, &
             7.1827D-10, 1.0814D-09, 3.1553D-09, 4.9538D-09, 5.3600D-09, &
             2.8912D-09, 1.7910D-11, 1.6223D-11, 3.3349D-10, 4.1767D-09, &
             6.7411D-10, 3.3799D-09, 4.1403D-09, 1.5558D-09, 1.2832D-09, &
             1.2515D-09, 1.5652D-11, 1.5125D-11, 3.6946D-10, 1.0108D-09 /)

         solx(181:225)  = (/ &
             1.2144D-09, 1.7466D-09, 1.1240D-08, 1.3858D-12, 1.5681D-09, &
             7.4306D-12, 9.9136D-12, 3.5767D-09, 4.5258D-10, 5.9562D-10, &
             8.0817D-10, 3.6533D-10, 7.1757D-10, 2.5198D-10, 5.2441D-10, &
             1.7857D-10, 1.7719D-10, 2.9140D-11, 1.4390D-10, 1.0931D-10, &
             1.3417D-10, 7.2470D-11, 2.6491D-10, 2.2827D-10, 1.7761D-10, &
             1.9660D-10, 2.5376D-12, 2.8008D-11, 1.9133D-10, 2.6675D-10, &
             2.0492D-10, 3.2772D-10, 2.9180D-10, 2.8274D-10, 8.6812D-13, &
             1.4787D-12, 3.7315D-11, 3.0340D-10, 4.1387D-10, 4.0489D-10, &
             4.6047D-10, 3.7104D-10, 1.4342D-12, 1.6759D-11, 3.5397D-10 /)

         solx(226:270)  = (/ &
             2.4332D-10, 2.8557D-10, 1.6082D-10, 1.6159D-10, 1.3599D-12, &
             3.2509D-11, 1.5312D-10, 2.3624D-10, 1.7504D-10, 3.4682D-10, &
             1.4023D-10, 1.5803D-10, 4.2293D-12, 1.0783D-12, 3.4992D-11, &
             1.2581D-10, 1.8550D-10, 9.3272D-11, 2.4131D-10, 1.1292D-14, &
             9.4772D-11, 7.8768D-13, 1.6113D-10, 8.7950D-11, 1.8989D-10, &
             1.7878D-10, 9.0315D-11, 1.5326D-10, 5.6782D-13, 5.0342D-11, &
             5.1086D-11, 4.2704D-10, 5.2110D-10, 8.5547D-10, 1.3453D-09, &
             1.1933D-09, 2.0211D-09, 8.1702D-13, 5.0994D-11, 2.1641D-09, &
             2.2344D-09, 1.6757D-09, 4.8231D-10, 9.3184D-10, 2.3797D-12 /)

         solx(271:286)  = (/ &
             1.7079D-10, 2.8843D-10, 3.9764D-10, 2.2828D-10, 5.1607D-10, &
             1.2023D-10, 2.7882D-10, 6.7411D-10, 3.1529D-10, 3.1369D-09, &
             3.4034D-09, 9.6809D-09, 7.6127D-10, 1.9659D-10, 3.8519D-13, &
             5.3760D-11 /)
                                                          

!..charge of the stable isotopes

         izsol(1:117)  = (/ &
         1,   1,   2,   2,   3,   3,   4,   5,   5,   6,   6,   7,   7, &
         8,   8,   8,   9,  10,  10,  10,  11,  12,  12,  12,  13,  14, &
        14,  14,  15,  16,  16,  16,  16,  17,  17,  18,  18,  18,  19, &
        19,  19,  20,  20,  20,  20,  20,  20,  21,  22,  22,  22,  22, &
        22,  23,  23,  24,  24,  24,  24,  25,  26,  26,  26,  26,  27, &
        28,  28,  28,  28,  28,  29,  29,  30,  30,  30,  30,  30,  31, &
        31,  32,  32,  32,  32,  32,  33,  34,  34,  34,  34,  34,  34, &
        35,  35,  36,  36,  36,  36,  36,  36,  37,  37,  38,  38,  38, &
        38,  39,  40,  40,  40,  40,  40,  41,  42,  42,  42,  42,  42 /)
 
         izsol(118:234)  = (/ &
        42,  42,  44,  44,  44,  44,  44,  44,  44,  45,  46,  46,  46, &
        46,  46,  46,  47,  47,  48,  48,  48,  48,  48,  48,  48,  48, &
        49,  49,  50,  50,  50,  50,  50,  50,  50,  50,  50,  50,  51, &
        51,  52,  52,  52,  52,  52,  52,  52,  52,  53,  54,  54,  54, &
        54,  54,  54,  54,  54,  54,  55,  56,  56,  56,  56,  56,  56, &
        56,  57,  57,  58,  58,  58,  58,  59,  60,  60,  60,  60,  60, &
        60,  60,  62,  62,  62,  62,  62,  62,  62,  63,  63,  64,  64, &
        64,  64,  64,  64,  64,  65,  66,  66,  66,  66,  66,  66,  66, &
        67,  68,  68,  68,  68,  68,  68,  69,  70,  70,  70,  70,  70 /)

         izsol(235:286)  = (/ &
        70,  70,  71,  71,  72,  72,  72,  72,  72,  72,  73,  73,  74, &
        74,  74,  74,  74,  75,  75,  76,  76,  76,  76,  76,  76,  76, &
        77,  77,  78,  78,  78,  78,  78,  78,  79,  80,  80,  80,  80, &
        80,  80,  80,  81,  81,  82,  82,  82,  82,  83,  90,  92,  92 /)


!..number of nucleons (protons and neutrons) in the stable isotopes

         iasol(1:117)  = (/ &
         1,   2,   3,   4,   6,   7,   9,  10,  11,  12,  13,  14,  15, &
        16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28, &
        29,  30,  31,  32,  33,  34,  36,  35,  37,  36,  38,  40,  39, &
        40,  41,  40,  42,  43,  44,  46,  48,  45,  46,  47,  48,  49, &
        50,  50,  51,  50,  52,  53,  54,  55,  54,  56,  57,  58,  59, &
        58,  60,  61,  62,  64,  63,  65,  64,  66,  67,  68,  70,  69, &
        71,  70,  72,  73,  74,  76,  75,  74,  76,  77,  78,  80,  82, &
        79,  81,  78,  80,  82,  83,  84,  86,  85,  87,  84,  86,  87, &
        88,  89,  90,  91,  92,  94,  96,  93,  92,  94,  95,  96,  97 /)

         iasol(118:234)  = (/ &
        98, 100,  96,  98,  99, 100, 101, 102, 104, 103, 102, 104, 105, &
       106, 108, 110, 107, 109, 106, 108, 110, 111, 112, 113, 114, 116, &
       113, 115, 112, 114, 115, 116, 117, 118, 119, 120, 122, 124, 121, &
       123, 120, 122, 123, 124, 125, 126, 128, 130, 127, 124, 126, 128, &
       129, 130, 131, 132, 134, 136, 133, 130, 132, 134, 135, 136, 137, &
       138, 138, 139, 136, 138, 140, 142, 141, 142, 143, 144, 145, 146, &
       148, 150, 144, 147, 148, 149, 150, 152, 154, 151, 153, 152, 154, &
       155, 156, 157, 158, 160, 159, 156, 158, 160, 161, 162, 163, 164, &
       165, 162, 164, 166, 167, 168, 170, 169, 168, 170, 171, 172, 173 /)

         iasol(235:286)  = (/ &
       174, 176, 175, 176, 174, 176, 177, 178, 179, 180, 180, 181, 180, &
       182, 183, 184, 186, 185, 187, 184, 186, 187, 188, 189, 190, 192, &
       191, 193, 190, 192, 194, 195, 196, 198, 197, 196, 198, 199, 200, &
       201, 202, 204, 203, 205, 204, 206, 207, 208, 209, 232, 235, 238 /)


! jcode tells the type progenitors each stable species can have.
! jcode = 0 if the species is the only stable one of that a
!       = 1 if the species can have proton-rich progenitors
!       = 2 if the species can have neutron-rich progenitors
!       = 3 if the species can only be made as itself (eg k40)

        jcode(1:286) = (/ &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, &
         0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, &
         0,   0,   0,   0,   0,   0,   2,   0,   0,   1,   0,   2,   0, &
         3,   0,   1,   0,   0,   0,   2,   2,   0,   1,   0,   1,   0, &
         2,   3,   0,   1,   0,   0,   2,   0,   1,   0,   0,   2,   0, &
         1,   0,   0,   0,   2,   0,   0,   1,   0,   0,   0,   2,   0, &
         0,   1,   0,   0,   2,   2,   0,   1,   1,   0,   2,   2,   2, &
         0,   0,   1,   1,   1,   0,   2,   2,   0,   2,   1,   1,   1, &
         0,   0,   0,   0,   2,   2,   2,   0,   1,   1,   0,   3,   0, &
         2,   2,   1,   1,   0,   1,   0,   2,   2,   0,   1,   1,   0, &
         2,   2,   2,   0,   0,   1,   1,   1,   0,   2,   2,   2,   2, &
         1,   2,   1,   1,   1,   1,   0,   0,   0,   2,   2,   2,   0, &
         2,   1,   1,   1,   3,   0,   2,   2,   2,   0,   1,   1,   1, &
         0,   3,   0,   2,   2,   2,   0,   1,   1,   1,   0,   3,   0, &
         2,   3,   0,   1,   1,   0,   2,   0,   1,   0,   2,   0,   0, &
         2,   2,   1,   0,   1,   0,   1,   2,   2,   0,   0,   1,   1, &
         0,   2,   0,   2,   2,   0,   1,   1,   1,   0,   2,   0,   2, &
         0,   1,   1,   0,   0,   2,   2,   0,   1,   1,   0,   0,   0, &
         2,   2,   0,   3,   1,   1,   0,   0,   0,   2,   3,   0,   1, &
         0,   0,   2,   2,   0,   2,   1,   1,   1,   0,   0,   2,   2, &
         0,   0,   1,   1,   0,   0,   2,   2,   0,   1,   1,   0,   0, &
         0,   0,   2,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0/)

      
   ! get sum and stuff residual into hydrogen
        sum = 0.0d0
        do j=1,solsiz
           sum = sum + solx(j)
        enddo
        sum = 1.0d0 - sum
        solx(1) = solx(1) + sum
        
        sum  = 0.0d0
        do j=1,solsiz
           if (izsol(j) .ge. 3) then
              sum  = sum + solx(j)
           endif
        enddo
        zsol = sum
        
        sum = 0.0d0
        do j=1,solsiz
           sum = sum + dble(izsol(j))/dble(iasol(j))*solx(j)
        enddo
        yesol = sum
        
        nullify(Xsol_names_dict)
        ierr = 0
        do i=1,solsiz
           call integer_dict_define(Xsol_names_dict, namsol(i), i, ierr)
           if (ierr /= 0) then
              write(*,*) 'FATAL ERROR: chem_init, init_ag_data failed in integer_dict_define'
              flush(6)
              call mesa_error(__FILE__,__LINE__)
           end if
        end do
        call integer_dict_define(Xsol_names_dict, 'p', 1, ierr)
        if (ierr /= 0) then
           write(*,*) 'FATAL ERROR: chem_init, init_ag_data failed in integer_dict_define'
           flush(6)
           call mesa_error(__FILE__,__LINE__)
        end if
        call integer_dict_define(Xsol_names_dict, 'd', 2, ierr)
        if (ierr /= 0) then
           write(*,*) 'FATAL ERROR: chem_init, init_ag_data failed in integer_dict_define'
           flush(6)
           call mesa_error(__FILE__,__LINE__)
        end if
        call integer_dict_create_hash(Xsol_names_dict, ierr)
        if (ierr /= 0) then
           write(*,*) 'FATAL ERROR: chem_init, init_ag_data failed in integer_dict_create_hash'
           flush(6)
           call mesa_error(__FILE__,__LINE__)
        end if
         
      end subroutine init_ag_data
            
      
      subroutine init_chem_element_names
         integer :: i
         chem_element_name(:) = ''

         chem_element_name(1:num_chem_elements) = el_name(1:num_chem_elements)

         do i=1,num_chem_elements
            if (len_trim(chem_element_name(i)) == 0) then
               write(*,*)'missing chem_element_name(i)', i
               flush(6)
               call mesa_error(__FILE__,__LINE__) 
            end if
         end do

      end subroutine init_chem_element_names

      
      
      subroutine init_chem_element_main_iso_names
         ! the iso with the largest number abundance according to Lodders03
         integer :: i
         chem_element_main_iso_name(:) = ''

         !periodic table, row 1
         chem_element_main_iso_name(e_h) = 'h1'
         chem_element_main_iso_name(e_he) = 'he4'

         !periodic table, row 2
         chem_element_main_iso_name(e_li) = 'li7'
         chem_element_main_iso_name(e_be) = 'be9'
         chem_element_main_iso_name(e_b) = 'b11'
         chem_element_main_iso_name(e_c) = 'c12'
         chem_element_main_iso_name(e_n) = 'n14'
         chem_element_main_iso_name(e_o) = 'o16'
         chem_element_main_iso_name(e_f) = 'f19'
         chem_element_main_iso_name(e_ne) = 'ne20'
         
         !periodic table, row 3
         chem_element_main_iso_name(e_na) = 'na23'
         chem_element_main_iso_name(e_mg) = 'mg24'
         chem_element_main_iso_name(e_al) = 'al27'
         chem_element_main_iso_name(e_si) = 'si28'
         chem_element_main_iso_name(e_p) = 'p31'
         chem_element_main_iso_name(e_s) = 's32'
         chem_element_main_iso_name(e_cl) = 'cl35'
         chem_element_main_iso_name(e_ar) = 'ar36'

         !periodic table, row 4
         chem_element_main_iso_name(e_k) = 'k39'
         chem_element_main_iso_name(e_ca) = 'ca40'
         chem_element_main_iso_name(e_sc) = 'sc45'
         chem_element_main_iso_name(e_ti) = 'ti48'
         chem_element_main_iso_name(e_v) = 'v51'
         chem_element_main_iso_name(e_cr) = 'cr52'
         chem_element_main_iso_name(e_mn) = 'mn55'
         chem_element_main_iso_name(e_fe) = 'fe56'
         chem_element_main_iso_name(e_co) = 'co59'
         chem_element_main_iso_name(e_ni) = 'ni58'
         chem_element_main_iso_name(e_cu) = 'cu63'
         chem_element_main_iso_name(e_zn) = 'zn64'
         chem_element_main_iso_name(e_ga) = 'ga69'
         chem_element_main_iso_name(e_ge) = 'ge74'
         chem_element_main_iso_name(e_as) = 'as75'
         chem_element_main_iso_name(e_se) = 'se80'
         chem_element_main_iso_name(e_br) = 'br79'
         chem_element_main_iso_name(e_kr) = 'kr84'

         !periodic table, row 5
         chem_element_main_iso_name(e_rb) = 'rb85'
         chem_element_main_iso_name(e_sr) = 'sr88'
         chem_element_main_iso_name(e_y) = 'y89'
         chem_element_main_iso_name(e_zr) = 'zr90'
         chem_element_main_iso_name(e_nb) = 'nb93'
         chem_element_main_iso_name(e_mo) = 'mo98'
         chem_element_main_iso_name(e_tc) = 'tc97'
         chem_element_main_iso_name(e_ru) = 'ru102'
         chem_element_main_iso_name(e_rh) = 'rh103'
         chem_element_main_iso_name(e_pd) = 'pd106'
         chem_element_main_iso_name(e_ag) = 'ag107'
         chem_element_main_iso_name(e_cd) = 'cd114'
         chem_element_main_iso_name(e_in) = 'in115'
         chem_element_main_iso_name(e_sn) = 'sn120'
         chem_element_main_iso_name(e_sb) = 'sb121'
         chem_element_main_iso_name(e_te) = 'te130'
         chem_element_main_iso_name(e_i ) = 'i127'
         chem_element_main_iso_name(e_xe) = 'xe132'

         !periodic table, row 6
         chem_element_main_iso_name(e_cs) = 'cs133'
         chem_element_main_iso_name(e_ba) = 'ba138'
         chem_element_main_iso_name(e_la) = 'la139'
         chem_element_main_iso_name(e_ce) = 'ce140'
         chem_element_main_iso_name(e_pr) = 'pr141'
         chem_element_main_iso_name(e_nd) = 'nd142'
         chem_element_main_iso_name(e_pm) = 'pm145'
         chem_element_main_iso_name(e_sm) = 'sm152'
         chem_element_main_iso_name(e_eu) = 'eu153'
         chem_element_main_iso_name(e_gd) = 'gd158'
         chem_element_main_iso_name(e_tb) = 'tb159'
         chem_element_main_iso_name(e_dy) = 'dy164'
         chem_element_main_iso_name(e_ho) = 'ho165'
         chem_element_main_iso_name(e_er) = 'er166'
         chem_element_main_iso_name(e_tm) = 'tm169'
         chem_element_main_iso_name(e_yb) = 'yb174'
         chem_element_main_iso_name(e_lu) = 'lu175'
         chem_element_main_iso_name(e_hf) = 'hf180'
         chem_element_main_iso_name(e_ta) = 'ta181'
         chem_element_main_iso_name(e_w ) = 'w184'
         chem_element_main_iso_name(e_re) = 're187'
         chem_element_main_iso_name(e_os) = 'os192'
         chem_element_main_iso_name(e_ir) = 'ir193'
         chem_element_main_iso_name(e_pt) = 'pt195'
         chem_element_main_iso_name(e_au) = 'au197'
         chem_element_main_iso_name(e_hg) = 'hg202'
         chem_element_main_iso_name(e_tl) = 'tl205'
         chem_element_main_iso_name(e_pb) = 'pb208'
         chem_element_main_iso_name(e_bi) = 'bi209'
         chem_element_main_iso_name(e_po) = 'po210'
         chem_element_main_iso_name(e_at) = 'at210'
         ! need a rule here to continue for the unstable nuclei -- longest half-life?
         chem_element_main_iso_name(e_rn) = 'rn222'
         chem_element_main_iso_name(e_fr) = 'fr223'
         chem_element_main_iso_name(e_ra) = 'ra226'
         chem_element_main_iso_name(e_ac) = 'ac227'
         chem_element_main_iso_name(e_th) = 'th232'
         chem_element_main_iso_name(e_pa) = 'pa231'
         chem_element_main_iso_name(e_u)  = 'u238'
         chem_element_main_iso_name(e_np) = 'np237'
         chem_element_main_iso_name(e_pu) = 'pu244'
         chem_element_main_iso_name(e_am) = 'am243'
         chem_element_main_iso_name(e_cm) = 'cm247'
         chem_element_main_iso_name(e_bk) = 'bk247'
         chem_element_main_iso_name(e_cf) = 'cf251'
         chem_element_main_iso_name(e_es) = 'es252'
         chem_element_main_iso_name(e_fm) = 'fm257'
         chem_element_main_iso_name(e_md) = 'md258'
         chem_element_main_iso_name(e_no) = 'no259'
         chem_element_main_iso_name(e_lr) = 'lr262'
         chem_element_main_iso_name(e_rf) = 'rf261'
         chem_element_main_iso_name(e_db) = 'db268'
         chem_element_main_iso_name(e_sg) = 'sg271'
         chem_element_main_iso_name(e_bh) = 'bh274'
         chem_element_main_iso_name(e_hs) = 'hs270'
         chem_element_main_iso_name(e_mt) = 'mt278'
         chem_element_main_iso_name(e_ds) = 'ds281'
         chem_element_main_iso_name(e_rg) = 'rg281'
         chem_element_main_iso_name(e_cn) = 'cn285'


         do i=1,num_chem_elements
            if (len_trim(chem_element_main_iso_name(i)) == 0) then
               write(*,*)'missing chem_element_main_iso_name', i
               flush(6)
               call mesa_error(__FILE__,__LINE__) 
            end if
         end do
         
      end subroutine init_chem_element_main_iso_names
      
      
      subroutine init_element_atomic_weights
         use utils_lib, only: integer_dict_lookup
         integer :: i, isotope_index, ierr
         element_atomic_weight(:) = 0d0

         do i = e_h, e_cn
            call integer_dict_lookup(chem_isos_dict, \
               chem_element_main_iso_name(i), isotope_index, ierr)
            if (ierr /= 0) then
               stop 'init_element_atomic_weights'
               return
            end if
            element_atomic_weight(i) = chem_isos% W(isotope_index)
         end do
      end subroutine init_element_atomic_weights
      
      
      subroutine init_AG89_data ! fraction by mass of total Z
         ! anders & grevesse 1989, paper not available on ADS
         integer :: i
         real(dp) :: z_sum  

         AG89_element_zfrac(:) = 0d0
         
         AG89_element_zfrac(e_li) =     5.2663415161043120D-07
         AG89_element_zfrac(e_be) =     8.7533846996258026D-09
         AG89_element_zfrac(e_b)  =     3.0535981584981396D-07
         AG89_element_zfrac(e_c)  =     1.6164192224602159D-01
         AG89_element_zfrac(e_n)  =     5.8425953868553440D-02
         AG89_element_zfrac(e_o)  =     5.0655547566525427D-01
         AG89_element_zfrac(e_f)  =     2.1339634220190109D-05
         AG89_element_zfrac(e_ne) =     9.2345201069972446D-02
         AG89_element_zfrac(e_na) =     1.7588936076737712D-03
         AG89_element_zfrac(e_mg) =     3.4766459385626718D-02
         AG89_element_zfrac(e_al) =     3.0576538214253392D-03
         AG89_element_zfrac(e_si) =     3.7438035164761560D-02
         AG89_element_zfrac(e_p)  =     4.2953684074804963D-04
         AG89_element_zfrac(e_s)  =     2.2000396518735602D-02
         AG89_element_zfrac(e_cl) =     1.7836963566662123D-04
         AG89_element_zfrac(e_ar) =     4.8868631434862610D-03
         AG89_element_zfrac(e_k)  =     1.8289986382724957D-04
         AG89_element_zfrac(e_ca) =     3.2987013574392004D-03
         AG89_element_zfrac(e_sc) =     2.0504272999081344D-06
         AG89_element_zfrac(e_ti) =     1.5319766333496369D-04
         AG89_element_zfrac(e_v)  =     1.9889381301661414D-05
         AG89_element_zfrac(e_cr) =     9.3528485499287748D-04
         AG89_element_zfrac(e_mn) =     6.9978620325668458D-04
         AG89_element_zfrac(e_fe) =     6.7005139944813982D-02
         AG89_element_zfrac(e_co) =     1.7686377328884702D-04
         AG89_element_zfrac(e_ni) =     3.8650578305194534D-03
         AG89_element_zfrac(e_cu) =     4.4243068859971588D-05
         AG89_element_zfrac(e_zn) =     1.0994428157112287D-04
         
         z_sum = sum(AG89_element_zfrac(:))
         do i = e_li, e_zn
            AG89_element_zfrac(i) = AG89_element_zfrac(i) / z_sum
         end do
         
      end subroutine init_AG89_data


      subroutine init_GN93_data ! fraction by mass of total Z
         ! Grevesse and Noels 1993
         integer :: i
         real(dp) :: z_sum  
         include 'formats.dek'         

         GN93_element_zfrac(:) = -20.0d0
         
        !GN93_element_zfrac(e_H)=12.00d0
        !GN93_element_zfrac(e_He)=10.99d0
         GN93_element_zfrac(e_Li)=3.31d0 !meteor
         GN93_element_zfrac(e_Be)=1.42d0 !meteor
         GN93_element_zfrac(e_B )=2.79d0  !meteor
         GN93_element_zfrac(e_C )=8.55d0
         GN93_element_zfrac(e_N )=7.97d0
         GN93_element_zfrac(e_O )=8.87d0
         GN93_element_zfrac(e_F )=4.48d0
         GN93_element_zfrac(e_Ne)=8.08d0
         GN93_element_zfrac(e_Na)=6.33d0
         GN93_element_zfrac(e_Mg)=7.58d0
         GN93_element_zfrac(e_Al)=6.47d0
         GN93_element_zfrac(e_Si)=7.55d0
         GN93_element_zfrac(e_P )=5.45d0
         GN93_element_zfrac(e_S )=7.20d0
         GN93_element_zfrac(e_Cl)=5.28d0
         GN93_element_zfrac(e_Ar)=6.52d0
         GN93_element_zfrac(e_K )=5.12d0
         GN93_element_zfrac(e_Ca)=6.36d0
         GN93_element_zfrac(e_Sc)=3.17d0
         GN93_element_zfrac(e_Ti)=5.02d0
         GN93_element_zfrac(e_V )=4.00d0
         GN93_element_zfrac(e_Cr)=5.67d0
         GN93_element_zfrac(e_Mn)=5.39d0
         GN93_element_zfrac(e_Fe)=7.50d0
         GN93_element_zfrac(e_Co)=4.92d0
         GN93_element_zfrac(e_Ni)=6.25d0
         GN93_element_zfrac(e_Cu)=4.21d0
         GN93_element_zfrac(e_Zn)=4.60d0
         GN93_element_zfrac(e_Ga)=3.13d0
         GN93_element_zfrac(e_Ge)=3.41d0
         GN93_element_zfrac(e_As)=2.37d0
         GN93_element_zfrac(e_Se)=3.38d0
         GN93_element_zfrac(e_Br)=2.63d0
         GN93_element_zfrac(e_Kr)=3.23d0
         GN93_element_zfrac(e_Rb)=2.41d0
         GN93_element_zfrac(e_Sr)=2.97d0
         GN93_element_zfrac(e_Y )=2.24d0
         GN93_element_zfrac(e_Zr)=2.60d0
         GN93_element_zfrac(e_Nb)=1.42d0
         GN93_element_zfrac(e_Mo)=1.92d0
         GN93_element_zfrac(e_Ru)=1.84d0
         GN93_element_zfrac(e_Rh)=1.12d0
         GN93_element_zfrac(e_Pd)=1.69d0
         GN93_element_zfrac(e_Ag)=1.24d0
         GN93_element_zfrac(e_Cd)=1.77d0
         GN93_element_zfrac(e_In)=0.82d0
         GN93_element_zfrac(e_Sn)=2.14d0
         GN93_element_zfrac(e_Sb)=1.03d0
         GN93_element_zfrac(e_Te)=2.24d0
         GN93_element_zfrac(e_I )=1.51d0
         GN93_element_zfrac(e_Xe)=2.23d0
         GN93_element_zfrac(e_Cs)=1.13d0
         GN93_element_zfrac(e_Ba)=2.13d0
         GN93_element_zfrac(e_La)=1.17d0
         GN93_element_zfrac(e_Ce)=1.58d0
         GN93_element_zfrac(e_Pr)=0.71d0
         GN93_element_zfrac(e_Nd)=1.50d0
         GN93_element_zfrac(e_Sm)=1.01d0
         GN93_element_zfrac(e_Eu)=0.51d0
         GN93_element_zfrac(e_Gd)=1.12d0
         GN93_element_zfrac(e_Tb)=0.35d0
         GN93_element_zfrac(e_Dy)=1.14d0
         GN93_element_zfrac(e_Ho)=0.51d0
         GN93_element_zfrac(e_Er)=0.93d0
         GN93_element_zfrac(e_Tm)=0.15d0
         GN93_element_zfrac(e_Yb)=1.08d0
         GN93_element_zfrac(e_Lu)=0.13d0
         GN93_element_zfrac(e_Hf)=0.88d0
         GN93_element_zfrac(e_Ta)=-0.13d0
         GN93_element_zfrac(e_W )=0.69d0
         GN93_element_zfrac(e_Re)=0.28d0
         GN93_element_zfrac(e_Os)=1.45d0
         GN93_element_zfrac(e_Ir)=1.37d0
         GN93_element_zfrac(e_Pt)=1.69d0
         GN93_element_zfrac(e_Au)=0.87d0
         GN93_element_zfrac(e_Hg)=1.17d0
         GN93_element_zfrac(e_Tl)=0.83d0
         GN93_element_zfrac(e_Pb)=2.06d0
         GN93_element_zfrac(e_Bi)=0.71d0
         GN93_element_zfrac(e_Th)=0.09d0
         GN93_element_zfrac(e_U)=-0.50d0

         ! convert to fraction of Z by mass
         z_sum = 0
         do i = e_li, e_u
            GN93_element_zfrac(i) = &
               exp10(GN93_element_zfrac(i))*element_atomic_weight(i)
            z_sum = z_sum + GN93_element_zfrac(i)
         end do
         do i = e_li, e_u
            GN93_element_zfrac(i) = GN93_element_zfrac(i) / z_sum
         end do

      end subroutine init_GN93_data
      
      
      subroutine init_GS98_data ! fraction by mass of total Z
         ! Grevesse and Sauval 1998, Table 1
         integer :: i
         real(dp) :: z_sum  
         include 'formats.dek'         

         GS98_element_zfrac(:) = -20.0d0

         !GS98_element_zfrac(e_H)=12.00d0
         !GS98_element_zfrac(e_He)=10.93d0
         GS98_element_zfrac(e_Li)=3.31d0 !meteor
         GS98_element_zfrac(e_Be)=1.42d0 !meteor
         GS98_element_zfrac(e_B)=2.79d0  !meteor
         GS98_element_zfrac(e_C)=8.52d0
         GS98_element_zfrac(e_N)=7.92d0
         GS98_element_zfrac(e_O)=8.83d0
         GS98_element_zfrac(e_F)=4.48d0
         GS98_element_zfrac(e_Ne)=8.08d0
         GS98_element_zfrac(e_Na)=6.32d0
         GS98_element_zfrac(e_Mg)=7.58d0
         GS98_element_zfrac(e_Al)=6.49d0
         GS98_element_zfrac(e_Si)=7.56d0
         GS98_element_zfrac(e_P)=5.56d0
         GS98_element_zfrac(e_S)=7.20d0
         GS98_element_zfrac(e_Cl)=5.28d0
         GS98_element_zfrac(e_Ar)=6.40d0
         GS98_element_zfrac(e_K)=5.13d0
         GS98_element_zfrac(e_Ca)=6.35d0
         GS98_element_zfrac(e_Sc)=3.10d0
         GS98_element_zfrac(e_Ti)=4.94d0
         GS98_element_zfrac(e_V)=4.02d0
         GS98_element_zfrac(e_Cr)=5.69d0
         GS98_element_zfrac(e_Mn)=5.53d0
         GS98_element_zfrac(e_Fe)=7.50d0
         GS98_element_zfrac(e_Co)=4.91d0
         GS98_element_zfrac(e_Ni)=6.25d0
         GS98_element_zfrac(e_Cu)=4.29d0
         GS98_element_zfrac(e_Zn)=4.67d0
         GS98_element_zfrac(e_Ga)=3.13d0
         GS98_element_zfrac(e_Ge)=3.63d0
         GS98_element_zfrac(e_As)=2.37d0
         GS98_element_zfrac(e_Se)=3.41d0
         GS98_element_zfrac(e_Br)=2.63d0
         GS98_element_zfrac(e_Kr)=3.31d0
         GS98_element_zfrac(e_Rb)=2.41d0
         GS98_element_zfrac(e_Sr)=2.92d0
         GS98_element_zfrac(e_Y)=2.23d0
         GS98_element_zfrac(e_Zr)=2.61d0
         GS98_element_zfrac(e_Nb)=1.40d0
         GS98_element_zfrac(e_Mo)=1.97d0
         GS98_element_zfrac(e_Ru)=1.83d0
         GS98_element_zfrac(e_Rh)=1.10d0
         GS98_element_zfrac(e_Pd)=1.70d0
         GS98_element_zfrac(e_Ag)=1.24d0
         GS98_element_zfrac(e_Cd)=1.76d0
         GS98_element_zfrac(e_In)=0.82d0
         GS98_element_zfrac(e_Sn)=2.14d0
         GS98_element_zfrac(e_Sb)=1.03d0
         GS98_element_zfrac(e_Te)=2.24d0
         GS98_element_zfrac(e_I)=1.51d0
         GS98_element_zfrac(e_Xe)=2.17d0
         GS98_element_zfrac(e_Cs)=1.13d0
         GS98_element_zfrac(e_Ba)=2.22d0
         GS98_element_zfrac(e_La)=1.22d0
         GS98_element_zfrac(e_Ce)=1.63d0
         GS98_element_zfrac(e_Pr)=0.80d0
         GS98_element_zfrac(e_Nd)=1.49d0
         GS98_element_zfrac(e_Sm)=0.98d0
         GS98_element_zfrac(e_Eu)=0.55d0
         GS98_element_zfrac(e_Gd)=1.09d0
         GS98_element_zfrac(e_Tb)=0.35d0
         GS98_element_zfrac(e_Dy)=1.17d0
         GS98_element_zfrac(e_Ho)=0.51d0
         GS98_element_zfrac(e_Er)=0.97d0
         GS98_element_zfrac(e_Tm)=0.15d0
         GS98_element_zfrac(e_Yb)=0.96d0
         GS98_element_zfrac(e_Lu)=0.13d0
         GS98_element_zfrac(e_Hf)=0.75d0
         GS98_element_zfrac(e_Ta)=-0.13d0
         GS98_element_zfrac(e_W)=0.69d0
         GS98_element_zfrac(e_Re)=0.28d0
         GS98_element_zfrac(e_Os)=1.39d0
         GS98_element_zfrac(e_Ir)=1.37d0
         GS98_element_zfrac(e_Pt)=1.69d0
         GS98_element_zfrac(e_Au)=0.85d0
         GS98_element_zfrac(e_Hg)=1.13d0
         GS98_element_zfrac(e_Tl)=0.83d0
         GS98_element_zfrac(e_Pb)=2.06d0
         GS98_element_zfrac(e_Bi)=0.71d0
         GS98_element_zfrac(e_Th)=0.09d0
         GS98_element_zfrac(e_U)=-0.50d0
         
         ! convert to fraction of Z by mass
         z_sum = 0d0
         do i = e_li, e_u
            GS98_element_zfrac(i) = &
               exp10(GS98_element_zfrac(i))*element_atomic_weight(i)
            z_sum = z_sum + GS98_element_zfrac(i)
         end do
         do i = e_li, e_u
            GS98_element_zfrac(i) = GS98_element_zfrac(i) / z_sum
         end do
         
      end subroutine init_GS98_data
            
      subroutine init_L03_data ! fraction by mass of total Z
         ! Lodders 2003, ApJ, Table 1 recommended abundance
         integer :: i
         real(dp) :: z_sum  
         include 'formats.dek'         

         L03_element_zfrac(:) = -20.0d0
         
        !L03_element_zfrac(e_H)=12d0
        !L03_element_zfrac(e_He)=10.89d0
        L03_element_zfrac(e_Li)=3.28d0
        L03_element_zfrac(e_Be)=1.41d0
        L03_element_zfrac(e_B)=2.78d0
        L03_element_zfrac(e_C)=8.39d0
        L03_element_zfrac(e_N)=7.83d0
        L03_element_zfrac(e_O)=8.69d0
        L03_element_zfrac(e_F)=4.46d0
        L03_element_zfrac(e_Ne)=7.87d0
        L03_element_zfrac(e_Na)=6.30d0
        L03_element_zfrac(e_Mg)=7.55d0
        L03_element_zfrac(e_Al)=6.46d0
        L03_element_zfrac(e_Si)=7.54d0
        L03_element_zfrac(e_P)=5.46d0
        L03_element_zfrac(e_S)=7.19d0
        L03_element_zfrac(e_Cl)=5.26d0
        L03_element_zfrac(e_Ar)=6.55d0
        L03_element_zfrac(e_K)=5.11d0
        L03_element_zfrac(e_Ca)=6.34d0
        L03_element_zfrac(e_Sc)=3.07d0
        L03_element_zfrac(e_Ti)=4.92d0
        L03_element_zfrac(e_V)=4.00d0
        L03_element_zfrac(e_Cr)=5.65d0
        L03_element_zfrac(e_Mn)=5.50d0
        L03_element_zfrac(e_Fe)=7.47d0
        L03_element_zfrac(e_Co)=4.91d0
        L03_element_zfrac(e_Ni)=6.22d0
        L03_element_zfrac(e_Cu)=4.26d0
        L03_element_zfrac(e_Zn)=4.63d0
        L03_element_zfrac(e_Ga)=3.10d0
        L03_element_zfrac(e_Ge)=3.62d0
        L03_element_zfrac(e_As)=2.32d0
        L03_element_zfrac(e_Se)=3.36d0
        L03_element_zfrac(e_Br)=2.59d0
        L03_element_zfrac(e_Kr)=3.28d0
        L03_element_zfrac(e_Rb)=2.36d0
        L03_element_zfrac(e_Sr)=2.91d0
        L03_element_zfrac(e_Y)=2.20d0
        L03_element_zfrac(e_Zr)=2.60d0
        L03_element_zfrac(e_Nb)=1.42d0
        L03_element_zfrac(e_Mo)=1.96d0
        L03_element_zfrac(e_Ru)=1.82d0
        L03_element_zfrac(e_Rh)=1.11d0
        L03_element_zfrac(e_Pd)=1.70d0
        L03_element_zfrac(e_Ag)=1.23d0
        L03_element_zfrac(e_Cd)=1.74d0
        L03_element_zfrac(e_In)=0.80d0
        L03_element_zfrac(e_Sn)=2.11d0
        L03_element_zfrac(e_Sb)=1.06d0
        L03_element_zfrac(e_Te)=2.22d0
        L03_element_zfrac(e_I)=1.54d0
        L03_element_zfrac(e_Xe)=2.27d0
        L03_element_zfrac(e_Cs)=1.10d0
        L03_element_zfrac(e_Ba)=2.18d0
        L03_element_zfrac(e_La)=1.18d0
        L03_element_zfrac(e_Ce)=1.61d0
        L03_element_zfrac(e_Pr)=0.78d0
        L03_element_zfrac(e_Nd)=1.46d0
        L03_element_zfrac(e_Sm)=0.95d0
        L03_element_zfrac(e_Eu)=0.52d0
        L03_element_zfrac(e_Gd)=1.06d0
        L03_element_zfrac(e_Tb)=0.31d0
        L03_element_zfrac(e_Dy)=1.13d0
        L03_element_zfrac(e_Ho)=0.49d0
        L03_element_zfrac(e_Er)=0.95d0
        L03_element_zfrac(e_Tm)=0.11d0
        L03_element_zfrac(e_Yb)=0.94d0
        L03_element_zfrac(e_Lu)=0.09d0
        L03_element_zfrac(e_Hf)=0.77d0
        L03_element_zfrac(e_Ta)=-0.14d0
        L03_element_zfrac(e_W)=0.65d0
        L03_element_zfrac(e_Re)=0.26d0
        L03_element_zfrac(e_Os)=1.37d0
        L03_element_zfrac(e_Ir)=1.35d0
        L03_element_zfrac(e_Pt)=1.67d0
        L03_element_zfrac(e_Au)=0.83d0
        L03_element_zfrac(e_Hg)=1.16d0
        L03_element_zfrac(e_Tl)=0.81d0
        L03_element_zfrac(e_Pb)=2.05d0
        L03_element_zfrac(e_Bi)=0.68d0
        L03_element_zfrac(e_Th)=0.09d0
        L03_element_zfrac(e_U)=-0.49d0

         ! convert to fraction of Z by mass
         z_sum = 0d0
         do i = e_li, e_u
            L03_element_zfrac(i) = &
               exp10(L03_element_zfrac(i))*element_atomic_weight(i)
            z_sum = z_sum + L03_element_zfrac(i)
         end do
         do i = e_li, e_u
            L03_element_zfrac(i) = L03_element_zfrac(i) / z_sum
         end do
         
      end subroutine init_L03_data
      
      
      subroutine init_AGS05_data ! fraction by mass of total Z
         ! Asplund, Grevesse and Sauval 2005      
         integer :: i
         real(dp) :: z_sum  
         include 'formats.dek'         

         AGS05_element_zfrac(:) = -20.0d0
                  
         ! first store log abundances from the paper (photosphere unless otherwise noted)
         ! relative to log abundance of H = 12.00d0
         AGS05_element_zfrac(e_li) = 3.25d0 !meteor
         AGS05_element_zfrac(e_be) = 1.38d0
         AGS05_element_zfrac(e_b ) = 2.70d0
         AGS05_element_zfrac(e_c ) = 8.39d0
         AGS05_element_zfrac(e_n ) = 7.78d0
         AGS05_element_zfrac(e_o ) = 8.66d0
         AGS05_element_zfrac(e_f ) = 4.56d0
         AGS05_element_zfrac(e_ne) = 7.84d0 !indirect
         AGS05_element_zfrac(e_na) = 6.17d0
         AGS05_element_zfrac(e_mg) = 7.53d0
         AGS05_element_zfrac(e_al) = 6.37d0
         AGS05_element_zfrac(e_si) = 7.51d0
         AGS05_element_zfrac(e_p ) = 5.36d0
         AGS05_element_zfrac(e_s ) = 7.14d0
         AGS05_element_zfrac(e_cl) = 5.50d0
         AGS05_element_zfrac(e_ar) = 6.18d0 !indirect
         AGS05_element_zfrac(e_k ) = 5.08d0
         AGS05_element_zfrac(e_ca) = 6.31d0
         AGS05_element_zfrac(e_sc) = 3.05d0
         AGS05_element_zfrac(e_ti) = 4.90d0
         AGS05_element_zfrac(e_v ) = 4.00d0
         AGS05_element_zfrac(e_cr) = 5.64d0
         AGS05_element_zfrac(e_mn) = 5.39d0
         AGS05_element_zfrac(e_fe) = 7.45d0
         AGS05_element_zfrac(e_co) = 4.92d0
         AGS05_element_zfrac(e_ni) = 6.23d0
         AGS05_element_zfrac(e_cu) = 4.21d0
         AGS05_element_zfrac(e_zn) = 4.60d0
         AGS05_element_zfrac(e_ga) = 2.88d0
         AGS05_element_zfrac(e_ge) = 3.58d0
         AGS05_element_zfrac(e_as) = 2.29d0 !meteor
         AGS05_element_zfrac(e_se) = 3.33d0 !meteor
         AGS05_element_zfrac(e_br) = 2.56d0 !meteor
         AGS05_element_zfrac(e_kr) = 3.28d0 !indirect
         AGS05_element_zfrac(e_rb) = 2.60d0
         AGS05_element_zfrac(e_sr) = 2.92d0
         AGS05_element_zfrac(e_y ) = 2.21d0
         AGS05_element_zfrac(e_zr) = 2.59d0
         AGS05_element_zfrac(e_nb) = 1.42d0
         AGS05_element_zfrac(e_mo) = 1.92d0
         AGS05_element_zfrac(e_Ru) = 1.84d0
         AGS05_element_zfrac(e_Rh) = 1.12d0
         AGS05_element_zfrac(e_Pd) = 1.69d0
         AGS05_element_zfrac(e_Ag) = 0.94d0
         AGS05_element_zfrac(e_Cd) = 1.77d0
         AGS05_element_zfrac(e_In) = 1.60d0
         AGS05_element_zfrac(e_Sn) = 2.00d0
         AGS05_element_zfrac(e_Sb) = 1.00d0
         AGS05_element_zfrac(e_Te) = 2.19d0 !meteor
         AGS05_element_zfrac(e_I ) = 1.51d0 !meteor
         AGS05_element_zfrac(e_Xe) = 2.27d0 !indirect
         AGS05_element_zfrac(e_Cs) = 1.07d0 !meteor
         AGS05_element_zfrac(e_Ba) = 2.17d0
         AGS05_element_zfrac(e_La) = 1.13d0
         AGS05_element_zfrac(e_Ce) = 1.58d0
         AGS05_element_zfrac(e_Pr) = 0.71d0
         AGS05_element_zfrac(e_Nd) = 1.45d0
         AGS05_element_zfrac(e_Sm) = 1.01d0
         AGS05_element_zfrac(e_Eu) = 0.52d0
         AGS05_element_zfrac(e_Gd) = 1.12d0
         AGS05_element_zfrac(e_Tb) = 0.28d0
         AGS05_element_zfrac(e_Dy) = 1.14d0
         AGS05_element_zfrac(e_Ho) = 0.51d0
         AGS05_element_zfrac(e_Er) = 0.93d0
         AGS05_element_zfrac(e_Tm) = 0.00d0
         AGS05_element_zfrac(e_Yb) = 1.08d0
         AGS05_element_zfrac(e_Lu) = 0.06d0
         AGS05_element_zfrac(e_Hf) = 0.88d0
         AGS05_element_zfrac(e_Ta) = -0.17d0 !meteor
         AGS05_element_zfrac(e_W ) = 1.11d0
         AGS05_element_zfrac(e_Re) = 0.23d0 !meteor
         AGS05_element_zfrac(e_Os) = 1.45d0
         AGS05_element_zfrac(e_Ir) = 1.38d0
         AGS05_element_zfrac(e_Pt) = 1.64d0 !meteor
         AGS05_element_zfrac(e_Au) = 1.01d0
         AGS05_element_zfrac(e_Hg) = 1.13d0 !meteor
         AGS05_element_zfrac(e_Tl) = 0.90d0
         AGS05_element_zfrac(e_Pb) = 2.00d0
         AGS05_element_zfrac(e_Bi) = 0.65d0 !meteor
         AGS05_element_zfrac(e_Th) = 0.06d0 !meteor
         AGS05_element_zfrac(e_U) = -0.52d0

         ! convert to fraction of Z by mass
         z_sum = 0d0
         do i = e_li, e_u
            AGS05_element_zfrac(i) = &
               exp10(AGS05_element_zfrac(i))*element_atomic_weight(i)
            z_sum = z_sum + AGS05_element_zfrac(i)
         end do
         do i = e_li, e_u
            AGS05_element_zfrac(i) = AGS05_element_zfrac(i) / z_sum
         end do
         
      end subroutine init_AGS05_data
      
      
      subroutine init_AGSS09_data ! fraction by mass of total Z
         ! Asplund, Grevesse, Sauval, and Scott 2009 abundances
         ! Annu. Rev. Astron. Astrophys. 2009. 47:481–522
         integer :: i
         real(dp) :: z_sum  
         include 'formats.dek'         
         
         AGSS09_element_zfrac(:) = -20.0d0
         
         ! first store log abundances from the paper
         AGSS09_element_zfrac(e_li) = 3.26d0
         AGSS09_element_zfrac(e_be) = 1.38d0
         AGSS09_element_zfrac(e_b ) = 2.70d0
         AGSS09_element_zfrac(e_c ) = 8.43d0
         AGSS09_element_zfrac(e_n ) = 7.83d0
         AGSS09_element_zfrac(e_o ) = 8.69d0
         AGSS09_element_zfrac(e_f ) = 4.56d0
         AGSS09_element_zfrac(e_ne) = 7.93d0
         AGSS09_element_zfrac(e_na) = 6.24d0
         AGSS09_element_zfrac(e_mg) = 7.60d0
         AGSS09_element_zfrac(e_al) = 6.45d0
         AGSS09_element_zfrac(e_si) = 7.51d0
         AGSS09_element_zfrac(e_p ) = 5.41d0
         AGSS09_element_zfrac(e_s ) = 7.12d0
         AGSS09_element_zfrac(e_cl) = 5.50d0
         AGSS09_element_zfrac(e_ar) = 6.40d0
         AGSS09_element_zfrac(e_k ) = 5.03d0
         AGSS09_element_zfrac(e_ca) = 6.34d0
         AGSS09_element_zfrac(e_sc) = 3.15d0
         AGSS09_element_zfrac(e_ti) = 4.95d0
         AGSS09_element_zfrac(e_v ) = 3.93d0
         AGSS09_element_zfrac(e_cr) = 5.64d0
         AGSS09_element_zfrac(e_mn) = 5.43d0
         AGSS09_element_zfrac(e_fe) = 7.50d0
         AGSS09_element_zfrac(e_co) = 4.99d0
         AGSS09_element_zfrac(e_ni) = 6.22d0
         AGSS09_element_zfrac(e_cu) = 4.19d0
         AGSS09_element_zfrac(e_zn) = 4.56d0
         AGSS09_element_zfrac(e_ga) = 3.04d0
         AGSS09_element_zfrac(e_ge) = 3.65d0
         AGSS09_element_zfrac(e_as) = 2.30d0 !meteor
         AGSS09_element_zfrac(e_se) = 3.34d0 !meteor
         AGSS09_element_zfrac(e_br) = 2.54d0 !meteor
         AGSS09_element_zfrac(e_kr) = 3.25d0 !indirect
         AGSS09_element_zfrac(e_rb) = 2.52d0
         AGSS09_element_zfrac(e_sr) = 2.87d0
         AGSS09_element_zfrac(e_y ) = 2.21d0
         AGSS09_element_zfrac(e_zr) = 2.58d0
         AGSS09_element_zfrac(e_nb) = 1.46d0
         AGSS09_element_zfrac(e_mo) = 1.88d0
         AGSS09_element_zfrac(e_Ru) = 1.75d0
         AGSS09_element_zfrac(e_Rh) = 0.91d0
         AGSS09_element_zfrac(e_Pd) = 1.57d0
         AGSS09_element_zfrac(e_Ag) = 0.94d0
         AGSS09_element_zfrac(e_Cd) = 1.71d0
         AGSS09_element_zfrac(e_In) = 0.80d0
         AGSS09_element_zfrac(e_Sn) = 2.04d0
         AGSS09_element_zfrac(e_Sb) = 1.01d0
         AGSS09_element_zfrac(e_Te) = 2.18d0
         AGSS09_element_zfrac(e_I ) = 1.55d0
         AGSS09_element_zfrac(e_Xe) = 2.24d0
         AGSS09_element_zfrac(e_Cs) = 1.08d0
         AGSS09_element_zfrac(e_Ba) = 2.18d0
         AGSS09_element_zfrac(e_La) = 1.10d0
         AGSS09_element_zfrac(e_Ce) = 1.58d0
         AGSS09_element_zfrac(e_Pr) = 0.72d0
         AGSS09_element_zfrac(e_Nd) = 1.42d0
         AGSS09_element_zfrac(e_Sm) = 0.96d0
         AGSS09_element_zfrac(e_Eu) = 0.52d0
         AGSS09_element_zfrac(e_Gd) = 1.07d0
         AGSS09_element_zfrac(e_Tb) = 0.30d0
         AGSS09_element_zfrac(e_Dy) = 1.10d0
         AGSS09_element_zfrac(e_Ho) = 0.48d0
         AGSS09_element_zfrac(e_Er) = 0.92d0
         AGSS09_element_zfrac(e_Tm) = 0.10d0
         AGSS09_element_zfrac(e_Yb) = 0.84d0
         AGSS09_element_zfrac(e_Lu) = 0.10d0
         AGSS09_element_zfrac(e_Hf) = 0.85d0
         AGSS09_element_zfrac(e_Ta) = -0.12d0
         AGSS09_element_zfrac(e_W ) = 0.85d0
         AGSS09_element_zfrac(e_Re) = 0.26d0
         AGSS09_element_zfrac(e_Os) = 1.40d0
         AGSS09_element_zfrac(e_Ir) = 1.38d0
         AGSS09_element_zfrac(e_Pt) = 1.62d0
         AGSS09_element_zfrac(e_Au) = 0.92d0
         AGSS09_element_zfrac(e_Hg) = 1.17d0
         AGSS09_element_zfrac(e_Tl) = 0.90d0
         AGSS09_element_zfrac(e_Pb) = 1.75d0
         AGSS09_element_zfrac(e_Bi) = 0.65d0
         AGSS09_element_zfrac(e_Th) = 0.02d0
         AGSS09_element_zfrac(e_U) = -0.54d0
         
         ! convert to fraction of Z by mass
         z_sum = 0
         do i = e_li, e_u
            AGSS09_element_zfrac(i) = &
               exp10(AGSS09_element_zfrac(i))*element_atomic_weight(i)
            z_sum = z_sum + AGSS09_element_zfrac(i)
         end do
         do i = e_li, e_u
            AGSS09_element_zfrac(i) = AGSS09_element_zfrac(i) / z_sum
         end do
         
      end subroutine init_AGSS09_data


      subroutine init_A09_Przybilla_data ! fraction by mass of total Z
         ! provided by Ehsan Moravveji, Oct 12, 2013.
         ! The mass fraction is taken from Asplund et al. (2009), and modified by the
         ! B-star measurement of Nieva & Przybilla 2012, A&A, 539, 143 and
         ! Przybilla et al. (2013), EAS proceeding to be published
         ! The modified elements are: he, c, n, o, ne, mg, al, si, s, ar, fe
         integer :: i
         real(dp) :: z_sum  
         include 'formats.dek'
         
         A09_Prz_zfrac(:)  = -20.0d0

!         A09_Prz_zfrac(e_h ) = 12.00d0
!         A09_Prz_zfrac(e_he) = 10.99d0
         A09_Prz_zfrac(e_li) = 3.26d0
         A09_Prz_zfrac(e_be) = 1.38d0
         A09_Prz_zfrac(e_b ) = 2.70d0
         A09_Prz_zfrac(e_c ) = 8.33d0
         A09_Prz_zfrac(e_n ) = 7.79d0
         A09_Prz_zfrac(e_o ) = 8.76d0
         A09_Prz_zfrac(e_f ) = 4.56d0
         A09_Prz_zfrac(e_ne) = 8.09d0
         A09_Prz_zfrac(e_na) = 6.24d0
         A09_Prz_zfrac(e_mg) = 7.56d0
         A09_Prz_zfrac(e_al) = 6.30d0
         A09_Prz_zfrac(e_si) = 7.50d0
         A09_Prz_zfrac(e_p ) = 5.41d0
         A09_Prz_zfrac(e_s ) = 7.14d0
         A09_Prz_zfrac(e_cl) = 5.50d0
         A09_Prz_zfrac(e_ar) = 6.50d0
         A09_Prz_zfrac(e_k ) = 5.03d0
         A09_Prz_zfrac(e_ca) = 6.34d0
         A09_Prz_zfrac(e_sc) = 3.15d0
         A09_Prz_zfrac(e_ti) = 4.95d0
         A09_Prz_zfrac(e_v ) = 3.93d0
         A09_Prz_zfrac(e_cr) = 5.64d0
         A09_Prz_zfrac(e_mn) = 5.43d0
         A09_Prz_zfrac(e_fe) = 7.52d0
         A09_Prz_zfrac(e_co) = 4.99d0
         A09_Prz_zfrac(e_ni) = 6.22d0
         A09_Prz_zfrac(e_cu) = 4.19d0
         A09_Prz_zfrac(e_zn) = 4.56d0
         A09_Prz_zfrac(e_ga) = 3.04d0
         A09_Prz_zfrac(e_ge) = 3.65d0
         A09_Prz_zfrac(e_as) = 2.30d0
         A09_Prz_zfrac(e_se) = 3.34d0
         A09_Prz_zfrac(e_br) = 2.54d0
         A09_Prz_zfrac(e_kr) = 3.25d0
         A09_Prz_zfrac(e_rb) = 2.52d0
         A09_Prz_zfrac(e_sr) = 2.87d0
         A09_Prz_zfrac(e_y ) = 2.21d0
         A09_Prz_zfrac(e_zr) = 2.58d0
         A09_Prz_zfrac(e_nb) = 1.46d0
         A09_Prz_zfrac(e_mo) = 1.88d0
         A09_Prz_zfrac(e_ru) = 1.75d0
         A09_Prz_zfrac(e_rh) = 0.91d0
         A09_Prz_zfrac(e_pd) = 1.57d0
         A09_Prz_zfrac(e_ag) = 0.94d0
         A09_Prz_zfrac(e_cd) = 1.71d0
         A09_Prz_zfrac(e_in) = 0.80d0
         A09_Prz_zfrac(e_sn) = 2.04d0
         A09_Prz_zfrac(e_sb) = 1.01d0
         A09_Prz_zfrac(e_te) = 2.18d0
         A09_Prz_zfrac(e_i ) = 1.55d0
         A09_Prz_zfrac(e_xe) = 2.24d0
         A09_Prz_zfrac(e_cs) = 1.08d0
         A09_Prz_zfrac(e_ba) = 2.18d0
         A09_Prz_zfrac(e_la) = 1.10d0
         A09_Prz_zfrac(e_ce) = 1.58d0
         A09_Prz_zfrac(e_pr) = 0.72d0
         A09_Prz_zfrac(e_nd) = 1.42d0
         A09_Prz_zfrac(e_sm) = 0.96d0
         A09_Prz_zfrac(e_eu) = 0.52d0
         A09_Prz_zfrac(e_gd) = 1.07d0
         A09_Prz_zfrac(e_tb) = 0.30d0
         A09_Prz_zfrac(e_dy) = 1.10d0
         A09_Prz_zfrac(e_ho) = 0.48d0
         A09_Prz_zfrac(e_er) = 0.92d0
         A09_Prz_zfrac(e_tm) = 0.10d0
         A09_Prz_zfrac(e_yb) = 0.84d0
         A09_Prz_zfrac(e_lu) = 0.10d0
         A09_Prz_zfrac(e_hf) = 0.85d0
         A09_Prz_zfrac(e_ta) = -0.12d0
         A09_Prz_zfrac(e_w ) = 0.85d0
         A09_Prz_zfrac(e_re) = 0.26d0
         A09_Prz_zfrac(e_os) = 1.40d0
         A09_Prz_zfrac(e_ir) = 1.38d0
         A09_Prz_zfrac(e_pt) = 1.62d0
         A09_Prz_zfrac(e_au) = 0.92d0
         A09_Prz_zfrac(e_hg) = 1.17d0
         A09_Prz_zfrac(e_tl) = 0.90d0
         A09_Prz_zfrac(e_pb) = 1.75d0
         A09_Prz_zfrac(e_bi) = 0.65d0
         A09_Prz_zfrac(e_th) = 0.02d0
         A09_Prz_zfrac(e_u ) = -0.54d0
                           
         ! convert to fraction of Z by mass
         z_sum = 0d0
         do i = e_li, e_u
            A09_Prz_zfrac(i) = &
               exp10(A09_Prz_zfrac(i))*element_atomic_weight(i)
            z_sum = z_sum + A09_Prz_zfrac(i)
         end do
         do i = e_li, e_u
            A09_Prz_zfrac(i) = A09_Prz_zfrac(i) / z_sum
         end do
         
      end subroutine init_A09_Przybilla_data
      
      
      subroutine init_L09_data ! fraction by mass of total Z
         ! Lodders 09
         integer :: i
         real(dp) :: z_sum  
         include 'formats.dek'
         
         L09_element_zfrac(:) = 0
         
         ! mass fractions
         L09_element_zfrac(e_li) = 1.054594933683d-08
         L09_element_zfrac(e_be) = 1.5087555571d-10
         L09_element_zfrac(e_b ) = 5.5633306761d-09
         L09_element_zfrac(e_c ) = 0.002365544090848D0
         L09_element_zfrac(e_n ) = 0.0008161934768925D0
         L09_element_zfrac(e_o ) = 0.0068991682737478D0
         L09_element_zfrac(e_f ) = 4.1844135602D-07
         L09_element_zfrac(e_ne) = 0.0018162023028794D0
         L09_element_zfrac(e_na) = 3.6352024324D-05
         L09_element_zfrac(e_mg) = 0.000683514477408D0
         L09_element_zfrac(e_al) = 6.2568980455D-05
         L09_element_zfrac(e_si) = 0.000769722819625D0
         L09_element_zfrac(e_p ) = 7.0479812061D-06
         L09_element_zfrac(e_s ) = 0.000370123705833609D0
         L09_element_zfrac(e_cl) = 5.0250763788D-06
         L09_element_zfrac(e_ar) = 9.2220355099056D-05
         L09_element_zfrac(e_k ) = 4.0297305059079D-06
         L09_element_zfrac(e_ca) = 6.63135868398494D-05
         L09_element_zfrac(e_sc) = 4.2402933957D-08
         L09_element_zfrac(e_ti) = 3.24207135484D-06
         L09_element_zfrac(e_v ) = 4.0049954611698D-07
         L09_element_zfrac(e_cr) = 1.870484358513D-05
         L09_element_zfrac(e_mn) = 1.3890521841D-05
         L09_element_zfrac(e_fe) = 0.0012986862725768D0
         L09_element_zfrac(e_co) = 3.7979113651D-06
         L09_element_zfrac(e_ni) = 7.901833309378D-05
         L09_element_zfrac(e_cu) = 9.4275308656D-07
         L09_element_zfrac(e_zn) = 2.324135489336D-06
         L09_element_zfrac(e_ga) = 6.9975797859D-08
         L09_element_zfrac(e_ge) = 2.27918509226D-07
         L09_element_zfrac(e_as) = 1.2531874861D-08
         L09_element_zfrac(e_se) = 1.460613983153D-07
         L09_element_zfrac(e_br) = 2.1909826069D-07
         L09_element_zfrac(e_kr) = 1.2837104766341D-07
         L09_element_zfrac(e_rb) = 1.69493949899D-08
         L09_element_zfrac(e_sr) = 5.58119030984D-08
         L09_element_zfrac(e_y ) = 1.1287452840D-08
         L09_element_zfrac(e_zr) = 2.696974516888D-08
         L09_element_zfrac(e_nb) = 1.9870212075D-09
         L09_element_zfrac(e_mo) = 6.70675811282D-09
         L09_element_zfrac(e_Ru) = 4.935148192679D-09
         L09_element_zfrac(e_Rh) = 1.0439120240D-09
         L09_element_zfrac(e_Pd) = 3.958834334296D-09
         L09_element_zfrac(e_Ag) = 1.44909561511D-09
         L09_element_zfrac(e_Cd) = 4.850725813727D-09
         L09_element_zfrac(e_In) = 5.60277526584D-10
         L09_element_zfrac(e_Sn) = 1.1748571050167D-08
         L09_element_zfrac(e_Sb) = 1.04476117832D-09
         L09_element_zfrac(e_Te) = 1.64381492798D-08
         L09_element_zfrac(e_I ) = 3.8266730451D-09
         L09_element_zfrac(e_Xe) = 1.9631791443838D-08
         L09_element_zfrac(e_Cs) = 1.3516072159D-09
         L09_element_zfrac(e_Ba) = 1.6848783894182D-08
         L09_element_zfrac(e_La) = 1.7400268564D-09
         L09_element_zfrac(e_Ce) = 4.5166246603333D-09
         L09_element_zfrac(e_Pr) = 6.6431263199D-10
         L09_element_zfrac(e_Nd) = 5.82366496835D-09
         L09_element_zfrac(e_Sm) = 1.10017534847D-09
         L09_element_zfrac(e_Eu) = 4.1023195079D-10
         L09_element_zfrac(e_Gd) = 1.5505832574086D-09
         L09_element_zfrac(e_Tb) = 2.7612856334D-10
         L09_element_zfrac(e_Dy) = 1.79997715438722D-09
         L09_element_zfrac(e_Ho) = 4.1129202414D-10
         L09_element_zfrac(e_Er) = 1.2035968713737D-09
         L09_element_zfrac(e_Tm) = 1.8794799164D-10
         L09_element_zfrac(e_Yb) = 1.2153809425605D-09
         L09_element_zfrac(e_Lu) = 1.826667993422D-10
         L09_element_zfrac(e_Hf) = 7.619544268076D-10
         L09_element_zfrac(e_Ta) = 1.04130101121661D-10
         L09_element_zfrac(e_W ) = 6.9059258989518D-10
         L09_element_zfrac(e_Re) = 2.9647265833D-10
         L09_element_zfrac(e_Os) = 3.52828020217507D-09
         L09_element_zfrac(e_Ir) = 3.5336600059D-09
         L09_element_zfrac(e_Pt) = 6.8100262401586D-09
         L09_element_zfrac(e_Au) = 1.0522666072D-09
         L09_element_zfrac(e_Hg) = 2.5169209215851D-09
         L09_element_zfrac(e_Tl) = 1.02465539439D-09
         L09_element_zfrac(e_Pb) = 1.879940103275D-08
         L09_element_zfrac(e_Bi) = 7.9004226175D-10
         L09_element_zfrac(e_Th) = 2.7961831384D-10
         L09_element_zfrac(e_U) = 1.546830543D-10
         
         ! convert from mass fraction to fraction of Z by mass
         z_sum = sum(L09_element_zfrac(e_li:e_u))
         do i = e_li, e_u
            L09_element_zfrac(i) = L09_element_zfrac(i) / z_sum
         end do
                  
      end subroutine init_L09_data


      subroutine allocate_nuclide_data(d,n,ierr)
         type(nuclide_data), intent(out) :: d
         integer, intent(in) :: n
         integer, intent(out) :: ierr  
         ierr = 0
         allocate(d% name(n), d% W(n), d% Z(n), d% N(n), d% Z_plus_N(n),  &
               d% spin(N), d% binding_energy(n), d% Z53(n), &
               d% isomeric_state(n), d% mass_excess(n), d% pfcn(npart,n),  &
               d% chem_id(n), d% nuclide(n), stat=ierr)
         if (ierr /= 0) return
         d% nnuclides = n
      end subroutine allocate_nuclide_data
      

      subroutine free_nuclide_data(n)
         type(nuclide_data), intent(inout) :: n
         if (associated(n% name)) &
            deallocate( &
               n% name, n% W, n% Z, n% N, n% Z_plus_N, n% spin, n% binding_energy, n% Z53, &
               n% isomeric_state, n% mass_excess, n% pfcn, n% chem_id, n% nuclide)
         n% nnuclides = 0
      end subroutine free_nuclide_data
      
      
      subroutine free_lodders03_table()
            use utils_lib, only : integer_dict_free
            deallocate(lodders03_tab6% isotopic_percent)
            call integer_dict_free(lodders03_tab6% name_dict)
            nullify(lodders03_tab6% name_dict)
      end subroutine free_lodders03_table


      ! returns the index of a particular nuclide in the full chem_isos set
      ! returns nuclide_not_found if name not found
      function get_nuclide_index(nuclei) result(indx)
         use utils_lib, only: integer_dict_lookup
         character(len=*), intent(in) :: nuclei
         integer :: indx, ierr
         if (.not. chem_has_been_initialized) then
            write(*,*) 'must call chem_init before calling any other routine in chem_lib'
            indx = nuclide_not_found
            return
         end if
         ierr = 0
         call integer_dict_lookup(chem_isos_dict, nuclei, indx, ierr)
         if (ierr /= 0) indx = nuclide_not_found
      end function get_nuclide_index

      
      subroutine set_some_isos
         ih1 = get_nuclide_index('h1')
         ih2 = get_nuclide_index('h2')
         ih3 = get_nuclide_index('h3')
         ihe3 = get_nuclide_index('he3')
         ihe4 = get_nuclide_index('he4')
         ili6 = get_nuclide_index('li6')
         ili7 = get_nuclide_index('li7')
         ili8 = get_nuclide_index('li8')
         ibe7 = get_nuclide_index('be7')
         ibe8 = get_nuclide_index('be8')
         ibe9 = get_nuclide_index('be9')
         ibe10 = get_nuclide_index('be10')
         ibe11 = get_nuclide_index('be11')
         ib8 = get_nuclide_index('b8')
         ib10 = get_nuclide_index('b10')
         ib11 = get_nuclide_index('b11')
         ib12 = get_nuclide_index('b12')
         ib13 = get_nuclide_index('b13')
         ib14 = get_nuclide_index('b14')
         ic9 = get_nuclide_index('c9')
         ic10 = get_nuclide_index('c10')
         ic11 = get_nuclide_index('c11')
         ic12 = get_nuclide_index('c12')
         ic13 = get_nuclide_index('c13')
         ic14 = get_nuclide_index('c14')
         ic15 = get_nuclide_index('c15')
         ic16 = get_nuclide_index('c16')
         in12 = get_nuclide_index('n12')
         in13 = get_nuclide_index('n13')
         in14 = get_nuclide_index('n14')
         in15 = get_nuclide_index('n15')
         in16 = get_nuclide_index('n16')
         in17 = get_nuclide_index('n17')
         in18 = get_nuclide_index('n18')
         in19 = get_nuclide_index('n19')
         in20 = get_nuclide_index('n20')
         io13 = get_nuclide_index('o13')
         io14 = get_nuclide_index('o14')
         io15 = get_nuclide_index('o15')
         io16 = get_nuclide_index('o16')
         io17 = get_nuclide_index('o17')
         io18 = get_nuclide_index('o18')
         io19 = get_nuclide_index('o19')
         io20 = get_nuclide_index('o20')
         if15 = get_nuclide_index('f15')
         if16 = get_nuclide_index('f16')
         if17 = get_nuclide_index('f17')
         if18 = get_nuclide_index('f18')
         if19 = get_nuclide_index('f19')
         if20 = get_nuclide_index('f20')
         if21 = get_nuclide_index('f21')
         if22 = get_nuclide_index('f22')
         if23 = get_nuclide_index('f23')
         if24 = get_nuclide_index('f24')
         ine17 = get_nuclide_index('ne17')
         ine18 = get_nuclide_index('ne18')
         ine19 = get_nuclide_index('ne19')
         ine20 = get_nuclide_index('ne20')
         ine21 = get_nuclide_index('ne21')
         ine22 = get_nuclide_index('ne22')
         ine23 = get_nuclide_index('ne23')
         ine24 = get_nuclide_index('ne24')
         ine25 = get_nuclide_index('ne25')
         ine26 = get_nuclide_index('ne26')
         ine27 = get_nuclide_index('ne27')
         ine28 = get_nuclide_index('ne28')
         ina20 = get_nuclide_index('na20')
         ina21 = get_nuclide_index('na21')
         ina22 = get_nuclide_index('na22')
         ina23 = get_nuclide_index('na23')
         ina24 = get_nuclide_index('na24')
         ina25 = get_nuclide_index('na25')
         ina26 = get_nuclide_index('na26')
         ina27 = get_nuclide_index('na27')
         ina28 = get_nuclide_index('na28')
         ina29 = get_nuclide_index('na29')
         ina30 = get_nuclide_index('na30')
         ina31 = get_nuclide_index('na31')
         img20 = get_nuclide_index('mg20')
         img21 = get_nuclide_index('mg21')
         img22 = get_nuclide_index('mg22')
         img23 = get_nuclide_index('mg23')
         img24 = get_nuclide_index('mg24')
         img25 = get_nuclide_index('mg25')
         img26 = get_nuclide_index('mg26')
         img27 = get_nuclide_index('mg27')
         img28 = get_nuclide_index('mg28')
         img29 = get_nuclide_index('mg29')
         img30 = get_nuclide_index('mg30')
         img31 = get_nuclide_index('mg31')
         img32 = get_nuclide_index('mg32')
         img33 = get_nuclide_index('mg33')
         ial22 = get_nuclide_index('al22')
         ial23 = get_nuclide_index('al23')
         ial24 = get_nuclide_index('al24')
         ial25 = get_nuclide_index('al25')
         ial26 = get_nuclide_index('al26')
         ial27 = get_nuclide_index('al27')
         ial28 = get_nuclide_index('al28')
         ial29 = get_nuclide_index('al29')
         ial30 = get_nuclide_index('al30')
         ial31 = get_nuclide_index('al31')
         ial32 = get_nuclide_index('al32')
         ial33 = get_nuclide_index('al33')
         ial34 = get_nuclide_index('al34')
         ial35 = get_nuclide_index('al35')
         isi22 = get_nuclide_index('si22')
         isi23 = get_nuclide_index('si23')
         isi24 = get_nuclide_index('si24')
         isi25 = get_nuclide_index('si25')
         isi26 = get_nuclide_index('si26')
         isi27 = get_nuclide_index('si27')
         isi28 = get_nuclide_index('si28')
         isi29 = get_nuclide_index('si29')
         isi30 = get_nuclide_index('si30')
         isi31 = get_nuclide_index('si31')
         isi32 = get_nuclide_index('si32')
         isi33 = get_nuclide_index('si33')
         isi34 = get_nuclide_index('si34')
         isi35 = get_nuclide_index('si35')
         isi36 = get_nuclide_index('si36')
         isi37 = get_nuclide_index('si37')
         isi38 = get_nuclide_index('si38')
         ip26 = get_nuclide_index('p26')
         ip27 = get_nuclide_index('p27')
         ip28 = get_nuclide_index('p28')
         ip29 = get_nuclide_index('p29')
         ip30 = get_nuclide_index('p30')
         ip31 = get_nuclide_index('p31')
         ip32 = get_nuclide_index('p32')
         ip33 = get_nuclide_index('p33')
         ip34 = get_nuclide_index('p34')
         ip35 = get_nuclide_index('p35')
         ip36 = get_nuclide_index('p36')
         ip37 = get_nuclide_index('p37')
         ip38 = get_nuclide_index('p38')
         ip39 = get_nuclide_index('p39')
         ip40 = get_nuclide_index('p40')
         is27 = get_nuclide_index('s27')
         is28 = get_nuclide_index('s28')
         is29 = get_nuclide_index('s29')
         is30 = get_nuclide_index('s30')
         is31 = get_nuclide_index('s31')
         is32 = get_nuclide_index('s32')
         is33 = get_nuclide_index('s33')
         is34 = get_nuclide_index('s34')
         is35 = get_nuclide_index('s35')
         is36 = get_nuclide_index('s36')
         is37 = get_nuclide_index('s37')
         is38 = get_nuclide_index('s38')
         is39 = get_nuclide_index('s39')
         is40 = get_nuclide_index('s40')
         is41 = get_nuclide_index('s41')
         is42 = get_nuclide_index('s42')
         icl31 = get_nuclide_index('cl31')
         icl32 = get_nuclide_index('cl32')
         icl33 = get_nuclide_index('cl33')
         icl34 = get_nuclide_index('cl34')
         icl35 = get_nuclide_index('cl35')
         icl36 = get_nuclide_index('cl36')
         icl37 = get_nuclide_index('cl37')
         icl38 = get_nuclide_index('cl38')
         icl39 = get_nuclide_index('cl39')
         icl40 = get_nuclide_index('cl40')
         icl41 = get_nuclide_index('cl41')
         icl42 = get_nuclide_index('cl42')
         icl43 = get_nuclide_index('cl43')
         icl44 = get_nuclide_index('cl44')
         iar31 = get_nuclide_index('ar31')
         iar32 = get_nuclide_index('ar32')
         iar33 = get_nuclide_index('ar33')
         iar34 = get_nuclide_index('ar34')
         iar35 = get_nuclide_index('ar35')
         iar36 = get_nuclide_index('ar36')
         iar37 = get_nuclide_index('ar37')
         iar38 = get_nuclide_index('ar38')
         iar39 = get_nuclide_index('ar39')
         iar40 = get_nuclide_index('ar40')
         iar41 = get_nuclide_index('ar41')
         iar42 = get_nuclide_index('ar42')
         iar43 = get_nuclide_index('ar43')
         iar44 = get_nuclide_index('ar44')
         iar45 = get_nuclide_index('ar45')
         iar46 = get_nuclide_index('ar46')
         iar47 = get_nuclide_index('ar47')
         ik35 = get_nuclide_index('k35')
         ik36 = get_nuclide_index('k36')
         ik37 = get_nuclide_index('k37')
         ik38 = get_nuclide_index('k38')
         ik39 = get_nuclide_index('k39')
         ik40 = get_nuclide_index('k40')
         ik41 = get_nuclide_index('k41')
         ik42 = get_nuclide_index('k42')
         ik43 = get_nuclide_index('k43')
         ik44 = get_nuclide_index('k44')
         ik45 = get_nuclide_index('k45')
         ik46 = get_nuclide_index('k46')
         ik47 = get_nuclide_index('k47')
         ica35 = get_nuclide_index('ca35')
         ica36 = get_nuclide_index('ca36')
         ica37 = get_nuclide_index('ca37')
         ica38 = get_nuclide_index('ca38')
         ica39 = get_nuclide_index('ca39')
         ica40 = get_nuclide_index('ca40')
         ica41 = get_nuclide_index('ca41')
         ica42 = get_nuclide_index('ca42')
         ica43 = get_nuclide_index('ca43')
         ica44 = get_nuclide_index('ca44')
         ica45 = get_nuclide_index('ca45')
         ica46 = get_nuclide_index('ca46')
         ica47 = get_nuclide_index('ca47')
         ica48 = get_nuclide_index('ca48')
         ica49 = get_nuclide_index('ca49')
         ica50 = get_nuclide_index('ca50')
         ica51 = get_nuclide_index('ca51')
         ica52 = get_nuclide_index('ca52')
         ica53 = get_nuclide_index('ca53')
         isc40 = get_nuclide_index('sc40')
         isc41 = get_nuclide_index('sc41')
         isc42 = get_nuclide_index('sc42')
         isc43 = get_nuclide_index('sc43')
         isc44 = get_nuclide_index('sc44')
         isc45 = get_nuclide_index('sc45')
         isc46 = get_nuclide_index('sc46')
         isc47 = get_nuclide_index('sc47')
         isc48 = get_nuclide_index('sc48')
         isc49 = get_nuclide_index('sc49')
         isc50 = get_nuclide_index('sc50')
         isc51 = get_nuclide_index('sc51')
         isc52 = get_nuclide_index('sc52')
         isc53 = get_nuclide_index('sc53')
         iti39 = get_nuclide_index('ti39')
         iti40 = get_nuclide_index('ti40')
         iti41 = get_nuclide_index('ti41')
         iti42 = get_nuclide_index('ti42')
         iti43 = get_nuclide_index('ti43')
         iti44 = get_nuclide_index('ti44')
         iti45 = get_nuclide_index('ti45')
         iti46 = get_nuclide_index('ti46')
         iti47 = get_nuclide_index('ti47')
         iti48 = get_nuclide_index('ti48')
         iti49 = get_nuclide_index('ti49')
         iti50 = get_nuclide_index('ti50')
         iti51 = get_nuclide_index('ti51')
         iti52 = get_nuclide_index('ti52')
         iti53 = get_nuclide_index('ti53')
         iti54 = get_nuclide_index('ti54')
         iti55 = get_nuclide_index('ti55')
         iv43 = get_nuclide_index('v43')
         iv44 = get_nuclide_index('v44')
         iv45 = get_nuclide_index('v45')
         iv46 = get_nuclide_index('v46')
         iv47 = get_nuclide_index('v47')
         iv48 = get_nuclide_index('v48')
         iv49 = get_nuclide_index('v49')
         iv50 = get_nuclide_index('v50')
         iv51 = get_nuclide_index('v51')
         iv52 = get_nuclide_index('v52')
         iv53 = get_nuclide_index('v53')
         iv54 = get_nuclide_index('v54')
         iv55 = get_nuclide_index('v55')
         iv56 = get_nuclide_index('v56')
         iv57 = get_nuclide_index('v57')
         icr43 = get_nuclide_index('cr43')
         icr44 = get_nuclide_index('cr44')
         icr45 = get_nuclide_index('cr45')
         icr46 = get_nuclide_index('cr46')
         icr47 = get_nuclide_index('cr47')
         icr48 = get_nuclide_index('cr48')
         icr49 = get_nuclide_index('cr49')
         icr50 = get_nuclide_index('cr50')
         icr51 = get_nuclide_index('cr51')
         icr52 = get_nuclide_index('cr52')
         icr53 = get_nuclide_index('cr53')
         icr54 = get_nuclide_index('cr54')
         icr55 = get_nuclide_index('cr55')
         icr56 = get_nuclide_index('cr56')
         icr57 = get_nuclide_index('cr57')
         icr58 = get_nuclide_index('cr58')
         icr59 = get_nuclide_index('cr59')
         icr60 = get_nuclide_index('cr60')
         icr61 = get_nuclide_index('cr61')
         icr62 = get_nuclide_index('cr62')
         icr63 = get_nuclide_index('cr63')
         icr64 = get_nuclide_index('cr64')
         icr65 = get_nuclide_index('cr65')
         icr66 = get_nuclide_index('cr66')
         imn46 = get_nuclide_index('mn46')
         imn47 = get_nuclide_index('mn47')
         imn48 = get_nuclide_index('mn48')
         imn49 = get_nuclide_index('mn49')
         imn50 = get_nuclide_index('mn50')
         imn51 = get_nuclide_index('mn51')
         imn52 = get_nuclide_index('mn52')
         imn53 = get_nuclide_index('mn53')
         imn54 = get_nuclide_index('mn54')
         imn55 = get_nuclide_index('mn55')
         imn56 = get_nuclide_index('mn56')
         imn57 = get_nuclide_index('mn57')
         imn58 = get_nuclide_index('mn58')
         imn59 = get_nuclide_index('mn59')
         imn60 = get_nuclide_index('mn60')
         imn61 = get_nuclide_index('mn61')
         imn62 = get_nuclide_index('mn62')
         imn63 = get_nuclide_index('mn63')
         ife46 = get_nuclide_index('fe46')
         ife47 = get_nuclide_index('fe47')
         ife48 = get_nuclide_index('fe48')
         ife49 = get_nuclide_index('fe49')
         ife50 = get_nuclide_index('fe50')
         ife51 = get_nuclide_index('fe51')
         ife52 = get_nuclide_index('fe52')
         ife53 = get_nuclide_index('fe53')
         ife54 = get_nuclide_index('fe54')
         ife55 = get_nuclide_index('fe55')
         ife56 = get_nuclide_index('fe56')
         ife57 = get_nuclide_index('fe57')
         ife58 = get_nuclide_index('fe58')
         ife59 = get_nuclide_index('fe59')
         ife60 = get_nuclide_index('fe60')
         ife61 = get_nuclide_index('fe61')
         ife62 = get_nuclide_index('fe62')
         ife63 = get_nuclide_index('fe63')
         ife64 = get_nuclide_index('fe64')
         ife65 = get_nuclide_index('fe65')
         ife66 = get_nuclide_index('fe66')
         ife68 = get_nuclide_index('fe68')
         ico50 = get_nuclide_index('co50')
         ico51 = get_nuclide_index('co51')
         ico52 = get_nuclide_index('co52')
         ico53 = get_nuclide_index('co53')
         ico54 = get_nuclide_index('co54')
         ico55 = get_nuclide_index('co55')
         ico56 = get_nuclide_index('co56')
         ico57 = get_nuclide_index('co57')
         ico58 = get_nuclide_index('co58')
         ico59 = get_nuclide_index('co59')
         ico60 = get_nuclide_index('co60')
         ico61 = get_nuclide_index('co61')
         ico62 = get_nuclide_index('co62')
         ico63 = get_nuclide_index('co63')
         ico64 = get_nuclide_index('co64')
         ico65 = get_nuclide_index('co65')
         ico66 = get_nuclide_index('co66')
         ico67 = get_nuclide_index('co67')
         ini50 = get_nuclide_index('ni50')
         ini51 = get_nuclide_index('ni51')
         ini52 = get_nuclide_index('ni52')
         ini53 = get_nuclide_index('ni53')
         ini54 = get_nuclide_index('ni54')
         ini55 = get_nuclide_index('ni55')
         ini56 = get_nuclide_index('ni56')
         ini57 = get_nuclide_index('ni57')
         ini58 = get_nuclide_index('ni58')
         ini59 = get_nuclide_index('ni59')
         ini60 = get_nuclide_index('ni60')
         ini61 = get_nuclide_index('ni61')
         ini62 = get_nuclide_index('ni62')
         ini63 = get_nuclide_index('ni63')
         ini64 = get_nuclide_index('ni64')
         ini65 = get_nuclide_index('ni65')
         ini66 = get_nuclide_index('ni66')
         ini67 = get_nuclide_index('ni67')
         ini68 = get_nuclide_index('ni68')
         ini69 = get_nuclide_index('ni69')
         ini70 = get_nuclide_index('ni70')
         ini71 = get_nuclide_index('ni71')
         ini72 = get_nuclide_index('ni72')
         ini73 = get_nuclide_index('ni73')
         icu56 = get_nuclide_index('cu56')
         icu57 = get_nuclide_index('cu57')
         icu58 = get_nuclide_index('cu58')
         icu59 = get_nuclide_index('cu59')
         icu60 = get_nuclide_index('cu60')
         icu61 = get_nuclide_index('cu61')
         icu62 = get_nuclide_index('cu62')
         icu63 = get_nuclide_index('cu63')
         icu64 = get_nuclide_index('cu64')
         icu65 = get_nuclide_index('cu65')
         icu66 = get_nuclide_index('cu66')
         icu67 = get_nuclide_index('cu67')
         icu68 = get_nuclide_index('cu68')
         icu69 = get_nuclide_index('cu69')
         icu70 = get_nuclide_index('cu70')
         icu71 = get_nuclide_index('cu71')
         icu72 = get_nuclide_index('cu72')
         izn55 = get_nuclide_index('zn55')
         izn56 = get_nuclide_index('zn56')
         izn57 = get_nuclide_index('zn57')
         izn58 = get_nuclide_index('zn58')
         izn59 = get_nuclide_index('zn59')
         izn60 = get_nuclide_index('zn60')
         izn61 = get_nuclide_index('zn61')
         izn62 = get_nuclide_index('zn62')
         izn63 = get_nuclide_index('zn63')
         izn64 = get_nuclide_index('zn64')
         izn65 = get_nuclide_index('zn65')
         izn66 = get_nuclide_index('zn66')
         izn67 = get_nuclide_index('zn67')
         izn68 = get_nuclide_index('zn68')
         izn69 = get_nuclide_index('zn69')
         izn70 = get_nuclide_index('zn70')
         izn71 = get_nuclide_index('zn71')
         izn72 = get_nuclide_index('zn72')
         izn73 = get_nuclide_index('zn73')
         izn74 = get_nuclide_index('zn74')
         
         iga60 = get_nuclide_index('ga60')
         iga61 = get_nuclide_index('ga61')
         iga62 = get_nuclide_index('ga62')
         iga63 = get_nuclide_index('ga63')
         iga64 = get_nuclide_index('ga64')
         iga65 = get_nuclide_index('ga65')
         iga66 = get_nuclide_index('ga66')
         iga67 = get_nuclide_index('ga67')
         iga68 = get_nuclide_index('ga68')
         iga69 = get_nuclide_index('ga69')
         iga70 = get_nuclide_index('ga70')
         iga71 = get_nuclide_index('ga71')
         iga72 = get_nuclide_index('ga72')
         iga73 = get_nuclide_index('ga73')
         iga74 = get_nuclide_index('ga74')
         iga75 = get_nuclide_index('ga75')
         
         ige59 = get_nuclide_index('ge59')
         ige60 = get_nuclide_index('ge60')
         ige61 = get_nuclide_index('ge61')
         ige62 = get_nuclide_index('ge62')
         ige63 = get_nuclide_index('ge63')
         ige64 = get_nuclide_index('ge64')
         ige65 = get_nuclide_index('ge65')
         ige66 = get_nuclide_index('ge66')
         ige67 = get_nuclide_index('ge67')
         ige68 = get_nuclide_index('ge68')
         ige69 = get_nuclide_index('ge69')
         ige70 = get_nuclide_index('ge70')
         ige71 = get_nuclide_index('ge71')
         ige72 = get_nuclide_index('ge72')
         ige73 = get_nuclide_index('ge73')
         ige74 = get_nuclide_index('ge74')
         ige75 = get_nuclide_index('ge75')
         ige76 = get_nuclide_index('ge76')
         
         ias71 = get_nuclide_index('as71')
         ias72 = get_nuclide_index('as72')
         ias73 = get_nuclide_index('as73')
         ias74 = get_nuclide_index('as74')
         ias75 = get_nuclide_index('as75')
         ias76 = get_nuclide_index('as76')
         ias77 = get_nuclide_index('as77')
         ias78 = get_nuclide_index('as78')
         ias79 = get_nuclide_index('as79')
         
         ise68 = get_nuclide_index('se68')
         ise69 = get_nuclide_index('se69')
         ise70 = get_nuclide_index('se70')
         ise71 = get_nuclide_index('se71')
         ise72 = get_nuclide_index('se72')
         ise73 = get_nuclide_index('se73')
         ise74 = get_nuclide_index('se74')
         ise75 = get_nuclide_index('se75')
         ise76 = get_nuclide_index('se76')

         ikr70 = get_nuclide_index('kr70')
         ikr72 = get_nuclide_index('kr72')
         isr74 = get_nuclide_index('sr74')
         isr76 = get_nuclide_index('sr76')
         izr77 = get_nuclide_index('zr77')
         izr80 = get_nuclide_index('zr80')
         imo82 = get_nuclide_index('mo82')
         imo84 = get_nuclide_index('mo84')
         iru88 = get_nuclide_index('ru88')
         ipd92 = get_nuclide_index('pd92')
         isr75 = get_nuclide_index('sr75')
         iru86 = get_nuclide_index('ru86')
         iru87 = get_nuclide_index('ru87')
         ipd89 = get_nuclide_index('pd89')
         ipd91 = get_nuclide_index('pd91')
         icd93 = get_nuclide_index('cd93')
         icd96 = get_nuclide_index('cd96')
         isn98 = get_nuclide_index('sn98')
         isn100 = get_nuclide_index('sn100')
         isn102 = get_nuclide_index('sn102')
         isn104 = get_nuclide_index('sn104')
         ineut = get_nuclide_index('neut')
         iprot = get_nuclide_index('prot')
         
      end subroutine set_some_isos
         
         
      integer function category_id(cname)
         character (len=*), intent(in)  :: cname 
         ! returns id for the category if there is a matching name
         ! returns 0 otherwise.
         integer :: i, len
         character (len=maxlen_category_name) :: nam
         len = len_trim(cname)
         do i = 1, maxlen_category_name
            if (i <= len) then
               nam(i:i) = cname(i:i)
            else
               nam(i:i) = ' '
            end if
         end do
         do i = 1, num_categories
            if (category_name(i)==nam) then
               category_id = i
               return
            end if
         end do
         category_id = 0
      end function category_id
      
      
      subroutine set_category_names
         integer :: i
         category_name(:) = ''

         category_name(ipp) = 'pp'
         category_name(icno) = 'cno'
         category_name(i3alf) = 'tri_alfa'

         category_name(i_burn_c) = 'burn_c'
         category_name(i_burn_n) = 'burn_n'
         category_name(i_burn_o) = 'burn_o'
         category_name(i_burn_ne) = 'burn_ne'
         category_name(i_burn_na) = 'burn_na'
         category_name(i_burn_mg) = 'burn_mg'
         category_name(i_burn_si) = 'burn_si'
         category_name(i_burn_s) = 'burn_s'
         category_name(i_burn_ar) = 'burn_ar'
         category_name(i_burn_ca) = 'burn_ca'
         category_name(i_burn_ti) = 'burn_ti'
         category_name(i_burn_cr) = 'burn_cr'
         category_name(i_burn_fe) = 'burn_fe'

         category_name(icc) = 'c12_c12'
         category_name(ico) = 'c12_o16'
         category_name(ioo) = 'o16_o16'

         category_name(iphoto) = 'photo'
         category_name(ipnhe4) = 'pnhe4'
         
         category_name(i_ni56_co56) = 'ni56_co56'
         category_name(i_co56_fe56) = 'co56_fe56'

         category_name(iother) = 'other'
         
         do i=1,num_categories
            if (len_trim(category_name(i)) == 0) then
               write(*,*) 'missing name for category', i
               if (i > 1) write(*,*) 'following ' // trim(category_name(i-1))
               flush(6)
               stop 'set_category_names'
            end if
         end do
         
      end subroutine set_category_names


      end module chem_def

