! ***********************************************************************
!
!   Copyright (C) 2010,2020  Aaron Dotter
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

module free_eos_table

   use math_lib
   use const_def
   use eos_def
   use eos_lib

   implicit none

   logical, parameter :: debug = .false.

   !FreeEOS supports 20 elements
   integer, parameter :: Neps = 20
   integer :: ifopt, ifmod, ifion
   integer, parameter :: i_lnRho = num_eos_basic_results + 1
   integer, parameter :: i_dpe = i_lnRho + 1
   integer, parameter :: i_dsp = i_dpe + 1
   integer, parameter :: i_dse = i_dsp + 1
   integer, parameter :: num_results = i_dse

   integer, parameter :: kif = 2 !for P(Rho,T)

   !for MESA 
   integer, parameter ::   h1 =  1
   integer, parameter ::  he4 =  2
   integer, parameter ::  c12 =  3
   integer, parameter ::  n14 =  4
   integer, parameter ::  o16 =  5
   integer, parameter :: ne20 =  6
   integer, parameter :: na23 =  7
   integer, parameter :: mg24 =  8
   integer, parameter :: Al27 =  9
   integer, parameter :: Si28 = 10
   integer, parameter ::  P31 = 11
   integer, parameter ::  S32 = 12
   integer, parameter :: Cl35 = 13
   integer, parameter :: Ar40 = 14
   integer, parameter :: Ca40 = 15
   integer, parameter :: Ti48 = 16
   integer, parameter :: Cr52 = 17
   integer, parameter :: Mn55 = 18
   integer, parameter :: Fe56 = 19
   integer, parameter :: Ni58 = 20

!!!!!

contains 

   !for the 4 basic EOS options
   subroutine free_eos_set_version(eos_version)
      !integer, parameter :: ifopt = 3, ifmod = 1, ifion = -2 !EOS1
      !integer, parameter :: ifopt = 3, ifmod = 1, ifion = -1 !EOS1a
      !integer, parameter :: ifopt = 2, ifmod = 1, ifion = -1 !EOS2
      !integer, parameter :: ifopt = 1, ifmod = 1, ifion =  0 !EOS3
      !integer, parameter :: ifopt = 1, ifmod=101, ifion =  0 !EOS4
      integer, intent(in) :: eos_version
      integer :: my_ifopt, my_ifmod, my_ifion
      select case (eos_version)
      case (1) !EOS1 could also try EOS1a with ifion=-1
         my_ifopt = 3; my_ifmod = 1; my_ifion = -2
      case (2) !EOS2
         my_ifopt = 2; my_ifmod = 1; my_ifion = -1
      case (3) !EOS3
         my_ifopt = 1; my_ifmod = 1; my_ifion =  0
      case (4) !EOS4
         my_ifopt = 1; my_ifmod=101; my_ifion =  0
      case (5) !EOS1a
         my_ifopt = 3; my_ifmod = 1; my_ifion = -1
      case default !default to EOS4
         my_ifopt = 1; my_ifmod=101; my_ifion =  0
      end select
      call free_eos_set_options(my_ifopt,my_ifmod,my_ifion)
   end subroutine free_eos_set_version

   !for complete control over EOS options
   subroutine free_eos_set_options(my_ifopt,my_ifmod,my_ifion)
      integer, intent(in) :: my_ifopt, my_ifmod, my_ifion
      ifopt = my_ifopt; ifmod = my_ifmod; ifion = my_ifion
   end subroutine free_eos_set_options

   subroutine free_eos_eval(logRho,logT,mass_frac,result)
      implicit none
      real(dp), intent(inout) :: logRho
      real(dp), intent(in)  :: logT,mass_frac(Neps)
      real(dp), intent(out) :: result(num_results)
      integer :: iter
      real(dp) :: logf, T, Rho, P, logP, Cf, Cp, Qf, Qp, Sf, St, grada, RTP, Rmue, &
         fh2, fhe2, fhe3, xmu1, xmu3, eta, gamma1, Prad, gamma2, gamma3, h2rat, h2plusrat, &
         lambda, gamma_e, chiRho, chiT, dse, dsp, dpe, dp_dt_constRho, dS_dRho_constT
      real(dp), dimension(Neps) :: atom_wgt, eps
      real(dp), dimension(3) :: degeneracy,pressure,density,energy,enthalpy,entropy

      ! FreeEOS uses an abundance array called eps:
      ! EPS(:) = ( H,He,C,N,O,Ne,Na,Mg,Al,Si,P,S,Cl,Ar,Ca,Ti,Cr,Mn,Fe,Ni)
      ! consisting of 20 elements, each entry of eps is the mass fraction 
      ! of that element divided by its atomic weight
      
      !these are FreeEOS masses, not MESA values
      atom_wgt = [ 1.007825_dp, 4.0026_dp,  12.0111_dp, 14.0067_dp, 15.9994_dp, &
         20.179_dp,   22.98977_dp, 24.305_dp,  26.9815_dp, 28.086_dp, &
         30.9738_dp,  32.06_dp,    35.453_dp,  39.948_dp,  40.08_dp, &
         47.9_dp,     51.996_dp,   54.938_dp,  55.847_dp,  58.71_dp ]

      eps = mass_frac / atom_wgt

      call free_eos( ifopt, ifmod, ifion, kif, eps, Neps, &
         logRho, logT, logf, T, Rho, logRho, P, logP, &
         Cf, Cp, Qf, Qp, Sf, St, grada, RTP, Rmue, fh2, &
         fhe2, fhe3, xmu1, xmu3, eta, gamma1, gamma2, gamma3, &
         h2rat, h2plusrat, lambda, gamma_e, degeneracy, &
         pressure, density, energy, enthalpy, entropy, iter )

      chiRho = pressure(2)
      chiT   = pressure(3)
      Prad  = radiation_pressure(T) !mesa/eos function

      if(entropy(1) < 0._dp) then
         if(debug) write(*,*) 'T, Rho, S', T, Rho, entropy(1)
         grada = -1.0E10_dp
      endif

      dP_dT_constRho = (p/T)*pressure(3)
      dS_dRho_constT = entropy(2)/Rho
      
      dpe = (rho/p)*energy(2) + chiT - 1._dp         !good
      dse = T*(entropy(3)/energy(3)) - 1._dp         !good
      dsp = -rho*rho*(dS_dRho_constT/dP_dT_constRho) - 1._dp !good
      
      result(i_lnPgas) = log10(P - Prad)
      result(i_lnE)    = log10(energy(1))
      result(i_lnS)    = log10(entropy(1))
      result(i_grad_ad)= grada
      result(i_chiRho) = chiRho
      result(i_chiT)   = chiT
      result(i_Cp)     = Cp
      result(i_Cv)     = energy(3)/T !(1/T)*dE/dlnT
      result(i_dE_dRho)= energy(2)/Rho !(1/Rho)*dE/dlnRho
      result(i_dS_dT)  = entropy(3)/T !(1/T)*dS/dlnT
      result(i_dS_dRho)= entropy(2)/Rho !(1/Rho)*dS/dlnRho
      result(i_mu)     = 1_dp/xmu1
      result(i_lnfree_e)= log(xmu3)
      result(i_gamma1) = gamma1
      result(i_gamma3) = gamma3
      result(i_eta)    = eta
      result(i_lnRho)  = log10(Rho)
      result(i_dpe)    = dpe
      result(i_dsp)    = dsp
      result(i_dse)    = dse
   end subroutine free_eos_eval

end module free_eos_table


program make_free_eos_table
   use free_eos_table
   use chem_def
   use chem_lib
   use const_lib
   use utils_lib

   implicit none

   character(len=64) :: eosDT_file, mass_list, data_dir, table_file, arg
   integer :: io, num_logQs, num_logTs, num_logWs, table_version, i
   integer :: eos_version, ierr, num_DT, num_FreeEOS
   real(dp) :: log10T, log10Rho, logT, logRho, mass_frac(Neps)
   real(dp) :: dlog10T, dlog10Q, log10Qmin, log10Qmax, log10Tmin, log10Tmax
   real(dp) :: log10Pgas, logPgas, X, Y, Z

   !for MESA EOS
   integer :: eos_handle
   integer, target :: chem_id_array(neps)
   integer, pointer :: chem_id(:)
   integer, pointer :: net_iso(:) => NULL()
   real(dp) :: abar, zbar, z2bar, z53bar, ye, mass_correction, sumx

   namelist / eos_table / mass_list, table_version, eos_version, &
      log10Tmin, log10Tmax, dlog10T, log10Qmin, log10Qmax, dlog10Q

   data_dir = 'data/eosFreeEOS_data/'

   call read_namelist

   call set_mass_fractions

   call setup_mesa

   if(command_argument_count()>2.and.command_argument_count()<5)then
      call get_command_argument(1,arg)
      read(arg,*) Z

      call get_command_argument(2,arg)
      read(arg,*) X

      call get_command_argument(3,eosDT_file)

      !set mass fractions according to XYZ
      Y = 1._dp - X - Z
      if(debug) write(*,*) 'X,Y,Z=', X, Y, Z
      mass_frac(1) = X*mass_frac(1)
      mass_frac(2) = Y*mass_frac(2)
      mass_frac(3:Neps) = Z*mass_frac(3:Neps)
   else
      write(*,*) './free_eos_table Z X eosDT_file'
      write(*,*) '    Z (real) is mass fraction of metals       '
      write(*,*) '    X (real) is mass fraction of hydrogen     '
      write(*,*) '    eosDT_file is the name of a rho,T eos table; written to data/ dir'
      stop
   endif

   call free_eos_set_version( eos_version )
   write(*,'(a5,a32,3(a3,f5.2))') &
      'file=',trim(eosDT_file),' X=',X,' Y=',Y,' Z=',Z
   table_file = trim(data_dir) // trim(eosDT_file)
   open(newunit=io,file=trim(table_file))
   call write_table(io)
   close(io)

   !final check
   call num_eos_files_loaded(num_DT, num_FreeEOS)

   write(*,*) ' final counts: '
   write(*,*) ' num DT tables loaded = ', num_DT
   write(*,*) ' num FreeEOS tables loaded = ', num_FreeEOS
   
contains

   subroutine read_namelist
      integer :: io_unit

      ! set defaults
      mass_list = 'mass_frac.txt'
      table_version = 51
      eos_version = 1

      !set T range
      log10Tmin = 3d0
      log10Tmax = 8.2d0
      dlog10T = 0.02d0        !default 0.02

      !for eosDT
      log10Qmin = -9d0
      log10Qmax =  4.5d0
      dlog10Q = 0.03d0        !default 0.03

      !now, read namelist
      open(newunit=io_unit,file='inlist',action='read',delim='quote',status='old',iostat=ierr)
      if(ierr/=0) stop 'free_eos_table: problem opening inlist file'
      read(io_unit, nml=eos_table, iostat=ierr) !'
      if(ierr/=0) stop 'free_eos_table: problem reading inlist file'
      close(io_unit)

      num_logTs = 1 + int( (log10Tmax - log10Tmin) / dlog10T )
      num_logQs = 1 + int( (log10Qmax - log10Qmin) / dlog10Q )      

      if(debug)then
         write(*,*) 'dlog10T = ', dlog10T
         write(*,*) 'num_logTs = ', num_logTs
         write(*,*) 'dlog10Q = ', dlog10Q
         write(*,*) 'num_logQs = ', num_logQs
      endif
   end subroutine read_namelist

   subroutine set_mass_fractions
      character(len=2) :: element
      integer :: io_unit
      mass_frac(1) = 1d0 !H
      mass_frac(2) = 1d0 !He
      open(newunit=io_unit,file=trim(mass_list),action='read',status='old',iostat=ierr)
      if(ierr/=0) then 
         write(*,*) 'free_eos_table: problem opening mass fractions list: ', trim(mass_list)
         stop
      endif
      read(io_unit,*,iostat=ierr) !header line
      if(ierr/=0) then
         write(*,*) 'free_eos_table: problem reading mass fractions list: ', trim(mass_list)
         stop
      endif
      do i=3,Neps !read mass fractions of C - Ni
         read(io_unit,*,iostat=ierr) element, mass_frac(i)
         if(ierr/=0) then
            write(*,*) 'free_eos_table: problem reading mass fractions list: ', trim(mass_list)
            stop
         endif
      enddo

      close(io_unit)

      if( abs(sum(mass_frac(3:Neps))-1_dp) > 1.0E-12_dp ) then
         write(*,*) 'free_eos_table: WARNING! Mass fractions of Z do not sum to 1'
      endif

      if(debug)then
         do i=1,Neps
            write(*,*) i, mass_frac(i)
         enddo
      endif
   end subroutine set_mass_fractions

   subroutine write_table(io_unit)
      integer, intent(in) :: io_unit
      real(dp) :: log10Q, mesa_frac, eos_result(num_results)
      real(dp), parameter :: tiny = 1.0E-5_dp
      integer :: mesa_count, iT
      real(dp), parameter :: min_log10T_for_FreeEOS=2.0_dp
      real(dp), allocatable :: results(:,:), logTs(:), logRhos(:), mesa_fracs(:)

      log10Q = 0._dp
      mesa_count = 0
      eos_result = 0._dp
      allocate(results(num_results,num_logTs))
      allocate(logTs(num_logTs), logRhos(num_logTs), mesa_fracs(num_logTs))
      results = 0._dp
      logTs=0._dp
      logRhos=0._dp
      mesa_fracs=0._dp
      
      !write header      
      write(io_unit,'(99(a14))') 'version', 'X', 'Z', 'num logTs', 'logT min', &
         'logT max', 'del logT', 'num logQs', 'logQ min', 'logQ max', 'del logQ'

      write(io_unit,'(i14,2f14.4,2(i10,4x,3(f14.4)))') table_version, X, Z, &
         num_logTs, log10Tmin, log10Tmax, dlog10T, num_logQs, log10Qmin, log10Qmax, dlog10Q

      log10Q = log10Qmin

      do while (log10Q <= log10Qmax)
         log10T = log10Tmax
         iT = num_logTs
         
         !write sub-header
         write(io_unit,'(/,7x,a)') 'logQ = logRho - 2*logT + 12'            
         write(io_unit,'(2x,f14.6/)') log10Q

         !original  '(99(a40,1x))'
         write(io_unit,'(99(a22))') 'logT', &
            'logPgas', 'logE', 'logS', 'chiRho', 'chiT', 'Cp', 'Cv', 'dE_dRho', &
            'dS_dT', 'dS_dRho', 'mu', 'log10_free_e', 'gamma1', 'gamma3', 'grad_ad', &
            'eta', 'MESA', 'logRho', 'dpe', 'dsp', 'dse'
         
         do while (log10T >= log10Tmin)
                  
            if(debug) write(*,*) 'log10Q, log10T=', log10Q, log10T

            log10Rho = log10Q + 2d0*log10T - 12.0d0
            logRho = log10Rho*ln10
            logT   = log10T*ln10

            if(log10T > min_log10T_for_FreeEOS)then
               call free_eos_eval(logRho,logT,mass_frac,eos_result)
            else
               call mesa_eos_eval(logRho,logT,mass_frac,eos_result)
               mesa_frac = 1._dp
            endif

            if(check_for_bad(eos_result))then
               log10Rho = log10Q + 2d0*log10T - 12.0d0
               logRho = log10Rho*ln10
               logT   = log10T*ln10
               if(debug) write(*,*) 'logT   = ', logT
               if(debug) write(*,*) 'logRho = ', logRho
               if(debug) write(*,*) 'free: ', eos_result
               call mesa_eos_eval(logRho, logT, mass_Frac, eos_result)
               if(debug) write(*,*) 'mesa: ', eos_result
               mesa_count = mesa_count + 1
               mesa_frac = 1.0_dp
            else
               mesa_frac = 0.0_dp
            endif

            results(:,iT) = eos_result
            mesa_fracs(iT) = mesa_frac
            logTs(iT) = log10T
            logRhos(iT) = log10Rho
                        
            log10T = log10T - dlog10T
            iT = iT - 1

         enddo

         do iT=1,num_logTs
               !write(io_unit,'(99(1pes40.16e3, 1x))') &
               !write(io_unit,'(f4.2,1p,99(e13.5))') &
            write(io_unit,'(1p,99(e22.14))') &
               logTs(iT),              &
               results(i_lnPgas,iT),   &
               results(i_lnE,iT),      &
               results(i_lnS,iT),      &
               results(i_chiRho,iT),   &
               results(i_chiT,iT),     &
               results(i_Cp,iT),       &
               results(i_Cv,iT),       &
               results(i_dE_dRho,iT),  &
               results(i_dS_dT,iT),    &
               results(i_dS_dRho,iT),  &
               results(i_mu,iT),       &
               results(i_lnfree_e,iT)/ln10, &  !MESA tables are based on OPAL tables, which
               results(i_gamma1,iT),        &  !list  log10(free_e) rather than ln(free_e)
               results(i_gamma3,iT),        & 
               results(i_grad_ad,iT),  &       
               results(i_eta,iT),      &
               mesa_fracs(iT),         &
               logRhos(iT),            &
               results(i_dpe,iT),      &
               results(i_dsp,iT),      &
               results(i_dse,iT)
         enddo
         log10Q = log10Q + dlog10Q
         write(io_unit,*)
      enddo

      if(mesa_count > 0) then
         write(*,*) trim(table_file), ' had this many MESA EOS calls: ', mesa_count
      endif
   end subroutine write_table

   subroutine Setup_MESA
      !..allocate and load the eos tables
      character (len=256) :: eos_file_prefix, my_mesa_dir
      integer :: info
      double precision :: logT_all_HELM, logT_all_OPAL
      logical :: use_cache

      eos_file_prefix = 'mesa'

      info = 0

      my_mesa_dir = '../..'
      call const_init(my_mesa_dir,info)
      if (info /= 0) then
         write(*,*) 'const_init failed'
         stop 1
      end if

      call math_init()

      call chem_init('isotopes.data', info)
      if (info /= 0) then
         write(*,*) 'chem_init failed'
         stop 1
      end if

      use_cache = .true.

      call eos_init( ' ', use_cache, info)
      if (info /= 0) then
         write(*,*) 'failed in eos_init'
         stop 1
      end if

      !eos_handle = alloc_eos_handle(info)
      eos_handle = alloc_eos_handle_using_inlist('inlist', info)
      if (info /= 0) then
         write(*,*) 'failed in alloc_eos_handle'
         stop 1
      end if

      chem_id => chem_id_array
      allocate(net_iso(num_chem_elements))
      net_iso = 0

      chem_id(  h1) =   ih1; net_iso(  ih1) =   h1
      chem_id( he4) =  ihe4; net_iso( ihe4) =  he4
      chem_id( c12) =  ic12; net_iso( ic12) =  c12
      chem_id( n14) =  in14; net_iso( in14) =  n14
      chem_id( o16) =  io16; net_iso( io16) =  o16
      chem_id(ne20) = ine20; net_iso(ine20) = ne20
      chem_id(na23) = ina23; net_iso(ina23) = na23
      chem_id(mg24) = img24; net_iso(img24) = mg24
      chem_id(al27) = ial27; net_iso(ial27) = al27
      chem_id(si28) = isi28; net_iso(isi28) = si28
      chem_id( p31) =  ip31; net_iso( ip31) =  p31
      chem_id( s32) =  is32; net_iso( is32) =  s32
      chem_id(cl35) = icl35; net_iso(icl35) = cl35
      chem_id(ar40) = iar40; net_iso(iar40) = ar40
      chem_id(ca40) = ica40; net_iso(ica40) = ca40
      chem_id(ti48) = iti48; net_iso(iti48) = ti48
      chem_id(cr52) = icr52; net_iso(icr52) = cr52
      chem_id(mn55) = imn55; net_iso(imn55) = mn55
      chem_id(fe56) = ife56; net_iso(ife56) = fe56
      chem_id(ni58) = ini58; net_iso(ini58) = ni58

      call basic_composition_info( Neps, chem_id, mass_frac, X, Y, Z, abar, zbar, z2bar, z53bar, ye, mass_correction, sumx)
   end subroutine Setup_MESA

   logical function check_for_bad(array)
      real(dp), intent(in) :: array(:)
      integer :: i
      do i=1,size(array)
         if(is_bad(array(i)))then
            check_for_bad = .true.
            return
         endif
      enddo
      check_for_bad = array(i_grad_ad) < 0._dp
   end function check_for_bad

   subroutine mesa_eos_eval( logRho0, logT, mass_Frac, eos_result)
      real(dp), intent(in) :: logRho0, logT
      real(dp), intent(in) :: mass_frac(neps)
      real(dp), intent(out) :: eos_result(num_results)

      real(dp) :: T, Pgas, Rho, log10Rho, log10Pgas, log10T, logRho
      real(dp) :: dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas
      real(dp), dimension(num_eos_basic_results) :: &
         res, d_dlnRho_const_T, d_dlnT_const_Rho
      real(dp) :: d_dxa_const_TRho(num_eos_d_dxa_results,neps)
      logical :: off_table
      real(dp), parameter :: logRho_min = -32.23619130191664_dp !-14 * ln10
      integer :: ierr      

      T = exp(logT)
      log10T = logT/ln10

      logRho = max(logRho_min, logRho0)
      Rho = exp(logRho)
      log10Rho = logRho/ln10
         
      call eosDT_get(eos_handle, &
         Neps, chem_id, net_iso, mass_frac, &
         Rho, log10Rho, T, log10T, &
         res, d_dlnRho_const_T, d_dlnT_const_Rho, &
         d_dxa_const_TRho, ierr)

      if (ierr/=0) then !bail to HELM

         logRho = max(logRho_min, logRho0)
         Rho = exp(logRho)
         log10Rho = logRho/ln10
         
         call eosDT_get_component(eos_handle, i_eos_HELM, &
            Neps, chem_id, net_iso, mass_frac, &
            Rho, log10Rho, T, log10T, &
            res, d_dlnRho_const_T, d_dlnT_const_Rho, &
            d_dxa_const_TRho, ierr)
      endif

      if(ierr/=0) then
         write(*,*) 'X = ', X
         write(*,*) 'Z = ', Z
         write(*,*) 'logT = ', log10T
         write(*,*) 'logRho=', log10Rho
         write(*,*) 'res = ', res
         write(*,*) 'ierr= ', ierr
         stop
      endif
         
      eos_result(1:num_eos_basic_results) = res
      eos_result(i_lnRho) = logRho
      eos_result(i_dpe) = 0._dp
      eos_result(i_dsp) = 0._dp
      eos_result(i_dse) = 0._dp

      eos_result(i_lnfree_e) = res(i_lnfree_e)/ln10
   end subroutine mesa_eos_eval

end program make_free_eos_table
