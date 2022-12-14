! ***********************************************************************
!
!   Copyright (C) 2012-2019  Phil Arras & The MESA Team
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

      module create_initial_model

      use star_private_def
      use const_def
      use chem_def

      implicit none

      private
      public :: build_initial_model


      integer :: eos_handle, kap_handle, species

      real(dp) :: X, Z, Y, G, abar, zbar, z2bar, z53bar, ye
      real (dp) :: Tmin,eps,R_try


      integer, parameter :: max_species=1000
      real(dp) :: xa(max_species)
      integer, pointer, dimension(:) :: net_iso, chem_id

      integer, parameter :: nzmax=100000

      type create_star_info
         integer :: nz
         real(dp), dimension(nzmax) :: rg, mg, Pg, Tg, taug, Lg, rhog, intdmTg
         real(dp) :: mass, radius, Teff, luminosity
      end type create_star_info

      integer, parameter :: max_create_star_handles = 3
      type (create_star_info), target, save :: &
         create_star_handles(max_create_star_handles)


      contains


      subroutine get_create_star_ptr(id,cs,ierr)
         integer, intent(in) :: id
         type (create_star_info), pointer :: cs
         integer, intent(out) :: ierr
         include 'formats'
         if (id < 1 .or. id > max_create_star_handles) then
            ierr = -1
            write(*,2) 'bad id for get_create_star_ptr', id
            return
         end if
         cs => create_star_handles(id)
         ierr = 0
      end subroutine get_create_star_ptr



      subroutine build_initial_model(s, ierr)
         use chem_lib, only: basic_composition_info, chem_Xsol
         use adjust_xyz, only: get_xa_for_standard_metals
         use alloc, only: allocate_star_info_arrays
         use star_utils, only: store_r_in_xh, store_rho_in_xh, store_T_in_xh

         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         integer :: initial_zfracs, i_lum, i, j, k, itry, max_try, id0, id1, id2
         real(dp) :: M, R, initial_y, initial_h1, initial_h2, initial_he3, initial_he4, &
            S0, Pc0, rhoc0, e0(2), S1, Pc1, e1(2), S2, Pc2, e2(2), det, dPc, dS, safefac, &
            initial_z, xsol_he3, xsol_he4, mass_correction, mat(2,2), minv(2,2), sumx
         type (create_star_info), pointer :: cs

         include 'formats'

         ierr = 0

         if (s% use_other_build_initial_model) then
            call s% other_build_initial_model(s% id, ierr)
            return
         end if

         id0 = 1; id1 = 2; id2 = 3
         call get_create_star_ptr(id0, cs, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'build_initial_model')

         if (s% nvar_hydro > 4) then
            write(*,*) 'sorry, current build_initial_model only supports the basic 4 vars.'
            ierr = -1
            return
         end if

         M = s% mass_in_gm_for_create_initial_model
         R = s% radius_in_cm_for_create_initial_model

         if (M <= 0) then
            write(*,1) 'must set mass_in_gm_for_create_initial_model', M
            ierr = -1
            return
         end if

         if (R <= 0) then
            write(*,1) 'must set radius_in_cm_for_create_initial_model', R
            ierr = -1
            return
         end if

         initial_z = s% initial_z
         initial_y = s% initial_y
         initial_h1 = max(0d0, min(1d0, 1d0 - (initial_z + initial_y)))
         initial_h2 = 0d0
         xsol_he3 = chem_Xsol('he3')
         xsol_he4 = chem_Xsol('he4')
         initial_he3 = initial_y*xsol_he3/(xsol_he3 + xsol_he4)
         initial_he4 = initial_y*xsol_he4/(xsol_he3 + xsol_he4)
         initial_zfracs = s% initial_zfracs_for_create_initial_model

         chem_id => s% chem_id
         net_iso => s% net_iso
         eos_handle = s% eos_handle
         kap_handle = s% kap_handle
         species = s% species

         call get_xa_for_standard_metals(s, &
            species, s% chem_id, s% net_iso, &
            initial_h1, initial_h2, initial_he3, initial_he4, &
            initial_zfracs, &
            s% initial_dump_missing_heaviest, xa, ierr)
         if (ierr /= 0) then
            write(*,*) 'create initial model failed in get_xa_for_standard_metals'
            return
         end if

         call basic_composition_info( &
            species, s% chem_id, xa(:), x, y, z, abar, zbar, z2bar, z53bar, ye, &
            mass_correction, sumx)

         G = standard_cgrav
         Tmin = 0.d0 ! sets surface isotherm
         eps = s% initial_model_eps ! integration accuracy
         R_try = R   ! used to set grid spacing near the center

         ! init guess
         Pc0=exp10(s% center_logP_1st_try_for_create_initial_model)
         S0=boltzm/mp *  s% entropy_1st_try_for_create_initial_model

         max_try = s% max_tries_for_create_initial_model
         do itry = 1, max_try+1

            if (itry > max_try) then
               write(*,'(a,i10)') &
                  'create initial model failed to converge in allowed number of tries', max_try
               ierr = -1
               return
            end if

            write(*,2) 'create initial model iteration', itry

            Pc1=Pc0*1.01d0
            S1=S0

            Pc2=Pc0
            S2=S0*1.01d0

!$OMP PARALLEL DO PRIVATE(i) SCHEDULE(dynamic,2)
            do i=0,2
               select case(i)
               case(0)
                  call PSerrfunc(id0,Pc0,S0,M,R,e0)
               case(1)
                  call PSerrfunc(id1,Pc1,S1,M,R,e1)
               case(2)
                  call PSerrfunc(id2,Pc2,S2,M,R,e2)
               end select
            end do
!$OMP END PARALLEL DO

            if (abs(e0(1)) < s% abs_e01_tolerance_for_create_initial_model .and. &
                abs(e0(2)) < s% abs_e02_tolerance_for_create_initial_model) exit

            mat(1,1)=(e1(1)-e0(1))/(0.01d0*Pc0)
            mat(1,2)=(e2(1)-e0(1))/(0.01d0*S0)
            mat(2,1)=(e1(2)-e0(2))/(0.01d0*Pc0)
            mat(2,2)=(e2(2)-e0(2))/(0.01d0*S0)
            det = mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)
            minv(1,1)=mat(2,2)/det
            minv(2,2)=mat(1,1)/det
            minv(1,2)=-mat(1,2)/det
            minv(2,1)=-mat(2,1)/det
            dPc = - (minv(1,1)*e0(1)+minv(1,2)*e0(2))
            dS = - (minv(2,1)*e0(1)+minv(2,2)*e0(2))

            safefac = 1.d0 + 10.d0*max( abs(dPc/Pc0), abs(dS/S0) )
            Pc0 = Pc0 + dPc/safefac
            S0 = S0 + dS/safefac

         end do

         s% nz = cs% nz - 1 ! skip center point
         call allocate_star_info_arrays(s, ierr)
         if (ierr /= 0) then
            return
         end if

         s% M_center = 0d0
         s% L_center = 0d0
         s% R_center = 0d0
         s% v_center = 0d0

         s% mstar = cs% mass
         s% star_mass = cs% mass/Msun
         s% xmstar = cs% mass

         i_lum = s% i_lum

         do k=1, s% nz
            i = s% nz - k + 2 ! skip center point
            call store_rho_in_xh(s, k, cs% rhog(i))
            call store_T_in_xh(s, k, cs% Tg(i))
            call store_r_in_xh(s, k, cs% rg(i))
            if (i_lum /= 0) s% xh(i_lum, k) = cs% Lg(i)
            do j=1,species
               s% xa(j,k) = xa(j)
            end do
            s% q(k) = cs% mg(i)/s% xmstar
         end do
         s% dq(s% nz) = s% q(s% nz)
         do k=1, s% nz - 1
            s% dq(k) = s% q(k) - s% q(k+1)
         end do

         write(*,*) 'done build_initial_model'

      end subroutine build_initial_model



      ! output
      !errvec(1)=(mass-M)/M
      !errvec(2)=(radius-R)/R
      subroutine PSerrfunc(id,Pcguess,Sguess,M,R,errvec)
         integer, intent(in) :: id
         real(dp) :: Pcguess,Sguess,M,R,errvec(2)

         integer :: ierr
         real(dp) :: rhoc,Tc,Pc,S
         real(dp) :: P,T,rho,y(3),dydP(3),r1,m1,intdmT1,dy1(3),dy2(3),dy3(3),dy4(3),d_P
         real(dp) :: kap,tau,grav
         logical :: exitnow
         integer :: i, nz, k
         type (create_star_info), pointer :: cs

         call get_create_star_ptr(id, cs, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'PSerrfunc')

         Pc=Pcguess
         S=Sguess


         ! PA. given S and Pc, integrate adiabat outward

         !S = 13.0 * boltzm / mp
         !Pc = 1.0d13

         call get_TRho_from_PS(cs,Pc,S,Tc,rhoc)

         ! initialize grids to zero
         cs% rg=0.d0
         cs% mg=0.d0
         cs% Pg=0.d0
         cs% Tg=0.d0
         cs% rhog=0.d0
         cs% intdmTg=0.d0
         cs% taug=0.d0
         cs% Lg=0.d0

         ! record central point
         i=1
         cs% rg(i)=0.d0
         cs% mg(i)=0.d0
         cs% Pg(i) = Pc
         cs% Tg(i) = Tc
         cs% rhog(i) = rhoc
         cs% intdmTg(i)=0.d0
         cs% taug(i) = 1.d20

         ! compute first point off center
         r1 = 1.d-2 * R_try * min(0.1d0,eps)

         m1 = four_thirds_pi * rhoc * r1*r1*r1
         P = Pc - two_thirds*pi * G*rhoc*rhoc*r1*r1
         intdmT1=m1*Tc
         y=(/r1,m1,intdmT1/)
         call get_TRho_from_PS(cs,P,S,T,rho)

         ! record first point off center
         i=2
         cs% rg(i)=r1
         cs% mg(i)=m1
         cs% Pg(i) = P
         cs% Tg(i) = T
         cs% rhog(i) = rho
         cs% intdmTg(i) = intdmT1
         cs% taug(i) = 1.d20
         cs% Lg(i)=0.d0

         exitnow=.false.
         do

            ! what should stepsize be
            call derivs(cs,P,S,y,dydP)
            d_P=-eps*min( P , abs(y(1)/dydP(1)) , abs(y(2)/dydP(2)) )

            ! 4th order rk step
            dy1=dydP*d_P
            call derivs(cs,P+dP/2d0,S,y+dy1/2d0,dydP)
            dy2=dydP*d_P
            call derivs(cs,P+dP/2d0,S,y+dy2/2d0,dydP)
            dy3=dydP*d_P
            call derivs(cs,P+dP,S,y+dy3,dydP)
            dy4=dydP*d_P
            y=y+dy1/6.d0+dy2/3.d0+dy3/3.d0+dy4/6.d0
            P=P+d_P

            call get_TRho_from_PS(cs,P,S,T,rho)

            call get_kap_from_rhoT(cs,log10(rho),log10(T),kap)
            grav = G*y(2)/(y(1)*y(1))
            tau = kap * P / grav

            i=i+1
            cs% rg(i)=y(1)
            cs% mg(i)=y(2)
            cs% Pg(i) = P
            cs% Tg(i) = T
            cs% rhog(i) = rho
            cs% taug(i) = tau
            cs% intdmTg(i) = y(3)

            if (T<150.d0) then
               print *,"temp too low in integration"
               stop
            endif

            if (tau < 2.d0/3.d0) then
               exitnow=.true.
            endif

            if (exitnow) exit

         end do

         ! interpolate last point to tau=2/3
         nz=i
         cs% nz = nz
         P = cs% Pg(i-1) + (cs% Pg(i)-cs% Pg(i-1))*(2.d0/3.d0-cs% taug(i-1))/(cs% taug(i)-cs% taug(i-1))
         cs% rg(i) = cs% rg(i-1) + (cs% rg(i)-cs% rg(i-1))*(P-cs% Pg(i-1))/(cs% Pg(i)-cs% Pg(i-1))
         cs% mg(i) = cs% mg(i-1) + (cs% mg(i)-cs% mg(i-1))*(P-cs% Pg(i-1))/(cs% Pg(i)-cs% Pg(i-1))
         cs% Tg(i) = cs% Tg(i-1) + (cs% Tg(i)-cs% Tg(i-1))*(P-cs% Pg(i-1))/(cs% Pg(i)-cs% Pg(i-1))
         cs% rhog(i) = cs% rhog(i-1) + (cs% rhog(i)-cs% rhog(i-1))*(P-cs% Pg(i-1))/(cs% Pg(i)-cs% Pg(i-1))
         cs% intdmTg(i) = cs% intdmTg(i-1) + (cs% intdmTg(i)-cs% intdmTg(i-1))*(P-cs% Pg(i-1))/(cs% Pg(i)-cs% Pg(i-1))
         cs% Pg(i)=P
         cs% taug(i)=2.d0/3.d0

         cs% mass = cs% mg(nz)
         cs% radius = cs% rg(nz)
         cs% Teff = cs% Tg(nz)
         cs% luminosity = pi4*pow2(cs% radius)*boltz_sigma*pow4(cs% Teff)
         do k=1,nz
            cs% Lg(k)=cs% luminosity*cs% intdmTg(k)/cs% intdmTg(nz)
         end do
         write(*,"(a20,2x,es15.8)") "log10(L/Lsun)=", log10(cs%luminosity/Lsun)
         write(*,"(a20,2x,es15.8)") "Mass=", cs% mass
         write(*,"(a20,2x,es15.8)") "Radius=", cs% radius
         write(*,*) ''

         errvec(1)=(cs% mass-M)/M
         errvec(2)=(cs% radius-R)/R

      end subroutine PSerrfunc



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      ! y(1)=r
      ! y(2)=m
      ! y(3)=\int_0^m dm T(m)

      subroutine derivs(cs,P,S,y,dydP)
         type (create_star_info), pointer :: cs
         real(dp) :: P,S,y(3),dydP(3)
         real(dp) :: r,m,T,rho,intdmT

         r=y(1)
         m=y(2)
         intdmT=y(3)
         call get_TRho_from_PS(cs,P,S,T,rho)
         !write(*,"(a24,10(2x,es15.8))") "S*mp/kb,lgPc,lgTc,rhoc=",S*mp/boltzm,log10(P),log10(T),rho
         dydP(1)=-r*r/(G*m*rho)
         dydP(2)=-pi4*r*r*r*r/(G*m)
         dydP(3)=dydP(2)*T

      end subroutine


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine get_kap_from_rhoT(cs,logrho,logT,kap)
         use kap_def, only: num_kap_fracs
         use kap_lib
         type (create_star_info), pointer :: cs
         real(dp) :: logrho,logT,kap
         real(dp) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
         real(dp) :: eta, d_eta_dlnRho, d_eta_dlnT
         real(dp) :: dlnkap_dlnRho, dlnkap_dlnT
         real(dp) :: kap_fracs(num_kap_fracs), dlnkap_dxa(species)
         integer :: ierr

         ierr=0

         ! this ignores lnfree and eta
         lnfree_e=0; d_lnfree_e_dlnRho=0; d_lnfree_e_dlnT=0
         eta=0; d_eta_dlnRho=0; d_eta_dlnT=0
         
         call kap_get( &
              kap_handle, species, chem_id, net_iso, xa, &
              logRho, logT, &
              lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
              eta, d_eta_dlnRho, d_eta_dlnT, &
              kap_fracs, kap, dlnkap_dlnRho, dlnkap_dlnT, dlnkap_dxa, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in kap_get get_kap_from_rhoT'
            return
         end if

         if(ierr/=0) then
            write(*,*) 'kap_get failed'
            stop
         endif

      end subroutine get_kap_from_rhoT


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! ignore radiation pressure: P=Pgas here

      subroutine get_TRho_from_PS(cs,P,S,T,rho)
         use eos_lib
         use eos_def
         type (create_star_info), pointer :: cs
         real(dp) :: P,S,T,rho

         real(dp) :: logT_result,log10Rho,dlnRho_dlnPgas_const_T,dlnRho_dlnT_const_Pgas
         real(dp), dimension(num_eos_basic_results) :: &
            res, d_dlnRho_const_T, d_dlnT_const_Rho
         real(dp), dimension(num_eos_d_dxa_results, species) :: d_dxa_const_TRho
         integer, parameter :: max_iter = 100
         integer :: eos_calls,ierr
         real(dp), parameter :: logT_tol = 1.d-6, other_tol = 1.d-6, logT_guess = 4.d0, &
            logT_bnd1= arg_not_provided, logT_bnd2= arg_not_provided, &
            other_at_bnd1= arg_not_provided, other_at_bnd2= arg_not_provided

         call eosPT_get_T( &
            eos_handle, &
            species, chem_id, net_iso, xa, &
            log10(P), i_lnS, log(S), &
            logT_tol, other_tol, max_iter, logT_guess, &
            logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
            logT_result, rho, log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
            res, d_dlnRho_const_T, d_dlnT_const_Rho, &
            d_dxa_const_TRho, &
            eos_calls, ierr)
         if (ierr /=0) then
            print *,"failure in eosPT_get_T"
            stop
         endif
         T = exp10(logT_result)

    ! don't let T get below min temp. use to make a surface isotherm.
         if (T < Tmin) T=Tmin

      end subroutine get_TRho_from_PS


      end module create_initial_model

