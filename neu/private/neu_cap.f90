! ***********************************************************************
!
!   Copyright (C) 2022  The MESA Team
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

   module neu_cap
      use const_def
      use neu_def
      use num_lib
      use num_def
      use math_lib
      use utils_lib, only: mesa_error, is_bad
      use interp_2d_lib_db

      implicit none

      ! Based of Lippuner & Roberts 2017 ApJS https://ui.adsabs.harvard.edu/abs/2017ApJS..233...18L/abstract


      real(dp), parameter :: Qec = (mn-mp) * clight2 ! ergs
      real(dp), parameter :: Qpc = -Qec ! ergs
      real(dp), parameter :: me_erg = me * clight2 ! ergs

      type neu_cap_rates
         real(dp),allocatable,dimension(:) :: logT, mu
         integer :: ntemp, nmu, ntot ! Num of temperature, mu, and total
         real(dp),pointer,dimension(:) :: lecf1=>null(), lpcf1=>null(), lneuf1=>null(), laneuf1=>null()
         real(dp),pointer,dimension(:) :: Qecf1=>null(), Qpcf1=>null(), Qneuf1=>null(), Qaneuf1=>null()
         logical :: init=.false.
      end type neu_cap_rates

      type(neu_cap_rates) :: rate_data

      integer,dimension(6),parameter :: ict = (/1,1,1,0,0,0/)


      contains
  

      subroutine  lambda_ec(logT, mu, rate, dratedt, dratedmu, Q, Qdt, Qdmu, ierr)
         real(dp),intent(in) :: logT, mu ! Temp in K mu electron chemical potential
         real(dp), intent(out) :: rate, dratedt, dratedmu, Q, Qdt, Qdmu
         integer, intent(inout) :: ierr
         real(dp) :: fval(6)

         ierr = 0

         if(.not. rate_data% init) then
            write(*,*) 'Neutrino capture data not initialized'
            ierr=-1
            return
         end if

         call interp_evbicub_db( mu, logT, &
                                 rate_data% mu, rate_data% nMu,&
                                 rate_data% logT, rate_data% nTemp,&
                                 1,1, rate_data% lecf1,rate_data% nMu, &
                                 ict, fval, ierr)
         if(ierr/=0) return

         rate = fval(1)
         dratedt = fval(3)! Fix log to T
         dratedmu = fval(2)


         fval = 0
         call interp_evbicub_db(mu, logT, &
                                 rate_data% mu, rate_data% nMu,&
                                 rate_data% logT, rate_data% nTemp,&
                                 1,1, rate_data% Qecf1,rate_data% nMu, &
                                 ict, fval, ierr)
         if(ierr/=0) return

         Q = fval(1)
         Qdt = fval(3)! Fix log to T
         Qdmu = fval(2)

      end subroutine lambda_ec

      subroutine  lambda_pc(logT, mu, rate, dratedt, dratedmu, Q, Qdt, Qdmu, ierr)
         real(dp),intent(in) :: logT, mu ! Temp in K mu electron chemical potential
         real(dp), intent(out) :: rate, dratedt, dratedmu, Q, Qdt, Qdmu
         integer, intent(inout) :: ierr
         real(dp) :: fval(6)

         ierr = 0
         if(.not. rate_data% init) then
            write(*,*) 'Neutrino capture data not initialized'
            ierr=-1
            return
         end if

         call interp_evbicub_db(mu, logT, &
                                 rate_data% mu, rate_data% nMu,&
                                 rate_data% logT, rate_data% nTemp,&
                                 1,1, rate_data% lpcf1, rate_data% nMu, &
                                 ict, fval, ierr)
         if(ierr/=0) return

         rate = fval(1)
         dratedt = fval(3)! Fix log to T
         dratedmu = fval(2)


         fval = 0
         call interp_evbicub_db(mu, logT, &
                                 rate_data% mu, rate_data% nMu,&
                                 rate_data% logT, rate_data% nTemp,&
                                 1,1, rate_data% Qpcf1,rate_data% nMu, &
                                 ict, fval, ierr)
         if(ierr/=0) return

         Q = fval(1)
         Qdt = fval(3)! Fix log to T
         Qdmu = fval(2)


      end subroutine lambda_pc

      subroutine  lambda_neu(logT, mu, rate, dratedt, dratedmu, Q, Qdt, Qdmu, ierr)
         real(dp),intent(in) :: logT, mu ! Temp in K mu electron chemical potential
         real(dp), intent(out) :: rate, dratedt, dratedmu, Q, Qdt, Qdmu
         integer, intent(inout) :: ierr
         real(dp) :: fval(6)

         ierr = 0
         if(.not. rate_data% init) then
            write(*,*) 'Neutrino capture data not initialized'
            ierr=-1
            return
         end if

         call interp_evbicub_db(mu, logT, &
                                 rate_data% mu, rate_data% nMu,&
                                 rate_data% logT, rate_data% nTemp,&
                                 1,1, rate_data% lneuf1,rate_data% nMu, &
                                 ict, fval, ierr)
         if(ierr/=0) return

         rate = fval(1)
         Qdt = fval(3)! Fix log to T
         Qdmu = fval(2)

         fval = 0
         call interp_evbicub_db(mu, logT, &
                                 rate_data% mu, rate_data% nMu,&
                                 rate_data% logT, rate_data% nTemp,&
                                 1,1, rate_data% Qneuf1,rate_data% nMu, &
                                 ict, fval, ierr)
         if(ierr/=0) return

         Q = fval(1)
         Qdt = fval(3)! Fix log to T
         Qdmu = fval(2)




      end subroutine lambda_neu

      subroutine  lambda_aneu(logT, mu, rate, dratedt, dratedmu, Q, Qdt, Qdmu, ierr)
         real(dp),intent(in) :: logT, mu ! Temp in K mu electron chemical potential
         real(dp), intent(out) :: rate, dratedt, dratedmu, Q, Qdt, Qdmu
         integer, intent(inout) :: ierr
         real(dp) :: fval(6)

         ierr = 0
         if(.not. rate_data% init) then
            write(*,*) 'Neutrino capture data not initialized'
            ierr=-1
            return
         end if

         call interp_evbicub_db(mu, logT, &
                                 rate_data% mu, rate_data% nMu,&
                                 rate_data% logT, rate_data% nTemp,&
                                 1,1, rate_data% laneuf1, rate_data% nMu,&
                                 ict, fval, ierr)
         if(ierr/=0) return

         rate = fval(1)
         dratedt = fval(3)! Fix log to T
         dratedmu = fval(2)


         fval = 0
         call interp_evbicub_db(mu, logT, &
                                 rate_data% mu, rate_data% nMu,&
                                 rate_data% logT, rate_data% nTemp,&
                                 1,1, rate_data% Qaneuf1,rate_data% nMu, &
                                 ict, fval, ierr)
         if(ierr/=0) return

         Q = fval(1)
         Qdt = fval(3)! Fix log to T
         Qdmu = fval(2)

      end subroutine lambda_aneu

  
      subroutine neu_cap_init(filename, ierr)
         use interp_2d_lib_db
         character(len=*), intent(in) :: filename
         integer, intent(inout) :: ierr
         integer :: io,ioerr,ntot, i,j, ind, nMu, nTemp
         character(len=1) :: tmp
         integer :: ibcMumin, ibcMumax, ibcTempmin, ibcTempmax, ilinx, iliny
         real(dp), allocatable :: bcMumin(:),bcMumax(:), bcTempmin(:),  bcTempmax(:)

         integer :: count

         ierr = 0

         open(newunit=io,file=trim(filename),iostat=ioerr)
         if(ioerr/=0) then
            ierr = -1
            write(*,*) "Could not open file ",trim(filename)
            return
         end if

         ! First line is number of lines of data
         ! Second is a header
         ! Third onwards is the data

         read(io,*) tmp, rate_data% ntot, rate_data% ntemp, rate_data% nmu 

         ntot = rate_data% ntot
         nTemp = rate_data% ntemp
         nMu = rate_data% nmu
         
         ! Read empty line
         read(io,'(a)')

         if(allocated(rate_data% logT)) deallocate(rate_data% logT)
         if(allocated(rate_data% mu)) deallocate(rate_data% mu)

         allocate(rate_data% logT(1:nTemp), rate_data% mu(1:nMu),&
                  rate_data% lecf1(1:ntot*4),rate_data% lpcf1(1:ntot*4),&
                  rate_data% lneuf1(1:ntot*4),rate_data% laneuf1(1:ntot*4),&
                  rate_data% Qecf1(1:ntot*4),rate_data% Qpcf1(1:ntot*4),&
                  rate_data% Qneuf1(1:ntot*4),rate_data% Qaneuf1(1:ntot*4)&
                  )

         ! Read rest of data
         do i=1,nTemp
            do j=1, nMu
               ind = 1 + ((j-1) + (i-1)*nMu) * 4
               read(io,*) rate_data% logT(i), rate_data% mu(j),&
                       rate_data% lecf1(ind),rate_data% lpcf1(ind),&
                       rate_data% lneuf1(ind),rate_data% laneuf1(ind),&
                       rate_data% Qecf1(ind),rate_data% Qpcf1(ind),&
                       rate_data% Qneuf1(ind),rate_data% Qaneuf1(ind)
            end do
         end do

         ! use "not a knot" bc's
         allocate(bcMumin(nMu),bcMumax(nMu), bcTempmin(nTemp),  bcTempmax(nTemp))

         ibcMumin = 0; bcMumin(:) = 0d0
         ibcMumax = 0; bcMumax(:) = 0d0
         ibcTempmin = 0; bcTempmin(:) = 0d0
         ibcTempmax = 0; bcTempmax(:) = 0d0

         ! Setup the interpolants
         call interp_mkbicub_db( rate_data% mu(1:nMu), nMu,&
                                 rate_data% logT(1:nTemp), nTemp,&
                                 rate_data% lecf1, nMu, &
                                 ibcMumin, bcMumin, ibcMumax, bcMumax,  &
                                 ibcTempmin, bcTempmin, ibcTempmax, bcTempmax, &
                                 ilinx, iliny, ierr &
         )

         if(ierr/=0) call mesa_error(__FILE__,__LINE__)

         call interp_mkbicub_db( rate_data% mu(1:nMu), nMu,&
                                 rate_data% logT(1:nTemp), nTemp,&
                                 rate_data% lpcf1, nMu, &
                                 ibcMumin, bcMumin, ibcMumax, bcMumax,  &
                                 ibcTempmin, bcTempmin, ibcTempmax, bcTempmax, &
                                 ilinx, iliny, ierr &
         )
         if(ierr/=0) call mesa_error(__FILE__,__LINE__)

         call interp_mkbicub_db( rate_data% mu(1:nMu), nMu,&
                                 rate_data% logT(1:nTemp), nTemp,&
                                 rate_data% lneuf1, nMu, &
                                 ibcMumin, bcMumin, ibcMumax, bcMumax,  &
                                 ibcTempmin, bcTempmin, ibcTempmax, bcTempmax, &
                                 ilinx, iliny, ierr &
         )
         if(ierr/=0) call mesa_error(__FILE__,__LINE__)


         call interp_mkbicub_db( rate_data% mu(1:nMu), nMu,&
                                 rate_data% logT(1:nTemp), nTemp,&
                                 rate_data% laneuf1, nMu, &
                                 ibcMumin, bcMumin, ibcMumax, bcMumax,  &
                                 ibcTempmin, bcTempmin, ibcTempmax, bcTempmax, &
                                 ilinx, iliny, ierr &
         )
         if(ierr/=0) call mesa_error(__FILE__,__LINE__)

         

         call interp_mkbicub_db( rate_data% mu(1:nMu), nMu,&
                                 rate_data% logT(1:nTemp), nTemp,&
                                 rate_data% Qecf1, nMu, &
                                 ibcMumin, bcMumin, ibcMumax, bcMumax,  &
                                 ibcTempmin, bcTempmin, ibcTempmax, bcTempmax, &
                                 ilinx, iliny, ierr &
         )

         if(ierr/=0) call mesa_error(__FILE__,__LINE__)

         call interp_mkbicub_db( rate_data% mu(1:nMu), nMu,&
                                 rate_data% logT(1:nTemp), nTemp,&
                                 rate_data% Qpcf1, nMu, &
                                 ibcMumin, bcMumin, ibcMumax, bcMumax,  &
                                 ibcTempmin, bcTempmin, ibcTempmax, bcTempmax, &
                                 ilinx, iliny, ierr &
         )
         if(ierr/=0) call mesa_error(__FILE__,__LINE__)

         call interp_mkbicub_db( rate_data% mu(1:nMu), nMu,&
                                 rate_data% logT(1:nTemp), nTemp,&
                                 rate_data% Qneuf1, nMu, &
                                 ibcMumin, bcMumin, ibcMumax, bcMumax,  &
                                 ibcTempmin, bcTempmin, ibcTempmax, bcTempmax, &
                                 ilinx, iliny, ierr &
         )
         if(ierr/=0) call mesa_error(__FILE__,__LINE__)


         call interp_mkbicub_db( rate_data% mu(1:nMu), nMu,&
                                 rate_data% logT(1:nTemp), nTemp,&
                                 rate_data% Qaneuf1, nMu, &
                                 ibcMumin, bcMumin, ibcMumax, bcMumax,  &
                                 ibcTempmin, bcTempmin, ibcTempmax, bcTempmax, &
                                 ilinx, iliny, ierr &
         )
         if(ierr/=0) call mesa_error(__FILE__,__LINE__)

         rate_data% init = .true.

      end subroutine neu_cap_init

      subroutine neu_cap_shutdown()

         if(rate_data% init) then

            if(allocated(rate_data% logT)) deallocate(rate_data% logT)
            if(allocated(rate_data% mu)) deallocate(rate_data% mu)

            deallocate(rate_data% lecf1)
            deallocate(rate_data% lpcf1)
            deallocate(rate_data% lneuf1)
            deallocate(rate_data% laneuf1)
            deallocate(rate_data% Qecf1)
            deallocate(rate_data% Qpcf1)
            deallocate(rate_data% Qneuf1)
            deallocate(rate_data% Qaneuf1)

            nullify(rate_data% lecf1, rate_data% lpcf1, rate_data% lneuf1, rate_data% laneuf1)
            nullify(rate_data% Qecf1, rate_data% Qpcf1, rate_data% Qneuf1, rate_data% Qaneuf1)

            rate_data% init = .false.
         end if
      end subroutine neu_cap_shutdown


   end module neu_cap
