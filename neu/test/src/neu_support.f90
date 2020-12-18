      module neu_support
      use neu_def
      use neu_lib
      use math_lib
      use const_def, only: dp
      use utils_lib, only: mkdir, mesa_error

      implicit none
      
      contains

!..tests the neutrino loss rate routine

!..ionmax  = number of isotopes in the network
!..xmass   = mass fractions
!..ymass   = molar fractions
!..aion    = number of nucleons
!..zion    = number of protons



      subroutine do_test_neutrinos()
         real(dp),parameter :: logT_start=6.d0,logT_end=10.5d0
         real(dp),parameter :: logRho_start=6.d0,logRho_end=10.5d0
         real(dp),parameter :: abar_start=1.d0,abar_end=60d0
         real(dp),parameter :: zbar_start=1.d0,zbar_end=26d0

         integer,parameter :: num_temps=10, num_rhos=10, num_abars=10, num_zbars=10

         real(dp) :: T, logT, Rho, logRho, abar, zbar

         integer :: i,j,k,l, info

         logical :: flags(num_neu_types) ! true if should include the type

         real(dp) :: loss(num_neu_rvs) ! total from all sources
         real(dp) :: sources(num_neu_types, num_neu_rvs)

         flags = .true.


         do i=1,num_temps

            logT = logT_start + (i-1)*(logT_end-logT_start)/num_temps
            T = exp10(logT)

            do j=1,num_rhos
               logRho = logRho_start + (j-1)*(logRho_end-logRho_start)/num_rhos
               rho = exp10(logRho)
               do k=1,num_abars
                  abar = abar_start + (k-1)*(abar_end-abar_start)/num_abars

                  do l=1,num_zbars
                     zbar = zbar_start + (l-1)*(zbar_end-zbar_start)/num_zbars

                     call neu_get( T, logT, Rho, logRho,abar,zbar,7.5d0,flags,loss,sources,info)

                     write(*,*) logT,logRho,abar,zbar,loss,sources
                  end do
               end do
            end do
         end do

      end subroutine do_test_neutrinos


!       subroutine do_test_neutrinos(plotting_flag)
!       logical, intent(in) :: plotting_flag
      
!       integer, parameter :: ionmax=7
!       integer, parameter :: ih1 = 1, ihe4 = 2, ic12 = 3, in14 = 4, io16 = 5, ine20 = 6, img24 = 7
!       real(dp) xmass(ionmax),ymass(ionmax), &
!                       aion(ionmax),zion(ionmax),temp,den,abar,zbar,z2bar, &
!                       snu,snudt,snudd,snuda,snudz

! !..local variables

!       real(dp), parameter :: log10_Tlim = 7.9d0
!       real(dp) t, rho, tmin, tmax, snu_alt, plas_alt, phot_alt, pair_alt
!       real(dp) rho_min,rho_max,dt,drho,splas_fxt,spair_fxt,sphot_fxt,sbrem_fxt,sreco_fxt
!       real(dp) :: spair, splas, sphot, sbrem, sreco, snu_fxt, dsnudt,dsnudd,dsnuda,dsnudz
!       integer tpoints,rho_points
!       integer io,io_rho,io_tmp,io_first,io_last
      
!       integer j, i, info

!       logical :: flags(num_neu_types) ! true if should include the type

!       real(dp) :: loss(num_neu_rvs) ! total from all sources
!       real(dp) :: sources(num_neu_types, num_neu_rvs)

!       logical, parameter :: performance_testing = .false.
      
!       logical, parameter :: plot_values = .true.
      
!       logical :: do_files
      
!       do_files = plotting_flag .and. (.not. performance_testing)

!       if (performance_testing) write(*,*) 'doing performance testing'

!       flags = .true.
      
! !..set the mass fractions, z's and a's of the composition
! !..hydrogen
!       aion(ih1)  = 1.0d0
!       zion(ih1)  = 1.0d0

! !..helium
!       aion(ihe4)  = 4.0d0
!       zion(ihe4)  = 2.0d0

! !..carbon 12
!       aion(ic12)  = 12.0d0
!       zion(ic12)  = 6.0d0

! !..nitrogen 14
!       aion(in14)  = 14.0d0
!       zion(in14)  = 7.0d0

! !..oxygen 16
!       aion(io16)  = 16.0d0
!       zion(io16)  = 8.0d0

! !..neon 20
!       aion(ine20)  = 20.0d0
!       zion(ine20)  = 10.0d0

! !..magnesium 24
!       aion(img24)  = 24.0d0
!       zion(img24)  = 12.0d0

      
!       ! mainly H
!       xmass(ih1)   = 0.98d0
!       xmass(ihe4)  = 0.00d0
!       xmass(ic12)  = 0.00d0
!       xmass(in14)  = 0.00d0
!       xmass(io16)  = 0.00d0
!       xmass(ine20) = 0.00d0
!       xmass(img24) = 0.02d0
      
!       ! mainly C
!       xmass(ih1)   = 0.00d0
!       xmass(ihe4)  = 0.00d0
!       xmass(ic12)  = 0.98d0
!       xmass(in14)  = 0.00d0
!       xmass(io16)  = 0.00d0
!       xmass(ine20) = 0.00d0
!       xmass(img24) = 0.02d0

!       ! mainly He
!       xmass(ih1)   = 0.00d0
!       xmass(ihe4)  = 0.98d0
!       xmass(ic12)  = 0.00d0
!       xmass(in14)  = 0.00d0
!       xmass(io16)  = 0.00d0
!       xmass(ine20) = 0.00d0
!       xmass(img24) = 0.02d0

!       ! pure He
!       xmass = 0
!       xmass(ihe4)  = 1

!       !  He, C, and O
!       xmass(ih1)   = 0.00d0
!       xmass(ihe4)  = 0.50d0
!       xmass(ic12)  = 0.25d0
!       xmass(in14)  = 0.00d0
!       xmass(io16)  = 0.25d0
!       xmass(ine20) = 0.00d0
!       xmass(img24) = 0.00d0
      

! !..get abar and zbar 
!       call azbar(xmass,aion,zion,ionmax,ymass,abar,zbar,z2bar)
      
! !..set the ranges
!       tmax = 9.0d0
!       tmin = 6.9d0
!       rho_min = -1.0d0
!       rho_max = 11.0d0
!       tpoints = 101
!       rho_points = 101

!       if (performance_testing) then
!          tpoints = 1501
!          rho_points = 1501
!       end if
     
!       dt = (tmax-tmin)/(tpoints-1)
!       drho = (rho_max-rho_min)/(rho_points-1)

!       io_rho = 41
!       io_tmp = 42
!       io_first = 43
!       io = io_first

!       call mkdir('plot_data')     
!       if (do_files) then   
!          open(unit=io_rho,file='plot_data/rho.data')
!          open(unit=io_tmp,file='plot_data/tmp.data')
!          if (plot_values) then
!             write(*,*) 'plot values'
!             open(unit=io,file='plot_data/snu.data')
!             io = io+1; open(unit=io,file='plot_data/pair.data')
!             io = io+1; open(unit=io,file='plot_data/plas.data')
!             io = io+1; open(unit=io,file='plot_data/phot.data')      
!             io = io+1; open(unit=io,file='plot_data/brem.data')
!             io = io+1; open(unit=io,file='plot_data/reco.data')   
!             io = io+1; open(unit=io,file='plot_data/snudt.data')
!             io = io+1; open(unit=io,file='plot_data/snudd.data')
!             io = io+1; open(unit=io,file='plot_data/snuda.data')
!             io = io+1; open(unit=io,file='plot_data/snudz.data')    
!          else    
!             write(*,*) 'plot differences'
!             open(unit=io,file='plot_data/snu_diff.data')
!             io = io+1; open(unit=io,file='plot_data/pair_diff.data')
!             io = io+1; open(unit=io,file='plot_data/plas_diff.data')
!             io = io+1; open(unit=io,file='plot_data/phot_diff.data')      
!             io = io+1; open(unit=io,file='plot_data/brem_diff.data')      
!             io = io+1; open(unit=io,file='plot_data/reco_diff.data')      
!          end if
!       end if

!       io_last = io

!  01   format(e26.10)
      
!       do i=1,tpoints
!          t = tmin + dt*(i-1)
!          if (do_files) write(io_tmp,01) t
!          temp = exp10(T)
!          do j=1,rho_points
!             rho = rho_min + drho*(j-1)
!             if (i .eq. 1 .and. do_files) write(io_rho,01) rho
!             den  = exp10(rho)

! !..get the neutrino losses
!             call neu_get(
!      >            temp,t,den,rho,abar,zbar,z2bar,log10_Tlim,flags,loss,sources,info)
!             if (info == 0) then
!                snu   = loss(ineu)
!                snudt = loss(idneu_dT)
!                snudd = loss(idneu_dRho)
!                snuda = loss(idneu_dabar)
!                snudz = loss(idneu_dzbar)
!                spair = sources(pair_neu_type, ineu)
!                splas = sources(plas_neu_type, ineu)
!                sphot = sources(phot_neu_type, ineu)
!                sbrem = sources(brem_neu_type, ineu)
!                sreco = sources(reco_neu_type, ineu)
!             else
!                if (.not. plotting_flag) call mesa_error(__FILE__,__LINE__)
!                snu=0d0; snudt=0d0; snudd=0d0; snuda=0d0; snudz=0d0
!                spair=0d0; splas=0d0; sphot=0d0; sbrem=0d0; sreco=0d0
!             end if

!             if (do_files) then

!                io = io_first; write(io,01) safe_log10(snu)
!                io = io+1; write(io,01) safe_log10(spair)
!                io = io+1; write(io,01) safe_log10(splas)
!                io = io+1; write(io,01) safe_log10(sphot)
!                io = io+1; write(io,01) safe_log10(sbrem)
!                io = io+1; write(io,01) safe_log10(sreco)
!                io = io+1; write(io,01) snudt
!                io = io+1; write(io,01) snudd
!                io = io+1; write(io,01) snuda
!                io = io+1; write(io,01) snudz
               
!             end if
            
!             if (.not. plotting_flag) then
   
!                if (i == 3*tpoints/4 .and. j == 3*rho_points/4) then
!  2             format(a20,1pe24.16)
!                   write(*,*)
!                   write(*,2) 'log10 temp', t
!                   write(*,2) 'log10 dens', rho
!                   write(*,*)
!                   write(*,'(a26)') 'ergs/g/sec'
!                   write(*,2) 'log10 plas', safe_log10(splas)
!                   write(*,2) 'log10 brem', safe_log10(sbrem)
!                   write(*,2) 'log10 phot', safe_log10(sphot)
!                   write(*,2) 'log10 pair', safe_log10(spair)
!                   write(*,*)
!                   write(*,'(a26)') 'derivatives'
!                   write(*,2) 'snudt', snudt
!                   write(*,2) 'snudd', snudd
!                   write(*,*)
!                end if

!             end if
            
!          enddo
!       enddo

!       if (.not. do_files) return
      
!       close(io_rho)
!       close(io_tmp)
!       do io=io_first,io_last
!          close(io)
!       end do
      
      
!       end subroutine do_test_neutrinos


!       subroutine azbar(xmass,aion,zion,ionmax, ymass,abar,zbar,z2bar)

! !..this routine calculates composition variables for an eos routine

! !..input:
! !..mass fractions     = xmass(1:ionmax)
! !..number of nucleons = aion(1:ionmax)
! !..charge of nucleus  = zion(1:ionmax)
! !..number of isotopes = ionmax

! !..output:
! !..molar abundances        = ymass(1:ionmax), 
! !..mean number of nucleons = abar
! !..mean nucleon charge     = zbar


! !..declare
!       integer          i,ionmax
!       real(dp) xmass(ionmax),aion(ionmax),zion(ionmax), &
!               ymass(ionmax),abar,zbar,z2bar,zbarxx,z2barxx,ytot1

!          zbarxx  = 0.0d0
!          z2barxx  = 0.0d0
!          ytot1   = 0.0d0
!          do i=1,ionmax
!             ymass(i) = xmass(i)/aion(i)
!             ytot1    = ytot1 + ymass(i)
!             zbarxx   = zbarxx + zion(i) * ymass(i)
!             z2barxx   = z2barxx + zion(i)**2 * ymass(i)
!          enddo
!          abar   = 1.0d0/ytot1
!          zbar   = zbarxx * abar
!          z2bar   = z2barxx * abar

!       end subroutine azbar


      end module neu_support

