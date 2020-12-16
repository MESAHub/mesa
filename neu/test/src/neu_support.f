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


      subroutine do_test_neutrinos(plotting_flag)
      logical, intent(in) :: plotting_flag
      
      integer, parameter :: ionmax=7
      integer, parameter :: ih1 = 1, ihe4 = 2, ic12 = 3, in14 = 4, io16 = 5, ine20 = 6, img24 = 7
      double precision xmass(ionmax),ymass(ionmax),
     1                 aion(ionmax),zion(ionmax),temp,den,abar,zbar,z2bar,
     2                 snu,snudt,snudd,snuda,snudz

!..local variables

      double precision, parameter :: log10_Tlim = 7.9d0
      double precision t, rho, tmin, tmax, snu_alt, plas_alt, phot_alt, pair_alt
      double precision rho_min,rho_max,dt,drho,splas_fxt,spair_fxt,sphot_fxt,sbrem_fxt,sreco_fxt
      double precision :: spair, splas, sphot, sbrem, sreco, snu_fxt, dsnudt,dsnudd,dsnuda,dsnudz
      integer tpoints,rho_points
      integer io,io_rho,io_tmp,io_first,io_last
      
      integer j, i, info

      logical :: flags(num_neu_types) ! true if should include the type

      double precision :: loss(num_neu_rvs) ! total from all sources
      double precision :: sources(num_neu_types, num_neu_rvs)

      logical, parameter :: performance_testing = .false.
      
      logical, parameter :: plot_values = .true.
      
      logical :: do_files
      
      do_files = plotting_flag .and. (.not. performance_testing)

      if (performance_testing) write(*,*) 'doing performance testing'

      flags = .true.
      
!..set the mass fractions, z's and a's of the composition
!..hydrogen
      aion(ih1)  = 1.0d0
      zion(ih1)  = 1.0d0

!..helium
      aion(ihe4)  = 4.0d0
      zion(ihe4)  = 2.0d0

!..carbon 12
      aion(ic12)  = 12.0d0
      zion(ic12)  = 6.0d0

!..nitrogen 14
      aion(in14)  = 14.0d0
      zion(in14)  = 7.0d0

!..oxygen 16
      aion(io16)  = 16.0d0
      zion(io16)  = 8.0d0

!..neon 20
      aion(ine20)  = 20.0d0
      zion(ine20)  = 10.0d0

!..magnesium 24
      aion(img24)  = 24.0d0
      zion(img24)  = 12.0d0

      
      ! mainly H
      xmass(ih1)   = 0.98d0
      xmass(ihe4)  = 0.00d0
      xmass(ic12)  = 0.00d0
      xmass(in14)  = 0.00d0
      xmass(io16)  = 0.00d0
      xmass(ine20) = 0.00d0
      xmass(img24) = 0.02d0
      
      ! mainly C
      xmass(ih1)   = 0.00d0
      xmass(ihe4)  = 0.00d0
      xmass(ic12)  = 0.98d0
      xmass(in14)  = 0.00d0
      xmass(io16)  = 0.00d0
      xmass(ine20) = 0.00d0
      xmass(img24) = 0.02d0

      ! mainly He
      xmass(ih1)   = 0.00d0
      xmass(ihe4)  = 0.98d0
      xmass(ic12)  = 0.00d0
      xmass(in14)  = 0.00d0
      xmass(io16)  = 0.00d0
      xmass(ine20) = 0.00d0
      xmass(img24) = 0.02d0

      ! pure He
      xmass = 0
      xmass(ihe4)  = 1

      !  He, C, and O
      xmass(ih1)   = 0.00d0
      xmass(ihe4)  = 0.50d0
      xmass(ic12)  = 0.25d0
      xmass(in14)  = 0.00d0
      xmass(io16)  = 0.25d0
      xmass(ine20) = 0.00d0
      xmass(img24) = 0.00d0
      

!..get abar and zbar 
      call azbar(xmass,aion,zion,ionmax,ymass,abar,zbar,z2bar)
      
!..set the ranges
      tmax = 9.0d0
      tmin = 6.9d0
      rho_min = -1.0d0
      rho_max = 11.0d0
      tpoints = 101
      rho_points = 101

      if (performance_testing) then
         tpoints = 1501
         rho_points = 1501
      end if
     
      dt = (tmax-tmin)/(tpoints-1)
      drho = (rho_max-rho_min)/(rho_points-1)

      io_rho = 41
      io_tmp = 42
      io_first = 43
      io = io_first

      call mkdir('plot_data')     
      if (do_files) then   
         open(unit=io_rho,file='plot_data/rho.data')
         open(unit=io_tmp,file='plot_data/tmp.data')
         if (plot_values) then
            write(*,*) 'plot values'
            open(unit=io,file='plot_data/snu.data')
            io = io+1; open(unit=io,file='plot_data/pair.data')
            io = io+1; open(unit=io,file='plot_data/plas.data')
            io = io+1; open(unit=io,file='plot_data/phot.data')      
            io = io+1; open(unit=io,file='plot_data/brem.data')
            io = io+1; open(unit=io,file='plot_data/reco.data')   
            io = io+1; open(unit=io,file='plot_data/snudt.data')
            io = io+1; open(unit=io,file='plot_data/snudd.data')
            io = io+1; open(unit=io,file='plot_data/snuda.data')
            io = io+1; open(unit=io,file='plot_data/snudz.data')    
         else    
            write(*,*) 'plot differences'
            open(unit=io,file='plot_data/snu_diff.data')
            io = io+1; open(unit=io,file='plot_data/pair_diff.data')
            io = io+1; open(unit=io,file='plot_data/plas_diff.data')
            io = io+1; open(unit=io,file='plot_data/phot_diff.data')      
            io = io+1; open(unit=io,file='plot_data/brem_diff.data')      
            io = io+1; open(unit=io,file='plot_data/reco_diff.data')      
         end if
      end if

      io_last = io

 01   format(e26.10)
      
      do i=1,tpoints
         t = tmin + dt*(i-1)
         if (do_files) write(io_tmp,01) t
         temp = exp10(T)
         do j=1,rho_points
            rho = rho_min + drho*(j-1)
            if (i .eq. 1 .and. do_files) write(io_rho,01) rho
            den  = exp10(rho)

!..get the neutrino losses
            call neu_get(
     >            temp,t,den,rho,abar,zbar,z2bar,log10_Tlim,flags,loss,sources,info)
            if (info == 0) then
               snu   = loss(ineu)
               snudt = loss(idneu_dT)
               snudd = loss(idneu_dRho)
               snuda = loss(idneu_dabar)
               snudz = loss(idneu_dzbar)
               spair = sources(pair_neu_type, ineu)
               splas = sources(plas_neu_type, ineu)
               sphot = sources(phot_neu_type, ineu)
               sbrem = sources(brem_neu_type, ineu)
               sreco = sources(reco_neu_type, ineu)
            else
               if (.not. plotting_flag) call mesa_error(__FILE__,__LINE__)
               snu=0d0; snudt=0d0; snudd=0d0; snuda=0d0; snudz=0d0
               spair=0d0; splas=0d0; sphot=0d0; sbrem=0d0; sreco=0d0
            end if

            if (do_files) then

               io = io_first; write(io,01) safe_log10(snu)
               io = io+1; write(io,01) safe_log10(spair)
               io = io+1; write(io,01) safe_log10(splas)
               io = io+1; write(io,01) safe_log10(sphot)
               io = io+1; write(io,01) safe_log10(sbrem)
               io = io+1; write(io,01) safe_log10(sreco)
               io = io+1; write(io,01) snudt
               io = io+1; write(io,01) snudd
               io = io+1; write(io,01) snuda
               io = io+1; write(io,01) snudz
               
            end if
            
            if (.not. plotting_flag) then
   
               if (i == 3*tpoints/4 .and. j == 3*rho_points/4) then
 2             format(a20,1pe24.16)
                  write(*,*)
                  write(*,2) 'log10 temp', t
                  write(*,2) 'log10 dens', rho
                  write(*,*)
                  write(*,'(a26)') 'ergs/g/sec'
                  write(*,2) 'log10 plas', safe_log10(splas)
                  write(*,2) 'log10 brem', safe_log10(sbrem)
                  write(*,2) 'log10 phot', safe_log10(sphot)
                  write(*,2) 'log10 pair', safe_log10(spair)
                  write(*,*)
                  write(*,'(a26)') 'derivatives'
                  write(*,2) 'snudt', snudt
                  write(*,2) 'snudd', snudd
                  write(*,*)
               end if

            end if
            
         enddo
      enddo

      if (.not. do_files) return
      
      close(io_rho)
      close(io_tmp)
      do io=io_first,io_last
         close(io)
      end do
      
      
      end subroutine do_test_neutrinos


      subroutine azbar(xmass,aion,zion,ionmax,
     1                 ymass,abar,zbar,z2bar)
      implicit none

!..this routine calculates composition variables for an eos routine

!..input:
!..mass fractions     = xmass(1:ionmax)
!..number of nucleons = aion(1:ionmax)
!..charge of nucleus  = zion(1:ionmax)
!..number of isotopes = ionmax

!..output:
!..molar abundances        = ymass(1:ionmax), 
!..mean number of nucleons = abar
!..mean nucleon charge     = zbar


!..declare
      integer          i,ionmax
      double precision xmass(ionmax),aion(ionmax),zion(ionmax),
     1                 ymass(ionmax),abar,zbar,z2bar,zbarxx,z2barxx,ytot1

      zbarxx  = 0.0d0
      z2barxx  = 0.0d0
      ytot1   = 0.0d0
      do i=1,ionmax
       ymass(i) = xmass(i)/aion(i)
       ytot1    = ytot1 + ymass(i)
       zbarxx   = zbarxx + zion(i) * ymass(i)
       z2barxx   = z2barxx + zion(i)**2 * ymass(i)
      enddo
      abar   = 1.0d0/ytot1
      zbar   = zbarxx * abar
      z2bar   = z2barxx * abar
      return
      end subroutine azbar



      double precision function ifermi12(f)
      implicit none

c..this routine applies a rational function expansion to get the inverse
c..fermi-dirac integral of order 1/2 when it is equal to f.
c..maximum error is 4.19d-9.   reference: antia apjs 84,101 1993

c..declare
      integer          i,m1,k1,m2,k2
      double precision f,an,a1(12),b1(12),a2(12),b2(12),rn,den,ff


c..load the coefficients of the expansion
      data  an,m1,k1,m2,k2 /0.5d0, 4, 3, 6, 5/
      data  (a1(i),i=1,5)/ 1.999266880833d4,   5.702479099336d3,
     1     6.610132843877d2,   3.818838129486d1,
     2     1.0d0/
      data  (b1(i),i=1,4)/ 1.771804140488d4,  -2.014785161019d3,
     1     9.130355392717d1,  -1.670718177489d0/
      data  (a2(i),i=1,7)/-1.277060388085d-2,  7.187946804945d-2, 
     1                    -4.262314235106d-1,  4.997559426872d-1,
     2                    -1.285579118012d0,  -3.930805454272d-1,
     3     1.0d0/
      data  (b2(i),i=1,6)/-9.745794806288d-3,  5.485432756838d-2,
     1                    -3.299466243260d-1,  4.077841975923d-1,
     2                    -1.145531476975d0,  -6.067091689181d-2/


      if (f .lt. 4.0d0) then
       rn  = f + a1(m1)
       do i=m1-1,1,-1
        rn  = rn*f + a1(i)
       enddo
       den = b1(k1+1)
       do i=k1,1,-1
        den = den*f + b1(i)
       enddo
       ifermi12 = log(f * rn/den)

      else
       ff = 1.0d0/pow(f,1.0d0/(1.0d0 + an))
       rn = ff + a2(m2)
       do i=m2-1,1,-1
        rn = rn*ff + a2(i)
       enddo
       den = b2(k2+1)
       do i=k2,1,-1
        den = den*ff + b2(i)
       enddo
       ifermi12 = rn/(den*ff)
      end if
      return
      end function ifermi12






      double precision function zfermim12(x)
      implicit none

c..this routine applies a rational function expansion to get the fermi-dirac
c..integral of order -1/2 evaluated at x. maximum error is 1.23d-12.
c..reference: antia apjs 84,101 1993

c..declare
      integer          i,m1,k1,m2,k2
      double precision x,an,a1(12),b1(12),a2(12),b2(12),rn,den,xx

c..load the coefficients of the expansion
      data  an,m1,k1,m2,k2 /-0.5d0, 7, 7, 11, 11/
      data  (a1(i),i=1,8)/ 1.71446374704454d7,    3.88148302324068d7,
     1                     3.16743385304962d7,    1.14587609192151d7,
     2                     1.83696370756153d6,    1.14980998186874d5,
     3                     1.98276889924768d3,    1.0d0/
      data  (b1(i),i=1,8)/ 9.67282587452899d6,    2.87386436731785d7,
     1                     3.26070130734158d7,    1.77657027846367d7,
     2                     4.81648022267831d6,    6.13709569333207d5,
     3                     3.13595854332114d4,    4.35061725080755d2/
      data (a2(i),i=1,12)/-4.46620341924942d-15, -1.58654991146236d-12,
     1                    -4.44467627042232d-10, -6.84738791621745d-8,
     2                    -6.64932238528105d-6,  -3.69976170193942d-4,
     3                    -1.12295393687006d-2,  -1.60926102124442d-1,
     4                    -8.52408612877447d-1,  -7.45519953763928d-1,
     5                     2.98435207466372d0,    1.0d0/
      data (b2(i),i=1,12)/-2.23310170962369d-15, -7.94193282071464d-13,
     1                    -2.22564376956228d-10, -3.43299431079845d-8,
     2                    -3.33919612678907d-6,  -1.86432212187088d-4,
     3                    -5.69764436880529d-3,  -8.34904593067194d-2,
     4                    -4.78770844009440d-1,  -4.99759250374148d-1,
     5                     1.86795964993052d0,    4.16485970495288d-1/


      if (x .lt. 2.0d0) then
       xx = exp(x)
       rn = xx + a1(m1)
       do i=m1-1,1,-1
        rn = rn*xx + a1(i)
       enddo
       den = b1(k1+1)
       do i=k1,1,-1
        den = den*xx + b1(i)
       enddo
       zfermim12 = xx * rn/den
c..
      else
       xx = 1.0d0/(x*x)
       rn = xx + a2(m2)
       do i=m2-1,1,-1
        rn = rn*xx + a2(i)
       enddo
       den = b2(k2+1)
       do i=k2,1,-1
        den = den*xx + b2(i)
       enddo
       zfermim12 = sqrt(x)*rn/den
      end if
      return
      end function zfermim12


      end module neu_support

