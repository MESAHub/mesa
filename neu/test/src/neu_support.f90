module neu_support
   use neu_def
   use neu_lib
   use math_lib
   use const_def, only: dp,mev_to_ergs
   use utils_lib, only: mkdir, mesa_error

   implicit none
   
   contains
   
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

                  write(*,'(99(1pe26.16))') logT,logRho,abar,zbar,loss,sources
               end do
            end do
         end do
      end do

   end subroutine do_test_neutrinos


   subroutine do_test_neu_captures()
      real(dp),parameter :: logT_start=9.d0, logT_end=11.0d0
      real(dp),parameter :: mu_start=-5d0, mu_end = 0d0

      integer,parameter :: num_temps=21,  num_mu=21

      real(dp) :: T, logT, mu

      real(dp) :: lec,lpc,lneu,laneu
      real(dp) :: Qec,Qpc,Qneu,Qaneu
      real(dp) :: dratedt, dratedrho, qdt, qdmu
      integer :: ierr,esteps

      integer :: i,j,k


      write(*,'(A)')
      write(*,'(A)')
      write(*,'(A)')
      write(*,*) "Testing neutrino and electron captures"
      write(*,'(A)')
      write(*,'(A)')
      write(*,'(A)')

      call neu_lib_init('../../data/neu_data/neutrino_captures.txt',ierr)
      if(ierr/=0) then
         write(*,*) 'neu capture init failed',ierr
         return
      end if

      do i=1,num_temps

         logT = logT_start + (i-1)*(logT_end-logT_start)/num_temps
         T = exp10(logT)

         do j=1,num_mu
            mu = mu_start + (j-1) * (mu_end-mu_start)/num_mu

            call neu_prot_ec(logT, mu, lec,  dratedt, dratedrho, Qec, Qdt, Qdmu,  ierr)
            call neu_neut_pc(logT, mu, lpc,  dratedt, dratedrho, Qpc, Qdt, Qdmu,  ierr)
            call neu_neut_neu_cap(logT, mu, lneu,  dratedt, dratedrho, Qneu, Qdt, Qdmu,  ierr)
            call neu_prot_aneu_cap(logT, mu, laneu,  dratedt, dratedrho, Qaneu, Qdt, Qdmu,  ierr)

            write(*,'(99(1pe26.16))') logT, mu, lec, lpc, lneu, laneu, Qec, Qpc, Qneu, Qaneu

         end do
      end do

      call neu_lib_shutdown(ierr)


   end subroutine do_test_neu_captures


end module neu_support

