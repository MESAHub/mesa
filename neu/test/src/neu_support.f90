module neu_support
   use neu_def
   use neu_lib
   use math_lib
   use const_def, only: dp
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

end module neu_support

