program neu_plotter

   use eos_def
   use eos_lib, only: eosDT_get_new
   use neu_def
   use neu_lib
   use chem_def
   use chem_lib
   use const_lib
   use math_lib
   use num_lib, only : dfridr

   implicit none

   real(dp) :: Rho, T, log10Rho, log10T, Zbar, Abar
   integer :: ierr
   character (len=32) :: my_mesa_dir

   real(dp), parameter :: log10_Tlim = 7.5d0 ! this is what the neu/test uses?
   logical :: flags(num_neu_types) ! true if should include the type of loss
   real(dp) :: loss(num_neu_rvs) ! total from all sources
   real(dp) :: sources(num_neu_types, num_neu_rvs)

   real(dp) :: res1

   integer :: nT, nRho, nZbar, nAbar
   real(dp) :: logT_center, delta_logT, logRho_center, delta_logRho
   real(dp) :: logT_min, logT_max, logRho_min, logRho_max

   real(dp) :: Zbar_center, delta_Zbar, Abar_center, delta_Abar
   real(dp) :: Zbar_min, Zbar_max, Abar_min, Abar_max

   real(dp) :: logT_step, logRho_step, Zbar_step, Abar_step

   integer :: iounit

   integer :: i_var
   integer :: j,k, njs,nks
   real(dp) :: jval, kval

   real(dp) :: var, dvardx_0, dvardx, err, dx_0, xdum
   logical :: doing_partial, doing_dfridr, doing_d_dlnd

   character(len=4) :: xname, yname

   real(dp), parameter :: UNSET = -999
   real(dp), parameter :: min_derivative_error = 1d-4

   namelist /plotter/ &
      nT, nRho, nZbar, nAbar, &
      logT_center, delta_logT, logRho_center, delta_logRho, &
      logT_min, logT_max, logRho_min, logRho_max, &
      Zbar_center, delta_Zbar, Abar_center, delta_Abar, &
      Zbar_min, Zbar_max, Abar_min, Abar_max, &
      xname, yname, doing_partial, doing_dfridr, doing_d_dlnd, &
      i_var


   include 'formats'

   ierr = 0

   my_mesa_dir = '../..'
   call const_init(my_mesa_dir,ierr)
   if (ierr /= 0) then
      write(*,*) 'const_init failed'
      stop 1
   end if

   call math_init()

   logRho_center = UNSET
   logT_center = UNSET
   Zbar_center = UNSET
   Abar_center = UNSET

   delta_logRho = UNSET
   delta_logT = UNSET
   delta_Zbar = UNSET
   delta_Abar = UNSET

   doing_dfridr = .false.
   doing_d_dlnd = .true.

   ! check i_var
   if ((i_var < 0) .or. (i_var > num_neu_types)) then
      write(*,*) 'invalid value of i_var'
      stop
   end if
   
   ! get info from namelist
   open(newunit=iounit, file='inlist_plotter')
   read(iounit, nml=plotter)
   close(iounit)

   if (trim(xname) == trim(yname)) then
      write(*,*) 'xname == yname'
      stop
   end if

   ! file for output
   open(newunit=iounit, file='neu_plotter.dat')

   if (doing_partial) then
      if (doing_d_dlnd) then
         write(*,*) 'plotting partial of log(neu) w.r.t. logRho'
         write(iounit,*) 'log(neu) (partial w.r.t. logRho)'
      else
         write(*,*) 'plotting partial of log(neu) w.r.t. logT'
         write(iounit,*) 'log(neu) (partial w.r.t. logT)'
      end if
   else
      write(*,*) 'plotting log(neu)'
      write(iounit,*) 'log(neu)'
   end if

   select case(xname)
   case('T')
      njs = nT
      write(iounit,*) 'log10(T)'
   case('Rho')
      njs = nRho
      write(iounit,*) 'log10(Rho)'
   case('Zbar')
      njs = nZbar
      write(iounit,*) 'Zbar'
   case('Abar')
      njs = nAbar
      write(iounit,*) 'Abar'
   case default
      write(*,*) 'invalid xname'
      stop
   end select

   select case(yname)
   case('T')
      nks = nT
      write(iounit,*) 'log10(T)'
   case('Rho')
      nks = nRho
      write(iounit,*) 'log10(Rho)'
   case('Zbar')
      nks = nZbar
      write(iounit,*) 'Zbar'
   case('Abar')
      nks = nAbar
      write(iounit,*) 'Abar'
   case default
      write(*,*) 'invalid yname'
      stop
   end select


   if ((logT_center == UNSET) .or. (delta_logT == UNSET)) then
      logT_center = 0.5d0 * (logT_max + logT_min)
      delta_logT = (logT_max - logT_min)
   else
      logT_min = logT_center - delta_logT * 0.5d0
      logT_max = logT_center - delta_logT * 0.5d0
   end if

   if ((logRho_center == UNSET) .or. (delta_logRho == UNSET)) then
      logRho_center = 0.5d0 * (logRho_max + logRho_min)
      delta_logRho = (logRho_max - logRho_min)
   else
      logRho_min = logRho_center - delta_logRho * 0.5d0
      logRho_max = logRho_center - delta_logRho * 0.5d0
   end if

   if ((Zbar_center == UNSET) .or. (delta_Zbar == UNSET)) then
      Zbar_center = 0.5d0 * (Zbar_max + Zbar_min)
      delta_Zbar = (Zbar_max - Zbar_min)
   else
      Zbar_min = Zbar_center - delta_Zbar * 0.5d0
      Zbar_max = Zbar_center - delta_Zbar * 0.5d0
   end if

   if ((Abar_center == UNSET) .or. (delta_Abar == UNSET)) then
      Abar_center = 0.5d0*(Abar_max + Abar_min)
      delta_Abar = (Abar_max - Abar_min)
   else
      Abar_min = Abar_center - delta_Abar * 0.5d0
      Abar_max = Abar_center - delta_Abar * 0.5d0
   end if


   if (nT .gt. 1) then
      logT_step = delta_logT / (nT-1d0)
   else
      logT_step = 0
   end if

   if (nRho .gt. 1) then
      logRho_step = delta_logRho / (nRho-1d0)
   else
      logRho_step = 0
   end if

   if (nZbar .gt. 1) then
      Zbar_step = delta_Zbar / (nZbar-1d0)
   else
      Zbar_step = 0
   end if

   if (nAbar .gt. 1) then
      Abar_step = delta_Abar / (nAbar-1d0)
   else
      Abar_step = 0
   end if


   write(iounit,*) nks, njs

   log10T = logT_center
   T = exp10(log10T)
   log10Rho = logRho_center
   Rho = exp10(log10Rho)
   Zbar = Zbar_center
   Abar = Abar_center

   do j=1,njs !x
      do k=1,nks !y

         select case(xname)
         case('T')
            log10T = logT_min + logT_step*(j - 1)
            T = exp10(log10T)
            jval = log10T
         case('Rho')
            log10Rho = logRho_min + logRho_step*(j - 1)
            rho = exp10(log10Rho)
            jval = log10Rho
         case('Zbar')
            Zbar = Zbar_min + Zbar_step*(j - 1)
            jval = Zbar
         case('Abar')
            Abar = Abar_min + Abar_step*(j - 1)
            jval = Abar
         end select

         select case(yname)
         case('T')
            log10T = logT_min + logT_step*(k - 1)
            T = exp10(log10T)
            kval = log10T
         case('Rho')
            log10Rho = logRho_min + logRho_step*(k - 1)
            rho = exp10(log10Rho)
            kval = log10Rho
         case('Zbar')
            Zbar = Zbar_min + Zbar_step*(k - 1)
            kval = Zbar
         case('Abar')
            Abar = Abar_min + Abar_step*(k - 1)
            kval = Abar
         end select

         ! all flags on
         flags = .true.
         call neu_get(T, log10T, Rho, log10Rho, abar, zbar, log10_Tlim, flags, &
            loss, sources, ierr)

         if (doing_partial) then
            if (doing_d_dlnd) then
               if (i_var == 0) then
                  res1 = Rho * loss(idneu_dRho) / loss(ineu)
               else
                  res1 = Rho * sources(i_var, idneu_dRho) / loss(ineu)
               endif
            else
               if (i_var == 0) then
                  res1 = T * loss(idneu_dT) / loss(ineu)
               else
                  res1 = T * sources(i_var, idneu_dT) / loss(ineu)
               endif
            end if
         else
            if (i_var == 0) then
               res1 = log(loss(ineu))
            else
               res1 = log(sources(i_var, ineu))
            endif
         end if

         if (doing_dfridr) then
            if (i_var == 0) then
               res1 = log(loss(ineu))
            else
               res1 = log(sources(i_var, ineu))
            endif
            if (doing_d_dlnd) then
               if (i_var == 0) then
                  dvardx_0 = Rho * loss(idneu_dRho) / loss(ineu)
               else
                  dvardx_0 = Rho * sources(i_var, idneu_dRho) / loss(ineu)
               endif
            else
               if (i_var == 0) then
                  dvardx_0 = T * loss(idneu_dT) / loss(ineu)
               else
                  dvardx_0 = T * sources(i_var, idneu_dT) / loss(ineu)
               endif
            end if

            dx_0 = 1d-3
            err = 0d0
            dvardx = dfridr(dx_0,dfridr_func,err)
            xdum = (dvardx - dvardx_0)/max(abs(dvardx_0),min_derivative_error)
            res1 = safe_log10(xdum)
         end if


         write(iounit,*) kval, jval, res1
      end do
   end do


   if (ierr /= 0) then
      write(*,*) 'bad result from neu_get'
      stop 1
   end if

contains

   real(dp) function dfridr_func(delta_x) result(val)
      real(dp), intent(in) :: delta_x
      integer :: ierr
      real(dp) :: var, log_var, lnT, lnd
      include 'formats'
      ierr = 0

      lnT = log10T*ln10
      lnd = log10Rho*ln10

      ! must call eos to get new lnfree_e info
      
      if (doing_d_dlnd) then
         log_var = (lnd + delta_x)/ln10
         var = exp10(log_var)
         ! use zbar**2 as z2bar
         call neu_get(T, log10T, var, log_var, abar, zbar, log10_Tlim, flags, &
            loss, sources, ierr)
      else
         log_var = (lnT + delta_x)/ln10
         var = exp10(log_var)
         ! use zbar**2 as z2bar
         call neu_get(var, log_var, Rho, log10Rho, abar, zbar, log10_Tlim, flags, &
            loss, sources, ierr)
      end if

      if (i_var == 0) then
         val = log(loss(ineu))
      else
         val = log(sources(i_var, ineu))
      endif

   end function dfridr_func

end program neu_plotter
