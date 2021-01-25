! ***********************************************************************
!
!   Copyright (C) 2011-2019  Bill Paxton & The MESA Team
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


      module condint
      
      use const_def, only: dp
      use math_lib
      use utils_lib, only: mesa_error
      
      implicit none


      integer, parameter :: num_logTs=29, num_logRhos=71, num_logzs=15
      !!! NB: These parameters must be consistent with the table "condtabl.d"!
      logical :: initialized = .false.
      
      real(dp) :: logTs(num_logTs), logRhos(num_logRhos), logzs(num_logzs)
      real(dp), target :: f_ary(4*num_logRhos*num_logTs*num_logzs) ! for bicubic splines
      real(dp), pointer :: f(:,:,:,:)
      integer :: ilinx(num_logzs), iliny(num_logzs)
      
      
      contains
      
      
      subroutine init_potekhin(ierr)
         use kap_def, only: kap_dir
         use interp_2d_lib_db, only: interp_mkbicub_db
         integer, intent(out) :: ierr
         
         character (len=256) :: filename
         integer :: read_err, iz, it, ir, shift
         integer :: ibcxmin                   ! bc flag for x=xmin
         real(dp) :: bcxmin(num_logTs)               ! bc data vs. y at x=xmin
         integer :: ibcxmax                   ! bc flag for x=xmax
         real(dp) :: bcxmax(num_logTs)               ! bc data vs. y at x=xmax
         integer :: ibcymin                   ! bc flag for y=ymin
         real(dp) :: bcymin(num_logRhos)               ! bc data vs. x at y=ymin
         integer :: ibcymax                   ! bc flag for y=ymax
         real(dp) :: bcymax(num_logRhos)               ! bc data vs. x at y=ymax
         real(dp) :: Z
         real(dp), pointer :: f1(:)
         
         include 'formats'
         
         ierr = 0
         if (initialized) return
         
         shift = 4*num_logRhos*num_logTs
         f(1:4,1:num_logRhos,1:num_logTs,1:num_logzs) => f_ary(1:shift*num_logzs)
                  
         filename = trim(kap_dir) // '/condtabl.data'
         open(1,file=trim(filename),status='OLD',iostat=ierr)
         if (ierr /= 0) then
            write(*,*)
            write(*,*)
            write(*,*)
            write(*,*)
            write(*,*)
            write(*,*) 'NOTICE: missing ' // trim(filename)
            write(*,*) 'Please remove the directory mesa/data/kap_data,'
            write(*,*) 'and rerun the mesa ./install script.'
            write(*,*)
            call mesa_error(__FILE__,__LINE__)
         end if
         !print*,'Reading thermal conductivity data...'
         read_err = 0
         read(1,'(A)',iostat=read_err) ! skip the first line
         if (read_err /= 0) ierr = read_err
         do iz = 1, num_logzs
            read (1,*,iostat=read_err) z, (logTs(it),it=1,num_logTs)
            if (read_err /= 0) ierr = read_err
            if (z .eq. 1d0) then
               logzs(iz) = 0d0
            else
               logzs(iz) = log10(z)
            endif
            do ir = 1, num_logRhos
               read(1,*,iostat=read_err) logRhos(ir), (f(1,ir,it,iz),it=1,num_logTs)
               if (read_err /= 0) ierr = read_err
            end do
         end do
         close(1)
         if (ierr /= 0) then
            write(*,*)
            write(*,*)
            write(*,*)
            write(*,*)
            write(*,*)
            write(*,*) 'NOTICE: error trying to read ' // trim(filename)
            write(*,*) 'Please remove the directory mesa/data/kap_data,'
            write(*,*) 'and rerun the mesa ./install script.'
            write(*,*)
            call mesa_error(__FILE__,__LINE__)
         end if
         ! just use "not a knot" bc's at edges of tables
         ibcxmin = 0; bcxmin(1:num_logTs) = 0d0
         ibcxmax = 0; bcxmax(1:num_logTs) = 0d0
         ibcymin = 0; bcymin(1:num_logRhos) = 0d0
         ibcymax = 0; bcymax(1:num_logRhos) = 0d0
         do iz = 1, num_logzs
            f1(1:shift) => f_ary(1+(iz-1)*shift:iz*shift) 
            call interp_mkbicub_db( &
               logRhos, num_logRhos, logTs, num_logTs, f1, num_logRhos, &
               ibcxmin, bcxmin, ibcxmax, bcxmax, &
               ibcymin, bcymin, ibcymax, bcymax, &
               ilinx(iz), iliny(iz), ierr)
            if (ierr /= 0) then
               write(*,*)
               write(*,*)
               write(*,*)
               write(*,*)
               write(*,*)
               write(*,*) 'NOTICE: error in ' // trim(filename)
               write(*,*) 'Please report the problem.'
               write(*,*)
               call mesa_error(__FILE__,__LINE__)
            end if
         end do
         initialized = .true.
      end subroutine init_potekhin
      
      
      subroutine do_electron_conduction( &
            zbar, logRho_in, logT_in, kap, dlogK_dlogRho, dlogK_dlogT, ierr)
         real(dp), intent(in) :: zbar, logRho_in, logT_in
         real(dp), intent(out) :: kap, dlogK_dlogRho, dlogK_dlogT
         integer, intent(out) :: ierr
         
         integer :: iz, iz1, iz2, shift
         real(dp) :: zlog, logRho, logT
         real(dp) :: alfa, beta, logK, &
            logK1, kap1, dlogK1_dlogRho, dlogK1_dlogT, &
            logK2, kap2, dlogK2_dlogRho, dlogK2_dlogT
         real(dp), pointer :: f1(:)
            
         include 'formats'
         
         ierr = 0
         shift = 4*num_logRhos*num_logTs

         logRho = max(logRhos(1),min(logRhos(num_logRhos),logRho_in))
         logT = max(logTs(1),min(logTs(num_logTs),logT_in))
         zlog = max(logzs(1),min(logzs(num_logzs),log10(max(1d-30,zbar))))
         
         if (zlog <= logzs(1)) then ! use 1st
            call get1(1, logK, dlogK_dlogRho, dlogK_dlogT, ierr)
            kap = exp10(logK)
            return
         end if
         
         if (zlog >= logzs(num_logzs)) then ! use last
            call get1(num_logzs, logK, dlogK_dlogRho, dlogK_dlogT, ierr)
            kap = exp10(logK)
            return
         end if
         
         iz1 = -1
         do iz = 2, num_logzs
            if (zlog >= logzs(iz-1) .and. zlog <= logzs(iz)) then
               iz1 = iz-1; iz2 = iz; exit
            end if
         end do
         if (iz1 < 0) then
            write(*,2) 'num_logzs', num_logzs
            do iz = 1, num_logzs
               write(*,2) 'logzs(iz)', iz, logzs(iz)
            end do
            write(*,1) 'zlog', zlog
            write(*,*) 'confusion in do_electron_conduction'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         call get1(iz1, logK1, dlogK1_dlogRho, dlogK1_dlogT, ierr)
         if (ierr /= 0) then
            write(*,*) 'interp failed for iz1 in do_electron_conduction', iz1, logRho, logT
            call mesa_error(__FILE__,__LINE__)
         end if
         
         call get1(iz2, logK2, dlogK2_dlogRho, dlogK2_dlogT, ierr)
         if (ierr /= 0) then
            write(*,*) 'interp failed for iz2 in do_electron_conduction', iz2, logRho, logT
            call mesa_error(__FILE__,__LINE__)
         end if
         
         ! linear interpolation in zlog
         alfa = (zlog - logzs(iz1)) / (logzs(iz2) - logzs(iz1))
         beta = 1d0-alfa
         logK = alfa*logK2 + beta*logK1
         dlogK_dlogRho = alfa*dlogK2_dlogRho + beta*dlogK1_dlogRho
         dlogK_dlogT = alfa*dlogK2_dlogT + beta*dlogK1_dlogT

         kap = exp10(logK)
         
         
         contains
         
         
         subroutine get1(iz, logK, dlogK_dlogRho, dlogK_dlogT, ierr)
            use kap_eval_support, only: Do_Kap_Interpolations
            use const_def
            integer, intent(in) :: iz
            real(dp), intent(out) :: logK, dlogK_dlogRho, dlogK_dlogT
            integer, intent(out) :: ierr
            logical, parameter :: dbg = .false.
            real(dp) :: fval, df_dx, df_dy, &
               logRho0, logRho1, logT0, logT1
            integer :: i_logRho, j_logT, k
            include 'formats'
            ierr = 0
            f1(1:shift) => f_ary(1+(iz-1)*shift:iz*shift) 
            
            if (logRho < logRhos(2)) then
               i_logRho = 1
            else
               i_logRho = num_logRhos-1
               do k=2,num_logRhos-1
                  if (logRho >= logRhos(k) .and. logRho < logRhos(k+1)) then
                     i_logRho = k; exit
                  end if
               end do
            end if
            logRho0 = logRhos(i_logRho)
            logRho1 = logRhos(i_logRho+1)
            
            if (logT < logTs(2)) then
               j_logT = 1
            else
               j_logT = num_logTs-1
               do k=2,num_logTs-1
                  if (logT >= logTs(k) .and. logT < logTs(k+1)) then
                     j_logT = k; exit
                  end if
               end do
            end if
            logT0 = logTs(j_logT)
            logT1 = logTs(j_logT+1)            
            
            call Do_Kap_Interpolations( &
               f1, num_logRhos, num_logTs, i_logRho, j_logT, logRho0, &
               logRho, logRho1, logT0, logT, logT1, fval, df_dx, df_dy)
            if (ierr /= 0) return
            
            ! fval(1) = CK; fval(2) = DRK = dCK/dlogRho; fval(3) = DTK = dCK/dlogT
            ! chi = thermal conductivity, = 10**CK (cgs units)
            ! conduction opacity kappa = 16*boltz_sigma*T^3 / (3*rho*chi)
            ! logK = 3*logT - logRho - CK + log10(16*boltz_sigma/3)
            logK = 3d0*logT - logRho - fval + log10(16d0 * boltz_sigma / 3d0)
            if (dbg) write(*,2) 'do_electron_conduction', iz, logK
            dlogK_dlogRho = -1d0 - df_dx
            dlogK_dlogT = 3d0 - df_dy   
            
         end subroutine get1


      end subroutine do_electron_conduction


      end module condint
