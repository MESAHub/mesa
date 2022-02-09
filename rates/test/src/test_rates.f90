! ***********************************************************************
!
!   Copyright (C) 2011-2020  Bill Paxton & The MESA Team
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
   module test_rates_support
   
   use rates_def
   use rates_lib
   use chem_lib
   use const_lib
   use const_def, only: missing_value
   use math_lib
   use utils_lib, only: mesa_error
   
   implicit none


   contains


   subroutine setup
      use chem_def
      
      integer :: ierr, iso_A, iso_B, iso_C, iso_D
      real(dp) :: W_A, W_B, W_C, W_D, a, b, g_A, g_B, g_C, g_D, Q
      character (len=32) :: my_mesa_dir
      
      include 'formats'
            
      ierr = 0
      
      my_mesa_dir = '../..'         
      call const_init(my_mesa_dir,ierr)     
      if (ierr /= 0) then
         write(*,*) 'const_init failed'
         call mesa_error(__FILE__,__LINE__)
      end if        
      
      call math_init()
      
      call chem_init('isotopes.data', ierr)
      if (ierr /= 0) then
         write(*,*) 'chem_init failed'
         call mesa_error(__FILE__,__LINE__)
      end if

      ! use special weak reaction data in test directory
      
      call rates_init('reactions.list', '', 'rate_tables', &
         .true., &
         .true., 'test_special.states', 'test_special.transitions', &
         '', ierr) 
      if (ierr /= 0) then
         write(*,*) 'rates_init failed'
         call mesa_error(__FILE__,__LINE__)
      end if
      
      call rates_warning_init(.true., 10d0)
      
      call read_raw_rates_records(ierr)
      if (ierr /= 0) then
         write(*,*) 'read_raw_rates_records failed'
         call mesa_error(__FILE__,__LINE__)
      end if

   end subroutine setup


   subroutine do_test_rates()

      integer :: ierr
      type (T_Factors), target :: tf_rec
      type (T_Factors), pointer :: tf
      real(dp) :: logT, temp, raw_rate
      integer :: i, t,ir
      
      integer :: nrates_to_eval
      integer, allocatable :: irs(:)
      real(dp), allocatable :: raw_rates(:)
      real(dp),dimension(9) :: temps
      
      logical, parameter :: dbg = .false.
      
      include 'formats'
      
      write(*,'(A)')
      
      temps = (/6.0d0,6.5d0,7.0d0,7.5d0,8.0d0,8.5d0,9.0d0,9.5d0,10.0d0/)
      
      tf => tf_rec
         
      if (dbg) then
      
         nrates_to_eval = 1
         allocate(irs(nrates_to_eval),raw_rates(nrates_to_eval))
      
         irs(1:nrates_to_eval) = (/ &
            ir_s32_ga_si28     &
            /)
      
      else
      
         nrates_to_eval = num_predefined_reactions
         allocate(irs(nrates_to_eval),raw_rates(nrates_to_eval))
         do i=1, nrates_to_eval
            irs(i) = i    
         end do

      end if
         
      do t=1,size(temps)
         logT = temps(t)
         temp = exp10(logT)
         call eval_tfactors(tf, logT, temp)

         write(*,1) 'logT', logT         
         write(*,1) 'temp', temp
         write(*,'(A)')
         
         raw_rates = missing_value
                     
         call get_raw_rates(nrates_to_eval, irs, temp, tf, raw_rates, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         
         do i=1,nrates_to_eval
            if (raw_rates(i) == missing_value) then
               write(*,*) 'missing value for ' // trim(reaction_Name(irs(i)))
               call mesa_error(__FILE__,__LINE__)
            end if
            write(*,1) trim(reaction_Name(irs(i))), raw_rates(i)
         end do
         write(*,'(A)')
      end do

      write(*,*) 'done'
      write(*,'(A)')
      
   end subroutine do_test_rates 
   

   subroutine test1
      integer :: ierr
      type (T_Factors), target :: tf_rec
      type (T_Factors), pointer :: tf
      real(dp) :: logT, temp, raw_rate, raw_rate1, raw_rate2
      integer :: i, ir
      logical, parameter :: dbg = .false.
      
      include 'formats'
      
      write(*,'(A)')
      write(*,*) 'test1'
      
      tf => tf_rec
      
      temp = 3d9
      logT = log10(temp)
      call eval_tfactors(tf, logT, temp)
      
      write(*,1) 'logT', logT         
      write(*,1) 'temp', temp
      write(*,'(A)')
      
      ir = rates_reaction_id('r_ni56_wk_co56')
      if (ir == 0) then
         write(*,*) 'failed to find rate id'
         call mesa_error(__FILE__,__LINE__)
      end if
      
      call run1 
      raw_rate1 = raw_rate  
      
      temp = 3.000000001d9
      logT = log10(temp)
      call eval_tfactors(tf, logT, temp)
      
      write(*,1) 'logT', logT         
      write(*,1) 'temp', temp
      write(*,1) 'raw_rate1', raw_rate1
      write(*,'(A)')
      
      stop
      
      ir = rates_reaction_id('r_s32_ga_si28')
      call run1   
      raw_rate2 = raw_rate
      
      write(*,1) 'raw_rate2', raw_rate2
      write(*,1) 'raw_rate2-raw_rate1', raw_rate2-raw_rate1

      write(*,*) 'done'
      write(*,'(A)')
      
      contains
      
      subroutine run1
         include 'formats'
         call get_raw_rate(ir, temp, tf, raw_rate, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         write(*,1) trim(reaction_Name(ir)), raw_rate
         write(*,'(A)')
      end subroutine run1
      
   end subroutine test1 
      
   subroutine do_test_rate_table
      integer :: ierr
      type (T_Factors), target :: tf_rec
      type (T_Factors), pointer :: tf
      real(dp) :: logT, temp, raw_rate
      integer :: i, ir
      logical, parameter :: dbg = .false.
      
      include 'formats'
      
      
      write(*,'(A)')
      write(*,*) 'do_test_rate_table'
      
      tf => tf_rec
      
      temp = 9.0d8
      logT = log10(temp)
      call eval_tfactors(tf, logT, temp)
      
      write(*,1) 'logT', logT         
      write(*,1) 'temp', temp
      write(*,'(A)')
      
      ir = rates_reaction_id('r3')
      call run1   

      write(*,*) 'done'
      write(*,'(A)')
      
      contains
      
      subroutine run1
         include 'formats'
         call get_raw_rate(ir, temp, tf, raw_rate, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         write(*,1) trim(reaction_Name(ir)), raw_rate
         write(*,'(A)')
      end subroutine run1
   
   end subroutine do_test_rate_table
   
   subroutine teardown
      call rates_shutdown
   end subroutine teardown

   end module test_rates_support




program test_rates
   
   use test_screen
   use test_weak
   use test_ecapture
   use test_rates_support
   
   implicit none

   call setup
   
   !call do_test_rates(rates_JR_if_available); stop
   !call test1; stop
   
   call do_test_screen
   
   call do_test_weak
   
   call do_test_ecapture
         
   call do_test_rates()
   call do_test_rate_table

   call teardown

end program
