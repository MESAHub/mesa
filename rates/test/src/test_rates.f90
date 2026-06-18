! ***********************************************************************
!
!   Copyright (C) 2011-2020  Bill Paxton & The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! ***********************************************************************

module test_rates_support

   use rates_def
   use rates_lib
   use chem_lib
   use const_def, only: missing_value
   use const_lib, only: const_init
   use math_lib
   use utils_lib, only: mesa_error

   implicit none

contains

   subroutine setup
      use chem_def

      integer :: ierr
      character(len=32) :: my_mesa_dir

      include 'formats'

      ierr = 0

      my_mesa_dir = '../..'
      call const_init(my_mesa_dir, ierr)
      if (ierr /= 0) then
         write (*, *) 'const_init failed'
         call mesa_error(__FILE__, __LINE__)
      end if

      call math_init()

      call chem_init('isotopes.data', ierr)
      if (ierr /= 0) then
         write (*, *) 'chem_init failed'
         call mesa_error(__FILE__, __LINE__)
      end if

      ! use special weak reaction data in test directory

      call rates_init('reactions.list', '', 'rate_tables', &
                      .true., &
                      .true., 'test_special.states', 'test_special.transitions', &
                      '', ierr)
      if (ierr /= 0) then
         write (*, *) 'rates_init failed'
         call mesa_error(__FILE__, __LINE__)
      end if

      call rates_warning_init(.true., 10d0)

      call read_raw_rates_records(ierr)
      if (ierr /= 0) then
         write (*, *) 'read_raw_rates_records failed'
         call mesa_error(__FILE__, __LINE__)
      end if

   end subroutine setup

   subroutine do_test_rates()

      integer :: ierr
      type(T_Factors), target :: tf_rec
      type(T_Factors), pointer :: tf
      real(dp) :: logT, temp
      integer :: i, t

      integer :: nrates_to_eval
      integer, allocatable :: irs(:)
      real(dp), allocatable :: raw_rates(:)
      real(dp), dimension(9) :: temps

      logical, parameter :: dbg = .false.

      include 'formats'

      write (*, '(A)')

      temps = [6.0d0, 6.5d0, 7.0d0, 7.5d0, 8.0d0, 8.5d0, 9.0d0, 9.5d0, 10.0d0]

      tf => tf_rec

      if (dbg) then

         nrates_to_eval = 1
         allocate (irs(nrates_to_eval), raw_rates(nrates_to_eval))

         irs(1:nrates_to_eval) = [ &
                                 ir_s32_ga_si28 &
                                 ]

      else

         nrates_to_eval = num_predefined_reactions
         allocate (irs(nrates_to_eval), raw_rates(nrates_to_eval))
         do i = 1, nrates_to_eval
            irs(i) = i
         end do

      end if

      do t = 1, size(temps)
         logT = temps(t)
         temp = exp10(logT)
         call eval_tfactors(tf, logT, temp)

         write (*, 1) 'logT', logT
         write (*, 1) 'temp', temp
         write (*, '(A)')

         raw_rates = missing_value

         call get_raw_rates(nrates_to_eval, irs, temp, tf, raw_rates, ierr)
         if (ierr /= 0) call mesa_error(__FILE__, __LINE__)

         do i = 1, nrates_to_eval
            if (raw_rates(i) == missing_value) then
               write (*, *) 'missing value for '//trim(reaction_Name(irs(i)))
               call mesa_error(__FILE__, __LINE__)
            end if
            write (*, 1) trim(reaction_Name(irs(i))), raw_rates(i)
         end do
         write (*, '(A)')
      end do

      write (*, *) 'done'
      write (*, '(A)')

   end subroutine do_test_rates

   subroutine test1
      integer :: ierr
      type(T_Factors), target :: tf_rec
      type(T_Factors), pointer :: tf
      real(dp) :: logT, temp, raw_rate, raw_rate1, raw_rate2
      integer :: ir
      logical, parameter :: dbg = .false.

      include 'formats'

      write (*, '(A)')
      write (*, *) 'test1'

      tf => tf_rec

      temp = 3d9
      logT = log10(temp)
      call eval_tfactors(tf, logT, temp)

      write (*, 1) 'logT', logT
      write (*, 1) 'temp', temp
      write (*, '(A)')

      ir = rates_reaction_id('r_ni56_wk_co56')
      if (ir == 0) then
         write (*, *) 'failed to find rate id'
         call mesa_error(__FILE__, __LINE__)
      end if

      call run1
      raw_rate1 = raw_rate

      temp = 3.000000001d9
      logT = log10(temp)
      call eval_tfactors(tf, logT, temp)

      write (*, 1) 'logT', logT
      write (*, 1) 'temp', temp
      write (*, 1) 'raw_rate1', raw_rate1
      write (*, '(A)')

      stop

      ir = rates_reaction_id('r_s32_ga_si28')
      call run1
      raw_rate2 = raw_rate

      write (*, 1) 'raw_rate2', raw_rate2
      write (*, 1) 'raw_rate2-raw_rate1', raw_rate2 - raw_rate1

      write (*, *) 'done'
      write (*, '(A)')

   contains

      subroutine run1
         include 'formats'
         call get_raw_rate(ir, temp, tf, raw_rate, ierr)
         if (ierr /= 0) call mesa_error(__FILE__, __LINE__)
         write (*, 1) trim(reaction_Name(ir)), raw_rate
         write (*, '(A)')
      end subroutine run1

   end subroutine test1

   subroutine do_test_FL_epsnuc_3alf
      real(dp) :: T  ! temperature
      real(dp) :: Rho  ! density
      real(dp) :: Y  ! helium mass fraction
      real(dp) :: UE  ! electron molecular weight
      real(dp) :: eps_nuc  ! eps_nuc in ergs/g/sec
      real(dp) :: deps_nuc_dT  ! partial wrt temperature
      real(dp) :: deps_nuc_dRho  ! partial wrt density
      include 'formats'
      T = 1d7
      Rho = 1d10
      Y = 1
      UE = 2
      call eval_FL_epsnuc_3alf(T, Rho, Y, UE, eps_nuc, deps_nuc_dT, deps_nuc_dRho)
      write (*, 1) 'FL_epsnuc_3alf', eps_nuc
      write (*, '(A)')
   end subroutine do_test_FL_epsnuc_3alf

   subroutine do_test_rate_table
      integer :: ierr
      type(T_Factors), target :: tf_rec
      type(T_Factors), pointer :: tf
      real(dp) :: logT, temp, raw_rate
      integer :: ir
      logical, parameter :: dbg = .false.

      include 'formats'

      write (*, '(A)')
      write (*, *) 'do_test_rate_table'

      tf => tf_rec

      temp = 9.0d8
      logT = log10(temp)
      call eval_tfactors(tf, logT, temp)

      write (*, 1) 'logT', logT
      write (*, 1) 'temp', temp
      write (*, '(A)')

      ir = rates_reaction_id('r3')
      call run1

      write (*, *) 'done'
      write (*, '(A)')

   contains

      subroutine run1
         include 'formats'
         call get_raw_rate(ir, temp, tf, raw_rate, ierr)
         if (ierr /= 0) call mesa_error(__FILE__, __LINE__)
         write (*, 1) trim(reaction_Name(ir)), raw_rate
         write (*, '(A)')
      end subroutine run1

   end subroutine do_test_rate_table

   subroutine do_test_reverse_coefficients
      use chem_def, only: nuclide_data
      use const_def, only: pi, kB=>boltzm, NA=>avo, hbar

      integer, parameter :: nchecks = 6

      integer :: i
      integer, dimension(nchecks) :: chapters, inverse_exps
      type(nuclide_data), pointer :: nuclides
      real(dp) :: fac

      include 'formats'

      ! Check representative reverse rates with Ni-No > 0, = 0, and < 0.
      ! Include both 1 -> 1 and 2 -> 2 cases for Ni-No = 0.
      ! For each case, use the first non-weak reaclib reaction with the requested
      ! chapter and inverse_exp. The printed handle shows which reaction was used.
      ! If this test network does not include a representative for one case,
      ! print a skip line and continue.
      ! The mass ratio always carries a fixed 3/2 power. Ni-No enters through fac and inverse_exp.
      chapters = [r_one_one, r_three_one, r_two_one, r_two_two, r_one_two, r_one_three]
      inverse_exps = [0, 2, 1, 0, -1, -2]

      nuclides => reaclib_rates% nuclides
      if (.not.associated(nuclides)) call mesa_error(__FILE__, __LINE__)

      fac = pow(1d9*kB/(2d0*pi*hbar*hbar*NA),1.5d0)/NA

      write (*, '(A)')
      write (*, *) 'do_test_reverse_coefficients'
      write (*, *) 'inverse_coefficients(1) = log[(m_in/m_out)^(3/2)*(g_in/g_out)*fac^(Ni-No)]'
      write (*, '(a32, 2x, a3, 2x, a3, 2x, a5, 2x, a26, 2x, a)') &
         'rate', 'Ni', 'No', 'Ni-No', 'inverse_coefficients(1)', 'note'

      do i = 1, nchecks
         call check1(chapters(i), inverse_exps(i))
      end do

      write (*, *) 'done'
      write (*, '(A)')

   contains

      subroutine check1(chapter, expected_inverse_exp)
         integer, intent(in) :: chapter, expected_inverse_exp

         integer :: ierr, i, j, lo, hi
         integer :: Ni, No, Nt, inverse_exp
         integer, dimension(max_species_per_reaction) :: ps
         character(len=max_id_length) :: handle
         real(dp), dimension(max_species_per_reaction) :: g, mass
         real(dp) :: tmp, expected_log_coeff, stored_log_coeff, tol

         include 'formats'

         ! Select the first non-weak reaclib reaction in this chapter with the
         ! requested inverse_exp. This keeps the test independent of branch-specific
         ! reaction handles while still checking one representative for each Ni, No case.
         lo = 0
         do i = 1, reaclib_rates% nreactions
            if (reaclib_rates% chapter(i) /= chapter) cycle
            if (reaclib_rates% reaction_flag(i) == 'w' .or. reaclib_rates% reaction_flag(i) == 'e') cycle
            if (reaclib_rates% inverse_exp(i) /= expected_inverse_exp) cycle
            lo = i
            exit
         end do

         if (lo <= 0) then
            Ni = Nin(chapter)
            No = Nout(chapter)
            write (*, '(a32, 2x, i3, 2x, i3, 2x, i5, 2x, a26, 2x, a)') &
               '(none)', Ni, No, expected_inverse_exp, '', 'skipped'
            return
         end if

         handle = trim(reaclib_rates% reaction_handle(lo))
         call reaclib_indices_for_reaction(handle, reaclib_rates, lo, hi, ierr)
         if (ierr /= 0 .or. lo <= 0 .or. hi < lo) then
            write (*, *) 'failed to find reaction '//trim(handle)
            call mesa_error(__FILE__, __LINE__)
         end if

         Ni = Nin(reaclib_rates% chapter(lo))
         No = Nout(reaclib_rates% chapter(lo))
         Nt = Ni + No
         inverse_exp = Ni - No

         ps(1:Nt) = reaclib_rates% pspecies(1:Nt, lo)
         g(1:Nt) = 2d0*nuclides% spin(ps(1:Nt)) + 1d0
         mass(1:Nt) = nuclides% W(ps(1:Nt))

         tmp = product(mass(1:Ni))/product(mass(Ni+1:Nt))
         expected_log_coeff = log(pow(tmp, 1.5d0)*(product(g(1:Ni))/product(g(Ni+1:Nt))))
         if (inverse_exp /= 0) expected_log_coeff = expected_log_coeff + dble(inverse_exp)*log(fac)

         tol = 1d-12*max(1d0, abs(expected_log_coeff))

         do j = lo, hi
            if (reaclib_rates% inverse_exp(j) /= inverse_exp) then
               write (*, *) 'bad inverse_exp for '//trim(handle)
               call mesa_error(__FILE__, __LINE__)
            end if

            stored_log_coeff = reaclib_rates% inverse_coefficients(1, j)
            if (abs(stored_log_coeff - expected_log_coeff) > tol) then
               write (*, *) 'bad inverse coefficient for '//trim(handle)
               write (*, 1) 'stored', stored_log_coeff
               write (*, 1) 'expected', expected_log_coeff
               call mesa_error(__FILE__, __LINE__)
            end if
         end do

         write (*, '(a32, 2x, i3, 2x, i3, 2x, i5, 2x, 1pe26.16, 2x, a)') &
            trim(handle), Ni, No, inverse_exp, expected_log_coeff, ''

      end subroutine check1

   end subroutine do_test_reverse_coefficients

   subroutine do_test2_FL_epsnuc_3alf
      real(dp) :: T  ! temperature
      real(dp) :: Rho  ! density
      real(dp) :: Y  ! helium mass fraction
      real(dp) :: UE  ! electron molecular weight
      real(dp) :: eps_nuc1, eps_nuc2  ! eps_nuc in ergs/g/sec
      real(dp) :: deps_nuc_dT  ! partial wrt temperature
      real(dp) :: deps_nuc_dRho  ! partial wrt density
      real(dp) :: dT, dRho
      include 'formats'
      T = 7.9432823472428218d+07
      dT = T*1d-8
      Rho = 3.1622776601683793d+09
      dRho = Rho*1d-8
      Y = 1
      UE = 2
      call eval_FL_epsnuc_3alf(T, Rho + dRho, Y, UE, eps_nuc1, deps_nuc_dT, deps_nuc_dRho)
      write (*, 1) 'FL_epsnuc_3alf 1', eps_nuc1
      write (*, '(A)')
      call eval_FL_epsnuc_3alf(T, Rho, Y, UE, eps_nuc2, deps_nuc_dT, deps_nuc_dRho)
      write (*, 1) 'FL_epsnuc_3alf 2', eps_nuc2
      write (*, '(A)')
      write (*, 1) 'analytic deps_nuc_dRho', deps_nuc_dRho
      write (*, 1) 'numerical deps_nuc_dRho', (eps_nuc1 - eps_nuc2)/dRho
      write (*, '(A)')
      write (*, 1) 'analytic dlneps_nuc_dlnRho', deps_nuc_dRho*Rho/eps_nuc2
      write (*, 1) 'numerical dlneps_nuc_dlnRho', (eps_nuc1 - eps_nuc2)/dRho*Rho/eps_nuc2
      write (*, '(A)')
   end subroutine do_test2_FL_epsnuc_3alf

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
   call do_test_FL_epsnuc_3alf()
   call do_test_rate_table
   call do_test_reverse_coefficients()

   call teardown

end program test_rates
