program test_turb
   use math_lib
   use auto_diff
   use const_def
   use turb

   write(*,*) 'CURRENTLY UNTESTED.'

   contains

   subroutine header(text)
      character(len=*), intent(in) :: text

      write(*,'(a)') ' ----------------------------------------------------------------'
      write(*,'(a)') ' '
      write(*,'(a)') ' '//text
      write(*,'(a)') ' '
      write(*,'(a)') ' ----------------------------------------------------------------'

   end subroutine header

   subroutine should_print0(affix, a, z)
      character(len=*), intent(in) :: affix ! to insert ' approximately'
      real(dp), intent(in) :: a, z

      write(*,'(2(a),1(1pd26.16),a,99(1pd26.16))') &
         ' Should print', affix, a, '  :  ', z
      write(*,'(a)') ''
   end subroutine should_print0

   subroutine should_print1(affix, a, b, z)
      character(len=*), intent(in) :: affix ! to insert ' approximately'
      real(dp), intent(in) :: a, b
      type(auto_diff_real_1var_order1), intent(in) :: z

      write(*,'(2(a),2(1pd26.16),a,99(1pd26.16))') &
         ' Should print', affix, a, b, '  :  ', z
      write(*,'(a)') ''
   end subroutine should_print1


   subroutine should_print2(affix, a, b, c, z)
      character(len=*), intent(in) :: affix ! to insert ' approximately'
      real(dp), intent(in) :: a, b, c
      type(auto_diff_real_2var_order1) :: z
      write(*,'(2(a),3(1pd26.16),a,99(1pd26.16))') &
         ' Should print', affix, a, b, c, '  :  ', z
      write(*,'(a)') ''
   end subroutine should_print2


end program test_turb
