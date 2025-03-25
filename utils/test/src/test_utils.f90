
program test_utils

   use utils_def
   use utils_lib
   use const_def, only: dp

   implicit none

   call test_dict

   call test_idict

   call test_token_read

contains

   subroutine test_dict
      type(integer_dict), pointer :: dict

      integer :: value, ierr
      logical :: duplicate

      write (*, '(A)')
      write (*, *) 'test_dict'

      nullify (dict)

      call integer_dict_define_and_report_duplicates(dict, 'c', 3, duplicate, ierr)
      if (ierr /= 0 .or. duplicate) call mesa_error(__FILE__, __LINE__)
      call integer_dict_define_and_report_duplicates(dict, 'a', 1, duplicate, ierr)
      if (ierr /= 0 .or. duplicate) call mesa_error(__FILE__, __LINE__)
      call integer_dict_define_and_report_duplicates(dict, 'd', 4, duplicate, ierr)
      if (ierr /= 0 .or. duplicate) call mesa_error(__FILE__, __LINE__)
      call integer_dict_define_and_report_duplicates(dict, 'b', 0, duplicate, ierr)
      if (ierr /= 0 .or. duplicate) call mesa_error(__FILE__, __LINE__)
      ! redefine some
      call integer_dict_define_and_report_duplicates(dict, 'b', 2, duplicate, ierr)
      if (ierr /= 0 .or. .not. duplicate) call mesa_error(__FILE__, __LINE__)
      call integer_dict_define_and_report_duplicates(dict, 'd', 4, duplicate, ierr)
      if (ierr /= 0 .or. .not. duplicate) call mesa_error(__FILE__, __LINE__)
      call integer_dict_define_and_report_duplicates(dict, 'c', 3, duplicate, ierr)
      if (ierr /= 0 .or. .not. duplicate) call mesa_error(__FILE__, __LINE__)

      call integer_dict_create_hash(dict, ierr)
      if (ierr /= 0) call mesa_error(__FILE__, __LINE__)

      call integer_dict_lookup(dict, 'b', value, ierr)
      if (ierr /= 0) call mesa_error(__FILE__, __LINE__)
      if (value /= 2) call mesa_error(__FILE__, __LINE__)

      call integer_dict_lookup(dict, 'a', value, ierr)
      if (ierr /= 0) call mesa_error(__FILE__, __LINE__)
      if (value /= 1) call mesa_error(__FILE__, __LINE__)

      call integer_dict_lookup(dict, 'd', value, ierr)
      if (ierr /= 0) call mesa_error(__FILE__, __LINE__)
      if (value /= 4) call mesa_error(__FILE__, __LINE__)

      call integer_dict_lookup(dict, 'bogus', value, ierr)
      if (ierr == 0) call mesa_error(__FILE__, __LINE__)
      ierr = 0

      call integer_dict_free(dict)

      write (*, *) 'okay'
      write (*, '(A)')

   end subroutine test_dict

   subroutine test_idict
      type(integer_idict), pointer :: idict

      integer :: value, ierr
      logical :: duplicate

      write (*, '(A)')
      write (*, *) 'test_idict'

      nullify (idict)

      call integer_idict_define_and_report_duplicates(idict, 196, 48, 3, duplicate, ierr)
      if (ierr /= 0 .or. duplicate) call mesa_error(__FILE__, __LINE__)
      call integer_idict_define_and_report_duplicates(idict, 1547, 974, 1, duplicate, ierr)
      if (ierr /= 0 .or. duplicate) call mesa_error(__FILE__, __LINE__)
      call integer_idict_define_and_report_duplicates(idict, 592, 8, 4, duplicate, ierr)
      if (ierr /= 0 .or. duplicate) call mesa_error(__FILE__, __LINE__)
      call integer_idict_define_and_report_duplicates(idict, -51, 885, 0, duplicate, ierr)
      if (ierr /= 0 .or. duplicate) call mesa_error(__FILE__, __LINE__)
      ! redefine some
      call integer_idict_define_and_report_duplicates(idict, -51, 885, 2, duplicate, ierr)
      if (ierr /= 0 .or. .not. duplicate) call mesa_error(__FILE__, __LINE__)
      call integer_idict_define_and_report_duplicates(idict, 592, 8, 4, duplicate, ierr)
      if (ierr /= 0 .or. .not. duplicate) call mesa_error(__FILE__, __LINE__)
      call integer_idict_define_and_report_duplicates(idict, 196, 48, 3, duplicate, ierr)
      if (ierr /= 0 .or. .not. duplicate) call mesa_error(__FILE__, __LINE__)

      call integer_idict_create_hash(idict, ierr)
      if (ierr /= 0) call mesa_error(__FILE__, __LINE__)

      call integer_idict_lookup(idict, -51, 885, value, ierr)
      if (ierr /= 0) call mesa_error(__FILE__, __LINE__)
      if (value /= 2) call mesa_error(__FILE__, __LINE__)

      call integer_idict_lookup(idict, 1547, 974, value, ierr)
      if (ierr /= 0) call mesa_error(__FILE__, __LINE__)
      if (value /= 1) call mesa_error(__FILE__, __LINE__)

      call integer_idict_lookup(idict, 592, 8, value, ierr)
      if (ierr /= 0) call mesa_error(__FILE__, __LINE__)
      if (value /= 4) call mesa_error(__FILE__, __LINE__)

      call integer_idict_lookup(idict, 0, 18888888, value, ierr)
      if (ierr == 0) call mesa_error(__FILE__, __LINE__)
      ierr = 0

      call integer_idict_free(idict)

      write (*, *) 'okay'
      write (*, '(A)')

   end subroutine test_idict

   subroutine test_token_read
      integer :: iounit, n, i, t, ierr
      character(len=256) :: buffer, string, filename

      write (*, '(A)')
      write (*, *) 'test_token_read'
      write (*, '(A)')

      filename = 'token.txt'
      ierr = 0
      iounit = alloc_iounit(ierr)
      if (ierr /= 0) call mesa_error(__FILE__, __LINE__)
      open (unit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
      if (ierr /= 0) call mesa_error(__FILE__, __LINE__)
      n = 0
      i = 0

      do
         t = token(iounit, n, i, buffer, string)
         select case (t)
         case (string_token)
            write (*, *) 'string_token', len_trim(string), trim(string)
         case (name_token)
            write (*, *) 'name_token', len_trim(string), trim(string)
         case (left_paren_token)
            write (*, *) 'left_paren_token'
         case (right_paren_token)
            write (*, *) 'right_paren_token'
         case (comma_token)
            write (*, *) 'comma_token'
         case (eof_token)
            write (*, *) 'eof_token'
            exit
         case default
         end select

      end do

      close (iounit)
      call free_iounit(iounit)

      write (*, '(A)')
      write (*, *) 'done test_token_read'
      write (*, '(A)')

   end subroutine test_token_read

end program test_utils
