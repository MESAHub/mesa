! ***********************************************************************
!
!   Copyright (C) 2010  Rich Townsend
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

module utils_lib

  ! Uses
  
  use utils_def, only: max_io_unit
  use const_def, only: dp, qp, strlen

  use utils_nan

  ! No implicit typing

  implicit none

  ! Module variables

  logical :: assigned(max_io_unit) = .false.

  character(*), private, parameter :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
  CHARACTER(*), private, parameter :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  ! Procedures

contains

  integer function utils_OMP_GET_THREAD_NUM()
    use utils_openmp, only: eval_OMP_GET_THREAD_NUM
    utils_OMP_GET_THREAD_NUM = eval_OMP_GET_THREAD_NUM()
  end function utils_OMP_GET_THREAD_NUM

  !****

  integer function utils_OMP_GET_MAX_THREADS()
    use utils_openmp, only: eval_OMP_GET_MAX_THREADS
    utils_OMP_GET_MAX_THREADS = eval_OMP_GET_MAX_THREADS()
  end function utils_OMP_GET_MAX_THREADS

  !****

  subroutine utils_OMP_SET_NUM_THREADS(threads)
    use utils_openmp, only: eval_OMP_SET_NUM_THREADS
    integer :: threads
    call eval_OMP_SET_NUM_THREADS(threads)
  end subroutine utils_OMP_SET_NUM_THREADS

  !****

  subroutine get_compiler_version(compiler_name,compiler_version_name)
    character(len=*) :: compiler_name, compiler_version_name
    integer :: intel_compiler_build_date = 0, gcc_major = 0, gcc_minor = 0, gcc_patch = 0
    character(len=3) :: gcc_major_string, gcc_minor_string, gcc_patch_string
#ifdef __INTEL_COMPILER
    compiler_name = "ifort"
    intel_compiler_build_date = __INTEL_COMPILER_BUILD_DATE
    write(compiler_version_name,'(i8)') intel_compiler_build_date
#elif __GFORTRAN__
    compiler_name = "gfortran"
    !compiler_version_name = __VERSION__
    gcc_major = __GNUC__
    gcc_minor = __GNUC_MINOR__
    gcc_patch = __GNUC_PATCHLEVEL__
    write(gcc_major_string,'(i3)') gcc_major
    write(gcc_minor_string,'(i3)') gcc_minor
    write(gcc_patch_string,'(i3)') gcc_patch
    compiler_version_name = trim(adjustl(gcc_major_string)) // '.' &
                         // trim(adjustl(gcc_minor_string)) // '.' &
                         // trim(adjustl(gcc_patch_string))
#else
    compiler_name = "unknown"
    compiler_version_name = "unknown"
#endif
  end subroutine get_compiler_version

  !****

  subroutine get_mesasdk_version(version, ierr)
    use iso_fortran_env
    implicit none
    character(len=*), intent(out) :: version
    integer, intent(out) :: ierr
    character(len=strlen) :: mesasdk_root, filename
    integer :: unit, root_len, name_len

    ierr = 0
    version = 'unknown' !set here in case there is a problem below

    call get_environment_variable(name='MESASDK_VERSION', value=version, length=name_len, status=ierr)
    if (ierr /= 0) then
       ierr=0
       return
    endif

    call get_environment_variable(name='MESASDK_ROOT', value=mesasdk_root, length=root_len, status=ierr)
    if (ierr /= 0 .or. root_len==0) then
       ierr=0
       return
    endif

    filename=trim(mesasdk_root) // '/bin/mesasdk_version'
    open(newunit=unit, file=filename, status='old', action='read', iostat=ierr)
    if (ierr /= 0)then
       ierr=0
       return
    endif

    read(unit,'(A)')
    read(unit,'(A)')
    read(unit,'(A)')
    read(unit,'(A)')
    read(unit,'(5X, A)') version

    close(unit)

    call strip(version,'"')

  end subroutine get_mesasdk_version

  !****

  elemental subroutine strip(string,set)
    character(len=*), intent(inout) :: string
    character(len=*), intent(in)    :: set
    integer :: old, new, stride
    old=1; new=1
    do
       stride=scan(string(old:),set)
       if(stride>0)then
          string(new:new+stride-2)=string(old:old+stride-2)
          old=old+stride
          new=new+stride-1
       else
          string(new:)=string(old:)
          return
       end if
    end do
  end subroutine strip

  !****

  subroutine realloc_double(ptr,new_size,ierr)
    real(dp), pointer :: ptr(:)
    integer, intent(in) :: new_size
    integer, intent(out) :: ierr
    real(dp), pointer :: new_ptr(:)
    integer :: i
    ierr = 0
    allocate(new_ptr(new_size),stat=ierr)
    if (ierr /= 0) return
    if (associated(ptr)) then
       do i = 1, min(new_size,size(ptr,dim=1))
          new_ptr(i) = ptr(i)
       end do
       deallocate(ptr)
    end if
    ptr => new_ptr
  end subroutine realloc_double

  !****
  
  subroutine realloc_double2(ptr,new_size1,new_size2,ierr)
    real(dp), pointer :: ptr(:,:)
    integer, intent(in) :: new_size1,new_size2
    integer, intent(out) :: ierr
    real(dp), pointer :: new_ptr(:,:)
    integer :: i1,i2, i,j
    ierr = 0
    allocate(new_ptr(new_size1,new_size2),stat=ierr)
    if (ierr /= 0) return
    if (associated(ptr)) then
       i1 = min(new_size1,size(ptr,dim=1))
       i2 = min(new_size2,size(ptr,dim=2))
       ! ifort uses stack for array copy temp storage
       ! for large copies, this can produce seg faults
       ! doing the explicit loops seems to be safe
       !new_ptr(1:i1,1:i2) = ptr(1:i1,1:i2)
       do i=1,i1
          do j=1,i2
             new_ptr(i,j) = ptr(i,j)
          end do
       end do
       deallocate(ptr)
    end if
    ptr => new_ptr
  end subroutine realloc_double2

  !****
  
  subroutine realloc_quad(ptr,new_size,ierr)
    real(qp), pointer :: ptr(:)
    integer, intent(in) :: new_size
    integer, intent(out) :: ierr
    real(qp), pointer :: new_ptr(:)
    integer :: i
    ierr = 0
    allocate(new_ptr(new_size),stat=ierr)
    if (ierr /= 0) return
    if (associated(ptr)) then
       do i = 1, min(new_size,size(ptr,dim=1))
          new_ptr(i) = ptr(i)
       end do
       deallocate(ptr)
    end if
    ptr => new_ptr
  end subroutine realloc_quad

  !****

  subroutine realloc_quad2(ptr,new_size1,new_size2,ierr)
    real(qp), pointer :: ptr(:,:)
    integer, intent(in) :: new_size1,new_size2
    integer, intent(out) :: ierr
    real(qp), pointer :: new_ptr(:,:)
    integer :: i1,i2, i,j
    ierr = 0
    allocate(new_ptr(new_size1,new_size2),stat=ierr)
    if (ierr /= 0) return
    if (associated(ptr)) then
       i1 = min(new_size1,size(ptr,dim=1))
       i2 = min(new_size2,size(ptr,dim=2))
       ! ifort uses stack for array copy temp storage
       ! for large copies, this can produce seg faults
       ! doing the explicit loops seems to be safe
       !new_ptr(1:i1,1:i2) = ptr(1:i1,1:i2)
       do i=1,i1
          do j=1,i2
             new_ptr(i,j) = ptr(i,j)
          end do
       end do
       deallocate(ptr)
    end if
    ptr => new_ptr
  end subroutine realloc_quad2

  !****
  
  subroutine realloc_double3(ptr,new_size1,new_size2,new_size3,ierr)
    real(dp), pointer :: ptr(:,:,:)
    integer, intent(in) :: new_size1,new_size2,new_size3
    integer, intent(out) :: ierr
    real(dp), pointer :: new_ptr(:,:,:)
    integer :: i1,i2,i3, i,j,k
    ierr = 0
    allocate(new_ptr(new_size1,new_size2,new_size3),stat=ierr)
    if (ierr /= 0) return
    if (associated(ptr)) then
       i1 = min(new_size1,size(ptr,dim=1))
       i2 = min(new_size2,size(ptr,dim=2))
       i3 = min(new_size3,size(ptr,dim=3))
       ! ifort uses stack for array copy temp storage
       ! for large copies, this can produce seg faults
       ! doing the explicit loops seems to be safe
       !new_ptr(1:i1,1:i2,1:i3) = ptr(1:i1,1:i2,1:i3)
       do i=1,i1
          do j=1,i2
             do k=1,i3
                new_ptr(i,j,k) = ptr(i,j,k)
             end do
          end do
       end do
       deallocate(ptr)
    end if
    ptr => new_ptr
  end subroutine realloc_double3

  !****

  subroutine realloc_real(ptr,new_size,ierr)
    real, pointer :: ptr(:)
    integer, intent(in) :: new_size
    integer, intent(out) :: ierr
    real, pointer :: new_ptr(:)
    integer :: i
    ierr = 0
    allocate(new_ptr(new_size),stat=ierr)
    if (ierr /= 0) return
    if (associated(ptr)) then
       do i=1,min(new_size,size(ptr,dim=1))
          new_ptr(i) = ptr(i)
       end do
       deallocate(ptr)
    end if
    ptr => new_ptr
  end subroutine realloc_real

  !****

  subroutine realloc_integer(ptr,new_size,ierr)
    integer, pointer :: ptr(:)
    integer, intent(in) :: new_size
    integer, intent(out) :: ierr
    integer, pointer :: new_ptr(:)
    integer :: i
    ierr = 0
    allocate(new_ptr(new_size),stat=ierr)
    if (ierr /= 0) return
    if (associated(ptr)) then
       do i = 1, min(new_size,size(ptr,dim=1))
          new_ptr(i) = ptr(i)
       end do
       deallocate(ptr)
    end if
    ptr => new_ptr
  end subroutine realloc_integer

  !****

  subroutine realloc_integer2(ptr,new_size1,new_size2,ierr)
    integer, pointer :: ptr(:,:)
    integer, intent(in) :: new_size1,new_size2
    integer, intent(out) :: ierr
    integer, pointer :: new_ptr(:,:)
    integer :: i1,i2, i,j
    ierr = 0
    allocate(new_ptr(new_size1,new_size2),stat=ierr)
    if (ierr /= 0) return
    if (associated(ptr)) then
       i1 = min(new_size1,size(ptr,dim=1))
       i2 = min(new_size2,size(ptr,dim=2))
       ! ifort uses stack for array copy temp storage
       ! for large copies, this can produce seg faults
       ! doing the explicit loops seems to be safe
       !new_ptr(1:i1,1:i2) = ptr(1:i1,1:i2)
       do i=1,i1
          do j=1,i2
             new_ptr(i,j) = ptr(i,j)
          end do
       end do
       deallocate(ptr)
    end if
    ptr => new_ptr
  end subroutine realloc_integer2

  !****

  subroutine realloc_logical(ptr,new_size,ierr)
    logical, pointer :: ptr(:)
    integer, intent(in) :: new_size
    integer, intent(out) :: ierr
    logical, pointer :: new_ptr(:)
    integer :: i
    ierr = 0
    allocate(new_ptr(new_size),stat=ierr)
    if (ierr /= 0) return
    if (associated(ptr)) then
       do i = 1, min(new_size,size(ptr,dim=1))
          new_ptr(i) = ptr(i)
       end do
       deallocate(ptr)
    end if
    ptr => new_ptr
  end subroutine realloc_logical

  !****

  subroutine do1D(ptr,sz,dealloc,ierr)
    real(dp),dimension(:),pointer::ptr
    integer, intent(in) :: sz
    logical, intent(in) :: dealloc
    integer, intent(out) :: ierr
    if (dealloc) then
       deallocate(ptr,stat=ierr)
    else
       allocate(ptr(sz),stat=ierr)
    end if
  end subroutine do1D

  !****

  subroutine do2D(ptr,sz1,sz2,dealloc,ierr)
    real(dp),dimension(:,:),pointer::ptr
    integer, intent(in) :: sz1,sz2
    logical, intent(in) :: dealloc
    integer, intent(out) :: ierr
    if (dealloc) then
       deallocate(ptr,stat=ierr)
    else
       allocate(ptr(sz1,sz2),stat=ierr)
    end if
  end subroutine do2D

  !****

  subroutine do3D(ptr,sz1,sz2,sz3,dealloc,ierr)
    real(dp),dimension(:,:,:),pointer::ptr
    integer, intent(in) :: sz1,sz2,sz3
    logical, intent(in) :: dealloc
    integer, intent(out) :: ierr
    if (dealloc) then
       deallocate(ptr,stat=ierr)
    else
       allocate(ptr(sz1,sz2,sz3),stat=ierr)
    end if
  end subroutine do3D

  !****

  subroutine do4D(ptr,sz1,sz2,sz3,sz4,dealloc,ierr)
    real(dp),dimension(:,:,:,:),pointer::ptr
    integer, intent(in) :: sz1,sz2,sz3,sz4
    logical, intent(in) :: dealloc
    integer, intent(out) :: ierr
    if (dealloc) then
       deallocate(ptr,stat=ierr)
    else
       allocate(ptr(sz1,sz2,sz3,sz4),stat=ierr)
    end if
  end subroutine do4D

  !****

  subroutine do1D_integer(ptr,sz,dealloc,ierr)
    integer,dimension(:),pointer::ptr
    integer, intent(in) :: sz
    logical, intent(in) :: dealloc
    integer, intent(out) :: ierr
    if (dealloc) then
       deallocate(ptr,stat=ierr)
    else
       allocate(ptr(sz),stat=ierr)
    end if
  end subroutine do1D_integer

  !***

  subroutine do2D_integer(ptr,sz1,sz2,dealloc,ierr)
    integer,dimension(:,:),pointer::ptr
    integer, intent(in) :: sz1,sz2
    logical, intent(in) :: dealloc
    integer, intent(out) :: ierr
    if (dealloc) then
       deallocate(ptr,stat=ierr)
    else
       allocate(ptr(sz1,sz2),stat=ierr)
    end if
  end subroutine do2D_integer

  !****

  subroutine do1D_logical(ptr,sz,dealloc,ierr)
    logical,dimension(:),pointer::ptr
    integer, intent(in) :: sz
    logical, intent(in) :: dealloc
    integer, intent(out) :: ierr
    if (dealloc) then
       deallocate(ptr,stat=ierr)
    else
       allocate(ptr(sz),stat=ierr)
    end if
  end subroutine do1D_logical

  !****

  subroutine alloc1(sz,a,ierr)
    real(dp), dimension(:), pointer :: a
    integer, intent(in) :: sz
    integer, intent(out) :: ierr
    allocate(a(sz),stat=ierr); if (ierr /= 0) return
  end subroutine alloc1

  !****
  
  subroutine alloc2(sz1,sz2,a,ierr)
    real(dp), dimension(:,:), pointer :: a
    integer, intent(in) :: sz1,sz2
    integer, intent(out) :: ierr
    allocate(a(sz1,sz2),stat=ierr); if (ierr /= 0) return
  end subroutine alloc2

  !****

  subroutine alloc3(sz1,sz2,sz3,a,ierr)
    real(dp), dimension(:,:,:), pointer :: a
    integer, intent(in) :: sz1,sz2,sz3
    integer, intent(out) :: ierr
    allocate(a(sz1,sz2,sz3),stat=ierr); if (ierr /= 0) return
  end subroutine alloc3

  !****

  subroutine realloc_if_needed_1(ptr,sz,extra,ierr)
    real(dp), pointer :: ptr(:)
    integer, intent(in) :: sz, extra
    integer, intent(out) :: ierr
    ierr = 0
    if (associated(ptr)) then
       if (size(ptr,dim=1) >= sz) return
    end if
    call realloc_double(ptr, sz + extra, ierr)
  end subroutine realloc_if_needed_1

  !****

  subroutine quad_realloc_if_needed_1(ptr,sz,extra,ierr)
    real(qp), pointer :: ptr(:)
    integer, intent(in) :: sz, extra
    integer, intent(out) :: ierr
    ierr = 0
    if (associated(ptr)) then
       if (size(ptr,dim=1) >= sz) return
    end if
    call realloc_quad(ptr, sz + extra, ierr)
  end subroutine quad_realloc_if_needed_1

  !****

  subroutine realloc_integer_if_needed_1(ptr,sz,extra,ierr)
    integer, pointer :: ptr(:)
    integer, intent(in) :: sz, extra
    integer, intent(out) :: ierr
    ierr = 0
    if (associated(ptr)) then
       if (size(ptr,dim=1) >= sz) return
    end if
    call realloc_integer(ptr, sz + extra, ierr)
  end subroutine realloc_integer_if_needed_1

  !****

  subroutine enlarge_if_needed_1(ptr,sz,extra,ierr)
    real(dp), pointer :: ptr(:)
    integer, intent(in) :: sz, extra
    integer, intent(out) :: ierr
    ierr = 0
    if (associated(ptr)) then
       if (size(ptr,dim=1) >= sz) return
       deallocate(ptr)
    end if
    allocate(ptr(sz + extra), stat=ierr)
  end subroutine enlarge_if_needed_1

  !****

  subroutine enlarge_if_needed_2(ptr,sz1,sz2,extra,ierr)
    real(dp), pointer :: ptr(:,:)
    integer, intent(in) :: sz1, sz2, extra
    integer, intent(out) :: ierr
    ierr = 0
    if (associated(ptr)) then
       if (size(ptr,1) == sz1 .and. size(ptr,2) >= sz2) return
       deallocate(ptr)
    end if
    allocate(ptr(sz1, sz2 + extra), stat=ierr)
  end subroutine enlarge_if_needed_2

  !****

  subroutine quad_enlarge_if_needed_1(ptr,sz,extra,ierr)
    real(qp), pointer :: ptr(:)
    integer, intent(in) :: sz, extra
    integer, intent(out) :: ierr
    ierr = 0
    if (associated(ptr)) then
       if (size(ptr,dim=1) >= sz) return
       deallocate(ptr)
    end if
    allocate(ptr(sz + extra), stat=ierr)
  end subroutine quad_enlarge_if_needed_1

  !****

  subroutine enlarge_integer_if_needed_1(ptr,sz,extra,ierr)
    integer, pointer :: ptr(:)
    integer, intent(in) :: sz, extra
    integer, intent(out) :: ierr
    ierr = 0
    if (associated(ptr)) then
       if (size(ptr,dim=1) >= sz) return
       deallocate(ptr)
    end if
    allocate(ptr(sz + extra), stat=ierr)
  end subroutine enlarge_integer_if_needed_1

  !****

  subroutine remove_underbars(str, name)
    character (len=*), intent(in) :: str
    character (len=*), intent(out) :: name
    integer :: i, len
    len = len_trim(str)
    name = ''
    do i=1,len
       if (str(i:i) == '_') then
          name(i:i) = ' '
       else
          name(i:i) = str(i:i)
       end if
    end do
  end subroutine remove_underbars

  !****

  integer function token(iounit, n, i, buffer, string)
    use utils_def
    integer, intent(in) :: iounit
    integer, intent(inout) :: n ! number of characters currently in buffer
    integer, intent(inout) :: i ! number of characters already read from buffer
    character (len=*), intent(inout) :: buffer ! line of text from input file
    character (len=*), intent(inout) :: string ! holds string or name for string or name token
    character (len=1) :: tab_str

    integer :: info, j, j1, j2, l, str_len

    token = 0
    info = 0

    line_loop: do
       do while (i >= n)
          read(iounit,fmt='(a)',iostat=info) buffer
          if (info /= 0) then
             token = eof_token
             return 
          end if
          n = len_trim(buffer)
          i = 0
          !write(*,'(i6,3x,a)') n, trim(buffer)
       end do
       token_loop: do while (i < n) ! have non-empty buffer
          i = i+1
          if (buffer(i:i) == char(9)) cycle token_loop ! skip tabs
          select case(buffer(i:i))
          case ('!')
             i = n
             cycle line_loop
          case (' ')
             cycle token_loop
          case ('&') ! ignore &'s
             cycle token_loop
          case ('(')
             token = left_paren_token; return
          case (')')
             token = right_paren_token; return
          case (',')
             token = comma_token; return
          case ('"')
             j = 1; str_len = len(string)
             do
                i = i+1
                if (i > n) exit
                if (buffer(i:i) == '"') exit
                if (j > str_len) exit
                string(j:j) = buffer(i:i)
                j = j+1
             end do
             do while (j <= str_len)
                string(j:j) = ' '
                j = j+1
             end do
             token = string_token
             return
          case ('''')
             j = 1; str_len = len(string)
             do
                i = i+1
                if (i > n) exit
                if (buffer(i:i) == '''') exit
                if (j > str_len) exit
                string(j:j) = buffer(i:i)
                j = j+1
             end do
             do while (j <= str_len)
                string(j:j) = ' '
                j = j+1
             end do
             token = string_token
             return
          case default
             j1 = i; j2 = i
             name_loop: do
                if (i+1 > n) exit
                if (buffer(i+1:i+1) == ' ') exit
                if (buffer(i+1:i+1) == '(') exit
                if (buffer(i+1:i+1) == ')') exit
                if (buffer(i+1:i+1) == ',') exit
                i = i+1
                j2 = i
             end do name_loop
             str_len = len(string)
             l = j2-j1+1
             if (l > str_len) then
                l = str_len
                j2 = l+j1-1
             end if
             string(1:l) = buffer(j1:j2)
             do j = l+1, str_len
                string(j:j) = ' '
             end do
             token = name_token
             return
          end select
       end do token_loop
    end do line_loop

  end function token

  !****

  subroutine integer_dict_define_and_report_duplicates(dict, key, value, duplicate, ierr)
    use utils_dict
    type (integer_dict), pointer :: dict ! pass null for empty dict
    character (len=*), intent(in) :: key
    integer, intent(in) :: value
    logical, intent(out) :: duplicate ! true if key was already defined
    ! if already defined, old value is replaced by new one.
    integer, intent(out) :: ierr ! error if len_trim(key) > maxlen_key_string
    call do_integer_dict_define(dict, key, value, duplicate, ierr)
  end subroutine integer_dict_define_and_report_duplicates

  !****

  subroutine integer_dict_define(dict, key, value, ierr)
    use utils_def, only: integer_dict
    type (integer_dict), pointer :: dict ! pass null for empty dict
    character (len=*), intent(in) :: key
    integer, intent(in) :: value
    integer, intent(out) :: ierr ! error if len_trim(key) > maxlen_key_string
    logical :: duplicate
    call integer_dict_define_and_report_duplicates(dict, key, value, duplicate, ierr)
  end subroutine integer_dict_define

  !****

  subroutine integer_dict_create_hash(dict, ierr)
    use utils_dict
    type (integer_dict), pointer :: dict
    integer, intent(out) :: ierr
    call do_integer_dict_create_hash(dict, ierr)
  end subroutine integer_dict_create_hash

  !****

  subroutine integer_dict_lookup(dict, key, value, ierr)
    use utils_dict
    type (integer_dict), pointer :: dict
    character (len=*), intent(in) :: key
    integer, intent(out) :: value
    integer, intent(out) :: ierr ! 0 if found key in dict, -1 if didn't
    call do_integer_dict_lookup(dict, key, value, ierr)
  end subroutine integer_dict_lookup

  !****

  integer function integer_dict_size(dict) ! number of entries
    use utils_dict
    type (integer_dict), pointer :: dict
    integer_dict_size = size_integer_dict(dict)
  end function integer_dict_size

  !****

  subroutine integer_dict_map(dict, fcn)
    use utils_dict
    type (integer_dict), pointer :: dict
    interface
       subroutine fcn(key, value, ierr)
         character (len=*), intent(in) :: key
         integer, intent(in) :: value
         integer, intent(out) :: ierr ! /= 0 means terminate map calls
       end subroutine fcn
    end interface
    integer :: ierr
    call do_integer_dict_map(dict, fcn, ierr)
  end subroutine integer_dict_map

  !****

  subroutine get_dict_entries(dict, keys, values)
    use utils_dict
    type (integer_dict), pointer :: dict
    character (len=maxlen_key_string), pointer :: keys(:)
    integer, pointer :: values(:)
    call do_get_dict_entries(dict, keys, values)
  end subroutine get_dict_entries

  !****

  subroutine integer_dict_free(dict)
    use utils_dict
    type (integer_dict), pointer :: dict
    call do_integer_dict_free(dict)
  end subroutine integer_dict_free

  !****

  subroutine integer_idict_define_and_report_duplicates(idict, key1, key2, value, duplicate, ierr)
    use utils_idict
    type (integer_idict), pointer :: idict ! pass null for empty idict
    integer, intent(in) :: key1, key2, value
    logical, intent(out) :: duplicate ! true if key was already defined
    ! if already defined, old value is replaced by new one.
    integer, intent(out) :: ierr
    call do_integer_idict_define(idict, key1, key2, value, duplicate, ierr)
  end subroutine integer_idict_define_and_report_duplicates

  !****

  subroutine integer_idict_define(idict, key1, key2, value, ierr)
    use utils_def, only: integer_idict
    type (integer_idict), pointer :: idict ! pass null for empty idict
    integer, intent(in) :: key1, key2, value
    integer, intent(out) :: ierr
    logical :: duplicate
    call integer_idict_define_and_report_duplicates( &
         idict, key1, key2, value, duplicate, ierr)
  end subroutine integer_idict_define

  !****

  subroutine integer_idict_create_hash(idict, ierr)
    use utils_idict
    type (integer_idict), pointer :: idict
    integer, intent(out) :: ierr
    call do_integer_idict_create_hash(idict, ierr)
  end subroutine integer_idict_create_hash

  !****

  subroutine integer_idict_lookup(idict, key1, key2, value, ierr)
    use utils_idict
    type (integer_idict), pointer :: idict
    integer, intent(in) :: key1, key2
    integer, intent(out) :: value
    integer, intent(out) :: ierr ! 0 if found key in idict, -1 if didn't
    call do_integer_idict_lookup(idict, key1, key2, value, ierr)
  end subroutine integer_idict_lookup

  !****

  integer function integer_idict_size(idict) ! number of entries
    use utils_idict
    type (integer_idict), pointer :: idict
    integer_idict_size = size_integer_idict(idict)
  end function integer_idict_size

  !****

  subroutine integer_idict_map(idict, fcn)
    use utils_idict
    type (integer_idict), pointer :: idict
    interface
       subroutine fcn(key1, key2, value, ierr)
         integer, intent(in) :: key1, key2, value
         integer, intent(out) :: ierr ! /= 0 means terminate map calls
       end subroutine fcn
    end interface
    integer :: ierr
    call do_integer_idict_map(idict, fcn, ierr)
  end subroutine integer_idict_map

  !****

  subroutine get_idict_entries(idict, key1s, key2s, values)
    use utils_idict
    type (integer_idict), pointer :: idict
    integer, pointer, dimension(:) :: key1s, key2s, values
    call do_get_idict_entries(idict, key1s, key2s, values)
  end subroutine get_idict_entries

  !****
  
  subroutine integer_idict_free(idict)
    use utils_idict
    type (integer_idict), pointer :: idict
    call do_integer_idict_free(idict)
  end subroutine integer_idict_free

  !****

  function StrUpCase ( Input_String ) result ( Output_String )
    character(len=*), intent(in) :: Input_String
    character(len(Input_string)) :: Output_String
    integer :: i, n
    Output_String = Input_String
    do i = 1, len( Output_String )
       n = index( LOWER_CASE, Output_String( i:i ) )
       if ( n /= 0 ) Output_String( i:i ) = UPPER_CASE( n:n )
    end do
  end function StrUpCase

  !****

  function StrLowCase ( Input_String ) result ( Output_String )
    character(len=*), intent(in) :: Input_String
    character(len(Input_string)) :: Output_String
    integer :: i, n
    Output_String = Input_String
    do i = 1, len( Output_String )
       n = index( UPPER_CASE, Output_String( i:i ) )
       if ( n /= 0 ) Output_String( i:i ) = LOWER_CASE( n:n )
    end do
  end function StrLowCase

  !****

  subroutine mkdir(folder)
    use utils_system, only : mkdir_p
    character(len=*), intent(in) :: folder
    integer :: res
    
    res = mkdir_p(folder)
    
    if(res/=0)then
      write(*,*) "mkdir failed for ",trim(folder)
      write(*,*) "error code ",res
      call mesa_error(__FILE__,__LINE__)
    end if
    
  end subroutine mkdir
  
  subroutine mv(file_in,file_out,skip_errors)
    use utils_system, only: mv_c => mv
    character(len=*),intent(in) :: file_in,file_out
    logical, optional, intent(in) :: skip_errors
    integer res
    
    res = mv_c(file_in,file_out)
    
    if(res/=0)then
       if (present(skip_errors))then
          if (skip_errors) then
             write(*,*) "mv failed for '"//trim(file_in)//"' '"//trim(file_out)//"' skipping"
          else
            call error()
          end if
       else
          call error()
       end if
    end if
    
    contains 
    
      subroutine error()
            write(*,*) "mv failed for '"//trim(file_in)//"' '"//trim(file_out)//"'"
            write(*,*) "Error code: ",res
            call mesa_error(__FILE__,__LINE__)
      end subroutine error
  
  end subroutine mv
  
  subroutine cp_file(file_in,file_out,skip_errors)
    use utils_system, only: cp_c => cp
    character(len=*),intent(in) :: file_in,file_out
    logical, optional, intent(in) :: skip_errors
    integer res
    
    res = cp_c(file_in,file_out)
    
    if(res/=0)then
       if (present(skip_errors))then
          if (skip_errors) then
             write(*,*) "cp failed for '"//trim(file_in)//"' '"//trim(file_out)//"' skipping"
          else
            call error()
          end if
       else
          call error()
       end if
    end if
    
    contains 
    
      subroutine error()
            write(*,*) "cp failed for '"//trim(file_in)//"' '"//trim(file_out)//"'"
            write(*,*) "Error code: ",res
            call mesa_error(__FILE__,__LINE__)
      end subroutine error
  
  end subroutine cp_file

   logical function folder_exists(folder)
      use utils_system, only: is_dir
      character(len=*),intent(in) :: folder

      folder_exists = is_dir(folder)

   end function folder_exists

  integer function alloc_iounit(ierr)
    use utils_def
    integer, intent(out) :: ierr
    integer :: i
    !$omp critical (utils_alloc_io_unit)
    ierr = 0
    alloc_iounit = -1
    do i = min_io_unit, max_io_unit
       if (.not. assigned(i)) then
          assigned(i) = .true.
          alloc_iounit = i
          exit
       end if
    end do
    !$omp end critical (utils_alloc_io_unit)
    if (alloc_iounit == -1) then
       ierr = -1
    end if
  end function alloc_iounit

  !****

  integer function number_iounits_allocated()
    use utils_def
    integer :: i, cnt
    cnt = 0
    !$omp critical (utils_alloc_io_unit)
    do i = min_io_unit, max_io_unit
       if (assigned(i)) cnt = cnt + 1
    end do
    !$omp end critical (utils_alloc_io_unit)
    number_iounits_allocated = cnt
  end function number_iounits_allocated

  !****

  subroutine free_iounit(iounit)
    use utils_def
    integer, intent(in) :: iounit
    logical :: bad_iounit
    bad_iounit = .false.
    !$omp critical (utils_alloc_io_unit)
    if (iounit >= min_io_unit .and. iounit <= max_io_unit) then
       assigned(iounit) = .false.
    else
       bad_iounit = .true.
    end if
    !$omp end critical (utils_alloc_io_unit)
    if (bad_iounit) then
       write(*,*) 'called free_iounit with invalid arg', iounit
       stop 'free_iounit'
    end if
  end subroutine free_iounit

  !****

  subroutine append_line(n, arry, filename, format_str, initialize, ierr)
    integer, intent(in) :: n
    real(dp), intent(in) :: arry(n)
    character (len=256), intent(in) :: filename, format_str
    logical, intent(in) :: initialize
    integer, intent(out) :: ierr
    integer :: iounit
    ierr = 0
    iounit = alloc_iounit(ierr); if (ierr /= 0) return
    if (initialize) then
       open(unit=iounit, file=trim(filename), action='write', iostat=ierr)
    else
       open(unit=iounit, file=trim(filename), position='append', action='write', iostat=ierr)
    end if
    if (ierr /= 0) then
       call free_iounit(iounit)
       return
    end if
    write(unit=iounit, fmt=format_str) arry
    close(iounit)
    call free_iounit(iounit)
  end subroutine append_line

  !****

  subroutine append_data(n, arry, filename, format_str, initialize, ierr)
    integer, intent(in) :: n
    real(dp), intent(in) :: arry(n)
    character (len=256), intent(in) :: filename, format_str
    logical, intent(in) :: initialize
    integer, intent(out) :: ierr
    integer :: iounit, i
    ierr = 0
    iounit = alloc_iounit(ierr); if (ierr /= 0) return
    if (initialize) then
       open(unit=iounit, file=trim(filename), action='write', iostat=ierr)
    else
       open(unit=iounit, file=trim(filename), position='append', action='write', iostat=ierr)
    end if
    if (ierr /= 0) then
       call free_iounit(iounit)
       return
    end if
    do i=1,n
       write(unit=iounit, fmt=format_str) arry(i)
    end do
    close(iounit)
    call free_iounit(iounit)
  end subroutine append_data

  subroutine append_data_filter_badnums(n, arry, filename, format_str, initialize, ierr)
    integer, intent(in) :: n
    real(dp), intent(in) :: arry(n)
    character (len=256), intent(in) :: filename, format_str
    logical, intent(in) :: initialize
    integer, intent(out) :: ierr
    integer :: iounit, i
    ierr = 0
    iounit = alloc_iounit(ierr); if (ierr /= 0) return
    if (initialize) then
       open(unit=iounit, file=trim(filename), action='write', iostat=ierr)
    else
       open(unit=iounit, file=trim(filename), position='append', action='write', iostat=ierr)
    end if
    if (ierr /= 0) then
       call free_iounit(iounit)
       return
    end if
    do i=1,n
       if (is_bad(arry(i))) then
        write(unit=iounit, fmt=format_str) -1d99
       else
        write(unit=iounit, fmt=format_str) arry(i)
       end if
    end do
    close(iounit)
    call free_iounit(iounit)
  end subroutine append_data_filter_badnums


   subroutine mesa_error(file, line,msg)
      use iso_fortran_env, only: error_unit
      character(len=*), intent(in) :: file
      character(len=*), optional,intent(in) :: msg
      integer, intent(in) :: line
      character (len=strlen) :: bt_disable
      if(present(msg)) then
         write(error_unit,"(A, A, A, I4, A, A)") "File: ", file, ", Line: ",line,", Message: ",msg
      else
         write(error_unit,"(A, A, A, I4)") "File: ", file, ", Line: ", line
      end if
      call get_environment_variable("MESA_ERROR_BACKTRACE_DISABLE", bt_disable)
      if (len_trim(bt_disable) > 0) then
         stop 1
      else
         error stop 1
      end if
   end subroutine mesa_error


    ! If flag is true return str1 else return str2
    character(len=strlen) function switch_str(str1,str2,flag)
        character(len=*), intent(in) :: str1,str2
        logical, intent(in) :: flag
        logical, parameter :: dbg=.false.
        
        if(flag) then
            switch_str=str1(1:min(len_trim(str1),strlen))
            if(len_trim(str1) > strlen .and. dbg) & 
                write(*,*) "Warning ",trim(str1), "truncated to ",switch_str
        else
            switch_str=str2(1:min(len_trim(str2),strlen))
            if(len_trim(str2) > strlen .and. dbg) &
                write(*,*) "Warning ",trim(str2), "truncated to ",switch_str
        end if
    
    end function switch_str

   subroutine split_line(line, num, out)
      ! Given a string line, split on whitespace, into num sub-strings storing them in out
      character(len=*),intent(in) :: line
      integer, intent(in) :: num
      character(len=*),dimension(:) :: out

      integer :: i,kstart,k

      if(size(out)<num) call mesa_error(__FILE__,__LINE__,'out array not large enough for num sub-strings')

      out = ''

      k = 1
      kstart = 1
      outer: do i=1, num
         inner: do
            if(line(k:k)==' ' .and. line(k+1:k+1)==' ' .and. k < len(line)) then
               k = k+1
               cycle inner
            end if

            if(line(k:k)==' ' .or. k > len(line))then
               !write(*,*) '*',i,kstart,k,line(kstart:k-1)
               out(i) = line(kstart:k-1)
               k = k+1
               kstart = k
               cycle outer
            end if
            k = k + 1
         end do inner
      end do outer

      end subroutine split_line

      
   ! backward compatibility so Bill can debug older versions of files without changing these calls
      logical function is_bad_num(x)
         real(dp), intent(in) :: x
         is_bad_num = is_bad(x)
      end function is_bad_num
      
            
      logical function is_bad_real(x)
         real, intent(in) :: x
         is_bad_real = is_bad(x)
      end function is_bad_real
      
            
      logical function is_bad_quad(x)
         real(qp), intent(in) :: x
         is_bad_quad = is_bad(x)
      end function is_bad_quad
            
      subroutine fill_with_NaNs(ptr)
         real(dp) :: ptr(:)
         call set_nan(ptr)
      end subroutine fill_with_NaNs
      
      
      subroutine fill_with_NaNs_2D(ptr)
         real(dp) :: ptr(:,:)
         call set_nan(ptr)
      end subroutine fill_with_NaNs_2D
      
      
      subroutine fill_with_NaNs_3D(ptr)
         real(dp) :: ptr(:,:,:)
         call set_nan(ptr)
      end subroutine fill_with_NaNs_3D
      
      
      subroutine fill_with_NaNs_4D(ptr)
         real(dp) :: ptr(:,:,:,:)
         call set_nan(ptr)
      end subroutine fill_with_NaNs_4D
      
      subroutine set_to_NaN(x)
         real(dp) :: x
         real(dp) :: xa(1)
         call set_nan(xa)
         x = xa(1)
      end subroutine set_to_NaN
      

end module utils_lib

