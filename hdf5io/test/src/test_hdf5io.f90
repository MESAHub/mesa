program test_hdf5io
   
   ! Uses
   
   use const_lib
   use hdf5io_lib
   
   ! No implicit typing
   
   implicit none
   
   integer, parameter :: n = 3
   
   character(12) :: s_0
   integer(i4) :: is_0
   integer(i4) :: id_0
   real(sp) :: rs_0
   real(dp) :: rd_0
   real(sp) :: cs_0
   real(dp) :: cd_0
   logical :: l_0
   character(12) :: s_1(n)
   integer(i4) :: is_1(n)
   integer(i4) :: id_1(n)
   real(sp) :: rs_1(n)
   real(dp) :: rd_1(n)
   real(sp) :: cs_1(n)
   real(dp) :: cd_1(n)
   logical :: l_1(n)
   
   ! Run tests
   
   call write_file()
   call read_file()

contains
   
   subroutine write_file()
      
      type(hdf5io_t) :: hi
      
      ! Set up values
      
      s_0 = 'hello world!'
      s_1 = s_0
      
      is_0 = 1_i4
      is_1 = is_0
      
      id_0 = 2_i8
      id_1 = id_0
      
      rs_0 = 3._sp
      rs_1 = rs_0
      
      rd_0 = 4._dp
      rd_1 = rd_0
      
      cs_0 = CMPLX(rs_0, 5 * rs_0, sp)
      cs_1 = cs_0
      
      cd_0 = CMPLX(rd_0, 6 * rd_0, dp)
      cd_1 = cd_0
      
      l_0 = .TRUE.
      l_1 = [.FALSE., .TRUE., .FALSE.]
      
      ! Read the file
      
      hi = hdf5io_t('test.h5', CREATE_FILE)
      
      call hi%write_dset('s_0', s_0)
      call hi%write_dset('s_1', s_1)
      
      call hi%write_dset('is_0', is_0)
      call hi%write_dset('is_1', is_1)
      
      call hi%write_dset('id_0', id_0)
      call hi%write_dset('id_1', id_1)
      
      call hi%write_dset('rs_0', rs_0)
      call hi%write_dset('rs_1', rs_1)
      
      call hi%write_dset('rd_0', rd_0)
      call hi%write_dset('rd_1', rd_1)
      
      call hi%write_dset('cs_0', cs_0)
      call hi%write_dset('cs_1', cs_1)
      
      call hi%write_dset('cd_0', cd_0)
      call hi%write_dset('cd_1', cd_1)
      
      call hi%write_dset('l_0', l_0)
      call hi%write_dset('l_1', l_1)
      
      call hi%final()
      
      ! Finish
   
   end subroutine write_file
   
   subroutine read_file()
      
      type(hdf5io_t) :: hi
      
      character(12) :: s_0_chk
      integer(i4) :: is_0_chk
      integer(i4) :: id_0_chk
      real(sp) :: rs_0_chk
      real(dp) :: rd_0_chk
      real(sp) :: cs_0_chk
      real(dp) :: cd_0_chk
      logical :: l_0_chk
      character(12) :: s_1_chk(n)
      integer(i4) :: is_1_chk(n)
      integer(i4) :: id_1_chk(n)
      real(sp) :: rs_1_chk(n)
      real(dp) :: rd_1_chk(n)
      real(sp) :: cs_1_chk(n)
      real(dp) :: cd_1_chk(n)
      logical :: l_1_chk(n)
      
      100   format(A)
      110   format(2X, 'scalar', 1X, 1L)
      120   format(2X, 'array ', 1X, 1L)
      
      ! Open the file
      
      hi = hdf5io_t('test.h5', OPEN_FILE_RO)
      
      ! Read and check values
      
      call hi%read_dset('s_0', s_0_chk)
      call hi%read_dset('s_1', s_1_chk)
      
      write (*, 100) 'character dataset'
      write (*, 110) s_0 == s_0_chk
      write (*, 120) ALL(s_1 == s_1_chk)
      
      !
      
      call hi%read_dset('is_0', is_0_chk)
      call hi%read_dset('is_1', is_1_chk)
      
      write (*, 100) 'integer(i4) dataset'
      write (*, 110) is_0 == is_0_chk
      write (*, 120) ALL(is_1 == is_1_chk)
      
      !
      
      call hi%read_dset('id_0', id_0_chk)
      call hi%read_dset('id_1', id_1_chk)
      
      write (*, 100) 'integer(i8) dataset'
      write (*, 110) id_0 == id_0_chk
      write (*, 120) ALL(id_1 == id_1_chk)
      
      !
      
      call hi%read_dset('rs_0', rs_0_chk)
      call hi%read_dset('rs_1', rs_1_chk)
      
      write (*, 100) 'real(sp) dataset'
      write (*, 110) rs_0 == rs_0_chk
      write (*, 120) ALL(rs_1 == rs_1_chk)
      
      !
      
      call hi%read_dset('rd_0', rd_0_chk)
      call hi%read_dset('rd_1', rd_1_chk)
      
      write (*, 100) 'real(dp) dataset'
      write (*, 110) rd_0 == rd_0_chk
      write (*, 120) ALL(rd_1 == rd_1_chk)
      
      !
      
      call hi%read_dset('cs_0', cs_0_chk)
      call hi%read_dset('cs_1', cs_1_chk)
      
      write (*, 100) 'complex(sp) dataset'
      write (*, 110) cs_0 == cs_0_chk
      write (*, 120) ALL(cs_1 == cs_1_chk)
      
      !
      
      call hi%read_dset('cd_0', cd_0_chk)
      call hi%read_dset('cd_1', cd_1_chk)
      
      write (*, 100) 'complex(dp) dataset'
      write (*, 110) cd_0 == cd_0_chk
      write (*, 120) ALL(cd_1 == cd_1_chk)
      
      !
      
      call hi%read_dset('l_0', l_0_chk)
      call hi%read_dset('l_1', l_1_chk)
      
      write (*, 100) 'logical dataset'
      write (*, 110) l_0 .EQV. l_0_chk
      write (*, 120) ALL(l_1 .EQV. l_1_chk)
      
      call hi%final()
      
      ! Finish
      
      return
   
   end subroutine read_file

end program test_hdf5io
