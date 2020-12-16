module test_gyre_mod

  use const_def, only : dp, mesa_dir
  use gyre_lib

  implicit none

contains

  subroutine do_test ()

    real(dp), allocatable :: global_data(:)
    real(dp), allocatable :: point_data(:,:)
    integer               :: version

    integer  :: ipar(1)
    real(dp) :: rpar(1)
    integer  :: retcode

    ! Initialize

    call gyre_init('gyre.in')

    call gyre_set_constant('GYRE_DIR', TRIM(mesa_dir)//'/gyre/gyre')

    ! Read a model from file

    call gyre_read_model('model.dat')

    ! Find modes

    call gyre_get_modes(0, user_sub, ipar, rpar)
    call gyre_get_modes(1, user_sub, ipar, rpar)

    write(*,*) 'done file model'

    ! Load a model into memory

    call load_model('model.dat', global_data, point_data, version)

    call gyre_set_model(global_data, point_data, version)

    ! Find modes

    call gyre_get_modes(0, user_sub, ipar, rpar)
    call gyre_get_modes(1, user_sub, ipar, rpar)

    write(*,*) 'done memory model'

  end subroutine do_test

  !****

  subroutine user_sub (md, ipar, rpar, retcode)

    type(mode_t), intent(in) :: md
    integer, intent(inout)   :: ipar(:)
    real(dp), intent(inout)  :: rpar(:)
    integer, intent(out)     :: retcode

    integer :: n_p
    integer :: n_g
    integer :: n_pg

    ! Print out mode info

    write(*,*) md%md_p%l, md%n_p, md%n_g, md%n_pg, REAL(md%freq('UHZ')), md%E_norm()

    ! Finish

    retcode = 0

    return

  end subroutine user_sub

  !****

  subroutine load_model (file, global_data, point_data, version)

    character(LEN=*), intent(in)       :: file
    real(dp), allocatable, intent(out) :: global_data(:)
    real(dp), allocatable, intent(out) :: point_data(:,:)
    integer, intent(out)               :: version

    integer :: unit
    integer :: n
    integer :: k
    integer :: k_chk

    ! Read a model from the MESA-format file

    open(NEWUNIT=unit, FILE=file, STATUS='OLD')

    ! Read the header

    allocate(global_data(3))

    read(unit, *) n, global_data, version

    select case (version)
    case (1)
       backspace(unit)
    case (19)
    case (100)
    case default
       stop 'Unrecognized MESA file version'
    end select

    ! Read the data

    allocate(point_data(18,n))

    read_loop : do k = 1, n
       read(unit, *) k_chk, point_data(:,k)
       if(k /= k_chk) stop 'Index mismatch'
    end do read_loop

    close(unit)

    ! Finish

    return

  end subroutine load_model

end module test_gyre_mod

program test_gyre

  use test_gyre_mod

  implicit none

  call do_test()

end program test_gyre
