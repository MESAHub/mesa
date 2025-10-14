! This stub fakes enough of GYRE to make sure that test cases can
! optionally use GYRE without failing to compile on systems built with
! USE_GYRE=NO in utils/makefile_header.

module gyre_mesa_m  ! stub

   use forum_m, only: RD

   implicit none

   type :: point_t
      integer  :: s
      real(RD) :: x
   end type point_t

   type :: grid_t
      type(point_t), allocatable :: pt(:)
      integer                    :: n_k
   end type grid_t

   type :: wave_t
      type(grid_t), allocatable :: gr
      integer :: n
   contains
      procedure, public :: freq
      procedure, public :: grid
      procedure, public :: xi_r
      procedure, public :: xi_h
      procedure, public :: dW_dx
   end type wave_t

   type, extends (wave_t) :: mode_t
      integer :: n_pg
      integer :: n_p
      integer :: n_g
   end type mode_t

   interface set_constant
      module procedure set_constant_r_
      module procedure set_constant_c_
   end interface set_constant

contains

   function xi_r(this, k)

      class(wave_t), intent(in) :: this
      integer, intent(in)       :: k
      complex(RD)               :: xi_r

      xi_r = 0._RD

   end function xi_r


  function xi_h(this, k)

     class(wave_t), intent(in) :: this
     integer, intent(in)       :: k
     complex(RD)               :: xi_h

     xi_h = 0._RD

  end function xi_h


  function dW_dx(this, k)

     class(wave_t), intent(in) :: this
     integer, intent(in)       :: k
     real(RD)                  :: dW_dx

     dW_dx = 0._RD

  end function dW_dx


  function freq(this, freq_units, freq_frame)

     class(wave_t), intent(in)          :: this
     character(*), intent(in)           :: freq_units
     character(*), optional, intent(in) :: freq_frame
     complex(RD)                        :: freq

     freq = 0._RD

  end function freq


  function grid(this) result (gr)

     class(wave_t), intent(in) :: this
     type(grid_t)              :: gr

     ! Return the wave's grid

     gr = this%gr

     return

  end function grid


  subroutine init(file)

     character(*), intent(in) :: file
     write(*,*) ''
     write(*,*) 'WARNING: You are initializing the stub version of GYRE, which means'
     write(*,*) 'that MESA was installed without the embedded interfaces to GYRE.'
     write(*,*) 'Calls to GYRE in this run will not produce useful output.'
     write(*,*) ''

  end subroutine init


  subroutine final()
  end subroutine final


  subroutine read_model (file)

     character(LEN=*), intent(in) :: file

  end subroutine read_model


  subroutine set_model (global_data, point_data, version)

     real(RD), intent(in) :: global_data(:)
     real(RD), intent(in) :: point_data(:,:)
     integer, intent(in)  :: version

  end subroutine set_model


  subroutine get_modes (l, user_sub, ipar, rpar)

     integer, intent(in)     :: l
     interface
        subroutine user_sub (md, ipar, rpar, retcode)
           import mode_t
           import RD
           implicit none
           type(mode_t), intent(in) :: md
           integer, intent(inout)   :: ipar(:)
           real(RD), intent(inout)  :: rpar(:)
           integer, intent(out)     :: retcode
        end subroutine user_sub
     end interface
     integer, intent(inout)  :: ipar(:)
     real(RD), intent(inout) :: rpar(:)

  end subroutine get_modes


  subroutine set_constant_r_ (name, value)

     character(*), intent(in) :: name
     real(RD), intent(in)     :: value

  end subroutine set_constant_r_


  subroutine set_constant_c_ (name, value)

     character(*), intent(in) :: name
     character(*), intent(in) :: value

  end subroutine set_constant_c_

end module gyre_mesa_m  ! stub
