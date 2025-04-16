module mod_colors
  use colors_def
  use math_lib
  use const_def
  use utils_lib
  use colors_ctrls_io
  use colors_lib, only: colors_init, colors_shutdown
  
  implicit none
  
  private
  
  public :: init_colors, shutdown_colors, copy_colors_to_star
  
  contains
  
  !************************************
  !******** Initialize Colors *********
  !************************************
  
  subroutine init_colors(ierr)
    integer, intent(out) :: ierr
    
    ! Initialize colors system
    call colors_init(ierr)
    if (ierr /= 0) then
      write(*,*) 'Error initializing colors module'
      return
    end if
    
    ! Read the colors controls - using the renamed global instance
    call read_colors_controls(color_settings, ierr)
    if (ierr /= 0) then
      write(*,*) 'Error reading colors controls'
      return
    end if
    
    ! Output success message and settings
    write(*,*) 'Colors module initialized successfully'
    call write_colors_controls_info(color_settings, 6)  ! 6 is stdout
  end subroutine init_colors
  
  !************************************
  !******** Shutdown Colors **********
  !************************************
  
  subroutine shutdown_colors()
    ! Call the primary shutdown function
    call colors_shutdown()
    write(*,*) 'Colors module shut down successfully'
  end subroutine shutdown_colors
  
  !************************************
  !****** Copy to Star Info **********
  !************************************
  
  ! This will be called from run_star_extras.f90 or another file that has access to star_info
  subroutine copy_colors_to_star(id)
    integer, intent(in) :: id
    ! No implementation here - this is just an interface that will be implemented elsewhere
  end subroutine copy_colors_to_star
  
end module mod_colors