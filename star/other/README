star/other

The routines in this directory allow you to override/modify many
physical and numerical aspects of MESA.

In your working copy of run_star_extras, replace 
    include 'standard_run_star_extras.inc'
by the contents of the included file from star/job or mesa/include.
In your work directory, do ./mk and ./rn to check that it is okay.

Don't make any edits to any of the files in star/job or mesa/include,
or to the files in star/other. You do all of this in your private copy
of run_star_extras.

There are 3 main things to do.  Make sure to do them all.

1) Copy the template routine from the other_*.f file into your copy of
run_star_extras.f and then further modify it there. (The template
routines are usually named either null_other_* or default_other_*.)
You must place it within the module (between the contains statement
and the end module statement), at the same level as the other
subroutines in run_star_extras.f.

Rename the subroutine to something descriptive for your purposes.  For
example, if you were using other_energy, you might have

    subroutine my_energy_routine(id, ierr)
       use const_def, only: Rsun
       integer, intent(in) :: id
       integer, intent(out) :: ierr
       type (star_info), pointer :: s
       integer :: k
       ierr = 0
       call star_ptr(id, s, ierr)
       if (ierr /= 0) return
       s% extra_heat(:) = 1 ! erg/g/sec
       return
    end subroutine my_energy_routine


2) Edit extra_controls to set the appropriate procedure pointer to
point to your routine.

    subroutine extras_controls(s, ierr)
       type (star_info), pointer :: s
       integer, intent(out) :: ierr
       ierr = 0

       ! this is the place to set any procedure pointers you want to change
       ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)
       s% other_energy => energy_routine
       
       ...
    end subroutine extras_controls


3) In your controls namelist, you will need to set the appropriate
flag (e.g., use_other_energy = .true.) to enable the routine.


NOTE: if you'd like to have some inlist controls for your routine,
you can use the x_ctrl array of real(dp) variables that is in &controls
e.g., in the &controls inlist, you can set
    x_ctrl(1) = my_special_param
then in your routine, you can access that by
    s% x_ctrl(1)
Before you can use s, you need to get it using the id argument.
Here's an example of how to do that -- add these lines at the start of your routine:
        use star_lib, only: star_ptr
        type (star_info), pointer :: s
        call star_ptr(id, s, ierr)
        if (ierr /= 0) then ! OOPS
           return
        end if

for integer control values, you can use x_integer_ctrl
for logical control values, you can use x_logical_ctrl
