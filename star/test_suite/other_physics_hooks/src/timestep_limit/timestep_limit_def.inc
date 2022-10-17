! declarations for other_timestep_limit

real(dp) :: &
      delta_lgT_limit, &
      delta_lgT_hard_limit, &
      delta_lgT_limit_min_lgT, &
      delta_mdot_atol, &
      delta_mdot_rtol, &
      delta_mdot_limit, &
      delta_mdot_hard_limit

integer :: &
      rotation_steps_limit, &
      rotation_steps_hard_limit

namelist /other_timestep_limit/ &
      delta_lgT_limit, &
      delta_lgT_hard_limit, &
      delta_lgT_limit_min_lgT, &
      rotation_steps_limit, &
      rotation_steps_hard_limit, &
      delta_mdot_atol, &
      delta_mdot_rtol, &
      delta_mdot_limit, &
      delta_mdot_hard_limit

! end of declarations for other_timestep_limit
      
