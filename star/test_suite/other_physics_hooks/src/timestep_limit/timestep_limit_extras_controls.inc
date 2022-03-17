         s% ctrl% use_other_timestep_limit = .true.
         s% other_timestep_limit => set_other_timestep_limit

         ! setup for other_timestep_limit
         call read_inlist_other_timestep_limit(ierr)
