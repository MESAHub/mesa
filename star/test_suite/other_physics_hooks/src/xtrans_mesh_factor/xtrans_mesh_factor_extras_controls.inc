         
         s% ctrl% use_other_mesh_delta_coeff_factor = .true.
         s% other_mesh_delta_coeff_factor => other_mesh_delta_coeff_factor
         
         ! setup for xtrans_mesh_factor
         call read_inlist_xtrans_mesh_factor(ierr)
         