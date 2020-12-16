         
         s% use_other_mesh_functions = .true.
         s% how_many_other_mesh_fcns => how_many_other_mesh_fcns
         s% other_mesh_fcn_data => other_mesh_fcn_data

         ! setup for convective_bdy
         call read_inlist_convective_bdy(ierr)
