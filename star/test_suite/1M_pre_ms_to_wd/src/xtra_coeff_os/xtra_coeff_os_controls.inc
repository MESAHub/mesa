         if (s% use_other_mesh_delta_coeff_factor) then ! setup for xtra_coeff_os
            call read_inlist_xtra_coeff_os(ierr)
            if (ierr /= 0) return
            s% other_mesh_delta_coeff_factor => other_mesh_delta_coeff_factor
         end if
