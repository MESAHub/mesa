
         real(dp) :: xtra_coef_above_xtrans
         real(dp) :: xtra_coef_below_xtrans
         real(dp) :: xtra_dist_above_xtrans
         real(dp) :: xtra_dist_below_xtrans
         
         namelist /xtrans_mesh_factor/ &
            xtra_coef_above_xtrans, &
            xtra_coef_below_xtrans, &
            xtra_dist_above_xtrans, &
            xtra_dist_below_xtrans
            