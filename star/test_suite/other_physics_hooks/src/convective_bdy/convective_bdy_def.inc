      ! declarations for convective_bdy
      
         real(dp) :: &
            convective_bdy_weight, &
            convective_bdy_dq_limit, &
            convective_bdy_min_dt_yrs
         
         namelist /convective_bdy/ &
            convective_bdy_weight, &
            convective_bdy_dq_limit, &
            convective_bdy_min_dt_yrs
         
      ! end of declarations for convective_bdy
