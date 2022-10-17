! declarations for other_winds

character (len = strlen) :: wind_scheme

real(dp) :: Kudritzki_scaling_factor
real(dp) :: Grafener_scaling_factor
real(dp) :: Stern51_scaling_factor

namelist /other_winds/ &
      wind_scheme, &
      Kudritzki_scaling_factor, &
      Grafener_scaling_factor, &
      Stern51_scaling_factor

! end of declarations for other_winds
