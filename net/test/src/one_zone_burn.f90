program one_zone_burn
   use test_net_support, only : load_libs
   use mod_one_zone_burn, only : do_one_burn
   call load_libs
   call do_one_burn('inlist_one_zone_burn', .false.)
   !call do_one_burn('inlist_one_zone_burn_const_density',.false.)
end program one_zone_burn
