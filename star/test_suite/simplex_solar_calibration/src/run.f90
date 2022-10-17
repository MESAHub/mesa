program run
   use simplex_search_run_support, only : do_run_star_simplex
   use run_star_extras, only : extras_controls
   call do_run_star_simplex(extras_controls)
end program
