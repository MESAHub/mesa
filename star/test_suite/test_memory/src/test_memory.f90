program test_memory

  use star_def
  use star_lib
  use eos_lib
  use net_lib
  use net_def
  use kap_def

  integer :: id, ierr, i
  type(star_info), pointer :: s

  do i = 1, 2
     write(*,*) 'iteration', i
     call alloc_star(id, ierr)
     write(*,*) 'id =', id
     if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'error in alloc_star')
     call star_ptr(id, s, ierr)
     if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'error in star_ptr')
     s% inlist_fname = 'inlist'
     call read_star_job(s, '', ierr)
     if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'error in read_star_job')
     call starlib_init(s, ierr)
     if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'error in starlib_init')
     call star_set_kap_and_eos_handles(id, ierr)
     if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'error in star_set_kap_and_eos_handles')
     call star_setup(id, '', ierr)
     if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'error in star_setup')
     call star_create_pre_ms_model(id, &
          s%job%pre_ms_T_c, &
          s%job%pre_ms_guess_rho_c, &
          s%job%pre_ms_d_log10_P, &
          s%job%pre_ms_logT_surf_limit, &
          s%job%pre_ms_logP_surf_limit, &
          s%job%initial_zfracs, &
          s%job%dump_missing_metals_into_heaviest, &
          s%job%change_net .OR. s%job%change_initial_net, &
          s%job%new_net_name, &
          s%job%pre_ms_relax_num_steps, ierr)
     if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'error in star_create_pre_ms_model')
     call starlib_shutdown()
     call free_star(id, ierr)
     if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'error in free_star')
  end do

end program test_memory
