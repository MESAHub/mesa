      
      program plotter
      use create_eosPT_files, only: setup_eos
      use make_PTplots, only: make_plot_files, Do_Test
      call setup_eos
      !call Do_Test
      call make_plot_files
      end program plotter

