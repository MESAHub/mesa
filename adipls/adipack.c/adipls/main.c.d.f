      program runadi
c
c  main programme for adiabatic pulsations
c
c  NOTE: all setups have been moved to s/r setups_adi.
c
c  Dated: 10/3/90
c
      implicit double precision (a-h,o-z)
      dimension x_arg(1), aa_arg(1), data_arg(1)
c
c  unit numbers read as input parameters in evolution part of code.
c  Used here to suppress output to istdpr regardless of input value
c
      common/cstdio_in/ istdin_in, istdpr_in
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  initialize unit for input and output (lest they are used before
c  being set from parameter input)
c
      data istdin_in, istdpr_in /5, 6/
c
      call setups_adi
c
      i_paramset=0
      i_inout=1
      nn_arg=0
      iaa_arg=1
      ivarmd=1
c
      call adipls(i_paramset, ierr_param, i_inout, 
     *  x_arg, aa_arg, data_arg, nn_arg, ivarmd, iaa_arg)
c
      write(istdou,'(//'' return from s/r adipls with ierr_param ='',
     *  i3)') ierr_param
c
      stop
      end
